/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tvm/ImpulseFunction.h>

#include <mc_tvm/RobotFrame.h>

#include <numeric>

namespace mc_tvm
{

ImpulseFunction::ImpulseFunction(const mc_rbdyn::Robot & robot, const mc_rbdyn::RobotFrame & frame, const Eigen::Vector3d normal, double lambda)
: tvm::function::abstract::LinearFunction(robot.mb().nrDof()), robot_(robot), frame_(frame), normal_(normal), lambda(lambda)
  , jac_(frame.tvm_frame().rbdJacobian()), coriolis_calculator_(rbd::Coriolis(robot_.mb()))
{
  assert(frame_->robot().robotIndex() == robot_.robotIndex() && "ImpulseFunction frame must belong to the robot");
  // assert(axis_ >= 0 && axis_ < 3 && "Axis must be in [0, 2] as there are only three translational axis in the frame (1=x, 2=y, 3=z)");
  registerUpdates(Update::B, &ImpulseFunction::updateb);
  registerUpdates(Update::Jacobian, &ImpulseFunction::updateJacobian);
  addOutputDependency<ImpulseFunction>(Output::B, Update::B);
  addOutputDependency<ImpulseFunction>(Output::Jacobian, Update::Jacobian);
  auto & tvm_robot = robot.tvmRobot();
  addInputDependency<ImpulseFunction>(Update::Jacobian, tvm_robot, Robot::Output::H);
  addInputDependency<ImpulseFunction>(Update::Jacobian, tvm_robot, Robot::Output::C);
  addInputDependency<ImpulseFunction>(Update::B, tvm_robot, Robot::Output::H); // TODO check if needed
  q_ddot_var_ = tvm::dot(tvm_robot.q(), 2);
  addVariable(q_ddot_var_, true);
  jacobian_[q_ddot_var_.get()] = Eigen::MatrixXd::Identity(robot.mb().nrDof(), robot.mb().nrDof());
  jacobian_[q_ddot_var_.get()].properties(tvm::internal::MatrixProperties::IDENTITY);

  addVariable(tvm::dot(tvm_robot.q(), 1), true);
  velocity_.setZero();

  Eigen::MatrixXd P_n_sub = normal_ * normal_.transpose();
  P_n = Eigen::Matrix<double, 6, 6>::Zero();
  P_n.block<3,3>(0,0) = P_n_sub;

  tau_imp = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  tau_imp2 = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  tau_imp_const = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  tau_imp_act = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  tau_imp_deriv = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  tau_imp_deriv_term1 = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  tau_imp_deriv_term2 = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  tau_imp_deriv_term3 = Eigen::VectorXd::Zero(robot_.mb().nrDof());

  last_joint_velocities_ = Eigen::VectorXd::Zero(robot_.mb().nrDof());

  // Initialize the moving average filter window
  // for(int i = 0; i < acc_error_window_size_; ++i)
  // {
  //   acc_error_window_.emplace_back(Eigen::VectorXd::Ones(robot_.mb().nrDof()));
  // }

  alpha_s_ = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  q_d = Eigen::VectorXd::Zero(robot_.mb().nrDof());

}


void ImpulseFunction::updateb() // TODO possibly make this function dependent on updateJacobian and use the jacobians etc in class variables rather than local function variables
{
  const auto & robot = robot_.tvmRobot();
  auto M = robot.H();

  // Check the mass matrix before inversion (for better debugging)
  assert(M.rows() == robot_.mb().nrDof());
  assert(M.cols() == robot_.mb().nrDof());
  assert(std::abs(M.determinant()) > 1e-5 && "Mass matrix is singular");

  // Get jacobian of the robot frame
  const rbd::MultiBody robot_mb = robot_.mb();
  const rbd::MultiBodyConfig mbc = robot_.mbc();
  const auto & world_frame_jacobian = jac_.jacobian(robot_mb, mbc);
  Eigen::MatrixXd full_world_frame_jacobian(6, robot_.mb().nrDof());
  jac_.fullJacobian(robot_mb, world_frame_jacobian, full_world_frame_jacobian);

  const auto & world_frame_jacobian_dot = jac_.jacobianDot(robot_mb, mbc);
  Eigen::MatrixXd full_world_frame_jacobian_dot(6, robot_.mb().nrDof());
  jac_.fullJacobian(robot_mb, world_frame_jacobian_dot, full_world_frame_jacobian_dot);

  assert(full_world_frame_jacobian.cols() == robot_.mb().nrDof());
  
  // Generate necessary matrices
  Eigen::MatrixXd j_m_uninverted = full_world_frame_jacobian * M.inverse() * full_world_frame_jacobian.transpose();
  assert(std::abs(j_m_uninverted.determinant()) > 1e-5 && "j_m_uninverted matrix is singular");
  Eigen::MatrixXd j_m = full_world_frame_jacobian.transpose() * j_m_uninverted.inverse();

  Eigen::MatrixXd J_=full_world_frame_jacobian;
  Eigen::MatrixXd J_d_=full_world_frame_jacobian_dot;

  Eigen::MatrixXd C = coriolis_calculator_.coriolis(robot_.mb(), robot_.mbc());

  Eigen::MatrixXd M_d_ = C+C.transpose();

  Eigen::MatrixXd J_dq_new = (-1.f/lambda)*(J_d_.transpose()*j_m_uninverted.inverse() -
    j_m*(J_d_*M.inverse()*J_.transpose()-J_*M.inverse()*M_d_*M.inverse()*J_.transpose() +
      J_*M.inverse()*J_d_.transpose())*j_m_uninverted.inverse())*P_n*J_ + (-1.f/lambda)*j_m*P_n*J_d_ - j_m*P_n*J_;

  q_d = tvm::dot(robot.q(),1)->value();

  b_ = J_dq_new * q_d;

  // These are now used as constants but if used in a final version should be taken in initialization from input parameters
  double cres = 1;
  double delta_t_ = 0.001;
  double timestep = 0.002;

  // Take the numerical derivative of the joint velocities
  num_qdd = (q_d - last_joint_velocities_) / timestep;
  last_joint_velocities_ = q_d;

  alpha_s_ += robot_.tvmRobot().alphaD()->value()*0.002;

  tau_imp2 = -1.f*j_m*P_n*J_*alpha_s_;
  tau_imp = -1.f*j_m*P_n*J_*q_d;
  tau_imp_act = (-1.f*(cres+1)/delta_t_)*j_m*P_n*J_*q_d;
  tau_imp_deriv = (-1.f*(cres+1)/delta_t_)*((J_d_.transpose()*j_m_uninverted.inverse() -
    j_m*(J_d_*M.inverse()*J_.transpose()-J_*M.inverse()*M_d_*M.inverse()*J_.transpose() +
    J_*M.inverse()*J_d_.transpose())*j_m_uninverted.inverse())*P_n*J_ * tvm::dot(robot.q(), 1)->value()
    + j_m*P_n*(J_d_ * tvm::dot(robot.q(), 1)->value() + J_ * tvm::dot(robot.q(), 2)->value()));

  tau_imp_const += (delta_t_/(cres+1))*tau_imp_deriv * 0.002;

  tau_imp_deriv_term1 = (-1.f*(cres+1)/delta_t_)*((J_d_.transpose()*j_m_uninverted.inverse() -
    j_m*(J_d_*M.inverse()*J_.transpose()-J_*M.inverse()*M_d_*M.inverse()*J_.transpose() +
    J_*M.inverse()*J_d_.transpose())*j_m_uninverted.inverse())*P_n*J_ * q_d);
  tau_imp_deriv_term2 = (-1.f*(cres+1)/delta_t_)*j_m*P_n*J_d_ * q_d;
  tau_imp_deriv_term3 = (-1.f*(cres+1)/delta_t_)*j_m*P_n*J_ * tvm::dot(robot.q(), 2)->value();

  tau_imp_deriv_num = (-1.f*(cres+1)/delta_t_)*((J_d_.transpose()*j_m_uninverted.inverse() -
    j_m*(J_d_*M.inverse()*J_.transpose()-J_*M.inverse()*M_d_*M.inverse()*J_.transpose() +
    J_*M.inverse()*J_d_.transpose())*j_m_uninverted.inverse())*P_n*J_ * q_d
    + j_m*P_n*(J_d_ * q_d + J_ * num_qdd));
}

void ImpulseFunction::updateJacobian()
{
  const auto & robot = robot_.tvmRobot();
  auto M = robot.H();

  // Check the mass matrix before inversion (for better debugging)
  assert(M.rows() == robot_.mb().nrDof());
  assert(M.cols() == robot_.mb().nrDof());
  assert(std::abs(M.determinant()) > 1e-5 && "Mass matrix is singular");

  // Get jacobian of the robot frame
  const rbd::MultiBody robot_mb = robot_.mb();
  const rbd::MultiBodyConfig mbc = robot_.mbc();
  const auto & world_frame_jacobian = jac_.jacobian(robot_mb, mbc);
  Eigen::MatrixXd full_world_frame_jacobian(6, robot_.mb().nrDof());
  jac_.fullJacobian(robot_mb, world_frame_jacobian, full_world_frame_jacobian);

  assert(full_world_frame_jacobian.cols() == robot_.mb().nrDof());

  // // Generate necessary matrices
  Eigen::MatrixXd j_m_uninverted = /*full_world_frame_jacobian.transpose() * */(full_world_frame_jacobian * M.inverse() * full_world_frame_jacobian.transpose());
  assert(std::abs(j_m_uninverted.determinant()) > 1e-5 && "j_m_uninverted is singular");
  Eigen::MatrixXd j_m_before_premult_J = j_m_uninverted.inverse();

Eigen::MatrixXd j_m = full_world_frame_jacobian.transpose() * j_m_uninverted.inverse();

  // // Get the average acceleration error
  // Eigen::VectorXd current_acc = robot_.tvmRobot().alphaD()->value();
  // assert(current_acc.size() == num_qdd.size());
  //
  // Eigen::VectorXd acc_error = Eigen::VectorXd::Zero(robot_.mb().nrDof()); // num_qdd.cwiseQuotient(current_acc);
  // for (int i = 0; i < current_acc.size(); ++i)
  // {
  //   if(std::abs(current_acc(i)) > 1e-6)
  //   {
  //     acc_error(i) = num_qdd(i) / current_acc(i);
  //   }
  //   else
  //   {
  //     acc_error(i) = 1.0;
  //   }
  //   if(acc_error(i) < 0.3)
  //   {
  //     acc_error(i) = 0.3;
  //   } else if (acc_error(i) > 4.0)
  //   {
  //     acc_error(i) = 4.0;
  //   }
  // }
  // acc_error_window_.push_back(acc_error);
  // if (acc_error_window_.size() > acc_error_window_size_) {
  //   acc_error_window_.erase(acc_error_window_.begin());
  // }
  // assert(acc_error_window_.size() == acc_error_window_size_);
  // Eigen::VectorXd sum = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  // for (auto i : acc_error_window_)
  // {
  //   sum += i;
  // }
  // Eigen::VectorXd average_error = sum / acc_error_window_.size();

  J_ddq = /*average_error.asDiagonal()**/((-1.f/lambda) * j_m * P_n * full_world_frame_jacobian);// + 0.5*0.002*J_dq_new;              // multiplies ddq variable
  // J_ddq = -1.f * j_m * P_n * full_world_frame_jacobian;              // multiplies ddq variable


  splitJacobian(J_ddq, q_ddot_var_);
}

Eigen::VectorXd & ImpulseFunction::JointPos()
{
  q = robot_.tvmRobot().q()->value();
  return q;
}

Eigen::VectorXd & ImpulseFunction::JointVel()
{
  alpha = robot_.tvmRobot().alpha()->value();
  return alpha;
}

Eigen::VectorXd & ImpulseFunction::JointVels()
{
  return alpha_s_;
}

Eigen::VectorXd & ImpulseFunction::JointAcc()
{
  alpha_d = robot_.tvmRobot().alphaD()->value();
  return alpha_d;
}

Eigen::VectorXd & ImpulseFunction::JointVelUsed()
{
  return q_d;
}

Eigen::VectorXd & ImpulseFunction::JointAccUsed()
{
  q_dd = tvm::dot(robot_.tvmRobot().q(),2)->value();
  return q_dd;
}

Eigen::VectorXd & ImpulseFunction::JointAccNum()
{
  return num_qdd;
}

} // namespace mc_tvm


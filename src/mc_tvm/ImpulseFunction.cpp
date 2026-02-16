/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tvm/ImpulseFunction.h>

#include <mc_tvm/RobotFrame.h>

#include <numeric>

namespace mc_tvm
{

ImpulseFunction::ImpulseFunction(const mc_rbdyn::Robot & robot, const mc_rbdyn::RobotFrame & frame, const Eigen::Vector3d normal, double lambda_high, double lambda_low, double c_res, double delta_t, const Eigen::VectorXd & limit_high, const Eigen::VectorXd & limit_low, bool enforce_high_limit)
: tvm::function::abstract::LinearFunction(robot.mb().nrDof()), robot_(robot), frame_(frame), normal_(normal), lambda_high(lambda_high), lambda_low(lambda_low), c_res_(c_res), delta_t_(delta_t), limit_high_(limit_high), limit_low_(limit_low), enforce_high_limit_(enforce_high_limit)
  , jac_(frame.tvm_frame().rbdJacobian()), coriolis_calculator_(rbd::Coriolis(robot_.mb()))
{
  assert(frame_->robot().robotIndex() == robot_.robotIndex() && "ImpulseFunction frame must belong to the robot");
  auto & tvm_robot = robot.tvmRobot();
  registerUpdates(Update::B, &ImpulseFunction::updateb);
  registerUpdates(Update::Jacobian, &ImpulseFunction::updateJacobian);
  addOutputDependency<ImpulseFunction>(Output::B, Update::B);
  addOutputDependency<ImpulseFunction>(Output::Jacobian, Update::Jacobian);
  addInputDependency<ImpulseFunction>(Update::Jacobian, tvm_robot, Robot::Output::H);
  addInputDependency<ImpulseFunction>(Update::Jacobian, tvm_robot, Robot::Output::C);
  addInputDependency<ImpulseFunction>(Update::B, tvm_robot, Robot::Output::H); // TODO check if needed
  addVariable(tvm_robot.alphaD(), true);

  startParam = tvm_robot.qFloatingBase()->size();

  pre_multiplier_ = Eigen::MatrixXd::Identity(robot.mb().nrDof(), robot.mb().nrDof());
  pre_multiplier_.block<6, 6>(0, 0).setZero();

  b_ = Eigen::VectorXd::Zero(robot.mb().nrDof());

  Eigen::MatrixXd P_n_sub = normal_ * normal_.transpose();
  P_n = Eigen::Matrix<double, 6, 6>::Zero();
  P_n.block<3,3>(0,0) = P_n_sub;

  tau_imp_pred = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  tau_imp_act = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  tau_imp_deriv = Eigen::VectorXd::Zero(robot_.mb().nrDof());

  last_joint_velocities_ = Eigen::VectorXd::Zero(robot_.mb().nrDof());

  q_d = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  num_qd = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  num_qdd = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  lambda = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  diff_upper_ = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  diff_lower_ = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  high_lambda_latch_ = std::vector<bool>(robot_.mb().nrDof(), false);
  high_lambda_latch_count_ = Eigen::VectorXi::Zero(robot_.mb().nrDof());

  constraint_right_side_ = Eigen::VectorXd::Zero(robot_.mb().nrDof());
}


void ImpulseFunction::updateb() // TODO possibly make this function dependent on updateJacobian and use the jacobians etc in class variables rather than local function variables
{
  // mc_rtc::log::info("Landed in the updateb of ImpulseFunction");
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

  Eigen::MatrixXd J_dq_new = (-1.f) * ((c_res_+1.f)/(delta_t_/**lambda*/)) *(J_d_.transpose()*j_m_uninverted.inverse() -
    j_m*(J_d_*M.inverse()*J_.transpose()-J_*M.inverse()*M_d_*M.inverse()*J_.transpose() +
      J_*M.inverse()*J_d_.transpose())*j_m_uninverted.inverse())*P_n*J_ +  (-1.f) * ((c_res_+1.f)/(delta_t_/**lambda*/)) *j_m*P_n*J_d_;// - j_m*P_n*J_;


  q_d = tvm::dot(robot.q(),1)->value();
  q_dd = tvm::dot(robot.q(),2)->value();

  b_ = pre_multiplier_*J_dq_new * q_d;

  // These are now used as constants but if used in a final version should be taken in initialization from input parameters
  double timestep = 0.002;

  auto current_vel = robot_.encoderVelocities();

  num_qd(0) = 0.0;
  num_qd(1) = 0.0;
  num_qd(2) = 0.0;
  num_qd(3) = 0.0;
  num_qd(4) = 0.0;
  num_qd(5) = 0.0;
  num_qd(6)  = current_vel.at(6);
  num_qd(7)  = current_vel.at(7);
  num_qd(8)  = current_vel.at(8);
  num_qd(9)  = current_vel.at(9);
  num_qd(10) = current_vel.at(10);
  num_qd(11) = current_vel.at(11);
  num_qd(12) = current_vel.at(0);
  num_qd(13) = current_vel.at(1);
  num_qd(14) = current_vel.at(2);
  num_qd(15) = current_vel.at(3);
  num_qd(16) = current_vel.at(4);
  num_qd(17) = current_vel.at(5);
  num_qd(18) = current_vel.at(12);
  num_qd(19) = current_vel.at(13);
  num_qd(20) = current_vel.at(14);
  num_qd(21) = current_vel.at(15);
  num_qd(22) = current_vel.at(16);
  num_qd(23) = current_vel.at(35);
  num_qd(24) = current_vel.at(36);
  num_qd(25) = current_vel.at(37);
  num_qd(26) = current_vel.at(38);
  num_qd(27) = current_vel.at(39);
  num_qd(28) = current_vel.at(40);
  num_qd(29) = current_vel.at(41);
  num_qd(30) = current_vel.at(42);
  num_qd(31) = current_vel.at(43);
  num_qd(32) = current_vel.at(17);
  num_qd(33) = current_vel.at(18);
  num_qd(34) = current_vel.at(19);
  num_qd(35) = current_vel.at(20);
  num_qd(36) = current_vel.at(21);
  num_qd(37) = current_vel.at(22);
  num_qd(38) = current_vel.at(23);
  num_qd(39) = current_vel.at(24);
  num_qd(40) = current_vel.at(25);

  // Take the numerical derivative of the joint velocities
  for (int i = 0; i < robot_.mb().nrDof(); ++i)
  {
    double current_speed = num_qd(i);
    double past_speed = last_joint_velocities_(i);
    num_qdd(i) = (current_speed - past_speed) / timestep;
  }
  last_joint_velocities_ = num_qd;

  // now add the limits termwise, lambda*sng(tau_I_max - tau_I)*sqrt(tau_I_max - tau_I)
  tau_imp_act = (-1.f*(c_res_+1)/delta_t_)*j_m*P_n*J_*num_qd;
  tau_imp_pred = (-1.f*(c_res_+1)/delta_t_)*j_m*P_n*J_*q_d;

  assert(b_.size() == robot_.mb().nrDof());
  assert(limit_high_.size() == robot_.mb().nrDof());
  assert(limit_low_.size() == robot_.mb().nrDof());
  assert(tau_imp_pred.size() == robot_.mb().nrDof());
  getLambda();
  if(enforce_high_limit_)
  {
    for (int i = startParam; i < b_.size(); ++i)
    {
      if (diff_upper_(i) >= 0 /*|| diff2 <= 0*/)
      {
        b_(i) += /*std::sqrt*/(lambda(i))*/*std::sqrt*/(diff_upper_(i));
        constraint_right_side_(i) = -1.0*lambda(i)*/*std::sqrt*/(diff_upper_(i));
      } else if (diff_lower_(i) < 0)
      {
        b_(i) -= /*std::sqrt*/(lambda(i))*/*std::sqrt*/(-1.0*diff_upper_(i));
        constraint_right_side_(i) = -1.0*lambda(i)*diff_upper_(i);
      }else
      {
        b_(i) -= /*std::sqrt*/(lambda(i))*/*std::sqrt*/(-1.f*diff_upper_(i));
        constraint_right_side_(i) = -1.0*lambda(i)*diff_upper_(i);
      }
    }
  } else
  {
    for (int i = startParam; i < b_.size(); ++i)
    {
      if (diff_lower_(i) < 0/* || diff2 >= 0*/)
      {
        b_(i) -= /*std::sqrt*/(lambda(i))*/*std::sqrt*/(-1.0*diff_lower_(i));
        constraint_right_side_(i) = -1.0*lambda(i)*/*std::sqrt*/(diff_lower_(i));
      } else if (diff_upper_(i) >= 0)
      {
        b_(i) += /*std::sqrt*/(lambda(i))*/*std::sqrt*/(diff_lower_(i));
        constraint_right_side_(i) = -1.0*lambda(i)*diff_lower_(i);
      } else
      {
        b_(i) += /*std::sqrt*/(lambda(i))*/*std::sqrt*/(diff_lower_(i));
        constraint_right_side_(i) = -1.0*lambda(i)*diff_lower_(i);
      }
    }
  }


   // last_joint_velocities_ = q;

   // alpha_s_ += robot_.tvmRobot().alphaD()->value()*0.002;

   // tau_imp2 = -1.f*j_m*P_n*J_*alpha_s_;
   // tau_imp_act = (-1.f*(c_res_+1)/delta_t_)*j_m*P_n*J_*q_d;
   tau_imp_deriv = (-1.f*(c_res_+1)/delta_t_)*((J_d_.transpose()*j_m_uninverted.inverse() -
     j_m*(J_d_*M.inverse()*J_.transpose()-J_*M.inverse()*M_d_*M.inverse()*J_.transpose() +
     J_*M.inverse()*J_d_.transpose())*j_m_uninverted.inverse())*P_n*J_ * q_d
     + j_m*P_n*(J_d_ * q_d + J_ * q_dd/*q_ddot_var_->value()*/));

   tau_imp_deriv_num = (-1.f*(c_res_+1)/delta_t_)*((J_d_.transpose()*j_m_uninverted.inverse() -
     j_m*(J_d_*M.inverse()*J_.transpose()-J_*M.inverse()*M_d_*M.inverse()*J_.transpose() +
     J_*M.inverse()*J_d_.transpose())*j_m_uninverted.inverse())*P_n*J_ * q_d
     + j_m*P_n*(J_d_ * q_d + J_ * num_qdd));
}

void ImpulseFunction::updateJacobian()
{
  // mc_rtc::log::info("Landed in the updateJacobian of ImpulseFunction");
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

  J_ddq = -1.f * ((c_res_+1.f)/(delta_t_/**lambda*/)) * pre_multiplier_ * j_m * P_n * full_world_frame_jacobian;// + 0.5*0.002*J_dq_new;              // multiplies ddq variable
  // J_ddq = -1.f * j_m * P_n * full_world_frame_jacobian;              // multiplies ddq variable

  splitJacobian(J_ddq, robot.alphaD());
}

void ImpulseFunction::getLambda()
{
  double lambda_increment = (lambda_high - lambda_low)/static_cast<double>(lambda_growing_steps);
  for (int i = 0; i < robot_.mb().nrDof(); ++i)
  {
    // double tau_I = tau_imp_act(i);
    diff_upper_(i) = tau_imp_pred(i) - limit_multiplier_*limit_high_(i);
    diff_lower_(i) = tau_imp_pred(i) - limit_multiplier_*limit_low_(i);
    if(diff_lower_(i) < 0 || diff_upper_(i) > 0)
    {
      if(lambda(i) < lambda_high)
      {
        lambda(i) += lambda_increment;
      } else
      {
        lambda(i) = lambda_high;
      }
      high_lambda_latch_[i] = true;
    } else
    {
      if(high_lambda_latch_[i])
      {
        high_lambda_latch_count_(i) += 1;
        if(high_lambda_latch_count_(i) >= high_lambda_latch_max_)
        {
          high_lambda_latch_[i] = false;
          high_lambda_latch_count_(i) = 0;
        }
        if(lambda(i) < lambda_high)
        {
          lambda(i) += lambda_increment;
        } else
        {
          lambda(i) = lambda_high;
        }
      } else
      {
        if(lambda(i) > lambda_low)
        {
          lambda(i) -= lambda_increment;
        } else
        {
          lambda(i) = lambda_low;
        }
      }
    }
  }


}

// Eigen::VectorXd & ImpulseFunction::JointPos()
// {
//   q = robot_.tvmRobot().q()->value();
//   return q;
// }
//
// Eigen::VectorXd & ImpulseFunction::JointVel()
// {
//   alpha = robot_.tvmRobot().alpha()->value();
//   return alpha;
// }
//
// Eigen::VectorXd & ImpulseFunction::JointVels()
// {
//   return alpha_s_;
// }
//
Eigen::VectorXd & ImpulseFunction::JointAcc()
{
  alpha_d = robot_.tvmRobot().alphaD()->value();
  return alpha_d;
}
//
// Eigen::VectorXd & ImpulseFunction::JointVelUsed()
// {
//   return q_d;
// }
//
// Eigen::VectorXd & ImpulseFunction::JointAccUsed()
// {
//   q_dd = q_ddot_var_->value();
//   return q_dd;
// }
//
Eigen::VectorXd & ImpulseFunction::JointAccNum()
{
  return num_qdd;
}

Eigen::VectorXd & ImpulseFunction::JointVelNum()
{
  return num_qd;
}

} // namespace mc_tvm


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

  b_ = Eigen::VectorXd::Zero(robot.mb().nrDof());

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

  alpha_s_ = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  q_d = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  num_qdd = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  lambda = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  diff_upper_ = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  diff_lower_ = Eigen::VectorXd::Zero(robot_.mb().nrDof());

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

  q_d = (Eigen::VectorXd)tvm::dot(robot.q(),1)->value();
  q_dd = tvm::dot(robot.q(),2)->value();

  b_ = J_dq_new * q_d;

  // now add the limits termwise, lambda*sng(tau_I_max - tau_I)*sqrt(tau_I_max - tau_I)
  tau_imp_act = (-1.f*(c_res_+1)/delta_t_)*j_m*P_n*J_*q_d;
  assert(b_.size() == robot_.mb().nrDof());
  assert(limit_high_.size() == robot_.mb().nrDof());
  assert(limit_low_.size() == robot_.mb().nrDof());
  assert(tau_imp.size() == robot_.mb().nrDof());
  getLambda();
  if(enforce_high_limit_)
  {
    for (int i = 0; i < b_.size(); ++i)
    {
      if (diff_upper_(i) >= 0 /*|| diff2 <= 0*/)
      {
        b_(i) += /*std::sqrt*/(lambda(i))*/*std::sqrt*/(diff_upper_(i));
        constraint_right_side_(i) = -1.0*lambda(i)*diff_upper_(i);
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
    for (int i = 0; i < b_.size(); ++i)
    {
      if (diff_lower_(i) < 0/* || diff2 >= 0*/)
      {
        b_(i) -= /*std::sqrt*/(lambda(i))*/*std::sqrt*/(-1.0*diff_lower_(i));
        constraint_right_side_(i) = -1.0*lambda(i)*diff_lower_(i);
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

   // These are now used as constants but if used in a final version should be taken in initialization from input parameters
   double timestep = 0.002;

   // Take the numerical derivative of the joint velocities
  for (int i = 0; i < robot_.mb().nrDof(); ++i)
  {
    double current_speed = q_d(i);
    double past_speed = last_joint_velocities_(i);
    num_qdd(i) = (current_speed - past_speed) / timestep;
  }
   // num_qdd = (q_d - last_joint_velocities_) / timestep;
   last_joint_velocities_ = q_d;

   // alpha_s_ += robot_.tvmRobot().alphaD()->value()*0.002;

   // tau_imp2 = -1.f*j_m*P_n*J_*alpha_s_;
   tau_imp = -1.f*j_m*P_n*J_*q_d;
   // tau_imp_act = (-1.f*(c_res_+1)/delta_t_)*j_m*P_n*J_*q_d;
   tau_imp_deriv = (-1.f*(c_res_+1)/delta_t_)*((J_d_.transpose()*j_m_uninverted.inverse() -
     j_m*(J_d_*M.inverse()*J_.transpose()-J_*M.inverse()*M_d_*M.inverse()*J_.transpose() +
     J_*M.inverse()*J_d_.transpose())*j_m_uninverted.inverse())*P_n*J_ * q_d
     + j_m*P_n*(J_d_ * q_d + J_ * q_dd/*q_ddot_var_->value()*/));

   // tau_imp_const += (delta_t_/(c_res_+1))*tau_imp_deriv * 0.002;

   // tau_imp_deriv_term1 = (-1.f*(c_res_+1)/delta_t_)*((J_d_.transpose()*j_m_uninverted.inverse() -
   //   j_m*(J_d_*M.inverse()*J_.transpose()-J_*M.inverse()*M_d_*M.inverse()*J_.transpose() +
   //   J_*M.inverse()*J_d_.transpose())*j_m_uninverted.inverse())*P_n*J_ * q_d);
   // tau_imp_deriv_term2 = (-1.f*(c_res_+1)/delta_t_)*j_m*P_n*J_d_ * q_d;
   // tau_imp_deriv_term3 = (-1.f*(c_res_+1)/delta_t_)*j_m*P_n*J_ * q_dd/*q_ddot_var_->value()*/;

   tau_imp_deriv_num = (-1.f*(c_res_+1)/delta_t_)*((J_d_.transpose()*j_m_uninverted.inverse() -
     j_m*(J_d_*M.inverse()*J_.transpose()-J_*M.inverse()*M_d_*M.inverse()*J_.transpose() +
     J_*M.inverse()*J_d_.transpose())*j_m_uninverted.inverse())*P_n*J_ * q_d
     + j_m*P_n*(J_d_ * q_d + J_ * num_qdd));
     // mc_rtc::log::info("Exiting the updateb of ImpulseFunction");
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

  J_ddq = -1.f * ((c_res_+1.f)/(delta_t_/**lambda*/)) * j_m * P_n * full_world_frame_jacobian;// + 0.5*0.002*J_dq_new;              // multiplies ddq variable
  // J_ddq = -1.f * j_m * P_n * full_world_frame_jacobian;              // multiplies ddq variable


  splitJacobian(J_ddq, robot.alphaD());
}

void ImpulseFunction::getLambda()
{
  double lambda_increment = (lambda_high - lambda_low)/static_cast<double>(lambda_growing_steps);
  for (int i = 0; i < robot_.mb().nrDof(); ++i)
  {
    // double tau_I = tau_imp_act(i);
    diff_upper_(i) = tau_imp_act(i) - limit_multiplier_*limit_high_(i);
    diff_lower_(i) = tau_imp_act(i) - limit_multiplier_*limit_low_(i);
    if(diff_lower_(i) < 0 || diff_upper_(i) > 0)
    {
      if(lambda(i) < lambda_high)
      {
        lambda(i) += lambda_increment;
      } else
      {
        lambda(i) = lambda_high;
      }
      high_lambda_latch_ = true;
    } else
    {
      if(high_lambda_latch_)
      {
        high_lambda_latch_count_ += 1;
        if(high_lambda_latch_count_ >= high_lambda_latch_max_)
        {
          high_lambda_latch_ = false;
          high_lambda_latch_count_ = 0;
        }
        lambda(i) += lambda_increment;
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

} // namespace mc_tvm


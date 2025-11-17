/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tvm/ImpulseFunction.h>

#include <mc_tvm/Robot.h>
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
  q_ddot_var_ = tvm_robot.alphaD();
  addVariable(q_ddot_var_, true);
  velocity_.setZero();

  Eigen::MatrixXd P_n_sub = normal_ * normal_.transpose();
  P_n = Eigen::Matrix<double, 6, 6>::Zero();
  P_n.block<3,3>(0,0) = P_n_sub;

}


void ImpulseFunction::updateb() // TODO possibly make this function dependent on updateJacobian and use the jacobians etc in class variables rather than local function variables
{
  const auto & robot = robot_.tvmRobot();
  auto M = robot.H();
  // auto C = robot.C();

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

  // Get normal of hitting plane
  // auto normal = robot_.frame(frame_->name()).position().rotation().col(axis_);

  // Generate necessary matrices
  Eigen::MatrixXd j_m_uninverted = full_world_frame_jacobian * M.inverse() * full_world_frame_jacobian.transpose();
  assert(std::abs(j_m_uninverted.determinant()) > 1e-5 && "j_m_uninverted matrix is singular");
  Eigen::MatrixXd j_m = full_world_frame_jacobian.transpose() * j_m_uninverted.inverse();
  // Eigen::MatrixXd P_n_sub = normal_ * normal_.transpose();
  // Eigen::Matrix<double, 6, 6> P_n = Eigen::Matrix<double, 6, 6>::Zero();
  // P_n.block<3,3>(0,0) = P_n_sub;


  // J_dq  = j_m * P_n * full_world_frame_jacobian / delta_t_;

  Eigen::MatrixXd J_=full_world_frame_jacobian;
  Eigen::MatrixXd J_d_=full_world_frame_jacobian_dot;


  auto C = coriolis_calculator_.coriolis(robot_.mb(), robot_.mbc());

  Eigen::MatrixXd M_d_ = C+C.transpose();

  Eigen::MatrixXd J_dq_new = (J_d_.transpose()*j_m_uninverted.inverse() -
    j_m*(J_d_*M.inverse()*J_.transpose()-J_*M.inverse()*M_d_*M.inverse()*J_.transpose() +
      J_*M.inverse()*J_d_.transpose())*j_m_uninverted.inverse())*P_n*J_ + j_m*P_n*J_d_ + lambda*j_m*P_n*J_;

  b_ = J_dq_new * tvm::dot(robot.q(), 1)->value();
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

  // Get normal of hitting plane
  // auto normal = robot_.frame(frame_->name()).position().rotation().col(axis_);

  // // Generate necessary matrices
  Eigen::MatrixXd j_m_uninverted = /*full_world_frame_jacobian.transpose() * */(full_world_frame_jacobian * M.inverse() * full_world_frame_jacobian.transpose());
  assert(std::abs(j_m_uninverted.determinant()) > 1e-5 && "Mass matrix is singular");
  Eigen::MatrixXd j_m = full_world_frame_jacobian.transpose() * j_m_uninverted.inverse();
  // Eigen::MatrixXd P_n_sub = normal_ * normal_.transpose();
  // Eigen::Matrix<double, 6, 6> P_n = Eigen::Matrix<double, 6, 6>::Zero();
  // P_n.block<3,3>(0,0) = P_n_sub;

  J_ddq = j_m * P_n * full_world_frame_jacobian;              // multiplies ddq variable
  splitJacobian(J_ddq, q_ddot_var_);
}
} // namespace mc_tvm


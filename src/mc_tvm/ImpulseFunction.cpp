/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tvm/ImpulseFunction.h>

#include <mc_tvm/Robot.h>
#include <mc_tvm/RobotFrame.h>

namespace mc_tvm
{

ImpulseFunction::ImpulseFunction(const mc_rbdyn::Robot & robot, const mc_rbdyn::RobotFrame & frame, double delta_t, double c_res, double limit_multiplier, int axis)
: tvm::function::abstract::LinearFunction(robot.mb().nrDof()), robot_(robot), frame_(frame), delta_t_(delta_t), c_res_(c_res), limit_multiplier_(limit_multiplier), axis_(axis)
  , jac_(frame.tvm_frame().rbdJacobian())
{
  assert(frame_->robot().robotIndex() == robot_.robotIndex() && "ImpulseFunction frame must belong to the robot");
  assert(axis_ >= 0 && axis_ < 3 && "Axis must be in [0, 2] as there are only three translational axis in the frame (1=x, 2=y, 3=z)");
  registerUpdates(Update::B, &ImpulseFunction::updateb);
  registerUpdates(Update::Jacobian, &ImpulseFunction::updateJacobian);
  addOutputDependency<ImpulseFunction>(Output::B, Update::B);
  addOutputDependency<ImpulseFunction>(Output::Jacobian, Update::Jacobian);
  auto & tvm_robot = robot.tvmRobot();
  addInputDependency<ImpulseFunction>(Update::Jacobian, tvm_robot, Robot::Output::H);
  addInputDependency<ImpulseFunction>(Update::B, tvm_robot, Robot::Output::C);
  addVariable(tvm::dot(tvm_robot.q(), 2), true);
  addVariable(tvm::dot(tvm_robot.q(), 1), true);
  velocity_.setZero();
}


void ImpulseFunction::updateb()
{
  b_ = - robot_.tvmRobot().limits().tl * limit_multiplier_ / (c_res_+1);
}

void ImpulseFunction::updateJacobian()
{
  const auto & robot = robot_.tvmRobot();
  auto M = robot.H();

  // Get jacobian of the robot frame
  const rbd::MultiBody robot_mb = robot_.mb();
  const rbd::MultiBodyConfig mbc = robot_.mbc();
  const auto & world_frame_jacobian = jac_.jacobian(robot_mb, mbc);
  Eigen::MatrixXd full_world_frame_jacobian(6, robot_.mb().nrDof());
  jac_.fullJacobian(robot_mb, world_frame_jacobian, full_world_frame_jacobian);

  // Get normal of hitting plane
  auto normal = robot_.frame(frame_->name()).position().rotation().col(axis_);

  // Generate necessary matrices
  auto j_m = full_world_frame_jacobian.transpose() * (full_world_frame_jacobian * M.inverse() * full_world_frame_jacobian.transpose()).inverse();
  auto P_n = normal * normal.transpose();

  splitJacobian(j_m*P_n*full_world_frame_jacobian, robot.alphaD());
  splitJacobian(j_m*P_n*full_world_frame_jacobian/delta_t_, robot.alpha());

  splitJacobian(robot.H(), robot.alphaD());
}
} // namespace mc_tvm


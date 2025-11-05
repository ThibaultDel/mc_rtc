/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tvm/DynamicFunction.h>
#include <mc_tvm/ForceFunction.h>

#include <mc_rbdyn/Robot.h>
#include <mc_tvm/Robot.h>

namespace mc_tvm
{

ForceFunction::ForceFunction(const mc_rbdyn::Robot & robot,
                             const Eigen::MatrixXd jTransposePseudoInverse,
                             bool compensateExternalForces)
: tvm::function::abstract::LinearFunction(robot.mb().nrDof()), robot_(robot),
  jTransposePseudoInverse_(jTransposePseudoInverse), compensateExternalForces_(compensateExternalForces)
{
  target_ = Eigen::VectorXd::Zero(jTransposePseudoInverse.rows());
  registerUpdates(Update::B, &ForceFunction::updateb);
  registerUpdates(Update::Jacobian, &ForceFunction::updateJacobian);
  addOutputDependency<ForceFunction>(Output::B, Update::B);
  addOutputDependency<ForceFunction>(Output::Jacobian, Update::Jacobian);
  auto & tvm_robot = robot.tvmRobot();
  addInputDependency<ForceFunction>(Update::Jacobian, tvm_robot, Robot::Output::H);
  addInputDependency<ForceFunction>(Update::B, tvm_robot, Robot::Output::C);
  addInputDependency<ForceFunction>(Update::B, tvm_robot, Robot::Output::ExternalForces);
  addVariable(tvm::dot(tvm_robot.q(), 2), true);
}

void ForceFunction::updateb() // Ax + b = 0
{
  Eigen::VectorXd D = robot_.tvmRobot().C();
  assert(target_.size() == D.size());
  assert((target_.array() == target_.array()).all()); // not NaN

  if(!jTransposePseudoInverse_.allFinite()) mc_rtc::log::error("[ForceFunction] NaN in jTranspose_ before solve");

  if(!D.allFinite()) mc_rtc::log::error("[ForceFunction] NaN in dynamic bias C() before solve");

  if(!compensateExternalForces_)
  {
    Eigen::VectorXd extForces = robot_.tvmRobot().tauExternal();
    D -= extForces;
  }
  b_ = jTransposePseudoInverse_ * D - target_;
}

void ForceFunction::updateJacobian()
{
  const auto & robot = robot_.tvmRobot();
  if(!jTransposePseudoInverse_.allFinite()) mc_rtc::log::error("[ForceFunction] NaN in jTranspose_ before solve");
  splitJacobian(jTransposePseudoInverse_ * robot.H(), robot.alphaD());
}

} // namespace mc_tvm

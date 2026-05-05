/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tvm/DynamicFunction.h>
#include <mc_tvm/TorqueFunction.h>

#include <mc_rbdyn/Robot.h>
#include <mc_tvm/Robot.h>
#include <RBDyn/MultiBodyConfig.h>
#include <SpaceVecAlg/EigenTypedef.h>

namespace mc_tvm
{

TorqueFunction::TorqueFunction(const mc_rbdyn::Robot & robot, bool compensateExternalForces, bool compensateGravity)
: tvm::function::abstract::LinearFunction(robot.mb().nrDof()), robot_(robot),
  compensateExternalForces_(compensateExternalForces), compensateGravity_(compensateGravity),
  j0_(robot_.mb().joint(0).type() == rbd::Joint::Free ? 1 : 0)
{
  registerUpdates(Update::B, &TorqueFunction::updateb);
  addOutputDependency<TorqueFunction>(Output::B, Update::B);
  auto & tvm_robot = robot.tvmRobot();
  addInputDependency<TorqueFunction>(Update::B, tvm_robot, Robot::Output::tau);
  addVariable(tvm::dot(tvm_robot.q(), 2), true);
  addVariable(tvm_robot.tau(), true); // x
  jacobian_[tvm_robot.tau().get()] = Eigen::MatrixXd::Identity(robot_.mb().nrDof(), robot_.mb().nrDof());
  jacobian_[tvm_robot.tau().get()].properties(tvm::internal::MatrixProperties::IDENTITY); // A
  reset();
}

void TorqueFunction::updateb() // Ax + b = 0
{
  b_ = -torque_;

  torque_extForces_ = robot_.tvmRobot().tauExternal();
  torque_gravity_ = robot_.tvmRobot().C();
  if(compensateExternalForces_) { b_ += torque_extForces_; }
  if(compensateGravity_) { b_ -= torque_gravity_; }
}

void TorqueFunction::reset()
{
  torque_mc_rtc_ = robot_.mbc().jointTorque;
  torque_ = rbd::sDofToVector(robot_.mb(), torque_mc_rtc_);
  torque_extForces_ = robot_.tvmRobot().tauExternal();
  torque_gravity_ = robot_.tvmRobot().C();
}

void TorqueFunction::torque(const std::string & j, const std::vector<double> & tau)
{
  if(!robot_.hasJoint(j))
  {
    mc_rtc::log::error("[TorqueFunction] No joint named {} in {}", j, robot_.name());
    return;
  }
  auto jIndex = static_cast<size_t>(robot_.mb().jointIndexByName(j));
  if(torque_mc_rtc_[jIndex].size() != tau.size())
  {
    mc_rtc::log::error("[TorqueFunction] Wrong size for input target on joint {}, excepted {} got {}", j,
                       torque_mc_rtc_[jIndex].size(), tau.size());
    return;
  }
  torque_mc_rtc_[static_cast<size_t>(jIndex)] = tau;
  torque_ = rbd::sDofToVector(robot_.mb(), torque_mc_rtc_);
}

bool TorqueFunction::isValidTorque(const std::vector<std::vector<double>> & ref,
                                   const std::vector<std::vector<double>> & in)
{
  if(ref.size() != in.size()) { return false; }
  for(size_t i = 0; i < ref.size(); ++i)
  {
    if(ref[i].size() != in[i].size()) { return false; }
  }
  return true;
}

void TorqueFunction::torque(const std::vector<std::vector<double>> & tau)
{
  if(!isValidTorque(torque_mc_rtc_, tau))
  {
    mc_rtc::log::error("[TorqueFunction] Invalid torque provided for {}", robot_.name());
    return;
  }
  torque_mc_rtc_ = tau;
  torque_ = rbd::sDofToVector(robot_.mb(), torque_mc_rtc_);
}

} // namespace mc_tvm

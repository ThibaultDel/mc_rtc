/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tvm/DynamicFunction.h>
#include <mc_tvm/TorqueFunction.h>

#include <mc_rbdyn/Robot.h>
#include <mc_tvm/Robot.h>

namespace mc_tvm
{

TorqueFunction::TorqueFunction(const mc_rbdyn::Robot & robot, bool compensateExternalForces)
: tvm::function::abstract::LinearFunction(robot.mb().nrDof()), robot_(robot),
  compensateExternalForces_(compensateExternalForces), j0_(robot_.mb().joint(0).type() == rbd::Joint::Free ? 1 : 0)
{
  registerUpdates(Update::B, &TorqueFunction::updateb);
  registerUpdates(Update::Jacobian, &TorqueFunction::updateJacobian);
  addOutputDependency<TorqueFunction>(Output::B, Update::B);
  addOutputDependency<TorqueFunction>(Output::Jacobian, Update::Jacobian);
  auto & tvm_robot = robot.tvmRobot();
  addInputDependency<TorqueFunction>(Update::Jacobian, tvm_robot, Robot::Output::H);
  addInputDependency<TorqueFunction>(Update::B, tvm_robot, Robot::Output::C);
  addInputDependency<TorqueFunction>(Update::B, tvm_robot, Robot::Output::ExternalForces);
  addVariable(tvm::dot(tvm_robot.q(), 2), true);
  velocity_.setZero();

  reset();
}

void TorqueFunction::updateb() // Ax + b = 0
{
  b_ = robot_.tvmRobot().C() - torque_;
  if(!compensateExternalForces_)
  {
    Eigen::VectorXd extForces = robot_.tvmRobot().tauExternal();
    b_ -= extForces;
  }
}

void TorqueFunction::updateJacobian()
{
  const auto & robot = robot_.tvmRobot();
  splitJacobian(robot.H(), robot.alphaD());
}

void TorqueFunction::reset()
{
  torque_ = robot_.tvmRobot().tau()->value();
  torque_mc_rtc_ = robot_.mbc().jointTorque;
  mcrtcTorqueToEigen();
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
  mcrtcTorqueToEigen();
}

namespace
{
bool isValidTorque(const std::vector<std::vector<double>> & ref, const std::vector<std::vector<double>> & in)
{
  if(ref.size() != in.size()) { return false; }
  for(size_t i = 0; i < ref.size(); ++i)
  {
    if(ref[i].size() != in[i].size()) { return false; }
  }
  return true;
}
} // namespace

void TorqueFunction::torque(const std::vector<std::vector<double>> & tau)
{
  if(!isValidTorque(torque_mc_rtc_, tau))
  {
    mc_rtc::log::error("[TorqueFunction] Invalid torque provided for {}", robot_.name());
    return;
  }
  torque_mc_rtc_ = tau;
  mcrtcTorqueToEigen();
}

void TorqueFunction::eigenToMCrtcTorque()
{
  int pos = 0;
  if(robot_.mb().nrJoints() > 0 && robot_.mb().joint(0).type() == rbd::Joint::Free)
  {
    pos = 6; // Skip the floating base joints
  }
  for(int jI = j0_; jI < robot_.mb().nrJoints(); ++jI)
  {
    auto jIdx = static_cast<size_t>(jI);
    const auto & j = robot_.mb().joint(jI);
    if(j.dof() == 1) // prismatic or revolute
    {
      torque_mc_rtc_[jIdx][0] = torque_[pos];
      pos++;
    }
  }
}

void TorqueFunction::mcrtcTorqueToEigen()
{
  int pos = 0;
  if(robot_.mb().nrJoints() > 0 && robot_.mb().joint(0).type() == rbd::Joint::Free)
  {
    pos = 6; // Skip the floating base joints
  }
  for(int jI = j0_; jI < robot_.mb().nrJoints(); ++jI)
  {
    auto jIdx = static_cast<size_t>(jI);
    const auto & j = robot_.mb().joint(jI);
    if(j.dof() == 1) // prismatic or revolute
    {
      torque_[pos] = torque_mc_rtc_[jIdx][0];
      pos++;
    }
  }
}

} // namespace mc_tvm

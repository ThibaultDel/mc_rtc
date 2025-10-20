/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tvm/DynamicFunction.h>
#include <mc_tvm/ForceFunction.h>

#include <mc_rbdyn/Robot.h>
#include <mc_tvm/Robot.h>
#include "mc_rbdyn/RobotFrame.h"

namespace mc_tvm
{

ForceFunction::ForceFunction(const mc_rbdyn::Robot & robot,
                             const mc_rbdyn::RobotFrame & frame,
                             bool compensateExternalForces)
: tvm::function::abstract::LinearFunction(robot.mb().nrDof()), robot_(robot), frame_(frame),
  compensateExternalForces_(compensateExternalForces), j0_(robot_.mb().joint(0).type() == rbd::Joint::Free ? 1 : 0)
{

  rbd::Jacobian jac = rbd::Jacobian(robot_.mb(), frame_.name());
  jTranspose_ = jac.jacobian(robot_.mb(), robot_.mbc());
  jTranspose_.transposeInPlace();

  registerUpdates(Update::B, &ForceFunction::updateb);
  registerUpdates(Update::Jacobian, &ForceFunction::updateJacobian);
  addOutputDependency<ForceFunction>(Output::B, Update::B);
  addOutputDependency<ForceFunction>(Output::Jacobian, Update::Jacobian);
  auto & tvm_robot = robot.tvmRobot();
  addInputDependency<ForceFunction>(Update::Jacobian, tvm_robot, Robot::Output::H);
  addInputDependency<ForceFunction>(Update::B, tvm_robot, Robot::Output::C);
  addInputDependency<ForceFunction>(Update::B, tvm_robot, Robot::Output::ExternalForces);
  addVariable(tvm::dot(tvm_robot.q(), 2), true);
  velocity_.setZero();

  reset();
}

void ForceFunction::updateb() // Ax + b = 0
{
  Eigen::VectorXd D = robot_.tvmRobot().C();
  if(!compensateExternalForces_)
  {
    Eigen::VectorXd extForces = robot_.tvmRobot().tauExternal();
    D -= extForces;
  }
  b_ = jTranspose_.llt().solve(D) - force_;
}

void ForceFunction::updateJacobian()
{
  const auto & robot = robot_.tvmRobot();
  splitJacobian(jTranspose_.llt().solve(robot.H()), robot.alphaD());
}

void ForceFunction::reset()
{
  force_ = jTranspose_.llt().solve(robot_.tvmRobot().tau()->value());
  // force_mc_rtc_ = robot_.mbc().jointTorque;
  // mcrtcforceToEigen();
}

// void ForceFunction::force(const std::string & j, const std::vector<double> & tau)
// {
//   if(!robot_.hasJoint(j))
//   {
//     mc_rtc::log::error("[ForceFunction] No joint named {} in {}", j, robot_.name());
//     return;
//   }
//   auto jIndex = static_cast<size_t>(robot_.mb().jointIndexByName(j));
//   if(force_mc_rtc_[jIndex].size() != tau.size())
//   {
//     mc_rtc::log::error("[ForceFunction] Wrong size for input target on joint {}, excepted {} got {}", j,
//                        force_mc_rtc_[jIndex].size(), tau.size());
//     return;
//   }
//   force_mc_rtc_[static_cast<size_t>(jIndex)] = tau;
//   mcrtcforceToEigen();
// }

// namespace
// {
// bool isValidforce(const std::vector<std::vector<double>> & ref, const std::vector<std::vector<double>> & in)
// {
//   if(ref.size() != in.size()) { return false; }
//   for(size_t i = 0; i < ref.size(); ++i)
//   {
//     if(ref[i].size() != in[i].size()) { return false; }
//   }
//   return true;
// }
// } // namespace

// void ForceFunction::force(const std::vector<std::vector<double>> & tau)
// {
//   if(!isValidforce(force_mc_rtc_, tau))
//   {
//     mc_rtc::log::error("[ForceFunction] Invalid force provided for {}", robot_.name());
//     return;
//   }
//   force_mc_rtc_ = tau;
//   mcrtcforceToEigen();
// }

// void ForceFunction::eigenToMCrtcforce()
// {
//   int pos = 0;
//   if(robot_.mb().nrJoints() > 0 && robot_.mb().joint(0).type() == rbd::Joint::Free)
//   {
//     pos = 6; // Skip the floating base joints
//   }
//   for(int jI = j0_; jI < robot_.mb().nrJoints(); ++jI)
//   {
//     auto jIdx = static_cast<size_t>(jI);
//     const auto & j = robot_.mb().joint(jI);
//     if(j.dof() == 1) // prismatic or revolute
//     {
//       force_mc_rtc_[jIdx][0] = force_[pos];
//       pos++;
//     }
//   }
// }

// void ForceFunction::mcrtcforceToEigen()
// {
//   int pos = 0;
//   if(robot_.mb().nrJoints() > 0 && robot_.mb().joint(0).type() == rbd::Joint::Free)
//   {
//     pos = 6; // Skip the floating base joints
//   }
//   for(int jI = j0_; jI < robot_.mb().nrJoints(); ++jI)
//   {
//     auto jIdx = static_cast<size_t>(jI);
//     const auto & j = robot_.mb().joint(jI);
//     if(j.dof() == 1) // prismatic or revolute
//     {
//       force_[pos] = force_mc_rtc_[jIdx][0];
//       pos++;
//     }
//   }
// }

} // namespace mc_tvm

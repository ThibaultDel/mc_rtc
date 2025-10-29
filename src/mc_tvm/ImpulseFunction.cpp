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
  // jacobian_[tvm_robot.tau().get()] = -Eigen::MatrixXd::Identity(robot_.mb().nrDof(), robot_.mb().nrDof());
  // jacobian_[tvm_robot.tau().get()].properties(tvm::internal::MatrixProperties::MINUS_IDENTITY);
  velocity_.setZero();
}

// ImpulseFunction::ForceContact::ForceContact(const mc_rbdyn::RobotFrame & frame,
//                                             std::vector<sva::PTransformd> points,
//                                             double dir)
// : frame_(frame), points_(std::move(points)), dir_(dir), jac_(frame.tvm_frame().rbdJacobian()),
//   blocks_(jac_.compactPath(frame.robot().mb())), force_jac_(6, jac_.dof()), full_jac_(6, frame.robot().mb().nrDof())
// {
//   for(size_t i = 0; i < points_.size(); ++i) { forces_.add(tvm::Space(3).createVariable("force" + std::to_string(i))); }
//   forces_.setZero();
// }
//
// void ImpulseFunction::ForceContact::updateJacobians(ImpulseFunction & parent)
// {
//   const auto & robot = frame_->robot();
//   const auto & bodyJac = jac_.bodyJacobian(robot.mb(), robot.mbc());
//   for(int i = 0; i < forces_.numberOfVariables(); ++i)
//   {
//     const auto & force = forces_[i];
//     const auto & point = points_[static_cast<size_t>(i)];
//     jac_.translateBodyJacobian(bodyJac, robot.mbc(), point.translation(), force_jac_);
//     full_jac_.setZero();
//     jac_.addFullJacobian(blocks_, force_jac_, full_jac_);
//     parent.jacobian_[force.get()].noalias() = -dir_ * full_jac_.block(3, 0, 3, robot.mb().nrDof()).transpose();
//   }
// }
//
// sva::ForceVecd ImpulseFunction::ForceContact::force() const
// {
//   sva::ForceVecd ret = sva::ForceVecd::Zero();
//   for(int i = 0; i < forces_.numberOfVariables(); ++i)
//   {
//     const auto & force = forces_[i];
//     const auto & point = points_[static_cast<size_t>(i)];
//     ret += point.transMul(sva::ForceVecd(Eigen::Vector3d::Zero(), force->value()));
//   }
//   return ret;
// }
//
// const tvm::VariableVector & ImpulseFunction::addContact(const mc_rbdyn::RobotFrame & frame,
//                                                         std::vector<sva::PTransformd> points,
//                                                         double dir)
// {
//   if(frame.robot().name() != robot_.name())
//   {
//     mc_rtc::log::error_and_throw<std::runtime_error>(
//         "Attempted to add a contact for {} to dynamic function belonging to {}", frame.robot().name(), robot_.name());
//   }
//   auto & fc = contacts_.emplace_back(frame, std::move(points), dir);
//   for(const auto & var : fc.forces_) { addVariable(var, true); }
//   addInputDependency<ImpulseFunction>(Update::Jacobian, frame.tvm_frame(), mc_tvm::RobotFrame::Output::Jacobian);
//   return fc.forces_;
// }
//
// void ImpulseFunction::removeContact(const mc_rbdyn::RobotFrame & frame)
// {
//   auto it = findContact(frame);
//   if(it != contacts_.end())
//   {
//     for(const auto & var : it->forces_) { removeVariable(var); }
//     contacts_.erase(it);
//   }
// }
//
// sva::ForceVecd ImpulseFunction::contactForce(const mc_rbdyn::RobotFrame & frame) const
// {
//   auto it = findContact(frame);
//   if(it != contacts_.end()) { return (*it).force(); }
//   else
//   {
//     mc_rtc::log::error("No contact at frame {} in dynamic function for {}", frame.name(), robot_.name());
//     return sva::ForceVecd(Eigen::Vector6d::Zero());
//   }
// }

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
  // for(auto & c : contacts_) { c.updateJacobians(*this); }
}

// auto ImpulseFunction::findContact(const mc_rbdyn::RobotFrame & frame) const -> std::vector<ForceContact>::const_iterator
// {
//   return std::find_if(contacts_.begin(), contacts_.end(), [&](const auto & c) { return c.frame_.get() == &frame; });
// }

} // namespace mc_tvm


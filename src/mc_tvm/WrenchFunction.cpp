/*
 * Copyright 2015-2026 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tvm/DynamicFunction.h>
#include <mc_tvm/WrenchFunction.h>

#include <mc_tvm/Robot.h>
#include <mc_tvm/RobotFrame.h>

#include <mc_rbdyn/Robot.h>
#include <mc_rbdyn/RobotFrame.h>
#include <SpaceVecAlg/EigenTypedef.h>

#include <Eigen/Cholesky>
#include <Eigen/src/Core/Matrix.h>

namespace mc_tvm
{

WrenchFunction::WrenchFunction(const mc_rbdyn::RobotFrame & frame)
: tvm::function::abstract::LinearFunction(6), frame_(frame), tvm_frame_(frame.tvm_frame()), robot_(frame.robot()),
  tvm_robot_(robot_.tvmRobot()), frameJac_(tvm_frame_.rbdJacobian()), shortJacMat_(6, frameJac_.dof()),
  jacMat_(6, robot_.mb().nrDof()), dynamicJacMatTranspose_(6, robot_.mb().nrDof()),
  cartesianInertiaMat_(Eigen::Matrix6d::Zero()), wrench_(sva::ForceVecd::Zero())
{
  registerUpdates(Update::B, &WrenchFunction::updateb);
  registerUpdates(Update::Jacobian, &WrenchFunction::updateJacobian);
  addOutputDependency<WrenchFunction>(Output::B, Update::B);
  addOutputDependency<WrenchFunction>(Output::Jacobian, Update::Jacobian);

  auto & robot = frame_->robot();
  auto & tvm_robot = robot.tvmRobot();

  addInputDependency<WrenchFunction>(Update::Jacobian, tvm_robot, Robot::Output::H);
  addInternalDependency<WrenchFunction>(Update::Jacobian, Update::B);
  // addInputDependency<WrenchFunction>(Update::B, tvm_frame_, mc_tvm::RobotFrame::Output::Jacobian);
  addInputDependency<WrenchFunction>(Update::Jacobian, tvm_frame_, mc_tvm::RobotFrame::Output::Jacobian);
  addInputDependency<WrenchFunction>(Update::B, tvm_robot, Robot::Output::C);
  addInputDependency<WrenchFunction>(Update::B, tvm_robot, Robot::Output::tau);
  addVariable(tvm::dot(tvm_robot.q(), 2), true);
  addVariable(tvm_robot.tau(), true);
  jacobian_[tvm_robot.tau().get()] = Eigen::MatrixXd::Identity(robot.mb().nrDof(), robot.mb().nrDof());
  jacobian_[tvm_robot.tau().get()].properties(tvm::internal::MatrixProperties::IDENTITY);
  velocity_.setZero();
  reset();
}

void WrenchFunction::reset()
{
  wrench_ = currentWrench();
}

sva::ForceVecd WrenchFunction::currentWrench()
{
  computeDynamicJacobian();
  auto tau = robot_.mbc().jointTorque;
  Eigen::VectorXd tau_vec = rbd::sDofToVector(robot_.mb(), tau);
  Eigen::VectorXd wrenchVec_ = dynamicJacobianTranspose() * tau_vec;
  return sva::ForceVecd(wrenchVec_.head<3>(), wrenchVec_.tail<3>());
}

void WrenchFunction::computeDynamicJacobian()
{
  shortJacMat_ = frameJac_.jacobian(robot_.mb(), robot_.mbc(), tvm_frame_.position());
  frameJac_.fullJacobian(robot_.mb(), shortJacMat_, jacMat_);

  rbd::ForwardDynamics fd(robot_.mb());
  fd.computeH(robot_.mb(), robot_.mbc());

  const Eigen::MatrixXd & H = fd.H();

  // 1. Factorize H (SPD)
  Eigen::LDLT<Eigen::MatrixXd> H_ldlt(H);

  // 2. Compute M^{-1} J^T
  Eigen::MatrixXd MinvJt = H_ldlt.solve(jacMat_.transpose());

  // 3. Compute Cartesian Inertia Λ = (J M^{-1} J^T)^{-1} using LDLT
  Eigen::Matrix6d JMinvJT = jacMat_ * MinvJt;
  cartesianInertiaMat_ = JMinvJT.ldlt().solve(Eigen::Matrix6d::Identity());

  // 4. Compute (J^#)^T = Λ J M^{-1}
  Eigen::MatrixXd JMInv = MinvJt.transpose(); // = J M^{-1}
  dynamicJacMatTranspose_ = cartesianInertiaMat_ * JMInv;
}

void WrenchFunction::updateJacobian()
{
  computeDynamicJacobian();
  splitJacobian(dynamicJacobianTranspose(), tvm_robot_.tau());
}

void WrenchFunction::updateb() // Ax + b = 0
{
  b_ = -wrench_.vector();
}

} // namespace mc_tvm

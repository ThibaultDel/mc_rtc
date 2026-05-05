/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_tvm/api.h>
#include <mc_tvm/fwd.h>

#include <mc_rbdyn/fwd.h>

#include <tvm/function/abstract/LinearFunction.h>

#include <RBDyn/Jacobian.h>

#include <SpaceVecAlg/SpaceVecAlg>

namespace mc_tvm
{

/** This class implements a wrench function for a given frame */
class MC_TVM_DLLAPI WrenchFunction : public tvm::function::abstract::LinearFunction
{
public:
  using Output = tvm::function::abstract::LinearFunction::Output;
  DISABLE_OUTPUTS(Output::JDot)
  SET_UPDATES(WrenchFunction, B, Jacobian)

  /** Constructor
   *
   * Set the objective to the current frame wrench
   *
   */
  WrenchFunction(const mc_rbdyn::RobotFrame & frame);

  /** Set the target wrench to the current frame wrench */
  void reset();

  /** Get the current wrench estimated */
  sva::ForceVecd currentWrench();

  /** Get the current objective */
  sva::ForceVecd targetWrench() const { return wrench_; }

  /** Set the objective */
  void targetWrench(const sva::ForceVecd & wrench) noexcept { wrench_ = wrench; }

  /** Get the frame */
  const mc_rbdyn::RobotFrame & frame() const noexcept { return *frame_; }

  /** Get the current Transpose Dynamic Jacobian */
  const Eigen::MatrixXd & dynamicJacobianTranspose() const { return dynamicJacMatTranspose_; }

  const Eigen::Matrix6d & cartesianInertia() const { return cartesianInertiaMat_; }

protected:
  void updateb();
  void updateJacobian();
  void computeDynamicJacobian();

  mc_rbdyn::ConstRobotFramePtr frame_;
  mc_tvm::RobotFrame & tvm_frame_;
  const mc_rbdyn::Robot & robot_;
  const mc_tvm::Robot & tvm_robot_;

  /** Computation intermediate */
  rbd::Jacobian frameJac_;
  Eigen::MatrixXd shortJacMat_;
  Eigen::MatrixXd jacMat_;
  Eigen::MatrixXd dynamicJacMatTranspose_; // (J^#)^T = Λ J M^{-1}
  Eigen::Matrix6d cartesianInertiaMat_; // Λ = (J M^{-1} J^T)^{-1}

  /** Target */
  sva::ForceVecd wrench_;
};

} // namespace mc_tvm

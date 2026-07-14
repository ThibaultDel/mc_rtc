/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_tvm/api.h>

#include <mc_rbdyn/fwd.h>

#include <tvm/function/abstract/LinearFunction.h>

#include <RBDyn/Jacobian.h>

#include <SpaceVecAlg/EigenTypedef.h>
<<<<<<< HEAD
#include <SpaceVecAlg/ForceVec.h>
#include <SpaceVecAlg/SpaceVecAlg>

// #include <mc_rbdyn/VirtualforceSensor.h>
=======
>>>>>>> bastien/hrp5p

namespace mc_tvm
{

/** This class implements a force function for a given robot */
class MC_TVM_DLLAPI ForceFunction : public tvm::function::abstract::LinearFunction
{
public:
  using Output = tvm::function::abstract::LinearFunction::Output;
  DISABLE_OUTPUTS(Output::JDot)
  SET_UPDATES(ForceFunction, Jacobian, B)

  /** Constructor
   *
   * Set the objective to the current force of robot
   *
   */
  ForceFunction(const mc_rbdyn::Robot & robot,
                Eigen::MatrixXd jTransposePseudoInverse,
                bool compensateExternalForces = false);

  void forceTarget(const Eigen::VectorXd & tf) { target_ = tf; }
  Eigen::MatrixXd getJacobianT() const { return jTransposePseudoInverse_; }
  void setJacobianTPseudoInverse(const Eigen::MatrixXd & jTransposePseudoInverse)
  {
    jTransposePseudoInverse_ = jTransposePseudoInverse;
  }

  /** Access the full target force */
  Eigen::VectorXd forceTarget() const { return target_; }

  void compensateExternalForces(bool compensate) { compensateExternalForces_ = compensate; }

  bool isCompensatingExternalForces() const { return compensateExternalForces_; }

protected:
  // void updateValue();
  void updateb();
  void updateJacobian();

  const mc_rbdyn::Robot & robot_;
  Eigen::MatrixXd jTransposePseudoInverse_;

  bool compensateExternalForces_;

  Eigen::VectorXd target_;
};

} // namespace mc_tvm

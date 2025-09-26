/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_tvm/api.h>

#include <mc_rbdyn/fwd.h>

#include "mc_rbdyn/RobotFrame.h"
#include <tvm/function/abstract/LinearFunction.h>

#include <RBDyn/Jacobian.h>

#include <SpaceVecAlg/SpaceVecAlg>

// #include <mc_rbdyn/VirtualforceSensor.h>

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
                const mc_rbdyn::RobotFrame & frame,
                bool compensateExternalForces = false);

  /** Set the target force to the current robot's force */
  void reset();

  /** Set the target for a given joint
   *
   *  \param j Joint name
   *
   *  \param tau Target configuration
   *
   */
  void force(const std::string & j, const std::vector<double> & tau);

  /** Set the fully body force */
  void force(const std::vector<std::vector<double>> & tau);

  /** Access the full target force */
  const std::vector<std::vector<double>> & force() const noexcept { return force_mc_rtc_; }

  void compensateExternalForces(bool compensate) { compensateExternalForces_ = compensate; }

  bool isCompensatingExternalForces() const { return compensateExternalForces_; }

protected:
  // void updateValue();
  void updateb();
  void updateJacobian();

  const mc_rbdyn::Robot & robot_;
  /** Frame */
  const mc_rbdyn::RobotFrame & frame_;
  Eigen::MatrixXd jTranspose_;

  bool compensateExternalForces_;

  void eigenToMCrtcforce();
  void mcrtcforceToEigen();

  /** Target */
  Eigen::VectorXd force_;
  std::vector<std::vector<double>> force_mc_rtc_;
  /** Starting joint */
  int j0_;
};

} // namespace mc_tvm

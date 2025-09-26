/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_tvm/api.h>

#include <mc_rbdyn/fwd.h>

#include <tvm/function/abstract/LinearFunction.h>

#include <RBDyn/Jacobian.h>

#include <SpaceVecAlg/SpaceVecAlg>

#include <mc_rbdyn/VirtualTorqueSensor.h>

namespace mc_tvm
{

/** This class implements a torque function for a given robot */
class MC_TVM_DLLAPI TorqueFunction : public tvm::function::abstract::LinearFunction
{
public:
  using Output = tvm::function::abstract::LinearFunction::Output;
  DISABLE_OUTPUTS(Output::JDot)
  SET_UPDATES(TorqueFunction, Jacobian, B)

  /** Constructor
   *
   * Set the objective to the current torque of robot
   *
   */
  TorqueFunction(const mc_rbdyn::Robot & robot, bool compensateExternalForces = false);

  /** Set the target torque to the current robot's torque */
  void reset();

  /** Set the target for a given joint
   *
   *  \param j Joint name
   *
   *  \param tau Target configuration
   *
   */
  void torque(const std::string & j, const std::vector<double> & tau);

  /** Set the fully body torque */
  void torque(const std::vector<std::vector<double>> & tau);

  /** Access the full target torque */
  const std::vector<std::vector<double>> & torque() const noexcept { return torque_mc_rtc_; }

  void compensateExternalForces(bool compensate) { compensateExternalForces_ = compensate; }

  bool isCompensatingExternalForces() const { return compensateExternalForces_; }

protected:
  // void updateValue();
  void updateb();
  void updateJacobian();

  const mc_rbdyn::Robot & robot_;
  bool compensateExternalForces_;

  void eigenToMCrtcTorque();
  void mcrtcTorqueToEigen();

  /** Target */
  Eigen::VectorXd torque_;
  std::vector<std::vector<double>> torque_mc_rtc_;

  /** Starting joint */
  int j0_;
};

} // namespace mc_tvm

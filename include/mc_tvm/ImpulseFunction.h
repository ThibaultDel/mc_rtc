/*
* Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_tvm/api.h>

#include <mc_rbdyn/fwd.h>

#include <tvm/function/abstract/LinearFunction.h>

#include <RBDyn/Jacobian.h>

#include <SpaceVecAlg/SpaceVecAlg>

namespace mc_tvm
{

/** Implement the equation of motion for a given robot.
 *
 * It can be given contacts that will be integrated into the equation of
 * motion (\see ImpulseFunction::addContact).
 *
 * It manages the force variables related to these contacts.
 *
 * Notably, it does not take care of enforcing Newton 3rd law of motion when
 * two actuated robots are in contact.
 *
 */
struct MC_TVM_DLLAPI ImpulseFunction : public tvm::function::abstract::LinearFunction
{
public:
  using Output = tvm::function::abstract::LinearFunction::Output;
  DISABLE_OUTPUTS(Output::JDot)
  SET_UPDATES(ImpulseFunction, Jacobian, B)

  /** Construct the equation of motion for a given robot */
  ImpulseFunction(const mc_rbdyn::Robot & robot, const mc_rbdyn::RobotFrame & frame, const Eigen::Vector3d normal , double lambda/*, int axis*/);

protected:
  void updateb();

  void updateJacobian();

  const mc_rbdyn::Robot & robot_;
  mc_rbdyn::ConstRobotFramePtr frame_;
  // const double delta_t_;
  const double lambda;
  // const double c_res_;
  // const double limit_multiplier_;
  // const int axis_; // limited between 0 and 2, 0 = x-axis, 1 = y-axis and 2 = z-axis
  const Eigen::Vector3d normal_;

  Eigen::MatrixXd J_ddq;
  Eigen::MatrixXd J_dq;

  rbd::Jacobian jac_;

  /* Persisting tvm variables to avoid temporaries/dangling pointers */
  tvm::VariablePtr q_ddot_var_;
  tvm::VariablePtr q_dot_var_;
  /* Persisting vectors to avoid any temporary-vector lifetime issues */
  tvm::VariableVector q_ddot_vars_;
  tvm::VariableVector q_dot_vars_;
};

using ImpulseFunctionPtr = std::shared_ptr<ImpulseFunction>;

} // namespace mc_tvm

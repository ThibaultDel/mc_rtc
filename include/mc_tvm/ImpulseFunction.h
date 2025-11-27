/*
* Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_tvm/api.h>

#include <mc_rbdyn/fwd.h>

#include <tvm/function/abstract/LinearFunction.h>

#include <RBDyn/Jacobian.h>

#include <RBDyn/Coriolis.h>

#include <SpaceVecAlg/SpaceVecAlg>

#include <mc_tvm/Robot.h>

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
  ImpulseFunction(const mc_rbdyn::Robot & robot, const mc_rbdyn::RobotFrame & frame, const Eigen::Vector3d normal , double lambda_high, double lambda_low, double c_res, double delta_t, const Eigen::VectorXd & limit);

  Eigen::VectorXd & EffectiveLambda(){ return lambda; }

  Eigen::VectorXd & ImpulsiveTorquePrediction(){ return tau_imp; }

  Eigen::VectorXd & ImpulsiveTorquePrediction2(){ return tau_imp2; }

  Eigen::VectorXd & ImpulsiveTorquePredictionConstraint(){ return tau_imp_const; }

  Eigen::VectorXd & ActualImpulsiveTorquePrediction(){ return tau_imp_act; }

  Eigen::VectorXd & ImpulsiveTorquePredictionDerivative(){ return tau_imp_deriv; }

  Eigen::VectorXd & ImpulsiveTorquePredictionDerivativeNum(){ return tau_imp_deriv_num; }

  Eigen::VectorXd & ImpulsiveTorquePredictionDerivative_term1(){ return tau_imp_deriv_term1; }

  Eigen::VectorXd & ImpulsiveTorquePredictionDerivative_term2(){ return tau_imp_deriv_term2; }

  Eigen::VectorXd & ImpulsiveTorquePredictionDerivative_term3(){ return tau_imp_deriv_term3; }

  Eigen::VectorXd & JointPos();

  Eigen::VectorXd & JointVel();

  Eigen::VectorXd & JointVels();

  Eigen::VectorXd & JointVelUsed();

  Eigen::VectorXd & JointAcc();

  Eigen::VectorXd & JointAccUsed();

  Eigen::VectorXd & JointAccNum();

protected:
  void updateb();

  void updateJacobian();

  const mc_rbdyn::Robot & robot_;
  mc_rbdyn::ConstRobotFramePtr frame_;
  const Eigen::Vector3d normal_;
  // const double lambda;
  const double lambda_high;
  const double lambda_low;
  const double c_res_;
  const double delta_t_;
  const Eigen::VectorXd & limit_;

  Eigen::MatrixXd J_ddq;
  Eigen::MatrixXd J_dq;

  Eigen::VectorXd q;
  Eigen::VectorXd alpha;
  Eigen::VectorXd alpha_s_;
  Eigen::VectorXd alpha_d;
  Eigen::VectorXd q_d;
  Eigen::VectorXd q_dd;

  Eigen::VectorXd lambda;

  Eigen::VectorXd tau_imp;
  Eigen::VectorXd tau_imp2;
  Eigen::VectorXd tau_imp_const;
  Eigen::VectorXd tau_imp_act;
  Eigen::VectorXd tau_imp_deriv;
  Eigen::VectorXd tau_imp_deriv_num;
  Eigen::VectorXd tau_imp_deriv_term1;
  Eigen::VectorXd tau_imp_deriv_term2;
  Eigen::VectorXd tau_imp_deriv_term3;

  Eigen::VectorXd last_joint_velocities_;
  Eigen::VectorXd num_qdd;

  // Implement a moving average filter on the acceleration error
  // std::vector<Eigen::VectorXd> acc_error_window_;
  // int acc_error_window_size_ = 3;

  rbd::Jacobian jac_;

  Eigen::Matrix<double, 6, 6> P_n;

  rbd::Coriolis coriolis_calculator_;

  /* Persisting tvm variables to avoid temporaries/dangling pointers */
  tvm::VariablePtr q_ddot_var_;

};

using ImpulseFunctionPtr = std::shared_ptr<ImpulseFunction>;

} // namespace mc_tvm

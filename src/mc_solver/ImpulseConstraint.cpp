/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_solver/ImpulseConstraint.h>

#include <mc_solver/TasksQPSolver.h>

#include <Tasks/Bounds.h>
#include <Tasks/QPConstr.h>

#include "TVMImpulseConstraint.h"

namespace mc_solver
{

// static mc_rtc::void_ptr initialize_tasks(const mc_rbdyn::Robots & robots, unsigned int robotIndex, double timeStep)
// {
//   // TODO implement a Tasks backend version of the ImpulseConstraint, look at KinematicsConstraint for an example
// }

mc_tvm::ImpulseFunctionPtr initialize_imp_cstr(const mc_rbdyn::Robot & robot,
                                     const mc_rbdyn::RobotFrame & frame,
                                     const Eigen::Vector3d normal,
                                     double lambda_high,
                                     double lambda_low,
                                     double c_res,
                                     double delta_t,
                                     const Eigen::VectorXd & limit_high,
                                     const Eigen::VectorXd & limit_low,
                                     bool enforce_high)
{
  return /*mc_rtc::make_void_ptr<mc_tvm::ImpulseFunctionPtr>*/(std::make_shared<mc_tvm::ImpulseFunction>(robot, frame, normal, lambda_high, lambda_low, c_res, delta_t,
                                                                                                      limit_high, limit_low, enforce_high));
}

TVMImpulseConstraint::TVMImpulseConstraint(const mc_rbdyn::Robot & robot,
                                                 const mc_rbdyn::RobotFrame & frame,
                                                 const Eigen::Vector3d normal,
                                                 double lambda_high,
                                                 double lambda_low,
                                                 double delta_t,
                                                 double c_res,
                                                 double limit_multiplier)
: robot_(robot), frame_(frame), lambda_high_(lambda_high), lambda_low_(lambda_low), delta_t_(delta_t), c_res_(c_res),
  limit_multiplier_(limit_multiplier),
  const_upper_limit_(robot_.tvmRobot().limits().tu * limit_multiplier_),
  const_lower_limit_(robot_.tvmRobot().limits().tl * limit_multiplier_),
  imp_constr_lower_(initialize_imp_cstr(robot, frame, normal, lambda_high_, lambda_low_, c_res, delta_t, const_upper_limit_, const_lower_limit_, false)),
  imp_constr_upper_(initialize_imp_cstr(robot, frame, normal, lambda_high_, lambda_low_, c_res, delta_t, const_upper_limit_, const_lower_limit_, true)),
  upper_limit_(robot_.tvmRobot().limits().tu * limit_multiplier_ * (delta_t_/(c_res_+1.f))),
  lower_limit_(robot_.tvmRobot().limits().tl * limit_multiplier_ * (delta_t_/(c_res_+1.f)))
{
  mc_rtc::log::info("Size of upper_limit_ is {} x {}", upper_limit_.rows(), upper_limit_.cols());
  mc_rtc::log::info("Size of lower_limit_ is {} x {}", lower_limit_.rows(), lower_limit_.cols());
  lower = Eigen::VectorXd::Zero(robot_.mb().nrDof());
  upper = Eigen::VectorXd::Zero(robot_.mb().nrDof());
}

void TVMImpulseConstraint::addToSolver(mc_solver::TVMQPSolver & solver)
{
      auto & problem = tvm_solver(solver).problem();
      mc_rtc::log::info("Size of the lower limit is {}", lower_limit_.size());
      mc_tvm::ImpulseFunctionPtr imp_fn_up = /**static_cast<mc_tvm::ImpulseFunctionPtr *>*/(imp_constr_upper_/*.get()*/);
      mc_tvm::ImpulseFunctionPtr imp_fn_low = /**static_cast<mc_tvm::ImpulseFunctionPtr *>(imp_constr_lower_.get())*/imp_constr_lower_;
      mc_rtc::log::info("Created the TVM functions");
      auto imp_up = problem.add( imp_fn_up <= 0.0/*robot_.tvmRobot().limits().tu * limit_multiplier_ * lambda_low_*/, tvm::task_dynamics::None(), {tvm::requirements::PriorityLevel(0)});
      auto imp_low = problem.add( imp_fn_low >= 0.0/*robot_.tvmRobot().limits().tl * limit_multiplier_ * lambda_low_*/, tvm::task_dynamics::None(), {tvm::requirements::PriorityLevel(0)});
      mc_rtc::log::info("made the problem additions");
      constraints_.push_back(imp_up);
      constraints_.push_back(imp_low);
      mc_rtc::log::info("Pushed back the constraints");
}

void TVMImpulseConstraint::removeFromSolver(mc_solver::TVMQPSolver & solver)
{
  for(auto & c : mimics_constraints_)
  {
    solver.problem().removeSubstitutionFor(*solver.problem().constraint(*c));
    solver.problem().remove(*c);
  }
  for(auto & c : constraints_) { solver.problem().remove(*c); }
  constraints_.clear();
  mimics_constraints_.clear();
}

static mc_rtc::void_ptr initialize_tvm(const mc_rbdyn::Robot & robot, const mc_rbdyn::RobotFrame & frame, const Eigen::Vector3d normal, double lambda_high, double lambda_low, double delta_t, double c_res, double limit_multiplier)
{
  return mc_rtc::make_void_ptr<TVMImpulseConstraint>(robot, frame, normal, lambda_high, lambda_low, delta_t, c_res, limit_multiplier);
}

static mc_rtc::void_ptr initialize(QPSolver::Backend backend,
                                   const mc_rbdyn::Robots & robots,
                                   unsigned int robotIndex,
                                   const mc_rbdyn::RobotFrame & frame,
                                   const Eigen::Vector3d normal,
                                   double lambda_high,
                                   double lambda_low,
                                   double delta_t,
                                   double c_res,
                                   double limit_multiplier)
{
  switch(backend)
  {
    case QPSolver::Backend::Tasks:
      mc_rtc::log::error("No implementation for the ImpulseConstraint with the Tasks backend");
      assert(false);
    case QPSolver::Backend::TVM:
      return initialize_tvm(robots.robot(robotIndex), frame, normal, lambda_high, lambda_low, delta_t, c_res, limit_multiplier);
    default:
      mc_rtc::log::error_and_throw("[ImpulseConstraint] Not implemented for solver backend: {}", backend);
  }
}

ImpulseConstraint::ImpulseConstraint(const mc_rbdyn::Robots & robots, unsigned int robotIndex, const mc_rbdyn::RobotFrame & frame, const Eigen::Vector3d normal, double lambda_high, double lambda_low, double delta_t, double c_res, double limit_multiplier, mc_rtc::Logger & logger)
: constraint_(initialize(backend_, robots, robotIndex, frame, normal, lambda_high, lambda_low, delta_t, c_res, limit_multiplier)), logger_(logger)
{
  add_logs();
}

void ImpulseConstraint::addToSolverImpl(mc_solver::QPSolver & solver)
{
  switch(backend_)
  {
    case QPSolver::Backend::Tasks:
      mc_rtc::log::error("No implementation for the ImpulseConstraint with the Tasks backend");
      assert(false);
    case QPSolver::Backend::TVM:
      static_cast<TVMImpulseConstraint *>(constraint_.get())->addToSolver(tvm_solver(solver));
      break;
    default:
      break;
  }
}

void ImpulseConstraint::removeFromSolverImpl(mc_solver::QPSolver & solver)
{
  switch(backend_)
  {
    case QPSolver::Backend::Tasks:
      mc_rtc::log::error("No implementation for the ImpulseConstraint with the Tasks backend");
      assert(false);
    case QPSolver::Backend::TVM:
      // logger_.removeLogEntry("ImpulseConstraintEvaluation");
      // logger_.removeLogEntry("ImpulseConstraintLowerlimit");
      // logger_.removeLogEntry("ImpulseConstraintUpperlimit");
      // logger_.removeLogEntry("PredictedImpulsiveTorque");
      static_cast<TVMImpulseConstraint *>(constraint_.get())->removeFromSolver(tvm_solver(solver));
      break;
    default:
      break;
  }
}

void ImpulseConstraint::add_logs()
{
  logger_.addLogEntry("ImpulseConstraint_Evaluation_lower", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunctionLow()->value();});

  logger_.addLogEntry("ImpulseConstraint_Evaluation_upper", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunctionHigh()->value();});

  logger_.addLogEntry("ImpulseConstraint_Elementwise_lambda_low", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunctionLow()->EffectiveLambda();});

  logger_.addLogEntry("ImpulseConstraint_Elementwise_lambda_high", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunctionHigh()->EffectiveLambda();});

  logger_.addLogEntry("ImpulseConstraint_ImpulsiveTorqueLowerlimit", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->LowerLimit();});

  logger_.addLogEntry("ImpulseConstraint_ImpulsiveTorqueUpperlimit", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->UpperLimit();});
  //
  logger_.addLogEntry("ImpulseConstraint_PredictedImpulsiveTorque", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunctionLow()->ActualImpulsiveTorquePrediction();});
  //
  // logger_.addLogEntry("PredictedImpulsiveTorque2", this, [&, this]()
  // {return static_cast<TVMImpulseConstraint *>(constraint_.get())->ImpulsiveTorqures2();});
  //
  // logger_.addLogEntry("PredictedImpulsiveTorqueIntegrated", this, [&, this]()
  // {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunctionLow()->ImpulsiveTorquePredictionConstraint();});
  //
  // logger_.addLogEntry("ActualPredictedImpulsiveTorque", this, [&, this]()
  // {return static_cast<TVMImpulseConstraint *>(constraint_.get())->ActualImpulsiveTorqures();});
  //
  logger_.addLogEntry("ImpulseConstraint_PredictedImpulsiveTorqueDerivativeExpected", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunctionLow()->ImpulsiveTorquePredictionDerivative();});
  //
  logger_.addLogEntry("ImpulseConstraint_PredictedImpulsiveTorqueDerivativeNumerical", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunctionLow()->ImpulsiveTorquePredictionDerivativeNum();});
  //
  // logger_.addLogEntry("PredictedImpulsiveTorqueDerivative_term1", this, [&, this]()
  // {return static_cast<TVMImpulseConstraint *>(constraint_.get())->ImpulsiveTorquresDerivative_term1();});
  //
  // logger_.addLogEntry("PredictedImpulsiveTorqueDerivative_term2", this, [&, this]()
  // {return static_cast<TVMImpulseConstraint *>(constraint_.get())->ImpulsiveTorquresDerivative_term2();});
  //
  // logger_.addLogEntry("PredictedImpulsiveTorqueDerivative_term3", this, [&, this]()
  // {return static_cast<TVMImpulseConstraint *>(constraint_.get())->ImpulsiveTorquresDerivative_term3();});
  //
  logger_.addLogEntry("ImpulseConstraint_RightSideUpperLimit", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunctionHigh()->RightSide();});

  logger_.addLogEntry("ImpulseConstraint_RightSideLowerLimit", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunctionLow()->RightSide();});

  // equal to qIn and qOut
  // logger_.addLogEntry("TVM_q", this, [&, this]()
  // {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunction()->JointPos();});

  // Does not give any value for some reason
  // logger_.addLogEntry("TVM_alpha", this, [&, this]()
  // {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunction()->JointVel();});

  // the value used to determine s in the passivity plugin, also in integral_of_reference_acceleration in MyPluginPassivityHumanoid
  // logger_.addLogEntry("TVM_alphas", this, [&, this]()
  // {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunctionLow()->JointVels();});

  // equal to alphaDOut
  logger_.addLogEntry("ImpulseConstraint_CommandedAcceleration", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunctionLow()->JointAcc();});

  // Does not give any value for some reason
  // logger_.addLogEntry("TVM_qd", this, [&, this]()
  // {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunctionLow()->JointVelUsed();});

  // This is alphaDOut
  // logger_.addLogEntry("TVM_qdd", this, [&, this]()
  // {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunction()->JointAccUsed();});

  // this is the numerical derivative of joint velocities AlphaIn
  logger_.addLogEntry("ImpulseConstraint_numericalAcceleration", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunctionLow()->JointAccNum();});
}


} // namespace mc_solver

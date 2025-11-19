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

mc_rtc::void_ptr initialize_imp_cstr(const mc_rbdyn::Robot & robot,
                                     const mc_rbdyn::RobotFrame & frame,
                                     const Eigen::Vector3d normal,
                                     double lambda)
{
  return mc_rtc::make_void_ptr<mc_tvm::ImpulseFunctionPtr>(std::make_shared<mc_tvm::ImpulseFunction>(robot, frame, normal, lambda));
}

TVMImpulseConstraint::TVMImpulseConstraint(const mc_rbdyn::Robot & robot,
                                                 const mc_rbdyn::RobotFrame & frame,
                                                 const Eigen::Vector3d normal,
                                                 double lambda,
                                                 double delta_t,
                                                 double c_res,
                                                 double limit_multiplier)
: robot_(robot), frame_(frame), lambda_(lambda), delta_t_(delta_t), c_res_(c_res),
  limit_multiplier_(limit_multiplier), imp_constr_(initialize_imp_cstr(robot, frame, normal, lambda_)),
  upper_limit_(robot_.tvmRobot().limits().tu * limit_multiplier_ * (delta_t_/(c_res_+1.f))),
  lower_limit_(robot_.tvmRobot().limits().tl * limit_multiplier_ * (delta_t_/(c_res_+1.f)))
{
  mc_rtc::log::info("Size of upper_limit_ is {} x {}", upper_limit_.rows(), upper_limit_.cols());
  mc_rtc::log::info("Size of lower_limit_ is {} x {}", lower_limit_.rows(), lower_limit_.cols());
}

void TVMImpulseConstraint::addToSolver(mc_solver::TVMQPSolver & solver)
{
      auto & problem = tvm_solver(solver).problem();
      mc_rtc::log::info("Size of the lower limit is {}", lower_limit_.size());
      mc_tvm::ImpulseFunctionPtr imp_fn = *static_cast<mc_tvm::ImpulseFunctionPtr *>(imp_constr_.get());
      auto imp = problem.add( lower_limit_ <= imp_fn <= upper_limit_, tvm::task_dynamics::None(), {tvm::requirements::PriorityLevel(0)});
      constraints_.push_back(imp);
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

static mc_rtc::void_ptr initialize_tvm(const mc_rbdyn::Robot & robot, const mc_rbdyn::RobotFrame & frame, const Eigen::Vector3d normal, double lambda, double delta_t, double c_res, double limit_multiplier)
{
  return mc_rtc::make_void_ptr<TVMImpulseConstraint>(robot, frame, normal, lambda, delta_t, c_res, limit_multiplier);
}

static mc_rtc::void_ptr initialize(QPSolver::Backend backend,
                                   const mc_rbdyn::Robots & robots,
                                   unsigned int robotIndex,
                                   const mc_rbdyn::RobotFrame & frame,
                                   const Eigen::Vector3d normal,
                                   double lambda,
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
      return initialize_tvm(robots.robot(robotIndex), frame, normal, lambda, delta_t, c_res, limit_multiplier);
    default:
      mc_rtc::log::error_and_throw("[ImpulseConstraint] Not implemented for solver backend: {}", backend);
  }
}

ImpulseConstraint::ImpulseConstraint(const mc_rbdyn::Robots & robots, unsigned int robotIndex, const mc_rbdyn::RobotFrame & frame, const Eigen::Vector3d normal, double lambda, double delta_t, double c_res, double limit_multiplier, mc_rtc::Logger & logger)
: constraint_(initialize(backend_, robots, robotIndex, frame, normal, lambda, delta_t, c_res, limit_multiplier)), logger_(logger)
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
      logger_.removeLogEntry("ImpulseConstraintEvaluation");
      logger_.removeLogEntry("ImpulseConstraintLowerlimit");
      logger_.removeLogEntry("ImpulseConstraintUpperlimit");
      logger_.removeLogEntry("PredictedImpulsiveTorque");
      static_cast<TVMImpulseConstraint *>(constraint_.get())->removeFromSolver(tvm_solver(solver));
      break;
    default:
      break;
  }
}

void ImpulseConstraint::add_logs()
{
  logger_.addLogEntry("ImpulseConstraintEvaluation", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunction()->value();});

  logger_.addLogEntry("ImpulseConstraintLowerlimit", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->LowerLimit();});

  logger_.addLogEntry("ImpulseConstraintUpperlimit", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->UpperLimit();});

  logger_.addLogEntry("PredictedImpulsiveTorque", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->ImpulsiveTorqures();});

  logger_.addLogEntry("PredictedImpulsiveTorqueDerivative", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->ImpulsiveTorquresDerivative();});

  logger_.addLogEntry("ConstraintRightSideUpper", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->RightSideUpper();});

  logger_.addLogEntry("ConstraintRightSideLower", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->RightSideLower();});

  logger_.addLogEntry("Usable_Joint_Position", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunction()->JointPos();});

  logger_.addLogEntry("Usable_Joint_Velocity", this, [&, this]()
  {return static_cast<TVMImpulseConstraint *>(constraint_.get())->impFunction()->JointVel();});
}

Eigen::VectorXd & TVMImpulseConstraint::RightSideLower()
{
  Eigen::VectorXd lower_lim = robot_.tvmRobot().limits().tl*limit_multiplier_;
  Eigen::VectorXd impulsive_torques = ((c_res_+1.f)/delta_t_)*ImpulsiveTorqures();
  lower = lambda_ * (lower_lim - impulsive_torques);
  return lower; // lambda_*(robot_.tvmRobot().limits().tl-ImpulsiveTorqures());
}

Eigen::VectorXd & TVMImpulseConstraint::RightSideUpper()
{
  Eigen::VectorXd upper_lim = robot_.tvmRobot().limits().tu*limit_multiplier_;
  Eigen::VectorXd impulsive_torques = ((c_res_+1.f)/delta_t_)*ImpulsiveTorqures();
  upper = lambda_ * (upper_lim - impulsive_torques);
  return upper; // lambda_*(robot_.tvmRobot().limits().tl-ImpulsiveTorqures());
  Eigen::VectorXd acc_liit = robot_.tvmRobot().limits().tu;
}

} // namespace mc_solver

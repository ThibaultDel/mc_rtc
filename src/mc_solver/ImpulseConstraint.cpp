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
                                     double delta_t,
                                     double c_res,
                                     double limit_multiplier,
                                     int axis)
{
  return mc_rtc::make_void_ptr<mc_tvm::ImpulseFunctionPtr>(std::make_shared<mc_tvm::ImpulseFunction>(robot, frame, delta_t, c_res, limit_multiplier, axis));
}

TVMImpulseConstraint::TVMImpulseConstraint(const mc_rbdyn::Robot & robot,
                                                 const mc_rbdyn::RobotFrame & frame,
                                                 double delta_t,
                                                 double c_res,
                                                 double limit_multiplier,
                                                 int axis)
: robot_(robot), frame_(frame), delta_t_(delta_t), c_res_(c_res), limit_multiplier_(limit_multiplier), imp_constr_(initialize_imp_cstr(robot, frame, delta_t, c_res, limit_multiplier, axis)), upper_limit_(robot_.tvmRobot().limits().tu * limit_multiplier_ / (c_res_+1)), lower_limit_(robot_.tvmRobot().limits().tl * limit_multiplier_ / (c_res_+1))
{
}

void TVMImpulseConstraint::addToSolver(mc_solver::TVMQPSolver & solver)
{
      auto & problem = tvm_solver(solver).problem();
      mc_rtc::log::info("Size of the lower limit is {}", lower_limit_.size());
      mc_tvm::ImpulseFunctionPtr imp_fn = *static_cast<mc_tvm::ImpulseFunctionPtr *>(imp_constr_.get());
      const auto imp = problem.add( lower_limit_ <= imp_fn <= upper_limit_, tvm::task_dynamics::Constant(), {tvm::requirements::PriorityLevel(0)});
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

static mc_rtc::void_ptr initialize_tvm(const mc_rbdyn::Robot & robot, const mc_rbdyn::RobotFrame & frame, double delta_t, double c_res, double limit_multiplier, int axis)
{
  return mc_rtc::make_void_ptr<TVMImpulseConstraint>(robot, frame, delta_t, c_res, limit_multiplier, axis);
}

static mc_rtc::void_ptr initialize(QPSolver::Backend backend,
                                   const mc_rbdyn::Robots & robots,
                                   unsigned int robotIndex,
                                   const mc_rbdyn::RobotFrame & frame,
                                   double delta_t,
                                   double c_res,
                                   double limit_multiplier,
                                   int axis)
{
  switch(backend)
  {
    case QPSolver::Backend::Tasks:
      mc_rtc::log::error("No implementation for the ImpulseConstraint with the Tasks backend");
      assert(false);
    case QPSolver::Backend::TVM:
      return initialize_tvm(robots.robot(robotIndex), frame, delta_t, c_res, limit_multiplier, axis);
    default:
      mc_rtc::log::error_and_throw("[ImpulseConstraint] Not implemented for solver backend: {}", backend);
  }
}

ImpulseConstraint::ImpulseConstraint(const mc_rbdyn::Robots & robots, unsigned int robotIndex, const mc_rbdyn::RobotFrame & frame, double delta_t, double c_res, double limit_multiplier, int axis, mc_rtc::Logger & logger)
: constraint_(initialize(backend_, robots, robotIndex, frame, delta_t, c_res, limit_multiplier, axis)), logger_(logger)
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

}


} // namespace mc_solver

/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_solver/TVMQPSolver.h>
#include <mc_tvm/PeakTorqueFunction.h>
#include <tvm/ControlProblem.h>
#include <tvm/task_dynamics/None.h>
#include <vector>

namespace mc_solver
{

struct TVMPeakTorqueConstraint
{
  mc_tvm::PeakTorqueFunctionPtr fn_;
  std::vector<tvm::TaskWithRequirementsPtr> constraints_;

  TVMPeakTorqueConstraint(const mc_tvm::PeakTorqueFunctionPtr & fn)
    : fn_(fn)
  {
  }

  void addToSolver(mc_solver::TVMQPSolver & solver)
  {
    if (fn_)
    {
      constraints_.push_back(
          solver.problem().add(fn_ <= 0.0, tvm::task_dynamics::None(), {tvm::requirements::PriorityLevel(0)}));
    }
  }

  void removeFromSolver(mc_solver::TVMQPSolver & solver)
  {
    for (auto & c : constraints_)
    {
      solver.problem().remove(*c);
    }
    constraints_.clear();
  }
};

} // namespace mc_solver

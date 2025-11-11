/*
* Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_solver/TVMQPSolver.h>

#include <mc_tvm/Robot.h>

#include <tvm/ControlProblem.h>
#include <tvm/hint/internal/DiagonalCalculator.h>
#include <tvm/task_dynamics/Constant.h>
#include <tvm/task_dynamics/None.h>


namespace mc_solver
{

struct TVMImpulseConstraint
{
  const mc_rbdyn::Robot & robot_;
  mc_rbdyn::ConstRobotFramePtr frame_;
  const double delta_t_;
  const double c_res_;
  const double limit_multiplier_;
  std::vector<tvm::TaskWithRequirementsPtr> constraints_;
  std::vector<tvm::TaskWithRequirementsPtr> mimics_constraints_;

  TVMImpulseConstraint(const mc_rbdyn::Robot & robot, const mc_rbdyn::RobotFrame & frame, double delta_t, double c_res, double limit_multiplier, int axis);

  void addToSolver(mc_solver::TVMQPSolver & solver);

  void removeFromSolver(mc_solver::TVMQPSolver & solver);

  mc_tvm::ImpulseFunctionPtr impFunction() const{ return *static_cast<mc_tvm::ImpulseFunctionPtr *>(imp_constr_.get());}

  Eigen::VectorXd & LowerLimit() { return lower_limit_;}
  Eigen::VectorXd & UpperLimit() { return upper_limit_;}


protected:
  mc_rtc::void_ptr imp_constr_;

  Eigen::VectorXd upper_limit_;
  Eigen::VectorXd lower_limit_;

};

} // namespace mc_solver

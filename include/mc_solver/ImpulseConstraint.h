/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_rbdyn/Robots.h>
#include <mc_solver/ConstraintSet.h>
#include "mc_tvm/ImpulseFunction.h"

#include <mc_rtc/log/Logger.h>

#include <mc_rtc/void_ptr.h>

namespace mc_solver
{

/** \class ImpulseConstraint
 * Holds kinematic constraints (joint limits) for the robot
 */

struct MC_SOLVER_DLLAPI ImpulseConstraint : public ConstraintSet
{
public:
  /** Regular constraint constructor
   * Builds a regular joint limits constraint, prefer a damped joint limits
   * constraint in general
   * See tasks::qp::JointLimitsConstr for detail
   * \param robots The robots including the robot affected by this constraint
   * \param robotIndex The index of the robot affected by this constraint
   * \param timeStep Solver timestep
   */
  ImpulseConstraint(const mc_rbdyn::Robots & robots,
                        unsigned int robotIndex,
                        const mc_rbdyn::RobotFrame & frame,
                        double delta_t,
                        double c_res,
                        double limit_multiplier,
                        int axis,
                        mc_rtc::Logger & logger);

protected:
  /** Implementation of mc_solver::ConstraintSet::addToSolver */
  void addToSolverImpl(mc_solver::QPSolver & solver) override;
  /** Implementation of mc_solver::ConstraintSet::removeFromSolver */
  void removeFromSolverImpl(mc_solver::QPSolver & solver) override;
  /** Holds the constraint implementation
   *
   * In Tasks backend:
   * - tasks::qp::JointLimitsConstr for non-damped constraint
   * - tasks::qp::DampedJointLimitsConstr for damped-constraint
   *
   * In TVM backend:
   * - internal implementation
   *
   * The deleter carries the initial type of the constraint
   */
  mc_rtc::void_ptr constraint_;

  void add_logs();

  mc_rtc::Logger & logger_;

  // auto evaluation;
};

} // namespace mc_solver

/*
 * Copyright 2015-2026 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_tasks/WrenchTask.h>
#include <SpaceVecAlg/EigenTypedef.h>
#include <SpaceVecAlg/ForceVec.h>
#include <SpaceVecAlg/MotionVec.h>
#include <SpaceVecAlg/SpaceVecAlg>
#include <Eigen/src/Core/Matrix.h>

namespace mc_tasks
{
struct MC_TASKS_DLLAPI TorqueCartesianTask_test : public WrenchTask
{
public:
  /*! \brief Constructor
   *
   * \param solver QP solver instance
   *
   * \param frame Frame controlled by this task
   *
   * \param stiffness Task stiffness
   *
   * \param weight Task weight
   */
  TorqueCartesianTask_test(const mc_solver::QPSolver & solver,
                           const mc_rbdyn::RobotFrame & frame,
                           double stiffness = 100.0,
                           double weight = 500.0);

  /*! \brief Constructor
   *
   * \param solver QP solver instance
   *
   * \param bodyName Name of the body to control
   *
   * \param rIndex Index of the robot controlled by this task
   *
   * \param stiffness Task stiffness
   *
   * \param weight Task weight
   */
  TorqueCartesianTask_test(const mc_solver::QPSolver & solver,
                           const std::string & bodyName,
                           unsigned int rIndex,
                           double stiffness = 100.0,
                           double weight = 500.0);

  void reset() override;

  void setTorqueFeedforward(const Eigen::VectorXd & tau_ff);
  Eigen::VectorXd torqueGravity() const { return torqueGravityCompensation_; }
  void setCompensateGravity(bool compensate) { compensateGravity_ = compensate; }
  bool isCompensatingGravity() { return compensateGravity_; }

protected:
  void addToGUI(mc_rtc::gui::StateBuilder & gui) override;
  void addToLogger(mc_rtc::Logger & logger) override;
  void update(mc_solver::QPSolver & solver) override;

private:
  /** Robot handled by the task */
  const mc_rbdyn::Robots & robots_;
  unsigned int rIndex_;
  const int nbActuatedJoints; // Number of actuated joints (excluding floating base)
  Eigen::VectorXd torqueFeedforward_; // tau_ff

  /** True if the task is compensating gravity */
  bool compensateGravity_ = false;
  Eigen::VectorXd torqueGravityCompensation_; // tau_g
};

} // namespace mc_tasks

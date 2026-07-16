/*
 * Copyright 2015-2026 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_tasks/TorquePDCartesianTask.h>
#include <SpaceVecAlg/MotionVec.h>
#include <Eigen/src/Core/Matrix.h>

namespace mc_tasks
{
/*! \brief Cartesian-space PD task relative to a body, in torque control.
 *
 * The TorquePDRelativeCartesianTask computes desired joint torques using a proportional–
 * derivative (PD) control law in cartesian space. It relies on TorqueTask to find
 * the joint accelerations that minimize the error between the commanded joint
 * torques and the torques predicted by the robot dynamic model, while
 * respecting the constraints of the controller. As such, this task is intended
 * to be used in torque-control mode only.
 *
 * The commanded torque is computed as:
 *
 *   \f[
 *     \tau = J^T[K_p (x_d - x) + K_d (\dot{x}_d - \dot{x})]
 *   \f]
 *
 * where \f$x\f$ and \f$\dot{x}\f$ are the current cartesian position and velocity,
 * and \f$x_d\f$ and \f$\dot{x}_d\f$ are the desired cartesian position and velocity.
 * Cartesian position and velocity are expressed in the relative frame, and converted to the world frame to be handle by
 * the TorquePDCartesianTask.
 *
 * Additional torque components can be added to the command:
 *  - feedforward torques \f$\tau_{ff}\f$,
 *  - external torque compensation \f$\tau_{ext}\f$,
 *  - gravity compensation \f$\tau_{g}\f$.
 *
 * By default, all additional torque components are set to zero.
 *
 * The default targets are the current cartesian configuration for position
 * (\f$x_d = x\f$) and zero velocity for the desired velocity
 * (\f$\dot{x}_d = 0\f$).
 *
 */
struct MC_TASKS_DLLAPI TorquePDRelativeCartesianTask : public TorquePDCartesianTask
{
public:
  /*! \brief Constructor
   *
   * \param solver QP solver instance
   *
   * \param frame Frame controlled by this task
   *
   * \param relative Relative frame for the task.
   *
   * \param stiffness Task stiffness
   *
   * \param weight Task weight
   */
  TorquePDRelativeCartesianTask(const mc_solver::QPSolver & solver,
                                const mc_rbdyn::RobotFrame & frame,
                                const mc_rbdyn::Frame & relative,
                                double stiffness = 100.0,
                                double weight = 500.0);

  /*! \brief Constructor
   *
   * \param solver QP solver instance
   *
   * \param bodyName Name of the body to control
   *
   * \param relBodyName Name of the body relatively to which the end-effector
   * is controlled. If empty, defaults to the robot's base.
   *
   * \param robots Robots controlled by this task
   *
   * \param rIndex Index of the robot controlled by this task
   *
   * \param stiffness Task stiffness
   *
   * \param weight Task weight
   */
  TorquePDRelativeCartesianTask(const mc_solver::QPSolver & solver,
                                const std::string & bodyName,
                                const std::string & relBodyName,
                                const mc_rbdyn::Robots & robots,
                                unsigned int rIndex,
                                double stiffness = 100.0,
                                double weight = 500.0);

  /*! \brief Constructor that defaults to the robot's base as relative frame
   *
   * \param solver QP solver instance
   *
   * \param bodyName Name of the body to control
   *
   * \param robots Robots controlled by this task
   *
   * \param rIndex Index of the robot controlled by this task
   *
   * \param stiffness Task stiffness
   *
   * \param weight Task weight
   */
  TorquePDRelativeCartesianTask(const mc_solver::QPSolver & solver,
                                const std::string & bodyName,
                                const mc_rbdyn::Robots & robots,
                                unsigned int rIndex,
                                double stiffness = 100.0,
                                double weight = 500.0);

  void setPosTarget(const sva::PTransformd & xd) override; // xd is in the relative frame.
  void setVelTarget(const sva::MotionVecd & xd_dot) override; // xd_dot is in the relative frame.

  sva::PTransformd posTarget() override; // Return the position target in the relative frame.
  sva::MotionVecd velTarget() override; // Return the velocity target in the relative frame.

protected:
  void addToGUI(mc_rtc::gui::StateBuilder & gui) override;

private:
  void update(mc_solver::QPSolver & solver) override;
  void reset() override;
  mc_rbdyn::ConstFramePtr relative_;
  sva::PTransformd posTarget_rel_; // xd
  sva::MotionVecd velTarget_rel_; // xd_dot
};

} // namespace mc_tasks

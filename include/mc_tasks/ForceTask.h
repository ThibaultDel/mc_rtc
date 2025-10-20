/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_tasks/TrajectoryTaskGeneric.h>

#include <mc_rbdyn/RobotFrame.h>

namespace mc_tasks
{

/*! \brief Control a frame 6D force */
struct MC_TASKS_DLLAPI ForceTask : public TrajectoryTaskGeneric
{
public:
  /*! \brief Constructor
   *
   * \param frame Frame controlled by this task
   *
   * \param weight Task weight
   *
   * \param compensateExternalForces If true, the task will try to compensate the external forces acting on the robot
   *
   */
  ForceTask(const mc_rbdyn::RobotFrame & frame, double weight = 500.0, bool compensateExternalForces = false);

  /*! \brief Constructor
   *
   * Prefer the frame-based constructor
   *
   * \param frameName Name of the surface frame to control
   *
   * \param robots Robots controlled by this task
   *
   * \param robotIndex Index of the robot controlled by this task
   *
   * \param weight Task weight
   *
   * \param compensateExternalForces If true, the task will try to compensate for external forces acting on the robot
   *
   */
  ForceTask(const std::string & frameName,
            const mc_rbdyn::Robots & robots,
            unsigned int robotIndex,
            double weight = 500,
            bool compensateExternalForces = false);

  /*! \brief Reset the task
   *
   * Set the task target to the current frame position
   *
   * Reset its target velocity and acceleration to zero.
   *
   */
  void reset() override;

  /*! \brief Get the task's target */
  virtual sva::ForceVecd target() const;

  /*! \brief Set the task's target
   *
   * \param force Target in world frame
   *
   */
  virtual void target(const sva::ForceVecd & force);

  /*! \brief Retrieve the controlled frame name */
  inline const std::string & getFrameName() const noexcept { return frame_->name(); }

  /*! \brief Return the controlled frame (const) */
  const mc_rbdyn::RobotFrame & frame() const noexcept { return *frame_; }

  /** Returns the wrench of the frame in the inertial frame */
  inline sva::ForceVecd wrench() const noexcept { return frame_->wrench(); }

  /** Set if the task is compensating external forces */
  void compensateExternalForces(bool compensate);

  /** True if the task is compensating external forces */
  bool isCompensatingExternalForces() const;

  // /** Add support for the following criterias:
  //  *
  //  * - wrench: completed when the surface wrench reaches the given wrench, if
  //  *   some values are NaN, this direction is ignored. Only valid if the surface
  //  *   controlled by this task is attached to a force sensor, throws otherwise
  //  *
  //  * @throws If wrench is used but the surface is not attached to a force sensor
  //  */
  // std::function<bool(const mc_tasks::MetaTask & task, std::string &)> buildCompletionCriteria(
  //     double dt,
  //     const mc_rtc::Configuration & config) const override;

  void addToLogger(mc_rtc::Logger & logger) override;

  /*! \brief Load parameters from a Configuration object */
  void load(mc_solver::QPSolver & solver, const mc_rtc::Configuration & config) override;

protected:
  mc_rbdyn::ConstRobotFramePtr frame_;
  sva::ForceVecd curForce;
  bool compensateExternalForces_ = false;
  void addToGUI(mc_rtc::gui::StateBuilder & gui) override;
};

} // namespace mc_tasks

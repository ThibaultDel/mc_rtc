/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_tasks/TrajectoryTaskGeneric.h>

#include <mc_rbdyn/RobotFrame.h>
#include <SpaceVecAlg/SpaceVecAlg>

namespace mc_tasks
{

/*! \brief Control a frame 6D wrench */
struct MC_TASKS_DLLAPI WrenchTask : public TrajectoryTaskGeneric
{
public:
  /*! \brief Constructor
   *
   * \param frame Frame controlled by this task
   *
   * \param weight Task weight
   */
  WrenchTask(const mc_rbdyn::RobotFrame & frame, double weight = 500.0);

  /*! \brief Constructor
   *
   * Prefer the frame-based constructor
   *
   * \param surfaceName Name of the surface frame to control
   *
   * \param robots Robots controlled by this task
   *
   * \param robotIndex Index of the robot controlled by this task
   *
   * \param weight Task weight
   *
   */
  WrenchTask(const std::string & surfaceName,
             const mc_rbdyn::Robots & robots,
             unsigned int robotIndex,
             double weight = 500);

  /*! \brief Reset the task
   *
   * Set the task target to the current frame position
   *
   * Reset its target velocity and acceleration to zero.
   *
   */
  void reset() override;

  /*! \brief Get the task's target */
  sva::ForceVecd targetWrench() const;

  /*! \brief Set the task's target
   *
   * \param wrench Target in world frame
   *
   */
  void targetWrench(const sva::ForceVecd & worldWrench);

  /*! \brief Retrieve the controlled frame name */
  inline const std::string & frameName() const noexcept { return frame_->name(); }

  /*! \brief Return the controlled frame (const) */
  const mc_rbdyn::RobotFrame & frame() const noexcept { return *frame_; }

  /** Returns the wrench of the frame in the inertial frame */
  sva::ForceVecd currentWrench();

  ///** Returns the Transpose of the Jacobian mapping joint torques to the frame wrench */
  Eigen::MatrixXd dynamicJacobianTranspose() const noexcept;

  Eigen::Matrix6d cartesianInertia() const noexcept;

protected:
  mc_rbdyn::ConstRobotFramePtr frame_;
  void addToGUI(mc_rtc::gui::StateBuilder & gui) override;
  void addToLogger(mc_rtc::Logger & logger) override;
};

} // namespace mc_tasks

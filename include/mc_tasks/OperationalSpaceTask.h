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
/*! \brief Operational Space Controller (OSC) task.
 *
 * The OperationalSpaceTask computes desired wrench based on the Operational Space Formulation developed by Oussama
 * Khatib (1987) It delegates the resolution of joint accelerations to the WrenchTask, which computes accelerations that
 * best match the desired wrench while respecting controller constraints.
 *
 * The commanded wrench is computed as:
 *
 *   \f[
 *     F = Λ(q)\ddot{x} + μ(q, \dot{q}​) + p(q) - F_{ext}
 *   \f]
 *
 * This task provides interfaces to set the desired acceleration \f$\ddot{x}\f$ using a customizable PD control law, as
 * well as Integral control and feedforward terms:
 *
 *   \f[
 *     \ddot{x} = K_p (x_d - x) + K_d (\dot{x}_d - \dot{x})
 *          + K_i \int (x_d - x)\,dt
 *          + \ddot{x}_{ff}
 *   \f]
 *
 * where:
 *  - \f$K_i \int (x_d - x) dt\f$ is an optional integral term (with anti-windup),
 *  - \f$\ddot{x}_{ff}\f$ are feedforward accelerations,
 *  - \f$F_{ext}\f$ compensates external disturbances,
 *   - \f$μ(q, \dot{q}​) + p(q)\f$ is the Coriolis and gravity compensation in Cartesian space.
 *
 * By default, \f$K_i = 0\f$ and all additional wrench components are zero.
 *
 * By default, the desired position is the current cartesian configuration
 * (\f$x_d = x\f$) and the desired velocity to zero (\f$\dot{x}_d = 0\f$).
 *
 * Cartesian position, velocity and Jacobian are expressed in the world frame.
 */
struct MC_TASKS_DLLAPI OperationalSpaceTask : public WrenchTask
{
public:
  /*! \brief Constructor
   *
   * \param solver QP solver instance
   *
   * \param frame Frame controlled by this task
   *
   * \param weight Task weight
   */
  OperationalSpaceTask(const mc_solver::QPSolver & solver, const mc_rbdyn::RobotFrame & frame, double weight = 500.0);

  /*! \brief Constructor
   *
   * \param solver QP solver instance
   *
   * \param bodyName Name of the body to control
   *
   * \param rIndex Index of the robot controlled by this task
   *
   * \param weight Task weight
   */
  OperationalSpaceTask(const mc_solver::QPSolver & solver,
                       const std::string & bodyName,
                       unsigned int rIndex,
                       double weight = 500.0);

  void reset() override;

  // Setter
  void setStiffnessRotational(double kp)
  {
    setStiffness(sva::MotionVecd(Eigen::Vector3d::Constant(kp), stiffness_.linear()));
  }
  void setStiffnessTranslational(double kp)
  {
    setStiffness(sva::MotionVecd(stiffness_.angular(), Eigen::Vector3d::Constant(kp)));
  }
  void setDampingRotational(double kd)
  {
    setDamping(sva::MotionVecd(Eigen::Vector3d::Constant(kd), damping_.linear()));
  }
  void setDampingTranslational(double kd)
  {
    setDamping(sva::MotionVecd(damping_.angular(), Eigen::Vector3d::Constant(kd)));
  }
  void setStiffness(const sva::MotionVecd & stiffness) { stiffness_ = stiffness; }
  void setDamping(const sva::MotionVecd & damping) { damping_ = damping; }

  void enableIntegralTerm(bool enable);
  void setIntegralGainRotational(double ki)
  {
    setIntegralGain(sva::MotionVecd(Eigen::Vector3d::Constant(ki), integralGain_.angular()));
  }
  void setIntegralGainTranslational(double ki)
  {
    setIntegralGain(sva::MotionVecd(integralGain_.linear(), Eigen::Vector3d::Constant(ki)));
  }
  void setMaxIntegralWrench(double maxIntegralForce, double maxIntegralTorque)
  {
    setMaxIntegralWrench(
        sva::ForceVecd(Eigen::Vector3d::Constant(maxIntegralForce), Eigen::Vector3d::Constant(maxIntegralTorque)));
  }
  void setIntegralGain(const sva::MotionVecd & integralGain)
  {
    integralGain_ = integralGain;
    enableIntegralTerm(true);
  }
  void setMaxIntegralWrench(const sva::ForceVecd & maxIntegralWrench) { maxIntegralWrench_ = maxIntegralWrench; }
  void setMaxIntegralForce(double maxIntegralForce)
  {
    setMaxIntegralWrench(sva::ForceVecd(maxIntegralWrench_.couple(), Eigen::Vector3d::Constant(maxIntegralForce)));
  }
  void setMaxIntegralTorque(double maxIntegralTorque)
  {
    setMaxIntegralWrench(sva::ForceVecd(Eigen::Vector3d::Constant(maxIntegralTorque), maxIntegralWrench_.force()));
  }

  void setPosTarget(const sva::PTransformd & xd) { posTarget_ = xd; }
  void setVelTarget(const sva::MotionVecd & xd_dot) { velTarget_ = xd_dot; }
  void setAccelerationFeedforward(const sva::MotionVecd & xd_ddot) { accelerationFeedforward_ = xd_ddot; }
  void setDeriveVelocityTargetFromPosition(bool compute);
  void setCompensateExternalWrench(bool compensate) { compensateExternalWrench_ = compensate; }
  void setCompensateGravity(bool compensate) { compensateGravity_ = compensate; }

  // Getter
  const sva::MotionVecd & stiffness() const { return stiffness_; }
  const sva::MotionVecd & damping() const { return damping_; }
  const sva::PTransformd & posTarget() const { return posTarget_; }
  const sva::MotionVecd & velTarget() const { return velTarget_; }
  const sva::MotionVecd & accelerationFeedforward() const { return accelerationFeedforward_; }
  const sva::MotionVecd & integralGain() const { return integralGain_; }
  const sva::ForceVecd & maxIntegralWrench() const { return maxIntegralWrench_; }
  bool integralTermEnabled() const { return integralTermEnabled_; }
  bool deriveVelocityTargetFromPosition() const { return deriveVelocityTargetFromPosition_; }

  const sva::ForceVecd & externalWrench() const { return externalWrench_; }
  const sva::ForceVecd & gravityWrench() const { return gravityWrench_; }
  bool isCompensatingExternalWrench() const { return compensateExternalWrench_; }
  bool isCompensatingGravity() const { return compensateGravity_; }

protected:
  void addToGUI(mc_rtc::gui::StateBuilder & gui) override;
  void addToLogger(mc_rtc::Logger & logger) override;
  void update(mc_solver::QPSolver & solver) override;
  sva::PTransformd posTarget_; // xd
  sva::MotionVecd velTarget_; // xd_dot

private:
  /** Robot handled by the task */
  const mc_rbdyn::Robots & robots_;
  unsigned int rIndex_;

  const int nbActuatedJoints; // Number of actuated joints (excluding floating base)

  sva::MotionVecd stiffness_; // Kp - Unit: N/m (for position) or Nm/rad (for orientation).
  sva::MotionVecd damping_; // Kd - Unit: Ns/m (for position) or Nms/rad (for orientation).

  bool integralTermEnabled_ = false;
  sva::MotionVecd integralGain_; // Ki - Unit: N/(m*s) (for position) or Nm/(rad*s) (for orientation)
  sva::ForceVecd maxIntegralWrench_; // Anti-windup limits for the integral term
  sva::ForceVecd integralWrench_; // K_i * integralError_ - Unit: N (for force) or Nm (for torque)

  sva::MotionVecd accelerationFeedforward_; // \ddot{x}_ff
  /** True if the task is compensating external forces */
  bool compensateExternalWrench_ = false;
  sva::ForceVecd externalWrench_; // F_ext
  /** True if the task is compensating gravity */
  bool compensateGravity_ = true;
  sva::ForceVecd gravityWrench_; // μ + p (Cartesian Coriolis + gravity compensation)
  void computeForceGravityCompensation();

  sva::MotionVecd posError_;
  sva::MotionVecd velError_;
  sva::MotionVecd integralError_;

  sva::ForceVecd wrenchTarget_;

  bool deriveVelocityTargetFromPosition_ = false;
  sva::PTransformd prevPosTarget_;
};

} // namespace mc_tasks

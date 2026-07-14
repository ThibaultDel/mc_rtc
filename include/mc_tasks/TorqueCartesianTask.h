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
/*! \brief Cartesian-space PD task, in torque control.
 *
 * The TorqueCartesianTask computes desired wrench using a
 * proportional–derivative (PD) control law in cartesian space. It delegates the resolution of
 * joint accelerations to the WrenchTask, which computes accelerations that best
 * match the desired wrench while respecting controller constraints.
 *
 * This task is intended exclusively for torque-control mode.
 *
 * The commanded wrench is computed as:
 *
 *   \f[
 *     \tau = J^T[K_p (x_d - x) + K_d (\dot{x}_d - \dot{x})]
 *   \f]
 *
 *   \f[
 *     F = J#^T(\tau)
 *   \f]
 *
 * The full torque can include additional terms:
 *
 *   \f[
 *     \tau = K_p (x_d - x) + K_d (\dot{x}_d - \dot{x})
 *          + K_i \int (x_d - x)\,dt
 *          + \tau_{ff} + \tau_{ext} + \tau_{g}
 *   \f]
 *
 * where:
 *  - \f$K_i \int (x_d - x) dt\f$ is an optional integral term (with anti-windup),
 *  - \f$\tau_{ff}\f$ are feedforward torques,
 *  - \f$\tau_{ext}\f$ compensates external disturbances,
 *  - \f$\tau_{g}\f$ is gravity compensation.
 *
 * By default, \f$K_i = 0\f$ and all additional torque components are zero.
 *
 * By default, the desired position is the current cartesian configuration
 * (\f$x_d = x\f$) and the desired velocity to zero (\f$\dot{x}_d = 0\f$).
 *
 * Cartesian position, velocity and Jacobian are expressed in the world frame.
 */
struct MC_TASKS_DLLAPI TorqueCartesianTask : public WrenchTask
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
  TorqueCartesianTask(const mc_solver::QPSolver & solver,
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
  TorqueCartesianTask(const mc_solver::QPSolver & solver,
                      const std::string & bodyName,
                      unsigned int rIndex,
                      double stiffness = 100.0,
                      double weight = 500.0);

  void reset() override;

  void setStiffness(double stiffness) { setStiffness(Eigen::Vector6d::Constant(stiffness)); }
  void setDamping(double damping) { setDamping(Eigen::Vector6d::Constant(damping)); }
  void setStiffness(const Eigen::Vector6d & stiffness) { stiffness_ = stiffness; }
  void setDamping(const Eigen::Vector6d & damping) { damping_ = damping; }
  void enableIntegralTerm(bool enable);
  // Set the integral gain, activating the integral term if is not already enabled
  void setIntegralGain(double integralGain) { setIntegralGain(Eigen::Vector6d::Constant(integralGain)); }
  void setMaxIntegralWrench(double maxIntegralForce, double maxIntegralTorque)
  {
    setMaxIntegralWrench(
        sva::ForceVecd(Eigen::Vector3d::Constant(maxIntegralForce), Eigen::Vector3d::Constant(maxIntegralTorque)));
  }
  void setIntegralGain(const Eigen::Vector6d & integralGain)
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
  void setTorqueFeedforward(const Eigen::VectorXd & tau_ff);
  void setDeriveVelocityTargetFromPosition(bool compute);

  const Eigen::Vector6d & stiffness() const { return stiffness_; }
  const Eigen::Vector6d & damping() const { return damping_; }
  sva::PTransformd posTarget() { return posTarget_; }
  sva::MotionVecd velTarget() { return velTarget_; }
  const Eigen::VectorXd & torqueFeedforward() const { return torqueFeedforward_; }
  const Eigen::Vector6d & integralGain() const { return integralGain_; }
  const sva::ForceVecd & maxIntegralWrench() const { return maxIntegralWrench_; }
  bool integralTermEnabled() const { return integralTermEnabled_; }
  bool deriveVelocityTargetFromPosition() const { return deriveVelocityTargetFromPosition_; }

  Eigen::VectorXd torqueExternalForces() const { return torqueExtForcesCompensation_; }
  sva::ForceVecd gravityWrench() const { return wrenchGravityCompensation_; }
  void setCompensateExternalForces(bool compensate) { compensateExternalForces_ = compensate; }
  bool isCompensatingExternalForces() { return compensateExternalForces_; }
  void setCompensateGravity(bool compensate) { compensateGravity_ = compensate; }
  bool isCompensatingGravity() { return compensateGravity_; }

  bool useAutoDamping() const { return useAutoDamping_; }
  void setUseAutoDamping(bool use);
  Eigen::Vector6d autoDampingNaturalFrequency() const { return omega_n_; }
  void setAutoDampingNaturalFrequency(const Eigen::Vector6d & omega_n);
  Eigen::Vector6d autoDampingDampingRatio() const { return zeta_; }
  void setAutoDampingDampingRatio(const Eigen::Vector6d & zeta);

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

  Eigen::Vector6d stiffness_; // Kp
  Eigen::Vector6d damping_; // Kd
  bool integralTermEnabled_ = false;
  Eigen::Vector6d integralGain_; // Ki
  sva::ForceVecd maxIntegralWrench_; // Anti-windup limits for the integral term
  sva::ForceVecd integralWrench_; // K_i * integralError_

  Eigen::VectorXd torqueFeedforward_; // tau_ff
  /** True if the task is compensating external forces */
  bool compensateExternalForces_ = false;
  Eigen::VectorXd torqueExtForcesCompensation_; // tau_ext
  /** True if the task is compensating gravity */
  bool compensateGravity_ = false;
  sva::ForceVecd wrenchGravityCompensation_; // mu + p (Cartesian Coriolis + gravity compensation)
  void computeForceGravityCompensation();

  sva::MotionVecd posError_;
  sva::MotionVecd velError_;
  sva::MotionVecd integralError_;

  sva::ForceVecd wrenchTarget_;

  bool deriveVelocityTargetFromPosition_ = false;
  sva::PTransformd prevPosTarget_;

  bool useAutoDamping_ = false;
  Eigen::Vector6d omega_n_; // Natural frequency for critical damping (used if useAutoDamping_ is true)
  Eigen::Vector6d zeta_; // Damping ratio (used if useAutoDamping_ is true)
};

} // namespace mc_tasks

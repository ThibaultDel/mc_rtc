/*
 * Copyright 2015-2026 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_tasks/TorqueTask.h>
#include <Eigen/src/Core/Matrix.h>

namespace mc_tasks
{
/*! \brief Joint-space position task in torque control.
 *
 * The TorqueJointTask computes desired joint torques using a
 * proportional–derivative (PD) control law in joint space. It delegates the resolution of
 * joint accelerations to the TorqueTask, which computes accelerations that best
 * match the desired torques while respecting controller constraints.
 *
 * This task is intended exclusively for torque-control mode.
 *
 * The commanded torque is:
 *
 *   \f[
 *     \tau = K_p (q_d - q) + K_d (\dot{q}_d - \dot{q})
 *   \f]
 *
 * The full torque command can include additional terms:
 *
 *   \f[
 *     \tau = K_p (q_d - q) + K_d (\dot{q}_d - \dot{q})
 *          + K_i \int (q_d - q)\,dt
 *          + \tau_{ff} + \tau_{ext} + \tau_{g}
 *   \f]
 *
 * where:
 *  - \f$K_i \int (q_d - q) dt\f$ is an optional integral term (with anti-windup),
 *  - \f$\tau_{ff}\f$ are feedforward torques,
 *  - \f$\tau_{ext}\f$ compensates external disturbances,
 *  - \f$\tau_{g}\f$ is gravity compensation.
 *
 * By default, \f$K_i = 0\f$ and all additional torque components are zero.
 *
 * By default, the desired position is set to the current joint configuration
 * (\f$q_d = q\f$) and the desired velocity to zero (\f$\dot{q}_d = 0\f$).
 */
struct MC_TASKS_DLLAPI TorqueJointTask : public TorqueTask
{
public:
  TorqueJointTask(const mc_solver::QPSolver & solver, unsigned int rIndex, double stiffness, double weight);
  void reset() override;
  void setStiffness(double stiffness);
  void setDamping(double damping);
  void setStiffness(const Eigen::VectorXd & stiffness);
  void setDamping(const Eigen::VectorXd & damping);
  void enableIntegralTerm(bool enable);
  // Set the integral gain, activating the integral term if is not already enabled
  void setIntegralGain(double integralGain);
  void setMaxIntegralTorque(double maxIntegralTorque);
  void setIntegralGain(const Eigen::VectorXd & integralGain);
  void setMaxIntegralTorque(const Eigen::VectorXd & maxIntegralTorque);

  void setPosTarget(const Eigen::VectorXd & qd);
  void setVelTarget(const Eigen::VectorXd & qd_dot);
  void setTorqueFeedforward(const Eigen::VectorXd & tau_ff);
  void setDeriveVelocityTargetFromPosition(bool compute);

  const Eigen::VectorXd & stiffness() const { return stiffness_; }
  const Eigen::VectorXd & damping() const { return damping_; }
  const Eigen::VectorXd & posTarget() const { return posTarget_; }
  const Eigen::VectorXd & velTarget() const { return velTarget_; }
  const Eigen::VectorXd & torqueFeedforward() const { return torqueFeedforward_; }
  const Eigen::VectorXd & integralGain() const { return integralGain_; }
  const Eigen::VectorXd & maxIntegralTorque() const { return maxIntegralTorque_; }
  bool integralTermEnabled() const { return integralTermEnabled_; }
  bool deriveVelocityTargetFromPosition() const { return deriveVelocityTargetFromPosition_; }

protected:
  void addToGUI(mc_rtc::gui::StateBuilder & gui) override;
  void addToLogger(mc_rtc::Logger & logger) override;

private:
  void update(mc_solver::QPSolver & solver) override;

  /** Robot handled by the task */
  const mc_rbdyn::Robots & robots_;
  unsigned int rIndex_;

  const int nbActuatedJoints; // Number of actuated joints (excluding floating base)

  Eigen::VectorXd stiffness_; // Kp
  Eigen::VectorXd damping_; // Kd
  bool integralTermEnabled_;
  Eigen::VectorXd integralGain_; // Ki
  Eigen::VectorXd maxIntegralTorque_; // Anti-windup limits for the integral term in Nm
  Eigen::VectorXd integralTorque_; // K_i * integralError_

  Eigen::VectorXd posTarget_; // qd
  Eigen::VectorXd velTarget_; // qd_dot
  Eigen::VectorXd torqueFeedforward_; // tau_ff

  Eigen::VectorXd posError_;
  Eigen::VectorXd velError_;
  Eigen::VectorXd integralError_;

  bool deriveVelocityTargetFromPosition_;
  Eigen::VectorXd prevPosTarget_;
};

} // namespace mc_tasks

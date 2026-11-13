/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_solver/api.h>
#include <mc_solver/ConstraintSet.h>
#include <mc_tvm/PeakTorqueFunction.h>
#include <mc_rbdyn/Robot.h>
#include <mc_rbdyn/RobotFrame.h>

#include <memory>
#include <string>
#include <vector>

namespace mc_solver
{

struct MC_SOLVER_DLLAPI JointFrictionConfig
{
  double K = 5000.0;                   ///< Presliding joint elasticity [Nm/rad]
  double B = 50.0;                     ///< Presliding joint damping [Nms/rad]
  double F_S = 2.5;                    ///< Breakaway static friction torque [Nm]
  double F_C = 0.4;                    ///< Coulomb friction torque [Nm]
  double v_S = 0.06;                   ///< Stribeck velocity [rad/s]
  double D = 4.5;                      ///< Viscous friction coefficient [Nms/rad]
};

struct MC_SOLVER_DLLAPI PeakTorqueConfig
{
  double gamma = 25.0;                 ///< CBF decay gain gamma
  double limit_multiplier = 0.85;      ///< Safety limit scaling (e.g. 0.85 = 85% of rated limit)
  double worker_rate_hz = 50.0;        ///< Background worker thread rate in Hz
  double c_res = 1.0;                  ///< Coefficient of restitution
  double delta_t = 0.005;              ///< Collision shock duration [s]
  double horizon = 0.050;              ///< Forward sensitivity simulation horizon [s]

  // Default elastoplastic friction model parameters (fallback if joint not in joint_friction)
  double K = 5000.0;                   ///< Presliding joint elasticity [Nm/rad]
  double B = 50.0;                     ///< Presliding joint damping [Nms/rad]
  double F_S = 2.5;                    ///< Breakaway static friction torque [Nm]
  double F_C = 0.4;                    ///< Coulomb friction torque [Nm]
  double v_S = 0.06;                   ///< Stribeck velocity [rad/s]
  double D = 4.5;                      ///< Viscous friction coefficient [Nms/rad]

  // Optional per-joint friction parameters map (jointName -> JointFrictionConfig)
  std::map<std::string, JointFrictionConfig> joint_friction;

  // Optional manual contact surfaces (if empty, harvested dynamically from solver contacts)
  std::vector<std::string> contact_surfaces;
  // Optional custom rotor inertias J_m (if empty, default 0.05 kg*m^2 is used)
  Eigen::VectorXd rotor_inertias;
  // Optional custom torque limits (if empty, extracted from robot limits)
  Eigen::VectorXd custom_torque_limits;
};

struct TVMPeakTorqueConstraint;

/**
 * @brief Control Barrier Function (CBF) constraint bounding predicted peak drive-side reaction torques.
 *
 * Implements a dual-rate architecture:
 * - A background worker computes the contact nullspace projection and multi-rate elastoplastic
 *   shock sensitivity Jacobian at 50-100 Hz.
 * - The 1 kHz TVM QP evaluates the acceleration barrier constraint (J_v * alphaD <= gamma * slack)
 *   lock-free using double-buffering.
 */
class MC_SOLVER_DLLAPI PeakTorqueConstraint : public ConstraintSet
{
public:
  PeakTorqueConstraint(
      const mc_rbdyn::Robot & robot,
      const mc_rbdyn::RobotFrame & impact_frame,
      const Eigen::Vector3d & impact_normal,
      const PeakTorqueConfig & config = PeakTorqueConfig());

  PeakTorqueConstraint(
      const mc_rbdyn::Robot & robot,
      const std::string & impact_frame_name,
      const Eigen::Vector3d & impact_normal,
      const PeakTorqueConfig & config = PeakTorqueConfig());

  ~PeakTorqueConstraint() override;

  void addToSolverImpl(mc_solver::QPSolver & solver) override;
  void removeFromSolverImpl(mc_solver::QPSolver & solver) override;
  void update(mc_solver::QPSolver & solver) override;

  // Background worker lifecycle
  void startWorker();
  void stopWorker();

  // Telemetry accessors
  const mc_tvm::PeakTorqueFunctionPtr & function() const;
  const Eigen::VectorXd & predictedPeakTorques() const;
  const Eigen::VectorXd & torqueLimits() const;
  const Eigen::VectorXd & cbfMargins() const;
  const Eigen::MatrixXd & jacobianV() const;
  double workerRateHz() const;
  double lastSolveTimeMs() const;
  uint64_t updateSeq() const;

  // Configuration setters
  void gamma(double g);
  double gamma() const;
  void limitMultiplier(double m);
  double limitMultiplier() const;
  void impactNormal(const Eigen::Vector3d & n);
  const Eigen::Vector3d & impactNormal() const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace mc_solver

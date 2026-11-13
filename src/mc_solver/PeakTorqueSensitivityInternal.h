/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/SVD>
#include <vector>
#include <cmath>

namespace mc_solver::internal
{

struct JointFrictionParams
{
  double K = 1000.0;    // Presliding stiffness [Nm/rad]
  double B = 50.0;      // Presliding damping [Nms/rad]
  double F_S = 15.0;    // Static friction / breakaway torque [Nm]
  double F_C = 10.0;    // Coulomb friction [Nm]
  double v_S = 0.05;    // Stribeck velocity [rad/s]
  double D = 5.0;       // Viscous damping [Nms/rad]

  double h = 0.001;     // Timestep [s]
  double Z = 0.0;       // Z = 1 / (K * h + B)

  void initialize(double timestep)
  {
    h = timestep;
    Z = 1.0 / (K * h + B);
  }
};

class CisnerosFrictionElement
{
public:
  CisnerosFrictionElement() = default;
  explicit CisnerosFrictionElement(const JointFrictionParams & params) : params_(params) {}

  double evaluate_force(double v_star, int & regime) const
  {
    double f_trial = params_.K * e_ + params_.B * v_star;
    if (std::abs(f_trial) <= params_.F_S)
    {
      regime = 0; // STICK
      return f_trial;
    }
    if (f_trial > params_.F_S)
    {
      regime = 1; // SLIDE+
      return params_.F_S;
    }
    regime = -1; // SLIDE-
    return -params_.F_S;
  }

  double evaluate_derivative(double v_star) const
  {
    double f_trial = params_.K * e_ + params_.B * v_star;
    if (std::abs(f_trial) <= params_.F_S)
    {
      return params_.B;
    }
    return 0.0;
  }

  void update_deformation(double e_prev, double f_current)
  {
    e_ = params_.Z * (params_.B * e_prev + params_.h * f_current);
  }

  void reset() { e_ = 0.0; }
  double deformation() const { return e_; }
  void set_timestep(double ts) { params_.initialize(ts); }
  const JointFrictionParams & params() const { return params_; }

private:
  JointFrictionParams params_;
  double e_ = 0.0;
};

struct ReductionResults
{
  int N = 0;              // Total DOFs
  int n = 0;              // Actuated DOFs
  int m = 0;              // Contact constraints
  int rank = 0;           // Numerical rank of J_c
  int nr = 0;             // Reduced dimensions N - rank

  Eigen::MatrixXd S;      // Joint selection matrix (n x N)
  Eigen::MatrixXd M0;     // Mass matrix with reflected motor inertias (N x N)
  Eigen::MatrixXd Z;      // Contact null-space orthonormal basis (N x nr)
  Eigen::MatrixXd Mc;     // Reduced mass matrix (nr x nr)
  Eigen::MatrixXd A;      // Joint mapping matrix S * Z (n x nr)
  Eigen::MatrixXd Jm_A_Mc_inv; // Motor shock projector: diag(J_m) * A * Mc^-1 (n x nr)
  Eigen::MatrixXd P_tau;  // Presliding reaction projector: I - J_m * A * Mc^-1 * A^T (n x n)
};

class ModelReduction
{
public:
  explicit ModelReduction(double svd_threshold = 1e-10) : svd_threshold_(svd_threshold) {}

  bool reduce(
      const Eigen::MatrixXd & M_robot,
      const Eigen::VectorXd & J_m,
      const Eigen::MatrixXd & J_c,
      int total_dofs,
      int actuated_dofs,
      ReductionResults & out) const;

private:
  double svd_threshold_ = 1e-10;
};

struct SensitivityOutput
{
  Eigen::VectorXd tau_r_peak;    ///< (n x 1) Peak absolute reaction torques [Nm]
  Eigen::MatrixXd dTauMax_dvpre; ///< (n x N) Sensitivity Jacobian w.r.t pre-impact generalized velocity
  bool success = false;
};

class PeakTorqueSensitivitySolver
{
public:
  PeakTorqueSensitivitySolver(
      const ReductionResults & red,
      const Eigen::VectorXd & J_m,
      const std::vector<JointFrictionParams> & joint_params,
      double delta_t = 0.005,
      double c_res = 1.0);

  SensitivityOutput solve(
      const Eigen::VectorXd & v_pre,
      const Eigen::VectorXd & J_imp,
      double horizon = 0.050);

private:
  ReductionResults red_;
  Eigen::VectorXd J_m_;
  std::vector<JointFrictionParams> joint_params_;
  double delta_t_ = 0.005;
  double c_res_ = 1.0;
};

} // namespace mc_solver::internal

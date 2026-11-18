/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include "PeakTorqueSensitivityInternal.h"
#include <mc_rtc/logging.h>
#include <cmath>

namespace mc_solver::internal
{

bool ModelReduction::reduce(
    const Eigen::MatrixXd & M_robot,
    const Eigen::VectorXd & J_m,
    const Eigen::MatrixXd & J_c,
    int total_dofs,
    int actuated_dofs,
    ReductionResults & out) const
{
  out.N = total_dofs;
  out.n = actuated_dofs;
  out.m = static_cast<int>(J_c.rows());

  const int startParam = out.N - out.n;
  out.S = Eigen::MatrixXd::Zero(out.n, out.N);
  out.S.block(0, startParam, out.n, out.n).setIdentity();

  // 1. Reflected motor inertia mapping: M0 = M_robot + S^T * diag(J_m) * S
  out.M0 = M_robot + out.S.transpose() * J_m.asDiagonal() * out.S;

  if (out.m == 0)
  {
    // Fixed-base robot or no external contact constraints:
    // Contact null-space is the entire unconstrained joint space
    out.rank = 0;
    out.nr = out.N;
    out.Z = Eigen::MatrixXd::Identity(out.N, out.nr);
    out.Mc = out.M0;
    out.A = out.S;

    Eigen::LLT<Eigen::MatrixXd> llt(out.Mc);
    if (llt.info() != Eigen::Success)
    {
      mc_rtc::log::warning("[PeakTorqueSensitivityInternal] Unconstrained mass matrix Mc Cholesky LLT decomposition failed!");
      return false;
    }

    Eigen::MatrixXd Mc_inv = llt.solve(Eigen::MatrixXd::Identity(out.nr, out.nr));
    out.Jm_A_Mc_inv = J_m.asDiagonal() * out.A * Mc_inv;
    out.P_tau = Eigen::MatrixXd::Identity(out.n, out.n) - out.Jm_A_Mc_inv * out.A.transpose();
    return true;
  }

  // 2. SVD of contact Jacobian J_c: J_c = U * Sigma * V^T
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(J_c, Eigen::ComputeFullV);
  svd.setThreshold(svd_threshold_);

  out.rank = static_cast<int>(svd.rank());
  out.nr = static_cast<int>(out.N - out.rank);
  if (out.nr <= 0)
  {
    mc_rtc::log::warning("[PeakTorqueSensitivityInternal] Reduced DOF dimension nr <= 0 (rank={}, N={})", out.rank, out.N);
    return false;
  }

  // Right null-space basis: last nr columns of V
  out.Z = svd.matrixV().rightCols(out.nr);

  // 3. Reduced mass matrix Mc and joint mapping matrix A
  out.Mc = out.Z.transpose() * out.M0 * out.Z;
  out.A = out.S * out.Z;

  // 4. Cholesky factorization LLT check of Mc
  Eigen::LLT<Eigen::MatrixXd> llt(out.Mc);
  if (llt.info() != Eigen::Success)
  {
    mc_rtc::log::warning("[PeakTorqueSensitivityInternal] Reduced mass matrix Mc Cholesky LLT decomposition failed!");
    return false;
  }

  // 5. Precompute motor shock projector: Jm_A_Mc_inv = diag(J_m) * A * Mc^-1
  Eigen::MatrixXd Mc_inv = llt.solve(Eigen::MatrixXd::Identity(out.nr, out.nr));
  out.Jm_A_Mc_inv = J_m.asDiagonal() * out.A * Mc_inv;

  // 6. Reaction torque projector: P_tau = I - J_m * A * Mc^-1 * A^T
  out.P_tau = Eigen::MatrixXd::Identity(out.n, out.n) - out.Jm_A_Mc_inv * out.A.transpose();

  return true;
}

PeakTorqueSensitivitySolver::PeakTorqueSensitivitySolver(
    const ReductionResults & red,
    const Eigen::VectorXd & J_m,
    const std::vector<JointFrictionParams> & joint_params,
    double delta_t,
    double c_res)
  : red_(red),
    J_m_(J_m),
    joint_params_(joint_params),
    delta_t_(delta_t),
    c_res_(c_res)
{
}

SensitivityOutput PeakTorqueSensitivitySolver::solve(
    const Eigen::VectorXd & v_pre,
    const Eigen::VectorXd & J_imp,
    double horizon)
{
  SensitivityOutput out;
  const int n = red_.n;
  const int N = red_.N;
  const int nr = red_.nr;

  out.tau_r_peak = Eigen::VectorXd::Zero(n);
  out.dTauMax_dvpre = Eigen::MatrixXd::Zero(n, N);

  // 1. Constrained Operational Space Effective Mass
  Eigen::VectorXd w = red_.Z.transpose() * J_imp; // (nr x 1)
  Eigen::LLT<Eigen::MatrixXd> llt_Mc(red_.Mc);
  if (llt_Mc.info() != Eigen::Success)
  {
    mc_rtc::log::warning("[PeakTorqueSensitivityInternal] Cholesky LLT of Mc failed in solve!");
    return out;
  }

  Eigen::VectorXd Mc_inv_w = llt_Mc.solve(w);
  double denom = w.dot(Mc_inv_w);
  if (denom <= 1e-9)
  {
    mc_rtc::log::warning("[PeakTorqueSensitivityInternal] Non-positive operational space denominator in impact law ({:.4e})", denom);
    return out;
  }
  double M_eq = 1.0 / denom;

  // 2. Linear reduced impulse sensitivity: d(p_c) / dv_pre (nr x N)
  // p_c = -(1 + c_res) * M_eq * w * (J_imp^T * v_pre)
  Eigen::MatrixXd dpc_dvpre = -(1.0 + c_res_) * M_eq * (w * J_imp.transpose());

  // Momentum impulse magnitude
  double v_n = J_imp.dot(v_pre);
  double p_ext = -(1.0 + c_res_) * M_eq * v_n;
  Eigen::VectorXd p_c = w * p_ext;

  // 3. Multi-rate timestepping setup
  const double h_shock = 0.0001; // 0.1 ms fine resolution for dynamic shock pulse
  const int N_shock = static_cast<int>(std::round(delta_t_ / h_shock));

  const double h_ringdown = 0.001; // 1.0 ms nominal step for ringdown
  const double t_ringdown = std::max(0.0, horizon - delta_t_);
  const int N_ringdown = static_cast<int>(std::round(t_ringdown / h_ringdown));

  std::vector<CisnerosFrictionElement> friction_elements(n);
  for (int i = 0; i < n; ++i)
  {
    friction_elements[i] = CisnerosFrictionElement(joint_params_[i]);
  }

  // Initial states at t = 0 (pre-impact consistent state)
  Eigen::VectorXd z = red_.Z.transpose() * v_pre;
  Eigen::MatrixXd dz_dv = red_.Z.transpose(); // dz(0) / dv_pre = Z^T (nr x N)
  Eigen::MatrixXd de_dv = Eigen::MatrixXd::Zero(n, N);

  // Helper lambda for shock pulse shaping
  auto eval_pulse = [this](double t) -> double {
    if (t < 0.0 || t > delta_t_ || delta_t_ <= 1e-9) return 0.0;
    return (M_PI / (2.0 * delta_t_)) * std::sin(M_PI * t / delta_t_);
  };

  // Integration loop
  int total_steps = N_shock + N_ringdown;
  double t_current = 0.0;

  for (int step = 1; step <= total_steps; ++step)
  {
    bool in_shock = (step <= N_shock);
    double h_step = in_shock ? h_shock : h_ringdown;
    t_current += h_step;

    for (int i = 0; i < n; ++i)
    {
      friction_elements[i].set_timestep(h_step);
    }

    // External contact torque pulse
    Eigen::VectorXd tau_c_ext = Eigen::VectorXd::Zero(nr);
    Eigen::MatrixXd dtau_c_ext_dv = Eigen::MatrixXd::Zero(nr, N);

    if (in_shock)
    {
      double s_val = eval_pulse(t_current);
      tau_c_ext = s_val * p_c;
      dtau_c_ext_dv = s_val * dpc_dvpre;
    }

    // Newton solver for discrete Euler step:
    // Mc * (z_next - z) = h * (tau_c_ext - A^T * f_next)
    Eigen::VectorXd z_next = z;
    Eigen::VectorXd f_next = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd df_du = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd df_de = Eigen::VectorXd::Zero(n);

    // Save e_k before update
    Eigen::VectorXd e_k(n);
    for (int i = 0; i < n; ++i) e_k(i) = friction_elements[i].deformation();

    // 1-step / Newton solve
    for (int iter = 0; iter < 10; ++iter)
    {
      Eigen::VectorXd u_next = red_.A * z_next;
      for (int i = 0; i < n; ++i)
      {
        const auto & p = friction_elements[i].params();
        double vs = u_next(i) + p.Z * p.K * e_k(i);
        int regime = 0;
        f_next(i) = friction_elements[i].evaluate_force(vs, regime);
        double df_dvs = friction_elements[i].evaluate_derivative(vs);
        df_du(i) = df_dvs;
        df_de(i) = df_dvs * p.Z * p.K;
      }

      Eigen::VectorXd res_step = red_.Mc * (z_next - z) - h_step * (tau_c_ext - red_.A.transpose() * f_next);
      if (res_step.norm() < 1e-7) break;

      Eigen::MatrixXd J_newton = red_.Mc + h_step * red_.A.transpose() * df_du.asDiagonal() * red_.A;
      z_next -= J_newton.ldlt().solve(res_step);
    }

    // Reaction torque at timestep:
    // tau_r = Jm_A_Mc_inv * tau_c_ext + P_tau * f_next
    Eigen::VectorXd tau_shock_rotor = red_.Jm_A_Mc_inv * tau_c_ext;
    Eigen::VectorXd tau_r = tau_shock_rotor + red_.P_tau * f_next;

    // Sensitivity propagation via Implicit Function Theorem:
    // J_newton * (dz_next / dv) = Mc * (dz / dv) - h * A^T * diag(df/de) * (de / dv) + h * dtau_c_ext_dv
    Eigen::MatrixXd J_sens = red_.Mc + h_step * red_.A.transpose() * df_du.asDiagonal() * red_.A;
    Eigen::LLT<Eigen::MatrixXd> llt_J(J_sens);
    if (llt_J.info() != Eigen::Success)
    {
      mc_rtc::log::warning("[PeakTorqueSensitivityInternal] Cholesky LLT of J_sens failed at step {}", step);
      return out;
    }

    Eigen::MatrixXd RHS = red_.Mc * dz_dv - h_step * red_.A.transpose() * df_de.asDiagonal() * de_dv + h_step * dtau_c_ext_dv;
    Eigen::MatrixXd dz_next_dv = llt_J.solve(RHS);

    Eigen::MatrixXd df_dv = df_du.asDiagonal() * red_.A * dz_next_dv + df_de.asDiagonal() * de_dv;
    Eigen::MatrixXd dtau_r_dv = red_.Jm_A_Mc_inv * dtau_c_ext_dv + red_.P_tau * df_dv;

    // Peak tracking
    for (int i = 0; i < n; ++i)
    {
      double abs_r = std::abs(tau_r(i));
      if (abs_r >= out.tau_r_peak(i))
      {
        out.tau_r_peak(i) = abs_r;
        double sgn = (tau_r(i) >= 0.0) ? 1.0 : -1.0;
        out.dTauMax_dvpre.row(i) = sgn * dtau_r_dv.row(i);
      }
    }

    // Update deformation state and deformation sensitivity
    Eigen::MatrixXd de_next_dv(n, N);
    for (int i = 0; i < n; ++i)
    {
      const auto & p = friction_elements[i].params();
      friction_elements[i].update_deformation(e_k(i), f_next(i));
      de_next_dv.row(i) = p.Z * (p.B * de_dv.row(i) + h_step * df_dv.row(i));
    }

    z = z_next;
    dz_dv = dz_next_dv;
    de_dv = de_next_dv;
  }

  out.success = true;
  return out;
}

} // namespace mc_solver::internal

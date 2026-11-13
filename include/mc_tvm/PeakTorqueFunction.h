/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_tvm/api.h>
#include <mc_tvm/Robot.h>
#include <mc_rbdyn/Robot.h>

#include <tvm/function/abstract/LinearFunction.h>

#include <Eigen/Core>
#include <vector>
#include <memory>
#include <atomic>
#include <array>

namespace mc_tvm
{

/**
 * @brief Thread-safe double-buffer struct storing latest sensitivity calculation results.
 */
struct PeakTorqueBuffer
{
  Eigen::MatrixXd J_v;        ///< (n x N) full pre-impact velocity sensitivity (floating-base + actuated)
  Eigen::VectorXd tau_r_peak; ///< (n x 1) predicted peak absolute reaction torques [Nm]
  uint64_t seq = 0;
};

/**
 * @brief TVM LinearFunction enforcing Control Barrier Function (CBF) on predicted peak impact torques:
 *        J_v * alphaD <= gamma * (tau_limit - tau_r_peak(v))
 *
 * Implemented as a LinearFunction f(alphaD) = J_full * alphaD + b <= 0,
 * where J_full = J_v in R^(n x N) directly mapping to the full robot acceleration alphaD in R^N,
 * and b_i = -gamma * (tau_limit_i - tau_r_peak_i).
 */
class MC_TVM_DLLAPI PeakTorqueFunction : public tvm::function::abstract::LinearFunction
{
public:
  using Output = tvm::function::abstract::LinearFunction::Output;
  DISABLE_OUTPUTS(Output::JDot)
  SET_UPDATES(PeakTorqueFunction, Jacobian, B)

  PeakTorqueFunction(const mc_rbdyn::Robot & robot,
                     double gamma,
                     double limit_multiplier,
                     const std::shared_ptr<const std::array<PeakTorqueBuffer, 2>> & buffers,
                     const std::shared_ptr<const std::atomic<int>> & read_idx);

  void updateb();
  void updateJacobian();

  // Parameter setters
  void gamma(double g) { gamma_ = g; }
  double gamma() const { return gamma_; }

  void limitMultiplier(double m);
  double limitMultiplier() const { return limit_multiplier_; }

  void customTorqueLimits(const Eigen::VectorXd & limits);

  // Telemetry getters for logging and GUI
  const Eigen::VectorXd & predictedPeakTorques() const { return tau_r_peak_; }
  const Eigen::VectorXd & torqueLimits() const { return tau_limit_; }
  const Eigen::VectorXd & cbfMargins() const { return cbf_margin_; }
  const Eigen::VectorXd & rightHandSide() const { return rhs_; }
  const Eigen::MatrixXd & jacobianV() const { return J_v_; }

private:
  const mc_rbdyn::Robot & robot_;
  double gamma_ = 25.0;
  double limit_multiplier_ = 0.85;

  std::shared_ptr<const std::array<PeakTorqueBuffer, 2>> buffers_;
  std::shared_ptr<const std::atomic<int>> read_idx_;

  Eigen::MatrixXd J_full_;          // (n x N)
  Eigen::MatrixXd J_v_;             // (n x N) full pre-impact sensitivity
  Eigen::VectorXd tau_limit_;       // (n x 1)
  Eigen::VectorXd tau_r_peak_;      // (n x 1)
  Eigen::VectorXd cbf_margin_;      // (n x 1)
  Eigen::VectorXd rhs_;             // (n x 1)
  int startParam_ = 6;
};

using PeakTorqueFunctionPtr = std::shared_ptr<PeakTorqueFunction>;

} // namespace mc_tvm

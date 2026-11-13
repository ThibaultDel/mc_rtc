/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tvm/PeakTorqueFunction.h>
#include <mc_rtc/logging.h>
#include <cassert>

namespace mc_tvm
{

PeakTorqueFunction::PeakTorqueFunction(
    const mc_rbdyn::Robot & robot,
    double gamma,
    double limit_multiplier,
    const std::shared_ptr<const std::array<PeakTorqueBuffer, 2>> & buffers,
    const std::shared_ptr<const std::atomic<int>> & read_idx)
  : tvm::function::abstract::LinearFunction(robot.mb().nrDof() - (robot.tvmRobot().qFloatingBase() ? robot.tvmRobot().qFloatingBase()->space().tSize() : 0)),
    robot_(robot),
    gamma_(gamma),
    limit_multiplier_(limit_multiplier),
    buffers_(buffers),
    read_idx_(read_idx)
{
  auto & tvm_robot = robot_.tvmRobot();
  registerUpdates(Update::B, &PeakTorqueFunction::updateb);
  registerUpdates(Update::Jacobian, &PeakTorqueFunction::updateJacobian);
  addOutputDependency<PeakTorqueFunction>(Output::B, Update::B);
  addOutputDependency<PeakTorqueFunction>(Output::Jacobian, Update::Jacobian);
  addVariable(tvm_robot.alphaD(), true);

  startParam_ = tvm_robot.qFloatingBase() ? tvm_robot.qFloatingBase()->space().tSize() : 0;
  const int n = robot_.mb().nrDof() - startParam_;
  const int N = robot_.mb().nrDof();

  J_full_ = Eigen::MatrixXd::Zero(n, N);
  J_v_ = Eigen::MatrixXd::Zero(n, N);
  b_ = Eigen::VectorXd::Zero(n);
  tau_r_peak_ = Eigen::VectorXd::Zero(n);
  cbf_margin_ = Eigen::VectorXd::Zero(n);
  rhs_ = Eigen::VectorXd::Zero(n);

  limitMultiplier(limit_multiplier_);
  mc_rtc::log::info("[PeakTorqueFunction] Created TVM CBF constraint for {} actuated joints (floating base = {} DOFs)",
                    n, startParam_);
}

void PeakTorqueFunction::limitMultiplier(double m)
{
  limit_multiplier_ = m;
  const auto & tu = robot_.tvmRobot().limits().tu;
  const int n = robot_.mb().nrDof() - startParam_;
  if (tu.size() >= startParam_ + n)
  {
    tau_limit_ = limit_multiplier_ * tu.segment(startParam_, n);
  }
  else
  {
    tau_limit_ = Eigen::VectorXd::Constant(n, 100.0);
  }
}

void PeakTorqueFunction::customTorqueLimits(const Eigen::VectorXd & limits)
{
  const int n = robot_.mb().nrDof() - startParam_;
  if (limits.size() == n)
  {
    tau_limit_ = limit_multiplier_ * limits;
  }
  else if (limits.size() == robot_.mb().nrDof())
  {
    tau_limit_ = limit_multiplier_ * limits.tail(n);
  }
}

void PeakTorqueFunction::updateb()
{
  int idx = read_idx_->load(std::memory_order_acquire);
  tau_r_peak_ = (*buffers_)[idx].tau_r_peak;

  // CBF formulation: J_v * alphaD_act <= gamma * (tau_limit - T_max)
  // LinearFunction: J * alphaD + b <= 0  <=>  J * alphaD <= -b
  // Therefore, -b = gamma * (tau_limit - T_max)  =>  b = -gamma * (tau_limit - T_max)
  cbf_margin_ = tau_limit_ - tau_r_peak_;
  rhs_ = gamma_ * cbf_margin_;
  b_ = -rhs_;
}

void PeakTorqueFunction::updateJacobian()
{
  int idx = read_idx_->load(std::memory_order_acquire);
  J_v_ = (*buffers_)[idx].J_v;

  const int n = robot_.mb().nrDof() - startParam_;
  const int N = robot_.mb().nrDof();

  if (J_v_.rows() == n && J_v_.cols() == N)
  {
    // Full pre-impact generalized velocity sensitivity J_v in R^(n x N)
    J_full_ = J_v_;
  }
  else if (J_v_.rows() == n && J_v_.cols() == n)
  {
    // Fallback: actuated-only sensitivity padded with zeros for floating base
    J_full_.leftCols(startParam_).setZero();
    J_full_.rightCols(n) = J_v_;
  }
  else
  {
    J_full_.setZero();
  }

  // Bind full Jacobian into TVM linear function matrix for alphaD
  jacobian_[robot_.tvmRobot().alphaD().get()] = J_full_;
}

} // namespace mc_tvm

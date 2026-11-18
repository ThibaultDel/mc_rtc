/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_solver/PeakTorqueConstraint.h>
#include <mc_solver/TVMQPSolver.h>
#include <mc_tvm/RobotFrame.h>
#include <mc_rtc/logging.h>
#include <RBDyn/FD.h>
#include <RBDyn/Jacobian.h>

#include "TVMPeakTorqueConstraint.h"
#include "PeakTorqueSensitivityInternal.h"

#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <chrono>

namespace mc_solver
{

struct WorkerInput
{
  Eigen::VectorXd v_pre;
  Eigen::VectorXd q0;
  Eigen::MatrixXd M_robot;
  Eigen::MatrixXd J_c;
  Eigen::VectorXd J_n;
  bool has_data = false;
};

struct PeakTorqueConstraint::Impl
{
  const mc_rbdyn::Robot & robot_;
  std::string impact_frame_name_;
  Eigen::Vector3d impact_normal_;
  PeakTorqueConfig config_;

  int total_dofs_ = 0;
  int actuated_dofs_ = 0;
  int startParam_ = 0;
  Eigen::VectorXd J_m_;
  std::vector<internal::JointFrictionParams> joint_params_;

  // TVM Function and Constraint binding
  mc_tvm::PeakTorqueFunctionPtr fn_;
  std::unique_ptr<TVMPeakTorqueConstraint> tvm_constraint_;

  // Lock-free double buffer
  std::shared_ptr<std::array<mc_tvm::PeakTorqueBuffer, 2>> buffers_;
  std::shared_ptr<std::atomic<int>> read_idx_;

  // Background worker thread
  std::thread worker_thread_;
  std::atomic<bool> worker_running_{false};
  std::atomic<double> last_solve_time_ms_{0.0};

  WorkerInput input_snapshot_;
  std::mutex input_mutex_;
  std::condition_variable cv_;
  std::mutex cv_mutex_;

  Impl(const mc_rbdyn::Robot & robot,
       const std::string & impact_frame_name,
       const Eigen::Vector3d & impact_normal,
       const PeakTorqueConfig & config)
    : robot_(robot),
      impact_frame_name_(impact_frame_name),
      impact_normal_(impact_normal.normalized()),
      config_(config)
  {
    total_dofs_ = static_cast<int>(robot_.mb().nrDof());
    startParam_ = robot_.tvmRobot().qFloatingBase() ? robot_.tvmRobot().qFloatingBase()->space().tSize() : 0;
    actuated_dofs_ = total_dofs_ - startParam_;

    // Initialize double-buffers
    buffers_ = std::make_shared<std::array<mc_tvm::PeakTorqueBuffer, 2>>();
    read_idx_ = std::make_shared<std::atomic<int>>(0);

    for (int b = 0; b < 2; ++b)
    {
      (*buffers_)[b].J_v = Eigen::MatrixXd::Zero(actuated_dofs_, total_dofs_);
      (*buffers_)[b].tau_r_peak = Eigen::VectorXd::Zero(actuated_dofs_);
      (*buffers_)[b].seq = 0;
    }

    // Reflected rotor inertias J_m: use user-specified, or default fallback (0.05 kg*m^2)
    if (config_.rotor_inertias.size() == actuated_dofs_)
    {
      J_m_ = config_.rotor_inertias;
    }
    else if (config_.rotor_inertias.size() == total_dofs_)
    {
      J_m_ = config_.rotor_inertias.tail(actuated_dofs_);
    }
    else
    {
      J_m_ = Eigen::VectorXd::Constant(actuated_dofs_, 0.05);
    }

    // Initialize friction parameters per joint
    joint_params_.resize(actuated_dofs_);

    // Map actuated DOFs to joint names and assign per-joint or fallback friction params
    for (size_t jIdx = 0; jIdx < robot_.mb().joints().size(); ++jIdx)
    {
      int dofPos = robot_.mb().jointPosInDof(static_cast<int>(jIdx));
      int actDof = dofPos - startParam_;
      const auto & j = robot_.mb().joint(static_cast<int>(jIdx));
      if (actDof >= 0 && j.dof() > 0)
      {
        JointFrictionConfig f_cfg;
        if (config_.joint_friction.count(j.name()))
        {
          f_cfg = config_.joint_friction.at(j.name());
        }
        else
        {
          f_cfg.K = config_.K;
          f_cfg.B = config_.B;
          f_cfg.F_S = config_.F_S;
          f_cfg.F_C = config_.F_C;
          f_cfg.v_S = config_.v_S;
          f_cfg.D = config_.D;
        }

        for (int d = 0; d < j.dof(); ++d)
        {
          if (actDof + d < actuated_dofs_)
          {
            joint_params_[actDof + d].K = f_cfg.K;
            joint_params_[actDof + d].B = f_cfg.B;
            joint_params_[actDof + d].F_S = f_cfg.F_S;
            joint_params_[actDof + d].F_C = f_cfg.F_C;
            joint_params_[actDof + d].v_S = f_cfg.v_S;
            joint_params_[actDof + d].D = f_cfg.D;
            joint_params_[actDof + d].initialize(0.001);
          }
        }
      }
    }

    // Create TVM LinearFunction
    fn_ = std::make_shared<mc_tvm::PeakTorqueFunction>(
        robot_, config_.gamma, config_.limit_multiplier, buffers_, read_idx_);

    if (config_.custom_torque_limits.size() > 0)
    {
      fn_->customTorqueLimits(config_.custom_torque_limits);
    }

    tvm_constraint_ = std::make_unique<TVMPeakTorqueConstraint>(fn_);

    startWorker();
  }

  ~Impl()
  {
    stopWorker();
  }

  void startWorker()
  {
    if (worker_running_) return;
    worker_running_ = true;
    worker_thread_ = std::thread(&Impl::workerLoop, this);
    mc_rtc::log::info("[PeakTorqueConstraint] Started asynchronous sensitivity worker thread (target rate: {} Hz).",
                      config_.worker_rate_hz);
  }

  void stopWorker()
  {
    if (!worker_running_) return;
    worker_running_ = false;
    cv_.notify_all();
    if (worker_thread_.joinable())
    {
      worker_thread_.join();
    }
    mc_rtc::log::info("[PeakTorqueConstraint] Stopped asynchronous sensitivity worker thread.");
  }

  void workerLoop()
  {
    internal::ModelReduction reducer;

    while (worker_running_)
    {
      auto t_start = std::chrono::steady_clock::now();

      WorkerInput in;
      {
        std::lock_guard<std::mutex> lock(input_mutex_);
        in = input_snapshot_;
      }

      if (in.has_data && in.v_pre.size() == total_dofs_ && in.M_robot.rows() == total_dofs_)
      {
        try
        {
          auto t_solve_start = std::chrono::high_resolution_clock::now();

          // 1. Contact nullspace reduction
          internal::ReductionResults red;
          if (reducer.reduce(in.M_robot, J_m_, in.J_c, total_dofs_, actuated_dofs_, red))
          {
            // 2. Multi-rate dynamic impulse shock sensitivity
            internal::PeakTorqueSensitivitySolver solver(red, J_m_, joint_params_, config_.delta_t, config_.c_res);
            auto sens = solver.solve(in.v_pre, in.J_n, config_.horizon);

            if (sens.success)
            {
              auto t_solve_end = std::chrono::high_resolution_clock::now();
              double ms = std::chrono::duration<double, std::milli>(t_solve_end - t_solve_start).count();
              last_solve_time_ms_.store(ms, std::memory_order_relaxed);

              // 3. Publish to inactive double-buffer
              int cur_read = read_idx_->load(std::memory_order_relaxed);
              int next_write = 1 - cur_read;
              (*buffers_)[next_write].J_v = sens.dTauMax_dvpre;
              (*buffers_)[next_write].tau_r_peak = sens.tau_r_peak;
              (*buffers_)[next_write].seq++;
              read_idx_->store(next_write, std::memory_order_release);
            }
          }
        }
        catch (const std::exception & e)
        {
          mc_rtc::log::warning("[PeakTorqueConstraint] Worker exception: {}", e.what());
        }
      }

      auto t_elapsed = std::chrono::steady_clock::now() - t_start;
      auto target_period = std::chrono::duration<double>(1.0 / std::max(1.0, config_.worker_rate_hz));
      auto sleep_dur = target_period - t_elapsed;
      if (sleep_dur > std::chrono::milliseconds(1))
      {
        std::unique_lock<std::mutex> lock(cv_mutex_);
        cv_.wait_for(lock, sleep_dur, [this]() { return !worker_running_; });
      }
    }
  }

  void updateInput(mc_solver::QPSolver & solver)
  {
    const auto & mb = robot_.mb();
    const auto & mbc = robot_.mbc();

    WorkerInput in;
    in.v_pre = Eigen::VectorXd::Zero(total_dofs_);
    in.q0 = Eigen::VectorXd::Zero(actuated_dofs_);

    // 1. Generalized velocity v_pre
    int dof_idx = 0;
    for (size_t i = 0; i < mb.joints().size(); ++i)
    {
      int dof = mb.joints()[i].dof();
      if (dof == 6)
      {
        // Floating base velocity: [angular, linear]
        const auto & b_vel = mbc.alpha[i];
        if (b_vel.size() == 6)
        {
          in.v_pre.segment(dof_idx, 6) = Eigen::Vector6d(b_vel.data());
        }
      }
      else if (dof > 0)
      {
        for (int d = 0; d < dof; ++d)
        {
          if (mbc.alpha[i].size() > static_cast<size_t>(d))
          {
            in.v_pre(dof_idx + d) = mbc.alpha[i][d];
          }
        }
      }
      dof_idx += dof;
    }

    // 2. Actuated joint configurations q0
    for (size_t jIdx = 0; jIdx < mb.joints().size(); ++jIdx)
    {
      int dofPos = mb.jointPosInDof(static_cast<int>(jIdx));
      int actDof = dofPos - startParam_;
      int dof = mb.joint(static_cast<int>(jIdx)).dof();
      if (actDof >= 0 && dof > 0 && mbc.q[jIdx].size() > 0)
      {
        for (int d = 0; d < dof; ++d)
        {
          if (actDof + d < actuated_dofs_ && static_cast<size_t>(d) < mbc.q[jIdx].size())
          {
            in.q0(actDof + d) = mbc.q[jIdx][d];
          }
        }
      }
    }

    // 3. Mass matrix from ForwardDynamics
    rbd::ForwardDynamics fd(robot_.mb());
    fd.computeH(robot_.mb(), robot_.mbc());
    in.M_robot = fd.H();

    // 4. Bilateral support contact surfaces
    std::vector<std::string> surfaces = config_.contact_surfaces;
    if (surfaces.empty())
    {
      for (const auto & c : solver.contacts())
      {
        if (c.r1Index() == robot_.robotIndex() && c.r1Surface())
        {
          surfaces.push_back(c.r1Surface()->name());
        }
        else if (c.r2Index() == robot_.robotIndex() && c.r2Surface())
        {
          surfaces.push_back(c.r2Surface()->name());
        }
      }
    }

    if (!surfaces.empty())
    {
      in.J_c = Eigen::MatrixXd::Zero(6 * surfaces.size(), total_dofs_);
      for (size_t i = 0; i < surfaces.size(); ++i)
      {
        if (robot_.hasSurface(surfaces[i]))
        {
          const auto & surf = robot_.surface(surfaces[i]);
          rbd::Jacobian c_jac(mb, surf.bodyName());
          Eigen::MatrixXd J_full(6, total_dofs_);
          c_jac.fullJacobian(mb, c_jac.jacobian(mb, mbc), J_full);
          in.J_c.middleRows(i * 6, 6) = J_full;
        }
      }
    }

    // 5. Impact normal Jacobian J_n
    if (robot_.hasFrame(impact_frame_name_))
    {
      const auto & frame = robot_.frame(impact_frame_name_);
      auto jac = frame.tvm_frame().rbdJacobian();
      Eigen::MatrixXd J_h(6, total_dofs_);
      jac.fullJacobian(mb, jac.jacobian(mb, mbc), J_h);
      Eigen::MatrixXd J_lin = J_h.bottomRows(3);
      in.J_n = impact_normal_.transpose() * J_lin;
    }
    else
    {
      in.J_n = Eigen::VectorXd::Zero(total_dofs_);
    }

    in.has_data = true;

    {
      std::lock_guard<std::mutex> lock(input_mutex_);
      input_snapshot_ = in;
    }
    cv_.notify_one();
  }
};

PeakTorqueConstraint::PeakTorqueConstraint(
    const mc_rbdyn::Robot & robot,
    const mc_rbdyn::RobotFrame & impact_frame,
    const Eigen::Vector3d & impact_normal,
    const PeakTorqueConfig & config)
  : impl_(std::make_unique<Impl>(robot, impact_frame.name(), impact_normal, config))
{
}

PeakTorqueConstraint::PeakTorqueConstraint(
    const mc_rbdyn::Robot & robot,
    const std::string & impact_frame_name,
    const Eigen::Vector3d & impact_normal,
    const PeakTorqueConfig & config)
  : impl_(std::make_unique<Impl>(robot, impact_frame_name, impact_normal, config))
{
}

PeakTorqueConstraint::~PeakTorqueConstraint() = default;

void PeakTorqueConstraint::addToSolverImpl(mc_solver::QPSolver & solver)
{
  switch (solver.backend())
  {
    case QPSolver::Backend::TVM:
      if (impl_->tvm_constraint_)
      {
        impl_->tvm_constraint_->addToSolver(tvm_solver(solver));
      }
      break;
    default:
      mc_rtc::log::error_and_throw<std::runtime_error>(
          "[PeakTorqueConstraint] PeakTorqueConstraint is currently only implemented for the TVM solver backend.");
  }
}

void PeakTorqueConstraint::removeFromSolverImpl(mc_solver::QPSolver & solver)
{
  switch (solver.backend())
  {
    case QPSolver::Backend::TVM:
      if (impl_->tvm_constraint_)
      {
        impl_->tvm_constraint_->removeFromSolver(tvm_solver(solver));
      }
      break;
    default:
      break;
  }
}

void PeakTorqueConstraint::update(mc_solver::QPSolver & solver)
{
  impl_->updateInput(solver);
}

void PeakTorqueConstraint::startWorker()
{
  impl_->startWorker();
}

void PeakTorqueConstraint::stopWorker()
{
  impl_->stopWorker();
}

const mc_tvm::PeakTorqueFunctionPtr & PeakTorqueConstraint::function() const
{
  return impl_->fn_;
}

const Eigen::VectorXd & PeakTorqueConstraint::predictedPeakTorques() const
{
  return impl_->fn_->predictedPeakTorques();
}

const Eigen::VectorXd & PeakTorqueConstraint::torqueLimits() const
{
  return impl_->fn_->torqueLimits();
}

const Eigen::VectorXd & PeakTorqueConstraint::cbfMargins() const
{
  return impl_->fn_->cbfMargins();
}

const Eigen::MatrixXd & PeakTorqueConstraint::jacobianV() const
{
  return impl_->fn_->jacobianV();
}

double PeakTorqueConstraint::workerRateHz() const
{
  return impl_->config_.worker_rate_hz;
}

double PeakTorqueConstraint::lastSolveTimeMs() const
{
  return impl_->last_solve_time_ms_.load(std::memory_order_relaxed);
}

uint64_t PeakTorqueConstraint::updateSeq() const
{
  int idx = impl_->read_idx_->load(std::memory_order_relaxed);
  return (*impl_->buffers_)[idx].seq;
}

void PeakTorqueConstraint::gamma(double g)
{
  impl_->config_.gamma = g;
  impl_->fn_->gamma(g);
}

double PeakTorqueConstraint::gamma() const
{
  return impl_->config_.gamma;
}

void PeakTorqueConstraint::limitMultiplier(double m)
{
  impl_->config_.limit_multiplier = m;
  impl_->fn_->limitMultiplier(m);
}

double PeakTorqueConstraint::limitMultiplier() const
{
  return impl_->config_.limit_multiplier;
}

void PeakTorqueConstraint::impactNormal(const Eigen::Vector3d & n)
{
  impl_->impact_normal_ = n.normalized();
}

const Eigen::Vector3d & PeakTorqueConstraint::impactNormal() const
{
  return impl_->impact_normal_;
}

} // namespace mc_solver

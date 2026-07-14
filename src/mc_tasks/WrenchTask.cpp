/*
 * Copyright 2015-2026 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tasks/WrenchTask.h>

#include <mc_tasks/MetaTaskLoader.h>
#include <mc_tasks/TrajectoryTaskGeneric.h>

#include <mc_tvm/WrenchFunction.h>

#include <mc_tasks/MetaTaskLoader.h>

// #include <mc_solver/TasksQPSolver.h>

#include <mc_rbdyn/rpy_utils.h>

#include <mc_rbdyn/hat.h>
#include <mc_rtc/ConfigurationHelpers.h>
#include <mc_rtc/deprecated.h>
#include <mc_rtc/gui/Force.h>
#include <mc_rtc/gui/NumberInput.h>
#include <SpaceVecAlg/SpaceVecAlg>
#include "mc_tasks/MetaTask.h"

namespace mc_tasks
{

static inline mc_rtc::void_ptr_caster<mc_tvm::WrenchFunction> tvm_error{};

WrenchTask::WrenchTask(const mc_rbdyn::RobotFrame & frame, double weight)
: TrajectoryTaskGeneric(frame.robot().robots(), frame.robot().robotIndex(), 0, weight), frame_(frame)
{
  finalize<Backend::TVM, mc_tvm::WrenchFunction>(frame);
  type_ = "wrench";
  name_ = "wrench_" + frame.robot().name() + "_" + frame.name();
  isNoneTaskDynamics_ = true;
}

WrenchTask::WrenchTask(const std::string & surfaceName,
                       const mc_rbdyn::Robots & robots,
                       unsigned int robotIndex,
                       double weight)
: WrenchTask(robots.robot(robotIndex).frame(surfaceName), weight)
{
}

void WrenchTask::reset()
{
  TrajectoryTaskGeneric::reset();
  switch(backend_)
  {
    // case Backend::Tasks:
    //   tasks_error(errorT)->target(frame_->position());
    //   break;
    case Backend::TVM:
      tvm_error(errorT)->reset();
      break;
    default:
      break;
  }
}

sva::ForceVecd WrenchTask::targetWrench() const
{
  switch(backend_)
  {
    // case Backend::Tasks:
    //   return tasks_error(errorT)->target();
    case Backend::TVM:
      return tvm_error(errorT)->targetWrench();
    default:
      mc_rtc::log::error_and_throw("Not implemented");
  }
}

void WrenchTask::targetWrench(const sva::ForceVecd & worldWrench)
{
  switch(backend_)
  {
    // case Backend::Tasks:
    //   tasks_error(errorT)->target(worldwrench);
    //   break;
    case Backend::TVM:
      tvm_error(errorT)->targetWrench(worldWrench);
      break;
    default:
      mc_rtc::log::error("Not implemented");
      break;
  }
}

sva::ForceVecd WrenchTask::currentWrench()
{
  switch(backend_)
  {
    case Backend::TVM:
      return tvm_error(errorT)->currentWrench();
    default:
      mc_rtc::log::error("Not implemented");
      return sva::ForceVecd::Zero();
  }
}

Eigen::MatrixXd WrenchTask::dynamicJacobianTranspose() const noexcept
{
  switch(backend_)
  {
    case Backend::TVM:
      return tvm_error(errorT)->dynamicJacobianTranspose();
    default:
      mc_rtc::log::error("Not implemented");
      return Eigen::MatrixXd();
  }
}

Eigen::Matrix6d WrenchTask::cartesianInertia() const noexcept
{
  switch(backend_)
  {
    case Backend::TVM:
      return tvm_error(errorT)->cartesianInertia();
    default:
      mc_rtc::log::error("Not implemented");
      return Eigen::Matrix6d::Zero();
  }
}

void WrenchTask::addToLogger(mc_rtc::Logger & logger)
{
  logger.addLogEntry(name_ + "_wrench_current", this, [this]() { return currentWrench(); });
  logger.addLogEntry(name_ + "_wrench_target", this, [this]() { return targetWrench(); });
  logger.addLogEntry(name_ + "_eval", this, [this]() { return eval(); });
  logger.addLogEntry(name_ + "_speed", this, [this]() { return speed(); });
  logger.addLogEntry(name_ + "_weight", this, [this]() { return weight(); });
  logger.addLogEntry(name_ + "_dimWeight", this, [this]() { return dimWeight(); });
}

void WrenchTask::addToGUI(mc_rtc::gui::StateBuilder & gui)
{
  MetaTask::addToGUI(gui);
  auto fConf_wrench_target = mc_rtc::gui::ForceConfig();
  fConf_wrench_target.color = mc_rtc::gui::Color::Blue;
  fConf_wrench_target.force_scale = 0.01;

  auto fConf_wrench = mc_rtc::gui::ForceConfig();
  fConf_wrench.color = mc_rtc::gui::Color::Yellow;
  fConf_wrench.force_scale = 0.01;
  gui.addElement({"Tasks", name_, "Wrench"},
                 mc_rtc::gui::Force(
                     "wrench_target", fConf_wrench_target, [this]() { return this->targetWrench(); },
                     [this]() { return this->frame_->position(); }),
                 mc_rtc::gui::Force(
                     "wrench", fConf_wrench, [this]() { return this->currentWrench(); },
                     [this]() { return this->frame_->position(); }));
  gui.addElement({"Tasks", name_, "Gains"},
                 mc_rtc::gui::NumberInput(
                     "weight", [this]() { return this->weight(); }, [this](const double & w) { this->weight(w); }));
}

} // namespace mc_tasks

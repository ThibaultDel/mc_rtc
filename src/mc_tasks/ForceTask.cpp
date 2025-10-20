/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tasks/ForceTask.h>

#include <mc_tasks/MetaTaskLoader.h>
#include <mc_tvm/ForceFunction.h>

#include <mc_solver/TasksQPSolver.h>

#include <mc_rbdyn/rpy_utils.h>

#include <mc_rbdyn/hat.h>
#include <mc_rtc/ConfigurationHelpers.h>
#include <mc_rtc/deprecated.h>
#include <mc_rtc/gui/Force.h>
#include <mc_rtc/gui/NumberInput.h>

namespace mc_tasks
{
static inline mc_rtc::void_ptr_caster<mc_tvm::ForceFunction> tvm_error{};
namespace details
{
inline static mc_rtc::void_ptr_caster<mc_tvm::ForceFunction> tvm_error{};
struct TVMForceTask : public TrajectoryTaskGeneric
{
  TVMForceTask(const mc_rbdyn::Robots & robots,
               const mc_rbdyn::RobotFrame & frame,
               unsigned int robotIndex,
               double weight,
               bool compensateExternalForces = false)
  : TrajectoryTaskGeneric(robots, robotIndex, 0, weight)
  {
    finalize<Backend::TVM, mc_tvm::ForceFunction>(robots.robot(robotIndex), frame, compensateExternalForces);
    type_ = "force";
    name_ = std::string("force_") + robots.robot(robotIndex).name();
    isNoneTaskDynamics_ = true;
  }

  void compensateExternalForces(bool compensate) { tvm_error(errorT)->compensateExternalForces(compensate); }

  bool isCompensatingExternalForces() const { return tvm_error(errorT)->isCompensatingExternalForces(); }

  void update(mc_solver::QPSolver & solver) override { TrajectoryTaskGeneric::update(solver); }

  void force(const sva::ForceVecd & p) { tvm_error(errorT)->force(p); }
};
} // namespace details

ForceTask::ForceTask(const std::string & frameName,
                     const mc_rbdyn::Robots & robots,
                     unsigned int robotIndex,
                     double weight,
                     bool compensateExternalForces)
: ForceTask(robots.robot(robotIndex).frame(frameName), weight, compensateExternalForces)
{
}

void ForceTask::reset()
{
  TrajectoryTaskGeneric::reset();
  switch(backend_)
  {
    case Backend::TVM:
      tvm_error(errorT)->reset();
      break;
    default:
      break;
  }
}

/*! \brief Load parameters from a Configuration object */
void ForceTask::load(mc_solver::QPSolver & solver, const mc_rtc::Configuration & config)
{
  if(config.has("weight")) { weight(config("weight")); }
  if(config.has("compensateExternalForces")) { compensateExternalForces(config("compensateExternalForces")); }
  TrajectoryBase::load(solver, config);
}

sva::ForceVecd ForceTask::target() const
{
  switch(backend_)
  {
    case Backend::TVM:
      return tvm_error(errorT)->force();
    default:
      mc_rtc::log::error_and_throw("Not implemented");
  }
}

void ForceTask::target(const sva::ForceVecd & force)
{
  switch(backend_)
  {
    case Backend::TVM:
      tvm_error(errorT)->force(force);
      break;
    default:
      break;
  }
}

void ForceTask::compensateExternalForces(bool compensate)
{
  switch(backend_)
  {
    case Backend::TVM:
      tvm_error(errorT)->compensateExternalForces(compensate);
      break;
    default:
      mc_rtc::log::error_and_throw("Compensating external forces is only supported in TVM backend");
  }
}

bool ForceTask::isCompensatingExternalForces() const
{
  switch(backend_)
  {
    case Backend::TVM:
      return tvm_error(errorT)->isCompensatingExternalForces();
    default:
      mc_rtc::log::error_and_throw("Compensating external forces is only supported in TVM backend");
  }
}

void ForceTask::addToLogger(mc_rtc::Logger & logger)
{
  TrajectoryBase::addToLogger(logger);
  logger.addLogEntry(name_ + "_force", this, [this]() { return frame_->wrench(); });
  logger.addLogEntry(name_ + "_target_force", this, [this]() { return target(); });
}

// std::function<bool(const mc_tasks::MetaTask &, std::string &)> ForceTask::buildCompletionCriteria(
//     double dt,
//     const mc_rtc::Configuration & config) const
// {
//   if(config.has("wrench"))
//   {
//     if(!frame_->hasForceSensor())
//     {
//       mc_rtc::log::error_and_throw<std::invalid_argument>("[{}] Attempted to use \"wrench\" as completion criteria
//       but "
//                                                           "frame \"{}\" is not attached to a force sensor",
//                                                           name(), frame_->name());
//     }
//     sva::ForceVecd target_w = config("wrench");
//     Eigen::Vector6d target = target_w.vector();
//     Eigen::Vector6d dof = Eigen::Vector6d::Ones();
//     for(int i = 0; i < 6; ++i)
//     {
//       if(std::isnan(target(i)))
//       {
//         dof(i) = 0.;
//         target(i) = 0.;
//       }
//       else if(target(i) < 0) { dof(i) = -1.; }
//     }
//     return [dof, target](const mc_tasks::MetaTask & t, std::string & out)
//     {
//       const auto & self = static_cast<const mc_tasks::ForceTask &>(t);
//       Eigen::Vector6d w = self.robots.robot(self.rIndex).surfaceWrench(self.surface()).vector();
//       for(int i = 0; i < 6; ++i)
//       {
//         if(dof(i) * fabs(w(i)) < target(i)) { return false; }
//       }
//       out += "wrench";
//       return true;
//     };
//   }
//   return MetaTask::buildCompletionCriteria(dt, config);
// }

void ForceTask::addToGUI(mc_rtc::gui::StateBuilder & gui)
{
  TrajectoryTaskGeneric::addToGUI(gui);
  gui.addElement({"Tasks", name_}, mc_rtc::gui::Force(
                                       "Force Target", [this]() { return this->wrench(); },
                                       [this]() -> sva::PTransformd { return frame_->position(); }));
  gui.addElement({"Tasks", name_, "Gains"},
                 mc_rtc::gui::NumberInput(
                     "weight", [this]() { return this->weight(); }, [this](const double & w) { this->weight(w); }));
}

} // namespace mc_tasks

namespace
{

static mc_tasks::MetaTaskPtr loadForceTask(mc_solver::QPSolver & solver, const mc_rtc::Configuration & config)
{
  const auto robotIndex = robotIndexFromConfig(config, solver.robots(), "transform");
  const auto & robot = solver.robots().robot(robotIndex);
  const auto & frame = [&]() -> const mc_rbdyn::RobotFrame &
  {
    if(config.has("surface"))
    {
      mc_rtc::log::deprecated("ForceTask", "surface", "frame");
      return robot.frame(config("surface"));
    }
    else { return robot.frame(config("frame")); }
  }();
  auto t = std::make_shared<mc_tasks::ForceTask>(frame);
  t->load(solver, config);
  return t;
}

static auto reg_dep = mc_tasks::MetaTaskLoader::register_load_function(
    "surfaceTransform",
    [](mc_solver::QPSolver & solver, const mc_rtc::Configuration & config)
    {
      mc_rtc::log::deprecated("TaskLoading", "surfaceTransform", "transform");
      return loadForceTask(solver, config);
    });
static auto reg = mc_tasks::MetaTaskLoader::register_load_function("transform", &loadForceTask);

} // namespace

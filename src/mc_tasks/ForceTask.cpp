/*
 * Copyright 2015-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tasks/ForceTask.h>

#include <mc_tasks/MetaTaskLoader.h>
#include <mc_tasks/TrajectoryTaskGeneric.h>

#include <mc_solver/TVMQPSolver.h>
// #include <mc_solver/TasksQPSolver.h>

#include <mc_tvm/ForceFunction.h>
#include <mc_tvm/Robot.h>

#include <mc_rbdyn/configuration_io.h>

#include <mc_rtc/gui/Force.h>
#include <mc_rtc/gui/NumberInput.h>
#include <vector>

namespace mc_tasks
{

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

inline static mc_rtc::void_ptr_caster<details::TVMForceTask> tvm_error{};

inline static mc_rtc::void_ptr make_error(MetaTask::Backend backend,
                                          const mc_solver::QPSolver & solver,
                                          const mc_rbdyn::RobotFrame & frame,
                                          unsigned int rIndex,
                                          double weight)
{
  switch(backend)
  {
    // case MetaTask::Backend::Tasks:
    //   return mc_rtc::make_void_ptr<tasks::qp::ForceTask>(solver.robots().mbs(), static_cast<int>(rIndex),
    //                                                        solver.robot(rIndex).mbc().tau, weight);
    case MetaTask::Backend::TVM:
      return mc_rtc::make_void_ptr<details::TVMForceTask>(solver.robots(), frame, rIndex, weight);
    default:
      mc_rtc::log::error_and_throw("[ForceTask] Not implemented for solver backend: {}", backend);
  }
}

ForceTask::ForceTask(const mc_solver::QPSolver & solver,
                     const std::string & bodyName,
                     const mc_rbdyn::Robots & robots,
                     unsigned int robotIndex,
                     double weight,
                     bool compensateExternalForces)
: ForceTask(solver, robots.robot(robotIndex).frame(bodyName), weight, compensateExternalForces)
{
}

ForceTask::ForceTask(const mc_solver::QPSolver & solver,
                     const mc_rbdyn::RobotFrame & frame,
                     double weight,
                     bool compensateExternalForces)
: robots_(frame.robot().robots()), frame_(frame),
  pt_(make_error(backend_, solver, frame, frame.robot().robotIndex(), weight)), dt_(solver.dt())
{
  compensateExternalForces_ = compensateExternalForces;
  eval_ = this->eval();
  speed_ = Eigen::VectorXd::Zero(eval_.size());
  curForce = frame_.wrench();

  type_ = "force6d";
  name_ = "force6d_" + frame.robot().name() + "_" + frame.name();
  name(name_);
}

void ForceTask::update(mc_solver::QPSolver & solver)
{
  switch(backend_)
  {
    // case Backend::Tasks:
    // {
    //   const auto & pt = *tasks_error(pt_);
    //   speed_ = pt.dimWeight().asDiagonal() * (pt.eval() - eval_) / dt_;
    //   eval_ = pt.eval();
    //   break;
    // }
    case Backend::TVM:
    {
      auto & pt = *tvm_error(pt_);
      pt.update(solver);
      speed_ = (pt.eval() - eval_) / dt_;
      eval_ = pt.dimWeight().asDiagonal() * pt.eval();
      break;
    }
    default:
      break;
  }
}

void ForceTask::reset()
{
  curForce = frame_.wrench();
}

void ForceTask::removeFromSolver(mc_solver::QPSolver & solver)
{
  if(!inSolver_) { return; }
  inSolver_ = false;
  switch(backend_)
  {
    // case Backend::Tasks:
    //   tasks_solver(solver).removeTask(tasks_error(pt_));
    //   break;
    case Backend::TVM:
      MetaTask::removeFromSolver(*tvm_error(pt_), solver);
      break;
    default:
      break;
  }
}

void ForceTask::addToSolver(mc_solver::QPSolver & solver)
{
  if(inSolver_) { return; }
  inSolver_ = true;
  switch(backend_)
  {
    // case Backend::Tasks:
    //   tasks_solver(solver).addTask(tasks_error(pt_));
    //   break;
    case Backend::TVM:
      MetaTask::addToSolver(*tvm_error(pt_), solver);
      break;
    default:
      break;
  }
}

void ForceTask::selectActiveJoints(mc_solver::QPSolver & solver,
                                   const std::vector<std::string> & activeJointsName,
                                   const std::map<std::string, std::vector<std::array<int, 2>>> & activeDofs)
{
  switch(backend_)
  {
    case Backend::TVM:
      tvm_error(pt_)->selectActiveJoints(solver, activeJointsName, activeDofs);
      break;
    default:
      mc_rtc::log::error_and_throw("selectActiveJoints not implemented for backend {}", backend_);
  }
}

void ForceTask::selectUnactiveJoints(mc_solver::QPSolver & solver,
                                     const std::vector<std::string> & unactiveJointsName,
                                     const std::map<std::string, std::vector<std::array<int, 2>>> & unactiveDofs)
{
  switch(backend_)
  {
    case Backend::TVM:
      tvm_error(pt_)->selectUnactiveJoints(solver, unactiveJointsName, unactiveDofs);
      break;
    default:
      mc_rtc::log::error_and_throw("selectUnactiveJoints not implemented for backend {}", backend_);
  }
}

void ForceTask::resetJointsSelector(mc_solver::QPSolver & solver)
{
  switch(backend_)
  {
    case Backend::TVM:
      tvm_error(pt_)->resetJointsSelector(solver);
      break;
    default:
      mc_rtc::log::error_and_throw("resetJointsSelector not implemented for backend {}", backend_);
  }
}

void ForceTask::add_force(const sva::ForceVecd & dtr)
{
  set_force(curForce + dtr);
}

void ForceTask::set_force(const sva::ForceVecd & tf)
{
  curForce = tf;
  // std::vector<std::vector<double>> force = {{curForce.couple().x()}, {curForce.couple().y()},
  // {curForce.couple().z()},
  //                                           {curForce.force().x()},  {curForce.force().y()}, {curForce.force().z()}};
  switch(backend_)
  {
    // case Backend::Tasks:
    //   tasks_error(pt_)->force(p);
    //   break;
    case Backend::TVM:
      tvm_error(pt_)->force(tf);
      break;
    default:
      break;
  }
}

sva::ForceVecd ForceTask::get_force()
{
  return curForce;
}

void ForceTask::dimWeight(const Eigen::VectorXd & dimW)
{
  switch(backend_)
  {
    // case Backend::Tasks:
    //   tasks_error(pt_)->dimWeight(dimW);
    //   break;
    case Backend::TVM:
      tvm_error(pt_)->dimWeight(dimW);
      break;
    default:
      break;
  }
}

Eigen::VectorXd ForceTask::dimWeight() const
{
  switch(backend_)
  {
    // case Backend::Tasks:
    //   return tasks_error(pt_)->dimWeight();
    case Backend::TVM:
      return tvm_error(pt_)->dimWeight();
    default:
      mc_rtc::log::error_and_throw("Not implemented");
  }
}

Eigen::VectorXd ForceTask::eval() const
{
  switch(backend_)
  {
    // case Backend::Tasks:
    // {
    //   auto & pt = *tasks_error(pt_);
    //   return pt.dimWeight().asDiagonal() * pt.eval();
    // }
    case Backend::TVM:
      return tvm_error(pt_)->eval();
    default:
      mc_rtc::log::error_and_throw("Not implemented");
  }
}

Eigen::VectorXd ForceTask::speed() const
{
  return speed_;
}

void ForceTask::weight(double w)
{
  switch(backend_)
  {
    // case Backend::Tasks:
    //   tasks_error(pt_)->weight(w);
    //   break;
    case Backend::TVM:
      tvm_error(pt_)->weight(w);
      break;
    default:
      break;
  }
}

double ForceTask::weight() const
{
  switch(backend_)
  {
    // case Backend::Tasks:
    //   return tasks_error(pt_)->weight();
    case Backend::TVM:
      return tvm_error(pt_)->weight();
    default:
      mc_rtc::log::error_and_throw("Not implemented");
  }
}

bool ForceTask::inSolver() const
{
  return inSolver_;
}

void ForceTask::compensateExternalForces(bool compensate)
{
  switch(backend_)
  {
    case Backend::TVM:
      tvm_error(pt_)->compensateExternalForces(compensate);
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
      return tvm_error(pt_)->isCompensatingExternalForces();
    default:
      mc_rtc::log::error_and_throw("Compensating external forces is only supported in TVM backend");
  }
}

void ForceTask::load(mc_solver::QPSolver & solver, const mc_rtc::Configuration & config)
{
  MetaTask::load(solver, config);
  if(config.has("weight")) { weight(config("weight")); }
}

void ForceTask::addToLogger(mc_rtc::Logger & logger)
{
  MC_RTC_LOG_HELPER(name_ + "_target", curForce);
  logger.addLogEntry(name_, this, [this]() { return frame_.position(); });
}

void ForceTask::removeFromLogger(mc_rtc::Logger & logger)
{
  MetaTask::removeFromLogger(logger);
}

void ForceTask::addToGUI(mc_rtc::gui::StateBuilder & gui)
{
  MetaTask::addToGUI(gui);
  gui.addElement({"Tasks", name_}, mc_rtc::gui::Force(
                                       "Force Target", [this]() { return this->get_force(); },
                                       [this]() -> sva::PTransformd { return frame_.position(); }));
  gui.addElement({"Tasks", name_, "Gains"},
                 mc_rtc::gui::NumberInput(
                     "weight", [this]() { return this->weight(); }, [this](const double & w) { this->weight(w); }));
}

void ForceTask::name(const std::string & name)
{
  MetaTask::name(name);
}

} // namespace mc_tasks

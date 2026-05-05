#include <mc_tasks/TorquePDRelativeCartesianTask.h>

#include <mc_rtc/gui/ArrayInput.h>
#include <mc_rtc/gui/ArrayLabel.h>
#include <mc_rtc/gui/NumberInput.h>
#include <mc_rtc/gui/NumberSlider.h>
#include <mc_rtc/gui/Transform.h>
#include <SpaceVecAlg/SpaceVecAlg>
#include "mc_tasks/TorquePDCartesianTask.h"

namespace mc_tasks
{

TorquePDRelativeCartesianTask::TorquePDRelativeCartesianTask(const mc_solver::QPSolver & solver,
                                                             const mc_rbdyn::RobotFrame & frame,
                                                             const mc_rbdyn::Frame & relative,
                                                             double stiffness,
                                                             double weight)
: TorquePDCartesianTask(solver, frame, stiffness, weight), relative_(relative),
  posTarget_rel_(sva::PTransformd::Identity()), velTarget_rel_(sva::MotionVecd::Zero())
{
  if(backend_ == Backend::Tasks)
    mc_rtc::log::error_and_throw<std::runtime_error>(
        "[mc_tasks] Can't use TorquePDRelativeCartesianTask with {} backend, please use TVM backend", backend_);

  name_ = "pd_relative_cartesian_" + frame.name() + "_" + solver.robots().robot(frame.robot().robotIndex()).name();
  type_ = "pd_relative_cartesian";

  setStiffness(stiffness);
  setDamping(2.0 * sqrt(stiffness)); // Critical damping by default

  reset();
}

TorquePDRelativeCartesianTask::TorquePDRelativeCartesianTask(const mc_solver::QPSolver & solver,
                                                             const std::string & bodyName,
                                                             const std::string & relBodyName,
                                                             const mc_rbdyn::Robots & robots,
                                                             unsigned int rIndex,
                                                             double stiffness,
                                                             double weight)
: TorquePDRelativeCartesianTask(solver,
                                robots.robot(rIndex).frame(bodyName),
                                relBodyName.size() ? robots.robot(rIndex).frame(relBodyName)
                                                   : robots.robot(rIndex).frame(0),
                                stiffness,
                                weight)
{
}

TorquePDRelativeCartesianTask::TorquePDRelativeCartesianTask(const mc_solver::QPSolver & solver,
                                                             const std::string & bodyName,
                                                             const mc_rbdyn::Robots & robots,
                                                             unsigned int rIndex,
                                                             double stiffness,
                                                             double weight)
: TorquePDRelativeCartesianTask(solver,
                                robots.robot(rIndex).frame(bodyName),
                                robots.robot(rIndex).frame(0),
                                stiffness,
                                weight)
{
}

void TorquePDRelativeCartesianTask::reset()
{
  sva::PTransformd X_0_body = frame().position();
  sva::PTransformd X_0_rel = relative_->position();
  posTarget_rel_ = X_0_body * (X_0_rel.inv());
  velTarget_rel_.Zero();
}

void TorquePDRelativeCartesianTask::update(mc_solver::QPSolver & solver)
{
  posTarget_ = posTarget_rel_ * relative_->position();
  velTarget_ = relative_->position() * velTarget_rel_ + relative_->velocity();
  TorquePDCartesianTask::update(solver);
}

void TorquePDRelativeCartesianTask::setPosTarget(const sva::PTransformd & xd)
{
  posTarget_rel_ = xd;
}

void TorquePDRelativeCartesianTask::setVelTarget(const sva::MotionVecd & xd_dot)
{
  velTarget_rel_ = xd_dot;
}

sva::PTransformd TorquePDRelativeCartesianTask::posTarget()
{
  return posTarget_rel_;
}

sva::MotionVecd TorquePDRelativeCartesianTask::velTarget()
{
  return velTarget_rel_;
}

void TorquePDRelativeCartesianTask::addToGUI(mc_rtc::gui::StateBuilder & gui)
{
  TorquePDCartesianTask::addToGUI(gui);
  gui.removeElement({"Tasks", name_, "Position Target"}, "pos_target");
  gui.removeElement({"Tasks", name_, "Velocity Target"}, "vel_target");

  gui.addElement({"Tasks", name_, "Position Target"},
                 mc_rtc::gui::Transform(
                     "pos_target", [this]() { return posTarget_rel_ * relative_->position(); },
                     [this](const sva::PTransformd & X_0_target)
                     { posTarget_rel_ = X_0_target * relative_->position().inv(); }),
                 mc_rtc::gui::Transform("rel_pos", [this]() { return relative_->position(); }),
                 mc_rtc::gui::Transform("pos", [this]() { return frame().position(); }));

  gui.addElement({"Tasks", name_, "Velocity Target"},
                 mc_rtc::gui::ArrayInput(
                     "vel_target", [this]() { return velTarget_rel_; },
                     [this](const sva::MotionVecd & xd_dot) { velTarget_rel_ = xd_dot; }));
}

} // namespace mc_tasks

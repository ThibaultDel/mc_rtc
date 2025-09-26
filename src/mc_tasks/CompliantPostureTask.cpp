#include <mc_tasks/CompliantPostureTask.h>

#include <mc_rtc/gui/Checkbox.h>
#include <mc_tvm/Robot.h>
#include "mc_rtc/gui/ArrayInput.h"

namespace mc_tasks
{

CompliantPostureTask::CompliantPostureTask(const mc_solver::QPSolver & solver,
                                           unsigned int rIndex,
                                           double stiffness,
                                           double weight)
: PostureTask(solver, rIndex, stiffness, weight),
  gamma_(Eigen::VectorXd::Zero(solver.robots().robot(rIndex).mb().nrDof())),
  tvm_robot_(solver.robots().robot(rIndex).tvmRobot()),
  refAccel_(Eigen::VectorXd::Zero(solver.robots().robot(rIndex).mb().nrDof()))
{
  if(backend_ == Backend::Tasks)
    mc_rtc::log::error_and_throw<std::runtime_error>(
        "[mc_tasks] Can't use CompliantEndEffectorTask with {} backend, please use TVM or TVMHierarchical backend",
        backend_);
  name_ = std::string("compliant_posture_") + solver.robots().robot(rIndex).name();
  type_ = "compliant_posture";
}

void CompliantPostureTask::refAccel(const Eigen::VectorXd & refAccel) noexcept
{
  refAccel_ = refAccel;
}

void CompliantPostureTask::update(mc_solver::QPSolver & solver)
{
  Eigen::VectorXd disturbance = tvm_robot_.alphaDExternal();
  // mc_rtc::log::info("Ref accel from disturbance : {}", disturbance.transpose());
  Eigen::VectorXd disturbedAccel = refAccel_ + gamma_.asDiagonal() * disturbance;
  PostureTask::refAccel(disturbedAccel);
  PostureTask::update(solver);
}

void CompliantPostureTask::makeCompliant(bool compliance)
{
  if(compliance) { gamma_.setOnes(); }
  else { gamma_.setZero(); }
}

void CompliantPostureTask::makeCompliant(Eigen::VectorXd gamma)
{
  gamma_ = gamma;
}

bool CompliantPostureTask::isCompliant(void)
{
  return gamma_.norm() > 0;
}

void CompliantPostureTask::addToGUI(mc_rtc::gui::StateBuilder & gui)
{
  gui.addElement(
      {"Tasks", name_, "Compliance"},
      mc_rtc::gui::Checkbox(
          "Compliance is active", [this]() { return isCompliant(); }, [this]() { makeCompliant(!isCompliant()); }),
      mc_rtc::gui::ArrayInput("Gamma", {"Joint_1", "Joint_2", "Joint_3", "Joint_4", "Joint_5", "Joint_6", "Joint_7"},
                              gamma_));
  PostureTask::addToGUI(gui);
}

} // namespace mc_tasks

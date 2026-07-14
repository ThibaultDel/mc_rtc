/*
 * Copyright 2015-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tasks/CompliantEndEffectorTask.h>

<<<<<<< HEAD
#include <mc_rtc/gui/Checkbox.h>
#include <mc_tvm/Robot.h>
#include <SpaceVecAlg/EigenTypedef.h>
#include "mc_rtc/gui/ArrayInput.h"
=======
#include <mc_rtc/gui/ArrayInput.h>
#include <mc_rtc/gui/Checkbox.h>
#include <mc_tvm/Robot.h>
#include <SpaceVecAlg/EigenTypedef.h>
#include <Eigen/src/Core/Matrix.h>
>>>>>>> bastien/hrp5p

namespace mc_tasks
{

CompliantEndEffectorTask::CompliantEndEffectorTask(const std::string & bodyName,
                                                   const mc_rbdyn::Robots & robots,
                                                   unsigned int robotIndex,
                                                   double stiffness,
                                                   double weight)
: EndEffectorTask(robots.robot(robotIndex).frame(bodyName), stiffness, weight),
<<<<<<< HEAD
  compliant_matrix_(Eigen::Matrix6d::Zero()), tvm_robot_(nullptr), rIdx_(robotIndex), bodyName_(bodyName),
  refAccel_(Eigen::Vector6d::Zero())
{
  const mc_rbdyn::RobotFrame & frame = robots.robot(robotIndex).frame(bodyName);

  if(backend_ == Backend::Tasks)
    mc_rtc::log::error_and_throw<std::runtime_error>(
        "[mc_tasks] Can't use CompliantEndEffectorTask with {} backend, please use TVM or TVMHierarchical backend",
        backend_);

  type_ = "compliant_body6d";
  name_ = "compliant_body6d_" + frame.robot().name() + "_" + frame.name();
=======
  compliant_matrix_(Eigen::Matrix6d::Zero()), robot_(robots.robot(robotIndex)),
  tvm_robot_(robots.robot(robotIndex).tvmRobot()), rIdx_(robotIndex), bodyName_(bodyName),
  frame_(robots.robot(robotIndex).frame(bodyName)), refAccel_(Eigen::Vector6d::Zero())
{
  switch(backend_)
  {
    case Backend::Tasks:
    case Backend::TVM:
      break;
    default:
      mc_rtc::log::error_and_throw<std::runtime_error>(
          "[mc_tasks] Can't use CompliantEndEffectorTask with {} backend, please use Tasks or TVM backend", backend_);
      break;
  }

  type_ = "compliant_body6d";
  name_ = "compliant_body6d_" + frame_.robot().name() + "_" + frame_.name();
>>>>>>> bastien/hrp5p
  EndEffectorTask::name(name_);
}

void CompliantEndEffectorTask::refAccel(const Eigen::Vector6d & refAccel) noexcept
{
  refAccel_ = refAccel;
}

void CompliantEndEffectorTask::makeCompliant(bool compliance)
{
  if(compliance) { compliant_matrix_.diagonal().setOnes(); }
<<<<<<< HEAD
  else { compliant_matrix_.diagonal().setZero(); }
=======
  else
  {
    compliant_matrix_.diagonal().setZero();
  }
>>>>>>> bastien/hrp5p
}

void CompliantEndEffectorTask::setComplianceVector(Eigen::Vector6d gamma)
{
  compliant_matrix_.diagonal() = gamma;
}

bool CompliantEndEffectorTask::isCompliant(void)
{
  return compliant_matrix_.diagonal().norm() > 0;
}

Eigen::Vector6d CompliantEndEffectorTask::getComplianceVector(void)
{
  return compliant_matrix_.diagonal();
}

void CompliantEndEffectorTask::addToSolver(mc_solver::QPSolver & solver)
{
<<<<<<< HEAD
  tvm_robot_ = &solver.robots().robot(rIdx_).tvmRobot();
  jac_ = new rbd::Jacobian(tvm_robot_->robot().mb(), bodyName_);

  EndEffectorTask::addToSolver(solver);
=======
  EndEffectorTask::addToSolver(solver);
  jac_ = new rbd::Jacobian(robot_.mb(), frame_.body());
>>>>>>> bastien/hrp5p
}

void CompliantEndEffectorTask::update(mc_solver::QPSolver & solver)
{
<<<<<<< HEAD
  auto J = jac_->jacobian(tvm_robot_->robot().mb(), tvm_robot_->robot().mbc());
  Eigen::Vector6d disturbance = J * tvm_robot_->alphaDExternal();
=======
  Eigen::MatrixXd J = jac_->jacobian(solver.robot(rIdx_).mb(), solver.robot(rIdx_).mbc());
  Eigen::VectorXd acc;

  if(backend_ == Backend::Tasks)
  {
    if(robot_.compensationTorquesAcc()) { acc = robot_.compensationTorquesAcc().value(); }
    else
    {
      acc = robot_.externalTorquesAcc();
    }
  }
  else
  {
    if(tvm_robot_.alphaDCompensation()) { acc = tvm_robot_.alphaDCompensation().value(); }
    else
    {
      acc = tvm_robot_.alphaDExternal();
    }
  }
  Eigen::Vector6d disturbance = J * acc;

>>>>>>> bastien/hrp5p
  Eigen::Vector6d disturbedAccel = refAccel_ + compliant_matrix_ * disturbance;

  EndEffectorTask::positionTask->refAccel(disturbedAccel.tail(3));
  EndEffectorTask::orientationTask->refAccel(disturbedAccel.head(3));

  // EndEffectorTask::update(solver);
}

<<<<<<< HEAD
=======
void CompliantEndEffectorTask::load(mc_solver::QPSolver & solver, const mc_rtc::Configuration & config)
{
  MetaTask::load(solver, config);
  if(config.has("stiffness"))
  {
    auto s = config("stiffness");
    if(s.size())
    {
      Eigen::VectorXd stiff = s;
      positionTask->stiffness(stiff);
      orientationTask->stiffness(stiff);
    }
    else
    {
      double stiff = s;
      positionTask->stiffness(stiff);
      orientationTask->stiffness(stiff);
    }
  }
  if(config.has("damping"))
  {
    auto d = config("damping");
    if(d.size())
    {
      positionTask->setGains(positionTask->dimStiffness(), d);
      orientationTask->setGains(orientationTask->dimStiffness(), d);
    }
    else
    {
      positionTask->setGains(positionTask->stiffness(), d);
      orientationTask->setGains(orientationTask->stiffness(), d);
    }
  }
  if(config.has("compliance"))
  {
    auto g = config("compliance");
    if(g.size()) { setComplianceVector(g); }
    else
    {
      makeCompliant((double)g != 0);
    }
  }
  if(config.has("weight"))
  {
    double w = config("weight");
    positionTask->weight(w);
    orientationTask->weight(w);
  }
}

>>>>>>> bastien/hrp5p
void CompliantEndEffectorTask::addToGUI(mc_rtc::gui::StateBuilder & gui)
{
  gui.addElement({"Tasks", name_, "Compliance"}, mc_rtc::gui::Checkbox(
                                                     "Compliance is active", [this]() { return isCompliant(); },
                                                     [this]() { makeCompliant(!isCompliant()); }));
  gui.addElement({"Tasks", name_, "Compliance"},
                 mc_rtc::gui::ArrayInput(
                     "Compliance parameters", {"rx", "ry", "rz", "x", "y", "z"}, [this]()
                     { return getComplianceVector(); }, [this](Eigen::Vector6d v) { setComplianceVector(v); }));

  EndEffectorTask::addToGUI(gui);
}

} // namespace mc_tasks

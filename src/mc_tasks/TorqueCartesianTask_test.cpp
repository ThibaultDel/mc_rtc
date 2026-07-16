#include <mc_tasks/MetaTask.h>
#include <mc_tasks/TorqueCartesianTask_test.h>
#include <mc_tvm/WrenchFunction.h>

#include <mc_rtc/gui/ArrayInput.h>
#include <mc_rtc/gui/ArrayLabel.h>
#include <mc_rtc/gui/Checkbox.h>
#include <mc_rtc/gui/NumberInput.h>
#include <mc_rtc/gui/NumberSlider.h>
#include <mc_rtc/gui/Transform.h>
#include <mc_tvm/Robot.h>
#include <SpaceVecAlg/EigenTypedef.h>
#include <SpaceVecAlg/SpaceVecAlg>
#include "mc_rtc/logging.h"
#include "mc_tasks/WrenchTask.h"
#include <Eigen/src/Core/Matrix.h>

namespace mc_tasks
{

TorqueCartesianTask_test::TorqueCartesianTask_test(const mc_solver::QPSolver & solver,
                                                   const mc_rbdyn::RobotFrame & frame,
                                                   double stiffness,
                                                   double weight)
: WrenchTask(frame, weight), robots_(solver.robots()), rIndex_(frame.robot().robotIndex()),
  nbActuatedJoints(
      (robots_.robot(rIndex_).mb().nrJoints() > 0 && robots_.robot(rIndex_).mb().joint(0).type() == rbd::Joint::Free)
          ? robots_.robot(rIndex_).mb().nrDof() - 6
          : robots_.robot(rIndex_).mb().nrDof()),
  torqueFeedforward_(Eigen::VectorXd::Zero(nbActuatedJoints)),
  torqueGravityCompensation_(Eigen::VectorXd::Zero(nbActuatedJoints))
{
  if(backend_ == Backend::Tasks)
    mc_rtc::log::error_and_throw<std::runtime_error>(
        "[mc_tasks] Can't use TorqueCartesianTask_test with {} backend, please use TVM backend", backend_);

  name_ = "torque_cartesian_test_" + frame.name() + "_" + solver.robots().robot(frame.robot().robotIndex()).name();
  type_ = "torque_cartesian_test";
  // reset();
}

TorqueCartesianTask_test::TorqueCartesianTask_test(const mc_solver::QPSolver & solver,
                                                   const std::string & bodyName,
                                                   unsigned int rIndex,
                                                   double stiffness,
                                                   double weight)
: TorqueCartesianTask_test(solver, solver.robots().robot(rIndex).frame(bodyName), stiffness, weight)
{
}

void TorqueCartesianTask_test::reset()
{
  torqueFeedforward_.setZero();
  WrenchTask::reset();
}

void TorqueCartesianTask_test::update(mc_solver::QPSolver & solver)
{
  Eigen::Vector6d forceComponents = Eigen::Vector6d::Zero();
  Eigen::VectorXd torqueComponents = Eigen::VectorXd::Zero(nbActuatedJoints); // tau_ff + tau_ext_forces + tau_gravity

  auto & tvm_robot = robots_.robot(rIndex_).tvmRobot();
  torqueComponents += torqueFeedforward_;
  torqueGravityCompensation_ = tvm_robot.C();
  if(compensateGravity_) { torqueComponents += torqueGravityCompensation_; }
  forceComponents += WrenchTask::dynamicJacobianTranspose() * torqueComponents;
  WrenchTask::targetWrench(sva::ForceVecd(forceComponents.head<3>(), forceComponents.tail<3>()));
}

void TorqueCartesianTask_test::setTorqueFeedforward(const Eigen::VectorXd & tau_ff)
{
  if(tau_ff.size() != nbActuatedJoints)
  {
    mc_rtc::log::error_and_throw<std::runtime_error>(
        "[TorqueCartesianTask_test] Torque feedforward vector size should be {}, got {}", nbActuatedJoints,
        tau_ff.size());
  }
  torqueFeedforward_ = tau_ff;
}

void TorqueCartesianTask_test::addToGUI(mc_rtc::gui::StateBuilder & gui)
{
  WrenchTask::addToGUI(gui);
}

void TorqueCartesianTask_test::addToLogger(mc_rtc::Logger & logger)
{
  WrenchTask::addToLogger(logger);
}

} // namespace mc_tasks

#include <mc_tasks/TorquePDCartesianTask.h>

#include <RBDyn/Jacobian.h>

#include <mc_rtc/gui/ArrayInput.h>
#include <mc_rtc/gui/ArrayLabel.h>
#include <mc_rtc/gui/Checkbox.h>
#include <mc_rtc/gui/NumberInput.h>
#include <mc_rtc/gui/NumberSlider.h>
#include <mc_rtc/gui/Transform.h>
#include "mc_rtc/logging.h"

namespace mc_tasks
{

TorquePDCartesianTask::TorquePDCartesianTask(const mc_solver::QPSolver & solver,
                                             const mc_rbdyn::RobotFrame & frame,
                                             double stiffness,
                                             double weight)
: TorqueTask(solver, frame.robot().robotIndex(), weight), posTarget_(sva::PTransformd::Identity()),
  velTarget_(sva::MotionVecd::Zero()), robots_(solver.robots()), rIndex_(frame.robot().robotIndex()),
  nbActuatedJoints(
      (robots_.robot(rIndex_).mb().nrJoints() > 0 && robots_.robot(rIndex_).mb().joint(0).type() == rbd::Joint::Free)
          ? robots_.robot(rIndex_).mb().nrDof() - 6
          : robots_.robot(rIndex_).mb().nrDof()),
  torqueFeedforward_(Eigen::VectorXd::Zero(nbActuatedJoints)), stiffness_(Eigen::Vector6d::Zero()),
  damping_(Eigen::Vector6d::Zero()), posError_(sva::MotionVecd::Zero()), velError_(sva::MotionVecd::Zero()),
  torqueTarget_(Eigen::VectorXd::Zero(nbActuatedJoints)), frame_(frame)
{
  if(backend_ == Backend::Tasks)
    mc_rtc::log::error_and_throw<std::runtime_error>(
        "[mc_tasks] Can't use TorquePDCartesianTask with {} backend, please use TVM backend", backend_);

  name_ = "pd_cartesian_" + frame.name() + "_" + solver.robots().robot(frame.robot().robotIndex()).name();
  type_ = "pd_cartesian";

  setStiffness(stiffness);
  setDamping(2.0 * sqrt(stiffness)); // Critical damping by default

  reset();
}

TorquePDCartesianTask::TorquePDCartesianTask(const mc_solver::QPSolver & solver,
                                             const std::string & bodyName,
                                             unsigned int rIndex,
                                             double stiffness,
                                             double weight)
: TorquePDCartesianTask(solver, solver.realRobots().robot(rIndex).frame(bodyName), stiffness, weight)
{
}

void TorquePDCartesianTask::reset()
{
  posTarget_ = frame_->position();
  velTarget_.Zero();
  torqueFeedforward_.setZero();
}

void TorquePDCartesianTask::update(mc_solver::QPSolver & solver)
{

  // auto & realRobot = solver.realRobots().robot(rIndex_);
  // rbd::Jacobian jac_body(realRobot.mb(), frame_->body());
  // Eigen::MatrixXd J_body = jac_body.jacobian(realRobot.mb(), realRobot.mbc());

  // // sva::PTransformd X_0_body = frame_->position();
  // // posError_ = sva::transformError(X_0_body, posTarget_);

  // // sva::MotionVecd vel_0_body = frame_->velocity();
  // // velError_ = velTarget_ - vel_0_body;

  // posError_ = sva::transformVelocity(frame_->position() * posTarget_.inv());
  // sva::MotionVecd V_p_p = jac_body.velocity(frame_->robot().mb(), frame_->robot().mbc(), frame_->X_b_f());
  // velError_ = V_p_p - posError_.cross(V_p_p) - velTarget_;

  // // torqueTarget_ =
  // //     J_body.transpose() * (stiffness_.asDiagonal() * posError_.vector() + damping_.asDiagonal() *
  // velError_.vector())
  // //     + torqueFeedforward_;

  // torqueTarget_ =
  //     J_body.transpose() * (-1 * stiffness_.asDiagonal() * posError_.vector() - damping_.asDiagonal() *
  //     velError_.vector())
  //     + torqueFeedforward_;

  // --- Robot & Jacobian ---
  auto & realRobot = solver.realRobots().robot(rIndex_);
  rbd::Jacobian jac_body(realRobot.mb(), frame_->body());
  Eigen::MatrixXd J = jac_body.jacobian(realRobot.mb(), realRobot.mbc());

  // --- 1. Pose error (LEFT-invariant) ---
  sva::PTransformd X_err = frame_->position().inv() * posTarget_;
  // log map → se(3)
  sva::MotionVecd xi = sva::transformVelocity(X_err);

  // --- 2. Velocity error (GEOMETRICALLY CORRECT) ---
  // transport desired velocity into current frame
  sva::MotionVecd V_d_in_b = frame_->position().inv() * velTarget_;

  // error
  sva::MotionVecd xi_dot = V_d_in_b - frame_->velocity();

  // --- 3. Task-space wrench ---
  Eigen::VectorXd F = stiffness_.asDiagonal() * xi.vector() + damping_.asDiagonal() * xi_dot.vector();

  // --- 4. Joint torques ---
  torqueTarget_ = J.transpose() * F + torqueFeedforward_;

  std::vector<std::vector<double>> torque_vector = realRobot.mbc().jointTorque;
  Eigen::VectorXd torque_target_full = Eigen::VectorXd::Zero(realRobot.mb().nrDof());
  torque_target_full.tail(nbActuatedJoints) = torqueTarget_;

  torque_vector = rbd::sVectorToDof(realRobot.mb(), torque_target_full);

  TorqueTask::torqueTarget(torque_vector);
  TorqueTask::update(solver);
}

void TorquePDCartesianTask::setStiffness(double stiffness)
{
  stiffness_ = Eigen::Vector6d::Constant(stiffness);
}

void TorquePDCartesianTask::setDamping(double damping)
{
  damping_ = Eigen::Vector6d::Constant(damping);
}

void TorquePDCartesianTask::setStiffness(const Eigen::Vector6d & stiffness)
{
  stiffness_ = stiffness;
}

void TorquePDCartesianTask::setDamping(const Eigen::Vector6d & damping)
{
  damping_ = damping;
}

void TorquePDCartesianTask::setPosTarget(const sva::PTransformd & xd)
{
  posTarget_ = xd;
}

void TorquePDCartesianTask::setVelTarget(const sva::MotionVecd & xd_dot)
{
  velTarget_ = xd_dot;
}

void TorquePDCartesianTask::setTorqueFeedforward(const Eigen::VectorXd & tau_ff)
{
  if(tau_ff.size() != nbActuatedJoints)
  {
    mc_rtc::log::error_and_throw<std::runtime_error>(
        "[TorquePDCartesianTask] Torque feedforward vector size should be {}, got {}", nbActuatedJoints, tau_ff.size());
  }
  torqueFeedforward_ = tau_ff;
}

void TorquePDCartesianTask::addToGUI(mc_rtc::gui::StateBuilder & gui)
{
  gui.addElement(
      {"Tasks", name_, "Gains"}, mc_rtc::gui::ArrayInput("Stiffness", stiffness_),
      mc_rtc::gui::ArrayInput("Damping", damping_),
      mc_rtc::gui::NumberInput(
          "Constant Stiffness & Critical Damping", [this]() { return stiffness_[0]; },
          [this](const double & g)
          {
            setStiffness(g);
            setDamping(2.0 * sqrt(g));
          }),
      mc_rtc::gui::NumberInput(
          "Constant Stiffness", [this]() { return stiffness_[0]; }, [this](const double & s) { setStiffness(s); }),
      mc_rtc::gui::NumberInput(
          "Constant Damping", [this]() { return damping_[0]; }, [this](const double & d) { setDamping(d); }));

  gui.addElement({"Tasks", name_, "Details"}, mc_rtc::gui::ArrayLabel("Position Error", posError_),
                 mc_rtc::gui::ArrayLabel("Velocity Error", velError_));

  gui.addElement({"Tasks", name_, "Position Target"},
                 mc_rtc::gui::Transform(
                     "pos_target", [this]() { return this->posTarget(); },
                     [this](const sva::PTransformd & X_0_target) { this->setPosTarget(X_0_target); }),
                 mc_rtc::gui::Transform("pos", [this]() { return this->frame_->position(); }));

  gui.addElement({"Tasks", name_, "Velocity Target"}, mc_rtc::gui::ArrayInput("vel_target", velTarget_));

  std::vector<std::string> active_gripper_joints;
  const auto & robot = robots_.robot(rIndex_);

  for(const auto & g : robot.grippers())
  {
    for(const auto & n : g.get().activeJoints()) { active_gripper_joints.push_back(n); }
  }
  auto isActiveGripperJoint = [&](const std::string & j)
  { return std::find(active_gripper_joints.begin(), active_gripper_joints.end(), j) != active_gripper_joints.end(); };

  int i = 0;
  for(const auto & j : robot.mb().joints())
  {
    if(j.dof() != 1 || j.isMimic() || isActiveGripperJoint(j.name())) { continue; }
    auto jIndex = robot.jointIndexByName(j.name());

    auto updateTorqueFeedforward = [this](int i, double v)
    {
      this->torqueFeedforward_[i] = v;
      setTorqueFeedforward(torqueFeedforward_);
    };

    gui.addElement({"Tasks", name_, "Torque Feedforward"},
                   mc_rtc::gui::NumberSlider(
                       j.name(), [this, i]() { return this->torqueFeedforward_[i]; },
                       [i, updateTorqueFeedforward](double v) { updateTorqueFeedforward(i, v); },
                       -robot.tl()[jIndex][0], robot.tu()[jIndex][0]));
    i++;
  }
  TorqueTask::addToGUI(gui);
}

void TorquePDCartesianTask::addToLogger(mc_rtc::Logger & logger)
{
  TorqueTask::addToLogger(logger);
  logger.removeLogEntry(name_ + "_torque");
  logger.addLogEntry(name_ + "_stiffness", [this]() { return stiffness_; });
  logger.addLogEntry(name_ + "_damping", [this]() { return damping_; });
  logger.addLogEntry(name_ + "_posTarget", [this]() { return posTarget_; });
  logger.addLogEntry(name_ + "_velTarget", [this]() { return velTarget_; });
  logger.addLogEntry(name_ + "_torqueFeedforward", [this]() { return torqueFeedforward_; });
  logger.addLogEntry(name_ + "_posError", [this]() { return posError_; });
  logger.addLogEntry(name_ + "_velError", [this]() { return velError_; });
  logger.addLogEntry(name_ + "_torque", [this]() { return torqueTarget_; });
}

} // namespace mc_tasks

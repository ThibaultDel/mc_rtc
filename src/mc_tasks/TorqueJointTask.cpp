#include <mc_tasks/TorqueJointTask.h>

#include <mc_rtc/gui/ArrayInput.h>
#include <mc_rtc/gui/ArrayLabel.h>
#include <mc_rtc/gui/Checkbox.h>
#include <mc_rtc/gui/NumberInput.h>
#include <mc_rtc/gui/NumberSlider.h>

namespace mc_tasks
{

TorqueJointTask::TorqueJointTask(const mc_solver::QPSolver & solver,
                                 unsigned int rIndex,
                                 double stiffness,
                                 double weight)
: TorqueTask(solver, rIndex, weight), robots_(solver.robots()), rIndex_(rIndex),
  nbActuatedJoints(
      (robots_.robot(rIndex_).mb().nrJoints() > 0 && robots_.robot(rIndex_).mb().joint(0).type() == rbd::Joint::Free)
          ? robots_.robot(rIndex_).mb().nrDof() - 6
          : robots_.robot(rIndex_).mb().nrDof()),
  stiffness_(Eigen::VectorXd::Zero(nbActuatedJoints)), damping_(Eigen::VectorXd::Zero(nbActuatedJoints)),
  integralTermEnabled_(false), integralGain_(Eigen::VectorXd::Zero(nbActuatedJoints)),
  maxIntegralTorque_(Eigen::VectorXd::Zero(nbActuatedJoints)), integralTorque_(Eigen::VectorXd::Zero(nbActuatedJoints)),
  posTarget_(Eigen::VectorXd::Zero(nbActuatedJoints)), velTarget_(Eigen::VectorXd::Zero(nbActuatedJoints)),
  torqueFeedforward_(Eigen::VectorXd::Zero(nbActuatedJoints)), posError_(Eigen::VectorXd::Zero(nbActuatedJoints)),
  velError_(Eigen::VectorXd::Zero(nbActuatedJoints)), integralError_(Eigen::VectorXd::Zero(nbActuatedJoints)),
  deriveVelocityTargetFromPosition_(false), prevPosTarget_(Eigen::VectorXd::Zero(nbActuatedJoints))
{
  if(backend_ == Backend::Tasks)
    mc_rtc::log::error_and_throw<std::runtime_error>(
        "[mc_tasks] Can't use TorqueJointTask with {} backend, please use TVM backend", backend_);
  name_ = std::string("torque_joint_") + solver.robots().robot(rIndex_).name();
  type_ = "torque_joint";

  setStiffness(stiffness);
  setDamping(2.0 * sqrt(stiffness)); // Critical damping by default
  setIntegralGain(0.01 * stiffness); // Small integral gain by default
  setMaxIntegralTorque(10); // Anti-windup limits for the integral term in Nm
  enableIntegralTerm(false); // Integral term disabled by default
  reset();
}

void TorqueJointTask::reset()
{
  posTarget_ = rbd::sParamToVector(robots_.robot(rIndex_).mb(), robots_.robot(rIndex_).q()).tail(nbActuatedJoints);
  prevPosTarget_ = posTarget_;
  velTarget_.setZero();
  torqueFeedforward_.setZero();
  integralError_.setZero();
  integralTorque_.setZero();
}

void TorqueJointTask::update(mc_solver::QPSolver & solver)
{
  Eigen::VectorXd torqueTarget = Eigen::VectorXd::Zero(nbActuatedJoints);
  auto & robot = solver.robots().robot(rIndex_);

  const auto & q_mbc = robot.q(); // MBC order
  const auto & q_dot_mbc = robot.alpha(); // MBC order
  const auto & refOrder = robot.refJointOrder();

  std::vector<double> q_map(refOrder.size());
  std::vector<double> q_dot_map(refOrder.size());

  for(size_t i = 0; i < refOrder.size(); ++i)
  {
    int mbcIndex = robot.jointIndexInMBC(i);
    if(mbcIndex >= 0)
    {
      q_map[i] = q_mbc[size_t(mbcIndex)][0];
      q_dot_map[i] = q_dot_mbc[size_t(mbcIndex)][0];
    }
    else
    {
      mc_rtc::log::warning(
          "[TorqueJointTask] Joint '{}' is in refJointOrder but not in mbc, skipping it in the control computation",
          refOrder[i]);
    }
  }

  Eigen::VectorXd q = Eigen::VectorXd::Map(q_map.data(), int(q_map.size()));
  Eigen::VectorXd q_dot = Eigen::VectorXd::Map(q_dot_map.data(), int(q_dot_map.size()));

  posError_ = posTarget_ - q;

  if(deriveVelocityTargetFromPosition_)
  {
    velTarget_ = (posTarget_ - prevPosTarget_) / solver.dt();
    prevPosTarget_ = posTarget_;
  }
  velError_ = velTarget_ - q_dot;

  if(integralTermEnabled_)
  {
    integralError_ += posError_ * solver.dt();
    Eigen::VectorXd tau_i = integralGain_.cwiseProduct(integralError_);
    integralTorque_ = // Anti-windup
        tau_i.cwiseMax(-maxIntegralTorque_).cwiseMin(maxIntegralTorque_);
    integralError_ = integralTorque_.cwiseQuotient(integralGain_);
    torqueTarget += integralTorque_;
  }

  torqueTarget += stiffness_.cwiseProduct(posError_);
  torqueTarget += damping_.cwiseProduct(velError_);
  torqueTarget += torqueFeedforward_;

  std::vector<std::vector<double>> torque_target = robot.mbc().jointTorque;

  int i = 0;
  for(const auto & joint_name : robot.refJointOrder())
  {
    torque_target[robot.jointIndexByName(joint_name)][0] = torqueTarget[i];
    i++;
  }

  TorqueTask::torqueTarget(torque_target);
  TorqueTask::update(solver);
}

void TorqueJointTask::setStiffness(double stiffness)
{
  setStiffness(Eigen::VectorXd::Constant(nbActuatedJoints, stiffness));
}

void TorqueJointTask::setDamping(double damping)
{
  setDamping(Eigen::VectorXd::Constant(nbActuatedJoints, damping));
}

void TorqueJointTask::setIntegralGain(double integralGain)
{
  setIntegralGain(Eigen::VectorXd::Constant(nbActuatedJoints, integralGain));
}

void TorqueJointTask::setMaxIntegralTorque(double maxIntegralTorque)
{
  setMaxIntegralTorque(Eigen::VectorXd::Constant(nbActuatedJoints, maxIntegralTorque));
}

void TorqueJointTask::setStiffness(const Eigen::VectorXd & stiffness)
{
  if(stiffness.size() != nbActuatedJoints)
  {
    mc_rtc::log::error_and_throw<std::runtime_error>("[TorqueJointTask] Stiffness vector size should be {}, got {}",
                                                     nbActuatedJoints, stiffness.size());
  }
  stiffness_ = stiffness;
}

void TorqueJointTask::setDamping(const Eigen::VectorXd & damping)
{
  if(damping.size() != nbActuatedJoints)
  {
    mc_rtc::log::error_and_throw<std::runtime_error>("[TorqueJointTask] Damping vector size should be {}, got {}",
                                                     nbActuatedJoints, damping.size());
  }
  damping_ = damping;
}

void TorqueJointTask::enableIntegralTerm(bool enable)
{
  integralError_.setZero();
  integralTorque_.setZero();
  integralTermEnabled_ = enable;
}

void TorqueJointTask::setIntegralGain(const Eigen::VectorXd & integralGain)
{
  if(integralGain.size() != nbActuatedJoints)
  {
    mc_rtc::log::error_and_throw<std::runtime_error>("[TorqueJointTask] Integral gain vector size should be {}, got {}",
                                                     nbActuatedJoints, integralGain.size());
  }
  integralGain_ = integralGain;
  enableIntegralTerm(true);
}

void TorqueJointTask::setMaxIntegralTorque(const Eigen::VectorXd & maxIntegralTorque)
{
  if(maxIntegralTorque.size() != nbActuatedJoints)
  {
    mc_rtc::log::error_and_throw<std::runtime_error>(
        "[TorqueJointTask] Max integral torque vector size should be {}, got {}", nbActuatedJoints,
        maxIntegralTorque.size());
  }
  maxIntegralTorque_ = maxIntegralTorque;
}

void TorqueJointTask::setPosTarget(const Eigen::VectorXd & qd)
{
  if(qd.size() != nbActuatedJoints)
  {
    mc_rtc::log::error_and_throw<std::runtime_error>(
        "[TorqueJointTask] Position target vector size should be {}, got {}", nbActuatedJoints, qd.size());
  }
  posTarget_ = qd;
}

void TorqueJointTask::setVelTarget(const Eigen::VectorXd & qd_dot)
{
  if(qd_dot.size() != nbActuatedJoints)
  {
    mc_rtc::log::error_and_throw<std::runtime_error>(
        "[TorqueJointTask] Velocity target vector size should be {}, got {}", nbActuatedJoints, qd_dot.size());
  }
  if(deriveVelocityTargetFromPosition_)
  {
    mc_rtc::log::warning("[TorqueJointTask] Trying to set velocity target while deriveVelocityTargetFromPosition is "
                         "enabled, ignoring the command");
    return;
  }
  velTarget_ = qd_dot;
}

void TorqueJointTask::setTorqueFeedforward(const Eigen::VectorXd & tau_ff)
{
  if(tau_ff.size() != nbActuatedJoints)
  {
    mc_rtc::log::error_and_throw<std::runtime_error>(
        "[TorqueJointTask] Torque feedforward vector size should be {}, got {}", nbActuatedJoints, tau_ff.size());
  }
  torqueFeedforward_ = tau_ff;
}

void TorqueJointTask::setDeriveVelocityTargetFromPosition(bool compute)
{
  deriveVelocityTargetFromPosition_ = compute;
  if(compute)
    prevPosTarget_ = posTarget_; // Reset previous position target to avoid large velocity target on the first update
}

void TorqueJointTask::addToGUI(mc_rtc::gui::StateBuilder & gui)
{
  std::vector<std::string> active_gripper_joints;
  const auto & robot = robots_.robot(rIndex_);
  std::vector<std::string> jointNames = robot.refJointOrder();

  for(const auto & g : robot.grippers())
  {
    for(const auto & n : g.get().activeJoints()) { active_gripper_joints.push_back(n); }
  }
  auto isActiveGripperJoint = [&](const std::string & j)
  { return std::find(active_gripper_joints.begin(), active_gripper_joints.end(), j) != active_gripper_joints.end(); };

  // jointNames.reserve(size_t(nbActuatedJoints));
  // for(const auto & joint : robot.mb().joints())
  // {
  //   if(joint.dof() == 1 && !joint.isMimic() && !isActiveGripperJoint(joint.name()))
  //   {
  //     jointNames.push_back(joint.name());
  //   }
  // }

  gui.addElement(
      {"Tasks", name_, "Gains"}, mc_rtc::gui::ArrayInput("Stiffness", jointNames, stiffness_),
      mc_rtc::gui::ArrayInput("Damping", jointNames, damping_),
      mc_rtc::gui::NumberInput(
          "Constant Stiffness", [this]() { return stiffness_[0]; }, [this](const double & s) { setStiffness(s); }),
      mc_rtc::gui::NumberInput(
          "Constant Damping", [this]() { return damping_[0]; }, [this](const double & d) { setDamping(d); }),
      mc_rtc::gui::Checkbox(
          "Enable Integral Term", [this]() { return integralTermEnabled(); },
          [this]() { enableIntegralTerm(!integralTermEnabled()); }),
      mc_rtc::gui::Checkbox(
          "Derive Velocity Target from Position Target", [this]() { return deriveVelocityTargetFromPosition(); },
          [this]() { setDeriveVelocityTargetFromPosition(!deriveVelocityTargetFromPosition()); }));

  gui.addElement({"Tasks", name_, "Gains", "Integral Term"},
                 mc_rtc::gui::ArrayInput("Integral Gain", jointNames, integralGain_),
                 mc_rtc::gui::ArrayInput("Max Integral Torque", jointNames, maxIntegralTorque_),
                 mc_rtc::gui::NumberInput(
                     "Constant Integral Gain", [this]() { return integralGain_[0]; },
                     [this](const double & g) { setIntegralGain(g); }),
                 mc_rtc::gui::NumberInput(
                     "Constant Max Integral Torque", [this]() { return maxIntegralTorque_[0]; },
                     [this](const double & t) { setMaxIntegralTorque(t); }),
                 mc_rtc::gui::ArrayLabel("Integral Torque", jointNames, integralTorque_));

  gui.addElement({"Tasks", name_, "Details"}, mc_rtc::gui::ArrayLabel("Position Error", jointNames, posError_),
                 mc_rtc::gui::ArrayLabel("Velocity Error", jointNames, velError_));

  int i = 0;
  for(const auto & j : robot.mb().joints())
  {
    if(j.dof() != 1 || j.isMimic() || isActiveGripperJoint(j.name())) { continue; }
    auto jIndex = robot.jointIndexByName(j.name());
    bool isContinuous = robot.ql()[jIndex][0] == -std::numeric_limits<double>::infinity();
    auto updatePosTarget = [this](int i, double v)
    {
      this->posTarget_[i] = v;
      setPosTarget(posTarget_);
    };

    auto updateVelTarget = [this](int i, double v)
    {
      this->velTarget_[i] = v;
      setVelTarget(velTarget_);
    };

    auto updateTorqueFeedforward = [this](int i, double v)
    {
      this->torqueFeedforward_[i] = v;
      setTorqueFeedforward(torqueFeedforward_);
    };

    if(isContinuous)
    {
      gui.addElement({"Tasks", name_, "Position Target"},
                     mc_rtc::gui::NumberInput(
                         j.name(), [this, i]() { return this->posTarget_[i]; },
                         [i, updatePosTarget](double v) { updatePosTarget(i, v); }));
    }
    else
    {
      gui.addElement({"Tasks", name_, "Position Target"},
                     mc_rtc::gui::NumberSlider(
                         j.name(), [this, i]() { return this->posTarget_[i]; }, [i, updatePosTarget](double v)
                         { updatePosTarget(i, v); }, robot.ql()[jIndex][0], robot.qu()[jIndex][0]));
    }

    gui.addElement({"Tasks", name_, "Velocity Target"},
                   mc_rtc::gui::NumberSlider(
                       j.name(), [this, i]() { return this->velTarget_[i]; }, [i, updateVelTarget](double v)
                       { updateVelTarget(i, v); }, robot.vl()[jIndex][0], robot.vu()[jIndex][0]));

    gui.addElement({"Tasks", name_, "Torque Feedforward"},
                   mc_rtc::gui::NumberSlider(
                       j.name(), [this, i]() { return this->torqueFeedforward_[i]; },
                       [i, updateTorqueFeedforward](double v) { updateTorqueFeedforward(i, v); },
                       -robot.tl()[jIndex][0], robot.tu()[jIndex][0]));
    i++;
  }
  TorqueTask::addToGUI(gui);
}

void TorqueJointTask::addToLogger(mc_rtc::Logger & logger)
{
  logger.addLogEntry(name_ + "_stiffness", [this]() { return stiffness_; });
  logger.addLogEntry(name_ + "_damping", [this]() { return damping_; });
  logger.addLogEntry(name_ + "_posTarget", [this]() { return posTarget_; });
  logger.addLogEntry(name_ + "_velTarget", [this]() { return velTarget_; });
  logger.addLogEntry(name_ + "_torqueFeedforward", [this]() { return torqueFeedforward_; });
  logger.addLogEntry(name_ + "_posError", [this]() { return posError_; });
  logger.addLogEntry(name_ + "_velError", [this]() { return velError_; });
  logger.addLogEntry(name_ + "_deriveVelocityTargetFromPosition",
                     [this]() { return deriveVelocityTargetFromPosition_; });
  logger.addLogEntry(name_ + "_integralTermEnabled", [this]() { return integralTermEnabled_; });
  logger.addLogEntry(name_ + "_integralGain", [this]() { return integralGain_; });
  logger.addLogEntry(name_ + "_maxIntegralTorque", [this]() { return maxIntegralTorque_; });
  logger.addLogEntry(name_ + "_integralError", [this]() { return integralError_; });
  logger.addLogEntry(name_ + "_integralTorque", [this]() { return integralTorque_; });
  TorqueTask::addToLogger(logger);
}

} // namespace mc_tasks

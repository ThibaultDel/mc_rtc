#include <mc_tasks/MetaTask.h>
#include <mc_tasks/TorqueCartesianTask.h>
#include <mc_tvm/WrenchFunction.h>

#include <mc_rtc/gui/ArrayInput.h>
#include <mc_rtc/gui/ArrayLabel.h>
#include <mc_rtc/gui/Checkbox.h>
#include <mc_rtc/gui/NumberInput.h>
#include <mc_rtc/gui/NumberSlider.h>
#include <mc_rtc/gui/Transform.h>
#include <mc_tvm/Robot.h>
#include <mc_tvm/RobotFrame.h>
#include <SpaceVecAlg/EigenTypedef.h>
#include <SpaceVecAlg/SpaceVecAlg>
#include "mc_rtc/logging.h"
#include "mc_tasks/WrenchTask.h"
#include <Eigen/src/Core/Matrix.h>
#include <vector>

namespace mc_tasks
{

TorqueCartesianTask::TorqueCartesianTask(const mc_solver::QPSolver & solver,
                                         const mc_rbdyn::RobotFrame & frame,
                                         double stiffness,
                                         double weight)
: WrenchTask(frame, weight), posTarget_(sva::PTransformd::Identity()), velTarget_(sva::MotionVecd::Zero()),
  robots_(solver.robots()), rIndex_(frame.robot().robotIndex()),
  nbActuatedJoints(
      (robots_.robot(rIndex_).mb().nrJoints() > 0 && robots_.robot(rIndex_).mb().joint(0).type() == rbd::Joint::Free)
          ? robots_.robot(rIndex_).mb().nrDof() - 6
          : robots_.robot(rIndex_).mb().nrDof()),
  stiffness_(Eigen::Vector6d::Zero()), damping_(Eigen::Vector6d::Zero()), integralGain_(Eigen::Vector6d::Zero()),
  maxIntegralWrench_(sva::ForceVecd::Zero()), integralWrench_(sva::ForceVecd::Zero()),
  torqueFeedforward_(Eigen::VectorXd::Zero(nbActuatedJoints)),
  torqueExtForcesCompensation_(Eigen::VectorXd::Zero(nbActuatedJoints)),
  wrenchGravityCompensation_(sva::ForceVecd::Zero()), posError_(sva::MotionVecd::Zero()),
  velError_(sva::MotionVecd::Zero()), integralError_(sva::MotionVecd::Zero()), wrenchTarget_(sva::ForceVecd::Zero()),
  prevPosTarget_(sva::PTransformd::Identity()), omega_n_(Eigen::Vector6d::Constant(5.0)),
  zeta_(Eigen::Vector6d::Constant(1.0))
{
  if(backend_ == Backend::Tasks)
    mc_rtc::log::error_and_throw<std::runtime_error>(
        "[mc_tasks] Can't use TorqueCartesianTask with {} backend, please use TVM backend", backend_);

  name_ = "torque_cartesian_" + frame.name() + "_" + solver.robots().robot(frame.robot().robotIndex()).name();
  type_ = "torque_cartesian";

  setStiffness(stiffness);
  setDamping(2.0 * sqrt(stiffness)); // Critical damping by default
  setIntegralGain(0.01 * stiffness); // Small integral gain by default
  setMaxIntegralWrench(5, 10); // Anti-windup limits for the integral term
  enableIntegralTerm(false); // Integral term disabled by default
  reset();
}

TorqueCartesianTask::TorqueCartesianTask(const mc_solver::QPSolver & solver,
                                         const std::string & bodyName,
                                         unsigned int rIndex,
                                         double stiffness,
                                         double weight)
: TorqueCartesianTask(solver, solver.robots().robot(rIndex).frame(bodyName), stiffness, weight)
{
}

void TorqueCartesianTask::reset()
{
  posTarget_ = frame().position();
  prevPosTarget_ = posTarget_;
  velTarget_.Zero();
  torqueFeedforward_.setZero();
  integralError_.Zero();
  integralWrench_.Zero();
  // WrenchTask::reset();
}

void TorqueCartesianTask::update(mc_solver::QPSolver & solver)
{
  wrenchTarget_ = sva::ForceVecd::Zero();

  // tau_ff + tau_ext_forces
  Eigen::VectorXd torqueComponents = Eigen::VectorXd::Zero(nbActuatedJoints);

  // Compute Force components (PD + integral if enabled)
  posError_ = sva::transformError(frame().position(), posTarget_);

  if(deriveVelocityTargetFromPosition_)
  {
    velTarget_ = sva::transformError(prevPosTarget_, posTarget_) / solver.dt();
    prevPosTarget_ = posTarget_;
  }

  velError_ = velTarget_ - frame().velocity();

  if(integralTermEnabled_)
  {
    integralError_ += posError_ * solver.dt();
    Eigen::Vector6d force_i = integralGain_.cwiseProduct(integralError_.vector());
    integralWrench_ = // Anti-windup
        sva::ForceVecd(force_i.head<3>().cwiseMax(-maxIntegralWrench_.force()).cwiseMin(maxIntegralWrench_.force()),
                       force_i.tail<3>().cwiseMax(-maxIntegralWrench_.couple()).cwiseMin(maxIntegralWrench_.couple()));

    Eigen::Vector6d integralErrorVec = integralWrench_.vector().cwiseQuotient(integralGain_);
    integralError_ = sva::MotionVecd(integralErrorVec.head<3>(), integralErrorVec.tail<3>());
    wrenchTarget_ += integralWrench_;
  }

  // Now all quantities are in local frame, consistent with cartesianInertia
  // wrenchTarget_ += sva::ForceVecd(cartesianInertia() * stiffness_.cwiseProduct(posError_.vector()));
  // wrenchTarget_ += sva::ForceVecd(cartesianInertia() * damping_.cwiseProduct(velError_.vector()));

  Eigen::Vector6d pdWrenchVec =
      cartesianInertia() * (stiffness_.cwiseProduct(posError_.vector()) + damping_.cwiseProduct(velError_.vector()));
  wrenchTarget_ += sva::ForceVecd(pdWrenchVec);

  // Compute torque components (feedforward + compensation)
  auto & tvm_robot = robots_.robot(rIndex_).tvmRobot();
  torqueComponents += torqueFeedforward_;
  torqueExtForcesCompensation_ = tvm_robot.tauExternal();
  if(compensateExternalForces_) { torqueComponents += torqueExtForcesCompensation_; }

  computeForceGravityCompensation();
  if(compensateGravity_) { wrenchTarget_ += wrenchGravityCompensation_; }

  wrenchTarget_ += sva::ForceVecd(WrenchTask::dynamicJacobianTranspose() * torqueComponents);
  WrenchTask::targetWrench(wrenchTarget_);
}

void TorqueCartesianTask::computeForceGravityCompensation()
{
  auto & robot = robots_.robot(rIndex_);
  Eigen::VectorXd Cg = robot.tvmRobot().C(); // C*qdot + g
  mc_tvm::RobotFrame & tvm_frame = frame().tvm_frame();
  // rbd::Jacobian jacobian(robot.mb(), frame().body());
  rbd::Jacobian jacobian(tvm_frame.rbdJacobian());
  sva::MotionVecd jdotAlpha = jacobian.normalAcceleration(robot.mb(), robot.mbc(),
                                                          tvm_frame.position(), // X_b_p: same point as Jacobian
                                                          tvm_frame.velocity() // V_b_p: velocity of that frame point
  );
  // sva::MotionVecd jdotAlpha = jacobian.normalAcceleration(robot.mb(), robot.mbc());
  Eigen::Vector6d forceGravityVec =
      WrenchTask::dynamicJacobianTranspose() * Cg - WrenchTask::cartesianInertia() * jdotAlpha.vector();
  wrenchGravityCompensation_ = sva::ForceVecd(forceGravityVec.head<3>(), forceGravityVec.tail<3>());
}

void TorqueCartesianTask::setTorqueFeedforward(const Eigen::VectorXd & tau_ff)
{
  if(tau_ff.size() != nbActuatedJoints)
  {
    mc_rtc::log::error_and_throw<std::runtime_error>(
        "[TorqueCartesianTask] Torque feedforward vector size should be {}, got {}", nbActuatedJoints, tau_ff.size());
  }
  torqueFeedforward_ = tau_ff;
}

void TorqueCartesianTask::setDeriveVelocityTargetFromPosition(bool compute)
{
  deriveVelocityTargetFromPosition_ = compute;
  if(compute)
    prevPosTarget_ = posTarget_; // Reset previous position target to avoid large velocity target on the first update
}

void TorqueCartesianTask::enableIntegralTerm(bool enable)
{
  integralTermEnabled_ = enable;
  if(!enable)
  {
    integralError_.Zero();
    integralWrench_.Zero();
  }
}

void TorqueCartesianTask::setUseAutoDamping(bool use)
{
  useAutoDamping_ = use;
  if(useAutoDamping_)
  {
    stiffness_ = omega_n_.array().square();
    damping_ = (2 * zeta_.array() * omega_n_.array()).matrix();
  }
}
void TorqueCartesianTask::setAutoDampingNaturalFrequency(const Eigen::Vector6d & omega_n)
{
  omega_n_ = omega_n;
  if(useAutoDamping_)
  {
    stiffness_ = omega_n_.array().square();
    damping_ = (2 * zeta_.array() * omega_n_.array()).matrix();
  }
}
void TorqueCartesianTask::setAutoDampingDampingRatio(const Eigen::Vector6d & zeta)
{
  zeta_ = zeta;
  if(useAutoDamping_) { damping_ = (2 * zeta_.array() * omega_n_.array()).matrix(); }
}

void TorqueCartesianTask::addToGUI(mc_rtc::gui::StateBuilder & gui)
{

  std::vector<std::string> active_gripper_joints;
  std::vector<std::string> jointNames;
  std::vector<std::string> gainNames = {"Rx", "Ry", "Rz", "Tx", "Ty", "Tz"};
  const auto & robot = robots_.robot(rIndex_);

  for(const auto & g : robot.grippers())
  {
    for(const auto & n : g.get().activeJoints()) { active_gripper_joints.push_back(n); }
  }
  auto isActiveGripperJoint = [&](const std::string & j)
  { return std::find(active_gripper_joints.begin(), active_gripper_joints.end(), j) != active_gripper_joints.end(); };

  jointNames.reserve(size_t(nbActuatedJoints));
  for(const auto & joint : robot.mb().joints())
  {
    if(joint.dof() == 1 && !joint.isMimic() && !isActiveGripperJoint(joint.name()))
    {
      jointNames.push_back(joint.name());
    }
  }

  gui.addElement(
      {"Tasks", name_, "Gains"}, mc_rtc::gui::ArrayInput("Stiffness", gainNames, stiffness_),
      mc_rtc::gui::ArrayInput("Damping", gainNames, damping_),
      mc_rtc::gui::NumberInput(
          "Constant Stiffness", [this]() { return stiffness_[0]; }, [this](const double & s) { setStiffness(s); }),
      mc_rtc::gui::NumberInput(
          "Constant Damping", [this]() { return damping_[0]; }, [this](const double & d) { setDamping(d); }),
      mc_rtc::gui::Checkbox(
          "Enable Integral Term", [this]() { return integralTermEnabled(); },
          [this]() { enableIntegralTerm(!integralTermEnabled()); }),
      mc_rtc::gui::Checkbox(
          "Derive Velocity Target from Position Target", [this]() { return deriveVelocityTargetFromPosition(); },
          [this]() { setDeriveVelocityTargetFromPosition(!deriveVelocityTargetFromPosition()); }),
      mc_rtc::gui::NumberInput(
          "Natural Frequency (for auto damping) Rotation", [this]() { return omega_n_[0]; },
          [this](const double & f)
          {
            Eigen::Vector6d omega_n = autoDampingNaturalFrequency();
            omega_n.head<3>().setConstant(f);
            setAutoDampingNaturalFrequency(omega_n);
          }),
      mc_rtc::gui::NumberInput(
          "Damping Ratio (for auto damping) Rotation", [this]() { return zeta_[0]; },
          [this](const double & r)
          {
            Eigen::Vector6d zeta = autoDampingDampingRatio();
            zeta.head<3>().setConstant(r);
            setAutoDampingDampingRatio(zeta);
          }),
      mc_rtc::gui::NumberInput(
          "Natural Frequency (for auto damping) Translation", [this]() { return omega_n_[3]; },
          [this](const double & f)
          {
            Eigen::Vector6d omega_n = autoDampingNaturalFrequency();
            omega_n.tail<3>().setConstant(f);
            setAutoDampingNaturalFrequency(omega_n);
          }),
      mc_rtc::gui::NumberInput(
          "Damping Ratio (for auto damping) Translation", [this]() { return zeta_[3]; },
          [this](const double & r)
          {
            Eigen::Vector6d zeta = autoDampingDampingRatio();
            zeta.tail<3>().setConstant(r);
            setAutoDampingDampingRatio(zeta);
          }),
      mc_rtc::gui::Checkbox(
          "Use Automatic Damping (based on natural frequency and damping ratio)", [this]() { return useAutoDamping_; },
          [this]() { setUseAutoDamping(!useAutoDamping()); }));

  gui.addElement({"Tasks", name_, "Gains", "Integral Term"},
                 mc_rtc::gui::ArrayInput("Integral Gain", gainNames, integralGain_),
                 mc_rtc::gui::ArrayInput("Max Integral Wrench", maxIntegralWrench_),
                 mc_rtc::gui::NumberInput(
                     "Constant Integral Gain", [this]() { return integralGain_[0]; },
                     [this](const double & g) { setIntegralGain(g); }),
                 mc_rtc::gui::NumberInput(
                     "Constant Max Integral Force", [this]() { return maxIntegralWrench_.force()[0]; },
                     [this](const double & w) { setMaxIntegralForce(w); }),
                 mc_rtc::gui::NumberInput(
                     "Constant Max Integral Torque", [this]() { return maxIntegralWrench_.couple()[0]; },
                     [this](const double & w) { setMaxIntegralTorque(w); }));

  gui.addElement(
      {"Tasks", name_, "Additional Forces"},
      mc_rtc::gui::Checkbox(
          "Compensate External Forces", [this]() { return isCompensatingExternalForces(); },
          [this]() { setCompensateExternalForces(!isCompensatingExternalForces()); }),
      mc_rtc::gui::Checkbox(
          "Compensate Gravity (+ Coriolis)", [this]() { return isCompensatingGravity(); },
          [this]() { setCompensateGravity(!isCompensatingGravity()); }),
      mc_rtc::gui::ArrayLabel("Torque External Forces", jointNames, [this]() { return this->torqueExternalForces(); }),
      mc_rtc::gui::ArrayLabel("Gravity Wrench", [this]() { return this->gravityWrench(); }));

  gui.addElement({"Tasks", name_, "Details"}, mc_rtc::gui::ArrayLabel("Position Error", posError_),
                 mc_rtc::gui::ArrayLabel("Velocity Error", velError_));

  gui.addElement({"Tasks", name_, "Position Target"},
                 mc_rtc::gui::Transform(
                     "pos_target", [this]() { return this->posTarget(); },
                     [this](const sva::PTransformd & X_0_target) { this->setPosTarget(X_0_target); }),
                 mc_rtc::gui::Transform("pos", [this]() { return this->frame().position(); }));

  gui.addElement({"Tasks", name_, "Velocity Target"}, mc_rtc::gui::ArrayInput("vel_target", velTarget_));

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
  WrenchTask::addToGUI(gui);
}

void TorqueCartesianTask::addToLogger(mc_rtc::Logger & logger)
{
  WrenchTask::addToLogger(logger);
  logger.addLogEntry(name_ + "_stiffness", [this]() { return stiffness_; });
  logger.addLogEntry(name_ + "_damping", [this]() { return damping_; });
  logger.addLogEntry(name_ + "_posTarget", [this]() { return posTarget_; });
  logger.addLogEntry(name_ + "_velTarget", [this]() { return velTarget_; });
  logger.addLogEntry(name_ + "_torqueFeedforward", [this]() { return torqueFeedforward_; });
  logger.addLogEntry(name_ + "_posError", [this]() { return posError_; });
  logger.addLogEntry(name_ + "_velError", [this]() { return velError_; });
  logger.addLogEntry(name_ + "_deriveVelocityTargetFromPosition",
                     [this]() { return deriveVelocityTargetFromPosition_; });
  logger.addLogEntry(name_ + "_compensateExternalForces", [this]() { return compensateExternalForces_; });
  logger.addLogEntry(name_ + "_gravityWrenchCompensation", [this]() { return wrenchGravityCompensation_; });
  logger.addLogEntry(name_ + "_compensateGravity", [this]() { return compensateGravity_; });
  logger.addLogEntry(name_ + "_integralTermEnabled", [this]() { return integralTermEnabled_; });
  logger.addLogEntry(name_ + "_integralGain", [this]() { return integralGain_; });
  logger.addLogEntry(name_ + "_maxIntegralWrench", [this]() { return maxIntegralWrench_; });
  logger.addLogEntry(name_ + "_integralError", [this]() { return integralError_; });
  logger.addLogEntry(name_ + "_integralWrench", [this]() { return integralWrench_; });
}

} // namespace mc_tasks

/*
 * Copyright 2015-2026 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_tasks/MetaTask.h>

#include <mc_rtc/void_ptr.h>

namespace mc_tasks
{

/*! \brief Torque tracking task based on the robot dynamic model.
 *
 * The TorqueTask computes the joint accelerations that minimize the
 * error between the desired torque vector and the torque produced by the
 * dynamic model. This task is therefore intended to be used in torque-control
 * mode only.
 *
 * The torque target is expressed directly in joint space and represents the
 * torques that should be realized at the actuated joints.
 *
 * As an additional feature, the task can optionally compensate for external
 * torques acting on the system if such estimates are available.
 *
 * By default, the desired torque target is zero and no external torque
 * compensation is applied.
 *
 */
struct MC_TASKS_DLLAPI TorqueTask : public MetaTask
{
public:
  TorqueTask(const mc_solver::QPSolver & solver,
             unsigned int rIndex,
             double weight = 10,
             bool compensateExternalForces = false,
             bool compensateGravity = false);

  void reset() override;

  /*! \brief Load parameters from a Configuration object */
  void load(mc_solver::QPSolver & solver, const mc_rtc::Configuration & config) override;

  /*! \brief Set the task dimensional weight
   *
   * For simple cases (using 0/1 as weights) prefer \ref selectActiveJoints or \ref selectUnactiveJoints which are
   * simpler to use
   */
  void dimWeight(const Eigen::VectorXd & dimW) override;

  Eigen::VectorXd dimWeight() const override;

  /*! \brief Select active joints for this task
   *
   * Manipulate \ref dimWeight() to achieve its effect
   */
  void selectActiveJoints(mc_solver::QPSolver & solver,
                          const std::vector<std::string> & activeJointsName,
                          const std::map<std::string, std::vector<std::array<int, 2>>> & activeDofs = {}) override;

  /*! \brief Select inactive joints for this task
   *
   * Manipulate \ref dimWeight() to achieve its effect
   */
  void selectUnactiveJoints(mc_solver::QPSolver & solver,
                            const std::vector<std::string> & unactiveJointsName,
                            const std::map<std::string, std::vector<std::array<int, 2>>> & unactiveDofs = {}) override;

  /*! \brief Reset the joint selector effect
   *
   * Reset dimWeight to ones
   */
  void resetJointsSelector(mc_solver::QPSolver & solver) override;

  Eigen::VectorXd eval() const override;

  Eigen::VectorXd speed() const override;

  Eigen::VectorXd torqueTargetVector() const { return torque_vector_; }

  /** Get current torque objective */
  std::vector<std::vector<double>> torqueTarget() const;

  /** Set joint weights for the torque task */
  void jointWeights(const std::map<std::string, double> & jws);

  /** Set specific joint targets
   *
   * \param tau Map of joint's name to torque's configuration
   *
   */
  void target(const std::map<std::string, std::vector<double>> & tau);

  /** Set task's weight */
  void weight(double w);

  /** Get task's weight */
  double weight() const;

  /** True if the task is in the solver */
  bool inSolver() const;

  /** Set if the task is compensating external forces */
  void setCompensateExternalForces(bool compensate);

  /** True if the task is compensating external forces */
  bool isCompensatingExternalForces();

  /** Set if the task is compensating gravity */
  void setCompensateGravity(bool compensate);

  /** True if the task is compensating gravity */
  bool isCompensatingGravity();

  Eigen::VectorXd torqueExternalForces() const;
  Eigen::VectorXd torqueGravity() const;

protected:
  void addToSolver(mc_solver::QPSolver & solver) override;

  void removeFromSolver(mc_solver::QPSolver & solver) override;

  void update(mc_solver::QPSolver &) override;

  void addToGUI(mc_rtc::gui::StateBuilder &) override;

  void addToLogger(mc_rtc::Logger & logger) override;

  /** Change torque objective without precising joint names */
  void torqueTarget(const std::vector<std::vector<double>> & tau);

private:
  /** True if added to solver */
  bool inSolver_ = false;
  /** Robot handled by the task */
  const mc_rbdyn::Robots & robots_;
  unsigned int rIndex_;
  /** Holds the constraint implementation
   *
   * In Tasks backend:
   * - None
   *
   * In TVM backend:
   * - details::TVMTorqueTask
   */
  mc_rtc::void_ptr pt_;
  /** Solver timestep */
  double dt_;
  /** True if the task is compensating external forces */
  bool compensateExternalForces_ = false;
  /** True if the task is compensating gravity */
  bool compensateGravity_ = false;

  /** Store the target torque */
  std::vector<std::vector<double>> torque_;
  /** Store the target torque vector */
  Eigen::VectorXd torque_vector_;
  /** Store mimic information */
  std::unordered_map<std::string, std::vector<int>> mimics_;
  /** Store the previous eval vector */
  Eigen::VectorXd eval_;
  /** Store the task speed */
  Eigen::VectorXd speed_;
};

using TorqueTaskPtr = std::shared_ptr<TorqueTask>;

} // namespace mc_tasks

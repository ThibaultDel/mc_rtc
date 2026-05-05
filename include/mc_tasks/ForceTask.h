/*
 * Copyright 2015-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <SpaceVecAlg/SpaceVecAlg>
#include <Eigen/src/Core/Matrix.h>

#pragma once

#include <mc_tasks/MetaTask.h>

#include <mc_rtc/void_ptr.h>

#include <Tasks/QPTasks.h>

namespace mc_tasks
{

/*! \brief Controls an end-effector force
 *
 * This task is a thin wrapper around the appropriate tasks in Tasks.
 * The task objective is given in the world frame.
 */
struct MC_TASKS_DLLAPI ForceTask : public MetaTask
{
public:
  /*! \brief Constructor
   *
   * \param solver QP solver
   *
   * \param robots Robots controlled by this task
   *
   * \param robotIndex Index of the robot controlled by this task
   *
   * \param weight Task weight
   *
   * \param compensateExternalForces If true, the task will try to compensate for external forces acting on the robot
   *
   */
  ForceTask(const mc_solver::QPSolver & solver,
            const mc_rbdyn::Robots & robots,
            unsigned int robotIndex,
            Eigen::MatrixXd jTransposePseudoInverse,
            double weight = 1000.0,
            bool compensateExternalForces = false);

  /*! \brief Reset the task
   *
   * Set the task objective to the current end-effector position
   */
  void reset() override;

  /*! \brief Increment the target position
   *
   * \param dtr Change in target position
   *
   */
  void addForceTarget(const Eigen::VectorXd & dtr);

  /*! \brief Change the target position
   *
   * \param tf New target position
   *
   */
  void setForceTarget(const Eigen::VectorXd & tf);

  /*! \brief Returns the current target positions
   *
   * \returns Current target position
   *
   */
  Eigen::VectorXd getForceTarget() const;

  void dimWeight(const Eigen::VectorXd & dimW) override;

  Eigen::VectorXd dimWeight() const override;

  void selectActiveJoints(mc_solver::QPSolver & solver,
                          const std::vector<std::string> & activeJointsName,
                          const std::map<std::string, std::vector<std::array<int, 2>>> & activeDofs = {}) override;

  void selectUnactiveJoints(mc_solver::QPSolver & solver,
                            const std::vector<std::string> & unactiveJointsName,
                            const std::map<std::string, std::vector<std::array<int, 2>>> & unactiveDofs = {}) override;

  void resetJointsSelector(mc_solver::QPSolver & solver) override;

  // Eigen::VectorXd getCurrentForce() const;

  /** Set task's weight */
  void weight(double w);

  /** Get task's weight */
  double weight() const;

  /** True if the task is in the solver */
  bool inSolver() const;

  /** Set if the task is compensating external forces */
  void compensateExternalForces(bool compensate);

  /** True if the task is compensating external forces */
  bool isCompensatingExternalForces() const;

  void setJacobianTPseudoInverse(const Eigen::MatrixXd & jTransposePseudoInverse);

  Eigen::MatrixXd getJacobianTPseudoInverse() const;

  Eigen::VectorXd eval() const override;

  Eigen::VectorXd speed() const override;

  void load(mc_solver::QPSolver & solver, const mc_rtc::Configuration & config) override;

  using MetaTask::name;

  void name(const std::string & name) override;

private:
  Eigen::VectorXd target_;
  // Eigen::Vector6d targetVector_;
  // Eigen::VectorXd curForce_;
  // Eigen::Vector6d curForceVector_;

  /** True if added to solver */
  bool inSolver_ = false;
  /** Robot handled by the task */
  const mc_rbdyn::Robots & robots_;
  Eigen::MatrixXd jTransposePseudoInverse_;

  unsigned int rIndex_;
  /** Holds the constraint implementation
   *
   * In Tasks backend:
   * - None
   *
   * In TVM backend:
   * - details::TVMForceTask
   */
  mc_rtc::void_ptr pt_;
  /** Solver timestep */
  double dt_;
  /** True if the task is compensating external forces */
  bool compensateExternalForces_ = false;
  /** Store the previous eval vector */
  Eigen::VectorXd eval_;
  /** Store the task speed */
  Eigen::VectorXd speed_;

protected:
  void removeFromSolver(mc_solver::QPSolver & solver) override;

  void addToSolver(mc_solver::QPSolver & solver) override;

  void update(mc_solver::QPSolver &) override;

  void addToLogger(mc_rtc::Logger & logger) override;

  void removeFromLogger(mc_rtc::Logger & logger) override;

  void addToGUI(mc_rtc::gui::StateBuilder & gui) override;
};

} // namespace mc_tasks

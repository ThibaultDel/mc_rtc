/*
 * Copyright 2015-2022 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include <mc_tvm/ImpulseFunction.h>

#include <mc_tvm/Robot.h>
#include <mc_tvm/RobotFrame.h>

#include <numeric>

namespace mc_tvm
{

ImpulseFunction::ImpulseFunction(const mc_rbdyn::Robot & robot, const mc_rbdyn::RobotFrame & frame, double delta_t, double c_res, double limit_multiplier, int axis)
: tvm::function::abstract::LinearFunction(robot.mb().nrDof()), robot_(robot), frame_(frame), delta_t_(delta_t), c_res_(c_res), limit_multiplier_(limit_multiplier), axis_(axis)
  , jac_(frame.tvm_frame().rbdJacobian())
{
  assert(frame_->robot().robotIndex() == robot_.robotIndex() && "ImpulseFunction frame must belong to the robot");
  assert(axis_ >= 0 && axis_ < 3 && "Axis must be in [0, 2] as there are only three translational axis in the frame (1=x, 2=y, 3=z)");
  // registerUpdates(Update::B, &ImpulseFunction::updateb);
  registerUpdates(Update::Jacobian, &ImpulseFunction::updateJacobian);
  // addOutputDependency<ImpulseFunction>(Output::B, Update::B);
  addOutputDependency<ImpulseFunction>(Output::Jacobian, Update::Jacobian);
  auto & tvm_robot = robot.tvmRobot();
  addInputDependency<ImpulseFunction>(Update::Jacobian, tvm_robot, Robot::Output::H);
  // addInputDependency<ImpulseFunction>(Update::B, tvm_robot, Robot::Output::C); // TODO check if needed
  // addVariable(tvm::dot(tvm_robot.q(), 2), true);
  // addVariable(tvm::dot(tvm_robot.q(), 1), true);
  q_ddot_var_ = tvm::dot(tvm_robot.q(), 2);
  q_dot_var_  = tvm::dot(tvm_robot.q(), 1);
  addVariable(q_ddot_var_, true);
  addVariable(q_dot_var_, true);
  q_ddot_vars_.add(q_ddot_var_);
  q_dot_vars_.add(q_dot_var_);
  // b_ = - robot_.tvmRobot().limits().tl * limit_multiplier_ / (c_res_+1);
  b_ = Eigen::VectorXd::Zero(robot_.mb().nrDof());

  if(!q_ddot_var_ || !q_dot_var_)
  {
    mc_rtc::log::error_and_throw("[ImpulseFunction] created variable is null");
  }
  mc_rtc::log::info("[ImpulseFunction] q_ddot_var_ = {} use_count = {}", (void*)q_ddot_var_.get(), q_ddot_var_.use_count());
  mc_rtc::log::info("[ImpulseFunction] q_dot_var_  = {} use_count = {}", (void*)q_dot_var_.get(),  q_dot_var_.use_count());

  velocity_.setZero();
}


void ImpulseFunction::updateb()
{
  // b_ = - robot_.tvmRobot().limits().tl * limit_multiplier_ / (c_res_+1);
  // b_ = robot_.tvmRobot().limits().tl;

  // mc_rtc::log::info("the evaluation of the function:");
  // mc_rtc::log::info(value());

  b_ = Eigen::VectorXd::Zero(robot_.mb().nrDof());

  // mc_rtc::log::info(" actual valuation:");
  // mc_rtc::log::info(value());

}

void ImpulseFunction::updateJacobian()
{
  // // mc_rtc::log::info("the evaluation of the function:");
  // // mc_rtc::log::info(b());
  //
  const auto & robot = robot_.tvmRobot();
  auto M = robot.H();

  // fixed asserts (no chained ==)
  assert(M.rows() == robot_.mb().nrDof());
  assert(M.cols() == robot_.mb().nrDof());
  assert(std::abs(M.determinant()) > 1e-5 && "Mass matrix is singular");

  // Get jacobian of the robot frame
  const rbd::MultiBody robot_mb = robot_.mb();
  const rbd::MultiBodyConfig mbc = robot_.mbc();
  const auto & world_frame_jacobian = jac_.jacobian(robot_mb, mbc);
  Eigen::MatrixXd full_world_frame_jacobian(6, robot_.mb().nrDof());
  jac_.fullJacobian(robot_mb, world_frame_jacobian, full_world_frame_jacobian);

  assert(full_world_frame_jacobian.cols() == robot_.mb().nrDof());

  // Get normal of hitting plane
  auto normal = robot_.frame(frame_->name()).position().rotation().col(axis_);

  // // Generate necessary matrices
  Eigen::MatrixXd j_m_uninverted = /*full_world_frame_jacobian.transpose() * */(full_world_frame_jacobian * M.inverse() * full_world_frame_jacobian.transpose());
  assert(std::abs(j_m_uninverted.determinant()) > 1e-5 && "Mass matrix is singular");
  Eigen::MatrixXd j_m = full_world_frame_jacobian.transpose() * j_m_uninverted.inverse();
  Eigen::MatrixXd P_n_sub = normal * normal.transpose();
  Eigen::Matrix<double, 6, 6> P_n = Eigen::Matrix<double, 6, 6>::Zero();
  P_n.block<3,3>(0,0) = P_n_sub;

  // full matrix: nrDof x nrDof
  J_ddq = j_m * P_n * full_world_frame_jacobian;              // multiplies ddq variable
  J_dq  = j_m * P_n * full_world_frame_jacobian / delta_t_;  // multiplies dq variable

  auto evaluation1 = J_ddq * robot.alphaD()->value();
  auto evaluation2 = /*J_ddq * robot.alphaD()->value() + */J_dq * robot.alpha()->value() /*- b_*/;
  auto evaluation3 = /*J_ddq * robot.alphaD()->value() + *//*J_dq * robot.alpha()->value()*/ - b_;
  auto evaluation = evaluation1 + evaluation2 + evaluation3;
  // auto evaluation = J_ddq * robot.alphaD()->value()/* + J_dq * robot.alpha()->value()*//* - b_*/;


  // mc_rtc::log::info("evaluation1:");
  // mc_rtc::log::info(evaluation1.transpose());
  // mc_rtc::log::info("evaluation2:");
  // mc_rtc::log::info(evaluation2.transpose());
  // mc_rtc::log::info("evaluation3:");
  // mc_rtc::log::info(evaluation3.transpose());
  // mc_rtc::log::info("evaluation:");
  // mc_rtc::log::info(evaluation);


  // J_ddq = Eigen::MatrixXd::Identity(robot_mb.nrDof(), robot_mb.nrDof());              // multiplies ddq variable
  // J_dq  = Eigen::MatrixXd::Identity(robot_mb.nrDof(), robot_mb.nrDof());  // multiplies dq variable


  // Sanity checks: variables non-null and sizes match
  if(!q_ddot_var_ || !q_dot_var_)
  {
    mc_rtc::log::error_and_throw("[ImpulseFunction] variable ptr became null before splitJacobian");
  }

  auto sumSizes = [](const tvm::VariableVector & vars)
  {
    return std::accumulate(vars.begin(), vars.end(), 0,
      [](int acc, const tvm::VariablePtr & v){ return acc + static_cast<int>(v->space().tSize()); });
  };

  const int cols_ddq = static_cast<int>(J_ddq.cols());
  const int cols_dq  = static_cast<int>(J_dq.cols());
  const int size_ddq = sumSizes(q_ddot_vars_);
  const int size_dq  = sumSizes(q_dot_vars_);

  if(cols_ddq != size_ddq)
  {
    mc_rtc::log::error_and_throw("[ImpulseFunction] J_ddq cols {} != sum(vars) {}", cols_ddq, size_ddq);
  }
  if(cols_dq != size_dq)
  {
    mc_rtc::log::error_and_throw("[ImpulseFunction] J_dq cols {} != sum(vars) {}", cols_dq, size_dq);
  }

  //
  splitJacobian(J_ddq, q_ddot_vars_);
  splitJacobian(J_dq, q_dot_vars_);

  // mc_rtc::log::info(" actual valuation:");
  // mc_rtc::log::info(value());
  // mc_rtc::log::info("lower bounds are");
  // mc_rtc::log::info(robot_.tvmRobot().limits().tl * limit_multiplier_ / (c_res_+1));
  // mc_rtc::log::info("upper bounds are");
  // mc_rtc::log::info(robot_.tvmRobot().limits().tu * limit_multiplier_ / (c_res_+1));

  // splitJacobian(robot.H(), robot.alphaD());
  // splitJacobian(robot.H(), robot.alpha());
  // J_ddq = Eigen::MatrixXd::Identity(robot_.mb().nrDof(), robot_.mb().nrDof());
  // J_dq = Eigen::MatrixXd::Identity(robot_.mb().nrDof(), robot_.mb().nrDof());

  // auto n = static_cast<Eigen::DenseIndex>(tvm::dot(robot.q(), 2)->space().tSize());
  // mc_rtc::log::info("size of n is {} and there are {} available columns", n, J_ddq.cols());
  // auto n2 = static_cast<Eigen::DenseIndex>(tvm::dot(robot.q(), 1)->space().tSize());
  // mc_rtc::log::info("size of n is {} and there are {} available columns", n2, J_dq.cols());

  // splitJacobian(J_ddq, q_ddot_var_);
  // splitJacobian(J_dq, q_dot_var_);

  // mc_rtc::log::info("size of J is {}x{}", full_world_frame_jacobian.rows(), full_world_frame_jacobian.cols());
  // mc_rtc::log::info("size of M is {}x{}", M.rows(), M.cols());
  // mc_rtc::log::info("size of Pn is {}x{}", P_n.rows(), P_n.cols());
  // mc_rtc::log::info("size of j_m is {}x{}", j_m.rows(), j_m.cols());
  // mc_rtc::log::info("size of J_ddq is {}x{}", J_ddq.rows(), J_ddq.cols());
  // mc_rtc::log::info("size of J_dq is {}x{}", J_dq.rows(), J_dq.cols());

  // mc_rtc::log::info("the evaluation of the function:");
  // mc_rtc::log::info(value());

}
} // namespace mc_tvm


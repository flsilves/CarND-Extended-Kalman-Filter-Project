#include "kalman_filter.h"

#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  const MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

  // Predicted measurement error
  const VectorXd z_pred = H_ * x_;
  const VectorXd y = z - z_pred;

  // Kalman Gain
  const MatrixXd Ht = H_.transpose();
  const MatrixXd S = H_ * P_ * Ht + R_;
  const MatrixXd Si = S.inverse();
  const MatrixXd PHt = P_ * Ht;
  const MatrixXd K = PHt * Si;

  // Update
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  const MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

  // Predicted measurement error
  const VectorXd z_pred = Tools::ConvertCartesianToPolarCoords(x_);
  VectorXd y = z - z_pred;
  y(1) = Tools::NormalizeAngle(y(1));

  // Kalman Gain
  const MatrixXd Ht = H_.transpose();
  const MatrixXd S = H_ * P_ * Ht + R_;
  const MatrixXd Si = S.inverse();
  const MatrixXd PHt = P_ * Ht;
  const MatrixXd K = PHt * Si;

  // Update
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;
}

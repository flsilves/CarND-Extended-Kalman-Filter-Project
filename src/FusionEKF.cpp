#include "FusionEKF.h"

#include <iostream>

#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

FusionEKF::FusionEKF() {
  // measurement function - laser
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,  //
      0, 1, 0, 0;

  // measurement covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,  //
      0, 0.0225;

  // measurement function - radar (computed each step)
  H_radar_ = MatrixXd(3, 4);

  // measurement covariance matrix - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,  //
      0, 0.0009, 0,        //
      0, 0, 0.09;

  // State transition
  ekf_.F_ = Eigen::MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,  //
      0, 1, 0, 1,         //
      0, 0, 1, 0,         //
      0, 0, 0, 1;

  // State Covariance matrix
  // Initialization with low error for position and higher error for velocity
  ekf_.P_ = Eigen::MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,  //
      0, 1, 0, 0,         //
      0, 0, 1000, 0,      //
      0, 0, 0, 1000;      //

  // State transition noise
  ekf_.Q_ = Eigen::MatrixXd::Identity(4, 4);
}

FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  if (!is_initialized_) {
    InitializeProcess(measurement_pack);
  } else {
    Prediction(measurement_pack.timestamp_);
    Update(measurement_pack.raw_measurements_, measurement_pack.sensor_type_);
  }

  UpdateLastMeasurementTime(measurement_pack.timestamp_);

  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

void FusionEKF::UpdateLastMeasurementTime(long long timestamp) {
  previous_timestamp_ = timestamp;
}

void FusionEKF::Prediction(long long timestamp) {
  // elapsed time
  float dt = (timestamp - previous_timestamp_) / 1000000.0;
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  ekf_.Q_ = Eigen::MatrixXd(4, 4);
  ekf_.Q_ << dt_4 / 4 * noise_ax_, 0, dt_3 / 2 * noise_ax_, 0, 0,
      dt_4 / 4 * noise_ay_, 0, dt_3 / 2 * noise_ay_, dt_3 / 2 * noise_ax_, 0,
      dt_2 * noise_ax_, 0, 0, dt_3 / 2 * noise_ay_, 0, dt_2 * noise_ay_;

  ekf_.Predict();
  previous_timestamp_ = timestamp;
}

void FusionEKF::InitializeProcess(const MeasurementPackage &measurement_pack) {
  switch (measurement_pack.sensor_type_) {
    case MeasurementPackage::LASER:
      InitializeProcessWithLaserData(measurement_pack.raw_measurements_[0],
                                     measurement_pack.raw_measurements_[1]);
      break;
    case MeasurementPackage::RADAR:
      InitializeProcessWithRadarData(measurement_pack.raw_measurements_[0],
                                     measurement_pack.raw_measurements_[1],
                                     measurement_pack.raw_measurements_[2]);
      break;
  }
}

void FusionEKF::Update(Eigen::VectorXd z,
                       MeasurementPackage::SensorType sensor_type) {
  switch (sensor_type) {
    case MeasurementPackage::SensorType::RADAR:
      H_radar_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.H_ = H_radar_;
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(z);
      break;
    case MeasurementPackage::SensorType::LASER:
      ekf_.H_ = H_laser_;
      ekf_.R_ = R_laser_;
      ekf_.Update(z);
      break;
  }
}

void FusionEKF::InitializeProcessWithLaserData(float position_x,
                                               float position_y) {
  ekf_.x_ = Eigen::VectorXd(4);
  ekf_.x_ << position_x, position_y, 0, 0;

  is_initialized_ = true;
}

void FusionEKF::InitializeProcessWithRadarData(float range, float angle,
                                               float radial_velocity) {
  const float px = range * cosf(angle);
  const float py = range * sinf(angle);
  const float vx = radial_velocity * cosf(angle);
  const float vy = radial_velocity * sinf(angle);

  ekf_.x_ = Eigen::VectorXd(4);
  ekf_.x_ << px, py, vx, vy;

  is_initialized_ = true;
}

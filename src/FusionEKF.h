#ifndef FusionEKF_H_
#define FusionEKF_H_

#include <fstream>
#include <string>
#include <vector>

#include "Eigen/Dense"
#include "kalman_filter.h"
#include "measurement_package.h"
#include "tools.h"

class FusionEKF {
 public:
  FusionEKF();
  virtual ~FusionEKF();

  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  KalmanFilter ekf_;

 private:
  bool is_initialized_{false};

  void InitializeProcess(const MeasurementPackage &measurement_pack);

  void InitializeProcessWithLaserData(float px, float py);

  void InitializeProcessWithRadarData(float range, float angle,
                                      float radial_velocity);

  void Prediction(long long timestamp);

  void Update(Eigen::VectorXd z, MeasurementPackage::SensorType type);

  void UpdateLastMeasurementTime(long long timestamp);
  long long previous_timestamp_{0ll};

  Tools tools;
  Eigen::MatrixXd R_laser_;
  Eigen::MatrixXd R_radar_;
  Eigen::MatrixXd H_laser_;
  Eigen::MatrixXd H_radar_;

  float noise_ax_{9.0f};
  float noise_ay_{9.0f};
};

#endif  // FusionEKF_H_

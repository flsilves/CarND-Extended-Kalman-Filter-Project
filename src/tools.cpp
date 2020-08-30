#include "tools.h"

#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  const auto n_estimations = estimations.size();
  rmse << 0, 0, 0, 0;

  if (n_estimations == 0 || (n_estimations != ground_truth.size())) {
    std::cerr << "Invalid input data for RMSE calculation" << std::endl;
    return rmse;
  }

  for (size_t i = 0; i < n_estimations; ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }
  rmse = rmse / n_estimations;
  rmse = rmse.array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state) {
  MatrixXd Hj(3, 4);

  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float c1 = px * px + py * py;
  float c2 = sqrt(c1);
  float c3 = (c1 * c2);

  // check division by zero
  if (c1 < 0.000001) {
    c1 = 0.000001f;
  }

  float Hj_11 = px / c2;
  float Hj_12 = py / c2;
  float Hj_13 = 0;
  float Hj_14 = 0;

  float Hj_21 = -1 * py / c1;
  float Hj_22 = px / c1;
  float Hj_23 = 0;
  float Hj_24 = 0;

  float Hj_31 = py * (vx * py - vy * px) / c3;
  float Hj_32 = px * (px * vy - py * vx) / c3;
  float Hj_33 = Hj_11;
  float Hj_34 = Hj_12;

  Hj << Hj_11, Hj_12, Hj_13, Hj_14,  //
      Hj_21, Hj_22, Hj_23, Hj_24,    //
      Hj_31, Hj_32, Hj_33, Hj_34;    //

  return Hj;
}

float Tools::NormalizeAngle(float value) {
  double a = fmod(value + M_PI, 2 * M_PI);
  return a >= 0 ? (a - M_PI) : (a + M_PI);
}

Eigen::VectorXd Tools::ConvertCartesianToPolarCoords(const Eigen::VectorXd &z) {
  const float px = z[0];
  const float py = z[1];
  const float vx = z[2];
  const float vy = z[3];

  float rho = sqrtf(px * px + py * py);
  float phi = atan2f(py, px);
  float rho_dot = 0;

  if (abs(rho) < 0.000001) {  // Avoid dividing by zero
    rho = 0.000001;
  }

  rho_dot = (px * vx + py * vy) / rho;

  Eigen::VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot;
  return z_pred;
}
#include "tools.h"

#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd x(4);
  return x;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state) { return {}; }

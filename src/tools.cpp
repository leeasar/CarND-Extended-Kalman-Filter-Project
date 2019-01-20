#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // Check the estimation vector size is non-zero and that vector sizes match
  if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
    return rmse;
  }

  // Accumulate squared residuals
  for (int i = 0; i < estimations.size(); i++) {

    VectorXd residual = estimations[i] - ground_truth[i];

    // Residuals squared
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  // Calculating the mean
  rmse = rmse / estimations.size();

  // Calculating the squared root
  rmse = rmse.array().sqrt();

  // Return an array of RMSE results
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  MatrixXd Hj(3,4);

  // Get state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  // Pre-compute a set of terms to avoid repeated calculation
  double c1 = px * px + py * py;
  double c2 = sqrt(c1);
  double c3 = (c1 * c2);

  // Check division by zero or close to zero
  if (c1 >= 0 && c1 < 0.0001) { c1 = 0.0001; }
  if (c1 < 0 && c1 > -0.0001) { c1 = -0.0001; }

  // Compute the Jacobian matrix
  Hj << (px/c2), (py/c2), 0, 0,
      -(py/c1), (px/c1), 0, 0,
      py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;
}

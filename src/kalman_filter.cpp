#include "kalman_filter.h"
#include <iostream>


using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

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
  /**
   * Predict the state x
   */

   x_ = F_ * x_;
   MatrixXd Ft = F_.transpose();
   P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * Update function for laser (use KF)
   */
   VectorXd z_pred = H_ * x_;
   VectorXd y = z - z_pred;
   MatrixXd Ht = H_.transpose();
   MatrixXd S = H_ * P_ * Ht + R_;
   MatrixXd Si = S.inverse();
   MatrixXd PHt = P_ * Ht;
   MatrixXd K = PHt * Si;

   // New estimate for state x
   x_ = x_ + (K * y);

   // New state covariance matrix (uncertainty about position and velocity)
   int x_size = x_.size();
   MatrixXd I = MatrixXd::Identity(x_size, x_size);
   P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * Update function for radar (use EKF)
   * Otherwise this is same as KF, but the process to calculate y is more complex
   */

   float px = x_(0);
   float py = x_(1);
   float vx = x_(2);
   float vy = x_(3);

   float rho = sqrt(px * px + py * py);
   float theta = atan2(py, px);
   float rho_dot = (px * vx + py * vy) / rho;
   VectorXd h = VectorXd(3);

   // Use a theta between PI and -PI
   while (theta > M_PI || theta < -M_PI) {
     if (theta > M_PI) { theta -= M_PI; }
     else { theta += M_PI;}
    }

   h << rho, theta, rho_dot;
   VectorXd y = z - h;

   MatrixXd Ht = H_.transpose();
   MatrixXd S = H_ * P_ * Ht + R_;
   MatrixXd Si = S.inverse();
   MatrixXd PHt = P_ * Ht;
   MatrixXd K = PHt * Si;

   // New estimate for state x
   x_ = x_ + (K * y);

   // New state covariance matrix (uncertainty about position and velocity)
   int x_size = x_.size();
   MatrixXd I = MatrixXd::Identity(x_size, x_size);
   P_ = (I - K * H_) * P_;
}

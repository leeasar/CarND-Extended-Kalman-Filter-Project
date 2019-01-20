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
   * Otherwise this is same as KF, but the process to calculate z_pred is more complex
   */

   double px = x_(0);
   double py = x_(1);
   double vx = x_(2);
   double vy = x_(3);

   double rho = sqrt(px * px + py * py);
   if (rho >= 0.0 && rho < 0.0001) { rho = 0.0001; }
   if (rho < 0.0 && rho > -0.0001) { rho = -0.0001; }

   double phi = atan2(py, px);
   double rho_dot = (rho != 0 ? (px * vx + py * vy) / rho : 0);

   VectorXd z_pred = VectorXd(3);
   z_pred << rho, phi, rho_dot;
   VectorXd y = z - z_pred;

   // Use a phi between PI and -PI
   while (y(1) > M_PI || y(1) < -M_PI) {
     if (y(1) > M_PI) { y(1) -= 2*M_PI; }
     else { y(1) += 2*M_PI;}
    }

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

#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  //measurement matrix - laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */

  //state covariance matrix
  MatrixXd P_ = MatrixXd(4, 4);
  P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;

  //state transition matrix
  MatrixXd F_ = MatrixXd(4, 4);
  F_ << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;

  //covariance matrix (process noise)
  MatrixXd Q_ = MatrixXd(4, 4);
  Q_ << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;

  //state
  VectorXd x_ = VectorXd(4);
  x_ << 1, 1, 1, 1;

  //initializing the Kalman filter
  ekf_.Init(x_, P_, F_, H_laser_, R_laser_, Q_);

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize state
      ekf_.x_ << measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]),
                measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]),
                measurement_pack.raw_measurements_[2] * cos(measurement_pack.raw_measurements_[1]),
                measurement_pack.raw_measurements_[2] * sin(measurement_pack.raw_measurements_[1]);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state with the initial location and zero velocity
      ekf_.x_ << measurement_pack.raw_measurements_[0],
                measurement_pack.raw_measurements_[1],
                5.25, // approximate based on this specific data file
                0.25; // approximate based on this specific data file
    }

    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }


  /**
   * Prediction
   */

  // Using values 9 and 9 as instructed
  float noise_ax = 9.0;
  float noise_ay = 9.0;

  // Computing the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  // Computing temporary variables for computing Q
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  // Modify the state transition matrix F with the new time stamps
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // Setting the process covariance matrix Q based on the current dt
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
         0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
         dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
         0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_); // Use Jacobian measurement matrix for radar
    ekf_.R_ = R_radar_; // Use the measurement covariance matrix for radar
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
    ekf_.H_ = H_laser_; // Use the measurement matrix for laser
    ekf_.R_ = R_laser_; // Use the measurement covariaance matrix for laser
    ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  std::cout << "x_ = " << ekf_.x_ << std::endl;
  std::cout << "P_ = " << ekf_.P_ << std::endl;
}

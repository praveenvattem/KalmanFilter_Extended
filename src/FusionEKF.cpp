#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd::Zero(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
 
  //State Transition Matrix
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1,0,1,0,
             0,1,0,1,
             0,0,1,0,
             0,0,0,1;
  
  H_laser_<< 1,0,0,0,
	    0,1,0,0;

  //Process covariance matrix
  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ << 1,0,0,0,
	     0,1,0,0,
	     0,0,1000,0,
	     0,0,0,1000;

  //Process noise covariance matrix
  ekf_.Q_ = MatrixXd::Zero(4, 4);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      cout << "Initial Measurement: Radar"<< endl;

      // rho and phi are the range and angle of the radar measurement
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];

      double px = rho*cos(phi);
      double py = rho*sin(phi);

      if(px < 0.0001){px = 0.0001;}
      if(py < 0.0001){py = 0.0001;}

      ekf_.x_ << px, py, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      cout << "Initial Measurement: Laser."<< endl;
      ekf_.x_<< measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    //Timestamp of first measurement
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  //Calculate elapsed time
  double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  // Modify the state transition matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;


  //Acceleration noise components
  double noise_ax = 9.0;
  double noise_ay = 9.0;

  // Set the process noise covariance matrix Q
  ekf_.Q_(0,0) = 0.25*noise_ax*pow(dt,4);
  ekf_.Q_(0,2) = 0.5*pow(dt,3)*noise_ax;
  ekf_.Q_(1,1) = 0.25*pow(dt,4)*noise_ay;
  ekf_.Q_(1,3) = 0.5*pow(dt,3)*noise_ay;
  ekf_.Q_(2,0) = 0.5*pow(dt,3)*noise_ax;
  ekf_.Q_(2,2) = pow(dt,2)*noise_ax;
  ekf_.Q_(3,1) = 0.5*pow(dt,3)*noise_ay;
  ekf_.Q_(3,3) = pow(dt,2)*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }
  else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

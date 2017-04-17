#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = .8; //originally 30 now set to ~3 for bicycle acceleration

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODONE:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = n_x_ + 2;

  /// Sigma Point Matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  ///* Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_ + n_aug_);
  weights_(0) = weight_0;
  for(int i=1; i<2*n_aug_ + 1; i++) {
    double weight = 0.5/(n_aug_ + lambda_);
    weights_(i) = weight;
  }

  // ///* the current NIS for radar
  NIS_radar_ = 0.0;
  //
  // ///* the current NIS for laser
  NIS_laser_ = 0.0;
}

//Counter to keep track of iterations
int counter = 0;

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODONE:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  //std::cout << "Processing Measurement"<< std::endl;

  //ignore zero measurements
  // if (meas_package.raw_measurements_[0] == 0 || meas_package.raw_measurements_[1] == 0) {
  //   std::cout << "Ignoring measurement" << std::endl;
  //   return;
  // }

  if(!is_initialized_) {

    // first measurement
    std::cout << "UKF: " << std::endl;
    x_ << 1.0, 1.0, 1.0, 1.0, 1.0;
    //std::cout << "x_ initial" << x_ << std::endl;
    //state covariance matrix P
    P_ <<  1.0,0.0,0.0,0.0,0.0,
           0.0,1.0,0.0,0.0,0.0,
           0.0,0.0,1.0,0.0,0.0,
           0.0,0.0,0.0,1.0,0.0,
           0.0,0.0,0.0,0.0,1.0;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double px = 0.0;
      double py = 0.0;
      double rho = meas_package.raw_measurements_[0];
      double theta = meas_package.raw_measurements_[1];
      px = rho * cos(theta);
      py = rho * sin(theta);
      // std::cout << "px: " << px << std::endl;
      // std::cout << "py: " << py << std::endl;
      x_ << px, py, 0.0, 0.0, 0.0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      //std::cout << "Laser Measurement" << std::endl;
      double x = meas_package.raw_measurements_[0];
      double y = meas_package.raw_measurements_[1];
      // std::cout << "px: " << x << std::endl;
      // std::cout << "py: " << y << std::endl;
      x_ << x, y, 0.0, 0.0, 0.0;
    }
    // set the timestamp
    time_us_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    std::cout << "Initialized" << std::endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

   //compute the time elapsed between the current and previous measurements
	float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  //make sure dt is greater than zero and do update and Predict
  //std::cout << "dt " << dt << std::endl;
  if ( dt < 0.0001 )
  {
    //dt = 0.0001;
    return;
  }
  while (dt > 0.1)
  {
    const double delta_t = 0.05;
    Prediction(delta_t);
    dt -= delta_t;
  }
  Prediction(dt);

  //update the timestamp
  time_us_ = meas_package.timestamp_;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODONE:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
    // Laser updates
    UpdateLidar(meas_package);
  }

//   if(counter < 100) {
//   // print the output
//   std::cout << "x_ = " << x_ << std::endl;
//   std::cout << "P_ = " << P_ << std::endl;
// }

  //iterate counter
  counter++;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODONE:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //std::cout << "Predicting"<< std::endl;

  MatrixXd Xsig = MatrixXd(11, 5);
  Xsig.fill(0.0);
  GenerateSigmaPoints(&Xsig);
  MatrixXd Xsig_aug = MatrixXd(15, 7);
  Xsig_aug.fill(0.0);
  AugmentSigmaPoints(&Xsig_aug);
  SigmaPointPredict(Xsig_aug, delta_t);
  PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODONE:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  //std::cout << "Begin Lidar Update" << std::endl;
  int n_z = 2;

  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd::Zero(n_z, n_z);
  S.fill(0.0);
  MatrixXd Zsig = MatrixXd::Zero(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  PredictLidarMeasurement(&z_pred, &S, &Zsig);

  //create vector for incoming laser measurement
  VectorXd z = VectorXd(n_z);
  z <<
      meas_package.raw_measurements_[0],
      meas_package.raw_measurements_[1];

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //std::cout << "z_diff" << z_diff << std::endl;
    //angle normalization
    // if(z_diff(1) > M_PI){
    //     z_diff(1) = (int(z_diff(1) - M_PI)%int(2*M_PI)) - M_PI;
    //   }
    // if(z_diff(1) < -M_PI){
    //   z_diff(1) = (int(z_diff(1) + M_PI)%int(2*M_PI)) + M_PI;
    // }
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    // if(x_diff(3) > M_PI){
    //     x_diff(3) = (int(x_diff(3) - M_PI)%int(2*M_PI)) - M_PI;
    //   }
    // if(x_diff(3) < -M_PI){
    //   x_diff(3) = (int(x_diff(3) + M_PI)%int(2*M_PI)) + M_PI;
    // }
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  // std::cout << "Kalman Gain" << std::endl;
  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;
  // std::cout << "Angle Normalization" << std::endl;
  //angle normalization
  // if(z_diff(1) > M_PI){
  //   z_diff(1) = (int(z_diff(1) - M_PI)%int(2*M_PI)) - M_PI;
  // }
  // if(z_diff(1) < -M_PI){
  //   z_diff(1) = (int(z_diff(1) + M_PI)%int(2*M_PI)) + M_PI;
  // }
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  // std::cout << "Update State and Mean Covariance Matrix" << std::endl;
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  // std::cout << "Update NIS Lidar" << std::endl;
  //update NIS_Lidar
  UpdateNISLidar(meas_package, z_pred , S);
  if(counter == 1) {
  //print result
  //std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  //std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
}
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODONE:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //std::cout << "Begin Radar Update" << std::endl;
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  VectorXd z_pred = VectorXd(3);
  MatrixXd S = MatrixXd(3, 3);
  S.fill(0.0);
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  PredictRadarMeasurement(&z_pred, &S, &Zsig );

  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z <<
      meas_package.raw_measurements_[0],
      meas_package.raw_measurements_[1],
      meas_package.raw_measurements_[2];

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    // if(z_diff(1) > M_PI){
    //   z_diff(1) = (int(z_diff(1) - M_PI)%int(2*M_PI)) - M_PI;
    // }
    // if(z_diff(1) < -M_PI){
    //   z_diff(1) = (int(z_diff(1) + M_PI)%int(2*M_PI)) + M_PI;
    // }
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    // if(x_diff(3) > M_PI){
    //   x_diff(3) = (int(x_diff(3) - M_PI)%int(2*M_PI)) - M_PI;
    // }
    // if(x_diff(3) < -M_PI){
    //   x_diff(3) = (int(x_diff(3) + M_PI)%int(2*M_PI)) + M_PI;
    // }
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  // if(z_diff(1) > M_PI){
  //   z_diff(1) = (int(z_diff(1) - M_PI)%int(2*M_PI)) - M_PI;
  // }
  // if(z_diff(1) < -M_PI){
  //   z_diff(1) = (int(z_diff(1) + M_PI)%int(2*M_PI)) + M_PI;
  // }
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  //update NIS_Radar
  UpdateNISRadar(meas_package, z_pred , S);

//   if(counter < 3) {
//   //print result
//   std::cout << "Updated state x: " << std::endl << x_ << std::endl;
//   std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
// }
//counter++;
}

void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {
  if(counter == 1){
  //std::cout << "Begin Generate Sigma Points" << std::endl;
  //std::cout << x_ << std::endl;
}
  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
  Xsig.fill(0.0);
  if(counter == 0){
  std::cout << Xsig << std::endl;
}
  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();

  //set first column of sigma point matrix
  Xsig.col(0)  = x_;
//   if(counter == 1) {
//   std::cout << Xsig.col(0) << std::endl;
// }
  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig.col(i+1)     = x_ + sqrt(lambda_+ n_x_) * A.col(i);
    Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_+ n_x_) * A.col(i);
  }
  if(counter == 1) {
  //print result
  //std::cout << "Generated Sigma Points Xsig = " << std::endl << Xsig << std::endl;
}
  *Xsig_out = Xsig;
}

void UKF::AugmentSigmaPoints(MatrixXd* Xsig_out) {

  //std::cout << "Begin Augment Sigma Points " << std::endl;

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(n_x_,n_x_) = std_a_*std_a_;
  P_aug(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  if(counter == 1) {
  //print result
  //std::cout << "Augmented Sigma Pts Xsig_pred = " << std::endl << Xsig_aug << std::endl;
  }
  *Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPredict(MatrixXd Xsig_aug, double delta_t) {
  //std::cout << "Begin Sigma Point Predict" << std::endl;
  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //std::cout << "Writing values to Xsig" << std::endl;
    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  if(counter == 1) {
  //print result
  //std::cout << "Xsig_pred = " << std::endl << Xsig_pred_ << std::endl;
}
}

void UKF::PredictMeanAndCovariance() {
  //std::cout << "predict mean and covariance" << std::endl;
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
    //std::cout << x_ << std::endl;
  }
  //std::cout << "here 1" << n_aug_ << std::endl;
  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    //std::cout << "here 2" << std::endl;
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //std::cout << "x_diff" << x_diff << std::endl;
    //angle normalization
    if(x_diff(3) > M_PI){
        x_diff(3) = (int(x_diff(3) - M_PI)%int(2*M_PI)) - M_PI;
      }
    if(x_diff(3) < -M_PI){
      x_diff(3) = (int(x_diff(3) + M_PI)%int(2*M_PI)) + M_PI;
    }
    // while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    // while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    //std::cout << "here 2a" << std::endl;
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
  //std::cout << "here 3" << std::endl;

  if(counter == 1) {
  //print result
  // std::cout << "Predicted state" << std::endl;
  // std::cout << x_ << std::endl;
  // std::cout << "Predicted covariance matrix" << std::endl;
  // std::cout << P_ << std::endl;
}
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Z_sig_out) {
  //std::cout << "Begin Radar Measurement" << std::endl;
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    if (p_x == 0 && p_y == 0)
    {
      Zsig(0,i) = 0;
      Zsig(1,i) = 0;
      Zsig(2,i) = 0;
    }
    else
    {
      Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
      Zsig(1,i) = atan2(p_y,p_x);                                 //phi
      Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    if(z_diff(1) > M_PI){
      z_diff(1) = (int(z_diff(1) - M_PI)%int(2*M_PI)) - M_PI;
    }
    if(z_diff(1) < -M_PI){
      z_diff(1) = (int(z_diff(1) + M_PI)%int(2*M_PI)) + M_PI;
    }
    // while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    // while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

  //print result
  // std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  // std::cout << "S: " << std::endl << S << std::endl;

  //write result
  *z_out = z_pred;
  *S_out = S;
  *Z_sig_out = Zsig;
}

void UKF::PredictLidarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Z_sig_out) {

  //std::cout << "Begin Predict Lidar Measurement" << std::endl;
  //set measurement dimension, lidar can measure px,py
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;


      Zsig(0,i) = p_x;
      Zsig(1,i) = p_y;
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //std::cout << "z_diff" << z_diff << std::endl;
    //angle normalization
    if(z_diff(1) > M_PI){
      z_diff(1) = (int(z_diff(1) - M_PI)%int(2*M_PI)) - M_PI;
    }
    if(z_diff(1) < -M_PI){
      z_diff(1) = (int(z_diff(1) + M_PI)%int(2*M_PI)) + M_PI;
    }
    // while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    // while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  S = S + R;

//   if(counter < 2) {
//   //print result
//   std::cout << "z_pred: " << std::endl << z_pred << std::endl;
//   std::cout << "S: " << std::endl << S << std::endl;
//   std:cout << "Zsig: " << std::endl << Zsig << std::endl;
// } counter++;

  //write result
  *z_out = z_pred;
  *S_out = S;
  *Z_sig_out = Zsig;
}

void UKF::UpdateNISRadar(MeasurementPackage meas_package, VectorXd z_pred, MatrixXd S) {
  // update NIS
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  NIS_radar_ = z_diff.transpose()*S.inverse()*z_diff;}

void UKF::UpdateNISLidar(MeasurementPackage meas_package, VectorXd z_pred, MatrixXd S) {
  //std::cout << "Begin Update NIS Lidar" << endl;
  // update NIS
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  NIS_laser_ = z_diff.transpose()*S.inverse()*z_diff;
  //std::cout << "NIS Lidar Value" << result << endl;
}

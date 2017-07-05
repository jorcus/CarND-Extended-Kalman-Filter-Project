#include "kalman_filter.h"
#include <math.h>       /* pow */

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  // Kalman Filter Variable
  x_ = x_in;  // Object State
  P_ = P_in;  // Object covariance Matrix
  F_ = F_in;  // State Transition Matrix
  H_ = H_in;  // Measurement Matrix
  R_ = R_in;  // Measurement Covariance Matrix
  Q_ = Q_in;  // Process Covariance Matrix
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */

  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K =  PHt * Si;

  //new state
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  

  float rho = sqrt(px*px + py*py);
  float theta = atan2(py, px);
  float ro_dot = (x_(0)*x_(2) + x_(1)*x_(3))/rho;
  if (fabs(rho) < 0.001) {ro_dot = 0;}

  VectorXd z_pred(3);
  z_pred << rho, theta, ro_dot;

  VectorXd y = z - z_pred;
  
  //angle normalization
  while (y(1)> M_PI) y(1)-=2.*M_PI;
  while (y(1)<-M_PI) y(1)+=2.*M_PI;
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // New State
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

  /*
  //KF Prediction step
  x_ = F_ * x_ + u;
  MatrixXd Ft = F.transpose();
  P_ = F_ * P_ * Ft + Q_;
  */
}
#include <iostream>
#include "tools.h"
#include <cmath>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() == 0 || estimations.size() != ground_truth.size() || ground_truth.size() ==0){
    cout << "Estimation or ground truth vector size is incorrect";
    return rmse;
  }

  //accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); ++i){
    VectorXd residuals = (estimations[i] - ground_truth[i]);
    residuals = residuals.array()*residuals.array();
    rmse += residuals;
  }

  rmse /= estimations.size();
  rmse = rmse.array().sqrt(); 
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
  //recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  double rho_squared = px*px + py*py;
  double rho = sqrt(rho_squared);
  
  /*
  //check if px or py are small
  if (px < 0.0001){px = 0.0001;}
  
  if (py < 0.0001){py = 0.0001;}
  */

  //check division by zero
  if (rho_squared < 0.0001){
    cout << "Jacobian - error division by zero" << endl;
    return Hj;
   }

  //compute the Jacobian matrix
  Hj << px/rho, py/rho, 0, 0,
       -py/(rho*rho), px/(rho*rho), 0, 0,
        py*(vx*py-vy*px)/(rho*rho*rho), px*(vy*px-vx*py)/(rho*rho*rho), px/rho, py/rho;
    
  return Hj; 
}

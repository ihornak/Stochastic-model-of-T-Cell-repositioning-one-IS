//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS
#ifndef POINT_ABSCISSA_H
#define POINT_ABSCISSA_H
#include <iostream>
#include <stdio.h>     
#include <math.h>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/SparseLU>
using Eigen::MatrixXd;
using namespace Eigen;

using namespace std;

namespace Dynein
{
	Vector3d point;
	double abscissa;
}

#endif 

#ifndef IS_CORTICAL_SL_PARAMETER_H_
#define IS_CORTICAL_SL_PARAMETER_H_
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/SparseLU>

using Eigen::MatrixXd;
using namespace Eigen;

//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS


namespace IS_Cortical_Sl_parameter
{
	const unsigned int number_of_polygon = 12;
	const unsigned int number_of_mito = 5;    
	const double force_per_lenght = 5e-6;
	const double force = 50e-12;

    const double replace_boundary =2.5e-7;
    const double catching_radius = 6.8e-6;
    const double zero_force_boundary = 3e-7;
    const double ellipsoid_control = 0.05;
    const Vector3d orientation( 0.0 , 0.0 , - 1.0 );
    
    const double front_radius = 3.8e-6;
    const double rear_radius = 5.5e-6;
    const double radius = 2e-6;
    const double radius_inner = 0;

    
}




#endif /* IS_CAPTURE_SHRINKAGE_PARAMETERS_H_ */

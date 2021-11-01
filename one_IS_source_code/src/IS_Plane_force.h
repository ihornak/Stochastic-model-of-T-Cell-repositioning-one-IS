//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS

#ifndef IS_PLANE_FORCE_H_
#define IS_PLANE_FORCE_H_
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
using Eigen::MatrixXd;


namespace IS_Plane_force
{
	const unsigned int switch_plane = 1;
	const double plane_force_k1 = 1e-4;//1e-4
	const double plane_force_k2 = 4.0;
	const double z_values = - 3.8e-6;
	const double x = 0.0;
	const double y = 0.0;
	const double z = 1.0;
}



#endif /* IS_PLANE_FORCE_H_ */

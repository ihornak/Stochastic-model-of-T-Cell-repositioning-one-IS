//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS

#ifndef NUCLEUS_H_
#define NUCLEUS_H_
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include "Nucleus_parametres.h"
#include "Cell_parametres.h"


using Eigen::MatrixXd;
using namespace Eigen;

class Nucleus {
public:
	Nucleus();
	Nucleus( Vector3d center_argument );
	Nucleus& operator=( const Nucleus &tmp );
	double get_A_Axis() const;
	double get_B_Axis() const;
	Vector3d get_center();
	Vector3d force_position_nucleus_wall( Vector3d position );
    	Vector3d force_position_nucleus_wall_MTOC( Vector3d position );

	virtual ~Nucleus();

private:
	double a_axis;
	double b_axis;
	Vector3d center;
};

#endif /* NUCLEUS_H_ */

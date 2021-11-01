//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS
#ifndef ISCAPTSHRINKAGE_H_
#define ISCAPTSHRINKAGE_H_
#include <eigen3/Eigen/Dense>
using Eigen::MatrixXd;
using namespace Eigen;
#include <iostream>
#include <string>
#include "GeneralGeometry.h"

using namespace std;


class IS_Capt_Shrinkage {
public:
	IS_Capt_Shrinkage();
	IS_Capt_Shrinkage( Vector3d center_of_IS_front , double radius_argument , Vector3d center_of_IS_rear );
	IS_Capt_Shrinkage( Vector3d center_of_IS_front , double radius_argument , Vector3d center_of_IS_rear , Vector3d trap_arg );
	IS_Capt_Shrinkage( IS_Capt_Shrinkage& tmp );
	IS_Capt_Shrinkage operator=( const IS_Capt_Shrinkage& tmp );
	Vector3d get_center_of_IS_front() const;
	Vector3d get_center_of_IS_rear() const;
	double get_radius_of_IS() const;
	Vector3d get_axis_of_IS() const;
	void set_center_of_IS_front( Vector3d center_argument_front );
	void set_center_of_IS_rear( Vector3d center_argument_rear );
	void set_radius_of_IS( double radius_center );
	void set_axis_of_IS( Vector3d axis_argument );
	bool check_IS_capture_shrinkage_caught_Main( Vector3d position );
        bool check_IS_capture_segment_control( Vector3d first_point , Vector3d tangent );
	virtual ~IS_Capt_Shrinkage();

private:
	Vector3d center_of_IS_front;
	Vector3d center_of_IS_rear;
	Vector3d axis_of_IS;
	double radius_of_IS;
	
	//MatrixXd trap;





};

#endif /* ISPLANAR_H_ */

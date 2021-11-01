//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS

#ifndef GENERALGEOMETRY_H_
#define GENERALGEOMETRY_H_

#include <eigen3/Eigen/Dense>
using Eigen::MatrixXd;
using namespace Eigen;
#include "Cell_parametres.h"
#include <vector>

#include <iostream>
#include <string>
using namespace std;



double distance_two_lines( Vector3d a_vector , Vector3d b_vector , Vector3d c_vector , Vector3d d_vector );
double distance_two_line_segments( Vector3d a_vector , Vector3d b_vector , Vector3d c_vector , Vector3d d_vector );
double distance_two_line_segments( Vector3d a_vector , Vector3d b_vector , Vector3d c_vector , Vector3d d_vector , unsigned int interaction_index);
double distance_two_line_segments( Vector3d a_vector , Vector3d b_vector , Vector3d c_vector , Vector3d d_vector , Vector3d& p_1 , Vector3d& p_2 , double& t_return, double& v_return );
double distance_point_segment( Vector3d a_vector , Vector3d b_vector , Vector3d y_vector );
double distance_point_segment( Vector3d a_vector , Vector3d b_vector , Vector3d y_vector , Vector3d& closest_point_of_segment );
double distance_plane_point( Vector3d plane_axis , Vector3d point_on_plane , Vector3d point );
double distance_point_line( Vector3d first_line_point , Vector3d second_line_point , Vector3d y_vector );
Vector3d project_point_on_surface_of_elipse( Vector3d position );

double coordinate_norm( MatrixXd coordinates );




Vector3d projection_of_tangent_on_plane( Vector3d plane_axis , Vector3d tangent );
Vector3d projection_of_point_on_plane( Vector3d plane_axis , Vector3d point_on_plane , Vector3d point_to_be_projected );
Vector3d projection_of_point_on_eclipse_plane_in_point( double A_Axis , double B_axis , Vector3d point_on_surface, Vector3d point_to_be_projected );



double magnitude_3D_matrix( MatrixXd& matrix);



#endif /* GENERALGEOMETRY_H_ */

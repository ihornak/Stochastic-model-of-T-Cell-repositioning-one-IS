#ifndef ISCORTICALSL_H_
#define ISCORTICALSL_H_
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/IterativeLinearSolvers>
//#include <unsupported/Eigen/IterativeSolvers>
#include <eigen3/Eigen/SparseLU>
//#include <eigen3/Eigen/src/SparseLU>
using Eigen::MatrixXd;
using namespace Eigen;
#include <iostream>
#include <string>
#include "GeneralGeometry.h"
#include "Cell_parametres.h"
#include "KISS.h"
using namespace std;


//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS

class ISCorticalSl {
public:
    ISCorticalSl();
    ISCorticalSl( Vector3d axis_tmp , Vector3d plane_point_tmp , double layer_tmp );
    ISCorticalSl( Vector3d center_front_of_IS_tmp , Vector3d center_rear_of_IS_tmp , double radius_tmp  , double radius_inner_tmp);
    ISCorticalSl(const ISCorticalSl& other);
    ISCorticalSl& operator=( const ISCorticalSl &other );
    bool check_ISCorticalSl_caught( Vector3d position );
    double covered_distance_of_segment( Vector3d first_point , Vector3d second_point );
    Vector3d get_point_on_plane();
    Vector3d get_axis();

    Vector3d get_center_front_of_IS();
    Vector3d get_center_rear_of_IS();
    double get_radius_inner();	
    double get_radius();
    Vector3d project_point_on_surface( Vector3d position ); 

    std::vector<Vector3d> control_IS_Cortical_Sliding_points( std::vector<Vector3d> points_to_control , unsigned int& number_of_IS_points );
    std::vector<Vector3d> control_IS_Cortical_Sliding_points2( std::vector<Vector3d> points_to_control , unsigned int& number_in  , unsigned int& number_out );
    std::vector<Vector3d> control_IS_Cortical_Sliding_points3( std::vector<Vector3d> points_to_control , unsigned int& number_in  , unsigned int& number_out );


    std::vector<Vector3d> create_IS_Cortical_Sliding_points( unsigned int number_of_points );

    bool control_IS_Cortical_Sliding_points( Vector3d point );
    
    
    bool control_IS_Cortical_Sliding_one_point( Vector3d point );    
    std::vector<Vector3d> control_of_inside_and_margins_IS_1( std::vector<Vector3d> points_to_control , std::vector<Vector3d>& points_margins );    
    
    

	virtual ~ISCorticalSl();
private:
    Vector3d axis;
    Vector3d plane_point;
    double layer;
    double radius;
    double radius_inner;
    Vector3d center_front_of_IS;
    Vector3d center_rear_of_IS;	


};

#endif /* ISCORTICALSL_H_ */

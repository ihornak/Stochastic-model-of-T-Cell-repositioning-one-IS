//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS

#ifndef MTOC2_H_
#define MTOC2_H_


#include <iostream>
#include <eigen3/Eigen/Dense>
using Eigen::MatrixXd;
using namespace Eigen;
#include <vector>
#include <random>
#include <iostream>



#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/SparseLU>
using namespace std;

#include "simulationofCell.h"
#include "MTOCparam.h"
#include "GeneralGeometry.h"
#include "mersenne.h"
#include "simulationofCell.h"
//#include <eigen3/Eigen/src/SparseLU>





class MTOC2
{
public:
	MTOC2();
	MTOC2( unsigned int number_of_points_arg );
    MTOC2( unsigned int number_of_points_arg , Vector3d center_of_MTOC  );
    MTOC2( unsigned int number_of_points_arg , Vector3d center_of_MTOC , string plane );
    MTOC2( unsigned int number_of_points_arg , Vector3d center_of_MTOC , string plane  , string plane_2 );
    unsigned int get_number_of_points() const;

    void add_vector( Vector3d position );
    MTOC2& operator=( const MTOC2& tmp );
    MatrixXd get_coordinates() const;
    Vector3d get_point( unsigned int number ) const;
    Vector3d get_point_original( unsigned int number ) const;
    Vector3d get_original_orientation( unsigned int number ) const;
    void set_orientation( unsigned int index , Vector3d orientation  );
    Vector3d get_center() const;
    Vector3d get_center_original( ) const;
    MatrixXd get_original_coordinates() const;
    Vector3d get_center_of_gravity() const;
    double get_effectiveFriction() const { return this->effectiveFriction; }


    double get_radius() const;
    double set_radius(  );
    double set_polygon_angle(  );
    void set_friction();
    void set_coordinates( MatrixXd& coordinates_arg );
    void set_time_clock( double time_tmp ){ this->time_clock = time_tmp;  };

    Vector3d tangent_to_center( unsigned int number );
    Vector3d tangent_from_next_bead( unsigned int number );
    Vector3d tangent_between_beads(  unsigned int number_1 , unsigned int number_2 );
    void controlMTOC( );
    void adjustMTOC( );

    void getMatrixes_Sparse( MatrixXd &projection_Matrix );
    void getMatrixes_Sparse_2( MatrixXd &projection_Matrix );
    void getMatrixes_Sparse_3( MatrixXd &projection_Matrix );
    void getMatrixes_Sparse_4( MatrixXd &projection_Matrix );


    void print_norm_tangent_to_center();
    void print_points();
    Vector3d get_side_center( unsigned int side );
    unsigned int get_opposite_point_index( unsigned int point_index );
    unsigned int get_opposite_direct_point_index( unsigned int point_index );
    unsigned int get_random_point_on_opposite_side( unsigned int point_index );
    unsigned int get_random_point_on_opposite_side_without_bias( unsigned int point_index );





    MatrixXd rotate_MTOC_from_original_orientations();
    MatrixXd rotate_MTOC_from_original_orientations_2();
    MatrixXd set_coordinates_on_same_altitude();
    double magnitude_of_the_force( MatrixXd& Forces  );


    void Euler_algorithm( MatrixXd& Forces );
    void oneStepMidStepAlgorithm_1_half( MatrixXd& Forces );
    void oneStepMidStepAlgorithm_2_half( MatrixXd& Forces , MatrixXd& coordinates_arg);
    void oneStepMidStepAlgorithm( );
    void MidStepAlgorithm( double Time_0 , double Time_1 );

    double get_average_altitude();
    void set_axis_of_rotation();

    //original coordinates
    void resize_from_originals();


    void getRandomForces(  double timeStep , MatrixXd &randomForces );



	virtual ~MTOC2();

private:
	unsigned int number_of_points;
	MatrixXd coordinates;
	double radius_of_MTOC;
	double polygon_angle = 3.1415926535897 / 2.0 / 4.0;
  double effectiveFriction;
  MatrixXd original_orientations;
  Vector3d original_MTOC_orientation;
  std::vector<Vector3d> axis_of_rotation;
	MatrixXd original_coordinates;
	double time_clock;


};

#endif /* MTOC2_H_ */

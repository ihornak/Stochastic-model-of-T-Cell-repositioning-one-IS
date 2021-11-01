//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS

#ifndef SURFACE_H
#define SURFACE_H


#include "simulationofCell.h"
#include "Cell_parametres.h"
#include "KISS.h"
#include "IS_Dynein_Cell_surface.h"
#include "IS_Cortical_Sl_parameter.h"
#include "GeneralGeometry.h"

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include "ISCaptShrinkage.h"
#include "ISCorticalSl.h"



using Eigen::MatrixXd;
using namespace Eigen;
#include "Dynein.h"

using namespace std;


class Surface
{
public:
Surface();
Surface(  double density , double width );
Surface(  double density , double width , unsigned int capture_shrinkage );
Surface(  double density , double width , IS_Capt_Shrinkage IS_tmp );
Surface(  double density , double width , string cortical_Sliding , double density_in_IS );
Surface(  double density , double width , string cortical_Sliding , double density_in_IS , ISCorticalSl ISCorticalSl_arg );
Surface( string cortical_Sliding , ISCorticalSl ISCorticalSl_arg , unsigned int minus_integer );



Surface(const Surface& other);
Surface& operator=( const Surface &tmp );

double get_original_number_capt(){ return this->original_number_capt; };
double get_original_number_cort_1(){ return this->original_number_cort_1; };


double get_X_width(){ return this->x_width; };
double get_Y_width(){ return this->y_width; };
double get_Z_width(){ return this->z_width; };
unsigned int get_dynein_compartment_id_projected( Vector3d point_position );
unsigned int get_neccessary_dimension(){ return this->neccessary_dimension; };
std::map< unsigned int , std::vector<Vector3d> > get_map(){ return this->surface; };
bool dynein_surface_catch_control( Vector3d point_position , Vector3d& Catching_point  );
bool dynein_surface_catch_control_tangent( Vector3d point_position , Vector3d a_point , Vector3d tangent_b , Vector3d& Catching_point );
void erase_dynein_point( unsigned int map_key , unsigned int array_index );
void erase_vector_points_with_key( unsigned int map_key );


std::vector<Vector3d> get_dynein_points( unsigned int key );
void set_dynein_points( unsigned int key , std::vector<Vector3d> vectors_arg );
void add_dynein_point( unsigned int key , Vector3d vector_arg );

unsigned int get_dynein_motors_number();
std::vector<Vector3d> create_capture_shrinkage_points(  unsigned int number_of_points  , IS_Capt_Shrinkage IS_tmp );
void add_dynein_points_to_compartment_with_number( std::vector<Vector3d> vectors , unsigned int compartment_id );

Vector3d project_point_on_surface( Vector3d position );
void project_and_add_points_to_surface( std::vector< Vector3d > points );
std::vector<Vector3d> create_cortical_sliding_points(  unsigned int number_of_points  , ISCorticalSl IS_tmp );

std::vector<Vector3d> create_Cortical_Sliding_points_outside_IS( unsigned int number_of_points , ISCorticalSl ISCorticalSl_arg );
std::vector<Vector3d> get_all_dynein_point();
void erase_map( );
void change_part_of_IS_capture_points(   IS_Capt_Shrinkage IS_tmp);
void change_part_of_IS_points( ISCorticalSl IS_tmp );
void change_part_of_points_surface_and_IS( ISCorticalSl IS_tmp );
void control_points_outside_IS(  ISCorticalSl IS_tmp );
double triangular_distribution( double upper_boundary );
void print_inner_map();

~Surface();

private:
    double A_axis;
    double B_axis;

    double x_width;
    double y_width;
    double z_width;
    const unsigned int neccessary_dimension = 1000;
    unsigned int original_number_capt;
    unsigned int original_number_cort_1;   

    double density;
    std::map< unsigned int , std::vector<Vector3d> > surface;

};

#endif // SURFACE_H

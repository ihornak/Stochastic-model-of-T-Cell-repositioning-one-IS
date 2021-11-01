

#ifndef CELL_H_
#define CELL_H_

#include "simulationofCell.h"
#include "Microtubule.h"
#include "Nucleus.h"
#include "MTOC2.h"
#include <errno.h>
using Eigen::MatrixXd;

#include <omp.h>
#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include "ISCaptShrinkage.h"
#include "IS_Capture_shrinkage_parameters.h"
#include "ISCorticalSl.h"
//#include "ISCorticalSliding.h"
#include "MTOCparam.h"
#include "IS_Plane_force.h"
#include "CellShape.h"
#include "Cell_parametres.h"
#include "IS_Cortical_Sl_parameter.h"
#include "IS_Dynein_Cell_surface.h"
#include "GeneralGeometry.h"
#include "Surface.h"
#include "mtoc_simple.h"
#include "ideal_mtoc.h"
#include <fstream>               //mandatory for ofstreaming
#include <sstream>
#include <random>
#include <iostream>

using namespace std;

//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning

class Cell
{
public:
	//This is constructor for unconstrained cytoskeleton in 3D space
	Cell(  );
    	Cell( unsigned int number_Of_Microtubules , unsigned int number_Of_extra_Microtubules , unsigned int servant );
    	Cell( unsigned int number_Of_Microtubules , unsigned int number_Of_extra_Microtubules , unsigned int servant  , unsigned int servant_2 );



    	Microtubule getMicrotubule( unsigned int index );
    	unsigned int get_microtubule_number();
    	void resize_microtubules( );
    	unsigned int get_number_of_microtubule_with_dynein_index( unsigned int dynein_index );
    	IS_Capt_Shrinkage get_IS_capture_shrinkage();
    	double get_cytoskeleton_friction();
    	MTOC2 get_MTOC(   ) const { return this->MTOC; };
    	Vector3d get_MTOC_center(  ) const { return this->MTOC.get_center(); }; 
    	Vector3d get_IS_center(  );    
    	double get_MTOC_IS_distance();

    	void print_Cell( unsigned int index );
    	void print_Cell_parametres( );

    	ISCorticalSl get_IS_cortical_sliding_first(){ return this->IS_Cortic_Sl; };

    	//Microtubule wall interaction
        Vector3d force_position_cell_wall( Vector3d position );
        Vector3d force_position_cell_wall_simple_elipsoid( Vector3d position  );
        Vector3d force_position_cell_wall_two_elipsoid( Vector3d position  );
        Vector3d force_position_cell_wall_simple_elipsoid_MTOC( Vector3d position  );
        void force_whole_microtubule_cell_wall( MatrixXd &force , unsigned int index );
        void MTOC_microtubule_force( Vector3d& center_bead_force , Vector3d& point_bead_force , unsigned int microtubule_number );
        void MTOC_microtubule_two_points_force( Vector3d& first_point_force , Vector3d& second_point_force , unsigned int microtubule_number );
        void MTOC_microtubule_two_points_force_with_bending( Vector3d& first_point_force , Vector3d& second_point_force , Vector3d& bending_opposite_MTOC_point , Vector3d& bending_MTOC_point , unsigned int microtubule_number );
        Vector3d MTOC_microtubule_two_points(  unsigned int microtubule_number );
        Vector3d force_position_cell_wall_sphere( Vector3d position  );
        //MICRO NUCLEUS INTERACTION
        void force_whole_microtubule_nucleus( MatrixXd &force , unsigned int index );
        //MTOC WALL INTERACTION
        Vector3d MTOC_point_wall_interaction( Vector3d position );
        void force_MTOC_cell_wall( MatrixXd &force );
        //MTOC NUCLEUS INTERACTION
        void MTOC_nucleus_interaction( MatrixXd &force_on_MTOC );
        //MTOC IS INTERACTIONideal_MTOC();
        Vector3d MTOC_dynein_capture_shrinkage_force( Vector3d MTOC_bead ,  Vector3d Microtubule_catching  );
        Vector3d MTOC_dynein_capture_shrinkage_force( unsigned int microtubule );
        Vector3d MTOC_dynein_cortical_sliding_force( Vector3d MTOC_bead ,  Vector3d Microtubule_catching  );
        void MTOC_naive_circular_force(  MatrixXd &force );
        void MTOC_naive_force_to_center(  MatrixXd &force_MTOC );
        //MICRO IS_CAPTUTE_SHRINKAGE INTERACTION
        void IS_micro_control_one_micro_capture_shrinkage_shrinkage( unsigned int microtubule );
        //MICRO IS CAPTURE SHRINKAGE INTERACTION
        void IS_capture_shrinkage_shrinkage_of_microtubule();
        void IS_micro_control_one_micro_cortical_sliding( unsigned int microtubule );
        bool check_IS_cortical_sliding_caught( Vector3d position );
        bool check_IS_cortical_sliding_caught_Main( Vector3d position );
        void check_IS_cortical_sl_catching_points();

   	//MICRO IS_SL INTERACTION
        void IS_cortical_sl_all_micro_control_CAUGHT( );
        bool check_IS_cortical_sl( Vector3d position );
        bool IS_cortical_sl_one_micro_control_CAUGHT( unsigned int microtubule );
        void force_IS_cortical_sl( MatrixXd& force_micro , unsigned int microtubule );
        void force_IS_cortical_sl_basic( MatrixXd& force_micro , unsigned int microtubule );
        void force_IS_cortical_sl_basic_with_covered_micro( MatrixXd& force_micro , unsigned int microtubule );
        void force_IS_cortical_sl_basic_with_covered_micro_closest_point( MatrixXd& force_micro , unsigned int microtubule );
        void force_IS_cortical_sl_all_points( MatrixXd& force_micro , unsigned int microtubule );
        void construct_catching_points_4( unsigned int microtubule );
        void construct_catching_points_8( unsigned int microtubule );
        double get_covered_distance_by_IS( unsigned int microtubule );
        Vector3d dynein_force_one_bead(Vector3d position , Vector3d tangent_to_be_projected , double absolute_value_of_force );
        void whole_microtubule_dynein_force(unsigned int microtubule , MatrixXd& force_whole_microtubule );
        void whole_microtubule_dynein_force_distributed_micro_7(unsigned int microtubule , MatrixXd& force_whole_microtubule );
        void microtubule_dynein(unsigned int microtubule , MatrixXd& force_whole_microtubule );
        void get_Projection_Matrixes( unsigned int microtubule , MatrixXd &projection_matrix );
        void project_microtubules_7_on_elipsoid( );
        void check_IS_sliding_all_points_of_microtubule_3( unsigned int microtubule );
        void IS_sl_check_catching_points_all_micro_2( );


   	//MICRO DYNEIN CELL SURFACE INTERACTION
        void Dynein_on_surface_force_2(unsigned int microtubule , MatrixXd& force_whole_microtubule );
        void whole_microtubule_Surface_Dynein_force(unsigned int microtubule , MatrixXd& force_whole_microtubule );
        void set_density_surface_dynein( Surface& density );
        void set_density_IS_capture_shrinkage( Surface& density ){ this->Capture_Shrinkage_dynein = density; };


   	//MICRO DYNEIN CELL SURFACE INTERACTION std::pair < Vector3d ,double >
        void all_microtubules_random_dynein_on_surface_catch_pair_abscissa();
        void microtubule_random_dynein_on_surface_catch_pair_abscissa( unsigned int microtubule );
        std::vector<Vector3d> microtubule_point_attach_dynein( Vector3d point_position );
        void stepping_and_detachment_of_all_microtubule();
        void stepping_and_detachment_of_all_microtubule_projection();

        void all_microtubules_random_dynein_on_surface_catch_pair_abscissa_2();
        void microtubule_random_dynein_on_surface_catch_pair_abscissa_2( unsigned int microtubule );
        void microtubule_random_dynein_on_surface_catch_pair_abscissa_plane( unsigned int microtubule );
        void microtubule_random_dynein_on_surface_catch_pair_abscissa_space( unsigned int microtubule );
        void microtubule_random_dynein_on_surface_catch_pair_abscissa_space_line( unsigned int microtubule );
        void microtubule_random_dynein_on_surface_catch_pair_abscissa_aproximation_space_line( unsigned int microtubule );
        void microtubule_random_dynein_on_surface_catch_pair_abscissa_aproximation_space_line_2( unsigned int microtubule );
        void microtubule_random_dynein_on_surface_catch_pair_abscissa_aproximation_true( unsigned int microtubule );
        void attach_dynein_on_microtubule( Vector3d point , Vector3d tangent , unsigned int key , unsigned int microtubule );
        void catch_pair_abscissa();
        void microtubule_catch_pair_abscissa( unsigned int microtubule );

        std::vector<Vector3d> get_dynein_in_compartment( unsigned int ID_map_index );
        unsigned int get_number_of_dynein_motors();
	unsigned int get_number_of_dynein_motors_in_micro_capture();
        unsigned int get_number_of_dynein_motors_IS();

   	//MICRO REAL DYNEIN CELL SURFACE INTERACTION
        void catch_pair_abscissa_real_dynein();
        void microtubule_catch_pair_abscissa_real_dynein( unsigned int microtubule );
        void stepping_and_detachment_of_all_microtubule_projection_real_dynein();
        void detach_microtubule_NINE_and_project_points_on_surface( unsigned int microtubule );

	////////////////////////////////////////////////////////////////////////////////////
        //MICRO REAL DYNEIN CELL IS INTERACTION
        void check_caught_micro_IS_with_real_dynein( unsigned int microtubule );
        void check_caught_micro_IS_with_real_dynein_all_micro( );

        void microtubule_catch_pair_abscissa_real_dynein_in_IS( unsigned int microtubule );
        void microtubule_catch_pair_abscissa_real_dynein_in_IS_all_micro( );
        void control_length_of_micro_IS();
	void control_length_of_micro_IS_2();
	void control_microtubule_20_three_beads_dyneins();


       //general utilities
       Vector3d project_point_on_surface( Vector3d position );
       unsigned int get_dynein_compartment_id( Vector3d position );
       unsigned int get_dynein_compartment_id_projected( Vector3d position );
       void create_point_in_IS_cortical_sliding(  );
       void project_and_add_points_to_surface( std::vector< Vector3d > points );
       void add_points_to_surface( std::vector< Vector3d > points );
       Vector3d force_MTOC_bending_first_point( unsigned int microtubule );
       Vector3d force_MTOC_bending_first_point_2( unsigned int microtubule , Vector3d& force_2 );


       //bending mtoc force
       unsigned int get_opposite_micro_number( unsigned int microtubule );
       void bending_micro_force( unsigned int micro_id ,  Vector3d& first_micro , Vector3d& second_micro , Vector3d& mtoc_center );

       //abstract mtoc
       void bending_abstract_mtoc( unsigned int micro_id ,  Vector3d& first_micro , Vector3d& second_micro , Vector3d& mtoc_center );
       void resize_micro_MTOC(  );
       void resize_micro_with_different_first_segment_and_MTOC(  );
       void ultimate_resizing();

        //PRINT
        void BIG_PRINT( unsigned int time_Step );
        void BIG_PRINT_numerical( unsigned int time_Step , unsigned int key );
        void BIG_PRINT_numerical_micro( unsigned int time_Step , unsigned int key );
        void print_time( unsigned int time_Step );
        void print_center_MTOC();


        void print_center_MTOC_numerical( double time_Step , unsigned int key );
	void print_dynein_capture_shrinkage_numerical( double time_Step , unsigned int key );
        void print_dynein_cortical_sliding_numerical( double time_Step , unsigned int key );
	void print_microtubules_numerical( unsigned int time_Step , unsigned int key );
 	void print_dynein_numerical( double time_Step , unsigned int key );
	void print_dynein_numerical_unattached( double time_Step , unsigned int key );


        unsigned int get_number_of_dynein_motors_cortical_sliding_micro();
	void get_attached_dynein_motors_cortical_sliding_micro_both_surface_IS( double& concentration_in_IS ,  double& concentration_outside_IS );

	void print_center_numerical_parameters( unsigned int key );
        void print_center_MTOC_radius_numerical( unsigned int key );
        void print_center_MTOC_azimuth_numerical( unsigned int key );
        void print_center_MTOC_altitude_numerical( unsigned int key );
        void print_center_MTOC_z_numerical( unsigned int key );

        Vector3d force_simple_one_micro( unsigned int microtubule );
        Vector3d simple_dynein_capture_shrinkage_force( unsigned int microtubule );
        Vector3d force_nucleus( );
        Vector3d force_wall( );
        void MidStep_3();
        void Euler_algorithm( );

       unsigned int get_number_of_dynein_motors_capture_shrinkage_1();
       void timeDevelopment( double Time_0 , double Time_1 , unsigned int number_of_Pictures );
       void timeDevelopment( double Time_0 , double Time_1 );
       void timeDevelopment_2( double Time_0 , double Time_1 );
       void timeDevelopment_2( double Time_0 , double Time_1 , unsigned int number_of_Pictures );
       void Euler_timeDevelopment( double Time_0 , double Time_1 );
       void Euler_timeDevelopment( double Time_0 , double Time_1 , unsigned int number_of_Pictures );
       void timeDevelopment_numerical_results( double Time_0 , double Time_1, unsigned int key );
       void timeDevelopment_numerical_results_2( double Time_0 , double Time_1, unsigned int key );
       void test_dynein( string param );

	virtual ~Cell();
private:
	double a_axis;
	double b_axis;
	Nucleus nucleus;
        MTOC2 MTOC;

        unsigned int number_of_microtubules;
        unsigned int number_of_microtubules_extra;
        Microtubule *array_Of_Microtubules = NULL;
        IS_Capt_Shrinkage IS_Capture_Shrinkage;
        ISCorticalSl IS_Cortic_Sl;
        Surface density_surface_dynein;
        Surface Capture_Shrinkage_dynein;

        ideal_MTOC abstract_MTOC;
        MTOC_simple mtoc;
        std::random_device rd;
        unsigned int cell_id;


};

#endif /* UNCONSTRAINEDCELL_H_ */

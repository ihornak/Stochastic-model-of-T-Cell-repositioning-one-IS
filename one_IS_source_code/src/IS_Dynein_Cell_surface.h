#ifndef IS_DYNEIN_CELL_SURFACE_H_
#define IS_DYNEIN_CELL_SURFACE_H_

//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS
namespace IS_Dynein_Cell_surface
{
	const unsigned int number_of_polygon_lower = 0;
    	const unsigned int number_of_polygon_higher = 12;
	const unsigned int number_of_mito = 5;    
	const double force_per_lenght = 0.18e-6;//N/microM 0.16   0.12e-6
	const double force = 5e-12;//N/microM 0.45  10e-12
	const unsigned int lower_side = 2;
    	const double close_boundary = 0.3e-6;
    	const double dynein_on_surface_catch_radius = 12.5e-9;
    	const double treshold_cut_tangent = 10.0e-7;
    	const double new_radius = 0.075e-7;
    
}


#endif

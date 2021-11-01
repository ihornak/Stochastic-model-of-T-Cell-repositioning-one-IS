//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS

#ifndef CELL_PARAMETRES_H_
#define CELL_PARAMETRES_H_
#include <iostream>
#include <stdio.h>      /* printf */
#include <math.h>
#include "IS_Capture_shrinkage_parameters.h"
using namespace std;



namespace Cell_parametres
{
    	extern double time;
    
	const double wall_cell_k1 = 6.0e-5;
	const double wall_cell_k2 = 1.0;
    	const double wall_cell_k3 = 1e-14;
	const double A_AXIS = 5.0e-6;
	const double B_AXIS = 5.0e-6;
    	const double IS_area_procentage = 0.2;
	const double B_AXIS_UPPER = 5.3e-6;
    	const double B_AXIS_Lower = 5.0e-6;
    	const double projection_radius = 6.9e-6;
	const double width_of_dynein = 0.1e-6;
    	const double width_of_dynein_2 = 1.5e-6;
	extern double force_lower_ellipsoid;
	extern double number_micro_capture_shrinkage;
	extern double number_micro_cortical_sliding;
	extern double procentage_of_the_force;
	const double cylinder_width = 1.0e-6;
	const double cylinder_height_realistic = 1.5e-6;
	extern double z_cylinder;
	const double B_AXIS_below = 3.0e-6;
	const double cylinder_height = 3.0e-6;
	const double lower_ellipsoid_B_AXIS = 3.0e-6;
    	const double surface_density = 10;    
	double calculate_z_cylinder(  );
	void calculate_repulsive_force(  );

}


#endif /* CELL_PARAMETRES_H_ */

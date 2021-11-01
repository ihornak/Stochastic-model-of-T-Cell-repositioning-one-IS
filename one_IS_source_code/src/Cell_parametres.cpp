//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning
#include "Cell_parametres.h"
#include "simulationofCell.h"


namespace Cell_parametres {
	double z_cylinder = 0;
    double time = 0;
	double number_micro_capture_shrinkage;
	double number_micro_cortical_sliding;
	double procentage_of_the_force = 1.0;


	double force_lower_ellipsoid = 0.0;
	double calculate_z_cylinder(  )
	{
		if( cylinder_width >=  A_AXIS )
		{
			cout<<" Cell_parametres::cylinder_width >=  Cell_parametres::A_Axis in  calculate_z_cylinder(  )"	<<endl;
			cout<<" ERROR_ID = 687494596 "<<endl;
			throw("");
		}
		double tmp = 1.0 - ( cylinder_width / A_AXIS ) * ( cylinder_width / A_AXIS );
		double z_value = B_AXIS * sqrt( tmp );
		return z_value;
	}
	void calculate_repulsive_force(  )
	{
		double dynein_force = number_micro_capture_shrinkage * IS_Capture_shrinkage_param::dynein_Force_capture_shrinkage;
		dynein_force = dynein_force + number_micro_cortical_sliding * sim_of_Cell::dynein_Force_cortical_sliding;
		force_lower_ellipsoid = dynein_force * procentage_of_the_force;
	}




}

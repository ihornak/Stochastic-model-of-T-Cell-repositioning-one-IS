#ifndef SIMULATIONOFCELL_H_
#define SIMULATIONOFCELL_H_



namespace sim_of_Cell
{
        //Parameters of the simulations
        const unsigned int Microtubules = 100;
        const double Time_1 = 40;
        const double density_of_dynein_in_IS_Capture_Shrinkage = 300;
        const double density_of_dynein_motor_Cortical_Sliding = 300;
        const double density_of_dynein_motor_surface = 0;
        const unsigned int number_of_samples = 4;

	//Random forces switched of
        const bool random_force_switch = false;
        const bool random_force_switch_MTOC = false;
	//Physical constants
	const double resting_distance = 0.8e-6; 	
	const double Temperature = 300.0;									
	const double k_bending_analytical =  2.3e-23;	
	const double k_bending =  2.3e-23 / resting_distance;
	const double viscosity_water = 0.001002;
	const double viscosity = 0.09;
        const double boltzmann_constant = 1.3806503e-23;	
	const double deeper_viscosity = viscosity * 10.0;
	const double surface_width = 0.5e-6;
	
	
	//Important Model constants
	const unsigned int MicrotubulePoints = 20;	
	const unsigned int REDUCTION = 5;            			
        const double treshold_distance = 1e-6;
        const bool big_print_switch = false;
        const double time_Step = 0.7e-3; 
        const double time_Step_half = time_Step / 2.0;
        const double dynein_Force_cortical_sliding = 5e-12;
        const double PI = 3.1415926535897;
        const double angle = 0.0;
        const double multiply_friction_constant = 1;
        const bool paralelization_switch = 1;
        const unsigned int general_switch = 1;
        const bool cytoskeleton_asymetry_switch = true;
        const double surface_replacement_procentage = 0.1;

}


#endif /* SIMULATIONOFCELL_H_ */

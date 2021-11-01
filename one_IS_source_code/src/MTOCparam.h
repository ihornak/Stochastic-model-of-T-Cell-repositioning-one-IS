//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS
#ifndef MTOCPARAM_H_
#define MTOCPARAM_H_


namespace MTOCparam
{
	const double MTOC_MTOC_spring = 1e-11;  // interaction between polygon of MTOC
	const double polygon_side = 1.5e-7;
	const double distance = 3e-7;
	const double polygon_distance = 1e-7;
	const double distance_mito_MTOC = 3e-7;
	const double distance_micro_MTOC = 3e-7;
    	const double elastic_kappa_MTOC_microtubule = 1e-11;
	const double kappa_MTOC =  8e-18 ;//8e-18
	const double inter_final_MTOC_spring = 0.0;
	const double inter_MTOC_diagonal_spring = 1e-11;
	const double opposite_sides_resting_distance = 8e-7;
	const double elastic_kappa_MTOC_opposidde_sides = 2e-11;
	const double elastic_kappa_disssipatice_control = 1e-10;
    	const unsigned int boundary_MTOC_points = 20;
    	const unsigned int micro_in_polygon = 5;
    
    	const double MTOC2_center_micro_kappa = 0e-6;//
    	const double MTOC2_point_micro_kappa = 30e-6;//3
    	const double MTOC2_circle_param = 1.0e7;//
    	const double MTOC_radius = 0.5e-6;
    	const double MTOC_center_index = 2; 
    	const double MTOC_angle_rotation = 2.95;
    	const unsigned int number_of_opposite_side = 2;
}



#endif /* MTOCPARAM_H_ */

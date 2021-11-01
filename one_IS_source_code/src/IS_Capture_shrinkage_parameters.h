//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS

#ifndef IS_CAPTURE_SHRINKAGE_PARAMETERS_H_
#define IS_CAPTURE_SHRINKAGE_PARAMETERS_H_
#include "simulationofCell.h"

namespace IS_Capture_shrinkage_param
{
    const unsigned int number_of_polygon_lower = 0;
    const unsigned int number_of_polygon_upper = 1;
    const unsigned int number_of_mito = 1;
    const double cut_distance = 0.5e-6;
    const double dynein_Force_capture_shrinkage = 4e-12;
    const double cutting_distance_depol_micro = 0.3e-6; 
    const double depolimerization_distance_threshold = 2e-6;
    const double procentage_constant = 1.3;

    const double minimal_radius =  4.9999e-6;


    const double z_coordinate_front = 3.7e-6;
    const double z_coordinate_back = 6.3e-6;
    const double radius = 0.4e-6;
    const double azimutal_angle = 0.0 *  sim_of_Cell::PI;
    const double polar_angle = 0.75 *   sim_of_Cell::PI;   //



}






#endif 

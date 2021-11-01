

#ifndef DYNEIN_H
#define DYNEIN_H
#include <iostream>
#include <stdio.h>     
#include <math.h>
#include "simulationofCell.h"

using namespace std;

//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS
namespace Dynein
{
    const double k_d = 1.0; //s^{ -1 }
    const double F_D = 3e-12;
    
    double calculate_detach_probability( double Force );
    double calculate_detach_probability_per_time_step( double Force );
        
    
    const double F_stall = 6e-12;    
    double backward_stepping_probability_between_0_Stall_force( double Force );
    double backward_stepp_prob_smaller_than_Stall_f_per_step( double Force );//
        
    
    const double v_F = 1000.0e-9; //nm s^{ -1 }
    const double d = 8e-9;//m    
    const double forward_stepping_probability = v_F / d;
    const double forward_stepping_probability_per_step = v_F / d / ( 1.0 / sim_of_Cell::time_Step );
    
    const double v_b = 6e-9;//m * s^{-1}
    const double backward_prob_force_bigger_then_Stall_force = v_b / d;
    const double backward_prob_force_bigger_then_Stall_force_per_step = backward_prob_force_bigger_then_Stall_force / ( 1.0 / sim_of_Cell::time_Step );
    
    
    const double L_0 = 110e-9; //m
    const double tmp_cutting_distance = 200e-9;
    const double K_a = 5.0;
    const double attach_rate_per_step = 5.0 * sim_of_Cell::time_Step;
    const double multiplicator = 1.0;
    const double reduction_abscissa = 2e-9;
    const double L_0_perpen = 10e-9;
    
    const double alpha = 1e-4;
    
    double detachment_probability( double Force );
    double detachment_probability_per_step( double Force );
    
    const unsigned int number_of_segment_steps = 7;
    
}

#endif // DYNEIN_H

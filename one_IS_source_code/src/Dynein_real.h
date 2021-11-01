
#ifndef DYNEIN_REAL_H
#define DYNEIN_REAL_H
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "simulationofCell.h"

using namespace std;


//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS


namespace Dynein_real
{
    const double k_d = 1.0; 
    const double F_D = 2e-12;
    const double F_stall = 4e-12;
    const double step = 8e-9;
    const double v_F = 1000.0e-9; 
    const double forward_stepping_probability = v_F / step;
    const double forward_stepping_probability_per_step = v_F / step * sim_of_Cell::time_Step;
    const double v_b = 6e-9;
    const double backward_prob_force_bigger_then_Stall_force = v_b / step;
    const double backward_prob_force_bigger_then_Stall_force_per_step = backward_prob_force_bigger_then_Stall_force * sim_of_Cell::time_Step;

    const double K_a = 5.0;
    const double attach_rate_per_step = K_a * sim_of_Cell::time_Step;
    const double alpha = 4e-4;
    const double L_0 = 18e-9;

    double calculate_attachment_probability( double distance );
    double calculate_attachment_probability_per_time_step( double distance );

    double calculate_detach_probability( double Force );
    double calculate_detach_probability_per_time_step( double Force );



    double backward_stepping_probability_between_0_Stall_force( double Force );
    double backward_stepp_prob_smaller_than_Stall_f_per_step( double Force );//

    double detachment_probability( double Force );
    double detachment_probability_per_step( double Force );

}

#endif // DYNEIN_H

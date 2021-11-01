

#include "Dynein.h"
//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS
namespace Dynein
{
    //double k_d = 5;
    
    double calculate_detach_probability( double Force )
    {
        return k_d * exp ( fabs( Force ) / F_D );
    }
    
    double calculate_detach_probability_per_time_step( double Force )
    {
        double probability_per_second = calculate_detach_probability( Force );
        double constant_per_time_step = ( 1.0 / sim_of_Cell::time_Step );
        return probability_per_second / constant_per_time_step;
    }
    
    double backward_stepping_probability_between_0_Stall_force( double Force )
    {
        return v_F / d * ( 1.0 - Force / F_stall );        
    }
    
    double backward_stepp_prob_smaller_than_Stall_f_per_step( double Force )
    {
        double probability_per_second = backward_stepping_probability_between_0_Stall_force( Force );
        double constant_per_time_step = ( 1.0 / sim_of_Cell::time_Step );
        return probability_per_second / constant_per_time_step;
    }
    

    double detachment_probability( double Force )
    {
        return k_d * exp( fabs( Force ) / F_D );
        //return 0.0;
    }
    
    double detachment_probability_per_step( double Force )
    {
        //cout<<"Force = "<<Force<<endl;
        double prob = detachment_probability( Force ) / ( 1.0 / sim_of_Cell::time_Step );
        prob = 1;
        //cout<<"prob = "<<prob<<endl;
        return prob;    
    }
    
    
    
}

#include <iostream>
#include "src/MTOC2.h"
#include "src/Microtubule.h"
#include "src/Cell.h"
#include "ISCorticalSl.h"
#include "generalUtilities.h"
#include "Surface.h"
#include "Dynein.h"
#include "numerical_result.h"
#include "IS_Cortical_Sl_parameter.h"
#include "Generator.h"




int main(int argc, char **argv )
{
    cout.precision(17); 
 
    //This is the initialization of the random number generator located in mersenne.h/cpp files
    //We use Mersenne Twister and the generator is initialized with time
    //Every thread has its onwn generator, we use 28 threads(increase if the number of threads > 28)
    mersenne_twisters_64::initialize_generators( 28 );
    //Time measurment
    time_t start,end;
    time ( &start );

     //////////////////////////////////////////////////////////////////////////////////////////////// 
    //This part produces the files for the video depicting the process
    //The number of microtubules
    unsigned int number_Of_Microtubules = sim_of_Cell::Microtubules;
    //The cytoskeleton is symetrical - no extra microtubules are created at one side of the MTOC
    unsigned int number_Of_extra_Microtubules = 0;
    //Times of simulations
    double Time_0 = 0.0;
    double Time_1 = sim_of_Cell::Time_1;
    string cortical_Sliding = "cortical_Sliding";
    unsigned int a = 0;
    double surface_density = sim_of_Cell::density_of_dynein_motor_surface;
    double corti_sl_density = sim_of_Cell::density_of_dynein_motor_Cortical_Sliding;


    //Constructor of the cell
    //Creation of the cell with two IS 
    //The cytosleton is not in the configuration of the lowest energy.
    Cell cell( number_Of_Microtubules , number_Of_extra_Microtubules , a );
    //Cortical sliding IS of the Cell
    ISCorticalSl tmp_1 = cell.get_IS_cortical_sliding_first();

    //Creation of the dynein in the cortical sliding IS
    Surface first_Surface(  surface_density ,  sim_of_Cell::surface_width , cortical_Sliding , corti_sl_density , tmp_1 );
    //Setting of the dyneins to the IS    
    cell.set_density_surface_dynein( first_Surface );
    unsigned int capture_shrinkage_density = sim_of_Cell::density_of_dynein_in_IS_Capture_Shrinkage;
    IS_Capt_Shrinkage tmp;
    //Capture-shrinkage IS of the cell
    tmp = cell.get_IS_capture_shrinkage();
    //Creation of the capture-shrinkage dyneins
    Surface Surface_capture_shrinkage( capture_shrinkage_density, 1e-6 , tmp );
    //Passing the dyneins to the IS
    cell.set_density_IS_capture_shrinkage( Surface_capture_shrinkage );

    //Relaxation of the cytoskeleton
    cell.timeDevelopment_2( Time_0 , 2 );
    //The number of pictures for the video of the repositioning    
    unsigned int number_of_pictures = 10;
    //One simulation run    
    cell.timeDevelopment_2( Time_0 , Time_1 , number_of_pictures );

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    //This part performs parallel computation of multiple simulation runs and saves numerical results
    
    
    /*
    numerical_result::test_parallel( );

    //measurement of the time and printing of the duration of the simulation
    
    */
    time(&end);
    double dif = difftime( end , start );
    printf ("Elasped time is %10.6lf seconds.\n", dif );


    return 0;

}



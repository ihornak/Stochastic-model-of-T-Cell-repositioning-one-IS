#include "numerical_result.h"

namespace numerical_result 
{

    void repolarization_sample( unsigned int index )
    {
            //Index identifies the simulation run
    	    //Times of simulation
            double Time_0 = 0.0;
            double Time_1 =  sim_of_Cell::Time_1;

            Cell bunka( sim_of_Cell::Microtubules , 0 , index  ); 
            
            //Setting of cortical sliding dyneins
	    ISCorticalSl tmp_1 = bunka.get_IS_cortical_sliding_first();             
            Surface first_Surface(  sim_of_Cell::density_of_dynein_motor_surface , sim_of_Cell::surface_width , "cortical_Sliding" , sim_of_Cell::density_of_dynein_motor_Cortical_Sliding , tmp_1 );
            bunka.set_density_surface_dynein( first_Surface );
            
            //Setting of capture-shrinkage dyneins
    	    IS_Capt_Shrinkage tmp;	
            tmp = bunka.get_IS_capture_shrinkage();
            Surface Surface_capture_shrinkage( sim_of_Cell::density_of_dynein_in_IS_Capture_Shrinkage, 1e-6 , tmp );
            bunka.set_density_IS_capture_shrinkage( Surface_capture_shrinkage );
      
            bunka.timeDevelopment_2( Time_0 , 3 );
            bunka.timeDevelopment_numerical_results_2( Time_0 , Time_1, index ); 
    }	

   void test_parallel( )
   {
	//Generation of numbers identifying simulation runs
    	std::uniform_int_distribution<> distribution{ 0 , 10000 };
    	unsigned int number_of_generator = omp_get_thread_num();
	unsigned int lower_boundary = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] ); 
	lower_boundary = lower_boundary * 100;
	
	//Prints all parameters of the simulation
	numerical_result::create_document_file();
	int j;
	//Parallel computing of multiple simulation runs
	#pragma omp parallel for
	for ( int j = lower_boundary ; j < lower_boundary + sim_of_Cell::number_of_samples ; ++j )
	{
		int threadNum;
		threadNum = omp_get_thread_num();
		//One simulation run
		repolarization_sample( j );
	}
   }



    
//Creates the file with the parameters of the simulation
void  create_document_file()
{
      		std::string name_of_text_file = "./picturesVideos/numerical_results/dokument.txt";	
      		ofstream fout;
      		fout.open( name_of_text_file.c_str() , std::fstream::app );	
      		fout<<sim_of_Cell::Time_1<<endl;  //write to file 
		fout<<sim_of_Cell::Microtubules<<endl;
		fout<<sim_of_Cell::viscosity<<endl;
		fout<<sim_of_Cell::density_of_dynein_motor_Cortical_Sliding<<endl;
		fout<<sim_of_Cell::density_of_dynein_motor_surface<<endl;
		fout<<sim_of_Cell::density_of_dynein_in_IS_Capture_Shrinkage<<endl;
		fout<<Cell_parametres::A_AXIS<<endl;
		fout<<Nucleus_parametres::A_AXIS<<endl;
		fout<<Nucleus_parametres::z_coordinate<<endl;
		fout<<IS_Capture_shrinkage_param::polar_angle/sim_of_Cell::PI<<endl;		
		fout<<IS_Capture_shrinkage_param::polar_angle<<endl;
		fout<<sim_of_Cell::MicrotubulePoints<<endl;
		fout<<sim_of_Cell::resting_distance<<endl;
		fout<<sim_of_Cell::REDUCTION<<endl;
		fout<<sim_of_Cell::number_of_samples<<endl;
                fout<<MTOCparam::MTOC_radius<<endl;
                fout<<IS_Cortical_Sl_parameter::radius<<endl;		
      		fout.close();
}	


    
}




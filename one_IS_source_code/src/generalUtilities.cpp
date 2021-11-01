//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS
#include "generalUtilities.h"


void get_program_ID( unsigned int number_of_micro , unsigned int number_of_pictures , double TIME_1 , double A_axis , double B_axis , double plane )
{
	//cout<<"get_program_ID"<<endl;
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );

    unsigned int year = ( unsigned int ) (now->tm_year + 1900);
    unsigned int month = ( unsigned int ) (now->tm_mon + 1);
    unsigned int day = ( unsigned int ) (now->tm_mday);
	FILE *out;
	out = fopen( "program_ID.txt" , "w");
	fprintf(out," year = %d , month = %d , day = %d\n", year , month , day );
	fprintf(out," number of microtubules = %d , number of pictures = %d , time = %5.2f\n", number_of_micro ,  number_of_pictures , TIME_1 );
	fprintf(out," axisA = %d , number of pictures = %d , time = %5.2f\n", number_of_micro ,  number_of_pictures , TIME_1 );
	fprintf(out," Cell_parametres::viscosity = %8.6f Cell_parametres::capture_shrinkage_force = %20.18f \n", sim_of_Cell::viscosity , IS_Capture_shrinkage_param::dynein_Force_capture_shrinkage );
	fprintf(out," mito = %d , polygon = %d\n", IS_Capture_shrinkage_param::number_of_mito , IS_Capture_shrinkage_param::number_of_polygon_upper  );
	fprintf(out,"sim_of_Cell::multiply_friction_constant = %5.2f\n", sim_of_Cell::multiply_friction_constant );
	fclose( out );
}



double calculateEffectiveFriction_Mitochondria( double lenght , double radius )
{
	//it has cylindrical shape
	double Volume = ( sim_of_Cell::PI * radius * radius ) * lenght;
	double radius_sphere_same_volume = std::pow( ( Volume * 3.0 / 4.0 / sim_of_Cell::PI ) , 1.0 / 3.0 );
	double viscosity =  sim_of_Cell::viscosity;  // 9.3 * 10.0


        double tmp = 3.0 / 4.0 * radius * radius * lenght;
        double d_v = 2.0 * std::pow( tmp , 1.0 / 3.0);
        double d_n = 2.0 * sqrt( 2.0 * radius * lenght / 3.14159265 );
        double d_s = sqrt( 2.0 * radius * radius + radius * lenght );			//tady dvojka neni umyslne
        double K = 1.0 / 3.0 * d_n / d_v + 2.0 / 3.0 * d_n / d_v;

        double effectiveFriction = 3.0 * sim_of_Cell::PI  * radius_sphere_same_volume * K * 10; //* viscosity
	return effectiveFriction;
}



double calculateEffectiveFriction_all_Mitochondria( unsigned int microtubule_number )
{
	double effectiveFriction = calculateEffectiveFriction_Mitochondria( 1.5e-6 , 0.4e-6);
	return ( double )microtubule_number * effectiveFriction;
}


double calculateEffectiveFriction_Endoplasmatic_Reticulum( double radius_of_the_cell )
{
	double Volume_Cell = 4.0 / 3.0 * sim_of_Cell::PI * radius_of_the_cell * radius_of_the_cell * radius_of_the_cell;
	double Endoplasmatic_Reticulum_Volume = 0.1 * Volume_Cell;

	double Tmp_Endoplasmatic_Reticulum_Volume = 3.0 / 4.0 * Endoplasmatic_Reticulum_Volume / sim_of_Cell::PI;
	double d_v = 2.0 * std::pow( Tmp_Endoplasmatic_Reticulum_Volume , 1.0 / 3.0 );

	double Area_Surface_Cell = 4 * sim_of_Cell::PI * radius_of_the_cell * radius_of_the_cell;
	double Area_Surface_Endoplasmatic_Reticulum = 20.0 * Area_Surface_Cell;

	double Tmp_Area_Surface_Endoplasmatic_Reticulum = Area_Surface_Endoplasmatic_Reticulum / sim_of_Cell::PI / 4.0;
	double d_s = sqrt( Tmp_Area_Surface_Endoplasmatic_Reticulum );
	double d_n = d_s;
        double K = 1.0 / 3.0 * d_n / d_v + 2.0 / 3.0 * d_s / d_v;


        double effectiveFriction = 3.0 * sim_of_Cell::PI * sim_of_Cell::viscosity * d_v * K;
	return effectiveFriction;
}


double calculate_form_Friction_Rough_ER_Golgi( double radius_of_the_cell )
{
	double Volume_Cell = 4.0 / 3.0 * sim_of_Cell::PI * radius_of_the_cell * radius_of_the_cell * radius_of_the_cell;
	double Rough_Golgi = 0.09 * Volume_Cell;
	double Tmp_Rough_Golgi = 3.0 / 4.0 * Rough_Golgi / sim_of_Cell::PI;
	double d_v = 2.0 * std::pow( Tmp_Rough_Golgi , 1.0 / 3.0 ); 

/*
	double Tmp_Endoplasmatic_Reticulum_Volume = 3.0 / 4.0 * Endoplasmatic_Reticulum_Volume / sim_of_Cell::PI;
	double d_v = 2.0 * std::pow( Tmp_Endoplasmatic_Reticulum_Volume , 1.0 / 3.0 );

	double Area_Surface_Cell = 4 * sim_of_Cell::PI * radius_of_the_cell * radius_of_the_cell;
	double Area_Surface_Endoplasmatic_Reticulum = 20.0 * Area_Surface_Cell;

	double Tmp_Area_Surface_Endoplasmatic_Reticulum = Area_Surface_Endoplasmatic_Reticulum / sim_of_Cell::PI / 4.0;
	double d_s = sqrt( Tmp_Area_Surface_Endoplasmatic_Reticulum );
	double d_n = d_s;
        double K = 1.0 / 3.0 * d_n / d_v + 2.0 / 3.0 * d_s / d_v;


        double effectiveFriction = 3.0 * sim_of_Cell::PI * sim_of_Cell::viscosity * d_v * K;
*/
	return 0.0;


}





double calculate_Friction_Smooth_ER_Golgi( double radius_of_the_cell )
{
	double Volume_Cell = 4.0 / 3.0 * sim_of_Cell::PI * radius_of_the_cell * radius_of_the_cell * radius_of_the_cell;
	double Rough_Golgi = 0.06 * Volume_Cell;
	double Tmp_Rough_Golgi = 3.0 / 4.0 * Rough_Golgi / sim_of_Cell::PI;
	double d_v = 2.0 * std::pow( Tmp_Rough_Golgi , 1.0 / 3.0 ); 

	double Area_Surface_Cell = 4 * sim_of_Cell::PI * radius_of_the_cell * radius_of_the_cell;
	double Area_Surface_Rough_Golgi = 23.0 / 2.0 * Area_Surface_Cell;

	double Tmp_Area_Surface_Endoplasmatic_Reticulum = Area_Surface_Rough_Golgi / sim_of_Cell::PI / 4.0;
	double d_s = sqrt( Tmp_Area_Surface_Endoplasmatic_Reticulum );
	double d_n = d_v * 2.0;
        double K = 1.0 / 3.0 * d_n / d_v + 2.0 / 3.0 * d_s / d_v;


        double effectiveFriction = 3.0 * sim_of_Cell::PI * 10.0 * d_v * K; //* sim_of_Cell::deeper_viscosity
	return effectiveFriction;

}


double calculate_Friction_Rough_ER( double radius_of_the_cell )
{
	double Volume_Cell = 4.0 / 3.0 * sim_of_Cell::PI * radius_of_the_cell * radius_of_the_cell * radius_of_the_cell;
	double Rough_Golgi = 0.09 * Volume_Cell;
	double Tmp_Rough_Golgi = 3.0 / 4.0 * Rough_Golgi / sim_of_Cell::PI;
	double d_v = 2.0 * std::pow( Tmp_Rough_Golgi , 1.0 / 3.0 ); 

	double Area_Surface_Cell = 4 * sim_of_Cell::PI * radius_of_the_cell * radius_of_the_cell;
	double Area_Surface_Rough_Golgi = 35.0 /2.0 * Area_Surface_Cell;

	double Tmp_Area_Surface_Endoplasmatic_Reticulum = Area_Surface_Rough_Golgi / sim_of_Cell::PI / 4.0;
	double d_s = sqrt( Tmp_Area_Surface_Endoplasmatic_Reticulum );
	double d_n = d_v * 2.0;
        double K = 1.0 / 3.0 * d_n / d_v + 2.0 / 3.0 * d_s / d_v;


        double effectiveFriction = 3.0 * sim_of_Cell::PI * 10.0 * d_v * K; //* sim_of_Cell::deeper_viscosity
	return effectiveFriction;

}




double calculateEffectiveFriction_Golgi_apparatus( double radius_of_the_cell )
{
	return 2.0 * calculateEffectiveFriction_Endoplasmatic_Reticulum( radius_of_the_cell );
}




double Endoplasmatic_Reticulum_coef( double radius_of_the_cell )
{
	double Volume_Cell = 4.0 / 3.0 * sim_of_Cell::PI * radius_of_the_cell * radius_of_the_cell * radius_of_the_cell;
	double Endoplasmatic_Reticulum_Volume = 0.1 * Volume_Cell;

	double Tmp_Endoplasmatic_Reticulum_Volume = 3.0 / 4.0 * Endoplasmatic_Reticulum_Volume / sim_of_Cell::PI;
	double d_v = 2.0 * std::pow( Tmp_Endoplasmatic_Reticulum_Volume , 1.0 / 3.0 );

	double Area_Surface_Cell = 4 * sim_of_Cell::PI * radius_of_the_cell * radius_of_the_cell;
	double Area_Surface_Endoplasmatic_Reticulum = 20.0 * Area_Surface_Cell;

	double Tmp_Area_Surface_Endoplasmatic_Reticulum = Area_Surface_Endoplasmatic_Reticulum / sim_of_Cell::PI / 4.0;
	double d_s = sqrt( Tmp_Area_Surface_Endoplasmatic_Reticulum );
	double d_n = d_s;
        double K = 1.0 / 3.0 * d_n / d_v + 2.0 / 3.0 * d_s / d_v;


        double effectiveFriction = 3.0 * sim_of_Cell::PI * d_v * K;
	return effectiveFriction;
}















double calculate_number_of_mitochondria( double radius_of_cell , double procentage , double lenght_of_mitochondria , double radius_of_mitochondria )
{
	double Surface_of_mitochondria_cut = 2.0 * radius_of_mitochondria * lenght_of_mitochondria;
	double Surface_of_cell = 4.0 * sim_of_Cell::PI * radius_of_cell * radius_of_cell;

	unsigned int number_of_mitochondrion = ( unsigned int ) ( Surface_of_cell * procentage / Surface_of_mitochondria_cut );
	return number_of_mitochondrion;
}


double calculate_number_of_mitochondria_micro_lengths( unsigned int number_of_microtubule , double lenght_of_microtubule , double lenght_of_mitochondria , double procentage )
{
	double number_of_mitochondria = number_of_microtubule * lenght_of_microtubule * procentage / lenght_of_mitochondria;
	return number_of_mitochondria;
}


double mito_cell_volume( double radius_of_cell , double lenght_of_mitochondria , double radius_of_mitochondria , double number_of_mitochondria )
{
	double Volume_cell = sim_of_Cell::PI * radius_of_cell * radius_of_cell * radius_of_cell * 4.0 / 3.0;
	double Volume_mito = sim_of_Cell::PI * radius_of_mitochondria * radius_of_mitochondria * lenght_of_mitochondria * number_of_mitochondria;
	return Volume_mito / Volume_cell;
}


double compute_time( double radius_of_the_cell )
{
	double lenght_of_mito = 1.4e-6;
	double radius_of_mito = 0.4e-6;
	double number_of_microtubule = 180.0;
	double number_of_mitochondrion = 180.0;

	double trajectory_lenght = sim_of_Cell::PI * radius_of_the_cell;

	unsigned int number_of_motors = 100;
	double mean_force_of_motor = 1.5e-12;


	Microtubule tmp( 20 );
	double friction_micro = tmp.get_effective_friction_whole_microtubule();
	double mitochondrion_friction = calculateEffectiveFriction_Mitochondria( lenght_of_mito, radius_of_mito );
	double effective_friction_endo_ret = calculateEffectiveFriction_Endoplasmatic_Reticulum( radius_of_the_cell );
	double effective_friction_Golgi = effective_friction_endo_ret / 2.0;

	double total_friction = friction_micro * number_of_microtubule + mitochondrion_friction * number_of_mitochondrion + effective_friction_Golgi + effective_friction_endo_ret;
	double total_force = number_of_motors * mean_force_of_motor;
	double V = total_force / total_friction;
	double time = trajectory_lenght / V;
	return time;
}

void erase_distance_measurements()
{
	remove( "picturesVideos/distance_results/distance_MTOC_IS_cortical_sliding_rear_plane.txt" );
	remove( "picturesVideos/distance_results/distance_MTOC_center_IS_cortical_sliding_rear_plane.txt" );
	remove( "picturesVideos/distance_results/distance_MTOC_IS_cortical_sliding_front_plane.txt" );
	remove( "picturesVideos/distance_results/distance_MTOC_center_IS_cortical_sliding_front_plane.txt" );
	remove( "picturesVideos/distance_results/distance_MTOC_IS_cortical_sliding_center_front.txt" );
	remove( "picturesVideos/distance_results/distance_MTOC_center_IS_cortical_sliding_center_front.txt" );
	remove( "picturesVideos/distance_results/distance_MTOC_IS_cortical_sliding_rear.txt" );
	remove( "picturesVideos/distance_results/distance_MTOC_center_IS_cortical_sliding_center_rear.txt" );

	remove( "picturesVideos/distance_results/distance_MTOC_IS_capture_shrinkage_front_plane.txt" );
	remove( "picturesVideos/distance_results/distance_MTOC_center_IS_capture_shrinkage_center_of_front_plane.txt" );
	remove( "picturesVideos/distance_results/distance_MTOC_center_IS_capture_shrinkage_front_plane.txt" );
	remove( "picturesVideos/distance_results/center_MTOC.txt" );
    remove( "picturesVideos/distance_results/time.txt" );
}




//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS

#include "Cell.h"
#include <tuple>

Cell::Cell(  )
{
	this->a_axis = Cell_parametres::A_AXIS;
	this->b_axis = Cell_parametres::B_AXIS;
        this->number_of_microtubules_extra = 0;
	this->nucleus = Nucleus();
	//initializing number of microtubules
	this->number_of_microtubules = 16;
	//create array of microtubules
	this->array_Of_Microtubules = new Microtubule[ this->number_of_microtubules ];
	//this is unconstrained space - radius_of_Cell has no meaning
	//IS_position is in x z plane - angle is determined in sim_of_Cell::IS_angle
	//In unconstrained cell IS has no meaning
	this->IS_Capture_Shrinkage = IS_Capt_Shrinkage();
	//this->IS_Cortical_Sliding = ISCorticalSliding();
        this->IS_Cortic_Sl = ISCorticalSl();
        this->MTOC = MTOC2( this->number_of_microtubules );
        //this->abstract_MTOC = ideal_MTOC();
        this->cell_id = 0;
}













Cell::Cell( unsigned int number_Of_Microtubules , unsigned int number_Of_extra_Microtubules , unsigned int servant )
{
        this->cell_id = servant;
	std::uniform_real_distribution<> distribution{ 0 , 1 };
	std::uniform_real_distribution<> distribution_angles{ 0 , 2.0 *  sim_of_Cell::PI };
	unsigned int number_of_generator = omp_get_thread_num();

        this->a_axis = Cell_parametres::A_AXIS;
	this->b_axis = Cell_parametres::B_AXIS;

	Vector3d center_of_nucleus( 0.0 , 0.0 , Nucleus_parametres::z_coordinate ); //1.5e-6
	this->nucleus = Nucleus( center_of_nucleus );

	if( ( this->a_axis < this->nucleus.get_A_Axis() ) || ( this->b_axis < this->nucleus.get_B_Axis() ) )
	{
		cout<<"( this->a_axis < this->nucleus.get_A_Axis() ) || ( this->b_axis < this->nucleus.get_B_Axis() ) in Cell( unsigned int number_Of_Microtubules )"<<endl;
		cout<<"ERROR_ID = 976456864648645"<<endl;
		throw("");
	}
	if( number_Of_Microtubules % 4 != 0 )
    	{
		cout<<"number_Of_Microtubules % 4 != 0 "<<endl;
        	cout<<"Cell::Cell(  unsigned int number_Of_Microtubules , unsigned int number_Of_extra_Microtubules , unsigned int servant )"<<endl;
		cout<<"ERROR_ID = 8764356435646"<<endl;
		throw("");
    	}

    	this->number_of_microtubules = number_Of_Microtubules + number_Of_extra_Microtubules;
    	this->number_of_microtubules_extra = number_Of_extra_Microtubules;


	//////////////////////////////////////////////////////////Capture-Shrinkage////////////////////////////////////////////////////////
	//HERE EVERYTHING ABOUT IS - MICROTUBULE CATCHING WILL BE SET
 	double azimutal_angle = distribution_angles( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );

	azimutal_angle = IS_Capture_shrinkage_param::azimutal_angle;
	double planar_component = sin( IS_Capture_shrinkage_param::polar_angle );
	double vertical_component = cos( IS_Capture_shrinkage_param::polar_angle );
	double planar_component_x = planar_component * cos( azimutal_angle );
	double planar_component_y = planar_component * sin( azimutal_angle );

	Vector3d orientation( 0.0 , 0.0 , 0.0 );
	orientation( 0 ) = planar_component_x;
	orientation( 1 ) = planar_component_y;
	orientation( 2 ) = vertical_component;
	orientation = orientation / orientation.norm();
	
	//////////////////////////////////////////////////////////Capture-Shrinkage////////////////////////////////////////////////////////
	double radius_of_IS_argument = IS_Capture_shrinkage_param::radius;
	Vector3d  center_of_IS_cap_sh_fron_center = IS_Capture_shrinkage_param::z_coordinate_front * orientation;
	Vector3d  center_of_IS_cap_sh_rear_center = IS_Capture_shrinkage_param::z_coordinate_back * orientation;
	this->IS_Capture_Shrinkage = IS_Capt_Shrinkage( center_of_IS_cap_sh_fron_center , radius_of_IS_argument , center_of_IS_cap_sh_rear_center );


	//////////////////////////////////////////////////////////Cortical Sliding////////////////////////////////////////////////////////
    	Vector3d front_center = orientation * IS_Cortical_Sl_parameter::front_radius;
    	Vector3d rear_center = orientation * IS_Cortical_Sl_parameter::rear_radius;
    	double radius_outer = IS_Cortical_Sl_parameter::radius;
    	double radius_inner = IS_Cortical_Sl_parameter::radius_inner;
    	this->IS_Cortic_Sl = ISCorticalSl( front_center , rear_center , radius_outer , radius_inner );


        unsigned int number_Of_Microtubules_MTOC;
        if( number_Of_Microtubules < MTOCparam::boundary_MTOC_points )
        {
            number_Of_Microtubules_MTOC = number_Of_Microtubules;
        }
        else
        {
             number_Of_Microtubules_MTOC = MTOCparam::boundary_MTOC_points;
        }


	//CAREFULL ABOUT THIS MAGIC NUMBERS!!!!
   	double z_coordinate = this->nucleus.get_B_Axis() -0.2e-6;
   	Vector3d center_of_MTOC( 0.0 , 0.0 , z_coordinate );
   	this->MTOC = MTOC2( number_Of_Microtubules_MTOC , center_of_MTOC , "plane" );
   	this->mtoc = MTOC_simple( center_of_MTOC );


   	unsigned int side = 4;
   	unsigned int micro_per_side = number_Of_Microtubules / side;

        unsigned int polygon_per_side;
   	//polygon_per_side gives me the number of polygon at every end of the centrioles

   	if( micro_per_side % MTOCparam::micro_in_polygon == 0 )
   	{
       	polygon_per_side = micro_per_side / MTOCparam::micro_in_polygon;
   	}
   	else
   	{
       	polygon_per_side = micro_per_side / MTOCparam::micro_in_polygon + 1;
   	}

   unsigned int unfinished_polygon_number = micro_per_side % MTOCparam::micro_in_polygon;
   this->array_Of_Microtubules = new Microtubule[ this->number_of_microtubules ];

   Vector3d MTOC_center = this->MTOC.get_center();
   unsigned int counter_micro = 0;
   unsigned int polygon_counter = 0;


   double constant = 1;
   Vector3d center_MTOC_tmp = MTOC_center;


   Vector3d original_orientation = this->MTOC.get_center() / this->MTOC.get_center().norm();
   this->abstract_MTOC = ideal_MTOC( this->number_of_microtubules , original_orientation );
   this->abstract_MTOC.set_center( this->MTOC.get_center() );
   this->abstract_MTOC.set_orientation_center( this->MTOC.get_center() );
   this->abstract_MTOC.set_mtoc_center( this->MTOC.get_center() );
   this->abstract_MTOC.set_original_MTOC_center( this->MTOC.get_center() );
   this->abstract_MTOC.set_original_MTOC_1_point( this->MTOC.get_point( 1 ) );


   for( unsigned int strana = 0 ; strana < side ; strana ++ )
   {
	unsigned int beads_per_side = 0;
        for( unsigned int polygon = 0 ; polygon < polygon_per_side ; polygon ++ )
        {
            unsigned int meze;
            if( polygon != polygon_per_side - 1 )
            {
                meze = MTOCparam::micro_in_polygon;
            }
            else
            {
                if( unfinished_polygon_number == 0 )
                {
                    meze = MTOCparam::micro_in_polygon;
                }
                else
                {
                    meze = unfinished_polygon_number;
                }
            }
            for( unsigned int  microtubule = 0 ; microtubule < meze ; microtubule ++ )
            {
                   unsigned int point_number = strana * MTOCparam::micro_in_polygon  + microtubule + 1;
		   unsigned int opposite_point_number = this->MTOC.get_random_point_on_opposite_side_without_bias( point_number );
                   Vector3d MTOC_Point = this->MTOC.get_point( point_number );
                   Vector3d opposite_MTOC_Point = this->MTOC.get_point( opposite_point_number );

                   this->abstract_MTOC.set_mtoc_point(  counter_micro , MTOC_Point  );
                   Vector3d orientation = MTOC_Point - opposite_MTOC_Point;
                   orientation = orientation / orientation.norm();


                   Vector3d first_Point = opposite_MTOC_Point;
                   Vector3d second_Point = MTOC_Point;
		   std::uniform_int_distribution<int> distribution( 0 , sim_of_Cell::REDUCTION );
    		   unsigned int number_of_generator = omp_get_thread_num();
    		   unsigned int reduction = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );

                   unsigned int number_of_points = sim_of_Cell::MicrotubulePoints;
                   number_of_points = sim_of_Cell::MicrotubulePoints - reduction;

                   Microtubule tmp( first_Point , second_Point , orientation , microtubule , polygon_counter , strana , point_number , opposite_point_number , number_of_points );
                   this->array_Of_Microtubules[ counter_micro ] = tmp;
                   this->abstract_MTOC.set_orientation( counter_micro , this->array_Of_Microtubules[ counter_micro ].getPoint( 1 ) );


                   Vector3d mtoc_point = this->MTOC.get_center() - orientation * sim_of_Cell::resting_distance;
                   this->abstract_MTOC.set_point( counter_micro , mtoc_point );
                   counter_micro = counter_micro + 1;

            }


            polygon_counter = polygon_counter + 1;
        }
    }
        //Here I will create extra polygon
        //extra microtubules will be added to one side
        unsigned int full_polygon_number = number_Of_extra_Microtubules / MTOCparam::micro_in_polygon;
        unsigned int unfinished_polygon_mito = number_Of_extra_Microtubules - full_polygon_number * MTOCparam::micro_in_polygon;
        this->abstract_MTOC.set_true_mtoc_coordinates( this->MTOC.get_coordinates() );
}


Cell::Cell( unsigned int number_Of_Microtubules , unsigned int number_Of_extra_Microtubules , unsigned int servant  , unsigned int servant_2 )
{
        this->cell_id = servant;
	std::uniform_real_distribution<> distribution{ 0 , 1 };
	std::uniform_real_distribution<> distribution_angles{ 0 , 2.0 *  sim_of_Cell::PI };
	unsigned int number_of_generator = omp_get_thread_num();


        this->a_axis = Cell_parametres::A_AXIS;
	this->b_axis = Cell_parametres::B_AXIS;

	Vector3d center_of_nucleus( 0.0 , 0.0 , Nucleus_parametres::z_coordinate ); //1.5e-6
	this->nucleus = Nucleus( center_of_nucleus );

	if( ( this->a_axis < this->nucleus.get_A_Axis() ) || ( this->b_axis < this->nucleus.get_B_Axis() ) )
	{
		cout<<"( this->a_axis < this->nucleus.get_A_Axis() ) || ( this->b_axis < this->nucleus.get_B_Axis() ) in Cell( unsigned int number_Of_Microtubules )"<<endl;
		cout<<"ERROR_ID = 976456864648645"<<endl;
		throw("");
	}
	if( number_Of_Microtubules % 4 != 0 )
    {
		cout<<"number_Of_Microtubules % 4 != 0 "<<endl;
        cout<<"Cell::Cell(  unsigned int number_Of_Microtubules , unsigned int number_Of_extra_Microtubules , unsigned int servant )"<<endl;
		cout<<"ERROR_ID = 8764356435646"<<endl;
		throw("");
    }


    this->number_of_microtubules = number_Of_Microtubules + number_Of_extra_Microtubules;
    this->number_of_microtubules_extra = number_Of_extra_Microtubules;



 	double azimutal_angle = distribution_angles( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );


	azimutal_angle = 0;
	IS_Capture_shrinkage_param::azimutal_angle;
	Vector3d orientation( 0.0 , 0.0 , 0.0 );
	double planar_component = sin( IS_Capture_shrinkage_param::polar_angle );
	double vertical_component = cos( IS_Capture_shrinkage_param::polar_angle );
	double planar_component_x = planar_component * cos( azimutal_angle );
	double planar_component_y = planar_component * sin( azimutal_angle );


	orientation( 0 ) = planar_component_x;
	orientation( 1 ) = planar_component_y;
	orientation( 2 ) = vertical_component;
	orientation = orientation / orientation.norm();

	//////////////////////////////////////////////////////////Capture-Shrinkage////////////////////////////////////////////////////////
	double radius_of_IS_argument = IS_Capture_shrinkage_param::radius;
	Vector3d  center_of_IS_cap_sh_fron_center = IS_Capture_shrinkage_param::z_coordinate_front * orientation;

	Vector3d  center_of_IS_cap_sh_rear_center = IS_Capture_shrinkage_param::z_coordinate_back * orientation;
	this->IS_Capture_Shrinkage = IS_Capt_Shrinkage( center_of_IS_cap_sh_fron_center , radius_of_IS_argument , center_of_IS_cap_sh_rear_center );

	//////////////////////////////////////////////////////////Cortical Sliding////////////////////////////////////////////////////////

    	Vector3d front_center = orientation * IS_Cortical_Sl_parameter::front_radius;
    	Vector3d rear_center = orientation * IS_Cortical_Sl_parameter::rear_radius;
    	double radius_outer = IS_Cortical_Sl_parameter::radius;
    	double radius_inner = IS_Cortical_Sl_parameter::radius_inner;
    	this->IS_Cortic_Sl = ISCorticalSl( front_center , rear_center , radius_outer , radius_inner );


        unsigned int number_Of_Microtubules_MTOC;
        if( number_Of_Microtubules < MTOCparam::boundary_MTOC_points )
        {
            number_Of_Microtubules_MTOC = number_Of_Microtubules;
        }
        else
        {
             number_Of_Microtubules_MTOC = MTOCparam::boundary_MTOC_points;
        }


   double z_coordinate = this->nucleus.get_B_Axis() -0.2e-6;
   Vector3d center_of_MTOC( 0.0 , 0.0 , z_coordinate );
   this->MTOC = MTOC2( number_Of_Microtubules_MTOC , center_of_MTOC , "plane" , "plane");
   this->mtoc = MTOC_simple( center_of_MTOC );


   unsigned int side = 4;
   unsigned int micro_per_side = number_Of_Microtubules / side;
   unsigned int polygon_per_side;

   if( micro_per_side % MTOCparam::micro_in_polygon == 0 )
   {
       polygon_per_side = micro_per_side / MTOCparam::micro_in_polygon;
   }
   else
   {
       polygon_per_side = micro_per_side / MTOCparam::micro_in_polygon + 1;
   }

   unsigned int unfinished_polygon_number = micro_per_side % MTOCparam::micro_in_polygon;
   this->array_Of_Microtubules = new Microtubule[ this->number_of_microtubules ];

   Vector3d MTOC_center = this->MTOC.get_center();
   unsigned int counter_micro = 0;
   unsigned int polygon_counter = 0;


   double constant = 1;
   Vector3d center_MTOC_tmp = MTOC_center;


   Vector3d original_orientation = this->MTOC.get_center() / this->MTOC.get_center().norm();
   this->abstract_MTOC = ideal_MTOC( this->number_of_microtubules , original_orientation );
   this->abstract_MTOC.set_center( this->MTOC.get_center() );
   this->abstract_MTOC.set_orientation_center( this->MTOC.get_center() );
   this->abstract_MTOC.set_mtoc_center( this->MTOC.get_center() );
   this->abstract_MTOC.set_original_MTOC_center( this->MTOC.get_center() );
   this->abstract_MTOC.set_original_MTOC_1_point( this->MTOC.get_point( 1 ) );


   for( unsigned int strana = 0 ; strana < side ; strana ++ )
   {
	unsigned int beads_per_side = 0;
        for( unsigned int polygon = 0 ; polygon < polygon_per_side ; polygon ++ )
        {
            unsigned int meze;
            if( polygon != polygon_per_side - 1 )
            {
                meze = MTOCparam::micro_in_polygon;
            }
            else
            {
                if( unfinished_polygon_number == 0 )
                {
                    meze = MTOCparam::micro_in_polygon;
                }
                else
                {
                    meze = unfinished_polygon_number;
                }
            }
            for( unsigned int  microtubule = 0 ; microtubule < meze ; microtubule ++ )
            {
                   unsigned int point_number = strana * MTOCparam::micro_in_polygon  + microtubule + 1;
		   unsigned int opposite_point_number = this->MTOC.get_random_point_on_opposite_side_without_bias( point_number );
                   Vector3d MTOC_Point = this->MTOC.get_point( point_number );
                   Vector3d opposite_MTOC_Point = this->MTOC.get_point( opposite_point_number );

                   this->abstract_MTOC.set_mtoc_point(  counter_micro , MTOC_Point  );
                   Vector3d orientation = MTOC_Point - opposite_MTOC_Point;
                   orientation = orientation / orientation.norm();


                   Vector3d first_Point = opposite_MTOC_Point;
                   Vector3d second_Point = MTOC_Point;
		   std::uniform_int_distribution<int> distribution( 0 , sim_of_Cell::REDUCTION );
    		   unsigned int number_of_generator = omp_get_thread_num();
    		   unsigned int reduction = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );

                   unsigned int number_of_points = sim_of_Cell::MicrotubulePoints;
                   number_of_points = sim_of_Cell::MicrotubulePoints - reduction;
                   Microtubule tmp( first_Point , second_Point , orientation , microtubule , polygon_counter , strana , point_number , opposite_point_number , number_of_points );//true

                   this->array_Of_Microtubules[ counter_micro ] = tmp;
                   this->abstract_MTOC.set_orientation( counter_micro , this->array_Of_Microtubules[ counter_micro ].getPoint( 1 ) );

                   Vector3d mtoc_point = this->MTOC.get_center() - orientation * sim_of_Cell::resting_distance;
                   this->abstract_MTOC.set_point( counter_micro , mtoc_point );
                   counter_micro = counter_micro + 1;

            }

            polygon_counter = polygon_counter + 1;
        }
    }
        unsigned int full_polygon_number = number_Of_extra_Microtubules / MTOCparam::micro_in_polygon;
        unsigned int unfinished_polygon_mito = number_Of_extra_Microtubules - full_polygon_number * MTOCparam::micro_in_polygon;

        this->abstract_MTOC.set_true_mtoc_coordinates( this->MTOC.get_coordinates() );
}









Microtubule Cell::getMicrotubule( unsigned int index )
{
	if( index >= this->number_of_microtubules )
	{
		cout<<"index >= this->number_of_microtubules in Cell::getMicrotubule( unsigned int index )"<<endl;
		throw("");
	}

	return this->array_Of_Microtubules[ index ];
}


unsigned int Cell::get_microtubule_number()
{
    return this->number_of_microtubules;
}

double Cell::get_cytoskeleton_friction()
{

     double viscosity =  sim_of_Cell::viscosity;  // 9.3 * 10.0
     double radius = 1.25e-8;
     double microtubule_lenght = sim_of_Cell::resting_distance * (double)( sim_of_Cell::MicrotubulePoints - 1 );
     double nominator = 4.0 * sim_of_Cell::PI  * microtubule_lenght; //* viscosity
     double denominator = log( microtubule_lenght / ( 2.0 * radius ) ) + 0.84;
     double effectiveFriction_micro = nominator / denominator;
     return effectiveFriction_micro * this->get_microtubule_number();

}


void Cell::resize_microtubules( )
{
    for( unsigned int index = 0 ; index < this->number_of_microtubules ; index ++ )
	{
       	this->array_Of_Microtubules[ index ].resizeMicrotubule();
    }


}

unsigned int Cell::get_number_of_microtubule_with_dynein_index( unsigned int dynein_index )
{
    unsigned int number = 0;
    for( unsigned int microtubule = 0 ; microtubule < this->get_microtubule_number() ; microtubule ++  )
    {
        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() ==  dynein_index )
        {
            number = number + 1;
        }
    }
    return number;

}

IS_Capt_Shrinkage Cell::get_IS_capture_shrinkage()
{
	IS_Capt_Shrinkage tmp;
	tmp = this->IS_Capture_Shrinkage;
	return tmp;
}


Vector3d Cell::force_simple_one_micro( unsigned int microtubule )
{
    //force goes from MTOC to micro
    Vector3d position_micro = this->array_Of_Microtubules[ microtubule ].getPoint( 0 );
    Vector3d position_mtoc = this->mtoc.get_position();
    Vector3d tangent = position_micro - position_mtoc;
    double distance = tangent.norm();
    double absolute_value_of_the_force = distance * MTOCparam::MTOC2_center_micro_kappa;
    Vector3d orientation = tangent / tangent.norm();
    Vector3d force = absolute_value_of_the_force * orientation;
    return force;
}

Vector3d Cell::force_nucleus( )
{
    Vector3d overall_force( 0.0 , 0.0 , 0.0 );

    Vector3d position = this->mtoc.get_position();
    Vector3d force_center = this->nucleus.force_position_nucleus_wall_MTOC( position );
    return force_center;
}
Vector3d Cell::force_wall( )
{
    Vector3d position = this->mtoc.get_position();
    Vector3d force = force_position_cell_wall_two_elipsoid( position  );
    return force;
}



void Cell::MidStep_3()
{

//MICROTUBULES
    MatrixXd* random_Forces = new MatrixXd[ this->number_of_microtubules ];
    MatrixXd* original_Coordinates = new MatrixXd[ this->number_of_microtubules ];
    MatrixXd* external_Forces = new MatrixXd[ this->number_of_microtubules ];
    MatrixXd* dynein_Forces = new MatrixXd[ this->number_of_microtubules ];

//MTOC
    MatrixXd MTOC_original_coordinates = this->MTOC.get_coordinates();
    MatrixXd MTOC_forces = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );
    MatrixXd MTOC_wall_force = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );
    MatrixXd MTOC_nucleus_force = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );



for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
{
    	original_Coordinates[ microtubule ] = this->array_Of_Microtubules[ microtubule ].getCoordinates();
}




for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
{
	if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() < 1 )
	{
            MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

	    //CELL WALL INTERACTION
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;


            // NUCLEUS INTERACTION
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

            //MTOC INTERACTION

            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );

            this->MTOC_microtubule_two_points_force_with_bending( first_bead_force , second_bead_force , bending_force_opposite_MTOC , bending_force_MTOC , microtubule );

            if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() > 2  )
            {
                unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
                unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

                for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
                {
                    external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                    if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                    {
                          external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                    }
                }
                for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
                {
                    MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension ) + bending_force_MTOC( dimension )  ;//

                    MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
                }
            }

		external_Forces[ microtubule ] = external_Force;
	}
        else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9 )
        {
            MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

	    //CELL WALL INTERACTION
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

		// NUCLEUS INTERACTION
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

            //MTOC INTERACTION
            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );

            this->MTOC_microtubule_two_points_force_with_bending( first_bead_force , second_bead_force , bending_force_opposite_MTOC , bending_force_MTOC , microtubule );


            unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();


              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {
                  external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                  if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                  {
                        external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                  }
              }
              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {

                  MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension )+ bending_force_MTOC( dimension );//

                  MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
              }

            external_Forces[ microtubule ] = external_Force;

        MatrixXd force_dynein_IS_surface = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
        this->array_Of_Microtubules[ microtubule ].force_dynein_abscissa_real_dynein( force_dynein_IS_surface );
        external_Force = external_Force + force_dynein_IS_surface;



       external_Forces[ microtubule ] = external_Force;
        //DYNEIN INTERACTION

    }

    else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
	{
        MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

		//CELL WALL INTERACTION
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

		// NUCLEUS INTERACTION
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

        //MTOC INTERACTION
            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );


	    if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() > 3 )
            {
            	this->MTOC_microtubule_two_points_force_with_bending( first_bead_force , second_bead_force , bending_force_opposite_MTOC , bending_force_MTOC , microtubule );

            	if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() > 2  )
            	{
            		unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            		unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

              	for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              	{
                  	external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                  	if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                  	{
                     	   external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                  	}
              	}
              	for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              	{
                 	 MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension ) + bending_force_MTOC( dimension )  ;//

                 	 MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
              	}
           	}
            }
	    else
            {

                if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                {
            		unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            		Vector3d second_bead_force = this->MTOC_microtubule_two_points( microtubule );
              		for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              		{
                 		external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                 		MTOC_forces(3* number_MTOC_point + dimension , 0 ) = MTOC_forces( 3 * number_MTOC_point + dimension , 0 ) + (- 1) * second_bead_force( dimension );//
              		}

                }


            }



            MatrixXd force_dynein_IS_surface = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->array_Of_Microtubules[ microtubule ].force_dynein_abscissa_real_dynein( force_dynein_IS_surface );
            dynein_Forces[ microtubule ] = force_dynein_IS_surface;

            external_Force = external_Force + force_dynein_IS_surface;
            external_Forces[ microtubule ] = external_Force;
    }

}


this->force_MTOC_cell_wall( MTOC_wall_force );
this->MTOC_nucleus_interaction( MTOC_nucleus_force );
MTOC_forces = MTOC_forces + MTOC_wall_force + MTOC_nucleus_force;




this->MTOC.oneStepMidStepAlgorithm_1_half( MTOC_forces );
for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ ) //this->number_of_microtubules
{
	MatrixXd randomForce = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
	this->array_Of_Microtubules[ microtubule ].oneStepMidStepAlgorithm_1_half_producing_random( randomForce , external_Forces[ microtubule ]  );
	random_Forces[ microtubule ] = randomForce;

}









MTOC_forces = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );
MTOC_wall_force = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );
MTOC_nucleus_force = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );
this->resize_micro_with_different_first_segment_and_MTOC(  );



for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
{
	if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() < 1 )
	{
        MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

		//CELL WALL INTERACTION
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

		// NUCLEUS INTERACTION
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

        //MTOC INTERACTION

            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );

            this->MTOC_microtubule_two_points_force_with_bending( first_bead_force , second_bead_force , bending_force_opposite_MTOC , bending_force_MTOC , microtubule );

            if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() > 2  )
            {
            unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {
                  external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                  if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                  {
                        external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                  }
              }
              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {
                  MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension ) + bending_force_MTOC( dimension )  ;//

                  MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
              }
            }

		external_Forces[ microtubule ] = external_Force;
	}


    else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9 )
    {

        MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

		//CELL WALL INTERACTION
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

		// NUCLEUS INTERACTION
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

        //MTOC INTERACTION
            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );

            this->MTOC_microtubule_two_points_force_with_bending( first_bead_force , second_bead_force , bending_force_opposite_MTOC , bending_force_MTOC , microtubule );


            unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {
                  external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                  if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                  {
                        external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                  }
              }
              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {

                  MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension )+ bending_force_MTOC( dimension );//

                  MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
              }

            external_Forces[ microtubule ] = external_Force;

        MatrixXd force_dynein_IS_surface = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
        this->array_Of_Microtubules[ microtubule ].force_dynein_abscissa_real_dynein( force_dynein_IS_surface );
        external_Force = external_Force + force_dynein_IS_surface;

       external_Forces[ microtubule ] = external_Force;


    }

    else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
	{
        MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

		//CELL WALL INTERACTION
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

		// NUCLEUS INTERACTION
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

        //MTOC INTERACTION
            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );

            this->MTOC_microtubule_two_points_force_with_bending( first_bead_force , second_bead_force , bending_force_opposite_MTOC , bending_force_MTOC , microtubule );


            unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {
                  external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                  if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                  {
                        external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                  }
              }
              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {

                  MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension )+ bending_force_MTOC( dimension );//

                  MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
              }


            MatrixXd force_dynein_IS_surface = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->array_Of_Microtubules[ microtubule ].force_dynein_abscissa_real_dynein( force_dynein_IS_surface );

            external_Forces[ microtubule ] = external_Force;
    }

}


this->force_MTOC_cell_wall( MTOC_wall_force );
this->MTOC_nucleus_interaction( MTOC_nucleus_force );
MTOC_forces = MTOC_forces + MTOC_wall_force + MTOC_nucleus_force;
this->MTOC.oneStepMidStepAlgorithm_2_half( MTOC_forces , MTOC_original_coordinates );


for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ ) //this->number_of_microtubules
{
    if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
    {
        MatrixXd force = external_Forces[ microtubule ] + dynein_Forces[ microtubule ];
        this->array_Of_Microtubules[ microtubule ].oneStepMidStepAlgorithm_2_half( original_Coordinates[ microtubule ] , random_Forces[ microtubule ] , force );
    }
    else
    {
        this->array_Of_Microtubules[ microtubule ].oneStepMidStepAlgorithm_2_half( original_Coordinates[ microtubule ] , random_Forces[ microtubule ] , external_Forces[ microtubule ] );
    }
}

    delete[] random_Forces;
    random_Forces = NULL;
    delete[] original_Coordinates;
    original_Coordinates = NULL;
    delete[] external_Forces;
    external_Forces = NULL;
    delete[] dynein_Forces;
    dynein_Forces = NULL;



}


























void Cell::Euler_algorithm()
{

//MICROTUBULES
    MatrixXd* random_Forces = new MatrixXd[ this->number_of_microtubules ];
    MatrixXd* original_Coordinates = new MatrixXd[ this->number_of_microtubules ];
    MatrixXd* external_Forces = new MatrixXd[ this->number_of_microtubules ];

//MTOC
    MatrixXd MTOC_original_coordinates = this->MTOC.get_coordinates();
    MatrixXd MTOC_forces = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );
    MatrixXd MTOC_wall_force = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );
    MatrixXd MTOC_nucleus_force = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );






for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
{
	MatrixXd randomForce = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
	this->array_Of_Microtubules[ microtubule ].getRandomForces(  sim_of_Cell::time_Step , randomForce );
	random_Forces[ microtubule ] = randomForce;

	if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() < 1 )
	{
		original_Coordinates[ microtubule ] = this->array_Of_Microtubules[ microtubule ].getCoordinates();
        MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

		//CELL WALL INTERACTION
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

		// NUCLEUS INTERACTION
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

        //MTOC INTERACTION

            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );

            if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() > 2  )
            {
            unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {
                  external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                  if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                  {
                        external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                  }
              }

              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {

                  MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension ) + bending_force_MTOC( dimension )  ;//

                  MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
              }
            }

		external_Forces[ microtubule ] = external_Force;
	}
    else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 1 )
	{
		original_Coordinates[ microtubule ] = this->array_Of_Microtubules[ microtubule ].getCoordinates();

		//external forces will contain ALL( !!! ) forces from every structure in cell ( dynein , wall , MTOC )
		MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

		//CELL WALL INTERACTION
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;


		//NUCLEUS INTERACTION
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;



            Vector3d dyneinForceVector = this->array_Of_Microtubules[ microtubule ].dynein_force_capture_shrinkage( );
            unsigned int index_last_point = this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1;
            for( unsigned int i = 0 ; i < 3 ; i ++ )
            {
                external_Force( 3 * index_last_point + i , 0 ) = external_Force( 3 * index_last_point + i , 0 ) + dyneinForceVector( i );
            }

            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );


            unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {
                  external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                  if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                  {
                        external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                  }
              }
              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {

                  MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension ) + bending_force_MTOC( dimension );//

                  MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
              }
            external_Forces[ microtubule ] = external_Force;
	}
    else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 2 )
	{
        Vector3d force = MTOC_dynein_capture_shrinkage_force( microtubule );
        unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
        for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
        {
            MTOC_forces( 3 * number_MTOC_point + dimension , 0 ) = MTOC_forces( 3 * number_MTOC_point + dimension , 0 ) + force( dimension );
        }

	}

    else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9 )
    {
        original_Coordinates[ microtubule ] = this->array_Of_Microtubules[ microtubule ].getCoordinates();
        MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

		//CELL WALL INTERACTION
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

		// NUCLEUS INTERACTION
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

        //MTOC INTERACTION
            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );

            unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {
                  external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                  if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                  {
                        external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                  }
              }
              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {

                  MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension )+ bending_force_MTOC( dimension );//

                  MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
              }

            external_Forces[ microtubule ] = external_Force;

        MatrixXd force_dynein_IS_surface = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
        //this->array_Of_Microtubules[ microtubule ].force_dynein_abscissa( force_dynein_IS_surface );
        this->array_Of_Microtubules[ microtubule ].force_dynein_abscissa_real_dynein( force_dynein_IS_surface );
        external_Force = external_Force + force_dynein_IS_surface;

       external_Forces[ microtubule ] = external_Force;
        //DYNEIN INTERACTION

    }

    else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
	{
        original_Coordinates[ microtubule ] = this->array_Of_Microtubules[ microtubule ].getCoordinates();
        MatrixXd external_Force = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );

		//CELL WALL INTERACTION
            MatrixXd wall_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_cell_wall( wall_force_matrix , microtubule );
            external_Force = external_Force + wall_force_matrix;

		// NUCLEUS INTERACTION
            MatrixXd nucleus_force_matrix = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->force_whole_microtubule_nucleus( nucleus_force_matrix , microtubule );
            external_Force = external_Force + nucleus_force_matrix;

        //MTOC INTERACTION
            Vector3d first_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d second_bead_force( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_MTOC( 0.0 , 0.0 , 0.0 );
            Vector3d bending_force_opposite_MTOC( 0.0 , 0.0 , 0.0 );

            unsigned int number_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
            unsigned int num_opp_MTOC_point = this->array_Of_Microtubules[ microtubule ].get_MTOC_opposite_point();

              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {
                  external_Force( dimension , 0 ) = external_Force( dimension , 0 ) + first_bead_force( dimension );
                  if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() >= 2 )
                  {
                        external_Force( 3 + dimension , 0 ) = external_Force( 3 + dimension , 0 ) + second_bead_force( dimension );
                  }
              }
              for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
              {

                  MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( number_MTOC_point ) + dimension , 0 ) + ( - 1 ) * second_bead_force( dimension )+ bending_force_MTOC( dimension );//

                  MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) = MTOC_forces( 3 * ( num_opp_MTOC_point ) + dimension , 0 ) + ( - 1 ) * first_bead_force( dimension )+ bending_force_opposite_MTOC( dimension );//
              }


            MatrixXd force_dynein_IS_surface = MatrixXd::Zero( 3 * this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() , 1 );
            this->array_Of_Microtubules[ microtubule ].force_dynein_abscissa_real_dynein( force_dynein_IS_surface );
            external_Force = external_Force + force_dynein_IS_surface;
            external_Forces[ microtubule ] = external_Force;
    }

}

this->force_MTOC_cell_wall( MTOC_wall_force );
this->MTOC_nucleus_interaction( MTOC_nucleus_force );
MatrixXd MTOC_forces_naive = MatrixXd::Zero( 3 * ( this->MTOC.get_number_of_points() + 1 ) , 1 );

MTOC_forces = MTOC_forces + MTOC_wall_force + MTOC_nucleus_force;


this->MTOC.Euler_algorithm( MTOC_forces );


for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ ) //
{

	if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 2  )
	{
		continue;
	}
	else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 3 )
	{
		continue;
	}
	else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 5  )
	{
		continue;
	}
	else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 6 )
	{
		continue;
	}
	else
	{
			this->array_Of_Microtubules[ microtubule ].Euler_algorithm( random_Forces[ microtubule ] , external_Forces[ microtubule ] );
	}
}

    delete[] random_Forces;
    random_Forces = NULL;
    delete[] original_Coordinates;
    original_Coordinates = NULL;
    delete[] external_Forces;
    external_Forces = NULL;

}


void Cell::timeDevelopment_2( double Time_0 , double Time_1 )
{
	std::uniform_real_distribution<> distribution{0.0, 1.0};
        double numberOfSteps =  ( ( Time_1 - Time_0 ) / sim_of_Cell::time_Step );
	unsigned int numberOfSteps2 = nearbyint ( numberOfSteps ); //This is control for unprecise values of doubles in C++
	for( unsigned int step = 0 ; step < numberOfSteps2 ; step ++ )
	{

       		double time = step * sim_of_Cell::time_Step;
        	Cell_parametres::time = time;

        	if( step % 1400 == 0)
		{
            		cout<<"time = "<<time<<endl;
		}

		this->MidStep_3();
        	this->stepping_and_detachment_of_all_microtubule_projection_real_dynein();
		this->control_length_of_micro_IS_2();
        }

}



















void Cell::timeDevelopment_2( double Time_0 , double Time_1 , unsigned int number_of_Pictures )
{
	unsigned int key = 0;
        double numberOfSteps =  ( ( Time_1 - Time_0 ) / sim_of_Cell::time_Step );
	unsigned int numberOfSteps2 = nearbyint ( numberOfSteps ); //This is control for unprecise values of doubles in C++

	unsigned int picture_time_step = 1;

        cout<<"number_of_Pictures number_of_Pictures = "<<number_of_Pictures<<endl;
	if( numberOfSteps2 > number_of_Pictures )
	{
		picture_time_step = numberOfSteps2 / number_of_Pictures;
	}

	this->print_Cell_parametres();

	unsigned int printing_counter = 0;
	for( unsigned int step = 0 ; step < numberOfSteps2 ; step ++ )
	{
       		double time = step * sim_of_Cell::time_Step;
		this->MTOC.set_time_clock( time );

        	if( numberOfSteps2 >= number_of_Pictures )
		{
			if( step % picture_time_step == 0 )
			{
				this->print_Cell( step / picture_time_step );
			}
		}

        	if( step % 1000 == 0)
		{
            		cout<<"key = "<<key<<endl;
            		cout<<"time = "<<time<<endl;
        	}


      		if( step % 50 == 0 )
		{
			this->BIG_PRINT_numerical( step , key );	
        	}
        	
      		if( step % 50 == 0 )
		{
			//This is the distance control
                	double MTOC_IS_distance = this->get_MTOC_IS_distance();
                	if( MTOC_IS_distance < sim_of_Cell::treshold_distance )
                	{

                	}
              	}        	
        	




      		if( ( step %  100 ) == 0 ) //14285
		{
			if ( sim_of_Cell::density_of_dynein_motor_surface < 0.00001 )
			{	
				this->density_surface_dynein.change_part_of_IS_points( this->IS_Cortic_Sl );
			}
			else
			{
				throw("");
				this->density_surface_dynein.change_part_of_points_surface_and_IS( this->IS_Cortic_Sl );
			}
        	}



      		if( ( step %  7142 ) == 0 )
		{
			
			this->BIG_PRINT_numerical_micro( printing_counter , key );
			printing_counter = printing_counter + 1;
        	}

		//FUNCTIONS
        	//////////////////////////////////////
        	this->MidStep_3(); //

        	//////////////////////////////////////

        	this->stepping_and_detachment_of_all_microtubule_projection_real_dynein();
        	this->check_caught_micro_IS_with_real_dynein_all_micro( );
        	this->microtubule_catch_pair_abscissa_real_dynein_in_IS_all_micro( );
		this->catch_pair_abscissa_real_dynein();

		if( ( step %  20 ) == 0 )
		{
			this->control_length_of_micro_IS_2();		
		}
		
      		if( ( step %  100 ) == 0 ) //14285
		{
			this->test_dynein( "-111111111111111111111111111" );		
			this->density_surface_dynein.control_points_outside_IS(  this->IS_Cortic_Sl );		
			this->test_dynein( "000000000000000" );
		}
        }
        
        
}



void Cell::print_Cell( unsigned int index )
{

      char name_of_text_file [55];
	  sprintf ( name_of_text_file , "picturesVideos/textFiles/microtubule_%d.txt", index );

	  FILE *out;
	  out = fopen( name_of_text_file , "w");


	  for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
	  {
                if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 7 )
		{
			cout<<"this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 7"<<endl;
		}

		if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 0 )
         	{
                    Vector3d tangent_last;
                    try
                    {
                        tangent_last = this->array_Of_Microtubules[ microtubule ].getTangent2( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 );
                    }
                    catch( int e )
                    {
                        throw("");
                    }
		    if( tangent_last.norm() > 2 *  sim_of_Cell::resting_distance )
		    {
		    }
		    else
		    {

	                Vector3d center = this->MTOC.get_center();
	                double parametr = 0.0;

                	for( unsigned int point = 0 ; point < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() ; point ++  )
                	{
                    		parametr = parametr + 0.05;
                    		Vector3d coordinate_point = this->array_Of_Microtubules[ microtubule ].getPoint( point );
                    		fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  coordinate_point( 0 ) * 1e6 , coordinate_point( 1 ) * 1e6 , coordinate_point( 2 ) * 1e6  );
                	}
               		Vector3d tangent_last = this->array_Of_Microtubules[ microtubule ].getTangent2( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 );
                	Vector3d point_last = this->array_Of_Microtubules[ microtubule ].getPoint( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 );
                	Vector3d ill_point = point_last + tangent_last / 10.0;
                	parametr = parametr + 0.05;
                	fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  ill_point( 0 ) * 1e6 , ill_point( 1 ) * 1e6 , ill_point( 2 ) * 1e6  );
                	fprintf(out,"%u\n", this->array_Of_Microtubules[ microtubule ].get_dynein_index() );
                     	fprintf(out,"END\n" );
		   }

                }
	 else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 1 )
         {
                if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() > 1 )
                {
                    Vector3d center = this->MTOC.get_center();
                    double parametr = 0.0;

                    for( unsigned int point = 0 ; point < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() ; point ++  )
                    {
                        parametr = parametr + 0.05;
                            Vector3d coordinate_point = this->array_Of_Microtubules[ microtubule ].getPoint( point );
                        fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  coordinate_point( 0 ) * 1e6 , coordinate_point( 1 ) * 1e6 , coordinate_point( 2 ) * 1e6  );
                    }
                    Vector3d tangent_last = this->array_Of_Microtubules[ microtubule ].getTangent2( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 );
                    Vector3d point_last = this->array_Of_Microtubules[ microtubule ].getPoint( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 );
                    Vector3d ill_point = point_last + tangent_last / 10.0;
                    parametr = parametr + 0.05;
                    fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  ill_point( 0 ) * 1e6 , ill_point( 1 ) * 1e6 , ill_point( 2 ) * 1e6  );
                    fprintf(out,"%u\n", this->array_Of_Microtubules[ microtubule ].get_dynein_index() );
                    fprintf(out,"END\n" );
                }
                else
                {
                    Vector3d center = this->MTOC.get_center();
                    double parametr = 0.0;
                    for( unsigned int point = 0 ; point < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() ; point ++  )
                    {
                            parametr = parametr + 0.05;
                    Vector3d coordinate_point = this->array_Of_Microtubules[ microtubule ].getPoint( point );
                        fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  coordinate_point( 0 ) * 1e6 , coordinate_point( 1 ) * 1e6 , coordinate_point( 2 ) * 1e6  );
                    }

                    Vector3d catching_point = this->array_Of_Microtubules[ microtubule ].get_IS_position_catching();
                    Vector3d tang = catching_point - this->array_Of_Microtubules[ microtubule ].getPoint( 0 );

                    Vector3d tmp_position = this->array_Of_Microtubules[ microtubule ].getPoint( 0 ) + tang / 2.0;
                    parametr = parametr + 0.05;
                    fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  tmp_position( 0 ) * 1e6 , tmp_position( 1 ) * 1e6 , tmp_position( 2 ) * 1e6  );

                    parametr = parametr + 0.05;
                    fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  catching_point( 0 ) * 1e6 , catching_point( 1 ) * 1e6 , catching_point( 2 ) * 1e6  );
                    fprintf(out,"%u\n", this->array_Of_Microtubules[ microtubule ].get_dynein_index() );
                    fprintf(out,"END\n" );

                }

         }
         else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 2 )
         {
            Vector3d center = this->MTOC.get_center();
            double parametr = 0.0;
            fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  center( 0 ) * 1e6 , center( 1 ) * 1e6 , center( 2 ) * 1e6 );
            Vector3d catching_point = this->array_Of_Microtubules[ microtubule ].get_IS_position_catching();
            Vector3d tang = catching_point - center;

            Vector3d tmp_position = center + tang / 3.0;
            parametr = parametr + 0.05;
            fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  tmp_position( 0 ) * 1e6 , tmp_position( 1 ) * 1e6 , tmp_position( 2 ) * 1e6  );

            Vector3d tmp_position_2 = center + tang * 2.0 / 3.0;
            parametr = parametr + 0.05;
            fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  tmp_position_2( 0 ) * 1e6 , tmp_position_2( 1 ) * 1e6 , tmp_position_2( 2 ) * 1e6  );


            parametr = parametr + 0.05;
            fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  catching_point( 0 ) * 1e6 , catching_point( 1 ) * 1e6 , catching_point( 2 ) * 1e6  );
            fprintf(out,"%u\n", this->array_Of_Microtubules[ microtubule ].get_dynein_index() );
            fprintf(out,"END\n" );


         }
         else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9 )
         {

		    Vector3d tangent_last = this->array_Of_Microtubules[ microtubule ].getTangent2( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 );
		    if( tangent_last.norm() > 4 *  sim_of_Cell::resting_distance )
		    {
		    
		    }
		    else
		    {
                    Vector3d center = this->MTOC.get_center();
                    double parametr = 0.0;
                    fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  center( 0 )  * 1e6, center( 1 )  * 1e6, center( 2 )   * 1e6);

                    for( unsigned int point = 0 ; point < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() ; point ++  )
                    {
                            parametr = parametr + 0.05;
                            Vector3d coordinate_point = this->array_Of_Microtubules[ microtubule ].getPoint( point );
                            fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  coordinate_point( 0 ) * 1e6 , coordinate_point( 1 ) * 1e6 , coordinate_point( 2 ) * 1e6  );
                    }

                    Vector3d point_last = this->array_Of_Microtubules[ microtubule ].getPoint( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 );
                    Vector3d ill_point = point_last + tangent_last / 10.0;
                    parametr = parametr + 0.05;
                    fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  ill_point( 0 ) * 1e6 , ill_point( 1 ) * 1e6 , ill_point( 2 ) * 1e6  );
                    fprintf(out,"%u\n", this->array_Of_Microtubules[ microtubule ].get_dynein_index() );
                    fprintf(out,"END\n" );
                    }

         }

         else if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
         {
                    Vector3d tangent_last;
                    try
                    {
                        tangent_last = this->array_Of_Microtubules[ microtubule ].getTangent2( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 );
                    }
                    catch( int e )
                    {
                        cout<<"exception"<<endl;
                        cout<<"void void Cell::print_Cell()"<<endl;
                        throw("");
                    }
		    if( tangent_last.norm() > 4 *  sim_of_Cell::resting_distance )
		    {

		    }
		    else
		    {

                    Vector3d center = this->MTOC.get_center();
                    double parametr = 0.0;
                    fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  center( 0 )  * 1e6, center( 1 )  * 1e6, center( 2 )   * 1e6);

                    for( unsigned int point = 0 ; point < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() ; point ++  )
                    {
                        parametr = parametr + 0.05;
                            Vector3d coordinate_point = this->array_Of_Microtubules[ microtubule ].getPoint( point );
                        fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  coordinate_point( 0 ) * 1e6 , coordinate_point( 1 ) * 1e6 , coordinate_point( 2 ) * 1e6  );
                    }
                    Vector3d tangent_last;
                    try
                    {
                        tangent_last = this->array_Of_Microtubules[ microtubule ].getTangent2( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 );
                    }
                    catch( int e )
                    {
                        cout<<"exception"<<endl;
                        cout<<"void void Cell::print_Cell()"<<endl;
                        throw("");
                    }

                    Vector3d point_last = this->array_Of_Microtubules[ microtubule ].getPoint( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 );
                    Vector3d ill_point = point_last + tangent_last / 10.0;
                    parametr = parametr + 0.05;
                    fprintf(out,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  ill_point( 0 ) * 1e6 , ill_point( 1 ) * 1e6 , ill_point( 2 ) * 1e6  );
                    fprintf(out,"%u\n", this->array_Of_Microtubules[ microtubule ].get_dynein_index() );
                    fprintf(out,"END\n" );
		    }
         }
         else
	{
		cout<<endl;
		cout<<endl;
		cout<<endl;
		cout<<endl;
	}

	  }
	  fclose( out );



      //MTOC printing
	  char name_of_MTOC_file [55];
	  sprintf ( name_of_MTOC_file , "picturesVideos/textFiles/MTOC_%d.txt", index );
	  FILE *out_MTOC;
	  out_MTOC = fopen( name_of_MTOC_file , "w");
          Vector3d MTOC_point = this->MTOC.get_center();
          fprintf( out_MTOC , "%10.5f,%10.5lf,%10.5lf\n",  MTOC_point( 0 ) * 1e6 , MTOC_point( 1 ) * 1e6 , MTOC_point( 2 ) * 1e6 );
          fprintf( out_MTOC , "END\n" );


	  for( unsigned int poly_number = 1 ; poly_number <= this->MTOC.get_number_of_points() ; poly_number ++ )
	  {
                  Vector3d MTOC_point =   this->MTOC.get_point( poly_number );
		  fprintf( out_MTOC , "%10.5f,%10.5lf,%10.5lf\n",  MTOC_point( 0 ) * 1e6 , MTOC_point( 1 ) * 1e6 , MTOC_point( 2 ) * 1e6 );
		  fprintf( out_MTOC , "END\n" );
	  }
	  fclose( out_MTOC );


	  char name_of_IS_cathing_file [55];
	  sprintf ( name_of_IS_cathing_file , "picturesVideos/textFiles/IS_catching_%d.txt", index );
	  FILE *out_IS_catching;
	  out_IS_catching = fopen( name_of_IS_cathing_file , "w");

	  for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
	  {
            Vector3d IS_catching = this->array_Of_Microtubules[ microtubule ].get_IS_position_catching();
            if( IS_catching( 0 ) != 666.0 )
            {
                fprintf( out_IS_catching , "%10.5f,%10.5lf,%10.5lf\n",  IS_catching( 0 ) * 1e6 , IS_catching( 1 ) * 1e6 , IS_catching( 2 ) * 1e6 );
                fprintf( out_IS_catching , "END\n" );
            }
          }
	  fclose( out_IS_catching );


	  //!!!!!!!!!!!!!!!!!!!! IS capture shrinkage
	  //CAREFULL: povray enables cylinder scatch from two points and radius

	  char IS_capture_shrinkage [55];
	  sprintf ( IS_capture_shrinkage , "picturesVideos/textFiles/IS_capture_shrinkage_%d.txt", index );
	  FILE *out_IS_capture_shrinkage;
	  out_IS_capture_shrinkage = fopen( IS_capture_shrinkage , "w");
	  Vector3d IS_capture_shrinkage_center = this->IS_Capture_Shrinkage.get_center_of_IS_front();
	  double radius_IS_Capture_Shrinkage = this->IS_Capture_Shrinkage.get_radius_of_IS();
	  Vector3d IS_capture_shrinkage_rear_point = this->IS_Capture_Shrinkage.get_center_of_IS_rear();



	  fprintf( out_IS_capture_shrinkage , "%10.5f,%10.5lf,%10.5lf\n",  IS_capture_shrinkage_center( 0 ) * 1e6 , IS_capture_shrinkage_center( 1 ) * 1e6 , IS_capture_shrinkage_center( 2 ) * 1e6 );
	  fprintf( out_IS_capture_shrinkage , "END\n" );
	  fprintf( out_IS_capture_shrinkage , "%10.5f,%10.5lf,%10.5lf\n",  IS_capture_shrinkage_rear_point( 0 ) * 1e6 , IS_capture_shrinkage_rear_point( 1 ) * 1e6 , IS_capture_shrinkage_rear_point( 2 )  * 1e6 );
	  fprintf( out_IS_capture_shrinkage , "END\n" );
	  fprintf( out_IS_capture_shrinkage , "%10.5f\n",  radius_IS_Capture_Shrinkage * 1e6 );
      fclose( out_IS_capture_shrinkage );


          char IS_cortical_sl[55];
	  sprintf ( IS_cortical_sl , "picturesVideos/textFiles/IS_cortical_sl_%d.txt", index );
	  FILE *out_IS_cortical_sl;
	  out_IS_cortical_sl = fopen( IS_cortical_sl , "w");

	  Vector3d front_center1 = this->IS_Cortic_Sl.get_center_front_of_IS();
	  Vector3d rear_center1 = this->IS_Cortic_Sl.get_center_rear_of_IS();


	  Vector3d point_on_plane = this->IS_Cortic_Sl.get_point_on_plane();
	  Vector3d axis = this->IS_Cortic_Sl.get_axis();
          Vector3d tmp_rear_point = point_on_plane + ( - 1.0 ) * axis * 2e-6;
	  double radius = IS_Cortical_Sl_parameter::radius;
          fprintf( out_IS_cortical_sl , "%10.5f,%10.5lf,%10.5lf\n",  front_center1( 0 ) * 1e6 , front_center1( 1 ) * 1e6 , front_center1( 2 ) * 1e6 );
	  fprintf( out_IS_cortical_sl , "END\n" );
	  fprintf( out_IS_cortical_sl , "%10.5f,%10.5lf,%10.5lf\n",  rear_center1( 0 ) * 1e6 , rear_center1( 1 ) * 1e6 , rear_center1( 2 )  * 1e6 );
	  fprintf( out_IS_cortical_sl , "END\n" );
	  fprintf( out_IS_cortical_sl , "%10.5f\n",  radius * 1e6 );
	  fclose( out_IS_cortical_sl );






      char MTOC_center [55];
      sprintf ( MTOC_center , "picturesVideos/textFiles/MTOC_center_%d.txt", index );
      FILE *out_MTOC_center;
      out_MTOC_center = fopen( MTOC_center , "w");
      Vector3d center_of_MTOC = this->MTOC.get_center();
      fprintf( out_MTOC_center , "%10.5f,%10.5lf,%10.5lf\n",  center_of_MTOC( 0 ) * 1e6 , center_of_MTOC( 1 ) * 1e6 , center_of_MTOC( 2 ) * 1e6 );
      fclose( out_MTOC_center );


           char Dynein_surface_randomly_distributed [155];
           sprintf ( Dynein_surface_randomly_distributed , "picturesVideos/textFiles/Dynein_surface_randomly_distributed_%d.txt", index );
           FILE *out_Dynein_surface_randomly_distributed;
           out_Dynein_surface_randomly_distributed = fopen( Dynein_surface_randomly_distributed , "w");
           std::map< unsigned int , std::vector<Vector3d> > mapa = this->density_surface_dynein.get_map();


           for(  std::map< unsigned int , std::vector<Vector3d> >::iterator it=mapa.begin(); it != mapa.end(); ++it)
           {

                std::vector<Vector3d> vector_tmp = it->second;
                for( unsigned int i = 0 ; i < vector_tmp.size() ; i ++ )
                {
                    Vector3d vect = vector_tmp.at( i );
                    fprintf( out_Dynein_surface_randomly_distributed , "%10.5f,%10.5lf,%10.5lf\n",  vect( 0 ) * 1e6 , vect( 1 ) * 1e6 , vect( 2 ) * 1e6 );
                    fprintf( out_Dynein_surface_randomly_distributed , "END\n" );
                }
           }

           fclose( out_Dynein_surface_randomly_distributed );

           char Dynein_IS_capture[155];
           sprintf ( Dynein_IS_capture , "picturesVideos/textFiles/Dynein_IS_capture_%d.txt", index );
           FILE *out_Dynein_IS_capture;
           out_Dynein_IS_capture = fopen( Dynein_IS_capture , "w");
           std::map< unsigned int , std::vector<Vector3d> > mapa_capt = this->Capture_Shrinkage_dynein.get_map();
           for(  std::map< unsigned int , std::vector<Vector3d> >::iterator it=mapa_capt.begin(); it != mapa_capt.end(); ++it)
           {

                std::vector<Vector3d> vector_tmp = it->second;
                for( unsigned int i = 0 ; i < vector_tmp.size() ; i ++ )
                {
                    Vector3d vect = vector_tmp.at( i );
                    fprintf( out_Dynein_IS_capture , "%10.5f,%10.5lf,%10.5lf\n",  vect( 0 ) * 1e6 , vect( 1 ) * 1e6 , vect( 2 ) * 1e6 );
                    fprintf( out_Dynein_IS_capture , "END\n" );
                }
           }

           fclose( out_Dynein_IS_capture );



      char dynein_abscissa [155];
      sprintf ( dynein_abscissa , "picturesVideos/textFiles/dynein_abscissa_%d.txt", index );
      FILE *out_dynein_abscissa;
      out_dynein_abscissa = fopen( dynein_abscissa , "w");
      for( unsigned int micro = 0 ; micro < this->number_of_microtubules ; micro ++ )
      {
          if( ( this->array_Of_Microtubules[ micro ].get_dynein_index() == 9 ) || ( this->array_Of_Microtubules[ micro ].get_dynein_index() == 20 ) )
          {
                std::vector< std::pair < Vector3d ,double  > > vector_dynein_abscissa = this->array_Of_Microtubules[ micro ].get_Dynein_abscissa();
                if( vector_dynein_abscissa.size() < 1 )
                {
                    continue;
                }
                else
                {
                      for( unsigned int point = 0 ; point < vector_dynein_abscissa.size() ; point ++ )
                      {
                          std::pair < Vector3d ,double  > pair_tmp = vector_dynein_abscissa.at( point );
                          Vector3d position = std::get< 0 >( pair_tmp );
                          fprintf(out_dynein_abscissa,"%10.5f,%10.5lf,%10.5lf\n",position( 0 ) * 1e6 , position( 1 ) * 1e6 , position( 2 ) * 1e6 );
                          fprintf( out_dynein_abscissa , "END\n" );
                      }
                }
          }
      }
      fclose( out_dynein_abscissa );


      char dynein_abscissa_attachment [155];
      sprintf ( dynein_abscissa_attachment , "picturesVideos/textFiles/dynein_abscissa_attachment_%d.txt", index );
      FILE *out_dynein_abscissa_attachment;
      out_dynein_abscissa_attachment = fopen( dynein_abscissa_attachment , "w");
      for( unsigned int micro = 0 ; micro < this->number_of_microtubules ; micro ++ )
      {
          if( this->array_Of_Microtubules[ micro ].get_dynein_index() == 9  )
          {
                std::vector< std::pair < Vector3d ,double  > > vector_dynein_abscissa = this->array_Of_Microtubules[ micro ].get_Dynein_abscissa();
                for( unsigned int dynein = 0 ; dynein < vector_dynein_abscissa.size() ; dynein ++ )
                {
                    std::pair < Vector3d ,double  > one_pair = vector_dynein_abscissa.at( dynein );
                    double abscissa =   std::get<1>( one_pair );
                    unsigned int lower_bead_index =  this->array_Of_Microtubules[ micro ].get_index_according_to_abscissa( abscissa );

                    Vector3d lower_point = this->array_Of_Microtubules[ micro ].getPoint( lower_bead_index );
                    Vector3d tangent = this->array_Of_Microtubules[ micro ].getTangent2( lower_bead_index );


                    double abscissa_minus_lower_bead = abscissa - ( ( double ) lower_bead_index ) * this->array_Of_Microtubules[ micro ].getRestDist();
                    Vector3d p_of_att = lower_point + ( abscissa_minus_lower_bead / tangent.norm()  ) * tangent; // this->getRestDist()
                    fprintf(out_dynein_abscissa_attachment,"%10.5f,%10.5lf,%10.5lf\n",p_of_att( 0 )*1e6,p_of_att( 1 )*1e6,p_of_att( 2 )*1e6 );
                    fprintf( out_dynein_abscissa_attachment , "END\n" );

                }
          }

      }
      fclose( out_dynein_abscissa_attachment );


	  char MTOC_SPLINES_1[55];
	  sprintf ( name_of_MTOC_file , "picturesVideos/textFiles/MTOCSpline_%d.txt", index );
	  FILE *out_MTOC_spline_1;
	  out_MTOC_spline_1 = fopen( name_of_MTOC_file , "w");


	  int number_of_segments = 10;
	  for( unsigned int poly_number = 1 ; poly_number <= this->MTOC.get_number_of_points() ; poly_number ++ )
	  {
                  Vector3d MTOC_point =   this->MTOC.get_point( poly_number );
		  Vector3d tangent_to_center = center_of_MTOC - MTOC_point;
		  Vector3d mini_tangent = tangent_to_center / (double)number_of_segments;
	          double parametr = 0.0;
		  for( int index = 0 ; index <= number_of_segments  ; index ++ )
		  {
			Vector3d tmp_point = MTOC_point + ( double ) index * mini_tangent;
			parametr = parametr + 0.05;
			fprintf(out_MTOC_spline_1,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  tmp_point( 0 ) * 1e6 , tmp_point( 1 ) * 1e6 , tmp_point( 2 ) * 1e6  );

		  }
                  fprintf(out_MTOC_spline_1,"%u\n", 1 );
                  fprintf(out_MTOC_spline_1,"END\n" );

	  }


	  for( unsigned int poly_number = 1 ; poly_number <= this->MTOC.get_number_of_points() ; poly_number ++ )
	  {
                  Vector3d MTOC_point =   this->MTOC.get_point( poly_number );
		  Vector3d next_MTOC_point( 0.0 , 0.0 , 0.0 );
		  if( poly_number < this->MTOC.get_number_of_points() )
		  {
			next_MTOC_point =  this->MTOC.get_point( poly_number + 1 );
		  }
		  else
		  {
			next_MTOC_point =  this->MTOC.get_point( 1 );
		  }

		  Vector3d tangent_to_center = next_MTOC_point - MTOC_point;
		  Vector3d mini_tangent = tangent_to_center / (double)number_of_segments;
	          double parametr = 0.0;
		  for( int index = 0 ; index <= number_of_segments  ; index ++ )
		  {
			Vector3d tmp_point = MTOC_point + ( double ) index * mini_tangent;
			parametr = parametr + 0.05;
			fprintf(out_MTOC_spline_1,"%10.5f <%10.5f,%10.5lf,%10.5lf>\n", parametr ,  tmp_point( 0 ) * 1e6 , tmp_point( 1 ) * 1e6 , tmp_point( 2 ) * 1e6  );

		  }
                  fprintf(out_MTOC_spline_1,"%u\n", 0 );
                  fprintf(out_MTOC_spline_1,"END\n" );



	  }
	  fclose( out_MTOC_spline_1 );




}

void Cell::print_Cell_parametres( )
{
	FILE *out;
	out = fopen( "picturesVideos/cell_shapes/Cell_Parametres.txt" , "w");
	fprintf(out,"%15.12f , %15.12f\n", this->a_axis * 1e6 , this->b_axis * 1e6 );
	fclose( out );

	FILE *out2;
	out2 = fopen( "picturesVideos/cell_shapes/Nucleus_Parametres.txt" , "w");
	fprintf(out2,"%15.12f , %15.12f\n",  this->nucleus.get_A_Axis() * 1e6 , this->nucleus.get_B_Axis() * 1e6  );
	fprintf(out2,"END\n" );
	Vector3d nucleus_center = this->nucleus.get_center();
	fprintf(out2,"%15.12f , %15.12f , %15.12f", nucleus_center( 0 ) * 1e6 , nucleus_center( 1 ) * 1e6 , nucleus_center( 2 ) * 1e6 );
	fclose( out2 );
}


Vector3d Cell::force_position_cell_wall( Vector3d position )
{
    return this->force_position_cell_wall_two_elipsoid( position ); // force_position_cell_wall_simple_elipsoid
}



Vector3d Cell::force_position_cell_wall_simple_elipsoid( Vector3d position )
{
	Vector3d force( 0.0 , 0.0 , 0.0 );
	double confirm_ellipsoid_value = ( position( 0 ) * position( 0 ) + position( 1 ) * position( 1 ) ) / ( this->a_axis * this->a_axis );
	confirm_ellipsoid_value = confirm_ellipsoid_value + ( position( 2 ) * position( 2 ) ) / ( this->b_axis * this->b_axis );

	if( confirm_ellipsoid_value < 1.0 )
	{

	}
	else
	{
		double constant_of_intersection = sqrt( 1.0 / confirm_ellipsoid_value );
		Vector3d point_of_intersection = constant_of_intersection * position;

		double distance_behind_wall = ( position - point_of_intersection ).norm();
		double absolut_value_of_the_force = Cell_parametres::wall_cell_k1 * ( exp( Cell_parametres::wall_cell_k2 * distance_behind_wall ) - 1.0 );
		Vector3d orientation = ( -1.0 ) * ( position - point_of_intersection ) / ( position - point_of_intersection ).norm();
		force = orientation * absolut_value_of_the_force;
	}
	return force;
}


Vector3d Cell::force_position_cell_wall_two_elipsoid( Vector3d position  )
{
    Vector3d force( 0.0 , 0.0 , 0.0 );

    double element = 0.0;

    double BB;
    double B_multiplicator;
    if( position( 2 )  < 0 )
    {
        BB = Cell_parametres::B_AXIS_Lower + element;
        B_multiplicator = ( Cell_parametres::B_AXIS_Lower + element  ) * ( Cell_parametres::B_AXIS_Lower + element  );
    }
    else
    {
        BB =  ( Cell_parametres::B_AXIS + element );
        B_multiplicator = ( Cell_parametres::B_AXIS + element ) * ( Cell_parametres::B_AXIS + element );
    }
    double numerator = ( Cell_parametres::A_AXIS + element ) * BB;
    double A_multiplicator = ( Cell_parametres::A_AXIS + element )  * ( Cell_parametres::A_AXIS + element ) ;

    double denominator = ( position( 0 ) * position( 0 ) + position( 1 ) * position( 1 ) ) * B_multiplicator;
    denominator = denominator + position( 2 ) * position( 2 ) * A_multiplicator;
    denominator = sqrt( denominator );
    double c_multiplicator = numerator / denominator;
    if( c_multiplicator >= 1 )
    {
        return force;
    }
    else
    {
        Vector3d point_of_intersection = c_multiplicator * position;
        double distance_behind_wall = ( position - point_of_intersection ).norm();
        double absolut_value_of_the_force = Cell_parametres::wall_cell_k1 * ( exp( Cell_parametres::wall_cell_k2 * distance_behind_wall ) - 1.0 );
        Vector3d orientation = ( -1.0 ) * ( position - point_of_intersection ) / ( position - point_of_intersection ).norm();
        force = orientation * absolut_value_of_the_force;
        return force;
    }


}




Vector3d Cell::force_position_cell_wall_sphere( Vector3d position  )
{
    Vector3d force( 0.0 , 0.0 , 0.0 );
    double abs_val = position.norm();
    if( abs_val < this->a_axis )
    {

    }
    else
    {
        double distance = abs_val - this->a_axis;
        double absolut_value_of_the_force = Cell_parametres::wall_cell_k1 * ( exp( Cell_parametres::wall_cell_k2 * distance ) - 1.0 );
        Vector3d orientation = -1.0 * position / position.norm();
        force = orientation * absolut_value_of_the_force;
    }

    return force;
}



Vector3d Cell::force_position_cell_wall_simple_elipsoid_MTOC( Vector3d position )
{
	Vector3d force( 0.0 , 0.0 , 0.0 );
	double confirm_ellipsoid_value = ( position( 0 ) * position( 0 ) + position( 1 ) * position( 1 ) ) / ( this->a_axis * this->a_axis );
	confirm_ellipsoid_value = confirm_ellipsoid_value + ( position( 2 ) * position( 2 ) ) / ( this->b_axis * this->b_axis );

	if( confirm_ellipsoid_value < 1.0 )
	{
	}
	else
	{
		//first, I have to find the point of intersection
		double constant_of_intersection = sqrt( 1.0 / confirm_ellipsoid_value );
		Vector3d point_of_intersection = constant_of_intersection * position;

		double distance_behind_wall = ( position - point_of_intersection ).norm();

		double absolut_value_of_the_force = Cell_parametres::wall_cell_k1 * ( exp( Cell_parametres::wall_cell_k2 * distance_behind_wall ) - 1.0 );
		Vector3d orientation = ( -1.0 ) * ( position - point_of_intersection ) / ( position - point_of_intersection ).norm();
		force = orientation * absolut_value_of_the_force;
	}

    return force;
}





void Cell::MTOC_microtubule_force( Vector3d& center_bead_force , Vector3d& point_bead_force , unsigned int microtubule_number )
{
    if( microtubule_number >= this->number_of_microtubules )
    {
        cout<<"MTOC_microtubule_force( MatrixXd force , unsigned int microtubule_number )"<<endl;
        cout<<"microtubule_number >= this->number_of_microtubules"<<endl;
        cout<<"microtubule_number = "<<microtubule_number<<endl;
        cout<<"ERROR_ID Cell974463414665"<<endl;
        throw("");
    }


    Vector3d micro_first_point = this->array_Of_Microtubules[ microtubule_number ].getPoint( 0 );
    unsigned int MTOC_point_index = this->array_Of_Microtubules[ microtubule_number ].get_MTOC_point();

    Vector3d MTOC_center = this->MTOC.get_center();
    Vector3d MTOC_point = this->MTOC.get_point( MTOC_point_index );

    //force center first micro_bead
    double distance_center_micro = ( micro_first_point - MTOC_center ).norm();
    distance_center_micro = abs( distance_center_micro - this->MTOC.get_radius() );

    Vector3d orientation = ( MTOC_center - micro_first_point );
    orientation = orientation / orientation.norm();


    double center_bead_force_abs_value = distance_center_micro  * MTOCparam::MTOC2_center_micro_kappa;
    center_bead_force = center_bead_force_abs_value * orientation;

    //force MTOC_point bead
    double distance_point_micro = ( MTOC_point - micro_first_point ).norm();

    if( this->array_Of_Microtubules[ microtubule_number ].get_dynein_index() == 7 )
    {

    }

    Vector3d orientation_point_bead = MTOC_point - micro_first_point;
    orientation_point_bead = orientation_point_bead / orientation_point_bead.norm();
    double point_bead_force_abs_value = distance_point_micro * MTOCparam::MTOC2_point_micro_kappa;
    point_bead_force = point_bead_force_abs_value * orientation_point_bead;

}

void Cell::MTOC_microtubule_two_points_force( Vector3d& first_point_force , Vector3d& second_point_force , unsigned int microtubule_number )
{
    if( microtubule_number >= this->number_of_microtubules )
    {
        cout<<"MTOC_microtubule_two_points_force( MatrixXd force , unsigned int microtubule_number )"<<endl;
        cout<<"microtubule_number >= this->number_of_microtubules"<<endl;
        cout<<"microtubule_number = "<<microtubule_number<<endl;
        cout<<"ERROR_ID = 666719816168645"<<endl;
        throw("");
    }

    if( this->array_Of_Microtubules[ microtubule_number ].getNumberOfPoints() >= 2 )
    {

        unsigned int mtoc_point_index = this->array_Of_Microtubules[ microtubule_number ].get_MTOC_point();
        unsigned int opposite_mtoc_point_index = this->array_Of_Microtubules[ microtubule_number ].get_MTOC_opposite_point();
        Vector3d MTOC_opposite_point = this->MTOC.get_point( opposite_mtoc_point_index );
        Vector3d first_point_micro = this->array_Of_Microtubules[ microtubule_number ].getPoint(0);

        double distance_1 = ( MTOC_opposite_point - first_point_micro ).norm();
        if( distance_1 == 0 )
        {

        }
        else
        {
            Vector3d orientation_1 = ( MTOC_opposite_point - first_point_micro );
            orientation_1 = orientation_1 / orientation_1.norm();
            double force_1_abs = distance_1 * MTOCparam::MTOC2_point_micro_kappa;
            Vector3d force_1 = force_1_abs * orientation_1;
            first_point_force = force_1;
        }



        Vector3d MTOC_point = this->MTOC.get_point( mtoc_point_index );
        Vector3d second_point_micro = this->array_Of_Microtubules[ microtubule_number ].getPoint( 1 );

        double distance_2 = ( MTOC_point - second_point_micro ).norm();
        if( distance_2 == 0 )
        {

        }
        else
        {
            Vector3d orientation_2 = ( MTOC_point - second_point_micro );
            orientation_2 = orientation_2 / orientation_2.norm();
            double force_2_abs = distance_2 * MTOCparam::MTOC2_point_micro_kappa;
            Vector3d force_2 = force_2_abs * orientation_2;
            second_point_force = force_2;
        }


    }
    else if( this->array_Of_Microtubules[ microtubule_number ].getNumberOfPoints() == 1 )
    {
        unsigned int opposite_mtoc_point_index = this->array_Of_Microtubules[ microtubule_number ].get_MTOC_opposite_point();
        Vector3d MTOC_opposite_point = this->MTOC.get_point( opposite_mtoc_point_index );
        Vector3d first_point_micro = this->array_Of_Microtubules[ microtubule_number ].getPoint(0);

        double distance_1 = ( MTOC_opposite_point - first_point_micro ).norm();
        if( distance_1 == 0 )
        {

        }
        else
        {
            Vector3d orientation_1 = ( MTOC_opposite_point - first_point_micro );
            orientation_1 = orientation_1 / orientation_1.norm();
            double force_1_abs = distance_1 * MTOCparam::MTOC2_point_micro_kappa;
            Vector3d force_1 = force_1_abs * orientation_1;
            first_point_force = force_1;

        }
    }



}



void Cell::MTOC_microtubule_two_points_force_with_bending( Vector3d& first_point_force , Vector3d& second_point_force , Vector3d& bending_opposite_MTOC_point , Vector3d& bending_MTOC_point , unsigned int microtubule_number )
{

    if( microtubule_number >= this->number_of_microtubules )
    {
        cout<<"MTOC_microtubule_two_points_force( MatrixXd force , unsigned int microtubule_number )"<<endl;
        cout<<"microtubule_number >= this->number_of_microtubules"<<endl;
        cout<<"microtubule_number = "<<microtubule_number<<endl;
        cout<<"ERROR_ID = 666719816168645"<<endl;
        throw("");
    }

    if( this->array_Of_Microtubules[ microtubule_number ].getNumberOfPoints() >= 2 )
    {

        unsigned int mtoc_point_index = this->array_Of_Microtubules[ microtubule_number ].get_MTOC_point();
        unsigned int opposite_mtoc_point_index = this->array_Of_Microtubules[ microtubule_number ].get_MTOC_opposite_point();
        Vector3d MTOC_opposite_point = this->MTOC.get_point( opposite_mtoc_point_index );
        Vector3d first_point_micro = this->array_Of_Microtubules[ microtubule_number ].getPoint(0);

        double distance_1 = ( MTOC_opposite_point - first_point_micro ).norm();
        if( distance_1 == 0 )
        {

        }
        else
        {
            Vector3d orientation_1 = ( MTOC_opposite_point - first_point_micro );
            orientation_1 = orientation_1 / orientation_1.norm();
            double force_1_abs = distance_1 * MTOCparam::MTOC2_point_micro_kappa;
            Vector3d force_1 = force_1_abs * orientation_1;
            first_point_force = force_1;
        }



        Vector3d MTOC_point = this->MTOC.get_point( mtoc_point_index );
        Vector3d second_point_micro = this->array_Of_Microtubules[ microtubule_number ].getPoint( 1 );

        double distance_2 = ( MTOC_point - second_point_micro ).norm();
        if( distance_2 == 0 )
        {

        }
        else
        {
            Vector3d orientation_2 = ( MTOC_point - second_point_micro );
            orientation_2 = orientation_2 / orientation_2.norm();
            double force_2_abs = distance_2 * MTOCparam::MTOC2_point_micro_kappa;
            Vector3d force_2 = force_2_abs * orientation_2;
            second_point_force = force_2;
        }

	double stiffness_coeffcient = 1.0;

        Vector3d bending_force_one_bead_micro( 0.0 , 0.0 , 0.0 );

        Vector3d mic_first_segment = this->array_Of_Microtubules[ microtubule_number ].getTangent2( 0 );
        Vector3d MTOC_segment = MTOC_point - MTOC_opposite_point;

        Vector3d first = mic_first_segment;
	Vector3d second = MTOC_segment;


        double magThisVec = first.norm();
	first = first * ( 1.0 / magThisVec );

	double magNextVec = second.norm();
	second = second * ( 1.0 / magNextVec );		//ATTENTION - here, it will be divided by magNextVec - the last tangent

	bending_force_one_bead_micro( 0 ) = ( first( 0 ) - second( 0 ) * first.dot( second ) ) / magNextVec;
	bending_force_one_bead_micro( 1 ) = ( first( 1 ) - second( 1 ) * first.dot( second ) ) / magNextVec;
	bending_force_one_bead_micro( 2 ) = ( first( 2 ) - second( 2 ) * first.dot( second ) ) / magNextVec;
        bending_force_one_bead_micro = ( -1.0 ) * sim_of_Cell::k_bending_analytical / magThisVec * bending_force_one_bead_micro * stiffness_coeffcient;
        second_point_force = second_point_force + bending_force_one_bead_micro;


        Vector3d bending_force_MTOC_point( 0.0 , 0.0 , 0.0 );
        bending_force_MTOC_point( 0 ) = ( - second( 0 ) + first( 0 ) * first.dot( second ) ) / magThisVec;
        bending_force_MTOC_point( 1 ) = ( - second( 1 ) + first( 1 ) * first.dot( second ) ) / magThisVec;
	bending_force_MTOC_point( 2 ) = ( - second( 2 ) + first( 2 ) * first.dot( second ) ) / magThisVec;
        bending_force_MTOC_point = ( -1.0 ) * sim_of_Cell::k_bending_analytical / magNextVec * bending_force_MTOC_point * stiffness_coeffcient;
        bending_MTOC_point = bending_force_MTOC_point;

        Vector3d central_point_force = ( -1.0 ) * bending_force_one_bead_micro + ( -1.0 ) * bending_force_MTOC_point;
        first_point_force = first_point_force + central_point_force;

        bending_opposite_MTOC_point = central_point_force;


    }
    else if( this->array_Of_Microtubules[ microtubule_number ].getNumberOfPoints() == 1 )
    {
        unsigned int opposite_mtoc_point_index = this->array_Of_Microtubules[ microtubule_number ].get_MTOC_opposite_point();
        Vector3d MTOC_opposite_point = this->MTOC.get_point( opposite_mtoc_point_index );
        Vector3d first_point_micro = this->array_Of_Microtubules[ microtubule_number ].getPoint(0);

        double distance_1 = ( MTOC_opposite_point - first_point_micro ).norm();
        if( distance_1 == 0 )
        {

        }
        else
        {
            Vector3d orientation_1 = ( MTOC_opposite_point - first_point_micro );
            orientation_1 = orientation_1 / orientation_1.norm();
            double force_1_abs = distance_1 * MTOCparam::MTOC2_point_micro_kappa;
            Vector3d force_1 = force_1_abs * orientation_1;
            first_point_force = force_1;

        }
    }




}

Vector3d Cell::MTOC_microtubule_two_points(  unsigned int microtubule_number )
{
    if( microtubule_number >= this->number_of_microtubules )
    {
        cout<<"MTOC_microtubule_two_points_force_with_bending( MatrixXd force , unsigned int microtubule_number )"<<endl;
        cout<<"microtubule_number >= this->number_of_microtubules"<<endl;
        cout<<"microtubule_number = "<<microtubule_number<<endl;
        cout<<"ERROR_ID = 6667198641546865645"<<endl;
        throw("");
    }






    unsigned int mtoc_point_index = this->array_Of_Microtubules[ microtubule_number ].get_MTOC_point();
    Vector3d MTOC_point = this->MTOC.get_point( mtoc_point_index );
    Vector3d second_point_micro = this->array_Of_Microtubules[ microtubule_number ].getPoint( 1 );

    double distance_2 = ( MTOC_point - second_point_micro ).norm();
    if( distance_2 == 0 )
    {
	Vector3d force( 0.0 , 0.0 , 0.0 );
	return force;
    }
    else
    {
        Vector3d orientation_2 = ( MTOC_point - second_point_micro );
        orientation_2 = orientation_2 / orientation_2.norm();
        double force_2_abs = distance_2 * MTOCparam::MTOC2_point_micro_kappa;
        Vector3d force_2 = force_2_abs * orientation_2;
        return force_2;
    }



}





void Cell::force_whole_microtubule_nucleus( MatrixXd &force , unsigned int index )
{
	if( index >= this->number_of_microtubules )
	{
		cout<<"index >= this->number_of_microtubules in Cell::force_whole_microtubule_nucleus( MatrixXd &force , unsigned int index )"<<endl;
		cout<<"index = "<<index<<endl;
		throw("");
	}
	Microtubule temporary_microtubule = this->getMicrotubule( index );

	if( force.rows() != 3 * temporary_microtubule.getNumberOfPoints()  )
	{
		cout<<"force.rows() != 3 * temporary_microtubule.getNumberOfPoints() in force_whole_microtubule_nucleus"<<endl;
		cout<<"index = "<<index<<endl;
		throw("");
	}

	for( unsigned int bead_index = 0 ; bead_index < temporary_microtubule.getNumberOfPoints() ; bead_index ++ )
	{
		Vector3d position = temporary_microtubule.getPoint( bead_index );
		Vector3d force_one_bead = this->nucleus.force_position_nucleus_wall( position );

		for( unsigned int index_dimension = 0 ; index_dimension < 3 ; index_dimension ++ )
		{
			force( 3 * bead_index + index_dimension , 0 ) = force_one_bead( index_dimension );
		}
	}

}


Vector3d Cell::MTOC_point_wall_interaction( Vector3d position )
{
     return force_position_cell_wall_two_elipsoid( position );//force_position_cell_wall_simple_elipsoid_MTOC force_position_cell_wall_sphere
}




void Cell::force_MTOC_cell_wall( MatrixXd &force_MTOC )
{

    if( force_MTOC.rows() != 3 * ( this->MTOC.get_number_of_points() + 1 ) )
    {
        cout<<"force_MTOC.rows() != this->MTOC.get_number_of_points() * 3"<<endl;
        cout<<"Cell::force_MTOC_cell_wall( MatrixXd &force_MTOC )"<<endl;
        cout<<"ERROR_ID Cell 4646116911546945"<<endl;
        throw("");
    }


    Vector3d center = this->MTOC.get_center();
    Vector3d force = this->MTOC_point_wall_interaction( center );

    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        force_MTOC( dimension , 0 ) = force( dimension );
    }

}


void Cell::MTOC_nucleus_interaction( MatrixXd &force_on_MTOC )
{

    if( ( force_on_MTOC.rows() != ( this->MTOC.get_number_of_points() + 1 ) * 3 ) || ( force_on_MTOC.cols() != 1 ) )
    {
        cout<<"( force_on_MTOC.rows() != ( this->MTOC.get_number_of_points() + 1 ) * 3 ) || ( force_on_MTOC.cols() != 1 )"<<endl;
        cout<<" Cell::MTOC_nucleus_interaction( MatrixXd &force_on_MTOC ) "<<endl;
        cout<<"CELL ERROR_ID = 68464186441684"<<endl;
        throw("");
    }

    Vector3d center_position = this->MTOC.get_center();
    Vector3d force_center = this->nucleus.force_position_nucleus_wall( center_position );
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        force_on_MTOC( dimension , 0 ) = force_center( dimension , 0 );
    }

    for( unsigned int bod = 1 ; bod <= this->MTOC.get_number_of_points() ; bod ++ )
    {
        Vector3d position_point = this->MTOC.get_point( bod );
        Vector3d force_on_point = this->nucleus.force_position_nucleus_wall( position_point );
        for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
        {
            force_on_MTOC( 3 * bod + dimension , 0 ) = force_on_point( dimension , 0 );
        }
    }



}





void Cell::IS_capture_shrinkage_shrinkage_of_microtubule()
{
    for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ ) //
	{
		unsigned int dynein_index = this->array_Of_Microtubules[ microtubule ].get_dynein_index();
		if( ( dynein_index > 0 ) && ( dynein_index < 4 ) )
		{
			this->IS_micro_control_one_micro_capture_shrinkage_shrinkage( microtubule );
		}
    }

}






void Cell::IS_micro_control_one_micro_capture_shrinkage_shrinkage( unsigned int microtubule )
{
    if( ( this->array_Of_Microtubules[ microtubule ].get_dynein_index() < 1 ) && ( this->array_Of_Microtubules[ microtubule ].get_dynein_index() > 3 ) )
    {
        cout<<"( this->array_Of_Microtubules[ microtubule ].get_dynein_index()<1)&&(this->array_Of_Microtubules[ microtubule ].get_dynein_index() > 3 )"<<endl;
    }

    unsigned int dynein_index = this->array_Of_Microtubules[ microtubule ].get_dynein_index();

    if ( ( dynein_index == 1 ) && ( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() ) > 1 )
	{

		Vector3d bead_position = this->array_Of_Microtubules[ microtubule ].getPoint(  this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 );
		Vector3d catching_point = this->array_Of_Microtubules[ microtubule ].get_IS_position_catching();
        Vector3d front_center_IS = this->IS_Capture_Shrinkage.get_center_of_IS_front();
        Vector3d rear_center_IS = this->IS_Capture_Shrinkage.get_center_of_IS_rear();
        Vector3d axis_S = this->IS_Capture_Shrinkage.get_axis_of_IS();
        double vzdalenost_1 = distance_plane_point( axis_S , catching_point , bead_position );
        double vzdalenost_2 = distance_point_line( front_center_IS , rear_center_IS , bead_position );

        if( ( bead_position - catching_point ).norm() < IS_Capture_shrinkage_param::cut_distance )
		{
			this->array_Of_Microtubules[ microtubule ].Substract_one_Bead();

		}

	}


	else if ( ( dynein_index == 1 ) && ( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() ) == 1 )
	{
		Vector3d bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( 0 );
		Vector3d tangent_Center_IS_bead = ( bead_position - this->array_Of_Microtubules[ microtubule ].get_IS_position_catching() );


        Vector3d catching_point = this->array_Of_Microtubules[ microtubule ].get_IS_position_catching();
        Vector3d front_center_IS = this->IS_Capture_Shrinkage.get_center_of_IS_front();
        Vector3d rear_center_IS = this->IS_Capture_Shrinkage.get_center_of_IS_rear();
        Vector3d axis_S = this->IS_Capture_Shrinkage.get_axis_of_IS();
        double vzdalenost_1 = distance_plane_point( axis_S , catching_point , bead_position );
        double vzdalenost_2 = distance_point_line( front_center_IS , rear_center_IS , bead_position );

        if( ( bead_position - catching_point ).norm() < IS_Capture_shrinkage_param::cut_distance )
		{
            cout<<"destroy_Microtubule_capture_shrinkage() = "<<tangent_Center_IS_bead.norm()<<endl;
			this->array_Of_Microtubules[ microtubule ].destroy_Microtubule_capture_shrinkage();
        }

	}

}






bool Cell::check_IS_cortical_sl( Vector3d position )
{
    return this->IS_Cortic_Sl.check_ISCorticalSl_caught( position );
}





bool Cell::IS_cortical_sl_one_micro_control_CAUGHT( unsigned int microtubule )
{

    if( this->array_Of_Microtubules[ microtubule ].get_polygon_number() >= IS_Cortical_Sl_parameter::number_of_polygon )
    {
        return 0;
    }
    if( this->array_Of_Microtubules[ microtubule ].getID() >= IS_Cortical_Sl_parameter::number_of_mito )
    {
        return 0;
    }

    for( unsigned int i = 0 ; i < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() ; i ++ )
    {
        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 0 )
        {
            Vector3d bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( i );
            bool mark = this->check_IS_cortical_sl( bead_position );
            if( mark == true )
            {
                return true;
            }
        }
    }
    return false;
}








double Cell::get_covered_distance_by_IS( unsigned int microtubule )
{
    double distance_cumulative = 0;
    for( unsigned int bead_number = 0 ; bead_number < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 ; bead_number ++ )
    {
        Vector3d first_bead = this->array_Of_Microtubules[ microtubule ].getPoint( bead_number );
        Vector3d second_bead = this->array_Of_Microtubules[ microtubule ].getPoint( bead_number + 1 );
        double distance = this->IS_Cortic_Sl.covered_distance_of_segment( first_bead , second_bead );
        distance_cumulative = distance_cumulative + distance;
    }
    return distance_cumulative;
}


Vector3d Cell::dynein_force_one_bead( Vector3d position , Vector3d tangent_to_be_projected , double absolute_value_of_force )
{
    Vector3d projection = projection_of_tangent_on_plane( position , tangent_to_be_projected );
    Vector3d orientation = projection / projection.norm();
    Vector3d force = orientation * absolute_value_of_force;
    return force;
}







void Cell::Dynein_on_surface_force_2(unsigned int microtubule , MatrixXd& force_whole_microtubule )
{
     double covered_distance = this->get_covered_distance_by_IS( microtubule );
     double lenght_of_microtubule = ( double ) ( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 );
     lenght_of_microtubule = lenght_of_microtubule * this->array_Of_Microtubules[ microtubule ].getRestDist();
     double lenght_not_covered_by_IS = lenght_of_microtubule - covered_distance;
     double overal_force = lenght_not_covered_by_IS * IS_Dynein_Cell_surface::force_per_lenght;
     double force_on_one_bead = overal_force / ( (double) this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() );

     for( unsigned int point_number = 0 ; point_number < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 ; point_number ++ )
     {
         Vector3d tangent = this->array_Of_Microtubules[ microtubule ].getTangent2( point_number );
         Vector3d orientation = tangent / tangent.norm();
         Vector3d force = orientation * force_on_one_bead;
         for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
         {
            force_whole_microtubule( 3 * point_number + dimension ) = force( dimension );
         }
     }

     unsigned int number_of_points = this->array_Of_Microtubules[ microtubule ].getNumberOfPoints();
     Vector3d tangent = this->array_Of_Microtubules[ microtubule ].getTangent2( number_of_points - 2 );
     Vector3d orientation = tangent / tangent.norm();
     Vector3d force = orientation * force_on_one_bead;
     for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
     {
        force_whole_microtubule( 3 * ( number_of_points - 1 ) + dimension ) = force( dimension );
     }

}





void Cell::whole_microtubule_Surface_Dynein_force(unsigned int microtubule , MatrixXd& force_whole_microtubule )
{
     double covered_distance = this->get_covered_distance_by_IS( microtubule );
     double lenght_of_microtubule = ( double ) ( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 );
     lenght_of_microtubule = lenght_of_microtubule * this->array_Of_Microtubules[ microtubule ].getRestDist();
     double lenght_not_covered_by_IS = lenght_of_microtubule - covered_distance;
     double overal_force = lenght_not_covered_by_IS * IS_Dynein_Cell_surface::force_per_lenght;
     double force_on_one_bead = overal_force / ( (double) this->array_Of_Microtubules[ microtubule ].getNumberOfPoints()  );


    for( unsigned int i = 0; i < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints(); i ++ )
    {
        Vector3d tangent_to_be( 0.0 , 0.0 , 0.0 );
        if( i == 0)
        {
           tangent_to_be  = this->array_Of_Microtubules[ microtubule ].getTangent2( 0 );
        }
        else if( i == this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 )
        {
            unsigned int tangent_index = this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2;
            tangent_to_be  = this->array_Of_Microtubules[ microtubule ].getTangent2( tangent_index );
        }
        else
        {
            tangent_to_be = this->array_Of_Microtubules[ microtubule ].getPoint( i + 1 ) - this->array_Of_Microtubules[ microtubule ].getPoint( i - 1 );
            tangent_to_be = this->array_Of_Microtubules[ microtubule ].getTangent2( i );
        }

        Vector3d orientation = tangent_to_be;
        orientation = orientation / orientation.norm();
        Vector3d force = orientation * force_on_one_bead;
        for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++  )
        {
           force_whole_microtubule( 3 * i + dimension , 0 ) = force_whole_microtubule( 3 * i + dimension , 0 ) + force( dimension );
        }
    }
}


void Cell::set_density_surface_dynein( Surface& density )
{
    this->density_surface_dynein = density;
}



void Cell::all_microtubules_random_dynein_on_surface_catch_pair_abscissa()
{
unsigned int microtubule;
unsigned int chunk = 7;
#pragma omp parallel shared(  chunk ) private( microtubule  )
{
    for( microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
    {
        if( ( this->array_Of_Microtubules[ microtubule ].get_dynein_index() != 0 ) && ( this->array_Of_Microtubules[ microtubule ].get_dynein_index() != 9 ) )
        {
            continue;
        }
        if( this->array_Of_Microtubules[ microtubule ].get_polygon_number() >= IS_Dynein_Cell_surface::number_of_polygon_higher )
        {
            continue;
        }
        if( this->array_Of_Microtubules[ microtubule ].get_polygon_number() < IS_Dynein_Cell_surface::number_of_polygon_lower )
        {
            continue;
        }
        if( this->array_Of_Microtubules[ microtubule ].getID() >= IS_Dynein_Cell_surface::number_of_mito )
        {
            continue;
        }

        this->microtubule_random_dynein_on_surface_catch_pair_abscissa( microtubule );
    }
}
}


void Cell::microtubule_random_dynein_on_surface_catch_pair_abscissa( unsigned int microtubule )
{

    if( ( this->array_Of_Microtubules[  microtubule].get_dynein_index() != 0 ) && ( this->array_Of_Microtubules[  microtubule].get_dynein_index() != 9 ) )
    {
        cout<<"(this->array_Of_Microtubules[microtubule].get_dynein_index()!= 0)&&( this->array_Of_Microtubules[  microtubule].get_dynein_index() != 9 )"<<endl;
        cout<<"Cell ERROR_ID = 468168434875465"<<endl;
        throw("");
    }
    else
    {
        unsigned int number_of_points_in_micro = this->array_Of_Microtubules[ microtubule ].getNumberOfPoints();

        if( number_of_points_in_micro <= 1 )
        {
            return;
        }


        double resting_distance = this->array_Of_Microtubules[ microtubule ].getRestDist();
        unsigned int number_of_test_bead_in_segment = resting_distance / this->density_surface_dynein.get_X_width(); //this->density_surface_dynein.get_X_width()

        for( unsigned int i = 0 ; i < number_of_points_in_micro - 1 ; i ++ ) //NARAZNIK od 1 do ( number_of_points_in_micro - 1 ) / 2
        {

            Vector3d position_bead = this->array_Of_Microtubules[ microtubule ].getPoint( i );
            Vector3d tangent = this->array_Of_Microtubules[ microtubule ].getTangent2( i );
            Vector3d mini_tangent = tangent / ( double )( number_of_test_bead_in_segment + 1 );


            for( unsigned int tangent_counter = 0 ; tangent_counter <= number_of_test_bead_in_segment ; tangent_counter ++ )
            {
                Vector3d tmp_point_position = position_bead + ( ( double ) tangent_counter ) * mini_tangent;
                std::vector<Vector3d> catched_points = microtubule_point_attach_dynein( tmp_point_position );
                if( catched_points.size() > 0 )
                {

                    this->array_Of_Microtubules[ microtubule ].set_dynein_index( 9 );
                    for( unsigned int point_index = 0 ; point_index < catched_points.size() ; point_index ++ )
                    {
                        Vector3d catched_point = catched_points.at( point_index );
                        Vector3d closest_point_of_segment( 0.0 , 0.0 , 0.0 );
                        Vector3d a_vector = position_bead;
                        Vector3d b_vector = tangent;
                        Vector3d y_vector = catched_point;
                        distance_point_segment( a_vector , b_vector , y_vector , closest_point_of_segment );

                        double distance_bead_c_point = ( closest_point_of_segment - position_bead ).norm();
                        double abscissa = ( double ) i * this->array_Of_Microtubules[ microtubule ].getRestDist() + distance_bead_c_point;
                        std::pair < Vector3d , double > pair_point_abs( closest_point_of_segment , abscissa );
                        this->array_Of_Microtubules[ microtubule ].add_pair( pair_point_abs );

                    }
                }


            }

        }

    }
}


std::vector<Vector3d> Cell::microtubule_point_attach_dynein( Vector3d point_position  )
{
        std::uniform_real_distribution<> distribution{ 0 , 1 };
	unsigned int number_of_generator = omp_get_thread_num();


        std::vector< Vector3d > returned_points;
        double BB;
        double B_multiplicator;
        if( point_position( 2 ) < 0 )
        {
            BB = Cell_parametres::B_AXIS_Lower;
            B_multiplicator = Cell_parametres::B_AXIS_Lower * Cell_parametres::B_AXIS_Lower;
        }
        else
        {
            BB = Cell_parametres::B_AXIS;
            B_multiplicator = Cell_parametres::B_AXIS * Cell_parametres::B_AXIS;
        }
        double numerator = Cell_parametres::A_AXIS * BB;
        double A_multiplicator = Cell_parametres::A_AXIS * Cell_parametres::A_AXIS;

        double denominator = ( point_position( 0 ) * point_position( 0 ) + point_position( 1 ) * point_position( 1 ) ) * B_multiplicator;
        denominator = denominator + point_position( 2 ) * point_position( 2 ) * A_multiplicator;
        denominator = sqrt( denominator );
        double c_multiplicator = numerator / denominator;

        Vector3d position_updated = point_position * c_multiplicator;
        Vector3d add_to_get_segment( Cell_parametres::A_AXIS , Cell_parametres::A_AXIS , Cell_parametres::B_AXIS_Lower );

        Vector3d position_to_get_index = position_updated + add_to_get_segment;

        unsigned int x_index = position_to_get_index( 0 ) / this->density_surface_dynein.get_X_width();
        unsigned int y_index = position_to_get_index( 1 ) / this->density_surface_dynein.get_Y_width();
        unsigned int z_index = position_to_get_index( 2 ) / this->density_surface_dynein.get_Z_width();

        unsigned int neccessary_dimension = this->density_surface_dynein.get_neccessary_dimension();
        unsigned int ID_map_index = x_index * neccessary_dimension * neccessary_dimension + y_index * neccessary_dimension + z_index;

    try
    {
        std::vector<Vector3d> dynein_on_surface_points = this->density_surface_dynein.get_map().at( ID_map_index );
        for( unsigned int point_id = 0 ; point_id < dynein_on_surface_points.size() ; point_id ++ )
        {
            Vector3d dynein_position = dynein_on_surface_points.at( point_id );
	    double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );

            if( probability < Dynein::attach_rate_per_step )
            {
                this->density_surface_dynein.erase_dynein_point( ID_map_index , point_id );
                returned_points.push_back( dynein_position );
            }
        }
    }
    catch (const std::out_of_range& oor)
    {

    }
    return returned_points;
}



void Cell::stepping_and_detachment_of_all_microtubule_projection_real_dynein()
{
    this->MTOC.resize_from_originals(  );




//MICROTUBULES

    for( unsigned int micro_index = 0 ; micro_index < this->number_of_microtubules ; micro_index ++ )
    {
        if( this->array_Of_Microtubules[ micro_index ].get_dynein_index() == 20  )
        {
	    std::vector<Vector3d> motors_for_surface = this->array_Of_Microtubules[ micro_index ].control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_2();
	    IS_Capt_Shrinkage tmp = this->IS_Capture_Shrinkage;
            std::vector<Vector3d> vectors_to_be_added = this->Capture_Shrinkage_dynein.create_capture_shrinkage_points( motors_for_surface.size() , tmp );
            this->Capture_Shrinkage_dynein.add_dynein_points_to_compartment_with_number( vectors_to_be_added , 1 );
        }
        else if( this->array_Of_Microtubules[ micro_index ].get_dynein_index() == 9  )
        {
            this->array_Of_Microtubules[ micro_index ].resizeMicrotubule_with_different_tangent_lenghts();
            std::vector<Vector3d> motors_for_surface = this->array_Of_Microtubules[ micro_index ].stepping_detach_real_dynein_abscissa_projection();
            std::vector<Vector3d> points_margins;
    	    std::vector<Vector3d> vector_to_project = this->IS_Cortic_Sl.control_of_inside_and_margins_IS_1( motors_for_surface , points_margins );
 	    this->project_and_add_points_to_surface( vector_to_project );  
	    if( points_margins.size(  ) > 0 )
	    {
		std::vector<Vector3d> motors_for_surface_tmp =  this->density_surface_dynein.create_cortical_sliding_points(  points_margins.size(  )  , this->IS_Cortic_Sl );
                this->project_and_add_points_to_surface( motors_for_surface_tmp );
	    }
	    
	    if( motors_for_surface.size() != points_margins.size(  )  + vector_to_project.size() )
	    {
	    	cout<<"motors_for_surface.size() != points_margins.size(  )  + vector_to_project.size()"<<endl;
	    	cout<<"void Cell::stepping_and_detachment_of_all_microtubule_projection_real_dynein()"<<endl;
	    	cout<<"ERROR_ID = 16115151531"<<endl;
	        throw("");
	    }
 	      	   
        }
        else
        {
            this->array_Of_Microtubules[ micro_index ].resizeMicrotubule_with_different_tangent_lenghts();
        }
    }





}




unsigned int Cell::get_number_of_dynein_motors_capture_shrinkage_1()
{
    unsigned int sum_of_motor_in_micro = 0;
    for( unsigned int motor_id = 0 ; motor_id < this->number_of_microtubules; motor_id ++ )
    {
        if( this->array_Of_Microtubules[ motor_id ].get_dynein_index() == 20 )
        {
            unsigned int number_of_dynein_in_micro = this->array_Of_Microtubules[ motor_id ].get_Dynein_abscissa().size();
            sum_of_motor_in_micro = sum_of_motor_in_micro + number_of_dynein_in_micro;
        }
    }

    return sum_of_motor_in_micro;
}



void Cell::test_dynein( string param )
{

	unsigned int number_dynein_20 = this->get_number_of_dynein_motors_capture_shrinkage_1(  );
	unsigned int number_dynein_20_IS = this->Capture_Shrinkage_dynein.get_all_dynein_point(  ).size();
	unsigned int basic_number_capt_1 = this->Capture_Shrinkage_dynein.get_original_number_capt();





	if( number_dynein_20 + number_dynein_20_IS != basic_number_capt_1 )
	{
		cout<<param<<endl;
		cout<<"number_dynein_20 = "<<number_dynein_20<<endl;
		cout<<"number_dynein_20_IS = "<<number_dynein_20_IS<<endl;
		cout<<"basic_number_capt_1 = "<<basic_number_capt_1<<endl;
		cout<<"number_dynein_20 + number_dynein_20_IS ! = basic_number_capt_1 "<<endl;
		cout<<"void Cell_two_IS::test_dynein( string param )"<<endl;
		cout<<"cell_id = "<<this->cell_id<<endl;
		cout<<"ERROR_ID = 364156456132516"<<endl;
		throw("");
	}




	//cortical
	unsigned int vector_IS_number = this->density_surface_dynein.get_dynein_motors_number();
     	unsigned int original_cort_1 = this->density_surface_dynein.get_original_number_cort_1();

    	unsigned int IS_cort_1_counter = 0;
    	for( unsigned int micro_index = 0 ; micro_index < this->number_of_microtubules ; micro_index ++ )
    	{

        	if( this->array_Of_Microtubules[ micro_index ].get_dynein_index() == 9  )
        	{

	    		std::vector< Vector3d > motors_from_micro = this->array_Of_Microtubules[ micro_index ].get_Dynein_points_without_erasing();
			IS_cort_1_counter = IS_cort_1_counter + motors_from_micro.size();
       		}
     	}


     	if( vector_IS_number + IS_cort_1_counter !=  original_cort_1 )
     	{
		cout<<param<<endl;
		cout<<"IS_cort_1_counter = "<<IS_cort_1_counter<<endl;
		cout<<"vector_IS_number = "<<vector_IS_number<<endl;
		cout<<"original_cort_1 = "<<original_cort_1<<endl;
		cout<<"cell_id = "<<this->cell_id<<endl;
		cout<<"void Cell::test_dynein( string param )"<<endl;
		cout<<"ERROR_ID = 546413256"<<endl;
		throw("");
     	}


}




void Cell::detach_microtubule_NINE_and_project_points_on_surface( unsigned int microtubule )
{
    if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() != 9 )
    {
        cout<<"this->array_Of_Microtubules[ microtubule ].get_dynein_index() != 9"<<endl;
        cout<<"void Cell::detach_microtubule_NINE_and_project_points_on_surface( unsigned int microtubule )"<<endl;
        cout<<"microtubule = "<<microtubule<<endl;
        cout<<"ERROR_ID = 684635843548164"<<endl;
        throw("");
    }

    std::vector<Vector3d> motors_for_surface = this->array_Of_Microtubules[ microtubule ].get_Dynein_points_and_erase();

    if( motors_for_surface.size() > 0 )
    {
        this->project_and_add_points_to_surface( motors_for_surface );
    }

}




void Cell::check_caught_micro_IS_with_real_dynein( unsigned int microtubule )
{
    //this is the function for control, whether the microtubule 0 or 9 is caught
    // the microtubule has to determine the point, where it is depolimerized


    std::uniform_real_distribution<> distribution{ 0 , 1 };

    if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
    {
	cout<<"void Cell::check_caught_micro_IS_with_real_dynein( unsigned int microtubule )"<<endl;
	cout<<" this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 "<<endl;
	unsigned int ERROR_ID = 978161641;
	cout<<"ERROR_ID = "<<ERROR_ID<<endl;
	throw("");
    }

    unsigned int dynein_index_to_remember;

    for( unsigned int segment_id = 1 ; segment_id < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 ; segment_id ++ )
    {
	if( ( ( segment_id == 1 ) || ( segment_id == 2 ) ) && ( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() > 5 )   )
	{
		continue;
	}
        Vector3d bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id );
        Vector3d tangent;
        try
        {
            tangent = this->array_Of_Microtubules[ microtubule ].getTangent2( segment_id );
        }
        catch( int e )
        {
            cout<<"exception"<<endl;
            cout<<"void Cell::check_caught_micro_IS_with_real_dynein( unsigned int microtubule )"<<endl;
            throw("");
        }
        bool answer = this->IS_Capture_Shrinkage.check_IS_capture_segment_control( bead_position , tangent );

        if( answer == true  )
        {
            std::vector<Vector3d> dynein_points = this->Capture_Shrinkage_dynein.get_dynein_points( 1 );
            std::vector< Vector3d > replace_motors;

            for( unsigned int dynein_index = 0 ; dynein_index < dynein_points.size() ; dynein_index ++ )
            {

                Vector3d motor_position = dynein_points[ dynein_index ];
                Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
                double distance = distance_point_segment( bead_position , tangent , motor_position , cl_p_of_seg );

		unsigned int number_of_generator = omp_get_thread_num();
    		double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );

                if( probability < Dynein_real::calculate_attachment_probability_per_time_step( distance ) )
                {
			dynein_index_to_remember = dynein_index;
			//erasing all dynein from nine
                        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9  )
                        {
                            std::vector< Vector3d  > points_positions = this->array_Of_Microtubules[ microtubule ].get_Dynein_points_and_erase();
                            this->project_and_add_points_to_surface( points_positions );

                        }



			double radius_tmp =this->Capture_Shrinkage_dynein.triangular_distribution( IS_Capture_shrinkage_param::radius );
			std::uniform_real_distribution<> distribution_angle( 0.0 , 2.0 * sim_of_Cell::PI );
         	        unsigned int number_of_generator = omp_get_thread_num();
			double azimutal_angle_tmp = distribution_angle(  mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
			double x = radius_tmp * cos( azimutal_angle_tmp );
			double y = radius_tmp * sin( azimutal_angle_tmp );
                        double z =  this->a_axis;
                        Vector3d position( x , y , z );
			//rotation to IS
			//I compute angle between z_axis and IS axis and rotate the point
			Vector3d z_axis( 0.0 , 0.0 , 1.0 );
			double angle = acos( z_axis.dot( this->IS_Capture_Shrinkage.get_axis_of_IS() *  ( -1.0 ) ) / ( this->IS_Capture_Shrinkage.get_axis_of_IS().norm() ) );
			Vector3d axis_of_rotation = this->IS_Capture_Shrinkage.get_axis_of_IS().cross( z_axis );
			axis_of_rotation = axis_of_rotation / axis_of_rotation.norm();

			Quaternion<double> q;
			q = AngleAxis<double>( angle , axis_of_rotation );
			position = q * position;


                        Vector3d caught_position = this->project_point_on_surface( position );
                        this->array_Of_Microtubules[ microtubule ].set_IS_position_catching( caught_position  );
                        unsigned int number_of_beads = this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1;

                        for( unsigned int sub_id = number_of_beads ; sub_id > segment_id ; sub_id -- )
                        {
                        	this->array_Of_Microtubules[ microtubule ].Substract_one_Bead();
                        }
                        this->array_Of_Microtubules[ microtubule ].add_one_final_Bead_catching_position();

			this->array_Of_Microtubules[ microtubule ].set_dynein_index( 20 );
                        this->array_Of_Microtubules[ microtubule ].set_lenght_after_catching( this->array_Of_Microtubules[ microtubule ].get_lenght_of_microtubule() );
			this->array_Of_Microtubules[ microtubule ].set_lenght_of_tangents();
			break;
		}
	    }

	}
    }


    if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20  )
    {

    	std::vector<Vector3d> dynein_points = this->Capture_Shrinkage_dynein.get_dynein_points( 1 );
        Vector3d motor_position = dynein_points[ dynein_index_to_remember ];


    	unsigned int bead_index = this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2;
    	Vector3d bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( bead_index );
    	Vector3d tangent;
    	try
    	{
        	tangent = this->array_Of_Microtubules[ microtubule ].getTangent2( bead_index );
    	}
    	catch( int e )
    	{
        	cout<<"exception"<<endl;
        	cout<<"void Cell::check_caught_micro_IS_with_real_dynein( unsigned int microtubule )"<<endl;
        	throw("");
    	}


   	Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
   	distance_point_segment( bead_position , tangent , motor_position , cl_p_of_seg );
   	double distance_lower_bead = 0;
   	for( unsigned int i = 0 ; i < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 ; i ++ )
   	{
   		distance_lower_bead = distance_lower_bead + this->array_Of_Microtubules[ microtubule ].getTangent2( i ).norm();
   	}

   	double distance_point_attachment = ( bead_position - cl_p_of_seg ).norm();
   	double abscissa = ( distance_lower_bead + distance_point_attachment );
   	if( distance_point_attachment == tangent.norm() )
   	{
   		abscissa = abscissa - 1e-8;
   	}

   	if( abscissa > this->array_Of_Microtubules[ microtubule ].get_lenght_of_microtubule() )
   	{
   		throw("");
  	 }

  	std::pair < Vector3d , double > tmp_pair( cl_p_of_seg , abscissa );
   	this->array_Of_Microtubules[ microtubule ].add_pair( tmp_pair );

    	dynein_points.erase (dynein_points.begin() + dynein_index_to_remember );
    	this->Capture_Shrinkage_dynein.set_dynein_points( 1 , dynein_points );

    }


}





void Cell::check_caught_micro_IS_with_real_dynein_all_micro( )
{
    Vector3d center_of_MTOC = this->MTOC.get_center();
    Vector3d center_of_IS = this->IS_Capture_Shrinkage.get_center_of_IS_front();
    double distance = ( center_of_MTOC - center_of_IS ).norm();

    if( distance > 0.0e-6 )
    {
        for( unsigned int micro_id = 0 ; micro_id < this->get_microtubule_number() ; micro_id ++ )
        {
            if( this->array_Of_Microtubules[ micro_id ].get_dynein_index() != 20 )
            {
                this->check_caught_micro_IS_with_real_dynein( micro_id );
            }
        }
    }
}


void Cell::microtubule_catch_pair_abscissa_real_dynein_in_IS( unsigned int microtubule )
{

    std::uniform_real_distribution<> distribution{ 0 , 1 };
    unsigned int number_of_generator = omp_get_thread_num();


    for( unsigned int segment_id = 1 ; segment_id < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 ; segment_id ++ )
    {


        Vector3d bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id );
        Vector3d tangent;
        try
        {
            tangent = this->array_Of_Microtubules[ microtubule ].getTangent2( segment_id );
        }
        catch( int e )
        {
            cout<<"exception"<<endl;
            cout<<"void void Microtubule::set_lenght_of_tangents()"<<endl;
            throw("");
        }


        bool answer = this->IS_Capture_Shrinkage.check_IS_capture_segment_control( bead_position , tangent );
        if( ( answer == true ) && ( this->Capture_Shrinkage_dynein.get_dynein_points( 1 ).size() > 0 ) )
        {

            std::vector<Vector3d> replace_motors;
            std::vector<Vector3d> motors_to_attach;
            std::vector<Vector3d> dynein_points = this->Capture_Shrinkage_dynein.get_dynein_points( 1 );

            for( unsigned int dynein_index = 0 ; dynein_index < dynein_points.size() ; dynein_index ++ )
            {

                Vector3d motor_position = dynein_points[ dynein_index ];
                Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
                double distance = distance_point_segment( bead_position , tangent , motor_position , cl_p_of_seg );

		double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );




                if( probability < Dynein_real::calculate_attachment_probability_per_time_step( distance ) )
                {

                    motors_to_attach.push_back( motor_position );
                }
                else
                {
                    replace_motors.push_back( motor_position );
                }

            }


            for( unsigned int counter_id = 0 ; counter_id < motors_to_attach.size() ; counter_id ++ )
            {
                Vector3d motor_position = motors_to_attach[ counter_id ];
                Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
                distance_point_segment( bead_position , tangent , motor_position , cl_p_of_seg );

                double distance_lower_bead = this->array_Of_Microtubules[ microtubule ].get_distance_to_lower_bead_with_index( segment_id );
                double distance_point_attachment = ( bead_position - cl_p_of_seg ).norm();
                double abscissa = ( distance_lower_bead + distance_point_attachment );

                if( abs( abscissa - this->array_Of_Microtubules[ microtubule ].get_lenght_of_microtubule() ) < Dynein_real::step / 2.0 )
                {
                    abscissa = abscissa - Dynein_real::step;
                }

                std::pair < Vector3d , double > tmp_pair( cl_p_of_seg , abscissa );
                this->array_Of_Microtubules[ microtubule ].add_pair( tmp_pair );
            }


            if( replace_motors.size() == 0 )
            {
                this->Capture_Shrinkage_dynein.erase_vector_points_with_key( 1 );
            }
	    else
	    {
                this->Capture_Shrinkage_dynein.set_dynein_points( 1 , replace_motors );
            }
        }

    }


}


void Cell::microtubule_catch_pair_abscissa_real_dynein_in_IS_all_micro( )
{
    for( unsigned int microtubule = 0 ; microtubule < this->get_microtubule_number() ; microtubule ++ )
    {
            if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
            {
                this->microtubule_catch_pair_abscissa_real_dynein_in_IS( microtubule );
            }
    }

}

void Cell::control_length_of_micro_IS()
{


    for( unsigned int microtubule = 0 ; microtubule < this->get_microtubule_number() ; microtubule ++ )
    {
        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
        {

            double original_lenght = this->array_Of_Microtubules[ microtubule ].get_lenght_after_catching();
            double current_lenght = this->array_Of_Microtubules[ microtubule ].get_lenght_of_microtubule();

            if( current_lenght < original_lenght )
            {
                this->array_Of_Microtubules[ microtubule ].set_lenght_after_catching( current_lenght );
            }

            if( current_lenght > original_lenght * IS_Capture_shrinkage_param::procentage_constant )
            {

                std::vector<Vector3d> erased_vectors = this->array_Of_Microtubules[ microtubule ].get_dynein_points_in_IS_and_erase();
		IS_Capt_Shrinkage tmp = this->IS_Capture_Shrinkage;
                std::vector<Vector3d> vectors_to_be_added = this->Capture_Shrinkage_dynein.create_capture_shrinkage_points( erased_vectors.size() , tmp );
                this->Capture_Shrinkage_dynein.add_dynein_points_to_compartment_with_number( vectors_to_be_added , 1 );
                this->array_Of_Microtubules[ microtubule ].set_dynein_index( 0 );

		Vector3d posledni_tangent = this->array_Of_Microtubules[ microtubule ].get_last_Tangent();
		Vector3d position = this->array_Of_Microtubules[ microtubule ].getPoint( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 );
		position = position + posledni_tangent / posledni_tangent.norm() * sim_of_Cell::resting_distance;
		this->array_Of_Microtubules[ microtubule ].setPoint( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 , position );
		this->array_Of_Microtubules[ microtubule ].set_lenght_of_tangents();


            }

            if( this->array_Of_Microtubules[ microtubule ].get_last_Tangent().norm() > 4 * sim_of_Cell::resting_distance )
            {
                std::vector<Vector3d> erased_vectors = this->array_Of_Microtubules[ microtubule ].get_dynein_points_in_IS_and_erase();
		IS_Capt_Shrinkage tmp = this->IS_Capture_Shrinkage;
                std::vector<Vector3d> vectors_to_be_added = this->Capture_Shrinkage_dynein.create_capture_shrinkage_points( erased_vectors.size() , tmp );
                this->Capture_Shrinkage_dynein.add_dynein_points_to_compartment_with_number( vectors_to_be_added , 1 );
                this->array_Of_Microtubules[ microtubule ].set_dynein_index( 0 );

		Vector3d posledni_tangent = this->array_Of_Microtubules[ microtubule ].get_last_Tangent();
		Vector3d position = this->array_Of_Microtubules[ microtubule ].getPoint( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 2 );
		position = position + posledni_tangent / posledni_tangent.norm() * sim_of_Cell::resting_distance;
		this->array_Of_Microtubules[ microtubule ].setPoint( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 , position );
		this->array_Of_Microtubules[ microtubule ].set_lenght_of_tangents();

            }





        }
    }

}





void Cell::control_length_of_micro_IS_2()
{
    for( unsigned int microtubule = 0 ; microtubule < this->get_microtubule_number() ; microtubule ++ )
    {
        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
        {
            std::vector<Vector3d> erased_vectors = this->array_Of_Microtubules[ microtubule ].control_length_of_micro_IS_2();
	    if( erased_vectors.size() >= 1 )
	    {
	    	IS_Capt_Shrinkage tmp = this->IS_Capture_Shrinkage;
            	std::vector<Vector3d> vectors_to_be_added = this->Capture_Shrinkage_dynein.create_capture_shrinkage_points( erased_vectors.size() , tmp );
            	this->Capture_Shrinkage_dynein.add_dynein_points_to_compartment_with_number( vectors_to_be_added , 1 );
            }

        }
    }

}




void Cell::control_microtubule_20_three_beads_dyneins()
{
    for( unsigned int microtubule = 0 ; microtubule < this->get_microtubule_number() ; microtubule ++ )
    {
        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 20 )
        {
		if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() == 3 )
		{
			if( this->array_Of_Microtubules[ microtubule ].get_number_of_dynein_points_IS() == 0 )
			{
				this->array_Of_Microtubules[ microtubule ].set_dynein_index( 0 );
			}
		}
        }
    }
}







void Cell::microtubule_random_dynein_on_surface_catch_pair_abscissa_2( unsigned int microtubule )
{

    std::uniform_real_distribution<> distribution{ 0 , 1 };
    unsigned int number_of_generator = omp_get_thread_num();
    for( unsigned int segment_id = 0 ; segment_id < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints()  - 1 ; segment_id ++ )
    {

        //in this for cycle the points are projected to the surface ( lower_bead_position )
        Vector3d lower_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id );
        Vector3d lower_bead_position_TMP = lower_bead_position;

        lower_bead_position = this->project_point_on_surface( lower_bead_position );

        Vector3d upper_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id + 1 );
        Vector3d upper_bead_position_TMP = upper_bead_position;

        upper_bead_position = this->project_point_on_surface( upper_bead_position );

        Vector3d tangent = upper_bead_position - lower_bead_position;
        Vector3d tangent_TMP = this->array_Of_Microtubules[ microtubule ].getTangent2( segment_id );

        std::vector<double> final_vector;
        std::vector<double> x_tran_vec;
        std::vector<double> y_tran_vec;
        std::vector<double> z_tran_vec;


        int lower_index_X = ( int ) ( lower_bead_position[ 0 ] / this->density_surface_dynein.get_X_width() );
        int upper_index_X = ( int ) ( upper_bead_position[ 0 ] / this->density_surface_dynein.get_X_width() );
        unsigned int x_transitions = abs( upper_index_X - lower_index_X );

        int smer_x;
        if( tangent[ 0 ] > 0 )
        {
            smer_x = 1;
        }
        else
        {
            smer_x = - 1;
        }
        if( x_transitions > 0 )
        {

            for( int x_index = 0 ; x_index < x_transitions ; x_index ++ )
            {
                int index_boundary;
                if( smer_x * lower_bead_position[ 0 ] > 0 )
                {
                    index_boundary = lower_index_X + ( x_index + 1 ) * smer_x;
                }
                else
                {
                    index_boundary = lower_index_X + ( x_index ) * smer_x;
                }
                double boundary = ( double ) index_boundary * this->density_surface_dynein.get_X_width();
                double distance = boundary - lower_bead_position[ 0 ];
                double t_tmp = distance / tangent[ 0 ];
                x_tran_vec.push_back( t_tmp );

            }
        }

        int lower_index_Y = ( int ) ( lower_bead_position[ 1 ] / this->density_surface_dynein.get_Y_width() );
        int upper_index_Y = ( int ) ( upper_bead_position[ 1 ] / this->density_surface_dynein.get_Y_width() );
        unsigned int y_transitions = abs( upper_index_Y - lower_index_Y );

        int smer_y;
        if( tangent[ 1 ] > 0 )
        {
            smer_y = 1;
        }
        else
        {
            smer_y = - 1;
        }
        if( y_transitions > 0 )
        {

            for( int y_index = 0 ; y_index < y_transitions ; y_index ++ )
            {
                int index_boundary;
                if( smer_y * lower_bead_position[ 1 ] > 0 )
                {
                    index_boundary = lower_index_Y + ( y_index + 1 ) * smer_y;
                }
                else
                {
                    index_boundary = lower_index_Y + ( y_index ) * smer_y;
                }
                double boundary = ( double ) index_boundary * this->density_surface_dynein.get_Y_width();
                double distance = boundary - lower_bead_position[ 1 ];
                double t_tmp = distance / tangent[ 1 ];
                y_tran_vec.push_back( t_tmp );

            }
        }

        int lower_index_Z = ( int ) ( lower_bead_position[ 2 ] / this->density_surface_dynein.get_Z_width() );
        int upper_index_Z = ( int ) ( upper_bead_position[ 2 ] / this->density_surface_dynein.get_Z_width() );
        unsigned int z_transitions = abs( upper_index_Z - lower_index_Z );

        int smer_z;
        if( tangent[ 2 ] > 0 )
        {
            smer_z = 1;
        }
        else
        {
            smer_z = - 1;
        }
        if( z_transitions > 0 )
        {

            for( int z_index = 0 ; z_index < z_transitions ; z_index ++ )
            {
                int index_boundary;
                if( smer_z * lower_bead_position[ 2 ] > 0 )
                {
                    index_boundary = lower_index_Z + ( z_index + 1 ) * smer_z;
                }
                else
                {
                    index_boundary = lower_index_Z + ( z_index ) * smer_z;
                }
                double boundary = ( double ) index_boundary * this->density_surface_dynein.get_Z_width();
                double distance = boundary - lower_bead_position[ 2 ];
                double t_tmp = distance / tangent[ 2 ];
                z_tran_vec.push_back( t_tmp );

            }
        }

        final_vector.push_back( 0 );
        final_vector.insert( final_vector.end(), x_tran_vec.begin(), x_tran_vec.end() );
        final_vector.insert( final_vector.end(), y_tran_vec.begin(), y_tran_vec.end() );
        final_vector.insert( final_vector.end(), z_tran_vec.begin(), z_tran_vec.end() );

        std::vector<double> compartment_IDs;

        std::sort ( final_vector.begin(), final_vector.end() );

        double push_par = 1e-10;
        for( unsigned int parameter_index = 0 ; parameter_index < final_vector.size() ; parameter_index ++ )
        {
            double parameter = final_vector.at( parameter_index ) + push_par;

            Vector3d Vector_point = lower_bead_position + tangent * parameter;
            unsigned int compartment_id = get_dynein_compartment_id( Vector_point );
            compartment_IDs.push_back( compartment_id );
        }


        for( unsigned int comp_index = 0 ; comp_index < compartment_IDs.size() ; comp_index ++ )
        {
            unsigned int compartment_id = compartment_IDs.at( comp_index );
            std::vector<Vector3d> dynein_points = this->get_dynein_in_compartment( compartment_id );

            for( unsigned int bod_id  = 0 ; bod_id < dynein_points.size() ; bod_id ++ )
            {

                Vector3d motor_position = dynein_points.at( bod_id );
                Vector3d closest_point_of_segment( 0.0 , 0.0 , 0.0 );


                double distance = distance_point_segment( lower_bead_position , tangent , motor_position , closest_point_of_segment );
                distance_point_segment( lower_bead_position_TMP , tangent_TMP , motor_position , closest_point_of_segment );
                //this vector will contain motors that will be kept in map
                std::vector< Vector3d > replace_motors;

                if( abs( distance ) < Dynein::L_0 )
                {
		    double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );



                    if( probability < Dynein::attach_rate_per_step )
                    {

                        Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
                        distance_point_segment( lower_bead_position_TMP , tangent_TMP , closest_point_of_segment , cl_p_of_seg );
                        double distance_from_start_to_lower_bead = this->array_Of_Microtubules[ microtubule ].getRestDist() * ( double ) segment_id;
                        double distance_point_attachment = ( lower_bead_position_TMP - cl_p_of_seg ).norm();
                        double abscissa = ( distance_from_start_to_lower_bead + distance_point_attachment );

                        std::pair < Vector3d , double > tmp_pair( cl_p_of_seg , abscissa );
                        this->array_Of_Microtubules[ microtubule ].add_pair( tmp_pair );
                        this->array_Of_Microtubules[ microtubule ].set_dynein_index( 9 );

                    }
                    else
                    {
                        replace_motors.push_back( motor_position );
                    }
                }
                replace_motors.push_back( motor_position );
                this->density_surface_dynein.set_dynein_points( compartment_id , replace_motors );


            }
        }

    }
}


void Cell::microtubule_random_dynein_on_surface_catch_pair_abscissa_plane( unsigned int microtubule )
{

    std::uniform_real_distribution<> distribution{ 0 , 1 };
    unsigned int number_of_generator = omp_get_thread_num();
    for( unsigned int segment_id = 0 ; segment_id < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints()  - 1 ; segment_id ++ )
    {
        Vector3d lower_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id );
        Vector3d lower_bead_position_TMP = lower_bead_position;
        lower_bead_position = this->project_point_on_surface( lower_bead_position );

        Vector3d upper_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id + 1 );
        Vector3d upper_bead_position_TMP = upper_bead_position;

        upper_bead_position = this->project_point_on_surface( upper_bead_position );
        Vector3d tangent = upper_bead_position - lower_bead_position;
        Vector3d tangent_TMP = this->array_Of_Microtubules[ microtubule ].getTangent2( segment_id );

        std::vector<double> final_vector;
        std::vector<double> x_tran_vec;
        std::vector<double> y_tran_vec;
        std::vector<double> z_tran_vec;

        int lower_index_X = ( int ) ( lower_bead_position[ 0 ] / this->density_surface_dynein.get_X_width() );
        int upper_index_X = ( int ) ( upper_bead_position[ 0 ] / this->density_surface_dynein.get_X_width() );
        unsigned int x_transitions = abs( upper_index_X - lower_index_X );

        int smer_x;
        if( tangent[ 0 ] > 0 )
        {
            smer_x = 1;
        }
        else
        {
            smer_x = - 1;
        }
        if( x_transitions > 0 )
        {


            //je tam int kvuli pretypovani
            for( int x_index = 0 ; x_index < x_transitions ; x_index ++ )
            {
                int index_boundary;
                //tahle podminka mi rika, jestli prvni hranice je ta zakladni, nebo zakladni + 1
                if( smer_x * lower_bead_position[ 0 ] > 0 )
                {
                    index_boundary = lower_index_X + ( x_index + 1 ) * smer_x;
                }
                else
                {
                    index_boundary = lower_index_X + ( x_index ) * smer_x;
                }
                double boundary = ( double ) index_boundary * this->density_surface_dynein.get_X_width();
                double distance = boundary - lower_bead_position[ 0 ];
                double t_tmp = distance / tangent[ 0 ];
                x_tran_vec.push_back( t_tmp );

            }
        }

        int lower_index_Y = ( int ) ( lower_bead_position[ 1 ] / this->density_surface_dynein.get_Y_width() );
        int upper_index_Y = ( int ) ( upper_bead_position[ 1 ] / this->density_surface_dynein.get_Y_width() );
        unsigned int y_transitions = abs( upper_index_Y - lower_index_Y );

        int smer_y;
        if( tangent[ 1 ] > 0 )
        {
            smer_y = 1;
        }
        else
        {
            smer_y = - 1;
        }
        if( y_transitions > 0 )
        {

            for( int y_index = 0 ; y_index < y_transitions ; y_index ++ )
            {
                int index_boundary;
                if( smer_y * lower_bead_position[ 1 ] > 0 )
                {
                    index_boundary = lower_index_Y + ( y_index + 1 ) * smer_y;
                }
                else
                {
                    index_boundary = lower_index_Y + ( y_index ) * smer_y;
                }
                double boundary = ( double ) index_boundary * this->density_surface_dynein.get_Y_width();
                double distance = boundary - lower_bead_position[ 1 ];
                double t_tmp = distance / tangent[ 1 ];
                y_tran_vec.push_back( t_tmp );

            }
        }


        int lower_index_Z = ( int ) ( lower_bead_position[ 2 ] / this->density_surface_dynein.get_Z_width() );
        int upper_index_Z = ( int ) ( upper_bead_position[ 2 ] / this->density_surface_dynein.get_Z_width() );
        unsigned int z_transitions = abs( upper_index_Z - lower_index_Z );

        int smer_z;
        if( tangent[ 2 ] > 0 )
        {
            smer_z = 1;
        }
        else
        {
            smer_z = - 1;
        }
        if( z_transitions > 0 )
        {

            for( int z_index = 0 ; z_index < z_transitions ; z_index ++ )
            {
                int index_boundary;
                if( smer_z * lower_bead_position[ 2 ] > 0 )
                {
                    index_boundary = lower_index_Z + ( z_index + 1 ) * smer_z;
                }
                else
                {
                    index_boundary = lower_index_Z + ( z_index ) * smer_z;
                }
                double boundary = ( double ) index_boundary * this->density_surface_dynein.get_Z_width();
                double distance = boundary - lower_bead_position[ 2 ];
                double t_tmp = distance / tangent[ 2 ];
                z_tran_vec.push_back( t_tmp );

            }
        }

        final_vector.push_back( 0 );
        final_vector.insert( final_vector.end(), x_tran_vec.begin(), x_tran_vec.end() );
        final_vector.insert( final_vector.end(), y_tran_vec.begin(), y_tran_vec.end() );
        final_vector.insert( final_vector.end(), z_tran_vec.begin(), z_tran_vec.end() );

        std::vector<double> compartment_IDs;
        std::sort ( final_vector.begin(), final_vector.end() );

        double push_par = 1e-10;
        for( unsigned int parameter_index = 0 ; parameter_index < final_vector.size() ; parameter_index ++ )
        {
            double parameter = final_vector.at( parameter_index ) + push_par;

            Vector3d Vector_point = lower_bead_position + tangent * parameter;
            unsigned int compartment_id = this->density_surface_dynein.get_dynein_compartment_id_projected( Vector_point );
            if( compartment_id != 0 )
            {
                compartment_IDs.push_back( compartment_id );
            }

        }


        Vector3d axis_of_plane = lower_bead_position_TMP.cross( upper_bead_position_TMP );
        axis_of_plane = axis_of_plane / axis_of_plane.norm();
        for( unsigned int comp_index = 0 ; comp_index < compartment_IDs.size() ; comp_index ++ )
        {
            unsigned int compartment_id = compartment_IDs.at( comp_index );
            std::vector<Vector3d> dynein_points = this->get_dynein_in_compartment( compartment_id );

            std::vector< Vector3d > replace_motors;
            for( unsigned int bod_id  = 0 ; bod_id < dynein_points.size() ; bod_id ++ )
            {

                Vector3d motor_position = dynein_points.at( bod_id );
                double distance = distance_plane_point( axis_of_plane , lower_bead_position_TMP , motor_position );
                if( motor_position( 2 ) < -5.0e-6 ) //!!!!!!NARAZNIK
                {
                    replace_motors.push_back( motor_position );
                    continue;
                }

                if( abs( distance ) < Dynein::L_0 )
                {

		    double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
                    if( probability < Dynein::attach_rate_per_step )
                    {
                        Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
                        double distance_tmp = distance_point_segment( lower_bead_position_TMP , tangent_TMP , motor_position , cl_p_of_seg );

                        if( distance_tmp > IS_Dynein_Cell_surface::treshold_cut_tangent )
                        {
                            replace_motors.push_back( motor_position );
                            continue;
                        }

                        double distance_from_start_to_lower_bead = this->array_Of_Microtubules[ microtubule ].getRestDist() * ( double ) segment_id;

                        Vector3d distance_cl_p_of_seg_motor_position = ( motor_position - cl_p_of_seg );
                        Vector3d perpendicular_tangent = distance_cl_p_of_seg_motor_position / distance_cl_p_of_seg_motor_position.norm() * Dynein::L_0;
                        Vector3d motor_2 = cl_p_of_seg + perpendicular_tangent / 2.5;

                        double distance_point_attachment = ( lower_bead_position_TMP - cl_p_of_seg ).norm();
                        double abscissa = ( distance_from_start_to_lower_bead + distance_point_attachment );
                        double distance = (double) ( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 ) * this->array_Of_Microtubules[ microtubule ].getRestDist();


                        if( distance_point_attachment > this->array_Of_Microtubules[ microtubule ].getRestDist() )
                        {
                            replace_motors.push_back( motor_position );
                        }
                        else
                        {
                            Vector3d motor_position_2 = motor_position / motor_position.norm() * cl_p_of_seg.norm();
                            std::pair < Vector3d , double > tmp_pair( motor_position_2 , abscissa );
                            this->array_Of_Microtubules[ microtubule ].add_pair( tmp_pair );
                            this->array_Of_Microtubules[ microtubule ].set_dynein_index( 9 );
                        }
                    }
                    else
                    {
                        replace_motors.push_back( motor_position );
                    }
                }
                else
                {
                    replace_motors.push_back( motor_position );
                }
            }

            this->density_surface_dynein.set_dynein_points( compartment_id , replace_motors );

        }

    }
}


void Cell::microtubule_random_dynein_on_surface_catch_pair_abscissa_space( unsigned int microtubule )
{
    std::uniform_real_distribution<> distribution{ 0 , 1 };
    unsigned int number_of_generator = omp_get_thread_num();


    for( unsigned int segment_id = 0 ; segment_id < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints()  - 1 ; segment_id ++ )
    {
        Vector3d lower_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id );
        lower_bead_position = this->project_point_on_surface( lower_bead_position );
        Vector3d upper_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id + 1 );
        upper_bead_position = this->project_point_on_surface( upper_bead_position );
        Vector3d tangent = upper_bead_position - lower_bead_position;

        std::vector<double> final_vector;
        std::vector<double> x_tran_vec;
        std::vector<double> y_tran_vec;
        std::vector<double> z_tran_vec;

        int lower_index_X = ( int ) ( lower_bead_position[ 0 ] / this->density_surface_dynein.get_X_width() );
        int upper_index_X = ( int ) ( upper_bead_position[ 0 ] / this->density_surface_dynein.get_X_width() );
        unsigned int x_transitions = abs( upper_index_X - lower_index_X );

        int smer_x;
        if( tangent[ 0 ] > 0 )
        {
            smer_x = 1;
        }
        else
        {
            smer_x = - 1;
        }
        if( x_transitions > 0 )
        {

            for( int x_index = 0 ; x_index < x_transitions ; x_index ++ )
            {
                int index_boundary;
                if( smer_x * lower_bead_position[ 0 ] > 0 )
                {
                    index_boundary = lower_index_X + ( x_index + 1 ) * smer_x;
                }
                else
                {
                    index_boundary = lower_index_X + ( x_index ) * smer_x;
                }
                double boundary = ( double ) index_boundary * this->density_surface_dynein.get_X_width();
                double distance = boundary - lower_bead_position[ 0 ];
                double t_tmp = distance / tangent[ 0 ];
                x_tran_vec.push_back( t_tmp );

            }
        }

        int lower_index_Y = ( int ) ( lower_bead_position[ 1 ] / this->density_surface_dynein.get_Y_width() );
        int upper_index_Y = ( int ) ( upper_bead_position[ 1 ] / this->density_surface_dynein.get_Y_width() );
        unsigned int y_transitions = abs( upper_index_Y - lower_index_Y );

        int smer_y;
        if( tangent[ 1 ] > 0 )
        {
            smer_y = 1;
        }
        else
        {
            smer_y = - 1;
        }
        if( y_transitions > 0 )
        {

            for( int y_index = 0 ; y_index < y_transitions ; y_index ++ )
            {
                int index_boundary;
                if( smer_y * lower_bead_position[ 1 ] > 0 )
                {
                    index_boundary = lower_index_Y + ( y_index + 1 ) * smer_y;
                }
                else
                {
                    index_boundary = lower_index_Y + ( y_index ) * smer_y;
                }
                double boundary = ( double ) index_boundary * this->density_surface_dynein.get_Y_width();
                double distance = boundary - lower_bead_position[ 1 ];
                double t_tmp = distance / tangent[ 1 ];
                y_tran_vec.push_back( t_tmp );

            }
        }

        int lower_index_Z = ( int ) ( lower_bead_position[ 2 ] / this->density_surface_dynein.get_Z_width() );
        int upper_index_Z = ( int ) ( upper_bead_position[ 2 ] / this->density_surface_dynein.get_Z_width() );
        unsigned int z_transitions = abs( upper_index_Z - lower_index_Z );

        int smer_z;
        if( tangent[ 2 ] > 0 )
        {
            smer_z = 1;
        }
        else
        {
            smer_z = - 1;
        }
        if( z_transitions > 0 )
        {

            for( int z_index = 0 ; z_index < z_transitions ; z_index ++ )
            {
                int index_boundary;
                if( smer_z * lower_bead_position[ 2 ] > 0 )
                {
                    index_boundary = lower_index_Z + ( z_index + 1 ) * smer_z;
                }
                else
                {
                    index_boundary = lower_index_Z + ( z_index ) * smer_z;
                }
                double boundary = ( double ) index_boundary * this->density_surface_dynein.get_Z_width();
                double distance = boundary - lower_bead_position[ 2 ];
                double t_tmp = distance / tangent[ 2 ];
                z_tran_vec.push_back( t_tmp );

            }
        }

        final_vector.push_back( 0 );
        final_vector.insert( final_vector.end(), x_tran_vec.begin(), x_tran_vec.end() );
        final_vector.insert( final_vector.end(), y_tran_vec.begin(), y_tran_vec.end() );
        final_vector.insert( final_vector.end(), z_tran_vec.begin(), z_tran_vec.end() );

        std::vector<double> compartment_IDs;
        std::sort ( final_vector.begin(), final_vector.end() );

        double push_par = 1e-10;
        for( unsigned int parameter_index = 0 ; parameter_index < final_vector.size() ; parameter_index ++ )
        {
            double parameter = final_vector.at( parameter_index ) + push_par;
            Vector3d Vector_point = lower_bead_position + tangent * parameter;
            unsigned int compartment_id = this->density_surface_dynein.get_dynein_compartment_id_projected( Vector_point );
            if( compartment_id != 0 )
            {
                compartment_IDs.push_back( compartment_id );
            }

        }



        for( unsigned int comp_index = 0 ; comp_index < compartment_IDs.size() ; comp_index ++ )
        {
            unsigned int compartment_id = compartment_IDs.at( comp_index );
            std::vector<Vector3d> dynein_points = this->get_dynein_in_compartment( compartment_id );

            std::vector< Vector3d > replace_motors;
            for( unsigned int bod_id  = 0 ; bod_id < dynein_points.size() ; bod_id ++ )
            {
                Vector3d motor_position = dynein_points.at( bod_id );
                if( motor_position( 2 ) < -5.0e-6 ) //!!!!!!NARAZNIK
                {
                    replace_motors.push_back( motor_position );
                    continue;
                }

                Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
                double distance = distance_point_segment( lower_bead_position , tangent , motor_position , cl_p_of_seg );
                if( abs( distance ) < Dynein::L_0 )
                {
		    double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );


                    if( probability < Dynein::attach_rate_per_step )
                    {

                        double distance_from_start_to_lower_bead = this->array_Of_Microtubules[ microtubule ].getRestDist() * ( double ) segment_id;
                        double distance_point_attachment = ( lower_bead_position - cl_p_of_seg ).norm();
                        double abscissa = ( distance_from_start_to_lower_bead + distance_point_attachment );

                        if( distance_point_attachment > this->array_Of_Microtubules[ microtubule ].getRestDist() )
                        {
                            replace_motors.push_back( motor_position );
                        }
                        else
                        {
                            std::pair < Vector3d , double > tmp_pair( motor_position , abscissa );
                            this->array_Of_Microtubules[ microtubule ].add_pair( tmp_pair );
                            this->array_Of_Microtubules[ microtubule ].set_dynein_index( 9 );
                        }
                    }
                    else
                    {
                        replace_motors.push_back( motor_position );
                    }
                }
                else
                {
                    replace_motors.push_back( motor_position );
                }
            }

            this->density_surface_dynein.set_dynein_points( compartment_id , replace_motors );
        }

    }

}


void Cell::microtubule_random_dynein_on_surface_catch_pair_abscissa_space_line( unsigned int microtubule )
{


    std::uniform_real_distribution<> distribution{ 0 , 1 };
    unsigned int number_of_generator = omp_get_thread_num();

    for( unsigned int segment_id = 0 ; segment_id < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints()  - 1 ; segment_id ++ )
    {
        Vector3d lower_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id );
        lower_bead_position = this->project_point_on_surface( lower_bead_position );
        Vector3d upper_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id + 1 );
        upper_bead_position = this->project_point_on_surface( upper_bead_position );
        Vector3d tangent = upper_bead_position - lower_bead_position;


        std::vector<double> final_vector;
        std::vector<double> x_tran_vec;
        std::vector<double> y_tran_vec;
        std::vector<double> z_tran_vec;


        int lower_index_X = ( int ) ( lower_bead_position[ 0 ] / this->density_surface_dynein.get_X_width() );
        int upper_index_X = ( int ) ( upper_bead_position[ 0 ] / this->density_surface_dynein.get_X_width() );
        unsigned int x_transitions = abs( upper_index_X - lower_index_X );

        int smer_x;
        if( tangent[ 0 ] > 0 )
        {
            smer_x = 1;
        }
        else
        {
            smer_x = - 1;
        }
        if( x_transitions > 0 )
        {


            for( int x_index = 0 ; x_index < x_transitions ; x_index ++ )
            {
                int index_boundary;
                if( smer_x * lower_bead_position[ 0 ] > 0 )
                {
                    index_boundary = lower_index_X + ( x_index + 1 ) * smer_x;
                }
                else
                {
                    index_boundary = lower_index_X + ( x_index ) * smer_x;
                }
                double boundary = ( double ) index_boundary * this->density_surface_dynein.get_X_width();
                double distance = boundary - lower_bead_position[ 0 ];
                double t_tmp = distance / tangent[ 0 ];
                x_tran_vec.push_back( t_tmp );

            }
        }

        int lower_index_Y = ( int ) ( lower_bead_position[ 1 ] / this->density_surface_dynein.get_Y_width() );
        int upper_index_Y = ( int ) ( upper_bead_position[ 1 ] / this->density_surface_dynein.get_Y_width() );
        unsigned int y_transitions = abs( upper_index_Y - lower_index_Y );

        int smer_y;
        if( tangent[ 1 ] > 0 )
        {
            smer_y = 1;
        }
        else
        {
            smer_y = - 1;
        }
        if( y_transitions > 0 )
        {

            for( int y_index = 0 ; y_index < y_transitions ; y_index ++ )
            {
                int index_boundary;
                if( smer_y * lower_bead_position[ 1 ] > 0 )
                {
                    index_boundary = lower_index_Y + ( y_index + 1 ) * smer_y;
                }
                else
                {
                    index_boundary = lower_index_Y + ( y_index ) * smer_y;
                }
                double boundary = ( double ) index_boundary * this->density_surface_dynein.get_Y_width();
                double distance = boundary - lower_bead_position[ 1 ];
                double t_tmp = distance / tangent[ 1 ];
                y_tran_vec.push_back( t_tmp );

            }
        }

        int lower_index_Z = ( int ) ( lower_bead_position[ 2 ] / this->density_surface_dynein.get_Z_width() );
        int upper_index_Z = ( int ) ( upper_bead_position[ 2 ] / this->density_surface_dynein.get_Z_width() );
        unsigned int z_transitions = abs( upper_index_Z - lower_index_Z );

        int smer_z;
        if( tangent[ 2 ] > 0 )
        {
            smer_z = 1;
        }
        else
        {
            smer_z = - 1;
        }
        if( z_transitions > 0 )
        {

            for( int z_index = 0 ; z_index < z_transitions ; z_index ++ )
            {
                int index_boundary;
                if( smer_z * lower_bead_position[ 2 ] > 0 )
                {
                    index_boundary = lower_index_Z + ( z_index + 1 ) * smer_z;
                }
                else
                {
                    index_boundary = lower_index_Z + ( z_index ) * smer_z;
                }
                double boundary = ( double ) index_boundary * this->density_surface_dynein.get_Z_width();
                double distance = boundary - lower_bead_position[ 2 ];
                double t_tmp = distance / tangent[ 2 ];
                z_tran_vec.push_back( t_tmp );

            }
        }

        final_vector.push_back( 0 );
        final_vector.insert( final_vector.end(), x_tran_vec.begin(), x_tran_vec.end() );
        final_vector.insert( final_vector.end(), y_tran_vec.begin(), y_tran_vec.end() );
        final_vector.insert( final_vector.end(), z_tran_vec.begin(), z_tran_vec.end() );

        std::vector<double> compartment_IDs;
        std::sort ( final_vector.begin(), final_vector.end() );

        double push_par = 1e-10;
        for( unsigned int parameter_index = 0 ; parameter_index < final_vector.size() ; parameter_index ++ )
        {
            double parameter = final_vector.at( parameter_index ) + push_par;
            Vector3d Vector_point = lower_bead_position + tangent * parameter;
            unsigned int compartment_id = this->density_surface_dynein.get_dynein_compartment_id_projected( Vector_point );
            if( compartment_id != 0 )
            {
                compartment_IDs.push_back( compartment_id );
            }

        }


        Vector3d axis_of_plane = lower_bead_position.cross( upper_bead_position );
        axis_of_plane = axis_of_plane / axis_of_plane.norm();
        for( unsigned int comp_index = 0 ; comp_index < compartment_IDs.size() ; comp_index ++ )
        {
            unsigned int compartment_id = compartment_IDs.at( comp_index );
            cout<<"compartment_id = "<<compartment_id<<endl;
            std::vector<Vector3d> dynein_points = this->get_dynein_in_compartment( compartment_id );
            std::vector< Vector3d > replace_motors;
            for( unsigned int bod_id  = 0 ; bod_id < dynein_points.size() ; bod_id ++ )
            {
                Vector3d motor_position = dynein_points.at( bod_id );
                if( motor_position( 2 ) < -5.0e-6 )
                {
                    replace_motors.push_back( motor_position );
                    continue;
                }

                Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
                double distance = distance_point_segment( lower_bead_position , tangent , motor_position , cl_p_of_seg );
                double distance_plane = distance_plane_point( axis_of_plane , lower_bead_position , motor_position );
                if( ( abs( distance ) < Dynein::multiplicator * Dynein::L_0 ) ) //( abs( distance_plane ) < Dynein::L_0 ) &&
                {
		    double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
                    if( probability < Dynein::attach_rate_per_step )
                    {

                        double distance_from_start_to_lower_bead = this->array_Of_Microtubules[ microtubule ].getRestDist() * ( double ) segment_id;
                        double distance_point_attachment = ( lower_bead_position - cl_p_of_seg ).norm();
                        double abscissa = ( distance_from_start_to_lower_bead + distance_point_attachment );

                        if( distance_point_attachment > this->array_Of_Microtubules[ microtubule ].getRestDist() )
                        {
                            replace_motors.push_back( motor_position );
                        }
                        else
                        {
                            std::pair < Vector3d , double > tmp_pair( motor_position , abscissa );
                            this->array_Of_Microtubules[ microtubule ].add_pair( tmp_pair );
                            this->array_Of_Microtubules[ microtubule ].set_dynein_index( 9 );
                        }
                    }
                    else
                    {
                        replace_motors.push_back( motor_position );
                    }
                }
                else
                {
                    replace_motors.push_back( motor_position );
                }
            }

            this->density_surface_dynein.set_dynein_points( compartment_id , replace_motors );
        }

    }

}

void Cell::microtubule_random_dynein_on_surface_catch_pair_abscissa_aproximation_space_line( unsigned int microtubule )
{
    std::uniform_real_distribution<> distribution{ 0 , 1 };
    unsigned int number_of_generator = omp_get_thread_num();

    for( unsigned int segment_id = 1 ; segment_id < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints()  - 1 ; segment_id ++ )
    {

        Vector3d lower_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id );
        lower_bead_position = this->project_point_on_surface( lower_bead_position );
        Vector3d upper_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id + 1 );
        upper_bead_position = this->project_point_on_surface( upper_bead_position );
        Vector3d tangent = upper_bead_position - lower_bead_position;
        Vector3d mini_segment = tangent / ( double ) Dynein::number_of_segment_steps;

        std::vector< unsigned int > compartment_IDs;
        unsigned int last_index = 0;
        for( unsigned int counter = 0 ; counter < Dynein::number_of_segment_steps ; counter ++ )
        {
            Vector3d tmp_position = lower_bead_position + ( double ) counter * mini_segment;
            unsigned int compartment_id = this->density_surface_dynein.get_dynein_compartment_id_projected( tmp_position );
            if( counter == 0 )
            {

                if( compartment_id != 0 )
                {
                    compartment_IDs.push_back( compartment_id );
                }
            }
            else
            {
                if( compartment_id != 0 )
                {
                    if( compartment_id != last_index )
                    {
                        compartment_IDs.push_back( compartment_id );
                    }
                }
            }
            last_index = compartment_id;
        }

        Vector3d axis_of_plane = lower_bead_position.cross( upper_bead_position );
        axis_of_plane = axis_of_plane / axis_of_plane.norm();
        for( unsigned int comp_index = 0 ; comp_index < compartment_IDs.size() ; comp_index ++ )
        {
            unsigned int compartment_id = compartment_IDs.at( comp_index );
            std::vector<Vector3d> dynein_points = this->get_dynein_in_compartment( compartment_id );

            std::vector< Vector3d > replace_motors;
            for( unsigned int bod_id  = 0 ; bod_id < dynein_points.size() ; bod_id ++ )
            {
                Vector3d motor_position = dynein_points.at( bod_id );
                Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
                double distance = distance_point_segment( lower_bead_position , tangent , motor_position , cl_p_of_seg );
                double distance_plane = distance_plane_point( axis_of_plane , lower_bead_position , motor_position );


                if(( abs( distance_plane ) < Dynein::L_0 ) &&  ( abs( distance ) < Dynein::multiplicator * Dynein::L_0 ) ) //
                {
		    double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
                    if( probability < Dynein::attach_rate_per_step )
                    {


                        double distance_from_start_to_lower_bead = this->array_Of_Microtubules[ microtubule ].getRestDist() * ( double ) segment_id;
                        double distance_point_attachment = ( lower_bead_position - cl_p_of_seg ).norm();

                        double abscissa = ( distance_from_start_to_lower_bead + distance_point_attachment ) - Dynein::reduction_abscissa;

                        unsigned int index_lower_bead = (unsigned int) ( abscissa / sim_of_Cell::resting_distance );
                        Vector3d lower_point = this->array_Of_Microtubules[ microtubule ].getPoint( index_lower_bead );
                        Vector3d tangent = this->array_Of_Microtubules[ microtubule ].getTangent2( index_lower_bead );

                        double abscissa_minus_lower_bead = abscissa - ( ( double ) index_lower_bead ) * sim_of_Cell::resting_distance ;
                        Vector3d point_of_attachment = lower_bead_position + distance_point_attachment * tangent / tangent.norm();

                        Vector3d ancher_atachment = motor_position - point_of_attachment;
                        double distance_ancher_atachment = ancher_atachment.norm();



                        if( distance_point_attachment > this->array_Of_Microtubules[ microtubule ].getRestDist() )
                        {
                            replace_motors.push_back( motor_position );
                        }
                        else
                        {

                            if( segment_id == 9 )
                            {

                            }
                            std::pair < Vector3d , double > tmp_pair( motor_position , abscissa );
                            this->array_Of_Microtubules[ microtubule ].add_pair( tmp_pair );
                            this->array_Of_Microtubules[ microtubule ].set_dynein_index( 9 );
                        }
                    }
                    else
                    {
                        replace_motors.push_back( motor_position );
                    }
                }
                else
                {
                    replace_motors.push_back( motor_position );
                }
            }

            this->density_surface_dynein.set_dynein_points( compartment_id , replace_motors );
        }
    }


}




void Cell::microtubule_random_dynein_on_surface_catch_pair_abscissa_aproximation_space_line_2( unsigned int microtubule )
{

		    std::uniform_real_distribution<> distribution{ 0 , 1 };
		    unsigned int number_of_generator = omp_get_thread_num();


    for( unsigned int segment_id = 1 ; segment_id < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints()  - 1 ; segment_id ++ )
    {

        Vector3d lower_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id );
        Vector3d lower_bead_position_2 = this->project_point_on_surface( lower_bead_position );
        Vector3d upper_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id + 1 );
        Vector3d upper_bead_position_2 = this->project_point_on_surface( upper_bead_position );
        Vector3d tangent = this->array_Of_Microtubules[ microtubule ].getTangent2( segment_id );
        Vector3d tangent_2 = upper_bead_position_2 - lower_bead_position_2;
        Vector3d mini_segment = tangent_2 / ( double ) Dynein::number_of_segment_steps;

        std::vector< unsigned int > compartment_IDs;
        unsigned int last_index = 0;
        for( unsigned int counter = 0 ; counter < Dynein::number_of_segment_steps ; counter ++ )
        {
            Vector3d tmp_position = lower_bead_position_2 + ( double ) counter * mini_segment;
            unsigned int compartment_id = this->density_surface_dynein.get_dynein_compartment_id_projected( tmp_position );
            if( counter == 0 )
            {

                if( compartment_id != 0 )
                {
                    compartment_IDs.push_back( compartment_id );
                }
            }
            else
            {
                if( compartment_id != 0 )
                {
                    if( compartment_id != last_index )
                    {
                        compartment_IDs.push_back( compartment_id );
                    }
                }
            }
            last_index = compartment_id;
        }

        Vector3d axis_of_plane = lower_bead_position.cross( upper_bead_position );
        axis_of_plane = axis_of_plane / axis_of_plane.norm();
        for( unsigned int comp_index = 0 ; comp_index < compartment_IDs.size() ; comp_index ++ )
        {
            unsigned int compartment_id = compartment_IDs.at( comp_index );

            std::vector<Vector3d> dynein_points = this->get_dynein_in_compartment( compartment_id );

            std::vector< Vector3d > replace_motors;
            for( unsigned int bod_id  = 0 ; bod_id < dynein_points.size() ; bod_id ++ )
            {
                Vector3d motor_position = dynein_points.at( bod_id );

                Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
                double distance = distance_point_segment( lower_bead_position , tangent , motor_position , cl_p_of_seg );
                double distance_plane = distance_plane_point( axis_of_plane , lower_bead_position , motor_position );


                if(( abs( distance_plane ) < Dynein::L_0 ) &&  ( abs( distance ) < Dynein::multiplicator * Dynein::L_0 ) ) //
                {
		    double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
                    if( probability < Dynein::attach_rate_per_step )
                    {

                        double distance_lower_bead = this->array_Of_Microtubules[ microtubule ].getRestDist() * ( double ) segment_id;
                        distance_lower_bead = 0;
                        for( unsigned int i = 0 ; i < segment_id ; i ++ )
                        {
                            distance_lower_bead = distance_lower_bead + this->array_Of_Microtubules[ microtubule ].getTangent2( i ).norm();
                        }

                        double distance_point_attachment = ( lower_bead_position - cl_p_of_seg ).norm();
                        double abscissa = ( distance_lower_bead + distance_point_attachment ) ;//- Dynein::reduction_abscissa
                        double lenght_micro = this->array_Of_Microtubules[ microtubule ].getRestDist() *
                        ( double )( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 );
                        if( abscissa >= lenght_micro - 8.0e-9 )
                        {
                            replace_motors.push_back( motor_position );
                            continue;
                        }
                        else
                        {

                        Vector3d new_motor_position;
                        if( cl_p_of_seg.norm() > this->a_axis - Dynein::L_0 )
                        {
                            new_motor_position = cl_p_of_seg;
                        }
                        else
                        {
                            Vector3d perpendicular_tangent = ( motor_position - cl_p_of_seg );

                            perpendicular_tangent = perpendicular_tangent / perpendicular_tangent.norm();
                            new_motor_position = cl_p_of_seg + perpendicular_tangent * Dynein::L_0; // /* 3.0 / 4.0  + perpendicular_tangent * Dynein::L_0
                        }
                        std::pair < Vector3d , double > tmp_pair( new_motor_position , abscissa );
                        this->array_Of_Microtubules[ microtubule ].add_pair( tmp_pair );
                        this->array_Of_Microtubules[ microtubule ].set_dynein_index( 9 );
                        }
                    }
                    else
                    {
                        replace_motors.push_back( motor_position );
                    }
                }
                else
                {
                    replace_motors.push_back( motor_position );
                }
            }

            this->density_surface_dynein.set_dynein_points( compartment_id , replace_motors );
        }



    }


}





void Cell::microtubule_random_dynein_on_surface_catch_pair_abscissa_aproximation_true( unsigned int microtubule )
{
    for( unsigned int segment_id = 1 ; segment_id < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints()  ; segment_id ++ ) // - 1
    {

        Vector3d lower_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id );
        cout<<lower_bead_position.norm()<<endl;
    }


}


void Cell::all_microtubules_random_dynein_on_surface_catch_pair_abscissa_2()
{
    for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
    {

        if( this->array_Of_Microtubules[ microtubule ].get_polygon_number() >= IS_Dynein_Cell_surface::number_of_polygon_higher )
        {
            continue;
        }
        if( this->array_Of_Microtubules[ microtubule ].get_polygon_number() < IS_Dynein_Cell_surface::number_of_polygon_lower )
        {
            continue;
        }
        if( this->array_Of_Microtubules[ microtubule ].getID() >= IS_Dynein_Cell_surface::number_of_mito )
        {
            continue;
        }
        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 0 )
        {
            this->microtubule_random_dynein_on_surface_catch_pair_abscissa_aproximation_space_line_2( microtubule );
        }
        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9 )
        {
            this->microtubule_random_dynein_on_surface_catch_pair_abscissa_aproximation_space_line_2( microtubule );
        }
    }

}





Vector3d Cell::MTOC_dynein_capture_shrinkage_force( Vector3d MTOC_bead , Vector3d Microtubule_catching )
{
	Vector3d orientace = Microtubule_catching - MTOC_bead;
	orientace = orientace / orientace.norm();
	Vector3d force = orientace * IS_Capture_shrinkage_param::dynein_Force_capture_shrinkage;
	return force;

}

Vector3d Cell::MTOC_dynein_capture_shrinkage_force( unsigned int microtubule )
{
    unsigned int MTOC_point_index = this->array_Of_Microtubules[ microtubule ].get_MTOC_point();
    Vector3d MTOC_point = this->MTOC.get_point( MTOC_point_index );
    Vector3d Microtubule_catching = this->array_Of_Microtubules[ microtubule ].get_IS_position_catching();

   	Vector3d orientace = Microtubule_catching - MTOC_point;
	orientace = orientace / orientace.norm();
	Vector3d force = orientace * IS_Capture_shrinkage_param::dynein_Force_capture_shrinkage;
	return force;

}


Vector3d Cell::simple_dynein_capture_shrinkage_force( unsigned int microtubule )
{
    Vector3d MTOC_point = this->mtoc.get_position();
    Vector3d Microtubule_catching = this->array_Of_Microtubules[ microtubule ].get_IS_position_catching();

   	Vector3d orientace = Microtubule_catching - MTOC_point;
	orientace = orientace / orientace.norm();
	Vector3d force = orientace * IS_Capture_shrinkage_param::dynein_Force_capture_shrinkage;
	return force;

}





Vector3d Cell::MTOC_dynein_cortical_sliding_force( Vector3d MTOC_bead ,  Vector3d Microtubule_catching  )
{
	Vector3d orientace = Microtubule_catching - MTOC_bead;
	orientace = orientace / orientace.norm();
	Vector3d force = orientace * IS_Capture_shrinkage_param::dynein_Force_capture_shrinkage;
	return force;
}


void Cell::MTOC_naive_circular_force(  MatrixXd &force_MTOC )
{
    if( force_MTOC.rows() != 3 * ( this->MTOC.get_number_of_points() + 1 ) )
    {
        cout<<"force_MTOC.rows() != this->MTOC.get_number_of_points() * 3"<<endl;
        cout<<"Cell::MTOC_naive_circular_force( MatrixXd &force_MTOC )"<<endl;
        cout<<"ERROR_ID Cell 4646116911546945"<<endl;
        throw("");
    }

    Vector3d center = this->MTOC.get_center();


    Vector3d orientation( 0.0 , 0.0 , 0.0 );
    orientation( 0 ) = center( 2 );
    orientation( 2 ) = - center( 0 );
    orientation = orientation / orientation.norm();

    double overall_force = IS_Capture_shrinkage_param::dynein_Force_capture_shrinkage * 5.0;
    double force_one_bead = overall_force / (double) (  this->MTOC.get_number_of_points() + 1 );
    Vector3d force = force_one_bead * orientation;

    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        force_MTOC( dimension , 0 ) = force_MTOC( dimension , 0 ) + force( dimension );
    }
    for( unsigned int point = 0 ; point <= this->MTOC.get_number_of_points() ; point )
    {
        for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
        {
            force_MTOC( 3 * point + dimension , 0 ) = force_MTOC(  3 * point + dimension , 0 ) + force( dimension );
        }
    }




}

void Cell::MTOC_naive_force_to_center(  MatrixXd &force_MTOC )
{
    if( force_MTOC.rows() != 3 * ( this->MTOC.get_number_of_points() + 1 ) )
    {
        cout<<"force_MTOC.rows() != this->MTOC.get_number_of_points() * 3"<<endl;
        cout<<"Cell::MTOC_naive_force_to_center( MatrixXd &force_MTOC )"<<endl;
        cout<<"ERROR_ID Cell 4646116911546945"<<endl;
        throw("");
    }

    Vector3d center = this->MTOC.get_center();

    Vector3d orientation = ( -1.0 ) * center;
    orientation = orientation / orientation.norm();

    double overall_force = 0; 
    double force_one_bead = overall_force / (double) (  this->MTOC.get_number_of_points() + 1 );
    Vector3d force = force_one_bead * orientation;

    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        force_MTOC( dimension , 0 ) = force_MTOC( dimension , 0 ) + force( dimension );
    }
    for( unsigned int point = 0 ; point <= this->MTOC.get_number_of_points() ; point ++  )
    {
        for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
        {
            force_MTOC( 3 * point + dimension , 0 ) = force_MTOC(  3 * point + dimension , 0 ) + force( dimension );
        }
    }

}





Vector3d Cell::project_point_on_surface( Vector3d position )
{
    double A_multiplicator = Cell_parametres::A_AXIS * Cell_parametres::A_AXIS;
    double B_multiplicator;
    double numerator;
    if( position( 2 ) > 0 )
    {
        B_multiplicator = Cell_parametres::B_AXIS * Cell_parametres::B_AXIS;
        numerator = Cell_parametres::A_AXIS * Cell_parametres::B_AXIS;
    }
    else
    {
        B_multiplicator = Cell_parametres::B_AXIS_Lower * Cell_parametres::B_AXIS_Lower;
        numerator = Cell_parametres::A_AXIS * Cell_parametres::B_AXIS_Lower;
    }

    double denominator = ( position( 0 ) * position( 0 ) + position( 1 ) * position( 1 ) ) * B_multiplicator;
    denominator = denominator + position( 2 ) * position( 2 ) * A_multiplicator;
    denominator = sqrt( denominator );
    double c_multiplicator = numerator / denominator;

    Vector3d new_position = position * c_multiplicator;
    return new_position;
}


unsigned int Cell::get_dynein_compartment_id( Vector3d position )
{

    Vector3d add_to_get_segment( Cell_parametres::A_AXIS , Cell_parametres::A_AXIS , Cell_parametres::B_AXIS_Lower );
    Vector3d position_to_get_index = position + add_to_get_segment;

    unsigned int x_index = position_to_get_index( 0 ) / this->density_surface_dynein.get_X_width();
    unsigned int y_index = position_to_get_index( 1 ) / this->density_surface_dynein.get_Y_width();
    unsigned int z_index = position_to_get_index( 2 ) / this->density_surface_dynein.get_Z_width();

    unsigned int neccessary_dimension = this->density_surface_dynein.get_neccessary_dimension();
    unsigned int ID_map_index = x_index * neccessary_dimension * neccessary_dimension + y_index * neccessary_dimension + z_index;
    return ID_map_index;
}





void Cell::project_and_add_points_to_surface( std::vector< Vector3d > points )
{
    for( unsigned int point_id = 0 ; point_id < points.size() ; point_id ++ )
    {
        Vector3d point_position = points.at( point_id );
        Vector3d point_position_new = this->project_point_on_surface( point_position );
        unsigned int map_key = get_dynein_compartment_id( point_position_new );
        this->density_surface_dynein.add_dynein_point( map_key , point_position_new );
    }

}

void  Cell::add_points_to_surface( std::vector< Vector3d > points )
{
    for( unsigned int point_id = 0 ; point_id < points.size() ; point_id ++ )
    {
        Vector3d point_position = points.at( point_id );
        unsigned int map_key = get_dynein_compartment_id( point_position );
        this->density_surface_dynein.add_dynein_point( map_key , point_position );
    }
}




Vector3d Cell::force_MTOC_bending_first_point( unsigned int microtubule )
{

    if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() < 2 )
    {
        cout<<"this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() < 2"<<endl;
        cout<<" force_MTOC_bending_first_point( unsigned int microtubule ) "<<endl;
        cout<<"ERROR_ID = 57453496468448"<<endl;
        throw("");
    }

    if( microtubule >= this->get_microtubule_number()  )
    {
        cout<<" microtubule >= this->get_microtubule_number()  "<<endl;
        cout<<" force_MTOC_bending_first_point( unsigned int microtubule ) "<<endl;
        cout<<"ERROR_ID = 6846846316834685"<<endl;
        throw("");
    }


    double circle_radius = this->MTOC.get_radius() + sim_of_Cell::resting_distance;
    double sphere_surface = 2.0 * sim_of_Cell::PI * circle_radius;
    double equilibrium_distance = sphere_surface / ( double ) this->get_microtubule_number();

    Vector3d force( 0.0 , 0.0 , 0.0 );
    Vector3d position_first_point = this->array_Of_Microtubules[ microtubule ].getPoint( 1 );
    for( unsigned int micro_count = 0 ; micro_count < this->get_microtubule_number() ; micro_count ++ )
    {
        if( micro_count == microtubule )
        {
            continue;

        }
        if( this->array_Of_Microtubules[ micro_count ].getNumberOfPoints() < 2 )
        {
            continue;
        }


        Vector3d position_tmp = this->array_Of_Microtubules[ micro_count ].getPoint( 1 );
        double distance = ( position_first_point - position_tmp ).norm();

        if( distance == 0 )
        {
        }

        if( ( distance > equilibrium_distance ) )
        {

            continue;
        }
        else
        {

            double constant_2 = 1.20;
            double distance_difference = abs( equilibrium_distance - distance );
            double force_abs_value = 1.0 * Cell_parametres::wall_cell_k1 * ( exp( constant_2 * distance_difference ) - 1.0 );

            Vector3d orientation = position_first_point - position_tmp;

            orientation = orientation / orientation.norm();
            Vector3d force_tmp = orientation * force_abs_value;
            force = force + force_tmp;
        }

    }
    return force;
}




Vector3d Cell::force_MTOC_bending_first_point_2( unsigned int microtubule , Vector3d& force_2 )
{

    if( this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() < 2 )
    {
        cout<<"this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() < 2"<<endl;
        cout<<" force_MTOC_bending_first_point_2( unsigned int microtubule ) "<<endl;
        cout<<"ERROR_ID = 57453496468448"<<endl;
        throw("");
    }

    if( microtubule >= this->get_microtubule_number()  )
    {
        cout<<" microtubule >= this->get_microtubule_number()  "<<endl;
        cout<<" force_MTOC_bending_first_point_2( unsigned int microtubule ) "<<endl;
        cout<<"ERROR_ID = 6846846316834685"<<endl;
        throw("");
    }

    Vector3d force( 0.0 , 0.0 , 0.0 );
    Vector3d position_first_point = this->array_Of_Microtubules[ microtubule ].getPoint( 1 );

    Vector3d force_2_tmp( 0.0 , 0.0 , 0.0 );
    Vector3d position_second_point( 0.0 , 0.0 , 0.0 );
    if(  this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() > 2 )
    {
        position_second_point = this->array_Of_Microtubules[ microtubule ].getPoint( 2 );
    }



    for( unsigned int micro_count = 0 ; micro_count < this->get_microtubule_number() ; micro_count ++ )
    {
        if( micro_count == microtubule )
        {
            continue;
        }
        if( this->array_Of_Microtubules[ micro_count ].getNumberOfPoints() < 2 )
        {
            continue;
        }

        double constant_2 = 0;// 
        double constant_22 = 0;//
        Vector3d position_tmp = this->array_Of_Microtubules[ micro_count ].getPoint( 1 );
        double distance = ( position_first_point - position_tmp ).norm();

        if( distance == 0 )
        {
            continue;
        }
        else if( distance < 5e-8)
        {
            double force_abs_value = constant_2 / ( 5e-8 * 5e-8 );
            Vector3d orientation = position_first_point - position_tmp;
            orientation = orientation / orientation.norm();
            Vector3d force_tmp = orientation * force_abs_value;
            force = force + force_tmp;
        }
        else
        {
            double force_abs_value = constant_2 / ( distance * distance );
            Vector3d orientation = position_first_point - position_tmp;
            orientation = orientation / orientation.norm();
            Vector3d force_tmp = orientation * force_abs_value;
            force = force + force_tmp;
        }

        if(  this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() > 2 )
        {
            if( this->array_Of_Microtubules[ micro_count ].getNumberOfPoints() < 3 )
            {
                continue;
            }
            Vector3d position_tmp_2 = this->array_Of_Microtubules[ micro_count ].getPoint( 2 );
            double distance_2 = ( position_second_point - position_tmp_2 ).norm();
            if( distance == 0 )
            {
                continue;
            }
            else if( distance < 5e-8)
            {
                double force_abs_value_2 = constant_22 / ( 5e-8 * 5e-8 );
                Vector3d orientation_2 = position_second_point - position_tmp_2;
                orientation_2 = orientation_2 / orientation_2.norm();
                Vector3d ll = orientation_2 * force_abs_value_2;
                force_2_tmp = force_2_tmp + ll;
            }
            else
            {
                double force_abs_value_2 = constant_22 / ( distance_2 * distance_2 );
                Vector3d orientation_2 = position_second_point - position_tmp_2;
                orientation_2 = orientation_2 / orientation_2.norm();
                Vector3d ll = orientation_2 * force_abs_value_2;
                force_2_tmp = force_2_tmp + ll;
            }
        }
    }
    force_2 = force_2_tmp;
    return force;
}





unsigned int Cell::get_opposite_micro_number( unsigned int microtubule )
{
    unsigned int this_micro_polygon_number = this->array_Of_Microtubules[ microtubule ].get_polygon_number();
    unsigned int this_micro_id = this->array_Of_Microtubules[ microtubule ].getID();

    unsigned int number_of_mito_without_extra = this->number_of_microtubules - this->number_of_microtubules_extra;
    unsigned int number_of_poly_in_cell = number_of_mito_without_extra / MTOCparam::micro_in_polygon;
    unsigned int half_number_of_poly_in_cell = number_of_poly_in_cell / 2;

    unsigned int second_micro_poly_number = this_micro_polygon_number + half_number_of_poly_in_cell;
    unsigned int second_micro_number = second_micro_poly_number * MTOCparam::micro_in_polygon + this_micro_id;
    return second_micro_number;

}






void Cell::bending_micro_force( unsigned int micro_id , Vector3d& first_micro , Vector3d& second_micro , Vector3d& mtoc_center )
{

    unsigned int number_of_regular_micro = this->get_microtubule_number() - this->number_of_microtubules_extra;
    if( micro_id >= number_of_regular_micro / 2 )
    {
        cout<<"  micro_id >= number_of_regular_micro  / 2  "<<endl;
        cout<<"bending_micro_force( unsigned int micro_id , Vector3d first_micro , Vector3d second_micro , Vector3d mtoc_center )"<<endl;
        cout<<"ERROR_ID = 896794894739489484"<<endl;
        throw("");
    }
    unsigned int second_microtubule_number = this->get_opposite_micro_number( micro_id );

    Vector3d first_micro_point = this->array_Of_Microtubules[ micro_id ].getPoint( 1 );
    Vector3d second_micro_point = this->array_Of_Microtubules[ second_microtubule_number ].getPoint( 1 );
    Vector3d MTOC_center = this->MTOC.get_center();

    Vector3d first = MTOC_center - first_micro_point;
    first = first / first .norm();
    Vector3d second = second_micro_point - MTOC_center;
    second = second / second.norm();

    double dot_product = first.dot( second );


    first_micro( 0 ) = ( - second( 0 ) + first( 0 ) * dot_product ) / sim_of_Cell::resting_distance;
	first_micro( 1 ) = ( - second( 1 ) + first( 1 ) * dot_product ) / sim_of_Cell::resting_distance;
	first_micro( 2 ) = ( - second( 2 ) + first( 2 ) * dot_product ) / sim_of_Cell::resting_distance;

    second_micro( 0 ) = ( first( 0 ) - second( 0 ) * dot_product ) / sim_of_Cell::resting_distance;
	second_micro( 1 ) = ( first( 1 ) - second( 1 ) * dot_product ) / sim_of_Cell::resting_distance;
	second_micro( 2 ) = ( first( 2 ) - second( 2 ) * dot_product ) / sim_of_Cell::resting_distance;

    mtoc_center( 0 ) = - first_micro( 0 ) - second_micro( 0 );
    mtoc_center( 1 ) = - first_micro( 1 ) - second_micro( 1 );
    mtoc_center( 2 ) = - first_micro( 2 ) - second_micro( 2 );

    first_micro = first_micro * sim_of_Cell::k_bending;
    second_micro = second_micro * sim_of_Cell::k_bending;
    mtoc_center = mtoc_center * sim_of_Cell::k_bending;
}



void Cell::bending_abstract_mtoc( unsigned int micro_id ,  Vector3d& first_micro , Vector3d& second_micro , Vector3d& mtoc_center )
{

    if( micro_id >= this->get_microtubule_number() )
    {
        cout<<"micro_id >= this->get_microtubule_number()"<<endl;
        cout<<"ERROR_ID = 9877306816846348"<<endl;
        cout<<"bending_abstract_mtoc( unsigned int micro_id ,  Vector3d& first_micro , Vector3d& second_micro , Vector3d& mtoc_center )"<<endl;
        throw("");
    }

    if( this->array_Of_Microtubules[ micro_id ].getNumberOfPoints() < 2 )
    {
        return;
    }

    Vector3d first_tangent = this->array_Of_Microtubules[ micro_id ].getTangent2( 0 );
    Vector3d first_tangent_MTOC = this->abstract_MTOC.get_tangent( micro_id );

    Vector3d first = first_tangent_MTOC;
    first = first / first .norm();
    Vector3d second = first_tangent;
    second = second / second.norm();

    double dot_product = first.dot( second );


    first_micro( 0 ) = ( first( 0 ) - second( 0 ) * dot_product ) / sim_of_Cell::resting_distance;
	first_micro( 1 ) = ( first( 1 ) - second( 1 ) * dot_product ) / sim_of_Cell::resting_distance;
	first_micro( 2 ) = ( first( 2 ) - second( 2 ) * dot_product ) / sim_of_Cell::resting_distance;
    mtoc_center = - first_micro;// - second_micro


    first_micro = first_micro * sim_of_Cell::k_bending;
    second_micro = second_micro * sim_of_Cell::k_bending;
    mtoc_center = mtoc_center * sim_of_Cell::k_bending;

}


void Cell::resize_micro_MTOC(  )
{
    Vector3d center_real_MTOC = this->MTOC.get_center();

    MatrixXd original_orientation = this->abstract_MTOC.get_orientations();
    MatrixXd coordinates_arg = MatrixXd::Zero( 3 * ( this->get_microtubule_number() + 1 ) , 1 );
    this->abstract_MTOC.rotate_orientations( center_real_MTOC , coordinates_arg );

    MatrixXd original_mtoc_coordinates = this->abstract_MTOC.get_mtoc_coordinates();
    MatrixXd mtoc_rotated = MatrixXd::Zero( 3 * ( this->get_microtubule_number() + 1 ) , 1 );
    this->abstract_MTOC.rotate_MTOC( center_real_MTOC , mtoc_rotated );

    this->abstract_MTOC.set_orientations( coordinates_arg );
    this->abstract_MTOC.set_mtoc_coordinates( mtoc_rotated );



    for( unsigned int micro = 0 ; micro < this->number_of_microtubules ; micro ++ )
    {
        Vector3d orientation_tmp = this->abstract_MTOC.get_original_micro_orientation( micro );
        orientation_tmp = orientation_tmp / orientation_tmp.norm();
        this->array_Of_Microtubules[ micro ].resizeMicrotubule( orientation_tmp );
    }

    this->abstract_MTOC.set_orientations( original_orientation );
    this->abstract_MTOC.set_mtoc_coordinates( original_mtoc_coordinates );


}


void Cell::resize_micro_with_different_first_segment_and_MTOC(  )
{
    this->MTOC.resize_from_originals();

    for( unsigned int micro_index = 0 ; micro_index < this->number_of_microtubules ; micro_index ++ )
    {
        if( this->array_Of_Microtubules[ micro_index ].get_dynein_index() == 20  )
        {
            this->array_Of_Microtubules[ micro_index ].resizeMicrotubule_with_different_tangent_lenghts();
        }
        else
        {
            this->array_Of_Microtubules[ micro_index ].resizeMicrotubule_with_different_tangent_lenghts();
        }
    }

}


void Cell::ultimate_resizing()
{
    MatrixXd new_coordinates = this->MTOC.rotate_MTOC_from_original_orientations();
    this->MTOC.set_coordinates( new_coordinates );
}


void Cell::catch_pair_abscissa()
{
        for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
    {

        if( this->array_Of_Microtubules[ microtubule ].get_polygon_number() >= IS_Dynein_Cell_surface::number_of_polygon_higher )
        {
            continue;
        }
        if( this->array_Of_Microtubules[ microtubule ].get_polygon_number() < IS_Dynein_Cell_surface::number_of_polygon_lower )
        {
            continue;
        }
        if( this->array_Of_Microtubules[ microtubule ].getID() >= IS_Dynein_Cell_surface::number_of_mito )
        {
            continue;
        }
        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 0 )
        {
            this->microtubule_catch_pair_abscissa( microtubule );
        }
        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9 )
        {
            this->microtubule_catch_pair_abscissa( microtubule );
        }

    }
}

void Cell::catch_pair_abscissa_real_dynein()
{
    for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules ; microtubule ++ )
    {
        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 0 )
        {
            this->microtubule_catch_pair_abscissa_real_dynein( microtubule );
        }
        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9 )
        {
            this->microtubule_catch_pair_abscissa_real_dynein( microtubule );
        }

    }
}






void Cell::microtubule_catch_pair_abscissa( unsigned int microtubule )
{
    std::uniform_real_distribution<> distribution{ 0 , 1 };
    unsigned int number_of_generator = omp_get_thread_num();

    double lenght_micro = 0.0;
    for( unsigned int segment_id = 0 ; segment_id < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints() - 1 ; segment_id ++ )
    {
        lenght_micro = lenght_micro + this->array_Of_Microtubules[ microtubule ].getTangent2( segment_id ).norm();
    }


    for( unsigned int segment_id = 1 ; segment_id < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints()  - 1 ; segment_id ++ )
    {

        Vector3d lower_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id );
        Vector3d lower_bead_position_2 = this->project_point_on_surface( lower_bead_position );
        Vector3d upper_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id + 1 );
        Vector3d upper_bead_position_2 = this->project_point_on_surface( upper_bead_position );
        Vector3d tangent = this->array_Of_Microtubules[ microtubule ].getTangent2( segment_id );
        Vector3d tangent_2 = upper_bead_position_2 - lower_bead_position_2;
        Vector3d mini_segment = tangent_2 / ( double ) Dynein::number_of_segment_steps;

        std::vector< unsigned int > compartment_IDs;
        unsigned int last_index = 0;
        for( unsigned int counter = 0 ; counter < Dynein::number_of_segment_steps ; counter ++ )
        {
            Vector3d tmp_position = lower_bead_position_2 + ( double ) counter * mini_segment;
            unsigned int compartment_id = this->density_surface_dynein.get_dynein_compartment_id_projected( tmp_position );
            if( counter == 0 )
            {

                if( compartment_id != 0 )
                {
                    compartment_IDs.push_back( compartment_id );
                }
            }
            else
            {
                if( compartment_id != 0 )
                {
                    if( compartment_id != last_index )
                    {
                        compartment_IDs.push_back( compartment_id );
                    }
                }
            }
            last_index = compartment_id;
        }

        Vector3d axis_of_plane = lower_bead_position.cross( upper_bead_position );
        axis_of_plane = axis_of_plane / axis_of_plane.norm();
        for( unsigned int comp_index = 0 ; comp_index < compartment_IDs.size() ; comp_index ++ )
        {
            unsigned int compartment_id = compartment_IDs.at( comp_index );

            std::vector<Vector3d> dynein_points = this->get_dynein_in_compartment( compartment_id );
            std::vector< Vector3d > replace_motors;
            for( unsigned int bod_id  = 0 ; bod_id < dynein_points.size() ; bod_id ++ )
            {
                Vector3d motor_position = dynein_points.at( bod_id );
                Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
                double distance = distance_point_segment( lower_bead_position , tangent , motor_position , cl_p_of_seg );
                double distance_plane = distance_plane_point( axis_of_plane , lower_bead_position , motor_position );


                if(( abs( distance_plane ) < Dynein::L_0 ) &&  ( abs( distance ) < Dynein::multiplicator * Dynein::L_0 ) ) //
                {

		    double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );


                    if( probability < Dynein::attach_rate_per_step )
                    {
                        double distance_lower_bead = this->array_Of_Microtubules[ microtubule ].getRestDist() * ( double ) segment_id;
                        distance_lower_bead = 0;
                        for( unsigned int i = 0 ; i < segment_id ; i ++ )
                        {
                            distance_lower_bead = distance_lower_bead + this->array_Of_Microtubules[ microtubule ].getTangent2( i ).norm();
                        }

                        double distance_point_attachment = ( lower_bead_position - cl_p_of_seg ).norm();
                        double abscissa = ( distance_lower_bead + distance_point_attachment ) ;//- Dynein::reduction_abscissa

                        if( abscissa >= lenght_micro - 8.0e-9 )
                        {
                            replace_motors.push_back( motor_position );
                            continue;
                        }
                        else
                        {

                            Vector3d new_motor_position;
                            if( cl_p_of_seg.norm() > this->a_axis - Dynein::L_0 )
                            {
                                new_motor_position = cl_p_of_seg;
                            }
                            else
                            {
                                Vector3d perpendicular_tangent = ( motor_position - cl_p_of_seg );
                                perpendicular_tangent = perpendicular_tangent / perpendicular_tangent.norm();
                                new_motor_position = cl_p_of_seg + perpendicular_tangent * Dynein::L_0; // /* 3.0 / 4.0  + perpendicular_tangent * Dynein::L_0
                            }
                            std::pair < Vector3d , double > tmp_pair( new_motor_position , abscissa );
                            this->array_Of_Microtubules[ microtubule ].add_pair( tmp_pair );
                            this->array_Of_Microtubules[ microtubule ].set_dynein_index( 9 );
                        }

                    }
                    else
                    {
                        replace_motors.push_back( motor_position );
                    }
                }
                else
                {
                    replace_motors.push_back( motor_position );
                }
            }

            this->density_surface_dynein.set_dynein_points( compartment_id , replace_motors );
        }



    }



}


std::vector<Vector3d> Cell::get_dynein_in_compartment( unsigned int ID_map_index )
{

    std::vector<Vector3d> returned_points;
    std::vector<Vector3d> dynein_on_surface_points = this->density_surface_dynein.get_dynein_points( ID_map_index );
    returned_points = dynein_on_surface_points;
    return returned_points;
}

unsigned int Cell::get_number_of_dynein_motors()
{
    unsigned int number_of_motors = this->density_surface_dynein.get_dynein_motors_number();
    unsigned int sum_of_motor_in_micro = 0;
    for( unsigned int motor_id = 0 ; motor_id < this->number_of_microtubules; motor_id ++ )
    {
        if( this->array_Of_Microtubules[ motor_id ].get_dynein_index() == 9 )
        {
            unsigned int number_of_dynein_in_micro = this->array_Of_Microtubules[ motor_id ].get_Dynein_abscissa().size();
            sum_of_motor_in_micro = sum_of_motor_in_micro + number_of_dynein_in_micro;
        }
    }

    return number_of_motors + sum_of_motor_in_micro;
}











unsigned int Cell::get_number_of_dynein_motors_cortical_sliding_micro()
{

    unsigned int sum_of_motor_in_micro = 0;
    for( unsigned int motor_id = 0 ; motor_id < this->number_of_microtubules; motor_id ++ )
    {
        if( this->array_Of_Microtubules[ motor_id ].get_dynein_index() == 9 )
        {
            unsigned int number_of_dynein_in_micro = this->array_Of_Microtubules[ motor_id ].get_Dynein_abscissa().size();
            sum_of_motor_in_micro = sum_of_motor_in_micro + number_of_dynein_in_micro;
        }
    }
    return sum_of_motor_in_micro;
}














void Cell::get_attached_dynein_motors_cortical_sliding_micro_both_surface_IS( double& concentration_in_IS ,  double& concentration_outside_IS )
{

    unsigned int number_in_IS_tmp = 0;
    unsigned int number_outside_IS_tmp = 0;
    for( unsigned int microtubule = 0 ; microtubule < this->number_of_microtubules; microtubule ++ )
    {
        if( this->array_Of_Microtubules[ microtubule ].get_dynein_index() == 9 )
        {
            std::vector<Vector3d> microtubule_anchor_points = this->array_Of_Microtubules[ microtubule ].get_catching_points_IS_sliding( );
	    unsigned int number_in = 0;
	    unsigned int number_out = 0;


	    if( microtubule_anchor_points.size() > 0 )
	    {
	    	this->IS_Cortic_Sl.control_IS_Cortical_Sliding_points3( microtubule_anchor_points , number_in , number_out );
	    }

	    number_in_IS_tmp = number_in_IS_tmp + number_in;
	    number_outside_IS_tmp = number_outside_IS_tmp + number_out;
        }
    }

    double micro_2_to_meter_2 = 1e12;
    double virtual_sphere_radius = ( this->a_axis + this->b_axis ) / 2.0;

    double surface = 4.0 * sim_of_Cell::PI * virtual_sphere_radius * virtual_sphere_radius * micro_2_to_meter_2;
    double number_of_dynein_surface = ( surface * sim_of_Cell::density_of_dynein_motor_surface );


    double area_of_IS_cor_sl = sim_of_Cell::PI * IS_Cortical_Sl_parameter::radius * IS_Cortical_Sl_parameter::radius * micro_2_to_meter_2;
    double number_of_dynein_IS = ( area_of_IS_cor_sl * sim_of_Cell::density_of_dynein_motor_Cortical_Sliding  );


   double density_outside = 0;
   if(  number_of_dynein_surface > number_outside_IS_tmp )
   {
        density_outside =   ( number_of_dynein_surface - ( double ) number_outside_IS_tmp ) / surface;
   }


   double density_inside = 0;
   if(  number_of_dynein_IS > number_in_IS_tmp )
   {
   	density_inside =   ( number_of_dynein_IS - ( double ) number_in_IS_tmp ) / area_of_IS_cor_sl;
   }

   concentration_in_IS = density_inside;
   concentration_outside_IS = density_outside;


}









unsigned int Cell::get_number_of_dynein_motors_in_micro_capture()
{
    unsigned int number_of_motors = this->Capture_Shrinkage_dynein.get_dynein_motors_number();

    unsigned int sum_of_motor_in_micro = 0;
    unsigned int number_of_20_micro = 0;
    for( unsigned int motor_id = 0 ; motor_id < this->number_of_microtubules; motor_id ++ )
    {
        if( this->array_Of_Microtubules[ motor_id ].get_dynein_index() == 20 )
        {
            unsigned int number_of_dynein_in_micro = this->array_Of_Microtubules[ motor_id ].get_Dynein_abscissa().size();
            sum_of_motor_in_micro = sum_of_motor_in_micro + number_of_dynein_in_micro;
	    number_of_20_micro = number_of_20_micro + 1;
        }
    }

    return sum_of_motor_in_micro;
}




unsigned int Cell::get_number_of_dynein_motors_IS()
{
    unsigned int number_of_motors = this->Capture_Shrinkage_dynein.get_dynein_motors_number();

    unsigned int sum_of_motor_in_micro = 0;
    unsigned int number_of_20_micro = 0;
    for( unsigned int motor_id = 0 ; motor_id < this->number_of_microtubules; motor_id ++ )
    {
        if( this->array_Of_Microtubules[ motor_id ].get_dynein_index() == 20 )
        {
            unsigned int number_of_dynein_in_micro = this->array_Of_Microtubules[ motor_id ].get_Dynein_abscissa().size();
            sum_of_motor_in_micro = sum_of_motor_in_micro + number_of_dynein_in_micro;
	    number_of_20_micro = number_of_20_micro + 1;
        }
    }

    cout<<"IS: sum_of_motor_in_micro = "<<sum_of_motor_in_micro<<endl;
    cout<<"IS: whole number = "<<number_of_motors + sum_of_motor_in_micro<<endl;
    cout<<"number of 20 micro = "<<number_of_20_micro<<endl;
    return number_of_motors + sum_of_motor_in_micro;
}


void Cell::microtubule_catch_pair_abscissa_real_dynein( unsigned int microtubule )
{

    std::uniform_real_distribution<> distribution{ 0 , 1 };
    double lenght_micro = this->array_Of_Microtubules[ microtubule ].get_lenght_of_microtubule();


    for( unsigned int segment_id = 1 ; segment_id < this->array_Of_Microtubules[ microtubule ].getNumberOfPoints()  - 1 ; segment_id ++ )
    {

        Vector3d lower_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id );
        Vector3d lower_bead_position_2 = this->project_point_on_surface( lower_bead_position );
        Vector3d upper_bead_position = this->array_Of_Microtubules[ microtubule ].getPoint( segment_id + 1 );
        Vector3d upper_bead_position_2 = this->project_point_on_surface( upper_bead_position );
        Vector3d tangent = this->array_Of_Microtubules[ microtubule ].getTangent2( segment_id );
        Vector3d tangent_2 = upper_bead_position_2 - lower_bead_position_2;
        Vector3d mini_segment = tangent_2 / ( double ) Dynein::number_of_segment_steps;

        std::vector< unsigned int > compartment_IDs;
        unsigned int last_index = 0;
        for( unsigned int counter = 0 ; counter < Dynein::number_of_segment_steps ; counter ++ )
        {
            Vector3d tmp_position = lower_bead_position_2 + ( double ) counter * mini_segment;
            unsigned int compartment_id = this->density_surface_dynein.get_dynein_compartment_id_projected( tmp_position );
            if( counter == 0 )
            {

                if( compartment_id != 0 )
                {
                    compartment_IDs.push_back( compartment_id );
                }
            }
            else
            {
                if( compartment_id != 0 )
                {
                    if( compartment_id != last_index )
                    {
                        compartment_IDs.push_back( compartment_id );
                    }
                }
            }
            last_index = compartment_id;
        }

        Vector3d axis_of_plane = lower_bead_position.cross( upper_bead_position );
        axis_of_plane = axis_of_plane / axis_of_plane.norm();
        for( unsigned int comp_index = 0 ; comp_index < compartment_IDs.size() ; comp_index ++ )
        {
            unsigned int compartment_id = compartment_IDs.at( comp_index );

            std::vector<Vector3d> dynein_points = this->get_dynein_in_compartment( compartment_id );
            std::vector< Vector3d > replace_motors;
            for( unsigned int bod_id  = 0 ; bod_id < dynein_points.size() ; bod_id ++ )
            {
                Vector3d motor_position = dynein_points.at( bod_id );

                Vector3d cl_p_of_seg( 0.0 , 0.0 , 0.0 );
                double distance = distance_point_segment( lower_bead_position , tangent , motor_position , cl_p_of_seg );

		    unsigned int number_of_generator = omp_get_thread_num();
		    double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
                    if( probability < Dynein_real::calculate_attachment_probability_per_time_step( distance ) )
                    {
                        double distance_lower_bead = 0;
                        for( unsigned int i = 0 ; i < segment_id ; i ++ )
                        {
                            distance_lower_bead = distance_lower_bead + this->array_Of_Microtubules[ microtubule ].getTangent2( i ).norm();
                        }

                        double distance_point_attachment = ( lower_bead_position - cl_p_of_seg ).norm();
                        double abscissa = ( distance_lower_bead + distance_point_attachment ) ;

                        if( abscissa >= lenght_micro - 8.0e-9 )
                        {
                            replace_motors.push_back( motor_position );
                            continue;
                        }
                        else
                        {

                            Vector3d new_motor_position = cl_p_of_seg;
                            std::pair < Vector3d , double > tmp_pair( new_motor_position , abscissa );
                            this->array_Of_Microtubules[ microtubule ].add_pair( tmp_pair );
                            this->array_Of_Microtubules[ microtubule ].set_dynein_index( 9 );
                        }

                    }
                    else
                    {
                        replace_motors.push_back( motor_position );
                    }

            }
            this->density_surface_dynein.set_dynein_points( compartment_id , replace_motors );
        }



    }

}

Vector3d Cell::get_IS_center(  )
{
	Vector3d front_center = this->get_IS_capture_shrinkage( ).get_center_of_IS_front();
	Vector3d center_on_sphere = front_center / front_center.norm() * Cell_parametres::A_AXIS;	
	return center_on_sphere;
}



double Cell::get_MTOC_IS_distance()
{
	Vector3d center_of_MTOC = this->get_MTOC_center(  );
	Vector3d IS_center = this->get_IS_center(  );
	Vector3d distance = center_of_MTOC - IS_center;
	double distance_1 = distance.norm();

	return distance_1;

}







void Cell::BIG_PRINT( unsigned int time_Step )
{

    this->print_time( time_Step );
    this->print_center_MTOC( );
}





void Cell::BIG_PRINT_numerical( unsigned int time_Step , unsigned int key )
{
    this->print_center_MTOC_numerical( time_Step , key );
    this->print_dynein_capture_shrinkage_numerical( time_Step , key );
    this->print_dynein_cortical_sliding_numerical( time_Step , key );

}


void Cell::BIG_PRINT_numerical_micro( unsigned int time_Step , unsigned int key )
{
    this->print_microtubules_numerical( time_Step , key );
    this->print_dynein_numerical( time_Step , key );
}






void Cell::print_center_numerical_parameters( unsigned int key )
{
      std::string path = "./picturesVideos/numerical_results/";
      std::string number = std::to_string(key);
      std::string name_of_text_file = path + "parameters_" + number + ".txt";

      Vector3d center = this->IS_Capture_Shrinkage.get_center_of_IS_front();
      ofstream fout;
      fout.open( name_of_text_file.c_str() , std::fstream::app );
      fout<<center( 0 )<< " " << center( 1 ) <<" "<<center( 2 )<<endl;  //write to file
      fout.close();


}



void Cell::print_time( unsigned int time_Step )
{
	FILE *out;
	out = fopen( "./picturesVideos/distance_results/time.txt" , "a");
    double time = ( double )( time_Step ) * sim_of_Cell::time_Step;
    fprintf(out,"%15.12f\n", time );
    fclose( out );
}

void Cell::print_center_MTOC()
{
	  FILE *out;
	  out = fopen( "./picturesVideos/distance_results/MTOC_center.txt" , "a");
	  Vector3d center = this->MTOC.get_center();
	  fprintf(out,"%15.12f %15.12f %15.12f\n", center( 0 ) , center( 1 ) , center( 2 ) );
	  fclose( out );
}



void Cell::print_center_MTOC_numerical( double time_Step , unsigned int key )
{
      double time = ( double )( time_Step ) * sim_of_Cell::time_Step;
      Vector3d center = this->MTOC.get_center();
      std::string path = "./picturesVideos/numerical_results/";
      std::string number = std::to_string(key);
      std::string name_of_text_file = path + "MTOC_center_" + number + ".txt";


      ofstream fout;
      fout.open( name_of_text_file.c_str() , std::fstream::app );
      fout<<time<<" "<<center( 0 )<< " " << center( 1 ) <<" "<<center( 2 )<<endl;  //write to file
      fout.close();
}

void Cell::print_dynein_capture_shrinkage_numerical( double time_Step , unsigned int key )
{
      double time = ( double )( time_Step ) * sim_of_Cell::time_Step;
      Vector3d center = this->MTOC.get_center();
      std::string path = "./picturesVideos/numerical_results/";
      std::string number = std::to_string(key);
      std::string name_of_text_file = path + "IS_CAP_SH_" + number + ".txt";
      unsigned int number_of_dynein = this->get_number_of_dynein_motors_in_micro_capture();


      ofstream fout;
      fout.open( name_of_text_file.c_str() , std::fstream::app );
      fout<<number_of_dynein<<endl;  //write to file
      fout.close();
}


void Cell::print_dynein_cortical_sliding_numerical( double time_Step , unsigned int key )
{
      double time = ( double )( time_Step ) * sim_of_Cell::time_Step;
      Vector3d center = this->MTOC.get_center();
      std::string path = "./picturesVideos/numerical_results/";
      std::string number = std::to_string(key);
      std::string name_of_text_file = path + "IS_CORTICAL_SLIDING_" + number + ".txt";
      unsigned int number_of_dynein = this->get_number_of_dynein_motors_cortical_sliding_micro();


      ofstream fout;
      fout.open( name_of_text_file.c_str() , std::fstream::app );
      fout<<number_of_dynein<<endl;  //write to file
      fout.close();
}












void Cell::print_dynein_numerical( double time_Step , unsigned int key )
{
      double time = ( double )( time_Step ) * sim_of_Cell::time_Step;
      std::string path = "./picturesVideos/numerical_results/";
      std::string number_string = std::to_string( key );
      std::string time_string = std::to_string( time_Step );
      std::string name_of_text_file = path + "dynein_" + time_string + "_"+ number_string + ".txt";



      ofstream fout;
      fout.open( name_of_text_file.c_str() , std::fstream::app );
      std::map< unsigned int , std::vector<Vector3d> > mapa = this->density_surface_dynein.get_map();
      for(  std::map< unsigned int , std::vector<Vector3d> >::iterator it=mapa.begin(); it != mapa.end(); ++it)
      {

          std::vector<Vector3d> vector_tmp = it->second;
          for( unsigned int i = 0 ; i < vector_tmp.size() ; i ++ )
          {
                Vector3d vect = vector_tmp.at( i );
      		fout<<vect( 0 )<<" "<< vect( 1 ) <<" "<<vect( 2 )<<endl;  //write to file
          }
      }
      fout<<"END"<<endl;

      for( unsigned int micro = 0 ; micro < this->number_of_microtubules ; micro ++ )
      {
          if( this->array_Of_Microtubules[ micro ].get_dynein_index() == 9 )
          {
                std::vector< std::pair < Vector3d ,double  > > vector_dynein_abscissa = this->array_Of_Microtubules[ micro ].get_Dynein_abscissa();
                if( vector_dynein_abscissa.size() < 1 )
                {
                    continue;
                }
                else
                {
                      for( unsigned int point = 0 ; point < vector_dynein_abscissa.size() ; point ++ )
                      {
                          std::pair < Vector3d ,double  > pair_tmp = vector_dynein_abscissa.at( point );
                          Vector3d position = std::get< 0 >( pair_tmp );
                          fout<<position( 0 )<<" "<< position( 1 ) <<" "<<position( 2 )<<endl;
                      }
                }
          }
      }
      fout<<"END"<<endl;
      for( unsigned int micro = 0 ; micro < this->number_of_microtubules ; micro ++ )
      {
          if( this->array_Of_Microtubules[ micro ].get_dynein_index() == 20 )
          {
                std::vector< std::pair < Vector3d ,double  > > vector_dynein_abscissa = this->array_Of_Microtubules[ micro ].get_Dynein_abscissa();
                if( vector_dynein_abscissa.size() < 1 )
                {
                    continue;
                }
                else
                {
                      for( unsigned int point = 0 ; point < vector_dynein_abscissa.size() ; point ++ )
                      {
                          std::pair < Vector3d ,double  > pair_tmp = vector_dynein_abscissa.at( point );
                          Vector3d position = std::get< 0 >( pair_tmp );
                          fout<<position( 0 )<<" "<< position( 1 ) <<" "<<position( 2 )<<endl;
                      }
                }
          }
      }

       fout.close();

}





void Cell::print_dynein_numerical_unattached( double time_Step , unsigned int key )
{
      double time = ( double )( time_Step ) * sim_of_Cell::time_Step;
      std::string path = "./picturesVideos/numerical_results/";
      std::string number_string = std::to_string( key );
      unsigned int time_Step_2 = time_Step;

      std::string time_string = std::to_string( time_Step_2 );
      std::string name_of_text_file = path + "dynein_unattached_" + time_string + "_"+ number_string + ".txt";



      ofstream fout;
      fout.open( name_of_text_file.c_str() , std::fstream::app );
      std::map< unsigned int , std::vector<Vector3d> > mapa = this->density_surface_dynein.get_map();
      for(  std::map< unsigned int , std::vector<Vector3d> >::iterator it=mapa.begin(); it != mapa.end(); ++it)
      {

          std::vector<Vector3d> vector_tmp = it->second;
          for( unsigned int i = 0 ; i < vector_tmp.size() ; i ++ )
          {
                Vector3d vect = vector_tmp.at( i );
      		fout<<vect( 0 )<<" "<< vect( 1 ) <<" "<<vect( 2 )<<endl;  //write to file
          }
      }
      fout<<"END"<<endl;

			mapa = this->Capture_Shrinkage_dynein.get_map();
      for(  std::map< unsigned int , std::vector<Vector3d> >::iterator it=mapa.begin(); it != mapa.end(); ++it)
      {

          std::vector<Vector3d> vector_tmp = it->second;
          for( unsigned int i = 0 ; i < vector_tmp.size() ; i ++ )
          {
                Vector3d vect = vector_tmp.at( i );
      		fout<<vect( 0 )<<" "<< vect( 1 ) <<" "<<vect( 2 )<<endl;  //write to file
          }
      }


       fout.close();

}













void Cell::print_microtubules_numerical( unsigned int time_Step , unsigned int key )
{
      std::string path = "./picturesVideos/numerical_results/";
      std::string number_string = std::to_string( key );
      std::string time_string = std::to_string( time_Step );
      std::string name_of_text_file = path + "micro_" + time_string + "_"+ number_string + ".txt";


      ofstream fout;
      fout.open( name_of_text_file.c_str() , std::fstream::app );

      for( unsigned int microtubule_index = 0 ; microtubule_index < this->number_of_microtubules ; microtubule_index ++ )
      {
          for( unsigned int bead_index = 0 ; bead_index < this->array_Of_Microtubules[ microtubule_index ].getNumberOfPoints() ; bead_index ++ )
          {
                Vector3d micro_point = this->array_Of_Microtubules[ microtubule_index ].getPoint( bead_index );
                fout<<micro_point( 0 )<<" "<< micro_point( 1 ) <<" "<<micro_point( 2 )<<endl;
          }
          fout<<this->array_Of_Microtubules[ microtubule_index ].get_dynein_index()<<endl;
	  fout<<"END"<<endl;
      }



      fout.close();
}






void Cell::print_center_MTOC_radius_numerical( unsigned int key )
{
      Vector3d center = this->MTOC.get_center();
      double radius = center.norm();


      std::string path = "./picturesVideos/numerical_results/";
      std::string number = std::to_string(key);
      std::string name_of_text_file = path + "MTOC_radius_" + number + ".txt";


      ofstream fout;
      fout.open( name_of_text_file.c_str() , std::fstream::app );
      fout<<radius<<endl;  //write to file
      fout.close();

}


void Cell::print_center_MTOC_azimuth_numerical( unsigned int key )
{
      Vector3d referential_axis( 1.0 , 0.0 , 0.0);
      Vector3d referential_axis_2( 0.0 , 1.0 , 0.0);    //( 1.0 , 0.0 , 0.0)
      Vector3d center = this->MTOC.get_center();
      center( 2 ) = 0;
      center = center / center.norm();
      double cosinus_angle = center.dot( referential_axis );
      double angle = acos( cosinus_angle );

      double cos_tmp = center.dot( referential_axis_2 );
      if( cos_tmp > 0 )
      {

      }
      else
      {
            angle = 2.0 * sim_of_Cell::PI - angle;
      }


      std::string path = "./picturesVideos/numerical_results/";
      std::string number = std::to_string(key);
      std::string name_of_text_file = path + "MTOC_azimuth_" + number + ".txt";

      ofstream fout;
      fout.open( name_of_text_file.c_str() , std::fstream::app );
      fout<<angle<<endl;  //write to file
      fout.close();




}





void Cell::force_whole_microtubule_cell_wall( MatrixXd &force , unsigned int index )
{
	if( index >= this->number_of_microtubules )
	{
		cout<<"index >= this->number_of_microtubules in Cell::force_whole_microtubule( MatrixXd &force , unsigned int index )"<<endl;
		cout<<"index = "<<index<<endl;
		throw("");
	}
	Microtubule temporary_microtubule = this->getMicrotubule( index );

	if( force.rows() != 3 * temporary_microtubule.getNumberOfPoints()  )
	{
		cout<<"force.rows() != 3 * temporary_microtubule.getNumberOfPoints() in force_whole_microtubule"<<endl;
		cout<<"index = "<<index<<endl;
		throw("");
	}


	for( unsigned int bead_index = 0 ; bead_index < temporary_microtubule.getNumberOfPoints() ; bead_index ++ )
	{
		Vector3d position = temporary_microtubule.getPoint( bead_index );//force_position_cell_wall_two_elipsoid
        	Vector3d force_one_bead = this->force_position_cell_wall_two_elipsoid( position );//force_position_cell_wall_simple_elipsoid
       		for( unsigned int index_dimension = 0 ; index_dimension < 3 ; index_dimension ++ )
		{
			force( 3 * bead_index + index_dimension , 0 ) = force_one_bead( index_dimension );
		}
	}
	if( this->array_Of_Microtubules[ index ].getNumberOfPoints() == 20 )
    	{
        	for( unsigned int index_dimension = 0 ; index_dimension < 3 ; index_dimension ++ )
		{
			force( 3 * ( this->array_Of_Microtubules[ index ].getNumberOfPoints() - 1 ) + index_dimension , 0 ) = 0;
		}
    	}


}







void Cell::timeDevelopment_numerical_results( double Time_0 , double Time_1 , unsigned int key )
{

    	double numberOfSteps =  ( ( Time_1 - Time_0 ) / sim_of_Cell::time_Step );
	unsigned int numberOfSteps2 = nearbyint ( numberOfSteps ); //This is control for unprecise values of doubles in C++

	unsigned int print_integer = 0;

	this->print_center_numerical_parameters( key );
	unsigned int printing_counter = 0;

        this->BIG_PRINT_numerical_micro( printing_counter , key );
	printing_counter = printing_counter + 1;
	for( unsigned int step = 0 ; step < numberOfSteps2 ; step ++ )
	{

        	double time = step * sim_of_Cell::time_Step;
        	if( step % 10000 == 0)
		{

            		cout<<"key = "<<key<<endl;
            		cout<<"time = "<<time<<endl;
        	}


      		if( step % 50 == 0 )
		{
			this->BIG_PRINT_numerical( step , key );

        	}

      		if( step % 50 == 0 )
		{
			//This is the distance control
                	double MTOC_IS_distance = this->get_MTOC_IS_distance();
                	if( MTOC_IS_distance < sim_of_Cell::treshold_distance )
                	{
                		break;
                	}
              	}    


      		if( ( step %  7142 ) == 0 ) //14285
		{

			if( sim_of_Cell::density_of_dynein_motor_surface < 0.001  )
			{

				ISCorticalSl tmp_1 = this->get_IS_cortical_sliding_first();
				unsigned int number_caught_dynein = get_number_of_dynein_motors_cortical_sliding_micro();
				Surface surface( "Cortical_Sliding" , tmp_1 , number_caught_dynein );
				this->set_density_surface_dynein( surface );

			}
			else
			{
                                string cortical_Sliding = "cortical_Sliding";
				ISCorticalSl tmp_1 = this->get_IS_cortical_sliding_first();
				double concentration_in_IS = 0;
				double concentration_outside_IS = 0;
				this->get_attached_dynein_motors_cortical_sliding_micro_both_surface_IS( concentration_in_IS ,  concentration_outside_IS );
				Surface first_Surface(  concentration_outside_IS , sim_of_Cell::surface_width , cortical_Sliding , concentration_in_IS , tmp_1 );  //get_dynein_motors_number()
				this->set_density_surface_dynein( first_Surface );
			}


			print_integer = step + 1000;
        	}



      		if( step == print_integer )
		{
			this->BIG_PRINT_numerical_micro( printing_counter , key );
			printing_counter = printing_counter + 1;
        	}

		//FUNCTIONS
        	//////////////////////////////////////
        	this->MidStep_3();//
        	//////////////////////////////////////
        	this->stepping_and_detachment_of_all_microtubule_projection_real_dynein();

        	this->check_caught_micro_IS_with_real_dynein_all_micro( );
        	this->microtubule_catch_pair_abscissa_real_dynein_in_IS_all_micro( );

		this->catch_pair_abscissa_real_dynein();
		this->control_length_of_micro_IS_2();

	}

}






void Cell::timeDevelopment_numerical_results_2( double Time_0 , double Time_1 , unsigned int key )
{

    	double numberOfSteps =  ( ( Time_1 - Time_0 ) / sim_of_Cell::time_Step );
	unsigned int numberOfSteps2 = nearbyint ( numberOfSteps ); //This is control for unprecise values of doubles in C++

	unsigned int print_integer = 0;

	this->print_center_numerical_parameters( key );
	unsigned int printing_counter = 0;

	for( unsigned int step = 0 ; step < numberOfSteps2 ; step ++ )
	{

        	double time = step * sim_of_Cell::time_Step;
        	if( step % 1000 == 0)
		{
            		cout<<"key = "<<key<<endl;
            		cout<<"time = "<<time<<endl;
        	}


      		if( step % 50 == 0 )
		{
			this->BIG_PRINT_numerical( step , key );
        	}

      		if( step % 50 == 0 )
		{
                	double MTOC_IS_distance = this->get_MTOC_IS_distance();
                	if( MTOC_IS_distance < sim_of_Cell::treshold_distance )
                	{
                		break;
                	}
              	} 



      		if( ( step %  200 ) == 0 )
		{


			if ( sim_of_Cell::density_of_dynein_motor_surface < 0.00001 )
			{	//I fear double imprecisions
				
				this->density_surface_dynein.change_part_of_IS_points( this->IS_Cortic_Sl );
			}
			else
			{
				cout<<"sim_of_Cell::density_of_dynein_motor_surface > 0.00001"<<endl;
				cout<<"void Cell::timeDevelopment_numerical_results_2( double Time_0 , double Time_1 , unsigned int key )"<<endl;
				cout<<"ERROR_ID = 43616121561512"<<endl;
				throw("");
				this->density_surface_dynein.change_part_of_points_surface_and_IS( this->IS_Cortic_Sl );
			}
        	}

		if( sim_of_Cell::big_print_switch == true )
		{
	      		if( ( step %  7142 ) == 0 )
			{
				this->BIG_PRINT_numerical_micro( printing_counter , key );
				printing_counter = printing_counter + 1;
			}
		}


		//FUNCTIONS
        	//////////////////////////////////////
        	this->MidStep_3(); //

        	//////////////////////////////////////

        	this->stepping_and_detachment_of_all_microtubule_projection_real_dynein();
        	this->check_caught_micro_IS_with_real_dynein_all_micro( );
        	this->microtubule_catch_pair_abscissa_real_dynein_in_IS_all_micro( );
		this->catch_pair_abscissa_real_dynein();
      		if( ( step %  100 ) == 0 ) 
		{
			this->test_dynein( "000000000000000" );
			this->density_surface_dynein.control_points_outside_IS(  this->IS_Cortic_Sl );
		}
		this->control_length_of_micro_IS_2();
	}

}














Cell::~Cell()
{

	delete[] array_Of_Microtubules;
	array_Of_Microtubules = NULL;
	// TODO Auto-generated destructor stub
}

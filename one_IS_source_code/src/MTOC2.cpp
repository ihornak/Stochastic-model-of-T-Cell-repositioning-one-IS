//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS
#include "MTOC2.h"

MTOC2::MTOC2()
{
    this->number_of_points = 1;
    this->coordinates = MatrixXd::Zero( 3 * ( this->number_of_points + 1 ) , 1 );
    this->effectiveFriction = 1e-6;
    this->radius_of_MTOC = this->set_radius();
    this->polygon_angle = this->set_polygon_angle();
    this->original_orientations = MatrixXd::Zero( 3 * ( this->number_of_points + 1 ) , 1 );
    this->original_coordinates = this->coordinates;
    this->time_clock = 0;
}

MTOC2::MTOC2( unsigned int number_of_points_arg )
{

	if( ( number_of_points_arg % 4 != 0 ) && ( number_of_points_arg > 4 ) )
	{
		cout<<"number_of_points_arg % 4 != 0"<<endl;
        	cout<<"( number_of_points_arg % 4 != 0 ) && ( number_of_points_arg > 4 )"<<endl;
		cout<<"in MTOC2::MTOC2( unsigned int number_of_points_arg)"<<endl;
		throw("");
	}
	if( number_of_points_arg > MTOCparam::boundary_MTOC_points )
    	{
		cout<<"number_of_points_arg = "<<number_of_points_arg<<endl;
        	cout<<"ERROR_ID MTOC2468418946"<<endl;
		cout<<"in MTOC2::MTOC2( unsigned int number_of_points_arg)"<<endl;
		throw("");
    	}


	this->number_of_points = number_of_points_arg;
	this->coordinates = MatrixXd::Zero( 3 * ( this->number_of_points + 1 ) , 1 );
    	this->set_friction();
    	this->radius_of_MTOC = this->set_radius();
    	this->polygon_angle = this->set_polygon_angle();

    if( number_of_points_arg == 0 )
    {
        cout<<" number_of_points_arg == 0 "<<endl;
        cout<<"in MTOC2::MTOC2( unsigned int number_of_points_arg )"<<endl;
        throw("");
    }
	if( number_of_points_arg == 1 )
    {
        this->coordinates( 3 , 0 ) = this->radius_of_MTOC;
    }
	else if( number_of_points_arg == 2 )
    {
        this->coordinates( 3 , 0 ) = this->radius_of_MTOC;
        this->coordinates( 7 , 0 ) = this->radius_of_MTOC;
    }
    else if( number_of_points_arg == 3 )
    {
        this->coordinates( 3 , 0 ) = this->radius_of_MTOC;
        double length = this->radius_of_MTOC * cos( sim_of_Cell::PI / 4.0 );

        this->coordinates( 6 , 0 ) = - length;
        this->coordinates( 7 , 0 ) = length;

        this->coordinates( 9 , 0 ) = - length;
        this->coordinates( 10 , 0 ) = - length;
    }
	else if( number_of_points_arg == 4 )
	{
		this->coordinates( 3 , 0 ) = this->radius_of_MTOC;
		this->coordinates( 7 , 0 ) = this->radius_of_MTOC;
		this->coordinates( 9 , 0 ) = ( -1.0 ) * this->radius_of_MTOC;
		this->coordinates( 13 , 0 ) =  ( -1.0 ) * this->radius_of_MTOC;
	}
	else
    {
        Vector3d array_of_axis[ 4 ];
        array_of_axis[ 0 ] = Vector3d( 1.0 , 0.0 , 0.0 );
        array_of_axis[ 1 ] = Vector3d( 0.0 , 1.0 , 0.0 );
        array_of_axis[ 2 ] = Vector3d( - 1.0 , 0.0 , 0.0 );
        array_of_axis[ 3 ] = Vector3d( 0.0 , -1.0 , 0.0 );


        double by_the_ax = this->radius_of_MTOC * cos( this->polygon_angle );
        double vertical_to_the_ax = this->radius_of_MTOC * sin( this->polygon_angle );
        unsigned int number_of_points_in_polygon = number_of_points_arg / 4;
        double angle_in_polygon = 2.0 * sim_of_Cell::PI / ( double )number_of_points_in_polygon;

        for( unsigned int polygon = 0 ; polygon < 4 ; polygon ++ )
        {
            Vector3d axis = array_of_axis[ polygon ];
            Vector3d first_point = axis * by_the_ax;
            first_point( 2 ) = vertical_to_the_ax;
            for( unsigned int dim_j = 0 ; dim_j < 3 ; dim_j ++ )
            {
                this->coordinates( 3 * ( polygon * number_of_points_in_polygon + 1 ) + dim_j , 0 ) = first_point( dim_j );
            }
            for( unsigned int point = 1 ; point < number_of_points_in_polygon ; point ++ )
            {
                Quaternion<double> quaternion;
                quaternion = AngleAxis<double>( angle_in_polygon , axis );
                first_point = quaternion * first_point;
                for( unsigned int dim_j = 0 ; dim_j < 3 ; dim_j ++ )
                {
                    this->coordinates( 3 * ( polygon * number_of_points_in_polygon + point + 1 ) + dim_j , 0 ) = first_point( dim_j );
                }
            }
        }
    }
    this->original_coordinates = this->coordinates;
    this->time_clock = 0;
}

MTOC2::MTOC2( unsigned int number_of_points_arg , Vector3d center_of_MTOC  )
{
	if( number_of_points_arg % 4 != 0 )
	{
		cout<<"number_of_points_arg % 4 != 0"<<endl;
		cout<<"in MTOC2::MTOC2( unsigned int number_of_points_arg)"<<endl;
		throw("");
	}
	if( number_of_points_arg > MTOCparam::boundary_MTOC_points )
    {
		cout<<"number_of_points_arg = "<<number_of_points_arg<<endl;
        	cout<<"ERROR_ID MTOC2468418946"<<endl;
		cout<<"in MTOC2::MTOC2( unsigned int number_of_points_arg)"<<endl;
		throw("");
    }


	this->number_of_points = number_of_points_arg;
	this->coordinates = MatrixXd::Zero( 3 * ( this->number_of_points + 1 ) , 1 );
    	this->set_friction();
    	this->radius_of_MTOC = this->set_radius();
    	this->polygon_angle = this->set_polygon_angle();

    if( number_of_points_arg == 0 )
    {
        cout<<" number_of_points_arg == 0 "<<endl;
        cout<<"in MTOC2::MTOC2( unsigned int number_of_points_arg )"<<endl;
        throw("");
    }
	if( number_of_points_arg == 1 )
    {
        this->coordinates( 3 , 0 ) = this->radius_of_MTOC;
    }
	else if( number_of_points_arg == 2 )
    {
        //if I do it in x axis, the sparse matrix_2 will not be invertible
        this->coordinates( 3 , 0 ) = this->radius_of_MTOC;
        this->coordinates( 7 , 0 ) = this->radius_of_MTOC;
    }
    else if( number_of_points_arg == 3 )
    {
        //if I do it in x axis, the sparse matrix_2 will not be invertible
        this->coordinates( 3 , 0 ) = this->radius_of_MTOC;
        double length = this->radius_of_MTOC * cos( sim_of_Cell::PI / 4.0 );

        this->coordinates( 6 , 0 ) = - length;
        this->coordinates( 7 , 0 ) = length;

        this->coordinates( 9 , 0 ) = - length;
        this->coordinates( 10 , 0 ) = - length;
    }
	else if( number_of_points_arg == 4 )
	{
		this->coordinates( 3 , 0 ) = this->radius_of_MTOC;
		this->coordinates( 7 , 0 ) = this->radius_of_MTOC;
		this->coordinates( 9 , 0 ) = ( -1.0 ) * this->radius_of_MTOC;
		this->coordinates( 13 , 0 ) =  ( -1.0 ) * this->radius_of_MTOC;
	}
	else
    {
        Vector3d array_of_axis[ 4 ];
        array_of_axis[ 0 ] = Vector3d( 1.0 , 0.0 , 0.0 );
        array_of_axis[ 1 ] = Vector3d( 0.0 , 1.0 , 0.0 );
        array_of_axis[ 2 ] = Vector3d( - 1.0 , 0.0 , 0.0 );
        array_of_axis[ 3 ] = Vector3d( 0.0 , -1.0 , 0.0 );

        double by_the_ax = this->radius_of_MTOC * cos( this->polygon_angle );
        double vertical_to_the_ax = this->radius_of_MTOC * sin( this->polygon_angle );
        unsigned int number_of_points_in_polygon = number_of_points_arg / 4;
        double angle_in_polygon = 2.0 * sim_of_Cell::PI / ( double )number_of_points_in_polygon;

        for( unsigned int polygon = 0 ; polygon < 4 ; polygon ++ )
        {
            Vector3d axis = array_of_axis[ polygon ];
            Vector3d first_point = axis * by_the_ax;
            first_point( 2 ) = vertical_to_the_ax;
            for( unsigned int dim_j = 0 ; dim_j < 3 ; dim_j ++ )
            {
                this->coordinates( 3 * ( polygon * number_of_points_in_polygon + 1 ) + dim_j , 0 ) = first_point( dim_j );
            }
            for( unsigned int point = 1 ; point < number_of_points_in_polygon ; point ++ )
            {
                Quaternion<double> quaternion;
                quaternion = AngleAxis<double>( angle_in_polygon , axis );
                first_point = quaternion * first_point;
                for( unsigned int dim_j = 0 ; dim_j < 3 ; dim_j ++ )
                {
                    this->coordinates( 3 * ( polygon * number_of_points_in_polygon + point + 1 ) + dim_j , 0 ) = first_point( dim_j );
                }
            }
        }
    }

    for( unsigned int i = 0 ; i <= this->number_of_points ; i ++ )
    {
        for( unsigned int j = 0 ; j  < 3 ; j ++  )
        {
            this->coordinates( 3 * i + j , 0 ) = this->coordinates( 3 * i + j , 0 ) + center_of_MTOC( j );
        }
    }

    this->original_coordinates = this->coordinates;
    this->time_clock = 0;
}


MTOC2::MTOC2( unsigned int number_of_points_arg , Vector3d center_of_MTOC , string plane )
{

    if( number_of_points_arg > MTOCparam::boundary_MTOC_points )
    {
		cout<<"number_of_points_arg = "<<number_of_points_arg<<endl;
        	cout<<"ERROR_ID = 97442468418946"<<endl;
		cout<<"in MTOC2::MTOC2( unsigned int number_of_points_arg)"<<endl;
		throw("");
    }

    this->number_of_points = number_of_points_arg;
    this->coordinates = MatrixXd::Zero( 3 * ( this->number_of_points + 1 ) , 1 );
    this->set_friction();
    this->radius_of_MTOC = this->set_radius();
    this->polygon_angle = this->set_polygon_angle();
    this->original_MTOC_orientation = center_of_MTOC / center_of_MTOC.norm();
    this->original_orientations = MatrixXd::Zero( 3 * ( this->number_of_points ) , 1 );


    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        this->coordinates( dimension , 0 ) = center_of_MTOC( dimension );
    }
    //the axis is always the center of the polygon
    double delta_angle = 2.0 * sim_of_Cell::PI / ( double )number_of_points_arg;
    double initial_angle = - delta_angle * 2.0; //!!MAGIC NUMBER

    //construction on x-y plane + z of mtoc center
    for( unsigned int point_index = 0 ; point_index <  number_of_points_arg ; point_index ++ )
    {
        double angle = initial_angle + point_index * delta_angle;
        double x = cos( angle ) * this->radius_of_MTOC;
        double y = sin( angle ) * this->radius_of_MTOC;
        this->coordinates( 3 * ( point_index + 1 ) , 0 ) = x;
        this->coordinates( 3 * ( point_index + 1 ) + 1 , 0 ) = y;
        this->coordinates( 3 * ( point_index + 1 ) + 2 , 0 ) = this->get_center()( 2 );

    }


    Vector3d center = this->get_center();
    for( unsigned int point_index = 1 ; point_index <= this->get_number_of_points() ; point_index ++ )
    {
        Vector3d point = this->get_point( point_index );
        Vector3d orientation_tmp = point - center;
        orientation_tmp = orientation_tmp / orientation_tmp.norm();

        for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
        {
            this->original_orientations( 3 * ( point_index - 1 ) + dimension , 0 ) = orientation_tmp( dimension );
        }
    }

    this->original_coordinates = this->coordinates;
    this->time_clock = 0;


}





MTOC2::MTOC2( unsigned int number_of_points_arg , Vector3d center_of_MTOC , string plane  , string plane_2 )
{
    if( number_of_points_arg > MTOCparam::boundary_MTOC_points )
    {
		cout<<"number_of_points_arg = "<<number_of_points_arg<<endl;
        	cout<<"ERROR_ID = 97442468418946"<<endl;
		cout<<"in MTOC2::MTOC2( unsigned int number_of_points_arg)"<<endl;
		throw("");
    }

    this->number_of_points = number_of_points_arg;
	this->coordinates = MatrixXd::Zero( 3 * ( this->number_of_points + 1 ) , 1 );
    this->set_friction();
    this->radius_of_MTOC = this->set_radius();
    this->polygon_angle = this->set_polygon_angle();
    this->original_MTOC_orientation = center_of_MTOC / center_of_MTOC.norm();
    this->original_orientations = MatrixXd::Zero( 3 * ( this->number_of_points ) , 1 );




    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        this->coordinates( dimension , 0 ) = center_of_MTOC( dimension );
    }
    double delta_angle = 2.0 * sim_of_Cell::PI / ( double )number_of_points_arg;
    double initial_angle = - delta_angle * 2.0; //!!MAGIC NUMBER

    std::uniform_real_distribution<> distribution_angles{ 0 , delta_angle };
    unsigned int number_of_generator = omp_get_thread_num();
    double azimutal_angle_tmp = distribution_angles( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
    initial_angle = initial_angle + azimutal_angle_tmp;

    for( unsigned int point_index = 0 ; point_index <  number_of_points_arg ; point_index ++ )
    {
        double angle = initial_angle + point_index * delta_angle;
        double x = cos( angle ) * this->radius_of_MTOC;
        double y = sin( angle ) * this->radius_of_MTOC;
        this->coordinates( 3 * ( point_index + 1 ) , 0 ) = x;
        this->coordinates( 3 * ( point_index + 1 ) + 1 , 0 ) = y;
        this->coordinates( 3 * ( point_index + 1 ) + 2 , 0 ) = this->get_center()( 2 );

    }

    Vector3d center = this->get_center();
    for( unsigned int point_index = 1 ; point_index <= this->get_number_of_points() ; point_index ++ )
    {
        Vector3d point = this->get_point( point_index );
        Vector3d orientation_tmp = point - center;
        orientation_tmp = orientation_tmp / orientation_tmp.norm();

        for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
        {
            this->original_orientations( 3 * ( point_index - 1 ) + dimension , 0 ) = orientation_tmp( dimension );
        }
    }
}


unsigned int MTOC2::get_number_of_points() const
{
    return this->number_of_points;
}



void MTOC2::set_friction()
{
    double d_v = 0.5e-6;
    double K = 1.0;
    double effectiveF = sim_of_Cell::multiply_friction_constant * 3.0 * sim_of_Cell::PI * sim_of_Cell::viscosity * d_v * K;
    this->effectiveFriction = 3.58369e-07;
}



void MTOC2::add_vector( Vector3d position )
{
    for( unsigned int i = 0 ; i <= this->number_of_points ; i ++ )
    {
        for( unsigned int j = 0 ; j < 3 ; j ++ )
        {
            this->coordinates( 3 * i + j , 0 ) = this->coordinates( 3 * i + j , 0 ) + position( j );
        }
    }
}

MTOC2& MTOC2::operator=( const MTOC2& tmp )
{
    if( this == &tmp )
	{
		cout<<" this == &tmp  in MTOC2::operator=( const MTOC2 &tmp )"<<endl;
		throw("");
	}
    this->number_of_points = tmp.get_number_of_points();
    this->coordinates = tmp.get_coordinates();
    this->radius_of_MTOC = tmp.radius_of_MTOC;
    this->polygon_angle = tmp.polygon_angle;
    this->effectiveFriction = tmp.effectiveFriction;
    this->original_orientations = tmp.original_orientations;
    this->original_MTOC_orientation = tmp.original_MTOC_orientation;
    this->original_coordinates = tmp.original_coordinates;
    this->time_clock = tmp.time_clock;
    return *this;
}

double MTOC2::set_radius(  )
{
    return MTOCparam::MTOC_radius;
}

double MTOC2::set_polygon_angle(  )
{
    return 3.1415926535897 / 2.0 / 8.0;
}




MatrixXd MTOC2::get_coordinates() const
{
    return this->coordinates;
}


void MTOC2::set_coordinates( MatrixXd& coordinates_arg )
{
     this->coordinates = coordinates_arg;
}


MatrixXd MTOC2::get_original_coordinates() const
{
    return this->original_coordinates;
}


Vector3d MTOC2::get_point( unsigned int number ) const
{
	if( ( number == 0 ) || ( number > this->number_of_points ) )
	{
                cout<<"( number == 0 ) || ( number > this->number_of_points )"<<endl;
		cout<<" number = "<<number<<endl;
		cout<<" in Vector3d MTOC2::get_point( unsigned int number ) "<<endl;
                cout<<"ERROR_ID_MTOC2_2469464351"<<endl;
		throw("");
		return Vector3d( 0.0 , 0.0 , 0.0 );
	}
	else
	{
		Vector3d point_position( 0.0 , 0.0 , 0.0 );

		for( unsigned int i = 0 ; i < 3 ; i ++ )
		{
			point_position( i ) = this->coordinates( 3 * number + i , 0);
		}
		return point_position;
	}

}



Vector3d MTOC2::get_point_original( unsigned int number ) const
{
	if( ( number == 0 ) || ( number > this->number_of_points ) )
	{
                cout<<"( number == 0 ) || ( number > this->number_of_points )"<<endl;
		cout<<" number = "<<number<<endl;
		cout<<" in Vector3d MTOC2::get_point_original( unsigned int number ) "<<endl;
                cout<<"ERROR_ID_MTOC2_221561311"<<endl;
		throw("");
		return Vector3d( 0.0 , 0.0 , 0.0 );
	}
	else
	{
		Vector3d point_position( 0.0 , 0.0 , 0.0 );

		for( unsigned int i = 0 ; i < 3 ; i ++ )
		{
			point_position( i ) = this->original_coordinates( 3 * number + i , 0);
		}
		return point_position;
	}

}


Vector3d MTOC2::get_original_orientation( unsigned int number ) const
{
    if( ( number == 0 ) || ( number > this->number_of_points ) )
	{
        	cout<<"( number == 0 ) || ( number > this->number_of_points )"<<endl;
		cout<<" number = "<<number<<endl;
		cout<<" in Vector3d MTOC2::get_original_orientation( unsigned int number ) "<<endl;
        	cout<<"ERROR_ID = 58168321648687645"<<endl;
		throw("");
		return Vector3d( 0.0 , 0.0 , 0.0 );
	}
	else
	{
		Vector3d point_position( 0.0 , 0.0 , 0.0 );

		for( unsigned int i = 0 ; i < 3 ; i ++ )
		{
			point_position( i ) = this->original_orientations( 3 * ( number - 1 ) + i , 0 );
		}
		return point_position;
	}

}






Vector3d MTOC2::get_center( ) const
{
	Vector3d point_position( 0.0 , 0.0 , 0.0 );
	for( unsigned int i = 0 ; i < 3 ; i ++ )
	{
		point_position( i ) = this->coordinates( i , 0);
	}
	return point_position;

}


Vector3d MTOC2::get_center_of_gravity() const
{
	Vector3d center = this->get_center( );
        for( unsigned int point_id = 1 ; point_id <= this->get_number_of_points() ; point_id ++ )
        {
        	Vector3d point = this->get_point( point_id );
		center = center + point;
        }

	center = center / ( ( double ) ( get_number_of_points() + 1 )  );
	return center;

}



Vector3d MTOC2::get_center_original( ) const
{
	Vector3d point_position( 0.0 , 0.0 , 0.0 );
	for( unsigned int i = 0 ; i < 3 ; i ++ )
	{
		point_position( i ) = this->original_coordinates( i , 0);
	}
	return point_position;

}



double MTOC2::get_radius() const
{
    return this->radius_of_MTOC;
}



Vector3d MTOC2::tangent_to_center( unsigned int number )
{
	Vector3d tmp = ( this->get_point( number ) - this->get_center() );
	return tmp;
}


Vector3d MTOC2::tangent_from_next_bead( unsigned int number )
{
    if( ( number > this->number_of_points )  || (  number == 0 )  )
    {
        cout<<"( number >= this->number_of_points )  || (  number == 0 ) "<<endl;
        cout<<"number = "<<number<<endl;
        cout<<"Vector3d MTOC2::tangent_to_next_bead( unsigned int number )"<<endl;
        throw("");
    }

    Vector3d tangent( 1.0 ,  1.0 ,  1.0  );
    if( number < this->number_of_points )
    {
    	tangent = this->get_point( number ) - this->get_point( number + 1 );
    }
    else if ( number == this->number_of_points )
    {
	tangent = this->get_point( 1 ) - this->get_point( number );
    }


    return tangent;
}



Vector3d MTOC2::tangent_between_beads(  unsigned int number_1 , unsigned int number_2 )
{
    if( number_1 >= this->number_of_points )
    {
        cout<<"number_1 >= this->number_of_points"<<endl;
        cout<<"number_1 = "<<number_1<<endl;
        cout<<"Vector3d MTOC2::tangent_between_beads( unsigned int number_1 )"<<endl;
        throw("");
    }

    if( number_2 > this->number_of_points )
    {
        cout<<"number_2 >= this->number_of_points"<<endl;
        cout<<"number_2 = "<<number_2<<endl;
        cout<<"Vector3d MTOC2::tangent_between_beads( unsigned int number_2 )"<<endl;
        throw("");
    }
    if( number_2 == number_1 )
    {
        cout<<"number_2 == number_1"<<endl;
        cout<<"Vector3d MTOC2::tangent_between_beads( unsigned int number_2 )"<<endl;
        throw("");
    }


    Vector3d tangent = this->get_point( number_1 ) - this->get_point( number_2 );
    return tangent;
}




void MTOC2::controlMTOC( )
{
    for( unsigned int point = 1 ; point <= this->number_of_points ; point++ )
    {
        double distance = this->tangent_to_center( point ).norm();
        double percentage = abs( ( distance - this->radius_of_MTOC ) / this->radius_of_MTOC );


        if( percentage > 0.05 )
        {
            this->adjustMTOC();
            break;
        }

    }

}




void MTOC2::adjustMTOC( )
{
    for( unsigned int i = 0 ; i < this->number_of_points ; i ++ )
    {
        Vector3d tangent = this->tangent_to_center( i + 1 );
        Vector3d orientation = tangent / tangent.norm();
        tangent = orientation * this->radius_of_MTOC;
        for( unsigned int j = 0 ; j < 3 ; j ++ )
        {
            this->coordinates( 3 * ( i + 1 ) + j , 0 ) = this->coordinates( j , 0 ) + tangent( j );
        }
    }

}









void MTOC2::getMatrixes_Sparse( MatrixXd &projection_Matrix )
{

	unsigned int coordinates = 3.0 * ( this->number_of_points + 1);
    MatrixXd projection_Matrix_local = MatrixXd::Zero( 3.0 * ( this->number_of_points + 1) , 3.0 * ( this->number_of_points + 1) );
    for( unsigned int i = 0 ; i < 3.0 * ( this->number_of_points + 1) ; i ++ )
    {
        projection_Matrix_local( i , i ) = 1.0;
    }


	unsigned int bonds = this->number_of_points;
	SparseMatrix<double> n_i_mu( bonds  , coordinates );

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList_n_i_mu;
	for( unsigned int bond = 1 ; bond <= this->number_of_points; bond ++ )
	{
		Vector3d tangent = this->tangent_to_center( bond );
        tangent = tangent / tangent.norm();

		for( unsigned int j = 0 ; j < 3 ; j ++ )
		{
				tripletList_n_i_mu.push_back(T( bond - 1 , j , - tangent( j ) ));
		}

        for( unsigned int j = 0 ; j < 3 ; j ++ )
		{
				tripletList_n_i_mu.push_back(T( bond - 1 , 3 * ( bond ) + j , tangent( j ) ));
		}

	}
	n_i_mu.setFromTriplets( tripletList_n_i_mu.begin() , tripletList_n_i_mu.end() );
    n_i_mu.makeCompressed();
    SparseMatrix<double> n_i_mu2Transpose = n_i_mu.transpose();
    n_i_mu2Transpose.makeCompressed();
    SparseMatrix<double> G_uv2 = n_i_mu * n_i_mu2Transpose;


	SparseMatrix<double> Identity( bonds , bonds );
	std::vector<T> tripletList_Identity;

    for( unsigned int index = 0 ; index < bonds ; index ++ )
	{
		tripletList_Identity.push_back(T( index , index , 1.0 ));
	}
    Identity.setFromTriplets( tripletList_Identity.begin() , tripletList_Identity.end() );
    Identity.makeCompressed();
    G_uv2.makeCompressed();


    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
	solver.compute( G_uv2 );

    SparseMatrix<double>  G_uv2_inv = solver.solve( Identity );

    MatrixXd tmp = G_uv2_inv * n_i_mu;
	tmp = n_i_mu2Transpose * tmp;
	projection_Matrix = projection_Matrix_local - tmp;

}


void MTOC2::getMatrixes_Sparse_2( MatrixXd &projection_Matrix )
{
	unsigned int coordinates = 3.0 * ( this->number_of_points + 1);
    	MatrixXd projection_Matrix_local = MatrixXd::Zero( 3.0 * ( this->number_of_points + 1) , 3.0 * ( this->number_of_points + 1) );

    	for( unsigned int i = 0 ; i < 3.0 * ( this->number_of_points + 1) ; i ++ )
    	{
        	projection_Matrix_local( i , i ) = 1.0;
    	}


	unsigned int bonds_1 = this->number_of_points;
    	unsigned int bonds_2 = this->number_of_points - 1;
    	unsigned int bonds = bonds_1 + bonds_2;

	SparseMatrix<double> n_i_mu( bonds  , coordinates );
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList_n_i_mu;
	for( unsigned int bond = 1 ; bond <= bonds_1; bond ++ )
	{
		Vector3d tangent = this->tangent_to_center( bond );
        	tangent = tangent / tangent.norm();

		for( unsigned int j = 0 ; j < 3 ; j ++ )
		{
				tripletList_n_i_mu.push_back(T( bond - 1 , j , - tangent( j ) ));
		}

        	for( unsigned int j = 0 ; j < 3 ; j ++ )
		{
				tripletList_n_i_mu.push_back(T( bond - 1 , 3 * ( bond ) + j , tangent( j ) ));
		}
	}


	for( unsigned int bond = 1 ; bond <= bonds_2; bond ++ )
	{
       	 	Vector3d tangent = this->tangent_from_next_bead( bond );
        	tangent = tangent / tangent.norm();
		for( unsigned int j = 0 ; j < 3 ; j ++ )
		{
			tripletList_n_i_mu.push_back(T( bonds_1 - 1 + bond , 3 * ( bond ) + j ,  tangent( j ) ));
            		tripletList_n_i_mu.push_back(T( bonds_1 - 1 + bond , 3 * ( bond + 1 ) + j , - tangent( j ) ));// +
        	}
	}

	n_i_mu.setFromTriplets( tripletList_n_i_mu.begin() , tripletList_n_i_mu.end() );
    	n_i_mu.makeCompressed();

    	SparseMatrix<double> n_i_mu2Transpose = n_i_mu.transpose();
    	n_i_mu2Transpose.makeCompressed();

    	MatrixXd G_uv2 = n_i_mu * n_i_mu2Transpose;
    	MatrixXd G_uv2_inv = G_uv2.inverse();
    	MatrixXd tmp = G_uv2_inv * n_i_mu;
	tmp = n_i_mu2Transpose * tmp;
	projection_Matrix = projection_Matrix_local - tmp;

}


void MTOC2::getMatrixes_Sparse_3( MatrixXd &projection_Matrix )
{
	    unsigned int coordinates = 3.0 * ( this->number_of_points + 1);
	    MatrixXd projection_Matrix_local = MatrixXd::Zero( 3.0 * ( this->number_of_points + 1) , 3.0 * ( this->number_of_points + 1) );

	    for( unsigned int i = 0 ; i < 3.0 * ( this->number_of_points + 1) ; i ++ )
	    {
		projection_Matrix_local( i , i ) = 1.0;
	    }


	    unsigned int bonds_1 = this->number_of_points;
	    unsigned int bonds_2 = this->number_of_points;
	    unsigned int bonds_3 = this->number_of_points / 2;
	    unsigned int bonds = bonds_1 + bonds_2;
	    unsigned int bonds_22 = bonds + bonds_3;



	    SparseMatrix<double> n_i_mu( bonds_22  , coordinates );
	    typedef Eigen::Triplet<double> T;
	    std::vector<T> tripletList_n_i_mu;
	for( unsigned int bond = 1 ; bond <= bonds_1; bond ++ )
	{
		Vector3d tangent = this->tangent_to_center( bond );
        	tangent = tangent / tangent.norm();

		for( unsigned int j = 0 ; j < 3 ; j ++ )
		{
				tripletList_n_i_mu.push_back(T( bond - 1 , j , - tangent( j ) ));
		}

        	for( unsigned int j = 0 ; j < 3 ; j ++ )
		{
				tripletList_n_i_mu.push_back(T( bond - 1 , 3 * ( bond ) + j , tangent( j ) ));
		}

	}


	for( unsigned int bond = 1 ; bond <= bonds_2; bond ++ )
	{
        	Vector3d tangent = this->tangent_from_next_bead( bond );
        	tangent = tangent / tangent.norm();
		if( bond < bonds_2 )
		{
			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				tripletList_n_i_mu.push_back(T( bonds_1 - 1 + bond , 3 * ( bond ) + j ,  tangent( j ) ));
		    		tripletList_n_i_mu.push_back(T( bonds_1 - 1 + bond , 3 * ( bond + 1 ) + j , - tangent( j ) ));// +
			}
		}
		if( bond == bonds_2 )
		{
			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				tripletList_n_i_mu.push_back(T( bonds_1 - 1 + bond , 3 * ( bond ) + j ,  tangent( j ) ));
		    		tripletList_n_i_mu.push_back(T( bonds_1 - 1 + bond , 3 + j , - tangent( j ) ));// +
			}
		}
	}

	for( unsigned int bond = 1 ; bond <= bonds_3; bond ++ )
	{
		unsigned int index_naproti = get_opposite_direct_point_index( bond );
        	Vector3d tangent = this->tangent_between_beads( bond , index_naproti );
        	tangent = tangent / tangent.norm();

		for( unsigned int j = 0 ; j < 3 ; j ++ )
		{
			tripletList_n_i_mu.push_back(T( bonds - 1 + bond , 3 * ( bond ) + j ,  tangent( j ) ));
            		tripletList_n_i_mu.push_back(T( bonds - 1 + bond , 3 * (  index_naproti ) + j , - tangent( j ) ));// +
        	}
	}


	    n_i_mu.setFromTriplets( tripletList_n_i_mu.begin() , tripletList_n_i_mu.end() ); 
    	    n_i_mu.makeCompressed();

    	    SparseMatrix<double> n_i_mu2Transpose = n_i_mu.transpose();
    	    n_i_mu2Transpose.makeCompressed();

    	    MatrixXd G_uv2 = n_i_mu * n_i_mu2Transpose;
    	    MatrixXd G_uv2_inv = G_uv2.inverse();
    	    MatrixXd tmp = G_uv2_inv * n_i_mu;
    	    tmp = n_i_mu2Transpose * tmp;
    	    projection_Matrix = projection_Matrix_local - tmp;

}



unsigned int MTOC2::get_random_point_on_opposite_side_without_bias( unsigned int point_index )
{
    if( point_index == 0 )
    {
        cout<<"point_index == 0"<<endl;
        cout<<"unsigned int MTOC2::get_random_point_on_opposite_side_without_bias( unsigned int point_index )"<<endl;
        cout<<"ERROR_ID = 74989613114"<<endl;
        throw("");
    }
    if( point_index > this->number_of_points )
    {
        cout<<"point_index > this->number_of_points"<<endl;
        cout<<"unsigned int MTOC2::get_random_point_on_opposite_side_without_bias( unsigned int point_index )"<<endl;
        cout<<"ERROR_ID = 98410463015614"<<endl;
        throw("");
    }

    unsigned int exactly_opposite_point = 0;
    if( point_index <= this->number_of_points / 2 )
    {
	exactly_opposite_point = point_index + this->number_of_points / 2;
    }
    else
    {
	exactly_opposite_point = point_index - this->number_of_points / 2;
    }

    std::uniform_int_distribution<int> distribution( ( -1.0 ) * MTOCparam::number_of_opposite_side ,  MTOCparam::number_of_opposite_side );
    unsigned int number_of_generator = omp_get_thread_num();
    int chosen_random_point = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] ) + exactly_opposite_point ;

    if ( chosen_random_point < 1 )
    {
	chosen_random_point = chosen_random_point + this->number_of_points;
    }
    else if ( chosen_random_point > 20 )
    {
	chosen_random_point = chosen_random_point - this->number_of_points;
    }
    return chosen_random_point;
}









void MTOC2::getMatrixes_Sparse_4( MatrixXd &projection_Matrix )
{
    	unsigned int coordinates = 3.0 * ( this->number_of_points + 1);
    	MatrixXd projection_Matrix_local = MatrixXd::Zero( 3.0 * ( this->number_of_points + 1) , 3.0 * ( this->number_of_points + 1) );

        for( unsigned int i = 0 ; i < 3.0 * ( this->number_of_points + 1) ; i ++ )
        {
             projection_Matrix_local( i , i ) = 1.0;
        }


	unsigned int bonds_1 = this->number_of_points;
    	unsigned int bonds_2 = this->number_of_points - 1;
    	unsigned int bonds_3 = this->number_of_points / 2;
    	unsigned int bonds_4 = this->number_of_points / 4 * 3;
    	unsigned int bonds = bonds_1 + bonds_2;
    	unsigned int bonds_22 = bonds + bonds_3;
    	unsigned int bonds_44 = bonds_22 + bonds_4;

	SparseMatrix<double> n_i_mu( bonds_44  , coordinates );

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList_n_i_mu;
	for( unsigned int bond = 1 ; bond <= bonds_1; bond ++ )
	{
		Vector3d tangent = this->tangent_to_center( bond );
        	tangent = tangent / tangent.norm();

		for( unsigned int j = 0 ; j < 3 ; j ++ )
		{
				tripletList_n_i_mu.push_back(T( bond - 1 , j , - tangent( j ) ));
		}

        	for( unsigned int j = 0 ; j < 3 ; j ++ )
		{
				tripletList_n_i_mu.push_back(T( bond - 1 , 3 * ( bond ) + j , tangent( j ) ));
		}
	}


	for( unsigned int bond = 1 ; bond <= bonds_2; bond ++ )
	{
        	Vector3d tangent = this->tangent_from_next_bead( bond );
        	tangent = tangent / tangent.norm();
		for( unsigned int j = 0 ; j < 3 ; j ++ )
		{
			tripletList_n_i_mu.push_back(T( bonds_1 - 1 + bond , 3 * ( bond ) + j ,  tangent( j ) ));
            		tripletList_n_i_mu.push_back(T( bonds_1 - 1 + bond , 3 * ( bond + 1 ) + j , - tangent( j ) ));// +
        	}
	}

	for( unsigned int bond = 1 ; bond <= bonds_3; bond ++ )
	{

        	Vector3d tangent = this->tangent_between_beads( bond , bond + this->number_of_points / 2 );
        	tangent = tangent / tangent.norm();

		for( unsigned int j = 0 ; j < 3 ; j ++ )
		{
			tripletList_n_i_mu.push_back(T( bonds - 1 + bond , 3 * ( bond ) + j ,  tangent( j ) ));
            		tripletList_n_i_mu.push_back(T( bonds - 1 + bond , 3 * (  bond + this->number_of_points / 2 ) + j , - tangent( j ) ));// +
        	}
	}

    for( unsigned int bond = 1 ; bond <= this->number_of_points / 4 * 3; bond ++ )
    {
        	Vector3d tangent = this->tangent_between_beads( bond , bond + this->number_of_points / 4 - 1 );
        	tangent = tangent / tangent.norm();
		for( unsigned int j = 0 ; j < 3 ; j ++ )
		{
			tripletList_n_i_mu.push_back(T( bonds_22 - 1 + bond , 3 * ( bond ) + j ,  tangent( j ) ));
            		tripletList_n_i_mu.push_back(T( bonds_22 - 1 + bond , 3 * ( bond + this->number_of_points / 4 - 1 ) + j , - tangent( j ) ));// +
        	}
    }





	n_i_mu.setFromTriplets( tripletList_n_i_mu.begin() , tripletList_n_i_mu.end() );
    	n_i_mu.makeCompressed();

    	SparseMatrix<double> n_i_mu2Transpose = n_i_mu.transpose();
    	n_i_mu2Transpose.makeCompressed();

    	MatrixXd G_uv2 = n_i_mu * n_i_mu2Transpose;
    	MatrixXd G_uv2_inv = G_uv2.inverse();
    	MatrixXd tmp = G_uv2_inv * n_i_mu;
	tmp = n_i_mu2Transpose * tmp;
	projection_Matrix = projection_Matrix_local - tmp;


}



void MTOC2::Euler_algorithm( MatrixXd& Forces )
{
    	MatrixXd original = this->get_coordinates();

	//Kronecker delta
	MatrixXd projectionMatrix = MatrixXd::Zero( 3 * ( this->number_of_points + 1) , 3 * ( this->number_of_points + 1) );
    	this->getMatrixes_Sparse_3( projectionMatrix );

    	double timeHalfStep =  sim_of_Cell::time_Step; // / 2.0
	MatrixXd projected_forces = projectionMatrix * Forces;

    	MatrixXd V_0 = ( projected_forces ) * ( 1.0 / this->effectiveFriction ); // dynein_ForceP
	MatrixXd R_half = original + V_0 * timeHalfStep;
	this->set_coordinates( R_half );

}







void MTOC2::oneStepMidStepAlgorithm_1_half( MatrixXd& Forces )
{
    	MatrixXd original = this->get_coordinates();

	//Kronecker delta
	MatrixXd projectionMatrix = MatrixXd::Zero( 3 * ( this->number_of_points + 1) , 3 * ( this->number_of_points + 1) );
    	this->getMatrixes_Sparse_3( projectionMatrix );

    	double timeHalfStep =  sim_of_Cell::time_Step_half; // / 2.0
	MatrixXd projected_forces = projectionMatrix * Forces;

    	MatrixXd V_0 = ( projected_forces ) * ( 1.0 / this->effectiveFriction ); // dynein_ForceP
	MatrixXd R_half = original + V_0 * timeHalfStep;
	this->set_coordinates( R_half );

}

void MTOC2::oneStepMidStepAlgorithm_2_half( MatrixXd& Forces , MatrixXd& coordinates_arg)
{
    //MatrixXd original = this->get_coordinates();
    //Kronecker delta
    MatrixXd projectionMatrix = MatrixXd::Zero( 3 * ( this->number_of_points + 1) , 3 * ( this->number_of_points + 1) );
    this->getMatrixes_Sparse_3( projectionMatrix );

    MatrixXd projected_forces = projectionMatrix * Forces;
    MatrixXd V_0 = ( projected_forces ) / this->effectiveFriction;

    MatrixXd R_final =  coordinates_arg + V_0 * sim_of_Cell::time_Step;
    this->set_coordinates( R_final );
}

void MTOC2::oneStepMidStepAlgorithm( )
{
    MatrixXd Forces = MatrixXd::Zero( 3 * ( this->number_of_points + 1 ) , 1 );
    double aaa = 10e-12;
    for( unsigned int i = 0 ; i < 11  ;  i ++)
    {
        Forces( 3 * i , 0 ) = aaa;
    }


    MatrixXd original_coordinates = this->get_coordinates();
    this->oneStepMidStepAlgorithm_1_half( Forces );
    this->oneStepMidStepAlgorithm_2_half( Forces , original_coordinates );

}

void MTOC2::MidStepAlgorithm( double Time_0 , double Time_1 )
{
    unsigned int number_of_steps = ( Time_1 - Time_0 ) / sim_of_Cell::time_Step;
    for( unsigned int i = 0 ; i < number_of_steps ; i ++ )
    {

        if( i % 10 == 0 )
        {
            cout<<"i = "<<i<<endl;
            this->controlMTOC( );
        }
        this->oneStepMidStepAlgorithm();
    }
}

double MTOC2::get_average_altitude()
{
    Vector3d center_orientation = this->get_center() / this->get_center().norm();
    double sum_magnitudes = 0;
    for( unsigned int point_id = 1 ; point_id <= this->get_number_of_points() ; point_id ++ )
    {
        Vector3d point = this->get_point( point_id );
        double cosinus = point.dot( center_orientation ) / ( center_orientation.norm() * point.norm() );

        double magnitude = point.norm();
        double projected_magnitude = cosinus * magnitude;
        sum_magnitudes = sum_magnitudes + projected_magnitude;
    }

    sum_magnitudes = sum_magnitudes + this->get_center().norm();
    double average = sum_magnitudes / ( double )( this->get_number_of_points() + 1 );
    return average;
}

void MTOC2::set_axis_of_rotation()
{
    Vector3d current_center = this->get_center();
    Vector3d current_center_tmp = current_center / current_center.norm();
    double angle = current_center_tmp.dot( this->original_MTOC_orientation ) / ( current_center_tmp.norm() * this->original_MTOC_orientation.norm() );

    angle = acos( angle );
    if( abs( angle - MTOCparam::MTOC_angle_rotation ) < 0.2  )
    {
        if( this->axis_of_rotation.size() == 0 )
        {
            Vector3d axis_of_rotation_arg = this->original_MTOC_orientation.cross( current_center_tmp );
            axis_of_rotation_arg = axis_of_rotation_arg / axis_of_rotation_arg.norm();
            this->axis_of_rotation.push_back( axis_of_rotation_arg );
        }
    }

}


void MTOC2::print_norm_tangent_to_center()
{
    for( unsigned int i = 1 ; i <= this->number_of_points ; i ++ )
    {
        Vector3d tangent = this->tangent_to_center( 1 );
        cout<<" i =  "<<i<<endl;
        cout<<" tangent.norm() =  "<<tangent.norm()<<endl;
    }
}


void MTOC2::print_points()
{
    cout<<"center"<<endl;
    cout<<this->get_center()<<endl;
    for( unsigned int i = 1 ; i <= this->number_of_points ; i ++ )
    {
        cout<<" Point_id = "<<i<<endl;
        cout<<this->get_point( i )<<endl;
    }

}

Vector3d MTOC2::get_side_center( unsigned int side )
{
    if( side >= 4 )
    {
        cout<<"side >= 4"<<endl;
        cout<<" MTOC2::get_side_center( unsigned int side ) "<<endl;
        throw("");
    }

    unsigned int number_of_points_in_side = this->get_number_of_points() / 4;

    Vector3d center_of_side( 0.0 , 0.0 , 0.0 );
    for( unsigned int point_counter = 0 ; point_counter < number_of_points_in_side ; point_counter ++ )
    {
        Vector3d point_position = this->get_point( 1 + number_of_points_in_side * side + point_counter );
        center_of_side = center_of_side + point_position;
    }

    center_of_side = center_of_side / 4.0;
    return center_of_side;
}










MatrixXd MTOC2::rotate_MTOC_from_original_orientations()
{
    Vector3d axis_z( 0.0 , 0.0 , 1 );
    double average_magnitude = get_average_altitude();
    MatrixXd new_coordinates = this->get_coordinates();
    MatrixXd upper_MTOC_coordinates = MatrixXd::Zero( 3 * ( this->get_number_of_points() + 1 ) , 1 );


    Vector3d instant_center = this->get_center();
    Vector3d instant_center_norm = instant_center / instant_center.norm();

    double angle = instant_center_norm.dot( this->original_MTOC_orientation );
    angle = acos( angle );

    Vector3d axis_of_rotation =  instant_center_norm.cross( this->original_MTOC_orientation );
    axis_of_rotation = axis_of_rotation / axis_of_rotation.norm();

    if( axis_of_rotation.norm() == 0 )
    {
        return this->coordinates;
    }

    Quaternion<double> q;
    q = AngleAxis<double>( angle , axis_of_rotation );
    if( angle == 0 )
    {
	return this->get_coordinates();
    }

    Vector3d upper_center = q *	instant_center_norm;

    upper_center( 2 ) = average_magnitude;
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
    	upper_MTOC_coordinates( dimension , 0 ) = upper_center( dimension );
    }

    //rotace nahoru
    for( unsigned int index = 1 ; index <= this->get_number_of_points(); index ++ )
    {
	    Vector3d point = this->get_point( index );
            Vector3d point_upper = q * point;
            for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
            {
                upper_MTOC_coordinates( 3 * index + dimension , 0 ) = point_upper( dimension );
            }
     }

    //vyrovnani vysky
    for( unsigned int index = 0 ; index <= this->get_number_of_points(); index ++ )
    {
	    upper_MTOC_coordinates( 3 * index + 2 , 0 ) = average_magnitude;
    }


    double uhel_prumer = 0;
    for( unsigned int index = 1 ; index <= this->get_number_of_points(); index ++ )
    {
            Vector3d orientation_local = this->get_original_orientation( index );
            Vector3d point_upper( 0.0 , 0.0 , 0.0 );
            for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
            {
                 point_upper( dimension ) = upper_MTOC_coordinates( 3 * index + dimension , 0 );
            }
	    point_upper = point_upper - upper_center;
	    double cos_uhel = point_upper.dot( orientation_local ) / ( point_upper.norm() *  orientation_local.norm() );
	    if( abs( cos_uhel ) > 1.0000 )
	    {
		if( cos_uhel > 0 )
		{
			cos_uhel = 1.0000;
		}
		if( cos_uhel < 0 )
		{
			cos_uhel = -1.0000;
		}
	    }
	    double uhel = acos( cos_uhel );
 	    uhel_prumer = uhel_prumer + uhel;
    }

    uhel_prumer = uhel_prumer / (double) this->get_number_of_points();
    Vector3d Axis_of_plane_rotation_1( 0.0 , 0.0 , 1 );
    Vector3d Axis_of_plane_rotation_2 = ( - 1.0 ) * Axis_of_plane_rotation_1;

    Quaternion<double> q_2;
    q_2 = AngleAxis<double>( uhel_prumer , Axis_of_plane_rotation_1 );

    Quaternion<double> q_3;
    q_3 = AngleAxis<double>( uhel_prumer , Axis_of_plane_rotation_2 );

    MatrixXd intermediate_coordinates_1 = MatrixXd::Zero( 3 * ( this->get_number_of_points() + 1 ) , 1 );
    MatrixXd intermediate_coordinates_2 = MatrixXd::Zero( 3 * ( this->get_number_of_points() + 1 ) , 1 );

    MatrixXd final_intermediate_coordinates = MatrixXd::Zero( 3 * ( this->get_number_of_points() + 1 ) , 1 );

    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
    	intermediate_coordinates_1( dimension , 0 ) = upper_center( dimension );
	intermediate_coordinates_2( dimension , 0 ) = upper_center( dimension );
    }

    for( unsigned int index = 1 ; index <= this->get_number_of_points(); index ++ )
    {
            Vector3d orientation_local = this->get_original_orientation( index );
	    Vector3d point_tmp_1 = q_2 * orientation_local * this->radius_of_MTOC;
            Vector3d point_tmp_2 = q_3 * orientation_local * this->radius_of_MTOC;
            point_tmp_1( 2 ) = average_magnitude;
            point_tmp_2( 2 ) = average_magnitude;
            for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
            {
    	         intermediate_coordinates_1( 3 * index + dimension , 0 ) = point_tmp_1( dimension );
	         intermediate_coordinates_2( 3 * index +  dimension , 0 ) = point_tmp_2( dimension );
            }
    }



    double norm_1 = coordinate_norm( intermediate_coordinates_1 - upper_MTOC_coordinates );
    double norm_2 = coordinate_norm( intermediate_coordinates_2 - upper_MTOC_coordinates );

    if( norm_1 >= norm_2 )
    {
	final_intermediate_coordinates = intermediate_coordinates_2;
    }
    else
    {
	final_intermediate_coordinates = intermediate_coordinates_1;
    }

    Quaternion<double> q_final;
    q_final = AngleAxis<double>( angle , ( - 1 ) * axis_of_rotation );

    Vector3d new_center = q_final * upper_center;
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
	    new_coordinates( dimension , 0 ) = new_center( dimension );
    }

    std::uniform_int_distribution<int> distribution( 1 ,  2 );
    unsigned int number_of_generator = omp_get_thread_num();
    unsigned int switch_local = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
    if( switch_local == 1 )
    {
        for( unsigned int index = 1 ; index <= 10 ; index ++ )
        {
              Vector3d point( 0.0 , 0.0 , 0.0 );
    	      for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	      {
	               point( dimension ) = final_intermediate_coordinates( 3 *  index + dimension , 0 );
    	      }
	      point = q_final * point;
    	      for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	      {
	               new_coordinates( 3 *  index + dimension , 0 ) = point( dimension );
    	      }
        }
    }
    else
    {
        for( unsigned int index = 11 ; index <= 20 ; index ++ )
        {
              Vector3d point( 0.0 , 0.0 , 0.0 );
    	      for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	      {
	               point( dimension ) = final_intermediate_coordinates( 3 *  index + dimension , 0 );
    	      }
	      point = q_final * point;

    	      for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	      {
	               new_coordinates( 3 *  index + dimension , 0 ) = point( dimension );
    	      }
        }
    }

    return new_coordinates;


}


MatrixXd MTOC2::rotate_MTOC_from_original_orientations_2()
{
    Vector3d instant_center = this->get_center();
    Vector3d instant_center_norm = instant_center / instant_center.norm();

    double angle = instant_center_norm.dot( this->original_MTOC_orientation );
    angle = acos( angle );

    Vector3d axis_of_rotation =  instant_center_norm.cross( this->original_MTOC_orientation );
    axis_of_rotation = axis_of_rotation / axis_of_rotation.norm();

    if( ( axis_of_rotation.norm() == 0 ) || ( angle == 0 ) )
    {
        return this->coordinates;
    }


    Quaternion<double> q;
    q = AngleAxis<double>( angle , axis_of_rotation );

    double average_magnitude = get_average_altitude();
    MatrixXd upper_MTOC_coordinates = MatrixXd::Zero( 3 * ( this->get_number_of_points() + 1 ) , 1 );

    Vector3d upper_center = q *	instant_center_norm;
    upper_center( 2 ) = average_magnitude;
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
    	upper_MTOC_coordinates( dimension , 0 ) = upper_center( dimension );
    }

    for( unsigned int index = 1 ; index <= this->get_number_of_points(); index ++ )
    {
	    Vector3d point = this->get_point( index );
            Vector3d point_upper = q * point;
            for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
            {
                upper_MTOC_coordinates( 3 * index + dimension , 0 ) = point_upper( dimension );
            }
     }


    for( unsigned int index = 0 ; index <= this->get_number_of_points(); index ++ )
    {
	    upper_MTOC_coordinates( 3 * index + 2 , 0 ) = average_magnitude;
    }

   return upper_MTOC_coordinates;

}






void MTOC2::resize_from_originals()
{
     //return;

    Vector3d axis_z( 0.0 , 0.0 , 1 );
    double average_magnitude = get_average_altitude();
    MatrixXd new_coordinates = this->get_coordinates();
    MatrixXd upper_MTOC_coordinates = MatrixXd::Zero( 3 * ( this->get_number_of_points() + 1 ) , 1 );


    Vector3d instant_center = this->get_center();
    Vector3d instant_center_norm = instant_center / instant_center.norm();

    double angle = instant_center_norm.dot( this->original_MTOC_orientation );
    angle = acos( angle );

    Vector3d axis_of_rotation =  instant_center_norm.cross( this->original_MTOC_orientation );
    axis_of_rotation = axis_of_rotation / axis_of_rotation.norm();

    if( axis_of_rotation.norm() == 0 )
    {
	return;
    }

    Quaternion<double> q;
    q = AngleAxis<double>( angle , axis_of_rotation );
    if( angle == 0 )
    {
	return;
    }

    Vector3d upper_center = q *	instant_center_norm;

    upper_center( 2 ) = average_magnitude;
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
    	upper_MTOC_coordinates( dimension , 0 ) = upper_center( dimension );
    }

    for( unsigned int index = 1 ; index <= this->get_number_of_points(); index ++ )
    {
	    Vector3d point = this->get_point( index );
            Vector3d point_upper = q * point;
            for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
            {
                upper_MTOC_coordinates( 3 * index + dimension , 0 ) = point_upper( dimension );
            }
     }
    for( unsigned int index = 0 ; index <= this->get_number_of_points(); index ++ )
    {
	    upper_MTOC_coordinates( 3 * index + 2 , 0 ) = average_magnitude;
    }


    double uhel_prumer = 0;
    for( unsigned int index = 1 ; index <= this->get_number_of_points(); index ++ )
    {
            Vector3d orientation_local = this->get_original_orientation( index );
            Vector3d point_upper( 0.0 , 0.0 , 0.0 );
            for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
            {
                 point_upper( dimension ) = upper_MTOC_coordinates( 3 * index + dimension , 0 );
            }
	    point_upper = point_upper - upper_center;
	    double cos_uhel = point_upper.dot( orientation_local ) / ( point_upper.norm() *  orientation_local.norm() );
	    if( abs( cos_uhel ) > 1.0000 )
	    {
		if( cos_uhel > 0 )
		{
			cos_uhel = 1.0000;
		}
		if( cos_uhel < 0 )
		{
			cos_uhel = -1.0000;
		}
	    }
	    double uhel = acos( cos_uhel );

 	    uhel_prumer = uhel_prumer + uhel;
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////

    uhel_prumer = uhel_prumer / (double) this->get_number_of_points();
    Vector3d Axis_of_plane_rotation_1( 0.0 , 0.0 , 1 );
    Vector3d Axis_of_plane_rotation_2 = ( - 1.0 ) * Axis_of_plane_rotation_1;

    Quaternion<double> q_2;
    q_2 = AngleAxis<double>( uhel_prumer , Axis_of_plane_rotation_1 );

    Quaternion<double> q_3;
    q_3 = AngleAxis<double>( uhel_prumer , Axis_of_plane_rotation_2 );

    MatrixXd intermediate_coordinates_1 = MatrixXd::Zero( 3 * ( this->get_number_of_points() + 1 ) , 1 );
    MatrixXd intermediate_coordinates_2 = MatrixXd::Zero( 3 * ( this->get_number_of_points() + 1 ) , 1 );

    MatrixXd final_intermediate_coordinates = MatrixXd::Zero( 3 * ( this->get_number_of_points() + 1 ) , 1 );

    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
    	intermediate_coordinates_1( dimension , 0 ) = upper_center( dimension );
	intermediate_coordinates_2( dimension , 0 ) = upper_center( dimension );
    }

    for( unsigned int index = 1 ; index <= this->get_number_of_points(); index ++ )
    {
            Vector3d orientation_local = this->get_original_orientation( index );
	    Vector3d point_tmp_1 = q_2 * orientation_local * this->radius_of_MTOC;
            Vector3d point_tmp_2 = q_3 * orientation_local * this->radius_of_MTOC;
            point_tmp_1( 2 ) = average_magnitude;
            point_tmp_2( 2 ) = average_magnitude;
            for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
            {
    	         intermediate_coordinates_1( 3 * index + dimension , 0 ) = point_tmp_1( dimension );
	         intermediate_coordinates_2( 3 * index +  dimension , 0 ) = point_tmp_2( dimension );
            }
    }


    double norm_1 = coordinate_norm( intermediate_coordinates_1 - upper_MTOC_coordinates );
    double norm_2 = coordinate_norm( intermediate_coordinates_2 - upper_MTOC_coordinates );

    if( norm_1 >= norm_2 )
    {
	final_intermediate_coordinates = intermediate_coordinates_2;
    }
    else
    {
	final_intermediate_coordinates = intermediate_coordinates_1;
    }

    Quaternion<double> q_final;
    q_final = AngleAxis<double>( angle , ( - 1 ) * axis_of_rotation );

    Vector3d new_center = q_final * upper_center;
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
	    new_coordinates( dimension , 0 ) = new_center( dimension );
    }


    std::uniform_int_distribution<int> distribution( 1 ,  2 );
    unsigned int number_of_generator = omp_get_thread_num();
    unsigned int switch_local = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
    if( switch_local == 1 )
    {
        for( unsigned int index = 1 ; index <= 10 ; index ++ )
        {
              Vector3d point( 0.0 , 0.0 , 0.0 );
    	      for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	      {
	               point( dimension ) = final_intermediate_coordinates( 3 *  index + dimension , 0 );
    	      }
	      point = q_final * point;
    	      for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	      {
	               new_coordinates( 3 *  index + dimension , 0 ) = point( dimension );
    	      }
        }
    }
    else
    {
        for( unsigned int index = 11 ; index <= 20 ; index ++ )
        {
              Vector3d point( 0.0 , 0.0 , 0.0 );
    	      for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	      {
	               point( dimension ) = final_intermediate_coordinates( 3 *  index + dimension , 0 );
    	      }
	      point = q_final * point;

    	      for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	      {
	               new_coordinates( 3 *  index + dimension , 0 ) = point( dimension );
    	      }
        }
    }

    this->coordinates = new_coordinates;

}




void MTOC2::getRandomForces(  double timeStep , MatrixXd &randomForces )
{

	std::normal_distribution<> distribution{0,1};
        unsigned int number_of_generator = omp_get_thread_num();


	if( ( randomForces.rows() !=  3 * ( this->get_number_of_points() + 1 ) ) || ( ( randomForces.cols() != 1 ) ) )
	{
		cout<<"( bendingForce.getNN() != 3 * this->numberOfPoints ) || ( ( bendingForce.getNN() != 1 ) ) in number_of_points::getBendingForces( MatDoub bendingForce )"<<endl;
		throw("");
	}

	unsigned int coordinates = 3 * ( this->get_number_of_points() + 1 );
	double root = 1.0 * sqrt( 2.0 * sim_of_Cell::boltzmann_constant * sim_of_Cell::Temperature * this->get_effectiveFriction(  ) / timeStep );

	if( sim_of_Cell::random_force_switch_MTOC == true )
	{

	}
	else
	{
        	root = 0;
	}
	for( unsigned int i = 0 ; i < coordinates ; i ++ )
	{
		double normal_number = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
		randomForces( i , 0 ) = normal_number * root;
	}

}





MatrixXd MTOC2::set_coordinates_on_same_altitude()
{
    MatrixXd new_coordinates = MatrixXd::Zero( 3 * ( this->get_number_of_points() + 1 )  , 1 );
    double average_magnitude = get_average_altitude();

    Vector3d center = this->get_center();
    center = center / center.norm() * average_magnitude;
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
       new_coordinates( dimension , 0 ) = center( dimension );
    }


    for( unsigned int index = 1 ; index <= 20 ; index ++ )
    {
        Vector3d point = this->get_point( index );
        Vector3d new_point = point / point.norm() * average_magnitude;

        for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
        {
             new_coordinates( 3 * index + dimension , 0 ) = new_point( dimension );
        }
    }

   return new_coordinates;
}



double MTOC2::magnitude_of_the_force( MatrixXd& Forces  )
{
    Vector3d force3d( 0.0 , 0.0 , 0.0 );
    double sum = 0;
    for( unsigned int i = 0 ; i < Forces.rows() / 3 ; i ++ )
    {
        for( unsigned int j = 0 ; j < 3 ; j ++ )
        {
            force3d( j ) = Forces( 3 * i + j , 0 );
        }
        sum = sum + force3d.norm();
    }
    return sum;
}

unsigned int MTOC2::get_opposite_point_index( unsigned int point_index )
{
    if( this->get_number_of_points() < MTOCparam::boundary_MTOC_points )
    {
        cout<<"this->get_number_of_points() < MTOCparam::boundary_MTOC_points"<<endl;
        cout<<"unsigned int MTOC2::get_opposite_point_index( unsigned int point_index )"<<endl;
        cout<<"ERROR_ID = 31668938"<<endl;
        throw("");
    }
    if( point_index == 0 )
    {
        cout<<"point_index == 0"<<endl;
        cout<<"unsigned int MTOC2::get_opposite_point_index( unsigned int point_index )"<<endl;
        cout<<"ERROR_ID = 897616481648"<<endl;
        throw("");
    }
    if( point_index > this->number_of_points )
    {
        cout<<"point_index > this->number_of_points"<<endl;
        cout<<"unsigned int MTOC2::get_opposite_point_index( unsigned int point_index )"<<endl;
        cout<<"ERROR_ID = 8755816938472"<<endl;
        throw("");
    }


    unsigned int opposite_index;
    if( point_index <= this->number_of_points / 2 )
    {
        opposite_index = point_index + this->number_of_points / 2;
    }
    else
    {
       opposite_index = point_index - this->number_of_points / 2;
    }

    return opposite_index;
}

unsigned int MTOC2::get_opposite_direct_point_index( unsigned int point_index )
{
    if( this->get_number_of_points() < MTOCparam::boundary_MTOC_points )
    {
        cout<<"this->get_number_of_points() < MTOCparam::boundary_MTOC_points"<<endl;
        cout<<"unsigned int MTOC2::get_opposite_direct_point_index( unsigned int point_index )"<<endl;
        cout<<"ERROR_ID = 61568412615"<<endl;
        throw("");
    }
    if( point_index == 0 )
    {
        cout<<"point_index == 0"<<endl;
        cout<<"unsigned int MTOC2::get_opposite_direct_point_index( unsigned int point_index )"<<endl;
        cout<<"ERROR_ID = 898976450364"<<endl;
        throw("");
    }
    if( point_index > this->number_of_points )
    {
        cout<<"point_index > this->number_of_points"<<endl;
        cout<<"unsigned int MTOC2::get_opposite_direct_point_index( unsigned int point_index )"<<endl;
        cout<<"ERROR_ID = 145414754156875468545435"<<endl;
        throw("");
    }

    unsigned int opposite_index;
    unsigned int point_index_tmp = point_index - 1;

    if( point_index <= this->number_of_points / 2 )
    {
        unsigned int lower_quarter_index = point_index_tmp / ( this->get_number_of_points() / 4 );
        unsigned int central_index = lower_quarter_index * ( this->get_number_of_points() / 4 ) + MTOCparam::MTOC_center_index;
        unsigned int difference = central_index - point_index_tmp;
        unsigned int opposite_central_index = central_index + this->number_of_points / 2;
        opposite_index = opposite_central_index + difference + 1;
    }
    else
    {
        unsigned int upper_quarter_index = point_index_tmp / ( this->get_number_of_points() / 4 );
        unsigned int central_index = upper_quarter_index * ( this->get_number_of_points() / 4 ) + MTOCparam::MTOC_center_index;
        unsigned int difference = central_index - point_index_tmp;
        unsigned int opposite_central_index = central_index - this->number_of_points / 2;
        opposite_index = opposite_central_index + difference + 1;
    }

    return opposite_index;

}



unsigned int MTOC2::get_random_point_on_opposite_side( unsigned int point_index )
{
    if( point_index == 0 )
    {
        cout<<"point_index == 0"<<endl;
        cout<<"unsigned int MTOC2::get_random_point_on_opposite_side( unsigned int point_index )"<<endl;
        cout<<"ERROR_ID = 898976450364"<<endl;
        throw("");
    }
    if( point_index > this->number_of_points )
    {
        cout<<"point_index > this->number_of_points"<<endl;
        cout<<"unsigned int MTOC2::get_random_point_on_opposite_side( unsigned int point_index )"<<endl;
        cout<<"ERROR_ID = 145414754156875468545435"<<endl;
        throw("");
    }

    unsigned int point_index_2 =  point_index - 1;
    unsigned int side_id = ( int ) point_index_2 / 5;
    unsigned int opposite_side_id = 0;
    if( side_id < 2 )
    {
         opposite_side_id = side_id + 2;
    }
    else
    {
         opposite_side_id = side_id - 2;
    }

    unsigned int lower_boundary = ( int ) opposite_side_id * this->number_of_points / 4.0 + 1;
    unsigned int higher_boundary = lower_boundary + this->number_of_points / 4.0;

    std::uniform_int_distribution<int> distribution( lower_boundary ,  higher_boundary - 1 );
    unsigned int number_of_generator = omp_get_thread_num();
    unsigned int return_index = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );


    return return_index;
}









MTOC2::~MTOC2() {
	// TODO Auto-generated destructor stub
}

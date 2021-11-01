#include "Microtubule.h"

//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS


Microtubule::Microtubule()
{

	this->microtubule_id = 10000;
	this->polygon = 10000;
	this->restDistancePoints = sim_of_Cell::resting_distance;
        this->MTOC_point = 1;
	this->kappa =  2.2e-23 / restDistancePoints;// 2.2e-25 v pohode  2.2e-26
	this->numberOfPoints = sim_of_Cell::MicrotubulePoints;
	this->dydein_index = 0.0;
	this->effective_friction = 1e-7;
	this->IS_position_catching = Vector3d( 0.0 , 0.0 , 0.0 );


	if( this->numberOfPoints < 1  )
	{
		cout<<" this->numberOfPoints < 1 in Microtubule::Microtubule() "<<endl;
		throw("");
	}

	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
	for( unsigned int i = 0 ; i < this->numberOfPoints ; i ++ )
	{
		this->coordinates( 3 * i , 0 ) = i * this->restDistancePoints;
		this->coordinates( 3 * i + 1 , 0 ) = 0.0;
		this->coordinates( 3 * i + 2 , 0 ) = 0.0;
	}

    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }
    this->set_bending_matrix( );
}



Microtubule::Microtubule( unsigned int ID )
{
	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  2.2e-23 / restDistancePoints;// 2.2e-25 v pohode  2.2e-26
	this->numberOfPoints = sim_of_Cell::MicrotubulePoints;
	this->microtubule_id = ID;
        this->MTOC_point = 1;
	this->polygon = 10000;
	this->dydein_index = 0.0;
	this->effective_friction = this->calculateEffectiveFriction_Howard();
	this->IS_position_catching = Vector3d( 0.0 , 0.0 , 0.0 );
	if( this->numberOfPoints < 1  )
	{
		cout<<" this->numberOfPoints < 1 in Microtubule::Microtubule( unsigned int ID ) "<<endl;
		throw("");
	}

	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );

	for( unsigned int i = 0 ; i < this->numberOfPoints ; i ++ )
	{
		this->coordinates( 3 * i , 0 ) = i * this->restDistancePoints;
		this->coordinates( 3 * i + 1 , 0 ) = 0.0;
		this->coordinates( 3 * i + 2 , 0 ) = 0.0;
	}
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }
    this->set_bending_matrix( );
}



Microtubule::Microtubule( Vector3d MTOC , unsigned int ID )
{
	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  2.2e-23 / restDistancePoints;// 2.2e-25 v pohode  2.2e-26
	this->numberOfPoints = sim_of_Cell::MicrotubulePoints;
	this->microtubule_id = ID;
    this->MTOC_point = 1;
	this->polygon = 10000;
	this->dydein_index = 0.0;
	this->effective_friction = this->calculateEffectiveFriction_Howard();
	this->IS_position_catching = Vector3d( 0.0 , 0.0 , 0.0 );

	if( this->numberOfPoints < 1  )
	{
		cout<<" this->numberOfPoints < 1 in Microtubule::Microtubule( Vector3d MTOC , unsigned int ID )"<<endl;
		throw("");
	}

	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );

	for( unsigned int i = 0 ; i < this->numberOfPoints ; i ++ )
	{
		this->coordinates( 3 * i , 0 ) = MTOC( 0 ) +  i * this->restDistancePoints;
		this->coordinates( 3 * i + 1 , 0 ) = MTOC( 1 ) + 0.0;
		this->coordinates( 3 * i + 2 , 0 ) = MTOC( 2 ) + 0.0;
	}
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }
    this->set_bending_matrix( );
}



Microtubule::Microtubule( Vector3d MTOC , Vector3d first_Point , unsigned int ID )
{
	if( ( MTOC - first_Point ).norm() == 0 )
	{
		cout<<"( MTOC - first_Point ).norm() == 0 in Microtubule::Microtubule( Vector3d MTOC , Vector3d first_Point , unsigned int ID )"<<endl;
		throw("");
	}


	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  2.2e-23 / restDistancePoints;// 2.2e-25 v pohode  2.2e-26
	this->numberOfPoints = sim_of_Cell::MicrotubulePoints;
	this->microtubule_id = ID;
	this->polygon = 10000;
    this->MTOC_point = 1;
	this->IS_position_catching = Vector3d( 0.0 , 0.0 , 0.0 );
	this->dydein_index = 0.0;
	this->effective_friction = this->calculateEffectiveFriction_Howard();

	if( this->numberOfPoints < 1  )
	{
		cout<<" this->numberOfPoints < 1 in Microtubule::Microtubule( Vector3d MTOC , Vector3d first_Point , unsigned int ID )"<<endl;
		throw("");
	}



	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
	Vector3d direction( ( first_Point[ 0 ] - MTOC[ 0 ]  ) , ( first_Point[ 1 ]  - MTOC[ 1 ] )  , 0 );
	direction = direction / direction.norm();
	Vector3d axis( direction( 1 ) , - direction( 0 )  , 0 );

	//constructing Quaterion - object to rotate vectors
	//the same Quaterion will be used for all rotation
	Quaternion<double> q;
	q = AngleAxis<double>( sim_of_Cell::angle , axis );



	//coordinates of first point:: equals first_Point
	this->coordinates( 0 , 0 ) = first_Point( 0 );
	this->coordinates( 1 , 0 ) = first_Point( 1 );
	this->coordinates( 2 , 0 ) = first_Point( 2 );


	if( this->numberOfPoints > 1 )
	{
		//first tangent - between first_Point and second bead of microtubule
		Vector3d tangent( this->restDistancePoints * direction( 0 ) , this->restDistancePoints * direction( 1 ) , 0.0 );

		this->coordinates( 3 , 0 ) = first_Point( 0 ) + tangent( 0 );
		this->coordinates( 4 , 0 ) = first_Point( 1 ) + tangent( 1 );
		this->coordinates( 5 , 0 ) = first_Point( 2 ) + tangent( 2 );


		//Others tangents are created here- there magnitude is still the same
		//Only difference is orientation
		//that is because "spider-like" structure of microtubule cytoskeleton - see movies
		for( unsigned int i = 2 ; i < this->numberOfPoints ; i ++ )
		{
			tangent = q * tangent;
			this->coordinates( 3 * i , 0 ) = this->coordinates( 3 * ( i - 1 ) , 0 ) + tangent( 0 );
			this->coordinates( 3 * i + 1 , 0 ) = this->coordinates( 3 * ( i - 1 ) + 1 , 0 ) + tangent( 1 );
			this->coordinates( 3 * i + 2 , 0 ) = this->coordinates( 3 * ( i - 1 ) + 2 , 0 ) + tangent( 2 );
		}

	}
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }
    this->set_bending_matrix( );
}



Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID , unsigned int poly )
{
	//Creates microtubule using two first beads
	//It assumes that the microtubule will be in the direction of segment between two first points
	//third_Point = second_Point + ( second_Point - first_Point ) until this->number of points

	if( ( second_Point - first_Point ).norm() == 0 )
	{
		cout<<"( second_Point - first_Point ).norm() == 0 in Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID , 666 )"<<endl;
		throw("");
	}

	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  2.2e-23 / restDistancePoints;
	this->numberOfPoints = sim_of_Cell::MicrotubulePoints;
	this->microtubule_id = ID;
	this->polygon = poly;
    this->MTOC_point = 1;
	this->dydein_index = 0.0;
	this->effective_friction = this->calculateEffectiveFriction_Howard();
	this->IS_position_catching = Vector3d( 666.0 , 0.0 , 0.0 );

	if( abs( ( second_Point - first_Point ).norm() / this->restDistancePoints - 1.0 ) > 1e-8)
	{
		cout<<"abs( ( second_Point - first_Point ).norm() / this->restDistancePoints ) > 1e-8 in Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID , 666 )"<<endl;
		throw("");
	}


	if( this->numberOfPoints < 2  )
	{
		cout<<" this->numberOfPoints < 2 in Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID )"<<endl;
		throw("");
	}


	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );

	Vector3d direction( second_Point[ 0 ] - first_Point[ 0 ] , second_Point[ 1 ] - first_Point[ 1 ]  , second_Point[ 2 ] - first_Point[ 2 ] );
	direction = direction / direction.norm();

	this->coordinates( 0 , 0 ) = first_Point( 0 );
	this->coordinates( 1 , 0 ) = first_Point( 1 );
	this->coordinates( 2 , 0 ) = first_Point( 2 );

	this->coordinates( 3 , 0 ) = second_Point( 0 );
	this->coordinates( 4 , 0 ) = second_Point( 1 );
	this->coordinates( 5 , 0 ) = second_Point( 2 );


	if( this->numberOfPoints > 2 )
	{
		Vector3d tangent( this->restDistancePoints * direction( 0 ) , this->restDistancePoints * direction( 1 ) , 0.0 );
		//Others tangents are created here- there magnitude is still the same
		//Only difference is orientation
		for( unsigned int i = 2 ; i < this->numberOfPoints ; i ++ )
		{
			this->coordinates( 3 * i , 0 ) = this->coordinates( 3 * ( i - 1 ) , 0 ) + tangent( 0 );
			this->coordinates( 3 * i + 1 , 0 ) = this->coordinates( 3 * ( i - 1 ) + 1 , 0 ) + tangent( 1 );
			this->coordinates( 3 * i + 2 , 0 ) = this->coordinates( 3 * ( i - 1 ) + 2 , 0 ) + tangent( 2 );
		}
	}
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }
    this->set_bending_matrix( );
}


Microtubule::Microtubule( Vector3d first_Point , Vector3d orientation , unsigned int ID , unsigned int poly , double a_axis , double b_axis )
{
	//Creates microtubule using two first beads
	//It assumes that the microtubule will be in the direction of segment between two first points
	//third_Point = second_Point + ( second_Point - first_Point ) until this->number of points

	if( this->numberOfPoints < 2  )
	{
		cout<<" this->numberOfPoints < 2 in Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID )"<<endl;
		throw("");
	}

	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  2.2e-23 / restDistancePoints;
	this->numberOfPoints = sim_of_Cell::MicrotubulePoints;
	this->microtubule_id = ID;
	this->polygon = poly;
    this->MTOC_point = 1;
	this->dydein_index = 0.0;
	this->effective_friction = this->calculateEffectiveFriction_Howard();

	//IS: it has no meaning in unconstrained cell
	this->IS_position_catching = Vector3d( 666.0 , 0.0 , 0.0 );

	bool index = confirm_inner_ellipsoid( first_Point , a_axis , b_axis );




	if( index == 0  )
	{
		cout<<"a_axis = "<<a_axis<<endl;
		cout<<"b_axis = "<<b_axis<<endl;
		cout<<"first_Point = "<<first_Point<<endl;
		cout<<"index = "<<index<<endl;
		cout<<"first_Point is outside of ellipsod in Microtubule ... double a_axis , double b_axis  ... ERROR_ID=6789468416351684"<<endl;
		throw("");
	}

	Vector3d direction( orientation( 0 ) , orientation( 1 )  , 0 );
	direction = direction / direction.norm();
	//axis that I will use to rotate the vectors
	Vector3d axis_rotation( direction( 1 ) , - direction( 0 )  , 0 );


	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
	//coordinates of first point:: equals first_Point
	this->coordinates( 0 , 0 ) = first_Point( 0 );
	this->coordinates( 1 , 0 ) = first_Point( 1 );
	this->coordinates( 2 , 0 ) = first_Point( 2 );

	Vector3d last_point = first_Point;
	Vector3d last_orientation = orientation;
	for( unsigned int index = 1 ; index < this->numberOfPoints ; index ++ )
	{
		Vector3d tmp_point = last_point + orientation * this->restDistancePoints;
		if( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 1 )
		{
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
		}
		else
		{
			//I aplly quaternion to rotate
			double rotation_angle = - 3.14159265358979 / 180.0;
			Quaternion<double> q;
			q = AngleAxis<double>( rotation_angle , axis_rotation );
			do
			{
				//rotation will be done by one degree
				orientation = q * orientation;
				tmp_point = last_point + orientation * this->restDistancePoints;

			} while ( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 0 );
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
			last_orientation = orientation;
		}

	}
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }
    this->set_bending_matrix( );

}






Microtubule::Microtubule( Vector3d first_Point , Vector3d orientation , unsigned int ID , unsigned int poly , double a_axis , double b_axis , unsigned int number_of_points  )
{
	//Creates microtubule using two first beads
	//It assumes that the microtubule will be in the direction of segment between two first points
	//third_Point = second_Point + ( second_Point - first_Point ) until this->number of points

	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  2.2e-23 / restDistancePoints;
	this->numberOfPoints = number_of_points;
	this->microtubule_id = ID;
	this->polygon = poly;
    this->MTOC_point = 1;
	this->dydein_index = 0.0;
	//this->effective_friction = this->calculateEffectiveFriction( );
	this->effective_friction = this->calculateEffectiveFriction_Howard();

	if( this->numberOfPoints < 2  )
	{
		cout<<" this->numberOfPoints < 2 in Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID )"<<endl;
		throw("");
	}


	//IS: it has no meaning in unconstrained cell
	this->IS_position_catching = Vector3d( 666.0 , 0.0 , 0.0 );

	//confirm that first bead lies in the ellipsoid - if not throw error
	bool index = confirm_inner_ellipsoid( first_Point , a_axis , b_axis );




	if( index == 0  )
	{
		cout<<"a_axis = "<<a_axis<<endl;
		cout<<"b_axis = "<<b_axis<<endl;
		cout<<"first_Point = "<<first_Point<<endl;
		cout<<"index = "<<index<<endl;
		cout<<"first_Point is outside of ellipsod in Microtubule ... double a_axis , double b_axis  ... ERROR_ID=689451765"<<endl;
		throw("");
	}


	Vector3d direction( orientation( 0 ) , orientation( 1 )  , 0 );
	direction = direction / direction.norm();

	Vector3d axis_rotation( direction( 1 ) , - direction( 0 )  , 0 );


	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
	this->coordinates( 0 , 0 ) = first_Point( 0 );
	this->coordinates( 1 , 0 ) = first_Point( 1 );
	this->coordinates( 2 , 0 ) = first_Point( 2 );

	Vector3d last_point = first_Point;
	Vector3d last_orientation = orientation;
	for( unsigned int index = 1 ; index < this->numberOfPoints ; index ++ )
	{
		Vector3d tmp_point = last_point + orientation * this->restDistancePoints;
		if( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 1 )
		{
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
		}
		else
		{

			double rotation_angle = - 3.14159265358979 / 180.0;
			Quaternion<double> q;
			q = AngleAxis<double>( rotation_angle , axis_rotation );
			do
			{
				//rotation will be done by one degree
				orientation = q * orientation;
				tmp_point = last_point + orientation * this->restDistancePoints;

			} while ( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 0 );
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
			last_orientation = orientation;
		}

	}
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }
    this->set_bending_matrix( );
}


Microtubule::Microtubule( Vector3d first_Point , Vector3d orientation , unsigned int ID , unsigned int poly , unsigned int side_arg , unsigned int MTOC_point_arg , unsigned int number_of_points  )
{

    	//Creates microtubule using two first beads
	//It assumes that the microtubule will be in the direction of segment between two first points
	//third_Point = second_Point + ( second_Point - first_Point ) until this->number of points

	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  2.2e-23 / restDistancePoints;
	this->numberOfPoints = number_of_points;
	this->microtubule_id = ID;
	this->polygon = poly;
    this->side = side_arg;
    this->MTOC_point = MTOC_point_arg;
	this->dydein_index = 0.0;
	//this->effective_friction = this->calculateEffectiveFriction( );
this->effective_friction = this->calculateEffectiveFriction_Howard();
	if( this->numberOfPoints < 2  )
	{
		cout<<" this->numberOfPoints < 2 in Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID )"<<endl;
		throw("");
	}

	this->IS_position_catching = Vector3d( 666.0 , 0.0 , 0.0 );
    double a_axis = Cell_parametres::A_AXIS;
    double b_axis = Cell_parametres::B_AXIS;
	bool bool_index = confirm_inner_ellipsoid( first_Point , a_axis , b_axis );

	if( bool_index == 0  )
	{
		cout<<"a_axis = "<<a_axis<<endl;
		cout<<"b_axis = "<<b_axis<<endl;
		cout<<"first_Point = "<<first_Point<<endl;
		cout<<"bool_index = "<<bool_index<<endl;
		cout<<"first_Point is outside of ellipsod in Microtubule ... double a_axis , double b_axis  ... ERROR_ID=46461343154127"<<endl;
		throw("");
	}


	//direction from MTOC to first_point
	//this is direction of microtubule in x - y plane - does not change for entire microtubule
	Vector3d direction( orientation( 0 ) , orientation( 1 )  , 0 );
	direction = direction / direction.norm();
	//axis that I will use to rotate the vectors
	Vector3d axis_rotation( direction( 1 ) , - direction( 0 )  , 0 );



	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
	this->coordinates( 0 , 0 ) = first_Point( 0 );
	this->coordinates( 1 , 0 ) = first_Point( 1 );
	this->coordinates( 2 , 0 ) = first_Point( 2 );


	//last point and orientation will be changed as we go from
	Vector3d last_point = first_Point;
	Vector3d last_orientation = orientation;
	for( unsigned int index = 1 ; index < this->numberOfPoints ; index ++ )
	{
		Vector3d tmp_point = last_point + orientation * this->restDistancePoints;
		if( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 1 )
		{
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
		}
		else
		{
			//I aplly quaternion to rotate
			double rotation_angle = - 3.14159265358979 / 180.0;
			Quaternion<double> q;
			q = AngleAxis<double>( rotation_angle , axis_rotation );
			do
			{
				//rotation will be done by one degree
				orientation = q * orientation;
				tmp_point = last_point + orientation * this->restDistancePoints;

			} while ( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 0 );
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
			last_orientation = orientation;
		}

	}
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }
    this->set_bending_matrix( );

}


Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , Vector3d orientation , unsigned int ID , unsigned int poly , unsigned int side_arg , unsigned int MTOC_point_arg , unsigned int second_MTOC_point_arg , unsigned int number_of_points  )
{

    if( MTOC_point_arg == 0 )
    {
        cout<<"MTOC_point_arg == 0"<<endl;
        cout<<"ERROR_ID = 689645684316464"<<endl;
        throw("");
    }
    if( second_MTOC_point_arg == 0 )
    {
        cout<<"MTOC_point_arg == 0"<<endl;
        cout<<"ERROR_ID = 8987625148688"<<endl;
        throw("");
    }

    	//Creates microtubule using two first beads
	//It assumes that the microtubule will be in the direction of segment between two first points
	//third_Point = second_Point + ( second_Point - first_Point ) until this->number of points

	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  2.2e-23 / restDistancePoints;
	this->numberOfPoints = number_of_points;
	this->microtubule_id = ID;
	this->polygon = poly;
    this->side = side_arg;
    this->MTOC_point = MTOC_point_arg;
    this->second_MTOC_point = second_MTOC_point_arg;
	this->dydein_index = 0.0;
	//this->effective_friction = this->calculateEffectiveFriction( );
	this->effective_friction = this->calculateEffectiveFriction_Howard();
    this->restDistancePoints_first = ( first_Point - second_Point ).norm();

    if( abs( ( first_Point - second_Point ).norm() - this->restDistancePoints_first ) > 1e-10  )
    {
        cout<<"( first_Point - second_Point ).norm() = "<<( first_Point - second_Point ).norm()<<endl;
        cout<<"this->restDistancePoints_first = "<<this->restDistancePoints_first<<endl;
        cout<<"abs( ( first_Point - second_Point ).norm() - this->restDistancePoints_first ) > 1e-10"<<endl;
        cout<<"ERROR_ID = 61543554435"<<endl;
        throw("");
    }



	if( this->numberOfPoints < 2  )
	{
		cout<<" this->numberOfPoints < 2 in Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID )"<<endl;
		throw("");
	}

	this->IS_position_catching = Vector3d( 666.0 , 0.0 , 0.0 );
    double a_axis = Cell_parametres::A_AXIS;
    double b_axis = Cell_parametres::B_AXIS;
	bool bool_index = confirm_inner_ellipsoid( first_Point , a_axis , b_axis );

	if( bool_index == 0  )
	{
		cout<<"a_axis = "<<a_axis<<endl;
		cout<<"b_axis = "<<b_axis<<endl;
		cout<<"first_Point = "<<first_Point<<endl;
		cout<<"bool_index = "<<bool_index<<endl;
		cout<<"first_Point is outside of ellipsod in Microtubule ... double a_axis , double b_axis  ... ERROR_ID=46461343154127"<<endl;
		throw("");
	}

	Vector3d direction( orientation( 0 ) , orientation( 1 )  , 0 );
	direction = direction / direction.norm();
	Vector3d axis_rotation( direction( 1 ) , - direction( 0 )  , 0 );


	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
	this->coordinates( 0 , 0 ) = first_Point( 0 );
	this->coordinates( 1 , 0 ) = first_Point( 1 );
	this->coordinates( 2 , 0 ) = first_Point( 2 );

    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        	this->coordinates( 3 + dimension , 0 ) = second_Point( dimension );
    }

	Vector3d last_point = second_Point;
	Vector3d last_orientation = orientation;
	for( unsigned int index = 2 ; index < this->numberOfPoints ; index ++ )
	{
		Vector3d tmp_point = last_point + orientation * this->restDistancePoints;
		if( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 1 )
		{
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
		}
		else
		{
			//I aplly quaternion to rotate
			double rotation_angle = - 3.14159265358979 / 180.0;
			Quaternion<double> q;
			q = AngleAxis<double>( rotation_angle , axis_rotation );
			do
			{
				//rotation will be done by one degree
				orientation = q * orientation;
				tmp_point = last_point + orientation * this->restDistancePoints;

			} while ( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 0 );
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
			last_orientation = orientation;
		}

	}
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }
    this->set_bending_matrix( );

}


Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , Vector3d orientation , unsigned int ID , unsigned int poly , unsigned int side_arg , unsigned int MTOC_point_arg , unsigned int second_MTOC_point_arg , unsigned int number_of_points , bool rovna_mikrotubula )
{

    if( MTOC_point_arg == 0 )
    {
        cout<<"MTOC_point_arg == 0"<<endl;
        cout<<"ERROR_ID = 689645684316464"<<endl;
        throw("");
    }
    if( second_MTOC_point_arg == 0 )
    {
        cout<<"MTOC_point_arg == 0"<<endl;
        cout<<"ERROR_ID = 8987625148688"<<endl;
        throw("");
    }

    	//Creates microtubule using two first beads
	//It assumes that the microtubule will be in the direction of segment between two first points
	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  2.2e-23 / restDistancePoints;
	this->numberOfPoints = number_of_points;
	this->microtubule_id = ID;
	this->polygon = poly;
    this->side = side_arg;
    this->MTOC_point = MTOC_point_arg;
    this->second_MTOC_point = second_MTOC_point_arg;
    this->dydein_index = 0.0;
    this->effective_friction = this->calculateEffectiveFriction_Howard();
    this->restDistancePoints_first = 2.0 * MTOCparam::MTOC_radius;

    if( abs( ( first_Point - second_Point ).norm() - this->restDistancePoints_first ) > 1e-10  )
    {
        cout<<"( first_Point - second_Point ).norm() = "<<( first_Point - second_Point ).norm()<<endl;
        cout<<"this->restDistancePoints_first = "<<this->restDistancePoints_first<<endl;
        cout<<"abs( ( first_Point - second_Point ).norm() - this->restDistancePoints_first ) > 1e-10"<<endl;
        cout<<"ERROR_ID = 154355454435"<<endl;
        throw("");
    }



	if( this->numberOfPoints < 2  )
	{
		cout<<" this->numberOfPoints < 2 in Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID )"<<endl;
		throw("");
	}

	this->IS_position_catching = Vector3d( 666.0 , 0.0 , 0.0 );
    double a_axis = Cell_parametres::A_AXIS;
    double b_axis = Cell_parametres::B_AXIS;
	bool bool_index = confirm_inner_ellipsoid( first_Point , a_axis , b_axis );

	if( bool_index == 0  )
	{
		cout<<"a_axis = "<<a_axis<<endl;
		cout<<"b_axis = "<<b_axis<<endl;
		cout<<"first_Point = "<<first_Point<<endl;
		cout<<"bool_index = "<<bool_index<<endl;
		cout<<"first_Point is outside of ellipsod in Microtubule ... double a_axis , double b_axis  ... ERROR_ID=46461343154127"<<endl;
		throw("");
	}

	Vector3d direction( orientation( 0 ) , orientation( 1 )  , 0 );
	direction = direction / direction.norm();
	Vector3d axis_rotation( direction( 1 ) , - direction( 0 )  , 0 );


	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
	this->coordinates( 0 , 0 ) = first_Point( 0 );
	this->coordinates( 1 , 0 ) = first_Point( 1 );
	this->coordinates( 2 , 0 ) = first_Point( 2 );

    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        	this->coordinates( 3 + dimension , 0 ) = second_Point( dimension );
    }




	//last point and orientation will be changed as we go from
	Vector3d last_point = second_Point;
	Vector3d last_orientation = orientation;
	for( unsigned int index = 2 ; index < this->numberOfPoints ; index ++ )
	{
		Vector3d tmp_point = last_point + orientation * this->restDistancePoints;
		if( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 1 )
		{
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
		}
		else
		{

			tmp_point = last_point + orientation * this->restDistancePoints;

			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
			last_orientation = orientation;
		}

	}
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int segment_couter = 0 ; segment_couter < this->numberOfPoints - 1 ; segment_couter ++ )
    {
        this->lenght_of_tangents( segment_couter , 0 ) = this->getTangent2( segment_couter ).norm();
    }

    this->set_bending_matrix( );
}





Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , Vector3d orientation , unsigned int ID , unsigned int poly , unsigned int side_arg , unsigned int MTOC_point_arg , unsigned int second_MTOC_point_arg , unsigned int number_of_points , double first_bead_distance )
{



    if( MTOC_point_arg == 0 )
    {
        cout<<"MTOC_point_arg == 0"<<endl;
        cout<<"ERROR_ID = 689645684316464"<<endl;
        throw("");
    }
    if( second_MTOC_point_arg == 0 )
    {
        cout<<"MTOC_point_arg == 0"<<endl;
        cout<<"ERROR_ID = 8987625148688"<<endl;
        throw("");
    }

    	//Creates microtubule using two first beads
	//It assumes that the microtubule will be in the direction of segment between two first points
	//third_Point = second_Point + ( second_Point - first_Point ) until this->number of points
	this->restDistancePoints = sim_of_Cell::resting_distance;
	this->kappa =  2.2e-23 / restDistancePoints;
	this->numberOfPoints = number_of_points;
	this->microtubule_id = ID;
	this->polygon = poly;
    this->side = side_arg;
    this->MTOC_point = MTOC_point_arg;
    this->second_MTOC_point = second_MTOC_point_arg;
	this->dydein_index = 0.0;
	this->effective_friction = this->calculateEffectiveFriction_Howard();
    this->restDistancePoints_first = first_bead_distance;


    if( abs( ( first_Point - second_Point ).norm() - this->restDistancePoints_first ) > 1e-10  )
    {
        cout<<"( first_Point - second_Point ).norm() = "<<( first_Point - second_Point ).norm()<<endl;
        cout<<"this->restDistancePoints_first = "<<this->restDistancePoints_first<<endl;
        cout<<"abs( ( first_Point - second_Point ).norm() - this->restDistancePoints_first ) > 1e-10"<<endl;
        cout<<"ERROR_ID = 61543554545"<<endl;
        throw("");
    }



	if( this->numberOfPoints < 2  )
	{
		cout<<" this->numberOfPoints < 2 in Microtubule::Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID )"<<endl;
		throw("");
	}

	this->IS_position_catching = Vector3d( 666.0 , 0.0 , 0.0 );
    double a_axis = Cell_parametres::A_AXIS;
    double b_axis = Cell_parametres::B_AXIS;

	bool bool_index = confirm_inner_ellipsoid( first_Point , a_axis , b_axis );

	if( bool_index == 0  )
	{
		cout<<"a_axis = "<<a_axis<<endl;
		cout<<"b_axis = "<<b_axis<<endl;
		cout<<"first_Point = "<<first_Point<<endl;
		cout<<"bool_index = "<<bool_index<<endl;
		cout<<"first_Point is outside of ellipsod in Microtubule ... double a_axis , double b_axis  ... ERROR_ID=46461343154127"<<endl;
		throw("");
	}


	//direction from MTOC to first_point
	Vector3d direction( orientation( 0 ) , orientation( 1 )  , 0 );
	direction = direction / direction.norm();
	//axis that I will use to rotate the vectors
	Vector3d axis_rotation( direction( 1 ) , - direction( 0 )  , 0 );


	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
	this->coordinates( 0 , 0 ) = first_Point( 0 );
	this->coordinates( 1 , 0 ) = first_Point( 1 );
	this->coordinates( 2 , 0 ) = first_Point( 2 );

    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        	this->coordinates( 3 + dimension , 0 ) = second_Point( dimension );
    }




	//last point and orientation will be changed as we go from
	Vector3d last_point = second_Point;
	Vector3d last_orientation = orientation;
	for( unsigned int index = 2 ; index < this->numberOfPoints ; index ++ )
	{
		Vector3d tmp_point = last_point + orientation * this->restDistancePoints;
		if( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 1 )
		{
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}

			last_point = tmp_point;
		}
		else
		{
			//I aplly quaternion to rotate
			double rotation_angle = - 3.14159265358979 / 180.0;
			Quaternion<double> q;
			q = AngleAxis<double>( rotation_angle , axis_rotation );
			do
			{
				//rotation will be done by one degree
				orientation = q * orientation;
				tmp_point = last_point + orientation * this->restDistancePoints;


			} while ( confirm_inner_ellipsoid( tmp_point , a_axis , b_axis ) == 0 );
			for( unsigned int index_axis = 0 ; index_axis < 3 ; index_axis ++ )
			{
				this->coordinates( 3 * index + index_axis , 0 ) = tmp_point( index_axis );
			}
			last_point = tmp_point;
			last_orientation = orientation;
		}

	}


    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    this->set_lenght_of_tangents();
    this->set_bending_matrix( );

}




















Microtubule::Microtubule( const Microtubule & tmp )
{

	this->restDistancePoints = tmp.getRestDist();
	this->kappa =  tmp.getKappa();
	this->numberOfPoints = tmp.getNumberOfPoints();
	this->microtubule_id = tmp.microtubule_id;
	this->effective_friction = tmp.effective_friction;
	this->IS_position_catching = tmp.IS_position_catching;

	this->coordinates = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
	this->coordinates = tmp.getCoordinates();
	this->polygon = tmp.polygon;
	this->dydein_index = tmp.dydein_index;
	this->dydein_index = tmp.effective_friction;


    this->Dynein_motors_2 = tmp.Dynein_motors_2;
    this->MTOC_point = tmp.MTOC_point;
    this->second_MTOC_point = tmp.second_MTOC_point;
    this->restDistancePoints_first = tmp.restDistancePoints_first;
    this->lenght_of_tangents = tmp.lenght_of_tangents;
    this->b_MATRIX = tmp.b_MATRIX;
    this->Kronecker_Delta = tmp.Kronecker_Delta;

}




void Microtubule::setCoordinates( MatrixXd &coordinatesTmp )
{
	if( ( coordinatesTmp.rows() !=  3 * this->numberOfPoints ) || ( ( coordinatesTmp.cols() != 1 ) ) )
	{
		cout<<"( ( coordinatesTmp.rows() !=  3 * this->numberOfPoints ) || ( ( coordinatesTmp.cols() != 1 ) ) in Microtubule::setCoordinates( MatDoub setCoordinates )"<<endl;
		cout<<"this->numberOfPoints = "<<this->numberOfPoints<<endl;
		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
		throw("");
	}
	this->coordinates = coordinatesTmp;
}

void Microtubule::setPoint( unsigned int number_of_point ,  Vector3d position )
{

    if( number_of_point >= this->getNumberOfPoints() )
    {
        cout<<"( number_of_point >= this->getNumberOfPoints() )"<<endl;
        cout<<"Microtubule::setPoint( unsigned int number_of_point ,  Vector3d position )"<<endl;
        throw("");
    }

    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        this->coordinates( 3 * number_of_point + dimension , 0 ) = position( dimension );
    }
}







MatrixXd Microtubule::getCoordinates()  const
{
	return this->coordinates;
}

double Microtubule::getRestDist() const
{
	return this->restDistancePoints;
}

double Microtubule::getKappa()  const
{
	return this->kappa;
}









Vector3d Microtubule::get_IS_position_catching( ) const
{
	return this->IS_position_catching;
}

double Microtubule::get_effective_friction()
{
	return this->effective_friction;
}

double Microtubule::get_effective_friction_whole_microtubule()
{
	return this->effective_friction * ( (double) this->numberOfPoints );
}




unsigned int Microtubule::getID() const
{
	return this->microtubule_id;
}


unsigned int Microtubule::getSide()  const
{
    return this->side;
}



unsigned int Microtubule::get_polygon_number() const
{
	return this->polygon;
}

unsigned int Microtubule::get_MTOC_point()
{
    return this->MTOC_point;
}


unsigned int Microtubule::get_MTOC_opposite_point()
{
    return this->second_MTOC_point;
}


unsigned int Microtubule::get_dynein_index() const
{
	return this->dydein_index;
}


void Microtubule::alter_flexural_rigidity( double ratio )
{
    this->kappa = this->kappa * ratio;
}

void Microtubule::set_dynein_index( unsigned int new_value )
{
	this->dydein_index = new_value;
}


void Microtubule::set_IS_position_catching( Vector3d IS_vector )
{
	this->IS_position_catching = IS_vector;
}



Microtubule& Microtubule::operator=( const Microtubule &tmp )
{
	if( this == &tmp )
	{
		cout<<" this == &tmp  in Microtubule::operator=( const Microtubule &tmp )"<<endl;
		throw("");
	}

	this->restDistancePoints = tmp.getRestDist();
	this->kappa =  tmp.getKappa();// 2.2e-25 v pohode  2.2e-26
	this->numberOfPoints = tmp.getNumberOfPoints();
	this->coordinates = tmp.getCoordinates();
	this->dydein_index = tmp.get_dynein_index();
	this->IS_position_catching = tmp.get_IS_position_catching();
	this->effective_friction = tmp.effective_friction;

    this->microtubule_id = tmp.microtubule_id;
    this->side = tmp.side;
    this->polygon = tmp.get_polygon_number();
	this->MTOC_point = tmp.MTOC_point;
    this->second_MTOC_point = tmp.second_MTOC_point;
    this->restDistancePoints_first = tmp.restDistancePoints_first;
    this->lenght_of_tangents = tmp.lenght_of_tangents;
    this->b_MATRIX = tmp.b_MATRIX;
    this->Kronecker_Delta = tmp.Kronecker_Delta;
    return *this;
}



unsigned int Microtubule::getNumberOfPoints()  const
{
	return this->numberOfPoints;
}



double Microtubule::calculateEffectiveFriction(  )
{
	//radius of the cell is always 1.25e-8
	//Stokes law is used
	//FrictionForce = 3 * PI * viscosity * diametr * K
	//K - dynamic shape factor

	//Drag on Nonspherical Objects - David Leith - Aerosol Science and Technology

	double lenghtOfSegment = this->getRestDist();
	double viscosity =  sim_of_Cell::viscosity;  // 9.3 * 10.0

    double radius = 1.25e-8;
    double tmp = 3.0 / 4.0 * radius * radius * lenghtOfSegment;
    double d_v = 2.0 * std::pow( tmp , 1.0 / 3.0);
    double d_u = 2.0 * sqrt( 2.0 * radius * lenghtOfSegment / 3.14159265 );
    double d_s = sqrt( 2.0 * radius * radius + radius * lenghtOfSegment );			//tady dvojka neni umyslne
    double K = 1.0 / 3.0 * d_u / d_v + 2.0 / 3.0 * d_s / d_v;

    double effectiveFriction = 3.0 * sim_of_Cell::PI * sim_of_Cell::multiply_friction_constant * viscosity * d_v * K;
    return effectiveFriction;
}

double Microtubule::calculateEffectiveFriction_Howard(  )
{
    double lenghtOfSegment = this->getRestDist();
    double viscosity =  sim_of_Cell::viscosity;  // 9.3 * 10.0
    double radius = 1.25e-8;

    double nominator = 4.0 * sim_of_Cell::PI * viscosity * lenghtOfSegment;
    double denominator = log( lenghtOfSegment / ( 2.0 * radius ) ) + 0.84;

    double effectiveFriction = nominator / denominator;
    return effectiveFriction;
}





Vector3d Microtubule::getPoint( unsigned int index )  const
{
	if( index >= this->numberOfPoints  )
	{

		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
		cout<<" this->getNumberOfPoints() = "<<this->getNumberOfPoints()<<endl;
		cout<<"index = "<<index<<endl;
		cout<<"index >= this->numberOfPoints = in Vector3d Microtubule::getPoint( unsigned int index ) "<<endl;
		throw("");
	}


	Vector3d result( 0.0 , 0.0 , 0.0 );
	for( unsigned int i = 0 ; i < 3 ; i ++ )
	{
		result( i ) = this->coordinates( 3 * ( index ) + i , 0 );
	}

	return result;
}


Vector3d Microtubule::get_last_Point( )  const
{
    return this->getPoint( this->getNumberOfPoints() - 1 );
}






Vector3d Microtubule::getTangent2( unsigned int index )  const
{
	if( index >= this->numberOfPoints - 1 )
	{

        	cout<<"this->dydein_index = "<<this->dydein_index<<endl;
        	cout<<" number of points in microtubule = "<<this->getNumberOfPoints()<<endl;
        	cout<<"index = "<<index<<endl;
		cout<<"index >= this->numberOfPoints - 1 in Microtubule::getTangent2( unsigned int index )"<<endl;
        	int ERROR_ID = 1411;
		cout<<" Microtubule ERROR_ID = "<<ERROR_ID<<endl;
        	double lenght = this->get_lenght_of_microtubule();
        	cout<<"............................................."<<endl;
        	cout<<"lenght = "<<lenght<<endl;
        	if( ( this->dydein_index == 9 ) || ( this->dydein_index == 20 ) )
        	{
            		cout<<"condition"<<endl;
            		cout<<"this->Dynein_motors_2.size() = "<<this->Dynein_motors_2.size()<<endl;
            		for( unsigned int i = 0 ; i < this->Dynein_motors_2.size() ; i ++ )
            		{
                		std::pair < Vector3d , double > pair_tmp = Dynein_motors_2.at( i );
                		double abscissa = std::get<1>( pair_tmp );
                		if( abscissa > lenght )
                		{
                    			cout<<"abscissa = "<<abscissa<<endl;
                		}
           		}
        	}

        	throw ERROR_ID;

	}

	Vector3d tangent( 0.0 , 0.0 , 0.0 );
	for( unsigned int j = 0 ; j < 3 ; j ++ )
	{
		tangent( j ) = this->coordinates( 3 * ( index + 1 ) + j , 0 ) - this->coordinates( 3 * index + j , 0 );
	}
	return tangent;
}



Vector3d Microtubule::get_last_Tangent( )  const
{
	return getTangent2( this->getNumberOfPoints() - 2 );
}







double Microtubule::get_lenght_of_microtubule() const
{
    double lenght_of_microtubule = 0;
    for( unsigned int segment_id = 0 ; segment_id < this->getNumberOfPoints() - 1 ; segment_id ++ )
    {
        Vector3d tangent;
        try
        {
            tangent = this->getTangent2( segment_id );
        }
        catch( int e )
        {
            cout<<"double Microtubule::get_lenght_of_microtubule() const"<<endl;
        }
        lenght_of_microtubule = lenght_of_microtubule + tangent.norm();
    }
    return lenght_of_microtubule;
}


double Microtubule::get_lenght_of_microtubule_outside_MTOC() const
{
    if( this->getNumberOfPoints() <= 2 )
    {
	return 0;
    }
    double lenght_of_microtubule = 0;
    for( unsigned int segment_id = 1 ; segment_id < this->getNumberOfPoints() - 1 ; segment_id ++ )
    {
        Vector3d tangent;
        try
        {
            tangent = this->getTangent2( segment_id );
        }
        catch( int e )
        {
            cout<<"double Microtubule::get_lenght_of_microtubule() const"<<endl;
        }
        lenght_of_microtubule = lenght_of_microtubule + tangent.norm();
    }
    return lenght_of_microtubule;
}




double Microtubule::get_distance_to_lower_bead_with_index( unsigned int index ) const
{
    if( index < this->numberOfPoints - 1 )
    {
        double lenght = 0;
        for( unsigned int counter = 0 ; counter < index ; counter ++ )
        {
            lenght = lenght + this->getTangent2( counter ).norm();
        }
        return lenght;
    }
    else if( index == this->numberOfPoints - 1 )
    {
        return this->get_lenght_of_microtubule();
    }
    else
    {
        cout<<"index > this->numberOfPoints - 1 "<<endl;
        cout<<"double Microtubule::get_distance_to_lower_bead_with_index( unsigned int index ) const"<<endl;
        unsigned int ERROR_ID = 981468;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
        throw ERROR_ID;
    }



}


void Microtubule::set_lenght_after_catching( double lenght_after_catching_arg )
{
    this->lenght_after_catching = lenght_after_catching_arg;
}

double Microtubule::get_lenght_after_catching(  )
{
    return this->lenght_after_catching;
}



Vector3d Microtubule::der_scalar_index( unsigned int scalar , unsigned int index)
{
	Vector3d force( 0.0 , 0.0 , 0.0 );

	if( ( index < 0 ) || ( index >= this->numberOfPoints ) )
	{
		cout<<"( index < 0 ) || ( index >= this->getNumberOfPoints()) in der_scalar_index( unsigned int scalar , unsigned int index)"<<endl;
		throw("");
	}
	if( ( scalar < 0 ) || ( scalar > this->numberOfPoints - 3 ) )
	{
		cout<<"scalar = "<<scalar<<endl;
		cout<<"( scalar < 0 ) || ( scalar >= this->getNumberOfPoints() - 3 ) in der_scalar_index( unsigned int scalar , unsigned int index)"<<endl;
		throw("");
	}

	if( this->numberOfPoints < 3 )
	{
		return force;
	}

	if( ( index - scalar  > 2 ) &&  ( scalar > index ) )
	{
		return force;
	}


	if(  index == scalar  )										//this is border case of bead of minimum number
	{
		Vector3d first = this->getTangent2( index );
		double magThisVec = first.norm();
		first = first * ( 1.0 / magThisVec );

		Vector3d second = this->getTangent2( index + 1 );
		double magNextVec = second.norm();
		second = second * ( 1.0 / magNextVec );

		force( 0 ) = ( - second( 0 ) + first( 0 ) * first.dot( second ) ) / magThisVec;
		force( 1 ) = ( - second( 1 ) + first( 1 ) * first.dot( second ) ) / magThisVec;
		force( 2 ) = ( - second( 2 ) + first( 2 ) * first.dot( second ) ) / magThisVec;
	}
	else if ( index == scalar + 2 )
	{
		Vector3d first = this->getTangent2( index - 2 );		
		double magThisVec = first.norm();
		first = first * ( 1.0 / magThisVec );


		Vector3d second = this->getTangent2( index - 1);
		double magNextVec = second.norm();
		second = second * ( 1.0 / magNextVec );		

		force( 0 ) = ( first( 0 ) - second( 0 ) * first.dot( second ) ) / magNextVec;
		force( 1 ) = ( first( 1 ) - second( 1 ) * first.dot( second ) ) / magNextVec;
		force( 2 ) = ( first( 2 ) - second( 2 ) * first.dot( second ) ) / magNextVec;
	}


	else if ( index == scalar + 1)
	{
//		cout<<1<<endl;
		Vector3d first = this->getTangent2( index - 1 );
		double magThisVec = first.norm();
		first = first * ( 1.0 / magThisVec );


		Vector3d second = this->getTangent2( index );
		double magNextVec = second.norm();
		second = second * ( 1.0 / magNextVec );		

		force( 0 ) = ( second( 0 ) - first( 0 ) * first.dot( second ) ) / magThisVec - ( first( 0 ) - second( 0 ) * first.dot( second ) ) / magNextVec;;
		force( 1 ) = ( second( 1 ) - first( 1 ) * first.dot( second ) ) / magThisVec - ( first( 1 ) - second( 1 ) * first.dot( second ) ) / magNextVec;
		force( 2 ) = ( second( 2 ) - first( 2 ) * first.dot( second ) ) / magThisVec - ( first( 2 ) - second( 2 ) * first.dot( second ) ) / magNextVec;
	}

	return force;




}


void Microtubule::getN_i_mu( MatrixXd &n_i_mu )
{
	for( unsigned int i = 0 ; i < this->numberOfPoints - 1 ; i ++ )
	{
		Vector3d tangent = this->getTangent2( i );
		tangent = tangent * ( 1.0 / tangent.norm() );
		for( unsigned int j = 0 ; j < 6 ; j ++ )
		{
			int tmp = 3 * i;
			if( j < 3)
			{
				n_i_mu( tmp + j , i ) = - tangent( j );
			}
			else
			{
				n_i_mu( tmp + j , i ) =  tangent( j - 3 );
			}
		}
	}
}












void Microtubule::getMatrixes( MatrixXd &inv_Matrix , MatrixXd &projection_Matrix )
{
	//Brownian Dynamics algorithm for bead-rod semiflexible chain with anisotropic friction
	//Montesi Morse Pasquali
	//Projection matrix is expressed in equation 16
	//projectionMatrix has to be Kronecker Delta			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	MatrixXd scalars =  MatrixXd::Zero(this->numberOfPoints - 2 , 1 );
	for( unsigned int i = 0 ; i < this->numberOfPoints - 2 ; i ++ )
	{
		scalars( i , 0 ) = this->getTangent2( i ).dot(this->getTangent2( i + 1 ) ) / ( this->getTangent2( i ).norm() * this->getTangent2( i + 1 ).norm() );
	}

	unsigned int coordinates = 3.0 * this->numberOfPoints;
	unsigned int bonds = this->numberOfPoints - 1;
        MatrixXd n_i_mu2 = MatrixXd::Zero( coordinates , bonds );
        getN_i_mu( n_i_mu2 );
        MatrixXd n_i_mu2Transpose = n_i_mu2.transpose();


	MatrixXd G_uv2 = MatrixXd::Zero( bonds , bonds );
	for( unsigned int index = 0 ; index < bonds ; index ++ )
	{
		G_uv2( index , index ) = 2;
		if( index != 0 )
		{
			G_uv2( index - 1 , index ) = - scalars( index - 1 , 0 );
			G_uv2( index , index - 1 ) = G_uv2( index - 1 , index );
		}
	}

	inv_Matrix = G_uv2.inverse();
	MatrixXd tmp = inv_Matrix * n_i_mu2Transpose;
	tmp = n_i_mu2 * tmp;


	projection_Matrix = projection_Matrix - tmp;
}











void Microtubule::getMatrixes_2( MatrixXd &inv_Matrix , MatrixXd &projection_Matrix )
{
	//Brownian Dynamics algorithm for bead-rod semiflexible chain with anisotropic friction
	//Montesi Morse Pasquali
	//Projection matrix is expressed in equation 16
	//projectionMatrix has to be Kronecker Delta			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	MatrixXd scalars =  MatrixXd::Zero(this->numberOfPoints - 2 , 1 );
	for( unsigned int i = 0 ; i < this->numberOfPoints - 2 ; i ++ )
	{
		scalars( i , 0 ) = this->getTangent2( i ).dot(this->getTangent2( i + 1 ) ) / ( this->getTangent2( i ).norm() * this->getTangent2( i + 1 ).norm() );
	}

	unsigned int coordinates = 3.0 * this->numberOfPoints;
	unsigned int bonds = this->numberOfPoints - 1;
        MatrixXd n_i_mu2 = MatrixXd::Zero( coordinates , bonds );
        getN_i_mu( n_i_mu2 );

        MatrixXd n_i_mu2Transpose = n_i_mu2.transpose();


	SparseMatrix<double> G_uv2( bonds , bonds );
	std::vector<T> tripletList_G_uv2;



	SparseMatrix<double> Identity( bonds , bonds );
	std::vector<T> tripletList_Identity;
        Identity.setIdentity ();

	for( unsigned int index = 0 ; index < bonds ; index ++ )
	{
		tripletList_G_uv2.push_back(T( index , index , 2.0 ));
		if( index != 0 )
		{

			tripletList_G_uv2.push_back(T( index - 1 , index , - scalars( index - 1 , 0 ) ));
			tripletList_G_uv2.push_back(T( index , index - 1 , - scalars( index - 1 , 0 ) ));
		}
	}

	G_uv2.setFromTriplets( tripletList_G_uv2.begin() , tripletList_G_uv2.end() );




	Eigen::SimplicialLLT <Eigen::SparseMatrix<double> > solver;//
	solver.compute( G_uv2 );

	SparseMatrix<double>  G_uv2_inv = solver.solve( Identity );
	MatrixXd tmp = G_uv2 * n_i_mu2Transpose;
	tmp = n_i_mu2 * tmp;


	projection_Matrix = projection_Matrix - tmp;
}



void Microtubule::getMatrixes_3( MatrixXd &projection_Matrix )
{
	//Brownian Dynamics algorithm for bead-rod semiflexible chain with anisotropic friction
	//Montesi Morse Pasquali
	//Projection matrix is expressed in equation 16
	//projectionMatrix has to be Kronecker Delta			



	unsigned int coordinates = 3.0 * this->numberOfPoints;
	unsigned int bonds = this->numberOfPoints - 1;


	SparseMatrix<double> n_i_mu2( coordinates , bonds );
        std::vector<T> tripletList_n_i_mu;
	for( unsigned int i = 0 ; i < this->numberOfPoints - 1 ; i ++ )
	{
		Vector3d tangent = this->getTangent2( i );
		tangent = tangent * ( 1.0 / tangent.norm() );
		for( unsigned int j = 0 ; j < 6 ; j ++ )
		{
			int tmp = 3 * i;
			if( j < 3)
			{
				tripletList_n_i_mu.push_back(T( tmp + j , j , - tangent( j ) ) );
			}
			else
			{
				tripletList_n_i_mu.push_back(T( tmp + j , i , - tangent( j - 3 ) ) );
			}
		}
	}
	n_i_mu2.setFromTriplets( tripletList_n_i_mu.begin() , tripletList_n_i_mu.end() );
        SparseMatrix<double> n_i_mu2Transpose = n_i_mu2.transpose();



	MatrixXd scalars =  MatrixXd::Zero(this->numberOfPoints - 2 , 1 );
	for( unsigned int i = 0 ; i < this->numberOfPoints - 2 ; i ++ )
	{
		scalars( i , 0 ) = this->getTangent2( i ).dot(this->getTangent2( i + 1 ) ) / ( this->getTangent2( i ).norm() * this->getTangent2( i + 1 ).norm() );
	}

	SparseMatrix<double> G_uv2( bonds , bonds );
	std::vector<T> tripletList_G_uv2;



	SparseMatrix<double> Identity( bonds , bonds );
	std::vector<T> tripletList_Identity;
        Identity.setIdentity ();

	for( unsigned int index = 0 ; index < bonds ; index ++ )
	{
		tripletList_G_uv2.push_back(T( index , index , 2.0 ));
		if( index != 0 )
		{

			tripletList_G_uv2.push_back(T( index - 1 , index , - scalars( index - 1 , 0 ) ));
			tripletList_G_uv2.push_back(T( index , index - 1 , - scalars( index - 1 , 0 ) ));
		}
	}

	G_uv2.setFromTriplets( tripletList_G_uv2.begin() , tripletList_G_uv2.end() );


	Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver;//
	solver.compute( G_uv2 );

	MatrixXd  G_uv2_inv = solver.solve( Identity );
	MatrixXd tmp = G_uv2_inv * n_i_mu2Transpose;
	MatrixXd tmp_2 = n_i_mu2 * tmp;





	projection_Matrix = this->Kronecker_Delta - tmp_2;
}



void Microtubule::getMatrixes_4(MatrixXd &projection_Matrix )
{
	MatrixXd scalars =  MatrixXd::Zero(this->numberOfPoints - 2 , 1 );
	for( unsigned int i = 0 ; i < this->numberOfPoints - 2 ; i ++ )
	{
		scalars( i , 0 ) = this->getTangent2( i ).dot(this->getTangent2( i + 1 ) ) / ( this->getTangent2( i ).norm() * this->getTangent2( i + 1 ).norm() );
	}

	unsigned int coordinates = 3.0 * this->numberOfPoints;
	unsigned int bonds = this->numberOfPoints - 1;


	SparseMatrix<double> n_i_mu2( coordinates , bonds );
        std::vector<T> tripletList_n_i_mu;
	for( unsigned int i = 0 ; i < this->numberOfPoints - 1 ; i ++ )
	{
		Vector3d tangent = this->getTangent2( i );
		tangent = tangent * ( 1.0 / tangent.norm() );
		for( unsigned int j = 0 ; j < 6 ; j ++ )
		{
			int tmp = 3 * i;
			if( j < 3)
			{
				tripletList_n_i_mu.push_back(T( tmp + j , j , - tangent( j ) ) );
			}
			else
			{
				tripletList_n_i_mu.push_back(T( tmp + j , i , - tangent( j - 3 ) ) );
			}
		}
	}
	n_i_mu2.setFromTriplets( tripletList_n_i_mu.begin() , tripletList_n_i_mu.end() );



	MatrixXd G_uv2 = MatrixXd::Zero( bonds , bonds );
	for( unsigned int index = 0 ; index < bonds ; index ++ )
	{
		G_uv2( index , index ) = 2;
		if( index != 0 )
		{
			G_uv2( index - 1 , index ) = - scalars( index - 1 , 0 );
			G_uv2( index , index - 1 ) = G_uv2( index - 1 , index );
		}
	}
	MatrixXd tmp = G_uv2.inverse() * n_i_mu2.transpose();
	tmp = n_i_mu2 * tmp;


	projection_Matrix = this->Kronecker_Delta - tmp;


}






void Microtubule::getMatrixes_Sparse( MatrixXd &inv_Matrix , MatrixXd &projection_Matrix )
{
	//Brownian Dynamics algorithm for bead-rod semiflexible chain with anisotropic friction
	//Montesi Morse Pasquali
	//Projection matrix is expressed in equation 16
	//projectionMatrix has to be Kronecker Delta			

	MatrixXd scalars =  MatrixXd::Zero(this->numberOfPoints - 2 , 1 );
	for( unsigned int i = 0 ; i < this->numberOfPoints - 2 ; i ++ )
	{
		scalars( i , 0 ) = this->getTangent2( i ).dot(this->getTangent2( i + 1 ) ) / ( this->getTangent2( i ).norm() * this->getTangent2( i + 1 ).norm() );
	}


	unsigned int coordinates = 3.0 * this->numberOfPoints;
	unsigned int bonds = this->numberOfPoints - 1;
	SparseMatrix<double> n_i_mu( coordinates , bonds );





	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList_n_i_mu;
	for( unsigned int i = 0 ; i < this->numberOfPoints - 1 ; i ++ )
	{
		Vector3d tangent = this->getTangent2( i );
		tangent = tangent * ( 1.0 / tangent.norm() );
		for( unsigned int j = 0 ; j < 6 ; j ++ )
		{
			int tmp = 3 * i;
			if( j < 3)
			{
				tripletList_n_i_mu.push_back(T( tmp + j , i , - tangent( j ) ));
			}
			else
			{   tripletList_n_i_mu.push_back(T( tmp + j , i , tangent( j - 3 ) ));
			}
		}
	}
	n_i_mu.setFromTriplets( tripletList_n_i_mu.begin() , tripletList_n_i_mu.end() );


	SparseMatrix<double> G_uv2( bonds , bonds );
	std::vector<T> tripletList_G_uv2;

	SparseMatrix<double> Identity( bonds , bonds );
	std::vector<T> tripletList_Identity;

	for( unsigned int index = 0 ; index < bonds ; index ++ )
	{
		tripletList_Identity.push_back(T( index , index , 1.0 ));
		tripletList_G_uv2.push_back(T( index , index , 2.0 ));
		if( index != 0 )
		{

			tripletList_G_uv2.push_back(T( index - 1 , index , - scalars( index - 1 , 0 ) ));
			tripletList_G_uv2.push_back(T( index , index - 1 , - scalars( index - 1 , 0 ) ));
		}
	}
	G_uv2.setFromTriplets( tripletList_G_uv2.begin() , tripletList_G_uv2.end() );
	Identity.setFromTriplets( tripletList_Identity.begin() , tripletList_Identity.end() );

	Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
	solver.compute( G_uv2 );

	SparseMatrix<double>  G_uv2_inv = solver.solve( Identity );
	inv_Matrix = G_uv2_inv;
	SparseMatrix<double> n_i_mu2Transpose = n_i_mu.transpose();

	MatrixXd tmp = inv_Matrix * n_i_mu2Transpose;
	tmp = n_i_mu * tmp;
	projection_Matrix = projection_Matrix - tmp;

}


void Microtubule::get_Sparse_Projection_Matrix_micro_7_second( MatrixXd &projection_Matrix )
{

	MatrixXd scalars =  MatrixXd::Zero(this->numberOfPoints - 2 , 1 );
	for( unsigned int i = 0 ; i < this->numberOfPoints - 2 ; i ++ )
	{
		scalars( i , 0 ) = this->getTangent2( i ).dot(this->getTangent2( i + 1 ) ) / ( this->getTangent2( i ).norm() * this->getTangent2( i + 1 ).norm() );
	}

	unsigned int coordinates = 3.0 * this->numberOfPoints;
	unsigned int bonds = 2 * this->numberOfPoints - 1;//
	SparseMatrix<double> n_i_mu( coordinates , bonds );

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList_n_i_mu;
	for( unsigned int i = 0 ; i < this->numberOfPoints - 1 ; i ++ )
	{
		Vector3d tangent = this->getTangent2( i );
		tangent = tangent * ( 1.0 / tangent.norm() );
		for( unsigned int j = 0 ; j < 6 ; j ++ )
		{
			int tmp = 3 * i;
			if( j < 3)
			{
				tripletList_n_i_mu.push_back(T( tmp + j , i , - tangent( j ) ));
			}
			else
			{   tripletList_n_i_mu.push_back(T( tmp + j , i , tangent( j - 3 ) ));
			}
		}

	}

	unsigned int bond_index = this->numberOfPoints - 1;
	for( unsigned int i = 0 ; i < this->numberOfPoints ; i ++ )
	{
        Vector3d position = this->getPoint( i );
        double multiple_B;
        if( position( 2 ) > 0 )
        {
            multiple_B = Cell_parametres::B_AXIS * Cell_parametres::B_AXIS;
        }
        {
            multiple_B = Cell_parametres::projection_radius * Cell_parametres::projection_radius;
        }
        double multiple_A = Cell_parametres::A_AXIS * Cell_parametres::A_AXIS;
		Vector3d tangent( 2.0 * multiple_B * position( 0 ) , 2.0 * multiple_B * position( 1 ) , 2.0 * multiple_A  * position( 2 ) );
        tangent = tangent / tangent.norm();
        for( unsigned int j = 0 ; j < 3 ; j ++ )
		{
            tripletList_n_i_mu.push_back(T( 3 * i + j , bond_index , tangent( j ) ));
        }
        bond_index = bond_index + 1;
    }




	n_i_mu.setFromTriplets( tripletList_n_i_mu.begin() , tripletList_n_i_mu.end() );
    n_i_mu.setFromTriplets( tripletList_n_i_mu.begin() , tripletList_n_i_mu.end() );
    n_i_mu.makeCompressed();

    SparseMatrix<double> n_i_mu2Transpose = n_i_mu.transpose();
    n_i_mu2Transpose.makeCompressed();


    MatrixXd G_uv2 = n_i_mu2Transpose * n_i_mu;
    MatrixXd G_uv2_inv = G_uv2.inverse();
    MatrixXd tmp = G_uv2_inv * n_i_mu2Transpose;
	tmp = n_i_mu * tmp;

    MatrixXd projection_Matrix_local = MatrixXd::Zero( 3.0 * this->getNumberOfPoints() , 3.0 * this->getNumberOfPoints() );
    for( unsigned int i = 0 ; i < 3.0 * this->getNumberOfPoints() ; i ++ )
    {
        projection_Matrix_local( i , i ) = 1.0;
    }
	projection_Matrix = projection_Matrix_local - tmp;
}










void Microtubule::set_bending_matrix( )
{

	MatrixXd MATRIX = MatrixXd::Zero( 3 * this->numberOfPoints , 3 * this->numberOfPoints );


	if( this->numberOfPoints == 3 )
	{
		MATRIX( 0 , 0 ) = -1.0;
		MATRIX( 0 , 3 ) = 2.0;
                MATRIX( 0 , 6 ) = -1.0;

		MATRIX( 1 , 1 ) = -1.0;
		MATRIX( 1 , 4 ) = 2.0;
                MATRIX( 1 , 7 ) = -1.0;


		MATRIX( 2 , 2 ) = -1.0;
		MATRIX( 2 , 5 ) = 2.0;
                MATRIX( 2 , 8 ) = -1.0;

		//second bead
		MATRIX( 3 , 0 ) = 2.0;
		MATRIX( 3 , 3 ) = -4.0;
                MATRIX( 3 , 6 ) = 2.0;

		MATRIX( 4 , 1 ) = 2.0;
		MATRIX( 4 , 4 ) = -4.0;
                MATRIX( 4 , 7 ) = 2.0;

		MATRIX( 5 , 2 ) = 2.0;
		MATRIX( 5 , 5 ) = -4.0;
                MATRIX( 5 , 8 ) = 2.0;

		//third bead
		MATRIX( 6 , 0 ) = -1.0;
		MATRIX( 6, 3 ) = 2.0;
                MATRIX( 6 , 6 ) = -1.0;

		MATRIX( 7 , 1 ) = -1.0;
		MATRIX( 7, 4 ) = 2.0;
                MATRIX( 7 , 7 ) = -1.0;

		MATRIX( 8 , 2 ) = -1.0;
		MATRIX( 8, 5 ) = 2.0;
                MATRIX( 8 , 8 ) = -1.0;

		//this->b_MATRIX = MATRIX;


	}
	else if( this->numberOfPoints == 4 )
	{
		//MatrixXd bending_matrix_4 = MatrixXd::Zero( 12 , 12 );

		//first bead
		MATRIX( 0 , 0 ) = -1.0;
		MATRIX( 0 , 3 ) = 2.0;
                MATRIX( 0 , 6 ) = -1.0;


		MATRIX( 1 , 1 ) = -1.0;
		MATRIX( 1 , 4 ) = 2.0;
                MATRIX( 1 , 7 ) = -1.0;


		MATRIX( 2 , 2 ) = -1.0;
		MATRIX( 2 , 5 ) = 2.0;
                MATRIX( 2 , 8 ) = -1.0;

		//second bead
		MATRIX( 3 , 0 ) = 2.0;
		MATRIX( 3 , 3 ) = -5.0;
                MATRIX( 3 , 6 ) = 4.0;
                MATRIX( 3 , 9 ) = -1.0;

		MATRIX( 4 , 1 ) = 2.0;
		MATRIX( 4 , 4 ) = -5.0;
                MATRIX( 4 , 7 ) = 4.0;
                MATRIX( 4 , 10 ) = -1.0;


		MATRIX( 5 , 2 ) = 2.0;
		MATRIX( 5 , 5 ) = -5.0;
                MATRIX( 5 , 8 ) = 4.0;
                MATRIX( 5 , 11 ) = -1.0;


		//third bead
		MATRIX( 6 , 0 ) = -1.0;
		MATRIX( 6 , 3 ) = 4.0;
                MATRIX( 6 , 6 ) = -5.0;
                MATRIX( 6 , 9 ) = 2.0;

		MATRIX( 7 , 1 ) = -1.0;
		MATRIX( 7 , 4 ) = 4.0;
                MATRIX( 7 , 7 ) = -5.0;
                MATRIX( 7 , 10 ) = 2.0;


		MATRIX( 8 , 2 ) = -1.0;
		MATRIX( 8 , 5 ) = 4.0;
                MATRIX( 8 , 8 ) = -5.0;
                MATRIX( 8 , 11 ) = 2.0;

		//fourth bead
		MATRIX( 9 , 3 ) = -1.0;
		MATRIX( 9 , 6 ) = 2.0;
                MATRIX( 9 , 9 ) = -1.0;

		MATRIX( 10 , 4 ) = -1.0;
		MATRIX( 10 , 7 ) = 2.0;
                MATRIX( 10 , 10 ) = -1.0;

		MATRIX( 11 , 5 ) = -1.0;
		MATRIX( 11 , 8 ) = 2.0;
                MATRIX( 11 , 11 ) = -1.0;
		//this->b_MATRIX = MATRIX;


	}
	else if( this->numberOfPoints >= 5 )
	{
		unsigned int N_beads = this->numberOfPoints;
//		MatrixXd MATRIX = MatrixXd::Zero( 3 * N_beads , 3 * N_beads );
		MATRIX( 0 , 0 ) = -1.0;
		MATRIX( 0 , 3 ) = 2.0;
                MATRIX( 0 , 6 ) = -1.0;


		MATRIX( 1 , 1 ) = -1.0;
		MATRIX( 1 , 4 ) = 2.0;
                MATRIX( 1 , 7 ) = -1.0;


		MATRIX( 2 , 2 ) = -1.0;
		MATRIX( 2 , 5 ) = 2.0;
                MATRIX( 2 , 8 ) = -1.0;

		//second bead
		MATRIX( 3 , 0 ) = 2.0;
		MATRIX( 3 , 3 ) = -5.0;
                MATRIX( 3 , 6 ) = 4.0;
                MATRIX( 3 , 9 ) = -1.0;

		MATRIX( 4 , 1 ) = 2.0;
		MATRIX( 4 , 4 ) = -5.0;
                MATRIX( 4 , 7 ) = 4.0;
                MATRIX( 4 , 10 ) = -1.0;


		MATRIX( 5 , 2 ) = 2.0;
		MATRIX( 5 , 5 ) = -5.0;
                MATRIX( 5 , 8 ) = 4.0;
                MATRIX( 5 , 11 ) = -1.0;


		//last bead
		MATRIX( 3 * ( N_beads - 1 ) , 3 * ( N_beads - 3 ) ) = -1.0;
		MATRIX( 3 * ( N_beads - 1 ) , 3 * ( N_beads - 2 )  ) = 2.0;
		MATRIX( 3 * ( N_beads - 1 ) , 3 * ( N_beads - 1 )  ) = -1.0;

		MATRIX( 3 * ( N_beads - 1 ) + 1 , 3 * ( N_beads - 3 ) + 1 ) = -1.0;
		MATRIX( 3 * ( N_beads - 1 ) + 1 , 3 * ( N_beads - 2 ) + 1 ) = 2.0;
		MATRIX( 3 * ( N_beads - 1 ) + 1 , 3 * ( N_beads - 1 ) + 1 ) = -1.0;

		MATRIX( 3 * ( N_beads - 1 ) + 2 , 3 * ( N_beads - 3 ) + 2 ) = -1.0;
		MATRIX( 3 * ( N_beads - 1 ) + 2 , 3 * ( N_beads - 2 ) + 2 ) = 2.0;
		MATRIX( 3 * ( N_beads - 1 ) + 2 , 3 * ( N_beads - 1 ) + 2 ) = -1.0;

		for( unsigned int bead_id = 0 ; bead_id < N_beads ; bead_id ++ )
		{

			if( ( bead_id > 1 ) && ( bead_id < N_beads - 2 ) )
			{

				MATRIX( 3 * ( bead_id ) , 3 * ( bead_id - 2 ) ) = - 1.0;
				MATRIX( 3 * ( bead_id ) , 3 * ( bead_id - 1 ) ) = 4.0;
				MATRIX( 3 * ( bead_id ) , 3 * ( bead_id ) ) = -6.0;
				MATRIX( 3 * ( bead_id ) , 3 * ( bead_id + 1 ) ) = 4.0;
				MATRIX( 3 * ( bead_id ) , 3 * ( bead_id + 2 ) ) = -1.0;

				MATRIX( 3 * ( bead_id ) + 1 , 3 * ( bead_id - 2 ) + 1 ) = - 1.0;
				MATRIX( 3 * ( bead_id ) + 1 , 3 * ( bead_id - 1 ) + 1 ) = 4.0;
				MATRIX( 3 * ( bead_id ) + 1 , 3 * ( bead_id ) + 1 ) = -6.0;
				MATRIX( 3 * ( bead_id ) + 1 , 3 * ( bead_id + 1 ) + 1 ) = 4.0;
				MATRIX( 3 * ( bead_id ) + 1 , 3 * ( bead_id + 2 ) + 1 ) = -1.0;

				MATRIX( 3 * ( bead_id ) + 2 , 3 * ( bead_id - 2 ) + 2 ) = - 1.0;
				MATRIX( 3 * ( bead_id ) + 2 , 3 * ( bead_id - 1 ) + 2 ) = 4.0;
				MATRIX( 3 * ( bead_id ) + 2 , 3 * ( bead_id ) + 2 ) = -6.0;
				MATRIX( 3 * ( bead_id ) + 2 , 3 * ( bead_id + 1 ) + 2 ) = 4.0;
				MATRIX( 3 * ( bead_id ) + 2 , 3 * ( bead_id + 2 ) + 2 ) = -1.0;


			}



		}


		MATRIX( 3 * ( N_beads - 2 ) , 3 * ( N_beads - 4 ) ) = -1.0;
		MATRIX( 3 * ( N_beads - 2 ) , 3 * ( N_beads - 3 )  ) = 4.0;
		MATRIX( 3 * ( N_beads - 2 ) , 3 * ( N_beads - 2 )  ) = -5.0;
		MATRIX( 3 * ( N_beads - 2 ) , 3 * ( N_beads - 1 )  ) = 2.0;

		MATRIX( 3 * ( N_beads - 2 ) + 1 , 3 * ( N_beads - 4 ) + 1 ) = -1.0;
		MATRIX( 3 * ( N_beads - 2 ) + 1 , 3 * ( N_beads - 3 ) + 1  ) = 4.0;
		MATRIX( 3 * ( N_beads - 2 ) + 1 , 3 * ( N_beads - 2 ) + 1  ) = -5.0;
		MATRIX( 3 * ( N_beads - 2 ) + 1 , 3 * ( N_beads - 1 ) + 1  ) = 2.0;

		MATRIX( 3 * ( N_beads - 2 ) + 2 , 3 * ( N_beads - 4 ) + 2 ) = -1.0;
		MATRIX( 3 * ( N_beads - 2 ) + 2 , 3 * ( N_beads - 3 ) + 2  ) = 4.0;
		MATRIX( 3 * ( N_beads - 2 ) + 2 , 3 * ( N_beads - 2 ) + 2  ) = -5.0;
		MATRIX( 3 * ( N_beads - 2 ) + 2 , 3 * ( N_beads - 1 ) + 2  ) = 2.0;

	}



	SparseMatrix<double> b_SPARSE( 3 * this->numberOfPoints , 3 * this->numberOfPoints );
        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList_n_i_mu;
	for( unsigned int row = 0 ; row < 3 * this->numberOfPoints ; row ++ )
	{
		for( unsigned int colum = 0 ; colum < 3 * this->numberOfPoints ; colum ++ )
		{
			if( MATRIX( row , colum ) != 0 )
			{
				tripletList_n_i_mu.push_back(T( row , colum , MATRIX( row , colum ) ));
			}
    		}

    	}

        b_SPARSE.setFromTriplets( tripletList_n_i_mu.begin() , tripletList_n_i_mu.end() );
 	this->b_MATRIX = b_SPARSE;


	MatrixXd projectionMatrix = MatrixXd::Zero( 3 * this->numberOfPoints , 3 * this->numberOfPoints );
	for( unsigned int index = 0 ; index < 3 * this->numberOfPoints ; index ++ )
	{
		projectionMatrix( index , index ) = 1.0;
	}



	//


	SparseMatrix<double> kronecker( 3 * this->numberOfPoints , 3 * this->numberOfPoints );
        std::vector<T> triplet_kronecker;
	for( unsigned int row = 0 ; row < 3 * this->numberOfPoints ; row ++ )
	{
		for( unsigned int colum = 0 ; colum < 3 * this->numberOfPoints ; colum ++ )
		{
			if( MATRIX( row , colum ) != 0 )
			{
				triplet_kronecker.push_back(T( row , colum , projectionMatrix( row , colum ) ));
			}
    		}

    	}
	this->Kronecker_Delta = kronecker;




}



void Microtubule::getBendingForces_2( MatrixXd &bendingForce )
{
		if( ( bendingForce.rows() !=  3 * this->numberOfPoints ) || ( ( bendingForce.cols() != 1 ) ) )
		{
			cout<<"( bendingForce.getNN() != 3 * this->numberOfPoints ) || ( ( bendingForce.getNN() != 1 ) ) in Microtubule::getBendingForces( MatDoub bendingForce )"<<endl;
			cout<<this->numberOfPoints <<endl;
			cout<<this->get_dynein_index()<<endl;

			throw("");
		}

		if( this->dydein_index == 20 )
		{
			cout<<"Microtubule ERROR_ID = "<<6466156156135<<endl;
			cout<<" this->dydein_index == 20 "<<endl;
			cout<<"void Microtubule::getBendingForces_2( MatrixXd &bendingForce )"<<endl;
			throw("");
		}

		bendingForce = this->b_MATRIX * this->coordinates * this->kappa / ( this->restDistancePoints * this->restDistancePoints);

}








void Microtubule::getBendingForces( MatrixXd &bendingForce )
{


		if( ( bendingForce.rows() !=  3 * this->numberOfPoints ) || ( ( bendingForce.cols() != 1 ) ) )
		{
			cout<<"( bendingForce.getNN() != 3 * this->numberOfPoints ) || ( ( bendingForce.getNN() != 1 ) ) in Microtubule::getBendingForces( MatDoub bendingForce )"<<endl;
			cout<<this->numberOfPoints <<endl;
			cout<<this->get_dynein_index()<<endl;

			throw("");
		}

		if( ( this->getNumberOfPoints() <= 2 )  )
        	{
	    		//bendingForce = bendingForce * this->kappa / 10.0;
            		return;
        	}

		if( ( this->getNumberOfPoints() == 3 )  )
        	{

			for( unsigned int index = 0 ; index < this->getNumberOfPoints() ; index ++ )
			{

				Vector3d first = this->getTangent2( 0 );
				Vector3d second = this->getTangent2( 1 );

				double coefficient = second.norm() / ( this->restDistancePoints );
				double treshold = 0.1;
				if( coefficient < treshold )
				{

				}
				else
				{
					coefficient = treshold;
				}


				first = first * ( 1.0 / first.norm() );
				second = second * ( 1.0 / second.norm() );



				double a_1 = ( - second( 0 ) + first( 0 ) * first.dot( second ) ) / this->restDistancePoints * coefficient;
				double a_2 = ( - second( 1 ) + first( 1 ) * first.dot( second ) ) / this->restDistancePoints * coefficient;
				double a_3 = ( - second( 2 ) + first( 2 ) * first.dot( second ) ) / this->restDistancePoints * coefficient;

				double b_1 = ( first( 0 ) - second( 0 ) * first.dot( second ) ) / this->restDistancePoints * coefficient;
				double b_2 = ( first( 1 ) - second( 1 ) * first.dot( second ) ) / this->restDistancePoints * coefficient;
				double b_3 = ( first( 2 ) - second( 2 ) * first.dot( second ) ) / this->restDistancePoints * coefficient;


				if( index == 0 )										
				{
					bendingForce( 3 * index +  0 , 0 ) = a_1;
					bendingForce( 3 * index +  1 , 0 ) = a_2;
					bendingForce( 3 * index +  2 , 0 ) = a_3;
				}

				else if( index == 2 )		
				{
					bendingForce( 3 * index + 0 , 0 ) = b_1;
					bendingForce( 3 * index + 1 , 0 ) = b_2;
					bendingForce( 3 * index + 2 , 0 ) = b_3;

				}

				else if(  index == 1  )
				{
					bendingForce( 3 * index +  0 , 0 ) = ( -1.0 ) * ( a_1 + b_1 );
					bendingForce( 3 * index +  1 , 0 ) = ( -1.0 ) * ( a_2 + b_2 );
					bendingForce( 3 * index +  2 , 0 ) = ( -1.0 ) * ( a_3 + b_3 );

				}
			}
        	}

		else
		{

		    for( unsigned int index = 0 ; index < this->numberOfPoints ; index ++ )
		    {

			if( index == 0 )										//this is border case of bead of minimum number
			{
				Vector3d first( 0.0 , 0.0 , 0.0 );
				Vector3d second( 0.0 , 0.0 , 0.0 );

				for( unsigned int j = 0 ; j < 3 ; j ++ )
				{
					first( j ) = this->coordinates( 3 * ( index + 1 ) + j , 0 ) - this->coordinates( j , 0 );
					second( j ) =  this->coordinates( 3 * ( index + 2 ) + j , 0 ) - this->coordinates( 3 * ( index + 1 ) + j , 0 );
				}

				double magThisVec = first.norm();
				first = first * ( 1.0 / magThisVec );

				double magNextVec = second.norm();
				second = second * ( 1.0 / magNextVec );

				bendingForce( 3 * index +  0 , 0 ) = ( - second( 0 ) + first( 0 ) * first.dot( second ) ) / magThisVec;
				bendingForce( 3 * index +  1 , 0 ) = ( - second( 1 ) + first( 1 ) * first.dot( second ) ) / magThisVec;
				bendingForce( 3 * index +  2 , 0 ) = ( - second( 2 ) + first( 2 ) * first.dot( second ) ) / magThisVec;
			}

			else if( index == this->numberOfPoints - 1 )		
			{

				Vector3d first( 0.0 , 0.0 , 0.0 );
				Vector3d second( 0.0 , 0.0 , 0.0 );
				for( unsigned int j = 0 ; j < 3 ; j ++ )
				{

					first( j ) = this->coordinates( 3 * ( index - 1 ) + j , 0 ) - this->coordinates( 3 * ( index - 2 ) + j , 0 );
					second( j ) = this->coordinates( 3 * index + j , 0 ) - this->coordinates( 3 * ( index - 1 ) + j , 0 );
				}

				double magThisVec = first.norm();
				first = first * ( 1.0 / magThisVec );

				double magNextVec = second.norm();
				second = second * ( 1.0 / magNextVec );		

				bendingForce( 3 * index + 0 , 0 ) = ( first( 0 ) - second( 0 ) * first.dot( second ) ) / magNextVec;
				bendingForce( 3 * index + 1 , 0 ) = ( first( 1 ) - second( 1 ) * first.dot( second ) ) / magNextVec;
				bendingForce( 3 * index + 2 , 0 ) = ( first( 2 ) - second( 2 ) * first.dot( second ) ) / magNextVec;

			}
			else
			{
				if( this->numberOfPoints > 3 )
				{
					Vector3d tangent_j = this->getTangent2( index );
                                        Vector3d tangent_j_M_1 = this->getTangent2( index - 1 );

					double tangent_jNorm = tangent_j.norm();
					tangent_j = tangent_j * ( 1.0 / tangent_jNorm );

					double tangent_j_M_1Norm = tangent_j_M_1.norm();
					tangent_j_M_1 = tangent_j_M_1 * ( 1.0 / tangent_j_M_1Norm );


					bendingForce( 3 * index + 0 , 0 ) = ( tangent_j( 0 ) - tangent_j_M_1( 0 ) * tangent_j.dot( tangent_j_M_1 ) ) / tangent_j_M_1Norm;
					bendingForce( 3 * index + 0 , 0 ) = bendingForce( 3 * index + 0 , 0 ) + ( - tangent_j_M_1( 0 ) + tangent_j( 0 ) * tangent_j.dot( tangent_j_M_1 ) ) / tangent_jNorm;

					bendingForce( 3 * index + 1 , 0 ) = ( tangent_j( 1 ) - tangent_j_M_1( 1 ) * tangent_j.dot( tangent_j_M_1 ) ) / tangent_j_M_1Norm;
					bendingForce( 3 * index + 1 , 0 ) = bendingForce( 3 * index + 1 , 0 ) + ( - tangent_j_M_1( 1 ) + tangent_j( 1 ) * tangent_j.dot( tangent_j_M_1 ) ) / tangent_jNorm;

					bendingForce( 3 * index + 2 , 0 ) = ( tangent_j( 2 ) - tangent_j_M_1( 2 ) * tangent_j.dot( tangent_j_M_1 ) ) / tangent_j_M_1Norm;
					bendingForce( 3 * index + 2 , 0 ) = bendingForce( 3 * index + 2 , 0 ) + ( - tangent_j_M_1( 2 ) + tangent_j( 2 ) * tangent_j.dot( tangent_j_M_1 ) ) / tangent_jNorm;


					if ( index == this->numberOfPoints - 2 )
					{

					}

					else
					{

                  				Vector3d first = this->getTangent2( index );
                                                Vector3d second = this->getTangent2( index + 1 );

						double magThisVec = first.norm();
						first = first * ( 1.0 / magThisVec );

						double magNextVec = second.norm();
						second = second * ( 1.0 / magNextVec );		//ATTENTION - here, it will be divided by magNextVec - the last tangent


						bendingForce( 3 * index + 0 , 0 ) = bendingForce( 3 * index + 0 , 0 ) + ( - second( 0 ) + first( 0 ) * tangent_j.dot( second ) ) / magThisVec;
						bendingForce( 3 * index + 1 , 0 ) = bendingForce( 3 * index + 1 , 0 ) + ( - second( 1 ) + first( 1 ) * tangent_j.dot( second ) ) / magThisVec;
						bendingForce( 3 * index + 2 , 0 ) = bendingForce( 3 * index + 2 , 0 ) + ( - second( 2 ) + first( 2 ) * tangent_j.dot( second ) ) / magThisVec;
					}


					if ( index == 1 )
					{

					}
					else
					{
                  				Vector3d first = this->getTangent2( index - 1 );
                                                Vector3d second = this->getTangent2( index - 2 );

						double magThisVec = first.norm();
						first = first * ( 1.0 / magThisVec );

						double magNextVec = second.norm();
						second = second * ( 1.0 / magNextVec );

						bendingForce( 3 * index + 0 , 0 ) = bendingForce( 3 * index + 0 , 0 ) + ( second( 0 ) - first( 0 ) * first.dot( second ) ) / magThisVec;
						bendingForce( 3 * index + 1 , 0 ) = bendingForce( 3 * index + 1 , 0 ) + ( second( 1 ) - first( 1 ) * first.dot( second ) ) / magThisVec;
						bendingForce( 3 * index + 2 , 0 ) = bendingForce( 3 * index + 2 , 0 ) + ( second( 2 ) - first( 2 ) * first.dot( second ) ) / magThisVec;

					}
				}
			}

		    }

		}


		bendingForce = bendingForce * this->kappa;
}









void Microtubule::getRandomForces(  double timeStep , MatrixXd &randomForces )
{

	std::normal_distribution<> distribution{0,1};
        unsigned int number_of_generator = omp_get_thread_num();


	if( ( randomForces.rows() !=  3 * this->numberOfPoints ) || ( ( randomForces.cols() != 1 ) ) )
	{
		cout<<"( bendingForce.getNN() != 3 * this->numberOfPoints ) || ( ( bendingForce.getNN() != 1 ) ) in Microtubule::getBendingForces( MatDoub bendingForce )"<<endl;
		throw("");
	}

	unsigned int coordinates = 3 * this->numberOfPoints;
	double root = 1.0 * sqrt( 2.0 * sim_of_Cell::boltzmann_constant * sim_of_Cell::Temperature * this->effective_friction / timeStep );

	if( sim_of_Cell::random_force_switch == true )
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







void Microtubule::metricForceWLC( MatrixXd &metricForces , MatrixXd &G_uv_Inverse )
{

	if( ( metricForces.rows() !=  3 * this->numberOfPoints ) || ( ( metricForces.cols() != 1 ) ) )
	{
		cout<<"( bendingForce.getNN() != 3 * this->numberOfPoints ) || ( ( bendingForce.getNN() != 1 ) ) in metricForceWLC( MatrixXd &metricForces , MatrixXd &G_uv_Inverse )"<<endl;
		throw("");
	}

	if( ( G_uv_Inverse.rows() !=  this->numberOfPoints - 1 ) || ( ( G_uv_Inverse.cols() != G_uv_Inverse.rows() ) ) )
	{
		cout<<" G_uv_Inverse.rows() = "<< G_uv_Inverse.rows()<<endl;
		cout<<"( ( G_uv_Inverse.rows() !=  this->numberOfPoints - 1 ) || ( ( G_uv_Inverse.cols() != G_uv_Inverse.rows() ) ) in Microtubule::metricForceWLC( MatrixXd &metricForces , MatrixXd &G_uv_Inverse )"<<endl;
		throw("");
	}

	if( this->numberOfPoints < 3)
	{
		return;
	}

	MatrixXd scalars =  MatrixXd::Zero( this->numberOfPoints - 2 , 1 );
	for( unsigned int i = 0 ; i < this->numberOfPoints - 2 ; i ++ )
	{
		scalars( i , 0 ) = this->getTangent2( i ).dot(this->getTangent2( i + 1 ) ) / ( this->getTangent2( i ).norm() * this->getTangent2( i + 1 ).norm() );
	}


	for( unsigned int i = 0 ; i < this->numberOfPoints ; i ++ )
	{
		if( i == 0 )
		{
			Vector3d force = der_scalar_index( 0 , 0 );
			force = force * ( ( - 1.0 ) * 2.0 * sim_of_Cell::Temperature * sim_of_Cell::boltzmann_constant );
			double tmp =  G_uv_Inverse( 0 , 1 );
			force = force * ( tmp );

			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				metricForces( j , 0 ) = force( j );
			}
		}

		else if( i == this->numberOfPoints - 1 )
		{
			Vector3d force = der_scalar_index( this->numberOfPoints - 3 , this->numberOfPoints - 1 );
			double tmp =  G_uv_Inverse( i - 2 , i - 1 );
			force = force * ( tmp );
			force = force * ( ( - 1.0 ) * 2.0 * sim_of_Cell::Temperature * sim_of_Cell::boltzmann_constant );
			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				metricForces( 3 * i + j , 0 ) = force( j );
			}
		}

		else
		{
			if( this->numberOfPoints == 3 )
			{
				Vector3d force = der_scalar_index( 0 , 1 );
				double tmp =  G_uv_Inverse( 0 , 1 );
				force = force * ( tmp );
				force = force * ( ( - 1.0 ) * 2.0 * sim_of_Cell::Temperature * sim_of_Cell::boltzmann_constant );
				for( unsigned int j = 0 ; j < 3 ; j ++ )
				{
					metricForces( 3 * i + j , 0 ) = force( j );
				}
			}

			else if( i == 1 )
			{
				Vector3d force1 = der_scalar_index( 0 , 1 );
				double tmp1 = G_uv_Inverse( 0 , 1 );
				force1 = force1 * ( tmp1 );
				Vector3d force2 = der_scalar_index( 1 , 1 );
				double tmp2 = G_uv_Inverse( 1 , 2 );
				force2 = force2 * ( tmp2 );
				Vector3d force = ( force1 + force2 ) * ( ( -2.0 ) * sim_of_Cell::boltzmann_constant * sim_of_Cell::Temperature );

				for( unsigned int j = 0 ; j < 3 ; j ++ )
				{
					metricForces( 3 * ( i ) + j , 0 ) = force( j );
				}
			}

			else if( i == this->numberOfPoints - 2 )
			{
				Vector3d force1 = this->der_scalar_index( i - 2 , i );
				double tmp1 = G_uv_Inverse(  i - 2 , i - 1 );
				force1 = force1 * ( tmp1 );

				Vector3d force2 = this->der_scalar_index( i - 1 , i );
				double tmp2 = G_uv_Inverse( i - 1 , i );
				force2 = force2 * ( tmp2 );

				Vector3d force = ( force1 + force2 ) * ( ( -2.0 ) * sim_of_Cell::boltzmann_constant * sim_of_Cell::Temperature );

				for( unsigned int j = 0 ; j < 3 ; j ++ )
				{
					metricForces( 3 * ( i ) + j , 0 ) = force( j );
				}
			}
			else
			{
				Vector3d force1 = this->der_scalar_index( i - 2 , i );
				double tmp1 = G_uv_Inverse( i - 2 , i - 1 );
				force1 = force1 * ( tmp1 );

				Vector3d force2 = this->der_scalar_index( i - 1 , i );
				double tmp2 = G_uv_Inverse( i - 1 , i );
				force2 = force2 * ( tmp2 );

				Vector3d force3 = this->der_scalar_index( i , i );
				double tmp3 = G_uv_Inverse( i , i + 1 );
				force3 = force3 * ( tmp3 );

				Vector3d force = ( force1 + force2 + force3 ) * (  ( -2.0 ) * sim_of_Cell::boltzmann_constant * sim_of_Cell::Temperature );

				for( unsigned int j = 0 ; j < 3 ; j ++ )
				{
					metricForces( 3 * ( i ) + j , 0 ) = force( j );
				}
			}
		}
	}
}

void Microtubule::metricForceConsLenght( MatrixXd &metrFor_C_L , MatrixXd &G_uv_Inverse )
{	

	if( ( metrFor_C_L.rows() !=  3 * this->numberOfPoints ) || ( ( metrFor_C_L.cols() != 1 ) ) )
	{
		cout<<"( bendingForce.getNN() != 3 * this->numberOfPoints ) || ( ( bendingForce.getNN() != 1 ) ) in metricForceWLC( MatrixXd &metricForces , MatrixXd &G_uv_Inverse )"<<endl;
		throw("");
	}

	if( ( G_uv_Inverse.rows() !=  this->numberOfPoints - 1 ) || ( ( G_uv_Inverse.cols() != G_uv_Inverse.rows() ) ) )
	{
		cout<<" G_uv_Inverse.rows() = "<< G_uv_Inverse.rows()<<endl;
		cout<<"( ( G_uv_Inverse.rows() !=  this->numberOfPoints - 1 ) || ( ( G_uv_Inverse.cols() != G_uv_Inverse.rows() ) ) in Microtubule::metricForceWLC( MatrixXd &metricForces , MatrixXd &G_uv_Inverse )"<<endl;
		throw("");
	}


	if( this->numberOfPoints < 3)
	{
		return;
	}


	MatrixXd scalars =  MatrixXd::Zero( this->numberOfPoints - 2 , 1 );
	for( unsigned int i = 0 ; i < this->numberOfPoints - 2 ; i ++ )
	{
		scalars( i , 0 ) = this->getTangent2( i ).dot(this->getTangent2( i + 1 ) ) / ( this->getTangent2( i ).norm() * this->getTangent2( i + 1 ).norm() );
	}


	for( unsigned int i = 0 ; i < this->numberOfPoints ; i ++ )
	{
		//this for loop calculates metric forces for every point and writes into 3 * N coordinates
		if( i == 0 )
		{
			Vector3d force = der_scalar_index( 0 , 0 );
			force = force * ( sim_of_Cell::Temperature * sim_of_Cell::boltzmann_constant );
			double tmp =  G_uv_Inverse( 0 , 1 );
			force = force * ( tmp );

			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				metrFor_C_L( j , 0 ) = force( j );
			}
		}

		else if( i == this->numberOfPoints - 1 )
		{
			Vector3d force = der_scalar_index( this->numberOfPoints - 3 , this->numberOfPoints - 1 );
			double tmp =  G_uv_Inverse( i - 2 , i - 1 );
			force = force * ( tmp );
			force = force * ( sim_of_Cell::Temperature * sim_of_Cell::boltzmann_constant );
			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				metrFor_C_L( 3 * i + j , 0 ) = force( j );
			}
		}

		else
		{
			if( this->numberOfPoints == 3 )
			{
				Vector3d force = der_scalar_index( 0 , 1 );
				double tmp =  G_uv_Inverse( 0 , 1 );
				force = force * ( tmp );
				force = force * ( sim_of_Cell::Temperature * sim_of_Cell::boltzmann_constant );
				for( unsigned int j = 0 ; j < 3 ; j ++ )
				{
					metrFor_C_L( 3 * i + j , 0 ) = force( j );
				}
			}

			else if( i == 1 )
			{
				Vector3d force1 = der_scalar_index( 0 , 1 );
				double tmp1 = G_uv_Inverse( 0 , 1 );
				force1 = force1 * ( tmp1 );
				Vector3d force2 = der_scalar_index( 1 , 1 );
				double tmp2 = G_uv_Inverse( 1 , 2 );
				force2 = force2 * ( tmp2 );
				Vector3d force = ( force1 + force2 ) * ( sim_of_Cell::boltzmann_constant * sim_of_Cell::Temperature );

				for( unsigned int j = 0 ; j < 3 ; j ++ )
				{
					metrFor_C_L( 3 * ( i ) + j , 0 ) = force( j );
				}
			}

			else if( i == this->numberOfPoints - 2 )
			{
				Vector3d force1 = this->der_scalar_index( i - 2 , i );
				double tmp1 = G_uv_Inverse(  i - 2 , i - 1 );
				force1 = force1 * ( tmp1 );

				Vector3d force2 = this->der_scalar_index( i - 1 , i );
				double tmp2 = G_uv_Inverse( i - 1 , i );
				force2 = force2 * ( tmp2 );

				Vector3d force = ( force1 + force2 ) * ( sim_of_Cell::boltzmann_constant * sim_of_Cell::Temperature );

				for( unsigned int j = 0 ; j < 3 ; j ++ )
				{
					metrFor_C_L( 3 * ( i ) + j , 0 ) = force( j );
				}
			}
			else
			{
				Vector3d force1 = this->der_scalar_index( i - 2 , i );
				double tmp1 = G_uv_Inverse( i - 2 , i - 1 );
				force1 = force1 * ( tmp1 );

				Vector3d force2 = this->der_scalar_index( i - 1 , i );
				double tmp2 = G_uv_Inverse( i - 1 , i );
				force2 = force2 * ( tmp2 );

				Vector3d force3 = this->der_scalar_index( i , i );
				double tmp3 = G_uv_Inverse( i , i + 1 );
				force3 = force3 * ( tmp3 );

				Vector3d force = ( force1 + force2 + force3 ) * ( sim_of_Cell::boltzmann_constant * sim_of_Cell::Temperature );

				for( unsigned int j = 0 ; j < 3 ; j ++ )
				{
					metrFor_C_L( 3 * ( i ) + j , 0 ) = force( j );
				}
			}
		}
	}

}







Vector3d Microtubule::dynein_force_capture_shrinkage( )
{
	Vector3d position_last_bead = this->getPoint( this->getNumberOfPoints() - 1 );
	Vector3d orientace = this->get_IS_position_catching() - position_last_bead;
	orientace = orientace / orientace.norm();
	Vector3d force = orientace * IS_Capture_shrinkage_param::dynein_Force_capture_shrinkage;
	return force;
}










Vector3d Microtubule::dynein_force_cortical_sliding( )
{

	Vector3d position_last_bead = this->getPoint( this->getNumberOfPoints() - 1 );
	Vector3d orientace = this->get_IS_position_catching() - position_last_bead;

	orientace = orientace / orientace.norm();
	Vector3d force = orientace * sim_of_Cell::dynein_Force_cortical_sliding;
	return force;
}





void Microtubule::oneStepMidStepAlgorithm()
{

	//MICROTUBULES
	MatrixXd random_Forces;
	MatrixXd original_Coordinates;
	MatrixXd external_Forces;


	if( this->get_dynein_index() <= 1 )
	{
		original_Coordinates = this->getCoordinates();
		random_Forces = MatrixXd::Zero( 3 * this->getNumberOfPoints() , 1 );
		this->getRandomForces( sim_of_Cell::time_Step , random_Forces );

		external_Forces = MatrixXd::Zero( 3 * this->getNumberOfPoints() , 1 );

		//The dynein is added here
		//dynein forces - they apply only when microtubule is cought by dynein motors -
		if( this->get_dynein_index() == 1 )
		{
			Vector3d dyneinForceVector = this->dynein_force_capture_shrinkage( );

			for( unsigned int i = 0 ; i < 3 ; i ++ )
			{
				external_Forces( 3 * ( this->getNumberOfPoints() - 1 ) + i , 0 ) = dyneinForceVector( i );
			}
		}
		this->oneStepMidStepAlgorithm_1_half( random_Forces , external_Forces );
	}


	if( this->get_dynein_index() <= 1 )
	{
		external_Forces = MatrixXd::Zero( 3 * this->getNumberOfPoints() , 1 );

		this->oneStepMidStepAlgorithm_2_half( original_Coordinates , random_Forces , external_Forces );
	}

}





void Microtubule::oneStepMidStepAlgorithm_1_half( MatrixXd& randomForces , MatrixXd& extrenal_Forces  )
{

	if( ( randomForces.rows() != this->numberOfPoints * 3 ) || ( randomForces.cols() != 1 )  )
	{
		cout<<"( randomForces.rows() != this->numberOfPoints * 3 ) || ( randomForces.cols() != 1 )  in oneStepMidStepAlgorithm_1_half "<<endl;
		cout<<" randomForces.rows() = "<<randomForces.rows()<<endl;
		cout<<"randomForces.cols() = "<<randomForces.cols()<<endl;
		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
		throw("");
	}

	unsigned int rows_tmp = extrenal_Forces.rows();
    unsigned int arg_tmp = this->numberOfPoints * 3;

    if( ( rows_tmp != arg_tmp )  )
	{
        cout<<"ERROR_ID Microtubule6461651685152"<<endl;
		cout<<"( extrenal_Forces.rows() != this->numberOfPoints * 3 ) || ( extrenal_Forces.cols() != 1 )  in oneStepMidStepAlgorithm_1_half "<<endl;
		cout<<" extrenal_Forces.rows() = "<<extrenal_Forces.rows()<<endl;
		cout<<"extrenal_Forces.cols() = "<<randomForces.cols()<<endl;
        cout<<"this->numberOfPoints = "<<this->numberOfPoints<<endl;
		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
        cout<<"this->microtubule_id = "<<this->microtubule_id<<endl;
		throw("");
	}

	if( this->numberOfPoints == 0 )
	{
		cout<<"this->numberOfPoints == 0 in Microtubule::oneStepMidStepAlgorithm_1_half( MatrixXd& randomForces , MatrixXd& extrenal_Forces  )"<<endl;
		throw("");
	}
	//Coordinates before step
	MatrixXd original = this->getCoordinates();

	if( this->numberOfPoints == 1 )
	{

		MatrixXd V_0 = (  extrenal_Forces + randomForces ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original + V_0 * sim_of_Cell::time_Step_half;
		this->setCoordinates( R_half );
	}
	else
	{
		MatrixXd projectionMatrix = MatrixXd::Zero( 3 * this->numberOfPoints , 3 * this->numberOfPoints );
		for( unsigned int index = 0 ; index < 3 * this->numberOfPoints ; index ++ )
		{
			projectionMatrix( index , index ) = 1.0;
		}
		MatrixXd G_uvINV = MatrixXd::Zero( this->numberOfPoints - 1 , this->numberOfPoints - 1 );

        	this->getMatrixes( G_uvINV , projectionMatrix );

     		double timeHalfStep =  sim_of_Cell::time_Step_half; /// 2.0

		MatrixXd bendingForce = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
		if( this->get_dynein_index() == 20 )
		{
			this->getBendingForces( bendingForce );
		}
		else
		{
			this->getBendingForces_2( bendingForce );
		}

		MatrixXd forces = bendingForce + extrenal_Forces + randomForces;


		MatrixXd forces_Projection = projectionMatrix * forces;


		MatrixXd V_0 = ( forces_Projection ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original + V_0 * timeHalfStep;
		this->setCoordinates( R_half );

	}

}



void Microtubule::oneStepMidStepAlgorithm_1_half_producing_random( MatrixXd& randomForces , MatrixXd& extrenal_Forces )
{
	if( ( randomForces.rows() != this->numberOfPoints * 3 ) || ( randomForces.cols() != 1 )  )
	{
		cout<<"( randomForces.rows() != this->numberOfPoints * 3 ) || ( randomForces.cols() != 1 )  in oneStepMidStepAlgorithm_1_half "<<endl;
		cout<<" randomForces.rows() = "<<randomForces.rows()<<endl;
		cout<<"randomForces.cols() = "<<randomForces.cols()<<endl;
		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
		throw("");
	}

	unsigned int rows_tmp = extrenal_Forces.rows();
    unsigned int arg_tmp = this->numberOfPoints * 3;

    if( ( rows_tmp != arg_tmp )  )
	{
        cout<<"ERROR_ID Microtubule6461651685152"<<endl;
		cout<<"( extrenal_Forces.rows() != this->numberOfPoints * 3 ) || ( extrenal_Forces.cols() != 1 )  in oneStepMidStepAlgorithm_1_half "<<endl;
		cout<<" extrenal_Forces.rows() = "<<extrenal_Forces.rows()<<endl;
		cout<<"extrenal_Forces.cols() = "<<randomForces.cols()<<endl;
        cout<<"this->numberOfPoints = "<<this->numberOfPoints<<endl;
		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
        cout<<"this->microtubule_id = "<<this->microtubule_id<<endl;
		throw("");
	}

	if( this->numberOfPoints == 0 )
	{
		cout<<"this->numberOfPoints == 0 in Microtubule::oneStepMidStepAlgorithm_1_half( MatrixXd& randomForces , MatrixXd& extrenal_Forces  )"<<endl;
		throw("");
	}
	MatrixXd original = this->getCoordinates();

	if( this->numberOfPoints == 1 )
	{
		//Physically speaking, no changes in inertia caused by configuration
		MatrixXd V_0 = (  extrenal_Forces + randomForces ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original + V_0 * sim_of_Cell::time_Step_half;
        	this->setCoordinates( R_half );
	}
	else
	{

		MatrixXd projectionMatrix = MatrixXd::Zero( 3 * this->numberOfPoints , 3 * this->numberOfPoints );
		for( unsigned int index = 0 ; index < 3 * this->numberOfPoints ; index ++ )
		{
			projectionMatrix( index , index ) = 1.0;
		}
		MatrixXd G_uvINV = MatrixXd::Zero( this->numberOfPoints - 1 , this->numberOfPoints - 1 );
        	this->getMatrixes( G_uvINV , projectionMatrix );


     		double timeHalfStep =  sim_of_Cell::time_Step_half; /// 2.0

		MatrixXd bendingForce = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
    this->getBendingForces( bendingForce );


		MatrixXd forces = bendingForce + extrenal_Forces;
		MatrixXd forces_Projection = projectionMatrix * forces  + randomForces;


		MatrixXd V_0 = ( forces_Projection ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original + V_0 * timeHalfStep;


		this->setCoordinates( R_half );

	}


}

















void Microtubule::oneStepMidStepAlgorithm_2_half( MatrixXd& original_Coordinates , MatrixXd& randomForces , MatrixXd& extrenal_Forces  )
{

	if( this->numberOfPoints == 1 )
	{

		MatrixXd V_0 = (  extrenal_Forces + randomForces ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original_Coordinates + V_0 * sim_of_Cell::time_Step;
		this->setCoordinates( R_half );
	}
	else
	{
		//Kronecker delta
		MatrixXd projectionMatrix = MatrixXd::Zero( 3 * this->numberOfPoints , 3 * this->numberOfPoints );
		for( unsigned int index = 0 ; index < 3 * this->numberOfPoints ; index ++ )
		{
			projectionMatrix( index , index ) = 1.0;
		}

		// Inversin matrix
		MatrixXd G_uvINV = MatrixXd::Zero( this->numberOfPoints - 1 , this->numberOfPoints - 1 );
        	this->getMatrixes( G_uvINV , projectionMatrix );



		// ------------------------------------------------------------------------------------------ Inter Microtubule forces
		//Bending Forces
		MatrixXd bendingForce = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );

    this->getBendingForces( bendingForce );


		//

		// ------------------------------------------------------------------------------------------ Projections
		//Projection
		MatrixXd forces = bendingForce + extrenal_Forces;
		MatrixXd forces_Projection = projectionMatrix * forces + randomForces;

		MatrixXd V_0 = ( forces_Projection ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original_Coordinates + V_0 * sim_of_Cell::time_Step;

        	this->setCoordinates( R_half );

	}
}



void Microtubule::project_force( MatrixXd& original_force , MatrixXd& projected_forces )
{
	if( ( original_force.rows() != this->numberOfPoints * 3 ) || ( original_force.cols() != 1 )  )
	{
		cout<<"( original_force.rows() != this->original_force * 3 ) || ( randomForces.cols() != 1 )  in oneStepMidStepAlgorithm_1_half "<<endl;
		cout<<" original_force.rows() = "<<original_force.rows()<<endl;
		cout<<"original_force.cols() = "<<original_force.cols()<<endl;
		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
		cout<<"Microtubule::project_force( MatrixXd& original_force , MatrixXd& projected_forces )"<<endl;
		throw("");
	}
	if( ( projected_forces.rows() != this->numberOfPoints * 3 ) || ( projected_forces.cols() != 1 )  )
	{
		cout<<"( projected_forces.rows() != this->numberOfPoints * 3 ) || ( projected_forces.cols() != 1 )  in oneStepMidStepAlgorithm_1_half "<<endl;
		cout<<" projected_forces.rows() = "<<projected_forces.rows()<<endl;
		cout<<"projected_forces.cols() = "<<projected_forces.cols()<<endl;
		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
		cout<<"Microtubule::project_force( MatrixXd& original_force , MatrixXd& projected_forces )"<<endl;
		throw("");
	}

	MatrixXd projectionMatrix = MatrixXd::Zero( 3 * this->numberOfPoints , 3 * this->numberOfPoints );
	for( unsigned int index = 0 ; index < 3 * this->numberOfPoints ; index ++ )
	{
		projectionMatrix( index , index ) = 1.0;
	}

		// Inversin matrix
	MatrixXd G_uvINV = MatrixXd::Zero( this->numberOfPoints - 1 , this->numberOfPoints - 1 );
        this->getMatrixes( G_uvINV , projectionMatrix );
	projected_forces = projectionMatrix * original_force;
}





void Microtubule::oneStepMidStepAlgorithm_2_half_projected( MatrixXd& original_Coordinates , MatrixXd& randomForces , MatrixXd& other_Forces  )
{
		if( ( randomForces.rows() != this->numberOfPoints * 3 ) || ( randomForces.cols() != 1 )  )
		{
			cout<<"( randomForces.rows() != this->numberOfPoints * 3 ) || ( randomForces.cols() != 1 )  in oneStepMidStepAlgorithm_1_half "<<endl;
			cout<<" randomForces.rows() = "<<randomForces.rows()<<endl;
			cout<<"randomForces.cols() = "<<randomForces.cols()<<endl;
			cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
			throw("");
		}
		if( ( other_Forces.rows() != this->numberOfPoints * 3 ) || ( other_Forces.cols() != 1 )  )
		{
			cout<<"( extrenal_Forces.rows() != this->numberOfPoints * 3 ) || ( extrenal_Forces.cols() != 1 )  in oneStepMidStepAlgorithm_1_half "<<endl;
			cout<<" extrenal_Forces.rows() = "<<other_Forces.rows()<<endl;
			cout<<"extrenal_Forces.cols() = "<<other_Forces.cols()<<endl;
			cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
			throw("");
		}



		if( this->numberOfPoints == 0 )
		{
			cout<<"this->numberOfPoints == 0 in Microtubule::oneStepMidStepAlgorithm_1_half( MatrixXd& randomForces , MatrixXd& extrenal_Forces  )"<<endl;
			throw("");
		}

	if( this->numberOfPoints == 1 )
	{
		//if this->numberOfPoints == 1 there are no geometrical problems! No projection is needed and also metric forces do not exist.

		MatrixXd V_0 = (  other_Forces + randomForces ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original_Coordinates + V_0 * sim_of_Cell::time_Step;
		this->setCoordinates( R_half );
	}
	else
	{
		MatrixXd forces_Projection = randomForces + other_Forces;
		MatrixXd V_0 = ( forces_Projection ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original_Coordinates + V_0 * sim_of_Cell::time_Step;
        	this->setCoordinates( R_half );

	}
}















void Microtubule::Euler_algorithm( MatrixXd& randomForces , MatrixXd& extrenal_Forces  )
{

    unsigned int rows_tmp = extrenal_Forces.rows();
    unsigned int arg_tmp = this->numberOfPoints * 3;


    if( ( randomForces.rows() != this->numberOfPoints * 3 ) || ( randomForces.cols() != 1 )  )
	{
		cout<<"( randomForces.rows() != this->numberOfPoints * 3 ) || ( randomForces.cols() != 1 )  in Euler_algorithm "<<endl;
		cout<<" randomForces.rows() = "<<randomForces.rows()<<endl;
		cout<<"randomForces.cols() = "<<randomForces.cols()<<endl;
		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
		throw("");
	}


    if( ( rows_tmp != arg_tmp )  )
	{
        cout<<"Microtubule ERROR_ID = 13561456464568"<<endl;
		cout<<"( extrenal_Forces.rows() != this->numberOfPoints * 3 ) || ( extrenal_Forces.cols() != 1 )  in Euler_algorithm "<<endl;
		cout<<" extrenal_Forces.rows() = "<<extrenal_Forces.rows()<<endl;
        cout<<"this->numberOfPoints = "<<this->numberOfPoints<<endl;
		cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
        cout<<"this->microtubule_id = "<<this->microtubule_id<<endl;
		throw("");
	}

	if( this->numberOfPoints == 0 )
	{
        cout<<"ERROR_ID = 132561497813613"<<endl;
		cout<<"this->numberOfPoints == 0 in Microtubule::Euler_algorithm( MatrixXd& randomForces , MatrixXd& extrenal_Forces  )"<<endl;
		throw("");
	}
	//Coordinates before step
	MatrixXd original = this->getCoordinates();

	if( this->numberOfPoints == 1 )
	{
		//if this->numberOfPoints == 1 there are no geometrical problems! No projection is needed and also metric forces do not exist.

		MatrixXd V_0 = (  extrenal_Forces + randomForces ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original + V_0 * sim_of_Cell::time_Step;
		this->setCoordinates( R_half );
	}
	else
	{
		//Kronecker delta
		MatrixXd projectionMatrix = MatrixXd::Zero( 3 * this->numberOfPoints , 3 * this->numberOfPoints );
		for( unsigned int index = 0 ; index < 3 * this->numberOfPoints ; index ++ )
		{
			projectionMatrix( index , index ) = 1.0;
		}

		// Inversin matrix
		MatrixXd G_uvINV = MatrixXd::Zero( this->numberOfPoints - 1 , this->numberOfPoints - 1 );
                this->getMatrixes_Sparse( G_uvINV , projectionMatrix );
		// ------------------------------------------------------------------------------------------ Inter Microtubule forces
		//Bending Forces
		MatrixXd bendingForce = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
		this->getBendingForces( bendingForce );

		// ------------------------------------------------------------------------------------------ Projections
		//Projection
		MatrixXd forces = bendingForce + extrenal_Forces + randomForces;
		MatrixXd forces_Projection = projectionMatrix * forces; //
		//cout<<"..............................................................."<<endl;
		MatrixXd V_0 = ( forces_Projection ) * ( 1.0 / this->effective_friction );
		MatrixXd R_half = original + V_0 * sim_of_Cell::time_Step;

        this->setCoordinates( R_half );

	}
}



bool Microtubule::distanceControl()
{
	unsigned int control = 0;
	if( this->getNumberOfPoints() > 1  )
	{

		for( unsigned int i = 0 ; i < this->numberOfPoints - 1 ; i ++ )
		{
			Vector3d tangent = this->getTangent2( i );
			double percentage = abs( ( tangent.norm() - this->restDistancePoints ) / this->restDistancePoints );
			if( percentage > 0.005 )
			{
				control = 1;
				break;
			}
		}
	}
	return control;
}


void Microtubule::resizeMicrotubule()
{
	if( this->getNumberOfPoints() > 2 )
	{
		Vector3d tangents[ this->numberOfPoints - 1 ];
		for( unsigned int i = 0 ; i < this->numberOfPoints - 1 ; i ++ )
		{
			Vector3d tangent = this->getTangent2( i );
			tangents[ i ] = tangent * this->restDistancePoints / tangent.norm();
		}

		//Regrow of microtubule
		//Bead 0 stays at the place - in MTOC - other beads are resized
		for( unsigned int i = 1 ; i< this->numberOfPoints ; i ++ )
		{

			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				coordinates( 3 * i + j , 0 ) = coordinates( 3 * ( i - 1 ) + j , 0 ) + tangents[ i - 1 ]( j );
			}
		}
	}

}


void Microtubule::resizeMicrotubule( Vector3d orientation )
{
    if( this->getNumberOfPoints() > 2 )
	{
		Vector3d tangents[ this->numberOfPoints - 1 ];

        tangents[ 0 ] = orientation / orientation.norm() * this->restDistancePoints ;
		for( unsigned int i = 1 ; i < this->numberOfPoints - 1 ; i ++ )
		{
			Vector3d tangent = this->getTangent2( i );
			tangents[ i ] = tangent / tangent.norm() * this->restDistancePoints ;
		}

		//Regrow of microtubule
		//Bead 0 stays at the place - in MTOC - other beads are resized
		for( unsigned int i = 1 ; i< this->numberOfPoints ; i ++ )
		{
			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				coordinates( 3 * i + j , 0 ) = coordinates( 3 * ( i - 1 ) + j , 0 ) + tangents[ i - 1 ]( j );
			}
		}
	}

}


void Microtubule::resizeMicrotubule_first_segment_different( )
{
    //Tangents array contains resized vector
	if( this->getNumberOfPoints() > 1 )
	{
		Vector3d tangents[ this->numberOfPoints - 1 ];
        Vector3d orientation_1 = this->getTangent2( 0 ) / this->getTangent2( 0 ).norm();

        tangents[ 0 ] = orientation_1 * this->restDistancePoints_first;
		for( unsigned int i = 1 ; i < this->numberOfPoints - 1 ; i ++ )
		{
			Vector3d tangent = this->getTangent2( i );
            Vector3d orientation = tangent / tangent.norm();
			tangents[ i ] = orientation * this->restDistancePoints;
		}

		//Regrow of microtubule
		//Bead 0 stays at the place - in MTOC - other beads are resized
		for( unsigned int i = 1 ; i< this->numberOfPoints ; i ++ )
		{
			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				coordinates( 3 * i + j , 0 ) = coordinates( 3 * ( i - 1 ) + j , 0 ) + tangents[ i - 1 ]( j );
			}
		}
	}
}

void Microtubule::resizeMicrotubule_with_different_tangent_lenghts( )
{
	if( this->getNumberOfPoints() > 1 )
	{
		Vector3d tangents[ this->numberOfPoints - 1 ];
        	Vector3d orientation_1 = this->getTangent2( 0 ) / this->getTangent2( 0 ).norm();


        	tangents[ 0 ] = orientation_1 * this->restDistancePoints_first;
		for( unsigned int i = 0 ; i < this->numberOfPoints - 1 ; i ++ )
		{
			Vector3d tangent = this->getTangent2( i );
           		Vector3d orientation = tangent / tangent.norm();
			tangents[ i ] = orientation * this->lenght_of_tangents( i , 0 );
		}

		//Regrow of microtubule
		//Bead 0 stays at the place - in MTOC - other beads are resized
		for( unsigned int i = 1 ; i< this->numberOfPoints ; i ++ )
		{
			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				coordinates( 3 * i + j , 0 ) = coordinates( 3 * ( i - 1 ) + j , 0 ) + tangents[ i - 1 ]( j );
			}
		}
	}
}


void Microtubule::midStepAlgorithm( double startTime , double endTime )
{
	unsigned int numberOfSteps = ( unsigned int ) ( ( endTime - startTime ) / sim_of_Cell::time_Step );
	for( unsigned int i = 0 ; i < numberOfSteps ; i ++ )
	{
		this->oneStepMidStepAlgorithm();

		if( i % 100 == 0)
		{
			if( distanceControl() == 1 )
			{
                this->print_Microtubule_Tangent();
				this->resizeMicrotubule();
			}
		}

	}
}





void Microtubule::midStepAlgorithm( double startTime , double endTime , IS_Capt_Shrinkage  IS_Capt_Shrinkage_argument  )
{
	unsigned int numberOfSteps = ( unsigned int ) ( ( endTime - startTime ) / sim_of_Cell::time_Step );


	for( unsigned int i = 0 ; i < numberOfSteps ; i ++ )
	{

		if( i % 100 == 0)
		{

			this->IS_micro_control_one_micro_capture_shrinkage( IS_Capt_Shrinkage_argument );
			if( distanceControl() == 1 )
			{
				this->resizeMicrotubule();
			}
		}
		this->oneStepMidStepAlgorithm();

	}


}







MatrixXd Microtubule::relaxationTimeMidStep( double startTime , double endTime , double DELTATT )
{

	Vector3d reference( 0.0 , 0.0 , 0.0 );
	reference( 0 ) = 1.0;
	unsigned int numSteps = ( ( endTime - startTime ) / DELTATT );
	MatrixXd equipartitionTheorem = MatrixXd::Zero( this->numberOfPoints - 2 , 1 );

	cout<<"numSteps = "<<numSteps<<endl;

	for( unsigned int i = 0 ; i < numSteps ; i ++ )
	{
		cout<<" i = "<<i<<endl;
		this->midStepAlgorithm( 0.0 , DELTATT  );
		for( unsigned int index = 0 ; index < this->numberOfPoints - 2; index ++ )
		{
			Vector3d first = this->getTangent2( index );
			double uhel1 = acos( first.dot( reference ) / first.norm() );		//I take two angles and measure difference
			Vector3d second = this->getTangent2( index + 1 );
			double uhel2 = acos( second.dot( reference ) / second.norm() );
			double uhel3 = (uhel1 - uhel2);
			equipartitionTheorem( index , 0 ) = equipartitionTheorem( index , 0 ) + uhel3 * uhel3;
		}
	}

	equipartitionTheorem = equipartitionTheorem * ( 1.0 / ( double ) numSteps );
	return equipartitionTheorem;

}


double Microtubule::equipartition_Theorem_analytical()
{
	return sim_of_Cell::Temperature * sim_of_Cell::boltzmann_constant / this->kappa;
}





void Microtubule::Substract_one_Bead()
{
	if( this->numberOfPoints == 1  )
	{
		cout<<"this->numberOfPoints == 1  and I call Microtubule::Substract_one_Bead()"<<endl;
		throw("");
	}

    if( this->numberOfPoints > 2 )
    {
        this->numberOfPoints = this->numberOfPoints - 1;
        if( this->numberOfPoints == 1 )
        {
            cout<<".............SUBSTRACTION....."<<endl;
        }

        MatrixXd tmp = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
        for( unsigned int i = 0 ; i < 3 * this->numberOfPoints ; i ++ )
        {
            tmp( i , 0 ) = this->coordinates( i , 0 );
        }

        this->coordinates = tmp;
    }
    else
    {
        this->set_dynein_index( 2 );
    }
}


void Microtubule::add_one_final_Bead_catching_position()
{
        this->numberOfPoints = this->numberOfPoints + 1;
        if( this->numberOfPoints == 1 )
        {
            cout<<".............SUBSTRACTION....."<<endl;
        }
        MatrixXd tmp = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
        for( unsigned int i = 0 ; i < 3 * ( this->numberOfPoints - 1 ) ; i ++ )
        {
            tmp( i , 0 ) = this->coordinates( i , 0 );
        }

        for( unsigned int i = 0 ; i < 3 ; i ++ )
        {
            tmp( 3 * ( this->numberOfPoints - 1 ) + i , 0 ) = this->IS_position_catching( i );
        }

        this->coordinates = tmp;

}





void Microtubule::Substract_All_Bead()
{
	this->numberOfPoints = 1;
	MatrixXd tmp = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
	for( unsigned int i = 0 ; i < 3 ; i ++ )
	{
		tmp( i , 0 ) = this->coordinates( i , 0 );
	}

	this->coordinates = tmp;

}

void Microtubule::smoother_Microtubule()
{
	if( this->numberOfPoints == 1  )
	{
		cout<<" Microtubule::smoother_Microtubule(): this->numberOfPoints == 1 "<<endl;
		throw("");
	}
	unsigned int new_number_of_bead = 2 * this->numberOfPoints - 1;
	MatrixXd tmp = MatrixXd::Zero( 3 * new_number_of_bead , 1 );

	for( unsigned int i = 0 ; i < this->numberOfPoints - 1 ; i ++ )
	{
		Vector3d point = this->getPoint( i );
		Vector3d tangent = this->getTangent2( i );
		for( unsigned int j = 0 ; j < 3 ; j ++ )
		{
			tmp( 3 * 2 * i + j , 0 ) = point( j );
			tmp( 3 * ( 2 * i + 1 ) + j , 0 ) = point( j ) + tangent( j ) / 2.0;
		}
	}

	Vector3d point = this->getPoint( this->numberOfPoints - 1 );
	for( unsigned int j = 0 ; j < 3 ; j ++ )
	{
		tmp( 3 * 2 * ( this->numberOfPoints - 1 ) + j , 0  ) = point( j );
	}

	this->numberOfPoints = new_number_of_bead;
	this->coordinates = tmp;
	this->restDistancePoints = this->restDistancePoints / 2.0;

}

void Microtubule::destroy_Microtubule_capture_shrinkage()
{
	if( this->numberOfPoints != 1 )
	{
		cout<<"this->numberOfPoints != 1 in Microtubule::destroy_Microtubule_capture_shrinkage()"<<endl;
		throw("");
	}
	this->set_dynein_index( 2 );
	this->numberOfPoints = 0;

}

void Microtubule::destroy_Microtubule_cortical_sliding()
{
	if( this->numberOfPoints != 1 )
	{
		cout<<"this->numberOfPoints != 1 in Microtubule::destroy_Microtubule_cortical_sliding()"<<endl;
		cout<<"this->numberOfPoints = "<<this->numberOfPoints <<endl;
		throw("");
	}
	this->set_dynein_index( 5 );
	this->numberOfPoints = 0;
}






void Microtubule::IS_micro_control_one_micro_capture_shrinkage( IS_Capt_Shrinkage IS_Capt_Shrinkage_argument )
{
	unsigned int dynein_index = this->get_dynein_index();
	if( dynein_index == 0 )
	{
		for( unsigned int bead = 0 ; bead < this->getNumberOfPoints() ; bead ++ )
		{

			Vector3d bead_position = this->getPoint( bead );

			Vector3d tangent_Center_IS_bead = ( bead_position - IS_Capt_Shrinkage_argument.get_center_of_IS_front() );
			double cos_two_lines = tangent_Center_IS_bead.dot( IS_Capt_Shrinkage_argument.get_axis_of_IS() ) / ( tangent_Center_IS_bead.norm() * IS_Capt_Shrinkage_argument.get_axis_of_IS().norm() );
			if( cos_two_lines < 0 )
			{
				continue;
			}

			else
			{

				if( tangent_Center_IS_bead.norm() > IS_Capt_Shrinkage_argument.get_radius_of_IS() )
				{
					continue;
				}

				else
				{

					this->set_dynein_index( 1 );
					unsigned int number_of_substraction = this->getNumberOfPoints() - bead;
					if( number_of_substraction < this->getNumberOfPoints() )
					{
						for( unsigned int substraction = 0 ; substraction < number_of_substraction ; substraction ++ )
						{
							this->set_IS_position_catching( bead_position );
							this->Substract_one_Bead();
						}
						if( this->getNumberOfPoints() > 1  )
						{
							this->smoother_Microtubule();
						}
						break;
					}
					if( number_of_substraction == this->getNumberOfPoints() )
					{
						for( unsigned int substraction = 0 ; substraction < number_of_substraction - 1 ; substraction ++ )
						{
							this->set_IS_position_catching( bead_position );
							this->Substract_one_Bead();
						}
					}

				}
			}
		}

	}
	else if ( ( dynein_index == 1 ) && ( this->getNumberOfPoints() ) > 1 )
	{
		Vector3d bead_position = this->getPoint(  this->getNumberOfPoints() - 1 );

		Vector3d catching_point = this->get_IS_position_catching();
		if( ( bead_position - catching_point ).norm() < 1e-8 )
		{
			this->Substract_one_Bead();
		}
	}

	else if ( ( dynein_index == 1 ) && ( this->getNumberOfPoints() ) == 1 )
	{

		Vector3d bead_position = this->getPoint( 0 );
		Vector3d tangent_Center_IS_bead = ( bead_position - IS_Capt_Shrinkage_argument.get_center_of_IS_front() );

		if( tangent_Center_IS_bead.norm() > 0.2e-7 )
		{

		}
		else
		{
			cout<<"this->array_Of_Microtubules[ microtubule ].destroy_Microtubule();"<<endl;
			this->destroy_Microtubule_capture_shrinkage();
			cout<<"this->array_Of_Microtubules[ microtubule ].get_dynein_index() = "<<this->get_dynein_index()<<endl;
		}

	}

}





void Microtubule::print_Microtubule()
{

	for( unsigned int i = 0 ; i < this->numberOfPoints ; i ++ )
	{
		cout<<" i = "<<i<<endl;
		cout<<this->getPoint( i )<<endl;
	}

}

void Microtubule::print_Microtubule_Tangent()
{
    for( unsigned int i = 0 ; i < this->numberOfPoints - 1 ; i ++ )
	{
		cout<<" i = "<<i<<endl;
		cout<<this->getTangent2( i ).norm()<<endl;
	}
	cout<<"resting_distance = "<<this->getRestDist()<<endl;

}




bool Microtubule::confirm_inner_ellipsoid( Vector3d position , double a_axis , double b_axis )
{

	if( position( 2 ) > 0 )
	{
		double confirm_ellipsoid_value = ( position( 0 ) * position( 0 ) + position( 1 ) * position( 1 ) ) / ( a_axis * a_axis );
		confirm_ellipsoid_value = confirm_ellipsoid_value + ( position( 2 ) * position( 2 ) ) / ( b_axis * b_axis );
		if( confirm_ellipsoid_value < 1.0 )
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	if( position( 2 ) <= 0 )
	{
		b_axis = b_axis - 1.0e-6;
		double confirm_ellipsoid_value = ( position( 0 ) * position( 0 ) + position( 1 ) * position( 1 ) ) / ( a_axis * a_axis );
		confirm_ellipsoid_value = confirm_ellipsoid_value + ( position( 2 ) * position( 2 ) ) / ( b_axis * b_axis );
		if( confirm_ellipsoid_value < 1.0 )
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	return false;

}



void Microtubule::distribute_force_on_beads( MatrixXd& force_on_microtubule , Vector3d force , unsigned int bead_segment , double t )
{
    if( force.norm() < 1e-16 )
    {
        return;
    }

    if( bead_segment >= this->numberOfPoints )
    {
        cout<<"bead_segment >= this->numberOfPoints"<<endl;
        cout<<"Microtubule ERROR_ID = 646161654346468"<<endl;
        cout<<"bead_segment = "<<bead_segment<<endl;
        cout<<"this->numberOfPoints = "<<this->numberOfPoints<<endl;
        cout<<"this->get_dynein_index() = "<<this->get_dynein_index()<<endl;
        throw("");
    }

        Vector3d tangent = this->getTangent2( bead_segment );
        double cosinus = tangent.dot( force ) / ( tangent.norm() * force.norm() );

    //PARALLEL FORCES
        Vector3d parallel_forces = ( force.norm() * cosinus ) * ( tangent / tangent.norm() );


        for( unsigned int i = 0 ; i < 3 ; i ++ )
        {
            force_on_microtubule( 3 * bead_segment + i , 0 ) = force_on_microtubule( 3 * bead_segment + i , 0 ) + parallel_forces( i ) / 2.0;
            force_on_microtubule( 3 * ( bead_segment + 1 ) + i , 0 ) = force_on_microtubule( 3 * bead_segment + i , 0 ) + parallel_forces( i ) / 2.0;
        }

    //PERPENDICULAR FORCES
        Vector3d perpendicullar_forces = force - parallel_forces;
        for( unsigned int i = 0 ; i < 3 ; i ++ )
        {
            force_on_microtubule( 3 * bead_segment + i , 0 ) = force_on_microtubule( 3 * bead_segment + i , 0 ) + perpendicullar_forces( i ) / 2.0;
            force_on_microtubule( 3 * ( bead_segment + 1 ) + i , 0 ) = force_on_microtubule( 3 * bead_segment + i , 0 ) + perpendicullar_forces( i ) / 2.0;
        }

    //MOMENT FORCES
        double ratio = abs( 0.5 - t );
        Vector3d moment_forces = perpendicullar_forces * ratio;

        if( t < 0.5 )
        {
            for( unsigned int i = 0 ; i < 3 ; i ++ )
            {
                force_on_microtubule( 3 * bead_segment + i , 0 ) = force_on_microtubule( 3 * bead_segment + i , 0 ) + moment_forces( i ) / 2.0;
                force_on_microtubule( 3 * ( bead_segment + 1 ) + i , 0 ) = force_on_microtubule( 3 * bead_segment + i , 0 ) - moment_forces( i ) / 2.0;
            }
        }
        else
        {
            for( unsigned int i = 0 ; i < 3 ; i ++ )
            {
                force_on_microtubule( 3 * bead_segment + i , 0 ) = force_on_microtubule( 3 * bead_segment + i , 0 ) - moment_forces( i ) / 2.0;
                force_on_microtubule( 3 * ( bead_segment + 1 ) + i , 0 ) = force_on_microtubule( 3 * bead_segment + i , 0 ) + moment_forces( i ) / 2.0;
            }
        }

}






//MICRO INTERACTION WITH CORTICAL SLIDING ON THE SURFACE with abscissa
void Microtubule::add_pair( std::pair < Vector3d ,double  > point_abscissa )
{
    this->Dynein_motors_2.push_back( point_abscissa );
}



std::vector< Vector3d  > Microtubule::get_Dynein_points_and_erase()
{
    if( this->get_dynein_index() != 9 )
    {
        cout<<"this->get_dynein_index() != 9"<<endl;
        cout<<"std::vector< Vector3d  > Microtubule::get_Dynein_points_and_erase()"<<endl;
        cout<<"ERROR_ID = 6835461681364813"<<endl;
        throw("");
    }

    std::vector< Vector3d  > point_coordinates;
    for( unsigned int pair_index = 0 ; pair_index < this->Dynein_motors_2.size() ; pair_index ++ )
    {
        std::pair < Vector3d , double > pair_tmp = this->Dynein_motors_2.at( pair_index );
        Vector3d anchor_position = std::get<0>( pair_tmp );
        point_coordinates.push_back( anchor_position );
    }
    this->Dynein_motors_2.clear();
    return point_coordinates;
}





Vector3d Microtubule::force_real_dynein_one_pair(  std::pair < Vector3d , double > pair_tmp   )
{
    Vector3d force( 0.0 , 0.0 , 0.0 );
    Vector3d anchor_position = std::get<0>( pair_tmp );
    double abscissa = std::get<1>( pair_tmp );
    Vector3d point_of_attachment = get_attachment_point_according_to_abscissa( abscissa );
    Vector3d ancher_atachment = anchor_position - point_of_attachment;

    double distance_ancher_atachment = ancher_atachment.norm();

    if( distance_ancher_atachment < Dynein_real::L_0 )
    {

    }
    else
    {
        double distance_spring = distance_ancher_atachment - Dynein_real::L_0;
        double force_absolute_value = Dynein_real::alpha * distance_spring;
        force = force_absolute_value * ancher_atachment / ancher_atachment.norm();

    }
    return force;


}







std::vector< Vector3d > Microtubule::stepping_detach_real_dynein_abscissa_projection( )
{
    std::uniform_real_distribution<> distribution{ 0 , 1 };
    unsigned int number_of_generator = omp_get_thread_num();






    std::vector< Vector3d > return_vectors;
    std::vector< std::pair < Vector3d , double > > new_vectors;


    for( unsigned int pair_index = 0 ; pair_index < this->Dynein_motors_2.size() ; pair_index ++ )
    {
        std::pair < Vector3d ,double > pair_tmp = this->Dynein_motors_2[ pair_index ];
        Vector3d anchor_position = std::get<0>( pair_tmp );
        double abscissa = std::get<1>( pair_tmp );

        Vector3d force_one_pair = this->force_real_dynein_one_pair(  pair_tmp  );
	double force = force_one_pair.norm();

        unsigned int lower_bead_index = 0;
        if( abscissa < this->getTangent2( 0 ).norm() )
        {
            lower_bead_index = 0;
        }
        else
        {
            unsigned int lower_bead_index;
            try
            {
                lower_bead_index = this->get_index_according_to_abscissa( abscissa );
            }
            catch( unsigned int error_id )
            {
                cout<<"std::vector< Vector3d > Microtubule::stepping_detach_real_dynein_abscissa_projection( )"<<endl;
                throw("");
            }
        }


        Vector3d tangent( 0.0 , 0.0 , 0.0 );
        try
        {
            if( abscissa < this->get_lenght_of_microtubule() )
            {
                tangent = this->getTangent2( lower_bead_index );
            }
            else if( abscissa == this->get_lenght_of_microtubule() )
            {
                tangent = ( -1.0 ) * this->getTangent2( this->getNumberOfPoints() - 2 );
            }
        }
        catch ( string e)
        {
            cout<<"error in std::vector< Vector3d > Microtubule::stepping_detach_real_dynein_abscissa_projection( )"<<endl;
        }

        double detach_prob = Dynein_real::detachment_probability_per_step( force );//NARAZNIK

	double probability = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
        if( probability < detach_prob )
        {
            return_vectors.push_back( anchor_position );
            continue;
        }

        if( force == 0 )
        {
 	    double probability_2 = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
            if( probability_2 < Dynein_real::forward_stepping_probability_per_step )
            {
                abscissa = abscissa - Dynein_real::step;
            }
        }
        else
        {
            double cosinus = tangent.dot( force_one_pair ) / ( force_one_pair.norm() * tangent.norm() );
            if( cosinus < 0 )
            {

                double probability_2 = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
                if( probability_2 < Dynein_real::forward_stepping_probability_per_step )
                {
                    abscissa = abscissa - Dynein_real::step;
                }
            }
            else
            {
                if( force < Dynein_real::F_stall )
                {
                    double prav_tmp = Dynein_real::backward_stepp_prob_smaller_than_Stall_f_per_step( force );
                    double probability_2 = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
                    if( probability_2 < prav_tmp )
                    {
                        abscissa = abscissa - Dynein_real::step;
                    }

                }
                else
                {
                    double tmp = Dynein_real::backward_prob_force_bigger_then_Stall_force_per_step;
                    double probability_2 = distribution( mersenne_twisters_64::vector_of_mersenne_twisters[ number_of_generator ] );
                    if( probability_2 < tmp )
                    {
                        abscissa = abscissa + Dynein_real::step;
                    }

                }
            }
        }

        if( abscissa <= this->getTangent2( 0 ).norm()  )
        {
            return_vectors.push_back( anchor_position );
        }
        else if( abscissa > this->get_lenght_of_microtubule() )
        {
            return_vectors.push_back( anchor_position );
        }
        else
        {
            std::pair < Vector3d ,double > pair_2( anchor_position , abscissa );
            new_vectors.push_back( pair_2 );
        }

    }

    this->Dynein_motors_2 = new_vectors;
    if( new_vectors.size() == 0 )
    {
        if( this->get_dynein_index() == 9 )
        {
            this->set_dynein_index( 0 );
        }
    }

    return return_vectors;
}

















void Microtubule::force_dynein_abscissa_real_dynein( MatrixXd& force_dynein )
{

    if( ( force_dynein.rows() != 3 * this->getNumberOfPoints() ) || ( force_dynein.cols() != 1 )  )
    {
        cout<<"( force_dynein.rows() != 3 * this->getNumberOfPoints() ) || ( force_dynein.cols() != 1 )"<<endl;
        cout<<"Microtubule::force_dynein_abscissa_real_dynein( MatrixXd& force_dynein )"<<endl;
        cout<<"Microtubule ERROR_ID = 3516145415461115626"<<endl;
        throw("");
    }


    for( unsigned int pair_index = 0 ; pair_index < this->Dynein_motors_2.size() ; pair_index ++ )
    {
        std::pair < Vector3d , double > pair_tmp = this->Dynein_motors_2.at( pair_index );
        Vector3d force = this->force_real_dynein_one_pair( pair_tmp );

        Vector3d anchor_position = std::get<0>( pair_tmp );
        double abscissa = std::get<1>( pair_tmp );

        unsigned int index_lower_bead = this->get_index_according_to_abscissa( abscissa );


        if( index_lower_bead < this->getNumberOfPoints() - 1 )
        {

            Vector3d tangent = this->getTangent2( index_lower_bead );
            double abscissa_minus_lower_bead = abscissa - get_distance_to_lower_bead_with_index( index_lower_bead );
            double t_parameter = abscissa_minus_lower_bead / tangent.norm();
            this->distribute_force_on_beads(  force_dynein , force , index_lower_bead , t_parameter );
        }
        else if( index_lower_bead == this->getNumberOfPoints() - 1 )
        {
            double t_parameter = 1.0;
            this->distribute_force_on_beads(  force_dynein , force ,  index_lower_bead - 1 , t_parameter );
        }

    }

}

void Microtubule::set_to_number_of_points( unsigned int number_of_point )
{
    if( number_of_point > this->numberOfPoints )
    {
        cout<<"void Microtubule::set_to_number_of_points( unsigned int number_of_point )"<<endl;
        cout<<"number_of_point = "<<number_of_point<<endl;
        cout<<"this->getNumberOfPoints() = "<<this->getNumberOfPoints()<<endl;
        cout<<"ERROR_ID = 136841246816161"<<endl;
        throw("");
    }


    this->numberOfPoints = number_of_point;
    MatrixXd tmp = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
    for( unsigned int i = 0 ; i < 3 * this->numberOfPoints ; i ++ )
    {
        tmp( i , 0 ) = this->coordinates( i , 0 );
    }

    this->coordinates = tmp;
}

unsigned int Microtubule::get_index_according_to_abscissa( double abscissa )
{

    if( abscissa > this->get_lenght_of_microtubule() )
    {
        cout<<"............................................"<<endl;
        cout<<this->get_polygon_number()<<endl;
        cout<<this->getID()<<endl;
        cout<<"dynein index = "<<this->get_dynein_index()<<endl;
        cout<<"abscissa >= this->get_lenght_of_microtubule()"<<endl;
        double diffference = abscissa - this->get_lenght_of_microtubule();
        printf( "%.50f", diffference );
        cout<<endl;


        cout<<"abscissa = "<<endl;
        printf( "%.50f\n", abscissa );
        cout<<"lenght = "<<endl;
        printf( "%.50f\n", this->get_lenght_of_microtubule() );
        cout<<"time = "<<Cell_parametres::time<<endl;


        unsigned int ERROR_ID = 81346;
        cout<<"unsigned int Microtubule::get_index_according_to_abscissa( double abscissa )"<<endl;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
        throw ERROR_ID;
    }
    else if( abscissa == this->get_lenght_of_microtubule() )
    {
        return ( this->numberOfPoints - 1);
    }
    else
    {



        double lower_length = 0;
        double upper_length = 0;
        for( unsigned int segment_id = 0 ; segment_id < this->getNumberOfPoints() - 1 ; segment_id ++ )
        {
            double length_of_segment = this->getTangent2( segment_id ).norm();
            upper_length = upper_length + length_of_segment;
            if( ( abscissa < upper_length ) && ( abscissa >= lower_length ) )
            {
                return segment_id;
            }

            lower_length = lower_length + length_of_segment;
        }
    }
    return 0;
}


Vector3d Microtubule::get_attachment_point_according_to_abscissa( double abscissa )
{
    if( abscissa == this->get_lenght_of_microtubule() )
    {
        return this->getPoint( this->numberOfPoints - 1 );
    }
    else
    {
        unsigned int lower_index;
        try
        {
            lower_index = this->get_index_according_to_abscissa( abscissa );
        }
        catch( unsigned int error_id )
        {
            cout<<"std::vector< Vector3d > Microtubule::get_attachment_point_according_to_abscissa( )"<<endl;
            throw("");
        }

        double lower_length = 0;
        for( unsigned int segment = 0 ; segment < lower_index ; segment ++ )
        {
            lower_length = lower_length + this->getTangent2( segment ).norm();
        }

        double abscissa_minus_lower_bead = abscissa - lower_length;
        Vector3d micro_bead = this->getPoint( lower_index );
        Vector3d tangent = this->getTangent2( lower_index );
        Vector3d attachment_point = micro_bead + tangent / tangent.norm() * abscissa_minus_lower_bead;
        return attachment_point;
    }
}



void Microtubule::set_lenght_of_tangents()
{
    this->lenght_of_tangents = MatrixXd::Zero( this->numberOfPoints - 1 , 1 );
    for( unsigned int counter = 0 ; counter < this->numberOfPoints - 1 ; counter ++ )
    {
        Vector3d tangent;
        try
        {
            tangent = this->getTangent2( counter );
        }
        catch( int e )
        {
            cout<<"exception"<<endl;
            cout<<"void void Microtubule::set_lenght_of_tangents()"<<endl;
            throw("");
        }
        this->lenght_of_tangents( counter , 0 ) = this->getTangent2( counter ).norm();
    }
}


std::vector<Vector3d> Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule()
{
    if( this->get_dynein_index() != 20 )
    {
        cout<<"Microtubule::control_and_resizing_of_depolimerizing_microtubule()"<<endl;
        cout<<" this->get_dynein_index() != 20"<<endl;
        unsigned int ERROR_ID = 5614681;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
        throw ERROR_ID;
    }




    //first
    this->resizeMicrotubule_with_different_tangent_lenghts();

    //second
    std::vector<Vector3d> detached_motors = this->stepping_detach_real_dynein_abscissa_projection();

    if( this->getNumberOfPoints() > 3 )
    {
        if( this->getTangent2( this->numberOfPoints - 2 ).norm() < this->restDistancePoints / 2.0   )
        {
            //one bead is removed - numberOfPoints, coordinates and lenghtofSegments have to be redone
            this->numberOfPoints = this->numberOfPoints - 1;
            MatrixXd tmp_new = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
            for( unsigned int counter = 0 ; counter < 3 * ( this->numberOfPoints - 1 ) ; counter ++ )
            {
                tmp_new( counter , 0 ) = this->coordinates( counter , 0 );
            }
            for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
            {
                tmp_new( 3 * ( this->numberOfPoints - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
            }
            this->coordinates = tmp_new;

            cout<<"this->numberOfPoints = "<<this->numberOfPoints<<endl;




        }
        else
        {
            //setting the last point
            for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
            {
                this->coordinates( 3 * ( this->getNumberOfPoints() - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );

            }
        }
        this->set_lenght_of_tangents();
    }

    else if( this->getNumberOfPoints() == 3 )
    {
            //setting the last point
        for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
        {
            this->coordinates( 3 * ( this->getNumberOfPoints() - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
        }

        this->set_lenght_of_tangents();





    }



    std::vector<Vector3d> vectors_to_be_returned_to_IS;

    for( unsigned int counter = 0 ; counter < detached_motors.size() ; counter ++ )
    {
        vectors_to_be_returned_to_IS.push_back( detached_motors[ counter ] );
    }

    std::vector<Vector3d> kicked_due_to_repolarization = this->control_motor_detachment_with_depolimerization();
    for( unsigned int counter = 0 ; counter < kicked_due_to_repolarization.size() ; counter ++ )
    {
        vectors_to_be_returned_to_IS.push_back( kicked_due_to_repolarization[ counter ] );
    }


    if( this->Dynein_motors_2.size() == 0 )
    {
        //cout<<"SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP SETUP"<<endl;
        //this->set_dynein_index( 0 );
    }
    if( this->getNumberOfPoints() <= 2 )
    {
        this->set_dynein_index( 0 );
    }


    return vectors_to_be_returned_to_IS;


}







std::vector<Vector3d>  Microtubule::control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_2()
{

    if( this->get_dynein_index() != 20 )
    {
        cout<<"Microtubule::control_and_resizing_of_depolimerizing_microtubule()"<<endl;
        cout<<" this->get_dynein_index() != 20"<<endl;
        unsigned int ERROR_ID = 5614681;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
        throw ERROR_ID;
    }

    this->resizeMicrotubule_with_different_tangent_lenghts();
    std::vector<Vector3d> detached_motors = this->stepping_detach_real_dynein_abscissa_projection();


    if( this->numberOfPoints > 3 )
    {

       if( this->getTangent2( this->numberOfPoints - 2 ).norm() < 0.9 * sim_of_Cell::resting_distance  )
       {
	    double lenght_micro_outside_MTOC_2 = this->get_lenght_of_microtubule_outside_MTOC();
	    unsigned int number_of_points_outside_MTOC_2 = this->getNumberOfPoints() - 2;
	    double a_a = abs( lenght_micro_outside_MTOC_2 / ( double ) number_of_points_outside_MTOC_2 - sim_of_Cell::resting_distance  );
	    double b_b = abs( lenght_micro_outside_MTOC_2 / ( double ) ( number_of_points_outside_MTOC_2 - 1.0 ) - sim_of_Cell::resting_distance );

	   if( b_b <= a_a )
	   {

		std::vector< Vector3d > orientations;
	    	for( unsigned int point_id = 1 ; point_id < this->getNumberOfPoints() - 1 ; point_id ++ )
	    	{
			orientations.push_back( this->getTangent2( point_id ) / this->getTangent2( point_id ).norm() );
	    	}
	    	double lenght_micro_outside_MTOC = this->get_lenght_of_microtubule_outside_MTOC();

	    	//zmenim pocet beadu
	   	this->numberOfPoints = this->numberOfPoints - 1;
	    	unsigned int number_of_points_outside_MTOC = this->getNumberOfPoints() - 3;
		if( this->getNumberOfPoints() - 3 > 0 )
		{
	    		this->restDistancePoints = lenght_micro_outside_MTOC / ( double ) ( number_of_points_outside_MTOC );
		}
		else
		{
			this->restDistancePoints = lenght_micro_outside_MTOC;
	    	}

            	MatrixXd tmp_new = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
    	    	for( unsigned int counter = 0 ; counter < 6 ; counter ++ )
    	    	{
			tmp_new( counter , 0 ) = this->coordinates( counter , 0 );
    	    	}



    	   	for( unsigned int i = 2 ; i < this->numberOfPoints ; i ++ )
    	   	{
			Vector3d orient = this->restDistancePoints * orientations[ i - 2 ];
			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				tmp_new( 3 * i + j , 0 ) = tmp_new( 3 * ( i - 1 ) + j , 0 ) + orient( j );
			}
    	   	}

    	  	for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	  	{
			tmp_new( 3 * ( this->numberOfPoints - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
    	  	}

    	  	this->coordinates = tmp_new;


	   }
           else
	   {

	
		std::vector< Vector3d > orientations;
	    	for( unsigned int point_id = 1 ; point_id < this->getNumberOfPoints() - 1 ; point_id ++ )
	    	{
			orientations.push_back( this->getTangent2( point_id ) / this->getTangent2( point_id ).norm() );
	    	}
	    	double lenght_micro_outside_MTOC = this->get_lenght_of_microtubule_outside_MTOC();

	   	this->numberOfPoints = this->numberOfPoints;
	    	unsigned int number_of_points_outside_MTOC = this->getNumberOfPoints() - 2;
	    	this->restDistancePoints = lenght_micro_outside_MTOC / ( double ) ( number_of_points_outside_MTOC );

            	MatrixXd tmp_new = MatrixXd::Zero( 3 * this->numberOfPoints , 1 );
    	    	for( unsigned int counter = 0 ; counter < 6 ; counter ++ )
    	    	{
			tmp_new( counter , 0 ) = this->coordinates( counter , 0 );
    	    	}



    	   	for( unsigned int i = 2 ; i < this->numberOfPoints ; i ++ )
    	   	{
			Vector3d orient = this->restDistancePoints * orientations[ i - 2 ];
			for( unsigned int j = 0 ; j < 3 ; j ++ )
			{
				tmp_new( 3 * i + j , 0 ) = tmp_new( 3 * ( i - 1 ) + j , 0 ) + orient( j );
			}
    	   	}

    	  	for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	  	{
			tmp_new( 3 * ( this->numberOfPoints - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
    	  	}

    	  	this->coordinates = tmp_new;


    	       for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	       {
		    this->coordinates( 3 * ( this->numberOfPoints - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
    	       }


	   }


       }
       else
       {
    	   for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    	   {
		this->coordinates( 3 * ( this->numberOfPoints - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
    	   }
       }
   }
   else
   {
       for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
       {
	   this->coordinates( 3 * ( this->numberOfPoints - 1 ) + dimension , 0 ) = this->IS_position_catching( dimension );
       }
   }



    this->set_lenght_of_tangents();




    std::vector<Vector3d> vectors_to_be_returned_to_IS;

    for( unsigned int counter = 0 ; counter < detached_motors.size() ; counter ++ )
    {
        vectors_to_be_returned_to_IS.push_back( detached_motors[ counter ] );
    }

    std::vector<Vector3d> kicked_due_to_repolarization = this->control_motor_detachment_with_depolimerization();


    for( unsigned int counter = 0 ; counter < kicked_due_to_repolarization.size() ; counter ++ )
    {
        vectors_to_be_returned_to_IS.push_back( kicked_due_to_repolarization[ counter ] );
    }


    if( this->Dynein_motors_2.size() == 0 )
    {

    }
    if( this->getNumberOfPoints() <= 2 )
    {
        this->set_dynein_index( 0 );
    }

    return vectors_to_be_returned_to_IS;



}


























std::vector<Vector3d> Microtubule::control_motor_detachment_with_depolimerization()
{


    std::vector< Vector3d > kicked_vector;
    std::vector< std::pair < Vector3d , double > > remaining_vector;
    for( unsigned int pair_index = 0 ; pair_index < this->Dynein_motors_2.size() ; pair_index ++ )
    {
        std::pair < Vector3d ,double > pair_tmp = this->Dynein_motors_2[ pair_index ];
        double abscissa = std::get< 1 >( pair_tmp );
        Vector3d point = std::get< 0 >( pair_tmp );

        if( abs( abscissa - this->get_lenght_of_microtubule() ) < Dynein_real::step / 10.0 )
        {
            abscissa = abscissa - Dynein_real::step / 10.0;
        }

        std::pair < Vector3d ,double > pair_tmp_2( point , abscissa );


        if( abscissa > this->get_lenght_of_microtubule() )
        {
            kicked_vector.push_back( point );
        }
        else
        {
            remaining_vector.push_back( pair_tmp_2 );
        }
    }
    this->Dynein_motors_2 = remaining_vector;
    return kicked_vector;
}




std::vector<Vector3d> Microtubule::get_dynein_points_in_IS_and_erase()
{
    if( this->get_dynein_index() != 20 )
    {
        cout<<"this->get_dynein_index() != 20"<<endl;
        cout<<"std::vector< Vector3d  > Microtubule::get_dynein_points_in_IS_and_erase()"<<endl;
        unsigned int ERROR_ID = 781646134;
        cout<<"ERROR_ID = "<<ERROR_ID<<endl;
        throw ERROR_ID;
    }

    std::vector< Vector3d  > point_coordinates;
    for( unsigned int pair_index = 0 ; pair_index < this->Dynein_motors_2.size() ; pair_index ++ )
    {
        std::pair < Vector3d , double > pair_tmp = this->Dynein_motors_2.at( pair_index );
        Vector3d anchor_position = std::get<0>( pair_tmp );
        point_coordinates.push_back( anchor_position );
    }
    this->Dynein_motors_2.clear();
    return point_coordinates;

}

unsigned int Microtubule::get_number_of_dynein_points_IS()
{
    return this->Dynein_motors_2.size();
}





std::vector<Vector3d> Microtubule::control_length_of_micro_IS_2()
{

    if( this->dydein_index != 20 )
    {
	cout<<"this->dydein_index != 20"<<endl;
	cout<<"ERROR_ID = 5468435615646"<<endl;
	throw("");
    }

    double distance = ( this->getPoint( 0 ) - this->IS_position_catching ).norm();

    for(  unsigned int bead_id = 2 ; bead_id < this->getNumberOfPoints() ;  bead_id ++ )
    {

	Vector3d bead_position = this->getPoint( bead_id );
	double distance_2 = ( bead_position - this->IS_position_catching ).norm();
	if( distance_2 > distance * 1.1 )
	{
	     std::vector<Vector3d> erased_vectors = this->get_dynein_points_in_IS_and_erase();

	     Vector3d posledni_tangent = this->get_last_Tangent();
	     Vector3d position = this->getPoint( this->getNumberOfPoints() - 2 );
	     position = position + posledni_tangent / posledni_tangent.norm() * sim_of_Cell::resting_distance;
	     this->setPoint( this->getNumberOfPoints() - 1 , position );
	     this->set_lenght_of_tangents();


	     this->dydein_index	= 0;
	     return erased_vectors;
	}
    }

    std::vector<Vector3d> empty_vector;
    return empty_vector;
}



std::vector<Vector3d> Microtubule::get_catching_points_IS_sliding( )
{
    if( this->get_dynein_index() != 9 )
    {
        cout<<"this->get_dynein_index() != 9"<<endl;
        cout<<"std::vector< Vector3d  > Microtubule::get_Dynein_points_and_erase()"<<endl;
        cout<<"ERROR_ID = 6835461681364813"<<endl;
        throw("");
    }

    std::vector< Vector3d  > point_coordinates;
    for( unsigned int pair_index = 0 ; pair_index < this->Dynein_motors_2.size() ; pair_index ++ )
    {
        std::pair < Vector3d , double > pair_tmp = this->Dynein_motors_2.at( pair_index );
        Vector3d anchor_position = std::get<0>( pair_tmp );
        point_coordinates.push_back( anchor_position );
    }
    return point_coordinates;
}




std::vector< Vector3d  > Microtubule::get_Dynein_points_without_erasing()
{

    if( this->get_dynein_index() < 3  )
    {
        cout<<"this->get_dynein_index() < 3"<<endl;
        cout<<"std::vector< Vector3d  > Microtubule::get_Dynein_points_without_erasing()"<<endl;
        cout<<"ERROR_ID = 97964651643515"<<endl;
        throw("");
    }

    std::vector< Vector3d  > point_coordinates;
    for( unsigned int pair_index = 0 ; pair_index < this->Dynein_motors_2.size() ; pair_index ++ )
    {
        std::pair < Vector3d , double > pair_tmp = this->Dynein_motors_2.at( pair_index );
        Vector3d anchor_position = std::get<0>( pair_tmp );
        point_coordinates.push_back( anchor_position );
    }
    return point_coordinates;
}








Microtubule::~Microtubule()
{
	// TODO Auto-generated destructor stub
}


#include "GeneralGeometry.h"

//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS


double distance_two_lines( Vector3d a_vector , Vector3d b_vector , Vector3d c_vector , Vector3d d_vector )
{
	Vector3d q_vector = a_vector - c_vector;
	Vector3d cross_product = b_vector.cross( d_vector );
	double nominator = q_vector.dot( cross_product );
	double distance = abs( nominator ) / cross_product.norm();

	return distance;
}

double distance_two_line_segments( Vector3d a_vector , Vector3d b_vector , Vector3d c_vector , Vector3d d_vector )
{
	double result = 0.0;
	//this function is thoroughly documented on papers in general notes
	//I will define basic substitutions
	//sunstitutions defined when expressing \frac{ \par{ r } }{ \par{ t } }
	Vector3d q_vector = a_vector - c_vector;
	double P_scalar = b_vector.dot( q_vector );
	double B_scalar = b_vector.dot( b_vector );
	double R_scalar = b_vector.dot( d_vector );

	//substitutions defined when expressing \frac{ \par{ r } }{ \par{ v } }
	double S_scalar = d_vector.dot( q_vector );

	double E_scalar = d_vector.dot( d_vector );
	double v = ( R_scalar * P_scalar - B_scalar * S_scalar ) / ( R_scalar * R_scalar - B_scalar * E_scalar );
	double t = ( R_scalar * v - P_scalar ) / B_scalar;

	//if the lines are parallel, v will be nan, because R_scalar * R_scalar - B_scalar * E_scalar == 0

	if( ( ( v > 0 ) && ( v < 1 ) ) && ( ( t > 0 ) && ( t < 1 ) ) )
	{
		//minimum is in the square
		Vector3d point_1 = a_vector + b_vector * t;
		Vector3d point_2 = c_vector + d_vector * v;
		result = ( point_1 - point_2 ).norm();
	}
	else
	{
		MatrixXd result_array = MatrixXd::Zero( 8 , 3 );
		//Corner investigation
		//1: v == 0 ^ t == 0
		double result_corner_1 = ( q_vector ).norm();
		result_array( 0 , 0 ) = result_corner_1;
		result_array( 0 , 1 ) = 0;
		result_array( 0 , 2 ) = 0;

		//2: t == 1 ^ v = 0
		double result_corner_2 = ( a_vector + b_vector - c_vector ).norm();
		result_array( 1 , 0 ) = result_corner_2;
		result_array( 1 , 1 ) = 1.0;
		result_array( 1 , 2 ) = 0.0;

		//3: t == 1 ^ v = 1
		double result_corner_3 = ( a_vector + b_vector - c_vector - d_vector ).norm();
		result_array( 2 , 0 ) = result_corner_3;
		result_array( 2 , 1 ) = 1.0;
		result_array( 2 , 2 ) = 1.0;


		//4: t == 1 ^ v = 1
		double result_corner_4 = ( a_vector - c_vector - d_vector ).norm();
		result_array( 3 , 0 ) = result_corner_4;
		result_array( 3 , 1 ) = 0.0;
		result_array( 3 , 2 ) = 1.0;



		//Investigation of Edges
		v = 0;
		t = ( - 1.0 ) * P_scalar / B_scalar;


		if( ( t > 0 ) && ( t < 1 ) )
		{
			Vector3d point_1_edge_1 = a_vector + b_vector * t;
			Vector3d point_2_edge_1 = c_vector + d_vector * v;
			double result_corner = ( point_1_edge_1 - point_2_edge_1 ).norm();
			result_array( 4 , 0 ) = result_corner;
			result_array( 4 , 1 ) = t;
			result_array( 4 , 2 ) = v;
		}
		else
		{
			result_array( 4 , 0 ) = 1.1 * result_array( 0 , 0 );
		}

		v = 1;
		t = ( R_scalar - P_scalar ) / B_scalar;
		if( ( t > 0 ) && ( t < 1 ) )
		{
			Vector3d point_1_edge_3 = a_vector + b_vector * t;
			Vector3d point_2_edge_3 = c_vector + d_vector * v;
			double result_corner = ( point_1_edge_3 - point_2_edge_3 ).norm();
			result_array( 5 , 0 ) = result_corner;
			result_array( 5 , 1 ) = t;
			result_array( 5 , 2 ) = v;
		}
		else
		{
			result_array( 5 , 0 ) = 1.1 * result_array( 0 , 0 );
		}


		t = 0;
		v = S_scalar / E_scalar;
		if( ( v > 0 ) && ( v < 1 ) )
		{
			Vector3d point_1_edge_2 = a_vector + b_vector * t;
			Vector3d point_2_edge_2 = c_vector + d_vector * v;
			double result_corner = ( point_1_edge_2 - point_2_edge_2 ).norm();
			result_array( 6 , 0 ) = result_corner;
			result_array( 6 , 1 ) = t;
			result_array( 6 , 2 ) = v;
		}
		else
		{
			result_array( 6 , 0 ) = 1.1 * result_array( 0 , 0 );
		}

		t = 1;
		v = ( S_scalar + R_scalar ) / E_scalar;
		if( ( v > 0 ) && ( v < 1 ) )
		{
	        Vector3d point_1_edge_4 = a_vector;
	        point_1_edge_4 = point_1_edge_4 + ( b_vector * t );
	        Vector3d point_2_edge_4 = c_vector;
			point_2_edge_4 = point_2_edge_4 + ( d_vector * v );
			double result_corner = ( point_1_edge_4 - point_2_edge_4 ).norm();
			result_array( 7 , 0 ) = result_corner;
			result_array( 7 , 1 ) = t;
			result_array( 7 , 2 ) = v;
		}
		else
		{
			result_array( 7 , 0 ) = 1.1 * result_array( 0 , 0 );
		}

		double tmp = result_array( 0 , 0 );
		for( unsigned int i = 1 ; i < result_array.rows() ; i ++ )
		{
			if( result_array( i , 0 ) < tmp )
			{
				tmp = result_array( i , 0 );
			}
		}

		result = tmp;
	}
	return result;
}



double distance_two_line_segments( Vector3d a_vector , Vector3d b_vector , Vector3d c_vector , Vector3d d_vector , unsigned int interaction_index )
{
	double result = 0.0;
	//this function is thoroughly documented on papers in general notes
	//I will define basic substitutions
	//sunstitutions defined when expressing \frac{ \par{ r } }{ \par{ t } }
	Vector3d q_vector = a_vector - c_vector;
	double P_scalar = b_vector.dot( q_vector );
	double B_scalar = b_vector.dot( b_vector );
	double R_scalar = b_vector.dot( d_vector );

	//substitutions defined when expressing \frac{ \par{ r } }{ \par{ v } }
	double S_scalar = d_vector.dot( q_vector );


	double E_scalar = d_vector.dot( d_vector );

	double v = ( R_scalar * P_scalar - B_scalar * S_scalar ) / ( R_scalar * R_scalar - B_scalar * E_scalar );
	double t = ( R_scalar * v - P_scalar ) / B_scalar;




	if( ( ( v > 0 ) && ( v < 1 ) ) && ( ( t > 0 ) && ( t < 1 ) ) )
	{
		//minimum is in the square
		Vector3d point_1 = a_vector + b_vector * t;
		Vector3d point_2 = c_vector + d_vector * v;
		result = ( point_1 - point_2 ).norm();
		interaction_index = 0;
	}
	else
	{
		MatrixXd result_array = MatrixXd::Zero( 8 , 1 );


		//Corner investigation
		//1: v == 0 ^ t == 0
		double result_corner_1 = ( q_vector ).norm();
		result_array( 0 , 0 ) = result_corner_1;
		//2: t == 1 ^ v = 0
		double result_corner_2 = ( a_vector + b_vector - c_vector ).norm();
		result_array( 1 , 0 ) = result_corner_2;
		//3: t == 1 ^ v = 1
		double result_corner_3 = ( a_vector + b_vector - c_vector - d_vector ).norm();
		result_array( 2 , 0 ) = result_corner_3;
		//4: t == 1 ^ v = 1
		double result_corner_4 = ( a_vector - c_vector - d_vector ).norm();
		result_array( 3 , 0 ) = result_corner_4;


		//Investigation of Edges
		// Edge 1: v == 0
		v = 0;
		t = ( - 1.0 ) * P_scalar / B_scalar;
		if( ( t > 0 ) && ( t < 1 ) )
		{
			Vector3d point_1_edge_1 = a_vector + b_vector * t;
			Vector3d point_2_edge_1 = c_vector + d_vector * v;
			double result_corner = ( point_1_edge_1 - point_2_edge_1 ).norm();
			result_array( 4 , 0 ) = result_corner;
		}
		else
		{
			result_array( 4 , 0 ) = 1.1 * result_array( 0 , 0 );
		}





		// Edge 3: v == 0
		v = 1;
		t = ( R_scalar - P_scalar ) / B_scalar;
		if( ( t > 0 ) && ( t < 1 ) )
		{
			Vector3d point_1_edge_3 = a_vector + b_vector * t;
			Vector3d point_2_edge_3 = c_vector + d_vector * v;
			double result_corner = ( point_1_edge_3 - point_2_edge_3 ).norm();
			result_array( 5 , 0 ) = result_corner;
		}
		else
		{
			result_array( 5 , 0 ) = 1.1 * result_array( 0 , 0 );
		}




		// Edge 2: v == 0
		t = 0;
		v = S_scalar / E_scalar;
		if( ( v > 0 ) && ( v < 1 ) )
		{
			Vector3d point_1_edge_2 = a_vector + b_vector * t;
			Vector3d point_2_edge_2 = c_vector + d_vector * v;
			double result_corner = ( point_1_edge_2 - point_2_edge_2 ).norm();
			result_array( 6 , 0 ) = result_corner;
		}
		else
		{
			result_array( 6 , 0 ) = 1.1 * result_array( 0 , 0 );
		}




		// Edge 4: v == 0
		t = 1;
		v = ( S_scalar + R_scalar ) / E_scalar;
		if( ( v > 0 ) && ( v < 1 ) )
		{
	        Vector3d point_1_edge_4 = a_vector;
	        point_1_edge_4 = point_1_edge_4 + ( b_vector * t );
	        Vector3d point_2_edge_4 = c_vector;
			point_2_edge_4 = point_2_edge_4 + ( d_vector * v );
			double result_corner = ( point_1_edge_4 - point_2_edge_4 ).norm();
			result_array( 7 , 0 ) = result_corner;
		}
		else
		{
			result_array( 7 , 0 ) = 1.1 * result_array( 0 , 0 );
		}

		double tmp = result_array( 0 , 0 );

		for( unsigned int i = 1 ; i < result_array.rows() ; i ++ )
		{
			if( result_array( i , 0 ) < tmp )
			{
				interaction_index = i;
				tmp = result_array( i , 0 );

			}
		}

		result = tmp;
	}
	return result;
}


double distance_two_line_segments( Vector3d a_vector , Vector3d b_vector , Vector3d c_vector , Vector3d d_vector , Vector3d& p_1 , Vector3d& p_2 , double& t_return, double& v_return )
{
        //this function is the same as the preceding one, but it returns the points and parametres 
    
    
	double result = 0.0;
	//this function is thoroughly documented on papers in general notes
	//I will define basic substitutions
	//sunstitutions defined when expressing \frac{ \par{ r } }{ \par{ t } }
	Vector3d q_vector = a_vector - c_vector;
	double P_scalar = b_vector.dot( q_vector );
	double B_scalar = b_vector.dot( b_vector );
	double R_scalar = b_vector.dot( d_vector );

	//substitutions defined when expressing \frac{ \par{ r } }{ \par{ v } }
	double S_scalar = d_vector.dot( q_vector );

	double E_scalar = d_vector.dot( d_vector );
	double v = ( R_scalar * P_scalar - B_scalar * S_scalar ) / ( R_scalar * R_scalar - B_scalar * E_scalar );
	double t = ( R_scalar * v - P_scalar ) / B_scalar;

	//if the lines are parallel, v will be nan, because R_scalar * R_scalar - B_scalar * E_scalar == 0


	if( ( ( v > 0 ) && ( v < 1 ) ) && ( ( t > 0 ) && ( t < 1 ) ) )
	{
		//minimum is in the square
		Vector3d point_1 = a_vector + b_vector * t;
		Vector3d point_2 = c_vector + d_vector * v;
		result = ( point_1 - point_2 ).norm();
		p_1 = point_1;
		p_2 = point_2;
		t_return = t;
		v_return = v;
	}
	else
	{
		MatrixXd result_array = MatrixXd::Zero( 8 , 3 );
		//Corner investigation
		//1: v == 0 ^ t == 0
		double result_corner_1 = ( q_vector ).norm();
		result_array( 0 , 0 ) = result_corner_1;
		result_array( 0 , 1 ) = 0;
		result_array( 0 , 2 ) = 0;

		//2: t == 1 ^ v = 0
		double result_corner_2 = ( a_vector + b_vector - c_vector ).norm();
		result_array( 1 , 0 ) = result_corner_2;
		result_array( 1 , 1 ) = 1.0;
		result_array( 1 , 2 ) = 0.0;

		//3: t == 1 ^ v = 1
		double result_corner_3 = ( a_vector + b_vector - c_vector - d_vector ).norm();
		result_array( 2 , 0 ) = result_corner_3;
		result_array( 2 , 1 ) = 1.0;
		result_array( 2 , 2 ) = 1.0;


		//4: t == 1 ^ v = 1
		double result_corner_4 = ( a_vector - c_vector - d_vector ).norm();
		result_array( 3 , 0 ) = result_corner_4;
		result_array( 3 , 1 ) = 0.0;
		result_array( 3 , 2 ) = 1.0;



		//Investigation of Edges
		// Edge 1: v == 0
		v = 0;
		t = ( - 1.0 ) * P_scalar / B_scalar;

		if( ( t > 0 ) && ( t < 1 ) )
		{
			Vector3d point_1_edge_1 = a_vector + b_vector * t;
			Vector3d point_2_edge_1 = c_vector + d_vector * v;
			double result_corner = ( point_1_edge_1 - point_2_edge_1 ).norm();
			result_array( 4 , 0 ) = result_corner;
			result_array( 4 , 1 ) = t;
			result_array( 4 , 2 ) = v;
		}
		else
		{
			result_array( 4 , 0 ) = 1.1 * result_array( 0 , 0 );
		}


		// Edge 3: v == 0
		v = 1;
		t = ( R_scalar - P_scalar ) / B_scalar;
		if( ( t > 0 ) && ( t < 1 ) )
		{
			Vector3d point_1_edge_3 = a_vector + b_vector * t;
			Vector3d point_2_edge_3 = c_vector + d_vector * v;
			double result_corner = ( point_1_edge_3 - point_2_edge_3 ).norm();
			result_array( 5 , 0 ) = result_corner;
			result_array( 5 , 1 ) = t;
			result_array( 5 , 2 ) = v;
		}
		else
		{
			result_array( 5 , 0 ) = 1.1 * result_array( 0 , 0 );
		}




		// Edge 2: v == 0
		t = 0;
		v = S_scalar / E_scalar;
		if( ( v > 0 ) && ( v < 1 ) )
		{
			Vector3d point_1_edge_2 = a_vector + b_vector * t;
			Vector3d point_2_edge_2 = c_vector + d_vector * v;
			double result_corner = ( point_1_edge_2 - point_2_edge_2 ).norm();
			result_array( 6 , 0 ) = result_corner;
			result_array( 6 , 1 ) = t;
			result_array( 6 , 2 ) = v;
		}
		else
		{
			result_array( 6 , 0 ) = 1.1 * result_array( 0 , 0 );
		}

		t = 1;
		v = ( S_scalar + R_scalar ) / E_scalar;
		if( ( v > 0 ) && ( v < 1 ) )
		{
	        Vector3d point_1_edge_4 = a_vector;
	        point_1_edge_4 = point_1_edge_4 + ( b_vector * t );
	        Vector3d point_2_edge_4 = c_vector;
			point_2_edge_4 = point_2_edge_4 + ( d_vector * v );
			double result_corner = ( point_1_edge_4 - point_2_edge_4 ).norm();
			result_array( 7 , 0 ) = result_corner;
			result_array( 7 , 1 ) = t;
			result_array( 7 , 2 ) = v;
		}
		else
		{
			result_array( 7 , 0 ) = 1.1 * result_array( 0 , 0 );
		}

		double tmp = result_array( 0 , 0 );
		for( unsigned int i = 1 ; i < result_array.rows() ; i ++ )
		{
			if( result_array( i , 0 ) < tmp )
			{
				tmp = result_array( i , 0 );
				t_return = result_array( i , 1 );
				v_return = result_array( i , 2 );
				p_1 = a_vector + t_return * b_vector;
				p_2 = c_vector + v_return * d_vector;
			}
		}

		result = tmp;
	}
	return result;
}


double coordinate_norm( MatrixXd coordinates )
{
	if( ( coordinates.rows() % 3 != 0 ) || ( coordinates.cols() != 1 ) )
	{
                cout<<"double coordinate_norm( MatrixXd coordinates )"<<endl;
		cout<<" ( coordinates.rows() % 3 != 0 ) || ( coordinates.colums() != 1 )"<<endl;
		cout<<" coordinates.rows() = "<< coordinates.rows()<<endl;
                cout<<" coordinates.cols = "<< coordinates.cols()<<endl;
		throw( "" );
	}

	double cumulative = 0;

	for( unsigned int i = 0 ; i < coordinates.rows() / 3 ; i ++ )
	{
		Vector3d point( 0.0 , 0.0 , 0.0 );
		for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
		{
			point( dimension ) = coordinates( 3 * i + dimension , 0 );
		}
		cumulative = cumulative + point.norm();
	}

	return cumulative;
}





double distance_point_segment( Vector3d a_vector , Vector3d b_vector , Vector3d y_vector )
{
	//segment: \vec{ x } = \vec{ a } + t * \vec{ b }
	//the result can be computed using easy differential equation
	double result = 0.0;
	double AAA = b_vector.dot( y_vector );
	double BBB = b_vector.dot( a_vector );
	double CCC = b_vector.norm() * b_vector.norm();
	double t = ( AAA - BBB ) / CCC;

	if( t > 1 )
	{
		Vector3d closest_point = a_vector + b_vector;
		double distance = ( y_vector - closest_point ).norm();
		return distance;
	}
	else if( t < 0 )
	{
		Vector3d closest_point = a_vector;
		double distance = ( y_vector - closest_point ).norm();
		return distance;
	}
	else
	{
		Vector3d closest_point = a_vector + t * b_vector;
		double distance = ( y_vector - closest_point ).norm();
		return distance;
	}

	return result;
}


double distance_point_segment( Vector3d a_vector , Vector3d b_vector , Vector3d y_vector , Vector3d& closest_point_of_segment )
{
	//segment: \vec{ x } = \vec{ a } + t * \vec{ b }
	//the result can be computed using easy differential equation
	double result = 0.0;
	double AAA = b_vector.dot( y_vector );
	double BBB = b_vector.dot( a_vector );
	double CCC = b_vector.norm() * b_vector.norm();
	double t = ( AAA - BBB ) / CCC;
    Vector3d closest_point;
	if( t > 1 )
	{
		closest_point = a_vector + b_vector;
		closest_point_of_segment = closest_point;
	}
	else if( t < 0 )
	{
		closest_point = a_vector;		
		closest_point_of_segment = closest_point;
	}
	else
	{
		closest_point = a_vector + t * b_vector;	
		closest_point_of_segment = closest_point;

	}
    double distance = ( y_vector - closest_point ).norm();
	return distance;

}

double distance_plane_point( Vector3d plane_axis , Vector3d point_on_plane , Vector3d point )
{
				double D = ( - 1.0 ) * plane_axis.dot( point_on_plane );
				double nominator = plane_axis.dot( point ) + D;
				double denominator = plane_axis.norm();
				double distance = nominator / denominator;
				return distance;
}

double distance_point_line( Vector3d first_line_point , Vector3d second_line_point , Vector3d y_vector )
{
    Vector3d first_line = y_vector - first_line_point;
    Vector3d second_line = y_vector - second_line_point;
    double numerator = ( first_line.cross( second_line ) ).norm();
    double denominator = ( second_line_point - first_line_point ).norm();
    
    return numerator / denominator;
}


Vector3d project_point_on_surface_of_elipse( Vector3d position )
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




Vector3d projection_of_tangent_on_plane( Vector3d plane_axis , Vector3d tangent )
{
    //vector plane axis - axis vector of the plane
    //tangent that will be projeted
    Vector3d projection = tangent - tangent.dot( plane_axis ) / ( plane_axis.norm() * plane_axis.norm() ) * plane_axis;
    return projection;
}

Vector3d projection_of_point_on_plane( Vector3d plane_axis , Vector3d point_on_plane , Vector3d point_to_be_projected )
{
    Vector3d plane_axis_tmp = plane_axis / plane_axis.norm();
    Vector3d projected_point = point_to_be_projected - plane_axis_tmp.dot( point_to_be_projected - point_on_plane ) * plane_axis_tmp;
    return projected_point;
}


double magnitude_3D_matrix( MatrixXd& matrix)
{
    if( ( matrix.rows() % 3 != 0 ) && ( matrix.cols() != 1 ) )
    {
        cout<<"( matrix.rows() % 3 != 0 ) && ( matrix.cols() != 1 ) in double magnitude_3D_matrix( MatrixXd& matrix)"<<endl;
        cout<<"General Geometry ERROR_ID = 9646451574794646841"<<endl;
        throw("");
    }
    
    unsigned int number_of_vectors = matrix.rows() / 3;
    double overall_force = 0;
    for( unsigned int number_bead = 0 ; number_bead < number_of_vectors ; number_bead++ )
    {
        Vector3d force( 0 , 0 , 0 );
        for( unsigned int i = 0 ; i < 3 ; i ++ )
        {
            force( i ) = matrix( 3 * number_bead + i , 0 );
        }
        overall_force = overall_force + force.norm();
    }
    
    return overall_force;
}











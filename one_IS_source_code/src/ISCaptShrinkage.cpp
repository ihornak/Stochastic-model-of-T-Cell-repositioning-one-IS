//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS

#include "ISCaptShrinkage.h"

IS_Capt_Shrinkage::IS_Capt_Shrinkage() {
	this->center_of_IS_front = Vector3d( 0.0 , 0.0 , 0.0 );
	this->center_of_IS_rear = Vector3d( -1.0 , 0.0 , 0.0 );
	this->radius_of_IS = 1.0;
	this->axis_of_IS = Vector3d( 1.0 , 0.0 , 0.0 );

}

IS_Capt_Shrinkage::IS_Capt_Shrinkage( Vector3d center_of_IS_front , double radius_argument , Vector3d center_of_IS_rear )
{
	if( ( center_of_IS_front - center_of_IS_rear ).norm() == 0 )
	{
		cout<<"( center_of_IS_front - center_of_IS_rear ).norm() == 0 "<<endl;
		cout<<" ERROR_ID = 5978134268 "<<endl;
		throw("");
	}

	this->axis_of_IS = ( center_of_IS_front - center_of_IS_rear ) / ( center_of_IS_front - center_of_IS_rear ).norm();
	this->center_of_IS_front = center_of_IS_front;
	this->center_of_IS_rear = center_of_IS_rear;
	this->radius_of_IS = radius_argument;

}

IS_Capt_Shrinkage::IS_Capt_Shrinkage( Vector3d center_of_IS_front , double radius_argument , Vector3d center_of_IS_rear , Vector3d trap_arg )
{
	if( ( center_of_IS_front - center_of_IS_rear ).norm() == 0 )
	{
		cout<<"( center_of_IS_front - center_of_IS_rear ).norm() == 0 "<<endl;
		cout<<" ERROR_ID = 5978134268 "<<endl;
		throw("");
	}

	this->axis_of_IS = ( center_of_IS_front - center_of_IS_rear ) / ( center_of_IS_front - center_of_IS_rear ).norm();

	// TODO Auto-generated constructor stub
	this->center_of_IS_front = center_of_IS_front;
	this->center_of_IS_rear = center_of_IS_rear;
	this->radius_of_IS = radius_argument;

}




IS_Capt_Shrinkage::IS_Capt_Shrinkage( IS_Capt_Shrinkage& tmp )
{
	this->center_of_IS_front = tmp.get_center_of_IS_front();
	this->center_of_IS_rear = tmp.center_of_IS_rear;
	this->radius_of_IS = tmp.get_radius_of_IS();
	this->axis_of_IS = tmp.get_axis_of_IS();

}


IS_Capt_Shrinkage IS_Capt_Shrinkage::operator=( const IS_Capt_Shrinkage& tmp )
{
	if( this == &tmp )
	{
		cout<<" *this == &tmp "<<endl;
		throw("");
	}


	this->center_of_IS_front = tmp.get_center_of_IS_front();
	this->center_of_IS_rear = tmp.center_of_IS_rear;
	this->radius_of_IS = tmp.get_radius_of_IS();
	this->axis_of_IS = tmp.get_axis_of_IS();
	
	return *this;
}



Vector3d IS_Capt_Shrinkage::get_center_of_IS_front() const
{
	return this->center_of_IS_front;
}


Vector3d IS_Capt_Shrinkage::get_center_of_IS_rear() const
{
	return this->center_of_IS_rear;
}



double IS_Capt_Shrinkage::get_radius_of_IS() const
{
	return this->radius_of_IS;
}

Vector3d IS_Capt_Shrinkage::get_axis_of_IS() const
{
	return this->axis_of_IS;
}




void IS_Capt_Shrinkage::set_center_of_IS_front( Vector3d center_argument_front )
{
	this->center_of_IS_front = center_argument_front;
}


void IS_Capt_Shrinkage::set_center_of_IS_rear( Vector3d center_argument_rear )
{
	this->center_of_IS_rear = center_argument_rear;
}




void IS_Capt_Shrinkage::set_radius_of_IS( double radius_center )
{
	this->radius_of_IS = radius_center;
}


void IS_Capt_Shrinkage::set_axis_of_IS( Vector3d axis_argument )
{
	this->axis_of_IS = axis_argument;
}





bool IS_Capt_Shrinkage::check_IS_capture_shrinkage_caught_Main( Vector3d position )
{

	double radius_front = this->get_center_of_IS_front().norm();
	double radius_back = this->get_center_of_IS_rear().norm();

	if( ( position.norm() > radius_front ) &&  ( position.norm() < radius_back ) )
	{
		double distance = distance_point_line( this->get_center_of_IS_front() , this->get_center_of_IS_rear() , position );
		if( distance < this->get_radius_of_IS() )
		{
			return true;
		}
		else
		{
			return false;
		}	
	
	}
	else
	{
		return false;
	}

}



bool IS_Capt_Shrinkage::check_IS_capture_segment_control( Vector3d first_point , Vector3d tangent )
{


    Vector3d c_vector = this->get_center_of_IS_rear();
    Vector3d d_vector = this->get_center_of_IS_front() - this->get_center_of_IS_rear();
    double height_of_IS = d_vector.norm();
    Vector3d p_1( 0.0 , 0.0 , 0.0 );
    Vector3d p_2( 0.0 , 0.0 , 0.0 );
    double t_return = 0;
    double v_return = 0;
    double distance = distance_two_line_segments( first_point , tangent , c_vector , d_vector , p_1 , p_2 , t_return, v_return );
    if( distance <= this->radius_of_IS )
    {
        Vector3d axis_of_IS = this->get_center_of_IS_front() - this->get_center_of_IS_rear();
        axis_of_IS = axis_of_IS / axis_of_IS.norm();
        
        Vector3d point_on_central_plane = c_vector + d_vector / 2.0;
        double distance_plane_segment = distance_plane_point( axis_of_IS , point_on_central_plane , p_2 ); 
        if( abs( distance_plane_segment ) < height_of_IS / 2.0 )
        {
            return true;
        }
        else
        {
            return false;
        }
        
    }
    else
    {
        return false;
    }

}











IS_Capt_Shrinkage::~IS_Capt_Shrinkage() {
	// TODO Auto-generated destructor stub
}

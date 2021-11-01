//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS

#include "Nucleus.h"

Nucleus::Nucleus()
{
	this->a_axis = Nucleus_parametres::A_AXIS;
	this->b_axis = Nucleus_parametres::B_AXIS;
	this->center = Vector3d( 0.0 , 0.0 , 0.0 );
}

Nucleus::Nucleus( Vector3d center_argument )
{
	this->a_axis = Nucleus_parametres::A_AXIS;
	this->b_axis = Nucleus_parametres::B_AXIS;
	this->center = center_argument;
}


Nucleus& Nucleus::operator=( const Nucleus &tmp )
{
	this->a_axis = tmp.get_A_Axis();
	this->b_axis = tmp.get_B_Axis();
	this->center = tmp.center;
	return *this;
}



double Nucleus::get_A_Axis() const
{
	return this->a_axis;
}
double Nucleus::get_B_Axis() const
{
	return this->b_axis;
}


Vector3d Nucleus::force_position_nucleus_wall( Vector3d position )
{
    Vector3d force( 0.0 , 0.0 , 0.0 );
    double B_multiplicator = this->b_axis * this->b_axis;
    double numerator = this->a_axis * this->b_axis;                    
    double A_multiplicator = this->a_axis * this->a_axis;
    position = position - this->get_center();  
    
    double denominator = ( position( 0 ) * position( 0 ) + position( 1 ) * position( 1 ) ) * B_multiplicator;         
    denominator = denominator + position( 2 ) * position( 2 ) * A_multiplicator;
    denominator = sqrt( denominator );
    double c_multiplicator = numerator / denominator;
    if( c_multiplicator <= 1 )
    {
        
    }       
    else
    {
        Vector3d point_of_intersection = c_multiplicator * position;
        double distance_behind_wall = ( position - point_of_intersection ).norm();
        double abs_val_force = Nucleus_parametres::wall_nucleus_k1 * ( exp( Nucleus_parametres::wall_nucleus_k2 * distance_behind_wall ) - 1.0 );
        Vector3d orientation = ( point_of_intersection - position );
        orientation = orientation / orientation.norm();
        force = orientation * abs_val_force;
       
    }
    return force;
}


Vector3d Nucleus::force_position_nucleus_wall_MTOC( Vector3d position )
{
    Vector3d force( 0.0 , 0.0 , 0.0 );
    double B_multiplicator = this->b_axis * this->b_axis;

    double numerator = this->a_axis * this->b_axis;                    
    double A_multiplicator = this->a_axis * this->a_axis;
    position = position - this->get_center();  
    
    double denominator = ( position( 0 ) * position( 0 ) + position( 1 ) * position( 1 ) ) * B_multiplicator;         
    denominator = denominator + position( 2 ) * position( 2 ) * A_multiplicator;
    denominator = sqrt( denominator );
    double c_multiplicator = numerator / denominator;
    if( c_multiplicator <= 1 )
    {
        return force;
    }       
    else
    {
        Vector3d point_of_intersection = c_multiplicator * position;
        double distance_behind_wall = ( position - point_of_intersection ).norm();
        
        double abs_val_force = Nucleus_parametres::wall_nucleus_k1 * ( exp( Nucleus_parametres::wall_nucleus_k2 * distance_behind_wall ) - 1.0 );
        Vector3d orientation = ( point_of_intersection - position ) / ( position - point_of_intersection ).norm();
        force = orientation * abs_val_force;
        return force;
    }
  
}

Vector3d Nucleus::get_center()
{
	return this->center;
}




Nucleus::~Nucleus() {
	// TODO Auto-generated destructor stub
}


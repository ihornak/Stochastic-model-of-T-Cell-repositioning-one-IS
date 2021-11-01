#include "CellShape.h"

//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS

Cell_Shape::Cell_Shape()
{
	// TODO Auto-generated constructor stub
	this->A_Axis = 4.5e-6;
	this->B_Axis = 5.3e-6;
	this->cylinder_width = 1.5e-6;
	this->z_value = ( -1.0 ) * this->B_Axis + this->cylinder_width;

}

Cell_Shape::Cell_Shape( double A_Axis_arg , double B_Axis_arg , double cylinder_width_arg )
{
	this->A_Axis = A_Axis_arg;
	this->B_Axis = B_Axis_arg;
	this->cylinder_width = 1.5e-6;
	this->z_value = this->calculate_z_cylinder();
}

double Cell_Shape::calculate_z_cylinder(  )
{
	if( this->cylinder_width >=  this->A_Axis )
	{
		cout<<" this->cylinder_width >=  this->A_Axis in Cell_Shape::calculate_z_cylinder(  )"	<<endl;
		cout<<" ERROR_ID = 687494596 "<<endl;
		throw("");
	}
	double tmp = 1.0 - ( this->cylinder_width / this->A_Axis ) * ( this->cylinder_width / this->A_Axis );
	double z_value = this->B_Axis * sqrt( tmp );
	return z_value;
}


double Cell_Shape::get_A_Axis()
{
	return this->A_Axis;
}

double Cell_Shape::get_B_Axis()
{
	return this->B_Axis;
}


double Cell_Shape::get_cylinder_width()
{
	return this->cylinder_width;
}

double Cell_Shape::get_z_value()
{
	return this->z_value;
}

















Cell_Shape::~Cell_Shape() {
	// TODO Auto-generated destructor stub
}


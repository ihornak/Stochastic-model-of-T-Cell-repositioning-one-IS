//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS
//Unused - for  future simulations
#include "ideal_mtoc.h"

ideal_MTOC::ideal_MTOC()
{
    this->number_of_points = 1;
    this->number_of_true_mtoc_points = this->set_number_of_true_mtoc_points();
    this->coordinates =  MatrixXd::Zero( 3 , 1 );
    this->original_orientation = Vector3d( 0.0 , 0.0 , 1.0 );
    this->original_orientation = this->original_orientation / this->original_orientation.norm();
    this->orientations_micro = MatrixXd::Zero( 3 * ( this->number_of_points + 1 ) , 1 );
    this->mtoc_coordinates = MatrixXd::Zero( 3 * ( this->number_of_points + 1 ) , 1 );
    this->true_mtoc_coordinates = MatrixXd::Zero( 3 * ( this->number_of_points + 1 ) , 1 );
    this->original_MTOC_center = Vector3d( 0.0 , 0.0 , 1.0 );
    this->original_MTOC_1_point = Vector3d( 0.0 , 0.0 , 1.0 );
}

//TADY POKRACUJ!!!!!!!!!!!!!!!!
ideal_MTOC::ideal_MTOC( unsigned int number_of_MTOC_points , Vector3d original_orientation )
{
    this->number_of_points = number_of_MTOC_points;
    this->number_of_true_mtoc_points = this->set_number_of_true_mtoc_points();
    this->original_orientation = original_orientation;
    this->original_orientation = this->original_orientation / this->original_orientation.norm();
    this->coordinates =  MatrixXd::Zero( 3 * ( this->number_of_points + 1 ) , 1 );
    this->orientations_micro = MatrixXd::Zero( 3 * ( this->number_of_points + 1 ) , 1 );
    this->mtoc_coordinates = MatrixXd::Zero( 3 * ( this->number_of_points + 1 ) , 1 );
    this->true_mtoc_coordinates = MatrixXd::Zero(  3 * ( 21 ) , 1  );
    this->original_MTOC_center = Vector3d( 0.0 , 0.0 , 1.0 );
    this->original_MTOC_1_point = Vector3d( 0.0 , 0.0 , 1.0 );
}



ideal_MTOC::ideal_MTOC(const ideal_MTOC& other)
{
    this->number_of_points = other.number_of_points;
    this->number_of_true_mtoc_points = other.number_of_true_mtoc_points;
    this->number_of_true_mtoc_points = this->set_number_of_true_mtoc_points();
    this->original_orientation = other.original_orientation;
    this->coordinates = other.coordinates;
    this->orientations_micro = other.orientations_micro;
    this->mtoc_coordinates = other.mtoc_coordinates;
}


ideal_MTOC& ideal_MTOC::operator=(const ideal_MTOC& other)
{
    this->number_of_points = other.number_of_points;
    this->number_of_true_mtoc_points = other.number_of_true_mtoc_points;
    this->original_orientation = other.original_orientation;
    this->coordinates = other.coordinates;
    this->orientations_micro = other.orientations_micro;
    this->mtoc_coordinates = other.mtoc_coordinates;
    return *this;
}


unsigned int ideal_MTOC::set_number_of_true_mtoc_points()
{
    return MTOCparam::boundary_MTOC_points;
}


unsigned int ideal_MTOC::get_number_of_true_mtoc_points()
{
    return this->number_of_true_mtoc_points;
}



void ideal_MTOC::set_center(Vector3d point )
{
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        this->coordinates( dimension , 0  ) = point( dimension );
    }
}

void ideal_MTOC::set_mtoc_center(Vector3d point )
{
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        this->mtoc_coordinates( dimension , 0  ) = point( dimension );
    }
}

void ideal_MTOC::set_true_mtoc_center(Vector3d point )
{
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        this->true_mtoc_coordinates( dimension , 0  ) = point( dimension );
    }      
}


void ideal_MTOC::set_original_MTOC_center( Vector3d argument )
{
    this->original_MTOC_center = argument;
}

void ideal_MTOC::set_original_MTOC_1_point( Vector3d argument )
{
    this->original_MTOC_1_point = argument;
}




void ideal_MTOC::set_orientation_center( Vector3d point )
{
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        this->orientations_micro( dimension , 0  ) = point( dimension );
    }  
}



Vector3d ideal_MTOC::get_center( )
{
    Vector3d point( 0.0 , 0.0 , 0.0 );
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        point( dimension ) = this->coordinates( dimension , 0  );
    }
    return point; 
}


Vector3d ideal_MTOC::get_orientation_center( )
{
    Vector3d point( 0.0 , 0.0 , 0.0 );
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        point( dimension ) = this->coordinates( dimension , 0  );
    }
    return point; 
}

Vector3d ideal_MTOC::get_moc_center( )
{
    Vector3d point( 0.0 , 0.0 , 0.0 );
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        point( dimension ) = this->mtoc_coordinates( dimension , 0  );
    }
    return point;     
    
}

Vector3d ideal_MTOC::get_true_moc_center( )
{
    Vector3d point( 0.0 , 0.0 , 0.0 );
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        point( dimension ) = this->true_mtoc_coordinates( dimension , 0  );
    }
    return point;    
}

void ideal_MTOC::set_point( unsigned int point_index , Vector3d point )
{
    if( point_index >= this->number_of_points )
    {
        cout<<"point_index >= this->number_of_points"<<endl;
        cout<<"void ideal_MTOC::set_point( unsigned int point_index , Vector3d point )"<<endl;
        cout<<"ERROR_ID = 8976469464684346854"<<endl;
        throw("");
    }

    

    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        this->coordinates( 3 * ( point_index + 1 ) + dimension , 0  ) = point( dimension );
    }
    
}


void ideal_MTOC::set_orientation( unsigned int point_index , Vector3d point )
{
    
    if( point_index >= this->number_of_points )
    {
        cout<<"point_index >= this->number_of_points"<<endl;
        cout<<"void ideal_MTOC::set_orientation( unsigned int point_index , Vector3d point )"<<endl;
        cout<<"ERROR_ID = 312648472524"<<endl;
        throw("");
    }

    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        this->orientations_micro( 3 * ( point_index + 1 ) + dimension , 0  ) = point( dimension );
    }
}


void ideal_MTOC::set_mtoc_point( unsigned int point_index , Vector3d point )
{
    if( point_index >= this->number_of_points )
    {
        cout<<"point_index >= this->number_of_points"<<endl;
        cout<<"void ideal_MTOC::set_mtoc_point( unsigned int point_index , Vector3d point )"<<endl;
        cout<<"ERROR_ID = 8791251691358"<<endl;
        throw("");
    }

    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        this->mtoc_coordinates( 3 * ( point_index + 1 ) + dimension , 0  ) = point( dimension );
    }    
    
}

void ideal_MTOC::set_true_mtoc_point( unsigned int point_index , Vector3d point )
{
    if( point_index >= this->number_of_true_mtoc_points )
    {
        cout<<"point_index >= this->number_of_points"<<endl;
        cout<<"void ideal_MTOC::set_true_mtoc_point( unsigned int point_index , Vector3d point )"<<endl;
        cout<<"ERROR_ID = 369178452"<<endl;
        throw("");
    }

    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        this->true_mtoc_coordinates( 3 * ( point_index + 1 ) + dimension , 0  ) = point( dimension );
    } 
    
}


Vector3d ideal_MTOC::get_point( unsigned int point_index )
{
    if( point_index >= this->number_of_points )
    {
        cout<<"point_index >= this->number_of_points"<<endl;
        cout<<"void ideal_MTOC::get_point( unsigned int point_index , Vector3d point )"<<endl;
        cout<<"ERROR_ID = 8798964364888"<<endl;
        throw("");
    }    
    
    Vector3d point( 0.0 , 0.0 , 0.0 );
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        point( dimension ) = this->coordinates( 3 * ( point_index + 1 ) + dimension , 0  );
    }
    return point;    
}


Vector3d ideal_MTOC::get_orientation( unsigned int point_index )
{
    
    
    if( point_index >= this->number_of_points )
    {
        cout<<"point_index >= this->number_of_points"<<endl;
        cout<<"void ideal_MTOC::get_orientation( unsigned int point_index , Vector3d point )"<<endl;
        cout<<"ERROR_ID = 8791324634126918431"<<endl;
        throw("");
    } 
    
    Vector3d point( 0.0 , 0.0 , 0.0 );
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        point( dimension ) = this->orientations_micro( 3 * ( point_index + 1 ) + dimension , 0  );
    }
    return point; 
    
}


Vector3d ideal_MTOC::get_mtoc_point( unsigned int point_index )
{
    if( point_index >= this->number_of_points )
    {
        cout<<"point_index >= this->number_of_points"<<endl;
        cout<<"void ideal_MTOC::get_mtoc_point( unsigned int point_index , Vector3d point )"<<endl;
        cout<<"ERROR_ID = 75891813445984"<<endl;
        throw("");
    }    
    
    Vector3d point( 0.0 , 0.0 , 0.0 );
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        point( dimension ) = this->mtoc_coordinates( 3 * ( point_index + 1 ) + dimension , 0  );
    }
    return point;      
    
}

Vector3d ideal_MTOC::get_true_mtoc_point( unsigned int point_index )
{
    if( point_index >= this->get_number_of_true_mtoc_points() )
    {
        cout<<"point_index >= this->number_of_points"<<endl;
        cout<<"void ideal_MTOC::get_true_mtoc_point( unsigned int point_index , Vector3d point )"<<endl;
        cout<<"ERROR_ID = 84916548544844"<<endl;
        throw("");
    }    
    
    Vector3d point( 0.0 , 0.0 , 0.0 );
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        point( dimension ) = this->true_mtoc_coordinates( 3 * ( point_index + 1 ) + dimension , 0  );
    }
    return point;     
    
}


Vector3d ideal_MTOC::get_tangent( unsigned int tangent_index )
{
    //remember that the tangent will return point - center
    if( tangent_index >= this->number_of_points )
    {
        cout<<" tangent_index >= this->number_of_points "<<endl;
        cout<<"void ideal_MTOC::get_tangent( unsigned int point_index , Vector3d point )"<<endl;
        cout<<"ERROR_ID = 8798964364888"<<endl;
        throw("");
    }
    
    Vector3d point = this->get_point( tangent_index );
    Vector3d center = this->get_center( );
    Vector3d tangent = center - point;  
    return tangent;   
}

Vector3d ideal_MTOC::get_orientation_tangent( unsigned int tangent_index )
{
    //remember that the tangent will return point - center
    if( tangent_index >= this->number_of_points )
    {
        cout<<" tangent_index >= this->number_of_points "<<endl;
        cout<<"void ideal_MTOC::get_orientation_tangent( unsigned int point_index , Vector3d point )"<<endl;
        cout<<"ERROR_ID = 145594683544"<<endl;
        throw("");
    }
    
    Vector3d point = this->get_orientation( tangent_index );
    Vector3d center = this->get_orientation_center( );
    Vector3d tangent = point - center;  
    return tangent; 
}




MatrixXd ideal_MTOC::get_coordinates(  )
{
    return this->coordinates;    
}

MatrixXd ideal_MTOC::get_orientations()
{
    return this->orientations_micro;    
}

MatrixXd ideal_MTOC::get_mtoc_coordinates()
{
    return this->mtoc_coordinates;    
}

MatrixXd ideal_MTOC::get_true_mtoc_coordinates()
{
    return this->true_mtoc_coordinates;
}




void ideal_MTOC::set_coordinates(  MatrixXd coordinates_tmp  )
{
    if( coordinates_tmp.rows() != 3 * ( this->number_of_points + 1 ) )
    {
        cout<<"coordinates_tmp.rows() != 3 * ( this->number_of_points + 1 )"<<endl;
        cout<<"void ideal_MTOC::set_coordinates(  MatrixXd coordinates_tmp  )"<<endl;
        cout<<"ERROR_ID = 78941346489465"<<endl;
        throw("");        
    }
    this->coordinates = coordinates_tmp;
}


void ideal_MTOC::set_orientations(  MatrixXd coordinates_tmp  )
{
    if( coordinates_tmp.rows() != 3 * ( this->number_of_points + 1 ) )
    {
        cout<<"coordinates_tmp.rows() != 3 * ( this->number_of_points + 1 )"<<endl;
        cout<<"void ideal_MTOC::set_coordinates(  MatrixXd coordinates_tmp  )"<<endl;
        cout<<"ERROR_ID = 78941346489465"<<endl;
        throw("");        
    }
    this->orientations_micro = coordinates_tmp;
}

void ideal_MTOC::set_mtoc_coordinates(  MatrixXd coordinates_tmp )
{
    if( coordinates_tmp.rows() != 3 * ( this->number_of_points + 1 ) )
    {
        cout<<"coordinates_tmp.rows() != 3 * ( this->number_of_points + 1 )"<<endl;
        cout<<"void ideal_MTOC::set_mtoc_orientations(  MatrixXd coordinates_tmp  )"<<endl;
        cout<<"ERROR_ID = 326164648316489"<<endl;
        throw("");        
    }
    this->mtoc_coordinates = coordinates_tmp;  
        
}


void ideal_MTOC::set_true_mtoc_coordinates(  MatrixXd coordinates_tmp )
{
    if( coordinates_tmp.rows() != 3 * ( this->number_of_true_mtoc_points + 1 ) )
    {
        cout<<"coordinates_tmp.rows() != 3 * ( this->number_of_points + 1 )"<<endl;
        cout<<"void ideal_MTOC::set_true_mtoc_coordinates(  MatrixXd coordinates_tmp  )"<<endl;
        cout<<"ERROR_ID = 164568416431564"<<endl;
        throw("");        
    }
    this->true_mtoc_coordinates = coordinates_tmp;
    
}






void ideal_MTOC::rotate_MTOC( Vector3d center_real_MTOC , MatrixXd& coordinates_arg )
{
    if( this->number_of_points != ( coordinates_arg.rows() / 3 - 1 ) )
    {
        cout<<"this->number_of_points != ( coordinates_arg.rows() / 3 - 1 )"<<endl;
        cout<<" void ideal_MTOC::rotate_MTOC( Vector3d center_real_MTOC , MatrixXd& coordinates_arg ) "<<endl;
        cout<<"ERROR_ID = 987733412648315"<<endl;
        throw("");
    }
        
    //construction of quaterniaon
    Vector3d center = this->get_center(); 
    double cosinus_angle = center_real_MTOC.dot( center ) / ( center.norm() * center_real_MTOC.norm() );
    double rotation_angle = acos( cosinus_angle );
    if( rotation_angle < 1e-4 )
    {
        coordinates_arg = this->coordinates;
        return;
    }
    else
    {
     
        Vector3d axis_rotation = center.cross( center_real_MTOC );
        axis_rotation = axis_rotation / axis_rotation.norm();
        Quaternion<double> q;
        q = AngleAxis<double>( rotation_angle , axis_rotation );
        Vector3d new_center =  q * center;
        
    
        for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
        {
            coordinates_arg( dimension , 0 ) = new_center( dimension );
        }

    
        for( unsigned int point = 0 ; point  < this->number_of_points ; point ++ )
        {
            Vector3d point_tmp = this->get_mtoc_point( point );
            Vector3d new_point = q * point_tmp;
            //cout<<new_point - point_tmp<<endl;
            for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
            {
                coordinates_arg( 3 * ( point + 1 ) + dimension , 0 ) = new_point( dimension );
            }           
        }
    }
}



void ideal_MTOC::rotate_orientations( Vector3d center_real_MTOC , MatrixXd& coordinates_arg )
{
    if( this->number_of_points != ( coordinates_arg.rows() / 3 - 1 ) )
    {
        cout<<"this->number_of_points != ( coordinates_arg.rows() / 3 - 1 )"<<endl;
        cout<<" void ideal_MTOC::rotate_MTOC( Vector3d center_real_MTOC , MatrixXd& coordinates_arg ) "<<endl;
        cout<<"ERROR_ID = 987733412648315"<<endl;
        throw("");
    }
        
    //construction of quaterniaon
    Vector3d center = this->get_center(); 
    double cosinus_angle = center_real_MTOC.dot( center ) / ( center.norm() * center_real_MTOC.norm() );
    double rotation_angle = acos( cosinus_angle );
    if( rotation_angle < 1e-4 )
    {
        coordinates_arg = this->orientations_micro;
        return;
    }
    else
    {
     
        Vector3d axis_rotation = center.cross( center_real_MTOC );
        axis_rotation = axis_rotation / axis_rotation.norm();
        Quaternion<double> q;
        q = AngleAxis<double>( rotation_angle , axis_rotation );
        Vector3d new_center =  q * center;
        
    
        for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
        {
            coordinates_arg( dimension , 0 ) = new_center( dimension );
        }

    
        for( unsigned int point = 0 ; point  < this->number_of_points ; point ++ )
        {
            Vector3d point_tmp = this->get_orientation( point );
            Vector3d new_point = q * point_tmp;
            for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
            {
                coordinates_arg( 3 * ( point + 1 ) + dimension , 0 ) = new_point( dimension );
            }           
        }
    }
    
    
}


void ideal_MTOC::rotate_mtoc_coordinates( Vector3d center_real_MTOC , MatrixXd& coordinates_arg )
{
    if( this->number_of_points != ( coordinates_arg.rows() / 3 - 1 ) )
    {
        cout<<"this->number_of_points != ( coordinates_arg.rows() / 3 - 1 )"<<endl;
        cout<<" void ideal_MTOC::rotate_mtoc_coordinates( Vector3d center_real_MTOC , MatrixXd& coordinates_arg ) "<<endl;
        cout<<"ERROR_ID = 8761489164"<<endl;
        throw("");
    }
        
    //construction of quaterniaon
    Vector3d center = this->get_moc_center(); 
    double cosinus_angle = center_real_MTOC.dot( center ) / ( center.norm() * center_real_MTOC.norm() );
    double rotation_angle = acos( cosinus_angle );
    if( rotation_angle < 1e-4 )
    {
        coordinates_arg = this->mtoc_coordinates;
        return;
    }
    else
    {
     
        Vector3d axis_rotation = center.cross( center_real_MTOC );
        axis_rotation = axis_rotation / axis_rotation.norm();
        Quaternion<double> q;
        q = AngleAxis<double>( rotation_angle , axis_rotation );
        Vector3d new_center =  q * center;
        
    
        for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
        {
            coordinates_arg( dimension , 0 ) = new_center( dimension );
        }

    
        for( unsigned int point = 0 ; point  < this->number_of_points ; point ++ )
        {
            Vector3d point_tmp = this->get_mtoc_point( point );
            Vector3d new_point = q * point_tmp;
            for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
            {
                coordinates_arg( 3 * ( point + 1 ) + dimension , 0 ) = new_point( dimension );
            }           
        }
    }
    
    
}

MatrixXd ideal_MTOC::translate_true_mtoc_coordinates_to_radius( double radius_arg )
{
    if( radius_arg > Cell_parametres::A_AXIS + 1e-6 )
    {       
        cout<<"radius_arg > Cell_parametres::A_AXIS"<<endl;
        cout<<"void ideal_MTOC::translate_true_mtoc_coordinates_to_radius( double radius_arg )"<<endl;
        cout<<"ERROR_ID = 5445454354"<<endl;
        throw("");
    }
    
    MatrixXd translated = MatrixXd::Zero( 3 * ( this->get_number_of_true_mtoc_points() + 1 ) , 1 );
    Vector3d center_new = this->get_true_moc_center();
    center_new( 2 ) = radius_arg;

    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        translated( dimension , 0 ) = center_new( dimension );
    }
    
    for( unsigned int point_index = 0 ; point_index < this->get_number_of_true_mtoc_points() ; point_index ++ )
    {
        Vector3d point = this->get_true_mtoc_point( point_index );  
        //cout<<"point"<<endl;
        //cout<<point<<endl;
        point( 2 ) = radius_arg;
        for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
        {
            translated( 3 * ( point_index + 1 ) + dimension , 0 ) = point( dimension );
        }
    }
    //cout<<translated<<endl;
    return translated;
}





void ideal_MTOC::rotate_and_translate_true_mtoc_coordinates( Vector3d center_real_MTOC , MatrixXd& coordinates_arg )
{
    if( coordinates_arg.rows() != 3 * ( this->number_of_true_mtoc_points + 1 ) )
    {
        cout<<"this->number_of_true_mtoc_points != ( coordinates_arg.rows() / 3 - 1 )"<<endl;
        cout<<" void ideal_MTOC::rotate_and_translate_true_mtoc_coordinates( Vector3d center_real_MTOC , MatrixXd& coordinates_arg ) "<<endl;
        cout<<"ERROR_ID = 3126456873164"<<endl;
        throw("");
    }
    //construction of quaterniaon
    Vector3d center = this->get_true_moc_center(); 
    double cosinus_angle = center_real_MTOC.dot( center ) / ( center.norm() * center_real_MTOC.norm() );
    double rotation_angle = acos( cosinus_angle );
    Vector3d axis_rotation = center.cross( center_real_MTOC );
    //cout<<"rotation_angle = "<<rotation_angle<<endl;
    //cout<<axis_rotation / axis_rotation.norm()<<endl;
    
    
    if( axis_rotation.norm() == 0 )
    {
        coordinates_arg = this->true_mtoc_coordinates;
        double radius_of_mtoc = center_real_MTOC.norm();
        //cout<<"radius_of_mtoc = "<<radius_of_mtoc<<endl; 
        MatrixXd new_matrix = translate_true_mtoc_coordinates_to_radius( radius_of_mtoc );
        
        coordinates_arg = new_matrix;
        return;
    }
    else
    {
        MatrixXd original_Matrix = this->get_true_mtoc_coordinates();        
        double radius_of_mtoc = center_real_MTOC.norm();
        MatrixXd new_matrix = translate_true_mtoc_coordinates_to_radius( radius_of_mtoc );
        
        
        //cout<<original_Matrix - new_matrix<<endl;
        coordinates_arg = new_matrix;
        this->set_true_mtoc_coordinates( new_matrix );
        Vector3d center_2 = this->get_true_moc_center();

        
        axis_rotation = center_2.cross( center_real_MTOC );
        axis_rotation = axis_rotation / axis_rotation.norm();
        Quaternion<double> q;
        q = AngleAxis<double>( rotation_angle , axis_rotation );
        Vector3d new_center =  q * center_2;

        
    
        for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
        {
            coordinates_arg( dimension , 0 ) = new_center( dimension );
        }

    
        for( unsigned int point = 0 ; point  < this->number_of_true_mtoc_points ; point ++ )
        {
            //cout<<"point = "<<point<<endl;
            Vector3d point_tmp = this->get_true_mtoc_point( point );
            Vector3d new_point = q * point_tmp;
            

            for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
            {
                coordinates_arg( 3 * ( point + 1 ) + dimension , 0 ) = new_point( dimension );
            }           
        }
        
        this->set_true_mtoc_coordinates( original_Matrix );   
        
    }
    
    
}

MatrixXd ideal_MTOC::rotate_and_translate_true_mtoc_coordinates_2( Vector3d center_real_MTOC , MatrixXd coordinates_of_mtoc )
{
    MatrixXd coordinates_return = MatrixXd::Zero( 3 * ( MTOCparam::boundary_MTOC_points + 1 ) , 1 );
    //cout<<"LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL"<<endl;    
    double radius = center_real_MTOC.norm();
    

    Vector3d original_center_tmp( 0.0 , 0.0 , 1.0 );    
    Vector3d center_real_MTOC_tmp = center_real_MTOC / center_real_MTOC.norm();
    //center_real_MTOC_tmp = Vector3d( 1 , 0.0 , 0.0 );
    
    double cosinus_angle = center_real_MTOC_tmp.dot( original_center_tmp ) / ( original_center_tmp.norm() * center_real_MTOC_tmp.norm() );
    double rotation_angle = acos( cosinus_angle );
    Vector3d axis_rotation = original_center_tmp.cross( center_real_MTOC_tmp ); 
    if( axis_rotation.norm() < 1e-10 )
    {
        return coordinates_of_mtoc;
    }
    
    axis_rotation = axis_rotation / axis_rotation.norm();
    Quaternion<double> q;
    q = AngleAxis<double>( rotation_angle , axis_rotation );   
    Vector3d new_center =  q * original_center_tmp;
    //cout<<"radius = "<<new_center<<endl;
    
    new_center = new_center / new_center.norm() * radius;
    
    
    Vector3d first_point = this->original_MTOC_1_point;
    first_point( 2 ) = radius;
    
    double first_point_radius = first_point.norm();
    Vector3d first_point_tmp = first_point / first_point.norm();
    Vector3d new_first_point =  q * first_point_tmp;
    new_first_point = new_first_point / new_first_point.norm() * first_point_radius;
    

    
    
    
    for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
    {
        coordinates_return( dimension , 0 ) = new_center( dimension );
        coordinates_return( 3 + dimension , 0 ) = new_first_point( dimension );
    }
    
    double delta_angle = 2.0 * sim_of_Cell::PI / ( double ) MTOCparam::boundary_MTOC_points;
    Vector3d axis_for = new_center;// / new_center.norm();

    
    
    for( unsigned int point = 2 ; point <= MTOCparam::boundary_MTOC_points ; point ++ )
    {
        double local_angle = delta_angle * (double) ( point - 1 );   
        Quaternion<double> q_for;
        q_for = AngleAxis<double>( local_angle , axis_for ); 
        Vector3d new_point =  q_for * new_first_point;
        
        for( unsigned int dimension = 0 ; dimension < 3 ; dimension ++ )
        {
            coordinates_return( 3 * ( point ) + dimension , 0 ) = new_point( dimension );
        }
    }
    
    return coordinates_return;
    
}




Vector3d ideal_MTOC::get_original_micro_orientation( unsigned int index )
{
    Vector3d tangent = this->get_orientation( index ) - this->get_mtoc_point( index );  
    return tangent / tangent.norm();    
}




ideal_MTOC::~ideal_MTOC()
{

}

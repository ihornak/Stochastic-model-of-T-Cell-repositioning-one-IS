//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS
//Unused - for  future simulations
#ifndef IDEAL_MTOC_H
#define IDEAL_MTOC_H
#include <iostream>
#include "MTOCparam.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/SparseLU>
#include "Cell_parametres.h"
#include "simulationofCell.h"
using Eigen::MatrixXd;
using namespace Eigen;

using namespace std;

class ideal_MTOC
{
public:
ideal_MTOC();
ideal_MTOC(  unsigned int number_of_MTOC_points , Vector3d original_orientation );

ideal_MTOC(  unsigned int number_of_MTOC_points , Vector3d original_orientation , MatrixXd coordinates_tmp );
ideal_MTOC(const ideal_MTOC& other);
ideal_MTOC(  unsigned int number_of_MTOC_points , Vector3d original_orientation , MatrixXd coordinates_tmp ,  MatrixXd orientation_tmp );

ideal_MTOC& operator=(const ideal_MTOC& other);

//setters and getters
    unsigned int set_number_of_true_mtoc_points();
    unsigned int get_number_of_true_mtoc_points();
    
    
    void set_center(Vector3d point );
    void set_orientation_center( Vector3d point );
    void set_mtoc_center( Vector3d point );    
    void set_true_mtoc_center(Vector3d point );
    void set_original_MTOC_center( Vector3d argument );
    void set_original_MTOC_1_point( Vector3d argument );

    
    Vector3d get_center( );
    Vector3d get_orientation_center( );
    Vector3d get_moc_center( );
    Vector3d get_true_moc_center( );
    
    
    void set_point( unsigned int point_index , Vector3d point );
    void set_orientation( unsigned int point_index , Vector3d point );
    void set_mtoc_point( unsigned int point_index , Vector3d point );
    void set_true_mtoc_point( unsigned int point_index , Vector3d point );
    
    
    
    Vector3d get_point( unsigned int point_index );
    Vector3d get_orientation( unsigned int point_index );
    Vector3d get_mtoc_point( unsigned int point_index );
    Vector3d get_true_mtoc_point( unsigned int point_index );
    
    
    Vector3d get_tangent( unsigned int tangent_index );
    Vector3d get_orientation_tangent( unsigned int tangent_index );
    
    
    MatrixXd get_coordinates();
    MatrixXd get_orientations();
    MatrixXd get_mtoc_coordinates();
    MatrixXd get_true_mtoc_coordinates();
    
    
    void set_coordinates(  MatrixXd coordinates_tmp  );
    void set_orientations(  MatrixXd coordinates_tmp );
    void set_mtoc_coordinates(  MatrixXd coordinates_tmp );
    void set_true_mtoc_coordinates(  MatrixXd coordinates_tmp );
    
//ROTATION    
    void rotate_MTOC( Vector3d center_real_MTOC , MatrixXd& coordinates_arg );
    void rotate_orientations( Vector3d center_real_MTOC , MatrixXd& coordinates_arg );
    void rotate_mtoc_coordinates( Vector3d center_real_MTOC , MatrixXd& coordinates_arg );
    MatrixXd translate_true_mtoc_coordinates_to_radius( double radius_arg );
    void rotate_and_translate_true_mtoc_coordinates( Vector3d center_real_MTOC , MatrixXd& coordinates_arg );   
    Vector3d get_original_micro_orientation( unsigned int index );
    
    
    
    MatrixXd rotate_and_translate_true_mtoc_coordinates_2( Vector3d center_real_MTOC , MatrixXd coordinates_of_mtoc );
    
    
    
~ideal_MTOC();    
    
    
private:
    unsigned int number_of_points;
    unsigned int number_of_true_mtoc_points;
    MatrixXd mtoc_coordinates;
    MatrixXd coordinates;
    MatrixXd orientations_micro;       
    MatrixXd true_mtoc_coordinates;
    Vector3d original_orientation;

    Vector3d original_MTOC_center;
    Vector3d original_MTOC_1_point;
    
    
    
};

#endif // IDEAL_MTOC_H

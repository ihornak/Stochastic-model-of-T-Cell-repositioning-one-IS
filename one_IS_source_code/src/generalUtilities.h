//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS

#ifndef GENERALUTILITIES_H_
#define GENERALUTILITIES_H_


#include <iostream>
#include <string>
#include <vector>
#include <eigen3/Eigen/Dense>
#include "MTOCparam.h"
#include "simulationofCell.h"
#include "Cell.h"
#include "GeneralGeometry.h"
#include <time.h>
#include <ctime>
#include "Cell_parametres.h"
#include "IS_Capture_shrinkage_parameters.h"


#include <limits>
typedef std::numeric_limits< double > dbl;
using namespace std;


#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/SparseLU>
using Eigen::MatrixXd;
using namespace Eigen;


void get_program_ID( unsigned int number_of_micro , unsigned int number_of_pictures , double TIME_1 , double A_axis , double B_axis , double plane );
double calculate_Friction_Smooth_ER_Golgi( double radius_of_the_cell );
double calculateEffectiveFriction_Mitochondria( double lenght , double radius );
double calculateEffectiveFriction_Endoplasmatic_Reticulum( double radius_of_the_cell );
double calculateEffectiveFriction_Golgi_apparatus( double radius_of_the_cell );
double calculate_Friction_Rough_ER( double radius_of_the_cell );

double Endoplasmatic_Reticulum_coef( double radius_of_the_cell );
double Mitochondria_coef( double lenght , double radius );
double calculate_number_of_mitochondria( double radius_of_cell , double procentage , double lenght_of_mitochondria , double radius_of_mitochondria );
double calculate_number_of_mitochondria_micro_lengths( unsigned int number_of_microtubule , double lenght_of_microtubule , double lenght_of_mitochondria , double procentage );
double mito_cell_volume( unsigned int number_of_microtubule , double lenght_of_microtubule , double lenght_of_mitochondria , double procentage );
double compute_time( double radius_of_the_cell );
double calculateEffectiveFriction_all_Mitochondria( unsigned int microtubule_number );
void erase_distance_measurements();








#endif /* GENERALUTILITIES_H_ */

//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS
//Unused, for future simulations
#ifndef MTOC_SIMPLE_H
#define MTOC_SIMPLE_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/SparseLU>
using Eigen::MatrixXd;
using namespace Eigen;
#include "simulationofCell.h"
#include <iostream>
using namespace std;


class MTOC_simple
{
public:
MTOC_simple();
MTOC_simple( Vector3d position_tmp );
MTOC_simple(const MTOC_simple& other);

Vector3d get_position();
void set_position( Vector3d position_tmp );
double get_friction();
void set_friction( double friction_tmp );

void oneStepMidStepAlgorithm_1_half( Vector3d force );
void oneStepMidStepAlgorithm_2_half( Vector3d force , Vector3d position_tmp );



~MTOC_simple();
MTOC_simple& operator=(const MTOC_simple& other);
bool operator==(const MTOC_simple& other) const;

private:
    Vector3d position;
    double friction;
};

#endif // MTOC_SIMPLE_H

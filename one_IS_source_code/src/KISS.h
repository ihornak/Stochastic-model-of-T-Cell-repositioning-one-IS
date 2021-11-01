//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS

#ifndef KISS_H_
#define KISS_H_
#include <stdio.h>
#include <climits>
#include <math.h>
#include <iostream>
#include "mersenne.h"
#include "simulationofCell.h"
#include "IS_Cortical_Sl_parameter.h"

using namespace std;


static unsigned long long


x = 9781641336351611ULL , c = 816813135002ULL,
y = 9361448128545ULL , z = 981641ULL , t;


#define MWC (t=(x<<58)+c, c=(x>>6), x+=t, c+=(x<t), x)
#define XSH ( y^=(y<<13), y^=(y>>17), y^=(y<<43) )
#define CNG ( z=6906969069LL*z+1234567 )
#define KISS (MWC+XSH+CNG)


double KISS_function();
double randx();
double rand_x( double min ,  double max );
double normalRandom();
double normal_Random( double mu, double sigma);
int RAND_INT( int min , int max );
void seed_generator_with_time();
double triangular_distribution();




#endif 

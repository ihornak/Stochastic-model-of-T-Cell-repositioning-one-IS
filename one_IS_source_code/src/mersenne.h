//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS
#include <stdio.h>
#include <stdio.h>
#include <climits>
#include <math.h>
#include <iostream>
#include <random>
#include <chrono>
#ifndef MERSENNE_H_
#define MERSENNE_H_
using namespace std;

#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL 
#define LM 0x7FFFFFFFULL 



namespace mersenne_twisters_64
{
	extern std::vector< std::mt19937_64 > vector_of_mersenne_twisters;
	void initialize_generators( unsigned int number_of_threads );
	void initialize_generators_stable_seed( unsigned int number_of_threads );
	double triangular_distribution( double upper_boundary );
}



static unsigned long long mt[NN];
static int mti=NN+1;
void init_genrand64(unsigned long long seed);
void init_by_array64(unsigned long long init_key[], unsigned long long key_length);
unsigned long long genrand64_int64(void);
double randx_mersenne();
double rand_x_mersenne( double min ,  double max );
double normalRandom_mersenne();
int RAND_INT_mersenne( int min , int max );

#endif

//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS
#include "KISS.h"


void seed_generator_with_time()
{

        unsigned int random_variable_1 = genrand64_int64();
	unsigned int random_variable_2 = genrand64_int64();
        unsigned int random_variable_3 = genrand64_int64();
	unsigned int random_variable_4 = genrand64_int64();

	x = random_variable_1;
	c = random_variable_2;
	y = random_variable_3;
	z = random_variable_4;

}



double KISS_function()
{
	return  ( double )KISS;
}

double randx()
{
	return ( double )KISS / ( double ) ( ULLONG_MAX  + 1.0);
}

double rand_x( double min ,  double max )
{
	return randx() * ( max - min ) + min;
}

double normalRandom()
{
  double u1=randx();
  double u2=randx();
  return cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1));
}

double normal_Random( double mu, double sigma)
{
  double ran = mu + sigma * normalRandom();
  return ran;
}

int RAND_INT( int min , int max )
{
	double lower = ( double ) min;
	double upper = ( double ) max - 1e-12; // boudary
	double tmp = rand_x( lower ,  upper );
	int final = ( unsigned int )  tmp;
	return final;
}



double triangular_distribution()
{
	double R = IS_Cortical_Sl_parameter::radius;
	double x = randx();
	
	double result = sqrt( x ) * R;
	return result;
}




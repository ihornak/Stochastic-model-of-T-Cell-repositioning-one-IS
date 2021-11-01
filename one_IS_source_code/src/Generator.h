//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS

#ifndef GENERATOR_H_
#define GENERATOR_H_
#include <random>
#include <iostream>






class Generator {
public:

	//Microtubules in unconstrained 3D space
	Generator( unsigned int seed ); 
	double get_probability();
	


	virtual ~Generator();


private:

	std::mt19937_64 rd;
	std::uniform_real_distribution<> dist_uniform;
    
};

#endif /* MICROTUBULE_H_ */

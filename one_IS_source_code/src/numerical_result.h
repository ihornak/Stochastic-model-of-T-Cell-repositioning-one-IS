#ifndef NUMERICAL_RESULT_H
#define NUMERICAL_RESULT_H


#include <iostream>
#include <random>

#include "MTOC2.h"
#include "Microtubule.h"
#include "Cell.h"
#include "ISCorticalSl.h"
#include "generalUtilities.h"
#include "Surface.h"
#include "Dynein.h"
#include "IS_Capture_shrinkage_parameters.h"
#include "MTOCparam.h"
#include "mersenne.h"
#include <random>
#include <iostream>





namespace numerical_result
{

    void numerical_results_f1( unsigned int number_of_samples );
    void repolarization_sample( unsigned int index );	
    void test_parallel( );


    //Creates the file with the parameters of the simulation
    void create_document_file();
}






#endif // NUMERICAL_RESULT_H



#ifndef MICROTUBULE_H_
#define MICROTUBULE_H_
#include "MTOCparam.h"
#include "KISS.h"
#include <fstream>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <string>
#include <map>
#include "mersenne.h"
#include <random>
#include <iostream>


#include "simulationofCell.h"
#include "ISCaptShrinkage.h"
#include "Cell_parametres.h"
#include "IS_Cortical_Sl_parameter.h"
#include "IS_Capture_shrinkage_parameters.h"
#include "GeneralGeometry.h"
#include "IS_Dynein_Cell_surface.h"
#include "Dynein.h"
#include "Dynein_real.h"


#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>	
#include <Eigen/SparseQR>
#include <eigen3/Eigen/IterativeLinearSolvers>
using Eigen::MatrixXd;
using namespace Eigen;
typedef Eigen::Triplet<double> T;

//FOR THE COMMENTARY OF THE FUNCTIONS SEE https://github.com/ihornak/Stochastic-model-of-T-cell-repositioning-two-IS

class Microtubule {
public:

	//Microtubules in unconstrained 3D space
	Microtubule();
	Microtubule( unsigned int ID );
	Microtubule( Vector3d MTOC , unsigned int ID );
	Microtubule( Vector3d MTOC , Vector3d first_Point , unsigned int ID );
	//microtubule ID is number of microtubule in polygon
	Microtubule( Vector3d first_Point , Vector3d second_Point , unsigned int ID , unsigned int poly);   //PREDELAT POZDEJI!!!!!
	Microtubule( Vector3d first_Point , Vector3d orientation , unsigned int ID , unsigned int poly , double a_axis , double b_axis );
	Microtubule( Vector3d first_Point , Vector3d orientation , unsigned int ID , unsigned int poly , double a_axis , double b_axis , unsigned int number_of_points  );
    Microtubule( Vector3d first_Point , Vector3d orientation , unsigned int ID , unsigned int poly , unsigned int side_arg , unsigned int MTOC_point_arg , unsigned int number_of_points  );
    Microtubule( Vector3d first_Point , Vector3d second_Point , Vector3d orientation , unsigned int ID , unsigned int poly , unsigned int side_arg , unsigned int MTOC_point_arg , unsigned int second_MTOC_point_arg , unsigned int number_of_points  );
	Microtubule( Vector3d first_Point , Vector3d second_Point , Vector3d orientation , unsigned int ID , unsigned int poly , unsigned int side_arg , unsigned int MTOC_point_arg , unsigned int second_MTOC_point_arg , unsigned int number_of_points , bool rovna_mikrotubula );
    Microtubule( Vector3d first_Point , Vector3d second_Point , Vector3d orientation , unsigned int ID , unsigned int poly , unsigned int side_arg , unsigned int MTOC_point_arg , unsigned int second_MTOC_point_arg , unsigned int number_of_points , double first_bead_distance );

	Microtubule( const char* filename , double rad_of_Cell , Vector3d IS , unsigned int ID );
	Microtubule( const Microtubule & tmp );


        //Setters and getters
        void setCoordinates( MatrixXd &coordinatesTmp );
        void setPoint( unsigned int number_of_point ,  Vector3d position );
        MatrixXd getCoordinates() const;
        double getRestDist() const;
        double getKappa() const;
        double get_effective_friction();
        double get_effective_friction_whole_microtubule();
        Vector3d get_IS_position_catching( ) const;
        unsigned int getID()  const;
        unsigned int getSide()  const;
        unsigned int get_polygon_number()  const;
        unsigned int get_MTOC_point();
        unsigned int get_MTOC_opposite_point();
        unsigned int get_dynein_index() const;
        void alter_flexural_rigidity( double ratio );
        void set_dynein_index( unsigned int new_value );
        void set_IS_position_catching( Vector3d IS_vector );
        Vector3d getPoint( unsigned int index )  const;
        Vector3d get_last_Point( )  const;
        Vector3d getTangent2( unsigned int index ) const;
	Vector3d get_last_Tangent( )  const;
        double get_lenght_of_microtubule() const; 
	double get_lenght_of_microtubule_outside_MTOC() const;
        double get_distance_to_lower_bead_with_index( unsigned int index ) const;
        MatrixXd get_lenght_of_tangents(){ return this->lenght_of_tangents; };
        void set_lenght_after_catching( double lenght_after_catching_arg );
        double get_lenght_after_catching( );
	Microtubule& operator=( const Microtubule &tmp );
	unsigned int getNumberOfPoints()  const;
	double calculateEffectiveFriction( );
	double calculateEffectiveFriction_Howard( );
	
    
	Vector3d der_scalar_index( unsigned int scalar , unsigned int index);
	void getN_i_mu( MatrixXd &n_i_mu );
	void getMatrixes( MatrixXd &inv_Matrix , MatrixXd &projection_Matrix );
	void getMatrixes_2( MatrixXd &inv_Matrix , MatrixXd &projection_Matrix );
        void getMatrixes_3(MatrixXd &projection_Matrix );
	void getMatrixes_4(MatrixXd &projection_Matrix );	
	void getMatrixes_Sparse( MatrixXd &inv_Matrix , MatrixXd &projection_Matrix );
        
    	void get_Sparse_Projection_Matrix_micro_7(  MatrixXd &projection_Matrix );
    	void get_Sparse_Projection_Matrix_micro_7_second( MatrixXd &projection_Matrix );
    	void get_Sparse_Projection_Matrix_micro_7_third( MatrixXd &projection_Matrix );
    	void get_Sparse_Projection_Matrix_micro_7_fourth( MatrixXd &projection_Matrix );
    	void getBendingForces( MatrixXd &bendingForce );
    	void set_bending_matrix( );
    	void getBendingForces_2( MatrixXd &bendingForce );
	void getRandomForces(  double timeStep , MatrixXd &randomForces );

	//Returns metric Forces
	void metricForceWLC( MatrixXd &metricForces , MatrixXd &G_uv_Inverse );
	void metricForceConsLenght( MatrixXd &metricForces , MatrixXd &G_uv_Inverse );

	//Forces cause by dynein in the cell - Constrained by geometry of the cell
	Vector3d dynein_force_capture_shrinkage( );
    	Vector3d dynein_force_capture_shrinkage_2();
	Vector3d dynein_force_cortical_sliding( );

	bool distanceControl();
	void resizeMicrotubule();
    	void resizeMicrotubule( Vector3d orientation );
    	void resizeMicrotubule_first_segment_different( );
    	void resizeMicrotubule_with_different_tangent_lenghts( );
    	bool distanceEllipsoid( );

	void midStepAlgorithm( double startTime , double endTime );
	void midStepAlgorithm( double startTime , double endTime , IS_Capt_Shrinkage  IS_Capt_Shrinkage_argument  );
	MatrixXd relaxationTimeMidStep( double startTime , double endTime , double DELTATT );
	double equipartition_Theorem_analytical();


	void Substract_one_Bead();
    	void add_one_final_Bead_catching_position();
	void Substract_All_Bead();
	void smoother_Microtubule();
	void destroy_Microtubule_capture_shrinkage();
	void destroy_Microtubule_cortical_sliding();
	void IS_micro_control_one_micro_capture_shrinkage( IS_Capt_Shrinkage IS_Capt_Shrinkage_argument );
	//Printing of microtubule
        void print_Microtubule();
        void print_Microtubule_Tangent();
        void print_Dynein_surface_points(); 
        
        bool confirm_inner_ellipsoid( Vector3d position , double a_axis , double b_axis );

    	void add_catching_point_IS_sliding( unsigned int bead_number , Vector3d position );
    	void set_IS_sliding_catching_points( std::map < unsigned int , Vector3d> mapa_tmp );
    	std::map < unsigned int , Vector3d> get_plane_axis_point();
    	void print_IS_sliding_catching_points( );
    	void check_IS_catching_points();
    	void check_IS_catching_points_2();
    	void check_IS_sliding_all_points();
    	void check_IS_sliding_all_points_2();    
	

        //MICRO interactions with cortical sliding
        void add_pair( std::pair < Vector3d ,double  > point_abscissa );
        std::vector< std::pair < Vector3d ,double  > > get_Dynein_abscissa(){ return this->Dynein_motors_2; };
        std::vector< Vector3d  > get_Dynein_points_and_erase();

       //Interaction with dynein
        Vector3d force_real_dynein_one_pair(  std::pair < Vector3d , double > pair_tmp   );
        std::vector< Vector3d > stepping_detach_real_dynein_abscissa_projection( );
        void force_dynein_abscissa_real_dynein( MatrixXd& force_dynein );    
        void set_to_number_of_points( unsigned int number_of_point );
        unsigned int get_index_according_to_abscissa( double abscissa );
        Vector3d get_attachment_point_according_to_abscissa( double abscissa ); 
        void set_lenght_of_tangents();
	std::vector<Vector3d> control_length_of_micro_IS_2();
	std::vector<Vector3d> get_catching_points_IS_sliding( );

       //Interaction with real dynein
        std::vector<Vector3d>control_and_resizing_and_motor_adjustment_depolimerizing_microtubule();
	std::vector<Vector3d>  control_and_resizing_and_motor_adjustment_depolimerizing_microtubule_2();
        std::vector<Vector3d> control_motor_detachment_with_depolimerization();
        std::vector<Vector3d> get_dynein_points_in_IS_and_erase();
        unsigned int get_number_of_dynein_points_IS();  
        
        //General utilities
        void distribute_force_on_beads( MatrixXd& force_on_microtubule , Vector3d force , unsigned int bead_segment , double t );
        
        
        //Performs one step using Mid-Step algorithm
	void oneStepMidStepAlgorithm();
	void oneStepMidStepAlgorithm_1_half( MatrixXd& randomForces , MatrixXd& extrenal_Forces );
	void oneStepMidStepAlgorithm_1_half_producing_random( MatrixXd& randomForces , MatrixXd& extrenal_Forces );
	void oneStepMidStepAlgorithm_2_half( MatrixXd& original_Coordinates , MatrixXd& randomForces , MatrixXd& extrenal_Forces );
	void oneStepMidStepAlgorithm_2_half_projected( MatrixXd& original_Coordinates , MatrixXd& randomForces , MatrixXd& other_Forces );
	void project_force( MatrixXd& original_force , MatrixXd& projected_forces );
    	void Euler_algorithm( MatrixXd& randomForces , MatrixXd& extrenal_Forces );    
  	std::vector< Vector3d  > get_Dynein_points_without_erasing(  );              
	virtual ~Microtubule();


private:
	unsigned int numberOfPoints;    
    	unsigned int side;
    	unsigned int polygon;
	unsigned int microtubule_id;
    	unsigned int MTOC_point;
    	unsigned int second_MTOC_point;     
	double restDistancePoints;
    	double restDistancePoints_first;
	double kappa;
	double effective_friction;
	MatrixXd coordinates;
	Vector3d IS_position_catching;
	unsigned int dydein_index;
        std::vector< std::pair < Vector3d , double > > Dynein_motors_2;
        MatrixXd lenght_of_tangents;
        double lenght_after_catching;
	SparseMatrix<double> b_MATRIX;
	SparseMatrix<double> Kronecker_Delta;

    
};

#endif /* MICROTUBULE_H_ */

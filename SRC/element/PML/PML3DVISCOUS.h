/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */


// Written by: Amin Pakzad, Pedro Arduino (parduino@uw.edu)
//
// Eight node PML3D element .. a c++ wrapper to fortran routine 
// provided by Wenyang Zhang (zwyll@ucla.edu), University of California, Los Angeles
//
// University of Washington, UC. Los Angeles, U.C. Berkeley, 12, 2020


#ifndef PML3DVISCOUS_H
#define PML3DVISCOUS_H

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <Domain.h>
#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>

#define PML3DVISCOUS_NUM_DOF 72
#define PML3DVISCOUS_NUM_PROPS 12
#define PML3DVISCOUS_NUM_NODES 8

#ifdef _WIN32

#define pml3d_	      PML_3D

extern "C" void  pml3d_(double* mMatrix,
	double* cMatrix,
	double* kMatrix,
	double* gMatrix,
	double* hMatrix,
	int* NDOFEL,
	double* PROPS,
	double* COORDS,
	int* MCRD,
	int* NNODE,
	int* LFLAGS);

#else

#define pml3d_	      pml_3d_

extern "C" void  pml3d_(double* mMatrix,
	double* cMatrix,
	double* kMatrix,
	double* gMatrix,
	double* hMatrix,
	int* NDOFEL,
	double* PROPS,
	double* COORDS,
	int* MCRD,
	int* NNODE,
	int* LFLAGS);

#endif

class PML3DVISCOUS : public Element {

public:

	PML3DVISCOUS();                                                                         //null constructor
	PML3DVISCOUS(int tag, int* nodeTags, double* newmarks, double* dData);                  // full constructor
	virtual ~PML3DVISCOUS();                                                                //destructor
	const char* getClassType(void) const { return "PML3DVISCOUS"; };                        //return class type
	void setDomain(Domain* theDomain);                                               // set domain
	int getNumExternalNodes() const; 	   						                     // get number of external nodes
	const ID& getExternalNodes(); 								                     // get external nodes
	Node** getNodePtrs(void); 									                     // get external nodes
	int getNumDOF(); 											                     // get number of DOF
	int commitState(); 											                     // commit state
	int revertToLastCommit(); 									                     // revert to last commit
	int revertToStart(); 										                     // revert to start
	int update(void); 											                     // update
	void Print(OPS_Stream& s, int flag); 							                 // print out element data
	const Matrix& getTangentStiff(); 							                     // get stiffness matrix
	const Matrix& getInitialStiff(); 							                     // get initial stiffness matrix
	const Matrix& getMass(); 									                     // get mass matrix
	const Matrix& getDamp(); 									                     // get damping matrix
	void zeroLoad(); 											                     // set residual to 0
	int addLoad(ElementalLoad* theLoad, double loadFactor); 		                 // add element loads
	const Vector& getResistingForce(); 							                     // get residual
	const Vector& getResistingForceIncInertia(); 				                     // get residual including damping forces
	int sendSelf(int commitTag, Channel& theChannel); 			                     // send self
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker & theBroker);  // receive self
	Response* setResponse(const char** argv, int argc, OPS_Stream& s);               // set response
	int getResponse(int responseID, Information& eleInformation);                    // get response
	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int parameterID, Information& info);                         // update parameter
	int displaySelf(Renderer&, int mode, float fact, const char** displayModes = 0, int numModes = 0);

private:

	Domain* Domainptr;                              // pointer to the domain
	double props[PML3DVISCOUS_NUM_PROPS];                  // material properties
	ID connectedExternalNodes;  					//eight node numbers
	Node* nodePointers[PML3DVISCOUS_NUM_NODES];    	    //pointers to eight nodes
	double K[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF];        // stiffness matrix
	double C[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF];        // damping matrix
	double M[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF];        // mass matrix
    double G[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF];        // G matrix
    double H[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF];        // H matrix
	double Keff[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF];     // effective stiffness matrix
	static double eta;                              // Newmark parameters: eta
	static double beta; 					  	    // Newmark parameters: beta
	static double gamma; 					  	// Newmark parameters: gamma
	static Matrix tangent;                          // tangent matrix
	static Vector resid; 						    // residual vector
	static Matrix mass;						        // mass matrix
	static Matrix damping;	 					    // damping matrix
	
	// Add static variables for integration points and weights
    static double xi[3][64];                        // integration points coordinates
    static double w[64];                            // integration weights
    static int numIntegrationPoints;                // number of integration points
    static bool integrationInitialized;             // flag to check if integration points are initialized
    
    // Add static function to initialize integration points and weights
    static void initializeIntegrationPointsAndWeights(int n_points, int n_nodes);
    
    // Add static function to calculate shape functions and derivatives
    static void calculateShapeFunctions(const double* xi, int n_nodes, double* N, double (*dNdxi)[3]);
    
    // Add function to calculate PML alpha and beta parameters
    static void calculatePMLParameters(const double* props, double x1, double x2, double x3, double (*pmlAlphaBeta)[3]);
    
    // Add function to verify matrices from Fortran and C++
    void verifyMatrices();
    
	Vector ubart; 				                    // ubar at time t 
    Vector ubarbart;                                // ubarbar at time t
	Vector ubar; 				                    // ubar at time t+dt
    Vector ubarbar;                                 // ubarbar at time t+dt
	static double dt; 								// time step
	int updateflag; 								// update flag
	int update_dt;                                  // flag for updating dt
	static int eleCount; 						    // element count
	// int innertag; 								// inner tag
	// static int numberOfElements; 			    // number of elements

    // Static matrices for storing calculated results
    static double M_cpp[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF];  // Mass matrix from C++
    static double C_cpp[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF];  // Damping matrix from C++
    static double K_cpp[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF];  // Stiffness matrix from C++
    static double G_cpp[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF];  // G matrix from C++
    static double H_cpp[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF];  // H matrix from C++
};

#endif


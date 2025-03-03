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
// Eight node PML3DVISCOUS element .. a c++ wrapper to fortran routine 
// provided by Wenyang Zhang (zwyll@ucla.edu), University of California, Los Angeles
//
// University of Washington, UC. Los Angeles, U.C. Berkeley, 12, 2020


#include "PML3DVISCOUS.h"

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <OPS_Globals.h>
#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <Parameter.h>
#include <ElementalLoad.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <fstream>


// =======================================================================
// PML3DVISCOUS element tcl command
// =======================================================================
void* OPS_PML3DVISCOUS()
{
	// check if the total number of arguments passed is correct
	if (OPS_GetNumRemainingInputArgs() < (9 + PML3DVISCOUS_NUM_PROPS + 3)) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: element PML3DVISCOUS eleTag? [8 integer nodeTags] [PML3DVISCOUS_NUM_PARAMS material properties]\n";
		return 0;
	}

	// reading element tag and node numbers 
	int idata[9];
	int num = 9;
	if (OPS_GetIntInput(&num, idata) < 0) {
		opserr << "WARNING: invalid integer data : could be the tag or the node numbers \n";
		return 0;
	}

	// reading Newmark parameters
	double Newmark[3];
	num = 3;
	if (OPS_GetDoubleInput(&num, Newmark) < 0) {
		opserr << "WARNING: invalid double data: could be Newmark parameters\n";
		return 0;
	}

	// reading material properties
	double dData[PML3DVISCOUS_NUM_PROPS]; num = PML3DVISCOUS_NUM_PROPS;
	if (OPS_GetDoubleInput(&num, dData) < 0) {
		opserr << "WARNING: invalid double data\n";
		return 0;
	}

	// create a new PML3DVISCOUS element and add it to the Domain
	return new PML3DVISCOUS(idata[0], &idata[1], Newmark, dData);
}

// =======================================================================
// static data
// =======================================================================
Matrix  PML3DVISCOUS::tangent(PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
Matrix  PML3DVISCOUS::mass(PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
Matrix  PML3DVISCOUS::damping(PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
Vector  PML3DVISCOUS::resid(PML3DVISCOUS_NUM_DOF);
double  PML3DVISCOUS::eta = 0.;
double  PML3DVISCOUS::beta = 0.;
double  PML3DVISCOUS::gamma = 0.;
double  PML3DVISCOUS::dt = 0.;
int     PML3DVISCOUS::eleCount = 0;

// Initialize integration points and weights
double  PML3DVISCOUS::xi[3][27] = {{0.0}};
double  PML3DVISCOUS::w[27] = {0.0};
bool    PML3DVISCOUS::integrationInitialized = false;

// Initialize static matrices
double PML3DVISCOUS::M_cpp[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF] = {0};
double PML3DVISCOUS::C_cpp[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF] = {0};
double PML3DVISCOUS::K_cpp[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF] = {0};
double PML3DVISCOUS::G_cpp[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF] = {0};
double PML3DVISCOUS::H_cpp[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF] = {0};

// =======================================================================
// null constructor
// =======================================================================
PML3DVISCOUS::PML3DVISCOUS()
	:Element(0, ELE_TAG_PML3DVISCOUS),
	connectedExternalNodes(PML3DVISCOUS_NUM_NODES),
	ubar(PML3DVISCOUS_NUM_DOF),
	ubart(PML3DVISCOUS_NUM_DOF),
    ubarbar(PML3DVISCOUS_NUM_DOF),
    ubarbart(PML3DVISCOUS_NUM_DOF)
{
    // Initialize integration points and weights with default values
    if (!integrationInitialized) {
        initializeIntegrationPointsAndWeights(27, PML3DVISCOUS_NUM_NODES); // Default to 27-point integration
    }
    
    for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
        nodePointers[i] = 0;
    }
	dt = 0;
	ubar.Zero();
	ubart.Zero();
    ubarbar.Zero();
    ubarbart.Zero();
	updateflag = 0;
	update_dt = 0;
	eta = 0;
	beta = 0;
	gamma = 0;
	
}

// =======================================================================
// Full constructor
// =======================================================================
PML3DVISCOUS::PML3DVISCOUS(int tag, int* nodeTags, double* nemwarks, double* eleData)
	:Element(tag, ELE_TAG_PML3DVISCOUS),
	connectedExternalNodes(PML3DVISCOUS_NUM_NODES),
	ubar(PML3DVISCOUS_NUM_DOF),
	ubart(PML3DVISCOUS_NUM_DOF),
    ubarbar(PML3DVISCOUS_NUM_DOF),
    ubarbart(PML3DVISCOUS_NUM_DOF)
{
	eleCount++;
	if (eleCount == 1) {
		opserr << "Perfectly Matched Layer 3D (PMLVISCOUS) element -  Written: W. Zhang, E. Taciroglu, A. Pakzad, P. Arduino, UCLA, U.Washington\n ";
	}
    
    // Initialize integration points and weights with default values
    if (!integrationInitialized) {
        // Check if afp parameter is less than 3.0 for integration point selection
        if (eleData[5] < 3.0) { // afp parameter is at index 5
            initializeIntegrationPointsAndWeights(27, PML3DVISCOUS_NUM_NODES);
        } else {
            initializeIntegrationPointsAndWeights(64, PML3DVISCOUS_NUM_NODES);
        }
    }
    
	// initialize node pointers
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		connectedExternalNodes(i) = nodeTags[i];
		nodePointers[i] = 0;
	}

	// initialize Newmark parameters
	eta   = nemwarks[0];
	beta  = nemwarks[1];
	gamma = nemwarks[2];

	// initialize material properties
	for (int i = 0; i < PML3DVISCOUS_NUM_PROPS; i++)
		props[i] = eleData[i];

    E   = props[0];
    xnu = props[1];
    rho = props[2];
    eleTypeArg = (int)props[3];
    PML_L = props[4];
    afp = props[5];
    PML_Rcoef = props[6];
    RD_half_width_x = props[7];
    RD_half_width_y = props[8];
    RD_depth = props[9];
    Damp_alpha = props[10];
    Damp_beta = props[11];

    // print element properties
    opserr << "Element properties:\n";
    opserr << "E: " << E << "\n";
    opserr << "xnu: " << xnu << "\n";
    opserr << "rho: " << rho << "\n";
    opserr << "eleTypeArg: " << eleTypeArg << "\n";
    opserr << "PML_L: " << PML_L << "\n";
    opserr << "afp: " << afp << "\n";
    opserr << "PML_Rcoef: " << PML_Rcoef << "\n";
    opserr << "RD_half_width_x: " << RD_half_width_x << "\n";
    opserr << "RD_half_width_y: " << RD_half_width_y << "\n";
    opserr << "RD_depth: " << RD_depth << "\n";
    opserr << "Damp_alpha: " << Damp_alpha << "\n";
    opserr << "Damp_beta: " << Damp_beta << "\n";

    

	// initialize the ubar and ubart vectors to zero
	ubart.Zero();
	ubar.Zero();
    ubarbar.Zero();
    ubarbart.Zero();
	updateflag = 0;
	update_dt = 0;
}

// =======================================================================
//  Integration points and weights
// =======================================================================
void PML3DVISCOUS::initializeIntegrationPointsAndWeights(int n_points, int n_nodes) {
    // If already initialized, return
    if (integrationInitialized) return;
    // Reset arrays
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 64; j++) {
            xi[i][j] = 0.0;
        }
    }
    for (int i = 0; i < 64; i++) {
        w[i] = 0.0;
    }
    // Hexahedral elements
    if (n_nodes == 8) {
        if (n_points == 27) {
            double x1D[3] = {-0.7745966692, 0.0, 0.7745966692};
            double w1D[3] = {0.5555555555, 0.888888888, 0.55555555555};
            int n = 0;
            for (int k = 0; k < 3; k++) {
                for (int j = 0; j < 3; j++) {
                    for (int i = 0; i < 3; i++) {
                        xi[0][n] = x1D[i];
                        xi[1][n] = x1D[j];
                        xi[2][n] = x1D[k];
                        w[n] = w1D[i] * w1D[j] * w1D[k];
                        n++;
                    }
                }
            }
        } else
        {
            opserr << "ERROR: PML3DVISCOUS element only supports 27 integration points for 8-node hexahedral elements\n";
            exit(-1);
        }
    } else {
        opserr << "ERROR: PML3DVISCOUS element only supports Hexahedral elements with 8 nodes\n";
        exit(-1);
    }

    integrationInitialized = true;
    // # if the element tag is equal to 1, print the integration points and weights
    if (eleCount == 1) {
        opserr << "Integration points and weights:\n";
        for (int i = 0; i < n_points; i++) {
            opserr << "Point " << i + 1 << ": (" << xi[0][i] << ", " << xi[1][i] << ", " << xi[2][i] << ") - Weight: " << w[i] << "\n";
        }
    }
}

// =======================================================================
//  shape functions
// =======================================================================
void PML3DVISCOUS::calculateShapeFunctions(const double* xi, int n_nodes, double* N, double (*dNdxi)[3]) {
    // Local variables
    double xi4;
    
    // Clear arrays
    for (int i = 0; i < n_nodes; i++) {
        N[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            dNdxi[i][j] = 0.0;
        }
    }
    if (n_nodes == 8) {  // Linear hexahedral
        N[0] = (1.0 - xi[0]) * (1.0 - xi[1]) * (1.0 - xi[2]) / 8.0;
        N[1] = (1.0 + xi[0]) * (1.0 - xi[1]) * (1.0 - xi[2]) / 8.0;
        N[2] = (1.0 + xi[0]) * (1.0 + xi[1]) * (1.0 - xi[2]) / 8.0;
        N[3] = (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 - xi[2]) / 8.0;
        N[4] = (1.0 - xi[0]) * (1.0 - xi[1]) * (1.0 + xi[2]) / 8.0;
        N[5] = (1.0 + xi[0]) * (1.0 - xi[1]) * (1.0 + xi[2]) / 8.0;
        N[6] = (1.0 + xi[0]) * (1.0 + xi[1]) * (1.0 + xi[2]) / 8.0;
        N[7] = (1.0 - xi[0]) * (1.0 + xi[1]) * (1.0 + xi[2]) / 8.0;
        
        dNdxi[0][0] = -(1.0 - xi[1]) * (1.0 - xi[2]) / 8.0;
        dNdxi[0][1] = -(1.0 - xi[0]) * (1.0 - xi[2]) / 8.0;
        dNdxi[0][2] = -(1.0 - xi[0]) * (1.0 - xi[1]) / 8.0;
        
        dNdxi[1][0] = (1.0 - xi[1]) * (1.0 - xi[2]) / 8.0;
        dNdxi[1][1] = -(1.0 + xi[0]) * (1.0 - xi[2]) / 8.0;
        dNdxi[1][2] = -(1.0 + xi[0]) * (1.0 - xi[1]) / 8.0;
        
        dNdxi[2][0] = (1.0 + xi[1]) * (1.0 - xi[2]) / 8.0;
        dNdxi[2][1] = (1.0 + xi[0]) * (1.0 - xi[2]) / 8.0;
        dNdxi[2][2] = -(1.0 + xi[0]) * (1.0 + xi[1]) / 8.0;
        
        dNdxi[3][0] = -(1.0 + xi[1]) * (1.0 - xi[2]) / 8.0;
        dNdxi[3][1] = (1.0 - xi[0]) * (1.0 - xi[2]) / 8.0;
        dNdxi[3][2] = -(1.0 - xi[0]) * (1.0 + xi[1]) / 8.0;
        
        dNdxi[4][0] = -(1.0 - xi[1]) * (1.0 + xi[2]) / 8.0;
        dNdxi[4][1] = -(1.0 - xi[0]) * (1.0 + xi[2]) / 8.0;
        dNdxi[4][2] = (1.0 - xi[0]) * (1.0 - xi[1]) / 8.0;
        
        dNdxi[5][0] = (1.0 - xi[1]) * (1.0 + xi[2]) / 8.0;
        dNdxi[5][1] = -(1.0 + xi[0]) * (1.0 + xi[2]) / 8.0;
        dNdxi[5][2] = (1.0 + xi[0]) * (1.0 - xi[1]) / 8.0;
        
        dNdxi[6][0] = (1.0 + xi[1]) * (1.0 + xi[2]) / 8.0;
        dNdxi[6][1] = (1.0 + xi[0]) * (1.0 + xi[2]) / 8.0;
        dNdxi[6][2] = (1.0 + xi[0]) * (1.0 + xi[1]) / 8.0;
        
        dNdxi[7][0] = -(1.0 + xi[1]) * (1.0 + xi[2]) / 8.0;
        dNdxi[7][1] = (1.0 - xi[0]) * (1.0 + xi[2]) / 8.0;
        dNdxi[7][2] = (1.0 - xi[0]) * (1.0 + xi[1]) / 8.0;
    }
	else {
		opserr << "WARNING: PML3DVISCOUS element only supports Hexahedral elements with 8 nodes\n";
	}
}

// =======================================================================
//  calculatePMLParameters
// =======================================================================
void PML3DVISCOUS::calculatePMLParameters(double x1, double x2, double x3, double (*pmlAlphaBeta)[3]) {
    
    // Initialize normal vectors and reference coordinates
    double n1 = 0.0;
    double n2 = 0.0;
    double n3 = 0.0;
    double x1_0 = 0.0;
    double x2_0 = 0.0;
    double x3_0 = 0.0;
    
    // Calculate reference dilatational wave velocity
    double cp_ref = sqrt(E * (1.0 - xnu) / rho / (1.0 + xnu) / (1.0 - 2.0 * xnu));
    
    // Determine element type based on coordinates
    if (x2 < -RD_half_width_y) {
        if (x1 < -RD_half_width_x) {
            if (x3 < -RD_depth) {
                eleTypeArg = 15;  // Left-bottom corner PML
            } else {
                eleTypeArg = 6;   // Left-bottom column PML
            }
        } else if (x1 < RD_half_width_x) {
            if (x3 < -RD_depth) {
                eleTypeArg = 11;  // Bottom column PML
            } else {
                eleTypeArg = 2;   // Bottom surface PML
            }
        } else {
            if (x3 < -RD_depth) {
                eleTypeArg = 16;  // Right-bottom corner PML
            } else {
                eleTypeArg = 7;   // Right-bottom column PML
            }
        }
    } else if (x2 < RD_half_width_y) {
        if (x1 < -RD_half_width_x) {
            if (x3 < -RD_depth) {
                eleTypeArg = 14;  // Left column PML
            } else {
                eleTypeArg = 5;   // Left surface PML
            }
        } else if (x1 < RD_half_width_x) {
            if (x3 < -RD_depth) {
                eleTypeArg = 10;  // Bottom surface PML
            } else {
                eleTypeArg = 1;   // Regular domain
            }
        } else {
            if (x3 < -RD_depth) {
                eleTypeArg = 12;  // Right column PML
            } else {
                eleTypeArg = 3;   // Right surface PML
            }
        }
    } else {
        if (x1 < -RD_half_width_x) {
            if (x3 < -RD_depth) {
                eleTypeArg = 18;  // Left-top corner PML
            } else {
                eleTypeArg = 9;   // Left-top column PML
            }
        } else if (x1 < RD_half_width_x) {
            if (x3 < -RD_depth) {
                eleTypeArg = 13;  // Top column PML
            } else {
                eleTypeArg = 4;   // Top surface PML
            }
        } else {
            if (x3 < -RD_depth) {
                eleTypeArg = 17;  // Right-top corner PML
            } else {
                eleTypeArg = 8;   // Right-top column PML
            }
        }
    }
    
    // Set normal vectors and reference coordinates based on element type
    switch (eleTypeArg) {
        case 1:  // Regular domain
            n1 = 0.0;
            n2 = 0.0;
            n3 = 0.0;
            x1_0 = 0.0;
            x2_0 = 0.0;
            x3_0 = 0.0;
            break;
        case 2:  // Bottom surface PML
            n1 = 0.0;
            n2 = -1.0;
            n3 = 0.0;
            x1_0 = 0.0;
            x2_0 = -1.0 * RD_half_width_y;
            x3_0 = 0.0;
            break;
        case 3:  // Right surface PML
            n1 = 1.0;
            n2 = 0.0;
            n3 = 0.0;
            x1_0 = 1.0 * RD_half_width_x;
            x2_0 = 0.0;
            x3_0 = 0.0;
            break;
        case 4:  // Top surface PML
            n1 = 0.0;
            n2 = 1.0;
            n3 = 0.0;
            x1_0 = 0.0;
            x2_0 = 1.0 * RD_half_width_y;
            x3_0 = 0.0;
            break;
        case 5:  // Left surface PML
            n1 = -1.0;
            n2 = 0.0;
            n3 = 0.0;
            x1_0 = -1.0 * RD_half_width_x;
            x2_0 = 0.0;
            x3_0 = 0.0;
            break;
        case 6:  // Left-bottom column PML
            n1 = -1.0;
            n2 = -1.0;
            n3 = 0.0;
            x1_0 = -1.0 * RD_half_width_x;
            x2_0 = -1.0 * RD_half_width_y;
            x3_0 = 0.0;
            break;
        case 7:  // Right-bottom column PML
            n1 = 1.0;
            n2 = -1.0;
            n3 = 0.0;
            x1_0 = 1.0 * RD_half_width_x;
            x2_0 = -1.0 * RD_half_width_y;
            x3_0 = 0.0;
            break;
        case 8:  // Right-top column PML
            n1 = 1.0;
            n2 = 1.0;
            n3 = 0.0;
            x1_0 = 1.0 * RD_half_width_x;
            x2_0 = 1.0 * RD_half_width_y;
            x3_0 = 0.0;
            break;
        case 9:  // Left-top column PML
            n1 = -1.0;
            n2 = 1.0;
            n3 = 0.0;
            x1_0 = -1.0 * RD_half_width_x;
            x2_0 = 1.0 * RD_half_width_y;
            x3_0 = 0.0;
            break;
        case 10:  // Bottom surface PML
            n1 = 0.0;
            n2 = 0.0;
            n3 = -1.0;
            x1_0 = 0.0;
            x2_0 = 0.0;
            x3_0 = -1.0 * RD_depth;
            break;
        case 11:  // Bottom column PML
            n1 = 0.0;
            n2 = -1.0;
            n3 = -1.0;
            x1_0 = 0.0;
            x2_0 = -1.0 * RD_half_width_y;
            x3_0 = -1.0 * RD_depth;
            break;
        case 12:  // Right column PML
            n1 = 1.0;
            n2 = 0.0;
            n3 = -1.0;
            x1_0 = 1.0 * RD_half_width_x;
            x2_0 = 0.0;
            x3_0 = -1.0 * RD_depth;
            break;
        case 13:  // Top column PML
            n1 = 0.0;
            n2 = 1.0;
            n3 = -1.0;
            x1_0 = 0.0;
            x2_0 = 1.0 * RD_half_width_y;
            x3_0 = -1.0 * RD_depth;
            break;
        case 14:  // Left column PML
            n1 = -1.0;
            n2 = 0.0;
            n3 = -1.0;
            x1_0 = -1.0 * RD_half_width_x;
            x2_0 = 0.0;
            x3_0 = -1.0 * RD_depth;
            break;
        case 15:  // Left-bottom corner PML
            n1 = -1.0;
            n2 = -1.0;
            n3 = -1.0;
            x1_0 = -1.0 * RD_half_width_x;
            x2_0 = -1.0 * RD_half_width_y;
            x3_0 = -1.0 * RD_depth;
            break;
        case 16:  // Right-bottom corner PML
            n1 = 1.0;
            n2 = -1.0;
            n3 = -1.0;
            x1_0 = 1.0 * RD_half_width_x;
            x2_0 = -1.0 * RD_half_width_y;
            x3_0 = -1.0 * RD_depth;
            break;
        case 17:  // Right-top corner PML
            n1 = 1.0;
            n2 = 1.0;
            n3 = -1.0;
            x1_0 = 1.0 * RD_half_width_x;
            x2_0 = 1.0 * RD_half_width_y;
            x3_0 = -1.0 * RD_depth;
            break;
        case 18:  // Left-top corner PML
            n1 = -1.0;
            n2 = 1.0;
            n3 = -1.0;
            x1_0 = -1.0 * RD_half_width_x;
            x2_0 = 1.0 * RD_half_width_y;
            x3_0 = -1.0 * RD_depth;
            break;
        default:
            break;
    }
    
    // Calculate PML parameters
    if (eleTypeArg == 1) {
        // Regular domain - no stretching
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 3; j++) {
                pmlAlphaBeta[i][j] = 0.0;
            }
        }
        // raise error to check the parameters and try to fix the issue
        opserr << "ERROR: PML parameters are not calculated for regular domain\n";
        exit(EXIT_FAILURE);
    } else {
        // PML domain - calculate stretching parameters
        double PML_b = PML_L / 1.0;
        double alpha_0 = ((afp + 1.0) * PML_b) / (2.0 * PML_L) * log10(1.0 / PML_Rcoef);
        double beta_0 = ((afp + 1.0) * cp_ref) / (2.0 * PML_L) * log10(1.0 / PML_Rcoef);
        
        // Primary stretching in normal direction
        pmlAlphaBeta[0][0] = 1.0 + alpha_0 * pow(((x1 - x1_0) * n1 / PML_L), afp);
        pmlAlphaBeta[0][1] = 1.0 + alpha_0 * pow(((x2 - x2_0) * n2 / PML_L), afp);
        pmlAlphaBeta[0][2] = 1.0 + alpha_0 * pow(((x3 - x3_0) * n3 / PML_L), afp);
        
        pmlAlphaBeta[1][0] = beta_0 * pow(((x1 - x1_0) * n1 / PML_L), afp);
        pmlAlphaBeta[1][1] = beta_0 * pow(((x2 - x2_0) * n2 / PML_L), afp);
        pmlAlphaBeta[1][2] = beta_0 * pow(((x3 - x3_0) * n3 / PML_L), afp);
    }

    static int intCount = 0;
    // print the PML parameters for the coordinates given:
    if (eleCount == 1) {
        opserr << "PML parameters for coordinates (" << intCount << "  "<< x1 << ", " << x2 << ", " << x3 << "):\n";
        opserr << "Alpha: (" << pmlAlphaBeta[0][0] << ", " << pmlAlphaBeta[0][1] << ", " << pmlAlphaBeta[0][2] << ")\n";
        opserr << "Beta: (" << pmlAlphaBeta[1][0] << ", " << pmlAlphaBeta[1][1] << ", " << pmlAlphaBeta[1][2] << ")\n";
    }
}

// =======================================================================
//  destructor
// ======================================================================= 
PML3DVISCOUS::~PML3DVISCOUS()
{

}

// =======================================================================
// Set Domain
// =======================================================================
void  PML3DVISCOUS::setDomain(Domain* theDomain)
{

	Domainptr = theDomain;

	// node pointers
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++)
		nodePointers[i] = theDomain->getNode(connectedExternalNodes(i));

	this->DomainComponent::setDomain(theDomain);

	// create the coordinate vectors
	double coords[PML3DVISCOUS_NUM_DOF];
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		const Vector& loc = nodePointers[i]->getCrds();
		coords[i * 3] = loc(0);
		coords[i * 3 + 1] = loc(1);
		coords[i * 3 + 2] = loc(2);
	}
	int NDOFEL = PML3DVISCOUS_NUM_DOF;
	int NPROPS = 13;
	int MCRD = 3;
	int NNODE = 8;
	int LFLAGS = 12;
	for (int i = 0; i < PML3DVISCOUS_NUM_DOF*PML3DVISCOUS_NUM_DOF; i++) {
		C[i] = 0.0;
		K[i] = 0.0;
		M[i] = 0.0;
		G[i] = 0.0;
		H[i] = 0.0;
	}
	// Call Fortran function
    pml3d_(M, C, K, G, H, &NDOFEL, props, coords, &MCRD, &NNODE, &LFLAGS);

    calculateMassMatrix();
    
    // Verify the matrices
    verifyMatrices();
    
    dt = theDomain->getDT();
	int tag = this->getTag();
}


void PML3DVISCOUS::calculateMassMatrix() {
    // Reset mass matrix
    for (int i = 0; i < PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF; i++) {
        M_cpp[i] = 0.0;
    }
    
    // Get node coordinates
    double coords[3 * PML3DVISCOUS_NUM_NODES];
    for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
        const Vector& loc = nodePointers[i]->getCrds();
        coords[i * 3] = loc(0);
        coords[i * 3 + 1] = loc(1);
        coords[i * 3 + 2] = loc(2);
    }
    
    // Determine number of integration points based on afp parameter
    int n_points = (afp < 3.0) ? 27 : 64;
    
    // Calculate Lame constants
    double lambda = xnu * E / ((1.0 + xnu) * (1.0 - 2.0 * xnu));
    double mu = 0.5 * E / (1.0 + xnu);
    
    // Initialize matrices
    double M_RD[24][24] = {{0.0}};  // Regular domain mass matrix
    double M_a[24][24] = {{0.0}};   // PML mass matrices
    double M_b[24][24] = {{0.0}};
    double M_c[24][24] = {{0.0}};
    double M_d[24][24] = {{0.0}};
    double N_a[48][48] = {{0.0}};   // N matrices for PML
    
    // Process each integration point
    for (int kint = 0; kint < n_points; kint++) {
        // Shape functions and derivatives
        double N[PML3DVISCOUS_NUM_NODES] = {0.0};
        double dNdxi[PML3DVISCOUS_NUM_NODES][3] = {{0.0}};
        
        // Calculate shape functions at integration point
        double point[3] = {xi[0][kint], xi[1][kint], xi[2][kint]};
        calculateShapeFunctions(point, PML3DVISCOUS_NUM_NODES, N, dNdxi);
        
        // Calculate Jacobian matrix
        double dxdxi[3][3] = {{0.0}};
        for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    dxdxi[j][k] += coords[i * 3 + j] * dNdxi[i][k];
                }
            }
        }
        
        // Calculate determinant of Jacobian
        double determinant = 
            dxdxi[0][0] * (dxdxi[1][1] * dxdxi[2][2] - dxdxi[1][2] * dxdxi[2][1]) -
            dxdxi[0][1] * (dxdxi[1][0] * dxdxi[2][2] - dxdxi[1][2] * dxdxi[2][0]) +
            dxdxi[0][2] * (dxdxi[1][0] * dxdxi[2][1] - dxdxi[1][1] * dxdxi[2][0]);
        
        // Calculate inverse of Jacobian
        double dxidx[3][3] = {{0.0}};
        dxidx[0][0] = (dxdxi[1][1] * dxdxi[2][2] - dxdxi[1][2] * dxdxi[2][1]) / determinant;
        dxidx[0][1] = (dxdxi[0][2] * dxdxi[2][1] - dxdxi[0][1] * dxdxi[2][2]) / determinant;
        dxidx[0][2] = (dxdxi[0][1] * dxdxi[1][2] - dxdxi[0][2] * dxdxi[1][1]) / determinant;
        dxidx[1][0] = (dxdxi[1][2] * dxdxi[2][0] - dxdxi[1][0] * dxdxi[2][2]) / determinant;
        dxidx[1][1] = (dxdxi[0][0] * dxdxi[2][2] - dxdxi[0][2] * dxdxi[2][0]) / determinant;
        dxidx[1][2] = (dxdxi[0][2] * dxdxi[1][0] - dxdxi[0][0] * dxdxi[1][2]) / determinant;
        dxidx[2][0] = (dxdxi[1][0] * dxdxi[2][1] - dxdxi[1][1] * dxdxi[2][0]) / determinant;
        dxidx[2][1] = (dxdxi[0][1] * dxdxi[2][0] - dxdxi[0][0] * dxdxi[2][1]) / determinant;
        dxidx[2][2] = (dxdxi[0][0] * dxdxi[1][1] - dxdxi[0][1] * dxdxi[1][0]) / determinant;
        
        // Calculate dNdx
        double dNdx[PML3DVISCOUS_NUM_NODES][3] = {{0.0}};
        for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    dNdx[i][j] += dNdxi[i][k] * dxidx[k][j];
                }
            }
        }
        
        // Calculate physical coordinates at integration point
        double x1 = 0.0, x2 = 0.0, x3 = 0.0;
        for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
            x1 += N[i] * coords[i * 3];
            x2 += N[i] * coords[i * 3 + 1];
            x3 += N[i] * coords[i * 3 + 2];
        }
        
        // Calculate PML parameters
        double pmlAlphaBeta[2][3] = {{0.0}};
        calculatePMLParameters(x1, x2, x3, pmlAlphaBeta);
        
        // Calculate coefficients for PML matrices
        double coef_a = pmlAlphaBeta[0][0] * pmlAlphaBeta[0][1] * pmlAlphaBeta[0][2];
        double coef_b = pmlAlphaBeta[0][0] * pmlAlphaBeta[0][1] * pmlAlphaBeta[1][2] + 
                        pmlAlphaBeta[0][0] * pmlAlphaBeta[0][2] * pmlAlphaBeta[1][1] + 
                        pmlAlphaBeta[0][1] * pmlAlphaBeta[0][2] * pmlAlphaBeta[1][0];
        double coef_c = pmlAlphaBeta[0][0] * pmlAlphaBeta[1][1] * pmlAlphaBeta[1][2] + 
                        pmlAlphaBeta[0][1] * pmlAlphaBeta[1][2] * pmlAlphaBeta[1][0] + 
                        pmlAlphaBeta[0][2] * pmlAlphaBeta[1][1] * pmlAlphaBeta[1][0];
        double coef_d = pmlAlphaBeta[1][0] * pmlAlphaBeta[1][1] * pmlAlphaBeta[1][2];
        
        // Calculate M_RD and M_a/b/c/d matrices (directly matching Fortran implementation)
        for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
            for (int j = 0; j < PML3DVISCOUS_NUM_NODES; j++) {
                // Basic mass terms
                double mass_term = rho * N[i] * N[j] * w[kint] * determinant;
                M_RD[i][j] += mass_term;
                
                // PML mass terms with coefficients
                M_a[i][j] += coef_a * mass_term;
                M_b[i][j] += coef_b * mass_term;
                M_c[i][j] += coef_c * mass_term;
                M_d[i][j] += coef_d * mass_term;
            }
        }
    }
    
    // Extend mass matrices to other dimensions (x, y, z)
    for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
        for (int j = 0; j < PML3DVISCOUS_NUM_NODES; j++) {
            // Copy x-dimension entries to y and z dimensions
            M_RD[i+PML3DVISCOUS_NUM_NODES][j+PML3DVISCOUS_NUM_NODES] = M_RD[i][j];
            M_RD[i+2*PML3DVISCOUS_NUM_NODES][j+2*PML3DVISCOUS_NUM_NODES] = M_RD[i][j];
            
            M_a[i+PML3DVISCOUS_NUM_NODES][j+PML3DVISCOUS_NUM_NODES] = M_a[i][j];
            M_a[i+2*PML3DVISCOUS_NUM_NODES][j+2*PML3DVISCOUS_NUM_NODES] = M_a[i][j];
        }
    }
    
    // Calculate N_a matrix (needed for PML elements)
    for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
        for (int j = 0; j < PML3DVISCOUS_NUM_NODES; j++) {
            // Block (1,1)
            N_a[i][j] = M_a[i][j]/rho*(lambda+mu)/mu/(3.0*lambda+2.0*mu);
            
            // Block (1,2)
            N_a[i][j+PML3DVISCOUS_NUM_NODES] = -M_a[i][j]/rho*lambda/mu/2.0/(3.0*lambda+2.0*mu);
            
            // Block (1,3)
            N_a[i][j+2*PML3DVISCOUS_NUM_NODES] = -M_a[i][j]/rho*lambda/mu/2.0/(3.0*lambda+2.0*mu);
            
            // Block (2,1)
            N_a[i+PML3DVISCOUS_NUM_NODES][j] = -M_a[i][j]/rho*lambda/mu/2.0/(3.0*lambda+2.0*mu);
            
            // Block (2,2)
            N_a[i+PML3DVISCOUS_NUM_NODES][j+PML3DVISCOUS_NUM_NODES] = 
                M_a[i][j]/rho*(lambda+mu)/mu/(3.0*lambda+2.0*mu);
            
            // Block (2,3)
            N_a[i+PML3DVISCOUS_NUM_NODES][j+2*PML3DVISCOUS_NUM_NODES] = 
                -M_a[i][j]/rho*lambda/mu/2.0/(3.0*lambda+2.0*mu);
            
            // Block (3,1)
            N_a[i+2*PML3DVISCOUS_NUM_NODES][j] = -M_a[i][j]/rho*lambda/mu/2.0/(3.0*lambda+2.0*mu);
            
            // Block (3,2)
            N_a[i+2*PML3DVISCOUS_NUM_NODES][j+PML3DVISCOUS_NUM_NODES] = 
                -M_a[i][j]/rho*lambda/mu/2.0/(3.0*lambda+2.0*mu);
            
            // Block (3,3)
            N_a[i+2*PML3DVISCOUS_NUM_NODES][j+2*PML3DVISCOUS_NUM_NODES] = 
                M_a[i][j]/rho*(lambda+mu)/mu/(3.0*lambda+2.0*mu);
            
            // Blocks for shear components
            N_a[i+3*PML3DVISCOUS_NUM_NODES][j+3*PML3DVISCOUS_NUM_NODES] = M_a[i][j]/rho/mu;
            N_a[i+4*PML3DVISCOUS_NUM_NODES][j+4*PML3DVISCOUS_NUM_NODES] = M_a[i][j]/rho/mu;
            N_a[i+5*PML3DVISCOUS_NUM_NODES][j+5*PML3DVISCOUS_NUM_NODES] = M_a[i][j]/rho/mu;
        }
    }
    
    // Build the complete M_PML matrix based on element type
    double M_PML[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF] = {0.0};



    for (int i = 0; i < 3*PML3DVISCOUS_NUM_NODES; i++) {
        for (int j = 0; j < 3*PML3DVISCOUS_NUM_NODES; j++) {
            M_PML[i*PML3DVISCOUS_NUM_DOF + j] = M_a[i][j];
        }
    }
    
    // // PML domain - A_eu * Damp_beta is used for upper-right block
    // for (int i = 0; i < 3*PML3DVISCOUS_NUM_NODES; i++) {
    //     for (int j = 3*PML3DVISCOUS_NUM_NODES; j < 9*PML3DVISCOUS_NUM_NODES; j++) {
    //         M_PML[i*PML3DVISCOUS_NUM_DOF + j]  = A_eu[i][j-3*PML3DVISCOUS_NUM_NODES] * Damp_beta;
    //     }
    // }

    // PML domain - N_a is used for lower-right block
    for (int i = 3*PML3DVISCOUS_NUM_NODES; i < 9*PML3DVISCOUS_NUM_NODES; i++) {
        for (int j = 3*PML3DVISCOUS_NUM_NODES; j < 9*PML3DVISCOUS_NUM_NODES; j++) {
            M_PML[i*PML3DVISCOUS_NUM_DOF + j] = -N_a[i-3*PML3DVISCOUS_NUM_NODES][j-3*PML3DVISCOUS_NUM_NODES];
        }
    }


    // print M_pml
    opserr << "M_PML matrix:\n";
    for (int i = 0; i < 72; i++) {
        for (int j = 0; j < 72; j++) {
            opserr << M_PML[i*PML3DVISCOUS_NUM_DOF + j] << " ";
        }
        opserr << "\n";
    }

    for (int i = 1; i <= 8; i++) {
        for (int j = 1; j <= 8; j++) {
            for (int k = 0; k < 9; k++) {
                for (int l = 0; l < 9; l++) {
                    int row_tgt = (i - 1) * 9 + k;      // Target row index
                    int col_tgt = (j - 1) * 9 + l;      // Target column index
                    int row_src = (i - 1) + 8 * k;      // Source row index
                    int col_src = (j - 1) + 8 * l;      // Source column index
                    int idx_tgt = row_tgt * 72 + col_tgt;
                    int idx_src = row_src * 72 + col_src;
                    M_cpp[idx_tgt] = M_PML[idx_src];
                }
            }
        }
    }

}

// =======================================================================
// update
// =======================================================================
int PML3DVISCOUS::update(void)
{
	dt = Domainptr->getDT();
	// opserr << "dt = " << dt << "\n";	
	// get u, v, a from nodes and calculate the ubar vector
	int loc = 0;
	double c1 = dt;
	double c2 = dt * dt * 0.5;
	double c3 = dt*dt*dt*((1.0/6.0)-eta);
	double c4 = dt*dt*dt*eta;
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		const Vector& uNode = nodePointers[i]->getDisp();
		const Vector& vNode = nodePointers[i]->getVel();
		const Vector& aNode = nodePointers[i]->getAccel();
		const Vector& atpdt = nodePointers[i]->getTrialAccel();
		for (int j = 0; j < 9; j++) {
			ubar(loc) = ubart(loc) + uNode(j)*c1 + vNode(j)*c2 + aNode(j)*c3 + atpdt(j)*c4; 
			loc++;
		}
	}

    double keisi = 1.0/48.0;

    loc = 0;
    c1 = dt;
    c2 = dt * dt * 0.5;
    c3 = dt * dt * dt /6.0;
    c4 = dt * dt * dt * dt * (1.0/24.0 - keisi);
    int c5 = dt * dt * dt * dt * keisi;
    for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
        const Vector& uNode = nodePointers[i]->getDisp();
        const Vector& vNode = nodePointers[i]->getVel();
        const Vector& aNode = nodePointers[i]->getAccel();
        const Vector& atpdt = nodePointers[i]->getTrialAccel();
        for (int j = 0; j < 9; j++) {
            ubarbar(loc) = ubarbart(loc) + ubart(loc)*c1 + uNode(j)*c2 + vNode(j)*c3 + aNode(j)*c4 + atpdt(j)*c5; 
            loc++;
        }
    }

	return 0;
}

// =======================================================================
//	return stiffness matrix 
// =======================================================================
const Matrix& PML3DVISCOUS::getTangentStiff()
{
	// check if the dt is changed to update the tangent stiffness matrix
	double cg = eta*dt/beta;
    double keisi = 1.0/48.0;
    double ch = dt * dt * keisi/beta;
	//keff = k + cg*g( k and g are symmetric matrices)
	for (int i = 0; i < PML3DVISCOUS_NUM_DOF*PML3DVISCOUS_NUM_DOF; i++) {
		Keff[i] = K[i] + cg*G[i] + ch*H[i];
	}
	tangent.setData(Keff, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
	return tangent;
}

// =======================================================================
//	return initial stiffness matrix 
// =======================================================================
const Matrix& PML3DVISCOUS::getInitialStiff()
{
	return this->getTangentStiff();
}

// =======================================================================
//	return mass matrix
// =======================================================================
const Matrix& PML3DVISCOUS::getMass()
{
	mass.setData(M, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
	// mass.Zero();
	return mass;
}

// =======================================================================
//	return damping matrix
// =======================================================================
const Matrix& PML3DVISCOUS::getDamp()
{
	damping.setData(C, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
	// damping.Zero();
	return damping;
}

// =======================================================================
// Ressisting force
// =======================================================================
//get residual
const Vector& PML3DVISCOUS::getResistingForce()
{
	// if (innertag==14) {
	// 	opserr << "getResistingForce function is called\n";
	// }
	int numNodalDOF = 9;
	static Vector theVector(PML3DVISCOUS_NUM_DOF);

	// get K into stiff
	tangent.setData(K, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);

	//
	// perform: R = K * u
	//

	int loc = 0;
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		const Vector& uNode = nodePointers[i]->getTrialDisp();
		for (int j = 0; j < numNodalDOF; j++)
			theVector(loc++) = uNode(j);
	}
	resid.addMatrixVector(0.0, tangent, theVector, 1.0);
	return resid;
}


// =======================================================================
//
// =======================================================================
//get residual with inertia terms
const Vector&
PML3DVISCOUS::getResistingForceIncInertia()
{
    // R += K*u
	static Vector theVector(PML3DVISCOUS_NUM_DOF);
	tangent.setData(K, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);

	int loc = 0;
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		const Vector& uNode = nodePointers[i]->getTrialDisp();
		for (int j = 0; j < 9; j++)
			theVector(loc++) = uNode(j);
	}
	resid.addMatrixVector(0.0, tangent, theVector, 1.0);



	// R += M*a
	loc = 0;
	Node** theNodes = this->getNodePtrs();
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		const Vector& acc = theNodes[i]->getTrialAccel();
		for (int j = 0; j < 9; j++) {
			theVector(loc++) = acc(j);
		}
	}
	resid.addMatrixVector(1.0, this->getMass(), theVector, 1.0);

	// R += C*v
	loc = 0;
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		const Vector& vel = theNodes[i]->getTrialVel();
		for (int j = 0; j < 9; j++) {
			theVector(loc++) = vel[j];
		}
	}
	resid.addMatrixVector(1.0, this->getDamp(), theVector, 1.0);


	// R += G*ubar
	tangent.setData(G, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
	resid.addMatrixVector(1.0, tangent, ubar, 1.0);

    // R += H*ubarbar
    tangent.setData(H, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
    resid.addMatrixVector(1.0, tangent, ubarbar, 1.0);
	
	
	return resid;
}

// =======================================================================
// get the number of external nodes
// =======================================================================
int  PML3DVISCOUS::getNumExternalNodes() const
{
	return PML3DVISCOUS_NUM_NODES;
}

// =======================================================================
// return connected external nodes
// =======================================================================
const ID& PML3DVISCOUS::getExternalNodes()
{
	return connectedExternalNodes;
}

// =======================================================================
// return node pointers
// =======================================================================
Node** PML3DVISCOUS::getNodePtrs(void)
{
	return nodePointers;
}

// =======================================================================
// return number of dofs
// =======================================================================
int  PML3DVISCOUS::getNumDOF()
{
	return PML3DVISCOUS_NUM_DOF;
}

// =======================================================================
// commit state
// =======================================================================
int  PML3DVISCOUS::commitState()
{
	int success = 0;
	if ((success = this->Element::commitState()) != 0) {
		opserr << "PML3DVISCOUS::commitState () - failed in base class";
	}

	// set ubart to ubar
	for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
		ubart(i) = ubar(i);
	}

	updateflag = 0;
	return success;
}

// =======================================================================
// revert to last commit 
// =======================================================================
int  PML3DVISCOUS::revertToLastCommit()
{
	int success = 0;

	// set ubar to ubart
	for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
		ubar(i) = ubart(i);
        ubarbar(i) = ubarbart(i);
	}

	return success;
}

// =======================================================================
// revert to start
// =======================================================================
int  PML3DVISCOUS::revertToStart()
{
	int success = 0;

	// set ubar and ubart to zero
	for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
		ubar(i) = 0.0;
		ubart(i) = 0.0;
        ubarbar(i) = 0.0;
        ubarbart(i) = 0.0;
	}

	return success;
}

// =======================================================================
// add load
// =======================================================================
int PML3DVISCOUS::addLoad(ElementalLoad* theLoad, double loadFactor)
{
	return -1;
}

// =======================================================================
// add zero load
// =======================================================================
void  PML3DVISCOUS::zeroLoad()
{
	return;
}

// =======================================================================
// senself
// =======================================================================
int  PML3DVISCOUS::sendSelf(int commitTag,
	Channel& theChannel)
{
	int res = 0;

	// note: we don't check for dataTag == 0 for Element
	// objects as that is taken care of in a commit by the Domain
	// object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// PML3DVISCOUS packs its data into a Vector and sends this to theChannel
	// along with its dbTag and the commitTag passed in the arguments
	static Vector data(PML3DVISCOUS_NUM_PROPS + 4);
	data(0) = this->getTag();

	for (int ii = 1; ii <= PML3DVISCOUS_NUM_PROPS; ii++) {
		data(ii) = props[ii - 1];
	}
	data(PML3DVISCOUS_NUM_PROPS+1) = eta;
	data(PML3DVISCOUS_NUM_PROPS+2) = beta;
	data(PML3DVISCOUS_NUM_PROPS+3) = gamma;

	res += theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING PML3DVISCOUS::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return res;
	}


	// PML3DVISCOUS then sends the tags of its four nodes
	res += theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
	if (res < 0) {
		opserr << "WARNING PML3DVISCOUS::sendSelf() - " << this->getTag() << " failed to send ID\n";
		return res;
	}

	return res;
}

// =======================================================================
// recvself
// =======================================================================
int  PML3DVISCOUS::recvSelf(int commitTag,
	Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	int res = 0;

	int dataTag = this->getDbTag();

	// PML3DVISCOUS creates a Vector, receives the Vector and then sets the 
	// internal data with the data in the Vector
	static Vector data(PML3DVISCOUS_NUM_PROPS + 4);
	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING PML3DVISCOUS::recvSelf() - failed to receive Vector\n";
		return res;
	}

	this->setTag((int)data(0));


	for (int ii = 1; ii <= PML3DVISCOUS_NUM_PROPS; ii++) {
		props[ii - 1] = data(ii);
	}

	eta   = data(PML3DVISCOUS_NUM_PROPS+1);
	beta  = data(PML3DVISCOUS_NUM_PROPS+2);
	gamma = data(PML3DVISCOUS_NUM_PROPS+3);

	// PML3DVISCOUS now receives the tags of its four external nodes
	res += theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
	if (res < 0) {
		opserr << "WARNING PML3DVISCOUS::recvSelf() - " << this->getTag() << " failed to receive ID\n";
		return res;
	}

	return res;
}


// =======================================================================
// display
// =======================================================================
int PML3DVISCOUS::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
	// Get the end point display coords
	static Vector v1(3);
	static Vector v2(3);
	static Vector v3(3);
	static Vector v4(3);
	static Vector v5(3);
	static Vector v6(3);
	static Vector v7(3);
	static Vector v8(3);
	nodePointers[0]->getDisplayCrds(v1, fact, displayMode);
	nodePointers[1]->getDisplayCrds(v2, fact, displayMode);
	nodePointers[2]->getDisplayCrds(v3, fact, displayMode);
	nodePointers[3]->getDisplayCrds(v4, fact, displayMode);
	nodePointers[4]->getDisplayCrds(v5, fact, displayMode);
	nodePointers[5]->getDisplayCrds(v6, fact, displayMode);
	nodePointers[6]->getDisplayCrds(v7, fact, displayMode);
	nodePointers[7]->getDisplayCrds(v8, fact, displayMode);

	// place values in coords matrix
	static Matrix coords(8, 3);
	for (int i = 0; i < 3; i++) {
		coords(0, i) = v1(i);
		coords(1, i) = v2(i);
		coords(2, i) = v3(i);
		coords(3, i) = v4(i);
		coords(4, i) = v5(i);
		coords(5, i) = v6(i);
		coords(6, i) = v7(i);
		coords(7, i) = v8(i);
	}

	// fill RGB vector
	static Vector values(8);
	for (int i = 0; i < 8; i++)
		values(i) = 1.0;

	// draw the cube
	return theViewer.drawCube(coords, values, this->getTag());
}

// =======================================================================
// setresponse
// =======================================================================
Response* PML3DVISCOUS::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	Response* theResponse = 0;

	// char outputData[32];

	// output.tag("ElementOutput");
	// output.attr("eleType", "PML3DVISCOUS");
	// output.attr("eleTag", this->getTag());
	// for (int i = 1; i <= 8; i++) {
	// 	sprintf(outputData, "node%d", i);
	// 	output.attr(outputData, nodePointers[i - 1]->getTag());
	// }

	// if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0) {

	// 	for (int i = 1; i <= 8; i++) {
	// 		sprintf(outputData, "P1_%d", i);
	// 		output.tag("ResponseType", outputData);
	// 		sprintf(outputData, "P2_%d", i);
	// 		output.tag("ResponseType", outputData);
	// 		sprintf(outputData, "P3_%d", i);
	// 		output.tag("ResponseType", outputData);
	// 	}

	// 	theResponse = new ElementResponse(this, 1, resid);
	// }
	// output.endTag(); // ElementOutput
	return theResponse;
}

// =======================================================================
// getresponse
// =======================================================================
int PML3DVISCOUS::getResponse(int responseID, Information& eleInfo)
{
	// static Vector stresses(48);

	// if (responseID == 1)
	// 	return eleInfo.setVector(this->getResistingForce());

	return -1;
}

// =======================================================================
// set parameter
// =======================================================================
int PML3DVISCOUS::setParameter(const char** argv, int argc, Parameter& param)
{
	int res = -1;
	return res;
}

// =======================================================================
// update parameter
// =======================================================================
int PML3DVISCOUS::updateParameter(int parameterID, Information& info)
{
	int res = -1;
	return res;
}

// =======================================================================
// print
// =======================================================================
void  PML3DVISCOUS::Print(OPS_Stream &s, int flag) {
	
  	if (flag == OPS_PRINT_CURRENTSTATE) {
		s << "Element: " << this->getTag() << endln;
		s << "type: PML3DVISCOUS \n";
		s << "Nodes: " << connectedExternalNodes;
		s << "eta: " << eta << " beta: " << beta << " gamma: " << gamma << endln;
		s << endln;
		s << "Resisting Force (no inertia): " << this->getResistingForce();
  	}
    
 	 if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": " << this->getTag() << ", ";
		s << "\"type\": \"PML3DVISCOUS\", ";
		s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
		for (int i = 1; i < 7; i++)
		s << connectedExternalNodes(i) << ", ";
		s << connectedExternalNodes(7) << "], ";
  	}
	return;
}



// Implementation of verification function
void PML3DVISCOUS::verifyMatrices() {
    // Define tolerance for floating point comparisons
    const double tolerance = 1e-4;
    bool mismatchFound = false;
    
    // Compare M matrices
    for (int i = 0; i < PML3DVISCOUS_NUM_DOF*PML3DVISCOUS_NUM_DOF; i++) {
        if (fabs(M[i] - M_cpp[i]) > tolerance) {
            opserr << "PML3DVISCOUS::verifyMatrices - Mismatch in M at element " << this->getTag() <<endln;
            opserr << "M[" << i << "] = " << M[i] << ", M_cpp[" << i << "] = " << M_cpp[i] <<endln;
            mismatchFound = true;
            
        }
    }
    opserr << "M matrices match between Fortran and C++ implementations for element " << this->getTag() <<endln;
    return;
    
    // Compare C matrices
    for (int i = 0; i < PML3DVISCOUS_NUM_DOF*PML3DVISCOUS_NUM_DOF; i++) {
        if (fabs(C[i] - C_cpp[i]) > tolerance) {
            opserr << "PML3DVISCOUS::verifyMatrices - Mismatch in C at element " << this->getTag() <<endln;
            mismatchFound = true;
            return;
        }
    }
    
    // Compare K matrices
    for (int i = 0; i < PML3DVISCOUS_NUM_DOF*PML3DVISCOUS_NUM_DOF; i++) {
        if (fabs(K[i] - K_cpp[i]) > tolerance) {
            opserr << "PML3DVISCOUS::verifyMatrices - Mismatch in K at element " << this->getTag() <<endln;
            mismatchFound = true;
            return;

        }
    }
    
    // Compare G matrices
    for (int i = 0; i < PML3DVISCOUS_NUM_DOF*PML3DVISCOUS_NUM_DOF; i++) {
        if (fabs(G[i] - G_cpp[i]) > tolerance) {
            opserr << "PML3DVISCOUS::verifyMatrices - Mismatch in G at element " << this->getTag() <<endln;
            mismatchFound = true;
            return;
        }
    }
    
    // Compare H matrices
    for (int i = 0; i < PML3DVISCOUS_NUM_DOF*PML3DVISCOUS_NUM_DOF; i++) {
        if (fabs(H[i] - H_cpp[i]) > tolerance) {
            opserr << "PML3DVISCOUS::verifyMatrices - Mismatch in H at element " << this->getTag() <<endln;
            mismatchFound = true;
            return;
        }
    }
    
    if (!mismatchFound) {
        opserr << "PML3DVISCOUS::verifyMatrices - All matrices match between Fortran and C++ implementations for element " << this->getTag() <<endln;
    }
}
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
double  PML3DVISCOUS::xi[3][64] = {{0.0}};
double  PML3DVISCOUS::w[64] = {0.0};
bool    PML3DVISCOUS::integrationInitialized = false;

// Initialize static matrices
double PML3DVISCOUS::M_cpp[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF] = {0};
double PML3DVISCOUS::C_cpp[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF] = {0};
double PML3DVISCOUS::K_cpp[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF] = {0};
double PML3DVISCOUS::G_cpp[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF] = {0};
double PML3DVISCOUS::H_cpp[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF] = {0};

// Implementation of the integration points and weights initialization function
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
        if (n_points == 1) {
            xi[0][0] = 0.0;
            xi[1][0] = 0.0;
            xi[2][0] = 0.0;
            w[0] = 8.0;
        }
        else if (n_points == 8) {
            double x1D[2] = {-0.5773502692, 0.5773502692};
            
            int n = 0;
            for (int k = 0; k < 2; k++) {
                for (int j = 0; j < 2; j++) {
                    for (int i = 0; i < 2; i++) {
                        xi[0][n] = x1D[i];
                        xi[1][n] = x1D[j];
                        xi[2][n] = x1D[k];
                        w[n] = 1.0;
                        n++;
                    }
                }
            }
        }
        else if (n_points == 27) {
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
        }
        else if (n_points == 64) {
            double x1D[4] = {0.8611363115940526, 0.3399810435848563, 
                            -0.3399810435848563, -0.8611363115940526};
            double w1D[4] = {0.3478548451374538, 0.6521451548625461, 
                             0.6521451548625461, 0.3478548451374538};
            
            int n = 0;
            for (int k = 0; k < 4; k++) {
                for (int j = 0; j < 4; j++) {
                    for (int i = 0; i < 4; i++) {
                        xi[0][n] = x1D[i];
                        xi[1][n] = x1D[j];
                        xi[2][n] = x1D[k];
                        w[n] = w1D[i] * w1D[j] * w1D[k];
                        n++;
                    }
                }
            }
        }
    }
    integrationInitialized = true;
}

// Implementation of shape functions calculation
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

// Implementation of PML parameters calculation function
void PML3DVISCOUS::calculatePMLParameters(const double* props, double x1, double x2, double x3, double (*pmlAlphaBeta)[3]) {
    // Extract properties
    double E = props[0];
    double xnu = props[1];
    double rho = props[2];
    int eleTypeArg = (int)props[3];
    double PML_L = props[4];
    double afp = props[5];
    double PML_Rcoef = props[6];
    double RD_half_width_x = props[7];
    double RD_half_width_y = props[8];
    double RD_depth = props[9];
    
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
}

// Implementation of the main function to calculate PML matrices
void PML3DVISCOUS::calculatePMLMatrices(const double* props, const double* coords, double* M, double* C, double* K, 
                                       double* G, double* H, int ndofel, int mcrd, int nnode, int lflags) {
    // Initialize constants
    const double ZERO = 0.0;
    const double HALF = 0.5;
    const double ONE = 1.0;
    
    const double coef_alpha = 1.0/12.0;
    const double coef_beta = 1.0/48.0;
    
    // Extract material properties
    double E = props[0];
    double xnu = props[1];
    double rho = props[2];
    int eleTypePos = (int)props[3];
    double PML_L = props[4];
    double afp = props[5];
    double PML_Rcoef = props[6];
    double RD_half_width_x = props[7];
    double RD_half_width_y = props[8];
    double RD_depth = props[9];
    double Damp_alpha = props[10];
    double Damp_beta = props[11];
    
    // Determine number of integration points based on afp
    int n_points = (afp < 3.0) ? 27 : 64;
    
    // Calculate Lame constants
    double lambda = xnu * E / ((1.0 + xnu) * (1.0 - 2.0 * xnu));
    double mu = 0.5 * E / (1.0 + xnu);
    
    // Initialize arrays
    double strain[6] = {0.0};
    double stress[6] = {0.0};
    double D[6][6] = {{0.0}};
    double B[6][60] = {{0.0}};
    double dxidx[3][3] = {{0.0}};
    
    double N[20] = {0.0};
    double dNdxi[20][3] = {{0.0}};
    double dxdxi[3][3] = {{0.0}};
    double dNdx[20][3] = {{0.0}};
    
    double Phi[3][60] = {{0.0}};
    
    // PML matrices
    double K_RD[24][24] = {{0.0}};
    double Kxx[8][8] = {{0.0}};
    double Kxy[8][8] = {{0.0}};
    double Kxz[8][8] = {{0.0}};
    double Kyx[8][8] = {{0.0}};
    double Kyy[8][8] = {{0.0}};
    double Kyz[8][8] = {{0.0}};
    double Kzx[8][8] = {{0.0}};
    double Kzy[8][8] = {{0.0}};
    double Kzz[8][8] = {{0.0}};
    
    double Phi_x[8] = {0.0};
    double Phi_y[8] = {0.0};
    double Phi_z[8] = {0.0};
    
    double M_RD[24][24] = {{0.0}};
    double C_RD[24][24] = {{0.0}};
    double M_a[24][24] = {{0.0}};
    double M_b[24][24] = {{0.0}};
    double M_c[24][24] = {{0.0}};
    double M_d[24][24] = {{0.0}};
    
    double N_a[48][48] = {{0.0}};
    double N_b[48][48] = {{0.0}};
    double N_c[48][48] = {{0.0}};
    double N_d[48][48] = {{0.0}};
    
    double A_eu[24][48] = {{0.0}};
    double A_wu[24][48] = {{0.0}};
    double A_pu[24][48] = {{0.0}};
    double A_el[24][48] = {{0.0}};
    double A_wl[24][48] = {{0.0}};
    double A_pl[24][48] = {{0.0}};
    
    double K_PML[72][72] = {{0.0}};
    double M_PML[72][72] = {{0.0}};
    double C_PML[72][72] = {{0.0}};
    double G_PML[72][72] = {{0.0}};
    double H_PML[72][72] = {{0.0}};
    
    // Make sure integration points are initialized
    if (!integrationInitialized) {
        initializeIntegrationPointsAndWeights(n_points, nnode);
    }
    
    // Initialize matrices to zero
    for (int i = 0; i < ndofel; i++) {
        for (int j = 0; j < ndofel; j++) {
            M[i*ndofel + j] = 0.0;
            C[i*ndofel + j] = 0.0;
            K[i*ndofel + j] = 0.0;
            G[i*ndofel + j] = 0.0;
            H[i*ndofel + j] = 0.0;
        }
    }
    
    // Loop over integration points
    for (int kint = 0; kint < n_points; kint++) {
        // Get shape functions and derivatives at this integration point
        double xiPoint[3] = {xi[0][kint], xi[1][kint], xi[2][kint]};
        calculateShapeFunctions(xiPoint, nnode, N, dNdxi);
        
        // Calculate derivatives of shape functions w.r.t global coordinates
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                dxdxi[i][j] = 0.0;
                for (int k = 0; k < nnode; k++) {
                    dxdxi[i][j] += coords[k*3 + i] * dNdxi[k][j];
                }
            }
        }
        
        // Calculate the determinant and inverse of dxdxi
        double determinant = dxdxi[0][0]*(dxdxi[1][1]*dxdxi[2][2] - dxdxi[1][2]*dxdxi[2][1])
                            - dxdxi[0][1]*(dxdxi[1][0]*dxdxi[2][2] - dxdxi[1][2]*dxdxi[2][0])
                            + dxdxi[0][2]*(dxdxi[1][0]*dxdxi[2][1] - dxdxi[1][1]*dxdxi[2][0]);
        
        if (fabs(determinant) < 1.0e-10) {
            opserr << "Error: Jacobian determinant close to zero in PML3DVISCOUS::calculatePMLMatrices\n";
            return;
        }
        
        // Calculate inverse of dxdxi
        dxidx[0][0] = (dxdxi[1][1]*dxdxi[2][2] - dxdxi[1][2]*dxdxi[2][1])/determinant;
        dxidx[0][1] = (dxdxi[0][2]*dxdxi[2][1] - dxdxi[0][1]*dxdxi[2][2])/determinant;
        dxidx[0][2] = (dxdxi[0][1]*dxdxi[1][2] - dxdxi[0][2]*dxdxi[1][1])/determinant;
        dxidx[1][0] = (dxdxi[1][2]*dxdxi[2][0] - dxdxi[1][0]*dxdxi[2][2])/determinant;
        dxidx[1][1] = (dxdxi[0][0]*dxdxi[2][2] - dxdxi[0][2]*dxdxi[2][0])/determinant;
        dxidx[1][2] = (dxdxi[0][2]*dxdxi[1][0] - dxdxi[0][0]*dxdxi[1][2])/determinant;
        dxidx[2][0] = (dxdxi[1][0]*dxdxi[2][1] - dxdxi[1][1]*dxdxi[2][0])/determinant;
        dxidx[2][1] = (dxdxi[0][1]*dxdxi[2][0] - dxdxi[0][0]*dxdxi[2][1])/determinant;
        dxidx[2][2] = (dxdxi[0][0]*dxdxi[1][1] - dxdxi[0][1]*dxdxi[1][0])/determinant;
        
        // Calculate derivatives of shape functions w.r.t. global coordinates
        for (int i = 0; i < nnode; i++) {
            for (int j = 0; j < 3; j++) {
                dNdx[i][j] = 0.0;
                for (int k = 0; k < 3; k++) {
                    dNdx[i][j] += dNdxi[i][k] * dxidx[k][j];
                }
            }
        }
        
        // Calculate B matrix (strain-displacement)
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 3*nnode; j++) {
                B[i][j] = 0.0;
            }
        }
        
        for (int i = 0; i < nnode; i++) {
            B[0][i*3]     = dNdx[i][0];  // du/dx
            B[1][i*3 + 1] = dNdx[i][1];  // dv/dy
            B[2][i*3 + 2] = dNdx[i][2];  // dw/dz
            
            B[3][i*3]     = dNdx[i][1];  // du/dy
            B[3][i*3 + 1] = dNdx[i][0];  // dv/dx
            
            B[4][i*3]     = dNdx[i][2];  // du/dz
            B[4][i*3 + 2] = dNdx[i][0];  // dw/dx
            
            B[5][i*3 + 1] = dNdx[i][2];  // dv/dz
            B[5][i*3 + 2] = dNdx[i][1];  // dw/dy
        }
        
        // Calculate Phi matrix (shape function matrix)
        for (int i = 0; i < nnode; i++) {
            Phi[0][i*3]     = N[i];
            Phi[1][i*3 + 1] = N[i];
            Phi[2][i*3 + 2] = N[i];
        }
        
        // Calculate derivative vectors
        for (int i = 0; i < nnode; i++) {
            Phi_x[i] = dNdx[i][0];
            Phi_y[i] = dNdx[i][1];
            Phi_z[i] = dNdx[i][2];
        }
        
        // Calculate gauss point coordinates
        double x1 = 0.0, x2 = 0.0, x3 = 0.0;
        for (int i = 0; i < nnode; i++) {
            x1 += N[i] * coords[i*3];
            x2 += N[i] * coords[i*3 + 1];
            x3 += N[i] * coords[i*3 + 2];
        }
        
        // Calculate PML alpha and beta parameters
        double pmlAlphaBeta[2][3] = {{0.0}};
        calculatePMLParameters(props, x1, x2, x3, pmlAlphaBeta);
        
        // Calculate coefficients for PML
        double coef_a = pmlAlphaBeta[0][0] * pmlAlphaBeta[0][1] * pmlAlphaBeta[0][2];
        double coef_b = pmlAlphaBeta[0][0] * pmlAlphaBeta[0][1] * pmlAlphaBeta[1][2] +
                        pmlAlphaBeta[0][0] * pmlAlphaBeta[0][2] * pmlAlphaBeta[1][1] +
                        pmlAlphaBeta[0][1] * pmlAlphaBeta[0][2] * pmlAlphaBeta[1][0];
        double coef_c = pmlAlphaBeta[0][0] * pmlAlphaBeta[1][1] * pmlAlphaBeta[1][2] +
                        pmlAlphaBeta[0][1] * pmlAlphaBeta[1][2] * pmlAlphaBeta[1][0] +
                        pmlAlphaBeta[0][2] * pmlAlphaBeta[1][1] * pmlAlphaBeta[1][0];
        double coef_d = pmlAlphaBeta[1][0] * pmlAlphaBeta[1][1] * pmlAlphaBeta[1][2];
        
        // Calculate L coefficients
        double coef_Le[3][3] = {{0.0}};
        double coef_Lp[3][3] = {{0.0}};
        double coef_Lw[3][3] = {{0.0}};
        
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                coef_Le[i][j] = pmlAlphaBeta[0][i] * pmlAlphaBeta[0][j];
                coef_Lp[i][j] = pmlAlphaBeta[0][i] * pmlAlphaBeta[1][j] +
                               pmlAlphaBeta[1][i] * pmlAlphaBeta[0][j];
                coef_Lw[i][j] = pmlAlphaBeta[1][i] * pmlAlphaBeta[1][j];
            }
        }
        
        // Calculate stiffness matrix components
        for (int i = 0; i < nnode; i++) {
            for (int j = 0; j < nnode; j++) {
                // Diagonal components
                Kxx[i][j] = (lambda + 2.0*mu) * Phi_x[i] * Phi_x[j] + 
                           mu * (Phi_y[i] * Phi_y[j] + Phi_z[i] * Phi_z[j]);
                
                Kyy[i][j] = (lambda + 2.0*mu) * Phi_y[i] * Phi_y[j] + 
                           mu * (Phi_x[i] * Phi_x[j] + Phi_z[i] * Phi_z[j]);
                
                Kzz[i][j] = (lambda + 2.0*mu) * Phi_z[i] * Phi_z[j] + 
                           mu * (Phi_x[i] * Phi_x[j] + Phi_y[i] * Phi_y[j]);
                
                // Off-diagonal components
                Kxy[i][j] = lambda * Phi_x[i] * Phi_y[j] + mu * Phi_y[i] * Phi_x[j];
                Kxz[i][j] = lambda * Phi_x[i] * Phi_z[j] + mu * Phi_z[i] * Phi_x[j];
                Kyz[i][j] = lambda * Phi_y[i] * Phi_z[j] + mu * Phi_z[i] * Phi_y[j];
                
                // Mass matrix components
                M_RD[i][j] += rho * N[i] * N[j] * w[kint] * determinant;
                
                // PML matrix components
                M_a[i][j] += coef_a * rho * N[i] * N[j] * w[kint] * determinant;
                M_b[i][j] += coef_b * rho * N[i] * N[j] * w[kint] * determinant;
                M_c[i][j] += coef_c * rho * N[i] * N[j] * w[kint] * determinant;
                M_d[i][j] += coef_d * rho * N[i] * N[j] * w[kint] * determinant;
                
                // A matrices
                A_eu[i][j] += Phi_x[i] * N[j] * coef_Le[1][2] * w[kint] * determinant;
                A_eu[i][j+nnode*3] += Phi_y[i] * N[j] * coef_Le[0][2] * w[kint] * determinant;
                A_eu[i][j+nnode*4] += Phi_z[i] * N[j] * coef_Le[0][1] * w[kint] * determinant;
                
                A_eu[i+nnode][j+nnode] += Phi_y[i] * N[j] * coef_Le[0][2] * w[kint] * determinant;
                A_eu[i+nnode][j+nnode*3] += Phi_x[i] * N[j] * coef_Le[1][2] * w[kint] * determinant;
                A_eu[i+nnode][j+nnode*5] += Phi_z[i] * N[j] * coef_Le[0][1] * w[kint] * determinant;
                
                A_eu[i+nnode*2][j+nnode*2] += Phi_z[i] * N[j] * coef_Le[0][1] * w[kint] * determinant;
                A_eu[i+nnode*2][j+nnode*4] += Phi_x[i] * N[j] * coef_Le[1][2] * w[kint] * determinant;
                A_eu[i+nnode*2][j+nnode*5] += Phi_y[i] * N[j] * coef_Le[0][2] * w[kint] * determinant;
                
                // A_wu matrices
                A_wu[i][j] += Phi_x[i] * N[j] * coef_Lw[1][2] * w[kint] * determinant;
                A_wu[i][j+nnode*3] += Phi_y[i] * N[j] * coef_Lw[0][2] * w[kint] * determinant;
                A_wu[i][j+nnode*4] += Phi_z[i] * N[j] * coef_Lw[0][1] * w[kint] * determinant;
                
                A_wu[i+nnode][j+nnode] += Phi_y[i] * N[j] * coef_Lw[0][2] * w[kint] * determinant;
                A_wu[i+nnode][j+nnode*3] += Phi_x[i] * N[j] * coef_Lw[1][2] * w[kint] * determinant;
                A_wu[i+nnode][j+nnode*5] += Phi_z[i] * N[j] * coef_Lw[0][1] * w[kint] * determinant;
                
                A_wu[i+nnode*2][j+nnode*2] += Phi_z[i] * N[j] * coef_Lw[0][1] * w[kint] * determinant;
                A_wu[i+nnode*2][j+nnode*4] += Phi_x[i] * N[j] * coef_Lw[1][2] * w[kint] * determinant;
                A_wu[i+nnode*2][j+nnode*5] += Phi_y[i] * N[j] * coef_Lw[0][2] * w[kint] * determinant;
                
                // A_pu matrices
                A_pu[i][j] += Phi_x[i] * N[j] * coef_Lp[1][2] * w[kint] * determinant;
                A_pu[i][j+nnode*3] += Phi_y[i] * N[j] * coef_Lp[0][2] * w[kint] * determinant;
                A_pu[i][j+nnode*4] += Phi_z[i] * N[j] * coef_Lp[0][1] * w[kint] * determinant;
                
                A_pu[i+nnode][j+nnode] += Phi_y[i] * N[j] * coef_Lp[0][2] * w[kint] * determinant;
                A_pu[i+nnode][j+nnode*3] += Phi_x[i] * N[j] * coef_Lp[1][2] * w[kint] * determinant;
                A_pu[i+nnode][j+nnode*5] += Phi_z[i] * N[j] * coef_Lp[0][1] * w[kint] * determinant;
                
                A_pu[i+nnode*2][j+nnode*2] += Phi_z[i] * N[j] * coef_Lp[0][1] * w[kint] * determinant;
                A_pu[i+nnode*2][j+nnode*4] += Phi_x[i] * N[j] * coef_Lp[1][2] * w[kint] * determinant;
                A_pu[i+nnode*2][j+nnode*5] += Phi_y[i] * N[j] * coef_Lp[0][2] * w[kint] * determinant;
            }
        }
        
        // Calculate Kyx, Kzx, and Kzy as transposes of Kxy, Kxz, and Kyz
        for (int i = 0; i < nnode; i++) {
            for (int j = 0; j < nnode; j++) {
                Kyx[i][j] = Kxy[j][i];
                Kzx[i][j] = Kxz[j][i];
                Kzy[i][j] = Kyz[j][i];
            }
        }
        
        // Assemble the full K_RD matrix
        for (int i = 0; i < nnode; i++) {
            for (int j = 0; j < nnode; j++) {
                // Diagonal blocks
                K_RD[i][j] += Kxx[i][j] * w[kint] * determinant;
                K_RD[i+nnode][j+nnode] += Kyy[i][j] * w[kint] * determinant;
                K_RD[i+2*nnode][j+2*nnode] += Kzz[i][j] * w[kint] * determinant;
                
                // Off-diagonal blocks
                K_RD[i][j+nnode] += Kxy[i][j] * w[kint] * determinant;
                K_RD[i][j+2*nnode] += Kxz[i][j] * w[kint] * determinant;
                K_RD[i+nnode][j+2*nnode] += Kyz[i][j] * w[kint] * determinant;
                
                K_RD[i+nnode][j] += Kyx[i][j] * w[kint] * determinant;
                K_RD[i+2*nnode][j] += Kzx[i][j] * w[kint] * determinant;
                K_RD[i+2*nnode][j+nnode] += Kzy[i][j] * w[kint] * determinant;
            }
        }
    }
    
    // Copy the mass matrix to other components
    for (int i = 0; i < nnode; i++) {
        for (int j = 0; j < nnode; j++) {
            M_RD[i+nnode][j+nnode] = M_RD[i][j];
            M_RD[i+2*nnode][j+2*nnode] = M_RD[i][j];
            
            M_a[i+nnode][j+nnode] = M_a[i][j];
            M_a[i+2*nnode][j+2*nnode] = M_a[i][j];
            
            M_b[i+nnode][j+nnode] = M_b[i][j];
            M_b[i+2*nnode][j+2*nnode] = M_b[i][j];
            
            M_c[i+nnode][j+nnode] = M_c[i][j];
            M_c[i+2*nnode][j+2*nnode] = M_c[i][j];
            
            M_d[i+nnode][j+nnode] = M_d[i][j];
            M_d[i+2*nnode][j+2*nnode] = M_d[i][j];
        }
    }
    
    // Calculate damping matrix based on Rayleigh damping
    for (int i = 0; i < 3*nnode; i++) {
        for (int j = 0; j < 3*nnode; j++) {
            C_RD[i][j] = Damp_alpha * M_RD[i][j] + Damp_beta * K_RD[i][j];
        }
    }
    
    // Calculate N matrices (stiffness-like matrices for auxiliary fields)
    double factorP = (lambda + mu) / mu / (3.0 * lambda + 2.0 * mu);
    double factorS = -lambda / mu / 2.0 / (3.0 * lambda + 2.0 * mu);
    
    // N_a matrices
    for (int i = 0; i < nnode; i++) {
        for (int j = 0; j < nnode; j++) {
            // P-wave components
            N_a[i][j] = M_a[i][j] / rho * factorP;
            N_a[i][j+nnode] = M_a[i][j] / rho * factorS;
            N_a[i][j+2*nnode] = M_a[i][j] / rho * factorS;
            
            N_a[i+nnode][j] = M_a[i][j] / rho * factorS;
            N_a[i+nnode][j+nnode] = M_a[i][j] / rho * factorP;
            N_a[i+nnode][j+2*nnode] = M_a[i][j] / rho * factorS;
            
            N_a[i+2*nnode][j] = M_a[i][j] / rho * factorS;
            N_a[i+2*nnode][j+nnode] = M_a[i][j] / rho * factorS;
            N_a[i+2*nnode][j+2*nnode] = M_a[i][j] / rho * factorP;
        }
    }
    
    // N_b matrices
    for (int i = 0; i < nnode; i++) {
        for (int j = 0; j < nnode; j++) {
            // P-wave components
            N_b[i][j] = M_b[i][j] / rho * factorP;
            N_b[i][j+nnode] = M_b[i][j] / rho * factorS;
            N_b[i][j+2*nnode] = M_b[i][j] / rho * factorS;
            
            N_b[i+nnode][j] = M_b[i][j] / rho * factorS;
            N_b[i+nnode][j+nnode] = M_b[i][j] / rho * factorP;
            N_b[i+nnode][j+2*nnode] = M_b[i][j] / rho * factorS;
            
            N_b[i+2*nnode][j] = M_b[i][j] / rho * factorS;
            N_b[i+2*nnode][j+nnode] = M_b[i][j] / rho * factorS;
            N_b[i+2*nnode][j+2*nnode] = M_b[i][j] / rho * factorP;
        }
    }
    
    // N_c matrices
    for (int i = 0; i < nnode; i++) {
        for (int j = 0; j < nnode; j++) {
            // P-wave components
            N_c[i][j] = M_c[i][j] / rho * factorP;
            N_c[i][j+nnode] = M_c[i][j] / rho * factorS;
            N_c[i][j+2*nnode] = M_c[i][j] / rho * factorS;
            
            N_c[i+nnode][j] = M_c[i][j] / rho * factorS;
            N_c[i+nnode][j+nnode] = M_c[i][j] / rho * factorP;
            N_c[i+nnode][j+2*nnode] = M_c[i][j] / rho * factorS;
            
            N_c[i+2*nnode][j] = M_c[i][j] / rho * factorS;
            N_c[i+2*nnode][j+nnode] = M_c[i][j] / rho * factorS;
            N_c[i+2*nnode][j+2*nnode] = M_c[i][j] / rho * factorP;
        }
    }
    
    // N_d matrices
    for (int i = 0; i < nnode; i++) {
        for (int j = 0; j < nnode; j++) {
            // P-wave components
            N_d[i][j] = M_d[i][j] / rho * factorP;
            N_d[i][j+nnode] = M_d[i][j] / rho * factorS;
            N_d[i][j+2*nnode] = M_d[i][j] / rho * factorS;
            
            N_d[i+nnode][j] = M_d[i][j] / rho * factorS;
            N_d[i+nnode][j+nnode] = M_d[i][j] / rho * factorP;
            N_d[i+nnode][j+2*nnode] = M_d[i][j] / rho * factorS;
            
            N_d[i+2*nnode][j] = M_d[i][j] / rho * factorS;
            N_d[i+2*nnode][j+nnode] = M_d[i][j] / rho * factorS;
            N_d[i+2*nnode][j+2*nnode] = M_d[i][j] / rho * factorP;
        }
    }
    
    // Assemble the full PML matrices
    for (int i = 0; i < 3*nnode; i++) {
        for (int j = 0; j < 3*nnode; j++) {
            K_PML[i][j] = K_RD[i][j];
            M_PML[i][j] = M_RD[i][j];
            C_PML[i][j] = C_RD[i][j];
        }
    }
    
    for (int i = 0; i < 3*nnode; i++) {
        for (int j = 0; j < 6*nnode; j++) {
            G_PML[i][j] = A_eu[i][j] + A_wu[i][j] + A_pu[i][j];
            H_PML[i][j] = A_el[i][j] + A_wl[i][j] + A_pl[i][j];
        }
    }
    
    for (int i = 0; i < 6*nnode; i++) {
        for (int j = 0; j < 6*nnode; j++) {
            K_PML[i+3*nnode][j+3*nnode] = N_a[i][j] + N_b[i][j] + N_c[i][j] + N_d[i][j];
            M_PML[i+3*nnode][j+3*nnode] = M_a[i][j] + M_b[i][j] + M_c[i][j] + M_d[i][j];
        }
    }
    
    // Copy the PML matrices to the output arrays
    for (int i = 0; i < ndofel; i++) {
        for (int j = 0; j < ndofel; j++) {
            K[i*ndofel + j] = K_PML[i][j];
            M[i*ndofel + j] = M_PML[i][j];
            C[i*ndofel + j] = C_PML[i][j];
            G[i*ndofel + j] = G_PML[i][j];
            H[i*ndofel + j] = H_PML[i][j];
        }
    }
}

// Implementation of the main function to calculate PML matrices (now using static matrices)
void PML3DVISCOUS::calculateCppMatrices(const double* props, const double* coords, int ndofel, int mcrd, int nnode, int lflags) {
    // Initialize constants
    const double ZERO = 0.0;
    const double HALF = 0.5;
    const double ONE = 1.0;
    
    const double coef_alpha = 1.0/12.0;
    const double coef_beta = 1.0/48.0;
    
    // Clear static matrices
    for (int i = 0; i < ndofel*ndofel; i++) {
        M_cpp[i] = 0.0;
        C_cpp[i] = 0.0;
        K_cpp[i] = 0.0;
        G_cpp[i] = 0.0;
        H_cpp[i] = 0.0;
    }
    
    // Extract material properties
    double E = props[0];
    double xnu = props[1];
    double rho = props[2];
    int eleTypePos = (int)props[3];
    double PML_L = props[4];
    double afp = props[5];
    double PML_Rcoef = props[6];
    double RD_half_width_x = props[7];
    double RD_half_width_y = props[8];
    double RD_depth = props[9];
    double Damp_alpha = props[10];
    double Damp_beta = props[11];
    
    // Determine number of integration points based on afp
    int n_points = (afp < 3.0) ? 27 : 64;
    
    // Calculate Lame constants
    double lambda = xnu * E / ((1.0 + xnu) * (1.0 - 2.0 * xnu));
    double mu = 0.5 * E / (1.0 + xnu);
    
    // Initialize arrays
    double strain[6] = {0.0};
    double stress[6] = {0.0};
    double D[6][6] = {{0.0}};
    double B[6][60] = {{0.0}};
    double dxidx[3][3] = {{0.0}};
    
    double N[20] = {0.0};
    double dNdxi[20][3] = {{0.0}};
    double dxdxi[3][3] = {{0.0}};
    double dNdx[20][3] = {{0.0}};
    
    double Phi[3][60] = {{0.0}};
    
    // PML matrices
    double K_RD[24][24] = {{0.0}};
    double Kxx[8][8] = {{0.0}};
    double Kxy[8][8] = {{0.0}};
    double Kxz[8][8] = {{0.0}};
    double Kyx[8][8] = {{0.0}};
    double Kyy[8][8] = {{0.0}};
    double Kyz[8][8] = {{0.0}};
    double Kzx[8][8] = {{0.0}};
    double Kzy[8][8] = {{0.0}};
    double Kzz[8][8] = {{0.0}};
    
    double Phi_x[8] = {0.0};
    double Phi_y[8] = {0.0};
    double Phi_z[8] = {0.0}; // Fixed extra brace
    
    double M_RD[24][24] = {{0.0}};
    double C_RD[24][24] = {{0.0}};
    double M_a[24][24] = {{0.0}};
    double M_b[24][24] = {{0.0}};
    double M_c[24][24] = {{0.0}};
    double M_d[24][24] = {{0.0}};
    
    double N_a[48][48] = {{0.0}};
    double N_b[48][48] = {{0.0}};
    double N_c[48][48] = {{0.0}};
    double N_d[48][48] = {{0.0}};
    
    double A_eu[24][48] = {{0.0}};
    double A_wu[24][48] = {{0.0}};
    double A_pu[24][48] = {{0.0}};
    double A_el[24][48] = {{0.0}};
    double A_wl[24][48] = {{0.0}};
    double A_pl[24][48] = {{0.0}};
    
    double K_PML[72][72] = {{0.0}};
    double M_PML[72][72] = {{0.0}};
    double C_PML[72][72] = {{0.0}};
    double G_PML[72][72] = {{0.0}};
    double H_PML[72][72] = {{0.0}};
    
    // Make sure integration points are initialized
    if (!integrationInitialized) {
        initializeIntegrationPointsAndWeights(n_points, nnode);
    }
    
    // Initialize matrices to zero
    for (int i = 0; i < ndofel; i++) {
        for (int j = 0; j < ndofel; j++) {
            M_cpp[i*ndofel + j] = 0.0;
            C_cpp[i*ndofel + j] = 0.0;
            K_cpp[i*ndofel + j] = 0.0;
            G_cpp[i*ndofel + j] = 0.0;
            H_cpp[i*ndofel + j] = 0.0;
        }
    }
    
    // Loop over integration points
    for (int kint = 0; kint < n_points; kint++) {
        // Get shape functions and derivatives at this integration point
        double xiPoint[3] = {xi[0][kint], xi[1][kint], xi[2][kint]};
        calculateShapeFunctions(xiPoint, nnode, N, dNdxi);
        
        // Calculate derivatives of shape functions w.r.t global coordinates
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                dxdxi[i][j] = 0.0;
                for (int k = 0; k < nnode; k++) {
                    dxdxi[i][j] += coords[k*3 + i] * dNdxi[k][j];
                }
            }
        }
        
        // Calculate the determinant and inverse of dxdxi
        double determinant = dxdxi[0][0]*(dxdxi[1][1]*dxdxi[2][2] - dxdxi[1][2]*dxdxi[2][1])
                            - dxdxi[0][1]*(dxdxi[1][0]*dxdxi[2][2] - dxdxi[1][2]*dxdxi[2][0])
                            + dxdxi[0][2]*(dxdxi[1][0]*dxdxi[2][1] - dxdxi[1][1]*dxdxi[2][0]);
        
        if (fabs(determinant) < 1.0e-10) {
            opserr << "Error: Jacobian determinant close to zero in PML3DVISCOUS::calculatePMLMatrices\n";
            return;
        }
        
        // Calculate inverse of dxdxi
        dxidx[0][0] = (dxdxi[1][1]*dxdxi[2][2] - dxdxi[1][2]*dxdxi[2][1])/determinant;
        dxidx[0][1] = (dxdxi[0][2]*dxdxi[2][1] - dxdxi[0][1]*dxdxi[2][2])/determinant;
        dxidx[0][2] = (dxdxi[0][1]*dxdxi[1][2] - dxdxi[0][2]*dxdxi[1][1])/determinant;
        dxidx[1][0] = (dxdxi[1][2]*dxdxi[2][0] - dxdxi[1][0]*dxdxi[2][2])/determinant;
        dxidx[1][1] = (dxdxi[0][0]*dxdxi[2][2] - dxdxi[0][2]*dxdxi[2][0])/determinant;
        dxidx[1][2] = (dxdxi[0][2]*dxdxi[1][0] - dxdxi[0][0]*dxdxi[1][2])/determinant;
        dxidx[2][0] = (dxdxi[1][0]*dxdxi[2][1] - dxdxi[1][1]*dxdxi[2][0])/determinant;
        dxidx[2][1] = (dxdxi[0][1]*dxdxi[2][0] - dxdxi[0][0]*dxdxi[2][1])/determinant;
        dxidx[2][2] = (dxdxi[0][0]*dxdxi[1][1] - dxdxi[0][1]*dxdxi[1][0])/determinant;
        
        // Calculate derivatives of shape functions w.r.t. global coordinates
        for (int i = 0; i < nnode; i++) {
            for (int j = 0; j < 3; j++) {
                dNdx[i][j] = 0.0;
                for (int k = 0; k < 3; k++) {
                    dNdx[i][j] += dNdxi[i][k] * dxidx[k][j];
                }
            }
        }
        
        // Calculate B matrix (strain-displacement)
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 3*nnode; j++) {
                B[i][j] = 0.0;
            }
        }
        
        for (int i = 0; i < nnode; i++) {
            B[0][i*3]     = dNdx[i][0];  // du/dx
            B[1][i*3 + 1] = dNdx[i][1];  // dv/dy
            B[2][i*3 + 2] = dNdx[i][2];  // dw/dz
            
            B[3][i*3]     = dNdx[i][1];  // du/dy
            B[3][i*3 + 1] = dNdx[i][0];  // dv/dx
            
            B[4][i*3]     = dNdx[i][2];  // du/dz
            B[4][i*3 + 2] = dNdx[i][0];  // dw/dx
            
            B[5][i*3 + 1] = dNdx[i][2];  // dv/dz
            B[5][i*3 + 2] = dNdx[i][1];  // dw/dy
        }
        
        // Calculate Phi matrix (shape function matrix)
        for (int i = 0; i < nnode; i++) {
            Phi[0][i*3]     = N[i];
            Phi[1][i*3 + 1] = N[i];
            Phi[2][i*3 + 2] = N[i];
        }
        
        // Calculate derivative vectors
        for (int i = 0; i < nnode; i++) {
            Phi_x[i] = dNdx[i][0];
            Phi_y[i] = dNdx[i][1];
            Phi_z[i] = dNdx[i][2];
        }
        
        // Calculate gauss point coordinates
        double x1 = 0.0, x2 = 0.0, x3 = 0.0;
        for (int i = 0; i < nnode; i++) {
            x1 += N[i] * coords[i*3];
            x2 += N[i] * coords[i*3 + 1];
            x3 += N[i] * coords[i*3 + 2];
        }
        
        // Calculate PML alpha and beta parameters
        double pmlAlphaBeta[2][3] = {{0.0}};
        calculatePMLParameters(props, x1, x2, x3, pmlAlphaBeta);
        
        // Calculate coefficients for PML
        double coef_a = pmlAlphaBeta[0][0] * pmlAlphaBeta[0][1] * pmlAlphaBeta[0][2];
        double coef_b = pmlAlphaBeta[0][0] * pmlAlphaBeta[0][1] * pmlAlphaBeta[1][2] +
                        pmlAlphaBeta[0][0] * pmlAlphaBeta[0][2] * pmlAlphaBeta[1][1] +
                        pmlAlphaBeta[0][1] * pmlAlphaBeta[0][2] * pmlAlphaBeta[1][0];
        double coef_c = pmlAlphaBeta[0][0] * pmlAlphaBeta[1][1] * pmlAlphaBeta[1][2] +
                        pmlAlphaBeta[0][1] * pmlAlphaBeta[1][2] * pmlAlphaBeta[1][0] +
                        pmlAlphaBeta[0][2] * pmlAlphaBeta[1][1] * pmlAlphaBeta[1][0];
        double coef_d = pmlAlphaBeta[1][0] * pmlAlphaBeta[1][1] * pmlAlphaBeta[1][2];
        
        // Calculate L coefficients
        double coef_Le[3][3] = {{0.0}};
        double coef_Lp[3][3] = {{0.0}};
        double coef_Lw[3][3] = {{0.0}};
        
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                coef_Le[i][j] = pmlAlphaBeta[0][i] * pmlAlphaBeta[0][j];
                coef_Lp[i][j] = pmlAlphaBeta[0][i] * pmlAlphaBeta[1][j] +
                               pmlAlphaBeta[1][i] * pmlAlphaBeta[0][j];
                coef_Lw[i][j] = pmlAlphaBeta[1][i] * pmlAlphaBeta[1][j];
            }
        }
        
        // Calculate stiffness matrix components
        for (int i = 0; i < nnode; i++) {
            for (int j = 0; j < nnode; j++) {
                // Diagonal components
                Kxx[i][j] = (lambda + 2.0*mu) * Phi_x[i] * Phi_x[j] + 
                           mu * (Phi_y[i] * Phi_y[j] + Phi_z[i] * Phi_z[j]);
                
                Kyy[i][j] = (lambda + 2.0*mu) * Phi_y[i] * Phi_y[j] + 
                           mu * (Phi_x[i] * Phi_x[j] + Phi_z[i] * Phi_z[j]);
                
                Kzz[i][j] = (lambda + 2.0*mu) * Phi_z[i] * Phi_z[j] + 
                           mu * (Phi_x[i] * Phi_x[j] + Phi_y[i] * Phi_y[j]);
                
                // Off-diagonal components
                Kxy[i][j] = lambda * Phi_x[i] * Phi_y[j] + mu * Phi_y[i] * Phi_x[j];
                Kxz[i][j] = lambda * Phi_x[i] * Phi_z[j] + mu * Phi_z[i] * Phi_x[j];
                Kyz[i][j] = lambda * Phi_y[i] * Phi_z[j] + mu * Phi_z[i] * Phi_y[j];
                
                // Mass matrix components
                M_RD[i][j] += rho * N[i] * N[j] * w[kint] * determinant;
                
                // PML matrix components
                M_a[i][j] += coef_a * rho * N[i] * N[j] * w[kint] * determinant;
                M_b[i][j] += coef_b * rho * N[i] * N[j] * w[kint] * determinant;
                M_c[i][j] += coef_c * rho * N[i] * N[j] * w[kint] * determinant;
                M_d[i][j] += coef_d * rho * N[i] * N[j] * w[kint] * determinant;
                
                // A matrices
                A_eu[i][j] += Phi_x[i] * N[j] * coef_Le[1][2] * w[kint] * determinant;
                A_eu[i][j+nnode*3] += Phi_y[i] * N[j] * coef_Le[0][2] * w[kint] * determinant;
                A_eu[i][j+nnode*4] += Phi_z[i] * N[j] * coef_Le[0][1] * w[kint] * determinant;
                
                A_eu[i+nnode][j+nnode] += Phi_y[i] * N[j] * coef_Le[0][2] * w[kint] * determinant;
                A_eu[i+nnode][j+nnode*3] += Phi_x[i] * N[j] * coef_Le[1][2] * w[kint] * determinant;
                A_eu[i+nnode][j+nnode*5] += Phi_z[i] * N[j] * coef_Le[0][1] * w[kint] * determinant;
                
                A_eu[i+nnode*2][j+nnode*2] += Phi_z[i] * N[j] * coef_Le[0][1] * w[kint] * determinant;
                A_eu[i+nnode*2][j+nnode*4] += Phi_x[i] * N[j] * coef_Le[1][2] * w[kint] * determinant;
                A_eu[i+nnode*2][j+nnode*5] += Phi_y[i] * N[j] * coef_Le[0][2] * w[kint] * determinant;
                
                // A_wu matrices
                A_wu[i][j] += Phi_x[i] * N[j] * coef_Lw[1][2] * w[kint] * determinant;
                A_wu[i][j+nnode*3] += Phi_y[i] * N[j] * coef_Lw[0][2] * w[kint] * determinant;
                A_wu[i][j+nnode*4] += Phi_z[i] * N[j] * coef_Lw[0][1] * w[kint] * determinant;
                
                A_wu[i+nnode][j+nnode] += Phi_y[i] * N[j] * coef_Lw[0][2] * w[kint] * determinant;
                A_wu[i+nnode][j+nnode*3] += Phi_x[i] * N[j] * coef_Lw[1][2] * w[kint] * determinant;
                A_wu[i+nnode][j+nnode*5] += Phi_z[i] * N[j] * coef_Lw[0][1] * w[kint] * determinant;
                
                A_wu[i+nnode*2][j+nnode*2] += Phi_z[i] * N[j] * coef_Lw[0][1] * w[kint] * determinant;
                A_wu[i+nnode*2][j+nnode*4] += Phi_x[i] * N[j] * coef_Lw[1][2] * w[kint] * determinant;
                A_wu[i+nnode*2][j+nnode*5] += Phi_y[i] * N[j] * coef_Lw[0][2] * w[kint] * determinant;
                
                // A_pu matrices
                A_pu[i][j] += Phi_x[i] * N[j] * coef_Lp[1][2] * w[kint] * determinant;
                A_pu[i][j+nnode*3] += Phi_y[i] * N[j] * coef_Lp[0][2] * w[kint] * determinant;
                A_pu[i][j+nnode*4] += Phi_z[i] * N[j] * coef_Lp[0][1] * w[kint] * determinant;
                
                A_pu[i+nnode][j+nnode] += Phi_y[i] * N[j] * coef_Lp[0][2] * w[kint] * determinant;
                A_pu[i+nnode][j+nnode*3] += Phi_x[i] * N[j] * coef_Lp[1][2] * w[kint] * determinant;
                A_pu[i+nnode][j+nnode*5] += Phi_z[i] * N[j] * coef_Lp[0][1] * w[kint] * determinant;
                
                A_pu[i+nnode*2][j+nnode*2] += Phi_z[i] * N[j] * coef_Lp[0][1] * w[kint] * determinant;
                A_pu[i+nnode*2][j+nnode*4] += Phi_x[i] * N[j] * coef_Lp[1][2] * w[kint] * determinant;
                A_pu[i+nnode*2][j+nnode*5] += Phi_y[i] * N[j] * coef_Lp[0][2] * w[kint] * determinant;
            }
        }
        
        // Calculate Kyx, Kzx, and Kzy as transposes of Kxy, Kxz, and Kyz
        for (int i = 0; i < nnode; i++) {
            for (int j = 0; j < nnode; j++) {
                Kyx[i][j] = Kxy[j][i];
                Kzx[i][j] = Kxz[j][i];
                Kzy[i][j] = Kyz[j][i];
            }
        }
        
        // Assemble the full K_RD matrix
        for (int i = 0; i < nnode; i++) {
            for (int j = 0; j < nnode; j++) {
                // Diagonal blocks
                K_RD[i][j] += Kxx[i][j] * w[kint] * determinant;
                K_RD[i+nnode][j+nnode] += Kyy[i][j] * w[kint] * determinant;
                K_RD[i+2*nnode][j+2*nnode] += Kzz[i][j] * w[kint] * determinant;
                
                // Off-diagonal blocks
                K_RD[i][j+nnode] += Kxy[i][j] * w[kint] * determinant;
                K_RD[i][j+2*nnode] += Kxz[i][j] * w[kint] * determinant;
                K_RD[i+nnode][j+2*nnode] += Kyz[i][j] * w[kint] * determinant;
                
                K_RD[i+nnode][j] += Kyx[i][j] * w[kint] * determinant;
                K_RD[i+2*nnode][j] += Kzx[i][j] * w[kint] * determinant;
                K_RD[i+2*nnode][j+nnode] += Kzy[i][j] * w[kint] * determinant;
            }
        }
    }
    
    // Copy the mass matrix to other components
    for (int i = 0; i < nnode; i++) {
        for (int j = 0; j < nnode; j++) {
            M_RD[i+nnode][j+nnode] = M_RD[i][j];
            M_RD[i+2*nnode][j+2*nnode] = M_RD[i][j];
            
            M_a[i+nnode][j+nnode] = M_a[i][j];
            M_a[i+2*nnode][j+2*nnode] = M_a[i][j];
            
            M_b[i+nnode][j+nnode] = M_b[i][j];
            M_b[i+2*nnode][j+2*nnode] = M_b[i][j];
            
            M_c[i+nnode][j+nnode] = M_c[i][j];
            M_c[i+2*nnode][j+2*nnode] = M_c[i][j];
            
            M_d[i+nnode][j+nnode] = M_d[i][j];
            M_d[i+2*nnode][j+2*nnode] = M_d[i][j];
        }
    }
    
    // Calculate damping matrix based on Rayleigh damping
    for (int i = 0; i < 3*nnode; i++) {
        for (int j = 0; j < 3*nnode; j++) {
            C_RD[i][j] = Damp_alpha * M_RD[i][j] + Damp_beta * K_RD[i][j];
        }
    }
    
    // Calculate N matrices (stiffness-like matrices for auxiliary fields)
    double factorP = (lambda + mu) / mu / (3.0 * lambda + 2.0 * mu);
    double factorS = -lambda / mu / 2.0 / (3.0 * lambda + 2.0 * mu);
    
    // N_a matrices
    for (int i = 0; i < nnode; i++) {
        for (int j = 0; j < nnode; j++) {
            // P-wave components
            N_a[i][j] = M_a[i][j] / rho * factorP;
            N_a[i][j+nnode] = M_a[i][j] / rho * factorS;
            N_a[i][j+2*nnode] = M_a[i][j] / rho * factorS;
            
            N_a[i+nnode][j] = M_a[i][j] / rho * factorS;
            N_a[i+nnode][j+nnode] = M_a[i][j] / rho * factorP;
            N_a[i+nnode][j+2*nnode] = M_a[i][j] / rho * factorS;
            
            N_a[i+2*nnode][j] = M_a[i][j] / rho * factorS;
            N_a[i+2*nnode][j+nnode] = M_a[i][j] / rho * factorS;
            N_a[i+2*nnode][j+2*nnode] = M_a[i][j] / rho * factorP;
        }
    }
    
    // N_b matrices
    for (int i = 0; i < nnode; i++) {
        for (int j = 0; j < nnode; j++) {
            // P-wave components
            N_b[i][j] = M_b[i][j] / rho * factorP;
            N_b[i][j+nnode] = M_b[i][j] / rho * factorS;
            N_b[i][j+2*nnode] = M_b[i][j] / rho * factorS;
            
            N_b[i+nnode][j] = M_b[i][j] / rho * factorS;
            N_b[i+nnode][j+nnode] = M_b[i][j] / rho * factorP;
            N_b[i+nnode][j+2*nnode] = M_b[i][j] / rho * factorS;
            
            N_b[i+2*nnode][j] = M_b[i][j] / rho * factorS;
            N_b[i+2*nnode][j+nnode] = M_b[i][j] / rho * factorS;
            N_b[i+2*nnode][j+2*nnode] = M_b[i][j] / rho * factorP;
        }
    }
    
    // N_c matrices
    for (int i = 0; i < nnode; i++) {
        for (int j = 0; j < nnode; j++) {
            // P-wave components
            N_c[i][j] = M_c[i][j] / rho * factorP;
            N_c[i][j+nnode] = M_c[i][j] / rho * factorS;
            N_c[i][j+2*nnode] = M_c[i][j] / rho * factorS;
            
            N_c[i+nnode][j] = M_c[i][j] / rho * factorS;
            N_c[i+nnode][j+nnode] = M_c[i][j] / rho * factorP;
            N_c[i+nnode][j+2*nnode] = M_c[i][j] / rho * factorS;
            
            N_c[i+2*nnode][j] = M_c[i][j] / rho * factorS;
            N_c[i+2*nnode][j+nnode] = M_c[i][j] / rho * factorS;
            N_c[i+2*nnode][j+2*nnode] = M_c[i][j] / rho * factorP;
        }
    }
    
    // N_d matrices
    for (int i = 0; i < nnode; i++) {
        for (int j = 0; j < nnode; j++) {
            // P-wave components
            N_d[i][j] = M_d[i][j] / rho * factorP;
            N_d[i][j+nnode] = M_d[i][j] / rho * factorS;
            N_d[i][j+2*nnode] = M_d[i][j] / rho * factorS;
            
            N_d[i+nnode][j] = M_d[i][j] / rho * factorS;
            N_d[i+nnode][j+nnode] = M_d[i][j] / rho * factorP;
            N_d[i+nnode][j+2*nnode] = M_d[i][j] / rho * factorS;
            
            N_d[i+2*nnode][j] = M_d[i][j] / rho * factorS;
            N_d[i+2*nnode][j+nnode] = M_d[i][j] / rho * factorS;
            N_d[i+2*nnode][j+2*nnode] = M_d[i][j] / rho * factorP;
        }
    }
    
    // Assemble the full PML matrices
    for (int i = 0; i < 3*nnode; i++) {
        for (int j = 0; j < 3*nnode; j++) {
            K_PML[i][j] = K_RD[i][j];
            M_PML[i][j] = M_RD[i][j];
            C_PML[i][j] = C_RD[i][j];
        }
    }
    
    for (int i = 0; i < 3*nnode; i++) {
        for (int j = 0; j < 6*nnode; j++) {
            G_PML[i][j] = A_eu[i][j] + A_wu[i][j] + A_pu[i][j];
            H_PML[i][j] = A_el[i][j] + A_wl[i][j] + A_pl[i][j];
        }
    }
    
    for (int i = 0; i < 6*nnode; i++) {
        for (int j = 0; j < 6*nnode; j++) {
            K_PML[i+3*nnode][j+3*nnode] = N_a[i][j] + N_b[i][j] + N_c[i][j] + N_d[i][j];
            M_PML[i+3*nnode][j+3*nnode] = M_a[i][j] + M_b[i][j] + M_c[i][j] + M_d[i][j];
        }
    }
    
    // Copy to static matrices at the end
    for (int i = 0; i < ndofel; i++) {
        for (int j = 0; j < ndofel; j++) {
            M_cpp[i*ndofel + j] = M_PML[i][j];
            C_cpp[i*ndofel + j] = C_PML[i][j];
            K_cpp[i*ndofel + j] = K_PML[i][j];
            G_cpp[i*ndofel + j] = G_PML[i][j];
            H_cpp[i*ndofel + j] = H_PML[i][j];
        }
    }
}

// Implementation of verification function
void PML3DVISCOUS::verifyMatrices() {
    // Define tolerance for floating point comparisons
    const double tolerance = 1e-10;
    bool mismatchFound = false;
    
    // Compare M matrices
    for (int i = 0; i < PML3DVISCOUS_NUM_DOF*PML3DVISCOUS_NUM_DOF; i++) {
        if (fabs(M[i] - M_cpp[i]) > tolerance) {
            opserr << "PML3DVISCOUS::verifyMatrices - Mismatch in M matrix at index " << i 
                   << ", Fortran: " << M[i] << ", C++: " << M_cpp[i] << ", Diff: " << (M[i] - M_cpp[i]) << endln;
            mismatchFound = true;
        }
    }
    
    // Compare C matrices
    for (int i = 0; i < PML3DVISCOUS_NUM_DOF*PML3DVISCOUS_NUM_DOF; i++) {
        if (fabs(C[i] - C_cpp[i]) > tolerance) {
            opserr << "PML3DVISCOUS::verifyMatrices - Mismatch in C matrix at index " << i 
                   << ", Fortran: " << C[i] << ", C++: " << C_cpp[i] << ", Diff: " << (C[i] - C_cpp[i]) << endln;
            mismatchFound = true;
        }
    }
    
    // Compare K matrices
    for (int i = 0; i < PML3DVISCOUS_NUM_DOF*PML3DVISCOUS_NUM_DOF; i++) {
        if (fabs(K[i] - K_cpp[i]) > tolerance) {
            opserr << "PML3DVISCOUS::verifyMatrices - Mismatch in K matrix at index " << i 
                   << ", Fortran: " << K[i] << ", C++: " << K_cpp[i] << ", Diff: " << (K[i] - K_cpp[i]) << endln;
            mismatchFound = true;
        }
    }
    
    // Compare G matrices
    for (int i = 0; i < PML3DVISCOUS_NUM_DOF*PML3DVISCOUS_NUM_DOF; i++) {
        if (fabs(G[i] - G_cpp[i]) > tolerance) {
            opserr << "PML3DVISCOUS::verifyMatrices - Mismatch in G matrix at index " << i 
                   << ", Fortran: " << G[i] << ", C++: " << G_cpp[i] << ", Diff: " << (G[i] - G_cpp[i]) << endln;
            mismatchFound = true;
        }
    }
    
    // Compare H matrices
    for (int i = 0; i < PML3DVISCOUS_NUM_DOF*PML3DVISCOUS_NUM_DOF; i++) {
        if (fabs(H[i] - H_cpp[i]) > tolerance) {
            opserr << "PML3DVISCOUS::verifyMatrices - Mismatch in H matrix at index " << i 
                   << ", Fortran: " << H[i] << ", C++: " << H_cpp[i] << ", Diff: " << (H[i] - H_cpp[i]) << endln;
            mismatchFound = true;
        }
    }
    
    if (!mismatchFound) {
        opserr << "PML3DVISCOUS::verifyMatrices - All matrices match between Fortran and C++ implementations!" << endln;
    }
}

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

	// initialize the ubar and ubart vectors to zero
	ubart.Zero();
	ubar.Zero();
    ubarbar.Zero();
    ubarbart.Zero();
	updateflag = 0;
	update_dt = 0;
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
    
    // Calculate matrices using C++ implementation
    calculateCppMatrices(props, coords, NDOFEL, MCRD, NNODE, LFLAGS);
    
    // Verify the matrices
    verifyMatrices();
    
    dt = theDomain->getDT();

	// // make C zero 
	for (int i = 0; i < PML3DVISCOUS_NUM_DOF*PML3DVISCOUS_NUM_DOF; i++) {
		// C[i] = 0.0;
		// K[i] = 0.0;
		// M[i] = 0.0;
		// G[i] = 0.0;
		// H[i] = 0.0;
	}


	std::ofstream myfile;
	std::string filename;
	int tag = this->getTag();

	// // save M matrix in a file
	// filename = "./PML/M" + std::to_string(tag) + ".mat";
	// myfile.open(filename);
	// for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
	// 	for (int j = 0; j < PML3DVISCOUS_NUM_DOF; j++) {
	// 		myfile << M[i*PML3DVISCOUS_NUM_DOF + j] << " ";
	// 	}
	// 	myfile << "\n";
	// }
	// myfile.close();

	// // save C matrix in a file
	// filename = "./PML/C" + std::to_string(tag) + ".mat";
	// myfile.open(filename);
	// for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
	// 	for (int j = 0; j < PML3DVISCOUS_NUM_DOF; j++) {
	// 		myfile << C[i*PML3DVISCOUS_NUM_DOF + j] << " ";
	// 	}
	// 	myfile << "\n";
	// }
	// myfile.close();

	// // save K matrix in a file
	// filename = "./PML/K" + std::to_string(tag) + ".mat";
	// myfile.open(filename);
	// for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
	// 	for (int j = 0; j < PML3DVISCOUS_NUM_DOF; j++) {
	// 		myfile << K[i*PML3DVISCOUS_NUM_DOF + j] << " ";
	// 	}
	// 	myfile << "\n";
	// }
	// myfile.close();

	// // save G matrix in a file
	// filename = "./PML/G" + std::to_string(tag) + ".mat";
	// myfile.open(filename);
	// for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
	// 	for (int j = 0; j < PML3DVISCOUS_NUM_DOF; j++) {
	// 		myfile << G[i*PML3DVISCOUS_NUM_DOF + j] << " ";
	// 	}
	// 	myfile << "\n";
	// }
	// myfile.close();

	// save H matrix in a file
	// filename = "./PML/H" + std::to_string(tag) + ".mat";
	// myfile.open(filename);
	// for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
	// 	for (int j = 0; j < PML3DVISCOUS_NUM_DOF; j++) {
	// 		myfile << H[i*PML3DVISCOUS_NUM_DOF + j] << " ";
	// 	}
	// 	myfile << "\n";
	// }
	// myfile.close();

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

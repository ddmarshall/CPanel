/*******************************************************************************
 * Copyright (c) 2014 Chris Satterwhite
 * Copyright (c) 2018 David D. Marshall <ddmarsha@calpoly.edu>
 *
 * This program and the accompanying materials are made
 * available under the terms of the Eclipse Public License 2.0
 * which is available at https://www.eclipse.org/legal/epl-2.0/
 *
 * See LICENSE.md file in the project root for full license information.
 *
 * SPDX-License-Identifier: EPL-2.0
 *
 * Contributors:
 *    Chris Satterwhite - initial code and implementation
 *    David D. Marshall - misc. changes
 ******************************************************************************/

#include <stdio.h>
#include <iostream>
#include "CPanelMgr.h"
#include "inputParams.h"
#include "geometry.h"
#include "cpFile.h"
#include "cpNode.h"


void usage(const char * argv[])
{
    printf("\n");
    printf("CPanel is an unstructured panel code developed by students at California Polytechnic State University - San Luis Obispo.\n");
    printf("\n");
    printf("Usage: %s <infile> \n",argv[0]);
    printf("   <infile>   Required parameter to specify settings input file.\n");
    printf("\n");
    printf("The input file contains the geometry filename, reference values,freestream conditions, and solver options.\n\n");
    printf("Currently supported geometry file formats are:\n");
    printf("\t.tri (Cart3D format)\n");
    printf("\t.tricp (Modified .tri file including normal vectors from underlying bezier surfaces)\n\n");
    printf("This example below may be used as a template for the input file.\n");
    printf("\n");
    printf("%%%% CPanel Input File %%%%\n\n");
    printf("%% Reference Geometry %%\n");
    printf("GeomFile = wing.tri\n");
    printf("S_ref =  6.0\n");
    printf("b_ref =  6.0\n");
    printf("c_ref =  1.0\n");
    printf("X_cg =  0.25\n");
    printf("Y_cg =  0.0\n");
    printf("Z_cg =  0.0\n\n");
    printf("%% Cases %%\n");
    printf("Velocity (ft/s)\n");
    printf("3\n");
    printf("120\n");
    printf("160\n");
    printf("200\n");
    printf("Angle_of_Attack (degrees)\n");
    printf("6\n");
    printf("0\n");
    printf("2\n");
    printf("4\n");
    printf("6\n");
    printf("8\n");
    printf("10\n");
    printf("Angle_of_Sideslip (degrees)\n");
    printf("1\n");
    printf("0\n");
    printf("Mach_Number\n");
    printf("1\n");
    printf("0.3\n");
    printf("\n");
    printf("%% Solver Options (0 = OFF, 1 = ON) %%\n");
    printf("Surface_Streamlines\n");
    printf("1\n");
    printf("Stability_Derivatives\n");
    printf("1\n");
    printf("Write_Influence_Coefficients\n");
    printf("0\n");
    printf("Vortex_Particle_Wake\n");
    printf("0\n");
    printf("Accelerate_Code\n");
    printf("1\n");
    printf("Volume_Mesh (Xo Xf Yo Yf Zo Zf nX nY nZ)\n");
    printf("1\n");
    printf("-0.5 6.5 -8 8 -0.5 1.25 40 1 30\n");
    printf("\n");
    printf("%% Vortex Particle Wake Options %%\n");
    printf("Time_Step\n");
    printf("0\n");
    printf("Number_of_Timesteps\n");
    printf("40\n");
    printf("Unsteady_Mode\n");
    printf("0\n");
    printf("/filepath/kinematicFile.bKin\n\n");


    printf("The number following Case variables indicates the number of different values for those variables, followed by the actual values\n\n");
    printf("Write_Influence_Coefficients option writes dense matrix to file for future use. Saves time of recomputing influence coefficients if running the geometry multiple times. For fine meshes (> 10000 panels), can take significant amount of time for I/O and file can become very large.\n\n");
    printf("\n");
    printf("For questions, contact :\n");
    printf("\tDavid D. Marshall\n");
    printf("\tCal Poly, San Luis Obispo\n");
    printf("\tddmarsha@calpoly.edu\n\n");
    printf("v. 1.0 - 1/24/2014\n");
    printf("\n");
    exit(1);
}

int main(int argc, const char * argv[])
{
    // Start Timer
    time_t ts,tf;
    time(&ts);

    // Check arguments
    if (argc != 2)
    {
        printf("ERROR : Incorrect Usage");
        usage(argv);
        exit(EXIT_FAILURE);
    }

    // Check for file existence
    cpFile inFile(argv[1]);

    inputParams inData(&inFile);

    if (!inData.set())
    {
        usage(argv);
        exit(EXIT_FAILURE);
    }

    std::cout << "Running CPanel with the following inputs...\n" << std::endl;
    inData.print(std::cout);
    std::cout << std::endl;
    geometry geom(&inData);
    caseMgr cm(&inData,&geom);


    time(&tf);
    std::cout << "Elapsed time for program execution : " << difftime(tf,ts) << " seconds" << std::endl;
}

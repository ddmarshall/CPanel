/*******************************************************************************
 * Copyright (c) 2017 Connor Sousa
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
 *    Connor Sousa - initial code and implementation
 *    David D. Marshall - misc. changes
 ******************************************************************************/

#ifndef __CPanel__cpCaseVP__
#define __CPanel__cpCaseVP__

#include <iostream>
#include <string>
#include "cpCase.h"
#include "particle.h"
#include "vortexFil.h"
#include "particleFMM.h"

#include "octreeFile.h"



class cpCaseVP : public cpCase
{

	using particles_type = std::vector<particle *>;
	using particles_index_type = particles_type::size_type;
	using filaments_type = std::vector<vortexFil *>;
	using filaments_index_type = filaments_type::size_type;

    particles_type particles;
    filaments_type filaments;
    wakePanels_type * w2panels; // Second row of buffer wake
    
    Eigen::VectorXd wake2Doublets;
    
    void convectBufferWake();
    void setSourceStrengthsVP();
    bool solutionConvergence();
    
    bool solveMatrixEq();
    bool solveVPmatrixEq();
    bool edgeIsUsed(edge* thisEdge, std::vector<edge*> pEdges);
    void collapseBufferWake();
    void compVelocity();
    void trefftzPlaneAnalysisVP();
    void stabilityDerivativesVP();
    void particleStrengthUpdate();
    void convectParticles();
    
    Eigen::Vector3d Vinf(Eigen::Vector3d POI);
    Eigen::Vector3d VinfPlusVecPot(Eigen::Vector3d POI);
    
    Eigen::Vector3d rungeKuttaStepper(Eigen::Vector3d POI);
    Eigen::Vector3d velocityInflFromEverything(Eigen::Vector3d POI);
    Eigen::Vector3d velocityInflFromEverything(particle* part);

    void writeFilesVP();
    void writeBodyDataVP(const std::string & path,const Eigen::MatrixXd &nodeMat);
    void writeWakeDataVP(const std::string & path, const Eigen::MatrixXd &nodeMat);
    void writeBuffWake2Data(const std::string & path, const Eigen::MatrixXd &nodeMat);
    void writeFilamentData(const std::string & path);
    void writeParticleData(const std::string & path);

    void readBodyKinFile();
    void populateVolMesh();
    
    int timestep = 1;
    double dt;
    int numSteps;
    bool accelerate;
    bool unsteady;
    
    
    particleOctree partOctree;
    particleFMM FMM;
    Eigen::MatrixXd bodyKin;
        
    //    void moveGeometry();
    
    
public:
    cpCaseVP( geometry *ggeom, double V, double aalpha, double bbeta, double mmach, inputParams* inParams )
      : cpCase( ggeom, V, aalpha, bbeta, mmach, inParams )
    {
        bPanels = geom->getBodyPanels();
        wPanels = geom->getWakePanels();
        w2panels = geom->getBufferWake2Panels();
        dt = geom->getDt();
        PG = sqrt(1-pow(mach,2));
        
        numSteps = inParams->numSteps;
        accelerate = inParams->accelerateCode;
        unsteady = inParams->unsteady;
    }
    
    
    void run(bool printFlag, bool surfStreamFlag, bool stabDerivFlag);
    
    
};





#endif /* defined(__CPanel__cpCaseVP__) */

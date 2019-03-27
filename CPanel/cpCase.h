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

#ifndef __CPanel__runCase__
#define __CPanel__runCase__

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <Eigen/Dense>
#include <vector>
#include <Eigen/IterativeLinearSolvers>
#include "geometry.h"
#include "VTUfile.h"
#include "bodyStreamline.h"
#include "cpFile.h"
#include "inputParams.h"

#include "octreeFile.h"



class cpCase
{
protected:
    geometry *geom;
    
    inputParams* params;
    double Vmag;
    double mach;
    double PG; // Prandtl-Glauert Correction - (1-M^2)^(1/2)
    double alpha;
    double beta;
    
    Eigen::Vector3d Vinf;
    Eigen::Matrix3d transform;
    
    using bodyPanels_type = std::vector<bodyPanel *>;
    using bodyPanels_index_type = bodyPanels_type::size_type;
    using wakePanels_type = std::vector<wakePanel *>;
    using wakePanels_index_type = wakePanels_type::size_type;
    std::vector<bodyPanel*>* bPanels;
    std::vector<wakePanel*>* wPanels;
    Eigen::VectorXd sigmas;
    
    double CL_trefftz = 0;
    double CD_trefftz = 0;
    Eigen::Vector3d Fbody;
    Eigen::Vector3d Fwind;
    Eigen::Vector3d CM; //[roll,pitch,yaw]
    Eigen::VectorXd spanLoc;
    Eigen::VectorXd Cl;
    Eigen::VectorXd Cd;
    
    Eigen::Vector3d dF_dAlpha;
    Eigen::Vector3d dF_dBeta;
    Eigen::Vector3d dM_dAlpha;
    Eigen::Vector3d dM_dBeta;
    
    using streamlines_type = std::vector<bodyStreamline *>;
    using streamlines_index_type = streamlines_type::size_type;

    streamlines_type bStreamlines;
    Eigen::Vector3d windToBody(double V,double alpha,double beta);
    
    Eigen::Vector3d bodyToWind(const Eigen::Vector3d &vec);
    void setSourceStrengths();
    bool solveMatrixEq();
    void compVelocity();
    void trefftzPlaneAnalysis();
    void createStreamlines();
    void stabilityDerivatives();
    void writeVTU(std::string filename);
    void writeFiles();
    void writeBodyData(const std::string & path, const Eigen::MatrixXd &nodeMat);
    void writeWakeData(const std::string &, const Eigen::MatrixXd &nodeMat);
    void writeSpanwiseData(const std::string & path);
    void writeBodyStreamlines(const std::string & path);




    struct meshDat {
        std::vector<Eigen::Vector3d> velocity;
        std::vector<double> coef_press;
        std::vector<Eigen::Vector3d> cellCenter;
    } volMeshDat;
    
    using cells_type = std::vector<Eigen::VectorXi>;
    using cells_index_type = cells_type::size_type;
    std::vector<Eigen::Matrix<size_t, Eigen::Dynamic, 1>> cells;
    Eigen::MatrixXd pts;

    Eigen::Vector3d velocityAtPoint(Eigen::Vector3d POI);
    void createVolMesh();
    void populateVolMesh();
    void writeVolMeshData(const std::string & path, Eigen::MatrixXd &nodeMat, std::vector<Eigen::Matrix<size_t, Eigen::Dynamic, 1>> cells);

    Eigen::MatrixXd solnMat; // For unsteady sims, but needs to be in parent class for simple output



public:
    cpCase(geometry *ggeom, double V, double aalpha, double bbeta, double mmach, inputParams* inParams)
      : geom(ggeom), params(inParams), Vmag(V), mach(mmach), alpha(aalpha), beta(bbeta)
    {
        Vinf = windToBody(V,alpha,beta);
        bPanels = geom->getBodyPanels();
        wPanels = geom->getWakePanels();
        PG = sqrt(1-pow(mach,2));
    }
    
    virtual ~cpCase();
    
    void run(bool printFlag, bool surfStreamFlag, bool stabDerivFlag);
    
    double getMach() {return mach;}
    double getV() {return Vmag;}
    double getAlpha() {return alpha;}
    double getBeta() {return beta;}
    double getCL() {return CL_trefftz;}
    double getCD() {return CD_trefftz;}
    Eigen::Vector3d getMoment() {return CM;}
    Eigen::Vector3d getBodyForces() {return Fbody;}
    Eigen::Vector3d getWindForces() {return Fwind;}
    Eigen::Vector3d get_dF_dAlpha() {return dF_dAlpha;}
    Eigen::Vector3d get_dF_dBeta() {return dF_dBeta;}
    Eigen::Vector3d get_dM_dAlpha() {return dM_dAlpha;}
    Eigen::Vector3d get_dM_dBeta() {return dM_dBeta;}
    Eigen::MatrixXd get_soln_mat() {return solnMat;}
    
};
#endif /* defined(__CPanel__runCase__) */

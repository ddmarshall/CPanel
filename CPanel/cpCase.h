//
//  runCase.h
//  CPanel
//
//  Created by Chris Satterwhite on 10/13/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__runCase__
#define __CPanel__runCase__

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <boost/filesystem/operations.hpp>
#include "geometry.h"
#include "VTUfile.h"
#include "bodyStreamline.h"
#include "cpFile.h"
#include "inputParams.h"



class cpCase
{
    geometry *geom;
//    std::string path;
//    std::string name;
//    cpFile* inFile;
//    cpFile* gFile;
    inputParams* params;
    double Vmag;
    double mach;
    double PG; // Prandtl-Glauert Correction
    double alpha;
    double beta;
//    double Sref;
//    double bref;
//    double cref;
//    Eigen::Vector3d cg;
    Eigen::Vector3d Vinf;
    Eigen::Matrix3d transform; 
    
    
    std::vector<bodyPanel*>* bPanels;
    std::vector<wakePanel*>* wPanels;
    Eigen::VectorXd sigmas;
    
    double CL_trefftz;
    double CD_trefftz;
    Eigen::Vector3d Fbody;
    Eigen::Vector3d Fwind;
    Eigen::Vector3d CM; //[roll,pitch,yaw]
    Eigen::VectorXd spanLoc;
    Eigen::VectorXd Cl;
    Eigen::VectorXd Cd;
    
    std::vector<bodyStreamline*> bStreamlines;
    Eigen::Vector3d windToBody(double V,double alpha,double beta);
    
    Eigen::Vector3d bodyToWind(const Eigen::Vector3d &vec);
    void setSourceStrengths();
    bool solveMatrixEq();
    void compVelocity();
    void trefftzPlaneAnalysis();
    void createStreamlines();
    void writeVTU(std::string filename);
    void writeFiles();
    void writeBodyData(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat);
    void writeWakeData(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat);
    void writeSpanwiseData(boost::filesystem::path path);
    void writeBodyStreamlines(boost::filesystem::path path);
    
public:
//    cpCase(geometry *geom,double V,double alpha,double beta, double mach, cpFile* inputFile, cpFile* geometryFile, std::string path,std::string name) : geom(geom), Vmag(V), alpha(alpha), beta(beta), mach(mach), path(path), name(name),inFile(inputFile),gFile(geometryFile)
//    {
//        Vinf = windToBody(V,alpha,beta);
//        bPanels = geom->getBodyPanels();
//        wPanels = geom->getWakePanels();
//        Sref = geom->getSref();
//        bref = geom->getbref();
//        cref = geom->getcref();
//        cg = geom->getCG();
//        PG = sqrt(1-pow(mach,2));
//    }
    cpCase(geometry *geom,double V, double alpha, double beta, double mach, inputParams* inParams) : geom(geom), Vmag(V), alpha(alpha), beta(beta), mach(mach), params(inParams)
    {
        Vinf = windToBody(V,alpha,beta);
        bPanels = geom->getBodyPanels();
        wPanels = geom->getWakePanels();
        PG = sqrt(1-pow(mach,2));
    }
    
    virtual ~cpCase();
    
    void run();
    
    double getMach() {return mach;}
    double getV() {return Vmag;}
    double getAlpha() {return alpha;}
    double getBeta() {return beta;}
    double getCL() {return CL_trefftz;}
    double getCD() {return CD_trefftz;}
    Eigen::Vector3d getMoment() {return CM;}
    Eigen::Vector3d getBodyForces() {return Fbody;}
    Eigen::Vector3d getWindForces() {return Fwind;}
    
};
#endif /* defined(__CPanel__runCase__) */

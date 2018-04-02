/*******************************************************************************
 * Copyright (c) 2015 Chris Satterwhite
 * Copyright (c) 2016 Connor Sousa
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
 *    Connor Sousa - particle code and implementation
 *    David D. Marshall - misc. improvements
 ******************************************************************************/

#include "influenceTests.h"

void influenceTests::panelVConstTest(Eigen::MatrixXd testingPoints, Eigen::MatrixXd velocityData, bodyPanel* tPan)
{
    // Test panel vertical velocity survey
    Eigen::MatrixXd cpanelVelSurv(testingPoints.rows(),3), survDifference(testingPoints.rows(),3);
    for(int i=0; i<cpanelVelSurv.rows(); i++)
    {
        cpanelVelSurv.row(i) = tPan->panelV(testingPoints.row(i))*4*M_PI;
//        Eigen::Vector3d survpnt = testingPoints.row(i);
//        Eigen::Vector3d pntVel =  tPan->pntVInf(testingPoints.row(i))*4*M_PI;
//        Eigen::Vector3d constVel = tPan->panelV(testingPoints.row(i))*4*M_PI;
//        Eigen::Vector3d expected = velocityData.row(i);
        
        for(int j=0; j<testingPoints.cols(); j++)
        {
            survDifference(i,j) = std::abs(cpanelVelSurv(i,j)-velocityData(i,j))*2/(velocityData(i,j)+cpanelVelSurv(i,j))*100;
        }
    }
    
//    std::cout << "\n\ntesting pts " << testingPoints << std::endl;
//    std::cout << "\n\nvelocity dat" << velocityData  << std::endl;
//    std::cout << "\n\ncpanel Vel  " << cpanelVelSurv << std::endl;
//    
//    std::cout << "\n\nsurv difference" << survDifference << std::endl;
//    
    std::cout << maxCoefficient(survDifference);

    

}

void influenceTests::panelVPntTest(Eigen::MatrixXd testingPoints, Eigen::MatrixXd velocityData, bodyPanel* tPan)
{
    // Test panel vertical velocity survey
    Eigen::MatrixXd cpanelVelSurv(testingPoints.rows(),3), survDifference(testingPoints.rows(),3);
    for(int i=0; i<cpanelVelSurv.rows(); i++)
    {
        cpanelVelSurv.row(i) = tPan->pntVInf(testingPoints.row(i))*4*M_PI;
        Eigen::Vector3d survpnt = testingPoints.row(i);
        Eigen::Vector3d pntVel = tPan->pntVInf(testingPoints.row(i))*4*M_PI;
        Eigen::Vector3d constVel = tPan->panelV(testingPoints.row(i))*4*M_PI;
        Eigen::Vector3d expected = velocityData.row(i);
        
        for(int j=0; j<testingPoints.cols(); j++)
        {
            survDifference(i,j) = std::abs(cpanelVelSurv(i,j)-velocityData(i,j))*2/(velocityData(i,j)+cpanelVelSurv(i,j))*100;
        }
    }
//    std::cout << "\n\ntesting pts " << testingPoints << std::endl;
//    std::cout << "\n\nvelocity dat" << velocityData  << std::endl;
//    std::cout << "\n\ncpanel Vel  " << cpanelVelSurv << std::endl;
//
//    std::cout << "\n\nsurv difference" << survDifference << std::endl;
//

    std::cout << maxCoefficient(survDifference);
    
    
}


double influenceTests::maxCoefficient(Eigen::MatrixXd mat)
{
    double maxCoef = 0;
    for(int i=0; i< mat.rows(); i++){
        for (int j=0; j<mat.cols(); j++) {
            if (isnormal(mat(i,j))) {
                if(std::abs(mat(i,j)) > maxCoef){
                    maxCoef = std::abs(mat(i,j));
                }
            }
        }
    }
    
    return maxCoef;
}


void influenceTests::velocityComparer(DoubVelData* dat)
{
    
    

influenceTests* test = new influenceTests();

// Constant doublet vertical survey
std::cout << "\nVertical const doublet velocity max difference in percent %: " << std::flush;
test->panelVConstTest(dat->vertSurveyPts, dat->vertVelConst, dat->testPan);

// Constant doublet diagonal survey
std::cout << "\nDiagonal const doublet velocity max difference in percent %: " << std::flush;
test->panelVConstTest(dat->diagSurveyPts, dat->diagVelConst, dat->testPan);

// Constant doublet median survey
std::cout << "\nMedian const doublet velocity max difference in percent %: " << std::flush;
test->panelVConstTest(dat->MedSurveyPts, dat->medVelConst, dat->testPan);

// Point doublet vertical survey
std::cout << "\nVertical pnt doublet velocity max difference in percent %: " << std::flush;
test->panelVPntTest(dat->vertSurveyPts, dat->vertVelPnt, dat->testPan);

// Point doublet diagonal survey
std::cout << "\nDiagonal pnt doublet velocity max difference in percent %: " << std::flush;
test->panelVPntTest(dat->diagSurveyPts, dat->diagVelPnt, dat->testPan);

// Point doublet median survey
std::cout << "\nMedian pnt doublet velocity max difference in percent %: " << std::flush;
test->panelVPntTest(dat->MedSurveyPts, dat->medVelPnt, dat->testPan);




}


std::vector<particle*> influenceTests::randomParticleGenerator(int numParts){
    // Function generates random particles with position and strength within 0 and 1 for each component
    std::vector<particle*> particles;
    
    double radius = .5/numParts; // ~ avg. dist btw parts
    for (int i=0; i<numParts; i++) {
        
        // Create random position and strength
        Eigen::Vector3d pos, strength;
        for (int j=0; j<pos.size(); j++)
        {
            pos(j) = pick_a_number(0,1);
            strength(j) = pick_a_number(0,1);
        }
        
        particle* p = new particle(pos,strength,radius, {0,0,0}, {0,0,0});
        particles.push_back(p);
    }
    
    return particles;
}


void influenceTests::FMMtests(){
    
    std::vector<particle*> testParts = randomParticleGenerator(750);
    
//    BarnesHutSpeedTest(testParts);
//    barnesHutTest(testParts);
//    interactionListTest(testParts);
//    speedTestGraphComparison();
    BarnesHutCarpetPlotData(testParts);
}

void influenceTests::barnesHutTest(std::vector<particle*> testParts){
    // This test creates an octree with random particles and then computes the velocity influnce by both the brute force method and then the Barnes Hut algorithm
    
    
    
    // Test params
    Eigen::Vector3d POI = {0.5,0.5,0.5};

    
    // Manually sum up particle influences
    std::cout << "Brute force method..." << std::endl;
    Eigen::Vector3d allPartsVel = Eigen::Vector3d::Zero();
    
    clock_t partBegin = clock();
    for (int i=0; i<testParts.size(); i++)
    {
        allPartsVel+=testParts[i]->partVelInfl(POI);
    }
    clock_t partEnd = clock();
    std::cout << "Individual Particle Calc Time = " <<  (partEnd-partBegin) << " clocks" << std::endl;

    
    
    // Barnes-Hut
    std::cout << "Barnes-Hut method:" << std::endl;
    
    std::cout << "Building octree..." << std::endl;
    clock_t octBuild_start = clock();
        particleOctree partTree;
        partTree.setMaxMembers(1);
        partTree.setMaxTheta(0.1);
        partTree.addData(testParts);
    clock_t octBuild_end = clock();
    std::cout << "Treecode build time = " <<  (octBuild_end-octBuild_start) << " clocks" << std::endl;

    // Writing octree output (not inlcuded in time)
    std::string file_name = "/Users/C_Man/Desktop/CPanelDevelopment/FMM/testTree.txt";
    octreeFile* oct;
    oct = new octreeFile(file_name,&partTree);
    
    
    std::cout << "Calculating expansions..." << std::endl;
    clock_t expBuild_start = clock();
        particleFMM FMM;
        FMM.build(&partTree);
    clock_t expBuild_end = clock();
    std::cout << "Particle Expansions build time = " <<  (expBuild_end-expBuild_start) << " clocks" << std::endl;

    
    std::cout << "Velocity using Barnes-Hut..." << std::endl;
    clock_t bhBegin = clock();
        Eigen::Vector3d barnesHutVel = FMM.barnesHutVel(POI);
    clock_t bhEnd = clock();
    std::cout << "Barnes Hut time = " <<  (bhEnd-bhBegin) << " clocks" << std::endl;


    std::cout << "Brute force method vel = {"<<allPartsVel.x()<<", "<<allPartsVel.y()<<","<<allPartsVel.z()<<"}"<<std::endl;
    std::cout << "Barnes-Hut velocity    = {"<<barnesHutVel.x()<<", "<<barnesHutVel.y()<<","<<barnesHutVel.z()<<"}"<<std::endl;

    double percDiffX = (allPartsVel.x()-barnesHutVel.x())/allPartsVel.x()*100;
    double percDiffY = (allPartsVel.y()-barnesHutVel.y())/allPartsVel.y()*100;
    double percDiffZ = (allPartsVel.z()-barnesHutVel.z())/allPartsVel.z()*100;

    std::cout << "% diff = {" << percDiffX<<"%, "<<percDiffY<<"%, "<<percDiffZ<<"%} "<<std::endl;
    
    std::cout <<  std::endl;
    
}

void influenceTests::BarnesHutCarpetPlotData(std::vector<particle*> testParts)
{
    // Test conditions:
    std::vector<int> maxMembers{1,3,5,7,10,12,15};
    std::vector<double> maxThetas{0.05, 0.1, 0.15, 0.2, 0.25, 0.3};
    
    Eigen::MatrixXd errorMag(maxMembers.size(),maxThetas.size());
    Eigen::MatrixXd boost(maxMembers.size(),maxThetas.size());
    
    
    // ========== Error test ========== //
    Eigen::Vector3d POI = {0.5,0.5,0};

    // Brute force
    Eigen::Vector3d allPartsVel = Eigen::Vector3d::Zero();
    for (int i=0; i<testParts.size(); i++)
    {
        allPartsVel+=testParts[i]->partVelInfl(POI);
    }
    
    // Barnes Hut
    for (int i=0; i<maxMembers.size(); i++)
    {
        for (int j=0; j<maxThetas.size(); j++)
        {
            particleOctree partTree;
            partTree.setMaxMembers(maxMembers[i]);
            partTree.setMaxTheta(maxThetas[j]);
            partTree.addData(testParts);
            
            particleFMM FMM;
            FMM.build(&partTree);
            Eigen::Vector3d barnesHutVel = FMM.barnesHutVel(POI);
            
            Eigen::Vector3d error;
            error.x() = (allPartsVel.x()-barnesHutVel.x())/allPartsVel.x()*100;
            error.y() = (allPartsVel.y()-barnesHutVel.y())/allPartsVel.y()*100;
            error.z() = (allPartsVel.z()-barnesHutVel.z())/allPartsVel.z()*100;
            
//            std::cout << "% diff = {" << error.x()<<"%, "<<error.y()<<"%, "<<error.z()<<"%} "<<std::endl;

            errorMag(i,j) = error.norm();
        }
    }

    
    std::cout << "error = \n" << errorMag << std::endl;

    // ============ Speed Test ============== //
    
    // Brute Force
    clock_t bruteForceBegin = clock();
    std::vector<Eigen::Vector3d> velOnVec;
    for (int i=0; i<testParts.size(); i++)
    {
        Eigen::Vector3d velOn = Eigen::Vector3d::Zero();
        for (int j=0; j<testParts.size(); j++)
        {
            velOn += testParts[j]->partVelInfl(testParts[i]->pos);
        }
        velOnVec.push_back(velOn);
    }
    clock_t bfTime = (clock()-bruteForceBegin);
    
    // Converting to double so can do math with it
    std::string bfTimeString = std::to_string(bfTime);
    std::string::size_type sz;     // alias of size_t
    double bfTimeD = std::stod (bfTimeString, &sz);
    
    
    

    // Barnes Hut
    
    for (int i=0; i<maxMembers.size(); i++)
    {
        for (int j=0; j<maxThetas.size(); j++)
        {
            clock_t barnesHutBegin = clock();
            
            particleOctree partTree;
            partTree.setMaxMembers(maxMembers[i]);
            partTree.setMaxTheta(maxThetas[j]);
            partTree.addData(testParts);
            
            particleFMM FMM;
            FMM.build(&partTree);
            
            std::vector<Eigen::Vector3d> velOnVecBH;
            for (int k=0; k<testParts.size(); k++)
            {
                velOnVecBH.push_back(FMM.barnesHutVel(testParts[k]->pos));
            }
            clock_t bhTime = (clock()-barnesHutBegin);
            
            // Converting clock time to double so can do math with it
            std::string bhTimeString = std::to_string(bhTime);
            double bhTimeD = std::stod (bhTimeString, &sz);
            
            boost(i,j) = (bfTimeD/bhTimeD)-1;
        }
    }
    
    
    std::cout << "Boost = \n" << boost << std::endl;
    
    
    
}

void influenceTests::BarnesHutSpeedTest(std::vector<particle*> testParts){
    // Function goes through each interaction like what would happen to advance a cloud of particles.
    std::cout << "------------Number of test particles = " << testParts.size() << "--------------" << std::endl;
    
    
    // Brute Force
//    std::cout << "Brute force method..." << std::endl;
    clock_t bruteForceBegin = clock();
    std::vector<Eigen::Vector3d> velOnVec;
    for (int i=0; i<testParts.size(); i++)
    {
        Eigen::Vector3d velOn = Eigen::Vector3d::Zero();
        for (int j=0; j<testParts.size(); j++)
        {
            velOn += testParts[j]->partVelInfl(testParts[i]->pos);
        }
        velOnVec.push_back(velOn);
    }
    clock_t bfTime = (clock()-bruteForceBegin);
    std::cout << "Brute force method = " <<  bfTime << " clocks" << std::endl;

    
    
    // Barnes Hut
    clock_t barnesHutBegin = clock();
    
    // Build octree
    particleOctree partTree;
    partTree.setMaxMembers(1);
    partTree.addData(testParts);
    
    // Build Expansions
    particleFMM FMM;
    FMM.build(&partTree);
    
    std::vector<Eigen::Vector3d> velOnVecBH;
//    clock_t barnesHutBegin = clock();
    for (int i=0; i<testParts.size(); i++)
    {
//        clock_t oneBH = clock();
        velOnVecBH.push_back(FMM.barnesHutVel(testParts[i]->pos));
//        clock_t oneBHend = clock();
//        std::cout << "one BH calc: " << oneBHend-oneBH << std::endl;
    }
    clock_t bhTime = (clock()-barnesHutBegin);
    std::cout << "Barnes Hut method  = " << bhTime  << " clocks" << std::endl;

    // Converting clock time to double so can do math with it
    std::string bhTimeString = std::to_string(bhTime);
    std::string::size_type sz;     // alias of size_t
    double bhTimeD = std::stod (bhTimeString, &sz);
    
    std::string bfTimeString = std::to_string(bfTime);
    double bfTimeD = std::stod (bfTimeString, &sz);
    
    std::cout << "Speed up  = " << (bfTimeD/bhTimeD)-1 << " times faster" << std::endl;
}



void influenceTests::speedTestGraphComparison()
{
    int a = 500;
    int b = 10000;
    int n = 20;
    std::vector<int> numParts = linspace(a, b, n);
    
    for (int i=0; i<numParts.size(); i++)
    {
        std::vector<particle*> testParts = randomParticleGenerator(numParts[i]);
        BarnesHutSpeedTest(testParts);
    }

}


std::vector<int> influenceTests::linspace(int a, int b, int n)
{
    std::vector<int> linSpaced;
    
    int step = (b-a)/(n-1);
    
    for (int i=0; i<n; i++) {
        linSpaced.push_back(a+step*i);
    }
    return linSpaced;
}

std::vector<double> influenceTests::linspace(double a, double b, int n)
{
    std::vector<double> linSpaced;
    
    double step = (b-a)/(n-1);
    
    for (int i=0; i<n; i++) {
        linSpaced.push_back(a+step*i);
    }
    return linSpaced;
}

std::default_random_engine & global_urng( )
{
    static std::default_random_engine u{};
    return u;
}

double influenceTests::pick_a_number( double from, double upto )
{
    static std::uniform_real_distribution<> d{};
    using parm_t = decltype(d)::param_type;
    
    return d( global_urng(), parm_t{from, upto} );
}

void influenceTests::potentialComparer(bodyPanel* pan)
{
    // comparing 4 sided panel potential in CPanel to Katz
    
    Eigen::Vector3d POI;
    POI << 1.2, 0.5 , 0.75;
    
    double phi = pan->dubPhiInf(POI);
    
    std::cout << "CPanel Phi = "  << phi << std::endl;
    
    cpNode* n1 = pan->getNodes()[0];
    cpNode* n2 = pan->getNodes()[1];
    cpNode* n3 = pan->getNodes()[2];
    cpNode* n4 = pan->getNodes()[3];
    
    double x,y,z;
    x = POI.x();
    y = POI.y();
    z = POI.z();
    
    double x1,x2,x3,x4;
    x1 = n1->getPnt().x();
    x2 = n2->getPnt().x();
    x3 = n3->getPnt().x();
    x4 = n4->getPnt().x();
    
    double y1,y2,y3,y4;
    y1 = n1->getPnt().y();
    y2 = n2->getPnt().y();
    y3 = n3->getPnt().y();
    y4 = n4->getPnt().y();
    
    double z1,z2,z3,z4;
    z1 = n1->getPnt().z();
    z2 = n2->getPnt().z();
    z3 = n3->getPnt().z();
    z4 = n4->getPnt().z();
    
    double m12, m23, m34, m41;
    m12 = (y2 - y1)/(x2 - x1);
    m23 = (y3 - y2)/(x3 - x2);
    m34 = (y4 - y3)/(x4 - x3);
    m41 = (y1 - y4)/(x1 - x4);
    
    double r1, r2, r3, r4;
    r1 = pow((x-x1)*(x-x1) + (y-y1)*(y-y1) + (z)*(z) , 0.5);
    r2 = pow((x-x2)*(x-x2) + (y-y2)*(y-y2) + (z)*(z) , 0.5);
    r3 = pow((x-x3)*(x-x3) + (y-y3)*(y-y3) + (z)*(z) , 0.5);
    r4 = pow((x-x4)*(x-x4) + (y-y4)*(y-y4) + (z)*(z) , 0.5);
    
    double e1, e2, e3, e4;
    e1 = (x-x1)*(x-x1) + z*z;
    e2 = (x-x2)*(x-x2) + z*z;
    e3 = (x-x3)*(x-x3) + z*z;
    e4 = (x-x4)*(x-x4) + z*z;
    
    double h1, h2, h3, h4;
    h1 = (x-x1)*(y-y1);
    h2 = (x-x2)*(y-y2);
    h3 = (x-x3)*(y-y3);
    h4 = (x-x4)*(y-y4);
    
    
    double KatzPhi = pan->getMu()/4/M_PI * (atan((m12*e1 - h1)/(z*r1)) - atan((m12*e2 - h2)/(z*r2)) +
                                     atan((m23*e2 - h2)/(z*r2)) - atan((m23*e3 - h3)/(z*r3)) +
                                     atan((m34*e3 - h3)/(z*r3)) - atan((m34*e4 - h4)/(z*r4)) +
                                     atan((m41*e4 - h4)/(z*r4)) - atan((m41*e1 - h1)/(z*r1)));

    
    std::cout << "Katz Phi = "  << KatzPhi << std::endl;

    
    
    
}




























//influenceTests::influenceTests()
//{
//    createTestPanel();
//        TEST_ADD(influenceTests::testPntClose)
//        TEST_ADD(influenceTests::testPntFar)
////        TEST_ADD(influenceTests::testCollocationPnt)
//    //    TEST_ADD(influenceTests::testDoublet)
//    //    TEST_ADD(influenceTests::testSource)
//}
//
//void influenceTests::createTestPanel()  //Connor: I think this is old and doesn't work with new panels
//{
//    testNodes.resize(3,3);
//    testNodes(0,0) = 1;
//    testNodes(0,1) = 0;
//    testNodes(0,2) = 0;
//    testNodes(1,0) = -1;
//    testNodes(1,1) = 1;
//    testNodes(1,2) = 0;
//    testNodes(2,0) = -1;
//    testNodes(2,1) = -1;
//    testNodes(2,2) = 0;
//    Eigen::Vector3i verts;
//    verts << 0,1,2;
//    Eigen::Vector3d bezNorm;
//    bezNorm << 0,0,1;
//    std::vector<edge*> edges;
//    testPan = new bodyPanel(&testNodes,edges,bezNorm,1,false);
//}
//
//
//void influenceTests::testPntClose()
//{
//    Eigen::Vector3d POI;
//    POI << 2,2,2; // |POI|/(side length) ~ 1.2 so panel formulation will be used;
//
//    double phiSrc,phiDub;
//    Eigen::Vector3d vSrc,vDub;
//    phiSrc = 0;
//    phiDub = 0;
//    vSrc = Eigen::Vector3d::Zero();
//    vDub = Eigen::Vector3d::Zero();
//    testPan->panelPhiInf(POI, phiSrc, phiDub);
//    testPan->panelVInf(POI, vSrc, vDub);
//    TEST_ASSERT(isEqual(phiSrc,0.043411,5))
//    TEST_ASSERT(isEqual(vSrc(0),0.007363,5))
//    TEST_ASSERT(isEqual(vSrc(1),0.006451,5))
//    TEST_ASSERT(isEqual(vSrc(2),0.006645,5))
//
//    TEST_ASSERT(isEqual(phiDub,0.006645,5))
//    TEST_ASSERT(isEqual(vDub(0),-0.003381,5))
//    TEST_ASSERT(isEqual(vDub(1),-0.003026,5))
//    TEST_ASSERT(isEqual(vDub(2),0.000158,5))
//}
//
//void influenceTests::testPntFar()
//{
//    Eigen::Vector3d POI;
//    POI << 7,7,7; // |POI|/(side length) ~ 5.5 so far field approximation will be used;
//
//    double phiSrc,phiDub;
//    Eigen::Vector3d vSrc,vDub;
//    phiSrc = 0;
//    phiDub = 0;
//    vSrc = Eigen::Vector3d::Zero();
//    vDub = Eigen::Vector3d::Zero();
//    testPan->panelPhiInf(POI, phiSrc, phiDub);
//    testPan->panelVInf(POI, vSrc, vDub);
//    TEST_ASSERT(isEqual(phiSrc,0.012919,5))
//    TEST_ASSERT(isEqual(vSrc(0),0.000625,5))
//    TEST_ASSERT(isEqual(vSrc(1),0.000596,5))
//    TEST_ASSERT(isEqual(vSrc(2),0.000596,5))
//    
//    TEST_ASSERT(isEqual(phiDub,0.0005958,5))
//    TEST_ASSERT(isEqual(vDub(0),-0.000086,6))
//    TEST_ASSERT(isEqual(vDub(1),-0.000082,6))
//    TEST_ASSERT(isEqual(vDub(2),0.0000027,6))
//}

//void influenceTests::testCollocationPnt()
//{
//    Eigen::Vector3d POI;
//    POI << -1.0/3,0,0;
//
//    double phiSource;
//    Eigen::Vector3d vSource;
//    phiSource = p.sourcePhi(1,POI,vertsLocal);
//    vSource = p.sourceV(1,POI,vertsLocal);
//    TEST_ASSERT(phiSource == 0)
//    TEST_ASSERT(vSource(0) == 0)
//    TEST_ASSERT(vSource(1) == 0)
//    TEST_ASSERT(vSource(2) == 0.5)
//
//    double phiDoublet;
//    Eigen::Vector3d vDoublet;
//    phiDoublet = p.doubletPhi(1,POI,vertsLocal);
//    vDoublet = p.doubletV(1,POI,vertsLocal);
//    TEST_ASSERT(phiDoublet == -0.5)
//    TEST_ASSERT(vDoublet(0) == 0)
//    TEST_ASSERT(vDoublet(1) == 0)
//    TEST_ASSERT(isEqual(vDoublet(2),1.432394487,9))
//}

//bool influenceTests::isEqual(double var1, double var2, int decimalPrecision)
//{
//    double eps = pow(10,-decimalPrecision);
//    return (std::abs(var1-var2)<eps);
//    
//}

//void influenceTests::testDoublet()
//{
//    std::string path = "/Users/Chris/Desktop/Thesis/Code/Geometry and Solution Files/";
//    std::string filename = "doubletCheck.txt";
//    std::ifstream fid;
//    fid.open(path+filename);
//    if (fid.is_open())
//    {
//        Eigen::MatrixXd nodes;
//        nodes.resize(3,3);
//        Eigen::Vector3i verts;
//        verts << 0,1,2;
//        Eigen::Vector3d col;
//        Eigen::MatrixXd pnts(9,3);
//        Eigen::VectorXd pot(9);
//        Eigen::MatrixXd vel(9,3);
//
//        for (int i=0; i<3; i++)
//        {
//            fid >> nodes(i,0) >> nodes(i,1) >> nodes(i,2);
//        }
//        fid >> col(0) >> col(1) >> col(2);
//        for (int i=0; i<9; i++)
//        {
//            fid >> pnts(i,0) >> pnts(i,1) >> pnts(i,2);
//        }
//        for (int i=0; i<9; i++)
//        {
//            fid >> pot(i);
//        }
//        for (int i=0; i<9; i++)
//        {
//            fid >> vel(i,0) >> vel(i,1) >> vel(i,2);
//        }
//        vel *= -1; //Different sign convention
//        Eigen::VectorXd phiCalc(9);
//        Eigen::MatrixXd Vcalc(9,3);
//        bodyPanel pan(verts,&nodes,1);
//        Eigen::Matrix3d localVerts = pan.getLocalVerts();
//
//        for (int i=0; i<9; i++)
//        {
//            phiCalc(i) = pan.doubletPhi(1,pnts.row(i),localVerts);
//            TEST_ASSERT(isEqual(phiCalc(i),pot(i),5))
//            Vcalc.row(i) = pan.doubletV(true, pnts.row(i));
//            for (int j=0; j<3; j++)
//            {
//                std::cout << Vcalc(i,j) << "," << vel(i,j) << std::endl;
//                TEST_ASSERT(isEqual(Vcalc(i,j),vel(i,j),5))
//            }
//        }
//    }
//    else
//    {
//        std::cout << "ERROR : File not found" << std::endl;
//        exit(EXIT_FAILURE);
//    }
//
//}
//
//void influenceTests::testSource()
//{
//    std::string path = "/Users/Chris/Desktop/Thesis/Code/Geometry and Solution Files/";
//    std::string filename = "sourceCheck.txt";
//    std::ifstream fid;
//    fid.open(path+filename);
//    if (fid.is_open())
//    {
//        Eigen::MatrixXd nodes;
//        nodes.resize(3,3);
//        Eigen::Vector3i verts;
//        verts << 0,1,2;
//        Eigen::Vector3d col;
//        Eigen::MatrixXd pnts(9,3);
//        Eigen::VectorXd pot(9);
//        Eigen::MatrixXd vel(9,3);
//
//        for (int i=0; i<3; i++)
//        {
//            fid >> nodes(i,0) >> nodes(i,1) >> nodes(i,2);
//        }
//        fid >> col(0) >> col(1) >> col(2);
//        for (int i=0; i<9; i++)
//        {
//            fid >> pnts(i,0) >> pnts(i,1) >> pnts(i,2);
//        }
//        for (int i=0; i<9; i++)
//        {
//            fid >> pot(i);
//        }
//        for (int i=0; i<9; i++)
//        {
//            fid >> vel(i,0) >> vel(i,1) >> vel(i,2);
//        }
//        pot *= -1; //Different sign convention;
//        Eigen::VectorXd phiCalc(9);
//        Eigen::MatrixXd Vcalc(9,3);
//        bodyPanel pan(verts,&nodes,1);
//        Eigen::Matrix3d localVerts = pan.getLocalVerts();
//
//        for (int i=0; i<9; i++)
//        {
//            phiCalc(i) = pan.sourcePhi(1,pnts.row(i),localVerts);
//            TEST_ASSERT(isEqual(phiCalc(i),pot(i),5))
//            Vcalc.row(i) = pan.sourceV(true, pnts.row(i));
//            for (int j=0; j<3; j++)
//            {
//                TEST_ASSERT(isEqual(-Vcalc(i,j),vel(i,j),5))
//            }
//        }
//    }
//    else
//    {
//        std::cout << "ERROR : File not found" << std::endl;
//        exit(EXIT_FAILURE);
//    }
//
//}

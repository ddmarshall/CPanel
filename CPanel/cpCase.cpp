//
//  runCase.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 10/13/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "cpCase.h"

cpCase::~cpCase()
{
    for (int i=0; i<bStreamlines.size(); i++)
    {
        delete bStreamlines[i];
    }
}

void cpCase::run(bool printFlag, bool surfStreamFlag, bool stabDerivFlag)
{
    bool converged;
    std::string check = "\u2713";
    setSourceStrengths();
    converged = solveMatrixEq();
    if (printFlag)
    {
        std::cout << std::setw(17) << std::left << check << std::flush;
    }
    
    compVelocity();
    
    if (printFlag)
    {
        std::cout << std::setw(17) << std::left << check << std::flush;
    }
    
    trefftzPlaneAnalysis();
    
    if (printFlag)
    {
        std::cout << std::setw(18) << std::left << check << std::flush;
    }
    
    if (surfStreamFlag)
    {
        createStreamlines();
        if (printFlag)
        {
            std::cout << std::setw(14) << std::left << check << std::flush;
        }
    }
    else
    {
        if (printFlag)
        {
            std::cout << std::setw(14) << std::left << "X" << std::flush;
        }
    }
    
    if (stabDerivFlag)
    {
        stabilityDerivatives();
        if (printFlag)
        {
            std::cout << std::setw(24) << std::left << check << std::flush;
        }
    }
    else
    {
        if (printFlag)
        {
            std::cout << std::setw(24) << std::left << "X" << std::flush;
        }
    }
    
    if (params->volMeshFlag) {
        createVolMesh();
        populateVolMesh();
        if (printFlag)
        {
            std::cout << std::setw(23) << std::left << check << std::endl;
        }
    }
    else
    {
        if (printFlag)
        {
            std::cout << std::setw(23) << std::left << "X" << std::endl;
        }
    }
    
    
    if (!converged && printFlag)
    {
        std::cout << "*** Warning : Solution did not converge ***" << std::endl;;
    }
    
    if (printFlag)
    {
        writeFiles();
    }
}

Eigen::Vector3d cpCase::windToBody(double V, double alpha, double beta)
{
    alpha *= M_PI/180;
    beta *= M_PI/180;
    
    Eigen::Matrix3d T;
    T(0,0) = cos(alpha)*cos(beta);
    T(0,1) = cos(alpha)*sin(beta);
    T(0,2) = -sin(alpha);
    T(1,0) = -sin(beta);
    T(1,1) = cos(beta);
    T(1,2) = 0;
    T(2,0) = sin(alpha)*cos(beta);
    T(2,1) = sin(alpha)*sin(beta);
    T(2,2) = cos(alpha);
    Eigen::Vector3d Vt;
    Eigen::Vector3d Vel;
    Vt << V,0,0;
    transform = T;
    
    Vel = transform*Vt;
    
    return Vel;
}

Eigen::Vector3d cpCase::bodyToWind(const Eigen::Vector3d &vec)
{
    Eigen::Matrix3d T = transform.transpose();
    return T*vec;
}


void cpCase::setSourceStrengths()
{
    sigmas.resize(bPanels->size());
    for (int i=0; i<bPanels->size(); i++)
    {
        (*bPanels)[i]->setSigma(Vinf,0);
        sigmas(i) = (*bPanels)[i]->getSigma();
    }
}


bool cpCase::solveMatrixEq()
{
    bool converged = true;
    
    // Solve matrix equations and set potential for all panels;
    Eigen::MatrixXd* A = geom->getA();
    Eigen::MatrixXd* B = geom->getB();
    Eigen::VectorXd RHS = -(*B)*sigmas;
    Eigen::VectorXd doubletStrengths(bPanels->size());
    
    
    Eigen::BiCGSTAB<Eigen::MatrixXd> res;
    res.compute((*A));
    doubletStrengths = res.solve(RHS);
    if (res.error() > pow(10,-10))
    {
        converged = false;
    }
    
    for (int i=0; i<bPanels->size(); i++)
    {
        (*bPanels)[i]->setMu(doubletStrengths(i));
        (*bPanels)[i]->setPotential(Vinf);
    }
    for (int i=0; i<wPanels->size(); i++)
    {
        (*wPanels)[i]->setMu();
        (*wPanels)[i]->setPotential(Vinf);
    }
    return converged;
}



void cpCase::compVelocity()
{
    //  Velocity Survey with known doublet and source strengths
    CM.setZero();
    Eigen::Vector3d moment;
    Fbody = Eigen::Vector3d::Zero();
    bodyPanel* p;
    for (int i=0; i<bPanels->size(); i++)
    {
        p = (*bPanels)[i];
        p->computeVelocity(PG,Vinf);
        p->computeCp(Vmag);
        Fbody += -p->getCp()*p->getArea()*p->getBezNormal()/params->Sref;
        moment = p->computeMoments(params->cg);
        CM(0) += moment(0)/(params->Sref*params->bref);
        CM(1) += moment(1)/(params->Sref*params->cref);
        CM(2) += moment(2)/(params->Sref*params->bref);
    }
    Fwind = bodyToWind(Fbody);
}


void cpCase::trefftzPlaneAnalysis()
{
    std::vector<wake*> wakes = geom->getWakes();
    CL_trefftz = 0;
    CD_trefftz = 0;
    for (int i=0; i<wakes.size(); i++)
    {
        wakes[i]->trefftzPlane(Vmag,params->Sref);
        CL_trefftz += wakes[i]->getCL()/PG;
        CD_trefftz += wakes[i]->getCD()/pow(PG,2);
    }
}


void cpCase::createStreamlines()
{
    // Gather TE Panels
    std::vector<surface*> surfs = geom->getSurfaces();
    std::vector<std::pair<Eigen::Vector3d,bodyPanel*>> streamPnts;
    bodyStreamline* s;
    
    for (int i=0; i<surfs.size(); i++)
    {
        streamPnts = surfs[i]->getStreamlineStartPnts(Vinf,PG); // will need to modify this in order to find the stagnation points for a moving body...
        for (int j=0; j<streamPnts.size(); j++)
        {
            s = new bodyStreamline(std::get<0>(streamPnts[j]),std::get<1>(streamPnts[j]),Vinf,PG,geom,3,false);
            bStreamlines.push_back(s);
        }
    }
    
    //    // Off Body Streamline (Temporary)
    //    std::vector<streamline*> sLines;
    //    streamline* sline;
    //    Eigen::Vector3d start;
    //    int ny = 1;
    //    int nz = 6;
    //    double dz = 0.1;
    //    for (int i=0; i<ny; i++)
    //    {
    //        for (int j=0; j<nz; j++)
    //        {
    ////            start << -1,-0.5*params->bref+params->bref*i/(ny-1),-0.35+j*dz;
    //            start << -1,-2.85,-0.45+j*dz;
    //            sline = new streamline(start,20,0.0001,Vinf,PG,geom);
    //            sLines.push_back(sline);
    //        }
    //    }
    //
    //    piece p;
    //    std::vector<piece> pieces;
    //    std::vector<pntDataArray> data;
    //    pntDataArray vel("Velocity");
    //    Eigen::MatrixXi con;
    //    Eigen::MatrixXd pntMat;
    //    std::vector<Eigen::Vector3d> pnts,velocities;
    //
    //    for (int i=0; i<sLines.size(); i++)
    //    {
    //        pnts = sLines[i]->getPnts();
    //        velocities = sLines[i]->getVelocities();
    //        vel.data.resize(velocities.size(),3);
    //        pntMat.resize(pnts.size(),3);
    //        con.resize(pnts.size()-1,2);
    //        for (int j=0; j<pnts.size(); j++)
    //        {
    //            pntMat.row(j) = pnts[j];
    //            vel.data.row(j) = velocities[j];
    //            if (j<con.rows())
    //            {
    //                con(j,0) = j;
    //                con(j,1) = j+1;
    //            }
    //        }
    //        data.push_back(vel);
    //        p.pnts = pntMat;
    //        p.connectivity = con;
    //        p.pntData = data;
    //
    //        pieces.push_back(p);
    //        data.clear();
    //    }
    //
    //    boost::filesystem::path path = boost::filesystem::current_path();
    //
    //    std::string fname = path.string()+"/OffBodyStreamlines.vtu";
    //    VTUfile wakeFile(fname,pieces);
    //
    //    for (int i=0; i<sLines.size(); i++)
    //    {
    //        delete sLines[i];
    //    }
}

void cpCase::stabilityDerivatives()
{
    double delta = 0.5;
    double dRad = delta*M_PI/180;
    cpCase dA(geom,Vmag,alpha+delta,beta,mach,params);
    cpCase dB(geom,Vmag,alpha,beta+delta,mach,params);
    
    dA.run(false,false,false);
    dB.run(false,false,false);
    
    Eigen::Vector3d FA = dA.getWindForces();
    FA(2) = dA.getCL();
    FA(0) = dA.getCD();
    
    Eigen::Vector3d FB = dB.getWindForces();
    FB(2) = dB.getCL();
    FB(0) = dB.getCD();
    
    Eigen::Vector3d F = Fwind;
    F(2) = CL_trefftz;
    F(0) = CD_trefftz;
    
    // Finish calculating stability derivatives
    dF_dAlpha = (FA-F)/dRad;
    dF_dBeta = (FB-F)/dRad;
    
    dM_dAlpha = (dA.getMoment()-CM)/dRad;
    dM_dBeta = (dB.getMoment()-CM)/dRad;
    
    
}

void cpCase::writeFiles()
{
    std::stringstream caseLabel;
    caseLabel << "/V" << Vmag << "_Mach" << mach << "_alpha" << alpha << "_beta" << beta;
    boost::filesystem::path subdir = boost::filesystem::current_path().string()+caseLabel.str();
    if (!boost::filesystem::exists(subdir))
    {
        boost::filesystem::create_directories(subdir);
    }
    Eigen::MatrixXd nodeMat = geom->getNodePnts();
    writeBodyData(subdir,nodeMat);
    if (geom->getWakes().size() > 0)
    {
        writeWakeData(subdir,nodeMat);
        writeSpanwiseData(subdir);
    }
    
    
    if (params->volMeshFlag)
    {
        writeVolMeshData(subdir, pts, cells);
    }
    
    
    if (params->surfStreamFlag)
    {
        writeBodyStreamlines(subdir);
    }
    
}

void cpCase::writeBodyData(boost::filesystem::path path,const Eigen::MatrixXd &nodeMat)
{
    std::vector<cellDataArray> data;
    cellDataArray mu("Doublet Strengths"),sigma("Source Strengths"),pot("Velocity Potential"),V("Velocity"),Cp("Cp"),bN("bezNormals"),x("xPosition"),y("yPosition"),z("zPostition");
    Eigen::MatrixXi con(bPanels->size(),3);
    mu.data.resize(bPanels->size(),1);
    sigma.data.resize(bPanels->size(),1);
    pot.data.resize(bPanels->size(),1);
    V.data.resize(bPanels->size(),3);
    Cp.data.resize(bPanels->size(),1);
    bN.data.resize(bPanels->size(),3);
    x.data.resize(bPanels->size(),1);
    y.data.resize(bPanels->size(),1);
    z.data.resize(bPanels->size(),1);
    
    for (int i=0; i<bPanels->size(); i++)
    {
        mu.data(i,0) = (*bPanels)[i]->getMu();
        sigma.data(i,0) = (*bPanels)[i]->getSigma();
        pot.data(i,0) = (*bPanels)[i]->getPotential();
        V.data.row(i) = (*bPanels)[i]->getGlobalV();
        Cp.data(i,0) = (*bPanels)[i]->getCp();
        con.row(i) = (*bPanels)[i]->getVerts();
        bN.data.row(i) = (*bPanels)[i]->getBezNormal();
        x.data(i,0) = (*bPanels)[i]->getCenter().x();
        y.data(i,0) = (*bPanels)[i]->getCenter().y();
        z.data(i,0) = (*bPanels)[i]->getCenter().z();
    }
    
    data.push_back(mu);
    data.push_back(sigma);
    data.push_back(pot);
    data.push_back(V);
    data.push_back(Cp);
    data.push_back(bN);
    data.push_back(x);
    data.push_back(y);
    data.push_back(z);
    
    piece body;
    body.pnts = nodeMat;
    body.connectivity = con;
    body.cellData = data;
    
    std::string fname = path.string()+"/surfaceData.vtu";
    VTUfile bodyFile(fname,body);
}

void cpCase::writeWakeData(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat)
{
    std::vector<cellDataArray> data;
    cellDataArray mu("Doublet Strengths"),pot("Velocity Potential");
    Eigen::MatrixXi con;
    con.resize(wPanels->size(),3);
    mu.data.resize(wPanels->size(),1);
    pot.data.resize(wPanels->size(),1);
    for (int i=0; i<wPanels->size(); i++)
    {
        mu.data(i,0) = (*wPanels)[i]->getMu();
        pot.data(i,0) = (*wPanels)[i]->getPotential();
        con.row(i) = (*wPanels)[i]->getVerts();
    }
    data.push_back(mu);
    data.push_back(pot);
    
    piece wake;
    wake.pnts = nodeMat;
    wake.connectivity = con;
    wake.cellData = data;
    std::string fname = path.string()+"/wakeData.vtu";
    VTUfile wakeFile(fname,wake);
}


void cpCase::writeSpanwiseData(boost::filesystem::path path)
{
    std::vector<wake*> wakes = geom->getWakes();
    for (int i=0; i<wakes.size(); i++)
    {
        Eigen::VectorXd spanLoc,Cl,Cd;
        spanLoc = 2*wakes[i]->getSpanwisePnts()/params->bref;
        Cl = wakes[i]->getSpanwiseCl()/PG;
        Cd = wakes[i]->getSpanwiseCd()/pow(PG,2);
        
        std::stringstream ss;
        ss << path.string() << "/spanwiseData_Wake" << i+1 << ".csv";
        std::string fname = ss.str();
        std::ofstream fout;
        fout.open(fname);
        if (fout)
        {
            fout << "2y/b,Cl,Cdi" << std::endl;
            for (int i=0; i<spanLoc.size(); i++)
            {
                fout << spanLoc(i) << "," << Cl(i) << "," << Cd(i) << std::endl;
            }
        }
        fout.close();
    }
}

void cpCase::writeBodyStreamlines(boost::filesystem::path path)
{
    piece p;
    std::vector<piece> pieces;
    std::vector<pntDataArray> data;
    pntDataArray vel("Velocity");
    Eigen::MatrixXi con;
    Eigen::MatrixXd pntMat;
    std::vector<Eigen::Vector3d> pnts,velocities;
    
    for (int i=0; i<bStreamlines.size(); i++)
    {
        pnts = bStreamlines[i]->getPnts();
        velocities = bStreamlines[i]->getVelocities();
        vel.data.resize(velocities.size(),3);
        pntMat.resize(pnts.size(),3);
        con.resize(pnts.size()-1,2);
        for (int j=0; j<pnts.size(); j++)
        {
            pntMat.row(j) = pnts[j];
            vel.data.row(j) = velocities[j];
            if (j<con.rows())
            {
                con(j,0) = j;
                con(j,1) = j+1;
            }
        }
        data.push_back(vel);
        p.pnts = pntMat;
        p.connectivity = con;
        p.pntData = data;
        
        pieces.push_back(p);
        data.clear();
    }
    
    std::string fname = path.string()+"/streamlines.vtu";
    VTUfile wakeFile(fname,pieces);
}


void cpCase::writeVolMeshData(boost::filesystem::path path, Eigen::MatrixXd &nodeMat, std::vector<Eigen::VectorXi> cells){
    int nCells = (int)cells.size();
    
    std::vector<cellDataArray> data;
    cellDataArray vel("Velocity"), Cp("Cp");
    Eigen::MatrixXi con;
    con.resize(nCells,8);
    vel.data.resize(nCells,3);
    Cp.data.resize(nCells,1);
    for (int i=0; i<nCells; i++)
    {
        vel.data.row(i) = volMeshDat.velocity[i];
        Cp.data(i,0) = volMeshDat.coef_press[i];
        con.row(i) = cells[i];
    }
    data.push_back(vel);
    data.push_back(Cp);
    
    piece volMesh;
    volMesh.pnts = nodeMat;
    volMesh.connectivity = con;
    volMesh.cellData = data;
    std::string fname = path.string()+"/volumeMesh.vtu";
    VTUfile wakeFile(fname,volMesh);
}


Eigen::Vector3d cpCase::velocityAtPoint(Eigen::Vector3d POI){
    
    Eigen::Vector3d vel = Eigen::Vector3d::Zero();
    
    for(int i=0; i<(*bPanels).size(); i++){
        vel += (*bPanels)[i]->panelV(POI);
    }
    
    for (int i=0; i<(*wPanels).size(); i++) {
        vel += (*wPanels)[i]->panelV(POI);
    }
    
    return vel;
    
}


void cpCase::createVolMesh(){
    
    // Retrieve mesh limits
    double x0, xf, y0, yf, z0, zf;
    x0 = params->volMeshBounds[0];
    xf = params->volMeshBounds[1];
    y0 = params->volMeshBounds[2];
    yf = params->volMeshBounds[3];
    z0 = params->volMeshBounds[4];
    zf = params->volMeshBounds[5];
    
    // Retrieve mesh resolution
    int nX, nY, nZ;
    nX = params->volMeshRes[0];
    nY = params->volMeshRes[1];
    nZ = params->volMeshRes[2];
    
    double hx, hy, hz;
    hx = (xf-x0)/nX;
    hy = (yf-y0)/nY;
    hz = (zf-z0)/nZ;
    
    // Add one for number of points (nodes)
    int nXp = nX+1;
    int nYp = nY+1;
    int nZp = nZ+1;
    int numPts = nXp * nYp * nZp;
    
    // Creating pnt vector
    pts.resize(numPts, 3);
    int count = 0;
    for (int k=0; k<nZp; k++) {
        for (int j=0; j<nYp; j++) {
            for (int i=0; i<nXp; i++) {
                pts(count,0) = x0 + hx*i;
                pts(count,1) = y0 + hy*j;
                pts(count,2) = z0 + hz*k;
                count++;
            }
        }
    }
    
    
    // Create cells
    for (int i=0; i<nX; i++) {
        for (int j=0; j<nY; j++) {
            for (int k=0; k<nZ; k++) {
                Eigen::VectorXi cell_pts(8);
                cell_pts[0] = (k*(nXp*nYp) + j*(nXp) + i);
                cell_pts[1] = (k*(nXp*nYp) + j*(nXp) + i) + 1;
                
                cell_pts[2] = (k*(nXp*nYp) + (j+1)*(nXp) + i);
                cell_pts[3] = (k*(nXp*nYp) + (j+1)*(nXp) + i) + 1;
                
                cell_pts[4] = ((k+1)*(nXp*nYp) + j*(nXp) + i);
                cell_pts[5] = ((k+1)*(nXp*nYp) + j*(nXp) + i) + 1;
                
                cell_pts[6] = ((k+1)*(nXp*nYp) + (j+1)*(nXp) + i);
                cell_pts[7] = ((k+1)*(nXp*nYp) + (j+1)*(nXp) + i) + 1;
                
                cells.push_back(cell_pts);
                
                // Find cell center
                Eigen::MatrixXd cellCorners(8,3);
                for (int m=0; m<cellCorners.rows(); m++) {
                    cellCorners.row(m) = pts.row(cell_pts[m]);
                }
                
                Eigen::Vector3d center;
                for (int m=0; m<cellCorners.cols(); m++) {
                    center(m) = cellCorners.col(m).mean();
                }
                volMeshDat.cellCenter.push_back(center);
            }
        }
    }
    
}

void cpCase::populateVolMesh(){
    
    // Clear Past Timestep Mesh
    volMeshDat.velocity.clear();
    volMeshDat.coef_press.clear();
    
    for (int i=0; i<cells.size(); i++) {
        Eigen::Vector3d velInCell = velocityAtPoint(volMeshDat.cellCenter[i]) + Vinf;
        volMeshDat.velocity.push_back(velInCell);
    }
    
    for (int i=0; i<cells.size(); i++) {
        // Incompressible Bernoulli equation
        double Cp = pow( volMeshDat.velocity[i].norm() / Vmag , 2 );
        volMeshDat.coef_press.push_back(Cp);
    }
    
}

}


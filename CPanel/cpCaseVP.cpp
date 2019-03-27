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

#include "cpCaseVP.h"

#include <iostream>
#include <iomanip>

#include <boost/filesystem/operations.hpp>


void cpCaseVP::run(bool printFlag, bool surfStreamFlag, bool stabDerivFlag){

    if(unsteady){
        readBodyKinFile();
        solnMat.resize(numSteps, 10);
    }

    bool matrix_convergence = false;
    std::string check = "\u2713";
    std::setprecision(5);
    if(printFlag) std::cout << std::setw(11) << " Time step " << std::to_string(timestep) + "/" + std::to_string(numSteps) + ".  " << "Flow time = " << std::setw(5) << std::to_string(timestep*dt) + ". " << std::endl;

    // First time step (one row of wake panels)
    setSourceStrengths();
    matrix_convergence = solveMatrixEq();
    if (unsteady) compVelocity();
    writeFilesVP();
    //    moveGeometry();


    // Initially convect buffer wake and solve
    timestep++;
    if(printFlag) std::cout << std::setw(11) << " Time step " << std::to_string(timestep) + "/" + std::to_string(numSteps) + ".  " << "Flow time = " << std::setw(5) << std::to_string(timestep*dt) + ". " << std::endl;

    convectBufferWake();
    matrix_convergence = solveVPmatrixEq();
    if (unsteady) compVelocity();
    writeFilesVP();
    //    moveGeometry();


    // For rest of timesteps
    while (timestep < numSteps){

        timestep++;
        if(printFlag) std::cout << std::setw(11) << " Time step " << (std::to_string(timestep) + "/" + std::to_string(numSteps) + ".  ") << "Flow time = " << std::setw(8) << std::to_string(timestep*dt) + ".  " << particles.size() << " particles" << std::endl;

        convectParticles();

        // Collapse second row of buffer wake
        collapseBufferWake();

//        moveGeometry();

        // Advance first row of wake panels to second row
        convectBufferWake();

        // Influence of body and particles onto particles
        particleStrengthUpdate();

        // Build/update particle octree
        if(accelerate){
            partOctree.removeData();
            partOctree.setMaxMembers(10); //Barnes Hut
            partOctree.addData(particles);

            FMM.build(&partOctree);
        }

        // Influence of particles onto body
        setSourceStrengthsVP();

        // Solve system of equations
        matrix_convergence = solveVPmatrixEq();

        // Check for simulation convergence
        if(!matrix_convergence){
            std::cout << "*** Warning : Matrix solution did not converge ***" << std::endl;
            std::exit(0);
        }else{
            bool soln_conv = solutionConvergence();
            if (soln_conv) {
                if(printFlag) std::cout << "Converged..." << std::endl;
                break;
            }
        }

        if(unsteady) compVelocity();

        writeFilesVP();

    }

    if(printFlag) std::cout << "\n (\u2713 - Complete, X - Not Requested)\n"<< std::setw(17) << std::left << " Solve System" << std::setw(17) << std::left << "Surface Data" << std::setw(18) << std::left << "Trefftz Plane" <<  std::setw(14) << std::left << "Streamlines" << std::setw(24) << std::left << "Stability Derivatives" << std::setw(23) << std::left << "Volume Mesh" << std::endl << std::setw(19) << std::left << " " + check << std::flush;

    if(!unsteady) compVelocity(); // calculated each time step in unsteady
    if(printFlag) std::cout << std::setw(19) << std::left << check << std::flush;

    trefftzPlaneAnalysisVP();
    if(printFlag) std::cout << std::setw(20) << std::left << check << std::flush;


    if (surfStreamFlag)
    {
        createStreamlines();
        if(printFlag) std::cout << std::setw(14) << std::left << check << std::flush;
    }
    else
    {
        if(printFlag) std::cout << std::setw(14) << std::left << "X" << std::flush;
    }

    if (stabDerivFlag)
    {
        stabilityDerivativesVP();
        if(printFlag) std::cout << std::setw(24) << std::left << check << std::flush;
    }
    else
    {
        if(printFlag) std::cout << std::setw(24) << std::left << "X" << std::flush;
    }


    if (params->volMeshFlag) {
        createVolMesh();
        populateVolMesh();
        if (printFlag) std::cout << std::setw(23) << std::left << check << std::endl;
    }
    else
    {
        if (printFlag) std::cout << std::setw(23) << std::left << "X" << std::endl;
    }


    if (!matrix_convergence && printFlag)
    {
        std::cout << "*** Warning : Solution did not converge ***" << std::endl;
    }

    if (printFlag) writeFilesVP();


}

void cpCaseVP::convectBufferWake(){
    wake2Doublets.resize(static_cast<Eigen::VectorXd::Index>(w2panels->size()));
    for (wakePanels_index_type i=0; i<w2panels->size(); i++)
    {
        (*w2panels)[i]->setPrevStrength((*w2panels)[i]->getMu());
        double parentMu = (*w2panels)[i]->getBufferParent()->getMu();
        (*w2panels)[i]->setMu(parentMu);
        wake2Doublets[static_cast<Eigen::VectorXd::Index>(i)] = parentMu;
    }

}


void cpCaseVP::setSourceStrengthsVP(){
    sigmas.resize(static_cast<Eigen::VectorXd::Index>(bPanels->size()));

    for (bodyPanels_index_type i=0; i<bPanels->size(); i++)
    {
        Eigen::Vector3d sumVelInfl = Eigen::Vector3d::Zero();

        // Influence from particles
        if (accelerate){
            sumVelInfl += FMM.barnesHutVel((*bPanels)[i]->getCenter());
        }
        else
        {
            for(particles_index_type j=0; j<particles.size(); j++)
            {
                sumVelInfl += particles[j]->velInfl((*bPanels)[i]->getCenter());
            }
        }

        // Influence from vortex filaments
        for(filaments_index_type j=0; j<filaments.size(); j++)
        {
            sumVelInfl += filaments[j]->velInfl( (*bPanels)[i]->getCenter() );
        }

        (*bPanels)[i]->setSigma( Vinf((*bPanels)[i]->getCenter()) + sumVelInfl , 0 );

        sigmas(static_cast<Eigen::VectorXd::Index>(i)) = (*bPanels)[i]->getSigma();
    }
}

bool cpCaseVP::solutionConvergence(){

    if(unsteady | (timestep < 5)){
        return false;
    }

    // Check for simulation convergence. Will use the Trefftz plane because it is so cheap. Medium size mesh takes <1% of time step time.

    // Store previous time step values
    double prevCD = CD_trefftz;
    double prevCL = CL_trefftz;

    trefftzPlaneAnalysisVP();

    double changeCD = std::abs((CD_trefftz - prevCD)/prevCD);
    double changeCL = std::abs((CL_trefftz - prevCL)/prevCL);

    double val = 0.0001; // 0.01%

    if((changeCL < val) && (changeCD < val)){
        return true;
    }

    return false;

}



bool cpCaseVP::solveMatrixEq(){
    bool converged = true;

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

    for (bodyPanels_index_type i=0; i<bPanels->size(); i++)
    {
        (*bPanels)[i]->setMu(doubletStrengths(static_cast<Eigen::VectorXd::Index>(i)));
        (*bPanels)[i]->setPotential(Vinf((*bPanels)[i]->getCenter()));
    }

    for (wakePanels_index_type i=0; i<wPanels->size(); i++)
    {
        (*wPanels)[i]->setMu();
        (*wPanels)[i]->setPotential(Vinf((*wPanels)[i]->getCenter()));
    }

    return converged;
}

bool cpCaseVP::solveVPmatrixEq(){
    bool converged = true;

    Eigen::MatrixXd* A = geom->getA();
    Eigen::MatrixXd* B = geom->getB();
    Eigen::MatrixXd* C = geom->getC();
    Eigen::VectorXd RHS = -(*B)*(sigmas) - (*C)*(wake2Doublets); // For first timestep, doublet strength for second row of buffer wake is zero, which is equivalent to not including their influence.
    Eigen::VectorXd doubletStrengths(bPanels->size());


    Eigen::BiCGSTAB<Eigen::MatrixXd> res;
    res.compute((*A));
    doubletStrengths = res.solve(RHS);
    if (res.error() > pow(10,-10))
    {
        converged = false;
    }

    for (bodyPanels_index_type i=0; i<bPanels->size(); i++)
    {
        (*bPanels)[i]->setMu(doubletStrengths(static_cast<Eigen::VectorXd::Index>(i)));
        (*bPanels)[i]->setPotential(VinfPlusVecPot((*bPanels)[i]->getCenter()));

    }

    for (wakePanels_index_type i=0; i<wPanels->size(); i++)
    {
        (*wPanels)[i]->setPrevStrength((*wPanels)[i]->getMu());
        (*wPanels)[i]->setMu();

        (*wPanels)[i]->setPotential(VinfPlusVecPot((*wPanels)[i]->getCenter()));

        (*w2panels)[i]->setPotential(VinfPlusVecPot((*w2panels)[i]->getCenter())); // Included in this loop because there are always the same number of w2pans as w1pans
    }

    return converged;
}


bool cpCaseVP::edgeIsUsed(edge* thisEdge, std::vector<edge*> pEdges){

    for(size_t i=0; i<pEdges.size(); i++){
        if(thisEdge == pEdges[i]){
            return true;
        }
    }

    return false;
}

void cpCaseVP::collapseBufferWake(){

    std::vector<edge*> usedEdges;

    for(wakePanels_index_type i=0; i<w2panels->size(); i++)
    {
        std::vector<edge*> pEdges = (*w2panels)[i]->edgesInOrder();

        Eigen::Vector3d strength = Eigen::Vector3d::Zero();
        for (size_t j=1; j<4; j++) //edges 2-4 only
        {
            if (!edgeIsUsed(pEdges[j],usedEdges))
            {
                usedEdges.push_back(pEdges[j]);
                strength += (*w2panels)[i]->edgeStrength( pEdges[j], static_cast<int>(j) ); // Don't need to pass in pEdges...
            }
        }

        Eigen::Vector3d pos = rungeKuttaStepper((*w2panels)[i]->getCenter());
//        Eigen::Vector3d pos = (*w2panels)[i]->getCenter(); // For moving geometry
        Eigen::Vector3d ptVel = Vinf(pos);
        double radius = (*w2panels)[i]->getPartRadius(ptVel,dt); // VinfLocal

        particle* p = new particle(pos, strength, radius, {0,0,0}, {0,0,0}, timestep); // Zeroes are previous pos and strength values used for Adams-Bashforth scheme
        p->parentPanel = (*w2panels)[i];
        particles.push_back(p);

    }

    // Create filament
    if(filaments.size() == 0)
    {
        for(wakePanels_index_type i=0; i<w2panels->size(); i++)
        {
            vortexFil* fil;
            Eigen::Vector3d p1,p2;

            p1 = (*w2panels)[i]->pointsInOrder()[2]->getPnt();
            p2 = (*w2panels)[i]->pointsInOrder()[3]->getPnt();

            fil = new vortexFil(p1, p2,-(*w2panels)[i]->getMu(), (*w2panels)[i]); // Negative strength is because filament is actually the upstream edge being convected which is oriented the opposite direction as downstream edge

            filaments.push_back(fil);
            (*w2panels)[i]->setVortFil(fil);
        }
    }
    else
    {
        for(wakePanels_index_type i=0; i<w2panels->size(); i++){
            filaments[i]->setStrength(-(*w2panels)[i]->getMu()); // Negative strength is because filament is actually the upstream edge being convected which is oriented the opposite direction as downstream edge
        }

    }
}

void cpCaseVP::compVelocity(){
    //  Velocity Survey with known doublet and source strengths
    CM.setZero();
    Eigen::Vector3d moment;
    Fbody = Eigen::Vector3d::Zero();
    bodyPanel* p;

    for (bodyPanels_index_type i=0; i<bPanels->size(); i++)
    {
        p = (*bPanels)[i];
        p->computeVelocity(PG,Vinf(p->getCenter()));
        if (unsteady) {
            p->computeCp( Vinf((*bPanels)[i]->getCenter()).norm() , dt ); // Katz 13.169 ///
        }else{
            p->computeCp( Vmag );
        }
        Fbody += -p->getCp()*p->getArea()*p->getBezNormal()/params->Sref;
        moment = p->computeMoments(params->cg);
        CM(0) += moment(0)/(params->Sref*params->bref);
        CM(1) += moment(1)/(params->Sref*params->cref);
        CM(2) += moment(2)/(params->Sref*params->bref);
    }

    Fwind = bodyToWind(Fbody);

    if(unsteady){
        //        trefftzPlaneAnalysisVP(); // Trefftz doesn't mean anything right now.

        //   flow_time | CN | CA | CY | CL | CD | CY | Cm | Cl | Cn |
        solnMat(timestep-1,0) = (timestep-1) * dt;

        solnMat(timestep-1,1) = Fbody.z();
        solnMat(timestep-1,2) = Fbody.x();
        solnMat(timestep-1,3) = Fbody.y();

        solnMat(timestep-1,4) = Fwind.z(); // Need to modify for unsteady sims
        solnMat(timestep-1,5) = Fwind.x();
        solnMat(timestep-1,6) = Fwind.y();

        solnMat(timestep-1,7) = CM.y();
        solnMat(timestep-1,8) = CM.x();
        solnMat(timestep-1,9) = CM.z();
    }


}


void cpCaseVP::trefftzPlaneAnalysisVP(){

    std::vector<wake*> wakes = geom->getWakes();
    CL_trefftz = 0;
    CD_trefftz = 0;
    for (size_t i=0; i<wakes.size(); i++)
    {
        wakes[i]->trefftzPlaneVP(Vmag,params->Sref, &particles, timestep);
        CL_trefftz += wakes[i]->getCL()/PG;
        CD_trefftz += wakes[i]->getCD()/pow(PG,2);
    }
}

void cpCaseVP::stabilityDerivativesVP()
{
    double delta = 0.5;
    double dRad = delta*M_PI/180;
    cpCaseVP dA(geom,Vmag,alpha+delta,beta,mach,params);
    cpCaseVP dB(geom,Vmag,alpha,beta+delta,mach,params);

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


void cpCaseVP::particleStrengthUpdate(){
    // This function uses the update equations found in 'Vortex Methods for DNS of...' by Plouhmhans. It uses the particle strength exchange for the viscous diffusion

    std::vector<Eigen::Vector3d> stretchDiffVec; // Creating vector values because the strength change needs to be set after all particle influences have been calculated

    for(particles_index_type i=0; i<particles.size(); i++)
    {
        Eigen::Vector3d dAlpha_diff = Eigen::Vector3d::Zero();
        Eigen::Vector3d dAlpha_stretch = Eigen::Vector3d::Zero();

        // FarField2 condition
        if( std::abs( particles[i]->pos.x() ) > 20*params->cref)
        {
            stretchDiffVec.push_back(Eigen::Vector3d::Zero());
        }
        // FarField1 condition
        else if( std::abs( particles[i]->pos.x() ) > 12*params->cref)
        {
            if(accelerate)
            {
                dAlpha_stretch = FMM.barnesHutStretch(particles[i]); ///
                dAlpha_diff = FMM.barnesHutDiff(particles[i]);
            }
            else
            {
                // Stretching from particles
                for(particles_index_type j=0; j<particles.size(); j++)
                {
                    if(i!=j) // Kroneger Delta Function
                    {
                        dAlpha_diff += particles[i]->viscousDiffusion(particles[j]);
                        dAlpha_stretch += particles[i]->vortexStretching(particles[j]);
                    }
                }
            }
            stretchDiffVec.push_back( dAlpha_diff + dAlpha_stretch );
        }
        else // Treat normally
        {
            // Influence from Particles
            if(accelerate && timestep > 3){
                dAlpha_stretch += FMM.barnesHutStretch(particles[i]);
                dAlpha_diff += FMM.barnesHutDiff(particles[i]);
            }else
            {
                for(particles_index_type j=0; j<particles.size(); j++)
                {
                    if(i!=j){ // Kroneger Delta Function
                        dAlpha_diff += particles[i]->viscousDiffusion(particles[j]);
                        dAlpha_stretch += particles[i]->vortexStretching(particles[j]);
                    }
                }
            }

            // Stretching from body panels
            for(bodyPanels_index_type j=0; j<bPanels->size(); j++)
            {
                dAlpha_diff += (*bPanels)[j]->partStretching(particles[i]);
            }

            // Stretching from wake panels
            for(wakePanels_index_type j=0;j<wPanels->size(); j++)
            {
                dAlpha_diff += (*wPanels)[j]->partStretching(particles[i]);
            }

            stretchDiffVec.push_back( dAlpha_diff + dAlpha_stretch );
        }
    }
    for(particles_index_type i=0;i<particles.size();i++)
    {
        Eigen::Vector3d newStrength;
        if(particles[i]->getprevStrengthUpdate().isZero())
        {
            newStrength = particles[i]->strength + stretchDiffVec[i]*dt;
        }else
        {   // Adams bashforth
            newStrength = particles[i]->strength + dt*(1.5*stretchDiffVec[i] - 0.5*particles[i]->getprevStrengthUpdate());
        }
        particles[i]->setprevStrengthUpdate(stretchDiffVec[i]);
        particles[i]->setStrength(newStrength);
    }
}

void cpCaseVP::convectParticles(){
    std::vector<Eigen::Vector3d> newPartPositions;

    for(particles_index_type i=0;i<particles.size();i++){

        Eigen::Vector3d newPos;

        // Adams bashforth scheme
        Eigen::Vector3d velOnPart = velocityInflFromEverything(particles[i]);

        if(particles[i]->getPrevVelInfl().isZero())
        {
            newPos = particles[i]->pos + dt*velOnPart;
        }
        else{
            newPos = particles[i]->pos + dt*(1.5*velOnPart - 0.5*particles[i]->getPrevVelInfl());
        }

        particles[i]->setPrevVelInfl(velOnPart);
        newPartPositions.push_back(newPos);
        //        }
    }

    for(particles_index_type i=0;i<particles.size();i++){
        particles[i]->setPos(newPartPositions[i]);

    }
}


Eigen::Vector3d cpCaseVP::Vinf(Eigen::Vector3d POI){

    if (!unsteady)
    {
        return windToBody( Vmag , alpha , beta );
    }
    else
    {
        Eigen::Vector3d localVel;

        // U = U3 + (-q*z + r*y)
        localVel.x() = bodyKin(timestep-1, 0) - bodyKin(timestep-1, 4)*POI.z() + bodyKin(timestep-1, 5)*POI.y();

        // V = V3 + (-r*x + p*z)
        localVel.y() = bodyKin(timestep-1, 1) - bodyKin(timestep-1, 5)*POI.x() + bodyKin(timestep-1, 3)*POI.z();

        // W = W3 + (-p*y + q*x)
        localVel.z() = bodyKin(timestep-1, 2) - bodyKin(timestep-1, 3)*POI.y() + bodyKin(timestep-1, 4)*POI.x();

        return localVel;
    }
}

Eigen::Vector3d cpCaseVP::VinfPlusVecPot(Eigen::Vector3d POI){
    Eigen::Vector3d vInfluence = Vinf(POI);

    // Particles
    if(accelerate && timestep > 2)
    {
        vInfluence += FMM.barnesHutVel(POI);
    }else
    {
        for (particles_index_type i=0; i<particles.size(); i++) {
            vInfluence += particles[i]->velInfl(POI);
        }
    }


    // Filaments
    for (filaments_index_type i=0; i<filaments.size(); i++)
    {
        vInfluence += filaments[i]->velInfl(POI);
    }

    return vInfluence;
}



Eigen::Vector3d cpCaseVP::rungeKuttaStepper( Eigen::Vector3d POI ){

    // RK4 algorithm makes symmeterized influences VERY impractical

    Eigen::Vector3d k1, k2, k3, k4;

    k1 = velocityInflFromEverything(POI);
    k2 = velocityInflFromEverything(POI+k1*dt/2);
    k3 = velocityInflFromEverything(POI+k2*dt/2);
    k4 = velocityInflFromEverything(POI+k3*dt);

    return POI + dt*(k1/6 + k2/3 + k3/3 + k4/6);

}

Eigen::Vector3d cpCaseVP::velocityInflFromEverything( Eigen::Vector3d POI ){
    // Freestream influence
    Eigen::Vector3d velOnPart = Vinf(POI);

    // Particle influence
    if(accelerate && timestep > 3){
        velOnPart += FMM.barnesHutVel(POI);
    }
    else{
        for(particles_index_type j=0;j<particles.size();j++)
        {
            velOnPart += particles[j]->velInfl(POI);
        }
    }


    // Body panel influence
    for(bodyPanels_index_type j=0;j<bPanels->size();j++){
        velOnPart += (*bPanels)[j]->panelV(POI);
    }


    // Buffer wake influence
    for(wakePanels_index_type j=0;j<wPanels->size();j++){
        velOnPart += (*wPanels)[j]->panelV(POI);
    }

    for(wakePanels_index_type j=0;j<w2panels->size();j++){
        velOnPart += (*w2panels)[j]->panelV(POI);
    }

    // Vortex Filament influence
    for(filaments_index_type i=0; i<filaments.size(); i++){
        velOnPart +=filaments[i]->velInfl(POI);
    }

    return velOnPart;
}

Eigen::Vector3d cpCaseVP::velocityInflFromEverything(particle* part){

    // Function is overloaded for symmeterized velocity influence for particles
    Eigen::Vector3d pos = part->pos;

    // Freestream influence
    Eigen::Vector3d velOnPart = Vinf(pos);

    // If particle reaches very far field condition
    if ((pos-Eigen::Vector3d::Zero()).norm() > 20*params->cref)
    {
        return velOnPart;
    }

    // Particle influence
    if(accelerate && timestep > 2){
        velOnPart += FMM.barnesHutVel(part);
    }
    else{
        for(particles_index_type j=0;j<particles.size();j++)
        {
            velOnPart += particles[j]->velInfl(pos);
        }
    }

    // If particle reaches first farfield condition
    if((pos-Eigen::Vector3d::Zero()).norm() > 12*params->cref)
    {
        return velOnPart;
    }


    // Body panel influence
    for(bodyPanels_index_type j=0;j<bPanels->size();j++){
        velOnPart += (*bPanels)[j]->panelV(pos);
    }


    // Buffer wake influence
    for(wakePanels_index_type j=0;j<wPanels->size();j++){
        velOnPart += (*wPanels)[j]->panelV(pos);
    }
    for(wakePanels_index_type j=0;j<w2panels->size();j++){
        velOnPart += (*w2panels)[j]->panelV(pos);
    }


    // Vortex Filament influence
    for(filaments_index_type i=0; i<filaments.size(); i++){
        velOnPart +=filaments[i]->velInfl(pos);
    }

    return velOnPart;
}

void cpCaseVP::writeFilesVP(){
    std::stringstream caseLabel;
    caseLabel << "V" << Vmag << "_Mach" << mach << "_alpha" << alpha << "_beta" << beta;
    boost::filesystem::path subdir = params->inputFile->path + caseLabel.str();
    if (!boost::filesystem::exists(subdir))
    {
        boost::filesystem::create_directories(subdir);
    }
    Eigen::MatrixXd nodeMat = geom->getNodePnts();

    writeBodyDataVP(subdir.string(),nodeMat);
    writeWakeDataVP(subdir.string(),nodeMat);
    writeBuffWake2Data(subdir.string(),nodeMat);
    writeSpanwiseData(subdir.string());
    writeParticleData(subdir.string());
    writeFilamentData(subdir.string());

    if (params->volMeshFlag && cells.size() > 0) //
    {
        writeVolMeshData(subdir.string(), pts, cells);
    }

    if (params->surfStreamFlag)
    {
        writeBodyStreamlines(subdir.string());
    }

}

void cpCaseVP::writeBodyDataVP(const std::string & path,const Eigen::MatrixXd &nodeMat){
    std::vector<cellDataArray> data;
    cellDataArray mu("Doublet Strengths"),sigma("Source Strengths"),pot("Velocity Potential"),V("Velocity"),Cp("Cp"),bN("bezNormals");
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> con(bPanels->size(),3);
    mu.data.resize(static_cast<Eigen::MatrixXd::Index>(bPanels->size()),1);
    sigma.data.resize(static_cast<Eigen::MatrixXd::Index>(bPanels->size()),1);
    pot.data.resize(static_cast<Eigen::MatrixXd::Index>(bPanels->size()),1);
    V.data.resize(static_cast<Eigen::MatrixXd::Index>(bPanels->size()),3);
    Cp.data.resize(static_cast<Eigen::MatrixXd::Index>(bPanels->size()),1);
    bN.data.resize(static_cast<Eigen::MatrixXd::Index>(bPanels->size()),3);
    for (bodyPanels_index_type i=0; i<bPanels->size(); i++)
    {
    	Eigen::MatrixXd::Index ii(static_cast<Eigen::MatrixXd::Index>(i));
        mu.data(ii,0) = (*bPanels)[i]->getMu();
        sigma.data(ii,0) = (*bPanels)[i]->getSigma();
        pot.data(ii,0) = (*bPanels)[i]->getPotential();
        V.data.row(ii) = (*bPanels)[i]->getGlobalV();
        Cp.data(ii,0) = (*bPanels)[i]->getCp();
        con.row(ii) = (*bPanels)[i]->getVerts();
        bN.data.row(ii) = (*bPanels)[i]->getBezNormal();
    }

    data.push_back(mu);
    data.push_back(sigma);
    data.push_back(pot);
    data.push_back(V);
    data.push_back(Cp);
    data.push_back(bN);

    piece body;
    body.pnts = nodeMat;
    body.connectivity = con;
    body.cellData = data;

    std::string fname = path+"/surfaceData-"+std::to_string(timestep)+".vtu";
    VTUfile bodyFile(fname,body);
}

void cpCaseVP::writeWakeDataVP(const std::string & path, const Eigen::MatrixXd &nodeMat){
    std::vector<cellDataArray> data;
    cellDataArray mu("Doublet Strengths"),pot("Velocity Potential");
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> con(wPanels->size(),(*wPanels)[0]->getVerts().size()); // Assumes wake won't mix tris and quads
    mu.data.resize(static_cast<Eigen::MatrixXd::Index>(wPanels->size()),1);
    pot.data.resize(static_cast<Eigen::MatrixXd::Index>(wPanels->size()),1);
    for (wakePanels_index_type i=0; i<wPanels->size(); i++)
    {
    	Eigen::MatrixXd::Index ii(static_cast<Eigen::MatrixXd::Index>(i));
        mu.data(ii,0) = (*wPanels)[i]->getMu();
        pot.data(ii,0) = (*wPanels)[i]->getPotential();
        con.row(ii) = (*wPanels)[i]->getVerts();
    }

    data.push_back(mu);
    data.push_back(pot);

    piece wake;
    wake.pnts = nodeMat;
    wake.connectivity = con;
    wake.cellData = data;

    std::string fname = path+"/wakeData-"+std::to_string(timestep)+".vtu";
    VTUfile wakeFile(fname,wake);
}


void cpCaseVP::writeBuffWake2Data(const std::string & path, const Eigen::MatrixXd &nodeMat){
    std::vector<cellDataArray> data;
    cellDataArray mu("Doublet Strengths"),pot("Velocity Potential");
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> con;
    con.resize(static_cast<Eigen::MatrixXd::Index>(w2panels->size()),4);
    mu.data.resize(static_cast<Eigen::MatrixXd::Index>(w2panels->size()),1);
    pot.data.resize(static_cast<Eigen::MatrixXd::Index>(w2panels->size()),1);
    for (wakePanels_index_type i=0; i<w2panels->size(); i++)
    {
    	Eigen::MatrixXd::Index ii(static_cast<Eigen::MatrixXd::Index>(i));
        mu.data(ii,0) = (*w2panels)[i]->getMu();
        pot.data(ii,0) = (*w2panels)[i]->getPotential();
        con.row(ii) = (*w2panels)[i]->getVerts();
    }
    data.push_back(mu);
    data.push_back(pot);

    piece wake;
    wake.pnts = nodeMat;
    wake.connectivity = con;
    wake.cellData = data;
    std::string fname = path+"/bufferWake2Data-"+std::to_string(timestep)+".vtu";
    VTUfile wakeFile(fname,wake);
}

void cpCaseVP::writeFilamentData(const std::string & path){
    std::vector<cellDataArray> data;
    cellDataArray mu("Gamma");
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> con;

    Eigen::MatrixXd nodeMat(2*filaments.size(),3);
    for(filaments_index_type i=0; i<filaments.size(); i++){
    	Eigen::MatrixXd::Index ii(static_cast<Eigen::MatrixXd::Index>(i));
        nodeMat.row(2*ii) = filaments[i]->getP1();
        nodeMat.row(2*ii+1) = filaments[i]->getP2();
    }


    con.resize(static_cast<Eigen::MatrixXd::Index>(filaments.size()),2); // 2 pnts per filament

    mu.data.resize(static_cast<Eigen::MatrixXd::Index>(filaments.size()),1);

    size_t nodeCounter = 0; // Used to make filament end points
    for (filaments_index_type i=0; i<filaments.size(); i++)
    {
    	Eigen::MatrixXd::Index ii(static_cast<Eigen::MatrixXd::Index>(i));
        mu.data(ii,0) = filaments[i]->getStrength();
        con(ii,0) = nodeCounter;
        con(ii,1) = nodeCounter+1;
        nodeCounter = nodeCounter+2;
    }

    data.push_back(mu);

    piece fils;
    fils.pnts = nodeMat;
    fils.connectivity = con;
    fils.cellData = data;
    std::string fname = path+"/filaments-"+std::to_string(timestep)+".vtu";
    VTUfile filFile(fname,fils);
}

void cpCaseVP::writeParticleData(const std::string & path){

    Eigen::MatrixXd partMat(particles.size(),3);
    for (particles_index_type i=0; i<particles.size(); i++)
    {
        partMat.row(static_cast<Eigen::MatrixXd::Index>(i)) = particles[i]->pos;
    }

    std::vector<cellDataArray> data;
    cellDataArray strength("Strength"), shedTime("Time Step Shed");
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> con(particles.size(),1);
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> shed(particles.size(),1);

    strength.data.resize(static_cast<Eigen::MatrixXd::Index>(particles.size()),3);
    shedTime.data.resize(static_cast<Eigen::MatrixXd::Index>(particles.size()),1);
    for (particles_index_type i=0; i<particles.size(); i++)
    {
    	Eigen::MatrixXd::Index ii(static_cast<Eigen::MatrixXd::Index>(i));
        strength.data.row(ii) = particles[i]->strength;
        shedTime.data(ii,0) = particles[i]->shedTime;
        con(ii) = i;
    }
    data.push_back(strength);
    data.push_back(shedTime);

    piece parts;
    parts.pnts = partMat;
    parts.connectivity = con;
    parts.cellData = data;

    std::string fname = path+"/particleData-"+std::to_string(timestep)+".vtu";
    VTUfile partFile(fname,parts);
}


void cpCaseVP::readBodyKinFile(){
    std::ifstream fid;
    fid.open(params->bodyKinFileLoc);

    bodyKin.resize(numSteps, 6); // U, V, W, p, q, r

    for (int i=0; i<numSteps; i++)
    {
        for(int j=0; j < 6; j++)
        {
            fid >> bodyKin(i,j);
        }
    }

    fid.close();
}


void cpCaseVP::populateVolMesh(){

    // Clear Past Timestep Mesh
    volMeshDat.velocity.clear();
    volMeshDat.coef_press.clear();

    for (cells_index_type i=0; i<cells.size(); i++) {
        Eigen::Vector3d velInCell = velocityInflFromEverything(volMeshDat.cellCenter[i]);
        volMeshDat.velocity.push_back(velInCell);
    }

    for (cells_index_type i=0; i<cells.size(); i++) {
        // Incompressible Bernoulli equation
        double Cp = pow( volMeshDat.velocity[i].norm() / Vmag , 2 );
        volMeshDat.coef_press.push_back(Cp);
    }

}




//void cpCaseVP::moveGeometry(){
//    // Steps to continue development include (but are not limited to)
//        // In the 'velocityInflFromEverything' functions, initialize with zero velocity instead of freestream
//        // Uncomment all 'moveGeometry' calls
//        // The 'collapseBufferWake' function must seed particle at the panel center instead of tracking it with freestream.
//        // Instead of the 'bodyMovement' vector, you'll probably want to use Eigen::MatrixXd and just slice with .row
//
//    std::vector<double> bodyMovement = { -99.619, 0, -8.7156, 0, 0, 0 };
//
//    geom->moveGeom(bodyMovement);
//
//    // Filaments aren't part of geom... but they should be.
//    for (int i=0; i<filaments.size(); i++) {
//        filaments[i]->moveFilament(bodyMovement, dt);
//    }
//
//
//}






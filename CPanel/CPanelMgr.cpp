/*******************************************************************************
 * Copyright (c) 2015 Chris Satterwhite
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

#include <iostream>
#include <iomanip>

#include "CPanelMgr.h"

void caseMgr::setCases()
{
    cpCase* c;
    cpCaseVP* cVP;
    for (int v=0; v<p->velocities.rows(); v++)
    {
        for (int a=0; a<p->alphas.rows(); a++)
        {
            for (int b=0; b<p->betas.rows(); b++)
            {
                for (int m=0; m<p->machs.rows(); m++)
                {
                    if (p->vortexParticles) {
                        cVP = new cpCaseVP(geom,p->velocities(v),p->alphas(a),p->betas(b),p->machs(m),p);
                        casesVP.push_back(cVP);
                    }
                    else{
                        c = new cpCase(geom,p->velocities(v),p->alphas(a),p->betas(b),p->machs(m),p);
                        cases.push_back(c);
                    }
                }
            }
        }
    }
}

void caseMgr::runCases()
{

    if( p->vortexParticles)
    {
        std::cout << "\nRunning " << casesVP.size() << " Cases...\n" << std::flush;

        for (cases_index_type i=0; i<casesVP.size(); i++)
        {
            std::string out;
            std::stringstream outstream;
            outstream << "\n\nCase # " << i+1 << "/" << casesVP.size() << ":" << std::endl;
            out = outstream.str();
            std::cout << out << std::flush;

            casesVP[i]->run(true,p->surfStreamFlag,p->stabDerivFlag);
        }
    }
    else
    {
        std::cout << "\nRunning " << cases.size() << " Cases... (\u2713 - Complete, X - Not Requested)\n" << std::endl;
        std::cout << std::setw(10) << std::left << "Case #" << std::setw(15) << std::left << "Solve System" << std::setw(15) << std::left << "Surface Data" << std::setw(16) << std::left << "Trefftz Plane" <<  std::setw(14) << std::left << "Streamlines" << std::setw(24) << std::left << "Stability Derivatives" << std::setw(15) << std::left << "Volume Mesh" << std::endl;

        for (cases_index_type i=0; i<cases.size(); i++)
        {
            std::string out;
            std::stringstream outstream;
            outstream << i+1 << "/" << cases.size();
            out = outstream.str();
            std::cout << std::setw(10) << std::left << out << std::flush;

            cases[i]->run(true,p->surfStreamFlag,p->stabDerivFlag);
        }
    }



    std::cout << "Complete" << std::endl;

}

void caseMgr::writeSummary()
{
    std::ofstream out;
    std::string outFile = p->inputFile->path+p->inputFile->name+".CPout";
    out.open(outFile);

    if (out)
    {
        out << "-----------INPUTS-----------\n\n";
        p->print(out);
        out << "\n\n";
        out << "-----------RESULTS-----------\n\n";
        out << std::left; // Left justify results
        out << std::fixed;

        outSpacing.resize(9);
        outSpacing << 11,7,8,8,8,8,12,12,12;
        for (cases_index_type i=0; i<cases.size(); i++)
        {
            writeCase(i+1, cases[i], out);
        }
        for (cases_index_type i=0; i<casesVP.size(); i++) {
            writeCase(i+1, casesVP[i], out);
        }
    }
    std::cout << "Case data written to " << outFile << std::endl;

}

void caseMgr::writeCase(size_t caseNumber, cpCase* c, std::ofstream &outStream)
{
    outStream << std::endl;
    outStream << std::setprecision(3);
    outStream << "---Case #" << caseNumber << "---" << std::endl;
    outStream << "\n\t--Flow Conditions--\n" << std::endl;
    outStream << "\t\t" << std::setw(outSpacing(0)) << "Velocity" << std::setw(outSpacing(1)) << "Mach" << std::setw(outSpacing(2)) << "Alpha" << std::setw(outSpacing(3)) << "Beta" << "\n";
    outStream << "\t\t" << std::setw(outSpacing(0)) << c->getV() << std::setw(outSpacing(1)) << c->getMach() << std::setw(outSpacing(2)) << c->getAlpha() << std::setw(outSpacing(3)) << c->getBeta() << "\n";

    outStream << std::setprecision(5);
    outStream << "\n\t--Force Coefficients--\n" << std::endl;
    outStream << "\t\t-Trefftz Plane-" << std::endl;
    outStream << "\t\t\tCL = " << c->getCL() << "\tCDi = " << c->getCD() << std::endl;
    outStream << "\n\t\t-Surface Integrated Force Coefficients-" << std::endl;
    outStream << "\t\t\t" << std::setw(15) << "Body Axis" << std::setw(10) << "CN" << std::setw(10) << "CA" << std::setw(10) << "CY" << std::endl;
    outStream << "\t\t\t" << std::setw(15) << " " << std::setw(10) << c->getBodyForces()(2) << std::setw(10) << c->getBodyForces()(0) << std::setw(10) << c->getBodyForces()(1) << std::endl;
    outStream << "\t\t\t" << std::setw(15) << "Wind Axis" << std::setw(10) << "CL" << std::setw(10) << "CD" << std::setw(10) << "CY" << std::endl;
    outStream << "\t\t\t" << std::setw(15) << " " << std::setw(10) << c->getWindForces()(2) << std::setw(10) << c->getWindForces()(0) << std::setw(10) << c->getWindForces()(1) << std::endl;

    outStream << "\n\t--Moment Coefficients--\n" << std::endl;
    outStream << "\t\t" << std::setw(12) << "Cm (pitch)" << std::setw(12) << "Cl (roll)" << std::setw(12) << "Cn (yaw)" << std::endl;
    outStream << "\t\t" << std::setw(12) << c->getMoment()(1) << std::setw(12) << c->getMoment()(0) << std::setw(12) << c->getMoment()(2) << std::endl;

    if (p->stabDerivFlag)
    {
        outStream << "\n\t--Stability Derivatives--\n" << std::endl;
        outStream << "\t\t" << std::setw(12) << "CL_alpha" << std::setw(12) << "CY_alpha" << std::setw(12) << "Cm_alpha" << std::setw(12) << "Cl_alpha" << std::setw(12) << "Cn_alpha" << std::endl;
        outStream << "\t\t" << std::setw(12) << c->get_dF_dAlpha()(2) << std::setw(12) << c->get_dF_dAlpha()(1) << std::setw(12) << c->get_dM_dAlpha()(1) << std::setw(12) << c->get_dM_dAlpha()(0) << std::setw(12) << c->get_dM_dAlpha()(2) << std::endl;

        outStream << "\t\t" << std::setw(12) << "CL_beta" << std::setw(12) << "CY_beta" << std::setw(12) << "Cm_beta" << std::setw(12) << "Cl_beta" << std::setw(12) << "Cn_beta" << std::endl;
        outStream << "\t\t" << std::setw(12) << c->get_dF_dBeta()(2) << std::setw(12) << c->get_dF_dBeta()(1) << std::setw(12) << c->get_dM_dBeta()(1) << std::setw(12) << c->get_dM_dBeta()(0) << std::setw(12) << c->get_dM_dBeta()(2) << std::endl;
        outStream << "\n\n" << std::endl;
    }

    if (p->unsteady) {
        outStream << "\n\t--Time Dependent Forces and Moments--\n" << std::endl;

        outStream << "\t\t" << std::setw(9) << "Time" << std::setw(10) << "CN_body" << std::setw(10) << "CA_body" << std::setw(10) << "CY_body" <<  std::setw(10) << "CL_wind" << std::setw(10) << "CD_wind" << std::setw(10) << "CY_wind" << std::setw(12) << "Cm (pitch)" << std::setw(12) << "Cl (roll)" << std::setw(12) << "Cn (yaw)" << std::endl;


        Eigen::MatrixXd resMat = c->get_soln_mat();
        for (int i=0; i<resMat.rows(); i++) {
            outStream << "\t\t" << std::setw(9) << resMat(i,0) << std::setw(10) << resMat(i,1) << std::setw(10) << resMat(i,2) << std::setw(10) << resMat(i,3) << std::setw(10) << resMat(i,4) << std::setw(10) << resMat(i,5) <<  std::setw(10) << resMat(i,6) << std::setw(12) << resMat(i,7) << std::setw(12) << resMat(i,8) << std::setw(12) << resMat(i,9) << std::endl;
        }

        //  timestep | CN | CA | CY | CL | CD | CY | Cm | Cl | Cn |
    }

}

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
 *    Connor Sousa - Vortex particle implementation
 *    David D. Marshall - misc. changes
 ******************************************************************************/

#ifndef __CPanel__edge__
#define __CPanel__edge__

#include <stdio.h>
#include <Eigen/Dense>
#include <vector>

class bodyPanel;
class wakePanel;
class cpNode;
class secondBufferWake;

class edge
{
	using bodyPanels_type = std::vector<bodyPanel *>;
	using bodyPanels_index_type = bodyPanels_type::size_type;
	using wakePanels_type = std::vector<wakePanel *>;
	using wakePanels_index_type = wakePanels_type::size_type;
	using secondBufferWakes_type = std::vector<secondBufferWake *>;
	using secondBufferWakes_index_type = secondBufferWakes_type::size_type;

    cpNode* n1;
    cpNode* n2;
    bodyPanels_type bodyPans;
    wakePanels_type wakePans;
    secondBufferWakes_type SBW2;
    bool TE; //Edge is at surface-wake junction
    
    void checkTE();
    
public:
    edge(cpNode* n1,cpNode* n2);
    
    void addBodyPan(bodyPanel* b);
    void addWakePan(wakePanel* w);
    void setNeighbors();
    bool sameEdge(cpNode* node1, cpNode* node2);
    void flipDir();
    
    edge* nextTE();
    
    bool isTE();
    double length();
    double distToEdge(const Eigen::Vector3d &pnt);
    Eigen::Vector3d edgeVelocity();
    
    Eigen::Vector3d TEgamma();
    
    std::vector<cpNode*> getNodes();
    cpNode* getN1() {return n1;}
    cpNode* getN2() {return n2;}
    cpNode* getOtherNode(cpNode* current);

    
    std::vector<bodyPanel*> getBodyPans() {return bodyPans;}
    std::vector<wakePanel*> getWakePans() {return wakePans;}
    bodyPanel* getOtherBodyPan(bodyPanel* currentPan);
    wakePanel* getOtherWakePan(wakePanel* currentPan);
    Eigen::Vector3d getVector();
    Eigen::Vector3d getMidPoint();
    Eigen::Vector3d getNormal();
    
    bool containsNode(cpNode* node);
};

#endif /* defined(__CPanel__edge__) */

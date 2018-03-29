//
//  edge.h
//  CPanel
//
//  Created by Chris Satterwhite on 1/25/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

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
    cpNode* n1;
    cpNode* n2;
    std::vector<bodyPanel*> bodyPans;
    std::vector<wakePanel*> wakePans;
    std::vector<secondBufferWake*> SBW2;
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

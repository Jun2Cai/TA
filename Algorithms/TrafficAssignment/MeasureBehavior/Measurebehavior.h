//
// Created by junjun on 30.09.22.
//

#ifndef ROUTINGFRAMEWORK_MEASUREBEHAVIOR_H
#define ROUTINGFRAMEWORK_MEASUREBEHAVIOR_H


#include "Algorithms/TrafficAssignment/Measurement/CompareTravelTime.h"



class Measurebehavior {


public:
    template <typename InputGraph>
    void measures(std::string anaFileName, std::map<std::string, std::vector<int32_t>> odPairPath, InputGraph inputGraph, AlignedVector<int> trafficFlows, int numIterations) {
    ;

    }


};


#endif //ROUTINGFRAMEWORK_MEASUREBEHAVIOR_H

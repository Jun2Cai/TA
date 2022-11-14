//
// Created by junjun on 30.09.22.
//

#ifndef ROUTINGFRAMEWORK_MEASUREBEHAVIOR_H
#define ROUTINGFRAMEWORK_MEASUREBEHAVIOR_H

#include "Algorithms/TrafficAssignment/Measurement/PercentageOfDifferentSat.h"
#include "Algorithms/TrafficAssignment/Measurement/LongestSubPathMeasure.h"


class Measurebehavior {


public:
    template <typename InputGraph>
    void measures(std::string anaFileName,std::map<std::string, std::vector<int32_t>> odPairPath, InputGraph inputGraph, AlignedVector<int> trafficFlows, int numIterations) {
        PercentageOfDifferentSat percentageOfDifferentSat;
        percentageOfDifferentSat.measure(anaFileName, odPairPath, inputGraph, trafficFlows, numIterations);
        LongestSubPathMeasure longestSubPathMeasure;
        longestSubPathMeasure.measure(anaFileName, odPairPath, inputGraph, trafficFlows, numIterations);
    }


};


#endif //ROUTINGFRAMEWORK_MEASUREBEHAVIOR_H

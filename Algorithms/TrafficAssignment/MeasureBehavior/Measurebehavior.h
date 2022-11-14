//
// Created by junjun on 30.09.22.
//

#ifndef ROUTINGFRAMEWORK_MEASUREBEHAVIOR_H
#define ROUTINGFRAMEWORK_MEASUREBEHAVIOR_H

#include "Algorithms/TrafficAssignment/Measurement/PercentageOfDifferentSat.h"
#include "Algorithms/TrafficAssignment/Measurement/LongestSubPathMeasure.h"
#include "Algorithms/TrafficAssignment/Measurement/CompareTravelTime.h"
#include "Algorithms/TrafficAssignment/Measurement/IncreasedFlow.h"


class Measurebehavior {


public:
    template <typename InputGraph>
    void measures(std::string anaFileName, InputGraph inputGraph, AlignedVector<int> trafficFlows, int numIterations) {
/*        PercentageOfDifferentSat percentageOfDifferentSat;
        percentageOfDifferentSat.setTrafficFlows(trafficFlows);
        percentageOfDifferentSat.measure(anaFileName, odPairPath, inputGraph,  numIterations);
        LongestSubPathMeasure longestSubPathMeasure;
        longestSubPathMeasure.setTrafficFlows(trafficFlows);
        longestSubPathMeasure.measure(anaFileName, odPairPath, inputGraph, trafficFlows, numIterations);
        CompareTravelTime compareTravelTime;
        compareTravelTime.measure(anaFileName, odPairPath, inputGraph,  numIterations);*/
        ;
        IncreasedFlow increasedFlow;
        increasedFlow.setTrafficFlows(trafficFlows);
        increasedFlow.measure(anaFileName, inputGraph, numIterations);

    }


};


#endif //ROUTINGFRAMEWORK_MEASUREBEHAVIOR_H

//
// Created by junjun on 30.09.22.
//

#ifndef ROUTINGFRAMEWORK_MEASUREBEHAVIOR_H
#define ROUTINGFRAMEWORK_MEASUREBEHAVIOR_H


#include "Algorithms/TrafficAssignment/Measurement/CompareTravelTime.h"
#include "Algorithms/TrafficAssignment/Measurement/QuantileRatio.h"


template<typename InputGraph>
class MeasureBehavior {

public:

    explicit MeasureBehavior(const InputGraph& inputGraph, const std::vector<ClusteredOriginDestination>& odPairs) : quantileRatio(inputGraph, odPairs) {}

    void measureFirstIteration(const std::vector<std::vector<int32_t>>& odPairPaths, const AlignedVector<int>& trafficFlows) {
        quantileRatio.measureFirstIteration(odPairPaths, trafficFlows);
    }

    void measureLastIterationAndWriteOutput(const std::vector<std::vector<int32_t>>& odPairPaths, const AlignedVector<int>& trafficFlows, const std::string& anaFileName) {
        quantileRatio.measureLastIterationAndWriteOutput(odPairPaths, trafficFlows, anaFileName);
    }

private:

    QuantileRatio<InputGraph> quantileRatio;

};


#endif //ROUTINGFRAMEWORK_MEASUREBEHAVIOR_H

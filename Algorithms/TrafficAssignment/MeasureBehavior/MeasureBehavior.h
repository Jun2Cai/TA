//
// Created by junjun on 30.09.22.
//

#ifndef ROUTINGFRAMEWORK_MEASUREBEHAVIOR_H
#define ROUTINGFRAMEWORK_MEASUREBEHAVIOR_H


#include "Algorithms/TrafficAssignment/Measurement/TravelTime.h"
#include "Algorithms/TrafficAssignment/Measurement/Quantile.h"
#include "Algorithms/TrafficAssignment/Measurement/LongestSubPath.h"
#include "Algorithms/TrafficAssignment/Measurement/QuantileDistribution.h"
#include "Algorithms/TrafficAssignment/Measurement/TravelTime.h"

template<typename InputGraph>
class MeasureBehavior {

public:

    explicit MeasureBehavior(const InputGraph& inputGraph, const std::vector<ClusteredOriginDestination>& odPairs) : quantile(inputGraph, odPairs), longestSubPath(inputGraph, odPairs), quantileDistribution(inputGraph, odPairs), travelTime(inputGraph, odPairs) {}

    void measureFirstIteration(const std::vector<std::vector<int32_t>>& odPairPaths, const AlignedVector<int>& trafficFlows) {
        quantile.measureFirstIteration(odPairPaths, trafficFlows);
        longestSubPath.measureFirstIteration(odPairPaths, trafficFlows);
        quantileDistribution.measureFirstIteration(odPairPaths, trafficFlows);
        travelTime.measureFirstIteration(odPairPaths);
    }

    void measureLastIterationAndWriteOutput(const std::vector<std::vector<int32_t>>& odPairPaths, const AlignedVector<int>& trafficFlows, const std::string& anaFileName) {
        quantile.measureLastIterationAndWriteOutput(odPairPaths, trafficFlows, anaFileName);
        longestSubPath.measureLastIterationAndWriteOutput(odPairPaths, trafficFlows, anaFileName);
        quantileDistribution.measureLastIterationAndWriteOutput(odPairPaths, trafficFlows, anaFileName);
        travelTime.measureLastIterationAndWriteOutput(odPairPaths, anaFileName);
    }

private:

    Quantile<InputGraph> quantile;
    LongestSubPath<InputGraph> longestSubPath;
    QuantileDistribution<InputGraph> quantileDistribution;
    TravelTime<InputGraph> travelTime;

};


#endif //ROUTINGFRAMEWORK_MEASUREBEHAVIOR_H

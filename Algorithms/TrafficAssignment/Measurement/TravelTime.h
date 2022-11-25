//
// Created by junjun on 30.10.22.
//

#ifndef ROUTINGFRAMEWORK_TRAVELTIME_H
#define ROUTINGFRAMEWORK_TRAVELTIME_H




#include "Tools/Simd/AlignedVector.h"
#include <fstream>
#include <vector>
#include <map>
#include <iostream>
#include "DataStructures/Utilities/OriginDestination.h"

template <typename InputGraph>
class TravelTime  {

public :
    explicit TravelTime(const InputGraph &inputGraph, const std::vector<ClusteredOriginDestination> &odPairs)
            : inputGraph(inputGraph), odPairs(odPairs), travelTimeInFstIter(odPairs.size()) {}

    void measureFirstIteration(const std::vector<std::vector<int32_t>> &odPairPaths) {
        for (int i = 0; i < odPairs.size(); ++i) {
            const auto &path = odPairPaths[i];
            travelTimeInFstIter[i] = getTotalTravelTime(path);
        }
    }

    void measureLastIterationAndWriteOutput(const std::vector<std::vector<int32_t>> &odPairPaths, const std::string &anaFileName) {
    auto odAnaFileName = anaFileName + "_travelTime.csv";
    std::ofstream odAnaFile;
        odAnaFile.open(odAnaFileName);
        if (!odAnaFile.good())
            throw std::invalid_argument("file cannot be opened -- '" + odAnaFileName + "'");
        odAnaFile << "origin,destination,,first iteration,,last iteration,,ratio" << std::endl ;
        double sumRatio = 0;
        for (int i = 0; i < odPairs.size(); ++i) {
            const auto &path = odPairPaths[i];
            assert((path.empty() && quanInFstIter[i] == StaticalFunction::Quantiles()) ||
                   (!path.empty() && quanInFstIter[i] != StaticalFunction::Quantiles()));
            odAnaFile << std::to_string(odPairs[i].origin) << "," << std::to_string(odPairs[i].destination) << ",,";
            int travelTimeInFirstIter = travelTimeInFstIter[i];
            int travelTimeInLastIter = getTotalTravelTime(path);
            double ratio = (travelTimeInFirstIter == 0)?1:(static_cast<double>(travelTimeInLastIter)/static_cast<double>(travelTimeInFirstIter));
            odAnaFile << travelTimeInFirstIter << ",," << travelTimeInLastIter << ",," << ratio << std::endl;
            sumRatio += ratio;

        }

        odAnaFile << ",,,,,,mean," << sumRatio/static_cast<double>(odPairs.size());
}


private :
            const InputGraph &inputGraph;
            const std::vector<ClusteredOriginDestination> &odPairs;
            std::vector<int> travelTimeInFstIter;
            int getTotalTravelTime(const std::vector<int32_t> &path) {
                int totalTime = 0;
                for (auto e:path) {
                    totalTime += inputGraph.travelTime(e);
                }
                return totalTime;
            }


};

#endif //ROUTINGFRAMEWORK_TRAVELTIME_H

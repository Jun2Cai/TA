//
// Created by junjun on 25.11.22.
//

#ifndef ROUTINGFRAMEWORK_LONGESTSUBPATH_H
#define ROUTINGFRAMEWORK_LONGESTSUBPATH_H
#include "Algorithms/TrafficAssignment/Measurement/IMeasure.h"
#include "Algorithms/TrafficAssignment/Measurement/StaticalFunction.h"
#include "Tools/Simd/AlignedVector.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include <fstream>
#include <vector>
#include <iostream>

template<typename InputGraph>
class LongestSubPath : public IMeasure {

public :
    struct LongestSubPathOfDifferentQuantiles {
        int lq0 = 0;//min
        int lq25 = 0;
        int lq50 = 0;//median
        int lq75 = 0;
        int lq100 = 0;//max
    };
    explicit LongestSubPath(const InputGraph &inputGraph, const std::vector<ClusteredOriginDestination> &odPairs)
            : inputGraph(inputGraph), odPairs(odPairs), quanInFstIter(odPairs.size()), longestSubPathInFstIer(odPairs.size()){}

    void measureFirstIteration(const std::vector<std::vector<int32_t>> &odPairPaths,
                               const AlignedVector<int> &trafficFlows) {
        for (int i = 0; i < odPairs.size(); ++i) {
            const auto &path = odPairPaths[i];
            quanInFstIter[i] = StaticalFunction::weightedQuantile(path, inputGraph, trafficFlows);
            longestSubPathInFstIer[i].lq0 = getLongestSubPath(path, trafficFlows, quanInFstIter[i].q0);
            longestSubPathInFstIer[i].lq25 = getLongestSubPath(path, trafficFlows, quanInFstIter[i].q25);
            longestSubPathInFstIer[i].lq50 = getLongestSubPath(path, trafficFlows, quanInFstIter[i].q50);
            longestSubPathInFstIer[i].lq75 = getLongestSubPath(path, trafficFlows, quanInFstIter[i].q75);
            longestSubPathInFstIer[i].lq100 = getLongestSubPath(path, trafficFlows, quanInFstIter[i].q100);
        }
    }

    void measureLastIterationAndWriteOutput(const std::vector<std::vector<int32_t>> &odPairPaths,
                                            const AlignedVector<int> &trafficFlows, const std::string &anaFileName) {
        double sumMin = 0;
        double sum25 = 0;
        double sum50 = 0;
        double sum75 = 0;
        double sumMax = 0;

        auto odAnaFileName = anaFileName + "_longestSubPath.csv";
        std::ofstream odAnaFile;
        odAnaFile.open(odAnaFileName);
        if (!odAnaFile.good())
            throw std::invalid_argument("file cannot be opened -- '" + odAnaFileName + "'");
        odAnaFile << "origin,destination,,first iteration,,,,,last iteration,,,,,ratio" << std::endl ;
        LongestSubPathOfDifferentQuantiles longestSubPathInLast;

        for (int i = 0; i < odPairs.size(); ++i) {
            const auto &path = odPairPaths[i];
            assert((path.empty() && quanInFstIter[i] == StaticalFunction::Quantiles()) ||
                   (!path.empty() && quanInFstIter[i] != StaticalFunction::Quantiles()));
            odAnaFile << std::to_string(odPairs[i].origin) << "," << std::to_string(odPairs[i].destination) << ",";
            const auto &longestSubPathInFst = longestSubPathInFstIer[i];
            longestSubPathInLast.lq0 = getLongestSubPath(path, trafficFlows, quanInFstIter[i].q0);
            longestSubPathInLast.lq25 = getLongestSubPath(path, trafficFlows, quanInFstIter[i].q25);
            longestSubPathInLast.lq50 = getLongestSubPath(path, trafficFlows, quanInFstIter[i].q50);
            longestSubPathInLast.lq75 = getLongestSubPath(path, trafficFlows, quanInFstIter[i].q75);
            longestSubPathInLast.lq100 = getLongestSubPath(path, trafficFlows, quanInFstIter[i].q100);
            odAnaFile <<"," << longestSubPathInFst.lq0 << "," << longestSubPathInFst.lq25 <<"," <<longestSubPathInFst.lq50 << ","
                      << longestSubPathInFst.lq75 << "," << longestSubPathInFst.lq100 << ",,"
                      <<  longestSubPathInLast.lq0 << ", " << longestSubPathInLast.lq25 << ", " <<longestSubPathInLast.lq50 << ", "
                      <<longestSubPathInLast.lq75 << ", " <<longestSubPathInLast.lq100 << ",,";

            if (!path.empty()) {

                double ratioMin =
                        getRatio(static_cast<double>(longestSubPathInFst.lq0), static_cast<double>(longestSubPathInLast.lq0));
                double ratio25 =
                        getRatio(static_cast<double>(longestSubPathInFst.lq25), static_cast<double>(longestSubPathInLast.lq25));
                double ratio50 =
                        getRatio(static_cast<double>(longestSubPathInFst.lq50), static_cast<double>(longestSubPathInLast.lq50));
                double ratio75 =
                        getRatio(static_cast<double>(longestSubPathInFst.lq75), static_cast<double>(longestSubPathInLast.lq75));
                double ratioMax =
                        getRatio(static_cast<double>(longestSubPathInFst.lq100), static_cast<double>(longestSubPathInLast.lq100));
                sumMin += ratioMin;
                sum25 += ratio25;
                sum50 += ratio50;
                sum75 += ratio75;
                sumMax += ratioMax;
                odAnaFile << ratioMin << ", " << ratio25 << ", " << ratio50 << ", " << ratio75 << ", " << ratioMax
                          << "\n";
            } else {
                odAnaFile << "1.0,1.0,1.0,1.0,1.0\n";
            }
        }
        odAnaFile.flush();

        const auto count = static_cast<double>(odPairs.size());
        odAnaFile << ",,,,,,,,,,,,mean,,,,,,,,,,,,,," << sumMin / count * 1.0 << "," << sum25 / count * 1.0 << "," << sum50 / count * 1.0 << ","
                  << sum75 / count * 1.0 << "," << sumMax / count * 1.0 << std::endl;
    }

private:

    const InputGraph &inputGraph;
    const std::vector<ClusteredOriginDestination> &odPairs;
    std::vector<StaticalFunction::Quantiles> quanInFstIter;
    std::vector<LongestSubPathOfDifferentQuantiles> longestSubPathInFstIer;

    int getLongestSubPath(const std::vector<int32_t> &path, const AlignedVector<int> &trafficFLows, const int threshold) {
        int max = 0;
        int tmpTraTime = 0;
        if (!path.empty()) {
            for (auto e : path) {
                if (trafficFLows[e] >= threshold) {
                    tmpTraTime += inputGraph.travelTime(e);
                    if (tmpTraTime > max) {
                        max  = tmpTraTime;
                    } else {
                        tmpTraTime = 0;
                    }
                }
            }
        }

        return max;

    }

    double getRatio(double valueInFst, double valueInLast) {
        if (valueInFst == 0) {
            return 1;
        }
        return valueInLast/valueInFst;
    }
};
#endif //ROUTINGFRAMEWORK_LONGESTSUBPATH_H

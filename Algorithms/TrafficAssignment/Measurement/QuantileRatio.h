//
// Created by junjun on 19.11.22.
//

#ifndef ROUTINGFRAMEWORK_QUANTILERATIO_H
#define ROUTINGFRAMEWORK_QUANTILERATIO_H

#include "Algorithms/TrafficAssignment/Measurement/IMeasure.h"
#include "Algorithms/TrafficAssignment/Measurement/StaticalFunction.h"
#include "Tools/Simd/AlignedVector.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include <fstream>
#include <vector>
#include <iostream>

template<typename InputGraph>
class QuantileRatio : public IMeasure {
public :
    explicit QuantileRatio(const InputGraph &inputGraph, const std::vector<ClusteredOriginDestination> &odPairs)
            : inputGraph(inputGraph), odPairs(odPairs), quanInFstIter(odPairs.size()) {}

    void measureFirstIteration(const std::vector<std::vector<int32_t>> &odPairPaths,
                               const AlignedVector<int> &trafficFlows) {
        for (int i = 0; i < odPairs.size(); ++i) {
            const auto &path = odPairPaths[i];
            quanInFstIter[i] = StaticalFunction::weightedQuantile(path, inputGraph, trafficFlows);

        }
    }

    void measureLastIterationAndWriteOutput(const std::vector<std::vector<int32_t>> &odPairPaths,
                                            const AlignedVector<int> &trafficFlows, const std::string &anaFileName) {
        double sumMin = 0;
        double sum25 = 0;
        double sum50 = 0;
        double sum75 = 0;
        double sumMax = 0;

        auto odAnaFileName = anaFileName + "_QuantileRatio.csv";
        std::ofstream odAnaFile;
        odAnaFile.open(odAnaFileName);
        if (!odAnaFile.good())
            throw std::invalid_argument("file cannot be opened -- '" + odAnaFileName + "'");
        odAnaFile << origin,destination,,first iteration,,,,,last iteration,,,,,ratio" <<std::endl;

        for (int i = 0; i < odPairs.size(); ++i) {
            const auto &path = odPairPaths[i];
            assert((path.empty() && quanInFstIter[i] == StaticalFunction::Quantiles()) ||
                   (!path.empty() && quanInFstIter[i] != StaticalFunction::Quantiles()));
            odAnaFile << std::to_string(odPairs[i].origin) << "," << std::to_string(odPairs[i].destination) << ",";

            const auto &quantilesAtFirstIter = quanInFstIter[i];
            const auto &quantilesAtLastIter = StaticalFunction::weightedQuantile(path, inputGraph,
                                                                                 trafficFlows);
            odAnaFile << static_cast<double>(quantilesAtFirstIter.q0) << ", " << static_cast<double>(quantilesAtFirstIter.q25) << ", "
                      << static_cast<double>(quantilesAtFirstIter.q50) << ", "<< static_cast<double>(quantilesAtFirstIter.q75) << ", "<< static_cast<double>(quantilesAtFirstIter.q100) << ", , "
                      << static_cast<double>(quantilesAtLastIter.q0) << ", " << static_cast<double>(quantilesAtLastIter.q25) << ", "
                      << static_cast<double>(quantilesAtLastIter.q50) << ", "<< static_cast<double>(quantilesAtLastIter.q75) << ", "<< static_cast<double>(quantilesAtLastIter.q100) << ", , ";
            if (!path.empty()) {
                double ratioMin =
                        static_cast<double>(quantilesAtLastIter.q0) / static_cast<double>(quantilesAtFirstIter.q0);
                double ratio25 = static_cast<double>(quantilesAtLastIter.q25) /
                                 static_cast<double>(quantilesAtFirstIter.q25);
                double ratio50 = static_cast<double>(quantilesAtLastIter.q50) /
                                 static_cast<double>(quantilesAtFirstIter.q50);
                double ratio75 = static_cast<double>(quantilesAtLastIter.q75) /
                                 static_cast<double>(quantilesAtFirstIter.q75);
                double ratioMax = static_cast<double>(quantilesAtLastIter.q100) /
                                  static_cast<double>(quantilesAtFirstIter.q100);
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

};


#endif //ROUTINGFRAMEWORK_QUANTILERATIO_H
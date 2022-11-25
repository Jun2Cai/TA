//
// Created by junjun on 25.11.22.
//

#ifndef ROUTINGFRAMEWORK_QUANTILEDISTRIBUTION_H
#define ROUTINGFRAMEWORK_QUANTILEDISTRIBUTION_H
#include "Algorithms/TrafficAssignment/Measurement/StaticalFunction.h"
#include "Tools/Simd/AlignedVector.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include <fstream>
#include <vector>
#include <iostream>

template<typename InputGraph>
class QuantileDistribution  {
public :
    struct DistributionOnQuantiles{
    double dq0 = 0;
    double dq25 = 0;
    double dq50 = 0;
    double dq75 = 0;
    double dq100 = 0;
};
    explicit QuantileDistribution(const InputGraph &inputGraph, const std::vector<ClusteredOriginDestination> &odPairs)
            : inputGraph(inputGraph), odPairs(odPairs), quanInFstIter(odPairs.size()), disQuanInFstIter(odPairs.size()) {}

    void measureFirstIteration(const std::vector<std::vector<int32_t>> &odPairPaths,
                               const AlignedVector<int> &trafficFlows) {
        for (int i = 0; i < odPairs.size(); ++i) {
            const auto &path = odPairPaths[i];
            quanInFstIter[i] = StaticalFunction::weightedQuantile(path, inputGraph, trafficFlows);
            disQuanInFstIter[i] = {
                    getDistribution(path, trafficFlows,  quanInFstIter[i].q0),
                    getDistribution(path, trafficFlows,  quanInFstIter[i].q25),
                    getDistribution(path, trafficFlows,  quanInFstIter[i].q50),
                    getDistribution(path, trafficFlows,  quanInFstIter[i].q75),
                    getDistribution(path, trafficFlows,  quanInFstIter[i].q100),
            };

        }
    }

    void measureLastIterationAndWriteOutput(const std::vector<std::vector<int32_t>> &odPairPaths,
                                            const AlignedVector<int> &trafficFlows, const std::string &anaFileName) {
        double sumMin = 0;
        double sum25 = 0;
        double sum50 = 0;
        double sum75 = 0;
        double sumMax = 0;

        auto odAnaFileName = anaFileName + "_distribution.csv";
        std::ofstream odAnaFile;
        odAnaFile.open(odAnaFileName);
        if (!odAnaFile.good())
            throw std::invalid_argument("file cannot be opened -- '" + odAnaFileName + "'");
        odAnaFile << "origin,destination,,first iteration,,,,,last iteration,,,,,,ratio" <<std::endl;

        for (int i = 0; i < odPairs.size(); ++i) {
            const auto &path = odPairPaths[i];
            assert((path.empty() && quanInFstIter[i] == StaticalFunction::Quantiles()) ||
                   (!path.empty() && quanInFstIter[i] != StaticalFunction::Quantiles()));
            odAnaFile << std::to_string(odPairs[i].origin) << "," << std::to_string(odPairs[i].destination) << ",";

            const auto &quantilesAtFirstIter = quanInFstIter[i];
            const auto distributionOnQuantilesAtFst = disQuanInFstIter[i];
            DistributionOnQuantiles distributionOnQuantilesAtLast = {
                    getDistribution(path, trafficFlows,  quantilesAtFirstIter.q0),
                    getDistribution(path, trafficFlows,  quantilesAtFirstIter.q25),
                    getDistribution(path, trafficFlows,  quantilesAtFirstIter.q50),
                    getDistribution(path, trafficFlows,  quantilesAtFirstIter.q75),
                    getDistribution(path, trafficFlows,  quantilesAtFirstIter.q100)};

            odAnaFile <<"," << static_cast<double>(distributionOnQuantilesAtFst.dq0) << ", " << static_cast<double>(distributionOnQuantilesAtFst.dq25) << ", "
                      << static_cast<double>(distributionOnQuantilesAtFst.dq50) << ", "<< static_cast<double>(distributionOnQuantilesAtFst.dq75) << ", "<< static_cast<double>(distributionOnQuantilesAtFst.dq100) << ", , "
                    << static_cast<double>(distributionOnQuantilesAtLast.dq0) << ", " << static_cast<double>(distributionOnQuantilesAtLast.dq25) << ", "
                    << static_cast<double>(distributionOnQuantilesAtLast.dq50) << ", "<< static_cast<double>(distributionOnQuantilesAtLast.dq75) << ", "<< static_cast<double>(distributionOnQuantilesAtLast.dq100) << ", , ";
            if (!path.empty()) {
                double ratioMin =
                        static_cast<double>(distributionOnQuantilesAtLast.dq0)
                        / static_cast<double>(distributionOnQuantilesAtFst.dq0);
                double ratio25 =  static_cast<double>(distributionOnQuantilesAtLast.dq25)
                                  / static_cast<double>(distributionOnQuantilesAtFst.dq25);
                double ratio50 =  static_cast<double>(distributionOnQuantilesAtLast.dq50)
                                  / static_cast<double>(distributionOnQuantilesAtFst.dq50);
                double ratio75 =  static_cast<double>(distributionOnQuantilesAtLast.dq75)
                                  / static_cast<double>(distributionOnQuantilesAtFst.dq75);
                double ratioMax =  static_cast<double>(distributionOnQuantilesAtLast.dq100)
                                   / static_cast<double>(distributionOnQuantilesAtFst.dq100);
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
        odAnaFile << ",,,,,,,,,,,,,,mean," << sumMin / count * 1.0 << "," << sum25 / count * 1.0 << "," << sum50 / count * 1.0 << ","
                  << sum75 / count * 1.0 << "," << sumMax / count * 1.0 << std::endl;
    }

private:

    const InputGraph &inputGraph;
    const std::vector<ClusteredOriginDestination> &odPairs;
    std::vector<StaticalFunction::Quantiles> quanInFstIter;
    std::vector<DistributionOnQuantiles> disQuanInFstIter;

    double getDistribution(const std::vector<int32_t> &path,  const AlignedVector<int> &trafficFlows, int threshold ) {
        int totalFlows = 0;
        int sumOfFlows = 0;
        for (auto e : path) {
            if (trafficFlows[e] >= threshold) {
                sumOfFlows += trafficFlows[e];
            }
            totalFlows += trafficFlows[e];
        }
        if (totalFlows == 0) {
            return 1;
        }
        return sumOfFlows*1.0/totalFlows;
    }

};

#endif //ROUTINGFRAMEWORK_QUANTILEDISTRIBUTION_H

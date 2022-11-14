//
// Created by junjun on 30.10.22.
//

#ifndef ROUTINGFRAMEWORK_COMPARETRAVELTIME_H
#define ROUTINGFRAMEWORK_COMPARETRAVELTIME_H



#include "Algorithms/TrafficAssignment/Measurement/IMeasure.h"
#include "Tools/Simd/AlignedVector.h"
#include <fstream>
#include <vector>
#include <map>
#include <iostream>

        class CompareTravelTime : public IMeasure {

public :

    template <typename InputGraph>
    void measure(std::string anaFileName,std::map<std::string, std::vector<int32_t>> odPairPath, InputGraph inputGraph,  int numIterations) {

        //use static values to store values from the first iteration
        static std::vector<double> travelTimeInFstIter;
        int i = 0;
        double totalRatio = 0;

        std::ofstream odAnaFile;


        if (numIterations > 1) {
            odAnaFile.open(anaFileName + "_compareTravelTieme.csv");
            odAnaFile << "origin,destination,travel time in first iteration, travel time in last iteration,ratio"
                      << std::endl;
        }
        for (auto iter = odPairPath.begin(); iter != odPairPath.end(); iter++) {
            auto pair = iter->first;
            auto path = iter->second;
            if (numIterations > 1) {
                odAnaFile << pair << ",";
            }
            double travelTime = 0;
            for (auto e: path) {
                    travelTime += inputGraph.travelTime(e);

            }

            if (numIterations == 1) {
                travelTimeInFstIter.push_back(travelTime);
            } else {
                odAnaFile << travelTimeInFstIter.at(i)  << "," << travelTime << ",";
                double ratio = 1;
                if (travelTimeInFstIter.at(i) > 0) {
                    ratio = travelTime / travelTimeInFstIter.at(i);
                }
                odAnaFile << ratio << std::endl;
                i++;
                totalRatio += ratio;
            }
        }

        if (numIterations > 1) {
            odAnaFile << "mean" << ",,,," << totalRatio / (i + 1);
        }


    }


};

#endif //ROUTINGFRAMEWORK_COMPARETRAVELTIME_H

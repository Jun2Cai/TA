//
// Created by junjun on 05.10.22.
//

#ifndef ROUTINGFRAMEWORK_MEANMEASURE_H
#define ROUTINGFRAMEWORK_MEANMEASURE_H

#include <fstream>
#include <vector>
#include <map>
#include "Algorithms/TrafficAssignment/Measurement/IMeasure.h"
#include "Tools/Simd/AlignedVector.h"


class MeanMeasure : public IMeasure {
public :
    template <typename InputGraph>
    void measure(std::string anaFileName,std::map<std::string, std::vector<int32_t>> odPairPath, InputGraph inputGraph, AlignedVector<int> trafficFlows, int numIterations) {
        //use static values to store the information of the first iteration
        static std::vector<std::vector<double>>  storeInfoFromFirstIteration;
        static std::vector<double>  sumInfo;
        std::ofstream odAnaFile;
        if (numIterations > 1) {
            odAnaFile.open(anaFileName + "_mean.csv");
            odAnaFile << "origin,destination,first iteration,0_20,20_40,40_60,60_80,80_100,last iteration,0_20,20_40,40_60,60_80,80_100"<<std::endl;
        }

}


};


#endif //ROUTINGFRAMEWORK_MEANMEASURE_H

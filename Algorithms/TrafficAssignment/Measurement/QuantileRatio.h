//
// Created by junjun on 19.11.22.
//

#ifndef ROUTINGFRAMEWORK_QUANTILERATIO_H
#define ROUTINGFRAMEWORK_QUANTILERATIO_H

#include "Algorithms/TrafficAssignment/Measurement/IMeasure.h"
#include "Algorithms/TrafficAssignment/Measurement/StaticalFunction.h"
#include "Tools/Simd/AlignedVector.h"
#include <fstream>
#include <vector>
#include <map>
#include <iostream>
    class QuantileRatio : public IMeasure {
    public :

        template <typename InputGraph>
        void measure(std::string anaFileName, std::map<std::string, std::vector<int32_t>> odPairPath, InputGraph inputGraph, AlignedVector<int> trafficFlows, int numIterations) {
        static std::map<std::string, std::vector<int>> quanInFstIter;
        if (numIterations == 1) {
            for (auto iter = odPairPath.begin(); iter!= odPairPath.end(); iter++) {
                auto pair = iter->first;
                auto path = iter->second;
                if (path.size() > 0) {

                     quanInFstIter[pair] = StaticalFunction::weightedQuantile(path, inputGraph, trafficFlows) ;

                }
            }
        } else {

            double sumMin = 0;
            double sum25 = 0;
            double sum50 = 0;
            double sum75 = 0;
            double sumMax =0;
            int count = 0;
            std::ofstream odAnaFile;
            odAnaFile.open(anaFileName+"_QuantileRatio.csv");
            for (auto iter = odPairPath.begin(); iter != odPairPath.end(); iter++) {
                auto pair = iter->first;
                auto path = iter->second;
                if (path.size() > 0 && (quanInFstIter.find(pair) != quanInFstIter.end())) {
                    odAnaFile << pair << ",";
                    double ratioMin = StaticalFunction::weightedQuantile(path, inputGraph, trafficFlows).at(0)/quanInFstIter[pair].at(0);
                    double ratio25 = StaticalFunction::weightedQuantile(path, inputGraph, trafficFlows).at(1)/quanInFstIter[pair].at(1);
                    double ratio50 = StaticalFunction::weightedQuantile(path, inputGraph, trafficFlows).at(2)/quanInFstIter[pair].at(2);
                    double ratio75 = StaticalFunction::weightedQuantile(path, inputGraph, trafficFlows).at(3)/quanInFstIter[pair].at(3);
                    double ratioMax = StaticalFunction::weightedQuantile(path, inputGraph, trafficFlows).at(4)/quanInFstIter[pair].at(4);
                    sumMin += ratioMin;
                    sum25 += ratio25;
                    sum50 += ratio50;
                    sum75 += ratio75;
                    sumMax += ratioMax;
                    count++;
                    odAnaFile << ratioMin << ", " <<ratio25 << ", " <<ratio50 << ", " <<ratio75 << ", " << ratioMax << std::endl;

                }
            }

            odAnaFile << "," << std::endl;
            odAnaFile << sumMin/count*1.0 <<sum25/count*1.0 <<sum50/count*1.0 <<sum75/count*1.0 <<sumMax/count*1.0 ;
        }
    }



};



#endif //ROUTINGFRAMEWORK_QUANTILERATIO_H

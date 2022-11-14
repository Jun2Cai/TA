//
// Created by junjun on 05.10.22.
//

#ifndef ROUTINGFRAMEWORK_PERCENTAGEOFDIFFERENTSAT_H
#define ROUTINGFRAMEWORK_PERCENTAGEOFDIFFERENTSAT_H

#include <fstream>
#include <vector>
#include <map>
#include "Algorithms/TrafficAssignment/Measurement/IMeasure.h"
#include "Tools/Simd/AlignedVector.h"


class PercentageOfDifferentSat : public IMeasure {
public :
    template <typename InputGraph>
    void measure(std::string anaFileName,std::map<std::string, std::vector<int32_t>> odPairPath, InputGraph inputGraph, int numIterations) {
        //use static values to store the information of the first iteration
        static std::vector<std::vector<double>>  storeInfoFromFirstIteration;
        static std::vector<double>  sumInfo;
        std::ofstream odAnaFile;
        int i = 0;
        if (numIterations > 1) {
            odAnaFile.open(anaFileName + "_percentageOfDifferentSat.csv");
            odAnaFile << "origin,destination,first iteration,0_20,20_40,40_60,60_80,80_100,last iteration,0_20,20_40,40_60,60_80,80_100"<<std::endl;
        }
        double sumOfSatInLength0 = 0;
        double sumOfSatInLength1 = 0;
        double sumOfSatInLength2 = 0;
        double sumOfSatInLength3 = 0;
        double sumOfSatInLength4 = 0;
        int sumOfMeaningfulPair = 0;
        for(auto iter = odPairPath.begin(); iter != odPairPath.end(); iter++) {
            auto pair = iter -> first;
            auto path = iter -> second;


            int totalLength = 0;
            int len0 = 0;
            int len1 = 0;
            int len2 = 0;
            int len3 = 0;
            int len4 = 0;
            for(auto e : path) {
                if (path.size() > 0) {
                    auto vol = trafficFlows[e];
                    auto sat = (double)vol / (double)inputGraph.capacity(e);
                    if (sat >=0 && sat <0.2) {
                        len0 += inputGraph.length(e);

                    } else if (sat >= 0.2 && sat < 0.4) {
                        len1 += inputGraph.length(e);
                    } else if (sat >= 0.4 && sat < 0.6) {
                        len2 += inputGraph.length(e);
                    } else if (sat >= 0.6 && sat < 0.8) {
                        len3 += inputGraph.length(e);
                    } else if (sat > 0.8) {
                        len4 += inputGraph.length(e);

                    }
                    totalLength += inputGraph.length(e);
                }
            }


            if (path.size() > 0 && totalLength!=0 ) {
                if (numIterations > 1) {
                    odAnaFile << pair << ",," ;
                }
            }



            double size = (double)path.size();
            double length = (double)totalLength ;

            if (size > 0 && length  > 0 ) {

                sumOfMeaningfulPair++;
                sumOfSatInLength0 += len0/length;
                sumOfSatInLength1 += len1/length;
                sumOfSatInLength2 += len2/length;
                sumOfSatInLength3 += len3/length;
                sumOfSatInLength4 += len4/length;
                if (numIterations > 1) {

                    auto edge = storeInfoFromFirstIteration.at(i);

                    odAnaFile << edge.at(0) << "," <<  edge.at(1) << "," <<  edge.at(2) << "," << edge.at(3) << "," << edge.at(4) << ",," <<
                              len0/length << ","<< len1/length  << ","<< len2/length << ","<< len3/length << ","<< len4/length <<std::endl;
                    i++;
                } else {
                    std::vector<double> sat5;
                    sat5.push_back(len0/length);
                    sat5.push_back(len1/length);
                    sat5.push_back(len2/length);
                    sat5.push_back(len3/length);
                    sat5.push_back(len4/length);
                    storeInfoFromFirstIteration.push_back(sat5);
                }
            }
        }

        if (numIterations == 1) {
            sumInfo.push_back(sumOfSatInLength0/sumOfMeaningfulPair);
            sumInfo.push_back(sumOfSatInLength1/sumOfMeaningfulPair);
            sumInfo.push_back(sumOfSatInLength2/sumOfMeaningfulPair);
            sumInfo.push_back(sumOfSatInLength3/sumOfMeaningfulPair);
            sumInfo.push_back(sumOfSatInLength4/sumOfMeaningfulPair);
        } else {
            odAnaFile << std::endl;
            odAnaFile << "mean,,,"<< sumInfo.at(0) << "," << sumInfo.at(1) << "," << sumInfo.at(2) << "," << sumInfo.at(3) << ","<< sumInfo.at(4) << ", ," << sumOfSatInLength0/sumOfMeaningfulPair<<","<< sumOfSatInLength1/sumOfMeaningfulPair<<","<< sumOfSatInLength2/sumOfMeaningfulPair<<","<< sumOfSatInLength3/sumOfMeaningfulPair<<","<< sumOfSatInLength4/sumOfMeaningfulPair<<std::endl;
            odAnaFile << "increase percentage,,,"<< (sumOfSatInLength0/sumOfMeaningfulPair) / sumInfo.at(0)
            << "," << (sumOfSatInLength1/sumOfMeaningfulPair) / sumInfo.at(1)
            << "," <<(sumOfSatInLength2/sumOfMeaningfulPair) / sumInfo.at(2)
            << "," << (sumOfSatInLength3/sumOfMeaningfulPair) / sumInfo.at(3)
            << "," << (sumOfSatInLength4/sumOfMeaningfulPair) / sumInfo.at(4)
            << std::endl;
        }


    }
    void setTrafficFlows(AlignedVector<int> flow) {
        trafficFlows = flow;
    }

private:
    AlignedVector<int> trafficFlows;
};


#endif //ROUTINGFRAMEWORK_PERCENTAGEOFDIFFERENTSAT_H

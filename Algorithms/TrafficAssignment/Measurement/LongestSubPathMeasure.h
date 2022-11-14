//
// Created by junjun on 05.10.22.
//

#ifndef ROUTINGFRAMEWORK_LONGESTSUBPATHMEASURE_H
#define ROUTINGFRAMEWORK_LONGESTSUBPATHMEASURE_H


#include "Algorithms/TrafficAssignment/Measurement/IMeasure.h"
#include "Tools/Simd/AlignedVector.h"
#include <fstream>
#include <vector>
#include <map>
#include <iostream>

class LongestSubPathMeasure : public IMeasure {

public :

    template <typename InputGraph>
    void measure(std::string anaFileName,std::map<std::string, std::vector<int32_t>> odPairPath, InputGraph inputGraph, AlignedVector<int> trafficFlows, int numIterations) {
        //use static values to store values from the first iteration
        static std::vector<int> longestSubPath20FromFirstIter;
        static std::vector<int> longestSubPath40FromFirstIter;
        static std::vector<int> longestSubPath60FromFirstIter;
        static std::vector<int> longestSubPath80FromFirstIter;
        static std::vector<int> longestSubPath100FromFirstIter;


        std::ofstream odAnaFile;

        if (numIterations > 1) {
            odAnaFile.open(anaFileName+"_longestSubPath.csv");
            odAnaFile << "origin,destination,first iteration, >20, >40, >60, >80, >100,last iteration,>20, >40, >60, >80, >100" << std::endl;
        }


        double sumOfPath20 = 0;
        double sumOfPath40 = 0;
        double sumOfPath60 = 0;
        double sumOfPath80 = 0;
        double sumOfPath100 = 0;

        int i = 0;

        for(auto iter = odPairPath.begin(); iter != odPairPath.end(); iter++) {

            auto pair = iter->first;
            auto path = iter->second;



                int longestSubPath20 = 0;
                int longestSubPath40 = 0;
                int longestSubPath60 = 0;
                int longestSubPath80 = 0;
                int longestSubPath100 = 0;



                for (int j = 1; j <= 5 ; j++) {
                    int longestSubPath = 0;
                    int temp = 0;
                    for(auto e : path) {
                        auto vol = trafficFlows[e];
                        auto sat = (double)vol /  (double)inputGraph.capacity(e);
                        if (sat > j*0.2) {
                            temp += inputGraph.length(e);
                            if (temp > longestSubPath) {
                                longestSubPath = temp;
                            }
                        } else {
                            temp = 0;
                        }
                    }

                    switch (j) {
                        case 1:
                            longestSubPath20 = longestSubPath;
                            break;
                        case 2:
                            longestSubPath40 = longestSubPath;
                            break;
                        case 3:
                            longestSubPath60 = longestSubPath;
                            break;
                        case 4:
                            longestSubPath80 = longestSubPath;
                            break;
                        case 5:
                            longestSubPath100 = longestSubPath;
                            break;
                    }


                }

            if (path.size() > 0) {
                if (numIterations > 1) {
                    odAnaFile << pair << ",," ;
                }
            }











                    if (numIterations == 1) {
                        longestSubPath20FromFirstIter.push_back(longestSubPath20);
                        longestSubPath40FromFirstIter.push_back(longestSubPath40);
                        longestSubPath60FromFirstIter.push_back(longestSubPath60);
                        longestSubPath80FromFirstIter.push_back(longestSubPath80);
                        longestSubPath100FromFirstIter.push_back(longestSubPath100);
                    } else {

                            odAnaFile  <<longestSubPath20FromFirstIter.at(i) << "," << longestSubPath40FromFirstIter.at(i) << ","
                         << longestSubPath60FromFirstIter.at(i) << "," << longestSubPath80FromFirstIter.at(i) << "," << longestSubPath100FromFirstIter.at(i) << ",," <<
                           longestSubPath20 << "," << longestSubPath40 << "," <<longestSubPath60 << "," <<longestSubPath80 << "," <<longestSubPath100 << "," << std::endl;

                            sumOfPath20 += longestSubPath20FromFirstIter.at(i) == 0 ? 1 : (double) longestSubPath20/longestSubPath20FromFirstIter.at(i);
                            sumOfPath40 += longestSubPath40FromFirstIter.at(i) == 0 ? 1 :(double) longestSubPath40/longestSubPath40FromFirstIter.at(i);
                            sumOfPath60 += longestSubPath60FromFirstIter.at(i) == 0 ? 1 :(double) longestSubPath60/longestSubPath60FromFirstIter.at(i);
                            sumOfPath80 += longestSubPath80FromFirstIter.at(i) == 0 ? 1 :(double) longestSubPath80/longestSubPath80FromFirstIter.at(i);
                            sumOfPath100 += longestSubPath100FromFirstIter.at(i) == 0 ? 1 :(double) longestSubPath100/longestSubPath100FromFirstIter.at(i);
                            i++;

                    }








        }



        if (numIterations > 1)  {
            odAnaFile << std::endl;
            odAnaFile << " ,,increase percentage: sum(last/first)/size," << sumOfPath20 / (i+1) << ","
                    << sumOfPath40 / (i+1) << ","
                    << sumOfPath60 /(i+1) << ","
                    << sumOfPath80 / (i+1) << ","
                    << sumOfPath100 / (i+1) << "," << std::endl;

        }





}

void setTrafficFlows(AlignedVector<int> flow) {
    this->trafficFlows = flow;
}
private:
    AlignedVector<int> trafficFlows;




};
#endif //ROUTINGFRAMEWORK_LONGESTSUBPATHMEASURE_H

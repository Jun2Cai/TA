//
// Created by junjun on 13.11.22.
//

#ifndef ROUTINGFRAMEWORK_INCREASEDFLOW_H
#define ROUTINGFRAMEWORK_INCREASEDFLOW_H

#include "Tools/Simd/AlignedVector.h"
#include <fstream>
#include <vector>
#include <map>
#include <iostream>
#include "Algorithms/TrafficAssignment/Measurement/IMeasure.h"
class IncreasedFlow : IMeasure {
public :

    template <typename InputGraph>
    void measure(std::string anaFileName, InputGraph graph, int numIterations) {
        //use static values to store values from the first iteration


        std::ofstream odAnaFile;
        static AlignedVector<int> trafficFlowInFstIter;
        int i = 0;
        double sum = 0;
        if (numIterations > 1) {

            odAnaFile.open(anaFileName + "_increased.csv");
            odAnaFile << ",ratio" << std::endl;

            for (int u = 0; u < graph.numVertices() - 1 ; u++) {
                for (int e = graph.firstEdge(u); e < graph.lastEdge(u); ++e) {

                    if (trafficFlowInFstIter[e] > 0 ) {
                        double ratio = trafficFlows[e]*(1.0)/trafficFlowInFstIter[e];
                        if (ratio > 1) {
                            odAnaFile << "," << ratio << std::endl;
                            i++;
                            sum += ratio;
                        }

                    }




                }
            }

            odAnaFile  << std::endl;
            odAnaFile << "mean," << (sum/i);
        } else {
            trafficFlowInFstIter = trafficFlows;

        }






    }



    void setTrafficFlows(AlignedVector<int> flow) {
        this->trafficFlows = flow;
    }
private:
    AlignedVector<int> trafficFlows;




};


#endif //ROUTINGFRAMEWORK_INCREASEDFLOW_H

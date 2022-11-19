//
// Created by junjun on 18.11.22.
//

#ifndef ROUTINGFRAMEWORK_STATICALFUNCTION_H
#define ROUTINGFRAMEWORK_STATICALFUNCTION_H

#include <vector>
#include "Tools/Simd/AlignedVector.h"
#include <map>
#include <algorithm>
#include <math.h>
class StaticalFunction {

public :
     template <typename InputGraph>
      static std::vector<int32_t> weightedQuantile(std::vector<int32_t> path, InputGraph inputGraph, AlignedVector<int> trafficFlows) {
      std::vector<int32_t> quantile;
      quantile.resize(100);
      //flow, travelTime, edge
      std::vector<std::pair<int,int>> infoOfEachE;
      int totalTime = 0;
      for (auto e : path) {
          std::pair<int,int> pair(trafficFlows[e], inputGraph.travelTime(e));
          infoOfEachE.push_back(pair);
          totalTime += inputGraph.travelTime(e);
      }


      std::sort(infoOfEachE.begin(), infoOfEachE.end());
      int quantile25 = ceil(totalTime * 0.25);
      int quantile50 = ceil(totalTime * 0.5);
      int quantile75 = ceil(totalTime * 0.75);
      quantile.push_back(std::get<0>(infoOfEachE.front()));
      quantile.push_back(getQuantile(infoOfEachE,quantile25));
      quantile.push_back(getQuantile(infoOfEachE,quantile50));
      quantile.push_back(getQuantile(infoOfEachE,quantile75));
      quantile.push_back(std::get<0>(infoOfEachE.back()));
      quantile.push_back(1);




      return quantile;
}

private:
    static int getQuantile(std::vector<std::pair<int,int>> infoOfEachE, int quantileValue) {
        int tmp = 0;
        for (auto sorted: infoOfEachE) {
            tmp += std::get<1>(sorted);
            if (tmp > quantileValue) {
            return std::get<0>(sorted);
            }
         }
        return -999;
      }

    };
#endif //ROUTINGFRAMEWORK_STATICALFUNCTION_H

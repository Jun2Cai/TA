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

    struct Quantiles {

        friend bool operator==(const Quantiles &quan1, const Quantiles &quan2) {
            return quan1.q0 == quan2.q0 && quan1.q25 == quan2.q25 && quan1.q50 == quan2.q50 && quan1.q75 == quan2.q75 &&
                   quan1.q100 == quan2.q100;
        };

        friend bool operator!=(const Quantiles& q1, const Quantiles& q2) {return !(q1 == q2);};

        int q0 = 0; // min
        int q25 = 0;
        int q50 = 0; // median
        int q75 = 0;
        int q100 = 0; // max
    };

    template<typename InputGraph>
    static Quantiles weightedQuantile(const std::vector<int32_t> &path, const InputGraph &inputGraph,
                                      const AlignedVector<int> &trafficFlows) {

        if (path.empty())
            return {};

        //flow, travelTime, edge
        std::vector<std::pair<int, int>> infoOfEachE;
        int totalTime = 0;
        for (auto e: path) {
            std::pair<int, int> pair(trafficFlows[e], inputGraph.travelTime(e));
            infoOfEachE.push_back(pair);
            totalTime += inputGraph.travelTime(e);
        }


        std::sort(infoOfEachE.begin(), infoOfEachE.end());
        int quantile25 = ceil(totalTime * 0.25);
        int quantile50 = ceil(totalTime * 0.5);
        int quantile75 = ceil(totalTime * 0.75);

        Quantiles quantiles;
        quantiles.q0 = std::get<0>(infoOfEachE.front());
        assert(quantiles.q0 > 0);
        quantiles.q25 = getQuantile(infoOfEachE, quantile25);
        quantiles.q50 = getQuantile(infoOfEachE, quantile50);
        quantiles.q75 = getQuantile(infoOfEachE, quantile75);
        quantiles.q100 = std::get<0>(infoOfEachE.back());

        return quantiles;
    }

private:
    static int getQuantile(std::vector<std::pair<int, int>> infoOfEachE, int quantileValue) {
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

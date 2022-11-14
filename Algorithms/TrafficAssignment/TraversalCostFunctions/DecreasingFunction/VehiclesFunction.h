//
// Created by junjun on 21.10.22.
//

#ifndef ROUTINGFRAMEWORK_VEHICLESFUNCTION_H
#define ROUTINGFRAMEWORK_VEHICLESFUNCTION_H
#pragma once

#include <cassert>
#include <cmath>

#include <vectorclass/vectorclass.h>

template <typename GraphT>
class VehiclesFunction {
public:
    // Constructs an inverse travel cost function.
    VehiclesFunction(const GraphT& graph, double al, double be, double sc) : graph(graph) {
        alpha = al;
        beta = be;
        scale = sc;
        setAB(3,2);
    }

    // Returns the travel cost on edge e, given the flow x on e.
    double operator()(const int e, const double x) const {
        //a>b

        assert(x >= 0);
        return scale * graph.length(e) * (beta *pow(b,log(x+1)/log(a)) / (x + 1) + alpha) + (1-scale) * graph.travelTime(e);
    }

    // Returns the derivative of e's travel cost function at x.
    double derivative(const int e, const double x) const {
        assert(x >= 0);
        const double tmp = x + 1;
        return scale * graph.length(e) * beta *(log(b)- log(a))*pow(2,log(tmp)/log(a))/(log(a)*tmp*tmp);
    }

    // Returns the second derivative of e's travel cost function at x.
    double secondDerivative(const int e, const double x) const {
        assert(x >= 0);
        const double tmp = x + 1;
        return scale * graph.length(e) * beta * (log(b)- log(a)) * (log(b)-2*log(a)) * pow(b, log(tmp)/log(a))/(log(a)*log(a)*tmp*tmp*tmp);
    }

    void setAB(int al,int be) {
        a = al;
        b = be;
    }


private:
    const GraphT& graph; // The graph on whose edges we operate.
    double alpha;
    double beta;
    double scale;
    int a;
    int b;
};

#endif //ROUTINGFRAMEWORK_VEHICLESFUNCTION_H

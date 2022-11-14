#pragma once

#include <cassert>
#include <cmath>

#include <vectorclass/vectorclass.h>

template <typename GraphT>
class DelayedFunction{
public:

    DelayedFunction(const GraphT& graph, double al, double be, double sc) : graph(graph) {
        alpha = al;
        beta = be;
        scale = sc;
        a = 500;
    }

    // Returns the travel cost on edge e, given the flow x on e.
    //0.45
    double operator()(const int e, const double x) const {
        assert(x >= 0);
        if (x == 0) {
            return 1;
        }
        return scale * graph.length(e) * (beta * (1-exp(-a/x)) +  alpha)+ (1-scale) * graph.travelTime(e);



    }

    // Returns the derivative of e's travel cost function at x.
    double derivative(const int e, const double x) const {
        assert(x >= 0);
        if (x == 0) {
            return 0;
        }
        return scale * graph.length(e) * beta * (-a * exp(-a/x)) / (x * x);
    }

    // Returns the second derivative of e's travel cost function at x.
    double secondDerivative(const int e, const double x) const {
        assert(x >= 0);
        if (x == 0) {
            return 0;
        }
        return scale * graph.length(e) * beta * a * (2 * x - a) * exp(-a / x) / (x * x * x * x);
    }


private:
    const GraphT& graph; // The graph on whose edges we operate.
    double alpha;
    double beta;
    double scale;
    double a;
};

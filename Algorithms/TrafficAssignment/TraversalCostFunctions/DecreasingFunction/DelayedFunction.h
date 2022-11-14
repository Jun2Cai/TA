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
        a = 1;
    }

    // Returns the travel cost on edge e, given the flow x on e.
    //0.45
    double operator()(const int e, const double x) const {
        assert(x >= 0);
        double y = x/graph.capacity(e);
        return scale * graph.length(e) * (beta * (1-exp(-a/(y+0.01))) + alpha) + (1-scale) * graph.travelTime(e);
    }

    // Returns the derivative of e's travel cost function at x.
    double derivative(const int e, const double x) const {
        assert(x >= 0);
        double y = x/graph.capacity(e);
        double tmp = y + 0.05;
        return scale * graph.length(e) * beta * (-a * exp(-a/tmp)/(tmp*tmp));
    }

    // Returns the second derivative of e's travel cost function at x.
    double secondDerivative(const int e, const double x) const {
        assert(x >= 0);
        double y = x/graph.capacity(e);
        double tmp = y + 0.05;
        return scale * graph.length(e) * beta * (tmp * 2 * a * exp(-a / tmp) - a * a * exp(-a / tmp)) / (tmp*tmp*tmp*tmp);
    }


private:
    const GraphT& graph; // The graph on whose edges we operate.
    double alpha;
    double beta;
    double scale;
    double a;
};

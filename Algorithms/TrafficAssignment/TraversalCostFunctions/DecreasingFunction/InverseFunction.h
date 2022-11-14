#pragma once

#include <cassert>
#include <cmath>

#include <vectorclass/vectorclass.h>

// An inverse travel cost function, which decreases (rather than increases) as the flow increases.
// It should model the costs for operating public transit. The more people use public transit, the
// lower the costs per person are. However, a traffic assignment minimizing solely the operational
// cost may result in long detours for some passengers. Therefore, the (static) travel time also
// contributes to the travel cost.
template <typename GraphT>
class InverseFunction {
public:
    // Constructs an inverse travel cost function.
    InverseFunction(const GraphT& graph, double al, double be, double sc) : graph(graph) {
        alpha = al;
        beta = be;
        scale = sc;
    }

    // Returns the travel cost on edge e, given the flow x on e.
    double operator()(const int e, const double x) const {
        assert(x >= 0);
        double y = x / graph.capacity(e);
        return scale * graph.length(e) * (beta / (y + 1) + alpha) + (1-scale) * graph.travelTime(e);
    }

    // Returns the derivative of e's travel cost function at x.
    double derivative(const int e, const double x) const {
        assert(x >= 0);
        double y = x / graph.capacity(e);
        const double tmp = y + 1;
        return scale * graph.length(e) * -beta / (tmp * tmp);
    }

    // Returns the second derivative of e's travel cost function at x.
    double secondDerivative(const int e, const double x) const {
        assert(x >= 0);
        double y = x / graph.capacity(e);
        const double tmp = y + 1;
        return scale * graph.length(e) * 2 * beta / (tmp * tmp * tmp);
    }


private:
    const GraphT& graph; // The graph on whose edges we operate.
    double alpha;
    double beta;
    double scale;
};

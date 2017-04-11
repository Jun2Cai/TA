#pragma once

#include <array>
#include <type_traits>

#include "Algorithms/Dijkstra/BiDijkstra.h"
#include "Algorithms/Dijkstra/Dijkstra.h"
#include "DataStructures/Graph/Graph.h"
#include "DataStructures/Labels/Containers/StampedDistanceLabelContainer.h"
#include "DataStructures/Queues/Heap.h"

// Implementation of a CH query. Depending on the used label set, it keeps track of parent vertices
// and/or edges, and computes multiple shortest paths simultaneously, possibly using SSE or AVX
// instructions. The algorithm can be used with different distance label containers and queues.
template <
    typename CH, template <typename> class DistanceLabelContT, typename LabelSetT, typename QueueT,
    bool useStallOnDemand = true>
class CHQuery {
 private:
  using Graph = typename CH::SearchGraph; // The graph type on which we compute shortest paths.
  using Weight = typename CH::Weight;     // The attribute used as edge weight.

  static constexpr int K = LabelSetT::K; // The number of simultaneous shortest-path computations.

 public:
  // Constructs a CH search instance that uses stall-on-demand.
  template <bool cond = useStallOnDemand>
  CHQuery(std::enable_if_t<cond, const CH&> ch)
      : chSearch(ch.getUpwardGraph(), ch.getDownwardGraph(),
                 CHQueryPruningCriterion(ch.getDownwardGraph()),
                 CHQueryPruningCriterion(ch.getUpwardGraph())) {}

  // Constructs a CH search instance that does not use stall-on-demand.
  template <bool cond = useStallOnDemand>
  CHQuery(std::enable_if_t<!cond, const CH&> ch)
      : chSearch(ch.getUpwardGraph(), ch.getDownwardGraph()) {}

  // Ensures that the internal data structures fit for the size of the graph.
  void resize() {
    chSearch.resize();
  }

  // Runs a CH search from s to t.
  void run(const int s, const int t) {
    chSearch.run(s, t);
  }

  // Runs a CH search that computes multiple shortest paths simultaneously.
  void run(const std::array<int, K>& sources, const std::array<int, K>& targets) {
    chSearch.run(sources, targets);
  }

  // Returns the length of the i-th shortest path.
  int getDistance(const int i = 0) {
    return chSearch.getDistance(i);
  }

  // Returns the vertices on the i-th packed shortest path.
  std::vector<int> getPackedPath(const int i = 0) {
    return chSearch.getPath(i);
  }

  // Returns the edges on the i-th packed shortest path.
  std::vector<int> getPackedEdgePath(const int i = 0) {
    return chSearch.getEdgePath(i);
  }

 private:
  // The pruning criterion for a CH query, also known as stall-on-demand.
  struct CHQueryPruningCriterion {
    using LabelMask = typename LabelSetT::LabelMask;     // The label mask type we use.
    using DistLabel = typename LabelSetT::DistanceLabel; // The distance label type we use.

    // Constructs a pruning criterion for a CH query.
    CHQueryPruningCriterion(const Graph& oppositeGraph) : oppositeGraph(oppositeGraph) {}

    // Returns true if the search can be pruned at u.
    template <typename DistLabelContainerT>
    bool operator()(const int u, const DistLabel& distToU, DistLabelContainerT& distLabels) const {
      LabelMask continueSearch = true;
      FORALL_INCIDENT_EDGES(oppositeGraph, u, e)
        continueSearch &=
            distToU < distLabels[oppositeGraph.edgeHead(e)] + oppositeGraph.template get<Weight>(e);
      return !continueSearch;
    }

    const Graph& oppositeGraph; // The downward graph if we prune the upward search or vice versa.
  };

  // The stopping criterion for a CH query computing k shortest paths simultaneously. We can stop
  // the forward search as soon as mu_i <= Qf.minKey for all i = 1, ..., k. The reverse search is
  // stopped in the same way.
  template <typename>
  struct CHQueryStoppingCriterion {
    // Constructs a stopping criterion for a CH query.
    CHQueryStoppingCriterion(const QueueT& forwardQueue, const QueueT& reverseQueue,
                             const int& maxTentativeDistance)
        : forwardQueue(forwardQueue),
          reverseQueue(reverseQueue),
          maxTentativeDistance(maxTentativeDistance) {}

    // Returns true if we can stop the forward search.
    bool stopForwardSearch() const {
      return forwardQueue.empty() || maxTentativeDistance <= forwardQueue.minKey();
    }

    // Returns true if we can stop the reverse search.
    bool stopReverseSearch() const {
      return reverseQueue.empty() || maxTentativeDistance <= reverseQueue.minKey();
    }

    const QueueT& forwardQueue;      // The priority queue of the forward search.
    const QueueT& reverseQueue;      // The priority queue of the reverse search.
    const int& maxTentativeDistance; // The largest of all k tentative distances.
  };

  using CHDijkstra = std::conditional_t<
      useStallOnDemand,
      Dijkstra<Graph, Weight, DistanceLabelContT, LabelSetT, QueueT, CHQueryPruningCriterion>,
      Dijkstra<Graph, Weight, DistanceLabelContT, LabelSetT, QueueT>>;

  BiDijkstra<CHDijkstra, CHQueryStoppingCriterion> chSearch; // The modified bidirectional search.
};

// An alias template for a standard CH search.
template <typename CH, typename LabelSetT, bool useStallOnDemand = true>
using StandardCHQuery = CHQuery<
    CH, StampedDistanceLabelContainer, LabelSetT, QuadHeap, useStallOnDemand>;
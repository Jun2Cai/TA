#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <queue>
#include <vector>

#include "Algorithms/CCH/CCH.h"
#include "Algorithms/CCH/EliminationTreeQuery.h"
#include "Algorithms/CH/CH.h"
#include "DataStructures/Labels/BasicLabelSet.h"
#include "DataStructures/Labels/ParentInfo.h"

// This class implements a closest-POI query in a customizable contraction hierarchy, based on a
// traversal of its separator decomposition tree. It works in two phases. Given a set of POIs, the
// caller must first build a POI index by invoking the method buildPoiIndexFor(). The index is then
// used to run closest-POI queries by invoking the method findClosestPois().
class CCHClosestPoiQuery {
 public:
  // Precomputed auxiliary information to accelerate POI queries.
  struct PoiIndex {
    friend class CCHClosestPoiQuery;
   private:
    // Builds the POI index for the specified set of POI vertices.
    PoiIndex(const std::vector<int32_t>& pointsOfInterest, const int numVertices)
        : pointsOfInterest(pointsOfInterest),
          numPoisAmongFirst(numVertices + 1, pointsOfInterest.size()) {
      assert(!pointsOfInterest.empty());
      assert(std::is_sorted(pointsOfInterest.begin(), pointsOfInterest.end()));
      assert(numVertices > pointsOfInterest.back());
      for (auto i = 0, j = 0; i < pointsOfInterest.size(); ++i)
        for (; j <= pointsOfInterest[i]; ++j)
          numPoisAmongFirst[j] = i;
    }

    const std::vector<int32_t>& pointsOfInterest; // The POI vertices in increasing order of ID.
    std::vector<int32_t> numPoisAmongFirst;       // The number of POIs among the first i vertices.
  };

  // A point of interest returned by a POI query.
  struct Poi {
    // Returns true if this POI is closer to the source than the specified POI.
    bool operator<(const Poi& rhs) const noexcept {
      return dist < rhs.dist;
    }

    int vertex; // The vertex that contains this POI.
    int dist;   // The shortest-path distance from the source to this POI.
  };

  // Creates an instance of a closest-POI query in the specified hierarchy.
  CCHClosestPoiQuery(const CCH& cch, const CH& minWeightedCH)
      : cch(cch), elimTreeQuery(minWeightedCH, cch.getEliminationTree()) {
    assert(cch.getUpwardGraph().numVertices() == minWeightedCH.upwardGraph().numVertices());
  }

  // Builds the POI index for the specified set of POI vertices.
  PoiIndex buildPoiIndexFor(const std::vector<int32_t>& pointsOfInterest) const {
    return {pointsOfInterest, cch.getUpwardGraph().numVertices()};
  }

  // Returns the k closest POI vertices to s.
  const std::vector<Poi>& findClosestPois(const int s, const PoiIndex& idx, const int k = 1) {
    assert(idx.numPoisAmongFirst.size() == cch.getUpwardGraph().numVertices() + 1);
    const auto& decomp = cch.getSeparatorDecomposition();
    elimTreeQuery.pinForwardSearch(s);
    recursionStack.emplace_back();
    currentlyClosestPois.push({INVALID_VERTEX, INFTY});

    while (!recursionStack.empty()) {
      const auto subgraph = recursionStack.back();
      recursionStack.pop_back();

      if (subgraph.dist >= currentlyClosestPois.top().dist)
        continue;

      const auto firstPoi = idx.numPoisAmongFirst[subgraph.firstVertex];
      const auto lastPoi = idx.numPoisAmongFirst[decomp.lastSeparatorVertex(subgraph.node)];
      if (lastPoi - firstPoi <= RECURSION_THRESHOLD) {
        handleBaseCase(firstPoi, lastPoi, idx, k);
        continue;
      }

      const auto firstPoiInSep = idx.numPoisAmongFirst[decomp.firstSeparatorVertex(subgraph.node)];
      handleBaseCase(firstPoiInSep, lastPoi, idx, k);

      const auto stackSize = recursionStack.size();
      auto child = decomp.leftChild(subgraph.node), firstVertex = subgraph.firstVertex;
      while (child != 0) {
        const auto firstPoiInChild = idx.numPoisAmongFirst[firstVertex];
        const auto lastPoiInChild = idx.numPoisAmongFirst[decomp.lastSeparatorVertex(child)];
        if (lastPoiInChild - firstPoiInChild > 0) {
          const auto childContainsSrc = firstVertex <= s && s < decomp.lastSeparatorVertex(child);
          const auto dist = childContainsSrc ? 0 : computeDistToSubgraph(child);
          recursionStack.push_back({child, dist, firstVertex});
        }
        firstVertex = decomp.lastSeparatorVertex(child);
        child = decomp.rightSibling(child);
      }
      assert(firstVertex == decomp.firstSeparatorVertex(subgraph.node));
      std::sort(recursionStack.rbegin(), recursionStack.rend() - stackSize);
    }

    if (currentlyClosestPois.top().vertex == INVALID_VERTEX)
      currentlyClosestPois.pop();
    if (currentlyClosestPois.size() < closestPois.size())
      closestPois.erase(closestPois.begin() + currentlyClosestPois.size(), closestPois.end());
    else
      closestPois.insert(closestPois.end(), currentlyClosestPois.size() - closestPois.size(), {});
    for (auto iter = closestPois.rbegin(); iter < closestPois.rend(); ++iter) {
      *iter = currentlyClosestPois.top();
      currentlyClosestPois.pop();
    }
    assert(currentlyClosestPois.empty());
    assert(std::is_sorted(closestPois.begin(), closestPois.end()));
    return closestPois;
  }

 private:
  // A subgraph of the hierarchy that has not been examined yet.
  struct ActiveSubgraph {
    // Returns true if this subgraph is closer to the source than the specified subgraph.
    bool operator<(const ActiveSubgraph& rhs) const noexcept {
      return dist < rhs.dist;
    }

    int node;        // The separator decomposition node that corresponds to this subgraph.
    int dist;        // A lower bound on the distance to any vertex in this subgraph.
    int firstVertex; // The first vertex in this subgraph.
  };

  static constexpr int RECURSION_THRESHOLD = 8; // Used to stop the recursion during queries.

  // Checks for each specified POI whether it improves the farthest of the k currently closest POIs.
  void handleBaseCase(const int firstPoi, const int lastPoi, const PoiIndex& idx, const int k) {
    assert(firstPoi >= 0);
    assert(lastPoi <= idx.pointsOfInterest.size());
    for (auto i = firstPoi; i < lastPoi; ++i) {
      const auto poi = idx.pointsOfInterest[i];
      elimTreeQuery.runReverseSearch(poi);
      const auto distToPoi = elimTreeQuery.getDistance();
      if (currentlyClosestPois.size() < k) {
        currentlyClosestPois.push({poi, distToPoi});
      } else if (distToPoi < currentlyClosestPois.top().dist) {
        currentlyClosestPois.pop();
        currentlyClosestPois.push({poi, distToPoi});
      }
      assert(currentlyClosestPois.size() <= k);
    }
  }

  // Returns a lower bound on the distance from the source to any vertex in the specified subgraph.
  int computeDistToSubgraph(const int node) {
    const auto& upwardGraph = cch.getUpwardGraph();
    const auto highestVertex = cch.getSeparatorDecomposition().lastSeparatorVertex(node) - 1;
    const auto firstBoundaryVertex = &upwardGraph.edgeHead(upwardGraph.firstEdge(highestVertex));
    const auto lastBoundaryVertex = &upwardGraph.edgeHead(upwardGraph.lastEdge(highestVertex));
    elimTreeQuery.runReverseSearch(firstBoundaryVertex, lastBoundaryVertex);
    return elimTreeQuery.getDistance();
  }

  using ElimTreeQuery = EliminationTreeQuery<BasicLabelSet<0, ParentInfo::NO_PARENT_INFO>>;

  const CCH& cch;              // The metric-independent CCH.
  ElimTreeQuery elimTreeQuery; // A standard elimination tree query.

  std::vector<ActiveSubgraph> recursionStack;    // The stack of active subgraphs.
  std::priority_queue<Poi> currentlyClosestPois; // The closest POIs, with the farthest one on top.
  std::vector<Poi> closestPois;                  // The result of the last query, from close to far.
};

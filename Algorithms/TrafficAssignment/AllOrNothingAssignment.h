#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <vector>
#include "Algorithms/TrafficAssignment/MeasureBehavior/Measurebehavior.h"
#include "Algorithms/TrafficAssignment/Measurement/StaticalFunction.h"
#include "DataStructures/Utilities/OriginDestination.h"
#include "Stats/TrafficAssignment/AllOrNothingAssignmentStats.h"
#include "Tools/CommandLine/ProgressBar.h"
#include "Tools/Simd/AlignedVector.h"
#include "Tools/Timer.h"

// Implementation of an iterative all-or-nothing traffic assignment. Each OD pair is processed in
// turn and the corresponding OD flow (in our case always a single flow unit) is assigned to each
// edge on the shortest path between O and D. Other O-D paths are not assigned any flow. The
// procedure can be used with different shortest-path algorithms.
template <typename ShortestPathAlgoT>
class AllOrNothingAssignment {
 private:
  using InputGraph = typename ShortestPathAlgoT::InputGraph;

 public:
  // Constructs an all-or-nothing assignment instance.
  AllOrNothingAssignment(const InputGraph& graph,
                         const std::vector<ClusteredOriginDestination>& odPairs,
                         const bool verbose = true)
      : stats(odPairs.size()),
        shortestPathAlgo(graph),
        inputGraph(graph),
        odPairs(odPairs),
        verbose(verbose) {
    Timer timer;
    shortestPathAlgo.preprocess();
    stats.totalPreprocessingTime = timer.elapsed();
    stats.lastRoutingTime = stats.totalPreprocessingTime;
    stats.totalRoutingTime = stats.totalPreprocessingTime;
    if (verbose) std::cout << "  Prepro: " << stats.totalPreprocessingTime << "ms" << std::endl;
  }

  // Assigns all OD flows to their currently shortest paths.
  void run(const int skipInterval = 1, int origin = 0, int destination = 0, int output = 0,std::string fileName = "") {
    Timer timer;
    ++stats.numIterations;
    if (verbose) std::cout << "Iteration " << stats.numIterations << ": " << std::flush;
    shortestPathAlgo.customize();
    stats.lastCustomizationTime = timer.elapsed();
    std::map<std::string, std::vector<int32_t>> odPairPath;
    timer.restart();
    ProgressBar bar(std::ceil(1.0 * odPairs.size() / (K * skipInterval)), verbose);
    trafficFlows.assign(inputGraph.numEdges(), 0);
    stats.startIteration();
    auto totalNumPairsSampledBefore = 0;
    #pragma omp parallel
    {
      auto queryAlgo = shortestPathAlgo.getQueryAlgoInstance();
      auto checksum = int64_t{0};
      auto prevMinPathCost = int64_t{0};
      auto avgChange = 0.0;
      auto maxChange = 0.0;
      auto numPairsSampledBefore = 0;

      #pragma omp for schedule(dynamic, 64) nowait
      for (auto i = 0; i < odPairs.size(); i += K * skipInterval) {
        // Run multiple shortest-path computations simultaneously.
        std::array<int, K> sources;
        std::array<int, K> targets;
        sources.fill(odPairs[i].origin);
        targets.fill(odPairs[i].destination);
        auto k = 1;
        for (; k < K && i + k * skipInterval < odPairs.size(); ++k) {
          sources[k] = odPairs[i + k * skipInterval].origin;
          targets[k] = odPairs[i + k * skipInterval].destination;
        }
        queryAlgo.run(sources, targets, k, odPairPath);









        for (auto j = 0; j < k; ++j) {
          // Maintain the avg and max change in the OD distances between the last two iterations.
          const auto dst = odPairs[i + j * skipInterval].destination;
          const auto dist = queryAlgo.getDistance(dst, j);
          const auto prevDist = stats.lastDistances[i + j * skipInterval];
          const auto change = 1.0 * std::abs(dist - prevDist) / prevDist;
          numPairsSampledBefore += prevDist != -1;
          checksum += dist;
          prevMinPathCost += prevDist != -1 ? dist : 0;
          stats.lastDistances[i + j * skipInterval] = dist;
          avgChange += std::max(0.0, change);
          maxChange = std::max(maxChange, change);
        }
        ++bar;
      }

      #pragma omp critical (combineResults)
      {
        queryAlgo.addLocalToGlobalFlows();
        stats.lastChecksum += checksum;
        stats.prevMinPathCost += prevMinPathCost;
        stats.avgChangeInDistances += avgChange;
        stats.maxChangeInDistances = std::max(stats.maxChangeInDistances, maxChange);
        totalNumPairsSampledBefore += numPairsSampledBefore;
      }

    }
    bar.finish();

    shortestPathAlgo.propagateFlowsToInputEdges(trafficFlows); 
    std::for_each(trafficFlows.begin(), trafficFlows.end(), [&](int& f) { f *= skipInterval; });
    stats.lastQueryTime = timer.elapsed();
    stats.avgChangeInDistances /= totalNumPairsSampledBefore;
    stats.finishIteration();
    output++;
    output--;
    //test
/*    std::vector<int> test = {0,1,3,2,4,5,6,7,8,9,10};
    auto testValue = StaticalFunction::weightedQuantile(test, inputGraph, trafficFlows);
    for (int i = 0; i < 11; i++) {
        std::cout << "flow " << trafficFlows[i] << std::endl;
    }

      for (int i = 0; i < 11; i++) {
          std::cout << "time " << inputGraph.travelTime(i) << std::endl;
      }
    for(auto v : testValue) {
        std::cout << "quantia" << v << std::endl;
    }*/

     //junjunjun
     Measurebehavior measures;
     if (output == 1 ) {
             std::string iteration;
         if (stats.numIterations == 1) {
             iteration = "first";
         } else {
             iteration = "last";
         }
         auto anaFileName = fileName + "_measure";

         measures.measures(anaFileName, odPairPath,inputGraph, trafficFlows, stats.numIterations);


     }


    auto odFileName = fileName + "_" + std::to_string(origin) + "_" + std::to_string(destination) + "_" + std::to_string(stats.numIterations);
      outputOdPairPath(origin, destination,odFileName, odPairPath);






    if (verbose) {
      std::cout << " done.\n";
      std::cout << "  Checksum: " << stats.lastChecksum;
      std::cout << "  Custom: " << stats.lastCustomizationTime << "ms";
      std::cout << "  Queries: " << stats.lastQueryTime << "ms";
      std::cout << "  Routing: " << stats.lastRoutingTime << "ms\n";
      std::cout << std::flush;
    }





  }

/*
  void meanMeasure(std::string anaFileName,std::map<std::string, std::vector<int32_t>> odPairPath) {

          std::ofstream odAnaFile;
          odAnaFile.open(anaFileName+"_mean.csv");
          odAnaFile << "origin,destination,0_20,20_40,40_60,60_80,80_100,,0_20,20_40,40_60,60_80,80_100"<<std::endl;
          double sumOfSatInCount0 = 0;
          double sumOfSatInLength0 = 0;
          double sumOfSatInCount1 = 0;
          double sumOfSatInLength1 = 0;
          double sumOfSatInCount2 = 0;
          double sumOfSatInLength2 = 0;
          double sumOfSatInCount3 = 0;
          double sumOfSatInLength3 = 0;
          double sumOfSatInCount4 = 0;
          double sumOfSatInLength4 = 0;
          int sumOfMeaningfulPair = 0;

          for(auto iter = odPairPath.begin(); iter != odPairPath.end(); iter++) {
              auto pair = iter -> first;
              auto path = iter -> second;
              odAnaFile << pair << ",";

              int sat0 = 0;
              int sat1 = 0;
              int sat2 = 0;
              int sat3 = 0;
              int sat4 = 0;

              int totalLength = 0;
              int len0 = 0;
              int len1 = 0;
              int len2 = 0;
              int len3 = 0;
              int len4 = 0;



              for(auto e : path) {
                  if (path.size() > 0) {
                  auto vol = trafficFlows[e];
                  auto sat = (double)vol / (double)inputGraph.capacity(e);
                  if (sat >=0 && sat <0.2) {
                      sat0++;
                      len0 += inputGraph.length(e);
                      totalLength += inputGraph.length(e);

                  } else if (sat >= 0.2 && sat < 0.4) {
                      sat1++;
                      len1 += inputGraph.length(e);
                      totalLength += inputGraph.length(e);
                  } else if (sat >= 0.4 && sat < 0.6) {
                      sat2++;
                      len2 += inputGraph.length(e);
                      totalLength += inputGraph.length(e);
                  } else if (sat >= 0.6 && sat < 0.8) {
                      sat3++;
                      len3 += inputGraph.length(e);
                      totalLength += inputGraph.length(e);
                  } else if (sat > 0.8) {
                      sat4++;
                      len4 += inputGraph.length(e);
                      totalLength += inputGraph.length(e);
                  }
              }}



              double size = (double)path.size();
              double length = (double)totalLength ;


              if (size > 0 && length  > 0) {
                  sumOfMeaningfulPair++;
                  odAnaFile << sat0/size << "," << sat1/size << "," <<sat2/size << ","<<sat3/size << ","<<sat4/size << ",,"<<
                            len0/length << ","<< len1/length  << ","<< len2/length << ","<< len3/length << ","<< len4/length <<std::endl;
                  sumOfSatInCount0 += sat0/size;
                  sumOfSatInCount1 += sat1/size;
                  sumOfSatInCount2 += sat2/size;
                  sumOfSatInCount3 += sat3/size;
                  sumOfSatInCount4 += sat4/size;
                  sumOfSatInLength0 += len0/length;
                  sumOfSatInLength1 += len1/length;
                  sumOfSatInLength2 += len2/length;
                  sumOfSatInLength3 += len3/length;
                  sumOfSatInLength4 += len4/length;
              }


          }

          odAnaFile << "mean,," << sumOfSatInCount0/sumOfMeaningfulPair <<","<< sumOfSatInCount1/sumOfMeaningfulPair <<","<< sumOfSatInCount2/sumOfMeaningfulPair << "," <<  sumOfSatInCount3/sumOfMeaningfulPair << "," << sumOfSatInCount4/sumOfMeaningfulPair<<
                    ","<< "mean" << "," << sumOfSatInLength0/sumOfMeaningfulPair<<","<< sumOfSatInLength1/sumOfMeaningfulPair<<","<< sumOfSatInLength2/sumOfMeaningfulPair<<","<< sumOfSatInLength3/sumOfMeaningfulPair<<","<< sumOfSatInLength4/sumOfMeaningfulPair<<std::endl;

      }

    void longestSubPathMeasure(std::string anaFileName,std::map<std::string, std::vector<int32_t>> odPairPath)  {
        std::ofstream odAnaFile;
        odAnaFile.open(anaFileName+"_longestSubPath.csv");
        odAnaFile << "origin,destination,longest continuous sub-path"<<std::endl;
        for(auto iter = odPairPath.begin(); iter != odPairPath.end(); iter++) {
            auto pair = iter->first;
            auto path = iter->second;



            if (path.size() > 0) {
                odAnaFile << pair << ",";
                int longestCount= 0;
                int temp = 0;
                for(auto e : path) {
                    auto vol = trafficFlows[e];
                    auto sat = (double)vol / (double)inputGraph.capacity(e);

                    if (sat > 0.8) {
                        temp++;
                        if (temp > longestCount) {
                            longestCount = temp;
                        }
                    } else {
                        temp = 0;
                    }


                }


                odAnaFile  << longestCount << std::endl;
            }





        }
  }*/

  void outputOdPairPath(int origin, int destination, std::string fileName, std::map<std::string, std::vector<int32_t>> odPairPath) {
      auto key = std::to_string(origin)+","+std::to_string(destination);
      if (odPairPath.count(key) > 0) {
          std::vector<int32_t> path = odPairPath[key];
          if (path.size() > 0) {

              //outputs origin-destination path in json-file
              std::ofstream odFile;
              auto odFileName =
                      fileName + ".json";

              odFile.open(odFileName);
              odFile << "{\"type\" : \"FeatureCollection\" ," << std::endl;
              odFile << "\"features\" : [{" << std::endl;
              odFile << "\"type\" : \"Feature\"," << std::endl;
              odFile << "\"properties\" : {\"stroke\" : \"#" << "FF0000" << "\" } ," << std::endl;

              odFile << "\"geometry\" : {\"type\":\"MultiLineString\",\"coordinates\":" << std::endl;
              odFile << "[";
              odFile << "[[" << inputGraph.LatLngAttribute::latLng(origin).lngInDeg() << ","
                     << inputGraph.LatLngAttribute::latLng(origin).latInDeg() << "],";
              odFile <<"["<<inputGraph.LatLngAttribute::latLng(inputGraph.edgeHead(path.front())).lngInDeg() << ","
                     << inputGraph.LatLngAttribute::latLng(inputGraph.edgeHead(path.front())).latInDeg() << "]]";
              for (int i = 0; i < path.size() - 1; i++) {
                  int tail = inputGraph.edgeHead(path.at(i));
                  int head = inputGraph.edgeHead(path.at(i+1));
                  odFile <<", " << "[[" << inputGraph.LatLngAttribute::latLng(tail).lngInDeg() << ","
                         << inputGraph.LatLngAttribute::latLng(tail).latInDeg() << "], ";
                  odFile <<"[" << inputGraph.LatLngAttribute::latLng(head).lngInDeg() << ","
                         << inputGraph.LatLngAttribute::latLng(head).latInDeg() << "]] ";

              }

              odFile << "] }} , { ";
              //mark origin
              odFile << "\"type\" : \"Feature\"," << std::endl;
              odFile << "\"properties\" : {\"stroke\" : \"#" << "000000" << "\" } ," << std::endl;

              odFile << "\"geometry\" : {\"type\":\"Point\",\"coordinates\":" << std::endl;
              odFile << "["<< inputGraph.LatLngAttribute::latLng(origin).lngInDeg() << ","
                     << inputGraph.LatLngAttribute::latLng(origin).latInDeg() << "]";

              odFile << "}}]}";

          }
          }
      }







  // Returns the traffic flow on edge e.
  const int& trafficFlowOn(const int e) const {
    assert(e >= 0); assert(e < inputGraph.numEdges());
    return trafficFlows[e];
  }

  AllOrNothingAssignmentStats stats; // Statistics about the execution.

 private:
  // The maximum number of simultaneous shortest-path computations.
  static constexpr int K = ShortestPathAlgoT::K;

  using ODPairs = std::vector<ClusteredOriginDestination>;

  ShortestPathAlgoT shortestPathAlgo; // Algorithm computing shortest paths between OD pairs.
  const InputGraph& inputGraph;       // The input graph.
  const ODPairs& odPairs;             // The OD pairs to be assigned onto the graph.
  AlignedVector<int> trafficFlows;    // The traffic flows on the edges.
  const bool verbose;                 // Should informative messages be displayed?
};

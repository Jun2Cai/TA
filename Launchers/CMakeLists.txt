add_executable(AssignTraffic AssignTraffic.cc ../Algorithms/TrafficAssignment/Measurement/LongestSubPath.h ../Algorithms/TrafficAssignment/Measurement/QuantileDistribution.h)
target_compile_definitions(AssignTraffic PRIVATE CSV_IO_NO_THREAD)
target_link_libraries(AssignTraffic ${RoutingKit_LIBRARY})
if(OpenMP_FOUND)
  target_link_libraries(AssignTraffic OpenMP::OpenMP_CXX)
endif()

option(TA_USE_CFW "Make each descent direction conjugate to the last direction." ON)
option(TA_USE_SIMD_SEARCH "Use SIMD optimizations for the centralized shortest-path search." ON)
set(TA_LOG_K "" CACHE STRING "Choose the number of simultaneous shortest-path computations.")

if(NOT TA_USE_CFW)
  target_compile_definitions(AssignTraffic PRIVATE TA_NO_CFW)
endif()
if(NOT TA_USE_SIMD_SEARCH)
  target_compile_definitions(AssignTraffic PRIVATE TA_NO_SIMD_SEARCH)
endif()
if(NOT TA_LOG_K STREQUAL "")
  target_compile_definitions(AssignTraffic PRIVATE TA_LOG_K=${TA_LOG_K})
endif()

add_executable(
    CalculateDemand CalculateDemand.cc
    ../External/randomc/mersenne.cpp ../External/randomc/userintf.cpp ../External/stocc/stoc1.cpp)
target_compile_definitions(CalculateDemand PRIVATE CSV_IO_NO_THREAD)
target_link_libraries(CalculateDemand ${PROJ_LIBRARY} ${RoutingKit_LIBRARY})
if(OpenMP_FOUND)
  target_link_libraries(CalculateDemand OpenMP::OpenMP_CXX)
endif()

set(DC_MAX_NUM_SOURCES "" CACHE STRING "Choose the number of sources to be substituted in FRAD.")
if(NOT DC_MAX_NUM_SOURCES STREQUAL "")
  target_compile_definitions(CalculateDemand PRIVATE DC_MAX_NUM_SOURCES=${DC_MAX_NUM_SOURCES})
endif()

add_executable(RunP2PAlgo RunP2PAlgo.cc)
target_compile_definitions(RunP2PAlgo PRIVATE CSV_IO_NO_THREAD)
target_link_libraries(RunP2PAlgo ${RoutingKit_LIBRARY})
if(OpenMP_FOUND)
  target_link_libraries(RunP2PAlgo OpenMP::OpenMP_CXX)
endif()

option(USE_FAST_ELIMINATION_TREE_QUERY "Use the fast elimination tree query." ON)
if(NOT USE_FAST_ELIMINATION_TREE_QUERY)
  target_compile_definitions(RunP2PAlgo PRIVATE NO_FAST_ELIMINATION_TREE_QUERY)
endif()

add_executable(RunPoiAlgo RunPoiAlgo.cc)
target_link_libraries(RunPoiAlgo ${RoutingKit_LIBRARY})
if(OpenMP_FOUND)
  target_link_libraries(RunPoiAlgo OpenMP::OpenMP_CXX)
endif()

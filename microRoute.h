#pragma once
#include <stdlib.h>
#include <list>
#include <map>
#include <math.h>
namespace facebook {
namespace terragraph {
struct RouteIndex {
  uint8_t tx_idx; // transmit beamforming index of the SNR report
  uint8_t rx_idx; // receive beamforming index of the SNR report
  RouteIndex(const uint8_t tx_idx_temp, const uint8_t rx_idx_temp)
      : tx_idx(tx_idx_temp), rx_idx(rx_idx_temp) {}
  bool operator<(const RouteIndex& A) const {
    return tx_idx < A.tx_idx || (tx_idx == A.tx_idx && rx_idx < A.rx_idx);
  }
  bool operator==(const RouteIndex& A) const {
    return tx_idx == A.tx_idx && rx_idx == A.rx_idx;
  }
};

struct clusterReport {
  std::multimap<uint8_t, RouteIndex> txRxClusterMap;
  uint8_t bfPeakTx;
  uint8_t bfPeakRx;
};

typedef std::multimap<uint8_t, RouteIndex>::iterator multimap_iterator;
typedef std::multimap<uint8_t, RouteIndex>::reverse_iterator
    multimap_rev_iterator;
typedef std::map<RouteIndex, uint8_t>::iterator map_iterator;

class MicroRouteDetection {

 public:
  MicroRouteDetection(
      uint8_t bfSampling,
      uint8_t clusterMaxNum,
      uint8_t beamIdxMaxNum,
      const std::multimap<uint8_t, RouteIndex>& bfSnrReport,
      const std::map<RouteIndex, uint8_t>& bfTxRxSnrReport,
      const std::map<RouteIndex, uint8_t>& bfTxRxClusterIdxReport);

  // initial clustering to obtain micro-route candidate clusters
  void initialClustering();

  // load the best micro-route beamforming information from FW
  void bestMicroRouteFromFW(RouteIndex routeIndexFromFW);

  // micro-route identification
  std::multimap<uint8_t, RouteIndex>& microRouteDiscovery();

 private:
  // This evaluates whether a beamforming combination pair
  // establishes a new micro-route or not;
  // It drops the <tx,rx> pair which lies in the beamwidth
  // of the identified micro-route peaks
  bool validNewMicroRoute(
      const RouteIndex& routeIndex,
      uint8_t Snr, uint8_t MaxSnr);

  // This evaluates whether there are multiple <tx, rx> pairs
  // sharing the same maxSnr values
  bool equalPeakValue(uint8_t maxSnr);

  // This compares the <tx, rx> pairs with same Snr values and
  // find the optimal pair based on grid average
  RouteIndex compareEqualPeak(uint8_t maxSnr);

  // Add beam combination pair routeIndex with a delta value
  bool bfIndexSum(RouteIndex& routeIndex, const RouteIndex& delta);

  // Construct a cluster based on the cluster peak RouteIndexPeak as well as the
  // calculated beamwidth
  void clusterConstruction(
    std::multimap<uint8_t, RouteIndex> &tempBfSnrReport,
    const RouteIndex& routeIndexPeak, const uint8_t clusterIdx);

  // During clustering, if the beamforming pair routeIndex exists in previous
  // clusters,
  // it shall be removed from bfTxRxSnrReport/bfSnrReport
  void eliminateEntryFromClusterReport();

  // Computes the half beamwidth
  void halfBeamWidthInitialCalculation();

  // The following multi-map and map are used to
  // store input from E2E (firmware report).
  // This multimap is used for matrix clustering (data compression)
  std::multimap<uint8_t, RouteIndex> bfSnrReport_;

  // This map is used for matrix clustering (data compression)
  std::map<RouteIndex, uint8_t> bfTxRxSnrReport_;

  // This is used to record the cluster index
  // for each <tx, rx> beam combination,
  // where second element is the cluster index (set to 0 at initialization)
  std::map<RouteIndex, uint8_t> bfTxRxClusterIdxReport_;

  // Use to store a list of cluster reports
  std::list<clusterReport> bfTxRxClusters_;

  // this multimap is used for micro-route discovery
  std::multimap<uint8_t, RouteIndex> bfSnrReportMicroRoute_;

  // this multimap is used to store final micro-route indices
  std::multimap<uint8_t, RouteIndex> microRoute_;

  std::list<RouteIndex> gridDeltaValues_;

  // the min index distance between two beamforming combination pairs
  // 1 for initial scan: 1, 2, 3, 4... (64x64)
  // 2 for periodical scan: 2, 4, 6... (31x31)
  uint8_t bfSampling_;

  // max num of clusters
  uint8_t clusterMaxNum_;

  // grid size of the grid average, used to evaluate equal peaks in the
  // compareEqualPeak function
  uint8_t gridSize_;

  // beam index max value = 63
  uint8_t beamIdxMaxNum_;

  // to calculate the beam index separation for limited/full nulling,
  // map input beam index x (in the set: 31,30,...,1,0,32,33,...,62,63) to beam
  // index y (in the set: 0,1,...,31,32,33,...,62,63). beamIndexMapping[x] = y
  uint8_t beamIndexMapping_[64];

  // half beamwidth values for all beamforming pairs are pre-calculated
  // to avoid repeating expensive calculations
  uint8_t halfBeamWidth_[64];
  
  // this flag indicates whether we obtain the best 
  // micro-route beamforming pair from FW report or not
  bool bestRouteFromFWFlag_;

  // best micro-route beamforming pair from FW
  RouteIndex bestRouteIndexFromFW_;
};

} // namespace terragraph
} // namespace facebook

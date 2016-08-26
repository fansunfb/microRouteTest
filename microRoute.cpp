#include "microRoute.h"
//using namespace std;
namespace facebook {
namespace terragraph {

MicroRouteDetection::MicroRouteDetection(
    uint8_t bfSampling,
    uint8_t clusterMaxNum,
    uint8_t beamIdxMaxNum,
    const std::multimap<uint8_t, RouteIndex>& bfSnrFwReport,
    const std::map<RouteIndex, uint8_t>& bfTxRxSnrFwReport,
    const std::map<RouteIndex, uint8_t>& bfTxRxClusterIdxReport)
    : gridSize_(9),
      bfSampling_(bfSampling),
      clusterMaxNum_(clusterMaxNum),
      beamIdxMaxNum_(beamIdxMaxNum),
      bestRouteFromFWFlag_(false),
      bestRouteIndexFromFW_(-1,-1),
      bfSnrReport_(bfSnrFwReport),
      bfSnrReportMicroRoute_(bfSnrFwReport),
      bfTxRxSnrReport_(bfTxRxSnrFwReport),
      bfTxRxClusterIdxReport_(bfTxRxClusterIdxReport) {
  // build 3x3 grid to compare the <tx, rx> pairs with same Snr values
  RouteIndex RouteIndexTemp1(-bfSampling_, -bfSampling_);
  RouteIndex RouteIndexTemp2(-bfSampling_, 0);
  RouteIndex RouteIndexTemp3(-bfSampling_, bfSampling_);
  RouteIndex RouteIndexTemp4(0, -bfSampling_);
  RouteIndex RouteIndexTemp5(0, bfSampling_);
  RouteIndex RouteIndexTemp6(bfSampling_, -bfSampling_);
  RouteIndex RouteIndexTemp7(bfSampling_, 0);
  RouteIndex RouteIndexTemp8(bfSampling_, bfSampling_);

  // any better way to initialize the following std::list?
  gridDeltaValues_.push_back(RouteIndexTemp1);
  gridDeltaValues_.push_back(RouteIndexTemp2);
  gridDeltaValues_.push_back(RouteIndexTemp3);
  gridDeltaValues_.push_back(RouteIndexTemp4);
  gridDeltaValues_.push_back(RouteIndexTemp5);
  gridDeltaValues_.push_back(RouteIndexTemp6);
  gridDeltaValues_.push_back(RouteIndexTemp7);
  gridDeltaValues_.push_back(RouteIndexTemp8);

  for (uint8_t tempIdx = 0; tempIdx <= beamIdxMaxNum_; tempIdx++) {
    if (tempIdx <= 31)
      beamIndexMapping_[tempIdx] = (31 - tempIdx);
    else
      beamIndexMapping_[tempIdx] = tempIdx;
  }
  halfBeamWidthInitialCalculation();
}

void
MicroRouteDetection::initialClustering() {
  uint8_t clusterIdx = 1;

  // When obtain a micro-route candidate, clusterIdx increases, the peak/center
  // of the cluster will be removed from bfSnrReport_ and bfTxRxSnrReport_
  while ((bfSnrReport_.size() != 0) && (clusterIdx <= clusterMaxNum_)) {
    uint8_t maxSnr = 0;
    RouteIndex routeIndexSelectedPeak(-1, -1);
    std::multimap<uint8_t, RouteIndex> tempBfSnrReport;
    clusterReport txRxSnrCluster;

    // if the potential peak for an additional cluster lie in previous clusters,
    // remove it
    if (bfTxRxClusters_.size() != 0) {
      eliminateEntryFromClusterReport();
    }
    if (bfSnrReport_.empty())
      break;

    if (equalPeakValue(maxSnr)) {
      uint8_t peakNum = bfSnrReport_.count(maxSnr);

      // check the 3x3 grid via search over bfTxRxSnrReport_ map
      routeIndexSelectedPeak = compareEqualPeak(maxSnr);

      std::pair<multimap_iterator, multimap_iterator> iterPair =
          bfSnrReport_.equal_range(maxSnr);

      for (multimap_iterator it = iterPair.first; it != iterPair.second; ++it) {
        if (it->second == routeIndexSelectedPeak) {
          if (bfTxRxClusterIdxReport_.at(routeIndexSelectedPeak) == 0) {
            tempBfSnrReport.insert(*it);
            bfTxRxClusterIdxReport_.at(routeIndexSelectedPeak) = clusterIdx;
          }
          bfTxRxSnrReport_.erase(it->second);
          bfSnrReport_.erase(it);
          break;
        }
      }
    } else {
      multimap_rev_iterator rit = bfSnrReport_.rbegin();
      routeIndexSelectedPeak = rit->second;
      tempBfSnrReport.insert(*(rit));
      bfTxRxClusterIdxReport_.at(routeIndexSelectedPeak) = clusterIdx;
      bfTxRxSnrReport_.erase(rit->second);
      bfSnrReport_.erase(--(rit.base()));
    }
    clusterConstruction(tempBfSnrReport, routeIndexSelectedPeak, clusterIdx);
    txRxSnrCluster.txRxClusterMap = tempBfSnrReport;
    txRxSnrCluster.bfPeakTx = routeIndexSelectedPeak.tx_idx;
    txRxSnrCluster.bfPeakRx = routeIndexSelectedPeak.rx_idx;
    bfTxRxClusters_.push_back(txRxSnrCluster);
    clusterIdx++;
  }
}

void
MicroRouteDetection::bestMicroRouteFromFW(
    bool microRouteFromFWFlag,
    RouteIndex routeIndexFromFW) {
  bestRouteFromFWFlag_ = microRouteFromFWFlag;
  bestRouteIndexFromFW_ = routeIndexFromFW;
}

std::multimap<uint8_t, RouteIndex>&
MicroRouteDetection::microRouteDiscovery() {
  // if equal peaks exist, maybe we shall consider equalPeakValue(maxSnr) using
  // compareEqualPeak(maxSnr)
  uint8_t maxSnr = bfSnrReportMicroRoute.rbegin()->first;
  for (multimap_rev_iterator rit = bfSnrReportMicroRoute.rbegin();
      rit != bfSnrReportMicroRoute.rend(); ) {
    uint8_t snrTemp = rit->first;
    RouteIndex routeIndexTemp(-1,-1);
    if (!bestRouteFromFWFlag)
      routeIndexTemp = rit->second;
    else
      routeIndexTemp = bestRouteIndexFromFW;

    if (microRoute.empty()) {
      microRoute.insert(
          std::make_pair<uint8_t, RouteIndex>(snrTemp, routeIndexTemp));
      bfSnrReportMicroRoute.erase(--(rit.base()));
    } else {
      if (validNewMicroRoute(routeIndexTemp, snrTemp, maxSnr)) {
        microRoute.insert(
            std::make_pair<uint8_t, RouteIndex>(snrTemp, routeIndexTemp));
        bfSnrReportMicroRoute.erase(--(rit.base()));
      } else {
        bfSnrReportMicroRoute.erase(--(rit.base()));
      }
    }
  }
  return microRoute;
}

bool
MicroRouteDetection::validNewMicroRoute(
    const RouteIndex& routeIndex,
    uint8_t Snr, uint8_t MaxSnr) {
  bool validFlag = false;
  uint8_t snrThreshold = 2;
  for (multimap_iterator it = microRoute.begin(); it != microRoute.end();
      ++it) {
    // drop the <tx,rx> pair which lies
    // in the beamwidth of an identified micro-route
    uint8_t txHalfBeamwidth = halfBeamWidth[it->second.tx_idx];
    uint8_t rxHalfBeamwidth = halfBeamWidth[it->second.rx_idx];
    uint8_t mappedTxIdxRef = beamIndexMapping[it->second.tx_idx];
    uint8_t mappedRxIdxRef = beamIndexMapping[it->second.rx_idx];
    uint8_t mappedTxIdx = beamIndexMapping[routeIndex.tx_idx];
    uint8_t mappedRxIdx = beamIndexMapping[routeIndex.rx_idx];
    uint8_t txSeperation = abs(mappedTxIdx - mappedTxIdxRef);
    uint8_t rxSeperation = abs(mappedRxIdx - mappedRxIdxRef);
    if ((txSeperation <= txHalfBeamwidth) &&
        (rxSeperation <= rxHalfBeamwidth)) {
      validFlag = false;
      break;
    }

    if ((txSeperation <= txHalfBeamwidth) ||
        (rxSeperation <= rxHalfBeamwidth)) {
      if (Snr < (int)(MaxSnr - snrThreshold)) {
        validFlag = false;
        break;
      } else
        validFlag = true;
    }
    if ((txSeperation >= txHalfBeamwidth) && (rxSeperation >= rxHalfBeamwidth))
      validFlag = true;
  }
  return validFlag;
}

bool
MicroRouteDetection::equalPeakValue(uint8_t maxSnr) {
  multimap_rev_iterator rit;
  rit = bfSnrReport_.rbegin();
  maxSnr = rit->first;
  if (bfSnrReport_.count(maxSnr) > 1) {
    return true;
  } else
    return false;
}

RouteIndex
MicroRouteDetection::compareEqualPeak(uint8_t maxSnr) {
  std::pair<multimap_iterator, multimap_iterator> iterPair =
      bfSnrReport_.equal_range(maxSnr);
  RouteIndex routeIndexTemp(-1, -1);
  RouteIndex routeIndexSelectedPeak(-1, -1);
  uint16_t snrSum;
  uint8_t validNeighborNum;
  uint8_t averageSnr = 0;
  uint8_t averageSnrSelectedPeak = 0;
  uint8_t indexTemp = 0;

  for (multimap_iterator it = iterPair.first; it != iterPair.second;) {
    if (bfTxRxClusterIdxReport.at(it->second) > 0) {
      bfTxRxSnrReport_.erase(it->second);
      bfSnrReport_.erase(it++);
    } else
      ++it;
  }
  iterPair = bfSnrReport_.equal_range(maxSnr);
  routeIndexSelectedPeak = iterPair.first->second;
  for (multimap_iterator it = iterPair.first; it != iterPair.second; ++it) {
    if (bfTxRxSnrReport_.find(it->second) != bfTxRxSnrReport_.end()) {
      snrSum = 0;
      validNeighborNum = 0;
      // check the nearby 3x3 grid based on bfSampling_
      for (std::list<RouteIndex>::iterator iter = gridDeltaValues_.begin();
          iter != gridDeltaValues_.end(); ++iter) {
        routeIndexTemp = it->second;
        if (bfIndexSum(routeIndexTemp, *iter)) {
          map_iterator mit = bfTxRxSnrReport_.find(routeIndexTemp);
          if (mit != bfTxRxSnrReport_.end()) {
            validNeighborNum++;
            snrSum += mit->second;
          }
        }
      }
      if (validNeighborNum > 0)
        averageSnr = snrSum / validNeighborNum;
      indexTemp++;
    }
    // update RouteIndexSelectedPeak
    if (averageSnr > averageSnrSelectedPeak) {
      routeIndexSelectedPeak = it->second;
      averageSnrSelectedPeak = averageSnr;
    }
  }
  return routeIndexSelectedPeak;
}

bool
MicroRouteDetection::bfIndexSum(
    RouteIndex& routeIndex, const RouteIndex& delta) {
  // this might be negative as delta can be delta can be negative
  int8_t tx_idx_temp = routeIndex.tx_idx + delta.tx_idx;
  int8_t rx_idx_temp = routeIndex.rx_idx + delta.rx_idx;
  if ((tx_idx_temp > beamIdxMaxNum_) || (rx_idx_temp > beamIdxMaxNum_)
      || (tx_idx_temp < 0) || (rx_idx_temp < 0))
    return false;
  else {
    routeIndex.tx_idx = tx_idx_temp;
    routeIndex.rx_idx = rx_idx_temp;
    return true;
  }
}

void
MicroRouteDetection::clusterConstruction(
    std::multimap<uint8_t, RouteIndex>& tempBfSnrReport,
    const RouteIndex& RouteIndexPeak, const uint8_t clusterIdx) {
  uint8_t temp = 0;
  for (multimap_iterator it = bfSnrReport_.begin(); it != bfSnrReport_.end();
      ++it) {
    uint8_t txHalfBeamwidth = halfBeamWidth[RouteIndexPeak.tx_idx];
    uint8_t rxHalfBeamwidth = halfBeamWidth[RouteIndexPeak.rx_idx];
    uint8_t mappedTxIdxPeak = beamIndexMapping[RouteIndexPeak.tx_idx];
    uint8_t mappedRxIdxPeak = beamIndexMapping[RouteIndexPeak.rx_idx];
    uint8_t mappedTxIdx = beamIndexMapping[it->second.tx_idx];
    uint8_t mappedRxIdx = beamIndexMapping[it->second.rx_idx];

    uint8_t txSeperation = abs(mappedTxIdxPeak - mappedTxIdx);
    uint8_t rxSeperation = abs(mappedRxIdxPeak - mappedRxIdx);
    temp++;
    // check BeamWidth
    if ((txSeperation <= txHalfBeamwidth) &&
        (rxSeperation <= rxHalfBeamwidth)) {
      tempBfSnrReport.insert(*it);
      bfTxRxClusterIdxReport_.at(it->second) = clusterIdx;
    }
  }
}

void
MicroRouteDetection::eliminateEntryFromClusterReport() {
  for (multimap_rev_iterator rit = bfSnrReport_.rbegin();
       rit != bfSnrReport_.rend();) {
    RouteIndex routeIndex = rit->second;

    // if the routeIndex exists in clusters, remove it from
    // bfTxRxSnrReport_/bfSnrReport_
    if (bfTxRxClusterIdxReport_.at(routeIndex) > 0) {
      bfTxRxSnrReport_.erase(routeIndex);
      bfSnrReport_.erase(--(rit.base()));
    } else {
      break;
    }
  }
}

void
MicroRouteDetection::halfBeamWidthInitialCalculation() {
  double cbResolution = 90 / (beamIdxMaxNum_+1);
  uint16_t broadSideBeamWidth = 8;
  double pi = 3.1415;
  uint8_t bfIndex;
  for (bfIndex=0; bfIndex <= beamIdxMaxNum_; ++bfIndex)
  {
    if (bfIndex > 31)
      halfBeamWidth_[bfIndex] = round(0.5 * (1 / cbResolution)
          * broadSideBeamWidth
          / (cos((pi / 180) * abs(bfIndex - 32) * cbResolution)));
    else
      halfBeamWidth_[bfIndex] = round(0.5 * (1 / cbResolution)
          * broadSideBeamWidth
          / (cos((pi / 180) * bfIndex * cbResolution)));
  }
}

} // namespace terragraph
} // namespace facebook

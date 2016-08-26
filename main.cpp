#include <fstream>
//#include <cstdlib>
#include <sstream>
#include <iostream>
//using std::make_pair;
#include "microRoute.h"
using namespace facebook::terragraph; //using namespace std;
int main() {
  std::multimap<uint8_t, RouteIndex> bfSnrReport;
  std::map<RouteIndex, uint8_t> bfTxRxSnrReport; // needs to define the compare for the key (struct RouteIndex)
  std::map<RouteIndex, uint8_t> bfTxRxClusterInitialIdxReport;
  uint8_t bfSampling = 1;
  // parse from a firmware measurement based csv file
  std::ifstream file("parsed_chanSurveyFW_syslog_L8.csv");
  std::string csvLine;
  std::getline(file, csvLine); //pass the first row with only words, i.e. tx, rx, SNR, ...
  while (std::getline(file, csvLine))
  {
    std::istringstream csvStream(csvLine);
    std::string csvElement;
    std::getline(csvStream, csvElement, ',');
    uint8_t tx_idx_temp = atoi(csvElement.c_str());
    std::getline(csvStream, csvElement, ',');
    uint8_t rx_idx_temp = atoi(csvElement.c_str());
    std::getline(csvStream, csvElement, ',');
    if (csvElement != "nan") {
      // this snr has to be int; uint will not work
      int8_t snr_temp = atoi(csvElement.c_str());
      if (snr_temp > 0) {
        RouteIndex routeIndex(tx_idx_temp, rx_idx_temp);
        bfSnrReport.insert(std::make_pair<uint8_t, RouteIndex>(snr_temp, routeIndex));    
      }
    }
  }
  uint8_t clusterMaxNum = 8; // max num of candidate clusters in the initial clustering stage 
  uint8_t snr_threshold = 0; // this threshold is used to perform pre-filtering, and it is based on SNR for MCS 3
  uint8_t beamIdxMaxNum = 63;
  for(multimap_iterator it = bfSnrReport.begin(); it != bfSnrReport.end(); ) {
    uint8_t snrTemp = it->first;
    RouteIndex routeIndex = it->second;
    // pre-filtering based on snr_threshold
    if (snrTemp > snr_threshold) {
      bfTxRxSnrReport.insert(std::make_pair<RouteIndex, int8_t>(routeIndex, snrTemp));
      bfTxRxClusterInitialIdxReport.insert(std::make_pair<RouteIndex, int8_t>(routeIndex, 0)); // cluster index of each tx, rx pair
      it ++;
    }
    else {
      bfSnrReport.erase(it ++); //important step
    }
  }
  std::cout << "After pre-filtering, size of bfSnrReport=" << bfSnrReport.size() << ",size of bfTxRxSnrReport=" << bfTxRxSnrReport.size() << std::endl;

  MicroRouteDetection microRoute(bfSampling, clusterMaxNum, beamIdxMaxNum, 
      bfSnrReport, bfTxRxSnrReport, bfTxRxClusterInitialIdxReport);
  microRoute.initialClustering();

  std::multimap<uint8_t, RouteIndex> microRoutes;

  RouteIndex bestMicroRouteFW(19,49);
  // For micro-route discovery, insert best micro-route beamforming pair from FW report
  microRoute.bestMicroRouteFromFW(bestMicroRouteFW);
  microRoutes = microRoute.microRouteDiscovery();

  std::cout << "" << std::endl;
  std::cout << "Final microRoute size=" << microRoutes.size() << std::endl;
  for(multimap_iterator it = microRoutes.begin(); it != microRoutes.end(); ) {
    std::cout << "peak tx=" << (int)it->second.tx_idx << ",peak rx=" << (int)it->second.rx_idx << ",SNR=" << (int)it->first << std::endl;
    it ++;
  }
}
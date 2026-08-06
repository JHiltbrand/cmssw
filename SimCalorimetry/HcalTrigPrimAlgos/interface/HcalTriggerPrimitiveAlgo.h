#ifndef HcalSimAlgos_HcalTriggerPrimitiveAlgo_h
#define HcalSimAlgos_HcalTriggerPrimitiveAlgo_h

#include "CalibCalorimetry/HcalTPGAlgos/interface/HcaluLUTTPGCoder.h"
#include "CalibFormats/CaloTPG/interface/HcalTPGCompressor.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CondFormats/HcalObjects/interface/HcalElectronicsMap.h"
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "SimCalorimetry/HcalTrigPrimAlgos/interface/HcalFeatureHFEMBit.h"
#include "SimCalorimetry/HcalTrigPrimAlgos/interface/HcalFinegrainBit.h"

#include <map>
#include <vector>

class CaloGeometry;
class IntegerCaloSamples;

class HcalTriggerPrimitiveAlgo {
public:
  HcalTriggerPrimitiveAlgo(const std::vector<uint32_t>& FG_HF_thresholds,
                           bool useTDCfromDB,
                           int numberOfSamples,
                           int numberOfPresamples,
                           int numberOfFilterPresamplesHBQIE11,
                           int numberOfSamplesHF,
                           int numberOfPresamplesHF,
                           bool useTDCInMinBiasBits);
  ~HcalTriggerPrimitiveAlgo();

  template <typename... Digis>
  void run(const HcalTPGCoder* incoder,
           const HcalTPGCompressor* outcoder,
           const HcalDbService* conditions,
           HcalTrigPrimDigiCollection& result,
           const HcalTrigTowerGeometry* trigTowerGeometry,
           float rctlsb,
           const HcalFeatureBit* LongvrsShortCut,
           const Digis&... digis);

  template <typename T, typename... Args>
  void addDigis(const T& collection, const Args&... digis) {
    addDigis(collection);
    addDigis(digis...);
  };

  template <typename T>
  void addDigis(const T& collection) {
    for (const auto& digi : collection) {
      addSignal(digi);
    }
  };

  template <typename D>
  void addDigis(const HcalDataFrameContainer<D>& collection) {
    for (auto i = collection.begin(); i != collection.end(); ++i) {
      D digi(*i);
      addSignal(digi);
    }
  };

  void setWeightsQIE11(const edm::ParameterSet& weightsQIE11);
  void setWeightQIE11(int aieta, int weight);
  void setCodedVetoThresholds(const edm::ParameterSet& codedVetoThresholds);
  void setCodedVetoThreshold(int aieta, int codedVetoThreshold);
  void setNCTScaleShift(int);
  void setRCTScaleShift(int);

  void setNumFilterPresamplesHBQIE11(int presamples) { numberOfFilterPresamplesHBQIE11_ = presamples; }

  void setFixSaturationFlag(bool fix_saturation);
  void overrideParameters(const edm::ParameterSet& ps);

private:
  /// adds the signal to the map
  void addSignal(const QIE10DataFrame& frame);
  void addSignal(const QIE11DataFrame& frame);
  void addSignal(const IntegerCaloSamples& samples);
  void addUpgradeFG(const HcalTrigTowerDetId& id, int depth, const std::vector<std::bitset<2>>& bits);
  void addUpgradeTDCFG(const HcalTrigTowerDetId& id, const QIE11DataFrame& frame);

  bool passTDC(const QIE10DataFrame& digi, int ts) const;
  bool validChannel(const QIE10DataFrame& digi, int ts) const;

  // 2017 and later: QIE11
  void analyzeQIE11(IntegerCaloSamples& samples,
                    std::vector<bool> sample_saturation,
                    HcalTriggerPrimitiveDigi& result,
                    const HcalFinegrainBit& fg_algo);
  // With dual anode readout
  void analyzeHFQIE10(const IntegerCaloSamples& SAMPLES,
                      HcalTriggerPrimitiveDigi& result,
                      const int HF_LUMI_SHIFT,
                      const HcalFeatureBit* HCALFEM);

  void analyzeZDC(IntegerCaloSamples& samples, HcalTriggerPrimitiveDigi& result);

  // Member initialized by constructor
  const HcaluLUTTPGCoder* incoder_;
  const HcalTPGCompressor* outcoder_;
  const HcalDbService* conditions_;
  double theThreshold;
  std::array<std::array<int, 2>, 17> weightsQIE11_;
  std::array<int, 17> codedVetoThresholds_;
  std::vector<uint32_t> FG_HF_thresholds_;
  bool useTDCfromDB_;
  int numberOfSamples_;
  int numberOfPresamples_;
  int numberOfFilterPresamplesHBQIE11_;
  int numberOfSamplesHF_;
  int numberOfPresamplesHF_;
  bool useTDCInMinBiasBits_;
  int NCTScaleShift;
  int RCTScaleShift;

  const HcalTrigTowerGeometry* theTrigTowerGeometry;

  typedef std::map<HcalTrigTowerDetId, IntegerCaloSamples> SumMap;
  SumMap theSumMap;

  typedef std::map<HcalTrigTowerDetId, std::vector<bool>> SatMap;
  SatMap theSatMap;

  struct HFUpgradeDetails {
    IntegerCaloSamples samples;
    QIE10DataFrame digi;
    std::vector<bool> validity;
    std::vector<std::bitset<2>> fgbits;
    std::vector<bool> passTDC;
  };
  typedef std::map<HcalTrigTowerDetId, std::map<uint32_t, std::array<HFUpgradeDetails, 4>>> HFUpgradeDetailMap;
  HFUpgradeDetailMap theHFUpgradeDetailMap;

  typedef std::vector<IntegerCaloSamples> SumFGContainer;
  typedef std::map<HcalTrigTowerDetId, SumFGContainer> TowerMapFGSum;
  TowerMapFGSum theTowerMapFGSum;

  // ==============================
  // =  HF Veto
  // ==============================
  // Sum = Long + Short;" // intermediate calculation.
  //  if ((Short < MinSignalThresholdET OR Long  < MinSignalThresholdET)
  //     AND Sum > PMTNoiseThresholdET) VetoedSum = 0;
  //  else VetoedSum = Sum;
  // ==============================
  // Map from FG id to veto booleans
  HcalFeatureBit* LongvrsShortCut;

  typedef std::map<HcalTrigTowerDetId, std::vector<bool>> FGbitMap;
  FGbitMap fgMap_;

  typedef std::vector<HcalFinegrainBit::Tower> FGUpgradeContainer;
  typedef std::map<HcalTrigTowerDetId, FGUpgradeContainer> FGUpgradeMap;
  FGUpgradeMap fgUpgradeMap_;

  typedef std::vector<HcalFinegrainBit::TowerTDC> FGUpgradeTDCContainer;
  typedef std::map<HcalTrigTowerDetId, FGUpgradeTDCContainer> FGUpgradeTDCMap;
  FGUpgradeTDCMap fgUpgradeTDCMap_;

  bool fix_saturation_ = false;

  edm::ParameterSet override_parameters_;

  bool override_adc_hf_ = false;
  uint32_t override_adc_hf_value_;
  bool override_tdc_hf_ = false;
  unsigned long long override_tdc_hf_value_;

  // Maximum valid TDC value for setting timing bits
  static constexpr int tdcmax_ = 49;

  // Fine-grain in HF ignores tower 29, and starts with 30
  static const int FIRST_FINEGRAIN_TOWER = 30;

  static const int QIE10_LINEARIZATION_ET = HcaluLUTTPGCoder::QIE10_LUT_BITMASK;
  static const int QIE11_LINEARIZATION_ET = HcaluLUTTPGCoder::QIE11_LUT_BITMASK;
  // Consider CaloTPGTranscoderULUT.h for values
  static const int QIE10_MAX_LINEARIZATION_ET = 0x7FF;
  static const int QIE11_MAX_LINEARIZATION_ET = 0x7FF;
  static const int QIE10_ZDC_MAX_LINEARIZATION_ET = 0x3FF;
};

template <typename... Digis>
void HcalTriggerPrimitiveAlgo::run(const HcalTPGCoder* incoder,
                                   const HcalTPGCompressor* outcoder,
                                   const HcalDbService* conditions,
                                   HcalTrigPrimDigiCollection& result,
                                   const HcalTrigTowerGeometry* trigTowerGeometry,
                                   float rctlsb,
                                   const HcalFeatureBit* LongvrsShortCut,
                                   const Digis&... digis) {
  theTrigTowerGeometry = trigTowerGeometry;

  incoder_ = dynamic_cast<const HcaluLUTTPGCoder*>(incoder);
  outcoder_ = outcoder;
  conditions_ = conditions;

  theSumMap.clear();
  theSatMap.clear();
  theTowerMapFGSum.clear();
  fgUpgradeMap_.clear();
  fgUpgradeTDCMap_.clear();
  theHFUpgradeDetailMap.clear();

  // Add all digi collections
  addDigis(digis...);

  // Prepare the fine-grain calculation algorithm for HB
  int version = 0;
  version = conditions_->getHcalTPParameters()->getFGVersionHBHE();
  HcalFinegrainBit fg_algo(version);

  for (auto& item : theSumMap) {
    result.push_back(HcalTriggerPrimitiveDigi(item.first));
    HcalTrigTowerDetId detId(item.second.id());
    if (detId.ietaAbs() >= theTrigTowerGeometry->firstHFTower(detId.version()) && detId.ietaAbs() < 42) {
      analyzeHFQIE10(item.second, result.back(), NCTScaleShift, LongvrsShortCut);
    } else if (detId.ietaAbs() >= theTrigTowerGeometry->firstHFTower(detId.version()) && detId.ietaAbs() > 41) {
      analyzeZDC(item.second, result.back());
    } else {
      SatMap::iterator item_sat = theSatMap.find(detId);
      if (item_sat == theSatMap.end())
        analyzeQIE11(item.second, std::vector<bool>(), result.back(), fg_algo);
      else
        analyzeQIE11(item.second, item_sat->second, result.back(), fg_algo);
    }
  }

  // Free up some memory
  theSumMap.clear();
  theSatMap.clear();
  theTowerMapFGSum.clear();
  fgUpgradeMap_.clear();
  fgUpgradeTDCMap_.clear();
  theHFUpgradeDetailMap.clear();

  return;
}

#endif

#include "CalibFormats/CaloObjects/interface/IntegerCaloSamples.h"
#include "CondFormats/HcalObjects/interface/HcalTPParameters.h"
#include "CondFormats/HcalObjects/interface/HcalTPChannelParameters.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "SimCalorimetry/HcalTrigPrimAlgos/interface/HcalTriggerPrimitiveAlgo.h"

#include <iostream>

HcalTriggerPrimitiveAlgo::HcalTriggerPrimitiveAlgo(const std::vector<uint32_t>& FG_HF_thresholds,
                                                   bool useTDCfromDB,
                                                   int numberOfSamples,
                                                   int numberOfPresamples,
                                                   int numberOfFilterPresamplesHBQIE11,
                                                   int numberOfSamplesHF,
                                                   int numberOfPresamplesHF,
                                                   bool useTDCInMinBiasBits)
    : incoder_(nullptr),
      outcoder_(nullptr),
      theThreshold(0),
      FG_HF_thresholds_(FG_HF_thresholds),
      useTDCfromDB_(useTDCfromDB),
      numberOfSamples_(numberOfSamples),
      numberOfPresamples_(numberOfPresamples),
      numberOfFilterPresamplesHBQIE11_(numberOfFilterPresamplesHBQIE11),
      numberOfSamplesHF_(numberOfSamplesHF),
      numberOfPresamplesHF_(numberOfPresamplesHF),
      useTDCInMinBiasBits_(useTDCInMinBiasBits),
      NCTScaleShift(0),
      RCTScaleShift(0),
      override_parameters_() {}

HcalTriggerPrimitiveAlgo::~HcalTriggerPrimitiveAlgo() {}

void HcalTriggerPrimitiveAlgo::setFixSaturationFlag(bool fix_saturation) { fix_saturation_ = fix_saturation; }

void HcalTriggerPrimitiveAlgo::overrideParameters(const edm::ParameterSet& ps) {
  override_parameters_ = ps;

  if (override_parameters_.exists("ADCThresholdHF")) {
    override_adc_hf_ = true;
    override_adc_hf_value_ = override_parameters_.getParameter<uint32_t>("ADCThresholdHF");
  }
  if (override_parameters_.exists("TDCMaskHF")) {
    override_tdc_hf_ = true;
    override_tdc_hf_value_ = override_parameters_.getParameter<unsigned long long>("TDCMaskHF");
  }
}

void HcalTriggerPrimitiveAlgo::addSignal(const QIE10DataFrame& frame) {
  DetId detId = DetId(frame.detid());
  if (detId.det() == DetId::Hcal) {
    HcalDetId detId = frame.detid();
    // prevent QIE10 calibration channels from entering TP emulation
    if (detId.subdet() != HcalForward)
      return;

    auto ids = theTrigTowerGeometry->towerIds(frame.id());
    for (const auto& id : ids) {
      if (id.version() == 0) {
        edm::LogError("HcalTPAlgo") << "Encountered QIE10 data frame mapped to TP version 0:" << id;
        continue;
      }
      int nsamples = frame.samples();

      IntegerCaloSamples samples(id, nsamples);
      samples.setPresamples(frame.presamples());
      incoder_->adc2Linear(frame, samples, false);

      // Don't add to final collection yet
      // HF PMT veto sum is calculated in analyzerHF()
      IntegerCaloSamples zero_samples(id, nsamples);
      zero_samples.setPresamples(frame.presamples());
      addSignal(zero_samples);

      auto fid = HcalDetId(frame.id());
      auto& details = theHFUpgradeDetailMap[id][fid.maskDepth()];
      auto& detail = details[fid.depth() - 1];
      detail.samples = samples;
      detail.digi = frame;
      detail.validity.resize(nsamples);
      detail.passTDC.resize(nsamples);
      incoder_->lookupMSB(frame, detail.fgbits);
      for (int idx = 0; idx < nsamples; ++idx) {
        detail.validity[idx] = validChannel(frame, idx);
        detail.passTDC[idx] = passTDC(frame, idx);
      }
    }
  } else if (detId.det() == DetId::Calo && detId.subdetId() == HcalZDCDetId::SubdetectorId) {
    HcalZDCDetId detId = frame.detid();
    // skip RPD Channels
    if (detId.section() != HcalZDCDetId::EM && detId.section() != HcalZDCDetId::HAD) {
      return;
    }
    // skip FSC and dummy channels
    if (detId.section() == HcalZDCDetId::EM && detId.channel() > 5) {
      return;
    }

    auto ids = theTrigTowerGeometry->towerIds_ZDC(frame.id());
    for (const auto& id : ids) {
      int nsamples = frame.samples();

      IntegerCaloSamples samples(id, nsamples);
      IntegerCaloSamples samples_PUsub(id, nsamples);

      samples.setPresamples(frame.presamples());
      samples_PUsub.setPresamples(frame.presamples());

      incoder_->adc2Linear(frame, samples, false);
      incoder_->adc2Linear(frame, samples_PUsub, true);

      for (int i = 1; i < samples.size(); ++i) {
        if (samples_PUsub[i - 1] > samples[i])
          samples[i] = 0;
        else
          samples[i] -= samples_PUsub[i - 1];
      }

      addSignal(samples);
    }
  }
}

void HcalTriggerPrimitiveAlgo::addSignal(const QIE11DataFrame& frame) {
  HcalDetId detId(frame.id());
  // prevent QIE11 calibration channels from entering TP emulation
  if (detId.subdet() != HcalEndcap && detId.subdet() != HcalBarrel)
    return;

  std::vector<HcalTrigTowerDetId> ids = theTrigTowerGeometry->towerIds(detId);
  assert(ids.size() == 1 || ids.size() == 2);
  IntegerCaloSamples samples1(ids[0], int(frame.samples()));

  samples1.setPresamples(frame.presamples());
  incoder_->adc2Linear(frame, samples1);

  std::vector<std::bitset<2>> msb(frame.samples(), 0);
  incoder_->lookupMSB(frame, msb);

  if (ids.size() == 2) {
    // make a second trigprim for the other one, and share the energy
    IntegerCaloSamples samples2(ids[1], samples1.size());
    for (int i = 0; i < samples1.size(); ++i) {
      samples1[i] = uint32_t(samples1[i]);
      samples2[i] = samples1[i];
    }
    samples2.setPresamples(frame.presamples());
    addSignal(samples2);
    addUpgradeFG(ids[1], detId.depth(), msb);
    addUpgradeTDCFG(ids[1], frame);
  }
  addSignal(samples1);
  addUpgradeFG(ids[0], detId.depth(), msb);
  addUpgradeTDCFG(ids[0], frame);
}

void HcalTriggerPrimitiveAlgo::addSignal(const IntegerCaloSamples& samples) {
  HcalTrigTowerDetId id(samples.id());
  SumMap::iterator itr = theSumMap.find(id);

  if (itr == theSumMap.end()) {
    theSumMap.insert(std::make_pair(id, samples));
  } else {
    // wish CaloSamples had a +=
    for (int i = 0; i < samples.size(); ++i) {
      (itr->second)[i] += samples[i];
    }
  }

  // if fix_saturation == true, keep track of tower with saturated input LUT
  if (fix_saturation_) {
    SatMap::iterator itr_sat = theSatMap.find(id);

    assert((itr == theSumMap.end()) == (itr_sat == theSatMap.end()));

    if (itr_sat == theSatMap.end()) {
      std::vector<bool> check_sat;
      for (int i = 0; i < samples.size(); ++i) {
        if (!(samples[i] < QIE11_LINEARIZATION_ET)) {
          check_sat.push_back(true);
        } else
          check_sat.push_back(false);
      }
      theSatMap.insert(std::make_pair(id, check_sat));
    } else {
      for (int i = 0; i < samples.size(); ++i) {
        if (!(samples[i] < QIE11_LINEARIZATION_ET))
          (itr_sat->second)[i] = true;
      }
    }
  }
}

void HcalTriggerPrimitiveAlgo::analyzeQIE11(IntegerCaloSamples& samples,
                                            std::vector<bool> sample_saturation,
                                            HcalTriggerPrimitiveDigi& result,
                                            const HcalFinegrainBit& fg_algo) {
  HcalDetId detId(samples.id());

  // Get the |ieta| for current sample
  int theIeta = detId.ietaAbs();

  unsigned int dgSamples = samples.size();
  unsigned int dgPresamples = samples.presamples();

  unsigned int tpSamples = numberOfSamples_;
  unsigned int tpPresamples = numberOfPresamples_;

  unsigned int filterSamples = weightsQIE11_[theIeta].size();
  unsigned int filterPresamples = numberOfFilterPresamplesHBQIE11_;

  unsigned int shift = dgPresamples - tpPresamples;

  // shrink keeps the FIR filter from going off the end of the 8TS vector
  unsigned int shrink = filterSamples - 1;

  auto& msb = fgUpgradeMap_[samples.id()];
  auto& timingTDC = fgUpgradeTDCMap_[samples.id()];
  IntegerCaloSamples sum(samples.id(), samples.size());

  std::vector<HcalTrigTowerDetId> ids = theTrigTowerGeometry->towerIds(detId);

  // keep track of tower with saturated energy and force the total TP saturated
  bool force_saturation[samples.size()];
  for (int i = 0; i < samples.size(); i++) {
    force_saturation[i] = false;
  }

  //slide algo window
  for (unsigned int ibin = 0; ibin < dgSamples - shrink; ++ibin) {
    int algosumvalue = 0;
    bool check_sat = false;
    //TP energy calculation for PFA2
    if (weightsQIE11_[theIeta][0] == 255) {
      for (unsigned int i = 0; i < filterSamples; i++) {
        //add up value * scale factor
        // In addition, divide by two in the 10 degree phi segmentation region
        // to mimic 5 degree segmentation for the trigger
        unsigned int sample = samples[ibin + i];

        if (fix_saturation_ && (sample_saturation.size() > ibin + i))
          check_sat = (check_sat | sample_saturation[ibin + i] | (sample > QIE11_MAX_LINEARIZATION_ET));

        if (sample > QIE11_MAX_LINEARIZATION_ET)
          sample = QIE11_MAX_LINEARIZATION_ET;

        // Usually use a segmentation factor of 1.0 but for ieta >= 21 use 2
        int segmentationFactor = 1;
        if (ids.size() == 2) {
          segmentationFactor = 2;
        }

        algosumvalue += int(sample / segmentationFactor);
      }
      if (algosumvalue < 0)
        sum[ibin] = 0;  // low-side
                        //high-side
      //else if (algosumvalue>QIE11_LINEARIZATION_ET) sum[ibin]=QIE11_LINEARIZATION_ET;
      else
        sum[ibin] = algosumvalue;  //assign value to sum[]

      if (check_sat)
        force_saturation[ibin] = true;
      //TP energy calculation for PFA1' and PFA1
    } else {
      //add up value * scale factor
      // In addition, divide by two in the 10 degree phi segmentation region
      // to mimic 5 degree segmentation for the trigger
      int sampleTS = samples[ibin + 1];
      int sampleTSminus1 = samples[ibin];

      if (fix_saturation_ && (sample_saturation.size() > ibin + 1))
        check_sat |= sample_saturation[ibin + 1] | (sampleTS >= QIE11_MAX_LINEARIZATION_ET);

      if (sampleTS > QIE11_MAX_LINEARIZATION_ET)
        sampleTS = QIE11_MAX_LINEARIZATION_ET;

      if (sampleTSminus1 > QIE11_MAX_LINEARIZATION_ET || sample_saturation[ibin])
        sampleTSminus1 = QIE11_MAX_LINEARIZATION_ET;

      // Usually use a segmentation factor of 1.0 but for ieta >= 21 use 2
      int segmentationFactor = 1;
      if (ids.size() == 2) {
        segmentationFactor = 2;
      }

      // Based on the |ieta| of the sample, retrieve the correct region weight
      int theWeight = weightsQIE11_[theIeta][0];

      algosumvalue = ((sampleTS << 8) - (sampleTSminus1 * theWeight)) / 256 / segmentationFactor;

      if (algosumvalue < 0)
        sum[ibin] = 0;  // low-side
                        //high-side
      //else if (algosumvalue>QIE11_LINEARIZATION_ET) sum[ibin]=QIE11_LINEARIZATION_ET;
      else
        sum[ibin] = algosumvalue;  //assign value to sum[]

      if (check_sat)
        force_saturation[ibin] = true;
    }
  }

  std::vector<int> finegrain(tpSamples, false);

  IntegerCaloSamples output(samples.id(), tpSamples);
  output.setPresamples(tpPresamples);

  // Based on the |ieta| of the sample, retrieve the correct region "coded" veto threshold
  // where two of the possible values have special meaning
  unsigned int codedVetoThreshold = codedVetoThresholds_[theIeta];

  // Anything in range (1, 2048) inclusive shall activate the veto
  unsigned int actualVetoThreshold = codedVetoThreshold;
  bool applyVetoThreshold = codedVetoThreshold > 0 && codedVetoThreshold <= 2048;

  // Special value to disable vetoing in the PFA1' algo is 0
  if (codedVetoThreshold > 0) {
    if (codedVetoThreshold <= 2048) {
      // Special value to run the veto in PFA1' with no threshold
      if (codedVetoThreshold == 2048)
        actualVetoThreshold = 0;
    } else {
      edm::LogWarning("HcalTPAlgo") << "Specified veto threshold value " << codedVetoThreshold
                                    << " is not in range (1, 2048) ! Vetoing in PFA1' will not be enabled !";
    }
  }

  for (unsigned int ibin = 0; ibin < tpSamples; ++ibin) {
    // ibin - index for output TP
    // idx - index for samples + shift - filterPresamples
    int idx = ibin + shift - filterPresamples;

    // When idx is <= 0 peakfind would compare out-of-bounds of the vector. Avoid this ambiguity
    if (idx <= 0) {
      output[ibin] = 0;
      continue;
    }

    //Only run the peak-finder when the PFA2 FIR filter is running, which corresponds to weights = 1
    if (weightsQIE11_[theIeta][0] == 255) {
      bool isPeak = (sum[idx] > sum[idx - 1] && sum[idx] >= sum[idx + 1] && sum[idx] > theThreshold);
      if (isPeak) {
        output[ibin] = std::min<unsigned int>(sum[idx], QIE11_MAX_LINEARIZATION_ET);

        if (fix_saturation_ && force_saturation[idx] && ids.size() == 2)
          output[ibin] = QIE11_MAX_LINEARIZATION_ET / 2;
        else if (fix_saturation_ && force_saturation[idx])
          output[ibin] = QIE11_MAX_LINEARIZATION_ET;

      } else {
        // Not a peak
        output[ibin] = 0;
      }
    } else {
      // Only if the sum for the future time sample is above the veto
      // threshold and the now sum is not a peak and the now sum is not
      // saturated does the current sum get zeroed
      if (applyVetoThreshold && sum[idx + 1] >= actualVetoThreshold &&
          (sum[idx] < sum[idx + 1] || force_saturation[idx + 1]) && !force_saturation[idx])
        output[ibin] = 0;
      else {
        // Here, either the "now" sum is a peak or the vetoing criteria are not satisfied
        // so assign the appropriate sum to the output
        output[ibin] = std::min<unsigned int>(sum[idx], QIE11_MAX_LINEARIZATION_ET);
        if (fix_saturation_ && force_saturation[idx]) {
          output[ibin] = QIE11_MAX_LINEARIZATION_ET;
          if (ids.size() == 2)
            output[ibin] /= 2;
        }
      }
    }
    // peak-finding is not applied for FG bits
    // compute(msb) returns two bits (MIP). compute(timingTDC,ids) returns 6 bits (1 depth, 1 prompt, 1 delayed 01, 1 delayed 10, 2 reserved)
    finegrain[ibin] = fg_algo.compute(timingTDC[idx + filterPresamples], ids[0]).to_ulong() |
                      fg_algo.compute(msb[idx + filterPresamples]).to_ulong() << 4;
    if (ibin == tpPresamples && (idx + filterPresamples) != dgPresamples)
      edm::LogError("HcalTriggerPritimveAlgo")
          << "TP SOI (tpPresamples = " << tpPresamples
          << ") is not aligned with digi SOI (dgPresamples = " << dgPresamples << ")";
  }
  outcoder_->compress(output, finegrain, result);
}

void HcalTriggerPrimitiveAlgo::analyzeZDC(IntegerCaloSamples& samples, HcalTriggerPrimitiveDigi& result) {
  HcalTrigTowerDetId detId(samples.id());

  unsigned int tpSamples;
  unsigned int tpPresamples;

  tpSamples = samples.size();
  tpPresamples = samples.presamples();
  result.setSize(tpSamples);
  result.setPresamples(tpPresamples);

  IntegerCaloSamples output(samples.id(), tpSamples);
  output.setPresamples(tpPresamples);

  for (int i = 0; i < samples.size(); i++) {
    if (samples[i] > QIE10_ZDC_MAX_LINEARIZATION_ET)
      output[i] = QIE10_ZDC_MAX_LINEARIZATION_ET;
    else
      output[i] = samples[i];
    HcalTriggerPrimitiveSample zdcSample(output[i]);
    result.setSample(i, zdcSample);
  }
}

bool HcalTriggerPrimitiveAlgo::passTDC(const QIE10DataFrame& digi, int ts) const {
  auto parameters = conditions_->getHcalTPParameters();
  auto adc_threshold = parameters->getADCThresholdHF();
  auto tdc_mask = parameters->getTDCMaskHF();

  if (override_adc_hf_)
    adc_threshold = override_adc_hf_value_;
  if (override_tdc_hf_)
    tdc_mask = override_tdc_hf_value_;

  if (digi[ts].adc() < adc_threshold)
    return true;

  return (1ul << digi[ts].le_tdc()) & tdc_mask;
}

bool HcalTriggerPrimitiveAlgo::validChannel(const QIE10DataFrame& digi, int ts) const {
  // channels with invalid data should not contribute to the sum
  if (digi.linkError() || ts >= digi.samples() || !digi[ts].ok())
    return false;

  auto mask = conditions_->getHcalTPChannelParameter(HcalDetId(digi.id()))->getMask();
  if (mask)
    return false;

  return true;
}

void HcalTriggerPrimitiveAlgo::analyzeHFQIE10(const IntegerCaloSamples& samples,
                                              HcalTriggerPrimitiveDigi& result,
                                              const int hf_lumi_shift,
                                              const HcalFeatureBit* embit) {
  // Align digis and TP
  const int shift = samples.presamples() - numberOfPresamplesHF_;
  assert(shift >= 0);
  assert((shift + numberOfSamplesHF_) <= samples.size());
  assert(hf_lumi_shift >= 2);

  // Try to find the HFDetails from the map corresponding to our samples
  const HcalTrigTowerDetId detId(samples.id());
  auto it = theHFUpgradeDetailMap.find(detId);
  // Missing values will give an empty digi
  if (it == theHFUpgradeDetailMap.end()) {
    return;
  }

  std::vector<std::bitset<2>> finegrain(numberOfSamplesHF_, false);

  // Set up out output of IntergerCaloSamples
  IntegerCaloSamples output(samples.id(), numberOfSamplesHF_);
  output.setPresamples(numberOfPresamplesHF_);

  for (const auto& item : it->second) {
    auto& details = item.second;
    for (int ibin = 0; ibin < numberOfSamplesHF_; ++ibin) {
      const int idx = ibin + shift;

      int long_fiber_val = 0;
      int long_fiber_count = 0;
      int short_fiber_val = 0;
      int short_fiber_count = 0;

      bool saturated = false;

      for (auto i : {0, 2}) {
        if (idx < details[i].samples.size() and details[i].validity[idx] and details[i].passTDC[idx]) {
          long_fiber_val += details[i].samples[idx];
          saturated = saturated || (details[i].samples[idx] == QIE10_LINEARIZATION_ET);
          ++long_fiber_count;
        }
      }
      for (auto i : {1, 3}) {
        if (idx < details[i].samples.size() and details[i].validity[idx] and details[i].passTDC[idx]) {
          short_fiber_val += details[i].samples[idx];
          saturated = saturated || (details[i].samples[idx] == QIE10_LINEARIZATION_ET);
          ++short_fiber_count;
        }
      }

      if (saturated) {
        output[ibin] = QIE10_MAX_LINEARIZATION_ET;
      } else {
        // For details of the energy handling, see:
        // https://cms-docdb.cern.ch/cgi-bin/DocDB/ShowDocument?docid=12306
        // If both readouts are valid, average of the two energies is taken
        // division by 2 is compensated by adjusting the total scale shift in the end
        if (long_fiber_count == 2)
          long_fiber_val >>= 1;
        if (short_fiber_count == 2)
          short_fiber_val >>= 1;

        auto sum = long_fiber_val + short_fiber_val;
        // Similar to above, if both channels are valid,
        // average of the two energies is calculated
        // division by 2 here is also compensated by adjusting the total scale shift in the end
        if (long_fiber_count > 0 and short_fiber_count > 0)
          sum >>= 1;

        output[ibin] += sum;
      }

      for (const auto& detail : details) {
        if (idx < int(detail.digi.size()) and detail.validity[idx] and
            HcalDetId(detail.digi.id()).ietaAbs() >= FIRST_FINEGRAIN_TOWER) {
          if (useTDCInMinBiasBits_ && !detail.passTDC[idx])
            continue;
          finegrain[ibin][1] = finegrain[ibin][1] or detail.fgbits[idx][0];
          // what is commonly called the "second" HF min-bias bit is
          // actually the 0-th bit, which can also be used instead for the EM bit
          // (called finegrain[ibin][0] below) in non-HI running
          finegrain[ibin][0] = finegrain[ibin][0] or detail.fgbits[idx][1];
        }
      }
      // the EM bit is only used if the "second" FG bit is disabled
      if (embit != nullptr and FG_HF_thresholds_.at(1) != 255) {
        finegrain[ibin][0] = embit->fineGrainbit(details[1].digi,
                                                 details[3].digi,
                                                 details[0].digi,
                                                 details[2].digi,
                                                 details[1].validity[idx],
                                                 details[3].validity[idx],
                                                 details[0].validity[idx],
                                                 details[2].validity[idx],
                                                 idx);
      }
    }
  }

  for (int bin = 0; bin < numberOfSamplesHF_; ++bin) {
    output[bin] = std::min({(unsigned int)QIE10_MAX_LINEARIZATION_ET, output[bin] >> (hf_lumi_shift - 2)});
  }
  std::vector<int> finegrain_converted;
  finegrain_converted.reserve(finegrain.size());
  for (const auto& fg : finegrain)
    finegrain_converted.push_back(fg.to_ulong());
  outcoder_->compress(output, finegrain_converted, result);
}

void HcalTriggerPrimitiveAlgo::addUpgradeFG(const HcalTrigTowerDetId& id,
                                            int depth,
                                            const std::vector<std::bitset<2>>& bits) {
  auto it = fgUpgradeMap_.find(id);
  if (it == fgUpgradeMap_.end()) {
    FGUpgradeContainer element;
    element.resize(bits.size());
    it = fgUpgradeMap_.insert(std::make_pair(id, element)).first;
  }
  for (unsigned int i = 0; i < bits.size(); ++i) {
    it->second[i][0][depth - 1] = bits[i][0];
    it->second[i][1][depth - 1] = bits[i][1];
  }
}

void HcalTriggerPrimitiveAlgo::addUpgradeTDCFG(const HcalTrigTowerDetId& id, const QIE11DataFrame& frame) {
  HcalDetId detId(frame.id());
  if (detId.subdet() != HcalEndcap && detId.subdet() != HcalBarrel)
    return;

  std::vector<HcalTrigTowerDetId> ids = theTrigTowerGeometry->towerIds(detId);
  assert(ids.size() == 1 || ids.size() == 2);
  IntegerCaloSamples samples1(ids[0], int(frame.samples()));
  samples1.setPresamples(frame.presamples());
  incoder_->adc2Linear(frame, samples1);                                  // use linearization LUT
  std::vector<unsigned short> bits12_15 = incoder_->group0FGbits(frame);  // get 4 energy bits (12-15) from group 0 LUT

  bool is_compressed = false;
  if (detId.subdet() == HcalBarrel) {
    is_compressed = (frame.flavor() == 3);
    // 0 if frame.flavor is 0 (uncompressed), 1 if frame.flavor is 3 (compressed)
  }

  // Access tdc1, tdc2 thresholds from first two bytes of auxi1 field of HcalTPChannelParameters, if dedicated switch is on
  int tdc1 = 0;
  int tdc2 = 0;
  bool hasTDCThr = false;

  if (useTDCfromDB_) {
    auto tpchparams = conditions_->getHcalTPChannelParameter(detId, false);
    const uint32_t auxi1 = tpchparams->getauxi1();
    const int tp_tdc1 = auxi1 & 0xFF;
    const int tp_tdc2 = (auxi1 >> 8) & 0xFF;

    // accept only if sane and nonzero
    const bool empty = (tp_tdc1 == 0 && tp_tdc2 == 0);
    const bool bad = (tp_tdc2 < tp_tdc1) || (tp_tdc2 > tdcmax_);

    if (empty || bad) {
      throw cms::Exception("HBTDCThresholds") << "Missing/invalid HB TDC thresholds from HcalTPChannelParameters for "
                                              << detId << " (auxi1: tdc1/tdc2 not set or out of range).";
    }

    tdc1 = tp_tdc1;
    tdc2 = tp_tdc2;
    hasTDCThr = true;
  }

  // Pack bits12_15 with tdc1, tdc2, and hasTDCThr
  auto packBits = [&](unsigned short b1215) -> int {
    int packed = (int(b1215) & 0xF);        // bits 0..3   : original bits12-15
    packed |= ((tdc1 & 0xFF) << 4);         // bits 4..11  : tdc1
    packed |= ((tdc2 & 0xFF) << 12);        // bits 12..19 : tdc2
    packed |= ((hasTDCThr ? 1 : 0) << 20);  // bit  20     : hasTDCThr flag
    return packed;
  };

  auto it = fgUpgradeTDCMap_.find(id);
  if (it == fgUpgradeTDCMap_.end()) {
    FGUpgradeTDCContainer element;
    element.resize(frame.samples());
    it = fgUpgradeTDCMap_.insert(std::make_pair(id, element)).first;
  }
  for (int i = 0; i < frame.samples(); i++) {
    it->second[i][detId.depth() - 1] = std::make_pair(std::make_pair(packBits(bits12_15[i]), is_compressed),
                                                      std::make_pair(frame[i].tdc(), samples1[i]));
  }
}

void HcalTriggerPrimitiveAlgo::setWeightsQIE11(const edm::ParameterSet& weightsQIE11) {
  // Names are just abs(ieta) for HB
  std::vector<std::string> ietaStrs = weightsQIE11.getParameterNames();
  for (auto& ietaStr : ietaStrs) {
    // Strip off "ieta" part of key and just use integer value in map
    auto const& v = weightsQIE11.getParameter<std::vector<int>>(ietaStr);
    weightsQIE11_[std::stoi(ietaStr.substr(4))] = {{v[0], v[1]}};
  }
}

void HcalTriggerPrimitiveAlgo::setWeightQIE11(int aieta, int weight) {
  // Simple map of |ieta| in HB to weight
  // Only one weight for SOI-1 TS
  weightsQIE11_[aieta] = {{weight, 255}};
}

void HcalTriggerPrimitiveAlgo::setCodedVetoThresholds(const edm::ParameterSet& codedVetoThresholds) {
  // Names are just abs(ieta) for HB
  std::vector<std::string> ietaStrs = codedVetoThresholds.getParameterNames();
  for (auto& ietaStr : ietaStrs) {
    // Strip off "ieta" part of key and just use integer value in map
    auto const& v = codedVetoThresholds.getParameter<int>(ietaStr);
    codedVetoThresholds_[std::stoi(ietaStr.substr(4))] = {v};
  }
}

void HcalTriggerPrimitiveAlgo::setCodedVetoThreshold(int aieta, int codedVetoThreshold) {
  // Simple map of |ieta| in HB to veto threshold
  codedVetoThresholds_[aieta] = {codedVetoThreshold};
}

void HcalTriggerPrimitiveAlgo::setNCTScaleShift(int shift) { NCTScaleShift = shift; }

void HcalTriggerPrimitiveAlgo::setRCTScaleShift(int shift) { RCTScaleShift = shift; }

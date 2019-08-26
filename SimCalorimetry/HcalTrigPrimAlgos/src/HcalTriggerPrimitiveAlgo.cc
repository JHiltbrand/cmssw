#include "SimCalorimetry/HcalTrigPrimAlgos/interface/HcalTriggerPrimitiveAlgo.h"

#include "CalibFormats/CaloObjects/interface/IntegerCaloSamples.h"
#include "CondFormats/HcalObjects/interface/HcalTPParameters.h"
#include "CondFormats/HcalObjects/interface/HcalTPChannelParameters.h"

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"

#include "EventFilter/HcalRawToDigi/interface/HcalDCCHeader.h"
#include "EventFilter/HcalRawToDigi/interface/HcalHTRData.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"

#include <iostream>
#include <numeric>

using namespace std;

std::pair<double, double>
timing(const QIE11DataFrame& frame) {
  // edm::DataFrame::size_type 
  int n = (int)frame.size();
  double ft = -999.;
  double rt = -999.;
  int sig_bx = (int)frame.presamples();     

  int dir = -1; int step = 1;  
  int i = sig_bx;

  int nbins = 500;

  while ((i > 2) && (i < (int)frame.size() - 2) && (i < n) && ((rt < -998.) || (ft < 998.))) {

    unsigned rise = frame[i].tdc() % 100;
    unsigned fall = frame[i].tdc() / 100;

    //    printf("runnin i=%d with frame_size=%d, rise=%u and fall=%u\n",i,n,rise,fall);  

    if (rt < -998. && rise != 62 && rise != 63) {
      rt = rise * 25. / nbins + (i - sig_bx) * 25.;
    }
    if (((ft < -998.) || (ft < rt))  && 
	(fall != 62) && (fall != 63)) {
      ft = fall * 25. / nbins + (i - sig_bx) * 25.;
    }

    i += dir * step;
    ++step;
    dir *= -1;
  }
  /*
  if (rt > -998 or ft > -998) 
  std::cout << "rise " << rt << " fall " << ft << std::endl; 
  */
  return {rt, ft};
}

void
update(HcalUpgradeTriggerPrimitiveDigi& digi, const Sample& sample, int soi, std::vector<int> lineardepthdata, std::vector<int> linearized)
{

  //   std::vector<int> linearized;

  std::vector<double> rise_avg;
  std::vector<double> rise_rms;
  std::vector<double> fall_avg;
  std::vector<double> fall_rms;
  std::vector<int> oot;

  for (int i = 0; i < 8; ++i) {
    //     linearized.push_back(sample[i][soi]);

    auto rise = sample.rise(i);
    if (rise.size() > 0) {
      double avg = std::accumulate(rise.begin(), rise.end(), 0.) / rise.size();
      double sqrs = std::accumulate(rise.begin(), rise.end(), 0., [](double x, double y) { return x + y * y; });

      rise_avg.push_back(avg);
      rise_rms.push_back(sqrt(sqrs / rise.size() - avg * avg));
    } else {
      rise_avg.push_back(-1e6);
      rise_rms.push_back(-1e6);
    }

    auto fall = sample.fall(i);
    if (fall.size() > 0) {
      double avg = std::accumulate(fall.begin(), fall.end(), 0.) / fall.size();
      double sqrs = std::accumulate(fall.begin(), fall.end(), 0., [](double x, double y) { return x + y * y; });

      fall_avg.push_back(avg);
      fall_rms.push_back(sqrt(sqrs / fall.size() - avg * avg));
    } else {
      fall_avg.push_back(-1e6);
      fall_rms.push_back(-1e6);
    }
    oot.push_back(sample(i)[soi]);

  }

  digi.setDepthData(lineardepthdata);

  digi.setSampleData(linearized);

  digi.setOOTData(oot);
  digi.setTimingData(rise_avg, rise_rms, fall_avg, fall_rms);
}

void
update_legacy(HcalUpgradeTriggerPrimitiveDigi& digi, const Sample& sample, int soi)
{

  std::vector<int> linearized;

  for (int i = 0; i < 3; ++i) {
    linearized.push_back(sample[i][soi]);
  }

  digi.setDepthData(linearized);

}

HcalTriggerPrimitiveAlgo::HcalTriggerPrimitiveAlgo( bool pf, int latency,
                                                    uint32_t FG_threshold, const std::vector<uint32_t>& FG_HF_thresholds, uint32_t ZS_threshold,
                                                    int numberOfSamples, int numberOfPresamples,
                                                    int numberOfSamplesHF, int numberOfPresamplesHF, bool useTDCInMinBiasBits,
                                                    uint32_t minSignalThreshold, uint32_t PMT_NoiseThreshold,
						    bool upgrade
                                                    )
                                                   : incoder_(nullptr), outcoder_(nullptr),
                                                   theThreshold(0), peakfind_(pf), latency_(latency),
                                                   FG_threshold_(FG_threshold), FG_HF_thresholds_(FG_HF_thresholds), ZS_threshold_(ZS_threshold),
                                                   numberOfSamples_(numberOfSamples),
                                                   numberOfPresamples_(numberOfPresamples),
                                                   numberOfSamplesHF_(numberOfSamplesHF),
                                                   numberOfPresamplesHF_(numberOfPresamplesHF),
                                                   useTDCInMinBiasBits_(useTDCInMinBiasBits),
                                                   minSignalThreshold_(minSignalThreshold),
                                                   PMT_NoiseThreshold_(PMT_NoiseThreshold),
                                                   NCTScaleShift(0), RCTScaleShift(0),
                                                   upgrade_(upgrade),  
                                                   peak_finder_algorithm_(2),
                                                   peak_finder_algorithm_name_("PFA2p"),
                                                   override_parameters_()
{
   //No peak finding setting (for Fastsim)
   if (!peakfind_){
      numberOfSamples_ = 1; 
      numberOfPresamples_ = 0;
      numberOfSamplesHF_ = 1; 
      numberOfPresamplesHF_ = 0;
   }
   // Switch to integer for comparisons - remove compiler warning
   ZS_threshold_I_ = ZS_threshold_;
}


HcalTriggerPrimitiveAlgo::~HcalTriggerPrimitiveAlgo() {

}

void
HcalTriggerPrimitiveAlgo::setUpgradeFlags(bool hb, bool he, bool hf)
{
   upgrade_hb_ = hb;
   upgrade_he_ = he;
   upgrade_hf_ = hf;
}


void
HcalTriggerPrimitiveAlgo::overrideParameters(const edm::ParameterSet& ps)
{
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


void HcalTriggerPrimitiveAlgo::addSignal(const HBHEDataFrame & frame) {
   // TODO: Need to add support for seperate 28, 29 in HE
   //Hack for 300_pre10, should be removed.
   if (frame.id().depth()==5) return;

   std::vector<HcalTrigTowerDetId> ids = theTrigTowerGeometry->towerIds(frame.id());
   assert(ids.size() == 1 || ids.size() == 2);
   IntegerCaloSamples samples1(ids[0], int(frame.size()));

   samples1.setPresamples(frame.presamples());
   incoder_->adc2Linear(frame, samples1);

   std::vector<bool> msb;
   incoder_->lookupMSB(frame, msb);

   if(ids.size() == 2) {
      // make a second trigprim for the other one, and split the energy
      IntegerCaloSamples samples2(ids[1], samples1.size());
      for(int i = 0; i < samples1.size(); ++i) {
         samples1[i] = uint32_t(samples1[i]*0.5);
         samples2[i] = samples1[i];
      }
      samples2.setPresamples(frame.presamples());
      addSignal(samples2);
      addFG(ids[1], msb);
   }
   addSignal(samples1);
   addFG(ids[0], msb);
}


void HcalTriggerPrimitiveAlgo::addSignal(const HFDataFrame & frame) {
   if(frame.id().depth() == 1 || frame.id().depth() == 2) {
      std::vector<HcalTrigTowerDetId> ids = theTrigTowerGeometry->towerIds(frame.id());
      std::vector<HcalTrigTowerDetId>::const_iterator it;
      for (it = ids.begin(); it != ids.end(); ++it) {
         HcalTrigTowerDetId trig_tower_id = *it;
         IntegerCaloSamples samples(trig_tower_id, frame.size());
         samples.setPresamples(frame.presamples());
         incoder_->adc2Linear(frame, samples);

         // Don't add to final collection yet
         // HF PMT veto sum is calculated in analyzerHF()
         IntegerCaloSamples zero_samples(trig_tower_id, frame.size());
         zero_samples.setPresamples(frame.presamples());
         addSignal(zero_samples);

         // Pre-LS1 Configuration
         if (trig_tower_id.version() == 0) {
            // Mask off depths: fgid is the same for both depths
            uint32_t fgid = (frame.id().maskDepth());

            if ( theTowerMapFGSum.find(trig_tower_id) == theTowerMapFGSum.end() ) {
               SumFGContainer sumFG;
               theTowerMapFGSum.insert(std::pair<HcalTrigTowerDetId, SumFGContainer>(trig_tower_id, sumFG));
            }

            SumFGContainer& sumFG = theTowerMapFGSum[trig_tower_id];
            SumFGContainer::iterator sumFGItr;
            for ( sumFGItr = sumFG.begin(); sumFGItr != sumFG.end(); ++sumFGItr) {
               if (sumFGItr->id() == fgid) { break; }
            }
            // If find
            if (sumFGItr != sumFG.end()) {
               for (int i=0; i<samples.size(); ++i) {
                  (*sumFGItr)[i] += samples[i];
               }
            }
            else {
               //Copy samples (change to fgid)
               IntegerCaloSamples sumFGSamples(DetId(fgid), samples.size());
               sumFGSamples.setPresamples(samples.presamples());
               for (int i=0; i<samples.size(); ++i) {
                  sumFGSamples[i] = samples[i];
               }
               sumFG.push_back(sumFGSamples);
            }

            // set veto to true if Long or Short less than threshold
            if (HF_Veto.find(fgid) == HF_Veto.end()) {
               vector<bool> vetoBits(samples.size(), false);
               HF_Veto[fgid] = vetoBits;
            }
            for (int i=0; i<samples.size(); ++i) {
               if (samples[i] < minSignalThreshold_) {
                  HF_Veto[fgid][i] = true;
               }
            }
         }
         // HF 1x1
         else if (trig_tower_id.version() == 1) {
            uint32_t fgid = (frame.id().maskDepth());
            HFDetails& details = theHFDetailMap[trig_tower_id][fgid];
            // Check the frame type to determine long vs short
            if (frame.id().depth() == 1) { // Long
               details.long_fiber = samples;
               details.LongDigi = frame;
            } else if (frame.id().depth() == 2) { // Short
               details.short_fiber = samples;
               details.ShortDigi = frame;
            } else {
                // Neither long nor short... So we have no idea what to do
                edm::LogWarning("HcalTPAlgo") << "Unable to figure out what to do with data frame for " << frame.id();
                return;
            }
         }
         // Uh oh, we are in a bad/unknown state! Things will start crashing.
         else {
             return;
         }
      }
   }
}

void
HcalTriggerPrimitiveAlgo::addSignal(const QIE10DataFrame& frame)
{
   HcalDetId detId = frame.detid();
   // prevent QIE10 calibration channels from entering TP emulation
   if(detId.subdet() != HcalForward) return;

   auto ids = theTrigTowerGeometry->towerIds(frame.id());
   for (const auto& id: ids) {
      if (id.version() == 0) {
         edm::LogError("HcalTPAlgo") << "Encountered QIE10 data frame mapped to TP version 0:" << id;
         continue;
      }

      int nsamples=frame.samples();

      IntegerCaloSamples samples(id, nsamples);
      samples.setPresamples(frame.presamples());
      incoder_->adc2Linear(frame, samples);

      // Don't add to final collection yet
      // HF PMT veto sum is calculated in analyzerHF()
      IntegerCaloSamples zero_samples(id, nsamples);
      zero_samples.setPresamples(frame.presamples());
      addSignal(zero_samples);

      auto fid = HcalDetId(frame.id());
      auto& details = theHFUpgradeDetailMap[id][fid.maskDepth()];
      auto& detail = details[fid.depth()-1];
      detail.samples = samples;
      detail.digi = frame;
      detail.validity.resize(nsamples);
      detail.passTDC.resize(nsamples);
      incoder_->lookupMSB(frame, detail.fgbits);
      for (int idx = 0; idx < nsamples; ++idx){
         detail.validity[idx] = validChannel(frame, idx);
         detail.passTDC[idx] = passTDC(frame, idx);
      }
   }
}

void
HcalTriggerPrimitiveAlgo::addSignal(const QIE11DataFrame& frame)
{
   HcalDetId detId(frame.id());

   int depth = HcalDetId(frame.id()).depth();
   // prevent QIE11 calibration channels from entering TP emulation
   if(detId.subdet() != HcalEndcap && detId.subdet() != HcalBarrel) return;

   std::vector<HcalTrigTowerDetId> ids = theTrigTowerGeometry->towerIds(detId);
   assert(ids.size() == 1 || ids.size() == 2);
   IntegerCaloSamples samples1(ids[0], int(frame.samples()));

   samples1.setPresamples(frame.presamples());
   incoder_->adc2Linear(frame, samples1);

   std::vector<std::bitset<2>> msb(frame.samples(), 0);
   incoder_->lookupMSB(frame, msb);

   if(ids.size() == 2) {
      // make a second trigprim for the other one, and share the energy
      IntegerCaloSamples samples2(ids[1], samples1.size());
      for(int i = 0; i < samples1.size(); ++i) {
         samples1[i] = uint32_t(samples1[i]);
         samples2[i] = samples1[i];
      }
      samples2.setPresamples(frame.presamples());
      addSignal(samples2, depth, timing(frame));
      //      addSignal(samples2);
      addUpgradeFG(ids[1], detId.depth(), msb);
   }
   //   addSignal(samples1);
   addSignal(samples1, depth, timing(frame));
   addUpgradeFG(ids[0], detId.depth(), msb);
}

//void HcalTriggerPrimitiveAlgo::addSignal(const IntegerCaloSamples & samples) {
void HcalTriggerPrimitiveAlgo::addSignal(const IntegerCaloSamples & samples, int depth, const std::pair<double, double>& tdc) {
   HcalTrigTowerDetId id(samples.id());
   SumMap::iterator itr = theSumMap.find(id);
   if(itr == theSumMap.end()) {
      theSumMap.insert(std::make_pair(id, samples));
      theDepthMap.insert(std::make_pair(id, Sample()));
   }
   else {
      // wish CaloSamples had a +=
      for(int i = 0; i < samples.size(); ++i) {
         (itr->second)[i] += samples[i];
      }
   }
   theDepthMap[id].add(depth, samples, tdc);
}


void HcalTriggerPrimitiveAlgo::analyze(IntegerCaloSamples & samples, HcalTriggerPrimitiveDigi & result) {
   int shrink = numberOfSamples_ - 1;
   std::vector<bool>& msb = fgMap_[samples.id()];
   IntegerCaloSamples sum(samples.id(), samples.size());

   //slide algo window
   for(int ibin = 0; ibin < int(samples.size())- shrink; ++ibin) {
      int algosumvalue = 0;
      for(int i = 0; i < numberOfSamples_; i++) {
         //add up value * scale factor
         algosumvalue += int(samples[ibin+i]);
      }
      if (algosumvalue<0) sum[ibin]=0;            // low-side
                                                  //high-side
      //else if (algosumvalue>QIE8_LINEARIZATION_ET) sum[ibin]=QIE8_LINEARIZATION_ET;
      else sum[ibin] = algosumvalue;              //assign value to sum[]
   }

   // Align digis and TP
   int dgPresamples=samples.presamples(); 
   int tpPresamples=numberOfPresamples_;
   int shift = dgPresamples - tpPresamples;
   int dgSamples=samples.size();
   int tpSamples=numberOfSamples_;
   if(peakfind_){
       if((shift<shrink) || (shift + tpSamples + shrink > dgSamples - (peak_finder_algorithm_ - 1) )   ){
	    edm::LogInfo("HcalTriggerPrimitiveAlgo::analyze") << 
		"TP presample or size from the configuration file is out of the accessible range. Using digi values from data instead...";
	    shift=shrink;
	    tpPresamples=dgPresamples-shrink;
	    tpSamples=dgSamples-(peak_finder_algorithm_-1)-shrink-shift;
       }
   }

   std::vector<int> finegrain(tpSamples,false);

   IntegerCaloSamples output(samples.id(), tpSamples);
   output.setPresamples(tpPresamples);

   for (int ibin = 0; ibin < tpSamples; ++ibin) {
      // ibin - index for output TP
      // idx - index for samples + shift
      int idx = ibin + shift;

      //Peak finding
      if (peakfind_) {
         bool isPeak = false;
         switch (peak_finder_algorithm_) {
            case 1 :
               isPeak = (samples[idx] > samples[idx-1] && samples[idx] >= samples[idx+1] && samples[idx] > theThreshold);
               break;
            case 2:
               isPeak = (sum[idx] > sum[idx-1] && sum[idx] >= sum[idx+1] && sum[idx] > theThreshold);
               break;
            default:
               break;
         }

         if (isPeak){
            output[ibin] = std::min<unsigned int>(sum[idx],QIE8_LINEARIZATION_ET);
            finegrain[ibin] = msb[idx];
         }
         // Not a peak
         else output[ibin] = 0;
      }
      else { // No peak finding, just output running sum
         output[ibin] = std::min<unsigned int>(sum[idx],QIE8_LINEARIZATION_ET);
         finegrain[ibin] = msb[idx];
      }

      // Only Pegged for 1-TS algo.
      if (peak_finder_algorithm_ == 1) {
         if (samples[idx] >= QIE8_LINEARIZATION_ET)
            output[ibin] = QIE8_LINEARIZATION_ET;
      }
   }
   outcoder_->compress(output, finegrain, result);
}

void HcalTriggerPrimitiveAlgo::analyze(IntegerCaloSamples & samples, HcalUpgradeTriggerPrimitiveDigi & result) {
   int shrink = numberOfSamples_ - 1;
   std::vector<bool>& msb = fgMap_[samples.id()];
   IntegerCaloSamples sum(samples.id(), samples.size());

   //slide algo window
   for(int ibin = 0; ibin < int(samples.size())- shrink; ++ibin) {
      int algosumvalue = 0;
      for(int i = 0; i < numberOfSamples_; i++) {
         //add up value * scale factor
         algosumvalue += int(samples[ibin+i]);
      }
      if (algosumvalue<0) sum[ibin]=0;            // low-side
                                                  //high-side
      //else if (algosumvalue>QIE8_LINEARIZATION_ET) sum[ibin]=QIE8_LINEARIZATION_ET;
      else sum[ibin] = algosumvalue;              //assign value to sum[]
   }

   // Align digis and TP
   int dgPresamples=samples.presamples(); 
   int tpPresamples=numberOfPresamples_;
   int shift = dgPresamples - tpPresamples;
   int dgSamples=samples.size();
   int tpSamples=numberOfSamples_;
   if(peakfind_){
       if((shift<shrink) || (shift + tpSamples + shrink > dgSamples - (peak_finder_algorithm_ - 1) )   ){
	    edm::LogInfo("HcalTriggerPrimitiveAlgo::analyze") << 
		"TP presample or size from the configuration file is out of the accessible range. Using digi values from data instead...";
	    shift=shrink;
	    tpPresamples=dgPresamples-shrink;
	    tpSamples=dgSamples-(peak_finder_algorithm_-1)-shrink-shift;
       }
   }

   std::vector<int> finegrain(tpSamples,false);

   IntegerCaloSamples output(samples.id(), tpSamples);
   output.setPresamples(tpPresamples);

    std::vector<int> depth_sums(3, 0);
   
   for (int ibin = 0; ibin < tpSamples; ++ibin) {
      // ibin - index for output TP
      // idx - index for samples + shift
      int idx = ibin + shift;

      //Peak finding
      if (peakfind_) {
         bool isPeak = false;
         switch (peak_finder_algorithm_) {
            case 1 :
               isPeak = (samples[idx] > samples[idx-1] && samples[idx] >= samples[idx+1] && samples[idx] > theThreshold);
               break;
            case 2:
               isPeak = (sum[idx] > sum[idx-1] && sum[idx] >= sum[idx+1] && sum[idx] > theThreshold);
               break;
            default:
               break;
         }

         if (isPeak){
            output[ibin] = std::min<unsigned int>(sum[idx],QIE8_LINEARIZATION_ET);
            finegrain[ibin] = msb[idx];

	    // Only provide depth information for the SOI.  This is the
            // energy value used downstream, even if the peak is found for
            // another sample.
            if (ibin == numberOfPresamples_) {
               for (int d = 0; d < 3; ++d) {
                  auto algosumvalue = 0;
                  for (int i = 0; i < tpSamples; ++i)
                   algosumvalue   += int(theDepthMap[samples.id()][d][idx + i]);
                  depth_sums[d] += std::min<unsigned int>(algosumvalue, QIE8_LINEARIZATION_ET);
               }
            }
	    
         }
         // Not a peak
         else output[ibin] = 0;
      }
      else { // No peak finding, just output running sum
         output[ibin] = std::min<unsigned int>(sum[idx],QIE8_LINEARIZATION_ET);
         finegrain[ibin] = msb[idx];

	  // See comment above.
         if (ibin == numberOfPresamples_) {
            for (int d = 0; d < 3; ++d) {
               auto algosumvalue = 0;
               for (int i = 0; i < tpSamples; ++i)
                  algosumvalue += int(theDepthMap[samples.id()][d][idx + i]);
               depth_sums[d] += std::min<unsigned int>(algosumvalue, QIE8_LINEARIZATION_ET);
            }
         }
      }

      // Only Pegged for 1-TS algo.
      if (peak_finder_algorithm_ == 1) {
         if (samples[idx] >= QIE8_LINEARIZATION_ET)
            output[ibin] = QIE8_LINEARIZATION_ET;
      }
   }
   outcoder_->compress(output, finegrain, result);

   update_legacy(result, theDepthMap[samples.id()], numberOfPresamples_);
}


void
HcalTriggerPrimitiveAlgo::analyzeQIE11(IntegerCaloSamples& samples, HcalTriggerPrimitiveDigi& result, const HcalFinegrainBit& fg_algo)
{
   int shrink = numberOfSamples_ - 1;
   auto& msb = fgUpgradeMap_[samples.id()];
   IntegerCaloSamples sum(samples.id(), samples.size());

   HcalDetId detId(samples.id());
   std::vector<HcalTrigTowerDetId> ids = theTrigTowerGeometry->towerIds(detId);
   //slide algo window
   for(int ibin = 0; ibin < int(samples.size())- shrink; ++ibin) {
      int algosumvalue = 0;
      for(int i = 0; i < numberOfSamples_; i++) {
	//add up value * scale factor
	// In addition, divide by two in the 10 degree phi segmentation region
	// to mimic 5 degree segmentation for the trigger
	unsigned int sample = samples[ibin+i];
	if(sample>QIE11_MAX_LINEARIZATION_ET) sample = QIE11_MAX_LINEARIZATION_ET;
	if(ids.size()==2) algosumvalue += int(sample * 0.5);
	else algosumvalue += int(sample);
      }
      if (algosumvalue<0) sum[ibin]=0;            // low-side
                                                  //high-side
      //else if (algosumvalue>QIE11_LINEARIZATION_ET) sum[ibin]=QIE11_LINEARIZATION_ET;
      else sum[ibin] = algosumvalue;              //assign value to sum[]
   }

   // Align digis and TP
   int dgPresamples=samples.presamples(); 
   int tpPresamples=numberOfPresamples_;
   int shift = dgPresamples - tpPresamples;
   int dgSamples=samples.size();
   int tpSamples=numberOfSamples_;

   if((shift<shrink) || (shift + tpSamples + shrink > dgSamples - (peak_finder_algorithm_ - 1) )   ){
      edm::LogInfo("HcalTriggerPrimitiveAlgo::analyze") << 
         "TP presample or size from the configuration file is out of the accessible range. Using digi values from data instead...";
      shift=shrink;
      tpPresamples=dgPresamples-shrink;
      tpSamples=dgSamples-(peak_finder_algorithm_-1)-shrink-shift;
   }

   std::vector<int> finegrain(tpSamples,false);

   IntegerCaloSamples output(samples.id(), tpSamples);
   output.setPresamples(tpPresamples);

   for (int ibin = 0; ibin < tpSamples; ++ibin) {
      // ibin - index for output TP
      // idx - index for samples + shift
      int idx = ibin + shift;
      bool isPeak = (sum[idx] > sum[idx-1] && sum[idx] >= sum[idx+1] && sum[idx] > theThreshold);

      if (isPeak){
         output[ibin] = std::min<unsigned int>(sum[idx],QIE11_MAX_LINEARIZATION_ET);
      } else {
         // Not a peak
         output[ibin] = 0;
      }
      // peak-finding is not applied for FG bits
      finegrain[ibin] = fg_algo.compute(msb[idx]).to_ulong();
   }
   outcoder_->compress(output, finegrain, result);
}

void
HcalTriggerPrimitiveAlgo::analyzeQIE11(IntegerCaloSamples& samples, HcalUpgradeTriggerPrimitiveDigi& result, const HcalFinegrainBit& fg_algo)
{
   auto& msb = fgUpgradeMap_[samples.id()];
   IntegerCaloSamples sum(samples.id(), samples.size());

   HcalDetId detId(samples.id());
   std::vector<HcalTrigTowerDetId> ids = theTrigTowerGeometry->towerIds(detId);

   int theIeta = detId.ietaAbs();
   int theIphi = detId.iphi();
   int shrink = numberOfSamples_ - 1;

   // By default assume we do not weight and use presamples.
   unsigned int weightedPresamples = 0;
   if (!weights_.empty()) { weightedPresamples = (weights_.begin()->second).size(); }

   // store individual sample data
   // std::map<int, std::vector<int>> thesampledata; //first=ibin, second=sample
   // Save the 8TS digi as the TP for weight studies
   std::vector<int> thesampledata(8, 0);
   for(int ibin = 0; ibin < int(samples.size()); ++ibin) {
       thesampledata[ibin] = samples[ibin];
   }

   double segmentationFactor = 1.0;
   if (ids.size() == 2) segmentationFactor = 0.5;

   //slide algo window
   for(int ibin = 0; ibin < int(samples.size()-shrink); ++ibin) {
       int algosumvalue = 0;
       for(int i = 0; i < numberOfSamples_; i++) {
	       //add up value * scale factor
	       // In addition, divide by two in the 10 degree phi segmentation region
	       // to mimic 5 degree segmentation for the trigger
	       unsigned int sample = samples[ibin+i]; //samples[ibin+i];
	       //	printf("ibin= %d and i=%d, so ibin+i=%d and sample[ibin+i]=%u\n",ibin,i,ibin+i,sample);    

	       if(sample>QIE11_MAX_LINEARIZATION_ET) sample = QIE11_MAX_LINEARIZATION_ET;

	       algosumvalue += int(sample * segmentationFactor);

       }
       // Use the presamples to get an additional factor to add to algosum
       // Since we are in nested loop, add weighted presample to sum once when i == 0 only
       if (weightedPresamples == 2) {

           int presamplevalue = 0;

           if (ibin >= 1) {presamplevalue += int(samples[ibin-1] * segmentationFactor * weights_[theIeta][1]);}
           if (ibin >= 2) {presamplevalue += int(samples[ibin-2] * segmentationFactor * weights_[theIeta][0]);}
           
           algosumvalue += presamplevalue;
       }

       else if (weightedPresamples == 1) {

           int presamplevalue = 0;

           if (ibin >= 1) {presamplevalue += int(samples[ibin-1] * segmentationFactor * weights_[theIeta][0]);}

           algosumvalue += presamplevalue; 
   
       }
       if (algosumvalue<0) {
	       sum[ibin]=0;            // low-side
       //} else if (algosumvalue>QIE11_LINEARIZATION_ET) { sum[ibin]=QIE11_LINEARIZATION_ET;
       } else {
	       sum[ibin] = algosumvalue;              //assign value to sum[]
       }
   }

   // Align digis and TP
   int dgPresamples=samples.presamples(); //3
   int tpPresamples=numberOfPresamples_;//2
   int shift = dgPresamples - tpPresamples; // shift=1 and shrink=4weights-1 = 3
   int dgSamples=samples.size(); //8
   int tpSamples=numberOfSamples_; //4

   // vector data to store per sample info for reconstructing pulse shape
   //   std::vector<int> sampledata = {0,0,0,0};
   // vector data to store the depth info for the SOI sample only
   std::vector<int> depth_sums(8, 0);
   
   // To include or not to include, that is the question
   //if((shift<shrink) || (shift + tpSamples + shrink > dgSamples - (peak_finder_algorithm_ - 1) )   ){
   //  edm::LogInfo("HcalTriggerPrimitiveAlgo::analyze") << 
   //    "TP presample or size from the configuration file is out of the accessible range. Using digi values from data instead...";
   //  //  shift=shrink; // shift=3
   //  tpPresamples=dgPresamples-shrink;// =0
   //  tpSamples=dgSamples-(peak_finder_algorithm_-1)-shrink-shift; //8-2-3-1=2
   //}
   std::vector<int> finegrain(tpSamples,false);

   IntegerCaloSamples output(samples.id(), tpSamples);
   output.setPresamples(tpPresamples);

   for (int ibin = 0; ibin < tpSamples; ++ibin) {
      // ibin - index for output TP
      // idx - index for samples + shift
     int idx = ibin + shift;

     bool isPeak = (sum[idx] > sum[idx-1] && sum[idx] >= sum[idx+1] && sum[idx] > theThreshold);

     if (isPeak){
         output[ibin] = std::min<unsigned int>(sum[idx],QIE11_MAX_LINEARIZATION_ET);

         //       printf("++++++++InPeak, sum is=%u\n",sum[idx]);
         /*
         std::map<int, std::vector<int>>::iterator it;
         for ( it = thesampledata.begin(); it != thesampledata.end(); it++ ){
	 if  ((it->first)==idx) {
	     sampledata=it->second;
	 }
         }
         */
         // Only provide depth information for the SOI.  This is the
         // energy value used downstream, even if the peak is found for
         // another sample.
         if (ibin == numberOfPresamples_) {
	         for (int d = 0; d < 8; ++d) {
	             auto algosumvalue = 0;
	             for (int i = 0; i < tpSamples; ++i)
	                 algosumvalue += int(theDepthMap[samples.id()][d][idx + i]);
	             depth_sums[d] += std::min<unsigned int>(algosumvalue, QIE11_MAX_LINEARIZATION_ET);
	         }
         }  
     } else {
         // Not a peak
         output[ibin] = 0;
         
         // See comment above.
         if (ibin == numberOfPresamples_) {
	         for (int d = 0; d < 8; ++d) {
	             auto algosumvalue = 0;
	             for (int i = 0; i < tpSamples; ++i)
	                 algosumvalue += int(theDepthMap[samples.id()][d][idx + i]);
	             depth_sums[d] += std::min<unsigned int>(algosumvalue, QIE11_MAX_LINEARIZATION_ET);
	         }
         }
     }
     // peak-finding is not applied for FG bits
     finegrain[ibin] = fg_algo.compute(msb[idx]).to_ulong();
   }

   outcoder_->compress(output, finegrain, result);
   update(result, theDepthMap[samples.id()], numberOfPresamples_,depth_sums, thesampledata); 
}

void HcalTriggerPrimitiveAlgo::analyzeHF(IntegerCaloSamples & samples, HcalTriggerPrimitiveDigi & result, const int hf_lumi_shift) {
   HcalTrigTowerDetId detId(samples.id());

   // Align digis and TP
   int dgPresamples=samples.presamples(); 
   int tpPresamples=numberOfPresamplesHF_;
   int shift = dgPresamples - tpPresamples;
   int dgSamples=samples.size();
   int tpSamples=numberOfSamplesHF_;
   if(shift<0 || shift+tpSamples>dgSamples){
	edm::LogInfo("HcalTriggerPrimitiveAlgo::analyzeHF") << 
	    "TP presample or size from the configuration file is out of the accessible range. Using digi values from data instead...";
	tpPresamples=dgPresamples;
	shift=0;
	tpSamples=dgSamples;
   }

   std::vector<int> finegrain(tpSamples, false);

   TowerMapFGSum::const_iterator tower2fg = theTowerMapFGSum.find(detId);
   assert(tower2fg != theTowerMapFGSum.end());

   const SumFGContainer& sumFG = tower2fg->second;
   // Loop over all L+S pairs that mapped from samples.id()
   // Note: 1 samples.id() = 6 x (L+S) without noZS
   for (SumFGContainer::const_iterator sumFGItr = sumFG.begin(); sumFGItr != sumFG.end(); ++sumFGItr) {
      const std::vector<bool>& veto = HF_Veto[sumFGItr->id().rawId()];
      for (int ibin = 0; ibin < tpSamples; ++ibin) {
         int idx = ibin + shift;
         // if not vetod, add L+S to total sum and calculate FG
	 bool vetoed = idx<int(veto.size()) && veto[idx];
         if (!(vetoed && (*sumFGItr)[idx] > PMT_NoiseThreshold_)) {
            samples[idx] += (*sumFGItr)[idx];
            finegrain[ibin] = (finegrain[ibin] || (*sumFGItr)[idx] >= FG_threshold_);
         }
      }
   }

   IntegerCaloSamples output(samples.id(), tpSamples);
   output.setPresamples(tpPresamples);

   for (int ibin = 0; ibin < tpSamples; ++ibin) {
      int idx = ibin + shift;
      output[ibin] = samples[idx] >> hf_lumi_shift;
      static const int MAX_OUTPUT = QIE8_LINEARIZATION_ET;  // QIE8_LINEARIZATION_ET = 1023
      if (output[ibin] > MAX_OUTPUT) output[ibin] = MAX_OUTPUT;
   }
   outcoder_->compress(output, finegrain, result);
}

void HcalTriggerPrimitiveAlgo::analyzeHF2016(
        const IntegerCaloSamples& samples,
        HcalTriggerPrimitiveDigi& result,
        const int hf_lumi_shift,
        const HcalFeatureBit* embit
        ) {
    // Align digis and TP
    const int SHIFT = samples.presamples() - numberOfPresamplesHF_;
    assert(SHIFT >= 0);
    assert((SHIFT + numberOfSamplesHF_) <= samples.size());

    // Try to find the HFDetails from the map corresponding to our samples
    const HcalTrigTowerDetId detId(samples.id());
    HFDetailMap::const_iterator it = theHFDetailMap.find(detId);
    // Missing values will give an empty digi
    if (it == theHFDetailMap.end()) {
        return;
    }

    std::vector<std::bitset<2>> finegrain(numberOfSamplesHF_, false);

    // Set up out output of IntergerCaloSamples
    IntegerCaloSamples output(samples.id(), numberOfSamplesHF_);
    output.setPresamples(numberOfPresamplesHF_);

    for (const auto& item: it->second) {
        auto& details = item.second;
        for (int ibin = 0; ibin < numberOfSamplesHF_; ++ibin) {
            const int IDX = ibin + SHIFT;
            int long_fiber_val = 0;
            if (IDX < details.long_fiber.size()) {
                long_fiber_val = details.long_fiber[IDX];
            }
            int short_fiber_val = 0;
            if (IDX < details.short_fiber.size()) {
                short_fiber_val = details.short_fiber[IDX];
            }
            output[ibin] += (long_fiber_val + short_fiber_val);

            uint32_t ADCLong = details.LongDigi[ibin].adc();
            uint32_t ADCShort = details.ShortDigi[ibin].adc();

            if (details.LongDigi.id().ietaAbs() >= FIRST_FINEGRAIN_TOWER) {
               finegrain[ibin][1] = (ADCLong > FG_HF_thresholds_[0] || ADCShort > FG_HF_thresholds_[0]);

               if (embit != nullptr)
                  finegrain[ibin][0] = embit->fineGrainbit(details.ShortDigi, details.LongDigi, ibin);
            }
        }
    }

    for (int bin = 0; bin < numberOfSamplesHF_; ++bin) {
       static const unsigned int MAX_OUTPUT = QIE8_LINEARIZATION_ET;  // QIE8_LINEARIZATION_ET = 1023
       output[bin] = min({MAX_OUTPUT, output[bin] >> hf_lumi_shift});
    }

    std::vector<int> finegrain_converted;
    for (const auto& fg: finegrain)
       finegrain_converted.push_back(fg.to_ulong());
    outcoder_->compress(output, finegrain_converted, result);
}

bool
HcalTriggerPrimitiveAlgo::passTDC(const QIE10DataFrame& digi, int ts) const
{
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

bool
HcalTriggerPrimitiveAlgo::validChannel(const QIE10DataFrame& digi, int ts) const
{
   // channels with invalid data should not contribute to the sum
   if(digi.linkError() || ts>=digi.samples() || !digi[ts].ok()) return false;

   auto mask = conditions_->getHcalTPChannelParameter(HcalDetId(digi.id()))->getMask();
   if (mask)
      return false;

   return true;
}

void HcalTriggerPrimitiveAlgo::analyzeHFQIE10(
        const IntegerCaloSamples& samples, HcalTriggerPrimitiveDigi& result,
        const int hf_lumi_shift, const HcalFeatureBit* embit)
{
    // Align digis and TP
    const int shift = samples.presamples() - numberOfPresamplesHF_;
    assert(shift >= 0);
    assert((shift + numberOfSamplesHF_) <= samples.size());
    assert(hf_lumi_shift>=2);

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

    for (const auto& item: it->second) {
        auto& details = item.second;
        for (int ibin = 0; ibin < numberOfSamplesHF_; ++ibin) {
            const int idx = ibin + shift;

            int long_fiber_val = 0;
            int long_fiber_count = 0;
            int short_fiber_val = 0;
            int short_fiber_count = 0;

            bool saturated = false;

            for (auto i: {0, 2}) {
               if (idx < details[i].samples.size() and details[i].validity[idx] and details[i].passTDC[idx]) {
                  long_fiber_val += details[i].samples[idx];
                  saturated = saturated || (details[i].samples[idx] == QIE10_LINEARIZATION_ET);
                  ++long_fiber_count;
               }
            }
            for (auto i: {1, 3}) {
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
		if (long_fiber_count == 2) long_fiber_val >>=1;
		if (short_fiber_count == 2) short_fiber_val >>=1;
		
		auto sum = long_fiber_val + short_fiber_val;
		// Similar to above, if both channels are valid, 
		// average of the two energies is calculated
		// division by 2 here is also compensated by adjusting the total scale shift in the end
		if (long_fiber_count > 0 and short_fiber_count > 0) sum >>=1;

		output[ibin] += sum;
            }

            for (const auto& detail: details) {
               if (idx < int(detail.digi.size()) and detail.validity[idx] and HcalDetId(detail.digi.id()).ietaAbs() >= FIRST_FINEGRAIN_TOWER) {
		 if(useTDCInMinBiasBits_ && !detail.passTDC[idx]) continue;
		 finegrain[ibin][1] = finegrain[ibin][1] or detail.fgbits[idx][0];
		 // what is commonly called the "second" HF min-bias bit is
		 // actually the 0-th bit, which can also be used instead for the EM bit
		 // (called finegrain[ibin][0] below) in non-HI running
		 finegrain[ibin][0] = finegrain[ibin][0] or detail.fgbits[idx][1];
               }
            }
	    // the EM bit is only used if the "second" FG bit is disabled
            if (embit != nullptr and FG_HF_thresholds_.at(1) != 255) {
               finegrain[ibin][0] = embit->fineGrainbit(
                     details[1].digi, details[3].digi,
                     details[0].digi, details[2].digi,
                     details[1].validity[idx], details[3].validity[idx],
                     details[0].validity[idx], details[2].validity[idx],
                     idx
               );
            }
        }
    }

    for (int bin = 0; bin < numberOfSamplesHF_; ++bin) {
       output[bin] = min({(unsigned int) QIE10_MAX_LINEARIZATION_ET, output[bin] >> (hf_lumi_shift-2)});
    }
    std::vector<int> finegrain_converted;
    for (const auto& fg: finegrain)
       finegrain_converted.push_back(fg.to_ulong());
    outcoder_->compress(output, finegrain_converted, result);
}

void HcalTriggerPrimitiveAlgo::analyzeHFQIE10(
        const IntegerCaloSamples& samples, HcalUpgradeTriggerPrimitiveDigi& result,
        const int hf_lumi_shift, const HcalFeatureBit* embit)
{
    // Align digis and TP
    const int shift = samples.presamples() - numberOfPresamplesHF_;
    assert(shift >= 0);
    assert((shift + numberOfSamplesHF_) <= samples.size());
    assert(hf_lumi_shift>=2);

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

    for (const auto& item: it->second) {
        auto& details = item.second;
        for (int ibin = 0; ibin < numberOfSamplesHF_; ++ibin) {
            const int idx = ibin + shift;

            int long_fiber_val = 0;
            int long_fiber_count = 0;
            int short_fiber_val = 0;
            int short_fiber_count = 0;

            bool saturated = false;

            for (auto i: {0, 2}) {
               if (idx < details[i].samples.size() and details[i].validity[idx] and details[i].passTDC[idx]) {
                  long_fiber_val += details[i].samples[idx];
                  saturated = saturated || (details[i].samples[idx] == QIE10_LINEARIZATION_ET);
                  ++long_fiber_count;
               }
            }
            for (auto i: {1, 3}) {
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
		if (long_fiber_count == 2) long_fiber_val >>=1;
		if (short_fiber_count == 2) short_fiber_val >>=1;
		
		auto sum = long_fiber_val + short_fiber_val;
		// Similar to above, if both channels are valid, 
		// average of the two energies is calculated
		// division by 2 here is also compensated by adjusting the total scale shift in the end
		if (long_fiber_count > 0 and short_fiber_count > 0) sum >>=1;

		output[ibin] += sum;
            }

            for (const auto& detail: details) {
               if (idx < int(detail.digi.size()) and detail.validity[idx] and HcalDetId(detail.digi.id()).ietaAbs() >= FIRST_FINEGRAIN_TOWER) {
		 if(useTDCInMinBiasBits_ && !detail.passTDC[idx]) continue;
		 finegrain[ibin][1] = finegrain[ibin][1] or detail.fgbits[idx][0];
		 // what is commonly called the "second" HF min-bias bit is
		 // actually the 0-th bit, which can also be used instead for the EM bit
		 // (called finegrain[ibin][0] below) in non-HI running
		 finegrain[ibin][0] = finegrain[ibin][0] or detail.fgbits[idx][1];
               }
            }
	    // the EM bit is only used if the "second" FG bit is disabled
            if (embit != nullptr and FG_HF_thresholds_.at(1) != 255) {
               finegrain[ibin][0] = embit->fineGrainbit(
                     details[1].digi, details[3].digi,
                     details[0].digi, details[2].digi,
                     details[1].validity[idx], details[3].validity[idx],
                     details[0].validity[idx], details[2].validity[idx],
                     idx
               );
            }
        }
    }

    for (int bin = 0; bin < numberOfSamplesHF_; ++bin) {
       output[bin] = min({(unsigned int) QIE10_MAX_LINEARIZATION_ET, output[bin] >> (hf_lumi_shift-2)});
    }
    std::vector<int> finegrain_converted;
    for (const auto& fg: finegrain)
       finegrain_converted.push_back(fg.to_ulong());
    outcoder_->compress(output, finegrain_converted, result);
}

void HcalTriggerPrimitiveAlgo::runZS(HcalTrigPrimDigiCollection & result){
   for (HcalTrigPrimDigiCollection::iterator tp = result.begin(); tp != result.end(); ++tp){
      bool ZS = true;
      for (int i=0; i<tp->size(); ++i) {
         if (tp->sample(i).compressedEt()  > ZS_threshold_I_) {
            ZS=false;
            break;
         }
      }
      if (ZS) tp->setZSInfo(false,true);
      else tp->setZSInfo(true,false);
   }
}

void HcalTriggerPrimitiveAlgo::runFEFormatError(const FEDRawDataCollection* rawraw,
                                                const HcalElectronicsMap *emap,
                                                HcalTrigPrimDigiCollection & result
                                                ){
  std::set<uint32_t> FrontEndErrors;

  for(int i=FEDNumbering::MINHCALFEDID; i<=FEDNumbering::MAXHCALFEDID; ++i) {
    const FEDRawData& raw = rawraw->FEDData(i);
    if (raw.size()<12) continue;
    const HcalDCCHeader* dccHeader=(const HcalDCCHeader*)(raw.data());
    if(!dccHeader) continue;
    HcalHTRData htr;
    for (int spigot=0; spigot<HcalDCCHeader::SPIGOT_COUNT; spigot++) {
      if (!dccHeader->getSpigotPresent(spigot)) continue;
      dccHeader->getSpigotData(spigot,htr,raw.size());
      int dccid = dccHeader->getSourceId();
      int errWord = htr.getErrorsWord() & 0x1FFFF;
      bool HTRError = (!htr.check() || htr.isHistogramEvent() || (errWord & 0x800)!=0);

      if(HTRError) {
        bool valid =false;
        for(int fchan=0; fchan<3 && !valid; fchan++) {
          for(int fib=0; fib<9 && !valid; fib++) {
            HcalElectronicsId eid(fchan,fib,spigot,dccid-FEDNumbering::MINHCALFEDID);
            eid.setHTR(htr.readoutVMECrateId(),htr.htrSlot(),htr.htrTopBottom());
            DetId detId = emap->lookup(eid);
            if(detId.null()) continue;
            HcalSubdetector subdet=(HcalSubdetector(detId.subdetId()));
            if (detId.det()!=4||
              (subdet!=HcalBarrel && subdet!=HcalEndcap &&
              subdet!=HcalForward )) continue;
            std::vector<HcalTrigTowerDetId> ids = theTrigTowerGeometry->towerIds(detId);
            for (std::vector<HcalTrigTowerDetId>::const_iterator triggerId=ids.begin(); triggerId != ids.end(); ++triggerId) {
              FrontEndErrors.insert(triggerId->rawId());
            }
            //valid = true;
          }
        }
      }
    }
  }

  // Loop over TP collection
  // Set TP to zero if there is FE Format Error
  HcalTriggerPrimitiveSample zeroSample(0);
  for (HcalTrigPrimDigiCollection::iterator tp = result.begin(); tp != result.end(); ++tp){
    if (FrontEndErrors.find(tp->id().rawId()) != FrontEndErrors.end()) {
      for (int i=0; i<tp->size(); ++i) tp->setSample(i, zeroSample);
    }
  }
}

void HcalTriggerPrimitiveAlgo::addFG(const HcalTrigTowerDetId& id, std::vector<bool>& msb){
   FGbitMap::iterator itr = fgMap_.find(id);
   if (itr != fgMap_.end()){
      std::vector<bool>& _msb = itr->second;
      for (size_t i=0; i<msb.size(); ++i)
         _msb[i] = _msb[i] || msb[i];
   }
   else fgMap_[id] = msb;
}

bool
HcalTriggerPrimitiveAlgo::validUpgradeFG(const HcalTrigTowerDetId& id, int depth) const
{
   if (depth > LAST_FINEGRAIN_DEPTH)
      return false;
   if (id.ietaAbs() > LAST_FINEGRAIN_TOWER)
      return false;
   if (id.ietaAbs() == HBHE_OVERLAP_TOWER and not upgrade_hb_)
      return false;
   return true;
}

bool
HcalTriggerPrimitiveAlgo::needLegacyFG(const HcalTrigTowerDetId& id) const
{
   // This tower (ietaAbs == 16) does not accept upgraded FG bits,
   // but needs pseudo legacy ones to ensure that the tower is processed
   // even when the QIE8 depths in front of it do not have energy deposits.
   if (id.ietaAbs() == HBHE_OVERLAP_TOWER and not upgrade_hb_)
      return true;
   return false;
}

void
HcalTriggerPrimitiveAlgo::addUpgradeFG(const HcalTrigTowerDetId& id, int depth, const std::vector<std::bitset<2>>& bits)
{
   if (not validUpgradeFG(id, depth)) {
      if (needLegacyFG(id)) {
         std::vector<bool> pseudo(bits.size(), false);
         addFG(id, pseudo);
      }
      return;
   }

   auto it = fgUpgradeMap_.find(id);
   if (it == fgUpgradeMap_.end()) {
      FGUpgradeContainer element;
      element.resize(bits.size());
      it = fgUpgradeMap_.insert(std::make_pair(id, element)).first;
   }
   for (unsigned int i = 0; i < bits.size(); ++i) {
      it->second[i][0][depth-1] = bits[i][0];
      it->second[i][1][depth-1] = bits[i][1];
   }
}

void HcalTriggerPrimitiveAlgo::setPeakFinderAlgorithm(int algo){
   if (algo <=0 && algo>2)
      throw cms::Exception("ERROR: Only algo 1 & 2 are supported.") << std::endl;
   peak_finder_algorithm_ = algo;
}

void HcalTriggerPrimitiveAlgo::setPeakFinderAlgorithm(std::string algo) {

    peak_finder_algorithm_name_ = algo;
    if (peak_finder_algorithm_name_ == "PFA1p_AVE") {
        // average weights when fixing w3 = 1 (1 PS) 
        weights_ =
        {{1, {-0.42}},
        {2,  {-0.42}},
        {3,  {-0.42}},
        {4,  {-0.42}},
        {5,  {-0.42}},
        {6,  {-0.42}},
        {7,  {-0.42}},
        {8,  {-0.42}},
        {9,  {-0.42}},
        {10, {-0.42}},
        {11, {-0.42}},
        {12, {-0.42}},
        {13, {-0.42}},
        {14, {-0.42}},
        {15, {-0.42}},
        {16, {-0.42}},
        {17, {-0.45}},
        {18, {-0.45}},
        {19, {-0.45}},
        {20, {-0.45}},
        {21, {-0.83}},
        {22, {-0.83}},
        {23, {-0.83}},
        {24, {-0.83}},
        {25, {-0.83}},
        {26, {-0.83}},
        {27, {-0.83}},
        {28, {-0.83}}};

    } else if (peak_finder_algorithm_name_ == "PFA1pp_AVE") {
        // per-eta weights when fixing w3 = 1 (2 PS) 
        weights_ =
        {{1, {0.21,-0.61}},
        {2,  {0.21,-0.61}},
        {3,  {0.21,-0.61}},
        {4,  {0.21,-0.61}},
        {5,  {0.21,-0.61}},
        {6,  {0.21,-0.61}},
        {7,  {0.21,-0.61}},
        {8,  {0.21,-0.61}},
        {9,  {0.21,-0.61}},
        {10, {0.21,-0.61}},
        {11, {0.21,-0.61}},
        {12, {0.21,-0.61}},
        {13, {0.21,-0.61}},
        {14, {0.21,-0.61}},
        {15, {0.21,-0.61}},
        {16, {0.21,-0.61}},
        {17, {0.20,-0.59}},
        {18, {0.20,-0.59}},
        {19, {0.20,-0.59}},
        {20, {0.20,-0.59}},
        {21, {-0.07,-0.56}},
        {22, {-0.07,-0.56}},
        {23, {-0.07,-0.56}},
        {24, {-0.07,-0.56}},
        {25, {-0.07,-0.56}},
        {26, {-0.07,-0.56}},
        {27, {-0.07,-0.56}},
        {28, {-0.07,-0.56}}};

    } else if (peak_finder_algorithm_name_ == "PFA2p") {
    
        // per-ieta weights when fixing w3 = w4 = 1.0 (2 PS)
        weights_ =
        {{1, {-0.09, -0.39}},
         {2, {-0.09, -0.38}},
         {3, {-0.10, -0.38}},
         {4, {-0.13, -0.39}},
         {5, {-0.09, -0.39}},
         {6, {-0.10, -0.39}},
         {7, {-0.08, -0.37}},
         {8, {-0.09, -0.38}},
         {9, {-0.11, -0.38}},
         {10,{ -0.09, -0.40}},
         {11,{ -0.12, -0.39}},
         {12,{ -0.09, -0.39}},
         {13,{ -0.11, -0.38}},
         {14,{ -0.13, -0.40}},
         {15,{ -0.13, -0.41}},
         {16,{ -0.12, -0.40}},
         {17,{ -0.07, -0.42}},
         {18,{ -0.11, -0.40}},
         {19,{ -0.11, -0.42}},
         {20,{ -0.09, -0.43}},
         {21,{ -0.18, -0.83}},
         {22,{ -0.22, -0.85}},
         {23,{ -0.24, -0.88}},
         {24,{ -0.28, -0.91}},
         {25,{ -0.32, -0.94}},
         {26,{ -0.34, -0.96}},
         {27,{ -0.34, -0.96}},
         {28,{ -0.40, -1.04}}};

    } else if (peak_finder_algorithm_name_ == "PFA2p_AVE") {
        // Average HBHE weights when fixing w3 = w4 = 1.0 
        weights_ =
        {{1, {-0.11,-0.39}},
        {2,  {-0.11,-0.39}},
        {3,  {-0.11,-0.39}},
        {4,  {-0.11,-0.39}},
        {5,  {-0.11,-0.39}},
        {6,  {-0.11,-0.39}},
        {7,  {-0.11,-0.39}},
        {8,  {-0.11,-0.39}},
        {9,  {-0.11,-0.39}},
        {10, {-0.11,-0.39}},
        {11, {-0.11,-0.39}},
        {12, {-0.11,-0.39}},
        {13, {-0.11,-0.39}},
        {14, {-0.11,-0.39}},
        {15, {-0.11,-0.39}},
        {16, {-0.11,-0.39}},
        {17, {-0.10,-0.42}},
        {18, {-0.10,-0.42}},
        {19, {-0.10,-0.42}},
        {20, {-0.10,-0.42}},
        {21, {-0.32,-0.95}},
        {22, {-0.32,-0.95}},
        {23, {-0.32,-0.95}},
        {24, {-0.32,-0.95}},
        {25, {-0.32,-0.95}},
        {26, {-0.32,-0.95}},
        {27, {-0.32,-0.95}},
        {28, {-0.32,-0.95}}};

    } else if (peak_finder_algorithm_name_ == "PFA3") {
        // HBHE weights when fixing w3 = w4 = 1.0 
        weights_ =
        {{1, {-2.0}},
        {2,  {-2.0}},
        {3,  {-2.0}},
        {4,  {-2.0}},
        {5,  {-2.0}},
        {6,  {-2.0}},
        {7,  {-2.0}},
        {8,  {-2.0}},
        {9,  {-2.0}},
        {10, {-2.0}},
        {11, {-2.0}},
        {12, {-2.0}},
        {13, {-2.0}},
        {14, {-2.0}},
        {15, {-2.0}},
        {16, {-2.0}},
        {17, {-2.0}},
        {18, {-2.0}},
        {19, {-2.0}},
        {20, {-2.0}},
        {21, {-2.0}},
        {22, {-2.0}},
        {23, {-2.0}},
        {24, {-2.0}},
        {25, {-2.0}},
        {26, {-2.0}},
        {27, {-2.0}},
        {28, {-2.0}}};

    } else if (peak_finder_algorithm_name_ == "PFA3p") {
        // per-ieta weights when fixing w1 = w2 = 1 (1 PS)
        weights_ =
        {{1, {-0.42}},
        {2,  {-0.42}},
        {3,  {-0.41}},
        {4,  {-0.44}},
        {5,  {-0.46}},
        {6,  {-0.44}},
        {7,  {-0.45}},
        {8,  {-0.45}},
        {9,  {-0.46}},
        {10, {-0.45}},
        {11, {-0.48}},
        {12, {-0.47}},
        {13, {-0.50}},
        {14, {-0.51}},
        {15, {-0.54}},
        {16, {-0.54}},
        {17, {-0.52}},
        {18, {-0.53}},
        {19, {-0.51}},
        {20, {-0.52}},
        {21, {-0.63}},
        {22, {-0.69}},
        {23, {-0.78}},
        {24, {-0.89}},
        {25, {-1.10}},
        {26, {-1.27}},
        {27, {-1.22}},
        {28, {-1.75}}};

    } else if (peak_finder_algorithm_name_ == "PFA3p_AVE") {
        // Average weights when fixing w1 = w2 = 1 (1 PS)
        weights_ =
        {{1, {-0.48}},
        {2,  {-0.48}},
        {3,  {-0.48}},
        {4,  {-0.48}},
        {5,  {-0.48}},
        {6,  {-0.48}},
        {7,  {-0.48}},
        {8,  {-0.48}},
        {9,  {-0.48}},
        {10, {-0.48}},
        {11, {-0.48}},
        {12, {-0.48}},
        {13, {-0.48}},
        {14, {-0.48}},
        {15, {-0.48}},
        {16, {-0.48}},
        {17, {-0.52}},
        {18, {-0.52}},
        {19, {-0.52}},
        {20, {-0.52}},
        {21, {-1.19}},
        {22, {-1.19}},
        {23, {-1.19}},
        {24, {-1.19}},
        {25, {-1.19}},
        {26, {-1.19}},
        {27, {-1.19}},
        {28, {-1.19}}};
    }
}

void HcalTriggerPrimitiveAlgo::setNCTScaleShift(int shift){
   NCTScaleShift = shift;
}

void HcalTriggerPrimitiveAlgo::setRCTScaleShift(int shift){
   RCTScaleShift = shift;
}

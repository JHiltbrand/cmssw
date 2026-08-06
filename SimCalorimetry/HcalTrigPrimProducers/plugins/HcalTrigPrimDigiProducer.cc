#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/HcalTPGCompressor.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGCoder.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGRecord.h"
#include "CondFormats/HcalObjects/interface/HcalElectronicsMap.h"
#include "CondFormats/HcalObjects/interface/HcalLutMetadata.h"
#include "CondFormats/HcalObjects/interface/HcalTPChannelParameters.h"
#include "CondFormats/DataRecord/interface/HcalLutMetadataRcd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "SimCalorimetry/HcalTrigPrimAlgos/interface/HcalTriggerPrimitiveAlgo.h"

#include <algorithm>
#include <vector>

class HcalTrigPrimDigiProducer : public edm::stream::EDProducer<edm::stream::WatchRuns> {
public:
  explicit HcalTrigPrimDigiProducer(const edm::ParameterSet& ps);

  /**Produces the EDM products,*/
  void beginRun(const edm::Run& r, const edm::EventSetup& c) override;
  void produce(edm::Event& e, const edm::EventSetup& c) override;

private:
  HcalTriggerPrimitiveAlgo theAlgo_;

  /// input tags for HCAL digis
  std::vector<edm::InputTag> inputUpgradeLabel_;
  // this seems a strange way of doing things
  edm::EDGetTokenT<QIE11DigiCollection> tok_hbhe_up_;
  edm::EDGetTokenT<QIE10DigiCollection> tok_hf_up_;

  bool overrideDBvetoThresholdsHB_;
  bool overrideDBweightsAndFilterHB_;

  double MinLongEnergy_, MinShortEnergy_, LongShortSlope_, LongShortOffset_;

  bool HFEMB_;
  edm::ParameterSet LongShortCut_;
  edm::ESGetToken<HcalTPGCoder, HcalTPGRecord> tok_tpgCoder_;
  edm::ESGetToken<CaloTPGTranscoder, CaloTPGRecord> tok_tpgTranscoder_;
  edm::ESGetToken<HcalLutMetadata, HcalLutMetadataRcd> tok_lutMetadata_;
  edm::ESGetToken<HcalTrigTowerGeometry, CaloGeometryRecord> tok_trigTowerGeom_;
  edm::ESGetToken<HcalTopology, HcalRecNumberingRecord> tok_hcalTopo_;
  edm::ESGetToken<HcalDbService, HcalDbRecord> tok_dbService_;
  edm::ESGetToken<HcalDbService, HcalDbRecord> tok_dbService_beginRun_;
};

HcalTrigPrimDigiProducer::HcalTrigPrimDigiProducer(const edm::ParameterSet& ps)
    : theAlgo_(ps.getParameter<std::vector<uint32_t> >("FG_HF_thresholds"),
               ps.getParameter<bool>("useTDCfromDB"),
               ps.getParameter<int>("numberOfSamples"),
               ps.getParameter<int>("numberOfPresamples"),
               ps.getParameter<int>("numberOfFilterPresamplesHBQIE11"),
               ps.getParameter<int>("numberOfSamplesHF"),
               ps.getParameter<int>("numberOfPresamplesHF"),
               ps.getParameter<bool>("useTDCInMinBiasBits")),
      inputUpgradeLabel_(ps.getParameter<std::vector<edm::InputTag> >("inputUpgradeLabel")) {
  overrideDBvetoThresholdsHB_ = ps.getParameter<bool>("overrideDBvetoThresholdsHB");
  overrideDBweightsAndFilterHB_ = ps.getParameter<bool>("overrideDBweightsAndFilterHB");

  theAlgo_.setWeightsQIE11(ps.getParameter<edm::ParameterSet>("weightsQIE11"));
  theAlgo_.setCodedVetoThresholds(ps.getParameter<edm::ParameterSet>("codedVetoThresholds"));

  if (ps.exists("parameters")) {
    auto pset = ps.getUntrackedParameter<edm::ParameterSet>("parameters");
    theAlgo_.overrideParameters(pset);
  }
  theAlgo_.setFixSaturationFlag(ps.getParameter<bool>("applySaturationFix"));

  HFEMB_ = false;
  if (ps.exists("LSConfig")) {
    LongShortCut_ = ps.getUntrackedParameter<edm::ParameterSet>("LSConfig");
    HFEMB_ = LongShortCut_.getParameter<bool>("HcalFeatureHFEMBit");
    MinLongEnergy_ = LongShortCut_.getParameter<double>("Min_Long_Energy");    //minimum long energy
    MinShortEnergy_ = LongShortCut_.getParameter<double>("Min_Short_Energy");  //minimum short energy
    LongShortSlope_ =
        LongShortCut_.getParameter<double>("Long_vrs_Short_Slope");  //slope of the line that cuts are based on
    LongShortOffset_ = LongShortCut_.getParameter<double>("Long_Short_Offset");  //offset of line
  }
  tok_tpgCoder_ = esConsumes<HcalTPGCoder, HcalTPGRecord>();
  tok_tpgTranscoder_ = esConsumes<CaloTPGTranscoder, CaloTPGRecord>();
  tok_lutMetadata_ = esConsumes<HcalLutMetadata, HcalLutMetadataRcd>();
  tok_trigTowerGeom_ = esConsumes<HcalTrigTowerGeometry, CaloGeometryRecord>();
  tok_hcalTopo_ = esConsumes<HcalTopology, HcalRecNumberingRecord, edm::Transition::BeginRun>();

  tok_hbhe_up_ = consumes<QIE11DigiCollection>(inputUpgradeLabel_[0]);
  tok_hf_up_ = consumes<QIE10DigiCollection>(inputUpgradeLabel_[1]);

  tok_dbService_ = esConsumes<HcalDbService, HcalDbRecord>();
  tok_dbService_beginRun_ = esConsumes<HcalDbService, HcalDbRecord, edm::Transition::BeginRun>();

  produces<HcalTrigPrimDigiCollection>();

  edm::ParameterSet hfSS = ps.getParameter<edm::ParameterSet>("tpScales").getParameter<edm::ParameterSet>("HF");

  theAlgo_.setNCTScaleShift(hfSS.getParameter<int>("NCTShift"));
  theAlgo_.setRCTScaleShift(hfSS.getParameter<int>("RCTShift"));
}

void HcalTrigPrimDigiProducer::beginRun(const edm::Run& run, const edm::EventSetup& eventSetup) {
  edm::ESHandle<HcalDbService> db = eventSetup.getHandle(tok_dbService_beginRun_);
  const HcalTopology* topo = &eventSetup.getData(tok_hcalTopo_);

  const HcalElectronicsMap* emap = db->getHcalMapping();

  int lastHBRing = topo->lastHBRing();

  std::vector<HcalElectronicsId> vIds = emap->allElectronicsIdTrigger();
  for (std::vector<HcalElectronicsId>::const_iterator eId = vIds.begin(); eId != vIds.end(); eId++) {
    HcalTrigTowerDetId hcalTTDetId(emap->lookupTrigger(*eId));
    if (hcalTTDetId.null())
      continue;

    int aieta = std::abs(hcalTTDetId.ieta());

    // Filter weight represented in fixed point 8 bit
    int fixedPointWeight = 255;
    // Coded veto threshold in range (0, 2048)
    // Default, special value of 0 will disable vetoing
    int codedVetoThreshold = 0;
    // Number of filter presamples
    int presamples = 0;

    // The absence of TT channels in the HcalTPChannelParameters
    // is intepreted as to not use the new filter
    auto tpParam = db->getHcalTPChannelParameter(hcalTTDetId, false);
    if (tpParam) {
      // Fix number of filter presamples to one if we are using DB weights
      // Size of filter is already known when using DB weights
      // Weight from DB represented as 8-bit integer
      fixedPointWeight = tpParam->getauxi1();
      codedVetoThreshold = tpParam->getauxi2();
      presamples = 1;
    }

    // If the aieta already has a weight in the map, then move on
    if (aieta <= lastHBRing) {
      if (!overrideDBvetoThresholdsHB_) {
        theAlgo_.setCodedVetoThreshold(aieta, codedVetoThreshold);
      }
      if (!overrideDBweightsAndFilterHB_) {
        theAlgo_.setNumFilterPresamplesHBQIE11(presamples);
        theAlgo_.setWeightQIE11(aieta, fixedPointWeight);
      }
    }
  }
}

void HcalTrigPrimDigiProducer::produce(edm::Event& iEvent, const edm::EventSetup& eventSetup) {
  edm::ESHandle<HcalTPGCoder> inputCoder = eventSetup.getHandle(tok_tpgCoder_);

  edm::ESHandle<CaloTPGTranscoder> outTranscoder = eventSetup.getHandle(tok_tpgTranscoder_);

  edm::ESHandle<HcalLutMetadata> lutMetadata = eventSetup.getHandle(tok_lutMetadata_);
  float rctlsb = lutMetadata->getRctLsb();

  edm::ESHandle<HcalTrigTowerGeometry> pG = eventSetup.getHandle(tok_trigTowerGeom_);

  std::unique_ptr<HcalTrigPrimDigiCollection> result(new HcalTrigPrimDigiCollection());

  edm::Handle<QIE11DigiCollection> hbheUpDigis;
  edm::Handle<QIE10DigiCollection> hfUpDigis;

  iEvent.getByToken(tok_hbhe_up_, hbheUpDigis);
  iEvent.getByToken(tok_hf_up_, hfUpDigis);

  if (!hbheUpDigis.isValid()) {
    edm::LogInfo("HcalTrigPrimDigiProducer")
        << "\nWarning: Upgrade HBHEDigiCollection with input tag " << inputUpgradeLabel_[0]
        << "\nrequested in configuration, but not found in the event."
        << "\nQuit returning empty product." << std::endl;

    // put empty HcalTrigPrimDigiCollection in the event
    iEvent.put(std::move(result));

    return;
  }

  if (!hfUpDigis.isValid()) {
    edm::LogInfo("HcalTrigPrimDigiProducer") << "\nWarning: HFDigiCollection with input tag " << inputUpgradeLabel_[1]
                                             << "\nrequested in configuration, but not found in the event."
                                             << "\nQuit returning empty product." << std::endl;

    // put empty HcalTrigPrimDigiCollection in the event
    iEvent.put(std::move(result));

    return;
  }

  edm::ESHandle<HcalDbService> pSetup = eventSetup.getHandle(tok_dbService_);

  HcalFeatureBit* hfembit = nullptr;

  if (HFEMB_) {
    hfembit = new HcalFeatureHFEMBit(MinShortEnergy_,
                                     MinLongEnergy_,
                                     LongShortSlope_,
                                     LongShortOffset_,
                                     *pSetup);  //inputs values that cut will be based on
  }

  theAlgo_.run(inputCoder.product(),
               outTranscoder->getHcalCompressor().get(),
               pSetup.product(),
               *result,
               &(*pG),
               rctlsb,
               hfembit,
               *hbheUpDigis,
               *hfUpDigis);

  edm::LogInfo("HcalTrigPrimDigiProducer") << "HcalTrigPrims: " << result->size();

  iEvent.put(std::move(result));
}

DEFINE_FWK_MODULE(HcalTrigPrimDigiProducer);

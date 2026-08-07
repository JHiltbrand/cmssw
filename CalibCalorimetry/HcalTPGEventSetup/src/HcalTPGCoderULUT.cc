// -*- C++ -*-
//
// Package:    HcalTPGCoderULUT
// Class:      HcalTPGCoderULUT
//
/**\class HcalTPGCoderULUT HcalTPGCoderULUT.h src/HcalTPGCoderULUT/interface/HcalTPGCoderULUT.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jeremiah Mans
//         Created:  Fri Sep 15 11:49:44 CDT 2006
//
//

// system include files
#include <memory>
#include <string>

// user include files

#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ESProductHost.h"
#include "FWCore/Utilities/interface/ReusableObjectHolder.h"

#include "CalibCalorimetry/HcalTPGAlgos/interface/HcaluLUTTPGCoder.h"
#include "CalibFormats/HcalObjects/interface/HcalTPGRecord.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CondFormats/HcalObjects/interface/HcalElectronicsMap.h"

class HcalTPGCoderULUT : public edm::ESProducer {
public:
  HcalTPGCoderULUT(const edm::ParameterSet&);
  ~HcalTPGCoderULUT() override;

  typedef std::shared_ptr<HcalTPGCoder> ReturnType;

  ReturnType produce(const HcalTPGRecord&);

private:
  using HostType = edm::ESProductHost<HcaluLUTTPGCoder, HcalDbRecord>;

  void buildCoder(const HcalTopology*, const HcalElectronicsMap*, const HcalTimeSlew*, HcaluLUTTPGCoder*);

  edm::ReusableObjectHolder<HostType> holder_;
  edm::ESGetToken<HcalTopology, HcalRecNumberingRecord> topoToken_;
  edm::ESGetToken<HcalTimeSlew, HcalTimeSlewRecord> delayToken_;
  edm::ESGetToken<HcalDbService, HcalDbRecord> serviceToken_;
  bool LUTGenerationMode_, linearLUTs_;
  bool contain1TSHB_;
  double containPhaseNSHB_;
  bool applyFixPCC_;
  bool overrideDBweightsAndFilterHB_;
  double nPedWidthsForZS_;
  bool overrideDBnPedWidthsForZS_;
  double linearLSB_QIE11_;
  int maskBit_;
  bool overrideFGHF_;
  std::array<uint32_t, 2> FG_HF_thresholds_;
  bool overrideHBLLP_;
  std::array<uint32_t, 4> HB_LLP_thresholds_;
};

HcalTPGCoderULUT::HcalTPGCoderULUT(const edm::ParameterSet& iConfig) {
  contain1TSHB_ = iConfig.getParameter<bool>("contain1TSHB");
  containPhaseNSHB_ = iConfig.getParameter<double>("containPhaseNSHB");
  overrideDBweightsAndFilterHB_ = iConfig.getParameter<bool>("overrideDBweightsAndFilterHB");
  nPedWidthsForZS_ = iConfig.getParameter<double>("nPedWidthsForZS");
  overrideDBnPedWidthsForZS_ = iConfig.getParameter<bool>("overrideDBnPedWidthsForZS");
  applyFixPCC_ = iConfig.getParameter<bool>("applyFixPCC");

  auto cc = setWhatProduced(this);
  topoToken_ = cc.consumes();
  delayToken_ = cc.consumes(edm::ESInputTag{"", "HBHE"});
  serviceToken_ = cc.consumes();

  LUTGenerationMode_ = iConfig.getParameter<bool>("LUTGenerationMode");
  auto scales = iConfig.getParameter<edm::ParameterSet>("tpScales").getParameter<edm::ParameterSet>("HBHE");
  linearLSB_QIE11_ = scales.getParameter<double>("LSBQIE11");
  maskBit_ = iConfig.getParameter<int>("MaskBit");
  overrideFGHF_ = iConfig.getParameter<bool>("overrideFGHF");
  FG_HF_thresholds_ = iConfig.getParameter<std::array<uint32_t, 2> >("FG_HF_thresholds");
  overrideHBLLP_ = iConfig.getParameter<bool>("overrideHBLLP");
  HB_LLP_thresholds_ = iConfig.getParameter<std::array<uint32_t, 4> >("HB_LLP_thresholds");
}

void HcalTPGCoderULUT::buildCoder(const HcalTopology* topo,
                                  const HcalElectronicsMap* emap,
                                  const HcalTimeSlew* delay,
                                  HcaluLUTTPGCoder* theCoder) {
  using namespace edm::es;

  theCoder->init(topo, emap, delay);

  theCoder->setOverrideDBweightsAndFilterHB(overrideDBweightsAndFilterHB_);
  theCoder->set1TSContainHB(contain1TSHB_);
  theCoder->setContainPhaseHB(containPhaseNSHB_);
  theCoder->setNpedWidthsForZS(nPedWidthsForZS_);
  theCoder->setApplyFixPCC(applyFixPCC_);
  theCoder->setAllLinear(linearLSB_QIE11_);
  theCoder->setLUTGenerationMode(LUTGenerationMode_);
  theCoder->setMaskBit(maskBit_);
  theCoder->setOverrideFGHF(overrideFGHF_);
  theCoder->setFGHFthresholds(FG_HF_thresholds_);
  theCoder->setOverrideHBLLP(overrideHBLLP_);
  theCoder->setHBLLPthresholds(HB_LLP_thresholds_);
}

HcalTPGCoderULUT::~HcalTPGCoderULUT() {}

HcalTPGCoderULUT::ReturnType HcalTPGCoderULUT::produce(const HcalTPGRecord& iRecord) {
  auto host = holder_.makeOrGet([]() { return new HostType; });

  const auto& topo = iRecord.get(topoToken_);
  const auto& delayRcd = iRecord.getRecord<HcalDbRecord>();
  const auto& dbServ = iRecord.get(serviceToken_);
  const auto* emap = dbServ.getHcalMapping();
  const auto& delay = delayRcd.get(delayToken_);

  host->ifRecordChanges<HcalDbRecord>(iRecord, [this, &topo, emap, &delay, h = host.get()](auto const& rec) {
    buildCoder(&topo, emap, &delay, h);
    h->update(rec.get(serviceToken_));
  });

  return host;
}

//define this as a plug-in
DEFINE_FWK_EVENTSETUP_MODULE(HcalTPGCoderULUT);

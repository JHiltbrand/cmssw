import FWCore.ParameterSet.Config as cms

from CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi import CaloTPGTranscoder
from CalibCalorimetry.CaloTPG.tpScales_cff import tpScales
from CalibCalorimetry.HcalPlugins.Hcal_PCCUpdate_cff import PCCUpdate
from SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cfi import simHcalTriggerPrimitiveDigis

HcalTPGCoderULUT = cms.ESProducer("HcalTPGCoderULUT",
    LUTGenerationMode = cms.bool(True),
    contain1TSHB = cms.bool(True),
    containPhaseNSHB = cms.double(6.0),
    applyFixPCC = PCCUpdate.applyFixPCC,
    overrideDBweightsAndFilterHB = cms.bool(False),
    nPedWidthsForZS = cms.double(0.0),
    overrideDBnPedWidthsForZS = cms.bool(False),
    tpScales = tpScales,
    MaskBit = cms.int32(0x8000),
    overrideFGHF = cms.bool(False),
    FG_HF_thresholds = cms.vuint32(17, 255),
    overrideHBLLP = cms.bool(False),
    HB_LLP_thresholds = cms.vuint32(16, 80, 64, 64)
)

HcalTrigTowerGeometryESProducer = cms.ESProducer("HcalTrigTowerGeometryESProducer")

#Starting from CMSSW_15_1_X, the FG_HF_thresholds are now controlled by the TPParameter tag in the global tag, and further modification of these settings is no longer needed


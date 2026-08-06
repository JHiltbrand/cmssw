import FWCore.ParameterSet.Config as cms

from CalibCalorimetry.CaloTPG.tpScales_cff import tpScales

LSParameter =cms.untracked.PSet(
    HcalFeatureHFEMBit= cms.bool(False),
    Min_Long_Energy= cms.double(10),#makes a cut based on energy deposited in short vrs long
    Min_Short_Energy= cms.double(10),
    Long_vrs_Short_Slope= cms.double(100.2),
    Long_Short_Offset= cms.double(10.1)
)

simHcalTriggerPrimitiveDigis = cms.EDProducer("HcalTrigPrimDigiProducer",
    weightsQIE11 = cms.PSet(
        ieta1 = cms.vint32(255, 255),
        ieta2 = cms.vint32(255, 255),
        ieta3 = cms.vint32(255, 255),
        ieta4 = cms.vint32(255, 255),
        ieta5 = cms.vint32(255, 255),
        ieta6 = cms.vint32(255, 255),
        ieta7 = cms.vint32(255, 255),
        ieta8 = cms.vint32(255, 255),
        ieta9 = cms.vint32(255, 255),
        ieta10 = cms.vint32(255, 255),
        ieta11 = cms.vint32(255, 255),
        ieta12 = cms.vint32(255, 255),
        ieta13 = cms.vint32(255, 255),
        ieta14 = cms.vint32(255, 255),
        ieta15 = cms.vint32(255, 255),
        ieta16 = cms.vint32(255, 255)
    ),

    FG_HF_thresholds = cms.vuint32(17, 255), ## thresholds for setting fine grain bit

    # To be used when overriding the CondDB, default is with vetoing off ("coded" threshold = 0)
    # To run PFA1' + vetoing with no threshold, use 2048
    # All other values (1, 2047) are interpreted literally as the PFA1' veto threshold 
    codedVetoThresholds = cms.PSet(
        ieta1  = cms.int32(0),
        ieta2  = cms.int32(0),
        ieta3  = cms.int32(0),
        ieta4  = cms.int32(0),
        ieta5  = cms.int32(0),
        ieta6  = cms.int32(0),
        ieta7  = cms.int32(0),
        ieta8  = cms.int32(0),
        ieta9  = cms.int32(0),
        ieta10 = cms.int32(0),
        ieta11 = cms.int32(0),
        ieta12 = cms.int32(0),
        ieta13 = cms.int32(0),
        ieta14 = cms.int32(0),
        ieta15 = cms.int32(0),
        ieta16 = cms.int32(0)
    ),

    overrideHBLLP = cms.bool(False), ## switch: False = read thresholds from TPParameters (default), True = override with HB_LLP_thresholds                                                         
    ## defaults for energy requirement for bits 12-15 are high / low to avoid FG bit 0-4 being set when not intended                                                                                              
    HB_LLP_thresholds = cms.vuint32(0, 0, 999, 999),  ## default energy thresholds for setting HB LLP bit                                                                                           
                                                      ## depths 1,2 max energy, depths 3+ min energy, prompt min energy, delayed min energy                                            

    useTDCfromDB = cms.bool(False), ## switch: False = read HB TDC LUT thresholds from hard-coded values, True = read them from TPChannelParameters

    overrideDBvetoThresholdsHB = cms.bool(False),
    numberOfSamples = cms.int32(4),
    numberOfPresamples = cms.int32(2),
    numberOfSamplesHF = cms.int32(2),
    numberOfPresamplesHF = cms.int32(1),
    numberOfFilterPresamplesHBQIE11 = cms.int32(0),
    useTDCInMinBiasBits = cms.bool(False), # TDC information not used in MB fine grain bits
    LSConfig=LSParameter,

    applySaturationFix = cms.bool(True), # Apply the TP energy saturation fix for Peak Finder Algorithm only for >= Run-3

    # Input digi label (_must_ be without zero-suppression!)
    inputUpgradeLabel = cms.VInputTag(
        cms.InputTag('simHcalUnsuppressedDigis:HBHEQIE11DigiCollection'),
        cms.InputTag('simHcalUnsuppressedDigis:HFQIE10DigiCollection')),
    overrideDBweightsAndFilterHB = cms.bool(False),

    tpScales = tpScales,
)


import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("RESP", eras.Phase2C9)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
            wantSummary = cms.untracked.bool(True),
            numberOfThreads = cms.untracked.uint32(2),
            numberOfStreams = cms.untracked.uint32(1),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/Phase2HLTTDRWinter20DIGI/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW/PU200_110X_mcRun4_realistic_v3-v2/110000/005E74D6-B50E-674E-89E6-EAA9A617B476.root',)
)

process.load('Configuration.Geometry.GeometryExtendedRun4D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D49_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '123X_mcRun4_realistic_v3', '')

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load("L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff") 
process.load("L1Trigger.TrackerDTC.DTC_cff") 
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

process.l1tLayer1Barrel9 = process.l1tLayer1Barrel.clone()
process.l1tLayer1Barrel9.puAlgo.nFinalSort = 32
process.l1tLayer1Barrel9.regions[0].etaBoundaries = [ -1.5, -0.5, 0.5, 1.5 ] 
process.l1tLayer1Barrel9.boards=cms.VPSet(
        cms.PSet(
            regions=cms.vuint32(*[0+9*ie+i for ie in range(3) for i in range(3)])),
        cms.PSet(
            regions=cms.vuint32(*[3+9*ie+i for ie in range(3) for i in range(3)])),
        cms.PSet(
            regions=cms.vuint32(*[6+9*ie+i for ie in range(3) for i in range(3)])),
    )

process.PFInputsTask = cms.Task(
    process.TTClustersFromPhase2TrackerDigis,
    process.TTStubsFromPhase2TrackerDigis,
    process.ProducerDTC,
    process.offlineBeamSpot,
    process.l1tTTTracksFromTrackletEmulation,
    process.SimL1EmulatorTask
)
process.RunPF = cms.Path( 
        process.l1tLayer1Barrel +
        process.l1tLayer1Barrel9 +
        process.l1tLayer1HGCal +
        process.l1tLayer1HGCalNoTK +
        process.l1tLayer1HF
)

process.RunPF.associate(process.PFInputsTask)
process.schedule = cms.Schedule(process.RunPF)

for det in "Barrel", "Barrel9", "HGCal", "HGCalNoTK", "HF":
    l1pf = getattr(process, 'l1tLayer1'+det)
    l1pf.dumpFileName = cms.untracked.string("TTbar_PU200_123X_"+det+".dump")

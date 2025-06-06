import FWCore.ParameterSet.Config as cms

#Primary vertex selection
PVCut = "abs(z) < 24. & position.rho < 2. & ndof > 4 & !isFake"

#Jet selection
looseJetCut = "(chargedHadronEnergyFraction()>0 && chargedMultiplicity()>0 && chargedEmEnergyFraction()<0.99 && neutralHadronEnergyFraction()<0.99 && neutralEmEnergyFraction()<0.99 && (chargedMultiplicity()+neutralMultiplicity())>1) && abs(eta)<=2.4 "

tightJetCut = "(chargedHadronEnergyFraction()>0 && chargedMultiplicity()>0 && chargedEmEnergyFraction()<0.99 && neutralHadronEnergyFraction()<0.90 && neutralEmEnergyFraction()<0.90 && (chargedMultiplicity()+neutralMultiplicity())>1) && abs(eta)<=2.4 "

#Loose muon selection
looseMuonCut  = "((isGlobalMuon || isTrackerMuon) && isPFMuon)"
looseIsoCut   = "((pfIsolationR04.sumChargedHadronPt + max(0., pfIsolationR04.sumNeutralHadronEt + pfIsolationR04.sumPhotonEt - 0.5 * pfIsolationR04.sumPUPt) ) / pt < 0.25)"

#Medium muon selection. Also requires either good global muon or tight segment compatibility
mediumMuonCut = looseMuonCut + " innerTrack.validFraction > 0.8"

#Tight muon selection. Lacks distance to primary vertex variables, dz<0.5, dxy < 0.2. Now done at .cc
tightMuonCut  = "isGlobalMuon && isPFMuon && globalTrack.normalizedChi2 < 10. && globalTrack.hitPattern.numberOfValidMuonHits > 0 && " + \
               "numberOfMatchedStations > 1 && innerTrack.hitPattern.numberOfValidPixelHits > 0 && innerTrack.hitPattern.trackerLayersWithMeasurement > 5 "
tightIsoCut   = "(pfIsolationR04.sumChargedHadronPt + max(0., pfIsolationR04.sumNeutralHadronEt + pfIsolationR04.sumPhotonEt - 0.5 * pfIsolationR04.sumPUPt) ) / pt < 0.15"

#Electron selections
looseEleCut = "((full5x5_sigmaIetaIeta() < 0.011 && superCluster().isNonnull() && superCluster().seed().isNonnull() && (deltaEtaSuperClusterTrackAtVtx() - superCluster().eta() + superCluster().seed().eta()) < 0.00477 && abs(deltaPhiSuperClusterTrackAtVtx()) < 0.222 && hadronicOverEm() < 0.298 && abs(1.0 - eSuperClusterOverP())*1.0/ecalEnergy() < 0.241 && gsfTrack.hitPattern().numberOfLostHits('MISSING_INNER_HITS') <= 1 && abs(superCluster().eta()) < 1.479) || (full5x5_sigmaIetaIeta() < 0.0314 && superCluster().isNonnull() && superCluster().seed().isNonnull() && (deltaEtaSuperClusterTrackAtVtx() - superCluster().eta() + superCluster().seed().eta()) < 0.00868 && abs(deltaPhiSuperClusterTrackAtVtx()) < 0.213 && hadronicOverEm() < 0.101  && abs(1.0 - eSuperClusterOverP())*1.0/ecalEnergy() < 0.14 && gsfTrack.hitPattern().numberOfLostHits('MISSING_INNER_HITS') <= 1 && abs(superCluster().eta()) > 1.479))"

tightEleCut = "((full5x5_sigmaIetaIeta() < 0.00998 && superCluster().isNonnull() && superCluster().seed().isNonnull() && (deltaEtaSuperClusterTrackAtVtx() - superCluster().eta() + superCluster().seed().eta()) < 0.00308  && abs(deltaPhiSuperClusterTrackAtVtx()) < 0.0816 && hadronicOverEm() < 0.0414 && abs(1.0 - eSuperClusterOverP())*1.0/ecalEnergy() < 0.0129 && gsfTrack.hitPattern().numberOfLostHits('MISSING_INNER_HITS') <= 1 && abs(superCluster().eta()) < 1.479) ||  (full5x5_sigmaIetaIeta() < 0.0292 && superCluster().isNonnull() && superCluster().seed().isNonnull() && (deltaEtaSuperClusterTrackAtVtx() - superCluster().eta() + superCluster().seed().eta()) < 0.00605 && abs(deltaPhiSuperClusterTrackAtVtx()) < 0.0394  && hadronicOverEm() < 0.0641  && abs(1.0 - eSuperClusterOverP())*1.0/ecalEnergy() <	0.0129 && gsfTrack.hitPattern().numberOfLostHits('MISSING_INNER_HITS') <= 1 && abs(superCluster().eta()) > 1.479))"


from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer

topSingleMuonMediumDQM = DQMEDAnalyzer('TopSingleLeptonDQM',
  ## ------------------------------------------------------
  ## SETUP
  ##
  ## configuration of the MonitoringEnsemble(s)
  ## [mandatory] : optional PSets may be omitted
  ##
  setup = cms.PSet(
    directory = cms.string("Physics/Top/TopSingleMuonMediumDQM/"),
    sources = cms.PSet(
      muons = cms.InputTag("muons"),
      elecs = cms.InputTag("gedGsfElectrons"),
      jets  = cms.InputTag("ak4PFJetsCHS"),
      mets  = cms.VInputTag("pfMet"),
      pvs   = cms.InputTag("offlinePrimaryVertices")

    ),
    monitoring = cms.PSet(
      verbosity = cms.string("DEBUG")
    ),
    pvExtras = cms.PSet(
      select = cms.string(PVCut)
    ),
    elecExtras = cms.PSet(
      select     = cms.string(tightEleCut + "& pt>20 & abs(eta)<2.5 & (abs(superCluster().eta()) <= 1.4442 || abs(superCluster().eta()) >= 1.5660)"),
			rho = cms.InputTag("fixedGridRhoFastjetAll"),
    ),                                     
    muonExtras = cms.PSet(
      select    = cms.string(tightMuonCut + " && pt>20 & abs(eta)<2.4"),                                               
      #select    = cms.string(looseMuonCut + " && pt>20 & abs(eta)<2.4"),                                              
      isolation = cms.string(looseIsoCut)
    ),
    jetExtras = cms.PSet(
      jetCorrector = cms.InputTag("dqmAk4PFCHSL1FastL2L3Corrector"), #Use pak4PFCHSL1FastL2L3Residual for data!!!                                            
      select = cms.string("pt>30 & abs(eta)<2.4"),                                                                                               
      jetBTaggers = cms.PSet(
				cvsVertex = cms.PSet(
            label = cms.InputTag("pfDeepCSVJetTags:probb"),
	          workingPoint = cms.double(0.4168)
	          # deepCSV Medium from https://btv-wiki.docs.cern.ch/ScaleFactors/UL2018/
        )
      ),                                                                                               
    ),
    massExtras = cms.PSet(
      lowerEdge = cms.double( 70.),
      upperEdge = cms.double(110.)
    ),
  ),

  preselection = cms.PSet(
    vertex = cms.PSet(
      src    = cms.InputTag("offlinePrimaryVertices"),
      select = cms.string(PVCut)
    )
  ),

  selection = cms.VPSet(
    cms.PSet(
      label  = cms.string("muons:step0"),
      src    = cms.InputTag("muons"),
      select = cms.string(tightMuonCut + " && pt>20 & abs(eta)<2.4"),     
      min    = cms.int32(1),
    ),
    cms.PSet(
      label  = cms.string("jets/pf:step1"),
      src    = cms.InputTag("ak4PFJetsCHS"),
      select = cms.string("pt>30 & abs(eta)<2.4"),
      min = cms.int32(4),
    ), 
    cms.PSet(
      label  = cms.string("met:step2"),
      src    = cms.InputTag("pfMet"),
      select = cms.string("pt>30"),                                             
    ),
  )
)


topSingleElectronMediumDQM = DQMEDAnalyzer('TopSingleLeptonDQM',
  ## ------------------------------------------------------
  ## SETUP
  ##
  ## configuration of the MonitoringEnsemble(s)
  ## [mandatory] : optional PSets may be omitted
  ##
  setup = cms.PSet(
    directory = cms.string("Physics/Top/TopSingleElectronMediumDQM/"),
    sources = cms.PSet(
      muons = cms.InputTag("muons"),
      elecs = cms.InputTag("gedGsfElectrons"),
      jets  = cms.InputTag("ak4PFJetsCHS"),
      mets  = cms.VInputTag("pfMet"),
      pvs   = cms.InputTag("offlinePrimaryVertices")
    ),
    monitoring = cms.PSet(
      verbosity = cms.string("DEBUG")
    ),
    pvExtras = cms.PSet(                                                                                           
      select = cms.string(PVCut)
    ),
    elecExtras = cms.PSet(
      select   = cms.string(tightEleCut + "& pt>20 & abs(eta)<2.5 & (abs(superCluster().eta()) <= 1.4442 || abs(superCluster().eta()) >= 1.5660)"),
			rho = cms.InputTag("fixedGridRhoFastjetAll"),
    ),
    muonExtras = cms.PSet(
      select     = cms.string(tightMuonCut + " & pt>20 & abs(eta)<2.4"),
      isolation  = cms.string(looseIsoCut),
    ),
    jetExtras = cms.PSet(
			jetCorrector = cms.InputTag("dqmAk4PFCHSL1FastL2L3Corrector"), #Use pak4PFCHSL1FastL2L3Residual for data!!!
      select = cms.string("pt>30 & abs(eta)<2.4"),
      jetBTaggers = cms.PSet(
			  cvsVertex = cms.PSet(
          label = cms.InputTag("pfDeepCSVJetTags:probb"),
	        workingPoint = cms.double(0.4168)
	          # deepCSV Medium from https://btv-wiki.docs.cern.ch/ScaleFactors/UL2018/
        )
      )                                                
    ),
    massExtras = cms.PSet(
      lowerEdge = cms.double( 70.),
      upperEdge = cms.double(110.)
    ),
  ),
  preselection = cms.PSet(
    vertex = cms.PSet(
      src    = cms.InputTag("offlinePrimaryVertices"),
      select = cms.string(PVCut)
    )
  ),
  selection = cms.VPSet(
    cms.PSet(
      label = cms.string("elecs:step0"),
      src   = cms.InputTag("gedGsfElectrons"),
      select = cms.string("pt>20 & abs(eta)<2.5 & (abs(superCluster().eta()) <= 1.4442 || abs(superCluster().eta()) >= 1.5660) &&" + tightEleCut),
 #    select = cms.string("pt>30 & abs(eta)<2.5 & abs(gsfTrack.d0)<0.02 & gsfTrack.hitPattern().numberOfLostHits('MISSING_INNER_HITS') <= 0 & (abs(superCluster.eta) <= 1.4442 || abs(superCluster.eta) >= 1.5660) & " + EletightIsoCut),
      min = cms.int32(1),
    ),
    cms.PSet(
      label = cms.string("jets/pf:step1"),
      src   = cms.InputTag("ak4PFJetsCHS"),
      select = cms.string("pt>30 & abs(eta)<2.4"),
      min = cms.int32(4),
    ), 
    cms.PSet(
      label  = cms.string("met:step2"),
      src    = cms.InputTag("pfMet"),
      select = cms.string("pt>30"),
    ),
  )
)


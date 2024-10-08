/*! \class   TTClusterAssociator
 *  \brief   Plugin to create the MC truth for TTClusters.
 *  \details After moving from SimDataFormats to DataFormats,
 *           the template structure of the class was maintained
 *           in order to accomodate any types other than PixelDigis
 *           in case there is such a need in the future.
 *
 *  \author Nicola Pozzobon
 *  \date   2013, Jul 19
 *  (tidy up: Ian Tomalin, 2020)
 *
 */

#ifndef L1_TRACK_TRIGGER_CLUSTER_ASSOCIATOR_H
#define L1_TRACK_TRIGGER_CLUSTER_ASSOCIATOR_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "L1Trigger/TrackTrigger/interface/classNameFinder.h"
#include "SimDataFormats/Associations/interface/TTClusterAssociationMap.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"

#include "L1Trigger/TrackTrigger/interface/TTStubAlgorithm.h"
#include "L1Trigger/TrackTrigger/interface/TTStubAlgorithmRecord.h"

#include "Geometry/CommonTopologies/interface/Topology.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include <memory>
#include <map>
#include <vector>

template <typename T>
class TTClusterAssociator : public edm::stream::EDProducer<> {
  /// NOTE since pattern hit correlation must be performed within a stacked module, one must store
  /// Clusters in a proper way, providing easy access to them in a detector/member-wise way
public:
  /// Constructors
  explicit TTClusterAssociator(const edm::ParameterSet& iConfig);

private:
  /// Data members
  edm::Handle<edm::DetSetVector<PixelDigiSimLink> > thePixelDigiSimLinkHandle_;
  edm::Handle<std::vector<TrackingParticle> > trackingParticleHandle_;

  std::vector<edm::InputTag> ttClustersInputTags_;

  edm::EDGetTokenT<edm::DetSetVector<PixelDigiSimLink> > digisimLinkToken_;
  edm::EDGetTokenT<std::vector<TrackingParticle> > tpToken_;
  std::vector<edm::EDGetTokenT<edmNew::DetSetVector<TTCluster<T> > > > ttClustersTokens_;

  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> theTrackerGeometryToken_;

  /// Mandatory methods
  void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

};  /// Close class

/*! \brief   Implementation of methods
 *  \details Here, in the header file, the methods which do not depend
 *           on the specific type <T> that can fit the template.
 *           Other methods, with type-specific features, are implemented
 *           in the source file.
 */

/// Constructors
template <typename T>
TTClusterAssociator<T>::TTClusterAssociator(const edm::ParameterSet& iConfig) {
  digisimLinkToken_ =
      consumes<edm::DetSetVector<PixelDigiSimLink> >(iConfig.getParameter<edm::InputTag>("digiSimLinks"));
  tpToken_ = consumes<std::vector<TrackingParticle> >(iConfig.getParameter<edm::InputTag>("trackingParts"));

  ttClustersInputTags_ = iConfig.getParameter<std::vector<edm::InputTag> >("TTClusters");

  for (const auto& iTag : ttClustersInputTags_) {
    ttClustersTokens_.push_back(consumes<edmNew::DetSetVector<TTCluster<T> > >(iTag));

    produces<TTClusterAssociationMap<T> >(iTag.instance());
  }

  theTrackerGeometryToken_ = esConsumes();

  /// Print some information when loaded
  edm::LogInfo("TTClusterAssociator< ") << templateNameFinder<T>() << " > loaded.";
}

/// Implement the producer
template <>
void TTClusterAssociator<Ref_Phase2TrackerDigi_>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup);

#endif

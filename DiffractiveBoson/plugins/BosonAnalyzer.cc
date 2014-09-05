// -*- C++ -*-
//
// Package:    BosonAnalyzer
// Class:      BosonAnalyzer
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/LorentzVector.h"

// dataformats
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

// Math Tools
#include "DataFormats/Math/interface/Point3D.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

// root
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TTree.h"

// c++
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <map>

#define NELEMS(x)  (sizeof(x) / sizeof(x[0]));

using namespace edm;
using namespace std;
using namespace HepMC;

//
// class declaration
//

class BosonAnalyzer : public edm::EDAnalyzer {
  public:
    explicit BosonAnalyzer(const edm::ParameterSet&);
    ~BosonAnalyzer();


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    void Init();

    // ----------member data ---------------------------

    void FillCollections(const edm::Event&, const edm::EventSetup&, bool debug);
    void SavingInformation();
    void PrintOrder();

    template <class T, class W>
      math::XYZTLorentzVector DiSystem(T obj1, W obj2);

    template <class T, class W>
      double MassT(T lepton, W met);

    template <class T, class W>
      double InvariantMass(T lepton1, W lepton2);

    edm::InputTag particleFlowTag_;
    edm::InputTag VertexTag_;
    double pTPFThresholdCharged_;
    double energyPFThresholdBar_;
    double energyPFThresholdEnd_;
    double energyPFThresholdHF_;
    edm::InputTag GenParticleTag_;
    double EBeam_;
    bool debug_;
    edm::InputTag electronTag_;
    edm::InputTag muonTag_;
    edm::InputTag metTag_;

    std::vector<const reco::Vertex*> VertexVector;
    std::vector< math::XYZVector > VertexPosition;
    std::vector<const reco::Track*> TracksVector;
    std::vector<const reco::PFCandidate*> PFVector;
    std::vector<const reco::GenParticle*> protonVector;
    std::vector<const reco::Muon*> MuonVector;
    std::vector<const reco::GsfElectron*> ElectronVector;
    std::vector<const reco::PFMET*> NeutrinoVector;

    std::vector< math::XYZTLorentzVector > protonLorentzVector;
    std::vector< math::XYZTLorentzVector > LeadingElectronsP4;
    std::vector< math::XYZTLorentzVector > LeadingMuonsP4;
    std::vector< math::XYZTLorentzVector >  METP4;

    /*
    // Possible to use both cases
    std::vector <ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > protonLorentzVector;
    std::vector <ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > LeadingElectronsP4;
    std::vector <ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > LeadingMuonsP4;
    std::vector <ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > METP4;
     */

    double METSumEt;

    int nTracks;
    int nVertex;
    double ZMassDiElectron;
    double WMassElectron;
    double ZMassDiMuon;
    double WMassMuon;

    double ZEtaDiElectron;
    double ZEtaDiMuon;
    double WEtaElectron;
    double WEtaMuon;

    std::vector<double> LeadingElectronsIsolation;
    std::vector<double> LeadingElectronsIsoEcal;
    std::vector<double> LeadingElectronsIsoHcal;
    std::vector<double> LeadingElectronsInnerHits;
    std::vector<double> LeadingElectronsDCot;
    std::vector<double> LeadingElectronsDist;
    std::vector<double> LeadingElectronsDeltaEtaTkClu;
    std::vector<double> LeadingElectronsDeltaPhiTkClu;
    std::vector<double> LeadingElectronsSigmaIeIe;
    std::vector<double> LeadingElectronsHE;
    std::vector<double> LeadingMuonsIsolation;

    bool ZMassE;
    bool ZMassMu;
    bool WMassE;
    bool WMassMu;

    bool ZFillE;
    bool ZFillMu;
    bool WFillE;
    bool WFillMu;

    TTree* eventTree_;

    struct orderPT
    {
      template <class T, class W>
	inline bool operator() (T vec1, W vec2)
	{
	  return (vec1->pt() > vec2->pt());
	}
    };

    struct orderETA
    {
      template <class T, class W>
	inline bool operator() (T vec1, W vec2)
	{
	  return (vec1->eta() > vec2->eta());
	}
    };

    struct orderAbsolutPZ
    {
      template <class T, class W>
	inline bool operator() (T vec1, W vec2)
	{
	  return (fabs(vec1->pz()) > fabs(vec2->pz()));
	}
    };

};


//
// constructors and destructor
//

BosonAnalyzer::BosonAnalyzer(const edm::ParameterSet& iConfig):
  particleFlowTag_(iConfig.getParameter<edm::InputTag>("ParticleFlowTag")),
  VertexTag_(iConfig.getParameter<edm::InputTag>("VertexTag")),
  pTPFThresholdCharged_(iConfig.getParameter<double>("pTPFThresholdCharged")),
  energyPFThresholdBar_(iConfig.getParameter<double>("energyPFThresholdBar")),
  energyPFThresholdEnd_(iConfig.getParameter<double>("energyPFThresholdEnd")),
  energyPFThresholdHF_(iConfig.getParameter<double>("energyPFThresholdHF")),
  GenParticleTag_(iConfig.getParameter<edm::InputTag>("GenParticleTag")),
  EBeam_(iConfig.getParameter<double>("EBeam")),
  debug_(iConfig.getParameter<bool>("debug")),
  electronTag_(iConfig.getParameter<edm::InputTag>("electronTag")),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonTag")),
  metTag_(iConfig.getParameter<edm::InputTag>("metTag"))
{

  edm::Service<TFileService> fs;
  eventTree_ = fs->make<TTree>("Event","Event");
  eventTree_->Branch("VertexPosition",&VertexPosition);
  eventTree_->Branch("ProtonsP4",&protonLorentzVector);
  eventTree_->Branch("ElectronsP4",&LeadingElectronsP4);
  eventTree_->Branch("MuonsP4",&LeadingMuonsP4);
  eventTree_->Branch("METP4",&METP4);
  eventTree_->Branch("METSumEt",&METSumEt,"METSumEt/D");
  eventTree_->Branch("LeadingElectronsIsolation",&LeadingElectronsIsolation);
  eventTree_->Branch("LeadingElectronsIsoEcal",&LeadingElectronsIsoEcal);
  eventTree_->Branch("LeadingElectronsIsoHcal",&LeadingElectronsIsoHcal);
  eventTree_->Branch("LeadingElectronsInnerHits",&LeadingElectronsInnerHits);
  eventTree_->Branch("LeadingElectronsDCot",&LeadingElectronsDCot);
  eventTree_->Branch("LeadingElectronsDist",&LeadingElectronsDist);
  eventTree_->Branch("LeadingElectronsDeltaEtaTkClu",&LeadingElectronsDeltaEtaTkClu);
  eventTree_->Branch("LeadingElectronsDeltaPhiTkClu",&LeadingElectronsDeltaPhiTkClu);
  eventTree_->Branch("LeadingElectronsSigmaIeIe",&LeadingElectronsSigmaIeIe);
  eventTree_->Branch("LeadingElectronsHE",&LeadingElectronsHE);
  eventTree_->Branch("LeadingMuonsIsolation",&LeadingMuonsIsolation);
  eventTree_->Branch("ZMassDiElectron",&ZMassDiElectron,"ZMassDiElectron/D");
  eventTree_->Branch("WMassElectron",&WMassElectron,"WMassElectron/D");
  eventTree_->Branch("ZMassDiMuon",&ZMassDiMuon,"ZMassDiMuon/D");
  eventTree_->Branch("WMassMuon",&WMassMuon,"WMassMuon/D");
  eventTree_->Branch("ZEtaDiElectron",&ZEtaDiElectron,"ZEtaDiElectron/D");
  eventTree_->Branch("ZEtaDiMuon",&ZEtaDiMuon,"ZEtaDiMuon/D");
  eventTree_->Branch("WEtaElectron",&WEtaElectron,"WEtaElectron/D");
  eventTree_->Branch("WEtaMuon",&WEtaMuon,"WEtaMuon/D");
  eventTree_->Branch("nVertex",&nVertex,"nVertex/I");
  eventTree_->Branch("nTracks",&nTracks,"nTracks/I");

}

BosonAnalyzer::~BosonAnalyzer()
{
}


// ------------ method called once each job just before starting event loop  ------------
void BosonAnalyzer::beginJob()
{
  cout << "\n--- B O S O N   A N A L Y Z E R---\n" << endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void BosonAnalyzer::endJob()
{

  cout << "\n--- S U M M A R Y---" << endl;
  cout << "\n--- E N D---\n" << endl;

}

// ------------ method called for each event  ------------
void BosonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  Init(); // Clean Variables
  FillCollections(iEvent, iSetup, false); // Debug option, true = show Boson Mass, false = do not show Boson Mass.
  if(debug_) PrintOrder();
  SavingInformation();
  eventTree_->Fill();

}

// ------------ Clean Variables --------------

void BosonAnalyzer::Init(){

  VertexVector.clear();
  VertexPosition.clear();
  TracksVector.clear();
  PFVector.clear();
  protonVector.clear();
  MuonVector.clear();
  ElectronVector.clear();
  NeutrinoVector.clear();
  protonLorentzVector.clear();
  LeadingElectronsP4.clear();
  LeadingMuonsP4.clear();
  LeadingElectronsIsolation.clear();
  LeadingElectronsIsoEcal.clear();
  LeadingElectronsIsoHcal.clear();
  LeadingElectronsInnerHits.clear();
  LeadingElectronsDCot.clear();
  LeadingElectronsDist.clear();
  LeadingElectronsDeltaEtaTkClu.clear();
  LeadingElectronsDeltaPhiTkClu.clear();
  LeadingElectronsSigmaIeIe.clear();
  LeadingElectronsHE.clear();
  LeadingMuonsIsolation.clear();
  METP4.clear();

  nTracks = -999;
  nVertex = -999;
  ZMassDiElectron = -999.;
  WMassElectron = -999.;
  ZMassDiMuon = -999.;
  WMassMuon = -999.;
  ZEtaDiElectron = -999.;
  ZEtaDiMuon = -999.;
  WEtaElectron = -999.;
  WEtaMuon = -999.;
  METSumEt = -999.;

  ZMassE = false;
  ZMassMu = false;
  WMassE = false;
  WMassMu = false;

  ZFillE = false;
  ZFillMu = false;
  WFillE = false;
  WFillMu = false;

}


// ------------ Fill Vectors, All Handles  ------------
void BosonAnalyzer::FillCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool debug)
{

  // Fill Vertex
  Handle<edm::View<reco::Vertex> > vertex;
  iEvent.getByLabel(VertexTag_, vertex);

  int vertexsize = vertex->size();
  int itVertex;

  if(vertex->size()>0){
    for(itVertex=0; itVertex < vertexsize; ++itVertex){
      const reco::Vertex* vertexAll = &((*vertex)[itVertex]);
      if (vertexAll->isValid()==0) continue; 
      VertexVector.push_back(vertexAll);
      math::XYZVector coord(vertexAll->x(),vertexAll->y(),vertexAll->z());
      VertexPosition.push_back(coord);
    }
  }

  // Fill Tracks
  Handle<edm::View<reco::Track> > tracks;
  iEvent.getByLabel("generalTracks", tracks);

  int trackssize = tracks->size();
  int itTracks;

  if(tracks->size()>0){
    for(itTracks=0; itTracks < trackssize; ++itTracks){
      const reco::Track* tracksAll = &((*tracks)[itTracks]);
      TracksVector.push_back(tracksAll);
    }
  }

  nTracks = TracksVector.size();
  nVertex = VertexVector.size();

  // Fill Particle Flow
  Handle <reco::PFCandidateCollection> PFCandidates;
  iEvent.getByLabel(particleFlowTag_,PFCandidates);
  reco::PFCandidateCollection::const_iterator iter;

  int pfsize = PFCandidates->size();
  int itPF;

  if(PFCandidates->size()>0){

    for(itPF=0; itPF < pfsize; ++itPF){
      const reco::PFCandidate* pfAll = &((*PFCandidates)[itPF]);
      double energy=pfAll->energy();
      double pt=pfAll->pt();
      double eta=pfAll->eta();
      double charge=pfAll->charge();
      if (fabs(eta)>4.7) continue;
      if ( (fabs(charge) >0 && pt > pTPFThresholdCharged_ ) ||
	  (fabs(charge) == 0 && ( (fabs(eta) <= 1.5 && energy > energyPFThresholdBar_) ||
				  (fabs(eta) > 1.5 && fabs(eta) <= 3 && energy > energyPFThresholdEnd_) ||
				  (fabs(eta) > 3 && energy >energyPFThresholdHF_) ) ) )
      { 
	PFVector.push_back(pfAll);
      }
    }

  }

  // Fill Protons
  Handle<reco::GenParticleCollection> genParticle;
  iEvent.getByLabel(GenParticleTag_, genParticle);
  int gensize = genParticle->size();
  int itGen;

  if(genParticle->size()>0){
    for(itGen=0; itGen < gensize; ++itGen){
      const reco::GenParticle* genAll = &((*genParticle)[itGen]);
      if (genAll->status() != 1) continue;
      if (fabs(genAll->pdgId()) != 2212) continue;
      protonVector.push_back(genAll);
    }
  }


  // Fill reco::Muon
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByLabel(muonTag_,muons);

  int muonsize = muons->size();
  int itMuon;
  MuonVector.clear();

  if(muons->size()>0){
    for(itMuon=0; itMuon < muonsize; ++itMuon){
      const reco::Muon* muonAll = &((*muons)[itMuon]);
      MuonVector.push_back(muonAll);
    }
  }

  // Fill reco::GsfElectron
  edm::Handle<reco::GsfElectronCollection> electrons;
  iEvent.getByLabel(electronTag_,electrons);

  int electronsize = electrons->size();
  int itElectron;
  ElectronVector.clear();

  if(electrons->size()>0){
    for(itElectron=0; itElectron < electronsize; ++itElectron){
      const reco::GsfElectron* electronAll = &((*electrons)[itElectron]);
      ElectronVector.push_back(electronAll);
    }
  }

  // Fill pfMET
  edm::Handle<reco::PFMETCollection> met;
  iEvent.getByLabel(metTag_,met);

  NeutrinoVector.clear();
  int neutrinosize = met->size();
  int itmet;

  if(met->size()>0){
    for(itmet=0; itmet < neutrinosize; ++itmet){
      const reco::PFMET* neutrinoAll = &((*met)[itmet]);
      NeutrinoVector.push_back(neutrinoAll);
    }

  }

  // Sorting Vectors, pT order
  std::sort(PFVector.begin(), PFVector.end(), orderPT());
  std::sort(protonVector.begin(), protonVector.end(), orderAbsolutPZ());
  std::sort(ElectronVector.begin(), ElectronVector.end(), orderPT());
  std::sort(MuonVector.begin(), MuonVector.end(), orderPT());
  std::sort(NeutrinoVector.begin(), NeutrinoVector.end(), orderPT());

  // Saving One or Two Leading Protons
  if (protonVector.size()==1){
    protonLorentzVector.push_back(protonVector[0]->p4());
  }

  if (protonVector.size()>1){
    protonLorentzVector.push_back(protonVector[0]->p4());
    protonLorentzVector.push_back(protonVector[1]->p4());
  }
  //--->

  if (ElectronVector.size()>1){
    ZMassDiElectron = InvariantMass(ElectronVector[0],ElectronVector[1]);
    ZEtaDiElectron = DiSystem(ElectronVector[0],ElectronVector[1]).eta();
  }

  if (MuonVector.size()>1){  
    ZMassDiMuon = InvariantMass(MuonVector[0],MuonVector[1]);
    ZEtaDiMuon = DiSystem(MuonVector[0],MuonVector[1]).eta();
  }

  if (ElectronVector.size()>0){
    WMassElectron = MassT(ElectronVector[0],NeutrinoVector[0]);
    WEtaElectron = DiSystem(ElectronVector[0],NeutrinoVector[0]).eta();
  }

  if (MuonVector.size()>0){
    WMassMuon = MassT(MuonVector[0],NeutrinoVector[0]);
    WEtaMuon = DiSystem(MuonVector[0],NeutrinoVector[0]).eta();
  }

  if(debug){
    std::cout << "M A S S" << std::endl;
    std::cout << "Z Mass (ee): " << ZMassDiElectron << " | eta: " << ZEtaDiElectron<<  std::endl;
    std::cout << "Z Mass (MuMu): " << ZMassDiMuon << " | eta: " << ZEtaDiMuon <<  std::endl;
    std::cout << "W Mass (eNu): " << WMassElectron << " | eta: " << WEtaElectron <<  std::endl;
    std::cout << "W Mass (MuNu): " << WMassMuon << " | eta: " << WEtaMuon <<  std::endl;
  }

  if(ZMassDiElectron > 60. && ZMassDiElectron < 110.) ZMassE = true;
  if(ZMassDiMuon > 60. && ZMassDiMuon < 110.) ZMassMu = true;
  if(WMassElectron > 60. && WMassElectron < 110.) WMassE = true;
  if(WMassMuon > 60. && WMassMuon < 110.) WMassMu = true;

  if(ZMassE && !ZMassMu && !WMassE && !WMassMu) ZFillE = true;
  if(ZMassMu && !ZMassE && !WMassE && !WMassMu) ZFillMu = true;
  if(WMassE && !ZMassE && !ZMassMu && !WMassMu) WFillE = true;
  if(WMassMu && !ZMassE && !ZMassMu && !WMassE) WFillMu = true;

}

void BosonAnalyzer::SavingInformation(){

  if(ZFillE){
    LeadingElectronsP4.push_back(ElectronVector[0]->p4());
    LeadingElectronsP4.push_back(ElectronVector[1]->p4());
    LeadingElectronsIsolation.push_back(ElectronVector[0]->dr03TkSumPt()/ElectronVector[0]->pt());
    LeadingElectronsIsolation.push_back(ElectronVector[1]->dr03TkSumPt()/ElectronVector[1]->pt());
    LeadingElectronsIsoEcal.push_back(ElectronVector[0]->dr03EcalRecHitSumEt()/ElectronVector[0]->pt());
    LeadingElectronsIsoEcal.push_back(ElectronVector[1]->dr03EcalRecHitSumEt()/ElectronVector[1]->pt());
    LeadingElectronsIsoHcal.push_back(ElectronVector[0]->dr03HcalTowerSumEt()/ElectronVector[0]->pt());
    LeadingElectronsIsoHcal.push_back(ElectronVector[1]->dr03HcalTowerSumEt()/ElectronVector[1]->pt());
    LeadingElectronsInnerHits.push_back(ElectronVector[0]->gsfTrack()->trackerExpectedHitsInner().numberOfHits());
    LeadingElectronsInnerHits.push_back(ElectronVector[1]->gsfTrack()->trackerExpectedHitsInner().numberOfHits());
    LeadingElectronsDCot.push_back(ElectronVector[0]->convDcot());
    LeadingElectronsDCot.push_back(ElectronVector[1]->convDcot());
    LeadingElectronsDist.push_back(ElectronVector[0]->convDist());
    LeadingElectronsDist.push_back(ElectronVector[1]->convDist());
    LeadingElectronsDeltaEtaTkClu.push_back(ElectronVector[0]->deltaEtaSuperClusterTrackAtVtx());
    LeadingElectronsDeltaEtaTkClu.push_back(ElectronVector[1]->deltaEtaSuperClusterTrackAtVtx());
    LeadingElectronsDeltaPhiTkClu.push_back(ElectronVector[0]->deltaPhiSuperClusterTrackAtVtx());
    LeadingElectronsDeltaPhiTkClu.push_back(ElectronVector[1]->deltaPhiSuperClusterTrackAtVtx());
    LeadingElectronsSigmaIeIe.push_back(ElectronVector[0]->sigmaIetaIeta());
    LeadingElectronsSigmaIeIe.push_back(ElectronVector[1]->sigmaIetaIeta());
    LeadingElectronsHE.push_back(ElectronVector[0]->hadronicOverEm());
    LeadingElectronsHE.push_back(ElectronVector[1]->hadronicOverEm());
  }

  if(ZFillMu){
    LeadingMuonsP4.push_back(MuonVector[0]->p4());
    LeadingMuonsP4.push_back(MuonVector[1]->p4());
    LeadingMuonsIsolation.push_back(MuonVector[0]->isolationR03().sumPt);
    LeadingMuonsIsolation.push_back(MuonVector[1]->isolationR03().sumPt);
  }

  if(WFillE){
    LeadingElectronsP4.push_back(ElectronVector[0]->p4());
    METP4.push_back(NeutrinoVector[0]->p4());
    METSumEt = NeutrinoVector[0]->sumEt();
    LeadingElectronsIsolation.push_back(ElectronVector[0]->dr03TkSumPt()/ElectronVector[0]->pt());
    LeadingElectronsIsoEcal.push_back(ElectronVector[0]->dr03EcalRecHitSumEt()/ElectronVector[0]->pt());
    LeadingElectronsIsoHcal.push_back(ElectronVector[0]->dr03HcalTowerSumEt()/ElectronVector[0]->pt());
    LeadingElectronsInnerHits.push_back(ElectronVector[0]->gsfTrack()->trackerExpectedHitsInner().numberOfHits());
    LeadingElectronsDCot.push_back(ElectronVector[0]->convDcot());
    LeadingElectronsDist.push_back(ElectronVector[0]->convDist());
    LeadingElectronsDeltaEtaTkClu.push_back(ElectronVector[0]->deltaEtaSuperClusterTrackAtVtx());
    LeadingElectronsDeltaPhiTkClu.push_back(ElectronVector[0]->deltaPhiSuperClusterTrackAtVtx());
    LeadingElectronsSigmaIeIe.push_back(ElectronVector[0]->sigmaIetaIeta());
    LeadingElectronsHE.push_back(ElectronVector[0]->hadronicOverEm());
  }

  if(WFillMu){
    LeadingMuonsP4.push_back(MuonVector[0]->p4());
    METP4.push_back(NeutrinoVector[0]->p4());
    METSumEt = NeutrinoVector[0]->sumEt();
    LeadingMuonsIsolation.push_back(MuonVector[0]->isolationR03().sumPt);
  }

}

template <class T, class W>
math::XYZTLorentzVector BosonAnalyzer::DiSystem(T obj1, W obj2){

  math::XYZTLorentzVector DiObj(0.,0.,0.,0.);
  DiObj += obj1->p4();
  DiObj += obj2->p4();  

  return DiObj;

}


template <class T, class W>
double BosonAnalyzer::MassT(T lepton, W met){

  double tmassWBoson_=-999.;
  double metT = sqrt(pow(met->px(),2) + pow(met->py(),2));
  double lepT = sqrt(pow(lepton->px(),2) + pow(lepton->py(),2));
  reco::Particle::LorentzVector WT = lepton->p4() + met->p4();
  tmassWBoson_ = sqrt(pow(metT+lepT,2) - (WT.px()*WT.px()) - (WT.py()*WT.py()));

  // Defense
  if (!isfinite(tmassWBoson_)) {
    tmassWBoson_ = -999.;
  }
  return tmassWBoson_;
}


template <class T, class W>
double BosonAnalyzer::InvariantMass(T lepton1, W lepton2){

  double massZBoson_=-999.;
  math::XYZTLorentzVector DiSystem(0.,0.,0.,0.);
  DiSystem += lepton1->p4();
  DiSystem += lepton2->p4();
  massZBoson_ = DiSystem.M();

  // Defense
  if (!isfinite(massZBoson_)) {
    massZBoson_ = -999.;
  }
  return massZBoson_;
}


void BosonAnalyzer::PrintOrder(){

  if(PFVector.size()>0){
    for (unsigned int i=0;i<PFVector.size();i++){
      cout << "reco::PF[" << i << "]\t---> pT [GeV]: " << PFVector[i]->pt() << " | eT [GeV]: " << PFVector[i]->et() << " | eta: " << PFVector[i]->eta() << " | phi: " << PFVector[i]->phi() << endl;
    }
  }

  if(protonVector.size()>0){
    for (unsigned int i=0;i<protonVector.size();i++){
      cout << "proton[" << i << "]\t---> pz [GeV]: " << protonVector[i]->pz() << " | pT [GeV]: " << protonVector[i]->pt() << " | eT [GeV]: " << protonVector[i]->et() << " | eta: " << protonVector[i]->eta() << " | phi: " << protonVector[i]->phi() << endl;
    }
  }

  if(ElectronVector.size()>0){
    for (unsigned int i=0;i<ElectronVector.size();i++){
      cout << "reco::Electron[" << i << "]\t---> pT [GeV]: " << ElectronVector[i]->pt() << " | eT [GeV]: " << ElectronVector[i]->et() << " | eta: " << ElectronVector[i]->eta() << " | phi: " << ElectronVector[i]->phi() << endl;
    }
  }

  if(MuonVector.size()>0){
    for (unsigned int i=0;i<MuonVector.size();i++){
      cout << "reco::Muon[" << i << "]\t---> pT [GeV]: " << MuonVector[i]->pt() << " | eT [GeV]: " << MuonVector[i]->et() << " | eta: " << MuonVector[i]->eta() << " | phi: " << MuonVector[i]->phi() << endl;
    }
  }

  if(NeutrinoVector.size()>0){
    for (unsigned int i=0;i<NeutrinoVector.size();i++){
      cout << "reco::Neutrino[" << i << "]\t---> pT [GeV]: " << NeutrinoVector[i]->pt() << " | eT [GeV]: " << NeutrinoVector[i]->et() << " | eta: " << NeutrinoVector[i]->eta() << " | phi: " << NeutrinoVector[i]->phi() << endl;
    }
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(BosonAnalyzer);

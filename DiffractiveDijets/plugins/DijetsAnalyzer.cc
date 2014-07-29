// -*- C++ -*-
//
// Package:    DijetsAnalyzer
// Class:      DijetsAnalyzer
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
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

// Tracks Associated with Jets
#include "DataFormats/JetReco/interface/JetTrackMatch.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/JetReco/interface/TrackJet.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"

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

class DijetsAnalyzer : public edm::EDAnalyzer {
  public:
    explicit DijetsAnalyzer(const edm::ParameterSet&);
    ~DijetsAnalyzer();


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
      double InvariantMass(T jet1, W jet2);

    edm::InputTag particleFlowTag_;
    edm::InputTag VertexTag_;
    double pTPFThresholdCharged_;
    double energyPFThresholdBar_;
    double energyPFThresholdEnd_;
    double energyPFThresholdHF_;
    edm::InputTag GenParticleTag_;
    double EBeam_;
    bool debug_;
    edm::InputTag JetTag_;

    std::vector<const reco::Vertex*> VertexVector;
    std::vector< math::XYZVector > VertexPosition;
    std::vector<const reco::Track*> TracksVector;
    std::vector<const reco::PFCandidate*> PFVector;
    std::vector<const reco::GenParticle*> protonVector;
    std::vector<const reco::PFJet*> JetsVector;

    std::vector< math::XYZTLorentzVector > protonLorentzVector;
    std::vector< math::XYZTLorentzVector > LeadingJetsP4;
    std::vector< math::XYZTLorentzVector > PFP4;

    int nTracks;
    int nVertex;
    double DijetMass;
    double DijetEta;
    double Mpf;

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

DijetsAnalyzer::DijetsAnalyzer(const edm::ParameterSet& iConfig):
  particleFlowTag_(iConfig.getParameter<edm::InputTag>("ParticleFlowTag")),
  VertexTag_(iConfig.getParameter<edm::InputTag>("VertexTag")),
  pTPFThresholdCharged_(iConfig.getParameter<double>("pTPFThresholdCharged")),
  energyPFThresholdBar_(iConfig.getParameter<double>("energyPFThresholdBar")),
  energyPFThresholdEnd_(iConfig.getParameter<double>("energyPFThresholdEnd")),
  energyPFThresholdHF_(iConfig.getParameter<double>("energyPFThresholdHF")),
  GenParticleTag_(iConfig.getParameter<edm::InputTag>("GenParticleTag")),
  EBeam_(iConfig.getParameter<double>("EBeam")),
  debug_(iConfig.getParameter<bool>("debug")),
  JetTag_(iConfig.getParameter<edm::InputTag>("JetTag"))
{

  edm::Service<TFileService> fs;
  eventTree_ = fs->make<TTree>("Event","Event");
  eventTree_->Branch("VertexPosition",&VertexPosition);
  eventTree_->Branch("ProtonsP4",&protonLorentzVector);
  eventTree_->Branch("JetsP4",&LeadingJetsP4);
  eventTree_->Branch("DijetEta",&DijetEta,"DijetEta/D");
  eventTree_->Branch("Mjj",&DijetMass,"DijetMass/D");
  eventTree_->Branch("Mpf",&Mpf,"Mpf/D");
  eventTree_->Branch("nVertex",&nVertex,"nVertex/I");
  eventTree_->Branch("nTracks",&nTracks,"nTracks/I");

}

DijetsAnalyzer::~DijetsAnalyzer()
{
}


// ------------ method called once each job just before starting event loop  ------------
void DijetsAnalyzer::beginJob()
{
  cout << "\n--- C E P   A N A L Y Z E R---\n" << endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void DijetsAnalyzer::endJob()
{

  cout << "\n--- S U M M A R Y---" << endl;
  cout << "\n--- E N D---\n" << endl;

}

// ------------ method called for each event  ------------
void DijetsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  Init(); // Clean Variables
  FillCollections(iEvent, iSetup, debug_);
  if(debug_) PrintOrder();
  SavingInformation();
  eventTree_->Fill();

}

// ------------ Clean Variables --------------

void DijetsAnalyzer::Init(){

  VertexVector.clear();
  VertexPosition.clear();
  TracksVector.clear();
  PFVector.clear();
  protonVector.clear();
  JetsVector.clear();
  protonLorentzVector.clear();
  LeadingJetsP4.clear();
  PFP4.clear();
  nTracks = -999;
  nVertex = -999;
  DijetMass = -999.;
  DijetEta = -999.;
  Mpf = -999.;

}


// ------------ Fill Vectors, All Handles  ------------
void DijetsAnalyzer::FillCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool debug)
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
    math::XYZTLorentzVector allCands(0.,0.,0.,0.);
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
	allCands+=pfAll->p4();
	PFVector.push_back(pfAll);
      }
    }
    Mpf = allCands.M();
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
      if (genAll->pdgId() != 2212) continue;
      protonVector.push_back(genAll);
    }
  }

  // Fill Jets
  Handle<edm::View<reco::PFJet> > jets;
  iEvent.getByLabel(JetTag_,jets);

  int jetsize = jets->size();
  int itJets;

  if(jets->size()>0){
    for(itJets=0; itJets < jetsize; ++itJets){
      const reco::PFJet* jetAll = &((*jets)[itJets]);
      JetsVector.push_back(jetAll);
    }
  }

  // Sorting Vectors, pT order
  std::sort(PFVector.begin(), PFVector.end(), orderETA());
  std::sort(JetsVector.begin(), JetsVector.end(), orderPT());
  std::sort(protonVector.begin(), protonVector.end(), orderAbsolutPZ());

  // Saving One or Two Leading Protons
  if (protonVector.size()==1){
    protonLorentzVector.push_back(protonVector[0]->p4());
  }

  if (protonVector.size()>1){
    protonLorentzVector.push_back(protonVector[0]->p4());
    protonLorentzVector.push_back(protonVector[1]->p4());
  }

  if (JetsVector.size()>1){
    DijetMass = InvariantMass(JetsVector[0],JetsVector[1]);
    DijetEta = DiSystem(JetsVector[0],JetsVector[1]).eta();
  }

  if(debug){
    std::cout << "D I J E T   I N F O" << std::endl;
    std::cout << "Mass: " << DijetMass << " [GeV] | eta: " << DijetEta <<  std::endl;
    if (PFVector.size()>1){
      std::cout << "\nP F   I N F O" << std::endl;
      std::cout << "Eta, max: " << PFVector[0]->eta() << " | Eta, min: " << PFVector[PFVector.size()-1]->eta() <<  std::endl;
    }
  }

}

void DijetsAnalyzer::SavingInformation(){

  if(JetsVector.size()>1 && PFVector.size()>1){
    LeadingJetsP4.push_back(JetsVector[0]->p4());
    LeadingJetsP4.push_back(JetsVector[1]->p4());
    PFP4.push_back(PFVector[0]->p4());
    PFP4.push_back(PFVector[PFVector.size()-1]->p4());
  }

}

template <class T, class W>
math::XYZTLorentzVector DijetsAnalyzer::DiSystem(T obj1, W obj2){

  math::XYZTLorentzVector DiObj(0.,0.,0.,0.);
  DiObj += obj1->p4();
  DiObj += obj2->p4();  

  return DiObj;

}

template <class T, class W>
double DijetsAnalyzer::InvariantMass(T jet1, W jet2){

  double mass_=-999.;
  math::XYZTLorentzVector DiSystem(0.,0.,0.,0.);
  DiSystem += jet1->p4();
  DiSystem += jet2->p4();
  mass_ = DiSystem.M();

  // Defense
  if (!isfinite(mass_)) {
    mass_ = -999.;
  }
  return mass_;
}


void DijetsAnalyzer::PrintOrder(){

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

  if(JetsVector.size()>0){
    for (unsigned int i=0;i<JetsVector.size();i++){
      cout << "reco::PFJet[" << i << "]\t---> pT [GeV]: " << JetsVector[i]->pt() << " | eT [GeV]: " << JetsVector[i]->et() << " | eta: " << JetsVector[i]->eta() << " | phi: " << JetsVector[i]->phi() << endl;
    }
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(DijetsAnalyzer);

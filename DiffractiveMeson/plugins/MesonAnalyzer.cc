// -*- C++ -*-
//
// Package:    MesonAnalyzer
// Class:      MesonAnalyzer
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

class MesonAnalyzer : public edm::EDAnalyzer {
  public:
    explicit MesonAnalyzer(const edm::ParameterSet&);
    ~MesonAnalyzer();


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
    edm::InputTag muonTag_;


    std::vector<const reco::Vertex*> VertexVector;
    std::vector< math::XYZVector > VertexPosition;
    std::vector<const reco::Track*> TracksVector;
    std::vector<const reco::PFCandidate*> PFVector;
    std::vector<const reco::GenParticle*> protonVector;
    std::vector<const reco::Muon*> MuonVector;


    std::vector< math::XYZTLorentzVector > protonLorentzVector;
    std::vector< math::XYZTLorentzVector > LeadingMuonsP4;
 

    /*
    // Possible to use both cases
    std::vector <ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > protonLorentzVector;
    std::vector <ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > LeadingElectronsP4;
    std::vector <ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > LeadingMuonsP4;
    std::vector <ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > METP4;
     */

    int nTracks;
    int nVertex;

    double JpsiMassDiMuon;

    double JpsiEtaDiMuon;


    std::vector<double> LeadingMuonsIsolation;


    bool JpsiMassMu;


    bool JpsiFillMu;


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

MesonAnalyzer::MesonAnalyzer(const edm::ParameterSet& iConfig):
  particleFlowTag_(iConfig.getParameter<edm::InputTag>("ParticleFlowTag")),
  VertexTag_(iConfig.getParameter<edm::InputTag>("VertexTag")),
  pTPFThresholdCharged_(iConfig.getParameter<double>("pTPFThresholdCharged")),
  energyPFThresholdBar_(iConfig.getParameter<double>("energyPFThresholdBar")),
  energyPFThresholdEnd_(iConfig.getParameter<double>("energyPFThresholdEnd")),
  energyPFThresholdHF_(iConfig.getParameter<double>("energyPFThresholdHF")),
  GenParticleTag_(iConfig.getParameter<edm::InputTag>("GenParticleTag")),
  EBeam_(iConfig.getParameter<double>("EBeam")),
  debug_(iConfig.getParameter<bool>("debug")),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonTag"))
 
{

  edm::Service<TFileService> fs;
  eventTree_ = fs->make<TTree>("Event","Event");
  eventTree_->Branch("VertexPosition",&VertexPosition);
  eventTree_->Branch("ProtonsP4",&protonLorentzVector);
  eventTree_->Branch("MuonsP4",&LeadingMuonsP4);
  eventTree_->Branch("LeadingMuonsIsolation",&LeadingMuonsIsolation);
  eventTree_->Branch("JpsiMassDiMuon",&JpsiMassDiMuon,"JpsiMassDiMuon/D");
  eventTree_->Branch("JpsiEtaDiMuon",&JpsiEtaDiMuon,"JpsiEtaDiMuon/D");
  eventTree_->Branch("nVertex",&nVertex,"nVertex/I");
  eventTree_->Branch("nTracks",&nTracks,"nTracks/I");

}

MesonAnalyzer::~MesonAnalyzer()
{
}


// ------------ method called once each job just before starting event loop  ------------
void MesonAnalyzer::beginJob()
{
  cout << "\n--- M E S O N   A N A L Y Z E R---\n" << endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void MesonAnalyzer::endJob()
{

  cout << "\n--- S U M M A R Y---" << endl;
  cout << "\n--- E N D---\n" << endl;

}

// ------------ method called for each event  ------------
void MesonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  Init(); // Clean Variables
  FillCollections(iEvent, iSetup, true); // Debug option, true = show Meson Mass, false = do not show Meson Mass.
  if(debug_) PrintOrder();
  SavingInformation();
  eventTree_->Fill();

}

// ------------ Clean Variables --------------

void MesonAnalyzer::Init(){

  VertexVector.clear();
  VertexPosition.clear();
  TracksVector.clear();
  PFVector.clear();
  protonVector.clear();
  MuonVector.clear();
  protonLorentzVector.clear();
  LeadingMuonsP4.clear();
  LeadingMuonsIsolation.clear();
  

  nTracks = -999;
  nVertex = -999;
  JpsiMassDiMuon = -999.;
  JpsiEtaDiMuon = -999.;

  JpsiMassMu = false;

  JpsiFillMu = false;
 

}


// ------------ Fill Vectors, All Handles  ------------
void MesonAnalyzer::FillCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool debug)
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
      if (genAll->pdgId() != 2212) continue;
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



  // Sorting Vectors, pT order
  std::sort(PFVector.begin(), PFVector.end(), orderPT());
  std::sort(protonVector.begin(), protonVector.end(), orderAbsolutPZ());
  std::sort(MuonVector.begin(), MuonVector.end(), orderPT());
 

  // Saving One or Two Leading Protons
  if (protonVector.size()==1){
    protonLorentzVector.push_back(protonVector[0]->p4());
  }

  if (protonVector.size()>1){
    protonLorentzVector.push_back(protonVector[0]->p4());
    protonLorentzVector.push_back(protonVector[1]->p4());
  }
  //--->


  if (MuonVector.size()>1){  
    JpsiMassDiMuon = InvariantMass(MuonVector[0],MuonVector[1]);
    JpsiEtaDiMuon = DiSystem(MuonVector[0],MuonVector[1]).eta();
  }


  if(debug){
    std::cout << "M A S S" << std::endl;
    std::cout << "Jpsi Mass (MuMu): " << JpsiMassDiMuon << " | eta: " << JpsiEtaDiMuon <<  std::endl;

  }

  if(JpsiMassDiMuon > 2.8 && JpsiMassDiMuon < 3.2) JpsiMassMu = true;

  if(JpsiMassMu) JpsiFillMu = true;


}

void MesonAnalyzer::SavingInformation(){

  if(JpsiFillMu){
    LeadingMuonsP4.push_back(MuonVector[0]->p4());
    LeadingMuonsP4.push_back(MuonVector[1]->p4());
    LeadingMuonsIsolation.push_back(MuonVector[0]->isolationR03().sumPt);
    LeadingMuonsIsolation.push_back(MuonVector[1]->isolationR03().sumPt);
  }

 

}

template <class T, class W>
math::XYZTLorentzVector MesonAnalyzer::DiSystem(T obj1, W obj2){

  math::XYZTLorentzVector DiObj(0.,0.,0.,0.);
  DiObj += obj1->p4();
  DiObj += obj2->p4();  

  return DiObj;

}


template <class T, class W>
double MesonAnalyzer::InvariantMass(T lepton1, W lepton2){

  double massJpsiMeson_=-999.;
  math::XYZTLorentzVector DiSystem(0.,0.,0.,0.);
  DiSystem += lepton1->p4();
  DiSystem += lepton2->p4();
  massJpsiMeson_ = DiSystem.M();

  // Defense
  if (!isfinite(massJpsiMeson_)) {
    massJpsiMeson_ = -999.;
  }
  return massJpsiMeson_;
}


void MesonAnalyzer::PrintOrder(){

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

  if(MuonVector.size()>0){
    for (unsigned int i=0;i<MuonVector.size();i++){
      cout << "reco::Muon[" << i << "]\t---> pT [GeV]: " << MuonVector[i]->pt() << " | eT [GeV]: " << MuonVector[i]->et() << " | eta: " << MuonVector[i]->eta() << " | phi: " << MuonVector[i]->phi() << endl;
    }
  }



}

//define this as a plug-in
DEFINE_FWK_MODULE(MesonAnalyzer);

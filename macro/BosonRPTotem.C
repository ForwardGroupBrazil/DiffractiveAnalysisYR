/*

// Y E L L O W   R E P O R T   M A C R O 
// Author: D. Figueiredo

Goal": obtain some extrapolations

 */

// C++
#include <stdio.h>
#include <math.h> 
#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <cmath>

// root
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TLegend.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TMatrixT.h>
#include "TMath.h"
#include "TLorentzVector.h"

#define PI 3.141592653589793

void BosonRPTotem()
{
  double lumi = 100; // pb
  double weight = pow(lumi,-1);

  bool debug = false; // print code output

  //------------ FWLite libraries ------
  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable();
  gROOT->ProcessLine("#include<vector>");

  //------------ files -----------------
  TFile *inf  = TFile::Open("out.root");
  TFile *out = TFile::Open("OutHisto.root","RECREATE");
  TTree *tr = (TTree*)inf->Get("BosonAnalyzer/Event");

  TFile *RPFileCMSMinus = TFile::Open("matrix_xi_vs_t_right.root");
  TFile *RPFileCMSPlus  = TFile::Open("matrix_xi_vs_t_left.root");

  TH2F *HistoRPCMSMinus = (TH2F*)RPFileCMSMinus->Get("accep_xi_vs_t_right");
  TH2F *HistoRPCMSPlus = (TH2F*)RPFileCMSPlus->Get("accep_xi_vs_t_left");

  //----------- define the tree branch --------
  std::vector <ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *protonLorentzVector;
  std::vector <ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *LeadingElectronsP4;
  std::vector <ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *LeadingMuonsP4;
  std::vector <ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *METP4;
  std::vector<double> *LeadingElectronsIsolation;
  std::vector<double> *LeadingElectronsIsoEcal;
  std::vector<double> *LeadingElectronsIsoHcal;
  std::vector<double> *LeadingElectronsInnerHits;
  std::vector<double> *LeadingElectronsDCot;
  std::vector<double> *LeadingElectronsDist;
  std::vector<double> *LeadingElectronsDeltaEtaTkClu;
  std::vector<double> *LeadingElectronsDeltaPhiTkClu;
  std::vector<double> *LeadingElectronsSigmaIeIe;
  std::vector<double> *LeadingElectronsHE;
  std::vector<double> *LeadingMuonsIsolation;

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

  TBranch *bprotonLorentzVector;
  TBranch *bLeadingElectronsP4;
  TBranch *bLeadingMuonsP4;
  TBranch *bMETP4;

  TBranch *bLeadingElectronsIsolation;
  TBranch *bLeadingElectronsIsoEcal;
  TBranch *bLeadingElectronsIsoHcal;
  TBranch *bLeadingElectronsInnerHits;
  TBranch *bLeadingElectronsDCot;
  TBranch *bLeadingElectronsDist;
  TBranch *bLeadingElectronsDeltaEtaTkClu;
  TBranch *bLeadingElectronsDeltaPhiTkClu;
  TBranch *bLeadingElectronsSigmaIeIe;
  TBranch *bLeadingElectronsHE;
  TBranch *bLeadingMuonsIsolation;

  TBranch *bMETSumEt;
  TBranch *bnTracks;
  TBranch *bnVertex;
  TBranch *bZMassDiElectron;
  TBranch *bWMassElectron;
  TBranch *bZMassDiMuon;
  TBranch *bWMassMuon;
  TBranch *bZEtaDiElectron;
  TBranch *bZEtaDiMuon;
  TBranch *bWEtaElectron;
  TBranch *bWEtaMuon;

  tr->SetBranchAddress("ProtonsP4",&protonLorentzVector,&bprotonLorentzVector);
  tr->SetBranchAddress("ElectronsP4",&LeadingElectronsP4,&bLeadingElectronsP4);
  tr->SetBranchAddress("MuonsP4",&LeadingMuonsP4,&bLeadingMuonsP4);
  tr->SetBranchAddress("METP4",&METP4,&bMETP4);
  tr->SetBranchAddress("LeadingElectronsIsolation",&LeadingElectronsIsolation,&bLeadingElectronsIsolation);
  tr->SetBranchAddress("LeadingElectronsIsoEcal",&LeadingElectronsIsoEcal,&bLeadingElectronsIsoEcal);
  tr->SetBranchAddress("LeadingElectronsIsoHcal",&LeadingElectronsIsoHcal,&bLeadingElectronsIsoHcal);
  tr->SetBranchAddress("LeadingElectronsInnerHits",&LeadingElectronsInnerHits,&bLeadingElectronsInnerHits);
  tr->SetBranchAddress("LeadingElectronsDCot",&LeadingElectronsDCot,&bLeadingElectronsDCot);
  tr->SetBranchAddress("LeadingElectronsDist",&LeadingElectronsDist,&bLeadingElectronsDist);
  tr->SetBranchAddress("LeadingElectronsDeltaEtaTkClu",&LeadingElectronsDeltaEtaTkClu,&bLeadingElectronsDeltaEtaTkClu);
  tr->SetBranchAddress("LeadingElectronsDeltaPhiTkClu",&LeadingElectronsDeltaPhiTkClu,&bLeadingElectronsDeltaPhiTkClu);
  tr->SetBranchAddress("LeadingElectronsSigmaIeIe",&LeadingElectronsSigmaIeIe,&bLeadingElectronsSigmaIeIe);
  tr->SetBranchAddress("LeadingElectronsHE",&LeadingElectronsHE,&bLeadingElectronsHE);
  tr->SetBranchAddress("LeadingMuonsIsolation",&LeadingMuonsIsolation,&bLeadingMuonsIsolation);
  tr->SetBranchAddress("ZMassDiElectron",&ZMassDiElectron,&bZMassDiElectron);
  tr->SetBranchAddress("WMassElectron",&WMassElectron,&bWMassElectron);
  tr->SetBranchAddress("ZMassDiMuon",&ZMassDiMuon,&bZMassDiMuon);
  tr->SetBranchAddress("WMassMuon",&WMassMuon,&bWMassMuon);
  tr->SetBranchAddress("ZEtaDiElectron",&ZEtaDiElectron,&bZEtaDiElectron);
  tr->SetBranchAddress("ZEtaDiMuon",&ZEtaDiMuon,&bZEtaDiMuon);
  tr->SetBranchAddress("WEtaElectron",&WEtaElectron,&bWEtaElectron);
  tr->SetBranchAddress("WEtaMuon",&WEtaMuon,&bWEtaMuon);
  tr->SetBranchAddress("nVertex",&nVertex,&bnVertex);
  tr->SetBranchAddress("nTracks",&nTracks,&bnTracks);

  // Creating Histograms
  std::vector<TH1F*> hVectorProtonEta;
  std::vector<TH1F*> hVectorProtonEnergy;
  std::vector<TH1F*> hVectorProtonPz;
  std::vector<TH1F*> hVectorProtonXi;
  std::vector<TH1F*> hVectorProtonT;
  std::vector<TH1F*> hVectorBosonM;
  std::vector<TH1F*> hVectorBosonEta;
  std::vector<TH1F*> hVectorLeptonsEta;
  std::vector<TH1F*> hVectorLeptonsPhi;
  std::vector<TH1F*> hVectorLeptonsPt;
  std::vector<TH1F*> hVectorCrossSection;

  std::string step0 = "no_accept_RP";
  std::string step1 = "accept_RP_boson";

  std::vector <std::string> GroupHisto;

  GroupHisto.push_back(step0);
  GroupHisto.push_back(step1);

  for (std::vector<std::string>::size_type j=0; j<GroupHisto.size(); j++){
    char name[300];
    sprintf(name,"ProtonEta_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonEta = new TH1F(name,"#eta; N events",100,-15,15);
    hVectorProtonEta.push_back(hProtonEta);

    sprintf(name,"ProtonPz_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonPz = new TH1F(name,"#P_{Z} [GeV]; N events",14000,-7000.,7000.);
    hVectorProtonPz.push_back(hProtonPz);

    sprintf(name,"ProtonEnergy_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonEnergy = new TH1F(name,"Energy [GeV]; N events",7000,0.,7000.);
    hVectorProtonEnergy.push_back(hProtonEnergy);

    sprintf(name,"ProtonXi_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonXi = new TH1F(name,"#xi; N events",100,0.,1.);
    hVectorProtonXi.push_back(hProtonXi);

    sprintf(name,"ProtonT_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonT = new TH1F(name,"|t|; N events",100,0.,1.);
    hVectorProtonT.push_back(hProtonT);

    sprintf(name,"BosonM_%s",GroupHisto.at(j).c_str());
    TH1F *hBosonM = new TH1F(name,"Mass; N events",500,0.,500.);
    hVectorBosonM.push_back(hBosonM);

    sprintf(name,"BosonEta_%s",GroupHisto.at(j).c_str());
    TH1F *hBosonEta = new TH1F(name,"#eta; N events",100,-6.,6.);
    hVectorBosonEta.push_back(hBosonEta);

    sprintf(name,"LeptonsEta_%s",GroupHisto.at(j).c_str());
    TH1F *hLeptonsEta = new TH1F(name,"#eta; N events",100,-6.,6.);
    hVectorLeptonsEta.push_back(hLeptonsEta);

    sprintf(name,"LeptonsPhi_%s",GroupHisto.at(j).c_str());
    TH1F *hLeptonsPhi = new TH1F(name,"#phi; N events",100,-4.,4.);
    hVectorLeptonsPhi.push_back(hLeptonsPhi);

    sprintf(name,"LeptonsPt_%s",GroupHisto.at(j).c_str());
    TH1F *hLeptonsPt = new TH1F(name,"p_{T}; N events",400,0.,400.);
    hVectorLeptonsPt.push_back(hLeptonsPt);

    sprintf(name,"CrossSection_%s",GroupHisto.at(j).c_str());
    TH1F *hCrossSection = new TH1F(name,"d#sigma/dM; N events",400,0.,400.);
    hVectorCrossSection.push_back(hCrossSection);

  }

  unsigned NEntries = tr->GetEntries();
  cout << "Reading TREE: "<<NEntries<<" events"<<endl;

  for(int unsigned i=0; i<NEntries; i++) {

    tr->GetEntry(i);

    double perc = -999;
    double xi_proton_plus = -999.;
    double xi_proton_minus = -999.;
    double t_proton_plus = -999.;
    double t_proton_minus = -999.;
    double accept = -999.;

    bool eleBarrel1 = false;
    bool eleEndCap1 = false;
    bool eleBarrel2 = false;
    bool eleEndCap2 = false;
    bool candSel = false;
    bool isoBarrel1 = false;
    bool isoEndCap1 = false;
    bool isoBarrel2 = false;
    bool isoEndCap2 = false;
    bool isolation = false;

    bool ZMassMu = false;
    bool WMassMu = false;
    bool ZMassE = false;
    bool WMassE = false;

    bool ZFillMu = false;
    bool WFillMu = false;
    bool ZFillE = false;
    bool WFillE = false;

    if (debug) std::cout << ">> EVENT " << i << " <<"  <<endl;

    if (protonLorentzVector->size() == 1){
      if (debug) cout << "--> proton(1) pZ [GeV]: " << protonLorentzVector->at(0).pz() << endl;
      perc = protonLorentzVector->at(0).pz()/6500.;
    }

    if (protonLorentzVector->size() > 1){
      if(debug){
	cout << "--> proton(1) pZ [GeV]: " << protonLorentzVector->at(0).pz() << endl;
	cout << "--> proton(2) pZ [GeV]: " << protonLorentzVector->at(1).pz() << endl;
      }
      perc = protonLorentzVector->at(0).pz()/6500.;
    }

    if (protonLorentzVector->size() > 0 && protonLorentzVector->at(0).pz() > 0. && fabs(perc) > 0.75){
      xi_proton_plus =  ( 1. - (protonLorentzVector->at(0).pz()/6500.) );

      TLorentzVector vec_pi(0.,0.,6500.,6500.);
      TLorentzVector vec_pf(protonLorentzVector->at(0).px(),protonLorentzVector->at(0).py(),protonLorentzVector->at(0).pz(),protonLorentzVector->at(0).energy());
      TLorentzVector vec_t = (vec_pf - vec_pi);

      t_proton_plus = vec_t.Mag2();
      accept = HistoRPCMSPlus->GetBinContent(HistoRPCMSPlus->GetXaxis()->FindBin(fabs(t_proton_plus)),HistoRPCMSPlus->GetYaxis()->FindBin(xi_proton_plus));

      if(debug){
	cout << "--> xi, plus: " << xi_proton_plus << endl;
	cout << "-->  t, plus: " << t_proton_plus << endl;
	cout << "--> RP Acceptance: " << accept << endl;
      }

      hVectorProtonEta.at(0)->Fill(protonLorentzVector->at(0).eta());
      hVectorProtonPz.at(0)->Fill(protonLorentzVector->at(0).pz());
      hVectorProtonEnergy.at(0)->Fill(protonLorentzVector->at(0).energy());
      hVectorProtonXi.at(0)->Fill(xi_proton_plus);
      hVectorProtonT.at(0)->Fill(fabs(t_proton_plus));

    }

    if (protonLorentzVector->size() > 0 && protonLorentzVector->at(0).pz() < 0. && fabs(perc) > 0.75){
      xi_proton_minus =  ( 1. + (protonLorentzVector->at(0).pz()/6500.) );

      TLorentzVector vec_pi(0.,0.,-6500.,6500.);
      TLorentzVector vec_pf(protonLorentzVector->at(0).px(),protonLorentzVector->at(0).py(),protonLorentzVector->at(0).pz(),protonLorentzVector->at(0).energy());
      TLorentzVector vec_t = (vec_pf - vec_pi);

      t_proton_minus = vec_t.Mag2();

      accept = HistoRPCMSMinus->GetBinContent(HistoRPCMSMinus->GetXaxis()->FindBin(fabs(t_proton_minus)),HistoRPCMSMinus->GetYaxis()->FindBin(xi_proton_minus));

      if(debug){
	cout << "--> xi, minus: " << xi_proton_minus << endl;
	cout << "-->  t, minus: " << t_proton_minus << endl;
	cout << "--> RP Acceptance: " << accept << endl;
      }

      hVectorProtonEta.at(0)->Fill(protonLorentzVector->at(0).eta());
      hVectorProtonPz.at(0)->Fill(protonLorentzVector->at(0).pz());
      hVectorProtonEnergy.at(0)->Fill(protonLorentzVector->at(0).energy());
      hVectorProtonXi.at(0)->Fill(xi_proton_minus);
      hVectorProtonT.at(0)->Fill(fabs(t_proton_minus));
    } 

    if (debug) cout << "-- END --\n" << endl;

    if (ZMassDiElectron > 60. && ZMassDiElectron < 110. ) ZMassE = true;
    if (ZMassDiMuon > 60. && ZMassDiMuon < 110. ) ZMassMu = true;
    if (WMassElectron > 60. && WMassElectron < 110.) WMassE = true;
    if (WMassMuon > 60. && WMassMuon < 110.) WMassMu = true;

    if(ZMassE && !ZMassMu && !WMassE && !WMassMu) ZFillE = true;
    if(ZMassMu && !ZMassE && !WMassE && !WMassMu) ZFillMu = true;
    if(WMassE && !ZMassE && !ZMassMu && !WMassMu) WFillE = true;
    if(WMassMu && !ZMassE && !ZMassMu && !WMassE) WFillMu = true;

    // Isolation One Leading Electrons
    if(WFillE){

      if ((fabs (LeadingElectronsP4->at(0).eta()) <= 1.4442) ){
	if (LeadingElectronsIsolation->at(0)<0.09 && LeadingElectronsIsoEcal->at(0)<0.07 && LeadingElectronsIsoHcal->at(0)<0.10) isoBarrel1 = true;
      }

      if ((fabs (LeadingElectronsP4->at(0).eta()) >= 1.5660) && (fabs (LeadingElectronsP4->at(0).eta()) <= 2.5)){
	if (LeadingElectronsIsolation->at(0)<0.04 && LeadingElectronsIsoEcal->at(0)<0.05 && LeadingElectronsIsoHcal->at(0)<0.025) isoEndCap1 = true;
      }

      if ((fabs (LeadingElectronsP4->at(0).eta()) <= 1.4442) ){
	if (LeadingElectronsInnerHits->at(0) == 0 && (fabs (LeadingElectronsDCot[0]) >= 0.02 || fabs (LeadingElectronsDist->at(0)) >= 0.02 ) && fabs (LeadingElectronsDeltaEtaTkClu->at(0)) < 0.004 && fabs (LeadingElectronsDeltaPhiTkClu->at(0)) < 0.06 && LeadingElectronsSigmaIeIe->at(0) < 0.01 && LeadingElectronsHE->at(0) < 0.04 ) eleBarrel1 = true;
      }

      if ((fabs (LeadingElectronsP4->at(0).eta()) >= 1.5660) && (fabs (LeadingElectronsP4->at(0).eta()) <= 2.5)){
	if (LeadingElectronsInnerHits->at(0) == 0 && (fabs (LeadingElectronsDCot[0]) >= 0.02 || fabs (LeadingElectronsDist->at(0)) >= 0.02 ) && fabs (LeadingElectronsDeltaEtaTkClu->at(0)) < 0.007 && fabs (LeadingElectronsDeltaPhiTkClu->at(0)) < 0.03 && LeadingElectronsSigmaIeIe->at(0) < 0.03 && LeadingElectronsHE->at(0) < 0.025) eleEndCap1 = true;
      }

      if (eleEndCap1 || eleBarrel1) candSel = true;

    }

    // Isolation Two Leading Electrons
    if(ZFillE){

      if ((fabs (LeadingElectronsP4->at(0).eta()) <= 1.4442) ){
	if (LeadingElectronsIsolation->at(0)<0.09 && LeadingElectronsIsoEcal->at(0)<0.07 && LeadingElectronsIsoHcal->at(0)<0.10) isoBarrel1 = true;
      }

      if ((fabs (LeadingElectronsP4->at(0).eta()) >= 1.5660) && (fabs (LeadingElectronsP4->at(0).eta()) <= 2.5)){
	if (LeadingElectronsIsolation->at(0)<0.04 && LeadingElectronsIsoEcal->at(0)<0.05 && LeadingElectronsIsoHcal->at(0)<0.025) isoEndCap1 = true;
      }

      if ((fabs (LeadingElectronsP4->at(1).eta()) <= 1.4442) ){
	if (LeadingElectronsIsolation->at(1)<0.09 && LeadingElectronsIsoEcal->at(1)<0.07 && LeadingElectronsIsoHcal->at(1)<0.10) isoBarrel2 = true;
      }

      if ((fabs (LeadingElectronsP4->at(1).eta()) >= 1.5660) && (fabs (LeadingElectronsP4->at(1).eta()) <= 2.5)){
	if (LeadingElectronsIsolation->at(1)<0.04 && LeadingElectronsIsoEcal->at(1)<0.05 && LeadingElectronsIsoHcal->at(1)<0.025) isoEndCap2 = true;
      }

      if ((isoEndCap1 || isoBarrel1) && (isoEndCap2 || isoBarrel2)) isolation = true;

      if ((fabs (LeadingElectronsP4->at(0).eta()) <= 1.4442) ){
	if (LeadingElectronsInnerHits->at(0) == 0 && (fabs (LeadingElectronsDCot[0]) >= 0.02 || fabs (LeadingElectronsDist->at(0)) >= 0.02 ) && fabs (LeadingElectronsDeltaEtaTkClu->at(0)) < 0.004 && fabs (LeadingElectronsDeltaPhiTkClu->at(0)) < 0.06 && LeadingElectronsSigmaIeIe->at(0) < 0.01 && LeadingElectronsHE->at(0) < 0.04 ) eleBarrel1 = true;
      }

      if ((fabs (LeadingElectronsP4->at(0).eta()) >= 1.5660) && (fabs (LeadingElectronsP4->at(0).eta()) <= 2.5)){
	if (LeadingElectronsInnerHits->at(0) == 0 && (fabs (LeadingElectronsDCot[0]) >= 0.02 || fabs (LeadingElectronsDist->at(0)) >= 0.02 ) && fabs (LeadingElectronsDeltaEtaTkClu->at(0)) < 0.007 && fabs (LeadingElectronsDeltaPhiTkClu->at(0)) < 0.03 && LeadingElectronsSigmaIeIe->at(0) < 0.03 && LeadingElectronsHE->at(0) < 0.025) eleEndCap1 = true;
      }

      if ((fabs (LeadingElectronsP4->at(1).eta()) <= 1.4442) ){
	if (LeadingElectronsInnerHits->at(1) == 0 && (fabs (LeadingElectronsDCot->at(1)) >= 0.02 || fabs (LeadingElectronsDist->at(1)) >= 0.02 ) && fabs (LeadingElectronsDeltaEtaTkClu->at(1)) < 0.004 && fabs (LeadingElectronsDeltaPhiTkClu->at(1)) < 0.06 && LeadingElectronsSigmaIeIe->at(1) < 0.01 && LeadingElectronsHE->at(1) < 0.04 ) eleBarrel2 = true;
      }

      if ((fabs (LeadingElectronsP4->at(1).eta()) >= 1.5660) && (fabs (LeadingElectronsP4->at(1).eta()) <= 2.5)){
	if (LeadingElectronsInnerHits->at(1) == 0 && (fabs (LeadingElectronsDCot->at(1)) >= 0.02 || fabs (LeadingElectronsDist->at(1)) >= 0.02 ) && fabs (LeadingElectronsDeltaEtaTkClu->at(1)) < 0.007 && fabs (LeadingElectronsDeltaPhiTkClu->at(1)) < 0.03 && LeadingElectronsSigmaIeIe->at(1) < 0.03 && LeadingElectronsHE->at(1) < 0.025) eleEndCap2 = true;
      }

      if ((eleEndCap1 || eleBarrel1) && (eleEndCap2 || eleBarrel2)) candSel = true;
    }

    // Isolation One Leading Muon
    if(WFillMu){
      if (LeadingMuonsIsolation->at(0) < 3) {
	isolation = true;
	candSel = true;
      }
    }

    // Isolation Two Leading Muons
    if(ZFillMu){
      if (LeadingMuonsIsolation->at(0) < 3 && LeadingMuonsIsolation->at(1) < 3 ) {
	isolation = true;
	candSel = true;
      }
    }

    if(accept > 0 && isolation && candSel && (ZMassMu || ZMassE) && LeadingElectronsP4->size()>1){
      hVectorProtonEta.at(1)->Fill(protonLorentzVector->at(0).eta(),accept);
      hVectorProtonPz.at(1)->Fill(protonLorentzVector->at(0).pz(),accept);
      hVectorProtonEnergy.at(1)->Fill(protonLorentzVector->at(0).energy(),accept);
      if(xi_proton_minus > -999.) hVectorProtonXi.at(1)->Fill(xi_proton_minus,accept);
      if(xi_proton_plus > -999.) hVectorProtonXi.at(1)->Fill(xi_proton_plus,accept);
      if(t_proton_minus > -999.) hVectorProtonT.at(1)->Fill(fabs(t_proton_minus),accept);
      if(t_proton_plus > -999.) hVectorProtonT.at(1)->Fill(fabs(t_proton_plus),accept);
      if (ZMassDiElectron > ZMassDiMuon) hVectorBosonM.at(1)->Fill(ZMassDiElectron,accept);
      if (ZMassDiElectron < ZMassDiMuon)  hVectorBosonM.at(1)->Fill(ZMassDiMuon,accept);
      hVectorBosonEta.at(1)->Fill(ZEtaDiElectron,accept);
      hVectorBosonEta.at(1)->Fill(ZEtaDiMuon,accept);
      hVectorLeptonsEta.at(1)->Fill(LeadingElectronsP4->at(0).eta(),accept);
      hVectorLeptonsEta.at(1)->Fill(LeadingElectronsP4->at(1).eta(),accept);
      hVectorLeptonsPhi.at(1)->Fill(LeadingElectronsP4->at(0).phi(),accept);
      hVectorLeptonsPhi.at(1)->Fill(LeadingElectronsP4->at(1).phi(),accept);
      hVectorLeptonsPt.at(1)->Fill(LeadingElectronsP4->at(0).pt(),accept);
      hVectorLeptonsPt.at(1)->Fill(LeadingElectronsP4->at(1).pt(),accept);
      if (ZMassDiElectron > ZMassDiMuon) hVectorCrossSection.at(1)->Fill(ZMassDiElectron,accept);
      if (ZMassDiElectron < ZMassDiMuon) hVectorCrossSection.at(1)->Fill(ZMassDiMuon,accept);
    }
  }

  hVectorCrossSection[1]->Scale(weight,"width");

  out->cd();
  hVectorProtonEta[0]->Write();
  hVectorProtonPz[0]->Write();
  hVectorProtonEnergy[0]->Write();
  hVectorProtonXi[0]->Write();
  hVectorProtonT[0]->Write();
  hVectorProtonEta[1]->Write();
  hVectorProtonPz[1]->Write();
  hVectorProtonEnergy[1]->Write();
  hVectorProtonXi[1]->Write();
  hVectorProtonT[1]->Write();
  hVectorBosonM[1]->Write();
  hVectorBosonEta[1]->Write();
  hVectorLeptonsEta[1]->Write();
  hVectorLeptonsPhi[1]->Write();
  hVectorLeptonsPt[1]->Write();
  hVectorCrossSection[1]->Write();

  out->Close();
  inf->Close();
  RPFileCMSMinus->Close();
  RPFileCMSPlus->Close();

}

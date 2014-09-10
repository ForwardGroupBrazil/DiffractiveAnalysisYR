/*

// Y E L L O W   R E P O R T   M A C R O 

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
#include <fstream>

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

void MesonRPTotem(string inputfile, string outputfile, double XSmcJpsi, double lumi)
{

  //Parameters
  // Lumi = luminosity [pb-1]
  // XSmcJpsi = cross section Jpsi = XSmc*BR*<s2>*factor flux; BR=6%, <s2>=10%, ff=7.7

  double weight = pow(lumi,-1);
  double EBeam = 6500.;

  bool debug = false; // print code output

  //------------ FWLite libraries ------
  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable();
  gROOT->ProcessLine("#include<vector>");

  //------------ files -----------------
  TFile *inf  = TFile::Open(inputfile.c_str());
  TFile *out = TFile::Open(outputfile.c_str(),"RECREATE");
  TTree *tr = (TTree*)inf->Get("mesonAnalyzer/Event");

  TString outtxt = outputfile;
  outtxt.ReplaceAll("root","txt");
  std::ofstream outstring(outtxt);

  TFile *RPFileCMSMinus = TFile::Open("matrix_xi_vs_t_right.root");
  TFile *RPFileCMSPlus  = TFile::Open("matrix_xi_vs_t_left.root");

  TH2F *HistoRPCMSMinus = (TH2F*)RPFileCMSMinus->Get("accep_xi_vs_t_right");
  TH2F *HistoRPCMSPlus = (TH2F*)RPFileCMSPlus->Get("accep_xi_vs_t_left");

  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  //----------- define the tree branch --------
  std::vector <ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *protonLorentzVector;
  std::vector <ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *LeadingMuonsP4;
  std::vector<double> *LeadingMuonsIsolation;


  int nTracks;
  int nVertex;


  double JpsiMassDiMuon;
  double JpsiEtaDiMuon;
  double AcceptJpsi = 0.;
  double selectedJpsi = 0.;

  TBranch *bprotonLorentzVector;
  TBranch *bLeadingMuonsP4;

  TBranch *bLeadingMuonsIsolation;

  TBranch *bnTracks;
  TBranch *bnVertex;
  TBranch *bJpsiMassDiMuon;
  TBranch *bJpsiEtaDiMuon;

  tr->SetBranchAddress("ProtonsP4",&protonLorentzVector,&bprotonLorentzVector);
  tr->SetBranchAddress("MuonsP4",&LeadingMuonsP4,&bLeadingMuonsP4);
  tr->SetBranchAddress("LeadingMuonsIsolation",&LeadingMuonsIsolation,&bLeadingMuonsIsolation);
  tr->SetBranchAddress("JpsiMassDiMuon",&JpsiMassDiMuon,&bJpsiMassDiMuon);
  tr->SetBranchAddress("JpsiEtaDiMuon",&JpsiEtaDiMuon,&bJpsiEtaDiMuon);
  tr->SetBranchAddress("nVertex",&nVertex,&bnVertex);
  tr->SetBranchAddress("nTracks",&nTracks,&bnTracks);

  // Creating Histograms
  std::vector<TH1F*> hVectorProtonEta;
  std::vector<TH1F*> hVectorProtonEnergy;
  std::vector<TH1F*> hVectorProtonPz;
  std::vector<TH1F*> hVectorProtonXi;
  std::vector<TH1F*> hVectorProtonT;
  std::vector<TH1F*> hVectorMesonM;
  std::vector<TH1F*> hVectorMesonEta;
  std::vector<TH1F*> hVectorLeptonsEta;
  std::vector<TH1F*> hVectorLeptonsPhi;
  std::vector<TH1F*> hVectorLeptonsPt;
  std::vector<TH1F*> hVectorCrossSection;
  std::vector<TH1F*> hVectorProtonAcceptanceXiMinus;
  std::vector<TH1F*> hVectorProtonAcceptanceXiPlus;
  std::vector<TH1F*> hVectorProtonAcceptanceTMinus;
  std::vector<TH1F*> hVectorProtonAcceptanceTPlus;
  std::vector<TH1F*> hVectorAccept;
  std::vector<TH1F*> hVectorVertex;

  std::string step0 = "no_accept_RP";
  std::string step1 = "accept_RP_meson";

  std::vector <std::string> GroupHisto;

  GroupHisto.push_back(step0);
  GroupHisto.push_back(step1);

  for (std::vector<std::string>::size_type j=0; j<GroupHisto.size(); j++){
    char name[300];
    sprintf(name,"ProtonEta_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonEta = new TH1F(name,";#eta; N events",100,-15,15);
    hVectorProtonEta.push_back(hProtonEta);

    sprintf(name,"ProtonPz_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonPz = new TH1F(name,";p_{z} [GeV]; N events",14000,-7000.,7000.);
    hVectorProtonPz.push_back(hProtonPz);

    sprintf(name,"ProtonEnergy_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonEnergy = new TH1F(name,";Energy [GeV]; N events",7000,0.,7000.);
    hVectorProtonEnergy.push_back(hProtonEnergy);

    sprintf(name,"ProtonXi_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonXi = new TH1F(name,";#xi; N events",50,0.,0.2);
    hVectorProtonXi.push_back(hProtonXi);

    sprintf(name,"ProtonT_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonT = new TH1F(name,";|t|; N events",50,0.,1);
    hVectorProtonT.push_back(hProtonT);

    sprintf(name,"MesonM_%s",GroupHisto.at(j).c_str());
    TH1F *hMesonM = new TH1F(name,";Mass; N events",1000,0.,10.);
    hVectorMesonM.push_back(hMesonM);

    sprintf(name,"MesonEta_%s",GroupHisto.at(j).c_str());
    TH1F *hMesonEta = new TH1F(name,";#eta; N events",100,-6.,6.);
    hVectorMesonEta.push_back(hMesonEta);

    sprintf(name,"LeptonsEta_%s",GroupHisto.at(j).c_str());
    TH1F *hLeptonsEta = new TH1F(name,";#eta; N events",100,-6.,6.);
    hVectorLeptonsEta.push_back(hLeptonsEta);

    sprintf(name,"LeptonsPhi_%s",GroupHisto.at(j).c_str());
    TH1F *hLeptonsPhi = new TH1F(name,";#phi; N events",100,-4.,4.);
    hVectorLeptonsPhi.push_back(hLeptonsPhi);

    sprintf(name,"LeptonsPt_%s",GroupHisto.at(j).c_str());
    TH1F *hLeptonsPt = new TH1F(name,";p_{T} [GeV]; N events",400,0.,400.);
    hVectorLeptonsPt.push_back(hLeptonsPt);

    sprintf(name,"CrossSection_%s",GroupHisto.at(j).c_str());
    TH1F *hCrossSection = new TH1F(name,";d#sigma/dM; N events",100,0.,10.);
    hVectorCrossSection.push_back(hCrossSection);

    sprintf(name,"ProtonAcceptanceXiMinus_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonAccXiMinus = new TH1F(name,";#xi; N events",50,0.,1.);
    hVectorProtonAcceptanceXiMinus.push_back(hProtonAccXiMinus);

    sprintf(name,"ProtonAcceptanceXiPlus_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonAccXiPlus = new TH1F(name,";#xi; N events",50,0.,1.);
    hVectorProtonAcceptanceXiPlus.push_back(hProtonAccXiPlus);

    sprintf(name,"ProtonAcceptanceTMinus_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonAccTMinus = new TH1F(name,";|t|; N events",50,0.,1.);
    hVectorProtonAcceptanceTMinus.push_back(hProtonAccTMinus);

    sprintf(name,"ProtonAcceptanceTPlus_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonAccTPlus = new TH1F(name,";|t|; N events",50,0.,1.);
    hVectorProtonAcceptanceTPlus.push_back(hProtonAccTPlus);

    sprintf(name,"Acceptance_%s",GroupHisto.at(j).c_str());
    TH1F *hAccept = new TH1F(name,";N_{RP}/N_{GEN}; N events",100,0.,1.);
    hVectorAccept.push_back(hAccept);

    sprintf(name,"VertexMultiplicity_%s",GroupHisto.at(j).c_str());
    TH1F *hVertex = new TH1F(name,";Vertex Multiplicity; N events",10,0,10);
    hVectorVertex.push_back(hVertex);

  }

  //unsigned NEntries = tr->GetEntries();
  unsigned NEntries = 100000;

  cout << "\nR U N N I N G" << endl;
  cout << "Reading TREE: "<<NEntries<<" events"<<endl;
  cout << "Input filename: " << inputfile << endl;
  cout << "Output filename: " << outputfile << endl;
  cout << "Integrated Luminosity: " << lumi << endl;
  cout << "Cross Section MC Jpsi(mumu): " << XSmcJpsi << "\n" << endl;

  hVectorProtonAcceptanceTPlus.at(0)=HistoRPCMSPlus->ProjectionX();
  hVectorProtonAcceptanceTMinus.at(0)=HistoRPCMSMinus->ProjectionX();
  hVectorProtonAcceptanceXiPlus.at(0)=HistoRPCMSPlus->ProjectionY();
  hVectorProtonAcceptanceXiMinus.at(0)=HistoRPCMSMinus->ProjectionY();

  hVectorProtonAcceptanceTPlus.at(0)->SetTitle("");
  hVectorProtonAcceptanceTMinus.at(0)->SetTitle("");
  hVectorProtonAcceptanceXiPlus.at(0)->SetTitle("");
  hVectorProtonAcceptanceXiMinus.at(0)->SetTitle("");

  hVectorProtonAcceptanceTPlus.at(0)->GetXaxis()->SetTitle("|t|, RP^{CMS, +} acceptance");
  hVectorProtonAcceptanceTMinus.at(0)->GetXaxis()->SetTitle("|t|, RP^{CMS, -} acceptance");
  hVectorProtonAcceptanceXiPlus.at(0)->GetXaxis()->SetTitle("#xi, RP^{CMS, +} acceptance");
  hVectorProtonAcceptanceXiMinus.at(0)->GetXaxis()->SetTitle("#xi, RP^{CMS, -} acceptance");

  hVectorProtonAcceptanceTPlus.at(0)->SetName("TMatrixAcceptanceRPPlus");
  hVectorProtonAcceptanceTMinus.at(0)->SetName("TMatrixAcceptanceRPMinus");
  hVectorProtonAcceptanceXiPlus.at(0)->SetName("XiMatrixAcceptanceRPPlus");
  hVectorProtonAcceptanceXiMinus.at(0)->SetName("XiMatrixAcceptanceRPMinus");

  int decade = 0;

  for(int unsigned i=0; i<NEntries; i++) {

    tr->GetEntry(i);

    //cout << "Event " << i << endl;

    double progress = 10.0*i/(1.0*NEntries);
    Int_t k = TMath::FloorNint(progress);
    if (k > decade)
      cout<<10*k<<" %"<<endl;
    decade = k; 

    double perc = -999;
    double xi_proton_plus = -999.;
    double xi_proton_minus = -999.;
    double t_proton_plus = -999.;
    double t_proton_minus = -999.;
    double accept = -999.;


    bool isolation = false;
    bool candSel = false;

    bool JpsiMassMu = false;

    bool JpsiFillMu = false;


    bool protonMinus = false;
    bool protonPlus = false;


    bool acceptJpsiMu = false;


    bool SingleVertex = false;
    bool accPT = false;
    bool accETA = false;

    if (debug) std::cout << ">> EVENT " << i << " <<"  <<endl;

    if (protonLorentzVector->size() == 1){
      if (debug) cout << "--> proton(1) pZ [GeV]: " << protonLorentzVector->at(0).pz() << endl;
      perc = protonLorentzVector->at(0).pz()/EBeam;
    }

    if (protonLorentzVector->size() > 1){
      if(debug){
	cout << "--> proton(1) pZ [GeV]: " << protonLorentzVector->at(0).pz() << endl;
	cout << "--> proton(2) pZ [GeV]: " << protonLorentzVector->at(1).pz() << endl;
      }
      perc = protonLorentzVector->at(0).pz()/EBeam;
    }

    if (protonLorentzVector->size() > 0 && protonLorentzVector->at(0).pz() > 0. && fabs(perc) > 0.75){
      protonPlus = true;
      xi_proton_plus =  ( 1. - (protonLorentzVector->at(0).pz()/EBeam) );

      TLorentzVector vec_pi(0.,0.,EBeam,EBeam);
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
      hVectorVertex.at(0)->Fill(nVertex);

    }

    if (protonLorentzVector->size() > 0 && protonLorentzVector->at(0).pz() < 0. && fabs(perc) > 0.75){
      protonMinus = true;
      xi_proton_minus =  ( 1. + (protonLorentzVector->at(0).pz()/EBeam) );

      TLorentzVector vec_pi(0.,0.,-EBeam,EBeam);
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
      hVectorVertex.at(0)->Fill(nVertex);

    } 

    if (debug) cout << "-- END --\n" << endl;

    if (nVertex == 1) SingleVertex = true;
    

   
    if (JpsiMassDiMuon > 3.05 && JpsiMassDiMuon < 3.15) JpsiMassMu = true;

    if(SingleVertex && JpsiMassMu ) JpsiFillMu = true;

    // Isolation Two Leading Muons
    if(JpsiFillMu){
      if (LeadingMuonsIsolation->at(0) < 3 && LeadingMuonsIsolation->at(1) < 3 ) {
	isolation = true;
	candSel = true;
      }
    }



    // Fill Acceptance Histogram
    hVectorAccept.at(0)->Fill(accept);
  
    if (JpsiFillMu && LeadingMuonsP4->size() > 1){
       if(fabs(LeadingMuonsP4->at(0).eta()) < 2.45 && fabs(LeadingMuonsP4->at(1).eta()) > 2.45) accETA = true;
       if(LeadingMuonsP4->at(0).pt() > 1.0 && LeadingMuonsP4->at(1).pt() > 1.0) accPT = true;
      
        acceptJpsiMu = true;
    }
    
    // Jpsi Meson, Proton Minus
    if(protonMinus && accept > 0 && isolation && candSel && (acceptJpsiMu)){

      AcceptJpsi+=accept;
      ++selectedJpsi;

      hVectorVertex.at(1)->Fill(nVertex);
      hVectorAccept.at(1)->Fill(accept);
      hVectorProtonEta.at(1)->Fill(protonLorentzVector->at(0).eta(),accept);
      hVectorProtonPz.at(1)->Fill(protonLorentzVector->at(0).pz(),accept);
      hVectorProtonEnergy.at(1)->Fill(protonLorentzVector->at(0).energy(),accept);
      if(xi_proton_minus > -999.) hVectorProtonXi.at(1)->Fill(xi_proton_minus,accept);
      if(xi_proton_plus > -999.) hVectorProtonXi.at(1)->Fill(xi_proton_plus,accept);
      if(t_proton_minus > -999.) hVectorProtonT.at(1)->Fill(fabs(t_proton_minus),accept);
      if(t_proton_plus > -999.) hVectorProtonT.at(1)->Fill(fabs(t_proton_plus),accept);

      if(acceptJpsiMu){
	hVectorMesonM.at(1)->Fill(JpsiMassDiMuon,accept);
	hVectorMesonEta.at(1)->Fill(JpsiEtaDiMuon,accept);
	hVectorLeptonsEta.at(1)->Fill(LeadingMuonsP4->at(0).eta(),accept);
	hVectorLeptonsEta.at(1)->Fill(LeadingMuonsP4->at(1).eta(),accept);
	hVectorLeptonsPhi.at(1)->Fill(LeadingMuonsP4->at(0).phi(),accept);
	hVectorLeptonsPhi.at(1)->Fill(LeadingMuonsP4->at(1).phi(),accept);
	hVectorLeptonsPt.at(1)->Fill(LeadingMuonsP4->at(0).pt(),accept);
	hVectorLeptonsPt.at(1)->Fill(LeadingMuonsP4->at(1).pt(),accept);
	hVectorCrossSection.at(1)->Fill(JpsiMassDiMuon,accept);
      }
    }

    // Jpsi Meson, Proton Plus
    if(protonPlus && accept > 0 && isolation && candSel && (acceptJpsiMu) ){

      AcceptJpsi+=accept;
      ++selectedJpsi;

      hVectorVertex.at(1)->Fill(nVertex);
      hVectorAccept.at(1)->Fill(accept);
      hVectorProtonEta.at(1)->Fill(protonLorentzVector->at(0).eta(),accept);
      hVectorProtonPz.at(1)->Fill(protonLorentzVector->at(0).pz(),accept);
      hVectorProtonEnergy.at(1)->Fill(protonLorentzVector->at(0).energy(),accept);
      if(xi_proton_minus > -999.) hVectorProtonXi.at(1)->Fill(xi_proton_minus,accept);
      if(xi_proton_plus > -999.) hVectorProtonXi.at(1)->Fill(xi_proton_plus,accept);
      if(t_proton_minus > -999.) hVectorProtonT.at(1)->Fill(fabs(t_proton_minus),accept);
      if(t_proton_plus > -999.) hVectorProtonT.at(1)->Fill(fabs(t_proton_plus),accept);

      if(acceptJpsiMu){
	hVectorMesonM.at(1)->Fill(JpsiMassDiMuon,accept);
	hVectorMesonEta.at(1)->Fill(JpsiEtaDiMuon,accept);
	hVectorLeptonsEta.at(1)->Fill(LeadingMuonsP4->at(0).eta(),accept);
	hVectorLeptonsEta.at(1)->Fill(LeadingMuonsP4->at(1).eta(),accept);
	hVectorLeptonsPhi.at(1)->Fill(LeadingMuonsP4->at(0).phi(),accept);
	hVectorLeptonsPhi.at(1)->Fill(LeadingMuonsP4->at(1).phi(),accept);
	hVectorLeptonsPt.at(1)->Fill(LeadingMuonsP4->at(0).pt(),accept);
	hVectorLeptonsPt.at(1)->Fill(LeadingMuonsP4->at(1).pt(),accept);
	hVectorCrossSection.at(1)->Fill(JpsiMassDiMuon,accept);
      }
    }

  }

  hVectorCrossSection[1]->Scale(weight,"width");

  out->cd();
  hVectorVertex[0]->Write();
  hVectorVertex[1]->Write();
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
  hVectorMesonM[1]->Write();
  hVectorMesonEta[1]->Write();
  hVectorLeptonsEta[1]->Write();
  hVectorLeptonsPhi[1]->Write();
  hVectorLeptonsPt[1]->Write();
  hVectorCrossSection[1]->Write();
  hVectorProtonAcceptanceTPlus[0]->Write();
  hVectorProtonAcceptanceTMinus[0]->Write();
  hVectorProtonAcceptanceXiPlus[0]->Write();
  hVectorProtonAcceptanceXiMinus[0]->Write();
  hVectorAccept[0]->Write();
  hVectorAccept[1]->Write();

  out->Close();
  inf->Close();
  RPFileCMSMinus->Close();
  RPFileCMSPlus->Close();

  cout << "\nS U M M A R Y" << endl;
  cout << "\n--> Meson Jpsi: " << endl;
  cout << "Total Selected Events Jpsi (not normalized): " << selectedJpsi << endl;
  if(selectedJpsi > 0) cout << "Average Acceptance for Jpsi: " << AcceptJpsi/selectedJpsi << endl;
  cout << "Visible Cross Section for Jpsi, weighted: " << (XSmcJpsi*AcceptJpsi)/NEntries << "\n" << endl;

  // Saving Text File

  outstring << "\nI N P U T S" << endl;
  outstring << "Reading TREE: "<<NEntries<<" events"<<endl;
  outstring << "Input filename: " << inputfile << endl;
  outstring << "Output filename: " << outputfile << endl;
  outstring << "Integrated Luminosity: " << lumi << endl;
  outstring << "Cross Section MC Jpsi: " << XSmcJpsi << endl;
  outstring << "\nS U M M A R Y" << endl;
  outstring << "\n--> Meson Jpsi: " << endl;
  outstring << "Total Selected Events Jpsi (not normalized): " << selectedJpsi << endl;
  outstring << "Total Selected Events Jpsi (normalized): " << AcceptJpsi << endl;
  if(selectedJpsi > 0) outstring << "Average Acceptance for Jpsi: " << AcceptJpsi/selectedJpsi << endl;
  outstring << "Visible Cross Section for Jpsi, weighted: " << (XSmcJpsi*AcceptJpsi)/NEntries << "\n" << endl;
  outstring << "Visible Cross Section for Jpsi: " << (XSmcJpsi*selectedJpsi)/NEntries << endl;

}


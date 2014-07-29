/*

// Y E L L O W   R E P O R T   M A C R O 
// Author: D. Figueiredo and E. Melo

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

void CEPRPTotem(string inputfile, string outputfile,double XSmc, double lumi)
{

  //Parameters
  // Lumi = luminosity [pb-1]
  // XSmc = cross section MC

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
  TTree *tr = (TTree*)inf->Get("DijetsAnalyzer/Event");

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
  std::vector <ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *LeadingJetsP4;
  std::vector <ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *PFP4;

  int nTracks;
  int nVertex;
  double Mjj;
  double Mx;
  double Mpf;
  double RjjMx;
  double RjjMpf;
  double etamax;
  double etamin;
  double accept = 0.;
  int selected = 0;

  TBranch *bprotonLorentzVector;
  TBranch *bLeadingJetsP4;
  TBranch *bPFP4;

  TBranch *bnTracks;
  TBranch *bnVertex;
  TBranch *bMjj;
  TBranch *bMpf;

  tr->SetBranchAddress("ProtonsP4",&protonLorentzVector,&bprotonLorentzVector);
  tr->SetBranchAddress("JetsP4",&LeadingJetsP4,&bLeadingJetsP4);
  tr->SetBranchAddress("PFP4",&PFP4,&bPFP4);
  tr->SetBranchAddress("Mjj",&Mjj,&bMjj);
  tr->SetBranchAddress("Mpf",&Mpf,&bMpf);
  tr->SetBranchAddress("nVertex",&nVertex,&bnVertex);
  tr->SetBranchAddress("nTracks",&nTracks,&bnTracks);

  // Creating Histograms
  std::vector<TH1F*> hVectorProtonAcceptanceXiMinus;
  std::vector<TH1F*> hVectorProtonAcceptanceXiPlus;
  std::vector<TH1F*> hVectorProtonAcceptanceTMinus;
  std::vector<TH1F*> hVectorProtonAcceptanceTPlus;
  std::vector<TH1F*> hVectorProtonEta;
  std::vector<TH1F*> hVectorProtonEnergy;
  std::vector<TH1F*> hVectorProtonPz;
  std::vector<TH1F*> hVectorProtonXiminus;
  std::vector<TH1F*> hVectorProtonTminus;
  std::vector<TH1F*> hVectorProtonXiplus;
  std::vector<TH1F*> hVectorProtonTplus;
  std::vector<TH1F*> hVectorVertex;
  std::vector<TH1F*> hVectorAccept;
  std::vector<TH1F*> hVectorDijetsM;
  std::vector<TH1F*> hVectorJetsEta;
  std::vector<TH1F*> hVectorJetsPt;
  std::vector<TH1F*> hVectorMpf;
  std::vector<TH1F*> hVectorMx;
  std::vector<TH1F*> hVectorRjjMpf;
  std::vector<TH1F*> hVectorRjjMx;

  std::string step0 = "no_accept_RP";
  std::string step1 = "accept_RP_boson";

  std::vector <std::string> GroupHisto;

  GroupHisto.push_back(step0);
  GroupHisto.push_back(step1);

  for (std::vector<std::string>::size_type j=0; j<GroupHisto.size(); j++){
    char name[300];

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

    sprintf(name,"ProtonEta_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonEta = new TH1F(name,";#eta; N events",100,-15,15);
    hVectorProtonEta.push_back(hProtonEta);

    sprintf(name,"ProtonPz_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonPz = new TH1F(name,";p_{z} [GeV]; N events",14000,-7000.,7000.);
    hVectorProtonPz.push_back(hProtonPz);

    sprintf(name,"ProtonEnergy_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonEnergy = new TH1F(name,";Energy [GeV]; N events",7000,0.,7000.);
    hVectorProtonEnergy.push_back(hProtonEnergy);

    sprintf(name,"ProtonXiMinus_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonXiminus = new TH1F(name,";#xi^{-}; N events",50,0.,0.2);
    hVectorProtonXiminus.push_back(hProtonXiminus);

    sprintf(name,"ProtonTMinus_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonTminus = new TH1F(name,";|t|^{-}; N events",50,0.,1);
    hVectorProtonTminus.push_back(hProtonTminus);

    sprintf(name,"ProtonXiPlus_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonXiplus = new TH1F(name,";#xi^{+}; N events",50,0.,0.2);
    hVectorProtonXiplus.push_back(hProtonXiplus);

    sprintf(name,"ProtonTPlus_%s",GroupHisto.at(j).c_str());
    TH1F *hProtonTplus = new TH1F(name,";|t|^{+}; N events",50,0.,1);
    hVectorProtonTplus.push_back(hProtonTplus);

    sprintf(name,"Mjj_%s",GroupHisto.at(j).c_str());
    TH1F *hDijetsM = new TH1F(name,";Mjj; N events",1000,0.,1000.);
    hVectorDijetsM.push_back(hDijetsM);

    sprintf(name,"Mpf_%s",GroupHisto.at(j).c_str());
    TH1F *hMpf = new TH1F(name,";Mpf; N events",1000,0.,1000.);
    hVectorMpf.push_back(hMpf);

    sprintf(name,"Mx_%s",GroupHisto.at(j).c_str());
    TH1F *hMx = new TH1F(name,";Mx; N events",1000,0.,1000.);
    hVectorMx.push_back(hMx);

    sprintf(name,"RjjMpf_%s",GroupHisto.at(j).c_str());
    TH1F *hRjjMpf = new TH1F(name,";Rjj = Mjj/Mpf; N events",500,0.,5.);
    hVectorRjjMpf.push_back(hRjjMpf);

    sprintf(name,"RjjMx_%s",GroupHisto.at(j).c_str());
    TH1F *hRjjMx = new TH1F(name,";Rjj = Mjj/Mx; N events",500,0.,5.);
    hVectorRjjMx.push_back(hRjjMx);

    sprintf(name,"JetsEta_%s",GroupHisto.at(j).c_str());
    TH1F *hJetsEta = new TH1F(name,";#eta; N events",100,-6.,6.);
    hVectorJetsEta.push_back(hJetsEta);

    sprintf(name,"JetsPt_%s",GroupHisto.at(j).c_str());
    TH1F *hJetsPt = new TH1F(name,";p_{T} [GeV]; N events",1500,0.,1500.);
    hVectorJetsPt.push_back(hJetsPt);


  }

  //unsigned NEntries = tr->GetEntries();
  unsigned NEntries = 500;

  cout << "\nR U N N I N G" << endl;
  cout << "Reading TREE: "<<NEntries<<" events"<<endl;
  cout << "Input filename: " << inputfile << endl;
  cout << "Output filename: " << outputfile << endl;
  cout << "Integrated Luminosity: " << lumi << endl;
  cout << "Cross Section MC CEP: " << XSmc << endl;

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

    double progress = 10.0*i/(1.0*NEntries);
    Int_t k = TMath::FloorNint(progress);
    if (k > decade)
      cout<<10*k<<" %"<<endl;
    decade = k; 

    double perc1 = -999;
    double perc2 = -999;
    double xi_proton_plus = -999.;
    double xi_proton_minus = -999.;
    double t_proton_plus = -999.;
    double t_proton_minus = -999.;
    double acceptMinus = -999.;
    double acceptPlus = -999.;
    double Mx = -999.;

    bool accPT = false;
    bool accETA = false;
    bool SingleVertex = false;

    bool Hminus = false;
    bool Hplus = false;
    bool genSel = false;

    if (debug) std::cout << ">> EVENT " << i << " <<"  <<endl;

    if (protonLorentzVector->size() == 1){
      if (debug) cout << "--> proton(1) pZ [GeV]: " << protonLorentzVector->at(0).pz() << endl;
      perc1 = protonLorentzVector->at(0).pz()/EBeam;
    }

    if (protonLorentzVector->size() > 1){
      if(debug){
	cout << "--> proton(1) pZ [GeV]: " << protonLorentzVector->at(0).pz() << endl;
	cout << "--> proton(2) pZ [GeV]: " << protonLorentzVector->at(1).pz() << endl;
      }
      perc1 = protonLorentzVector->at(0).pz()/EBeam;
      perc2 = protonLorentzVector->at(1).pz()/EBeam;

      if(protonLorentzVector->at(0).pz() > 0. && protonLorentzVector->at(1).pz() < 0.) Hplus = true; 
      if(protonLorentzVector->at(0).pz() < 0. && protonLorentzVector->at(1).pz() > 0.) Hminus = true;

    }

    if (protonLorentzVector->size()>1 && fabs(perc1) > 0.75 && fabs(perc2) > 0.75) genSel = true;


    if(genSel){

      if(Hplus){

	xi_proton_plus =  ( 1. - (protonLorentzVector->at(0).pz()/EBeam) );
	TLorentzVector vec_pi_plus(0.,0.,EBeam,EBeam);
	TLorentzVector vec_pf_plus(protonLorentzVector->at(0).px(),protonLorentzVector->at(0).py(),protonLorentzVector->at(0).pz(),protonLorentzVector->at(0).energy());
	TLorentzVector vec_t_plus = (vec_pf_plus - vec_pi_plus);
	t_proton_plus = vec_t_plus.Mag2();
	acceptPlus = HistoRPCMSPlus->GetBinContent(HistoRPCMSPlus->GetXaxis()->FindBin(fabs(t_proton_plus)),HistoRPCMSPlus->GetYaxis()->FindBin(xi_proton_plus));

	xi_proton_minus =  ( 1. + (protonLorentzVector->at(1).pz()/EBeam) );
	TLorentzVector vec_pi_minus(0.,0.,-EBeam,EBeam);
	TLorentzVector vec_pf_minus(protonLorentzVector->at(1).px(),protonLorentzVector->at(1).py(),protonLorentzVector->at(1).pz(),protonLorentzVector->at(1).energy());
	TLorentzVector vec_t_minus = (vec_pf_minus - vec_pi_minus);
	t_proton_minus = vec_t_minus.Mag2();
	acceptMinus = HistoRPCMSMinus->GetBinContent(HistoRPCMSMinus->GetXaxis()->FindBin(fabs(t_proton_minus)),HistoRPCMSMinus->GetYaxis()->FindBin(xi_proton_minus));
	Mx = EBeam*TMath::Sqrt(xi_proton_minus*xi_proton_plus);

	if(debug){
	  cout << "\nProton(0)+ and Proton(1)-" << endl;
	  cout << "--> xi, plus: " << xi_proton_plus << endl;
	  cout << "-->  t, plus: " << t_proton_plus << endl;
	  cout << "--> RP Acceptance, plus: " << acceptPlus << endl;
	  cout << "--> xi, minus: " << xi_proton_minus << endl;
	  cout << "-->  t, minus: " << t_proton_minus << endl;
	  cout << "--> RP Acceptance, minus: " << acceptMinus << endl;
	  cout << "--> Mx: " << Mx << endl;
	}

      }

      if(Hminus){

	xi_proton_plus =  ( 1. - (protonLorentzVector->at(1).pz()/EBeam) );
	TLorentzVector vec_pi_plus(0.,0.,EBeam,EBeam);
	TLorentzVector vec_pf_plus(protonLorentzVector->at(1).px(),protonLorentzVector->at(1).py(),protonLorentzVector->at(1).pz(),protonLorentzVector->at(1).energy());
	TLorentzVector vec_t_plus = (vec_pf_plus - vec_pi_plus);
	t_proton_plus = vec_t_plus.Mag2();
	acceptPlus = HistoRPCMSPlus->GetBinContent(HistoRPCMSPlus->GetXaxis()->FindBin(fabs(t_proton_plus)),HistoRPCMSPlus->GetYaxis()->FindBin(xi_proton_plus));

	xi_proton_minus =  ( 1. + (protonLorentzVector->at(0).pz()/EBeam) );
	TLorentzVector vec_pi_minus(0.,0.,-EBeam,EBeam);
	TLorentzVector vec_pf_minus(protonLorentzVector->at(0).px(),protonLorentzVector->at(0).py(),protonLorentzVector->at(0).pz(),protonLorentzVector->at(0).energy());
	TLorentzVector vec_t_minus = (vec_pf_minus - vec_pi_minus);
	t_proton_minus = vec_t_minus.Mag2();
	acceptMinus = HistoRPCMSMinus->GetBinContent(HistoRPCMSMinus->GetXaxis()->FindBin(fabs(t_proton_minus)),HistoRPCMSMinus->GetYaxis()->FindBin(xi_proton_minus));
	Mx = EBeam*TMath::Sqrt(xi_proton_minus*xi_proton_plus);

	if(debug){
	  cout << "\nProton(0)- and Proton(1)+" << endl;
	  cout << "--> xi, plus: " << xi_proton_plus << endl;
	  cout << "-->  t, plus: " << t_proton_plus << endl;
	  cout << "--> RP Acceptance, plus: " << acceptPlus << endl;
	  cout << "--> xi, minus: " << xi_proton_minus << endl;
	  cout << "-->  t, minus: " << t_proton_minus << endl;
	  cout << "--> RP Acceptance, minus: " << acceptMinus << endl;
	  cout << "--> Mx: " << Mx << endl;
	}

      }

      if (Hminus || Hplus){
	hVectorVertex.at(0)->Fill(nVertex);
	hVectorAccept.at(0)->Fill((acceptMinus + acceptPlus)/2.);
	hVectorProtonEta.at(0)->Fill(protonLorentzVector->at(0).eta());
	hVectorProtonPz.at(0)->Fill(protonLorentzVector->at(0).pz());
	hVectorProtonEnergy.at(0)->Fill(protonLorentzVector->at(0).energy());
	hVectorProtonEta.at(0)->Fill(protonLorentzVector->at(1).eta());
	hVectorProtonPz.at(0)->Fill(protonLorentzVector->at(1).pz());
	hVectorProtonEnergy.at(0)->Fill(protonLorentzVector->at(1).energy());
	hVectorProtonXiminus.at(0)->Fill(xi_proton_minus);
	hVectorProtonXiplus.at(0)->Fill(xi_proton_plus);
	hVectorProtonTminus.at(0)->Fill(fabs(t_proton_minus));
	hVectorProtonTplus.at(0)->Fill(fabs(t_proton_plus));
	if(LeadingJetsP4->size()>1){
	  hVectorDijetsM.at(0)->Fill(Mjj);
	  hVectorJetsPt.at(0)->Fill(LeadingJetsP4->at(0).pt());
	  hVectorJetsPt.at(0)->Fill(LeadingJetsP4->at(1).pt());
	  hVectorJetsEta.at(0)->Fill(LeadingJetsP4->at(0).eta());
	  hVectorJetsEta.at(0)->Fill(LeadingJetsP4->at(1).eta());
	  hVectorMpf.at(0)->Fill(Mpf);
	  hVectorMx.at(0)->Fill(Mx);
	  hVectorRjjMpf.at(0)->Fill(Mjj/Mpf);
	  hVectorRjjMx.at(0)->Fill(Mjj/Mx);
	}
      }

    }

    if (debug) cout << "-- END --\n" << endl;

    if(LeadingJetsP4->size()>1){
      if(fabs(LeadingJetsP4->at(0).eta()) < 3. && fabs(LeadingJetsP4->at(1).eta()) < 3. ) accETA = true;
      if(LeadingJetsP4->at(0).pt() > 30. && LeadingJetsP4->at(1).pt() > 30. ) accPT = true;
    }

    if (nVertex == 1) SingleVertex = true;

    if(SingleVertex && accPT && accETA && genSel && (Hplus || Hminus)) {
      //accept += (acceptMinus + acceptPlus)/2.;
      accept += acceptMinus*acceptPlus;
      ++selected;

      hVectorVertex.at(1)->Fill(nVertex,acceptMinus*acceptPlus);
      hVectorAccept.at(1)->Fill((acceptMinus + acceptPlus)/2.);
      hVectorProtonEta.at(1)->Fill(protonLorentzVector->at(0).eta(),acceptMinus*acceptPlus);
      hVectorProtonPz.at(1)->Fill(protonLorentzVector->at(0).pz(),acceptMinus*acceptPlus);
      hVectorProtonEnergy.at(1)->Fill(protonLorentzVector->at(0).energy(),acceptMinus*acceptPlus);
      hVectorProtonEta.at(1)->Fill(protonLorentzVector->at(1).eta(),acceptMinus*acceptPlus);
      hVectorProtonPz.at(1)->Fill(protonLorentzVector->at(1).pz(),acceptMinus*acceptPlus);
      hVectorProtonEnergy.at(1)->Fill(protonLorentzVector->at(1).energy(),acceptMinus*acceptPlus);
      hVectorProtonXiminus.at(1)->Fill(xi_proton_minus,acceptMinus*acceptPlus);
      hVectorProtonXiplus.at(1)->Fill(xi_proton_plus,acceptMinus*acceptPlus);
      hVectorProtonTminus.at(1)->Fill(fabs(t_proton_minus),acceptMinus*acceptPlus);
      hVectorProtonTplus.at(1)->Fill(fabs(t_proton_plus),acceptMinus*acceptPlus);
      hVectorJetsPt.at(1)->Fill(LeadingJetsP4->at(0).pt(),acceptMinus*acceptPlus);
      hVectorJetsPt.at(1)->Fill(LeadingJetsP4->at(1).pt(),acceptMinus*acceptPlus);
      hVectorJetsEta.at(1)->Fill(LeadingJetsP4->at(0).eta(),acceptMinus*acceptPlus);
      hVectorJetsEta.at(1)->Fill(LeadingJetsP4->at(1).eta(),acceptMinus*acceptPlus);
      hVectorDijetsM.at(1)->Fill(Mjj,acceptMinus*acceptPlus);
      hVectorMpf.at(1)->Fill(Mpf,acceptMinus*acceptPlus);
      hVectorMx.at(1)->Fill(Mx,acceptMinus*acceptPlus);
      hVectorRjjMpf.at(1)->Fill(Mjj/Mpf,acceptMinus*acceptPlus);
      hVectorRjjMx.at(1)->Fill(Mjj/Mx,acceptMinus*acceptPlus);
    }


  }// end loop events

  out->cd();

  hVectorProtonAcceptanceTPlus[0]->Write();
  hVectorProtonAcceptanceTMinus[0]->Write();
  hVectorProtonAcceptanceXiPlus[0]->Write();
  hVectorProtonAcceptanceXiMinus[0]->Write();

  hVectorVertex[0]->Write();
  hVectorAccept[0]->Write();
  hVectorProtonEta[0]->Write();
  hVectorProtonPz[0]->Write();
  hVectorProtonEnergy[0]->Write();
  hVectorProtonXiminus[0]->Write();
  hVectorProtonTminus[0]->Write();
  hVectorProtonXiplus[0]->Write();
  hVectorProtonTplus[0]->Write();
  hVectorJetsPt[0]->Write();
  hVectorJetsEta[0]->Write();
  hVectorDijetsM[0]->Write();
  hVectorMpf[0]->Write();
  hVectorMx[0]->Write();
  hVectorRjjMpf[0]->Write();
  hVectorRjjMx[0]->Write();

  hVectorVertex[1]->Write();
  hVectorAccept[1]->Write();
  hVectorProtonEta[1]->Write();
  hVectorProtonPz[1]->Write();
  hVectorProtonEnergy[1]->Write();
  hVectorProtonXiminus[1]->Write();
  hVectorProtonTminus[1]->Write();
  hVectorProtonXiplus[1]->Write();
  hVectorProtonTplus[1]->Write();
  hVectorJetsPt[1]->Write();
  hVectorJetsEta[1]->Write();
  hVectorDijetsM[1]->Write();
  hVectorMpf[1]->Write();
  hVectorMx[1]->Write();
  hVectorRjjMpf[1]->Write();
  hVectorRjjMx[1]->Write();

  out->Close();
  inf->Close();
  RPFileCMSMinus->Close();
  RPFileCMSPlus->Close();

  cout << "\nS U M M A R Y" << endl;
  cout << "--> CEP: " << endl;
  cout << "Total Selected CEP events      : " << selected << endl;
  cout << "Normalized Selected CEP events : " << accept << endl;
  cout << "Visible Cross Section for CEP  : " << (XSmc*accept)/NEntries << "\n" << endl;

  // Saving Text File
  outstring << "\nI N P U T S" << endl;
  outstring << "Reading TREE: "<<NEntries<<" events"<<endl;
  outstring << "Input filename: " << inputfile << endl;
  outstring << "Output filename: " << outputfile << endl;
  outstring << "Integrated Luminosity: " << lumi << endl;
  outstring << "Cross Section CEP MC: " << XSmc << endl;
  outstring << "\nS U M M A R Y" << endl;
  outstring << "--> CEP: " << endl;
  outstring << "Total Selected CEP events      : " << selected << endl;
  outstring << "Normalized Selected CEP events : " << accept << endl;
  outstring << "Visible Cross Section for CEP  : " << (XSmc*accept)/NEntries << "\n" << endl;

}


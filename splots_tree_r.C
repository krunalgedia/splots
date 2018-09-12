#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include <cstring>
#include <sstream>
#include <stdlib.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include "TColor.h"
#include "TAxis.h"
#include "TColor.h"
#include "TAxis.h"
#include "TLorentzVector.h"
#include <TMath.h>
#include <TLegend.h>

using namespace std;

void splots_tree_r()
{
    TFile* data = new TFile("r_TMVA_output.root");
    TTree* datatree = (TTree*)(data->Get("tree"));

    vector<double>* im_nl_j = 0;
    vector<double>* iBDT_response = 0;
    vector<double>* ijet_pt = 0;
    vector<double>* ijet_indices = 0;
    double injet;

    UInt_t          irun;
    UInt_t          ilumi;
    ULong64_t       ievt;
    Float_t         iJet_leptonDeltaR[23];   //[nJet]
    Float_t         iJet_pt[23];   //[nJet]
    Float_t         iJet_eta[23];   //[nJet]
    Float_t         iJet_phi[23];   //[nJet]
    Float_t         iJet_mass[23];   //[nJet]
    Float_t         iJet_energyRing_dR0_em[23];   //[nJet]
    Float_t         iJet_energyRing_dR1_em[23];   //[nJet]
    Float_t         iJet_energyRing_dR2_em[23];   //[nJet]
    Float_t         iJet_energyRing_dR3_em[23];   //[nJet]
    Float_t         iJet_energyRing_dR4_em[23];   //[nJet]
    Float_t         iJet_energyRing_dR0_mu[23];   //[nJet]
    Float_t         iJet_energyRing_dR1_mu[23];   //[nJet]
    Float_t         iJet_energyRing_dR2_mu[23];   //[nJet]
    Float_t         iJet_energyRing_dR3_mu[23];   //[nJet]
    Float_t         iJet_energyRing_dR4_mu[23];   //[nJet]
    Float_t         iJet_energyRing_dR0_ch[23];   //[nJet]
    Float_t         iJet_energyRing_dR1_ch[23];   //[nJet]
    Float_t         iJet_energyRing_dR2_ch[23];   //[nJet]
    Float_t         iJet_energyRing_dR3_ch[23];   //[nJet]
    Float_t         iJet_energyRing_dR4_ch[23];   //[nJet]
    Float_t         iJet_energyRing_dR0_neut[23];   //[nJet]
    Float_t         iJet_energyRing_dR1_neut[23];   //[nJet]
    Float_t         iJet_energyRing_dR2_neut[23];   //[nJet]
    Float_t         iJet_energyRing_dR3_neut[23];   //[nJet]
    Float_t         iJet_energyRing_dR4_neut[23];   //[nJet]
    Int_t           iJet_numDaughters_pt03[23];   //[nJet]
    Float_t         iJet_rawEnergy[23];   //[nJet]
    Int_t           iJet_numberOfDaughters[23];   //[nJet]
    Float_t         irho;
    Float_t         iJet_rawPt[23];   //[nJet]
    Float_t         iJet_chHEF[23];   //[nJet]
    Float_t         iJet_neHEF[23];   //[nJet]
    Float_t         iJet_chEmEF[23];   //[nJet]
    Float_t         iJet_neEmEF[23];   //[nJet]
    Float_t         iJet_muEF[23];   //[nJet]
    Float_t         iJet_leadTrackPt[23];   //[nJet]
    Float_t         iJet_leptonPt[23];   //[nJet]
    Float_t         iJet_leptonPtRel[23];   //[nJet]
    Float_t         iJet_leptonPdgId[23];   //[nJet]
    Float_t         iJet_leptonPtRelInv[23];   //[nJet]
    Float_t         iJet_vtxMass[23];   //[nJet]
    Float_t         iJet_vtxNtracks[23];   //[nJet]
    Float_t         iJet_vtxPt[23];   //[nJet]
    Float_t         iJet_vtx3DSig[23];   //[nJet]
    Float_t         iJet_vtx3DVal[23];   //[nJet]
    Float_t         iJet_ptd[23];   //[nJet]

    TBranch        *ib_run;   //!
    TBranch        *ib_lumi;   //!
    TBranch        *ib_evt;   //!
    TBranch        *ib_Jet_leptonDeltaR;   //!
    TBranch        *ib_Jet_pt;   //!
    TBranch        *ib_Jet_eta;   //!
    TBranch        *ib_Jet_phi;   //!
    TBranch        *ib_Jet_mass;   //!
    TBranch        *ib_Jet_energyRing_dR0_em;   //!
    TBranch        *ib_Jet_energyRing_dR1_em;   //!
    TBranch        *ib_Jet_energyRing_dR2_em;   //!
    TBranch        *ib_Jet_energyRing_dR3_em;   //!
    TBranch        *ib_Jet_energyRing_dR4_em;   //!
    TBranch        *ib_Jet_energyRing_dR0_mu;   //!
    TBranch        *ib_Jet_energyRing_dR1_mu;   //!
    TBranch        *ib_Jet_energyRing_dR2_mu;   //!
    TBranch        *ib_Jet_energyRing_dR3_mu;   //!
    TBranch        *ib_Jet_energyRing_dR4_mu;   //!
    TBranch        *ib_Jet_energyRing_dR0_ch;   //!
    TBranch        *ib_Jet_energyRing_dR1_ch;   //!
    TBranch        *ib_Jet_energyRing_dR2_ch;   //!
    TBranch        *ib_Jet_energyRing_dR3_ch;   //!
    TBranch        *ib_Jet_energyRing_dR4_ch;   //!
    TBranch        *ib_Jet_energyRing_dR0_neut;   //!
    TBranch        *ib_Jet_energyRing_dR1_neut;   //!
    TBranch        *ib_Jet_energyRing_dR2_neut;   //!
    TBranch        *ib_Jet_energyRing_dR3_neut;   //!
    TBranch        *ib_Jet_energyRing_dR4_neut;   //!
    TBranch        *ib_Jet_numDaughters_pt03;   //!
    TBranch        *ib_Jet_numberOfDaughters;   //!
    TBranch        *ib_rho;   //!
    TBranch        *ib_Jet_rawEnergy;   //!
    TBranch        *ib_Jet_rawPt;   //!
    TBranch        *ib_Jet_chHEF;   //!
    TBranch        *ib_Jet_neHEF;   //!
    TBranch        *ib_Jet_chEmEF;   //!
    TBranch        *ib_Jet_neEmEF;   //!
    TBranch        *ib_Jet_muEF;   //!
    TBranch        *ib_Jet_leadTrackPt;   //!
    TBranch        *ib_Jet_leptonPt;   //!
    TBranch        *ib_Jet_leptonPtRel;   //!
    TBranch        *ib_Jet_leptonPtRelInv;   //!
    TBranch        *ib_Jet_leptonPdgId;   //!
    TBranch        *ib_Jet_vtxMass;   //!
    TBranch        *ib_Jet_vtxNtracks;   //!
    TBranch        *ib_Jet_vtxPt;   //!
    TBranch        *ib_Jet_vtx3DSig;   //!
    TBranch        *ib_Jet_vtx3DVal;   //!
    TBranch        *ib_Jet_ptd;   //!


    datatree->SetBranchAddress("r_BDT_response",&iBDT_response);
    datatree->SetBranchAddress("jet_pt",&ijet_pt);
    datatree->SetBranchAddress("njet",&injet);
    datatree->SetBranchAddress("jet_indices",&ijet_indices);
    datatree->SetBranchAddress("m_nl_j",&im_nl_j );

    datatree->SetBranchAddress("run", &irun, &ib_run);
    datatree->SetBranchAddress("lumi", &ilumi, &ib_lumi);
    datatree->SetBranchAddress("evt", &ievt, &ib_evt);
    datatree->SetBranchAddress("rho", &irho, &ib_rho);
    datatree->SetBranchAddress("Jet_leptonDeltaR", iJet_leptonDeltaR, &ib_Jet_leptonDeltaR);
    datatree->SetBranchAddress("Jet_pt", iJet_pt, &ib_Jet_pt);
    datatree->SetBranchAddress("Jet_eta", iJet_eta, &ib_Jet_eta);
    datatree->SetBranchAddress("Jet_phi", iJet_phi, &ib_Jet_phi);
    datatree->SetBranchAddress("Jet_mass", iJet_mass, &ib_Jet_mass);
    datatree->SetBranchAddress("Jet_energyRing_dR0_em", iJet_energyRing_dR0_em, &ib_Jet_energyRing_dR0_em);
    datatree->SetBranchAddress("Jet_energyRing_dR1_em", iJet_energyRing_dR1_em, &ib_Jet_energyRing_dR1_em);
    datatree->SetBranchAddress("Jet_energyRing_dR2_em", iJet_energyRing_dR2_em, &ib_Jet_energyRing_dR2_em);
    datatree->SetBranchAddress("Jet_energyRing_dR3_em", iJet_energyRing_dR3_em, &ib_Jet_energyRing_dR3_em);
    datatree->SetBranchAddress("Jet_energyRing_dR4_em", iJet_energyRing_dR4_em, &ib_Jet_energyRing_dR4_em);
    datatree->SetBranchAddress("Jet_energyRing_dR0_mu", iJet_energyRing_dR0_mu, &ib_Jet_energyRing_dR0_mu);
    datatree->SetBranchAddress("Jet_energyRing_dR1_mu", iJet_energyRing_dR1_mu, &ib_Jet_energyRing_dR1_mu);
    datatree->SetBranchAddress("Jet_energyRing_dR2_mu", iJet_energyRing_dR2_mu, &ib_Jet_energyRing_dR2_mu);
    datatree->SetBranchAddress("Jet_energyRing_dR3_mu", iJet_energyRing_dR3_mu, &ib_Jet_energyRing_dR3_mu);
    datatree->SetBranchAddress("Jet_energyRing_dR4_mu", iJet_energyRing_dR4_mu, &ib_Jet_energyRing_dR4_mu);
    datatree->SetBranchAddress("Jet_energyRing_dR0_ch", iJet_energyRing_dR0_ch, &ib_Jet_energyRing_dR0_ch);
    datatree->SetBranchAddress("Jet_energyRing_dR1_ch", iJet_energyRing_dR1_ch, &ib_Jet_energyRing_dR1_ch);
    datatree->SetBranchAddress("Jet_energyRing_dR2_ch", iJet_energyRing_dR2_ch, &ib_Jet_energyRing_dR2_ch);
    datatree->SetBranchAddress("Jet_energyRing_dR3_ch", iJet_energyRing_dR3_ch, &ib_Jet_energyRing_dR3_ch);
    datatree->SetBranchAddress("Jet_energyRing_dR4_ch", iJet_energyRing_dR4_ch, &ib_Jet_energyRing_dR4_ch);
    datatree->SetBranchAddress("Jet_energyRing_dR0_neut", iJet_energyRing_dR0_neut, &ib_Jet_energyRing_dR0_neut);
    datatree->SetBranchAddress("Jet_energyRing_dR1_neut", iJet_energyRing_dR1_neut, &ib_Jet_energyRing_dR1_neut);
    datatree->SetBranchAddress("Jet_energyRing_dR2_neut", iJet_energyRing_dR2_neut, &ib_Jet_energyRing_dR2_neut);
    datatree->SetBranchAddress("Jet_energyRing_dR3_neut", iJet_energyRing_dR3_neut, &ib_Jet_energyRing_dR3_neut);
    datatree->SetBranchAddress("Jet_energyRing_dR4_neut", iJet_energyRing_dR4_neut, &ib_Jet_energyRing_dR4_neut);
    datatree->SetBranchAddress("Jet_rawEnergy", iJet_rawEnergy, &ib_Jet_rawEnergy);
    datatree->SetBranchAddress("Jet_numDaughters_pt03", iJet_numDaughters_pt03, &ib_Jet_numDaughters_pt03);
    datatree->SetBranchAddress("Jet_numberOfDaughters", iJet_numberOfDaughters, &ib_Jet_numberOfDaughters);
    datatree->SetBranchAddress("Jet_rawPt", iJet_rawPt, &ib_Jet_rawPt);
    datatree->SetBranchAddress("Jet_chHEF", iJet_chHEF, &ib_Jet_chHEF);
    datatree->SetBranchAddress("Jet_neHEF", iJet_neHEF, &ib_Jet_neHEF);
    datatree->SetBranchAddress("Jet_chEmEF", iJet_chEmEF, &ib_Jet_chEmEF);
    datatree->SetBranchAddress("Jet_neEmEF", iJet_neEmEF, &ib_Jet_neEmEF);
    datatree->SetBranchAddress("Jet_muEF", iJet_muEF, &ib_Jet_muEF);
    datatree->SetBranchAddress("Jet_leadTrackPt", iJet_leadTrackPt, &ib_Jet_leadTrackPt);
    datatree->SetBranchAddress("Jet_leptonPdgId", iJet_leptonPdgId, &ib_Jet_leptonPdgId);
    datatree->SetBranchAddress("Jet_leptonPt", iJet_leptonPt, &ib_Jet_leptonPt);
    datatree->SetBranchAddress("Jet_leptonPtRel", iJet_leptonPtRel, &ib_Jet_leptonPtRel);
    datatree->SetBranchAddress("Jet_leptonPtRelInv", iJet_leptonPtRelInv, &ib_Jet_leptonPtRelInv);
    datatree->SetBranchAddress("Jet_leptonDeltaR", iJet_leptonDeltaR, &ib_Jet_leptonDeltaR);
    datatree->SetBranchAddress("Jet_vtxMass", iJet_vtxMass, &ib_Jet_vtxMass);
    datatree->SetBranchAddress("Jet_vtxNtracks", iJet_vtxNtracks, &ib_Jet_vtxNtracks);
    datatree->SetBranchAddress("Jet_vtxPt", iJet_vtxPt, &ib_Jet_vtxPt);
    datatree->SetBranchAddress("Jet_vtx3DSig", iJet_vtx3DSig, &ib_Jet_vtx3DSig);
    datatree->SetBranchAddress("Jet_vtx3DVal", iJet_vtx3DVal, &ib_Jet_vtx3DVal);
    datatree->SetBranchAddress("Jet_ptd", iJet_ptd, &ib_Jet_ptd);


    TFile *outputFile = TFile::Open("r_TMVA_output_splots.root", "RECREATE" );
    TTree *splots_tree = new TTree("tree","tree");

    double BDT_response;
    double BDT_response_20;
    double BDT_response_30;
    double BDT_response_40;
    double BDT_response_60;
    double BDT_response_80;
    double BDT_response_100;
    double BDT_response_120;

    double jet_pt;
    int jet_indices;
    double m_nl_j;
//    TBranch *nBDT_response = splots_tree->Branch("r_BDT_response",&BDT_response);
    double         rho;
    double         rho_copy;
    double         run;
    double         lumi;
    double         evt;
    double         run_copy;
    double         lumi_copy;
    double         evt_copy;
      
    double         Jet_leptonDeltaR;   //[nJet]
    double         Jet_pt;   //[nJet]
    double         Jet_eta;   //[nJet]
    double         Jet_phi;   //[nJet]
    double         Jet_mass;   //[nJet]
    double         Jet_energyRing_dR0_em;   //[nJet]
    double         Jet_energyRing_dR1_em;   //[nJet]
    double         Jet_energyRing_dR2_em;   //[nJet]
    double         Jet_energyRing_dR3_em;   //[nJet]
    double         Jet_energyRing_dR4_em;   //[nJet]
    double         Jet_energyRing_dR0_mu;   //[nJet]
    double         Jet_energyRing_dR1_mu;   //[nJet]
    double         Jet_energyRing_dR2_mu;   //[nJet]
    double         Jet_energyRing_dR3_mu;   //[nJet]
    double         Jet_energyRing_dR4_mu;   //[nJet]
    double         Jet_energyRing_dR0_ch;   //[nJet]
    double         Jet_energyRing_dR1_ch;   //[nJet]
    double         Jet_energyRing_dR2_ch;   //[nJet]
    double         Jet_energyRing_dR3_ch;   //[nJet]
    double         Jet_energyRing_dR4_ch;   //[nJet]
    double         Jet_energyRing_dR0_neut;   //[nJet]
    double         Jet_energyRing_dR1_neut;   //[nJet]
    double         Jet_energyRing_dR2_neut;   //[nJet]
    double         Jet_energyRing_dR3_neut;   //[nJet]
    double         Jet_energyRing_dR4_neut;   //[nJet]
    int            Jet_numDaughters_pt03;   //[nJet]
    double         Jet_rawEnergy;   //[nJet]
    int            Jet_numberOfDaughters;   //[nJet]
    double         Jet_rawPt;   //[nJet]
    double         Jet_chHEF;   //[nJet]
    double         Jet_neHEF;   //[nJet]
    double         Jet_chEmEF;   //[nJet]
    double         Jet_neEmEF;   //[nJet]
    double         Jet_muEF;   //[nJet]
    double         Jet_leadTrackPt;   //[nJet]
    double         Jet_leptonPt;   //[nJet]
    double         Jet_leptonPtRel;   //[nJet]
    double         Jet_leptonPdgId;   //[nJet]
    double         Jet_leptonPtRelInv;   //[nJet]
    double         Jet_vtxMass;   //[nJet]
    double         Jet_vtxNtracks;   //[nJet]
    double         Jet_vtxPt;   //[nJet]
    double         Jet_vtx3DSig;   //[nJet]
    double         Jet_vtx3DVal;   //[nJet]
    double         Jet_ptd;   //[nJet]


    TBranch *nrho_copy = splots_tree->Branch("rho_copy",&rho_copy);
    TBranch *nrun = splots_tree->Branch("run",&run);
    TBranch *nlumi = splots_tree->Branch("lumi",&lumi);
    TBranch *nevt = splots_tree->Branch("evt",&evt);
    TBranch *nrun_copy = splots_tree->Branch("run_copy",&run_copy);
    TBranch *nlumi_copy = splots_tree->Branch("lumi_copy",&lumi_copy);
    TBranch *nevt_copy = splots_tree->Branch("evt_copy",&evt_copy);

    TBranch *nBDT_response = splots_tree->Branch("BDT_response",&BDT_response);
    TBranch *nBDT_response_20 = splots_tree->Branch("r_BDT_response_20",&BDT_response_20);
    TBranch *nBDT_response_30 = splots_tree->Branch("r_BDT_response_30",&BDT_response_30);
    TBranch *nBDT_response_40 = splots_tree->Branch("r_BDT_response_40",&BDT_response_40);
    TBranch *nBDT_response_60 = splots_tree->Branch("r_BDT_response_60",&BDT_response_60);
    TBranch *nBDT_response_80 = splots_tree->Branch("r_BDT_response_80",&BDT_response_80);
    TBranch *nBDT_response_100 = splots_tree->Branch("r_BDT_response_100",&BDT_response_100);
    TBranch *nBDT_response_120 = splots_tree->Branch("r_BDT_response_120",&BDT_response_120);
    TBranch *njet_pt = splots_tree->Branch("jet_pt",&jet_pt);
    TBranch *njet_indices = splots_tree->Branch("jet_indices",&jet_indices);
    TBranch *nm_nl_j = splots_tree->Branch("m_nl_j",&m_nl_j);
    TBranch *nJet_leptonDeltaR  = splots_tree->Branch("Jet_leptonDeltaR",&Jet_leptonDeltaR);
    TBranch *nJet_pt = splots_tree->Branch("Jet_pt",&Jet_pt);
    TBranch *nJet_eta = splots_tree->Branch("Jet_eta",&Jet_eta);
    TBranch *nJet_phi = splots_tree->Branch("Jet_phi",&Jet_phi);
    TBranch *nJet_mass = splots_tree->Branch("Jet_mass",&Jet_mass);
    TBranch *nJet_energyRing_dR0_em = splots_tree->Branch("Jet_energyRing_dR0_em",&Jet_energyRing_dR0_em);
    TBranch *nJet_energyRing_dR1_em = splots_tree->Branch("Jet_energyRing_dR1_em",&Jet_energyRing_dR1_em);
    TBranch *nJet_energyRing_dR2_em = splots_tree->Branch("Jet_energyRing_dR2_em",&Jet_energyRing_dR2_em);
    TBranch *nJet_energyRing_dR3_em = splots_tree->Branch("Jet_energyRing_dR3_em",&Jet_energyRing_dR3_em);
    TBranch *nJet_energyRing_dR4_em = splots_tree->Branch("Jet_energyRing_dR4_em",&Jet_energyRing_dR4_em);
    TBranch *nJet_energyRing_d0R_mu = splots_tree->Branch("Jet_energyRing_dR0_mu",&Jet_energyRing_dR0_mu);
    TBranch *nJet_energyRing_dR1_mu = splots_tree->Branch("Jet_energyRing_dR1_mu",&Jet_energyRing_dR1_mu);
    TBranch *nJet_energyRing_dR2_mu = splots_tree->Branch("Jet_energyRing_dR2_mu",&Jet_energyRing_dR2_mu);
    TBranch *nJet_energyRing_dR3_mu = splots_tree->Branch("Jet_energyRing_dR3_mu",&Jet_energyRing_dR3_mu);
    TBranch *nJet_energyRing_dR4_mu = splots_tree->Branch("Jet_energyRing_dR4_mu",&Jet_energyRing_dR4_mu);
    TBranch *nJet_energyRing_d0R_ch = splots_tree->Branch("Jet_energyRing_dR0_ch",&Jet_energyRing_dR0_ch);
    TBranch *nJet_energyRing_dR1_ch = splots_tree->Branch("Jet_energyRing_dR1_ch",&Jet_energyRing_dR1_ch);
    TBranch *nJet_energyRing_dR2_ch = splots_tree->Branch("Jet_energyRing_dR2_ch",&Jet_energyRing_dR2_ch);
    TBranch *nJet_energyRing_dR3_ch = splots_tree->Branch("Jet_energyRing_dR3_ch",&Jet_energyRing_dR3_ch);
    TBranch *nJet_energyRing_dR4_ch = splots_tree->Branch("Jet_energyRing_dR4_ch",&Jet_energyRing_dR4_ch);
    TBranch *nJet_energyRing_d0R_neut = splots_tree->Branch("Jet_energyRing_dR0_neut",&Jet_energyRing_dR0_neut);
    TBranch *nJet_energyRing_dR1_neut = splots_tree->Branch("Jet_energyRing_dR1_neut",&Jet_energyRing_dR1_neut);
    TBranch *nJet_energyRing_dR2_neut = splots_tree->Branch("Jet_energyRing_dR2_neut",&Jet_energyRing_dR2_neut);
    TBranch *nJet_energyRing_dR3_neut = splots_tree->Branch("Jet_energyRing_dR3_neut",&Jet_energyRing_dR3_neut);
    TBranch *nJet_energyRing_dR4_neut = splots_tree->Branch("Jet_energyRing_dR4_neut",&Jet_energyRing_dR4_neut);
    TBranch *nJet_numDaughters_pt03 = splots_tree->Branch("Jet_numDaughters_pt03",&Jet_numDaughters_pt03);
    TBranch *nJet_rawEnergy = splots_tree->Branch("Jet_rawEnergy",&Jet_rawEnergy);
    TBranch *nJet_numberOfDaughters = splots_tree->Branch("Jet_numberOfDaughters",&Jet_numberOfDaughters);
    TBranch *nrho = splots_tree->Branch("rho",&rho);
    TBranch *nJet_rawPt = splots_tree->Branch("Jet_rawPt",&Jet_rawPt);
    TBranch *nJet_chHEF = splots_tree->Branch("Jet_chHEF",&Jet_chHEF);
    TBranch *nJet_neHEF = splots_tree->Branch("Jet_neHEF",&Jet_neHEF);
    TBranch *nJet_chEmEF = splots_tree->Branch("Jet_chEmEF",&Jet_chEmEF);
    TBranch *nJet_neEmEF = splots_tree->Branch("Jet_neEmEF",&Jet_neEmEF);
    TBranch *nJet_muEF = splots_tree->Branch("Jet_muEF",&Jet_muEF);
    TBranch *nJet_leadTrackPt = splots_tree->Branch("Jet_leadTrackPt",&Jet_leadTrackPt);
    TBranch *nJet_leptonPt = splots_tree->Branch("Jet_leptonPt",&Jet_leptonPt);
    TBranch *nJet_leptonPtRel = splots_tree->Branch("Jet_leptonPtRel",&Jet_leptonPtRel);
    TBranch *nJet_leptonPdgId = splots_tree->Branch("Jet_leptonPdgId",&Jet_leptonPdgId);
    TBranch *nJet_leptonPtRelInv = splots_tree->Branch("Jet_leptonPtRelInv",&Jet_leptonPtRelInv);
    TBranch *nJet_vtxMass = splots_tree->Branch("Jet_vtxMass",&Jet_vtxMass);
    TBranch *nJet_vtxNtracks = splots_tree->Branch("Jet_vtxNtracks",&Jet_vtxNtracks);
    TBranch *nJet_vtxPt = splots_tree->Branch("Jet_vtxPt",&Jet_vtxPt);
    TBranch *nJet_vtx3DSig = splots_tree->Branch("Jet_vtx3DSig",&Jet_vtx3DSig);
    TBranch *nJet_vtx3DVal = splots_tree->Branch("Jet_vtx3DVal",&Jet_vtx3DVal);
    TBranch *nJet_ptd = splots_tree->Branch("Jet_ptd",&Jet_ptd);


    TH1F *p_bdt = new TH1F("r_p_bdt" ,"jet BDT" ,100 , -1, 1);
    TH1F *p_20_bdt = new TH1F("r_p_20_bdt" ,"jet BDT pt_20_30" ,100 , -1, 1);
    TH1F *p_30_bdt = new TH1F("r_p_30_bdt" ,"jet BDT pt_30_40" ,100 , -1, 1);
    TH1F *p_40_bdt = new TH1F("r_p_40_bdt" ,"jet BDT pt_40_60" ,100 , -1, 1);
    TH1F *p_60_bdt = new TH1F("r_p_60_bdt" ,"jet BDT pt_60_80" ,100 , -1, 1);
    TH1F *p_80_bdt = new TH1F("r_p_80_bdt" ,"jet BDT pt_80_100" ,100 , -1, 1);
    TH1F *p_100_bdt = new TH1F("r_p_100_bdt" ,"jet BDT pt_100_120" ,100 , -1, 1);
    TH1F *p_120_bdt = new TH1F("r_p_120_bdt" ,"jet BDT pt_>120" ,100 , -1, 1);

    Long64_t nentries = datatree->GetEntries();
    cout<<nentries<<endl;
    for (Long64_t ientry=0; ientry < nentries; ientry++)
    {
        datatree->GetEntry(ientry);
           
        for (Int_t i = 0; i < injet;i++)
        {
            BDT_response = (*iBDT_response)[i];
            jet_pt = (*ijet_pt)[i];
            jet_indices = (*ijet_indices)[i];
            m_nl_j = (*im_nl_j)[i];

            p_bdt->Fill((*iBDT_response)[i]);
   
            if (i==0)
            {
                rho = irho;
                run = irun ;
                lumi = ilumi;
                evt = ievt;             
            }
            else
            {
                rho = -100;
                run = -100;
                lumi = -100;
                evt = -100;
            }
            rho_copy = irho;        
            run_copy = irun ;
            lumi_copy = ilumi;
            evt_copy = ievt;

            Jet_leptonDeltaR = iJet_leptonDeltaR[jet_indices];
            Jet_mass = iJet_mass[jet_indices];
            Jet_pt = iJet_pt[jet_indices];
            Jet_eta = iJet_eta[jet_indices];
            Jet_phi = iJet_phi[jet_indices];
            Jet_rawEnergy = iJet_rawEnergy[jet_indices];
            Jet_energyRing_dR0_neut = iJet_energyRing_dR0_neut[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR1_neut = iJet_energyRing_dR1_neut[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR2_neut = iJet_energyRing_dR2_neut[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR3_neut = iJet_energyRing_dR3_neut[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR4_neut = iJet_energyRing_dR4_neut[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR0_ch = iJet_energyRing_dR0_ch[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR1_ch = iJet_energyRing_dR1_ch[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR2_ch = iJet_energyRing_dR2_ch[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR3_ch = iJet_energyRing_dR3_ch[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR4_ch = iJet_energyRing_dR4_ch[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR0_em = iJet_energyRing_dR0_em[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR1_em = iJet_energyRing_dR1_em[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR2_em = iJet_energyRing_dR2_em[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR3_em = iJet_energyRing_dR3_em[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR4_em = iJet_energyRing_dR4_em[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR0_mu = iJet_energyRing_dR0_mu[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR1_mu = iJet_energyRing_dR1_mu[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR2_mu = iJet_energyRing_dR2_mu[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR3_mu = iJet_energyRing_dR3_mu[jet_indices]/Jet_rawEnergy;
            Jet_energyRing_dR4_mu = iJet_energyRing_dR4_mu[jet_indices]/Jet_rawEnergy;
            Jet_numDaughters_pt03 = iJet_numDaughters_pt03[jet_indices];
            Jet_numberOfDaughters = iJet_numberOfDaughters[jet_indices];
            Jet_rawPt = iJet_rawPt[jet_indices];
            Jet_chHEF = iJet_chHEF[jet_indices];
            Jet_muEF = iJet_muEF[jet_indices];
            Jet_neHEF = iJet_neHEF[jet_indices];
            Jet_chEmEF = iJet_chEmEF[jet_indices];
            Jet_neEmEF = iJet_neEmEF[jet_indices];
            Jet_leadTrackPt = iJet_leadTrackPt[jet_indices];
            Jet_leptonPt = iJet_leptonPt[jet_indices];
            Jet_leptonPtRel = iJet_leptonPtRel[jet_indices];
            Jet_leptonDeltaR = iJet_leptonDeltaR[jet_indices];
            Jet_leptonPdgId = iJet_leptonPdgId[jet_indices];
            Jet_leptonPtRelInv = iJet_leptonPtRelInv[jet_indices];
            Jet_vtxMass = iJet_vtxMass[jet_indices];
            Jet_vtxNtracks = iJet_vtxNtracks[jet_indices];
            Jet_vtxPt = iJet_vtxPt[jet_indices];
            Jet_vtx3DSig = iJet_vtx3DSig[jet_indices];
            Jet_vtx3DVal = iJet_vtx3DVal[jet_indices];
            Jet_ptd = iJet_ptd[jet_indices];


            if(((*ijet_pt)[i]>=20) && ((*ijet_pt)[i]<30))
            {  
                p_20_bdt->Fill((*iBDT_response)[i]);  
                BDT_response_20 = (*iBDT_response)[i];                             
            }
    
            if(((*ijet_pt)[i]>=30) && ((*ijet_pt)[i]<40))
            {
                p_30_bdt->Fill((*iBDT_response)[i]);
                BDT_response_30 = (*iBDT_response)[i];
            }

            if(((*ijet_pt)[i]>=40) && ((*ijet_pt)[i]<60))
            {
                p_40_bdt->Fill((*iBDT_response)[i]);
                BDT_response_40 = (*iBDT_response)[i];
            }

            if(((*ijet_pt)[i]>=60) && ((*ijet_pt)[i]<80))
            {
                p_60_bdt->Fill((*iBDT_response)[i]);
                BDT_response_60 = (*iBDT_response)[i];
            }

            if(((*ijet_pt)[i]>=80) && ((*ijet_pt)[i]<100))
            {
                p_80_bdt->Fill((*iBDT_response)[i]);
                BDT_response_80 = (*iBDT_response)[i];
            }

           if(((*ijet_pt)[i]>=100) && ((*ijet_pt)[i]<120))
            {
                p_100_bdt->Fill((*iBDT_response)[i]);
                BDT_response_100 = (*iBDT_response)[i];
            }      
 
           if((*ijet_pt)[i]>=120)
            {
                p_120_bdt->Fill((*iBDT_response)[i]);
                BDT_response_120 = (*iBDT_response)[i];
            }
           splots_tree->Fill();
        }
    } 


p_bdt->Write();
p_20_bdt->Write();
p_30_bdt->Write();
p_40_bdt->Write();
p_60_bdt->Write();
p_80_bdt->Write();
p_100_bdt->Write();
p_120_bdt->Write();

splots_tree->Write();
outputFile->Close();
//splots_tree->Scan("BDT_response");
//splots_tree->Scan("jet_pt");

}



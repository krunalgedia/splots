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
#include "THStack.h"

using namespace std;
using namespace TMVA;

// Copy histogram from one root file to other :https://root.cern.ch/root/roottalk/roottalk04/2904.html

void copytree()
{

    TFile* r_data = new TFile("MuonEG_49.root");
    TTree* r_datatree = (TTree*)(r_data->Get("tree"));

    TFile *r_outputFile = TFile::Open("r_TMVA_output.root", "RECREATE" );
    TTree *r_tree = r_datatree->CloneTree();

    r_tree->Write();
    r_tree->Print();

    r_data->Close();
    r_outputFile->Close();

    delete r_data;
    delete r_outputFile;
  
}

void r_tmva()
{

    copytree();  

    TFile* r_data = new TFile("r_TMVA_output.root","UPDATE");
    TTree* r_datatree = (TTree*)(r_data->Get("tree"));
      
//---------------------------------------------------------------
// This loads the library
    TMVA::Tools::Instance();

// Create the Reader object

    TMVA::Reader *r_reader_l = new TMVA::Reader();
    TMVA::Reader *r_reader_s = new TMVA::Reader();
    TMVA::Reader *r_reader_o = new TMVA::Reader();

    vector<double>* r_m_nl_j=0;
    vector<double>* r_m_fl_j=0;

    vector<double>* r_eta_nl_j=0;
    vector<double>* r_eta_fl_j=0;
    vector<double>* r_eta_nlj_ll=0;
    vector<double>* r_eta_flj_ll=0;
    vector<double>* r_eta_j_ll=0;

    vector<double>* r_phi_nl_j=0;
    vector<double>* r_phi_fl_j=0;
    vector<double>* r_phi_nlj_ll=0;
    vector<double>* r_phi_flj_ll=0;
    vector<double>* r_phi_j_ll=0;

    Int_t r_Jet_mcFlavour[23];   //[nJet]

    TBranch *r_b_Jet_mcFlavour;

    vector<double>* r_jet_indices=0;

    double r_njet;

    r_datatree->SetBranchAddress("m_nl_j",&r_m_nl_j);
    r_datatree->SetBranchAddress("phi_nl_j",&r_phi_nl_j);
    r_datatree->SetBranchAddress("eta_nl_j",&r_eta_nl_j);
    r_datatree->SetBranchAddress("phi_nlj_ll",&r_phi_nlj_ll);
    r_datatree->SetBranchAddress("eta_nlj_ll",&r_eta_nlj_ll);
    r_datatree->SetBranchAddress("m_fl_j",&r_m_fl_j);
    r_datatree->SetBranchAddress("phi_fl_j",&r_phi_fl_j);
    r_datatree->SetBranchAddress("eta_fl_j",&r_eta_fl_j);
    r_datatree->SetBranchAddress("phi_flj_ll",&r_phi_flj_ll);
    r_datatree->SetBranchAddress("eta_flj_ll",&r_eta_flj_ll);
    r_datatree->SetBranchAddress("phi_j_ll",&r_phi_j_ll);
    r_datatree->SetBranchAddress("eta_j_ll",&r_eta_j_ll);
    r_datatree->SetBranchAddress("njet",&r_njet);

    r_datatree->SetBranchAddress("Jet_mcFlavour", r_Jet_mcFlavour, &r_b_Jet_mcFlavour);
    r_datatree->SetBranchAddress("jet_indices",&r_jet_indices);


    float r_var1,r_var2,r_var3,r_var4,r_var5,r_var6,r_var7,r_var8,r_var9,r_var10, r_var11, r_var12;
    float r_var1s,r_var2s,r_var3s,r_var4s,r_var5s,r_var6s,r_var7s,r_var8s,r_var9s,r_var10s, r_var11s, r_var12s;
    float r_var1o,r_var2o,r_var3o,r_var4o,r_var5o,r_var6o,r_var7o,r_var8o,r_var9o,r_var10o, r_var11o, r_var12o;

    r_reader_l->AddVariable("close_mlj[0]",&r_var1);
    r_reader_l->AddVariable("close_dphi", &r_var2);
    r_reader_l->AddVariable("close_deta", &r_var3);
    r_reader_l->AddVariable("close_lj2ll_dphi", &r_var4);
    r_reader_l->AddVariable("close_lj2ll_deta", &r_var5);
    r_reader_l->AddVariable("far_mlj", &r_var6);
    r_reader_l->AddVariable("far_dphi", &r_var7);
    r_reader_l->AddVariable("far_deta", &r_var8);
    r_reader_l->AddVariable("far_lj2ll_dphi", &r_var9);
    r_reader_l->AddVariable("far_lj2ll_deta", &r_var10);
    r_reader_l->AddVariable("j2ll_dphi", &r_var11);
    r_reader_l->AddVariable("j2ll_deta", &r_var12);

    r_reader_s->AddVariable("close_mlj[0]",&r_var1s);
    r_reader_s->AddVariable("close_dphi", &r_var2s);
    r_reader_s->AddVariable("close_deta", &r_var3s);
    r_reader_s->AddVariable("close_lj2ll_dphi", &r_var4s);
    r_reader_s->AddVariable("close_lj2ll_deta", &r_var5s);
    r_reader_s->AddVariable("far_mlj", &r_var6s);
    r_reader_s->AddVariable("far_dphi", &r_var7s);
    r_reader_s->AddVariable("far_deta", &r_var8s);
    r_reader_s->AddVariable("far_lj2ll_dphi", &r_var9s);
    r_reader_s->AddVariable("far_lj2ll_deta", &r_var10s);
    r_reader_s->AddVariable("j2ll_dphi", &r_var11s);
    r_reader_s->AddVariable("j2ll_deta", &r_var12s);

    r_reader_o->AddVariable("close_mlj[0]",&r_var1o);
    r_reader_o->AddVariable("close_dphi", &r_var2o);
    r_reader_o->AddVariable("close_deta", &r_var3o);
    r_reader_o->AddVariable("close_lj2ll_dphi", &r_var4o);
    r_reader_o->AddVariable("close_lj2ll_deta", &r_var5o);
    r_reader_o->AddVariable("far_mlj", &r_var6o);
    r_reader_o->AddVariable("far_dphi", &r_var7o);
    r_reader_o->AddVariable("far_deta", &r_var8o);
    r_reader_o->AddVariable("far_lj2ll_dphi", &r_var9o);
    r_reader_o->AddVariable("far_lj2ll_deta", &r_var10o);
    r_reader_o->AddVariable("j2ll_dphi", &r_var11o);
    r_reader_o->AddVariable("j2ll_deta", &r_var12o);


    r_reader_l->BookMVA("BDT method","/afs/cern.ch/work/f/fromeo/public/BTagging/KIN_Nominal/leading/TMVAClassification_BDT.weights.xml");
    r_reader_s->BookMVA("BDT method","/afs/cern.ch/work/f/fromeo/public/BTagging/KIN_Nominal/subleading/TMVAClassification_BDT.weights.xml");
    r_reader_o->BookMVA("BDT method","/afs/cern.ch/work/f/fromeo/public/BTagging/KIN_Nominal/others/TMVAClassification_BDT.weights.xml");


    double r_BDT_response_l, r_BDT_response_s;
    vector<double> r_BDT_response_o;
    double r_bdt_response_o_row;
    vector<double> r_BDT_response;
    double r_bdt_response_row;

    TBranch *r_bdt_l = r_datatree->Branch("r_BDT_response_l",&r_BDT_response_l);
    TBranch *r_bdt_s = r_datatree->Branch("r_BDT_response_s",&r_BDT_response_s);
    TBranch *r_bdt_o = r_datatree->Branch("r_BDT_response_o",&r_BDT_response_o);
    TBranch *r_bdt = r_datatree->Branch("r_BDT_response",&r_BDT_response);


    TH1F *r_p_l_bdt = new TH1F("r_p_l_bdt" ,"Leading jet BDT" ,100 , -1, 1);
    TH1F *r_p_s_bdt = new TH1F("r_p_s_bdt" ,"Sub-leading jet BDT" ,100 , -1, 1);
    TH1F *r_p_o_bdt = new TH1F("r_p_o_bdt" ,"Other jets BDT" ,100 , -1, 1);

    Long64_t r_nentries = r_datatree->GetEntries();
    for (Long64_t ientry=0; ientry < r_nentries; ientry++)
    {
        r_datatree->GetEntry(ientry);

        r_BDT_response_l = -2, r_BDT_response_s = -2, r_bdt_response_o_row = -2, r_bdt_response_row = -2;


        for (Int_t i = 0; i < r_njet;i++)
        {

            if ( i==0 )
            {
                r_var1 = fabs((*r_m_nl_j)[0]);
                r_var2 = fabs((*r_phi_nl_j)[0]);
                r_var3 = fabs((*r_eta_nl_j)[0]);
                r_var4 = fabs((*r_phi_nlj_ll)[0]);
                r_var5 = fabs((*r_eta_nlj_ll)[0]);
                r_var6 = fabs((*r_m_fl_j)[0]);
                r_var7 = fabs((*r_phi_fl_j)[0]);
                r_var8 = fabs((*r_eta_fl_j)[0]);
                r_var9 = fabs((*r_phi_flj_ll)[0]);
                r_var10 = fabs((*r_eta_flj_ll)[0]);
                r_var11 = fabs((*r_phi_j_ll)[0]);
                r_var12 = fabs((*r_eta_j_ll)[0]);

                r_BDT_response_l = r_reader_l->EvaluateMVA("BDT method");
                r_bdt_response_row = r_BDT_response_l;
                r_BDT_response.push_back(r_bdt_response_row);
                r_p_l_bdt->Fill(r_BDT_response_l);

            }

            if ( i==1 )
            {
                r_var1s = fabs((*r_m_nl_j)[1]);
                r_var2s = fabs((*r_phi_nl_j)[1]);
                r_var3s = fabs((*r_eta_nl_j)[1]);
                r_var4s = fabs((*r_phi_nlj_ll)[1]);
                r_var5s = fabs((*r_eta_nlj_ll)[1]);
                r_var6s = fabs((*r_m_fl_j)[1]);
                r_var7s = fabs((*r_phi_fl_j)[1]);
                r_var8s = fabs((*r_eta_fl_j)[1]);
                r_var9s = fabs((*r_phi_flj_ll)[1]);
                r_var10s = fabs((*r_eta_flj_ll)[1]);
                r_var11s = fabs((*r_phi_j_ll)[1]);
                r_var12s = fabs((*r_eta_j_ll)[1]);

                r_BDT_response_s = r_reader_s->EvaluateMVA("BDT method");
                r_bdt_response_row = r_BDT_response_s;
                r_BDT_response.push_back(r_bdt_response_row);
                r_p_s_bdt->Fill(r_BDT_response_s);

            }

            if (i>1)
            {
                r_var1o = fabs((*r_m_nl_j)[i]);
                r_var2o = fabs((*r_phi_nl_j)[i]);
                r_var3o = fabs((*r_eta_nl_j)[i]);
                r_var4o = fabs((*r_phi_nlj_ll)[i]);
                r_var5o = fabs((*r_eta_nlj_ll)[i]);
                r_var6o = fabs((*r_m_fl_j)[i]);
                r_var7o = fabs((*r_phi_fl_j)[i]);
                r_var8o = fabs((*r_eta_fl_j)[i]);
                r_var9o = fabs((*r_phi_flj_ll)[i]);
                r_var10o = fabs((*r_eta_flj_ll)[i]);
                r_var11o = fabs((*r_phi_j_ll)[i]);
                r_var12o = fabs((*r_eta_j_ll)[i]);

                r_bdt_response_o_row = r_reader_o->EvaluateMVA("BDT method");
                r_BDT_response_o.push_back(r_bdt_response_o_row);

                r_bdt_response_row = r_bdt_response_o_row;
                r_BDT_response.push_back(r_bdt_response_row);

                r_p_o_bdt->Fill(r_bdt_response_o_row);

            }
        }

        r_bdt_l->Fill();
        r_bdt_s->Fill();
        r_bdt_o->Fill();
        r_bdt->Fill();

        r_BDT_response.erase(r_BDT_response.begin(),r_BDT_response.begin()+r_njet);

        if (r_njet>1)
        {
            r_BDT_response_o.erase(r_BDT_response_o.begin(),r_BDT_response_o.begin()+r_njet-2);
        }
    }


    r_p_l_bdt->Write();
    r_p_s_bdt->Write();
    r_p_o_bdt->Write();
 
    r_datatree->Write();
    r_datatree->Print();


    TString f0("l");
    TString f1("ALP");

    TString r_yaxis("#events");	

    TString yaxis("#events");

    Double_t scale;

/*
    TCanvas *r_b_l_bdt = new TCanvas ("r_b_l_bdt","r_b_l_bdt");
    r_p_l_bdt->SetLineColor(kGreen + 2);
    r_p_l_bdt->SetLineWidth(1);
    r_p_l_bdt->GetXaxis()->SetTitle("BDT Response");
    r_p_l_bdt->GetXaxis()->CenterTitle();
    r_p_l_bdt->GetYaxis()->SetTitle(yaxis);
    r_p_l_bdt->GetYaxis()->CenterTitle();
    r_p_l_bdt->SetStats(0);
    scale = 1/r_p_l_bdt->Integral();
    r_p_l_bdt->Scale(scale);
    r_p_l_bdt->Draw("hist");

    p_l_bdt->SetLineColor(kRed);
    p_l_bdt->SetStats(0);
    scale = 1/p_l_bdt->Integral();
    p_l_bdt->Scale(scale);
    p_l_bdt->Draw("hist same");

    TLegend *r_s_l_bdt = new TLegend(0.6,0.6,0.8,0.8);
    r_s_l_bdt->SetFillColor(kWhite);
    r_s_l_bdt->AddEntry(p_l_bdt,r_mc,f0);
    r_s_l_bdt->AddEntry(r_p_l_bdt,r_r,f1);
    r_s_l_bdt->SetBorderSize(0);
    r_s_l_bdt->Draw("same");
    r_b_l_bdt->Draw();


    TCanvas *r_b_s_bdt = new TCanvas ("r_b_s_bdt","r_b_s_bdt");
    r_p_s_bdt->SetLineColor(kGreen + 2);
    r_p_s_bdt->SetLineWidth(1);
    r_p_s_bdt->GetXaxis()->SetTitle("BDT Response");
    r_p_s_bdt->GetXaxis()->CenterTitle();
    r_p_s_bdt->GetYaxis()->SetTitle(yaxis);
    r_p_s_bdt->GetYaxis()->CenterTitle();
    r_p_s_bdt->SetStats(0);
    scale = 1/r_p_s_bdt->Integral();
    r_p_s_bdt->Scale(scale);
    r_p_s_bdt->Draw("hist");

    p_s_bdt->SetLineColor(kRed);
    p_s_bdt->SetStats(0);
    scale = 1/p_s_bdt->Integral();
    p_s_bdt->Scale(scale);
    p_s_bdt->Draw("hist same");

    TLegend *r_s_s_bdt = new TLegend(0.6,0.6,0.8,0.8);
    r_s_s_bdt->SetFillColor(kWhite);
    r_s_s_bdt->AddEntry(p_s_bdt,r_mc,f0);
    r_s_s_bdt->AddEntry(r_p_s_bdt,r_r,f1);
    r_s_s_bdt->SetBorderSize(0);
    r_s_s_bdt->Draw("same");
    r_b_s_bdt->Draw();



    TCanvas *r_b_o_bdt = new TCanvas ("r_b_o_bdt","r_b_o_bdt");
    r_p_o_bdt->SetLineColor(kGreen + 2);
    r_p_o_bdt->SetLineWidth(1);
    r_p_o_bdt->GetXaxis()->SetTitle("BDT Response");
    r_p_o_bdt->GetXaxis()->CenterTitle();
    r_p_o_bdt->GetYaxis()->SetTitle(yaxis);
    r_p_o_bdt->GetYaxis()->CenterTitle();
    r_p_o_bdt->SetStats(0);
    scale = 1/r_p_o_bdt->Integral();
    r_p_o_bdt->Scale(scale);
    r_p_o_bdt->Draw("hist");

    p_o_bdt->SetLineColor(kRed);
    p_o_bdt->SetStats(0);
    scale = 1/p_o_bdt->Integral();
    p_o_bdt->Scale(scale);
    p_o_bdt->Draw("hist same");

    TLegend *r_s_o_bdt = new TLegend(0.6,0.6,0.8,0.8);
    r_s_o_bdt->SetFillColor(kWhite);
    r_s_o_bdt->AddEntry(p_o_bdt,r_mc,f0);
    r_s_o_bdt->AddEntry(r_p_o_bdt,r_r,f1);
    r_s_o_bdt->SetBorderSize(0);
    r_s_o_bdt->Draw("same");
    r_b_s_bdt->Draw();

*/
///////////////////////////////////////////////////////////////
/*
    TCanvas *b_l_bdt = new TCanvas ("b_l_bdt","b_l_bdt");
    p_ln_bdt->SetLineColor(kBlue);
    p_ln_bdt->SetLineWidth(1);
    p_ln_bdt->GetXaxis()->SetTitle("BDT Response");
    p_ln_bdt->GetXaxis()->CenterTitle();
    p_ln_bdt->GetYaxis()->SetTitle(yaxis);
    p_ln_bdt->GetYaxis()->CenterTitle();
    TString bin_width_l = p_ln_bdt->GetYaxis()->GetBinWidth(0);
    p_ln_bdt->SetStats(0);
    scale = 1/p_ln_bdt->Integral();
    p_ln_bdt->Scale(scale);
    p_ln_bdt->Draw("hist");

    p_lb_bdt->SetLineColor(kRed);
    p_lb_bdt->SetStats(0);
    scale = 1/p_lb_bdt->Integral();
    p_lb_bdt->Scale(scale);
    p_lb_bdt->Draw("hist same");

    TLegend *s_l_bdt = new TLegend(0.6,0.6,0.8,0.8);
    s_l_bdt->SetFillColor(kWhite);
    s_l_bdt->AddEntry(p_lb_bdt,blead,f0);
    s_l_bdt->AddEntry(p_ln_bdt,nlead,f0);
    s_l_bdt->SetBorderSize(0);
    s_l_bdt->Draw("same");
    b_l_bdt->Draw();

///////////////////////////////////////////////////////////////

    TCanvas *b_s_bdt = new TCanvas ("b_s_bdt","b_s_bdt");
    p_sn_bdt->SetLineColor(kBlue);
    p_sn_bdt->SetLineWidth(1);
    p_sn_bdt->GetXaxis()->SetTitle("BDT Response");
    p_sn_bdt->GetXaxis()->CenterTitle();
    p_sn_bdt->GetYaxis()->SetTitle(yaxis);
    p_sn_bdt->GetYaxis()->CenterTitle();
  //  p_sn_bdt->GetYaxis()->SetRangeUser(0,150);
    scale = 1/p_sn_bdt->Integral();
    p_sn_bdt->Scale(scale);
    p_sn_bdt->SetStats(0);
    p_sn_bdt->Draw("hist");

    p_sb_bdt->SetLineColor(kRed);
    p_sb_bdt->SetStats(0);
    scale = 1/p_sb_bdt->Integral();
    p_sb_bdt->Scale(scale);
    p_sb_bdt->Draw("hist same");

    TLegend *s_s_bdt = new TLegend(0.6,0.6,0.8,0.8);
    s_s_bdt->SetFillColor(kWhite);
    s_s_bdt->AddEntry(p_sb_bdt,bsub,f0);
    s_s_bdt->AddEntry(p_sn_bdt,nsub,f0);
    s_s_bdt->SetBorderSize(0);
    s_s_bdt->Draw("same");
    b_s_bdt->Draw();

//////////////////////////////////////////////////////////////////

    TCanvas *b_o_bdt = new TCanvas ("b_o_bdt","b_o_bdt");
    p_ob_bdt->SetLineColor(kRed);
    p_ob_bdt->SetLineWidth(1);
    p_ob_bdt->GetXaxis()->SetTitle("BDT Response");
    p_ob_bdt->GetXaxis()->CenterTitle();
    p_ob_bdt->GetYaxis()->SetTitle(yaxis);
    p_ob_bdt->GetYaxis()->CenterTitle();
    p_ob_bdt->GetYaxis()->SetRangeUser(0,60);
    p_ob_bdt->SetStats(0);
    scale = 1/p_ob_bdt->Integral();
    p_ob_bdt->Scale(scale);
    p_ob_bdt->Draw("hist");

    p_on_bdt->SetLineColor(kBlue);
    p_on_bdt->SetStats(0);
    scale = 1/p_on_bdt->Integral();
    p_on_bdt->Scale(scale);
    p_on_bdt->Draw("hist same");

    TLegend *s_o_bdt = new TLegend(0.6,0.6,0.8,0.8);
    s_o_bdt->SetFillColor(kWhite);
    s_o_bdt->AddEntry(p_ob_bdt,bother,f0);
    s_o_bdt->AddEntry(p_on_bdt,nother,f0);
    s_o_bdt->SetBorderSize(0);
    s_o_bdt->Draw("same");
    b_o_bdt->Draw();

*/
/////////////////////////////////////////////////////////////////
/*
    TCanvas *b_lb_bdt = new TCanvas ("b_lb_bdt","b_lb_bdt");
    p_lb_bdt->Fit("gaus","L");
    p_lb_bdt->Draw();
    b_lb_bdt->Draw();

    TCanvas *b_ln_bdt = new TCanvas ("b_ln_bdt","b_ln_bdt");
    p_ln_bdt->Fit("gaus","L");//,"",-0.95,1.0);
//    p_l_bdt->Fit("pol2","L","",-1.0,-0.95);
    p_ln_bdt->Draw();
    b_ln_bdt->Draw();
*/

} 

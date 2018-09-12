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
    TFile* data = new TFile("TTbar_14.root");
    TTree* datatree = (TTree*)(data->Get("tree"));
  
    TFile *outputFile = TFile::Open("TMVA_output.root", "RECREATE" );
    TTree *tree = datatree->CloneTree();
   
    tree->Write();
    tree->Print();
  
    data->Close();
    outputFile->Close();
   
    delete data;
    delete outputFile;
  
}

void tmva()
{

    copytree();  

    TFile* data = new TFile("TMVA_output.root","UPDATE");
    TTree* datatree = (TTree*)(data->Get("tree"));
      
//---------------------------------------------------------------
// This loads the library
    TMVA::Tools::Instance();

// Create the Reader object
    TMVA::Reader *reader_l = new TMVA::Reader();
    TMVA::Reader *reader_s = new TMVA::Reader();
    TMVA::Reader *reader_o = new TMVA::Reader();

    vector<double>* m_nl_j=0;
    vector<double>* m_fl_j=0;

    vector<double>* eta_nl_j=0;
    vector<double>* eta_fl_j=0;
    vector<double>* eta_nlj_ll=0;
    vector<double>* eta_flj_ll=0;
    vector<double>* eta_j_ll=0;

    vector<double>* phi_nl_j=0;
    vector<double>* phi_fl_j=0;
    vector<double>* phi_nlj_ll=0;
    vector<double>* phi_flj_ll=0;
    vector<double>* phi_j_ll=0;

    Int_t Jet_mcFlavour[23];   //[nJet]

    TBranch *b_Jet_mcFlavour;   

    vector<double>* jet_indices=0;

    double njet;


    datatree->SetBranchAddress("m_nl_j",&m_nl_j);
    datatree->SetBranchAddress("phi_nl_j",&phi_nl_j);
    datatree->SetBranchAddress("eta_nl_j",&eta_nl_j);
    datatree->SetBranchAddress("phi_nlj_ll",&phi_nlj_ll);
    datatree->SetBranchAddress("eta_nlj_ll",&eta_nlj_ll);
    datatree->SetBranchAddress("m_fl_j",&m_fl_j);
    datatree->SetBranchAddress("phi_fl_j",&phi_fl_j);
    datatree->SetBranchAddress("eta_fl_j",&eta_fl_j);
    datatree->SetBranchAddress("phi_flj_ll",&phi_flj_ll);
    datatree->SetBranchAddress("eta_flj_ll",&eta_flj_ll);
    datatree->SetBranchAddress("phi_j_ll",&phi_j_ll);
    datatree->SetBranchAddress("eta_j_ll",&eta_j_ll);
    datatree->SetBranchAddress("njet",&njet);

    datatree->SetBranchAddress("Jet_mcFlavour", Jet_mcFlavour, &b_Jet_mcFlavour);
    datatree->SetBranchAddress("jet_indices",&jet_indices);


    float var1,var2,var3,var4,var5,var6,var7,var8,var9,var10, var11, var12;
    float var1s,var2s,var3s,var4s,var5s,var6s,var7s,var8s,var9s,var10s, var11s, var12s;
    float var1o,var2o,var3o,var4o,var5o,var6o,var7o,var8o,var9o,var10o, var11o, var12o;

    float r_var1,r_var2,r_var3,r_var4,r_var5,r_var6,r_var7,r_var8,r_var9,r_var10, r_var11, r_var12;
    float r_var1s,r_var2s,r_var3s,r_var4s,r_var5s,r_var6s,r_var7s,r_var8s,r_var9s,r_var10s, r_var11s, r_var12s;
    float r_var1o,r_var2o,r_var3o,r_var4o,r_var5o,r_var6o,r_var7o,r_var8o,r_var9o,r_var10o, r_var11o, r_var12o;


    reader_l->AddVariable("close_mlj[0]",&var1); 
    reader_l->AddVariable("close_dphi", &var2);
    reader_l->AddVariable("close_deta", &var3);
    reader_l->AddVariable("close_lj2ll_dphi", &var4);
    reader_l->AddVariable("close_lj2ll_deta", &var5);
    reader_l->AddVariable("far_mlj", &var6);
    reader_l->AddVariable("far_dphi", &var7);
    reader_l->AddVariable("far_deta", &var8);
    reader_l->AddVariable("far_lj2ll_dphi", &var9);
    reader_l->AddVariable("far_lj2ll_deta", &var10);
    reader_l->AddVariable("j2ll_dphi", &var11);
    reader_l->AddVariable("j2ll_deta", &var12);

    reader_s->AddVariable("close_mlj[0]",&var1s);
    reader_s->AddVariable("close_dphi", &var2s);
    reader_s->AddVariable("close_deta", &var3s);
    reader_s->AddVariable("close_lj2ll_dphi", &var4s);
    reader_s->AddVariable("close_lj2ll_deta", &var5s);
    reader_s->AddVariable("far_mlj", &var6s);
    reader_s->AddVariable("far_dphi", &var7s);
    reader_s->AddVariable("far_deta", &var8s);
    reader_s->AddVariable("far_lj2ll_dphi", &var9s);
    reader_s->AddVariable("far_lj2ll_deta", &var10s);
    reader_s->AddVariable("j2ll_dphi", &var11s);
    reader_s->AddVariable("j2ll_deta", &var12s);

    reader_o->AddVariable("close_mlj[0]",&var1o);
    reader_o->AddVariable("close_dphi", &var2o);
    reader_o->AddVariable("close_deta", &var3o);
    reader_o->AddVariable("close_lj2ll_dphi", &var4o);
    reader_o->AddVariable("close_lj2ll_deta", &var5o);
    reader_o->AddVariable("far_mlj", &var6o);
    reader_o->AddVariable("far_dphi", &var7o);
    reader_o->AddVariable("far_deta", &var8o);
    reader_o->AddVariable("far_lj2ll_dphi", &var9o);
    reader_o->AddVariable("far_lj2ll_deta", &var10o);
    reader_o->AddVariable("j2ll_dphi", &var11o);
    reader_o->AddVariable("j2ll_deta", &var12o);


    reader_l->BookMVA("BDT method","/afs/cern.ch/work/f/fromeo/public/BTagging/KIN_Nominal/leading/TMVAClassification_BDT.weights.xml");
    reader_s->BookMVA("BDT method","/afs/cern.ch/work/f/fromeo/public/BTagging/KIN_Nominal/subleading/TMVAClassification_BDT.weights.xml");
    reader_o->BookMVA("BDT method","/afs/cern.ch/work/f/fromeo/public/BTagging/KIN_Nominal/others/TMVAClassification_BDT.weights.xml");


    double BDT_response_l, BDT_response_s;
    double BDT_response_lb, BDT_response_ln;
    double BDT_response_sb, BDT_response_sn;
    vector<double> BDT_response_o;
    double bdt_response_o_row;
    vector<double> BDT_response_ob;
    double bdt_response_ob_row;
    vector<double> BDT_response_on;
    double bdt_response_on_row;
    vector<double> BDT_response;
    double bdt_response_row;

    TBranch *bdt_l = datatree->Branch("BDT_response_l",&BDT_response_l);
    TBranch *bdt_lb = datatree->Branch("BDT_response_lb",&BDT_response_lb);
    TBranch *bdt_ln = datatree->Branch("BDT_response_ln",&BDT_response_ln);

    TBranch *bdt_s = datatree->Branch("BDT_response_s",&BDT_response_s);
    TBranch *bdt_sb = datatree->Branch("BDT_response_sb",&BDT_response_sb);
    TBranch *bdt_sn = datatree->Branch("BDT_response_sn",&BDT_response_sn);

    TBranch *bdt_o = datatree->Branch("BDT_response_o",&BDT_response_o);
    TBranch *bdt_ob = datatree->Branch("BDT_response_ob",&BDT_response_ob);
    TBranch *bdt_on = datatree->Branch("BDT_response_on",&BDT_response_on);

    TBranch *bdt = datatree->Branch("BDT_response",&BDT_response);

    TH1F *p_l_bdt = new TH1F("p_l_bdt" ,"Leading jet BDT" ,100 , -1, 1);
    TH1F *p_lb_bdt = new TH1F("p_lb_bdt" ,"Leading bjet BDT" ,100 , -1, 1);
    TH1F *p_ln_bdt = new TH1F("p_ln_bdt" ,"Leading non-bjet BDT" ,100 , -1, 1);
    TH1F *p_s_bdt = new TH1F("p_s_bdt" ,"Sub-leading jet BDT" ,100 , -1, 1);
    TH1F *p_sb_bdt = new TH1F("p_sb_bdt" ,"Sub-leading jet BDT" ,100 , -1, 1);
    TH1F *p_sn_bdt = new TH1F("p_sn_bdt" ,"Sub-leading jet BDT" ,100 , -1, 1);
    TH1F *p_o_bdt = new TH1F("p_o_bdt" ,"Other jets BDT" ,100 , -1, 1);
    TH1F *p_ob_bdt = new TH1F("p_ob_bdt" ,"Other bjets BDT" ,100 , -1, 1);
    TH1F *p_on_bdt = new TH1F("p_on_bdt" ,"Other non-bjets BDT" ,100 , -1, 1);


    Long64_t nentries = datatree->GetEntries();
    int check_ob = 0, check_on = 0;
   
    for (Long64_t ientry=0; ientry < nentries; ientry++)
    {
        datatree->GetEntry(ientry);

        BDT_response_l = -2, BDT_response_s = -2, BDT_response_lb = -2, BDT_response_ln = -2, BDT_response_sb = -2, BDT_response_sn = -2, bdt_response_o_row = -2, bdt_response_ob_row = -2, bdt_response_on_row = -2, bdt_response_row = -2;
       
        check_ob = 0, check_on =0;        

        for (Int_t i = 0; i < njet;i++)
        {
           
            if ( i==0 )
            {
                var1 = fabs((*m_nl_j)[0]);
                var2 = fabs((*phi_nl_j)[0]);
                var3 = fabs((*eta_nl_j)[0]); 
       	        var4 = fabs((*phi_nlj_ll)[0]);
                var5 = fabs((*eta_nlj_ll)[0]);
                var6 = fabs((*m_fl_j)[0]);
                var7 = fabs((*phi_fl_j)[0]);
                var8 = fabs((*eta_fl_j)[0]);
                var9 = fabs((*phi_flj_ll)[0]);
                var10 = fabs((*eta_flj_ll)[0]);
                var11 = fabs((*phi_j_ll)[0]);
                var12 = fabs((*eta_j_ll)[0]);

                BDT_response_l = reader_l->EvaluateMVA("BDT method");
                bdt_response_row = BDT_response_l;
                BDT_response.push_back(bdt_response_row);
                p_l_bdt->Fill(BDT_response_l);

                int t = (*jet_indices)[0];

                if ((fabs(Jet_mcFlavour[t])) == 5)
                {   
                    p_lb_bdt->Fill(BDT_response_l);
                    BDT_response_lb = BDT_response_l;
                }
                else
                {               
                    p_ln_bdt->Fill(BDT_response_l);
                    BDT_response_ln = BDT_response_l;
                }

            }       
           
           
            if (i==1)
            {
                var1s = (*m_nl_j)[1];
                var2s = fabs((*phi_nl_j)[1]);
                var3s = fabs((*eta_nl_j)[1]);
                var4s = fabs((*phi_nlj_ll)[1]);
                var5s = fabs((*eta_nlj_ll)[1]);
                var6s = fabs((*m_fl_j)[1]);
                var7s = fabs((*phi_fl_j)[1]);
                var8s = fabs((*eta_fl_j)[1]);
                var9s = fabs((*phi_flj_ll)[1]);
                var10s = fabs((*eta_flj_ll)[1]);
                var11s = fabs((*phi_j_ll)[1]);
                var12s = fabs((*eta_j_ll)[1]);

		BDT_response_s = reader_s->EvaluateMVA("BDT method");
                bdt_response_row = BDT_response_s;
                BDT_response.push_back(bdt_response_row);
                p_s_bdt->Fill(BDT_response_s);

                int t = (*jet_indices)[1];
                if (((fabs(Jet_mcFlavour[t])) == 5))
                {
                    p_sb_bdt->Fill(BDT_response_s);
                    BDT_response_sb = BDT_response_s;
                }
                else
                {
                    p_sn_bdt->Fill(BDT_response_s);
                    BDT_response_sn = BDT_response_s;
                }

            }

            if (i>1)
            {
                var1o = fabs((*m_nl_j)[i]);
                var2o = fabs((*phi_nl_j)[i]);
                var3o = fabs((*eta_nl_j)[i]);
                var4o = fabs((*phi_nlj_ll)[i]);
                var5o = fabs((*eta_nlj_ll)[i]);
                var6o = fabs((*m_fl_j)[i]);
                var7o = fabs((*phi_fl_j)[i]);
                var8o = fabs((*eta_fl_j)[i]);
                var9o = fabs((*phi_flj_ll)[i]);
                var10o = fabs((*eta_flj_ll)[i]);         
                var11o = fabs((*phi_j_ll)[i]);
                var12o = fabs((*eta_j_ll)[i]);         
         
                bdt_response_o_row = reader_o->EvaluateMVA("BDT method");
                BDT_response_o.push_back(bdt_response_o_row);
                   
                bdt_response_row = bdt_response_o_row;
                BDT_response.push_back(bdt_response_row);

                p_o_bdt->Fill(bdt_response_o_row);

                int t = (*jet_indices)[i];

                if (((fabs(Jet_mcFlavour[t])) == 5))
                {
                    p_ob_bdt->Fill(bdt_response_o_row);
                    BDT_response_ob.push_back(bdt_response_o_row);
                    check_ob = check_ob + 1;
                }
                else
                {
                    p_on_bdt->Fill(bdt_response_o_row);
                    BDT_response_on.push_back(bdt_response_o_row);
                    check_on = check_on + 1;
                }
            }                 
        }
                
        bdt_l->Fill();
        bdt_lb->Fill();
        bdt_ln->Fill();
        bdt_s->Fill();
        bdt_sb->Fill();
        bdt_sn->Fill();
        bdt_o->Fill();
        bdt_ob->Fill();
        bdt_on->Fill();
        bdt->Fill();
       // datatree->Fill(); 

        BDT_response.erase(BDT_response.begin(),BDT_response.begin()+njet);

        if (njet>1)
        { 
        BDT_response_o.erase(BDT_response_o.begin(),BDT_response_o.begin()+njet-2); 
        if (check_ob >0)
            {
            BDT_response_ob.erase(BDT_response_ob.begin(),BDT_response_ob.begin()+check_ob);
            }
        if (check_ob >0)
            {
            BDT_response_on.erase(BDT_response_on.begin(),BDT_response_on.begin()+check_on);
            }
        }
    }

    p_l_bdt->Write();   
    p_lb_bdt->Write();
    p_ln_bdt->Write();
    p_s_bdt->Write();
    p_sb_bdt->Write();
    p_sn_bdt->Write();
    p_o_bdt->Write();
    p_ob_bdt->Write();
    p_on_bdt->Write();

    datatree->Write();
    datatree->Print();


    TString blead("leading bjet"),f0("l");
    TString nlead("leading non-bjet");

    TString bsub("subleading bjet");
    TString nsub("subleading non-bjet");

    TString bother("other bjets");
    TString nother("other non-bjets");

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

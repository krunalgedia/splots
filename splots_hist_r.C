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

void splots_hist_r()
{
    TFile* data = new TFile("r_TMVA_output.root","UPDATE");
    TTree* datatree = (TTree*)(data->Get("tree"));
 
    double r_BDT_response_l=0;
    double r_BDT_response_s=0;

    vector<double>* r_BDT_response_o=0;
    vector<double>* r_BDT_response=0;

    vector<double>* r_jet_indices=0;

    double r_njet;

    datatree->SetBranchAddress("r_BDT_response",&r_BDT_response);
    datatree->SetBranchAddress("r_BDT_response_l",&r_BDT_response_l);
    datatree->SetBranchAddress("r_BDT_response_s",&r_BDT_response_s);
    datatree->SetBranchAddress("r_BDT_response_o",&r_BDT_response_o);

    datatree->SetBranchAddress("njet",&r_njet);
    datatree->SetBranchAddress("jet_indices",&r_jet_indices);

    TH1F *r_splot_p_l_bdt = new TH1F("r_splot_p_l_bdt" ,"r_Leading jet BDT" ,300 , -1, 1);
    TH1F *r_splot_p_s_bdt = new TH1F("r_splot_p_s_bdt" ,"r_Sub-leading jet BDT" ,300 , -1, 1);
    TH1F *r_splot_p_o_bdt = new TH1F("r_splot_p_o_bdt" ,"r_Other jets BDT" ,300 , -1, 1);

    Long64_t nentries = datatree->GetEntries();

    for (Long64_t ientry=0; ientry < nentries; ientry++)
    {
        datatree->GetEntry(ientry);

        for (Int_t i = 0; i < r_njet;i++)
        {
       
            if ( i==0 )
            {
                r_splot_p_l_bdt->Fill((*r_BDT_response)[0]);
            }

            if (i==1)
            {
                r_splot_p_s_bdt->Fill((*r_BDT_response)[1]);
            }

            if (i>1)
            {
                r_splot_p_o_bdt->Fill((*r_BDT_response)[i]);
            }
        }
    }

    r_splot_p_l_bdt->Write();
    r_splot_p_s_bdt->Write();
    r_splot_p_o_bdt->Write();

    datatree->Write();
//    datatree->Print();
/*
    TString blead("leading bjet"),f0("l");
    TString nlead("leading non-bjet");

    TString bsub("subleading bjet");
    TString nsub("subleading non-bjet");

    TString bother("other bjets");
    TString nother("other non-bjets");

    TString yaxis("#events");

    Double_t scale;

    TCanvas *r_splot_b_l_bdt = new TCanvas ("r_splot_b_l_bdt","r_splot_b_l_bdt");
    splot_p_ln_bdt->SetLineColor(kBlue);
    splot_p_ln_bdt->SetLineWidth(1);
    splot_p_ln_bdt->GetXaxis()->SetTitle("BDT Response");
    splot_p_ln_bdt->GetXaxis()->CenterTitle();
    splot_p_ln_bdt->GetYaxis()->SetTitle(yaxis);
    splot_p_ln_bdt->GetYaxis()->CenterTitle();
//    TString bin_width_l = p_ln_bdt->GetYaxis()->GetBinWidth(0);
    splot_p_ln_bdt->SetStats(0);
    scale = 1/splot_p_ln_bdt->Integral();
    splot_p_ln_bdt->Scale(scale);
    splot_p_ln_bdt->Draw("hist");

    splot_p_lb_bdt->SetLineColor(kRed);
    splot_p_lb_bdt->SetStats(0);
    scale = 1/splot_p_lb_bdt->Integral();
    splot_p_lb_bdt->Scale(scale);
    splot_p_lb_bdt->Draw("hist same");

    TLegend *splot_s_l_bdt = new TLegend(0.6,0.6,0.8,0.8);
    splot_s_l_bdt->SetFillColor(kWhite);
    splot_s_l_bdt->AddEntry(splot_p_lb_bdt,blead,f0);
    splot_s_l_bdt->AddEntry(splot_p_ln_bdt,nlead,f0);
    splot_s_l_bdt->SetBorderSize(0);
    splot_s_l_bdt->Draw("same");
    splot_b_l_bdt->Draw();



    TCanvas *splot_b_s_bdt = new TCanvas ("splot_b_s_bdt","splot_b_s_bdt");
    splot_p_sn_bdt->SetLineColor(kBlue);
    splot_p_sn_bdt->SetLineWidth(1);
    splot_p_sn_bdt->GetXaxis()->SetTitle("BDT Response");
    splot_p_sn_bdt->GetXaxis()->CenterTitle();
    splot_p_sn_bdt->GetYaxis()->SetTitle(yaxis);
    splot_p_sn_bdt->GetYaxis()->CenterTitle();
    scale = 1/splot_p_sn_bdt->Integral();
    splot_p_sn_bdt->Scale(scale);
    splot_p_sn_bdt->SetStats(0);
    splot_p_sn_bdt->Draw("hist");

    splot_p_sb_bdt->SetLineColor(kRed);
    splot_p_sb_bdt->SetStats(0);
    scale = 1/splot_p_sb_bdt->Integral();
    splot_p_sb_bdt->Scale(scale);
    splot_p_sb_bdt->Draw("hist same");

    TLegend *splot_s_s_bdt = new TLegend(0.6,0.6,0.8,0.8);
    splot_s_s_bdt->SetFillColor(kWhite);
    splot_s_s_bdt->AddEntry(splot_p_sb_bdt,bsub,f0);
    splot_s_s_bdt->AddEntry(splot_p_sn_bdt,nsub,f0);
    splot_s_s_bdt->SetBorderSize(0);
    splot_s_s_bdt->Draw("same");
    splot_b_s_bdt->Draw();




    TCanvas *splot_b_o_bdt = new TCanvas ("splot_b_o_bdt","splot_b_o_bdt");
    splot_p_ob_bdt->SetLineColor(kRed);
    splot_p_ob_bdt->SetLineWidth(1);
    splot_p_ob_bdt->GetXaxis()->SetTitle("BDT Response");
    splot_p_ob_bdt->GetXaxis()->CenterTitle();
    splot_p_ob_bdt->GetYaxis()->SetTitle(yaxis);
    splot_p_ob_bdt->GetYaxis()->CenterTitle();
//    splot_p_ob_bdt->GetYaxis()->SetRangeUser(0,60);
    splot_p_ob_bdt->SetStats(0);
    scale = 1/splot_p_ob_bdt->Integral();
    splot_p_ob_bdt->Scale(scale);
    splot_p_ob_bdt->Draw("hist");

    splot_p_on_bdt->SetLineColor(kBlue);
    splot_p_on_bdt->SetStats(0);
    scale = 1/splot_p_on_bdt->Integral();
    splot_p_on_bdt->Scale(scale);
    splot_p_on_bdt->Draw("hist same");

    TLegend *splot_s_o_bdt = new TLegend(0.6,0.6,0.8,0.8);
    splot_s_o_bdt->SetFillColor(kWhite);
    splot_s_o_bdt->AddEntry(splot_p_ob_bdt,bother,f0);
    splot_s_o_bdt->AddEntry(splot_p_on_bdt,nother,f0);
    splot_s_o_bdt->SetBorderSize(0);
    splot_s_o_bdt->Draw("same");
    splot_b_o_bdt->Draw();

*/
}

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

void splots_hist_mc()
{
    TFile* data = new TFile("TMVA_output.root","UPDATE");
    TTree* datatree = (TTree*)(data->Get("tree"));
 
    double BDT_response_l=0;
    double BDT_response_lb=0;
    double BDT_response_ln=0;

    double BDT_response_s=0;
    double BDT_response_sb=0;
    double BDT_response_sn=0;

    vector<double>* BDT_response_o=0;
    vector<double>* BDT_response_ob=0;
    vector<double>* BDT_response_on=0;

    vector<double>* BDT_response=0;

    Int_t Jet_mcFlavour[23];   //[nJet]

    TBranch *b_Jet_mcFlavour;

    vector<double>* jet_indices=0;

    double njet;

    datatree->SetBranchAddress("BDT_response",&BDT_response);
    datatree->SetBranchAddress("BDT_response_l",&BDT_response_l);
    datatree->SetBranchAddress("BDT_response_lb",&BDT_response_lb);
    datatree->SetBranchAddress("BDT_response_ln",&BDT_response_ln);
    datatree->SetBranchAddress("BDT_response_s",&BDT_response_s);
    datatree->SetBranchAddress("BDT_response_sb",&BDT_response_sb);
    datatree->SetBranchAddress("BDT_response_sn",&BDT_response_sn);
    datatree->SetBranchAddress("BDT_response_o",&BDT_response_o);
    datatree->SetBranchAddress("BDT_response_ob",&BDT_response_ob);
    datatree->SetBranchAddress("BDT_response_on",&BDT_response_on);

    datatree->SetBranchAddress("njet",&njet);
    datatree->SetBranchAddress("Jet_mcFlavour", Jet_mcFlavour, &b_Jet_mcFlavour);
    datatree->SetBranchAddress("jet_indices",&jet_indices);

    TH1F *splot_p_l_bdt = new TH1F("splot_p_l_bdt" ,"Leading jet BDT" ,300 , -1, 1);
    TH1F *splot_p_lb_bdt = new TH1F("splot_p_lb_bdt" ,"Leading jet BDT" ,300 , -1, 1);
    TH1F *splot_p_ln_bdt = new TH1F("splot_p_ln_bdt" ,"Leading jet BDT" ,300 , -1, 1);
    TH1F *splot_p_s_bdt = new TH1F("splot_p_s_bdt" ,"Sub-leading jet BDT" ,300 , -1, 1);
    TH1F *splot_p_sb_bdt = new TH1F("splot_p_sb_bdt" ,"Sub-leading jet BDT" ,300 , -1, 1);
    TH1F *splot_p_sn_bdt = new TH1F("splot_p_sn_bdt" ,"Sub-leading jet BDT" ,300 , -1, 1);
    TH1F *splot_p_o_bdt = new TH1F("splot_p_o_bdt" ,"Other jets BDT" ,300 , -1, 1);
    TH1F *splot_p_ob_bdt = new TH1F("splot_p_ob_bdt" ,"Other bjets BDT" ,300 , -1, 1);
    TH1F *splot_p_on_bdt = new TH1F("splot_p_on_bdt" ,"Other non-bjets BDT" ,300 , -1, 1);

    Long64_t nentries = datatree->GetEntries();

    for (Long64_t ientry=0; ientry < nentries; ientry++)
    {
        datatree->GetEntry(ientry);

        for (Int_t i = 0; i < njet;i++)
        {
       
            if ( i==0 )
            {
            splot_p_l_bdt->Fill((*BDT_response)[0]);

                int t = (*jet_indices)[0];

                if ((fabs(Jet_mcFlavour[t])) == 5)
                {
                    splot_p_lb_bdt->Fill((*BDT_response)[0]);
                }
                else
                {
                    splot_p_ln_bdt->Fill((*BDT_response)[0]);
                }
            }

            if (i==1)
            {
       
                splot_p_s_bdt->Fill((*BDT_response)[1]);

                int t = (*jet_indices)[1];
                if (((fabs(Jet_mcFlavour[t])) == 5))
                {
                    splot_p_sb_bdt->Fill((*BDT_response)[1]);
                }
                else
                {
                    splot_p_sn_bdt->Fill((*BDT_response)[1]);
                }

            }

            if (i>1)
            {

                splot_p_o_bdt->Fill((*BDT_response)[i]);

                int t = (*jet_indices)[i];
                if (((fabs(Jet_mcFlavour[t])) == 5))
                {
                    splot_p_ob_bdt->Fill((*BDT_response)[i]);
                }
                else
                {
                    splot_p_on_bdt->Fill((*BDT_response)[i]);
                }

            }
        }
    }

    splot_p_l_bdt->Write();
    splot_p_lb_bdt->Write();
    splot_p_ln_bdt->Write();
    splot_p_s_bdt->Write();
    splot_p_sb_bdt->Write();
    splot_p_sn_bdt->Write();
    splot_p_o_bdt->Write();
    splot_p_ob_bdt->Write();
    splot_p_on_bdt->Write();

    datatree->Write();
//    datatree->Print();


    TString blead("leading bjet"),f0("l");
    TString nlead("leading non-bjet");

    TString bsub("subleading bjet");
    TString nsub("subleading non-bjet");

    TString bother("other bjets");
    TString nother("other non-bjets");

    TString yaxis("#events");

    Double_t scale;

    TCanvas *splot_b_l_bdt = new TCanvas ("splot_b_l_bdt","splot_b_l_bdt");
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


}

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
#include "TLine.h"
#include "TImage.h"
#include <TPaveLabel.h>
#include <RooAbsDataStore.h>
#include <TImage.h>
#include <TAttText.h>
#include <TTF.h>
#include <THStack.h>
#include <TPave.h>
#include <TPaveText.h>
#include "TLatex.h"

using namespace std;

void plots(TString label, TString canvas, TH1* h1, TH1* h2);
void plot_dataMC(TString label, TString canvas, TH1* h1, TH1* h2);
void sweights(TH1* f1,TH1* g1,TString title,TString save);
void plot_rho(TString label, TString canvas, TH1* h1, TH1* h2, TH1* h3);
void sig_bg_plots()
{

    TFile* data = new TFile("TMVA_output.root");
    TTree* datatree = (TTree*)(data->Get("tree"));

    TFile* sw_data = new TFile("qreg_pt_norm.root");//has sweights
    TTree* sw_datatree = (TTree*)(sw_data->Get("dataset"));
    TTree* r_sw_datatree = (TTree*)(sw_data->Get("r_dataset"));

    TFile* bdt_data = new TFile("TMVA_output_splots.root");//will give you rho weight
    TTree* bdt_datatree = (TTree*)(bdt_data->Get("tree"));

    TFile* r_bdt_data = new TFile("r_TMVA_output_splots.root");//will give you rho weight
    TTree* r_bdt_datatree = (TTree*)(r_bdt_data->Get("tree"));

    TFile* r_data = new TFile("r_TMVA_output.root");
    TTree* r_datatree = (TTree*)(r_data->Get("tree"));

    double bdt_response;
    double r_bdt_response;
    double yield_b_sw;
    double yield_n_sw;
    double r_yield_b_sw;
    double r_yield_n_sw;
    double p_sw_mc_norm;
    double p_sw_data_norm;
    double rho_weight;
    double rho;
    double r_rho;
 
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
    vector<double>* BDT_response=0;
    double njet;


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

    //Int_t Jet_mcFlavour[23];   //[nJet]
    //TBranch *b_Jet_mcFlavour;

    vector<double>* r_jet_indices=0;
    vector<double>* r_BDT_response=0;
    double r_njet;

    datatree->SetBranchAddress("BDT_response",&BDT_response);
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

    r_datatree->SetBranchAddress("r_BDT_response",&r_BDT_response);
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

    r_datatree->SetBranchAddress("jet_indices",&r_jet_indices);

    sw_datatree->SetBranchAddress("BDT_response",&bdt_response);
    bdt_datatree->SetBranchAddress("rho",&rho);
    bdt_datatree->SetBranchAddress("rho_weight",&rho_weight);
    r_sw_datatree->SetBranchAddress("BDT_response",&r_bdt_response);
    sw_datatree->SetBranchAddress("yield_b_sw",&yield_b_sw);
    r_sw_datatree->SetBranchAddress("r_yield_b_sw",&r_yield_b_sw);
    sw_datatree->SetBranchAddress("yield_n_sw",&yield_n_sw);
    r_sw_datatree->SetBranchAddress("r_yield_n_sw",&r_yield_n_sw);
    r_bdt_datatree->SetBranchAddress("rho",&r_rho);

    TH1F *p_sig_m_nl_j = new TH1F("p_sig_m_nl_j", "m_{(nl,j)} GeV",100,0,600);
    TH1F *p_bg_m_nl_j = new TH1F("p_bg_m_nl_j", "m_{(nl,j)} GeV",100,0,600);
    TH1F *p_m_nl_j = new TH1F("p_m_nl_j", "m_{(nl,j)}",100,0,600);
    TH1F *r_p_m_nl_j = new TH1F("r_p_m_nl_j", "m_{(nl,j)}",100,0,600);

    TH1F *p_sig_m_fl_j = new TH1F("p_sig_m_fl_j", "m_{(fl,j)} GeV",100,0,600);
    TH1F *p_bg_m_fl_j = new TH1F("p_bg_m_fl_j", "m_{(fl,j)} GeV",100,0,600);
    TH1F *p_m_fl_j = new TH1F("p_m_fl_j", "m_{(fl,j)}",100,0,600);
    TH1F *r_p_m_fl_j = new TH1F("r_p_m_fl_j", "m_{(fl,j)}",100,0,600);

    TH1F *p_sig_eta_nl_j = new TH1F("p_sig_eta_nl_j", "#Delta#eta_{(nl,j)}",100,0,5);
    TH1F *p_bg_eta_nl_j = new TH1F("p_bg_eta_nl_j", "#Delta#eta_{(nl,j)}",100,0,5);
    TH1F *p_eta_nl_j = new TH1F("p_eta_nl_j", "#Delta#eta_{(nl,j)}",100,0,5);
    TH1F *r_p_eta_nl_j = new TH1F("r_p_eta_nl_j", "#Delta#eta_{(nl,j)}",100,0,5);

    TH1F *p_sig_eta_fl_j = new TH1F("p_sig_eta_fl_j", "#Delta#eta_{(fl,j)}",100,0,5);
    TH1F *p_bg_eta_fl_j = new TH1F("p_bg_eta_fl_j", "#Delta#eta_{(fl,j)}",100,0,5);
    TH1F *p_eta_fl_j = new TH1F("p_eta_fl_j", "#Delta#eta_{(fl,j)}",100,0,5);
    TH1F *r_p_eta_fl_j = new TH1F("r_p_eta_fl_j", "#Delta#eta_{(fl,j)}",100,0,5);

    TH1F *p_sig_eta_nlj_ll = new TH1F("p_sig_eta_nlj_ll", "#Delta#eta_{(nlj,ll)}",100,0,6);
    TH1F *p_bg_eta_nlj_ll = new TH1F("p_bg_eta_nlj_ll", "#Delta#eta_{(nlj,ll)}",100,0,6);
    TH1F *p_eta_nlj_ll = new TH1F("p_eta_nlj_ll", "#Delta#eta_{(nlj,ll)}",100,0,6);
    TH1F *r_p_eta_nlj_ll = new TH1F("r_p_s_eta_nlj_ll", "#Delta#eta_{(nlj,ll)}",100,0,6);

    TH1F *p_sig_eta_flj_ll = new TH1F("p_sig_eta_flj_ll", "#Delta#eta_{(flj,ll)}",100,0,6);
    TH1F *p_bg_eta_flj_ll = new TH1F("p_bg_eta_flj_ll", "#Delta#eta_{(flj,ll)}",100,0,6);
    TH1F *p_eta_flj_ll = new TH1F("p_eta_flj_ll", "#Delta#eta_{(flj,ll)}",100,0,6);
    TH1F *r_p_eta_flj_ll = new TH1F("r_p_eta_flj_ll", "#Delta#eta_{(flj,ll)}",100,0,6);

    TH1F *p_sig_eta_j_ll = new TH1F("p_sig_eta_j_ll", "#Delta#eta_{(j,ll)}",100,0,6);
    TH1F *p_bg_eta_j_ll = new TH1F("p_bg_eta_j_ll", "#Delta#eta_{(j,ll)}",100,0,6);
    TH1F *p_eta_j_ll = new TH1F("p_eta_j_ll", "#Delta#eta_{(j,ll)}",100,0,6);
    TH1F *r_p_eta_j_ll = new TH1F("r_p_eta_j_ll", "#Delta#eta_{(j,ll)}",100,0,6);

    TH1F *p_sig_phi_nl_j = new TH1F("p_sig_phi_nl_j", "#Delta#phi_{(nl,j)}",100,0,3.5);
    TH1F *p_bg_phi_nl_j = new TH1F("p_bg_phi_nl_j", "#Delta#phi_{(nl,j)}",100,0,3.5);
    TH1F *p_phi_nl_j = new TH1F("p_phi_nl_j", "#Delta#phi_{(nl,j)}",100,0,3.5);
    TH1F *r_p_phi_nl_j = new TH1F("r_p_phi_nl_j", "#Delta#phi_{(nl,j)}",100,0,3.5);

    TH1F *p_sig_phi_fl_j = new TH1F("p_sig_phi_fl_j", "#Delta#phi_{(fl,j)}",100,0,3.5);
    TH1F *p_bg_phi_fl_j = new TH1F("p_bg_phi_fl_j", "#Delta#phi_{(fl,j)}",100,0,3.5);
    TH1F *p_phi_fl_j = new TH1F("p_phi_fl_j", "#Delta#phi_{(fl,j)}",100,0,3.5);
    TH1F *r_p_phi_fl_j = new TH1F("r_p_phi_fl_j", "#Delta#phi_{(fl,j)}",100,0,3.5);

    TH1F *p_sig_phi_nlj_ll = new TH1F("p_sig_phi_nlj_ll", "#Delta#phi_{(nlj,ll)}",100,0,3.5);
    TH1F *p_bg_phi_nlj_ll = new TH1F("p_bg_phi_nlj_ll", "#Delta#phi_{(nlj,ll)}",100,0,3.5);
    TH1F *p_phi_nlj_ll = new TH1F("p_phi_nlj_ll", "#Delta#phi_{(nlj,ll)}",100,0,3.5);
    TH1F *r_p_phi_nlj_ll = new TH1F("r_p_phi_nlj_ll", "#Delta#phi_{(nlj,ll)}",100,0,3.5);

    TH1F *p_sig_phi_flj_ll = new TH1F("p_sig_phi_flj_ll", "#Delta#phi_{(flj,ll)}",100,0,3.5);
    TH1F *p_bg_phi_flj_ll = new TH1F("p_bg_phi_flj_ll", "#Delta#phi_{(flj,ll)}",100,0,3.5);
    TH1F *p_phi_flj_ll = new TH1F("p_phi_flj_ll", "#Delta#phi_{(flj,ll)}",100,0,3.5);
    TH1F *r_p_phi_flj_ll = new TH1F("r_p_phi_flj_ll", "#Delta#phi_{(flj,ll)}",100,0,3.5);

    TH1F *p_sig_phi_j_ll = new TH1F("p_sig_phi_j_ll", "#Delta#phi_{(j,ll)}",100,0,3.5);
    TH1F *p_bg_phi_j_ll = new TH1F("p_bg_phi_j_ll", "#Delta#phi_{(j,ll)}",100,0,3.5);
    TH1F *p_phi_j_ll = new TH1F("p_phi_j_ll", "#Delta#phi_{(j,ll)}",100,0,3.5);
    TH1F *r_p_phi_j_ll = new TH1F("r_p_phi_j_ll", "#Delta#phi_{(j,ll)}",100,0,3.5);

    TH1F *p_bdt = new TH1F("p_bdt" ,"jet BDT" ,100 , -1, 1);
    TH1F *r_p_bdt = new TH1F("r_p_bdt" ,"jet BDT" ,100 , -1, 1);
    TH1F *p_l_bdt = new TH1F("p_l_bdt" ,"Leading jet BDT" ,100 , -1, 1);
    TH1F *r_p_l_bdt = new TH1F("r_p_l_bdt" ,"Leading jet BDT" ,100 , -1, 1);
    TH1F *p_lb_bdt = new TH1F("p_lb_bdt" ,"Leading jet BDT" ,100 , -1, 1);
    TH1F *p_ln_bdt = new TH1F("p_ln_bdt" ,"Leading jet BDT" ,100 , -1, 1);
    TH1F *p_s_bdt = new TH1F("p_s_bdt" ,"Sub-leading jet BDT" ,100 , -1, 1);
    TH1F *r_p_s_bdt = new TH1F("r_p_s_bdt" ,"Sub-leading jet BDT" ,100 , -1, 1);
    TH1F *p_sb_bdt = new TH1F("p_sb_bdt" ,"Sub-leading jet BDT" ,100 , -1, 1);
    TH1F *p_sn_bdt = new TH1F("p_sn_bdt" ,"Sub-leading jet BDT" ,100 , -1, 1);
    TH1F *p_o_bdt = new TH1F("p_o_bdt" ,"Other jets BDT" ,100 , -1, 1);
    TH1F *r_p_o_bdt = new TH1F("r_p_o_bdt" ,"Other jets BDT" ,100 , -1, 1);
    TH1F *p_ob_bdt = new TH1F("p_ob_bdt" ,"Other bjets BDT" ,100 , -1, 1);
    TH1F *p_on_bdt = new TH1F("p_on_bdt" ,"Other non-bjets BDT" ,100 , -1, 1);

    TH1F *bdt_rho_yield_b = new TH1F("bdt_rho_yield_b" ,"rho_bdt_yield_b" ,100 , -1, 1);
    TH1F *bdt_rho_yield_n = new TH1F("bdt_rho_yield_n" ,"rho_bdt_yield_n" ,100 , -1, 1);
    TH1F *bdt_rho_r_yield_b = new TH1F("bdt_rho_r_yield_b" ,"rho_bdt_r_yield_b" ,100 , -1, 1);
    TH1F *bdt_rho_r_yield_n = new TH1F("bdt_rho_r_yield_n" ,"rho_bdt_r_yield_n" ,100 , -1, 1);
    TH1F *n_yield_b_sw = new TH1F("n_yield_b_sw" ,"yield_b_sw" ,100 , -5, 5);
    TH1F *n_yield_n_sw = new TH1F("n_yield_n_sw" ,"yield_n_sw" ,100 , -5, 5); 
    TH1F *n_r_yield_b_sw = new TH1F("n_r_yield_b_sw" ,"r_yield_b_sw" ,100 , -5, 5);
    TH1F *n_r_yield_n_sw = new TH1F("n_r_yield_n_sw" ,"r_yield_n_sw" ,100 , -5, 5);


    TH1F *p_rho_unwgt = new TH1F("p_rho_unwgt" ,"p_rho_unwgt" ,100 , 0, 80);
    TH1F *p_rho_wgt = new TH1F("p_rho_wgt" ,"p_rho_wgt" ,100 , 0, 80);
    TH1F *p_r_rho = new TH1F("p_r_rho" ,"p_r_rho" ,100 , 0, 80);

/*
    Long64_t nentries = datatree->GetEntries();
    cout<<nentries<<"is the number of entries"<<endl;

    for (Long64_t ientry=0; ientry < nentries; ientry++)
    {
        datatree->GetEntry(ientry);

        for (Int_t i = 0; i < njet;i++)
        {
            int t = (*jet_indices)[i];
            
            p_m_nl_j->Fill((*m_nl_j)[i]);
            p_m_fl_j->Fill((*m_fl_j)[i]);
            p_eta_nl_j->Fill(abs((*eta_nl_j)[i]));
            p_eta_fl_j->Fill(abs((*eta_fl_j)[i]));
            p_eta_nlj_ll->Fill(abs((*eta_nlj_ll)[i]));
            p_eta_flj_ll->Fill(abs((*eta_flj_ll)[i]));
            p_eta_j_ll->Fill(abs((*eta_j_ll)[i]));
            p_phi_nl_j->Fill(abs((*phi_nl_j)[i]));
            p_phi_fl_j->Fill(abs((*phi_fl_j)[i]));
            p_phi_nlj_ll->Fill(abs((*phi_nlj_ll)[i]));
            p_phi_flj_ll->Fill(abs((*phi_flj_ll)[i]));
            p_phi_j_ll->Fill(abs((*phi_j_ll)[i]));
            
            if ((abs(Jet_mcFlavour[t])) == 5)
            {
                p_sig_m_nl_j->Fill((*m_nl_j)[i]);
	        p_sig_m_fl_j->Fill((*m_fl_j)[i]);
	        p_sig_eta_nl_j->Fill(abs((*eta_nl_j)[i]));
	        p_sig_eta_fl_j->Fill(abs((*eta_fl_j)[i]));
	        p_sig_eta_nlj_ll->Fill(abs((*eta_nlj_ll)[i]));
	        p_sig_eta_flj_ll->Fill(abs((*eta_flj_ll)[i]));
	        p_sig_eta_j_ll->Fill(abs((*eta_j_ll)[i]));
	        p_sig_phi_nl_j->Fill(abs((*phi_nl_j)[i]));
	        p_sig_phi_fl_j->Fill(abs((*phi_fl_j)[i]));
	        p_sig_phi_nlj_ll->Fill(abs((*phi_nlj_ll)[i]));
	        p_sig_phi_flj_ll->Fill(abs((*phi_flj_ll)[i]));
	        p_sig_phi_j_ll->Fill(abs((*phi_j_ll)[i]));
            }            
            
            else
            {
                p_bg_m_nl_j->Fill(abs((*m_nl_j)[i]));
                p_bg_m_fl_j->Fill(abs((*m_fl_j)[i]));
                p_bg_eta_nl_j->Fill(abs((*eta_nl_j)[i]));
                p_bg_eta_fl_j->Fill(abs((*eta_fl_j)[i]));
                p_bg_eta_nlj_ll->Fill(abs((*eta_nlj_ll)[i]));
                p_bg_eta_flj_ll->Fill(abs((*eta_flj_ll)[i]));
                p_bg_eta_j_ll->Fill(abs((*eta_j_ll)[i]));
                p_bg_phi_nl_j->Fill(abs((*phi_nl_j)[i]));
                p_bg_phi_fl_j->Fill(abs((*phi_fl_j)[i]));
                p_bg_phi_nlj_ll->Fill(abs((*phi_nlj_ll)[i]));
                p_bg_phi_flj_ll->Fill(abs((*phi_flj_ll)[i]));
                p_bg_phi_j_ll->Fill(abs((*phi_j_ll)[i]));
            }           

            p_bdt->Fill((*BDT_response)[i]);
            if ( i==0 )
            {
                p_l_bdt->Fill((*BDT_response)[0]);
                //int t = (*jet_indices)[0];
                if ((fabs(Jet_mcFlavour[t])) == 5)
                {
                    p_lb_bdt->Fill((*BDT_response)[0]);
                }
                else
                {
                    p_ln_bdt->Fill((*BDT_response)[0]);
                }
            }

            if (i==1)
            {

                p_s_bdt->Fill((*BDT_response)[1]);

                //int t = (*jet_indices)[1];
                if (((fabs(Jet_mcFlavour[t])) == 5))
                {
                    p_sb_bdt->Fill((*BDT_response)[1]);
                }
                else
                {
                    p_sn_bdt->Fill((*BDT_response)[1]);
                }

            }

            if (i>1)
            {

                p_o_bdt->Fill((*BDT_response)[i]);

                //int t = (*jet_indices)[i];
                if (((fabs(Jet_mcFlavour[t])) == 5))
                {
                    p_ob_bdt->Fill((*BDT_response)[i]);
                }
                else
                {
                    p_on_bdt->Fill((*BDT_response)[i]);
                }

            }

        }

    }


    Long64_t r_nentries = r_datatree->GetEntries();
    cout<<r_nentries<<"is the number of entries in data"<<endl;

    for (Long64_t r_ientry=0; r_ientry < r_nentries; r_ientry++)
    {
        r_datatree->GetEntry(r_ientry);

        for (Int_t i = 0; i < r_njet;i++)
        {

            r_p_m_nl_j->Fill((*r_m_nl_j)[i]);
            r_p_m_fl_j->Fill((*r_m_fl_j)[i]);
            r_p_eta_nl_j->Fill(abs((*r_eta_nl_j)[i]));
            r_p_eta_fl_j->Fill(abs((*r_eta_fl_j)[i]));
            r_p_eta_nlj_ll->Fill(abs((*r_eta_nlj_ll)[i]));
            r_p_eta_flj_ll->Fill(abs((*r_eta_flj_ll)[i]));
            r_p_eta_j_ll->Fill(abs((*r_eta_j_ll)[i]));
            r_p_phi_nl_j->Fill(abs((*r_phi_nl_j)[i]));
            r_p_phi_fl_j->Fill(abs((*r_phi_fl_j)[i]));
            r_p_phi_nlj_ll->Fill(abs((*r_phi_nlj_ll)[i]));
            r_p_phi_flj_ll->Fill(abs((*r_phi_flj_ll)[i]));
            r_p_phi_j_ll->Fill(abs((*r_phi_j_ll)[i]));

            r_p_bdt->Fill((*r_BDT_response)[i]);
            if ( i==0 )
            {
                r_p_l_bdt->Fill((*r_BDT_response)[0]);
            }

            if (i==1)
            {
                r_p_s_bdt->Fill((*r_BDT_response)[1]);
            }

            if (i>1)
            {
                r_p_o_bdt->Fill((*r_BDT_response)[i]);
            }

        }
    }

*/

    Long64_t sw_nentries = sw_datatree->GetEntries();
    cout<<sw_nentries<<"is the number of entries in sw_mc"<<endl;

    for (Long64_t sw_ientry=0; sw_ientry < sw_nentries; sw_ientry++)
    {
        sw_datatree->GetEntry(sw_ientry);
        bdt_datatree->GetEntry(sw_ientry);
  
        n_yield_b_sw->Fill(yield_b_sw);
        n_yield_n_sw->Fill(yield_n_sw);   

        if(yield_b_sw>0)
        {
            bdt_rho_yield_b->Fill(bdt_response);
        }
        if(yield_n_sw>0)
        {
            bdt_rho_yield_n->Fill(bdt_response);
        }

        if(rho>0)
        {
            p_rho_unwgt->Fill(rho);
            p_rho_wgt->Fill(rho,rho_weight);
        }
    }  

    Long64_t r_sw_nentries = r_sw_datatree->GetEntries();
    cout<<r_sw_nentries<<"is the number of entries in sw_mc"<<endl;

    for (Long64_t r_sw_ientry=0; r_sw_ientry < r_sw_nentries; r_sw_ientry++)
    {
        r_sw_datatree->GetEntry(r_sw_ientry);
        r_bdt_datatree->GetEntry(r_sw_ientry);

        n_r_yield_b_sw->Fill(r_yield_b_sw);
        n_r_yield_n_sw->Fill(r_yield_n_sw);

        if(r_yield_b_sw>0)
        {
             bdt_rho_r_yield_b->Fill(r_bdt_response);
        }
        if(r_yield_n_sw>0)
        {
             bdt_rho_r_yield_n->Fill(r_bdt_response);
        }

        if(r_rho>0)
        {
            p_r_rho->Fill(r_rho);
        }
    }

//    plot_rho(TString("rho"), TString("rho"), p_r_rho, p_rho_unwgt, p_rho_wgt);


    sweights(n_yield_b_sw,n_yield_n_sw,"#scale[4]{#font[12]{sWeights for MC}}","MC_sweights");
    sweights(n_r_yield_b_sw,n_r_yield_n_sw,"#scale[4]{#font[12]{sWeights for data}}","data_sweights");

    sweights(bdt_rho_yield_b,bdt_rho_yield_n,"#scale[4]{#font[12]{BDT score for +ve sWeighted MC}}","bdt_MC_sweights");
    sweights(bdt_rho_r_yield_b,bdt_rho_r_yield_n,"#scale[4]{#font[12]{BDT score for +ve sWeighted data}}","bdt_data_sweights");

/*  
    
    plots("m_{(nl,j)} GeV","m_nl_j",p_sig_m_nl_j,p_bg_m_nl_j);
    plots("m_{(fl,j)} GeV","m_fl_j",p_sig_m_fl_j,p_bg_m_fl_j);
    
    plots("#Delta#eta_{(nl,j)}","eta_nl_j",p_sig_eta_nl_j,p_bg_eta_nl_j);
    plots("#Delta#eta_{(fl,j)}","eta_fl_j",p_sig_eta_fl_j,p_bg_eta_fl_j);

    plots("#Delta#eta_{(nlj,ll)}","eta_nlj_ll",p_sig_eta_nlj_ll,p_bg_eta_nlj_ll);
    plots("#Delta#eta_{(flj,ll)}","eta_flj_ll",p_sig_eta_flj_ll,p_bg_eta_flj_ll);

    plots("#Delta#eta_{(j,ll)}","eta_j_ll",p_sig_eta_j_ll,p_bg_eta_j_ll);

    plots("#Delta#phi_{(nl,j)} rad","phi_nl_j",p_sig_phi_nl_j,p_bg_phi_nl_j);
    plots("#Delta#phi_{(fl,j)} rad","phi_fl_j",p_sig_phi_fl_j,p_bg_phi_fl_j);

    plots("#Delta#phi_{(nlj,ll)} rad","phi_nlj_ll",p_sig_phi_nlj_ll,p_bg_phi_nlj_ll);
    plots("#Delta#phi_{(flj,ll)} rad","phi_flj_ll",p_sig_phi_flj_ll,p_bg_phi_flj_ll);

    plots("#Delta#phi_{(j,ll)} rad","phi_j_ll",p_sig_phi_j_ll,p_bg_phi_j_ll);
    
    plots("Leading jet BDT response","p_l_bdt",p_lb_bdt,p_ln_bdt);
    plots("Sub-leading jet BDT response","p_s_bdt",p_sb_bdt,p_sn_bdt);
    plots("Other jet BDT response","p_o_bdt",p_ob_bdt,p_on_bdt);
    
    plot_dataMC("m_{(nl,j)} GeV","r_m_nl_j",r_p_m_nl_j,p_m_nl_j);
    plot_dataMC("m_{(fl,j)} GeV","r_m_fl_j",r_p_m_fl_j,p_m_fl_j);

    plot_dataMC("#Delta#eta_{(nl,j)}","r_eta_nl_j",r_p_eta_nl_j,p_eta_nl_j);
    plot_dataMC("#Delta#eta_{(fl,j)}","r_eta_fl_j",r_p_eta_fl_j,p_eta_fl_j);

    plot_dataMC("#Delta#eta_{(nlj,ll)}","r_eta_nlj_ll",r_p_eta_nlj_ll,p_eta_nlj_ll);
    plot_dataMC("#Delta#eta_{(flj,ll)}","r_eta_flj_ll",r_p_eta_flj_ll,p_eta_flj_ll);

    plot_dataMC("#Delta#eta_{(j,ll)}","r_eta_j_ll",r_p_eta_j_ll,p_eta_j_ll);

    plot_dataMC("#Delta#phi_{(nl,j)} rad","r_phi_nl_j",r_p_phi_nl_j,p_phi_nl_j);
    plot_dataMC("#Delta#phi_{(fl,j)} rad","r_phi_fl_j",r_p_phi_fl_j,p_phi_fl_j);

    plot_dataMC("#Delta#phi_{(nlj,ll)} rad","r_phi_nlj_ll",r_p_phi_nlj_ll,p_phi_nlj_ll);
    plot_dataMC("#Delta#phi_{(flj,ll)} rad","r_phi_flj_ll",r_p_phi_flj_ll,p_phi_flj_ll);

    plot_dataMC("#Delta#phi_{(j,ll)} rad","r_phi_j_ll",r_p_phi_j_ll,p_phi_j_ll);
    
    plot_dataMC("Leading jet BDT response","r_p_l_bdt",r_p_l_bdt,p_l_bdt);
    plot_dataMC("Sub-leading jet BDT response","r_p_s_bdt",r_p_s_bdt,p_s_bdt);
    plot_dataMC("Other jet BDT response","r_p_o_bdt",r_p_o_bdt,p_o_bdt);
    plot_dataMC("BDT response","r_p_bdt",r_p_bdt,p_bdt);
*/   
}
void plots(TString label, TString canvas, TH1* h1, TH1* h2)
    
{
    TH1F *g1 = (TH1F*)h1->Clone("g1");
    TH1F *g2 = (TH1F*)h2->Clone("g2");
    TCanvas* c1 = new TCanvas(canvas,canvas,65,52,560,403);
    c1->Range(0,0,1,1);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetFrameBorderMode(0);
    gStyle->SetPadBorderMode(0);
    c1->Divide(1,2);
    c1->cd(1);
    gPad->SetPad(0,0.33,1,0.9378806);
    gPad->Range(-69.03766,-0.0002509835,628.8703,0.06681207);
    gPad->SetGridx();
    gPad->SetTopMargin(0.006574267);
    gPad->SetBottomMargin(0.001912046);

    h1->Scale(1/(h1->GetSumOfWeights()));
    h2->Scale(1/(h2->GetSumOfWeights()));
    h1->SetLineColor(kRed);
    h2->SetLineColor(kBlue);
    THStack* hs = new THStack();
    hs->Add(h1,"HIST");
    hs->Add(h2,"HIST");
    hs->Draw("nostack");
    
    hs->GetHistogram()->GetXaxis()->SetLabelOffset(999);
    hs->GetHistogram()->GetXaxis()->SetTitle("");
    hs->GetHistogram()->GetYaxis()->SetTitle("a.u");
    hs->GetHistogram()->GetYaxis()->SetTitleOffset(0.85);
    hs->GetHistogram()->GetYaxis()->SetTitleSize(0.058);
    hs->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
    hs->GetHistogram()->GetYaxis()->CenterTitle();
    TLegend* l1a;
    if ((canvas == "phi_fl_j") || (canvas == "phi_j_ll"))
    {
      TLegend* la = new TLegend(0.1888489,0.7661965,0.4046763,0.9887377,NULL,"brNDC");
      l1a = la;
    }
    else
    {
      TLegend* la = new TLegend(0.6816547,0.7688758,0.897482,0.9892674,NULL,"brNDC");
      l1a = la;
    }
    l1a->AddEntry(h1,"signal","l");
    l1a->AddEntry(h2,"background","l");
    l1a->SetTextFont(12);
    l1a->SetTextSize(0.06237499);
    l1a->SetBorderSize(0);
    l1a->SetFillColor(0);
    l1a->Draw("same");
    
    
    c1->cd(2);
    gPad->SetPad(0.0,0.0,1.0,0.33);
    gPad->SetTickx();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTopMargin(0);
    gPad->SetBottomMargin(0.3536693);
    g2->Scale(g1->GetSumOfWeights()/(g2->GetSumOfWeights()));
    g1->Sumw2();
    g1->Divide(g2);     
    g1->SetStats(0);
    g1->SetLineColor(kBlack);
    g1->SetTitle("");
    g1->GetXaxis()->SetLabelSize(0.12);
    g1->GetXaxis()->SetTitle(label);
    g1->GetXaxis()->SetTitleSize(0.15);
    g1->GetXaxis()->SetTitleOffset(0.99);
    //g1->GetXaxis()->CenterTitle();
    g1->GetXaxis()->SetTickLength(0.03);
    g1->GetYaxis()->SetTitle("sig/bg");
    g1->GetYaxis()->CenterTitle();
    g1->GetYaxis()->SetTitleSize(0.115);
    g1->GetYaxis()->SetLabelSize(0.105);
    g1->GetYaxis()->SetTitleOffset(0.35);
    g1->GetYaxis()->SetNdivisions(406);
    g1->GetYaxis()->SetRangeUser(0.5,1.5);
    g1->Draw("ep");
    Double_t low = g1->GetXaxis()->GetXmin();
    Double_t high = g1->GetXaxis()->GetXmax();
    TLine *line = new TLine(low, 1, high, 1);
    line->SetLineColor(2);
    line->SetLineStyle(7);
    line->Draw(); 
    gPad->Modified();
    c1->Update();
    
    gSystem->ProcessEvents();

    TString addr = "/mnt/t3nfs01/data01/shome/krgedia/Workspace/lxplus/qreg_pt_norm/dataMC_pt/sig_bg_KIN/";

    TImage *img = TImage::Create();
    img->FromPad(c1);
    //img->WriteImage((addr+label+".png").c_str());
    img->WriteImage(addr+canvas+".png");
    
}


void plot_dataMC(TString label, TString canvas, TH1* h1, TH1* h2)
{
    TH1F *g1 = (TH1F*)h1->Clone("g1");
    TH1F *g2 = (TH1F*)h2->Clone("g2");
    TCanvas* c1 = new TCanvas(canvas,canvas,65,52,560,403);
    c1->Range(0,0,1,1);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetFrameBorderMode(0);
    gStyle->SetPadBorderMode(0);
    c1->Divide(1,2);
    c1->cd(1);
    gPad->SetPad(0,0.33,1,0.9378806);
    gPad->Range(-69.03766,-0.0002509835,628.8703,0.06681207);
    gPad->SetGridx();
    gPad->SetTopMargin(0.006574267);
    gPad->SetBottomMargin(0.001912046);
    h2->Scale((h1->GetSumOfWeights())/(h2->GetSumOfWeights()));
    h1->SetLineColor(kBlack);
    h2->SetLineColor(kRed);
    THStack* hs = new THStack();
    hs->Add(h1,"E0");
    hs->Add(h2,"HIST");
    hs->Draw("nostack");

    hs->GetHistogram()->GetXaxis()->SetLabelOffset(999);
    hs->GetHistogram()->GetXaxis()->SetTitle("");
    double bin_width = h1->GetXaxis()->GetBinWidth(0);
    //string str = to_string(bin_width);
    //std::string str = boost::lexical_cast<std::string>(bin_width);
    std::ostringstream strs;
    strs << bin_width;
    std::string str = strs.str();
    str.erase (str.find_last_not_of('0') + 1, std::string::npos);
    TString yaxis = "jets/" + str;
    hs->GetHistogram()->GetYaxis()->SetTitle(yaxis);
    hs->GetHistogram()->GetYaxis()->SetTitleOffset(0.85);
    hs->GetHistogram()->GetYaxis()->SetTitleSize(0.058);
    hs->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
    hs->GetHistogram()->GetYaxis()->CenterTitle();
    TLegend* l1a;
    if ((canvas == "phi_fl_j") || (canvas == "phi_j_ll"))
    {
      TLegend* la = new TLegend(0.1888489,0.7661965,0.4046763,0.9887377,NULL,"brNDC");
      l1a = la;
    }
    else
    {
      TLegend* la = new TLegend(0.6816547,0.7688758,0.897482,0.9892674,NULL,"brNDC");
      l1a = la;
    }
    l1a->AddEntry(h1,"data","LE");
    l1a->AddEntry(h2,"MC Truth","l");

    l1a->SetTextFont(12);
    l1a->SetTextSize(0.06237499);
    l1a->SetBorderSize(0);
    l1a->SetFillColor(0);
    l1a->Draw("same");


    c1->cd(2);
    gPad->SetPad(0.0,0.0,1.0,0.33);
    gPad->SetTickx();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTopMargin(0);
    gPad->SetBottomMargin(0.3536693);
    g2->Scale((g1->GetSumOfWeights())/(g2->GetSumOfWeights()));
    g1->Sumw2();
    g1->Divide(g2);
    g1->SetStats(0);
    g1->SetLineColor(kBlack);
    g1->SetTitle("");
    g1->GetXaxis()->SetLabelSize(0.12);
    g1->GetXaxis()->SetTitle(label);
    g1->GetXaxis()->SetTitleSize(0.15);
    g1->GetXaxis()->SetTitleOffset(0.99);
    g1->GetXaxis()->SetTickLength(0.03);
    g1->GetYaxis()->SetTitle("data/MC");
    g1->GetYaxis()->CenterTitle();
    g1->GetYaxis()->SetTitleSize(0.115);
    g1->GetYaxis()->SetLabelSize(0.105);
    g1->GetYaxis()->SetTitleOffset(0.35);
    g1->GetYaxis()->SetNdivisions(406);
    g1->GetYaxis()->SetRangeUser(0.5,1.5);
    g1->Draw();
    Double_t low = g1->GetXaxis()->GetXmin();
    Double_t high = g1->GetXaxis()->GetXmax();
    TLine *line = new TLine(low, 1, high, 1);
    line->SetLineColor(2);
    line->SetLineStyle(7);
    line->Draw();
    gPad->Modified();
    c1->Update();

    
    gSystem->ProcessEvents();

    TString addr = "/mnt/t3nfs01/data01/shome/krgedia/Workspace/lxplus/qreg_pt_norm/dataMC_pt/data_MC_KIN/";

    TImage *img = TImage::Create();
    img->FromPad(c1);
    //img->WriteImage((addr+label+".png").c_str());
    img->WriteImage(addr+canvas+".png");
    
}


void sweights(TH1* f1,TH1* g1,TString title,TString save)
{

    TCanvas* canvas3 = new TCanvas(title,title,121,146,641,493);
    canvas3->Divide(1,2);
    canvas3->cd(1);
    gPad->SetPad(0,0.4811206,1,0.9317905);
    gPad->SetGridx();
    gPad->SetBottomMargin(0);
    double bin_width = f1->GetXaxis()->GetBinWidth(0);
    std::ostringstream strs;
    strs << bin_width;
    std::string str = strs.str();
    str.erase (str.find_last_not_of('0') + 1, std::string::npos);
    TString yaxis = "jets/" + str;
    f1->SetLineColor(kGreen+3);
    f1->SetStats(0);
    f1->SetTitle("");
    f1->GetXaxis()->SetLabelOffset(999);
    f1->GetYaxis()->SetTitle(yaxis);
    f1->GetYaxis()->SetTitleOffset(0.7);
    f1->GetYaxis()->SetTitleSize(0.064);
    f1->GetYaxis()->SetLabelSize(0.062);
    f1->Draw("HIST");
    gPad->Update();
    TLatex *text1 = new TLatex(0.5,0.5,"#scale[3.5]{#color[2]{#font[12]{Signal}}}");
    text1->SetTextSize(0.03375527);
    text1->Draw();
    TLatex *text2 = new TLatex(0.0,0.0,title);
    text2->SetTextSize(0.03375527);
    text2->Draw();
    Double_t high = gPad->GetUymax();
    TLine *line = new TLine(0, 0, 0, high);
    line->SetLineColor(2);
    line->SetLineStyle(1);
    line->Draw();


    canvas3->cd(2);
    gPad->SetPad(0.001113586,0.002436054,0.9988864,0.4835566);
    gPad->Range(-6.236102,-110.8116,6.21657,512.4);
    gPad->SetGridx();
    gPad->SetTopMargin(0);
    gPad->SetBottomMargin(0.1778074);
    g1->SetLineColor(kBlue);
    g1->SetStats(0);
    g1->SetTitle("");
    if (save == "MC_sweights" || save == "data_sweights")
    {
        g1->GetXaxis()->SetTitle("sWeights");
    }
    else
    {
        g1->GetXaxis()->SetTitle("BDT response");
    }
    g1->GetXaxis()->SetTitleSize(0.071);
    g1->GetXaxis()->SetTitleOffset(0.96);
    g1->GetXaxis()->SetLabelSize(0.070);
    g1->GetYaxis()->SetTitle(yaxis);
    g1->GetYaxis()->SetTitleOffset(0.75);
    g1->GetYaxis()->SetTitleSize(0.062);
    g1->GetYaxis()->SetLabelSize(0.062);
    g1->Draw("HIST");
    gPad->Update();
    TLatex *text3 = new TLatex(0.5,0.5,"#scale[3.5]{#color[4]{#font[12]{Background}}}");
    text3->SetTextSize(0.03375527);
    text3->Draw();
    
    Double_t ghigh = gPad->GetUymax();
    TLine *gline = new TLine(0, 0, 0, ghigh);
    gline->SetLineColor(2);
    gline->SetLineStyle(1);
    gline->Draw();
   
    canvas3->Update();
/*
    TString addr = "/mnt/t3nfs01/data01/shome/krgedia/Workspace/lxplus/qreg_pt_norm/dataMC_pt/sweights/";

    gSystem->ProcessEvents();

    TImage *img = TImage::Create();
    img->FromPad(canvas3);
    img->WriteImage(addr+save+".png");
*/
}

void plot_rho(TString label, TString canvas, TH1* h1, TH1* h2, TH1* h3)
{
    TH1F *g1 = (TH1F*)h1->Clone("g1");
    TH1F *g1c = (TH1F*)h1->Clone("g1c");
    TH1F *g2 = (TH1F*)h2->Clone("g2");
    TH1F* g3 = (TH1F*)h3->Clone("g3");
    TCanvas* c1 = new TCanvas(canvas,canvas,65,52,560,403);
    c1->Range(0,0,1,1);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetFrameBorderMode(0);
    gStyle->SetPadBorderMode(0);
    c1->Divide(1,2);
    c1->cd(1);
    gPad->SetPad(0,0.33,1,0.9378806);
    gPad->Range(-69.03766,-0.0002509835,628.8703,0.06681207);
    gPad->SetGridx();
    gPad->SetTopMargin(0.006574267);
    gPad->SetBottomMargin(0.001912046);
    h2->Scale((h1->GetSumOfWeights())/(h2->GetSumOfWeights()));
    h3->Scale((h1->GetSumOfWeights())/(h3->GetSumOfWeights()));
    h1->SetLineColor(kBlack);
    h2->SetLineColor(kGreen+2);
    h3->SetLineColor(kRed);
    THStack* hs = new THStack();
    hs->Add(h1,"E0");
    hs->Add(h2,"HIST");
    hs->Add(h3,"HIST");
    hs->Draw("nostack");

    hs->GetHistogram()->GetXaxis()->SetLabelOffset(999);
    hs->GetHistogram()->GetXaxis()->SetTitle("");
    double bin_width = h1->GetXaxis()->GetBinWidth(0);
    std::ostringstream strs;
    strs << bin_width;
    std::string str = strs.str();
    str.erase (str.find_last_not_of('0') + 1, std::string::npos);
    TString yaxis = "events/" + str;
    hs->GetHistogram()->GetYaxis()->SetTitle(yaxis);
    hs->GetHistogram()->GetYaxis()->SetTitleOffset(0.85);
    hs->GetHistogram()->GetYaxis()->SetTitleSize(0.058);
    hs->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
    //hs->GetHistogram()->GetYaxis()->CenterTitle();
    TLegend* l1a = new TLegend(0.6816547,0.7688758,0.897482,0.9892674,NULL,"brNDC");
    l1a->AddEntry(h1,"data","LE");
    l1a->AddEntry(h2,"MC Truth (unweighted)","l");
    l1a->AddEntry(h3,"MC Truth (rho_reweighted)","l");
    l1a->SetTextFont(12);
    l1a->SetTextSize(0.06237499);
    l1a->SetBorderSize(0);
    l1a->SetFillColor(0);
    l1a->Draw("same");


    c1->cd(2);
    gPad->SetPad(0.0,0.0,1.0,0.33);
    gPad->SetTickx();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTopMargin(0);
    gPad->SetBottomMargin(0.3536693);
    g2->Scale((g1->GetSumOfWeights())/(g2->GetSumOfWeights()));
    g3->Scale((g1c->GetSumOfWeights())/(g3->GetSumOfWeights()));     
    g1->Sumw2();
    g1c->Sumw2();
    g1->Divide(g2);
    g1c->Divide(g3);
    g1->SetStats(0);
    g1->SetLineColor(kGreen+2);
    g1c->SetLineColor(kRed);
    g1->SetTitle("");
    g1->GetXaxis()->SetLabelSize(0.12);
    g1->GetXaxis()->SetTitle(label);
    g1->GetXaxis()->SetTitleSize(0.15);
    g1->GetXaxis()->SetTitleOffset(0.99);
    g1->GetXaxis()->SetTickLength(0.03);
    g1->GetYaxis()->SetTitle("data/MC");
    g1->GetYaxis()->CenterTitle();
    g1->GetYaxis()->SetTitleSize(0.115);
    g1->GetYaxis()->SetLabelSize(0.105);
    g1->GetYaxis()->SetTitleOffset(0.35);
    g1->GetYaxis()->SetNdivisions(406);
    g1->GetYaxis()->SetRangeUser(0.5,1.5);
    g1->Draw();
    g1c->Draw("same");
    Double_t low = g1->GetXaxis()->GetXmin();
    Double_t high = g1->GetXaxis()->GetXmax();
    TLine *line = new TLine(low, 1, high, 1);
    line->SetLineColor(1);
    line->SetLineStyle(7);
    line->Draw();
    gPad->Modified();
    c1->Update();


    gSystem->ProcessEvents();

    //TString addr = "/mnt/t3nfs01/data01/shome/krgedia/Workspace/lxplus/qreg_pt_norm/dataMC_pt/data_MC_KIN/";

    TImage *img = TImage::Create();
    img->FromPad(c1);
    //img->WriteImage((addr+label+".png").c_str());
    //img->WriteImage(addr+canvas+".png");
    
    }
       































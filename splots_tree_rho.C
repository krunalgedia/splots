#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include <vector>
#include "RooPolynomial.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooFormulaVar.h"
#include "RooDecay.h"
#include "RooGaussModel.h"
#include "RooDataHist.h"
#include "RooProdPdf.h"
#include "RooHistPdf.h"

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
#include <TF1.h>
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
#include <TPaveLabel.h>

// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;
using namespace std;

void splots_tree_rho()
{
    TFile* r_file1 = new TFile("r_TMVA_output_splots.root");
    TTree* r_datatree = (TTree*)(r_file1->Get("tree"));

    TFile* file1 = new TFile("TMVA_output_splots.root","UPDATE");
    TTree* datatree = (TTree*)(file1->Get("tree"));

    double irho;
    double iJet_mass;
    double iJet_chHEF;
    datatree->SetBranchAddress("rho", &irho);
    datatree->SetBranchAddress("Jet_mass", &iJet_mass);
    datatree->SetBranchAddress("Jet_chHEF", &iJet_chHEF);
/*
    TFile* r_file1 = new TFile("r_TMVA_output_splots.root");
    TTree* r_datatree = (TTree*)(r_file1->Get("tree"));
*/
    RooRealVar rho("rho","rho",-200,100);
    RooRealVar Jet_mass("Jet_mass","Jet_mass",30000);
    RooRealVar Jet_chHEF("Jet_chHEF","Jet_chHEF",15);

    RooDataSet data_set("data_set","data_set",RooArgSet(rho,Jet_mass,Jet_chHEF),Import(*datatree));
    RooDataSet r_data_set("r_data_set","r_data_set",RooArgSet(rho,Jet_mass,Jet_chHEF),Import(*r_datatree));

    TH1* data_p_rho = (TH1*) data_set.createHistogram("data_p_rho",rho,Binning(100,0,100));
    TH1* r_data_p_rho = (TH1*) r_data_set.createHistogram("r_data_p_rho",rho,Binning(100,0,100));

    TH1* data_p_Jet_mass = (TH1*) data_set.createHistogram("data_p_Jet_mass",Jet_mass,Binning(100,0,60));
    TH1* r_data_p_Jet_mass = (TH1*) r_data_set.createHistogram("r_data_p_Jet_mass",Jet_mass,Binning(100,0,60));

    TH1* data_p_Jet_chHEF = (TH1*) data_set.createHistogram("data_p_Jet_chHEF",Jet_chHEF,Binning(100,0,1.2));
    TH1* r_data_p_Jet_chHEF = (TH1*) r_data_set.createHistogram("r_data_p_Jet_chHEF",Jet_chHEF,Binning(100,0,1.2));

    data_p_rho->Scale(1/data_p_rho->Integral());
    r_data_p_rho->Scale(1/r_data_p_rho->Integral());

    data_p_Jet_mass->Scale(1/data_p_Jet_mass->Integral());
    r_data_p_Jet_mass->Scale(1/r_data_p_Jet_mass->Integral());

    data_p_Jet_chHEF->Scale(1/data_p_Jet_chHEF->Integral());
    r_data_p_Jet_chHEF->Scale(1/r_data_p_Jet_chHEF->Integral());

    TH1* rho_ratio = new TH1F("rho_ratio" ,"rho_ratio" ,100 , 0, 100);
    rho_ratio->Divide(r_data_p_rho,data_p_rho);

    TF1 *func = new TF1("func","pol9",0,100);

    rho_ratio->Fit(func);

    func->FixParameter(0,func->GetParameter(0));
    func->FixParameter(1,func->GetParameter(1));
    func->FixParameter(2,func->GetParameter(2));
    func->FixParameter(3,func->GetParameter(3));
    func->FixParameter(4,func->GetParameter(4));
    func->FixParameter(5,func->GetParameter(5));
    func->FixParameter(6,func->GetParameter(6));
    func->FixParameter(7,func->GetParameter(7));
    func->FixParameter(8,func->GetParameter(8));
    func->FixParameter(9,func->GetParameter(9));

    TH1* p_new_rho = new TH1F("p_new_rho" ," " ,100 , 0, 100);
    TH1* p_new_Jet_mass = new TH1F("p_new_Jet_mass" ," " ,100 , 0, 60);
    TH1* p_new_Jet_chHEF = new TH1F("p_new_Jet_chHEF" ," " ,100 , 0, 1.2);

    double rho_weight;
    TBranch *nrho = datatree->Branch("rho_weight",&rho_weight);

    Long64_t nentries = datatree->GetEntries();

    double rho_w = 0;
    double rho_func = 0;

    for (Long64_t ientry=0; ientry < nentries; ientry++)
    {
        datatree->GetEntry(ientry);
        
        if (irho != (-100))        
        {
            rho_w = irho;
        }

        rho_func = func->Eval(rho_w);
        
        if (rho_func < 0)
        {
            rho_func = 0;
        }    

        p_new_rho->Fill(rho_w,rho_func); 
        p_new_Jet_mass->Fill(iJet_mass,rho_func); 
        p_new_Jet_chHEF->Fill(iJet_chHEF,rho_func);
  
        rho_weight = rho_func;
        nrho->Fill();
    } 

    datatree->Write();
    datatree->Print();
    
/*
    
    TCanvas* c01 = new TCanvas ("c01","c01");
    TPaveLabel* ttitle = new TPaveLabel(0.1,0.962,0.9,0.995,"rho","BL");
    ttitle->SetFillColor(kRed-7);
    ttitle->Draw();
    c01->Divide(1,2);
    c01->cd(1);
    gPad->SetPad(0.0,0.0,1.0,0.8);
    gPad->SetTickx();
    gPad->SetTopMargin(0);
    p_new_rho->Scale(1/p_new_rho->Integral());
    p_new_rho->SetLineColor(kBlue);
    p_new_rho->Draw("hist");
    data_p_rho->SetLineColor(kRed);    
    data_p_rho->Draw("hist same");
    r_data_p_rho->SetLineColor(kBlack);
    r_data_p_rho->Draw("same");
    TLegend* leg_rho = new TLegend(0.6,0.6,0.8,0.8);
    leg_rho->SetTextSize(0.03);
    leg_rho->AddEntry(p_new_rho,"reweighted MC","l");
    leg_rho->AddEntry(data_p_rho,"unweighted MC","l");
    leg_rho->AddEntry(r_data_p_rho,"data","LEP");
    leg_rho->SetBorderSize(0); 
    leg_rho->Draw(); 
    c01->cd(2);
    gPad->SetPad(0.0,0.805,1.0,0.96);
    gPad->SetTickx();
    gPad->SetGridy();
    gPad->SetGridx();
    gPad->SetBottomMargin(0);
    TH1* p_data_mc_rho = new TH1F("" ,"" ,100, 0, 100);
    p_data_mc_rho->Divide(r_data_p_rho,p_new_rho);
    p_data_mc_rho->GetYaxis()->SetNdivisions(404);
    p_data_mc_rho->SetStats(0);
    p_data_mc_rho->GetYaxis()->SetLabelSize(0.13);
    p_data_mc_rho->Draw();    
    c01->Draw();
    
    TCanvas* c1 = new TCanvas ("c1","c1");
    TPaveLabel* title = new TPaveLabel(0.1,0.962,0.9,0.995,"Jet_mass","BL");
    title->SetFillColor(kRed-7);
    title->Draw();
    c1->Divide(1,2);
    c1->cd(1);
    gPad->SetPad(0.0,0.0,1.0,0.8);
    gPad->SetTickx();
    gPad->SetTopMargin(0);
    p_new_Jet_mass->Scale(1/p_new_Jet_mass->Integral());
    p_new_Jet_mass->SetLineColor(kBlue);
    p_new_Jet_mass->Draw("hist");
    data_p_Jet_mass->SetLineColor(kRed);
    data_p_Jet_mass->Draw("hist same");
    r_data_p_Jet_mass->SetLineColor(kBlack);
    r_data_p_Jet_mass->Draw("same");
    TLegend* leg_Jet_mass = new TLegend(0.6,0.6,0.8,0.8);    
    leg_Jet_mass->SetTextSize(0.03);
    leg_Jet_mass->AddEntry(p_new_Jet_mass,"reweighted MC","l");
    leg_Jet_mass->AddEntry(data_p_Jet_mass,"unweighted MC","l");
    leg_Jet_mass->AddEntry(r_data_p_Jet_mass,"data","LEP");
    leg_Jet_mass->SetBorderSize(0);  
    leg_Jet_mass->Draw();
    c1->cd(2);
    gPad->SetPad(0.0,0.805,1.0,0.96);
    gPad->SetTickx();
    gPad->SetGridy();
    gPad->SetGridx();
    gPad->SetBottomMargin(0);
    TH1* p_data_mc_Jet_mass = new TH1F("" ,"" ,100, 0, 60);
    p_data_mc_Jet_mass->Divide(r_data_p_Jet_mass,p_new_Jet_mass);
    p_data_mc_Jet_mass->GetYaxis()->SetNdivisions(404);
    p_data_mc_Jet_mass->SetStats(0);
    p_data_mc_Jet_mass->GetYaxis()->SetLabelSize(0.13);
    p_data_mc_Jet_mass->Draw(); 
    c1->Draw();

    TCanvas* c2 = new TCanvas ("c2","c2");
    TPaveLabel* title1 = new TPaveLabel(0.1,0.962,0.9,0.995,"Jet_chHEF","BL");
    title1->SetFillColor(kRed-7);
    title1->Draw();
    c2->Divide(1,2);
    c2->cd(1);
    gPad->SetPad(0.0,0.0,1.0,0.8);
    gPad->SetTickx();
    gPad->SetTopMargin(0);
    p_new_Jet_chHEF->Scale(1/p_new_Jet_chHEF->Integral());
    p_new_Jet_chHEF->SetLineColor(kBlue);
    p_new_Jet_chHEF->Draw("hist");
    data_p_Jet_chHEF->SetLineColor(kRed);
    data_p_Jet_chHEF->Draw("hist same");
    r_data_p_Jet_chHEF->SetLineColor(kBlack);
    r_data_p_Jet_chHEF->Draw("same");
    TLegend* leg_Jet_chHEF = new TLegend(0.6,0.6,0.8,0.8);
    leg_Jet_chHEF->SetTextSize(0.03);
    leg_Jet_chHEF->AddEntry(p_new_Jet_chHEF,"reweighted MC","l");
    leg_Jet_chHEF->AddEntry(data_p_Jet_chHEF,"unweighted MC","l");
    leg_Jet_chHEF->AddEntry(r_data_p_Jet_chHEF,"data","LEP");
    leg_Jet_chHEF->SetBorderSize(0);
    leg_Jet_chHEF->Draw();
    c2->cd(2);
    gPad->SetPad(0.0,0.805,1.0,0.96);
    gPad->SetTickx();
    gPad->SetGridy();
    gPad->SetGridx();
    gPad->SetBottomMargin(0);
    TH1* p_data_mc_Jet_chHEF = new TH1F("" ,"" ,100, 0, 1.2);
    p_data_mc_Jet_chHEF->Divide(r_data_p_Jet_chHEF,p_new_Jet_chHEF);
    p_data_mc_Jet_chHEF->GetYaxis()->SetNdivisions(404);
    p_data_mc_Jet_chHEF->SetStats(0);
    p_data_mc_Jet_chHEF->GetYaxis()->SetLabelSize(0.13);
    p_data_mc_Jet_chHEF->Draw();
    c2->Draw();
*/
}
/*

    TH1F *p_ratio_rho = new TH1F("p_ratio_rho" ,"DATA/MC_rho" ,100 , 0, 100);
    p_ratio_rho->Divide(r_data_p_rho,data_p_rho);
    

    TF1 *func = new TF1("func","pol9",0,100);

    p_ratio_rho->Fit(func);

    func->FixParameter(0,func->GetParameter(0));
    func->FixParameter(1,func->GetParameter(1));
    func->FixParameter(2,func->GetParameter(2));
    func->FixParameter(3,func->GetParameter(3));
    func->FixParameter(4,func->GetParameter(4));
    func->FixParameter(5,func->GetParameter(5));
    func->FixParameter(6,func->GetParameter(6));
    func->FixParameter(7,func->GetParameter(7));
    func->FixParameter(8,func->GetParameter(8));
    func->FixParameter(9,func->GetParameter(9));

    double p0 = func->GetParameter(0);
    double p1 = func->GetParameter(1);
    double p2 = func->GetParameter(2);
    double p3 = func->GetParameter(3);
    double p4 = func->GetParameter(4);
    double p5 = func->GetParameter(5);
    double p6 = func->GetParameter(6);  
    double p7 = func->GetParameter(7);
    double p8 = func->GetParameter(8);
    double p9 = func->GetParameter(9);
 
    TH1F *p_new_rho = new TH1F("p_new_rho" ,"new rho" ,100 , 0, 100);
    TH1F *p_new1_rho = new TH1F("p_new1_rho" ,"new1 rho" ,100 , 0, 100);

    Long64_t nentries = datatree->GetEntries();

    for (Long64_t ientry=0; ientry < nentries; ientry++)
    {
        datatree->GetEntry(ientry);
        p_new_rho->Fill(rho,func->Eval(rho)); 
        p_new1_rho->Fill(rho);//func->Eval(rho));
    }



    TCanvas* histcanvas = new TCanvas("histcanvas","histcanvas");
//    p_new_rho->Draw("hist");
//    p_new_rho->Scale(1/p_new_rho->Integral());
    p_new1_rho->Draw("same");
//    r_data_p_rho->Draw("same");
//    data_p_rho->Draw("same");
//    p_ratio_rho->Draw();
//    hist1->Draw("same");
    histcanvas->Draw();

//    hist->Divide(hist1);
   
    TF1 *func = new TF1("func","pol9",-100,100);

    hist->Fit(func);

    func->FixParameter(0,func->GetParameter(0));
    func->FixParameter(1,func->GetParameter(1));
    func->FixParameter(2,func->GetParameter(2));
    func->FixParameter(3,func->GetParameter(3));
    func->FixParameter(4,func->GetParameter(4));
    func->FixParameter(5,func->GetParameter(5));
    func->FixParameter(6,func->GetParameter(6));
    func->FixParameter(7,func->GetParameter(7));
    func->FixParameter(8,func->GetParameter(8));
    func->FixParameter(9,func->GetParameter(9));

    double p0 = func->GetParameter(0);
    double p1 = func->GetParameter(1);
    double p2 = func->GetParameter(2);
    double p3 = func->GetParameter(3);
    double p4 = func->GetParameter(4);
    double p5 = func->GetParameter(5);
    double p6 = func->GetParameter(6);  
    double p7 = func->GetParameter(7);
    double p8 = func->GetParameter(8);
    double p9 = func->GetParameter(9);

    string str = "1";
    
//    RooFormulaVar wfunc("w","event weight","(5+10*rho)",rho);
    RooFormulaVar wfunc("w","event weight","(0.746827-0.329127*rho+0.138822*rho*rho-0.0216198*rho*rho*rho+0.00182586*rho*rho*rho*rho-0.0000927829*rho*rho*rho*rho*rho+0.00000289744*rho*rho*rho*rho*rho*rho-0.0000000542824*rho*rho*rho*rho*rho*rho*rho+0.000000000559129*rho*rho*rho*rho*rho*rho*rho*rho-0.00000000000243266*rho*rho*rho*rho*rho*rho*rho*rho*rho)",rho);

    RooRealVar* w = (RooRealVar*) data_set.addColumn(wfunc) ;

    TH1* hh_data = data_set.createHistogram("rho,w",100,100) ;

    TCanvas* histcanvas = new TCanvas("histcanvas","histcanvas");
    hh_data->Draw("hist");
//    hist1->Draw("same");
    histcanvas->Draw();

 
    RooDataSet* jj = &data_set;
    RooDataSet wdata_set(jj->GetName(),jj->GetTitle(),jj,*jj->get(),0,w->GetName());

    wdata_set.Print() ;

   
    RooDataHist fit_hist("fitted rho","fitted rho",rho,hist,1.0);

    RooRealVar poly_c1("poly_c1","coefficient of x^1 term",0,-10,10);
    RooRealVar poly_c2("poly_c2","coefficient of x^2 term",0,-10,10);
    RooRealVar poly_c3("poly_c3","coefficient of x^3 term",0,-10,10);
    RooRealVar poly_c4("poly_c4","coefficient of x^4 term",0,-10,10);
    RooRealVar poly_c5("poly_c5","coefficient of x^5 term",0,-10,10);
    RooRealVar poly_c6("poly_c6","coefficient of x^6 term",0,-10,10);
    RooRealVar poly_c7("poly_c7","coefficient of x^7 term",0,-10,10);
    RooRealVar poly_c8("poly_c8","coefficient of x^8 term",0,-10,10);
    RooRealVar poly_c9("poly_c9","coefficient of x^9 term",0,-10,10);

    RooPolynomial model_hist("model_hist","model_hist",rho,RooArgList(poly_c1,poly_c2,poly_c3,poly_c4,poly_c5,poly_c6,poly_c7,poly_c8));
    model_hist.fitTo(fit_hist);

   

    TCanvas* canvas = new TCanvas("title","0");
    RooPlot* wframe_rho = w->frame(-2,2,100);
    data_set.plotOn(wframe_rho,MarkerSize(0.3),MarkerColor(kBlack),LineColor(kBlack));
//    wdata_set.plotOn(wframe_rho,MarkerSize(0.3),MarkerColor(kRed),LineColor(kRed));
//    r_data_set.plotOn(wframe_rho,MarkerSize(0.3),MarkerColor(kGreen),LineColor(kGreen));
    wframe_rho->Draw();
    canvas->Draw();
*/
//}
    
 






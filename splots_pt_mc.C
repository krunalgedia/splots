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
#include <TObject.h>
#include "TLatex.h"

// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;

// see below for implementation
void AddModel(RooWorkspace*);
void plot(RooPlot* frame_var_name_pt,RooRealVar var_name,double x_low, double x_high,double bins, TString title,RooDataHist r_p_pt_bdt,RooAddPdf totalPdf_pt,RooHistPdf pdf_b_pt,RooHistPdf pdf_n_pt,TCanvas* TCanvas_bdt_pt,TLegend* leg_pt);

void splots_pt_mc()
{
    RooWorkspace* wspace = new RooWorkspace("mywspace");
   // add the signal and background models to the workspace.
    AddModel(wspace);
}


void AddModel(RooWorkspace* ws)
{

    TFile* file1 = new TFile("TMVA_output_splots.root");
    TTree* datatree = (TTree*)(file1->Get("tree"));

    TFile* r_file1 = new TFile("r_TMVA_output_splots.root");
    TTree* r_datatree = (TTree*)(r_file1->Get("tree"));
  
    TH1F* r_data_p_bdt = (TH1F*)(r_file1->Get("r_p_bdt"));
    TH1F* data_p_bdt = (TH1F*)(file1->Get("p_bdt"));
    TH1F* data_p_b_bdt = (TH1F*)(file1->Get("p_b_bdt"));
    TH1F* data_p_n_bdt = (TH1F*)(file1->Get("p_n_bdt"));

    TH1F* r_data_p_20_bdt = (TH1F*)(r_file1->Get("r_p_20_bdt"));
    TH1F* data_p_20_bdt = (TH1F*)(file1->Get("p_20_bdt"));
    TH1F* data_p_b_20_bdt = (TH1F*)(file1->Get("p_b_20_bdt"));
    TH1F* data_p_n_20_bdt = (TH1F*)(file1->Get("p_n_20_bdt"));
    TH1F* r_data_p_30_bdt = (TH1F*)(r_file1->Get("r_p_30_bdt"));
    TH1F* data_p_b_30_bdt = (TH1F*)(file1->Get("p_b_30_bdt"));
    TH1F* data_p_30_bdt = (TH1F*)(file1->Get("p_30_bdt"));
    TH1F* data_p_n_30_bdt = (TH1F*)(file1->Get("p_n_30_bdt"));
    TH1F* r_data_p_40_bdt = (TH1F*)(r_file1->Get("r_p_40_bdt"));
    TH1F* data_p_40_bdt = (TH1F*)(file1->Get("p_40_bdt"));
    TH1F* data_p_b_40_bdt = (TH1F*)(file1->Get("p_b_40_bdt"));
    TH1F* data_p_n_40_bdt = (TH1F*)(file1->Get("p_n_40_bdt"));
    TH1F* r_data_p_60_bdt = (TH1F*)(r_file1->Get("r_p_60_bdt"));
    TH1F* data_p_60_bdt = (TH1F*)(file1->Get("p_60_bdt"));
    TH1F* data_p_b_60_bdt = (TH1F*)(file1->Get("p_b_60_bdt"));
    TH1F* data_p_n_60_bdt = (TH1F*)(file1->Get("p_n_60_bdt"));
    TH1F* r_data_p_80_bdt = (TH1F*)(r_file1->Get("r_p_80_bdt"));
    TH1F* data_p_80_bdt = (TH1F*)(file1->Get("p_80_bdt"));
    TH1F* data_p_b_80_bdt = (TH1F*)(file1->Get("p_b_80_bdt"));
    TH1F* data_p_n_80_bdt = (TH1F*)(file1->Get("p_n_80_bdt"));
    TH1F* r_data_p_100_bdt = (TH1F*)(r_file1->Get("r_p_100_bdt"));
    TH1F* data_p_100_bdt = (TH1F*)(file1->Get("p_100_bdt"));
    TH1F* data_p_b_100_bdt = (TH1F*)(file1->Get("p_b_100_bdt"));
    TH1F* data_p_n_100_bdt = (TH1F*)(file1->Get("p_n_100_bdt"));
    TH1F* r_data_p_120_bdt = (TH1F*)(r_file1->Get("r_p_120_bdt"));
    TH1F* data_p_120_bdt = (TH1F*)(file1->Get("p_120_bdt"));
    TH1F* data_p_b_120_bdt = (TH1F*)(file1->Get("p_b_120_bdt"));
    TH1F* data_p_n_120_bdt = (TH1F*)(file1->Get("p_n_120_bdt"));

    RooRealVar BDT_response("BDT_response", "BDT_response",-1,1);
    RooRealVar BDT_response_20("BDT_response_20", "BDT_response_20",-1,1);
    RooRealVar BDT_response_30("BDT_response_30", "BDT_response_30",-1,1);
    RooRealVar BDT_response_40("BDT_response_40", "BDT_response_40",-1,1);
    RooRealVar BDT_response_60("BDT_response_60", "BDT_response_60",-1,1);
    RooRealVar BDT_response_80("BDT_response_80", "BDT_response_80",-1,1);
    RooRealVar BDT_response_100("BDT_response_100", "BDT_response_100",-1,1);
    RooRealVar BDT_response_120("BDT_response_120", "BDT_response_120",-1,1);

    double max = std::numeric_limits<double>::infinity();

    RooDataHist r_p_bdt("r_p_bdt","r_p_bdt",RooArgList(BDT_response),Import(*r_data_p_bdt));
    RooDataHist p_bdt("p_bdt","p_bdt",RooArgList(BDT_response),Import(*data_p_bdt));
    RooDataHist p_b_bdt("p_b_bdt","p_b_bdt",RooArgList(BDT_response),Import(*data_p_b_bdt));
    RooDataHist p_n_bdt("p_n_bdt","p_n_bdt",RooArgList(BDT_response),Import(*data_p_n_bdt));

    RooDataHist r_p_20_bdt("r_p_20_bdt","r_p_20_bdt",RooArgList(BDT_response),Import(*r_data_p_20_bdt));
    RooDataHist p_20_bdt("p_20_bdt","p_20_bdt",RooArgList(BDT_response),Import(*data_p_20_bdt));
    RooDataHist p_b_20_bdt("p_b_20_bdt","p_b_20_bdt",RooArgList(BDT_response),Import(*data_p_b_20_bdt));
    RooDataHist p_n_20_bdt("p_n_20_bdt","p_n_20_bdt",RooArgList(BDT_response),Import(*data_p_n_20_bdt));
    RooDataHist r_p_30_bdt("r_p_30_bdt","r_p_30_bdt",RooArgList(BDT_response),Import(*r_data_p_30_bdt));
    RooDataHist p_30_bdt("p_30_bdt","p_30_bdt",RooArgList(BDT_response),Import(*data_p_30_bdt));
    RooDataHist p_b_30_bdt("p_b_30_bdt","p_b_30_bdt",RooArgList(BDT_response),Import(*data_p_b_30_bdt));
    RooDataHist p_n_30_bdt("p_n_30_bdt","p_n_30_bdt",RooArgList(BDT_response),Import(*data_p_n_30_bdt));
    RooDataHist r_p_40_bdt("r_p_40_bdt","r_p_40_bdt",RooArgList(BDT_response),Import(*r_data_p_40_bdt));
    RooDataHist p_40_bdt("p_40_bdt","p_40_bdt",RooArgList(BDT_response),Import(*data_p_40_bdt));
    RooDataHist p_b_40_bdt("p_b_40_bdt","p_b_40_bdt",RooArgList(BDT_response),Import(*data_p_b_40_bdt));
    RooDataHist p_n_40_bdt("p_n_40_bdt","p_n_40_bdt",RooArgList(BDT_response),Import(*data_p_n_40_bdt));
    RooDataHist r_p_60_bdt("r_p_60_bdt","r_p_60_bdt",RooArgList(BDT_response),Import(*r_data_p_60_bdt));
    RooDataHist p_60_bdt("p_60_bdt","p_60_bdt",RooArgList(BDT_response),Import(*data_p_60_bdt));
    RooDataHist p_b_60_bdt("p_b_60_bdt","p_b_60_bdt",RooArgList(BDT_response),Import(*data_p_b_60_bdt));
    RooDataHist p_n_60_bdt("p_n_60_bdt","p_n_60_bdt",RooArgList(BDT_response),Import(*data_p_n_60_bdt));
    RooDataHist r_p_80_bdt("r_p_80_bdt","r_p_80_bdt",RooArgList(BDT_response),Import(*r_data_p_80_bdt));
    RooDataHist p_80_bdt("p_80_bdt","p_80_bdt",RooArgList(BDT_response),Import(*data_p_80_bdt));
    RooDataHist p_b_80_bdt("p_b_80_bdt","p_b_80_bdt",RooArgList(BDT_response),Import(*data_p_b_80_bdt));
    RooDataHist p_n_80_bdt("p_n_80_bdt","p_n_80_bdt",RooArgList(BDT_response),Import(*data_p_n_80_bdt));
    RooDataHist r_p_100_bdt("r_p_100_bdt","r_p_100_bdt",RooArgList(BDT_response),Import(*r_data_p_100_bdt));
    RooDataHist p_100_bdt("p_100_bdt","p_100_bdt",RooArgList(BDT_response),Import(*data_p_100_bdt));
    RooDataHist p_b_100_bdt("p_b_100_bdt","p_b_100_bdt",RooArgList(BDT_response),Import(*data_p_b_100_bdt));
    RooDataHist p_n_100_bdt("p_n_100_bdt","p_n_100_bdt",RooArgList(BDT_response),Import(*data_p_n_100_bdt));
    RooDataHist r_p_120_bdt("r_p_120_bdt","r_p_120_bdt",RooArgList(BDT_response),Import(*r_data_p_120_bdt));
    RooDataHist p_120_bdt("p_120_bdt","p_120_bdt",RooArgList(BDT_response),Import(*data_p_120_bdt));
    RooDataHist p_b_120_bdt("p_b_120_bdt","p_b_120_bdt",RooArgList(BDT_response),Import(*data_p_b_120_bdt));
    RooDataHist p_n_120_bdt("p_n_120_bdt","p_n_120_bdt",RooArgList(BDT_response),Import(*data_p_n_120_bdt));

    RooRealVar yield_b("yield_b","yield_b",100,0,max);
    RooRealVar yield_n("yield_n","yield_n",100,0,max);

    RooRealVar yield_b_20("yield_b_20","yield_b_20",100,0,max);
    RooRealVar yield_n_20("yield_n_20","yield_n_20",100,0,max);
    RooRealVar yield_b_30("yield_b_30","yield_b_30",100,0,max);
    RooRealVar yield_n_30("yield_n_30","yield_n_30",100,0,max);
    RooRealVar yield_b_40("yield_b_40","yield_b_40",100,0,max);
    RooRealVar yield_n_40("yield_n_40","yield_n_40",100,0,max);
    RooRealVar yield_b_60("yield_b_60","yield_b_60",100,0,max);
    RooRealVar yield_n_60("yield_n_60","yield_n_60",100,0,max);
    RooRealVar yield_b_80("yield_b_80","yield_b_80",100,0,max);
    RooRealVar yield_n_80("yield_n_80","yield_n_80",100,0,max);
    RooRealVar yield_b_100("yield_b_100","yield_b_100",100,0,max);
    RooRealVar yield_n_100("yield_n_100","yield_n_100",100,0,max);
    RooRealVar yield_b_120("yield_b_120","yield_b_120",100,0,max);
    RooRealVar yield_n_120("yield_n_120","yield_n_120",100,0,max);

    RooHistPdf pdf_b("pdf_b","",BDT_response,p_b_bdt);
    RooHistPdf pdf_n("pdf_n","",BDT_response,p_n_bdt);

    RooHistPdf pdf_b_20("pdf_b_20","",BDT_response,p_b_20_bdt);
    RooHistPdf pdf_n_20("pdf_n_20","",BDT_response,p_n_20_bdt);
    RooHistPdf pdf_b_30("pdf_b_30","",BDT_response,p_b_30_bdt);
    RooHistPdf pdf_n_30("pdf_n_30","",BDT_response,p_n_30_bdt);
    RooHistPdf pdf_b_40("pdf_b_40","",BDT_response,p_b_40_bdt);
    RooHistPdf pdf_n_40("pdf_n_40","",BDT_response,p_n_40_bdt);
    RooHistPdf pdf_b_60("pdf_b_60","",BDT_response,p_b_60_bdt);
    RooHistPdf pdf_n_60("pdf_n_60","",BDT_response,p_n_60_bdt);
    RooHistPdf pdf_b_80("pdf_b_80","",BDT_response,p_b_80_bdt);
    RooHistPdf pdf_n_80("pdf_n_80","",BDT_response,p_n_80_bdt);
    RooHistPdf pdf_b_100("pdf_b_100","",BDT_response,p_b_100_bdt);
    RooHistPdf pdf_n_100("pdf_n_100","",BDT_response,p_n_100_bdt);
    RooHistPdf pdf_b_120("pdf_b_120","",BDT_response,p_b_120_bdt);
    RooHistPdf pdf_n_120("pdf_n_120","",BDT_response,p_n_120_bdt);

    RooArgList shapes_30;
    RooArgList yields_30;

    shapes_30.add(pdf_b_30);
    yields_30.add(yield_b_30);
    shapes_30.add(pdf_n_30);
    yields_30.add(yield_n_30);

    RooAddPdf totalPdf_30("totalPdf_30","totalPdf_30",shapes_30,yields_30);
    totalPdf_30.fitTo(p_30_bdt,Extended());

    RooAddPdf totalPdf_40("totalPdf_40","totalPdf_40",RooArgList(pdf_b_40,pdf_n_40),RooArgList(yield_b_40,yield_n_40));
    totalPdf_40.fitTo(p_40_bdt,Extended());

    RooAddPdf totalPdf_60("totalPdf_60","totalPdf_60",RooArgList(pdf_b_60,pdf_n_60),RooArgList(yield_b_60,yield_n_60));
    totalPdf_60.fitTo(p_60_bdt,Extended());

    RooAddPdf totalPdf_80("totalPdf_80","totalPdf_80",RooArgList(pdf_b_80,pdf_n_80),RooArgList(yield_b_80,yield_n_80));
    totalPdf_80.fitTo(p_80_bdt,Extended());

    RooAddPdf totalPdf_100("totalPdf_100","totalPdf_100",RooArgList(pdf_b_100,pdf_n_100),RooArgList(yield_b_100,yield_n_100));
    totalPdf_100.fitTo(p_100_bdt,Extended());

    RooAddPdf totalPdf_120("totalPdf_120","totalPdf_120",RooArgList(pdf_b_120,pdf_n_120),RooArgList(yield_b_120,yield_n_120));
    totalPdf_120.fitTo(p_120_bdt,Extended());

    RooAddPdf totalPdf("totalPdf","totalPdf",RooArgList(pdf_b,pdf_n),RooArgList(yield_b,yield_n));
    totalPdf.fitTo(p_bdt,Extended());

    RooPlot frame_bdt;
    TCanvas TCanvas_bdt;
    TLegend leg_bdt;

    plot(&frame_bdt,BDT_response,-1,1,100, TString("#scale[1.2]{#font[12]{BDT response for p_{T}-inclusive sample}}"),p_bdt,totalPdf,pdf_b,pdf_n,&TCanvas_bdt,&leg_bdt);


    RooPlot frame_bdt_30;
    TCanvas TCanvas_bdt_30;
    TLegend leg_30;

    plot(&frame_bdt_30,BDT_response,-1,1,100, TString("#scale[1.2]{#font[12]{BDT response for 30<jet p_{T}<40}}"),p_30_bdt,totalPdf_30,pdf_b_30,pdf_n_30,&TCanvas_bdt_30,&leg_30);

    RooPlot frame_bdt_40;
    TCanvas TCanvas_bdt_40;
    TLegend leg_40;

    plot(&frame_bdt_40,BDT_response,-1,1,100, TString("#scale[1.2]{#font[12]{BDT response for 40<jet p_{T}<60}}"),p_40_bdt,totalPdf_40,pdf_b_40,pdf_n_40,&TCanvas_bdt_40,&leg_40);

    RooPlot frame_bdt_60;
    TCanvas TCanvas_bdt_60;
    TLegend leg_60;

    plot(&frame_bdt_60,BDT_response,-1,1,100, TString("#scale[1.2]{#font[12]{BDT response for 60<jet p_{T}<80}}"),p_60_bdt,totalPdf_60,pdf_b_60,pdf_n_60,&TCanvas_bdt_60,&leg_60);

    RooPlot frame_bdt_80;
    TCanvas TCanvas_bdt_80;
    TLegend leg_80;

    plot(&frame_bdt_80,BDT_response,-1,1,100, TString("#scale[1.2]{#font[12]{BDT response for 80<jet p_{T}<100}}"),p_80_bdt,totalPdf_80,pdf_b_80,pdf_n_80,&TCanvas_bdt_80,&leg_80);

    RooPlot frame_bdt_100;
    TCanvas TCanvas_bdt_100;
    TLegend leg_100;

    plot(&frame_bdt_100,BDT_response,-1,1,100, TString("#scale[1.2]{#font[12]{BDT response for 100<jet p_{T}<120}}"),p_100_bdt,totalPdf_100,pdf_b_100,pdf_n_100,&TCanvas_bdt_100,&leg_100);

    RooPlot frame_bdt_120;
    TCanvas TCanvas_bdt_120;
    TLegend leg_120;

    plot(&frame_bdt_120,BDT_response,-1,1,100, TString("#scale[1.2]{#font[12]{BDT response for jet p_{T}>120}}"),p_120_bdt,totalPdf_120,pdf_b_120,pdf_n_120,&TCanvas_bdt_120,&leg_120);

}

    void plot(RooPlot* frame_var_name_pt,RooRealVar var_name,double x_low, double x_high,double bins, TString title,RooDataHist r_p_pt_bdt,RooAddPdf totalPdf_pt,RooHistPdf pdf_b_pt,RooHistPdf pdf_n_pt,TCanvas* TCanvas_bdt_pt,TLegend* leg_pt)
{
    frame_var_name_pt = var_name.frame(x_low, x_high, bins);
    frame_var_name_pt->SetTitle("");
    frame_var_name_pt->GetXaxis()->CenterTitle();
    r_p_pt_bdt.plotOn(frame_var_name_pt,MarkerSize(0.7),Name("data"));
    totalPdf_pt.paramOn(frame_var_name_pt,Format("NEU",AutoPrecision(2)),Layout(0.2,0.5,0.7));
    frame_var_name_pt->getAttText()->SetTextSize(0.03);
    frame_var_name_pt->getAttText()->SetTextColor(kRed+3);
    totalPdf_pt.plotOn(frame_var_name_pt,Components(pdf_b_pt),LineColor(kGreen+2),LineStyle(kDashed),Name("signal"));
    totalPdf_pt.plotOn(frame_var_name_pt,Components(pdf_n_pt),LineColor(kBlue),LineStyle(kDashed),Name("bcg"));
    totalPdf_pt.plotOn(frame_var_name_pt,Components(totalPdf_pt),LineColor(kRed),LineWidth(1),Name("MC"));
    TCanvas_bdt_pt = new TCanvas(title,title);
    double bin_width = frame_var_name_pt->GetXaxis()->GetBinWidth(0);
    std::ostringstream strs;
    strs << bin_width;
    std::string str = strs.str();
    str.erase (str.find_last_not_of('0') + 1, std::string::npos);
    TString yaxis = "jets/" + str;
    frame_var_name_pt->GetYaxis()->SetTitle(yaxis);    
    frame_var_name_pt->GetYaxis()->SetTitleOffset(1.5);
    frame_var_name_pt->GetXaxis()->SetTitleOffset(1.1);
    leg_pt = new TLegend(0.6,0.6,0.8,0.8);
    TObject* signal_leg = (RooPlot*)frame_var_name_pt->findObject("signal");
    TObject* bcg_leg = (RooPlot*)frame_var_name_pt->findObject("bcg");
    TObject* MC_leg = (RooPlot*)frame_var_name_pt->findObject("MC");
    leg_pt->SetFillColor(kWhite);
    leg_pt->SetTextSize(0.03);
    leg_pt->AddEntry(signal_leg,"bjet","l");
    leg_pt->AddEntry(bcg_leg,"non-bjet","l");
    leg_pt->AddEntry(MC_leg,"total MC","l");
    leg_pt->AddEntry("data","MC Truth","LEP");
    leg_pt->SetBorderSize(0);
    TLatex *texf = new TLatex(0.5,0.5,title);
    texf->SetTextSize(0.03375527);
    frame_var_name_pt->Draw();
    leg_pt->Draw();
    texf->Draw();
    TCanvas_bdt_pt->Draw();
}









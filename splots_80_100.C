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
#include <TPaveLabel.h>
#include <RooAbsDataStore.h>
#include <TImage.h>
#include <TAttText.h>
#include <TTF.h>
#include <THStack.h>
#include <TPave.h>
#include <TPaveText.h>

// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;
using namespace std;

// see below for implementation
void AddModel(RooWorkspace*);
void weight_plot(RooPlot* frame_var_name,RooRealVar var_name,double x_low,double x_high,double bins,TString title,RooDataSet* data_set,TCanvas* TCanvas_var_name,TLegend* leg_var_name);
void plots(RooRealVar var_name,RooDataSet* data_set_sig_bg,RooDataSet* r_data_set_sig_bg,RooDataSet* p_data_set_sig_bg,RooDataSet* p_r_data_set_sig_bg,RooDataSet* mc_truth_sig_bg,RooDataSet* p_data_set_b_n,double bins,double x_low,double x_high);
void sweights(RooRealVar var_name_b,RooRealVar var_name_n,RooDataSet* data_set_mc_r);
void sweights_bdt(RooRealVar var_name_b,RooRealVar var_name_n,RooRealVar var_name_bdt,RooDataSet* data_set_mc_r);


void splots_80_100()
{
    RooWorkspace* wspace = new RooWorkspace("mywspace");
    AddModel(wspace);
}


void AddModel(RooWorkspace* ws)
{

    TFile* file1 = new TFile("TMVA_output_splots.root");
    TTree* datatree = (TTree*)(file1->Get("tree"));

    TFile* r_file1 = new TFile("r_TMVA_output_splots.root");
    TTree* r_datatree = (TTree*)(r_file1->Get("tree"));

    RooRealVar BDT_response("BDT_response", "BDT_response",-1,1);

    RooRealVar m_nl_j("m_nl_j", "m_nl_j",0,30000);
    RooRealVar jet_pt("jet_pt", "jet_pt",0,30000);
    RooRealVar Jet_mcFlavour("Jet_mcFlavour","Jet_mcFlavour",-30,30);
    RooRealVar jet_indices("jet_indices", "jet_indices",-100,100);

    RooRealVar rho_weight("rho_weight","rho_weight",-100,100);
    RooRealVar rho_copy("rho_copy","rho_copy",-200,500);
    RooRealVar Jet_leptonDeltaR("Jet_leptonDeltaR","Jet_leptonDeltaR",-105,5);
    RooRealVar Jet_mass("Jet_mass","Jet_mass",30000);
    RooRealVar Jet_pt("Jet_pt","Jet_pt",80,100);
    RooRealVar Jet_eta("Jet_eta","Jet_eta",-3,3);
    RooRealVar Jet_phi("Jet_phi","Jet_phi",-4,4);
    RooRealVar Jet_energyRing_dR0_neut("Jet_energyRing_dR0_neut","Jet_energyRing_dR0_neut",0,1);
    RooRealVar Jet_energyRing_dR1_neut("Jet_energyRing_dR1_neut","Jet_energyRing_dR1_neut",0,1);
    RooRealVar Jet_energyRing_dR2_neut("Jet_energyRing_dR2_neut","Jet_energyRing_dR2_neut",0,1);
    RooRealVar Jet_energyRing_dR3_neut("Jet_energyRing_dR3_neut","Jet_energyRing_dR3_neut",0,1);
    RooRealVar Jet_energyRing_dR4_neut("Jet_energyRing_dR4_neut","Jet_energyRing_dR4_neut",0,1);
    RooRealVar Jet_energyRing_dR0_ch("Jet_energyRing_dR0_ch","Jet_energyRing_dR0_ch",0,1);
    RooRealVar Jet_energyRing_dR1_ch("Jet_energyRing_dR1_ch","Jet_energyRing_dR1_ch",0,1);
    RooRealVar Jet_energyRing_dR2_ch("Jet_energyRing_dR2_ch","Jet_energyRing_dR2_ch",0,1);
    RooRealVar Jet_energyRing_dR3_ch("Jet_energyRing_dR3_ch","Jet_energyRing_dR3_ch",0,1);
    RooRealVar Jet_energyRing_dR4_ch("Jet_energyRing_dR4_ch","Jet_energyRing_dR4_ch",0,1);
    RooRealVar Jet_energyRing_dR0_em("Jet_energyRing_dR0_em","Jet_energyRing_dR0_em",0,1);
    RooRealVar Jet_energyRing_dR1_em("Jet_energyRing_dR1_em","Jet_energyRing_dR1_em",0,1);
    RooRealVar Jet_energyRing_dR2_em("Jet_energyRing_dR2_em","Jet_energyRing_dR2_em",0,1);
    RooRealVar Jet_energyRing_dR3_em("Jet_energyRing_dR3_em","Jet_energyRing_dR3_em",0,1);
    RooRealVar Jet_energyRing_dR4_em("Jet_energyRing_dR4_em","Jet_energyRing_dR4_em",0,1);
    RooRealVar Jet_energyRing_dR0_mu("Jet_energyRing_dR0_mu","Jet_energyRing_dR0_mu",0,1);
    RooRealVar Jet_energyRing_dR1_mu("Jet_energyRing_dR1_mu","Jet_energyRing_dR1_mu",0,1);
    RooRealVar Jet_energyRing_dR2_mu("Jet_energyRing_dR2_mu","Jet_energyRing_dR2_mu",0,1);
    RooRealVar Jet_energyRing_dR3_mu("Jet_energyRing_dR3_mu","Jet_energyRing_dR3_mu",0,1);
    RooRealVar Jet_energyRing_dR4_mu("Jet_energyRing_dR4_mu","Jet_energyRing_dR4_mu",0,1);
    RooRealVar Jet_rawEnergy("Jet_rawEnergy","Jet_rawEnergy",0,4000);
    RooRealVar Jet_numDaughters_pt03("Jet_numDaughters_pt03","Jet_numDaughters_pt03",0,100);
    RooRealVar Jet_numberOfDaughters("Jet_numberOfDaughters","Jet_numberOfDaughters",0,1200000000);
    RooRealVar rho("rho","rho",-200,100);
    RooRealVar Jet_rawPt("Jet_rawPt","Jet_rawPt",0,1600);
    RooRealVar Jet_chHEF("Jet_chHEF","Jet_chHEF",0,12);
    RooRealVar Jet_muEF("Jet_muEF","Jet_muEF",0,12);
    RooRealVar Jet_neHEF("Jet_neHEF","Jet_neHEF",0,12);
    RooRealVar Jet_chEmEF("Jet_chEmEF","Jet_chEmEF",0,12);
    RooRealVar Jet_neEmEF("Jet_neEmEF","Jet_neEmEF",0,12);
    RooRealVar Jet_leadTrackPt("Jet_leadTrackPt","Jet_leadTrackPt",0,800);
    RooRealVar Jet_leptonPt("Jet_leptonPt","Jet_leptonPt",-200,1000);
    RooRealVar Jet_leptonPtRel("Jet_leptonPtRel","Jet_leptonPtRel",-200,100);
    RooRealVar Jet_leptonPdgId("Jet_leptonPdgId","Jet_leptonPdgId",-200,100);
    RooRealVar Jet_leptonPtRelInv("Jet_leptonPtRelInv","Jet_leptonPtRelInv",-200,200);
    RooRealVar Jet_vtxMass("Jet_vtxMass","Jet_vtxMass",0,10);
    RooRealVar Jet_vtxNtracks("Jet_vtxNtracks","Jet_vtxNtracks",0,16);
    RooRealVar Jet_vtxPt("Jet_vtxPt","Jet_vtxPt",0,600);
    RooRealVar Jet_vtx3DSig("Jet_vtx3DSig","Jet_vtx3DSig",0,400);
    RooRealVar Jet_vtx3DVal("Jet_vtx3DVal","Jet_vtx3DVal",0,18);
    RooRealVar Jet_ptd("Jet_ptd","Jet_ptd",-30,20);

    RooArgSet var_set;
    var_set.add(rho_weight);
    var_set.add(rho_copy);
    var_set.add(Jet_mcFlavour);
    var_set.add(BDT_response);
    var_set.add(Jet_leptonDeltaR);
    var_set.add(Jet_mass);
    var_set.add(Jet_pt);
    var_set.add(Jet_eta);
    var_set.add(Jet_phi);
    var_set.add(Jet_rawEnergy);
    var_set.add(Jet_energyRing_dR0_neut);   
    var_set.add(Jet_energyRing_dR1_neut);
    var_set.add(Jet_energyRing_dR2_neut);
    var_set.add(Jet_energyRing_dR3_neut);
    var_set.add(Jet_energyRing_dR4_neut);
    var_set.add(Jet_energyRing_dR0_ch);
    var_set.add(Jet_energyRing_dR1_ch);
    var_set.add(Jet_energyRing_dR2_ch);
    var_set.add(Jet_energyRing_dR3_ch);
    var_set.add(Jet_energyRing_dR4_ch);
    var_set.add(Jet_energyRing_dR0_em);
    var_set.add(Jet_energyRing_dR1_em);
    var_set.add(Jet_energyRing_dR2_em);
    var_set.add(Jet_energyRing_dR3_em);
    var_set.add(Jet_energyRing_dR4_em);
    var_set.add(Jet_energyRing_dR0_mu);
    var_set.add(Jet_energyRing_dR1_mu);
    var_set.add(Jet_energyRing_dR2_mu);
    var_set.add(Jet_energyRing_dR3_mu);
    var_set.add(Jet_energyRing_dR4_mu);
    var_set.add(Jet_numDaughters_pt03);
    var_set.add(Jet_numberOfDaughters);
    var_set.add(rho);
    var_set.add(Jet_rawPt);
    var_set.add(Jet_chHEF);
    var_set.add(Jet_muEF);
    var_set.add(Jet_neHEF);
    var_set.add(Jet_chEmEF);
    var_set.add(Jet_neEmEF);
    var_set.add(Jet_leadTrackPt);
    var_set.add(Jet_leptonPt);
    var_set.add(Jet_leptonPtRel);
    var_set.add(Jet_leptonPdgId);
    var_set.add(Jet_leptonPtRelInv);
    var_set.add(Jet_vtxMass);
    var_set.add(Jet_vtxNtracks);
    var_set.add(Jet_vtxPt);
    var_set.add(Jet_vtx3DSig);
    var_set.add(Jet_vtx3DVal);
    var_set.add(Jet_ptd);

    RooArgSet r_var_set;
    r_var_set.add(rho_copy);
    r_var_set.add(BDT_response);
    r_var_set.add(Jet_mass);
    r_var_set.add(Jet_leptonDeltaR);
    r_var_set.add(Jet_pt);
    r_var_set.add(Jet_eta);
    r_var_set.add(Jet_phi);
    r_var_set.add(Jet_rawEnergy);
    r_var_set.add(Jet_energyRing_dR0_neut);
    r_var_set.add(Jet_energyRing_dR1_neut);
    r_var_set.add(Jet_energyRing_dR2_neut);
    r_var_set.add(Jet_energyRing_dR3_neut);
    r_var_set.add(Jet_energyRing_dR4_neut);
    r_var_set.add(Jet_energyRing_dR0_ch);
    r_var_set.add(Jet_energyRing_dR1_ch);
    r_var_set.add(Jet_energyRing_dR2_ch);
    r_var_set.add(Jet_energyRing_dR3_ch);
    r_var_set.add(Jet_energyRing_dR4_ch);
    r_var_set.add(Jet_energyRing_dR0_em);
    r_var_set.add(Jet_energyRing_dR1_em);
    r_var_set.add(Jet_energyRing_dR2_em);
    r_var_set.add(Jet_energyRing_dR3_em);
    r_var_set.add(Jet_energyRing_dR4_em);
    r_var_set.add(Jet_energyRing_dR0_mu);
    r_var_set.add(Jet_energyRing_dR1_mu);
    r_var_set.add(Jet_energyRing_dR2_mu);
    r_var_set.add(Jet_energyRing_dR3_mu);
    r_var_set.add(Jet_energyRing_dR4_mu);
    r_var_set.add(Jet_numDaughters_pt03);
    r_var_set.add(Jet_numberOfDaughters);
    r_var_set.add(rho);
    r_var_set.add(Jet_rawPt);
    r_var_set.add(Jet_chHEF);
    r_var_set.add(Jet_muEF);
    r_var_set.add(Jet_neHEF);
    r_var_set.add(Jet_chEmEF);
    r_var_set.add(Jet_neEmEF);
    r_var_set.add(Jet_leadTrackPt);
    r_var_set.add(Jet_leptonPt);
    r_var_set.add(Jet_leptonPtRel);
    r_var_set.add(Jet_leptonPdgId);
    r_var_set.add(Jet_leptonPtRelInv);
    r_var_set.add(Jet_vtxMass);
    r_var_set.add(Jet_vtxNtracks);
    r_var_set.add(Jet_vtxPt);
    r_var_set.add(Jet_vtx3DSig);
    r_var_set.add(Jet_vtx3DVal);
    r_var_set.add(Jet_ptd);

    RooDataSet* x_data_set = new RooDataSet("x_data_set","x_data_set",var_set,Import(*datatree));
    RooDataSet* data_set = new RooDataSet("data_set","data_set",x_data_set,*x_data_set->get(),"",rho_weight.GetName()); 

    RooDataSet* t_mc_truth_sig = new RooDataSet("t_mc_truth_sig","t_mc_truth_sig",x_data_set,*x_data_set->get(),"abs(Jet_mcFlavour)==5");
    RooDataSet* t_mc_truth_bg = new RooDataSet("t_mc_truth_bg","t_mc_truth_bg",x_data_set,*x_data_set->get(),"abs(Jet_mcFlavour)!=5");
  
    RooDataSet* r_data_set = new RooDataSet("r_data_set","r_data_set",r_var_set,Import(*r_datatree));

    RooDataSet* mc_truth_sig = (RooDataSet*) data_set->reduce("abs(Jet_mcFlavour)==5");
    RooDataSet* mc_truth_bg = (RooDataSet*) data_set->reduce("abs(Jet_mcFlavour)!=5");

    TH1* data_p_b_bdt = (TH1*) mc_truth_sig->createHistogram("data_p_b_bdt",BDT_response,Binning(100,-1,1));
    TH1* data_p_n_bdt = (TH1*) mc_truth_bg->createHistogram("data_p_n_bdt",BDT_response,Binning(100,-1,1));

    RooDataHist p_b_bdt("p_b_bdt","p_b_bdt",RooArgList(BDT_response),data_p_b_bdt);
    RooDataHist p_n_bdt("p_n_bdt","p_n_bdt",RooArgList(BDT_response),data_p_n_bdt);

    RooHistPdf pdf_b("pdf_b","",BDT_response,p_b_bdt);
    RooHistPdf pdf_n("pdf_n","",BDT_response,p_n_bdt);
  
    double max = std::numeric_limits<double>::infinity();

    RooRealVar yield_b("yield_b","yield_b",100,0,max);
    RooRealVar yield_n("yield_n","yield_n",100,0,max);
    RooRealVar r_yield_b("r_yield_b","r_yield_b",100,0,max);
    RooRealVar r_yield_n("r_yield_n","r_yield_n",100,0,max);

    RooAddPdf totalPdf("totalPdf","totalPdf",RooArgList(pdf_b,pdf_n),RooArgList(yield_b,yield_n));
//    totalPdf.fitTo(data_set,Extended());

    RooAddPdf r_totalPdf("r_totalPdf","r_totalPdf",RooArgList(pdf_b,pdf_n),RooArgList(r_yield_b,r_yield_n));
//    r_totalPdf.fitTo(r_data_set,Extended());

    RooAddPdf* ptr_totalPdf = &totalPdf;
    RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",*data_set,ptr_totalPdf,RooArgList(yield_b,yield_n));

    RooAddPdf* r_ptr_totalPdf = &r_totalPdf;
    RooStats::SPlot* r_sData = new RooStats::SPlot("r_sData","An r_SPlot",*r_data_set,r_ptr_totalPdf,RooArgList(r_yield_b,r_yield_n));

    cout<<yield_b.getVal()/mc_truth_sig->sumEntries()<<" yield_b.getVal()/mc_truth_sig.sumEntries()"<<endl;
    cout<<yield_n.getVal()/mc_truth_bg->sumEntries()<<" yield_n.getVal()/mc_truth_bg.sumEntries()"<<endl;

    RooRealVar yield_b_sw("yield_b_sw","yield_b_sw",-10,10);
    RooRealVar yield_n_sw("yield_n_sw","yield_n_sw",-10,10);

    RooDataSet* data_set_sig = new RooDataSet("data_set_sig","data_set_sig",data_set,*data_set->get(),"","yield_b_sw");
    RooDataSet* data_set_bg = new RooDataSet("data_set_bg","data_set_bg",data_set,*data_set->get(),"","yield_n_sw");

    RooDataSet* p_data_set_b = (RooDataSet*) data_set->reduce("yield_b_sw>=0");
    RooDataSet* p_data_set_n = (RooDataSet*) data_set->reduce("yield_n_sw>=0");    

    RooDataSet* p_data_set_sig = new RooDataSet("p_data_set_sig","p_data_set_sig",p_data_set_b,*p_data_set_b->get(),"","yield_b_sw");
    RooDataSet* p_data_set_bg = new RooDataSet("p_data_set_bg","p_data_set_bg",p_data_set_n,*p_data_set_n->get(),"","yield_n_sw");

    Double_t p_sw_mc_evt = p_data_set_sig->sumEntries();

    RooRealVar p_sw_mc_events("p_sw_mc_events","p_sw_mc_events",p_sw_mc_evt);
    p_sw_mc_events.setConstant();

    RooFormulaVar p_sw_mc_norm = RooFormulaVar("p_sw_mc_norm","@0/@1",RooArgList(yield_b,p_sw_mc_events));
    p_data_set_b->addColumn(p_sw_mc_norm);

    RooDataSet* t_data_set_sig = new RooDataSet("t_data_set_sig","t_data_set_sig",p_data_set_b,*p_data_set_b->get());
    RooDataSet* t_data_set_bg = new RooDataSet("t_data_set_bg","t_data_set_bg",p_data_set_n,*p_data_set_n->get());

    RooRealVar r_yield_b_sw("r_yield_b_sw","r_yield_b_sw",-10,10);
    RooRealVar r_yield_n_sw("r_yield_n_sw","r_yield_n_sw",-10,10);

    RooDataSet* r_data_set_sig = new RooDataSet("r_data_set_sig","r_data_set_sig",r_data_set,*r_data_set->get(),"","r_yield_b_sw");
    RooDataSet* r_data_set_bg = new RooDataSet("r_data_set_bg","r_data_set_bg",r_data_set,*r_data_set->get(),"","r_yield_n_sw");

    RooDataSet* p_r_data_set_b = (RooDataSet*) r_data_set->reduce("r_yield_b_sw>=0");
    RooDataSet* p_r_data_set_n = (RooDataSet*) r_data_set->reduce("r_yield_n_sw>=0");

    RooDataSet* p_r_data_set_sig = new RooDataSet("p_r_data_set_sig","p_r_data_set_sig",p_r_data_set_b,*p_r_data_set_b->get(),"","r_yield_b_sw");
    RooDataSet* p_r_data_set_bg = new RooDataSet("p_r_data_set_bg","p_r_data_set_bg",p_r_data_set_n,*p_r_data_set_n->get(),"","r_yield_n_sw");

    Double_t p_sw_data_evt = p_r_data_set_sig->sumEntries();

    RooRealVar p_sw_data_events("p_sw_data_events","p_sw_data_events",p_sw_data_evt);
    p_sw_data_events.setConstant();

    RooFormulaVar p_sw_data_norm = RooFormulaVar("p_sw_data_norm","@0/@1",RooArgList(r_yield_b,p_sw_data_events));
    p_r_data_set_b->addColumn(p_sw_data_norm);

    RooDataSet* t_r_data_set_sig = new RooDataSet("t_r_data_set_sig","t_r_data_set_sig",p_r_data_set_b,*p_r_data_set_b->get());
    RooDataSet* t_r_data_set_bg = new RooDataSet("t_r_data_set_bg","t_r_data_set_bg",p_r_data_set_n,*p_r_data_set_n->get());

/*
    sweights(yield_b_sw,yield_n_sw,data_set);
    sweights(r_yield_b_sw,r_yield_n_sw,r_data_set);

    plots(rho,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,60);
    plots(rho,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,60);
    plots(Jet_leptonDeltaR,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,0.2);
    plots(Jet_leptonDeltaR,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,0.2);
    plots(Jet_mass,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,50);
    plots(Jet_mass,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,50);
    plots(Jet_pt,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,300);
    plots(Jet_pt,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,300);
    plots(Jet_eta,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,-2.6,2.6);
    plots(Jet_eta,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,-2.6,2.6);
    plots(Jet_phi,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,-3.2,3.2);
    plots(Jet_phi,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,-3.2,3.2);
    plots(Jet_rawEnergy,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,600);
    plots(Jet_rawEnergy,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,600);
    plots(Jet_energyRing_dR0_neut,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR0_neut,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR1_neut,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR1_neut,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR2_neut,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR2_neut,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR3_neut,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR3_neut,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR4_neut,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR4_neut,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR0_ch,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR0_ch,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR1_ch,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR1_ch,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR2_ch,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR2_ch,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR3_ch,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR3_ch,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR4_ch,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR4_ch,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR0_em,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR0_em,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR1_em,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR1_em,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR2_em,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR2_em,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR3_em,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR3_em,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR4_em,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR4_em,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR0_mu,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR0_mu,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR1_mu,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR1_mu,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR2_mu,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR2_mu,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR3_mu,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR3_mu,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_energyRing_dR4_mu,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_energyRing_dR4_mu,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_numDaughters_pt03,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,100);
    plots(Jet_numDaughters_pt03,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,100);
    plots(Jet_numberOfDaughters,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,80,0,80);
    plots(Jet_numberOfDaughters,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,80,0,80);
    plots(Jet_rawPt,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,400);
    plots(Jet_rawPt,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,400);
    plots(Jet_chHEF,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_chHEF,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_muEF,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_muEF,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_neHEF,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_neHEF,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_chEmEF,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_chEmEF,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_neEmEF,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_neEmEF,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);
    plots(Jet_leadTrackPt,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,140);
    plots(Jet_leadTrackPt,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,140);
    plots(Jet_leptonPt,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,160);
    plots(Jet_leptonPt,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,160);
    plots(Jet_leptonPtRel,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,30);
    plots(Jet_leptonPtRel,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,30);
    plots(Jet_leptonPdgId,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,40,-20,20);
    plots(Jet_leptonPdgId,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,40,-20,20);
    plots(Jet_leptonPtRelInv,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,60);
    plots(Jet_leptonPtRelInv,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,60);
    plots(Jet_vtxMass,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,5);
    plots(Jet_vtxMass,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,5);
    plots(Jet_vtxNtracks,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,15,0,15);
    plots(Jet_vtxNtracks,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,15,0,15);
    plots(Jet_vtxPt,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,100);
    plots(Jet_vtxPt,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,100);
    plots(Jet_vtx3DSig,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,110,0,110);
    plots(Jet_vtx3DSig,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,110,0,110);
    plots(Jet_vtx3DVal,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,5);
    plots(Jet_vtx3DVal,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,5);
    plots(Jet_ptd,data_set_sig,r_data_set_sig,p_data_set_sig,p_r_data_set_sig,mc_truth_sig,p_data_set_b,100,0,1);
    plots(Jet_ptd,data_set_bg,r_data_set_bg,p_data_set_bg,p_r_data_set_bg,mc_truth_bg,p_data_set_n,100,0,1);


    sweights_bdt(yield_b_sw,yield_n_sw,BDT_response,data_set);
    sweights_bdt(r_yield_b_sw,r_yield_n_sw,BDT_response,r_data_set);
*/

//Make TTrees of the datasets used



    RooAbsData::setDefaultStorageType(RooAbsData::Tree);
    TFile* qreg_file = new TFile("qreg_file_80_100.root","RECREATE");

    RooDataSet* dataset = new RooDataSet("dataset","dataset",data_set,*data_set->get()); 
    RooDataSet* mc_sig = new RooDataSet("mc_sig","mc_sig",t_mc_truth_sig,*t_mc_truth_sig->get());
    RooDataSet* mc_bg = new RooDataSet("mc_bg","mc_bg",t_mc_truth_bg,*t_mc_truth_bg->get());    
    RooDataSet* t_dataset_sig = new RooDataSet("t_dataset_sig","t_dataset_sig",t_data_set_sig,*t_data_set_sig->get());
    RooDataSet* t_dataset_bg = new RooDataSet("t_dataset_bg","t_dataset_bg",t_data_set_bg,*t_data_set_bg->get());
    RooDataSet* r_dataset = new RooDataSet("r_dataset","r_dataset",r_data_set,*r_data_set->get());
    RooDataSet* t_r_dataset_sig = new RooDataSet("t_r_dataset_sig","t_r_dataset_sig",t_r_data_set_sig,*t_r_data_set_sig->get());
    RooDataSet* t_r_dataset_bg = new RooDataSet("t_r_dataset_bg","t_r_dataset_bg",t_r_data_set_bg,*t_r_data_set_bg->get());

    const TTree *dataset_tree = dataset->tree();
    const TTree *mc_sig_tree = mc_sig->tree();
    const TTree *mc_bg_tree = mc_bg->tree();
    const TTree *t_dataset_sig_tree = t_dataset_sig->tree();
    const TTree *t_dataset_bg_tree = t_dataset_bg->tree();
    const TTree *r_dataset_tree = r_dataset->tree();
    const TTree *t_r_dataset_sig_tree = t_r_dataset_sig->tree();
    const TTree *t_r_dataset_bg_tree = t_r_dataset_bg->tree();

    dataset_tree->Write();
    mc_sig_tree->Write();
    mc_bg_tree->Write();
    t_dataset_sig_tree->Write();
    t_dataset_bg_tree->Write();
    r_dataset_tree->Write();
    t_r_dataset_sig_tree->Write();
    t_r_dataset_bg_tree->Write();

    qreg_file->Close();  

}

void plots(RooRealVar var_name,RooDataSet* data_set_sig_bg,RooDataSet* r_data_set_sig_bg,RooDataSet* p_data_set_sig_bg,RooDataSet* p_r_data_set_sig_bg,RooDataSet* mc_truth_sig_bg,RooDataSet* p_data_set_b_n,double bins,double x_low,double x_high)


{
    string variable_names = std::string(var_name.GetName());
    string dataset_names = std::string(data_set_sig_bg->GetName());
    string comps;
    if (dataset_names.find("sig")==9)
    {
        comps = " (signal)";
    }
    else
    {
        comps = " (background)";
    }    

    TString variable_name = variable_names;
    TString dataset_name = dataset_names;
    TString comp = comps;  

    TCanvas* c1_var_name = new TCanvas("c1_" + variable_name + comp,"c1_" + variable_name + comp,700,52,900,848);
    c1_var_name->Range(0,0,1,1);
    c1_var_name->SetFillColor(0);
    c1_var_name->SetBorderMode(0);
    c1_var_name->SetBorderSize(2);
    c1_var_name->SetFrameBorderMode(0);
    gStyle->SetPadBorderMode(0);
    TPaveLabel* gtitle = new TPaveLabel(0.09966592,0.956151,0.9003341,0.9890378,variable_name + comp,"BL");
    if (dataset_names.find("sig")==9)
    {
        gtitle->SetFillColor(46);
    }
    else
    {
        gtitle->SetFillColor(47);
    }
    gtitle->Draw();
    c1_var_name->Divide(1,3);
    c1_var_name->cd(1);
    gPad->SetPad(0.0,0.60,1.0,0.95);
    gPad->SetTickx();
    gPad->SetGridx();
    TH1* h1a = (TH1*) data_set_sig_bg->createHistogram("",var_name,Binning(bins,x_low,x_high));
    TH1* r_h1a = (TH1*) r_data_set_sig_bg->createHistogram("",var_name,Binning(bins,x_low,x_high));
    double s1a = h1a->GetSumOfWeights();
    r_h1a->Scale(s1a/r_h1a->GetSumOfWeights());
    double bin_width = h1a->GetXaxis()->GetBinWidth(0);
    std::string str = std::to_string(bin_width);
    str.erase (str.find_last_not_of('0') + 1, std::string::npos );
    TString jets = "jets/";
    h1a->SetLineColor(kBlue);
    r_h1a->SetLineColor(kBlack);

    THStack* f1a = new THStack();
    f1a->Add(h1a,"HIST");
    f1a->Add(r_h1a,"E0");
    f1a->Draw("nostack");
    f1a->GetHistogram()->GetXaxis()->SetLabelOffset(999);
    f1a->GetHistogram()->GetXaxis()->SetTitle("");
    f1a->GetHistogram()->GetYaxis()->SetTitle(jets+str);
    f1a->GetHistogram()->GetYaxis()->SetTitleOffset(0.95);
    f1a->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
    f1a->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
    f1a->GetHistogram()->GetYaxis()->CenterTitle();
    TLegend* l1a = new TLegend(0.7327394,0.728728,0.8830735,0.8783713,NULL,"brNDC");
    l1a->AddEntry(r_h1a,"sweighted data","LEP");
    l1a->AddEntry(h1a,"sweighted MC","l");
    l1a->SetTextFont(12);
    l1a->SetTextSize(0.05);
    l1a->SetBorderSize(0);
    l1a->Draw("same");

    c1_var_name->cd(2);
    gPad->SetPad(0.0,0.20,1.0,0.60);
    gPad->SetTickx();
    gPad->SetGridx();
    gPad->SetBottomMargin(0);
    gPad->SetTopMargin(0);
    TH1* h1b = (TH1*) p_data_set_sig_bg->createHistogram("",var_name,Binning(bins,x_low,x_high));
    TH1* r_h1b = (TH1*) p_r_data_set_sig_bg->createHistogram("",var_name,Binning(bins,x_low,x_high));
    double s1b = h1b->GetSumOfWeights();
    r_h1b->Scale(s1b/r_h1b->Integral());
    h1b->SetLineColor(kGreen+3);
    r_h1b->SetLineColor(kBlack);
    THStack* f1b = new THStack();
    f1b->Add(h1b,"HIST");
    f1b->Add(r_h1b,"E0");
    f1b->Draw("nostack");
    f1b->GetHistogram()->GetXaxis()->SetLabelOffset(999);
    f1b->GetHistogram()->GetXaxis()->SetTitle("");
    TString log = "log(jets)/";

    if (variable_names.find("energyRing")==4 || variable_names.find("rawPt")==4 || variable_names.find("muEF")==4 ||variable_names.find("neHEF")==4 || variable_names.find("chEmEF")==4 || variable_names.find("neEmEF")==4 || variable_names.find("leadTrackPt")==4 || variable_names.find("leptonPt")==4 || variable_names.find("vtxMass")==4 || variable_names.find("vtx3DSig")==4 ||variable_names.find("vtx3DVal")==4 ||variable_names.find("ptd")==4 || variable_names.find("vtxPt")==4)
    {
        gPad->SetLogy();
        f1b->GetHistogram()->GetYaxis()->SetTitle(log+str);
    }
    else 
    {
        f1b->GetHistogram()->GetYaxis()->SetTitle(jets+str);
    }
    f1b->GetHistogram()->GetYaxis()->SetTitleOffset(0.95);
    f1b->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
    f1b->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
    f1b->GetHistogram()->GetYaxis()->CenterTitle();
    TLegend* l1b = new TLegend(0.7048998,0.8337393,0.8552339,0.9829476,NULL,"brNDC");
    l1b->AddEntry(r_h1b,"+ve sweighted data","LEP");
    l1b->AddEntry(h1b,"+ve sweighted MC","l");
    l1b->SetTextFont(12);
    l1b->SetTextSize(0.05);
    l1b->SetBorderSize(0);
    l1b->Draw("same");

    c1_var_name->cd(3);
    gPad->SetPad(0.0,0.0,1.0,0.20);
//    gPad->SetTickx();
    gPad->SetTicky();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTopMargin(0);
    gPad->SetBottomMargin(0.2984166);
    TH1* h1c = (TH1*) data_set_sig_bg->createHistogram("",var_name,Binning(bins,x_low,x_high));
    TH1* r_h1c = (TH1*) r_data_set_sig_bg->createHistogram("",var_name,Binning(bins,x_low,x_high));
    double s1c = h1c->GetSumOfWeights();
    r_h1c->Scale(s1c/r_h1c->GetSumOfWeights());
    r_h1c->Divide(h1c);
    r_h1c->SetStats(0);
    r_h1c->SetLineColor(kBlack);
    r_h1c->SetTitle("");
    r_h1c->GetXaxis()->SetTitle(variable_name);
    r_h1c->GetXaxis()->SetTitleSize(0.12);
    r_h1c->GetXaxis()->SetTitleOffset(0.89);
    r_h1c->GetXaxis()->CenterTitle();
    r_h1c->GetXaxis()->SetLabelSize(0.11);
    r_h1c->GetXaxis()->SetTickLength(0.03);
    r_h1c->GetYaxis()->SetTitle("Data/MC");
    r_h1c->GetYaxis()->CenterTitle();
    r_h1c->GetYaxis()->SetTitleSize(0.08);
    r_h1c->GetYaxis()->SetLabelSize(0.08);
    r_h1c->GetYaxis()->SetTitleOffset(0.40);
    r_h1c->GetYaxis()->SetNdivisions(406);    
//    gPad->SetTicky();
    r_h1c->GetYaxis()->SetRangeUser(0.5,1.5);
    r_h1c->Draw("LEP");
    c1_var_name->Update();


    TCanvas* c2_var_name = new TCanvas("c2_" + variable_name + comp,"c2_" + variable_name + comp,640,52,900,848);
    c2_var_name->Range(0,0,1,1);
    c2_var_name->SetFillColor(0);
    c2_var_name->SetBorderMode(0);
    c2_var_name->SetBorderSize(2);
    c2_var_name->SetFrameBorderMode(0);
    gStyle->SetPadBorderMode(0);
    TPaveLabel* title = new TPaveLabel(0.09966592,0.9488429,0.9003341,0.9817296,variable_name + comp,"BL");
    if (dataset_names.find("sig")==9)
    {
        title->SetFillColor(42);
    }
    else
    {
        title->SetFillColor(41);
    }
    title->Draw();
    c2_var_name->Divide(1,2);
    c2_var_name->cd(1);
    gPad->SetPad(0,0.3,1,0.9378806);
    gPad->SetGridx();
    gPad->SetTopMargin(0.04971319);
    gPad->SetBottomMargin(0.001912046);

    TH1* h1 = (TH1*) mc_truth_sig_bg->createHistogram("",var_name,Binning(bins,x_low,x_high));
    TH1* h2 = (TH1*) data_set_sig_bg->createHistogram("",var_name,Binning(bins,x_low,x_high));
    TH1* h3 = (TH1*) p_data_set_sig_bg->createHistogram("",var_name,Binning(bins,x_low,x_high));
    TH1* h4 = (TH1*) p_data_set_b_n->createHistogram("",var_name,Binning(bins,x_low,x_high));
    double s2a = h1->GetSumOfWeights();
    h2->Scale(s2a/h2->GetSumOfWeights());
    h3->Scale(s2a/h3->GetSumOfWeights());
    h4->Scale(s2a/h4->GetSumOfWeights());
    h1->SetLineColor(kBlue);
    h2->SetLineColor(kRed);
    h3->SetLineColor(kBlack);
    h4->SetLineColor(kGreen+3);
    h1->SetLineStyle(1);
    h2->SetLineStyle(7);
    h3->SetLineStyle(1);
    h4->SetLineStyle(2);
    THStack* hs = new THStack();
    hs->Add(h1,"HIST");
    hs->Add(h2,"HIST");
    hs->Add(h3,"E0");
    hs->Add(h4,"HIST");
    hs->Draw("nostack");
    hs->GetHistogram()->GetXaxis()->SetLabelOffset(999);
    hs->GetHistogram()->GetXaxis()->SetTitle("");
    hs->GetHistogram()->GetYaxis()->SetTitle(jets+str);
    hs->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
    hs->GetHistogram()->GetYaxis()->SetTitleSize(0.03);
    hs->GetHistogram()->GetYaxis()->SetLabelSize(0.03);
    hs->GetHistogram()->GetYaxis()->CenterTitle();
    TLegend* l2a = new TLegend(0.7004454,0.7954111,0.8507795,0.9445507,NULL,"brNDC");
    l2a->AddEntry(h3,"+ve sweighted MC","LEP");
    l2a->AddEntry(h1,"MC Truth","l");
    l2a->AddEntry(h2,"sweighted MC","l");
    l2a->AddEntry(h4,"unwgt (+ve sw) MC","l");
    l2a->SetTextFont(12);
    l2a->SetTextSize(0.03);
    l2a->SetBorderSize(0);
    l2a->Draw("same");


    c2_var_name->cd(2);
    gPad->SetPad(0.0,0.0,1.0,0.3);
    gPad->SetTickx();
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTopMargin(0);
    gPad->SetBottomMargin(0.227365);
    TH1* d1 = (TH1*) p_data_set_sig_bg->createHistogram("",var_name,Binning(bins,x_low,x_high));
    TH1* d2 = (TH1*) mc_truth_sig_bg->createHistogram("",var_name,Binning(bins,x_low,x_high));
    double s2b = d1->GetSumOfWeights();
    d2->Scale(s2b/d2->GetSumOfWeights());
    d1->Divide(d2);
    d1->SetStats(0);
    d1->SetLineColor(kBlack);
    d1->SetTitle("");
    d1->GetXaxis()->SetLabelSize(0.06);
    d1->GetXaxis()->SetTitle(variable_name);
    d1->GetXaxis()->SetTitleSize(0.07);
    d1->GetXaxis()->SetTitleOffset(0.88);
    d1->GetXaxis()->CenterTitle();
    d1->GetXaxis()->SetTickLength(0.03);
    d1->GetYaxis()->SetTitle("+ve sweighted MC/ MC Truth");
    d1->GetYaxis()->CenterTitle();
    d1->GetYaxis()->SetTitleSize(0.06);
    d1->GetYaxis()->SetLabelSize(0.06);
    d1->GetYaxis()->SetTitleOffset(0.70);
    d1->GetYaxis()->SetNdivisions(406);
//    gPad->SetTicky();
    d1->GetYaxis()->SetRangeUser(0.5,1.5);
    d1->Draw("LEP");  
    c2_var_name->Update();

    gSystem->ProcessEvents();

    string addr = "/mnt/t3nfs01/data01/shome/krgedia/Workspace/lxplus/splots_pt_80_100/splots/";
    string ud1 = "_1";
    string ud2 = "_2";
    string png = ".png";
    string compo;
    if (comps == " (signal)")
    {
        compo = "_sig";
    }
    else  
    {
        compo = "_bg";
    }

    TImage *img = TImage::Create();
    img->FromPad(c1_var_name);
    img->WriteImage((addr+variable_names+compo+ud1+png).c_str());

    TImage *img2 = TImage::Create();
    img2->FromPad(c2_var_name);
    img2->WriteImage((addr+variable_names+compo+ud2+png).c_str());

}


void sweights(RooRealVar var_name_b,RooRealVar var_name_n,RooDataSet* data_set_mc_r)

{
    string dataset_names = std::string(data_set_mc_r->GetName());
    string comps;
    if (dataset_names.find("r")==0)
    {
        comps = "Data sweights";
    }
    else
    {
        comps = "MC sweights";
    }

    TString dataset_name = dataset_names;
    TString comp = comps;

    TCanvas* canvas3 = new TCanvas(comp,comp,900,900);
    TPaveLabel* title3 = new TPaveLabel(0.09966592,0.9488429,0.9003341,0.9817296,comp,"BL");
    if (comps == "MC sweights")
    {
        title3->SetFillColor(40);
    }
    else
    {
        title3->SetFillColor(46);
    }
    title3->Draw();
    canvas3->Divide(1,2);
    canvas3->cd(1);
    gPad->SetPad(0,0.4811206,1,0.9317905);
    gPad->SetGridx();
    gPad->SetBottomMargin(0);
    TH1* f1 = (TH1*) data_set_mc_r->createHistogram("",var_name_b,Binning(100,-5,5));
    double bin_width = f1->GetXaxis()->GetBinWidth(0);
    std::string str = std::to_string(bin_width);
    str.erase (str.find_last_not_of('0') + 1, std::string::npos );
    TString jets = "jets/";
    
    if (comps == "MC sweights")
    {
        f1->SetLineColor(kBlue);
    }
    else
    {
        f1->SetLineColor(kRed);
    }
    f1->SetStats(0);
    f1->SetTitle("");
    f1->GetXaxis()->SetLabelOffset(999);
    f1->GetYaxis()->SetTitle(jets+str);
    f1->GetYaxis()->SetTitleOffset(0.8);
    f1->GetYaxis()->SetTitleSize(0.05);
    f1->GetYaxis()->SetLabelSize(0.03);

    TPaveText *myText1 = new TPaveText(0.2,0.7,0.4,0.85,"NDC");
    myText1->SetTextSize(0.04);
    myText1->SetFillColor(0); 
    myText1->SetTextFont(62);
    myText1->AddText("Signal");
    f1->Draw("HIST");
    myText1->Draw();

    canvas3->cd(2);
    gPad->SetPad(0.001113586,0.002436054,0.9988864,0.4835566);
    gPad->SetGridx();
    gPad->SetTopMargin(0);
    gPad->SetBottomMargin(0.1367089);
    TH1* g1 = (TH1*) data_set_mc_r->createHistogram("",var_name_n,Binning(100,-5,5));
    if (comps == "MC sweights")
    {
        g1->SetLineColor(kBlue);
    }
    else
    {
        g1->SetLineColor(kRed);
    }
    g1->SetStats(0);
    g1->SetTitle("");
    g1->GetXaxis()->SetTitle("sweights");
    g1->GetXaxis()->SetTitleSize(0.05);
    g1->GetXaxis()->SetTitleOffset(0.75);
    g1->GetYaxis()->SetTitle(jets+str);
    g1->GetYaxis()->SetTitleOffset(0.8);
    g1->GetYaxis()->SetTitleSize(0.05);
    g1->GetYaxis()->SetLabelSize(0.03);
    TPaveText *myText2 = new TPaveText(0.2,0.7,0.4,0.85,"NDC");
    myText2->SetTextSize(0.04);
    myText2->SetFillColor(0);
    myText2->SetTextFont(62);
    myText2->AddText("background");
    g1->Draw("HIST");
    myText2->Draw();

    canvas3->Draw();

    string addr = "/mnt/t3nfs01/data01/shome/krgedia/Workspace/lxplus/splots_pt_80_100/sweights/";
    string png = ".png";
    string compo;

    if (comps == "MC sweights")
    {
        compo = "MC_sweights";
    }
    else
    {
        compo = "Data_sweights";
    }

    gSystem->ProcessEvents();

    TImage *img = TImage::Create();
    img->FromPad(canvas3);
    img->WriteImage((addr+compo+png).c_str());

}


void sweights_bdt(RooRealVar var_name_b,RooRealVar var_name_n,RooRealVar var_name_bdt,RooDataSet* data_set_mc_r)

{
    string dataset_names = std::string(data_set_mc_r->GetName());
    string comps;
    if (dataset_names.find("r")==0)
    {
        comps = "BDT score VS sweights (DATA)";
    }
    else
    {
        comps = "BDT score VS sweights (MC)";
    }

    TString dataset_name = dataset_names;
    TString comp = comps;

    TCanvas* canvas4 = new TCanvas(comp,comp,900,900);
    TPaveLabel* title4 = new TPaveLabel(0.09966592,0.9488429,0.9003341,0.9817296,comp,"BL");
    if (comps == "BDT score VS sweights (MC)")
    {
        title4->SetFillColor(40);
    }
    else
    {
        title4->SetFillColor(46);
    }
    title4->Draw();
    canvas4->Divide(1,2);
    canvas4->cd(1);
    gPad->SetPad(0,0.4811206,1,0.9317905);
    gPad->SetGridx();
    gPad->SetBottomMargin(0);
    TH1* f1 = (TH1*) data_set_mc_r->createHistogram("",var_name_b,Binning(100,-5,5),YVar(var_name_bdt,Binning(100,-1,1)));
    double bin_width = f1->GetXaxis()->GetBinWidth(0);
    std::string str = std::to_string(bin_width);
    str.erase (str.find_last_not_of('0') + 1, std::string::npos );
    TString jets = "jets/";

    if (comps == "BDT score VS sweights (MC)")
    {
        f1->SetLineColor(kBlue);
    }
    else
    {
        f1->SetLineColor(kRed);
    }
    f1->SetStats(0);
    f1->SetTitle("");
    f1->GetXaxis()->SetLabelOffset(999);
    f1->GetYaxis()->SetTitle("BDT_score");
    f1->GetYaxis()->SetTitleOffset(0.8);
    f1->GetYaxis()->SetTitleSize(0.05);
    f1->GetYaxis()->SetLabelSize(0.03);

    TPaveText *myText1 = new TPaveText(0.2,0.7,0.4,0.85,"NDC");
    myText1->SetTextSize(0.04);
    myText1->SetFillColor(0);
    myText1->SetTextFont(62);
    myText1->AddText("Signal");
    f1->Draw("HIST");
    myText1->Draw();

    canvas4->cd(2);
    gPad->SetPad(0.001113586,0.002436054,0.9988864,0.4835566);
    gPad->SetGridx();
    gPad->SetTopMargin(0);
    gPad->SetBottomMargin(0.1367089);
    TH1* g1 = (TH1*) data_set_mc_r->createHistogram("",var_name_n,Binning(100,-5,5),YVar(var_name_bdt,Binning(100,-1,1)));
    if (comps == "BDT score VS sweights (MC)")
    {
        g1->SetLineColor(kBlue);
    }
    else
    {
        g1->SetLineColor(kRed);
    }
    g1->SetStats(0);
    g1->SetTitle("");
    g1->GetXaxis()->SetTitle("sweights");
    g1->GetXaxis()->SetTitleSize(0.05);
    g1->GetXaxis()->SetTitleOffset(0.75);
    g1->GetYaxis()->SetTitle("BDT_score");
    g1->GetYaxis()->SetTitleOffset(0.8);
    g1->GetYaxis()->SetTitleSize(0.05);
    g1->GetYaxis()->SetLabelSize(0.03);
    TPaveText *myText2 = new TPaveText(0.2,0.7,0.4,0.85,"NDC");
    myText2->SetTextSize(0.04);
    myText2->SetFillColor(0);
    myText2->SetTextFont(62);
    myText2->AddText("background");
    g1->Draw("HIST");
    myText2->Draw();

    canvas4->Draw();

    string addr = "/mnt/t3nfs01/data01/shome/krgedia/Workspace/lxplus/splots_pt_80_100/bdt_sweights/";
    string png = ".png";
    string compo;

    if (comps == "BDT score VS sweights (MC)")
    {
        compo = "BDT_MC_sweights";
    }
    else
    {
        compo = "BDT_Data_sweights";
    }

    gSystem->ProcessEvents();

    TImage *img = TImage::Create();
    img->FromPad(canvas4);
    img->WriteImage((addr+compo+png).c_str());

}




















/*
    TCanvas* canvas3 = new TCanvas("canvas3","canvas3",900,900);
    TPaveLabel* title3 = new TPaveLabel(0.1,0.962,0.9,0.995,"title","BL");
    title3->SetFillColor(kRed-7);
    title3->Draw();
    canvas3->Divide(1,2);
    canvas3->cd(1);
    gPad->SetGridx();
    gPad->SetBottomMargin(0);
    TH1* f1 = (TH1*) data_set->createHistogram("",yield_b_sw,Binning(100,-5,5));
    TH1* f2 = (TH1*) r_data_set->createHistogram("",r_yield_b_sw,Binning(100,-5,5));
    double scale4 = f1->GetSumOfWeights();
    f2->Scale(scale4/f2->GetSumOfWeights());
    cout<<" f1 "<<f1->GetSumOfWeights()<<endl;
    cout<<" f2 "<<f2->GetSumOfWeights()<<endl;
    f1->SetLineColor(kBlue);
    f2->SetLineColor(kRed);
    f1->SetLineWidth(2);
    f2->SetLineWidth(2);
    f1->SetLineStyle(1);
    f2->SetLineStyle(1);
    THStack* fs = new THStack("fs","");
    fs->Add(f1,"HIST");
    fs->Add(f2,"E0");
    fs->Draw("nostack");
    fs->GetHistogram()->GetXaxis()->SetLabelOffset(999);
    fs->GetHistogram()->GetYaxis()->SetTitle(jets+str);
    fs->GetHistogram()->GetYaxis()->SetTitleOffset(0.8);
    fs->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
    fs->GetHistogram()->GetYaxis()->SetLabelSize(0.03);
    TLegend* leg4 = new TLegend(0.65,0.79,0.80,0.94);
    leg4->AddEntry(h1,"sweighted MC","l");
    leg4->AddEntry(h2,"sweighted data","LEP");
    leg4->SetTextFont(12);
    leg4->SetTextSize(0.03);
    leg4->SetBorderSize(0);
    leg4->Draw("same");
    
    canvas3->cd(2);
    gPad->SetGridx();
    gPad->SetTopMargin(0);
    TH1* g1 = (TH1*) data_set->createHistogram("",yield_n_sw,Binning(100,-5,5));
    TH1* g2 = (TH1*) r_data_set->createHistogram("",r_yield_n_sw,Binning(100,-5,5));
    double scale5 = g1->GetSumOfWeights();
    g2->Scale(scale5/g2->GetSumOfWeights());
    cout<<" g1 "<<g1->GetSumOfWeights()<<endl;
    cout<<" g2 "<<g2->GetSumOfWeights()<<endl;
    g1->SetLineColor(kBlue);
    g2->SetLineColor(kRed);
    g1->SetLineWidth(2);
    g2->SetLineWidth(2);
    g1->SetLineStyle(1);
    g2->SetLineStyle(1);
    THStack* gs = new THStack("gs","");
    gs->Add(g1,"HIST");
    gs->Add(g2,"E0");
    gs->Draw("nostack");
    gs->GetHistogram()->GetYaxis()->SetTitle(jets+str);
    gs->GetHistogram()->GetYaxis()->SetTitleOffset(0.8);
    gs->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
    gs->GetHistogram()->GetYaxis()->SetLabelSize(0.03);
    TLegend* leg5 = new TLegend(0.65,0.79,0.80,0.94);
    leg5->AddEntry(g1,"sweighted MC","l");
    leg5->AddEntry(g2,"sweighted data","LEP");
    leg5->SetTextFont(12);
    leg5->SetTextSize(0.03);
    leg5->SetBorderSize(0);
    leg5->Draw("same");

    canvas3->Draw();




*/



















    void weight_plot(RooPlot* frame_var_name,RooRealVar var_name,double x_low,double x_high,double bins,TString title,RooDataSet* data_set,TCanvas* TCanvas_var_name,TLegend* leg_var_name)
{
    frame_var_name = var_name.frame(x_low, x_high, bins);
    frame_var_name->SetTitle(title);
    frame_var_name->GetXaxis()->CenterTitle();
    frame_var_name->GetXaxis()->SetTitleOffset(1.1);
    frame_var_name->GetYaxis()->SetTitleOffset(1.5);
    data_set->plotOn(frame_var_name,MarkerSize(0.7),MarkerColor(kBlack),LineColor(kBlack),Name("sweighted MC"));

    TCanvas_var_name = new TCanvas(title,title);
    leg_var_name = new TLegend(0.6,0.6,0.8,0.8);
    TObject* sMC_leg = (RooPlot*)frame_var_name->findObject("sweighted MC");
    leg_var_name->SetFillColor(kWhite);
    leg_var_name->SetTextSize(0.03);
    leg_var_name->AddEntry(sMC_leg,title,"LEP");
    leg_var_name->SetBorderSize(0);
    frame_var_name->Draw();
    leg_var_name->Draw();
    TCanvas_var_name->Draw();
}







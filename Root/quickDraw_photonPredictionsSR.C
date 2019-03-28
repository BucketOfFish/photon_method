#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>

#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPad.h"
#include "TCut.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"

using namespace std;

void quickDraw_photonPredictionsSR( string period = "data15-16" , string channel  = "mm" , string var = "MET", string smearing_mode = "NoSmear", string data = "zjets" ) {


  bool DF = TString(channel).EqualTo("em");

  //gStyle->SetOptStat(0);
  gStyle->SetOptStat(kFALSE);

  //-----------------------------------------------
  // define filenames
  //-----------------------------------------------

  // set up labels
  string mcdir      = "";
  string gdatalabel = "";

  if     ( TString(period).Contains("data15-16") ){
    mcdir    = "ZMC16a/";
    gdatalabel = "data15-16";
  }
  else if( TString(period).Contains("data17")    ){
    mcdir    = "ZMC16cd/";
    gdatalabel = "data17";
  }
  else if( TString(period).Contains("data18")    ){
    mcdir    = "ZMC16cd/";
    gdatalabel = "data18";
  }

  string outPath = "../OutputNtuples/v1.5/";

  string gfilename      = outPath + "gdata/" + period + "_merged_processed_" + channel + "_" + smearing_mode + ".root";
  string tt_filename    = outPath + mcdir + "ttbar_merged_processed.root";
  string vv_filename    = outPath + mcdir + "diboson_merged_processed.root";
  string z_filename     = outPath + mcdir + "Zjets_merged_processed.root";

  cout << "period               " << period        << endl;
  cout << "channel              " << channel       << endl;
  cout << "smearing mode        " << smearing_mode << endl;
  cout << "g filename           " << gfilename     << endl;

  TCut Zweight("totalWeight");
  TCut lumi1516("36100");
  TCut lumi17("44000");

  //TCut gselection("lep_pT[0]>25 && lep_pT[1]>25 && jet_n>=2  && bjet_n==0");
  //TCut VR("jet_n>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lep_pT[0]>25 && lep_pT[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && nJet20==2 && mll>80 && mll<100 && mjj<60.0 && mjj>100.");

  //-----------------------------------------------
  // Define selections
  //-----------------------------------------------

  TCut SR("jet_n>=2 && lep_pT[0]>25 && lep_pT[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && nBJet20_MV2c10_FixedCutBEff_77==0 && (mjj<60.0 || mjj>100.)");

  TCut weight_g    = "totalWeight";
  TCut weight_g_rw = "totalWeight*ptreweight3*ptreweight5"; //2-STEP

  //-----------------------------------------------
  // define and draw histograms
  //-----------------------------------------------

  const unsigned int nmetbins =  20;
  float metmax = 400;

  int   nbins =  20;
  float xmin  =   0;
  float xmax  = 400;

  string xtitle = var;

  if( TString(var).EqualTo("MET") ){
    xtitle = "E_{T}^{miss} [GeV]";
    nbins = nmetbins;
    xmin  = 0.0;
    xmax  = metmax;
  }

  double metbins[nmetbins+1] = {10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400};

  string metLoose = "MET_loose";
  string metTight = "MET";
  string metTighter = "MET_tighter";
  string metTenacious = "MET_tenacious";

  if( TString(data).EqualTo("zjets") ){

    TChain* chz = new TChain("BaselineTree");

    chz->Add( z_filename.c_str() );

    cout << "Z+jets entries, loose      " << chz->GetEntries()    << endl;

    TH1F* hz_loose     = new TH1F("hz_loose"     ,"",nmetbins,metbins)  ;
    TH1F* hz_tight     = new TH1F("hz_tight"     ,"",nmetbins,metbins)  ;
    TH1F* hz_tighter   = new TH1F("hz_tighter"     ,"",nmetbins,metbins)    ;
    TH1F* hz_tenacious = new TH1F("hz_tenacious"     ,"",nmetbins,metbins)   ;

    chz->  Draw("min(MET_loose,400)>>hz_loose", SR*Zweight      , "goff");
    chz->  Draw("min(MET,400)>>hz_tight", SR*Zweight      , "goff");
    chz->  Draw("min(MET_tighter,400)>>hz_tighter", SR*Zweight      , "goff");
    chz->  Draw("min(MET_tenacious,400)>>hz_tenacious", SR*Zweight      , "goff");

    cout << "MET total" << endl;
    cout << "Z+jets MC integral, loose   " << hz_loose->Integral()  << endl;
    cout << "Z+jets MC integral, tight       " << hz_tight->Integral()  << endl;
    cout << "Z+jets MC integral, tighter   " << hz_tighter->Integral()  << endl;
    cout << "Z+jets MC integral, tenacious       " << hz_tenacious->Integral()  << endl;

    cout << "MET 200-400" << endl;
    cout << "Z+jets MC integral, loose   " << hz_loose->Integral(11,20)  << endl;
    cout << "Z+jets MC integral, tight       " << hz_tight->Integral(11,20)  << endl;
    cout << "Z+jets MC integral, tighter   " << hz_tighter->Integral(11,20)  << endl;
    cout << "Z+jets MC integral, tenacious       " << hz_tenacious->Integral(11,20)  << endl;

    //-----------------------------------------------
    // make pretty plots
    //-----------------------------------------------

    // make canvas 

    TCanvas *can = new TCanvas("can","can",600,600);
    can->cd();

    TPad* mainpad = new TPad("mainpad","mainpad",0.05,0.0,1.0,0.85);
    mainpad->Draw();
    mainpad->cd();

    gPad->SetLogy();

    // MET loose
    hz_loose->SetLineColor(1);
    hz_loose->SetLineWidth(2);
    hz_loose->SetLineStyle(1);
    
    hz_loose->GetXaxis()->SetTitle(xtitle.c_str());
    hz_loose->GetYaxis()->SetTitle("entries / bin");
    hz_loose->Draw("hist");
    can->Update();

    // Other MET WP
    hz_tight->SetLineColor(4);
    hz_tight->SetLineWidth(2);
    hz_tight->SetLineStyle(1);
    hz_tight->Draw("samehist");

    hz_tighter->SetLineColor(3);
    hz_tighter->SetLineWidth(2);
    hz_tighter->SetLineStyle(1);
    hz_tighter->Draw("samehist");

    hz_tenacious->SetLineColor(2);
    hz_tenacious->SetLineWidth(2);
    hz_tenacious->SetLineStyle(1);
    hz_tenacious->Draw("samehist");

    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);

    leg->AddEntry(hz_loose,"Z+jets (from MC, MET loose)","f");
    leg->AddEntry(hz_tight,"Z+jets (from MC, MET tight)","f");
    leg->AddEntry(hz_tighter,"Z+jets (from MC, MET tighter)","f");
    leg->AddEntry(hz_tenacious,"Z+jets (from MC, MET tenacious)","f");

    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.03);
    tex->DrawLatex(0.6,0.65,"ATLAS Internal");
    if(TString(period).Contains("data15-16") ) tex->DrawLatex(0.6,0.61,"36 fb^{-1} 2015-2016 data");
    if(TString(period).Contains("data17")    ) tex->DrawLatex(0.6,0.61,"44 fb^{-1} 2017 data");
    if(TString(channel).Contains("ee")       ) tex->DrawLatex(0.6,0.57,"ee events");
    if(TString(channel).Contains("em")       ) tex->DrawLatex(0.6,0.57,"e#mu events");
    if(TString(channel).Contains("mm")       ) tex->DrawLatex(0.6,0.57,"#mu#mu events");

    // Tight/loose ratio

    TCanvas *can2 = new TCanvas("can2","can2",600,600);
    can2->cd();
    
    TPad* respad = new TPad("respad","respad",0.05,0.7,1.0,0.95);
    respad->Draw();
    respad->cd();
    respad->SetGridy();

    TH1F* hratio = (TH1F*) hz_tight->Clone("hratio");
    TH1F* hrLoose = (TH1F*) hz_loose->Clone("hrLoose");

    for( int ibin = 1 ; ibin <= hrLoose->GetXaxis()->GetNbins() ; ibin++ ) hrLoose->SetBinError(ibin,0.0);
    hratio->Divide(hrLoose);

    hratio->GetXaxis()->SetTitle("");
    hratio->GetXaxis()->SetLabelSize(0.);
    hratio->GetYaxis()->SetNdivisions(5);
    hratio->GetYaxis()->SetTitle("E^{tight}_{T,miss}/E^{loose}_{T,miss}");
    hratio->GetYaxis()->SetTitleSize(0.1);
    hratio->GetYaxis()->SetTitleOffset(0.5);
    hratio->GetYaxis()->SetLabelSize(0.1);
    hratio->SetMinimum(0.0);
    hratio->SetMaximum(2.0);
    hratio->GetYaxis()->SetRangeUser(0.,2.0);
    hratio->Draw("E1");

    can2->cd();

    // Tighter/Loose ratio
    
    TPad* respad1 = new TPad("respad1","respad1",0.05,0.35,1.0,0.6);
    respad1->Draw();
    respad1->cd();
    respad1->SetGridy();

    TH1F* hratio1 = (TH1F*) hz_tighter->Clone("hratio1");
    //TH1F* hrTighter = (TH1F*) hz_Loose->Clone("hrTighter");

    //for( int ibin = 1 ; ibin <= hrTighter->GetXaxis()->GetNbins() ; ibin++ ) hrTighter->SetBinError(ibin,0.0);
    hratio1->Divide(hrLoose);

    hratio1->GetXaxis()->SetTitle("");
    hratio1->GetXaxis()->SetLabelSize(0.);
    hratio1->GetYaxis()->SetNdivisions(5);
    hratio1->GetYaxis()->SetTitle("E^{tighter}_{T,miss}/E^{loose}_{T,miss}");
    hratio1->GetYaxis()->SetTitleSize(0.1);
    hratio1->GetYaxis()->SetTitleOffset(0.5);
    hratio1->GetYaxis()->SetLabelSize(0.1);
    hratio1->SetMinimum(0.0);
    hratio1->SetMaximum(2.0);
    hratio1->GetYaxis()->SetRangeUser(0.0,2.0);
    hratio1->Draw("E1");
    
    // Tenacious/Loose ratio
    can2->cd();

    TPad* respad2 = new TPad("respad2","respad2",0.05,0.05,1.0,0.3);
    respad2->Draw();
    respad2->cd();
    respad2->SetGridy();

    TH1F* hratio2 = (TH1F*) hz_tenacious->Clone("hratio2");
    //TH1F* hrTenacious = (TH1F*) hz_tenacious->Clone("hrTenacious");

    //for( int ibin = 1 ; ibin <= hrTenacious->GetXaxis()->GetNbins() ; ibin++ ) hrTenacious->SetBinError(ibin,0.0);
    hratio2->Divide(hrLoose);

    hratio2->GetXaxis()->SetTitle("");
    hratio2->GetXaxis()->SetLabelSize(0.);
    hratio2->GetYaxis()->SetNdivisions(5);
    hratio2->GetYaxis()->SetTitle("E^{tenacious}_{T,miss}/E^{loose}_{T,miss}");
    hratio2->GetYaxis()->SetTitleSize(0.1);
    hratio2->GetYaxis()->SetTitleOffset(0.5);
    hratio2->GetYaxis()->SetLabelSize(0.1);
    hratio2->SetMinimum(0.0);
    hratio2->SetMaximum(2.0);
    hratio2->GetYaxis()->SetRangeUser(0.,2.0);
    hratio2->Draw("E1");

    //can2->Update();
    
    can->Print(Form("plots/v1.5/VR_SR_studies/quickData_Data_%s_%s_%s_%s_SR_zmc_metWPs.pdf",period.c_str(),channel.c_str(),var.c_str(),smearing_mode.c_str()));

    can2->Print(Form("plots/v1.5/VR_SR_studies/ratioPlots_%s_%s_%s_%s_SR_zmc_metWPs.pdf",period.c_str(),channel.c_str(),var.c_str(),smearing_mode.c_str()));
  }

  
  else if( TString(data).EqualTo("gdata") ){

    TChain* gtree = new TChain("BaselineTree");

    gtree->Add( gfilename.c_str() );

    cout << "g+jets entries   " << gtree->GetEntries()    << endl;

    TH1F* hg_rw_loose     = new TH1F("hg_rw_loose"     ,"",nmetbins,metbins)  ;
    TH1F* hg_rw_tight     = new TH1F("hg_rw_tight"     ,"",nmetbins,metbins)  ;
    TH1F* hg_rw_tighter   = new TH1F("hg_rw_tighter"     ,"",nmetbins,metbins)    ;
    TH1F* hg_rw_tenacious = new TH1F("hg_rw_tenacious"     ,"",nmetbins,metbins)   ;

    gtree->  Draw("min(MET_loose,400)>>hg_rw_loose",         SR*weight_g_rw      , "goff");
    gtree->  Draw("min(MET,400)>>hg_rw_tight",         SR*weight_g_rw      , "goff");
    gtree->  Draw("min(MET_tighter,400)>>hg_rw_tighter",     SR*weight_g_rw  , "goff");
    gtree->  Draw("min(MET_tenacious,400)>>hg_rw_tenacious", SR*weight_g_rw , "goff");

    cout << "MET total" << endl;
    cout << "g+jets data integral, loose   "   << hg_rw_loose->Integral()  << endl;
    cout << "g+jets data integral, tight   "   << hg_rw_tight->Integral()  << endl;
    cout << "g+jets data integral, tighter "   << hg_rw_tighter->Integral()  << endl;
    cout << "g+jets data integral, tenacious " << hg_rw_tenacious->Integral()  << endl;

    cout << "MET 200-400" << endl;
    cout << "g+jets data integral, loose     " << hg_rw_loose->Integral(11,20)  << endl;
    cout << "g+jets data integral, tight     " << hg_rw_tight->Integral(11,20)  << endl;
    cout << "g+jets data integral, tighter   " << hg_rw_tighter->Integral(11,20)  << endl;
    cout << "g+jets data integral, tenacious " << hg_rw_tenacious->Integral(11,20)  << endl;

    //-----------------------------------------------
    // make pretty plots
    //-----------------------------------------------

    // make canvas

    TCanvas *can = new TCanvas("can","can",600,600);
    can->cd();

    TPad* mainpad = new TPad("mainpad","mainpad",0.05,0.0,1.0,0.85);
    mainpad->Draw();
    mainpad->cd();

    gPad->SetLogy();

    // MET loose
    hg_rw_loose->SetLineColor(1);
    hg_rw_loose->SetLineWidth(2);
    hg_rw_loose->SetLineStyle(1);

    hg_rw_loose->GetXaxis()->SetTitle(xtitle.c_str());
    hg_rw_loose->GetYaxis()->SetTitle("entries / bin");
    hg_rw_loose->Draw("hist");
    can->Update();

    // Other MET WP
    hg_rw_tight->SetLineColor(4);
    hg_rw_tight->SetLineWidth(2);
    hg_rw_tight->SetLineStyle(1);
    hg_rw_tight->Draw("samehist");

    hg_rw_tighter->SetLineColor(3);
    hg_rw_tighter->SetLineWidth(2);
    hg_rw_tighter->SetLineStyle(1);
    hg_rw_tighter->Draw("samehist");

    hg_rw_tenacious->SetLineColor(2);
    hg_rw_tenacious->SetLineWidth(2);
    hg_rw_tenacious->SetLineStyle(1);
    hg_rw_tenacious->Draw("samehist");

    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);

    leg->AddEntry(hg_rw_loose,"Z+jets (from #gamma+jets, MET loose)","f");
    leg->AddEntry(hg_rw_tight,"Z+jets (from #gamma+jets, MET tight)","f");
    leg->AddEntry(hg_rw_tighter,"Z+jets (from #gamma+jets, MET tighter)","f");
    leg->AddEntry(hg_rw_tenacious,"Z+jets (from #gamma+jets, MET tenacious)","f");

    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.03);
    tex->DrawLatex(0.6,0.65,"ATLAS Internal");
    if(TString(period).Contains("data15-16") ) tex->DrawLatex(0.6,0.61,"36 fb^{-1} 2015-2016 data");
    if(TString(period).Contains("data17")    ) tex->DrawLatex(0.6,0.61,"44 fb^{-1} 2017 data");
    if(TString(channel).Contains("ee")       ) tex->DrawLatex(0.6,0.57,"ee events");
    if(TString(channel).Contains("em")       ) tex->DrawLatex(0.6,0.57,"e#mu events");
    if(TString(channel).Contains("mm")       ) tex->DrawLatex(0.6,0.57,"#mu#mu events");

    // Tight/loose ratio

    TCanvas *can2 = new TCanvas("can2","can2",600,600);
    can2->cd();

    TPad* respad = new TPad("respad","respad",0.05,0.7,1.0,0.95);
    respad->Draw();
    respad->cd();
    respad->SetGridy();

    TH1F* hratio = (TH1F*) hg_rw_tight->Clone("hratio");
    TH1F* hrLoose = (TH1F*) hg_rw_loose->Clone("hrLoose");

    for( int ibin = 1 ; ibin <= hrLoose->GetXaxis()->GetNbins() ; ibin++ ) hrLoose->SetBinError(ibin,0.0);
    hratio->Divide(hrLoose);

    hratio->GetXaxis()->SetTitle("");
    hratio->GetXaxis()->SetLabelSize(0.);
    hratio->GetYaxis()->SetNdivisions(5);
    hratio->GetYaxis()->SetTitle("E^{tight}_{T,miss}/E^{loose}_{T,miss}");
    hratio->GetYaxis()->SetTitleSize(0.1);
    hratio->GetYaxis()->SetTitleOffset(0.5);
    hratio->GetYaxis()->SetLabelSize(0.1);
    hratio->SetMinimum(0.0);
    hratio->SetMaximum(2.0);
    hratio->GetYaxis()->SetRangeUser(0.,2.0);
    hratio->Draw("E1");

    can2->cd();

    // Tighter ratio

    TPad* respad1 = new TPad("respad1","respad1",0.05,0.35,1.0,0.6);
    respad1->Draw();
    respad1->cd();
    respad1->SetGridy();

    TH1F* hratio1 = (TH1F*) hg_rw_tighter->Clone("hratio1");
    //TH1F* hrTighter = (TH1F*) hz_Loose->Clone("hrTighter");

    //for( int ibin = 1 ; ibin <= hrTighter->GetXaxis()->GetNbins() ; ibin++ ) hrTighter->SetBinError(ibin,0.0);
    hratio1->Divide(hrLoose);

    hratio1->GetXaxis()->SetTitle("");
    hratio1->GetXaxis()->SetLabelSize(0.);
    hratio1->GetYaxis()->SetNdivisions(5);
    hratio1->GetYaxis()->SetTitle("E^{tighter}_{T,miss}/E^{loose}_{T,miss}");
    hratio1->GetYaxis()->SetTitleSize(0.1);
    hratio1->GetYaxis()->SetTitleOffset(0.5);
    hratio1->GetYaxis()->SetLabelSize(0.1);
    hratio1->SetMinimum(0.0);
    hratio1->SetMaximum(2.0);
    hratio1->GetYaxis()->SetRangeUser(0.0,2.0);
    hratio1->Draw("E1");

    can2->cd();

    //Tenacious ratio
    TPad* respad2 = new TPad("respad2","respad2",0.05,0.05,1.0,0.3);
    respad2->Draw();
    respad2->cd();
    respad2->SetGridy();

    TH1F* hratio2 = (TH1F*) hg_rw_tenacious->Clone("hratio2");
    //TH1F* hrTenacious = (TH1F*) hz_tenacious->Clone("hrTenacious");

    //for( int ibin = 1 ; ibin <= hrTenacious->GetXaxis()->GetNbins() ; ibin++ ) hrTenacious->SetBinError(ibin,0.0);
    hratio2->Divide(hrLoose);

    hratio2->GetXaxis()->SetTitle("");
    hratio2->GetXaxis()->SetLabelSize(0.);
    hratio2->GetYaxis()->SetNdivisions(5);
    hratio2->GetYaxis()->SetTitle("E^{tenacious}_{T,miss}/E^{loose}_{T,miss}");
    hratio2->GetYaxis()->SetTitleSize(0.1);
    hratio2->GetYaxis()->SetTitleOffset(0.5);
    hratio2->GetYaxis()->SetLabelSize(0.1);
    hratio2->SetMinimum(0.0);
    hratio2->SetMaximum(2.0);
    hratio2->GetYaxis()->SetRangeUser(0.,2.0);
    hratio2->Draw("E1");

    //can2->Update();

    can2->Print(Form("plots/v1.5/VR_SR_studies/ratioPlots_%s_%s_%s_%s_SR_gdata_metWPs.pdf",period.c_str(),channel.c_str(),var.c_str(),smearing_mode.c_str()));

    can->Print(Form("plots/v1.5/VR_SR_studies/quickDraw_%s_%s_%s_%s_SR_gdata_metWPs.pdf",period.c_str(),channel.c_str(),var.c_str(),smearing_mode.c_str()));
  }

  else if( TString(data).EqualTo("ttbar") ){

    TChain* chtt = new TChain("BaselineTree");

    chtt->Add( tt_filename.c_str() );

    cout << "ttbar entries,      " << chtt->GetEntries()    << endl;

    TH1F* htt_loose     = new TH1F("htt_loose"     ,"",nmetbins,metbins)  ;
    TH1F* htt_tight     = new TH1F("htt_tight"     ,"",nmetbins,metbins)  ;
    TH1F* htt_tighter   = new TH1F("htt_tighter"     ,"",nmetbins,metbins)    ;
    TH1F* htt_tenacious = new TH1F("htt_tenacious"     ,"",nmetbins,metbins)   ;

    chtt->  Draw("min(MET_loose,400)>>htt_loose", SR*Zweight      , "goff");
    chtt->  Draw("min(MET,400)>>htt_tight", SR*Zweight      , "goff");
    chtt->  Draw("min(MET_tighter,400)>>htt_tighter", SR*Zweight      , "goff");
    chtt->  Draw("min(MET_tenacious,400)>>htt_tenacious", SR*Zweight      , "goff");

    cout << "MET total" << endl;
    cout << "ttbar MC integral, loose     " << htt_loose->Integral()  << endl;
    cout << "ttbar MC integral, tight     " << htt_tight->Integral()  << endl;
    cout << "ttbar MC integral, tighter   " << htt_tighter->Integral()  << endl;
    cout << "ttbar MC integral, tenacious " << htt_tenacious->Integral()  << endl;

    cout << "MET 200-400" << endl;
    cout << "ttbar MC integral, loose    " << htt_loose->Integral(11,20)  << endl;
    cout << "ttbar MC integral, tight    " << htt_tight->Integral(11,20)  << endl;
    cout << "ttbar MC integral, tighter  " << htt_tighter->Integral(11,20)  << endl;
    cout << "ttbar MC integral, tenacious" << htt_tenacious->Integral(11,20)  << endl;

    //-----------------------------------------------
    // make pretty plots
    //-----------------------------------------------

    // make canvas

    TCanvas *can = new TCanvas("can","can",600,600);
    can->cd();

    TPad* mainpad = new TPad("mainpad","mainpad",0.05,0.0,1.0,0.85);
    mainpad->Draw();
    mainpad->cd();

    gPad->SetLogy();

    // MET loose
    htt_loose->SetLineColor(1);
    htt_loose->SetLineWidth(2);
    htt_loose->SetLineStyle(1);

    htt_loose->GetXaxis()->SetTitle(xtitle.c_str());
    htt_loose->GetYaxis()->SetTitle("entries / bin");
    htt_loose->Draw("hist");
    can->Update();

    // Other MET WP
    htt_tight->SetLineColor(4);
    htt_tight->SetLineWidth(2);
    htt_tight->SetLineStyle(1);
    htt_tight->Draw("samehist");

    htt_tighter->SetLineColor(3);
    htt_tighter->SetLineWidth(2);
    htt_tighter->SetLineStyle(1);
    htt_tighter->Draw("samehist");

    htt_tenacious->SetLineColor(2);
    htt_tenacious->SetLineWidth(2);
    htt_tenacious->SetLineStyle(1);
    htt_tenacious->Draw("samehist");

    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);

    leg->AddEntry(htt_loose,"t#bar{t} (from MC, MET loose)","f");
    leg->AddEntry(htt_tight,"t#bar{t} (from MC, MET tight)","f");
    leg->AddEntry(htt_tighter,"t#bar{t} (from MC, MET tighter)","f");
    leg->AddEntry(htt_tenacious,"t#bar{t} (from MC, MET tenacious)","f");

    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.03);
    tex->DrawLatex(0.6,0.65,"ATLAS Internal");
    if(TString(period).Contains("data15-16") ) tex->DrawLatex(0.6,0.61,"36 fb^{-1} 2015-2016 data");
    if(TString(period).Contains("data17")    ) tex->DrawLatex(0.6,0.61,"44 fb^{-1} 2017 data");
    if(TString(channel).Contains("ee")       ) tex->DrawLatex(0.6,0.57,"ee events");
    if(TString(channel).Contains("em")       ) tex->DrawLatex(0.6,0.57,"e#mu events");
    if(TString(channel).Contains("mm")       ) tex->DrawLatex(0.6,0.57,"#mu#mu events");


    // Tight/loose ratio

    TCanvas *can2 = new TCanvas("can2","can2",600,600);
    can2->cd();

    TPad* respad = new TPad("respad","respad",0.05,0.7,1.0,0.95);
    respad->Draw();
    respad->cd();
    respad->SetGridy();

    TH1F* hratio = (TH1F*) htt_tight->Clone("hratio");
    TH1F* hrLoose = (TH1F*) htt_loose->Clone("hrLoose");

    for( int ibin = 1 ; ibin <= hrLoose->GetXaxis()->GetNbins() ; ibin++ ) hrLoose->SetBinError(ibin,0.0);
    hratio->Divide(hrLoose);

    hratio->GetXaxis()->SetTitle("");
    hratio->GetXaxis()->SetLabelSize(0.);
    hratio->GetYaxis()->SetNdivisions(5);
    hratio->GetYaxis()->SetTitle("E^{tight}_{T,miss}/E^{loose}_{T,miss}");
    hratio->GetYaxis()->SetTitleSize(0.1);
    hratio->GetYaxis()->SetTitleOffset(0.5);
    hratio->GetYaxis()->SetLabelSize(0.1);
    hratio->SetMinimum(0.0);
    hratio->SetMaximum(2.0);
    hratio->GetYaxis()->SetRangeUser(0.,2.0);
    hratio->Draw("E1");

    can2->cd();

    // Tighter ratio

    TPad* respad1 = new TPad("respad1","respad1",0.05,0.35,1.0,0.6);
    respad1->Draw();
    respad1->cd();
    respad1->SetGridy();

    TH1F* hratio1 = (TH1F*) htt_tighter->Clone("hratio1");
    //TH1F* hrTighter = (TH1F*) hz_Loose->Clone("hrTighter");

    //for( int ibin = 1 ; ibin <= hrTighter->GetXaxis()->GetNbins() ; ibin++ ) hrTighter->SetBinError(ibin,0.0);
    hratio1->Divide(hrLoose);

    hratio1->GetXaxis()->SetTitle("");
    hratio1->GetXaxis()->SetLabelSize(0.);
    hratio1->GetYaxis()->SetNdivisions(5);
    hratio1->GetYaxis()->SetTitle("E^{tighter}_{T,miss}/E^{loose}_{T,miss}");
    hratio1->GetYaxis()->SetTitleSize(0.1);
    hratio1->GetYaxis()->SetTitleOffset(0.5);
    hratio1->GetYaxis()->SetLabelSize(0.1);
    hratio1->SetMinimum(0.0);
    hratio1->SetMaximum(2.0);
    hratio1->GetYaxis()->SetRangeUser(0.0,2.0);
    hratio1->Draw("E1");

    can2->cd();

    //Tenacious ratio
    TPad* respad2 = new TPad("respad2","respad2",0.05,0.05,1.0,0.3);
    respad2->Draw();
    respad2->cd();
    respad2->SetGridy();

    TH1F* hratio2 = (TH1F*) htt_tenacious->Clone("hratio2");
    //TH1F* hrTenacious = (TH1F*) hz_tenacious->Clone("hrTenacious");

    //for( int ibin = 1 ; ibin <= hrTenacious->GetXaxis()->GetNbins() ; ibin++ ) hrTenacious->SetBinError(ibin,0.0);
    hratio2->Divide(hrLoose);

    hratio2->GetXaxis()->SetTitle("");
    hratio2->GetXaxis()->SetLabelSize(0.);
    hratio2->GetYaxis()->SetNdivisions(5);
    hratio2->GetYaxis()->SetTitle("E^{tenacious}_{T,miss}/E^{loose}_{T,miss}");
    hratio2->GetYaxis()->SetTitleSize(0.1);
    hratio2->GetYaxis()->SetTitleOffset(0.5);
    hratio2->GetYaxis()->SetLabelSize(0.1);
    hratio2->SetMinimum(0.0);
    hratio2->SetMaximum(2.0);
    hratio2->GetYaxis()->SetRangeUser(0.,2.0);
    hratio2->Draw("E1");

    //can2->Update();

    can2->Print(Form("plots/v1.5/VR_SR_studies/ratioPlots_%s_%s_%s_%s_SR_ttbar_metWPs.pdf",period.c_str(),channel.c_str(),var.c_str(),smearing_mode.c_str()));

    can->Print(Form("plots/v1.5/VR_SR_studies/quickDraw_%s_%s_%s_%s_SR_ttbar_metWPs.pdf",period.c_str(),channel.c_str(),var.c_str(),smearing_mode.c_str()));
  }

  else if( TString(data).EqualTo("vv") ){

    TChain* chvv = new TChain("BaselineTree");

    chvv->Add( vv_filename.c_str() );

    cout << "VV entries, loose      " << chvv->GetEntries()    << endl;

    TH1F* hvv_loose     = new TH1F("hvv_loose"     ,"",nmetbins,metbins)  ;
    TH1F* hvv_tight     = new TH1F("hvv_tight"     ,"",nmetbins,metbins)  ;
    TH1F* hvv_tighter   = new TH1F("hvv_tighter"     ,"",nmetbins,metbins)    ;
    TH1F* hvv_tenacious = new TH1F("hvv_tenacious"     ,"",nmetbins,metbins)   ;

    chvv->  Draw("min(MET_loose,400)>>hvv_loose", SR*Zweight      , "goff");
    chvv->  Draw("min(MET,400)>>hvv_tight", SR*Zweight      , "goff");
    chvv->  Draw("min(MET_tighter,400)>>hvv_tighter", SR*Zweight      , "goff");
    chvv->  Draw("min(MET_tenacious,400)>>hvv_tenacious", SR*Zweight      , "goff");

    cout << "MET total" << endl;
    cout << "VV MC integral, loose     " << hvv_loose->Integral()  << endl;
    cout << "VV MC integral, tight     " << hvv_tight->Integral()  << endl;
    cout << "VV MC integral, tighter   " << hvv_tighter->Integral()  << endl;
    cout << "VV MC integral, tenacious " << hvv_tenacious->Integral()  << endl;

    cout << "MET 200-400" << endl;
    cout << "VV MC integral, loose     " << hvv_loose->Integral(11,20)  << endl;
    cout << "VV MC integral, tight     " << hvv_tight->Integral(11,20)  << endl;
    cout << "VV MC integral, tighter   " << hvv_tighter->Integral(11,20)  << endl;
    cout << "VV MC integral, tenacious " << hvv_tenacious->Integral(11,20)  << endl;

    //-----------------------------------------------
    // make pretty plots
    //-----------------------------------------------

    // make canvas

    TCanvas *can = new TCanvas("can","can",600,600);
    can->cd();

    TPad* mainpad = new TPad("mainpad","mainpad",0.05,0.0,1.0,0.85);
    mainpad->Draw();
    mainpad->cd();

    gPad->SetLogy();

    // MET loose
    hvv_loose->SetLineColor(1);
    hvv_loose->SetLineWidth(2);
    hvv_loose->SetLineStyle(1);

    hvv_loose->GetXaxis()->SetTitle(xtitle.c_str());
    hvv_loose->GetYaxis()->SetTitle("entries / bin");
    hvv_loose->Draw("hist");
    can->Update();

    // Other MET WP
    hvv_tight->SetLineColor(4);
    hvv_tight->SetLineWidth(2);
    hvv_tight->SetLineStyle(1);
    hvv_tight->Draw("samehist");

    hvv_tighter->SetLineColor(3);
    hvv_tighter->SetLineWidth(2);
    hvv_tighter->SetLineStyle(1);
    hvv_tighter->Draw("samehist");

    hvv_tenacious->SetLineColor(2);
    hvv_tenacious->SetLineWidth(2);
    hvv_tenacious->SetLineStyle(1);
    hvv_tenacious->Draw("samehist");

    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);

    leg->AddEntry(hvv_loose,"VV (from MC, MET loose)","f");
    leg->AddEntry(hvv_tight,"VV (from MC, MET tight)","f");
    leg->AddEntry(hvv_tighter,"VV (from MC, MET tighter)","f");
    leg->AddEntry(hvv_tenacious,"VV (from MC, MET tenacious)","f");

    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.03);
    tex->DrawLatex(0.6,0.65,"ATLAS Internal");
    if(TString(period).Contains("data15-16") ) tex->DrawLatex(0.6,0.61,"36 fb^{-1} 2015-2016 data");
    if(TString(period).Contains("data17")    ) tex->DrawLatex(0.6,0.61,"44 fb^{-1} 2017 data");
    if(TString(channel).Contains("ee")       ) tex->DrawLatex(0.6,0.57,"ee events");
    if(TString(channel).Contains("em")       ) tex->DrawLatex(0.6,0.57,"e#mu events");
    if(TString(channel).Contains("mm")       ) tex->DrawLatex(0.6,0.57,"#mu#mu events");

    // Tight ratio

    TCanvas *can2 = new TCanvas("can2","can2",600,600);
    can2->cd();

    TPad* respad = new TPad("respad","respad",0.05,0.7,1.0,0.95);
    respad->Draw();
    respad->cd();
    respad->SetGridy();

    TH1F* hratio = (TH1F*) hvv_tight->Clone("hratio");
    TH1F* hrLoose = (TH1F*) hvv_loose->Clone("hrLoose");

    for( int ibin = 1 ; ibin <= hrLoose->GetXaxis()->GetNbins() ; ibin++ ) hrLoose->SetBinError(ibin,0.0);
    hratio->Divide(hrLoose);

    hratio->GetXaxis()->SetTitle("");
    hratio->GetXaxis()->SetLabelSize(0.);
    hratio->GetYaxis()->SetNdivisions(5);
    hratio->GetYaxis()->SetTitle("E^{tight}_{T,miss}/E^{loose}_{T,miss}");
    hratio->GetYaxis()->SetTitleSize(0.1);
    hratio->GetYaxis()->SetTitleOffset(0.5);
    hratio->GetYaxis()->SetLabelSize(0.1);
    hratio->SetMinimum(0.0);
    hratio->SetMaximum(2.0);
    hratio->GetYaxis()->SetRangeUser(0.,2.0);
    hratio->Draw("E1");

    can2->cd();

    // Tighter ratio

    TPad* respad1 = new TPad("respad1","respad1",0.05,0.35,1.0,0.6);
    respad1->Draw();
    respad1->cd();
    respad1->SetGridy();

    TH1F* hratio1 = (TH1F*) hvv_tighter->Clone("hratio1");
    //TH1F* hrTighter = (TH1F*) hz_Loose->Clone("hrTighter");

    //for( int ibin = 1 ; ibin <= hrTighter->GetXaxis()->GetNbins() ; ibin++ ) hrTighter->SetBinError(ibin,0.0);
    hratio1->Divide(hrLoose);

    hratio1->GetXaxis()->SetTitle("");
    hratio1->GetXaxis()->SetLabelSize(0.);
    hratio1->GetYaxis()->SetNdivisions(5);
    hratio1->GetYaxis()->SetTitle("E^{tighter}_{T,miss}/E^{loose}_{T,miss}");
    hratio1->GetYaxis()->SetTitleSize(0.1);
    hratio1->GetYaxis()->SetTitleOffset(0.5);
    hratio1->GetYaxis()->SetLabelSize(0.1);
    hratio1->SetMinimum(0.0);
    hratio1->SetMaximum(2.0);
    hratio1->GetYaxis()->SetRangeUser(0.0,2.0);
    hratio1->Draw("E1");

    can2->cd();

    //Tenacious ratio
    TPad* respad2 = new TPad("respad2","respad2",0.05,0.05,1.0,0.3);
    respad2->Draw();
    respad2->cd();
    respad2->SetGridy();

    TH1F* hratio2 = (TH1F*) hvv_tenacious->Clone("hratio2");
    //TH1F* hrTenacious = (TH1F*) hz_tenacious->Clone("hrTenacious");

    //for( int ibin = 1 ; ibin <= hrTenacious->GetXaxis()->GetNbins() ; ibin++ ) hrTenacious->SetBinError(ibin,0.0);
    hratio2->Divide(hrLoose);

    hratio2->GetXaxis()->SetTitle("");
    hratio2->GetXaxis()->SetLabelSize(0.);
    hratio2->GetYaxis()->SetNdivisions(5);
    hratio2->GetYaxis()->SetTitle("E^{tenacious}_{T,miss}/E^{loose}_{T,miss}");
    hratio2->GetYaxis()->SetTitleSize(0.1);
    hratio2->GetYaxis()->SetTitleOffset(0.5);
    hratio2->GetYaxis()->SetLabelSize(0.1);
    hratio2->SetMinimum(0.0);
    hratio2->SetMaximum(2.0);
    hratio2->GetYaxis()->SetRangeUser(0.,2.0);
    hratio2->Draw("E1");

    can2->Print(Form("plots/v1.5/VR_SR_studies/ratioPlots_%s_%s_%s_%s_SR_VV_metWPs.pdf",period.c_str(),channel.c_str(),var.c_str(),smearing_mode.c_str()));
    can->Print(Form("plots/v1.5/VR_SR_studies/quickDraw_%s_%s_%s_%s_SR_VV_metWPs.pdf",period.c_str(),channel.c_str(),var.c_str(),smearing_mode.c_str()));
  }

}

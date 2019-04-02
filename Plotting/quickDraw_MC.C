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

#include "../BasicSetting.C"

using namespace std;

void quickDraw_MC(string period = "data15-16" , string channel  = "mm" , string var = "HT", string smearing_mode = "NoSmear" ) {
  
  bool DF = TString(channel).EqualTo("em");

  bool normalize = true;
  if ((TString(var).Contains("pt")) || (DF)) normalize = false;
  
  gStyle->SetOptStat(0);

  //-----------------------------------------------
  // define filenames
  //-----------------------------------------------

  string Zfilename      = outputPath + "ZMC16a/Zjets_merged_processed.root";
  string gfilename      = outputPath + "gmc/gmc_" + channel + "_" + smearing_mode + ".root";
  
  cout << "channel              " << channel       << endl;
  cout << "smearing mode        " << smearing_mode << endl;
  cout << "Z datafilename       " << Zfilename     << endl;
  cout << "g filename           " << gfilename     << endl;
  cout << "DF?                  " << DF            << endl;
  cout << "" << endl;
  
  //-----------------------------------------------
  // add files to TChain
  //-----------------------------------------------

  TChain * gtree = new TChain("BaselineTree");
  if( !DF ) gtree->Add( gfilename.c_str() );
	
  TChain * Ztree = new TChain("BaselineTree");
  Ztree->Add( Zfilename.c_str() );

  cout << "g entries            " << gtree->GetEntries()  << endl;
  cout << "Z data entries       " << Ztree->GetEntries()  << endl;

  //-----------------------------------------------
  // Define selections
  //-----------------------------------------------
	
  TCut Zselection("mll>81 && mll<101 && jet_n >= 2 && is_OS && lep_pT[0]>25.0 && lep_pT[1]>25.0 && bjet_n == 0");
  //TCut Zweight("totalWeight");
  //TCut lumi1516("36100");
  //TCut lumi17("44000");
  
  TCut gselection("jet_n>=2 && lep_pT[0]>25 && lep_pT[1]>25 && bjet_n == 0");

  TCut CR("MET<60.0");

  TCut ee("channel==1");
  TCut mm("channel==0");
  TCut em("channel==2 || channel==3");
  if     ( TString(channel).EqualTo("ee") ) Zselection += ee;
  else if( TString(channel).EqualTo("mm") ) Zselection += mm;
  else if( TString(channel).EqualTo("em") ) Zselection += em;
  else{
    cout << "Unrecognized channel! quitting   " << channel << endl;
    exit(0);
  }

  //-----------------------------------------------
  // Set weights 
  //-----------------------------------------------
  TCut lumi("1.0");
  if( TString(period).EqualTo("data15-16") ) lumi = TCut("36100");
  if( TString(period).EqualTo("data17")    ) lumi = TCut("44100");
  if( TString(period).EqualTo("data18")    ) lumi = TCut("64610");

  TCut Zweight = "totalWeight*36100";//change luminosity according to data year
  TCut weight_g    = "totalWeight*36100";
  //TCut weight_g_rw = "totalWeight*ptreweight9*36100";//For checking step1 of reweighting
  TCut weight_g_rw = "totalWeight*ptreweight9*ptreweight10*36100";

  cout << "Z selection          " << Zselection.GetTitle()  << endl;  
  cout << "Z weight             " << Zweight.GetTitle()     << endl;
  cout << "g selection          " << gselection.GetTitle()  << endl;
  cout << "g weight             " << weight_g.GetTitle()    << endl;
  cout << "g weight (reweight)  " << weight_g_rw.GetTitle() << endl;
  cout << "luminosity           " << lumi.GetTitle()       << endl;

  //-----------------------------------------------
  // define and draw histograms
  //-----------------------------------------------
  
  TH1F* hZ    ;//= new TH1F();    
  TH1F* hg    ;//= new TH1F();    
  TH1F* hg_rw ;//= new TH1F(); 

  string xtitle = var;

  int   nbins =  20;
  float xmin  =   0;
  float xmax  = 100;

  if( TString(var).EqualTo("MET") ){
    xtitle = "E_{T}^{miss} [GeV]";
    nbins = 20;
    xmin  = 0.0;
    xmax  = 400;
  }
  
  else if( TString(var).EqualTo("METl") || TString(var).EqualTo("METt") ){
    if( TString(var).EqualTo("METl") ) xtitle = "E_{T,||}^{miss} [GeV]";
    if( TString(var).EqualTo("METt") ) xtitle = "E_{T,#perp}^{miss} [GeV]";
    nbins =   20;
    xmin  = -200;
    xmax  =  200;
  }

  else if( TString(var).EqualTo("Z_pt") ){
    xtitle = "p_{T} [GeV]";
  }

  else if( TString(var).EqualTo("jet_n") ){
    xtitle = "n_{jets}";
    nbins = 6;
    xmin  = 2;
    xmax  = 8;
  }
  
  else if( TString(var).EqualTo("bjet_n") ){
    xtitle = "n_{b-jets}";
    nbins = 4;
    xmin  = 0;
    xmax  = 4;
  }

  else if( TString(var).EqualTo("HT") ){
    xtitle = "H_{T}";
    //nbins =   20;
    //xmin  =    0;
    //xmax  = 1000;
  }

  else if( TString(var).EqualTo("mll") ){
    xtitle = "m_{ll} [GeV]";
    nbins =   30;
    xmin  =    0;
    xmax  =  300;
  }

  else if( TString(var).EqualTo("MT2W") ){
    xtitle = "m_{T2}^{W} [GeV]";
    nbins =   20;
    xmin  =    0;
    xmax  =  200;
  }

  else if( TString(var).EqualTo("lep_pT[0]") ){
    xtitle = "1^{st} lepton p_{T} [GeV]";
    nbins =   20;
    xmin  =    0;
    xmax  =  200;
  }

  else if( TString(var).EqualTo("lep_pT[1]") ){
    xtitle = "2^{nd} lepton p_{T} [GeV]";
    nbins =   20;
    xmin  =    0;
    xmax  =  100;
  }

  else if( TString(var).EqualTo("DPhi_METJetLeading") ){
    xtitle = "#Delta#phi(jet_{1},E_{T}^{miss})";
    nbins =   20;
    xmin  =    0;
    xmax  =  3.14;
  }

  else if( TString(var).EqualTo("DPhi_METJetSecond") ){
    xtitle = "#Delta#phi(jet_{2},E_{T}^{miss})";
    nbins =   20;
    xmin  =    0;
    xmax  =  3.14;
  }

  else{
    cout << "Error! unrecognized variable, need to set binning, quitting! " << var << endl;
    exit(0);
  }
  
  const unsigned int nZptbins = 16;
  double Zptbins[nZptbins+1] = {40, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500, 600, 700, 850, 1000};

  const unsigned int nHTbins = 16;
  double HTbins[nHTbins+1] = {40, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500, 600, \
700, 850, 1000};
  // Initialize histograms

  if( TString(var).EqualTo("Z_pt") ){
    hZ     = new TH1F("hZ"     , "" , nZptbins , Zptbins );
    //hg     = new TH1F("hg"     , "" , nZptbins , Zptbins );
    hg_rw  = new TH1F("hg_rw"  , "" , nZptbins , Zptbins );
  }

  else if( TString(var).EqualTo("HT") ){
    hZ     = new TH1F("hZ"     , "" , nHTbins , HTbins );
    //hg     = new TH1F("hg"     , "" , nZptbins , Zptbins );
    hg_rw  = new TH1F("hg_rw"  , "" , nHTbins , HTbins );
  }

  else{
    hZ     = new TH1F("hZ"     , "" , nbins , xmin , xmax );
    //hg     = new TH1F("hg"     , "" , nbins , xmin , xmax );
    hg_rw  = new TH1F("hg_rw"  , "" , nbins , xmin , xmax );
  }

  Ztree->Draw(Form("%s>>hZ",var.c_str())       , Zselection*Zweight      , "goff");

  if( !DF ){
    //gtree->Draw(Form("%s>>hg",var.c_str())     , gselection*weight_g     , "goff");
    gtree->Draw(Form("%s>>hg_rw",var.c_str())  , gselection*weight_g_rw  , "goff");
  }
  
  //DEBUGGING
  /*
  std::cout << "CHECKPOINT: Checking histogram bins after reweighting" << std::endl;
  std::cout << "reweighted Z" << std::endl;
  for (int i=0; i<=nZptbins; i++)
    std::cout << hZ->GetBinContent(i) << std::endl;
  std::cout << "reweighted G" << std::endl;
  for (int i=0; i<=nZptbins; i++)
    std::cout << hg_rw->GetBinContent(i) << std::endl;
  std::cout << "CHECKPOINT: Finished checking histogram bins after reweighting" << std::endl;
  */
  cout << "" << endl;
  cout << "Z MC integral      " << hZ->Integral()  << endl;
  //cout << "g MC raw integral  " << hg->Integral()  << endl;
  cout << "g MC rw integral   " << hg_rw->Integral()  << endl;
  cout << "" << endl;


  //-----------------------------------------------
  // make pretty plots
  //-----------------------------------------------

  cout << "MET100-150" << endl;
  cout << "Z                      " << hZ->Integral(11,15)     << endl;
  cout << "g MC (reweighted)    " << hg_rw->Integral(11,15)  << endl;
  //cout << "g MC (raw)           " << hg->Integral(11,15)     << endl;

  cout << "MET150-200" << endl;
  cout << "Z                      " << hZ->Integral(16,21)     << endl;
  cout << "g MC (reweighted)    " << hg_rw->Integral(16,21)  << endl;
  //cout << "g MC (raw)           " << hg->Integral(16,21)     << endl;

  //-------------------------------------------
  //--------- Normalize histograms ------------ 
  //-------------------------------------------
  
  float integralG; 
  float integralZ; 
  float normalization; 
  /*
  integralG = hg_rw->Integral(1,3);
  integralZ = hZ->Integral(1,3);
  //IF USING HT VAR 
  //integralG = hg_rw->Integral(10,15);
  //integralZ = hZ->Integral(10,15);
  normalization = integralZ/integralG;
  hg_rw->Scale(normalization);
  */

  TCanvas *can = new TCanvas("can","can",600,600);
  can->cd();

  TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);
  mainpad->Draw();
  mainpad->cd();
  
  gPad->SetLogy();
  
  hZ->SetLineColor(1);
  hZ->SetLineWidth(2);
  hZ->SetMarkerStyle(20);

  hZ->GetXaxis()->SetTitle(xtitle.c_str());
  hZ->GetYaxis()->SetTitle("entries / bin");
  hZ->Draw("E1");
  
  //hg->SetLineColor(2);
  hg_rw->SetLineColor(4);

  THStack *mcstack = new THStack("mcstack","mcstack");
  if( !DF ) mcstack->Add(hg_rw);
  mcstack->Draw("samehist");
  hZ->Draw("sameE1");
  //if( !DF ) hg->Draw("samehist");
  hZ->Draw("axissame");

  TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
  //leg->AddEntry(hZ,"data","lp");
  if( !DF ){
    leg->AddEntry(hg_rw,"Z+jets (from #gamma+jets, reweighted)","f");
    //leg->AddEntry(hg,"Z+jets (from #gamma+jets, raw)","f");
    leg->AddEntry(hZ,"Z+jets (from MC)","f");
  }
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.03);
  tex->DrawLatex(0.6,0.65,"ATLAS Internal");
  if(TString(channel).Contains("ee")       ) tex->DrawLatex(0.6,0.57,"ee events");
  if(TString(channel).Contains("em")       ) tex->DrawLatex(0.6,0.57,"e#mu events");
  if(TString(channel).Contains("mm")       ) tex->DrawLatex(0.6,0.57,"#mu#mu events");

  can->cd();

  TPad* respad = new TPad("respad","respad",0.0,0.8,1.0,1.0);
  respad->Draw();
  respad->cd();
  respad->SetGridy();
  
  TH1F* hratio = (TH1F*) hZ->Clone("hratio");
  TH1F* hmctot = (TH1F*) hg_rw->Clone("hmctot");
  for( int ibin = 1 ; ibin <= hmctot->GetXaxis()->GetNbins() ; ibin++ ) hmctot->SetBinError(ibin,0.0);
  hratio->Divide(hmctot);

  hratio->GetXaxis()->SetTitle("");
  hratio->GetXaxis()->SetLabelSize(0.);
  hratio->GetYaxis()->SetNdivisions(5);
  hratio->GetYaxis()->SetTitle("data/bkg");
  hratio->GetYaxis()->SetTitleSize(0.15);
  hratio->GetYaxis()->SetTitleOffset(0.3);
  hratio->GetYaxis()->SetLabelSize(0.15);
  hratio->SetMinimum(0.0);
  hratio->SetMaximum(2.0);
  hratio->GetYaxis()->SetRangeUser(0.0,2.0);
  hratio->Draw("E1");

  can->Print(Form(plotsPath + "mm_NoSmear_HT_MC_ZptHTreweigh.pdf"));
	
}

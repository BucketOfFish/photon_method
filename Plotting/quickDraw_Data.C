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

void quickDraw_Data( string period = "data15-16" , string channel  = "ee" , string var = "MET_tenacious" , bool VgSubtracted = false , string smearing_mode = "NoSmear" , bool normalize = true ) {
  
  if( TString(var).Contains("pt") ) normalize = false;

  if( TString(var).Contains("HT") ) normalize = false;

  bool DF = TString(channel).EqualTo("em");

  if( DF ) normalize = false;

  
  //  cout << "smearing mode " << smearing_mode << endl;
  //  exit(0);
  
  gStyle->SetOptStat(0);

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

  // string path2l         = "/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.1/";  
  // string Zfilename      = path2l + "Data/" + period + "_merged_processed.root";
  // string tt_filename    = path2l + mcdir + "ttbar_merged_processed.root";
  // string vv_filename    = path2l + mcdir + "diboson_merged_processed.root";

  // for v1.1 ntuples
  // string Zfilename      = "../OutputNtuples/Zdata/"       + period     + "_merged_processed.root";
  // string tt_filename    = "../OutputNtuples/"             + mcdir      + "ttbar_merged_processed.root";
  // string vv_filename    = "../OutputNtuples/"             + mcdir      + "diboson_merged_processed.root";
  // string gfilename      = "../OutputNtuples/gdata/JETM4_" + gdatalabel + "_" + channel + "_" + smearing_mode + ".root";
  // string vg_filename    = "../OutputNtuples/Vg/Vgall.root";

  //string VgString = "";
  //if( VgSubtracted ) VgString = "_VgSubtracted";
  
  string Zfilename      = ntuple_path + "zdata/data15-16_merged_processed.root";
  string tt_filename    = ntuple_path + mcdir + "ttbar_merged_processed.root";
  string vv_filename    = ntuple_path + mcdir + "diboson_merged_processed.root";
  string z_filename     = ntuple_path + mcdir + "Zjets_merged_processed.root";
  //string o_filename     = ntuple_path + mcdir + "triboson_higgs_topOther_merged_processed.root";
  //string gfilename      = ntuple_path + "gdata/" + period + VgString + "_" + channel + "_" + smearing_mode + ".root";
  string gfilename      = ntuple_path + "gdata/" + period + "_merged_processed_" + channel + "_" + smearing_mode + ".root";
  //string vg_filename    = ntuple_path + "gmc/Vgamma_merged_processed.root";
  
  cout << "period               " << period        << endl;
  cout << "channel              " << channel       << endl;
  cout << "smearing mode        " << smearing_mode << endl;
  cout << "Z datafilename       " << Zfilename     << endl;
  cout << "tt filename          " << tt_filename   << endl;
  cout << "vv filename          " << vv_filename   << endl;
  cout << "Z+jets filename      " << z_filename    << endl;
  //cout << "other filename       " << o_filename    << endl;
  //cout << "Vg filename          " << vg_filename   << endl;
  cout << "g filename           " << gfilename     << endl;
  cout << "DF?                  " << DF            << endl;
  
  //-----------------------------------------------
  // add files to TChain
  //-----------------------------------------------

  TChain * gtree = new TChain("BaselineTree");
  if( !DF ) gtree->Add( gfilename.c_str() );
	
  TChain * Ztree = new TChain("BaselineTree");
  Ztree->Add( Zfilename.c_str() );

  TChain* chtt = new TChain("BaselineTree");
  chtt->Add( tt_filename.c_str() );

  TChain* chz = new TChain("BaselineTree");
  chz->Add( z_filename.c_str() );

  //TChain* cho = new TChain("BaselineTree");
  //cho->Add( o_filename.c_str() );

  TChain* chvv = new TChain("BaselineTree");
  chvv->Add( vv_filename.c_str() );

  // TChain* chvg = new TChain("BaselineTree");
  //chvg->Add( vg_filename.c_str() );

  cout << "g entries            " << gtree->GetEntries()  << endl;
  cout << "Z data entries       " << Ztree->GetEntries()  << endl;
  cout << "ttbar entries        " << chtt->GetEntries()   << endl;
  cout << "diboson entries      " << chvv->GetEntries()   << endl;
  cout << "Z+jets entries       " << chz->GetEntries()    << endl;
  //cout << "other entries        " << cho->GetEntries()    << endl;
  //cout << "Vg entries           " << chvg->GetEntries()   << endl;

  //-----------------------------------------------
  // define selections
  //-----------------------------------------------
	
  //TCut Zselection("mll>81 && mll<101 && jet_n >= 2 && MET<200 && is_OS && lep_pT[0]>25.0 && lep_pT[1]>25.0 && bjet_n==0");
  TCut Zselection("jet_n >= 2 && MET<200 && is_OS && lep_pT[0]>25.0 && lep_pT[1]>25.0 && mjj<60 || mjj>100 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && jet_pT[0]>30 && jet_pT[1]>30 && mll>80 && mll<100 && HT>800");
  TCut Zweight("totalWeight");
  TCut lumi1516("36100");
  TCut lumi17("44000");
  
  //TCut gselection("lep_pT[0]>25 && lep_pT[1]>25 && jet_n>=2  && bjet_n==0");
  TCut gselection("jet_n>=2 && mjj<60 || mjj>100 && lep_pT[0]>25 && lep_pT[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && HT>800");

  //TCut vgselection("jet_n>=2  && bjet_n==0");
  TCut vgselection("jet_n>=2");
  //TCut ZCR("MET<60.0");
  TCut CR("MET<60.0");

  TCut VR("jet_n>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lep_pT[0]>25 && lep_pT[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && mll>80 && mll<100 && mjj<60 && mjj>100");  

  // TCut ee("trigMatch_1L2LTrig && (lepFlavor[0] == 1 && lepFlavor[1] == 1)");
  // TCut mm("lepFlavor[0] == 2 && lepFlavor[1] == 2");
  // TCut em("trigMatch_1L2LTrig && ( (lepFlavor[0] == 1 && lepFlavor[1] == 2) || (lepFlavor[0] == 2 && lepFlavor[1] == 1) )");

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

  if( TString(period).EqualTo("data15-16") ) Zweight *= lumi1516;
  if( TString(period).EqualTo("data17")    ) Zweight *= lumi17;

  TCut weight_g    = "totalWeight";
  //TCut weight_g_rw = "totalWeight*ptreweight3";
  TCut weight_g_rw = "totalWeight*ptreweight3*ptreweight5"; //2-STEP

  cout << "VR selection         " << VR.GetTitle()  << endl;
  cout << "Z selection          " << Zselection.GetTitle()  << endl;  
  cout << "Z weight             " << Zweight.GetTitle()     << endl;
  cout << "g selection          " << gselection.GetTitle()  << endl;
  cout << "g weight             " << weight_g.GetTitle()    << endl;
  cout << "g weight (reweight)  " << weight_g_rw.GetTitle() << endl;
  
  //-----------------------------------------------
  // define and draw histograms
  //-----------------------------------------------

  const unsigned int nptbins = 16;
  double ptbins[nptbins+1] = {40, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500, 600, 700, 850, 1000};

  int nmetbins =  20;
  float metmax = 200;
  
  TH1F* hZ    ;//= new TH1F();    
  TH1F* hg    ;//= new TH1F();    
  TH1F* htt   ;//= new TH1F();   
  TH1F* hvv   ;//= new TH1F();
  TH1F* hz    ;//= new TH1F();
  //TH1F* ho    ;//= new TH1F();
  //TH1F* hvg   ;//= new TH1F();
  //TH1F* hvg_rw;//= new TH1F();   
  TH1F* hg_rw ;//= new TH1F(); 

  string xtitle = var;

  /*
  if( TString(var).EqualTo("met") ){
    xtitle = "E_{T}^{miss} [GeV]";

    hZ    = new TH1F("hZ",   "" , nmetbins,0,metmax);
    hg    = new TH1F("hg",   "" , nmetbins,0,metmax);
    htt   = new TH1F("htt",  "" , nmetbins,0,metmax);
    hvv   = new TH1F("hvv",  "" , nmetbins,0,metmax);
    hg_rw = new TH1F("hg_rw","" , nmetbins,0,metmax);

    Ztree->Draw("met_Et>>hZ"    , Zselection              , "goff");
    chtt-> Draw("met_Et>>htt"   , Zselection*Zweight      , "goff");
    chvv-> Draw("met_Et>>hvv"   , Zselection*Zweight      , "goff");
    if( !DF ){
      gtree->Draw("MET>>hg"       , gselection*weight_g     , "goff");
      gtree->Draw("MET>>hg_rw"    , gselection*weight_g_rw  , "goff");
    }
  }

  else if( TString(var).EqualTo("pt") ){
    xtitle = "p_{T} [GeV]";

    hZ    = new TH1F("hZ",   "" , nptbins , ptbins );
    hg    = new TH1F("hg",   "" , nptbins , ptbins );
    htt   = new TH1F("htt",  "" , nptbins , ptbins );
    hvv   = new TH1F("hvv",  "" , nptbins , ptbins );
    hg_rw = new TH1F("hg_rw","" , nptbins , ptbins );

    Ztree->Draw("Ptll>>hZ"     , Zselection              , "goff");
    chtt-> Draw("Ptll>>htt"    , Zselection*Zweight      , "goff");
    chvv-> Draw("Ptll>>hvv"    , Zselection*Zweight      , "goff");
    if( !DF ){
      gtree->Draw("Z_pt>>hg"     , gselection*weight_g     , "goff");
      gtree->Draw("Z_pt>>hg_rw"  , gselection*weight_g_rw  , "goff");
    }
  }
  */

  int   nbins =  20;
  float xmin  =   0;
  float xmax  = 100;
  
  if( TString(var).EqualTo("MET") ){
    xtitle = "E_{T}^{miss} [GeV]";
    nbins = nmetbins;
    xmin  = 0.0;
    xmax  = metmax;
  }
  
  else if( TString(var).EqualTo("METl") || TString(var).EqualTo("METt") ){
    if( TString(var).EqualTo("METl") ) xtitle = "E_{T,||}^{miss} [GeV]";
    if( TString(var).EqualTo("METt") ) xtitle = "E_{T,#perp}^{miss} [GeV]";
    nbins =   20;
    xmin  = -200;
    xmax  =  200;
  }

  else if( TString(var).EqualTo("MET_loose") || TString(var).EqualTo("MET_tight") || TString(var).EqualTo("MET_tighter") || TString(var).EqualTo("MET_tenacious") ){
    if( TString(var).EqualTo("MET_loose") ) xtitle = "E_{T,loose}^{miss} [GeV]";
    if( TString(var).EqualTo("MET_tight") ) xtitle = "E_{T,tight}^{miss} [GeV]";
    if( TString(var).EqualTo("MET_tighter") ) xtitle = "E_{T,tighter}^{miss} [GeV]";
    if( TString(var).EqualTo("MET_tenacious") ) xtitle = "E_{T,tenacious}^{miss} [GeV]";
    nbins = 20;
    xmin  = 0;
    xmax  = 200;
  }

  else if( TString(var).EqualTo("Z_pt") ){
    xtitle = "p_{T} [GeV]";
  }

  else if( TString(var).EqualTo("HT") ){
    xtitle = "H_{T}";
    //nbins =   20;
    //xmin  =    0;
    //xmax  = 1000;
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

  // initialize histograms
  
  if( TString(var).EqualTo("Z_pt") ){
    hZ     = new TH1F("hZ"     , "" , nptbins , ptbins );
    hg     = new TH1F("hg"     , "" , nptbins , ptbins );
    htt    = new TH1F("htt"    , "" , nptbins , ptbins );
    hvv    = new TH1F("hvv"    , "" , nptbins , ptbins );
    hz     = new TH1F("hz"     , "" , nptbins , ptbins );
    //ho     = new TH1F("ho"     , "" , nptbins , ptbins );
    //hvg    = new TH1F("hvg"    , "" , nptbins , ptbins );
    //hvg_rw = new TH1F("hvg_rw" , "" , nptbins , ptbins );
    hg_rw  = new TH1F("hg_rw"  , "" , nptbins , ptbins );
  }

  else if( TString(var).EqualTo("HT") ){
    hZ     = new TH1F("hZ"     , "" , nptbins , ptbins );
    hg     = new TH1F("hg"     , "" , nptbins , ptbins );
    htt    = new TH1F("htt"    , "" , nptbins , ptbins );
    hvv    = new TH1F("hvv"    , "" , nptbins , ptbins );
    hz     = new TH1F("hz"     , "" , nptbins , ptbins );
    //ho     = new TH1F("ho"     , "" , nptbins , ptbins );
    //hvg    = new TH1F("hvg"    , "" , nptbins , ptbins );
    //hvg_rw = new TH1F("hvg_rw" , "" , nptbins , ptbins );
    hg_rw  = new TH1F("hg_rw"  , "" , nptbins , ptbins );
  }

  else{
    hZ     = new TH1F("hZ"     , "" , nbins , xmin , xmax );
    hg     = new TH1F("hg"     , "" , nbins , xmin , xmax );
    htt    = new TH1F("htt"    , "" , nbins , xmin , xmax );
    hvv    = new TH1F("hvv"    , "" , nbins , xmin , xmax );
    hz     = new TH1F("hz"     , "" , nbins , xmin , xmax );
    //ho     = new TH1F("ho"     , "" , nbins , xmin , xmax );
    //hvg    = new TH1F("hvg"    , "" , nbins , xmin , xmax );
    //hvg_rw = new TH1F("hvg_rw" , "" , nbins , xmin , xmax );
    hg_rw  = new TH1F("hg_rw"  , "" , nbins , xmin , xmax );
  }


  
  
  Ztree->Draw(Form("%s>>hZ",var.c_str())       , Zselection              , "goff");
  chtt-> Draw(Form("%s>>htt",var.c_str())      , Zselection*Zweight      , "goff");
  chz->  Draw(Form("%s>>hz",var.c_str())       , Zselection*Zweight      , "goff");
  //cho->  Draw(Form("%s>>ho",var.c_str())       , Zselection*Zweight      , "goff");
  chvv-> Draw(Form("%s>>hvv",var.c_str())      , Zselection*Zweight      , "goff");

  if( !DF ){
    gtree->Draw(Form("%s>>hg",var.c_str())     , gselection*weight_g     , "goff");
    gtree->Draw(Form("%s>>hg_rw",var.c_str())  , gselection*weight_g_rw  , "goff");
  }
  
  /*
  //-------------------------------
  //----------VR Studies-----------
  //-------------------------------
  Ztree->Draw(Form("%s>>hZ",var.c_str())       , VR              , "goff");
  chtt-> Draw(Form("%s>>htt",var.c_str())      , VR*Zweight      , "goff");
  chz->  Draw(Form("%s>>hz",var.c_str())       , VR*Zweight      , "goff");
  //cho->  Draw(Form("%s>>ho",var.c_str())       , Zselection*Zweight      , "goff");
  chvv-> Draw(Form("%s>>hvv",var.c_str())      , VR*Zweight      , "goff");

  if( !DF ){
    gtree->Draw(Form("%s>>hg",var.c_str())     , VR*weight_g     , "goff");
    gtree->Draw(Form("%s>>hg_rw",var.c_str())  , VR*weight_g_rw  , "goff");
  }
*/

  //chvg-> Draw("MET_raw>>hvg"      , vgselection*weight_g*lumi1516     , "goff");
  //chvg-> Draw("MET_raw>>hvg_rw"   , vgselection*weight_g_rw*lumi1516  , "goff");
  
  cout << "Z data integral      " << hZ->Integral()  << endl;
  cout << "tt MC integral       " << htt->Integral()  << endl;
  cout << "Z+jets MC integral   " << hz->Integral()  << endl;
  //cout << "other MC integral    " << ho->Integral()  << endl;
  cout << "VV MC integral       " << hvv->Integral()  << endl;
  //cout << "Vg MC integral       " << hvg->Integral()  << endl;
  //cout << "Vg MC integral (rw)  " << hvg_rw->Integral()  << endl;
  cout << "g data raw integral  " << hg->Integral()  << endl;
  cout << "g data rw integral   " << hg_rw->Integral()  << endl;

  //tree1->Draw("Z_pt>>h1(20,0,500)",mycut*weight_z);
  //tree2->Draw("gamma_pt>>h2(20,0,500)",mycut*weight_g);

  //-----------------------------------------------
  // normalize Z to MET<60 GeV region
  //-----------------------------------------------

  float SF   = 1.0;
  float SFrw = 1.0;
  
  if( normalize ){
    
    cout << "normalize to CR    " << CR.GetTitle()         << endl;

    TH1F* hZnorm    = new TH1F("hZnorm",   "",1,0,1);
    TH1F* hgnorm    = new TH1F("hgnorm",   "",1,0,1);
    TH1F* hgrwnorm  = new TH1F("hgrwnorm", "",1,0,1);
    TH1F* httnorm   = new TH1F("httnorm",  "",1,0,1);
    TH1F* hvvnorm   = new TH1F("hvvnorm",  "",1,0,1);

    Ztree->Draw("0.5>>hZnorm"     , Zselection+CR                , "goff");
    chtt-> Draw("0.5>>httnorm"    , (Zselection+CR)*Zweight      , "goff");
    chvv-> Draw("0.5>>hvvnorm"    , (Zselection+CR)*Zweight      , "goff");

    gtree->Draw("0.5>>hgnorm"     , (gselection+CR)*weight_g     , "goff");
    gtree->Draw("0.5>>hgrwnorm"   , (gselection+CR)*weight_g_rw  , "goff");

    SF   = ( hZnorm->Integral() - httnorm->Integral() - hvvnorm->Integral() ) / hgnorm->Integral();
    SFrw = ( hZnorm->Integral() - httnorm->Integral() - hvvnorm->Integral() ) / hgrwnorm->Integral();

    cout << "Scale reweighted Z by    " << SFrw << endl;
    cout << "Scale raw Z by           " << SF   << endl;
    
    hg->Scale(SF);
    hg_rw->Scale(SFrw);

  }

  if( TString(var).EqualTo("Z_pt") ) hg->Scale( hg_rw->Integral() / hg->Integral() );
      
  //-----------------------------------------------
  // make pretty plots
  //-----------------------------------------------

  cout << "MET100-150" << endl;
  cout << "2L data                " << hZ->Integral(11,15)     << endl;
  cout << "g data (reweighted)    " << hg_rw->Integral(11,15)  << endl;
  cout << "g data (raw)           " << hg->Integral(11,15)     << endl;
  //cout << "Vg MC                  " << hvg->Integral(11,15)    << endl;
  cout << "VV MC                  " << hvv->Integral(11,15)    << endl;
  cout << "tt MC                  " << htt->Integral(11,15)    << endl;
  //cout << "other MC               " << ho->Integral(11,15)     << endl;
  cout << "Z+jets MC              " << hz->Integral(11,15)     << endl;

  cout << "MET150-200" << endl;
  cout << "2L data                " << hZ->Integral(16,21)     << endl;
  cout << "g data (reweighted)    " << hg_rw->Integral(16,21)  << endl;
  cout << "g data (raw)           " << hg->Integral(16,21)     << endl;
  //cout << "Vg MC                  " << hvg->Integral(16,21)    << endl;
  cout << "VV MC                  " << hvv->Integral(16,21)    << endl;
  cout << "tt MC                  " << htt->Integral(16,21)    << endl;
  //cout << "other MC               " << ho->Integral(16,21)     << endl;
  cout << "Z+jets MC              " << hz->Integral(16,21)     << endl;

  
  // make canvas and draw 2L data vs. MC plot
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

  htt->SetLineColor(1);
  htt->SetFillColor(kRed-2);

  //ho->SetLineColor(1);
  //ho->SetFillColor(kMagenta-2);

  hvv->SetLineColor(1);
  hvv->SetFillColor(kGreen-2);

  hg_rw->SetLineColor(1);
  hg_rw->SetFillColor(kOrange-2);

  hg->Add(htt);
  hg->Add(hvv);
  //hg->Add(ho);

  hg->SetLineColor(4);
  hg->SetLineWidth(1);
  hg->SetLineStyle(2);

  hz->Add(htt);
  hz->Add(hvv);
  //hz->Add(ho);

  hz->SetLineColor(2);
  hz->SetLineWidth(1);
  hz->SetLineStyle(7);

  THStack *mcstack = new THStack("mcstack","mcstack");
  //mcstack->Add(ho);
  mcstack->Add(htt);
  mcstack->Add(hvv);
  if( !DF ) mcstack->Add(hg_rw);
  mcstack->Draw("samehist");
  hZ->Draw("sameE1");
  if( !DF ) hg->Draw("samehist");
  hZ->Draw("axissame");

  TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
  leg->AddEntry(hZ,"data","lp");
  if( !DF ){
    leg->AddEntry(hg_rw,"Z+jets (from #gamma+jets, reweighted)","f");
    leg->AddEntry(hg,"Z+jets (from #gamma+jets, raw)","f");
    leg->AddEntry(hz,"Z+jets (from MC)","f");
  }
  leg->AddEntry(hvv,"VV","f");
  leg->AddEntry(htt,"t#bar{t}+tW","f");
  //leg->AddEntry(ho,"Other MC","f");
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

  can->cd();

  TPad* respad = new TPad("respad","respad",0.0,0.8,1.0,1.0);
  respad->Draw();
  respad->cd();
  respad->SetGridy();
  
  TH1F* hratio = (TH1F*) hZ->Clone("hratio");
  TH1F* hmctot = (TH1F*) hg_rw->Clone("hmctot");
  hmctot->Add(htt);
  hmctot->Add(hvv);
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
  
  //can->Print(Form("plots/quickData_Data_%s_%s_%s_%s_VgSubtraction.pdf",period.c_str(),channel.c_str(),var.c_str(),smearing_mode.c_str()));

  can->Print(Form("%sVR_SR_studies/quickData_Data_%s_%s_%s_%s_VR_ht800cut.pdf",plots_path.c_str(),period.c_str(),channel.c_str(),var.c_str(),smearing_mode.c_str()));


  // TFile* outfile = TFile::Open( Form("plots/quickData_Data_%s_%s_%s_%s.root",period.c_str(),channel.c_str(),var.c_str(),smearing_mode.c_str()) , "RECREATE" );
  // outfile->cd();
  // hratio->Write();
  // outfile->Close();

  
  
  // TCanvas* can = new TCanvas("can","can",600,600);
  // can->cd();

  // gPad->SetLogy();
	
  // hZ->SetLineColor(1);
  // hZ->SetLineWidth(2);

  // hg->SetLineColor(4);
  // hg->SetLineWidth(1);
  // hg->SetLineStyle(9);

  // hg_rw->SetLineColor(2);
  // hg_rw->SetLineWidth(2);
  // hg_rw->SetLineStyle(2);

  // hZ->Draw("hist");
  // hg->Draw("samehist");
  // hg_rw->Draw("samehist");
	 
  // // h1->SetLineColor(4);
  // // h1->SetStats(0);
  // // h2->SetLineColor(2);
  // // h2->SetStats(0);
  // // h1->Draw();
  // // h2->Draw("same");
  // // c1->SetLogy();

  // can->Print( Form( "quickDraw_Data_normalized_%s.pdf",var.c_str() ) );
	
}

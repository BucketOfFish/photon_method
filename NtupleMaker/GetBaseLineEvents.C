//-----------------------------------------------------------------------------------------------
// this script takes the non-photon ntuples from QuickAna and generates smaller ntuples with information required by photon method
// the parameters of the function GetBaseLineEvents(string sampleID, string outputName, string pathToNtuples, bool isData, stringtreename = "outputTree" ) are:
// 	sampleID: DSID of the MC sample, or simply the input file name
//      outputName: name of output directory 
// 	pathToNtuples: the path to the input ntuple
// 	isData: put "1" if you are running a data sample, "0" if MC 
//      stringtreename: check tree, usually Zjets_NoSys
// example of code running command: root -l -b -q 'GetBaseLineEvents.C+("Zjets_merged_processed","zmc","/afs/cern.ch/user/b/benhoob/SusySkim2LJets/v1.2/SUSY2/SUSY2_Bkgs_mc16a/",0,"Zjets_NoSys")'
//-----------------------------------------------------------------------------------------------


#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <string>
#include <iomanip> 

#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TObject.h"

#include "../BasicSetting.C"
#include "../CommonFunctions/CommonFunctions.C"
#include "InputVariables.C"
#include "../PhotonSmearing/GetDijetVariables.C"

using namespace std;

void RebinHistogram(TH1D* hist) {
    float negative_yield = 0.;
    float positive_yield = 0.;
    for (int bin=1;bin<=hist->GetNbinsX();bin++) {
        if (hist->GetBinContent(bin)>=0) positive_yield += hist->GetBinContent(bin);
        else negative_yield += hist->GetBinContent(bin);
    }
    while (abs(negative_yield/positive_yield)>0.05) {
        hist->Rebin(2);
        negative_yield = 0.;
        positive_yield = 0.;
        for (int bin=1;bin<=hist->GetNbinsX();bin++) {
            if (hist->GetBinContent(bin)>=0) positive_yield += hist->GetBinContent(bin);
            else negative_yield += hist->GetBinContent(bin);
        }
    }
}

void GetBaseLineEvents(string sampleID, string outputName, string pathToNtuples, string isData, string treename = "outputTree" ) {

    //---------------------------------------------
    // open file, get Tree and EventCountHist
    //---------------------------------------------

    TH1::SetDefaultSumw2();

    TH1D* hist_EventCount = new TH1D("hist_EventCount","",3,0,3);
    float N_passMET100 = 0.;
    string  filename       = Form("%s%s.root",pathToNtuples.c_str(),sampleID.c_str()); 
    TFile*  inputFile      = TFile::Open(filename.c_str());
    Float_t _nGenEvents = 1.;
    if (isData == "MC") {
        cout << "Setting _nGenEvents = 1 for now NEED TO FIX" << endl;
        _nGenEvents    = 1.0;
        hist_EventCount->SetBinContent(1,1.0);
    }
    TTree*  inputTree              = (TTree*)inputFile->Get( treename.c_str() );

    std::cout << inputTree << std::endl;
    cout << endl;
    cout << "Opening file           : " << filename        << endl;
    cout << "Events in ntuple       : " << inputTree->GetEntries() << endl;
    if (isData == "MC") {
        cout << "Total generated events : " << _nGenEvents     << endl;
    }

    //-----------------------------
    // access existing branches
    //-----------------------------
    inputTree->SetBranchStatus("*", 0);
    bool trigMatch_1L2LTrig; SetInputBranch(inputTree, "trigMatch_1L2LTrig", &trigMatch_1L2LTrig);
    ULong64_t EventNumber; SetInputBranch(inputTree, "EventNumber", &EventNumber);
    Int_t RunNumber; SetInputBranch(inputTree, "RunNumber", &RunNumber);
    Float_t Mu; SetInputBranch(inputTree, "mu", &Mu);
    Int_t nVtx; SetInputBranch(inputTree, "nVtx", &nVtx);
    Float_t mll; SetInputBranch(inputTree, "mll", &mll);
    float MET; SetInputBranch(inputTree, "met_Et", &MET);
    float MET_phi; SetInputBranch(inputTree, "met_Phi", &MET_phi);
    float MET_softTerm; SetInputBranch(inputTree, "TST_Et", &MET_softTerm);
    float MET_softPhi; SetInputBranch(inputTree, "TST_Phi", &MET_softPhi);
    Bool_t trigMatch_2LTrigOR; SetInputBranch(inputTree, "trigMatch_2LTrigOR", &trigMatch_2LTrigOR);
    float MET_loose; SetInputBranch(inputTree, "met_Et_loose", &MET_loose);
    float MET_tight; SetInputBranch(inputTree, "met_Et_tight", &MET_tight);
    float MET_tighter; SetInputBranch(inputTree, "met_Et_tighter", &MET_tighter);
    float MET_tenacious; SetInputBranch(inputTree, "met_Et_tenacious", &MET_tenacious);
    Bool_t is2Lep2Jet; SetInputBranch(inputTree, "is2Lep2Jet", &is2Lep2Jet);
    Bool_t is2L2JInt; SetInputBranch(inputTree, "is2L2JInt", &is2L2JInt);
    int nBJet20_MV2c10_FixedCutBEff_77; SetInputBranch(inputTree, "nBJet20_MV2c10_FixedCutBEff_77", &nBJet20_MV2c10_FixedCutBEff_77);
    float mjj; SetInputBranch(inputTree, "mjj", &mjj);
    Double_t mll_RJ; SetInputBranch(inputTree, "mll_RJ", &mll_RJ);
    Double_t R_minH2P_minH3P; SetInputBranch(inputTree, "R_minH2P_minH3P", &R_minH2P_minH3P);
    Double_t RPT_HT5PP; SetInputBranch(inputTree, "RPT_HT5PP", &RPT_HT5PP);
    Double_t dphiVP; SetInputBranch(inputTree, "dphiVP", &dphiVP);
    Double_t H2PP; SetInputBranch(inputTree, "H2PP", &H2PP);
    Double_t H5PP; SetInputBranch(inputTree, "H5PP", &H5PP);
    int nJet20; SetInputBranch(inputTree, "nJet20", &nJet20);
    Double_t minDphi; SetInputBranch(inputTree, "minDphi", &minDphi);
    Double_t MZ; SetInputBranch(inputTree, "MZ", &MZ);
    Double_t NjS; SetInputBranch(inputTree, "NjS", &NjS);
    Double_t NjISR; SetInputBranch(inputTree, "NjISR", &NjISR);
    Double_t dphiISRI; SetInputBranch(inputTree, "dphiISRI", &dphiISRI);
    Double_t RISR; SetInputBranch(inputTree, "RISR", &RISR);
    Double_t PTISR; SetInputBranch(inputTree, "PTISR", &PTISR);
    Double_t PTI; SetInputBranch(inputTree, "PTI", &PTI);
    Double_t PTCM; SetInputBranch(inputTree, "PTCM", &PTCM);
    Double_t MJ; SetInputBranch(inputTree, "MJ", &MJ);
    Bool_t is3Lep3Jet; SetInputBranch(inputTree, "is3Lep3Jet", &is3Lep3Jet);
    Bool_t is4Lep3Jet; SetInputBranch(inputTree, "is4Lep3Jet", &is4Lep3Jet);
    Double_t lept1sign_VR; SetInputBranch(inputTree, "lept1sign_VR", &lept1sign_VR);
    Double_t lept2sign_VR; SetInputBranch(inputTree, "lept2sign_VR", &lept2sign_VR);
    Double_t lept1Pt_VR; SetInputBranch(inputTree, "lept1Pt_VR", &lept1Pt_VR);
    Double_t lept2Pt_VR; SetInputBranch(inputTree, "lept2Pt_VR", &lept2Pt_VR);
    Double_t MZ_VR; SetInputBranch(inputTree, "MZ_VR", &MZ_VR);
    Double_t MJ_VR; SetInputBranch(inputTree, "MJ_VR", &MJ_VR);
    Double_t RISR_VR; SetInputBranch(inputTree, "RISR_VR", &RISR_VR);
    Double_t PTISR_VR; SetInputBranch(inputTree, "PTISR_VR", &PTISR_VR);
    Double_t PTI_VR; SetInputBranch(inputTree, "PTI_VR", &PTI_VR);
    Double_t PTCM_VR; SetInputBranch(inputTree, "PTCM_VR", &PTCM_VR);
    Double_t dphiISRI_VR; SetInputBranch(inputTree, "dphiISRI_VR", &dphiISRI_VR);
    float DPhi_METJetLeading; SetInputBranch(inputTree, "DPhiJ1Met", &DPhi_METJetLeading);
    float DPhi_METJetSecond; SetInputBranch(inputTree, "DPhiJ2Met", &DPhi_METJetSecond);
    float HT; SetInputBranch(inputTree, "Ht30", &HT);
    Float_t Z_pt; SetInputBranch(inputTree, "Ptll", &Z_pt);
    Int_t jet_n; SetInputBranch(inputTree, "nJet30", &jet_n);
    Int_t bjet_n; SetInputBranch(inputTree, "nBJet30_MV2c10_FixedCutBEff_77", &bjet_n);
    Int_t nLep_signal; SetInputBranch(inputTree, "nLep_signal", &nLep_signal);
    Int_t nLep_base; SetInputBranch(inputTree, "nLep_base", &nLep_base);
    //std::vector<int>* lepFlavor = new std::vector<int>(10); SetInputBranch(inputTree, "lepFlavor", &lepFlavor);
    //std::vector<int>* lepCharge = new std::vector<int>(10); SetInputBranch(inputTree, "lepCharge", &lepCharge);
    //std::vector<float>* lep_pT = new std::vector<float>(10); SetInputBranch(inputTree, "lepPt", &lep_pT);
    //std::vector<float>* lep_eta = new std::vector<float>(10); SetInputBranch(inputTree, "lepEta", &lep_eta);
    //std::vector<float>* lep_phi = new std::vector<float>(10); SetInputBranch(inputTree, "lepPhi", &lep_phi);
    //std::vector<float>* jet_pT = new std::vector<float>(10); SetInputBranch(inputTree, "jetPt", &jet_pT);
    //std::vector<float>* jet_eta = new std::vector<float>(10); SetInputBranch(inputTree, "jetEta", &jet_eta);
    //std::vector<float>* jet_phi = new std::vector<float>(10); SetInputBranch(inputTree, "jetPhi", &jet_phi);
    //std::vector<float>* jet_m = new std::vector<float>(10); SetInputBranch(inputTree, "jetM", &jet_m);
    SetInputBranch(inputTree, "lepFlavor", &lepFlavor);
    SetInputBranch(inputTree, "lepCharge", &lepCharge);
    SetInputBranch(inputTree, "lepPt", &lep_pT);
    SetInputBranch(inputTree, "lepEta", &lep_eta);
    SetInputBranch(inputTree, "lepPhi", &lep_phi);
    SetInputBranch(inputTree, "jetPt", &jet_pT);
    SetInputBranch(inputTree, "jetEta", &jet_eta);
    SetInputBranch(inputTree, "jetPhi", &jet_phi);
    SetInputBranch(inputTree, "jetM", &jet_m);

    Double_t genWeight;
    Double_t eventWeight;
    Double_t leptonWeight;
    Double_t jvtWeight;
    Double_t bTagWeight;
    Double_t pileupWeight;
    Double_t FFWeight;
    if (isData == "MC") {
        SetInputBranch(inputTree, "genWeight", &genWeight);
        SetInputBranch(inputTree, "eventWeight", &eventWeight);
        SetInputBranch(inputTree, "leptonWeight", &leptonWeight);
        SetInputBranch(inputTree, "jvtWeight", &jvtWeight);
        SetInputBranch(inputTree, "bTagWeight", &bTagWeight);
        SetInputBranch(inputTree, "pileupWeight", &pileupWeight);
        SetInputBranch(inputTree, "FFWeight", &FFWeight);
    }

    //-----------------------------
    // add new branches
    //-----------------------------

    string outfilename = ntuple_path + "/" + outputName + "/" + sampleID.c_str() + ".root";
    cout << "Writing to : " << outfilename << endl;

    TFile   outputFile( outfilename.c_str() , "recreate" );

    TTree* BaselineTree;
    BaselineTree = new TTree("BaselineTree","baseline tree");

    TH1D* hist_cutflow_raw = new TH1D("hist_cutflow_raw","",8,0,8);
    TH1D* hist_cutflow_weight = new TH1D("hist_cutflow_weight","",8,0,8);

    TH1D* hist_dPt_Pt[bin_size];
    TH1D* hist_dPhi_Pt[bin_size];
    TH1D* hist_METl_Pt[bin_size];
    TH1D* hist_METt_Pt[bin_size];
    TH1D* hist_JetMETl_Pt[bin_size];
    TH1D* hist_2LPt_Pt[bin_size];
    TH1D* hist_Mll_dPt[dpt_bin_size];
    TH1D* hist_fsee_METl_Pt[bin_size];
    TH1D* hist_fsmm_METl_Pt[bin_size];

    Float_t MT2_max= 0;
    TBranch *b_MT2_max = BaselineTree->Branch("MT2_max",&MT2_max,"MT2_max/F");
    Float_t boost_phi= 0;
    TBranch *b_boost_phi = BaselineTree->Branch("boost_phi",&boost_phi,"boost_phi/F");
    Float_t boost_eta= 0;
    TBranch *b_boost_eta = BaselineTree->Branch("boost_eta",&boost_eta,"boost_eta/F");
    Float_t boost_pt= 0;
    TBranch *b_boost_pt = BaselineTree->Branch("boost_pt",&boost_pt,"boost_pt/F");

    BaselineTree->Branch("DPhi_METNonWJet",&DPhi_METNonWJet,"DPhi_METNonWJet/F");
    BaselineTree->Branch("NonWJet_pT",&NonWJet_pT,"NonWJet_pT/F");
    BaselineTree->Branch("DPhi_METNonWminJet",&DPhi_METNonWminJet,"DPhi_METNonWminJet/F");
    BaselineTree->Branch("NonWminJet_pT",&NonWminJet_pT,"NonWminJet_pT/F");
    BaselineTree->Branch("lepFlavor","std::vector<int>",&lepFlavor);
    BaselineTree->Branch("lepCharge","std::vector<int>",&lepCharge);
    BaselineTree->Branch("lep_pT","std::vector<float>",&lep_pT);
    BaselineTree->Branch("lep_phi","std::vector<float>",&lep_phi);
    BaselineTree->Branch("lep_eta","std::vector<float>",&lep_eta);
    BaselineTree->Branch("jet_pT","std::vector<float>",&jet_pT);
    BaselineTree->Branch("jet_phi","std::vector<float>",&jet_phi);
    BaselineTree->Branch("jet_eta","std::vector<float>",&jet_eta);
    BaselineTree->Branch("jet_m","std::vector<float>",&jet_m);
    BaselineTree->Branch("Mu",&Mu,"Mu/F");
    BaselineTree->Branch("nVtx",&nVtx,"nVtx/I");
    BaselineTree->Branch("mll",&mll,"mll/F");
    BaselineTree->Branch("MET",&MET,"MET/F");
    // 2019 RJR analysis variables -------------------------
    BaselineTree->Branch("MET_loose",&MET_loose,"MET_loose/F");
    BaselineTree->Branch("MET_tight",&MET_tight,"MET_tight/F");
    BaselineTree->Branch("MET_tighter",&MET_tighter,"MET_tighter/F");
    BaselineTree->Branch("MET_tenacious",&MET_tenacious,"MET_tenacious/F");
    BaselineTree->Branch("trigMatch_2LTrigOR",&trigMatch_2LTrigOR,"trigMatch_2LTrigOR/I");
    BaselineTree->Branch("is2Lep2Jet",&is2Lep2Jet,"is2Lep2Jet/I");
    BaselineTree->Branch("is2L2JInt",&is2L2JInt,"is2L2JInt/I");
    BaselineTree->Branch("nBJet20_MV2c10_FixedCutBEff_77",&nBJet20_MV2c10_FixedCutBEff_77,"nBJet20_MV2c10_FixedCutBEff_77/I");
    BaselineTree->Branch("mjj",&mjj,"mjj/F");
    BaselineTree->Branch("mll_RJ",&mll_RJ,"mll_RJ/F");
    BaselineTree->Branch("R_minH2P_minH3P",&R_minH2P_minH3P,"R_minH2P_minH3P/F");
    BaselineTree->Branch("RPT_HT5PP",&RPT_HT5PP,"RPT_HT5PP/F");
    BaselineTree->Branch("dphiVP",&dphiVP,"dphiVP/F");
    BaselineTree->Branch("H2PP",&H2PP,"H2PP/F");
    BaselineTree->Branch("H5PP",&H5PP,"H5PP/F");
    BaselineTree->Branch("nJet20",&nJet20,"nJet20/I");
    BaselineTree->Branch("minDphi",&minDphi,"minDphi/F");
    BaselineTree->Branch("MZ",&MZ,"MZ/F");
    BaselineTree->Branch("NjS",&NjS,"NjS/I");
    BaselineTree->Branch("NjISR",&NjISR,"NjISR/I");
    BaselineTree->Branch("dphiISRI",&dphiISRI,"dphiISRI/F");
    BaselineTree->Branch("RISR",&RISR,"RISR/F");
    BaselineTree->Branch("PTISR",&PTISR,"PTISR/F");
    BaselineTree->Branch("PTI",&PTI,"PTI/F");
    BaselineTree->Branch("PTCM",&PTCM,"PTCM/F");
    BaselineTree->Branch("MJ",&MJ,"MJ/F");
    BaselineTree->Branch("is3Lep3Jet",&is3Lep3Jet,"is3Lep3Jet/I");
    BaselineTree->Branch("is4Lep3Jet",&is4Lep3Jet,"is4Lep3Jet/I");
    BaselineTree->Branch("lept1sign_VR",&lept1sign_VR,"lept1sign_VR/I");
    BaselineTree->Branch("lept2sign_VR",&lept2sign_VR,"lept2sign_VR/I");
    BaselineTree->Branch("lept1Pt_VR",&lept1Pt_VR,"lept1Pt_VR/F");
    BaselineTree->Branch("lept2Pt_VR",&lept2Pt_VR,"lept2Pt_VR/F");
    BaselineTree->Branch("MZ_VR",&MZ_VR,"MZ_VR/F");
    BaselineTree->Branch("MJ_VR",&MJ_VR,"MJ_VR/F");
    BaselineTree->Branch("RISR_VR",&RISR_VR,"RISR_VR/F");
    BaselineTree->Branch("PTISR_VR",&PTISR_VR,"PTISR_VR/F");
    BaselineTree->Branch("PTI_VR",&PTI_VR,"PTI_VR/F");
    BaselineTree->Branch("PTCM_VR",&PTCM_VR,"PTCM_VR/F");
    BaselineTree->Branch("dphiISRI_VR",&dphiISRI_VR,"dphiISRI_VR/F");
    BaselineTree->Branch("METl",&METl,"METl/F");
    BaselineTree->Branch("METt",&METt,"METt/F");
    BaselineTree->Branch("MET_phi",&MET_phi,"MET_phi/F");
    BaselineTree->Branch("MET_softTerm",&MET_softTerm,"MET_softTerm/F");
    BaselineTree->Branch("MET_softPhi",&MET_softPhi,"MET_softPhi/F");
    BaselineTree->Branch("DPhi_2Lep",&DPhi_2Lep,"DPhi_2Lep/F");
    BaselineTree->Branch("DR_2Lep",&DR_2Lep,"DR_2Lep/F");
    BaselineTree->Branch("DR_Wmin2Jet",&DR_Wmin2Jet,"DR_Wmin2Jet/F");
    BaselineTree->Branch("DR_J0J1",&DR_J0J1,"DR_J0J1/F");
    BaselineTree->Branch("DPhi_METJetLeading",&DPhi_METJetLeading,"DPhi_METJetLeading/F");
    BaselineTree->Branch("DPhi_METJetSecond",&DPhi_METJetSecond,"DPhi_METJetSecond/F");
    BaselineTree->Branch("DPhi_METPhoton",&DPhi_METPhoton,"DPhi_METPhoton/F");
    BaselineTree->Branch("DPhi_METLepLeading",&DPhi_METLepLeading,"DPhi_METLepLeading/F");
    BaselineTree->Branch("DPhi_METLepSecond",&DPhi_METLepSecond,"DPhi_METLepSecond/F");
    BaselineTree->Branch("DPhi_METLepMin",&DPhi_METLepMin,"DPhi_METLepMin/F");
    BaselineTree->Branch("MinDPhi_PhotonJet",&MinDPhi_PhotonJet,"MinDPhi_PhotonJet/F");
    BaselineTree->Branch("HT",&HT,"HT/F");
    BaselineTree->Branch("Z_pt",&Z_pt,"Z_pt/F");
    BaselineTree->Branch("Z_eta",&Z_eta,"Z_eta/F");
    BaselineTree->Branch("channel",&channel,"channel/I");
    BaselineTree->Branch("is_OS",&is_OS,"is_OS/I");
    BaselineTree->Branch("Z_eta",&Z_eta,"Z_eta/F");
    BaselineTree->Branch("Z_phi",&Z_phi,"Z_phi/F");
    BaselineTree->Branch("bjet_n",&bjet_n,"bjet_n/I");
    BaselineTree->Branch("jet_n",&jet_n,"jet_n/I");
    BaselineTree->Branch("mj0j1",&mj0j1,"mj0j1/F");
    BaselineTree->Branch("W01_pt",&W01_pt,"W01_pt/F");
    BaselineTree->Branch("DPhi_METW01",&DPhi_METW01,"DPhi_METW01/F");
    BaselineTree->Branch("DPhi_W01Z",&DPhi_W01Z,"DPhi_W01Z/F");
    BaselineTree->Branch("mWmin",&mWmin,"mWmin/F");
    BaselineTree->Branch("Wmin_pt",&Wmin_pt,"Wmin_pt/F");
    BaselineTree->Branch("Wmin_eta",&Wmin_eta,"Wmin_eta/F");
    BaselineTree->Branch("DPhi_METWmin",&DPhi_METWmin,"DPhi_METWmin/F");
    BaselineTree->Branch("DPhi_WminZ",&DPhi_WminZ,"DPhi_WminZ/F");
    BaselineTree->Branch("EventNumber",&EventNumber,"EventNumber/I");
    BaselineTree->Branch("RunNumber",&RunNumber,"RunNumber/I");
    BaselineTree->Branch("totalWeight",&totalWeight,"totalWeight/D");

    hist_cutflow_raw->SetStats(0);
    hist_cutflow_raw->GetXaxis()->SetBinLabel(1, "2lep");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(2, "flavor");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(3, "trigger");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(4, "OS");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(5, "lep pT");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(6, "njet");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(7, "mll");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(8, "prompt");

    hist_cutflow_weight->SetStats(0);
    hist_cutflow_weight->GetXaxis()->SetBinLabel(1, "2lep");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(2, "flavor");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(3, "trigger");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(4, "OS");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(5, "lep pT");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(6, "njet");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(7, "mll");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(8, "prompt");

    for (int bin=0;bin<bin_size;bin++) {
        hist_METl_Pt[bin] = new TH1D(TString("hist_METl_Pt_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        hist_METl_Pt[bin]->SetStats(0);
        hist_METt_Pt[bin] = new TH1D(TString("hist_METt_Pt_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        hist_METt_Pt[bin]->SetStats(0);
    }

    TH1D* hist_low_njet = new TH1D("hist_low_njet","",bin_size,njet_bin);
    TH1D* hist_low_nbjet = new TH1D("hist_low_nbjet","",bin_size,njet_bin);
    TH1D* hist_low_pt = new TH1D("hist_low_pt","",bin_size,pt_bin);
    TH1D* hist_sm_pt = new TH1D("hist_sm_pt","",bin_size,sm_pt_bin);
    TH1D* hist_low_ht = new TH1D("hist_low_ht","",bin_size,ht_bin);

    //-----------------------------
    // these variables do not go to output
    //-----------------------------
    float n_3jet = 0;
    float Wtruth_corr = 0;
    float Wmin_corr = 0;
    float W12_corr = 0;
    float W80_corr = 0;
    float W01_corr = 0;
    float Wjigsaw_corr = 0;

    //-----------------------------
    // loop over events
    //-----------------------------

    Long64_t nentries = inputTree->GetEntries();

    //nentries = 1000;
    for (Long64_t i=0;i<nentries;i+=event_interval) {

        if (fmod(i,1e5)==0) std::cout << i << " events processed." << std::endl;
        inputTree->GetEntry(i);

        if (MET>100) N_passMET100 += 1; 
        if (jet_n>=3 && MET>150) n_3jet += 1;

        if ( nLep_signal  != 2                 ) continue; // exactly 2 signal leptons
        if ( nLep_base    != 2                 ) continue; // exactly 2 baseline leptons
        if ( lep_pT->at(0) < leading_lep_pt_cut ) continue; // 1st lep pT > 25 GeV
        if ( lep_pT->at(1) < second_lep_pt_cut  ) continue; // 2nd lep pT > 25 GeV

        //--- evaluate weight
        totalWeight = 1;
        if (isData == "MC") totalWeight = genWeight * eventWeight * leptonWeight * jvtWeight * bTagWeight * pileupWeight * FFWeight;

        //--- determine channel
        channel = -1;
        if( ( lepFlavor->at(0) == 1 && lepFlavor->at(1) == 1 ) && trigMatch_1L2LTrig  ) channel = 1; // ee
        if( ( lepFlavor->at(0) == 2 && lepFlavor->at(1) == 2 ) && trigMatch_1L2LTrig  ) channel = 0; // mumu
        if( ( lepFlavor->at(0) == 1 && lepFlavor->at(1) == 2 ) && trigMatch_1L2LTrig  ) channel = 2; // em
        if( ( lepFlavor->at(0) == 2 && lepFlavor->at(1) == 1 ) && trigMatch_1L2LTrig  ) channel = 3; // me

        //--- determine OS / SS
        is_OS = -1;
        if( lepCharge->at(0) != lepCharge->at(1) ) is_OS = 1;
        if( lepCharge->at(0) == lepCharge->at(1) ) is_OS = 0;

        if( channel < 0 ) continue; // require exactly 2 signal leptons and corresponding triggers
        if( is_OS != 1  ) continue; // require opposite-sign
        if( jet_n < 1   ) continue; // require at least 1 pT > 30 GeV jets

        int njet = hist_low_njet->FindBin(jet_n)-1;
        if (jet_n>njet_bin[bin_size]) njet = bin_size-1;
        int nbjet = hist_low_nbjet->FindBin(bjet_n)-1;
        if (bjet_n>njet_bin[bin_size]) nbjet = bin_size-1;
        int pt = hist_low_pt->FindBin(Z_pt)-1;
        int smpt = hist_sm_pt->FindBin(Z_pt)-1;
        if (Z_pt>pt_bin[bin_size]) pt = bin_size-1;
        int ht = hist_low_ht->FindBin(HT)-1;
        if (HT>ht_bin[bin_size]) ht = bin_size-1;

        //---------------------------------------------
        // here we compute the MET parallel and perpendicular components
        // and DR between photon and nearby jet
        // and Oslo's MET_rel
        //---------------------------------------------

        TLorentzVector lep0vec;
        TLorentzVector lep1vec;

        lep0vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
        lep1vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);

        Z_eta = ( lep0vec + lep1vec ).Eta();
        Z_phi = ( lep0vec + lep1vec ).Phi();

        METt = MET*TMath::Sin(MET_phi-Z_phi);
        METl = MET*TMath::Cos(MET_phi-Z_phi);
        TLorentzVector met_4vec;
        met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
        TLorentzVector z_4vec;
        z_4vec.SetPtEtaPhiM(Z_pt,Z_eta,Z_phi,0);
        DPhi_METPhoton = fabs(met_4vec.DeltaPhi(z_4vec));
        TLorentzVector lep0_4vec;
        lep0_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
        TLorentzVector lep1_4vec;
        lep1_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);
        DPhi_2Lep = fabs(lep0_4vec.DeltaPhi(lep1_4vec));
        DR_2Lep = lep0_4vec.DeltaR(lep1_4vec);
        DPhi_METLepLeading = fabs(met_4vec.DeltaPhi(lep0_4vec));
        DPhi_METLepSecond = fabs(met_4vec.DeltaPhi(lep1_4vec));
        DPhi_METLepMin = min(DPhi_METLepLeading,DPhi_METLepSecond);
        TLorentzVector tst_4vec;
        tst_4vec.SetPtEtaPhiM(MET_softTerm,0,MET_softPhi,0);
        DPhi_TSTLepLeading = fabs(tst_4vec.DeltaPhi(lep0_4vec));
        DPhi_TSTLepSecond = fabs(tst_4vec.DeltaPhi(lep1_4vec));
        DPhi_TSTLepMin = min(DPhi_TSTLepLeading,DPhi_TSTLepSecond);
        MinDR_Lep0Jet = 1000.;
        MinDR_Lep1Jet = 1000.;
        MinDR_PhotonJet = 1000.;
        MinDPhi_PhotonJet = 1000.;
        TLorentzVector jet_4vec;
        for (unsigned int j=0;j<jet_pT->size();j++) {
            jet_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),jet_m->at(j));
            float DR_Lep0Jet = jet_4vec.DeltaR(lep0_4vec);
            float DR_Lep1Jet = jet_4vec.DeltaR(lep1_4vec);
            if (MinDR_Lep0Jet>DR_Lep0Jet) MinDR_Lep0Jet = DR_Lep0Jet;
            if (MinDR_Lep1Jet>DR_Lep1Jet) MinDR_Lep1Jet = DR_Lep1Jet;
            float DR_PhotonJet = jet_4vec.DeltaR(z_4vec);
            float DPhi_PhotonJet = jet_4vec.DeltaPhi(z_4vec);
            if (MinDR_PhotonJet>DR_PhotonJet) MinDR_PhotonJet = DR_PhotonJet;
            if (MinDPhi_PhotonJet>DPhi_PhotonJet) MinDPhi_PhotonJet = DPhi_PhotonJet;
        }
        float min_DPhi_MET_LepJet = 1000.;
        float DPhi_MET_LepJet = 1000.;
        for (unsigned int j=0;j<jet_pT->size();j++) {
            jet_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),jet_m->at(j));
            DPhi_MET_LepJet = jet_4vec.DeltaR(met_4vec);
            if (min_DPhi_MET_LepJet>DPhi_MET_LepJet) min_DPhi_MET_LepJet = DPhi_MET_LepJet;
        }
        DPhi_MET_LepJet = lep0_4vec.DeltaR(met_4vec);
        if (min_DPhi_MET_LepJet>DPhi_MET_LepJet) min_DPhi_MET_LepJet = DPhi_MET_LepJet;
        DPhi_MET_LepJet = lep1_4vec.DeltaR(met_4vec);
        if (min_DPhi_MET_LepJet>DPhi_MET_LepJet) min_DPhi_MET_LepJet = DPhi_MET_LepJet;
        MET_rel = MET;
        if (min_DPhi_MET_LepJet<TMath::Pi()/2.) MET_rel = MET*TMath::Sin(min_DPhi_MET_LepJet);

        //---------------------------------------------
        // compute dijet system variables, m80jj, W pT, DR(2jet), etc.
        //---------------------------------------------
        z_4vec.SetPtEtaPhiM(Z_pt,Z_eta,Z_phi,0);
        met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
        GetDijetVariables(z_4vec,met_4vec);

        BaselineTree->Fill();     

    }

    std::cout << "write output..." << std::endl;
    BaselineTree->Write();

    hist_cutflow_raw->Write();
    hist_cutflow_weight->Write();
    for (int bin=0;bin<bin_size;bin++) {
        hist_METl_Pt[bin]->Write();
        hist_METt_Pt[bin]->Write();
    }
    if (isData == "MC") {
        hist_EventCount->SetBinContent(2,nentries);
        hist_EventCount->SetBinContent(3,N_passMET100);
        hist_EventCount->Write();
    }

    std::cout << "done." << std::endl;
    outputFile.Close();
    delete inputFile;

}

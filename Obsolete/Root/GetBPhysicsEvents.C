//-----------------------------------------------------------------------------------------------
// this script takes the photon ntuples from QuickAna and generates smaller ntuples with information required by photon method
// the parameters of the function GetPhotonEvents(string sampleID, string outputName, string pathToNtuples, bool isData) are:
// 	sampleID: DSID of the MC sample
// 	pathToNtuples: the path to the input ntuple
// 	isData: put "true" if you are running a data sample
// example of code running command: root -l -b -q 'GetPhotonEvents.C+("361039","gmc","root://eosatlas//eos/atlas/user/r/rshang/MET_template/v21/photon_mc/",false)'
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
#include "TObject.h"

#include "BasicSetting.C"
#include "InputVariables.C"

using namespace std;

vector<string> noSampleWeight;
vector<string> noEventWeight;

void getPhotonSmearingFunction(TString file, TString histname, TH1D* hist) {
	TFile *fData = TFile::Open( file );
	fData->cd("");
	TH1D* temp = (TH1D*)fData->Get(histname);
	for (int bin=1;bin<=temp->GetNbinsX();bin++) {
		hist->SetBinContent(bin,max(temp->GetBinContent(bin),0.));
	}
	fData->Close();
}
void GetBPhysicsEvents(string sampleID, string outputName, string pathToNtuples, int isData) {

	//---------------------------------------------
	// open file, get Tree and EventCountHist
	//---------------------------------------------

	TH1::SetDefaultSumw2();

	string  filename       = Form("%s%s.root",pathToNtuples.c_str(),sampleID.c_str()); 
	TFile*  inputFile      = TFile::Open(filename.c_str());
	TH1D*   EventCountHist = (TH1D*) inputFile->Get("EventCountHist");
	Float_t _nGenEvents    = EventCountHist->GetBinContent(2);
	TTree*  T              = (TTree*)inputFile->Get("outputTree");

	cout << endl;
	cout << "Opening file           : " << filename        << endl;
	cout << "Events in ntuple       : " << T->GetEntries() << endl;
	cout << "Total generated events : " << _nGenEvents     << endl;

	//-----------------------------
	// access existing branches
	//-----------------------------
	
	
	Long64_t EventNumber;
	Int_t RunNumber;
	//std::vector<double>* HLT_g10_loose = new std::vector<double>(10);
	//std::vector<double>* HLT_g15_loose_L1EM7 = new std::vector<double>(10);
	//std::vector<double>* HLT_g20_loose_L1EM12 = new std::vector<double>(10);
	//std::vector<double>* HLT_g20_loose_L1EM15 = new std::vector<double>(10);
	//std::vector<double>* HLT_g25_loose_L1EM15 = new std::vector<double>(10);
	//std::vector<double>* HLT_g35_loose_L1EM15 = new std::vector<double>(10);
	//std::vector<double>* HLT_g40_loose_L1EM15 = new std::vector<double>(10);
	//std::vector<double>* HLT_g45_loose_L1EM15 = new std::vector<double>(10);
	//std::vector<double>* HLT_g50_loose_L1EM15 = new std::vector<double>(10);
	//std::vector<double>* HLT_g20_loose_larpeb_L1EM15 = new std::vector<double>(10);
	//std::vector<double>* HLT_g40_loose_larpeb = new std::vector<double>(10);
	//std::vector<double>* HLT_g60_loose = new std::vector<double>(10);
	//std::vector<double>* HLT_g70_loose = new std::vector<double>(10);
	//std::vector<double>* HLT_g80_loose = new std::vector<double>(10);
	//std::vector<double>* HLT_g100_loose = new std::vector<double>(10);
	//std::vector<double>* HLT_g120_loose = new std::vector<double>(10);
	//std::vector<double>* HLT_g140_loose = new std::vector<double>(10);
	Float_t Mu;
	double MET;
	double MET_phi;
	double MET_softTerm;
	double MET_softPhi;
	Float_t truthMET;
	Float_t truthMET_Phi;
	double DPhi_METJetLeading;
	double DPhi_METJetSecond;
	double HT;
	Int_t jet_n;
	Int_t bjet_n;
	Int_t lep_n;
	std::vector<double>* lep_pT = new std::vector<double>(10);
	std::vector<double>* lep_eta = new std::vector<double>(10);
	std::vector<double>* lep_phi = new std::vector<double>(10);
	std::vector<double>* jet_btag = new std::vector<double>(10);
	std::vector<double>* jet_m = new std::vector<double>(10);
	std::vector<double>* jet_pT = new std::vector<double>(10);
	std::vector<double>* jet_eta = new std::vector<double>(10);
	std::vector<double>* jet_phi = new std::vector<double>(10);
	Int_t Jpsi_is_OS;
	double Jpsi_mll;
	double Jpsi_pt;
	double Jpsi_eta;
	double Jpsi_phi;
    T->SetBranchStatus("*", 0);
    T->SetBranchStatus("EventNumber"              ,1); 
    T->SetBranchStatus("RunNumber"              ,1); 
    T->SetBranchStatus("Mu"              ,1); 
    T->SetBranchStatus("MET"             ,1); 
    T->SetBranchStatus("MET_phi"         ,1); 
    T->SetBranchStatus("MET_softTerm"         ,1); 
    T->SetBranchStatus("MET_softPhi"         ,1); 
    T->SetBranchStatus("truthMET"             ,1); 
    T->SetBranchStatus("truthMET_Phi"         ,1); 
    T->SetBranchStatus("DPhi_METJetLeading"         ,1); 
    T->SetBranchStatus("DPhi_METJetSecond"          ,1); 
    T->SetBranchStatus("HT"              ,1); 
    T->SetBranchStatus("bjet_n"          ,1); 
    T->SetBranchStatus("lep_n"           ,1); 
    T->SetBranchStatus("lep_pT"          ,1); 
    T->SetBranchStatus("lep_eta"         ,1); 
    T->SetBranchStatus("lep_phi"         ,1); 
    T->SetBranchStatus("jet_n"           ,1); 
    T->SetBranchStatus("jet_btag"          ,1); 
    T->SetBranchStatus("jet_m"          ,1); 
    T->SetBranchStatus("jet_pT"          ,1); 
    T->SetBranchStatus("jet_eta"         ,1); 
    T->SetBranchStatus("jet_phi"         ,1); 
    //T->SetBranchStatus("HLT_g10_loose", 1); 
    //T->SetBranchStatus("HLT_g15_loose_L1EM7", 1); 
    //T->SetBranchStatus("HLT_g20_loose_L1EM12", 1); 
    //T->SetBranchStatus("HLT_g20_loose_L1EM15", 1); 
    //T->SetBranchStatus("HLT_g25_loose_L1EM15", 1); 
    //T->SetBranchStatus("HLT_g35_loose_L1EM15", 1); 
    //T->SetBranchStatus("HLT_g40_loose_L1EM15", 1); 
    //T->SetBranchStatus("HLT_g45_loose_L1EM15", 1); 
    //T->SetBranchStatus("HLT_g50_loose_L1EM15", 1); 
    //T->SetBranchStatus("HLT_g20_loose_larpeb_L1EM15", 1); 
    //T->SetBranchStatus("HLT_g40_loose_larpeb", 1); 
    //T->SetBranchStatus("HLT_g60_loose", 1); 
    //T->SetBranchStatus("HLT_g70_loose", 1); 
    //T->SetBranchStatus("HLT_g80_loose", 1); 
    //T->SetBranchStatus("HLT_g100_loose", 1); 
    //T->SetBranchStatus("HLT_g120_loose", 1); 
    //T->SetBranchStatus("HLT_g140_loose", 1); 
    T->SetBranchStatus("is_OS", 1); 
    T->SetBranchStatus("mll", 1); 
    T->SetBranchStatus("Z_pt", 1); 
    T->SetBranchStatus("Z_eta", 1); 
    T->SetBranchStatus("Z_phi", 1); 

	T->SetBranchAddress("EventNumber"              ,&EventNumber               );
	T->SetBranchAddress("RunNumber"              ,&RunNumber               );
	T->SetBranchAddress("Mu"              ,&Mu               );
	T->SetBranchAddress("MET"             ,&MET              );
	T->SetBranchAddress("MET_phi"         ,&MET_phi          );
	T->SetBranchAddress("MET_softTerm"         ,&MET_softTerm          );
	T->SetBranchAddress("MET_softPhi"         ,&MET_softPhi          );
	T->SetBranchAddress("truthMET"             ,&truthMET              );
	T->SetBranchAddress("truthMET_Phi"         ,&truthMET_Phi          );
	T->SetBranchAddress("DPhi_METJetLeading"         ,&DPhi_METJetLeading          );
	T->SetBranchAddress("DPhi_METJetSecond"          ,&DPhi_METJetSecond           );
	T->SetBranchAddress("HT"              ,&HT               );
	T->SetBranchAddress("bjet_n"          ,&bjet_n            );
	T->SetBranchAddress("lep_n"           ,&lep_n            );
	T->SetBranchAddress("lep_pT"          ,&lep_pT           );
	T->SetBranchAddress("lep_eta"         ,&lep_eta          );
	T->SetBranchAddress("lep_phi"         ,&lep_phi          );
	T->SetBranchAddress("jet_n"           ,&jet_n            );
	T->SetBranchAddress("jet_btag"          ,&jet_btag           );
	T->SetBranchAddress("jet_m"          ,&jet_m           );
	T->SetBranchAddress("jet_pT"          ,&jet_pT           );
	T->SetBranchAddress("jet_eta"         ,&jet_eta          );
	T->SetBranchAddress("jet_phi"         ,&jet_phi          );
	//T->SetBranchAddress("HLT_g10_loose", &HLT_g10_loose);
	//T->SetBranchAddress("HLT_g15_loose_L1EM7", &HLT_g15_loose_L1EM7);
	//T->SetBranchAddress("HLT_g20_loose_L1EM12", &HLT_g20_loose_L1EM12);
	//T->SetBranchAddress("HLT_g20_loose_L1EM15", &HLT_g20_loose_L1EM15);
	//T->SetBranchAddress("HLT_g25_loose_L1EM15", &HLT_g25_loose_L1EM15);
	//T->SetBranchAddress("HLT_g35_loose_L1EM15", &HLT_g35_loose_L1EM15);
	//T->SetBranchAddress("HLT_g40_loose_L1EM15", &HLT_g40_loose_L1EM15);
	//T->SetBranchAddress("HLT_g45_loose_L1EM15", &HLT_g45_loose_L1EM15);
	//T->SetBranchAddress("HLT_g50_loose_L1EM15", &HLT_g50_loose_L1EM15);
	//T->SetBranchAddress("HLT_g20_loose_larpeb_L1EM15", &HLT_g20_loose_larpeb_L1EM15);
	//T->SetBranchAddress("HLT_g40_loose_larpeb", &HLT_g40_loose_larpeb);
	//T->SetBranchAddress("HLT_g60_loose", &HLT_g60_loose);
	//T->SetBranchAddress("HLT_g70_loose", &HLT_g70_loose);
	//T->SetBranchAddress("HLT_g80_loose", &HLT_g80_loose);
	//T->SetBranchAddress("HLT_g100_loose", &HLT_g100_loose);
	//T->SetBranchAddress("HLT_g120_loose", &HLT_g120_loose);
	//T->SetBranchAddress("HLT_g140_loose", &HLT_g140_loose);
	T->SetBranchAddress("is_OS", &Jpsi_is_OS);
	T->SetBranchAddress("mll", &Jpsi_mll);
	T->SetBranchAddress("Z_pt", &Jpsi_pt);
	T->SetBranchAddress("Z_eta", &Jpsi_eta);
	T->SetBranchAddress("Z_phi", &Jpsi_phi);

	double sampleWeight;
	Float_t eventWeight;
	if (isData!=1) {
	        T->SetBranchStatus("sampleWeight",1);
	        T->SetBranchStatus("eventWeight",1);
		T->SetBranchAddress("sampleWeight",&sampleWeight);
		T->SetBranchAddress("eventWeight",&eventWeight);
	}

	//-----------------------------
	// add new branches
	//-----------------------------

	TFile   outputFile(TString(outputPath)+"/"+TString(outputName)+"/"+TString(sampleID.c_str())+".root","recreate");
	TTree BaselineTree("BaselineTree","baseline tree");
	double METl = 0.;
	double METt = 0.;
	double DPhi_METJpsi = 0.;
	double MinDR_JpsiJet = 0.;
	double MinDPhi_JpsiJet = 0.;
	int pt = 0;
	int ht = 0;
	int njet = 0;
	int nbjet = 0;
	double MT;
	BaselineTree.Branch("pt",&pt,"pt/I");
	BaselineTree.Branch("ht",&ht,"ht/I");
	BaselineTree.Branch("njet",&njet,"njet/I");
	BaselineTree.Branch("nbjet",&nbjet,"nbjet/I");
	BaselineTree.Branch("Mu",&Mu,"Mu/Float_t");
	BaselineTree.Branch("MET_raw",&MET,"MET_raw/D");
	BaselineTree.Branch("METl_raw",&METl,"METl_raw/D");
	BaselineTree.Branch("METt_raw",&METt,"METt_raw/D");
	BaselineTree.Branch("MET_phi_raw",&MET_phi,"MET_phi_raw/D");
	//BaselineTree.Branch("MET_softTerm",&MET_softTerm,"MET_softTerm/D");
	//BaselineTree.Branch("MET_softPhi",&MET_softPhi,"MET_softPhi/D");
	//BaselineTree.Branch("truthMET",&truthMET,"truthMET/Float_t");
	//BaselineTree.Branch("truthMETl",&truthMETl,"truthMETl/D");
	//BaselineTree.Branch("truthMETt",&truthMETt,"truthMETt/D");
	//BaselineTree.Branch("truthMET_Phi",&MET_phi,"MET_phi/Float_t");
	//BaselineTree.Branch("DPhi_METJetLeading_raw",&DPhi_METJetLeading,"DPhi_METJetLeading_raw/D");
	//BaselineTree.Branch("DPhi_METJetSecond_raw",&DPhi_METJetSecond,"DPhi_METJetSecond_raw/D");
	//BaselineTree.Branch("DPhi_METPhoton_raw",&DPhi_METPhoton,"DPhi_METPhoton_raw/D");
	//BaselineTree.Branch("MinDR_PhotonJet",&MinDR_PhotonJet,"MinDR_PhotonJet/D");
	//BaselineTree.Branch("MinDPhi_PhotonJet",&MinDPhi_PhotonJet,"MinDPhi_PhotonJet/D");
	BaselineTree.Branch("HT",&HT,"HT/F");
	BaselineTree.Branch("MT",&MT,"MT/D");
	BaselineTree.Branch("Jpsi_mll",&Jpsi_mll,"Jpsi_mll/D");
	BaselineTree.Branch("Jpsi_pt",&Jpsi_pt,"Jpsi_pt/D");
	BaselineTree.Branch("Jpsi_eta",&Jpsi_eta,"Jpsi_eta/D");
	BaselineTree.Branch("Jpsi_phi",&Jpsi_phi,"Jpsi_phi/D");
	BaselineTree.Branch("bjet_n",&bjet_n,"bjet_n/Int_t");
	BaselineTree.Branch("jet_n",&jet_n,"jet_n/Int_t");
	BaselineTree.Branch("jet_btag","std::vector<double>",&jet_btag);
	BaselineTree.Branch("jet_m","std::vector<double>",&jet_m);
	BaselineTree.Branch("jet_pT","std::vector<double>",&jet_pT);
	BaselineTree.Branch("jet_phi","std::vector<double>",&jet_phi);
	BaselineTree.Branch("jet_eta","std::vector<double>",&jet_eta);
	BaselineTree.Branch("EventNumber",&EventNumber,"EventNumber/Long64_t");
	BaselineTree.Branch("RunNumber",&RunNumber,"RunNumber/Int_t");
	BaselineTree.Branch("lep_n_raw",&lep_n,"lep_n_raw/Int_t");
	BaselineTree.Branch("lep_pT_raw","std::vector<double>",&lep_pT);
	BaselineTree.Branch("lep_phi_raw","std::vector<double>",&lep_phi);
	BaselineTree.Branch("lep_eta_raw","std::vector<double>",&lep_eta);
	double totalWeight = 0.;
	BaselineTree.Branch("totalWeight",&totalWeight,"totalWeight/D");

	TH1D* hist_low_njet = new TH1D("hist_low_njet","",bin_size,njet_bin);
	hist_low_njet->SetStats(0);
	TH1D* hist_low_nbjet = new TH1D("hist_low_nbjet","",bin_size,njet_bin);
	hist_low_nbjet->SetStats(0);
	TH1D* hist_low_pt = new TH1D("hist_low_pt","",bin_size,bphys_pt_bin);
	hist_low_pt->SetStats(0);
	TH1D* hist_sm_pt = new TH1D("hist_sm_pt","",bin_size,bphys_pt_bin);
	hist_sm_pt->SetStats(0);
	TH1D* hist_low_et = new TH1D("hist_low_et","",bin_size,bphys_pt_bin);
	hist_low_et->SetStats(0);
	TH1D* hist_low_ht = new TH1D("hist_low_ht","",bin_size,ht_bin);
	hist_low_ht->SetStats(0);
	TH1D* hist_medium_njet = new TH1D("hist_medium_njet","",bin_size,njet_bin);
	hist_medium_njet->SetStats(0);
	TH1D* hist_medium_nbjet = new TH1D("hist_medium_nbjet","",bin_size,njet_bin);
	hist_medium_nbjet->SetStats(0);
	TH1D* hist_medium_pt = new TH1D("hist_medium_pt","",bin_size,bphys_pt_bin);
	hist_medium_pt->SetStats(0);
	TH1D* hist_medium_et = new TH1D("hist_medium_et","",bin_size,bphys_pt_bin);
	hist_medium_et->SetStats(0);
	TH1D* hist_medium_ht = new TH1D("hist_medium_ht","",bin_size,ht_bin);
	hist_medium_ht->SetStats(0);
	TH1D* hist_high_njet = new TH1D("hist_high_njet","",bin_size,njet_bin);
	hist_high_njet->SetStats(0);
	TH1D* hist_high_nbjet = new TH1D("hist_high_nbjet","",bin_size,njet_bin);
	hist_high_nbjet->SetStats(0);
	TH1D* hist_high_pt = new TH1D("hist_high_pt","",bin_size,bphys_pt_bin);
	hist_high_pt->SetStats(0);
	TH1D* hist_high_et = new TH1D("hist_high_et","",bin_size,bphys_pt_bin);
	hist_high_et->SetStats(0);
	TH1D* hist_high_ht = new TH1D("hist_high_ht","",bin_size,ht_bin);
	hist_high_ht->SetStats(0);
	TH1D* hist_zmet_njet = new TH1D("hist_zmet_njet","",bin_size,njet_bin);
	hist_zmet_njet->SetStats(0);
	TH1D* hist_zmet_nbjet = new TH1D("hist_zmet_nbjet","",bin_size,njet_bin);
	hist_zmet_nbjet->SetStats(0);
	TH1D* hist_zmet_pt = new TH1D("hist_zmet_pt","",bin_size,bphys_pt_bin);
	hist_zmet_pt->SetStats(0);
	TH1D* hist_zmet_et = new TH1D("hist_zmet_et","",bin_size,bphys_pt_bin);
	hist_zmet_et->SetStats(0);
	TH1D* hist_zmet_ht = new TH1D("hist_zmet_ht","",bin_size,ht_bin);
	hist_zmet_ht->SetStats(0);
	TH1D* hist_bveto_njet = new TH1D("hist_bveto_njet","",bin_size,njet_bin);
	hist_bveto_njet->SetStats(0);
	TH1D* hist_bveto_nbjet = new TH1D("hist_bveto_nbjet","",bin_size,njet_bin);
	hist_bveto_nbjet->SetStats(0);
	TH1D* hist_bveto_pt = new TH1D("hist_bveto_pt","",bin_size,bphys_pt_bin);
	hist_bveto_pt->SetStats(0);
	TH1D* hist_bveto_et = new TH1D("hist_bveto_et","",bin_size,bphys_pt_bin);
	hist_bveto_et->SetStats(0);
	TH1D* hist_bveto_ht = new TH1D("hist_bveto_ht","",bin_size,ht_bin);
	hist_bveto_ht->SetStats(0);


	TLorentzVector obj_4vec;
	TLorentzVector isr_4vec;
	TLorentzVector z_4vec;
	//-----------------------------
	// loop over events
	//-----------------------------

	Long64_t nentries = T->GetEntries();
	for (Long64_t i=0;i<nentries;i+=event_interval) {
	//for (Long64_t i=0;i<100000;i++) {

		if (fmod(i,1e5)==0) std::cout << i << " events processed." << std::endl;
		T->GetEntry(i);

		if (lep_n<2) continue;
		if (jet_n==0) continue;
		if (channel!=0 && channel!=1) continue;
		if (Jpsi_is_OS!=1) continue;
		if (abs(Jpsi_mll-3.)>0.5 && abs(Jpsi_mll-9.5)>1.5) continue;

		// find the trigger and prescale
		double trigWeight = 0;

		//if (HLT_g15_loose_L1EM7->at(1)==1 && photon_pT->at(0)>(15+5) && photon_pT->at(0)<(20+5)) trigWeight = HLT_g15_loose_L1EM7->at(0);
		//if (HLT_g20_loose_L1EM12->at(1)==1 && photon_pT->at(0)>(20+5) && photon_pT->at(0)<(25+5)) trigWeight = HLT_g20_loose_L1EM12->at(0);
		//if (HLT_g25_loose_L1EM15->at(1)==1 && photon_pT->at(0)>(25+5) && photon_pT->at(0)<(35+5)) trigWeight = HLT_g25_loose_L1EM15->at(0);
		//if (HLT_g35_loose_L1EM15->at(1)==1 && photon_pT->at(0)>(35+5) && photon_pT->at(0)<(40+5)) trigWeight = HLT_g35_loose_L1EM15->at(0);
		//if (HLT_g40_loose_L1EM15->at(1)==1 && photon_pT->at(0)>(40+5) && photon_pT->at(0)<(45+5)) trigWeight = HLT_g40_loose_L1EM15->at(0);
		//if (HLT_g45_loose_L1EM15->at(1)==1 && photon_pT->at(0)>(45+5) && photon_pT->at(0)<(50+5)) trigWeight = HLT_g45_loose_L1EM15->at(0);
		//if (HLT_g50_loose_L1EM15->at(1)==1 && photon_pT->at(0)>(50+5) && photon_pT->at(0)<(60+5)) trigWeight = HLT_g50_loose_L1EM15->at(0);
		//if (HLT_g60_loose->at(1)==1 && photon_pT->at(0)>(60+5) && photon_pT->at(0)<(70+5)) trigWeight = HLT_g60_loose->at(0);
		//if (HLT_g70_loose->at(1)==1 && photon_pT->at(0)>(70+5) && photon_pT->at(0)<(80+5)) trigWeight = HLT_g70_loose->at(0);
		//if (HLT_g80_loose->at(1)==1 && photon_pT->at(0)>(80+5) && photon_pT->at(0)<(100+5)) trigWeight = HLT_g80_loose->at(0);
		//if (HLT_g100_loose->at(1)==1 && photon_pT->at(0)>(100+5) && photon_pT->at(0)<(140+5)) trigWeight = HLT_g100_loose->at(0);
		//if (HLT_g140_loose->at(1)==1 && photon_pT->at(0)>(140+5)) trigWeight = HLT_g140_loose->at(0);
		//if (trigWeight==0) continue;

		njet = hist_low_njet->FindBin(jet_n)-1;
		if (jet_n>njet_bin[bin_size]) njet = bin_size-1;
		nbjet = hist_low_nbjet->FindBin(bjet_n)-1;
		if (bjet_n>njet_bin[bin_size]) nbjet = bin_size-1;
		pt = hist_low_pt->FindBin(Jpsi_pt)-1;
		if (Jpsi_pt>bphys_pt_bin[bin_size]) pt = bin_size-1;
		ht = hist_low_ht->FindBin(HT)-1;
		if (HT>ht_bin[bin_size]) ht = bin_size-1;

		//totalWeight = trigWeight;
		totalWeight = 1;

		// here we compute the MET parallel and perpendicular components
		METt = MET*TMath::Sin(MET_phi-Jpsi_phi);
		METl = MET*TMath::Cos(MET_phi-Jpsi_phi);

		TLorentzVector Jpsi_4vec;
		Jpsi_4vec.SetPtEtaPhiM(Jpsi_pt,Jpsi_eta,Jpsi_phi,0);
		// hist_dPt_Pt histogram, i.e. the photon truth-reco response function
		if (isData!=1) {

			totalWeight = (sampleWeight*eventWeight)/_nGenEvents;

		}
		totalWeight = totalWeight*event_interval;

		DPhi_METJpsi = fabs(TMath::ATan2(METt,METl));
		MinDR_JpsiJet = 1000.;
		MinDPhi_JpsiJet = 1000.;
		TLorentzVector jet_4vec;
		for (int j=0;j<jet_pT->size();j++) {
			jet_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),jet_m->at(j));
			double DR_JpsiJet = jet_4vec.DeltaR(Jpsi_4vec);
			double DPhi_JpsiJet = jet_4vec.DeltaPhi(Jpsi_4vec);
			if (MinDR_JpsiJet>DR_JpsiJet) MinDR_JpsiJet = DR_JpsiJet;
			if (MinDPhi_JpsiJet>DPhi_JpsiJet) MinDPhi_JpsiJet = DPhi_JpsiJet;
		}

		TLorentzVector met_4vec;
		met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
		TLorentzVector lep_4vec;
		if (lep_pT->size()>0) lep_4vec.SetPtEtaPhiM(lep_pT->at(0),0,lep_phi->at(0),0);  // only transverse component
		else lep_4vec.SetPtEtaPhiM(0,0,0,0);
		if (lep_pT->size()>0) MT = (met_4vec+lep_4vec).M();
		else MT = 0;



		double resample_size = 1.;
		if (MET>200.) {
			if (totalWeight>1 && isData==1) {
				totalWeight = totalWeight/resample_size;
				for (int r=0;r<resample_size;r++) BaselineTree.Fill();
			}
			else if (Jpsi_pt>200. && outputName=="gmc") {
				totalWeight = totalWeight/resample_size;
				for (int r=0;r<resample_size;r++) BaselineTree.Fill();
			}
			else if (outputName=="Vg") {
				totalWeight = totalWeight/resample_size;
				for (int r=0;r<resample_size;r++) BaselineTree.Fill();
			}
			else BaselineTree.Fill();
		}
		else BaselineTree.Fill();

	}

	BaselineTree.Write();

	std::cout << "done." << std::endl;
	outputFile.Close();
	delete inputFile;


}

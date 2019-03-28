
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
#include "TSpectrum.h"
#include "TVirtualFFT.h"

TString tag;


TTree* inputTree;
Long64_t EventNumber;
Int_t RunNumber;
double totalWeight = 0.;
double mll;
double MET;
double MET_phi;
double METl = 0.;
double METt = 0.;
double lep0_truthPt= 0.;
double lep0_truthEta = 0.;
double lep0_truthPhi = 0.;
double lep1_truthPt= 0.;
double lep1_truthEta = 0.;
double lep1_truthPhi = 0.;
Int_t lep_n;
Int_t jet_n;
Int_t bjet_n;
double HT = 0.;
double Z_truthPt = 0.;
double Z_truthEta = 0.;
double Z_truthPhi = 0.;
double gamma_pt = 0.;
double gamma_eta = 0.;
double gamma_phi = 0.;
double Z_pt = 0.;
Float_t MT2W = 0.;
double Z_eta = 0.;
double Z_phi = 0.;
std::vector<double>* jet_pT = new std::vector<double>(10);
std::vector<double>* jet_eta = new std::vector<double>(10);
std::vector<double>* jet_phi = new std::vector<double>(10);
std::vector<double>* jet_m = new std::vector<double>(10);
std::vector<double>* jet_btag = new std::vector<double>(10);
std::vector<double>* lep_pT = new std::vector<double>(10);
std::vector<double>* lep_eta = new std::vector<double>(10);
std::vector<double>* lep_phi = new std::vector<double>(10);
double ptrw = 0;


void SetBasicBranch() {
	inputTree->SetBranchStatus("EventNumber"		,1);
	inputTree->SetBranchAddress("EventNumber"		,&EventNumber);
	inputTree->SetBranchStatus("RunNumber"			,1);
	inputTree->SetBranchAddress("RunNumber"			,&RunNumber);
	inputTree->SetBranchStatus("totalWeight"		,1);
	inputTree->SetBranchAddress("totalWeight"		,&totalWeight);
	inputTree->SetBranchStatus("MET"			,1);
	inputTree->SetBranchAddress("MET"			,&MET);
	inputTree->SetBranchStatus("MET_phi"			,1);
	inputTree->SetBranchAddress("MET_phi"			,&MET_phi);
	inputTree->SetBranchStatus("METl"			,1);
	inputTree->SetBranchAddress("METl"			,&METl);
	inputTree->SetBranchStatus("METt"			,1);
	inputTree->SetBranchAddress("METt"			,&METt);
	inputTree->SetBranchStatus("mll"			,1);
	inputTree->SetBranchAddress("mll"			,&mll);
	inputTree->SetBranchStatus("HT"				,1);
	inputTree->SetBranchAddress("HT"			,&HT);
	inputTree->SetBranchStatus("lep_n"			,1);
	inputTree->SetBranchAddress("lep_n"			,&lep_n);
	inputTree->SetBranchStatus("jet_n"			,1);
	inputTree->SetBranchAddress("jet_n"			,&jet_n);
	inputTree->SetBranchStatus("bjet_n"			,1);
	inputTree->SetBranchAddress("bjet_n"			,&bjet_n);
	inputTree->SetBranchStatus("gamma_pt"			,1);
	inputTree->SetBranchAddress("gamma_pt"			,&gamma_pt);
	inputTree->SetBranchStatus("gamma_eta"			,1);
	inputTree->SetBranchAddress("gamma_eta"			,&gamma_eta);
	inputTree->SetBranchStatus("gamma_phi"			,1);
	inputTree->SetBranchAddress("gamma_phi"			,&gamma_phi);
	inputTree->SetBranchStatus("MT2W"			,1);
	inputTree->SetBranchAddress("MT2W"			,&MT2W);
	inputTree->SetBranchStatus("Z_pt"			,1);
	inputTree->SetBranchAddress("Z_pt"			,&Z_pt);
	inputTree->SetBranchStatus("Z_eta"			,1);
	inputTree->SetBranchAddress("Z_eta"			,&Z_eta);
	inputTree->SetBranchStatus("Z_phi"			,1);
	inputTree->SetBranchAddress("Z_phi"			,&Z_phi);
	inputTree->SetBranchStatus("jet_pT"			,1);
	inputTree->SetBranchAddress("jet_pT"			,&jet_pT);
	inputTree->SetBranchStatus("jet_eta"			,1);
	inputTree->SetBranchAddress("jet_eta"			,&jet_eta);
	inputTree->SetBranchStatus("jet_phi"			,1);
	inputTree->SetBranchAddress("jet_phi"			,&jet_phi);
	inputTree->SetBranchStatus("jet_m"			,1);
	inputTree->SetBranchAddress("jet_m"			,&jet_m);
	inputTree->SetBranchStatus("lep_pT"			,1);
	inputTree->SetBranchAddress("lep_pT"			,&lep_pT);
	inputTree->SetBranchStatus("lep_eta"			,1);
	inputTree->SetBranchAddress("lep_eta"			,&lep_eta);
	inputTree->SetBranchStatus("lep_phi"			,1);
	inputTree->SetBranchAddress("lep_phi"			,&lep_phi);
	inputTree->SetBranchStatus("Z_truthPt"			,1);
	inputTree->SetBranchAddress("Z_truthPt"			,&Z_truthPt);
	inputTree->SetBranchStatus("Z_truthPhi"			,1);
	inputTree->SetBranchAddress("Z_truthPhi"		,&Z_truthPhi);
	inputTree->SetBranchStatus("Z_truthEta"			,1);
	inputTree->SetBranchAddress("Z_truthEta"		,&Z_truthEta);
}
double JetEtMissDPhi(double unit, int j) {
	if (j>=jet_pT->size()) return 1e10;
	TLorentzVector met_4vec;
	met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
	TLorentzVector jet_4vec;
	double DPhiMETJet = 0;
	jet_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),0);
	DPhiMETJet = abs(jet_4vec.DeltaPhi(met_4vec));
	return DPhiMETJet;
}
bool MonoJetSelection(double unit) {
	if (jet_n!=1) return false;
	if (bjet_n>0) return false;
	if (lep_pT->at(0)<25.*unit) return false;
	if (lep_pT->at(1)<25.*unit) return false;
	if (mll/unit<12) return false;
	if (Z_pt/unit<50) return false;
	if (MET/unit<50) return false;
	return true;
}
bool Strong2LSelection(double unit) {
	//if (lep_n!=2) return false;
	if (jet_n<2) return false;
	if (lep_pT->at(0)<50.*unit) return false;
	if (lep_pT->at(1)<25.*unit) return false;
	if (mll/unit<12) return false;
	return true;
}
bool Strong2LSRhigh(double unit) {
	double DPhi_METJet0 = JetEtMissDPhi(unit,0);
	double DPhi_METJet1 = JetEtMissDPhi(unit,1);
	double DPhi_METJetMin = min(DPhi_METJet0,DPhi_METJet1);
	if (DPhi_METJetMin<0.4) return false;
	if (HT/unit<1200) return false;
	return true;
}
bool Strong2LVRhigh(double unit) {
	double DPhi_METJet0 = JetEtMissDPhi(unit,0);
	double DPhi_METJet1 = JetEtMissDPhi(unit,1);
	double DPhi_METJetMin = min(DPhi_METJet0,DPhi_METJet1);
	if (DPhi_METJetMin>0.4) return false;
	if (HT/unit<1200) return false;
	return true;
}
bool Strong2LSRmedium(double unit) {
	double DPhi_METJet0 = JetEtMissDPhi(unit,0);
	double DPhi_METJet1 = JetEtMissDPhi(unit,1);
	double DPhi_METJetMin = min(DPhi_METJet0,DPhi_METJet1);
	if (DPhi_METJetMin<0.4) return false;
	if (HT/unit<400) return false;
	if (MT2W/unit<25) return false;
	return true;
}
bool Strong2LVRmedium(double unit) {
	double DPhi_METJet0 = JetEtMissDPhi(unit,0);
	double DPhi_METJet1 = JetEtMissDPhi(unit,1);
	double DPhi_METJetMin = min(DPhi_METJet0,DPhi_METJet1);
	if (DPhi_METJetMin>0.4) return false;
	if (HT/unit<400) return false;
	if (MT2W/unit<25) return false;
	return true;
}
bool Strong2LSRlow(double unit) {
	double DPhi_METJet0 = JetEtMissDPhi(unit,0);
	double DPhi_METJet1 = JetEtMissDPhi(unit,1);
	double DPhi_METJetMin = min(DPhi_METJet0,DPhi_METJet1);
	if (DPhi_METJetMin<0.4) return false;
	if (HT/unit<200) return false;
	if (MT2W/unit<70) return false;
	return true;
}
bool Strong2LVRlow(double unit) {
	double DPhi_METJet0 = JetEtMissDPhi(unit,0);
	double DPhi_METJet1 = JetEtMissDPhi(unit,1);
	double DPhi_METJetMin = min(DPhi_METJet0,DPhi_METJet1);
	if (DPhi_METJetMin>0.4) return false;
	if (HT/unit<200) return false;
	if (MT2W/unit<70) return false;
	return true;
}
double JetEtMissSig(double unit, int j) {
	if (j>=jet_pT->size()) return 1e10;
	TLorentzVector met_4vec;
	met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
	TLorentzVector jet_4vec;
	double DPhiMETJet = 0;
	double sig = 0.;
	jet_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),0);
	DPhiMETJet = abs(jet_4vec.DeltaPhi(met_4vec));
	sig = abs(jet_pT->at(j)/unit*TMath::Cos(DPhiMETJet));
	sig = pow(sig,0.5);
	return MET/sig;
}
double PhotonEtMissDPhi(double unit) {
	TLorentzVector met_4vec;
	met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
	TLorentzVector lept1_4vec;
	lept1_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
	TLorentzVector lept2_4vec;
	lept2_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);
	TLorentzVector z_4vec = lept1_4vec+lept2_4vec;
	double DPhiMETPhoton = 0;
	DPhiMETPhoton = abs(z_4vec.DeltaPhi(met_4vec));
	return DPhiMETPhoton;
}
double PhotonEtMissSig(double unit) {
	TLorentzVector met_4vec;
	met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
	TLorentzVector lept1_4vec;
	lept1_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
	TLorentzVector lept2_4vec;
	lept2_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);
	TLorentzVector z_4vec = lept1_4vec+lept2_4vec;
	double DPhiMETPhoton = 0;
	double sig = 0.;
	DPhiMETPhoton = abs(z_4vec.DeltaPhi(met_4vec));
	sig = abs(z_4vec.Pt()/unit*TMath::Cos(DPhiMETPhoton));
	return MET/sig;
}
double PhotonEtMissRatio(double unit) {
	TLorentzVector lept1_4vec;
	lept1_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
	TLorentzVector lept2_4vec;
	lept2_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);
	TLorentzVector z_4vec = lept1_4vec+lept2_4vec;
	return MET/z_4vec.Pt();
}
bool BadMuon(double unit) {
	if (MET/unit<100) return false;
	if (Z_pt/unit<400) return false;
	TLorentzVector met_4vec;
	met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
	TLorentzVector jet_4vec;
	double DPhiMETJet = 0;
	double RatioMETJet = 0;
	for (int j=0;j<jet_pT->size();j++) {
		jet_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),0);
		DPhiMETJet = abs(jet_4vec.DeltaPhi(met_4vec));
		RatioMETJet = (MET/unit)/pow(abs((jet_pT->at(j)/unit)*TMath::Cos(DPhiMETJet)),0.5);
		if (DPhiMETJet<0.4) return false;
		//if (RatioMETJet<15) return false;
	}
	TLorentzVector lept1_4vec;
	lept1_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
	TLorentzVector lept2_4vec;
	lept2_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);
	TLorentzVector z_4vec = lept1_4vec+lept2_4vec;
	double DPhiMETPhoton = abs(z_4vec.DeltaPhi(met_4vec));
	double RatioMETPhoton = (MET/unit)/abs((z_4vec.Pt()/unit)*TMath::Cos(DPhiMETPhoton));
	if (DPhiMETPhoton>0.4) return false;
	//if (DPhiMETPhoton>0.4 && DPhiMETPhoton<2.8) return false;
	//if (RatioMETPhoton<0.5) return false;
	//if (RatioMETPhoton>1.0) return false;
	return true;
}
void Strong2L(TString pathToNtuple, TString pathToOutput, TString TreeName, bool DoReweighting, double scale, double unit) {

	TH1::SetDefaultSumw2();

	std::cout << "open file " << pathToNtuple << std::endl;
	TFile*  inputFile      = TFile::Open(pathToNtuple);
	inputTree              = (TTree*)inputFile->Get(TreeName);
	inputTree->SetBranchStatus("*", 0);
	SetBasicBranch();
	ptrw = 1.;
	Z_truthPt = 0;
	Z_truthEta = 0;
	Z_truthPhi = 0;
	if (DoReweighting) {
		if (tag.Contains("PTRW")) {
			inputTree->SetBranchStatus("ptrw_bveto"			,1);
			inputTree->SetBranchAddress("ptrw_bveto"		,&ptrw);
			if (tag.Contains("low")) inputTree->SetBranchStatus("ptrw_ht200"			,1);
			if (tag.Contains("low")) inputTree->SetBranchAddress("ptrw_ht200"		,&ptrw);
			if (tag.Contains("medium")) inputTree->SetBranchStatus("ptrw_ht400"			,1);
			if (tag.Contains("medium")) inputTree->SetBranchAddress("ptrw_ht400"		,&ptrw);
			if (tag.Contains("high")) inputTree->SetBranchStatus("ptrw_ht1200"		,1);
			if (tag.Contains("high")) inputTree->SetBranchAddress("ptrw_ht1200"		,&ptrw);
		}
		else if (tag.Contains("HTRW")) {
			inputTree->SetBranchStatus("htrw_bveto"			,1);
			inputTree->SetBranchAddress("htrw_bveto"		,&ptrw);
		}
	}
	else {
		ptrw = 1.;
	}

	Long64_t nentries = inputTree->GetEntries();

	std::cout << "output file " << pathToOutput << std::endl;
	TFile   outputFile(pathToOutput,"recreate");

	Float_t bins[] = { 0,20,40,60,80,100,150,200,300,400,500 };
	Float_t bins_mll[] = { 12,41,61,81,101,151,201,251,301,401,501,1001 };
	if (tag.Contains("low")) bins_mll = { 12,41,61,81,101,151,201,251,301,401,501,1001 };
	if (tag.Contains("high")) bins_mll = { 12,41,61,81,101,151,201,251,301,401,1001 };
	if (tag.Contains("high")) bins_mll = { 12,41,61,81,101,151,201,251,301,401,501,1001 };
	Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1;
	Int_t  binnum_mll = sizeof(bins_mll)/sizeof(Float_t) - 1;
	TH1D Hist_MET		= TH1D(TString("Hist_MET"),"",binnum,bins);
	TH1D Hist_MET_20GeV	= TH1D(TString("Hist_MET_20GeV"),"",binnum,bins);
	TH1D Hist_HT		= TH1D(TString("Hist_HT"),"",15,0,1500);
	TH1D Hist_ZPt		= TH1D(TString("Hist_ZPt"),"",15,0,1500);
	//TH1D Hist_mll		= TH1D(TString("Hist_mll"),"",5,41,141);
	TH1D Hist_mll		= TH1D(TString("Hist_mll"),"",binnum_mll,bins_mll);
	TH1D Hist_nlep		= TH1D(TString("Hist_nlep"),"",5,0,5);
	TH1D Hist_DPhiMETJet1st	= TH1D(TString("Hist_DPhiMETJet1st"),"",16,0,3.2);
	TH1D Hist_DPhiMETJet2nd	= TH1D(TString("Hist_DPhiMETJet2nd"),"",16,0,3.2);
	TH1D Hist_DPhiMETJetMin	= TH1D(TString("Hist_DPhiMETJetMin"),"",16,0,3.2);
	TH1D Hist_LowMET_HT		= TH1D(TString("Hist_LowMET_HT"),"",15,0,1500);
	TH1D Hist_LowMET_mll		= TH1D(TString("Hist_LowMET_mll"),"",binnum_mll,bins_mll);
	TH1D Hist_LowMET_nlep		= TH1D(TString("Hist_LowMET_nlep"),"",5,0,5);
	TH1D Hist_LowMET_DPhiMETJet1st	= TH1D(TString("Hist_LowMET_DPhiMETJet1st"),"",16,0,3.2);
	TH1D Hist_LowMET_DPhiMETJet2nd	= TH1D(TString("Hist_LowMET_DPhiMETJet2nd"),"",16,0,3.2);
	TH1D Hist_LowMET_DPhiMETJetMin	= TH1D(TString("Hist_LowMET_DPhiMETJetMin"),"",16,0,3.2);

	int step = 1;
	for (Long64_t i=0;i<nentries;i+=step) {
		if (fmod(i,1e5)==0) std::cout << pathToNtuple << " " << i << "/" << nentries << " events processed." << std::endl;
		inputTree->GetEntry(i);
		if (!Strong2LSelection(unit)) continue;
		if (tag.Contains("Strong2LSRlow")) {if (!Strong2LSRlow(unit)) continue;}
		else if (tag.Contains("Strong2LVRlow")) {if (!Strong2LVRlow(unit)) continue;}
		else if (tag.Contains("Strong2LSRmedium")) {if (!Strong2LSRmedium(unit)) continue;}
		else if (tag.Contains("Strong2LVRmedium")) {if (!Strong2LVRmedium(unit)) continue;}
		else if (tag.Contains("Strong2LSRhigh")) {if (!Strong2LSRhigh(unit)) continue;}
		else if (tag.Contains("Strong2LVRhigh")) {if (!Strong2LVRhigh(unit)) continue;}

		
		if (tag.Contains("VgOff")) {
			if (RunNumber>=301890 && RunNumber<=301910) continue;
		}

		double DPhiMETJet1st = JetEtMissDPhi(unit,0);
		double DPhiMETJet2nd = JetEtMissDPhi(unit,1);
		double DPhiMETJetMin = min(DPhiMETJet1st,DPhiMETJet2nd);

		int met_bin = Hist_MET_20GeV.FindBin(MET/unit);
		double bin_width = Hist_MET_20GeV.GetBinLowEdge(met_bin+1)-Hist_MET_20GeV.GetBinLowEdge(met_bin);
		double bin_scale = 20./bin_width;
		//double bin_scale = 1;
		Hist_MET_20GeV.Fill(TMath::Min(MET/unit,Hist_MET_20GeV.GetBinCenter(Hist_MET_20GeV.GetNbinsX())),totalWeight*ptrw*scale*bin_scale);
		Hist_MET.Fill(TMath::Min(MET/unit,Hist_MET.GetBinCenter(Hist_MET.GetNbinsX())),totalWeight*ptrw*scale);

		if (MET/unit>100. && MET/unit<200.) {
			//std::cout << "totalWeight = " << totalWeight << ", ptrw = " << ptrw << ", scale = " << scale << std::endl;
			Hist_LowMET_mll.Fill(TMath::Min(mll/unit,Hist_LowMET_mll.GetBinCenter(Hist_LowMET_mll.GetNbinsX())),totalWeight*ptrw*scale);
			if (mll/unit>81 && mll/unit<101) {
				Hist_LowMET_HT.Fill(HT/unit,totalWeight*ptrw*scale);
				Hist_LowMET_nlep.Fill(lep_n,totalWeight*ptrw*scale);
				Hist_LowMET_DPhiMETJet1st.Fill(DPhiMETJet1st,totalWeight*ptrw*scale);
				Hist_LowMET_DPhiMETJet2nd.Fill(DPhiMETJet2nd,totalWeight*ptrw*scale);
				Hist_LowMET_DPhiMETJetMin.Fill(DPhiMETJetMin,totalWeight*ptrw*scale);
			}
		} 

		if (tag.Contains("low")) if (MET/unit<250.) continue;
		if (tag.Contains("medium")) if (MET/unit<400.) continue;
		if (tag.Contains("high")) if (MET/unit<200.) continue;
		Hist_mll.Fill(TMath::Min(mll/unit,Hist_mll.GetBinCenter(Hist_mll.GetNbinsX())),totalWeight*ptrw*scale);
		if (mll/unit>81 && mll/unit<101) {
			Hist_HT.Fill(HT/unit,totalWeight*ptrw*scale);
			Hist_nlep.Fill(lep_n,totalWeight*ptrw*scale);
			Hist_DPhiMETJet1st.Fill(DPhiMETJet1st,totalWeight*ptrw*scale);
			Hist_DPhiMETJet2nd.Fill(DPhiMETJet2nd,totalWeight*ptrw*scale);
			Hist_DPhiMETJetMin.Fill(DPhiMETJetMin,totalWeight*ptrw*scale);
		}
	}

	Hist_MET.Write();
	Hist_MET_20GeV.Write();
	Hist_HT.Write();
	Hist_mll.Write();
	Hist_nlep.Write();
	Hist_DPhiMETJet1st.Write();
	Hist_DPhiMETJet2nd.Write();
	Hist_DPhiMETJetMin.Write();
	Hist_LowMET_HT.Write();
	Hist_LowMET_mll.Write();
	Hist_LowMET_nlep.Write();
	Hist_LowMET_DPhiMETJet1st.Write();
	Hist_LowMET_DPhiMETJet2nd.Write();
	Hist_LowMET_DPhiMETJetMin.Write();
	std::cout << "done." << std::endl;
	outputFile.Close();
	delete inputFile;

}
void RunStrong2L() {

	bool DoReweighting;
	double Scale = 36.1*1000.;
	TString pathToNtuple;
	TString pathToOutput;
	TString OutputNtuples_dir = TString("/eos/atlas/user/r/rshang/OutputNtuples/");

	pathToNtuple = OutputNtuples_dir+TString("zjets_tt/ztt_ee.root");
	pathToOutput = TString("../OutputHistogram/Ztt_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);

	pathToNtuple = OutputNtuples_dir+TString("zjets_tt/ztt_mm.root");
	pathToOutput = TString("../OutputHistogram/Ztt_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);

	pathToNtuple = OutputNtuples_dir+TString("tt/ttee.root");
	pathToOutput = TString("../OutputHistogram/Top_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);

	pathToNtuple = OutputNtuples_dir+TString("vv/vvee.root");
	pathToOutput = TString("../OutputHistogram/VV_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);

	pathToNtuple = OutputNtuples_dir+TString("data/data_ee.root");
	pathToOutput = TString("../OutputHistogram/Data_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,1.0,1.0);

	pathToNtuple = OutputNtuples_dir+TString("data/data_em.root");
	pathToOutput = TString("../OutputHistogram/DataFS_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,0.5,1.0);

	pathToNtuple = OutputNtuples_dir+TString("tt/ttmm.root");
	pathToOutput = TString("../OutputHistogram/Top_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);

	pathToNtuple = OutputNtuples_dir+TString("vv/vvmm.root");
	pathToOutput = TString("../OutputHistogram/VV_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);

	pathToNtuple = OutputNtuples_dir+TString("data/data_mm.root");
	pathToOutput = TString("../OutputHistogram/Data_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,1.0,1.0);

	pathToNtuple = OutputNtuples_dir+TString("data/data_em.root");
	pathToOutput = TString("../OutputHistogram/DataFS_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,0.5,1.0);

	pathToNtuple = OutputNtuples_dir+TString("gdata/gdata_ee_NoSmear.root");
	pathToOutput = TString("../OutputHistogram/GData_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = true;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,1.0,1.0);

	if (tag.Contains("McSmear"))
		pathToNtuple = OutputNtuples_dir+TString("gdata/gdata_mm_McSmear.root");
	else if (tag.Contains("NoSmear"))
		pathToNtuple = OutputNtuples_dir+TString("gdata/gdata_mm_NoSmear.root");
	else if (tag.Contains("DataSmear"))
		pathToNtuple = OutputNtuples_dir+TString("gdata/gdata_mm_DataSmear.root");
	pathToOutput = TString("../OutputHistogram/GData_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = true;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,1.0,1.0);

	pathToNtuple = OutputNtuples_dir+TString("zjets_221/zjets_ee.root");
	pathToOutput = TString("../OutputHistogram/ZMC221_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);

	pathToNtuple = OutputNtuples_dir+TString("zjets_221/zjets_mm.root");
	pathToOutput = TString("../OutputHistogram/ZMC221_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);


	pathToNtuple = OutputNtuples_dir+TString("gmc/gmc_ee_NoSmear.root");
	pathToOutput = TString("../OutputHistogram/GMC_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = true;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	
	if (tag.Contains("McSmear"))
		pathToNtuple = OutputNtuples_dir+TString("gmc/gmc_mm_McSmear.root");
	else if (tag.Contains("NoSmear"))
		pathToNtuple = OutputNtuples_dir+TString("gmc/gmc_mm_NoSmear.root");
	pathToOutput = TString("../OutputHistogram/GMC_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = true;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	
	pathToNtuple = OutputNtuples_dir+TString("zjets/zjets_ee.root");
	pathToOutput = TString("../OutputHistogram/ZMC_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	
	pathToNtuple = OutputNtuples_dir+TString("zjets/zjets_mm.root");
	pathToOutput = TString("../OutputHistogram/ZMC_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	

}
void MakeSelectionHistogram_Strong2LPaper2017() {

	//RunStrong2L();
	//
	tag = TString("Strong2LVRhigh_PTRW_McSmear");  // the nominal prediction (pT-reweighting and MC smearing) in VR-high
	RunStrong2L();
	tag = TString("Strong2LSRhigh_PTRW_McSmear");  // the nominal prediction in SR-high
	RunStrong2L();
	tag = TString("Strong2LSRhigh_PTRW_McSmear_VgOff");  // the prediction without V-gamma subtraction in SR-high
	RunStrong2L();
	tag = TString("Strong2LVRhigh_PTRW_McSmear_VgOff");
	RunStrong2L();
	tag = TString("Strong2LVRhigh_HTRW_McSmear");  // the prediction with HT-reweighting in VR-high
	RunStrong2L();
	tag = TString("Strong2LVRhigh_PTRW_DataSmear"); // the prediction with data-driven smearing in VR-high
	RunStrong2L();

	tag = TString("Strong2LSRlow_PTRW_McSmear");
	RunStrong2L();
	tag = TString("Strong2LSRlow_PTRW_McSmear_VgOff");
	RunStrong2L();
	tag = TString("Strong2LVRlow_PTRW_McSmear");
	RunStrong2L();
	tag = TString("Strong2LVRlow_PTRW_McSmear_VgOff");
	RunStrong2L();
	tag = TString("Strong2LVRlow_HTRW_McSmear");
	RunStrong2L();
	tag = TString("Strong2LVRlow_PTRW_DataSmear");
	RunStrong2L();
	
	tag = TString("Strong2LSRmedium_PTRW_McSmear");
	RunStrong2L();
	tag = TString("Strong2LSRmedium_PTRW_McSmear_VgOff");
	RunStrong2L();
	tag = TString("Strong2LVRmedium_PTRW_McSmear");
	RunStrong2L();
	tag = TString("Strong2LVRmedium_PTRW_McSmear_VgOff");
	RunStrong2L();
	tag = TString("Strong2LVRmedium_HTRW_McSmear");
	RunStrong2L();
	tag = TString("Strong2LVRmedium_PTRW_DataSmear");
	RunStrong2L();
	
}

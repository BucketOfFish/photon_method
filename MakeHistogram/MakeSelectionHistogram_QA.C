
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
#include "TProfile.h"
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
#include "MT2_ROOT.h"
#include "TRandom.h"

TString tag;
//TString tag = TString("MET0");
//TString tag = TString("MET50");
//TString tag = TString("MET100");
//TString tag = TString("MET150");
//TString tag = TString("MET200");
//TString tag = TString("GoodJet");
//TString tag = TString("BadMuon");
//TString tag = TString("CompressWino");
//TString tag = TString("LowWino");
//TString tag = TString("LowPtLep");
//TString tag = TString("RJ2CA_Compress");
//TString tag = TString("RJ2CA_Low");
//TString tag = TString("RJ2CA_All");

TRandom myRandom;
TTree* inputTree;
Long64_t EventNumber;
Int_t RunNumber;
double totalWeight = 0.;
double mll;
double MET;
double MET_phi;
double MET_softTerm;
double MET_softPhi;
double METl = 0.;
double METt = 0.;
double lep0_truthPt= 0.;
double lep0_truthEta = 0.;
double lep0_truthPhi = 0.;
double lep1_truthPt= 0.;
double lep1_truthEta = 0.;
double lep1_truthPhi = 0.;
Int_t jet_n;
Int_t lep_n;
Int_t bjet_n;
double HT = 0.;
double Z_truthPt = 0.;
double Z_truthEta = 0.;
double Z_truthPhi = 0.;
double gamma_pt = 0.;
double gamma_eta = 0.;
double gamma_phi = 0.;
double Z_pt = 0.;
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
	inputTree->SetBranchStatus("MET_softTerm"		,1);
	inputTree->SetBranchAddress("MET_softTerm"		,&MET_softTerm);
	inputTree->SetBranchStatus("MET_softPhi"		,1);
	inputTree->SetBranchAddress("MET_softPhi"		,&MET_softPhi);
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
double GetMT2LorentzTransform(TLorentzVector boost_4vec, TLorentzVector lep0_4vec, TLorentzVector lep1_4vec, TLorentzVector met_4vec) {

	TVector3 boost_vec = boost_4vec.BoostVector();
	lep0_4vec.Boost(boost_vec);
	lep1_4vec.Boost(boost_vec);
	met_4vec.Boost(boost_vec);
	double MT2_new = ComputeMT2(lep0_4vec, lep1_4vec, met_4vec, 0, 0).Compute();
	return MT2_new;

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
	if (abs(lep_eta->at(0))>2.5) return false;
	if (abs(lep_eta->at(1))>2.5) return false;
	if (mll/unit<12) return false;
	if (Z_pt/unit<50) return false;
	if (HT/unit<100) return false;
	//if (HT/unit<200) return false;
	double DPhiMETJet1st = JetEtMissDPhi(unit,0);
	if (DPhiMETJet1st<0.4) return false;
	return true;
}
bool LowMllExcess2JetSelection(double unit) {
	if (jet_n!=2) return false;
	if (bjet_n>0) return false;
	if (lep_pT->at(0)<25.*unit) return false;
	if (lep_pT->at(1)<25.*unit) return false;
	if (abs(lep_eta->at(0))>2.5) return false;
	if (abs(lep_eta->at(1))>2.5) return false;
	if (mll/unit<12) return false;
	if (Z_pt/unit<50) return false;
	if (jet_pT->at(0)/unit<100) return false;
	double DPhiMETJet1st = JetEtMissDPhi(unit,0);
	if (DPhiMETJet1st<0.4) return false;
	return true;
}
bool LowMllExcess3JetSelection(double unit) {
	if (jet_n!=3) return false;
	if (bjet_n>0) return false;
	if (lep_pT->at(0)<25.*unit) return false;
	if (lep_pT->at(1)<25.*unit) return false;
	if (abs(lep_eta->at(0))>2.5) return false;
	if (abs(lep_eta->at(1))>2.5) return false;
	if (mll/unit<12) return false;
	if (Z_pt/unit<50) return false;
	if (jet_pT->at(0)/unit<100) return false;
	double DPhiMETJet1st = JetEtMissDPhi(unit,0);
	if (DPhiMETJet1st<0.4) return false;
	return true;
}
bool LowMllExcess4JetSelection(double unit) {
	if (jet_n!=4) return false;
	if (bjet_n>0) return false;
	if (lep_pT->at(0)<25.*unit) return false;
	if (lep_pT->at(1)<25.*unit) return false;
	if (abs(lep_eta->at(0))>2.5) return false;
	if (abs(lep_eta->at(1))>2.5) return false;
	if (mll/unit<12) return false;
	if (Z_pt/unit<50) return false;
	if (jet_pT->at(0)/unit<100) return false;
	double DPhiMETJet1st = JetEtMissDPhi(unit,0);
	if (DPhiMETJet1st<0.4) return false;
	return true;
}
bool LowMllExcess5JetSelection(double unit) {
	if (jet_n!=5) return false;
	if (bjet_n>0) return false;
	if (lep_pT->at(0)<25.*unit) return false;
	if (lep_pT->at(1)<25.*unit) return false;
	if (abs(lep_eta->at(0))>2.5) return false;
	if (abs(lep_eta->at(1))>2.5) return false;
	if (mll/unit<12) return false;
	if (Z_pt/unit<50) return false;
	if (jet_pT->at(0)/unit<100) return false;
	double DPhiMETJet1st = JetEtMissDPhi(unit,0);
	if (DPhiMETJet1st<0.4) return false;
	return true;
}
bool LowMllExcessBJetSelection(double unit) {
	if (jet_n!=1) return false;
	if (bjet_n!=1) return false;
	if (lep_pT->at(0)<25.*unit) return false;
	if (lep_pT->at(1)<25.*unit) return false;
	if (abs(lep_eta->at(0))>2.5) return false;
	if (abs(lep_eta->at(1))>2.5) return false;
	if (mll/unit<12) return false;
	if (Z_pt/unit<50) return false;
	if (jet_pT->at(0)/unit<100) return false;
	double DPhiMETJet1st = JetEtMissDPhi(unit,0);
	if (DPhiMETJet1st<0.4) return false;
	return true;
}
bool Strong2LSelection(double unit) {
	if (jet_n<2) return false;
	if (bjet_n>0) return false;
	if (lep_pT->at(0)<25.*unit) return false;
	if (lep_pT->at(1)<25.*unit) return false;
	if (abs(lep_eta->at(0))>2.5) return false;
	if (abs(lep_eta->at(1))>2.5) return false;
	if (jet_pT->at(0)<30.*unit) return false;
	if (jet_pT->at(1)<30.*unit) return false;
	if (mll/unit<12) return false;
	if (Z_pt/unit<50) return false;
	if (HT/unit<200) return false;
	return true;
}
bool LowPtZSelection(double unit) {
	if (jet_n<2) return false;
	if (bjet_n>0) return false;
	if (lep_pT->at(0)<25.*unit) return false;
	if (lep_pT->at(1)<25.*unit) return false;
	if (abs(lep_eta->at(0))>2.5) return false;
	if (abs(lep_eta->at(1))>2.5) return false;
	if (jet_pT->at(0)<30.*unit) return false;
	if (jet_pT->at(1)<30.*unit) return false;
	if (mll/unit<12) return false;
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
	//sig = abs(jet_pT->at(j)/unit*TMath::Cos(DPhiMETJet));
	sig = jet_pT->at(j)/unit;
	sig = pow(sig,0.5);
	return (MET/unit)/sig;
}
double MTWino(double unit, int j1, int j2) {
	if (j1>=jet_pT->size()) return 1e10;
	if (j2>=jet_pT->size()) return 1e10;
	TLorentzVector met_4vec;
	met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
	TLorentzVector jet1_4vec;
	TLorentzVector jet2_4vec;
	TLorentzVector w_4vec;
	jet1_4vec.SetPtEtaPhiM(jet_pT->at(j1),jet_eta->at(j1),jet_phi->at(j1),jet_m->at(j1));
	jet2_4vec.SetPtEtaPhiM(jet_pT->at(j2),jet_eta->at(j2),jet_phi->at(j2),jet_m->at(j2));
	w_4vec = jet1_4vec + jet2_4vec;
	TLorentzVector lept1_4vec;
	lept1_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
	TLorentzVector lept2_4vec;
	lept2_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);
	TLorentzVector z_4vec = lept1_4vec+lept2_4vec;
	//w_4vec.SetPtEtaPhiM(w_4vec.Pt(),w_4vec.Eta(),w_4vec.Phi(),0.);
	//z_4vec.SetPtEtaPhiM(z_4vec.Pt(),z_4vec.Eta(),z_4vec.Phi(),0.);
	double mt = ComputeMT2(z_4vec, w_4vec, met_4vec, 0, 0).Compute();
	return mt/unit;
}
double MT2W(double unit) {
	TLorentzVector met_4vec;
	TLorentzVector lep1_4vec;
	TLorentzVector lep2_4vec;
	lep1_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
	lep2_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);
	met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
	double mt = ComputeMT2(lep1_4vec, lep2_4vec, met_4vec, 0, 0).Compute();
	return mt/unit;
}
double JetInvMass(double unit, int j1, int j2) {
	if (j1>=jet_pT->size()) return 1e10;
	if (j2>=jet_pT->size()) return 1e10;
	TLorentzVector jet1_4vec;
	TLorentzVector jet2_4vec;
	TLorentzVector w_4vec;
	jet1_4vec.SetPtEtaPhiM(jet_pT->at(j1),jet_eta->at(j1),jet_phi->at(j1),jet_m->at(j1));
	jet2_4vec.SetPtEtaPhiM(jet_pT->at(j2),jet_eta->at(j2),jet_phi->at(j2),jet_m->at(j2));
	w_4vec = jet1_4vec + jet2_4vec;
	return w_4vec.M()/unit;
}
double JetInvMassMin(double unit) {
	TLorentzVector jet1_4vec;
	TLorentzVector jet2_4vec;
	TLorentzVector w_4vec;
	double min_mass = 1e10;
	for (int j1=0;j1<jet_pT->size()-1;j1++) {
		jet1_4vec.SetPtEtaPhiM(jet_pT->at(j1),jet_eta->at(j1),jet_phi->at(j1),jet_m->at(j1));
		for (int j2=j1+1;j2<jet_pT->size();j2++) {
			jet2_4vec.SetPtEtaPhiM(jet_pT->at(j2),jet_eta->at(j2),jet_phi->at(j2),jet_m->at(j2));
			w_4vec = jet1_4vec + jet2_4vec;
			if (min_mass>w_4vec.M()/unit) min_mass = w_4vec.M()/unit;
		}
	}
	return min_mass;
}
double JetInvMassMax(double unit) {
	TLorentzVector jet1_4vec;
	TLorentzVector jet2_4vec;
	TLorentzVector w_4vec;
	double max_mass = 0;
	for (int j1=0;j1<jet_pT->size()-1;j1++) {
		jet1_4vec.SetPtEtaPhiM(jet_pT->at(j1),jet_eta->at(j1),jet_phi->at(j1),jet_m->at(j1));
		for (int j2=j1+1;j2<jet_pT->size();j2++) {
			jet2_4vec.SetPtEtaPhiM(jet_pT->at(j2),jet_eta->at(j2),jet_phi->at(j2),jet_m->at(j2));
			w_4vec = jet1_4vec + jet2_4vec;
			if (max_mass<w_4vec.M()/unit) max_mass = w_4vec.M()/unit;
		}
	}
	return max_mass;
}
std::pair<int, int> FindJetPairCloseTo80GeV(double unit)
{
	TLorentzVector jet1_4vec;
	TLorentzVector jet2_4vec;
	TLorentzVector w_4vec;
	double min_diff  =1e10;
	int jet1 = 0;
	int jet2 = 0;
	for (int j1=0;j1<jet_pT->size()-1;j1++) {
		jet1_4vec.SetPtEtaPhiM(jet_pT->at(j1),jet_eta->at(j1),jet_phi->at(j1),jet_m->at(j1));
		for (int j2=j1+1;j2<jet_pT->size();j2++) {
			jet2_4vec.SetPtEtaPhiM(jet_pT->at(j2),jet_eta->at(j2),jet_phi->at(j2),jet_m->at(j2));
			w_4vec = jet1_4vec + jet2_4vec;
			if (min_diff>abs(w_4vec.M()/unit-80.)) {
				min_diff = abs(w_4vec.M()/unit-80.);
				jet1 = j1;
				jet2 = j2;
			}
		}
	}
	return std::make_pair(jet1,jet2);
}
TLorentzVector W4vecWithJetPairCloseTo80GeV(double unit)
{
	TLorentzVector jet1_4vec;
	TLorentzVector jet2_4vec;
	TLorentzVector w_4vec;
	double min_diff  =1e10;
	int jet1 = 0;
	int jet2 = 0;
	for (int j1=0;j1<jet_pT->size()-1;j1++) {
		jet1_4vec.SetPtEtaPhiM(jet_pT->at(j1),jet_eta->at(j1),jet_phi->at(j1),jet_m->at(j1));
		for (int j2=j1+1;j2<jet_pT->size();j2++) {
			jet2_4vec.SetPtEtaPhiM(jet_pT->at(j2),jet_eta->at(j2),jet_phi->at(j2),jet_m->at(j2));
			w_4vec = jet1_4vec + jet2_4vec;
			if (min_diff>abs(w_4vec.M()/unit-80.)) {
				min_diff = abs(w_4vec.M()/unit-80.);
				jet1 = j1;
				jet2 = j2;
			}
		}
	}
	jet1_4vec.SetPtEtaPhiM(jet_pT->at(jet1),jet_eta->at(jet1),jet_phi->at(jet1),jet_m->at(jet1));
	jet2_4vec.SetPtEtaPhiM(jet_pT->at(jet2),jet_eta->at(jet2),jet_phi->at(jet2),jet_m->at(jet2));
	return jet1_4vec+jet2_4vec;
}
double MTWinoUsing80GeVPair(double unit) {
	std::pair<int, int> jets;
	jets = FindJetPairCloseTo80GeV(unit);
	return MTWino(unit, jets.first, jets.second);
}
double MTWinoMax(double unit, int j1, int j2, double inv_mass) {
	if (j1>=jet_pT->size()) return 1e10;
	if (j2>=jet_pT->size()) return 1e10;
	TLorentzVector met_4vec;
	met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
	TLorentzVector jet1_4vec;
	TLorentzVector jet2_4vec;
	TLorentzVector w_4vec;
	jet1_4vec.SetPtEtaPhiM(jet_pT->at(j1),jet_eta->at(j1),jet_phi->at(j1),jet_m->at(j1));
	jet2_4vec.SetPtEtaPhiM(jet_pT->at(j2),jet_eta->at(j2),jet_phi->at(j2),jet_m->at(j2));
	w_4vec = jet1_4vec + jet2_4vec;
	TLorentzVector lept1_4vec;
	lept1_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
	TLorentzVector lept2_4vec;
	lept2_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);
	TLorentzVector z_4vec = lept1_4vec+lept2_4vec;
	w_4vec.SetPtEtaPhiM(w_4vec.Pt(),w_4vec.Eta(),w_4vec.Phi(),0.);
	z_4vec.SetPtEtaPhiM(z_4vec.Pt(),z_4vec.Eta(),z_4vec.Phi(),0.);
	double boost_pt;
	double boost_phi;
	TLorentzVector boost_4vec;
	if (inv_mass==0.) inv_mass = (w_4vec+z_4vec).M();
	//double inv_mass = (w_4vec+z_4vec).M();
	double beta_limit = 100.;
	double boost_pt_final = 0;
	double boost_phi_final = 0;
	double MT2_this = 0;
	double MT2_max = 0;
	for (int i=0;i<100;i++) {
		boost_phi = myRandom.Rndm()*2.*TMath::Pi()-TMath::Pi()+met_4vec.Phi();
		boost_4vec.SetPtEtaPhiM(beta_limit/2.,0.,boost_phi,1.);
		met_4vec.SetPtEtaPhiM(met_4vec.Pt(),-(w_4vec+z_4vec).Eta(),met_4vec.Phi(),inv_mass);
		MT2_this = GetMT2LorentzTransform(boost_4vec, w_4vec, z_4vec, met_4vec);
		if (MT2_max<MT2_this) {
			MT2_max = MT2_this;
			boost_phi_final = boost_phi;
		}
	}
	for (int i=0;i<100;i++) {
		boost_pt = myRandom.Rndm()*beta_limit;
		boost_4vec.SetPtEtaPhiM(boost_pt,0.,boost_phi_final,1.);
		met_4vec.SetPtEtaPhiM(met_4vec.Pt(),-(w_4vec+z_4vec).Eta(),met_4vec.Phi(),inv_mass);
		MT2_this = GetMT2LorentzTransform(boost_4vec, w_4vec, z_4vec, met_4vec);
		if (MT2_max<MT2_this) {
			MT2_max = MT2_this;
			boost_pt_final = boost_pt;
		}
	}
	return MT2_max/unit;
}
double MTWinoMaxUsing80GeVPair(double unit, int n_itr) {
	std::pair<int, int> jets;
	jets = FindJetPairCloseTo80GeV(unit);
	double result = MTWinoMax(unit, jets.first, jets.second,0.);
	if (n_itr>=1) max(result,MTWinoMax(unit, jets.first, jets.second,result));
	if (n_itr>=2) max(result,MTWinoMax(unit, jets.first, jets.second,result));
	if (n_itr>=3) max(result,MTWinoMax(unit, jets.first, jets.second,result));
	return result;
}
double JetInvMassUsing80GeVPair(double unit) {
	std::pair<int, int> jets;
	jets = FindJetPairCloseTo80GeV(unit);
	return JetInvMass(unit, jets.first, jets.second);
}
double IsrJetPtUsing80GeVPair(double unit) {
	std::pair<int, int> jets;
	TLorentzVector jet_4vec;
	TLorentzVector isr_4vec;
	jets = FindJetPairCloseTo80GeV(unit);
	for (int j1=0;j1<jet_pT->size()-1;j1++) {
		if (j1!=jets.first && j1!=jets.second) {
			jet_4vec.SetPtEtaPhiM(jet_pT->at(j1),jet_eta->at(j1),jet_phi->at(j1),jet_m->at(j1));
			isr_4vec += jet_4vec;
		}
	}
	return isr_4vec.Pt()/unit;;
}
double WJetPtUsing80GeVPair(double unit) {
	std::pair<int, int> jets;
	TLorentzVector jet1_4vec;
	TLorentzVector jet2_4vec;
	TLorentzVector w_4vec;
	jets = FindJetPairCloseTo80GeV(unit);
	jet1_4vec.SetPtEtaPhiM(jet_pT->at(jets.first),jet_eta->at(jets.first),jet_phi->at(jets.first),jet_m->at(jets.first));
	jet2_4vec.SetPtEtaPhiM(jet_pT->at(jets.second),jet_eta->at(jets.second),jet_phi->at(jets.second),jet_m->at(jets.second));
	w_4vec = jet1_4vec+jet2_4vec;
	return w_4vec.Pt()/unit;;
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
double WEtMissRatio(double unit, int j1, int j2) {
	if (j1>=jet_pT->size()) return 1e10;
	if (j2>=jet_pT->size()) return 1e10;
	TLorentzVector jet1_4vec;
	TLorentzVector jet2_4vec;
	TLorentzVector w_4vec;
	jet1_4vec.SetPtEtaPhiM(jet_pT->at(j1),jet_eta->at(j1),jet_phi->at(j1),jet_m->at(j1));
	jet2_4vec.SetPtEtaPhiM(jet_pT->at(j2),jet_eta->at(j2),jet_phi->at(j2),jet_m->at(j2));
	w_4vec = jet1_4vec + jet2_4vec;
	return MET/w_4vec.Pt();
}
double RJ2CA_HT5PP(double unit) {
	double ht5pp = 0;
	ht5pp = lep_pT->at(0)+lep_pT->at(1)/unit;
	ht5pp += jet_pT->at(0)+jet_pT->at(1)/unit;
	ht5pp += jet_pT->at(0)+jet_pT->at(1)/unit;
	ht5pp += MET/unit;
	return ht5pp;
}
double RJ2CA_PTPP(double unit) {
	double ptpp = 0;
	TLorentzVector all_4vec;
	TLorentzVector lep1_4vec;
	TLorentzVector lep2_4vec;
	TLorentzVector jet1_4vec;
	TLorentzVector jet2_4vec;
	TLorentzVector met_4vec;
	lep1_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
	lep2_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);
	jet1_4vec.SetPtEtaPhiM(jet_pT->at(0),jet_eta->at(0),jet_phi->at(0),jet_m->at(0));
	jet2_4vec.SetPtEtaPhiM(jet_pT->at(1),jet_eta->at(1),jet_phi->at(1),jet_m->at(1));
	met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
	all_4vec = lep1_4vec+lep2_4vec+jet1_4vec+jet2_4vec+met_4vec;
	ptpp = all_4vec.Pt()/unit;
	return ptpp;
}
double RJ2CA_H2PP(double unit) {
	double h2pp = 0;
	TLorentzVector lep1_4vec;
	TLorentzVector lep2_4vec;
	TLorentzVector jet1_4vec;
	TLorentzVector jet2_4vec;
	TLorentzVector met_4vec;
	lep1_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
	lep2_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);
	jet1_4vec.SetPtEtaPhiM(jet_pT->at(0),jet_eta->at(0),jet_phi->at(0),jet_m->at(0));
	jet2_4vec.SetPtEtaPhiM(jet_pT->at(1),jet_eta->at(1),jet_phi->at(1),jet_m->at(1));
	met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
	h2pp = (lep1_4vec+lep2_4vec+jet1_4vec+jet2_4vec).Pt()+met_4vec.Pt();
	h2pp = h2pp/unit;
	return h2pp;
}
double RJ2CA_H5PP(double unit) {
	double h5pp = 0;
	TLorentzVector lep1_4vec;
	TLorentzVector lep2_4vec;
	TLorentzVector jet1_4vec;
	TLorentzVector jet2_4vec;
	TLorentzVector met_4vec;
	lep1_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
	lep2_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);
	jet1_4vec.SetPtEtaPhiM(jet_pT->at(0),jet_eta->at(0),jet_phi->at(0),jet_m->at(0));
	jet2_4vec.SetPtEtaPhiM(jet_pT->at(1),jet_eta->at(1),jet_phi->at(1),jet_m->at(1));
	met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
	h5pp = lep1_4vec.P()+lep2_4vec.P()+jet1_4vec.P()+jet2_4vec.P()+met_4vec.Pt();
	h5pp = h5pp/unit;
	return h5pp;
}
double RJ2CA_PTCM(double unit) {
	double ptcm = 0;
	TLorentzVector lep1_4vec;
	TLorentzVector lep2_4vec;
	TLorentzVector jet1_4vec;
	TLorentzVector jet2_4vec;
	TLorentzVector jet3_4vec;
	TLorentzVector jet4_4vec;
	TLorentzVector met_4vec;
	lep1_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
	lep2_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);
	jet1_4vec.SetPtEtaPhiM(jet_pT->at(0),jet_eta->at(0),jet_phi->at(0),jet_m->at(0));
	jet2_4vec.SetPtEtaPhiM(jet_pT->at(1),jet_eta->at(1),jet_phi->at(1),jet_m->at(1));
	if (jet_pT->size()>=3) jet3_4vec.SetPtEtaPhiM(jet_pT->at(2),jet_eta->at(2),jet_phi->at(2),jet_m->at(2));
	if (jet_pT->size()>=4) jet4_4vec.SetPtEtaPhiM(jet_pT->at(3),jet_eta->at(3),jet_phi->at(3),jet_m->at(3));
	met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
	ptcm = (lep1_4vec+lep2_4vec+jet1_4vec+jet2_4vec+jet3_4vec+jet4_4vec+met_4vec).Pt();
	ptcm = ptcm/unit;
	return ptcm;
}
std::pair<int, int> RJ2CA_FindJetPairFromW(double unit)
{
	TLorentzVector jet0_4vec;
	TLorentzVector jet1_4vec;
	int jet1 = 0;
	int jet2 = 0;
	double max_dphi = 0;
	jet0_4vec.SetPtEtaPhiM(jet_pT->at(0),jet_eta->at(0),jet_phi->at(0),jet_m->at(0));
	for (int j1=1;j1<jet_pT->size();j1++) {
		jet1_4vec.SetPtEtaPhiM(jet_pT->at(j1),jet_eta->at(j1),jet_phi->at(j1),jet_m->at(j1));
		if (max_dphi<abs(jet0_4vec.DeltaPhi(jet1_4vec))) {
			max_dphi = abs(jet0_4vec.DeltaPhi(jet1_4vec));
			jet1 = j1;
		}
	}
	max_dphi = 0;
	for (int j1=1;j1<jet_pT->size();j1++) {
		if (j1==jet1) continue;
		jet1_4vec.SetPtEtaPhiM(jet_pT->at(j1),jet_eta->at(j1),jet_phi->at(j1),jet_m->at(j1));
		if (max_dphi<abs(jet0_4vec.DeltaPhi(jet1_4vec))) {
			max_dphi = abs(jet0_4vec.DeltaPhi(jet1_4vec));
			jet2 = j1;
		}
	}
	return std::make_pair(jet1,jet2);
}
TLorentzVector RJ2CA_W4vec(double unit)
{
	TLorentzVector w_4vec;
	TLorentzVector jet0_4vec;
	TLorentzVector jet1_4vec;
	std::pair<int, int> jets;
	jets = RJ2CA_FindJetPairFromW(unit);
	jet0_4vec.SetPtEtaPhiM(jet_pT->at(jets.first),jet_eta->at(jets.first),jet_phi->at(jets.first),jet_m->at(jets.first));	
	jet1_4vec.SetPtEtaPhiM(jet_pT->at(jets.second),jet_eta->at(jets.second),jet_phi->at(jets.second),jet_m->at(jets.second));	
	w_4vec = jet0_4vec+jet1_4vec;
	return w_4vec;
}
double RJ2CA_WjetsDR(double unit)
{
	TLorentzVector jet0_4vec;
	TLorentzVector jet1_4vec;
	std::pair<int, int> jets;
	jets = RJ2CA_FindJetPairFromW(unit);
	jet0_4vec.SetPtEtaPhiM(jet_pT->at(jets.first),jet_eta->at(jets.first),jet_phi->at(jets.first),jet_m->at(jets.first));	
	jet1_4vec.SetPtEtaPhiM(jet_pT->at(jets.second),jet_eta->at(jets.second),jet_phi->at(jets.second),jet_m->at(jets.second));	
	return jet0_4vec.DeltaR(jet1_4vec);
}
TLorentzVector RJ2CA_ISR4vec(double unit)
{
	TLorentzVector isr_4vec;
	TLorentzVector jet_4vec;
	std::pair<int, int> jets;
	jets = RJ2CA_FindJetPairFromW(unit);
	for (int j1=0;j1<jet_pT->size();j1++) {
		if (j1==jets.first) continue;
		if (j1==jets.second) continue;
		jet_4vec.SetPtEtaPhiM(jet_pT->at(j1),jet_eta->at(j1),jet_phi->at(j1),jet_m->at(j1));	
		isr_4vec += jet_4vec;
	}
	return isr_4vec;
}
bool RJ2CA_Low(double unit) {
	if (lep_n!=2) return false;
	if (jet_n!=2) return false;
	if (bjet_n>0) return false;
	if (lep_pT->at(0)<25.*unit) return false;
	if (lep_pT->at(1)<25.*unit) return false;
	if (jet_pT->at(0)<30.*unit) return false;
	if (jet_pT->at(1)<30.*unit) return false;
	if (mll/unit<80) return false;
	if (mll/unit>100) return false;
	if (MET/unit<100) return false;
	TLorentzVector jet0_4vec;
	TLorentzVector jet1_4vec;
	jet0_4vec.SetPtEtaPhiM(jet_pT->at(0),jet_eta->at(0),jet_phi->at(0),jet_m->at(0));
	jet1_4vec.SetPtEtaPhiM(jet_pT->at(1),jet_eta->at(1),jet_phi->at(1),jet_m->at(1));
	double mjj = (jet0_4vec+jet1_4vec).M();
	if (mjj/unit<70) return false;
	if (mjj/unit>90) return false;
	double HT5PP = RJ2CA_HT5PP(unit);
	double H5PP = RJ2CA_H5PP(unit);
	double H2PP = RJ2CA_H2PP(unit);
	double PTPP = RJ2CA_PTPP(unit);
	double RPT_HT5PP = PTPP/(PTPP+HT5PP);
	//if (H5PP<600) return false;
	if (H5PP<400) return false;
	if (H2PP/H5PP<0.35) return false;
	if (H2PP/H5PP>0.6) return false;
	if (RPT_HT5PP>0.05) return false;
	double DPhiMETJet1st = JetEtMissDPhi(unit,0);
	double DPhiMETJet2nd = JetEtMissDPhi(unit,1);
	double MinDPhi = min(DPhiMETJet1st,DPhiMETJet2nd);
	if (MinDPhi<2.4) return false;
	return true;
}
bool RJ2CA_Compress(double unit) {
	if (lep_n!=2) return false;
	if (jet_n<3) return false;
	if (jet_n>4) return false;
	if (bjet_n>0) return false;
	if (lep_pT->at(0)<25.*unit) return false;
	if (lep_pT->at(1)<25.*unit) return false;
	if (jet_pT->at(0)<30.*unit) return false;
	if (jet_pT->at(1)<30.*unit) return false;
	if (mll/unit<80) return false;
	if (mll/unit>100) return false;
	if (MET/unit<100) return false;
	if (jet_pT->at(0)<150.*unit) return false;
	TLorentzVector met_4vec;
	met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
	TLorentzVector w_4vec = RJ2CA_W4vec(unit);
	if (w_4vec.M()/unit<50) return false;
	if (w_4vec.M()/unit>110) return false;
	TLorentzVector isr_4vec = RJ2CA_ISR4vec(unit);
	if (isr_4vec.Pt()/unit<150) return false;
	//if (isr_4vec.Pt()/unit<180) return false;
	if (abs(isr_4vec.DeltaPhi(met_4vec))<2.5) return false;
	double RISR = MET/isr_4vec.Pt();
	if (RISR<0.5) return false;
	if (RISR>1.0) return false;
	//if (RISR<0.4) return false;
	//if (RISR>1.0) return false;
	if (RJ2CA_PTCM(unit)>30) return false;
	if (abs(met_4vec.DeltaPhi(w_4vec))>1.4) return false;
	if (RJ2CA_WjetsDR(unit)>2.0) return false;
	return true;
}
bool JetFakeEtMissCut(double unit) {
	double met_sig = JetEtMissSig(unit,0);
	if (met_sig<15) return false;
	return true;
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
bool CompressWino(double unit) {
	if (jet_n<3) return false;
	if (bjet_n>0) return false;
	if (lep_pT->at(0)<25.*unit) return false;
	if (lep_pT->at(1)<25.*unit) return false;
	if (mll/unit<86) return false;
	if (mll/unit>96) return false;
	if (MET/unit<150) return false;
	double mjj_80 = JetInvMassUsing80GeVPair(unit);
	if (mjj_80<70) return false;
	if (mjj_80>100) return false;
	if (JetInvMassMin(unit)>100) return false;
	if (JetInvMassMax(unit)<70) return false;
	if (jet_pT->at(0)<150.*unit) return false;  // not a good cut for large mass splitting models (but good for small dm)
	if (PhotonEtMissDPhi(unit)<0.2) return false;
	//if (MTWinoUsing80GeVPair(unit)>150) return false;
	return true;
}
bool LowWino(double unit) {
	if (jet_n<2) return false;
	if (bjet_n>0) return false;
	if (lep_pT->at(0)<25.*unit) return false;
	if (lep_pT->at(1)<25.*unit) return false;
	if (jet_pT->at(0)<50.*unit) return false;
	//if (jet_pT->at(0)<30.*unit) return false;
	//if (jet_pT->at(1)<30.*unit) return false;
	if (mll/unit<81) return false;
	if (mll/unit>101) return false;
	if (MET/unit<150) return false;
	//if (MET/unit<125) return false;
	double mjj_80 = JetInvMassUsing80GeVPair(unit);
	if (mjj_80<70) return false;
	if (mjj_80>100) return false;
	//if (JetInvMassMin(unit)>100) return false;
	//if (JetInvMassMax(unit)<70) return false;
	double DPhiMETJet1st = JetEtMissDPhi(unit,0);
	if (DPhiMETJet1st<0.4) return false;
	if (MTWinoMaxUsing80GeVPair(unit,0)>60.) return false;
	return true;
}
bool LowPtLep(double unit) {
	if (jet_n<2) return false;
	if (lep_pT->at(0)<7.*unit) return false;
	if (lep_pT->at(1)<7.*unit) return false;
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
		inputTree->SetBranchStatus("ptsmrw_bveto"			,1);
		inputTree->SetBranchAddress("ptsmrw_bveto"		,&ptrw);
		//inputTree->SetBranchStatus("ptsmrw_ht200"			,1);
		//inputTree->SetBranchAddress("ptsmrw_ht200"		,&ptrw);
		//inputTree->SetBranchStatus("ptrw_ht400"			,1);
		//inputTree->SetBranchAddress("ptrw_ht400"		,&ptrw);
		//inputTree->SetBranchStatus("ptrw_ht1200"		,1);
		//inputTree->SetBranchAddress("ptrw_ht1200"		,&ptrw);
	}
	else {
		ptrw = 1.;
	}

	Long64_t nentries = inputTree->GetEntries();

	std::cout << "output file " << pathToOutput << std::endl;
	TFile   outputFile(pathToOutput,"recreate");

	TH1D Hist_MET = TH1D(TString("Hist_MET"),"",15,0,300);
	TH1D Hist_METl = TH1D(TString("Hist_METl"),"",30,-300,300);
	TH1D Hist_METt = TH1D(TString("Hist_METt"),"",30,-300,300);
	TH1D Hist_ZPt = TH1D(TString("Hist_ZPt"),"",30,0,1500);
	TH1D Hist_ZPt2 = TH1D(TString("Hist_ZPt2"),"",20,0,200);
	TH1D Hist_mll = TH1D(TString("Hist_mll"),"",28,21,161);
	TH1D Hist_mll2 = TH1D(TString("Hist_mll2"),"",12,61,121);
	TH1D Hist_mjj80 = TH1D(TString("Hist_mjj80"),"",20,0,200);
	TH1D Hist_mjjmin = TH1D(TString("Hist_mjjmin"),"",20,0,200);
	TH1D Hist_mjjmax = TH1D(TString("Hist_mjjmax"),"",20,0,200);
	TH1D Hist_MT2 = TH1D(TString("Hist_MT2"),"",20,0,200);
	TH1D Hist_MTWino_80 = TH1D(TString("Hist_MTWino_80"),"",16,50,850);
	TH1D Hist_MTWino_80_fine = TH1D(TString("Hist_MTWino_80_fine"),"",10,50,250);
	TH1D Hist_MTWinoMax_80 = TH1D(TString("Hist_MTWinoMax_80"),"",16,0,800);
	TH1D Hist_MTWinoMax_80_fine = TH1D(TString("Hist_MTWinoMax_80_fine"),"",10,0,200);
	TH1D Hist_RatioMETJet1st = TH1D(TString("Hist_RatioMETJet1st"),"",16,0,40);
	TH1D Hist_RatioMETJet2nd = TH1D(TString("Hist_RatioMETJet2nd"),"",16,0,40);
	TH1D Hist_RatioMETJet3rd = TH1D(TString("Hist_RatioMETJet3rd"),"",16,0,40);
	TH1D Hist_RatioMETPhoton = TH1D(TString("Hist_RatioMETPhoton"),"",20,0,5);
	TH1D Hist_RatioMETW = TH1D(TString("Hist_RatioMETW"),"",20,0,5);
	TH1D Hist_DPhiMETPhoton = TH1D(TString("Hist_DPhiMETPhoton"),"",16,0,3.2);
	TH1D Hist_DPhiMETJet1st = TH1D(TString("Hist_DPhiMETJet1st"),"",16,0,3.2);
	TH1D Hist_DPhiMETJet2nd = TH1D(TString("Hist_DPhiMETJet2nd"),"",16,0,3.2);
	TH1D Hist_njet = TH1D(TString("Hist_njet"),"",15,0,15);
	TH1D Hist_HT = TH1D(TString("Hist_HT"),"",40,0,2000);
	TH1D Hist_Lep1stPt = TH1D(TString("Hist_Lep1stPt"),"",20,0,1000);
	TH1D Hist_Lep2ndPt = TH1D(TString("Hist_Lep2ndPt"),"",20,0,1000);
	TH1D Hist_Jet1stPt = TH1D(TString("Hist_Jet1stPt"),"",40,0,2000);
	TH1D Hist_Jet2ndPt = TH1D(TString("Hist_Jet2ndPt"),"",20,0,1000);
	TH1D Hist_Jet3rdPt = TH1D(TString("Hist_Jet3rdPt"),"",20,0,200);
	TH1D Hist_JetIsrPt = TH1D(TString("Hist_JetIsrPt"),"",20,0,1000);
	TH1D Hist_RatioMETJetIsr = TH1D(TString("Hist_RatioMETJetIsr"),"",20,0,2);
	TProfile Prof_ZPtResolution = TProfile(TString("Prof_ZPtResolution"),"",7,0,1400,0,500);
	TProfile Prof_JetPtResolution = TProfile(TString("Prof_JetPtResolution"),"",7,0,1400,0,500);
	TProfile Prof_TSTResolution = TProfile(TString("Prof_TSTResolution"),"",7,0,1400,0,500);
	TProfile Prof_ZPtResolution2 = TProfile(TString("Prof_ZPtResolution2"),"",8,0,400,0,500);
	TProfile Prof_JetPtResolution2 = TProfile(TString("Prof_JetPtResolution2"),"",8,0,400,0,500);
	TProfile Prof_TSTResolution2 = TProfile(TString("Prof_TSTResolution2"),"",8,0,400,0,500);
	TH1D Hist_RJ2CA_Wmass = TH1D(TString("Hist_RJ2CA_Wmass"),"",3,50,110);
	TH1D Hist_RJ2CA_RISR = TH1D(TString("Hist_RJ2CA_RISR"),"",10,0,1);

	int step = 1;
	for (Long64_t i=0;i<nentries;i+=step) {
		if (fmod(i,1e5)==0) std::cout << pathToNtuple << " " << i << "/" << nentries << " events processed." << std::endl;
		inputTree->GetEntry(i);
		if (tag.Contains("MonoJet")) {if (!MonoJetSelection(unit)) continue;}
		else if (tag.Contains("LowMllExcess2Jet")) {if (!LowMllExcess2JetSelection(unit)) continue;}
		else if (tag.Contains("LowMllExcess3Jet")) {if (!LowMllExcess3JetSelection(unit)) continue;}
		else if (tag.Contains("LowMllExcess4Jet")) {if (!LowMllExcess4JetSelection(unit)) continue;}
		else if (tag.Contains("LowMllExcess5Jet")) {if (!LowMllExcess5JetSelection(unit)) continue;}
		else if (tag.Contains("LowMllExcessBJet")) {if (!LowMllExcessBJetSelection(unit)) continue;}
		else if (tag.Contains("CompressWino")) {if (!CompressWino(unit)) continue;}
		else if (tag.Contains("LowWino")) {if (!LowWino(unit)) continue;}
		else if (tag.Contains("LowPtLep")) {if (!LowPtLep(unit)) continue;}
		else if (tag.Contains("RJ2CA_Compress")) {if (!RJ2CA_Compress(unit)) continue;}
		else if (tag.Contains("RJ2CA_Low")) {if (!RJ2CA_Low(unit)) continue;}
		else if (tag.Contains("RJ2CA_All")) {if (!RJ2CA_Low(unit) && !RJ2CA_Compress(unit)) continue;}
		else if (tag.Contains("TwoJets")) {if (!Strong2LSelection(unit)) continue;}
		else if (tag.Contains("LowPtZ")) {if (!LowPtZSelection(unit)) continue;}
		if (tag.Contains("BadMuon")) if (!BadMuon(unit)) continue;
		if (tag.Contains("OnZ")) if (mll/unit<61.||mll/unit>121.) continue;
		if (tag.Contains("MET50")) if (MET/unit<50.) continue;
		if (tag.Contains("MET80")) if (MET/unit<80.) continue;
		if (tag.Contains("MET100")) if (MET/unit<100.) continue;
		if (tag.Contains("MET120")) if (MET/unit<120.) continue;
		if (tag.Contains("MET140")) if (MET/unit<140.) continue;
		if (tag.Contains("MET160")) if (MET/unit<160.) continue;
		if (tag.Contains("MET150")) if (MET/unit<150.) continue;
		if (tag.Contains("MET200")) if (MET/unit<200.) continue;

		double DPhiMETJet1st = JetEtMissDPhi(unit,0);
		double DPhiMETJet2nd = JetEtMissDPhi(unit,1);
		double DPhiMETPhoton = PhotonEtMissDPhi(unit);
		double RatioMETPhoton = PhotonEtMissRatio(unit);
		double RatioMETW01 = (MET/unit)/WJetPtUsing80GeVPair(unit);
		double RatioMETJet1st = JetEtMissSig(unit,0);
		double RatioMETJet2nd = JetEtMissSig(unit,1);
		double RatioMETJet3rd = JetEtMissSig(unit,2);

		TLorentzVector w_4vec = RJ2CA_W4vec(unit);
		double RJ2CA_Wmass = w_4vec.M()/unit;
		TLorentzVector isr_4vec = RJ2CA_ISR4vec(unit);
		double RISR = MET/isr_4vec.Pt();

		Hist_MET.Fill(MET/unit,totalWeight*ptrw*scale);
		Hist_METl.Fill(METl/unit,totalWeight*ptrw*scale);
		Hist_METt.Fill(METt/unit,totalWeight*ptrw*scale);
		Hist_ZPt.Fill(Z_pt/unit,totalWeight*ptrw*scale);
		Hist_ZPt2.Fill(Z_pt/unit,totalWeight*ptrw*scale);
		Hist_mll.Fill(mll/unit,totalWeight*ptrw*scale);
		Hist_mll2.Fill(mll/unit,totalWeight*ptrw*scale);
		//Hist_mjj80.Fill(JetInvMassUsing80GeVPair(unit),totalWeight*ptrw*scale);
		//Hist_mjjmin.Fill(JetInvMassMin(unit),totalWeight*ptrw*scale);
		//Hist_mjjmax.Fill(JetInvMassMax(unit),totalWeight*ptrw*scale);
		Hist_mjjmax.Fill(JetInvMass(unit,1,2),totalWeight*ptrw*scale);
		Hist_MT2.Fill(MT2W(unit),totalWeight*ptrw*scale);
		//Hist_MTWino_80.Fill(MTWinoUsing80GeVPair(unit),totalWeight*ptrw*scale);
		//Hist_MTWino_80_fine.Fill(MTWinoUsing80GeVPair(unit),totalWeight*ptrw*scale);
		//Hist_MTWinoMax_80.Fill(MTWinoMaxUsing80GeVPair(unit,0),totalWeight*ptrw*scale);
		//Hist_MTWinoMax_80_fine.Fill(MTWinoMaxUsing80GeVPair(unit,0),totalWeight*ptrw*scale);
		//Hist_RatioMETPhoton.Fill(RatioMETPhoton,totalWeight*ptrw*scale);
		//Hist_RatioMETW.Fill(RatioMETW01,totalWeight*ptrw*scale);
		Hist_DPhiMETPhoton.Fill(DPhiMETPhoton,totalWeight*ptrw*scale);
		Hist_DPhiMETJet1st.Fill(DPhiMETJet1st,totalWeight*ptrw*scale);
		Hist_DPhiMETJet2nd.Fill(DPhiMETJet2nd,totalWeight*ptrw*scale);
		Hist_njet.Fill(jet_n,totalWeight*ptrw*scale);
		Hist_HT.Fill(HT/unit,totalWeight*ptrw*scale);
		Hist_Lep1stPt.Fill(lep_pT->at(0)/unit,totalWeight*ptrw*scale);
		Hist_Lep2ndPt.Fill(lep_pT->at(1)/unit,totalWeight*ptrw*scale);
		Hist_Jet1stPt.Fill(jet_pT->at(0)/unit,totalWeight*ptrw*scale);
		if (jet_pT->size()>=2) Hist_Jet2ndPt.Fill(jet_pT->at(1)/unit,totalWeight*ptrw*scale);
		if (jet_pT->size()>=3) Hist_Jet3rdPt.Fill(jet_pT->at(2)/unit,totalWeight*ptrw*scale);
		//Hist_JetIsrPt.Fill(IsrJetPtUsing80GeVPair(unit),totalWeight*ptrw*scale);
		//Hist_RatioMETJetIsr.Fill((MET/unit)/IsrJetPtUsing80GeVPair(unit),totalWeight*ptrw*scale);
		//Hist_RatioMETJet1st.Fill(RatioMETJet1st,totalWeight*ptrw*scale);
		//Hist_RatioMETJet2nd.Fill(RatioMETJet2nd,totalWeight*ptrw*scale);
		//Hist_RatioMETJet3rd.Fill(RatioMETJet3rd,totalWeight*ptrw*scale);
		//Hist_RJ2CA_Wmass.Fill(RJ2CA_Wmass,totalWeight*ptrw*scale);
		//Hist_RJ2CA_RISR.Fill(RISR,totalWeight*ptrw*scale);

		Prof_ZPtResolution.Fill(Z_pt/unit,pow(pow(Z_pt-Z_truthPt,2),0.5)/unit,totalWeight*ptrw*scale);
		Prof_TSTResolution.Fill(jet_pT->at(0)/unit,pow(pow(MET_softTerm,2),0.5)/unit,totalWeight*ptrw*scale);
		Prof_ZPtResolution2.Fill(MET/unit,pow(pow(Z_pt-Z_truthPt,2),0.5)/unit,totalWeight*ptrw*scale);
		Prof_TSTResolution2.Fill(MET/unit,pow(pow(MET_softTerm,2),0.5)/unit,totalWeight*ptrw*scale);
		TLorentzVector jet_4vec;
		jet_4vec.SetPtEtaPhiM(jet_pT->at(0),jet_eta->at(0),jet_phi->at(0),0);
		TLorentzVector lept1_4vec;
		lept1_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
		TLorentzVector lept2_4vec;
		lept2_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);
		TLorentzVector z_4vec = lept1_4vec+lept2_4vec;
		TLorentzVector zTruth_4vec;
		zTruth_4vec.SetPtEtaPhiM(Z_truthPt,z_4vec.Eta(),z_4vec.Phi(),91.*unit);
		TLorentzVector met_4vec;
		met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
		TLorentzVector tst_4vec;
		tst_4vec.SetPtEtaPhiM(MET_softTerm,0,MET_softPhi,0);
		TLorentzVector jetTruth_4vec;
		jetTruth_4vec = met_4vec+jet_4vec+z_4vec-zTruth_4vec+tst_4vec;
		Prof_JetPtResolution.Fill(jet_pT->at(0)/unit,pow(pow((jet_4vec-jetTruth_4vec).Pt(),2),0.5)/unit,totalWeight*ptrw*scale);
		Prof_JetPtResolution2.Fill(MET/unit,pow(pow((jet_4vec-jetTruth_4vec).Pt(),2),0.5)/unit,totalWeight*ptrw*scale);
	}

	Hist_MET.Write();
	Hist_METl.Write();
	Hist_METt.Write();
	Hist_ZPt.Write();
	Hist_ZPt2.Write();
	Hist_mll.Write();
	Hist_mll2.Write();
	Hist_mjj80.Write();
	Hist_mjjmin.Write();
	Hist_mjjmax.Write();
	Hist_MT2.Write();
	Hist_MTWino_80.Write();
	Hist_MTWino_80_fine.Write();
	Hist_MTWinoMax_80.Write();
	Hist_MTWinoMax_80_fine.Write();
	Hist_RatioMETJet1st.Write();
	Hist_RatioMETJet2nd.Write();
	Hist_RatioMETJet3rd.Write();
	Hist_RatioMETPhoton.Write();
	Hist_RatioMETW.Write();
	Hist_DPhiMETPhoton.Write();
	Hist_DPhiMETJet1st.Write();
	Hist_DPhiMETJet2nd.Write();
	Hist_njet.Write();
	Hist_HT.Write();
	Hist_Lep1stPt.Write();
	Hist_Lep2ndPt.Write();
	Hist_Jet1stPt.Write();
	Hist_Jet2ndPt.Write();
	Hist_Jet3rdPt.Write();
	Hist_JetIsrPt.Write();
	Hist_RatioMETJetIsr.Write();
	Prof_TSTResolution.Write();
	Prof_ZPtResolution.Write();
	Prof_JetPtResolution.Write();
	Prof_TSTResolution2.Write();
	Prof_ZPtResolution2.Write();
	Prof_JetPtResolution2.Write();
	Hist_RJ2CA_Wmass.Write();
	Hist_RJ2CA_RISR.Write();
	std::cout << "done." << std::endl;
	outputFile.Close();
	delete inputFile;

}
void Strong2LSignal(TString pathToNtuple, TString pathToOutput, TString TreeName, bool DoReweighting, double scale, double unit) {

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
		inputTree->SetBranchStatus("ptrw_bveto"			,1);
		inputTree->SetBranchAddress("ptrw_bveto"		,&ptrw);
		//inputTree->SetBranchStatus("ptrw_ht200"			,1);
		//inputTree->SetBranchAddress("ptrw_ht200"		,&ptrw);
		//inputTree->SetBranchStatus("ptrw_ht400"			,1);
		//inputTree->SetBranchAddress("ptrw_ht400"		,&ptrw);
		//inputTree->SetBranchStatus("ptrw_ht1200"		,1);
		//inputTree->SetBranchAddress("ptrw_ht1200"		,&ptrw);
	}
	else {
		ptrw = 1.;
	}

	Long64_t nentries = inputTree->GetEntries();

	std::cout << "output file " << pathToOutput << std::endl;
	TFile   outputFile(pathToOutput,"recreate");

	TH1D Hist_MET = TH1D(TString("Hist_MET"),"",30,0,600);
	TH1D Hist_mll = TH1D(TString("Hist_mll"),"",20,4,104);
	TH1D Hist_MT2 = TH1D(TString("Hist_MT2"),"",20,0,200);
	TH1D Hist_HT = TH1D(TString("Hist_HT"),"",20,0,200);
	TH1D Hist_Jet1stPt = TH1D(TString("Hist_Jet1stPt"),"",20,0,200);
	TH1D Hist_Jet2ndPt = TH1D(TString("Hist_Jet2ndPt"),"",20,0,200);
	TH1D Hist_nlep = TH1D(TString("Hist_nlep"),"",10,0,10);
	TH1D Hist_njet = TH1D(TString("Hist_njet"),"",10,0,10);

	int step = 1;
	for (Long64_t i=0;i<nentries;i+=step) {
		if (fmod(i,1e5)==0) std::cout << pathToNtuple << " " << i << "/" << nentries << " events processed." << std::endl;
		inputTree->GetEntry(i);
		if (tag.Contains("MonoJet")) {if (!MonoJetSelection(unit)) continue;}
		else if (tag.Contains("CompressWino")) {if (!CompressWino(unit)) continue;}
		else if (tag.Contains("LowWino")) {if (!LowWino(unit)) continue;}
		else if (tag.Contains("LowPtLep")) {if (!LowPtLep(unit)) continue;}
		else {if (!Strong2LSelection(unit)) continue;}
		if (tag.Contains("BadMuon")) if (!BadMuon(unit)) continue;
		if (tag.Contains("MET50")) if (MET/unit<50.) continue;
		if (tag.Contains("MET100")) if (MET/unit<100.) continue;
		if (tag.Contains("MET150")) if (MET/unit<150.) continue;
		if (tag.Contains("MET200")) if (MET/unit<200.) continue;

		double DPhiMETJet1st = JetEtMissDPhi(unit,0);
		double DPhiMETJet2nd = JetEtMissDPhi(unit,1);

		Hist_MET.Fill(MET/unit,totalWeight*ptrw*scale);
		Hist_mll.Fill(mll/unit,totalWeight*ptrw*scale);
		Hist_MT2.Fill(MT2W(unit),totalWeight*ptrw*scale);
		Hist_HT.Fill(HT/unit,totalWeight*ptrw*scale);
		Hist_Jet1stPt.Fill(jet_pT->at(0)/unit,totalWeight*ptrw*scale);
		Hist_Jet2ndPt.Fill(jet_pT->at(1)/unit,totalWeight*ptrw*scale);
		Hist_nlep.Fill(lep_n,totalWeight*ptrw*scale);
		Hist_njet.Fill(jet_n,totalWeight*ptrw*scale);
	}

	Hist_MET.Write();
	Hist_mll.Write();
	Hist_MT2.Write();
	Hist_HT.Write();
	Hist_Jet1stPt.Write();
	Hist_Jet2ndPt.Write();
	Hist_nlep.Write();
	Hist_njet.Write();
	std::cout << "done." << std::endl;
	outputFile.Close();
	delete inputFile;

}
void RunStrong2L() {

	bool DoReweighting;
	double Scale = 36.1*1000.;
	TString pathToNtuple;
	TString pathToOutput;

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/susy/392330_ee.root");
	pathToOutput = TString("../OutputHistogram/Susy_392330_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/susy/392330_mm.root");
	pathToOutput = TString("../OutputHistogram/Susy_392330_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/susy/392354_ee.root");
	pathToOutput = TString("../OutputHistogram/Susy_392354_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/susy/392354_mm.root");
	pathToOutput = TString("../OutputHistogram/Susy_392354_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/susy/392304_ee.root");
	pathToOutput = TString("../OutputHistogram/Susy_392304_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/susy/392304_mm.root");
	pathToOutput = TString("../OutputHistogram/Susy_392304_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/wjets/wjets_ee.root");
	pathToOutput = TString("../OutputHistogram/WMC_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	
	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/wjets/wjets_mm.root");
	pathToOutput = TString("../OutputHistogram/WMC_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/tt/ttee.root");
	pathToOutput = TString("../OutputHistogram/Top_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*0.5,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/vv/vvee.root");
	pathToOutput = TString("../OutputHistogram/VV_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*0.5,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/data/data_ee.root");
	pathToOutput = TString("../OutputHistogram/Data_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,1.0,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/tt/ttmm.root");
	pathToOutput = TString("../OutputHistogram/Top_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*0.5,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/vv/vvmm.root");
	pathToOutput = TString("../OutputHistogram/VV_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*0.5,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/data/data_mm.root");
	pathToOutput = TString("../OutputHistogram/Data_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,1.0,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/gdata/gdata_ee_NoSmear.root");
	pathToOutput = TString("../OutputHistogram/GData_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = true;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,1.3,1.0);

	if (tag.Contains("McSmear"))
		pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/gdata/gdata_mm_McSmear.root");
	else if (tag.Contains("NoSmear"))
		pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/gdata/gdata_mm_NoSmear.root");
	pathToOutput = TString("../OutputHistogram/GData_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = true;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,1.3,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/zjets_221/zjets_ee.root");
	pathToOutput = TString("../OutputHistogram/ZMC221_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*0.5,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/zjets_221/zjets_mm.root");
	pathToOutput = TString("../OutputHistogram/ZMC221_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*0.5,1.0);


	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/gmc/gmc_ee_NoSmear.root");
	pathToOutput = TString("../OutputHistogram/GMC_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = true;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*1.3,1.0);
	
	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/zjets/zjets_ee.root");
	pathToOutput = TString("../OutputHistogram/ZMC_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	
	if (tag.Contains("McSmear"))
		pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/gmc/gmc_mm_McSmear.root");
	else if (tag.Contains("NoSmear"))
		pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/gmc/gmc_mm_NoSmear.root");
	pathToOutput = TString("../OutputHistogram/GMC_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = true;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*1.3,1.0);
	
	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples/zjets/zjets_mm.root");
	pathToOutput = TString("../OutputHistogram/ZMC_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	

}
void RunStrong2LSignal() {

	bool DoReweighting;
	double Scale = 36.1*1000.;
	TString pathToNtuple;
	TString pathToOutput;

	pathToNtuple = TString("../OutputNtuples/susy/374226_ee.root");
	pathToOutput = TString("../OutputHistogram/Susy_374226_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2LSignal(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	pathToNtuple = TString("../OutputNtuples/susy/374231_ee.root");
	pathToOutput = TString("../OutputHistogram/Susy_374231_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2LSignal(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	pathToNtuple = TString("../OutputNtuples/susy/374238_ee.root");
	pathToOutput = TString("../OutputHistogram/Susy_374238_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2LSignal(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	pathToNtuple = TString("../OutputNtuples/susy/374239_ee.root");
	pathToOutput = TString("../OutputHistogram/Susy_374239_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2LSignal(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	pathToNtuple = TString("../OutputNtuples/susy/374241_ee.root");
	pathToOutput = TString("../OutputHistogram/Susy_374241_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2LSignal(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	pathToNtuple = TString("../OutputNtuples/susy/372449_ee.root");
	pathToOutput = TString("../OutputHistogram/Susy_372449_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2LSignal(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	pathToNtuple = TString("../OutputNtuples/susy/372449_ee.root");
	pathToOutput = TString("../OutputHistogram/Susy_372449_ee_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2LSignal(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);

	pathToNtuple = TString("../OutputNtuples/susy/374226_mm.root");
	pathToOutput = TString("../OutputHistogram/Susy_374226_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2LSignal(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	pathToNtuple = TString("../OutputNtuples/susy/374231_mm.root");
	pathToOutput = TString("../OutputHistogram/Susy_374231_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2LSignal(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	pathToNtuple = TString("../OutputNtuples/susy/374238_mm.root");
	pathToOutput = TString("../OutputHistogram/Susy_374238_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2LSignal(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	pathToNtuple = TString("../OutputNtuples/susy/374239_mm.root");
	pathToOutput = TString("../OutputHistogram/Susy_374239_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2LSignal(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	pathToNtuple = TString("../OutputNtuples/susy/374241_mm.root");
	pathToOutput = TString("../OutputHistogram/Susy_374241_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2LSignal(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	pathToNtuple = TString("../OutputNtuples/susy/372449_mm.root");
	pathToOutput = TString("../OutputHistogram/Susy_372449_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2LSignal(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	pathToNtuple = TString("../OutputNtuples/susy/372449_mm.root");
	pathToOutput = TString("../OutputHistogram/Susy_372449_mm_Strong2L_")+tag+TString(".root");
	DoReweighting = false;
	Strong2LSignal(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale,1.0);
	
}
void RunStrong2LSS() {

	bool DoReweighting;
	double Scale = 36.1*1000.;
	TString pathToNtuple;
	TString pathToOutput;

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples_SS/data/data_ee.root");
	pathToOutput = TString("../OutputHistogram/Data_ee_Strong2LSS_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,1.0,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples_SS/data/data_mm.root");
	pathToOutput = TString("../OutputHistogram/Data_mm_Strong2LSS_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,1.0,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples_SS/wjets/wjets_ee.root");
	pathToOutput = TString("../OutputHistogram/WMC_ee_Strong2LSS_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*1.0,1.0);
	
	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples_SS/wjets/wjets_mm.root");
	pathToOutput = TString("../OutputHistogram/WMC_mm_Strong2LSS_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*1.0,1.0);
	
	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples_SS/tt/ttee.root");
	pathToOutput = TString("../OutputHistogram/Top_ee_Strong2LSS_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*1.0,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples_SS/tt/ttmm.root");
	pathToOutput = TString("../OutputHistogram/Top_mm_Strong2LSS_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*1.0,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples_SS/vv/vvee.root");
	pathToOutput = TString("../OutputHistogram/VV_ee_Strong2LSS_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*1.0,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples_SS/vv/vvmm.root");
	pathToOutput = TString("../OutputHistogram/VV_mm_Strong2LSS_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*1.0,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples_SS/zjets_221/zjets_ee.root");
	pathToOutput = TString("../OutputHistogram/ZMC221_ee_Strong2LSS_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*1.0,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples_SS/zjets_221/zjets_mm.root");
	pathToOutput = TString("../OutputHistogram/ZMC221_mm_Strong2LSS_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*1.0,1.0);


}
void RunStrong2LDF() {

	bool DoReweighting;
	double Scale = 36.1*1000.;
	TString pathToNtuple;
	TString pathToOutput;

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples_DF/tt/ttem.root");
	pathToOutput = TString("../OutputHistogram/Top_em_Strong2LDF_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*0.5,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples_DF/vv/vvem.root");
	pathToOutput = TString("../OutputHistogram/VV_em_Strong2LDF_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*0.5,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples_DF/zjets_221/zjets_em.root");
	pathToOutput = TString("../OutputHistogram/ZMC221_em_Strong2LDF_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,Scale*1.0,1.0);

	pathToNtuple = TString("/eos/atlas/user/r/rshang/OutputNtuples_DF/data/data_em.root");
	pathToOutput = TString("../OutputHistogram/Data_em_Strong2LDF_")+tag+TString(".root");
	DoReweighting = false;
	Strong2L(pathToNtuple,pathToOutput,"BaselineTree",DoReweighting,1.0,1.0);


}
void MakeSelectionHistogram_QA() {

	//tag = TString("MonoJetOnZ_PTRW_McSmear");
	//RunStrong2L();
	//tag = TString("MonoJetOnZMET50_PTRW_McSmear");
	//RunStrong2L();
	//tag = TString("MonoJetOnZMET80_PTRW_McSmear");
	//RunStrong2L();
	//tag = TString("MonoJetOnZMET100_PTRW_McSmear");
	//RunStrong2L();
	//tag = TString("MonoJetOnZ_PTRW_NoSmear");
	//RunStrong2L();

	//tag = TString("MonoJet_PTRW_NoSmear");
	//RunStrong2L();
	//tag = TString("MonoJet_PTRW_McSmear");
	//RunStrong2L();
	//tag = TString("MonoJetMET80_PTRW_McSmear");
	//RunStrong2L();
	//tag = TString("MonoJetMET100_PTRW_McSmear");
	//RunStrong2L();
	//tag = TString("MonoJetMET120_PTRW_McSmear");
	//RunStrong2L();
	//tag = TString("MonoJetMET140_PTRW_McSmear");
	//RunStrong2L();
	//tag = TString("MonoJetMET160_PTRW_McSmear");
	//RunStrong2L();
	//tag = TString("LowMllExcess2JetMET100_PTRW_McSmear");
	//RunStrong2L();
	//tag = TString("LowMllExcess3JetMET100_PTRW_McSmear");
	//RunStrong2L();
	//tag = TString("LowMllExcess4JetMET100_PTRW_McSmear");
	//RunStrong2L();
	//tag = TString("LowMllExcess5JetMET100_PTRW_McSmear");
	//RunStrong2L();
	//tag = TString("LowMllExcessBJetMET100_PTRW_McSmear");
	//RunStrong2L();

	//tag = TString("MonoJet_PTRW_McSmear");
	//RunStrong2LSS();
	//tag = TString("MonoJetMET100_PTRW_McSmear");
	//RunStrong2LSS();
	tag = TString("MonoJetMET100_PTRW_McSmear");
	RunStrong2LDF();

	//tag = TString("TwoJetsOnZ_PTRW_McSmear");
	//RunStrong2L();
	//tag = TString("TwoJetsOnZMET100_PTRW_McSmear");
	//RunStrong2L();

	//tag = TString("LowPtZOnZ_PTRW_McSmear");
	//RunStrong2L();

	//RunStrong2LSignal();

}

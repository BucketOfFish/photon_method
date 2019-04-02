//-----------------------------------------------------------------------------------------------
// this script takes the non-photon ntuples from QuickAna and generates smaller ntuples with information required by photon method
// the parameters of the function GetBaseLineEvents(string sampleID, string pathToNtuples, string ch, bool isData) are:
// 	sampleID: DSID of the MC sample
// 	pathToNtuples: the path to the input ntuple
// 	ch: which dilepton channel to use (ee,mm,em)
// 	isData: put "true" if you are running a data sample
// example of code running command: root -l -b -q 'GetBaseLineEvents.C+("361063","vv","root://eosatlas//eos/atlas/user/l/longjon/Ntuples/v00-21/MC/","ee",false)'
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
#include "TH2.h"
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

#include "BasicSetting.C"
#include "InputVariables.C"
#include "OutputVariables.C"
#include "GetDijetVariables.C"
#include "GetJigsawVariables.C"
#include "MT2_ROOT.h"
double beta_limit = 100.;
TH2D* hist_mt2max_map = new TH2D("hist_mt2max_map","",100,0,beta_limit,100,-4,4);
TH1D* hist_mt2max_beta = new TH1D("hist_mt2max_beta","",100,0,beta_limit);
#include "GetMT2Max.C"

using namespace std;

bool doTruthJetMatch = true;
bool doJigsaw = false;

void RebinHistogram(TH1D* hist) {
	double negative_yield = 0.;
	double positive_yield = 0.;
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
void GetSingleEvent(string sampleID, string outputName, string pathToNtuples, string ch, bool isData) {

	if (isData) doTruthJetMatch = false;

	
	//---------------------------------------------
	// open file, get Tree and EventCountHist
	//---------------------------------------------

	TH1::SetDefaultSumw2();

	TH1D* hist_EventCount = new TH1D("hist_EventCount","",3,0,3);
	double N_passMET100 = 0.;
	string  filename       = Form("%s%s.root",pathToNtuples.c_str(),sampleID.c_str()); 
	TFile*  inputFile      = TFile::Open(filename.c_str());
	Float_t _nGenEvents = 0.;
	if (!isData) {
		TH1D*   EventCountHist = (TH1D*) inputFile->Get("EventCountHist");
		_nGenEvents    = EventCountHist->GetBinContent(2);
		hist_EventCount->SetBinContent(1,EventCountHist->GetBinContent(1));
	}
	TTree*  inputTree              = (TTree*)inputFile->Get("outputTree");

	cout << endl;
	cout << "Opening file           : " << filename        << endl;
	cout << "Events in ntuple       : " << inputTree->GetEntries() << endl;
	if (!isData) {
		cout << "Total generated events : " << _nGenEvents     << endl;
	}

	//-----------------------------
	// access existing branches
	//-----------------------------
	
	GetBranches(inputTree, isData);

	//-----------------------------
	// add new branches
	//-----------------------------

	TFile   outputFile(TString(outputPath)+"/"+TString(outputName.c_str())+"/"+TString(sampleID.c_str())+"_"+TString(ch.c_str())+".root","recreate");
	TTree* BaselineTree;
	BaselineTree = new TTree("BaselineTree","baseline tree");
	AddBranches(BaselineTree, isData);

	Double_t MT2_max= 0;
	TBranch *b_MT2_max = BaselineTree->Branch("MT2_max",&MT2_max,"MT2_max/Double_t");
	Double_t boost_phi= 0;
	TBranch *b_boost_phi = BaselineTree->Branch("boost_phi",&boost_phi,"boost_phi/Double_t");
	Double_t boost_eta= 0;
	TBranch *b_boost_eta = BaselineTree->Branch("boost_eta",&boost_eta,"boost_eta/Double_t");
	Double_t boost_pt= 0;
	TBranch *b_boost_pt = BaselineTree->Branch("boost_pt",&boost_pt,"boost_pt/Double_t");

	Double_t MT2_truth= 0;
	TBranch *b_MT2_truth = BaselineTree->Branch("MT2_truth",&MT2_truth,"MT2_truth/Double_t");

	//-----------------------------
	// these variables do not go to output
	//-----------------------------
	double n_3jet = 0;
	double Wtruth_corr = 0;
	double Wmin_corr = 0;
	double W12_corr = 0;
	double W80_corr = 0;
	double W01_corr = 0;
	double Wjigsaw_corr = 0;

	//-----------------------------
	// loop over events
	//-----------------------------

	Long64_t nentries = inputTree->GetEntries();

	bool stop = false;
	for (Long64_t i=0;i<nentries;i+=event_interval) {

		if (stop) continue;

		if (fmod(i,1e5)==0) std::cout << i << " events processed." << std::endl;
		//std::cout << i << " events processed." << std::endl;
		inputTree->GetEntry(i);

		std::cout << MT2W << std::endl;
		//if (MT2W>70) continue;
		//if (MT2W<30) continue;
		//if (MT2W>0.1) continue;

		if (MET>100) N_passMET100 += 1; 
		if (jet_n>=3 && MET>150) n_3jet += 1;

		//if (MET<100) continue;
		//if (isData && MET>150) continue;
		
		if (lep_n<2) continue;
		if (lep_pT->at(0)<leading_lep_pt_cut) continue;
		if (lep_pT->at(1)<second_lep_pt_cut) continue;

		totalWeight = 1;
		if (!isData) totalWeight = (sampleWeight*eventWeight*lep_weight->at(0)*lep_weight->at(1)*emtrigweight*eetrigweight*mmtrigweight)/_nGenEvents;
		//if (!isData) totalWeight = (sampleWeight*eventWeight*lep_weight->at(0)*lep_weight->at(1)*trig_sf)/_nGenEvents;
		//if (!isData) totalWeight = (sampleWeight*eventWeight)/_nGenEvents;

		hist_cutflow_raw->Fill(0.);
		hist_cutflow_weight->Fill(0.,totalWeight);
		if (ch=="ee" && channel!=1) continue;
		if (ch=="mm" && channel!=0) continue;
		if (ch=="em" && (channel!=2 && channel!=3)) continue;
		hist_cutflow_raw->Fill(1.);
		hist_cutflow_weight->Fill(1.,totalWeight);
		if (!useMETtrig) {
			if (ch=="ee" && eetrig!=1) continue;
			if (ch=="mm" && mmtrig!=1) continue;
			if (ch=="em" && emtrig!=1) continue;
		}
		else {
			if (pass_trig_MET->at(0)!=1) continue;
		}
		hist_cutflow_raw->Fill(2.);
		hist_cutflow_weight->Fill(2.,totalWeight);
		if (is_OS!=1) continue;
		hist_cutflow_raw->Fill(3.);
		hist_cutflow_weight->Fill(3.,totalWeight);
		if (jet_n==0) continue;
		hist_cutflow_raw->Fill(4.);
		hist_cutflow_weight->Fill(4.,totalWeight);
		if (mll<12.) continue;
		hist_cutflow_raw->Fill(5.);
		hist_cutflow_weight->Fill(5.,totalWeight);
		//if (!isData) {if (lep_isPrompt->at(0)==0 || lep_isPrompt->at(1)==0) continue;}
		hist_cutflow_raw->Fill(6.);
		hist_cutflow_weight->Fill(6.,totalWeight);

		int njet = hist_low_njet->FindBin(jet_n)-1;
		if (jet_n>njet_bin[bin_size]) njet = bin_size-1;
		int nbjet = hist_low_nbjet->FindBin(bjet_n)-1;
		if (bjet_n>njet_bin[bin_size]) nbjet = bin_size-1;
		int pt = hist_low_pt->FindBin(Z_pt)-1;
		int smpt = hist_sm_pt->FindBin(Z_pt)-1;
		if (Z_pt>pt_bin[bin_size]) pt = bin_size-1;
		int ht = hist_low_ht->FindBin(HT)-1;
		if (HT>ht_bin[bin_size]) ht = bin_size-1;

		totalWeight = 1.;


		if (fmod(i,1e5)==0) std::cout << i << " process MC info." << std::endl;

		double Z_truthPt_dilep = 0;
		if (!isData) {

			// here we compute truth Z pT and adjust MC weights

			totalWeight = (sampleWeight*eventWeight*lep_weight->at(0)*lep_weight->at(1)*emtrigweight*eetrigweight*mmtrigweight)/_nGenEvents;
			//totalWeight = (sampleWeight*eventWeight*lep_weight->at(0)*lep_weight->at(1)*trig_sf)/_nGenEvents;
			//totalWeight = (sampleWeight*eventWeight)/_nGenEvents;
			if (totalWeight!=totalWeight) totalWeight = 0.;
			if (sampleID.find("361381")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361382")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361383")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361384")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361385")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361386")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361387")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361388")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361389")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361390")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361391")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361392")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361393")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361394")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361395")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361405")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361406")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361407")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361408")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361409")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361410")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361411")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361412")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361413")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361414")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361415")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361416")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361417")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361418")!=std::string::npos) totalWeight = totalWeight*1.1;
			if (sampleID.find("361419")!=std::string::npos) totalWeight = totalWeight*1.1;

			TLorentzVector lep0_truth4vec;
			lep0_truthPt = lep_truthPt->at(0);
			lep0_truthEta = lep_truthEta->at(0);
			lep0_truthPhi = lep_truthPhi->at(0);
			lep0_truth4vec.SetPtEtaPhiM(lep_truthPt->at(0),lep_truthEta->at(0),lep_truthPhi->at(0),0);
			TLorentzVector lep1_truth4vec;
			lep1_truthPt = lep_truthPt->at(1);
			lep1_truthEta = lep_truthEta->at(1);
			lep1_truthPhi = lep_truthPhi->at(1);
			lep1_truth4vec.SetPtEtaPhiM(lep_truthPt->at(1),lep_truthEta->at(1),lep_truthPhi->at(1),0);
			TLorentzVector dilep_truth4vec;
			dilep_truth4vec = lep0_truth4vec + lep1_truth4vec;
			Z_truthPt = dilep_truth4vec.Pt();
			Z_truthPt_dilep = dilep_truth4vec.Pt();
			Z_truthEta = dilep_truth4vec.Eta();
			Z_truthPhi = dilep_truth4vec.Phi();
			lep0_truthGammaPt = lep_truthGammaPt->at(0);
			lep1_truthGammaPt = lep_truthGammaPt->at(1);
			if (ch=="ee") Z_truthPt = truth_ZpT/1e3;  // for Sherpa ee channel, Z truth pT includes vertex photon radiation
			//if (ch=="ee") Z_truthPhi = truth_Zphi;  // for Sherpa ee channel, Z truth pT includes vertex photon radiation

		}
		totalWeight = totalWeight*event_interval;

		//---------------------------------------------
		// here we compute the MET parallel and perpendicular components
		// and DR between photon and nearby jet
		// and Oslo's MET_rel
		//---------------------------------------------
		if (fmod(i,1e5)==0) std::cout << i << " compute the MET parallel" << std::endl;

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
		for (int j=0;j<jet_pT->size();j++) {
			jet_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),jet_m->at(j));
			double DR_Lep0Jet = jet_4vec.DeltaR(lep0_4vec);
			double DR_Lep1Jet = jet_4vec.DeltaR(lep1_4vec);
			if (MinDR_Lep0Jet>DR_Lep0Jet) MinDR_Lep0Jet = DR_Lep0Jet;
			if (MinDR_Lep1Jet>DR_Lep1Jet) MinDR_Lep1Jet = DR_Lep1Jet;
			double DR_PhotonJet = jet_4vec.DeltaR(z_4vec);
			double DPhi_PhotonJet = jet_4vec.DeltaPhi(z_4vec);
			if (MinDR_PhotonJet>DR_PhotonJet) MinDR_PhotonJet = DR_PhotonJet;
			if (MinDPhi_PhotonJet>DPhi_PhotonJet) MinDPhi_PhotonJet = DPhi_PhotonJet;
		}
		double min_DPhi_MET_LepJet = 1000.;
		double DPhi_MET_LepJet = 1000.;
		for (int j=0;j<jet_pT->size();j++) {
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
		if (fmod(i,1e5)==0) std::cout << i << " compute dijet variables" << std::endl;
		z_4vec.SetPtEtaPhiM(Z_pt,Z_eta,Z_phi,0);
		met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
		GetDijetVariables(z_4vec,met_4vec);

		if (lep_truthPt->size()<2) continue;
		if (truthJet_pT->size()<2) continue;
		TLorentzVector met_4vec_truth;
		met_4vec_truth.SetPtEtaPhiM(truthMET,0,truthMET_Phi,0);
		double lep0_pt_truth = lep_truthPt->at(0);
		double lep0_eta_truth = lep_truthEta->at(0);
		double lep0_phi_truth = lep_truthPhi->at(0);
		double lep1_pt_truth = lep_truthPt->at(1);
		double lep1_eta_truth = lep_truthEta->at(1);
		double lep1_phi_truth = lep_truthPhi->at(1);
		if (lep0_pt_truth<0) {
			lep0_pt_truth = 0;
			lep0_eta_truth = 0;
			lep0_phi_truth = 0;
		}
		if (lep1_pt_truth<0) {
			lep1_pt_truth = 0;
			lep1_eta_truth = 0;
			lep1_phi_truth = 0;
		}
		double jet0_pt_truth = truthJet_pT->at(0);
		double jet0_eta_truth = truthJet_eta->at(0);
		double jet0_phi_truth = truthJet_phi->at(0);
		double jet1_pt_truth = truthJet_pT->at(1);
		double jet1_eta_truth = truthJet_eta->at(1);
		double jet1_phi_truth = truthJet_phi->at(1);
		if (jet0_pt_truth<0) {
			jet0_pt_truth = 0;
			jet0_eta_truth = 0;
			jet0_phi_truth = 0;
		}
		if (jet1_pt_truth<0) {
			jet1_pt_truth = 0;
			jet1_eta_truth = 0;
			jet1_phi_truth = 0;
		}
		TLorentzVector lep0_4vec_truth;
		lep0_4vec_truth.SetPtEtaPhiM(lep0_pt_truth,lep0_eta_truth,lep0_phi_truth,0);
		TLorentzVector lep1_4vec_truth;
		lep1_4vec_truth.SetPtEtaPhiM(lep1_pt_truth,lep1_eta_truth,lep1_phi_truth,0);
		TLorentzVector z_4vec_truth = lep0_4vec_truth+lep1_4vec_truth;
		TLorentzVector jet0_4vec_truth;
		jet0_4vec_truth.SetPtEtaPhiM(jet0_pt_truth,jet0_eta_truth,jet0_phi_truth,0);
		TLorentzVector jet1_4vec_truth;
		jet1_4vec_truth.SetPtEtaPhiM(jet1_pt_truth,jet1_eta_truth,jet1_phi_truth,0);
		TLorentzVector w_4vec_truth = jet0_4vec_truth+jet1_4vec_truth;
		MT2_truth = ComputeMT2(z_4vec_truth, w_4vec_truth, met_4vec_truth, 0, 0).Compute();
		hist_mt2max_map->Reset();
		hist_mt2max_beta->Reset();
		MT2_max = GetMT2Max(MT2_truth,z_4vec_truth, w_4vec_truth, met_4vec_truth,&boost_pt,&boost_phi);

		//MT2_max = GetMT2Max(mll,lep0_4vec, lep1_4vec, met_4vec,&boost_pt,&boost_phi);
		//std::cout << "use mll, mt2_max = " << MT2_max << std::endl;
		//MT2_max = GetMT2Max(80,lep0_4vec, lep1_4vec, met_4vec,&boost_pt,&boost_phi);
		//std::cout << "use 80, mt2_max = " << MT2_max << std::endl;
		//hist_mt2max_map->Reset();
		//hist_mt2max_beta->Reset();
		//MT2_max = GetMT2Max(MT2W,lep0_4vec, lep1_4vec, met_4vec,&boost_pt,&boost_phi);
		//std::cout << "use MT2W, mt2_max = " << MT2_max << std::endl;
		//MT2_max = GetMT2Max(MT2_max,lep0_4vec, lep1_4vec, met_4vec,&boost_pt,&boost_phi);
		//std::cout << "use MT2_max, mt2_max = " << MT2_max << std::endl;
		//MT2_max = GetMT2Max(MT2_max,lep0_4vec, lep1_4vec, met_4vec,&boost_pt,&boost_phi);
		//std::cout << "use MT2_max, mt2_max = " << MT2_max << std::endl;
		//MT2_max = GetMT2Max(MT2_max,lep0_4vec, lep1_4vec, met_4vec,&boost_pt,&boost_phi);
		//std::cout << "use MT2_max, mt2_max = " << MT2_max << std::endl;
		//MT2_max = GetMT2Max(MT2_max,lep0_4vec, lep1_4vec, met_4vec,&boost_pt,&boost_phi);
		//std::cout << "use MT2_max, mt2_max = " << MT2_max << std::endl;

		//---------------------------------------------
		// prompt jet truth matching
		//---------------------------------------------
		if (fmod(i,1e5)==0) std::cout << i << " prompt jet truth matching" << std::endl;
		TLorentzVector truthJet_4vec;
		TLorentzVector recoJet_4vec;
		jet_isPrompt->clear();
		for (int j1=0;j1<jet_pT->size();j1++) {
			jet_isPrompt->push_back(0);
		}
		if (doTruthJetMatch) {
			if (!isData) {
				for (int j0=0;j0<truthJet_pT->size();j0++) {
					truthJet_4vec.SetPtEtaPhiM(truthJet_pT->at(j0),truthJet_eta->at(j0),truthJet_phi->at(j0),0);
					double min_DR = 1000.;
					int reco_this = 0;
					for (int j1=0;j1<jet_pT->size();j1++) {
						recoJet_4vec.SetPtEtaPhiM(jet_pT->at(j1),jet_eta->at(j1),jet_phi->at(j1),0);
						double DR_this = truthJet_4vec.DeltaR(recoJet_4vec);
						if (DR_this<min_DR) {
							min_DR = DR_this;
							reco_this = j1;
						}
					}
					if (min_DR<0.3) jet_isPrompt->at(reco_this) = 1;
				}
			}
		}
		int njet_turth_prompt = 0;
		for (int j=0;j<jet_isPrompt->size();j++) {
			if (jet_isPrompt->at(j)==1) njet_turth_prompt += 1;
		}

		//---------------------------------------------
		// apply Jigsaw rules to find ISR jets, without knowing the flavor of the objects
		// this is for the compressed region
		//---------------------------------------------
		if (fmod(i,1e5)==0) std::cout << i << " compute jigsaw variables" << std::endl;
		if (doJigsaw) {
			GetJigsawVariables();
		}

		//---------------------------------------------
		// Fill METl histograms for smearing, and mll histograms for mll modeling
		// hist_dPt_Pt histogram, i.e. the Z truth-reco response function
		// as well as the hist_2LPt_Pt histogram, which translates Z pT to dilepton sum pT
		//---------------------------------------------
		if (fmod(i,1e5)==0) std::cout << i << " Fill METl histograms" << std::endl;
		double pt37_cut = 37.;
		if (!isData) {
			if (Z_pt>pt37_cut && abs(MinDPhi_PhotonJet)>0.0) {
				truthMETt = truthMET*TMath::Sin(truthMET_Phi-Z_truthPhi);
				truthMETl = truthMET*TMath::Cos(truthMET_Phi-Z_truthPhi);
				int pt_truth = hist_low_pt->FindBin(Z_truthPt)-1;
				int dpt = hist_low_dpt->FindBin((Z_pt-Z_truthPt))-1;
				int dpt_dilep = hist_low_dpt->FindBin((Z_pt-Z_truthPt_dilep))-1;
				if (Z_truthPt>pt_bin[bin_size]) pt_truth = bin_size-1;
				if (pt>=0) hist_METl_Pt[pt]->Fill(METl,totalWeight);
				if (pt>=0) hist_METt_Pt[pt]->Fill(METt,totalWeight);
			}
		}

		stop = true;

		//BaselineTree.Fill();     
		BaselineTree->Fill();     

	}

	std::cout << "write output..." << std::endl;
	//T->Print();
	//BaselineTree.Write();
	BaselineTree->Write();

	hist_mt2max_map->Write();
	hist_mt2max_beta->Write();
	hist_cutflow_raw->Write();
	hist_cutflow_weight->Write();
	for (int bin=0;bin<bin_size;bin++) {
		hist_METl_Pt[bin]->Write();
		hist_METt_Pt[bin]->Write();
	}
	if (!isData) {
		hist_EventCount->SetBinContent(2,nentries);
		hist_EventCount->SetBinContent(3,N_passMET100);
		hist_EventCount->Write();
	}

	std::cout << "done." << std::endl;
	outputFile.Close();
	delete inputFile;


}

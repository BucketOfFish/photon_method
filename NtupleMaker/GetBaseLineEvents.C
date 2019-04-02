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
#include "InputVariables.C"
#include "OutputVariables.C"
#include "../PhotonSmearing/GetDijetVariables.C"
float beta_limit = 10.;

using namespace std;

bool doTruthJetMatch = false;
bool doJigsaw = false;

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
void GetBaseLineEvents(string sampleID, string outputName, string pathToNtuples, bool isData, string treename = "outputTree" ) {

        //if (isData) doTruthJetMatch = false;
	doTruthJetMatch = false;
	
	//---------------------------------------------
	// open file, get Tree and EventCountHist
	//---------------------------------------------

	TH1::SetDefaultSumw2();

	TH1D* hist_EventCount = new TH1D("hist_EventCount","",3,0,3);
	float N_passMET100 = 0.;
	string  filename       = Form("%s%s.root",pathToNtuples.c_str(),sampleID.c_str()); 
	TFile*  inputFile      = TFile::Open(filename.c_str());
	Float_t _nGenEvents = 1.;
	if (!isData) {
	        //TH1D*   EventCountHist = (TH1D*) inputFile->Get("EventCountHist");
		//_nGenEvents    = EventCountHist->GetBinContent(2);
		//hist_EventCount->SetBinContent(1,EventCountHist->GetBinContent(1));

		cout << "Setting _nGenEvents = 1 for now NEED TO FIX" << endl;
		_nGenEvents    = 1.0;
		hist_EventCount->SetBinContent(1,1.0);
}
	TTree*  inputTree              = (TTree*)inputFile->Get( treename.c_str() );

	std::cout << inputTree << std::endl;
	cout << endl;
	cout << "Opening file           : " << filename        << endl;
	cout << "Events in ntuple       : " << inputTree->GetEntries() << endl;
	if (!isData) {
		cout << "Total generated events : " << _nGenEvents     << endl;
	}
	
	//-----------------------------
	// access existing branches
	//-----------------------------
        //std::cout << "testing1" << std::endl; 
	GetBranches(inputTree, isData);
	//std::cout << "testing2" << std::endl;
	//-- see code for branches, cut and pasted below
	
	//-----------------------------
	// add new branches
	//-----------------------------

	string outfilename = ntuple_path + "/" + outputName + "/" + sampleID.c_str() + ".root";
	cout << "Writing to : " << outfilename << endl;
	
	TFile   outputFile( outfilename.c_str() , "recreate" );

	TTree* BaselineTree;
	BaselineTree = new TTree("BaselineTree","baseline tree");
	AddBranches(BaselineTree, isData);

	Float_t MT2_max= 0;
	TBranch *b_MT2_max = BaselineTree->Branch("MT2_max",&MT2_max,"MT2_max/F");
	Float_t boost_phi= 0;
	TBranch *b_boost_phi = BaselineTree->Branch("boost_phi",&boost_phi,"boost_phi/F");
	Float_t boost_eta= 0;
	TBranch *b_boost_eta = BaselineTree->Branch("boost_eta",&boost_eta,"boost_eta/F");
	Float_t boost_pt= 0;
	TBranch *b_boost_pt = BaselineTree->Branch("boost_pt",&boost_pt,"boost_pt/F");
	//Float_t MT2_truth= 0;
	//TBranch *b_MT2_truth = BaselineTree->Branch("MT2_truth",&MT2_truth,"MT2_truth/Float_t");
	
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
		//std::cout << i << " events processed." << std::endl;
		inputTree->GetEntry(i);

		if (MET>100) N_passMET100 += 1; 
		if (jet_n>=3 && MET>150) n_3jet += 1;

		//if (MET<100) continue;
		//if (isData && MET>150) continue;
		
		if ( nLep_signal  != 2                 ) continue; // exactly 2 signal leptons
		if ( nLep_base    != 2                 ) continue; // exactly 2 baseline leptons
		if ( lep_pT->at(0) < leading_lep_pt_cut ) continue; // 1st lep pT > 25 GeV
		if ( lep_pT->at(1) < second_lep_pt_cut  ) continue; // 2nd lep pT > 25 GeV

		//--- evaluate weight
		totalWeight = 1;
		//if (!isData) totalWeight = (sampleWeight*eventWeight*lep_weight->at(0)*lep_weight->at(1)*emtrigweight*eetrigweight*mmtrigweight)/_nGenEvents;
		//if (!isData) totalWeight = (sampleWeight*eventWeight*lep_weight->at(0)*lep_weight->at(1)*trig_sf)/_nGenEvents;
		//if (!isData) totalWeight = (sampleWeight*eventWeight)/_nGenEvents;
		if (!isData) totalWeight = genWeight * eventWeight * leptonWeight * jvtWeight * bTagWeight * pileupWeight * FFWeight;

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

		//hist_cutflow_raw->Fill(0.);
		//hist_cutflow_weight->Fill(0.,totalWeight);

		if( channel < 0 ) continue; // require exactly 2 signal leptons and corresponding triggers
		if( is_OS != 1  ) continue; // require opposite-sign
		if( jet_n < 1   ) continue; // require at least 1 pT > 30 GeV jets
		
		//if (ch=="ee" && channel!=1) continue;
		//if (ch=="mm" && channel!=0) continue;
		//if (ch=="em" && (channel!=2 && channel!=3)) continue;

		/*
		//--- commenting out for now
		  
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
		//if (is_OS==1) continue;
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
		*/
		
		int njet = hist_low_njet->FindBin(jet_n)-1;
		if (jet_n>njet_bin[bin_size]) njet = bin_size-1;
		int nbjet = hist_low_nbjet->FindBin(bjet_n)-1;
		if (bjet_n>njet_bin[bin_size]) nbjet = bin_size-1;
		int pt = hist_low_pt->FindBin(Z_pt)-1;
		int smpt = hist_sm_pt->FindBin(Z_pt)-1;
		if (Z_pt>pt_bin[bin_size]) pt = bin_size-1;
		int ht = hist_low_ht->FindBin(HT)-1;
		if (HT>ht_bin[bin_size]) ht = bin_size-1;

		/*
		totalWeight = 1.;

		if (fmod(i,1e5)==0) std::cout << i << " process MC info." << std::endl;

		float Z_truthPt_dilep = 0;
		if (!isData) {

			// here we compute truth Z pT and adjust MC weights

			//totalWeight = (sampleWeight*eventWeight*lep_weight->at(0)*lep_weight->at(1)*emtrigweight*eetrigweight*mmtrigweight)/_nGenEvents;
			//totalWeight = (sampleWeight*eventWeight*lep_weight->at(0)*lep_weight->at(1)*trig_sf)/_nGenEvents;
			totalWeight = (sampleWeight*eventWeight)/_nGenEvents;
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
			Z_truthPt       = dilep_truth4vec.Pt();
			Z_truthPt_dilep = dilep_truth4vec.Pt();
			Z_truthEta = dilep_truth4vec.Eta();
			Z_truthPhi = dilep_truth4vec.Phi();
			lep0_truthGammaPt = lep_truthGammaPt->at(0);
			lep1_truthGammaPt = lep_truthGammaPt->at(1);
			if (ch=="ee") Z_truthPt = truth_ZpT/1e3;  // for Sherpa ee channel, Z truth pT includes vertex photon radiation
			//if (ch=="ee") Z_truthPhi = truth_Zphi;  // for Sherpa ee channel, Z truth pT includes vertex photon radiation

		}
		totalWeight = totalWeight*event_interval;
		*/
		
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
		
		//if (fmod(i,1e5)==0) std::cout << i << " compute the MET parallel" << std::endl;

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
	if (!isData) {
		hist_EventCount->SetBinContent(2,nentries);
		hist_EventCount->SetBinContent(3,N_passMET100);
		hist_EventCount->Write();
	}

	std::cout << "done." << std::endl;
	outputFile.Close();
	delete inputFile;


}

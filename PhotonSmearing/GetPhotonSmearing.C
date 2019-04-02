//-----------------------------------------------------------------------------------------------
// this script takes the outputs from GetBaseLineEvents.C and GetPhotonEvents.C, and makes smeared photon information.
// the parameters of the function GetPhotonSmearing(string label, string ch, int isData, string period, int smearing_method) are:
//  label: filename
// 	ch: which dilepton channel (ee,mm)
// 	isData: 0 (MC) or 1 (data) 
//  period: data15-16 (input) -> ZMC16a (source file), data17 -> ZMC16cd, data18 -> ZMC16cd
//  smearing_method: 0 (no smear), 1 (MC smear), 2 (data smear), 3 (truth smear), 4 (MC smear)
// example: root -l -b 'GetPhotonSmearing.C+("ee-channel","ee",0,"data15-16",0)'  //for no smear
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
#include "TSpectrum.h"
#include "TVirtualFFT.h"

#include "../BasicSetting.C"
#include "PhotonVariables.C"
#include "GetDijetVariables.C"
#include "GetSmearingHistogram.C"
#include "GetMllHistogram.C"
#include "GetIndividualLeptonInfo.C"
#include "MT2_ROOT.h"

using namespace std;

vector<string> noSampleWeight;
vector<string> noEventWeight;


int RebinHistogram(TH1D* hist, int rebin) {
	float negative_yield = 0.;
	float positive_yield = 0.;
	float core_yield = 0.;
	int remainder = 0;
	if (rebin==0) {
		rebin = 1;
		for (int bin=1;bin<=hist->GetNbinsX();bin++) {
			if (hist->GetBinContent(bin)>=0) positive_yield += hist->GetBinContent(bin);
			else negative_yield += hist->GetBinContent(bin);
		}
		remainder = hist->GetNbinsX() % 2;
		while ((abs(negative_yield/positive_yield)>0.005 || core_yield/positive_yield<0.4) && remainder==0 && rebin<=32) {
			hist->Rebin(2);
			rebin = rebin*2;
			remainder = hist->GetNbinsX() % 2;
			negative_yield = 0.;
			positive_yield = 0.;
			core_yield = 0.;
			for (int bin=1;bin<=hist->GetNbinsX();bin++) {
				if (hist->GetBinContent(bin)>=0) positive_yield += hist->GetBinContent(bin);
				else negative_yield += hist->GetBinContent(bin);
				if (abs(hist->GetBinCenter(bin)-hist->GetMean())<hist->GetRMS()) {
					core_yield += hist->GetBinContent(bin); // core_yield = 68% for a perfect Guassian
				}
			}
		}
	}
	else {
		hist->Rebin(rebin);
	}
	//std::cout  << "negative_yield = " << negative_yield << ", positive_yield = " << positive_yield << ", remainder = " << remainder << std::endl;
	for (int bin=1;bin<=hist->GetNbinsX();bin++) {
		hist->SetBinContent(bin,max(hist->GetBinContent(bin),0.));
	}
	return rebin;
}
void getPhotonSmearingFunction(TString file, TString histname, TH1D* hist, int rebin) {
	TFile *fData = TFile::Open( file );
	fData->cd("");
	TH1D* temp = (TH1D*)fData->Get(histname);
	for (int bin=1;bin<=temp->GetNbinsX();bin++) {
		hist->SetBinContent(bin,temp->GetBinContent(bin));
	}
	fData->Close();
}
void getPhotonSmearingFunction2(TString file, TString histname, TH1D* hist, int rebin) {
	TFile *fData = TFile::Open( file );
	fData->cd("");
	TH1D* temp = (TH1D*)fData->Get(histname);
	//temp->Rebin(rebin);
	//hist->Rebin(rebin);
	for (int bin=1;bin<=temp->GetNbinsX();bin++) {
		hist->SetBinContent(bin,max(temp->GetBinContent(bin),0.));
	}
	fData->Close();
}
void GetPhotonSmearing(string label, string ch, int isData, string period, int smearing_method) {

        cout << "channel         " << ch              << endl;
	cout << "period          " << period          << endl;
	cout << "isData?         " << isData          << endl;
	cout << "smearing path   " << smearingPath    << endl;
	cout << "smearing method " << smearing_method << endl;

	//-----------------------------
	// prepare lep pT functions, to convert photon pT to dilepton sum pT
	//-----------------------------
	
	std::cout << "Prepare Mll histograms..." << std::endl;

	GetMllHistogram(ch,period);
	//for (int bin=0;bin<dpt_bin_size;bin++) {
	//	int rebin = RebinHistogram(hist_Mll_dPt[bin],0);
	//}
	for (int bin0=0;bin0<bin_size;bin0++) {
		for (int bin1=0;bin1<dpt_bin_size;bin1++) {
			int rebin = RebinHistogram(hist_Mll_dPt[bin0][bin1],0);
		}
	}

	//-----------------------------
	// prepare smearing functions
	//-----------------------------
	
	std::cout << "Prepare smearing histograms..." << std::endl;

	if (smearing_method == 0) photon_tag = "_NoSmear";
	if (smearing_method == 1) photon_tag = "_McSmear";
	if (smearing_method == 2) photon_tag = "_DataSmear";
	if (smearing_method == 3) photon_tag = "_TruthSmear";
	if (smearing_method == 4) photon_tag = "_McSmear";
	if (smearing_method == 5) photon_tag = "_DataSmear";

	cout << "smearing_method    " << smearing_method << endl;
	cout << "photon_tag         " << photon_tag      << endl;
	
	TSpectrum pfinder;
	TH1D* g_resp[bin_size];
	TH1D* z_resp[bin_size];
	TH1D* smear_raw[bin_size];
	TH1D* smear_fft_re[bin_size];
	TH1D* smear_fft_im[bin_size];
	TH1D* smear_fft_amp[bin_size];
	TH1D* smear_final[bin_size];
	TH1D* smear_final_phi[bin_size];
	float shift[bin_size];
	TH1D* g_metl_smear[bin_size];
	TH1D* g_metl_smear_2j[bin_size];
	for (int bin=0;bin<bin_size;bin++) {
		g_metl_smear[bin] = new TH1D(TString("g_metl_smear_")+TString::Itoa(bin,10),"",40000,-30000,10000);
		g_metl_smear[bin]->SetStats(0);
		g_metl_smear_2j[bin] = new TH1D(TString("g_metl_smear_2j_")+TString::Itoa(bin,10),"",40000,-30000,10000);
		g_metl_smear_2j[bin]->SetStats(0);
	}


	GetSmearingHistogram(ch, lumi, photon_tag,period,smearing_method);

	if (smearing_method >= 0) {  // if you want to use the deconvolution method to smear photon events. To enable this method, set "bool useDeconvolution = true" in BasicSetting.C
		for (int bin=0;bin<bin_size;bin++) {
			int rebin = 1;
			rebin = RebinHistogram(z_metl[bin],0);
			rebin = RebinHistogram(z_metl_2j[bin],rebin);
			rebin = RebinHistogram(g_metl[bin],rebin);
			rebin = RebinHistogram(g_metl_smear[bin],rebin);
			rebin = RebinHistogram(g_metl_smear_2j[bin],rebin);
			float gmetl_mean = g_metl[bin]->GetMean();
			float gmetl_rms = g_metl[bin]->GetRMS();
			float zmetl_mean = z_metl[bin]->GetMean();
			float zmetl_rms = z_metl[bin]->GetRMS();
			//std::cout << "--------------------------------------------------------" << std::endl;
			//std::cout << "bin " << bin << ", rebin " << rebin << std::endl;
			//std::cout << "gmetl_mean " << gmetl_mean << ", gmetl_rms " << gmetl_rms << std::endl;
			//std::cout << "zmetl_mean " << zmetl_mean << ", zmetl_rms " << zmetl_rms << std::endl;
			int newbin = 40000/rebin;
			Float_t *smear = new Float_t[2*((newbin+1)/2+1)];
			Float_t *fft_re = new Float_t[newbin];
			Float_t *fft_im = new Float_t[newbin];
			Double_t *z_smear_in = new Double_t[newbin];
			Double_t g_smear_in[newbin];
			Double_t j_resp_in[newbin];
		        Double_t *z_resp_in = new Double_t[newbin];
			Double_t *g_resp_in = new Double_t[newbin];
			g_resp[bin] = new TH1D(TString("g_resp_")+TString::Itoa(bin,10),"",newbin,-30000,10000);
			z_resp[bin] = new TH1D(TString("z_resp_")+TString::Itoa(bin,10),"",newbin,-30000,10000);
			smear_raw[bin] = new TH1D(TString("smear_raw_")+TString::Itoa(bin,10),"",newbin,-30000,10000);
			smear_fft_re[bin] = new TH1D(TString("smear_fft_re_")+TString::Itoa(bin,10),"",newbin,-30000,10000);
			smear_fft_im[bin] = new TH1D(TString("smear_fft_im_")+TString::Itoa(bin,10),"",newbin,-30000,10000);
			smear_fft_amp[bin] = new TH1D(TString("smear_fft_amp_")+TString::Itoa(bin,10),"",newbin,-30000,10000);
			for (int i=0;i<newbin;i++) {
				if (smearing_method == 3) {
				  cout << "smearing_method 3 not yet implemented, quitting!" << endl;
				  exit(0);
					// z_smear_in[i] = max(z_dpt[bin]->GetBinContent(i+1),0.);
					// if (i<newbin/2) g_smear_in[i] = max(g_dpt[bin]->GetBinContent(i+1+newbin/2),0.);
					// else g_smear_in[i] = 0.;
				}
				else {
					z_smear_in[i] = max(z_metl[bin]->GetBinContent(i+1),0.);
					if (i<newbin/2) g_smear_in[i] = max(g_metl[bin]->GetBinContent(i+1+newbin/2),0.);
					else g_smear_in[i] = 0.;
				}
				z_resp_in[i] = max(z_metl[bin]->GetBinContent(i+1),0.);
				g_resp_in[i] = max(g_metl[bin]->GetBinContent(i+1),0.);
				if (i<newbin/2) j_resp_in[i] = max(z_jetmetl[bin]->GetBinContent(i+1+newbin/2),0.);
				else j_resp_in[i] = 0.;
			}
			pfinder.Deconvolution(z_smear_in,g_smear_in,newbin,1000,1,1.0);
			pfinder.Deconvolution(z_resp_in,j_resp_in,newbin,1000,1,1.0);
			pfinder.Deconvolution(g_resp_in,j_resp_in,newbin,1000,1,1.0);
			for (int i=0;i<newbin;i++) {
				smear_raw[bin]->SetBinContent(i+1,z_smear_in[i]);
				g_resp[bin]->SetBinContent(i+1,g_resp_in[i]);
				z_resp[bin]->SetBinContent(i+1,z_resp_in[i]);
				smear[i] = z_smear_in[i];
			}
			// rebin = RebinHistogram(z_dpt[bin],0);
			// rebin = RebinHistogram(g_dpt[bin],rebin);
			// float gdpt_mean = g_dpt[bin]->GetMean();
			// float gdpt_rms = g_dpt[bin]->GetRMS();
			// float zdpt_mean = z_dpt[bin]->GetMean();
			// float zdpt_rms = z_dpt[bin]->GetRMS();
			float smear_mean = smear_raw[bin]->GetMean();
			float smear_rms = smear_raw[bin]->GetRMS();
			//std::cout << "gdpt_mean " << gdpt_mean << ", gdpt_rms " << gdpt_rms << std::endl;
			//std::cout << "zdpt_mean " << zdpt_mean << ", zdpt_rms " << zdpt_rms << std::endl;

			//std::cout << "smear_mean " << smear_mean << ", smear_rms " << smear_rms << std::endl;
			for (int i=0;i<newbin;i++) {
				if (gmetl_rms/zmetl_rms > 1.0) {
					smear_raw[bin]->SetBinContent(i+1,0.);
					smear[i] = 0.;
				}
				float smear_cut = 6.;
				if (ch=="mm" && sm_pt_bin[bin]>=0) smear_cut = 7.;
				if (abs(smear_raw[bin]->GetBinCenter(i+1)-smear_mean)/smear_rms>smear_cut) {
					smear_raw[bin]->SetBinContent(i+1,0.);
					smear[i] = 0.;
				}
			}
			// smearing_method 3 isn't implemented
			//if (smearing_method == 3) shift[bin] = -g_dpt[bin]->GetMean();
			//else shift[bin] = -g_metl[bin]->GetMean();

			shift[bin] = -g_metl[bin]->GetMean();
		}
	}
	for (int bin=0;bin<bin_size;bin++) {
		smear_final[bin] = new TH1D(TString("smear_final_")+TString::Itoa(bin,10),"",500,-1000,1000);
 		for (int i=0;i<500;i++) {
			int which_bin = smear_raw[bin]->FindBin(smear_final[bin]->GetBinCenter(i+1));
			smear_final[bin]->SetBinContent(i+1,smear_raw[bin]->GetBinContent(which_bin));
	 	}
	}


	//---------------------------------------------
	// open file, get Tree and EventCountHist
	//---------------------------------------------

	TH1::SetDefaultSumw2();

	string  infilename;
	//if (isData==1) infilename = TString(TString(smearingPath)+"gdata/gdata_raw.root"); 
	//if (isData==0) infilename = TString(TString(smearingPath)+"gmc/gmc_raw.root"); 

	if (isData==1){
	  infilename = outputPath + "gdata/" + label + ".root";
	  cout << __FILE__ << " " << __LINE__ << endl;
	}
	if (isData==0){
	  infilename = outputPath + "gmc/" + label + ".root"; //"../OutputNtuples_v1.2/gmc/SinglePhoton211_merged_processed.root";
	  cout << __FILE__ << " " << __LINE__ << endl;
	}

	TChain* T = new TChain("BaselineTree");

	T->Add( infilename.c_str() );
	
	cout << endl;
	cout << "Opening file           : " << infilename        << endl;
	cout << "Events in ntuple       : " << T->GetEntries() << endl;

	//---------------------------------------------
	// open file, get Tree and EventCountHist
	//---------------------------------------------

	TH1::SetDefaultSumw2();

	string outfilename;
	if (isData==1) outfilename = TString(TString(outputPath)+"gdata/" + label + "_"+TString(ch)+TString(photon_tag)+".root"); 
	if (isData==0){
	  outfilename = TString(TString(outputPath)+"gmc/gmc_"+TString(ch)+TString(photon_tag)+".root"); 

	}
	TFile*  f              = new TFile(outfilename.c_str(),"recreate");          
	TTree* BaselineTree = new TTree("BaselineTree","baseline tree");

	cout << endl;
	cout << "Create file           : " << outfilename        << endl;


	//-----------------------------
	// access existing branches
	//-----------------------------
	
	T->SetBranchAddress("totalWeight"     ,&totalWeight              );
	T->SetBranchAddress("pt"              ,&pt               );
	T->SetBranchAddress("ht"              ,&ht               );
	T->SetBranchAddress("njet"            ,&njet             );
	T->SetBranchAddress("nbjet"            ,&nbjet             );
	T->SetBranchAddress("HT"              ,&HT               );
	T->SetBranchAddress("jet_n"           ,&jet_n            );
	T->SetBranchAddress("bjet_n"           ,&bjet_n            );
	// if (isData!=1) T->SetBranchAddress("truthGamma_pt", &truthGamma_pt);
	// if (isData!=1) T->SetBranchAddress("truthGamma_phi", &truthGamma_phi);
	// if (isData!=1) T->SetBranchAddress("truthGamma_eta", &truthGamma_eta);
	T->SetBranchAddress("gamma_pt", &gamma_pt);
	T->SetBranchAddress("gamma_eta", &gamma_eta);
	T->SetBranchAddress("gamma_phi", &gamma_phi);
	if (isData!=1) T->SetBranchAddress("gamma_dR", &gamma_dR);
	T->SetBranchAddress("MET_raw"             ,&MET              );
        //------------------------------------------------------
        // 2019 RJR analysis variables -------------------------
        //------------------------------------------------------
        T->SetBranchAddress("MET_loose"             ,&MET_loose              );
        T->SetBranchAddress("MET_tight"             ,&MET_tight              );
        T->SetBranchAddress("MET_tighter"             ,&MET_tighter              );
        T->SetBranchAddress("MET_tenacious"             ,&MET_tenacious              );
        T->SetBranchAddress("trigMatch_2LTrigOR"         ,&trigMatch_2LTrigOR          );
        T->SetBranchAddress("is2Lep2Jet"         ,&is2Lep2Jet         );
        T->SetBranchAddress("is2L2JInt"         ,&is2L2JInt         );
        T->SetBranchAddress("nBJet20_MV2c10_FixedCutBEff_77"    ,&nBJet20_MV2c10_FixedCutBEff_77 );
        T->SetBranchAddress("mjj"         ,&mjj          );
        T->SetBranchAddress("mll_RJ"          ,&mll_RJ          );
        T->SetBranchAddress("R_minH2P_minH3P"          ,&R_minH2P_minH3P           );
        T->SetBranchAddress("RPT_HT5PP"              ,&RPT_HT5PP               );
        T->SetBranchAddress("dphiVP"          ,&dphiVP            );
        T->SetBranchAddress("H2PP"           ,&H2PP           );
        T->SetBranchAddress("H5PP"           ,&H5PP           );
        T->SetBranchAddress("nJet20"             ,&nJet20            );
        T->SetBranchAddress("minDphi"          ,&minDphi           );
        T->SetBranchAddress("MZ"         ,&MZ          );
        T->SetBranchAddress("NjS"         ,&NjS        );
        T->SetBranchAddress("NjISR"           ,&NjISR           );
        T->SetBranchAddress("dphiISRI"           ,&dphiISRI           );
        T->SetBranchAddress("RISR"          ,&RISR           );
        T->SetBranchAddress("PTISR"          ,&PTISR           );
        T->SetBranchAddress("PTI"         ,&PTI          );
        T->SetBranchAddress("PTCM"         ,&PTCM          );
        T->SetBranchAddress("MJ",            &MJ   );
        T->SetBranchAddress("is3Lep3Jet",           &is3Lep3Jet  );
        T->SetBranchAddress("is4Lep3Jet",           &is4Lep3Jet  );
        T->SetBranchAddress("lept1sign_VR",           &lept1sign_VR  );
        T->SetBranchAddress("lept2sign_VR",           &lept2sign_VR   );
        T->SetBranchAddress("lept1Pt_VR",           &lept1Pt_VR );
        T->SetBranchAddress("lept2Pt_VR",            &lept2Pt_VR );
        T->SetBranchAddress("MZ_VR",                  &MZ_VR );
        T->SetBranchAddress("MJ_VR",                 &MJ_VR );
        T->SetBranchAddress("RISR_VR",           &RISR_VR  );
        T->SetBranchAddress("PTISR_VR",           &PTISR_VR  );
        T->SetBranchAddress("PTI_VR",           &PTI_VR  );
        T->SetBranchAddress("PTCM_VR",           &PTCM_VR   );
        T->SetBranchAddress("dphiISRI_VR",           &dphiISRI_VR );
        T->SetBranchAddress("lepFlavor",                  &lepFlavor );
        T->SetBranchAddress("lepCharge",                 &lepCharge );
        //------------------------------------------------------
	T->SetBranchAddress("METl_raw"             ,&METl              );
	T->SetBranchAddress("METt_raw"             ,&METt              );
	T->SetBranchAddress("MET_phi_raw"         ,&MET_phi          );
	T->SetBranchAddress("DPhi_METJetLeading_raw"         ,&DPhi_METJetLeading          );
	T->SetBranchAddress("DPhi_METJetSecond_raw"          ,&DPhi_METJetSecond           );
	T->SetBranchAddress("MinDPhi_PhotonJet"          ,&MinDPhi_PhotonJet           );
	T->SetBranchAddress("lep_pT_raw"          ,&lep_pT           );
	T->SetBranchAddress("jet_pT"          ,&jet_pT           );
	T->SetBranchAddress("jet_eta"         ,&jet_eta          );
	T->SetBranchAddress("jet_phi"         ,&jet_phi          );
	T->SetBranchAddress("jet_m"         ,&jet_m          );

	//-----------------------------
	// add new branches
	//-----------------------------
	
	TBranch *b_totalWeight = BaselineTree->Branch("totalWeight",&totalWeight,"totalWeight/D");
	TBranch *b_HT = BaselineTree->Branch("HT",&HT,"HT/F");
	TBranch *b_jet_n = BaselineTree->Branch("jet_n",&jet_n,"jet_n/I");
	TBranch *b_bjet_n = BaselineTree->Branch("bjet_n",&bjet_n,"bjet_n/I");
	TBranch *b_lep_n = BaselineTree->Branch("lep_n",&lep_n,"lep_n/I");
	TBranch *b_mll = BaselineTree->Branch("mll",&mll,"mll/F");
	TBranch *b_gamma_pt = BaselineTree->Branch("gamma_pt",&gamma_pt,"gamma_pt/F");
	TBranch *b_gamma_pt_smear = BaselineTree->Branch("Z_pt",&gamma_pt_smear,"Z_pt/F");
	TBranch *b_gamma_phi_smear = BaselineTree->Branch("Z_phi",&gamma_phi_smear,"Z_phi/F");
	TBranch *b_gamma_eta_smear = BaselineTree->Branch("Z_eta",&gamma_eta,"Z_eta/F");
	// TBranch *b_truthGamma_pt = BaselineTree->Branch("Z_truthPt",&truthGamma_pt,"Z_truthPt/F");
	// TBranch *b_truthGamma_eta = BaselineTree->Branch("Z_truthEta",&truthGamma_eta,"Z_truthEta/F");
	// TBranch *b_truthGamma_phi = BaselineTree->Branch("Z_truthPhi",&truthGamma_phi,"Z_truthPhi/F");
	//-----------EDIT HERE TO ADD VARS FOR 2019 RJR ANALYSIS
	TBranch *b_MET_smear = BaselineTree->Branch("MET",&MET_smear,"MET/F");
        TBranch *b_MET_loose = BaselineTree->Branch("MET_loose",&MET_loose,"MET_loose/F");
        TBranch *b_MET_tight = BaselineTree->Branch("MET_tight",&MET_tight,"MET_tight/F");
        TBranch *b_MET_tighter = BaselineTree->Branch("MET_tighter",&MET_tighter,"MET_tighter/F");
        TBranch *b_MET_tenacious = BaselineTree->Branch("MET_tenacious",&MET_tenacious,"MET_tenacious/F");
        TBranch *b_is2Lep2Jet =BaselineTree->Branch("is2Lep2Jet",&is2Lep2Jet,"is2Lep2Jet/I");
        TBranch *b_is2L2JInt =BaselineTree->Branch("is2L2JInt",&is2L2JInt,"is2L2JInt/I");
        TBranch *b_nBJet20_MV2c10_FixedCutBEff_77 =BaselineTree->Branch("nBJet20_MV2c10_FixedCutBEff_77",&nBJet20_MV2c10_FixedCutBEff_77,"nBJet20_MV2c10_FixedCutBEff_77/I");
        TBranch *b_mjj =BaselineTree->Branch("mjj",&mjj,"mjj/F");
        TBranch *b_mll_RJ =BaselineTree->Branch("mll_RJ",&mll_RJ,"mll_RJ/F");
        TBranch *b_R_minH2P_minH3P =BaselineTree->Branch("R_minH2P_minH3P",&R_minH2P_minH3P,"R_minH2P_minH3P/F");
        TBranch *b_RPT_HT5PP =BaselineTree->Branch("RPT_HT5PP",&RPT_HT5PP,"RPT_HT5PP/F");
        TBranch *b_dphiVP =BaselineTree->Branch("dphiVP",&dphiVP,"dphiVP/F");
        TBranch *b_H2PP =BaselineTree->Branch("H2PP",&H2PP,"H2PP/F");
        TBranch *b_H5PP =BaselineTree->Branch("H5PP",&H5PP,"H5PP/F");
        TBranch *b_nJet20 =BaselineTree->Branch("nJet20",&nJet20,"nJet20/I");
        TBranch *b_minDphi =BaselineTree->Branch("minDphi",&minDphi,"minDphi/F");
        TBranch *b_MZ =BaselineTree->Branch("MZ",&MZ,"MZ/F");
        TBranch *b_NjS =BaselineTree->Branch("NjS",&NjS,"NjS/I");
        TBranch *b_NjISR =BaselineTree->Branch("NjISR",&NjISR,"NjISR/I");
        TBranch *b_dphiISRI =BaselineTree->Branch("dphiISRI",&dphiISRI,"dphiISRI/F");
        TBranch *b_RISR =BaselineTree->Branch("RISR",&RISR,"RISR/F");
        TBranch *b_PTISR =BaselineTree->Branch("PTISR",&PTISR,"PTISR/F");
        TBranch *b_PTI =BaselineTree->Branch("PTI",&PTI,"PTI/F");
        TBranch *b_PTCM =BaselineTree->Branch("PTCM",&PTCM,"PTCM/F");
        TBranch *b_MJ =BaselineTree->Branch("MJ",&MJ,"MJ/F");
        TBranch *b_is3Lep3Jet =BaselineTree->Branch("is3Lep3Jet",&is3Lep3Jet,"is3Lep3Jet/I");
        TBranch *b_is4Lep3Jet =BaselineTree->Branch("is4Lep3Jet",&is4Lep3Jet,"is4Lep3Jet/I");
        TBranch *b_lept1sign_VR =BaselineTree->Branch("lept1sign_VR",&lept1sign_VR,"lept1sign_VR/I");
        TBranch *b_lept2sign_VR =BaselineTree->Branch("lept2sign_VR",&lept2sign_VR,"lept2sign_VR/I");
        TBranch *b_lept1Pt_VR =BaselineTree->Branch("lept1Pt_VR",&lept1Pt_VR,"lept1Pt_VR/F");
        TBranch *b_lept2Pt_VR =BaselineTree->Branch("lept2Pt_VR",&lept2Pt_VR,"lept2Pt_VR/F");
        TBranch *b_MZ_VR =BaselineTree->Branch("MZ_VR",&MZ_VR,"MZ_VR/F");
        TBranch *b_MJ_VR =BaselineTree->Branch("MJ_VR",&MJ_VR,"MJ_VR/F");
        TBranch *b_RISR_VR =BaselineTree->Branch("RISR_VR",&RISR_VR,"RISR_VR/F");
        TBranch *b_PTISR_VR =BaselineTree->Branch("PTISR_VR",&PTISR_VR,"PTISR_VR/F");
        TBranch *b_PTI_VR =BaselineTree->Branch("PTI_VR",&PTI_VR,"PTI_VR/F");
        TBranch *b_PTCM_VR =BaselineTree->Branch("PTCM_VR",&PTCM_VR,"PTCM_VR/F");
        TBranch *b_dphiISRI_VR =BaselineTree->Branch("dphiISRI_VR",&dphiISRI_VR,"dphiISRI_VR/F");
        TBranch *b_lepFlavor =BaselineTree->Branch("lepFlavor","std::vector<int>",&lepFlavor);
        TBranch *b_lepCharge =BaselineTree->Branch("lepCharge","std::vector<int>",&lepCharge);
	//---------------------------------------------------------
	TBranch *b_METl_smear = BaselineTree->Branch("METl",&METl_smear,"METl/F");
	TBranch *b_METt_smear = BaselineTree->Branch("METt",&METt_smear,"METt/F");
	TBranch *b_MET_phi_smear = BaselineTree->Branch("MET_phi",&MET_phi_smear,"MET_phi/F");
	TBranch *b_DPhi_METJetLeading_smear = BaselineTree->Branch("DPhi_METJetLeading",&DPhi_METJetLeading_smear,"DPhi_METJetLeading/F");
	TBranch *b_DPhi_METJetSecond_smear = BaselineTree->Branch("DPhi_METJetSecond",&DPhi_METJetSecond_smear,"DPhi_METJetSecond/F");
	TBranch *b_DPhi_METLepLeading_smear = BaselineTree->Branch("DPhi_METLepLeading",&DPhi_METLepLeading_smear,"DPhi_METLepLeading/F");
	TBranch *b_DPhi_METLepSecond_smear = BaselineTree->Branch("DPhi_METLepSecond",&DPhi_METLepSecond_smear,"DPhi_METLepSecond/F");
	TBranch *b_DPhi_METPhoton_smear = BaselineTree->Branch("DPhi_METPhoton",&DPhi_METPhoton_smear,"DPhi_METPhoton/F");
	TBranch *b_DR_Wmin2Jet = BaselineTree->Branch("DR_Wmin2Jet",&DR_Wmin2Jet,"DR_Wmin2Jet/F");
	TBranch *b_DR_J0J1 = BaselineTree->Branch("DR_J0J1",&DR_J0J1,"DR_J0J1/F");
	TBranch *b_mWmin = BaselineTree->Branch("mWmin",&mWmin,"mWmin/F");
	TBranch *b_Wmin_pt = BaselineTree->Branch("Wmin_pt",&Wmin_pt,"Wmin_pt/F");
	TBranch *b_Wmin_eta = BaselineTree->Branch("Wmin_eta",&Wmin_eta,"Wmin_eta/F");
	TBranch *b_DPhi_METWmin = BaselineTree->Branch("DPhi_METWmin",&DPhi_METWmin,"DPhi_METWmin/F");
	TBranch *b_DPhi_WminZ = BaselineTree->Branch("DPhi_WminZ",&DPhi_WminZ,"DPhi_WminZ/F");
	TBranch *b_mj0j1 = BaselineTree->Branch("mj0j1",&mj0j1,"mj0j1/F");
	TBranch *b_W01_pt = BaselineTree->Branch("W01_pt",&W01_pt,"W01_pt/F");
	TBranch *b_DPhi_METW01 = BaselineTree->Branch("DPhi_METW01",&DPhi_METW01,"DPhi_METW01/F");
	TBranch *b_DPhi_W01Z = BaselineTree->Branch("DPhi_W01Z",&DPhi_W01Z,"DPhi_W01Z/F");
	TBranch *b_DPhi_METNonWJet = BaselineTree->Branch("DPhi_METNonWJet",&DPhi_METNonWJet,"DPhi_METNonWJet/F");
	TBranch *b_NonWJet_pT = BaselineTree->Branch("NonWJet_pT",&NonWJet_pT,"NonWJet_pT/F");
	TBranch *b_DPhi_METNonWminJet = BaselineTree->Branch("DPhi_METNonWminJet",&DPhi_METNonWminJet,"DPhi_METNonWminJet/F");
	TBranch *b_NonWminJet_pT = BaselineTree->Branch("NonWminJet_pT",&NonWminJet_pT,"NonWminJet_pT/F");

	TBranch *b_lep_pT = BaselineTree->Branch("lep_pT","std::vector<float>",&lep_pT);
	TBranch *b_lep_phi = BaselineTree->Branch("lep_phi","std::vector<float>",&lep_phi);
	TBranch *b_lep_eta = BaselineTree->Branch("lep_eta","std::vector<float>",&lep_eta);
	TBranch *b_jet_pT = BaselineTree->Branch("jet_pT","std::vector<float>",&jet_pT);
	TBranch *b_jet_phi = BaselineTree->Branch("jet_phi","std::vector<float>",&jet_phi);
	TBranch *b_jet_eta = BaselineTree->Branch("jet_eta","std::vector<float>",&jet_eta);
	TBranch *b_jet_m = BaselineTree->Branch("jet_m","std::vector<float>",&jet_m);

	TBranch *b_MT2W = BaselineTree->Branch("MT2W",&MT2W,"MT2W/F");
	TBranch *b_DR_2Lep = BaselineTree->Branch("DR_2Lep",&DR_2Lep,"DR_2Lep/F");

	TH1D* hist_low_njet = new TH1D("hist_low_njet","",bin_size,njet_bin);
	hist_low_njet->SetStats(0);
	TH1D* hist_low_nbjet = new TH1D("hist_low_nbjet","",bin_size,njet_bin);
	hist_low_nbjet->SetStats(0);
	TH1D* hist_low_pt = new TH1D("hist_low_pt","",bin_size,sm_pt_bin);
	hist_low_pt->SetStats(0);
	TH1D* hist_sm_pt = new TH1D("hist_sm_pt","",bin_size,sm_pt_bin);
	hist_sm_pt->SetStats(0);
	TH1D* hist_low_et = new TH1D("hist_low_et","",bin_size,et_bin);
	hist_low_et->SetStats(0);
	TH1D* hist_low_ht = new TH1D("hist_low_ht","",bin_size,ht_bin);
	hist_low_ht->SetStats(0);
	TH1D* hist_medium_njet = new TH1D("hist_medium_njet","",bin_size,njet_bin);
	hist_medium_njet->SetStats(0);
	TH1D* hist_medium_nbjet = new TH1D("hist_medium_nbjet","",bin_size,njet_bin);
	hist_medium_nbjet->SetStats(0);
	TH1D* hist_medium_pt = new TH1D("hist_medium_pt","",bin_size,sm_pt_bin);
	hist_medium_pt->SetStats(0);
	TH1D* hist_medium_et = new TH1D("hist_medium_et","",bin_size,et_bin);
	hist_medium_et->SetStats(0);
	TH1D* hist_medium_ht = new TH1D("hist_medium_ht","",bin_size,ht_bin);
	hist_medium_ht->SetStats(0);
	TH1D* hist_high_njet = new TH1D("hist_high_njet","",bin_size,njet_bin);
	hist_high_njet->SetStats(0);
	TH1D* hist_high_nbjet = new TH1D("hist_high_nbjet","",bin_size,njet_bin);
	hist_high_nbjet->SetStats(0);
	TH1D* hist_high_pt = new TH1D("hist_high_pt","",bin_size,sm_pt_bin);
	hist_high_pt->SetStats(0);
	TH1D* hist_high_et = new TH1D("hist_high_et","",bin_size,et_bin);
	hist_high_et->SetStats(0);
	TH1D* hist_high_ht = new TH1D("hist_high_ht","",bin_size,ht_bin);
	hist_high_ht->SetStats(0);
	TH1D* hist_zmet_njet = new TH1D("hist_zmet_njet","",bin_size,njet_bin);
	hist_zmet_njet->SetStats(0);
	TH1D* hist_zmet_nbjet = new TH1D("hist_zmet_nbjet","",bin_size,njet_bin);
	hist_zmet_nbjet->SetStats(0);
	TH1D* hist_zmet_pt = new TH1D("hist_zmet_pt","",bin_size,sm_pt_bin);
	hist_zmet_pt->SetStats(0);
	TH1D* hist_zmet_et = new TH1D("hist_zmet_et","",bin_size,et_bin);
	hist_zmet_et->SetStats(0);
	TH1D* hist_zmet_ht = new TH1D("hist_zmet_ht","",bin_size,ht_bin);
	hist_zmet_ht->SetStats(0);
	TH1D* hist_bveto_njet = new TH1D("hist_bveto_njet","",bin_size,njet_bin);
	hist_bveto_njet->SetStats(0);
	TH1D* hist_bveto_nbjet = new TH1D("hist_bveto_nbjet","",bin_size,njet_bin);
	hist_bveto_nbjet->SetStats(0);
	TH1D* hist_bveto_pt = new TH1D("hist_bveto_pt","",bin_size,sm_pt_bin);
	hist_bveto_pt->SetStats(0);
	TH1D* hist_bveto_et = new TH1D("hist_bveto_et","",bin_size,et_bin);
	hist_bveto_et->SetStats(0);
	TH1D* hist_bveto_ht = new TH1D("hist_bveto_ht","",bin_size,ht_bin);
	hist_bveto_ht->SetStats(0);

	TH1D* hist_low_dpt = new TH1D("hist_low_dpt","",dpt_bin_size,dpt_bin);
	hist_low_dpt->SetStats(0);
	TH1D* hist_low_pt_smear = new TH1D("hist_low_pt_smear","",bin_size,sm_pt_bin);
	hist_low_pt_smear->SetStats(0);
	TH1D* hist_low_htincl = new TH1D("hist_low_htincl","",bin_size,ht_bin);
	hist_low_htincl->SetStats(0);
	TH1D* hist_medium_pt_smear = new TH1D("hist_medium_pt_smear","",bin_size,sm_pt_bin);
	hist_medium_pt_smear->SetStats(0);
	TH1D* hist_medium_htincl = new TH1D("hist_medium_htincl","",bin_size,ht_bin);
	hist_medium_htincl->SetStats(0);
	TH1D* hist_high_pt_smear = new TH1D("hist_high_pt_smear","",bin_size,sm_pt_bin);
	hist_high_pt_smear->SetStats(0);
	TH1D* hist_high_htincl = new TH1D("hist_high_htincl","",bin_size,ht_bin);
	hist_high_htincl->SetStats(0);
	TH1D* hist_zmet_pt_smear = new TH1D("hist_zmet_pt_smear","",bin_size,sm_pt_bin);
	hist_zmet_pt_smear->SetStats(0);
	TH1D* hist_zmet_htincl = new TH1D("hist_zmet_htincl","",bin_size,ht_bin);
	hist_zmet_htincl->SetStats(0);
	TH1D* hist_bveto_pt_smear = new TH1D("hist_bveto_pt_smear","",bin_size,sm_pt_bin);
	hist_bveto_pt_smear->SetStats(0);
	TH1D* hist_bveto_htincl = new TH1D("hist_bveto_htincl","",bin_size,ht_bin);
	hist_bveto_htincl->SetStats(0);

	TH1D* hist_low_met = new TH1D("hist_low_met","",bin_size,met_bin);
	hist_low_met->SetStats(0);
	TH1D* hist_low_dphi = new TH1D("hist_low_dphi","",bin_size,dphi_bin);
	hist_low_dphi->SetStats(0);

	//-----------------------------
	// loop over events
	//-----------------------------

	Long64_t nentries = T->GetEntries();

	for (Long64_t i=0;i<nentries;i++) {

		if (fmod(i,1e5)==0) std::cout << i << " events processed." << std::endl;
		T->GetEntry(i);

		if( jet_n == 0 ) continue;

		// use the smearing function to smear MET and pT in photon events
		float photon_smear = 0;
		float photon_smear_phi = 0;
		int smpt = hist_sm_pt->FindBin(gamma_pt)-1;
		if (smpt>=0) {
			if (smearing_method != 0) {
				if (smear_final[smpt]->Integral()>0) photon_smear = smear_final[smpt]->GetRandom() + shift[smpt];
				else photon_smear = shift[smpt];
				smear_shift = shift[smpt];
			}
			else {
				photon_smear = 0;
				smear_shift = 0;
			}
		}

		gamma_pt_smear = gamma_pt-photon_smear; // sign of photon_smear is important!!!
		gamma_phi_smear = gamma_phi-photon_smear_phi; // sign of photon_smear is important!!!

		float photon_smear_l = gamma_pt-gamma_pt_smear*TMath::Cos(photon_smear_phi);
		float photon_smear_t = -gamma_pt_smear*TMath::Sin(photon_smear_phi);
		METl_smear = METl + photon_smear_l;  // sign of photon_smear is important!!!
		METt_smear = METt + photon_smear_t;  // sign of photon_smear is important!!!
		MET_smear = pow(METl_smear*METl_smear+METt_smear*METt_smear,0.5);

		pt_smear = hist_sm_pt->FindBin(gamma_pt_smear)-1;
		if (gamma_pt_smear>sm_pt_bin[bin_size]) pt_smear = bin_size-1;
		met_smear = hist_low_met->FindBin(MET_smear)-1;
		if (met_smear>met_bin[bin_size]) met_smear = bin_size-1;

		// recompute DPhi after smearing
		float METtx = METt*TMath::Cos(gamma_phi_smear+TMath::Pi()/2.);
		float METty = METt*TMath::Sin(gamma_phi_smear+TMath::Pi()/2.);
		float METlx_smear = METl_smear*TMath::Cos(gamma_phi_smear);
		float METly_smear = METl_smear*TMath::Sin(gamma_phi_smear);

		TLorentzVector met_4vec_smear;
		met_4vec_smear.SetXYZM(METtx+METlx_smear,METty+METly_smear,0,0);
		MET_phi_smear = met_4vec_smear.Phi();

		TLorentzVector jet0_4vec;
		if (jet_n<1) jet0_4vec.SetPtEtaPhiM(0,0,0,0);
		else jet0_4vec.SetPtEtaPhiM(jet_pT->at(0),jet_eta->at(0),jet_phi->at(0),jet_m->at(0));

		//cout << __LINE__ << endl;
		//cout << "jet_n " << jet_n << endl;
		
		TLorentzVector jet1_4vec;
		if (jet_n<2)jet1_4vec.SetPtEtaPhiM(0,0,0,0);
		else jet1_4vec.SetPtEtaPhiM(jet_pT->at(1),jet_eta->at(1),jet_phi->at(1),jet_m->at(1));

		//		TLorentzVector jet1_4vec;
		//cout << __LINE__ << endl;
		if (jet_n<2){
		  //cout << __LINE__ << endl;
		  jet1_4vec.SetPtEtaPhiM(0,0,0,0);
		  
		}
		else{
		  //cout << __LINE__ << endl;
		  jet1_4vec.SetPtEtaPhiM(jet_pT->at(1),jet_eta->at(1),jet_phi->at(1),jet_m->at(1));
		}
		
		//cout << __LINE__ << endl;
		DPhi_METJetLeading_smear = fabs(met_4vec_smear.DeltaPhi(jet0_4vec));
		DPhi_METJetSecond_smear = fabs(met_4vec_smear.DeltaPhi(jet1_4vec));
		
		//DPhi_METJetLeading_smear = -999;
		//DPhi_METJetSecond_smear  = -999;

		DPhi_METPhoton_smear = fabs(TMath::ATan2(METt,METl_smear));
		dphi_smear = hist_low_dphi->FindBin(DPhi_METPhoton_smear)-1;
		if (dphi_smear>dphi_bin[bin_size]) dphi_smear = bin_size-1;

		if (gamma_pt>50. && jet_n==1) g_metl_smear[smpt]->Fill(METl_smear,totalWeight);
		if (gamma_pt>50. && jet_n>=2) g_metl_smear_2j[smpt]->Fill(METl_smear,totalWeight);

		// translate photon pT to dilepton sum pT, and compute HTincl for photon events
		float photon_2LPt = 0;
		//if (pt_smear>=0) if (lep_sumpt[pt_smear]->Integral()>0) photon_2LPt = lep_sumpt[pt_smear]->GetRandom();
		HTincl = HT + photon_2LPt;
		//dpt = hist_low_dpt->FindBin(gamma_pt_smear-gamma_pt)-1;
		dpt = hist_low_dpt->FindBin(METl_smear)-1;
		if (dpt>=0 && pt_smear>=0) if (hist_Mll_dPt[pt_smear][dpt]->Integral()>0) mll = hist_Mll_dPt[pt_smear][dpt]->GetRandom();

		//---------------------------------------------
		// compute two lepton kinematics
		//---------------------------------------------
		TLorentzVector z_4vec;
		z_4vec.SetPtEtaPhiM(gamma_pt,gamma_eta,gamma_phi,mll);
		//z_4vec.SetPtEtaPhiM(gamma_pt_smear,gamma_eta,gamma_phi,mll);
		GetDijetVariables(z_4vec,met_4vec_smear);

		lep_pT->clear();
		lep_eta->clear();
		lep_phi->clear();
		lep_pT->push_back(0);
		lep_eta->push_back(0);
		lep_phi->push_back(0);
		lep_pT->push_back(0);
		lep_eta->push_back(0);
		lep_phi->push_back(0);
		int ntry = 0;
		while ((lep_pT->at(0)<leading_lep_pt_cut || lep_pT->at(1)<second_lep_pt_cut) && ntry<100) {
			ntry += 1;
			GetIndividualLeptonInfo(z_4vec);
		}

		TLorentzVector lep0_4vec;
		TLorentzVector lep1_4vec;
		lep_n = 2;
		lep0_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
		lep1_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);
		MT2W = ComputeMT2(lep0_4vec, lep1_4vec, met_4vec_smear, 0, 0).Compute();
		DPhi_METLepLeading_smear = fabs(met_4vec_smear.DeltaPhi(lep0_4vec));
		DPhi_METLepSecond_smear = fabs(met_4vec_smear.DeltaPhi(lep1_4vec));
/////////////////////////////////////////////////////////////////////////////////////////
//		cout << "DEBUGGING - leading: " << DPhi_METLepLeading_smear << " second: " << DPhi_METLepSecond_smear << endl;
/////////////////////////////////////////////////////////////////////////////////////////
		DR_2Lep = lep0_4vec.DeltaR(lep1_4vec);

		BaselineTree->Fill();

	}

	BaselineTree->Write();

	for (int bin=0;bin<bin_size;bin++) {
		g_metl[bin]->Write();
		z_metl[bin]->Write();
		z_metl_2j[bin]->Write();
		g_metl_smear[bin]->Write();
		g_metl_smear_2j[bin]->Write();
		smear_final[bin]->Write();
		if (smearing_method != 0) {
			g_resp[bin]->Write();
			z_resp[bin]->Write();
			smear_raw[bin]->Write();
			smear_fft_re[bin]->Write();
			smear_fft_im[bin]->Write();
			smear_fft_amp[bin]->Write();
		}
	}

	std::cout << "done." << std::endl;
	delete f;


}

//-----------------------------------------------------------------------------------------------
// this script takes the outputs from GetBaseLineEvents.C and GetPhotonEvents.C, and makes smeared photon information.
// the parameters of the function GetPhotonSmearing(string ch, bool isData) are:
// 	ch: which dilepton channel you are smearing the photon events to (ee,mm)
// 	isData: put "true" if you are running a data sample
// example of code running command: root -l -b -q 'GetPhotonSmearing.C+("ee",true)'
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

#include "BasicSetting.C"
#include "PhotonVariables.C"
#include "GetDijetVariables.C"
#include "GetJpsiSmearingHistogram.C"
#include "GetJpsiMllHistogram.C"
#include "GetIndividualLeptonInfo.C"
#include "MT2_ROOT.h"

using namespace std;

vector<string> noSampleWeight;
vector<string> noEventWeight;


int RebinHistogram(TH1D* hist, int rebin) {
	double negative_yield = 0.;
	double positive_yield = 0.;
	double core_yield = 0.;
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
	std::cout  << "negative_yield = " << negative_yield << ", positive_yield = " << positive_yield << ", remainder = " << remainder << std::endl;
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
void GetJpsiSmearing(string ch, int isData, string year) {


	//-----------------------------
	// prepare lep pT functions, to convert photon pT to dilepton sum pT
	//-----------------------------
	
	std::cout << "Prepare Mll histograms..." << std::endl;

	GetJpsiMllHistogram(ch);
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

	TSpectrum pfinder;
	TH1D* g_resp[bin_size];
	TH1D* z_resp[bin_size];
	TH1D* smear_raw[bin_size];
	TH1D* smear_fft_re[bin_size];
	TH1D* smear_fft_im[bin_size];
	TH1D* smear_fft_amp[bin_size];
	TH1D* smear_final[bin_size];
	TH1D* smear_final_phi[bin_size];
	double shift[bin_size];
	TH1D* g_metl_smear[bin_size];
	TH1D* g_metl_smear_2j[bin_size];
	for (int bin=0;bin<bin_size;bin++) {
		g_metl_smear[bin] = new TH1D(TString("g_metl_smear_")+TString::Itoa(bin,10),"",40000,-30000,10000);
		g_metl_smear[bin]->SetStats(0);
		g_metl_smear_2j[bin] = new TH1D(TString("g_metl_smear_2j_")+TString::Itoa(bin,10),"",40000,-30000,10000);
		g_metl_smear_2j[bin]->SetStats(0);
	}


	if (smearing_method != 0) GetJpsiSmearingHistogram(ch,lumi, photon_tag);

	if (smearing_method != 0) {  // if you want to use the deconvolution method to smear photon events. To enable this method, set "bool useDeconvolution = true" in BasicSetting.C
		for (int bin=0;bin<bin_size;bin++) {
			int rebin = 1;
			rebin = RebinHistogram(z_metl[bin],0);
			rebin = RebinHistogram(z_metl_2j[bin],rebin);
			rebin = RebinHistogram(g_metl[bin],rebin);
			rebin = RebinHistogram(g_metl_smear[bin],rebin);
			rebin = RebinHistogram(g_metl_smear_2j[bin],rebin);
			double gmetl_mean = g_metl[bin]->GetMean();
			double gmetl_rms = g_metl[bin]->GetRMS();
			double zmetl_mean = z_metl[bin]->GetMean();
			double zmetl_rms = z_metl[bin]->GetRMS();
			std::cout << "--------------------------------------------------------" << std::endl;
			std::cout << "bin " << bin << ", rebin " << rebin << std::endl;
			std::cout << "gmetl_mean " << gmetl_mean << ", gmetl_rms " << gmetl_rms << std::endl;
			std::cout << "zmetl_mean " << zmetl_mean << ", zmetl_rms " << zmetl_rms << std::endl;
			const int newbin = 40000/rebin;
			Double_t *smear = new Double_t[2*((newbin+1)/2+1)];
			Double_t *fft_re = new Double_t[newbin];
			Double_t *fft_im = new Double_t[newbin];
			Double_t *z_smear_in = new Double_t[newbin];
			Double_t *g_smear_in = new Double_t[newbin];
			Double_t *j_resp_in = new Double_t[newbin];
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
					z_smear_in[i] = max(z_dpt[bin]->GetBinContent(i+1),0.);
					if (i<newbin/2) g_smear_in[i] = max(g_dpt[bin]->GetBinContent(i+1+newbin/2),0.);
					else g_smear_in[i] = 0.;
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
			pfinder.Deconvolution(z_smear_in,g_smear_in,newbin,1000,1,1);
			pfinder.Deconvolution(z_resp_in,j_resp_in,newbin,1000,1,1);
			pfinder.Deconvolution(g_resp_in,j_resp_in,newbin,1000,1,1);
			for (int i=0;i<newbin;i++) {
				smear_raw[bin]->SetBinContent(i+1,z_smear_in[i]);
				g_resp[bin]->SetBinContent(i+1,g_resp_in[i]);
				z_resp[bin]->SetBinContent(i+1,z_resp_in[i]);
				smear[i] = z_smear_in[i];
			}
			rebin = RebinHistogram(z_dpt[bin],0);
			rebin = RebinHistogram(g_dpt[bin],rebin);
			double gdpt_mean = g_dpt[bin]->GetMean();
			double gdpt_rms = g_dpt[bin]->GetRMS();
			double zdpt_mean = z_dpt[bin]->GetMean();
			double zdpt_rms = z_dpt[bin]->GetRMS();
			double smear_mean = smear_raw[bin]->GetMean();
			double smear_rms = smear_raw[bin]->GetRMS();
			std::cout << "gdpt_mean " << gdpt_mean << ", gdpt_rms " << gdpt_rms << std::endl;
			std::cout << "zdpt_mean " << zdpt_mean << ", zdpt_rms " << zdpt_rms << std::endl;
			std::cout << "smear_mean " << smear_mean << ", smear_rms " << smear_rms << std::endl;
			for (int i=0;i<newbin;i++) {
				if (gmetl_rms/zmetl_rms > 1.0) {
					smear_raw[bin]->SetBinContent(i+1,0.);
					smear[i] = 0.;
				}
				double smear_cut = 6.;
				if (ch=="mm" && bphys_pt_bin[bin]>=0) smear_cut = 7.;
				if (abs(smear_raw[bin]->GetBinCenter(i+1)-smear_mean)/smear_rms>smear_cut) {
					smear_raw[bin]->SetBinContent(i+1,0.);
					smear[i] = 0.;
				}
			}
			if (smearing_method == 3) shift[bin] = -g_dpt[bin]->GetMean();
			else shift[bin] = -g_metl[bin]->GetMean();
		}
	}
	for (int bin=0;bin<bin_size;bin++) {
		smear_final[bin] = new TH1D(TString("smear_final_")+TString::Itoa(bin,10),"",500,-1000,1000);
		if (smearing_method != 0) {
	 		for (int i=0;i<500;i++) {
				int which_bin = smear_raw[bin]->FindBin(smear_final[bin]->GetBinCenter(i+1));
				smear_final[bin]->SetBinContent(i+1,smear_raw[bin]->GetBinContent(which_bin));
	 		}
		}
	}

	std::cout << "Finished preparing smearing histograms..." << std::endl;

	//---------------------------------------------
	// open file, get Tree and EventCountHist
	//---------------------------------------------

	TH1::SetDefaultSumw2();

	string  filename;
	if (isData==0) filename = TString(TString(smearingPath)+"bphys/bdata_raw.root"); 
	if (isData==1) filename = TString(TString(smearingPath)+"bphys/bdata_raw.root"); 
	if (isData==2) filename = TString(TString(smearingPath)+"bphys/bdata_raw.root"); 
	TFile*  inputFile      = TFile::Open(filename.c_str());
	TTree*  T              = (TTree*)inputFile->Get("BaselineTree");

	cout << endl;
	cout << "Opening file           : " << filename        << endl;
	cout << "Events in ntuple       : " << T->GetEntries() << endl;

	//---------------------------------------------
	// open file, get Tree and EventCountHist
	//---------------------------------------------

	TH1::SetDefaultSumw2();

	if (isData==0) filename = TString(TString(outputPath)+"bphys/bdata_"+TString(ch)+TString(photon_tag)+".root"); 
	if (isData==1) filename = TString(TString(outputPath)+"bphys/jpsidata_"+TString(ch)+TString(photon_tag)+".root"); 
	if (isData==2) filename = TString(TString(outputPath)+"bphys/upsidata_"+TString(ch)+TString(photon_tag)+".root"); 
	TFile*  f              = new TFile(filename.c_str(),"recreate");          
	TTree* BaselineTree = new TTree("BaselineTree","baseline tree");

	cout << endl;
	cout << "Create file           : " << filename        << endl;


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
	//if (isData!=1) T->SetBranchAddress("truthGamma_pt", &truthGamma_pt);
	//if (isData!=1) T->SetBranchAddress("truthGamma_phi", &truthGamma_phi);
	T->SetBranchAddress("Jpsi_pt", &gamma_pt);
	T->SetBranchAddress("Jpsi_eta", &gamma_eta);
	T->SetBranchAddress("Jpsi_phi", &gamma_phi);
	T->SetBranchAddress("MET_raw"             ,&MET              );
	T->SetBranchAddress("METl_raw"             ,&METl              );
	T->SetBranchAddress("METt_raw"             ,&METt              );
	T->SetBranchAddress("MET_phi_raw"         ,&MET_phi          );
	T->SetBranchAddress("DPhi_METJetLeading_raw"         ,&DPhi_METJetLeading          );
	T->SetBranchAddress("DPhi_METJetSecond_raw"          ,&DPhi_METJetSecond           );
	T->SetBranchAddress("jet_pT"          ,&jet_pT           );
	T->SetBranchAddress("jet_eta"         ,&jet_eta          );
	T->SetBranchAddress("jet_phi"         ,&jet_phi          );
	T->SetBranchAddress("jet_m"         ,&jet_m          );
	double Jpsi_mll = 0;
	T->SetBranchAddress("Jpsi_mll", &Jpsi_mll);

	//-----------------------------
	// add new branches
	//-----------------------------

	TBranch *b_totalWeight = BaselineTree->Branch("totalWeight",&totalWeight,"totalWeight/D");
	TBranch *b_HT = BaselineTree->Branch("HT",&HT,"HT/D");
	TBranch *b_jet_n = BaselineTree->Branch("jet_n",&jet_n,"jet_n/I");
	TBranch *b_bjet_n = BaselineTree->Branch("bjet_n",&bjet_n,"bjet_n/I");
	TBranch *b_lep_n = BaselineTree->Branch("lep_n",&lep_n,"lep_n/I");
	TBranch *b_mll = BaselineTree->Branch("mll",&mll,"mll/D");
	TBranch *b_Jpsi_mll = BaselineTree->Branch("Jpsi_mll",&Jpsi_mll,"Jpsi_mll/D");
	TBranch *b_Jpsi_pt = BaselineTree->Branch("Jpsi_pt",&gamma_pt,"Jpsi_pt/D");
	TBranch *b_Jpsi_pt_smear = BaselineTree->Branch("Z_pt",&gamma_pt_smear,"Z_pt/D");
	TBranch *b_Jpsi_phi_smear = BaselineTree->Branch("Z_phi",&gamma_phi_smear,"Z_phi/D");
	TBranch *b_Jpsi_eta_smear = BaselineTree->Branch("Z_eta",&gamma_eta,"Z_eta/D");
	TBranch *b_MET_smear = BaselineTree->Branch("MET",&MET_smear,"MET/D");
	TBranch *b_METl_smear = BaselineTree->Branch("METl",&METl_smear,"METl/D");
	TBranch *b_METt_smear = BaselineTree->Branch("METt",&METt_smear,"METt/D");
	TBranch *b_MET_phi_smear = BaselineTree->Branch("MET_phi",&MET_phi_smear,"MET_phi/D");
	TBranch *b_DPhi_METJetLeading_smear = BaselineTree->Branch("DPhi_METJetLeading",&DPhi_METJetLeading_smear,"DPhi_METJetLeading/D");
	TBranch *b_DPhi_METJetSecond_smear = BaselineTree->Branch("DPhi_METJetSecond",&DPhi_METJetSecond_smear,"DPhi_METJetSecond/D");
	TBranch *b_DPhi_METLepLeading_smear = BaselineTree->Branch("DPhi_METLepLeading",&DPhi_METLepLeading_smear,"DPhi_METLepLeading/D");
	TBranch *b_DPhi_METLepSecond_smear = BaselineTree->Branch("DPhi_METLepSecond",&DPhi_METLepSecond_smear,"DPhi_METLepSecond/D");
	TBranch *b_DPhi_METPhoton_smear = BaselineTree->Branch("DPhi_METPhoton",&DPhi_METPhoton_smear,"DPhi_METPhoton/D");
	TBranch *b_DR_Wmin2Jet = BaselineTree->Branch("DR_Wmin2Jet",&DR_Wmin2Jet,"DR_Wmin2Jet/D");
	TBranch *b_DR_J0J1 = BaselineTree->Branch("DR_J0J1",&DR_J0J1,"DR_J0J1/D");
	TBranch *b_mWmin = BaselineTree->Branch("mWmin",&mWmin,"mWmin/D");
	TBranch *b_Wmin_pt = BaselineTree->Branch("Wmin_pt",&Wmin_pt,"Wmin_pt/D");
	TBranch *b_Wmin_eta = BaselineTree->Branch("Wmin_eta",&Wmin_eta,"Wmin_eta/D");
	TBranch *b_DPhi_METWmin = BaselineTree->Branch("DPhi_METWmin",&DPhi_METWmin,"DPhi_METWmin/D");
	TBranch *b_DPhi_WminZ = BaselineTree->Branch("DPhi_WminZ",&DPhi_WminZ,"DPhi_WminZ/D");
	TBranch *b_mj0j1 = BaselineTree->Branch("mj0j1",&mj0j1,"mj0j1/D");
	TBranch *b_W01_pt = BaselineTree->Branch("W01_pt",&W01_pt,"W01_pt/D");
	TBranch *b_DPhi_METW01 = BaselineTree->Branch("DPhi_METW01",&DPhi_METW01,"DPhi_METW01/D");
	TBranch *b_DPhi_W01Z = BaselineTree->Branch("DPhi_W01Z",&DPhi_W01Z,"DPhi_W01Z/D");
	TBranch *b_DPhi_METNonWJet = BaselineTree->Branch("DPhi_METNonWJet",&DPhi_METNonWJet,"DPhi_METNonWJet/D");
	TBranch *b_NonWJet_pT = BaselineTree->Branch("NonWJet_pT",&NonWJet_pT,"NonWJet_pT/D");
	TBranch *b_DPhi_METNonWminJet = BaselineTree->Branch("DPhi_METNonWminJet",&DPhi_METNonWminJet,"DPhi_METNonWminJet/D");
	TBranch *b_NonWminJet_pT = BaselineTree->Branch("NonWminJet_pT",&NonWminJet_pT,"NonWminJet_pT/D");

	TBranch *b_lep_pT = BaselineTree->Branch("lep_pT","std::vector<double>",&lep_pT);
	TBranch *b_lep_phi = BaselineTree->Branch("lep_phi","std::vector<double>",&lep_phi);
	TBranch *b_lep_eta = BaselineTree->Branch("lep_eta","std::vector<double>",&lep_eta);
	TBranch *b_jet_pT = BaselineTree->Branch("jet_pT","std::vector<double>",&jet_pT);
	TBranch *b_jet_phi = BaselineTree->Branch("jet_phi","std::vector<double>",&jet_phi);
	TBranch *b_jet_eta = BaselineTree->Branch("jet_eta","std::vector<double>",&jet_eta);
	TBranch *b_jet_m = BaselineTree->Branch("jet_m","std::vector<double>",&jet_m);

	TBranch *b_MT2W = BaselineTree->Branch("MT2W",&MT2W,"MT2W/Float_t");
	TBranch *b_DR_2Lep = BaselineTree->Branch("DR_2Lep",&DR_2Lep,"DR_2Lep/D");


	//-----------------------------
	// loop over events
	//-----------------------------

	Long64_t nentries = T->GetEntries();

	for (Long64_t i=0;i<nentries;i++) {

		if (fmod(i,1e5)==0) std::cout << i << " events processed." << std::endl;
		T->GetEntry(i);

		if (isData==1 && Jpsi_mll>5.) continue; 
		if (isData==2 && Jpsi_mll<5.) continue; 

		// use the smearing function to smear MET and pT in photon events
		if (fmod(i,1e5)==0) std::cout << i << " smear photon." << std::endl;
		double photon_smear = 0;
		double photon_smear_phi = 0;
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

		if (fmod(i,1e5)==0) std::cout << i << " smear MET" << std::endl;
		double photon_smear_l = gamma_pt-gamma_pt_smear*TMath::Cos(photon_smear_phi);
		double photon_smear_t = -gamma_pt_smear*TMath::Sin(photon_smear_phi);
		METl_smear = METl + photon_smear_l;  // sign of photon_smear is important!!!
		METt_smear = METt + photon_smear_t;  // sign of photon_smear is important!!!
		MET_smear = pow(METl_smear*METl_smear+METt_smear*METt_smear,0.5);

		pt_smear = hist_sm_pt->FindBin(gamma_pt_smear)-1;
		if (gamma_pt_smear>bphys_sm_pt_bin[bin_size]) pt_smear = bin_size-1;

		// recompute DPhi after smearing
		if (fmod(i,1e5)==0) std::cout << i << " compute DPhi" << std::endl;
		double METtx = METt*TMath::Cos(gamma_phi_smear+TMath::Pi()/2.);
		double METty = METt*TMath::Sin(gamma_phi_smear+TMath::Pi()/2.);
		double METlx_smear = METl_smear*TMath::Cos(gamma_phi_smear);
		double METly_smear = METl_smear*TMath::Sin(gamma_phi_smear);
		TLorentzVector met_4vec_smear;
		met_4vec_smear.SetXYZM(METtx+METlx_smear,METty+METly_smear,0,0);
		MET_phi_smear = met_4vec_smear.Phi();
		TLorentzVector jet0_4vec;
		if (jet_n<1) jet0_4vec.SetPtEtaPhiM(0,0,0,0);
		else jet0_4vec.SetPtEtaPhiM(jet_pT->at(0),jet_eta->at(0),jet_phi->at(0),jet_m->at(0));
		TLorentzVector jet1_4vec;
		if (jet_n<2) jet1_4vec.SetPtEtaPhiM(0,0,0,0);
		else jet1_4vec.SetPtEtaPhiM(jet_pT->at(1),jet_eta->at(1),jet_phi->at(1),jet_m->at(1));
		DPhi_METJetLeading_smear = fabs(met_4vec_smear.DeltaPhi(jet0_4vec));
		DPhi_METJetSecond_smear = fabs(met_4vec_smear.DeltaPhi(jet1_4vec));
		DPhi_METPhoton_smear = fabs(TMath::ATan2(METt,METl_smear));

		// translate photon pT to dilepton sum pT, and compute HTincl for photon events
		double photon_2LPt = 0;
		//if (pt_smear>=0) if (lep_sumpt[pt_smear]->Integral()>0) photon_2LPt = lep_sumpt[pt_smear]->GetRandom();
		HTincl = HT + photon_2LPt;
		//dpt = hist_low_dpt->FindBin(gamma_pt_smear-gamma_pt)-1;
		if (fmod(i,1e5)==0) std::cout << i << " compute mll." << std::endl;
		dpt = hist_low_dpt->FindBin(METl_smear)-1;
		if (dpt>=0 && pt_smear>=0) if (hist_Mll_dPt[pt_smear][dpt]->Integral()>0) mll = hist_Mll_dPt[pt_smear][dpt]->GetRandom();

		//---------------------------------------------
		// compute dijet system variables, m80jj, W pT, DR(2jet), etc.
		//---------------------------------------------
		if (fmod(i,1e5)==0) std::cout << i << " compute dilepton variables" << std::endl;
		TLorentzVector z_4vec;
		z_4vec.SetPtEtaPhiM(gamma_pt_smear,gamma_eta,gamma_phi,mll);
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

		if (fmod(i,1e5)==0) std::cout << i << " compute MT2 variable" << std::endl;
		TLorentzVector lep0_4vec;
		TLorentzVector lep1_4vec;
		lep_n = 2;
		lep0_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
		lep1_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);
		MT2W = ComputeMT2(lep0_4vec, lep1_4vec, met_4vec_smear, 0, 0).Compute();
		DPhi_METLepLeading_smear = fabs(met_4vec_smear.DeltaPhi(lep0_4vec));
		DPhi_METLepSecond_smear = fabs(met_4vec_smear.DeltaPhi(lep1_4vec));
		DR_2Lep = lep0_4vec.DeltaR(lep1_4vec);

		if (fmod(i,1e5)==0) std::cout << i << " finish this event." << std::endl;
		BaselineTree->Fill();

	}

	BaselineTree->Write();

	if (smearing_method != 0) {
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
	}

	std::cout << "done." << std::endl;
	delete f;


}

//-----------------------------------------------------------------------------------------------
// this script takes the outputs from GetBaseLineEvents.C and GetPhotonEvents.C and GetPhotonSmearing.C, and makes photon reweighting factors.
// the parameters of the function GetPhotonReweighting(string ch, bool isData) are:
// 	ch: which dilepton channel you are smearing the photon events to (ee,mm)
// 	isData: put "true" if you are running a data sample
// example of code running command: root -l -b -q 'GetPhotonReweighting.C+("ee",true)'
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
#include "TGraphAsymmErrors.h"

#include "BasicSetting.C"
#include "GetJpsiReweightingHistogram.C"

using namespace std;

vector<string> noSampleWeight;
vector<string> noEventWeight;

void GetJpsiReweighting(string ch, int isData, string year) {

	if (smearing_method == 0) photon_tag = "_NoSmear";
	if (smearing_method == 1) photon_tag = "_McSmear";
	if (smearing_method == 2) photon_tag = "_DataSmear";
	if (smearing_method == 3) photon_tag = "_TruthSmear";

	// retrieve the histograms for pT, Njets, HT-reweighting and on-Z rescaling
	std::cout << "Prepare reweighting histograms..." << std::endl;
	GetJpsiReweightingHistogram(ch, isData,lumi, photon_tag);

	//---------------------------------------------
	// open file, get Tree and EventCountHist
	//---------------------------------------------

	TH1::SetDefaultSumw2();

	string  filename;
	if (isData==0) filename = TString(TString(outputPath)+"bphys/bdata_"+TString(ch)+TString(photon_tag)+".root"); 
	if (isData==1) filename = TString(TString(outputPath)+"bphys/jpsidata_"+TString(ch)+TString(photon_tag)+".root"); 
	if (isData==2) filename = TString(TString(outputPath)+"bphys/upsidata_"+TString(ch)+TString(photon_tag)+".root"); 
	TFile*  f              = new TFile(filename.c_str(),"update");          
	TTree*  outputTree              = (TTree*)f->Get("BaselineTree");

	cout << endl;
	cout << "Opening file           : " << filename        << endl;
	cout << "Events in ntuple       : " << outputTree->GetEntries() << endl;

	//if( _nGenEvents < 1.0 ) cout << "ERROR 0 EVENTS -----------------------------------------------------" << endl;

	//-----------------------------
	// access existing branches
	//-----------------------------

	double totalWeight = 0.;
	int pt = 0;
	int ptsm = 0;
	int metsm = 0;
	int ht = 0;
	int njet = 0;
	int nbjet = 0;
	Int_t jet_n = 0;
	double mll = 0.;
	double HT = 0.;
	double Jpsi_pt = 0.;
	double Jpsi_pt_smear = 0.;
	double MET_smear = 0.;
	double Jpsi_dR = 0.;
	std::vector<double>* lep_pT = new std::vector<double>(10);
	outputTree->SetBranchAddress("totalWeight"     ,&totalWeight              );
	outputTree->SetBranchAddress("pt"              ,&pt               );
	outputTree->SetBranchAddress("pt_smear"   ,&ptsm               );
	outputTree->SetBranchAddress("met_smear"   ,&metsm               );
	outputTree->SetBranchAddress("ht"              ,&ht               );
	outputTree->SetBranchAddress("njet"            ,&njet             );
	outputTree->SetBranchAddress("nbjet"            ,&nbjet             );
	outputTree->SetBranchAddress("mll"              ,&mll               );
	outputTree->SetBranchAddress("HT"              ,&HT               );
	outputTree->SetBranchAddress("jet_n"           ,&jet_n            );
	outputTree->SetBranchAddress("Jpsi_pt", &Jpsi_pt);
	outputTree->SetBranchAddress("Z_pt", &Jpsi_pt_smear);
	outputTree->SetBranchAddress("MET", &MET_smear);
	outputTree->SetBranchAddress("lep_pT"          ,&lep_pT           );

	//-----------------------------
	// add new branches
	//-----------------------------

	double ptrw_bveto = 0;
	double ptsmrw_bveto = 0;
	double etrw_bveto = 0;
	double etsmrw_bveto = 0;
	double htrw_bveto = 0;
	double njetrw_bveto = 0;
	double nbjetrw_bveto = 0;
	double ptrw_err_bveto = 0;
	double ptsmrw_err_bveto = 0;
	double etrw_err_bveto = 0;
	double etsmrw_err_bveto = 0;
	double htrw_err_bveto = 0;
	double njetrw_err_bveto = 0;
	double nbjetrw_err_bveto = 0;
	TBranch *b_ptrw_bveto = outputTree->Branch("ptrw_bveto",&ptrw_bveto,"ptrw_bveto/D");
	TBranch *b_ptsmrw_bveto = outputTree->Branch("ptsmrw_bveto",&ptsmrw_bveto,"ptsmrw_bveto/D");
	//TBranch *b_etrw_bveto = outputTree->Branch("etrw_bveto",&etrw_bveto,"etrw/D");
	//TBranch *b_etsmrw_bveto = outputTree->Branch("etsmrw_bveto",&etsmrw_bveto,"etsmrw_bveto/D");
	TBranch *b_htrw_bveto = outputTree->Branch("htrw_bveto",&htrw_bveto,"htrw_bveto/D");
	//TBranch *b_njetrw_bveto = outputTree->Branch("njetrw_bveto",&njetrw_bveto,"njetrw_bveto/D");
	//TBranch *b_nbjetrw_bveto = outputTree->Branch("nbjetrw_bveto",&nbjetrw_bveto,"nbjetrw_bveto/D");
	//TBranch *b_ptrw_err_bveto = outputTree->Branch("ptrw_err_bveto",&ptrw_err_bveto,"ptrw_err_bveto/D");
	//TBranch *b_ptsmrw_err_bveto = outputTree->Branch("ptsmrw_err_bveto",&ptsmrw_err_bveto,"ptsmrw_err_bveto/D");
	//TBranch *b_etrw_err_bveto = outputTree->Branch("etrw_err_bveto",&etrw_err_bveto,"etrw_err/D");
	//TBranch *b_etsmrw_err_bveto = outputTree->Branch("etsmrw_err_bveto",&etsmrw_err_bveto,"etsmrw_err_bveto/D");
	//TBranch *b_htrw_err_bveto = outputTree->Branch("htrw_err_bveto",&htrw_err_bveto,"htrw_err_bveto/D");
	//TBranch *b_njetrw_err_bveto = outputTree->Branch("njetrw_err_bveto",&njetrw_err_bveto,"njetrw_err_bveto/D");
	//TBranch *b_nbjetrw_err_bveto = outputTree->Branch("nbjetrw_err_bveto",&nbjetrw_err_bveto,"nbjetrw_err_bveto/D");

	//-----------------------------
	// build TGraph reweighting factors
	//-----------------------------

	TGraph * gr_bveto_ptrw_correction = new TGraph(hist_bveto_ptrw_correction);
	TGraph * gr_bveto_ptsmrw_correction = new TGraph(hist_bveto_ptsmrw_correction);
	TGraph * gr_bveto_etrw_correction = new TGraph(hist_bveto_etrw_correction);
	TGraph * gr_bveto_etsmrw_correction = new TGraph(hist_bveto_etsmrw_correction);
	TGraph * gr_bveto_htrw_correction = new TGraph(hist_bveto_htrw_correction);

	//-----------------------------
	// loop over events
	//-----------------------------

	Long64_t nentries = outputTree->GetEntries();

	for (Long64_t i=0;i<nentries;i++) {

		if (fmod(i,1e5)==0) std::cout << i << " events processed." << std::endl;
		outputTree->GetEntry(i);

		// assess the reweighting factor according to the pT/smeared pT/Njets/HT value
		double Jpsi_et_smear = pow(Jpsi_pt_smear*Jpsi_pt_smear+mll*mll,0.5);
		//double Jpsi_et_smear = Jpsi_pt_smear;
		ptrw_bveto = gr_bveto_ptrw_correction->Eval(Jpsi_pt);
		ptsmrw_bveto = gr_bveto_ptsmrw_correction->Eval(Jpsi_pt_smear);
		etrw_bveto = gr_bveto_etrw_correction->Eval(Jpsi_pt);
		etsmrw_bveto = gr_bveto_etsmrw_correction->Eval(Jpsi_et_smear);
		htrw_bveto = gr_bveto_htrw_correction->Eval(HT);
		njetrw_bveto = hist_bveto_njetrw_correction->GetBinContent(njet+1);
		nbjetrw_bveto = hist_bveto_nbjetrw_correction->GetBinContent(nbjet+1);
		ptrw_err_bveto = hist_bveto_ptrw_correction->GetBinError(pt+1);
		ptsmrw_err_bveto = hist_bveto_ptsmrw_correction->GetBinError(ptsm+1);
		etrw_err_bveto = hist_bveto_etrw_correction->GetBinError(pt+1);
		etsmrw_err_bveto = hist_bveto_etsmrw_correction->GetBinError(ptsm+1);
		htrw_err_bveto = hist_bveto_htrw_correction->GetBinError(ht+1);
		njetrw_err_bveto = hist_bveto_njetrw_correction->GetBinError(njet+1);
		nbjetrw_err_bveto = hist_bveto_nbjetrw_correction->GetBinError(nbjet+1);

		b_ptrw_bveto->Fill();     
		b_ptsmrw_bveto->Fill();     
		//b_etrw_bveto->Fill();     
		//b_etsmrw_bveto->Fill();     
		b_htrw_bveto->Fill();     
		//b_njetrw_bveto->Fill();     
		//b_nbjetrw_bveto->Fill();     
		//b_ptrw_err_bveto->Fill();     
		//b_ptsmrw_err_bveto->Fill();     
		//b_etrw_err_bveto->Fill();     
		//b_etsmrw_err_bveto->Fill();     
		//b_htrw_err_bveto->Fill();     
		//b_njetrw_err_bveto->Fill();     
		//b_nbjetrw_err_bveto->Fill();     

	}

	outputTree->Write();
	hist_bveto_ptrw_correction->Write();
	gr_bveto_ptrw_correction->Write();

	std::cout << "done." << std::endl;
	delete f;


}

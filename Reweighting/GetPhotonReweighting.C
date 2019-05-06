/*
-----------------------------------------------------------------------------------------------
This script takes the outputs from GetBaseLineEvents.C and GetPhotonEvents.C and GetPhotonSmearing.C, and makes photon reweighting factors.
The parameters of the function GetPhotonReweighting(string label, string ch, int isData, int smearing_method, int step) are:
label: takes as an input the dataset year as data15-16, data17 or data18
ch: which dilepton channel you are smearing the photon events to (ee,mm)
isData: 0 (MC) or 1 (data)
step: the latest version of this code is split into two steps (1) reweight with HT, then (2) reweight with Pt_z. First run with step = 1, then step = 2
example of code running command: root -l -b 'GetPhotonReweighting.C+("data15-16","mm",0,0)'
-----------------------------------------------------------------------------------------------
*/

#include "../Settings.C"
#include "../CommonFunctions/CommonLibraries.C"
#include "../CommonFunctions/CommonFunctions.C"
#include "GetSimpleReweightingHistograms.C"

using namespace std;

vector<string> noSampleWeight;
vector<string> noEventWeight;

void GetPhotonReweighting(string periodlabel, string ch, string isData, string smearing_mode, int step) {

    //---------------------------------------------
    // standard 1-d reweighting array 
    //---------------------------------------------

    TH1F* hreweight = GetSimpleReweightingHistograms(periodlabel, ch, smearing_mode, step);
    cout << "Got reweighting histogram hratio with integral " << hreweight->Integral() << endl;

    //---------------------------------------------
    // open file, get Tree and EventCountHist
    //---------------------------------------------

    cout << "Getting smeared photon data" << endl;

    TH1::SetDefaultSumw2();

    string  filename;
    if (isData == "Data") {
        filename = TString(TString(reweighting_path) + "gdata/" + periodlabel + "_merged_processed"  + "_" +TString(ch) + "_" + TString(smearing_mode) + ".root");
        cout << "opening data file" << endl;
    }
    if (isData == "MC") {
        filename = TString(TString(reweighting_path) + "gmc/gmc_" + TString(ch) + "_" + TString(smearing_mode) + ".root");
        cout << "bypassing gdata dir" <<  endl; 
        cout << "opening MC file" << endl;
    }

    TFile*  smeared_file              = new TFile(filename.c_str(),"update");          
    TTree*  outputTree              = (TTree*)smeared_file->Get("BaselineTree");

    cout << endl;
    cout << "Opening file           : " << filename        << endl;
    cout << "Events in ntuple       : " << outputTree->GetEntries() << endl;

    //-----------------------------
    // access existing branch and add new branch
    //-----------------------------

    float gamma_pt = 0.; SetInputBranch(outputTree, "gamma_pt", &gamma_pt);

    Float_t ptreweight = 0.;
    TBranch *b_ptreweight;
    if (step == 1)
        b_ptreweight = outputTree->Branch("ptreweight_step1",&ptreweight,"ptreweight_step1/F");
    else if (step == 2)
        b_ptreweight = outputTree->Branch("ptreweight_step2",&ptreweight,"ptreweight_step2/F");

    //-----------------------------
    // loop over events and fill new branch
    //-----------------------------

    Long64_t nentries = outputTree->GetEntries();

    for (Long64_t i=0; i<nentries; i++) {

        if (fmod(i,1e5)==0) std::cout << i << " events processed." << std::endl;
        outputTree->GetEntry(i);

        float gamma_pt_truncated = gamma_pt;
        if( gamma_pt_truncated < 40   ) gamma_pt_truncated = 40;
        if( gamma_pt_truncated > 1000 ) gamma_pt_truncated = 1000;

        int ptbin = hreweight->FindBin( gamma_pt_truncated );

        ptreweight = hreweight->GetBinContent(ptbin);

        b_ptreweight->Fill();
    }

    outputTree->Write();

    std::cout << "done." << std::endl;
    delete smeared_file;
}

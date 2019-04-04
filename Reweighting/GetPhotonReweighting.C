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
    // Standard 1-d reweighting array 
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
    // access existing branches
    //-----------------------------

    //double totalWeight = 0.; SetInputBranch(outputTree, "totalWeight", &totalWeight);
    //int pt = 0; SetInputBranch(outputTree, "pt", &pt);
    //int ptsm = 0; SetInputBranch(outputTree, "pt_smear", &ptsm);
    //int metsm = 0; SetInputBranch(outputTree, "met_smear", &metsm);
    //int ht = 0; SetInputBranch(outputTree, "ht", &ht);
    //int njet = 0; SetInputBranch(outputTree, "njet", &njet);
    //int nbjet = 0; SetInputBranch(outputTree, "nbjet", &nbjet);
    //Int_t jet_n = 0; SetInputBranch(outputTree, "jet_n", &jet_n);
    //float mll = 0.; SetInputBranch(outputTree, "mll", &mll);
    ////float HT = 0.; SetInputBranch(outputTree, "ht", &HT);
    float gamma_pt = 0.; SetInputBranch(outputTree, "gamma_pt", &gamma_pt);
    //float gamma_pt_smear = 0.; SetInputBranch(outputTree, "Z_pt", &gamma_pt_smear);
    //float gamma_ht = 0.; SetInputBranch(outputTree, "gamma_ht", &gamma_ht);
    //float gamma_ht_smear = 0.; SetInputBranch(outputTree, "HT", &gamma_ht_smear);
    //float MET_smear = 0.; SetInputBranch(outputTree, "MET", &MET_smear);
    //// variables for 2019 RJR analysis
    //float MET_loose = 0.; SetInputBranch(outputTree, "MET_loose", &MET_loose);
    //float MET_tight = 0.; SetInputBranch(outputTree, "MET_tight", &MET_tight);
    //float MET_tighter = 0.; SetInputBranch(outputTree, "MET_tighter", &MET_tighter);
    //float MET_tenacious = 0.; SetInputBranch(outputTree, "MET_tenacious", &MET_tenacious);
    //Int_t trigMatch_2LTrigOR; SetInputBranch(outputTree, "trigMatch_2LTrigOR", &trigMatch_2LTrigOR);
    //Int_t is2Lep2Jet; SetInputBranch(outputTree, "is2Lep2Jet", &is2Lep2Jet);
    //Int_t is2L2JInt; SetInputBranch(outputTree, "is2L2JInt", &is2L2JInt);
    //int nBJet20_MV2c10_FixedCutBEff_77; SetInputBranch(outputTree, "nBJet20_MV2c10_FixedCutBEff_77", &nBJet20_MV2c10_FixedCutBEff_77);
    //float mjj; SetInputBranch(outputTree, "mjj", &mjj);
    //Float_t mll_RJ; SetInputBranch(outputTree, "mll_RJ", &mll_RJ);
    //Float_t R_minH2P_minH3P; SetInputBranch(outputTree, "R_minH2P_minH3P", &R_minH2P_minH3P);
    //Float_t RPT_HT5PP; SetInputBranch(outputTree, "RPT_HT5PP", &RPT_HT5PP);
    //Float_t dphiVP; SetInputBranch(outputTree, "dphiVP", &dphiVP);
    //Float_t H2PP; SetInputBranch(outputTree, "H2PP", &H2PP);
    //Float_t H5PP; SetInputBranch(outputTree, "H5PP", &H5PP);
    //int nJet20; SetInputBranch(outputTree, "nJet20", &nJet20);
    //Float_t minDphi; SetInputBranch(outputTree, "minDphi", &minDphi);
    //Float_t MZ; SetInputBranch(outputTree, "MZ", &MZ);
    //Int_t NjS; SetInputBranch(outputTree, "NjS", &NjS);
    //Int_t NjISR; SetInputBranch(outputTree, "NjISR", &NjISR);
    //Float_t dphiISRI; SetInputBranch(outputTree, "dphiISRI", &dphiISRI);
    //Float_t RISR; SetInputBranch(outputTree, "RISR", &RISR);
    //Float_t PTISR; SetInputBranch(outputTree, "PTISR", &PTISR);
    //Float_t PTI; SetInputBranch(outputTree, "PTI", &PTI);
    //Float_t PTCM; SetInputBranch(outputTree, "PTCM", &PTCM);
    //Float_t MJ; SetInputBranch(outputTree, "MJ", &MJ);
    //Int_t is3Lep3Jet; SetInputBranch(outputTree, "is3Lep3Jet", &is3Lep3Jet);
    //Int_t is4Lep3Jet; SetInputBranch(outputTree, "is4Lep3Jet", &is4Lep3Jet);
    //Int_t lept1sign_VR; SetInputBranch(outputTree, "lept1sign_VR", &lept1sign_VR);
    //Int_t lept2sign_VR; SetInputBranch(outputTree, "lept2sign_VR", &lept2sign_VR);
    //Float_t lept1Pt_VR; SetInputBranch(outputTree, "lept1Pt_VR", &lept1Pt_VR);
    //Float_t lept2Pt_VR; SetInputBranch(outputTree, "lept2Pt_VR", &lept2Pt_VR);
    //Float_t MZ_VR; SetInputBranch(outputTree, "MZ_VR", &MZ_VR);
    //Float_t MJ_VR; SetInputBranch(outputTree, "MJ_VR", &MJ_VR);
    //Float_t RISR_VR; SetInputBranch(outputTree, "RISR_VR", &RISR_VR);
    //Float_t PTISR_VR; SetInputBranch(outputTree, "PTISR_VR", &PTISR_VR);
    //Float_t PTI_VR; SetInputBranch(outputTree, "PTI_VR", &PTI_VR);
    //Float_t PTCM_VR; SetInputBranch(outputTree, "PTCM_VR", &PTCM_VR);
    //Float_t dphiISRI_VR; SetInputBranch(outputTree, "dphiISRI_VR", &dphiISRI_VR);

    //std::vector<int>* lepFlavor = new std::vector<int>(10); SetInputBranch(outputTree, "lepFlavor", &lepFlavor);
    //std::vector<int>* lepCharge = new std::vector<int>(10); SetInputBranch(outputTree, "lepCharge", &lepCharge);
    //std::vector<float>* jet_pT = new std::vector<float>(10); SetInputBranch(outputTree, "jet_pT", &jet_pT);
    //std::vector<float>* lep_pT = new std::vector<float>(10); SetInputBranch(outputTree, "lep_pT", &lep_pT);

    //if (isData == "MC") {
        //float gamma_dR = 0.; SetInputBranch(outputTree, "gamma_dR", &gamma_dR);
    //}

    //-----------------------------
    // add new branches
    //-----------------------------

    Float_t ptreweight = 0.;
    TBranch *b_ptreweight;
    if (step == 1)
        b_ptreweight = outputTree->Branch("ptreweight3",&ptreweight,"ptreweight3/F");
    else if (step == 2)
        b_ptreweight = outputTree->Branch("ptreweight5",&ptreweight,"ptreweight5/F");

    //-----------------------------
    // loop over events
    //-----------------------------

    Long64_t nentries = outputTree->GetEntries();

    for (Long64_t i=0;i<nentries;i++) {

        if (fmod(i,1e5)==0) std::cout << i << " events processed." << std::endl;
        outputTree->GetEntry(i);

        float gamma_pt_truncated = gamma_pt;
        if( gamma_pt_truncated < 40   ) gamma_pt_truncated = 41;
        if( gamma_pt_truncated > 1000 ) gamma_pt_truncated = 999;

        int ptbin = hreweight->FindBin( gamma_pt_truncated );

        ptreweight = hreweight->GetBinContent(ptbin);

        b_ptreweight->Fill();
    }

    outputTree->Write();

    std::cout << "done." << std::endl;
    delete smeared_file;
}

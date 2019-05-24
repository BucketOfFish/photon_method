#include "../Common/Settings.C"
#include "../Common/CommonLibraries.C"
#include "../Common/CommonFunctions.C"
#include "../Common/CommonCuts.C"

using namespace std;

TH1F* GetSimpleReweightingHistograms(string period, string channel, string isData, string smearing_mode, int step ){

    cout << "Making reweighting histograms for period and year " << period << " " << channel << endl;
    gStyle->SetOptStat(0);

    //--- open files and create TChains
    string mc_folder = "";
    if (TString(period).Contains("data15-16")) mc_folder = "ZMC16a/";
    else if (TString(period).Contains("data17")) mc_folder = "ZMC16cd/";
    else if (TString(period).Contains("data18")) mc_folder = "ZMC16cd/";

    string data_filename = ntuple_path + "zdata/" + period + "_merged_processed.root";
    string tt_filename = ntuple_path + mc_folder + "ttbar_merged_processed.root";
    string vv_filename = ntuple_path + mc_folder + "diboson_merged_processed.root";
    string zjets_filename = ntuple_path + mc_folder + "Zjets_merged_processed.root";
    string photon_filename;
    if (isData == "Data") photon_filename = reweighting_path + "gdata/" + period + "_merged_processed" + "_" + channel + "_" + smearing_mode + ".root"; //Vg subtracted 
    else photon_filename = reweighting_path + "gmc/" "gmc" + "_" + channel + "_" + smearing_mode + ".root"; //Vg subtracted 

    cout << "Opening data file    " << data_filename << endl;
    cout << "Opening ttbar file   " << tt_filename << endl;
    cout << "Opening diboson file " << vv_filename << endl;
    cout << "Opening Z+jets file  " << zjets_filename << endl;
    cout << "Opening photon file  " << photon_filename << endl;

    TChain* tch_data = new TChain("BaselineTree"); tch_data->Add(data_filename.c_str());
    TChain* tch_tt = new TChain("BaselineTree"); tch_tt->Add(tt_filename.c_str());
    TChain* tch_vv = new TChain("BaselineTree"); tch_vv->Add(vv_filename.c_str());
    TChain* tch_zjets = new TChain("BaselineTree"); tch_zjets->Add(zjets_filename.c_str());
    TChain* tch_photon = new TChain("BaselineTree"); tch_photon->Add(photon_filename.c_str());

    cout << "data entries         " << tch_data->GetEntries() << endl;
    cout << "ttbar entries        " << tch_tt->GetEntries() << endl;
    cout << "diboson entries      " << tch_vv->GetEntries() << endl;
    cout << "Z+jets entries       " << tch_zjets->GetEntries() << endl;
    cout << "photon entries       " << tch_photon->GetEntries() << endl;

    //--- modify event selections and weights
    if (TString(channel).EqualTo("ee")) cuts::Zselection += cuts::ee;
    else if (TString(channel).EqualTo("mm")) cuts::Zselection += cuts::mm;
    else {
        cout << "Unrecognized channel! quitting   " << channel << endl;
        exit(0);
    }

    if( TString(period).EqualTo("data17")    ){
        TCut RunRange = TCut("RunNumber < 348000");  
        cout << "Data17! adding cut " << RunRange.GetTitle() << endl;
        cuts::Zselection *= RunRange;
    }

    cout << "Z selection          " << cuts::Zselection.GetTitle() << endl;
    cout << "Z weight             " << cuts::Zweight.GetTitle() << endl;
    cout << "g selection          " << cuts::gselection.GetTitle() << endl;
    cout << "g weight             " << cuts::weight_g.GetTitle() << endl;

    //--- fill reweighting histograms
    TH1F* hdata  = new TH1F("hdata", "", bins::nptbins, bins::ptbins);
    TH1F* htt    = new TH1F("htt", "", bins::nptbins, bins::ptbins);
    TH1F* hvv    = new TH1F("hvv", "", bins::nptbins, bins::ptbins);
    TH1F* hz     = new TH1F("hz", "", bins::nptbins, bins::ptbins);
    TH1F* histoG = new TH1F("histoG", "", bins::nptbins, bins::ptbins);    

    // step 1: HT
    if (step == 1) {
        tch_data->Draw("min(HT,999)>>hdata", cuts::Zselection, "goff");
        tch_tt->Draw("min(HT,999)>>htt", cuts::Zselection*cuts::Zweight, "goff");
        tch_vv->Draw("min(HT,999)>>hvv", cuts::Zselection*cuts::Zweight, "goff");
        tch_zjets->Draw("min(HT,999)>>hz", cuts::Zselection*cuts::Zweight, "goff");
        tch_photon->Draw("min(HT,999)>>histoG", cuts::gselection*cuts::weight_g, "goff");
    }

    // step 2: Z_pt
    else if (step == 2) {
        TCut g_rw("ptreweight_step1"); // from step 1
        tch_data->Draw("min(Ptll,999)>>hdata", cuts::Zselection, "goff");
        tch_tt->Draw("min(Ptll,999)>>htt", cuts::Zselection*cuts::Zweight, "goff");
        tch_vv->Draw("min(Ptll,999)>>hvv", cuts::Zselection*cuts::Zweight, "goff");
        tch_zjets->Draw("min(Ptll,999)>>hz", cuts::Zselection*cuts::Zweight, "goff");
        tch_photon->Draw("min(Ptll,999)>>histoG", cuts::gselection*cuts::weight_g*g_rw, "goff");
    }

    cout << "data integral        " << hdata->Integral() << endl;
    cout << "ttbar integral       " << htt->Integral() << endl;
    cout << "diboson integral     " << hvv->Integral() << endl;
    cout << "Z+jets integral      " << hz->Integral() << endl;
    cout << "photon integral      " << histoG->Integral() << endl;

    //--- calculate reweighting ratios
    TH1F* histoZ = (TH1F*) hdata->Clone("histoZ");
    histoZ->Add(htt, -1.0);
    histoZ->Add(hvv, -1.0);

    TH1F* hratio = (TH1F*) histoZ->Clone("hratio");
    hratio->Divide(histoG);

    cout << "histoG->Integral() " << histoG->Integral() << endl;
    cout << "histoZ->Integral() " << histoZ->Integral() << endl;
    cout << "hratio->Integral() " << hratio->Integral() << endl;

    return hratio;
}

void GetPhotonReweighting(string periodlabel, string ch, string isData, string smearing_mode, int step) {

    //---------------------------------------------
    // 1-d reweighting histogram 
    //---------------------------------------------

    TH1F* hreweight = GetSimpleReweightingHistograms(periodlabel, ch, isData, smearing_mode, step);
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

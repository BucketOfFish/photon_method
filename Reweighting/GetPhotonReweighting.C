#include "../Common/Settings.C"
#include "../Common/CommonLibraries.C"
#include "../Common/CommonFunctions.C"

using namespace std;

TH1F* GetSimpleReweightingHistograms(string period, string channel, string data_or_mc, string photon_filename, string smearing_mode, string reweight_var){

    cout << "Making reweighting histograms for period and year " << period << " " << channel << endl;
    gStyle->SetOptStat(0);

    //--- open files and create TChains
    string mc_period = "";
    if (TString(period).Contains("data15-16")) mc_period = "mc16a";
    else if (TString(period).Contains("data17")) mc_period = "mc16cd";
    else if (TString(period).Contains("data18")) mc_period = "mc16e";

    string data_filename = ntuple_path + "bkg_data/" + period + "_bkg.root";
    string tt_filename = ntuple_path + "bkg_mc/" + mc_period + "_ttbar.root";
    string vv_filename = ntuple_path + "bkg_mc/" + mc_period + "_diboson.root";
    string zjets_filename = ntuple_path + "bkg_mc/" + mc_period + "_Zjets.root";

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
    if (TString(channel).EqualTo("ee")) cuts::bkg_baseline += cuts::ee;
    else if (TString(channel).EqualTo("mm")) cuts::bkg_baseline += cuts::mm;
    else {
        cout << "Unrecognized channel! quitting   " << channel << endl;
        exit(0);
    }

    if (TString(period).EqualTo("data17")){
        TCut RunRange = TCut("RunNumber < 348000");  
        cout << "Data17! adding cut " << RunRange.GetTitle() << endl;
        cuts::bkg_baseline *= RunRange;
    }

    cout << "Z selection          " << cuts::bkg_baseline.GetTitle() << endl;
    cout << "Z weight             " << cuts::bkg_weight.GetTitle() << endl;
    cout << "g selection          " << cuts::photon_baseline.GetTitle() << endl;
    cout << "g weight             " << cuts::photon_weight.GetTitle() << endl;

    //--- fill reweighting histograms
    TH1F* hdata  = new TH1F("hdata", "", bins::n_reweighting_bins, bins::reweighting_bins);
    TH1F* htt    = new TH1F("htt", "", bins::n_reweighting_bins, bins::reweighting_bins);
    TH1F* hvv    = new TH1F("hvv", "", bins::n_reweighting_bins, bins::reweighting_bins);
    TH1F* hz     = new TH1F("hz", "", bins::n_reweighting_bins, bins::reweighting_bins);
    TH1F* histoG = new TH1F("histoG", "", bins::n_reweighting_bins, bins::reweighting_bins);    

    //--- reweighting variable histograms
    tch_data->Draw(("min("+reweight_var+",999)>>hdata").c_str(), cuts::bkg_baseline, "goff");
    tch_tt->Draw(("min("+reweight_var+",999)>>htt").c_str(), cuts::bkg_baseline*cuts::bkg_weight, "goff");
    tch_vv->Draw(("min("+reweight_var+",999)>>hvv").c_str(), cuts::bkg_baseline*cuts::bkg_weight, "goff");
    tch_zjets->Draw(("min("+reweight_var+",999)>>hz").c_str(), cuts::bkg_baseline*cuts::bkg_weight, "goff");
    tch_photon->Draw(("min("+reweight_var+",999)>>histoG").c_str(), cuts::photon_baseline*cuts::photon_weight, "goff");

    cout << "data integral        " << hdata->Integral() << endl;
    cout << "ttbar integral       " << htt->Integral() << endl;
    cout << "diboson integral     " << hvv->Integral() << endl;
    cout << "Z+jets integral      " << hz->Integral() << endl;
    cout << "photon integral      " << histoG->Integral() << endl;

    //--- calculate reweighting ratios
    TH1F* histoZ = (TH1F*) hdata->Clone("histoZ");
    histoZ->Add(htt, -1.0);
    histoZ->Add(hvv, -1.0);

    if (data_or_mc == "MC")
        histoZ = (TH1F*) hz->Clone("histoZ");

    TH1F* hratio = (TH1F*) histoZ->Clone("hratio");
    hratio->Divide(histoG);

    cout << "histoG->Integral() " << histoG->Integral() << endl;
    cout << "histoZ->Integral() " << histoZ->Integral() << endl;
    cout << "hratio->Integral() " << hratio->Integral() << endl;

    return hratio;
}

void GetPhotonReweighting(string period, string channel, string data_or_mc, string smearing_mode, string reweight_var) {

    //---------------------------------------------
    // open file, get Tree and EventCountHist
    //---------------------------------------------

    cout << "Getting smeared photon data" << endl;

    TH1::SetDefaultSumw2();

    string mc_period = "";
    if (TString(period).Contains("data15-16")) mc_period = "mc16a";
    else if (TString(period).Contains("data17")) mc_period = "mc16cd";
    else if (TString(period).Contains("data18")) mc_period = "mc16e";

    string photon_filename;
    if (data_or_mc == "Data") {
        photon_filename = TString(reweighting_path+"g_data/"+period+"_photon_"+channel+"_"+smearing_mode+".root");
        cout << "opening data file" << endl;
    }
    if (data_or_mc == "MC") {
        photon_filename = TString(reweighting_path+"g_mc/"+mc_period+"_SinglePhoton222_"+channel+"_"+smearing_mode+".root");
        cout << "opening MC file" << endl;
    }

    TFile* smeared_file = new TFile(photon_filename.c_str(),"update");          
    TTree* outputTree = (TTree*)smeared_file->Get("BaselineTree");

    cout << endl;
    cout << "Opening file           : " << photon_filename << endl;
    cout << "Events in ntuple       : " << outputTree->GetEntries() << endl;

    //---------------------------------------------
    // 1-d reweighting histogram 
    //---------------------------------------------

    TH1F* h_reweight = GetSimpleReweightingHistograms(period, channel, data_or_mc, photon_filename, smearing_mode, reweight_var);
    cout << "Got reweighting histogram hratio with integral " << h_reweight->Integral() << endl;

    //-----------------------------
    // loop over events and fill new branch
    //-----------------------------

    float gamma_var = 0.; SetInputBranch(outputTree, reweight_var, &gamma_var);

    Float_t reweight = 0.;
    TBranch *b_reweight;
    b_reweight = outputTree->Branch(("reweight_"+reweight_var).c_str(), &reweight, ("reweight_"+reweight_var+"/F").c_str());

    Long64_t nentries = outputTree->GetEntries();
    for (Long64_t i=0; i<nentries; i++) {

        if (fmod(i,1e5)==0) std::cout << i << " events processed." << std::endl;
        outputTree->GetEntry(i);

        float gamma_var_truncated = gamma_var;
        if(gamma_var_truncated < bins::reweighting_bins[0]) gamma_var_truncated = bins::reweighting_bins[0];
        if(gamma_var_truncated > bins::reweighting_bins[bins::n_reweighting_bins]) gamma_var_truncated = bins::reweighting_bins[bins::n_reweighting_bins];
        int ptbin = h_reweight->FindBin( gamma_var_truncated );

        int var_bin = h_reweight->FindBin(gamma_var);
        reweight = h_reweight->GetBinContent(var_bin);
        b_reweight->Fill();
    }

    outputTree->Write();

    std::cout << "done." << std::endl;
    delete smeared_file;
}

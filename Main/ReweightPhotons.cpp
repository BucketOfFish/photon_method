#include "Settings.cpp"

using namespace std;

TH1F* GetSimpleReweightingHistograms(ReweightingOptions options) {
    gStyle->SetOptStat(0);

    //--- open files and create TChains
    string data_filename = options.reduction_folder + options.period + "_data_bkg.root";
    string tt_filename = options.reduction_folder + options.mc_period + "_ttbar.root";
    string vv_filename = options.reduction_folder + options.mc_period + "_diboson.root";
    string zjets_filename = options.reduction_folder + options.mc_period + "_Zjets.root";

    TChain *tch_data, *tch_tt, *tch_vv, *tch_zjets, *tch_photon;

    if (options.is_data) {
        cout << "Opening data file      : " << data_filename << endl;
        tch_data = new TChain(options.out_tree_name.c_str()); tch_data->Add(data_filename.c_str());
        cout << "data entries           : " << tch_data->GetEntries() << endl;
        cout << "Opening ttbar file     : " << tt_filename << endl;
        tch_tt = new TChain(options.out_tree_name.c_str()); tch_tt->Add(tt_filename.c_str());
        cout << "ttbar entries          : " << tch_tt->GetEntries() << endl;
        cout << "Opening diboson file   : " << vv_filename << endl;
        tch_vv = new TChain(options.out_tree_name.c_str()); tch_vv->Add(vv_filename.c_str());
        cout << "diboson entries        : " << tch_vv->GetEntries() << endl;
    }
    else {
        cout << "Opening Z+jets file    : " << zjets_filename << endl;
        tch_zjets = new TChain(options.out_tree_name.c_str()); tch_zjets->Add(zjets_filename.c_str());
        cout << "Z+jets entries         : " << tch_zjets->GetEntries() << endl;
    }
    cout << "Opening photon file    : " << options.in_file_name << endl;
    tch_photon = new TChain(options.out_tree_name.c_str()); tch_photon->Add(options.in_file_name.c_str());
    cout << "photon entries         : " << tch_photon->GetEntries() << endl;
    cout << endl;

    //--- modify event selections and weights
    TCut reweight_region = cuts::reweight_region;
    if (TString(options.channel).EqualTo("ee")) reweight_region += cuts::ee;
    else if (TString(options.channel).EqualTo("mm")) reweight_region += cuts::mm;
    else failTest("Unrecognized channel " + options.channel);

    cout << "bkg selection          : " << reweight_region.GetTitle() << endl;
    cout << "bkg weight             : " << cuts::bkg_weight.GetTitle() << endl;
    cout << "photon selection       : " << reweight_region.GetTitle() << endl;
    cout << "photon weight          : " << cuts::photon_weight.GetTitle() << endl;
    cout << endl;

    //--- set reweighting variable
    int n_reweighting_bins = bins::n_pt_bins;
    double *reweighting_bins = bins::pt_bins;

    //--- fill reweighting histograms
    TH1F* hdata  = new TH1F("hdata", "", n_reweighting_bins, reweighting_bins);
    TH1F* htt    = new TH1F("htt", "", n_reweighting_bins, reweighting_bins);
    TH1F* hvv    = new TH1F("hvv", "", n_reweighting_bins, reweighting_bins);
    TH1F* hz     = new TH1F("hz", "", n_reweighting_bins, reweighting_bins);
    TH1F* histoG = new TH1F("histoG", "", n_reweighting_bins, reweighting_bins);    

    if (options.is_data) {
        tch_data->Draw((options.reweight_var+">>hdata").c_str(), reweight_region, "goff");
        tch_tt->Draw((options.reweight_var+">>htt").c_str(), reweight_region*cuts::bkg_weight, "goff");
        tch_vv->Draw((options.reweight_var+">>hvv").c_str(), reweight_region*cuts::bkg_weight, "goff");
        cout << "data integral          : " << hdata->Integral() << endl;
        cout << "ttbar integral         : " << htt->Integral() << endl;
        cout << "diboson integral       : " << hvv->Integral() << endl;
    }
    else {
        tch_zjets->Draw((options.reweight_var+">>hz").c_str(), reweight_region*cuts::bkg_weight, "goff");
        cout << "Z+jets integral        : " << hz->Integral() << endl;
    }

    tch_photon->Draw((options.reweight_var+">>histoG").c_str(), reweight_region*cuts::photon_weight, "goff");
    cout << "photon integral      " << histoG->Integral() << endl;
    cout << endl;

    //--- calculate reweighting ratios
    TH1F* histoZ;
    if (options.is_data) {
        histoZ = (TH1F*) hdata->Clone("histoZ");
        histoZ->Add(htt, -1.0);
        histoZ->Add(hvv, -1.0);
    }
    else histoZ = (TH1F*) hz->Clone("histoZ");

    TH1F* hratio = (TH1F*) histoZ->Clone("hratio");
    hratio->Divide(histoG);

    cout << "photon integral        : " << histoG->Integral() << endl;
    cout << "bkg integral           : " << histoZ->Integral() << endl;
    cout << "scaling factor         : " << hratio->Integral() << endl;
    cout << endl;

    return hratio;
}

void ReweightPhotons(ReweightingOptions options) {

    //---------------------------------------------
    // open file, get Tree and EventCountHist
    //---------------------------------------------
    TH1::SetDefaultSumw2();

    cout << BOLD(PBLU("Reweighting photon events")) << endl;
    cout << endl;

    cout << "Period                 : " << options.period << endl;
    cout << "Channel                : " << options.channel << endl;
    cout << "Opening file           : " << options.in_file_name << endl;
    cout << "Reading tree           : " << options.out_tree_name << endl;

    TFile* smeared_file = new TFile(options.in_file_name.c_str(),"update");          
    TTree* outputTree = (TTree*)smeared_file->Get(options.out_tree_name.c_str());
    cout << "Events in ntuple       : " << outputTree->GetEntries() << endl;
    cout << endl;

    //---------------------------------------------
    // 1-d reweighting histogram 
    //---------------------------------------------
    TH1F* h_reweight = GetSimpleReweightingHistograms(options);

    //-----------------------------
    // loop over events and fill new branch
    //-----------------------------
    float gamma_var = 0.; SetInputBranch(outputTree, options.reweight_var, &gamma_var);

    Float_t reweight = 0.;
    TBranch *b_reweight;
    b_reweight = outputTree->Branch(("reweight_"+options.reweight_var).c_str(), &reweight, ("reweight_"+options.reweight_var+"/F").c_str());

    Long64_t nentries = outputTree->GetEntries();
    for (Long64_t i=0; i<nentries; i++) {

        if (fmod(i,1e5)==0) cout << i << " events processed.\r" << flush;
        outputTree->GetEntry(i);

        //float gamma_var_truncated = gamma_var;
        //if(gamma_var_truncated < reweighting_bins[0]) gamma_var_truncated = reweighting_bins[0];
        //if(gamma_var_truncated > reweighting_bins[n_reweighting_bins]) gamma_var_truncated = reweighting_bins[n_reweighting_bins];
        int var_bin = h_reweight->FindBin(gamma_var);
        reweight = h_reweight->GetBinContent(var_bin);

        b_reweight->Fill();
    }
    cout << endl;

    outputTree->Write();

    cout << PBLU("Finished reweighting") << endl;
    delete smeared_file;
}

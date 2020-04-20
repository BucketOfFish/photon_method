#include "Settings.cpp"

using namespace std;

map<string, TH1F*> GetSimpleReweightingHistograms(ReweightingOptions options) {
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

    //--- and calculate reweighting ratios
    TH1F* histoZ;
    map<string, TH1F*> hratios;

    for (auto reweight_var : options.reweight_vars) {
        if (options.is_data) {
            tch_data->Draw((reweight_var+">>hdata").c_str(), reweight_region, "goff");
            tch_tt->Draw((reweight_var+">>htt").c_str(), reweight_region*cuts::bkg_weight, "goff");
            tch_vv->Draw((reweight_var+">>hvv").c_str(), reweight_region*cuts::bkg_weight, "goff");
            cout << "data integral          : " << hdata->Integral(0, n_reweighting_bins+1) << endl;
            cout << "ttbar integral         : " << htt->Integral(0, n_reweighting_bins+1) << endl;
            cout << "diboson integral       : " << hvv->Integral(0, n_reweighting_bins+1) << endl;
        }
        else {
            tch_zjets->Draw((reweight_var+">>hz").c_str(), reweight_region*cuts::bkg_weight, "goff");
            cout << "Z+jets integral        : " << hz->Integral(0, n_reweighting_bins+1) << endl;
        }

        tch_photon->Draw((reweight_var+">>histoG").c_str(), reweight_region*cuts::photon_weight, "goff");
        cout << "photon integral      " << histoG->Integral(0, n_reweighting_bins+1) << endl;
        cout << endl;

        if (options.is_data) {
            histoZ = (TH1F*) hdata->Clone("histoZ");
            histoZ->Add(htt, -1.0);
            histoZ->Add(hvv, -1.0);
        }
        else histoZ = (TH1F*) hz->Clone("histoZ");

        hratios[reweight_var] = (TH1F*) histoZ->Clone(("hratio_" + reweight_var).c_str());
        hratios[reweight_var]->Divide(histoG);

        cout << "photon integral        : " << histoG->Integral(0, n_reweighting_bins+1) << endl;
        cout << "bkg integral           : " << histoZ->Integral(0, n_reweighting_bins+1) << endl;
        cout << "scaling factor         : " << hratios[reweight_var]->Integral(0, n_reweighting_bins+1) << endl;
        cout << endl;
    }

    return hratios;
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
    map<string, TH1F*> reweight_hists = GetSimpleReweightingHistograms(options);

    //-----------------------------
    // loop over events and fill new branch
    //-----------------------------
    map<string, float> var_vals;
    map<string, Float_t> reweight_vals;
    map<string, TBranch*> reweight_branches;

    for (auto reweight_var : options.reweight_vars) {
        var_vals[reweight_var] = 0.0;
        reweight_vals[reweight_var] = 0.0;
        SetInputBranch(outputTree, reweight_var, &var_vals[reweight_var]);
        reweight_branches[reweight_var] = outputTree->Branch(("reweight_"+reweight_var).c_str(),
            &reweight_vals[reweight_var], ("reweight_"+reweight_var+"/F").c_str());
    }

    Long64_t nentries = outputTree->GetEntries();
    for (Long64_t i=0; i<nentries; i++) {
        if (fmod(i,1e5)==0) cout << i << " events processed.\r" << flush;
        outputTree->GetEntry(i);

        //float gamma_var_truncated = var_vals;
        //if(gamma_var_truncated < reweighting_bins[0]) gamma_var_truncated = reweighting_bins[0];
        //if(gamma_var_truncated > reweighting_bins[n_reweighting_bins]) gamma_var_truncated = reweighting_bins[n_reweighting_bins];
        for (auto reweight_var : options.reweight_vars) {
            int var_bin = reweight_hists[reweight_var]->FindBin(var_vals[reweight_var]);
            reweight_vals[reweight_var] = reweight_hists[reweight_var]->GetBinContent(var_bin);
            reweight_branches[reweight_var]->Fill();
        }
    }
    cout << endl;

    outputTree->Write();

    cout << PBLU("Finished reweighting") << endl;
    delete smeared_file;
}

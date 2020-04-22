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

    //--- store reweighting histograms
    TH1F* histoZ;
    map<string, TH1F*> hratios;

    for (auto reweight_var : options.reweight_vars) {
        //--- split variable by ":"
        std::stringstream unsplit_vars(reweight_var);
        std::string segment;
        std::vector<std::string> split_vars;
        while(std::getline(unsplit_vars, segment, ':')) {
           split_vars.push_back(segment);
        }

        //--- single variable reweighting
        reweight_var = split_vars[0];

        //--- set reweighting properties based on variable
        int n_reweighting_bins = bins::n_reweighting_bins.at(reweight_var);
        vector<double> reweighting_bins = bins::reweighting_bins.at(reweight_var);

        TH1F* hdata  = new TH1F("hdata", "", n_reweighting_bins, &reweighting_bins[0]);
        TH1F* htt    = new TH1F("htt", "", n_reweighting_bins, &reweighting_bins[0]);
        TH1F* hvv    = new TH1F("hvv", "", n_reweighting_bins, &reweighting_bins[0]);
        TH1F* hz     = new TH1F("hz", "", n_reweighting_bins, &reweighting_bins[0]);
        TH1F* histoG = new TH1F("histoG", "", n_reweighting_bins, &reweighting_bins[0]);    

        //--- fill reweighting histograms
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

        //--- get ratio
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

        //--- delete pointers
        delete hdata, htt, hvv, hz, histoG;
    }

    return hratios;
}

map<string, TH2F*> GetSimple2DReweightingHistograms(ReweightingOptions options) {
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

    //--- store reweighting histograms
    TH2F* histoZ;
    map<string, TH2F*> hratios;

    for (auto reweight_var : options.reweight_vars) {
        //--- split variable by ":"
        std::stringstream unsplit_vars(reweight_var);
        std::string segment;
        std::vector<std::string> split_vars;
        while(std::getline(unsplit_vars, segment, ':')) {
           split_vars.push_back(segment);
        }

        //--- set reweighting properties based on variables
        vector<int> n_reweighting_bins;
        vector<vector<double>> reweighting_bins;
        int var_n = 0;
        for (auto reweight_var : split_vars) {
            n_reweighting_bins[var_n] = bins::n_reweighting_bins.at(reweight_var);
            reweighting_bins[var_n] = bins::reweighting_bins.at(reweight_var);
            var_n++;
        }

        TH2F* hdata  = new TH2F("hdata", "", n_reweighting_bins[0], &(reweighting_bins[0])[0],
                                n_reweighting_bins[1], &(reweighting_bins[1])[0]);
        TH2F* htt    = new TH2F("htt", "", n_reweighting_bins[0], &(reweighting_bins[0])[0],
                                n_reweighting_bins[1], &(reweighting_bins[1])[0]);
        TH2F* hvv    = new TH2F("hvv", "", n_reweighting_bins[0], &(reweighting_bins[0])[0],
                                n_reweighting_bins[1], &(reweighting_bins[1])[0]);
        TH2F* hz     = new TH2F("hz", "", n_reweighting_bins[0], &(reweighting_bins[0])[0],
                                n_reweighting_bins[1], &(reweighting_bins[1])[0]);
        TH2F* histoG = new TH2F("histoG", "", n_reweighting_bins[0], &(reweighting_bins[0])[0],
                                n_reweighting_bins[1], &(reweighting_bins[1])[0]);

        //--- fill reweighting histograms
        vector<tuple<string, TChain*, TH2F*>> histogram_fills;
        if (options.is_data) {
            histogram_fills.push_back(make_tuple("data", tch_data, hdata));
            histogram_fills.push_back(make_tuple("ttbar", tch_tt, htt));
            histogram_fills.push_back(make_tuple("diboson", tch_vv, hvv));
        }
        else {
            histogram_fills.push_back(make_tuple("Z+jets", tch_zjets, hz));
        }
        histogram_fills.push_back(make_tuple("photon", tch_photon, histoG));

        for (auto histogram_fill : histogram_fills) {
            string hist_name = get<0>(histogram_fill);
            TChain *tchain = get<1>(histogram_fill);
            TH2F *hist = get<2>(histogram_fill);
            Long64_t nentries = tchain->GetEntries();
            vector<float> var_val = {0, 0};
            float weight;
            tchain->SetBranchAddress(split_vars[0].c_str(), &var_val[0]);
            tchain->SetBranchAddress(split_vars[1].c_str(), &var_val[1]);
            tchain->SetBranchAddress(cuts::bkg_weight, &weight);
            for (Long64_t i=0; i<nentries; i++) {
                tchain->GetEntry(i);
                hist->Fill(var_val[0], var_val[1], weight);
            }
            cout << padString(hist_name + " integral") << ": " << hist->Integral() << endl;
        }
        cout << endl;

        //--- get ratio
        if (options.is_data) {
            histoZ = (TH2F*) hdata->Clone("histoZ");
            histoZ->Add(htt, -1.0);
            histoZ->Add(hvv, -1.0);
        }
        else histoZ = (TH2F*) hz->Clone("histoZ");

        hratios[reweight_var] = (TH2F*) histoZ->Clone(("hratio_" + reweight_var).c_str());
        hratios[reweight_var]->Divide(histoG);

        //--- delete pointers
        delete hdata, htt, hvv, hz, histoG;
    }

    return hratios;
}

void ReweightSample(ReweightingOptions options) {
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
    // 1D reweighting histogram 
    //---------------------------------------------
    map<string, TH1F*> reweight_hists = GetSimpleReweightingHistograms(options);

    //---------------------------------------------
    // 2D reweighting histogram 
    //---------------------------------------------
    //map<string, TH2F*> reweight_hists = GetSimple2DReweightingHistograms(options);

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

void ReweightPhotons(ReweightingOptions options) {
    if (options.is_data) {
        //options.in_file_name = options.reweighting_folder + options.data_period + "_data_photon_" + options.channel + ".root";
        //options.out_file_name = options.in_file_name;
        //ReweightSample(options);

        options.in_file_name = options.reweighting_folder + options.mc_period + "_Vgamma_" + options.channel + ".root";
        options.out_file_name = options.in_file_name;
        ReweightSample(options);
    }
    else {
        options.in_file_name = options.reweighting_folder + options.mc_period + "_SinglePhoton222_" + options.channel + ".root";
        options.out_file_name = options.in_file_name;
        ReweightSample(options);
    }
}

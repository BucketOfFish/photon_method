#include "Settings.cpp"

using namespace std;

map<string, TChain*> getTChains(ReweightingOptions options) {
    //--- open files and create TChains
    map<string, string> filenames;
    filenames["data"] = options.reduction_folder + options.period + "_data_bkg.root";
    filenames["tt"] = options.reduction_folder + options.mc_period + "_ttbar.root";
    filenames["vv"] = options.reduction_folder + options.mc_period + "_diboson.root";
    filenames["zjets"] = options.reduction_folder + options.mc_period + "_Zjets.root";
    filenames["photon"] = options.in_file_name;

    map<string, TChain*> tchains;

    if (options.is_data) vector<string> processes = {"data", "tt", "vv", "photon"};
    else vector<string> processes = {"zjets", "photon"};
    for (auto process : processes) {
        cout << padString("Opening " + process + "file") << ": " << filenames[process] << endl;
        tchains[process] = new TChain(options.out_tree_name.c_str());
        tchains[process]->Add(filenames[process].c_str());
        cout << padString(process + " entries") << ": " << tchains[process]->GetEntries() << endl;
    }
    cout << endl;

    return tchains;
}

TCut getReweightRegion(ReweightingOptions options) {
    TCut reweight_region = cuts::reweight_region;
    if (TString(options.channel).EqualTo("ee")) reweight_region += cuts::ee;
    else if (TString(options.channel).EqualTo("mm")) reweight_region += cuts::mm;
    else failTest("Unrecognized channel " + options.channel);

    cout << "bkg selection          : " << reweight_region.GetTitle() << endl;
    cout << "bkg weight             : " << cuts::bkg_weight.GetTitle() << endl;
    cout << "photon selection       : " << reweight_region.GetTitle() << endl;
    cout << "photon weight          : " << cuts::photon_weight.GetTitle() << endl;
    cout << endl;

    return reweight_region;
}

map<string, TH1F*> GetSimpleReweightingHistograms(ReweightingOptions options) {
    //--- get TChains and modify event selections and weights
    gStyle->SetOptStat(0);
    map<string, TChain*> tchains = getTChains(options);
    TCut reweight_region = getReweightRegion(options);

    //--- store reweighting histograms
    TH1F* histoZ;
    map<string, TH1F*> hratios;
    vector<string> processes = {"data", "tt", "vv", "zjets", "photon"};

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

        map<string, TH1F*> hists;
        for (auto process : processes)
            hists[process] = new TH1F(process.c_str(), "", n_reweighting_bins, &reweighting_bins[0]);

        //--- fill reweighting histograms
        if (options.is_data) {
            tchains["data"]->Draw((reweight_var+">>data").c_str(), reweight_region, "goff");
            tchains["tt"]->Draw((reweight_var+">>tt").c_str(), reweight_region*cuts::bkg_weight, "goff");
            tchains["vv"]->Draw((reweight_var+">>vv").c_str(), reweight_region*cuts::bkg_weight, "goff");
            cout << "data integral          : " << hists["data"]->Integral(0, n_reweighting_bins+1) << endl;
            cout << "ttbar integral         : " << hists["tt"]->Integral(0, n_reweighting_bins+1) << endl;
            cout << "diboson integral       : " << hists["vv"]->Integral(0, n_reweighting_bins+1) << endl;
        }
        else {
            tchains["zjets"]->Draw((reweight_var+">>zjets").c_str(), reweight_region*cuts::bkg_weight, "goff");
            cout << "Z+jets integral        : " << hists["zjets"]->Integral(0, n_reweighting_bins+1) << endl;
        }

        tchains["photon"]->Draw((reweight_var+">>photon").c_str(), reweight_region*cuts::photon_weight, "goff");
        cout << "photon integral      " << hists["photon"]->Integral(0, n_reweighting_bins+1) << endl;
        cout << endl;

        //--- get ratio
        if (options.is_data) {
            histoZ = (TH1F*) hists["data"]->Clone("zjets");
            histoZ->Add(hists["tt"], -1.0);
            histoZ->Add(hists["vv"], -1.0);
        }
        else histoZ = (TH1F*) hists["zjets"]->Clone("zjets");

        hratios[reweight_var] = (TH1F*) histoZ->Clone(("hratio_" + reweight_var).c_str());
        hratios[reweight_var]->Divide(hists["photon"]);

        cout << "photon integral        : " << hists["photon"]->Integral(0, n_reweighting_bins+1) << endl;
        cout << "bkg integral           : " << histoZ->Integral(0, n_reweighting_bins+1) << endl;
        cout << "scaling factor         : " << hratios[reweight_var]->Integral(0, n_reweighting_bins+1) << endl;
        cout << endl;

        //--- delete pointers
        for (auto &[key, hist] : hists) delete hist;
    }

    return hratios;
}

map<string, TH2F*> GetSimple2DReweightingHistograms(ReweightingOptions options) {
    //--- get TChains and modify event selections and weights
    gStyle->SetOptStat(0);
    map<string, TChain*> tchains = getTChains(options);
    TCut reweight_region = getReweightRegion(options);

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
            histogram_fills.push_back(make_tuple("data", tchains["data"], hdata));
            histogram_fills.push_back(make_tuple("ttbar", tchains["tt"], htt));
            histogram_fills.push_back(make_tuple("diboson", tchains["vv"], hvv));
        }
        else {
            histogram_fills.push_back(make_tuple("Z+jets", tchains["zjets"], hz));
        }
        histogram_fills.push_back(make_tuple("photon", tchains["photon"], histoG));

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

#include "Settings.cpp"

using namespace std;

//------------------
// HELPER FUNCTIONS
//------------------

map<string, TChain*> getTChains(ReweightingOptions options) {
    //--- open files and create TChains
    map<string, string> filenames;
    filenames["data"] = options.reduction_folder + options.period + "_data_bkg.root";
    filenames["tt"] = options.reduction_folder + options.mc_period + "_ttbar.root";
    filenames["vv"] = options.reduction_folder + options.mc_period + "_diboson.root";
    filenames["zjets"] = options.reduction_folder + options.mc_period + "_Zjets.root";
    filenames["photon"] = options.in_file_name;

    map<string, TChain*> tchains;

    for (auto process : options.processes) {
        cout << padString("Opening " + process + "file") << filenames[process] << endl;
        tchains[process] = new TChain(options.out_tree_name.c_str());
        tchains[process]->Add(filenames[process].c_str());
        cout << padString(process + " entries") << tchains[process]->GetEntries() << endl;
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

TH1F* get1DReweightingRatioHist(ReweightingOptions options, map<string, TChain*> tchains, string reweight_var) {
    //--- set reweighting properties based on variable
    int n_reweighting_bins = bins::n_reweighting_bins.at(reweight_var);
    vector<double> reweighting_bins = bins::reweighting_bins.at(reweight_var);

    map<string, TH1F*> hists;
    for (auto process : options.processes)
        hists[process] = new TH1F(process.c_str(), "", n_reweighting_bins, &reweighting_bins[0]);

    //--- fill reweighting histograms
    if (options.is_data) {
        tchains["data"]->Draw((reweight_var+">>data").c_str(), options.reweight_region, "goff");
        tchains["tt"]->Draw((reweight_var+">>tt").c_str(), options.reweight_region*cuts::bkg_weight, "goff");
        tchains["vv"]->Draw((reweight_var+">>vv").c_str(), options.reweight_region*cuts::bkg_weight, "goff");
        cout << "data integral          : " << hists["data"]->Integral(0, n_reweighting_bins+1) << endl;
        cout << "ttbar integral         : " << hists["tt"]->Integral(0, n_reweighting_bins+1) << endl;
        cout << "diboson integral       : " << hists["vv"]->Integral(0, n_reweighting_bins+1) << endl;
    }
    else {
        tchains["zjets"]->Draw((reweight_var+">>zjets").c_str(), options.reweight_region*cuts::bkg_weight, "goff");
        cout << "Z+jets integral        : " << hists["zjets"]->Integral(0, n_reweighting_bins+1) << endl;
    }

    tchains["photon"]->Draw((reweight_var+">>photon").c_str(), options.reweight_region*cuts::photon_weight, "goff");
    cout << "photon integral      " << hists["photon"]->Integral(0, n_reweighting_bins+1) << endl;
    cout << endl;

    //--- get ratio
    TH1F* hratio;
    if (options.is_data) {
        hratio = (TH1F*) hists["data"]->Clone(("hratio_" + reweight_var).c_str());
        hratio->Add(hists["tt"], -1.0);
        hratio->Add(hists["vv"], -1.0);
    }
    else hratio = (TH1F*) hists["zjets"]->Clone(("hratio_" + reweight_var).c_str());

    cout << "photon integral        : " << hists["photon"]->Integral(0, n_reweighting_bins+1) << endl;
    cout << "bkg integral           : " << hratio->Integral(0, n_reweighting_bins+1) << endl;

    hratio->Divide(hists["photon"]);
    cout << "scaling factor         : " << hratio->Integral(0, n_reweighting_bins+1) << endl;
    cout << endl;

    //--- delete pointers
    for (auto &[key, hist] : hists) delete hist;

    return hratio;
}

//-----------------------
// REWEIGHTING FUNCTIONS
//-----------------------

void ReweightSample(ReweightingOptions options) {
    //---------------------------------------------
    // open files and make TChains
    //---------------------------------------------
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    cout << BOLD(PBLU("Reweighting photon events")) << endl;
    cout << endl;

    cout << padString("Period") << options.period << endl;
    cout << padString("Channel") << options.channel << endl;
    cout << padString("Opening file") << options.in_file_name << endl;
    cout << padString("Reading tree") << options.out_tree_name << endl;

    TFile* smeared_file = new TFile(options.in_file_name.c_str(),"update");          
    TTree* outputTree = (TTree*)smeared_file->Get(options.out_tree_name.c_str());
    cout << padString("Events in ntuple") << outputTree->GetEntries() << endl;
    cout << endl;

    if (options.is_data) options.processes = {"data", "tt", "vv", "photon"};
    else options.processes = {"zjets", "photon"};

    map<string, TChain*> tchains = getTChains(options);
    options.reweight_region = getReweightRegion(options);

    //---------------------------------------------
    // get reweighting histograms
    //---------------------------------------------
    struct ReweightHist {
        int dim;
        TH1F *h1d;
        TH2F *h2d;
    };
    
    map<string, ReweightHist> reweight_hists;

    for (auto reweight_var : options.reweight_vars) {
        //--- split variable by ":"
        std::stringstream unsplit_vars(reweight_var);
        std::string segment;
        std::vector<std::string> split_vars;
        while(std::getline(unsplit_vars, segment, ':')) {
           split_vars.push_back(segment);
        }

        //--- get 1D or 2D reweighting histogram
        if (split_vars.size() == 1) {
            TH1F* hist_ratio = get1DReweightingRatioHist(options, tchains, split_vars[0]);
            reweight_hists[reweight_var] = {1, hist_ratio, 0};
        }
        else {
            //TH2F* hist_ratio = get2DReweightingRatioHist(options, {split_vars[0], split_vars[1]});
            //reweight_hists[reweight_var] = {2, 0, hist_ratio};
        }
    }

    //---------------------------------------------
    // fill reweighting branches
    //---------------------------------------------
    //map<string, float> var_vals;
    //map<string, Float_t> reweight_vals;
    //map<string, TBranch*> reweight_branches;

    //for (auto reweight_var : options.reweight_vars) {
        //var_vals[reweight_var] = 0.0;
        //reweight_vals[reweight_var] = 0.0;
        //SetInputBranch(outputTree, reweight_var, &var_vals[reweight_var]);
        //reweight_branches[reweight_var] = outputTree->Branch(("reweight_"+reweight_var).c_str(),
            //&reweight_vals[reweight_var], ("reweight_"+reweight_var+"/F").c_str());

        //Long64_t nentries = outputTree->GetEntries();
        //for (Long64_t i=0; i<nentries; i++) {
            //if (fmod(i,1e5)==0) cout << i << " events processed.\r" << flush;
            //outputTree->GetEntry(i);

            ////float gamma_var_truncated = var_vals;
            ////if(gamma_var_truncated < reweighting_bins[0]) gamma_var_truncated = reweighting_bins[0];
            ////if(gamma_var_truncated > reweighting_bins[n_reweighting_bins]) gamma_var_truncated = reweighting_bins[n_reweighting_bins];
            //for (auto reweight_var : options.reweight_vars) {
                //int var_bin = reweight_hists[reweight_var]->FindBin(var_vals[reweight_var]);
                //reweight_vals[reweight_var] = reweight_hists[reweight_var]->GetBinContent(var_bin);
                //reweight_branches[reweight_var]->Fill();
            //}
        //}
    //}
    cout << endl;

    outputTree->Write();

    cout << PBLU("Finished reweighting") << endl;
    delete smeared_file;
}

void ReweightPhotons(ReweightingOptions options) {
    if (options.is_data) {
        options.in_file_name = options.reweighting_folder + options.data_period + "_data_photon_" + options.channel + ".root";
        options.out_file_name = options.in_file_name;
        ReweightSample(options);

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

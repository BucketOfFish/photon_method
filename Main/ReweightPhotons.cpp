#include "Settings.cpp"

using namespace std;

//------------
// DATA TYPES
//------------

struct BranchType {
    branch_type type;
    Int_t int_val;
    Float_t float_val;
};

struct ReweightHist {
    int dim;
    TH1F *h1d;
    TH2F *h2d;
};

//------------------
// HELPER FUNCTIONS
//------------------

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

vector<string> splitVars(string reweight_var) {
    /// Split variable by "+".
    stringstream unsplit_vars(reweight_var);
    string segment;
    vector<string> split_vars;
    while(std::getline(unsplit_vars, segment, '+')) {
       split_vars.push_back(segment);
    }

    return split_vars;
}

string combineVars(vector<string> split_vars) {
    /// Combine variables into a single string.
    string reweight_var = "";
    for (auto var : split_vars)
        reweight_var = reweight_var + "+" + var;
    reweight_var.erase(0);

    return reweight_var;
}

//------------------
// PREP FOR READING
//------------------

tuple<TFile*, TTree*> cloneTree(ReweightingOptions options) {
    //--- open files and make TChains
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    cout << BOLD(PBLU("Reweighting photon events")) << endl;
    cout << endl;

    cout << padString("Period") << options.period << endl;
    cout << padString("Channel") << options.channel << endl;
    cout << padString("Opening smeared file") << options.in_file_name << endl;
    cout << padString("Reading tree") << options.out_tree_name << endl;
    cout << padString("Saving reweighted file") << options.out_file_name << endl;

    TFile* input_file = new TFile(options.in_file_name.c_str(), "read");          
    TTree* input_tree = (TTree*)input_file->Get(options.in_tree_name.c_str());

    TFile* output_file = new TFile(options.out_file_name.c_str(), "recreate");          
    TTree* output_tree = input_tree->CloneTree();

    cout << padString("Events in ntuple") << output_tree->GetEntries() << endl;
    cout << endl;

    return make_tuple(output_file, output_tree);
}

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
        cout << padString("Opening " + process + " file") << filenames[process] << endl;
        tchains[process] = new TChain(options.out_tree_name.c_str());
        tchains[process]->Add(filenames[process].c_str());
        cout << padString(process + " entries") << tchains[process]->GetEntries() << endl;
    }
    cout << endl;

    return tchains;
}

//------------------------
// GET FEATURE HISTOGRAMS
//------------------------

map<string, ReweightHist> getReweightingHists(ReweightingOptions options, map<string, TChain*> tchains, string unsplit_vars) {
    /// Fill histograms for all relevant processes.
    vector<string> reweight_vars = splitVars(unsplit_vars);

    //--- set reweighting properties for each reweighting variable
    vector<int> n_reweighting_bins;
    vector<vector<double>> reweighting_bins;
    for (auto reweight_var : reweight_vars) {
        n_reweighting_bins.push_back(bins::n_reweighting_bins.at(reweight_var));
        reweighting_bins.push_back(bins::reweighting_bins.at(reweight_var));
    }

    //--- initialize, fill, and save reweighting histograms
    map<string, ReweightHist> hists;
    for (auto process : options.processes) {
        hists[process].dim = reweight_vars.size();
        TCut weight = cuts::bkg_weight;
        if (process == "data") weight = "1";
        else if (process == "photon") weight = cuts::photon_weight;
        if (hists[process].dim == 1) {
            string name = unsplit_vars + "_" + process;
            hists[process].h1d = new TH1F(name.c_str(), "", n_reweighting_bins[0], &(reweighting_bins[0][0]));
            tchains[process]->Draw((reweight_vars[0] + ">>" + name).c_str(), options.reweight_region * weight, "goff");
            hists[process].h1d->Write(("reweight_" + name).c_str());
        }
        else if (hists[process].dim == 2) {
            string name = unsplit_vars + "_" + process;
            // for some reason, the first variable given goes on the y axis
            hists[process].h2d = new TH2F(name.c_str(), "", n_reweighting_bins[1], &(reweighting_bins[1][0]),
                                    n_reweighting_bins[0], &(reweighting_bins[0][0]));
            tchains[process]->Draw((reweight_vars[1] + ":" + reweight_vars[0] + ">>" + name).c_str(),
                                    options.reweight_region * weight, "goff");
            hists[process].h2d->Write(("reweight_" + name).c_str());
        }
    }

    return hists;
}

//----------------------
// GET RATIO HISTOGRAMS
//----------------------

ReweightHist getReweightingRatioHist(ReweightingOptions options, map<string, ReweightHist> hists, string reweight_var) {
    /// Given filled histograms for various processes, return a histogram for the Z/photon ratio.
    ReweightHist hratio;
    vector<string> split_vars = splitVars(reweight_var);
    hratio.dim = split_vars.size();

    if (options.is_data) {
        if (hratio.dim == 1) {
            hratio.h1d = (TH1F*) hists["data"].h1d->Clone(("hratio_" + reweight_var).c_str());
            hratio.h1d->Add(hists["tt"].h1d, -1.0);
            hratio.h1d->Add(hists["vv"].h1d, -1.0);
            //hists["data"].h1d->Write("data_2D.png");
            //hists["tt"].h1d->Write("tt_2D.png");
            //hists["vv"].h1d->Write("vv_2D.png");
            //hratio.h1d->Write("total_2D.png");
            //hists["photon"].h1d->Write("photon_2D.png");
        }
        else if (hratio.dim == 2) {
            hratio.h2d = (TH2F*) hists["data"].h2d->Clone(("hratio_" + reweight_var).c_str());
            hratio.h2d->Add(hists["tt"].h2d, -1.0);
            hratio.h2d->Add(hists["vv"].h2d, -1.0);
        }
    }
    else {
        if (hratio.dim == 1)
            hratio.h1d = (TH1F*) hists["zjets"].h1d->Clone(("hratio_" + reweight_var).c_str());
        else if (hratio.dim == 2)
            hratio.h2d = (TH2F*) hists["zjets"].h2d->Clone(("hratio_" + reweight_var).c_str());
    }

    if (hratio.dim == 1) {
        cout << "\t" << padString("photon integral") << hists["photon"].h1d->Integral(0, hists["photon"].h1d->GetNbinsX()+1) << endl;
        cout << "\t" << padString("bkg integral") << hratio.h1d->Integral(0, hratio.h1d->GetNbinsX()+1) << endl;
        hratio.h1d->Divide(hists["photon"].h1d);
        hratio.h1d->Write(("ratio_" + reweight_var).c_str());
    }
    else if (hratio.dim == 2) {
        cout << "\t" << padString("photon integral") << hists["photon"].h2d->Integral(0, hists["photon"].h2d->GetNbinsX()+1,
                                                            0, hists["photon"].h2d->GetNbinsY()+1) << endl;
        cout << "\t" << padString("bkg integral") << hratio.h2d->Integral(0, hratio.h2d->GetNbinsX()+1,
                                                            0, hratio.h2d->GetNbinsY()+1) << endl;
        hratio.h2d->Divide(hists["photon"].h2d);
        hratio.h2d->Write(("ratio_" + reweight_var).c_str());
    }
    cout << endl;


    return hratio;
}

//-------------------------
// FILL REWEIGHTING BRANCH
//-------------------------

void fillReweightingBranches(ReweightingOptions options, TTree* output_tree, map<string, ReweightHist> reweight_hists) {
    map<string, BranchType> rw_feature_vals;
    map<string, Float_t> rw_weight_vals;
    map<string, TBranch*> rw_branches;

    //--- set variables to read features from input branches, and write weights to new branches
    for (auto unsplit_var : options.reweight_vars) {
        vector<string> split_vars = splitVars(unsplit_var);
        for (auto reweight_var : split_vars) {
            BranchType newBranch;
            rw_feature_vals[reweight_var] = newBranch;
            if (bins::reweighting_type.at(reweight_var) == INT) {
                rw_feature_vals[reweight_var].type = INT;
                rw_feature_vals[reweight_var].int_val = 0;
                SetInputBranch(output_tree, reweight_var, &rw_feature_vals[reweight_var].int_val);
            }
            else if (bins::reweighting_type.at(reweight_var) == FLOAT) {
                rw_feature_vals[reweight_var].type = FLOAT;
                rw_feature_vals[reweight_var].float_val = 0.0;
                SetInputBranch(output_tree, reweight_var, &rw_feature_vals[reweight_var].float_val);
            }
        }
        rw_weight_vals[unsplit_var] = 0.0;
        string branch_name = "reweight_" + unsplit_var;
        rw_branches[unsplit_var] = output_tree->Branch(branch_name.c_str(), &rw_weight_vals[unsplit_var], (branch_name + "/F").c_str());
    }

    //--- fill new branches
    Long64_t nentries = output_tree->GetEntries();
    for (Long64_t i=0; i<nentries; i++) {
        if (fmod(i,1e5)==0) cout << i << " events processed.\r" << flush;
        output_tree->GetEntry(i);

        for (auto unsplit_var : options.reweight_vars) {
            vector<string> split_vars = splitVars(unsplit_var);
            vector<int> feature_bins;
            for (auto reweight_var : split_vars) {
                int feature_bin = 0;
                while (true) {
                    if (bins::reweighting_type.at(reweight_var) == INT) {
                        if (bins::reweighting_bins.at(reweight_var)[feature_bin] > rw_feature_vals[reweight_var].int_val) break;
                    }
                    else if (bins::reweighting_type.at(reweight_var) == FLOAT) {
                        if (bins::reweighting_bins.at(reweight_var)[feature_bin] > rw_feature_vals[reweight_var].float_val) break;
                    }
                    feature_bin++;
                    if (feature_bin >= bins::n_reweighting_bins.at(reweight_var)) break;
                }
                feature_bins.push_back(feature_bin);
            }
            if (split_vars.size() == 1)
                rw_weight_vals[unsplit_var] = reweight_hists[unsplit_var].h1d->GetBinContent(feature_bins[0]);
            else if (split_vars.size() == 2)
                rw_weight_vals[unsplit_var] = reweight_hists[unsplit_var].h2d->GetBinContent(feature_bins[0], feature_bins[1]);
            //float gamma_var_truncated = rw_branches[unsplit_var];
            //if(gamma_var_truncated < reweighting_bins[0]) gamma_var_truncated = reweighting_bins[0];
            //if(gamma_var_truncated > reweighting_bins[n_reweighting_bins]) gamma_var_truncated = reweighting_bins[n_reweighting_bins];
            rw_branches[unsplit_var]->Fill();

            cout << unsplit_var << ", " << rw_feature_vals[unsplit_var].int_val << ", " << rw_feature_vals[unsplit_var].float_val <<
                ", " << rw_weight_vals[unsplit_var] << endl;
        }
    }
    cout << endl << endl;

    output_tree->Write();
}

//------------
// UNIT TESTS
//------------

void RunUnitTests(ReweightingOptions options) {
    options.period = "data15-16";
    options.data_period = DataPeriod(options.period);
    options.mc_period = getMCPeriod(options.period);
    options.channel = "ee";
    options.in_file_name = options.unit_test_folder + "SmearedNtuples/" + options.data_period + "_data_photon_" + options.channel + ".root";
    options.out_file_name = "test.root";
    options.reduction_folder = options.unit_test_folder + "ReducedNtuples/";
    options.reweight_vars = {"Ptll", "nJet30", "Ptll+Ht30"};

    //--- tree cloning test
    auto [output_file, output_tree] = cloneTree(options);
    map<string, vector<float>> check_branch_vals;
    check_branch_vals["Ptll"] = {243.606, 138.736, 40.529, 70.656, 257.859, 82.166, 107.491, 476.676, 55.058, 134.012};
    check_branch_vals["Ht30"] = {224.067, 136.098, 50.745, 84.949, 203.888, 84.017, 123.265, 914.663, 39.261, 194.367};
    vector<string> branches_to_check = {"Ptll", "Ht30"};
    map<string, float> branch_vals;
    for (auto branch : branches_to_check) {
        branch_vals[branch] = 0;
        SetInputBranch(output_tree, branch, &branch_vals[branch]);
    }
    bool all_match = true;
    for (int i = 0; i < 10; i++) {
        output_tree->GetEntry(i);
        for (auto branch : branches_to_check) {
            if (fabs(branch_vals[branch] - check_branch_vals[branch][i]) > 0.001) all_match = false;
        }
    }
    if (all_match) passTest("Passed cloning test");
    else failTest("Failed cloning test");
    cout << endl;

    //--- TChain creation test
    map<string, TChain*> tchains = getTChains(options);
    if (tchains["data"]->GetEntries() == 1000000 && tchains["tt"]->GetEntries() == 1000000) passTest("Passed TChain test");
    else failTest("Failed TChain test");
    cout << endl;

    //--- reweighting region with channel test
    options.reweight_region = getReweightRegion(options);
    if (options.reweight_region == "(nJet30>=2 && lepPt[0]>25.0 && lepPt[1]>25.0 && nLep_signal==2)&&(channel==1)")
        passTest("Passed reweighting region test");
    else failTest("Failed reweighting region test");
    cout << endl;

    //--- reweighting histograms test
    map<string, map<string, ReweightHist>> hists;
    for (auto vars : options.reweight_vars) {
        hists[vars] = getReweightingHists(options, tchains, vars);
    }
    if (abs(hists["Ptll"]["tt"].h1d->GetBinContent(1) - 3089.1179) < 0.001 &&
        abs(hists["Ptll"]["tt"].h1d->GetBinContent(10) - 4412.0952) < 0.001 &&
        abs(hists["Ptll"]["vv"].h1d->GetBinContent(8) - 248.56230) < 0.001 &&
        abs(hists["Ptll"]["vv"].h1d->GetBinContent(23) - 3.4948041) < 0.001 &&
        abs(hists["nJet30"]["data"].h1d->GetBinContent(3) - 84967.000) < 1 &&
        abs(hists["Ptll+Ht30"]["tt"].h2d->GetBinContent(1, 8) - 32.0249) < 0.001 &&
        abs(hists["Ptll+Ht30"]["photon"].h2d->GetBinContent(5, 9) - 205209.000) < 1)
        passTest("Passed feature histogram test");
    else failTest("Failed feature histogram test");
    cout << endl;

    //--- reweighting histogram ratios test
    map<string, ReweightHist> reweight_hists;
    for (auto vars : options.reweight_vars) {
        reweight_hists[vars] = getReweightingRatioHist(options, hists[vars], vars);
    }
    if (abs(reweight_hists["Ptll"].h1d->GetBinContent(2) - 0.000830079) < 0.00001 &&
        abs(reweight_hists["Ptll"].h1d->GetBinContent(5) - 0.0023039) < 0.00001 &&
        abs(reweight_hists["nJet30"].h1d->GetBinContent(5) - 0.00282705) < 0.00001 &&
        abs(reweight_hists["Ptll+Ht30"].h2d->GetBinContent(2, 8) - 0.00088203) < 0.00001 &&
        abs(reweight_hists["Ptll+Ht30"].h2d->GetBinContent(5, 9) - 0.00260129) < 0.00001)
        passTest("Passed histogram ratio test");
    else failTest("Failed histogram ratio test");
    cout << endl;

    //--- branch filling test
    fillReweightingBranches(options, output_tree, reweight_hists);

    map<string, vector<float>> rw_weight_checks;
    rw_weight_checks["Ptll"] = {0.0373285, 0.0152829, 0.0012275, 0.0052744, 0.0373285, 0.0070094, 0.0106032, 0.0405640,
                                    0.0031431, 0.0152829};
    rw_weight_checks["nJet30"] = {0, 0, 0, 0, 0, 0, 0, 0.00315326, 0, 0};
    rw_weight_checks["Ptll+Ht30"] = {0.056162, 0.0196737, 0, 0.0077757, 0.0518713, 0.0114504, 0.0123719, 0.0355054, 0, 0.0139699};

    map<string, Float_t> rw_weights;
    for (auto unsplit_vars : options.reweight_vars) {
        rw_weights[unsplit_vars] = 0.0;
        SetInputBranch(output_tree, "reweight_"+unsplit_vars, &rw_weights[unsplit_vars]);
    }
    
    all_match = true;
    for (Long64_t i=0; i<10; i++) {
        output_tree->GetEntry(i);

        for (auto unsplit_vars : options.reweight_vars) {
            if (abs(rw_weights[unsplit_vars] - rw_weight_checks[unsplit_vars][i]) > 0.01) all_match = false;
            //cout << unsplit_vars << ": " << rw_weights[unsplit_vars] << " " << rw_weight_checks[unsplit_vars][i] << endl;
        }
    }
    if (all_match) passTest("Passed branch filling test");
    else failTest("Failed branch filling test");
    cout << endl;

    ////--- save diagnostic histograms
    //for (auto unsplit_vars : options.reweight_vars) {
        //TCanvas *can = new TCanvas("can","can",600,600);
        //can->cd();
        //TPad* namepad = new TPad("namepad","namepad",0.0,0.0,1.0,1.0);
        //namepad->Draw();
        //namepad->cd();

        //THStack *hs_bkg = new THStack("hs_bkg","");
        //if (hists[unsplit_vars]["photon"].dim == 1) {
            //hs_bkg->Add(hists[unsplit_vars]["data"].h1d);
            ////hists[unsplit_vars]["tt"].h1d->Scale(-1);
            ////hs_bkg->Add(hists[unsplit_vars]["tt"].h1d);
            ////hists[unsplit_vars]["vv"].h1d->Scale(-1);
            ////hs_bkg->Add(hists[unsplit_vars]["vv"].h1d);
        //}
        //else if (hists[unsplit_vars]["photon"].dim == 2) {
            //hs_bkg->Add(hists[unsplit_vars]["data"].h2d);
            ////hists[unsplit_vars]["tt"].h2d->Scale(-1);
            ////hs_bkg->Add(hists[unsplit_vars]["tt"].h2d);
            ////hists[unsplit_vars]["vv"].h2d->Scale(-1);
            ////hs_bkg->Add(hists[unsplit_vars]["vv"].h2d);
        //}
        //hs_bkg->Draw("hist");

        //string name = unsplit_vars + "_photon";
        //output_tree->Draw((unsplit_vars + ">>" + name).c_str(), options.reweight_region * cuts::photon_weight *
            //("reweight_"+unsplit_vars).c_str(), "goff");
        //if (hists[unsplit_vars]["photon"].dim == 1)
            //hists[unsplit_vars]["photon"].h1d->Draw("samehist");
        //else if (hists[unsplit_vars]["photon"].dim == 2)
            //hists[unsplit_vars]["photon"].h2d->Draw("samehist");

        //can->Write((unsplit_vars + "_compare.eps").c_str());
        //delete can, namepad, hs_bkg;
    //}
    //passTest("Saved diagnostic histograms");
    //cout << endl;

    passTest("All reweighting tests passed");
    cout << endl;
    output_file->Close();
    remove(options.out_file_name.c_str());
}

//----------------------
// REWEIGHTING FUNCTION
//----------------------

void ReweightSample(ReweightingOptions options) {
    auto [output_file, output_tree] = cloneTree(options);
    map<string, TChain*> tchains = getTChains(options);

    options.reweight_region = getReweightRegion(options);

    map<string, ReweightHist> reweight_hists;
    for (auto vars : options.reweight_vars) {
        map<string, ReweightHist> hists = getReweightingHists(options, tchains, vars);
        reweight_hists[vars] = getReweightingRatioHist(options, hists, vars);
    }

    fillReweightingBranches(options, output_tree, reweight_hists);
    output_file->Close();

    cout << PBLU("Finished reweighting") << endl;
}

void ReweightPhotons(ReweightingOptions options) {
    if (options.unit_testing)
        RunUnitTests(options);
    else {
        if (options.is_data) {
            options.in_file_name = options.smearing_folder + options.data_period + "_data_photon_" + options.channel + ".root";
            options.out_file_name = options.reweighting_folder + options.data_period + "_data_photon_" + options.channel + ".root";
            ReweightSample(options);

            options.in_file_name = options.smearing_folder + options.mc_period + "_Vgamma_" + options.channel + ".root";
            options.out_file_name = options.reweighting_folder + options.mc_period + "_Vgamma_" + options.channel + ".root";
            ReweightSample(options);
        }
        else {
            options.in_file_name = options.reweighting_folder + options.mc_period + "_SinglePhoton222_" + options.channel + ".root";
            options.out_file_name = options.in_file_name;
            ReweightSample(options);
        }
    }
}

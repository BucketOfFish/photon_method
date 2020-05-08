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

TCut getReweightRegion(Options options) {
    TCut reweight_region = cuts::selections["reweight"];
    reweight_region += cuts::selections[options.channel];

    cout << "bkg selection          : " << reweight_region.GetTitle() << endl;
    cout << "bkg weight             : " << cuts::bkg_weight.GetTitle() << endl;
    //cout << "photon selection       : " << reweight_region.GetTitle() << endl;
    //cout << "photon weight          : " << cuts::photon_weight.GetTitle() << endl;
    cout << endl;

    return reweight_region;
}

bool isInZWindow(string reweight_var) {
    return (reweight_var.size() > 7 && reweight_var.rfind("Zwindow") == reweight_var.size()-7);
}

vector<string> splitVars(string reweight_var) {
    vector<string> split_vars;
    string delimiter = "__";
    string segment;
    size_t pos = 0;
    while ((pos = reweight_var.find(delimiter)) != string::npos) {
        segment = reweight_var.substr(0, pos);
        split_vars.push_back(segment);
        reweight_var.erase(0, pos + delimiter.length());
    }
    if (reweight_var != "Zwindow") split_vars.push_back(reweight_var);

    return split_vars;
}

//------------------
// PREP FOR READING
//------------------

tuple<TFile*, TTree*> cloneTree(Options options) {
    //--- open files and make TChains
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);

    cout << BOLD(PBLU("Reweighting photon events")) << endl;
    cout << endl;

    cout << padString("Period") << options.period << endl;
    cout << padString("Channel") << options.channel << endl;
    cout << padString("Opening smeared file") << options.smearing_file_name << endl;
    cout << padString("Reading tree") << options.tree_name << endl;
    cout << padString("Saving reweighted file") << options.reweighting_file_name << endl;

    TFile* input_file = new TFile(options.smearing_file_name.c_str(), "read");          
    TTree* input_tree = (TTree*)input_file->Get(options.tree_name.c_str());

    TFile* output_file = new TFile(options.reweighting_file_name.c_str(), "recreate");          
    TTree* output_tree = input_tree->CloneTree();

    cout << padString("Events in ntuple") << output_tree->GetEntries() << endl;
    cout << endl;

    return make_tuple(output_file, output_tree);
}

map<string, TChain*> getTChains(Options options) {
    //--- open files and create TChains
    map<string, string> filenames;
    filenames["data"] = options.reduction_folder + options.period + "_data_bkg.root";
    filenames["tt"] = options.reduction_folder + options.mc_period + "_ttbar.root";
    filenames["vv"] = options.reduction_folder + options.mc_period + "_diboson.root";
    filenames["zjets"] = options.reduction_folder + options.mc_period + "_Zjets.root";
    filenames["photon"] = options.smearing_file_name;

    map<string, TChain*> tchains;

    for (auto process : options.processes) {
        cout << padString("Opening " + process + " file") << filenames[process] << endl;
        tchains[process] = new TChain(options.tree_name.c_str());
        tchains[process]->Add(filenames[process].c_str());
        cout << padString(process + " entries") << tchains[process]->GetEntries() << endl;
    }
    cout << endl;

    return tchains;
}

//------------------------
// GET FEATURE HISTOGRAMS
//------------------------

map<string, ReweightHist> getReweightingHists(Options options, map<string, TChain*> tchains, string unsplit_vars) {
    /// Fill histograms for all relevant processes.
    vector<string> reweight_vars = splitVars(unsplit_vars);
    bool in_z_window = isInZWindow(unsplit_vars);
    TCut reweight_region = options.reweight_region;
    if (in_z_window) reweight_region += "(mll>81 && mll<101)";

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
            tchains[process]->Draw((reweight_vars[0] + ">>" + name).c_str(), reweight_region * weight, "goff");
            hists[process].h1d->Write(("reweight_" + name).c_str());
        }
        else if (hists[process].dim == 2) {
            string name = unsplit_vars + "_" + process;
            // for some reason, the first variable given goes on the y axis
            hists[process].h2d = new TH2F(name.c_str(), "", n_reweighting_bins[1], &(reweighting_bins[1][0]),
                                    n_reweighting_bins[0], &(reweighting_bins[0][0]));
            tchains[process]->Draw((reweight_vars[1] + ":" + reweight_vars[0] + ">>" + name).c_str(), reweight_region * weight, "goff");
            hists[process].h2d->Write(("reweight_" + name).c_str());
        }
    }

    return hists;
}

//----------------------
// GET RATIO HISTOGRAMS
//----------------------

ReweightHist getReweightingRatioHist(Options options, map<string, ReweightHist> hists, string reweight_var) {
    /// Given filled histograms for various processes, return a histogram for the Z/photon ratio.
    ReweightHist hratio;
    vector<string> split_vars = splitVars(reweight_var);
    hratio.dim = split_vars.size();

    if (options.is_data) {
        if (hratio.dim == 1) {
            hratio.h1d = (TH1F*) hists["data"].h1d->Clone(("hratio_" + reweight_var).c_str());
            for (auto process : options.processes) {
                //hists[process].h1d->Write((process+"_1D.png").c_str());
                if (process != "data" && process != "photon")
                    hratio.h1d->Add(hists[process].h1d, -1.0);
            }
            //hratio.h1d->Write("total_1D.png");
        }
        else if (hratio.dim == 2) {
            hratio.h2d = (TH2F*) hists["data"].h2d->Clone(("hratio_" + reweight_var).c_str());
            for (auto process : options.processes) {
                //hists[process].h2d->Write((process+"_2D.png").c_str());
                if (process != "data" && process != "photon")
                    hratio.h2d->Add(hists[process].h2d, -1.0);
            }
            //hratio.h2d->Write("total_2D.png");
        }
    }
    else {
        if (hratio.dim == 1)
            hratio.h1d = (TH1F*) hists["zjets"].h1d->Clone(("hratio_" + reweight_var).c_str());
        else if (hratio.dim == 2)
            hratio.h2d = (TH2F*) hists["zjets"].h2d->Clone(("hratio_" + reweight_var).c_str());
    }

    if (hratio.dim == 1) {
        //float photon_integral = hists["photon"].h1d->Integral(0, hists["photon"].h1d->GetNbinsX()+1);
        //float bkg_integral = hratio.h1d->Integral(0, hratio.h1d->GetNbinsX()+1);
        //cout << "\t" << padString("photon integral") << photon_integral << endl;
        //cout << "\t" << padString("bkg integral") << bkg_integral << endl;
        //hratio.h1d->Scale(photon_integral/bkg_integral);
        hratio.h1d->Divide(hists["photon"].h1d);
        hratio.h1d->Write(("ratio_" + reweight_var).c_str());
    }
    else if (hratio.dim == 2) {
        //float photon_integral = hists["photon"].h2d->Integral(0, hists["photon"].h2d->GetNbinsX()+1, 0, hists["photon"].h2d->GetNbinsY()+1);
        //float bkg_integral = hratio.h2d->Integral(0, hratio.h2d->GetNbinsX()+1, 0, hratio.h2d->GetNbinsY()+1);
        //cout << "\t" << padString("photon integral") << photon_integral << endl;
        //cout << "\t" << padString("bkg integral") << bkg_integral << endl;
        //hratio.h2d->Scale(photon_integral/bkg_integral);
        hratio.h2d->Divide(hists["photon"].h2d);
        hratio.h2d->Write(("ratio_" + reweight_var).c_str());
    }
    cout << endl;


    return hratio;
}

//-------------------------
// FILL REWEIGHTING BRANCH
//-------------------------

void fillReweightingBranches(Options options, TTree* output_tree, map<string, ReweightHist> reweight_hists) {
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
        rw_branches[unsplit_var] = output_tree->Branch(branch_name.c_str(), &rw_weight_vals[unsplit_var],
            (branch_name + "/F").c_str());
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
                        if (bins::reweighting_bins.at(reweight_var)[feature_bin] >
                            rw_feature_vals[reweight_var].int_val) break;
                    }
                    else if (bins::reweighting_type.at(reweight_var) == FLOAT) {
                        if (bins::reweighting_bins.at(reweight_var)[feature_bin] >
                            rw_feature_vals[reweight_var].float_val) break;
                    }
                    feature_bin++;
                    if (feature_bin > bins::n_reweighting_bins.at(reweight_var)) break;
                }
                feature_bins.push_back(feature_bin);
            }
            if (split_vars.size() == 1)
                rw_weight_vals[unsplit_var] = reweight_hists[unsplit_var].h1d->GetBinContent(feature_bins[0]);
            else if (split_vars.size() == 2)
                rw_weight_vals[unsplit_var] = reweight_hists[unsplit_var].h2d->GetBinContent(feature_bins[0],
                    feature_bins[1]);
            //float gamma_var_truncated = rw_branches[unsplit_var];
            //if(gamma_var_truncated < reweighting_bins[0]) gamma_var_truncated = reweighting_bins[0];
            //if(gamma_var_truncated > reweighting_bins[n_reweighting_bins]) gamma_var_truncated = reweighting_bins[n_reweighting_bins];
            rw_branches[unsplit_var]->Fill();
        }
    }
    cout << endl << endl;

    output_tree->Write();
}

//------------
// UNIT TESTS
//------------

void performReweightingUnitTests(Options options) {
    options.period = "data15-16";
    options.data_period = DataPeriod(options.period);
    options.mc_period = getMCPeriod(options.period);
    options.is_data = true;
    options.processes = {"data", "tt", "vv", "photon"};
    options.channel = "ee";
    options.smearing_file_name = options.unit_test_folder + "SmearedNtuples/" + options.data_period + "_data_photon_" + options.channel + ".root";
    options.reweighting_file_name = "test.root";
    options.reduction_folder = options.unit_test_folder + "ReducedNtuples/";
    options.reweight_vars = {"Ptll", "nJet30", "Ptll__Ht30", "Ptll__Zwindow"};

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
    if (options.reweight_region == "((nJet30>=2 && lepPt[0]>25.0 && lepPt[1]>25.0 && nLep_signal==2)&&(abs(lepFlavor[0])==abs(lepFlavor[1])))&&(channel==1)")
        passTest("Passed reweighting region test");
    else failTest("Failed reweighting region test");
    cout << endl;

    //--- Zwindow detection test
    if (!isInZWindow("nJet30") && isInZWindow("Ptll__Zwindow") && !isInZWindow("Ptll__Ht30") && isInZWindow("Ptll__Ht30__Zwindow"))
        passTest("Passed Z window detection test");
    else failTest("Failed Z window detection test");
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
        abs(hists["Ptll__Ht30"]["tt"].h2d->GetBinContent(1, 8) - 32.0249) < 0.001 &&
        abs(hists["Ptll__Ht30"]["photon"].h2d->GetBinContent(5, 9) - 205209.000) < 1 &&
        abs(hists["Ptll__Zwindow"]["vv"].h1d->GetBinContent(23) - 3.1094) < 0.001)
        passTest("Passed feature histogram test");
    else failTest("Failed feature histogram test");
    cout << endl;

    //--- reweighting histogram ratios test
    map<string, ReweightHist> reweight_hists;
    for (auto vars : options.reweight_vars) {
        reweight_hists[vars] = getReweightingRatioHist(options, hists[vars], vars);
    }
    if (abs(reweight_hists["Ptll"].h1d->GetBinContent(2) - 0.262335) < 0.00001 &&
        abs(reweight_hists["Ptll"].h1d->GetBinContent(5) - 0.728115) < 0.00001 &&
        abs(reweight_hists["nJet30"].h1d->GetBinContent(5) - 0.892541) < 0.00001 &&
        abs(reweight_hists["Ptll__Ht30"].h2d->GetBinContent(2, 8) - 0.278756) < 0.00001 &&
        abs(reweight_hists["Ptll__Ht30"].h2d->GetBinContent(5, 9) - 0.822111) < 0.00001 &&
        abs(reweight_hists["Ptll__Zwindow"].h1d->GetBinContent(5) - 0.768134) < 0.00001)
        passTest("Passed histogram ratio test");
    else failTest("Failed histogram ratio test");
    cout << endl;

    //--- branch filling test
    fillReweightingBranches(options, output_tree, reweight_hists);

    map<string, vector<float>> rw_weight_checks;
    rw_weight_checks["Ptll"] = {12.6804, 6.38236, 0.573366, 2.21522, 12.6804, 3.35102, 4.82996, 15.2658, 1.24056, 6.38236};
    rw_weight_checks["nJet30"] = {0.99553, 0, 0, 0, 0, 0, 0, 1.09558, 0, 0};
    rw_weight_checks["Ptll__Ht30"] = {16.9965, 8.06162, 0, 2.99799, 19.043, 4.94791, 5.8982, 17.0231, 0, 5.4536};
    rw_weight_checks["Ptll__Zwindow"] = {11.6063, 6.82505, 0.59245, 2.48084, 11.6063, 3.75864, 5.35801, 15.0041, 1.2998, 6.82505};

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
    remove(options.reweighting_file_name.c_str());
}

//----------------------
// REWEIGHTING FUNCTION
//----------------------

void ReweightSample(Options options) {
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

void ReweightPhotons(Options options) {
    if (options.unit_testing)
        performReweightingUnitTests(options);
    else {
        if (options.is_data) {
            options.smearing_file_name = options.smearing_folder + options.data_period + "_data_photon_" + options.channel + ".root";
            options.reweighting_file_name = options.reweighting_folder + options.data_period + "_data_photon_" + options.channel + ".root";
            ReweightSample(options);

            options.smearing_file_name = options.smearing_folder + options.mc_period + "_Vgamma_" + options.channel + ".root";
            options.reweighting_file_name = options.reweighting_folder + options.mc_period + "_Vgamma_" + options.channel + ".root";
            ReweightSample(options);
        }
        else {
            options.smearing_file_name = options.reweighting_folder + options.mc_period + "_SinglePhoton222_" + options.channel + ".root";
            options.reweighting_file_name = options.smearing_file_name;
            ReweightSample(options);
        }
    }
}

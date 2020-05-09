#include "Settings.cpp"

using namespace std;

//-----------------
// DATA STRUCTURES
//-----------------

using dataFrameMap = map<string, ROOT::RDataFrame*>;
using weightedDataFrameMap = map<string,
                             unique_ptr<ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void>>>;
using filteredDataFrameMap = map<string, // region
                             map<string, // feature
                             map<string, // SR/CR
                             unique_ptr<ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>>>>>;
using histMap = map<string, // region
                map<string, // feature
                map<string, // process
                ROOT::RDF::RResultPtr<TH1D>>>>;
struct Result {
    map<string, map<string, TH1D*>> histograms; // [feature][process] histogram
    map<string, float> process_yields;
    float photon_SF;
};
struct resultsMap {
    vector<string> plot_regions;
    vector<string> plot_features;
    vector<string> plot_channels;
    string period;
    map<string, Result> results; // results by [region]
};

//------------------
// HELPER FUNCTIONS
//------------------

 vector<string> splitStringBy(string input, string splitter) {
    vector<string> output;
    boost::algorithm::split(output, input, boost::is_any_of(splitter));
    return output;
}

tuple<string, string, string> getPlotRegionInfo(Options options, string channel, string region) {
    if (cuts::selections.count(region) == 0) {
        cout << "Unrecognized region! Exiting." << endl;
        exit(0);
    }
    TCut plot_region = cuts::selections[region];

    plot_region += cuts::selections[channel];
    plot_region += options.additional_plot_cut;

    TCut plot_CR = plot_region + cuts::CR;
    if (cuts::plot_region_met_portions.count(region) > 0) plot_region += cuts::plot_region_met_portions[region];

    string region_name = region + " " + channel;

    return make_tuple(region_name, string(plot_region), string(plot_CR));
}

ROOT::RDF::TH1DModel getHistogramInfo(string plot_feature) {
    map<string, ROOT::RDF::TH1DModel> plot_settings;
    plot_settings["met_Et"] = ROOT::RDF::TH1DModel("", "E_{T}^{miss} [GeV]", 20, 0, 300);
    plot_settings["METl"] = ROOT::RDF::TH1DModel("", "E_{T,||}^{miss} [GeV]", 30, -150, 150);
    plot_settings["METt"] = ROOT::RDF::TH1DModel("", "E_{T,#perp}^{miss} [GeV]", 30, -150, 150);
    plot_settings["met_Sign"] = ROOT::RDF::TH1DModel("", "E_{T}^{miss} significance", 20, 0, 50);
    plot_settings["MET_loose"] = ROOT::RDF::TH1DModel("", "E_{T,loose}^{miss} [GeV]", 20, 0, 200);
    plot_settings["MET_tight"] = ROOT::RDF::TH1DModel("", "E_{T,tight}^{miss} [GeV]", 20, 0, 200);
    plot_settings["MET_tighter"] = ROOT::RDF::TH1DModel("", "E_{T,tighter}^{miss} [GeV]", 20, 0, 200);
    plot_settings["MET_tenacious"] = ROOT::RDF::TH1DModel("", "E_{T,tenacious}^{miss} [GeV]", 20, 0, 200);
    plot_settings["mt2leplsp_0"] = ROOT::RDF::TH1DModel("", "m_{T2}^{0} [GeV]", 20, 0, 500);
    //plot_settings["Ptll"] = ROOT::RDF::TH1DModel("", "p_{T} [GeV]", 25, 0, 1000);
    //plot_settings["Ptll_reweight"] = ROOT::RDF::TH1DModel("", "p_{T} [GeV]", bins::n_reweighting_bins.at("Ptll"),
        //&bins::reweighting_bins.at("Ptll")[0]);
    plot_settings["Ptll"] = ROOT::RDF::TH1DModel("", "p_{T} [GeV]", bins::n_reweighting_bins.at("Ptll"),
        &bins::reweighting_bins.at("Ptll")[0]);
    plot_settings["Z_pt"] = ROOT::RDF::TH1DModel("", "p_{T} [GeV]", 20, 0, 100);
    plot_settings["nJet30"] = ROOT::RDF::TH1DModel("", "n_{jets}", 6, 2, 8);
    plot_settings["jet_n"] = ROOT::RDF::TH1DModel("", "n_{jets}", 6, 2, 8);
    plot_settings["jet_eta"] = ROOT::RDF::TH1DModel("", "jet_{#eta}", 30, -3, 3);
    plot_settings["jet_phi"] = ROOT::RDF::TH1DModel("", "jet_{#phi}", 20, 0, 3.14);
    plot_settings["jetPt"] = ROOT::RDF::TH1DModel("", "jet_{p_{T}} [GeV]", 20, 0, 300);
    plot_settings["jetPt[0]"] = ROOT::RDF::TH1DModel("", "jet_{p_{T},1} [GeV]", 20, 0, 300);
    plot_settings["jetPt[1]"] = ROOT::RDF::TH1DModel("", "jet_{p_{T},2} [GeV]", 20, 0, 300);
    plot_settings["bjet_n"] = ROOT::RDF::TH1DModel("", "n_{b-jets}", 4, 0, 4);
    //plot_settings["Ht30"] = ROOT::RDF::TH1DModel("", "H_{T}", 15, 0, 1500);
    //plot_settings["Ht30_reweight"] = ROOT::RDF::TH1DModel("", "H_{T} [GeV]", bins::n_reweighting_bins.at("Ht30"),
        //&bins::reweighting_bins.at("Ht30")[0]);
    plot_settings["Ht30"] = ROOT::RDF::TH1DModel("", "H_{T} [GeV]", bins::n_reweighting_bins.at("Ht30"),
        &bins::reweighting_bins.at("Ht30")[0]);
    plot_settings["mll"] = ROOT::RDF::TH1DModel("", "m_{ll} [GeV]", 30, 0, 300);
    plot_settings["MT2"] = ROOT::RDF::TH1DModel("", "m_{T2} [GeV]", 20, 0, 200);
    plot_settings["MT2W"] = ROOT::RDF::TH1DModel("", "m_{T2}^{W} [GeV]", 20, 0, 200);
    plot_settings["lepEta"] = ROOT::RDF::TH1DModel("", "lep_{#eta}", 30, -3, 3);
    plot_settings["lepPhi"] = ROOT::RDF::TH1DModel("", "lep_{#phi}", 20, 0, 3.14);
    plot_settings["lepPt"] = ROOT::RDF::TH1DModel("", "lep_{p_{T}} [GeV]", 20, 0, 300);
    plot_settings["lepPt[0]"] = ROOT::RDF::TH1DModel("", "lep_{p_{T},1} [GeV]", 20, 0, 300);
    plot_settings["lepPt[1]"] = ROOT::RDF::TH1DModel("", "lep_{p_{T},2} [GeV]", 20, 0, 200);
    plot_settings["lepEta[0]"] = ROOT::RDF::TH1DModel("", "lep_{#eta,1}", 30, -3, 3);
    plot_settings["lepEta[1]"] = ROOT::RDF::TH1DModel("", "lep_{#eta,2}", 30, -3, 3);
    plot_settings["DPhi_METLepLeading"] = ROOT::RDF::TH1DModel("", "#Delta#phi(lep_{1},E_{T}^{miss})", 20, 0, 3.14);
    plot_settings["DPhi_METLepSecond"] = ROOT::RDF::TH1DModel("", "#Delta#phi(lep_{2},E_{T}^{miss})", 20, 0, 3.14);
    plot_settings["dPhiMetJet1"] = ROOT::RDF::TH1DModel("", "#Delta#phi(jet_{1},E_{T}^{miss})", 20, 0, 3.14);
    plot_settings["dPhiMetJet2"] = ROOT::RDF::TH1DModel("", "#Delta#phi(jet_{2},E_{T}^{miss})", 20, 0, 3.14);
    plot_settings["dPhiMetJet12Min"] = ROOT::RDF::TH1DModel("", "#Delta#phi(jet_{min(1,2)},E_{T}^{miss})", 20, 0, 3.14);
    plot_settings["dPhiPllMet"] = ROOT::RDF::TH1DModel("", "#Delta#phi(p_{T},E_{T}^{miss})", 20, 0, 3.14);

    ROOT::RDF::TH1DModel hist_model = plot_settings[plot_feature];
    return hist_model;
}

//--------------------
// SET UP RDATAFRAMES
//--------------------

dataFrameMap getRDataFrames(Options options) {
    //--- load files
    string ntuple_path = options.reduction_folder;
    string photon_path = options.reweighting_folder;

    vector<string> periods = {"data15-16", "data17", "data18"};
    if (options.data_period != "all")
        periods = {options.data_period};

    map<string, TChain*> tchains;
    for (auto process : options.processes) {
        tchains[process] = new TChain("BaselineTree");
        cout << process << " filename(s)" << endl;
        for (auto period : periods) {
            if (!options.is_data || ((process != "photon") && (process != "data_bkg"))) period = getMCPeriod(period);
            vector<string> filenames = {ntuple_path + period + "_" + process + ".root"};
            if (process == "photon") {
                if (options.is_data)
                    filenames = {photon_path + period + "_data_photon_ee.root",
                                 photon_path + period + "_data_photon_mm.root"};
                else
                    filenames = {photon_path + getMCPeriod(period) + "_SinglePhoton222_ee.root",
                                 photon_path + getMCPeriod(period) + "_SinglePhoton222_mm.root"};
            }
            else if (process == "Vgamma") {
                filenames = {photon_path + getMCPeriod(period) + "_Vgamma_ee.root",
                             photon_path + getMCPeriod(period) + "_Vgamma_mm.root"};
            }
            for (auto filename : filenames) {
                cout << "\t" << filename << endl;
                tchains[process]->Add(filename.c_str());
            }
        }
    }
    cout << endl;

    //--- add files to RDataFrame
    dataFrameMap RDataFrames;
    for (auto process : options.processes) {
        RDataFrames[process] = new ROOT::RDataFrame(*tchains[process]);
        cout << padString(process + " entries") << *(RDataFrames[process]->Count()) << endl;
    }
    cout << endl;

    return RDataFrames;
}

weightedDataFrameMap weightRDataFrames(dataFrameMap dataframes, Options options) {
    map<string, string> plot_weights;
    for (auto process : options.processes)
        plot_weights[process] = cuts::bkg_weight;
    plot_weights["data_bkg"] = "1";
    plot_weights["photon_raw"] = cuts::photon_weight;
    plot_weights["photon_reweighted"] = cuts::photon_weight * options.reweight_branch.c_str();

    weightedDataFrameMap weighted_dataframes;
    for (auto process : options.processes) {
        if (process == "photon" || process == "Vgamma") {
            auto weighted_dataframe = dataframes[process]->Define("plot_raw_weight", plot_weights["photon_raw"])
                                                .Define("plot_reweighted_weight", plot_weights["photon_reweighted"]);
            using DFType = decltype(weighted_dataframe);
            weighted_dataframes[process] = std::make_unique<DFType>(weighted_dataframe);  
        }
        else {
            auto weighted_dataframe = dataframes[process]->Define("plot_weight", plot_weights[process]);
            using DFType = decltype(weighted_dataframe);
            weighted_dataframes[process] = std::make_unique<DFType>(weighted_dataframe);  
        }
    }
    return weighted_dataframes;
}

filteredDataFrameMap filterRDataFrames(weightedDataFrameMap* WRDataFramesPtr, Options options) {
    filteredDataFrameMap filtered_dataframes;

    for (auto const& [process, weighted_dataframe] : *WRDataFramesPtr) {
        for (string region : options.plot_regions) {
            for (string channel : options.plot_channels) {
                auto [region_name, plot_region, plot_CR] = getPlotRegionInfo(options, channel, region);
                auto filtered_df_plot = weighted_dataframe->Filter(plot_region);
                auto filtered_df_CR = weighted_dataframe->Filter(plot_CR);
                using DFType = decltype(filtered_df_plot);

                for (auto plot_feature : options.plot_features) {
                    filtered_dataframes[region_name][plot_feature]["plot_region"] =
                                        std::make_unique<DFType>(filtered_df_plot);  
                    filtered_dataframes[region_name][plot_feature]["control_region"] =
                                        std::make_unique<DFType>(filtered_df_CR);  
                }
            }
        }
    }

    return filtered_dataframes;
};

//-----------------
// FILL HISTOGRAMS
//-----------------

tuple<histMap, histMap> setUpHistograms(weightedDataFrameMap* WRDataFramesPtr, Options options) {
    cout << PBLU("setting up histograms") << endl;
    cout << endl;

    histMap plot_region_histograms;
    histMap control_region_histograms;

    for (string region : options.plot_regions) {
        for (string channel : options.plot_channels) {
            auto [region_name, plot_region, plot_CR] = getPlotRegionInfo(options, channel, region);
            cout << "\t" << padString("region name") << region_name << endl;
            cout << "\t" << padString("plot region") << plot_region << endl;
            cout << "\t" << padString("normalization reg.") << plot_CR << endl;
            cout << endl;

            for (auto plot_feature : options.plot_features) {
                for (auto const& [process, weighted_dataframe] : *WRDataFramesPtr) {
                    ROOT::RDF::TH1DModel hist_model = getHistogramInfo(plot_feature);
                    if (process == "photon" || process == "Vgamma") {
                        plot_region_histograms[region_name][plot_feature][process+"_raw"] =
                                            weighted_dataframe->Filter(plot_region)
                                            .Histo1D(hist_model, plot_feature, "plot_raw_weight");
                        plot_region_histograms[region_name][plot_feature][process+"_reweighted"] =
                                            weighted_dataframe->Filter(plot_region)
                                            .Histo1D(hist_model, plot_feature, "plot_reweighted_weight");
                        control_region_histograms[region_name][plot_feature][process+"_raw"] =
                                            weighted_dataframe->Filter(plot_CR)
                                            .Histo1D(hist_model, plot_feature, "plot_raw_weight");
                        control_region_histograms[region_name][plot_feature][process+"_reweighted"] =
                                            weighted_dataframe->Filter(plot_CR)
                                            .Histo1D(hist_model, plot_feature, "plot_reweighted_weight");
                    }
                    else {
                        plot_region_histograms[region_name][plot_feature][process] =
                                            weighted_dataframe->Filter(plot_region)
                                            .Histo1D(hist_model, plot_feature, "plot_weight");
                        control_region_histograms[region_name][plot_feature][process] =
                                            weighted_dataframe->Filter(plot_CR)
                                            .Histo1D(hist_model, plot_feature, "plot_weight");
                    }
                }
            }
        }
    }

    return make_tuple(plot_region_histograms, control_region_histograms);
};

resultsMap fillHistograms(tuple<histMap, histMap> region_hists, Options options) {
    cout << PBLU("making histograms") << endl;
    cout << endl;

    auto [plot_region_histograms, control_region_histograms] = region_hists;
    resultsMap results_map;
    results_map.plot_regions = options.plot_regions;
    results_map.plot_features = options.plot_features;
    results_map.plot_channels = options.plot_channels;

    for (string region : options.plot_regions) {
        for (string channel : options.plot_channels) {
            auto [region_name, plot_region, plot_CR] = getPlotRegionInfo(options, channel, region);
            cout << "\t" << padString("region name") << region_name << endl;

            auto prh = plot_region_histograms[region_name];
            auto crh = control_region_histograms[region_name];
            string first_feature = options.plot_features[0];
            auto prh0 = prh[first_feature];
            auto crh0 = crh[first_feature];

            int nbins = getHistogramInfo(first_feature).fNbinsX;
            auto fullInt = [nbins](TH1D* hist) {return hist->Integral(0, nbins+1);};

            //--- Vgamma subtraction
            if (options.do_vgamma_subtraction) {
                //prh0["photon_raw"]->Add(prh0["Vgamma_raw"], -1.0);
                //crh0["photon_raw"]->Add(crh0["Vgamma_raw"], -1.0);
                //prh0["photon_reweighted"]->Add(prh0["Vgamma_reweighted"], -1.0);
                //crh0["photon_reweighted"]->Add(crh0["Vgamma_reweighted"], -1.0);
            }

            //--- scaling
            map<string, float> CR_integrals;
            map<string, float> PR_integrals;
            float mc_bkg_CR_integral = 0.0;
            float mc_bkg_PR_integral = 0.0;
            for (auto process : options.processes) {
                if (process == "photon") {
                    cout << "\t" << padString("photon raw integral");
                    CR_integrals["photon_raw"] = crh0["photon_raw"]->Integral(0, nbins+1);
                    PR_integrals["photon_raw"] = prh0["photon_raw"]->Integral(0, nbins+1);
                    cout << PR_integrals["photon_raw"] << endl;

                    cout << "\t" << padString("photon reweight int.");
                    CR_integrals["photon_reweighted"] = crh0["photon_reweighted"]->Integral(0, nbins+1);
                    PR_integrals["photon_reweighted"] = prh0["photon_reweighted"]->Integral(0, nbins+1);
                    cout << PR_integrals["photon_reweighted"] << endl;
                }
                else {
                    cout << "\t" << padString(process + " integral");
                    CR_integrals[process] = crh0[process]->Integral(0, nbins+1);
                    PR_integrals[process] = prh0[process]->Integral(0, nbins+1);
                    cout << PR_integrals[process] << endl;
                }
                if (process != "data_bkg" && process != "photon" && process != "Zjets") {
                    mc_bkg_CR_integral += CR_integrals[process];
                    mc_bkg_PR_integral += PR_integrals[process];
                }
            }

            float zdata_integral = CR_integrals["data_bkg"] - mc_bkg_CR_integral;
            if (!options.is_data) zdata_integral = CR_integrals["Zjets"];
            float SF = CR_integrals["photon_raw"] == 0 ? 0.0 : zdata_integral / CR_integrals["photon_raw"];
            float SFrw = CR_integrals["photon_reweighted"] == 0 ? 0.0 : zdata_integral / CR_integrals["photon_reweighted"];

            cout << "\tScaling raw photon data by " << SF << endl;
            cout << "\tScaling reweighted photon data by " << SFrw << endl;

            for (auto plot_feature : options.plot_features) {
                prh[plot_feature]["photon_raw"]->Scale(SF);
                prh[plot_feature]["photon_reweighted"]->Scale(SFrw);
            }
            PR_integrals["photon_raw"] *= SF;
            PR_integrals["photon_reweighted"] *= SFrw;

            PR_integrals["bkg MC"] = mc_bkg_PR_integral;
            PR_integrals["photon + bkg MC"] = PR_integrals["photon_reweighted"] + mc_bkg_PR_integral;
            PR_integrals["Z + bkg MC"] = PR_integrals["Zjets"] + mc_bkg_PR_integral;

            cout << "\tScaled photon yield of " << PR_integrals["photon_reweighted"] << endl;
            cout << endl;

            map<string, map<string, TH1D*>> region_hists;
            for (auto plot_feature : options.plot_features) {
                for (auto [process, histogram] : prh[plot_feature]) {
                    TH1D* new_histogram = new TH1D;
                    histogram->Copy(*new_histogram);
                    region_hists[plot_feature][process] = new_histogram;
                }
            }
            results_map.results[region_name] = Result{region_hists, PR_integrals, SFrw};
        }
    }

    return results_map;
}

//-------------
// MAKE TABLES
//-------------

string toString(float val) {
    std::ostringstream out;
    out << std::setprecision(3) << std::fixed << val; // set printouts to 3 sig figs
    return out.str();
}

string getChannelString(map<string, string> channel_values, vector<string> plot_channels) {
    string channel_string = "";
    for (auto channel : plot_channels) {
        if (channel_string.length() > 0) channel_string += " / ";
        channel_string += channel_values[channel];
    }
    return channel_string;
}

void printPhotonYieldTables(Options options, resultsMap results_map, string save_name, bool blinded) {
    //--- [region_name][feature], with a dictionary of hists by process and a photon yield value
    //--- region_name can further be split into [region][channel]
    ofstream out_file;
    out_file.open(save_name);

    out_file << "\\documentclass{article}" << endl;
    out_file << endl;
    out_file << "\\usepackage[utf8]{inputenc}" << endl;
    out_file << "\\usepackage{color, colortbl}" << endl;
    out_file << endl;
    out_file << "\\begin{document}" << endl;
    out_file << endl;
    out_file << "\\definecolor{Gray}{gray}{0.9}" << endl;
    out_file << "\\newcolumntype{g}{>{\\columncolor{Gray}}c}" << endl;
    out_file << endl;
    out_file << "\\begin{table}" << endl;
    map<string, string> plot_channels = {{"ee", "ee"}, {"mm", "mm"}, {"SF", "SF"}};
    string channel_string = getChannelString(plot_channels, options.plot_channels);
    out_file << "\\caption{Photon Method Yields (" << channel_string << ")}" << endl;
    out_file << "\\begin{center}" << endl;
    out_file << "\\begin{tabular}{c|gcccgg}" << endl;
    vector<string> table_processes = {"data_bkg", "bkg MC", "photon", "Zjets", "photon + bkg MC", "Z + bkg MC"};
    out_file << "region & data & bkg$_{MC}$ & photon & Zjets & photon$_{tot}$ & Zjets$_{tot}$ \\\\" << endl;
    out_file << "\\hline" << endl;

    for (auto region : results_map.plot_regions) {
        out_file << region;
        for (auto process : table_processes) {
            if (process == "photon") process = "photon_reweighted";
            float yield_ee = results_map.results[region + " ee"].process_yields[process];
            float yield_mm = results_map.results[region + " mm"].process_yields[process];
            float yield_SF = results_map.results[region + " SF"].process_yields[process];
            map<string, string> channel_yields = {
                {"ee", toString(yield_ee)},
                {"mm", toString(yield_mm)},
                {"SF", toString(yield_SF)}
            };
            if ((process == "data_bkg") && blinded && (region.find("SR") != std::string::npos))
                channel_yields = {{"ee", "-"}, {"mm", "-"}, {"SF", "-"}};
            out_file << " & " << getChannelString(channel_yields, options.plot_channels);
        }
        out_file << " \\\\" << endl;
    }

    out_file << "\\end{tabular}" << endl;
    out_file << "\\end{center}" << endl;
    out_file << "\\end{table}" << endl;
    out_file << endl;
    out_file << "\\end{document}" << endl;

    out_file.close();
}

void printPhotonScaleFactorTables(Options options, resultsMap results_map, string save_name) {
    //--- [region_name][feature], with a dictionary of hists by process and a photon yield value
    //--- region_name can further be split into [region][channel]
    ofstream out_file;
    out_file.open(save_name);

    out_file << "\\documentclass{article}" << endl;
    out_file << "\\usepackage[utf8]{inputenc}" << endl;
    out_file << endl;
    out_file << "\\begin{document}" << endl;
    out_file << endl;
    out_file << "\\begin{table}" << endl;
    map<string, string> plot_channels = {{"ee", "ee"}, {"mm", "mm"}, {"SF", "SF"}};
    string channel_string = getChannelString(plot_channels, results_map.plot_channels);
    out_file << "\\caption{Photon Method Scale Factors (" << channel_string << ")}" << endl;
    out_file << "\\begin{center}" << endl;
    out_file << "\\begin{tabular}{c|c}" << endl;
    out_file << "region & photon scale factor \\\\" << endl;
    out_file << "\\hline" << endl;

    for (auto region : results_map.plot_regions) {
        float photon_ee_sf = results_map.results[region + " ee"].photon_SF;
        float photon_mm_sf = results_map.results[region + " mm"].photon_SF;
        float photon_SF_sf = results_map.results[region + " SF"].photon_SF;
        map<string, string> channel_photons = {{"ee", toString(photon_ee_sf)}, {"mm", toString(photon_mm_sf)}, {"SF", toString(photon_SF_sf)}};
        out_file << region << " & " << getChannelString(channel_photons, results_map.plot_channels) << " \\\\" << endl;
    }

    out_file << "\\end{tabular}" << endl;
    out_file << "\\end{center}" << endl;
    out_file << "\\end{table}" << endl;
    out_file << endl;
    out_file << "\\end{document}" << endl;

    out_file.close();
}

//------------
// MAKE PLOTS
//------------

tuple<THStack*, THStack*, THStack*> createStacks(map<string, TH1D*> histograms, TString formatted_feature, Options options) {
    THStack *data_stack = new THStack("data_stack", "");
    THStack *raw_g_stack = new THStack("raw_g_stack", "");
    THStack *reweight_g_stack = new THStack("reweight_g_stack", "");

    for (auto process_ptr = options.processes.rbegin(); process_ptr != options.processes.rend(); ++process_ptr) {
        auto process = *process_ptr;

        //--- set plotting options
        if (process == "photon") {
            histograms["photon_reweighted"]->SetLineColor(1);
            histograms["photon_reweighted"]->SetLineWidth(0);
            histograms["photon_reweighted"]->SetFillColor(options.process_colors["photon_reweighted"]);
            if (options.is_data) {
                histograms["photon_raw"]->SetLineColor(kRed);
                histograms["photon_raw"]->SetLineStyle(9);
                histograms["photon_raw"]->SetLineWidth(3);
                histograms["photon_raw"]->SetFillStyle(0);
            }
            else {
                histograms["photon_raw"]->SetLineColor(4);
                histograms["photon_raw"]->SetLineStyle(2);
                histograms["photon_raw"]->SetFillStyle(0);
            }
        }
        else if (process=="Zjets" && !options.is_data) {
            histograms[process]->SetFillColor(42);
            histograms[process]->SetLineStyle(1);
        }
        else if (process == "data_bkg")
            histograms[process]->SetMarkerStyle(20);
        else {
            histograms[process]->SetLineColor(1);
            histograms[process]->SetLineWidth(0);
            histograms[process]->SetFillColor(options.process_colors[process]);
        }

        //--- turn on overflow bin
        if (process == "photon") {
            histograms["photon_raw"]->GetXaxis()->SetRange(0, histograms["photon_raw"]->GetNbinsX() + 1);
            histograms["photon_reweighted"]->GetXaxis()->SetRange(0, histograms["photon_reweighted"]->GetNbinsX() + 1);
        }
        else
            histograms[process]->GetXaxis()->SetRange(0, histograms[process]->GetNbinsX() + 1);

        //--- add to relevant stack
        if (process == "data_bkg") {
            if (options.is_data) data_stack->Add(histograms[process]);
        }
        else if (process == "Zjets") {
            if (!options.is_data) data_stack->Add(histograms[process]);
        }
        else if (process == "photon") {
            raw_g_stack->Add(histograms["photon_raw"]);
            reweight_g_stack->Add(histograms["photon_reweighted"]);
        }
        else if (options.is_data) {
            raw_g_stack->Add(histograms[process]);
            reweight_g_stack->Add(histograms[process]);
        }
    }

    return make_tuple(data_stack, raw_g_stack, reweight_g_stack);
}

TString getPlotSaveName(string period, string channel, string plot_feature, string reweight_branch, bool is_data, string region, string plots_path) {
    TString plot_name;
    if (is_data)
        plot_name = Form("%s/%s_%s_%s_%s_%s", plots_path.c_str(), reweight_branch.c_str(), period.c_str(), channel.c_str(), plot_feature.c_str(), region.c_str());
    else
        plot_name = Form("%s/%s_%s_%s_%s_%s", plots_path.c_str(), reweight_branch.c_str(), getMCPeriod(period).c_str(), channel.c_str(), plot_feature.c_str(), region.c_str());
    plot_name += ".eps";
    return plot_name;
}

TLegend* getLegend(Options options, map<string, TH1D*> histograms) {
    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
    if (options.is_data) {
        leg->AddEntry(histograms["data_bkg"], options.process_latex["data_bkg"].c_str(), "lp");
        for (auto process : options.processes) {
            if ((process != "data_bkg") && (process != "Zjets")) {
                if (process == "photon") {
                    leg->AddEntry(histograms["photon_reweighted"],
                        options.process_latex["photon_reweighted"].c_str(), "f");
                    if (options.plot_unreweighted_photons)
                        leg->AddEntry(histograms["photon_raw"], options.process_latex["photon_raw"].c_str(), "f");
                }
                else
                    leg->AddEntry(histograms[process], options.process_latex[process].c_str(), "f");
            }
        }
    }
    else {
        leg->AddEntry(histograms["Zjets"], options.process_latex["Zjets"].c_str(), "f");
        if (options.plot_unreweighted_photons)
            leg->AddEntry(histograms["photon_raw"], options.process_latex["photon_raw"].c_str(), "f");
        leg->AddEntry(histograms["photon_reweighted"], options.process_latex["photon_reweighted"].c_str(), "f");
    }

    leg->SetBorderSize(0);
    leg->SetFillColor(0);

    return leg;
}

string getPlotTex(Options options) {
    string tex_string;
    if (TString(options.period).Contains("all")) tex_string = "139 fb^{-1} 2015-2018 data";
    if (options.is_data) {
        if(TString(options.period).Contains("data15-16")) tex_string = "36 fb^{-1} 2015-2016 data";
        if(TString(options.period).Contains("data17")) tex_string = "44 fb^{-1} 2017 data";
        if(TString(options.period).Contains("data18")) tex_string = "60 fb^{-1} 2018 data";
    }
    else {
        string mc_period = getMCPeriod(options.period);
        if(TString(mc_period).Contains("mc16a")) tex_string = "MC16a";
        if(TString(mc_period).Contains("mc16cd")) tex_string = "MC16cd";
        if(TString(mc_period).Contains("mc16e")) tex_string = "MC16e";
    }

    return tex_string;
}

tuple<TH1D*, TH1D*> getRatioPlots(Options options, map<string, TH1D*> histograms) {
    TH1D *hratio, *hratio_unreweighted, *hmctot, *hmctot_unreweighted;

    if (options.is_data) {
        hratio = (TH1D*) histograms["data_bkg"]->Clone("hratio");
        hratio_unreweighted = (TH1D*) histograms["data_bkg"]->Clone("hratio");
        hmctot = (TH1D*) histograms["photon_reweighted"]->Clone("hmctot");
        hmctot_unreweighted = (TH1D*) histograms["photon_raw"]->Clone("hmctot");
        for (auto process : options.processes) {
            if ((process != "data_bkg") && (process != "Zjets") && (process != "photon")) {
                hmctot->Add(histograms[process]);
                hmctot_unreweighted->Add(histograms[process]);
            }
        }
    }
    else {
        hratio = (TH1D*) histograms["Zjets"]->Clone("hratio");
        hratio_unreweighted = (TH1D*) histograms["Zjets"]->Clone("hratio");
        hmctot = (TH1D*) histograms["photon_reweighted"]->Clone("hmctot");
        hmctot_unreweighted = (TH1D*) histograms["photon_raw"]->Clone("hmctot");
    }

    for (int ibin=1; ibin <= hmctot->GetXaxis()->GetNbins(); ibin++) {
        hmctot->SetBinError(ibin, 0.0);
        hmctot_unreweighted->SetBinError(ibin, 0.0);
    }

    hratio_unreweighted->Divide(hmctot_unreweighted);
    hratio_unreweighted->SetMarkerStyle(20);
    auto unreweighted_color = histograms["photon_raw"]->GetLineColor();
    hratio_unreweighted->SetMarkerColor(unreweighted_color);
    hratio_unreweighted->SetLineColor(unreweighted_color);
    hratio_unreweighted->GetXaxis()->SetTitle("");
    hratio_unreweighted->GetXaxis()->SetLabelSize(0.);
    hratio_unreweighted->GetYaxis()->SetNdivisions(5);
    hratio_unreweighted->GetYaxis()->SetTitle("");
    hratio_unreweighted->GetYaxis()->SetTitleSize(0.15);
    hratio_unreweighted->GetYaxis()->SetTitleOffset(0.3);
    hratio_unreweighted->GetYaxis()->SetLabelSize(0.15);
    hratio_unreweighted->SetMinimum(0.0);
    hratio_unreweighted->SetMaximum(2.0);
    hratio_unreweighted->GetYaxis()->SetRangeUser(0.0,2.0);

    hratio->Divide(hmctot);
    hratio->SetMarkerStyle(20);
    auto reweighted_color = histograms["photon_reweighted"]->GetLineColor();
    hratio->SetMarkerColor(reweighted_color);
    hratio->SetLineColor(reweighted_color);
    hratio->GetXaxis()->SetTitle("");
    hratio->GetXaxis()->SetLabelSize(0.);
    hratio->GetYaxis()->SetNdivisions(5);
    if (options.is_data)
        hratio->GetYaxis()->SetTitle("data/bkg");
    else
        hratio->GetYaxis()->SetTitle("Z/#gamma MC");
    hratio->GetYaxis()->SetTitleSize(0.15);
    hratio->GetYaxis()->SetTitleOffset(0.3);
    hratio->GetYaxis()->SetLabelSize(0.15);
    hratio->SetMinimum(0.0);
    hratio->SetMaximum(2.0);
    hratio->GetYaxis()->SetRangeUser(0.0,2.0);

    hratio->SetTitle("");
    hratio_unreweighted->SetTitle("");

    return make_tuple(hratio, hratio_unreweighted);
}

void makePlot(resultsMap results_map, Options options) {
    for (auto region : results_map.plot_regions) {
        for (auto channel : options.plot_channels) {
            string region_name = region + " " + channel;
            for (auto feature : results_map.plot_features) {
                //--- draw title
                TCanvas *can = new TCanvas("can","can",600,600);
                can->cd();
                TPad* namepad = new TPad("namepad","namepad",0.0,0.0,1.0,1.0);
                namepad->Draw();
                namepad->cd();
                ROOT::RDF::TH1DModel plot_info = getHistogramInfo(feature);
                TString formatted_feature = plot_info.fTitle;
                TString plot_title = formatted_feature + " in " + region.c_str();
                //TH1D *h_name = new TH1D("h_name", plot_title, nbins, xmin, xmax);
                TH1D *h_name = new TH1D("h_name", plot_title, 1, 0, 1);
                h_name->Draw();

                //--- create comparison stacks
                auto hist_map = results_map.results[region_name].histograms[feature]; // map of histograms by [process]
                auto [data_stack, raw_g_stack, reweight_g_stack] = createStacks(hist_map, formatted_feature, options);

                //--- draw plot
                TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);
                mainpad->Draw();
                mainpad->cd();

                float max_y = max(reweight_g_stack->GetMaximum(), data_stack->GetMaximum()) * 1.5;
                float min_y = 0;
                vector<string> log_features = {"met_Et", "METl", "METt", "met_Sign", "Ptll", "Ht30"};
                if (find(log_features.begin(), log_features.end(), feature) != log_features.end()) {
                    mainpad->SetLogy();
                    float max_y = max(reweight_g_stack->GetMaximum(), data_stack->GetMaximum()) * 50;
                    min_y = pow(10.0, -2);
                }

                bool applicable_blinded = (options.is_data)
                                          && (options.blinded && (region.find("SR") != std::string::npos));

                if (options.is_data) {
                    reweight_g_stack->Draw("hist");
                    reweight_g_stack->SetMaximum(max_y);
                    reweight_g_stack->SetMinimum(min_y);
                    if (options.plot_unreweighted_photons)
                        raw_g_stack->GetStack()->Last()->Draw("samehist");
                    if (!applicable_blinded)
                        data_stack->Draw("sameE1");
                    reweight_g_stack->GetXaxis()->SetTitle(formatted_feature);
                    reweight_g_stack->GetYaxis()->SetTitle("entries / bin");
                }
                else {
                    data_stack->Draw("hist");
                    data_stack->SetMaximum(max_y);
                    data_stack->SetMinimum(min_y);
                    if (options.plot_unreweighted_photons)
                        raw_g_stack->GetStack()->Last()->Draw("samehist");
                    reweight_g_stack->Draw("samehist");
                    data_stack->GetXaxis()->SetTitle(formatted_feature);
                    data_stack->GetYaxis()->SetTitle("entries / bin");
                }

                //--- draw legend and labels
                TLegend *leg = getLegend(options, hist_map);
                leg->Draw();

                //--- write info
                TLatex *tex = new TLatex();
                tex->SetNDC();
                tex->SetTextSize(0.025);
                tex->DrawLatex(0.6,0.65,"ATLAS Internal");
                tex->DrawLatex(0.6,0.61,getPlotTex(options).c_str());
                if(TString(channel).Contains("ee")) tex->DrawLatex(0.6,0.57,"ee events");
                if(TString(channel).Contains("mm")) tex->DrawLatex(0.6,0.57,"#mu#mu events");
                if(TString(channel).Contains("em")) tex->DrawLatex(0.6,0.57,"e#mu events");
                if(TString(channel).Contains("me")) tex->DrawLatex(0.6,0.57,"#mu e events");
                if(TString(channel).Contains("SF")) tex->DrawLatex(0.6,0.57,"SF events");
                if(TString(channel).Contains("DF")) tex->DrawLatex(0.6,0.57,"DF events");

                //--- draw ratio
                can->cd();
                TPad* ratio_pad = new TPad("ratio_pad","ratio_pad",0.0,0.75,1.0,0.905);
                ratio_pad->Draw();
                ratio_pad->cd();
                ratio_pad->SetGridy();

                auto [hratio, hratio_unreweighted] = getRatioPlots(options, hist_map);
                if (applicable_blinded) {
                    TH1D *empty_hist = new TH1D("", "", 1, 0, 1);
                    empty_hist->Draw();
                    tex->SetTextSize(0.3);
                    tex->DrawLatex(0.42,0.42,"BLINDED");
                }
                else {
                    if (options.plot_unreweighted_photons)
                        hratio_unreweighted->Draw("E1");
                    hratio->Draw("sameE1");
                }

                //--- save plot
                TString plot_name = getPlotSaveName(options.data_period, channel, feature, options.reweight_branch, options.is_data, region, options.plots_folder);
                can->Print(plot_name);

                //--- clean up
                delete can;
                delete h_name;
            }
        }
    }
}

//----------------
// TEST FUNCTIONS
//----------------

void run_quickDraw(Options options);

void performPlottingUnitTests(Options options) {
    cout << BOLD(PBLU("Performing unit testing on plotting step")) << endl;
    cout << endl;

    options.data_period = "data15-16";
    options.is_data = true;
    options.reduction_folder = options.unit_test_folder + "ReducedNtuples/";
    options.reweighting_folder = options.unit_test_folder + "ReweightedNtuples/";
    options.plots_folder = "Diagnostics/Plots/";

    options.plot_regions = vector<string>{"VRDPhiLow6"};
    options.plot_features = vector<string>{"mll", "Ptll", "met_Et", "met_Sign", "mt2leplsp_0", "Ht30"};

    run_quickDraw(options);

    passTest("Check output folders for sample tables and plots");
}

//---------------
// MAIN FUNCTION
//---------------

void run_quickDraw(Options options) {
    cout << BOLD(PBLU("Making plots")) << endl;
    cout << endl;

    cout << padString("period") << options.data_period << endl;
    cout << padString("is data?") << options.is_data << endl;
    cout << endl;

    //--- get input data in the form of RDataFrames
    dataFrameMap RDataFrames = getRDataFrames(options);
    weightedDataFrameMap WRDataFrames = weightRDataFrames(RDataFrames, options);

    //--- set up all histograms and fill simulaneously
    auto region_hists = setUpHistograms(&WRDataFrames, options);
    resultsMap results_map = fillHistograms(region_hists, options);

    //--- print photon yield tables
    TString plot_name = getPlotSaveName(options.data_period, "yields", "allFeatures", options.reweight_branch,
                            options.is_data, "allRegions", options.plots_folder);
    string type = options.is_data ? "Data" : "MC";
    printPhotonYieldTables(options, results_map, options.plots_folder + options.reweight_branch + "_" +
            options.data_period + "_" + type + "_yields.txt", options.blinded);
    printPhotonScaleFactorTables(options, results_map, options.plots_folder + options.reweight_branch + "_" +
            options.data_period + "_" + type + "_scale_factors.txt");

    //--- draw and save plot
    if (!options.print_photon_yield_only)
        makePlot(results_map, options);
}

void MakePlots(Options options) {
    //--- set global options
    //ROOT::EnableImplicitMT();
    gStyle->SetOptStat(0);

    //--- either perform unit tests or run code
    if (options.unit_testing)
        performPlottingUnitTests(options);
    else
        run_quickDraw(options);
}

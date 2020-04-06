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
    vector<string> regions;
    vector<string> features;
    vector<string> channels;
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

tuple<string, string, string> getPlotRegionInfo(string channel, string region) {
    if (cuts::plot_regions.count(region) == 0) {
        cout << "Unrecognized region! Exiting." << endl;
        exit(0);
    }
    TCut plot_region = cuts::plot_regions[region];

    if (TString(channel).EqualTo("ee")) plot_region += cuts::ee;
    else if (TString(channel).EqualTo("mm")) plot_region += cuts::mm;
    else if (TString(channel).EqualTo("em")) plot_region += cuts::em;
    else if (TString(channel).EqualTo("me")) plot_region += cuts::me;
    else if (TString(channel).EqualTo("SF")) plot_region += cuts::SF;
    else if (TString(channel).EqualTo("DF")) plot_region += cuts::DF;
    else {
        cout << "Unrecognized channel " << channel << "! Exiting." << endl;
        exit(0);
    }

    TCut plot_CR = plot_region + cuts::CR;
    plot_region += cuts::plot_region_met_portions[region];

    string region_name = region + " " + channel;

    return make_tuple(region_name, string(plot_region), string(plot_CR));
}

ROOT::RDF::TH1DModel getHistogramInfo(string plot_feature) {
    map<string, ROOT::RDF::TH1DModel> plot_settings;
    plot_settings["met_Et"] = ROOT::RDF::TH1DModel("", "E_{T}^{miss} [GeV]", 50, 0, 1000);
    plot_settings["METl"] = ROOT::RDF::TH1DModel("", "E_{T,||}^{miss} [GeV]", 30, -150, 150);
    plot_settings["METt"] = ROOT::RDF::TH1DModel("", "E_{T,#perp}^{miss} [GeV]", 30, -150, 150);
    plot_settings["met_Sign"] = ROOT::RDF::TH1DModel("", "E_{T}^{miss} significance", 50, 0, 50);
    plot_settings["MET_loose"] = ROOT::RDF::TH1DModel("", "E_{T,loose}^{miss} [GeV]", 20, 0, 200);
    plot_settings["MET_tight"] = ROOT::RDF::TH1DModel("", "E_{T,tight}^{miss} [GeV]", 20, 0, 200);
    plot_settings["MET_tighter"] = ROOT::RDF::TH1DModel("", "E_{T,tighter}^{miss} [GeV]", 20, 0, 200);
    plot_settings["MET_tenacious"] = ROOT::RDF::TH1DModel("", "E_{T,tenacious}^{miss} [GeV]", 20, 0, 200);
    plot_settings["mt2leplsp_0"] = ROOT::RDF::TH1DModel("", "m_{T2}^{0} [GeV]", 50, 0, 500);
    plot_settings["Ptll"] = ROOT::RDF::TH1DModel("", "p_{T} [GeV]", 25, 0, 1000);
    plot_settings["Z_pt"] = ROOT::RDF::TH1DModel("", "p_{T} [GeV]", 20, 0, 100);
    plot_settings["nJet30"] = ROOT::RDF::TH1DModel("", "n_{jets}", 6, 2, 8);
    plot_settings["jet_n"] = ROOT::RDF::TH1DModel("", "n_{jets}", 6, 2, 8);
    plot_settings["jet_eta"] = ROOT::RDF::TH1DModel("", "jet_{#eta}", 30, -3, 3);
    plot_settings["jet_phi"] = ROOT::RDF::TH1DModel("", "jet_{#phi}", 20, 0, 3.14);
    plot_settings["jetPt"] = ROOT::RDF::TH1DModel("", "jet_{p_{T}} [GeV]", 20, 0, 300);
    plot_settings["bjet_n"] = ROOT::RDF::TH1DModel("", "n_{b-jets}", 4, 0, 4);
    plot_settings["Ht30"] = ROOT::RDF::TH1DModel("", "H_{T}", 30, 0, 1500);
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

    ROOT::RDF::TH1DModel hist_model = plot_settings[plot_feature];
    return hist_model;
}

//--------------------
// SET UP RDATAFRAMES
//--------------------

dataFrameMap getRDataFrames(PlottingOptions options) {
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
        cout << padString(process + " entries") << ": " << *(RDataFrames[process]->Count()) << endl;
    }
    cout << endl;

    return RDataFrames;
}

weightedDataFrameMap weightRDataFrames(dataFrameMap dataframes, PlottingOptions options) {
    map<string, string> plot_weights;
    for (auto process : options.processes)
        plot_weights[process] = cuts::bkg_weight;
    plot_weights["data_bkg"] = "1";
    plot_weights["photon_raw"] = cuts::photon_weight;
    plot_weights["photon_reweighted"] = cuts::photon_weight_rw;

    weightedDataFrameMap weighted_dataframes;
    for (auto [process, dataframe] : dataframes) {
        if (process == "photon") {
            auto weighted_dataframe = dataframe->Define("plot_raw_weight", plot_weights["photon_raw"])
                                                .Define("plot_reweighted_weight", plot_weights["photon_reweighted"]);
            using DFType = decltype(weighted_dataframe);
            weighted_dataframes[process] = std::make_unique<DFType>(weighted_dataframe);  
        }
        else {
            auto weighted_dataframe = dataframe->Define("plot_weight", plot_weights[process]);
            using DFType = decltype(weighted_dataframe);
            weighted_dataframes[process] = std::make_unique<DFType>(weighted_dataframe);  
        }
    }
    return weighted_dataframes;
}

filteredDataFrameMap filterRDataFrames(weightedDataFrameMap* WRDataFramesPtr, PlottingOptions options) {
    filteredDataFrameMap filtered_dataframes;

    for (auto const& [process, weighted_dataframe] : *WRDataFramesPtr) {
        for (string region : options.regions) {
            for (string channel : options.channels) {
                auto [region_name, plot_region, plot_CR] = getPlotRegionInfo(channel, region);
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

tuple<histMap, histMap> setUpHistograms(weightedDataFrameMap* WRDataFramesPtr, PlottingOptions options) {
    cout << PBLU("setting up histograms") << endl;
    cout << endl;

    histMap plot_region_histograms;
    histMap control_region_histograms;

    for (string region : options.regions) {
        for (string channel : options.channels) {
            auto [region_name, plot_region, plot_CR] = getPlotRegionInfo(channel, region);
            cout << "\t" << padString("region name") << ": " << region_name << endl;
            cout << "\t" << padString("plot region") << ": " << plot_region << endl;
            cout << "\t" << padString("normalization reg.") << ": " << plot_CR << endl;
            cout << endl;

            for (auto plot_feature : options.plot_features) {
                for (auto const& [process, weighted_dataframe] : *WRDataFramesPtr) {
                    ROOT::RDF::TH1DModel hist_model = getHistogramInfo(plot_feature);
                    if (process == "photon") {
                        plot_region_histograms[region_name][plot_feature]["photon_raw"] =
                                            weighted_dataframe->Filter(plot_region)
                                            .Histo1D(hist_model, plot_feature, "plot_raw_weight");
                        plot_region_histograms[region_name][plot_feature]["photon_reweighted"] =
                                            weighted_dataframe->Filter(plot_region)
                                            .Histo1D(hist_model, plot_feature, "plot_reweighted_weight");
                        control_region_histograms[region_name][plot_feature]["photon_raw"] =
                                            weighted_dataframe->Filter(plot_CR)
                                            .Histo1D(hist_model, plot_feature, "plot_raw_weight");
                        control_region_histograms[region_name][plot_feature]["photon_reweighted"] =
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

resultsMap fillHistograms(tuple<histMap, histMap> region_hists, PlottingOptions options) {
    cout << PBLU("making histograms") << endl;
    cout << endl;

    auto [plot_region_histograms, control_region_histograms] = region_hists;
    resultsMap results_map;
    results_map.regions = options.regions;
    results_map.features = options.plot_features;
    results_map.channels = options.channels;

    for (string region : options.regions) {
        for (string channel : options.channels) {
            auto [region_name, plot_region, plot_CR] = getPlotRegionInfo(channel, region);
            cout << "\t" << padString("region name") << ": " << region_name << endl;

            auto prh = plot_region_histograms[region_name];
            auto crh = control_region_histograms[region_name];
            string first_feature = options.plot_features[0];
            auto prh0 = prh[first_feature];
            auto crh0 = crh[first_feature];

            int nbins = getHistogramInfo(first_feature).fNbinsX;
            auto fullInt = [nbins](TH1D* hist) {return hist->Integral(0, nbins+1);};

            map<string, float> CR_integrals;
            map<string, float> PR_integrals;
            float mc_bkg_CR_integral = 0.0;
            float mc_bkg_PR_integral = 0.0;
            for (auto process : options.processes) {
                if (process == "photon") {
                    cout << "\t" << padString("photon raw integral") << ": ";
                    CR_integrals["photon_raw"] = crh0["photon_raw"]->Integral(0, nbins+1);
                    PR_integrals["photon_raw"] = prh0["photon_raw"]->Integral(0, nbins+1);
                    cout << PR_integrals["photon_raw"] << endl;

                    cout << "\t" << padString("photon reweight int.") << ": ";
                    CR_integrals["photon_reweighted"] = crh0["photon_reweighted"]->Integral(0, nbins+1);
                    PR_integrals["photon_reweighted"] = prh0["photon_reweighted"]->Integral(0, nbins+1);
                    cout << PR_integrals["photon_reweighted"] << endl;
                }
                else {
                    cout << "\t" << padString(process + " integral") << ": ";
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

            cout << "\tPhoton yield of " << PR_integrals["photon_reweighted"] << endl;
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

string getChannelString(map<string, string> channel_values, vector<string> channels) {
    string channel_string = "";
    for (auto channel : channels) {
        if (channel_string.length() > 0) channel_string += " / ";
        channel_string += channel_values[channel];
    }
    return channel_string;
}

void printPhotonYieldTables(PlottingOptions options, resultsMap results_map, string save_name, bool blinded) {
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
    map<string, string> channels = {{"ee", "ee"}, {"mm", "mm"}, {"SF", "SF"}};
    string channel_string = getChannelString(channels, options.channels);
    out_file << "\\caption{Photon Method Yields (" << channel_string << ")}" << endl;
    out_file << "\\begin{center}" << endl;
    out_file << "\\begin{tabular}{c";
    for (auto process : options.processes)
        out_file << "|c";
    out_file << "}" << endl;
    out_file << "region";
    for (auto process : options.processes)
        out_file << " & " << process;
    out_file << " \\\\" << endl;
    out_file << "\\hline" << endl;

    for (auto region : results_map.regions) {
        out_file << region;
        for (auto process : options.processes) {
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
            out_file << " & " << getChannelString(channel_yields, options.channels);
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

void printPhotonScaleFactorTables(PlottingOptions options, resultsMap results_map, string save_name) {
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
    out_file << "\\caption{Photon Method Scale Factors (" << channel_string << ")}" << endl;
    out_file << "\\begin{center}" << endl;
    out_file << "\\begin{tabular}{c|c}" << endl;
    map<string, string> channels = {{"ee", "ee"}, {"mm", "mm"}, {"SF", "SF"}};
    string channel_string = getChannelString(channels, results_map.channels);
    out_file << "region & photon scale factor \\\\" << endl;
    out_file << "\\hline" << endl;

    for (auto region : results_map.regions) {
        float photon_ee_sf = results_map.results[region + " ee"].photon_SF;
        float photon_mm_sf = results_map.results[region + " mm"].photon_SF;
        float photon_SF_sf = results_map.results[region + " SF"].photon_SF;
        map<string, string> channel_photons = {{"ee", toString(photon_ee_sf)}, {"mm", toString(photon_mm_sf)}, {"SF", toString(photon_SF_sf)}};
        out_file << region << " & " << getChannelString(channel_photons, results_map.channels) << " \\\\" << endl;
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

tuple<THStack*, THStack*, THStack*> createStacks(map<string, TH1D*> histograms, TString formatted_feature, PlottingOptions options) {
    //--- set plotting options
    vector<int> colors = {kGreen-5, kGreen+2, kRed+1, kSpring+10, kCyan, kMagenta, kYellow+1, kAzure-9, kGreen-2,
                          kViolet+5};
    int color_count = 0;
    for (auto [process, histogram] : histograms) {
        histograms[process]->SetLineColor(1);
        histograms[process]->SetLineWidth(1);
        if (process != "data_bkg") histograms[process]->SetFillColor(colors[color_count++]);
        if (process == "photon_raw") {
            histograms[process]->SetLineColor(4);
            histograms[process]->SetLineStyle(2);
        }
        else if (process == "photon_reweighted") {
            if (options.is_data)
                histograms[process]->SetFillColor(kOrange-2);
            else {
                histograms[process]->SetLineColor(kRed-2);
                histograms[process]->SetFillStyle(0);
            }
        }
        else if (process == "Zjets") {
            if (options.is_data) {
                histograms[process]->SetLineColor(2);
                histograms[process]->SetLineStyle(7);
            }
            else {
                histograms[process]->SetFillColor(42);
                histograms[process]->SetLineStyle(1);
            }
        }
        else if (process == "data_bkg")
            if (options.is_data) histograms[process]->SetMarkerStyle(20);
    }

    //--- turn on overflow bin
    for (auto [process, histogram] : histograms)
        histograms[process]->GetXaxis()->SetRange(0, histograms[process]->GetNbinsX() + 1);

    //--- make stacks
    THStack *data_stack = new THStack("data_stack", "");
    THStack *raw_g_stack = new THStack("raw_g_stack", "");
    THStack *reweight_g_stack = new THStack("reweight_g_stack", "");

    for (auto [process, histogram] : histograms) {
        if (process == "data_bkg") {
            if (options.is_data) data_stack->Add(histograms[process]);
        }
        else if (process == "Zjets") {
            if (!options.is_data) data_stack->Add(histograms[process]);
        }
        else if (process == "photon_raw") {
            raw_g_stack->Add(histograms[process]);
        }
        else if (process == "photon_reweighted") {
            reweight_g_stack->Add(histograms[process]);
        }
        else if (options.is_data) {
            raw_g_stack->Add(histograms[process]);
            reweight_g_stack->Add(histograms[process]);
        }
    }

    return make_tuple(data_stack, raw_g_stack, reweight_g_stack);
}

TString getPlotSaveName(string period, string channel, string plot_feature, bool is_data, string region, string plots_path) {
    TString plot_name;
    if (is_data)
        plot_name = Form("%s/%s_%s_%s_%s", plots_path.c_str(), period.c_str(), channel.c_str(), plot_feature.c_str(), region.c_str());
    else
        plot_name = Form("%s/%s_%s_%s_%s", plots_path.c_str(), getMCPeriod(period).c_str(), channel.c_str(), plot_feature.c_str(), region.c_str());
    plot_name += ".eps";
    return plot_name;
}

TLegend* getLegend(PlottingOptions options, map<string, TH1D*> histograms) {
    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
    if (options.is_data) {
        for (auto [process, histogram] : histograms) {
            if (process == "data_bkg") {
                leg->AddEntry(histograms[process], options.process_latex[process].c_str(), "lp");
            }
            else if (process != "photon_raw") {
                leg->AddEntry(histograms[process], options.process_latex[process].c_str(), "f");
            }
        }
    }
    else {
        leg->AddEntry(histograms["Zjets"], options.process_latex["Zjets"].c_str(), "f");
        leg->AddEntry(histograms["photon_raw"], options.process_latex["photon_raw"].c_str(), "f");
        leg->AddEntry(histograms["photon_reweighted"], options.process_latex["photon_reweighted"].c_str(), "f");
    }

    leg->SetBorderSize(0);
    leg->SetFillColor(0);

    return leg;
}

string getPlotTex(string period, bool is_data) {
    string tex_string;
    if(TString(period).Contains("all")) tex_string = "139 fb^{-1} 2015-2018 data";
    if (is_data) {
        if(TString(period).Contains("data15-16")) tex_string = "36 fb^{-1} 2015-2016 data";
        if(TString(period).Contains("data17")) tex_string = "44 fb^{-1} 2017 data";
        if(TString(period).Contains("data18")) tex_string = "60 fb^{-1} 2018 data";
    }
    else {
        string mc_period = getMCPeriod(period);
        if(TString(mc_period).Contains("mc16a")) tex_string = "MC16a";
        if(TString(mc_period).Contains("mc16cd")) tex_string = "MC16cd";
        if(TString(mc_period).Contains("mc16e")) tex_string = "MC16e";
    }

    return tex_string;
}

tuple<TH1D*, TH1D*> getRatioPlots(map<string, TH1D*> histograms, bool is_data) {
    TH1D *hratio, *hratio_unreweighted, *hmctot, *hmctot_unreweighted;

    if (is_data) {
        hratio = (TH1D*) histograms["data_bkg"]->Clone("hratio");
        hratio_unreweighted = (TH1D*) histograms["data_bkg"]->Clone("hratio");
        hmctot = (TH1D*) histograms["photon_reweighted"]->Clone("hmctot");
        hmctot_unreweighted = (TH1D*) histograms["photon_raw"]->Clone("hmctot");
        for (auto [process, histogram] : histograms) {
            if ((process != "data_bkg") && (process != "Zjets") && (process != "photon_raw")
               && (process != "photon_reweighted")) {
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
    hratio_unreweighted->SetMarkerColor(kBlue);
    hratio_unreweighted->SetLineColor(kBlue);
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
    hratio->SetMarkerColor(kRed);
    hratio->SetLineColor(kRed);
    hratio->GetXaxis()->SetTitle("");
    hratio->GetXaxis()->SetLabelSize(0.);
    hratio->GetYaxis()->SetNdivisions(5);
    if (is_data)
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

void makePlot(resultsMap results_map, PlottingOptions options) {
    for (auto region : results_map.regions) {
        vector<string> channels = {"SF"};
        for (auto channel : channels) {
            string region_name = region + " " + channel;
            for (auto feature : results_map.features) {
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
                mainpad->SetLogy();

                bool applicable_blinded = (options.is_data) && (options.blinded && (region.find("SR") != std::string::npos));

                float max_y = pow(10.0, 4);
                float min_y = pow(10.0, -2);
                if (options.is_data) {
                    reweight_g_stack->Draw("hist");
                    reweight_g_stack->SetMaximum(max_y);
                    reweight_g_stack->SetMinimum(min_y);
                    if (!applicable_blinded)
                        data_stack->Draw("sameE1");
                    reweight_g_stack->GetXaxis()->SetTitle(formatted_feature);
                    reweight_g_stack->GetYaxis()->SetTitle("entries / bin");
                }
                else {
                    data_stack->Draw("hist");
                    data_stack->SetMaximum(max_y);
                    data_stack->SetMinimum(min_y);
                    raw_g_stack->Draw("samehist");
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
                tex->DrawLatex(0.6,0.61,getPlotTex(options.data_period, options.is_data).c_str());
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

                auto [hratio, hratio_unreweighted] = getRatioPlots(hist_map, options.is_data);
                if (applicable_blinded) {
                    TH1D *empty_hist = new TH1D("", "", 1, 0, 1);
                    empty_hist->Draw();
                    tex->SetTextSize(0.3);
                    tex->DrawLatex(0.42,0.42,"BLINDED");
                }
                else {
                    if (!options.is_data)
                        hratio_unreweighted->Draw("E1");
                    hratio->Draw("sameE1");
                }

                //--- save plot
                TString plot_name = getPlotSaveName(options.data_period, channel, feature, options.is_data, region, options.plots_folder);
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

void testTablePrintout(PlottingOptions options, resultsMap results_map, bool SF_only=false) {
    if (SF_only)
        printPhotonYieldTables(options, results_map, "Diagnostics/test_yield_table_SF_only_blindeded.txt", true);
    else {
        printPhotonYieldTables(options, results_map, "Diagnostics/test_yield_table_blindeded.txt", true);
        printPhotonYieldTables(options, results_map, "Diagnostics/test_yield_table_unblindeded.txt", false);
        printPhotonScaleFactorTables(options, results_map, "Diagnostics/test_scale_factor_table.txt");
    }
}

void testMakePlot(resultsMap results_map, string plots_folder) {
    PlottingOptions options;
    options.data_period = "all";
    options.is_data = true;
    options.blinded = true;
    options.plots_folder = plots_folder;
    makePlot(results_map, options);

    options.blinded = false;
    makePlot(results_map, options);
}

void run_quickDraw(PlottingOptions options);

void unit_tests(PlottingOptions options) {
    cout << BOLD(PBLU("Performing unit testing on plotting step")) << endl;
    cout << endl;

    options.data_period = "data15-16";
    options.is_data = true;
    options.reduction_folder = options.unit_test_folder + "ReducedNtuples/";
    options.reweighting_folder = options.unit_test_folder + "ReweightedNtuples/";
    options.plots_folder = "Diagnostics/Plots/";
    run_quickDraw(options);

    passTest("Check output folders for sample tables and plots");

    //resultsMap results_map;
    //vector<string> regions{"SR_test1", "SR_test2", "SR_test3", "VR_test1"};
    //vector<string> features{"METl", "METt", "met_Et"};
    //results_map.regions = regions;
    //results_map.features = features;
    //results_map.channels = vector<string>{"ee", "mm", "SF"};

    //map<string, map<string, TH1D*>> test_hists; // [feature][process] histogram
    //map<string, int> n_entries;
    //n_entries["data_bkg"] = 10000;
    //n_entries["ttbar"] = 6000;
    //n_entries["diboson"] = 3000;
    //n_entries["Zjets"] = 1000;
    //n_entries["photon_raw"] = 1000;
    //n_entries["photon_reweighted"] = 1000;
    //for (auto feature : features) {
        //for (auto process : options.processes) {
            //test_hists[feature][process] = new TH1D("", "", 100, -3, 3);
            //test_hists[feature][process]->FillRandom("gaus", n_entries[process]);
        //}
    //}

    //results_map.results["SR_test1 ee"] = Result{test_hists, 1.3, 1.5, 1.8, 0.4, 1};
    //results_map.results["SR_test1 mm"] = Result{test_hists, 1.2, 1.6, 1.9, 0.6, 1};
    //results_map.results["SR_test1 SF"] = Result{test_hists, 2.5, 3.1, 2.8, 0.1, 1};

    //results_map.results["SR_test2 ee"] = Result{test_hists, 3.6, 3.5, 3.9, 0.2, 3.1};
    //results_map.results["SR_test2 mm"] = Result{test_hists, 3.1, 3.2, 3.9, 0.5, 3.1};
    //results_map.results["SR_test2 SF"] = Result{test_hists, 6.7, 6.7, 7.4, 0.3, 3.1};

    //results_map.results["SR_test3 ee"] = Result{test_hists, 12.9, 13.3, 11.8, 0.6, 1.8};
    //results_map.results["SR_test3 mm"] = Result{test_hists, 11.6, 12.8, 12.7, 0.7, 1.3};
    //results_map.results["SR_test3 SF"] = Result{test_hists, 24.5, 26.1, 23.8, 0.6, 1.9};

    //results_map.results["VR_test1 ee"] = Result{test_hists, 112.9, 113.3, 111.4, 0.2 2.5};
    //results_map.results["VR_test1 mm"] = Result{test_hists, 111.6, 112.8, 112.3, 0.3, 5};
    //results_map.results["VR_test1 SF"] = Result{test_hists, 124.5, 126.1, 123.7, 0.5, 1.8};

    //testTablePrintout(options, results_map);
    //results_map.channels = vector<string>{"SF"};
    //testTablePrintout(options, results_map, true);
    //string plot_folder = "Diagnostics/Plotting/";
    //testMakePlot(results_map, plot_folder);

    //passTest("Produced sample yield table");
    //passTest("Produced sample plots");
    //passTest("Passed all unit tests");
    //cout << endl;
}

//---------------
// MAIN FUNCTION
//---------------

void run_quickDraw(PlottingOptions options) {
    cout << BOLD(PBLU("Making plots")) << endl;
    cout << endl;

    cout << padString("period") << ": " << options.data_period << endl;
    cout << padString("is data?") << ": " << options.is_data << endl;
    cout << endl;

    //--- get input data in the form of RDataFrames
    dataFrameMap RDataFrames = getRDataFrames(options);
    weightedDataFrameMap WRDataFrames = weightRDataFrames(RDataFrames, options);

    //--- set up all histograms and fill simulaneously
    auto region_hists = setUpHistograms(&WRDataFrames, options);
    resultsMap results_map = fillHistograms(region_hists, options);

    //--- print photon yield tables
    TString plot_name = getPlotSaveName(options.data_period, "yields", "allFeatures", options.is_data,
                                        "allRegions", options.plots_folder);
    string type = options.is_data ? "Data" : "MC";
    printPhotonYieldTables(options, results_map, "FinalOutputs/" + options.data_period + "_" +
                           type + "_yields.txt", options.blinded);
    printPhotonScaleFactorTables(options, results_map, "FinalOutputs/" + options.data_period + "_" +
                                 type + "_scale_factors.txt");

    //--- draw and save plot
    if (!options.print_photon_yield_only)
        makePlot(results_map, options);
}

void MakePlots(PlottingOptions options) {
    //--- set global options
    //ROOT::EnableImplicitMT();
    gStyle->SetOptStat(0);

    //--- either perform unit tests or run code
    if (options.unit_testing)
        unit_tests(options);
    else
        run_quickDraw(options);
}

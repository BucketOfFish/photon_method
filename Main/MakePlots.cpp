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
    map<string, double> CR_yields;
    map<string, double> CR_uncertainties;
    map<string, double> process_yields;
    map<string, double> uncertainties;
    map<string, double> scale_factors;
    map<string, double> sf_uncertainties;
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

    TCut plot_CR;
    if (options.scaling_method == "MET")
        plot_CR = NMinus1Cut(plot_region, "met_Et") + cuts::CR_MET;
    if (options.scaling_method == "minDPhi2JetsMet")
        plot_CR = NMinus1Cut(plot_region, "minDPhi2JetsMet") + cuts::minDPhi2JetsMet_anti0p4;

    string region_name = region + " " + channel;

    return make_tuple(region_name, string(plot_region), string(plot_CR));
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
    /*
     * Use RDataFrames to set up a bunch of plot region histograms and associated control region histograms.
     * All cuts are applied here, but none of these histograms will be actually filled until the filling function
     * is called.
     */
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
                        if (options.plot_unreweighted_photons) {
                            plot_region_histograms[region_name][plot_feature][process+"_raw"] =
                                            weighted_dataframe->Filter(plot_region)
                                            .Histo1D(hist_model, plot_feature, "plot_raw_weight");
                            control_region_histograms[region_name][plot_feature][process+"_raw"] =
                                            weighted_dataframe->Filter(plot_CR)
                                            .Histo1D(hist_model, plot_feature, "plot_raw_weight");
                        }
                        if (options.plot_reweighted_photons) {
                            plot_region_histograms[region_name][plot_feature][process+"_reweighted"] =
                                            weighted_dataframe->Filter(plot_region)
                                            .Histo1D(hist_model, plot_feature, "plot_reweighted_weight");
                            control_region_histograms[region_name][plot_feature][process+"_reweighted"] =
                                            weighted_dataframe->Filter(plot_CR)
                                            .Histo1D(hist_model, plot_feature, "plot_reweighted_weight");
                        }
                    }
                    else {
                        if (process == "Zjets" && (options.is_data && !options.plot_zmc)) continue;
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
    /*
     * Fill the plot region and control region histograms given in region_hists, then perform data-bkg subtraction
     * in control regions to get scaling factors. Apply scaling factors to the appropriate process.
     */
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

            //--- get yields and uncertainties for each process, and sum totals for bkg MC
            map<string, double> CR_integrals;
            map<string, double> PR_integrals;
            map<string, double> CR_uncertainties;
            map<string, double> uncertainties;
            CR_integrals["bkg MC"] = 0.0;
            CR_uncertainties["bkg MC"] = 0;
            PR_integrals["bkg MC"] = 0.0;
            uncertainties["bkg MC"] = 0;
            for (auto process : options.processes) {
                if (process == "photon") {
                    if (options.plot_unreweighted_photons) {
                        cout << "\t" << padString("photon raw integral");
                        CR_integrals["photon_raw"] = crh0["photon_raw"]
                            ->IntegralAndError(0, nbins+1, CR_uncertainties["photon_raw"]);
                        PR_integrals["photon_raw"] = prh0["photon_raw"]
                            ->IntegralAndError(0, nbins+1, uncertainties["photon_raw"]);
                        cout << PR_integrals["photon_raw"] << " +/- " << uncertainties["photon_raw"] << endl;
                    }
                    if (options.plot_reweighted_photons) {
                        cout << "\t" << padString("photon reweight int.");
                        CR_integrals["photon_reweighted"] = crh0["photon_reweighted"]
                            ->IntegralAndError(0, nbins+1, CR_uncertainties["photon_reweighted"]);
                        PR_integrals["photon_reweighted"] = prh0["photon_reweighted"]
                            ->IntegralAndError(0, nbins+1, uncertainties["photon_reweighted"]);
                        cout << PR_integrals["photon_reweighted"] << " +/- " <<
                            uncertainties["photon_reweighted"] << endl;
                    }
                }
                else {
                    if (process == "Zjets" && (options.is_data && !options.plot_zmc)) continue;
                    cout << "\t" << padString(process + " integral");
                    CR_integrals[process] = crh0[process]->IntegralAndError(0, nbins+1, CR_uncertainties[process]);
                    PR_integrals[process] = prh0[process]->IntegralAndError(0, nbins+1, uncertainties[process]);
                    cout << PR_integrals[process] << " +/- " << uncertainties[process] << endl;
                }
                if (process != "data_bkg" && process != "photon" && process != "Zjets") {
                    CR_integrals["bkg MC"] += CR_integrals[process];
                    PR_integrals["bkg MC"] += PR_integrals[process];
                    uncertainties["bkg MC"] = sqrt(pow(uncertainties["bkg MC"], 2) + pow(uncertainties[process], 2));
                    CR_uncertainties["bkg MC"] =
                        sqrt(pow(CR_uncertainties["bkg MC"], 2) + pow(CR_uncertainties[process], 2));
                }
            }
            cout << endl;

            //--- scale each Z-like process in the CR
            double zdata_integral = CR_integrals["data_bkg"] - CR_integrals["bkg MC"];
            double zdata_uncertainty = sqrt(pow(CR_uncertainties["data_bkg"], 2) +
                pow(CR_uncertainties["bkg MC"], 2));
            if (!options.is_data) {
                zdata_integral = CR_integrals["Zjets"];
                zdata_uncertainty = CR_uncertainties["Zjets"];
            }

            vector<string> scaled_processes;
            if (options.plot_unreweighted_photons) scaled_processes.push_back("photon_raw");
            if (options.plot_reweighted_photons) scaled_processes.push_back("photon_reweighted");
            if (options.plot_zmc) scaled_processes.push_back("Zjets");

            map<string, double> scale_factors;
            map<string, double> sf_uncertainties;
            cout << "\tScaling via " << options.scaling_method << " method" << endl;
            for (auto process : scaled_processes) {
                double SF = CR_integrals[process] == 0 ? 0.0 : zdata_integral / CR_integrals[process];
                if (SF < 0) SF = 0;
                double SF_uncertainty_percentage = SF == 0 ? 0.0 :
                    (zdata_uncertainty / zdata_integral) + (CR_uncertainties[process] / CR_integrals[process]);
                double SF_uncertainty = SF * SF_uncertainty_percentage;
                cout << "\tScaling " << process << " yield by " << SF << endl;
                for (auto plot_feature : options.plot_features) prh[plot_feature][process]->Scale(SF);
                PR_integrals[process] *= SF;
                uncertainties[process] *= SF;
                scale_factors[process] = SF;
                sf_uncertainties[process] = SF_uncertainty;
                cout << "\tScaled " << process << " yield of " << PR_integrals[process] << endl;
            }
            cout << endl;

            //--- calculate process + bkg totals and uncertainties
            CR_integrals["photon + bkg MC"] = CR_integrals["photon_reweighted"] + CR_integrals["bkg MC"];
            CR_uncertainties["photon + bkg MC"] = sqrt(pow(CR_uncertainties["photon_reweighted"], 2)
                + pow(CR_uncertainties["bkg MC"], 2));
            CR_integrals["Z + bkg MC"] = CR_integrals["Zjets"] + CR_integrals["bkg MC"];
            CR_uncertainties["Z + bkg MC"] = sqrt(pow(CR_uncertainties["Zjets"], 2)
                + pow(CR_uncertainties["bkg MC"], 2));

            PR_integrals["photon + bkg MC"] = PR_integrals["photon_reweighted"] + PR_integrals["bkg MC"];
            uncertainties["photon + bkg MC"] = sqrt(pow(uncertainties["photon_reweighted"], 2)
                + pow(uncertainties["bkg MC"], 2));
            PR_integrals["Z + bkg MC"] = PR_integrals["Zjets"] + PR_integrals["bkg MC"];
            uncertainties["Z + bkg MC"] = sqrt(pow(uncertainties["Zjets"], 2) + pow(uncertainties["bkg MC"], 2));

            //--- return histograms along with yields, scale factors, and uncertainties
            map<string, map<string, TH1D*>> region_hists;
            for (auto plot_feature : options.plot_features) {
                for (auto [process, histogram] : prh[plot_feature]) {
                    TH1D* new_histogram = new TH1D;
                    histogram->Copy(*new_histogram);
                    region_hists[plot_feature][process] = new_histogram;
                }
            }
            results_map.results[region_name] = Result{region_hists, CR_integrals, CR_uncertainties, PR_integrals,
                uncertainties, scale_factors, sf_uncertainties};
        }
    }

    return results_map;
}

//-------------
// MAKE TABLES
//-------------

string toString(double val, int places) {
    std::ostringstream out;
    if (places <= 0)
        out << int(val); // just return val
    else
        out << std::setprecision(places) << std::fixed << val; // set printouts to n numbers after the decimal place
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

void printPhotonCRYieldTables(Options options, resultsMap results_map, string save_name, bool blinded) {
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
    //out_file << "\\definecolor{Gray}{gray}{0.9}" << endl;
    //out_file << "\\newcolumntype{g}{>{\\columncolor{Gray}}c}" << endl;
    //out_file << endl;
    out_file << "\\begin{table}" << endl;
    map<string, string> plot_channels = {{"ee", "ee"}, {"mm", "mm"}, {"SF", "SF"}};
    string channel_string = getChannelString(plot_channels, options.plot_channels);
    out_file << "\\caption{CR Yields (" << channel_string << ")}" << endl;
    out_file << "\\begin{center}" << endl;

    vector<string> processes = {};
    if (options.plot_reweighted_photons) processes.push_back("photon");
    if (options.plot_zmc) processes.push_back("Zjets");
    processes.push_back("bkg MC");
    if (options.plot_reweighted_photons) processes.push_back("photon + bkg MC");
    if (options.plot_zmc) processes.push_back("Z + bkg MC");
    processes.push_back("data_bkg");

    if (processes.size() == 4) {
        out_file << "\\begin{tabular}{l|rr|rr}" << endl;
        if (options.plot_reweighted_photons)
            out_file << "region & photon & other bkgs & photon_{tot} & data \\\\" << endl;
        else out_file << "region & Z MC & other bkgs & Z MC_{tot} & data \\\\" << endl;
    }
    else {
        out_file << "\\begin{tabular}{l|rrr|rrr}" << endl;
        out_file << "region & photon & Z MC & other bkgs & photon_{tot} & Z MC_{tot} & data \\\\" << endl;
    }
    out_file << "\\hline" << endl;

    for (auto region : results_map.plot_regions) {
        out_file << boost::replace_all_copy(region, "_", "\\_");
        for (auto process : processes) {
            if (process == "photon") process = "photon_reweighted";
            double yield_ee = results_map.results[region + " ee"].CR_yields[process];
            double yield_mm = results_map.results[region + " mm"].CR_yields[process];
            double yield_SF = results_map.results[region + " SF"].CR_yields[process];
            double uncertainty_ee = results_map.results[region + " ee"].CR_uncertainties[process];
            double uncertainty_mm = results_map.results[region + " mm"].CR_uncertainties[process];
            double uncertainty_SF = results_map.results[region + " SF"].CR_uncertainties[process];
            map<string, string> channel_yields = {
                {"ee", "$" + toString(yield_ee, 1) + "\\pm" + toString(uncertainty_ee, 1) + "$"},
                {"mm", "$" + toString(yield_mm, 1) + "\\pm" + toString(uncertainty_mm, 1) + "$"},
                {"SF", "$" + toString(yield_SF, 1) + "\\pm" + toString(uncertainty_SF, 1) + "$"}
            };
            if (process == "data_bkg")
                channel_yields = {
                    {"ee", toString(yield_ee, 0)},
                    {"mm", toString(yield_mm, 0)},
                    {"SF", toString(yield_SF, 0)}
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
    //out_file << "\\definecolor{Gray}{gray}{0.9}" << endl;
    //out_file << "\\newcolumntype{g}{>{\\columncolor{Gray}}c}" << endl;
    //out_file << endl;
    out_file << "\\begin{table}" << endl;
    map<string, string> plot_channels = {{"ee", "ee"}, {"mm", "mm"}, {"SF", "SF"}};
    string channel_string = getChannelString(plot_channels, options.plot_channels);
    out_file << "\\caption{Yields (" << channel_string << ")}" << endl;
    out_file << "\\begin{center}" << endl;

    vector<string> processes = {};
    if (options.plot_reweighted_photons) processes.push_back("photon");
    if (options.plot_zmc) processes.push_back("Zjets");
    processes.push_back("bkg MC");
    if (options.plot_reweighted_photons) processes.push_back("photon + bkg MC");
    if (options.plot_zmc) processes.push_back("Z + bkg MC");
    processes.push_back("data_bkg");

    if (processes.size() == 4) {
        out_file << "\\begin{tabular}{l|rr|rr|r}" << endl;
        if (options.plot_reweighted_photons)
            out_file << "region & photon & other bkgs & photon_{tot} & data & significance \\\\" << endl;
        else out_file << "region & Z MC & other bkgs & Z MC_{tot} & data & significance \\\\" << endl;
    }
    else {
        out_file << "\\begin{tabular}{l|rrr|rrr|rr}" << endl;
        out_file << "region & photon & Z MC & other bkgs & photon_{tot} & Z MC_{tot} & data & photon significance & Z MC significance \\\\" << endl;
    }
    out_file << "\\hline" << endl;

    for (auto region : results_map.plot_regions) {
        out_file << boost::replace_all_copy(region, "_", "\\_");
        for (auto process : processes) {
            if (process == "photon") process = "photon_reweighted";
            double yield_ee = results_map.results[region + " ee"].process_yields[process];
            double yield_mm = results_map.results[region + " mm"].process_yields[process];
            double yield_SF = results_map.results[region + " SF"].process_yields[process];
            double uncertainty_ee = results_map.results[region + " ee"].uncertainties[process];
            double uncertainty_mm = results_map.results[region + " mm"].uncertainties[process];
            double uncertainty_SF = results_map.results[region + " SF"].uncertainties[process];
            map<string, string> channel_yields = {
                {"ee", "$" + toString(yield_ee, 1) + "\\pm" + toString(uncertainty_ee, 1) + "$"},
                {"mm", "$" + toString(yield_mm, 1) + "\\pm" + toString(uncertainty_mm, 1) + "$"},
                {"SF", "$" + toString(yield_SF, 1) + "\\pm" + toString(uncertainty_SF, 1) + "$"}
            };
            if (process == "data_bkg")
                channel_yields = {
                    {"ee", toString(yield_ee, 0)},
                    {"mm", toString(yield_mm, 0)},
                    {"SF", toString(yield_SF, 0)}
                };
            if ((process == "data_bkg") && blinded && (region.find("SR") != std::string::npos))
                channel_yields = {{"ee", "-"}, {"mm", "-"}, {"SF", "-"}};
            out_file << " & " << getChannelString(channel_yields, options.plot_channels);
        }
        if (options.plot_reweighted_photons) {
            if (blinded && (region.find("SR") != std::string::npos))
                out_file << " & -";
            else
                out_file << " & " << toString((results_map.results[region + " SF"].process_yields["data_bkg"] -
                    results_map.results[region + " SF"].process_yields["photon + bkg MC"]) /
                    sqrt(pow(results_map.results[region + " SF"].uncertainties["data_bkg"], 2) +
                    pow(results_map.results[region + " SF"].uncertainties["photon + bkg MC"], 2)), 3);
        }
        if (options.plot_zmc) {
            if (blinded && (region.find("SR") != std::string::npos))
                out_file << " & -";
            else
                out_file << " & " << toString((results_map.results[region + " SF"].process_yields["data_bkg"] -
                    results_map.results[region + " SF"].process_yields["Z + bkg MC"]) /
                    sqrt(pow(results_map.results[region + " SF"].uncertainties["data_bkg"], 2) +
                    pow(results_map.results[region + " SF"].uncertainties["Z + bkg MC"], 2)), 3);
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
    out_file << "\\caption{Scale Factors (" << channel_string << "), Scaling by "
        << options.scaling_method << "}" << endl;
    out_file << "\\begin{center}" << endl;

    vector<string> processes;    
    if (options.plot_reweighted_photons) processes.push_back("photon_reweighted");
    //if (options.plot_unreweighted_photons) processes.push_back("photon_raw");
    if (options.plot_zmc) processes.push_back("Zjets");
    if (processes.size() == 1) out_file << "\\begin{tabular}{c|c}" << endl;
    else if (processes.size() == 2) out_file << "\\begin{tabular}{c|cc}" << endl;
    out_file << "region";
    for (auto process : processes) {
        if (process == "Zjets") process = "Z MC";
        if (process == "photon_reweighted") process = "photon_{reweighted}";
        if (process == "photon_raw") process = "photon_{raw}";
        out_file << " & " << process;
    }
    out_file << " \\\\" << endl;
    out_file << "\\hline" << endl;

    for (auto region : results_map.plot_regions) {
        out_file << boost::replace_all_copy(region, "_", "\\_");
        for (auto process : processes) {
            double ee_sf = results_map.results[region + " ee"].scale_factors[process];
            double mm_sf = results_map.results[region + " mm"].scale_factors[process];
            double SF_sf = results_map.results[region + " SF"].scale_factors[process];
            double ee_sf_unc = results_map.results[region + " ee"].sf_uncertainties[process];
            double mm_sf_unc = results_map.results[region + " mm"].sf_uncertainties[process];
            double SF_sf_unc = results_map.results[region + " SF"].sf_uncertainties[process];
            map<string, string> channel_sfs = {
                {"ee", "$" + toString(ee_sf, 3) + "\\pm" + toString(ee_sf_unc, 3) + "$"},
                {"mm", "$" + toString(mm_sf, 3) + "\\pm" + toString(mm_sf_unc, 3) + "$"},
                {"SF", "$" + toString(SF_sf, 3) + "\\pm" + toString(SF_sf_unc, 3) + "$"}
            };
            out_file << " & " << getChannelString(channel_sfs, results_map.plot_channels);
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

//------------
// MAKE PLOTS
//------------

tuple<THStack*, THStack*, THStack*, THStack*> createStacks(map<string, TH1D*> histograms, TString formatted_feature, Options options) {
    THStack *data_stack = new THStack("data_stack", "");
    THStack *raw_g_stack = new THStack("raw_g_stack", "");
    THStack *reweight_g_stack = new THStack("reweight_g_stack", "");
    THStack *zmc_stack = new THStack("zmc_stack", "");

    for (auto process_ptr = options.processes.rbegin(); process_ptr != options.processes.rend(); ++process_ptr) {
        auto process = *process_ptr;

        //--- set plotting options
        if (process == "photon") {
            if (options.plot_reweighted_photons) {
                histograms["photon_reweighted"]->SetLineColor(options.process_colors["photon_reweighted"]);
                histograms["photon_reweighted"]->SetLineWidth(0);
                histograms["photon_reweighted"]->SetFillColor(options.process_colors["photon_reweighted"]);
            }
            if (options.plot_unreweighted_photons) {
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
        }
        else if (process=="Zjets") {
            if (options.is_data && (options.plot_reweighted_photons || options.plot_unreweighted_photons)) {
                if (options.plot_zmc) {
                    histograms[process]->SetLineColor(kBlue);
                    histograms[process]->SetLineStyle(9);
                    histograms[process]->SetLineWidth(3);
                    histograms[process]->SetFillStyle(0);
                }
            }
            else {
                histograms[process]->SetFillColor(42);
                histograms[process]->SetLineStyle(1);
            }
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
            if (options.plot_unreweighted_photons)
                histograms["photon_raw"]->GetXaxis()->SetRange(0, histograms["photon_raw"]->GetNbinsX() + 1);
            if (options.plot_reweighted_photons)
                histograms["photon_reweighted"]->GetXaxis()
                    ->SetRange(0, histograms["photon_reweighted"]->GetNbinsX() + 1);
        }
        else {
            if (process == "Zjets" && options.is_data && !options.plot_zmc) continue;
            histograms[process]->GetXaxis()->SetRange(0, histograms[process]->GetNbinsX() + 1);
        }

        //--- add to relevant stack
        if (process == "data_bkg") {
            if (options.is_data) data_stack->Add(histograms[process]);
        }
        else if ((process == "Zjets") && (!options.is_data)) {
            data_stack->Add(histograms[process]);
        }
        else if (process == "photon") {
            raw_g_stack->Add(histograms["photon_raw"]);
            reweight_g_stack->Add(histograms["photon_reweighted"]);
        }
        else if (options.is_data) {
            if (options.plot_zmc)
                zmc_stack->Add(histograms[process]);
            if (process != "Zjets") {
                if (options.plot_unreweighted_photons) raw_g_stack->Add(histograms[process]);
                if (options.plot_reweighted_photons) reweight_g_stack->Add(histograms[process]);
            }
        }
    }

    return make_tuple(data_stack, raw_g_stack, reweight_g_stack, zmc_stack);
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
            if (process == "photon") {
                if (options.plot_reweighted_photons)
                    leg->AddEntry(histograms["photon_reweighted"],
                        options.process_latex["photon_reweighted"].c_str(), "f");
                if (options.plot_unreweighted_photons)
                    leg->AddEntry(histograms["photon_raw"], options.process_latex["photon_raw"].c_str(), "f");
            }
            else if (process == "Zjets") {
                if (options.plot_zmc)
                    leg->AddEntry(histograms["Zjets"], options.process_latex["Zjets"].c_str(), "f");
            }
            else if (process != "data_bkg") {
                leg->AddEntry(histograms[process], options.process_latex[process].c_str(), "f");
            }
        }
    }
    else {
        leg->AddEntry(histograms["Zjets"], options.process_latex["Zjets"].c_str(), "f");
        if (options.plot_unreweighted_photons)
            leg->AddEntry(histograms["photon_raw"], options.process_latex["photon_raw"].c_str(), "f");
        if (options.plot_reweighted_photons)
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

tuple<TH1D*, TH1D*, TH1D*> getRatioPlots(Options options, map<string, TH1D*> histograms) {
    TH1D *hratio, *hratio_unreweighted, *hratio_zmc, *hmctot, *hmctot_unreweighted, *hmctot_zmc;

    if (options.is_data) {
        hratio = (TH1D*) histograms["data_bkg"]->Clone("hratio");
        hratio_unreweighted = (TH1D*) histograms["data_bkg"]->Clone("hratio");
        hratio_zmc = (TH1D*) histograms["data_bkg"]->Clone("hratio");
        if (options.plot_reweighted_photons)
            hmctot = (TH1D*) histograms["photon_reweighted"]->Clone("hmctot");
        if (options.plot_unreweighted_photons)
            hmctot_unreweighted = (TH1D*) histograms["photon_raw"]->Clone("hmctot");
        if (options.plot_zmc)
            hmctot_zmc = (TH1D*) histograms["Zjets"]->Clone("hmctot");
        for (auto process : options.processes) {
            if ((process != "data_bkg") && (process != "Zjets") && (process != "photon")) {
                if (options.plot_reweighted_photons) hmctot->Add(histograms[process]);
                if (options.plot_unreweighted_photons) hmctot_unreweighted->Add(histograms[process]);
                if (options.plot_zmc) hmctot_zmc->Add(histograms[process]);
            }
        }
    }
    else {
        hratio = (TH1D*) histograms["Zjets"]->Clone("hratio");
        hratio_unreweighted = (TH1D*) histograms["Zjets"]->Clone("hratio");
        if (options.plot_reweighted_photons)
            hmctot = (TH1D*) histograms["photon_reweighted"]->Clone("hmctot");
        if (options.plot_unreweighted_photons)
            hmctot_unreweighted = (TH1D*) histograms["photon_raw"]->Clone("hmctot");
    }

    for (int ibin=1; ibin <= hmctot->GetXaxis()->GetNbins(); ibin++) {
        if (options.plot_reweighted_photons)
            hmctot->SetBinError(ibin, 0.0);
        if (options.plot_unreweighted_photons)
            hmctot_unreweighted->SetBinError(ibin, 0.0);
    }

    if (options.plot_unreweighted_photons) {
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
        hratio_unreweighted->SetTitle("");
    }

    if (options.plot_reweighted_photons) {
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
    }

    if (options.is_data && options.plot_zmc) {
        hratio_zmc->Divide(hmctot_zmc);
        hratio_zmc->SetMarkerStyle(20);
        auto zmc_color = histograms["Zjets"]->GetLineColor();
        hratio_zmc->SetMarkerColor(zmc_color);
        hratio_zmc->SetLineColor(zmc_color);
        hratio_zmc->GetXaxis()->SetTitle("");
        hratio_zmc->GetXaxis()->SetLabelSize(0.);
        hratio_zmc->GetYaxis()->SetNdivisions(5);
        hratio_zmc->GetYaxis()->SetTitle("");
        hratio_zmc->GetYaxis()->SetTitleSize(0.15);
        hratio_zmc->GetYaxis()->SetTitleOffset(0.3);
        hratio_zmc->GetYaxis()->SetLabelSize(0.15);
        hratio_zmc->SetMinimum(0.0);
        hratio_zmc->SetMaximum(2.0);
        hratio_zmc->GetYaxis()->SetRangeUser(0.0,2.0);
        hratio_zmc->SetTitle("");
    }

    //--- don't plot underflow and overflow
    if (options.plot_reweighted_photons) hratio->GetXaxis()->SetRange(1, hratio->GetNbinsX()-1);
    if (options.is_data && options.plot_zmc) hratio_zmc->GetXaxis()->SetRange(1, hratio_zmc->GetNbinsX()-1);
    if (options.plot_unreweighted_photons)
        hratio_unreweighted->GetXaxis()->SetRange(1, hratio_unreweighted->GetNbinsX()-1);

    return make_tuple(hratio, hratio_unreweighted, hratio_zmc);
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
                auto [data_stack, raw_g_stack, reweight_g_stack, zmc_stack] = createStacks(hist_map, formatted_feature, options);

                //--- draw plot
                TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);
                mainpad->Draw();
                mainpad->cd();

                double max_y = data_stack->GetMaximum() * 1.5;
                if (options.plot_reweighted_photons) max_y = max(reweight_g_stack->GetMaximum() * 1.5, max_y);
                if (options.plot_unreweighted_photons) max_y = max(raw_g_stack->GetMaximum() * 1.5, max_y);
                if (options.plot_zmc) max_y = max(zmc_stack->GetMaximum() * 1.5, max_y);
                double min_y = 0;
                //max_y /= 8;
                vector<string> log_features = {"met_Et", "METl", "METt", "met_Sign", "Ptll", "Ht30"};
                if (find(log_features.begin(), log_features.end(), feature) != log_features.end()) {
                    mainpad->SetLogy();
                    max_y = max_y / 1.5 * 50;
                    min_y = max_y / pow(10, 5);
                }

                bool applicable_blinded = (options.is_data)
                                          && (options.blinded && (region.find("SR") != std::string::npos));

                if (options.is_data) {
                    vector<THStack*> stacks;
                    if (options.plot_reweighted_photons) stacks.push_back(reweight_g_stack);
                    if (options.plot_unreweighted_photons)
                        stacks.push_back((THStack*)raw_g_stack->GetStack()->Last());
                    if (options.plot_zmc) {
                        if (options.plot_reweighted_photons || options.plot_unreweighted_photons)
                            stacks.push_back((THStack*)zmc_stack->GetStack()->Last());
                        else stacks.push_back(zmc_stack);
                    }
                    bool first_stack = true;
                    for (auto stack : stacks) {
                        if (first_stack) {
                            stack->Draw("hist");
                            first_stack = false;
                            stack->SetMaximum(max_y);
                            stack->SetMinimum(min_y);
                            stack->GetXaxis()->SetTitle(formatted_feature);
                            stack->GetYaxis()->SetTitle("entries / bin");
                        }
                        else stack->Draw("samehist");
                    }
                    if (!applicable_blinded)
                        data_stack->Draw("sameE1");
                }
                else {
                    data_stack->Draw("hist");
                    data_stack->SetMaximum(max_y);
                    data_stack->SetMinimum(min_y);
                    if (options.plot_unreweighted_photons)
                        raw_g_stack->GetStack()->Last()->Draw("samehist");
                    if (options.plot_reweighted_photons)
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

                auto [hratio, hratio_unreweighted, hratio_zmc] = getRatioPlots(options, hist_map);
                if (applicable_blinded) {
                    TH1D *empty_hist = new TH1D("", "", 1, 0, 1);
                    empty_hist->Draw();
                    tex->SetTextSize(0.3);
                    tex->DrawLatex(0.42,0.42,"BLINDED");
                }
                else {
                    if (options.plot_unreweighted_photons)
                        hratio_unreweighted->Draw("E1");
                    if (options.plot_zmc)
                        hratio_zmc->Draw("sameE1");
                    if (options.plot_reweighted_photons)
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

    options.period = "data15-16";
    options.data_period = DataPeriod(options.period);
    options.mc_period = getMCPeriod(options.period);

    options.is_data = true;
    options.reduction_folder = options.unit_test_folder + "ReducedNtuples/";
    options.reweighting_folder = options.unit_test_folder + "ReweightedNtuples/";
    options.plots_folder = "DiagnosticPlots/";

    options.make_diagnostic_plots = true;
    options.plot_regions = vector<string>{"VRZ", "VRZ_MET0_50"};
    options.plot_channels = vector<string>{"SF"};
    options.plot_features = vector<string>{"mll", "Ptll", "met_Et", "met_Sign", "mt2leplsp_0", "Ht30"};
    options.processes = {"data_bkg", "photon", "Zjets", "ttbar", "diboson"};

    options.blinded = true;
    options.print_photon_yield_only = false;
    options.do_vgamma_subtraction = false;

    options.plot_reweighted_photons = true;
    options.plot_unreweighted_photons = true;
    options.plot_zmc = true;

    options.scale_zmc = true;
    //options.scaling_method = "MET";
    options.scaling_method = "minDPhi2JetsMet";

    options.reweight_branch = "reweight_Ptll";

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
    printPhotonCRYieldTables(options, results_map, options.plots_folder + options.reweight_branch + "_" +
            options.data_period + "_" + type + "_CR_yields.txt", options.blinded);
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

    options.process_colors = {{"data_bkg", kBlack},
                             {"photon_raw", kYellow+2},
                             {"photon_reweighted", kGreen-1},
                             {"Zjets", kYellow-1},
                             {"Wjets", kSpring+1},
                             {"ttbar", kBlue+3},
                             {"diboson", kOrange+1},
                             {"lowMassDY", kGreen+4},
                             {"topOther", kAzure-9},
                             {"singleTop", kAzure+2},
                             {"triboson", kOrange-9},
                             {"higgs", kRed-10}};
    options.process_latex = {{"data_bkg", "Data"},
                             {"photon_raw", "Z+jets (from #gamma+jets, raw)"},
                             {"photon_reweighted", "Z+jets (from #gamma+jets, reweighted)"},
                             {"Zjets", "Z+jets (from MC)"},
                             {"Wjets", "W+jets"},
                             {"ttbar", "t#bar{t}"},
                             {"diboson", "Diboson"},
                             {"lowMassDY", "Low mass DY"},
                             {"topOther", "Top other"},
                             {"singleTop", "Single top"},
                             {"triboson", "Triboson"},
                             {"higgs", "Higgs"}};

    //--- either perform unit tests or run code
    if (options.unit_testing)
        performPlottingUnitTests(options);
    else
        run_quickDraw(options);
}

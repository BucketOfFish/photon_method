#include "../Common/Settings.C"
#include <unordered_map>
#include <boost/algorithm/string.hpp>

using namespace std;
using weightedDataFrame = ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void>;
struct Selection {
    string region;
    string channel;
    string feature;
    string process;
};
using histMap = unordered_map<string, unordered_map<string, unordered_map<string, ROOT::RDF::RResultPtr<TH1D>>>>; // [region][feature][process]
struct Result {
    unordered_map<string, unordered_map<string, TH1D*>> histograms; // [feature][process] histogram
    float photon_yield;
    float zmc_yield;
    float data_yield;
};
struct resultsMap {
    vector<string> regions;
    vector<string> features;
    vector<string> channels;
    string period;
    string data_or_mc;
    unordered_map<string, Result> results; // results by [region]
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
        cout << "Unrecognized channel! quitting   " << channel << endl;
        exit(0);
    }

    TCut plot_CR = plot_region + cuts::CR;
    plot_region += cuts::plot_region_met_portions[region];

    string region_name = region + " " + channel;

    return make_tuple(region_name, string(plot_region), string(plot_CR));
}

ROOT::RDF::TH1DModel getHistogramInfo(string plot_feature) {
    unordered_map<string, ROOT::RDF::TH1DModel> plot_settings;
    plot_settings["met_Et"] = ROOT::RDF::TH1DModel("", "E_{T}^{miss} [GeV]", 30, 0, 300);
    plot_settings["MET"] = ROOT::RDF::TH1DModel("", "E_{T}^{miss} [GeV]", 30, 0, 300);
    plot_settings["METl"] = ROOT::RDF::TH1DModel("", "E_{T,||}^{miss} [GeV]", 30, -150, 150);
    plot_settings["METt"] = ROOT::RDF::TH1DModel("", "E_{T,#perp}^{miss} [GeV]", 30, -150, 150);
    plot_settings["MET_loose"] = ROOT::RDF::TH1DModel("", "E_{T,loose}^{miss} [GeV]", 20, 0, 200);
    plot_settings["MET_tight"] = ROOT::RDF::TH1DModel("", "E_{T,tight}^{miss} [GeV]", 20, 0, 200);
    plot_settings["MET_tighter"] = ROOT::RDF::TH1DModel("", "E_{T,tighter}^{miss} [GeV]", 20, 0, 200);
    plot_settings["MET_tenacious"] = ROOT::RDF::TH1DModel("", "E_{T,tenacious}^{miss} [GeV]", 20, 0, 200);
    plot_settings["Ptll"] = ROOT::RDF::TH1DModel("", "p_{T} [GeV]", 20, 0, 100);
    plot_settings["Z_pt"] = ROOT::RDF::TH1DModel("", "p_{T} [GeV]", 20, 0, 100);
    plot_settings["nJet30"] = ROOT::RDF::TH1DModel("", "n_{jets}", 6, 2, 8);
    plot_settings["jet_n"] = ROOT::RDF::TH1DModel("", "n_{jets}", 6, 2, 8);
    plot_settings["bjet_n"] = ROOT::RDF::TH1DModel("", "n_{b-jets}", 4, 0, 4);
    plot_settings["HT"] = ROOT::RDF::TH1DModel("", "H_{T}", 20, 0, 1000);
    plot_settings["mll"] = ROOT::RDF::TH1DModel("", "m_{ll} [GeV]", 30, 0, 300);
    plot_settings["MT2"] = ROOT::RDF::TH1DModel("", "m_{T2} [GeV]", 20, 0, 200);
    plot_settings["MT2W"] = ROOT::RDF::TH1DModel("", "m_{T2}^{W} [GeV]", 20, 0, 200);
    plot_settings["lepPt[0]"] = ROOT::RDF::TH1DModel("", "lep_{p_{T},1} [GeV]", 20, 0, 300);
    plot_settings["lepPt[1]"] = ROOT::RDF::TH1DModel("", "lep_{p_{T},2} [GeV]", 20, 0, 200);
    plot_settings["lep_eta[0]"] = ROOT::RDF::TH1DModel("", "lep_{#eta,1} [GeV]", 30, -3, 3);
    plot_settings["lep_eta[1]"] = ROOT::RDF::TH1DModel("", "lep_{#eta,2} [GeV]", 30, -3, 3);
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

unordered_map<string, ROOT::RDataFrame*> getRDataFrames(string period, string photon_data_or_mc) {
    //--- load files
    string mc_period = getMCPeriod(period);
    string data_filename= ntuple_path + "bkg_data/" + period + "_bkg.root";
    string tt_filename = ntuple_path + "bkg_mc/" + mc_period + "_ttbar.root";
    string vv_filename = ntuple_path + "bkg_mc/" + mc_period + "_diboson.root";
    string zmc_filename = ntuple_path + "bkg_mc/" + mc_period + "_Zjets.root";
    string photon_ee_filename, photon_mm_filename;
    if (photon_data_or_mc == "MC") {
        photon_ee_filename = reweighting_path + "g_mc/" + mc_period + "_SinglePhoton222_ee.root";
        photon_mm_filename = reweighting_path + "g_mc/" + mc_period + "_SinglePhoton222_mm.root";
    }
    else {
        photon_ee_filename = reweighting_path + "g_data/" + period + "_photon_ee.root";
        photon_mm_filename = reweighting_path + "g_data/" + period + "_photon_mm.root";
    }

    cout << "data filename        " << data_filename << endl;
    cout << "ttbar filename       " << tt_filename << endl;
    cout << "diboson filename     " << vv_filename << endl;
    cout << "Z MC filename        " << zmc_filename << endl;
    cout << "photon filename (ee) " << photon_ee_filename << endl;
    cout << "photon filename (mm) " << photon_mm_filename << endl;
    cout << endl;

    //--- combine photon files
    TChain *tch_photon = new TChain("BaselineTree");
    tch_photon->Add(photon_ee_filename.c_str());
    tch_photon->Add(photon_mm_filename.c_str());

    //--- add files to RDataFrame
    unordered_map<string, ROOT::RDataFrame*> RDataFrames = {
        {"data", new ROOT::RDataFrame("BaselineTree", data_filename)},
        {"tt", new ROOT::RDataFrame("BaselineTree", tt_filename)},
        {"vv", new ROOT::RDataFrame("BaselineTree", vv_filename)},
        {"zmc", new ROOT::RDataFrame("BaselineTree", zmc_filename)},
        {"photon", new ROOT::RDataFrame(*tch_photon)},
    };

    cout << "data entries         " << *(RDataFrames["data"]->Count()) << endl;
    cout << "ttbar entries        " << *(RDataFrames["tt"]->Count()) << endl;
    cout << "diboson entries      " << *(RDataFrames["vv"]->Count()) << endl;
    cout << "Z MC entries         " << *(RDataFrames["zmc"]->Count()) << endl;
    cout << "photon entries       " << *(RDataFrames["photon"]->Count()) << endl;
    cout << endl;

    return RDataFrames;
}

unordered_map<string, unique_ptr<weightedDataFrame>> weightRDataFrames(unordered_map<string, ROOT::RDataFrame*> dataframes) {
    unordered_map<string, string> plot_weights;
    plot_weights["data"] = "1";
    plot_weights["tt"] = cuts::bkg_weight;
    plot_weights["vv"] = cuts::bkg_weight;
    plot_weights["zmc"] = cuts::bkg_weight;
    plot_weights["photon_raw"] = cuts::photon_weight;
    plot_weights["photon_reweighted"] = cuts::photon_weight_rw;

    unordered_map<string, unique_ptr<weightedDataFrame>> weighted_dataframes;
    for (auto const& [process, dataframe] : dataframes) {
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

//-----------------
// FILL HISTOGRAMS
//-----------------

tuple<histMap, histMap> setUpHistograms(unordered_map<string, unique_ptr<weightedDataFrame>>* WRDataFramesPtr, vector<string> plot_features, vector<string> regions, vector<string> channels) {
    cout << "setting up histograms" << endl;
    cout << endl;

    histMap plot_region_histograms;
    histMap control_region_histograms;

    for (string region : regions) {
        for (string channel : channels) {
            auto [region_name, plot_region, plot_CR] = getPlotRegionInfo(channel, region);
            cout << "\tregion name          " << region_name << endl;
            cout << "\tplot region          " << plot_region << endl;
            cout << "\tnormalization reg.   " << plot_CR << endl;
            cout << endl;

            for (auto plot_feature : plot_features) {
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

resultsMap fillHistograms(tuple<histMap, histMap> region_hists, vector<string> plot_features, vector<string> regions, vector<string> channels, string photon_data_or_mc) {
    cout << "making histograms" << endl;
    cout << endl;

    auto [plot_region_histograms, control_region_histograms] = region_hists;
    resultsMap results_map;
    results_map.regions = regions;
    results_map.features = plot_features;
    results_map.channels = channels;
    results_map.data_or_mc = photon_data_or_mc;
    for (string region : regions) {
        for (string channel : channels) {
            auto [region_name, plot_region, plot_CR] = getPlotRegionInfo(channel, region);
            auto prh = plot_region_histograms[region_name];
            auto crh = control_region_histograms[region_name];
            string first_feature = plot_features[0];
            auto prh0 = prh[first_feature];
            auto crh0 = crh[first_feature];

            int nbins = getHistogramInfo(first_feature).fNbinsX;
            auto fullInt = [nbins](TH1D* hist) {return hist->Integral(0, nbins+1);};

            cout << "\tregion name          " << region_name << endl;
            cout << endl;
            cout << "\tdata integral        " << prh0["data"]->Integral(0, nbins+1) << endl;
            cout << "\tttbar integral       " << prh0["tt"]->Integral(0, nbins+1) << endl;
            cout << "\tdiboson integral     " << prh0["vv"]->Integral(0, nbins+1) << endl;
            cout << "\tZ MC integral        " << prh0["zmc"]->Integral(0, nbins+1) << endl;
            cout << "\tg raw integral       " << prh0["photon_raw"]->Integral(0, nbins+1) << endl;
            cout << "\tg reweighted int.    " << prh0["photon_reweighted"]->Integral(0, nbins+1) << endl;
            cout << endl;

            float zdata_integral = crh0["data"]->Integral(0, nbins+1) - crh0["tt"]->Integral(0, nbins+1) - crh0["vv"]->Integral(0, nbins+1);
            if (photon_data_or_mc == "MC")
                zdata_integral = crh0["zmc"]->Integral(0, nbins+1);
            float SF = zdata_integral / crh0["photon_raw"]->Integral(0, nbins+1);
            if (crh0["photon_raw"]->Integral(0, nbins+1) == 0) SF = 0;
            float SFrw = zdata_integral / crh0["photon_reweighted"]->Integral(0, nbins+1);
            if (crh0["photon_reweighted"]->Integral(0, nbins+1) == 0) SFrw = 0;

            cout << "\tScaling raw photon data by " << SF << endl;
            cout << "\tScaling reweighted photon data by " << SFrw << endl;

            for (auto plot_feature : plot_features) {
                prh[plot_feature]["photon_raw"]->Scale(SF);
                prh[plot_feature]["photon_reweighted"]->Scale(SFrw);
            }

            float photon_yield = prh0["photon_reweighted"]->Integral(0, nbins+1);
            float zmc_yield = prh0["zmc"]->Integral(0, nbins+1);
            float zdata_yield = (photon_data_or_mc == "MC") ? -999 : zdata_integral;
            cout << "\tPhoton yield of " << photon_yield << endl;
            cout << endl;

            unordered_map<string, unordered_map<string, TH1D*>> region_hists;
            for (auto plot_feature : plot_features) {
                for (auto [process, histogram] : prh[plot_feature]) {
                    TH1D* new_histogram = new TH1D;
                    histogram->Copy(*new_histogram);
                    region_hists[plot_feature][process] = new_histogram;
                }
            }
            results_map.results[region_name] = Result{region_hists, photon_yield, zmc_yield, zdata_yield};
        }
    }

    return results_map;
}

//-------------
// MAKE TABLES
//-------------

void printPhotonYieldTables(resultsMap results_map, string save_name, bool blinded) {
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
    out_file << "\\caption{Photon Method Yields}" << endl;
    out_file << "\\begin{center}" << endl;
    out_file << "\\begin{tabular}{c|c|c|c}" << endl;
    out_file << "region & photon ee / mm / SF & Z MC ee / mm / SF & data ee / mm / SF \\\\" << endl;
    out_file << "\\hline" << endl;

    for (auto region : results_map.regions) {
        float photon_ee = results_map.results[region + " ee"].photon_yield;
        float photon_mm = results_map.results[region + " mm"].photon_yield;
        float photon_SF = results_map.results[region + " SF"].photon_yield;
        float zmc_ee = results_map.results[region + " ee"].zmc_yield;
        float zmc_mm = results_map.results[region + " mm"].zmc_yield;
        float zmc_SF = results_map.results[region + " SF"].zmc_yield;
        float data_ee = results_map.results[region + " ee"].data_yield;
        float data_mm = results_map.results[region + " mm"].data_yield;
        float data_SF = results_map.results[region + " SF"].data_yield;
        if (blinded && (region.find("SR") != std::string::npos))
            out_file << region << " & " << photon_ee << " / " << photon_mm << " / " << photon_SF << " & " << zmc_ee << " / " << zmc_mm << " / " << zmc_SF << " & - / - / - \\\\" << endl;
        else
            out_file << region << " & " << photon_ee << " / " << photon_mm << " / " << photon_SF << " & " << zmc_ee << " / " << zmc_mm << " / " << zmc_SF << " & " << data_ee << " / " << data_mm << " / " << data_SF << " \\\\" << endl;
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

tuple<THStack*, THStack*, THStack*> createStacks(unordered_map<string, TH1D*> histograms, string photon_data_or_mc, TString formatted_feature) {
    //--- set plotting options
    histograms["tt"]->SetLineColor(1); histograms["tt"]->SetFillColor(kRed-2);
    histograms["vv"]->SetLineColor(1); histograms["vv"]->SetFillColor(kGreen-2);
    histograms["photon_raw"]->SetLineColor(4); histograms["photon_raw"]->SetLineWidth(1); histograms["photon_raw"]->SetLineStyle(2);
    if (photon_data_or_mc == "Data") {
        histograms["data"]->SetLineColor(1); histograms["data"]->SetLineWidth(2); histograms["data"]->SetMarkerStyle(20);
        histograms["zmc"]->SetLineColor(2); histograms["zmc"]->SetLineWidth(1); histograms["zmc"]->SetLineStyle(7);
        histograms["photon_reweighted"]->SetLineColor(1); histograms["photon_reweighted"]->SetFillColor(kOrange-2);
    }
    else {
        histograms["zmc"]->SetLineColor(1); histograms["zmc"]->SetFillColor(42); histograms["zmc"]->SetLineStyle(1);
        histograms["photon_reweighted"]->SetLineWidth(1); histograms["photon_reweighted"]->SetLineColor(kRed); histograms["photon_reweighted"]->SetFillStyle(0);
    }

    //--- turn on overflow bin
    vector<string> processes{"data", "tt", "vv", "zmc", "photon_raw", "photon_reweighted"};
    for (auto process : processes)
        histograms[process]->GetXaxis()->SetRange(0, histograms[process]->GetNbinsX() + 1);

    //--- make stacks
    THStack *data_stack = new THStack("data_stack", "");
    THStack *raw_g_stack = new THStack("raw_g_stack", "");
    THStack *reweight_g_stack = new THStack("reweight_g_stack", "");

    if (photon_data_or_mc == "Data")
        data_stack->Add(histograms["data"]);
    else
        data_stack->Add(histograms["zmc"]);

    if (photon_data_or_mc == "Data") {
        raw_g_stack->Add(histograms["tt"]);
        raw_g_stack->Add(histograms["vv"]);
    }
    raw_g_stack->Add(histograms["photon_raw"]);

    if (photon_data_or_mc == "Data") {
        reweight_g_stack->Add(histograms["tt"]);
        reweight_g_stack->Add(histograms["vv"]);
    }
    reweight_g_stack->Add(histograms["photon_reweighted"]);

    return make_tuple(data_stack, raw_g_stack, reweight_g_stack);
}

TString getPlotSaveName(string period, string channel, string plot_feature, string photon_data_or_mc, string region) {
    TString plot_name;
    if (photon_data_or_mc == "Data")
        plot_name = Form("%s/%s_%s_%s_%s", plots_path.c_str(), period.c_str(), channel.c_str(), plot_feature.c_str(), region.c_str());
    else
        plot_name = Form("%s/%s_%s_%s_%s", plots_path.c_str(), getMCPeriod(period).c_str(), channel.c_str(), plot_feature.c_str(), region.c_str());
    plot_name += ".eps";
    return plot_name;
}

TLegend* getLegend(string photon_data_or_mc, unordered_map<string, TH1D*> histograms) {
    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
    if (photon_data_or_mc == "Data") {
        leg->AddEntry(histograms["data"],"data","lp");
        //leg->AddEntry(histograms["zmc"], "Z+jets (from MC)", "f");
        //leg->AddEntry(histograms["photon_raw"], "Z+jets (from #gamma+jets, raw)", "f");
        leg->AddEntry(histograms["photon_reweighted"], "Z+jets (from #gamma+jets, reweighted)", "f");
        leg->AddEntry(histograms["vv"], "VV", "f");
        leg->AddEntry(histograms["tt"], "t#bar{t}+tW", "f");
    }
    else {
        leg->AddEntry(histograms["zmc"], "Z+jets (from MC)", "f");
        leg->AddEntry(histograms["photon_raw"], "Z+jets (from #gamma+jets, raw)", "f");
        leg->AddEntry(histograms["photon_reweighted"], "Z+jets (from #gamma+jets, reweighted)", "f");
    }

    leg->SetBorderSize(0);
    leg->SetFillColor(0);

    return leg;
}

string getPlotTex(string period, string photon_data_or_mc) {
    string tex_string;
    if (photon_data_or_mc == "Data") {
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

tuple<TH1D*, TH1D*> getRatioPlots(unordered_map<string, TH1D*> histograms, string photon_data_or_mc) {
    TH1D *hratio, *hratio_unreweighted, *hmctot, *hmctot_unreweighted;

    if (photon_data_or_mc == "MC") {
        hratio = (TH1D*) histograms["zmc"]->Clone("hratio");
        hratio_unreweighted = (TH1D*) histograms["zmc"]->Clone("hratio");
        hmctot = (TH1D*) histograms["photon_reweighted"]->Clone("hmctot");
        hmctot_unreweighted = (TH1D*) histograms["photon_raw"]->Clone("hmctot");
    }
    else {
        hratio = (TH1D*) histograms["data"]->Clone("hratio");
        hratio_unreweighted = (TH1D*) histograms["data"]->Clone("hratio");
        hmctot = (TH1D*) histograms["photon_reweighted"]->Clone("hmctot");
        hmctot->Add(histograms["tt"]);
        hmctot->Add(histograms["vv"]);
        hmctot_unreweighted = (TH1D*) histograms["photon_raw"]->Clone("hmctot");
        hmctot_unreweighted->Add(histograms["tt"]);
        hmctot_unreweighted->Add(histograms["vv"]);
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
    if (photon_data_or_mc == "Data")
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

void makePlot(resultsMap results_map, string period, bool blinded) {
    for (auto region : results_map.regions) {
        vector<string> channels = {"ee", "mm", "SF"};
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
                TString plot_title = formatted_feature + " in " + region;
                //TH1D *h_name = new TH1D("h_name", plot_title, nbins, xmin, xmax);
                TH1D *h_name = new TH1D("h_name", plot_title, 1, 0, 1);
                h_name->Draw();

                //--- create comparison stacks
                auto hist_map = results_map.results[region_name].histograms[feature]; // map of histograms by [process]
                auto [data_stack, raw_g_stack, reweight_g_stack] = createStacks(hist_map, results_map.data_or_mc, formatted_feature);

                //--- draw plot
                TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);
                mainpad->Draw();
                mainpad->cd();
                mainpad->SetLogy();

                bool applicable_blinded = (blinded && (region.find("SR") != std::string::npos));

                if (results_map.data_or_mc == "Data") {
                    reweight_g_stack->Draw("hist");
                    if (!applicable_blinded)
                        data_stack->Draw("sameE1");
                    reweight_g_stack->GetXaxis()->SetTitle(formatted_feature);
                    reweight_g_stack->GetYaxis()->SetTitle("entries / bin");
                }
                else {
                    if (!applicable_blinded)
                        data_stack->Draw("hist");
                    raw_g_stack->Draw("samehist");
                    reweight_g_stack->Draw("samehist");
                    data_stack->GetXaxis()->SetTitle(formatted_feature);
                    data_stack->GetYaxis()->SetTitle("entries / bin");
                }

                //--- draw legend and labels
                TLegend *leg = getLegend(results_map.data_or_mc, hist_map);
                leg->Draw();

                //--- write info
                TLatex *tex = new TLatex();
                tex->SetNDC();
                tex->SetTextSize(0.03);
                tex->DrawLatex(0.6,0.65,"ATLAS Internal");
                tex->DrawLatex(0.6,0.61,getPlotTex(period, results_map.data_or_mc).c_str());
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

                auto [hratio, hratio_unreweighted] = getRatioPlots(hist_map, results_map.data_or_mc);
                if (applicable_blinded) {
                    TH1D *empty_hist = new TH1D("", "", 1, 0, 1);
                    empty_hist->Draw();
                    tex->SetTextSize(0.3);
                    tex->DrawLatex(0.42,0.42,"BLINDED");
                }
                else {
                    if (results_map.data_or_mc == "MC")
                        hratio_unreweighted->Draw("E1");
                    hratio->Draw("sameE1");
                }

                //--- save plot
                TString plot_name = getPlotSaveName(period, channel, feature, results_map.data_or_mc, region);
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

void testTablePrintout(resultsMap results_map) {
    printPhotonYieldTables(results_map, "Output/test_table_blindeded.txt", true);
    printPhotonYieldTables(results_map, "Output/test_table_unblindeded.txt", false);
}

void testMakePlot(resultsMap results_map) {
    results_map.data_or_mc = "Data";
    makePlot(results_map, "data18", true);
    results_map.data_or_mc = "MC";
    makePlot(results_map, "data18", false);
}

void unit_tests() {
    resultsMap results_map;
    vector<string> regions{"SR_test1", "SR_test2", "SR_test3", "VR_test1"};
    vector<string> features{"METl", "METt", "met_Et"};
    results_map.regions = regions;
    results_map.features = features;

    unordered_map<string, unordered_map<string, TH1D*>> test_hists; // [feature][process] histogram
    vector<string> processes{"data", "tt", "vv", "zmc", "photon_raw", "photon_reweighted"};
    unordered_map<string, int> n_entries;
    n_entries["data"] = 10000;
    n_entries["tt"] = 6000;
    n_entries["vv"] = 3000;
    n_entries["zmc"] = 1000;
    n_entries["photon_raw"] = 1000;
    n_entries["photon_reweighted"] = 1000;
    for (auto feature : features) {
        for (auto process : processes) {
            test_hists[feature][process] = new TH1D("", "", 100, -3, 3);
            test_hists[feature][process]->FillRandom("gaus", n_entries[process]);
        }
    }

    results_map.results["SR_test1 ee"] = Result{test_hists, 1.3, 1.5, 1.4};
    results_map.results["SR_test1 mm"] = Result{test_hists, 1.2, 1.6, 1.3};
    results_map.results["SR_test1 SF"] = Result{test_hists, 2.5, 3.1, 2.7};

    results_map.results["SR_test2 ee"] = Result{test_hists, 3.6, 3.5, 3.7};
    results_map.results["SR_test2 mm"] = Result{test_hists, 3.1, 3.2, 3.4};
    results_map.results["SR_test2 SF"] = Result{test_hists, 6.7, 6.7, 7.1};

    results_map.results["SR_test3 ee"] = Result{test_hists, 12.9, 13.3, 11.2};
    results_map.results["SR_test3 mm"] = Result{test_hists, 11.6, 12.8, 12.0};
    results_map.results["SR_test3 SF"] = Result{test_hists, 24.5, 26.1, 23.2};

    results_map.results["VR_test1 ee"] = Result{test_hists, 112.9, 113.3, 111.2};
    results_map.results["VR_test1 mm"] = Result{test_hists, 111.6, 112.8, 112.0};
    results_map.results["VR_test1 SF"] = Result{test_hists, 124.5, 126.1, 123.2};

    testTablePrintout(results_map);
    testMakePlot(results_map);

    n_entries["data"] = 0;
    n_entries["tt"] = 0;
    n_entries["vv"] = 0;
    n_entries["zmc"] = 0;
    n_entries["photon_raw"] = 0;
    n_entries["photon_reweighted"] = 0;
    for (auto feature : features) {
        for (auto process : processes) {
            test_hists[feature][process] = new TH1D("", "", 100, -3, 3);
            test_hists[feature][process]->FillRandom("gaus", n_entries[process]);
        }
    }

    results_map.results["SR_test1 ee"] = Result{test_hists, 0, 0, 0};
    results_map.results["SR_test1 mm"] = Result{test_hists, 0, 0, 0};
    results_map.results["SR_test1 SF"] = Result{test_hists, 0, 0, 0};

    results_map.results["SR_test2 ee"] = Result{test_hists, 0, 0, 0};
    results_map.results["SR_test2 mm"] = Result{test_hists, 0, 0, 0};
    results_map.results["SR_test2 SF"] = Result{test_hists, 0, 0, 0};

    results_map.results["SR_test3 ee"] = Result{test_hists, 0, 0, 0};
    results_map.results["SR_test3 mm"] = Result{test_hists, 0, 0, 0};
    results_map.results["SR_test3 SF"] = Result{test_hists, 0, 0, 0};

    results_map.results["VR_test1 ee"] = Result{test_hists, 0, 0, 0};
    results_map.results["VR_test1 mm"] = Result{test_hists, 0, 0, 0};
    results_map.results["VR_test1 SF"] = Result{test_hists, 0, 0, 0};

    //testTablePrintout(results_map);
    //testMakePlot(results_map);
}

//---------------
// MAIN FUNCTION
//---------------

void run_quickDraw(string period, string photon_data_or_mc, string plot_feature_list, string region_list, bool blinded, bool print_photon_yield_only) {
    cout << "period               " << period << endl;
    cout << "photon data          " << photon_data_or_mc << endl;
    cout << endl;

    //--- parse arguments
    vector<string> plot_features = splitStringBy(plot_feature_list, " ");
    vector<string> regions = splitStringBy(region_list, " ");
    vector<string> channels{"ee", "mm", "SF", "DF"};

    //--- get input data in the form of RDataFrames
    unordered_map<string, ROOT::RDataFrame*> RDataFrames = getRDataFrames(period, photon_data_or_mc);
    unordered_map<string, unique_ptr<weightedDataFrame>> WRDataFrames = weightRDataFrames(RDataFrames);

    //--- set up all histograms and fill simulaneously
    auto region_hists = setUpHistograms(&WRDataFrames, plot_features, regions, channels);
    resultsMap results_map = fillHistograms(region_hists, plot_features, regions, channels, photon_data_or_mc);

    //--- print photon yield tables
    TString plot_name = getPlotSaveName(period, "yields", "allFeatures", results_map.data_or_mc, "allRegions");
    printPhotonYieldTables(results_map, "Output/" + period + "_" + results_map.data_or_mc + "_yields.txt", blinded);

    //--- draw and save plot
    if (!print_photon_yield_only)
        makePlot(results_map, period, blinded);
}

void quickDraw(string period, string photon_data_or_mc, string plot_feature_list, string region_list, bool blinded, bool print_photon_yield_only) {
    //--- set global options
    gStyle->SetOptStat(0);
    ROOT::EnableImplicitMT();

    //--- either perform unit tests or run code
    //unit_tests();
    run_quickDraw(period, photon_data_or_mc, plot_feature_list, region_list, blinded, print_photon_yield_only);
}

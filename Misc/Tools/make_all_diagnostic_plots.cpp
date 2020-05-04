#include "../../Main/Settings.cpp"
#include "../../Main/MakePlots.cpp"

map<string, ROOT::RDataFrame*> getDataFrames(map<string, string> folders) {
    map<string, ROOT::RDataFrame*> dataframes;
    map<string, TChain*> tchains;

    tchains["data_bkg"] = new TChain("BaselineTree");
    tchains["data_bkg"]->Add((folders["reduced"] + "data15-16_data_bkg.root").c_str());
    tchains["data_bkg"]->Add((folders["reduced"] + "data17_data_bkg.root").c_str());
    tchains["data_bkg"]->Add((folders["reduced"] + "data18_data_bkg.root").c_str());
    dataframes["data_bkg"] = new ROOT::RDataFrame(*tchains["data_bkg"]);

    vector<string> processes = {"ttbar", "Zjets"};
    for (auto process : processes) {
        tchains[process] = new TChain("BaselineTree");
        tchains[process]->Add((folders["reduced"] + "mc16a_" + process + ".root").c_str());
        tchains[process]->Add((folders["reduced"] + "mc16cd" + process + ".root").c_str());
        tchains[process]->Add((folders["reduced"] + "mc16e" + process + ".root").c_str());
        dataframes[process] = new ROOT::RDataFrame(*tchains[process]);
    }

    tchains["photon_smeared"] = new TChain("BaselineTree");
    tchains["photon_smeared"]->Add((folders["smeared"] + "data15-16_data_photon_ee.root").c_str());
    tchains["photon_smeared"]->Add((folders["smeared"] + "data15-16_data_photon_mm.root").c_str());
    tchains["photon_smeared"]->Add((folders["smeared"] + "data17_data_photon_ee.root").c_str());
    tchains["photon_smeared"]->Add((folders["smeared"] + "data18_data_photon_mm.root").c_str());
    tchains["photon_smeared"]->Add((folders["smeared"] + "data17_data_photon_ee.root").c_str());
    tchains["photon_smeared"]->Add((folders["smeared"] + "data18_data_photon_mm.root").c_str());
    dataframes["photon_smeared"] = new ROOT::RDataFrame(*tchains["photon_smeared"]);

    tchains["photon_reweighted"] = new TChain("BaselineTree");
    tchains["photon_reweighted"]->Add((folders["reweighted"] + "data15-16_data_photon_ee.root").c_str());
    tchains["photon_reweighted"]->Add((folders["reweighted"] + "data15-16_data_photon_mm.root").c_str());
    tchains["photon_reweighted"]->Add((folders["reweighted"] + "data17_data_photon_ee.root").c_str());
    tchains["photon_reweighted"]->Add((folders["reweighted"] + "data18_data_photon_mm.root").c_str());
    tchains["photon_reweighted"]->Add((folders["reweighted"] + "data17_data_photon_ee.root").c_str());
    tchains["photon_reweighted"]->Add((folders["reweighted"] + "data18_data_photon_mm.root").c_str());
    dataframes["photon_reweighted"] = new ROOT::RDataFrame(*tchains["photon_reweighted"]);

    return dataframes;
}

TH1F* makePlot(map<string, ROOT::RDataFrame*> dataframes,
    tuple<string, string, string, string> plot_info) {
    /// For a given plot region and feature, return a map of histograms by [process][channel].
    auto [region, feature, photon_type, data_or_mc] = plot_info;
    ROOT::RDF::TH1DModel hist_model = getHistogramInfo(feature);

    //--- draw histograms for each process
    vector<string> processes = {"ttbar", "photon"};
    dataframes["photon"] = dataframes["photon_" + photon_type];
    if (data_or_mc == "data") processes.push_back("data_bkg");
    else if (data_or_mc == "MC") processes.push_back("Zjets");
    else failTest("Choose either data or MC");

    map<string, float> scale_factors = {{"SRC", 0.000}, {"SRLow2", 1.852}, {"SRMed2", 1.517}, {"SRHigh2", 1.649},
        {"SRLow23", 1.854}, {"SRMed23", 1.571}, {"SRHigh23", 1.715}, {"SRLow4", 1.847}, {"SRMed4", 1.472},
        {"SRHigh4", 1.618}, {"SRLowZ4", 1.862}, {"SRMedZ4", 1.457}, {"SRHighZ4", 1.595}, {"SRLowZ6", 1.721},
        {"SRMedZ6", 1.288}, {"SRHighZ6", 1.340}, {"VRC", 0.000}, {"VRLow2", 1.852}, {"VRMed2", 1.517},
        {"VRHigh2", 1.649}, {"VRLow23", 1.854}, {"VRMed23", 1.571}, {"VRHigh23", 1.715}, {"VRLow4", 1.847},
        {"VRMed4", 1.472}, {"VRHigh4", 1.618}, {"VRLowZ4", 1.862}, {"VRMedZ4", 1.457}, {"VRHighZ4", 1.595},
        {"VRLowZ6", 1.721}, {"VRMedZ6", 1.288}, {"VRHighZ6", 1.340}, {"VRDPhiLow2", 1.634}, {"VRDPhiMed2", 1.402},
        {"VRDPhiHigh2", 1.659}, {"VRDPhiLow6", 1.517}, {"VRDPhiMed6", 1.063}, {"VRDPhiHigh6", 1.362},
        {"reweight", 1.29}};

    map<string, map<string, ROOT::RDF::RResultPtr<TH1D>>> hists;
    vector<string> channels = {"ee", "mm", "SF"};
    for (string process : processes) {
        string weight = "totalWeight";
        if (process == "data") weight = "1";
        if (process == "photon" && photon_type == "reweighted") {
            weight = "totalWeight * reweight_Ptll__Ht30";
            weight += (" * " + to_string(scale_factors[region]));
        }
        for (string channel : channels) {
            auto [region_name, plot_region, plot_CR] = getPlotRegionInfo(channel, region);
            hists[process][channel] = dataframes[process]->Filter(plot_region)
                                .Define("rw_weight", weight)
                                .Histo1D(hist_model, feature, "rw_weight");
        }
    }

    //--- combine histograms into plot
    for (string region : regions) {
        for (string feature : features) {
            for (string channel : {"ee", "mm"}) {
                auto [region_name, plot_region, plot_CR] = getPlotRegionInfo(channel, region);
                TCanvas *can = new TCanvas("can","can",600,600);
                hists[region_name][feature]["photon"]->Draw("hist");
                hists[region_name][feature]["data"]->Draw("same e1");
                can->Write((region_name+"_"+feature+".eps").c_str());
                //can->Print(("Plots/"+region_name+"_"+feature+".eps").c_str());
                delete can;
            }
        }
    }

    return plot;
}

void make_all_diagnostic_plots() {
    map<string, string> folders;
    folders["reduced"] = "/public/data/Photon/Samples/ReducedNtuples/";
    folders["smeared"] = "/public/data/Photon/Samples/SmearedNtuples/";
    folders["reweighted"] = "/public/data/Photon/Samples/ReweightedNtuples/";
    map<string, ROOT::RDataFrame*> dataframes = getDataFrames(folders);

    vector<tuple<string, string, string, string>> requested_plots;
    requested_plots.push_back(make_tuple("baseline", "met_Et", "smeared", "data"));
    requested_plots.push_back(make_tuple("baseline", "Ptll", "smeared", "data"));
    requested_plots.push_back(make_tuple("baseline", "Ht30", "smeared", "data"));
    requested_plots.push_back(make_tuple("reweight", "met_Et", "reweighted", "data"));
    requested_plots.push_back(make_tuple("reweight", "Ptll", "reweighted", "data"));
    requested_plots.push_back(make_tuple("reweight", "Ht30", "reweighted", "data"));

    for (auto plot_info : requested_plots) {
        auto [region, feature, photon_type, data_or_mc] = plot_info;
        makePlot(dataframes, plot_info);
    }
}

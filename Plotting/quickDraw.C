#include "../Common/Settings.C"
#include <unordered_map>
#include <boost/algorithm/string.hpp>

using namespace std;
using weightedDataFrame = ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void>;
using histDictionary = unordered_map<string, unordered_map<string, unordered_map<string, ROOT::RDF::RResultPtr<TH1D>>>>; // [region][feature][process]
using histResults = unordered_map<string, unordered_map<string, tuple<unordered_map<string, TH1D*>, float>>>; // [region][feature], with a dictionary of hists by process and a photon yield value

//------------------
// HELPER FUNCTIONS
//------------------

vector<string> splitStringBySpaces(string input) {
    vector<string> output;
    boost::algorithm::split(output, input, boost::is_any_of(" "));
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

tuple<histDictionary, histDictionary> setUpHistograms(unordered_map<string, unique_ptr<weightedDataFrame>>* WRDataFramesPtr, vector<string> plot_features, vector<string> regions, vector<string> channels) {
    cout << "setting up histograms" << endl;
    cout << endl;

    histDictionary plot_region_histograms;
    histDictionary control_region_histograms;

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

histResults fillHistograms(histDictionary plot_region_histograms, histDictionary control_region_histograms, vector<string> plot_features, vector<string> regions, vector<string> channels, string photon_data_or_mc) {
    cout << "making histograms" << endl;
    cout << endl;

    histResults hist_results;
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
            float SFrw = zdata_integral / crh0["photon_reweighted"]->Integral(0, nbins+1);

            cout << "\tScaling raw photon data by " << SF << endl;
            cout << "\tScaling reweighted photon data by " << SFrw << endl;

            for (auto plot_feature : plot_features) {
                prh[plot_feature]["photon_raw"]->Scale(SF);
                prh[plot_feature]["photon_reweighted"]->Scale(SFrw);
            }

            float photon_yield = prh0["photon_reweighted"]->Integral(0, nbins+1);
            cout << "\tPhoton yield of " << photon_yield << endl;
            cout << endl;

            for (auto plot_feature : plot_features) {
                unordered_map<string, TH1D*> histograms;
                for (auto [process, histogram] : prh[plot_feature]) {
                    histograms[process] = new TH1D;
                    histogram->Copy(*histograms[process]);
                }
                hist_results[region_name][plot_feature] = make_tuple(histograms, photon_yield);
            }
        }
    }

    return hist_results;
}

//-------------
// MAKE TABLES
//-------------

void printPhotonYieldTables(histResults filled_histograms) {
    // [region][feature], with a dictionary of hists by process and a photon yield value
    cout << "photon yields" << endl;
    cout << endl;
    for (auto [region, region_hists] : filled_histograms) {
        string first_feature = region_hists.begin()->first;
        cout << region << ": " << get<1>(region_hists[first_feature]) << endl;
    }
}

//------------
// MAKE PLOTS
//------------

THStack* createStacks(unordered_map<string, TH1D*> histograms, string photon_data_or_mc, bool DF) {
    //--- create MC stack
    THStack *mcstack = new THStack("mcstack", "");

    histograms["tt"]->SetLineColor(1); histograms["tt"]->SetFillColor(kRed-2);
    histograms["vv"]->SetLineColor(1); histograms["vv"]->SetFillColor(kGreen-2);
    histograms["photon_reweighted"]->SetLineColor(1); histograms["photon_reweighted"]->SetFillColor(kOrange-2);
    mcstack->Add(histograms["tt"]);
    mcstack->Add(histograms["vv"]);
    if(!DF) mcstack->Add(histograms["photon_reweighted"]);

    //--- create comparison "stacks"
    if (photon_data_or_mc == "Data") {
        histograms["photon_raw"]->Add(histograms["tt"]); histograms["photon_raw"]->Add(histograms["vv"]);
        histograms["zmc"]->Add(histograms["tt"]); histograms["zmc"]->Add(histograms["vv"]);
    }

    //--- set plotting options
    histograms["photon_raw"]->SetLineColor(4); histograms["photon_raw"]->SetLineWidth(1); histograms["photon_raw"]->SetLineStyle(2);
    histograms["zmc"]->SetLineColor(2); histograms["zmc"]->SetLineWidth(1); histograms["zmc"]->SetLineStyle(7);

    //--- turn on overflow bin
    histograms["zmc"]->GetXaxis()->SetRange(0, histograms["zmc"]->GetNbinsX() + 1);
    histograms["photon_raw"]->GetXaxis()->SetRange(0, histograms["photon_raw"]->GetNbinsX() + 1);
    histograms["photon_reweighted"]->GetXaxis()->SetRange(0, histograms["photon_reweighted"]->GetNbinsX() + 1);

    return mcstack;
}

TString getPlotName(string period, string channel, string plot_feature, string photon_data_or_mc, string region) {
    TString plot_name;
    if (photon_data_or_mc == "Data")
        plot_name = Form("%s/%s_%s_%s_%s", plots_path.c_str(), period.c_str(), channel.c_str(), plot_feature.c_str(), region.c_str());
    else
        plot_name = Form("%s/%s_%s_%s_%s", plots_path.c_str(), getMCPeriod(period).c_str(), channel.c_str(), plot_feature.c_str(), region.c_str());
    plot_name += ".eps";
    return plot_name;
}

void makePlot(unordered_map<string, TH1D*> histograms, THStack *mcstack, TString plot_name, TString formatted_feature, string period, string channel, string photon_data_or_mc, string region) {
    //--- draw title
    TCanvas *can = new TCanvas("can","can",600,600);
    can->cd();
    TPad* namepad = new TPad("namepad","namepad",0.0,0.0,1.0,1.0);
    namepad->Draw();
    namepad->cd();
    TString plot_title = formatted_feature + " in " + region;
    //TH1D *h_name = new TH1D("h_name", plot_title, nbins, xmin, xmax);
    TH1D *h_name = new TH1D("h_name", plot_title, 1, 0, 1);
    h_name->Draw();

    //--- draw plot
    TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);
    mainpad->Draw();
    mainpad->cd();
    mainpad->SetLogy();

    bool DF = TString(channel).EqualTo("em");
    if (photon_data_or_mc == "Data") {
        mcstack->Draw("hist");
        mcstack->GetXaxis()->SetTitle(formatted_feature);
        mcstack->GetYaxis()->SetTitle("entries / bin");

        if( !DF ) {
            histograms["zmc"]->Draw("samehist");
            histograms["photon_raw"]->Draw("samehist");
        }
        histograms["data"]->SetLineColor(1); histograms["data"]->SetLineWidth(2); histograms["data"]->SetMarkerStyle(20);
        histograms["data"]->Draw("sameE1");
    }
    else {
        histograms["zmc"]->SetLineColor(1); histograms["zmc"]->SetFillColor(42); histograms["zmc"]->SetLineStyle(1);
        histograms["zmc"]->GetXaxis()->SetTitle(formatted_feature);
        histograms["zmc"]->GetYaxis()->SetTitle("entries / bin");
        histograms["zmc"]->Draw("hist");
        histograms["photon_raw"]->Draw("samehist");
        histograms["photon_reweighted"]->SetLineWidth(1); histograms["photon_reweighted"]->SetLineColor(kRed); histograms["photon_reweighted"]->SetFillStyle(0);
        histograms["photon_reweighted"]->Draw("samehist");
    }

    //--- draw legend and labels
    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
    if (photon_data_or_mc == "Data") {
        leg->AddEntry(histograms["data"],"data","lp");
        if(!DF){
            leg->AddEntry(histograms["zmc"], "Z+jets (from MC)", "f");
            leg->AddEntry(histograms["photon_raw"], "Z+jets (from #gamma+jets, raw)", "f");
            leg->AddEntry(histograms["photon_reweighted"], "Z+jets (from #gamma+jets, reweighted)", "f");
        }
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
    leg->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.03);
    tex->DrawLatex(0.6,0.65,"ATLAS Internal");
    if (photon_data_or_mc == "Data") {
        if(TString(period).Contains("data15-16")) tex->DrawLatex(0.6,0.61,"36 fb^{-1} 2015-2016 data");
        if(TString(period).Contains("data17")) tex->DrawLatex(0.6,0.61,"44 fb^{-1} 2017 data");
        if(TString(period).Contains("data18")) tex->DrawLatex(0.6,0.61,"60 fb^{-1} 2018 data");
    }
    else {
        string mc_period = getMCPeriod(period);
        if(TString(mc_period).Contains("mc16a")) tex->DrawLatex(0.6,0.61,"MC16a");
        if(TString(mc_period).Contains("mc16cd")) tex->DrawLatex(0.6,0.61,"MC16cd");
        if(TString(mc_period).Contains("mc16e")) tex->DrawLatex(0.6,0.61,"MC16e");
    }
    if(TString(channel).Contains("ee")) tex->DrawLatex(0.6,0.57,"ee events");
    if(TString(channel).Contains("em")) tex->DrawLatex(0.6,0.57,"e#mu events");
    if(TString(channel).Contains("mm")) tex->DrawLatex(0.6,0.57,"#mu#mu events");

    //--- draw ratio
    can->cd();
    TPad* ratio_pad = new TPad("ratio_pad","ratio_pad",0.0,0.75,1.0,0.905);
    ratio_pad->Draw();
    ratio_pad->cd();
    ratio_pad->SetGridy();

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
    hratio_unreweighted->Draw("E1");

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
    hratio->Draw("sameE1");

    //--- save plot
    can->Print(plot_name);
}

//---------------
// MAIN FUNCTION
//---------------

void quickDraw(string period, string plot_feature_list, string photon_data_or_mc, string region_list, bool print_photon_yield_only) {
    cout << "period               " << period << endl;
    cout << "photon data          " << photon_data_or_mc << endl;
    cout << endl;

    //--- parse arguments
    vector<string> plot_features = splitStringBySpaces(plot_feature_list);
    vector<string> regions = splitStringBySpaces(region_list);
    vector<string> channels{"ee", "mm"};

    //--- set global options
    gStyle->SetOptStat(0);
    ROOT::EnableImplicitMT();

    //--- get input data in the form of RDataFrames
    unordered_map<string, ROOT::RDataFrame*> RDataFrames = getRDataFrames(period, photon_data_or_mc);
    unordered_map<string, unique_ptr<weightedDataFrame>> WRDataFrames = weightRDataFrames(RDataFrames);

    //--- set up all histograms and fill simulaneously
    auto [plot_region_histograms, control_region_histograms] = setUpHistograms(&WRDataFrames, plot_features, regions, channels);
    histResults filled_histograms = fillHistograms(plot_region_histograms, control_region_histograms, plot_features, regions, channels, photon_data_or_mc);

    //--- print photon yield tables
    printPhotonYieldTables(filled_histograms);

    //--- create histogram stacks and set their plotting options; draw and save plot
    //if (!print_photon_yield_only) {
        //bool DF = TString(channel).EqualTo("em");
        //THStack *mcstack = createStacks(filled_histograms, photon_data_or_mc, DF);
        //TString plot_name = getPlotName(period, channel, plot_feature, photon_data_or_mc, region);
        //makePlot(filled_histograms, mcstack, plot_name, formatted_feature, period, channel, photon_data_or_mc, region);
    //}
}

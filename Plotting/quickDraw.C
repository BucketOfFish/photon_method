#include "../Common/Settings.C"
#include <unordered_map>
#include <boost/algorithm/string.hpp>

using namespace std;

vector<string> splitStringBySpaces(string input) {
    vector<string> output;
    boost::algorithm::split(output, input, boost::is_any_of(" "));
    return output;
}

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

tuple<TString, unordered_map<string, TH1D*>> initializeHistograms(string plot_feature) {
    std::tuple<string, int, float, float> plot_settings;

    if (plot_feature == "met_Et") plot_settings = std::make_tuple("E_{T}^{miss} [GeV]", 30, 0, 300);
    else if (plot_feature == "MET") plot_settings = std::make_tuple("E_{T}^{miss} [GeV]", 30, 0, 300);
    else if (plot_feature == "METl") plot_settings = std::make_tuple("E_{T,||}^{miss} [GeV]", 30, -150, 150);
    else if (plot_feature == "METt") plot_settings = std::make_tuple("E_{T,#perp}^{miss} [GeV]", 30, -150, 150);
    else if (plot_feature == "MET_loose") plot_settings = std::make_tuple("E_{T,loose}^{miss} [GeV]", 20, 0, 200);
    else if (plot_feature == "MET_tight") plot_settings = std::make_tuple("E_{T,tight}^{miss} [GeV]", 20, 0, 200);
    else if (plot_feature == "MET_tighter") plot_settings = std::make_tuple("E_{T,tighter}^{miss} [GeV]", 20, 0, 200);
    else if (plot_feature == "MET_tenacious") plot_settings = std::make_tuple("E_{T,tenacious}^{miss} [GeV]", 20, 0, 200);
    else if (plot_feature == "Ptll") plot_settings = std::make_tuple("p_{T} [GeV]", 20, 0, 100);
    else if (plot_feature == "Z_pt") plot_settings = std::make_tuple("p_{T} [GeV]", 20, 0, 100);
    else if (plot_feature == "nJet30") plot_settings = std::make_tuple("n_{jets}", 6, 2, 8);
    else if (plot_feature == "jet_n") plot_settings = std::make_tuple("n_{jets}", 6, 2, 8);
    else if (plot_feature == "bjet_n") plot_settings = std::make_tuple("n_{b-jets}", 4, 0, 4);
    else if (plot_feature == "HT") plot_settings = std::make_tuple("H_{T}", 20, 0, 1000);
    else if (plot_feature == "mll") plot_settings = std::make_tuple("m_{ll} [GeV]", 30, 0, 300);
    else if (plot_feature == "MT2") plot_settings = std::make_tuple("m_{T2} [GeV]", 20, 0, 200);
    else if (plot_feature == "MT2W") plot_settings = std::make_tuple("m_{T2}^{W} [GeV]", 20, 0, 200);
    else if (plot_feature == "lepPt[0]") plot_settings = std::make_tuple("lep_{p_{T},1} [GeV]", 20, 0, 300);
    else if (plot_feature == "lepPt[1]") plot_settings = std::make_tuple("lep_{p_{T},2} [GeV]", 20, 0, 200);
    else if (plot_feature == "lep_eta[0]") plot_settings = std::make_tuple("lep_{#eta,1} [GeV]", 30, -3, 3);
    else if (plot_feature == "lep_eta[1]") plot_settings = std::make_tuple("lep_{#eta,2} [GeV]", 30, -3, 3);
    else if (plot_feature == "DPhi_METLepLeading") plot_settings = std::make_tuple("#Delta#phi(lep_{1},E_{T}^{miss})", 20, 0, 3.14);
    else if (plot_feature == "DPhi_METLepSecond") plot_settings = std::make_tuple("#Delta#phi(lep_{2},E_{T}^{miss})", 20, 0, 3.14);
    else if (plot_feature == "dPhiMetJet1") plot_settings = std::make_tuple("#Delta#phi(jet_{1},E_{T}^{miss})", 20, 0, 3.14);
    else if (plot_feature == "dPhiMetJet2") plot_settings = std::make_tuple("#Delta#phi(jet_{2},E_{T}^{miss})", 20, 0, 3.14);
    else if (plot_feature == "dPhiMetJet12Min") plot_settings = std::make_tuple("#Delta#phi(jet_{min(1,2)},E_{T}^{miss})", 20, 0, 3.14);
    else {
        cout << "Error! unrecognized variable, need to set binning, quitting! " << plot_feature << endl;
        exit(0);
    }

    TString formatted_feature = std::get<0>(plot_settings).c_str();
    int nbins = std::get<1>(plot_settings);
    float xmin = std::get<2>(plot_settings);
    float xmax = std::get<3>(plot_settings);

    //--- initialize histograms
    TH1D *h_data, *h_photon, *h_photon_reweighted, *h_tt, *h_vv, *h_zmc;

    h_data = new TH1D("h_data", "", nbins, xmin, xmax);
    h_tt = new TH1D("h_tt", "", nbins, xmin, xmax);
    h_vv = new TH1D("h_vv", "", nbins, xmin, xmax);
    h_zmc = new TH1D("h_zmc", "", nbins, xmin, xmax);
    h_photon = new TH1D("h_photon", "", nbins, xmin, xmax);
    h_photon_reweighted = new TH1D("h_photon_reweighted", "", nbins, xmin, xmax);

    unordered_map<string, TH1D*> histograms = {
        {"data", h_data},
        {"tt", h_tt},
        {"vv", h_vv},
        {"zmc", h_zmc},
        {"photon", h_photon},
        {"photon_reweighted", h_photon_reweighted},
    };

    return make_tuple(formatted_feature, histograms);
}

tuple<string, string> getPlotRegions(string channel, string region) {
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

    cout << "bkg selection        " << plot_region.GetTitle() << endl;  
    cout << "bkg weight           " << cuts::bkg_weight.GetTitle() << endl;
    cout << "g selection          " << plot_region.GetTitle() << endl;
    cout << "g weight             " << cuts::photon_weight.GetTitle() << endl;
    cout << "g weight (reweight)  " << cuts::photon_weight_rw.GetTitle() << endl;
    cout << endl;

    return make_tuple(string(plot_region), string(plot_CR));
}

tuple<unordered_map<string, TH1D*>, float> fillHistograms(unordered_map<string, ROOT::RDataFrame*> RDataFrames, unordered_map<string, TH1D*> histograms, string plot_feature, string channel, string region, string photon_data_or_mc, bool DF) {
    auto [plot_region, plot_CR] = getPlotRegions(channel, region);
    cout << "plot region          " << plot_region << endl;
    cout << "normalization reg.   " << plot_CR << endl;
    cout << endl;

    //--- add weight branches
    unordered_map<string, string> plot_weights;
    plot_weights["data"] = "1";
    plot_weights["tt"] = cuts::bkg_weight;
    plot_weights["vv"] = cuts::bkg_weight;
    plot_weights["zmc"] = cuts::bkg_weight;
    plot_weights["photon_raw"] = cuts::photon_weight;
    plot_weights["photon_reweighted"] = cuts::photon_weight_rw;

    unordered_map<string, ROOT::RDF::RResultPtr<TH1D>> plot_region_histograms;
    unordered_map<string, ROOT::RDF::RResultPtr<TH1D>> control_region_histograms;
    for (auto const& [process, dataframe] : RDataFrames) {
        ROOT::RDF::TH1DModel hist_model = ROOT::RDF::TH1DModel(*histograms[process]);
        if (process == "photon") {
            auto weighted_dataframe = dataframe->Define("plot_raw_weight", plot_weights["photon_raw"])
                                             .Define("plot_reweighted_weight", plot_weights["photon_reweighted"]);
            //for (auto plot_feature : plot_features) {
                //string ptn = getPlotTypeName(plot_region, plot_feature);
                plot_region_histograms["photon_raw"] =
                    weighted_dataframe.Filter(plot_region)
                                      .Histo1D(hist_model, plot_feature, "plot_raw_weight");
                plot_region_histograms["photon_reweighted"] =
                    weighted_dataframe.Filter(plot_region)
                                      .Histo1D(hist_model, plot_feature, "plot_reweighted_weight");
                control_region_histograms["photon_raw"] =
                    weighted_dataframe.Filter(plot_CR)
                                      .Histo1D(hist_model, plot_feature, "plot_raw_weight");
                control_region_histograms["photon_reweighted"] =
                    weighted_dataframe.Filter(plot_CR)
                                      .Histo1D(hist_model, plot_feature, "plot_reweighted_weight");
            //}
        }
        else {
            auto weighted_dataframe = dataframe->Define("plot_weight", plot_weights[process]);
            //for (auto plot_feature : plot_features) {
                //string ptn = getPlotTypeName(plot_region, plot_feature);
                plot_region_histograms[process] =
                    weighted_dataframe.Filter(plot_region)
                                      .Histo1D(hist_model, plot_feature, "plot_weight");
                control_region_histograms[process] =
                    weighted_dataframe.Filter(plot_CR)
                                      .Histo1D(hist_model, plot_feature, "plot_weight");
            //}
        }
    }

    cout << "data integral        " << plot_region_histograms["data"]->Integral() << endl;
    cout << "ttbar integral       " << plot_region_histograms["tt"]->Integral() << endl;
    cout << "diboson integral     " << plot_region_histograms["vv"]->Integral() << endl;
    cout << "Z MC integral        " << plot_region_histograms["zmc"]->Integral() << endl;
    cout << "g raw integral       " << plot_region_histograms["photon_raw"]->Integral() << endl;
    cout << "g reweighted int.    " << plot_region_histograms["photon_reweighted"]->Integral() << endl;
    cout << endl;

    float zdata_integral = control_region_histograms["data"]->Integral() - control_region_histograms["tt"]->Integral() - control_region_histograms["vv"]->Integral();
    if (photon_data_or_mc == "MC")
        zdata_integral = control_region_histograms["zmc"]->Integral();
    float SF = zdata_integral / control_region_histograms["photon_raw"]->Integral();
    float SFrw = zdata_integral / control_region_histograms["photon_reweighted"]->Integral();

    cout << "Scaling raw photon data by " << SF << endl;
    cout << "Scaling reweighted photon data by " << SFrw << endl;

    plot_region_histograms["photon_raw"]->Scale(SF);
    plot_region_histograms["photon_reweighted"]->Scale(SFrw);

    float photon_yield = plot_region_histograms["photon_reweighted"]->Integral();
    cout << "Photon yield of " << photon_yield << endl;
    cout << endl;

    for (auto [process, histogram] : plot_region_histograms) {
        if (process == "photon_raw" || process == "photon_reweighted")
            histograms[process] = new TH1D;
        histogram->Copy(*histograms[process]);
    }
    histograms.erase("photon");

    return make_tuple(histograms, photon_yield);
}

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

float quickDraw(string period, string channel, string plot_feature_list, string photon_data_or_mc, string region_list, bool return_photon_yield_only) {
    ROOT::EnableImplicitMT();

    cout << "period               " << period << endl;
    cout << "channel              " << channel << endl;
    cout << "photon data          " << photon_data_or_mc << endl;
    cout << endl;

    //--- parse arguments
    vector<string> plot_features = splitStringBySpaces(plot_feature_list);
    vector<string> regions = splitStringBySpaces(region_list);
    //--- CHECKPOINT: using just first arguments for now
    string plot_feature = plot_features[0];
    string region = regions[0];

    //--- get input data in the form of RDataFrames; define plotting regions; initialize histograms
    gStyle->SetOptStat(0);
    unordered_map<string, ROOT::RDataFrame*> RDataFrames = getRDataFrames(period, photon_data_or_mc);
    auto [formatted_feature, empty_histograms] = initializeHistograms(plot_feature);

    //--- fill histograms and get yields
    bool DF = TString(channel).EqualTo("em");
    auto [filled_histograms, photon_yield] = fillHistograms(RDataFrames, empty_histograms, plot_feature, channel, region, photon_data_or_mc, DF);

    //--- create histogram stacks and set their plotting options; draw and save plot
    if (!return_photon_yield_only) {
        THStack *mcstack = createStacks(filled_histograms, photon_data_or_mc, DF);
        TString plot_name = getPlotName(period, channel, plot_feature, photon_data_or_mc, region);
        makePlot(filled_histograms, mcstack, plot_name, formatted_feature, period, channel, photon_data_or_mc, region);
    }

    return photon_yield;
}

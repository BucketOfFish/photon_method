#include "../Common/Settings.C"
#include <unordered_map>
#include <boost/algorithm/string.hpp>

using namespace std;

vector<string> splitStringBySpaces(string input) {
    vector<string> output;
    boost::algorithm::split(output, input, boost::is_any_of(" "));
    return output;
}

unordered_map<string, TChain*> getTChains(string period, string channel, string photon_data_or_mc, bool DF) {
    //--- load files
    string mc_period = getMCPeriod(period);
    string data_filename= ntuple_path + "bkg_data/" + period + "_bkg.root";
    string tt_filename = ntuple_path + "bkg_mc/" + mc_period + "_ttbar.root";
    string vv_filename = ntuple_path + "bkg_mc/" + mc_period + "_diboson.root";
    string zmc_filename = ntuple_path + "bkg_mc/" + mc_period + "_Zjets.root";
    string photon_filename;
    if (photon_data_or_mc == "MC") photon_filename = reweighting_path + "g_mc/" + mc_period + "_SinglePhoton222_" + channel + ".root";
    else photon_filename = reweighting_path + "g_data/" + period + "_photon_" + channel + ".root";

    cout << "data filename        " << data_filename << endl;
    cout << "ttbar filename       " << tt_filename << endl;
    cout << "diboson filename     " << vv_filename << endl;
    cout << "Z MC filename        " << zmc_filename << endl;
    cout << "photon filename      " << photon_filename << endl;
    cout << endl;

    //--- add files to TChain
    TChain* tch_data = new TChain("BaselineTree"); tch_data->Add(data_filename.c_str());
    TChain* tch_tt = new TChain("BaselineTree"); tch_tt->Add(tt_filename.c_str());
    TChain* tch_vv = new TChain("BaselineTree"); tch_vv->Add(vv_filename.c_str());
    TChain* tch_zmc = new TChain("BaselineTree"); tch_zmc->Add(zmc_filename.c_str());
    TChain* tch_photon = new TChain("BaselineTree"); if (!DF) tch_photon->Add(photon_filename.c_str());

    cout << "Z data entries       " << tch_data->GetEntries() << endl;
    cout << "ttbar entries        " << tch_tt->GetEntries() << endl;
    cout << "diboson entries      " << tch_vv->GetEntries() << endl;
    cout << "Z MC entries         " << tch_zmc->GetEntries() << endl;
    cout << "photon entries       " << tch_photon->GetEntries() << endl;
    cout << endl;

    unordered_map<string, TChain*> TChains = {
        {"data", tch_data},
        {"tt", tch_tt},
        {"vv", tch_vv},
        {"zmc", tch_zmc},
        {"photon", tch_photon},
    };

    return TChains;
}

tuple<TCut, TCut> getPlotRegions(string channel, string region, string additional_cut) {
    if (cuts::plot_regions.count(region) == 0) {
        cout << "Unrecognized region! Exiting." << endl;
        exit(0);
    }
    TCut plot_region = cuts::plot_regions[region];

    if (additional_cut != "1") plot_region += TCut(additional_cut.c_str());

    if (TString(channel).EqualTo("ee")) plot_region += cuts::ee;
    else if (TString(channel).EqualTo("mm")) plot_region += cuts::mm;
    else if (TString(channel).EqualTo("em")) plot_region += cuts::em;
    else {
        cout << "Unrecognized channel! quitting   " << channel << endl;
        exit(0);
    }

    TCut plot_CR = plot_region + cuts::CR;
    plot_region += cuts::plot_region_met_portions[region];

    cout << "Z selection          " << plot_region.GetTitle() << endl;  
    cout << "Z weight             " << cuts::bkg_weight.GetTitle() << endl;
    cout << "g selection          " << plot_region.GetTitle() << endl;
    cout << "g weight             " << cuts::photon_weight.GetTitle() << endl;
    cout << "g weight (reweight)  " << cuts::photon_weight_rw.GetTitle() << endl;
    cout << endl;

    return make_tuple(plot_region, plot_CR);
}

tuple<TString, unordered_map<string, TH1F*>> initializeHistograms(string plot_feature) {
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
    TH1F *h_data, *h_photon, *h_photon_reweighted, *h_tt, *h_vv, *h_zmc;

    h_data = new TH1F("h_data", "", nbins, xmin, xmax);
    h_tt = new TH1F("h_tt", "", nbins, xmin, xmax);
    h_vv = new TH1F("h_vv", "", nbins, xmin, xmax);
    h_zmc = new TH1F("h_zmc", "", nbins, xmin, xmax);
    h_photon = new TH1F("h_photon", "", nbins, xmin, xmax);
    h_photon_reweighted = new TH1F("h_photon_reweighted", "", nbins, xmin, xmax);

    unordered_map<string, TH1F*> histograms = {
        {"data", h_data},
        {"tt", h_tt},
        {"vv", h_vv},
        {"zmc", h_zmc},
        {"photon", h_photon},
        {"photon_reweighted", h_photon_reweighted},
    };

    return make_tuple(formatted_feature, histograms);
}

float fillHistograms(unordered_map<string, TChain*> TChains, unordered_map<string, TH1F*> histograms, string plot_feature, TCut plot_region, TCut plot_CR, string photon_data_or_mc, bool DF) {
    TChains["data"]->Draw(Form("%s>>h_data", plot_feature.c_str()), plot_region, "goff");
    TChains["tt"]->Draw(Form("%s>>h_tt", plot_feature.c_str()), plot_region*cuts::bkg_weight, "goff");
    TChains["vv"]->Draw(Form("%s>>h_vv", plot_feature.c_str()), plot_region*cuts::bkg_weight, "goff");
    if (!DF) {
        TChains["zmc"]->Draw(Form("%s>>h_zmc", plot_feature.c_str()), plot_region*cuts::bkg_weight, "goff");
        TChains["photon"]->Draw(Form("%s>>h_photon", plot_feature.c_str()), plot_region*cuts::photon_weight, "goff");
        TChains["photon"]->Draw(Form("%s>>h_photon_reweighted", plot_feature.c_str()), plot_region*cuts::photon_weight_rw, "goff");
    }

    cout << "data integral        " << histograms["data"]->Integral() << endl;
    cout << "tt integral          " << histograms["tt"]->Integral() << endl;
    cout << "VV integral          " << histograms["vv"]->Integral() << endl;
    cout << "Z MC integral        " << histograms["zmc"]->Integral() << endl;
    cout << "g raw integral       " << histograms["photon"]->Integral() << endl;
    cout << "g reweighted int.    " << histograms["photon_reweighted"]->Integral() << endl;
    cout << endl;

    //--- normalize Z to CR
    cout << "normalize to CR " << cuts::CR.GetTitle() << endl;

    TH1F* h_data_cr = new TH1F("h_data_cr", "", 1, 0, 1);
    TH1F* h_tt_cr = new TH1F("h_tt_cr", "", 1, 0, 1);
    TH1F* h_vv_cr = new TH1F("h_vv_cr", "", 1, 0, 1);
    TH1F* h_zmc_cr = new TH1F("h_zmc_cr", "", 1, 0, 1);
    TH1F* h_photon_cr = new TH1F("h_photon_cr", "", 1, 0, 1);
    TH1F* h_photon_reweighted_cr = new TH1F("h_photon_reweighted_cr", "", 1, 0, 1);

    TChains["data"]->Draw("0.5>>h_data_cr", plot_CR*cuts::bkg_weight, "goff");
    TChains["tt"]-> Draw("0.5>>h_tt_cr", plot_CR*cuts::bkg_weight, "goff");
    TChains["vv"]-> Draw("0.5>>h_vv_cr", plot_CR*cuts::bkg_weight, "goff");
    TChains["zmc"]->Draw("0.5>>h_zmc_cr", plot_CR*cuts::bkg_weight, "goff");
    TChains["photon"]->Draw("0.5>>h_photon_cr", plot_CR*cuts::photon_weight, "goff");
    TChains["photon"]->Draw("0.5>>h_photon_reweighted_cr", plot_CR*cuts::photon_weight_rw, "goff");

    float SF = (h_data_cr->Integral() - h_tt_cr->Integral() - h_vv_cr->Integral()) / h_photon_cr->Integral();
    float SFrw = (h_data_cr->Integral() - h_tt_cr->Integral() - h_vv_cr->Integral()) / h_photon_reweighted_cr->Integral();

    if (photon_data_or_mc == "MC") {
        SF = h_zmc_cr->Integral() / h_photon_cr->Integral();
        SFrw = h_zmc_cr->Integral() / h_photon_reweighted_cr->Integral();
    }

    cout << "Scale raw photon data by " << SF << endl;
    cout << "Scale reweighted photon data by " << SFrw << endl;

    histograms["photon"]->Scale(SF);
    histograms["photon_reweighted"]->Scale(SFrw);

    //--- print MET integrals
    cout << "MET100-150" << endl;
    cout << "2L data                " << histograms["data"]->Integral(11,15) << endl;
    cout << "VV MC                  " << histograms["vv"]->Integral(11,15) << endl;
    cout << "tt MC                  " << histograms["tt"]->Integral(11,15) << endl;
    cout << "Z+jets MC              " << histograms["zmc"]->Integral(11,15) << endl;
    cout << "g data (reweighted)    " << histograms["photon_reweighted"]->Integral(11,15) << endl;
    cout << "g data (raw)           " << histograms["photon"]->Integral(11,15) << endl;

    cout << "MET150-200" << endl;
    cout << "2L data                " << histograms["data"]->Integral(16,21) << endl;
    cout << "VV MC                  " << histograms["vv"]->Integral(16,21) << endl;
    cout << "tt MC                  " << histograms["tt"]->Integral(16,21) << endl;
    cout << "Z+jets MC              " << histograms["zmc"]->Integral(16,21) << endl;
    cout << "g data (reweighted)    " << histograms["photon_reweighted"]->Integral(16,21) << endl;
    cout << "g data (raw)           " << histograms["photon"]->Integral(16,21) << endl;

    float photon_yield = histograms["photon_reweighted"]->Integral();
    cout << endl;
    cout << "Photon yield of " << photon_yield << endl;
    cout << endl;

    return photon_yield;
}

THStack* createStacks(unordered_map<string, TH1F*> histograms, string photon_data_or_mc, bool DF) {
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
        histograms["photon"]->Add(histograms["tt"]); histograms["photon"]->Add(histograms["vv"]);
        histograms["zmc"]->Add(histograms["tt"]); histograms["zmc"]->Add(histograms["vv"]);
    }

    //--- set plotting options
    histograms["photon"]->SetLineColor(4); histograms["photon"]->SetLineWidth(1); histograms["photon"]->SetLineStyle(2);
    histograms["zmc"]->SetLineColor(2); histograms["zmc"]->SetLineWidth(1); histograms["zmc"]->SetLineStyle(7);

    //--- turn on overflow bin
    histograms["zmc"]->GetXaxis()->SetRange(0, histograms["zmc"]->GetNbinsX() + 1);
    histograms["photon"]->GetXaxis()->SetRange(0, histograms["photon"]->GetNbinsX() + 1);
    histograms["photon_reweighted"]->GetXaxis()->SetRange(0, histograms["photon_reweighted"]->GetNbinsX() + 1);

    return mcstack;
}

TString getPlotName(string period, string channel, string plot_feature, string photon_data_or_mc, string additional_cut, string region) {
    TString plot_name;
    if (photon_data_or_mc == "Data")
        plot_name = Form("%s/%s_%s_%s_%s", plots_path.c_str(), period.c_str(), channel.c_str(), plot_feature.c_str(), region.c_str());
    else
        plot_name = Form("%s/%s_%s_%s_%s", plots_path.c_str(), getMCPeriod(period).c_str(), channel.c_str(), plot_feature.c_str(), region.c_str());
    if (additional_cut != "1")
        plot_name += ("_" + additional_cut).c_str();
    plot_name += ".eps";
    return plot_name;
}

void makePlot(unordered_map<string, TH1F*> histograms, THStack *mcstack, TString plot_name, TString formatted_feature, string period, string channel, string photon_data_or_mc, string additional_cut, string region, bool DF) {
    //--- draw title
    TCanvas *can = new TCanvas("can","can",600,600);
    can->cd();
    TPad* namepad = new TPad("namepad","namepad",0.0,0.0,1.0,1.0);
    namepad->Draw();
    namepad->cd();
    TString plot_title = formatted_feature + " in " + region;
    if (additional_cut != "1")
        plot_title += (", " + additional_cut).c_str();
    //TH1F *h_name = new TH1F("h_name", plot_title, nbins, xmin, xmax);
    TH1F *h_name = new TH1F("h_name", plot_title, 1, 0, 1);
    h_name->Draw();

    //--- draw plot
    TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);
    mainpad->Draw();
    mainpad->cd();
    mainpad->SetLogy();

    if (photon_data_or_mc == "Data") {
        mcstack->Draw("hist");
        mcstack->GetXaxis()->SetTitle(formatted_feature);
        mcstack->GetYaxis()->SetTitle("entries / bin");

        if( !DF ) {
            histograms["zmc"]->Draw("samehist");
            histograms["photon"]->Draw("samehist");
        }
        histograms["data"]->SetLineColor(1); histograms["data"]->SetLineWidth(2); histograms["data"]->SetMarkerStyle(20);
        histograms["data"]->Draw("sameE1");
    }
    else {
        histograms["zmc"]->SetLineColor(1); histograms["zmc"]->SetFillColor(42); histograms["zmc"]->SetLineStyle(1);
        histograms["zmc"]->GetXaxis()->SetTitle(formatted_feature);
        histograms["zmc"]->GetYaxis()->SetTitle("entries / bin");
        histograms["zmc"]->Draw("hist");
        histograms["photon"]->Draw("samehist");
        histograms["photon_reweighted"]->SetLineWidth(1); histograms["photon_reweighted"]->SetLineColor(kRed); histograms["photon_reweighted"]->SetFillStyle(0);
        histograms["photon_reweighted"]->Draw("samehist");
    }

    //--- draw legend and labels
    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
    if (photon_data_or_mc == "Data") {
        leg->AddEntry(histograms["data"],"data","lp");
        if(!DF){
            leg->AddEntry(histograms["zmc"], "Z+jets (from MC)", "f");
            leg->AddEntry(histograms["photon"], "Z+jets (from #gamma+jets, raw)", "f");
            leg->AddEntry(histograms["photon_reweighted"], "Z+jets (from #gamma+jets, reweighted)", "f");
        }
        leg->AddEntry(histograms["vv"], "VV", "f");
        leg->AddEntry(histograms["tt"], "t#bar{t}+tW", "f");
    }
    else {
        leg->AddEntry(histograms["zmc"], "Z+jets (from MC)", "f");
        leg->AddEntry(histograms["photon"], "Z+jets (from #gamma+jets, raw)", "f");
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

    TH1F *hratio, *hratio_unreweighted, *hmctot, *hmctot_unreweighted;

    if (photon_data_or_mc == "MC") {
        hratio = (TH1F*) histograms["zmc"]->Clone("hratio");
        hratio_unreweighted = (TH1F*) histograms["zmc"]->Clone("hratio");
        hmctot = (TH1F*) histograms["photon_reweighted"]->Clone("hmctot");
        hmctot_unreweighted = (TH1F*) histograms["photon"]->Clone("hmctot");
    }
    else {
        hratio = (TH1F*) histograms["data"]->Clone("hratio");
        hratio_unreweighted = (TH1F*) histograms["data"]->Clone("hratio");
        hmctot = (TH1F*) histograms["photon_reweighted"]->Clone("hmctot");
        hmctot->Add(histograms["tt"]);
        hmctot->Add(histograms["vv"]);
        hmctot_unreweighted = (TH1F*) histograms["photon"]->Clone("hmctot");
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

float quickDraw(string period, string channel, string plot_feature_list, string photon_data_or_mc, string region_list, string additional_cut, bool return_photon_yield_only) {
    ROOT::EnableImplicitMT();

    bool DF = TString(channel).EqualTo("em");
    gStyle->SetOptStat(0);

    cout << "period               " << period << endl;
    cout << "channel              " << channel << endl;
    cout << "DF?                  " << DF << endl;
    cout << "photon data          " << photon_data_or_mc << endl;
    cout << endl;

    //--- parse arguments
    vector<string> plot_features = splitStringBySpaces(plot_feature_list);
    vector<string> regions = splitStringBySpaces(region_list);

    //--- CHECKPOINT: using just first arguments for now
    string plot_feature = plot_features[0];
    string region = regions[0];

    //--- get TChains; define plotting regions; initialize histograms
    auto [plot_region, plot_CR] = getPlotRegions(channel, region, additional_cut);
    unordered_map<string, TChain*> TChains = getTChains(period, channel, photon_data_or_mc, DF);
    auto [formatted_feature, histograms] = initializeHistograms(plot_feature);

    //--- fill histograms and get yields
    float photon_yield = fillHistograms(TChains, histograms, plot_feature, plot_region, plot_CR, photon_data_or_mc, DF);
    if (return_photon_yield_only)
        return histograms["photon_reweighted"]->Integral();

    //--- create histogram stacks and set their plotting options; draw and save plot
    THStack *mcstack = createStacks(histograms, photon_data_or_mc, DF);
    TString plot_name = getPlotName(period, channel, plot_feature, photon_data_or_mc, additional_cut, region);
    makePlot(histograms, mcstack, plot_name, formatted_feature, period, channel, photon_data_or_mc, additional_cut, region, DF);

    return photon_yield;
}

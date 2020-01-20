#include "../../Common/Settings.C"
#include <boost/algorithm/string/replace.hpp>

using namespace std;

void Z_METl_asymmetry(string period="data18", string channel="ee" , string plot_feature="METl", string region="SR", string additional_cut="1") {

    gStyle->SetOptStat(0);

    cout << "period               " << period << endl;
    cout << "channel              " << channel << endl;

    //--- load files
    string mc_period = "";
    if (TString(period).Contains("data15-16")) mc_period = "mc16a";
    else if (TString(period).Contains("data17")) mc_period = "mc16cd";
    else if (TString(period).Contains("data18")) mc_period = "mc16e";

    string zmc_filename = ntuple_path + "bkg_mc/" + mc_period + "_Zjets.root";

    cout << "Z MC filename        " << zmc_filename << endl;
    cout << "" << endl;

    //--- add files to TChain
    TChain* tch_zmc = new TChain("BaselineTree"); tch_zmc->Add(zmc_filename.c_str());

    cout << "Z MC entries         " << tch_zmc->GetEntries() << endl;

    //--- define selections
    TCut plot_region;
    if (region == "CR") plot_region = cuts::CR;
    else if (region == "baseline") plot_region = cuts::baseline;
    else if (region == "reweight") plot_region = cuts::reweight_region;
    else if (region == "VR") plot_region = cuts::VR;
    else if (region == "SR") plot_region = cuts::SR;
    else if (region == "VRcom") plot_region = cuts::VRcom;
    else if (region == "SRZ2016") plot_region = cuts::SRZ2016;
    else if (region == "SRlow2016") plot_region = cuts::SRlow2016;
    else if (region == "SRmed2016") plot_region = cuts::SRmed2016;
    else if (region == "SRhigh2016") plot_region = cuts::SRhigh2016;
    else {
        cout << "Unrecognized region! Exiting." << endl;
        exit(0);
    }

    if (additional_cut != "1") plot_region += TCut(additional_cut.c_str());
    if (TString(channel).EqualTo("ee")) plot_region += cuts::ee;
    else if (TString(channel).EqualTo("mm")) plot_region += cuts::mm;
    else if (TString(channel).EqualTo("em")) plot_region += cuts::em;
    else {
        cout << "Unrecognized channel! quitting   " << channel << endl;
        exit(0);
    }
    plot_region += TCut(additional_cut.c_str());

    cout << "Z selection          " << plot_region.GetTitle() << endl;  
    cout << "Z weight             " << cuts::bkg_weight.GetTitle() << endl;
    cout << "g selection          " << plot_region.GetTitle() << endl;
    cout << "g weight             " << cuts::photon_weight.GetTitle() << endl;
    cout << "g weight (reweight)  " << cuts::photon_weight_rw.GetTitle() << endl;

    //--- set histogram binning
    std::tuple<string, int, float, float> plot_settings;

    if (plot_feature == "met_Et") plot_settings = std::make_tuple("E_{T}^{miss} [GeV]", 30, 0, 300);
    else if (plot_feature == "MET") plot_settings = std::make_tuple("E_{T}^{miss} [GeV]", 30, 0, 300);
    else if (plot_feature == "METl") plot_settings = std::make_tuple("E_{T,||}^{miss} [GeV]", 30, -150, 150);
    else if (plot_feature == "METt") plot_settings = std::make_tuple("E_{T,#perp}^{miss} [GeV]", 30, -150, 150);

    TString formatted_feature = std::get<0>(plot_settings).c_str();
    int nbins = std::get<1>(plot_settings);
    float xmin = std::get<2>(plot_settings);
    float xmax = std::get<3>(plot_settings);

    //--- initialize histograms
    TH1F *h_zmc = new TH1F("h_zmc", "", nbins, xmin, xmax);
    TH1F *h_zmc_flipped = new TH1F("h_zmc_flipped", "", nbins, xmin, xmax);

    //--- draw histograms
    tch_zmc->Draw(Form("%s>>h_zmc", plot_feature.c_str()), plot_region*cuts::bkg_weight, "goff");

    cout << "" << endl;
    cout << "Z MC integral        " << h_zmc->Integral() << endl;
    cout << "" << endl;

    //--- flip plot
    float overlap_integral = 0.0;
    float total_integral = 0.0;
    for (int i=0; i<=nbins+1; i++) {
        h_zmc_flipped->SetBinContent(i, h_zmc->GetBinContent(nbins+1-i));
        overlap_integral += min(h_zmc->GetBinContent(i), h_zmc->GetBinContent(nbins+1-i));
        total_integral += h_zmc->GetBinContent(i);
    }
    float overlap_ratio = overlap_integral / total_integral;

    //--- set plotting options
    h_zmc->SetLineColor(2); h_zmc->SetLineWidth(1); h_zmc->SetLineStyle(7);
    h_zmc_flipped->SetLineColor(kRed); h_zmc_flipped->SetFillStyle(0);

    //--- turn on overflow bin
    h_zmc->GetXaxis()->SetRange(0, h_zmc->GetNbinsX() + 1);
    h_zmc_flipped->GetXaxis()->SetRange(0, h_zmc_flipped->GetNbinsX() + 1);

    //--- draw title
    TCanvas *can = new TCanvas("can","can",600,600);
    can->cd();
    TPad* namepad = new TPad("namepad","namepad",0.0,0.0,1.0,1.0);
    namepad->Draw();
    namepad->cd();
    TString plot_title = formatted_feature + " in " + region;
    if (additional_cut != "1")
        plot_title += (", " + additional_cut).c_str();
    TH1F *h_name = new TH1F("h_name", plot_title, nbins, xmin, xmax);
    h_name->Draw();

    //--- draw plot
    TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,1.0);
    mainpad->Draw();
    mainpad->cd();
    mainpad->SetLogy();

    h_zmc->SetLineColor(1); h_zmc->SetFillColor(42); h_zmc->SetLineStyle(1);
    h_zmc->GetXaxis()->SetTitle(formatted_feature);
    h_zmc->GetYaxis()->SetTitle("entries / bin");
    h_zmc->Draw("hist");
    h_zmc_flipped->Draw("samehist");

    //--- draw legend and labels
    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
    leg->AddEntry(h_zmc, "Z+jets (from MC)", "f");
    leg->AddEntry(h_zmc_flipped, "Z+jets flipped", "f");

    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.03);
    tex->DrawLatex(0.6,0.65,"ATLAS Internal");
    if(TString(mc_period).Contains("mc16a")) tex->DrawLatex(0.6,0.61,"MC16a");
    if(TString(mc_period).Contains("mc16cd")) tex->DrawLatex(0.6,0.61,"MC16cd");
    if(TString(mc_period).Contains("mc16e")) tex->DrawLatex(0.6,0.61,"MC16e");
    if(TString(channel).Contains("ee")) tex->DrawLatex(0.6,0.57,"ee events");
    if(TString(channel).Contains("em")) tex->DrawLatex(0.6,0.57,"e#mu events");
    if(TString(channel).Contains("mm")) tex->DrawLatex(0.6,0.57,"#mu#mu events");
    tex->DrawLatex(0.6,0.53,TString("Overlap of " + to_string(overlap_ratio)));

    //--- save plot
    TString plot_name;
    plot_name = Form("Plots/Z_METl_asymmetry_%s_%s_%s_%s", mc_period.c_str(), channel.c_str(), plot_feature.c_str(), region.c_str());
    if (additional_cut != "1")
        plot_name += ("_" + additional_cut).c_str();
    plot_name += ".eps";
    can->Print(plot_name);
}

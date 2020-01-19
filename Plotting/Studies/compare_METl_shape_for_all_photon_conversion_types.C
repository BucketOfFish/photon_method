#include "../../Common/Settings.C"
#include <boost/algorithm/string/replace.hpp>

using namespace std;

void compare_METl_shape_for_all_photon_conversion_types(string period="data18", string channel="mm" , string plot_feature="met_Et", string region="SR", string additional_cut="1") {

    gStyle->SetOptStat(0);

    cout << "period               " << period << endl;
    cout << "channel              " << channel << endl;

    //--- load files
    string mc_period = "";
    if (TString(period).Contains("data15-16")) mc_period = "mc16a";
    else if (TString(period).Contains("data17")) mc_period = "mc16cd";
    else if (TString(period).Contains("data18")) mc_period = "mc16e";

    string zmc_filename = ntuple_path + "bkg_mc/" + mc_period + "_Zjets.root";
    string photon_filename = reweighting_path + "g_mc/" + mc_period + "_SinglePhoton222_" + channel + ".root";

    cout << "Z MC filename        " << zmc_filename << endl;
    cout << "photon filename      " << photon_filename << endl;
    cout << "" << endl;

    //--- add files to TChain
    TChain* tch_zmc = new TChain("BaselineTree"); tch_zmc->Add(zmc_filename.c_str());
    TChain* tch_photon = new TChain("BaselineTree"); tch_photon->Add(photon_filename.c_str());

    cout << "Z MC entries         " << tch_zmc->GetEntries() << endl;
    cout << "photon entries       " << tch_photon->GetEntries() << endl;

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
    TH1F *h_photon, *h_photon_reweighted, *h_zmc;

    h_zmc = new TH1F("h_zmc", "", nbins, xmin, xmax);
    h_photon = new TH1F("h_photon", "", nbins, xmin, xmax);
    h_photon_reweighted = new TH1F("h_photon_reweighted", "", nbins, xmin, xmax);

    TH1F *h_photon_rw_pct_0 = new TH1F("h_photon_rw_pct_0", "", nbins, xmin, xmax);
    TH1F *h_photon_rw_pct_1 = new TH1F("h_photon_rw_pct_1", "", nbins, xmin, xmax);
    TH1F *h_photon_rw_pct_2 = new TH1F("h_photon_rw_pct_2", "", nbins, xmin, xmax);
    TH1F *h_photon_rw_pct_3 = new TH1F("h_photon_rw_pct_3", "", nbins, xmin, xmax);
    TH1F *h_photon_rw_pct_4 = new TH1F("h_photon_rw_pct_4", "", nbins, xmin, xmax);
    TH1F *h_photon_rw_pct_5 = new TH1F("h_photon_rw_pct_5", "", nbins, xmin, xmax);

    //--- draw histograms
    tch_zmc->Draw(Form("%s>>h_zmc", plot_feature.c_str()), plot_region*cuts::bkg_weight, "goff");
    tch_photon->Draw(Form("%s>>h_photon", plot_feature.c_str()), plot_region*cuts::photon_weight, "goff");
    tch_photon->Draw(Form("%s>>h_photon_reweighted", plot_feature.c_str()), plot_region*cuts::photon_weight_rw, "goff");
    tch_photon->Draw(Form("%s>>h_photon_rw_pct_0", plot_feature.c_str()), plot_region*cuts::photon_weight_rw*"PhotonConversionType==0", "goff");
    tch_photon->Draw(Form("%s>>h_photon_rw_pct_1", plot_feature.c_str()), plot_region*cuts::photon_weight_rw*"PhotonConversionType==1", "goff");
    tch_photon->Draw(Form("%s>>h_photon_rw_pct_2", plot_feature.c_str()), plot_region*cuts::photon_weight_rw*"PhotonConversionType==2", "goff");
    tch_photon->Draw(Form("%s>>h_photon_rw_pct_3", plot_feature.c_str()), plot_region*cuts::photon_weight_rw*"PhotonConversionType==3", "goff");
    tch_photon->Draw(Form("%s>>h_photon_rw_pct_4", plot_feature.c_str()), plot_region*cuts::photon_weight_rw*"PhotonConversionType==4", "goff");
    tch_photon->Draw(Form("%s>>h_photon_rw_pct_5", plot_feature.c_str()), plot_region*cuts::photon_weight_rw*"PhotonConversionType==5", "goff");

    cout << "" << endl;
    cout << "Z MC integral        " << h_zmc->Integral() << endl;
    cout << "g raw integral       " << h_photon->Integral() << endl;
    cout << "g reweighted int.    " << h_photon_reweighted->Integral() << endl;
    cout << "" << endl;

    //--- normalize Z to MET<60 GeV region
    cout << "normalize to CR " << cuts::CR.GetTitle() << endl;

    TH1F* h_zmc_cr = new TH1F("h_zmc_cr", "", 1, 0, 1);
    TH1F* h_photon_cr = new TH1F("h_photon_cr", "", 1, 0, 1);
    TH1F* h_photon_reweighted_cr = new TH1F("h_photon_reweighted_cr", "", 1, 0, 1);
    TH1F* h_photon_rw_pct_0_cr = new TH1F("h_photon_rw_pct_0_cr", "", 1, 0, 1);
    TH1F* h_photon_rw_pct_1_cr = new TH1F("h_photon_rw_pct_1_cr", "", 1, 0, 1);
    TH1F* h_photon_rw_pct_2_cr = new TH1F("h_photon_rw_pct_2_cr", "", 1, 0, 1);
    TH1F* h_photon_rw_pct_3_cr = new TH1F("h_photon_rw_pct_3_cr", "", 1, 0, 1);
    TH1F* h_photon_rw_pct_4_cr = new TH1F("h_photon_rw_pct_4_cr", "", 1, 0, 1);
    TH1F* h_photon_rw_pct_5_cr = new TH1F("h_photon_rw_pct_5_cr", "", 1, 0, 1);

    TCut current_CR = plot_region + cuts::CR;
    tch_zmc->Draw("0.5>>h_zmc_cr", current_CR*cuts::bkg_weight, "goff");
    tch_photon->Draw("0.5>>h_photon_cr", current_CR*cuts::photon_weight, "goff");
    tch_photon->Draw("0.5>>h_photon_reweighted_cr", current_CR*cuts::photon_weight_rw, "goff");
    tch_photon->Draw("0.5>>h_photon_rw_pct_0_cr", current_CR*cuts::photon_weight_rw*"PhotonConversionType==0", "goff");
    tch_photon->Draw("0.5>>h_photon_rw_pct_1_cr", current_CR*cuts::photon_weight_rw*"PhotonConversionType==1", "goff");
    tch_photon->Draw("0.5>>h_photon_rw_pct_2_cr", current_CR*cuts::photon_weight_rw*"PhotonConversionType==2", "goff");
    tch_photon->Draw("0.5>>h_photon_rw_pct_3_cr", current_CR*cuts::photon_weight_rw*"PhotonConversionType==3", "goff");
    tch_photon->Draw("0.5>>h_photon_rw_pct_4_cr", current_CR*cuts::photon_weight_rw*"PhotonConversionType==4", "goff");
    tch_photon->Draw("0.5>>h_photon_rw_pct_5_cr", current_CR*cuts::photon_weight_rw*"PhotonConversionType==5", "goff");

    float SF = h_zmc_cr->Integral() / h_photon_cr->Integral();
    float SFrw = h_zmc_cr->Integral() / h_photon_reweighted_cr->Integral();
    float SFrw_pct_0 = h_zmc_cr->Integral() / h_photon_rw_pct_0_cr->Integral();
    float SFrw_pct_1 = h_zmc_cr->Integral() / h_photon_rw_pct_1_cr->Integral();
    float SFrw_pct_2 = h_zmc_cr->Integral() / h_photon_rw_pct_2_cr->Integral();
    float SFrw_pct_3 = h_zmc_cr->Integral() / h_photon_rw_pct_3_cr->Integral();
    float SFrw_pct_4 = h_zmc_cr->Integral() / h_photon_rw_pct_4_cr->Integral();
    float SFrw_pct_5 = h_zmc_cr->Integral() / h_photon_rw_pct_5_cr->Integral();

    cout << "Scale raw photon data by " << SF << endl;
    cout << "Scale reweighted photon data by " << SFrw << endl;

    h_photon->Scale(SF);
    h_photon_reweighted->Scale(SFrw);
    h_photon_rw_pct_0->Scale(SFrw_pct_0);
    h_photon_rw_pct_1->Scale(SFrw_pct_1);
    h_photon_rw_pct_2->Scale(SFrw_pct_2);
    h_photon_rw_pct_3->Scale(SFrw_pct_3);
    h_photon_rw_pct_4->Scale(SFrw_pct_4);
    h_photon_rw_pct_5->Scale(SFrw_pct_5);

    //--- print MET integrals
    cout << "MET100-150" << endl;
    cout << "Z+jets MC              " << h_zmc->Integral(11,15) << endl;
    cout << "g data (reweighted)    " << h_photon_reweighted->Integral(11,15) << endl;
    cout << "g data (raw)           " << h_photon->Integral(11,15) << endl;

    cout << "MET150-200" << endl;
    cout << "Z+jets MC              " << h_zmc->Integral(16,21) << endl;
    cout << "g data (reweighted)    " << h_photon_reweighted->Integral(16,21) << endl;
    cout << "g data (raw)           " << h_photon->Integral(16,21) << endl;

    //--- set plotting options
    h_zmc->SetLineColor(2); h_zmc->SetLineWidth(1); h_zmc->SetLineStyle(7);
    h_photon->SetLineColor(4); h_photon->SetLineWidth(1); h_photon->SetLineStyle(2);
    h_photon_reweighted->SetLineColor(8); h_photon_reweighted->SetFillStyle(0);
    h_photon_rw_pct_0->SetLineColor(2); h_photon_rw_pct_0->SetFillStyle(0);
    h_photon_rw_pct_1->SetLineColor(3); h_photon_rw_pct_1->SetFillStyle(0);
    h_photon_rw_pct_2->SetLineColor(4); h_photon_rw_pct_2->SetFillStyle(0);
    h_photon_rw_pct_3->SetLineColor(5); h_photon_rw_pct_3->SetFillStyle(0);
    h_photon_rw_pct_4->SetLineColor(6); h_photon_rw_pct_4->SetFillStyle(0);
    h_photon_rw_pct_5->SetLineColor(7); h_photon_rw_pct_5->SetFillStyle(0);

    //--- turn on overflow bin
    h_zmc->GetXaxis()->SetRange(0, h_zmc->GetNbinsX() + 1);
    h_photon->GetXaxis()->SetRange(0, h_photon->GetNbinsX() + 1);
    h_photon_reweighted->GetXaxis()->SetRange(0, h_photon_reweighted->GetNbinsX() + 1);
    h_photon_rw_pct_0->GetXaxis()->SetRange(0, h_photon_rw_pct_0->GetNbinsX() + 1);
    h_photon_rw_pct_1->GetXaxis()->SetRange(0, h_photon_rw_pct_1->GetNbinsX() + 1);
    h_photon_rw_pct_2->GetXaxis()->SetRange(0, h_photon_rw_pct_2->GetNbinsX() + 1);
    h_photon_rw_pct_3->GetXaxis()->SetRange(0, h_photon_rw_pct_3->GetNbinsX() + 1);
    h_photon_rw_pct_4->GetXaxis()->SetRange(0, h_photon_rw_pct_4->GetNbinsX() + 1);
    h_photon_rw_pct_5->GetXaxis()->SetRange(0, h_photon_rw_pct_5->GetNbinsX() + 1);

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
    h_photon->Draw("samehist");
    h_photon_reweighted->Draw("samehist");
    h_photon_rw_pct_0->Draw("samehist");
    h_photon_rw_pct_1->Draw("samehist");
    h_photon_rw_pct_2->Draw("samehist");
    h_photon_rw_pct_3->Draw("samehist");
    h_photon_rw_pct_4->Draw("samehist");
    h_photon_rw_pct_5->Draw("samehist");

    //--- draw legend and labels
    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
    leg->AddEntry(h_zmc, "Z+jets (from MC)", "f");
    leg->AddEntry(h_photon, "Z+jets (from #gamma+jets, raw)", "f");
    leg->AddEntry(h_photon_reweighted, "Z+jets (from #gamma+jets, reweighted)", "f");
    leg->AddEntry(h_photon_rw_pct_0, "Z+jets (from #gamma+jets, PCT 0)", "f");
    leg->AddEntry(h_photon_rw_pct_1, "Z+jets (from #gamma+jets, PCT 1)", "f");
    leg->AddEntry(h_photon_rw_pct_2, "Z+jets (from #gamma+jets, PCT 2)", "f");
    leg->AddEntry(h_photon_rw_pct_3, "Z+jets (from #gamma+jets, PCT 3)", "f");
    leg->AddEntry(h_photon_rw_pct_4, "Z+jets (from #gamma+jets, PCT 4)", "f");
    leg->AddEntry(h_photon_rw_pct_5, "Z+jets (from #gamma+jets, PCT 5)", "f");

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

    //--- save plot
    TString plot_name;
    plot_name = Form("Plots/photon_conversion_type_comparison_%s_%s_%s_%s", mc_period.c_str(), channel.c_str(), plot_feature.c_str(), region.c_str());
    if (additional_cut != "1")
        plot_name += ("_" + additional_cut).c_str();
    plot_name += ".eps";
    can->Print(plot_name);
}

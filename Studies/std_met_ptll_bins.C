#include "../Common/Settings.C"
#include "../Common/CommonLibraries.C"
#include "../Common/CommonFunctions.C"

using namespace std;

//--- run with:
//--- root -l -b -q 'std_met_ptll_bins.C("mc16a", "ee", "met_Et", "Ptll")'

//--- Plots standard deviation of MET in Ptll bins

void std_met_ptll_bins(string mc_period, string channel, string plot_feature, string bin_feature) {

    gStyle->SetOptStat(0);

    //--- load files
    string zmc_filename = ntuple_path + "bkg_mc/" + mc_period + "_Zjets.root";
    string photon_filename = reweighting_path + "g_mc/" + mc_period + "_SinglePhoton222_" + channel + "_NoSmear.root";

    if (TString(channel).EqualTo("ee")) cuts::bkg_baseline += cuts::ee;
    else if (TString(channel).EqualTo("mm")) cuts::bkg_baseline += cuts::mm;
    else if (TString(channel).EqualTo("em")) cuts::bkg_baseline += cuts::em;
    else {
        cout << "Unrecognized channel! quitting   " << channel << endl;
        exit(0);
    }

    //--- add files to TChain
    TChain* tch_zmc = new TChain("BaselineTree"); tch_zmc->Add(zmc_filename.c_str());
    TChain* tch_photon = new TChain("BaselineTree"); tch_photon->Add(photon_filename.c_str());

    //--- set plot info
    string plot_title;
    if (plot_feature == "met_Et" || plot_feature == "MET") plot_title = "E_{T}^{miss} [GeV]";
    else if (plot_feature == "METl") plot_title = "E_{T,||}^{miss} [GeV]";
    else if (plot_feature == "METt") plot_title = "E_{T,#perp}^{miss} [GeV]";

    int nbins = 50;
    float xmin = -300;
    float xmax = 300;

    //--- get standard deviations in bins
    int n_feature_bins = 20;
    float feature_bins[20+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200};
    vector<float> z_std, photon_std;

    for (int i = 0; i < n_feature_bins; i++) {
        string additional_cut = bin_feature + ">" + to_string(int(feature_bins[i])) + "&&" + bin_feature + "<" + to_string(int(feature_bins[i+1]));
        TCut extended_bkg_baseline = cuts::bkg_baseline + TCut(TString(additional_cut));
        TCut extended_photon_baseline = cuts::photon_baseline + TCut(TString(additional_cut));

        TH1F *h_zmc, *h_photon;
        h_zmc = new TH1F("h_zmc", "", nbins, xmin, xmax);
        h_photon = new TH1F("h_photon", "", nbins, xmin, xmax);
        h_zmc->SetLineColor(1); h_photon->SetLineColor(2);

        tch_zmc->Draw(Form("%s>>h_zmc", plot_feature.c_str()), extended_bkg_baseline*cuts::bkg_weight, "goff");
        tch_photon->Draw(Form("%s>>h_photon", plot_feature.c_str()), extended_photon_baseline*cuts::photon_weight, "goff");

        TCanvas *can = new TCanvas("can","can",600,600);
        can->cd();
        TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,1.0);
        mainpad->Draw();
        mainpad->cd();

        float max_y = max(h_zmc->GetMaximum(), h_photon->GetMaximum()) * 1.1;
        h_zmc->GetYaxis()->SetRange(0, max_y);
        h_photon->GetYaxis()->SetRange(0, max_y);

        h_zmc->SetTitle(TString(plot_feature + " " + additional_cut));
        h_zmc->Draw("hist");
        h_photon->Draw("samehist");
        can->Print(Form("%s_%s_%s_%s.eps", plot_feature.c_str(), mc_period.c_str(), channel.c_str(), additional_cut.c_str()));

        z_std.push_back(h_zmc->GetStdDev());
        photon_std.push_back(h_photon->GetStdDev());

        delete can;
        delete h_zmc;
        delete h_photon;
    }

    TH1F *h_z_std, *h_photon_std;
    h_z_std = new TH1F("h_z_std", "", n_feature_bins, feature_bins);
    h_photon_std = new TH1F("h_photon_std", "", n_feature_bins, feature_bins);
    h_z_std->SetLineColor(1); h_photon_std->SetLineColor(2);

    for (int i = 0; i < n_feature_bins; i++) {
        h_z_std->SetBinContent(i, z_std[i]);
        h_photon_std->SetBinContent(i, photon_std[i]);
    }

    //--- make plots
    TCanvas *can = new TCanvas("can","can",600,600);
    can->cd();
    TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);
    mainpad->Draw();
    mainpad->cd();

    float max_y = max(h_z_std->GetMaximum(), h_photon_std->GetMaximum()) * 1.1;
    h_z_std->GetYaxis()->SetRange(0, max_y);
    h_photon_std->GetYaxis()->SetRange(0, max_y);

    h_z_std->SetTitle(TString(plot_feature + " Standard Deviation"));
    h_z_std->Draw("hist");
    h_photon_std->Draw("samehist");

    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
    leg->AddEntry(h_z_std, "Z+jets (from MC)", "f");
    leg->AddEntry(h_photon_std, "Z+jets (from #gamma+jets, raw)", "f");

    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.03);
    tex->DrawLatex(0.6,0.65,"ATLAS Internal");
    if(TString(mc_period).Contains("mc16a")) tex->DrawLatex(0.6,0.61,"MC16a");
    if(TString(mc_period).Contains("mc16cd")) tex->DrawLatex(0.6,0.61,"MC16cd");
    if(TString(mc_period).Contains("mc16e")) tex->DrawLatex(0.6,0.61,"MC16cd");
    if(TString(channel).Contains("ee")) tex->DrawLatex(0.6,0.57,"ee events");
    if(TString(channel).Contains("em")) tex->DrawLatex(0.6,0.57,"e#mu events");
    if(TString(channel).Contains("mm")) tex->DrawLatex(0.6,0.57,"#mu#mu events");

    can->cd();
    TPad* ratio_pad = new TPad("ratio_pad","ratio_pad",0.0,0.8,1.0,1.0);
    ratio_pad->Draw();
    ratio_pad->cd();
    ratio_pad->SetGridy();

    TH1F* hratio;
    hratio = (TH1F*) h_z_std->Clone("hratio");
    hratio->SetTitle("");
    hratio->Divide(h_photon_std);
    hratio->SetMarkerStyle(20);

    hratio->GetXaxis()->SetTitle("");
    hratio->GetXaxis()->SetLabelSize(0.);
    hratio->GetYaxis()->SetNdivisions(5);
    hratio->GetYaxis()->SetTitle("Z/#gamma");
    hratio->GetYaxis()->SetTitleSize(0.15);
    hratio->GetYaxis()->SetTitleOffset(0.3);
    hratio->GetYaxis()->SetLabelSize(0.15);
    hratio->SetMinimum(0.0);
    hratio->SetMaximum(2.0);
    hratio->GetYaxis()->SetRangeUser(0.0,2.0);
    hratio->Draw("E1");

   can->Print(Form("%s_%s_%s_ratio.eps", plot_feature.c_str(), mc_period.c_str(), channel.c_str()));
}

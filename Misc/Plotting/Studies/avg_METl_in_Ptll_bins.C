#include "../../Common/Settings.C"

using namespace std;

void avg_METl_in_Ptll_bins(string mc_period, string channel) {

    gStyle->SetOptStat(0);

    //--- load files
    string zmc_filename = ntuple_path + "bkg_mc/" + mc_period + "_Zjets.root";
    string photon_filename = reweighting_path + "g_mc/" + mc_period + "_SinglePhoton222_" + channel + ".root";

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
    string plot_title = "E_{T,||}^{miss} in Bins of p_T";

    int n_metl_bins = 50;
    float metl_min = -300;
    float metl_max = 300;

    //--- get standard deviations in bins
    vector<float> z_mean, photon_mean, z_std, photon_std;

    int n_pt_bins = 15;
    double pt_bins[] = {0,30,45,60,80,100,120,140,160,180,200,220,260,280,300,350,400};

    for (int i=0; i<n_pt_bins; i++) {
        string additional_cut = "Ptll>" + to_string(int(pt_bins[i])) + "&&Ptll<" + to_string(int(pt_bins[i+1]));
        TCut extended_bkg_baseline = cuts::bkg_baseline + TCut(TString(additional_cut));
        TCut extended_photon_baseline = cuts::photon_baseline + TCut(TString(additional_cut));

        TH1F *h_zmc, *h_photon;
        h_zmc = new TH1F("h_zmc", "", n_metl_bins, metl_min, metl_max);
        h_photon = new TH1F("h_photon", "", n_metl_bins, metl_min, metl_max);
        h_zmc->SetLineColor(1); h_photon->SetLineColor(2);

        tch_zmc->Draw("METl>>h_zmc", extended_bkg_baseline*cuts::bkg_weight, "goff");
        tch_photon->Draw("METl>>h_photon", extended_photon_baseline*cuts::photon_weight, "goff");

        float SF = h_zmc->Integral() / h_photon->Integral();
        h_photon->Scale(SF);

        TCanvas *can = new TCanvas("can","can",600,600);
        can->cd();
        TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,1.0);
        mainpad->Draw();
        mainpad->cd();

        h_zmc->SetTitle(TString("METl " + additional_cut));
        h_zmc->Draw("hist");
        h_photon->Draw("samehist");
        can->Print(Form("Plots/METl_%s_%s_%s.eps", mc_period.c_str(), channel.c_str(), additional_cut.c_str()));

        z_mean.push_back(h_zmc->GetMean());
        photon_mean.push_back(h_photon->GetMean());
        z_std.push_back(h_zmc->GetStdDev());
        photon_std.push_back(h_photon->GetStdDev());

        delete can;
        delete h_zmc;
        delete h_photon;
    }

    TH1F *h_z_mean = new TH1F("h_z_mean", "", n_pt_bins, pt_bins);
    TH1F *h_photon_mean = new TH1F("h_photon_mean", "", n_pt_bins, pt_bins);
    h_z_mean->SetLineColor(1); h_photon_mean->SetLineColor(2);

    for (int i = 0; i < n_pt_bins; i++) {
        h_z_mean->SetBinContent(i, z_mean[i]);
        h_photon_mean->SetBinContent(i, photon_mean[i]);
        //h_z_std->SetBinContent(i, z_std[i]);
        //h_photon_std->SetBinContent(i, photon_std[i]);
    }

    //--- make plots
    TCanvas *can = new TCanvas("can","can",600,600);
    can->cd();
    TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);
    mainpad->Draw();
    mainpad->cd();

    float max_y = max(h_z_mean->GetMaximum(), h_photon_mean->GetMaximum()) * 1.1;
    h_z_mean->GetYaxis()->SetRange(0, max_y);
    h_photon_mean->GetYaxis()->SetRange(0, max_y);

    h_z_mean->SetTitle("METl Mean");
    h_z_mean->Draw("hist");
    h_photon_mean->Draw("samehist");

    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
    leg->AddEntry(h_z_mean, "Z+jets (from MC)", "f");
    leg->AddEntry(h_photon_mean, "Z+jets (from #gamma+jets, raw)", "f");

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
    hratio = (TH1F*) h_z_mean->Clone("hratio");
    hratio->SetTitle("");
    hratio->Divide(h_photon_mean);
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

   can->Print(Form("Plots/avg_METl_%s_%s_ratio.eps", mc_period.c_str(), channel.c_str()));
}

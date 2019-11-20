#include "../Common/Settings.C"

using namespace std;

// e.g. root -l -b -q 'plot_Z_lepton_cm_theta_distribution_vs_functions.C("mc16a", "ee", "baseline")'
void plot_Z_lepton_cm_theta_distribution_vs_functions(string mc_period, string channel, string selection) {

    gStyle->SetOptStat(0);

    //--- load file and add to TChain
    string zmc_filename = "/eos/user/m/mazhang/PhotonMethod/v1.6/LeptonDistributions/Ntuples/bkg_mc/" + mc_period + "_Zjets.root";
    TChain* tch_zmc = new TChain("BaselineTree"); tch_zmc->Add(zmc_filename.c_str());

    //--- set up event selections
    TCut Zselection;
    if (selection == "inclusive") Zselection = TCut("1");
    else if (selection == "baseline") Zselection = cuts::bkg_baseline;
    else Zselection = cuts::bkg_baseline + (TCut)(TString)selection;
    if (TString(channel).EqualTo("ee")) Zselection += cuts::ee;
    else if (TString(channel).EqualTo("mm")) Zselection += cuts::mm;
    else {
        cout << "Unrecognized channel! quitting   " << channel << endl;
        exit(0);
    }

    //--- set plot info
    string x_label = "#theta_{#ell} [GeV]";
    float xmin = 0; float xmax = 3.2; int nbins = 50;

    string plot_title = mc_period + " " + channel + " lepTheta";
    string save_title = "Plots/lepTheta_" + mc_period + "_" + channel + "_" + selection + "_distribution.eps";

    //--- create Z histogram
    TH1F *h_zmc;
    h_zmc = new TH1F("h_zmc", "", nbins, xmin, xmax);
    h_zmc->SetLineColor(1);
    tch_zmc->Draw(Form("%s>>h_zmc", "Z_cm_lep_theta"), Zselection*cuts::bkg_weight, "goff");

    TCanvas *can = new TCanvas("can","can",600,600);
    can->cd();
    h_zmc->SetTitle(plot_title.c_str());
    h_zmc->Draw("hist");

    //--- plot curves on top
    //TString scale = "25000";
    TString scale = "18600";
    TF1 *f1 = new TF1("f1", scale+"*pow(sin(x), 2)", 0, 3.2);
    //h_zmc->Fit(f1);
    f1->SetTitle("sin^{2}(#theta)");
    f1->SetLineWidth(3);
    f1->SetFillStyle(3005);
    f1->SetLineColor(3);
    f1->Draw("same");
    TF1 *f2 = new TF1("f2", scale+"*pow(sin(x), 3)", 0, 3.2);
    //h_zmc->Fit(f2);
    f2->SetTitle("sin^{3}(#theta)");
    f2->SetLineWidth(3);
    f2->SetFillStyle(3005);
    f2->SetLineColor(4);
    f2->Draw("same");
    TF1 *f3 = new TF1("f3", scale+"*pow(sin(x), 4)", 0, 3.2);
    //h_zmc->Fit(f3);
    f3->SetTitle("sin^{4}(#theta)");
    f3->SetLineWidth(3);
    f3->SetFillStyle(3005);
    f3->SetLineColor(5);
    f3->Draw("same");
    TF1 *f4 = new TF1("f4", scale+"*(1+pow(cos(x), 2))*sin(x)", 0, 3.2);
    //h_zmc->Fit(f4);
    f4->SetTitle("(1+cos^{2}(#theta))*sin(#theta)");
    f4->SetLineWidth(3);
    f4->SetFillStyle(3005);
    f4->SetLineColor(6);
    f4->Draw("same");
    gPad->BuildLegend();

    //--- save plot
    can->Print(save_title.c_str());
}

#include "../Common/CommonLibraries.C"
#include "../Common/Settings.C"

using namespace std;

// e.g. root -l -b -q 'Z_lepton_theta_distribution_plots.C("mc16a", "ee", "baseline")'
void Z_lepton_theta_distribution_plots(string mc_period, string channel, string selection) {

    gStyle->SetOptStat(0);

    //--- load file and add to TChain
    string zmc_filename = "/eos/user/m/mazhang/PhotonMethod/v1.6/LeptonDistributions/Ntuples/bkg_mc/" + mc_period + "_Zjets.root";
    TChain* tch_zmc = new TChain("BaselineTree"); tch_zmc->Add(zmc_filename.c_str());

    //--- set up event selections
    TCut Zselection;
    if (selection == "inclusive") Zselection = TCut("1");
    else if (selection == "baseline") Zselection = cuts::bkg_baseline;
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
    tch_zmc->Draw(Form("%s>>h_zmc", "2*TMath::ATan(exp(-lep_eta))"), Zselection*cuts::bkg_weight, "goff");

    TCanvas *can = new TCanvas("can","can",600,600);
    can->cd();
    h_zmc->SetTitle(plot_title.c_str());
    h_zmc->Draw("hist");

    //--- plot curves on top
    TF1 *f1 = new TF1("f1", "1 ++ pow(cos(x),2)", 0, 3.2);
    h_zmc->Fit(f1);
    f1->SetTitle("1+cos^{2}(#theta)");
    f1->SetLineWidth(3);
    f1->SetFillStyle(3005);
    f1->SetLineColor(2);
    f1->Draw("same");
    TF1 *f2 = new TF1("f2", "1 ++ pow(cos(x),4)", 0, 3.2);
    h_zmc->Fit(f2);
    f2->SetTitle("1+cos^{4}(#theta)");
    f2->SetLineWidth(3);
    f2->SetFillStyle(3005);
    f2->SetLineColor(3);
    f2->Draw("same");
    gPad->BuildLegend();

    //--- save plot
    can->Print(save_title.c_str());
}

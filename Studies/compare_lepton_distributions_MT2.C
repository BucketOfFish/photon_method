#include "../Common/CommonLibraries.C"
#include "../Common/Settings.C"

using namespace std;

void compare_lepton_distributions_MT2() {

    gStyle->SetOptStat(0);

    //--- load files
    TChain* tch_uniform = new TChain("BaselineTree");
    tch_uniform->Add("/eos/user/m/mazhang/PhotonMethod/v1.6/LeptonDistributions/SmearedNtuples/UniformSampling/g_mc/mc16a_SinglePhoton222_ee_NoSmear.root");
    TChain* tch_nonuniform = new TChain("BaselineTree");
    tch_nonuniform->Add("/eos/user/m/mazhang/PhotonMethod/v1.6/LeptonDistributions/SmearedNtuples/NonuniformSampling/g_mc/mc16a_SinglePhoton222_ee_NoSmear.root");
    TChain* tch_DY = new TChain("BaselineTree");
    tch_DY->Add("/eos/user/m/mazhang/PhotonMethod/v1.6/LeptonDistributions/SmearedNtuples/DrellYanSampling_CorrectedBoostAngle/g_mc/mc16a_SinglePhoton222_ee_NoSmear.root");

    //--- set up event selections
    TCut gselection = cuts::photon_baseline + cuts::ee;

    //--- set plot info
    string x_label = "MT2";
    float xmin = 0; float xmax = 1000; int nbins = 50;

    string plot_title = "MC16a ee MT2";
    string save_title = "Plots/lepton_MT2_comparison.eps";

    //--- plot
    TCanvas *can = new TCanvas("can","can",600,600);
    can->cd();
    tch_uniform->SetTitle(plot_title.c_str());
    tch_uniform->Draw("hist");
    tch_nonuniform->Draw("hist");
    tch_DY->Draw("hist");
    gPad->BuildLegend();

    //--- save plot
    can->Print(save_title.c_str());
}

#include "../Common/Settings.C"
#include "../Common/CommonLibraries.C"
#include "../Common/CommonFunctions.C"

using namespace std;

//--- run with:
//--- root -l -b -q 'compare_feature_for_Z_vs_photon.C("mc16a", "ee", "inclusive", "non-uniform", "lepPt")'
//--- 1. mc16a, mc16cd
//--- 2. ee, mm
//--- 3. inclusive, baseline
//--- 4. non-uniform, uniform, Drell-Yan, sin3
//--- 5. lepPt, lepEta, lepPhi, MT2, mll, Ptll

//--- Compare photon vs. Z feature distributions for different types of angular distributions

void compare_feature_for_Z_vs_photon(string mc_period, string channel, string selection, string distribution, string feature) {

    gStyle->SetOptStat(0);

    //--- load files
    string zmc_filename = "/eos/user/m/mazhang/PhotonMethod/v1.7/Default/Ntuples/bkg_mc/" + mc_period + "_Zjets.root";
    string distribution_folder;
    if (distribution == "non-uniform") distribution_folder = "NonuniformSampling";
    else if (distribution == "uniform") distribution_folder = "UniformSampling";
    else if (distribution == "Drell-Yan") distribution_folder = "DrellYanSampling_CorrectedBoostAngle";
    else if (distribution == "sin3") distribution_folder = "Sin3Sampling";
    string photon_filename = "/eos/user/m/mazhang/PhotonMethod/v1.7/Default/" + distribution_folder + "/ReweightedNtuples/g_mc/" + mc_period + "_SinglePhoton222_" + channel + "_NoSmear.root";

    //--- add files to TChain
    TChain* tch_zmc = new TChain("BaselineTree"); tch_zmc->Add(zmc_filename.c_str());
    TChain* tch_photon = new TChain("BaselineTree"); tch_photon->Add(photon_filename.c_str());

    //--- set up event selections
    TCut Zselection, gselection;
    if (selection == "inclusive") { Zselection = TCut("1"); gselection = TCut("1"); }
    else if (selection == "baseline") { Zselection = cuts::bkg_baseline; gselection = cuts::photon_baseline; }
    if (TString(channel).EqualTo("ee")) { Zselection += cuts::ee; gselection += cuts::ee; }
    else if (TString(channel).EqualTo("mm")) { Zselection += cuts::mm; gselection += cuts::mm; }
    else {
        cout << "Unrecognized channel! quitting   " << channel << endl;
        exit(0);
    }

    //--- set plot info
    string x_label;
    int nbins = 50;
    float xmin = 0;
    float xmax = 200;

    if (feature == "lepPt") { x_label = "p_{T,#ell} [GeV]"; }
    else if (feature == "lepEta") { x_label = "#eta_{#ell} [GeV]"; xmin = -3; xmax = 3; }
    else if (feature == "lepPhi") { x_label = "#phi_{#ell} [GeV]"; xmin = -TMath::Pi(); xmax = TMath::Pi(); }
    else if (feature == "mll") { x_label = "m_{#ell#ell} [GeV]"; xmin = 0; xmax = 200; }
    else if (feature == "MT2") { x_label = "MT2 [GeV]"; xmin = -500; xmax = 500; }
    else if (feature == "Ptll") { x_label = "Pt_{#ell#ell} [GeV]"; xmin = 0; xmax = 300; }

    string plot_title = mc_period + " " + channel + " " + feature;
    string save_title = "Plots/" + feature + "_" + mc_period + "_" + channel + "_" + selection + "_" + distribution + "_comparison.eps";

    //--- create photon and Z histograms
    TH1F *h_zmc, *h_photon, *h_photon_rw;
    h_zmc = new TH1F("h_zmc", "", nbins, xmin, xmax);
    h_photon = new TH1F("h_photon", "", nbins, xmin, xmax);
    h_photon_rw = new TH1F("h_photon_rw", "", nbins, xmin, xmax);
    h_zmc->SetLineColor(1); h_photon->SetLineColor(2); h_photon_rw->SetLineColor(3);

    string feature_string = feature;
    if (feature == "lepEta") feature_string = "lep_eta";
    if (feature == "lepPhi") feature_string = "lep_phi";

    tch_zmc->Draw(Form("%s>>h_zmc", feature_string.c_str()), Zselection*cuts::bkg_weight, "goff");
    //if (feature == "Ptll")
        //tch_photon->Draw(Form("%s>>h_photon", feature_string.c_str()), gselection*cuts::photon_weight, "goff");
    tch_photon->Draw(Form("%s>>h_photon_rw", feature_string.c_str()), gselection*cuts::photon_weight_rw, "goff");

    float ymax = max(max(h_zmc->GetMaximum(), h_photon->GetMaximum()), h_photon_rw->GetMaximum()) * 1.1;
    h_zmc->GetYaxis()->SetRangeUser(0, ymax);
    h_photon->GetYaxis()->SetRangeUser(0, ymax);
    h_photon_rw->GetYaxis()->SetRangeUser(0, ymax);

    h_zmc->SetTitle(plot_title.c_str());

    //--- draw histograms on top of each other
    TCanvas *can = new TCanvas("can","can",600,600);
    can->cd();
    TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);
    mainpad->Draw();
    mainpad->cd();

    h_zmc->Draw("hist");
    //if (feature == "Ptll")
        //h_photon->Draw("samehist");
    h_photon_rw->Draw("samehist");

    //--- draw legend
    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
    leg->AddEntry(h_zmc, "Z+jets (from MC)", "f");
    //if (feature == "Ptll")
        //leg->AddEntry(h_photon, "Z+jets (from #gamma+jets)", "f");
    leg->AddEntry(h_photon_rw, "Z+jets (from reweighted #gamma+jets)", "f");

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

    //--- draw ratio of distributions
    can->cd();
    TPad* ratio_pad = new TPad("ratio_pad","ratio_pad",0.0,0.8,1.0,1.0);
    ratio_pad->Draw();
    ratio_pad->cd();
    ratio_pad->SetGridy();

    TH1F* hratio;
    hratio = (TH1F*) h_zmc->Clone("hratio");
    hratio->SetTitle("");
    hratio->Divide(h_photon_rw);
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

    //--- save plot
    can->Print(save_title.c_str());
}

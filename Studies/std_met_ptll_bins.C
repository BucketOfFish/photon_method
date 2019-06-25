#include "../Common/Settings.C"
#include "../Common/CommonLibraries.C"
#include "../Common/CommonFunctions.C"

using namespace std;

//--- run with:
//--- root -l -b -q 'std_met_ptll_bins.C("mc16a", "ee", "met_Et", "Ptll")'

void std_met_ptll_bins(string mc_period, string channel, string plot_feature, string bin_feature) {

    gStyle->SetOptStat(0);

    //--- load files
    string zmc_filename = ntuple_path + "bkg_mc/" + mc_period + "_Zjets.root";
    string photon_filename = reweighting_path + "g_mc/" + mc_period + "_SinglePhoton222_" + channel + "_NoSmear.root";

    if (TString(channel).EqualTo("ee")) cuts::Zselection += cuts::ee;
    else if (TString(channel).EqualTo("mm")) cuts::Zselection += cuts::mm;
    else if (TString(channel).EqualTo("em")) cuts::Zselection += cuts::em;
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
    std::vector<double> feature_bins = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 1000};

    for (int i = 0; i < feature_bins.size()-2; i++) {
        TString additional_cut = TString(bin_feature + ">" + feature_bins[i] + "&&" + bin_feature + "<" + feature_bins[i+1]);
        TCut extended_Zselection = cuts::Zselection + TCut(additional_cut);
        TCut extended_gselection = cuts::gselection + TCut(additional_cut);

        TH1F *h_zmc, *h_photon;
        h_zmc = new TH1F("h_zmc", "", nbins, xmin, xmax);
        h_photon = new TH1F("h_photon", "", nbins, xmin, xmax);

        tch_zmc->Draw(Form("%s>>h_zmc", plot_feature.c_str()), extended_Zselection*cuts::Zweight, "goff");
        tch_photon->Draw(Form("%s>>h_photon", plot_feature.c_str()), extended_gselection*cuts::weight_g, "goff");

        cout << h_zmc->GetStdDev() << " | " << h_photon->GetStdDev() << endl;
        delete h_zmc;
        delete h_photon;
    }

    ////--- make plots
    //string plotName = "default cuts";
    //if (additional_cut != "1") plotName += ("&&" + additional_cut);
    //boost::replace_all(plotName, "&&", " && ");
    //THStack *mcstack = new THStack("mcstack", plotName.c_str());

    //h_tt->SetLineColor(1); h_tt->SetFillColor(kRed-2);
    //h_vv->SetLineColor(1); h_vv->SetFillColor(kGreen-2);
    //h_photon_reweighted->SetLineColor(1); h_photon_reweighted->SetFillColor(kOrange-2);
    //mcstack->Add(h_tt);
    //mcstack->Add(h_vv);
    //mcstack->Add(h_photon_reweighted);

    //h_photon->Add(h_tt); h_photon->Add(h_vv);
    //h_photon->SetLineColor(4); h_photon->SetLineWidth(1); h_photon->SetLineStyle(2);

    //h_zmc->Add(h_tt); h_zmc->Add(h_vv);
    //h_zmc->SetLineColor(2); h_zmc->SetLineWidth(1); h_zmc->SetLineStyle(7);

    //TCanvas *can = new TCanvas("can","can",600,600);
    //can->cd();
    //TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);
    //mainpad->Draw();
    //mainpad->cd();
    //mainpad->SetLogy();

    //if (photon_data_or_mc == "Data") {
        //mcstack->Draw("hist");
        //mcstack->GetXaxis()->SetTitle(xtitle.c_str());
        //mcstack->GetYaxis()->SetTitle("entries / bin");

        //h_photon->Draw("samehist");
        //h_zmc->Draw("samehist");
        //h_zdata->SetLineColor(1); h_zdata->SetLineWidth(2); h_zdata->SetMarkerStyle(20);
        //h_zdata->Draw("sameE1");
    //}
    //else {
        //h_zmc->SetLineColor(1); h_zmc->SetFillColor(42); h_zmc->SetLineStyle(1);
        //h_zmc->GetXaxis()->SetTitle(xtitle.c_str());
        //h_zmc->GetYaxis()->SetTitle("entries / bin");
        //h_zmc->SetTitle(plotName.c_str());
        //h_zmc->Draw("hist");
        //h_photon->Draw("samehist");
        //h_photon_reweighted->SetLineWidth(1); h_photon_reweighted->SetFillStyle(0);
        //h_photon_reweighted->Draw("samehist");
    //}

    ////--- draw legend and labels
    //TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
    //if (photon_data_or_mc == "Data") {
        //leg->AddEntry(h_zdata,"data","lp");
        //leg->AddEntry(h_photon, "Z+jets (from #gamma+jets, raw)", "f");
        //leg->AddEntry(h_zmc, "Z+jets (from MC)", "f");
        //leg->AddEntry(h_photon_reweighted, "Z+jets (from #gamma+jets, reweighted)", "f");
        //leg->AddEntry(h_vv, "VV", "f");
        //leg->AddEntry(h_tt, "t#bar{t}+tW", "f");
    //}
    //else {
        //leg->AddEntry(h_photon, "Z+jets (from #gamma+jets, raw)", "f");
        //leg->AddEntry(h_zmc, "Z+jets (from MC)", "f");
        //leg->AddEntry(h_photon_reweighted, "Z+jets (from #gamma+jets, reweighted)", "f");
    //}

    //leg->SetBorderSize(0);
    //leg->SetFillColor(0);
    //leg->Draw();

    //TLatex *tex = new TLatex();
    //tex->SetNDC();
    //tex->SetTextSize(0.03);
    //tex->DrawLatex(0.6,0.65,"ATLAS Internal");
    //if (photon_data_or_mc == "Data") {
        //if(TString(period).Contains("data15-16")) tex->DrawLatex(0.6,0.61,"36 fb^{-1} 2015-2016 data");
        //if(TString(period).Contains("data17")) tex->DrawLatex(0.6,0.61,"44 fb^{-1} 2017 data");
    //}
    //else {
        //if(TString(period).Contains("data15-16")) tex->DrawLatex(0.6,0.61,"MC16a");
        //if(TString(period).Contains("data17")) tex->DrawLatex(0.6,0.61,"MC16cd");
        //if(TString(period).Contains("data18")) tex->DrawLatex(0.6,0.61,"MC16cd");
    //}
    //if(TString(channel).Contains("ee")) tex->DrawLatex(0.6,0.57,"ee events");
    //if(TString(channel).Contains("em")) tex->DrawLatex(0.6,0.57,"e#mu events");
    //if(TString(channel).Contains("mm")) tex->DrawLatex(0.6,0.57,"#mu#mu events");

    ////--- draw ratio
    //can->cd();
    //TPad* ratio_pad = new TPad("ratio_pad","ratio_pad",0.0,0.8,1.0,1.0);
    //ratio_pad->Draw();
    //ratio_pad->cd();
    //ratio_pad->SetGridy();

    //TH1F* hratio;
    //hratio = (TH1F*) h_zdata->Clone("hratio");
    //TH1F* hmctot = (TH1F*) h_photon_reweighted->Clone("hmctot");
    //hmctot->Add(h_tt);
    //hmctot->Add(h_vv);

    //if (photon_data_or_mc == "MC") {
        //hratio = (TH1F*) h_zmc->Clone("hratio");
        //hratio->SetTitle("");
        //hmctot = (TH1F*) h_photon_reweighted->Clone("hmctot");
    //}

    //for (int ibin=1; ibin <= hmctot->GetXaxis()->GetNbins(); ibin++)
        //hmctot->SetBinError(ibin, 0.0);
    //hratio->Divide(hmctot);
    //hratio->SetMarkerStyle(20);

    //hratio->GetXaxis()->SetTitle("");
    //hratio->GetXaxis()->SetLabelSize(0.);
    //hratio->GetYaxis()->SetNdivisions(5);
    //if (photon_data_or_mc == "Data")
        //hratio->GetYaxis()->SetTitle("data/bkg");
    //else
        //hratio->GetYaxis()->SetTitle("Z/#gamma MC");
    //hratio->GetYaxis()->SetTitleSize(0.15);
    //hratio->GetYaxis()->SetTitleOffset(0.3);
    //hratio->GetYaxis()->SetLabelSize(0.15);
    //hratio->SetMinimum(0.0);
    //hratio->SetMaximum(2.0);
    //hratio->GetYaxis()->SetRangeUser(0.0,2.0);
    //hratio->Draw("E1");

    ////--- save plot
    //if (additional_cut == "1")
        //can->Print(Form("%s/%s_%s_%s_%s_%s.eps", plots_path.c_str(), period.c_str(), channel.c_str(), smearing_mode.c_str(), plot_feature.c_str(), ("photon-"+photon_data_or_mc).c_str()));
    //else
        //can->Print(Form("%s/%s_%s_%s_%s_%s_%s.eps", plots_path.c_str(), period.c_str(), channel.c_str(), smearing_mode.c_str(), plot_feature.c_str(), additional_cut.c_str(), ("photon-"+photon_data_or_mc).c_str()));
}

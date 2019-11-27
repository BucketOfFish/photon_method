#include "../Common/Settings.C"
#include <boost/algorithm/string/replace.hpp>

using namespace std;

void quickDraw(string period="data15-16", string channel="mm" , string plot_feature="HT", string smearing_mode="NoSmear", string photon_data_or_mc="Data", string region="SR", string additional_cut="1") {

    bool DF = TString(channel).EqualTo("em");
    gStyle->SetOptStat(0);

    cout << "period               " << period << endl;
    cout << "channel              " << channel << endl;
    cout << "smearing mode        " << smearing_mode << endl;
    cout << "DF?                  " << DF << endl;
    cout << "photon data          " << photon_data_or_mc << endl;

    //--- load files
    string mc_period = "";
    if (TString(period).Contains("data15-16")) mc_period = "mc16a";
    else if (TString(period).Contains("data17")) mc_period = "mc16cd";
    else if (TString(period).Contains("data18")) mc_period = "mc16e";

    string data_filename= ntuple_path + "bkg_data/" + period + "_bkg.root";
    string tt_filename = ntuple_path + "bkg_mc/" + mc_period + "_ttbar.root";
    string vv_filename = ntuple_path + "bkg_mc/" + mc_period + "_diboson.root";
    string zmc_filename = ntuple_path + "bkg_mc/" + mc_period + "_Zjets.root";
    string photon_filename;
    if (photon_data_or_mc == "MC") photon_filename = reweighting_path + "g_mc/" + mc_period + "_SinglePhoton222_" + channel + "_" + smearing_mode + ".root";
    else photon_filename = reweighting_path + "g_data/" + period + "_photon_" + channel + "_" + smearing_mode + ".root";

    cout << "data filename        " << data_filename << endl;
    cout << "ttbar filename       " << tt_filename << endl;
    cout << "diboson filename     " << vv_filename << endl;
    cout << "Z MC filename        " << zmc_filename << endl;
    cout << "photon filename      " << photon_filename << endl;
    cout << "" << endl;

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

    if (plot_feature == "met_Et") plot_settings = std::make_tuple("E_{T}^{miss} [GeV]", 20, 0, 200);
    else if (plot_feature == "MET") plot_settings = std::make_tuple("E_{T}^{miss} [GeV]", 20, 0, 200);
    else if (plot_feature == "METl") plot_settings = std::make_tuple("E_{T}^{miss}{||} [GeV]", 30, -150, 150);
    else if (plot_feature == "METt") plot_settings = std::make_tuple("E_{T}^{miss}{#perp} [GeV]", 30, -150, 150);
    else if (plot_feature == "MET_loose") plot_settings = std::make_tuple("E_{T}^{miss}{loose} [GeV]", 20, 0, 200);
    else if (plot_feature == "MET_tight") plot_settings = std::make_tuple("E_{T}^{miss}{tight} [GeV]", 20, 0, 200);
    else if (plot_feature == "MET_tighter") plot_settings = std::make_tuple("E_{T}^{miss}{tighter} [GeV]", 20, 0, 200);
    else if (plot_feature == "MET_tenacious") plot_settings = std::make_tuple("E_{T}^{miss}{tenacious} [GeV]", 20, 0, 200);
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

    //--- draw histograms
    tch_data->Draw(Form("%s>>h_data", plot_feature.c_str()), plot_region, "goff");
    tch_tt->Draw(Form("%s>>h_tt", plot_feature.c_str()), plot_region*cuts::bkg_weight, "goff");
    tch_vv->Draw(Form("%s>>h_vv", plot_feature.c_str()), plot_region*cuts::bkg_weight, "goff");
    if (!DF) {
        tch_zmc->Draw(Form("%s>>h_zmc", plot_feature.c_str()), plot_region*cuts::bkg_weight, "goff");
        tch_photon->Draw(Form("%s>>h_photon", plot_feature.c_str()), plot_region*cuts::photon_weight, "goff");
        tch_photon->Draw(Form("%s>>h_photon_reweighted", plot_feature.c_str()), plot_region*cuts::photon_weight_rw, "goff");
    }

    cout << "" << endl;
    cout << "data integral        " << h_data->Integral() << endl;
    cout << "tt integral          " << h_tt->Integral() << endl;
    cout << "VV integral          " << h_vv->Integral() << endl;
    cout << "Z MC integral        " << h_zmc->Integral() << endl;
    cout << "g raw integral       " << h_photon->Integral() << endl;
    cout << "g reweighted int.    " << h_photon_reweighted->Integral() << endl;
    cout << "" << endl;

    //--- normalize Z to MET<60 GeV region
    cout << "normalize to CR " << cuts::CR.GetTitle() << endl;

    TH1F* h_data_cr = new TH1F("h_data_cr", "", 1, 0, 1);
    TH1F* h_tt_cr = new TH1F("h_tt_cr", "", 1, 0, 1);
    TH1F* h_vv_cr = new TH1F("h_vv_cr", "", 1, 0, 1);
    TH1F* h_zmc_cr = new TH1F("h_zmc_cr", "", 1, 0, 1);
    TH1F* h_photon_cr = new TH1F("h_photon_cr", "", 1, 0, 1);
    TH1F* h_photon_reweighted_cr = new TH1F("h_photon_reweighted_cr", "", 1, 0, 1);

    TCut current_CR = plot_region + cuts::CR;
    tch_data->Draw("0.5>>h_data_cr", current_CR*cuts::bkg_weight, "goff");
    tch_tt-> Draw("0.5>>h_tt_cr", current_CR*cuts::bkg_weight, "goff");
    tch_vv-> Draw("0.5>>h_vv_cr", current_CR*cuts::bkg_weight, "goff");
    tch_zmc->Draw("0.5>>h_zmc_cr", current_CR*cuts::bkg_weight, "goff");
    tch_photon->Draw("0.5>>h_photon_cr", current_CR*cuts::photon_weight, "goff");
    tch_photon->Draw("0.5>>h_photon_reweighted_cr", current_CR*cuts::photon_weight_rw, "goff");

    float SF = (h_data_cr->Integral() - h_tt_cr->Integral() - h_vv_cr->Integral()) / h_photon_cr->Integral();
    float SFrw = (h_data_cr->Integral() - h_tt_cr->Integral() - h_vv_cr->Integral()) / h_photon_reweighted_cr->Integral();

    if (photon_data_or_mc == "MC") {
        SF = h_zmc_cr->Integral() / h_photon_cr->Integral();
        SFrw = h_zmc_cr->Integral() / h_photon_reweighted_cr->Integral();
    }

    cout << "Scale raw photon data by " << SF << endl;
    cout << "Scale reweighted photon data by " << SFrw << endl;

    h_photon->Scale(SF);
    h_photon_reweighted->Scale(SFrw);

    //--- print MET integrals
    cout << "MET100-150" << endl;
    cout << "2L data                " << h_data->Integral(11,15) << endl;
    cout << "VV MC                  " << h_vv->Integral(11,15) << endl;
    cout << "tt MC                  " << h_tt->Integral(11,15) << endl;
    cout << "Z+jets MC              " << h_zmc->Integral(11,15) << endl;
    cout << "g data (reweighted)    " << h_photon_reweighted->Integral(11,15) << endl;
    cout << "g data (raw)           " << h_photon->Integral(11,15) << endl;

    cout << "MET150-200" << endl;
    cout << "2L data                " << h_data->Integral(16,21) << endl;
    cout << "VV MC                  " << h_vv->Integral(16,21) << endl;
    cout << "tt MC                  " << h_tt->Integral(16,21) << endl;
    cout << "Z+jets MC              " << h_zmc->Integral(16,21) << endl;
    cout << "g data (reweighted)    " << h_photon_reweighted->Integral(16,21) << endl;
    cout << "g data (raw)           " << h_photon->Integral(16,21) << endl;

    //--- create MC stack
    THStack *mcstack = new THStack("mcstack", "");

    h_tt->SetLineColor(1); h_tt->SetFillColor(kRed-2);
    h_vv->SetLineColor(1); h_vv->SetFillColor(kGreen-2);
    h_photon_reweighted->SetLineColor(1); h_photon_reweighted->SetFillColor(kOrange-2);
    mcstack->Add(h_tt);
    mcstack->Add(h_vv);
    if(!DF) mcstack->Add(h_photon_reweighted);

    //--- create comparison "stacks"
    h_photon->Add(h_tt); h_photon->Add(h_vv);
    h_photon->SetLineColor(4); h_photon->SetLineWidth(1); h_photon->SetLineStyle(2);

    h_zmc->Add(h_tt); h_zmc->Add(h_vv);
    h_zmc->SetLineColor(2); h_zmc->SetLineWidth(1); h_zmc->SetLineStyle(7);

    //--- turn on overflow bin
    h_zmc->GetXaxis()->SetRange(0, h_zmc->GetNbinsX() + 1);
    h_photon->GetXaxis()->SetRange(0, h_photon->GetNbinsX() + 1);
    h_photon_reweighted->GetXaxis()->SetRange(0, h_photon_reweighted->GetNbinsX() + 1);

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
    TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);
    mainpad->Draw();
    mainpad->cd();
    mainpad->SetLogy();

    if (photon_data_or_mc == "Data") {
        mcstack->Draw("hist");
        mcstack->GetXaxis()->SetTitle(formatted_feature);
        mcstack->GetYaxis()->SetTitle("entries / bin");

        if( !DF ) {
            h_zmc->Draw("samehist");
            h_photon->Draw("samehist");
        }
        h_data->SetLineColor(1); h_data->SetLineWidth(2); h_data->SetMarkerStyle(20);
        h_data->Draw("sameE1");
    }
    else {
        h_zmc->SetLineColor(1); h_zmc->SetFillColor(42); h_zmc->SetLineStyle(1);
        h_zmc->GetXaxis()->SetTitle(formatted_feature);
        h_zmc->GetYaxis()->SetTitle("entries / bin");
        h_zmc->Draw("hist");
        h_photon->Draw("samehist");
        h_photon_reweighted->SetLineWidth(1); h_photon_reweighted->SetFillStyle(0);
        h_photon_reweighted->Draw("samehist");
    }

    //--- draw legend and labels
    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
    if (photon_data_or_mc == "Data") {
        leg->AddEntry(h_data,"data","lp");
        if(!DF){
            leg->AddEntry(h_zmc, "Z+jets (from MC)", "f");
            leg->AddEntry(h_photon, "Z+jets (from #gamma+jets, raw)", "f");
            leg->AddEntry(h_photon_reweighted, "Z+jets (from #gamma+jets, reweighted)", "f");
        }
        leg->AddEntry(h_vv, "VV", "f");
        leg->AddEntry(h_tt, "t#bar{t}+tW", "f");
    }
    else {
        leg->AddEntry(h_zmc, "Z+jets (from MC)", "f");
        leg->AddEntry(h_photon, "Z+jets (from #gamma+jets, raw)", "f");
        leg->AddEntry(h_photon_reweighted, "Z+jets (from #gamma+jets, reweighted)", "f");
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

    TH1F* hratio;
    hratio = (TH1F*) h_data->Clone("hratio");
    TH1F* hmctot = (TH1F*) h_photon_reweighted->Clone("hmctot");
    hmctot->Add(h_tt);
    hmctot->Add(h_vv);

    if (photon_data_or_mc == "MC") {
        hratio = (TH1F*) h_zmc->Clone("hratio");
        hmctot = (TH1F*) h_photon_reweighted->Clone("hmctot");
    }

    for (int ibin=1; ibin <= hmctot->GetXaxis()->GetNbins(); ibin++)
        hmctot->SetBinError(ibin, 0.0);
    hratio->Divide(hmctot);
    hratio->SetMarkerStyle(20);

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
    hratio->Draw("E1");

    //--- save plot
    TString plot_name;
    if (photon_data_or_mc == "Data")
        plot_name = Form("%s/%s_%s_%s_%s_%s", plots_path.c_str(), period.c_str(), channel.c_str(), smearing_mode.c_str(), plot_feature.c_str(), region.c_str());
    else
        plot_name = Form("%s/%s_%s_%s_%s_%s", plots_path.c_str(), mc_period.c_str(), channel.c_str(), smearing_mode.c_str(), plot_feature.c_str(), region.c_str());
    if (additional_cut != "1")
        plot_name += ("_" + additional_cut).c_str();
    plot_name += ".eps";
    can->Print(plot_name);
}

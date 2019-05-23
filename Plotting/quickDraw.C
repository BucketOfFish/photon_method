#include "../Common/Settings.C"
#include "../Common/CommonLibraries.C"
#include "../Common/CommonCuts.C"
#include "../Common/CommonFunctions.C"

using namespace std;

void quickDraw(string period="data15-16", string channel="mm" , string var="HT", string smearing_mode="NoSmear", string isData="Data", string region="SR") {

    bool DF = TString(channel).EqualTo("em");

    bool normalize = true; // whether or not to normalize Z to MET<60 GeV region
    if (DF || TString(var).Contains("pt") || TString(var).Contains("HT")) normalize = false;

    gStyle->SetOptStat(0);

    //--- load files
    string Zfilename, gfilename;
    if (isData == "MC") {
        string Zfilename = ntuple_path + "ZMC16a/Zjets_merged_processed.root";
        string gfilename = reweighting_path + "gmc/gmc_" + channel + "_" + smearing_mode + ".root";
    }
    else {
        string Zfilename = ntuple_path + "zdata/data15-16_merged_processed.root";
        string gfilename = reweighting_path + "gdata/" + period + "_merged_processed_" + channel + "_" + smearing_mode + ".root";
    }

    string mcdir = "";
    if (TString(period).Contains("data15-16")) mcdir = "ZMC16a/";
    else if (TString(period).Contains("data17")) mcdir = "ZMC16cd/";
    else if (TString(period).Contains("data18")) mcdir = "ZMC16cd/";
    string tt_filename = ntuple_path + mcdir + "ttbar_merged_processed.root";
    string vv_filename = ntuple_path + mcdir + "diboson_merged_processed.root";
    string z_filename = ntuple_path + mcdir + "Zjets_merged_processed.root";

    cout << "period               " << period        << endl;
    cout << "channel              " << channel       << endl;
    cout << "smearing mode        " << smearing_mode << endl;
    cout << "Z datafilename       " << Zfilename     << endl;
    if (isData == "MC") {
        cout << "tt filename          " << tt_filename   << endl;
        cout << "vv filename          " << vv_filename   << endl;
        cout << "Z+jets filename      " << z_filename    << endl;
    }
    cout << "g filename           " << gfilename     << endl;
    cout << "DF?                  " << DF            << endl;
    cout << "" << endl;

    //--- add files to TChain
    TChain* gtree = new TChain("BaselineTree"); if (!DF) gtree->Add(gfilename.c_str());
    TChain* Ztree = new TChain("BaselineTree"); Ztree->Add(Zfilename.c_str());
    TChain* chtt = new TChain("BaselineTree"); chtt->Add(tt_filename.c_str());
    TChain* chz = new TChain("BaselineTree"); chz->Add(z_filename.c_str());
    TChain* chvv = new TChain("BaselineTree"); chvv->Add(vv_filename.c_str());

    cout << "g entries            " << gtree->GetEntries()  << endl;
    cout << "Z data entries       " << Ztree->GetEntries()  << endl;
    if (isData == "MC") {
        cout << "ttbar entries        " << chtt->GetEntries()   << endl;
        cout << "diboson entries      " << chvv->GetEntries()   << endl;
        cout << "Z+jets entries       " << chz->GetEntries()    << endl;
    }

    //--- define selections
    if (TString(channel).EqualTo("ee")) cuts::Zselection += cuts::ee;
    else if (TString(channel).EqualTo("mm")) cuts::Zselection += cuts::mm;
    else if (TString(channel).EqualTo("em")) cuts::Zselection += cuts::em;
    else {
        cout << "Unrecognized channel! quitting   " << channel << endl;
        exit(0);
    }

    if (region == "VR") cout << "VR selection         " << cuts::VR.GetTitle()  << endl;
    cout << "Z selection          " << cuts::Zselection.GetTitle()  << endl;  
    cout << "Z weight             " << cuts::Zweight.GetTitle()     << endl;
    cout << "g selection          " << cuts::gselection.GetTitle()  << endl;
    cout << "g weight             " << cuts::weight_g.GetTitle()    << endl;
    cout << "g weight (reweight)  " << cuts::weight_g_rw.GetTitle() << endl;

    //--- set histogram binning
    std::tuple<string, int, float, float> plot_settings;

    if (var == "met_Et") plot_settings = std::make_tuple("E_{T}^{miss} [GeV]", 20, 0, 400);
    else if (var == "METl") plot_settings = std::make_tuple("E_{T,||}^{miss} [GeV]", 20, -200, 200);
    else if (var == "METt") plot_settings = std::make_tuple("E_{T,#perp}^{miss} [GeV]", 20, -200, 200);
    else if (var == "MET_loose") plot_settings = std::make_tuple("E_{T,loose}^{miss} [GeV]", 20, 0, 200);
    else if (var == "MET_tight") plot_settings = std::make_tuple("E_{T,tight}^{miss} [GeV]", 20, 0, 200);
    else if (var == "MET_tighter") plot_settings = std::make_tuple("E_{T,tighter}^{miss} [GeV]", 20, 0, 200);
    else if (var == "MET_tenacious") plot_settings = std::make_tuple("E_{T,tenacious}^{miss} [GeV]", 20, 0, 200);
    else if (var == "Z_pt") plot_settings = std::make_tuple("p_{T} [GeV]", 20, 0, 100);
    else if (var == "jet_n") plot_settings = std::make_tuple("n_{jets}", 6, 2, 8);
    else if (var == "bjet_n") plot_settings = std::make_tuple("n_{b-jets}", 4, 0, 4);
    else if (var == "HT") plot_settings = std::make_tuple("H_{T}", 20, 0, 1000);
    else if (var == "mll") plot_settings = std::make_tuple("m_{ll} [GeV]", 30, 0, 300);
    else if (var == "MT2W") plot_settings = std::make_tuple("m_{T2}^{W} [GeV]", 20, 0, 200);
    else if (var == "lep_pT[0]") plot_settings = std::make_tuple("1^{st} lepton p_{T} [GeV]", 20, 0, 200);
    else if (var == "lep_pT[1]") plot_settings = std::make_tuple("2^{nd} lepton p_{T} [GeV]", 20, 0, 100);
    else if (var == "DPhi_METJetLeading") plot_settings = std::make_tuple("#Delta#phi(jet_{1},E_{T}^{miss})", 20, 0, 3.14);
    else if (var == "DPhi_METJetSecond") plot_settings = std::make_tuple("#Delta#phi(jet_{2},E_{T}^{miss})", 20, 0, 3.14);
    else {
        cout << "Error! unrecognized variable, need to set binning, quitting! " << var << endl;
        exit(0);
    }

    string xtitle = std::get<0>(plot_settings);
    int nbins = std::get<1>(plot_settings);
    float xmin = std::get<2>(plot_settings);
    float xmax = std::get<3>(plot_settings);

    //--- initialize histograms
    TH1F *hZ, *hg, *hg_rw, *htt, *hvv, *hz;

    if ((var == "Z_pt") || (var == "HT")) {
        const unsigned int n_new_bins = 16;
        double new_bins[n_new_bins+1] = {40,75,100,125,150,175,200,250,300,350,400,450,500,600,700,850,1000};
        hZ = new TH1F("hZ", "", n_new_bins, new_bins);
        hg = new TH1F("hg", "", n_new_bins, new_bins);
        hg_rw = new TH1F("hg_rw", "", n_new_bins, new_bins);
        htt = new TH1F("htt", "", n_new_bins, new_bins);
        hvv = new TH1F("hvv", "", n_new_bins, new_bins);
        hz = new TH1F("hz", "", n_new_bins, new_bins);
    }
    else {
        hZ = new TH1F("hZ", "", nbins, xmin, xmax);
        hg = new TH1F("hg", "", nbins, xmin, xmax);
        hg_rw = new TH1F("hg_rw", "", nbins, xmin, xmax);
        htt = new TH1F("htt", "", nbins, xmin, xmax);
        hvv = new TH1F("hvv", "", nbins, xmin, xmax);
        hz = new TH1F("hz", "", nbins, xmin, xmax);
    }

    //--- draw histograms
    Ztree->Draw(Form("%s>>hZ",var.c_str()), cuts::Zselection, "goff");
    if (isData == "MC") {
        chtt-> Draw(Form("%s>>htt", var.c_str()), cuts::Zselection*cuts::Zweight, "goff");
        chz->  Draw(Form("%s>>hz", var.c_str()), cuts::Zselection*cuts::Zweight, "goff");
        chvv-> Draw(Form("%s>>hvv", var.c_str()), cuts::Zselection*cuts::Zweight, "goff");
    }
    if (!DF) {
        gtree->Draw(Form("%s>>hg", var.c_str()), cuts::gselection*cuts::weight_g, "goff");
        gtree->Draw(Form("%s>>hg_rw", var.c_str()), cuts::gselection*cuts::weight_g_rw, "goff");
    }

    cout << "" << endl;
    cout << "Z MC integral      " << hZ->Integral()  << endl;
    if (isData == "MC") {
        cout << "tt MC integral       " << htt->Integral()  << endl;
        cout << "Z+jets MC integral   " << hz->Integral()  << endl;
        cout << "VV MC integral       " << hvv->Integral()  << endl;
    }
    cout << "g data raw integral  " << hg->Integral()  << endl;
    cout << "g MC rw integral   " << hg_rw->Integral()  << endl;
    cout << "" << endl;

    //--- normalize Z to MET<60 GeV region
    if (normalize) {
        cout << "normalize to CR " << cuts::CR.GetTitle() << endl;

        TH1F* hZnorm    = new TH1F("hZnorm", "", 1, 0, 1);
        TH1F* hgnorm    = new TH1F("hgnorm", "", 1, 0, 1);
        TH1F* hgrwnorm  = new TH1F("hgrwnorm", "", 1, 0, 1);
        TH1F* httnorm   = new TH1F("httnorm", "", 1, 0, 1);
        TH1F* hvvnorm   = new TH1F("hvvnorm", "", 1, 0, 1);

        Ztree->Draw("0.5>>hZnorm", cuts::Zselection+cuts::CR, "goff");
        chtt-> Draw("0.5>>httnorm", (cuts::Zselection+cuts::CR)*cuts::Zweight, "goff");
        chvv-> Draw("0.5>>hvvnorm", (cuts::Zselection+cuts::CR)*cuts::Zweight, "goff");
        gtree->Draw("0.5>>hgnorm", (cuts::gselection+cuts::CR)*cuts::weight_g, "goff");
        gtree->Draw("0.5>>hgrwnorm", (cuts::gselection+cuts::CR)*cuts::weight_g_rw, "goff");

        float SF   = ( hZnorm->Integral() - httnorm->Integral() - hvvnorm->Integral() ) / hgnorm->Integral();
        float SFrw = ( hZnorm->Integral() - httnorm->Integral() - hvvnorm->Integral() ) / hgrwnorm->Integral();

        cout << "Scale reweighted Z by    " << SFrw << endl;
        cout << "Scale raw Z by           " << SF   << endl;

        hg->Scale(SF);
        hg_rw->Scale(SFrw);
    }

    //if ( TString(var).EqualTo("Z_pt") ) hg->Scale( hg_rw->Integral() / hg->Integral() );

    //--- print MET integrals
    cout << "MET100-150" << endl;
    cout << "2L data                " << hZ->Integral(11,15) << endl;
    cout << "g data (reweighted)    " << hg_rw->Integral(11,15) << endl;
    cout << "g data (raw)           " << hg->Integral(11,15) << endl;
    if (isData == "MC") {
        cout << "VV MC                  " << hvv->Integral(11,15) << endl;
        cout << "tt MC                  " << htt->Integral(11,15) << endl;
        cout << "Z+jets MC              " << hz->Integral(11,15) << endl;
    }

    cout << "MET150-200" << endl;
    cout << "2L data                " << hZ->Integral(16,21) << endl;
    cout << "g data (reweighted)    " << hg_rw->Integral(16,21) << endl;
    cout << "g data (raw)           " << hg->Integral(16,21) << endl;
    if (isData == "MC") {
        cout << "VV MC                  " << hvv->Integral(16,21) << endl;
        cout << "tt MC                  " << htt->Integral(16,21) << endl;
        cout << "Z+jets MC              " << hz->Integral(16,21) << endl;
    }

    //--- make plots
    if (isData == "Data") {
        float integralG, integralZ, normalization; 

        TCanvas *can = new TCanvas("can","can",600,600);
        can->cd();

        TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);
        mainpad->Draw();
        mainpad->cd();

        gPad->SetLogy();

        hZ->SetLineColor(1); hZ->SetLineWidth(2); hZ->SetMarkerStyle(20);

        hZ->GetXaxis()->SetTitle(xtitle.c_str());
        hZ->GetYaxis()->SetTitle("entries / bin");
        hZ->Draw("E1");

        hg_rw->SetLineColor(4);

        THStack *mcstack = new THStack("mcstack","mcstack");
        if( !DF ) mcstack->Add(hg_rw);
        mcstack->Draw("samehist");
        hZ->Draw("sameE1");
        hZ->Draw("axissame");

        TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
        if( !DF ){
            leg->AddEntry(hg_rw,"Z+jets (from #gamma+jets, reweighted)","f");
            leg->AddEntry(hZ,"Z+jets (from MC)","f");
        }
        leg->SetBorderSize(0);
        leg->SetFillColor(0);
        leg->Draw();

        TLatex *tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.03);
        tex->DrawLatex(0.6,0.65,"ATLAS Internal");
        if(TString(channel).Contains("ee")       ) tex->DrawLatex(0.6,0.57,"ee events");
        if(TString(channel).Contains("em")       ) tex->DrawLatex(0.6,0.57,"e#mu events");
        if(TString(channel).Contains("mm")       ) tex->DrawLatex(0.6,0.57,"#mu#mu events");

        can->cd();

        TPad* respad = new TPad("respad","respad",0.0,0.8,1.0,1.0);
        respad->Draw();
        respad->cd();
        respad->SetGridy();

        TH1F* hratio = (TH1F*) hZ->Clone("hratio");
        TH1F* hmctot = (TH1F*) hg_rw->Clone("hmctot");
        for( int ibin = 1 ; ibin <= hmctot->GetXaxis()->GetNbins() ; ibin++ ) hmctot->SetBinError(ibin,0.0);
        hratio->Divide(hmctot);

        hratio->GetXaxis()->SetTitle("");
        hratio->GetXaxis()->SetLabelSize(0.);
        hratio->GetYaxis()->SetNdivisions(5);
        hratio->GetYaxis()->SetTitle("data/bkg");
        hratio->GetYaxis()->SetTitleSize(0.15);
        hratio->GetYaxis()->SetTitleOffset(0.3);
        hratio->GetYaxis()->SetLabelSize(0.15);
        hratio->SetMinimum(0.0);
        hratio->SetMaximum(2.0);
        hratio->GetYaxis()->SetRangeUser(0.0,2.0);
        hratio->Draw("E1");

        can->Print(Form("%s", (plots_path + channel + "_NoSmear_HT_MC_ZptHTreweigh.pdf").c_str()));
    }
    else {
        TCanvas *can = new TCanvas("can","can",600,600);
        can->cd();

        TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);
        mainpad->Draw();
        mainpad->cd();

        gPad->SetLogy();

        hZ->SetLineColor(1); hZ->SetLineWidth(2); hZ->SetMarkerStyle(20);

        hZ->GetXaxis()->SetTitle(xtitle.c_str());
        hZ->GetYaxis()->SetTitle("entries / bin");
        hZ->Draw("E1");

        htt->SetLineColor(1); htt->SetFillColor(kRed-2);
        hvv->SetLineColor(1); hvv->SetFillColor(kGreen-2);
        hg_rw->SetLineColor(1); hg_rw->SetFillColor(kOrange-2);

        hg->Add(htt); hg->Add(hvv);
        hg->SetLineColor(4); hg->SetLineWidth(1); hg->SetLineStyle(2);

        hz->Add(htt); hz->Add(hvv);
        hz->SetLineColor(2); hz->SetLineWidth(1); hz->SetLineStyle(7);

        THStack *mcstack = new THStack("mcstack","mcstack");
        mcstack->Add(htt);
        mcstack->Add(hvv);
        if( !DF ) mcstack->Add(hg_rw);
        mcstack->Draw("samehist");
        hZ->Draw("sameE1");
        if( !DF ) hg->Draw("samehist");
        hZ->Draw("axissame");

        TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
        leg->AddEntry(hZ,"data","lp");
        if( !DF ){
            leg->AddEntry(hg_rw,"Z+jets (from #gamma+jets, reweighted)","f");
            leg->AddEntry(hg,"Z+jets (from #gamma+jets, raw)","f");
            leg->AddEntry(hz,"Z+jets (from MC)","f");
        }
        leg->AddEntry(hvv,"VV","f");
        leg->AddEntry(htt,"t#bar{t}+tW","f");
        leg->SetBorderSize(0);
        leg->SetFillColor(0);
        leg->Draw();

        TLatex *tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.03);
        tex->DrawLatex(0.6,0.65,"ATLAS Internal");
        if(TString(period).Contains("data15-16") ) tex->DrawLatex(0.6,0.61,"36 fb^{-1} 2015-2016 data");
        if(TString(period).Contains("data17")    ) tex->DrawLatex(0.6,0.61,"44 fb^{-1} 2017 data");
        if(TString(channel).Contains("ee")       ) tex->DrawLatex(0.6,0.57,"ee events");
        if(TString(channel).Contains("em")       ) tex->DrawLatex(0.6,0.57,"e#mu events");
        if(TString(channel).Contains("mm")       ) tex->DrawLatex(0.6,0.57,"#mu#mu events");

        can->cd();

        TPad* respad = new TPad("respad","respad",0.0,0.8,1.0,1.0);
        respad->Draw();
        respad->cd();
        respad->SetGridy();

        TH1F* hratio = (TH1F*) hZ->Clone("hratio");
        TH1F* hmctot = (TH1F*) hg_rw->Clone("hmctot");
        hmctot->Add(htt);
        hmctot->Add(hvv);
        for( int ibin = 1 ; ibin <= hmctot->GetXaxis()->GetNbins() ; ibin++ ) hmctot->SetBinError(ibin,0.0);
        hratio->Divide(hmctot);

        hratio->GetXaxis()->SetTitle("");
        hratio->GetXaxis()->SetLabelSize(0.);
        hratio->GetYaxis()->SetNdivisions(5);
        hratio->GetYaxis()->SetTitle("data/bkg");
        hratio->GetYaxis()->SetTitleSize(0.15);
        hratio->GetYaxis()->SetTitleOffset(0.3);
        hratio->GetYaxis()->SetLabelSize(0.15);
        hratio->SetMinimum(0.0);
        hratio->SetMaximum(2.0);
        hratio->GetYaxis()->SetRangeUser(0.0,2.0);
        hratio->Draw("E1");

        can->Print(Form("%s/quickData_Data_%s_%s_%s_%s_VR_ht800cut.pdf",plots_path.c_str(),period.c_str(),channel.c_str(),var.c_str(),smearing_mode.c_str()));
    }
}

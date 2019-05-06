#include "../Settings.C"
#include "../CommonFunctions/CommonLibraries.C"
#include "../CommonFunctions/CommonCuts.C"

using namespace std;

void quickDraw_photonPredictionsSR(string period = "data15-16", string channel = "mm", string var = "MET", string smearing_mode = "NoSmear", string data = "zjets") {

    bool DF = TString(channel).EqualTo("em");

    gStyle->SetOptStat(0);

    //-----------------------------------------------
    // load files
    //-----------------------------------------------

    // set up labels
    string mcdir      = "";
    string gdatalabel = "";

    if     ( TString(period).Contains("data15-16") ){
        mcdir    = "ZMC16a/";
        gdatalabel = "data15-16";
    }
    else if( TString(period).Contains("data17")    ){
        mcdir    = "ZMC16cd/";
        gdatalabel = "data17";
    }
    else if( TString(period).Contains("data18")    ){
        mcdir    = "ZMC16cd/";
        gdatalabel = "data18";
    }

    string gfilename      = reweighting_path + "gdata/" + period + "_merged_processed_" + channel + "_" + smearing_mode + ".root";
    string tt_filename    = ntuple_path + mcdir + "ttbar_merged_processed.root";
    string vv_filename    = ntuple_path + mcdir + "diboson_merged_processed.root";
    string z_filename     = ntuple_path + mcdir + "Zjets_merged_processed.root";

    cout << "period               " << period        << endl;
    cout << "channel              " << channel       << endl;
    cout << "smearing mode        " << smearing_mode << endl;
    cout << "g filename           " << gfilename     << endl;

    //-----------------------------------------------
    // define histograms
    //-----------------------------------------------

    const unsigned int nmetbins =  20;
    double metbins[nmetbins+1] = {10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400};

    string xtitle = var;
    int   nbins =  20;
    float xmin  =   0;
    float xmax  = 400;

    if( TString(var).EqualTo("MET") ){
        xtitle = "E_{T}^{miss} [GeV]";
        nbins = nmetbins;
        xmin  = 0.0;
        xmax  = 400;
    }

    TH1F* h_loose;
    TH1F* h_tight;
    TH1F* h_tighter;
    TH1F* h_tenacious;
    string label;

    if( TString(data).EqualTo("zjets") ){

        TChain* chz = new TChain("BaselineTree");

        chz->Add( z_filename.c_str() );

        cout << "Z+jets entries, loose      " << chz->GetEntries()    << endl;

        TH1F* hz_loose     = new TH1F("hz_loose"     ,"",nmetbins,metbins)  ;
        TH1F* hz_tight     = new TH1F("hz_tight"     ,"",nmetbins,metbins)  ;
        TH1F* hz_tighter   = new TH1F("hz_tighter"     ,"",nmetbins,metbins)    ;
        TH1F* hz_tenacious = new TH1F("hz_tenacious"     ,"",nmetbins,metbins)   ;

        chz->  Draw("min(MET_loose,400)>>hz_loose", cuts::SR*cuts::Zweight      , "goff");
        chz->  Draw("min(MET,400)>>hz_tight", cuts::SR*cuts::Zweight      , "goff");
        chz->  Draw("min(MET_tighter,400)>>hz_tighter", cuts::SR*cuts::Zweight      , "goff");
        chz->  Draw("min(MET_tenacious,400)>>hz_tenacious", cuts::SR*cuts::Zweight      , "goff");

        cout << "MET total" << endl;
        cout << "Z+jets MC integral, loose   " << hz_loose->Integral()  << endl;
        cout << "Z+jets MC integral, tight       " << hz_tight->Integral()  << endl;
        cout << "Z+jets MC integral, tighter   " << hz_tighter->Integral()  << endl;
        cout << "Z+jets MC integral, tenacious       " << hz_tenacious->Integral()  << endl;

        cout << "MET 200-400" << endl;
        cout << "Z+jets MC integral, loose   " << hz_loose->Integral(11,20)  << endl;
        cout << "Z+jets MC integral, tight       " << hz_tight->Integral(11,20)  << endl;
        cout << "Z+jets MC integral, tighter   " << hz_tighter->Integral(11,20)  << endl;
        cout << "Z+jets MC integral, tenacious       " << hz_tenacious->Integral(11,20)  << endl;

        h_loose = hz_loose;
        h_tight = hz_tight;
        h_tighter = hz_tighter;
        h_tenacious = hz_tenacious;
        label = "zmc";
    }
    else if( TString(data).EqualTo("gdata") ){

        TChain* gtree = new TChain("BaselineTree");

        gtree->Add( gfilename.c_str() );

        cout << "g+jets entries   " << gtree->GetEntries()    << endl;

        TH1F* hg_rw_loose     = new TH1F("hg_rw_loose"     ,"",nmetbins,metbins)  ;
        TH1F* hg_rw_tight     = new TH1F("hg_rw_tight"     ,"",nmetbins,metbins)  ;
        TH1F* hg_rw_tighter   = new TH1F("hg_rw_tighter"     ,"",nmetbins,metbins)    ;
        TH1F* hg_rw_tenacious = new TH1F("hg_rw_tenacious"     ,"",nmetbins,metbins)   ;

        gtree->  Draw("min(MET_loose,400)>>hg_rw_loose",         cuts::SR*cuts::weight_g_rw      , "goff");
        gtree->  Draw("min(MET,400)>>hg_rw_tight",         cuts::SR*cuts::weight_g_rw      , "goff");
        gtree->  Draw("min(MET_tighter,400)>>hg_rw_tighter",     cuts::SR*cuts::weight_g_rw  , "goff");
        gtree->  Draw("min(MET_tenacious,400)>>hg_rw_tenacious", cuts::SR*cuts::weight_g_rw , "goff");

        cout << "MET total" << endl;
        cout << "g+jets data integral, loose   "   << hg_rw_loose->Integral()  << endl;
        cout << "g+jets data integral, tight   "   << hg_rw_tight->Integral()  << endl;
        cout << "g+jets data integral, tighter "   << hg_rw_tighter->Integral()  << endl;
        cout << "g+jets data integral, tenacious " << hg_rw_tenacious->Integral()  << endl;

        cout << "MET 200-400" << endl;
        cout << "g+jets data integral, loose     " << hg_rw_loose->Integral(11,20)  << endl;
        cout << "g+jets data integral, tight     " << hg_rw_tight->Integral(11,20)  << endl;
        cout << "g+jets data integral, tighter   " << hg_rw_tighter->Integral(11,20)  << endl;
        cout << "g+jets data integral, tenacious " << hg_rw_tenacious->Integral(11,20)  << endl;

        h_loose = hg_rw_loose;
        h_tight = hg_rw_tight;
        h_tighter = hg_rw_tighter;
        h_tenacious = hg_rw_tenacious;
        label = "gdata";
    }
    else if( TString(data).EqualTo("ttbar") ){

        TChain* chtt = new TChain("BaselineTree");

        chtt->Add( tt_filename.c_str() );

        cout << "ttbar entries,      " << chtt->GetEntries()    << endl;

        TH1F* htt_loose     = new TH1F("htt_loose"     ,"",nmetbins,metbins)  ;
        TH1F* htt_tight     = new TH1F("htt_tight"     ,"",nmetbins,metbins)  ;
        TH1F* htt_tighter   = new TH1F("htt_tighter"     ,"",nmetbins,metbins)    ;
        TH1F* htt_tenacious = new TH1F("htt_tenacious"     ,"",nmetbins,metbins)   ;

        chtt->  Draw("min(MET_loose,400)>>htt_loose", cuts::SR*cuts::Zweight      , "goff");
        chtt->  Draw("min(MET,400)>>htt_tight", cuts::SR*cuts::Zweight      , "goff");
        chtt->  Draw("min(MET_tighter,400)>>htt_tighter", cuts::SR*cuts::Zweight      , "goff");
        chtt->  Draw("min(MET_tenacious,400)>>htt_tenacious", cuts::SR*cuts::Zweight      , "goff");

        cout << "MET total" << endl;
        cout << "ttbar MC integral, loose     " << htt_loose->Integral()  << endl;
        cout << "ttbar MC integral, tight     " << htt_tight->Integral()  << endl;
        cout << "ttbar MC integral, tighter   " << htt_tighter->Integral()  << endl;
        cout << "ttbar MC integral, tenacious " << htt_tenacious->Integral()  << endl;

        cout << "MET 200-400" << endl;
        cout << "ttbar MC integral, loose    " << htt_loose->Integral(11,20)  << endl;
        cout << "ttbar MC integral, tight    " << htt_tight->Integral(11,20)  << endl;
        cout << "ttbar MC integral, tighter  " << htt_tighter->Integral(11,20)  << endl;
        cout << "ttbar MC integral, tenacious" << htt_tenacious->Integral(11,20)  << endl;

        h_loose = htt_loose;
        h_tight = htt_tight;
        h_tighter = htt_tighter;
        h_tenacious = htt_tenacious;
        label = "ttbar";
    }
    else if( TString(data).EqualTo("vv") ){

        TChain* chvv = new TChain("BaselineTree");

        chvv->Add( vv_filename.c_str() );

        cout << "VV entries, loose      " << chvv->GetEntries()    << endl;

        TH1F* hvv_loose     = new TH1F("hvv_loose"     ,"",nmetbins,metbins)  ;
        TH1F* hvv_tight     = new TH1F("hvv_tight"     ,"",nmetbins,metbins)  ;
        TH1F* hvv_tighter   = new TH1F("hvv_tighter"     ,"",nmetbins,metbins)    ;
        TH1F* hvv_tenacious = new TH1F("hvv_tenacious"     ,"",nmetbins,metbins)   ;

        chvv->  Draw("min(MET_loose,400)>>hvv_loose", cuts::SR*cuts::Zweight      , "goff");
        chvv->  Draw("min(MET,400)>>hvv_tight", cuts::SR*cuts::Zweight      , "goff");
        chvv->  Draw("min(MET_tighter,400)>>hvv_tighter", cuts::SR*cuts::Zweight      , "goff");
        chvv->  Draw("min(MET_tenacious,400)>>hvv_tenacious", cuts::SR*cuts::Zweight      , "goff");

        cout << "MET total" << endl;
        cout << "VV MC integral, loose     " << hvv_loose->Integral()  << endl;
        cout << "VV MC integral, tight     " << hvv_tight->Integral()  << endl;
        cout << "VV MC integral, tighter   " << hvv_tighter->Integral()  << endl;
        cout << "VV MC integral, tenacious " << hvv_tenacious->Integral()  << endl;

        cout << "MET 200-400" << endl;
        cout << "VV MC integral, loose     " << hvv_loose->Integral(11,20)  << endl;
        cout << "VV MC integral, tight     " << hvv_tight->Integral(11,20)  << endl;
        cout << "VV MC integral, tighter   " << hvv_tighter->Integral(11,20)  << endl;
        cout << "VV MC integral, tenacious " << hvv_tenacious->Integral(11,20)  << endl;

        h_loose = hvv_loose;
        h_tight = hvv_tight;
        h_tighter = hvv_tighter;
        h_tenacious = hvv_tenacious;
        label = "VV";
    }

    //-----------------------------------------------
    // draw histograms
    //-----------------------------------------------

    // make canvas
    TCanvas *can = new TCanvas("can","can",600,600);
    can->cd();

    TPad* mainpad = new TPad("mainpad","mainpad",0.05,0.0,1.0,0.85);
    mainpad->Draw();
    mainpad->cd();

    gPad->SetLogy();

    // MET loose
    h_loose->SetLineColor(1);
    h_loose->SetLineWidth(2);
    h_loose->SetLineStyle(1);

    h_loose->GetXaxis()->SetTitle(xtitle.c_str());
    h_loose->GetYaxis()->SetTitle("entries / bin");
    h_loose->Draw("hist");
    can->Update();

    // Other MET WP
    h_tight->SetLineColor(4);
    h_tight->SetLineWidth(2);
    h_tight->SetLineStyle(1);
    h_tight->Draw("samehist");

    h_tighter->SetLineColor(3);
    h_tighter->SetLineWidth(2);
    h_tighter->SetLineStyle(1);
    h_tighter->Draw("samehist");

    h_tenacious->SetLineColor(2);
    h_tenacious->SetLineWidth(2);
    h_tenacious->SetLineStyle(1);
    h_tenacious->Draw("samehist");

    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);

    TString plot_label = "";
    if (label == "zmc")
        plot_label = "Z+jets (from MC";
    if (label == "gdata")
        plot_label = "Z+jets (from #gamma+jets";
    if (label == "ttbar")
        plot_label = "t#bar{t} (from MC";
    if (label == "VV")
        plot_label = "VV (from MC";

    leg->AddEntry(h_loose, plot_label + ", MET loose)","f");
    leg->AddEntry(h_tight, plot_label + ", MET tight)","f");
    leg->AddEntry(h_tighter, plot_label + ", MET tighter)","f");
    leg->AddEntry(h_tenacious, plot_label + ", MET tenacious)","f");

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

    // Tight/loose ratio

    TCanvas *can2 = new TCanvas("can2","can2",600,600);
    can2->cd();

    TPad* respad = new TPad("respad","respad",0.05,0.7,1.0,0.95);
    respad->Draw();
    respad->cd();
    respad->SetGridy();

    TH1F* hratio = (TH1F*) h_tight->Clone("hratio");
    TH1F* hrLoose = (TH1F*) h_loose->Clone("hrLoose");

    for( int ibin = 1 ; ibin <= hrLoose->GetXaxis()->GetNbins() ; ibin++ ) hrLoose->SetBinError(ibin,0.0);
    hratio->Divide(hrLoose);

    hratio->GetXaxis()->SetTitle("");
    hratio->GetXaxis()->SetLabelSize(0.);
    hratio->GetYaxis()->SetNdivisions(5);
    hratio->GetYaxis()->SetTitle("E^{tight}_{T,miss}/E^{loose}_{T,miss}");
    hratio->GetYaxis()->SetTitleSize(0.1);
    hratio->GetYaxis()->SetTitleOffset(0.5);
    hratio->GetYaxis()->SetLabelSize(0.1);
    hratio->SetMinimum(0.0);
    hratio->SetMaximum(2.0);
    hratio->GetYaxis()->SetRangeUser(0.,2.0);
    hratio->Draw("E1");

    can2->cd();

    // Tighter ratio

    TPad* respad1 = new TPad("respad1","respad1",0.05,0.35,1.0,0.6);
    respad1->Draw();
    respad1->cd();
    respad1->SetGridy();

    TH1F* hratio1 = (TH1F*) h_tighter->Clone("hratio1");

    hratio1->Divide(hrLoose);

    hratio1->GetXaxis()->SetTitle("");
    hratio1->GetXaxis()->SetLabelSize(0.);
    hratio1->GetYaxis()->SetNdivisions(5);
    hratio1->GetYaxis()->SetTitle("E^{tighter}_{T,miss}/E^{loose}_{T,miss}");
    hratio1->GetYaxis()->SetTitleSize(0.1);
    hratio1->GetYaxis()->SetTitleOffset(0.5);
    hratio1->GetYaxis()->SetLabelSize(0.1);
    hratio1->SetMinimum(0.0);
    hratio1->SetMaximum(2.0);
    hratio1->GetYaxis()->SetRangeUser(0.0,2.0);
    hratio1->Draw("E1");

    can->Print(Form("%sVR_SR_studies/quickDraw_%s_%s_%s_%s_SR_%s_metWPs.pdf",plots_path.c_str(),period.c_str(),channel.c_str(),var.c_str(),smearing_mode.c_str()),label.c_str());

    can2->cd();

    //Tenacious ratio
    TPad* respad2 = new TPad("respad2","respad2",0.05,0.05,1.0,0.3);
    respad2->Draw();
    respad2->cd();
    respad2->SetGridy();

    TH1F* hratio2 = (TH1F*) h_tenacious->Clone("hratio2");

    hratio2->Divide(hrLoose);

    hratio2->GetXaxis()->SetTitle("");
    hratio2->GetXaxis()->SetLabelSize(0.);
    hratio2->GetYaxis()->SetNdivisions(5);
    hratio2->GetYaxis()->SetTitle("E^{tenacious}_{T,miss}/E^{loose}_{T,miss}");
    hratio2->GetYaxis()->SetTitleSize(0.1);
    hratio2->GetYaxis()->SetTitleOffset(0.5);
    hratio2->GetYaxis()->SetLabelSize(0.1);
    hratio2->SetMinimum(0.0);
    hratio2->SetMaximum(2.0);
    hratio2->GetYaxis()->SetRangeUser(0.,2.0);
    hratio2->Draw("E1");

    can2->Print(Form("%sVR_SR_studies/ratioPlots_%s_%s_%s_%s_SR_%s_metWPs.pdf",plots_path.c_str(),period.c_str(),channel.c_str(),var.c_str(),smearing_mode.c_str()),label.c_str());
}

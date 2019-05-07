#include "../Common/Settings.C"
#include "../Common/CommonLibraries.C"
#include "../Common/CommonCuts.C"

using namespace std;

TH1F* GetSimpleReweightingHistograms(string period, string channel, string smearing_mode, int step ){

    cout << "Doing histograms for period and year " << period << " " << channel << endl;
    gStyle->SetOptStat(0);

    //--- set up labels
    string mcdir      = "";
    string gdatalabel = "";
    if     ( TString(period).Contains("data15-16") ){
        mcdir    = "MC16a/";
        gdatalabel = "data15-16";
    }
    else if( TString(period).Contains("data17")    ){
        mcdir    = "MC16cd/";
        gdatalabel = "data17";
    }
    else if( TString(period).Contains("data18")    ){
        mcdir    = "MC16cd_2018/";
        gdatalabel = "data18";
    }

    //--- set up MC period
    string mcperiod = "";
    if( TString(period).Contains("data15-16") ) mcperiod = "ZMC16a/";
    if( TString(period).Contains("data17")    ) mcperiod = "ZMC16cd/";
    if( TString(period).Contains("data18")    ) mcperiod = "ZMC16cd/";

    // set up filenames
    string data_filename   = ntuple_path + "zdata/" + period + "_merged_processed.root";
    string tt_filename     = ntuple_path + mcperiod + "ttbar_merged_processed.root";
    string vv_filename     = ntuple_path + mcperiod + "diboson_merged_processed.root";
    string zjets_filename  = ntuple_path + mcperiod + "Zjets_merged_processed.root";
    string photon_filename = reweighting_path + "gdata/" + period + "_merged_processed" + "_" + channel + "_" + smearing_mode + ".root"; //Vg subtracted 

    cout << "Opening data file    " << data_filename   << endl;
    cout << "Opening ttbar file   " << tt_filename     << endl;
    cout << "Opening diboson file " << vv_filename     << endl;
    cout << "Opening Z+jets file  " << zjets_filename << endl;
    cout << "Opening photon file  " << photon_filename << endl;

    // create TChains
    TChain* chdata  = new TChain("BaselineTree");
    TChain* chtt    = new TChain("BaselineTree");
    TChain* chvv    = new TChain("BaselineTree");
    TChain* chzjets = new TChain("BaselineTree");
    TChain* chg     = new TChain("BaselineTree");

    // add files to TChain and get entries
    chdata->Add( data_filename.c_str() );
    chtt->Add( tt_filename.c_str() );
    chvv->Add( vv_filename.c_str() );
    chzjets->Add( zjets_filename.c_str() );
    chg->Add( photon_filename.c_str() );

    cout << "data entries         " << chdata->GetEntries()   << endl;
    cout << "ttbar entries        " << chtt->GetEntries()     << endl;
    cout << "diboson entries      " << chvv->GetEntries()     << endl;
    cout << "Z+jets entries       " << chzjets->GetEntries()  << endl;
    cout << "photon entries       " << chg->GetEntries()      << endl;

    // define selections and weights
    if ( TString(channel).EqualTo("ee") ) cuts::Zselection += cuts::ee;
    else if ( TString(channel).EqualTo("mm") ) cuts::Zselection += cuts::mm;
    else {
        cout << "Unrecognized channel! quitting   " << channel << endl;
        exit(0);
    }

    float lumi = GetLumi(period);
    cout << "Z selection          " << cuts::Zselection.GetTitle() << endl;
    cout << "g selection          " << cuts::gselection.GetTitle() << endl;
    cout << "weight               " << cuts::Zweight.GetTitle()     << endl;
    cout << "lumi                 " << lumi       << endl;

    // define histograms
    TH1F* hdata  = new TH1F("hdata"  ,"",nptbins,ptbins);
    TH1F* htt    = new TH1F("htt"    ,"",nptbins,ptbins);
    TH1F* hvv    = new TH1F("hvv"    ,"",nptbins,ptbins);
    TH1F* hz     = new TH1F("hz"     ,"",nptbins,ptbins);
    TH1F* histoG = new TH1F("histoG" ,"",nptbins,ptbins);    

    TCut RunRange("");
    if( TString(period).EqualTo("data17")    ){
        RunRange = TCut("RunNumber < 348000");  
        cout << "Data17! adding cut " << RunRange.GetTitle() << endl;
    }

    TCut lumi_cut = TCut(TString(to_string(lumi)));

    // fill histograms: HT -step1
    if (step == 1) {
        chdata-> Draw("min(HT,999)>>hdata"  ,cuts::Zselection       ,"goff");
        chtt->   Draw("min(HT,999)>>htt"    ,cuts::Zselection*RunRange*cuts::Zweight*lumi_cut,"goff");
        chvv->   Draw("min(HT,999)>>hvv"    ,cuts::Zselection*RunRange*cuts::Zweight*lumi_cut,"goff");
        chzjets->Draw("min(HT,999)>>hz"     ,cuts::Zselection*RunRange*cuts::Zweight*lumi_cut,"goff");
        chg    ->Draw("min(HT,999)>>histoG" ,cuts::gselection*cuts::weight_g,"goff");
    }

    // fill histograms: Z_pt -step2
    else if (step == 2) {
        TCut g_rw("ptreweight_step1"); // from step 1
        chdata-> Draw("min(Z_pt,999)>>hdata"  ,cuts::Zselection       ,"goff");
        chtt->   Draw("min(Z_pt,999)>>htt"    ,cuts::Zselection*RunRange*cuts::Zweight*lumi_cut,"goff");
        chvv->   Draw("min(Z_pt,999)>>hvv"    ,cuts::Zselection*RunRange*cuts::Zweight*lumi_cut,"goff");
        chzjets->Draw("min(Z_pt,999)>>hz"     ,cuts::Zselection*RunRange*cuts::Zweight*lumi_cut,"goff");
        chg    ->Draw("min(Z_pt,999)>>histoG" ,cuts::gselection*cuts::weight_g*g_rw,"goff");
    }

    cout << "data integral        " << hdata->Integral()   << endl;
    cout << "ttbar integral       " << htt->Integral()     << endl;
    cout << "diboson integral     " << hvv->Integral()     << endl;
    cout << "Z+jets integral      " << hz->Integral()      << endl;
    cout << "photon integral      " << histoG->Integral()  << endl;

    // make canvas and draw 2L data vs. MC plot
    TCanvas *can = new TCanvas("can","can",600,600);
    can->cd();

    gPad->SetLogy();

    hdata->SetLineColor(1);
    hdata->SetLineWidth(2);
    hdata->SetMarkerStyle(20);

    hdata->GetXaxis()->SetTitle("Z p_{T} [GeV]");
    hdata->GetYaxis()->SetTitle("entries / bin");
    hdata->Draw("E1");

    htt->SetLineColor(1);
    htt->SetFillColor(kRed-2);

    hvv->SetLineColor(1);
    hvv->SetFillColor(kGreen-2);

    hz->SetLineColor(1);
    hz->SetFillColor(kOrange-2);

    THStack *mcstack = new THStack("mcstack","mcstack");
    mcstack->Add(htt);
    mcstack->Add(hvv);
    mcstack->Add(hz);
    mcstack->Draw("samehist");
    hdata->Draw("sameE1");
    hdata->Draw("axissame");

    TLegend* leg = new TLegend(0.65,0.65,0.88,0.88);
    leg->AddEntry(hdata,"data","lp");
    leg->AddEntry(hz,"Z+jets","f");
    leg->AddEntry(hvv,"VV","f");
    leg->AddEntry(htt,"t#bar{t}","f");
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->Draw();
    can->Print(Form("%sGetSimpleReweightingHistograms_%s_%s_%s_2L.pdf",plots_path.c_str(),period.c_str(),channel.c_str(),smearing_mode.c_str()));

    //--- take data and remove MC htt and MC hvv
    TH1F* histoZ = (TH1F*) hdata->Clone("histoZ");
    histoZ->Add( htt , -1.0 );
    histoZ->Add( hvv , -1.0 );
    //histoZ->Add( hz  , -1.0 );

    TCanvas *can2 = new TCanvas("can2","can2",600,600);
    can2->cd();

    // float nZ = histoZ->Integral();
    // histoZ->Scale( 1.0 / nZ );

    // float nG = histoG->Integral();
    // histoG->Scale( 1.0 / nG );

    TH1F* hratio = (TH1F*) histoZ->Clone("hratio");
    hratio->Divide( histoG );

    can2->Divide(1,2);
    can2->cd(2);

    histoG->GetXaxis()->SetTitle("Z p_{T} [GeV]");
    histoG->GetYaxis()->SetTitle("entries / bin");

    gPad->SetLogy();

    histoZ->SetLineColor(2);
    histoG->SetLineColor(4);

    histoG->Draw("hist");
    histoZ->Draw("samehist");

    TLegend* leg2 = new TLegend(0.7,0.7,0.88,0.88);
    leg2->AddEntry(histoZ,"2L data - t#bar{t} - VV","f");
    leg2->AddEntry(histoG,"photon","f");
    leg2->SetBorderSize(0);
    leg2->SetFillColor(0);
    leg2->Draw();

    can2->cd(1);
    hratio->SetLineColor(1);
    hratio->Draw("hist");

    can2->Print(Form("%sGetSimpleReweightingHistograms_%s_%s_%s_Z_vs_g.pdf",plots_path.c_str(),period.c_str(),channel.c_str(),smearing_mode.c_str()));

    cout << "histoG->Integral() " << histoG->Integral() << endl;
    cout << "histoZ->Integral() " << histoZ->Integral() << endl;
    cout << "hratio->Integral() " << hratio->Integral() << endl;

    return hratio;
}

#include "../Common/Settings.C"
#include "../Common/CommonLibraries.C"
#include "../Common/CommonCuts.C"

using namespace std;

void quickDraw_MC(string period = "data15-16" , string channel  = "mm" , string var = "HT", string smearing_mode = "NoSmear" ) {

    bool DF = TString(channel).EqualTo("em");

    gStyle->SetOptStat(0);

    //-----------------------------------------------
    // load files
    //-----------------------------------------------

    string Zfilename      = ntuple_path + "ZMC16a/Zjets_merged_processed.root";
    string gfilename      = reweighting_path + "gmc/gmc_" + channel + "_" + smearing_mode + ".root";

    cout << "channel              " << channel       << endl;
    cout << "smearing mode        " << smearing_mode << endl;
    cout << "Z datafilename       " << Zfilename     << endl;
    cout << "g filename           " << gfilename     << endl;
    cout << "DF?                  " << DF            << endl;
    cout << "" << endl;

    //-----------------------------------------------
    // add files to TChain
    //-----------------------------------------------

    TChain * gtree = new TChain("BaselineTree");
    if( !DF ) gtree->Add( gfilename.c_str() );

    TChain * Ztree = new TChain("BaselineTree");
    Ztree->Add( Zfilename.c_str() );

    cout << "g entries            " << gtree->GetEntries()  << endl;
    cout << "Z data entries       " << Ztree->GetEntries()  << endl;

    //-----------------------------------------------
    // Define selections
    //-----------------------------------------------

    if     ( TString(channel).EqualTo("ee") ) cuts::Zselection += cuts::ee;
    else if( TString(channel).EqualTo("mm") ) cuts::Zselection += cuts::mm;
    else if( TString(channel).EqualTo("em") ) cuts::Zselection += cuts::em;
    else{
        cout << "Unrecognized channel! quitting   " << channel << endl;
        exit(0);
    }

    //-----------------------------------------------
    // Set weights 
    //-----------------------------------------------

    TCut lumi("1.0");
    if( TString(period).EqualTo("data15-16") ) lumi = TCut("36100");
    if( TString(period).EqualTo("data17")    ) lumi = TCut("44100");
    if( TString(period).EqualTo("data18")    ) lumi = TCut("64610");

    cout << "Z selection          " << cuts::Zselection.GetTitle()  << endl;  
    cout << "Z weight             " << cuts::Zweight.GetTitle()     << endl;
    cout << "g selection          " << cuts::gselection.GetTitle()  << endl;
    cout << "g weight             " << cuts::weight_g.GetTitle()    << endl;
    cout << "g weight (reweight)  " << cuts::weight_g_rw.GetTitle() << endl;
    cout << "luminosity           " << lumi.GetTitle()       << endl;

    //-----------------------------------------------
    // define and draw histograms
    //-----------------------------------------------

    std::tuple<string, int, float, float> plot_settings;

    if (var == "MET") plot_settings = std::make_tuple("E_{T}^{miss} [GeV]", 20, 0, 400);
    else if (var == "METl") plot_settings = std::make_tuple("E_{T,||}^{miss} [GeV]", 20, -200, 200);
    else if (var == "METt") plot_settings = std::make_tuple("E_{T,#perp}^{miss} [GeV]", 20, -200, 200);
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

    // Initialize histograms

    TH1F* hZ;
    TH1F* hg;
    TH1F* hg_rw ;

    if(var == "Z_pt") {
        const unsigned int nZptbins = 16;
        double Zptbins[nZptbins+1] = {40,75,100,125,150,175,200,250,300,350,400,450,500,600,700,850,1000};
        hZ = new TH1F("hZ", "", nZptbins, Zptbins);
        hg_rw = new TH1F("hg_rw", "", nZptbins, Zptbins);
    }
    else if(var == "HT") {
        const unsigned int nHTbins = 16;
        double HTbins[nHTbins+1] = {40,75,100,125,150,175,200,250,300,350,400,450,500,600,700,850,1000};
        hZ = new TH1F("hZ", "", nHTbins, HTbins);
        hg_rw = new TH1F("hg_rw", "", nHTbins, HTbins);
    }
    else {
        hZ = new TH1F("hZ", "", nbins, xmin, xmax);
        hg_rw = new TH1F("hg_rw", "", nbins, xmin, xmax);
    }

    Ztree->Draw(Form("%s>>hZ",var.c_str())       , cuts::Zselection*cuts::Zweight*lumi      , "goff");

    if( !DF ){
        //gtree->Draw(Form("%s>>hg",var.c_str())     , cuts::gselection*cuts::weight_g*lumi     , "goff");
        gtree->Draw(Form("%s>>hg_rw",var.c_str())  , cuts::gselection*cuts::weight_g_rw*lumi  , "goff");
    }

    cout << "" << endl;
    cout << "Z MC integral      " << hZ->Integral()  << endl;
    //cout << "g MC raw integral  " << hg->Integral()  << endl;
    cout << "g MC rw integral   " << hg_rw->Integral()  << endl;
    cout << "" << endl;

    //-----------------------------------------------
    // make pretty plots
    //-----------------------------------------------

    cout << "MET100-150" << endl;
    cout << "Z                      " << hZ->Integral(11,15)     << endl;
    cout << "g MC (reweighted)    " << hg_rw->Integral(11,15)  << endl;
    //cout << "g MC (raw)           " << hg->Integral(11,15)     << endl;

    cout << "MET150-200" << endl;
    cout << "Z                      " << hZ->Integral(16,21)     << endl;
    cout << "g MC (reweighted)    " << hg_rw->Integral(16,21)  << endl;
    //cout << "g MC (raw)           " << hg->Integral(16,21)     << endl;

    //-------------------------------------------
    //--------- Normalize histograms ------------ 
    //-------------------------------------------

    float integralG; 
    float integralZ; 
    float normalization; 
    /*
       integralG = hg_rw->Integral(1,3);
       integralZ = hZ->Integral(1,3);
    //IF USING HT VAR 
    //integralG = hg_rw->Integral(10,15);
    //integralZ = hZ->Integral(10,15);
    normalization = integralZ/integralG;
    hg_rw->Scale(normalization);
    */

    TCanvas *can = new TCanvas("can","can",600,600);
    can->cd();

    TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);
    mainpad->Draw();
    mainpad->cd();

    gPad->SetLogy();

    hZ->SetLineColor(1);
    hZ->SetLineWidth(2);
    hZ->SetMarkerStyle(20);

    hZ->GetXaxis()->SetTitle(xtitle.c_str());
    hZ->GetYaxis()->SetTitle("entries / bin");
    hZ->Draw("E1");

    //hg->SetLineColor(2);
    hg_rw->SetLineColor(4);

    THStack *mcstack = new THStack("mcstack","mcstack");
    if( !DF ) mcstack->Add(hg_rw);
    mcstack->Draw("samehist");
    hZ->Draw("sameE1");
    hZ->Draw("axissame");

    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
    //leg->AddEntry(hZ,"data","lp");
    if( !DF ){
        leg->AddEntry(hg_rw,"Z+jets (from #gamma+jets, reweighted)","f");
        //leg->AddEntry(hg,"Z+jets (from #gamma+jets, raw)","f");
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

    can->Print(Form((plots_path + "mm_NoSmear_HT_MC_ZptHTreweigh.pdf").c_str()));

}

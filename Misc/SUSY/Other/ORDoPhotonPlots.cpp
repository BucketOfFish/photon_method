void ORDoPhotonPlots() {
    // Options
    TString path = "/eos/user/m/mazhang/SUSY/finished/";
    TString file = "Zee/_T_06_09/data-tree/mc16_13TeV.364121.Sherpa_221_NNPDF30NNLO_Zee_MAXHTPTV140_280_CFilterBVeto.deriv.DAOD_SUSY2.e5299_s3126_r10724_p4189.root";
    bool is_data = false;

    TFile *noORDoPhoton = new TFile(path + "noORDoPhoton/" + file);
    TFile *yesORDoPhoton = new TFile(path + "yesORDoPhoton/" + file);
    TTree *noTree = (TTree*)noORDoPhoton->Get("tree_NoSys");
    TTree *yesTree = (TTree*)yesORDoPhoton->Get("tree_NoSys");

    TCut lumi = "(RandomRunNumber<320000 ? 36200 : (RandomRunNumber>320000 && RandomRunNumber<348000) ? 44300 : 58500)";
    TCut MC_weight = "genWeight*eventWeight*leptonWeight*jvtWeight*bTagWeight*pileupWeight*globalDiLepTrigSF";

    TString feature = "met_Et";
    TString plotName = "Zee.png";

    // Make Plots
    TCanvas *can = new TCanvas("can","can",800,600);
    can->cd();
    TPad* namepad = new TPad("namepad","namepad",0.0,0.0,1.0,1.0);
    namepad->Draw();
    namepad->cd();

    TCut cut = lumi*MC_weight;
    if (is_data) cut = "1";

    TH1F* hist_yes = new TH1F("hist_yes","",10,0,3.14159);
    TH1F* hist_no = new TH1F("hist_no","",10,0,3.14159);
    yesTree->Draw(feature + ">>hist_yes", cut, "goff");
    noTree->Draw(feature + ">>hist_no", cut, "goff");

    hist_yes->SetLineColor(kRed);
    hist_yes->SetLineWidth(2);
    hist_no->SetLineColor(kBlue);
    hist_no->SetLineWidth(2);

    TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
    leg->AddEntry(hist_yes, "DoPhoton true", "lp");
    leg->AddEntry(hist_no, "DoPhoton false", "lp");

    namepad->SetLogy();
    hist_yes->SetTitle("Effect of OR.DoPhoton on MET");
    hist_yes->GetXaxis()->SetTitle(feature);
    hist_yes->GetYaxis()->SetTitle("entries / bin");
    hist_yes->Draw();
    hist_no->Draw();
    leg->Draw();

    can->Print(plotName);

    delete can, namepad;
}

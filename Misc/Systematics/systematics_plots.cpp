#include "../../Main/Settings.cpp"

void systematics_plots() {
    // Options
    TChain chain("Zjets_NoSys");
    chain.Add("/public/data/SUSY_Systematics/Skimmed/StrongPreselectionInclusive/mc16a_Zjets_merged_processed.root");
    chain.Add("/public/data/SUSY_Systematics/Skimmed/StrongPreselectionInclusive/mc16cd_Zjets_merged_processed.root");
    chain.Add("/public/data/SUSY_Systematics/Skimmed/StrongPreselectionInclusive/mc16e_Zjets_merged_processed.root");

    vector<string> regions = {"SRC", "SRLow", "SRMed"}; 
    vector<string> systematics = {"nominal", "scale", "PDF"};
    map<string, vector<string>> largest_systematics;
    string nominal_branch = "MUR1_MUF1_PDF261000";
    largest_systematics["SRC"] = {nominal_branch, "MUR0_5_MUF1_PDF261000", "MUR1_MUF1_PDF13000"};
    largest_systematics["SRLow"] = {nominal_branch, "MUR0_5_MUF0_5_PDF261000", "MUR1_MUF1_PDF25300"};
    largest_systematics["SRMed"] = {nominal_branch, "MUR0_5_MUF1_PDF261000", "MUR1_MUF1_PDF25300"};
    vector<int> colors = {kBlack, kRed, kBlue};

    TCut lumi = "(RandomRunNumber<320000 ? 36200 : (RandomRunNumber>320000 && RandomRunNumber<348000) ? 44300 : 58500)";
    TCut MC_weight = "genWeight*eventWeight*leptonWeight*jvtWeight*bTagWeight*pileupWeight*globalDiLepTrigSF";

    // Make Plots
    for (auto region : regions) {
        cout << "Making Plots for Region " << region << endl;
        TCut inclusive_cut = NMinus1Cut(cuts::selections[region], "minDPhi2JetsMet");

        TCanvas *can = new TCanvas("can","can",800,600);
        can->cd();
        TPad* namepad = new TPad("namepad","namepad",0.0,0.0,1.0,1.0);
        namepad->Draw();
        namepad->cd();

        THStack *stack = new THStack("dPhi stack", "");
        TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
        int color = 0;
        for (auto LHE3Weight : largest_systematics[region]) {
            TH1F* hist = new TH1F(("hist_" + LHE3Weight).c_str(),"",10,0,3.14159);
            TCut LHE_branch = ("LHE3Weight_" + LHE3Weight).c_str();
            chain.Draw(("minDPhi2JetsMet>>hist_" + LHE3Weight).c_str(), inclusive_cut*lumi*MC_weight*LHE_branch, "goff");
            hist->SetLineColor(colors[color++]);
            hist->SetLineWidth(2);
            hist->GetXaxis()->SetTitle("minDPhi2JetsMet");
            hist->GetYaxis()->SetTitle("entries / bin");
            stack->Add(hist);
            leg->AddEntry(hist, LHE_branch.GetTitle(), "lp");
        }

        namepad->SetLogy();
        stack->SetTitle(("dPhi Shapes in " + region).c_str());
        stack->Draw("hist nostack");
        leg->Draw();

        TString plot_name = (region + ".png").c_str();
        can->Print(plot_name);

        for (auto hist : *stack) delete hist;
        delete can, namepad, stack;
    }	
}

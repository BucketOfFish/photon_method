#include "../../PhotonMethod/Main/Settings.cpp"

void make_2D_eta_phi_plot() {
    // options
    string tree_name = "BaselineTree";
    string path = "/public/data/Photon/Samples/ReducedNtuples/";
    vector<string> trees = {
        path + "mc16a_Zjets.root",
        path + "mc16cd_Zjets.root",
        path + "mc16e_Zjets.root",
    };
    vector<vector<string>> vars = {
        {"lepEta[0]", "lepPhi[0]", "Leading Lepton Eta vs. Phi", "lep_eta_phi"},
        {"jet_eta[0]", "jet_phi[0]", "Leading Jet Eta vs. Phi", "jet_eta_phi"},
    };
    string weight = "totalWeight";
    vector<string> regions = {"SRC", "SRLow", "SRMed", "SRHigh", "SRLowZ", "SRMedZ", "SRHighZ", "VRC_CR", "VRC",
        "VRLow", "VRMed", "VRHigh", "VRLowZ", "VRMedZ", "VRHighZ",
        "SRHigh4", "SRLow2", "VRHigh4", "EWK_VRHigh", "VRHighR", "VRllbb", "VRInt",
        "EWK_VRLow", "VRLow2", "VROffShell", "CRZZ", "CRtt", "CRZ", "CRDY", "DRHigh", "DRllbb", "DRInt", "DRLow",
        "DROffShell", "SRHigh8", "SRHigh16", "SRInt", "SRLow", "SROffShell"};

    // plotting
    TChain *chain = new TChain(tree_name.c_str());
    for (auto tree : trees) {
        chain->Add(tree.c_str());
    }
    TCanvas *can = new TCanvas("can","can",800,600);
    can->cd();
    gStyle->SetOptStat(0);
    for (auto region : regions) {
        for (auto var : vars) {
            TH2D *hist = new TH2D("hist", (var[2] + " for Zjets MC").c_str(), 20, -3.14159, 3.14159, 20, -2.74, 2.74);
            chain->Draw((var[0]+":"+var[1]+">>hist").c_str(), cuts::selections[region] * TCut(weight.c_str()), "colz");
            hist->GetXaxis()->SetTitle(var[1].c_str());
            hist->GetYaxis()->SetTitle(var[0].c_str());
            can->Print((var[3] + "_" + region + ".png").c_str());
            delete hist;
        }
    }
}

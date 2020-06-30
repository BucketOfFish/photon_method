#include "../../../Main/Settings.cpp"

void make_plots(string process, vector<string> trees, vector<string> regions) {
    // options
    string tree_name = "BaselineTree";
    vector<vector<string>> vars = {
        {"lepEta[0]", "lepPhi[0]", "Leading Lepton Eta vs. Phi", "lep_eta_phi"},
        {"jet_eta[0]", "jet_phi[0]", "Leading Jet Eta vs. Phi", "jet_eta_phi"},
    };
    string weight = "totalWeight";

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
            string title = var[2] + " for " + process + " MC in " + region;
            if (process == "data") title = var[2] + " for data in " + region;
            TH2D *hist = new TH2D("hist", title.c_str(), 20, -3.14159, 3.14159, 20, -2.74, 2.74);
            if (process != "data") chain->Draw((var[0]+":"+var[1]+">>hist").c_str(), cuts::selections[region] * TCut(weight.c_str()), "colz");
            else chain->Draw((var[0]+":"+var[1]+">>hist").c_str(), cuts::selections[region], "colz");
            hist->GetXaxis()->SetTitle(var[1].c_str());
            hist->GetYaxis()->SetTitle(var[0].c_str());
            can->Print((var[3] + "_" + process + "_" + region + ".png").c_str());
            delete hist;
        }
    }
}

void make_2D_eta_phi_plot() {
    string mc_path = "/public/data/Photon/Samples/ReducedNtuples/";
    vector<string> processes = {"diboson", "higgs", "lowMassDY", "singleTop", "topOther", "triboson", "ttbar",
        "Vgamma", "Wjets", "Zjets"};
    vector<string> regions = {"SRC", "SRLow", "SRMed", "SRHigh", "SRLowZ", "SRMedZ", "SRHighZ", "VRC_CR", "VRC",
        "VRLow", "VRMed", "VRHigh", "VRLowZ", "VRMedZ", "VRHighZ",
        "CRC", "CRLow", "CRMed", "CRHigh", "CRLowZ", "CRMedZ", "CRHighZ",
        "SRHigh4", "SRLow2", "VRHigh4", "EWK_VRHigh", "VRHighR", "VRllbb", "VRInt",
        "EWK_VRLow", "VRLow2", "VROffShell", "CRZZ", "CRtt", "CRZ", "CRDY", "DRHigh", "DRllbb", "DRInt", "DRLow",
        "DROffShell", "SRHigh8", "SRHigh16", "SRInt", "SRLow", "SROffShell"};
    for (auto process : processes) {
        vector<string> trees = {
            mc_path + "mc16a_" + process + ".root",
            mc_path + "mc16cd_" + process + ".root",
            mc_path + "mc16e_" + process + ".root",
        };
        make_plots(process, trees, regions);
    }

    string data_path = "/public/data/Photon/Samples/ReducedNtuples/";
    vector<string> trees = {
        data_path + "data15-16_data_bkg.root",
        data_path + "data17_data_bkg.root",
        data_path + "data18_data_bkg.root",
    };
    regions = {"VRC_CR", "VRC",
        "VRLow", "VRMed", "VRHigh", "VRLowZ", "VRMedZ", "VRHighZ",
        "CRC", "CRLow", "CRMed", "CRHigh", "CRLowZ", "CRMedZ", "CRHighZ",
        "VRHigh4", "EWK_VRHigh", "VRHighR", "VRllbb", "VRInt",
        "EWK_VRLow", "VRLow2", "VROffShell", "CRZZ", "CRtt", "CRZ", "CRDY", "DRHigh", "DRllbb", "DRInt", "DRLow",
        "DROffShell"};
    make_plots("data", trees, regions);
}

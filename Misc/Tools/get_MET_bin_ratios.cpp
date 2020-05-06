#include "../../Main/Settings.cpp"

void get_MET_bin_ratios() {
    //--- get photon and ZMC files
    vector<string> periods = vector<string>{"data15-16", "data17", "data18"};
    TChain *tch_photon = new TChain("BaselineTree");
    TChain *tch_zmc = new TChain("BaselineTree");

    for (auto period : periods) {
        string data_period = DataPeriod(period);
        string photon_path = "/public/data/Photon/Samples/ReweightedNtuples/";
        string zmc_path = "/public/data/Photon/Samples/ReducedNtuples/";
        tch_photon->Add((photon_path + data_period + "_data_photon_ee.root").c_str());
        tch_photon->Add((photon_path + data_period + "_data_photon_mm.root").c_str());
        tch_zmc->Add((zmc_path + data_period + "_data_bkg.root").c_str());
    }

    //--- options
    string reweight_branch = "reweight_Ptll__Ht30";

    //--- get ratio
    TH1F *photon_yield = new TH1F("photon_yield", "", 1,0,1);
    TH1F *zmc_yield = new TH1F("zmc_yield", "", 1,0,1);
    tch_photon->Draw("mll>>photon_yield", ("(met_Et<100) * totalWeight * " +  reweight_branch).c_str(), "goff");
    tch_zmc->Draw("mll>>zmc_yield", "(met_Et<100) * totalWeight", "goff");
    cout << "Ratio below met_Et<100: " << photon_yield->Integral(0,2) / zmc_yield->Integral(0,2) << endl;
    tch_photon->Draw("mll>>photon_yield", ("(met_Et>100) * totalWeight * " + reweight_branch).c_str(), "goff");
    tch_zmc->Draw("mll>>zmc_yield", "(met_Et>100) * totalWeight", "goff");
    cout << "Ratio above met_Et>100: " << photon_yield->Integral(0,2) / zmc_yield->Integral(0,2) << endl;
}

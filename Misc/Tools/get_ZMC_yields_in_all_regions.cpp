#include "../../Main/Settings.cpp"
#include <fstream>

void get_ZMC_yields_in_all_regions() {
    //--- settings
    vector<string> regions = vector<string>{"SRC", "SRLow2", "SRMed2", "SRHigh2", "SRLowZ4", "SRMedZ4", "SRHighZ4",
                                    "VRC", "VRLow2", "VRMed2", "VRHigh2", "VRLowZ4", "VRMedZ4", "VRHighZ4"};

    //--- get ZMC files
    vector<string> periods = vector<string>{"mc16a", "mc16cd", "mc16e"};
    TChain *tch_zmc = new TChain("BaselineTree");
    for (auto period : periods) {
        string zmc_path = "/public/data/Photon/Samples/ReducedNtuples/";
        string zmc_filename = zmc_path + period + "_Zjets.root";
        tch_zmc->Add(zmc_filename.c_str());
    }

    //--- output file
    ofstream outfile;
    outfile.open("ZMC_yields.txt");

    //--- get yields
    TH1D *yield = new TH1D("yield", "", 1, 0, 1);
    for (string region_name : regions) {
        TCut region = TCut(cuts::selections[region_name]) + TCut(cuts::plot_region_met_portions[region_name]);
        outfile << region_name << ": " << region << endl;
        outfile << tch_zmc->Draw("mll>>yield", TCut("totalWeight") * region, "goff") << endl;
        outfile << region_name << " yield: " << yield->Integral(0,2) << endl;
        outfile << endl;
    }
    outfile.close();
}

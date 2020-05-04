#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include "../../Main/Settings.cpp"

string oldpath = "/public/data/Photon/Samples/ReducedNtuples/";
string newpath = "/public/data/Photon/Samples/ReducedNtuples/New/";
string treename = "BaselineTree";

string selection = cuts::selections["VRDPhi"];

map<string, vector<string>> filename_sets {
    {"reweighted",
        {"data15-16_data_photon_ee.root", "data18_data_photon_ee.root", "mc16cd_SinglePhoton222_ee.root",
        "data15-16_data_photon_mm.root", "data18_data_photon_mm.root", "mc16cd_SinglePhoton222_mm.root",
        "data17_data_photon_ee.root", "mc16a_SinglePhoton222_ee.root", "mc16e_SinglePhoton222_ee.root",
        "data17_data_photon_mm.root", "mc16a_SinglePhoton222_mm.root", "mc16e_SinglePhoton222_mm.root",}},
    {"reweighted_noMC",
        {"data15-16_data_photon_ee.root", "data18_data_photon_ee.root",
        "data15-16_data_photon_mm.root", "data18_data_photon_mm.root",
        "data17_data_photon_ee.root", "data17_data_photon_mm.root",}},
    {"reduced",
        {
        "data15-16_data_bkg.root", "mc16a_singleTop.root", "mc16e_lowMassDY.root",
        "mc16a_topOther.root", "mc16cd_singleTop.root",
        "data17_data_bkg.root", "mc16a_triboson.root", "mc16cd_topOther.root", "mc16e_singleTop.root",
        "mc16a_ttbar.root", "mc16cd_triboson.root", "mc16e_topOther.root",
        "data18_data_bkg.root", "mc16cd_ttbar.root", "mc16e_triboson.root",
        "mc16a_Wjets.root", "mc16e_ttbar.root",
        "mc16a_diboson.root", "mc16a_Zjets.root", "mc16cd_Wjets.root",
        "mc16a_higgs.root", "mc16cd_diboson.root", "mc16cd_Zjets.root", "mc16e_Wjets.root",
        "mc16a_lowMassDY.root", "mc16cd_higgs.root", "mc16e_diboson.root", "mc16e_Zjets.root",
        "mc16cd_lowMassDY.root", "mc16e_higgs.root",
        }},
};
vector<string> filenames = filename_sets["reweighted_noMC"];

void skim_ntuples_with_selection() {
    TreeCreator *reducer = new TreeCreator();

    for (auto filename : filenames) {
        reducer->read(oldpath + filename, treename);

        reducer->setCut(selection);

        reducer->write(newpath + filename, treename);
    }
}

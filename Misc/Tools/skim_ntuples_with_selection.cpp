#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include "../../Main/Settings.cpp"

string oldpath = "/public/data/Photon/Samples/ReducedNtuples/";
string newpath = "/public/data/Photon/VRDPhiSamples/ReducedNtuples/";
string treename = "BaselineTree";

map<string, string> selections = {
    {"VRDPhi", "nLep_signal==2 && nLep_base==2 && trigMatch_2LTrig && is_OS && Ptll>40 && lepPt[0]>25 && lepPt[1]>25 && nJet30>=2 && jetPt[0]>30 && jetPt[1]>30 && minDPhi2JetsMet<0.4 && Ht30>250 && mt2leplsp_0>75 && (mll>81 && mll<101)"},
};
string selection = selections["VRDPhi"];

map<string, vector<string>> filename_sets {
    {"reweighted",
        {"data15-16_data_photon_ee.root", "data18_data_photon_ee.root", "mc16cd_SinglePhoton222_ee.root",
        "data15-16_data_photon_mm.root", "data18_data_photon_mm.root", "mc16cd_SinglePhoton222_mm.root",
        "data17_data_photon_ee.root", "mc16a_SinglePhoton222_ee.root", "mc16e_SinglePhoton222_ee.root",
        "data17_data_photon_mm.root", "mc16a_SinglePhoton222_mm.root", "mc16e_SinglePhoton222_mm.root",}},
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
vector<string> filenames = filename_sets["reduced"];

void skim_ntuples_with_selection() {
    TreeCreator *reducer = new TreeCreator();

    for (auto filename : filenames) {
        reducer->read(oldpath + filename, treename);

        reducer->setCut(selection);

        reducer->write(newpath + filename, treename);
    }
}

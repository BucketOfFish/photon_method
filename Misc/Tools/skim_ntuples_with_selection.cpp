#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include "../../Main/Settings.cpp"

string oldpath = "/public/data/Photon/Samples/ReweightedNtuples/";
string newpath = "/public/data/Photon/VRDPhiSamples/ReweightedNtuples/";
string treename = "BaselineTree";
string selection = "nLep_signal==2 && nLep_base==2 && trigMatch_2LTrig && is_OS && Ptll>40 && lepPt[0]>25 && lepPt[1]>25 && nJet30>=2 && jetPt[0]>30 && jetPt[1]>30 && minDPhi2JetsMet<0.4 && Ht30>250 && mt2leplsp_0>75 && (mll>81 && mll<101)";

vector<string> filenames {
    "data15-16_data_photon_ee.root", "data18_data_photon_ee.root", "mc16cd_SinglePhoton222_ee.root",
    "data15-16_data_photon_mm.root", "data18_data_photon_mm.root", "mc16cd_SinglePhoton222_mm.root",
    "data17_data_photon_ee.root", "mc16a_SinglePhoton222_ee.root", "mc16e_SinglePhoton222_ee.root",
    "data17_data_photon_mm.root", "mc16a_SinglePhoton222_mm.root", "mc16e_SinglePhoton222_mm.root",
};

void skim_ntuples_with_selection() {
    TreeCreator *reducer = new TreeCreator();

    for (auto filename : filenames) {
        reducer->read(oldpath + filename, treename);

        reducer->setCut(selection);

        reducer->write(newpath + filename, treename);
    }
}

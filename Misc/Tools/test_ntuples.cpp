#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

string path = "/eos/user/m/mazhang/PhotonMethod/v1.7/Samples/HistogramSampling/ReweightedNtuples/";

vector<string> filenames{
    //"data15-16_data_bkg.root",
    //"data15-16_data_photon.root",
    //"mc16a_diboson.root",
    //"mc16a_higgs.root",
    //"mc16a_lowMassDY.root",
    //"mc16a_SinglePhoton222.root",
    //"mc16a_singleTop.root",
    //"mc16a_topOther.root",
    //"mc16a_triboson.root",
    //"mc16a_ttbar.root",
    //"mc16a_Wjets.root",
    //"mc16a_Zjets.root",

    //"data17_data_bkg.root",
    //"data17_data_photon.root",
    //"mc16cd_diboson.root",
    //"mc16cd_higgs.root",
    //"mc16cd_lowMassDY.root",
    //"mc16cd_SinglePhoton222.root",
    //"mc16cd_singleTop.root",
    //"mc16cd_topOther.root",
    //"mc16cd_triboson.root",
    //"mc16cd_ttbar.root",
    //"mc16cd_Wjets.root",
    //"mc16cd_Zjets.root",

    //"data18_data_bkg.root",
    //"data18_data_photon.root",
    //"mc16e_diboson.root",
    //"mc16e_higgs.root",
    //"mc16e_lowMassDY.root",
    //"mc16e_SinglePhoton222.root",
    //"mc16e_singleTop.root",
    //"mc16e_topOther.root",
    //"mc16e_triboson.root",
    //"mc16e_ttbar.root",
    //"mc16e_Wjets.root",
    //"mc16e_Zjets.root",

    "data15-16_data_photon_ee.root",
    "data15-16_data_photon_mm.root",
    "data17_data_photon_ee.root",
    "data17_data_photon_mm.root",
    "data18_data_photon_ee.root",
    "data18_data_photon_mm.root",
    "mc16a_SinglePhoton222_ee.root",
    "mc16a_SinglePhoton222_mm.root",
    "mc16cd_SinglePhoton222_ee.root",
    "mc16cd_SinglePhoton222_mm.root",
    "mc16e_SinglePhoton222_ee.root",
    "mc16e_SinglePhoton222_mm.root",
};

//vector<string> branches_to_check = {"MET_sig", "met_Sign"};
vector<string> branches_to_check = {"met_Sign"};
map<string, float> branch_values;

void test_ntuples() {
    for (string filename : filenames) {
        TFile *myfile = TFile::Open((path+filename).c_str());
        TTree *mytree = (TTree*)myfile->Get("BaselineTree");

        mytree->SetBranchStatus("*", 0);
        for (auto branch : branches_to_check) {
            mytree->SetBranchStatus(branch.c_str(), 1);
            branch_values[branch] = 0.0; mytree->SetBranchAddress(branch.c_str(), &branch_values[branch]);
        }

        mytree->GetEntry(0);
        for (auto branch : branches_to_check) {
            if (branch_values[branch] == 0) {
                cout << "Check " << branch << " for " << filename << endl;
            }
        }
    }
}

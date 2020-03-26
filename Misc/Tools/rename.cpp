#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

string oldpath = "/eos/user/m/mazhang/public/PhotonMethod/";
string newpath = "/eos/user/m/mazhang/public/PhotonMethod/Renamed/";

vector<string> filenames{
    "data15-16_data_photon_ee.root",
    "data15-16_data_photon_mm.root",
    "data17_data_photon_ee.root",
    "data17_data_photon_mm.root",
    "data18_data_photon_ee.root",
    "data18_data_photon_mm.root",
};

void rename() {
    for (string filename : filenames) {
        TFile *oldfile = TFile::Open((oldpath+filename).c_str());
        TTree *oldtree = (TTree*)oldfile->Get("BaselineTree");

        TFile *newfile = new TFile((newpath+filename).c_str(), "recreate");
        TTree *newtree = oldtree->CloneTree(0);

        newtree->SetBranchStatus("MET_sig", 0);
        float met_sign; oldtree->SetBranchAddress("MET_sig", &met_sign);
        newtree->Branch("met_Sign", &met_sign, "met_Sign/F");
        //newtree->SetMakeClass(1);
        //oldtree->CopyAddresses(newtree);

        for (int i=0; i<oldtree->GetEntries(); i++)
        { // Main loop
            oldtree->GetEntry(i); // get ith entry
            newtree->Fill();
        }

        newfile->Write();
        newfile->Close();
        delete oldfile, oldtree, newfile, newtree;
    }
}

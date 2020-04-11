#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

string oldpath = "/public/data/Photon/NewSamples/HistogramSampling/ReweightedNtuples/";
string newpath = "/public/data/Photon/NewSamples/UnitTests/ReweightedNtuples/";

vector<string> filenames {
    "data15-16_data_photon_ee.root", "data18_data_photon_ee.root", "mc16cd_SinglePhoton222_ee.root",
    "data15-16_data_photon_mm.root", "data18_data_photon_mm.root", "mc16cd_SinglePhoton222_mm.root",
    "data17_data_photon_ee.root", "mc16a_SinglePhoton222_ee.root", "mc16e_SinglePhoton222_ee.root",
    "data17_data_photon_mm.root", "mc16a_SinglePhoton222_mm.root", "mc16e_SinglePhoton222_mm.root",
};

void skim_ntuples() {
    for (string filename : filenames) {
        TFile *oldfile = TFile::Open((oldpath+filename).c_str());
        TTree *oldtree = (TTree*)oldfile->Get("BaselineTree");    

        TFile *newfile = new TFile((newpath+filename).c_str(), "recreate");  
        TTree *newtree = oldtree->CloneTree(0); 

        //newtree->SetMakeClass(1);
        //oldtree->CopyAddresses(newtree);

        for (Long64_t iEntry=0; iEntry<1000000; iEntry++)
        { // Main loop
            oldtree->GetEntry(iEntry); // get ith entry
            newtree->Fill();
        }

        newfile->Write();
        newfile->Close();
        delete oldfile, oldtree, newfile, newtree;
    }
}

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

string oldpath = "/public/data/Photon/Samples/ReducedNtuples/";
string newpath = "/public/data/Photon/UnitTests/ReducedNtuples/";

vector<string> filenames {
    "data15-16_data_bkg.root",
    "data15-16_data_photon.root",
    "mc16a_higgs.root",
    "mc16a_singleTop.root",
    "mc16a_ttbar.root",
    "mc16a_lowMassDY.root",
    "mc16a_topOther.root",
    "mc16a_Wjets.root",
    "mc16a_diboson.root",
    "mc16a_SinglePhoton222.root",
    "mc16a_triboson.root",
    "mc16a_Zjets.root"
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

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

string oldpath = "/public/data/Photon/Samples/ReweightedNtuples/";
string newpath = "/public/data/Photon/Samples/ReweightedNtuples/NewBranch/";

vector<string> filenames{
    "data15-16_data_photon_ee.root",
    "data15-16_data_photon_mm.root",
    "mc16a_Vgamma_ee.root",
    "mc16a_Vgamma_mm.root",
    "data17_data_photon_ee.root",
    "data17_data_photon_mm.root",
    "mc16cd_Vgamma_ee.root",
    "mc16cd_Vgamma_mm.root",
    "data18_data_photon_ee.root",
    "data18_data_photon_mm.root",
    "mc16e_Vgamma_ee.root",
    "mc16e_Vgamma_mm.root",
};

void add_branch() {
    for (string filename : filenames) {
        TFile *oldfile = TFile::Open((oldpath+filename).c_str());
        TTree *oldtree = (TTree*)oldfile->Get("BaselineTree");

        TFile *newfile = new TFile((newpath+filename).c_str(), "recreate");
        oldtree->SetBranchStatus("dPhiPllMet", 0);
        TTree *newtree = oldtree->CloneTree(0);

        float gamma_pt; oldtree->SetBranchAddress("gamma_pt", &gamma_pt);
        float Z_eta; oldtree->SetBranchAddress("Z_eta", &Z_eta);
        float Z_phi; oldtree->SetBranchAddress("Z_phi", &Z_phi);
        float MET; oldtree->SetBranchAddress("met_Et", &MET);
        float MET_phi; oldtree->SetBranchAddress("met_Phi", &MET_phi);
        float mll = 91.188;

        float dPhiPllMet; newtree->Branch("dPhiPllMet", &dPhiPllMet, "dPhiPllMet/F");

        for (int i=0; i<oldtree->GetEntries(); i++) {
            oldtree->GetEntry(i);

            TLorentzVector z_4vec, met_4vec;
            z_4vec.SetPtEtaPhiM(gamma_pt, Z_eta, Z_phi, mll);
            met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
            
            dPhiPllMet = fabs(met_4vec.DeltaPhi(z_4vec));

            newtree->Fill();
        }

        newfile->Write();
        newfile->Close();
        delete oldfile, oldtree, newfile, newtree;
    }
}

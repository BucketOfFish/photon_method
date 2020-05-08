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

        vector<float> *lepEta = 0; oldtree->SetBranchAddress("lepEta", &lepEta);
        vector<float> *lepPhi = 0; oldtree->SetBranchAddress("lepPhi", &lepPhi);
        vector<float> *lepPt = 0; oldtree->SetBranchAddress("lepPt", &lepPt);
        vector<float> *lepM = 0; oldtree->SetBranchAddress("lepM", &lepM);
        float MET; oldtree->SetBranchAddress("met_Et", &MET);
        float MET_phi; oldtree->SetBranchAddress("met_Phi", &MET_phi);

        float dPhiPllMet; newtree->Branch("dPhiPllMet", &dPhiPllMet, "dPhiPllMet/F");

        for (int i=0; i<oldtree->GetEntries(); i++) {
            oldtree->GetEntry(i);

            TLorentzVector lep1_4vec, lep2_4vec, met_4vec;
            lep1_4vec.SetPtEtaPhiM(lepPt->at(0),lepEta->at(0),lepPhi->at(0),lepM->at(0));
            lep2_4vec.SetPtEtaPhiM(lepPt->at(1),lepEta->at(1),lepPhi->at(1),lepM->at(1));
            met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
            
            TLorentzVector boson_4vec = lep1_4vec + lep2_4vec;
            dPhiPllMet = fabs(met_4vec.DeltaPhi(boson_4vec));

            newtree->Fill();
        }

        newfile->Write();
        newfile->Close();
        delete oldfile, oldtree, newfile, newtree;
    }
}

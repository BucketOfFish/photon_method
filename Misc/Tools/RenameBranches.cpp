#include "Settings.cpp"

using namespace std;

void RenameBranches() {
    //string in_file_name = "/public/data/Photon/OldSamples/HistogramSampling/ReweightedNtuples/g_mc/mc16cd_SinglePhoton222_ee.root";
    string in_file_name = "/public/data/Photon/OldSamples/HistogramSampling/ReweightedNtuples/g_mc/mc16cd_SinglePhoton222_mm.root";
    vector<string> branches_to_copy = vector<string> {
        "PhotonConversionType",
        "met_Phi",
        "nBJet20_MV2c10_FixedCutBEff_77", "nJet30", "jetM", "jet_pT", "Ht30",
        "dPhiMetJet1", // dPhiMetJet1 is not abs()
        "gamma_pt", "totalWeight",
        "DatasetNumber", "H2PP", "H5PP", "H5PP_VR", "METOverPtISR", "METOverPtW", "METOverPtZ",
        "MJ", "MJ_VR", "MZ", "MZ_VR", "NjISR", "NjS", "PTCM", "PTCM_VR",
        "PTI", "PTISR", "PTISR_VR", "PTI_VR", "RISR", "RISR_VR", "RPT_HT5PP", "RPT_HT5PP_VR", "R_minH2P_minH3P",
        "R_minH2P_minH3P_VR", "Rjj", "Rll", "dPhiMetISR", "dPhiPjjMet", "dPhiPllMet", "dphiISRI", "dphiISRI_VR", 
        "dphiVP", "dphiVP_VR", "lept1Pt_VR", "lept2Pt_VR", "mTl3", "minDphi", "mjj",
        "mll_RJ", "mll_RJ_VR",
        "bjet_n", "jet_eta", "jet_phi", "MET_sig",
        "dPhiMetJet2", "dPhiMetJet12Min",
        "met_Et", "mll", "Ptll", "Z_eta", "Z_phi", "METt", "METl",
        "DPhi_METLepLeading",
        "lepPt", "lepEta", "lepPhi", "lepFlavor", "lepCharge",
        "nLep_base",
        "nLep_signal",
        "lepFlavor",
        "lepCharge",
        "channel",
        "reweight_Ptll",
    };
    BranchRenameOptions branches_to_rename = BranchRenameOptions {
        make_tuple("lep_eta", "lepEta"),
        make_tuple("lep_phi", "lepPhi"),
        make_tuple("jet_pT", "jetPt"),
        make_tuple("HT", "Ht30"),
    };
    BranchAddOptions branches_to_add = BranchAddOptions {
    };
    string out_file_name = "renamed.root";

    TreeCreator *renamer = new TreeCreator();
    renamer->read(in_file_name, "BaselineTree");
    renamer->setBranchesToCopy(branches_to_copy);
    renamer->setBranchesToRename(branches_to_rename);
    renamer->setBranchesToAdd(branches_to_add);
    renamer->setCut("nJet30>=0");
    renamer->setFinalCut("nJet30>=0");
    renamer->write(out_file_name, "BaselineTree");
}

#ifndef COMMON_SETTINGS
#define COMMON_SETTINGS

#include "CommonLibraries.C"
#include "CommonCuts.C"
#include "CommonBins.C"
#include "CommonFunctions.C"

vector<string> histFitterBranches {"DatasetNumber/I", "Etall/F", "H2PP/D", "H5PP/D", "H5PP_VR/D",
    "METOverPtISR/F", "METOverPtW/F", "METOverPtZ/F", "MJ/D", "MJ_VR/D", "MZ/D", "MZ_VR/D", "NjISR/D",
    "NjS/D", "PTCM/D", "PTCM_VR/D", "PTI/D", "PTISR/D", "PTISR_VR/D", "PTI_VR/D", "RISR/D", "RISR_VR/D",
    "RPT_HT5PP/D", "RPT_HT5PP_VR/D", "R_minH2P_minH3P/D", "R_minH2P_minH3P_VR/D", "Rjj/F", "Rll/F",
    "dPhiMetISR/F", "dPhiMetJet1/F", "dPhiMetJet2/F", "dPhiMetJet12Min/F", "dPhiPjjMet/F", "dPhiPllMet/F",
    "dphiISRI/D", "dphiISRI_VR/D", "dphiVP/D", "dphiVP_VR/D", "lept1Pt_VR/D", "lept2Pt_VR/D", "mTl3/D",
    "met_Sign/F", "minDphi/D", "mll_RJ/D", "mll_RJ_VR/D", "mt2leplsp_0/F", "nJet20/I", "mjj/F"};

using BranchRenameOptions = vector<tuple<string, string>>;
using BranchAddOptions = vector<tuple<string, string>>;

struct GlobalOptions {
    bool is_photon; // vs. bkg
    bool is_data; // vs. MC
    string period;
    string sampleID;

    string photon_mc_path;
    string photon_data_path;
    string bkg_mc_path;
    string bkg_data_path;

    string my_samples_folder;
    string sampling_method;
    string reduction_folder;
    string smearing_folder;
    string reweighting_folder;
    string plots_folder;

    string out_tree_name;
};

struct ReductionOptions {
    string in_file_name;
    string in_tree_name;
    string out_file_name;
    string out_tree_name;

    vector<string> branches_to_copy;
    BranchRenameOptions branches_to_rename;
    BranchAddOptions branches_to_add;

    string cut;

    bool unit_testing;
};

#endif

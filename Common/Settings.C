#ifndef COMMON_SETTINGS
#define COMMON_SETTINGS

#include "CommonLibraries.C"
#include "CommonCuts.C"
#include "CommonBins.C"
#include "CommonFunctions.C"

//std::string sample_folder = "/eos/user/m/mazhang/PhotonMethod/v1.7/Samples/";
//std::string sampling_method = "HistogramSampling";
std::string sample_folder = "/public/data/Photon/";
std::string sampling_method = "HistogramSampling";

std::string ntuple_path =  sample_folder + "Ntuples/";
std::string smearing_path = sample_folder + "/" + sampling_method + "/SmearedNtuples/";
std::string reweighting_path = sample_folder + "/" + sampling_method + "/ReweightedNtuples/";
std::string plots_path = sample_folder + "/" + sampling_method + "/Plots/";

vector<string> histFitterBranches {"DatasetNumber/I", "Etall/F", "H2PP/D", "H5PP/D", "H5PP_VR/D",
    "METOverPtISR/F", "METOverPtW/F", "METOverPtZ/F", "MJ/D", "MJ_VR/D", "MZ/D", "MZ_VR/D", "NjISR/D",
    "NjS/D", "PTCM/D", "PTCM_VR/D", "PTI/D", "PTISR/D", "PTISR_VR/D", "PTI_VR/D", "RISR/D", "RISR_VR/D",
    "RPT_HT5PP/D", "RPT_HT5PP_VR/D", "R_minH2P_minH3P/D", "R_minH2P_minH3P_VR/D", "Rjj/F", "Rll/F",
    "dPhiMetISR/F", "dPhiMetJet1/F", "dPhiMetJet2/F", "dPhiMetJet12Min/F", "dPhiPjjMet/F", "dPhiPllMet/F",
    "dphiISRI/D", "dphiISRI_VR/D", "dphiVP/D", "dphiVP_VR/D", "lept1Pt_VR/D", "lept2Pt_VR/D", "mTl3/D",
    "met_Sign/F", "minDphi/D", "mll_RJ/D", "mll_RJ_VR/D", "mt2leplsp_0/F", "nJet20/I", "mjj/F"};

using BranchRenameOptions = vector<tuple<string, string>>;
using BranchAddOptions = vector<tuple<string, string>>;

struct Options {
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

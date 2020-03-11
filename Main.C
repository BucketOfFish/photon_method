#include "Common/Settings.C"
#include "ReduceNtuples.C"

using namespace std;

void Main() {
    Options options;

    options.in_file_name = "/public/data/Photon/Ntuples/bkg_data/data15-16_bkg.root";
    options.in_tree_name = "BaselineTree";
    options.out_file_name = "/public/data/Photon/SkimmedSamples/data15-16_bkg.root";
    options.out_tree_name = "BaselineTree";

    options.branches_to_copy = vector<string> {
        "lepIsoFCTight", "lepIsPR", "nLep_signal", "nLep_base",
        "lepEta", "lepPhi", "lepM", "lepFlavor", "lepCharge", "lepPt",
        "PhotonConversionType",
        "met_Phi", "met_Et",
        "mll",
        "Ptll",
        "nBJet20_MV2c10_FixedCutBEff_77", "nJet30", "jetM", "jetPt", "Ht30",
        "minDPhi2JetsMet",
        "trigMatch_2LTrig", "trigMatch_2LTrigOR", "globalDiLepTrigSF",
        "getWeight", "eventWeight", "leptonWeight", "jvtWeight", "bTagWeight", "pileupWeight",
        "RunNumber", "RandomRunNumber",
    };
    options.branches_to_rename = BranchRenameOptions {
        make_tuple("PhotonPt", "gamma_pt"),
        make_tuple("PhotonEta", "gamma_eta"),
        make_tuple("PhotonPhi", "gamma_phi"),
        make_tuple("nBJet30_MV2c10_FixedCutBEff_77", "bjet_n"),
        make_tuple("jetEta", "jet_eta"),
        make_tuple("jetPhi", "jet_phi"),
        make_tuple("met_Sign", "MET_sig"),
    };

    options.cut = "met_Et>300";

    ReduceNtuples(options);
}

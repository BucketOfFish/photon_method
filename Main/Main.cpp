#include "Settings.cpp"
#include "ReduceNtuples.cpp"
#include "SmearPhotons.cpp"
//#include "InProgress/SmearPhotons.cpp"
#include "ReweightPhotons.cpp"
#include "MakePlots.cpp"

using namespace std;
using rvecf = ROOT::VecOps::RVec<float>;

//----------------------------------
// REDUCTION AND SMEARING FUNCTIONS
//----------------------------------

void initFillingFunctions() {
    unordered_map<string, string> filling_functions;

    filling_functions["getDPhiMetJet"] =
        "vector<float> getDPhiMetJet(rvecf jet_pT, rvecf jet_eta, rvecf jet_phi, int jet_n, double MET, double MET_phi) {"
            "TLorentzVector jet_4vec, met_4vec;"
            "met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);"
            ""
            "vector<float> dPhiMetJet;"
            "for (int i=0; i<jet_n; i++) {"
            "   jet_4vec.SetPtEtaPhiM(jet_pT[i],jet_eta[i],jet_phi[i],0);"
            "   dPhiMetJet.push_back(fabs(met_4vec.DeltaPhi(jet_4vec)));"
            "}"
            "return dPhiMetJet;"
        "}";
    filling_functions["getPhotonDataWeight"] =
        "double getPhotonDataWeight(float trigMatch_HLT_g15_loose_L1EM7, float trigPrescale_HLT_g15_loose_L1EM7,"
        "float trigMatch_HLT_g25_loose_L1EM15, float trigPrescale_HLT_g25_loose_L1EM15,"
        "float trigMatch_HLT_g35_loose_L1EM15, float trigPrescale_HLT_g35_loose_L1EM15,"
        "float trigMatch_HLT_g40_loose_L1EM15, float trigPrescale_HLT_g40_loose_L1EM15,"
        "float trigMatch_HLT_g45_loose_L1EM15, float trigPrescale_HLT_g45_loose_L1EM15,"
        "float trigMatch_HLT_g50_loose_L1EM15, float trigPrescale_HLT_g50_loose_L1EM15,"
        "float trigMatch_HLT_g60_loose, float trigPrescale_HLT_g60_loose,"
        "float trigMatch_HLT_g70_loose, float trigPrescale_HLT_g70_loose,"
        "float trigMatch_HLT_g80_loose, float trigPrescale_HLT_g80_loose,"
        "float trigMatch_HLT_g100_loose, float trigPrescale_HLT_g100_loose,"
        "float trigMatch_HLT_g140_loose, float trigPrescale_HLT_g140_loose, float gamma_pt) {"
            "double totalWeight = 0;"
            ""
            "if (trigMatch_HLT_g15_loose_L1EM7==1 && gamma_pt>(15) && gamma_pt<(25+5)) totalWeight = trigPrescale_HLT_g15_loose_L1EM7;"
            "if (trigMatch_HLT_g25_loose_L1EM15==1 && gamma_pt>(25+5) && gamma_pt<(35+5)) totalWeight = trigPrescale_HLT_g25_loose_L1EM15;"
            "if (trigMatch_HLT_g35_loose_L1EM15==1 && gamma_pt>(35+5) && gamma_pt<(40+5)) totalWeight = trigPrescale_HLT_g35_loose_L1EM15;"
            "if (trigMatch_HLT_g40_loose_L1EM15==1 && gamma_pt>(40+5) && gamma_pt<(45+5)) totalWeight = trigPrescale_HLT_g40_loose_L1EM15;"
            "if (trigMatch_HLT_g45_loose_L1EM15==1 && gamma_pt>(45+5) && gamma_pt<(50+5)) totalWeight = trigPrescale_HLT_g45_loose_L1EM15;"
            "if (trigMatch_HLT_g50_loose_L1EM15==1 && gamma_pt>(50+5) && gamma_pt<(60+5)) totalWeight = trigPrescale_HLT_g50_loose_L1EM15;"
            "if (trigMatch_HLT_g60_loose==1 && gamma_pt>(60+5) && gamma_pt<(70+5)) totalWeight = trigPrescale_HLT_g60_loose;"
            "if (trigMatch_HLT_g70_loose==1 && gamma_pt>(70+5) && gamma_pt<(80+5)) totalWeight = trigPrescale_HLT_g70_loose;"
            "if (trigMatch_HLT_g80_loose==1 && gamma_pt>(80+5) && gamma_pt<(100+5)) totalWeight = trigPrescale_HLT_g80_loose;"
            "if (trigMatch_HLT_g100_loose==1 && gamma_pt>(100+5) && gamma_pt<(140+5)) totalWeight = trigPrescale_HLT_g100_loose;"
            "if (trigMatch_HLT_g140_loose==1 && gamma_pt>(140+5)) totalWeight = trigPrescale_HLT_g140_loose;"
            ""
            "if (totalWeight > 100000000000) totalWeight=0;" //--- fix for large photon sample spikes
            ""
            "return totalWeight;"
        "}";
    filling_functions["getPhotonMCWeight"] =
        "double getPhotonMCWeight(float lumi, float genWeight, float eventWeight, float jvtWeight, float bTagWeight, float pileupWeight) {"
            "double totalWeight = lumi*genWeight*eventWeight*jvtWeight*bTagWeight*pileupWeight;"
            ""
            "if (totalWeight > 100000000000) totalWeight=0;" //--- fix for large photon sample spikes
            ""
            "return totalWeight;"
        "}";
    filling_functions["getZEta"] =
        "float getZEta(rvecf lep_pT, rvecf lep_eta, rvecf lep_phi) {"
            "TLorentzVector l0_4vec, l1_4vec;"
            "l0_4vec.SetPtEtaPhiM(lep_pT[0],lep_eta[0],lep_phi[0],0);"
            "l1_4vec.SetPtEtaPhiM(lep_pT[1],lep_eta[1],lep_phi[1],0);"
            ""
            "float Z_eta = (l0_4vec+l1_4vec).Eta();"
            "return Z_eta;"
        "}";
    filling_functions["getZPhi"] =
        "float getZPhi(rvecf lep_pT, rvecf lep_eta, rvecf lep_phi) {"
            "TLorentzVector l0_4vec, l1_4vec;"
            "l0_4vec.SetPtEtaPhiM(lep_pT[0],lep_eta[0],lep_phi[0],0);"
            "l1_4vec.SetPtEtaPhiM(lep_pT[1],lep_eta[1],lep_phi[1],0);"
            ""
            "float Z_phi = (l0_4vec+l1_4vec).Phi();"
            "return Z_phi;"
        "}";
    filling_functions["getZCMLepTheta"] =
        "vector<float> getZCMLepTheta(rvecf lep_pT, rvecf lep_eta, rvecf lep_phi, float Z_pt, float Z_eta, float Z_phi) {"
            "TLorentzVector l0_cm_4vec, l1_cm_4vec;"
            "l0_cm_4vec.SetPtEtaPhiM(lep_pT[0],lep_eta[0],lep_phi[0],0);"
            "l1_cm_4vec.SetPtEtaPhiM(lep_pT[1],lep_eta[1],lep_phi[1],0);"
            ""
            "TLorentzVector z_4vec;"
            "float Z_m = 91.1876;"
            "z_4vec.SetPtEtaPhiM(Z_pt,Z_eta,Z_phi,Z_m);"
            "TVector3 boost_vec(0, 0, -z_4vec.BoostVector().Mag());"
            ""
            "l0_cm_4vec.RotateZ(-z_4vec.Phi());"
            "l0_cm_4vec.RotateY(-z_4vec.Theta());"
            "l0_cm_4vec.Boost(boost_vec);"
            "l1_cm_4vec.RotateZ(-z_4vec.Phi());"
            "l1_cm_4vec.RotateY(-z_4vec.Theta());"
            "l1_cm_4vec.Boost(boost_vec);"
            ""
            "vector<float> Z_cm_lep_theta;"
            "Z_cm_lep_theta.push_back(l0_cm_4vec.Theta());"
            "Z_cm_lep_theta.push_back(l1_cm_4vec.Theta());"
            "return Z_cm_lep_theta;"
        "}";
    filling_functions["getDR2Lep"] =
        "float getDR2Lep(rvecf lep_pT, rvecf lep_eta, rvecf lep_phi) {"
            "TLorentzVector l0_4vec, l1_4vec;"
            "l0_4vec.SetPtEtaPhiM(lep_pT[0],lep_eta[0],lep_phi[0],0);"
            "l1_4vec.SetPtEtaPhiM(lep_pT[1],lep_eta[1],lep_phi[1],0);"
            ""
            "return l0_4vec.DeltaR(l1_4vec);"
        "}";
    filling_functions["getDPhi2Lep"] =
        "float getDPhi2Lep(rvecf lep_pT, rvecf lep_eta, rvecf lep_phi) {"
            "TLorentzVector l0_4vec, l1_4vec;"
            "l0_4vec.SetPtEtaPhiM(lep_pT[0],lep_eta[0],lep_phi[0],0);"
            "l1_4vec.SetPtEtaPhiM(lep_pT[1],lep_eta[1],lep_phi[1],0);"
            ""
            "return fabs(l0_4vec.DeltaPhi(l1_4vec));"
        "}";
    filling_functions["getDPhiMETZPhoton"] =
        "float getDPhiMETZPhoton(float Z_pt, float Z_eta, float Z_phi, float MET, float MET_phi) {"
            "TLorentzVector z_4vec;"
            "float Z_m = 91.1876;"
            "z_4vec.SetPtEtaPhiM(Z_pt,Z_eta,Z_phi,Z_m);"
            ""
            "TLorentzVector met_4vec;"
            "met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);"
            ""
            "return fabs(met_4vec.DeltaPhi(z_4vec));"
        "}";
    filling_functions["getDPhiMETLepLeading"] =
        "float getDPhiMETLepLeading(rvecf lep_pT, rvecf lep_eta, rvecf lep_phi, float MET, float MET_phi) {"
            "TLorentzVector l0_4vec;"
            "l0_4vec.SetPtEtaPhiM(lep_pT[0],lep_eta[0],lep_phi[0],0);"
            ""
            "TLorentzVector met_4vec;"
            "met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);"
            ""
            "return fabs(met_4vec.DeltaPhi(l0_4vec));"
        "}";
    filling_functions["getDPhiMETLepSecond"] =
        "float getDPhiMETLepSecond(rvecf lep_pT, rvecf lep_eta, rvecf lep_phi, float MET, float MET_phi) {"
            "TLorentzVector l1_4vec;"
            "l1_4vec.SetPtEtaPhiM(lep_pT[1],lep_eta[1],lep_phi[1],0);"
            ""
            "TLorentzVector met_4vec;"
            "met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);"
            ""
            "return fabs(met_4vec.DeltaPhi(l1_4vec));"
        "}";

    for (auto const& [key, val] : filling_functions)
        gInterpreter->Declare(val.c_str());
}

void ReductionStep(Options options) {
    string in_folder;
    if (options.is_data)
        if (options.is_photon) in_folder = options.photon_data_path;
        else in_folder = options.bkg_data_path;
    else
        if (options.is_photon) in_folder = options.photon_mc_path;
        else in_folder = options.bkg_mc_path;

    if (options.is_data) options.in_file_name = Form("%s/%s_merged_processed.root", in_folder.c_str(), options.period.c_str()); 
    else options.in_file_name = Form("%s%s/%s_merged_processed.root", in_folder.c_str(), options.period.c_str(), options.sampleID.c_str()); 

    if (options.is_data) {
       if (options.is_photon) options.in_tree_name = options.period;
       else options.in_tree_name = "data";
    }
    else options.in_tree_name = options.sampleID + "_NoSys";

    if (options.is_data) {
       if (options.is_photon) options.out_file_name = options.reduction_folder + "/" + options.period.c_str() + "_data_photon.root";
       else options.out_file_name = options.reduction_folder + "/" + options.period.c_str() + "_data_bkg.root";
    }
    else options.out_file_name = options.reduction_folder + "/" + options.period.c_str() + "_" + options.sampleID.c_str() + ".root";

    //--- branches to copy from old tree to new tree
    options.branches_to_copy = vector<string> {
        "PhotonConversionType",
        "met_Phi",
        "nBJet20_MV2c10_FixedCutBEff_77", "nJet30", "jetM", "jetPt", "Ht30",
        "dPhiMetJet1", "minDPhi2JetsMet", // dPhiMetJet1 is not abs()
        "genWeight", "eventWeight", "leptonWeight", "jvtWeight", "bTagWeight", "pileupWeight", "globalDiLepTrigSF",
        "RunNumber", "RandomRunNumber",
        "DatasetNumber", "H2PP", "H5PP", "H5PP_VR", "METOverPtISR", "METOverPtW", "METOverPtZ",
        "MJ", "MJ_VR", "MZ", "MZ_VR", "NjISR", "NjS", "PTCM", "PTCM_VR",
        "PTI", "PTISR", "PTISR_VR", "PTI_VR", "RISR", "RISR_VR", "RPT_HT5PP", "RPT_HT5PP_VR", "R_minH2P_minH3P",
        "R_minH2P_minH3P_VR", "Rjj", "Rll", "dPhiMetISR", "dPhiPjjMet", "dPhiPllMet", "dphiISRI", "dphiISRI_VR", 
        "dphiVP", "dphiVP_VR", "lept1Pt_VR", "lept2Pt_VR", "mTl3", "minDphi", "mjj",
        "mll_RJ", "mll_RJ_VR", "nJet20", "met_Sign",
        "mjj_minDPhiZMET", "mbb", "PtISR",
    };

    //--- branches to rename and copy
    options.branches_to_rename = BranchRenameOptions {
        make_tuple("nBJet30_MV2c10_FixedCutBEff_77", "bjet_n"),
        make_tuple("jetEta", "jet_eta"),
        make_tuple("jetPhi", "jet_phi"),
    };

    //--- new branches to add
    options.branches_to_add = BranchAddOptions {
        make_tuple("dPhiMetJet", "getDPhiMetJet(jetPt, jetEta, jetPhi, nJet30, met_Et, met_Phi)"),
        make_tuple("dPhiMetJet2", "dPhiMetJet[1]"),
        make_tuple("dPhiMetJet12Min", "std::min(dPhiMetJet[0], dPhiMetJet[1])"),
        make_tuple("lumi", to_string(GetLumi(options.period))), //if (TString(options.sampleID).Contains("Vg")) lumi *= -1;
    };
    
    //--- photon/bkg specific branches
    vector<string> additional_copy;
    BranchRenameOptions additional_rename;
    BranchAddOptions additional_add;

    if (options.is_photon) {
        additional_rename = BranchRenameOptions {
            make_tuple("met_Et", "met_Et_unsmeared"),
            make_tuple("PhotonPt", "gamma_pt"),
            make_tuple("PhotonEta", "gamma_eta"),
            make_tuple("PhotonPhi", "gamma_phi"),
        };
        additional_add = BranchAddOptions {
            make_tuple("METt_unsmeared", "met_Et*sin(met_Phi-PhotonPt)"),
            make_tuple("METl_unsmeared", "met_Et*cos(met_Phi-PhotonPhi)"),
            make_tuple("trigMatch_2LTrig", "1"),
            make_tuple("trigMatch_2LTrigOR", "1"),
        };
        if (options.is_data)
            additional_add.push_back(make_tuple("totalWeight", "getPhotonDataWeight(trigMatch_HLT_g15_loose_L1EM7,"
                "trigPrescale_HLT_g15_loose_L1EM7, trigMatch_HLT_g25_loose_L1EM15, trigPrescale_HLT_g25_loose_L1EM15, trigMatch_HLT_g35_loose_L1EM15,"
                "trigPrescale_HLT_g35_loose_L1EM15, trigMatch_HLT_g40_loose_L1EM15, trigPrescale_HLT_g40_loose_L1EM15, trigMatch_HLT_g45_loose_L1EM15,"
                "trigPrescale_HLT_g45_loose_L1EM15, trigMatch_HLT_g50_loose_L1EM15, trigPrescale_HLT_g50_loose_L1EM15, trigMatch_HLT_g60_loose,"
                "trigPrescale_HLT_g60_loose, trigMatch_HLT_g70_loose, trigPrescale_HLT_g70_loose, trigMatch_HLT_g80_loose, trigPrescale_HLT_g80_loose,"
                "trigMatch_HLT_g100_loose, trigPrescale_HLT_g100_loose, trigMatch_HLT_g140_loose, trigPrescale_HLT_g140_loose, PhotonPt)"));
        else {
            additional_add.push_back(make_tuple("totalWeight", "getPhotonMCWeight(lumi, genWeight, eventWeight, jvtWeight, bTagWeight, pileupWeight)"));
        }
    }
    else {
        additional_copy = vector<string> {
            "lepIsoFCTight", "lepIsPR",
            "lepEta", "lepPhi", "lepM", "lepFlavor", "lepCharge", "lepPt",
            "channel",
            "trigMatch_2LTrig", "trigMatch_2LTrigOR", "nLep_signal", "nLep_base",
            "met_Et", "mt2leplsp_0",
            "mll", "Ptll",
        };
        additional_add = BranchAddOptions {
            make_tuple("is_OS", "lepCharge[0]!=lepCharge[1]"),
            make_tuple("Z_eta", "getZEta(lepPt, lepEta, lepPhi)"),
            make_tuple("Z_phi", "getZPhi(lepPt, lepEta, lepPhi)"),
            make_tuple("METt", "met_Et*sin(met_Phi-Z_phi)"),
            make_tuple("METl", "met_Et*cos(met_Phi-Z_phi)"),
            make_tuple("Z_cm_lep_theta", "getZCMLepTheta(lepPt, lepEta, lepPhi, Ptll, Z_eta, Z_phi)"),
            make_tuple("DR_2Lep", "getDR2Lep(lepPt, lepEta, lepPhi)"),
            make_tuple("DPhi_2Lep", "getDPhi2Lep(lepPt, lepEta, lepPhi)"),
            make_tuple("DPhi_METZPhoton", "getDPhiMETZPhoton(Ptll, Z_eta, Z_phi, met_Et, met_Phi)"),
            make_tuple("DPhi_METLepLeading", "getDPhiMETLepLeading(lepPt, lepEta, lepPhi, met_Et, met_Phi)"),
            make_tuple("DPhi_METLepSecond", "getDPhiMETLepSecond(lepPt, lepEta, lepPhi, met_Et, met_Phi)"),
            make_tuple("DPhi_METLepMin", "std::min(DPhi_METLepLeading, DPhi_METLepSecond)"),
        };
        if (options.is_data)
            additional_add.push_back(make_tuple("totalWeight", "1.0"));
        else
            additional_add.push_back(make_tuple("totalWeight", "lumi*genWeight*eventWeight*leptonWeight*jvtWeight*bTagWeight*pileupWeight*FFWeight"));
    }

    options.branches_to_copy.insert(options.branches_to_copy.end(), additional_copy.begin(), additional_copy.end());
    options.branches_to_rename.insert(options.branches_to_rename.end(), additional_rename.begin(), additional_rename.end());
    options.branches_to_add.insert(options.branches_to_add.end(), additional_add.begin(), additional_add.end());

    //--- set selection cut
    if (options.is_photon)
        options.cut = cuts::selections["photon_baseline_ntuples"] + "totalWeight < 1000000";
    else
        options.cut = cuts::selections["bkg_baseline"] + "totalWeight < 1000000";
    options.final_cut = "totalWeight!=0";

    //--- make reduced ntuples
    ReduceNtuples(options);
}

void initSmearingFunctions() {
    unordered_map<string, string> smearing_functions;

    smearing_functions["getLepFlavors"] =
        "vector<int> getLepFlavors(int flavor) {"
            "vector<int> lepFlavors{flavor, flavor};"
            "return lepFlavors;"
        "}";
    smearing_functions["getLepCharges"] =
        "vector<int> getLepCharges() {"
            "int charge = (std::rand()%2)*2-1;"
            "vector<int> lepCharges{charge, -charge};"
            "return lepCharges;"
        "}";

    for (auto const& [key, val] : smearing_functions)
        gInterpreter->Declare(val.c_str());
}

//---------------
// MAIN FUNCTION
//---------------

void Main() {
    //ROOT::EnableImplicitMT(); // enable parallelization to speed up RDataFrame
    gErrorIgnoreLevel = kWarning; // turn off info dumps
    bins::init_binning_histograms(); // prepare histograms
    TH1::SetDefaultSumw2(); // histogram summing option

    Options options;

    string SUSY_folder = "/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.7/";
    options.photon_mc_path = SUSY_folder + "JETM4/JETM4_";
    options.photon_data_path = SUSY_folder + "JETM4/JETM4_Data/";
    options.bkg_mc_path = SUSY_folder + "SUSY2/SUSY2_Bkgs_";
    options.bkg_data_path = SUSY_folder + "SUSY2/SUSY2_Data/SUSY2_Data_v1.7/merged/";
    //options.bkg_mc_path = "/eos/atlas/user/l/longjon/Ntuples/2L2J_skims/skim_slim_v1.7/2LTrigOR_nBaseLep25-ge-2_nJet30-ge-2_metEt-gt-200_Ht30-gt-200-if-mll-gt-81/SUSY2_Bkgs_"
    //options.bkg_data_path = "/eos/atlas/user/l/longjon/Ntuples/2L2J_skims/skim_slim_v1.7/2LTrigOR_nBaseLep25-ge-2_nJet30-ge-2_metEt-gt-200_Ht30-gt-200-if-mll-gt-81/SUSY2_Data/"

    options.my_samples_folder = "/public/data/Photon/Samples/";
    //options.my_samples_folder = "/eos/user/m/mazhang/PhotonMethod/v1.7/Samples/";

    options.reduction_folder = options.my_samples_folder + "ReducedNtuples/";
    options.smearing_folder = options.my_samples_folder + "SmearedNtuples/";
    //options.reweighting_folder = options.my_samples_folder + "ReweightedNtuples/";
    options.reweighting_folder = options.my_samples_folder + "SkimmedNtuples/StrongPreselectionNtuples/";
    options.plots_folder = options.my_samples_folder + "Plots/ZMC/";

    options.unit_test_folder = "/public/data/Photon/UnitTestSamples/";
    //options.unit_test_folder = "/eos/user/m/mazhang/PhotonMethod/v1.7/UnitTestSamples/";

    options.tree_name = "BaselineTree";

    options.unit_testing = true;
    bool do_reduction = false;
    bool do_smearing = false;
    bool do_reweighting = false;
    bool do_plotting = true;

    //--- unit testing
    if (options.unit_testing) {
        if (do_reduction) ReductionStep(options);
        if (do_smearing) SmearPhotons(options);
        if (do_reweighting) ReweightPhotons(options);
        if (do_plotting) MakePlots(options);
        return;
    }

    //--- reduce ntuples
    if (do_reduction) {
        initFillingFunctions(); // functions used for adding new branches

        vector<bool> is_datas{true, false};
        vector<string> periods{"data15-16", "data17", "data18"};
        for (auto is_data : is_datas) {
            for (auto period : periods) {
                vector<string> sampleIDs{"data", "photon"};
                if (!is_data) sampleIDs = vector<string>{"SinglePhoton222", "Vgamma", "Zjets", "ttbar", "diboson", "higgs",
                    "lowMassDY", "singleTop", "topOther", "triboson", "Wjets"};

                for (auto sampleID : sampleIDs) {
                    if (sampleID == "photon") {
                        options.sampleID = "data";
                        options.is_photon = true;
                    }
                    else {
                        options.sampleID = sampleID;
                        if (sampleID == "SinglePhoton222") options.is_photon = true;
                        else if (sampleID == "Vgamma") options.is_photon = true;
                        else options.is_photon = false;
                    }
                    options.is_data = is_data;
                    options.period = period;
                    if (!options.is_data) {
                        if (options.period == "data15-16") options.period = "mc16a";
                        else if (options.period == "data17") options.period = "mc16cd";
                        else if (options.period == "data18") options.period = "mc16e";
                    }

                    ReductionStep(options);
                }
            }
        }
    }

    //--- smear photons
    if (do_smearing) {
        initSmearingFunctions();
        options.turn_off_shifting_and_smearing = false;

        vector<string> periods{"data15-16", "data17", "data18"};
        for (auto period : periods) {
            options.period = period;
            options.data_period = DataPeriod(options.period);
            options.mc_period = getMCPeriod(options.period);

            //vector<bool> is_datas{true, false};
            vector<string> channels{"ee", "mm"};
            vector<bool> is_datas{true};
            //vector<string> channels{"mm"};
            for (auto is_data : is_datas) {
                for (auto channel : channels) {
                    options.channel = channel;
                    options.is_data = is_data;
                    SmearPhotons(options);
                }
            }
        }
    }

    //--- reweight photons
    if (do_reweighting) {
        //options.reweight_vars = {"Ptll", "nBJet20_MV2c10_FixedCutBEff_77", "nJet30", "Ht30", "Ptll__Ht30",
            //"Ptll__Zwindow", "nBJet20_MV2c10_FixedCutBEff_77__Zwindow", "nJet30__Zwindow", "Ht30__Zwindow",
            //"Ptll__Ht30__Zwindow"};
        options.reweight_vars = {"Ptll", "Ht30", "Ptll__Ht30",
            "Ptll__Zwindow", "Ht30__Zwindow", "Ptll__Ht30__Zwindow"};

        vector<string> periods{"data15-16", "data17", "data18"};
        for (auto period : periods) {
            options.period = period;
            options.data_period = DataPeriod(options.period);
            options.mc_period = getMCPeriod(options.period);

            //vector<bool> is_datas{true, false};
            vector<bool> is_datas{true};
            vector<string> channels{"ee", "mm"};
            for (auto is_data : is_datas) {
                for (auto channel : channels) {
                    options.channel = channel;
                    options.is_data = is_data;
                    if (options.is_data) options.processes = {"data", "tt", "vv", "photon"};
                    else options.processes = {"zjets", "photon"};
                    ReweightPhotons(options);
                }
            }
        }
    }

    //--- make plots
    if (do_plotting) {
        options.reweight_branch = "reweight_Ptll";

        //options.plot_regions = vector<string>{"SRC", "SRLow2", "SRMed2", "SRHigh2", "SRLowZ4", "SRMedZ4", "SRHighZ4",
                                        //"VRC", "VRLow2", "VRMed2", "VRHigh2", "VRLowZ4", "VRMedZ4", "VRHighZ4"};
        options.plot_regions = vector<string>{"VRZjets_noZwindow", "VRZjets_noZwindow_noMETcut"};
        options.plot_features = vector<string>{"met_Et", "METl", "METt", "mt2leplsp_0", "dPhiMetJet12Min", "dPhiPllMet"};
        //options.plot_features = vector<string>{"met_Et"};
        //options.plot_channels = vector<string>{"ee", "mm", "SF"};
        options.plot_channels = vector<string>{"SF"};

        //options.additional_plot_cut = "Ht30>250";

        options.make_diagnostic_plots = true;
        //options.diagnostic_plots = vector<string>{"lepEta", "METl", "mll"};

        //options.processes = {"data_bkg", "photon", "Zjets", "ttbar", "diboson", "higgs", "singleTop", "topOther",
                             //"Wjets", "triboson"};
        options.processes = {"data_bkg", "photon", "Zjets", "ttbar", "diboson"};

        options.blinded = true;
        options.print_photon_yield_only = false;
        options.do_vgamma_subtraction = false;

        options.plot_reweighted_photons = false;
        options.plot_unreweighted_photons = false;
        options.plot_zmc = true;

        options.scale_zmc = true;
        options.scaling_method = "MET";
        //options.scaling_method = "dPhiMetJet12Min";

        //vector<string> periods{"data15-16", "data17", "data18"};
        vector<string> periods{"all"};
        //vector<bool> is_datas{true, false};
        vector<bool> is_datas{true};
        for (auto period : periods) {
            options.period = period;
            options.data_period = DataPeriod(options.period);
            options.mc_period = getMCPeriod(options.period);
            for (auto is_data : is_datas) {
                options.is_data = is_data;
                //for (auto reweight_branch : {"reweight_Ptll", "reweight_Ptll__Ht30"}) {
                    //options.reweight_branch = reweight_branch;
                    MakePlots(options);
                //}
            }
        }
    }
}

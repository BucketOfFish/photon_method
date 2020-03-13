#include "Common/Settings.C"
#include "ReduceNtuples.C"

using namespace std;
using rvecf = ROOT::VecOps::RVec<float>;

//------------------
// NTUPLE REDUCTION
//------------------

unordered_map<string, string> getFillingFunctions(bool is_data) {
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
    if (is_data) {
        filling_functions["getPhotonWeight"] =
            "float getPhotonWeight(float trigMatch_HLT_g15_loose_L1EM7, float trigPrescale_HLT_g15_loose_L1EM7,"
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
                "float totalWeight = 0;"
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
    }
    else {
        filling_functions["getPhotonWeight"] =
            "float getPhotonWeight(float lumi, float genWeight, float eventWeight, float jvtWeight, float bTagWeight, float pileupWeight) {"
                "float totalWeight = lumi*genWeight*eventWeight*jvtWeight*bTagWeight*pileupWeight;"
                ""
                "if (totalWeight > 100000000000) totalWeight=0;" //--- fix for large photon sample spikes
                ""
                "return totalWeight;"
            "}";
    }
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
            "float Z_m = 91.1876"
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
            "Z_cm_lep_theta->push_back(l0_cm_4vec.Theta());"
            "Z_cm_lep_theta->push_back(l1_cm_4vec.Theta());"
            "return Z_cm_lep_theta;"
        "}";
    filling_functions["getDR2Lep"] =
        "vector<float> getDR2Lep(rvecf lep_pT, rvecf lep_eta, rvecf lep_phi) {"
            "TLorentzVector l0_4vec, l1_4vec;"
            "l0_4vec.SetPtEtaPhiM(lep_pT[0],lep_eta[0],lep_phi[0],0);"
            "l1_4vec.SetPtEtaPhiM(lep_pT[1],lep_eta[1],lep_phi[1],0);"
            ""
            "return lep0_4vec.DeltaR(lep1_4vec);"
        "}";
    filling_functions["getDPhi2Lep"] =
        "vector<float> getDPhi2Lep(rvecf lep_pT, rvecf lep_eta, rvecf lep_phi) {"
            "TLorentzVector l0_4vec, l1_4vec;"
            "l0_4vec.SetPtEtaPhiM(lep_pT[0],lep_eta[0],lep_phi[0],0);"
            "l1_4vec.SetPtEtaPhiM(lep_pT[1],lep_eta[1],lep_phi[1],0);"
            ""
            "return fabs(lep0_4vec.DeltaPhi(lep1_4vec));"
        "}";
    filling_functions["getDPhiMETZPhoton"] =
        "vector<float> getDPhiMETZPhoton(float Z_pt, float Z_eta, float Z_phi, float MET, float MET_phi) {"
            "TLorentzVector z_4vec;"
            "float Z_m = 91.1876"
            "z_4vec.SetPtEtaPhiM(Z_pt,Z_eta,Z_phi,Z_m);"
            ""
            "TLorentzVector met_4vec;"
            "met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);"
            ""
            "return fabs(met_4vec.DeltaPhi(z_4vec));"
        "}";
    filling_functions["getDPhiMETLepLeading"] =
        "vector<float> getDPhiMETLepLeading(rvecf lep_pT, rvecf lep_eta, rvecf lep_phi, float MET, float MET_phi) {"
            "TLorentzVector l0_4vec;"
            "l0_4vec.SetPtEtaPhiM(lep_pT[0],lep_eta[0],lep_phi[0],0);"
            ""
            "TLorentzVector met_4vec;"
            "met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);"
            ""
            "return fabs(met_4vec.DeltaPhi(lep0_4vec));"
        "}";
    filling_functions["getDPhiMETLepSecond"] =
        "vector<float> getDPhiMETLepSecond(rvecf lep_pT, rvecf lep_eta, rvecf lep_phi, float MET, float MET_phi) {"
            "TLorentzVector l1_4vec;"
            "l1_4vec.SetPtEtaPhiM(lep_pT[1],lep_eta[1],lep_phi[1],0);"
            ""
            "TLorentzVector met_4vec;"
            "met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);"
            ""
            "return fabs(met_4vec.DeltaPhi(lep1_4vec));"
        "}";

    return filling_functions;
}

void ReductionStep(bool unit_testing) {
    Options options;

    //--- input/output
    options.sampleID = "photon";
    options.is_photon = (photonOrBackground == "photon");
    options.is_data = (sampleID == "data");
    options.period = "data15-16";
    if (!options.is_data) {
        if (options.period == "data15-16") options.period = "mc16a";
        else if (options.period == "data17") options.period = "mc16cd";
        else if (options.period == "data18") options.period = "mc16e";
    }

    options.in_file_name = "/eos/atlas/atlascerngroupdisk/phys-susy/2L2J-ANA-SUSY-2018-05/SusySkim2LJets/v1.7/JETM4/JETM4_Data/data15-16_merged_processed.root";
    options.in_tree_name = "data15-16";
    options.out_file_name = "test.root";
    options.out_tree_name = "BaselineTree";

    //--- branches to copy from old tree to new tree
    options.branches_to_copy = vector<string> {
        "lepIsoFCTight", "lepIsPR", "nLep_signal", "nLep_base",
        "lepEta", "lepPhi", "lepM", "lepFlavor", "lepCharge", "lepPt",
        "channel",
        "PhotonConversionType",
        "met_Phi",
        "nBJet20_MV2c10_FixedCutBEff_77", "nJet30", "jetM", "jetPt", "Ht30",
        "dPhiMetJet1", "minDPhi2JetsMet", // dPhiMetJet1 is not in abs value, but that should be ok
        "genWeight", "eventWeight", "leptonWeight", "jvtWeight", "bTagWeight", "pileupWeight", "globalDiLepTrigSF",
        "RunNumber", "RandomRunNumber",
    };

    //--- branches to rename and copy
    options.branches_to_rename = BranchRenameOptions {
        make_tuple("nBJet30_MV2c10_FixedCutBEff_77", "bjet_n"),
        make_tuple("jetEta", "jet_eta"),
        make_tuple("jetPhi", "jet_phi"),
        make_tuple("met_Sign", "MET_sig"),
    };

    //--- functions used for adding new branches
    auto ff = getFillingFunctions(options.is_data);

    //--- new branches to add
    options.branches_to_add = BranchAddOptions {
        make_tuple("dPhiMetJet", ff["getDPhiMetJet"], "getDPhiMetJet(jetPt, jetEta, jetPhi, nJet30, met_Et, met_Phi)"),
        make_tuple("dPhiMetJet2", "", "dPhiMetJet[1]"),
        make_tuple("dPhiMetJet12Min", "", "std::min(dPhiMetJet[0], dPhiMetJet[1])"),
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
            make_tuple("METt_unsmeared", "", "met_Et*sin(met_Phi-PhotonPt)"),
            make_tuple("METl_unsmeared", "", "met_Et*cos(met_Phi-PhotonPhi)"),
            make_tuple("trigMatch_2LTrig", "", "1"),
            make_tuple("trigMatch_2LTrigOR", "", "1"),
        };
        if (options.is_data)
            additional_add.push_back(make_tuple("totalWeight", ff["getPhotonWeight"], "getPhotonWeight(trigMatch_HLT_g15_loose_L1EM7,"
                "trigPrescale_HLT_g15_loose_L1EM7, trigMatch_HLT_g25_loose_L1EM15, trigPrescale_HLT_g25_loose_L1EM15, trigMatch_HLT_g35_loose_L1EM15,"
                "trigPrescale_HLT_g35_loose_L1EM15, trigMatch_HLT_g40_loose_L1EM15, trigPrescale_HLT_g40_loose_L1EM15, trigMatch_HLT_g45_loose_L1EM15,"
                "trigPrescale_HLT_g45_loose_L1EM15, trigMatch_HLT_g50_loose_L1EM15, trigPrescale_HLT_g50_loose_L1EM15, trigMatch_HLT_g60_loose,"
                "trigPrescale_HLT_g60_loose, trigMatch_HLT_g70_loose, trigPrescale_HLT_g70_loose, trigMatch_HLT_g80_loose, trigPrescale_HLT_g80_loose,"
                "trigMatch_HLT_g100_loose, trigPrescale_HLT_g100_loose, trigMatch_HLT_g140_loose, trigPrescale_HLT_g140_loose, PhotonPt)"));
        else if (!options.is_data) {
            float lumi = GetLumi(options.period);
            if (TString(options.sampleID).Contains("Vg"))
                lumi *= -1;
            additional_add.push_back(make_tuple("totalWeight", ff["getPhotonWeight"], "getPhotonWeight(lumi, genWeight, eventWeight, jvtWeight, bTagWeight, pileupWeight"));
        }
    }
    if (!options.is_photon) {
        additional_copy = vector<string> {
            "trigMatch_2LTrig", "trigMatch_2LTrigOR",
            "met_Et",
            "mll", "Ptll",
        };
        additional_add = BranchAddOptions {
            make_tuple("is_OS", "", "lepCharge[0]!=lepCharge[1]"),
            make_tuple("Z_eta", ff["getZEta"], "getZEta(lepPt, lepEta, lepPhi)"),
            make_tuple("Z_phi", ff["getZPhi"], "getZPhi(lepPt, lepEta, lepPhi)"),
            make_tuple("METt", "", "met_Et*sin(met_Phi-Z_phi)"),
            make_tuple("METl", "", "met_Et*cos(met_Phi-Z_phi)"),
            make_tuple("Z_cm_lep_theta", ff["getZCMLepTheta"], "getZCMLepTheta(lepPt, lepEta, lepPhi, Ptll, Z_eta, Z_phi)"),
            make_tuple("DR_2Lep", ff["getDR2Lep"], "getZCMLepTheta(lepPt, lepEta, lepPhi)"),
            make_tuple("DPhi_2Lep", ff["getDPhi2Lep"], "getZCMLepTheta(lepPt, lepEta, lepPhi)"),
            make_tuple("DPhi_METZPhoton", ff["getDPhiMETZPhoton"], "getZCMLepTheta(Ptll, Z_eta, Z_phi, met_Et, met_Phi)"),
            make_tuple("DPhi_METLepLeading", ff["getDPhiMETLepLeading"], "getZCMLepTheta(lepPt, lepEta, lepPhi, met_Et, met_Phi)"),
            make_tuple("DPhi_METLepSecond", ff["getDPhiMETLepSecond"], "getZCMLepTheta(lepPt, lepEta, lepPhi, met_Et, met_Phi)"),
            make_tuple("DPhi_METLepMin", "", "std::min(DPhi_METLepLeading, DPhi_METLepSecond)"),
        };
        if (options.is_data)
            additional_add.push_back(make_tuple("totalWeight", "", "1"));
        else if (!options.is_data)
            additional_add.push_back(make_tuple("totalWeight", "", "lumi*genWeight*eventWeight*leptonWeight*jvtWeight*bTagWeight*pileupWeight*FFWeight"));
    }

    options.branches_to_copy.insert(options.branches_to_copy.end(), additional_copy.begin(), additional_copy.end());
    options.branches_to_rename.insert(options.branches_to_rename.end(), additional_rename.begin(), additional_rename.end());
    options.branches_to_add.insert(options.branches_to_add.end(), additional_add.begin(), additional_add.end());

    //--- set selection cut
    if (!options.is_photon)
        options.cut = cuts::bkg_baseline;
    if (options.is_photon)
        options.cut = cuts::photon_baseline_ntuples;

    //--- make reduced ntuples
    options.unit_testing = unit_testing;
    ReduceNtuples(options);
}

//---------------
// MAIN FUNCTION
//---------------

void Main() {
    ROOT::EnableImplicitMT(); // enable parallelization to speed up RDataFrame

    bool unit_testing = false;
    ReductionStep(unit_testing);
    //TH1::SetDefaultSumw2();

    ////--- open input and output files and make TTrees
    //auto [inTree, outTree, inFile, outFile] = openTTrees(inFolder, outFolder, period, sampleID, isData, is_photon);
//}

//tuple<TTree*, TTree*, TFile*, TFile*> openTTrees(string inFolder, string outFolder, string period, string sampleID, bool isData, bool is_photon) {
    ///// open input and output files, get TTrees
    //string infilename = Form("%s%s/%s_merged_processed.root", inFolder.c_str(), period.c_str(), sampleID.c_str()); 
    //if (isData) infilename = Form("%s/%s_merged_processed.root", inFolder.c_str(), period.c_str()); 

    //string outfilename = ntuple_path + "/" + outFolder + "/" + period.c_str() + "_" + sampleID.c_str() + ".root";
    //if (isData) {
       //if (is_photon) outfilename = ntuple_path + "/" + outFolder + "/" + period.c_str() + "_photon.root";
       //else outfilename = ntuple_path + "/" + outFolder + "/" + period.c_str() + "_bkg.root";
    //}

    //string treeName = sampleID + "_NoSys";
    //if (isData) {
       //if (is_photon) treeName = period;
       //else treeName = "data";
    //}

    //cout << "Opening file           : " << infilename << endl;
    //cout << "Tree name              : " << treeName << endl;

    //TFile* inFile = TFile::Open(infilename.c_str());
    //TTree* inTree = (TTree*)inFile->Get(treeName.c_str());

    //cout << "Events in tree         : " << inTree->GetEntries() << endl;
    //cout << "Writing to             : " << outfilename << endl;
    //cout << endl;

    //TFile* outFile = TFile::Open(outfilename.c_str(), "recreate");
    //TTree* outTree = new TTree("BaselineTree", "baseline tree");

    //return make_tuple(inTree, outTree, inFile, outFile);
//}
}

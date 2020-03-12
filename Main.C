#include "Common/Settings.C"
#include "ReduceNtuples.C"

using namespace std;

void ReductionStep(bool unit_testing) {
    Options options;

    //--- input/output
    options.in_file_name = "/public/data/Photon/Ntuples/bkg_data/data15-16_bkg.root";
    options.in_tree_name = "BaselineTree";
    options.out_file_name = "/public/data/Photon/SkimmedSamples/data15-16_bkg.root";
    options.out_tree_name = "BaselineTree";
    options.is_photon = false;

    //--- branches to copy from old tree to new tree
    options.branches_to_copy = vector<string> {
        "lepIsoFCTight", "lepIsPR", "nLep_signal", "nLep_base",
        "lepEta", "lepPhi", "lepM", "lepFlavor", "lepCharge", "lepPt",
        "channel",
        "PhotonConversionType",
        "met_Phi",
        "mll",
        "Ptll",
        "nBJet20_MV2c10_FixedCutBEff_77", "nJet30", "jetM", "jetPt", "Ht30",
        "minDPhi2JetsMet",
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

    //--- new branches to add, with functions for calculating them
    string getDPhiMetJet =
        "vector<double> myFunc(vector<double> jet_pT, vector<double> jet_eta, vector<double> jet_phi, int jet_n, double MET, double MET_phi) {"
            "TLorentzVector jet_4vec, met_4vec;"
            "met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);"
            ""
            "vector<double> dPhiMetJet;"
            "for (int i=0; i<jet_n; i++) {"
            "   jet_4vec.SetPtEtaPhiM(jet_pT->at(i),jet_eta->at(i),jet_phi->at(i),0);"
            "   dPhiMetJet.push_back(fabs(met_4vec.DeltaPhi(jet_4vec)));"
            "}"
            "return dPhiMetJet;"
        "}";

    options.branches_to_add = BranchAddOptions {
        make_tuple("dPhiMetJet", getDPhiMetJet, "getDPhiMetJet(jet_pT, jet_eta, jet_phi, jet_n, met_Et, MET_phi)"),
        make_tuple("dPhiMetJet1", "", "dPhiMetJet[0]"),
        make_tuple("dPhiMetJet2", "", "dPhiMetJet[1]"),
        make_tuple("dPhiMetJet12Min", "", "min(dPhiMetJet1, dPhiMetJet2)"),
        //make_tuple("Z_cm_lep_theta", ""),
    };

    //--- set selection cut
    if (!options.is_photon)
        options.cut = cuts::bkg_baseline;
    if (options.is_photon)
        options.cut = cuts::photon_baseline_ntuples;
    
    //--- photon/bkg specific branches
    vector<string> additional_copy;
    BranchRenameOptions additional_rename;
    BranchAddOptions additional_add;

    string getPhotonWeight =
        "double getWeight(double t) {"
            "totalWeight = 0;"
            "if (trigMatch_HLT_g15_loose_L1EM7 ==1 && gamma_pt>(15) && gamma_pt<(25+5)) totalWeight = trigPrescale_HLT_g15_loose_L1EM7;"
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
            "if (totalWeight==0) continue;"
            ""
            "if (!isData) {
                "totalWeight = lumi * genWeight * eventWeight * jvtWeight * bTagWeight * pileupWeight;"
                "//totalWeight = lumi * genWeight * eventWeight * jvtWeight * bTagWeight;"
                "if( TString(sampleID).Contains("Vg") ) totalWeight = -1.0 * totalWeight;"
            "}
            ""
            "//--- fix for large photon sample spikes
            "if (totalWeight > 100000000000) continue;"
        "}";

    string getBkgWeight =
        "totalWeight = 1;";
        //if (!isData) totalWeight = lumi * genWeight * eventWeight * leptonWeight * jvtWeight * bTagWeight * pileupWeight * FFWeight;
        ////if (!isData) totalWeight = lumi * genWeight * eventWeight * leptonWeight * jvtWeight * bTagWeight * FFWeight;

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
            //make_tuple("totalWeight", getPhotonWeight, ""),
        };
    }
    if (!options.is_photon) {
        additional_copy = vector<string> {
            "trigMatch_2LTrig", "trigMatch_2LTrigOR",
            "met_Et",
        };
        additional_add = BranchAddOptions {
            make_tuple("is_OS", "", "lepCharge[0]!=lepCharge[1]"),
            //make_tuple("totalWeight", getBkgWeight, ""),
        };
            ////--- compute 4-vectors of objects
            //TLorentzVector lep0_4vec, lep1_4vec;
            //lep0_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
            //lep1_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);

            //TLorentzVector z_4vec;
            //Z_eta = (lep0_4vec+lep1_4vec).Eta();
            //Z_phi = (lep0_4vec+lep1_4vec).Phi();
            //z_4vec.SetPtEtaPhiM(Z_pt,Z_eta,Z_phi,91.1876);

            //TLorentzVector met_4vec;
            //met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);

            ////--- compute DR and DPhi between objects
            //DR_2Lep = lep0_4vec.DeltaR(lep1_4vec);
            //DPhi_METPhoton = fabs(met_4vec.DeltaPhi(z_4vec));
            //DPhi_2Lep = fabs(lep0_4vec.DeltaPhi(lep1_4vec));
            //DPhi_METLepLeading = fabs(met_4vec.DeltaPhi(lep0_4vec));
            //DPhi_METLepSecond = fabs(met_4vec.DeltaPhi(lep1_4vec));
            //DPhi_METLepMin = min(DPhi_METLepLeading,DPhi_METLepSecond);

            ////--- compute MET parallel and perpendicular components
            //METt = MET*TMath::Sin(MET_phi-Z_phi);
            //METl = MET*TMath::Cos(MET_phi-Z_phi);

            ////--- lepton angles in Z rest frame
            //TLorentzVector l0_cm_4vec, l1_cm_4vec;
            //l0_cm_4vec = lep0_4vec;
            //l1_cm_4vec = lep1_4vec;
            //TVector3 boost_vec(0, 0, -z_4vec.BoostVector().Mag());
            //l0_cm_4vec.RotateZ(-z_4vec.Phi());
            //l0_cm_4vec.RotateY(-z_4vec.Theta());
            //l0_cm_4vec.Boost(boost_vec);
            //l1_cm_4vec.RotateZ(-z_4vec.Phi());
            //l1_cm_4vec.RotateY(-z_4vec.Theta());
            //l1_cm_4vec.Boost(boost_vec);
            ////cout << lep0_4vec.Theta() << ", " << lep1_4vec.Theta() << ", (" << l0_cm_4vec.Theta() << ", " << l0_cm_4vec.Phi() << "), (" << l1_cm_4vec.Theta() << ", " << l1_cm_4vec.Phi() << ")" << endl;
            //Z_cm_lep_theta->clear();
            //Z_cm_lep_theta->push_back(l0_cm_4vec.Theta());
            //Z_cm_lep_theta->push_back(l1_cm_4vec.Theta());
        //outTree->Branch("Z_eta",&Z_eta,"Z_eta/F");
        //outTree->Branch("Z_phi",&Z_phi,"Z_phi/F");
        //outTree->Branch("DR_2Lep",&DR_2Lep,"DR_2Lep/F");
        //outTree->Branch("DPhi_METPhoton",&DPhi_METPhoton,"DPhi_METPhoton/F");
        //outTree->Branch("DPhi_2Lep",&DPhi_2Lep,"DPhi_2Lep/F");
        //outTree->Branch("DPhi_METLepLeading",&DPhi_METLepLeading,"DPhi_METLepLeading/F");
        //outTree->Branch("DPhi_METLepSecond",&DPhi_METLepSecond,"DPhi_METLepSecond/F");
        //outTree->Branch("DPhi_METLepMin",&DPhi_METLepMin,"DPhi_METLepMin/F");
    }

    options.branches_to_copy.insert(options.branches_to_copy.end(), additional_copy.begin(), additional_copy.end());
    options.branches_to_rename.insert(options.branches_to_rename.end(), additional_rename.begin(), additional_rename.end());
    options.branches_to_add.insert(options.branches_to_add.end(), additional_add.begin(), additional_add.end());

    //--- make reduced ntuples
    options.unit_testing = unit_testing;
    ReduceNtuples(options);
}

void Main() {
    bool unit_testing = true;
    ReductionStep(unit_testing);
}

//void blah_main() {
    //TH1::SetDefaultSumw2();
    //bool isData = (sampleID == "data");
    //bool is_photon = (photonOrBackground == "photon");
    //if (!isData) {
        //if (period == "data15-16") period = "mc16a";
        //else if (period == "data17") period = "mc16cd";
        //else if (period == "data18") period = "mc16e";
    //}
    //float lumi = GetLumi(period);
	//cout << "Using luminosity       : " << lumi          << endl;
    //cout << endl;

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

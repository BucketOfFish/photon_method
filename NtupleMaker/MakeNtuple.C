#include "../Common/Settings.C"

using namespace std;

void MakeNtuple(string outputFolder, string period, string pathToNtuples, string sampleID, string photonOrBackground, int everyNEntries=10) {

    //---------------------------------------------
    // open input and output files, get TTrees
    //---------------------------------------------

    TH1::SetDefaultSumw2();

    bool isData = (sampleID == "data");
    bool isPhoton = (photonOrBackground == "photon");

    if (!isData) {
        if (period == "data15-16") period = "mc16a";
        else if (period == "data17") period = "mc16cd";
        else if (period == "data18") period = "mc16e";
    }
    float lumi = GetLumi(period);

    string filename = Form("%s%s/%s_merged_processed.root", pathToNtuples.c_str(), period.c_str(), sampleID.c_str()); 
    if (isData) filename = Form("%s/%s_merged_processed.root", pathToNtuples.c_str(), period.c_str()); 
    TFile* inputFile = TFile::Open(filename.c_str());
    string treeName = sampleID + "_NoSys";
    if (isData) {
       if (isPhoton) treeName = period;
       else treeName = "data";
    }
    cout << treeName << endl;
    TTree* inputTree = (TTree*)inputFile->Get(treeName.c_str());

    cout << endl;
    cout << "Opening file           : " << filename        << endl;
    cout << "Events in ntuple       : " << inputTree->GetEntries() << endl;
	cout << "using luminosity       : " << lumi          << endl;

    string outfilename = ntuple_path + "/" + outputFolder + "/" + period.c_str() + "_" + sampleID.c_str() + ".root";
    if (isData) {
       if (isPhoton) outfilename = ntuple_path + "/" + outputFolder + "/" + period.c_str() + "_photon.root";
       else outfilename = ntuple_path + "/" + outputFolder + "/" + period.c_str() + "_bkg.root";
    }
    cout << "Writing to : " << outfilename << endl;
    TFile outputFile( outfilename.c_str(), "recreate" );
    TTree* BaselineTree = new TTree("BaselineTree", "baseline tree");

    //-----------------------------
    // access, copy, and create branches
    //-----------------------------

    inputTree->SetBranchStatus("*", 0);

    //--- MC weights
    Double_t genWeight; SetInputBranch(inputTree, "genWeight", &genWeight);
    Double_t eventWeight; SetInputBranch(inputTree, "eventWeight", &eventWeight);
    Double_t jvtWeight; SetInputBranch(inputTree, "jvtWeight", &jvtWeight);
    Double_t bTagWeight; SetInputBranch(inputTree, "bTagWeight", &bTagWeight);
    Double_t pileupWeight; SetInputBranch(inputTree, "pileupWeight", &pileupWeight);
    double totalWeight; BaselineTree->Branch("totalWeight",&totalWeight,"totalWeight/D");
    //--- photon weights
	int trigMatch_HLT_g15_loose_L1EM7; SetInputBranch(inputTree, "trigMatch_HLT_g15_loose_L1EM7", &trigMatch_HLT_g15_loose_L1EM7);
	int trigMatch_HLT_g25_loose_L1EM15; SetInputBranch(inputTree, "trigMatch_HLT_g25_loose_L1EM15", &trigMatch_HLT_g25_loose_L1EM15);
	int trigMatch_HLT_g35_loose_L1EM15; SetInputBranch(inputTree, "trigMatch_HLT_g35_loose_L1EM15", &trigMatch_HLT_g35_loose_L1EM15);
	int trigMatch_HLT_g40_loose_L1EM15; SetInputBranch(inputTree, "trigMatch_HLT_g40_loose_L1EM15", &trigMatch_HLT_g40_loose_L1EM15);
	int trigMatch_HLT_g45_loose_L1EM15; SetInputBranch(inputTree, "trigMatch_HLT_g45_loose_L1EM15", &trigMatch_HLT_g45_loose_L1EM15);
	int trigMatch_HLT_g50_loose_L1EM15; SetInputBranch(inputTree, "trigMatch_HLT_g50_loose_L1EM15", &trigMatch_HLT_g50_loose_L1EM15);
	int trigMatch_HLT_g60_loose; SetInputBranch(inputTree, "trigMatch_HLT_g60_loose", &trigMatch_HLT_g60_loose);
	int trigMatch_HLT_g70_loose; SetInputBranch(inputTree, "trigMatch_HLT_g70_loose", &trigMatch_HLT_g70_loose);
	int trigMatch_HLT_g80_loose; SetInputBranch(inputTree, "trigMatch_HLT_g80_loose", &trigMatch_HLT_g80_loose);
	int trigMatch_HLT_g100_loose; SetInputBranch(inputTree, "trigMatch_HLT_g100_loose", &trigMatch_HLT_g100_loose);
	int trigMatch_HLT_g120_loose; SetInputBranch(inputTree, "trigMatch_HLT_g120_loose", &trigMatch_HLT_g120_loose);
	int trigMatch_HLT_g140_loose; SetInputBranch(inputTree, "trigMatch_HLT_g140_loose", &trigMatch_HLT_g140_loose);
	float trigPrescale_HLT_g15_loose_L1EM7; SetInputBranch(inputTree, "trigPrescale_HLT_g15_loose_L1EM7", &trigPrescale_HLT_g15_loose_L1EM7);
	float trigPrescale_HLT_g25_loose_L1EM15; SetInputBranch(inputTree, "trigPrescale_HLT_g25_loose_L1EM15", &trigPrescale_HLT_g25_loose_L1EM15);
	float trigPrescale_HLT_g35_loose_L1EM15; SetInputBranch(inputTree, "trigPrescale_HLT_g35_loose_L1EM15", &trigPrescale_HLT_g35_loose_L1EM15);
	float trigPrescale_HLT_g40_loose_L1EM15; SetInputBranch(inputTree, "trigPrescale_HLT_g40_loose_L1EM15", &trigPrescale_HLT_g40_loose_L1EM15);
	float trigPrescale_HLT_g45_loose_L1EM15; SetInputBranch(inputTree, "trigPrescale_HLT_g45_loose_L1EM15", &trigPrescale_HLT_g45_loose_L1EM15);
	float trigPrescale_HLT_g50_loose_L1EM15; SetInputBranch(inputTree, "trigPrescale_HLT_g50_loose_L1EM15", &trigPrescale_HLT_g50_loose_L1EM15);
	float trigPrescale_HLT_g60_loose; SetInputBranch(inputTree, "trigPrescale_HLT_g60_loose", &trigPrescale_HLT_g60_loose);
	float trigPrescale_HLT_g70_loose; SetInputBranch(inputTree, "trigPrescale_HLT_g70_loose", &trigPrescale_HLT_g70_loose);
	float trigPrescale_HLT_g80_loose; SetInputBranch(inputTree, "trigPrescale_HLT_g80_loose", &trigPrescale_HLT_g80_loose);
	float trigPrescale_HLT_g100_loose; SetInputBranch(inputTree, "trigPrescale_HLT_g100_loose", &trigPrescale_HLT_g100_loose);
	float trigPrescale_HLT_g120_loose; SetInputBranch(inputTree, "trigPrescale_HLT_g120_loose", &trigPrescale_HLT_g120_loose);
	float trigPrescale_HLT_g140_loose; SetInputBranch(inputTree, "trigPrescale_HLT_g140_loose", &trigPrescale_HLT_g140_loose);
    //--- non-photon weights
    Double_t leptonWeight; SetInputBranch(inputTree, "leptonWeight", &leptonWeight);
    Double_t FFWeight; SetInputBranch(inputTree, "FFWeight", &FFWeight);
    bool trigMatch_1L2LTrigOR; SetInputBranch(inputTree, "trigMatch_1L2LTrigOR", &trigMatch_1L2LTrigOR);

    //--- event selection
    vector<int>* lepIsoFCTight = new vector<int>; CopyBranch(inputTree, BaselineTree, "lepIsoFCTight", "lepIsoFCTight", &lepIsoFCTight, "vector<int>");
    vector<float>* lep_eta = new vector<float>; CopyBranch(inputTree, BaselineTree, "lepEta", "lep_eta", &lep_eta, "vector<float>");
    vector<float>* lep_phi = new vector<float>; CopyBranch(inputTree, BaselineTree, "lepPhi", "lep_phi", &lep_phi, "vector<float>");
    int channel, is_OS;
    if (!isPhoton) {
        BaselineTree->Branch("channel", &channel, "channel/I");
        BaselineTree->Branch("is_OS", &is_OS, "is_OS/I");
    }

    //--- photon conversion types
    // 0 = unconverted;
    // 1 = single track with silicon hit; 2 = single track, no silicon
    // 3 = double track with silicon hits; 4 = double track, no silicon; 5 = double track, one silicon hit
    int photon_conversion_type; CopyBranch(inputTree, BaselineTree, "PhotonConversionType", "PhotonConversionType", &photon_conversion_type, "I");

    //--- MET components, and DR and DPhi between objects
    float MET_phi; CopyBranch(inputTree, BaselineTree, "met_Phi", "met_Phi", &MET_phi, "F");
    float MET, METl, METt;
    if (!isPhoton) {
        CopyBranch(inputTree, BaselineTree, "met_Et", "met_Et", &MET, "F");
        BaselineTree->Branch("METl",&METl,"METl/F");
        BaselineTree->Branch("METt",&METt,"METt/F");
    }
    else {
        CopyBranch(inputTree, BaselineTree, "met_Et", "met_Et_raw", &MET, "F");
        BaselineTree->Branch("METl_raw",&METl,"METl/F");
        BaselineTree->Branch("METt_raw",&METt,"METt/F");
    }
    //--- non-photon
    float Z_eta, Z_phi, DR_2Lep, DPhi_METPhoton, DPhi_2Lep, DPhi_METLepLeading, DPhi_METLepSecond, DPhi_METLepMin;
    if (!isPhoton) {
        BaselineTree->Branch("Z_eta",&Z_eta,"Z_eta/F");
        BaselineTree->Branch("Z_phi",&Z_phi,"Z_phi/F");
        BaselineTree->Branch("DR_2Lep",&DR_2Lep,"DR_2Lep/F");
        BaselineTree->Branch("DPhi_METPhoton",&DPhi_METPhoton,"DPhi_METPhoton/F");
        BaselineTree->Branch("DPhi_2Lep",&DPhi_2Lep,"DPhi_2Lep/F");
        BaselineTree->Branch("DPhi_METLepLeading",&DPhi_METLepLeading,"DPhi_METLepLeading/F");
        BaselineTree->Branch("DPhi_METLepSecond",&DPhi_METLepSecond,"DPhi_METLepSecond/F");
        BaselineTree->Branch("DPhi_METLepMin",&DPhi_METLepMin,"DPhi_METLepMin/F");
    }
    //--- photon-only
	float gamma_pt; CopyBranch(inputTree, BaselineTree, "PhotonPt", "gamma_pt", &gamma_pt, "F");
	float gamma_eta; CopyBranch(inputTree, BaselineTree, "PhotonEta", "gamma_eta", &gamma_eta, "F");
	float gamma_phi; CopyBranch(inputTree, BaselineTree, "PhotonPhi", "gamma_phi", &gamma_phi, "F");

    //--- not used here?
    Int_t RunNumber; CopyBranch(inputTree, BaselineTree, "RunNumber", "RunNumber", &RunNumber, "I");
    Int_t bjet_n; CopyBranch(inputTree, BaselineTree, "nBJet30_MV2c10_FixedCutBEff_77", "bjet_n", &bjet_n, "I");
    float HT; CopyBranch(inputTree, BaselineTree, "Ht30", "HT", &HT, "F");
    vector<float>* jet_pT = new vector<float>; CopyBranch(inputTree, BaselineTree, "jetPt", "jet_pT", &jet_pT, "vector<float>");
    vector<float>* jet_eta = new vector<float>; CopyBranch(inputTree, BaselineTree, "jetEta", "jet_eta", &jet_eta, "vector<float>");
    vector<float>* jet_phi = new vector<float>; CopyBranch(inputTree, BaselineTree, "jetPhi", "jet_phi", &jet_phi, "vector<float>");

    //--- HistFitter branches
    CopyAllBranches(inputTree, BaselineTree, histFitterBranches);

    float mll; CopyBranch(inputTree, BaselineTree, "mll", "mll", &mll, "F");
    float Z_pt; CopyBranch(inputTree, BaselineTree, "Ptll", "Ptll", &Z_pt, "F");
    float dPhiMetJet1; BaselineTree->Branch("dPhiMetJet1",&dPhiMetJet1,"dPhiMetJet1/F");
    float dPhiMetJet2; BaselineTree->Branch("dPhiMetJet2",&dPhiMetJet2,"dPhiMetJet2/F");
    float dPhiMetJet12Min; BaselineTree->Branch("dPhiMetJet12Min",&dPhiMetJet12Min,"dPhiMetJet12Min/F");

    vector<float>* jet_m = new vector<float>; CopyBranch(inputTree, BaselineTree, "jetM", "jetM", &jet_m, "vector<float>");
    vector<int>* lepFlavor = new vector<int>; CopyBranch(inputTree, BaselineTree, "lepFlavor", "lepFlavor", &lepFlavor, "vector<int>");
    vector<int>* lepCharge = new vector<int>; CopyBranch(inputTree, BaselineTree, "lepCharge", "lepCharge", &lepCharge, "vector<int>");
    vector<float>* lep_pT = new vector<float>; CopyBranch(inputTree, BaselineTree, "lepPt", "lepPt", &lep_pT, "vector<float>");
    int nBJet20_MV2c10_FixedCutBEff_77; CopyBranch(inputTree, BaselineTree, "nBJet20_MV2c10_FixedCutBEff_77", "nBJet20_MV2c10_FixedCutBEff_77", &nBJet20_MV2c10_FixedCutBEff_77, "I");
    Int_t jet_n; CopyBranch(inputTree, BaselineTree, "nJet30", "nJet30", &jet_n, "I");
    Int_t nLep_signal; CopyBranch(inputTree, BaselineTree, "nLep_signal", "nLep_signal", &nLep_signal, "I");
    Int_t nLep_base; CopyBranch(inputTree, BaselineTree, "nLep_base", "nLep_base", &nLep_base, "I");
    bool trigMatch_2LTrigOR; CopyBranch(inputTree, BaselineTree, "trigMatch_2LTrigOR", "trigMatch_2LTrigOR", &trigMatch_2LTrigOR, "O");

    // lepton angular distribution
    vector<float> *Z_cm_lep_theta = new vector<float>; BaselineTree->Branch("Z_cm_lep_theta","vector<float>",&Z_cm_lep_theta);

    //-----------------------------
    // loop over events
    //-----------------------------

    Long64_t nentries = inputTree->GetEntries();

    for (Long64_t i=0; i<nentries; i+=everyNEntries) {

        if (fmod(i,1e5)==0) cout << i << " events processed." << endl;
        inputTree->GetEntry(i);

        //--- event selection
        if (isPhoton) {
            if (nLep_base > 0) continue;
            if (jet_n < 1) continue;
            if (gamma_pt<15.) continue;
        }
        else {
            if ( nLep_signal  != 2                 ) continue; // exactly 2 signal leptons
            if ( nLep_base    != 2                 ) continue; // exactly 2 baseline leptons
            if ( lep_pT->at(0) < cuts::leading_lep_pt_cut ) continue; // 1st lep pT > 25 GeV
            if ( lep_pT->at(1) < cuts::second_lep_pt_cut  ) continue; // 2nd lep pT > 25 GeV
            if ( jet_n < 1   ) continue; // require at least 1 pT > 30 GeV jets
            if ( !trigMatch_1L2LTrigOR ) continue; // need 2 lepton trigger

            //--- determine channel
            channel = -1;
            if ( lepFlavor->at(0) == 2 && lepFlavor->at(1) == 2 ) channel = 0; // mumu
            if ( lepFlavor->at(0) == 1 && lepFlavor->at(1) == 1 ) channel = 1; // ee
            if ( lepFlavor->at(0) == 1 && lepFlavor->at(1) == 2 ) channel = 2; // em
            if ( lepFlavor->at(0) == 2 && lepFlavor->at(1) == 1 ) channel = 3; // me
            if ( channel < 0 ) continue; // require exactly 2 signal leptons

            //--- determine OS / SS
            is_OS = -1;
            if ( lepCharge->at(0) != lepCharge->at(1) ) is_OS = 1;
            if ( lepCharge->at(0) == lepCharge->at(1) ) is_OS = 0;
            if ( is_OS != 1 ) continue; // require opposite-sign
        }

        //--- evaluate weight
        if (isPhoton) {
            totalWeight = 0;
            if (trigMatch_HLT_g15_loose_L1EM7 ==1 && gamma_pt>(15) && gamma_pt<(25+5)) totalWeight = trigPrescale_HLT_g15_loose_L1EM7;
            if (trigMatch_HLT_g25_loose_L1EM15==1 && gamma_pt>(25+5) && gamma_pt<(35+5)) totalWeight = trigPrescale_HLT_g25_loose_L1EM15;
            if (trigMatch_HLT_g35_loose_L1EM15==1 && gamma_pt>(35+5) && gamma_pt<(40+5)) totalWeight = trigPrescale_HLT_g35_loose_L1EM15;
            if (trigMatch_HLT_g40_loose_L1EM15==1 && gamma_pt>(40+5) && gamma_pt<(45+5)) totalWeight = trigPrescale_HLT_g40_loose_L1EM15;
            if (trigMatch_HLT_g45_loose_L1EM15==1 && gamma_pt>(45+5) && gamma_pt<(50+5)) totalWeight = trigPrescale_HLT_g45_loose_L1EM15;
            if (trigMatch_HLT_g50_loose_L1EM15==1 && gamma_pt>(50+5) && gamma_pt<(60+5)) totalWeight = trigPrescale_HLT_g50_loose_L1EM15;
            if (trigMatch_HLT_g60_loose==1 && gamma_pt>(60+5) && gamma_pt<(70+5)) totalWeight = trigPrescale_HLT_g60_loose;
            if (trigMatch_HLT_g70_loose==1 && gamma_pt>(70+5) && gamma_pt<(80+5)) totalWeight = trigPrescale_HLT_g70_loose;
            if (trigMatch_HLT_g80_loose==1 && gamma_pt>(80+5) && gamma_pt<(100+5)) totalWeight = trigPrescale_HLT_g80_loose;
            if (trigMatch_HLT_g100_loose==1 && gamma_pt>(100+5) && gamma_pt<(140+5)) totalWeight = trigPrescale_HLT_g100_loose;
            if (trigMatch_HLT_g140_loose==1 && gamma_pt>(140+5)) totalWeight = trigPrescale_HLT_g140_loose;
            if (totalWeight==0) continue;

            if (!isData) {
                totalWeight = lumi * genWeight * eventWeight * jvtWeight * bTagWeight * pileupWeight;
                //totalWeight = lumi * genWeight * eventWeight * jvtWeight * bTagWeight;
                if( TString(sampleID).Contains("Vg") ) totalWeight = -1.0 * totalWeight;
            }

            //--- fix for large photon sample spikes
            if (totalWeight > 100000000000) continue;
        }
        else {
            totalWeight = 1;
            if (!isData) totalWeight = lumi * genWeight * eventWeight * leptonWeight * jvtWeight * bTagWeight * pileupWeight * FFWeight;
            //if (!isData) totalWeight = lumi * genWeight * eventWeight * leptonWeight * jvtWeight * bTagWeight * FFWeight;
        }

        //--- compute additional features
        if (isPhoton) {
            //--- compute MET parallel and perpendicular components
            METt = MET*TMath::Sin(MET_phi-gamma_phi);
            METl = MET*TMath::Cos(MET_phi-gamma_phi);
        }
        else {
            //--- compute 4-vectors of objects
            TLorentzVector lep0_4vec, lep1_4vec;
            lep0_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
            lep1_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);

            TLorentzVector z_4vec;
            Z_eta = (lep0_4vec+lep1_4vec).Eta();
            Z_phi = (lep0_4vec+lep1_4vec).Phi();
            z_4vec.SetPtEtaPhiM(Z_pt,Z_eta,Z_phi,91.1876);

            TLorentzVector met_4vec;
            met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);

            //--- compute DR and DPhi between objects
            DR_2Lep = lep0_4vec.DeltaR(lep1_4vec);
            DPhi_METPhoton = fabs(met_4vec.DeltaPhi(z_4vec));
            DPhi_2Lep = fabs(lep0_4vec.DeltaPhi(lep1_4vec));
            DPhi_METLepLeading = fabs(met_4vec.DeltaPhi(lep0_4vec));
            DPhi_METLepSecond = fabs(met_4vec.DeltaPhi(lep1_4vec));
            DPhi_METLepMin = min(DPhi_METLepLeading,DPhi_METLepSecond);

            //--- compute MET parallel and perpendicular components
            METt = MET*TMath::Sin(MET_phi-Z_phi);
            METl = MET*TMath::Cos(MET_phi-Z_phi);

            //--- lepton angles in Z rest frame
            TLorentzVector l0_cm_4vec, l1_cm_4vec;
            l0_cm_4vec = lep0_4vec;
            l1_cm_4vec = lep1_4vec;
            TVector3 boost_vec(0, 0, -z_4vec.BoostVector().Mag());
            l0_cm_4vec.RotateZ(-z_4vec.Phi());
            l0_cm_4vec.RotateY(-z_4vec.Theta());
            l0_cm_4vec.Boost(boost_vec);
            l1_cm_4vec.RotateZ(-z_4vec.Phi());
            l1_cm_4vec.RotateY(-z_4vec.Theta());
            l1_cm_4vec.Boost(boost_vec);
            //cout << lep0_4vec.Theta() << ", " << lep1_4vec.Theta() << ", (" << l0_cm_4vec.Theta() << ", " << l0_cm_4vec.Phi() << "), (" << l1_cm_4vec.Theta() << ", " << l1_cm_4vec.Phi() << ")" << endl;
            Z_cm_lep_theta->clear();
            Z_cm_lep_theta->push_back(l0_cm_4vec.Theta());
            Z_cm_lep_theta->push_back(l1_cm_4vec.Theta());
        }

        TLorentzVector jet0_4vec, jet1_4vec, met_4vec;
        jet0_4vec.SetPtEtaPhiM(jet_pT->at(0),jet_eta->at(0),jet_phi->at(0),0);
        if (jet_n > 1) jet1_4vec.SetPtEtaPhiM(jet_pT->at(1),jet_eta->at(1),jet_phi->at(1),0);
        met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);

        dPhiMetJet1 = fabs(met_4vec.DeltaPhi(jet0_4vec));
        if (jet_n > 1) dPhiMetJet2 = fabs(met_4vec.DeltaPhi(jet1_4vec));
        else dPhiMetJet2 = 999;
        dPhiMetJet12Min = min(dPhiMetJet1, dPhiMetJet2);

        BaselineTree->Fill();     
    }

    cout << "write output..." << endl;
    BaselineTree->Write();

    cout << "done." << endl;
    outputFile.Close();
    delete inputFile;
}

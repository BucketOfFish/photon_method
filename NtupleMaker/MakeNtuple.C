#include "../Common/Settings.C"

using namespace std;

tuple<TTree*, TTree*, TFile*, TFile*> openTTrees(string inFolder, string outFolder, string period, string sampleID, bool isData, bool isPhoton) {
    /// open input and output files, get TTrees
    string infilename = Form("%s%s/%s_merged_processed.root", inFolder.c_str(), period.c_str(), sampleID.c_str()); 
    if (isData) infilename = Form("%s/%s_merged_processed.root", inFolder.c_str(), period.c_str()); 

    string outfilename = ntuple_path + "/" + outFolder + "/" + period.c_str() + "_" + sampleID.c_str() + ".root";
    if (isData) {
       if (isPhoton) outfilename = ntuple_path + "/" + outFolder + "/" + period.c_str() + "_photon.root";
       else outfilename = ntuple_path + "/" + outFolder + "/" + period.c_str() + "_bkg.root";
    }

    string treeName = sampleID + "_NoSys";
    if (isData) {
       if (isPhoton) treeName = period;
       else treeName = "data";
    }

    cout << "Opening file           : " << infilename << endl;
    cout << "Tree name              : " << treeName << endl;

    TFile* inFile = TFile::Open(infilename.c_str());
    TTree* inTree = (TTree*)inFile->Get(treeName.c_str());

    cout << "Events in tree         : " << inTree->GetEntries() << endl;
    cout << "Writing to             : " << outfilename << endl;
    cout << endl;

    TFile* outFile = TFile::Open(outfilename.c_str(), "recreate");
    TTree* outTree = new TTree("BaselineTree", "baseline tree");

    return make_tuple(inTree, outTree, inFile, outFile);
}

void MakeNtuple(string inFolder, string outFolder, string period, string sampleID, string photonOrBackground, int everyNEntries=10) {
    //--- global settings
    TH1::SetDefaultSumw2();
    bool isData = (sampleID == "data");
    bool isPhoton = (photonOrBackground == "photon");
    if (!isData) {
        if (period == "data15-16") period = "mc16a";
        else if (period == "data17") period = "mc16cd";
        else if (period == "data18") period = "mc16e";
    }
    float lumi = GetLumi(period);
	cout << "Using luminosity       : " << lumi          << endl;
    cout << endl;

    //--- open input and output files and make TTrees
    auto [inTree, outTree, inFile, outFile] = openTTrees(inFolder, outFolder, period, sampleID, isData, isPhoton);

    //--- access, copy, and create branches
    inTree->SetBranchStatus("*", 0);
    //copyBranches(inTree, outTree, copyBranchNames, renameBranchNames);

    //--- MC weights
    Double_t genWeight; SetInputBranch(inTree, "genWeight", &genWeight);
    Double_t eventWeight; SetInputBranch(inTree, "eventWeight", &eventWeight);
    Double_t jvtWeight; SetInputBranch(inTree, "jvtWeight", &jvtWeight);
    Double_t bTagWeight; SetInputBranch(inTree, "bTagWeight", &bTagWeight);
    Double_t pileupWeight; SetInputBranch(inTree, "pileupWeight", &pileupWeight);
    double totalWeight; outTree->Branch("totalWeight",&totalWeight,"totalWeight/D");
    //--- photon weights
	int trigMatch_HLT_g15_loose_L1EM7; SetInputBranch(inTree, "trigMatch_HLT_g15_loose_L1EM7", &trigMatch_HLT_g15_loose_L1EM7);
	int trigMatch_HLT_g25_loose_L1EM15; SetInputBranch(inTree, "trigMatch_HLT_g25_loose_L1EM15", &trigMatch_HLT_g25_loose_L1EM15);
	int trigMatch_HLT_g35_loose_L1EM15; SetInputBranch(inTree, "trigMatch_HLT_g35_loose_L1EM15", &trigMatch_HLT_g35_loose_L1EM15);
	int trigMatch_HLT_g40_loose_L1EM15; SetInputBranch(inTree, "trigMatch_HLT_g40_loose_L1EM15", &trigMatch_HLT_g40_loose_L1EM15);
	int trigMatch_HLT_g45_loose_L1EM15; SetInputBranch(inTree, "trigMatch_HLT_g45_loose_L1EM15", &trigMatch_HLT_g45_loose_L1EM15);
	int trigMatch_HLT_g50_loose_L1EM15; SetInputBranch(inTree, "trigMatch_HLT_g50_loose_L1EM15", &trigMatch_HLT_g50_loose_L1EM15);
	int trigMatch_HLT_g60_loose; SetInputBranch(inTree, "trigMatch_HLT_g60_loose", &trigMatch_HLT_g60_loose);
	int trigMatch_HLT_g70_loose; SetInputBranch(inTree, "trigMatch_HLT_g70_loose", &trigMatch_HLT_g70_loose);
	int trigMatch_HLT_g80_loose; SetInputBranch(inTree, "trigMatch_HLT_g80_loose", &trigMatch_HLT_g80_loose);
	int trigMatch_HLT_g100_loose; SetInputBranch(inTree, "trigMatch_HLT_g100_loose", &trigMatch_HLT_g100_loose);
	int trigMatch_HLT_g120_loose; SetInputBranch(inTree, "trigMatch_HLT_g120_loose", &trigMatch_HLT_g120_loose);
	int trigMatch_HLT_g140_loose; SetInputBranch(inTree, "trigMatch_HLT_g140_loose", &trigMatch_HLT_g140_loose);
	float trigPrescale_HLT_g15_loose_L1EM7; SetInputBranch(inTree, "trigPrescale_HLT_g15_loose_L1EM7", &trigPrescale_HLT_g15_loose_L1EM7);
	float trigPrescale_HLT_g25_loose_L1EM15; SetInputBranch(inTree, "trigPrescale_HLT_g25_loose_L1EM15", &trigPrescale_HLT_g25_loose_L1EM15);
	float trigPrescale_HLT_g35_loose_L1EM15; SetInputBranch(inTree, "trigPrescale_HLT_g35_loose_L1EM15", &trigPrescale_HLT_g35_loose_L1EM15);
	float trigPrescale_HLT_g40_loose_L1EM15; SetInputBranch(inTree, "trigPrescale_HLT_g40_loose_L1EM15", &trigPrescale_HLT_g40_loose_L1EM15);
	float trigPrescale_HLT_g45_loose_L1EM15; SetInputBranch(inTree, "trigPrescale_HLT_g45_loose_L1EM15", &trigPrescale_HLT_g45_loose_L1EM15);
	float trigPrescale_HLT_g50_loose_L1EM15; SetInputBranch(inTree, "trigPrescale_HLT_g50_loose_L1EM15", &trigPrescale_HLT_g50_loose_L1EM15);
	float trigPrescale_HLT_g60_loose; SetInputBranch(inTree, "trigPrescale_HLT_g60_loose", &trigPrescale_HLT_g60_loose);
	float trigPrescale_HLT_g70_loose; SetInputBranch(inTree, "trigPrescale_HLT_g70_loose", &trigPrescale_HLT_g70_loose);
	float trigPrescale_HLT_g80_loose; SetInputBranch(inTree, "trigPrescale_HLT_g80_loose", &trigPrescale_HLT_g80_loose);
	float trigPrescale_HLT_g100_loose; SetInputBranch(inTree, "trigPrescale_HLT_g100_loose", &trigPrescale_HLT_g100_loose);
	float trigPrescale_HLT_g120_loose; SetInputBranch(inTree, "trigPrescale_HLT_g120_loose", &trigPrescale_HLT_g120_loose);
	float trigPrescale_HLT_g140_loose; SetInputBranch(inTree, "trigPrescale_HLT_g140_loose", &trigPrescale_HLT_g140_loose);
    //--- non-photon weights
    Double_t leptonWeight; SetInputBranch(inTree, "leptonWeight", &leptonWeight);
    Double_t FFWeight; SetInputBranch(inTree, "FFWeight", &FFWeight);
    bool trigMatch_1L2LTrigOR; SetInputBranch(inTree, "trigMatch_1L2LTrigOR", &trigMatch_1L2LTrigOR);

    //--- event selection
    vector<int>* lepIsoFCTight = new vector<int>; CopyBranch(inTree, outTree, "lepIsoFCTight", "lepIsoFCTight", &lepIsoFCTight, "vector<int>");
    vector<float>* lep_eta = new vector<float>; CopyBranch(inTree, outTree, "lepEta", "lep_eta", &lep_eta, "vector<float>");
    vector<float>* lep_phi = new vector<float>; CopyBranch(inTree, outTree, "lepPhi", "lep_phi", &lep_phi, "vector<float>");
    int channel, is_OS;
    if (!isPhoton) {
        outTree->Branch("channel", &channel, "channel/I");
        outTree->Branch("is_OS", &is_OS, "is_OS/I");
    }

    //--- photon conversion types
    // 0 = unconverted;
    // 1 = single track with silicon hit; 2 = single track, no silicon
    // 3 = double track with silicon hits; 4 = double track, no silicon; 5 = double track, one silicon hit
    int photon_conversion_type; CopyBranch(inTree, outTree, "PhotonConversionType", "PhotonConversionType", &photon_conversion_type, "I");

    //--- MET components, and DR and DPhi between objects
    float MET_phi; CopyBranch(inTree, outTree, "met_Phi", "met_Phi", &MET_phi, "F");
    float MET, METl, METt;
    if (!isPhoton) {
        CopyBranch(inTree, outTree, "met_Et", "met_Et", &MET, "F");
        outTree->Branch("METl",&METl,"METl/F");
        outTree->Branch("METt",&METt,"METt/F");
    }
    else {
        CopyBranch(inTree, outTree, "met_Et", "met_Et_raw", &MET, "F");
        outTree->Branch("METl_raw",&METl,"METl/F");
        outTree->Branch("METt_raw",&METt,"METt/F");
    }
    //--- non-photon
    float Z_eta, Z_phi, DR_2Lep, DPhi_METPhoton, DPhi_2Lep, DPhi_METLepLeading, DPhi_METLepSecond, DPhi_METLepMin;
    if (!isPhoton) {
        outTree->Branch("Z_eta",&Z_eta,"Z_eta/F");
        outTree->Branch("Z_phi",&Z_phi,"Z_phi/F");
        outTree->Branch("DR_2Lep",&DR_2Lep,"DR_2Lep/F");
        outTree->Branch("DPhi_METPhoton",&DPhi_METPhoton,"DPhi_METPhoton/F");
        outTree->Branch("DPhi_2Lep",&DPhi_2Lep,"DPhi_2Lep/F");
        outTree->Branch("DPhi_METLepLeading",&DPhi_METLepLeading,"DPhi_METLepLeading/F");
        outTree->Branch("DPhi_METLepSecond",&DPhi_METLepSecond,"DPhi_METLepSecond/F");
        outTree->Branch("DPhi_METLepMin",&DPhi_METLepMin,"DPhi_METLepMin/F");
    }
    //--- photon-only
	float gamma_pt; CopyBranch(inTree, outTree, "PhotonPt", "gamma_pt", &gamma_pt, "F");
	float gamma_eta; CopyBranch(inTree, outTree, "PhotonEta", "gamma_eta", &gamma_eta, "F");
	float gamma_phi; CopyBranch(inTree, outTree, "PhotonPhi", "gamma_phi", &gamma_phi, "F");

    //--- not used here?
    Int_t RunNumber; CopyBranch(inTree, outTree, "RunNumber", "RunNumber", &RunNumber, "I");
    Int_t bjet_n; CopyBranch(inTree, outTree, "nBJet30_MV2c10_FixedCutBEff_77", "bjet_n", &bjet_n, "I");
    float HT; CopyBranch(inTree, outTree, "Ht30", "HT", &HT, "F");
    vector<float>* jet_pT = new vector<float>; CopyBranch(inTree, outTree, "jetPt", "jet_pT", &jet_pT, "vector<float>");
    vector<float>* jet_eta = new vector<float>; CopyBranch(inTree, outTree, "jetEta", "jet_eta", &jet_eta, "vector<float>");
    vector<float>* jet_phi = new vector<float>; CopyBranch(inTree, outTree, "jetPhi", "jet_phi", &jet_phi, "vector<float>");

    //--- HistFitter branches
    CopyAllBranches(inTree, outTree, histFitterBranches);

    float mll; CopyBranch(inTree, outTree, "mll", "mll", &mll, "F");
    float Z_pt; CopyBranch(inTree, outTree, "Ptll", "Ptll", &Z_pt, "F");
    float dPhiMetJet1; outTree->Branch("dPhiMetJet1",&dPhiMetJet1,"dPhiMetJet1/F");
    float dPhiMetJet2; outTree->Branch("dPhiMetJet2",&dPhiMetJet2,"dPhiMetJet2/F");
    float dPhiMetJet12Min; outTree->Branch("dPhiMetJet12Min",&dPhiMetJet12Min,"dPhiMetJet12Min/F");

    vector<float>* jet_m = new vector<float>; CopyBranch(inTree, outTree, "jetM", "jetM", &jet_m, "vector<float>");
    vector<int>* lepFlavor = new vector<int>; CopyBranch(inTree, outTree, "lepFlavor", "lepFlavor", &lepFlavor, "vector<int>");
    vector<int>* lepCharge = new vector<int>; CopyBranch(inTree, outTree, "lepCharge", "lepCharge", &lepCharge, "vector<int>");
    vector<float>* lep_pT = new vector<float>; CopyBranch(inTree, outTree, "lepPt", "lepPt", &lep_pT, "vector<float>");
    int nBJet20_MV2c10_FixedCutBEff_77; CopyBranch(inTree, outTree, "nBJet20_MV2c10_FixedCutBEff_77", "nBJet20_MV2c10_FixedCutBEff_77", &nBJet20_MV2c10_FixedCutBEff_77, "I");
    Int_t jet_n; CopyBranch(inTree, outTree, "nJet30", "nJet30", &jet_n, "I");
    Int_t nLep_signal; CopyBranch(inTree, outTree, "nLep_signal", "nLep_signal", &nLep_signal, "I");
    Int_t nLep_base; CopyBranch(inTree, outTree, "nLep_base", "nLep_base", &nLep_base, "I");
    bool trigMatch_2LTrigOR; CopyBranch(inTree, outTree, "trigMatch_2LTrigOR", "trigMatch_2LTrigOR", &trigMatch_2LTrigOR, "O");
    float MET_sig; CopyBranch(inTree, outTree, "met_Sign", "MET_sig", &MET_sig, "F");

    // lepton angular distribution
    vector<float> *Z_cm_lep_theta = new vector<float>; outTree->Branch("Z_cm_lep_theta","vector<float>",&Z_cm_lep_theta);

    //-----------------------------
    // loop over events
    //-----------------------------

    cout << "Applying baseline selections" << endl;
	cout << "Background baseline    : " << cuts::bkg_baseline << endl;
	cout << "Photon baseline        : " << cuts::photon_baseline_ntuples << endl;
    if (isPhoton) inTree->Draw(">>eventList", cuts::photon_baseline_ntuples, "goff");
    else inTree->Draw(">>eventList", cuts::bkg_baseline, "goff");
    TEventList *eventList = (TEventList*)gDirectory->Get("eventList");
	cout << "N events selected      : " << eventList->GetN() << endl;

    for (Long64_t i=0; i<eventList->GetN(); i+=everyNEntries) {
        if (fmod(i,1e5)==0) cout << i << " events processed.\r" << flush;
        inTree->GetEntry(eventList->GetEntry(i));

        if (!isPhoton) {
            //--- determine channel
            channel = -1;
            if ( lepFlavor->at(0) == 2 && lepFlavor->at(1) == 2 ) channel = 0; // mumu
            if ( lepFlavor->at(0) == 1 && lepFlavor->at(1) == 1 ) channel = 1; // ee
            if ( lepFlavor->at(0) == 1 && lepFlavor->at(1) == 2 ) channel = 2; // em
            if ( lepFlavor->at(0) == 2 && lepFlavor->at(1) == 1 ) channel = 3; // me

            //--- determine OS / SS
            is_OS = (lepCharge->at(0)!=lepCharge->at(1));
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

        outTree->Fill();     
    }
    cout << endl;

    cout << "Writing output..." << endl;
    outTree->Write();

    inFile->Close();
    outFile->Close();
    delete inFile;
}

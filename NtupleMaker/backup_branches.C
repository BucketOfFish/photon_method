    //-----------------------------
    // access and copy over existing branches
    //-----------------------------

    ULong64_t EventNumber; CopyBranch(inputTree, BaselineTree, "EventNumber", "EventNumber", &EventNumber, "I");
    Float_t Mu; CopyBranch(inputTree, BaselineTree, "mu", "Mu", &Mu, "F");
    Int_t nVtx; CopyBranch(inputTree, BaselineTree, "nVtx", "nVtx", &nVtx, "I");
    Bool_t trigMatch_2LTrigOR; CopyBranch(inputTree, BaselineTree, "trigMatch_2LTrigOR", "trigMatch_2LTrigOR", &trigMatch_2LTrigOR, "I");
    float MET_loose; CopyBranch(inputTree, BaselineTree, "met_Et_loose", "MET_loose", &MET_loose, "F");
    float MET_tight; CopyBranch(inputTree, BaselineTree, "met_Et", "MET_tight", &MET_tight, "F");
    float MET_tighter; CopyBranch(inputTree, BaselineTree, "met_Et_tighter", "MET_tighter", &MET_tighter, "F");
    float MET_tenacious; CopyBranch(inputTree, BaselineTree, "met_Et_tenacious", "MET_tenacious", &MET_tenacious, "F");
    Bool_t is2Lep2Jet; CopyBranch(inputTree, BaselineTree, "is2Lep2Jet", "is2Lep2Jet", &is2Lep2Jet, "I");
    Bool_t is2L2JInt; CopyBranch(inputTree, BaselineTree, "is2L2JInt", "is2L2JInt", &is2L2JInt, "I");
    Double_t mll_RJ; CopyBranch(inputTree, BaselineTree, "mll_RJ", "mll_RJ", &mll_RJ, "F");
    Double_t R_minH2P_minH3P; CopyBranch(inputTree, BaselineTree, "R_minH2P_minH3P", "R_minH2P_minH3P", &R_minH2P_minH3P, "F");
    Double_t RPT_HT5PP; CopyBranch(inputTree, BaselineTree, "RPT_HT5PP", "RPT_HT5PP", &RPT_HT5PP, "F");
    Double_t dphiVP; CopyBranch(inputTree, BaselineTree, "dphiVP", "dphiVP", &dphiVP, "F");
    Double_t H2PP; CopyBranch(inputTree, BaselineTree, "H2PP", "H2PP", &H2PP, "F");
    Double_t H5PP; CopyBranch(inputTree, BaselineTree, "H5PP", "H5PP", &H5PP, "F");
    int nJet20; CopyBranch(inputTree, BaselineTree, "nJet20", "nJet20", &nJet20, "I");
    Double_t minDphi; CopyBranch(inputTree, BaselineTree, "minDphi", "minDphi", &minDphi, "F");
    Double_t MZ; CopyBranch(inputTree, BaselineTree, "MZ", "MZ", &MZ, "F");
    Double_t NjS; CopyBranch(inputTree, BaselineTree, "NjS", "NjS", &NjS, "I");
    Double_t NjISR; CopyBranch(inputTree, BaselineTree, "NjISR", "NjISR", &NjISR, "I");
    Double_t dphiISRI; CopyBranch(inputTree, BaselineTree, "dphiISRI", "dphiISRI", &dphiISRI, "F");
    Double_t RISR; CopyBranch(inputTree, BaselineTree, "RISR", "RISR", &RISR, "F");
    Double_t PTISR; CopyBranch(inputTree, BaselineTree, "PTISR", "PTISR", &PTISR, "F");
    Double_t PTI; CopyBranch(inputTree, BaselineTree, "PTI", "PTI", &PTI, "F");
    Double_t PTCM; CopyBranch(inputTree, BaselineTree, "PTCM", "PTCM", &PTCM, "F");
    Double_t MJ; CopyBranch(inputTree, BaselineTree, "MJ", "MJ", &MJ, "F");
    Bool_t is3Lep3Jet; CopyBranch(inputTree, BaselineTree, "is3Lep3Jet", "is3Lep3Jet", &is3Lep3Jet, "I");
    Bool_t is4Lep3Jet; CopyBranch(inputTree, BaselineTree, "is4Lep3Jet", "is4Lep3Jet", &is4Lep3Jet, "I");
    Double_t lept1sign_VR; CopyBranch(inputTree, BaselineTree, "lept1sign_VR", "lept1sign_VR", &lept1sign_VR, "I");
    Double_t lept2sign_VR; CopyBranch(inputTree, BaselineTree, "lept2sign_VR", "lept2sign_VR", &lept2sign_VR, "I");
    Double_t lept1Pt_VR; CopyBranch(inputTree, BaselineTree, "lept1Pt_VR", "lept1Pt_VR", &lept1Pt_VR, "F");
    Double_t lept2Pt_VR; CopyBranch(inputTree, BaselineTree, "lept2Pt_VR", "lept2Pt_VR", &lept2Pt_VR, "F");
    Double_t MZ_VR; CopyBranch(inputTree, BaselineTree, "MZ_VR", "MZ_VR", &MZ_VR, "F");
    Double_t MJ_VR; CopyBranch(inputTree, BaselineTree, "MJ_VR", "MJ_VR", &MJ_VR, "F");
    Double_t RISR_VR; CopyBranch(inputTree, BaselineTree, "RISR_VR", "RISR_VR", &RISR_VR, "F");
    Double_t PTISR_VR; CopyBranch(inputTree, BaselineTree, "PTISR_VR", "PTISR_VR", &PTISR_VR, "F");
    Double_t PTI_VR; CopyBranch(inputTree, BaselineTree, "PTI_VR", "PTI_VR", &PTI_VR, "F");
    Double_t PTCM_VR; CopyBranch(inputTree, BaselineTree, "PTCM_VR", "PTCM_VR", &PTCM_VR, "F");
    Double_t dphiISRI_VR; CopyBranch(inputTree, BaselineTree, "dphiISRI_VR", "dphiISRI_VR", &dphiISRI_VR, "F");
    float DPhi_METJetLeading; CopyBranch(inputTree, BaselineTree, "DPhiJ1Met", "DPhiJ1Met", &DPhi_METJetLeading, "F");
    float DPhi_METJetSecond; CopyBranch(inputTree, BaselineTree, "DPhiJ2Met", "DPhiJ2Met", &DPhi_METJetSecond, "F");

    //-----------------------------
    // add new branches and histograms
    //-----------------------------

    // dijet variables
    BaselineTree->Branch("DPhi_METNonWJet",&DPhi_METNonWJet,"DPhi_METNonWJet/F");
    BaselineTree->Branch("NonWJet_pT",&NonWJet_pT,"NonWJet_pT/F");
    BaselineTree->Branch("DPhi_METNonWminJet",&DPhi_METNonWminJet,"DPhi_METNonWminJet/F");
    BaselineTree->Branch("NonWminJet_pT",&NonWminJet_pT,"NonWminJet_pT/F");
    BaselineTree->Branch("DR_Wmin2Jet",&DR_Wmin2Jet,"DR_Wmin2Jet/F");
    BaselineTree->Branch("DR_J0J1",&DR_J0J1,"DR_J0J1/F");
    BaselineTree->Branch("mj0j1",&mj0j1,"mj0j1/F");
    BaselineTree->Branch("W01_pt",&W01_pt,"W01_pt/F");
    BaselineTree->Branch("DPhi_METW01",&DPhi_METW01,"DPhi_METW01/F");
    BaselineTree->Branch("DPhi_W01Z",&DPhi_W01Z,"DPhi_W01Z/F");
    BaselineTree->Branch("mWmin",&mWmin,"mWmin/F");
    BaselineTree->Branch("Wmin_pt",&Wmin_pt,"Wmin_pt/F");
    BaselineTree->Branch("Wmin_eta",&Wmin_eta,"Wmin_eta/F");
    BaselineTree->Branch("DPhi_METWmin",&DPhi_METWmin,"DPhi_METWmin/F");
    BaselineTree->Branch("DPhi_WminZ",&DPhi_WminZ,"DPhi_WminZ/F");


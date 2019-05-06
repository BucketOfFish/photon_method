    //-----------------------------
    // access existing branches
    //-----------------------------

    float MET; inputTree->SetBranchAddress("MET_raw", &MET);
    float DPhi_METJetLeading; inputTree->SetBranchAddress("DPhi_METJetLeading_raw", &DPhi_METJetLeading);
    float DPhi_METJetSecond; inputTree->SetBranchAddress("DPhi_METJetSecond_raw", &DPhi_METJetSecond);
    float MinDPhi_PhotonJet; inputTree->SetBranchAddress("MinDPhi_PhotonJet", &MinDPhi_PhotonJet);
    float HT; CopyBranch(inputTree, BaselineTree, "HT", "HT", &HT, "F");
    int bjet_n; CopyBranch(inputTree, BaselineTree, "bjet_n", "bjet_n", &bjet_n, "I");
    float MET_loose; CopyBranch(inputTree, BaselineTree, "MET_loose", "MET_loose", &MET_loose, "F");
    float MET_tighter; CopyBranch(inputTree, BaselineTree, "MET_tighter", "MET_tighter", &MET_tighter, "F");
    float MET_tenacious; CopyBranch(inputTree, BaselineTree, "MET_tenacious", "MET_tenacious", &MET_tenacious, "F");
    int is2Lep2Jet; CopyBranch(inputTree, BaselineTree, "is2Lep2Jet", "is2Lep2Jet", &is2Lep2Jet, "I");
    int is2L2JInt; CopyBranch(inputTree, BaselineTree, "is2L2JInt", "is2L2JInt", &is2L2JInt, "I");
    int nBJet20_MV2c10_FixedCutBEff_77; CopyBranch(inputTree, BaselineTree, "nBJet20_MV2c10_FixedCutBEff_77", "nBJet20_MV2c10_FixedCutBEff_77", &nBJet20_MV2c10_FixedCutBEff_77, "I");
    float mjj; CopyBranch(inputTree, BaselineTree, "mjj", "mjj", &mjj, "F");
    float mll_RJ; CopyBranch(inputTree, BaselineTree, "mll_RJ", "mll_RJ", &mll_RJ, "F");
    float R_minH2P_minH3P; CopyBranch(inputTree, BaselineTree, "R_minH2P_minH3P", "R_minH2P_minH3P", &R_minH2P_minH3P, "F");
    float RPT_HT5PP; CopyBranch(inputTree, BaselineTree, "RPT_HT5PP", "RPT_HT5PP", &RPT_HT5PP, "F");
    float dphiVP; CopyBranch(inputTree, BaselineTree, "dphiVP", "dphiVP", &dphiVP, "F");
    float H2PP; CopyBranch(inputTree, BaselineTree, "H2PP", "H2PP", &H2PP, "F");
    float H5PP; CopyBranch(inputTree, BaselineTree, "H5PP", "H5PP", &H5PP, "F");
    int nJet20; CopyBranch(inputTree, BaselineTree, "nJet20", "nJet20", &nJet20, "I");
    float minDphi; CopyBranch(inputTree, BaselineTree, "minDphi", "minDphi", &minDphi, "F");
    float MZ; CopyBranch(inputTree, BaselineTree, "MZ", "MZ", &MZ, "F");
    int NjS; CopyBranch(inputTree, BaselineTree, "NjS", "NjS", &NjS, "I");
    int NjISR; CopyBranch(inputTree, BaselineTree, "NjISR", "NjISR", &NjISR, "I");
    float dphiISRI; CopyBranch(inputTree, BaselineTree, "dphiISRI", "dphiISRI", &dphiISRI, "F");
    float RISR; CopyBranch(inputTree, BaselineTree, "RISR", "RISR", &RISR, "F");
    float PTISR; CopyBranch(inputTree, BaselineTree, "PTISR", "PTISR", &PTISR, "F");
    float PTI; CopyBranch(inputTree, BaselineTree, "PTI", "PTI", &PTI, "F");
    float PTCM; CopyBranch(inputTree, BaselineTree, "PTCM", "PTCM", &PTCM, "F");
    float MJ; CopyBranch(inputTree, BaselineTree, "MJ", "MJ", &MJ, "F");
    int is3Lep3Jet; CopyBranch(inputTree, BaselineTree, "is3Lep3Jet", "is3Lep3Jet", &is3Lep3Jet, "I");
    int is4Lep3Jet; CopyBranch(inputTree, BaselineTree, "is4Lep3Jet", "is4Lep3Jet", &is4Lep3Jet, "I");
    int lept1sign_VR; CopyBranch(inputTree, BaselineTree, "lept1sign_VR", "lept1sign_VR", &lept1sign_VR, "I");
    int lept2sign_VR; CopyBranch(inputTree, BaselineTree, "lept2sign_VR", "lept2sign_VR", &lept2sign_VR, "I");
    float lept1Pt_VR; CopyBranch(inputTree, BaselineTree, "lept1Pt_VR", "lept1Pt_VR", &lept1Pt_VR, "F");
    float lept2Pt_VR; CopyBranch(inputTree, BaselineTree, "lept2Pt_VR", "lept2Pt_VR", &lept2Pt_VR, "F");
    float MZ_VR; CopyBranch(inputTree, BaselineTree, "MZ_VR", "MZ_VR", &MZ_VR, "F");
    float MJ_VR; CopyBranch(inputTree, BaselineTree, "MJ_VR", "MJ_VR", &MJ_VR, "F");
    float RISR_VR; CopyBranch(inputTree, BaselineTree, "RISR_VR", "RISR_VR", &RISR_VR, "F");
    float PTISR_VR; CopyBranch(inputTree, BaselineTree, "PTISR_VR", "PTISR_VR", &PTISR_VR, "F");
    float PTI_VR; CopyBranch(inputTree, BaselineTree, "PTI_VR", "PTI_VR", &PTI_VR, "F");
    float PTCM_VR; CopyBranch(inputTree, BaselineTree, "PTCM_VR", "PTCM_VR", &PTCM_VR, "F");
    float dphiISRI_VR; CopyBranch(inputTree, BaselineTree, "dphiISRI_VR", "dphiISRI_VR", &dphiISRI_VR, "F");
    float MET_phi; CopyBranch(inputTree, BaselineTree, "MET_phi_raw", "MET_phi", &MET_phi, "F");

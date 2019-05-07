    //---------------------------------------------
    // open input and output files, get TTrees
    //---------------------------------------------

    Float_t _nGenEvents = 1.;
    TH1D* hist_EventCount = new TH1D("hist_EventCount","",3,0,3);
    if (isData == "MC") {
        cout << "Setting _nGenEvents = 1 for now NEED TO FIX" << endl;
        _nGenEvents    = 1.0;
        hist_EventCount->SetBinContent(1,1.0);
    }

    if (isData == "MC") {
        cout << "Total generated events : " << _nGenEvents     << endl;
    }

    float N_passMET100 = 0.;
    //--- in loop
        if (MET>100) N_passMET100 += 1; 

    if (isData == "MC") {
        hist_EventCount->SetBinContent(2,nentries);
        hist_EventCount->SetBinContent(3,N_passMET100);
        hist_EventCount->Write();
    }

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
    float MET_softTerm; CopyBranch(inputTree, BaselineTree, "TST_Et", "MET_softTerm", &MET_softTerm, "F");
    float MET_softPhi; CopyBranch(inputTree, BaselineTree, "TST_Phi", "MET_softPhi", &MET_softPhi, "F");
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

    TH1D* hist_cutflow_raw = new TH1D("hist_cutflow_raw","",8,0,8);
    hist_cutflow_raw->SetStats(0);
    hist_cutflow_raw->GetXaxis()->SetBinLabel(1, "2lep");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(2, "flavor");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(3, "trigger");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(4, "OS");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(5, "lep pT");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(6, "njet");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(7, "mll");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(8, "prompt");

    TH1D* hist_cutflow_weight = new TH1D("hist_cutflow_weight","",8,0,8);
    hist_cutflow_weight->SetStats(0);
    hist_cutflow_weight->GetXaxis()->SetBinLabel(1, "2lep");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(2, "flavor");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(3, "trigger");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(4, "OS");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(5, "lep pT");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(6, "njet");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(7, "mll");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(8, "prompt");

    TH1D* hist_METl_Pt[bin_size];
    TH1D* hist_METt_Pt[bin_size];
    for (int bin=0;bin<bin_size;bin++) {
        hist_METl_Pt[bin] = new TH1D(TString("hist_METl_Pt_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        hist_METl_Pt[bin]->SetStats(0);
        hist_METt_Pt[bin] = new TH1D(TString("hist_METt_Pt_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        hist_METt_Pt[bin]->SetStats(0);
    }
    TH1D* hist_Mll_dPt[dpt_bin_size];
    TH1D* hist_low_njet = new TH1D("hist_low_njet","",bin_size,njet_bin);
    TH1D* hist_low_nbjet = new TH1D("hist_low_nbjet","",bin_size,njet_bin);
    TH1D* hist_low_pt = new TH1D("hist_low_pt","",bin_size,pt_bin);
    TH1D* hist_sm_pt = new TH1D("hist_sm_pt","",bin_size,sm_pt_bin);
    TH1D* hist_low_ht = new TH1D("hist_low_ht","",bin_size,ht_bin);

    //-----------------------------
    // unused calculations
    //-----------------------------

        //--- MET soft term
        TLorentzVector tst_4vec;
        tst_4vec.SetPtEtaPhiM(MET_softTerm,0,MET_softPhi,0);
        float DPhi_TSTLepLeading = fabs(tst_4vec.DeltaPhi(lep0_4vec));
        float DPhi_TSTLepSecond = fabs(tst_4vec.DeltaPhi(lep1_4vec));
        float DPhi_TSTLepMin = min(DPhi_TSTLepLeading,DPhi_TSTLepSecond);

        //--- Oslo's MET_rel
        float MinDR_Lep0Jet = 1000.;
        float MinDR_Lep1Jet = 1000.;
        float MinDR_PhotonJet = 1000.;
        MinDPhi_PhotonJet = 1000.;
        TLorentzVector jet_4vec;
        for (unsigned int j=0;j<jet_pT->size();j++) {
            jet_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),jet_m->at(j));
            float DR_Lep0Jet = jet_4vec.DeltaR(lep0_4vec);
            float DR_Lep1Jet = jet_4vec.DeltaR(lep1_4vec);
            if (MinDR_Lep0Jet>DR_Lep0Jet) MinDR_Lep0Jet = DR_Lep0Jet;
            if (MinDR_Lep1Jet>DR_Lep1Jet) MinDR_Lep1Jet = DR_Lep1Jet;
            float DR_PhotonJet = jet_4vec.DeltaR(z_4vec);
            float DPhi_PhotonJet = jet_4vec.DeltaPhi(z_4vec);
            if (MinDR_PhotonJet>DR_PhotonJet) MinDR_PhotonJet = DR_PhotonJet;
            if (MinDPhi_PhotonJet>DPhi_PhotonJet) MinDPhi_PhotonJet = DPhi_PhotonJet;
        }
        float min_DPhi_MET_LepJet = 1000.;
        float DPhi_MET_LepJet = 1000.;
        for (unsigned int j=0;j<jet_pT->size();j++) {
            jet_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),jet_m->at(j));
            DPhi_MET_LepJet = jet_4vec.DeltaR(met_4vec);
            if (min_DPhi_MET_LepJet>DPhi_MET_LepJet) min_DPhi_MET_LepJet = DPhi_MET_LepJet;
        }
        DPhi_MET_LepJet = lep0_4vec.DeltaR(met_4vec);
        if (min_DPhi_MET_LepJet>DPhi_MET_LepJet) min_DPhi_MET_LepJet = DPhi_MET_LepJet;
        DPhi_MET_LepJet = lep1_4vec.DeltaR(met_4vec);
        if (min_DPhi_MET_LepJet>DPhi_MET_LepJet) min_DPhi_MET_LepJet = DPhi_MET_LepJet;
        float MET_rel = MET;
        if (min_DPhi_MET_LepJet<TMath::Pi()/2.) MET_rel = MET*TMath::Sin(min_DPhi_MET_LepJet);

        //--- photon DPhi variables
		TLorentzVector gamma_4vec;
		gamma_4vec.SetPtEtaPhiM(gamma_pt,gamma_eta,gamma_phi,0);
		MinDPhi_PhotonJet = 1000.;
		TLorentzVector jet_4vec;
		for (unsigned int j=0;j<jet_pT->size();j++) {
			jet_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),jet_m->at(j));
			//float DR_PhotonJet = jet_4vec.DeltaR(gamma_4vec);
			float DPhi_PhotonJet = jet_4vec.DeltaPhi(gamma_4vec);
			if (MinDPhi_PhotonJet>DPhi_PhotonJet) MinDPhi_PhotonJet = DPhi_PhotonJet;
		}

        //--- MT
        float MT; BaselineTree->Branch("MT",&MT,"MT/F");
		TLorentzVector met_4vec;
		met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
		if (lep_pT->size()>0) {
            TLorentzVector lep_4vec;
            lep_4vec.SetPtEtaPhiM(lep_pT->at(0),0,lep_phi->at(0),0);  // only transverse component
            MT = (met_4vec+lep_4vec).M();
        }
		else MT = 0;

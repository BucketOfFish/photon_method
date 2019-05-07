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

    //-----------------------------
    // make new branches
    //-----------------------------

    float DR_Wmin2Jet; BaselineTree->Branch("DR_Wmin2Jet", &DR_Wmin2Jet, "DR_Wmin2Jet/F");
    float DR_J0J1; BaselineTree->Branch("DR_J0J1", &DR_J0J1, "DR_J0J1/F");
    float mWmin; BaselineTree->Branch("mWmin", &mWmin, "mWmin/F");
    float Wmin_pt; BaselineTree->Branch("Wmin_pt", &Wmin_pt, "Wmin_pt/F");
    float Wmin_eta; BaselineTree->Branch("Wmin_eta", &Wmin_eta, "Wmin_eta/F");
    float DPhi_METWmin; BaselineTree->Branch("DPhi_METWmin", &DPhi_METWmin, "DPhi_METWmin/F");
    float DPhi_WminZ; BaselineTree->Branch("DPhi_WminZ", &DPhi_WminZ, "DPhi_WminZ/F");
    float mj0j1; BaselineTree->Branch("mj0j1", &mj0j1, "mj0j1/F");
    float W01_pt; BaselineTree->Branch("W01_pt", &W01_pt, "W01_pt/F");
    float DPhi_METW01; BaselineTree->Branch("DPhi_METW01", &DPhi_METW01, "DPhi_METW01/F");
    float DPhi_W01Z; BaselineTree->Branch("DPhi_W01Z", &DPhi_W01Z, "DPhi_W01Z/F");
    float DPhi_METNonWJet; BaselineTree->Branch("DPhi_METNonWJet", &DPhi_METNonWJet, "DPhi_METNonWJet/F");
    float NonWJet_pT; BaselineTree->Branch("NonWJet_pT", &NonWJet_pT, "NonWJet_pT/F");
    float DPhi_METNonWminJet; BaselineTree->Branch("DPhi_METNonWminJet", &DPhi_METNonWminJet, "DPhi_METNonWminJet/F");
    float NonWminJet_pT; BaselineTree->Branch("NonWminJet_pT", &NonWminJet_pT, "NonWminJet_pT/F");
    float DR_2Lep; BaselineTree->Branch("DR_2Lep", &DR_2Lep, "DR_2Lep/F");

    //-----------------------------
    // make histograms
    //-----------------------------

    TH1D* hist_low_njet = new TH1D("hist_low_njet", "", bin_size, njet_bin); hist_low_njet->SetStats(0);
    TH1D* hist_low_nbjet = new TH1D("hist_low_nbjet", "", bin_size, njet_bin); hist_low_nbjet->SetStats(0);
    TH1D* hist_low_pt = new TH1D("hist_low_pt", "", bin_size, sm_pt_bin); hist_low_pt->SetStats(0);
    TH1D* hist_sm_pt = new TH1D("hist_sm_pt", "", bin_size, sm_pt_bin); hist_sm_pt->SetStats(0);
    TH1D* hist_low_et = new TH1D("hist_low_et", "", bin_size, et_bin); hist_low_et->SetStats(0);
    TH1D* hist_low_ht = new TH1D("hist_low_ht", "", bin_size, ht_bin); hist_low_ht->SetStats(0);
    TH1D* hist_medium_njet = new TH1D("hist_medium_njet", "", bin_size, njet_bin); hist_medium_njet->SetStats(0);
    TH1D* hist_medium_nbjet = new TH1D("hist_medium_nbjet", "", bin_size, njet_bin); hist_medium_nbjet->SetStats(0);
    TH1D* hist_medium_pt = new TH1D("hist_medium_pt", "", bin_size, sm_pt_bin); hist_medium_pt->SetStats(0);
    TH1D* hist_medium_et = new TH1D("hist_medium_et", "", bin_size, et_bin); hist_medium_et->SetStats(0);
    TH1D* hist_medium_ht = new TH1D("hist_medium_ht", "", bin_size, ht_bin); hist_medium_ht->SetStats(0);
    TH1D* hist_high_njet = new TH1D("hist_high_njet", "", bin_size, njet_bin); hist_high_njet->SetStats(0);
    TH1D* hist_high_nbjet = new TH1D("hist_high_nbjet", "", bin_size, njet_bin); hist_high_nbjet->SetStats(0);
    TH1D* hist_high_pt = new TH1D("hist_high_pt", "", bin_size, sm_pt_bin); hist_high_pt->SetStats(0);
    TH1D* hist_high_et = new TH1D("hist_high_et", "", bin_size, et_bin); hist_high_et->SetStats(0);
    TH1D* hist_high_ht = new TH1D("hist_high_ht", "", bin_size, ht_bin); hist_high_ht->SetStats(0);
    TH1D* hist_zmet_njet = new TH1D("hist_zmet_njet", "", bin_size, njet_bin); hist_zmet_njet->SetStats(0);
    TH1D* hist_zmet_nbjet = new TH1D("hist_zmet_nbjet", "", bin_size, njet_bin); hist_zmet_nbjet->SetStats(0);
    TH1D* hist_zmet_pt = new TH1D("hist_zmet_pt", "", bin_size, sm_pt_bin); hist_zmet_pt->SetStats(0);
    TH1D* hist_zmet_et = new TH1D("hist_zmet_et", "", bin_size, et_bin); hist_zmet_et->SetStats(0);
    TH1D* hist_zmet_ht = new TH1D("hist_zmet_ht", "", bin_size, ht_bin); hist_zmet_ht->SetStats(0);
    TH1D* hist_bveto_njet = new TH1D("hist_bveto_njet", "", bin_size, njet_bin); hist_bveto_njet->SetStats(0);
    TH1D* hist_bveto_nbjet = new TH1D("hist_bveto_nbjet", "", bin_size, njet_bin); hist_bveto_nbjet->SetStats(0);
    TH1D* hist_bveto_pt = new TH1D("hist_bveto_pt", "", bin_size, sm_pt_bin); hist_bveto_pt->SetStats(0);
    TH1D* hist_bveto_et = new TH1D("hist_bveto_et", "", bin_size, et_bin); hist_bveto_et->SetStats(0);
    TH1D* hist_bveto_ht = new TH1D("hist_bveto_ht", "", bin_size, ht_bin); hist_bveto_ht->SetStats(0);
    TH1D* hist_low_dpt = new TH1D("hist_low_dpt", "", dpt_bin_size, dpt_bin); hist_low_dpt->SetStats(0);
    TH1D* hist_low_pt_smear = new TH1D("hist_low_pt_smear", "", bin_size, sm_pt_bin); hist_low_pt_smear->SetStats(0);
    TH1D* hist_low_htincl = new TH1D("hist_low_htincl", "", bin_size, ht_bin); hist_low_htincl->SetStats(0);
    TH1D* hist_medium_pt_smear = new TH1D("hist_medium_pt_smear", "", bin_size, sm_pt_bin); hist_medium_pt_smear->SetStats(0);
    TH1D* hist_medium_htincl = new TH1D("hist_medium_htincl", "", bin_size, ht_bin); hist_medium_htincl->SetStats(0);
    TH1D* hist_high_pt_smear = new TH1D("hist_high_pt_smear", "", bin_size, sm_pt_bin); hist_high_pt_smear->SetStats(0);
    TH1D* hist_high_htincl = new TH1D("hist_high_htincl", "", bin_size, ht_bin); hist_high_htincl->SetStats(0);
    TH1D* hist_zmet_pt_smear = new TH1D("hist_zmet_pt_smear", "", bin_size, sm_pt_bin); hist_zmet_pt_smear->SetStats(0);
    TH1D* hist_zmet_htincl = new TH1D("hist_zmet_htincl", "", bin_size, ht_bin); hist_zmet_htincl->SetStats(0);
    TH1D* hist_bveto_pt_smear = new TH1D("hist_bveto_pt_smear", "", bin_size, sm_pt_bin); hist_bveto_pt_smear->SetStats(0);
    TH1D* hist_bveto_htincl = new TH1D("hist_bveto_htincl", "", bin_size, ht_bin); hist_bveto_htincl->SetStats(0);
    TH1D* hist_low_met = new TH1D("hist_low_met", "", bin_size, met_bin); hist_low_met->SetStats(0);
    TH1D* hist_low_dphi = new TH1D("hist_low_dphi", "", bin_size, dphi_bin); hist_low_dphi->SetStats(0);

    //-----------------------------
    // dijet
    //-----------------------------

float mj1j2 = 0.;
float DR_J1J2 = 0.;
float DR_Wmin2Jet = 0.;
float DR_J0J1 = 0.;
float DR_W80jj = 0.;
float mWmin = 0.;
float Wmin_pt = 0.;
float Wmin_eta = 0.;
float DPhi_METWmin = 0.;
float DPhi_WminZ = 0.;
float mj0j1 = 0.;
float m80jj = 0.;
float W01_pt = 0.;
float DPhi_METW01 = 0.;
float DPhi_W01Z = 0.;
float DPhi_W80Z = 0.;
float DPhi_METNonWminJet = 0.;
float NonWminJet_pT = 0.;
float W80_pt = 0.;
float DPhi_METW80 = 0.;
std::vector<int>* jet_isW80 = new std::vector<int>(10);
std::vector<int>* jet_isW01 = new std::vector<int>(10);
std::vector<int>* jet_isW12 = new std::vector<int>(10);
std::vector<int>* jet_isWmin = new std::vector<int>(10);
std::vector<int>* jet_isWminJ0 = new std::vector<int>(10);
std::vector<int>* jet_isWminJ1 = new std::vector<int>(10);
float DPhi_METNonWJet = 0.;
float NonWJet_pT = 0.;
float DPhi_METNonW12Jet = 0.;
float NonW12Jet_pT = 0.;
int W_j0 = 0;
int W_j1 = 0;
int Wmin_j0 = 0;
int Wmin_j1 = 0;

void GetDijetVariables(TLorentzVector z_4vec, TLorentzVector met_4vec, std::vector<float>* jet_pT, std::vector<float>* jet_eta, std::vector<float>* jet_phi, std::vector<float>* jet_m) {
	
	TLorentzVector obj_4vec;
	TLorentzVector isr_4vec;
	TLorentzVector jet0_4vec;
	TLorentzVector jet1_4vec;
	TLorentzVector dijet_4vec;
	TLorentzVector zmet_4vec;
	zmet_4vec = z_4vec + met_4vec;
	mj0j1 = -1;
	mj1j2 = -1;
	m80jj = -1;
	DR_W80jj = 0.;
	DR_Wmin2Jet = 0.;
	DR_J0J1 = 0.;
	DR_J1J2 = 0.;
	W80_pt = 0.;
	DPhi_METW80 = 0.;
	DPhi_W80Z = 0.;
	W01_pt = 0.;
	DPhi_METW01 = 0.;
	DPhi_W01Z = 0.;
	mWmin = -1;
	Wmin_pt = 0.;
	Wmin_eta = 0.;
	DPhi_METWmin = 0.;
	DPhi_WminZ = 0.;
	if (jet_pT->size()>1) {
		jet0_4vec.SetPtEtaPhiM(jet_pT->at(0),jet_eta->at(0),jet_phi->at(0),jet_m->at(0));
		jet1_4vec.SetPtEtaPhiM(jet_pT->at(1),jet_eta->at(1),jet_phi->at(1),jet_m->at(1));
		dijet_4vec = jet0_4vec + jet1_4vec;
		mj0j1 = dijet_4vec.M();
		DR_J0J1 = jet0_4vec.DeltaR(jet1_4vec);
		W01_pt = dijet_4vec.Pt();
		DPhi_METW01 = fabs(met_4vec.DeltaPhi(dijet_4vec));
		DPhi_W01Z = fabs(z_4vec.DeltaPhi(dijet_4vec));
	}
	if (jet_pT->size()>2) {
		jet0_4vec.SetPtEtaPhiM(jet_pT->at(1),jet_eta->at(1),jet_phi->at(1),jet_m->at(1));
		jet1_4vec.SetPtEtaPhiM(jet_pT->at(2),jet_eta->at(2),jet_phi->at(2),jet_m->at(2));
		dijet_4vec = jet0_4vec + jet1_4vec;
		mj1j2 = dijet_4vec.M();
		DR_J1J2 = jet0_4vec.DeltaR(jet1_4vec);
	}
	W_j0 = 0;
	W_j1 = 0;
	if (jet_pT->size()>1) {
		jet0_4vec.SetPtEtaPhiM(jet_pT->at(0),jet_eta->at(0),jet_phi->at(0),jet_m->at(0));
		jet1_4vec.SetPtEtaPhiM(jet_pT->at(1),jet_eta->at(1),jet_phi->at(1),jet_m->at(1));
		dijet_4vec = jet0_4vec + jet1_4vec;
		m80jj = dijet_4vec.M();
		W80_pt = dijet_4vec.Pt();
		DR_W80jj = jet0_4vec.DeltaR(jet1_4vec);
		DPhi_METW80 = fabs(met_4vec.DeltaPhi(dijet_4vec));
		DPhi_W80Z = fabs(z_4vec.DeltaPhi(dijet_4vec));
		for (unsigned int j0=0;j0<jet_pT->size();j0++) {
			for (unsigned int j1=j0+1;j1<jet_pT->size();j1++) {
				jet0_4vec.SetPtEtaPhiM(jet_pT->at(j0),jet_eta->at(j0),jet_phi->at(j0),jet_m->at(j0));
				jet1_4vec.SetPtEtaPhiM(jet_pT->at(j1),jet_eta->at(j1),jet_phi->at(j1),jet_m->at(j1));
				dijet_4vec = jet0_4vec + jet1_4vec;
				if (abs(m80jj-80.)>abs(dijet_4vec.M()-80.)) {
					m80jj = dijet_4vec.M();
					W80_pt = dijet_4vec.Pt();
					DR_W80jj = jet0_4vec.DeltaR(jet1_4vec);
					DPhi_METW80 = fabs(met_4vec.DeltaPhi(dijet_4vec));
					DPhi_W80Z = fabs(z_4vec.DeltaPhi(dijet_4vec));
					W_j0 = j0;
					W_j1 = j1;
				}
			}
		}
	}
	Wmin_j0 = 0;	
	Wmin_j1 = 0;
	if (jet_pT->size()>1) {
		double min_dPhi = 1000;
		for (unsigned int j0=0;j0<jet_pT->size();j0++) {
			jet0_4vec.SetPtEtaPhiM(jet_pT->at(j0),jet_eta->at(j0),jet_phi->at(j0),jet_m->at(j0));
			if (min_dPhi>abs(jet0_4vec.DeltaPhi(zmet_4vec))) {
				min_dPhi = abs(jet0_4vec.DeltaPhi(zmet_4vec));
				Wmin_j0 = j0;
			}
		}
		min_dPhi = 1000;
		for (unsigned int j1=0;j1<jet_pT->size();j1++) {
			jet1_4vec.SetPtEtaPhiM(jet_pT->at(j1),jet_eta->at(j1),jet_phi->at(j1),jet_m->at(j1));
			if (min_dPhi>abs(jet1_4vec.DeltaPhi(zmet_4vec))) {
			  if (j1!=(unsigned int) Wmin_j0) {
					min_dPhi = abs(jet1_4vec.DeltaPhi(zmet_4vec));
					Wmin_j1 = j1;
				}
			}
		}
		jet0_4vec.SetPtEtaPhiM(jet_pT->at(Wmin_j0),jet_eta->at(Wmin_j0),jet_phi->at(Wmin_j0),jet_m->at(Wmin_j0));
		jet1_4vec.SetPtEtaPhiM(jet_pT->at(Wmin_j1),jet_eta->at(Wmin_j1),jet_phi->at(Wmin_j1),jet_m->at(Wmin_j1));
		dijet_4vec = jet0_4vec + jet1_4vec;
		mWmin = dijet_4vec.M();
		Wmin_pt = dijet_4vec.Pt();
		Wmin_eta = dijet_4vec.Eta();
		DPhi_METWmin = fabs(met_4vec.DeltaPhi(dijet_4vec));
		DPhi_WminZ = fabs(z_4vec.DeltaPhi(dijet_4vec));
		DR_Wmin2Jet = jet0_4vec.DeltaR(jet1_4vec);
	}
	jet_isW01->clear();
	for (unsigned int j0=0;j0<jet_pT->size();j0++) {
		jet_isW01->push_back(0);
	}
	if (jet_pT->size()>1) {
		jet_isW01->at(0) = 1;
		jet_isW01->at(1) = 1;
	}
	jet_isW12->clear();
	for (unsigned int j0=0;j0<jet_pT->size();j0++) {
		jet_isW12->push_back(0);
	}
	if (jet_pT->size()>2) {
		jet_isW12->at(1) = 1;
		jet_isW12->at(2) = 1;
	}
	jet_isWmin->clear();
	jet_isWminJ0->clear();
	jet_isWminJ1->clear();
	for (unsigned int j0=0;j0<jet_pT->size();j0++) {
		jet_isWmin->push_back(0);
		jet_isWminJ0->push_back(0);
		jet_isWminJ1->push_back(0);
	}
	if (jet_pT->size()>1) {
		jet_isWmin->at(Wmin_j0) = 1;
		jet_isWminJ0->at(Wmin_j0) = 1;
		jet_isWmin->at(Wmin_j1) = 1;
		jet_isWminJ1->at(Wmin_j1) = 1;
	}
	isr_4vec.SetPtEtaPhiM(0,0,0,0);
	for (unsigned int j=0;j<jet_pT->size();j++) {
		if (jet_isW12->at(j)==0) {
			obj_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),jet_m->at(j));
			isr_4vec += obj_4vec;
		}
	}
	DPhi_METNonW12Jet = fabs(met_4vec.DeltaPhi(isr_4vec));
	NonW12Jet_pT = isr_4vec.Pt();
	isr_4vec.SetPtEtaPhiM(0,0,0,0);
	for (unsigned int j=0;j<jet_pT->size();j++) {
		if (jet_isWmin->at(j)==0) {
			obj_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),jet_m->at(j));
			isr_4vec += obj_4vec;
		}
	}
	DPhi_METNonWminJet = fabs(met_4vec.DeltaPhi(isr_4vec));
	NonWminJet_pT = isr_4vec.Pt();
	jet_isW80->clear();
	for (unsigned int j0=0;j0<jet_pT->size();j0++) {
		jet_isW80->push_back(0);
	}
	if (jet_pT->size()>1) {
		jet_isW80->at(W_j0) = 1;
		jet_isW80->at(W_j1) = 1;
	}
	isr_4vec.SetPtEtaPhiM(0,0,0,0);
	for (unsigned int j=0;j<jet_pT->size();j++) {
		if (jet_isW01->at(j)==0) {
			obj_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),jet_m->at(j));
			isr_4vec += obj_4vec;
		}
	}
	DPhi_METNonWJet = fabs(met_4vec.DeltaPhi(isr_4vec));
	NonWJet_pT = isr_4vec.Pt();

}

Int_t EventNumber;
Int_t RunNumber;
double totalWeight = 0.;
int pt = 0; // pT bin for reweighting
int smpt = 0; // pT bin for smearing
int ht = 0;
int njet = 0;
int nbjet = 0;
Int_t jet_n;
Int_t bjet_n;
float HT = 0.;
float truthGamma_pt = 0.;
float truthGamma_phi = 0.;
float truthGamma_eta = 0.;
Float_t gamma_pt = 0.;
Float_t gamma_ht = 0.;
float gamma_eta = 0.;
float gamma_phi = 0.;
float gamma_dR = 0.;
float MET;
//-------------------------------
//Variables for 2019 RJR analysis
//-------------------------------
float MET_loose;
float MET_tight;
float MET_tighter;
float MET_tenacious;
Int_t trigMatch_2LTrigOR;
Int_t is2Lep2Jet;
Int_t is2L2JInt;
int nBJet20_MV2c10_FixedCutBEff_77;
float mjj;
Float_t mll_RJ;
Float_t R_minH2P_minH3P;
Float_t RPT_HT5PP;
Float_t dphiVP;
Float_t H2PP;
Float_t H5PP;
int nJet20;
Float_t minDphi;
Float_t MZ;
Int_t NjS;
Int_t NjISR;
Float_t dphiISRI;
Float_t RISR;
Float_t PTISR;
Float_t PTI;
Float_t PTCM;
Float_t MJ;
Int_t is3Lep3Jet;
Int_t is4Lep3Jet;
Int_t lept1sign_VR;
Int_t lept2sign_VR;
Float_t lept1Pt_VR;
Float_t lept2Pt_VR;
Float_t MZ_VR;
Float_t MJ_VR;
Float_t RISR_VR;
Float_t PTISR_VR;
Float_t PTI_VR;
Float_t PTCM_VR;
Float_t dphiISRI_VR;
//-------------------------------------
float MET_phi;
float DPhi_METJetLeading;
float DPhi_METJetSecond;
float MinDPhi_PhotonJet;

// new variables to be added to ntuple
int lep_n = 0;
int dpt = 0;
int pt_smear = 0;
int met_smear = 0;
int dphi_smear = 0;
float HTincl = 0.;
Float_t mll = 0.;
float smear_shift = 0.;
float gamma_pt_smear = 0.;
float gamma_ht_smear = 0.;
float gamma_phi_smear = 0.;
float MET_smear;
float MET_phi_smear;
float DPhi_METJetLeading_smear;
float DPhi_METJetSecond_smear;
float DPhi_METLepLeading_smear;
float DPhi_METLepSecond_smear;
float METl_smear = 0.;
float METt_smear = 0.;
float DPhi_METPhoton_smear = 0.;
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

float MT2W;

int W_j0 = 0;
int W_j1 = 0;
int Wmin_j0 = 0;
int Wmin_j1 = 0;

int    channel = 0;
int    is_OS   = 0;
float  Z_eta   = 0.0;
float  Z_phi   = 0.0;
int    el_n    = 0;
int    mu_n    = 0;
float MET_rel = 0.;
float METl = 0.;
float METt = 0.;
float DPhi_2Lep = 0.;
float DR_2Lep = 0.;
float DPhi_METPhoton = 0.;
float DPhi_METLepLeading = 0.;
float DPhi_METLepSecond = 0.;
float DPhi_METLepMin = 0.;
float DPhi_TSTLepLeading = 0.;
float DPhi_TSTLepSecond = 0.;
float DPhi_TSTLepMin = 0.;
float MinDR_Lep0Jet = 0.;
float MinDR_Lep1Jet = 0.;
float MinDR_PhotonJet = 0.;
float FS_ee_weight = 0.;
float FS_mm_weight = 0.;
std::vector<int>* jet_isPrompt = new std::vector<int>(10);
// following are variables from Jigsaw
std::vector<int>* lep_MetHsph = new std::vector<int>(10);
std::vector<int>* jet_MetHsph = new std::vector<int>(10);
float DPhi_METISR = 0.;
float ISR_pT = 0.;
float ISR_eta = 0.;
float ISR_phi = 0.;
float Jigsaw_ZMass = 0.;
float Jigsaw_WMass = 0.;
float DPhi_METJigsawW = 0.;

void GetDijetVariables(TLorentzVector z_4vec, TLorentzVector met_4vec) {
	
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

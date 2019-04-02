
int W_j0 = 0;
int W_j1 = 0;
int Wmin_j0 = 0;
int Wmin_j1 = 0;
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
	//if (jet_pT->size()>1) {
	//	double min_dPhi = 0;
	//	isr_4vec.SetPtEtaPhiM(jet_pT->at(0),jet_eta->at(0),jet_phi->at(0),jet_m->at(0));
	//	for (int j0=1;j0<jet_pT->size();j0++) {
	//		jet0_4vec.SetPtEtaPhiM(jet_pT->at(j0),jet_eta->at(j0),jet_phi->at(j0),jet_m->at(j0));
	//		if (min_dPhi<abs(jet0_4vec.DeltaPhi(isr_4vec))) {
	//			min_dPhi = abs(jet0_4vec.DeltaPhi(isr_4vec));
	//			Wmin_j0 = j0;
	//		}
	//	}
	//	min_dPhi = 0;
	//	for (int j1=0;j1<jet_pT->size();j1++) {
	//		jet1_4vec.SetPtEtaPhiM(jet_pT->at(j1),jet_eta->at(j1),jet_phi->at(j1),jet_m->at(j1));
	//		if (min_dPhi<abs(jet1_4vec.DeltaPhi(isr_4vec))) {
	//			if (j1!=Wmin_j0) {
	//				min_dPhi = abs(jet1_4vec.DeltaPhi(isr_4vec));
	//				Wmin_j1 = j1;
	//			}
	//		}
	//	}
	//	jet0_4vec.SetPtEtaPhiM(jet_pT->at(Wmin_j0),jet_eta->at(Wmin_j0),jet_phi->at(Wmin_j0),jet_m->at(Wmin_j0));
	//	jet1_4vec.SetPtEtaPhiM(jet_pT->at(Wmin_j1),jet_eta->at(Wmin_j1),jet_phi->at(Wmin_j1),jet_m->at(Wmin_j1));
	//	dijet_4vec = jet0_4vec + jet1_4vec;
	//	mWmin = dijet_4vec.M();
	//	Wmin_pt = dijet_4vec.Pt();
	//	Wmin_eta = dijet_4vec.Eta();
	//	DPhi_METWmin = fabs(met_4vec.DeltaPhi(dijet_4vec));
	//	DPhi_WminZ = fabs(z_4vec.DeltaPhi(dijet_4vec));
	//	DR_Wmin2Jet = jet0_4vec.DeltaR(jet1_4vec);
	//	//if (abs(m80jj-80.)<10.) {
	//	//	Wmin_j0 = W_j0;
	//	//	Wmin_j1 = W_j1;
	//	//	mWmin = m80jj;
	//	//	Wmin_pt = W80_pt;
	//	//	DPhi_METWmin = DPhi_METW80;
	//	//	DPhi_WminZ = DPhi_W80Z;
	//	//}
	//}
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

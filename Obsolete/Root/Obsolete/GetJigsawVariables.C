#include "GetNThingsFromAGroup.C"
void GetJigsawVariables() {

	TLorentzVector obj_4vec;
	TLorentzVector vis_4vec;
	TLorentzVector inv_4vec;
	TLorentzVector isr_4vec;
	TLorentzVector jigsaw_z_4vec;
	TLorentzVector jigsaw_w_4vec;
	TLorentzVector thing1_4vec;
	TLorentzVector thing2_4vec;
	std::vector<int>* obj_flavor = new std::vector<int>(10);
	std::vector<int>* obj_idx = new std::vector<int>(10);
	TLorentzVector met_4vec;
	met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);

	lep_MetHsph->clear();
	jet_MetHsph->clear();
	for (int j1=0;j1<lep_pT->size();j1++) {
		lep_MetHsph->push_back(0);
	}
	for (int j1=0;j1<jet_pT->size();j1++) {
		jet_MetHsph->push_back(0);
	}
	obj_flavor->clear();
	obj_idx->clear();
	for (int l=0;l<lep_pT->size();l++) {
		obj_flavor->push_back(0);
		obj_idx->push_back(l);
	}
	for (int l=0;l<jet_pT->size();l++) {
		obj_flavor->push_back(1);
		obj_idx->push_back(l);
	}

	double minMass = 1E6;
	double maxPt = 0;
	Jigsaw_ZMass = 1E6;
	Jigsaw_WMass = 1E6;
	DPhi_METJigsawW = 1E6;
	vis_4vec.SetPtEtaPhiM(0,0,0,0);
	inv_4vec.SetPtEtaPhiM(0,0,0,0);
	for (int j1=0;j1<lep_pT->size();j1++) {
		obj_4vec.SetPtEtaPhiM(lep_pT->at(j1),0,lep_phi->at(j1),0);
		vis_4vec += obj_4vec;
	}
	for (int j1=0;j1<jet_pT->size();j1++) {
		obj_4vec.SetPtEtaPhiM(jet_pT->at(j1),0,jet_phi->at(j1),0);
		vis_4vec += obj_4vec;
	}
	int num_considered_obj = jet_pT->size()+lep_pT->size();
	if (num_considered_obj>9) num_considered_obj = 9;
	for (int num=0;num<=num_considered_obj;num++) {
		vector< vector<int> > not_methsph_list = GetNThingsFromAGroup(num,jet_pT->size()+lep_pT->size());
		for (int partition=0;partition<not_methsph_list.size();partition++) {
			// make sure leading lep is not in the ISR list
			int nlep_in_list = 0;
			int leadlep_in_list = 0;
			for (int obj=0;obj<not_methsph_list.at(partition).size();obj++) {
				int idx = obj_idx->at(not_methsph_list.at(partition).at(obj));
				int flavor = obj_flavor->at(not_methsph_list.at(partition).at(obj));
				if (flavor==0) nlep_in_list += 1;
				if (flavor==0 && idx==0) leadlep_in_list += 1;
			}
			//if (nlep_in_list>0) continue; // all lep are not in ISR list
			//if (leadlep_in_list>0) continue; // leading lep not in ISR list
			thing1_4vec.SetPtEtaPhiM(0,0,0,0);  // SUSY tree
			thing2_4vec.SetPtEtaPhiM(0,0,0,0);  // ISR tree
			for (int j1=0;j1<obj_idx->size();j1++) {
				bool in_not_methsph_list = false;
				int idx = obj_idx->at(j1);
				int flavor = obj_flavor->at(j1);
				if (flavor==0) obj_4vec.SetPtEtaPhiM(lep_pT->at(idx),lep_eta->at(idx),lep_phi->at(idx),0);
				else obj_4vec.SetPtEtaPhiM(jet_pT->at(idx),jet_eta->at(idx),jet_phi->at(idx),0);
				for (int obj=0;obj<not_methsph_list.at(partition).size();obj++) {
					if (j1==not_methsph_list.at(partition).at(obj)) {
						in_not_methsph_list = true;
						thing2_4vec += obj_4vec;
					}
				}
				if (in_not_methsph_list==false) {
					thing1_4vec += obj_4vec;
				}
			}
			inv_4vec.SetPtEtaPhiM(MET,-vis_4vec.Eta(),MET_phi,vis_4vec.M());
			//inv_4vec.SetPtEtaPhiM(MET,-vis_4vec.Eta(),MET_phi,thing1_4vec.M());
			thing1_4vec += inv_4vec;

			double this_mass = thing1_4vec.M() + thing2_4vec.M();
			double this_pt = thing1_4vec.Pt() + thing2_4vec.Pt();
			//double this_mass = pow(thing1_4vec.M(),2) + pow(thing2_4vec.M(),2);
			//this_mass = pow(this_mass,0.5);
			//double this_pt = pow(thing1_4vec.Pt(),2) + pow(thing2_4vec.Pt(),2);
			//this_pt = pow(this_pt,0.5);
			
			if (this_mass<minMass) {
				minMass = this_mass;
			//if (this_pt>maxPt) {
			//	maxPt = this_pt;
				lep_MetHsph->clear();
				jet_MetHsph->clear();
				for (int j1=0;j1<lep_pT->size();j1++) {
					lep_MetHsph->push_back(1);
				}
				for (int j1=0;j1<jet_pT->size();j1++) {
					jet_MetHsph->push_back(1);
				}
				for (int obj=0;obj<not_methsph_list.at(partition).size();obj++) {
					int j1 = not_methsph_list.at(partition).at(obj);
					int idx = obj_idx->at(j1);
					int flavor = obj_flavor->at(j1);
					if (flavor==0) {
						lep_MetHsph->at(idx) = 0;
					}
					else {
						jet_MetHsph->at(idx) = 0;
					}
				}
			}
		}
	}
	int n_methsph_lep = 0;
	int n_methsph_jet = 0;
	isr_4vec.SetPtEtaPhiM(0,0,0,0);
	jigsaw_z_4vec.SetPtEtaPhiM(0,0,0,0);
	jigsaw_w_4vec.SetPtEtaPhiM(0,0,0,0);
	for (int l=0;l<lep_MetHsph->size();l++) {
		if (lep_MetHsph->at(l)==1) {
			n_methsph_lep += 1;
			obj_4vec.SetPtEtaPhiM(lep_pT->at(l),lep_eta->at(l),lep_phi->at(l),0);
			jigsaw_z_4vec += obj_4vec;
		}
		else {
			obj_4vec.SetPtEtaPhiM(lep_pT->at(l),lep_eta->at(l),lep_phi->at(l),0);
			isr_4vec += obj_4vec;
		}
	}
	for (int j=0;j<jet_MetHsph->size();j++) {
		if (jet_MetHsph->at(j)==1) {
			n_methsph_jet += 1;
			obj_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),0);
			jigsaw_w_4vec += obj_4vec;
		}
		else {
			obj_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),0);
			isr_4vec += obj_4vec;
		}
	}
	DPhi_METISR = fabs(met_4vec.DeltaPhi(isr_4vec));
	ISR_pT = isr_4vec.Pt();
	ISR_eta = isr_4vec.Eta();
	ISR_phi = isr_4vec.Phi();
	Jigsaw_ZMass = jigsaw_z_4vec.M();
	Jigsaw_WMass = jigsaw_w_4vec.M();
	DPhi_METJigsawW = fabs(met_4vec.DeltaPhi(jigsaw_w_4vec));


}

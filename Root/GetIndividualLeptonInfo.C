
#include "TRandom.h"

TRandom myRandom;

void GetIndividualLeptonInfo(TLorentzVector z_4vec) {

	// Compute two lepton pT in Z boson rest (C.M.) frame.
	// lepton p = mll/2. in this frame.
	// mll is derived from the smearing code.
	// two leptons are back-to-back.
	// phi is taken randomly from 0-2pi
	// theta is taken randomly from 0-pi
	TLorentzVector lep0_cm_4vec;
	TLorentzVector lep1_cm_4vec;
	TVector3 boost_vec;
	TLorentzVector lep0_lab_4vec;
	TLorentzVector lep1_lab_4vec;
	if (z_4vec.M()>0.) {
		double lep_E_cm = z_4vec.M()/2.;
		double lep_phi_cm = myRandom.Rndm()*2.*TMath::Pi();
		double lep_theta_cm = myRandom.Rndm()*TMath::Pi()-0.5*TMath::Pi();
		lep0_cm_4vec.SetPxPyPzE(lep_E_cm*TMath::Cos(lep_theta_cm)*TMath::Cos(lep_phi_cm),lep_E_cm*TMath::Cos(lep_theta_cm)*TMath::Sin(lep_phi_cm),lep_E_cm*TMath::Sin(lep_theta_cm),lep_E_cm);
		lep1_cm_4vec.SetPxPyPzE(-lep_E_cm*TMath::Cos(lep_theta_cm)*TMath::Cos(lep_phi_cm),-lep_E_cm*TMath::Cos(lep_theta_cm)*TMath::Sin(lep_phi_cm),-lep_E_cm*TMath::Sin(lep_theta_cm),lep_E_cm);
	}
	// Now we boost the lepton vectors in the rest frame back to the Lab frame.
	// We use smeared photon pT, eta and phi, coded in z_4vec
	boost_vec = z_4vec.BoostVector();
	lep0_lab_4vec = lep0_cm_4vec;
	lep1_lab_4vec = lep1_cm_4vec;
	lep0_lab_4vec.Boost(boost_vec);
	lep1_lab_4vec.Boost(boost_vec);
	lep_pT->clear();
	lep_eta->clear();
	lep_phi->clear();
	if (lep0_lab_4vec.Pt()>lep1_lab_4vec.Pt()) {
		lep_pT->push_back(lep0_lab_4vec.Pt());
		lep_eta->push_back(lep0_lab_4vec.Eta());
		lep_phi->push_back(lep0_lab_4vec.Phi());
		lep_pT->push_back(lep1_lab_4vec.Pt());
		lep_eta->push_back(lep1_lab_4vec.Eta());
		lep_phi->push_back(lep1_lab_4vec.Phi());
	}
	else {
		lep_pT->push_back(lep1_lab_4vec.Pt());
		lep_eta->push_back(lep1_lab_4vec.Eta());
		lep_phi->push_back(lep1_lab_4vec.Phi());
		lep_pT->push_back(lep0_lab_4vec.Pt());
		lep_eta->push_back(lep0_lab_4vec.Eta());
		lep_phi->push_back(lep0_lab_4vec.Phi());
	}
	// sanity check
	//TLorentzVector twolep_cm_4vec;
	//twolep_cm_4vec = lep0_cm_4vec + lep1_cm_4vec;
	//TLorentzVector twolep_lab_4vec;
	//twolep_lab_4vec = lep0_lab_4vec + lep1_lab_4vec;
	//std::cout << "z_4vec pT = " << z_4vec.Pt() << ", eta = " << z_4vec.Eta() << ", phi = " << z_4vec.Phi() << ", m = " << z_4vec.M() << std::endl;
	//std::cout << "lep0_cm_4vec pT = " << lep0_cm_4vec.Pt() << ", eta = " << lep0_cm_4vec.Eta() << ", phi = " << lep0_cm_4vec.Phi() << std::endl;
	//std::cout << "lep1_cm_4vec pT = " << lep1_cm_4vec.Pt() << ", eta = " << lep1_cm_4vec.Eta() << ", phi = " << lep1_cm_4vec.Phi() << std::endl;
	//std::cout << "2lep_cm_4vec pT = " << twolep_cm_4vec.Pt() << ", eta = " << twolep_cm_4vec.Eta() << ", phi = " << twolep_cm_4vec.Phi() << ", M = " << twolep_cm_4vec.M() << std::endl;
	//std::cout << "lep0_lab_4vec pT = " << lep0_lab_4vec.Pt() << ", eta = " << lep0_lab_4vec.Eta() << ", phi = " << lep0_lab_4vec.Phi() << std::endl;
	//std::cout << "lep1_lab_4vec pT = " << lep1_lab_4vec.Pt() << ", eta = " << lep1_lab_4vec.Eta() << ", phi = " << lep1_lab_4vec.Phi() << std::endl;
	//std::cout << "2lep_lab_4vec pT = " << twolep_lab_4vec.Pt() << ", eta = " << twolep_lab_4vec.Eta() << ", phi = " << twolep_lab_4vec.Phi() << ", M = " << twolep_lab_4vec.M() << std::endl;
	//std::cout << "==================================================================================" << std::endl;

}

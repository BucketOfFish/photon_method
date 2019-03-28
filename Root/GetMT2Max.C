
#include "TRandom.h"

TRandom myRandom;

double GetMT2LorentzTransform(TLorentzVector boost_4vec, TLorentzVector lep0_4vec, TLorentzVector lep1_4vec, TLorentzVector met_4vec) {

	TVector3 boost_vec = boost_4vec.BoostVector();
	lep0_4vec.Boost(boost_vec);
	lep1_4vec.Boost(boost_vec);
	met_4vec.Boost(boost_vec);
	double MT2_new = ComputeMT2(lep0_4vec, lep1_4vec, met_4vec, 0, 0).Compute();
	//int binx = hist_mt2max_map->GetXaxis()->FindBin(boost_4vec.Pt());
	//int biny = hist_mt2max_map->GetYaxis()->FindBin(boost_4vec.Phi());
	//hist_mt2max_map->SetBinContent(binx,biny,max(MT2_new,hist_mt2max_map->GetBinContent(binx,biny)));
	//hist_mt2max_beta->SetBinContent(binx,max(MT2_new,hist_mt2max_beta->GetBinContent(binx)));
	return MT2_new;

}
double GetMT2Max(double inv_mass, TLorentzVector lep0_4vec, TLorentzVector lep1_4vec, TLorentzVector met_4vec, double *boost_pt_final, double *boost_phi_final) {

	lep0_4vec.SetPtEtaPhiM(lep0_4vec.Pt(),lep0_4vec.Eta(),lep0_4vec.Phi(),0.);
	lep1_4vec.SetPtEtaPhiM(lep1_4vec.Pt(),lep1_4vec.Eta(),lep1_4vec.Phi(),0.);
	TLorentzVector boost_4vec;
	double boost_pt = 0;
	double boost_phi = 0;
	double MT2_this = 0;
	double MT2_max = ComputeMT2(lep0_4vec, lep1_4vec, met_4vec, 0, 0).Compute();
	//for (int i=0;i<1000;i++) {
	//	boost_phi = myRandom.Rndm()*2.*TMath::Pi()-TMath::Pi()+met_4vec.Phi();
	//	for (int j=0;j<1000;j++) {
	//		boost_pt = myRandom.Rndm()*beta_limit;
	//		boost_4vec.SetPtEtaPhiM(boost_pt,0.,boost_phi,1.);
	//		double dm = (lep0_4vec+lep1_4vec).M()*(1.0*myRandom.Rndm()-0.5);
	//		met_4vec.SetPtEtaPhiM(met_4vec.Pt(),-(lep0_4vec+lep1_4vec).Eta(),met_4vec.Phi(),inv_mass);
	//		MT2_this = GetMT2LorentzTransform(boost_4vec, lep0_4vec, lep1_4vec, met_4vec);
	//		if (MT2_max<MT2_this) {
	//			MT2_max = MT2_this;
	//			*boost_phi_final = boost_phi-met_4vec.Phi();
	//			*boost_pt_final = boost_pt;
	//		}
	//	}
	//}
	inv_mass = (lep0_4vec+lep1_4vec).M();
	for (int i=0;i<100;i++) {
		boost_phi = myRandom.Rndm()*2.*TMath::Pi()-TMath::Pi()+met_4vec.Phi();
		boost_4vec.SetPtEtaPhiM(beta_limit/2.,0.,boost_phi,1.);
		met_4vec.SetPtEtaPhiM(met_4vec.Pt(),-(lep0_4vec+lep1_4vec).Eta(),met_4vec.Phi(),inv_mass);
		MT2_this = GetMT2LorentzTransform(boost_4vec, lep0_4vec, lep1_4vec, met_4vec);
		if (MT2_max<MT2_this) {
			MT2_max = MT2_this;
			*boost_phi_final = boost_phi;
		}
	}
	for (int i=0;i<100;i++) {
		boost_pt = myRandom.Rndm()*beta_limit;
		boost_4vec.SetPtEtaPhiM(boost_pt,0.,*boost_phi_final,1.);
		met_4vec.SetPtEtaPhiM(met_4vec.Pt(),-(lep0_4vec+lep1_4vec).Eta(),met_4vec.Phi(),inv_mass);
		MT2_this = GetMT2LorentzTransform(boost_4vec, lep0_4vec, lep1_4vec, met_4vec);
		if (MT2_max<MT2_this) {
			MT2_max = MT2_this;
			*boost_pt_final = boost_pt;
		}
	}
	return MT2_max;

}

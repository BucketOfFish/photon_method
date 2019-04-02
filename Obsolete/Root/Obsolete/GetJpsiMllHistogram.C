TH1D* hist_Mll_dPt[bin_size][dpt_bin_size];
TH1D* hist_low_dpt = new TH1D("hist_low_dpt","",dpt_bin_size,dpt_bin);
TH1D* hist_sm_pt = new TH1D("hist_sm_pt","",bin_size,bphys_pt_bin);

void GetJpsiMllHistogram(string ch) {

	double Z_truthPt_dilep = 0;
	//for (int bin=0;bin<dpt_bin_size;bin++) {
	//	hist_Mll_dPt[bin] = new TH1D(TString("hist_Mll_dPt_")+TString::Itoa(bin,10),"",mll_bin_size,mll_bin);
	//}
	for (int bin0=0;bin0<bin_size;bin0++) {
		for (int bin1=0;bin1<dpt_bin_size;bin1++) {
			hist_Mll_dPt[bin0][bin1] = new TH1D(TString("hist_Mll_dPt_")+TString::Itoa(bin0,10)+TString("_")+TString::Itoa(bin1,10),"",mll_bin_size,mll_bin);
		}
	}
	TFile fZ( TString(smearingPath)+"zjets_221/zjets_"+TString(ch)+".root" );
	TTree*  tZ              = (TTree*)fZ.Get("BaselineTree");
	tZ->SetBranchStatus("*", 0);
	tZ->SetBranchStatus("totalWeight", 1);
	tZ->SetBranchStatus("jet_n", 1);
	tZ->SetBranchStatus("bjet_n", 1);
	tZ->SetBranchStatus("Z_pt", 1);
	tZ->SetBranchStatus("Z_truthPt", 1);
	tZ->SetBranchStatus("HT", 1);
	tZ->SetBranchStatus("METl", 1);
	tZ->SetBranchStatus("mll", 1);
	tZ->SetBranchStatus("lep_pT" ,1);
	tZ->SetBranchAddress("totalWeight" ,&totalWeight);
	tZ->SetBranchAddress("jet_n" ,&jet_n);
	tZ->SetBranchAddress("bjet_n" ,&bjet_n);
	tZ->SetBranchAddress("Z_pt" ,&Z_pt);
	tZ->SetBranchAddress("Z_truthPt" ,&Z_truthPt);
	tZ->SetBranchAddress("HT" ,&HT);
	tZ->SetBranchAddress("METl" ,&METl);
	tZ->SetBranchAddress("mll" ,&mll);
	tZ->SetBranchAddress("lep_pT" ,&lep_pT           );
	for (int entry=0;entry<tZ->GetEntries();entry++) {
		tZ->GetEntry(entry);
		if (jet_n==0) continue;
		if (Z_pt>50) continue;
		if (lep_pT->at(0)<leading_lep_pt_cut) continue;
		if (lep_pT->at(1)<second_lep_pt_cut) continue;
		//int pt = hist_low_pt->FindBin(Z_pt)-1;
		int pt = hist_sm_pt->FindBin(Z_pt)-1;
		//int dpt = hist_low_dpt->FindBin((Z_pt-Z_truthPt))-1;
		int dpt = hist_low_dpt->FindBin(METl)-1;
		if (dpt>=0 && pt>=0) hist_Mll_dPt[pt][dpt]->Fill(mll,totalWeight);
	}
	fZ.Close();

}

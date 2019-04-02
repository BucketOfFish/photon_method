TH1D* z_dpt[bin_size];
TH1D* z_dphi[bin_size];
TH1D* g_dpt[bin_size];
TH1D* z_metl[bin_size];
TH1D* z_metl_2j[bin_size];
TH1D* tt_metl[bin_size];
TH1D* vv_metl[bin_size];
TH1D* g_metl[bin_size];
TH1D* z_jetmetl[bin_size];
TH1D* tt_jetmetl[bin_size];
TH1D* vv_jetmetl[bin_size];

TH1D* hist_low_pt = new TH1D("hist_low_pt","",bin_size,bphys_sm_pt_bin);

//double totalWeight = 0.;
//int jet_n = 0;
//int bjet_n = 0;
double Z_pt = 0.;
double Z_truthPt = 0.;
//double gamma_pt = 0.;
//double HT = 0.;
//double mll = 0.;
//double METl = 0.;

void GetJpsiSmearingHistogram(string ch, double lumi, string photon_tag) {

	for (int bin=0;bin<bin_size;bin++) {
		z_dpt[bin] = new TH1D(TString("z_dpt_")+TString::Itoa(bin,10),"",40000,-30000,10000);
		g_dpt[bin] = new TH1D(TString("g_dpt_")+TString::Itoa(bin,10),"",40000,-30000,10000);
		z_metl[bin] = new TH1D(TString("z_metl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
		z_metl_2j[bin] = new TH1D(TString("z_metl_2j_")+TString::Itoa(bin,10),"",40000,-30000,10000);
		z_jetmetl[bin] = new TH1D(TString("z_jetmetl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
		g_metl[bin] = new TH1D(TString("g_metl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
	}
	if (smearing_method == 2) {  // data-driven smearing function
		std::cout << "Get smearing function from data." << std::endl;
		cout << "Opening file           : " << TString(smearingPath)+"data/data_"+TString(ch)+".root"  << endl;
		TFile fZ( TString(smearingPath)+"data/data_"+TString(ch)+".root" );
		TTree*  tZ              = (TTree*)fZ.Get("BaselineTree");
		tZ->SetBranchStatus("*", 0);
		tZ->SetBranchStatus("totalWeight", 1);
		tZ->SetBranchStatus("jet_n", 1);
		tZ->SetBranchStatus("bjet_n", 1);
		tZ->SetBranchStatus("Z_pt", 1);
		tZ->SetBranchStatus("HT", 1);
		tZ->SetBranchStatus("mll", 1);
		tZ->SetBranchStatus("METl", 1);
		tZ->SetBranchAddress("totalWeight" ,&totalWeight);
		tZ->SetBranchAddress("jet_n" ,&jet_n);
		tZ->SetBranchAddress("bjet_n" ,&bjet_n);
		tZ->SetBranchAddress("Z_pt" ,&Z_pt);
		tZ->SetBranchAddress("HT" ,&HT);
		tZ->SetBranchAddress("mll" ,&mll);
		tZ->SetBranchAddress("METl" ,&METl);
		for (int entry=0;entry<tZ->GetEntries();entry++) {
			tZ->GetEntry(entry);
			if (Z_pt<50.) continue;
			int pt = hist_low_pt->FindBin(Z_pt)-1;
			if (jet_n!=1) continue;
			if (bjet_n!=0) continue;
			z_metl[pt]->Fill(METl,totalWeight);
			if (mll<90 || mll>92) continue;
			z_jetmetl[pt]->Fill(METl,totalWeight);
		}
		fZ.Close();
		//cout << "Opening file           : " << TString(smearingPath)+"tt/tt"+TString(ch)+".root"        << endl;
		//TFile ftt( TString(smearingPath)+"tt/tt"+TString(ch)+".root" );
		//TTree*  ttt              = (TTree*)ftt.Get("BaselineTree");
		//ttt->SetBranchStatus("*", 0);
		//ttt->SetBranchStatus("totalWeight", 1);
		//ttt->SetBranchStatus("jet_n", 1);
		//ttt->SetBranchStatus("bjet_n", 1);
		//ttt->SetBranchStatus("Z_pt", 1);
		//ttt->SetBranchStatus("HT", 1);
		//ttt->SetBranchStatus("mll", 1);
		//ttt->SetBranchStatus("METl", 1);
		//ttt->SetBranchAddress("totalWeight" ,&totalWeight);
		//ttt->SetBranchAddress("jet_n" ,&jet_n);
		//ttt->SetBranchAddress("bjet_n" ,&bjet_n);
		//ttt->SetBranchAddress("Z_pt" ,&Z_pt);
		//ttt->SetBranchAddress("HT" ,&HT);
		//ttt->SetBranchAddress("mll" ,&mll);
		//ttt->SetBranchAddress("METl" ,&METl);
		//for (int entry=0;entry<ttt->GetEntries();entry++) {
		//	ttt->GetEntry(entry);
		//	if (Z_pt<50.) continue;
		//	int pt = hist_low_pt->FindBin(Z_pt)-1;
		//	if (jet_n!=1) continue;
		//	if (bjet_n!=0) continue;
		//	z_metl[pt]->Fill(METl,-1.*lumi*totalWeight);
		//	if (mll<90 || mll>92) continue;
		//	z_jetmetl[pt]->Fill(METl,-1.*lumi*totalWeight);
		//}
		//ftt.Close();
		//cout << "Opening file           : " << TString(smearingPath)+"vv/vv"+TString(ch)+".root"       << endl;
		//TFile fvv( TString(smearingPath)+"vv/vv"+TString(ch)+".root" );
		//TTree*  tvv              = (TTree*)fvv.Get("BaselineTree");
		//tvv->SetBranchStatus("*", 0);
		//tvv->SetBranchStatus("totalWeight", 1);
		//tvv->SetBranchStatus("jet_n", 1);
		//tvv->SetBranchStatus("bjet_n", 1);
		//tvv->SetBranchStatus("Z_pt", 1);
		//tvv->SetBranchStatus("HT", 1);
		//tvv->SetBranchStatus("mll", 1);
		//tvv->SetBranchStatus("METl", 1);
		//tvv->SetBranchAddress("totalWeight" ,&totalWeight);
		//tvv->SetBranchAddress("jet_n" ,&jet_n);
		//tvv->SetBranchAddress("bjet_n" ,&bjet_n);
		//tvv->SetBranchAddress("Z_pt" ,&Z_pt);
		//tvv->SetBranchAddress("HT" ,&HT);
		//tvv->SetBranchAddress("mll" ,&mll);
		//tvv->SetBranchAddress("METl" ,&METl);
		//for (int entry=0;entry<tvv->GetEntries();entry++) {
		//	tvv->GetEntry(entry);
		//	if (Z_pt<50.) continue;
		//	int pt = hist_low_pt->FindBin(Z_pt)-1;
		//	if (jet_n!=1) continue;
		//	if (bjet_n!=0) continue;
		//	z_metl[pt]->Fill(METl,-1.*lumi*totalWeight);
		//	if (mll<90 || mll>92) continue;
		//	z_jetmetl[pt]->Fill(METl,-1.*lumi*totalWeight);
		//}
		//fvv.Close();
		TFile fPhoton( TString(smearingPath)+"bphys/bdata.root" );
		TTree*  tPhoton              = (TTree*)fPhoton.Get("BaselineTree");
		tPhoton->SetBranchStatus("*", 0);
		tPhoton->SetBranchStatus("totalWeight", 1);
		tPhoton->SetBranchStatus("jet_n", 1);
		tPhoton->SetBranchStatus("bjet_n", 1);
		tPhoton->SetBranchStatus("gamma_pt", 1);
		tPhoton->SetBranchStatus("HT", 1);
		tPhoton->SetBranchStatus("METl_raw", 1);
		tPhoton->SetBranchAddress("totalWeight" ,&totalWeight);
		tPhoton->SetBranchAddress("jet_n" ,&jet_n);
		tPhoton->SetBranchAddress("bjet_n" ,&bjet_n);
		tPhoton->SetBranchAddress("gamma_pt" ,&gamma_pt);
		tPhoton->SetBranchAddress("HT" ,&HT);
		tPhoton->SetBranchAddress("METl_raw" ,&METl);
		for (int entry=0;entry<tPhoton->GetEntries();entry++) {
			tPhoton->GetEntry(entry);
			if (gamma_pt<50.) continue;
			int pt = hist_low_pt->FindBin(gamma_pt)-1;
			if (jet_n!=1) continue;
			if (bjet_n!=0) continue;
			g_metl[pt]->Fill(METl,totalWeight);
		}
		fPhoton.Close();
	}
	else {
		std::cout << "Get smearing function from MC." << std::endl;
		cout << "Opening file           : " << TString(smearingPath)+"zjets/zjets_"+TString(ch)+".root"  << endl;
		TFile fZ( TString(smearingPath)+"zjets/zjets_"+TString(ch)+".root" );
		TTree*  tZ              = (TTree*)fZ.Get("BaselineTree");
		tZ->SetBranchStatus("*", 0);
		tZ->SetBranchStatus("totalWeight", 1);
		tZ->SetBranchStatus("jet_n", 1);
		tZ->SetBranchStatus("bjet_n", 1);
		tZ->SetBranchStatus("Z_pt", 1);
		tZ->SetBranchStatus("Z_truthPt", 1);
		tZ->SetBranchStatus("HT", 1);
		tZ->SetBranchStatus("mll", 1);
		tZ->SetBranchStatus("METl", 1);
		tZ->SetBranchStatus("RunNumber", 1);
		tZ->SetBranchStatus("EventNumber", 1);
		tZ->SetBranchAddress("totalWeight" ,&totalWeight);
		tZ->SetBranchAddress("jet_n" ,&jet_n);
		tZ->SetBranchAddress("bjet_n" ,&bjet_n);
		tZ->SetBranchAddress("Z_pt" ,&Z_pt);
		tZ->SetBranchAddress("Z_truthPt" ,&Z_truthPt);
		tZ->SetBranchAddress("HT" ,&HT);
		tZ->SetBranchAddress("mll" ,&mll);
		tZ->SetBranchAddress("METl" ,&METl);
		tZ->SetBranchAddress("RunNumber" ,&RunNumber);
		tZ->SetBranchAddress("EventNumber" ,&EventNumber);
		for (int entry=0;entry<tZ->GetEntries();entry++) {
			tZ->GetEntry(entry);
			if (Z_pt<50.) continue;
			int pt = hist_low_pt->FindBin(Z_pt)-1;
			if (jet_n==0) continue;
			if (jet_n>=2) z_metl_2j[pt]->Fill(METl,totalWeight);
			if (jet_n!=1) continue;
			if (EventNumber==62237 && RunNumber==361405) std::cout << "kick out large weight Ev 62237." << std::endl;
			if (EventNumber==62237 && RunNumber==361405) continue;  // this is a large weight event that screws the smearing function
			if (EventNumber==1416 && RunNumber==361414) std::cout << "kick out large weight Ev 1416." << std::endl;
			if (EventNumber==1416 && RunNumber==361414) continue;  // this is a large weight event that screws the smearing function
			if (EventNumber==3411 && RunNumber==361415) std::cout << "kick out large weight Ev 3411." << std::endl;
			if (EventNumber==3411 && RunNumber==361415) continue;  // this is a large weight event that screws the smearing function
			if (EventNumber==1449 && RunNumber==361415) std::cout << "kick out large weight Ev 1449." << std::endl;
			if (EventNumber==1449 && RunNumber==361415) continue;  // this is a large weight event that screws the smearing function
			z_metl[pt]->Fill(METl,totalWeight);
			z_dpt[pt]->Fill(Z_truthPt-Z_pt,totalWeight);
			if (mll<90 || mll>92) continue;
			z_jetmetl[pt]->Fill(METl,totalWeight);
		}
		fZ.Close();
		TFile fPhoton( TString(smearingPath)+"gmc/gmc_raw.root" );
		TTree*  tPhoton              = (TTree*)fPhoton.Get("BaselineTree");
		tPhoton->SetBranchStatus("*", 0);
		tPhoton->SetBranchStatus("totalWeight", 1);
		tPhoton->SetBranchStatus("jet_n", 1);
		tPhoton->SetBranchStatus("bjet_n", 1);
		tPhoton->SetBranchStatus("gamma_pt", 1);
		tPhoton->SetBranchStatus("truthGamma_pt", 1);
		tPhoton->SetBranchStatus("HT", 1);
		tPhoton->SetBranchStatus("METl_raw", 1);
		tPhoton->SetBranchAddress("totalWeight" ,&totalWeight);
		tPhoton->SetBranchAddress("jet_n" ,&jet_n);
		tPhoton->SetBranchAddress("bjet_n" ,&bjet_n);
		tPhoton->SetBranchAddress("gamma_pt" ,&gamma_pt);
		tPhoton->SetBranchAddress("truthGamma_pt" ,&truthGamma_pt);
		tPhoton->SetBranchAddress("HT" ,&HT);
		tPhoton->SetBranchAddress("METl_raw" ,&METl);
		for (int entry=0;entry<tPhoton->GetEntries();entry++) {
			tPhoton->GetEntry(entry);
			if (gamma_pt<50.) continue;
			int pt = hist_low_pt->FindBin(gamma_pt)-1;
			if (jet_n==0) continue;
			if (jet_n!=1) continue;
			g_metl[pt]->Fill(METl,totalWeight);
			g_dpt[pt]->Fill(truthGamma_pt-gamma_pt,totalWeight);
		}
		fPhoton.Close();
	}

}

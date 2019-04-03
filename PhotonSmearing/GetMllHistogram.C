TH1D* hist_Mll_dPt[bin_size][dpt_bin_size];

void GetMllHistogram(string ch,string period) {

    for (int bin0=0; bin0<bin_size; bin0++) {
        for (int bin1=0; bin1<dpt_bin_size; bin1++) {
            hist_Mll_dPt[bin0][bin1] = new TH1D(TString("hist_Mll_dPt_")+TString::Itoa(bin0,10)+TString("_")+TString::Itoa(bin1,10),"",mll_bin_size,mll_bin);
        }
    }

    cout << "Path is " << ntuple_path << endl;

    TH1D* hist_low_dpt = new TH1D("hist_low_dpt","",dpt_bin_size,dpt_bin);
    TH1D* hist_sm_pt = new TH1D("hist_sm_pt","",bin_size,sm_pt_bin);

    string filename = ntuple_path + "/ZMC16a/Zjets_merged_processed.root";
    cout << "Opening mll histo file : " << filename << endl;
    TFile fZ(filename.c_str());
    TTree* tZ = (TTree*)fZ.Get("BaselineTree");

    tZ->SetBranchStatus("*", 0);
    double totalWeight; SetInputBranch(tZ, "totalWeight", &totalWeight);
    float METl; SetInputBranch(tZ, "METl", &METl);
    int jet_n; SetInputBranch(tZ, "jet_n", &jet_n);
    int bjet_n; SetInputBranch(tZ, "bjet_n", &bjet_n);
    float Z_pt; SetInputBranch(tZ, "Z_pt", &Z_pt);
    float mll; SetInputBranch(tZ, "mll", &mll);
    std::vector<float>* lep_pT = new std::vector<float>(10); SetInputBranch(tZ, "lep_pT", &lep_pT);
    int channel; SetInputBranch(tZ, "channel", &channel);

    for (int entry=0; entry<tZ->GetEntries(); entry++) {
        tZ->GetEntry(entry);

        if( TString(ch).EqualTo("ee") && channel != 1 ) continue; // ee
        if( TString(ch).EqualTo("mm") && channel != 0 ) continue; // ee
        if (jet_n<2) continue;
        if (lep_pT->at(0)<leading_lep_pt_cut) continue;
        if (lep_pT->at(1)<second_lep_pt_cut) continue;
        int pt = hist_sm_pt->FindBin(Z_pt)-1;
        int dpt = hist_low_dpt->FindBin(METl)-1;
        if (dpt>=0 && pt>=0) hist_Mll_dPt[pt][dpt]->Fill(mll,totalWeight);
    }

    fZ.Close();
}

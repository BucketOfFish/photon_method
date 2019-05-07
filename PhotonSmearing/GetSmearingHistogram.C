// Takes input automatically from GetPhotonSmearing.C
TH1D* z_metl[bin_size];
TH1D* z_metl_2j[bin_size];
TH1D* g_metl[bin_size];
TH1D* z_jetmetl[bin_size];

TH1D* hist_low_pt = new TH1D("hist_low_pt","",bin_size,sm_pt_bin);

// SMEARING METHODS:
// 0 : no smearing
// 4 : R21 MC smearing
// 5 : R21 data smearing

using namespace std;

void GetSmearingHistogram(string ch, float lumi, string period, int smearing_method) {

    cout << "GetSmearingHistogram : smearing_method " << smearing_method << endl;

    for (int bin=0;bin<bin_size;bin++) {
        z_metl[bin] = new TH1D(TString("z_metl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        z_metl_2j[bin] = new TH1D(TString("z_metl_2j_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        z_jetmetl[bin] = new TH1D(TString("z_jetmetl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        g_metl[bin] = new TH1D(TString("g_metl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
    }

    //------------------------------------
    // SMEARING METHOD 5: for R21 data
    //------------------------------------

    if (smearing_method == 5) {  // R21 data-driven smearing function
        std::cout << "Get smearing function from R21 data." << std::endl;

        //--- smearing with R21 samples
        string datafilename = ntuple_path + "zdata/" + period + "_merged_processed.root";

        cout << "Opening data smearing file   : " << datafilename << endl;
        TFile fZ( datafilename.c_str() );

        TTree*  tZ              = (TTree*)fZ.Get("BaselineTree");
        tZ->SetBranchStatus("*", 0);
        double totalWeight; SetInputBranch(tZ, "totalWeight" ,&totalWeight);
        int jet_n; SetInputBranch(tZ, "jet_n" ,&jet_n);
        int bjet_n; SetInputBranch(tZ, "bjet_n" ,&bjet_n);
        float gZ_pt; SetInputBranch(tZ, "Z_pt" ,&gZ_pt);
        //float HT; SetInputBranch(tZ, "HT" ,&HT);
        float mll; SetInputBranch(tZ, "mll" ,&mll);
        float METl; SetInputBranch(tZ, "METl" ,&METl);
        int gchannel; SetInputBranch(tZ, "channel" ,&gchannel);
        for (int entry=0;entry<tZ->GetEntries();entry++) {
            tZ->GetEntry(entry);
            if( TString(ch).EqualTo("ee") && gchannel != 1 ) continue; // ee
            if( TString(ch).EqualTo("mm") && gchannel != 0 ) continue; // ee
            if (gZ_pt<50.) continue;
            int pt = hist_low_pt->FindBin(gZ_pt)-1;
            if (jet_n!=1) continue;
            if (bjet_n!=0) continue;
            z_metl[pt]->Fill(METl,totalWeight);
            if (mll<90 || mll>92) continue;
            z_jetmetl[pt]->Fill(METl,totalWeight);
        }

        //tZ->Close()
        fZ.Close();

        //--- smearing R21 samples
        string mcperiod = "";
        if( TString(period).EqualTo("data15-16") ) mcperiod = "ZMC16a/";
        if( TString(period).EqualTo("data17")    ) mcperiod = "ZMC16cd/";

        string ttfilename = ntuple_path + mcperiod + "ttbar_merged_processed.root";
        cout << "Opening tt smearing file   : " << ttfilename << endl;
        TFile ftt( ttfilename.c_str() );

        TTree*  ttt              = (TTree*)ftt.Get("BaselineTree");
        ttt->SetBranchStatus("*", 0);
        SetInputBranch(ttt, "totalWeight" ,&totalWeight);
        SetInputBranch(ttt, "jet_n" ,&jet_n);
        SetInputBranch(ttt, "bjet_n" ,&bjet_n);
        SetInputBranch(ttt, "Z_pt" ,&gZ_pt);
        //SetInputBranch(ttt, "HT" ,&HT);
        SetInputBranch(ttt, "mll" ,&mll);
        SetInputBranch(ttt, "METl" ,&METl);
        SetInputBranch(ttt, "channel" ,&gchannel);
        for (int entry=0;entry<ttt->GetEntries();entry++) {
            ttt->GetEntry(entry);
            if( TString(ch).EqualTo("ee") && gchannel != 1 ) continue; // ee
            if( TString(ch).EqualTo("mm") && gchannel != 0 ) continue; // ee			
            if (gZ_pt<50.) continue;
            int pt = hist_low_pt->FindBin(gZ_pt)-1;
            if (jet_n!=1) continue;
            if (bjet_n!=0) continue;
            z_metl[pt]->Fill(METl,-1.*lumi*totalWeight);
            if (mll<90 || mll>92) continue;
            z_jetmetl[pt]->Fill(METl,-1.*lumi*totalWeight);
        }
        //ttt->Close();
        ftt.Close();

        // smearing with R21 samples
        string vvfilename = ntuple_path + mcperiod + "diboson_merged_processed.root";
        cout << "Opening VV smearing file   : " << vvfilename << endl;
        TFile fvv( vvfilename.c_str() );

        TTree*  tvv              = (TTree*)fvv.Get("BaselineTree");
        tvv->SetBranchStatus("*", 0);
        SetInputBranch(tvv, "totalWeight" ,&totalWeight);
        SetInputBranch(tvv, "jet_n" ,&jet_n);
        SetInputBranch(tvv, "bjet_n" ,&bjet_n);
        SetInputBranch(tvv, "Z_pt" ,&gZ_pt);
        //SetInputBranch(tvv, "HT" ,&HT);
        SetInputBranch(tvv, "mll" ,&mll);
        SetInputBranch(tvv, "METl" ,&METl);
        SetInputBranch(tvv, "channel" ,&gchannel);
        for (int entry=0;entry<tvv->GetEntries();entry++) {
            tvv->GetEntry(entry);
            if( TString(ch).EqualTo("ee") && gchannel != 1 ) continue; // ee
            if( TString(ch).EqualTo("mm") && gchannel != 0 ) continue; // ee
            if (gZ_pt<50.) continue;
            int pt = hist_low_pt->FindBin(gZ_pt)-1;
            if (jet_n!=1) continue;
            if (bjet_n!=0) continue;
            z_metl[pt]->Fill(METl,-1.*lumi*totalWeight);
            if (mll<90 || mll>92) continue;
            z_jetmetl[pt]->Fill(METl,-1.*lumi*totalWeight);
        }
        //tvv->Close();
        fvv.Close();

        string gperiod = "";
        if( TString(period).EqualTo("data15-16") ) gperiod = "data15-16";
        if( TString(period).EqualTo("data17")    ) gperiod = "data17";

        string gfilename = ntuple_path + "gdata/" + gperiod + "_merged_processed.root";

        cout << "Opening photon smearing file   : " << gfilename << endl;
        TFile fPhoton( gfilename.c_str() );

        TTree*  tPhoton              = (TTree*)fPhoton.Get("BaselineTree");

        cout << "Setting photon branches" << endl;
        tPhoton->SetBranchStatus("*", 0);
        SetInputBranch(tPhoton, "totalWeight" ,&totalWeight);
        SetInputBranch(tPhoton, "jet_n" ,&jet_n);
        SetInputBranch(tPhoton, "bjet_n" ,&bjet_n);
        float gamma_pt; SetInputBranch(tPhoton, "gamma_pt" ,&gamma_pt);
        //float gamma_ht; SetInputBranch(tPhoton, "gamma_ht" ,&gamma_ht);
        //SetInputBranch(tPhoton, "HT" ,&HT);
        SetInputBranch(tPhoton, "METl_raw" ,&METl);
        cout << "Done setting photon branches" << endl;
        for (int entry=0;entry<tPhoton->GetEntries();entry++) {
            tPhoton->GetEntry(entry);
            if (gamma_pt<50.) continue;
            int pt = hist_low_pt->FindBin(gamma_pt)-1;
            if (jet_n!=1) continue;
            if (bjet_n!=0) continue;
            g_metl[pt]->Fill(METl,totalWeight);
        }
        //tPhoton->Close();
        fPhoton.Close();
    }

    //------------------------------------
    // SMEARING METHOD 4: for R21 MC
    //------------------------------------

    else if (smearing_method == 4) { // R21 MC-driven smearing function

        Float_t Z_ptD;
        Float_t mllD;

        std::cout << "Get smearing function from R21 MC." << std::endl;

        string mcperiod = "";
        if( TString(period).EqualTo("data15-16") ) mcperiod = "ZMC16a/";
        if( TString(period).EqualTo("data17")    ) mcperiod = "ZMC16cd/";

        string Zfilename = ntuple_path + mcperiod + "Zjets_merged_processed.root";

        cout << "Opening Z+jets MC smearing file           : " << Zfilename << endl;

        TFile fZ( Zfilename.c_str() );

        TTree*  tZ = (TTree*)fZ.Get("BaselineTree");
        tZ->SetBranchStatus("*", 0);
        double totalWeight; SetInputBranch(tZ, "totalWeight" ,&totalWeight);
        int jet_n; SetInputBranch(tZ, "jet_n" ,&jet_n);
        int bjet_n; SetInputBranch(tZ, "bjet_n" ,&bjet_n);
        float gZ_pt; SetInputBranch(tZ, "Z_pt" ,&gZ_pt);
        float HT; SetInputBranch(tZ, "HT" ,&HT);
        float mll; SetInputBranch(tZ, "mll" ,&mll);
        float METl; SetInputBranch(tZ, "METl" ,&METl);
        int gchannel; SetInputBranch(tZ, "channel" ,&gchannel);
        for (int entry=0;entry<tZ->GetEntries();entry++) {
            tZ->GetEntry(entry);

            if( TString(ch).EqualTo("ee") && gchannel != 1 ) continue; // ee
            if( TString(ch).EqualTo("mm") && gchannel != 0 ) continue; // ee			

            if (Z_ptD<50.) continue;
            int pt = hist_low_pt->FindBin(Z_ptD)-1;
            if (jet_n==0) continue;
            if (jet_n>=2) z_metl_2j[pt]->Fill(METl,totalWeight);
            if (jet_n!=1) continue;

            z_metl[pt]->Fill(METl,totalWeight);
            //z_dpt[pt]->Fill(Z_truthPt-Z_pt,totalWeight);
            if (mllD<90 || mllD>92) continue;
            z_jetmetl[pt]->Fill(METl,totalWeight);
        }
        fZ.Close();

        string gmcfilename = ntuple_path + "gmc/SinglePhoton222_merged_processed.root";

        cout << "Opening photon MC smearing file " << gmcfilename << endl;

        TFile fPhoton( gmcfilename.c_str() );

        TTree*  tPhoton              = (TTree*)fPhoton.Get("BaselineTree");
        tPhoton->SetBranchStatus("*", 0);
        SetInputBranch(tPhoton, "totalWeight" ,&totalWeight);
        SetInputBranch(tPhoton, "jet_n" ,&jet_n);
        SetInputBranch(tPhoton, "bjet_n" ,&bjet_n);
        float gamma_pt; SetInputBranch(tPhoton, "gamma_pt" ,&gamma_pt);
        //float gamma_ht; SetInputBranch(tPhoton, "gamma_ht" ,&gamma_ht);
        SetInputBranch(tPhoton, "HT" ,&HT);
        SetInputBranch(tPhoton, "METl_raw" ,&METl);
        for (int entry=0;entry<tPhoton->GetEntries();entry++) {
            tPhoton->GetEntry(entry);
            if (gamma_pt<50.) continue;
            int pt = hist_low_pt->FindBin(gamma_pt)-1;
            if (jet_n==0) continue;
            if (jet_n!=1) continue;
            g_metl[pt]->Fill(METl,totalWeight);
            //g_dpt[pt]->Fill(truthGamma_pt-gamma_pt,totalWeight);
        }
        fPhoton.Close();
    }

}

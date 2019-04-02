// Takes input automatically from GetPhotonSmearing.C
TH1D* z_dphi[bin_size];
TH1D* z_metl[bin_size];
TH1D* z_metl_2j[bin_size];
TH1D* tt_metl[bin_size];
TH1D* vv_metl[bin_size];
TH1D* g_metl[bin_size];
TH1D* z_jetmetl[bin_size];
TH1D* tt_jetmetl[bin_size];
TH1D* vv_jetmetl[bin_size];

TH1D* hist_low_pt = new TH1D("hist_low_pt","",bin_size,sm_pt_bin);

double totalWeightF = 0.;
float METlF = 0.;
float HTF = 0.;
float gZ_pt;
int gchannel;

// SMEARING METHODS:
// 0 : no smearing
// 4 : R21 MC smearing
// 5 : R21 data smearing

void GetSmearingHistogram(string ch, float lumi, string photon_tag,string period, int smearing_method) {

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
        //string datafilename = smearingPath + "Zdata/" + period + "_merged_processed.root";
        string datafilename = smearingPath + "zdata/" + period + "_merged_processed.root";

        cout << "Opening data smearing file   : " << datafilename << endl;
        TFile fZ( datafilename.c_str() );

        TTree*  tZ              = (TTree*)fZ.Get("BaselineTree");
        tZ->SetBranchStatus("*", 0);
        tZ->SetBranchStatus("totalWeight", 1);
        tZ->SetBranchStatus("jet_n", 1);
        tZ->SetBranchStatus("bjet_n", 1);
        tZ->SetBranchStatus("Z_pt", 1);
        //tZ->SetBranchStatus("HT", 1);
        tZ->SetBranchStatus("mll", 1);
        tZ->SetBranchStatus("METl", 1);
        tZ->SetBranchAddress("totalWeight" ,&totalWeightF);
        tZ->SetBranchAddress("jet_n" ,&jet_n);
        tZ->SetBranchAddress("bjet_n" ,&bjet_n);
        tZ->SetBranchAddress("Z_pt" ,&gZ_pt);
        //tZ->SetBranchAddress("HT" ,&HT);
        tZ->SetBranchAddress("mll" ,&mll);
        tZ->SetBranchAddress("METl" ,&METl);
        tZ->SetBranchAddress("channel" ,&gchannel);
        for (int entry=0;entry<tZ->GetEntries();entry++) {
            tZ->GetEntry(entry);
            if( TString(ch).EqualTo("ee") && gchannel != 1 ) continue; // ee
            if( TString(ch).EqualTo("mm") && gchannel != 0 ) continue; // ee
            if (gZ_pt<50.) continue;
            int pt = hist_low_pt->FindBin(gZ_pt)-1;
            if (jet_n!=1) continue;
            if (bjet_n!=0) continue;
            z_metl[pt]->Fill(METlF,totalWeightF);
            if (mll<90 || mll>92) continue;
            z_jetmetl[pt]->Fill(METlF,totalWeightF);
        }

        //tZ->Close()
        fZ.Close();

        //--- smearing R21 samples
        string mcperiod = "";
        if( TString(period).EqualTo("data15-16") ) mcperiod = "ZMC16a/";
        if( TString(period).EqualTo("data17")    ) mcperiod = "ZMC16cd/";

        string ttfilename = smearingPath + mcperiod + "ttbar_merged_processed.root";
        cout << "Opening tt smearing file   : " << ttfilename << endl;
        TFile ftt( ttfilename.c_str() );

        TTree*  ttt              = (TTree*)ftt.Get("BaselineTree");
        ttt->SetBranchStatus("*", 0);
        ttt->SetBranchStatus("totalWeight", 1);
        ttt->SetBranchStatus("jet_n", 1);
        ttt->SetBranchStatus("bjet_n", 1);
        ttt->SetBranchStatus("Z_pt", 1);
        //ttt->SetBranchStatus("HT", 1);
        ttt->SetBranchStatus("mll", 1);
        ttt->SetBranchStatus("METl", 1);
        ttt->SetBranchAddress("totalWeight" ,&totalWeightF);
        ttt->SetBranchAddress("jet_n" ,&jet_n);
        ttt->SetBranchAddress("bjet_n" ,&bjet_n);
        ttt->SetBranchAddress("Z_pt" ,&gZ_pt);
        //ttt->SetBranchAddress("HT" ,&HT);
        ttt->SetBranchAddress("mll" ,&mll);
        ttt->SetBranchAddress("METl" ,&METlF);
        ttt->SetBranchAddress("channel" ,&gchannel);
        for (int entry=0;entry<ttt->GetEntries();entry++) {
            ttt->GetEntry(entry);
            if( TString(ch).EqualTo("ee") && gchannel != 1 ) continue; // ee
            if( TString(ch).EqualTo("mm") && gchannel != 0 ) continue; // ee			
            if (gZ_pt<50.) continue;
            int pt = hist_low_pt->FindBin(gZ_pt)-1;
            if (jet_n!=1) continue;
            if (bjet_n!=0) continue;
            z_metl[pt]->Fill(METlF,-1.*lumi*totalWeightF);
            if (mll<90 || mll>92) continue;
            z_jetmetl[pt]->Fill(METlF,-1.*lumi*totalWeightF);
        }
        //ttt->Close();
        ftt.Close();

        // smearing with R21 samples
        string vvfilename = smearingPath + mcperiod + "diboson_merged_processed.root";
        cout << "Opening VV smearing file   : " << vvfilename << endl;
        TFile fvv( vvfilename.c_str() );

        TTree*  tvv              = (TTree*)fvv.Get("BaselineTree");
        tvv->SetBranchStatus("*", 0);
        tvv->SetBranchStatus("totalWeight", 1);
        tvv->SetBranchStatus("jet_n", 1);
        tvv->SetBranchStatus("bjet_n", 1);
        tvv->SetBranchStatus("Z_pt", 1);
        //tvv->SetBranchStatus("HT", 1);
        tvv->SetBranchStatus("mll", 1);
        tvv->SetBranchStatus("METl", 1);
        tvv->SetBranchAddress("totalWeight" ,&totalWeightF);
        tvv->SetBranchAddress("jet_n" ,&jet_n);
        tvv->SetBranchAddress("bjet_n" ,&bjet_n);
        tvv->SetBranchAddress("Z_pt" ,&gZ_pt);
        //tvv->SetBranchAddress("HT" ,&HT);
        tvv->SetBranchAddress("mll" ,&mll);
        tvv->SetBranchAddress("METl" ,&METlF);
        tvv->SetBranchAddress("channel" ,&gchannel);
        for (int entry=0;entry<tvv->GetEntries();entry++) {
            tvv->GetEntry(entry);
            if( TString(ch).EqualTo("ee") && gchannel != 1 ) continue; // ee
            if( TString(ch).EqualTo("mm") && gchannel != 0 ) continue; // ee
            if (gZ_pt<50.) continue;
            int pt = hist_low_pt->FindBin(gZ_pt)-1;
            if (jet_n!=1) continue;
            if (bjet_n!=0) continue;
            z_metl[pt]->Fill(METlF,-1.*lumi*totalWeightF);
            if (mll<90 || mll>92) continue;
            z_jetmetl[pt]->Fill(METlF,-1.*lumi*totalWeightF);
        }
        //tvv->Close();
        fvv.Close();

        string gperiod = "";
        if( TString(period).EqualTo("data15-16") ) gperiod = "data15-16";
        if( TString(period).EqualTo("data17")    ) gperiod = "data17";

        string gfilename = smearingPath + "gdata/" + gperiod + "_merged_processed.root";

        cout << "Opening photon smearing file   : " << gfilename << endl;
        TFile fPhoton( gfilename.c_str() );

        TTree*  tPhoton              = (TTree*)fPhoton.Get("BaselineTree");

        cout << "Setting photon branches" << endl;
        tPhoton->SetBranchStatus("*", 0);
        tPhoton->SetBranchStatus("totalWeight", 1);
        tPhoton->SetBranchStatus("jet_n", 1);
        tPhoton->SetBranchStatus("bjet_n", 1);
        tPhoton->SetBranchStatus("gamma_pt", 1);
        tPhoton->SetBranchStatus("gamma_ht", 1);
        //tPhoton->SetBranchStatus("HT", 1);
        tPhoton->SetBranchStatus("METl_raw", 1);
        tPhoton->SetBranchAddress("totalWeight" ,&totalWeight);
        tPhoton->SetBranchAddress("jet_n" ,&jet_n);
        tPhoton->SetBranchAddress("bjet_n" ,&bjet_n);
        tPhoton->SetBranchAddress("gamma_pt" ,&gamma_pt);
        tPhoton->SetBranchAddress("gamma_ht" ,&gamma_ht);
        //tPhoton->SetBranchAddress("HT" ,&HT);
        tPhoton->SetBranchAddress("METl_raw" ,&MET);
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

        string Zfilename = smearingPath + mcperiod + "Zjets_merged_processed.root";

        cout << "Opening Z+jets MC smearing file           : " << Zfilename << endl;

        TFile fZ( Zfilename.c_str() );

        TTree*  tZ = (TTree*)fZ.Get("BaselineTree");
        tZ->SetBranchStatus("*", 0);
        tZ->SetBranchStatus("totalWeight", 1);
        tZ->SetBranchStatus("jet_n", 1);
        tZ->SetBranchStatus("bjet_n", 1);
        tZ->SetBranchStatus("Z_pt", 1);
        //tZ->SetBranchStatus("Z_truthPt", 1);
        tZ->SetBranchStatus("HT", 1);
        tZ->SetBranchStatus("mll", 1);
        tZ->SetBranchStatus("METl", 1);
        tZ->SetBranchStatus("RunNumber", 1);
        tZ->SetBranchStatus("EventNumber", 1);
        tZ->SetBranchAddress("totalWeight" ,&totalWeight);
        tZ->SetBranchAddress("jet_n" ,&jet_n);
        tZ->SetBranchAddress("bjet_n" ,&bjet_n);
        tZ->SetBranchAddress("Z_pt" ,&Z_ptD);
        //tZ->SetBranchAddress("Z_truthPt" ,&Z_truthPt);
        tZ->SetBranchAddress("HT" ,&HT);
        tZ->SetBranchAddress("mll" ,&mllD);
        tZ->SetBranchAddress("METl" ,&METl);
        tZ->SetBranchAddress("RunNumber" ,&RunNumber);
        tZ->SetBranchAddress("EventNumber" ,&EventNumber);
        tZ->SetBranchAddress("channel" ,&gchannel);
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

        string gmcfilename = smearingPath + "gmc/SinglePhoton222_merged_processed.root";

        cout << "Opening photon MC smearing file " << gmcfilename << endl;

        TFile fPhoton( gmcfilename.c_str() );

        TTree*  tPhoton              = (TTree*)fPhoton.Get("BaselineTree");
        tPhoton->SetBranchStatus("*", 0);
        tPhoton->SetBranchStatus("totalWeight", 1);
        tPhoton->SetBranchStatus("jet_n", 1);
        tPhoton->SetBranchStatus("bjet_n", 1);
        tPhoton->SetBranchStatus("gamma_pt", 1);
        tPhoton->SetBranchStatus("gamma_ht", 1);
        //tPhoton->SetBranchStatus("truthGamma_pt", 1);
        tPhoton->SetBranchStatus("HT", 1);
        tPhoton->SetBranchStatus("METl_raw", 1);
        tPhoton->SetBranchAddress("totalWeight" ,&totalWeight);
        tPhoton->SetBranchAddress("jet_n" ,&jet_n);
        tPhoton->SetBranchAddress("bjet_n" ,&bjet_n);
        tPhoton->SetBranchAddress("gamma_pt" ,&gamma_pt);
        tPhoton->SetBranchAddress("gamma_ht" ,&gamma_ht);
        //tPhoton->SetBranchAddress("truthGamma_pt" ,&truthGamma_pt);
        tPhoton->SetBranchAddress("HT" ,&HT);
        tPhoton->SetBranchAddress("METl_raw" ,&METl);
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

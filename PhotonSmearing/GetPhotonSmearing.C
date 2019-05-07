#include "../Common/Settings.C"
#include "../Common/CommonLibraries.C"
#include "../Common/CommonFunctions.C"
#include "MT2.h"

using namespace std;

int RebinHistogram(TH1D* hist, int rebin) {

    // Does what?

    float negative_yield = 0.;
    float positive_yield = 0.;
    float core_yield = 0.;
    int remainder = 0;

    if (rebin==0) {
        rebin = 1;
        for (int bin=1;bin<=hist->GetNbinsX();bin++) {
            if (hist->GetBinContent(bin)>=0) positive_yield += hist->GetBinContent(bin);
            else negative_yield += hist->GetBinContent(bin);
        }
        remainder = hist->GetNbinsX() % 2;
        while ((abs(negative_yield/positive_yield)>0.005 || core_yield/positive_yield<0.4) && remainder==0 && rebin<=32) {
            hist->Rebin(2);
            rebin = rebin*2;
            remainder = hist->GetNbinsX() % 2;
            negative_yield = 0.;
            positive_yield = 0.;
            core_yield = 0.;
            for (int bin=1;bin<=hist->GetNbinsX();bin++) {
                if (hist->GetBinContent(bin)>=0) positive_yield += hist->GetBinContent(bin);
                else negative_yield += hist->GetBinContent(bin);
                if (abs(hist->GetBinCenter(bin)-hist->GetMean())<hist->GetRMS()) {
                    core_yield += hist->GetBinContent(bin); // core_yield = 68% for a perfect Guassian
                }
            }
        }
    }
    else {
        hist->Rebin(rebin);
    }

    for (int bin=1;bin<=hist->GetNbinsX();bin++) {
        hist->SetBinContent(bin,max(hist->GetBinContent(bin),0.));
    }

    return rebin;
}

TH1D* GetMllHistogram(string ch,string period) {

    TH1D* hist_Mll_dPt[bin_size][dpt_bin_size];

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

    for (int bin0=0; bin0<bin_size; bin0++) {
        for (int bin1=0; bin1<dpt_bin_size; bin1++) {
            int rebin = RebinHistogram(hist_Mll_dPt[bin0][bin1], 0);
        }
    }

    return hist_Mll_dPt;
}

TH1D* z_metl[bin_size];
TH1D* z_metl_2j[bin_size];
TH1D* g_metl[bin_size];
TH1D* z_jetmetl[bin_size];

void GetSmearingHistogram(string ch, float lumi, string period, int smearing_method) {

    // SMEARING METHODS:
    // 0 : no smearing
    // 4 : R21 MC smearing
    // 5 : R21 data smearing

    cout << "GetSmearingHistogram : smearing_method " << smearing_method << endl;

    for (int bin=0;bin<bin_size;bin++) {
        z_metl[bin] = new TH1D(TString("z_metl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        z_metl_2j[bin] = new TH1D(TString("z_metl_2j_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        z_jetmetl[bin] = new TH1D(TString("z_jetmetl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        g_metl[bin] = new TH1D(TString("g_metl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
    }

    TH1D* hist_low_pt = new TH1D("hist_low_pt","",bin_size,sm_pt_bin);

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

void GetIndividualLeptonInfo(TLorentzVector z_4vec) {

    TRandom myRandom;
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

void GetPhotonSmearing(string label, string ch, string isData, string period, int smearing_method) {

    cout << "channel         " << ch              << endl;
    cout << "period          " << period          << endl;
    cout << "isData?         " << isData          << endl;
    cout << "smearing path   " << smearing_path    << endl;
    cout << "smearing method " << smearing_method << endl;

    //-----------------------------
    // get and rebin mll histograms
    //-----------------------------

    std::cout << "Prepare Mll histograms..." << std::endl;
    TH1D* hist_Mll_dPt = GetMllHistogram(ch,period);

    //-----------------------------
    // prepare smearing functions
    //-----------------------------

    std::cout << "Prepare smearing histograms..." << std::endl;
    cout << "smearing_method    " << smearing_method << endl;

    TH1D* g_resp[bin_size];
    TH1D* z_resp[bin_size];
    TH1D* smear_raw[bin_size];
    TH1D* smear_fft_re[bin_size];
    TH1D* smear_fft_im[bin_size];
    TH1D* smear_fft_amp[bin_size];
    TH1D* smear_final[bin_size];
    float shift[bin_size];

    TH1D* g_metl_smear[bin_size];
    TH1D* g_metl_smear_2j[bin_size];
    for (int bin=0;bin<bin_size;bin++) {
        g_metl_smear[bin] = new TH1D(TString("g_metl_smear_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        g_metl_smear[bin]->SetStats(0);
        g_metl_smear_2j[bin] = new TH1D(TString("g_metl_smear_2j_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        g_metl_smear_2j[bin]->SetStats(0);
    }

    float lumi = GetLumi(period);
    GetSmearingHistogram(ch, lumi, period, smearing_method);

    TSpectrum pfinder;
    for (int bin=0;bin<bin_size;bin++) {
        int rebin = RebinHistogram(z_metl[bin],0);
        rebin = RebinHistogram(z_metl_2j[bin],rebin);
        rebin = RebinHistogram(g_metl[bin],rebin);
        rebin = RebinHistogram(g_metl_smear[bin],rebin);
        rebin = RebinHistogram(g_metl_smear_2j[bin],rebin);
        float gmetl_mean = g_metl[bin]->GetMean();
        float gmetl_rms = g_metl[bin]->GetRMS();
        float zmetl_mean = z_metl[bin]->GetMean();
        float zmetl_rms = z_metl[bin]->GetRMS();
        int newbin = 40000/rebin;
        Float_t *smear = new Float_t[2*((newbin+1)/2+1)];
        Float_t *fft_re = new Float_t[newbin];
        Float_t *fft_im = new Float_t[newbin];
        Double_t *z_smear_in = new Double_t[newbin];
        Double_t g_smear_in[newbin];
        Double_t j_resp_in[newbin];
        Double_t *z_resp_in = new Double_t[newbin];
        Double_t *g_resp_in = new Double_t[newbin];
        g_resp[bin] = new TH1D(TString("g_resp_")+TString::Itoa(bin,10),"",newbin,-30000,10000);
        z_resp[bin] = new TH1D(TString("z_resp_")+TString::Itoa(bin,10),"",newbin,-30000,10000);
        smear_raw[bin] = new TH1D(TString("smear_raw_")+TString::Itoa(bin,10),"",newbin,-30000,10000);
        smear_fft_re[bin] = new TH1D(TString("smear_fft_re_")+TString::Itoa(bin,10),"",newbin,-30000,10000);
        smear_fft_im[bin] = new TH1D(TString("smear_fft_im_")+TString::Itoa(bin,10),"",newbin,-30000,10000);
        smear_fft_amp[bin] = new TH1D(TString("smear_fft_amp_")+TString::Itoa(bin,10),"",newbin,-30000,10000);
        for (int i=0;i<newbin;i++) {
            z_smear_in[i] = max(z_metl[bin]->GetBinContent(i+1),0.);
            if (i<newbin/2) g_smear_in[i] = max(g_metl[bin]->GetBinContent(i+1+newbin/2),0.);
            else g_smear_in[i] = 0.;
            z_resp_in[i] = max(z_metl[bin]->GetBinContent(i+1),0.);
            g_resp_in[i] = max(g_metl[bin]->GetBinContent(i+1),0.);
            if (i<newbin/2) j_resp_in[i] = max(z_jetmetl[bin]->GetBinContent(i+1+newbin/2),0.);
            else j_resp_in[i] = 0.;
        }
        pfinder.Deconvolution(z_smear_in,g_smear_in,newbin,1000,1,1.0);
        pfinder.Deconvolution(z_resp_in,j_resp_in,newbin,1000,1,1.0);
        pfinder.Deconvolution(g_resp_in,j_resp_in,newbin,1000,1,1.0);
        for (int i=0;i<newbin;i++) {
            smear_raw[bin]->SetBinContent(i+1,z_smear_in[i]);
            g_resp[bin]->SetBinContent(i+1,g_resp_in[i]);
            z_resp[bin]->SetBinContent(i+1,z_resp_in[i]);
            smear[i] = z_smear_in[i];
        }
        float smear_mean = smear_raw[bin]->GetMean();
        float smear_rms = smear_raw[bin]->GetRMS();

        for (int i=0;i<newbin;i++) {
            if (gmetl_rms/zmetl_rms > 1.0) {
                smear_raw[bin]->SetBinContent(i+1,0.);
                smear[i] = 0.;
            }
            float smear_cut = 6.;
            if (ch=="mm" && sm_pt_bin[bin]>=0) smear_cut = 7.;
            if (abs(smear_raw[bin]->GetBinCenter(i+1)-smear_mean)/smear_rms>smear_cut) {
                smear_raw[bin]->SetBinContent(i+1,0.);
                smear[i] = 0.;
            }
        }

        shift[bin] = -g_metl[bin]->GetMean();
    }

    for (int bin=0;bin<bin_size;bin++) {
        smear_final[bin] = new TH1D(TString("smear_final_")+TString::Itoa(bin,10),"",500,-1000,1000);
        for (int i=0;i<500;i++) {
            int which_bin = smear_raw[bin]->FindBin(smear_final[bin]->GetBinCenter(i+1));
            smear_final[bin]->SetBinContent(i+1,smear_raw[bin]->GetBinContent(which_bin));
        }
    }

    //---------------------------------------------
    // get unsmeared input file
    //---------------------------------------------

    TH1::SetDefaultSumw2();

    string  infilename;
    if (isData == "MC") infilename = ntuple_path + "gmc/" + label + ".root";
    else if (isData == "Data") infilename = ntuple_path + "gdata/" + label + ".root";
    cout << __FILE__ << " " << __LINE__ << endl;

    TChain* inputTree = new TChain("BaselineTree");
    inputTree->Add( infilename.c_str() );

    cout << endl;
    cout << "Opening file           : " << infilename        << endl;
    cout << "Events in ntuple       : " << inputTree->GetEntries() << endl;

    //---------------------------------------------
    // create smeared output file
    //---------------------------------------------

    TH1::SetDefaultSumw2();

    string photon_tag = "";
    if (smearing_method == 0) photon_tag = "_NoSmear";
    if (smearing_method == 4) photon_tag = "_McSmear";
    if (smearing_method == 5) photon_tag = "_DataSmear";

    string outfilename;
    if (isData == "Data") outfilename = TString(TString(smearing_path)+"gdata/" + label + "_"+TString(ch)+TString(photon_tag)+".root"); 
    if (isData == "MC") outfilename = TString(TString(smearing_path)+"gmc/gmc_"+TString(ch)+TString(photon_tag)+".root"); 

    TFile* f = new TFile(outfilename.c_str(), "recreate");          
    TTree* BaselineTree = new TTree("BaselineTree", "baseline tree");

    cout << endl;
    cout << "Create file           : " << outfilename << endl;

    //-----------------------------
    // access existing branches
    //-----------------------------

    float gamma_phi; inputTree->SetBranchAddress("gamma_phi", &gamma_phi);

    double totalWeight; CopyBranch(inputTree, BaselineTree, "totalWeight", "totalWeight", &totalWeight, "D");
    int jet_n; CopyBranch(inputTree, BaselineTree, "jet_n", "jet_n", &jet_n, "I");
    float gamma_pt; CopyBranch(inputTree, BaselineTree, "gamma_pt", "gamma_pt",  &gamma_pt, "F");
    float gamma_eta; CopyBranch(inputTree, BaselineTree, "gamma_eta", "Z_eta",  &gamma_eta, "F");
    float METl; CopyBranch(inputTree, BaselineTree, "METl_raw", "METl", &METl, "F");
    float METt; CopyBranch(inputTree, BaselineTree, "METt_raw", "METt", &METt, "F");

    std::vector<int>* lepFlavor = new std::vector<int>(10); CopyBranch(inputTree, BaselineTree, "lepFlavor", "lepFlavor", &lepFlavor, "std::vector<int>");
    std::vector<int>* lepCharge = new std::vector<int>(10); CopyBranch(inputTree, BaselineTree, "lepCharge", "lepCharge", &lepCharge, "std::vector<int>");
    std::vector<float>* jet_pT = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jet_pT", "jet_pT", &jet_pT, "std::vector<float>");
    std::vector<float>* jet_eta = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jet_eta", "jet_eta", &jet_eta, "std::vector<float>");
    std::vector<float>* jet_phi = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jet_phi", "jet_phi", &jet_phi, "std::vector<float>");
    std::vector<float>* jet_m = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jet_m", "jet_m", &jet_m, "std::vector<float>");

    //-----------------------------
    // add new branches
    //-----------------------------

    float gamma_pt_smear; BaselineTree->Branch("Z_pt", &gamma_pt_smear, "Z_pt/F");
    float gamma_phi_smear; BaselineTree->Branch("Z_phi", &gamma_phi_smear, "Z_phi/F");
    int lep_n; BaselineTree->Branch("lep_n", &lep_n, "lep_n/I");
    float mll; BaselineTree->Branch("mll", &mll, "mll/F");
    float MET_smear; BaselineTree->Branch("MET", &MET_smear, "MET/F");
    float DPhi_METJetLeading_smear; BaselineTree->Branch("DPhi_METJetLeading", &DPhi_METJetLeading_smear, "DPhi_METJetLeading/F");
    float DPhi_METJetSecond_smear; BaselineTree->Branch("DPhi_METJetSecond", &DPhi_METJetSecond_smear, "DPhi_METJetSecond/F");
    float DPhi_METLepLeading_smear; BaselineTree->Branch("DPhi_METLepLeading", &DPhi_METLepLeading_smear, "DPhi_METLepLeading/F");
    float DPhi_METLepSecond_smear; BaselineTree->Branch("DPhi_METLepSecond", &DPhi_METLepSecond_smear, "DPhi_METLepSecond/F");
    float DPhi_METPhoton_smear; BaselineTree->Branch("DPhi_METPhoton", &DPhi_METPhoton_smear, "DPhi_METPhoton/F");
    float MT2W; BaselineTree->Branch("MT2W", &MT2W, "MT2W/F");

    //-----------------------------
    // loop over events
    //-----------------------------

    Long64_t nentries = inputTree->GetEntries();

    for (Long64_t i=0;i<nentries;i++) {

        if (fmod(i,1e5)==0) std::cout << i << " events processed." << std::endl;
        inputTree->GetEntry(i);

        if( jet_n == 0 ) continue;

        // use the smearing function to smear MET and pT in photon events
        float photon_smear = 0;
        float photon_smear_phi = 0;
        int smpt = hist_sm_pt->FindBin(gamma_pt)-1;
        if (smpt>=0) {
            if (smearing_method != 0) {
                if (smear_final[smpt]->Integral()>0) photon_smear = smear_final[smpt]->GetRandom() + shift[smpt];
                else photon_smear = shift[smpt];
            }
            else {
                photon_smear = 0;
            }
        }

        gamma_pt_smear = gamma_pt-photon_smear; // sign of photon_smear is important!!!
        gamma_phi_smear = gamma_phi-photon_smear_phi; // sign of photon_smear is important!!!

        float photon_smear_l = gamma_pt-gamma_pt_smear*TMath::Cos(photon_smear_phi);
        float photon_smear_t = -gamma_pt_smear*TMath::Sin(photon_smear_phi);
        float METl_smear = METl + photon_smear_l;  // sign of photon_smear is important!!!
        float METt_smear = METt + photon_smear_t;  // sign of photon_smear is important!!!
        float MET_smear = pow(METl_smear*METl_smear+METt_smear*METt_smear,0.5);

        int pt_smear = hist_sm_pt->FindBin(gamma_pt_smear)-1;
        if (gamma_pt_smear>sm_pt_bin[bin_size]) pt_smear = bin_size-1;
        int met_smear = hist_low_met->FindBin(MET_smear)-1;
        if (met_smear>met_bin[bin_size]) met_smear = bin_size-1;

        // recompute DPhi after smearing
        float METtx = METt*TMath::Cos(gamma_phi_smear+TMath::Pi()/2.);
        float METty = METt*TMath::Sin(gamma_phi_smear+TMath::Pi()/2.);
        float METlx_smear = METl_smear*TMath::Cos(gamma_phi_smear);
        float METly_smear = METl_smear*TMath::Sin(gamma_phi_smear);

        TLorentzVector met_4vec_smear;
        met_4vec_smear.SetXYZM(METtx+METlx_smear,METty+METly_smear,0,0);

        TLorentzVector jet0_4vec;
        if (jet_n<1) jet0_4vec.SetPtEtaPhiM(0,0,0,0);
        else jet0_4vec.SetPtEtaPhiM(jet_pT->at(0),jet_eta->at(0),jet_phi->at(0),jet_m->at(0));
        DPhi_METJetLeading_smear = fabs(met_4vec_smear.DeltaPhi(jet0_4vec));

        TLorentzVector jet1_4vec;
        if (jet_n<2) jet1_4vec.SetPtEtaPhiM(0,0,0,0);
        else jet1_4vec.SetPtEtaPhiM(jet_pT->at(1),jet_eta->at(1),jet_phi->at(1),jet_m->at(1));
        DPhi_METJetSecond_smear = fabs(met_4vec_smear.DeltaPhi(jet1_4vec));

        DPhi_METPhoton_smear = fabs(TMath::ATan2(METt,METl_smear));
        //int dphi_smear = hist_low_dphi->FindBin(DPhi_METPhoton_smear)-1;
        //if (dphi_smear>dphi_bin[bin_size]) dphi_smear = bin_size-1;

        if (gamma_pt>50. && jet_n==1) g_metl_smear[smpt]->Fill(METl_smear,totalWeight);
        if (gamma_pt>50. && jet_n>=2) g_metl_smear_2j[smpt]->Fill(METl_smear,totalWeight);

        // translate photon pT to dilepton sum pT, and compute HTincl for photon events
        float photon_2LPt = 0;
        //HTincl = HT + photon_2LPt;
        int dpt = hist_low_dpt->FindBin(METl_smear)-1;
        if (dpt>=0 && pt_smear>=0)
            if (hist_Mll_dPt[pt_smear][dpt]->Integral()>0)
                mll = hist_Mll_dPt[pt_smear][dpt]->GetRandom();

        //---------------------------------------------
        // compute two lepton kinematics
        //---------------------------------------------
        TLorentzVector z_4vec;
        z_4vec.SetPtEtaPhiM(gamma_pt,gamma_eta,gamma_phi,mll);
        GetDijetVariables(z_4vec, met_4vec_smear, jet_pT, jet_eta, jet_phi, jet_m);

        std::vector<float>* lep_phi = new std::vector<float>(10);
        std::vector<float>* lep_eta = new std::vector<float>(10);
        std::vector<float>* lep_pT = new std::vector<float>(10);
        lep_pT->push_back(0);
        lep_eta->push_back(0);
        lep_phi->push_back(0);
        lep_pT->push_back(0);
        lep_eta->push_back(0);
        lep_phi->push_back(0);
        int ntry = 0;
        while ((lep_pT->at(0)<leading_lep_pt_cut || lep_pT->at(1)<second_lep_pt_cut) && ntry<100) {
            ntry += 1;
            GetIndividualLeptonInfo(z_4vec);
        }

        TLorentzVector lep0_4vec;
        TLorentzVector lep1_4vec;
        lep_n = 2;
        lep0_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
        lep1_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);
        MT2W = ComputeMT2(lep0_4vec, lep1_4vec, met_4vec_smear, 0, 0).Compute();
        DPhi_METLepLeading_smear = fabs(met_4vec_smear.DeltaPhi(lep0_4vec));
        DPhi_METLepSecond_smear = fabs(met_4vec_smear.DeltaPhi(lep1_4vec));
        DR_2Lep = lep0_4vec.DeltaR(lep1_4vec);

        BaselineTree->Fill();
    }

    //-----------------------------
    // write tree and histograms
    //-----------------------------

    BaselineTree->Write();

    for (int bin=0;bin<bin_size;bin++) {
        g_metl[bin]->Write();
        z_metl[bin]->Write();
        z_metl_2j[bin]->Write();
        g_metl_smear[bin]->Write();
        g_metl_smear_2j[bin]->Write();
        smear_final[bin]->Write();
        if (smearing_method != 0) {
            g_resp[bin]->Write();
            z_resp[bin]->Write();
            smear_raw[bin]->Write();
            smear_fft_re[bin]->Write();
            smear_fft_im[bin]->Write();
            smear_fft_amp[bin]->Write();
        }
    }

    std::cout << "done." << std::endl;
    delete f;
}

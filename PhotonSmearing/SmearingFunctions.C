#include "../Common/Settings.C"
#include "../Common/CommonLibraries.C"
#include "../Common/CommonFunctions.C"
#include "MT2.h"

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

TH1D* z_metl[bins::smearing_bin_size];
TH1D* z_metl_2j[bins::smearing_bin_size];
TH1D* g_metl[bins::smearing_bin_size];
TH1D* z_jetmetl[bins::smearing_bin_size];

void GetSmearingHistogram(string ch, float lumi, string period, int smearing_method) {

    // SMEARING METHODS:
    // 0 : no smearing
    // 4 : R21 MC smearing
    // 5 : R21 data smearing

    cout << "GetSmearingHistogram : smearing_method " << smearing_method << endl;

    for (int bin=0;bin<bins::smearing_bin_size;bin++) {
        z_metl[bin] = new TH1D(TString("z_metl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        z_metl_2j[bin] = new TH1D(TString("z_metl_2j_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        z_jetmetl[bin] = new TH1D(TString("z_jetmetl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        g_metl[bin] = new TH1D(TString("g_metl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
    }

    TH1D* hist_low_pt = new TH1D("hist_low_pt","",bins::smearing_bin_size,bins::pt_bins);

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
        int jet_n; SetInputBranch(tZ, "nJet30" ,&jet_n);
        int bjet_n; SetInputBranch(tZ, "bjet_n" ,&bjet_n);
        float gZ_pt; SetInputBranch(tZ, "Ptll" ,&gZ_pt);
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
        SetInputBranch(ttt, "nJet30" ,&jet_n);
        SetInputBranch(ttt, "bjet_n" ,&bjet_n);
        SetInputBranch(ttt, "Ptll" ,&gZ_pt);
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
        SetInputBranch(tvv, "nJet30" ,&jet_n);
        SetInputBranch(tvv, "bjet_n" ,&bjet_n);
        SetInputBranch(tvv, "Ptll" ,&gZ_pt);
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
        int jet_n; SetInputBranch(tZ, "nJet30" ,&jet_n);
        int bjet_n; SetInputBranch(tZ, "bjet_n" ,&bjet_n);
        float gZ_pt; SetInputBranch(tZ, "Ptll" ,&gZ_pt);
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

std::vector<float>* lep_pT = new std::vector<float>(10); 
std::vector<float>* lep_eta = new std::vector<float>(10); 
std::vector<float>* lep_phi = new std::vector<float>(10); 
std::vector<float>* lep_flavor = new std::vector<int>(10); 
std::vector<float>* lep_charge = new std::vector<int>(10); 

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


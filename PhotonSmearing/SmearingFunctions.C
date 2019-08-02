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

std::vector<TLorentzVector>* SplitLeptons(TLorentzVector z_4vec) {

    double Z_pT = z_4vec.Pt();
    double Z_E = z_4vec.E();
    double Z_m = z_4vec.M();
    double e_m = 0.000511;
    double m_m = 0.106;
    double min_pT = min(cuts::leading_lep_pt_cut, cuts::second_lep_pt_cut);

    // Decay into electrons or muons?
    TRandomMixMax myRandom(0);
    int l_type = myRandom.Integer(1);
    double l_m;
    if (l_type == 0) l_m = e_m;
    else l_m = m_m;

    // Minimum decay angle in CM frame (in Z flight coordinate system)
    double l_pT_cm = sqrt(std::pow(Z_m, 2)/4 - std::pow(l_m, 2));
    double cos_min_theta = (Z_E - 2*sqrt(std::pow(l_m, 2) + std::pow(min_pT, 2))) / Z_pT;
    double min_theta = 0;
    if (cos_min_theta < 1) min_theta = acos(cos_min_theta);

    // Lepton decay angle sampling (in Z flight coordinate system)
    double l_phi_cm = myRandom.Rndm(0) * 2.*TMath::Pi();
    double l_theta_cm = myRandom.Rndm(0) * (TMath::Pi()/2-min_theta) + min_theta;
    double l_eta_cm = -log(tan(l_theta_cm/2));

    // Split γ/Z to two leptons in rest frame (in Z flight coordinate system)
    double l_px = l_pT_cm*TMath::Sin(l_theta_cm)*TMath::Cos(l_phi_cm);
    double l_py = l_pT_cm*TMath::Sin(l_theta_cm)*TMath::Sin(l_phi_cm);
    double l_pz = l_pT_cm*TMath::Cos(l_theta_cm);
    double l_E = Z_m/2;

    // Boost to lab frame
    double gamma = z_4vec.Gamma();
    double beta = z_4vec.Beta();
    double l_px_1 = gamma*(l_px - beta*l_E);
    double l_px_2 = gamma*(-l_px - beta*l_E);
    double l_E_1 = gamma*(l_E - beta*l_px);
    double l_E_2 = gamma*(l_E + beta*l_px);
	TLorentzVector l0_4vec, l1_4vec;
    l0_4vec.SetPxPyPzE(l_px_1, l_py, l_pz, l_E_1);
    l1_4vec.SetPxPyPzE(l_px_2, -l_py, -l_pz, l_E_2);

    // Rotate angle to coordinate system of detector

	// Checks
    TLorentzVector twol_4vec = l0_4vec + l1_4vec;
    std::cout << "z_4vec pT = " << z_4vec.Pt() << ", eta = " << z_4vec.Eta() << ", phi = " << z_4vec.Phi() << ", m = " << z_4vec.M() << std::endl;
    std::cout << "l0_4vec pT = " << l0_4vec.Pt() << ", eta = " << l0_4vec.Eta() << ", phi = " << l0_4vec.Phi() << ", m = " << l0_4vec.M() << std::endl;
    std::cout << "l1_4vec pT = " << l1_4vec.Pt() << ", eta = " << l1_4vec.Eta() << ", phi = " << l1_4vec.Phi() << ", m = " << l1_4vec.M() << std::endl;
    std::cout << "2l_4vec pT = " << twol_4vec.Pt() << ", eta = " << twol_4vec.Eta() << ", phi = " << twol_4vec.Phi() << ", m = " << twol_4vec.M() << std::endl;
    std::cout << "==================================================================================" << std::endl;

    // Return leptons
    std::vector<TLorentzVector>* leptons = new std::vector<TLorentzVector>;
    leptons->push_back(l0_4vec);
    leptons->push_back(l1_4vec);
    return leptons;
}

std::vector<TLorentzVector>* SplitLeptonsBackup(TLorentzVector z_4vec) {

    double Z_pT = z_4vec.Pt();
    double Z_E = z_4vec.E();
    double Z_m = z_4vec.M();
    double e_m = 0.000511;
    double m_m = 0.106;
    double min_pT = min(cuts::leading_lep_pt_cut, cuts::second_lep_pt_cut);

    // Decay into electrons or muons?
    TRandomMixMax myRandom(0);
    int l_type = myRandom.Integer(1);
    double l_m;
    if (l_type == 0) l_m = e_m;
    else l_m = m_m;

    // Minimum decay angle in CM frame
    double l_pT_cm = sqrt(std::pow(Z_m, 2)/4 - std::pow(l_m, 2));
    double cos_min_theta = (Z_E - 2*sqrt(std::pow(l_m, 2) + std::pow(min_pT, 2))) / Z_pT;
    double min_theta = 0;
    if (cos_min_theta < 1) min_theta = acos(cos_min_theta);
    min_theta = 1;
    //std::cout << Z_E << " " << Z_pT << " " << Z_m << " " << sqrt(pow(Z_m,2) + pow(Z_pT,2)) << " " << min_pT << " " << (Z_E - 2*sqrt(std::pow(l_m, 2) + std::pow(min_pT, 2))) / Z_pT << std::endl;

    // Lepton decay angle sampling (relative to path of Z flight)
    double l_phi_cm = myRandom.Rndm() * 2.*TMath::Pi();
    double l_theta_cm = myRandom.Rndm() * (TMath::Pi()/2-min_theta) + min_theta;
    double l_eta_cm = -log(tan(l_theta_cm/2));

    // Split γ/Z to two leptons in rest frame
	TLorentzVector l0_cm_4vec, l1_cm_4vec;
    //std::cout << z_4vec.Px() << " " << z_4vec.Py() << " " << z_4vec.Pz() << std::endl;
    std::cout << l_theta_cm << " " << l_eta_cm << " " << l_phi_cm << std::endl;
    double l_px = l_pT_cm*TMath::Sin(l_theta_cm)*TMath::Cos(l_phi_cm);
    double l_py = l_pT_cm*TMath::Sin(l_theta_cm)*TMath::Sin(l_phi_cm);
    double l_pz = l_pT_cm*TMath::Cos(l_theta_cm);
    l0_cm_4vec.SetPxPyPzE(l_px, l_py, l_pz, Z_m/2);
    l1_cm_4vec.SetPxPyPzE(-l_px, -l_py, -l_pz, Z_m/2);
    //l0_cm_4vec.SetPtEtaPhiM(l_pT_cm, l_eta_cm, l_phi_cm, l_m);
    //l1_cm_4vec.SetPtEtaPhiM(l_pT_cm, -l_eta_cm, l_phi_cm-TMath::Pi(), l_m);

    // Boost to lab frame along direction of Z flight
    //std::cout << boost_vec.X() << " " << boost_vec.Y() << " " << boost_vec.Z() << std::endl;
	TLorentzVector l0_lab_4vec = l0_cm_4vec;
    TLorentzVector l1_lab_4vec = l1_cm_4vec;
    TVector3 boost_vec = z_4vec.BoostVector();
	l0_lab_4vec.Boost(boost_vec);
	l1_lab_4vec.Boost(boost_vec);

    // Rotate angle to frame of detector
    std::cout << "boost = " << boost_vec.X() << ", " << boost_vec.Y() << ", " << boost_vec.Z() << " - Mag = " << boost_vec.Mag() << std::endl;
    TLorentzVector z_cm_4vec = z_4vec; 
    z_cm_4vec.Boost(-boost_vec);
    std::cout << "z_4vec pT = " << z_4vec.Pt() << ", eta = " << z_4vec.Eta() << ", phi = " << z_4vec.Phi() << ", m = " << z_4vec.M() << std::endl;
    std::cout << "z_cm_4vec pT = " << z_cm_4vec.Pt() << ", eta = " << z_cm_4vec.Eta() << ", phi = " << z_cm_4vec.Phi() << ", m = " << z_cm_4vec.M() << std::endl;
    std::cout << acos((cos(z_4vec.Phi())/z_4vec.Gamma())/sqrt(1 - pow(z_4vec.Beta()*cos(z_4vec.Phi()), 2))) << std::endl;
    std::cout << "----------------------------------------------------------------------------------" << std::endl;
    //double l_theta_cm = z_4vec.Theta() + l_dtheta*TMath::Cos(l_dphi);
    //double l_phi_cm = z_4vec.Phi() + l_dtheta*TMath::Sin(l_dphi);

	// Checks
    TLorentzVector twol_cm_4vec = l0_cm_4vec + l1_cm_4vec;
    TLorentzVector twol_lab_4vec = l0_lab_4vec + l1_lab_4vec;
    std::cout << "z_4vec pT = " << z_4vec.Pt() << ", eta = " << z_4vec.Eta() << ", phi = " << z_4vec.Phi() << ", m = " << z_4vec.M() << std::endl;
    //std::cout << "l_pT_cm = " << l_pT_cm << ", min_theta = " << min_theta << ", phi = " << l_phi_cm << ", theta = " << l_theta_cm << std::endl;
    std::cout << "l0_cm_4vec pT = " << l0_cm_4vec.Pt() << ", eta = " << l0_cm_4vec.Eta() << ", phi = " << l0_cm_4vec.Phi() << ", m = " << l0_cm_4vec.M() << std::endl;
    std::cout << "l1_cm_4vec pT = " << l1_cm_4vec.Pt() << ", eta = " << l1_cm_4vec.Eta() << ", phi = " << l1_cm_4vec.Phi() << ", m = " << l1_cm_4vec.M() << std::endl;
    std::cout << "2l_cm_4vec pT = " << twol_cm_4vec.Pt() << ", eta = " << twol_cm_4vec.Eta() << ", phi = " << twol_cm_4vec.Phi() << ", m = " << twol_cm_4vec.M() << std::endl;
    std::cout << "l0_lab_4vec pT = " << l0_lab_4vec.Pt() << ", eta = " << l0_lab_4vec.Eta() << ", phi = " << l0_lab_4vec.Phi() << ", m = " << l0_lab_4vec.M() << std::endl;
    std::cout << "l1_lab_4vec pT = " << l1_lab_4vec.Pt() << ", eta = " << l1_lab_4vec.Eta() << ", phi = " << l1_lab_4vec.Phi() << ", m = " << l1_lab_4vec.M() << std::endl;
    std::cout << "2l_lab_4vec pT = " << twol_lab_4vec.Pt() << ", eta = " << twol_lab_4vec.Eta() << ", phi = " << twol_lab_4vec.Phi() << ", m = " << twol_lab_4vec.M() << std::endl;
    std::cout << "==================================================================================" << std::endl;

    // Return leptons
    std::vector<TLorentzVector>* leptons = new std::vector<TLorentzVector>;
    leptons->push_back(l0_lab_4vec);
    leptons->push_back(l1_lab_4vec);
    return leptons;
}


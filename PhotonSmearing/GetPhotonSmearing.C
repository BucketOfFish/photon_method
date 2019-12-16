#include "../Common/Settings.C"
#include "SmearingFunctions.C"

using namespace std;

void GetPhotonSmearing(string period, string channel, string data_or_mc) {

    TString data_period = DataPeriod(period);
    TString mc_period = MCPeriod(period);

    cout << "channel         " << channel         << endl;
    cout << "period          " << period          << endl;
    cout << "isData?         " << data_or_mc          << endl;
    cout << "smearing path   " << smearing_path   << endl;

    //---------------------------------------------
    // get unsmeared input file
    //---------------------------------------------

    TH1::SetDefaultSumw2();

    TString infilename;
    if (data_or_mc == "MC") infilename = ntuple_path + "g_mc/" + mc_period + "_SinglePhoton222.root";
    else if (data_or_mc == "Data") infilename = ntuple_path + "g_data/" + data_period + "_photon.root";

    TChain* inputTree = new TChain("BaselineTree");
    inputTree->Add(infilename);

    cout << endl;
    cout << "Opening file           : " << infilename        << endl;
    cout << "Events in ntuple       : " << inputTree->GetEntries() << endl;

    //---------------------------------------------
    // create smeared output file
    //---------------------------------------------

    TH1::SetDefaultSumw2();

    TString outfilename;
    if (data_or_mc == "Data") outfilename = smearing_path+"g_data/"+data_period+"_photon_"+channel+".root"; 
    if (data_or_mc == "MC") outfilename = smearing_path+"g_mc/"+mc_period+"_SinglePhoton222_"+channel+".root";

    TFile* outputFile = new TFile(outfilename, "recreate");          
    TTree* BaselineTree = new TTree("BaselineTree", "baseline tree");
    BaselineTree->SetDirectory(outputFile);

    cout << endl;
    cout << "Create file           : " << outfilename << endl;

    //-----------------------------
    // access, copy, and create branches
    //-----------------------------

    float gamma_phi; inputTree->SetBranchAddress("gamma_phi", &gamma_phi);

    double totalWeight; CopyBranch(inputTree, BaselineTree, "totalWeight", "totalWeight", &totalWeight, "D");
    int bjet_n; CopyBranch(inputTree, BaselineTree, "bjet_n", "bjet_n", &bjet_n, "I");
    float gamma_pt; CopyBranch(inputTree, BaselineTree, "gamma_pt", "gamma_pt",  &gamma_pt, "F");
    float gamma_eta; CopyBranch(inputTree, BaselineTree, "gamma_eta", "Z_eta",  &gamma_eta, "F");
    float METl; CopyBranch(inputTree, BaselineTree, "METl_raw", "METl_raw", &METl, "F");
    float METt; CopyBranch(inputTree, BaselineTree, "METt_raw", "METt_raw", &METt, "F");
    float METl_smeared; BaselineTree->Branch("METl", &METl_smeared, "METl/F");
    float METt_smeared; BaselineTree->Branch("METt", &METt_smeared, "METt/F");
    float HT; CopyBranch(inputTree, BaselineTree, "HT", "HT", &HT, "F");
    float MET_raw; CopyBranch(inputTree, BaselineTree, "met_Et_raw", "met_Et_raw", &MET_raw, "F");
    float MET_phi; CopyBranch(inputTree, BaselineTree, "met_Phi", "met_Phi", &MET_phi, "F");

    vector<float>* jet_pT = new vector<float>(10); CopyBranch(inputTree, BaselineTree, "jet_pT", "jet_pT", &jet_pT, "vector<float>");
    vector<float>* jet_eta = new vector<float>(10); CopyBranch(inputTree, BaselineTree, "jet_eta", "jet_eta", &jet_eta, "vector<float>");
    vector<float>* jet_phi = new vector<float>(10); CopyBranch(inputTree, BaselineTree, "jet_phi", "jet_phi", &jet_phi, "vector<float>");

    float gamma_pt_smeared; BaselineTree->Branch("Ptll", &gamma_pt_smeared, "Ptll/F");
    float gamma_phi_smeared; BaselineTree->Branch("Z_phi", &gamma_phi_smeared, "Z_phi/F");
    int lep_n; BaselineTree->Branch("lep_n", &lep_n, "lep_n/I"); BaselineTree->Branch("nLep_signal", &lep_n, "nLep_signal/I"); BaselineTree->Branch("nLep_base", &lep_n, "nLep_base/I");
    float MET_smeared; BaselineTree->Branch("met_Et", &MET_smeared, "met_Et/F");
    float DPhi_METLepLeading_smeared; BaselineTree->Branch("DPhi_METLepLeading", &DPhi_METLepLeading_smeared, "DPhi_METLepLeading/F");
    float DPhi_METLepSecond_smeared; BaselineTree->Branch("DPhi_METLepSecond", &DPhi_METLepSecond_smeared, "DPhi_METLepSecond/F");
    float DPhi_METPhoton_smear; BaselineTree->Branch("DPhi_METPhoton", &DPhi_METPhoton_smear, "DPhi_METPhoton/F");
    float MT2; BaselineTree->Branch("mt2leplsp_0", &MT2, "mt2leplsp_0/F");
    float DR_2Lep; BaselineTree->Branch("DR_2Lep", &DR_2Lep, "DR_2Lep/F");
    int photon_conversion_type; CopyBranch(inputTree, BaselineTree, "PhotonConversionType", "PhotonConversionType", &photon_conversion_type, "I");
    float lep_theta_cm; BaselineTree->Branch("lep_theta_cm", &lep_theta_cm, "lep_theta_cm/F");

    //--- HistFitter branches
    int DatasetNumber; CopyBranch(inputTree, BaselineTree, "DatasetNumber", "DatasetNumber", &DatasetNumber, "I");
    float Etall; CopyBranch(inputTree, BaselineTree, "Etall", "Etall", &Etall, "F");
    double H2PP; CopyBranch(inputTree, BaselineTree, "H2PP", "H2PP", &H2PP, "D");
    double H5PP; CopyBranch(inputTree, BaselineTree, "H5PP", "H5PP", &H5PP, "D");
    double H5PP_VR; CopyBranch(inputTree, BaselineTree, "H5PP_VR", "H5PP_VR", &H5PP_VR, "D");
    float METOverPtISR; CopyBranch(inputTree, BaselineTree, "METOverPtISR", "METOverPtISR", &METOverPtISR, "F");
    float METOverPtW; CopyBranch(inputTree, BaselineTree, "METOverPtW", "METOverPtW", &METOverPtW, "F");
    float METOverPtZ; CopyBranch(inputTree, BaselineTree, "METOverPtZ", "METOverPtZ", &METOverPtZ, "F");
    double MJ; CopyBranch(inputTree, BaselineTree, "MJ", "MJ", &MJ, "D");
    double MJ_VR; CopyBranch(inputTree, BaselineTree, "MJ_VR", "MJ_VR", &MJ_VR, "D");
    double MZ; CopyBranch(inputTree, BaselineTree, "MZ", "MZ", &MZ, "D");
    double MZ_VR; CopyBranch(inputTree, BaselineTree, "MZ_VR", "MZ_VR", &MZ_VR, "D");
    double NjISR; CopyBranch(inputTree, BaselineTree, "NjISR", "NjISR", &NjISR, "D");
    double NjS; CopyBranch(inputTree, BaselineTree, "NjS", "NjS", &NjS, "D");
    double PTCM; CopyBranch(inputTree, BaselineTree, "PTCM", "PTCM", &PTCM, "D");
    double PTCM_VR; CopyBranch(inputTree, BaselineTree, "PTCM_VR", "PTCM_VR", &PTCM_VR, "D");
    double PTI; CopyBranch(inputTree, BaselineTree, "PTI", "PTI", &PTI, "D");
    double PTISR; CopyBranch(inputTree, BaselineTree, "PTISR", "PTISR", &PTISR, "D");
    double PTISR_VR; CopyBranch(inputTree, BaselineTree, "PTISR_VR", "PTISR_VR", &PTISR_VR, "D");
    double PTI_VR; CopyBranch(inputTree, BaselineTree, "PTI_VR", "PTI_VR", &PTI_VR, "D");
    //float Ptll; CopyBranch(inputTree, BaselineTree, "Ptll", "Ptll", &Ptll, "F");
    double RISR; CopyBranch(inputTree, BaselineTree, "RISR", "RISR", &RISR, "D");
    double RISR_VR; CopyBranch(inputTree, BaselineTree, "RISR_VR", "RISR_VR", &RISR_VR, "D");
    double RPT_HT5PP; CopyBranch(inputTree, BaselineTree, "RPT_HT5PP", "RPT_HT5PP", &RPT_HT5PP, "D");
    double RPT_HT5PP_VR; CopyBranch(inputTree, BaselineTree, "RPT_HT5PP_VR", "RPT_HT5PP_VR", &RPT_HT5PP_VR, "D");
    double R_minH2P_minH3P; CopyBranch(inputTree, BaselineTree, "R_minH2P_minH3P", "R_minH2P_minH3P", &R_minH2P_minH3P, "D");
    double R_minH2P_minH3P_VR; CopyBranch(inputTree, BaselineTree, "R_minH2P_minH3P_VR", "R_minH2P_minH3P_VR", &R_minH2P_minH3P_VR, "D");
    float Rjj; CopyBranch(inputTree, BaselineTree, "Rjj", "Rjj", &Rjj, "F");
    float Rll; CopyBranch(inputTree, BaselineTree, "Rll", "Rll", &Rll, "F");
    float dPhiMetISR; CopyBranch(inputTree, BaselineTree, "dPhiMetISR", "dPhiMetISR", &dPhiMetISR, "F");
    float dPhiMetJet1; CopyBranch(inputTree, BaselineTree, "dPhiMetJet1", "dPhiMetJet1", &dPhiMetJet1, "F");
    float dPhiMetJet2; CopyBranch(inputTree, BaselineTree, "dPhiMetJet2", "dPhiMetJet2", &dPhiMetJet2, "F");
    float dPhiMetJet12Min; CopyBranch(inputTree, BaselineTree, "dPhiMetJet12Min", "dPhiMetJet12Min", &dPhiMetJet12Min, "F");
    float dPhiPjjMet; CopyBranch(inputTree, BaselineTree, "dPhiPjjMet", "dPhiPjjMet", &dPhiPjjMet, "F");
    float dPhiPllMet; CopyBranch(inputTree, BaselineTree, "dPhiPllMet", "dPhiPllMet", &dPhiPllMet, "F");
    double dphiISRI; CopyBranch(inputTree, BaselineTree, "dphiISRI", "dphiISRI", &dphiISRI, "D");
    double dphiISRI_VR; CopyBranch(inputTree, BaselineTree, "dphiISRI_VR", "dphiISRI_VR", &dphiISRI_VR, "D");
    double dphiVP; CopyBranch(inputTree, BaselineTree, "dphiVP", "dphiVP", &dphiVP, "D");
    double dphiVP_VR; CopyBranch(inputTree, BaselineTree, "dphiVP_VR", "dphiVP_VR", &dphiVP_VR, "D");
    double lept1Pt_VR; CopyBranch(inputTree, BaselineTree, "lept1Pt_VR", "lept1Pt_VR", &lept1Pt_VR, "D");
    double lept2Pt_VR; CopyBranch(inputTree, BaselineTree, "lept2Pt_VR", "lept2Pt_VR", &lept2Pt_VR, "D");
    double mTl3; CopyBranch(inputTree, BaselineTree, "mTl3", "mTl3", &mTl3, "D");
    //float MET; CopyBranch(inputTree, BaselineTree, "met_Et", "met_Et", &MET, "F");
    float met_Sign; CopyBranch(inputTree, BaselineTree, "met_Sign", "met_Sign", &met_Sign, "F");
    double minDphi; CopyBranch(inputTree, BaselineTree, "minDphi", "minDphi", &minDphi, "D");
    double mll_RJ; CopyBranch(inputTree, BaselineTree, "mll_RJ", "mll_RJ", &mll_RJ, "D");
    double mll_RJ_VR; CopyBranch(inputTree, BaselineTree, "mll_RJ_VR", "mll_RJ_VR", &mll_RJ_VR, "D");
    float mt2leplsp_0; CopyBranch(inputTree, BaselineTree, "mt2leplsp_0", "mt2leplsp_0", &mt2leplsp_0, "F");
    int nJet20; CopyBranch(inputTree, BaselineTree, "nJet20", "nJet20", &nJet20, "I");
    float mjj; CopyBranch(inputTree, BaselineTree, "mjj", "mjj", &mjj, "F");
    float mll; CopyBranch(inputTree, BaselineTree, "mll", "mll", &mll, "F");
    int lepIsPR; CopyBranch(inputTree, BaselineTree, "lepIsPR", "lepIsPR", &lepIsPR, "I");

    vector<float>* jet_m = new vector<float>(10); CopyBranch(inputTree, BaselineTree, "jetM", "jetM", &jet_m, "vector<float>");
    Int_t jet_n; CopyBranch(inputTree, BaselineTree, "nJet30", "nJet30", &jet_n, "I");
    vector<float>* lep_pT = new vector<float>(10); BaselineTree->Branch("lepPt", "vector<float>", &lep_pT);
    vector<float>* lep_eta = new vector<float>(10); BaselineTree->Branch("lep_eta", "vector<float>", &lep_eta);
    vector<float>* lep_phi = new vector<float>(10); BaselineTree->Branch("lep_phi", "vector<float>", &lep_phi);
    vector<int>* lep_flavor = new vector<int>(10); BaselineTree->Branch("lepFlavor", "vector<int>", &lep_flavor);
    vector<int>* lep_charge = new vector<int>(10); BaselineTree->Branch("lepCharge", "vector<int>", &lep_charge);
    Int_t lepChannel; BaselineTree->Branch("channel", &lepChannel, "channel/I");
    int nBJet20_MV2c10_FixedCutBEff_77; CopyBranch(inputTree, BaselineTree, "nBJet20_MV2c10_FixedCutBEff_77", "nBJet20_MV2c10_FixedCutBEff_77", &nBJet20_MV2c10_FixedCutBEff_77, "I");
    bool trigMatch_2LTrigOR; CopyBranch(inputTree, BaselineTree, "trigMatch_2LTrigOR", "trigMatch_2LTrigOR", &trigMatch_2LTrigOR, "O");

    //---------------------------------------------
    // set channel and lepton flavors
    //---------------------------------------------

    int flavor;
    if (TString(channel).EqualTo("ee")) {
        flavor = 1;
        lepChannel = 1;
    }
    else if (TString(channel).EqualTo("mm")) {
        flavor = 2;
        lepChannel = 0;
    }

    lep_flavor->clear();
    lep_flavor->push_back(flavor);
    lep_flavor->push_back(flavor);

    //-----------------------------
    // get smearing histograms and perform smearing
    //-----------------------------

    bins::init_binning_histograms();
    vector<TH1D*> g_pt_smear_dist = GetSmearingDistribution(channel, period, data_or_mc);

    //-----------------------------
    // get Z lepton CM theta distribution
    //-----------------------------

    TH1F* h_lep_theta = GetLepThetaHistogram(period, channel, data_or_mc);

    std::vector<int> lep_theta_count;
    std::vector<float> cm_theta_bin_boundaries;
    for (int i=0; i<h_lep_theta->GetNbinsX(); i++) {
        lep_theta_count.push_back(h_lep_theta->GetBinContent(i));
        cm_theta_bin_boundaries.push_back(h_lep_theta->GetBinLowEdge(i));
    }
    cm_theta_bin_boundaries.push_back(h_lep_theta->GetBinLowEdge(h_lep_theta->GetNbinsX()) + h_lep_theta->GetBinWidth(h_lep_theta->GetNbinsX()));
    std::discrete_distribution<int> cm_theta_distribution (lep_theta_count.begin(),lep_theta_count.end());

    //-----------------------------
    // loop over events
    //-----------------------------

    unsigned random_seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine lep_theta_generator(random_seed);

    Long64_t nentries = inputTree->GetEntries();
    for (Long64_t i=0; i<nentries; i++) {

        if (fmod(i,1e5)==0) cout << i << " events processed." << endl;
        inputTree->GetEntry(i);

        //--- get random smearing values, then apply smearing (note signs)
        float gamma_pt_smear = 0;
        float gamma_phi_smear = 0;
        for (int i=0; i<3; i++) {
            int pt_smear_bin = bins::hist_pt_bins->FindBin(gamma_pt)-1;
            if (pt_smear_bin>=0 && g_pt_smear_dist[pt_smear_bin]->Integral()>0)
                gamma_pt_smear = g_pt_smear_dist[pt_smear_bin]->GetRandom();
        }

        gamma_pt_smeared = gamma_pt - gamma_pt_smear;
        gamma_phi_smeared = gamma_phi + gamma_phi_smear;

        TLorentzVector gamma_4vec, gamma_smeared_4vec, MET_4vec, MET_smeared_4vec;
        gamma_4vec.SetPtEtaPhiM(gamma_pt, 0, gamma_phi, 0);
        gamma_smeared_4vec.SetPtEtaPhiM(gamma_pt_smeared, 0, gamma_phi_smeared, 0);
        MET_4vec.SetPtEtaPhiM(MET_raw,0,MET_phi,0);
        MET_smeared_4vec = MET_4vec + gamma_4vec - gamma_smeared_4vec;

        MET_smeared = MET_smeared_4vec.Pt();
        DPhi_METPhoton_smear = gamma_phi_smeared - MET_smeared_4vec.Phi();
        METl_smeared = MET_smeared * TMath::Cos(DPhi_METPhoton_smear);
        METt_smeared = MET_smeared * TMath::Sin(DPhi_METPhoton_smear);

        int gamma_pt_smear_bin = bins::hist_pt_bins->FindBin(gamma_pt_smeared)-1;
        //if (gamma_pt_smeared>bins::pt_bins[bins::n_pt_bins]) gamma_pt_smear_bin = bins::n_pt_bins-1;
        int METl_bin = bins::hist_METl_bins->FindBin(METl_smeared)-1;
        mll = 91.1876;
        if (METl_bin>=0 && gamma_pt_smear_bin>=0) {
            if (hist_z_mll_bin_pt_metl[gamma_pt_smear_bin][METl_bin]->Integral()>0)
                mll = hist_z_mll_bin_pt_metl[gamma_pt_smear_bin][METl_bin]->GetRandom();
        }

        //---------------------------------------------
        // compute two lepton kinematics
        //---------------------------------------------
        TRandom1 myRandom;
        myRandom.SetSeed(0);

        TLorentzVector z_4vec;
        z_4vec.SetPtEtaPhiM(gamma_pt_smeared,gamma_eta,gamma_phi_smeared,mll);

        // boost along z axis (since we measure angles in CM relative to boost direction)
        TVector3 boost_vec_lab = z_4vec.BoostVector();
        TVector3 boost_vec(0, 0, boost_vec_lab.Mag());

        TLorentzVector l0_lab_4vec, l1_lab_4vec;
        while (true) {

            double lep_phi_cm = myRandom.Rndm()*2.*TMath::Pi();

            // Histogram sampling
            int lep_theta_bin = cm_theta_distribution(lep_theta_generator);
            float low_lep_theta = cm_theta_bin_boundaries[lep_theta_bin];
            float high_lep_theta = cm_theta_bin_boundaries[lep_theta_bin+1];
            lep_theta_cm = myRandom.Rndm()*(high_lep_theta-low_lep_theta) + low_lep_theta;

            // Split leptons in Z rest frame
            TLorentzVector l0_cm_4vec, l1_cm_4vec;
            double lep_E_cm = z_4vec.M()/2.;
            double lep_px_cm = lep_E_cm*TMath::Sin(lep_theta_cm)*TMath::Cos(lep_phi_cm);
            double lep_py_cm = lep_E_cm*TMath::Sin(lep_theta_cm)*TMath::Sin(lep_phi_cm);
            double lep_pz_cm = lep_E_cm*TMath::Cos(lep_theta_cm);
            l0_cm_4vec.SetPxPyPzE(lep_px_cm, lep_py_cm, lep_pz_cm, lep_E_cm);
            l1_cm_4vec.SetPxPyPzE(-lep_px_cm, -lep_py_cm, -lep_pz_cm, lep_E_cm);

            // Boost to lab frame using smeared photon pT, eta, and phi
            l0_lab_4vec = l0_cm_4vec;
            l1_lab_4vec = l1_cm_4vec;
            l0_lab_4vec.Boost(boost_vec);
            l1_lab_4vec.Boost(boost_vec);
            if (l0_lab_4vec.Pt() < l1_lab_4vec.Pt()) {
                TLorentzVector lep_placeholder = l1_lab_4vec;
                l1_lab_4vec = l0_lab_4vec;
                l0_lab_4vec = lep_placeholder;
            }

            // Rotate to lab coordinates
            l0_lab_4vec.RotateY(z_4vec.Theta());
            l0_lab_4vec.RotateZ(z_4vec.Phi());
            l1_lab_4vec.RotateY(z_4vec.Theta());
            l1_lab_4vec.RotateZ(z_4vec.Phi());

            // Select lepton flavor and charge
            int charge = myRandom.Integer(2)*2-1;

            // Add leptons to event
            lep_pT->clear();
            lep_pT->push_back(l0_lab_4vec.Pt());
            lep_pT->push_back(l1_lab_4vec.Pt());
            lep_eta->clear();
            lep_eta->push_back(l0_lab_4vec.Eta());
            lep_eta->push_back(l1_lab_4vec.Eta());
            lep_phi->clear();
            lep_phi->push_back(l0_lab_4vec.Phi());
            lep_phi->push_back(l1_lab_4vec.Phi());
            lep_charge->clear();
            lep_charge->push_back(charge);
            lep_charge->push_back(-charge);

            // Stop loop if we're ready
            if (lep_pT->at(0)>cuts::leading_lep_pt_cut && lep_pT->at(1)>cuts::second_lep_pt_cut) break;
        }

        lep_n = 2;
        MT2 = ComputeMT2(l0_lab_4vec, l1_lab_4vec, MET_smeared_4vec, 0, 0).Compute();
        DPhi_METLepLeading_smeared = fabs(MET_smeared_4vec.DeltaPhi(l0_lab_4vec));
        DPhi_METLepSecond_smeared = fabs(MET_smeared_4vec.DeltaPhi(l1_lab_4vec));
        DR_2Lep = l0_lab_4vec.DeltaR(l1_lab_4vec);

        BaselineTree->Fill();
    }

    //-----------------------------
    // write tree and histograms
    //-----------------------------

    outputFile->cd();
    BaselineTree->Write();

    cout << "done." << endl;
    delete outputFile;
}

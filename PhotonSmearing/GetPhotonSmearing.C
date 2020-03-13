#include "../Common/Settings.C"
#include "SmearingFunctions.C"

using namespace std;

void GetPhotonSmearing(string period, string channel, string data_or_mc, bool turn_off_shifting_and_smearing=false) {

    //---------------------------------------------
    // get unsmeared input file and create smeared output file
    //---------------------------------------------

    TString data_period = DataPeriod(period);
    TString mc_period = getMCPeriod(period);

    cout << "channel         " << channel         << endl;
    cout << "period          " << period          << endl;
    cout << "isData?         " << data_or_mc          << endl;
    cout << "smearing path   " << smearing_path   << endl;

    TH1::SetDefaultSumw2();

    TString infilename;
    if (data_or_mc == "MC") infilename = ntuple_path + "g_mc/" + mc_period + "_SinglePhoton222.root";
    else if (data_or_mc == "Data") infilename = ntuple_path + "g_data/" + data_period + "_photon.root";

    TChain* inputTree = new TChain("BaselineTree");
    inputTree->Add(infilename);

    cout << endl;
    cout << "Opening file           : " << infilename        << endl;
    cout << "Events in ntuple       : " << inputTree->GetEntries() << endl;

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
    float DPhi_METZPhoton_smear; BaselineTree->Branch("DPhi_METZPhoton", &DPhi_METZPhoton_smear, "DPhi_METZPhoton/F");
    float MT2; BaselineTree->Branch("mt2leplsp_0", &MT2, "mt2leplsp_0/F");
    float DR_2Lep; BaselineTree->Branch("DR_2Lep", &DR_2Lep, "DR_2Lep/F");
    int photon_conversion_type; CopyBranch(inputTree, BaselineTree, "PhotonConversionType", "PhotonConversionType", &photon_conversion_type, "I");
    float lep_theta_cm; BaselineTree->Branch("lep_theta_cm", &lep_theta_cm, "lep_theta_cm/F");

    //--- HistFitter branches
    CopyAllBranches(inputTree, BaselineTree, histFitterBranches);

    float mll; CopyBranch(inputTree, BaselineTree, "mll", "mll", &mll, "F");
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
    float MET_sig; CopyBranch(inputTree, BaselineTree, "MET_sig", "MET_sig", &MET_sig, "F");

    //---------------------------------------------
    // set diagnostics printing
    //---------------------------------------------

    bool diagnostics = false;

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
    // get smearing Gaussians
    //-----------------------------

    bins::init_binning_histograms();
    map<int, pair<float, float>> smearing_gaussians = GetSmearingDistribution(channel, period, data_or_mc, diagnostics);
    //vector<TH1D> smearing_hists = GetSmearingDistribution(channel, period, data_or_mc, diagnostics);

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
    std::default_random_engine random_generator(random_seed);

    TH1D* hist_g_smeared_metl_bin_pt[bins::n_pt_bins+2];
    for (int bin=0; bin<bins::n_pt_bins+2; bin++)
        hist_g_smeared_metl_bin_pt[bin] = new TH1D(TString("hist_g_smeared_metl_")+TString::Itoa(bin,10),"",bins::n_smearing_bins,bins::smearing_low,bins::smearing_high);

    Long64_t nentries = inputTree->GetEntries();
    for (Long64_t i=0; i<nentries; i++) {

        if (fmod(i,1e5)==0) cout << i << " events processed." << endl;
        inputTree->GetEntry(i);

        //--- get random smearing values, then apply smearing (note signs)
        float gamma_pt_smear = 0;
        float gamma_phi_smear = 0;
        gamma_pt_smeared = gamma_pt - gamma_pt_smear;
        gamma_phi_smeared = gamma_phi + gamma_phi_smear;

        int pt_bin = bins::hist_pt_bins->FindBin(gamma_pt);
        auto gaussian_vals = smearing_gaussians.find(pt_bin);
            
        if (turn_off_shifting_and_smearing)
            METl_smeared = METl;
        else {
            normal_distribution<float> smearing_gaussian(gaussian_vals->second.first, gaussian_vals->second.second);
            METl_smeared = METl + smearing_gaussian(random_generator);
        }
        METt_smeared = METt;
        MET_smeared = sqrt(pow(METl_smeared, 2) + pow(METt_smeared, 2));
        hist_g_smeared_metl_bin_pt[pt_bin]->Fill(METl_smeared);

        TLorentzVector MET_smeared_4vec;
        MET_smeared_4vec.SetPtEtaPhiM(MET_smeared,0,MET_phi,0);
        DPhi_METZPhoton_smear = gamma_phi_smeared - MET_smeared_4vec.Phi();

        int gamma_pt_smear_bin = bins::hist_pt_bins->FindBin(gamma_pt_smeared);
        int METl_bin = bins::hist_METl_bins->FindBin(METl_smeared);
        mll = 0;
        if (hist_z_mll_bin_pt_metl[gamma_pt_smear_bin][METl_bin]->Integral()>0)
            mll = hist_z_mll_bin_pt_metl[gamma_pt_smear_bin][METl_bin]->GetRandom();

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
        int n_tries = 0;
        while (n_tries++ < 10) {

            double lep_phi_cm = myRandom.Rndm()*2.*TMath::Pi();

            // Histogram sampling
            int lep_theta_bin = cm_theta_distribution(random_generator);
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
        if (n_tries==11) continue;

        lep_n = 2;
        MT2 = ComputeMT2(l0_lab_4vec, l1_lab_4vec, MET_smeared_4vec, 0, 0).Compute();
        DPhi_METLepLeading_smeared = fabs(MET_smeared_4vec.DeltaPhi(l0_lab_4vec));
        DPhi_METLepSecond_smeared = fabs(MET_smeared_4vec.DeltaPhi(l1_lab_4vec));
        DR_2Lep = l0_lab_4vec.DeltaR(l1_lab_4vec);

        BaselineTree->Fill();
    }

    //-----------------------------
    // Draw final METl distribution
    //-----------------------------

    if (diagnostics) {
        for (int pt_bin=0; pt_bin<bins::n_pt_bins+2; pt_bin++) {
            TCanvas *canvas = new TCanvas("canvas","canvas",600,600);
            canvas->cd();
            canvas->SetLogy();

            TH1D* g_smeared_hist = hist_g_smeared_metl_bin_pt[pt_bin];
            TString g_smeared_plot_name = Form("Plots/g_smeared_ptbin_%d.eps", pt_bin);
            g_smeared_hist->SetLineColor(1); g_smeared_hist->SetFillColor(42); g_smeared_hist->SetLineStyle(1);
            g_smeared_hist->GetXaxis()->SetTitle("METl");
            g_smeared_hist->GetYaxis()->SetTitle("entries / bin");
            g_smeared_hist->Draw("hist");
            canvas->Print(g_smeared_plot_name);

            delete canvas, g_smeared_hist;
        }
    }

    //-----------------------------
    // Write tree and histograms
    //-----------------------------

    outputFile->cd();
    BaselineTree->Write();

    cout << "done." << endl;
    delete outputFile;
}

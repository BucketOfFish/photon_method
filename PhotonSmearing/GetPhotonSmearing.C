#include "../Common/Settings.C"
#include "../Common/CommonLibraries.C"
#include "../Common/CommonFunctions.C"
#include "SmearingFunctions.C"

using namespace std;

TH1F* GetLepThetaHistogram(string period, string channel, string data_or_mc) {

    cout << "Getting Z lepton CM theta distribution histogram." << endl;
    gStyle->SetOptStat(0);

    //--- open files and create TChains
    string mc_period = period;
    if (mc_period == "data15-16") mc_period = "mc16a";
    else if (mc_period == "data17") mc_period = "mc16cd";
    else if (mc_period == "data18") mc_period = "mc16e";
    string data_filename = ntuple_path + "bkg_data/" + period + "_bkg.root";
    string tt_filename = ntuple_path + "bkg_mc/" + mc_period + "_ttbar.root";
    string vv_filename = ntuple_path + "bkg_mc/" + mc_period + "_diboson.root";
    string zjets_filename = ntuple_path + "bkg_mc/" + mc_period + "_Zjets.root";

    TChain *tch_data, *tch_tt, *tch_vv, *tch_zjets;

    if (data_or_mc == "Data") {
        cout << "Opening data file    " << data_filename << endl;
        tch_data = new TChain("BaselineTree"); tch_data->Add(data_filename.c_str());
        cout << "data entries         " << tch_data->GetEntries() << endl;
        cout << "Opening ttbar file   " << tt_filename << endl;
        tch_tt = new TChain("BaselineTree"); tch_tt->Add(tt_filename.c_str());
        cout << "ttbar entries        " << tch_tt->GetEntries() << endl;
        cout << "Opening diboson file " << vv_filename << endl;
        tch_vv = new TChain("BaselineTree"); tch_vv->Add(vv_filename.c_str());
        cout << "diboson entries      " << tch_vv->GetEntries() << endl;
    }
    else {
        cout << "Opening Z+jets file  " << zjets_filename << endl;
        tch_zjets = new TChain("BaselineTree"); tch_zjets->Add(zjets_filename.c_str());
        cout << "Z+jets entries       " << tch_zjets->GetEntries() << endl;
    }

    //--- modify event selections and weights
    if (TString(channel).EqualTo("ee")) cuts::bkg_baseline += cuts::ee;
    else if (TString(channel).EqualTo("mm")) cuts::bkg_baseline += cuts::mm;
    else {
        cout << "Unrecognized channel! quitting   " << channel << endl;
        exit(0);
    }

    if (TString(period).EqualTo("data17")){
        TCut RunRange = TCut("RunNumber < 348000");  
        cout << "Data17! adding cut " << RunRange.GetTitle() << endl;
        cuts::bkg_baseline *= RunRange;
    }

    cout << "Z selection          " << cuts::bkg_baseline.GetTitle() << endl;
    cout << "Z weight             " << cuts::bkg_weight.GetTitle() << endl;
    cout << "g selection          " << cuts::photon_baseline.GetTitle() << endl;
    cout << "g weight             " << cuts::photon_weight.GetTitle() << endl;

    //--- fill lep theta histograms
    TH1F* hdata  = new TH1F("hdata", "", 30, 0, 3);
    TH1F* htt    = new TH1F("htt", "", 30, 0, 3);
    TH1F* hvv    = new TH1F("hvv", "", 30, 0, 3);
    TH1F* hz     = new TH1F("hz", "", 30, 0, 3);
    TH1F* histoG = new TH1F("histoG", "", 30, 0, 3);

    if (data_or_mc == "Data") {
        tch_data->Draw("Z_cm_lep_theta>>hdata", cuts::bkg_baseline, "goff");
        tch_tt->Draw("Z_cm_lep_theta>>htt", cuts::bkg_baseline*cuts::bkg_weight, "goff");
        tch_vv->Draw("Z_cm_lep_theta>>hvv", cuts::bkg_baseline*cuts::bkg_weight, "goff");
        cout << "data integral        " << hdata->Integral() << endl;
        cout << "ttbar integral       " << htt->Integral() << endl;
        cout << "diboson integral     " << hvv->Integral() << endl;
    }
    else {
        tch_zjets->Draw("Z_cm_lep_theta>>hz", cuts::bkg_baseline*cuts::bkg_weight, "goff");
        cout << "Z+jets integral      " << hz->Integral() << endl;
    }

    //--- return lep theta histogram
    TH1F* histoZ;
    if (data_or_mc == "Data") {
        histoZ = (TH1F*) hdata->Clone("histoZ");
        histoZ->Add(htt, -1.0);
        histoZ->Add(hvv, -1.0);
    }
    else histoZ = (TH1F*) hz->Clone("histoZ");

    return histoZ;
}

void GetPhotonSmearing(string period, string channel, string data_or_mc, int smearing_method) {

    string label = "photon";
    if (data_or_mc == "MC") {
        label = "SinglePhoton222";
        if (period == "data15-16") period = "mc16a";
        else if (period == "data17") period = "mc16cd";
        else if (period == "data18") period = "mc16e";
    }

    cout << "channel         " << channel         << endl;
    cout << "period          " << period          << endl;
    cout << "isData?         " << data_or_mc          << endl;
    cout << "smearing path   " << smearing_path   << endl;
    cout << "smearing method " << smearing_method << endl;

    //---------------------------------------------
    // get unsmeared input file
    //---------------------------------------------

    TH1::SetDefaultSumw2();

    string infilename;
    if (data_or_mc == "MC") infilename = ntuple_path + "g_mc/" + period + "_" + label + ".root";
    else if (data_or_mc == "Data") infilename = ntuple_path + "g_data/" + period + "_" + label + ".root";

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
    if (smearing_method == 0) photon_tag = "NoSmear";
    if (smearing_method == 4) photon_tag = "McSmear";
    if (smearing_method == 5) photon_tag = "DataSmear";

    string outfilename;
    if (data_or_mc == "Data") outfilename = TString(smearing_path+"g_data/"+period+"_"+label+"_"+channel+"_"+photon_tag+".root"); 
    if (data_or_mc == "MC") outfilename = TString(smearing_path+"g_mc/"+period+"_"+label+"_"+channel+"_"+photon_tag+".root"); 

    TFile* f = new TFile(outfilename.c_str(), "recreate");          
    TTree* BaselineTree = new TTree("BaselineTree", "baseline tree");

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
    float METl_smear; BaselineTree->Branch("METl", &METl_smear, "METl/F");
    float METt_smear; BaselineTree->Branch("METt", &METt_smear, "METt/F");
    float HT; CopyBranch(inputTree, BaselineTree, "HT", "HT", &HT, "F");
    float MET_raw; CopyBranch(inputTree, BaselineTree, "met_Et_raw", "met_Et_raw", &MET_raw, "F");

    vector<float>* jet_pT = new vector<float>(10); CopyBranch(inputTree, BaselineTree, "jet_pT", "jet_pT", &jet_pT, "vector<float>");
    vector<float>* jet_eta = new vector<float>(10); CopyBranch(inputTree, BaselineTree, "jet_eta", "jet_eta", &jet_eta, "vector<float>");
    vector<float>* jet_phi = new vector<float>(10); CopyBranch(inputTree, BaselineTree, "jet_phi", "jet_phi", &jet_phi, "vector<float>");

    float gamma_pt_smear; BaselineTree->Branch("Ptll", &gamma_pt_smear, "Ptll/F");
    float gamma_phi_smear; BaselineTree->Branch("Z_phi", &gamma_phi_smear, "Z_phi/F");
    int lep_n; BaselineTree->Branch("lep_n", &lep_n, "lep_n/I"); BaselineTree->Branch("nLep_signal", &lep_n, "nLep_signal/I"); BaselineTree->Branch("nLep_base", &lep_n, "nLep_base/I");
    float MET_smear; BaselineTree->Branch("met_Et", &MET_smear, "met_Et/F");
    float DPhi_METLepLeading_smear; BaselineTree->Branch("DPhi_METLepLeading", &DPhi_METLepLeading_smear, "DPhi_METLepLeading/F");
    float DPhi_METLepSecond_smear; BaselineTree->Branch("DPhi_METLepSecond", &DPhi_METLepSecond_smear, "DPhi_METLepSecond/F");
    float DPhi_METPhoton_smear; BaselineTree->Branch("DPhi_METPhoton", &DPhi_METPhoton_smear, "DPhi_METPhoton/F");
    float MT2W; BaselineTree->Branch("MT2W", &MT2W, "MT2W/F");
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

    //-----------------------------
    // get smearing histograms
    //-----------------------------

    TH1D* hist_MET_bins = new TH1D("hist_MET_bins", "", bins::smearing_bin_size, bins::MET_bins); //hist_MET_bins->SetStats(0);
    TH1D* hist_METl_bins = new TH1D("hist_METl_bins", "", bins::METl_bin_size, bins::METl_bins);
    TH1D* hist_pt_bins = new TH1D("hist_pt_bins", "", bins::smearing_bin_size, bins::pt_bins);

    TH1D* hist_Mll_dPt[bins::smearing_bin_size][bins::METl_bin_size];
    for (int bin0=0; bin0<bins::smearing_bin_size; bin0++) {
        for (int bin1=0; bin1<bins::METl_bin_size; bin1++) {
            hist_Mll_dPt[bin0][bin1] = new TH1D(TString("hist_Mll_dPt_")+TString::Itoa(bin0,10)+TString("_")+TString::Itoa(bin1,10),"",bins::mll_bin_size,bins::mll_bin);
        }
    }
    if (smearing_method != 0) {

        cout << "Prepare Mll histograms..." << endl;
        cout << "Path is " << ntuple_path << endl;

        string filename = ntuple_path + "/ZMC16a/Zjets_merged_processed.root";
        cout << "Opening mll histo file : " << filename << endl;
        TFile fZ(filename.c_str());
        TTree* tZ = (TTree*)fZ.Get("BaselineTree");

        tZ->SetBranchStatus("*", 0);
        double totalWeight; SetInputBranch(tZ, "totalWeight", &totalWeight);
        float METl; SetInputBranch(tZ, "METl", &METl);
        int jet_n; SetInputBranch(tZ, "nJet30", &jet_n);
        vector<float>* tZ_lep_pT = new vector<float>(10); SetInputBranch(tZ, "lepPt", &tZ_lep_pT);
        int bjet_n; SetInputBranch(tZ, "bjet_n", &bjet_n);
        float Z_pt; SetInputBranch(tZ, "Ptll", &Z_pt);
        float mll; SetInputBranch(tZ, "mll", &mll);
        int channel; SetInputBranch(tZ, "channel", &channel);

        for (int entry=0; entry<tZ->GetEntries(); entry++) {
            tZ->GetEntry(entry);

            if( TString(channel).EqualTo("ee") && channel != 1 ) continue; // ee
            if( TString(channel).EqualTo("mm") && channel != 0 ) continue; // ee
            if (jet_n<2) continue;
            if (tZ_lep_pT->at(0) < cuts::leading_lep_pt_cut) continue;
            if (tZ_lep_pT->at(1) < cuts::second_lep_pt_cut) continue;
            int pt_bin = hist_pt_bins->FindBin(Z_pt)-1;
            int METl_bin = hist_METl_bins->FindBin(METl)-1;
            if (METl_bin>=0 && pt_bin>=0) hist_Mll_dPt[pt_bin][METl_bin]->Fill(mll,totalWeight);
        }

        fZ.Close();

        for (int bin0=0; bin0<bins::smearing_bin_size; bin0++) {
            for (int bin1=0; bin1<bins::METl_bin_size; bin1++) {
                int rebin = RebinHistogram(hist_Mll_dPt[bin0][bin1], 0);
            }
        }
    }

    if (smearing_method != 0) {
        cout << "Prepare smearing histograms..." << endl;
        cout << "smearing_method    " << smearing_method << endl;
        float lumi = GetLumi(period);
        GetSmearingHistogram(channel, lumi, period, smearing_method);
    }

    //-----------------------------
    // Get Z lepton CM theta distribution
    //-----------------------------

    TH1F* h_lep_theta = GetLepThetaHistogram(period, channel, data_or_mc);
    unsigned random_seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine lep_theta_generator(random_seed);
    std::vector<int> lep_theta_count;
    std::vector<float> lep_theta_boundaries;
    for (int i=0; i<h_lep_theta->GetNbinsX(); i++) {
        lep_theta_count.push_back(h_lep_theta->GetBinContent(i));
        lep_theta_boundaries.push_back(h_lep_theta->GetBinLowEdge(i));
    }
    lep_theta_boundaries.push_back(h_lep_theta->GetBinLowEdge(h_lep_theta->GetNbinsX()) + h_lep_theta->GetBinWidth(h_lep_theta->GetNbinsX()));
    std::discrete_distribution<int> lep_theta_distribution (lep_theta_count.begin(),lep_theta_count.end());

    //-----------------------------
    // something about convolution? who knows
    //-----------------------------

    TH1D* g_resp[bins::smearing_bin_size];
    TH1D* z_resp[bins::smearing_bin_size];
    TH1D* smear_fft_re[bins::smearing_bin_size];
    TH1D* smear_fft_im[bins::smearing_bin_size];
    TH1D* smear_fft_amp[bins::smearing_bin_size];
    TH1D* smear_raw[bins::smearing_bin_size];
    TH1D* smear_final[bins::smearing_bin_size];
    float shift[bins::smearing_bin_size];
    TH1D* g_metl_smear[bins::smearing_bin_size];
    TH1D* g_metl_smear_2j[bins::smearing_bin_size];
    for (int bin=0;bin<bins::smearing_bin_size;bin++) {
        g_metl_smear[bin] = new TH1D(TString("g_metl_smear_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        g_metl_smear[bin]->SetStats(0);
        g_metl_smear_2j[bin] = new TH1D(TString("g_metl_smear_2j_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        g_metl_smear_2j[bin]->SetStats(0);
    }

    if (smearing_method != 0) {
        TSpectrum pfinder;
        for (int bin=0;bin<bins::smearing_bin_size;bin++) {
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
                if (channel=="mm" && bins::pt_bins[bin]>=0) smear_cut = 7.;
                if (abs(smear_raw[bin]->GetBinCenter(i+1)-smear_mean)/smear_rms>smear_cut) {
                    smear_raw[bin]->SetBinContent(i+1,0.);
                    smear[i] = 0.;
                }
            }

            shift[bin] = -g_metl[bin]->GetMean();
        }

        for (int bin=0;bin<bins::smearing_bin_size;bin++) {
            smear_final[bin] = new TH1D(TString("smear_final_")+TString::Itoa(bin,10),"",500,-1000,1000);
            for (int i=0;i<500;i++) {
                int which_bin = smear_raw[bin]->FindBin(smear_final[bin]->GetBinCenter(i+1));
                smear_final[bin]->SetBinContent(i+1,smear_raw[bin]->GetBinContent(which_bin));
            }
        }
    }

    //-----------------------------
    // loop over events
    //-----------------------------

    Long64_t nentries = inputTree->GetEntries();
    for (Long64_t i=0; i<nentries; i++) {

        if (fmod(i,1e5)==0) cout << i << " events processed." << endl;
        inputTree->GetEntry(i);

        //--- get random photon smearing values
        float photon_pt_smear = 0;
        int pt_smear_bin = 0;
        if (smearing_method != 0) {
            pt_smear_bin = hist_pt_bins->FindBin(gamma_pt)-1;
            if (pt_smear_bin>=0) {
                if (smear_final[pt_smear_bin]->Integral()>0) photon_pt_smear = smear_final[pt_smear_bin]->GetRandom() + shift[pt_smear_bin];
                else photon_pt_smear = shift[pt_smear_bin];
            }
        }
        float photon_phi_smear = 0;

        //--- smear MET and photon variables (note signs in this section)
        gamma_pt_smear = gamma_pt - photon_pt_smear;
        gamma_phi_smear = gamma_phi - photon_phi_smear;
        float photon_pt_smear_l = gamma_pt-gamma_pt_smear*TMath::Cos(photon_phi_smear);
        float photon_pt_smear_t = -gamma_pt_smear*TMath::Sin(photon_phi_smear);
        METl_smear = METl + photon_pt_smear_l;
        METt_smear = METt + photon_pt_smear_t;
        MET_smear = pow(METl_smear*METl_smear+METt_smear*METt_smear,0.5);

        //--- recompute DPhi after smearing
        TLorentzVector met_4vec_smear;
        if (smearing_method != 0) {
            float METtx = METt*TMath::Cos(gamma_phi_smear+TMath::Pi()/2.);
            float METty = METt*TMath::Sin(gamma_phi_smear+TMath::Pi()/2.);
            float METlx_smear = METl_smear*TMath::Cos(gamma_phi_smear);
            float METly_smear = METl_smear*TMath::Sin(gamma_phi_smear);

            met_4vec_smear.SetXYZM(METtx+METlx_smear,METty+METly_smear,0,0);
            DPhi_METPhoton_smear = fabs(TMath::ATan2(METt,METl_smear));
            //int dphi_smear = hist_low_dphi->FindBin(DPhi_METPhoton_smear)-1;
            //if (dphi_smear>dphi_bin[smearing_bin_size]) dphi_smear = smearing_bin_size-1;

            if (gamma_pt>50. && jet_n==1) g_metl_smear[pt_smear_bin]->Fill(METl_smear,totalWeight);
            if (gamma_pt>50. && jet_n>=2) g_metl_smear_2j[pt_smear_bin]->Fill(METl_smear,totalWeight);
        }
        else {
            met_4vec_smear.SetXYZM(METl, METt, 0, 0);
            DPhi_METPhoton_smear = fabs(TMath::ATan2(METt, METl_smear));
        }

        //--- translate photon pT to dilepton sum pT, and compute HTincl for photon events
        if (smearing_method != 0) {
            int photon_pt_smear_bin = hist_pt_bins->FindBin(gamma_pt_smear)-1;
            if (gamma_pt_smear>bins::pt_bins[bins::smearing_bin_size]) photon_pt_smear_bin = bins::smearing_bin_size-1;
            int MET_bin = hist_MET_bins->FindBin(MET_smear)-1;
            if (MET_bin > bins::MET_bins[bins::smearing_bin_size]) MET_bin = bins::smearing_bin_size-1;

            float photon_2LPt = 0;
            //HTincl = HT + photon_2LPt;
            int METl_bin = hist_METl_bins->FindBin(METl_smear)-1;
            if (METl_bin>=0 && photon_pt_smear_bin>=0) {
                if (hist_Mll_dPt[photon_pt_smear_bin][METl_bin]->Integral()>0)
                    mll = hist_Mll_dPt[photon_pt_smear_bin][METl_bin]->GetRandom();
            }
        }
        else {
            mll = 91.1876;
        }

        //---------------------------------------------
        // set channel and lepton flavors
        //---------------------------------------------
        int flavor;
        if( TString(channel).EqualTo("ee")) {
            flavor = 1;
            lepChannel = 1;
        }
        else if( TString(channel).EqualTo("mm")) {
            flavor = 2;
            lepChannel = 0;
        }

        //---------------------------------------------
        // compute two lepton kinematics
        //---------------------------------------------
        TRandom1 myRandom;
        myRandom.SetSeed(0);
        bool create_pseudo_leptons = true;
        if (create_pseudo_leptons) {
            TLorentzVector z_4vec;
            z_4vec.SetPtEtaPhiM(gamma_pt,gamma_eta,gamma_phi,mll);
            //GetDijetVariables(z_4vec, met_4vec_smear, jet_pT, jet_eta, jet_phi, jet_m);

            // boost along z axis (since we measure angles in CM relative to boost direction)
            TVector3 boost_vec_lab = z_4vec.BoostVector();
            TVector3 boost_vec(0, 0, boost_vec_lab.Mag());

            TLorentzVector l0_lab_4vec, l1_lab_4vec;
            while (true) {

                double lep_phi_cm = myRandom.Rndm()*2.*TMath::Pi();

                //// Naive sampling (incorrect)
                //double lep_theta_cm = myRandom.Rndm()*TMath::Pi()-0.5*TMath::Pi();

                //// Uniform sampling
                //double lep_theta_cm = acos(1 - 2*myRandom.Rndm());

                //// Drell-Yan lepton angular distribution (with Z boost direction as +z)
                //double placeholder_1 = 4-8*myRandom.Rndm();
                //double placeholder_2 = pow(pow(placeholder_1,2)+4,1.0/2) + placeholder_1;
                //double numerator = pow(2.0,1/3)*pow(placeholder_2,2.0/3) - 2;
                //double denominator = pow(2.0,2/3)*pow(placeholder_2,1.0/3);
                //double lep_theta_cm = acos(numerator/denominator);

                //// Sin^3 sampling
                //lep_theta_cm = atanh(1.973926*(myRandom.Rndm() - 0.5))/1.6 + TMath::Pi()/2;

                // Histogram sampling
                int lep_theta_bin = lep_theta_distribution(lep_theta_generator);
                float low_lep_theta = lep_theta_boundaries[lep_theta_bin];
                float high_lep_theta = lep_theta_boundaries[lep_theta_bin+1];
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
                lep_flavor->clear();
                lep_flavor->push_back(flavor);
                lep_flavor->push_back(flavor);
                lep_charge->clear();
                lep_charge->push_back(charge);
                lep_charge->push_back(-charge);

                // Stop loop if we're ready
                if (lep_pT->at(0)>cuts::leading_lep_pt_cut && lep_pT->at(1)>cuts::second_lep_pt_cut) break;

                // Checks
                //TLorentzVector twol_cm_4vec = l0_cm_4vec + l1_cm_4vec;
                //TLorentzVector twol_lab_4vec = l0_lab_4vec + l1_lab_4vec;
                //cout << "z_4vec pT = " << z_4vec.Pt() << ", eta = " << z_4vec.Eta() << ", phi = " << z_4vec.Phi() << ", m = " << z_4vec.M() << endl;
                //cout << "l_pT_cm = " << l_pT_cm << ", min_theta = " << min_theta << ", phi = " << l_phi_cm << ", theta = " << l_theta_cm << endl;
                //cout << "l0_cm_4vec pT = " << l0_cm_4vec.Pt() << ", eta = " << l0_cm_4vec.Eta() << ", phi = " << l0_cm_4vec.Phi() << ", m = " << l0_cm_4vec.M() << endl;
                //cout << "l1_cm_4vec pT = " << l1_cm_4vec.Pt() << ", eta = " << l1_cm_4vec.Eta() << ", phi = " << l1_cm_4vec.Phi() << ", m = " << l1_cm_4vec.M() << endl;
                //cout << "2l_cm_4vec pT = " << twol_cm_4vec.Pt() << ", eta = " << twol_cm_4vec.Eta() << ", phi = " << twol_cm_4vec.Phi() << ", m = " << twol_cm_4vec.M() << endl;
                //cout << "l0_lab_4vec pT = " << l0_lab_4vec.Pt() << ", eta = " << l0_lab_4vec.Eta() << ", phi = " << l0_lab_4vec.Phi() << ", m = " << l0_lab_4vec.M() << endl;
                //cout << "l1_lab_4vec pT = " << l1_lab_4vec.Pt() << ", eta = " << l1_lab_4vec.Eta() << ", phi = " << l1_lab_4vec.Phi() << ", m = " << l1_lab_4vec.M() << endl;
                //cout << "2l_lab_4vec pT = " << twol_lab_4vec.Pt() << ", eta = " << twol_lab_4vec.Eta() << ", phi = " << twol_lab_4vec.Phi() << ", m = " << twol_lab_4vec.M() << endl;
                //cout << "==================================================================================" << endl;
            }

            lep_n = 2;
            MT2W = ComputeMT2(l0_lab_4vec, l1_lab_4vec, met_4vec_smear, 0, 0).Compute();
            DPhi_METLepLeading_smear = fabs(met_4vec_smear.DeltaPhi(l0_lab_4vec));
            DPhi_METLepSecond_smear = fabs(met_4vec_smear.DeltaPhi(l1_lab_4vec));
            DR_2Lep = l0_lab_4vec.DeltaR(l1_lab_4vec);
        }
        else {
            lep_n = 2;
            lep_pT->clear();
            lep_eta->clear();
            lep_phi->clear();
            lep_flavor->clear();
            lep_charge->clear();
            lep_pT->push_back(50);
            lep_pT->push_back(50);
            lep_eta->push_back(0);
            lep_eta->push_back(0);
            lep_phi->push_back(0);
            lep_phi->push_back(0);
            lep_flavor->push_back(11); // 11=electron, 13=muon
            lep_flavor->push_back(11);
            lep_charge->push_back(-1);
            lep_charge->push_back(1);
            MT2W = 0;
            DPhi_METLepLeading_smear = 0;
            DPhi_METLepSecond_smear = 0;
            DR_2Lep = 0;
        }

        BaselineTree->Fill();
    }

    //-----------------------------
    // write tree and histograms
    //-----------------------------

    BaselineTree->Write();

    cout << "done." << endl;
    delete f;
}

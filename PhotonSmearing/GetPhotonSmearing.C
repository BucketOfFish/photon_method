#include "../Common/Settings.C"
#include "../Common/CommonLibraries.C"
#include "../Common/CommonFunctions.C"
#include "SmearingFunctions.C"

using namespace std;

void GetPhotonSmearing(string label, string channel, string isData, string period, int smearing_method) {

    cout << "channel         " << channel         << endl;
    cout << "period          " << period          << endl;
    cout << "isData?         " << isData          << endl;
    cout << "smearing path   " << smearing_path   << endl;
    cout << "smearing method " << smearing_method << endl;

    //-----------------------------
    // get mll and smearing histograms
    //-----------------------------

    std::cout << "Prepare Mll histograms..." << std::endl;
    GetMllHistogram(channel,period); // hist_Mll_dPt 

    std::cout << "Prepare smearing histograms..." << std::endl;
    cout << "smearing_method    " << smearing_method << endl;
    float lumi = GetLumi(period);
    GetSmearingHistogram(channel, lumi, period, smearing_method);

    //---------------------------------------------
    // get unsmeared input file
    //---------------------------------------------

    TH1::SetDefaultSumw2();

    string  infilename;
    if (isData == "MC") infilename = ntuple_path + "gmc/" + label + ".root";
    else if (isData == "Data") infilename = ntuple_path + "gdata/" + label + ".root";

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
    if (isData == "Data") outfilename = TString(TString(smearing_path)+"gdata/" + label + "_"+TString(channel)+TString(photon_tag)+".root"); 
    if (isData == "MC") outfilename = TString(TString(smearing_path)+"gmc/gmc_"+TString(channel)+TString(photon_tag)+".root"); 

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
    float HT; CopyBranch(inputTree, BaselineTree, "HT", "HT", &HT, "F");
    float MET_raw; CopyBranch(inputTree, BaselineTree, "met_Et_raw", "met_Et_raw", &MET_raw, "F");

    std::vector<float>* jet_pT = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jet_pT", "jet_pT", &jet_pT, "std::vector<float>");
    std::vector<float>* jet_eta = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jet_eta", "jet_eta", &jet_eta, "std::vector<float>");
    std::vector<float>* jet_phi = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jet_phi", "jet_phi", &jet_phi, "std::vector<float>");

    float gamma_pt_smear; BaselineTree->Branch("Ptll", &gamma_pt_smear, "Ptll/F");
    float gamma_phi_smear; BaselineTree->Branch("Z_phi", &gamma_phi_smear, "Z_phi/F");
    int lep_n; BaselineTree->Branch("lep_n", &lep_n, "lep_n/I");
    float MET_smear; BaselineTree->Branch("met_Et", &MET_smear, "met_Et/F");
    float DPhi_METLepLeading_smear; BaselineTree->Branch("DPhi_METLepLeading", &DPhi_METLepLeading_smear, "DPhi_METLepLeading/F");
    float DPhi_METLepSecond_smear; BaselineTree->Branch("DPhi_METLepSecond", &DPhi_METLepSecond_smear, "DPhi_METLepSecond/F");
    float DPhi_METPhoton_smear; BaselineTree->Branch("DPhi_METPhoton", &DPhi_METPhoton_smear, "DPhi_METPhoton/F");
    float MT2W; BaselineTree->Branch("MT2W", &MT2W, "MT2W/F");
    float DR_2Lep; BaselineTree->Branch("DR_2Lep", &DR_2Lep, "DR_2Lep/F");

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

    std::vector<float>* jet_m = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jetM", "jetM", &jet_m, "std::vector<float>");
    Int_t jet_n; CopyBranch(inputTree, BaselineTree, "nJet30", "nJet30", &jet_n, "I");
    std::vector<int>* lepFlavor = new std::vector<int>(10); CopyBranch(inputTree, BaselineTree, "lepFlavor", "lepFlavor", &lepFlavor, "std::vector<int>");
    std::vector<int>* lepCharge = new std::vector<int>(10); CopyBranch(inputTree, BaselineTree, "lepCharge", "lepCharge", &lepCharge, "std::vector<int>");
    std::vector<float>* lep_pT = new std::vector<float>(10); BaselineTree->Branch("lepPt", "std::vector<float>", &lep_pT);
    int nBJet20_MV2c10_FixedCutBEff_77; CopyBranch(inputTree, BaselineTree, "nBJet20_MV2c10_FixedCutBEff_77", "nBJet20_MV2c10_FixedCutBEff_77", &nBJet20_MV2c10_FixedCutBEff_77, "I");
    Int_t nLep_signal; CopyBranch(inputTree, BaselineTree, "nLep_signal", "nLep_signal", &nLep_signal, "I");
    Int_t nLep_base; CopyBranch(inputTree, BaselineTree, "nLep_base", "nLep_base", &nLep_base, "I");
    bool trigMatch_2LTrigOR; CopyBranch(inputTree, BaselineTree, "trigMatch_2LTrigOR", "trigMatch_2LTrigOR", &trigMatch_2LTrigOR, "O");

    //-----------------------------
    // smearing histograms
    //-----------------------------

    TH1D* hist_low_met = new TH1D("hist_low_met", "", bin_size, met_bin); hist_low_met->SetStats(0);
    TH1D* g_resp[bin_size];
    TH1D* z_resp[bin_size];
    TH1D* smear_fft_re[bin_size];
    TH1D* smear_fft_im[bin_size];
    TH1D* smear_fft_amp[bin_size];
    TH1D* smear_raw[bin_size];
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
            if (channel=="mm" && sm_pt_bin[bin]>=0) smear_cut = 7.;
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
        MET_smear = pow(METl_smear*METl_smear+METt_smear*METt_smear,0.5);

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

        DPhi_METPhoton_smear = fabs(TMath::ATan2(METt,METl_smear));
        //int dphi_smear = hist_low_dphi->FindBin(DPhi_METPhoton_smear)-1;
        //if (dphi_smear>dphi_bin[bin_size]) dphi_smear = bin_size-1;

        if (gamma_pt>50. && jet_n==1) g_metl_smear[smpt]->Fill(METl_smear,totalWeight);
        if (gamma_pt>50. && jet_n>=2) g_metl_smear_2j[smpt]->Fill(METl_smear,totalWeight);

        // translate photon pT to dilepton sum pT, and compute HTincl for photon events
        float photon_2LPt = 0;
        //HTincl = HT + photon_2LPt;
        int dpt = hist_low_dpt->FindBin(METl_smear)-1;
        if (dpt>=0 && pt_smear>=0) {
            if (hist_Mll_dPt[pt_smear][dpt]->Integral()>0)
                mll = hist_Mll_dPt[pt_smear][dpt]->GetRandom();
        }

        // CHECKPOINT
        if (mll <= 0) // BIG HACK
            mll = 90; // BIG HACK

        //---------------------------------------------
        // compute two lepton kinematics
        //---------------------------------------------
        TLorentzVector z_4vec;
        z_4vec.SetPtEtaPhiM(gamma_pt,gamma_eta,gamma_phi,mll);
        //GetDijetVariables(z_4vec, met_4vec_smear, jet_pT, jet_eta, jet_phi, jet_m);

        lep_pT->clear();
        lep_eta->clear();
        lep_phi->clear();
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

    std::cout << "done." << std::endl;
    delete f;
}

#include "../Settings.C"
#include "../CommonFunctions/CommonLibraries.C"
#include "../CommonFunctions/CommonFunctions.C"
#include "GetDijetVariables.C"
#include "GetSmearingHistogram.C"
#include "GetMllHistogram.C"
//#include "GetIndividualLeptonInfo.C"
#include "MT2.h"

using namespace std;

vector<string> noSampleWeight;
vector<string> noEventWeight;

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

//  period: data15-16 (input) -> ZMC16a (source file), data17 -> ZMC16cd, data18 -> ZMC16cd
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

    GetMllHistogram(ch,period); // fill histogram hist_Mll_dPt

    for (int bin0=0; bin0<bin_size; bin0++) {
        for (int bin1=0; bin1<dpt_bin_size; bin1++) {
            int rebin = RebinHistogram(hist_Mll_dPt[bin0][bin1], 0);
        }
    }

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
    float MET; inputTree->SetBranchAddress("MET_raw", &MET);
    float DPhi_METJetLeading; inputTree->SetBranchAddress("DPhi_METJetLeading_raw", &DPhi_METJetLeading);
    float DPhi_METJetSecond; inputTree->SetBranchAddress("DPhi_METJetSecond_raw", &DPhi_METJetSecond);
    float MinDPhi_PhotonJet; inputTree->SetBranchAddress("MinDPhi_PhotonJet", &MinDPhi_PhotonJet);

    double totalWeight; CopyBranch(inputTree, BaselineTree, "totalWeight", "totalWeight", &totalWeight, "D");
    float HT; CopyBranch(inputTree, BaselineTree, "HT", "HT", &HT, "F");
    int jet_n; CopyBranch(inputTree, BaselineTree, "jet_n", "jet_n", &jet_n, "I");
    int bjet_n; CopyBranch(inputTree, BaselineTree, "bjet_n", "bjet_n", &bjet_n, "I");
    float gamma_pt; CopyBranch(inputTree, BaselineTree, "gamma_pt", "gamma_pt",  &gamma_pt, "F");
    float gamma_eta; CopyBranch(inputTree, BaselineTree, "gamma_eta", "Z_eta",  &gamma_eta, "F");
    float MET_loose; CopyBranch(inputTree, BaselineTree, "MET_loose", "MET_loose", &MET_loose, "F");
    float MET_tighter; CopyBranch(inputTree, BaselineTree, "MET_tighter", "MET_tighter", &MET_tighter, "F");
    float MET_tenacious; CopyBranch(inputTree, BaselineTree, "MET_tenacious", "MET_tenacious", &MET_tenacious, "F");
    int is2Lep2Jet; CopyBranch(inputTree, BaselineTree, "is2Lep2Jet", "is2Lep2Jet", &is2Lep2Jet, "I");
    int is2L2JInt; CopyBranch(inputTree, BaselineTree, "is2L2JInt", "is2L2JInt", &is2L2JInt, "I");
    int nBJet20_MV2c10_FixedCutBEff_77; CopyBranch(inputTree, BaselineTree, "nBJet20_MV2c10_FixedCutBEff_77", "nBJet20_MV2c10_FixedCutBEff_77", &nBJet20_MV2c10_FixedCutBEff_77, "I");
    float mjj; CopyBranch(inputTree, BaselineTree, "mjj", "mjj", &mjj, "F");
    float mll_RJ; CopyBranch(inputTree, BaselineTree, "mll_RJ", "mll_RJ", &mll_RJ, "F");
    float R_minH2P_minH3P; CopyBranch(inputTree, BaselineTree, "R_minH2P_minH3P", "R_minH2P_minH3P", &R_minH2P_minH3P, "F");
    float RPT_HT5PP; CopyBranch(inputTree, BaselineTree, "RPT_HT5PP", "RPT_HT5PP", &RPT_HT5PP, "F");
    float dphiVP; CopyBranch(inputTree, BaselineTree, "dphiVP", "dphiVP", &dphiVP, "F");
    float H2PP; CopyBranch(inputTree, BaselineTree, "H2PP", "H2PP", &H2PP, "F");
    float H5PP; CopyBranch(inputTree, BaselineTree, "H5PP", "H5PP", &H5PP, "F");
    int nJet20; CopyBranch(inputTree, BaselineTree, "nJet20", "nJet20", &nJet20, "I");
    float minDphi; CopyBranch(inputTree, BaselineTree, "minDphi", "minDphi", &minDphi, "F");
    float MZ; CopyBranch(inputTree, BaselineTree, "MZ", "MZ", &MZ, "F");
    int NjS; CopyBranch(inputTree, BaselineTree, "NjS", "NjS", &NjS, "I");
    int NjISR; CopyBranch(inputTree, BaselineTree, "NjISR", "NjISR", &NjISR, "I");
    float dphiISRI; CopyBranch(inputTree, BaselineTree, "dphiISRI", "dphiISRI", &dphiISRI, "F");
    float RISR; CopyBranch(inputTree, BaselineTree, "RISR", "RISR", &RISR, "F");
    float PTISR; CopyBranch(inputTree, BaselineTree, "PTISR", "PTISR", &PTISR, "F");
    float PTI; CopyBranch(inputTree, BaselineTree, "PTI", "PTI", &PTI, "F");
    float PTCM; CopyBranch(inputTree, BaselineTree, "PTCM", "PTCM", &PTCM, "F");
    float MJ; CopyBranch(inputTree, BaselineTree, "MJ", "MJ", &MJ, "F");
    int is3Lep3Jet; CopyBranch(inputTree, BaselineTree, "is3Lep3Jet", "is3Lep3Jet", &is3Lep3Jet, "I");
    int is4Lep3Jet; CopyBranch(inputTree, BaselineTree, "is4Lep3Jet", "is4Lep3Jet", &is4Lep3Jet, "I");
    int lept1sign_VR; CopyBranch(inputTree, BaselineTree, "lept1sign_VR", "lept1sign_VR", &lept1sign_VR, "I");
    int lept2sign_VR; CopyBranch(inputTree, BaselineTree, "lept2sign_VR", "lept2sign_VR", &lept2sign_VR, "I");
    float lept1Pt_VR; CopyBranch(inputTree, BaselineTree, "lept1Pt_VR", "lept1Pt_VR", &lept1Pt_VR, "F");
    float lept2Pt_VR; CopyBranch(inputTree, BaselineTree, "lept2Pt_VR", "lept2Pt_VR", &lept2Pt_VR, "F");
    float MZ_VR; CopyBranch(inputTree, BaselineTree, "MZ_VR", "MZ_VR", &MZ_VR, "F");
    float MJ_VR; CopyBranch(inputTree, BaselineTree, "MJ_VR", "MJ_VR", &MJ_VR, "F");
    float RISR_VR; CopyBranch(inputTree, BaselineTree, "RISR_VR", "RISR_VR", &RISR_VR, "F");
    float PTISR_VR; CopyBranch(inputTree, BaselineTree, "PTISR_VR", "PTISR_VR", &PTISR_VR, "F");
    float PTI_VR; CopyBranch(inputTree, BaselineTree, "PTI_VR", "PTI_VR", &PTI_VR, "F");
    float PTCM_VR; CopyBranch(inputTree, BaselineTree, "PTCM_VR", "PTCM_VR", &PTCM_VR, "F");
    float dphiISRI_VR; CopyBranch(inputTree, BaselineTree, "dphiISRI_VR", "dphiISRI_VR", &dphiISRI_VR, "F");
    float METl; CopyBranch(inputTree, BaselineTree, "METl_raw", "METl", &METl, "F");
    float METt; CopyBranch(inputTree, BaselineTree, "METt_raw", "METt", &METt, "F");
    float MET_phi; CopyBranch(inputTree, BaselineTree, "MET_phi_raw", "MET_phi", &MET_phi, "F");

    std::vector<int>* lepFlavor = new std::vector<int>(10); CopyBranch(inputTree, BaselineTree, "lepFlavor", "lepFlavor", &lepFlavor, "std::vector<int>");
    std::vector<int>* lepCharge = new std::vector<int>(10); CopyBranch(inputTree, BaselineTree, "lepCharge", "lepCharge", &lepCharge, "std::vector<int>");
    std::vector<float>* jet_pT = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jetPt", "jet_pT", &jet_pT, "std::vector<float>");
    std::vector<float>* jet_eta = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jetEta", "jet_eta", &jet_eta, "std::vector<float>");
    std::vector<float>* jet_phi = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jetPhi", "jet_phi", &jet_phi, "std::vector<float>");
    std::vector<float>* jet_m = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jetM", "jet_m", &jet_m, "std::vector<float>");

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
    float DR_Wmin2Jet; BaselineTree->Branch("DR_Wmin2Jet", &DR_Wmin2Jet, "DR_Wmin2Jet/F");
    float DR_J0J1; BaselineTree->Branch("DR_J0J1", &DR_J0J1, "DR_J0J1/F");
    float mWmin; BaselineTree->Branch("mWmin", &mWmin, "mWmin/F");
    float Wmin_pt; BaselineTree->Branch("Wmin_pt", &Wmin_pt, "Wmin_pt/F");
    float Wmin_eta; BaselineTree->Branch("Wmin_eta", &Wmin_eta, "Wmin_eta/F");
    float DPhi_METWmin; BaselineTree->Branch("DPhi_METWmin", &DPhi_METWmin, "DPhi_METWmin/F");
    float DPhi_WminZ; BaselineTree->Branch("DPhi_WminZ", &DPhi_WminZ, "DPhi_WminZ/F");
    float mj0j1; BaselineTree->Branch("mj0j1", &mj0j1, "mj0j1/F");
    float W01_pt; BaselineTree->Branch("W01_pt", &W01_pt, "W01_pt/F");
    float DPhi_METW01; BaselineTree->Branch("DPhi_METW01", &DPhi_METW01, "DPhi_METW01/F");
    float DPhi_W01Z; BaselineTree->Branch("DPhi_W01Z", &DPhi_W01Z, "DPhi_W01Z/F");
    float DPhi_METNonWJet; BaselineTree->Branch("DPhi_METNonWJet", &DPhi_METNonWJet, "DPhi_METNonWJet/F");
    float NonWJet_pT; BaselineTree->Branch("NonWJet_pT", &NonWJet_pT, "NonWJet_pT/F");
    float DPhi_METNonWminJet; BaselineTree->Branch("DPhi_METNonWminJet", &DPhi_METNonWminJet, "DPhi_METNonWminJet/F");
    float NonWminJet_pT; BaselineTree->Branch("NonWminJet_pT", &NonWminJet_pT, "NonWminJet_pT/F");
    float MT2W; BaselineTree->Branch("MT2W", &MT2W, "MT2W/F");
    float DR_2Lep; BaselineTree->Branch("DR_2Lep", &DR_2Lep, "DR_2Lep/F");

    TH1D* hist_low_njet = new TH1D("hist_low_njet", "", bin_size, njet_bin); hist_low_njet->SetStats(0);
    TH1D* hist_low_nbjet = new TH1D("hist_low_nbjet", "", bin_size, njet_bin); hist_low_nbjet->SetStats(0);
    TH1D* hist_low_pt = new TH1D("hist_low_pt", "", bin_size, sm_pt_bin); hist_low_pt->SetStats(0);
    TH1D* hist_sm_pt = new TH1D("hist_sm_pt", "", bin_size, sm_pt_bin); hist_sm_pt->SetStats(0);
    TH1D* hist_low_et = new TH1D("hist_low_et", "", bin_size, et_bin); hist_low_et->SetStats(0);
    TH1D* hist_low_ht = new TH1D("hist_low_ht", "", bin_size, ht_bin); hist_low_ht->SetStats(0);
    TH1D* hist_medium_njet = new TH1D("hist_medium_njet", "", bin_size, njet_bin); hist_medium_njet->SetStats(0);
    TH1D* hist_medium_nbjet = new TH1D("hist_medium_nbjet", "", bin_size, njet_bin); hist_medium_nbjet->SetStats(0);
    TH1D* hist_medium_pt = new TH1D("hist_medium_pt", "", bin_size, sm_pt_bin); hist_medium_pt->SetStats(0);
    TH1D* hist_medium_et = new TH1D("hist_medium_et", "", bin_size, et_bin); hist_medium_et->SetStats(0);
    TH1D* hist_medium_ht = new TH1D("hist_medium_ht", "", bin_size, ht_bin); hist_medium_ht->SetStats(0);
    TH1D* hist_high_njet = new TH1D("hist_high_njet", "", bin_size, njet_bin); hist_high_njet->SetStats(0);
    TH1D* hist_high_nbjet = new TH1D("hist_high_nbjet", "", bin_size, njet_bin); hist_high_nbjet->SetStats(0);
    TH1D* hist_high_pt = new TH1D("hist_high_pt", "", bin_size, sm_pt_bin); hist_high_pt->SetStats(0);
    TH1D* hist_high_et = new TH1D("hist_high_et", "", bin_size, et_bin); hist_high_et->SetStats(0);
    TH1D* hist_high_ht = new TH1D("hist_high_ht", "", bin_size, ht_bin); hist_high_ht->SetStats(0);
    TH1D* hist_zmet_njet = new TH1D("hist_zmet_njet", "", bin_size, njet_bin); hist_zmet_njet->SetStats(0);
    TH1D* hist_zmet_nbjet = new TH1D("hist_zmet_nbjet", "", bin_size, njet_bin); hist_zmet_nbjet->SetStats(0);
    TH1D* hist_zmet_pt = new TH1D("hist_zmet_pt", "", bin_size, sm_pt_bin); hist_zmet_pt->SetStats(0);
    TH1D* hist_zmet_et = new TH1D("hist_zmet_et", "", bin_size, et_bin); hist_zmet_et->SetStats(0);
    TH1D* hist_zmet_ht = new TH1D("hist_zmet_ht", "", bin_size, ht_bin); hist_zmet_ht->SetStats(0);
    TH1D* hist_bveto_njet = new TH1D("hist_bveto_njet", "", bin_size, njet_bin); hist_bveto_njet->SetStats(0);
    TH1D* hist_bveto_nbjet = new TH1D("hist_bveto_nbjet", "", bin_size, njet_bin); hist_bveto_nbjet->SetStats(0);
    TH1D* hist_bveto_pt = new TH1D("hist_bveto_pt", "", bin_size, sm_pt_bin); hist_bveto_pt->SetStats(0);
    TH1D* hist_bveto_et = new TH1D("hist_bveto_et", "", bin_size, et_bin); hist_bveto_et->SetStats(0);
    TH1D* hist_bveto_ht = new TH1D("hist_bveto_ht", "", bin_size, ht_bin); hist_bveto_ht->SetStats(0);
    TH1D* hist_low_dpt = new TH1D("hist_low_dpt", "", dpt_bin_size, dpt_bin); hist_low_dpt->SetStats(0);
    TH1D* hist_low_pt_smear = new TH1D("hist_low_pt_smear", "", bin_size, sm_pt_bin); hist_low_pt_smear->SetStats(0);
    TH1D* hist_low_htincl = new TH1D("hist_low_htincl", "", bin_size, ht_bin); hist_low_htincl->SetStats(0);
    TH1D* hist_medium_pt_smear = new TH1D("hist_medium_pt_smear", "", bin_size, sm_pt_bin); hist_medium_pt_smear->SetStats(0);
    TH1D* hist_medium_htincl = new TH1D("hist_medium_htincl", "", bin_size, ht_bin); hist_medium_htincl->SetStats(0);
    TH1D* hist_high_pt_smear = new TH1D("hist_high_pt_smear", "", bin_size, sm_pt_bin); hist_high_pt_smear->SetStats(0);
    TH1D* hist_high_htincl = new TH1D("hist_high_htincl", "", bin_size, ht_bin); hist_high_htincl->SetStats(0);
    TH1D* hist_zmet_pt_smear = new TH1D("hist_zmet_pt_smear", "", bin_size, sm_pt_bin); hist_zmet_pt_smear->SetStats(0);
    TH1D* hist_zmet_htincl = new TH1D("hist_zmet_htincl", "", bin_size, ht_bin); hist_zmet_htincl->SetStats(0);
    TH1D* hist_bveto_pt_smear = new TH1D("hist_bveto_pt_smear", "", bin_size, sm_pt_bin); hist_bveto_pt_smear->SetStats(0);
    TH1D* hist_bveto_htincl = new TH1D("hist_bveto_htincl", "", bin_size, ht_bin); hist_bveto_htincl->SetStats(0);
    TH1D* hist_low_met = new TH1D("hist_low_met", "", bin_size, met_bin); hist_low_met->SetStats(0);
    TH1D* hist_low_dphi = new TH1D("hist_low_dphi", "", bin_size, dphi_bin); hist_low_dphi->SetStats(0);

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
            //GetIndividualLeptonInfo(z_4vec);
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

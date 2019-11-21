#include "../Common/Settings.C"
#include "MT2.h"

TH1D* hist_z_metl[bins::smearing_bin_size];
TH1D* hist_g_metl[bins::smearing_bin_size];
TH1D* hist_z_onshell_metl[bins::smearing_bin_size];
TH1D* hist_z_mll_pt[bins::smearing_bin_size][bins::METl_bin_size];

void FillHistograms(string target_channel, TString period, int smearing_method) {

    /**
     * Fills global histograms hist_z_metl, hist_z_onshell_metl, hist_g_metl, and hist_z_mll_pt (2D)
     * by getting Z or photon data and entering the relevant features into binned histograms.
     */

    cout << "FillHistograms : smearing_method " << smearing_method << endl;

    for (int bin=0; bin<bins::smearing_bin_size; bin++) {
        hist_z_metl[bin] = new TH1D(TString("hist_z_metl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        hist_z_onshell_metl[bin] = new TH1D(TString("hist_z_jetmetl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        hist_g_metl[bin] = new TH1D(TString("hist_g_metl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
    }

    for (int bin0=0; bin0<bins::smearing_bin_size; bin0++) {
        for (int bin1=0; bin1<bins::METl_bin_size; bin1++) {
            hist_z_mll_pt[bin0][bin1] = new TH1D(TString("hist_z_Mll_dPt_")+TString::Itoa(bin0,10)+TString("_")+TString::Itoa(bin1,10),"",bins::mll_bin_size,bins::mll_bin);
        }
    }

    TString mc_period = MCPeriod(period);
    TString data_period = DataPeriod(period);

    //------------------------------------
    // SMEARING METHOD 0: no smearing
    //------------------------------------

    if (smearing_method == 0) return;

    //------------------------------------
    // SMEARING METHOD 4/5: for R21 MC/data
    //------------------------------------

    std::cout << "Getting Z histograms binned by pT and METl." << std::endl;

    TFile* data_file = new TFile(ntuple_path + "/bkg_data/" + data_period + "_bkg.root");
    TFile* ttbar_mc_file = new TFile(ntuple_path + "/bkg_mc/" + mc_period + "_ttbar.root");
    TFile* diboson_mc_file = new TFile(ntuple_path + "/bkg_mc/" + mc_period + "_diboson.root");
    TFile* Z_mc_file = new TFile(ntuple_path + "/bkg_mc/" + mc_period + "_Zjets.root");

    vector<TFile*> contributing_files;
    vector<int> file_weights;
    if (smearing_method == 4) {
        contributing_files = {Z_mc_file};
        file_weights = {1};
    }
    else if (smearing_method == 5) {
        contributing_files = {data_file, ttbar_mc_file, diboson_mc_file};
        file_weights = {1, -1, -1};
    }

    for (int i=0; i<contributing_files.size(); i++) {
        TTree* tree = (TTree*)contributing_files[i]->Get("BaselineTree");
        int fileWeight = file_weights[i];

        tree->SetBranchStatus("*", 0);
        double totalWeight; SetInputBranch(tree, "totalWeight", &totalWeight);
        int jet_n; SetInputBranch(tree, "nJet30", &jet_n);
        int bjet_n; SetInputBranch(tree, "bjet_n", &bjet_n);
        float ptll; SetInputBranch(tree, "Ptll", &ptll);
        float mll; SetInputBranch(tree, "mll", &mll);
        float METl; SetInputBranch(tree, "METl", &METl);
        int channel; SetInputBranch(tree, "channel", &channel);
        vector<float>* lep_pT = new vector<float>(10); SetInputBranch(tree, "lepPt", &lep_pT);

        for (int entry=0; entry<tree->GetEntries(); entry++) {
            tree->GetEntry(entry);
            if (TString(target_channel).EqualTo("ee") && channel != 1) continue;
            if (TString(target_channel).EqualTo("mm") && channel != 0) continue;
            if (ptll<50. || jet_n<2 || lep_pT->at(0)<cuts::leading_lep_pt_cut || lep_pT->at(1)<cuts::second_lep_pt_cut) continue;
            int pt_bin = bins::hist_pt_bins->FindBin(ptll)-1;
            int METl_bin = bins::hist_METl_bins->FindBin(METl)-1;
            hist_z_metl[pt_bin]->Fill(METl, fileWeight*totalWeight);
            if (mll>90 && mll<92) hist_z_onshell_metl[pt_bin]->Fill(METl, fileWeight*totalWeight);
            if (METl_bin>=0 && pt_bin>=0) hist_z_mll_pt[pt_bin][METl_bin]->Fill(mll, fileWeight*totalWeight);
        }
    }

    data_file->Close();
    ttbar_mc_file->Close();
    diboson_mc_file->Close();
    Z_mc_file->Close();

    //-- photon histogram

    std::cout << "Getting photon histogram binned by pT." << std::endl;

    TFile* photon_file;
    if (smearing_method == 4)
        photon_file = new TFile(ntuple_path + "/g_mc/" + mc_period + "_SinglePhoton222.root");
    else if (smearing_method == 5)
        photon_file = new TFile(ntuple_path + "/g_data/" + data_period + "_photon.root");

    TTree* tree = (TTree*)photon_file->Get("BaselineTree");
    tree->SetBranchStatus("*", 0);
    double totalWeight; SetInputBranch(tree, "totalWeight", &totalWeight);
    int jet_n; SetInputBranch(tree, "nJet30", &jet_n);
    int bjet_n; SetInputBranch(tree, "bjet_n", &bjet_n);
    float ptll; SetInputBranch(tree, "gamma_pt", &ptll);
    float METl; SetInputBranch(tree, "METl_raw", &METl);

    for (int entry=0; entry<tree->GetEntries(); entry++) {
        tree->GetEntry(entry);
        if (ptll<50. || jet_n!=1 || bjet_n!=0) continue;
        int pt_bin = bins::hist_pt_bins->FindBin(ptll)-1;
        hist_g_metl[pt_bin]->Fill(METl, totalWeight);
    }

    photon_file->Close();

}

TH1F* GetLepThetaHistogram(string period, string channel, string data_or_mc) {

    cout << "Getting Z lepton CM theta distribution histogram." << endl;
    gStyle->SetOptStat(0);

    //--- open files and create TChains
    string data_filename = ntuple_path + "bkg_data/" + period + "_bkg.root";
    if (period == "data15-16") period = "mc16a";
    else if (period == "data17") period = "mc16cd";
    else if (period == "data18") period = "mc16e";
    string tt_filename = ntuple_path + "bkg_mc/" + period + "_ttbar.root";
    string vv_filename = ntuple_path + "bkg_mc/" + period + "_diboson.root";
    string zjets_filename = ntuple_path + "bkg_mc/" + period + "_Zjets.root";

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

int RebinHistogram(TH1D* hist, int rebin) {

    /**
     * If rebin=0, rebin histogram by factors of 2, up to 32, until it looks reasonable or there are an odd number of bins.
     * This means having abs(negative_yield/positive_yield)<0.005 and core_yield/positive_yield<0.4.
     * Else, rebin by #(rebin).
     * Finally, set all negative bins of the histogram to 0.
     * Return the final rebinning factor for the histogram.
     */

    float negative_yield = 0.; // sum of bins with negative values
    float positive_yield = 0.; // sum of bins with positive values
    float core_yield = 0.; // sum of bins within one sigma of axis center

    if (rebin==0) {
        rebin = 1;
        for (int bin=1;bin<=hist->GetNbinsX();bin++) {
            if (hist->GetBinContent(bin)>=0) positive_yield += hist->GetBinContent(bin);
            else negative_yield += hist->GetBinContent(bin);
        }
        int has_odd_nbins = hist->GetNbinsX() % 2;
        while ((abs(negative_yield/positive_yield)>0.005 || core_yield/positive_yield<0.4) && has_odd_nbins==0 && rebin<=32) {
            hist->Rebin(2);
            rebin = rebin*2;
            has_odd_nbins = hist->GetNbinsX() % 2;
            negative_yield = 0.;
            positive_yield = 0.;
            core_yield = 0.;
            for (int bin=1;bin<=hist->GetNbinsX();bin++) {
                if (hist->GetBinContent(bin)>=0) positive_yield += hist->GetBinContent(bin);
                else negative_yield += hist->GetBinContent(bin);
                if (abs(hist->GetBinCenter(bin)-hist->GetMean())<hist->GetRMS()) { // if this bin is within a sigma of the center of the axis
                    core_yield += hist->GetBinContent(bin); // core_yield = 68% for a perfect Guassian
                }
            }
        }
    }
    else {
        // merge every #rebin bins into a single bin
        hist->Rebin(rebin);
    }

    for (int bin=1;bin<=hist->GetNbinsX();bin++) {
        hist->SetBinContent(bin,max(hist->GetBinContent(bin),0.));
    }

    return rebin;
}

TH1D* smear_final[bins::smearing_bin_size];
float shift[bins::smearing_bin_size];
TH1D* hist_g_metl_smear[bins::smearing_bin_size];
TH1D* hist_g_metl_smear_2j[bins::smearing_bin_size];

void ConvolveAndSmear(string channel, int smearing_method) {

    TH1D* smear_raw[bins::smearing_bin_size];
    for (int bin=0;bin<bins::smearing_bin_size;bin++) {
        hist_g_metl_smear[bin] = new TH1D(TString("hist_g_metl_smear_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        hist_g_metl_smear[bin]->SetStats(0);
        hist_g_metl_smear_2j[bin] = new TH1D(TString("hist_g_metl_smear_2j_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        hist_g_metl_smear_2j[bin]->SetStats(0);
    }

    //------------------------------------
    // SMEARING METHOD 0: no smearing
    //------------------------------------

    if (smearing_method == 0) return;

    //------------------------------------
    // SMEARING METHOD 4/5: for R21 MC/data
    //------------------------------------

    TSpectrum pfinder;
    for (int bin=0;bin<bins::smearing_bin_size;bin++) {

        int rebin_factor = RebinHistogram(hist_z_metl[bin],0);
        RebinHistogram(hist_g_metl[bin],rebin_factor);
        RebinHistogram(hist_g_metl_smear[bin],rebin_factor);
        RebinHistogram(hist_g_metl_smear_2j[bin],rebin_factor);
        for (int bin1=0; bin1<bins::METl_bin_size; bin1++) {
            RebinHistogram(hist_z_mll_pt[bin][bin1], 0);
        }
        int newbin = 40000/rebin_factor;

        Double_t *z_smear_in = new Double_t[newbin];
        Double_t g_smear_in[newbin];
        Double_t j_resp_in[newbin];
        Double_t *z_resp_in = new Double_t[newbin];
        Double_t *g_resp_in = new Double_t[newbin];
        smear_raw[bin] = new TH1D(TString("smear_raw_")+TString::Itoa(bin,10),"",newbin,-30000,10000);
        for (int i=0;i<newbin;i++) {
            z_smear_in[i] = max(hist_z_metl[bin]->GetBinContent(i+1),0.);
            if (i<newbin/2) g_smear_in[i] = max(hist_g_metl[bin]->GetBinContent(i+1+newbin/2),0.);
            else g_smear_in[i] = 0.;
            z_resp_in[i] = max(hist_z_metl[bin]->GetBinContent(i+1),0.);
            g_resp_in[i] = max(hist_g_metl[bin]->GetBinContent(i+1),0.);
            if (i<newbin/2) j_resp_in[i] = max(hist_z_onshell_metl[bin]->GetBinContent(i+1+newbin/2),0.);
            else j_resp_in[i] = 0.;
        }
        pfinder.Deconvolution(z_smear_in,g_smear_in,newbin,1000,1,1.0);
        pfinder.Deconvolution(z_resp_in,j_resp_in,newbin,1000,1,1.0);
        pfinder.Deconvolution(g_resp_in,j_resp_in,newbin,1000,1,1.0);
        for (int i=0;i<newbin;i++) {
            smear_raw[bin]->SetBinContent(i+1,z_smear_in[i]);
        }
        float smear_mean = smear_raw[bin]->GetMean();
        float smear_rms = smear_raw[bin]->GetRMS();

        float hist_gmetl_rms = hist_g_metl[bin]->GetRMS();
        float hist_zmetl_rms = hist_z_metl[bin]->GetRMS();
        for (int i=0;i<newbin;i++) {
            if (hist_gmetl_rms/hist_zmetl_rms > 1.0) {
                smear_raw[bin]->SetBinContent(i+1,0.);
            }
            float smear_cut = 6.;
            if (channel=="mm" && bins::pt_bins[bin]>=0) smear_cut = 7.;
            if (abs(smear_raw[bin]->GetBinCenter(i+1)-smear_mean)/smear_rms>smear_cut) {
                smear_raw[bin]->SetBinContent(i+1,0.);
            }
        }

        shift[bin] = -hist_g_metl[bin]->GetMean();
    }

    for (int bin=0;bin<bins::smearing_bin_size;bin++) {
        smear_final[bin] = new TH1D(TString("smear_final_")+TString::Itoa(bin,10),"",500,-1000,1000);
        for (int i=0;i<500;i++) {
            int which_bin = smear_raw[bin]->FindBin(smear_final[bin]->GetBinCenter(i+1));
            smear_final[bin]->SetBinContent(i+1,smear_raw[bin]->GetBinContent(which_bin));
        }
    }
}

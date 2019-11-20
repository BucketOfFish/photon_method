#include "../Common/Settings.C"
#include "../Common/CommonLibraries.C"
#include "../Common/CommonFunctions.C"
#include "MT2.h"

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

TH1D* hist_METl_bins = new TH1D("hist_METl_bins", "", bins::METl_bin_size, bins::METl_bins);
TH1D* hist_pt_bins = new TH1D("hist_pt_bins", "", bins::smearing_bin_size, bins::pt_bins);
TH1D* hist_MET_bins = new TH1D("hist_MET_bins", "", bins::smearing_bin_size, bins::MET_bins); //hist_MET_bins->SetStats(0);

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

    TH1D* hist_pt_bins = new TH1D("pt_bins","",bins::smearing_bin_size,bins::pt_bins);

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
            if (jet_n<2 || lep_pT->at(0)<cuts::leading_lep_pt_cut || lep_pT->at(1)<cuts::second_lep_pt_cut) continue;
            int pt_bin = hist_pt_bins->FindBin(ptll)-1;
            int METl_bin = hist_METl_bins->FindBin(METl)-1;
            std::cout << pt_bin << std::endl;
            hist_z_metl[pt_bin];
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
        TFile* photon_file = new TFile(ntuple_path + "/g_mc/" + mc_period + "_SinglePhoton222.root");
    else if (smearing_method == 5)
        TFile* photon_file = new TFile(ntuple_path + "/g_data/" + data_period + "_photon.root");

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
        int pt_bin = hist_pt_bins->FindBin(ptll)-1;
        hist_g_metl[pt_bin]->Fill(METl, totalWeight);
    }

    photon_file->Close();

}

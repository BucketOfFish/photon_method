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
TH1D* g_metl[bins::smearing_bin_size];
TH1D* z_onshell_metl[bins::smearing_bin_size];

void GetSmearingHistogram(string target_channel, float lumi, string period, int smearing_method) {

    cout << "GetSmearingHistogram : smearing_method " << smearing_method << endl;

    for (int bin=0;bin<bins::smearing_bin_size;bin++) {
        z_metl[bin] = new TH1D(TString("z_metl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        z_onshell_metl[bin] = new TH1D(TString("z_jetmetl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        g_metl[bin] = new TH1D(TString("g_metl_")+TString::Itoa(bin,10),"",40000,-30000,10000);
    }

    TH1D* pt_bins = new TH1D("pt_bins","",bins::smearing_bin_size,bins::pt_bins);

    //------------------------------------
    // SMEARING METHOD 0: no smearing
    //------------------------------------

    if (smearing_method == 0) return;

    //------------------------------------
    // SMEARING METHOD 4/5: for R21 MC/data
    //------------------------------------

    std::cout << "Getting Z METl binned by pT." << std::endl;

    TFile* data_file = new TFile((ntuple_path + "/bkg_data/" + period + "_bkg.root").c_str());
    TFile* ttbar_mc_file = new TFile((ntuple_path + "/bkg_mc/" + period + "_ttbar.root").c_str());
    TFile* diboson_mc_file = new TFile((ntuple_path + "/bkg_mc/" + period + "_diboson.root").c_str());
    TFile* Z_mc_file = new TFile((ntuple_path + "/bkg_mc/" + period + "_Zjets.root").c_str());

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
        int weight = file_weights[i];

        tree->SetBranchStatus("*", 0);
        double totalWeight; SetInputBranch(tree, "totalWeight", &totalWeight);
        int jet_n; SetInputBranch(tree, "nJet30", &jet_n);
        int bjet_n; SetInputBranch(tree, "bjet_n", &bjet_n);
        float ptll; SetInputBranch(tree, "Ptll", &ptll);
        float mll; SetInputBranch(tree, "mll", &mll);
        float METl; SetInputBranch(tree, "METl", &METl);
        int channel; SetInputBranch(tree, "channel", &channel);

        for (int entry=0; entry<tree->GetEntries(); entry++) {
            tree->GetEntry(entry);
            if (TString(target_channel).EqualTo("ee") && channel != 1) continue;
            if (TString(target_channel).EqualTo("mm") && channel != 0) continue;
            if (ptll<50. || jet_n!=1 || bjet_n!=0) continue;
            int pt_bin = pt_bins->FindBin(ptll)-1;
            z_metl[pt_bin]->Fill(METl, weight*totalWeight);
            if (mll<90 || mll>92) continue;
            z_onshell_metl[pt_bin]->Fill(METl, weight*totalWeight);
        }
    }

    data_file->Close();
    ttbar_mc_file->Close();
    diboson_mc_file->Close();
    Z_mc_file->Close();

    //-- photon histogram

    std::cout << "Getting photon METl binned by pT." << std::endl;

    TFile* photon_file;
    if (smearing_method == 4)
        TFile* photon_file = new TFile((ntuple_path + "/bkg_mc/" + period + "_SinglePhoton222.root").c_str());
    else if (smearing_method == 5)
        TFile* photon_file = new TFile((ntuple_path + "/bkg_data/" + period + "_photon.root").c_str());

    TTree* tree = (TTree*)photon_file->Get("BaselineTree");
    tree->SetBranchStatus("*", 0);
    double totalWeight; SetInputBranch(tree, "totalWeight", &totalWeight);
    int jet_n; SetInputBranch(tree, "nJet30", &jet_n);
    int bjet_n; SetInputBranch(tree, "bjet_n", &bjet_n);
    float ptll; SetInputBranch(tree, "Ptll", &ptll);
    float METl; SetInputBranch(tree, "METl", &METl);

    for (int entry=0; entry<tree->GetEntries(); entry++) {
        tree->GetEntry(entry);
        if (ptll<50. || jet_n!=1 || bjet_n!=0) continue;
        int pt_bin = pt_bins->FindBin(ptll)-1;
        g_metl[pt]->Fill(METl, totalWeight);
    }

    photon_file->Close();

}

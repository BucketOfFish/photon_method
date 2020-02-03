#include "../Common/Settings.C"
#include "MT2.h"

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

//TH1D* hist_z_mll_bin_pt_metl[bins::n_pt_bins+2][bins::n_METl_bins+2];

vector<TH1D> GetSmearingDistribution(string channel, TString period, string data_or_mc, bool diagnostics) {
//map<int, pair<float, float>> GetSmearingDistribution(string channel, TString period, string data_or_mc, bool diagnostics) {

    /**
     * Returns map with key = photon pt bin, value = (mean, std).
     * For a photon in a given pt bin, smear the event's MET using the given Gaussian numbers.
     */

    //------------------------------------
    // FILL METL HISTOGRAMS
    //------------------------------------

    TH1D* hist_z_metl_bin_pt[bins::n_pt_bins+2];  // 0 = underflow, n_pt_bins + 1 = overflow
    TH1D* hist_g_metl_bin_pt[bins::n_pt_bins+2];

    auto FillHistograms = [&hist_z_metl_bin_pt, &hist_g_metl_bin_pt] (string target_channel, TString period, string data_or_mc) {

        /**
         * Fills histogram arrays hist_z_metl_bin_pt, hist_g_metl_bin_pt, and hist_z_mll_bin_pt_metl (2D)
         * with Z/photon METl/mll distributions in binned regions.
         */

        cout << "Filling smearing histograms" << endl;

        for (int bin=0; bin<bins::n_pt_bins+2; bin++) {
            hist_z_metl_bin_pt[bin] = new TH1D(TString("hist_z_metl_")+TString::Itoa(bin,10),"",bins::n_smearing_bins,bins::smearing_low,bins::smearing_high);
            hist_g_metl_bin_pt[bin] = new TH1D(TString("hist_g_metl_")+TString::Itoa(bin,10),"",bins::n_smearing_bins,bins::smearing_low,bins::smearing_high);
            for (int bin1=0; bin1<bins::n_METl_bins+2; bin1++) {
                //hist_z_mll_bin_pt_metl[bin][bin1] = new TH1D(TString("hist_z_Mll_dPt_")+TString::Itoa(bin,10)+TString("_")+TString::Itoa(bin1,10),"",bins::n_mll_bins,bins::mll_bin);
            }
        }

        TString mc_period = MCPeriod(period);
        TString data_period = DataPeriod(period);

        cout << "Getting Z histograms binned by pT and METl." << endl;

        TFile* data_file = new TFile(ntuple_path + "/bkg_data/" + data_period + "_bkg.root");
        TFile* ttbar_mc_file = new TFile(ntuple_path + "/bkg_mc/" + mc_period + "_ttbar.root");
        TFile* diboson_mc_file = new TFile(ntuple_path + "/bkg_mc/" + mc_period + "_diboson.root");
        TFile* Z_mc_file = new TFile(ntuple_path + "/bkg_mc/" + mc_period + "_Zjets.root");

        vector<TFile*> contributing_files;
        vector<int> file_weights;

        if (data_or_mc == "MC") {
            contributing_files = {Z_mc_file};
            file_weights = {1};
        }
        else if (data_or_mc == "Data") {
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
            int nLep_signal; SetInputBranch(tree, "nLep_signal", &nLep_signal);
            vector<float>* lep_pT = new vector<float>(10); SetInputBranch(tree, "lepPt", &lep_pT);

            //tree->Draw(">>event_list", cuts::reweight_region);
            //tree->Draw(">>event_list", "nJet30>=2 && lepPt[0]>25 && lepPt[1]>25");
            tree->Draw(">>event_list", "nJet30>=2");
            auto event_list = (TEventList*) gDirectory->Get("event_list");
            for (int entry=0; entry<event_list->GetN(); entry++) {
                tree->GetEntry(event_list->GetEntry(entry));
                if (TString(target_channel).EqualTo("ee") && channel != 1) continue;
                if (TString(target_channel).EqualTo("mm") && channel != 0) continue;
                int pt_bin = bins::hist_pt_bins->FindBin(ptll);
                int METl_bin = bins::hist_METl_bins->FindBin(METl);
                //hist_z_mll_bin_pt_metl[pt_bin][METl_bin]->Fill(mll, fileWeight*totalWeight);
                //if (ptll<50. || jet_n<2 || lep_pT->at(0)<cuts::leading_lep_pt_cut || lep_pT->at(1)<cuts::second_lep_pt_cut) continue;
                if (jet_n<2 || lep_pT->at(0)<cuts::leading_lep_pt_cut || lep_pT->at(1)<cuts::second_lep_pt_cut) continue;
                hist_z_metl_bin_pt[pt_bin]->Fill(METl, fileWeight*totalWeight);
            }
        }

        data_file->Close();
        ttbar_mc_file->Close();
        diboson_mc_file->Close();
        Z_mc_file->Close();

        //-- photon histogram

        cout << "Getting photon histogram binned by pT." << endl;

        TFile* photon_file;
        if (data_or_mc == "MC")
            photon_file = new TFile(ntuple_path + "/g_mc/" + mc_period + "_SinglePhoton222.root");
        else if (data_or_mc == "Data")
            photon_file = new TFile(ntuple_path + "/g_data/" + data_period + "_photon.root");

        TTree* tree = (TTree*)photon_file->Get("BaselineTree");
        tree->SetBranchStatus("*", 0);
        double totalWeight; SetInputBranch(tree, "totalWeight", &totalWeight);
        int jet_n; SetInputBranch(tree, "nJet30", &jet_n);
        int bjet_n; SetInputBranch(tree, "bjet_n", &bjet_n);
        float ptll; SetInputBranch(tree, "gamma_pt", &ptll);
        int nLep_signal; SetInputBranch(tree, "nLep_signal", &nLep_signal);
        float METl; SetInputBranch(tree, "METl_raw", &METl);

        //tree->Draw(">>event_list", cuts::reweight_region);
        tree->Draw(">>event_list", "nJet30>=2");
        auto event_list = (TEventList*) gDirectory->Get("event_list");
        for (int entry=0; entry<event_list->GetN(); entry++) {
            tree->GetEntry(event_list->GetEntry(entry));
            //if (ptll<50. || jet_n!=1 || bjet_n!=0) continue;
            int pt_bin = bins::hist_pt_bins->FindBin(ptll);
            hist_g_metl_bin_pt[pt_bin]->Fill(METl, totalWeight);
        }

        photon_file->Close();
    };

    FillHistograms(channel, period, data_or_mc);

    //------------------------------------
    // REBINNING FUNCTION
    //------------------------------------

    auto RebinHistogram = [] (TH1D* hist, int rebin) {

        /**
         * If rebin=0, rebin histogram by factors of 2, up to 32, until it looks reasonable or there are an odd number of bins.
         * This means having abs(negative_yield/positive_yield)<0.005 and core_yield/positive_yield>0.4.
         * Else, rebin by #(rebin).
         * Finally, set all negative bins values in the histogram to 0.
         * Return the final rebinning factor for the histogram.
         */

        float negative_yield = 0.; // sum of bins with negative values
        float positive_yield = 0.001; // sum of bins with positive values
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
    };

    //------------------------------------
    // PERFORM SMEARING
    //------------------------------------

    map<int, pair<float, float>> smearing_gaussians;
    TSpectrum deconv_tool;
    vector<TH1D> metl_deconv_hists;
    for (int pt_bin=0;pt_bin<bins::n_pt_bins+2;pt_bin++) {

        //--- Rebin METl histograms to look reasonable
        int rebin_factor = RebinHistogram(hist_z_metl_bin_pt[pt_bin],0);
        RebinHistogram(hist_g_metl_bin_pt[pt_bin],rebin_factor);
        //for (int met_bin=0; met_bin<bins::n_METl_bins+2; met_bin++)
            //RebinHistogram(hist_z_mll_bin_pt_metl[pt_bin][met_bin], 0);
        int n_hist_bins = bins::n_smearing_bins/rebin_factor;

        //--- Save smearing value
        float smear_mean = hist_z_metl_bin_pt[pt_bin]->GetMean() - hist_g_metl_bin_pt[pt_bin]->GetMean();
        float smear_rms = sqrt(pow(hist_z_metl_bin_pt[pt_bin]->GetRMS(), 2) - pow(hist_g_metl_bin_pt[pt_bin]->GetRMS(), 2));
        if (isnan(smear_rms)) smear_rms = 0.0;
        
        pair<float, float> gaussian = make_pair(0.0, 0.0);
        if (channel=="ee") gaussian = make_pair(smear_mean, 0.0);
        else gaussian = make_pair(smear_mean, smear_rms);
        smearing_gaussians.insert(pair<int, pair<float, float>>(pt_bin, gaussian));

        //cout << "Pt bin " << pt_bin << endl;
        //cout << "Z_METl " << ": [" << endl;
        //for (int i=0;i<n_hist_bins+2;i++) {
            //cout << hist_z_metl_bin_pt[pt_bin]->GetBinContent(i) << ", ";
        //}
        //cout << "]" << endl;
        //cout << "g_METl " << ": [" << endl;
        //for (int i=0;i<n_hist_bins+2;i++) {
            //cout << hist_g_metl_bin_pt[pt_bin]->GetBinContent(i) << ", ";
        //}
        //cout << "]" << endl;

        //--- Deconvolution
        Double_t *z_metl_dist = new Double_t[n_hist_bins]; // pointer for passing to TSpectrum
        Double_t g_metl_dist[n_hist_bins];
        for (int i=0;i<n_hist_bins;i++) {
            z_metl_dist[i] = max(hist_z_metl_bin_pt[pt_bin]->GetBinContent(i+1),0.);
            g_metl_dist[i] = max(hist_g_metl_bin_pt[pt_bin]->GetBinContent(i+1),0.);
        }
        deconv_tool.Deconvolution(z_metl_dist,g_metl_dist,n_hist_bins,1000,1,1.0);
        TH1D new_metl_deconv_hist = TH1D(TString("metl_deconv_")+TString::Itoa(pt_bin, 10),"",n_hist_bins,-10000,10000);
        for (int i=1;i<=n_hist_bins;i++) {
            new_metl_deconv_hist.SetBinContent(i, z_metl_dist[i-1]);
        }
        metl_deconv_hists.push_back(new_metl_deconv_hist);

        //-- Print diagnostics and save histograms
        if (!diagnostics) continue;

        cout << "Pt bin " << pt_bin << " has a smearing mean, std of " << gaussian.first << ", " << gaussian.second << endl;

        TCanvas *canvas = new TCanvas("canvas","canvas",600,600);
        canvas->cd();
        canvas->SetLogy();

        TH1D* z_hist = hist_z_metl_bin_pt[pt_bin];
        TString z_plot_name = Form("Plots/z_ptbin_%d.eps", pt_bin);
        z_hist->SetLineColor(1); z_hist->SetFillColor(42); z_hist->SetLineStyle(1);
        z_hist->GetXaxis()->SetTitle("METl");
        z_hist->GetYaxis()->SetTitle("entries / bin");
        z_hist->Draw("hist");
        canvas->Print(z_plot_name);

        TH1D* g_hist = hist_g_metl_bin_pt[pt_bin];
        TString g_plot_name = Form("Plots/g_ptbin_%d.eps", pt_bin);
        g_hist->SetLineColor(1); g_hist->SetFillColor(42); g_hist->SetLineStyle(1);
        g_hist->GetXaxis()->SetTitle("METl");
        g_hist->GetYaxis()->SetTitle("entries / bin");
        g_hist->Draw("hist");
        canvas->Print(g_plot_name);

        //for (int metl_bin=0; metl_bin<bins::n_METl_bins+2; metl_bin++) {
            //TH1D* z_mll_hist = hist_z_mll_bin_pt_metl[pt_bin][metl_bin];
            //TString z_mll_plot_name = Form("Plots/z_mll_ptbin_%d_%d.eps", pt_bin, metl_bin);
            //z_mll_hist->SetLineColor(1); z_mll_hist->SetFillColor(42); z_mll_hist->SetLineStyle(1);
            //z_mll_hist->GetXaxis()->SetTitle("mll");
            //z_mll_hist->GetYaxis()->SetTitle("entries / bin");
            //z_mll_hist->Draw("hist");
            //canvas->Print(z_mll_plot_name);
        //}

        delete canvas;
    }

    //return smearing_gaussians;
    return metl_deconv_hists;
}

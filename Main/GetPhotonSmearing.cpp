#include "Settings.cpp"
#include "MT2.cpp"

using namespace std;

TH1F* GetLepThetaHistogram(SmearingOptions options, string period, string channel, string data_or_mc) {

    cout << "Getting Z lepton CM theta distribution histogram." << endl;
    gStyle->SetOptStat(0);

    //--- open files and create TChains
    string data_filename = options.in_file_path + period + "_data_bkg.root";
    if (period == "data15-16") period = "mc16a";
    else if (period == "data17") period = "mc16cd";
    else if (period == "data18") period = "mc16e";
    string tt_filename = options.in_file_path + period + "_ttbar.root";
    string vv_filename = options.in_file_path + period + "_diboson.root";
    string zjets_filename = options.in_file_path + period + "_Zjets.root";

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

vector<vector<TH1D*>> hist_z_mll_bin_pt_metl;

//vector<TH1D> GetSmearingDistribution(string channel, TString period, string data_or_mc, bool diagnostics) {
map<int, pair<float, float>> GetSmearingDistribution(SmearingOptions options, string channel, TString period, string data_or_mc, bool diagnostics) {

    /**
     * Returns map with key = photon pt bin, value = (mean, std).
     * For a photon in a given pt bin, smear the event's MET using the given Gaussian numbers.
     */

    //------------------------------------
    // FILL METL HISTOGRAMS
    //------------------------------------

    TH1D* hist_z_metl_bin_pt[bins::n_pt_bins+2];  // 0 = underflow, n_pt_bins + 1 = overflow
    TH1D* hist_g_metl_bin_pt[bins::n_pt_bins+2];

    auto FillHistograms = [&options, &hist_z_metl_bin_pt, &hist_g_metl_bin_pt] (string target_channel, TString period, string data_or_mc) {

        /**
         * Fills histogram arrays hist_z_metl_bin_pt, hist_g_metl_bin_pt, and hist_z_mll_bin_pt_metl (2D)
         * with Z/photon METl/mll distributions in binned regions.
         */

        cout << "Filling smearing histograms" << endl;

        for (int bin=0; bin<bins::n_pt_bins+2; bin++) {
            hist_z_metl_bin_pt[bin] = new TH1D(TString("hist_z_metl_")+TString::Itoa(bin,10),"",bins::n_smearing_bins,bins::smearing_low,bins::smearing_high);
            hist_g_metl_bin_pt[bin] = new TH1D(TString("hist_g_metl_")+TString::Itoa(bin,10),"",bins::n_smearing_bins,bins::smearing_low,bins::smearing_high);

            vector<TH1D*> hist_z_mll_bin_pt_metl_pt_bin;
            for (int bin1=0; bin1<bins::n_METl_bins+2; bin1++) {
                hist_z_mll_bin_pt_metl_pt_bin.push_back(new TH1D(TString("hist_z_Mll_dPt_")+TString::Itoa(bin,10)+TString("_")+TString::Itoa(bin1,10),"",bins::n_mll_bins,bins::mll_bin));
            }
            hist_z_mll_bin_pt_metl.push_back(hist_z_mll_bin_pt_metl_pt_bin);
        }

        string mc_period = getMCPeriod(period);
        string data_period = DataPeriod(period);

        cout << "Getting Z histograms binned by pT and METl." << endl;

        TFile* data_file = new TFile((options.in_file_path + data_period + "_data_bkg.root").c_str());
        TFile* ttbar_mc_file = new TFile((options.in_file_path + mc_period + "_ttbar.root").c_str());
        TFile* diboson_mc_file = new TFile((options.in_file_path + mc_period + "_diboson.root").c_str());
        TFile* Z_mc_file = new TFile((options.in_file_path + mc_period + "_Zjets.root").c_str());

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
            tree->Draw(">>event_list", "nJet30>=2 && lepPt[0]>25 && lepPt[1]>25");
            //tree->Draw(">>event_list", "nJet30>=2");
            auto event_list = (TEventList*) gDirectory->Get("event_list");
            for (int entry=0; entry<event_list->GetN(); entry++) {
                tree->GetEntry(event_list->GetEntry(entry));
                if (TString(target_channel).EqualTo("ee") && channel != 1) continue;
                if (TString(target_channel).EqualTo("mm") && channel != 0) continue;
                int pt_bin = bins::hist_pt_bins->FindBin(ptll);
                int METl_bin = bins::hist_METl_bins->FindBin(METl);
                hist_z_mll_bin_pt_metl[pt_bin][METl_bin]->Fill(mll, fileWeight*totalWeight);
                //if (ptll<50. || jet_n<2 || lep_pT->at(0)<cuts::leading_lep_pt_cut || lep_pT->at(1)<cuts::second_lep_pt_cut) continue;
                //if (jet_n<2 || lep_pT->at(0)<cuts::leading_lep_pt_cut || lep_pT->at(1)<cuts::second_lep_pt_cut) continue;
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
            photon_file = new TFile((options.in_file_path + mc_period + "_SinglePhoton222.root").c_str());
        else if (data_or_mc == "Data")
            photon_file = new TFile((options.in_file_path + data_period + "_data_photon.root").c_str());

        TTree* tree = (TTree*)photon_file->Get("BaselineTree");
        tree->SetBranchStatus("*", 0);
        double totalWeight; SetInputBranch(tree, "totalWeight", &totalWeight);
        int jet_n; SetInputBranch(tree, "nJet30", &jet_n);
        int bjet_n; SetInputBranch(tree, "bjet_n", &bjet_n);
        float ptll; SetInputBranch(tree, "gamma_pt", &ptll);
        int nLep_signal; SetInputBranch(tree, "nLep_signal", &nLep_signal);
        float METl; SetInputBranch(tree, "METl_unsmeared", &METl);

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
    //TSpectrum deconv_tool;
    //vector<TH1D> metl_deconv_hists;
    for (int pt_bin=0;pt_bin<bins::n_pt_bins+2;pt_bin++) {

        //--- Rebin METl histograms to look reasonable
        int rebin_factor = RebinHistogram(hist_z_metl_bin_pt[pt_bin],0);
        RebinHistogram(hist_g_metl_bin_pt[pt_bin],rebin_factor);
        for (int met_bin=0; met_bin<bins::n_METl_bins+2; met_bin++)
            RebinHistogram(hist_z_mll_bin_pt_metl[pt_bin][met_bin], 0);
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
        //Double_t *z_metl_dist = new Double_t[n_hist_bins]; // pointer for passing to TSpectrum
        //Double_t g_metl_dist[n_hist_bins];
        //for (int i=0;i<n_hist_bins;i++) {
            //z_metl_dist[i] = max(hist_z_metl_bin_pt[pt_bin]->GetBinContent(i+1),0.);
            //g_metl_dist[i] = max(hist_g_metl_bin_pt[pt_bin]->GetBinContent(i+1),0.);
        //}
        //deconv_tool.Deconvolution(z_metl_dist,g_metl_dist,n_hist_bins,1000,1,1.0);
        //TH1D new_metl_deconv_hist = TH1D(TString("metl_deconv_")+TString::Itoa(pt_bin, 10),"",n_hist_bins,-10000,10000);
        //for (int i=1;i<=n_hist_bins;i++) {
            //new_metl_deconv_hist.SetBinContent(i, z_metl_dist[i-1]);
        //}
        //metl_deconv_hists.push_back(new_metl_deconv_hist);

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

        for (int metl_bin=0; metl_bin<bins::n_METl_bins+2; metl_bin++) {
            TH1D* z_mll_hist = hist_z_mll_bin_pt_metl[pt_bin][metl_bin];
            TString z_mll_plot_name = Form("Plots/z_mll_ptbin_%d_%d.eps", pt_bin, metl_bin);
            z_mll_hist->SetLineColor(1); z_mll_hist->SetFillColor(42); z_mll_hist->SetLineStyle(1);
            z_mll_hist->GetXaxis()->SetTitle("mll");
            z_mll_hist->GetYaxis()->SetTitle("entries / bin");
            z_mll_hist->Draw("hist");
            canvas->Print(z_mll_plot_name);
        }

        delete canvas;
    }

    return smearing_gaussians;
    //return metl_deconv_hists;
}

void GetPhotonSmearing(SmearingOptions options, string period, string channel, string data_or_mc, bool turn_off_shifting_and_smearing=false) {

    //---------------------------------------------
    // get unsmeared input file and create smeared output file
    //---------------------------------------------

    string data_period = DataPeriod(period);
    string mc_period = getMCPeriod(period);

    cout << "channel         " << channel         << endl;
    cout << "period          " << period          << endl;
    cout << "isData?         " << data_or_mc          << endl;
    cout << "smearing output " << options.out_file_name   << endl;

    TH1::SetDefaultSumw2();

    string infilename;
    if (data_or_mc == "MC") infilename = options.in_file_path + mc_period + "_SinglePhoton222.root";
    else if (data_or_mc == "Data") infilename = options.in_file_path + data_period + "_data_photon.root";

    TChain* inputTree = new TChain("BaselineTree");
    inputTree->Add(infilename.c_str());

    cout << endl;
    cout << "Opening file           : " << infilename        << endl;
    cout << "Events in ntuple       : " << inputTree->GetEntries() << endl;

    TH1::SetDefaultSumw2();

    string outfilename = options.out_file_name;

    TFile* outputFile = new TFile(outfilename.c_str(), "recreate");          
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
    float METl; CopyBranch(inputTree, BaselineTree, "METl_unsmeared", "METl_unsmeared", &METl, "F");
    float METt; CopyBranch(inputTree, BaselineTree, "METt_unsmeared", "METt_unsmeared", &METt, "F");
    float METl_smeared; BaselineTree->Branch("METl", &METl_smeared, "METl/F");
    float METt_smeared; BaselineTree->Branch("METt", &METt_smeared, "METt/F");
    float HT; CopyBranch(inputTree, BaselineTree, "Ht30", "Ht30", &HT, "F");
    float MET_raw; CopyBranch(inputTree, BaselineTree, "met_Et_unsmeared", "met_Et_unsmeared", &MET_raw, "F");
    float MET_phi; CopyBranch(inputTree, BaselineTree, "met_Phi", "met_Phi", &MET_phi, "F");

    vector<float>* jet_pT = new vector<float>(10); CopyBranch(inputTree, BaselineTree, "jetPt", "jetPt", &jet_pT, "vector<float>");
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
    float lep_theta_cm; BaselineTree->Branch("Z_cm_lep_theta", &lep_theta_cm, "Z_cm_lep_theta/F");

    //--- HistFitter branches
    vector<string> histFitterBranches {"DatasetNumber/I", "H2PP/D", "H5PP/D", "H5PP_VR/D",
        "METOverPtISR/F", "METOverPtW/F", "METOverPtZ/F", "MJ/D", "MJ_VR/D", "MZ/D", "MZ_VR/D", "NjISR/D",
        "NjS/D", "PTCM/D", "PTCM_VR/D", "PTI/D", "PTISR/D", "PTISR_VR/D", "PTI_VR/D", "RISR/D", "RISR_VR/D",
        "RPT_HT5PP/D", "RPT_HT5PP_VR/D", "R_minH2P_minH3P/D", "R_minH2P_minH3P_VR/D", "Rjj/F", "Rll/F",
        "dPhiMetISR/F", "dPhiMetJet1/F", "dPhiMetJet2/F", "dPhiMetJet12Min/F", "dPhiPjjMet/F", "dPhiPllMet/F",
        "dphiISRI/D", "dphiISRI_VR/D", "dphiVP/D", "dphiVP_VR/D", "lept1Pt_VR/D", "lept2Pt_VR/D", "mTl3/D",
        "MET_sig/F", "minDphi/D", "minDPhi2JetsMet/F", "mll_RJ/D", "mll_RJ_VR/D", "nJet30/I", "nJet20/I", "mjj/F",
        "nBJet20_MV2c10_FixedCutBEff_77/I", "trigMatch_2LTrigOR/I", "genWeight/D", "eventWeight/D", "leptonWeight/D", "jvtWeight/D", "bTagWeight/D", "pileupWeight/D", "globalDiLepTrigSF/D", "RunNumber/I", "RandomRunNumber/I", "trigMatch_2LTrig/I", "lumi/D"};
    CopyAllBranches(inputTree, BaselineTree, histFitterBranches);

    vector<float>* dPhiMetJet = new vector<float>(10); CopyBranch(inputTree, BaselineTree, "dPhiMetJet", "dPhiMetJet", &dPhiMetJet, "vector<float>");
    vector<float>* jetM = new vector<float>(10); CopyBranch(inputTree, BaselineTree, "jetM", "jetM", &jetM, "vector<float>");
    float mll; BaselineTree->Branch("mll", &mll, "mll/F");
    int is_OS = 1; BaselineTree->Branch("is_OS", &is_OS, "is_OS/I");
    vector<float>* lep_pT = new vector<float>(10); BaselineTree->Branch("lepPt", "vector<float>", &lep_pT);
    vector<float>* lep_eta = new vector<float>(10); BaselineTree->Branch("lepEta", "vector<float>", &lep_eta);
    vector<float>* lep_phi = new vector<float>(10); BaselineTree->Branch("lepPhi", "vector<float>", &lep_phi);
    vector<int>* lep_flavor = new vector<int>(10); BaselineTree->Branch("lepFlavor", "vector<int>", &lep_flavor);
    vector<int>* lep_charge = new vector<int>(10); BaselineTree->Branch("lepCharge", "vector<int>", &lep_charge);
    vector<float>* lep_m = new vector<float>(10); BaselineTree->Branch("lepM", "vector<float>", &lep_m);
    vector<int>* lepIsoFCTight = new vector<int>{1,1}; BaselineTree->Branch("lepIsoFCTight", "vector<int>", &lepIsoFCTight);
    vector<int>* lepIsPR = new vector<int>{1,1}; BaselineTree->Branch("lepIsPR", "vector<int>", &lepIsPR);
    Int_t lepChannel; BaselineTree->Branch("channel", &lepChannel, "channel/I");

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
    map<int, pair<float, float>> smearing_gaussians = GetSmearingDistribution(options, channel, period, data_or_mc, diagnostics);
    //vector<TH1D> smearing_hists = GetSmearingDistribution(options, channel, period, data_or_mc, diagnostics);

    //-----------------------------
    // get Z lepton CM theta distribution
    //-----------------------------

    TH1F* h_lep_theta = GetLepThetaHistogram(options, period, channel, data_or_mc);

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
    nentries = 100;
    for (Long64_t i=0; i<nentries; i++) {

        if (fmod(i,1e5)==0) cout << i << " events processed.\r";
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
            lep_m->clear();
            if (lepChannel == 0) {
                lep_m->push_back(0.1056583);
                lep_m->push_back(0.1056583);
            }
            else if (lepChannel == 1) {
                lep_m->push_back(0.0005109);
                lep_m->push_back(0.0005109);
            }

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
    cout << endl;

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

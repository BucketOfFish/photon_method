#include "Settings.cpp"
#include "MT2.cpp"

using namespace std;

//------------------
// HELPER FUNCTIONS
//------------------

void initFeatureHists(map<string, map<int, TH1D*>> &plots) {
    //options.diagnostic_plots = vector<string>{"lepEta", "METl", "mll"};
    //for (auto feature : options.diagnostic_plots) {
    for (int i=0; i<bins::n_pt_bins+2; i++) {
        plots["lep_cm_theta"][i] = new TH1D("", "", 100, 0, 3);
        plots["lepEta"][i] = new TH1D("", "", 100, -3, 3);
        plots["lepPt"][i] = new TH1D("", "", 100, 0, 300);
        plots["met_Et_raw"][i] = new TH1D("", "", 100, 0, 300);
        plots["met_Et"][i] = new TH1D("", "", 100, 0, 300);
        plots["METl_raw"][i] = new TH1D("", "", 100, -100, 100);
        plots["METl"][i] = new TH1D("", "", 100, -150, 150);
        plots["METt"][i] = new TH1D("", "", 100, -150, 150);
        plots["mll"][i] = new TH1D("", "", 100, 0, 150);
    }
//------------------
}

//-----------------
// CONVERTOR CLASS
//-----------------

class PhotonToZConverter {
public:
    Options options;

    string in_file_name;
    TChain *inputTree;

    string out_file_name;
    TFile* outputFile;
    TTree *outputTree;

    std::default_random_engine *random_generator;
    TRandom1 myRandom;

    TH1D* hist_z_metl_bin_pt[bins::n_pt_bins+2];  // 0 = underflow, n_pt_bins + 1 = overflow
    TH1D* hist_g_metl_bin_pt[bins::n_pt_bins+2];
    TH1D* hist_g_smeared_metl_bin_pt[bins::n_pt_bins+2];
    vector<vector<TH1D*>> hist_z_mll_bin_pt_metl;

    TH1D* h_lep_cm_theta;
    std::vector<float> lep_cm_theta_bin_bounds;
    std::discrete_distribution<int> *lep_cm_theta_distribution;
    map<int, normal_distribution<float>> smearing_gaussians;

    //----------------
    // INITIALIZATION
    //----------------

    void openPhotonFile(Options options) {
        TH1::SetDefaultSumw2();

        if (options.run_vgamma) {
            this->in_file_name = options.reduction_folder + options.mc_period + "_Vgamma.root";
            this->out_file_name = options.smearing_folder + options.mc_period + "_Vgamma_" + options.channel + ".root";
        }
        else {
            if (options.is_data) {
                this->in_file_name = options.reduction_folder + options.data_period + "_data_photon.root";
                this->out_file_name = options.smearing_folder + options.data_period + "_data_photon_" + options.channel + ".root"; 
            }
            else {
                this->in_file_name = options.reduction_folder + options.mc_period + "_SinglePhoton222.root";
                this->out_file_name = options.smearing_folder + options.mc_period + "_SinglePhoton222_" + options.channel + ".root";
            }
        }

        cout << padString("channel") << options.channel         << endl;
        cout << padString("period") << options.period          << endl;
        cout << padString("data-based photons?") << options.is_data          << endl;
        cout << padString("smearing output") << this->out_file_name   << endl;

        // run on MC Vgamma if option is set
        if (options.run_vgamma) this->in_file_name = options.reduction_folder + options.mc_period + "_Vgamma.root";

        this->inputTree = new TChain("BaselineTree");
        this->inputTree->Add(this->in_file_name.c_str());

        cout << endl;
        cout << padString("Opening read file") << this->in_file_name << endl;
        cout << padString("Events in ntuple") << inputTree->GetEntries() << endl;
    }

    void openOutputFile(Options options) {
        this->outputFile = new TFile(this->out_file_name.c_str(), "recreate");          
        this->outputTree = new TTree("BaselineTree", "baseline tree");
        this->outputTree->SetDirectory(outputFile);

        cout << padString("Opening write file") << this->out_file_name << endl;
        cout << endl;
    }

    void initRandomGenerators() {
        unsigned random_seed = std::chrono::system_clock::now().time_since_epoch().count();
        this->random_generator = new default_random_engine(random_seed);
        myRandom.SetSeed(0);
    }

    void initHistograms() {
        /**
         * Set up photon and Z METl histograms, and smeared METl histograms, all binned in pT.
         * Set up Z mll histograms, binned in pT and METl.
         */
        for (int bin=0; bin<bins::n_pt_bins+2; bin++) {
            this->hist_z_metl_bin_pt[bin] = new TH1D(TString("hist_z_metl_")+TString::Itoa(bin,10),"",
                bins::n_smearing_bins,bins::smearing_low,bins::smearing_high);
            this->hist_g_metl_bin_pt[bin] = new TH1D(TString("hist_g_metl_")+TString::Itoa(bin,10),"",
                bins::n_smearing_bins,bins::smearing_low,bins::smearing_high);

            this->hist_g_smeared_metl_bin_pt[bin] = new TH1D(TString("hist_g_smeared_metl_")+TString::Itoa(bin,10),
                "",bins::n_smearing_bins,bins::smearing_low,bins::smearing_high);

            vector<TH1D*> hist_z_mll_bin_pt_metl_pt_bin;
            for (int bin1=0; bin1<bins::n_METl_bins+2; bin1++) {
                hist_z_mll_bin_pt_metl_pt_bin.push_back(new TH1D(TString("hist_z_Mll_dPt_")+TString::Itoa(bin,10)
                    +TString("_")+TString::Itoa(bin1,10),"",bins::n_mll_bins,bins::mll_bin));
            }
            this->hist_z_mll_bin_pt_metl.push_back(hist_z_mll_bin_pt_metl_pt_bin);
        }
    }

    Options getSmearingSelectionsWithChannel(Options options) {
        //--- modify event selections and weights
        options.bkg_smearing_selection = cuts::selections["bkg_baseline"] + cuts::selections[options.channel];
        options.photon_smearing_selection = cuts::selections["photon_baseline"];

        cout << padString("bkg selection") << options.bkg_smearing_selection.GetTitle() << endl;
        cout << padString("bkg weight") << cuts::bkg_weight.GetTitle() << endl;
        cout << padString("photon selection") << options.photon_smearing_selection.GetTitle() << endl;
        cout << padString("photon weight") << cuts::photon_weight.GetTitle() << endl;
        cout << endl;

        return options;
    }

    PhotonToZConverter(Options options) {
        this->options = options;

        //--- open files
        this->openPhotonFile(options);
        this->openOutputFile(options);

        //--- set up random generators
        this->initRandomGenerators();

        //--- set up histograms
        this->initHistograms();

        //--- get lepton histograms and smearing Gaussians
        options = getSmearingSelectionsWithChannel(options);
        this->getLepThetaHistograms(options);
        this->smearing_gaussians = this->getSmearingGaussians(options);
    }

    //---------------
    // Z MC FEATURES
    //---------------

    map<string, map<int, TH1D*>> zmc_hists;

    void fillZMCFeatureHists(Options options, TTree* ttree_zjets) {
        /// Compare lep eta, METl, and mll for photon vs. Z MC

        //--- initialize histograms
        initFeatureHists(zmc_hists);
        zmc_hists.erase("METl_raw");
        zmc_hists.erase("met_Et_raw");

        //--- fill histograms
        ttree_zjets->Draw(">>event_list", options.bkg_smearing_selection);
        auto event_list = (TEventList*) gDirectory->Get("event_list");

        ttree_zjets->SetBranchStatus("*", 0);
        vector<float> *lepEta = new vector<float>(10); SetInputBranch(ttree_zjets, "lepEta", &lepEta);
        vector<float> *lepPt = new vector<float>(10); SetInputBranch(ttree_zjets, "lepPt", &lepPt);
        float METl; SetInputBranch(ttree_zjets, "METl", &METl);
        float METt; SetInputBranch(ttree_zjets, "METt", &METt);
        float met_Et; SetInputBranch(ttree_zjets, "met_Et", &met_Et);
        float mll; SetInputBranch(ttree_zjets, "mll", &mll);
        float Ptll; SetInputBranch(ttree_zjets, "Ptll", &Ptll);
        double totalWeight; SetInputBranch(ttree_zjets, "totalWeight", &totalWeight);

        for (int i=0; i<event_list->GetN(); i++) {
            ttree_zjets->GetEntry(event_list->GetEntry(i));
            int pt_bin = bins::hist_pt_bins->FindBin(Ptll);
            //if (pt_bin == 0) cout << Ptll << endl;

            this->zmc_hists["lepEta"][pt_bin]->Fill(lepEta->at(0), totalWeight);
            this->zmc_hists["lepEta"][pt_bin]->Fill(lepEta->at(1), totalWeight);
            this->zmc_hists["lepPt"][pt_bin]->Fill(lepPt->at(0), totalWeight);
            this->zmc_hists["lepPt"][pt_bin]->Fill(lepPt->at(1), totalWeight);
            this->zmc_hists["met_Et"][pt_bin]->Fill(met_Et, totalWeight);
            this->zmc_hists["METl"][pt_bin]->Fill(METl, totalWeight);
            this->zmc_hists["METt"][pt_bin]->Fill(METt, totalWeight);
            this->zmc_hists["mll"][pt_bin]->Fill(mll, totalWeight);
        }

        //--- store histograms
        for (int pt_bin=0; pt_bin<bins::n_pt_bins+2; pt_bin++)
            this->zmc_hists["lep_cm_theta"][pt_bin] = this->h_lep_cm_theta;
    }

    map<string, map<int, TH1D*>> getZMCFeatureHists() {
        return this->zmc_hists;
    };

    //------------------
    // HELPER FUNCTIONS
    //------------------

    int rebinHistogram(TH1D* hist, int rebin) {
        /**
         * If rebin=0, rebin histogram by factors of 2, up to 32, until it looks reasonable or there are an odd
         * number of bins. This means having abs(negative_yield/positive_yield)<0.005 and
         * core_yield/positive_yield>0.4. Else, rebin by #(rebin).
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
            while ((abs(negative_yield/positive_yield)>0.005 || core_yield/positive_yield<0.4) &&
            has_odd_nbins==0 && rebin<=32) {
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

    //------------------
    // LEPTON SPLITTING
    //------------------

    void getLepThetaHistograms(Options options) {
        cout << PBLU("Getting Z lepton CM theta distribution histogram") << endl;
        cout << endl;
        gStyle->SetOptStat(0);

        //--- open files and create TChains
        TTree *ttree_zjets;
        string zjets_filename = options.reduction_folder + options.mc_period + "_Zjets.root";

        cout << padString("Opening Z+jets file") << zjets_filename << endl;
        ttree_zjets = (TTree*)(new TFile(zjets_filename.c_str()))->Get("BaselineTree");
        cout << padString("Z+jets entries") << ttree_zjets->GetEntries() << endl;
        cout << endl;

        //--- fill lep theta histogram
        this->h_lep_cm_theta = new TH1D("h_lep_cm_theta", "", 100, 0, 3.2);

        ttree_zjets->Draw("Z_cm_lep_theta>>h_lep_cm_theta", options.bkg_smearing_selection*cuts::bkg_weight, "goff");
        cout << padString("Z+jets integral") << h_lep_cm_theta->Integral(0, h_lep_cm_theta->GetNbinsX()+2) << endl;
        cout << endl;

        if (options.make_diagnostic_plots || options.unit_testing) {
            this->fillZMCFeatureHists(options, ttree_zjets);
        }

        //--- get lep theta histogram bin boundaries
        std::vector<int> h_lep_cm_bin_counts;
        for (int i=0; i<this->h_lep_cm_theta->GetNbinsX(); i++) {
            h_lep_cm_bin_counts.push_back(this->h_lep_cm_theta->GetBinContent(i));
            this->lep_cm_theta_bin_bounds.push_back(this->h_lep_cm_theta->GetBinLowEdge(i));
        }
        this->lep_cm_theta_bin_bounds.push_back(this->h_lep_cm_theta->GetBinLowEdge(this->h_lep_cm_theta->GetNbinsX()) +
            this->h_lep_cm_theta->GetBinWidth(this->h_lep_cm_theta->GetNbinsX()));
        this->lep_cm_theta_distribution = new std::discrete_distribution<int>(h_lep_cm_bin_counts.begin(),h_lep_cm_bin_counts.end());
    }

    float getRandomLepTheta() {
        int lep_theta_bin = (*this->lep_cm_theta_distribution)(*this->random_generator);
        float low_lep_theta = this->lep_cm_theta_bin_bounds[lep_theta_bin];
        float high_lep_theta = this->lep_cm_theta_bin_bounds[lep_theta_bin+1];
        float lep_theta_cm = this->myRandom.Rndm()*(high_lep_theta-low_lep_theta) + low_lep_theta;
        return lep_theta_cm;
    }

    tuple<bool, TLorentzVector, TLorentzVector> splitLeptons(TLorentzVector z_4vec, Options options) {
        // boost along z axis (since we measure angles in CM relative to boost direction)
        //TVector3 boost_vec(0, 0, z_4vec.BoostVector().Mag());
        TVector3 boost_vec = z_4vec.BoostVector();

        TRandom1 myRandom;
        myRandom.SetSeed(0);

        //float lep_M = 0;
        //if (options.channel == "mm") lep_M = 0.1056583;
        //else if (options.channel == "ee") lep_M = 0.0005109;
        //else failTest("Cannot split in this channel");

        TLorentzVector l0_lab_4vec, l1_lab_4vec;
        int n_tries = 0;
        bool good_event = false;
        while (n_tries++ < 10) {
            double lep_phi_cm = myRandom.Rndm()*2.*TMath::Pi();
            //float lep_theta_cm = this->getRandomLepTheta(); // Histogram sampling
            float lep_theta_cm = acos(1 - 2*myRandom.Rndm()); // uniform sampling on a sphere

            // Split leptons in Z rest frame
            TLorentzVector l0_cm_4vec, l1_cm_4vec;
            double lep_E_cm = z_4vec.M()/2.;
            double lep_px_cm = lep_E_cm*TMath::Sin(lep_theta_cm)*TMath::Cos(lep_phi_cm);
            double lep_py_cm = lep_E_cm*TMath::Sin(lep_theta_cm)*TMath::Sin(lep_phi_cm);
            double lep_pz_cm = lep_E_cm*TMath::Cos(lep_theta_cm);
            l0_cm_4vec.SetPxPyPzE(lep_px_cm, lep_py_cm, lep_pz_cm, lep_E_cm);
            l1_cm_4vec.SetPxPyPzE(-lep_px_cm, -lep_py_cm, -lep_pz_cm, lep_E_cm);
            //l0_cm_4vec.SetPtEtaPhiM(lep_E_cm, convertThetaToEta(lep_theta_cm), lep_phi_cm, lep_M);
            //l1_cm_4vec.SetPtEtaPhiM(lep_E_cm, convertThetaToEta(M_PI-lep_theta_cm), M_PI+lep_phi_cm, lep_M);

            // Boost to lab frame using smeared photon pT, eta, and phi
            l0_lab_4vec = l0_cm_4vec;
            l1_lab_4vec = l1_cm_4vec;
            l0_lab_4vec.Boost(boost_vec);
            l1_lab_4vec.Boost(boost_vec);

            // Rotate to lab coordinates
            //l0_lab_4vec.RotateZ(z_4vec.Phi());
            //l0_lab_4vec.RotateY(z_4vec.Theta());
            //l1_lab_4vec.RotateZ(z_4vec.Phi());
            //l1_lab_4vec.RotateY(z_4vec.Theta());

            if (l0_lab_4vec.Pt() < l1_lab_4vec.Pt()) {
                TLorentzVector lep_placeholder = l1_lab_4vec;
                l1_lab_4vec = l0_lab_4vec;
                l0_lab_4vec = lep_placeholder;
            }

            //cout << z_4vec.E() << " " << (l0_lab_4vec + l1_lab_4vec).E() << " " << l0_lab_4vec.E() << " " << l1_lab_4vec.E() << endl;
            //cout << z_4vec.Eta() << " " << (l0_lab_4vec + l1_lab_4vec).Eta() << " " << l0_lab_4vec.Eta() << " " << l1_lab_4vec.Eta() << endl;

            good_event = abs(l0_lab_4vec.Eta()) < 2.5 && abs(l1_lab_4vec.Eta()) < 2.5 &&
                  l0_lab_4vec.Pt() > cuts::leading_lep_pt_cut && l1_lab_4vec.Pt() > cuts::second_lep_pt_cut;
            if (good_event) break;
        }

        return make_tuple(good_event, l0_lab_4vec, l1_lab_4vec);
    }

    //--------------
    // MET SMEARING
    //--------------

    void getMETlAndMllHistograms() {
        /**
         * Fills histogram arrays with Z/photon METl distributions (and Z mll distribution) in binned regions.
         */

        cout << PBLU("Getting METl and mll histograms") << endl;
        cout << endl;

        //--- set up histograms and input files
        TFile* data_file = new TFile((this->options.reduction_folder + this->options.data_period + "_data_bkg.root").c_str());
        TFile* ttbar_mc_file = new TFile((this->options.reduction_folder + this->options.mc_period + "_ttbar.root").c_str());
        TFile* diboson_mc_file = new TFile((this->options.reduction_folder + this->options.mc_period + "_diboson.root").c_str());
        TFile* Z_mc_file = new TFile((this->options.reduction_folder + this->options.mc_period + "_Zjets.root").c_str());
        TFile* photon_file = new TFile(this->in_file_name.c_str());

        //--- z samples
        if (true) { // lazy scoping
            TTree* tree = (TTree*)Z_mc_file->Get("BaselineTree");
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

            tree->Draw(">>event_list", this->options.bkg_smearing_selection);
            //tree->Draw(">>event_list", "nJet30>=2 && lepPt[0]>25 && lepPt[1]>25");
            //tree->Draw(">>event_list", "nJet30>=2");
            auto event_list = (TEventList*) gDirectory->Get("event_list");
            for (int entry=0; entry<event_list->GetN(); entry++) {
                tree->GetEntry(event_list->GetEntry(entry));
                if (TString(this->options.channel).EqualTo("ee") && channel != 1) continue;
                if (TString(this->options.channel).EqualTo("mm") && channel != 0) continue;
                int pt_bin = bins::hist_pt_bins->FindBin(ptll);
                int METl_bin = bins::hist_METl_bins->FindBin(METl);
                this->hist_z_mll_bin_pt_metl[pt_bin][METl_bin]->Fill(mll, totalWeight);
                //if (ptll<50. || jet_n<2 || lep_pT->at(0)<cuts::leading_lep_pt_cut || lep_pT->at(1)<cuts::second_lep_pt_cut) continue;
                //if (jet_n<2 || lep_pT->at(0)<cuts::leading_lep_pt_cut || lep_pT->at(1)<cuts::second_lep_pt_cut) continue;
                this->hist_z_metl_bin_pt[pt_bin]->Fill(METl, totalWeight);
            }
        }

        //--- photon sample
        if (true) { // lazy scoping
            TTree* tree = (TTree*)photon_file->Get("BaselineTree");
            tree->SetBranchStatus("*", 0);
            double totalWeight; SetInputBranch(tree, "totalWeight", &totalWeight);
            int jet_n; SetInputBranch(tree, "nJet30", &jet_n);
            int bjet_n; SetInputBranch(tree, "bjet_n", &bjet_n);
            float ptll; SetInputBranch(tree, "gamma_pt", &ptll);
            int nLep_signal; SetInputBranch(tree, "nLep_signal", &nLep_signal);
            float METl; SetInputBranch(tree, "METl_unsmeared", &METl);

            tree->Draw(">>event_list", this->options.photon_smearing_selection);
            //tree->Draw(">>event_list", "nJet30>=2");
            auto event_list = (TEventList*) gDirectory->Get("event_list");
            for (int entry=0; entry<event_list->GetN(); entry++) {
                tree->GetEntry(event_list->GetEntry(entry));
                //if (ptll<50. || jet_n!=1 || bjet_n!=0) continue;
                int pt_bin = bins::hist_pt_bins->FindBin(ptll);
                this->hist_g_metl_bin_pt[pt_bin]->Fill(METl, totalWeight);
            }
        }

        //--- close files
        data_file->Close();
        ttbar_mc_file->Close();
        diboson_mc_file->Close();
        Z_mc_file->Close();
        photon_file->Close();
    };

    map<int, normal_distribution<float>> getSmearingGaussians(Options options) {
        /**
         * Returns map with key = photon pt bin, value = (mean, std).
         * For a photon in a given pt bin, smear the event's MET using the given Gaussian numbers.
         */

        //--- METl histograms
        this->getMETlAndMllHistograms();

        //--- perform smearing
        cout << PBLU("Calculating smearing Gaussians") << endl;

        map<int, normal_distribution<float>> smearing_gaussians;
        for (int pt_bin=0; pt_bin<bins::n_pt_bins+2; pt_bin++) {
            float smear_mean = 0;
            float smear_rms = 0;

            float z_int = this->hist_z_metl_bin_pt[pt_bin]->Integral(
                0, this->hist_z_metl_bin_pt[pt_bin]->GetNbinsX()+1);
            float g_int = this->hist_g_metl_bin_pt[pt_bin]->Integral(
                0, this->hist_g_metl_bin_pt[pt_bin]->GetNbinsX()+1);
            if (z_int > 0 && g_int > 0) {
                //--- Rebin METl histograms to look reasonable
                int rebin_factor = rebinHistogram(this->hist_z_metl_bin_pt[pt_bin], 0);
                rebinHistogram(this->hist_g_metl_bin_pt[pt_bin], rebin_factor);
                for (int met_bin=0; met_bin<bins::n_METl_bins+2; met_bin++)
                    rebinHistogram(this->hist_z_mll_bin_pt_metl[pt_bin][met_bin], 0);
                int n_hist_bins = bins::n_smearing_bins/rebin_factor;

                //--- Save smearing value
                this->hist_z_metl_bin_pt[pt_bin]->Fit("gaus", "Q");
                this->hist_g_metl_bin_pt[pt_bin]->Fit("gaus", "Q");
                TF1 *z_fit = (TF1*)this->hist_z_metl_bin_pt[pt_bin]->GetFunction("gaus");
                TF1 *g_fit = (TF1*)this->hist_g_metl_bin_pt[pt_bin]->GetFunction("gaus");
                smear_mean = z_fit->GetParameter(1) - g_fit->GetParameter(1);
                smear_rms = sqrt(pow(z_fit->GetParameter(2),2) - pow(g_fit->GetParameter(2),2));
                //smear_mean = z_fit->GetParameter(1);
                //smear_rms = z_fit->GetParameter(2);
                if (isnan(smear_rms)) smear_rms = 0;
            }
            
            if (options.channel=="ee") smearing_gaussians[pt_bin] = normal_distribution<float>(smear_mean, 0.0);
            else smearing_gaussians[pt_bin] = normal_distribution<float>(smear_mean, smear_rms);
        }
        cout << endl;

        return smearing_gaussians;
    }

    tuple<float, float, float> smearMETlAndMll(float METl, float gamma_pt) {
        int pt_bin = bins::hist_pt_bins->FindBin(gamma_pt);
        auto smearing_gaussian = this->smearing_gaussians[pt_bin];
            
        float METl_smeared;
        if (this->options.turn_off_shifting_and_smearing)
            METl_smeared = METl;
        else {
            METl_smeared = METl + smearing_gaussian(*random_generator);
            //METl_smeared = smearing_gaussian(*random_generator);
        }
        this->hist_g_smeared_metl_bin_pt[pt_bin]->Fill(METl_smeared);

        int METl_bin = bins::hist_METl_bins->FindBin(METl_smeared);
        float mll = 91.188;
        // fix for zero integral error
        //if (this->hist_z_mll_bin_pt_metl[pt_bin][METl_bin]->Integral(0, bins::n_mll_bins+1)>0)
        if (this->hist_z_mll_bin_pt_metl[pt_bin][METl_bin]->Integral()>0)
            mll = this->hist_z_mll_bin_pt_metl[pt_bin][METl_bin]->GetRandom();

        //float gamma_pt_smeared = gamma_pt - (METl_smeared-METl);
        float gamma_pt_smeared = gamma_pt;

        return make_tuple(METl_smeared, mll, gamma_pt_smeared);
    }

    //------------------
    // EVENT CONVERSION
    //------------------

    void convertEvents() {
        TTree *inputTree = this->inputTree;
        TTree *outputTree = this->outputTree;

        float gamma_phi; inputTree->SetBranchAddress("gamma_phi", &gamma_phi);

        double totalWeight; CopyBranch(inputTree, outputTree, "totalWeight", "totalWeight", &totalWeight, "D");
        int bjet_n; CopyBranch(inputTree, outputTree, "bjet_n", "bjet_n", &bjet_n, "I");
        float gamma_pt; CopyBranch(inputTree, outputTree, "gamma_pt", "gamma_pt",  &gamma_pt, "F");
        float gamma_eta; CopyBranch(inputTree, outputTree, "gamma_eta", "Z_eta",  &gamma_eta, "F");
        float METl; CopyBranch(inputTree, outputTree, "METl_unsmeared", "METl_unsmeared", &METl, "F");
        float METt; CopyBranch(inputTree, outputTree, "METt_unsmeared", "METt_unsmeared", &METt, "F");
        float METl_smeared; outputTree->Branch("METl", &METl_smeared, "METl/F");
        float METt_smeared; outputTree->Branch("METt", &METt_smeared, "METt/F");
        float HT; CopyBranch(inputTree, outputTree, "Ht30", "Ht30", &HT, "F");
        float MET_raw; CopyBranch(inputTree, outputTree, "met_Et_unsmeared", "met_Et_unsmeared", &MET_raw, "F");
        float MET_phi; CopyBranch(inputTree, outputTree, "met_Phi", "met_Phi", &MET_phi, "F");

        vector<float>* jet_pT = new vector<float>(10); CopyBranch(inputTree, outputTree, "jetPt", "jetPt", &jet_pT, "vector<float>");
        vector<float>* jet_eta = new vector<float>(10); CopyBranch(inputTree, outputTree, "jet_eta", "jet_eta", &jet_eta, "vector<float>");
        vector<float>* jet_phi = new vector<float>(10); CopyBranch(inputTree, outputTree, "jet_phi", "jet_phi", &jet_phi, "vector<float>");

        float gamma_pt_smeared; outputTree->Branch("Ptll", &gamma_pt_smeared, "Ptll/F");
        float gamma_phi_smeared; outputTree->Branch("Z_phi", &gamma_phi_smeared, "Z_phi/F");
        int lep_n; outputTree->Branch("lep_n", &lep_n, "lep_n/I"); outputTree->Branch("nLep_signal", &lep_n, "nLep_signal/I");
            outputTree->Branch("nLep_base", &lep_n, "nLep_base/I");
        float MET_smeared; outputTree->Branch("met_Et", &MET_smeared, "met_Et/F");
        float DPhi_METLepLeading_smeared; outputTree->Branch("DPhi_METLepLeading", &DPhi_METLepLeading_smeared, "DPhi_METLepLeading/F");
        float DPhi_METLepSecond_smeared; outputTree->Branch("DPhi_METLepSecond", &DPhi_METLepSecond_smeared, "DPhi_METLepSecond/F");
        float DPhi_METZPhoton_smear; outputTree->Branch("DPhi_METZPhoton", &DPhi_METZPhoton_smear, "DPhi_METZPhoton/F");
        float MT2; outputTree->Branch("mt2leplsp_0", &MT2, "mt2leplsp_0/F");
        float DR_2Lep; outputTree->Branch("DR_2Lep", &DR_2Lep, "DR_2Lep/F");
        int photon_conversion_type; CopyBranch(inputTree, outputTree, "PhotonConversionType", "PhotonConversionType", &photon_conversion_type, "I");
        float lep_theta_cm; outputTree->Branch("Z_cm_lep_theta", &lep_theta_cm, "Z_cm_lep_theta/F");
        float dPhiPllMet; outputTree->Branch("dPhiPllMet", &dPhiPllMet, "dPhiPllMet/F");

        //--- HistFitter branches
        vector<string> histFitterBranches {"DatasetNumber/I", "H2PP/D", "H5PP/D", "H5PP_VR/D",
            "METOverPtISR/F", "METOverPtW/F", "METOverPtZ/F", "MJ/D", "MJ_VR/D", "MZ/D", "MZ_VR/D", "NjISR/D",
            "NjS/D", "PTCM/D", "PTCM_VR/D", "PTI/D", "PTISR/D", "PTISR_VR/D", "PTI_VR/D", "RISR/D", "RISR_VR/D",
            "RPT_HT5PP/D", "RPT_HT5PP_VR/D", "R_minH2P_minH3P/D", "R_minH2P_minH3P_VR/D", "Rjj/F", "Rll/F",
            "dPhiMetISR/F", "dPhiMetJet1/F", "dPhiMetJet2/F", "dPhiMetJet12Min/F", "dPhiPjjMet/F", "dPhiPllMet/F",
            "dphiISRI/D", "dphiISRI_VR/D", "dphiVP/D", "dphiVP_VR/D", "lept1Pt_VR/D", "lept2Pt_VR/D", "mTl3/D",
            "met_Sign/F", "minDphi/D", "minDPhi2JetsMet/F", "mll_RJ/D", "mll_RJ_VR/D", "nJet30/I", "nJet20/I", "mjj/F",
            "nBJet20_MV2c10_FixedCutBEff_77/I", "trigMatch_2LTrigOR/I", "genWeight/D", "eventWeight/D", "leptonWeight/D",
            "jvtWeight/D", "bTagWeight/D", "pileupWeight/D", "globalDiLepTrigSF/D", "RunNumber/I", "RandomRunNumber/I",
            "trigMatch_2LTrig/I", "lumi/D", "mjj_minDPhiZMET/F", "mbb/F", "PtISR/F"};
        CopyAllBranches(inputTree, outputTree, histFitterBranches);

        vector<float>* dPhiMetJet = new vector<float>(10); CopyBranch(inputTree, outputTree, "dPhiMetJet", "dPhiMetJet", &dPhiMetJet, "vector<float>");
        vector<float>* jetM = new vector<float>(10); CopyBranch(inputTree, outputTree, "jetM", "jetM", &jetM, "vector<float>");
        float mll; outputTree->Branch("mll", &mll, "mll/F");
        int is_OS = 1; outputTree->Branch("is_OS", &is_OS, "is_OS/I");
        vector<float>* lep_pT = new vector<float>(10); outputTree->Branch("lepPt", "vector<float>", &lep_pT);
        vector<float>* lep_eta = new vector<float>(10); outputTree->Branch("lepEta", "vector<float>", &lep_eta);
        vector<float>* lep_phi = new vector<float>(10); outputTree->Branch("lepPhi", "vector<float>", &lep_phi);
        vector<int>* lep_charge = new vector<int>(10); outputTree->Branch("lepCharge", "vector<int>", &lep_charge);
        vector<float>* lep_m = new vector<float>(10); outputTree->Branch("lepM", "vector<float>", &lep_m);
        vector<int>* lepIsoFCTight = new vector<int>{1,1}; outputTree->Branch("lepIsoFCTight", "vector<int>", &lepIsoFCTight);
        vector<int>* lepIsPR = new vector<int>{1,1}; outputTree->Branch("lepIsPR", "vector<int>", &lepIsPR);
        Int_t lepChannel; outputTree->Branch("channel", &lepChannel, "channel/I");
        vector<int>* lep_flavor = new vector<int>(10); outputTree->Branch("lepFlavor", "vector<int>", &lep_flavor);

        //-----------------------------
        // set global lepton info
        //-----------------------------

        int flavor;
        if (TString(options.channel).EqualTo("ee")) {
            flavor = 1;
            lepChannel = 1;
        }
        else if (TString(options.channel).EqualTo("mm")) {
            flavor = 2;
            lepChannel = 0;
        }

        lep_flavor->clear();
        lep_flavor->push_back(flavor);
        lep_flavor->push_back(flavor);

        lep_m->clear();
        if (lepChannel == 0) {
            lep_m->push_back(0.1056583);
            lep_m->push_back(0.1056583);
        }
        else if (lepChannel == 1) {
            lep_m->push_back(0.0005109);
            lep_m->push_back(0.0005109);
        }

        //-----------------------------
        // loop over events
        //-----------------------------

        Long64_t nentries = inputTree->GetEntries();
        for (Long64_t i=0; i<nentries; i++) {
            if (fmod(i,1e5)==0) cout << i << " events processed.\r" << flush;
            inputTree->GetEntry(i);

            //--- smearing
            gamma_pt_smeared = gamma_pt;
            gamma_phi_smeared = gamma_phi;

            METt_smeared = METt;
            try {
                auto [METl_smeared_return, mll_return, gamma_pt_return] = this->smearMETlAndMll(METl, gamma_pt);
                METl_smeared = METl_smeared_return;
                mll = mll_return;
                gamma_pt = gamma_pt_return;
            }
            catch(...) {
                cout << PRED("Problem in event ") << i << endl;
                continue;
            }
            if (mll==0) continue;

            MET_smeared = sqrt(pow(METl_smeared, 2) + pow(METt_smeared, 2));
            TLorentzVector MET_smeared_4vec;
            MET_smeared_4vec.SetPtEtaPhiM(MET_smeared, 0, MET_phi, 0);
            DPhi_METZPhoton_smear = gamma_phi_smeared - MET_smeared_4vec.Phi();

            //--- lepton splitting
            TLorentzVector z_4vec;
            z_4vec.SetPtEtaPhiM(gamma_pt_smeared, gamma_eta, gamma_phi_smeared, mll);

            auto [good_event, l0_lab_4vec, l1_lab_4vec] = splitLeptons(z_4vec, options);
            if (!good_event) continue;

            lep_pT->clear();
            lep_pT->push_back(l0_lab_4vec.Pt());
            lep_pT->push_back(l1_lab_4vec.Pt());
            lep_eta->clear();
            lep_eta->push_back(l0_lab_4vec.Eta());
            lep_eta->push_back(l1_lab_4vec.Eta());
            lep_phi->clear();
            lep_phi->push_back(l0_lab_4vec.Phi());
            lep_phi->push_back(l1_lab_4vec.Phi());
            int charge = myRandom.Integer(2)*2-1; // random leading lepton charge
            lep_charge->clear();
            lep_charge->push_back(charge);
            lep_charge->push_back(-charge);

            lep_n = 2;
            MT2 = ComputeMT2(l0_lab_4vec, l1_lab_4vec, MET_smeared_4vec, 0, 0).Compute();
            DPhi_METLepLeading_smeared = fabs(MET_smeared_4vec.DeltaPhi(l0_lab_4vec));
            DPhi_METLepSecond_smeared = fabs(MET_smeared_4vec.DeltaPhi(l1_lab_4vec));
            DR_2Lep = l0_lab_4vec.DeltaR(l1_lab_4vec);
            dPhiPllMet = fabs(MET_smeared_4vec.DeltaPhi(l0_lab_4vec + l1_lab_4vec));

            outputTree->Fill();
        }
        cout << endl;

        //-----------------------------
        // write
        //-----------------------------

        this->outputFile->cd();
        this->outputTree->Write();

        cout << PBLU("Done with smearing") << endl;
        cout << endl;
        delete this->outputFile;
    }
};

//------------
// UNIT TESTS
//------------

Options setSmearingUnitTestOptions(Options options) {
    options.reduction_folder = options.unit_test_folder + "ReducedNtuples/";
    options.smearing_folder = "./";

    options.period = "data15-16";
    options.mc_period = getMCPeriod(options.period);
    options.data_period = DataPeriod(options.period);
    options.is_data = true;
    options.channel = "mm";

    options.turn_off_shifting_and_smearing = false;
    options.make_diagnostic_plots = false;

    options.unit_testing = true;
    options.run_vgamma = false;

    return options;
}

void performSmearingUnitTests(Options options) {
    cout << BOLD(PBLU("Performing unit testing on smearing step")) << endl;
    cout << endl;

    //--- get Z MC feature plots and set up photon feature plots for comparison
    PhotonToZConverter converter(options);
    map<string, map<int, TH1D*>> zmc_plots = converter.getZMCFeatureHists();

    TCanvas *can = new TCanvas("can","can",600,600);
    map<string, map<int, TH1D*>> photon_plots;

    TPad* mainpad = new TPad("mainpad","mainpad",0.0,0.0,1.0,0.8);
    mainpad->Draw();

    TPad* ratio_pad = new TPad("ratio_pad","ratio_pad",0.0,0.75,1.0,0.905);
    ratio_pad->Draw();
    ratio_pad->SetGridy();

    initFeatureHists(photon_plots);

    //--- open photon files
    string in_file_name;
    if (options.is_data) in_file_name = options.reduction_folder + options.data_period + "_data_photon.root";
    else in_file_name = options.reduction_folder + options.mc_period + "_SinglePhoton222.root";
    TFile *photon_file = new TFile(in_file_name.c_str());
    TTree *photon_tree = (TTree*)photon_file->Get("BaselineTree");

    float gamma_pt; photon_tree->SetBranchAddress("gamma_pt", &gamma_pt);
    float gamma_eta; photon_tree->SetBranchAddress("gamma_eta", &gamma_eta);
    float gamma_phi; photon_tree->SetBranchAddress("gamma_phi", &gamma_phi);
    double totalWeight; photon_tree->SetBranchAddress("totalWeight", &totalWeight);
    float METl_unsmeared; photon_tree->SetBranchAddress("METl_unsmeared", &METl_unsmeared);
    float METt; photon_tree->SetBranchAddress("METt_unsmeared", &METt);

    //--- perform photon splitting and smearing
    float Z_m = 91.188;
    for (int i=0; i<photon_tree->GetEntries(); i++) {
        photon_tree->GetEntry(i);
        auto [METl, mll, gamma_pt_return] = converter.smearMETlAndMll(METl_unsmeared, gamma_pt);
        gamma_pt = gamma_pt_return;
        int pt_bin = bins::hist_pt_bins->FindBin(gamma_pt);

        photon_plots["lep_cm_theta"][pt_bin]->Fill(converter.getRandomLepTheta());

        TLorentzVector z_4vec;
        z_4vec.SetPtEtaPhiM(gamma_pt, gamma_eta, gamma_phi, Z_m);
        auto [good_event, l0_lab_4vec, l1_lab_4vec] = converter.splitLeptons(z_4vec, options);
        if (good_event) {
            photon_plots["lepEta"][pt_bin]->Fill(l0_lab_4vec.Eta(), totalWeight);
            photon_plots["lepEta"][pt_bin]->Fill(l1_lab_4vec.Eta(), totalWeight);
            photon_plots["lepPt"][pt_bin]->Fill(l0_lab_4vec.Pt(), totalWeight);
            photon_plots["lepPt"][pt_bin]->Fill(l1_lab_4vec.Pt(), totalWeight);
        }

        photon_plots["METl_raw"][pt_bin]->Fill(METl_unsmeared, totalWeight);
        photon_plots["METl"][pt_bin]->Fill(METl, totalWeight);
        photon_plots["METt"][pt_bin]->Fill(METt, totalWeight);
        photon_plots["met_Et_raw"][pt_bin]->Fill(sqrt(pow(METt,2) + pow(METl_unsmeared,2)), totalWeight);
        photon_plots["met_Et"][pt_bin]->Fill(sqrt(pow(METt,2) + pow(METl,2)), totalWeight);
        photon_plots["mll"][pt_bin]->Fill(mll, totalWeight);
    }

    //--- draw and save comparison plots
    for (auto& [key, hist] : zmc_plots) {
        THStack *zmc_stack = new THStack("zmc_stack","");
        THStack *photon_stack = new THStack("gmc_stack","");;

        for (int i=0; i<bins::n_pt_bins+2; i++) {
            //--- add to stack
            zmc_stack->Add(zmc_plots[key][i]);
            photon_stack->Add(photon_plots[key][i]);

            //--- draw Z MC
            mainpad->cd();
            zmc_plots[key][i]->Draw("hist");

            //--- scale and draw photons
            float z_yield = zmc_plots[key][i]->Integral(0, zmc_plots[key][i]->GetNbinsX()+1);
            float g_yield = photon_plots[key][i]->Integral(0, photon_plots[key][i]->GetNbinsX()+1);

            //cout << key << " " << i << " " << z_yield << " " << g_yield << endl;
            float scale = g_yield != 0 ? z_yield/g_yield : 0;

            if (key == "METl") {
                float g_yield_raw = photon_plots["METl_raw"][i]->Integral(0, photon_plots["METl_raw"][i]->GetNbinsX()+1);
                float scale_raw = g_yield_raw != 0 ? z_yield/g_yield_raw : 0;

                photon_plots["METl_raw"][i]->Scale(scale_raw);
                photon_plots["METl_raw"][i]->SetLineColor(kGreen);
                photon_plots["METl_raw"][i]->Draw("hist same");
            }
            else if (key == "met_Et") {
                float g_yield_raw = photon_plots["met_Et_raw"][i]->Integral(0, photon_plots["met_Et_raw"][i]->GetNbinsX()+1);
                float scale_raw = g_yield_raw != 0 ? z_yield/g_yield_raw : 0;

                photon_plots["met_Et_raw"][i]->Scale(scale_raw);
                photon_plots["met_Et_raw"][i]->SetLineColor(kGreen);
                photon_plots["met_Et_raw"][i]->Draw("hist same");
            }

            TH1D *photon_plot_scaled = (TH1D*)photon_plots[key][i]->Clone();
            photon_plot_scaled->Scale(scale);

            //--- draw comparison
            photon_plot_scaled->SetLineColor(kRed);
            photon_plot_scaled->Draw("hist same");

            //--- draw legend
            TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
            leg->AddEntry(zmc_plots[key][i], "Zjets", "f");
            if (key == "METl")
                leg->AddEntry(photon_plots["METl_raw"][i], "photon raw", "f");
            leg->AddEntry(photon_plot_scaled, "photon corrected", "f");
            leg->SetBorderSize(0);
            leg->SetFillColor(0);
            leg->Draw();

            //--- draw ratio
            TH1D *hratio = (TH1D*) zmc_plots[key][i]->Clone(("hratio_"+to_string(i)).c_str());
            hratio->Divide(photon_plot_scaled);

            hratio->SetMarkerStyle(20);
            hratio->SetMarkerColor(kRed);
            hratio->SetLineColor(kRed);
            hratio->GetXaxis()->SetTitle("");
            hratio->GetXaxis()->SetLabelSize(0.);
            hratio->GetYaxis()->SetNdivisions(5);
            hratio->GetYaxis()->SetTitle("");
            hratio->GetYaxis()->SetTitleSize(0.15);
            hratio->GetYaxis()->SetTitleOffset(0.3);
            hratio->GetYaxis()->SetLabelSize(0.15);
            hratio->SetMinimum(0.0);
            hratio->SetMaximum(2.0);
            hratio->GetYaxis()->SetRangeUser(0.0,2.0);
            hratio->SetTitle("");

            ratio_pad->cd();
            hratio->Draw();

            //--- save
            converter.outputFile->cd();
            can->Write(key + "_pt_bin_" + i + ".eps");
        }

        //--- scale photons for plotting
        TH1D *zmc_total = (TH1D*)zmc_stack->GetStack()->Last();
        TH1D *photon_total = (TH1D*)photon_stack->GetStack()->Last();

        float z_yield = zmc_total->Integral(0, zmc_total->GetNbinsX()+1);
        float g_yield = photon_total->Integral(0, photon_total->GetNbinsX()+1);
        float scale = g_yield != 0 ? z_yield/g_yield : 0;
        photon_total->Scale(scale);

        //--- draw total plot
        mainpad->cd();
        zmc_total->Draw("hist");
        photon_total->SetLineColor(kRed);
        photon_total->Draw("hist same");

        //--- draw total legend
        TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
        leg->AddEntry(zmc_total, "Zjets", "f");
        leg->AddEntry(photon_total, "photon corrected", "f");
        leg->SetBorderSize(0);
        leg->SetFillColor(0);
        leg->Draw();

        //--- draw total ratio
        TH1D *hratio = (TH1D*) zmc_total->Clone("hratio_total");
        hratio->Divide((TH1D*) photon_total);

        hratio->SetMarkerStyle(20);
        hratio->SetMarkerColor(kRed);
        hratio->SetLineColor(kRed);
        hratio->GetXaxis()->SetTitle("");
        hratio->GetXaxis()->SetLabelSize(0.);
        hratio->GetYaxis()->SetNdivisions(5);
        hratio->GetYaxis()->SetTitle("");
        hratio->GetYaxis()->SetTitleSize(0.15);
        hratio->GetYaxis()->SetTitleOffset(0.3);
        hratio->GetYaxis()->SetLabelSize(0.15);
        hratio->SetMinimum(0.0);
        hratio->SetMaximum(2.0);
        hratio->GetYaxis()->SetRangeUser(0.0,2.0);
        hratio->SetTitle("");

        ratio_pad->cd();
        hratio->Draw();

        //--- save
        converter.outputFile->cd();
        can->Write((key + ".eps").c_str());

        passTest("Comparison plots for " + key + " produced");
    }

    //remove(converter.out_file_name.c_str());
    delete can;

    passTest("Passed all unit tests");
    cout << endl;
}

//---------------
// MAIN FUNCTION
//---------------

void SmearPhotons(Options options) {
    if (options.unit_testing) {
        options = setSmearingUnitTestOptions(options);
        performSmearingUnitTests(options);
    }
    else {
        cout << BOLD(PBLU("Performing smearing")) << endl;
        cout << endl;

        //--- smear photons
        options.run_vgamma = false;
        PhotonToZConverter photon_converter(options);
        photon_converter.convertEvents();

        //--- smear Vgamma - this MC component must be subtracted from photon data
        if (options.is_data) {
            options.run_vgamma = true;
            PhotonToZConverter vgamma_converter(options);
            vgamma_converter.convertEvents();
        }
    }
}

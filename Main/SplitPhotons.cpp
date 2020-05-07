#include "Settings.cpp"

using namespace std;

//----------------
// SPLITTER CLASS
//----------------

class PhotonSplitter {
public:
    std::default_random_engine *random_generator;
    TRandom3 myRandom;

    TFile* input_file;
    TChain* input_tree;
    TFile* output_file;
    TTree* output_tree;

    vector<vector<ROOT::RDF::RResultPtr<TH1D>>> z_mll_histptrs;
    vector<vector<ROOT::RDF::RResultPtr<TH1D>>> z_cmleptheta_histptrs;

    //----------------
    // INITIALIZATION
    //----------------

    void initRandomGenerators() {
        unsigned random_seed = std::chrono::system_clock::now().time_since_epoch().count();
        this->random_generator = new default_random_engine(random_seed);
        myRandom.SetSeed(0);
    }

    void openReadAndWriteFiles(Options options) {
        //--- open input and output files
        string in_file_name, out_file_name;

        if (options.run_vgamma) {
            in_file_name = options.reduction_folder + options.mc_period + "_Vgamma.root";
            out_file_name = options.splitting_folder + options.mc_period + "_Vgamma_" + options.channel + ".root";
        }
        else {
            if (options.is_data) {
                in_file_name = options.reduction_folder + options.data_period + "_data_photon.root";
                out_file_name = options.splitting_folder + options.data_period + "_data_photon_" + options.channel + ".root"; 
            }
            else {
                in_file_name = options.reduction_folder + options.mc_period + "_SinglePhoton222.root";
                out_file_name = options.splitting_folder + options.mc_period + "_SinglePhoton222_" + options.channel + ".root";
            }
        }

        cout << padString("channel") << options.channel << endl;
        cout << padString("period") << options.period << endl;
        cout << padString("data-based photons?") << options.is_data << endl;
        cout << padString("splitting output") << out_file_name << endl;

        this->input_tree = new TChain(options.tree_name.c_str());
        this->input_tree->Add(in_file_name.c_str());

        cout << endl;
        cout << padString("Opening read file") << in_file_name << endl;
        cout << padString("Events in ntuple") << this->input_tree->GetEntries() << endl;

        this->output_file = new TFile(out_file_name.c_str(), "recreate");          
        this->output_tree = new TTree("BaselineTree", "baseline tree");
        this->output_tree->SetDirectory(this->output_file);

        cout << padString("Opening write file") << out_file_name << endl;
        cout << endl;
    }

    PhotonSplitter(Options options) {
        TH1::SetDefaultSumw2();
        this->initRandomGenerators();
        this->openReadAndWriteFiles(options);
    }

    //---------------
    // Z MC FEATURES
    //---------------

    void fillZMCHistograms(Options options, bool save_plots) {
        /**
         * Fill Z MC mll and lepcm histograms, binned in pT and METl.
         */
        cout << PBLU("Getting Z lepton CM theta distribution histogram") << endl;
        cout << endl;
        gStyle->SetOptStat(0);

        //--- open Z MC file
        string zjets_filename = options.reduction_folder + options.mc_period + "_Zjets.root";
        TTree *ttree_zjets = (TTree*)(new TFile(zjets_filename.c_str()))->Get(options.tree_name.c_str());
        ROOT::RDataFrame *zjets_df = new ROOT::RDataFrame(*ttree_zjets);
        cout << padString("Opening Z+jets file") << zjets_filename << endl;
        cout << padString("Z+jets entries") << ttree_zjets->GetEntries() << endl;
        cout << endl;

        //--- select events
        TCut bkg_splitting_selection = (cuts::selections["bkg_baseline"] + cuts::selections[options.channel])
            * cuts::bkg_weight;
        auto reduced_df = zjets_df->Filter(bkg_splitting_selection.GetTitle());

        //--- prep for filling histograms
        for (int pt_bin=0; pt_bin<bins::n_pt_bins+2; pt_bin++) {
            vector<ROOT::RDF::RResultPtr<TH1D>> z_mll_histptrs_row;
            vector<ROOT::RDF::RResultPtr<TH1D>> z_cmleptheta_histptrs_row;
            for (int METl_bin=0; METl_bin<bins::n_METl_bins+2; METl_bin++) {
                float low_pt = pt_bin>0 ? bins::pt_bins[pt_bin-1] : -999999;
                float high_pt = pt_bin<bins::n_pt_bins+1 ? bins::pt_bins[pt_bin] : 999999;
                TCut pt_filter = ("Ptll > " + to_string(low_pt) + " && Ptll < " + to_string(high_pt)).c_str();

                float low_METl = METl_bin>0 ? bins::METl_bins[METl_bin-1] : -999999;
                float high_METl = METl_bin<bins::n_METl_bins+1 ? bins::METl_bins[METl_bin] : 999999;
                TCut METl_filter = ("METl > " + to_string(low_METl) + " && METl < " + to_string(high_METl)).c_str();

                TCut bin_filter = pt_filter + METl_filter;
                z_mll_histptrs_row.push_back(zjets_df->Filter(bin_filter.GetTitle())
                    .Histo1D("mll", cuts::bkg_weight.GetTitle()));
                z_cmleptheta_histptrs_row.push_back(zjets_df->Filter(bin_filter.GetTitle())
                    .Histo1D("Z_cm_lep_theta", cuts::bkg_weight.GetTitle()));
            }
            z_mll_histptrs.push_back(z_mll_histptrs_row);
            z_cmleptheta_histptrs.push_back(z_cmleptheta_histptrs_row);
        }

        //--- fill and save histograms
        for (int pt_bin=0; pt_bin<bins::n_pt_bins+2; pt_bin++) {
            for (int METl_bin=0; METl_bin<bins::n_METl_bins+2; METl_bin++) {
                z_mll_histptrs[pt_bin][METl_bin]->Draw();
                z_cmleptheta_histptrs[pt_bin][METl_bin]->Draw();

                //--- save plots
                if (save_plots) {
                    string mll_name = "hist_z_mll_dPt_"+to_string(pt_bin)+"_dMETl_"+to_string(METl_bin);
                    string cmleptheta_name = "hist_z_cmleptheta_dPt_"+to_string(pt_bin)+"_dMETl_"+to_string(METl_bin);

                    this->output_file->cd();
                    z_mll_histptrs[pt_bin][METl_bin]->Write((mll_name + ".eps").c_str());
                    z_cmleptheta_histptrs[pt_bin][METl_bin]->Write((cmleptheta_name + ".eps").c_str());
                }
            }
        }
    }

    //------------------
    // LEPTON SPLITTING
    //------------------

    bool getRandomSample(vector<vector<ROOT::RDF::RResultPtr<TH1D>>> hists, float ptll, float METl, float &val) {
        int pt_bin = bins::hist_pt_bins->FindBin(ptll);
        int METl_bin = bins::hist_METl_bins->FindBin(METl);
        if (hists[pt_bin][METl_bin]->Integral() <= 0) return false;
        val = hists[pt_bin][METl_bin]->GetRandom();
        cout << pt_bin << " " << METl_bin << " " << val << endl;
        return true;
    }

    //tuple<bool, TLorentzVector, TLorentzVector> splitLeptons(TLorentzVector z_4vec) {
        //// boost along z axis (since we measure angles in CM relative to boost direction)
        //TVector3 boost_vec(0, 0, z_4vec.BoostVector().Mag());

        //TRandom1 myRandom;
        //myRandom.SetSeed(0);

        //TLorentzVector l0_lab_4vec, l1_lab_4vec;
        ////int n_tries = 0;
        ////while (n_tries++ < 10) {

            //double lep_phi_cm = myRandom.Rndm()*2.*TMath::Pi();
            //float lep_theta_cm = this->getRandomLepTheta(); // Histogram sampling

            //// Split leptons in Z rest frame
            //TLorentzVector l0_cm_4vec, l1_cm_4vec;
            //double lep_E_cm = z_4vec.M()/2.;
            //double lep_px_cm = lep_E_cm*TMath::Sin(lep_theta_cm)*TMath::Cos(lep_phi_cm);
            //double lep_py_cm = lep_E_cm*TMath::Sin(lep_theta_cm)*TMath::Sin(lep_phi_cm);
            //double lep_pz_cm = lep_E_cm*TMath::Cos(lep_theta_cm);
            //l0_cm_4vec.SetPxPyPzE(lep_px_cm, lep_py_cm, lep_pz_cm, lep_E_cm);
            //l1_cm_4vec.SetPxPyPzE(-lep_px_cm, -lep_py_cm, -lep_pz_cm, lep_E_cm);

            //// Boost to lab frame using split photon pT, eta, and phi
            //l0_lab_4vec = l0_cm_4vec;
            //l1_lab_4vec = l1_cm_4vec;
            //l0_lab_4vec.Boost(boost_vec);
            //l1_lab_4vec.Boost(boost_vec);
            //if (l0_lab_4vec.Pt() < l1_lab_4vec.Pt()) {
                //TLorentzVector lep_placeholder = l1_lab_4vec;
                //l1_lab_4vec = l0_lab_4vec;
                //l0_lab_4vec = lep_placeholder;
            //}

            //// Rotate to lab coordinates
            //l0_lab_4vec.RotateY(z_4vec.Theta());
            //l0_lab_4vec.RotateZ(z_4vec.Phi());
            //l1_lab_4vec.RotateY(z_4vec.Theta());
            //l1_lab_4vec.RotateZ(z_4vec.Phi());

            //bool good_event = abs(l0_lab_4vec.Eta()) < 2.5 && abs(l1_lab_4vec.Eta()) < 2.5 &&
                              //l0_lab_4vec.Pt() > cuts::leading_lep_pt_cut && l1_lab_4vec.Pt() > cuts::second_lep_pt_cut;
        ////}
        ////if (n_tries==11) continue;

        //return make_tuple(good_event, l0_lab_4vec, l1_lab_4vec);
    //}

    //tuple<float, float> splitMETlAndMll(float METl, float gamma_pt) {

        //float mll = 91.188;
        //// fix for zero integral error
        ////if (this->z_mll_histptrs[pt_bin][METl_bin]->Integral(0, bins::n_mll_bins+1)>0)
        //if (this->z_mll_histptrs[pt_bin][METl_bin]->Integral()>0)
            //mll = this->z_mll_histptrs[pt_bin][METl_bin]->GetRandom();

        //return make_tuple(0, mll);
    //}

    //---------------
    // SPLIT PHOTONS
    //---------------

    //void splitPhotons(Options options) {
        //TTree* inputTree = this->input_tree;
        //TTree* outputTree = this->output_tree;

        ////--- fill branches
        //float gamma_phi; inputTree->SetBranchAddress("gamma_phi", &gamma_phi);

        //double totalWeight; CopyBranch(inputTree, outputTree, "totalWeight", "totalWeight", &totalWeight, "D");
        //int bjet_n; CopyBranch(inputTree, outputTree, "bjet_n", "bjet_n", &bjet_n, "I");
        //float gamma_pt; CopyBranch(inputTree, outputTree, "gamma_pt", "gamma_pt",  &gamma_pt, "F");
        //float gamma_eta; CopyBranch(inputTree, outputTree, "gamma_eta", "Z_eta",  &gamma_eta, "F");
        //float METl; CopyBranch(inputTree, outputTree, "METl_unsplit", "METl_unsplit", &METl, "F");
        //float METt; CopyBranch(inputTree, outputTree, "METt_unsplit", "METt_unsplit", &METt, "F");
        //float METl_split; outputTree->Branch("METl", &METl_split, "METl/F");
        //float METt_split; outputTree->Branch("METt", &METt_split, "METt/F");
        //float HT; CopyBranch(inputTree, outputTree, "Ht30", "Ht30", &HT, "F");
        //float MET_raw; CopyBranch(inputTree, outputTree, "met_Et_unsplit", "met_Et_unsplit", &MET_raw, "F");
        //float MET_phi; CopyBranch(inputTree, outputTree, "met_Phi", "met_Phi", &MET_phi, "F");

        //vector<float>* jet_pT = new vector<float>(10); CopyBranch(inputTree, outputTree, "jetPt", "jetPt", &jet_pT, "vector<float>");
        //vector<float>* jet_eta = new vector<float>(10); CopyBranch(inputTree, outputTree, "jet_eta", "jet_eta", &jet_eta, "vector<float>");
        //vector<float>* jet_phi = new vector<float>(10); CopyBranch(inputTree, outputTree, "jet_phi", "jet_phi", &jet_phi, "vector<float>");

        //float gamma_pt_split; outputTree->Branch("Ptll", &gamma_pt_split, "Ptll/F");
        //float gamma_phi_split; outputTree->Branch("Z_phi", &gamma_phi_split, "Z_phi/F");
        //int lep_n; outputTree->Branch("lep_n", &lep_n, "lep_n/I"); outputTree->Branch("nLep_signal", &lep_n, "nLep_signal/I");
            //outputTree->Branch("nLep_base", &lep_n, "nLep_base/I");
        //float MET_split; outputTree->Branch("met_Et", &MET_split, "met_Et/F");
        //float DPhi_METLepLeading_split; outputTree->Branch("DPhi_METLepLeading", &DPhi_METLepLeading_split, "DPhi_METLepLeading/F");
        //float DPhi_METLepSecond_split; outputTree->Branch("DPhi_METLepSecond", &DPhi_METLepSecond_split, "DPhi_METLepSecond/F");
        //float DPhi_METZPhoton_split; outputTree->Branch("DPhi_METZPhoton", &DPhi_METZPhoton_split, "DPhi_METZPhoton/F");
        //float DR_2Lep; outputTree->Branch("DR_2Lep", &DR_2Lep, "DR_2Lep/F");
        //int photon_conversion_type; CopyBranch(inputTree, outputTree, "PhotonConversionType", "PhotonConversionType", &photon_conversion_type, "I");
        //float lep_theta_cm; outputTree->Branch("Z_cm_lep_theta", &lep_theta_cm, "Z_cm_lep_theta/F");

        ////--- HistFitter branches
        //vector<string> histFitterBranches {"DatasetNumber/I", "H2PP/D", "H5PP/D", "H5PP_VR/D",
            //"METOverPtISR/F", "METOverPtW/F", "METOverPtZ/F", "MJ/D", "MJ_VR/D", "MZ/D", "MZ_VR/D", "NjISR/D",
            //"NjS/D", "PTCM/D", "PTCM_VR/D", "PTI/D", "PTISR/D", "PTISR_VR/D", "PTI_VR/D", "RISR/D", "RISR_VR/D",
            //"RPT_HT5PP/D", "RPT_HT5PP_VR/D", "R_minH2P_minH3P/D", "R_minH2P_minH3P_VR/D", "Rjj/F", "Rll/F",
            //"dPhiMetISR/F", "dPhiMetJet1/F", "dPhiMetJet2/F", "dPhiMetJet12Min/F", "dPhiPjjMet/F", "dPhiPllMet/F",
            //"dphiISRI/D", "dphiISRI_VR/D", "dphiVP/D", "dphiVP_VR/D", "lept1Pt_VR/D", "lept2Pt_VR/D", "mTl3/D",
            //"met_Sign/F", "minDphi/D", "minDPhi2JetsMet/F", "mll_RJ/D", "mll_RJ_VR/D", "nJet30/I", "nJet20/I", "mjj/F",
            //"nBJet20_MV2c10_FixedCutBEff_77/I", "trigMatch_2LTrigOR/I", "genWeight/D", "eventWeight/D", "leptonWeight/D",
            //"jvtWeight/D", "bTagWeight/D", "pileupWeight/D", "globalDiLepTrigSF/D", "RunNumber/I", "RandomRunNumber/I",
            //"trigMatch_2LTrig/I", "lumi/D", "mjj_minDPhiZMET/F", "mbb/F", "PtISR/F"};
        //CopyAllBranches(inputTree, outputTree, histFitterBranches);

        //vector<float>* dPhiMetJet = new vector<float>(10); CopyBranch(inputTree, outputTree, "dPhiMetJet", "dPhiMetJet", &dPhiMetJet, "vector<float>");
        //vector<float>* jetM = new vector<float>(10); CopyBranch(inputTree, outputTree, "jetM", "jetM", &jetM, "vector<float>");
        //float mll; outputTree->Branch("mll", &mll, "mll/F");
        //int is_OS = 1; outputTree->Branch("is_OS", &is_OS, "is_OS/I");
        //vector<float>* lep_pT = new vector<float>(10); outputTree->Branch("lepPt", "vector<float>", &lep_pT);
        //vector<float>* lep_eta = new vector<float>(10); outputTree->Branch("lepEta", "vector<float>", &lep_eta);
        //vector<float>* lep_phi = new vector<float>(10); outputTree->Branch("lepPhi", "vector<float>", &lep_phi);
        //vector<int>* lep_charge = new vector<int>(10); outputTree->Branch("lepCharge", "vector<int>", &lep_charge);
        //vector<float>* lep_m = new vector<float>(10); outputTree->Branch("lepM", "vector<float>", &lep_m);
        //vector<int>* lepIsoFCTight = new vector<int>{1,1}; outputTree->Branch("lepIsoFCTight", "vector<int>", &lepIsoFCTight);
        //vector<int>* lepIsPR = new vector<int>{1,1}; outputTree->Branch("lepIsPR", "vector<int>", &lepIsPR);
        //Int_t lepChannel; outputTree->Branch("channel", &lepChannel, "channel/I");
        //vector<int>* lep_flavor = new vector<int>(10); outputTree->Branch("lepFlavor", "vector<int>", &lep_flavor);

        ////-----------------------------
        //// set global lepton info
        ////-----------------------------

        //int flavor;
        //if (TString(options.channel).EqualTo("ee")) {
            //flavor = 1;
            //lepChannel = 1;
        //}
        //else if (TString(options.channel).EqualTo("mm")) {
            //flavor = 2;
            //lepChannel = 0;
        //}

        //lep_flavor->clear();
        //lep_flavor->push_back(flavor);
        //lep_flavor->push_back(flavor);

        //lep_m->clear();
        //if (lepChannel == 0) {
            //lep_m->push_back(0.1056583);
            //lep_m->push_back(0.1056583);
        //}
        //else if (lepChannel == 1) {
            //lep_m->push_back(0.0005109);
            //lep_m->push_back(0.0005109);
        //}

        ////-----------------------------
        //// loop over events
        ////-----------------------------

        //Long64_t nentries = inputTree->GetEntries();
        //for (Long64_t i=0; i<nentries; i++) {
            //if (fmod(i,1e5)==0) cout << i << " events processed.\r" << flush;
            //inputTree->GetEntry(i);

            ////--- splitting
            //gamma_pt_split = gamma_pt;
            //gamma_phi_split = gamma_phi;

            //METt_split = METt;
            //try {
                //auto [METl_split_return, mll_return] = this->splitMETlAndMll(METl, gamma_pt);
                //METl_split = METl_split_return;
                //mll = mll_return;
            //}
            //catch(...) {
                //cout << PRED("Problem in event ") << i << endl;
                //continue;
            //}
            //if (mll==0) continue;

            //MET_split = sqrt(pow(METl_split, 2) + pow(METt_split, 2));
            //TLorentzVector MET_split_4vec;
            //MET_split_4vec.SetPtEtaPhiM(MET_split, 0, MET_phi, 0);
            //DPhi_METZPhoton_split = gamma_phi_split - MET_split_4vec.Phi();

            ////--- lepton splitting
            //TLorentzVector z_4vec;
            //z_4vec.SetPtEtaPhiM(gamma_pt_split, gamma_eta, gamma_phi_split, mll);

            //auto [good_event, l0_lab_4vec, l1_lab_4vec] = splitLeptons(z_4vec);
            //if (!good_event) continue;

            //lep_pT->clear();
            //lep_pT->push_back(l0_lab_4vec.Pt());
            //lep_pT->push_back(l1_lab_4vec.Pt());
            //lep_eta->clear();
            //lep_eta->push_back(l0_lab_4vec.Eta());
            //lep_eta->push_back(l1_lab_4vec.Eta());
            //lep_phi->clear();
            //lep_phi->push_back(l0_lab_4vec.Phi());
            //lep_phi->push_back(l1_lab_4vec.Phi());
            //int charge = myRandom.Integer(2)*2-1; // random leading lepton charge
            //lep_charge->clear();
            //lep_charge->push_back(charge);
            //lep_charge->push_back(-charge);

            //lep_n = 2;
            //DPhi_METLepLeading_split = fabs(MET_split_4vec.DeltaPhi(l0_lab_4vec));
            //DPhi_METLepSecond_split = fabs(MET_split_4vec.DeltaPhi(l1_lab_4vec));
            //DR_2Lep = l0_lab_4vec.DeltaR(l1_lab_4vec);

            //outputTree->Fill();
        //}
        //cout << endl;

        ////-----------------------------
        //// write
        ////-----------------------------

        //this->output_file->cd();
        //outputTree->Write();

        //cout << PBLU("Done with splitting") << endl;
        //cout << endl;
        //delete this->output_file;
    //}
};

//------------
// UNIT TESTS
//------------

Options setSplittingUnitTestOptions(Options options) {
    options.reduction_folder = options.unit_test_folder + "ReducedNtuples/";
    options.splitting_folder = "./";

    options.period = "data15-16";
    options.mc_period = getMCPeriod(options.period);
    options.data_period = DataPeriod(options.period);
    options.channel = "mm";

    return options;
}

void performSplittingUnitTests(Options options) {
    cout << BOLD(PBLU("Performing unit testing on splitting step")) << endl;
    cout << endl;

    //--- get Z MC feature plots
    PhotonSplitter splitter(options);
    splitter.fillZMCHistograms(options, true);

    //--- random sampling test
    for (int i=0; i<10; i++) {
        float mll;
        if (splitter.getRandomSample(splitter.z_mll_histptrs, 38, 80, mll))
            cout << mll << endl;
    }
    //vector<vector<ROOT::RDF::RResultPtr<TH1D>>> z_mll_histptrs;
    //vector<vector<ROOT::RDF::RResultPtr<TH1D>>> z_cmleptheta_histptrs;

    //--- set up photon feature plots for comparison
    //map<string, map<int, TH1F*>> photon_plots;
    //for (int i=0; i<bins::n_pt_bins+2; i++) {
        //photon_plots["lepcm_theta"][i] = new TH1F("", "", 100, 0, 3);
        //photon_plots["lepEta"][i] = new TH1F("", "", 100, -3, 3);
        //photon_plots["METl_raw"][i] = new TH1F("", "", 100, -100, 100);
        //photon_plots["METl"][i] = new TH1F("", "", 100, -100, 100);
        //photon_plots["mll"][i] = new TH1F("", "", 100, 0, 150);
    //}

    ////--- open photon files
    //string in_file_name;
    //if (options.is_data) in_file_name = options.reduction_folder + options.data_period + "_data_photon.root";
    //else in_file_name = options.reduction_folder + options.mc_period + "_SinglePhoton222.root";
    //TFile *photon_file = new TFile(in_file_name.c_str());
    //TTree *photon_tree = (TTree*)photon_file->Get("BaselineTree");

    //float gamma_pt; photon_tree->SetBranchAddress("gamma_pt", &gamma_pt);
    //float gamma_eta; photon_tree->SetBranchAddress("gamma_eta", &gamma_eta);
    //float gamma_phi; photon_tree->SetBranchAddress("gamma_phi", &gamma_phi);
    //double totalWeight; photon_tree->SetBranchAddress("totalWeight", &totalWeight);

    ////--- perform photon splitting
    //float Z_m = 91.1876;
    //for (int i=0; i<photon_tree->GetEntries(); i++) {
        //photon_tree->GetEntry(i);
        //int pt_bin = bins::hist_pt_bins->FindBin(gamma_pt);

        //photon_plots["lepcm_theta"][pt_bin]->Fill(splitter.getRandomLepTheta());

        //TLorentzVector z_4vec;
        //z_4vec.SetPtEtaPhiM(gamma_pt, gamma_eta, gamma_phi, Z_m);
        //auto [good_event, l0_lab_4vec, l1_lab_4vec] = splitter.splitLeptons(z_4vec);
        //if (good_event) {
            //photon_plots["lepEta"][pt_bin]->Fill(l0_lab_4vec.Eta(), totalWeight);
            //photon_plots["lepEta"][pt_bin]->Fill(l0_lab_4vec.Eta(), totalWeight);
        //}

        //auto [METl, mll] = splitter.splitMETlAndMll(METl_unsplit, gamma_pt);
        //photon_plots["METl_raw"][pt_bin]->Fill(METl_unsplit, totalWeight);
        //photon_plots["METl"][pt_bin]->Fill(METl, totalWeight);
        //photon_plots["mll"][pt_bin]->Fill(mll, totalWeight);
    //}

    //--- draw and save comparison plots
    TCanvas *can = new TCanvas("can","can",600,600);

    //for (auto& [key, hist] : zmc_plots) {
        //THStack *zmc_stack = new THStack("zmc_stack","");
        //THStack *photon_stack = new THStack("gmc_stack","");;

        //for (int i=0; i<bins::n_pt_bins+2; i++) {
            //zmc_plots[key][i]->Draw("hist");
            //float z_yield = zmc_plots[key][i]->Integral(0, zmc_plots[key][i]->GetNbinsX()+1);
            //float g_yield = photon_plots[key][i]->Integral(0, photon_plots[key][i]->GetNbinsX()+1);
            ////cout << z_yield << " " << g_yield << endl;
            //float scale = g_yield != 0 ? z_yield/g_yield : 0;

            //if (key == "METl") {
                //float g_yield_raw = photon_plots[key][i]->Integral(0, photon_plots[key][i]->GetNbinsX()+1);
                //float scale_raw = g_yield_raw != 0 ? z_yield/g_yield_raw : 0;

                //photon_plots["METl_raw"][i]->Scale(scale_raw);
                //photon_plots["METl_raw"][i]->SetLineColor(kBlue);
                //photon_plots["METl_raw"][i]->Draw("hist same");
            //}

            //photon_plots[key][i]->Scale(scale);
            //photon_plots[key][i]->SetLineColor(kRed);
            //photon_plots[key][i]->Draw("hist same");

            //TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
            //leg->AddEntry(zmc_plots[key][i], "Zjets", "f");
            //if (key == "METl")
                //leg->AddEntry(photon_plots["METl_raw"][i], "photon raw", "f");
            //leg->AddEntry(photon_plots[key][i], "photon corrected", "f");
            //leg->SetBorderSize(0);
            //leg->SetFillColor(0);
            //leg->Draw();

            //can->Print("Diagnostics/Splitting/" + key + "_pt_bin_" + i + ".eps");

            //if (z_yield > 0) zmc_stack->Add(zmc_plots[key][i]);
            //if (g_yield > 0) photon_stack->Add(photon_plots[key][i]);
        //}

        //zmc_stack->GetStack()->Last()->Draw("hist");
        //photon_stack->GetStack()->Last()->Draw("hist same");
        ////photon_stack->GetStack()->Last()->Draw("hist");
        ////photon_stack->SetLineColor(kRed);

        //TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
        //leg->AddEntry(zmc_plots[key][0], "Zjets", "f");
        //leg->AddEntry(photon_plots[key][0], "photon corrected", "f");
        //leg->SetBorderSize(0);
        //leg->SetFillColor(0);
        //leg->Draw();

        //can->Print(("Diagnostics/Splitting/" + key + ".eps").c_str());

        //passTest("Comparison plots for " + key + " produced");
    //}

    delete can;

    passTest("Passed all unit tests");
    cout << endl;
}

//---------------
// MAIN FUNCTION
//---------------

void SplitPhotons(Options options) {
    if (options.unit_testing) {
        options = setSplittingUnitTestOptions(options);
        performSplittingUnitTests(options);
    }
    //else {
        //cout << BOLD(PBLU("Performing splitting")) << endl;
        //cout << endl;

        ////--- split photons
        //options.run_vgamma = false;
        //PhotonSplitter photon_splitter(options);
        //photon_splitter.splitPhotons(options);
        //if (options.diagnostic_plots) photon_splitter.drawMETlDistributions();

        ////--- split Vgamma - this MC component must be subtracted from photon data
        //if (options.is_data) {
            //options.run_vgamma = true;
            //PhotonSplitter vgamma_splitter(options);
            //vgamma_splitter.splitPhotons(options);
            //if (options.diagnostic_plots) vgamma_splitter.drawMETlDistributions();
        //}
    //}
}

#include "../Settings.C"
#include "../CommonFunctions/CommonLibraries.C"
#include "../CommonFunctions/CommonFunctions.C"
#include "../PhotonSmearing/GetDijetVariables.C"

using namespace std;

void RebinHistogram(TH1D* hist) {
    float negative_yield = 0.;
    float positive_yield = 0.;
    for (int bin=1;bin<=hist->GetNbinsX();bin++) {
        if (hist->GetBinContent(bin)>=0) positive_yield += hist->GetBinContent(bin);
        else negative_yield += hist->GetBinContent(bin);
    }
    while (abs(negative_yield/positive_yield)>0.05) {
        hist->Rebin(2);
        negative_yield = 0.;
        positive_yield = 0.;
        for (int bin=1;bin<=hist->GetNbinsX();bin++) {
            if (hist->GetBinContent(bin)>=0) positive_yield += hist->GetBinContent(bin);
            else negative_yield += hist->GetBinContent(bin);
        }
    }
}

void GetBaseLineEvents(string sampleID, string outputName, string pathToNtuples, string isData, string treename = "outputTree" ) {

    //---------------------------------------------
    // open input and output files, get TTrees
    //---------------------------------------------

    TH1::SetDefaultSumw2();

    string filename = Form("%s%s.root",pathToNtuples.c_str(),sampleID.c_str()); 
    TFile* inputFile = TFile::Open(filename.c_str());
    TTree* inputTree = (TTree*)inputFile->Get(treename.c_str());

    Float_t _nGenEvents = 1.;
    TH1D* hist_EventCount = new TH1D("hist_EventCount","",3,0,3);
    if (isData == "MC") {
        cout << "Setting _nGenEvents = 1 for now NEED TO FIX" << endl;
        _nGenEvents    = 1.0;
        hist_EventCount->SetBinContent(1,1.0);
    }

    std::cout << inputTree << std::endl;
    cout << endl;
    cout << "Opening file           : " << filename        << endl;
    cout << "Events in ntuple       : " << inputTree->GetEntries() << endl;
    if (isData == "MC") {
        cout << "Total generated events : " << _nGenEvents     << endl;
    }

    string outfilename = ntuple_path + "/" + outputName + "/" + sampleID.c_str() + ".root";
    cout << "Writing to : " << outfilename << endl;
    TFile outputFile( outfilename.c_str() , "recreate" );
    TTree* BaselineTree = new TTree("BaselineTree","baseline tree");

    //-----------------------------
    // access and copy over existing branches
    //-----------------------------

    inputTree->SetBranchStatus("*", 0);

    Double_t genWeight;
    Double_t eventWeight;
    Double_t leptonWeight;
    Double_t jvtWeight;
    Double_t bTagWeight;
    Double_t pileupWeight;
    Double_t FFWeight;
    if (isData == "MC") {
        SetInputBranch(inputTree, "genWeight", &genWeight);
        SetInputBranch(inputTree, "eventWeight", &eventWeight);
        SetInputBranch(inputTree, "leptonWeight", &leptonWeight);
        SetInputBranch(inputTree, "jvtWeight", &jvtWeight);
        SetInputBranch(inputTree, "bTagWeight", &bTagWeight);
        SetInputBranch(inputTree, "pileupWeight", &pileupWeight);
        SetInputBranch(inputTree, "FFWeight", &FFWeight);
    }

    Int_t RunNumber; CopyBranch(inputTree, BaselineTree, "RunNumber", "RunNumber", &RunNumber, "I");
    bool trigMatch_1L2LTrig; SetInputBranch(inputTree, "trigMatch_1L2LTrig", &trigMatch_1L2LTrig);
    Int_t nLep_signal; SetInputBranch(inputTree, "nLep_signal", &nLep_signal);
    Int_t nLep_base; SetInputBranch(inputTree, "nLep_base", &nLep_base);
    float MET; CopyBranch(inputTree, BaselineTree, "met_Et", "MET", &MET, "F");
    float MET_phi; CopyBranch(inputTree, BaselineTree, "met_Phi", "MET_Phi", &MET_phi, "F");
    float MET_softTerm; CopyBranch(inputTree, BaselineTree, "TST_Et", "MET_softTerm", &MET_softTerm, "F");
    float MET_softPhi; CopyBranch(inputTree, BaselineTree, "TST_Phi", "MET_softPhi", &MET_softPhi, "F");
    float HT; CopyBranch(inputTree, BaselineTree, "Ht30", "HT", &HT, "F");
    Float_t mll; CopyBranch(inputTree, BaselineTree, "mll", "mll", &mll, "F");
    float mjj; CopyBranch(inputTree, BaselineTree, "mjj", "mjj", &mjj, "F");
    Int_t jet_n; CopyBranch(inputTree, BaselineTree, "nJet30", "jet_n", &jet_n, "I");
    Int_t bjet_n; CopyBranch(inputTree, BaselineTree, "nBJet30_MV2c10_FixedCutBEff_77", "bjet_n", &bjet_n, "I");
    int nBJet20_MV2c10_FixedCutBEff_77; CopyBranch(inputTree, BaselineTree, "nBJet20_MV2c10_FixedCutBEff_77", "nBJet20_MV2c10_FixedCutBEff_77", &nBJet20_MV2c10_FixedCutBEff_77, "I");
    Float_t Z_pt; CopyBranch(inputTree, BaselineTree, "Ptll", "Z_pt", &Z_pt, "F");

    std::vector<int>* lepFlavor = new std::vector<int>(10); CopyBranch(inputTree, BaselineTree, "lepFlavor", "lepFlavor", &lepFlavor, "std::vector<int>");
    std::vector<int>* lepCharge = new std::vector<int>(10); CopyBranch(inputTree, BaselineTree, "lepCharge", "lepCharge", &lepCharge, "std::vector<int>");
    std::vector<int>* lepIsoFCTight = new std::vector<int>(10); CopyBranch(inputTree, BaselineTree, "lepIsoFCTight", "lepIsoFCTight", &lepCharge, "std::vector<int>");
    std::vector<float>* lep_pT = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "lepPt", "lep_pT", &lep_pT, "std::vector<float>");
    std::vector<float>* lep_eta = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "lepEta", "lep_eta", &lep_eta, "std::vector<float>");
    std::vector<float>* lep_phi = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "lepPhi", "lep_phi", &lep_phi, "std::vector<float>");
    std::vector<float>* jet_pT = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jetPt", "jet_pT", &jet_pT, "std::vector<float>");
    std::vector<float>* jet_eta = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jetEta", "jet_eta", &jet_eta, "std::vector<float>");
    std::vector<float>* jet_phi = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jetPhi", "jet_phi", &jet_phi, "std::vector<float>");
    std::vector<float>* jet_m = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jetM", "jet_m", &jet_m, "std::vector<float>");

    //-----------------------------
    // add new branches and histograms
    //-----------------------------

    float METl; BaselineTree->Branch("METl",&METl,"METl/F");
    float METt; BaselineTree->Branch("METt",&METt,"METt/F");
    float DPhi_2Lep; BaselineTree->Branch("DPhi_2Lep",&DPhi_2Lep,"DPhi_2Lep/F");
    float DR_2Lep; BaselineTree->Branch("DR_2Lep",&DR_2Lep,"DR_2Lep/F");
    float DPhi_METPhoton; BaselineTree->Branch("DPhi_METPhoton",&DPhi_METPhoton,"DPhi_METPhoton/F");
    float DPhi_METLepLeading; BaselineTree->Branch("DPhi_METLepLeading",&DPhi_METLepLeading,"DPhi_METLepLeading/F");
    float DPhi_METLepSecond; BaselineTree->Branch("DPhi_METLepSecond",&DPhi_METLepSecond,"DPhi_METLepSecond/F");
    float DPhi_METLepMin; BaselineTree->Branch("DPhi_METLepMin",&DPhi_METLepMin,"DPhi_METLepMin/F");
    float MinDPhi_PhotonJet; BaselineTree->Branch("MinDPhi_PhotonJet",&MinDPhi_PhotonJet,"MinDPhi_PhotonJet/F");
    int channel; BaselineTree->Branch("channel",&channel,"channel/I");
    int is_OS; BaselineTree->Branch("is_OS",&is_OS,"is_OS/I");
    float Z_eta; BaselineTree->Branch("Z_eta",&Z_eta,"Z_eta/F");
    float Z_phi; BaselineTree->Branch("Z_phi",&Z_phi,"Z_phi/F");
    double totalWeight; BaselineTree->Branch("totalWeight",&totalWeight,"totalWeight/D");

    TH1D* hist_cutflow_raw = new TH1D("hist_cutflow_raw","",8,0,8);
    hist_cutflow_raw->SetStats(0);
    hist_cutflow_raw->GetXaxis()->SetBinLabel(1, "2lep");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(2, "flavor");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(3, "trigger");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(4, "OS");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(5, "lep pT");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(6, "njet");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(7, "mll");
    hist_cutflow_raw->GetXaxis()->SetBinLabel(8, "prompt");

    TH1D* hist_cutflow_weight = new TH1D("hist_cutflow_weight","",8,0,8);
    hist_cutflow_weight->SetStats(0);
    hist_cutflow_weight->GetXaxis()->SetBinLabel(1, "2lep");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(2, "flavor");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(3, "trigger");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(4, "OS");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(5, "lep pT");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(6, "njet");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(7, "mll");
    hist_cutflow_weight->GetXaxis()->SetBinLabel(8, "prompt");

    TH1D* hist_METl_Pt[bin_size];
    TH1D* hist_METt_Pt[bin_size];
    for (int bin=0;bin<bin_size;bin++) {
        hist_METl_Pt[bin] = new TH1D(TString("hist_METl_Pt_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        hist_METl_Pt[bin]->SetStats(0);
        hist_METt_Pt[bin] = new TH1D(TString("hist_METt_Pt_")+TString::Itoa(bin,10),"",40000,-30000,10000);
        hist_METt_Pt[bin]->SetStats(0);
    }

    TH1D* hist_Mll_dPt[dpt_bin_size];

    TH1D* hist_low_njet = new TH1D("hist_low_njet","",bin_size,njet_bin);
    TH1D* hist_low_nbjet = new TH1D("hist_low_nbjet","",bin_size,njet_bin);
    TH1D* hist_low_pt = new TH1D("hist_low_pt","",bin_size,pt_bin);
    TH1D* hist_sm_pt = new TH1D("hist_sm_pt","",bin_size,sm_pt_bin);
    TH1D* hist_low_ht = new TH1D("hist_low_ht","",bin_size,ht_bin);

    //-----------------------------
    // loop over events
    //-----------------------------

    Long64_t nentries = inputTree->GetEntries();
    //nentries = 1000;

    float N_passMET100 = 0.;

    for (Long64_t i=0;i<nentries;i+=event_interval) {

        if (fmod(i,1e5)==0) std::cout << i << " events processed." << std::endl;
        inputTree->GetEntry(i);

        if (MET>100) N_passMET100 += 1; 

        if ( nLep_signal  != 2                 ) continue; // exactly 2 signal leptons
        if ( nLep_base    != 2                 ) continue; // exactly 2 baseline leptons
        if ( lep_pT->at(0) < leading_lep_pt_cut ) continue; // 1st lep pT > 25 GeV
        if ( lep_pT->at(1) < second_lep_pt_cut  ) continue; // 2nd lep pT > 25 GeV

        //--- evaluate weight
        totalWeight = 1;
        if (isData == "MC") totalWeight = genWeight * eventWeight * leptonWeight * jvtWeight * bTagWeight * pileupWeight * FFWeight;

        //--- determine channel
        channel = -1;
        if( ( lepFlavor->at(0) == 1 && lepFlavor->at(1) == 1 ) && trigMatch_1L2LTrig  ) channel = 1; // ee
        if( ( lepFlavor->at(0) == 2 && lepFlavor->at(1) == 2 ) && trigMatch_1L2LTrig  ) channel = 0; // mumu
        if( ( lepFlavor->at(0) == 1 && lepFlavor->at(1) == 2 ) && trigMatch_1L2LTrig  ) channel = 2; // em
        if( ( lepFlavor->at(0) == 2 && lepFlavor->at(1) == 1 ) && trigMatch_1L2LTrig  ) channel = 3; // me

        //--- determine OS / SS
        is_OS = -1;
        if( lepCharge->at(0) != lepCharge->at(1) ) is_OS = 1;
        if( lepCharge->at(0) == lepCharge->at(1) ) is_OS = 0;

        if( channel < 0 ) continue; // require exactly 2 signal leptons and corresponding triggers
        if( is_OS != 1  ) continue; // require opposite-sign
        if( jet_n < 1   ) continue; // require at least 1 pT > 30 GeV jets

        //---------------------------------------------
        // here we compute the MET parallel and perpendicular components
        // and DR between photon and nearby jet
        // and Oslo's MET_rel
        //---------------------------------------------

        TLorentzVector lep0vec;
        TLorentzVector lep1vec;

        lep0vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
        lep1vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);

        Z_eta = ( lep0vec + lep1vec ).Eta();
        Z_phi = ( lep0vec + lep1vec ).Phi();

        METt = MET*TMath::Sin(MET_phi-Z_phi);
        METl = MET*TMath::Cos(MET_phi-Z_phi);
        TLorentzVector met_4vec;
        met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
        TLorentzVector z_4vec;
        z_4vec.SetPtEtaPhiM(Z_pt,Z_eta,Z_phi,0);
        DPhi_METPhoton = fabs(met_4vec.DeltaPhi(z_4vec));
        TLorentzVector lep0_4vec;
        lep0_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
        TLorentzVector lep1_4vec;
        lep1_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);
        DPhi_2Lep = fabs(lep0_4vec.DeltaPhi(lep1_4vec));
        DR_2Lep = lep0_4vec.DeltaR(lep1_4vec);
        DPhi_METLepLeading = fabs(met_4vec.DeltaPhi(lep0_4vec));
        DPhi_METLepSecond = fabs(met_4vec.DeltaPhi(lep1_4vec));
        DPhi_METLepMin = min(DPhi_METLepLeading,DPhi_METLepSecond);
        TLorentzVector tst_4vec;
        tst_4vec.SetPtEtaPhiM(MET_softTerm,0,MET_softPhi,0);
        float DPhi_TSTLepLeading = fabs(tst_4vec.DeltaPhi(lep0_4vec));
        float DPhi_TSTLepSecond = fabs(tst_4vec.DeltaPhi(lep1_4vec));
        float DPhi_TSTLepMin = min(DPhi_TSTLepLeading,DPhi_TSTLepSecond);
        float MinDR_Lep0Jet = 1000.;
        float MinDR_Lep1Jet = 1000.;
        float MinDR_PhotonJet = 1000.;
        MinDPhi_PhotonJet = 1000.;
        TLorentzVector jet_4vec;
        for (unsigned int j=0;j<jet_pT->size();j++) {
            jet_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),jet_m->at(j));
            float DR_Lep0Jet = jet_4vec.DeltaR(lep0_4vec);
            float DR_Lep1Jet = jet_4vec.DeltaR(lep1_4vec);
            if (MinDR_Lep0Jet>DR_Lep0Jet) MinDR_Lep0Jet = DR_Lep0Jet;
            if (MinDR_Lep1Jet>DR_Lep1Jet) MinDR_Lep1Jet = DR_Lep1Jet;
            float DR_PhotonJet = jet_4vec.DeltaR(z_4vec);
            float DPhi_PhotonJet = jet_4vec.DeltaPhi(z_4vec);
            if (MinDR_PhotonJet>DR_PhotonJet) MinDR_PhotonJet = DR_PhotonJet;
            if (MinDPhi_PhotonJet>DPhi_PhotonJet) MinDPhi_PhotonJet = DPhi_PhotonJet;
        }
        float min_DPhi_MET_LepJet = 1000.;
        float DPhi_MET_LepJet = 1000.;
        for (unsigned int j=0;j<jet_pT->size();j++) {
            jet_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),jet_m->at(j));
            DPhi_MET_LepJet = jet_4vec.DeltaR(met_4vec);
            if (min_DPhi_MET_LepJet>DPhi_MET_LepJet) min_DPhi_MET_LepJet = DPhi_MET_LepJet;
        }
        DPhi_MET_LepJet = lep0_4vec.DeltaR(met_4vec);
        if (min_DPhi_MET_LepJet>DPhi_MET_LepJet) min_DPhi_MET_LepJet = DPhi_MET_LepJet;
        DPhi_MET_LepJet = lep1_4vec.DeltaR(met_4vec);
        if (min_DPhi_MET_LepJet>DPhi_MET_LepJet) min_DPhi_MET_LepJet = DPhi_MET_LepJet;
        float MET_rel = MET;
        if (min_DPhi_MET_LepJet<TMath::Pi()/2.) MET_rel = MET*TMath::Sin(min_DPhi_MET_LepJet);

        //---------------------------------------------
        // compute dijet system variables, m80jj, W pT, DR(2jet), etc.
        //---------------------------------------------
        z_4vec.SetPtEtaPhiM(Z_pt,Z_eta,Z_phi,0);
        met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
        GetDijetVariables(z_4vec, met_4vec, jet_pT, jet_eta, jet_phi, jet_m);

        BaselineTree->Fill();     
    }

    std::cout << "write output..." << std::endl;
    BaselineTree->Write();

    hist_cutflow_raw->Write();
    hist_cutflow_weight->Write();
    for (int bin=0;bin<bin_size;bin++) {
        hist_METl_Pt[bin]->Write();
        hist_METt_Pt[bin]->Write();
    }
    if (isData == "MC") {
        hist_EventCount->SetBinContent(2,nentries);
        hist_EventCount->SetBinContent(3,N_passMET100);
        hist_EventCount->Write();
    }

    std::cout << "done." << std::endl;
    outputFile.Close();
    delete inputFile;
}

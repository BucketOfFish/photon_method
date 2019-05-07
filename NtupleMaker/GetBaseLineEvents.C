#include "../Settings.C"
#include "../CommonFunctions/CommonLibraries.C"
#include "../CommonFunctions/CommonFunctions.C"

using namespace std;

void GetBaseLineEvents(string sampleID, string outputName, string pathToNtuples, string isData, string treeName = "outputTree" ) {

    //---------------------------------------------
    // open input and output files, get TTrees
    //---------------------------------------------

    TH1::SetDefaultSumw2();

    string filename = Form("%s%s.root",pathToNtuples.c_str(), sampleID.c_str()); 
    TFile* inputFile = TFile::Open(filename.c_str());
    TTree* inputTree = (TTree*)inputFile->Get(treeName.c_str());

    float lumi = GetLumi(pathToNtuples);

    cout << endl;
    cout << "Opening file           : " << filename        << endl;
    cout << "Events in ntuple       : " << inputTree->GetEntries() << endl;
	cout << "using luminosity       : " << lumi          << endl;

    string outfilename = ntuple_path + "/" + outputName + "/" + sampleID.c_str() + ".root";
    cout << "Writing to : " << outfilename << endl;
    TFile outputFile( outfilename.c_str(), "recreate" );
    TTree* BaselineTree = new TTree("BaselineTree", "baseline tree");

    //-----------------------------
    // access, copy, and create branches
    //-----------------------------

    inputTree->SetBranchStatus("*", 0);

    //--- triggers and weights
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
    bool trigMatch_1L2LTrig; SetInputBranch(inputTree, "trigMatch_1L2LTrig", &trigMatch_1L2LTrig);
    double totalWeight; BaselineTree->Branch("totalWeight",&totalWeight,"totalWeight/D");

    //--- event selection
    Int_t nLep_signal; SetInputBranch(inputTree, "nLep_signal", &nLep_signal);
    Int_t nLep_base; SetInputBranch(inputTree, "nLep_base", &nLep_base);
    Int_t jet_n; CopyBranch(inputTree, BaselineTree, "nJet30", "jet_n", &jet_n, "I");
    int channel; BaselineTree->Branch("channel",&channel,"channel/I");
    int is_OS; BaselineTree->Branch("is_OS",&is_OS,"is_OS/I");
    std::vector<int>* lepFlavor = new std::vector<int>(10); CopyBranch(inputTree, BaselineTree, "lepFlavor", "lepFlavor", &lepFlavor, "std::vector<int>");
    std::vector<int>* lepCharge = new std::vector<int>(10); CopyBranch(inputTree, BaselineTree, "lepCharge", "lepCharge", &lepCharge, "std::vector<int>");
    std::vector<int>* lepIsoFCTight = new std::vector<int>(10); CopyBranch(inputTree, BaselineTree, "lepIsoFCTight", "lepIsoFCTight", &lepCharge, "std::vector<int>");
    std::vector<float>* lep_pT = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "lepPt", "lep_pT", &lep_pT, "std::vector<float>");
    std::vector<float>* lep_eta = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "lepEta", "lep_eta", &lep_eta, "std::vector<float>");
    std::vector<float>* lep_phi = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "lepPhi", "lep_phi", &lep_phi, "std::vector<float>");

    //--- MET components, and DR and DPhi between objects
    float MET; CopyBranch(inputTree, BaselineTree, "met_Et", "MET", &MET, "F");
    float MET_phi; SetInputBranch(inputTree, "met_Phi", &MET_phi);
    Float_t Z_pt; CopyBranch(inputTree, BaselineTree, "Ptll", "Z_pt", &Z_pt, "F");
    float Z_eta; BaselineTree->Branch("Z_eta",&Z_eta,"Z_eta/F");
    float Z_phi; BaselineTree->Branch("Z_phi",&Z_phi,"Z_phi/F");
    float METl; BaselineTree->Branch("METl",&METl,"METl/F");
    float METt; BaselineTree->Branch("METt",&METt,"METt/F");
    float DR_2Lep; BaselineTree->Branch("DR_2Lep",&DR_2Lep,"DR_2Lep/F");
    float DPhi_METPhoton; BaselineTree->Branch("DPhi_METPhoton",&DPhi_METPhoton,"DPhi_METPhoton/F");
    float DPhi_2Lep; BaselineTree->Branch("DPhi_2Lep",&DPhi_2Lep,"DPhi_2Lep/F");
    float DPhi_METLepLeading; BaselineTree->Branch("DPhi_METLepLeading",&DPhi_METLepLeading,"DPhi_METLepLeading/F");
    float DPhi_METLepSecond; BaselineTree->Branch("DPhi_METLepSecond",&DPhi_METLepSecond,"DPhi_METLepSecond/F");
    float DPhi_METLepMin; BaselineTree->Branch("DPhi_METLepMin",&DPhi_METLepMin,"DPhi_METLepMin/F");

    //--- not used here?
    Int_t RunNumber; CopyBranch(inputTree, BaselineTree, "RunNumber", "RunNumber", &RunNumber, "I");
    float mjj; CopyBranch(inputTree, BaselineTree, "mjj", "mjj", &mjj, "F");
    Int_t bjet_n; CopyBranch(inputTree, BaselineTree, "nBJet30_MV2c10_FixedCutBEff_77", "bjet_n", &bjet_n, "I");
    float HT; CopyBranch(inputTree, BaselineTree, "Ht30", "HT", &HT, "F");
    Float_t mll; CopyBranch(inputTree, BaselineTree, "mll", "mll", &mll, "F");
    int nBJet20_MV2c10_FixedCutBEff_77; CopyBranch(inputTree, BaselineTree, "nBJet20_MV2c10_FixedCutBEff_77", "nBJet20_MV2c10_FixedCutBEff_77", &nBJet20_MV2c10_FixedCutBEff_77, "I");
    std::vector<float>* jet_pT = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jetPt", "jet_pT", &jet_pT, "std::vector<float>");
    std::vector<float>* jet_eta = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jetEta", "jet_eta", &jet_eta, "std::vector<float>");
    std::vector<float>* jet_phi = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jetPhi", "jet_phi", &jet_phi, "std::vector<float>");
    std::vector<float>* jet_m = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "jetM", "jet_m", &jet_m, "std::vector<float>");

    //-----------------------------
    // loop over events
    //-----------------------------

    Long64_t nentries = inputTree->GetEntries();

    for (Long64_t i=0;i<nentries;i+=event_interval) {

        if (fmod(i,1e5)==0) std::cout << i << " events processed." << std::endl;
        inputTree->GetEntry(i);

        //--- event selection
        if ( nLep_signal  != 2                 ) continue; // exactly 2 signal leptons
        if ( nLep_base    != 2                 ) continue; // exactly 2 baseline leptons
        if ( lep_pT->at(0) < leading_lep_pt_cut ) continue; // 1st lep pT > 25 GeV
        if ( lep_pT->at(1) < second_lep_pt_cut  ) continue; // 2nd lep pT > 25 GeV
        if ( jet_n < 1   ) continue; // require at least 1 pT > 30 GeV jets
        if ( !trigMatch_1L2LTrig ) continue; // need 2 lepton trigger

        //--- evaluate weight
        totalWeight = 1;
        if (isData == "MC") totalWeight = genWeight * eventWeight * leptonWeight * jvtWeight * bTagWeight * pileupWeight * FFWeight;

        //--- determine channel
        channel = -1;
        if ( lepFlavor->at(0) == 2 && lepFlavor->at(1) == 2 ) channel = 0; // mumu
        if ( lepFlavor->at(0) == 1 && lepFlavor->at(1) == 1 ) channel = 1; // ee
        if ( lepFlavor->at(0) == 1 && lepFlavor->at(1) == 2 ) channel = 2; // em
        if ( lepFlavor->at(0) == 2 && lepFlavor->at(1) == 1 ) channel = 3; // me

        if ( channel < 0 ) continue; // require exactly 2 signal leptons

        //--- determine OS / SS
        is_OS = -1;
        if ( lepCharge->at(0) != lepCharge->at(1) ) is_OS = 1;
        if ( lepCharge->at(0) == lepCharge->at(1) ) is_OS = 0;

        if ( is_OS != 1 ) continue; // require opposite-sign

        //--- compute 4-vectors of objects
        TLorentzVector lep0_4vec;
        TLorentzVector lep1_4vec;
        lep0_4vec.SetPtEtaPhiM(lep_pT->at(0),lep_eta->at(0),lep_phi->at(0),0);
        lep1_4vec.SetPtEtaPhiM(lep_pT->at(1),lep_eta->at(1),lep_phi->at(1),0);

        TLorentzVector z_4vec;
        Z_eta = (lep0vec+lep1vec).Eta();
        Z_phi = (lep0vec+lep1vec).Phi();
        z_4vec.SetPtEtaPhiM(Z_pt,Z_eta,Z_phi,0);

        TLorentzVector met_4vec;
        met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);

        //--- compute MET parallel and perpendicular components
        METt = MET*TMath::Sin(MET_phi-Z_phi);
        METl = MET*TMath::Cos(MET_phi-Z_phi);

        //--- compute DR and DPhi between objects
        DR_2Lep = lep0_4vec.DeltaR(lep1_4vec);
        DPhi_METPhoton = fabs(met_4vec.DeltaPhi(z_4vec));
        DPhi_2Lep = fabs(lep0_4vec.DeltaPhi(lep1_4vec));
        DPhi_METLepLeading = fabs(met_4vec.DeltaPhi(lep0_4vec));
        DPhi_METLepSecond = fabs(met_4vec.DeltaPhi(lep1_4vec));
        DPhi_METLepMin = min(DPhi_METLepLeading,DPhi_METLepSecond);

        BaselineTree->Fill();     
    }

    std::cout << "write output..." << std::endl;
    BaselineTree->Write();

    std::cout << "done." << std::endl;
    outputFile.Close();
    delete inputFile;
}

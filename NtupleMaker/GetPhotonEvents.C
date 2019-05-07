#include "../Common/Settings.C"
#include "../Common/CommonLibraries.C"
#include "../Common/CommonFunctions.C"

using namespace std;

void GetPhotonEvents(string sampleID, string outputName, string pathToNtuples, string isData, string treeName = "tree_NoSys" ) {

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
	int trigMatch_HLT_g15_loose_L1EM7; SetInputBranch(inputTree, "trigMatch_HLT_g15_loose_L1EM7", &trigMatch_HLT_g15_loose_L1EM7);
	int trigMatch_HLT_g25_loose_L1EM15; SetInputBranch(inputTree, "trigMatch_HLT_g25_loose_L1EM15", &trigMatch_HLT_g25_loose_L1EM15);
	int trigMatch_HLT_g35_loose_L1EM15; SetInputBranch(inputTree, "trigMatch_HLT_g35_loose_L1EM15", &trigMatch_HLT_g35_loose_L1EM15);
	int trigMatch_HLT_g40_loose_L1EM15; SetInputBranch(inputTree, "trigMatch_HLT_g40_loose_L1EM15", &trigMatch_HLT_g40_loose_L1EM15);
	int trigMatch_HLT_g45_loose_L1EM15; SetInputBranch(inputTree, "trigMatch_HLT_g45_loose_L1EM15", &trigMatch_HLT_g45_loose_L1EM15);
	int trigMatch_HLT_g50_loose_L1EM15; SetInputBranch(inputTree, "trigMatch_HLT_g50_loose_L1EM15", &trigMatch_HLT_g50_loose_L1EM15);
	int trigMatch_HLT_g60_loose; SetInputBranch(inputTree, "trigMatch_HLT_g60_loose", &trigMatch_HLT_g60_loose);
	int trigMatch_HLT_g70_loose; SetInputBranch(inputTree, "trigMatch_HLT_g70_loose", &trigMatch_HLT_g70_loose);
	int trigMatch_HLT_g80_loose; SetInputBranch(inputTree, "trigMatch_HLT_g80_loose", &trigMatch_HLT_g80_loose);
	int trigMatch_HLT_g100_loose; SetInputBranch(inputTree, "trigMatch_HLT_g100_loose", &trigMatch_HLT_g100_loose);
	int trigMatch_HLT_g120_loose; SetInputBranch(inputTree, "trigMatch_HLT_g120_loose", &trigMatch_HLT_g120_loose);
	int trigMatch_HLT_g140_loose; SetInputBranch(inputTree, "trigMatch_HLT_g140_loose", &trigMatch_HLT_g140_loose);
	float trigPrescale_HLT_g15_loose_L1EM7; SetInputBranch(inputTree, "trigPrescale_HLT_g15_loose_L1EM7", &trigPrescale_HLT_g15_loose_L1EM7);
	float trigPrescale_HLT_g25_loose_L1EM15; SetInputBranch(inputTree, "trigPrescale_HLT_g25_loose_L1EM15", &trigPrescale_HLT_g25_loose_L1EM15);
	float trigPrescale_HLT_g35_loose_L1EM15; SetInputBranch(inputTree, "trigPrescale_HLT_g35_loose_L1EM15", &trigPrescale_HLT_g35_loose_L1EM15);
	float trigPrescale_HLT_g40_loose_L1EM15; SetInputBranch(inputTree, "trigPrescale_HLT_g40_loose_L1EM15", &trigPrescale_HLT_g40_loose_L1EM15);
	float trigPrescale_HLT_g45_loose_L1EM15; SetInputBranch(inputTree, "trigPrescale_HLT_g45_loose_L1EM15", &trigPrescale_HLT_g45_loose_L1EM15);
	float trigPrescale_HLT_g50_loose_L1EM15; SetInputBranch(inputTree, "trigPrescale_HLT_g50_loose_L1EM15", &trigPrescale_HLT_g50_loose_L1EM15);
	float trigPrescale_HLT_g60_loose; SetInputBranch(inputTree, "trigPrescale_HLT_g60_loose", &trigPrescale_HLT_g60_loose);
	float trigPrescale_HLT_g70_loose; SetInputBranch(inputTree, "trigPrescale_HLT_g70_loose", &trigPrescale_HLT_g70_loose);
	float trigPrescale_HLT_g80_loose; SetInputBranch(inputTree, "trigPrescale_HLT_g80_loose", &trigPrescale_HLT_g80_loose);
	float trigPrescale_HLT_g100_loose; SetInputBranch(inputTree, "trigPrescale_HLT_g100_loose", &trigPrescale_HLT_g100_loose);
	float trigPrescale_HLT_g120_loose; SetInputBranch(inputTree, "trigPrescale_HLT_g120_loose", &trigPrescale_HLT_g120_loose);
	float trigPrescale_HLT_g140_loose; SetInputBranch(inputTree, "trigPrescale_HLT_g140_loose", &trigPrescale_HLT_g140_loose);
    Double_t genWeight;
    Double_t eventWeight;
    Double_t jvtWeight;
    Double_t bTagWeight;
    Double_t pileupWeight;
    if (isData == "MC") {
        SetInputBranch(inputTree, "genWeight", &genWeight);
        SetInputBranch(inputTree, "eventWeight", &eventWeight);
        SetInputBranch(inputTree, "jvtWeight", &jvtWeight);
        SetInputBranch(inputTree, "bTagWeight", &bTagWeight);
        SetInputBranch(inputTree, "pileupWeight", &pileupWeight);
    }
	double totalWeight = 0.; BaselineTree->Branch("totalWeight",&totalWeight,"totalWeight/D");

    //--- event selection
    Int_t nLep_signal; SetInputBranch(inputTree, "nLep_signal", &nLep_signal);
    Int_t nLep_base; SetInputBranch(inputTree, "nLep_base", &nLep_base);
    Int_t jet_n; CopyBranch(inputTree, BaselineTree, "nJet30", "jet_n", &jet_n, "I");
    std::vector<int>* lepFlavor = new std::vector<int>(10); CopyBranch(inputTree, BaselineTree, "lepFlavor", "lepFlavor", &lepFlavor, "std::vector<int>");
    std::vector<int>* lepCharge = new std::vector<int>(10); CopyBranch(inputTree, BaselineTree, "lepCharge", "lepCharge", &lepCharge, "std::vector<int>");
    std::vector<int>* lepIsoFCTight = new std::vector<int>(10); CopyBranch(inputTree, BaselineTree, "lepIsoFCTight", "lepIsoFCTight", &lepCharge, "std::vector<int>");
    std::vector<float>* lep_pT = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "lepPt", "lep_pT", &lep_pT, "std::vector<float>");
    std::vector<float>* lep_eta = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "lepEta", "lep_eta", &lep_eta, "std::vector<float>");
    std::vector<float>* lep_phi = new std::vector<float>(10); CopyBranch(inputTree, BaselineTree, "lepPhi", "lep_phi", &lep_phi, "std::vector<float>");

    //--- MET components
    float MET; CopyBranch(inputTree, BaselineTree, "met_Et", "MET_raw", &MET, "F");
    float MET_phi; SetInputBranch(inputTree, "met_Phi", &MET_phi);
	float gamma_pt; CopyBranch(inputTree, BaselineTree, "PhotonPt", "gamma_pt", &gamma_pt, "F");
	float gamma_eta; CopyBranch(inputTree, BaselineTree, "PhotonEta", "gamma_eta", &gamma_eta, "F");
	float gamma_phi; CopyBranch(inputTree, BaselineTree, "PhotonPhi", "gamma_phi", &gamma_phi, "F");
	float METl = 0.; BaselineTree->Branch("METl_raw",&METl,"METl_raw/F");
	float METt = 0.; BaselineTree->Branch("METt_raw",&METt,"METt_raw/F");

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

	TLorentzVector obj_4vec;
	TLorentzVector isr_4vec;
	TLorentzVector z_4vec;

	Long64_t nentries = inputTree->GetEntries();

	for (Long64_t i=0;i<nentries;i+=event_interval) {

		if (fmod(i,1e5)==0) std::cout << i << " events processed." << std::endl;
		inputTree->GetEntry(i);

        //--- event selection
		if (nLep_signal > 0) continue;
		if (jet_n < 1) continue;
		if (gamma_pt<15.) continue;

        //--- evaluate weight
        totalWeight = 0;
		if (trigMatch_HLT_g15_loose_L1EM7 ==1 && gamma_pt>(15) && gamma_pt<(25+5)) totalWeight = trigPrescale_HLT_g15_loose_L1EM7;
		if (trigMatch_HLT_g25_loose_L1EM15==1 && gamma_pt>(25+5) && gamma_pt<(35+5)) totalWeight = trigPrescale_HLT_g25_loose_L1EM15;
		if (trigMatch_HLT_g35_loose_L1EM15==1 && gamma_pt>(35+5) && gamma_pt<(40+5)) totalWeight = trigPrescale_HLT_g35_loose_L1EM15;
		if (trigMatch_HLT_g40_loose_L1EM15==1 && gamma_pt>(40+5) && gamma_pt<(45+5)) totalWeight = trigPrescale_HLT_g40_loose_L1EM15;
		if (trigMatch_HLT_g45_loose_L1EM15==1 && gamma_pt>(45+5) && gamma_pt<(50+5)) totalWeight = trigPrescale_HLT_g45_loose_L1EM15;
		if (trigMatch_HLT_g50_loose_L1EM15==1 && gamma_pt>(50+5) && gamma_pt<(60+5)) totalWeight = trigPrescale_HLT_g50_loose_L1EM15;
		if (trigMatch_HLT_g60_loose==1 && gamma_pt>(60+5) && gamma_pt<(70+5)) totalWeight = trigPrescale_HLT_g60_loose;
		if (trigMatch_HLT_g70_loose==1 && gamma_pt>(70+5) && gamma_pt<(80+5)) totalWeight = trigPrescale_HLT_g70_loose;
		if (trigMatch_HLT_g80_loose==1 && gamma_pt>(80+5) && gamma_pt<(100+5)) totalWeight = trigPrescale_HLT_g80_loose;
		if (trigMatch_HLT_g100_loose==1 && gamma_pt>(100+5) && gamma_pt<(140+5)) totalWeight = trigPrescale_HLT_g100_loose;
		if (trigMatch_HLT_g140_loose==1 && gamma_pt>(140+5)) totalWeight = trigPrescale_HLT_g140_loose;
		if (totalWeight==0) continue;

		if (isData == "MC") {
            totalWeight = lumi * genWeight * eventWeight * jvtWeight * bTagWeight * pileupWeight;
			if( TString(sampleID).Contains("Vg") ) totalWeight = -1.0 * totalWeight;
		}
		totalWeight = totalWeight*event_interval;

        //--- compute MET parallel and perpendicular components
		METt = MET*TMath::Sin(MET_phi-gamma_phi);
		METl = MET*TMath::Cos(MET_phi-gamma_phi);

		BaselineTree->Fill();
	}

	BaselineTree->Write();

	std::cout << "done." << std::endl;
	outputFile.Close();
	delete inputFile;
}

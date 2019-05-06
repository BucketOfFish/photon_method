#include "../Settings.C"
#include "../CommonFunctions/CommonLibraries.C"
#include "../CommonFunctions/CommonFunctions.C"

using namespace std;

vector<string> noSampleWeight;
vector<string> noEventWeight;

void GetPhotonEvents(string sampleID, string outputName, string pathToNtuples, string isData, string treeName = "tree_NoSys" ) {

    float mylumi = 1.0;
	if     ( TString(outputName).Contains("mc16cd_2018") ) mylumi =  6461;
	else if( TString(outputName).Contains("mc16cd")      ) mylumi = 44000;
   	else if( TString(outputName).Contains("mc16a")       ) mylumi = 36100;

    //---------------------------------------------
    // open input and output files, get TTrees
    //---------------------------------------------

	TH1::SetDefaultSumw2();

	string filename = Form("%s%s.root",pathToNtuples.c_str(), sampleID.c_str()); 
	TFile* inputFile = TFile::Open(filename.c_str());
	TTree* inputTree = (TTree*)inputFile->Get(treeName.c_str());

	Float_t _nGenEvents = 1.;

	cout << endl;
	cout << "Opening file           : " << filename        << endl;
	cout << "Events in ntuple       : " << inputTree->GetEntries() << endl;
	cout << "Total generated events : " << _nGenEvents     << endl;
	cout << "Output type            : " << outputName      << endl;
	cout << "Output path            : " << ntuple_path      << endl;
	cout << "updated event weights" << endl;
	cout << "using luminosity       : " << mylumi          << endl;
    string outfilename = ntuple_path + "/" + outputName + "/" + sampleID.c_str() + ".root";
    cout << "Writing to : " << outfilename << endl;
    TFile outputFile( outfilename.c_str(), "recreate" );
    TTree* BaselineTree = new TTree("BaselineTree", "baseline tree");
	
    //-----------------------------
    // access and copy over existing branches
    //-----------------------------
		
	inputTree->SetBranchStatus("*", 0);

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
    BaselineTree->Branch("genWeight",                &genWeight            ,"genWeight/D"                 );
    BaselineTree->Branch("eventWeight",              &eventWeight          ,"eventWeight/D"               );
    BaselineTree->Branch("jvtWeight",                &jvtWeight            ,"jvtWeight/D"                 );
    BaselineTree->Branch("bTagWeight",               &bTagWeight           ,"bTagWeight/D"                );
    BaselineTree->Branch("pileupWeight",             &pileupWeight         ,"pileupWeight/D"              );

    Int_t RunNumber; CopyBranch(inputTree, BaselineTree, "RunNumber", "RunNumber", &RunNumber, "I");
    Int_t nLep_signal; CopyBranch(inputTree, BaselineTree, "nLep_signal", "nLep_signal", &nLep_signal, "I");
    Int_t nLep_base; CopyBranch(inputTree, BaselineTree, "nLep_base", "nLep_base", &nLep_base, "I");
    float MET; CopyBranch(inputTree, BaselineTree, "met_Et", "MET_raw", &MET, "F");
    float MET_phi; CopyBranch(inputTree, BaselineTree, "met_Phi", "MET_Phi", &MET_phi, "F");
    float MET_softTerm; CopyBranch(inputTree, BaselineTree, "TST_Et", "MET_softTerm", &MET_softTerm, "F");
    float MET_softPhi; CopyBranch(inputTree, BaselineTree, "TST_Phi", "MET_softPhi", &MET_softPhi, "F");
    float HT; CopyBranch(inputTree, BaselineTree, "Ht30", "HT", &HT, "F");
    Float_t mll; CopyBranch(inputTree, BaselineTree, "mll", "mll", &mll, "F");
    float mjj; CopyBranch(inputTree, BaselineTree, "mjj", "mjj", &mjj, "F");
    Int_t jet_n; CopyBranch(inputTree, BaselineTree, "nJet30", "jet_n", &jet_n, "I");
    Int_t bjet_n; CopyBranch(inputTree, BaselineTree, "nBJet30_MV2c10_FixedCutBEff_77", "bjet_n", &bjet_n, "I");
    int nBJet20_MV2c10_FixedCutBEff_77; CopyBranch(inputTree, BaselineTree, "nBJet20_MV2c10_FixedCutBEff_77", "nBJet20_MV2c10_FixedCutBEff_77", &nBJet20_MV2c10_FixedCutBEff_77, "I");

	float PhotonPt; SetInputBranch(inputTree, "PhotonPt", &PhotonPt);
	float PhotonEta; SetInputBranch(inputTree, "PhotonEta", &PhotonEta);
	float PhotonPhi; SetInputBranch(inputTree, "PhotonPhi", &PhotonPhi);

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

	int trigMatch_HLT_g15_loose_L1EM7; CopyBranch(inputTree, BaselineTree, "trigMatch_HLT_g15_loose_L1EM7", "trigMatch_HLT_g15_loose_L1EM7", &trigMatch_HLT_g15_loose_L1EM7, "I");
	int trigMatch_HLT_g25_loose_L1EM15; CopyBranch(inputTree, BaselineTree, "trigMatch_HLT_g25_loose_L1EM15", "trigMatch_HLT_g25_loose_L1EM15", &trigMatch_HLT_g25_loose_L1EM15, "I");
	int trigMatch_HLT_g35_loose_L1EM15; CopyBranch(inputTree, BaselineTree, "trigMatch_HLT_g35_loose_L1EM15", "trigMatch_HLT_g35_loose_L1EM15", &trigMatch_HLT_g35_loose_L1EM15, "I");
	int trigMatch_HLT_g40_loose_L1EM15; CopyBranch(inputTree, BaselineTree, "trigMatch_HLT_g40_loose_L1EM15", "trigMatch_HLT_g40_loose_L1EM15", &trigMatch_HLT_g40_loose_L1EM15, "I");
	int trigMatch_HLT_g45_loose_L1EM15; CopyBranch(inputTree, BaselineTree, "trigMatch_HLT_g45_loose_L1EM15", "trigMatch_HLT_g45_loose_L1EM15", &trigMatch_HLT_g45_loose_L1EM15, "I");
	int trigMatch_HLT_g50_loose_L1EM15; CopyBranch(inputTree, BaselineTree, "trigMatch_HLT_g50_loose_L1EM15", "trigMatch_HLT_g50_loose_L1EM15", &trigMatch_HLT_g50_loose_L1EM15, "I");
	int trigMatch_HLT_g60_loose; CopyBranch(inputTree, BaselineTree, "trigMatch_HLT_g60_loose", "trigMatch_HLT_g60_loose", &trigMatch_HLT_g60_loose, "I");
	int trigMatch_HLT_g70_loose; CopyBranch(inputTree, BaselineTree, "trigMatch_HLT_g70_loose", "trigMatch_HLT_g70_loose", &trigMatch_HLT_g70_loose, "I");
	int trigMatch_HLT_g80_loose; CopyBranch(inputTree, BaselineTree, "trigMatch_HLT_g80_loose", "trigMatch_HLT_g80_loose", &trigMatch_HLT_g80_loose, "I");
	int trigMatch_HLT_g100_loose; CopyBranch(inputTree, BaselineTree, "trigMatch_HLT_g100_loose", "trigMatch_HLT_g100_loose", &trigMatch_HLT_g100_loose, "I");
	int trigMatch_HLT_g120_loose; CopyBranch(inputTree, BaselineTree, "trigMatch_HLT_g120_loose", "trigMatch_HLT_g120_loose", &trigMatch_HLT_g120_loose, "I");
	int trigMatch_HLT_g140_loose; CopyBranch(inputTree, BaselineTree, "trigMatch_HLT_g140_loose", "trigMatch_HLT_g140_loose", &trigMatch_HLT_g140_loose, "I");

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

	//-----------------------------
	// add new branches
	//-----------------------------
	
	float METl = 0.; BaselineTree->Branch("METl_raw",&METl,"METl_raw/F");
	float METt = 0.; BaselineTree->Branch("METt_raw",&METt,"METt_raw/F");
	float MinDPhi_PhotonJet = 0.; BaselineTree->Branch("MinDPhi_PhotonJet",&MinDPhi_PhotonJet,"MinDPhi_PhotonJet/F");
	float MT; BaselineTree->Branch("MT",&MT,"MT/F");
	float gamma_pt = 0.; BaselineTree->Branch("gamma_pt",&gamma_pt,"gamma_pt/F");
	float gamma_eta = 0.; BaselineTree->Branch("gamma_eta",&gamma_eta,"gamma_eta/F");
	float gamma_phi = 0.; BaselineTree->Branch("gamma_phi",&gamma_phi,"gamma_phi/F");

	double totalWeight = 0.; BaselineTree->Branch("totalWeight",&totalWeight,"totalWeight/D");

	//-----------------------------
	// loop over events
	//-----------------------------

	TLorentzVector obj_4vec;
	TLorentzVector isr_4vec;
	TLorentzVector z_4vec;
    std::vector<float>* photon_pT = new std::vector<float>(10);
    std::vector<float>* photon_eta = new std::vector<float>(10);
    std::vector<float>* photon_phi = new std::vector<float>(10);

	Long64_t nentries = inputTree->GetEntries();

	for (Long64_t i=0;i<nentries;i+=event_interval) {

		if (fmod(i,1e5)==0) std::cout << i << " events processed." << std::endl;
		inputTree->GetEntry(i);

		photon_pT->clear();
		photon_eta->clear();
		photon_phi->clear();
		photon_pT->push_back(PhotonPt);
		photon_eta->push_back(PhotonEta);
		photon_phi->push_back(PhotonPhi);

		if (nLep_signal > 0) continue;
		if (jet_n < 1) continue;
		if (photon_pT->size()==0) continue;
		if (photon_pT->at(0)<15.) continue;
		
		double trigWeight = 0;

		if (trigMatch_HLT_g15_loose_L1EM7 ==1 && photon_pT->at(0)>(15) && photon_pT->at(0)<(25+5)) trigWeight = trigPrescale_HLT_g15_loose_L1EM7;
		if (trigMatch_HLT_g25_loose_L1EM15==1 && photon_pT->at(0)>(25+5) && photon_pT->at(0)<(35+5)) trigWeight = trigPrescale_HLT_g25_loose_L1EM15;
		if (trigMatch_HLT_g35_loose_L1EM15==1 && photon_pT->at(0)>(35+5) && photon_pT->at(0)<(40+5)) trigWeight = trigPrescale_HLT_g35_loose_L1EM15;
		if (trigMatch_HLT_g40_loose_L1EM15==1 && photon_pT->at(0)>(40+5) && photon_pT->at(0)<(45+5)) trigWeight = trigPrescale_HLT_g40_loose_L1EM15;
		if (trigMatch_HLT_g45_loose_L1EM15==1 && photon_pT->at(0)>(45+5) && photon_pT->at(0)<(50+5)) trigWeight = trigPrescale_HLT_g45_loose_L1EM15;
		if (trigMatch_HLT_g50_loose_L1EM15==1 && photon_pT->at(0)>(50+5) && photon_pT->at(0)<(60+5)) trigWeight = trigPrescale_HLT_g50_loose_L1EM15;
		if (trigMatch_HLT_g60_loose==1 && photon_pT->at(0)>(60+5) && photon_pT->at(0)<(70+5)) trigWeight = trigPrescale_HLT_g60_loose;
		if (trigMatch_HLT_g70_loose==1 && photon_pT->at(0)>(70+5) && photon_pT->at(0)<(80+5)) trigWeight = trigPrescale_HLT_g70_loose;
		if (trigMatch_HLT_g80_loose==1 && photon_pT->at(0)>(80+5) && photon_pT->at(0)<(100+5)) trigWeight = trigPrescale_HLT_g80_loose;
		if (trigMatch_HLT_g100_loose==1 && photon_pT->at(0)>(100+5) && photon_pT->at(0)<(140+5)) trigWeight = trigPrescale_HLT_g100_loose;
		if (trigMatch_HLT_g140_loose==1 && photon_pT->at(0)>(140+5)) trigWeight = trigPrescale_HLT_g140_loose;
		if (trigWeight==0) continue;

		gamma_pt = photon_pT->at(0);
		gamma_eta = photon_eta->at(0);
		gamma_phi = photon_phi->at(0);

		totalWeight = trigWeight;

		// here we compute the MET parallel and perpendicular components
		METt = MET*TMath::Sin(MET_phi-gamma_phi);
		METl = MET*TMath::Cos(MET_phi-gamma_phi);

		TLorentzVector gamma_4vec;
		gamma_4vec.SetPtEtaPhiM(gamma_pt,gamma_eta,gamma_phi,0);

		if (isData == "MC") {
            totalWeight = mylumi * genWeight * eventWeight * jvtWeight * bTagWeight * pileupWeight;
			if( TString(sampleID).Contains("Vg") ) totalWeight = -1.0 * totalWeight;
		}
		totalWeight = totalWeight*event_interval;

		MinDPhi_PhotonJet = 1000.;
		TLorentzVector jet_4vec;
		for (unsigned int j=0;j<jet_pT->size();j++) {
			jet_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),jet_m->at(j));
			float DR_PhotonJet = jet_4vec.DeltaR(gamma_4vec);
			float DPhi_PhotonJet = jet_4vec.DeltaPhi(gamma_4vec);
			if (MinDPhi_PhotonJet>DPhi_PhotonJet) MinDPhi_PhotonJet = DPhi_PhotonJet;
		}

		TLorentzVector met_4vec;
		met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
		TLorentzVector lep_4vec;
		if (lep_pT->size()>0) lep_4vec.SetPtEtaPhiM(lep_pT->at(0),0,lep_phi->at(0),0);  // only transverse component
		else lep_4vec.SetPtEtaPhiM(0,0,0,0);
		if (lep_pT->size()>0) MT = (met_4vec+lep_4vec).M();
		else MT = 0;

		BaselineTree->Fill();

	}

	BaselineTree->Write();

	std::cout << "done." << std::endl;
	outputFile.Close();
	delete inputFile;

}

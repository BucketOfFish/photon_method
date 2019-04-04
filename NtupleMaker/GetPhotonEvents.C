//-----------------------------------------------------------------------------------------------
// this script takes the photon ntuples from QuickAna and generates smaller ntuples with information required by photon method
// Parameters for: GetPhotonEvents(string sampleID, string outputName, string pathToNtuples, int isData, string treeName = "tree_NoSys" ): 
// 	sampleID: DSID of the MC sample, input ntuple filename in most instances
//      outputName: pretty self-explanatory 
// 	pathToNtuples: the path to the input ntuples  
// 	isData:  "1" for data samples, "0" for MC 
//      treeName: Singlephoton2##_NoSys, you can find this out in the ntuple 
// Example usage: root -l -b -q 'GetPhotonEvents.C+("SinglePhoton211_merged_processed","gmc","/afs/cern.ch/user/b/benhoob/SusySkim2LJets/v1.2/JETM4/JETM4_mc16a/JETM4_mc16a_v1.2_v2/merged/",0,"SinglePhoton211_NoSys")'
//-----------------------------------------------------------------------------------------------

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
	// open file, get Tree and EventCountHist
	//---------------------------------------------

	TH1::SetDefaultSumw2();

	string  filename       = Form("%s%s.root",pathToNtuples.c_str(),sampleID.c_str()); 
	TFile*  inputFile      = TFile::Open(filename.c_str());
	//TH1D*   EventCountHist = (TH1D*) inputFile->Get("EventCountHist");
	//Float_t _nGenEvents    = EventCountHist->GetBinContent(2);
	Float_t _nGenEvents    = 1.;
	TTree*  inputTree              = (TTree*)inputFile->Get( treeName.c_str() );

	cout << endl;
	cout << "Opening file           : " << filename        << endl;
	cout << "Events in ntuple       : " << inputTree->GetEntries() << endl;
	cout << "Total generated events : " << _nGenEvents     << endl;
	cout << "Output type            : " << outputName      << endl;
	cout << "Output path            : " << ntuple_path      << endl;
	cout << "updated event weights" << endl;
	cout << "using luminosity       : " << mylumi          << endl;
	
	//-----------------------------
	// access existing branches
	//-----------------------------
		
	inputTree->SetBranchStatus("*", 0);
	ULong64_t EventNumber; SetInputBranch(inputTree, "EventNumber", &EventNumber);
	Int_t RunNumber; SetInputBranch(inputTree, "RunNumber", &RunNumber);
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
	Float_t Mu; SetInputBranch(inputTree, "mu", &Mu);
	Int_t nVtx; SetInputBranch(inputTree, "nVtx", &nVtx);
	float MET; SetInputBranch(inputTree, "met_Et", &MET); 
	//Variables for 2019 RJR analysis
	Bool_t trigMatch_2LTrigOR; SetInputBranch(inputTree, "trigMatch_2LTrigOR", &trigMatch_2LTrigOR);
	float MET_loose; SetInputBranch(inputTree, "met_Et_loose", &MET_loose);
	float MET_tight; SetInputBranch(inputTree, "met_Et_tight", &MET_tight);
	float MET_tighter; SetInputBranch(inputTree, "met_Et_tighter", &MET_tighter);
	float MET_tenacious; SetInputBranch(inputTree, "met_Et_tenacious", &MET_tenacious);
	Bool_t is2Lep2Jet; SetInputBranch(inputTree, "is2Lep2Jet", &is2Lep2Jet);
	Bool_t is2L2JInt; SetInputBranch(inputTree, "is2L2JInt", &is2L2JInt);
	int nBJet20_MV2c10_FixedCutBEff_77; SetInputBranch(inputTree, "nBJet20_MV2c10_FixedCutBEff_77", &nBJet20_MV2c10_FixedCutBEff_77);
	float mjj; SetInputBranch(inputTree, "mjj", &mjj);
	Double_t mll_RJ; SetInputBranch(inputTree, "mll_RJ", &mll_RJ);
	Double_t R_minH2P_minH3P; SetInputBranch(inputTree, "R_minH2P_minH3P", &R_minH2P_minH3P);
	Double_t RPT_HT5PP; SetInputBranch(inputTree, "RPT_HT5PP", &RPT_HT5PP);
	Double_t dphiVP; SetInputBranch(inputTree, "dphiVP", &dphiVP);
	Double_t H2PP; SetInputBranch(inputTree, "H2PP", &H2PP);
	Double_t H5PP; SetInputBranch(inputTree, "H5PP", &H5PP);
	int nJet20; SetInputBranch(inputTree, "nJet20", &nJet20);
	Double_t minDphi; SetInputBranch(inputTree, "minDphi", &minDphi);
	Double_t MZ; SetInputBranch(inputTree, "MZ", &MZ);
	Double_t NjS; SetInputBranch(inputTree, "NjS", &NjS);
	Double_t NjISR; SetInputBranch(inputTree, "NjISR", &NjISR);
	Double_t dphiISRI; SetInputBranch(inputTree, "dphiISRI", &dphiISRI);
	Double_t RISR; SetInputBranch(inputTree, "RISR", &RISR);
	Double_t PTISR; SetInputBranch(inputTree, "PTISR", &PTISR);
	Double_t PTI; SetInputBranch(inputTree, "PTI", &PTI);
	Double_t PTCM; SetInputBranch(inputTree, "PTCM", &PTCM);
	Double_t MJ; SetInputBranch(inputTree, "MJ", &MJ);
	Bool_t is3Lep3Jet; SetInputBranch(inputTree, "is3Lep3Jet", &is3Lep3Jet);
	Bool_t is4Lep3Jet; SetInputBranch(inputTree, "is4Lep3Jet", &is4Lep3Jet);
	Double_t lept1sign_VR; SetInputBranch(inputTree, "lept1sign_VR", &lept1sign_VR);
	Double_t lept2sign_VR; SetInputBranch(inputTree, "lept2sign_VR", &lept2sign_VR);
	Double_t lept1Pt_VR; SetInputBranch(inputTree, "lept1Pt_VR", &lept1Pt_VR);
	Double_t lept2Pt_VR; SetInputBranch(inputTree, "lept2Pt_VR", &lept2Pt_VR);
	Double_t MZ_VR; SetInputBranch(inputTree, "MZ_VR", &MZ_VR);
	Double_t MJ_VR; SetInputBranch(inputTree, "MJ_VR", &MJ_VR);
	Double_t RISR_VR; SetInputBranch(inputTree, "RISR_VR", &RISR_VR);
	Double_t PTISR_VR; SetInputBranch(inputTree, "PTISR_VR", &PTISR_VR);
	Double_t PTI_VR; SetInputBranch(inputTree, "PTI_VR", &PTI_VR);
	Double_t PTCM_VR; SetInputBranch(inputTree, "PTCM_VR", &PTCM_VR);
	Double_t dphiISRI_VR; SetInputBranch(inputTree, "dphiISRI_VR", &dphiISRI_VR);
	float MET_phi; SetInputBranch(inputTree, "met_Phi", &MET_phi);
	float MET_softTerm; SetInputBranch(inputTree, "TST_Et", &MET_softTerm);
	float MET_softPhi; SetInputBranch(inputTree, "TST_Phi", &MET_softPhi);
	Float_t truthMET = 0; SetInputBranch(inputTree, "truthMET", &truthMET);
	Float_t truthMET_Phi = 0; SetInputBranch(inputTree, "truthMET_Phi", &truthMET_Phi);
	float DPhi_METJetLeading; SetInputBranch(inputTree, "DPhiJ1Met", &DPhi_METJetLeading);
	float DPhi_METJetSecond; SetInputBranch(inputTree, "DPhiJ2Met", &DPhi_METJetSecond);
	float HT; SetInputBranch(inputTree, "Ht30", &HT);
	Int_t jet_n; SetInputBranch(inputTree, "nJet30", &jet_n);
	Int_t bjet_n; SetInputBranch(inputTree, "nBJet30_MV2c10_FixedCutBEff_77", &bjet_n);
	Int_t nLep_signal; SetInputBranch(inputTree, "nLep_signal", &nLep_signal);
	Int_t nLep_base; SetInputBranch(inputTree, "nLep_base", &nLep_base);
	float PhotonPt; SetInputBranch(inputTree, "PhotonPt", &PhotonPt);
	float PhotonEta; SetInputBranch(inputTree, "PhotonEta", &PhotonEta);
	float PhotonPhi; SetInputBranch(inputTree, "PhotonPhi", &PhotonPhi);
	std::vector<int>* lepFlavor = new std::vector<int>(10); SetInputBranch(inputTree, "lepFlavor", &lepFlavor);
	std::vector<int>* lepCharge = new std::vector<int>(10); SetInputBranch(inputTree, "lepCharge", &lepCharge);
	std::vector<int>* lepSignal = new std::vector<int>(10); SetInputBranch(inputTree, "lepSignal", &lepSignal);
	std::vector<float>* lep_pT = new std::vector<float>(10); SetInputBranch(inputTree, "lepPt", &lep_pT);
	std::vector<float>* lep_eta = new std::vector<float>(10); SetInputBranch(inputTree, "lepEta", &lep_eta);
	std::vector<float>* lep_phi = new std::vector<float>(10); SetInputBranch(inputTree, "lepPhi", &lep_phi);
	std::vector<float>* jet_m = new std::vector<float>(10); SetInputBranch(inputTree, "jetM", &jet_m);
	std::vector<float>* jet_pT = new std::vector<float>(10); SetInputBranch(inputTree, "jetPt", &jet_pT);
	std::vector<float>* jet_eta = new std::vector<float>(10); SetInputBranch(inputTree, "jetEta", &jet_eta);
	std::vector<float>* jet_phi = new std::vector<float>(10); SetInputBranch(inputTree, "jetPhi", &jet_phi);
	std::vector<float>* photon_pT = new std::vector<float>(10); SetInputBranch(inputTree, "photon_pT", &photon_pT);
	std::vector<float>* photon_eta = new std::vector<float>(10); SetInputBranch(inputTree, "photon_eta", &photon_eta);
	std::vector<float>* photon_phi = new std::vector<float>(10); SetInputBranch(inputTree, "photon_phi", &photon_phi);

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

	//-----------------------------
	// add new branches
	//-----------------------------

	string outputfilename = ntuple_path + "/" + outputName + "/" + sampleID.c_str() + ".root";
	cout << "Writing to output file : " << outputfilename << endl;
	TFile   outputFile( outputfilename.c_str() , "recreate" );
	
	TTree BaselineTree("BaselineTree","baseline tree");
	float gamma_pt = 0.;
	float gamma_eta = 0.;
	float gamma_phi = 0.;
	float gamma_dR = 999.;
	float METl = 0.;
	float METt = 0.;
	float DPhi_METPhoton = 0.;
	float MinDR_PhotonJet = 0.;
	float MinDPhi_PhotonJet = 0.;
	int pt = 0;
	int ht = 0;
	int njet = 0;
	int nbjet = 0;
	float MT;
	BaselineTree.Branch("pt",&pt,"pt/I");
	BaselineTree.Branch("ht",&ht,"ht/I");
	BaselineTree.Branch("njet",&njet,"njet/I");
	BaselineTree.Branch("nbjet",&nbjet,"nbjet/I");
	BaselineTree.Branch("Mu",&Mu,"Mu/F");
	BaselineTree.Branch("nVtx",&nVtx,"nVtx/I");
	BaselineTree.Branch("MET_raw",&MET,"MET_raw/F");
    //------------------------------------------------------
    // 2019 RJR analysis variables -------------------------
    //------------------------------------------------------
	BaselineTree.Branch("MET_loose",&MET_loose,"MET_loose/F");
	BaselineTree.Branch("MET_tight",&MET_tight,"MET_tight/F");
	BaselineTree.Branch("MET_tighter",&MET_tighter,"MET_tighter/F");
	BaselineTree.Branch("MET_tenacious",&MET_tenacious,"MET_tenacious/F"); 
    BaselineTree.Branch("trigMatch_2LTrigOR",&trigMatch_2LTrigOR,"trigMatch_2LTrigOR/I");
    BaselineTree.Branch("is2Lep2Jet",&is2Lep2Jet,"is2Lep2Jet/I");
    BaselineTree.Branch("is2L2JInt",&is2L2JInt,"is2L2JInt/I");
    BaselineTree.Branch("nBJet20_MV2c10_FixedCutBEff_77",&nBJet20_MV2c10_FixedCutBEff_77,"nBJet20_MV2c10_FixedCutBEff_77/I");
    BaselineTree.Branch("mjj",&mjj,"mjj/F");
    BaselineTree.Branch("mll_RJ",&mll_RJ,"mll_RJ/F");
    BaselineTree.Branch("R_minH2P_minH3P",&R_minH2P_minH3P,"R_minH2P_minH3P/F");
    BaselineTree.Branch("RPT_HT5PP",&RPT_HT5PP,"RPT_HT5PP/F");
    BaselineTree.Branch("dphiVP",&dphiVP,"dphiVP/F");
    BaselineTree.Branch("H2PP",&H2PP,"H2PP/F");
    BaselineTree.Branch("H5PP",&H5PP,"H5PP/F");
    BaselineTree.Branch("nJet20",&nJet20,"nJet20/I");
    BaselineTree.Branch("minDphi",&minDphi,"minDphi/F");
    BaselineTree.Branch("MZ",&MZ,"MZ/F");
    BaselineTree.Branch("NjS",&NjS,"NjS/I");
    BaselineTree.Branch("NjISR",&NjISR,"NjISR/I");
    BaselineTree.Branch("dphiISRI",&dphiISRI,"dphiISRI/F");
    BaselineTree.Branch("RISR",&RISR,"RISR/F");
    BaselineTree.Branch("PTISR",&PTISR,"PTISR/F");
    BaselineTree.Branch("PTI",&PTI,"PTI/F");
    BaselineTree.Branch("PTCM",&PTCM,"PTCM/F");
    BaselineTree.Branch("MJ",&MJ,"MJ/F");
    BaselineTree.Branch("is3Lep3Jet",&is3Lep3Jet,"is3Lep3Jet/I");
    BaselineTree.Branch("is4Lep3Jet",&is4Lep3Jet,"is4Lep3Jet/I");
    BaselineTree.Branch("lept1sign_VR",&lept1sign_VR,"lept1sign_VR/I");
    BaselineTree.Branch("lept2sign_VR",&lept2sign_VR,"lept2sign_VR/I");
    BaselineTree.Branch("lept1Pt_VR",&lept1Pt_VR,"lept1Pt_VR/F");
    BaselineTree.Branch("lept2Pt_VR",&lept2Pt_VR,"lept2Pt_VR/F");
    BaselineTree.Branch("MZ_VR",&MZ_VR,"MZ_VR/F");
    BaselineTree.Branch("MJ_VR",&MJ_VR,"MJ_VR/F");
    BaselineTree.Branch("RISR_VR",&RISR_VR,"RISR_VR/F");
    BaselineTree.Branch("PTISR_VR",&PTISR_VR,"PTISR_VR/F");
    BaselineTree.Branch("PTI_VR",&PTI_VR,"PTI_VR/F");
    BaselineTree.Branch("PTCM_VR",&PTCM_VR,"PTCM_VR/F");
    BaselineTree.Branch("dphiISRI_VR",&dphiISRI_VR,"dphiISRI_VR/F");
    BaselineTree.Branch("lepFlavor","std::vector<int>",&lepFlavor);
    BaselineTree.Branch("lepCharge","std::vector<int>",&lepCharge);
    //------------------------------------------------------
	BaselineTree.Branch("METl_raw",&METl,"METl_raw/F");
	BaselineTree.Branch("METt_raw",&METt,"METt_raw/F");
	BaselineTree.Branch("MET_phi_raw",&MET_phi,"MET_phi_raw/F");
	BaselineTree.Branch("DPhi_METJetLeading_raw",&DPhi_METJetLeading,"DPhi_METJetLeading_raw/F");
	BaselineTree.Branch("DPhi_METJetSecond_raw",&DPhi_METJetSecond,"DPhi_METJetSecond_raw/F");
	BaselineTree.Branch("MinDPhi_PhotonJet",&MinDPhi_PhotonJet,"MinDPhi_PhotonJet/F");
	BaselineTree.Branch("HT",&HT,"HT/F");
	BaselineTree.Branch("MT",&MT,"MT/F");
	BaselineTree.Branch("gamma_pt",&gamma_pt,"gamma_pt/F");
	BaselineTree.Branch("gamma_eta",&gamma_eta,"gamma_eta/F");
	BaselineTree.Branch("gamma_phi",&gamma_phi,"gamma_phi/F");
	BaselineTree.Branch("bjet_n",&bjet_n,"bjet_n/I");
	BaselineTree.Branch("jet_n",&jet_n,"jet_n/I");
	BaselineTree.Branch("jet_m","std::vector<float>",&jet_m);
	BaselineTree.Branch("jet_pT","std::vector<float>",&jet_pT);
	BaselineTree.Branch("jet_phi","std::vector<float>",&jet_phi);
	BaselineTree.Branch("jet_eta","std::vector<float>",&jet_eta);
	BaselineTree.Branch("EventNumber",&EventNumber,"EventNumber/I");
	BaselineTree.Branch("RunNumber",&RunNumber,"RunNumber/I");
	BaselineTree.Branch("nLep_signal",&nLep_signal,"nLep_signal/I");
	BaselineTree.Branch("nLep_base"  ,&nLep_base  ,"nLep_base/I");
	BaselineTree.Branch("lep_pT_raw","std::vector<float>",&lep_pT);
	BaselineTree.Branch("lep_phi_raw","std::vector<float>",&lep_phi);
	BaselineTree.Branch("lep_eta_raw","std::vector<float>",&lep_eta);

	BaselineTree.Branch("trigMatch_HLT_g15_loose_L1EM7", 		&trigMatch_HLT_g15_loose_L1EM7      ,"trigMatch_HLT_g15_loose_L1EM7/I" 	             );
	BaselineTree.Branch("trigMatch_HLT_g25_loose_L1EM15", 		&trigMatch_HLT_g25_loose_L1EM15     ,"trigMatch_HLT_g25_loose_L1EM15/I" 	     );
	BaselineTree.Branch("trigMatch_HLT_g35_loose_L1EM15", 		&trigMatch_HLT_g35_loose_L1EM15     ,"trigMatch_HLT_g35_loose_L1EM15/I" 	     );
	BaselineTree.Branch("trigMatch_HLT_g40_loose_L1EM15", 		&trigMatch_HLT_g40_loose_L1EM15     ,"trigMatch_HLT_g40_loose_L1EM15/I" 	     );
	BaselineTree.Branch("trigMatch_HLT_g45_loose_L1EM15", 		&trigMatch_HLT_g45_loose_L1EM15     ,"trigMatch_HLT_g45_loose_L1EM15/I" 	     );
	BaselineTree.Branch("trigMatch_HLT_g50_loose_L1EM15", 		&trigMatch_HLT_g50_loose_L1EM15     ,"trigMatch_HLT_g50_loose_L1EM15/I" 	     );
	BaselineTree.Branch("trigMatch_HLT_g60_loose", 			&trigMatch_HLT_g60_loose            ,"trigMatch_HLT_g60_loose/I" 		     );
	BaselineTree.Branch("trigMatch_HLT_g70_loose", 			&trigMatch_HLT_g70_loose            ,"trigMatch_HLT_g70_loose/I" 		     );
	BaselineTree.Branch("trigMatch_HLT_g80_loose", 			&trigMatch_HLT_g80_loose            ,"trigMatch_HLT_g80_loose/I" 		     );
	BaselineTree.Branch("trigMatch_HLT_g100_loose", 		&trigMatch_HLT_g100_loose           ,"trigMatch_HLT_g100_loose/I" 	             );
	BaselineTree.Branch("trigMatch_HLT_g120_loose", 		&trigMatch_HLT_g120_loose           ,"trigMatch_HLT_g120_loose/I" 	             );
	BaselineTree.Branch("trigMatch_HLT_g140_loose", 		&trigMatch_HLT_g140_loose           ,"trigMatch_HLT_g140_loose/I" 	             );

    BaselineTree.Branch("genWeight",                &genWeight            ,"genWeight/D"                 );
    BaselineTree.Branch("eventWeight",              &eventWeight          ,"eventWeight/D"               );
    BaselineTree.Branch("jvtWeight",                &jvtWeight            ,"jvtWeight/D"                 );
    BaselineTree.Branch("bTagWeight",               &bTagWeight           ,"bTagWeight/D"                );
    BaselineTree.Branch("pileupWeight",             &pileupWeight         ,"pileupWeight/D"              );

	double totalWeight = 0.;
	BaselineTree.Branch("totalWeight",&totalWeight,"totalWeight/D");

	TH1D* hist_low_njet = new TH1D("hist_low_njet","",bin_size,njet_bin);
	hist_low_njet->SetStats(0);
	TH1D* hist_low_nbjet = new TH1D("hist_low_nbjet","",bin_size,njet_bin);
	hist_low_nbjet->SetStats(0);
	TH1D* hist_low_pt = new TH1D("hist_low_pt","",bin_size,pt_bin);
	hist_low_pt->SetStats(0);
	TH1D* hist_sm_pt = new TH1D("hist_sm_pt","",bin_size,sm_pt_bin);
	hist_sm_pt->SetStats(0);
	TH1D* hist_low_et = new TH1D("hist_low_et","",bin_size,pt_bin);
	hist_low_et->SetStats(0);
	TH1D* hist_low_ht = new TH1D("hist_low_ht","",bin_size,ht_bin);
	hist_low_ht->SetStats(0);
	TH1D* hist_medium_njet = new TH1D("hist_medium_njet","",bin_size,njet_bin);
	hist_medium_njet->SetStats(0);
	TH1D* hist_medium_nbjet = new TH1D("hist_medium_nbjet","",bin_size,njet_bin);
	hist_medium_nbjet->SetStats(0);
	TH1D* hist_medium_pt = new TH1D("hist_medium_pt","",bin_size,pt_bin);
	hist_medium_pt->SetStats(0);
	TH1D* hist_medium_et = new TH1D("hist_medium_et","",bin_size,pt_bin);
	hist_medium_et->SetStats(0);
	TH1D* hist_medium_ht = new TH1D("hist_medium_ht","",bin_size,ht_bin);
	hist_medium_ht->SetStats(0);
	TH1D* hist_high_njet = new TH1D("hist_high_njet","",bin_size,njet_bin);
	hist_high_njet->SetStats(0);
	TH1D* hist_high_nbjet = new TH1D("hist_high_nbjet","",bin_size,njet_bin);
	hist_high_nbjet->SetStats(0);
	TH1D* hist_high_pt = new TH1D("hist_high_pt","",bin_size,pt_bin);
	hist_high_pt->SetStats(0);
	TH1D* hist_high_et = new TH1D("hist_high_et","",bin_size,pt_bin);
	hist_high_et->SetStats(0);
	TH1D* hist_high_ht = new TH1D("hist_high_ht","",bin_size,ht_bin);
	hist_high_ht->SetStats(0);
	TH1D* hist_zmet_njet = new TH1D("hist_zmet_njet","",bin_size,njet_bin);
	hist_zmet_njet->SetStats(0);
	TH1D* hist_zmet_nbjet = new TH1D("hist_zmet_nbjet","",bin_size,njet_bin);
	hist_zmet_nbjet->SetStats(0);
	TH1D* hist_zmet_pt = new TH1D("hist_zmet_pt","",bin_size,pt_bin);
	hist_zmet_pt->SetStats(0);
	TH1D* hist_zmet_et = new TH1D("hist_zmet_et","",bin_size,pt_bin);
	hist_zmet_et->SetStats(0);
	TH1D* hist_zmet_ht = new TH1D("hist_zmet_ht","",bin_size,ht_bin);
	hist_zmet_ht->SetStats(0);
	TH1D* hist_bveto_njet = new TH1D("hist_bveto_njet","",bin_size,njet_bin);
	hist_bveto_njet->SetStats(0);
	TH1D* hist_bveto_nbjet = new TH1D("hist_bveto_nbjet","",bin_size,njet_bin);
	hist_bveto_nbjet->SetStats(0);
	TH1D* hist_bveto_pt = new TH1D("hist_bveto_pt","",bin_size,pt_bin);
	hist_bveto_pt->SetStats(0);
	TH1D* hist_bveto_et = new TH1D("hist_bveto_et","",bin_size,pt_bin);
	hist_bveto_et->SetStats(0);
	TH1D* hist_bveto_ht = new TH1D("hist_bveto_ht","",bin_size,ht_bin);
	hist_bveto_ht->SetStats(0);

	TLorentzVector obj_4vec;
	TLorentzVector isr_4vec;
	TLorentzVector z_4vec;

	//-----------------------------
	// loop over events
	//-----------------------------

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

		if (trigMatch_HLT_g15_loose_L1EM7 ==1 && photon_pT->at(0)>(15)   && photon_pT->at(0)<(25+5)) trigWeight = trigPrescale_HLT_g15_loose_L1EM7;
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

		njet = hist_low_njet->FindBin(jet_n)-1;
		if (jet_n>njet_bin[bin_size]) njet = bin_size-1;
		nbjet = hist_low_nbjet->FindBin(bjet_n)-1;
		if (bjet_n>njet_bin[bin_size]) nbjet = bin_size-1;
		pt = hist_low_pt->FindBin(gamma_pt)-1;
		if (gamma_pt>pt_bin[bin_size]) pt = bin_size-1;
		ht = hist_low_ht->FindBin(HT)-1;
		if (HT>ht_bin[bin_size]) ht = bin_size-1;

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

		DPhi_METPhoton = fabs(TMath::ATan2(METt,METl));
		MinDR_PhotonJet = 1000.;
		MinDPhi_PhotonJet = 1000.;
		TLorentzVector jet_4vec;
		for (unsigned int j=0;j<jet_pT->size();j++) {
			jet_4vec.SetPtEtaPhiM(jet_pT->at(j),jet_eta->at(j),jet_phi->at(j),jet_m->at(j));
			float DR_PhotonJet = jet_4vec.DeltaR(gamma_4vec);
			float DPhi_PhotonJet = jet_4vec.DeltaPhi(gamma_4vec);
			if (MinDR_PhotonJet>DR_PhotonJet) MinDR_PhotonJet = DR_PhotonJet;
			if (MinDPhi_PhotonJet>DPhi_PhotonJet) MinDPhi_PhotonJet = DPhi_PhotonJet;
		}

		TLorentzVector met_4vec;
		met_4vec.SetPtEtaPhiM(MET,0,MET_phi,0);
		TLorentzVector lep_4vec;
		if (lep_pT->size()>0) lep_4vec.SetPtEtaPhiM(lep_pT->at(0),0,lep_phi->at(0),0);  // only transverse component
		else lep_4vec.SetPtEtaPhiM(0,0,0,0);
		if (lep_pT->size()>0) MT = (met_4vec+lep_4vec).M();
		else MT = 0;

		BaselineTree.Fill();

	}

	BaselineTree.Write();

	std::cout << "done." << std::endl;
	outputFile.Close();
	delete inputFile;

}

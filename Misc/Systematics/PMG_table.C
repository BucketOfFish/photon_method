#include "/home/matt/Projects/PhotonMethod/Main/Settings.cpp"

TCut getDPhiCR(string region) {
    return NMinus1Cut(cuts::selections[region], "minDPhi2JetsMet") + cuts::minDPhi2JetsMet_anti0p4;
}

void PMG_table() {
    // Options
    TChain chain("Zjets_NoSys");
    chain.Add("/public/data/SUSY_Systematics/Skimmed/StrongPreselectionInclusive/mc16a_Zjets_merged_processed.root");
    chain.Add("/public/data/SUSY_Systematics/Skimmed/StrongPreselectionInclusive/mc16cd_Zjets_merged_processed.root");
    chain.Add("/public/data/SUSY_Systematics/Skimmed/StrongPreselectionInclusive/mc16e_Zjets_merged_processed.root");

    //vector<string> regions = {"SRC", "SRLow", "SRMed", "SRHigh", "SRLowZ", "SRMedZ", "SRHighZ",
                              ////"CRC", "CRLow", "CRMed", "CRHigh", "CRLowZ", "CRMedZ", "CRHighZ", 
                              //"VRC", "VRLow", "VRMed", "VRHigh", "VRLowZ", "VRMedZ", "VRHighZ"}; 
    vector<string> regions = {"VRHighZ"}; 

    TCut lumi = "(RandomRunNumber<320000 ? 36200 : (RandomRunNumber>320000 && RandomRunNumber<348000) ? 44300 : 58500)";
    TCut MC_weight = "genWeight*eventWeight*leptonWeight*jvtWeight*bTagWeight*pileupWeight*globalDiLepTrigSF";

    //vector<string> systematics = {"scale", "PDF"};
    vector<string> systematics = {"PDF"};
    //systematics = {"ttbar_scale", "ttbar_PDF"};

    // Table Production
    for (auto systematic : systematics) {
        vector<string> LHE3Weights;
        if (systematic == "scale") {
            LHE3Weights.push_back("MUR1_MUF1_PDF261000");
            LHE3Weights.push_back("MUR0_5_MUF0_5_PDF261000");
            LHE3Weights.push_back("MUR0_5_MUF1_PDF261000");
            LHE3Weights.push_back("MUR1_MUF0_5_PDF261000");
            LHE3Weights.push_back("MUR1_MUF2_PDF261000");
            LHE3Weights.push_back("MUR2_MUF1_PDF261000");
            LHE3Weights.push_back("MUR2_MUF2_PDF261000");
        }
        else if (systematic == "PDF") {
            LHE3Weights.push_back("MUR1_MUF1_PDF261000");
            LHE3Weights.push_back("MUR1_MUF1_PDF13000");
            LHE3Weights.push_back("MUR1_MUF1_PDF25300");
            //LHE3Weights.push_back("MUR1_MUF1_PDF269000");
            //LHE3Weights.push_back("MUR1_MUF1_PDF270000");
            //for (int i=261001;i<261101;i++) LHE3Weights.push_back("MUR1_MUF1_PDF" + to_string(i));
        }
        else if (systematic == "ttbar_scale") {
            LHE3Weights.push_back("nominal");
            LHE3Weights.push_back("muR0p5,muF0p5");
            LHE3Weights.push_back("muR0p5,muF1");
            LHE3Weights.push_back("muR1,muF0p5");
            LHE3Weights.push_back("muR0p5,muF2");
            LHE3Weights.push_back("muR2,muF0p5");
            LHE3Weights.push_back("muR1,muF2");
            LHE3Weights.push_back("muR2,muF1");
            LHE3Weights.push_back("muR2,muF2");
        }
        else if (systematic == "ttbar_PDF") {
            LHE3Weights.push_back("nominal");
            LHE3Weights.push_back("PDFset265000");
            LHE3Weights.push_back("PDFset266000");
        }

        for (auto region : regions) {
            TCut SR_cut = cuts::selections[region];
            TCut CR_cut = getDPhiCR(region);
            cout << "Finding Systematics for Region " << region << endl;

            vector<float> ratios;
            for (auto LHE3Weight : LHE3Weights) {
                cout << "Working on Weight " << LHE3Weight << endl;
                TH1F* SR_hist = new TH1F(("hists_" + LHE3Weight + "_SR").c_str(),"",1,0,1);
                TH1F* CR_hist = new TH1F(("hists_" + LHE3Weight + "_CR").c_str(),"",1,0,1);
                TCut LHE_branch = ("LHE3Weight_" + LHE3Weight).c_str();
                chain.Draw(("mll>>hists_" + LHE3Weight + "_SR").c_str(), SR_cut*lumi*MC_weight*LHE_branch, "goff");
                chain.Draw(("mll>>hists_" + LHE3Weight + "_CR").c_str(), CR_cut*lumi*MC_weight*LHE_branch, "goff");
                float SR_yield = fabs(SR_hist->Integral(0,2));
                float CR_yield = fabs(CR_hist->Integral(0,2));
                ratios.push_back(SR_yield / CR_yield);
                cout << "SR yield: " << SR_yield << ", CR yield: " << CR_yield << ", ratio: " << ratios.back() << endl;
                delete SR_hist, CR_hist;
            }

            float uncertainty_max = 0;
            for (int i=0; i<ratios.size(); i++) {
                if (i==0) continue;
                float uncertainty = abs((ratios[i]-ratios[0]) / ratios[0]);
                if (uncertainty>uncertainty_max) uncertainty_max = uncertainty;
            }
            cout << "Region " << region << " " << systematic << " uncertainty: " << uncertainty_max << endl;
        }	
    }
}

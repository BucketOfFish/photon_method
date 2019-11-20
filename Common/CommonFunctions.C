#ifndef COMMON_FUNCTIONS
#define COMMON_FUNCTIONS

template<class variableType>
void SetInputBranch(TTree* inputTree, string branchName, variableType variablePointer) {
    // set variablePointer to reference value in branch if it exists
    for (auto branch : *(inputTree->GetListOfBranches())) {
        if (branch->GetName() == branchName) {
            inputTree->SetBranchStatus(branchName.c_str(), 1);
            inputTree->SetBranchAddress(branchName.c_str(), variablePointer);
        }
    }
}

template<class variableType>
void CopyBranch(TTree* inputTree, TTree* outputTree, string inputBranchName, string outputBranchName, variableType variablePointer, string varType) {
    // copy branch and set variablePointer to reference value in branch if it exists
    for (auto branch : *(inputTree->GetListOfBranches())) {
        if (branch->GetName() == inputBranchName) {
            inputTree->SetBranchStatus(inputBranchName.c_str(), 1);
            inputTree->SetBranchAddress(inputBranchName.c_str(), variablePointer);
            if (varType.compare(0, 6, "vector") == 0)
                outputTree->Branch(outputBranchName.c_str(), varType.c_str(), variablePointer);
            else
                outputTree->Branch(outputBranchName.c_str(), variablePointer, (outputBranchName+"/"+varType).c_str());
        }
    }
}

//  period: data15-16 (input) -> ZMC16a (source file), data17 -> ZMC16cd, data18 -> ZMC16e
float GetLumi(TString period) {
    float lumi = 1.0;
    if (period.Contains("mc16e") || period.Contains("data18")) lumi = 59900;
    else if (period.Contains("mc16cd") || period.Contains("data17")) lumi = 44000;
    else if (period.Contains("mc16a") || period.Contains("data15-16")) lumi = 36100;
    return lumi;
}

TString MCPeriod(TString period) {
    TString mc_period;
    if (period == "data15-16") mc_period = "mc16a";
    else if (period == "data17") mc_period = "mc16cd";
    else if (period == "data18") mc_period = "mc16e";
    else mc_period = period;
    return mc_period;
}

TString DataPeriod(TString period) {
    TString data_period;
    if (period == "mc16a") data_period = "data15-16";
    else if (period == "mc16cd") data_period = "data17";
    else if (period == "mc16e") data_period = "data18";
    else data_period = period;
    return data_period;
}

#endif

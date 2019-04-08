template<class variableType>
void SetInputBranch(TTree* inputTree, string branchName, variableType variablePointer) {
    inputTree->SetBranchStatus(branchName.c_str(), 1);
    inputTree->SetBranchAddress(branchName.c_str(), variablePointer);
}

template<class variableType>
void CopyBranch(TTree* inputTree, TTree* outputTree, string inputBranchName, string outputBranchName, variableType variablePointer, string varType) {
    inputTree->SetBranchStatus(inputBranchName.c_str(), 1);
    inputTree->SetBranchAddress(inputBranchName.c_str(), variablePointer);
    if (varType.compare(0, 11, "std::vector") == 0)
        outputTree->Branch(outputBranchName.c_str(), varType.c_str(), variablePointer);
    else
        outputTree->Branch(outputBranchName.c_str(), variablePointer, (outputBranchName+"/"+varType).c_str());
}

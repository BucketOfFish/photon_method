template<class variableType>
void SetInputBranch(TTree* inputTree, string branchName, variableType variablePointer) {
    inputTree->SetBranchStatus(branchName.c_str(), 1);
    inputTree->SetBranchAddress(branchName.c_str(), variablePointer);
}

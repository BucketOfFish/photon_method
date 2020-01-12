#ifndef COMMON_FUNCTIONS
#define COMMON_FUNCTIONS

using namespace std;

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

vector<int> int_copy_vars;
vector<double> double_copy_vars;
vector<float> float_copy_vars;

void CopyAllBranches(TTree* inputTree, TTree* outputTree, vector<string> branches) {
    int_copy_vars.clear(); int_copy_vars.reserve(100);
    double_copy_vars.clear(); double_copy_vars.reserve(100);
    float_copy_vars.clear(); float_copy_vars.reserve(100);
    for (string branch : branches) {
        string branch_name = branch.substr(0,branch.length()-2);
        char branch_type = branch.substr(branch.length()-1)[0];
        switch (branch_type) {
            case 'I':
                int int_copy_var; int_copy_vars.push_back(int_copy_var);
                CopyBranch(inputTree, outputTree, branch_name, branch_name, &int_copy_vars.back(), "I");
                break;
            case 'D':
                double double_copy_var; double_copy_vars.push_back(double_copy_var);
                CopyBranch(inputTree, outputTree, branch_name, branch_name, &double_copy_vars.back(), "D");
                break;
            case 'F':
                float float_copy_var; float_copy_vars.push_back(float_copy_var);
                CopyBranch(inputTree, outputTree, branch_name, branch_name, &float_copy_vars.back(), "F");
                break;
            default:
                cout << "Unknown branch type" << endl;
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

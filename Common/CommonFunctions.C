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
    int_copy_vars.clear(); int_copy_vars.reserve(branches.size());
    double_copy_vars.clear(); double_copy_vars.reserve(branches.size());
    float_copy_vars.clear(); float_copy_vars.reserve(branches.size());
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

TString getMCPeriod(TString period) {
    TString mc_period;
    if (period == "data15-16") mc_period = "mc16a";
    else if (period == "data17") mc_period = "mc16cd";
    else if (period == "data18") mc_period = "mc16e";
    else mc_period = period;
    return mc_period;
}

string getMCPeriod(string period) {
    return string(getMCPeriod(TString(period)));
}

TString DataPeriod(TString period) {
    TString data_period;
    if (period == "mc16a") data_period = "data15-16";
    else if (period == "mc16cd") data_period = "data17";
    else if (period == "mc16e") data_period = "data18";
    else data_period = period;
    return data_period;
}

template <typename T> constexpr string_view type_name() {
    // from https://stackoverflow.com/questions/81870/is-it-possible-to-print-a-variables-type-in-standard-c/56766138#56766138
    // e.g. cout << type_name<decltype(weighted_dataframe)>() << endl;
    std::string_view name, prefix, suffix;
#ifdef __clang__
    name = __PRETTY_FUNCTION__;
    prefix = "std::string_view type_name() [T = ";
    suffix = "]";
#elif defined(__GNUC__)
    name = __PRETTY_FUNCTION__;
    prefix = "constexpr std::string_view type_name() [with T = ";
    suffix = "; std::string_view = std::basic_string_view<char>]";
#elif defined(_MSC_VER)
    name = __FUNCSIG__;
    prefix = "class std::basic_string_view<char,struct std::char_traits<char> > __cdecl type_name<";
    suffix = ">(void)";
#endif
    name.remove_prefix(prefix.size());
    name.remove_suffix(suffix.size());
    return name;
}

#endif

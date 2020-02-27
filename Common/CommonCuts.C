namespace cuts {

    double leading_lep_pt_cut = 25.; // also used for smearing
    double second_lep_pt_cut = 25.; // also used for smearing

    TCut bkg_baseline("nJet30>=2 && is_OS && lepPt[0]>25.0 && lepPt[1]>25.0 && lepIsoFCTight[0] && lepIsoFCTight[1]");
    TCut photon_baseline("nJet30>=2 && lepPt[0]>25.0 && lepPt[1]>25.0");
    TCut baseline("nJet30>=2 && lepPt[0]>25.0 && lepPt[1]>25.0");

    TCut ee("channel==1");
    TCut mm("channel==0");
    TCut em("channel==2 || channel==3");

    TCut bkg_weight("totalWeight");
    TCut photon_weight("totalWeight");
    TCut photon_weight_rw("totalWeight*reweight_Ptll");

    TCut reweight_region("nJet30>=2 && lepPt[0]>25.0 && lepPt[1]>25.0 && nLep_signal==2");

    TCut CR("met_Et<100.0");

    unordered_map<string, TCut> plot_regions = {
        {"SR", "nJet30>=2 && lepPt[0]>25.0 && lepPt[1]>25.0 && nLep_signal==2"},
        {"VR", "nJet30>=2 && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && mll>80 && mll<100 && mjj<60 && mjj>100"},
        {"VRcom", "nJet30>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && met_Et>100 && met_Et<200"},

        // from https://arxiv.org/pdf/1611.05791.pdf
        {"SRZ2016", "nJet30>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && (mll>81 && mll<101) && dPhiMetJet12Min>0.4 && met_Et>225 && HT>600"},
        {"SRlow2016", "nJet30>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && mll>12 && dPhiMetJet12Min>0.4 && met_Et>200"},
        {"SRmed2016", "nJet30>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && mll>12 && dPhiMetJet12Min>0.4 && met_Et>200 && HT>400"},
        {"SRhigh2016", "nJet30>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && mll>12 && dPhiMetJet12Min>0.4 && met_Et>200 && HT>700"},

        // from https://indico.cern.ch/event/883484/contributions/3722767/attachments/1984038/3305237/20-02-10-2l.pdf
        {"SRC", "nLep_signal==2 && nJet30>=2 && met_Et>250 && mt2leplsp_0>90 && (Ptll>40 && Ptll<100) && MET_sig>10 && mll<50"},
        {"SRCZ", "nLep_signal==2 && nJet30>=4 && met_Et>250 && mt2leplsp_0>90 && (Ptll>40 && Ptll<100) && MET_sig>10 && (mll>81 && mll<101)"},
        {"SRLow4", "nLep_signal==2 && nJet30>=4 && HT>250 && met_Et>250 && mt2leplsp_0>100 && (Ptll>40 && Ptll<500) && (mll<150 && ~(mll>81 && mll<101))"},
        {"SRLowZ", "nLep_signal==2 && nJet30>=6 && HT>250 && met_Et>250 && mt2leplsp_0>100 && (Ptll>40 && Ptll<500) && (mll>81 && mll<101)"},
        {"SRMed4", "nLep_signal==2 && nJet30>=4 && HT>500 && met_Et>300 && mt2leplsp_0>75 && (Ptll>40 && Ptll<800) && (mll>100 && mll<550)"},
        {"SRMedZ", "nLep_signal==2 && nJet30>=6 && HT>500 && met_Et>300 && mt2leplsp_0>75 && (Ptll>40 && Ptll<800) && (mll>81 && mll<101)"},
        {"SRHigh4", "nLep_signal==2 && nJet30>=4 && HT>800 && met_Et>300 && mt2leplsp_0>75 && Ptll>40 && (mll>150 && mll<950)"},
        {"SRHighZ", "nLep_signal==2 && nJet30>=6 && HT>800 && met_Et>300 && mt2leplsp_0>75 && Ptll>40 && (mll>81 && mll<101)"},

        {"VRC", "nLep_signal==2 && nJet30>=2 && (met_Et>150 && met_Et<250) && mt2leplsp_0>90 && (Ptll>40 && Ptll<100) && MET_sig>10 && mll<50"},
        {"VRCZ", "nLep_signal==2 && nJet30>=4 && (met_Et>150 && met_Et<250) && mt2leplsp_0>90 && (Ptll>40 && Ptll<100) && MET_sig>10 && (mll>81 && mll<101)"},
        {"VRLow4", "nLep_signal==2 && nJet30>=4 && HT>250 && (met_Et>150 && met_Et<250) && mt2leplsp_0>100 && (Ptll>40 && Ptll<500) && (mll<150 && ~(mll>81 && mll<101))"},
        {"VRLowZ", "nLep_signal==2 && nJet30>=6 && HT>250 && (met_Et>150 && met_Et<250) && mt2leplsp_0>100 && (Ptll>40 && Ptll<500) && (mll>81 && mll<101)"},
        {"VRMed4", "nLep_signal==2 && nJet30>=4 && HT>500 && (met_Et>200 && met_Et<300) && mt2leplsp_0>75 && (Ptll>40 && Ptll<800) && (mll>100 && mll<550)"},
        {"VRMedZ", "nLep_signal==2 && nJet30>=6 && HT>500 && (met_Et>200 && met_Et<300) && mt2leplsp_0>75 && (Ptll>40 && Ptll<800) && (mll>81 && mll<101)"},
        {"VRHigh4", "nLep_signal==2 && nJet30>=4 && HT>800 && (met_Et>200 && met_Et<300) && mt2leplsp_0>75 && Ptll>40 && (mll>150 && mll<950)"},
        {"VRHighZ", "nLep_signal==2 && nJet30>=6 && HT>800 && (met_Et>200 && met_Et<300) && mt2leplsp_0>75 && Ptll>40 && (mll>81 && mll<101)"},
    };
}

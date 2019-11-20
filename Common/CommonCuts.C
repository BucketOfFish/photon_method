namespace cuts {

    double leading_lep_pt_cut = 25.; // also used for smearing
    double second_lep_pt_cut = 25.; // also used for smearing

    TCut bkg_baseline("nJet30>=2 && is_OS && lepPt[0]>25.0 && lepPt[1]>25.0 && lepIsoFCTight[0] && lepIsoFCTight[1]");
    TCut photon_baseline("nJet30>=2 && lepPt[0]>25.0 && lepPt[1]>25.0");
    TCut baseline("nJet30>=2 && lepPt[0]>25.0 && lepPt[1]>25.0");

    //TCut reweight_region("nJet30>=2 && met_Et>100 && met_Et<200");
    //TCut reweight_region("nJet30>=2 && met_Et<200");
    TCut reweight_region("nJet30>=2 && lepPt[0]>25.0 && lepPt[1]>25.0");

    //TCut CR("met_Et<100.0 && Ptll>25");
    TCut CR("met_Et<100.0");
    //TCut VR("nJet30>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && mll>80 && mll<100 && mjj<60 && mjj>100");  
    TCut VR("nJet30>=2 && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && mll>80 && mll<100 && mjj<60 && mjj>100");  
    TCut SR("nJet30>=2 && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && nBJet20_MV2c10_FixedCutBEff_77==0 && (mjj<60.0 || mjj>100.)");

    TCut VRcom("nJet30>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && met_Et>100 && met_Et<200");  
    // from https://arxiv.org/pdf/1611.05791.pdf
    TCut SRZ2016("nJet30>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && (mll>81 && mll<101) && dPhiMetJet12Min>0.4 && met_Et>225 && HT>600");  
    TCut SRlow2016("nJet30>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && mll>12 && dPhiMetJet12Min>0.4 && met_Et>200");  
    TCut SRmed2016("nJet30>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && mll>12 && dPhiMetJet12Min>0.4 && met_Et>200 && HT>400");  
    TCut SRhigh2016("nJet30>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && mll>12 && dPhiMetJet12Min>0.4 && met_Et>200 && HT>700");  

    TCut ee("channel==1");
    TCut mm("channel==0");
    TCut em("channel==2 || channel==3");

    TCut bkg_weight("totalWeight");
    TCut photon_weight("totalWeight");
    TCut photon_weight_rw("totalWeight*reweight_Ptll");
}

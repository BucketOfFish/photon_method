namespace cuts{
    //TCut Zselection("mll>81 && mll<101 && jet_n >= 2 && MET<200 && is_OS && lep_pT[0]>25.0 && lep_pT[1]>25.0 && bjet_n==0");
    TCut Zselection("jet_n >= 2 && MET<200 && is_OS && lep_pT[0]>25.0 && lep_pT[1]>25.0 && mjj<60 || mjj>100 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && jet_pT[0]>30 && jet_pT[1]>30 && mll>80 && mll<100 && HT>800");
    //TCut Zselection("mll>81 && mll<101 && jet_n >= 2 && is_OS && lep_pT[0]>25.0 && lep_pT[1]>25.0 && bjet_n == 0");

    //TCut gselection("lep_pT[0]>25 && lep_pT[1]>25 && jet_n>=2  && bjet_n==0");
    TCut gselection("jet_n>=2 && mjj<60 || mjj>100 && lep_pT[0]>25 && lep_pT[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && HT>800");
    //TCut gselection("jet_n>=2 && lep_pT[0]>25 && lep_pT[1]>25 && bjet_n == 0");

    //TCut vgselection("jet_n>=2  && bjet_n==0");
    TCut vgselection("jet_n>=2");

    //TCut ZCR("MET<60.0");
    TCut CR("MET<60.0");

    TCut VR("jet_n>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lep_pT[0]>25 && lep_pT[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && mll>80 && mll<100 && mjj<60 && mjj>100");  
    //TCut VR("jet_n>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lep_pT[0]>25 && lep_pT[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && nJet20==2 && mll>80 && mll<100 && mjj<60.0 && mjj>100.");
    TCut SR("jet_n>=2 && lep_pT[0]>25 && lep_pT[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && nBJet20_MV2c10_FixedCutBEff_77==0 && (mjj<60.0 || mjj>100.)");

    TCut ee("channel==1");
    TCut mm("channel==0");
    TCut em("channel==2 || channel==3");

    TCut Zweight("totalWeight");
    TCut weight_g    = "totalWeight";
    TCut weight_g_rw = "totalWeight*ptreweight_step1*ptreweight_step2"; //2-STEP

    TCut Zweight_MC = "totalWeight*36100"; //change luminosity according to data year
    TCut weight_g_MC    = "totalWeight*36100";
    TCut weight_g_rw_MC = "totalWeight*ptreweight_step1*ptreweight_step2*36100";
}

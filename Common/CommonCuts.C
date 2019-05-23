namespace cuts{

    TCut Zselection("nJet30 >= 2 && is_OS && lepPt[0]>25.0 && lepPt[1]>25.0 && bjet_n==0");

    TCut gselection("nJet30>=2 && bjet_n == 0");

    TCut vgselection("nJet30>=2");

    TCut CR("met_Et<60.0");

    TCut VR("nJet30>=2 && lepCharge[0]*lepCharge[1]<0 && abs(lepFlavor[0])==abs(lepFlavor[1]) && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && mll>80 && mll<100 && mjj<60 && mjj>100");  
    TCut SR("nJet30>=2 && lepPt[0]>25 && lepPt[1]>25 && jet_pT[0]>30 && jet_pT[1]>30 && nBJet20_MV2c10_FixedCutBEff_77==0 && (mjj<60.0 || mjj>100.)");

    TCut ee("channel==1");
    TCut mm("channel==0");
    TCut em("channel==2 || channel==3");

    TCut Zweight("totalWeight");
    TCut weight_g = "totalWeight";
    TCut weight_g_rw = "totalWeight*ptreweight_step1*ptreweight_step2";
}

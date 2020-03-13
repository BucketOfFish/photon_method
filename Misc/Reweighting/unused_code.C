    // make canvas and draw 2L data vs. MC plot
    TCanvas *can = new TCanvas("can","can",600,600);
    can->cd();

    gPad->SetLogy();

    hdata->SetLineColor(1);
    hdata->SetLineWidth(2);
    hdata->SetMarkerStyle(20);

    hdata->GetXaxis()->SetTitle("Z p_{T} [GeV]");
    hdata->GetYaxis()->SetTitle("entries / bin");
    hdata->Draw("E1");

    htt->SetLineColor(1);
    htt->SetFillColor(kRed-2);

    hvv->SetLineColor(1);
    hvv->SetFillColor(kGreen-2);

    hz->SetLineColor(1);
    hz->SetFillColor(kOrange-2);

    THStack *mcstack = new THStack("mcstack","mcstack");
    mcstack->Add(htt);
    mcstack->Add(hvv);
    mcstack->Add(hz);
    mcstack->Draw("samehist");
    hdata->Draw("sameE1");
    hdata->Draw("axissame");

    TLegend* leg = new TLegend(0.65,0.65,0.88,0.88);
    leg->AddEntry(hdata,"data","lp");
    leg->AddEntry(hz,"Z+jets","f");
    leg->AddEntry(hvv,"VV","f");
    leg->AddEntry(htt,"t#bar{t}","f");
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->Draw();
    can->Print(Form("%sGetSimpleReweightingHistograms_%s_%s_%s_2L.pdf",plots_path.c_str(),period.c_str(),channel.c_str(),smearing_mode.c_str()));

    //--- take data and remove MC htt and MC hvv
    TH1F* histoZ = (TH1F*) hdata->Clone("histoZ");
    histoZ->Add( htt , -1.0 );
    histoZ->Add( hvv , -1.0 );
    //histoZ->Add( hz  , -1.0 );

    TCanvas *can2 = new TCanvas("can2","can2",600,600);
    can2->cd();

    // float nZ = histoZ->Integral();
    // histoZ->Scale( 1.0 / nZ );

    // float nG = histoG->Integral();
    // histoG->Scale( 1.0 / nG );

    TH1F* hratio = (TH1F*) histoZ->Clone("hratio");
    hratio->Divide( histoG );

    can2->Divide(1,2);
    can2->cd(2);

    histoG->GetXaxis()->SetTitle("Z p_{T} [GeV]");
    histoG->GetYaxis()->SetTitle("entries / bin");

    gPad->SetLogy();

    histoZ->SetLineColor(2);
    histoG->SetLineColor(4);

    histoG->Draw("hist");
    histoZ->Draw("samehist");

    TLegend* leg2 = new TLegend(0.7,0.7,0.88,0.88);
    leg2->AddEntry(histoZ,"2L data - t#bar{t} - VV","f");
    leg2->AddEntry(histoG,"photon","f");
    leg2->SetBorderSize(0);
    leg2->SetFillColor(0);
    leg2->Draw();

    can2->cd(1);
    hratio->SetLineColor(1);
    hratio->Draw("hist");

    can2->Print(Form("%sGetSimpleReweightingHistograms_%s_%s_%s_Z_vs_g.pdf",plots_path.c_str(),period.c_str(),channel.c_str(),smearing_mode.c_str()));

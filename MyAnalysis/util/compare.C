int compare(){
	
    gROOT->SetStyle("Plain");

    // For the canvas:
    gStyle->SetCanvasColor(0);

    // For the Pad:
    gStyle->SetPadColor(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadBorderSize(2);

    // For the frame:
    gStyle->SetFrameBorderMode(0);


    // For the statistics box:
    gStyle->SetOptStat(0);

    // Margins:
    gStyle->SetPadBottomMargin(0.25);
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.1);

    // For the Global title:
    gStyle->SetOptTitle(0);
    gStyle->SetTitleColor(1);
    gStyle->SetTitleFillColor(10);
    gStyle->SetTitleTextColor(1);
    gStyle->SetTitleFont(42);
    gStyle->SetTitleFontSize(0.05);
    gStyle->SetTitleBorderSize(0);

    // For the axis
    gStyle->SetNdivisions(510, "X");
    gStyle->SetNdivisions(510, "Y");
    gStyle->SetTickLength(0.03);

    // For the axis titles:
    gStyle->SetTitleOffset(1.4, "X");
    gStyle->SetTitleOffset(1.2, "Y");
    gStyle->SetTitleOffset(0.5, "Z");
    gStyle->SetTitleSize(0.061, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");

    // For the axis labels:
    gStyle->SetLabelSize(0.04, "XYZ");
    gStyle->SetLabelOffset(0.01, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");

    // For the legend
    gStyle->SetLegendBorderSize(0);

    gROOT->ForceStyle();

    ////////////////////////////////////////

	TFile *f1 = new TFile("output_GetPrediction/prediction_histos_MyTest_data_METsoftSmeared_noAngSmear_N20_CR_v3.root", "READ", "", 0);
    TFile *f2 = new TFile("output_GetPrediction/bkg.root", "READ", "", 0);
	selection = (TH1F*) f1->FindObjectAny("VBF_MET_presel_4JV_dPhiSide_selection");
	prediction = (TH1F*) f1->FindObjectAny("VBF_MET_presel_4JV_dPhiSide_prediction_px");
	background2 = (TH1F*) f2->FindObjectAny("met_check_inVR_rebin");
	TH1F* background = new TH1F(*prediction);
	background->Reset();
	for (int i = 1; i <= background->GetXaxis()->GetNbins(); ++i) {
		float x = background->GetXaxis()->GetBinCenter(i);
		int j = background2->GetXaxis()->FindBin(x);
		float value = background2->GetBinContent(j);
		float error = background2->GetBinError(j);
		background->SetBinContent(i,value);
		background->SetBinError(i,error);
	}
	
	//double MinX = selection->GetXaxis()->GetBinLowEdge(1);
    //double MaxX = selection->GetXaxis()->GetBinUpEdge(selection->GetXaxis()->GetNbins());
	double MinX = 150;
    double MaxX = 500;
    double BinWidth = selection->GetXaxis()->GetBinWidth(selection->GetXaxis()->GetNbins());
    double MaxY = prediction->GetBinContent(prediction->GetMaximumBin());
    double MaxYsel = selection->GetBinContent(selection->GetMaximumBin());
    if (MaxY < MaxYsel) MaxY = MaxYsel;
    double YRangeMax = 2.*pow(10., int(log10(MaxY))+2);
    double MinY = prediction->GetBinContent(prediction->GetMinimumBin());
    double MinYsel = selection->GetBinContent(selection->GetMinimumBin());
    if (MinY > MinYsel) MinY = MinYsel;
    if (MinY < 0.001) MinY = 0.001;
    double YRangeMin = 0.5*pow(10., int(log10(MinY))-2);
    TString titlePrediction;
    TString titleSelection;
    TString titleBackground;
    TString RatioTitle;
    TString LumiTitle;
    TString Title;
    TString xTitle;
    TString yTitle;

	LumiTitle = "ATLAS internal, L = 36.1 fb^{  -1}, #sqrt{s} = 13 TeV";

    Title = "3 jets, 1.8<#Delta#phi(jj)<2.7, MET>150 GeV, M(jj)>0.6 TeV, p_{T}^{3rd}<50 GeV";
    xTitle = "#slash{E}_{T} (GeV)";
    yTitle = "Events";

    titlePrediction = "Data-driven Pred.";
    titleSelection = "Data";
    titleBackground = "non-QCD background";

    RatioTitle = "(Pred-Data)/Data";

    static Int_t c_LightBrown = TColor::GetColor( "#D9D9CC" );
    static Int_t c_LightGray  = TColor::GetColor( "#DDDDDD" );

    selection->SetAxisRange(MinX, MaxX, "X");
    selection->GetYaxis()->SetRangeUser(YRangeMin, YRangeMax);
    selection->SetMarkerStyle(20);
    selection->SetMarkerSize(0.9);
    selection->SetMarkerColor(kBlack);
    selection->SetXTitle(xTitle);
    selection->SetYTitle(yTitle);

    prediction->SetAxisRange(MinX, MaxX, "X");
    prediction->GetYaxis()->SetRangeUser(YRangeMin, YRangeMax);
    prediction->SetFillColor(c_LightGray);
    prediction->SetTitle("");
    prediction->SetXTitle(xTitle);
    prediction->SetYTitle(yTitle);

    background->SetAxisRange(MinX, MaxX, "X");
    background->GetYaxis()->SetRangeUser(YRangeMin, YRangeMax);
    background->SetTitle("");
    background->SetLineColor(kRed);
    background->SetLineWidth(2);
    background->SetXTitle(xTitle);
    background->SetYTitle(yTitle);

    TCanvas *c = new TCanvas("ca", "Comparison and ratio of two histos", 700, 700);

    TPad *pad1 = new TPad("pad1a", "pad1a", 0, 0.35, 1, 1);
    pad1->SetLogy();
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();

    prediction->DrawCopy("hist");
    selection->Draw("same");
    prediction->SetFillColor(kAzure-3);
    prediction->SetFillStyle(3354);
    prediction->DrawCopy("e2same");
    background->Scale(36.1/32.6);
    background->Draw("same");

    prediction->SetFillStyle(1001);
    //prediction->SetFillColor(c_LightBrown);
    prediction->SetFillColor(c_LightGray);

    //TLegend* leg1 = new TLegend(0.48, 0.63, 0.95, 0.83);
    TLegend* leg1 = new TLegend(0.44, 0.63, 0.91, 0.83);
    leg1->SetFillStyle(0);
    leg1->SetLineStyle(1);
    leg1->SetTextFont(42);
    //leg1->SetTextSize(0.04);
    leg1->SetTextSize(0.045);
    leg1->AddEntry(prediction, titlePrediction, "lf");
    leg1->AddEntry(selection, titleSelection, "lep");
    leg1->AddEntry(background, titleBackground, "l");
    leg1->Draw("same");

    TPaveText* pt = new TPaveText(0.11, 0.98, 0.95, 0.86, "NDC");
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.045);
    pt->AddText(Title);
    pt->AddText(LumiTitle);
    pt->Draw();

    c->cd();
    TPad *pad2 = new TPad("pad2a", "pad2a", 0, 0, 1, 0.35);
    pad2->SetTopMargin(0);
    pad2->Draw();
    pad2->cd();
    TH1F* r = new TH1F(*prediction);
    r->SetTitle("");
    r->SetLabelSize(0.08, "XYZ");
    r->SetLabelOffset(0.01, "XYZ");
    // r->SetTitleSize(0.09, "XYZ");
    r->SetTitleSize(0.125, "XYZ");
    r->SetTitleOffset(0.95, "X");
    r->SetTitleOffset(0.53, "Y");
    // r->SetTitleOffset(0.65, "Y");
    r->SetTickLength(0.05);
    r->SetYTitle(RatioTitle);
    r->SetStats(0);
    r->SetMarkerStyle(20);
    r->SetMarkerSize(0.9);
    r->SetMarkerColor(kBlack);
    r->Reset();
    r->Add(prediction, 1);
    r->Add(background, 1);
    r->Add(selection, -1);
    r->Divide(selection);
    r->SetMaximum(2.2);
    r->SetMinimum(-2.2);
    r->Draw("ep");
    TLine l;
    l.DrawLine(MinX, 0., MaxX+BinWidth, 0.);
    c->cd();
    
    c->SaveAs("compare.pdf");
	
	return 0;
}

int plot1D(vector<TH1F*> h, vector<string> t, vector<string> x, vector<string> y) {

    int N = h.size();

    TString LumiTitle;
    TString Title;
    TString xTitle;
    TString yTitle;

    //LumiTitle = "Simulation, L = 32.9 fb^{  -1}, #sqrt{s} = 13 TeV";
    LumiTitle = "Simulation, #sqrt{s} = 13 TeV";

    for (int i = 0; i < N; ++i) {

        TH1F* selection = h.at(i);
        Title = t.at(i).c_str();
        xTitle = x.at(i).c_str();
        yTitle = y.at(i).c_str();

        static Int_t c_LightBrown = TColor::GetColor( "#D9D9CC" );
        static Int_t c_LightGray  = TColor::GetColor( "#DDDDDD" );

        selection->SetMarkerStyle(20);
        selection->SetMarkerSize(0.9);
        selection->SetMarkerColor(kBlack);
        selection->SetXTitle(xTitle);
        selection->SetYTitle(yTitle);

        TCanvas *c = new TCanvas("c", "c", 800, 600);
        c->cd();
        c->SetLogy();

        selection->Draw();

        TPaveText* pt = new TPaveText(0.1, 0.98, 0.95, 0.87, "NDC");
        pt->SetBorderSize(0);
        pt->SetFillStyle(0);
        pt->SetTextAlign(12);
        pt->SetTextFont(42);
        pt->SetTextSize(0.04);
        pt->AddText(Title);
        pt->AddText(LumiTitle);
        pt->Draw();

        if (N == 0) {
            c->SaveAs("ABCD_mc_1D.pdf");
        } else {
            if ( i == 0 ) {
                c->SaveAs("ABCD_mc_1D.pdf(");
            } else if ( i == N - 1 ) {
                c->SaveAs("ABCD_mc_1D.pdf)");
            } else {
                c->SaveAs("ABCD_mc_1D.pdf");
            }
        }


    }

    return 0;
}

int plot2D(vector<TH2F*> h, vector<string> t, vector<string> x, vector<string> y, vector<string> z) {

    int N = h.size();

    TString LumiTitle;
    TString Title;
    TString xTitle;
    TString yTitle;
    TString zTitle;

    //LumiTitle = "Data, L = 32.9 fb^{  -1}, #sqrt{s} = 13 TeV";
    LumiTitle = "Simulation, #sqrt{s} = 13 TeV";

    for (int i = 0; i < N; ++i) {

        TH2F* selection = h.at(i);
        Title = t.at(i).c_str();
        xTitle = x.at(i).c_str();
        yTitle = y.at(i).c_str();
        zTitle = z.at(i).c_str();

        static Int_t c_LightBrown = TColor::GetColor( "#D9D9CC" );
        static Int_t c_LightGray  = TColor::GetColor( "#DDDDDD" );

        selection->SetMarkerStyle(20);
        selection->SetMarkerSize(0.9);
        selection->SetMarkerColor(kBlack);
        selection->SetXTitle(xTitle);
        selection->SetYTitle(yTitle);
        selection->SetZTitle(zTitle);

        TCanvas *c = new TCanvas("c", "c", 800, 600);
        c->cd();
        c->SetLogz();

        selection->Draw("COLZ");

        TPaveText* pt = new TPaveText(0.1, 0.98, 0.95, 0.87, "NDC");
        pt->SetBorderSize(0);
        pt->SetFillStyle(0);
        pt->SetTextAlign(12);
        pt->SetTextFont(42);
        pt->SetTextSize(0.04);
        pt->AddText(Title);
        pt->AddText(LumiTitle);
        pt->Draw();

        if (N == 0) {
            c->SaveAs("ABCD_mc_2D.pdf");
        } else {
            if ( i == 0 ) {
                c->SaveAs("ABCD_mc_2D.pdf(");
            } else if ( i == N - 1 ) {
                c->SaveAs("ABCD_mc_2D.pdf)");
            } else {
                c->SaveAs("ABCD_mc_2D.pdf");
            }
        }


    }

    return 0;
}

int plotABCDStudies() {

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
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.15);

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
    gStyle->SetTitleOffset(1.5, "X");
    gStyle->SetTitleOffset(1.5, "Y");
    gStyle->SetTitleOffset(1.0, "Z");
    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");

    // For the axis labels:
    gStyle->SetLabelSize(0.04, "XYZ");
    gStyle->SetLabelOffset(0.01, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");

    // For the legend
    gStyle->SetLegendBorderSize(0);

    gROOT->ForceStyle();

    ////////////////////////////////////////

    TFile *f = new TFile("ABCDStudiesOutput_mc.root", "READ", "", 0);
    TFile *f2 = new TFile("ABCDStudiesOutput_signal.root", "READ", "", 0);

    vector<TH1F*> h;
    vector<string> t;
    vector<string> x;
    vector<string> y;

	string SelTitle = "Jet1 p_{T} > 80 GeV, Jet2 p_{T} > 50 GeV & #Delta#eta > 4.8";

    h_Jet1_Pt  =  (TH1F*) f->FindObjectAny("h_Jet1_Pt");
    h.push_back(h_Jet1_Pt);
    t.push_back(SelTitle);
    x.push_back("Jet1 p_{T} (GeV)");
    y.push_back("Events");

    TH1F* h_Jet2_Pt  =  (TH1F*) f->FindObjectAny("h_Jet2_Pt");
    h.push_back(h_Jet2_Pt);
    t.push_back(SelTitle);
    x.push_back("Jet2 p_{T} (GeV)");
    y.push_back("Events");

    TH1F* h_Jet3_Pt  =  (TH1F*) f->FindObjectAny("h_Jet3_Pt");
    h.push_back(h_Jet3_Pt);
    t.push_back(SelTitle);
    x.push_back("Jet3 p_{T} (GeV)");
    y.push_back("Events");

    TH1F* h_Jet1_Eta  =  (TH1F*) f->FindObjectAny("h_Jet1_Eta");
    h.push_back(h_Jet1_Eta);
    t.push_back(SelTitle);
    x.push_back("Jet1 #eta");
    y.push_back("Events");

    TH1F* h_Jet2_Eta  =  (TH1F*) f->FindObjectAny("h_Jet2_Eta");
    h.push_back(h_Jet2_Eta);
    t.push_back(SelTitle);
    x.push_back("Jet2 #eta");
    y.push_back("Events");

    TH1F* h_Jet3_Eta  =  (TH1F*) f->FindObjectAny("h_Jet3_Eta");
    h.push_back(h_Jet3_Eta);
    t.push_back(SelTitle);
    x.push_back("Jet3 #eta");
    y.push_back("Events");

    TH1F* h_Jet1_Phi  =  (TH1F*) f->FindObjectAny("h_Jet1_Phi");
    h.push_back(h_Jet1_Phi);
    t.push_back(SelTitle);
    x.push_back("Jet1 #phi");
    y.push_back("Events");

    TH1F* h_Jet2_Phi  =  (TH1F*) f->FindObjectAny("h_Jet2_Phi");
    h.push_back(h_Jet2_Phi);
    t.push_back(SelTitle);
    x.push_back("Jet2 #phi");
    y.push_back("Events");

    TH1F* h_Jet3_Phi  =  (TH1F*) f->FindObjectAny("h_Jet3_Phi");
    h.push_back(h_Jet3_Phi);
    t.push_back(SelTitle);
    x.push_back("Jet3 #phi");
    y.push_back("Events");

    TH1F* h_Jet1_DeltaPhi  =  (TH1F*) f->FindObjectAny("h_Jet1_DeltaPhi");
    h.push_back(h_Jet1_DeltaPhi);
    t.push_back(SelTitle);
    x.push_back("#Delta#phi(MET,jet1)");
    y.push_back("Events");

    TH1F* h_Jet2_DeltaPhi  =  (TH1F*) f->FindObjectAny("h_Jet2_DeltaPhi");
    h.push_back(h_Jet2_DeltaPhi);
    t.push_back(SelTitle);
    x.push_back("#Delta#phi(MET,jet2)");
    y.push_back("Events");

    TH1F* h_Jet3_DeltaPhi  =  (TH1F*) f->FindObjectAny("h_Jet3_DeltaPhi");
    h.push_back(h_Jet3_DeltaPhi);
    t.push_back(SelTitle);
    x.push_back("#Delta#phi(MET,jet3)");
    y.push_back("Events");
                        
    TH1F* h_MET  =  (TH1F*) f->FindObjectAny("h_MET");
    h.push_back(h_MET);
    t.push_back(SelTitle);
    x.push_back("MET (GeV)");
    y.push_back("Events");

    TH1F* h_MET2  =  (TH1F*) f2->FindObjectAny("h_MET");

	//// overlay plotting MET
	
	h_MET->SetMarkerStyle(20);
	h_MET->SetMarkerSize(0.9);
	h_MET->SetMarkerColor(kBlack);
	h_MET->SetXTitle("MET (GeV)");
	h_MET->SetYTitle("Events");
	//h_DeltaPhijj->SetMinimum(1.);

	h_MET2->SetMarkerStyle(20);
	h_MET2->SetMarkerSize(0.9);
	h_MET2->SetMarkerColor(kRed);
	h_MET2->SetXTitle("MET (GeV)");
	h_MET2->SetYTitle("Events");
	//h_MET2->SetMinimum(1.);

	TCanvas *j = new TCanvas("c", "c", 800, 600);
	j->cd();
	j->SetLogy();

	h_MET->Scale(1./h_MET->Integral());
	h_MET->Draw();
	h_MET2->Scale(1./h_MET2->Integral());
	h_MET2->Draw("same");

	TPaveText* jt = new TPaveText(0.1, 0.98, 0.95, 0.87, "NDC");
	jt->SetBorderSize(0);
	jt->SetFillStyle(0);
	jt->SetTextAlign(12);
	jt->SetTextFont(42);
	jt->SetTextSize(0.04);
	jt->AddText(SelTitle.c_str());
	//jt->AddText("Data, L = 32.9 fb^{  -1}, #sqrt{s} = 13 TeV");
	jt->AddText("Simulation, #sqrt{s} = 13 TeV");
	jt->Draw();

	j->SaveAs("MET.pdf");

	//// end overlay plotting MET

	string SelTitle2 = "Jet1 p_{T} > 80 GeV, Jet2 p_{T} > 50 GeV & #Delta#eta > 3.0";
    TH1F* h_METsig  =  (TH1F*) f->FindObjectAny("h_METsig");
    h.push_back(h_METsig);
    t.push_back(SelTitle2);
    x.push_back("MET/#sqrt{HT} (GeV^{1/2})");
    y.push_back("Events");

    TH1F* h_METsig2  =  (TH1F*) f2->FindObjectAny("h_METsig");
    
	//// overlay plotting MET significance

	static Int_t c_LightBrown = TColor::GetColor( "#D9D9CC" );
	static Int_t c_LightGray  = TColor::GetColor( "#DDDDDD" );

	h_METsig->SetMarkerStyle(20);
	h_METsig->SetMarkerSize(0.9);
	h_METsig->SetMarkerColor(kBlack);
	h_METsig->SetXTitle("MET/#sqrt{HT} (GeV^{1/2})");
	h_METsig->SetYTitle("Events");

	h_METsig2->SetMarkerStyle(20);
	h_METsig2->SetMarkerSize(0.9);
	h_METsig2->SetMarkerColor(kRed);
	h_METsig2->SetXTitle("MET/#sqrt{HT} (GeV^{1/2})");
	h_METsig2->SetYTitle("Events");
	h_METsig2->SetMinimum(1.);

	TCanvas *c = new TCanvas("c", "c", 800, 600);
	c->cd();
	c->SetLogy();

	h_METsig2->Draw();
	h_METsig->Draw("same");

	TPaveText* pt = new TPaveText(0.1, 0.98, 0.95, 0.87, "NDC");
	pt->SetBorderSize(0);
	pt->SetFillStyle(0);
	pt->SetTextAlign(12);
	pt->SetTextFont(42);
	pt->SetTextSize(0.04);
	pt->AddText(SelTitle2.c_str());
	pt->AddText("Data, L = 32.9 fb^{  -1}, #sqrt{s} = 13 TeV");
	pt->Draw();

	c->SaveAs("METsig.pdf");

	//// overlay plotting MET significance

    TH1F* h_METsoft  =  (TH1F*) f->FindObjectAny("h_METsoft");
    h.push_back(h_METsoft);
    t.push_back(SelTitle2);
    x.push_back("MET_{soft} (GeV)");
    y.push_back("Events");

    TH1F* h_METsoft2  =  (TH1F*) f2->FindObjectAny("h_METsoft");

	//// overlay plotting METsoft
	
	h_METsoft->SetMarkerStyle(20);
	h_METsoft->SetMarkerSize(0.9);
	h_METsoft->SetMarkerColor(kBlack);
	h_METsoft->SetXTitle("MET_{soft} (GeV)");
	h_METsoft->SetYTitle("Events");

	h_METsoft2->SetMarkerStyle(20);
	h_METsoft2->SetMarkerSize(0.9);
	h_METsoft2->SetMarkerColor(kRed);
	h_METsoft2->SetXTitle("MET_{soft} (GeV)");
	h_METsoft2->SetYTitle("Events");
	h_METsoft2->SetMinimum(1.);

	TCanvas *d = new TCanvas("c", "c", 800, 600);
	d->cd();
	d->SetLogy();

	h_METsoft2->Draw();
	h_METsoft->Draw("same");

	TPaveText* bt = new TPaveText(0.1, 0.98, 0.95, 0.87, "NDC");
	bt->SetBorderSize(0);
	bt->SetFillStyle(0);
	bt->SetTextAlign(12);
	bt->SetTextFont(42);
	bt->SetTextSize(0.04);
	bt->AddText(SelTitle2.c_str());
	//bt->AddText("Data, L = 32.9 fb^{  -1}, #sqrt{s} = 13 TeV");
	bt->AddText("Simulation, #sqrt{s} = 13 TeV");
	bt->Draw();

	d->SaveAs("METsoft.pdf");

	//// end overlay plotting METsoft

    TH1F* h_DeltaPhijj  =  (TH1F*) f->FindObjectAny("h_DeltaPhijj");
    h.push_back(h_DeltaPhijj);
    t.push_back(SelTitle);
    x.push_back("#Delta#phi(jj)");
    y.push_back("Events");

    TH1F* h_DeltaPhijj2  =  (TH1F*) f2->FindObjectAny("h_DeltaPhijj");

	//// overlay plotting DeltaPhijj
	
	h_DeltaPhijj->SetMarkerStyle(20);
	h_DeltaPhijj->SetMarkerSize(0.9);
	h_DeltaPhijj->SetMarkerColor(kBlack);
	h_DeltaPhijj->SetXTitle("#Delta#phi(jj)");
	h_DeltaPhijj->SetYTitle("Events");
	//h_DeltaPhijj->SetMinimum(1.);

	h_DeltaPhijj2->SetMarkerStyle(20);
	h_DeltaPhijj2->SetMarkerSize(0.9);
	h_DeltaPhijj2->SetMarkerColor(kRed);
	h_DeltaPhijj2->SetXTitle("#Delta#phi(jj)");
	h_DeltaPhijj2->SetYTitle("Events");
	//h_DeltaPhijj2->SetMinimum(1.);

	TCanvas *i = new TCanvas("c", "c", 800, 600);
	i->cd();
	i->SetLogy();

	h_DeltaPhijj->Scale(1./h_DeltaPhijj->Integral());
	h_DeltaPhijj->Draw();
	h_DeltaPhijj2->Scale(1./h_DeltaPhijj2->Integral());
	h_DeltaPhijj2->Draw("same");

	TPaveText* it = new TPaveText(0.1, 0.98, 0.95, 0.87, "NDC");
	it->SetBorderSize(0);
	it->SetFillStyle(0);
	it->SetTextAlign(12);
	it->SetTextFont(42);
	it->SetTextSize(0.04);
	it->AddText(SelTitle.c_str());
	//it->AddText("Data, L = 32.9 fb^{  -1}, #sqrt{s} = 13 TeV");
	it->AddText("Simulation, #sqrt{s} = 13 TeV");
	it->Draw();

	i->SaveAs("dPhijj.pdf");

	//// end overlay plotting DeltaPhijj

    TH1F* h_DeltaEtajj  =  (TH1F*) f->FindObjectAny("h_DeltaEtajj");
    h.push_back(h_DeltaEtajj);
    t.push_back(SelTitle);
    x.push_back("#Delta#eta(jj)");
    y.push_back("Events");

    TH1F* h_DeltaEtajj2  =  (TH1F*) f2->FindObjectAny("h_DeltaEtajj");

	//// overlay plotting dEta
	
	h_DeltaEtajj->SetMarkerStyle(20);
	h_DeltaEtajj->SetMarkerSize(0.9);
	h_DeltaEtajj->SetMarkerColor(kBlack);
	h_DeltaEtajj->SetXTitle("#Delta#eta(jj)");
	h_DeltaEtajj->SetYTitle("Events");
	//h_DeltaEtajj->SetMinimum(1.);

	h_DeltaEtajj2->SetMarkerStyle(20);
	h_DeltaEtajj2->SetMarkerSize(0.9);
	h_DeltaEtajj2->SetMarkerColor(kRed);
	h_DeltaEtajj2->SetXTitle("#Delta#eta(jj)");
	h_DeltaEtajj2->SetYTitle("Events");
	//h_DeltaEtajj2->SetMinimum(1.);

	TCanvas *e = new TCanvas("c", "c", 800, 600);
	e->cd();
	e->SetLogy();

	h_DeltaEtajj->Scale(1./h_DeltaEtajj->Integral());
	h_DeltaEtajj->Draw();
	h_DeltaEtajj2->Scale(1./h_DeltaEtajj2->Integral());
	h_DeltaEtajj2->Draw("same");

	TPaveText* et = new TPaveText(0.1, 0.98, 0.95, 0.87, "NDC");
	et->SetBorderSize(0);
	et->SetFillStyle(0);
	et->SetTextAlign(12);
	et->SetTextFont(42);
	et->SetTextSize(0.04);
	et->AddText(SelTitle.c_str());
	//et->AddText("Data, L = 32.9 fb^{  -1}, #sqrt{s} = 13 TeV");
	et->AddText("Simulation, #sqrt{s} = 13 TeV");
	et->Draw();

	e->SaveAs("DeltaEtajj.pdf");

	//// end overlay plotting dEta

    TH1F* h_Mjj  =  (TH1F*) f->FindObjectAny("h_Mjj");
    h.push_back(h_Mjj);
    t.push_back(SelTitle);
    x.push_back("M(jj) (GeV)");
    y.push_back("Events");

    TH1F* h_Mjj2  =  (TH1F*) f2->FindObjectAny("h_Mjj");

	//// overlay plotting Mjj
	
	h_Mjj->SetMarkerStyle(20);
	h_Mjj->SetMarkerSize(0.9);
	h_Mjj->SetMarkerColor(kBlack);
	h_Mjj->SetXTitle("M(jj) (GeV)");
	h_Mjj->SetYTitle("Events");
	//h_Mjj->SetMinimum(1.);

	h_Mjj2->SetMarkerStyle(20);
	h_Mjj2->SetMarkerSize(0.9);
	h_Mjj2->SetMarkerColor(kRed);
	h_Mjj2->SetXTitle("M(jj) (GeV)");
	h_Mjj2->SetYTitle("Events");
	//h_Mjj2->SetMinimum(1.);

	TCanvas *g = new TCanvas("c", "c", 800, 600);
	g->cd();
	g->SetLogy();

	h_Mjj->Scale(1./h_Mjj->Integral());
	h_Mjj->Draw();
	h_Mjj2->Scale(1./h_Mjj2->Integral());
	h_Mjj2->Draw("same");

	TPaveText* gt = new TPaveText(0.1, 0.98, 0.95, 0.87, "NDC");
	gt->SetBorderSize(0);
	gt->SetFillStyle(0);
	gt->SetTextAlign(12);
	gt->SetTextFont(42);
	gt->SetTextSize(0.04);
	gt->AddText(SelTitle.c_str());
	//gt->AddText("Data, L = 32.9 fb^{  -1}, #sqrt{s} = 13 TeV");
	gt->AddText("Simulation, #sqrt{s} = 13 TeV");
	gt->Draw();

	g->SaveAs("Mjj.pdf");

	//// end overlay plotting Mjj

    plot1D(h,t,x,y);

    vector<TH2F*> h2;
    vector<string> t2;
    vector<string> x2;
    vector<string> y2;
    vector<string> z2;

	SelTitle = "M_{jj} > 600 GeV, #Delta#eta > 4.8";
    TH2F* h_MET_vs_dPhi_Incl  =  (TH2F*) f->FindObjectAny("h_MET_vs_dPhi_Incl");
    h2.push_back(h_MET_vs_dPhi_Incl);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("#Delta#phi(jj)");
    z2.push_back("Events");

	SelTitle = "600 GeV < M_{jj} < 1000 GeV, #Delta#eta > 4.8";
    TH2F* h_MET_vs_dPhi_Mjj600  =  (TH2F*) f->FindObjectAny("h_MET_vs_dPhi_Mjj600");
    h2.push_back(h_MET_vs_dPhi_Mjj600);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("#Delta#phi(jj)");
    z2.push_back("Events");

	SelTitle = "1000 GeV < M_{jj} < 1500 GeV, #Delta#eta > 4.8";
    TH2F* h_MET_vs_dPhi_Mjj1000  =  (TH2F*) f->FindObjectAny("h_MET_vs_dPhi_Mjj1000");
    h2.push_back(h_MET_vs_dPhi_Mjj1000);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("#Delta#phi(jj)");
    z2.push_back("Events");

	SelTitle = "1500 GeV < M_{jj} < 2000 GeV, #Delta#eta > 4.8";
    TH2F* h_MET_vs_dPhi_Mjj1500  =  (TH2F*) f->FindObjectAny("h_MET_vs_dPhi_Mjj1500");
    h2.push_back(h_MET_vs_dPhi_Mjj1500);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("#Delta#phi(jj)");
    z2.push_back("Events");

	SelTitle = "M_{jj} > 2000 GeV, #Delta#eta > 4.8";
    TH2F* h_MET_vs_dPhi_Mjj2000  =  (TH2F*) f->FindObjectAny("h_MET_vs_dPhi_Mjj2000");
    h2.push_back(h_MET_vs_dPhi_Mjj2000);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("#Delta#phi(jj)");
    z2.push_back("Events");

	SelTitle = "M_{jj} > 600 GeV, 1.8 < #Delta#phi < #pi";
    TH2F* h_MET_vs_dEta_Incl  =  (TH2F*) f->FindObjectAny("h_MET_vs_dEta_Incl");
    h2.push_back(h_MET_vs_dEta_Incl);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("#Delta#eta(jj)");
    z2.push_back("Events");

	SelTitle = "600 GeV < M_{jj} < 1000 GeV, 1.8 < #Delta#phi < #pi";
    TH2F* h_MET_vs_dEta_Mjj600  =  (TH2F*) f->FindObjectAny("h_MET_vs_dEta_Mjj600");
    h2.push_back(h_MET_vs_dEta_Mjj600);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("#Delta#eta(jj)");
    z2.push_back("Events");

	SelTitle = "1000 GeV < M_{jj} < 1500 GeV, 1.8 < #Delta#phi < #pi";
    TH2F* h_MET_vs_dEta_Mjj1000  =  (TH2F*) f->FindObjectAny("h_MET_vs_dEta_Mjj1000");
    h2.push_back(h_MET_vs_dEta_Mjj1000);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("#Delta#eta(jj)");
    z2.push_back("Events");

	SelTitle = "1500 GeV < M_{jj} < 2000 GeV, 1.8 < #Delta#phi < #pi";
    TH2F* h_MET_vs_dEta_Mjj1500  =  (TH2F*) f->FindObjectAny("h_MET_vs_dEta_Mjj1500");
    h2.push_back(h_MET_vs_dEta_Mjj1500);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("#Delta#eta(jj)");
    z2.push_back("Events");

	SelTitle = "M_{jj} > 2000 GeV, 1.8 < #Delta#phi < #pi";
    TH2F* h_MET_vs_dEta_Mjj2000  =  (TH2F*) f->FindObjectAny("h_MET_vs_dEta_Mjj2000");
    h2.push_back(h_MET_vs_dEta_Mjj2000);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("#Delta#eta(jj)");
    z2.push_back("Events");

    plot2D(h2,t2,x2,y2,z2);

    return 0;
}


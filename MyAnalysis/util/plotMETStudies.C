int plot1D(vector<TH1F*> h, vector<string> t, vector<string> x, vector<string> y) {

    int N = h.size();

    TString LumiTitle;
    TString Title;
    TString xTitle;
    TString yTitle;

    LumiTitle = "Simulation, L = 30 fb^{  -1}, #sqrt{s} = 13 TeV";

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
            c->SaveAs("c0LL.pdf");
        } else {
            if ( i == 0 ) {
                c->SaveAs("c0LL.pdf(");
            } else if ( i == N - 1 ) {
                c->SaveAs("c0LL.pdf)");
            } else {
                c->SaveAs("c0LL.pdf");
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

    LumiTitle = "Simulation, L = 30 fb^{  -1}, #sqrt{s} = 13 TeV";

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
            c->SaveAs("c0LL2.pdf");
        } else {
            if ( i == 0 ) {
                c->SaveAs("c0LL2.pdf(");
            } else if ( i == N - 1 ) {
                c->SaveAs("c0LL2.pdf)");
            } else {
                c->SaveAs("c0LL2.pdf");
            }
        }


    }

    return 0;
}

int plotMETStudies() {

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

    TFile *f = new TFile("METStudiesOutput_0L.root", "READ", "", 0);

    vector<TH1F*> h;
    vector<string> t;
    vector<string> x;
    vector<string> y;

    h_Jet1_Pt  =  (TH1F*) f->FindObjectAny("h_Jet1_Pt");
    h.push_back(h_Jet1_Pt);
    t.push_back("No selection");
    x.push_back("Jet1 p_{T} (GeV)");
    y.push_back("Events");

    TH1F* h_Jet2_Pt  =  (TH1F*) f->FindObjectAny("h_Jet2_Pt");
    h.push_back(h_Jet2_Pt);
    t.push_back("No selection");
    x.push_back("Jet2 p_{T} (GeV)");
    y.push_back("Events");

    TH1F* h_Jet3_Pt  =  (TH1F*) f->FindObjectAny("h_Jet3_Pt");
    h.push_back(h_Jet3_Pt);
    t.push_back("No selection");
    x.push_back("Jet3 p_{T} (GeV)");
    y.push_back("Events");

    TH1F* h_Jet1_Eta  =  (TH1F*) f->FindObjectAny("h_Jet1_Eta");
    h.push_back(h_Jet1_Eta);
    t.push_back("No selection");
    x.push_back("Jet1 #eta");
    y.push_back("Events");

    TH1F* h_Jet2_Eta  =  (TH1F*) f->FindObjectAny("h_Jet2_Eta");
    h.push_back(h_Jet2_Eta);
    t.push_back("No selection");
    x.push_back("Jet2 #eta");
    y.push_back("Events");

    TH1F* h_Jet3_Eta  =  (TH1F*) f->FindObjectAny("h_Jet3_Eta");
    h.push_back(h_Jet3_Eta);
    t.push_back("No selection");
    x.push_back("Jet3 #eta");
    y.push_back("Events");

    TH1F* h_Jet1_Phi  =  (TH1F*) f->FindObjectAny("h_Jet1_Phi");
    h.push_back(h_Jet1_Phi);
    t.push_back("No selection");
    x.push_back("Jet1 #phi");
    y.push_back("Events");

    TH1F* h_Jet2_Phi  =  (TH1F*) f->FindObjectAny("h_Jet2_Phi");
    h.push_back(h_Jet2_Phi);
    t.push_back("No selection");
    x.push_back("Jet2 #phi");
    y.push_back("Events");

    TH1F* h_Jet3_Phi  =  (TH1F*) f->FindObjectAny("h_Jet3_Phi");
    h.push_back(h_Jet3_Phi);
    t.push_back("No selection");
    x.push_back("Jet3 #phi");
    y.push_back("Events");

    TH1F* h_Jet1_DeltaPhi  =  (TH1F*) f->FindObjectAny("h_Jet1_DeltaPhi");
    h.push_back(h_Jet1_DeltaPhi);
    t.push_back("No selection");
    x.push_back("#Delta#phi(MET,jet1)");
    y.push_back("Events");

    TH1F* h_Jet2_DeltaPhi  =  (TH1F*) f->FindObjectAny("h_Jet2_DeltaPhi");
    h.push_back(h_Jet2_DeltaPhi);
    t.push_back("No selection");
    x.push_back("#Delta#phi(MET,jet2)");
    y.push_back("Events");

    TH1F* h_Jet3_DeltaPhi  =  (TH1F*) f->FindObjectAny("h_Jet3_DeltaPhi");
    h.push_back(h_Jet3_DeltaPhi);
    t.push_back("No selection");
    x.push_back("#Delta#phi(MET,jet3)");
    y.push_back("Events");

    TH1F* h_MHT  =  (TH1F*) f->FindObjectAny("h_MHT");
    h.push_back(h_MHT);
    t.push_back("No selection");
    x.push_back("MHT (GeV)");
    y.push_back("Events");

    TH1F* h_NoJVTMHT  =  (TH1F*) f->FindObjectAny("h_NoJVTMHT");
    h.push_back(h_NoJVTMHT);
    t.push_back("No selection");
    x.push_back("MHT no JVT (GeV)");
    y.push_back("Events");

    TH1F* h_TruthMHT  =  (TH1F*) f->FindObjectAny("h_TruthMHT");
    h.push_back(h_TruthMHT);
    t.push_back("No selection");
    x.push_back("Truth MHT (GeV)");
    y.push_back("Events");

    TH1F* h_MET  =  (TH1F*) f->FindObjectAny("h_MET");
    h.push_back(h_MET);
    t.push_back("No selection");
    x.push_back("MET (GeV)");
    y.push_back("Events");

    TH1F* h_TruthMET  =  (TH1F*) f->FindObjectAny("h_TruthMET");
    h.push_back(h_TruthMET);
    t.push_back("No selection");
    x.push_back("Truth MET (GeV)");
    y.push_back("Events");

    TH1F* h_METreplaced  =  (TH1F*) f->FindObjectAny("h_METreplaced");
    h.push_back(h_METreplaced);
    t.push_back("No selection");
    x.push_back("MET with truth jet term (GeV)");
    y.push_back("Events");

    plot1D(h,t,x,y);

    vector<TH2F*> h2;
    vector<string> t2;
    vector<string> x2;
    vector<string> y2;
    vector<string> z2;

    TH2F* h_MHTvsNoJVTMHT  =  (TH2F*) f->FindObjectAny("h_MHTvsNoJVTMHT");
    h2.push_back(h_MHTvsNoJVTMHT);
    t2.push_back("Jet1 p_{T} > 80 GeV, Jet2 p_{T} > 50 GeV & #Delta#eta > 3.5");
    x2.push_back("MHT (GeV)");
    y2.push_back("MHT no JVT (GeV)");
    z2.push_back("Events");

    TH2F* h_METvsMHT  =  (TH2F*) f->FindObjectAny("h_METvsMHT");
    h2.push_back(h_METvsMHT);
    t2.push_back("Jet1 p_{T} > 80 GeV, Jet2 p_{T} > 50 GeV & #Delta#eta > 3.5");
    x2.push_back("MET (GeV)");
    y2.push_back("MHT (GeV)");
    z2.push_back("Events");

    TH2F* h_METvsMETreplaced  =  (TH2F*) f->FindObjectAny("h_METvsMETreplaced");
    h2.push_back(h_METvsMETreplaced);
    t2.push_back("Jet1 p_{T} > 80 GeV, Jet2 p_{T} > 50 GeV & #Delta#eta > 3.5");
    x2.push_back("MET (GeV)");
    y2.push_back("MET with truth jet term (GeV)");
    z2.push_back("Events");

    TH2F* h_METvsNoJVTMHT  =  (TH2F*) f->FindObjectAny("h_METvsNoJVTMHT");
    h2.push_back(h_METvsNoJVTMHT);
    t2.push_back("Jet1 p_{T} > 80 GeV, Jet2 p_{T} > 50 GeV & #Delta#eta > 3.5");
    x2.push_back("MET (GeV)");
    y2.push_back("MHT no JVT (GeV)");
    z2.push_back("Events");

    TH2F* h_METvsElePt  =  (TH2F*) f->FindObjectAny("h_METvsElePt");
    h2.push_back(h_METvsElePt);
    t2.push_back("Jet1 p_{T} > 80 GeV, Jet2 p_{T} > 50 GeV & #Delta#eta > 3.5");
    x2.push_back("MET (GeV)");
    y2.push_back("Electron p_{T} (GeV)");
    z2.push_back("Events");

    TH2F* h_METvsMuPt  =  (TH2F*) f->FindObjectAny("h_METvsMuPt");
    h2.push_back(h_METvsMuPt);
    t2.push_back("Jet1 p_{T} > 80 GeV, Jet2 p_{T} > 50 GeV & #Delta#eta > 3.5");
    x2.push_back("MET (GeV)");
    y2.push_back("Muon p_{T} (GeV)");
    z2.push_back("Events");

    TH2F* h_METvsPhoPt  =  (TH2F*) f->FindObjectAny("h_METvsPhoPt");
    h2.push_back(h_METvsPhoPt);
    t2.push_back("Jet1 p_{T} > 80 GeV, Jet2 p_{T} > 50 GeV & #Delta#eta > 3.5");
    x2.push_back("MET (GeV)");
    y2.push_back("Photon p_{T} (GeV)");
    z2.push_back("Events");

    TH2F* h_METvsLepPt  =  (TH2F*) f->FindObjectAny("h_METvsLepPt");
    h2.push_back(h_METvsLepPt);
    t2.push_back("Jet1 p_{T} > 80 GeV, Jet2 p_{T} > 50 GeV & #Delta#eta > 3.5");
    x2.push_back("MET (GeV)");
    y2.push_back("Lepton p_{T} (GeV)");
    z2.push_back("Events");

    TH2F* h_METvsMETjet  =  (TH2F*) f->FindObjectAny("h_METvsMETjet");
    h2.push_back(h_METvsMETjet);
    t2.push_back("Jet1 p_{T} > 80 GeV, Jet2 p_{T} > 50 GeV & #Delta#eta > 3.5");
    x2.push_back("MET (GeV)");
    y2.push_back("MET jet term (GeV)");
    z2.push_back("Events");

    TH2F* h_METvsMETele  =  (TH2F*) f->FindObjectAny("h_METvsMETele");
    h2.push_back(h_METvsMETele);
    t2.push_back("Jet1 p_{T} > 80 GeV, Jet2 p_{T} > 50 GeV & #Delta#eta > 3.5");
    x2.push_back("MET (GeV)");
    y2.push_back("MET electron term (GeV)");
    z2.push_back("Events");

    TH2F* h_METvsMETmu  =  (TH2F*) f->FindObjectAny("h_METvsMETmu");
    h2.push_back(h_METvsMETmu);
    t2.push_back("Jet1 p_{T} > 80 GeV, Jet2 p_{T} > 50 GeV & #Delta#eta > 3.5");
    x2.push_back("MET (GeV)");
    y2.push_back("MET muon term (GeV)");
    z2.push_back("Events");

    TH2F* h_METvsMETgamma  =  (TH2F*) f->FindObjectAny("h_METvsMETgamma");
    h2.push_back(h_METvsMETgamma);
    t2.push_back("Jet1 p_{T} > 80 GeV, Jet2 p_{T} > 50 GeV & #Delta#eta > 3.5");
    x2.push_back("MET (GeV)");
    y2.push_back("MET photon term (GeV)");
    z2.push_back("Events");

    TH2F* h_METvsMETtrack  =  (TH2F*) f->FindObjectAny("h_METvsMETtrack");
    h2.push_back(h_METvsMETtrack);
    t2.push_back("Jet1 p_{T} > 80 GeV, Jet2 p_{T} > 50 GeV & #Delta#eta > 3.5");
    x2.push_back("MET (GeV)");
    y2.push_back("MET soft term (GeV)");
    z2.push_back("Events");

    plot2D(h2,t2,x2,y2,z2);

    return 0;
}


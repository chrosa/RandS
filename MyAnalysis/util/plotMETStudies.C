int plot1D(vector<TH1F*> h, vector<string> t, vector<string> x, vector<string> y) {

    int N = h.size();

    TString LumiTitle;
    TString Title;
    TString xTitle;
    TString yTitle;

    LumiTitle = "Simulation, L = 32.9 fb^{  -1}, #sqrt{s} = 13 TeV";

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
            c->SaveAs("c0LAddMu.pdf");
        } else {
            if ( i == 0 ) {
                c->SaveAs("c0LAddMu.pdf(");
            } else if ( i == N - 1 ) {
                c->SaveAs("c0LAddMu.pdf)");
            } else {
                c->SaveAs("c0LAddMu.pdf");
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

    LumiTitle = "Simulation, L = 32.9 fb^{  -1}, #sqrt{s} = 13 TeV";

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
            c->SaveAs("c0LAddMu2.pdf");
        } else {
            if ( i == 0 ) {
                c->SaveAs("c0LAddMu2.pdf(");
            } else if ( i == N - 1 ) {
                c->SaveAs("c0LAddMu2.pdf)");
            } else {
                c->SaveAs("c0LAddMu2.pdf");
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

	//string SelTitle = "No selection";
	//string SelTitle = "Baseline lepton veto; added all muons to MHT";
	string SelTitle = "Signal lepton veto; added all muons to MHT";
	//string SelTitle = "Baseline lepton veto; added passOR muons to MHT";
	//string SelTitle = "Signal lepton veto; added passOR muons to MHT";
	//string SelTitle = "Jet1 p_{T} > 80 GeV, Jet2 p_{T} > 50 GeV & #Delta#eta > 3.5";

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
                        
    TH1F* h_MHT  =  (TH1F*) f->FindObjectAny("h_MHT");
    h.push_back(h_MHT);
    t.push_back(SelTitle);
    x.push_back("MHT (GeV)");
    y.push_back("Events");

    TH1F* h_MHTnoOR  =  (TH1F*) f->FindObjectAny("h_MHTnoOR");
    h.push_back(h_MHTnoOR);
    t.push_back(SelTitle);
    x.push_back("MHT no OR (GeV)");
    y.push_back("Events");

    TH1F* h_MHTnoJVT  =  (TH1F*) f->FindObjectAny("h_MHTnoJVT");
    h.push_back(h_MHTnoJVT);
    t.push_back(SelTitle);
    x.push_back("MHT no JVT (GeV)");
    y.push_back("Events");

    TH1F* h_MHTnoJVTnoOR  =  (TH1F*) f->FindObjectAny("h_MHTnoJVTnoOR");
    h.push_back(h_MHTnoJVTnoOR);
    t.push_back(SelTitle);
    x.push_back("MHT no JVT no OR (GeV)");
    y.push_back("Events");

    TH1F* h_TruthMHT  =  (TH1F*) f->FindObjectAny("h_TruthMHT");
    h.push_back(h_TruthMHT);
    t.push_back(SelTitle);
    x.push_back("Truth MHT (GeV)");
    y.push_back("Events");

    TH1F* h_MET  =  (TH1F*) f->FindObjectAny("h_MET");
    h.push_back(h_MET);
    t.push_back(SelTitle);
    x.push_back("MET (GeV)");
    y.push_back("Events");

    TH1F* h_TruthMET  =  (TH1F*) f->FindObjectAny("h_TruthMET");
    h.push_back(h_TruthMET);
    t.push_back(SelTitle);
    x.push_back("Truth MET (GeV)");
    y.push_back("Events");

    TH1F* h_METreplaced  =  (TH1F*) f->FindObjectAny("h_METreplaced");
    h.push_back(h_METreplaced);
    t.push_back(SelTitle);
    x.push_back("MET with truth jet term (GeV)");
    y.push_back("Events");

    plot1D(h,t,x,y);

    vector<TH2F*> h2;
    vector<string> t2;
    vector<string> x2;
    vector<string> y2;
    vector<string> z2;


    TH2F* h_MET_vs_MHT  =  (TH2F*) f->FindObjectAny("h_MET_vs_MHT");
    h2.push_back(h_MET_vs_MHT);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("MHT (GeV)");
    z2.push_back("Events");

    TH2F* h_MET_vs_MHTnoOR  =  (TH2F*) f->FindObjectAny("h_MET_vs_MHTnoOR");
    h2.push_back(h_MET_vs_MHTnoOR);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("MHT no OR (GeV)");
    z2.push_back("Events");

    TH2F* h_MET_vs_MHTnoJVT  =  (TH2F*) f->FindObjectAny("h_MET_vs_MHTnoJVT");
    h2.push_back(h_MET_vs_MHTnoJVT);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("MHT no JVT (GeV)");
    z2.push_back("Events");

    TH2F* h_MET_vs_MHTnoJVTnoOR  =  (TH2F*) f->FindObjectAny("h_MET_vs_MHTnoJVTnoOR");
    h2.push_back(h_MET_vs_MHTnoJVTnoOR);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("MHT no JVT no OR (GeV)");
    z2.push_back("Events");

    TH2F* h_MHT_vs_MHTnoOR  =  (TH2F*) f->FindObjectAny("h_MHT_vs_MHTnoOR");
    h2.push_back(h_MHT_vs_MHTnoOR);
    t2.push_back(SelTitle);
    x2.push_back("MHT (GeV)");
    y2.push_back("MHT noOR (GeV)");
    z2.push_back("Events");

    TH2F* h_MHT_vs_MHTnoJVT  =  (TH2F*) f->FindObjectAny("h_MHT_vs_MHTnoJVT");
    h2.push_back(h_MHT_vs_MHTnoJVT);
    t2.push_back(SelTitle);
    x2.push_back("MHT (GeV)");
    y2.push_back("MHT no JVT (GeV)");
    z2.push_back("Events");

    TH2F* h_MHT_vs_MHTnoJVTnoOR  =  (TH2F*) f->FindObjectAny("h_MHT_vs_MHTnoJVTnoOR");
    h2.push_back(h_MHT_vs_MHTnoJVTnoOR);
    t2.push_back(SelTitle);
    x2.push_back("MHT (GeV)");
    y2.push_back("MHT no JVT no OR (GeV)");
    z2.push_back("Events");

    TH2F* h_JVTPhi_vs_METPhi  =  (TH2F*) f->FindObjectAny("h_JVTPhi_vs_METPhi");
    h2.push_back(h_JVTPhi_vs_METPhi);
    t2.push_back(SelTitle);
    x2.push_back("#Sum JVT jets #Phi");
    y2.push_back("MET #Phi");
    z2.push_back("Events");

    TH2F* h_JVTPhi_vs_MHTPhi  =  (TH2F*) f->FindObjectAny("h_JVTPhi_vs_MHTPhi");
    h2.push_back(h_JVTPhi_vs_MHTPhi);
    t2.push_back(SelTitle);
    x2.push_back("#Sum JVT jets #Phi");
    y2.push_back("MHT #Phi");
    z2.push_back("Events");

    TH2F* h_JVTPhi_vs_MHTnoORPhi  =  (TH2F*) f->FindObjectAny("h_JVTPhi_vs_MHTnoORPhi");
    h2.push_back(h_JVTPhi_vs_MHTnoORPhi);
    t2.push_back(SelTitle);
    x2.push_back("#Sum JVT jets #Phi");
    y2.push_back("MHT no OR #Phi");
    z2.push_back("Events");

    TH2F* h_MET_vs_METreplaced  =  (TH2F*) f->FindObjectAny("h_MET_vs_METreplaced");
    h2.push_back(h_MET_vs_METreplaced);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("MET with truth jet term (GeV)");
    z2.push_back("Events");

    TH2F* h_MHTnoOR_vs_MHTnoORreplaced  =  (TH2F*) f->FindObjectAny("h_MHTnoOR_vs_MHTnoORreplaced");
    h2.push_back(h_MHTnoOR_vs_MHTnoORreplaced);
    t2.push_back(SelTitle);
    x2.push_back("MHT no OR (GeV)");
    y2.push_back("MHT no OR with truth jet term (GeV)");
    z2.push_back("Events");

    TH2F* h_MHTnoJVTnoOR_vs_MHTnoJVTnoORreplaced  =  (TH2F*) f->FindObjectAny("h_MHTnoJVTnoOR_vs_MHTnoJVTnoORreplaced");
    h2.push_back(h_MHTnoJVTnoOR_vs_MHTnoJVTnoORreplaced);
    t2.push_back(SelTitle);
    x2.push_back("MHT no JVT no OR (GeV)");
    y2.push_back("MET no JVT no OR with truth jet term (GeV)");
    z2.push_back("Events");
   
    TH2F* h_MET_vs_ElePt  =  (TH2F*) f->FindObjectAny("h_MET_vs_ElePt");
    h2.push_back(h_MET_vs_ElePt);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("Electron p_{T} (GeV)");
    z2.push_back("Events");

    TH2F* h_MET_vs_MuPt  =  (TH2F*) f->FindObjectAny("h_MET_vs_MuPt");
    h2.push_back(h_MET_vs_MuPt);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("Muon p_{T} (GeV)");
    z2.push_back("Events");

    TH2F* h_MET_vs_PhoPt  =  (TH2F*) f->FindObjectAny("h_MET_vs_PhoPt");
    h2.push_back(h_MET_vs_PhoPt);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("Photon p_{T} (GeV)");
    z2.push_back("Events");

    TH2F* h_MET_vs_LepPt  =  (TH2F*) f->FindObjectAny("h_MET_vs_LepPt");
    h2.push_back(h_MET_vs_LepPt);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("Lepton p_{T} (GeV)");
    z2.push_back("Events");

    TH2F* h_MET_vs_METjet  =  (TH2F*) f->FindObjectAny("h_MET_vs_METjet");
    h2.push_back(h_MET_vs_METjet);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("MET jet term (GeV)");
    z2.push_back("Events");

    TH2F* h_MET_vs_METele  =  (TH2F*) f->FindObjectAny("h_MET_vs_METele");
    h2.push_back(h_MET_vs_METele);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("MET electron term (GeV)");
    z2.push_back("Events");

    TH2F* h_MET_vs_METmu  =  (TH2F*) f->FindObjectAny("h_MET_vs_METmu");
    h2.push_back(h_MET_vs_METmu);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("MET muon term (GeV)");
    z2.push_back("Events");

    TH2F* h_MET_vs_METgamma  =  (TH2F*) f->FindObjectAny("h_MET_vs_METgamma");
    h2.push_back(h_MET_vs_METgamma);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("MET photon term (GeV)");
    z2.push_back("Events");

    TH2F* h_MET_vs_METtrack  =  (TH2F*) f->FindObjectAny("h_MET_vs_METtrack");
    h2.push_back(h_MET_vs_METtrack);
    t2.push_back(SelTitle);
    x2.push_back("MET (GeV)");
    y2.push_back("MET soft term (GeV)");
    z2.push_back("Events");

    plot2D(h2,t2,x2,y2,z2);

    return 0;
}


int plot1D(vector<TH1F*> h, vector<string> t, vector<string> x, vector<string> y) {

    int N = h.size();

    TString LumiTitle;
    TString Title;
    TString xTitle;
    TString yTitle;

    LumiTitle = "ATLAS internal, L = 36.1 fb^{  -1}, #sqrt{s} = 13 TeV";

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

        if (N == 1) {
            c->SaveAs("TTO_data_incl_1D.pdf");
        } else {
            if ( i == 0 ) {
                c->SaveAs("TTO_data_incl_1D.pdf(");
            } else if ( i == N - 1 ) {
                c->SaveAs("TTO_data_incl_1D.pdf)");
            } else {
                c->SaveAs("TTO_data_incl_1D.pdf");
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

    LumiTitle = "ATLAS internal, L = 36.1 fb^{  -1}, #sqrt{s} = 13 TeV";

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

        if (N == 1) {
            c->SaveAs("TTO_data_incl_2D.pdf");
        } else {
            if ( i == 0 ) {
                c->SaveAs("TTO_data_incl_2D.pdf(");
            } else if ( i == N - 1 ) {
                c->SaveAs("TTO_data_incl_2D.pdf)");
            } else {
                c->SaveAs("TTO_data_incl_2D.pdf");
            }
        }


    }

    return 0;
}

int plotTriggerStudies() {

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

    TFile *f = new TFile("TriggerStudiesOutput_incl_data_v3.root", "READ", "", 0);

    vector<TH1F*> h;
    vector<string> t;
    vector<string> x;
    vector<string> y;

    h_MHT_all  =  (TH1F*) f->FindObjectAny("h_MHT_all");
    h.push_back(h_MHT_all);
    t.push_back("single jet triggers");
    x.push_back("MHT (GeV)");
    y.push_back("Events");

    h_MHT_triggered  =  (TH1F*) f->FindObjectAny("h_MHT_triggered");
    h.push_back(h_MHT_triggered);
    t.push_back("single jet triggers & xe70/xe90/xe110 trigger");
    x.push_back("MHT (GeV)");
    y.push_back("Events");

    TH1F* h_MHT_eff = new TH1F(*h_MHT_triggered);
    h_MHT_eff->Divide(h_MHT_all);
    h.push_back(h_MHT_eff);
    t.push_back("xe70/xe90/xe110 trigger");
    x.push_back("MHT (GeV)");
    y.push_back("Efficiency");
    
    //===

    h_MET_all  =  (TH1F*) f->FindObjectAny("h_MET_all");
    h.push_back(h_MET_all);
    t.push_back("single jet triggers");
    x.push_back("MET (GeV)");
    y.push_back("Events");

    h_MET_triggered  =  (TH1F*) f->FindObjectAny("h_MET_triggered");
    h.push_back(h_MET_triggered);
    t.push_back("single jet triggers & xe70/xe90/xe110 trigger");
    x.push_back("MET (GeV)");
    y.push_back("Events");

    TH1F* h_MET_eff = new TH1F(*h_MET_triggered);
    h_MET_eff->Divide(h_MET_all);
    h.push_back(h_MET_eff);
    t.push_back("xe70/xe90/xe110 trigger");
    x.push_back("MET (GeV)");
    y.push_back("Efficiency");
    
    //===

    h_MHT2jet_all  =  (TH1F*) f->FindObjectAny("h_MHT2jet_all");
    h.push_back(h_MHT2jet_all);
    t.push_back("single jet triggers, 2 jets");
    x.push_back("MHT (GeV)");
    y.push_back("Events");

    h_MHT2jet_triggered  =  (TH1F*) f->FindObjectAny("h_MHT2jet_triggered");
    h.push_back(h_MHT2jet_triggered);
    t.push_back("single jet triggers & xe70/xe90/xe110 trigger, 2 jets");
    x.push_back("MHT (GeV)");
    y.push_back("Events");

    TH1F* h_MHT2jet_eff = new TH1F(*h_MHT2jet_triggered);
    h_MHT2jet_eff->Divide(h_MHT2jet_all);
    h.push_back(h_MHT2jet_eff);
    t.push_back("xe70/xe90/xe110 trigger, 2 jets");
    x.push_back("MHT (GeV)");
    y.push_back("Efficiency");
    
    //===

    h_MET2jet_all  =  (TH1F*) f->FindObjectAny("h_MET2jet_all");
    h.push_back(h_MET2jet_all);
    t.push_back("single jet triggers, 2 jets");
    x.push_back("MET (GeV)");
    y.push_back("Events");

    h_MET2jet_triggered  =  (TH1F*) f->FindObjectAny("h_MET2jet_triggered");
    h.push_back(h_MET2jet_triggered);
    t.push_back("single jet triggers & xe70/xe90/xe110 trigger, 2 jets");
    x.push_back("MET (GeV)");
    y.push_back("Events");

    TH1F* h_MET2jet_eff = new TH1F(*h_MET2jet_triggered);
    h_MET2jet_eff->Divide(h_MET2jet_all);
    h.push_back(h_MET2jet_eff);
    t.push_back("xe70/xe90/xe110 trigger, 2 jets");
    x.push_back("MET (GeV)");
    y.push_back("Efficiency");


    plot1D(h,t,x,y);

    vector<TH2F*> h2;
    vector<string> t2;
    vector<string> x2;
    vector<string> y2;
    vector<string> z2;

    TH2F* h_MHTvsHT_all  =  (TH2F*) f->FindObjectAny("h_MHTvsHT_all");
    h2.push_back(h_MHTvsHT_all);
    t2.push_back("single jet triggers");
    x2.push_back("MHT (GeV)");
    y2.push_back("HT (GeV)");
    z2.push_back("Events");

    TH2F* h_MHTvsHT_triggered  =  (TH2F*) f->FindObjectAny("h_MHTvsHT_triggered");
    h2.push_back(h_MHTvsHT_triggered);
    t2.push_back("single jet triggers  & xe70/xe90/xe110 trigger");
    x2.push_back("MHT (GeV)");
    y2.push_back("HT (GeV)");
    z2.push_back("Events");

    TH2F* h_MHTvsHT_eff = new TH2F(*h_MHTvsHT_triggered);
    h_MHTvsHT_eff->Divide(h_MHTvsHT_all);
    h2.push_back(h_MHTvsHT_eff);
    t2.push_back("xe70/xe90/xe110 trigger");
    x2.push_back("MHT (GeV)");
    y2.push_back("HT (GeV)");
    z2.push_back("Efficiency");
    
    //===

    TH2F* h_METvsHT_all  =  (TH2F*) f->FindObjectAny("h_METvsHT_all");
    h2.push_back(h_METvsHT_all);
    t2.push_back("single jet triggers");
    x2.push_back("MET (GeV)");
    y2.push_back("HT (GeV)");
    z2.push_back("Events");

    TH2F* h_METvsHT_triggered  =  (TH2F*) f->FindObjectAny("h_METvsHT_triggered");
    h2.push_back(h_METvsHT_triggered);
    t2.push_back("single jet triggers  & xe70/xe90/xe110 trigger");
    x2.push_back("MET (GeV)");
    y2.push_back("HT (GeV)");
    z2.push_back("Events");

    TH2F* h_METvsHT_eff = new TH2F(*h_METvsHT_triggered);
    h_METvsHT_eff->Divide(h_METvsHT_all);
    h2.push_back(h_METvsHT_eff);
    t2.push_back("xe70/xe90/xe110 trigger");
    x2.push_back("MET (GeV)");
    y2.push_back("HT (GeV)");
    z2.push_back("Efficiency");

    //===
    
    TH2F* h_MHT2jetvsHT_all  =  (TH2F*) f->FindObjectAny("h_MHT2jetvsHT_all");
    h2.push_back(h_MHT2jetvsHT_all);
    t2.push_back("single jet triggers, 2 jets");
    x2.push_back("MHT (GeV)");
    y2.push_back("HT (GeV)");
    z2.push_back("Events");

    TH2F* h_MHT2jetvsHT_triggered  =  (TH2F*) f->FindObjectAny("h_MHT2jetvsHT_triggered");
    h2.push_back(h_MHT2jetvsHT_triggered);
    t2.push_back("single jet triggers  & xe70/xe90/xe110 trigger, 2 jets");
    x2.push_back("MHT (GeV)");
    y2.push_back("HT (GeV)");
    z2.push_back("Events");

    TH2F* h_MHT2jetvsHT_eff = new TH2F(*h_MHT2jetvsHT_triggered);
    h_MHT2jetvsHT_eff->Divide(h_MHT2jetvsHT_all);
    h2.push_back(h_MHT2jetvsHT_eff);
    t2.push_back("xe70/xe90/xe110 trigger, 2 jets");
    x2.push_back("MHT (GeV)");
    y2.push_back("HT (GeV)");
    z2.push_back("Efficiency");

	//===

    TH2F* h_MET2jetvsHT_all  =  (TH2F*) f->FindObjectAny("h_MET2jetvsHT_all");
    h2.push_back(h_MET2jetvsHT_all);
    t2.push_back("single jet triggers, 2 jets");
    x2.push_back("MET (GeV)");
    y2.push_back("HT (GeV)");
    z2.push_back("Events");

    TH2F* h_MET2jetvsHT_triggered  =  (TH2F*) f->FindObjectAny("h_MET2jetvsHT_triggered");
    h2.push_back(h_MET2jetvsHT_triggered);
    t2.push_back("single jet triggers  & xe70/xe90/xe110 trigger, 2 jets");
    x2.push_back("MET (GeV)");
    y2.push_back("HT (GeV)");
    z2.push_back("Events");

    TH2F* h_MET2jetvsHT_eff = new TH2F(*h_MET2jetvsHT_triggered);
    h_MET2jetvsHT_eff->Divide(h_MET2jetvsHT_all);
    h2.push_back(h_MET2jetvsHT_eff);
    t2.push_back("xe70/xe90/xe110 trigger, 2 jets");
    x2.push_back("MET (GeV)");
    y2.push_back("HT (GeV)");
    z2.push_back("Efficiency");

    plot2D(h2,t2,x2,y2,z2);

    return 0;
}


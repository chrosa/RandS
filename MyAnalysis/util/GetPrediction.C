#include <TROOT.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TChain.h>
#include <TPad.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TFile.h>
#include <TPostScript.h>
#include <TLegend.h>
#include <TMath.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "Prediction.h"

using namespace std;

TCanvas* DrawComparison(TH1F* prediction, TH1F* selection, TString Title, TString LumiTitle, TString xTitle, TString yTitle, bool isData)
{
    double MinX = selection->GetXaxis()->GetBinLowEdge(1);
    double MaxX = selection->GetXaxis()->GetBinUpEdge(selection->GetXaxis()->GetNbins());
    double BinWidth = selection->GetXaxis()->GetBinWidth(selection->GetXaxis()->GetNbins());
    double MaxY = prediction->GetBinContent(prediction->GetMaximumBin());
    double MaxYsel = selection->GetBinContent(selection->GetMaximumBin());
    if (MaxY < MaxYsel) MaxY = MaxYsel;
    double YRangeMax = 2.*pow(10., int(log10(MaxY))+3);
    double MinY = prediction->GetBinContent(prediction->GetMinimumBin());
    double MinYsel = selection->GetBinContent(selection->GetMinimumBin());
    if (MinY > MinYsel) MinY = MinYsel;
    if (MinY < 0.1) MinY = 0.1;
    double YRangeMin = 0.5*pow(10., int(log10(MinY))-0);
    TString titlePrediction;
    TString titleSelection;
    TString RatioTitle;

    if( isData ) {
        titlePrediction = "Pred. from Data";
        titleSelection = "Data";
        RatioTitle = "(Pred-Data)/Data";
    }
    else {
        titlePrediction = "Data-driven Pred. from MC";
        //titlePrediction = "Smeared Generator Jets";
        titleSelection = "MC Expectation";
        RatioTitle = "(Pred-MC)/MC";
        //RatioTitle = "(Gen-MC)/MC";
    }

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
    //prediction->SetFillColor(c_LightBrown);
    prediction->SetFillColor(c_LightGray);
    prediction->SetTitle("");
    prediction->SetXTitle(xTitle);
    prediction->SetYTitle(yTitle);

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
    r->Add(selection, -1);
    r->Divide(selection);
    r->SetMaximum(2.2);
    r->SetMinimum(-2.2);
    r->Draw("ep");
    TLine l;
    l.DrawLine(MinX, 0., MaxX+BinWidth, 0.);
    c->cd();
    return c;
}

int main()
{

    int debug = 0;

    // Somehow this does not work!!!
    //gROOT->Reset();
    gROOT->SetStyle("Plain");
    //gStyle->SetPalette(51, 0);
    //gStyle->SetHatchesLineWidth(1.2);

    // For the canvas:
    gStyle->SetCanvasColor(0);
    //gStyle->SetCanvasBorderMode(0);

    // For the Pad:
    gStyle->SetPadColor(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadBorderSize(2);
    //gStyle->SetPadBorderMode(0);

    // For the frame:
    gStyle->SetFrameBorderMode(0);

    // For the histo:
    // gStyle->SetMarkerSize(0.7);
    // gStyle->SetMarkerStyle(20);
    // gStyle->SetMarkerColor(1);

    // For the statistics box:
    gStyle->SetOptStat(0);
    //gStyle->SetOptFit(1011);

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
    //gStyle->SetTitleOffset(1.25, "Y");
    gStyle->SetTitleOffset(1.2, "Y");
    gStyle->SetTitleOffset(0.5, "Z");
    // gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetTitleSize(0.061, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    //gStyle->SetTitleX(0.15);
    //gStyle->SetTitleY(0.99);

    // For the axis labels:
    gStyle->SetLabelSize(0.04, "XYZ");
    gStyle->SetLabelOffset(0.01, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");

    // For the legend
    gStyle->SetLegendBorderSize(0);

    gROOT->ForceStyle();

    ////////////////////////////////////////

    string root_file;
    TChain* prediction = new TChain("PredictionTree");

    ifstream myfile1 ("filelist_RnS_data_all.txt");
    //ifstream myfile1 ("filelist_RnS_mc_all.txt");
    //ifstream myfile1 ("filelist_RnS.txt");

    if (myfile1.is_open()) {
        while( myfile1.good() ) {
            getline (myfile1,root_file);
            cout << root_file << endl;
            if (root_file.length() > 0 ) {
                TString path = root_file;
                prediction->Add(path);
            }

        }
        myfile1.close();
    }

    // ------------------------------------------------------------------- //

    // initialize new Prediction object
    Prediction *pred_;
    bool isData = true;
    bool VBF = true;
    bool VBFSR = true;
    bool HTMHT = false;
    TString postfix = "_MyTest_data_METsoftSmeared_muRes_noAngSmear_N20_SRtw_v1"; //CRlm: MET < 120; SRtw with trigger weight
    //TString postfix = "_test";

    pred_ = new Prediction(*prediction, postfix);

    cout << "after prediction in main" << endl;

    TString LumiTitle;
    if( isData ) LumiTitle = "ATLAS internal, L = 36.1 fb^{  -1}, #sqrt{s} = 13 TeV";
    else LumiTitle = "Simulation, L = 36.1 fb^{  -1}, #sqrt{s} = 13 TeV";

    vector<TString> xTitle_presel;
    xTitle_presel.push_back("H_{T} (GeV)");
    xTitle_presel.push_back("#slash{H}_{T} (GeV)");
    xTitle_presel.push_back("#slash{E}_{T} (GeV)");
    xTitle_presel.push_back("N_{Jets}");
    xTitle_presel.push_back("N_{b-Tags}");
    xTitle_presel.push_back("Jet1 p_{T} (GeV)");
    xTitle_presel.push_back("Jet2 p_{T} (GeV)");
    xTitle_presel.push_back("Jet1 #eta");
    xTitle_presel.push_back("Jet2 #eta");

    vector<TString> xTitle_deltaPhi;
    xTitle_deltaPhi.push_back("H_{T} (GeV)");
    xTitle_deltaPhi.push_back("#slash{H}_{T} (GeV)");
    xTitle_deltaPhi.push_back("#slash{E}_{T} (GeV)");
    xTitle_deltaPhi.push_back("Jet1 p_{T} (GeV)");
    xTitle_deltaPhi.push_back("Jet2 p_{T} (GeV)");
    xTitle_deltaPhi.push_back("Jet1 #eta");
    xTitle_deltaPhi.push_back("Jet2 #eta");

    vector<TString> xTitle_baseline_Bin1;
    xTitle_baseline_Bin1.push_back("Jet1 p_{T} (GeV)");
    xTitle_baseline_Bin1.push_back("Jet2 p_{T} (GeV)");
    xTitle_baseline_Bin1.push_back("Jet1 #eta");
    xTitle_baseline_Bin1.push_back("Jet2 #eta");
    xTitle_baseline_Bin1.push_back("#Delta#phi 1");
    xTitle_baseline_Bin1.push_back("#Delta#phi 2");

    vector<TString> xTitle_baseline_Bin2;
    xTitle_baseline_Bin2.push_back("Jet1 p_{T} (GeV)");
    xTitle_baseline_Bin2.push_back("Jet2 p_{T} (GeV)");
    xTitle_baseline_Bin2.push_back("Jet3 p_{T} (GeV)");
    xTitle_baseline_Bin2.push_back("Jet1 #eta");
    xTitle_baseline_Bin2.push_back("Jet2 #eta");
    xTitle_baseline_Bin2.push_back("Jet3 #eta");
    xTitle_baseline_Bin2.push_back("#Delta#phi 1");
    xTitle_baseline_Bin2.push_back("#Delta#phi 2");
    xTitle_baseline_Bin2.push_back("#Delta#phi 3");

    vector<TString> xTitle_baseline_Bin3;
    xTitle_baseline_Bin3.push_back("Jet1 p_{T} (GeV)");
    xTitle_baseline_Bin3.push_back("Jet2 p_{T} (GeV)");
    xTitle_baseline_Bin3.push_back("Jet3 p_{T} (GeV)");
    xTitle_baseline_Bin3.push_back("Jet1 #eta");
    xTitle_baseline_Bin3.push_back("Jet2 #eta");
    xTitle_baseline_Bin3.push_back("Jet3 #eta");
    xTitle_baseline_Bin3.push_back("#Delta#phi 1");
    xTitle_baseline_Bin3.push_back("#Delta#phi 2");
    xTitle_baseline_Bin3.push_back("#Delta#phi 3");

    vector<TString> xTitle_baseline_Bin4;
    xTitle_baseline_Bin4.push_back("Jet1 p_{T} (GeV)");
    xTitle_baseline_Bin4.push_back("Jet2 p_{T} (GeV)");
    xTitle_baseline_Bin4.push_back("Jet3 p_{T} (GeV)");
    xTitle_baseline_Bin4.push_back("Jet1 #eta");
    xTitle_baseline_Bin4.push_back("Jet2 #eta");
    xTitle_baseline_Bin4.push_back("Jet3 #eta");
    xTitle_baseline_Bin4.push_back("#Delta#phi 1");
    xTitle_baseline_Bin4.push_back("#Delta#phi 2");
    xTitle_baseline_Bin4.push_back("#Delta#phi 3");

    vector<TString> xTitle_baseline_withoutDeltaPhi_Bin1;
    xTitle_baseline_withoutDeltaPhi_Bin1.push_back("Jet1 p_{T} (GeV)");
    xTitle_baseline_withoutDeltaPhi_Bin1.push_back("Jet2 p_{T} (GeV)");
    xTitle_baseline_withoutDeltaPhi_Bin1.push_back("Jet1 #eta");
    xTitle_baseline_withoutDeltaPhi_Bin1.push_back("Jet2 #eta");
    xTitle_baseline_withoutDeltaPhi_Bin1.push_back("#Delta#phi 1");
    xTitle_baseline_withoutDeltaPhi_Bin1.push_back("#Delta#phi 2");

    vector<TString> xTitle_baseline_withoutDeltaPhi_Bin2;
    xTitle_baseline_withoutDeltaPhi_Bin2.push_back("Jet1 p_{T} (GeV)");
    xTitle_baseline_withoutDeltaPhi_Bin2.push_back("Jet2 p_{T} (GeV)");
    xTitle_baseline_withoutDeltaPhi_Bin2.push_back("Jet3 p_{T} (GeV)");
    xTitle_baseline_withoutDeltaPhi_Bin2.push_back("Jet1 #eta");
    xTitle_baseline_withoutDeltaPhi_Bin2.push_back("Jet2 #eta");
    xTitle_baseline_withoutDeltaPhi_Bin2.push_back("Jet3 #eta");
    xTitle_baseline_withoutDeltaPhi_Bin2.push_back("#Delta#phi 1");
    xTitle_baseline_withoutDeltaPhi_Bin2.push_back("#Delta#phi 2");
    xTitle_baseline_withoutDeltaPhi_Bin2.push_back("#Delta#phi 3");

    vector<TString> xTitle_baseline_withoutDeltaPhi_Bin3;
    xTitle_baseline_withoutDeltaPhi_Bin3.push_back("Jet1 p_{T} (GeV)");
    xTitle_baseline_withoutDeltaPhi_Bin3.push_back("Jet2 p_{T} (GeV)");
    xTitle_baseline_withoutDeltaPhi_Bin3.push_back("Jet3 p_{T} (GeV)");
    xTitle_baseline_withoutDeltaPhi_Bin3.push_back("Jet1 #eta");
    xTitle_baseline_withoutDeltaPhi_Bin3.push_back("Jet2 #eta");
    xTitle_baseline_withoutDeltaPhi_Bin3.push_back("Jet3 #eta");
    xTitle_baseline_withoutDeltaPhi_Bin3.push_back("#Delta#phi 1");
    xTitle_baseline_withoutDeltaPhi_Bin3.push_back("#Delta#phi 2");
    xTitle_baseline_withoutDeltaPhi_Bin3.push_back("#Delta#phi 3");

    vector<TString> xTitle_baseline_withoutDeltaPhi_Bin4;
    xTitle_baseline_withoutDeltaPhi_Bin4.push_back("Jet1 p_{T} (GeV)");
    xTitle_baseline_withoutDeltaPhi_Bin4.push_back("Jet2 p_{T} (GeV)");
    xTitle_baseline_withoutDeltaPhi_Bin4.push_back("Jet3 p_{T} (GeV)");
    xTitle_baseline_withoutDeltaPhi_Bin4.push_back("Jet1 #eta");
    xTitle_baseline_withoutDeltaPhi_Bin4.push_back("Jet2 #eta");
    xTitle_baseline_withoutDeltaPhi_Bin4.push_back("Jet3 #eta");
    xTitle_baseline_withoutDeltaPhi_Bin4.push_back("#Delta#phi 1");
    xTitle_baseline_withoutDeltaPhi_Bin4.push_back("#Delta#phi 2");
    xTitle_baseline_withoutDeltaPhi_Bin4.push_back("#Delta#phi 3");

    vector<TString> xTitle_VBF_presel;
    xTitle_VBF_presel.push_back("#Delta#phi (j_{1}, j_{2})");
    xTitle_VBF_presel.push_back("#Delta#eta (j_{1}, j_{2})");
    xTitle_VBF_presel.push_back("M(j_{1},j_{2}) (GeV)");
    xTitle_VBF_presel.push_back("Jet1 p_{T} (GeV)");
    xTitle_VBF_presel.push_back("Jet2 p_{T} (GeV)");
    xTitle_VBF_presel.push_back("Jet3 p_{T} (GeV)");
    xTitle_VBF_presel.push_back("Jet1 #eta");
    xTitle_VBF_presel.push_back("Jet2 #eta");
    xTitle_VBF_presel.push_back("Jet3 #eta");
    xTitle_VBF_presel.push_back("p_{T}(j_{1},j_{2}) (GeV)");
    xTitle_VBF_presel.push_back("MET (GeV)");
    xTitle_VBF_presel.push_back("METsoft (GeV)");
    xTitle_VBF_presel.push_back("METsig (GeV^{1/2})");
    xTitle_VBF_presel.push_back("MHTsig (GeV^{1/2})");
    xTitle_VBF_presel.push_back("#Delta#phi(MET,j_{1})");
    xTitle_VBF_presel.push_back("#Delta#phi(MET,j_{2})");
    xTitle_VBF_presel.push_back("#Delta#phi(MET,j_{3}))");

    vector<TString> xTitle_VBF_presel_4JV_dPhiSide;
    xTitle_VBF_presel_4JV_dPhiSide.push_back("#Delta#phi (j_{1}, j_{2})");
    xTitle_VBF_presel_4JV_dPhiSide.push_back("#Delta#eta (j_{1}, j_{2})");
    xTitle_VBF_presel_4JV_dPhiSide.push_back("M(j_{1},j_{2}) (GeV)");
    xTitle_VBF_presel_4JV_dPhiSide.push_back("Jet1 p_{T} (GeV)");
    xTitle_VBF_presel_4JV_dPhiSide.push_back("Jet2 p_{T} (GeV)");
    xTitle_VBF_presel_4JV_dPhiSide.push_back("Jet3 p_{T} (GeV)");
    xTitle_VBF_presel_4JV_dPhiSide.push_back("Jet1 #eta");
    xTitle_VBF_presel_4JV_dPhiSide.push_back("Jet2 #eta");
    xTitle_VBF_presel_4JV_dPhiSide.push_back("Jet3 #eta");
    xTitle_VBF_presel_4JV_dPhiSide.push_back("p_{T}(j_{1},j_{2}) (GeV)");
    xTitle_VBF_presel_4JV_dPhiSide.push_back("MET (GeV)");
    xTitle_VBF_presel_4JV_dPhiSide.push_back("METsoft (GeV)");
    xTitle_VBF_presel_4JV_dPhiSide.push_back("METsig (GeV^{1/2})");
    xTitle_VBF_presel_4JV_dPhiSide.push_back("MHTsig (GeV^{1/2})");
    xTitle_VBF_presel_4JV_dPhiSide.push_back("#Delta#phi(MET,j_{1})");
    xTitle_VBF_presel_4JV_dPhiSide.push_back("#Delta#phi(MET,j_{2})");
    xTitle_VBF_presel_4JV_dPhiSide.push_back("#Delta#phi(MET,j_{3})");

    vector<TString> xTitle_VBF_dEta;
    xTitle_VBF_dEta.push_back("#Delta#phi (j_{1}, j_{2})");
    xTitle_VBF_dEta.push_back("#Delta#eta (j_{1}, j_{2})");
    xTitle_VBF_dEta.push_back("M(j_{1},j_{2}) (GeV)");
    xTitle_VBF_dEta.push_back("Jet1 p_{T} (GeV)");
    xTitle_VBF_dEta.push_back("Jet2 p_{T} (GeV)");
    xTitle_VBF_dEta.push_back("Jet3 p_{T} (GeV)");
    xTitle_VBF_dEta.push_back("Jet1 #eta");
    xTitle_VBF_dEta.push_back("Jet2 #eta");
    xTitle_VBF_dEta.push_back("Jet3 #eta");
    xTitle_VBF_dEta.push_back("p_{T}(j_{1},j_{2}) (GeV)");
    xTitle_VBF_dEta.push_back("MET (GeV)");
    xTitle_VBF_dEta.push_back("METsoft (GeV)");
    xTitle_VBF_dEta.push_back("METsig (GeV^{1/2})");
    xTitle_VBF_dEta.push_back("MHTsig (GeV^{1/2})");
    xTitle_VBF_dEta.push_back("#Delta#phi(MET,j_{1})");
    xTitle_VBF_dEta.push_back("#Delta#phi(MET,j_{2})");
    xTitle_VBF_dEta.push_back("#Delta#phi(MET,j_{3})");

    vector<TString> xTitle_VBF_dEta_3JV;
    xTitle_VBF_dEta_3JV.push_back("#Delta#phi (j_{1}, j_{2})");
    xTitle_VBF_dEta_3JV.push_back("#Delta#eta (j_{1}, j_{2})");
    xTitle_VBF_dEta_3JV.push_back("M(j_{1},j_{2}) (GeV)");
    xTitle_VBF_dEta_3JV.push_back("Jet1 p_{T} (GeV)");
    xTitle_VBF_dEta_3JV.push_back("Jet2 p_{T} (GeV)");
    xTitle_VBF_dEta_3JV.push_back("Jet3 p_{T} (GeV)");
    xTitle_VBF_dEta_3JV.push_back("Jet1 #eta");
    xTitle_VBF_dEta_3JV.push_back("Jet2 #eta");
    xTitle_VBF_dEta_3JV.push_back("Jet3 #eta");
    xTitle_VBF_dEta_3JV.push_back("p_{T}(j_{1},j_{2}) (GeV)");
    xTitle_VBF_dEta_3JV.push_back("MET (GeV)");
    xTitle_VBF_dEta_3JV.push_back("METsoft (GeV)");
    xTitle_VBF_dEta_3JV.push_back("METsig (GeV^{1/2})");
    xTitle_VBF_dEta_3JV.push_back("MHTsig (GeV^{1/2})");
    xTitle_VBF_dEta_3JV.push_back("#Delta#phi(MET,j_{1})");
    xTitle_VBF_dEta_3JV.push_back("#Delta#phi(MET,j_{2})");
    xTitle_VBF_dEta_3JV.push_back("#Delta#phi(MET,j_{3})");

    vector<TString> xTitle_VBF_dEta_3JV_dPhiPTjj;
    xTitle_VBF_dEta_3JV_dPhiPTjj.push_back("#Delta#phi (j_{1}, j_{2})");
    xTitle_VBF_dEta_3JV_dPhiPTjj.push_back("#Delta#eta (j_{1}, j_{2})");
    xTitle_VBF_dEta_3JV_dPhiPTjj.push_back("M(j_{1},j_{2}) (GeV)");
    xTitle_VBF_dEta_3JV_dPhiPTjj.push_back("Jet1 p_{T} (GeV)");
    xTitle_VBF_dEta_3JV_dPhiPTjj.push_back("Jet2 p_{T} (GeV)");
    xTitle_VBF_dEta_3JV_dPhiPTjj.push_back("Jet3 p_{T} (GeV)");
    xTitle_VBF_dEta_3JV_dPhiPTjj.push_back("Jet1 #eta");
    xTitle_VBF_dEta_3JV_dPhiPTjj.push_back("Jet2 #eta");
    xTitle_VBF_dEta_3JV_dPhiPTjj.push_back("Jet3 #eta");
    xTitle_VBF_dEta_3JV_dPhiPTjj.push_back("p_{T}(j_{1},j_{2}) (GeV)");
    xTitle_VBF_dEta_3JV_dPhiPTjj.push_back("MET (GeV)");
    xTitle_VBF_dEta_3JV_dPhiPTjj.push_back("METsoft (GeV)");
    xTitle_VBF_dEta_3JV_dPhiPTjj.push_back("METsig (GeV^{1/2})");
    xTitle_VBF_dEta_3JV_dPhiPTjj.push_back("MHTsig (GeV^{1/2})");
    xTitle_VBF_dEta_3JV_dPhiPTjj.push_back("#Delta#phi(MET,j_{1})");
    xTitle_VBF_dEta_3JV_dPhiPTjj.push_back("#Delta#phi(MET,j_{2})");
    xTitle_VBF_dEta_3JV_dPhiPTjj.push_back("#Delta#phi(MET,j_{3})");

    vector<TString> xTitle_VBF_jj;
    xTitle_VBF_jj.push_back("#Delta#phi (j_{1}, j_{2})");
    xTitle_VBF_jj.push_back("#Delta#eta (j_{1}, j_{2})");
    xTitle_VBF_jj.push_back("M(j_{1},j_{2}) (GeV)");
    xTitle_VBF_jj.push_back("Jet1 p_{T} (GeV)");
    xTitle_VBF_jj.push_back("Jet2 p_{T} (GeV)");
    xTitle_VBF_jj.push_back("Jet3 p_{T} (GeV)");
    xTitle_VBF_jj.push_back("Jet1 #eta");
    xTitle_VBF_jj.push_back("Jet2 #eta");
    xTitle_VBF_jj.push_back("Jet3 #eta");
    xTitle_VBF_jj.push_back("p_{T}(j_{1},j_{2}) (GeV)");
    xTitle_VBF_jj.push_back("MET (GeV)");
    xTitle_VBF_jj.push_back("METsoft (GeV)");
    xTitle_VBF_jj.push_back("METsig (GeV^{1/2})");
    xTitle_VBF_jj.push_back("MHTsig (GeV^{1/2})");
    xTitle_VBF_jj.push_back("#Delta#phi(MET,j_{1})");
    xTitle_VBF_jj.push_back("#Delta#phi(MET,j_{2})");
    xTitle_VBF_jj.push_back("#Delta#phi(MET,j_{3})");

    vector<TString> hist_type_presel;
    hist_type_presel.push_back("HT_presel");
    hist_type_presel.push_back("MHT_presel");
    hist_type_presel.push_back("MET_presel");
    hist_type_presel.push_back("NJets_presel");
    hist_type_presel.push_back("NBJets_presel");
    hist_type_presel.push_back("Jet1Pt_presel");
    hist_type_presel.push_back("Jet2Pt_presel");
    hist_type_presel.push_back("Jet1Eta_presel");
    hist_type_presel.push_back("Jet2Eta_presel");

    vector<TString> hist_type_deltaPhi;
    hist_type_deltaPhi.push_back("HT_deltaPhi");
    hist_type_deltaPhi.push_back("MHT_deltaPhi");
    hist_type_deltaPhi.push_back("MET_deltaPhi");
    hist_type_deltaPhi.push_back("Jet1Pt_deltaPhi");
    hist_type_deltaPhi.push_back("Jet2Pt_deltaPhi");
    hist_type_deltaPhi.push_back("Jet1Eta_deltaPhi");
    hist_type_deltaPhi.push_back("Jet2Eta_deltaPhi");

    vector<TString> hist_type_baseline_Bin1;
    hist_type_baseline_Bin1.push_back("Jet1Pt_JetBin1_baseline");
    hist_type_baseline_Bin1.push_back("Jet2Pt_JetBin1_baseline");
    hist_type_baseline_Bin1.push_back("Jet1Eta_JetBin1_baseline");
    hist_type_baseline_Bin1.push_back("Jet2Eta_JetBin1_baseline");
    hist_type_baseline_Bin1.push_back("DeltaPhi1_JetBin1_baseline");
    hist_type_baseline_Bin1.push_back("DeltaPhi2_JetBin1_baseline");

    vector<TString> hist_type_baseline_Bin2;
    hist_type_baseline_Bin2.push_back("Jet1Pt_JetBin2_baseline");
    hist_type_baseline_Bin2.push_back("Jet2Pt_JetBin2_baseline");
    hist_type_baseline_Bin2.push_back("Jet3Pt_JetBin2_baseline");
    hist_type_baseline_Bin2.push_back("Jet1Eta_JetBin2_baseline");
    hist_type_baseline_Bin2.push_back("Jet2Eta_JetBin2_baseline");
    hist_type_baseline_Bin2.push_back("Jet3Eta_JetBin2_baseline");
    hist_type_baseline_Bin2.push_back("DeltaPhi1_JetBin2_baseline");
    hist_type_baseline_Bin2.push_back("DeltaPhi2_JetBin2_baseline");
    hist_type_baseline_Bin2.push_back("DeltaPhi3_JetBin2_baseline");

    vector<TString> hist_type_baseline_Bin3;
    hist_type_baseline_Bin3.push_back("Jet1Pt_JetBin3_baseline");
    hist_type_baseline_Bin3.push_back("Jet2Pt_JetBin3_baseline");
    hist_type_baseline_Bin3.push_back("Jet3Pt_JetBin3_baseline");
    hist_type_baseline_Bin3.push_back("Jet1Eta_JetBin3_baseline");
    hist_type_baseline_Bin3.push_back("Jet2Eta_JetBin3_baseline");
    hist_type_baseline_Bin3.push_back("Jet3Eta_JetBin3_baseline");
    hist_type_baseline_Bin3.push_back("DeltaPhi1_JetBin3_baseline");
    hist_type_baseline_Bin3.push_back("DeltaPhi2_JetBin3_baseline");
    hist_type_baseline_Bin3.push_back("DeltaPhi3_JetBin3_baseline");

    vector<TString> hist_type_baseline_Bin4;
    hist_type_baseline_Bin4.push_back("Jet1Pt_JetBin4_baseline");
    hist_type_baseline_Bin4.push_back("Jet2Pt_JetBin4_baseline");
    hist_type_baseline_Bin4.push_back("Jet3Pt_JetBin4_baseline");
    hist_type_baseline_Bin4.push_back("Jet1Eta_JetBin4_baseline");
    hist_type_baseline_Bin4.push_back("Jet2Eta_JetBin4_baseline");
    hist_type_baseline_Bin4.push_back("Jet3Eta_JetBin4_baseline");
    hist_type_baseline_Bin4.push_back("DeltaPhi1_JetBin4_baseline");
    hist_type_baseline_Bin4.push_back("DeltaPhi2_JetBin4_baseline");
    hist_type_baseline_Bin4.push_back("DeltaPhi3_JetBin4_baseline");

    vector<TString> hist_type_baseline_withoutDeltaPhi_Bin1;
    hist_type_baseline_withoutDeltaPhi_Bin1.push_back("Jet1Pt_JetBin1_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin1.push_back("Jet2Pt_JetBin1_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin1.push_back("Jet1Eta_JetBin1_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin1.push_back("Jet2Eta_JetBin1_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin1.push_back("DeltaPhi1_JetBin1_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin1.push_back("DeltaPhi2_JetBin1_baseline_withoutDeltaPhi");

    vector<TString> hist_type_baseline_withoutDeltaPhi_Bin2;
    hist_type_baseline_withoutDeltaPhi_Bin2.push_back("Jet1Pt_JetBin2_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin2.push_back("Jet2Pt_JetBin2_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin2.push_back("Jet3Pt_JetBin2_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin2.push_back("Jet1Eta_JetBin2_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin2.push_back("Jet2Eta_JetBin2_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin2.push_back("Jet3Eta_JetBin2_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin2.push_back("DeltaPhi1_JetBin2_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin2.push_back("DeltaPhi2_JetBin2_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin2.push_back("DeltaPhi3_JetBin2_baseline_withoutDeltaPhi");

    vector<TString> hist_type_baseline_withoutDeltaPhi_Bin3;
    hist_type_baseline_withoutDeltaPhi_Bin3.push_back("Jet1Pt_JetBin3_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin3.push_back("Jet2Pt_JetBin3_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin3.push_back("Jet3Pt_JetBin3_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin3.push_back("Jet1Eta_JetBin3_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin3.push_back("Jet2Eta_JetBin3_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin3.push_back("Jet3Eta_JetBin3_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin3.push_back("DeltaPhi1_JetBin3_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin3.push_back("DeltaPhi2_JetBin3_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin3.push_back("DeltaPhi3_JetBin3_baseline_withoutDeltaPhi");

    vector<TString> hist_type_baseline_withoutDeltaPhi_Bin4;
    hist_type_baseline_withoutDeltaPhi_Bin4.push_back("Jet1Pt_JetBin4_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin4.push_back("Jet2Pt_JetBin4_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin4.push_back("Jet3Pt_JetBin4_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin4.push_back("Jet1Eta_JetBin4_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin4.push_back("Jet2Eta_JetBin4_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin4.push_back("Jet3Eta_JetBin4_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin4.push_back("DeltaPhi1_JetBin4_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin4.push_back("DeltaPhi2_JetBin4_baseline_withoutDeltaPhi");
    hist_type_baseline_withoutDeltaPhi_Bin4.push_back("DeltaPhi3_JetBin4_baseline_withoutDeltaPhi");

    vector<TString> hist_type_VBF_presel;
    hist_type_VBF_presel.push_back("VBF_dPhi_presel");
    hist_type_VBF_presel.push_back("VBF_dEta_presel");
    hist_type_VBF_presel.push_back("VBF_Mjj_presel");
    hist_type_VBF_presel.push_back("VBF_Jet1Pt_presel");
    hist_type_VBF_presel.push_back("VBF_Jet2Pt_presel");
    hist_type_VBF_presel.push_back("VBF_Jet3Pt_presel");
    hist_type_VBF_presel.push_back("VBF_Jet1Eta_presel");
    hist_type_VBF_presel.push_back("VBF_Jet2Eta_presel");
    hist_type_VBF_presel.push_back("VBF_Jet3Eta_presel");
    hist_type_VBF_presel.push_back("VBF_PTjj_presel");
    hist_type_VBF_presel.push_back("VBF_MET_presel");
    hist_type_VBF_presel.push_back("VBF_METsoft_presel");
    hist_type_VBF_presel.push_back("VBF_METsig_presel");
    hist_type_VBF_presel.push_back("VBF_MHTsig_presel");
    hist_type_VBF_presel.push_back("VBF_minDeltaPhiPTj12_presel");
    hist_type_VBF_presel.push_back("VBF_maxDeltaPhiPTj12_presel");
    hist_type_VBF_presel.push_back("VBF_DeltaPhiPTj3_presel");

    vector<TString> hist_type_VBF_presel_4JV_dPhiSide;
    hist_type_VBF_presel_4JV_dPhiSide.push_back("VBF_dPhi_presel_4JV_dPhiSide");
    hist_type_VBF_presel_4JV_dPhiSide.push_back("VBF_dEta_presel_4JV_dPhiSide");
    hist_type_VBF_presel_4JV_dPhiSide.push_back("VBF_Mjj_presel_4JV_dPhiSide");
    hist_type_VBF_presel_4JV_dPhiSide.push_back("VBF_Jet1Pt_presel_4JV_dPhiSide");
    hist_type_VBF_presel_4JV_dPhiSide.push_back("VBF_Jet2Pt_presel_4JV_dPhiSide");
    hist_type_VBF_presel_4JV_dPhiSide.push_back("VBF_Jet3Pt_presel_4JV_dPhiSide");
    hist_type_VBF_presel_4JV_dPhiSide.push_back("VBF_Jet1Eta_presel_4JV_dPhiSide");
    hist_type_VBF_presel_4JV_dPhiSide.push_back("VBF_Jet2Eta_presel_4JV_dPhiSide");
    hist_type_VBF_presel_4JV_dPhiSide.push_back("VBF_Jet3Eta_presel_4JV_dPhiSide");
    hist_type_VBF_presel_4JV_dPhiSide.push_back("VBF_PTjj_presel_4JV_dPhiSide");
    hist_type_VBF_presel_4JV_dPhiSide.push_back("VBF_MET_presel_4JV_dPhiSide");
    hist_type_VBF_presel_4JV_dPhiSide.push_back("VBF_METsoft_presel_4JV_dPhiSide");
    hist_type_VBF_presel_4JV_dPhiSide.push_back("VBF_METsig_presel_4JV_dPhiSide");
    hist_type_VBF_presel_4JV_dPhiSide.push_back("VBF_MHTsig_presel_4JV_dPhiSide");
    hist_type_VBF_presel_4JV_dPhiSide.push_back("VBF_minDeltaPhiPTj12_presel_4JV_dPhiSide");
    hist_type_VBF_presel_4JV_dPhiSide.push_back("VBF_maxDeltaPhiPTj12_presel_4JV_dPhiSide");
    hist_type_VBF_presel_4JV_dPhiSide.push_back("VBF_DeltaPhiPTj3_presel_4JV_dPhiSide");

    vector<TString> hist_type_VBF_dEta;
    hist_type_VBF_dEta.push_back("VBF_dPhi_dEta");
    hist_type_VBF_dEta.push_back("VBF_dEta_dEta");
    hist_type_VBF_dEta.push_back("VBF_Mjj_dEta");
    hist_type_VBF_dEta.push_back("VBF_Jet1Pt_dEta");
    hist_type_VBF_dEta.push_back("VBF_Jet2Pt_dEta");
    hist_type_VBF_dEta.push_back("VBF_Jet3Pt_dEta");
    hist_type_VBF_dEta.push_back("VBF_Jet1Eta_dEta");
    hist_type_VBF_dEta.push_back("VBF_Jet2Eta_dEta");
    hist_type_VBF_dEta.push_back("VBF_Jet3Eta_dEta");
    hist_type_VBF_dEta.push_back("VBF_PTjj_dEta");
    hist_type_VBF_dEta.push_back("VBF_MET_dEta");
    hist_type_VBF_dEta.push_back("VBF_METsoft_dEta");
    hist_type_VBF_dEta.push_back("VBF_METsig_dEta");
    hist_type_VBF_dEta.push_back("VBF_MHTsig_dEta");
    hist_type_VBF_dEta.push_back("VBF_minDeltaPhiPTj12_dEta");
    hist_type_VBF_dEta.push_back("VBF_maxDeltaPhiPTj12_dEta");
    hist_type_VBF_dEta.push_back("VBF_DeltaPhiPTj3_dEta");

    vector<TString> hist_type_VBF_dEta_3JV;
    hist_type_VBF_dEta_3JV.push_back("VBF_dPhi_dEta_3JV");
    hist_type_VBF_dEta_3JV.push_back("VBF_dEta_dEta_3JV");
    hist_type_VBF_dEta_3JV.push_back("VBF_Mjj_dEta_3JV");
    hist_type_VBF_dEta_3JV.push_back("VBF_Jet1Pt_dEta_3JV");
    hist_type_VBF_dEta_3JV.push_back("VBF_Jet2Pt_dEta_3JV");
    hist_type_VBF_dEta_3JV.push_back("VBF_Jet3Pt_dEta_3JV");
    hist_type_VBF_dEta_3JV.push_back("VBF_Jet1Eta_dEta_3JV");
    hist_type_VBF_dEta_3JV.push_back("VBF_Jet2Eta_dEta_3JV");
    hist_type_VBF_dEta_3JV.push_back("VBF_Jet3Eta_dEta_3JV");
    hist_type_VBF_dEta_3JV.push_back("VBF_PTjj_dEta_3JV");
    hist_type_VBF_dEta_3JV.push_back("VBF_MET_dEta_3JV");
    hist_type_VBF_dEta_3JV.push_back("VBF_METsoft_dEta_3JV");
    hist_type_VBF_dEta_3JV.push_back("VBF_METsig_dEta_3JV");
    hist_type_VBF_dEta_3JV.push_back("VBF_MHTsig_dEta_3JV");
    hist_type_VBF_dEta_3JV.push_back("VBF_minDeltaPhiPTj12_dEta_3JV");
    hist_type_VBF_dEta_3JV.push_back("VBF_maxDeltaPhiPTj12_dEta_3JV");
    hist_type_VBF_dEta_3JV.push_back("VBF_DeltaPhiPTj3_dEta_3JV");

    vector<TString> hist_type_VBF_dEta_3JV_dPhiPTjj;
    hist_type_VBF_dEta_3JV_dPhiPTjj.push_back("VBF_dPhi_dEta_3JV_dPhiPTjj");
    hist_type_VBF_dEta_3JV_dPhiPTjj.push_back("VBF_dEta_dEta_3JV_dPhiPTjj");
    hist_type_VBF_dEta_3JV_dPhiPTjj.push_back("VBF_Mjj_dEta_3JV_dPhiPTjj");
    hist_type_VBF_dEta_3JV_dPhiPTjj.push_back("VBF_Jet1Pt_dEta_3JV_dPhiPTjj");
    hist_type_VBF_dEta_3JV_dPhiPTjj.push_back("VBF_Jet2Pt_dEta_3JV_dPhiPTjj");
    hist_type_VBF_dEta_3JV_dPhiPTjj.push_back("VBF_Jet3Pt_dEta_3JV_dPhiPTjj");
    hist_type_VBF_dEta_3JV_dPhiPTjj.push_back("VBF_Jet1Eta_dEta_3JV_dPhiPTjj");
    hist_type_VBF_dEta_3JV_dPhiPTjj.push_back("VBF_Jet2Eta_dEta_3JV_dPhiPTjj");
    hist_type_VBF_dEta_3JV_dPhiPTjj.push_back("VBF_Jet3Eta_dEta_3JV_dPhiPTjj");
    hist_type_VBF_dEta_3JV_dPhiPTjj.push_back("VBF_PTjj_dEta_3JV_dPhiPTjj");
    hist_type_VBF_dEta_3JV_dPhiPTjj.push_back("VBF_MET_dEta_3JV_dPhiPTjj");
    hist_type_VBF_dEta_3JV_dPhiPTjj.push_back("VBF_METsoft_dEta_3JV_dPhiPTjj");
    hist_type_VBF_dEta_3JV_dPhiPTjj.push_back("VBF_METsig_dEta_3JV_dPhiPTjj");
    hist_type_VBF_dEta_3JV_dPhiPTjj.push_back("VBF_MHTsig_dEta_3JV_dPhiPTjj");
    hist_type_VBF_dEta_3JV_dPhiPTjj.push_back("VBF_minDeltaPhiPTj12_dEta_3JV_dPhiPTjj");
    hist_type_VBF_dEta_3JV_dPhiPTjj.push_back("VBF_maxDeltaPhiPTj12_dEta_3JV_dPhiPTjj");
    hist_type_VBF_dEta_3JV_dPhiPTjj.push_back("VBF_DeltaPhiPTj3_dEta_3JV_dPhiPTjj");

    vector<TString> hist_type_VBF_jj;
    hist_type_VBF_jj.push_back("VBF_dPhi_jj");
    hist_type_VBF_jj.push_back("VBF_dEta_jj");
    hist_type_VBF_jj.push_back("VBF_Mjj_jj");
    hist_type_VBF_jj.push_back("VBF_Jet1Pt_jj");
    hist_type_VBF_jj.push_back("VBF_Jet2Pt_jj");
    hist_type_VBF_jj.push_back("VBF_Jet3Pt_jj");
    hist_type_VBF_jj.push_back("VBF_Jet1Eta_jj");
    hist_type_VBF_jj.push_back("VBF_Jet2Eta_jj");
    hist_type_VBF_jj.push_back("VBF_Jet3Eta_jj");
    hist_type_VBF_jj.push_back("VBF_PTjj_jj");
    hist_type_VBF_jj.push_back("VBF_MET_jj");
    hist_type_VBF_jj.push_back("VBF_METsoft_jj");
    hist_type_VBF_jj.push_back("VBF_METsig_jj");
    hist_type_VBF_jj.push_back("VBF_MHTsig_jj");
    hist_type_VBF_jj.push_back("VBF_minDeltaPhiPTj12_jj");
    hist_type_VBF_jj.push_back("VBF_maxDeltaPhiPTj12_jj");
    hist_type_VBF_jj.push_back("VBF_DeltaPhiPTj3_jj");

    // --------------------------------------------------------------------------------------------- //
    // plots for Mjj
    TString Title;
    TString yTitle = "Events";

    if (VBF) {

        if (VBFSR) {
            Title = "N_{j}>=3, M_{jj}>1.0 TeV, MET>150 GeV, #Delta#phi<2.7, #Delta#eta>2.5, p_{T}(j3)<50 GeV";
        } else {
            Title = "N_{j}>=3, M_{jj}>0.6 TeV, MET>100 GeV, #Delta#phi<2.7, #Delta#eta>2.5, p_{T}(j3)<50 GeV";
        }

        if( hist_type_VBF_presel.size() != xTitle_VBF_presel.size() ) cout << "Error: Missing xTitles VBF_presel!!" << endl;

        for(int i = 0; i < hist_type_VBF_presel.size(); i++ ) {
            TCanvas *c = DrawComparison( pred_->GetPredictionHisto(hist_type_VBF_presel.at(i)), pred_->GetSelectionHisto(hist_type_VBF_presel.at(i)), Title, LumiTitle, xTitle_VBF_presel.at(i), yTitle, isData);
            c->Print("output_GetPrediction/" + hist_type_VBF_presel.at(i) + postfix + ".pdf");
        }

        if (VBFSR) {
            Title = "N_{j}=3, M_{jj}>1.0 TeV, MET>150 GeV, 1.8<#Delta#phi<2.7, #Delta#eta>2.5, p_{T}(j3)<50 GeV";
        } else {
            Title = "N_{j}=3, M_{jj}>0.6 TeV, MET>100 GeV, 1.8<#Delta#phi<2.7, #Delta#eta>2.5, p_{T}(j3)<50 GeV";
        }

        if( hist_type_VBF_presel_4JV_dPhiSide.size() != xTitle_VBF_presel_4JV_dPhiSide.size() ) cout << "Error: Missing xTitles VBF_presel!!" << endl;

        for(int i = 0; i < hist_type_VBF_presel_4JV_dPhiSide.size(); i++ ) {
            TCanvas *c = DrawComparison( pred_->GetPredictionHisto(hist_type_VBF_presel_4JV_dPhiSide.at(i)), pred_->GetSelectionHisto(hist_type_VBF_presel_4JV_dPhiSide.at(i)), Title, LumiTitle, xTitle_VBF_presel_4JV_dPhiSide.at(i), yTitle, isData);
            c->Print("output_GetPrediction/" + hist_type_VBF_presel_4JV_dPhiSide.at(i) + postfix + ".pdf");
        }

        if (VBFSR) {
            Title = "N_{j}>=3, M_{jj}>1.0 TeV, MET>150 GeV, #Delta#phi<1.8, #Delta#eta>4.8, p_{T}(j3)<50 GeV";
        } else {
            Title = "N_{j}>=3, M_{jj}>0.6 TeV, MET>100 GeV, #Delta#phi<1.8, #Delta#eta>4.8, p_{T}(j3)<50 GeV";
        }

        if( hist_type_VBF_dEta.size() != xTitle_VBF_dEta.size() ) cout << "Error: Missing xTitles VBF_dEta!!" << endl;

        for(int i = 0; i < hist_type_VBF_dEta.size(); i++ ) {
            TCanvas *c = DrawComparison( pred_->GetPredictionHisto(hist_type_VBF_dEta.at(i)), pred_->GetSelectionHisto(hist_type_VBF_dEta.at(i)), Title, LumiTitle, xTitle_VBF_dEta.at(i), yTitle, isData);
            c->Print("output_GetPrediction/" + hist_type_VBF_dEta.at(i) + postfix + ".pdf");
        }

        if (VBFSR) {
            Title = "N_{j}=2, M_{jj}>1.0 TeV, MET>150 GeV, #Delta#phi<1.8, #Delta#eta>4.8, p_{T}(j3)<25 GeV";
        } else {
            Title = "N_{j}=2, M_{jj}>0.6 TeV, MET>100 GeV, #Delta#phi<1.8, #Delta#eta>4.8, p_{T}(j3)<25 GeV";
        }

        if( hist_type_VBF_dEta_3JV.size() != xTitle_VBF_dEta_3JV.size() ) cout << "Error: Missing xTitles VBF_dEta_3JV!!" << endl;

        for(int i = 0; i < hist_type_VBF_dEta_3JV.size(); i++ ) {
            TCanvas *c = DrawComparison( pred_->GetPredictionHisto(hist_type_VBF_dEta_3JV.at(i)), pred_->GetSelectionHisto(hist_type_VBF_dEta_3JV.at(i)), Title, LumiTitle, xTitle_VBF_dEta_3JV.at(i), yTitle, isData);
            c->Print("output_GetPrediction/" + hist_type_VBF_dEta_3JV.at(i) + postfix + ".pdf");
        }

        if (VBFSR) {
            Title = "N_{j}=2, M_{jj}>1.0 TeV, MET>150 GeV, #Delta#phi<1.8, #Delta#eta>2.5, p_{T}(j3)<25 GeV";
        } else {
            Title = "N_{j}=2, M_{jj}>0.6 TeV, MET>100 GeV, #Delta#phi<1.8, #Delta#eta>2.5, p_{T}(j3)<25 GeV";
        }

        if( hist_type_VBF_jj.size() != xTitle_VBF_jj.size() ) cout << "Error: Missing xTitles VBF_jj!!" << endl;

        for(int i = 0; i < hist_type_VBF_jj.size(); i++ ) {
            TCanvas *c = DrawComparison( pred_->GetPredictionHisto(hist_type_VBF_jj.at(i)), pred_->GetSelectionHisto(hist_type_VBF_jj.at(i)), Title, LumiTitle, xTitle_VBF_jj.at(i), yTitle, isData);
            c->Print("output_GetPrediction/" + hist_type_VBF_jj.at(i) + postfix + ".pdf");
        }

        if (VBFSR) {
            Title = "N_{j}=2, M_{jj}>1.0 TeV, MET>150 GeV, #Delta#phi<1.8, #Delta#eta>4.8, #Delta#phi(MET,j)>1.0, p_{T}(j3)<25 GeV";
        } else {
            Title = "N_{j}=2, M_{jj}>0.6 TeV, MET>100 GeV, #Delta#phi<1.8, #Delta#eta>4.8, #Delta#phi(MET,j)>1.0, p_{T}(j3)<25 GeV";
        }

        if( hist_type_VBF_dEta_3JV_dPhiPTjj.size() != xTitle_VBF_dEta_3JV_dPhiPTjj.size() ) cout << "Error: Missing xTitles VBF_dEta_3JV_dPhiPTjj!!" << endl;

        for(int i = 0; i < hist_type_VBF_dEta_3JV_dPhiPTjj.size(); i++ ) {
            TCanvas *c = DrawComparison( pred_->GetPredictionHisto(hist_type_VBF_dEta_3JV_dPhiPTjj.at(i)), pred_->GetSelectionHisto(hist_type_VBF_dEta_3JV_dPhiPTjj.at(i)), Title, LumiTitle, xTitle_VBF_dEta_3JV_dPhiPTjj.at(i), yTitle, isData);
            c->Print("output_GetPrediction/" + hist_type_VBF_dEta_3JV_dPhiPTjj.at(i) + postfix + ".pdf");
        }
    }

    if (HTMHT) {
        // --------------------------------------------------------------------------------------------- //
        // plots for preselection (2 jets)
        yTitle = "Events";
        Title = ">=2 jets";

        if( hist_type_presel.size() != xTitle_presel.size() ) cout << "Error: Missing xTitles preselection!!" << endl;

        for(int i = 0; i < hist_type_presel.size(); i++ ) {
            TCanvas *c = DrawComparison( pred_->GetPredictionHisto(hist_type_presel.at(i)), pred_->GetSelectionHisto(hist_type_presel.at(i)), Title, LumiTitle, xTitle_presel.at(i), yTitle, isData);
            c->Print("output_GetPrediction/" + hist_type_presel.at(i) + postfix + ".pdf");
        }

        // --------------------------------------------------------------------------------------------- //
        // deltaPhi after preselection
        Title = ">=2 jets";
        TCanvas *c =  DrawComparison( pred_->GetPredictionHisto("DeltaPhi1_presel"), pred_->GetSelectionHisto("DeltaPhi1_presel"), Title, LumiTitle,"#Delta#phi (jet1, MET)", yTitle, isData);
        c->Print("output_GetPrediction/DeltaPhi1_presel" + postfix + ".pdf");

        Title = ">=2 jets";
        c =  DrawComparison( pred_->GetPredictionHisto("DeltaPhi2_presel"), pred_->GetSelectionHisto("DeltaPhi2_presel"), Title, LumiTitle,"#Delta#phi (jet2, MET)", yTitle, isData);
        c->Print("output_GetPrediction/DeltaPhi2_presel" + postfix + ".pdf");

        Title = ">=2 jets";
        c =  DrawComparison( pred_->GetPredictionHisto("DeltaPhi3_presel"), pred_->GetSelectionHisto("DeltaPhi3_presel"), Title, LumiTitle,"#Delta#phi (jet3, MET)", yTitle, isData);
        c->Print("output_GetPrediction/DeltaPhi3_presel" + postfix + ".pdf");

        // --------------------------------------------------------------------------------------------- //
        // plots for preselection + deltaPhi cut
        Title = ">=4 jets, #Delta#phi(#slash{H}_{T}, jet1,2,3)";

        if( hist_type_deltaPhi.size() != xTitle_deltaPhi.size() ) cout << "Error: Missing xTitles preselection!!" << endl;

        for(int i = 0; i < hist_type_deltaPhi.size(); i++ ) {
            TCanvas *c = DrawComparison( pred_->GetPredictionHisto(hist_type_deltaPhi.at(i)), pred_->GetSelectionHisto(hist_type_deltaPhi.at(i)), Title, LumiTitle, xTitle_deltaPhi.at(i), yTitle, isData);
            c->Print("output_GetPrediction/" + hist_type_deltaPhi.at(i) + postfix + ".pdf");
        }

        // --------------------------------------------------------------------------------------------- //
        // plots for baseline
        Title = "2+3 jets, #Delta#phi, HT > 500 GeV, MET > 200 GeV";

        if( hist_type_baseline_Bin1.size() != xTitle_baseline_Bin1.size() ) cout << "Error: Missing xTitles baseline Bin1!!" << endl;

        for(int i = 0; i < hist_type_baseline_Bin1.size(); i++ ) {
            TCanvas *c = DrawComparison( pred_->GetPredictionHisto(hist_type_baseline_Bin1.at(i)), pred_->GetSelectionHisto(hist_type_baseline_Bin1.at(i)), Title, LumiTitle, xTitle_baseline_Bin1.at(i), yTitle, isData);
            c->Print("output_GetPrediction/" + hist_type_baseline_Bin1.at(i) + postfix + ".pdf");
        }

        ////////////////////////////////////////////////////////////////
        Title = "4-6 jets, #Delta#phi, HT > 500 GeV, MET > 200 GeV";

        if( hist_type_baseline_Bin2.size() != xTitle_baseline_Bin2.size() ) cout << "Error: Missing xTitles baseline Bin2!!" << endl;

        for(int i = 0; i < hist_type_baseline_Bin2.size(); i++ ) {
            TCanvas *c = DrawComparison( pred_->GetPredictionHisto(hist_type_baseline_Bin2.at(i)), pred_->GetSelectionHisto(hist_type_baseline_Bin2.at(i)), Title, LumiTitle, xTitle_baseline_Bin2.at(i), yTitle, isData);
            c->Print("output_GetPrediction/" + hist_type_baseline_Bin2.at(i) + postfix + ".pdf");
        }

        ////////////////////////////////////////////////////////////////
        Title = "7+8 jets, #Delta#phi, HT > 500 GeV, MET > 200 GeV";

        if( hist_type_baseline_Bin3.size() != xTitle_baseline_Bin3.size() ) cout << "Error: Missing xTitles baseline Bin3!!" << endl;

        for(int i = 0; i < hist_type_baseline_Bin3.size(); i++ ) {
            TCanvas *c = DrawComparison( pred_->GetPredictionHisto(hist_type_baseline_Bin3.at(i)), pred_->GetSelectionHisto(hist_type_baseline_Bin3.at(i)), Title, LumiTitle, xTitle_baseline_Bin3.at(i), yTitle, isData);
            c->Print("output_GetPrediction/" + hist_type_baseline_Bin3.at(i) + postfix + ".pdf");
        }

        ////////////////////////////////////////////////////////////////
        Title = ">=9 jets, #Delta#phi, HT > 500 GeV, MET > 200 GeV";

        if( hist_type_baseline_Bin4.size() != xTitle_baseline_Bin4.size() ) cout << "Error: Missing xTitles baseline Bin4!!" << endl;

        for(int i = 0; i < hist_type_baseline_Bin4.size(); i++ ) {
            TCanvas *c = DrawComparison( pred_->GetPredictionHisto(hist_type_baseline_Bin4.at(i)), pred_->GetSelectionHisto(hist_type_baseline_Bin4.at(i)), Title, LumiTitle, xTitle_baseline_Bin4.at(i), yTitle, isData);
            c->Print("output_GetPrediction/" + hist_type_baseline_Bin4.at(i) + postfix + ".pdf");
        }

        // --------------------------------------------------------------------------------------------- //
        // plots for baseline without deltaPhi
        Title = "2+3 jets, HT > 500 GeV, MET > 200 GeV";

        if( hist_type_baseline_withoutDeltaPhi_Bin1.size() != xTitle_baseline_withoutDeltaPhi_Bin1.size() ) cout << "Error: Missing xTitles baseline_withoutDeltaPhi Bin1!!" << endl;

        for(int i = 0; i < hist_type_baseline_withoutDeltaPhi_Bin1.size(); i++ ) {
            TCanvas *c = DrawComparison( pred_->GetPredictionHisto(hist_type_baseline_withoutDeltaPhi_Bin1.at(i)), pred_->GetSelectionHisto(hist_type_baseline_withoutDeltaPhi_Bin1.at(i)), Title, LumiTitle, xTitle_baseline_withoutDeltaPhi_Bin1.at(i), yTitle, isData);
            c->Print("output_GetPrediction/" + hist_type_baseline_withoutDeltaPhi_Bin1.at(i) + postfix + ".pdf");
        }

        ////////////////////////////////////////////////////////////////
        Title = "4-6 jets, HT > 500 GeV, MET > 200 GeV";

        if( hist_type_baseline_withoutDeltaPhi_Bin2.size() != xTitle_baseline_withoutDeltaPhi_Bin2.size() ) cout << "Error: Missing xTitles baseline_withoutDeltaPhi Bin2!!" << endl;

        for(int i = 0; i < hist_type_baseline_withoutDeltaPhi_Bin2.size(); i++ ) {
            TCanvas *c = DrawComparison( pred_->GetPredictionHisto(hist_type_baseline_withoutDeltaPhi_Bin2.at(i)), pred_->GetSelectionHisto(hist_type_baseline_withoutDeltaPhi_Bin2.at(i)), Title, LumiTitle, xTitle_baseline_withoutDeltaPhi_Bin2.at(i), yTitle, isData);
            c->Print("output_GetPrediction/" + hist_type_baseline_withoutDeltaPhi_Bin2.at(i) + postfix + ".pdf");
        }

        ////////////////////////////////////////////////////////////////
        Title = "7+8 jets, HT > 500 GeV, MET > 200 GeV";

        if( hist_type_baseline_withoutDeltaPhi_Bin3.size() != xTitle_baseline_withoutDeltaPhi_Bin3.size() ) cout << "Error: Missing xTitles baseline_withoutDeltaPhi Bin3!!" << endl;

        for(int i = 0; i < hist_type_baseline_withoutDeltaPhi_Bin3.size(); i++ ) {
            TCanvas *c = DrawComparison( pred_->GetPredictionHisto(hist_type_baseline_withoutDeltaPhi_Bin3.at(i)), pred_->GetSelectionHisto(hist_type_baseline_withoutDeltaPhi_Bin3.at(i)), Title, LumiTitle, xTitle_baseline_withoutDeltaPhi_Bin3.at(i), yTitle, isData);
            c->Print("output_GetPrediction/" + hist_type_baseline_withoutDeltaPhi_Bin3.at(i) + postfix + ".pdf");
        }

        ////////////////////////////////////////////////////////////////
        Title = ">=9 jets, HT > 500 GeV, MET > 200 GeV";

        if( hist_type_baseline_withoutDeltaPhi_Bin4.size() != xTitle_baseline_withoutDeltaPhi_Bin4.size() ) cout << "Error: Missing xTitles baseline_withoutDeltaPhi Bin4!!" << endl;

        for(int i = 0; i < hist_type_baseline_withoutDeltaPhi_Bin4.size(); i++ ) {
            TCanvas *c = DrawComparison( pred_->GetPredictionHisto(hist_type_baseline_withoutDeltaPhi_Bin4.at(i)), pred_->GetSelectionHisto(hist_type_baseline_withoutDeltaPhi_Bin4.at(i)), Title, LumiTitle, xTitle_baseline_withoutDeltaPhi_Bin4.at(i), yTitle, isData);
            c->Print("output_GetPrediction/" + hist_type_baseline_withoutDeltaPhi_Bin4.at(i) + postfix + ".pdf");
        }

        // --------------------------------------------------------------------------------------------- //
        // plots for HT inclusive
        //jet Bin 1
        Title = "2+3 jets, #Delta#phi cut, HT > 500 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("MHT_JetBin1_HTinclusive"), pred_->GetSelectionHisto("MHT_JetBin1_HTinclusive"), Title, LumiTitle,"#slash{H}_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/MHT_JetBin1_HTinclusive" + postfix + ".pdf");

        // jet Bin2
        Title = "4-6 jets, #Delta#phi cut, HT > 500 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("MHT_JetBin2_HTinclusive"), pred_->GetSelectionHisto("MHT_JetBin2_HTinclusive"), Title, LumiTitle,"#slash{H}_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/MHT_JetBin2_HTinclusive" + postfix + ".pdf");

        // jet Bin3
        Title = "7+8 jets, #Delta#phi cut, HT > 500 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("MHT_JetBin3_HTinclusive"), pred_->GetSelectionHisto("MHT_JetBin3_HTinclusive"), Title, LumiTitle,"#slash{H}_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/MHT_JetBin3_HTinclusive" + postfix + ".pdf");

        // jet Bin 4
        Title = ">=9 jets, #Delta#phi cut, HT > 500 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("MHT_JetBin4_HTinclusive"), pred_->GetSelectionHisto("MHT_JetBin4_HTinclusive"), Title, LumiTitle,"#slash{H}_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/MHT_JetBin4_HTinclusive" + postfix + ".pdf");

        //jet Bin 1
        Title = "2+3 jets, #Delta#phi cut, HT > 500 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("MET_JetBin1_HTinclusive"), pred_->GetSelectionHisto("MET_JetBin1_HTinclusive"), Title, LumiTitle,"#slash{E}_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/MET_JetBin1_HTinclusive" + postfix + ".pdf");

        // jet Bin2
        Title = "4-6 jets, #Delta#phi cut, HT > 500 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("MET_JetBin2_HTinclusive"), pred_->GetSelectionHisto("MET_JetBin2_HTinclusive"), Title, LumiTitle,"#slash{E}_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/MET_JetBin2_HTinclusive" + postfix + ".pdf");

        // jet Bin3
        Title = "7+8 jets, #Delta#phi cut, HT > 500 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("MET_JetBin3_HTinclusive"), pred_->GetSelectionHisto("MET_JetBin3_HTinclusive"), Title, LumiTitle,"#slash{E}_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/MET_JetBin3_HTinclusive" + postfix + ".pdf");

        // jet Bin 4
        Title = ">=9 jets, #Delta#phi cut, HT > 500 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("MET_JetBin4_HTinclusive"), pred_->GetSelectionHisto("MET_JetBin4_HTinclusive"), Title, LumiTitle,"#slash{E}_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/MET_JetBin4_HTinclusive" + postfix + ".pdf");

        // --------------------------------------------------------------------------------------------- //
        // baseline without deltaPhi (HT + MHT)
        //jet Bin 1
        Title = "2+3 jets, HT > 500 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("MHT_JetBin1_baseline_withoutDeltaPhi"), pred_->GetSelectionHisto("MHT_JetBin1_baseline_withoutDeltaPhi"), Title, LumiTitle,"#slash{H}_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/MHT_JetBin1_baseline_withoutDeltaPhi" + postfix + ".pdf");

        // jet Bin2
        Title = "4-6 jets, HT > 500 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("MHT_JetBin2_baseline_withoutDeltaPhi"), pred_->GetSelectionHisto("MHT_JetBin2_baseline_withoutDeltaPhi"), Title, LumiTitle,"#slash{H}_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/MHT_JetBin2_baseline_withoutDeltaPhi" + postfix + ".pdf");

        // jet Bin3
        Title = "7-8 jets, HT > 500 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("MHT_JetBin3_baseline_withoutDeltaPhi"), pred_->GetSelectionHisto("MHT_JetBin3_baseline_withoutDeltaPhi"), Title, LumiTitle,"#slash{H}_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/MHT_JetBin3_baseline_withoutDeltaPhi" + postfix + ".pdf");

        // jet Bin 4
        Title = ">=9 jets, HT > 500 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("MHT_JetBin4_baseline_withoutDeltaPhi"), pred_->GetSelectionHisto("MHT_JetBin4_baseline_withoutDeltaPhi"), Title, LumiTitle,"#slash{H}_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/MHT_JetBin4_baseline_withoutDeltaPhi" + postfix + ".pdf");

        //jet Bin 1
        Title = "2+3 jets, HT > 500 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("MET_JetBin1_baseline_withoutDeltaPhi"), pred_->GetSelectionHisto("MET_JetBin1_baseline_withoutDeltaPhi"), Title, LumiTitle,"#slash{E}_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/MET_JetBin1_baseline_withoutDeltaPhi" + postfix + ".pdf");

        // jet Bin2
        Title = "4-6 jets, HT > 500 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("MET_JetBin2_baseline_withoutDeltaPhi"), pred_->GetSelectionHisto("MET_JetBin2_baseline_withoutDeltaPhi"), Title, LumiTitle,"#slash{E}_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/MET_JetBin2_baseline_withoutDeltaPhi" + postfix + ".pdf");

        // jet Bin3
        Title = "7-8 jets, HT > 500 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("MET_JetBin3_baseline_withoutDeltaPhi"), pred_->GetSelectionHisto("MET_JetBin3_baseline_withoutDeltaPhi"), Title, LumiTitle,"#slash{E}_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/MET_JetBin3_baseline_withoutDeltaPhi" + postfix + ".pdf");

        // jet Bin 4
        Title = ">=9 jets, HT > 500 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("MET_JetBin4_baseline_withoutDeltaPhi"), pred_->GetSelectionHisto("MET_JetBin4_baseline_withoutDeltaPhi"), Title, LumiTitle,"#slash{E}_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/MET_JetBin4_baseline_withoutDeltaPhi" + postfix + ".pdf");

        //jet Bin 1
        Title = "2+3 jets, MET > 200 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("HT_JetBin1_baseline_withoutDeltaPhi"), pred_->GetSelectionHisto("HT_JetBin1_baseline_withoutDeltaPhi"), Title, LumiTitle,"H_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/HT_JetBin1_baseline_withoutDeltaPhi" + postfix + ".pdf");

        // jet Bin2
        Title = "4-6 jets, MET > 200 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("HT_JetBin2_baseline_withoutDeltaPhi"), pred_->GetSelectionHisto("HT_JetBin2_baseline_withoutDeltaPhi"), Title, LumiTitle,"H_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/HT_JetBin2_baseline_withoutDeltaPhi" + postfix + ".pdf");

        // jet Bin3
        Title = "7+8 jets, MET > 200 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("HT_JetBin3_baseline_withoutDeltaPhi"), pred_->GetSelectionHisto("HT_JetBin3_baseline_withoutDeltaPhi"), Title, LumiTitle,"H_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/HT_JetBin3_baseline_withoutDeltaPhi" + postfix + ".pdf");

        // jet Bin 4
        Title = ">=9 jets, MET > 200 GeV";
        c = DrawComparison( pred_->GetPredictionHisto("HT_JetBin4_baseline_withoutDeltaPhi"), pred_->GetSelectionHisto("HT_JetBin4_baseline_withoutDeltaPhi"), Title, LumiTitle,"H_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/HT_JetBin4_baseline_withoutDeltaPhi" + postfix + ".pdf");

        // --------------------------------------------------------------------------------------------- //
        // baseline plots
        Title = "#Delta#phi cut, HT > 500 GeV";
        c =  DrawComparison( pred_->GetPredictionHisto("NJets_baseline_withoutMET"), pred_->GetSelectionHisto("NJets_baseline_withoutMET"), Title, LumiTitle,"N_{Jets}", yTitle, isData);
        c->Print("output_GetPrediction/NJets_baseline_withoutMET" + postfix + ".pdf");

        Title = "#Delta#phi cut, HT > 500 GeV, MET > 200 GeV";
        c =  DrawComparison( pred_->GetPredictionHisto("NJets_baseline"), pred_->GetSelectionHisto("NJets_baseline"), Title, LumiTitle,"N_{Jets}", yTitle, isData);
        c->Print("output_GetPrediction/NJets_baseline" + postfix + ".pdf");

        Title = "HT > 500 GeV";
        c =  DrawComparison( pred_->GetPredictionHisto("NJets_baseline_withoutDeltaPhi_withoutMET"), pred_->GetSelectionHisto("NJets_baseline_withoutDeltaPhi_withoutMET"), Title, LumiTitle,"N_{Jets}", yTitle, isData);
        c->Print("output_GetPrediction/NJets_baseline_withoutDeltaPhi_withoutMET" + postfix + ".pdf");

        Title = "HT > 500 GeV, MET > 200 GeV";
        c =  DrawComparison( pred_->GetPredictionHisto("NJets_baseline_withoutDeltaPhi"), pred_->GetSelectionHisto("NJets_baseline_withoutDeltaPhi"), Title, LumiTitle,"N_{Jets}", yTitle, isData);
        c->Print("output_GetPrediction/NJets_baseline_withoutDeltaPhi" + postfix + ".pdf");

        Title = "#Delta#phi cut, HT > 500 GeV";
        c =  DrawComparison( pred_->GetPredictionHisto("NBJets_baseline_withoutMET"), pred_->GetSelectionHisto("NBJets_baseline_withoutMET"), Title, LumiTitle,"N_{b-Tags}", yTitle, isData);
        c->Print("output_GetPrediction/NBJets_baseline_withoutMET" + postfix + ".pdf");

        Title = "#Delta#phi cut, HT > 500 GeV, MET > 200 GeV";
        c =  DrawComparison( pred_->GetPredictionHisto("NBJets_baseline"), pred_->GetSelectionHisto("NBJets_baseline"), Title, LumiTitle,"N_{b-Tags}", yTitle, isData);
        c->Print("output_GetPrediction/NBJets_baseline" + postfix + ".pdf");

        Title = "HT > 500 GeV";
        c =  DrawComparison( pred_->GetPredictionHisto("NBJets_baseline_withoutDeltaPhi_withoutMET"), pred_->GetSelectionHisto("NBJets_baseline_withoutDeltaPhi_withoutMET"), Title, LumiTitle,"N_{b-Tags}", yTitle, isData);
        c->Print("output_GetPrediction/NBJets_baseline_withoutDeltaPhi_withoutMET" + postfix + ".pdf");

        Title = "HT > 500 GeV, MET > 200 GeV";
        c =  DrawComparison( pred_->GetPredictionHisto("NBJets_baseline_withoutDeltaPhi"), pred_->GetSelectionHisto("NBJets_baseline_withoutDeltaPhi"), Title, LumiTitle,"N_{b-Tags}", yTitle, isData);
        c->Print("output_GetPrediction/NBJets_baseline_withoutDeltaPhi" + postfix + ".pdf");

        Title = ">=4 jets, #Delta#phi cut, HT > 500 GeV";
        c =  DrawComparison( pred_->GetPredictionHisto("MHT_baseline"), pred_->GetSelectionHisto("MHT_baseline"), Title, LumiTitle,"#slash{H}_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/MHT_baseline" + postfix + ".pdf");

        Title = ">=4 jets, #Delta#phi cut, HT > 500 GeV";
        c =  DrawComparison( pred_->GetPredictionHisto("MET_baseline"), pred_->GetSelectionHisto("MET_baseline"), Title, LumiTitle,"#slash{E}_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/MET_baseline" + postfix + ".pdf");

        Title = ">=4 jets, #Delta#phi cut, MET > 200 GeV";
        c =  DrawComparison( pred_->GetPredictionHisto("HT_baseline"), pred_->GetSelectionHisto("HT_baseline"), Title, LumiTitle,"H_{T} (GeV)", yTitle, isData);
        c->Print("output_GetPrediction/HT_baseline" + postfix + ".pdf");
        // --------------------------------------------------------------------------------------------- //

    }

    return 1;
}

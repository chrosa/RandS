#include <TSystem.h>

#include <TProfile.h>
#include <TF1.h>
#include <TArrayF.h>

#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TLegend.h>

#include <TPostScript.h>
#include <TString.h>

#include <TMath.h>

#include <memory>
#include <string>
#include <cassert>
#include <cmath>
#include <iostream>

#include "Prediction.h"

using namespace std;

Prediction::Prediction(TChain& QCDPrediction, TString postfix)
{
    gROOT->ProcessLine("#include <vector>");

    // ------------- define all histos needed -------//
    // set histogram attributes
    int Npseudo = 20;
    int NbinsMHT = 100;
    int NbinsMHTsoft = 100;
    int NbinsMHTsig = 100;
    int NbinsHT = 100;
    int NbinsJetPt = 100;
    int NbinsJetEta = 100;
    int NbinsJetPhi = 100;
    double MHTsigmin = 0.;
    double MHTsigmax = 20.;
    double MHTsoftmin = 0.;
    double MHTsoftmax = 100.;
    double MHTmin = 0.;
    double MHTmax = 1000.;
    double HTmin = 0.;
    double HTmax = 5000.;
    double JetPtmin = 0.;
    double JetPtmax = 2500.;
    double JetVBFPtmin = 0.;
    double JetVBFPtmax = 500.;
    double JetEtamin = -5.;
    double JetEtamax = 5.;

    double HTSave = 0;
    double METSave = 100;
    double MHTjjSave = 9999.;
    double MHTSave = 9999.;
    double MjjSave = 0.;
    double dPhijjSave = 2.7;
    double dEtajjSave = 2.5;
    int NJetsSave = 0;

    double DEtaLoose = 3.0;
    double DEtaTight = 4.8;
    double DEtajj = 2.5;
    double DPhiSR = 1.8;
    double DPhiCRMax = 2.7;
    double DPhiCRMin = 1.8;
    double MjjCut = 1000.;
    double METCut = 150.;
    
    double MHTSigSeedMax = 5.;
    double METSigSeedMax = 999999.;
    double METSoftSeedMax = 30.;

    //double MHTSigSeedMax = 999999.;
    //double METSigSeedMax = 999999.;
    //double METSoftSeedMax = 999999.;

    bool blindSR = true;
    bool VBF = true;
    bool HTMHT = false;

    // define prediction histograms
    // preselection
    HT_presel_pred_raw = new TH2F("presel_HT_prediction", "presel_HT_prediction", NbinsHT, HTmin, HTmax, Npseudo, 0.5, Npseudo + 0.5);
    MHT_presel_pred_raw = new TH2F("presel_MHT_prediction", "presel_MHT_prediction", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    MET_presel_pred_raw = new TH2F("presel_MET_prediction", "presel_MET_prediction", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet1Pt_presel_pred_raw = new TH2F("presel_Jet1_Pt_prediction", "presel_Jet1_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Pt_presel_pred_raw = new TH2F("presel_Jet2_Pt_prediction", "presel_Jet2_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet3Pt_presel_pred_raw = new TH2F("presel_Jet3_Pt_prediction", "presel_Jet3_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet1Eta_presel_pred_raw = new TH2F("presel_Jet1_Eta_prediction", "presel_Jet1_Eta", NbinsJetEta,  JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Eta_presel_pred_raw = new TH2F("presel_Jet2_Eta_prediction", "presel_Jet2_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    Jet3Eta_presel_pred_raw = new TH2F("presel_Jet3_Eta_prediction", "presel_Jet3_Eta", NbinsJetEta,  JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi1_presel_pred_raw = new TH2F("presel_DeltaPhi1_prediction", "DeltaPhi1", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi2_presel_pred_raw = new TH2F("presel_DeltaPhi2_prediction", "DeltaPhi2", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi3_presel_pred_raw = new TH2F("presel_DeltaPhi3_prediction", "DeltaPhi3", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);

    // preselection + delta Phi
    HT_deltaPhi_pred_raw = new TH2F("deltaPhi_HT_prediction", "deltaPhi_HT_prediction", NbinsHT, HTmin, HTmax, Npseudo, 0.5, Npseudo + 0.5);
    MHT_deltaPhi_pred_raw = new TH2F("deltaPhi_MHT_prediction", "deltaPhi_MHT_prediction", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    MET_deltaPhi_pred_raw = new TH2F("deltaPhi_MET_prediction", "deltaPhi_MET_prediction", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet1Pt_deltaPhi_pred_raw = new TH2F("deltaPhi_Jet1_Pt_prediction", "deltaPhi_Jet1_Pt", NbinsJetPt,  JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Pt_deltaPhi_pred_raw = new TH2F("deltaPhi_Jet2_Pt_prediction", "deltaPhi_Jet2_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet3Pt_deltaPhi_pred_raw = new TH2F("deltaPhi_Jet3_Pt_prediction", "deltaPhi_Jet3_Pt", NbinsJetPt,  JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet1Eta_deltaPhi_pred_raw = new TH2F("deltaPhi_Jet1_Eta_prediction", "deltaPhi_Jet1_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Eta_deltaPhi_pred_raw = new TH2F("deltaPhi_Jet2_Eta_prediction", "deltaPhi_Jet2_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    Jet3Eta_deltaPhi_pred_raw = new TH2F("deltaPhi_Jet3_Eta_prediction", "deltaPhi_Jet3_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);

    // HT inclusive, baseline
    MHT_JetBin1_HTinclusive_pred_raw = new TH2F("MHT_JetBin1_HTinclusive_pred", "MHT_JetBin1_HTinclusive_pred", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    MHT_JetBin2_HTinclusive_pred_raw = new TH2F("MHT_JetBin2_HTinclusive_pred", "MHT_JetBin2_HTinclusive_pred", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    MHT_JetBin3_HTinclusive_pred_raw = new TH2F("MHT_JetBin3_HTinclusive_pred", "MHT_JetBin3_HTinclusive_pred", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    MHT_JetBin4_HTinclusive_pred_raw = new TH2F("MHT_JetBin4_HTinclusive_pred", "MHT_JetBin4_HTinclusive_pred", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);

    MET_JetBin1_HTinclusive_pred_raw = new TH2F("MET_JetBin1_HTinclusive_pred", "MET_JetBin1_HTinclusive_pred", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    MET_JetBin2_HTinclusive_pred_raw = new TH2F("MET_JetBin2_HTinclusive_pred", "MET_JetBin2_HTinclusive_pred", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    MET_JetBin3_HTinclusive_pred_raw = new TH2F("MET_JetBin3_HTinclusive_pred", "MET_JetBin3_HTinclusive_pred", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    MET_JetBin4_HTinclusive_pred_raw = new TH2F("MET_JetBin4_HTinclusive_pred", "MET_JetBin4_HTinclusive_pred", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);

    // baseline
    HT_baseline_pred_raw = new TH2F("HT_baseline_pred", "HT baseline", NbinsHT, HTmin, HTmax, Npseudo, 0.5, Npseudo + 0.5);
    MHT_baseline_pred_raw = new TH2F("MHT_baseline_pred", "MHT baseline", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    MET_baseline_pred_raw = new TH2F("MET_baseline_pred", "MET baseline", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);

    // baseline jet bin 1
    Jet1Pt_JetBin1_baseline_pred_raw = new TH2F("baseline_Jet1_Pt_JetBin1_prediction", "Jet1_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Pt_JetBin1_baseline_pred_raw = new TH2F("baseline_Jet2_Pt_JetBin1_prediction", "Jet2_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet1Eta_JetBin1_baseline_pred_raw = new TH2F("baseline_Jet1_Eta_JetBin1_prediction", "Jet1_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Eta_JetBin1_baseline_pred_raw = new TH2F("baseline_Jet2_Eta_JetBin1_prediction", "Jet2_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi1_JetBin1_baseline_pred_raw = new TH2F("baseline_DeltaPhi1_JetBin1_prediction", "DeltaPhi1", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi2_JetBin1_baseline_pred_raw = new TH2F("baseline_DeltaPhi2_JetBin1_prediction", "DeltaPhi2", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);

    // baseline jet bin 2
    Jet1Pt_JetBin2_baseline_pred_raw = new TH2F("baseline_Jet1_Pt_JetBin2_prediction", "Jet1_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Pt_JetBin2_baseline_pred_raw = new TH2F("baseline_Jet2_Pt_JetBin2_prediction", "Jet2_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet3Pt_JetBin2_baseline_pred_raw = new TH2F("baseline_Jet3_Pt_JetBin2_prediction", "Jet3_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet1Eta_JetBin2_baseline_pred_raw = new TH2F("baseline_Jet1_Eta_JetBin2_prediction", "Jet1_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Eta_JetBin2_baseline_pred_raw = new TH2F("baseline_Jet2_Eta_JetBin2_prediction", "Jet2_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    Jet3Eta_JetBin2_baseline_pred_raw = new TH2F("baseline_Jet3_Eta_JetBin2_prediction", "Jet3_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi1_JetBin2_baseline_pred_raw = new TH2F("baseline_DeltaPhi1_JetBin2_prediction", "DeltaPhi1", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi2_JetBin2_baseline_pred_raw = new TH2F("baseline_DeltaPhi2_JetBin2_prediction", "DeltaPhi2", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi3_JetBin2_baseline_pred_raw = new TH2F("baseline_DeltaPhi3_JetBin2_prediction", "DeltaPhi3", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);

    // baseline jet bin 3
    Jet1Pt_JetBin3_baseline_pred_raw = new TH2F("baseline_Jet1_Pt_JetBin3_prediction", "Jet1_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Pt_JetBin3_baseline_pred_raw = new TH2F("baseline_Jet2_Pt_JetBin3_prediction", "Jet2_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet3Pt_JetBin3_baseline_pred_raw = new TH2F("baseline_Jet3_Pt_JetBin3_prediction", "Jet3_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet1Eta_JetBin3_baseline_pred_raw = new TH2F("baseline_Jet1_Eta_JetBin3_prediction", "Jet1_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Eta_JetBin3_baseline_pred_raw = new TH2F("baseline_Jet2_Eta_JetBin3_prediction", "Jet2_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    Jet3Eta_JetBin3_baseline_pred_raw = new TH2F("baseline_Jet3_Eta_JetBin3_prediction", "Jet3_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi1_JetBin3_baseline_pred_raw = new TH2F("baseline_DeltaPhi1_JetBin3_prediction", "DeltaPhi1", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi2_JetBin3_baseline_pred_raw = new TH2F("baseline_DeltaPhi2_JetBin3_prediction", "DeltaPhi2", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi3_JetBin3_baseline_pred_raw = new TH2F("baseline_DeltaPhi3_JetBin3_prediction", "DeltaPhi3", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);

    // baseline jet bin 4
    Jet1Pt_JetBin4_baseline_pred_raw = new TH2F("baseline_Jet1_Pt_JetBin4_prediction", "Jet1_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Pt_JetBin4_baseline_pred_raw = new TH2F("baseline_Jet2_Pt_JetBin4_prediction", "Jet2_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet3Pt_JetBin4_baseline_pred_raw = new TH2F("baseline_Jet3_Pt_JetBin4_prediction", "Jet3_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet1Eta_JetBin4_baseline_pred_raw = new TH2F("baseline_Jet1_Eta_JetBin4_prediction", "Jet1_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Eta_JetBin4_baseline_pred_raw = new TH2F("baseline_Jet2_Eta_JetBin4_prediction", "Jet2_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    Jet3Eta_JetBin4_baseline_pred_raw = new TH2F("baseline_Jet3_Eta_JetBin4_prediction", "Jet3_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi1_JetBin4_baseline_pred_raw = new TH2F("baseline_DeltaPhi1_JetBin4_prediction", "DeltaPhi1", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi2_JetBin4_baseline_pred_raw = new TH2F("baseline_DeltaPhi2_JetBin4_prediction", "DeltaPhi2", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi3_JetBin4_baseline_pred_raw = new TH2F("baseline_DeltaPhi3_JetBin4_prediction", "DeltaPhi3", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);

    // baseline without delta Phi jet bin 1
    HT_JetBin1_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_HT_JetBin1_prediction", "HT_prediction", NbinsHT, HTmin, HTmax, Npseudo, 0.5, Npseudo + 0.5);
    MHT_JetBin1_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_MHT_JetBin1_prediction", "MHT_prediction", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    MET_JetBin1_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_MET_JetBin1_prediction", "MET_prediction", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet1Pt_JetBin1_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet1_Pt_JetBin1_prediction", "Jet1_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Pt_JetBin1_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet2_Pt_JetBin1_prediction", "Jet2_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet1Eta_JetBin1_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet1_Eta_JetBin1_prediction", "Jet1_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Eta_JetBin1_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet2_Eta_JetBin1_prediction", "Jet2_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi1_JetBin1_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_DeltaPhi1_JetBin1_prediction", "DeltaPhi1", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi2_JetBin1_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_DeltaPhi2_JetBin1_prediction", "DeltaPhi2", NbinsJetPhi,  0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);

    // baseline without delta Phi jet bin 2
    HT_JetBin2_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_HT_JetBin2_prediction", "HT_prediction", NbinsHT, HTmin, HTmax, Npseudo, 0.5, Npseudo + 0.5);
    MHT_JetBin2_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_MHT_JetBin2_prediction", "MHT_prediction", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    MET_JetBin2_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_MET_JetBin2_prediction", "MET_prediction", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet1Pt_JetBin2_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet1_Pt_JetBin2_prediction", "Jet1_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Pt_JetBin2_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet2_Pt_JetBin2_prediction", "Jet2_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet3Pt_JetBin2_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet3_Pt_JetBin2_prediction", "Jet3_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet1Eta_JetBin2_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet1_Eta_JetBin2_prediction", "Jet1_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Eta_JetBin2_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet2_Eta_JetBin2_prediction", "Jet2_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    Jet3Eta_JetBin2_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet3_Eta_JetBin2_prediction", "Jet3_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi1_JetBin2_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_DeltaPhi1_JetBin2_prediction", "DeltaPhi1", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi2_JetBin2_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_DeltaPhi2_JetBin2_prediction", "DeltaPhi2", NbinsJetPhi,  0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi3_JetBin2_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_DeltaPhi3_JetBin2_prediction", "DeltaPhi3", NbinsJetPhi,  0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);

    // baseline without delta Phi jet bin 3
    HT_JetBin3_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_HT_JetBin3_prediction", "HT_prediction", NbinsHT, HTmin, HTmax, Npseudo, 0.5, Npseudo + 0.5);
    MHT_JetBin3_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_MHT_JetBin3_prediction", "MHT_prediction", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    MET_JetBin3_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_MET_JetBin3_prediction", "MET_prediction", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet1Pt_JetBin3_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet1_Pt_JetBin3_prediction", "Jet1_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Pt_JetBin3_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet2_Pt_JetBin3_prediction", "Jet2_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet3Pt_JetBin3_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet3_Pt_JetBin3_prediction", "Jet3_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet1Eta_JetBin3_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet1_Eta_JetBin3_prediction", "Jet1_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Eta_JetBin3_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet2_Eta_JetBin3_prediction", "Jet2_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    Jet3Eta_JetBin3_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet3_Eta_JetBin3_prediction", "Jet3_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi1_JetBin3_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_DeltaPhi1_JetBin3_prediction", "DeltaPhi1", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi2_JetBin3_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_DeltaPhi2_JetBin3_prediction", "DeltaPhi2", NbinsJetPhi,  0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi3_JetBin3_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_DeltaPhi3_JetBin3_prediction", "DeltaPhi3", NbinsJetPhi,  0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);

    // baseline without delta Phi jet bin 4
    HT_JetBin4_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_HT_JetBin4_prediction", "HT_prediction", NbinsHT, HTmin, HTmax, Npseudo, 0.5, Npseudo + 0.5);
    MHT_JetBin4_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_MHT_JetBin4_prediction", "MHT_prediction", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    MET_JetBin4_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_MET_JetBin4_prediction", "MET_prediction", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet1Pt_JetBin4_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet1_Pt_JetBin4_prediction", "Jet1_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Pt_JetBin4_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet2_Pt_JetBin4_prediction", "Jet2_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet3Pt_JetBin4_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet3_Pt_JetBin4_prediction", "Jet3_Pt", NbinsJetPt, JetPtmin, JetPtmax, Npseudo, 0.5, Npseudo + 0.5);
    Jet1Eta_JetBin4_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet1_Eta_JetBin4_prediction", "Jet1_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    Jet2Eta_JetBin4_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet2_Eta_JetBin4_prediction", "Jet2_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    Jet3Eta_JetBin4_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_Jet3_Eta_JetBin4_prediction", "Jet3_Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi1_JetBin4_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_DeltaPhi1_JetBin4_prediction", "DeltaPhi1", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi2_JetBin4_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_DeltaPhi2_JetBin4_prediction", "DeltaPhi2", NbinsJetPhi,  0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    DeltaPhi3_JetBin4_baseline_withoutDeltaPhi_pred_raw = new TH2F("baseline_withoutDeltaPhi_DeltaPhi3_JetBin4_prediction", "DeltaPhi3", NbinsJetPhi,  0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);

    // Mjj
    VBF_dPhi_presel_pred_raw = new TH2F("VBF_dPhi_presel_prediction","dPhi", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    VBF_dEta_presel_pred_raw = new TH2F("VBF_dEta_presel_prediction","dEta", NbinsJetEta, 0., 10., Npseudo, 0.5, Npseudo + 0.5);
    VBF_Mjj_presel_pred_raw = new TH2F("VBF_Mjj_presel_prediction","Mjj", NbinsHT, HTmin, HTmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet1Pt_presel_pred_raw = new TH2F("VBF_Jet1Pt_presel_prediction","Jet1Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet2Pt_presel_pred_raw = new TH2F("VBF_Jet2Pt_presel_prediction","Jet2Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet3Pt_presel_pred_raw = new TH2F("VBF_Jet3Pt_presel_prediction","Jet3Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax/2., Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet1Eta_presel_pred_raw = new TH2F("VBF_Jet1Eta_presel_prediction","Jet1Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet2Eta_presel_pred_raw = new TH2F("VBF_Jet2Eta_presel_prediction","Jet2Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet3Eta_presel_pred_raw = new TH2F("VBF_Jet3Eta_presel_prediction","Jet3Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_PTjj_presel_pred_raw = new TH2F("VBF_PTjj_presel_prediction","PTjj", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_MET_presel_pred_raw = new TH2F("VBF_MET_presel_prediction","MET", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_METsoft_presel_pred_raw = new TH2F("VBF_METsoft_presel_prediction","METsoft", NbinsMHTsoft, MHTsoftmin, MHTsoftmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_METsig_presel_pred_raw = new TH2F("VBF_METsig_presel_prediction","METsig", NbinsMHTsig, MHTsigmin, MHTsigmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_MHTsig_presel_pred_raw = new TH2F("VBF_MHTsig_presel_prediction","MHTsig", NbinsMHTsig, MHTsigmin, MHTsigmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_minDeltaPhiPTj12_presel_pred_raw = new TH2F("VBF_minDeltaPhiPTj12_presel_prediction","minDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    VBF_maxDeltaPhiPTj12_presel_pred_raw = new TH2F("VBF_maxDeltaPhiPTj12_presel_prediction","maxDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    VBF_DeltaPhiPTj3_presel_pred_raw = new TH2F("VBF_DeltaPhiPTj3_presel_prediction","DeltaPhiPTj3", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);

    VBF_dPhi_presel_4JV_dPhiSide_pred_raw = new TH2F("VBF_dPhi_presel_4JV_dPhiSide_prediction","dPhi", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    VBF_dEta_presel_4JV_dPhiSide_pred_raw = new TH2F("VBF_dEta_presel_4JV_dPhiSide_prediction","dEta", NbinsJetEta, 0., 10., Npseudo, 0.5, Npseudo + 0.5);
    VBF_Mjj_presel_4JV_dPhiSide_pred_raw = new TH2F("VBF_Mjj_presel_4JV_dPhiSide_prediction","Mjj", NbinsHT, HTmin, HTmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet1Pt_presel_4JV_dPhiSide_pred_raw = new TH2F("VBF_Jet1Pt_presel_4JV_dPhiSide_prediction","Jet1Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet2Pt_presel_4JV_dPhiSide_pred_raw = new TH2F("VBF_Jet2Pt_presel_4JV_dPhiSide_prediction","Jet2Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet3Pt_presel_4JV_dPhiSide_pred_raw = new TH2F("VBF_Jet3Pt_presel_4JV_dPhiSide_prediction","Jet3Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax/2., Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet1Eta_presel_4JV_dPhiSide_pred_raw = new TH2F("VBF_Jet1Eta_presel_4JV_dPhiSide_prediction","Jet1Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet2Eta_presel_4JV_dPhiSide_pred_raw = new TH2F("VBF_Jet2Eta_presel_4JV_dPhiSide_prediction","Jet2Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet3Eta_presel_4JV_dPhiSide_pred_raw = new TH2F("VBF_Jet3Eta_presel_4JV_dPhiSide_prediction","Jet3Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_PTjj_presel_4JV_dPhiSide_pred_raw = new TH2F("VBF_PTjj_presel_4JV_dPhiSide_prediction","PTjj", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_MET_presel_4JV_dPhiSide_pred_raw = new TH2F("VBF_MET_presel_4JV_dPhiSide_prediction","MET", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_METsoft_presel_4JV_dPhiSide_pred_raw = new TH2F("VBF_METsoft_presel_4JV_dPhiSide_prediction","METsoft", NbinsMHTsoft, MHTsoftmin, MHTsoftmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_METsig_presel_4JV_dPhiSide_pred_raw = new TH2F("VBF_METsig_presel_4JV_dPhiSide_prediction","METsig", NbinsMHTsig, MHTsigmin, MHTsigmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_MHTsig_presel_4JV_dPhiSide_pred_raw = new TH2F("VBF_MHTsig_presel_4JV_dPhiSide_prediction","MHTsig", NbinsMHTsig, MHTsigmin, MHTsigmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_minDeltaPhiPTj12_presel_4JV_dPhiSide_pred_raw = new TH2F("VBF_minDeltaPhiPTj12_presel_4JV_dPhiSide_prediction","minDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    VBF_maxDeltaPhiPTj12_presel_4JV_dPhiSide_pred_raw = new TH2F("VBF_maxDeltaPhiPTj12_presel_4JV_dPhiSide_prediction","maxDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    VBF_DeltaPhiPTj3_presel_4JV_dPhiSide_pred_raw = new TH2F("VBF_DeltaPhiPTj3_presel_4JV_dPhiSide_prediction","DeltaPhiPTj3", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);

    VBF_dPhi_dEta_pred_raw = new TH2F("VBF_dPhi_dEta_prediction","dPhi", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    VBF_dEta_dEta_pred_raw = new TH2F("VBF_dEta_dEta_prediction","dEta", NbinsJetEta, 0., 10., Npseudo, 0.5, Npseudo + 0.5);
    VBF_Mjj_dEta_pred_raw = new TH2F("VBF_Mjj_dEta_prediction","Mjj", NbinsHT, HTmin, HTmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet1Pt_dEta_pred_raw = new TH2F("VBF_Jet1Pt_dEta_prediction","Jet1Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet2Pt_dEta_pred_raw = new TH2F("VBF_Jet2Pt_dEta_prediction","Jet2Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet3Pt_dEta_pred_raw = new TH2F("VBF_Jet3Pt_dEta_prediction","Jet3Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax/2., Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet1Eta_dEta_pred_raw = new TH2F("VBF_Jet1Eta_dEta_prediction","Jet1Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet2Eta_dEta_pred_raw = new TH2F("VBF_Jet2Eta_dEta_prediction","Jet2Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet3Eta_dEta_pred_raw = new TH2F("VBF_Jet3Eta_dEta_prediction","Jet3Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_PTjj_dEta_pred_raw = new TH2F("VBF_PTjj_dEta_prediction","PTjj", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_MET_dEta_pred_raw = new TH2F("VBF_MET_dEta_prediction","MET", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_METsoft_dEta_pred_raw = new TH2F("VBF_METsoft_dEta_prediction","METsoft", NbinsMHTsoft, MHTsoftmin, MHTsoftmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_METsig_dEta_pred_raw = new TH2F("VBF_METsig_dEta_prediction","METsig", NbinsMHTsig, MHTsigmin, MHTsigmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_MHTsig_dEta_pred_raw = new TH2F("VBF_MHTssig_dEta_prediction","MHTsig", NbinsMHTsig, MHTsigmin, MHTsigmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_minDeltaPhiPTj12_dEta_pred_raw = new TH2F("VBF_minDeltaPhiPTj12_dEta_prediction","minDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    VBF_maxDeltaPhiPTj12_dEta_pred_raw = new TH2F("VBF_maxDeltaPhiPTj12_dEta_prediction","maxDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    VBF_DeltaPhiPTj3_dEta_pred_raw = new TH2F("VBF_DeltaPhiPTj3_dEta_prediction","DeltaPhiPTj3", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);

    VBF_dPhi_jj_pred_raw = new TH2F("VBF_dPhi_jj_prediction","dPhi", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    VBF_dEta_jj_pred_raw = new TH2F("VBF_dEta_jj_prediction","dEta", NbinsJetEta, 0., 10., Npseudo, 0.5, Npseudo + 0.5);
    VBF_Mjj_jj_pred_raw = new TH2F("VBF_Mjj_jj_prediction","Mjj", NbinsHT, HTmin, HTmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet1Pt_jj_pred_raw = new TH2F("VBF_Jet1Pt_jj_prediction","Jet1Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet2Pt_jj_pred_raw = new TH2F("VBF_Jet2Pt_jj_prediction","Jet2Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet3Pt_jj_pred_raw = new TH2F("VBF_Jet3Pt_jj_prediction","Jet3Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax/2., Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet1Eta_jj_pred_raw = new TH2F("VBF_Jet1Eta_jj_prediction","Jet1Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet2Eta_jj_pred_raw = new TH2F("VBF_Jet2Eta_jj_prediction","Jet2Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet3Eta_jj_pred_raw = new TH2F("VBF_Jet3Eta_jj_prediction","Jet3Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_PTjj_jj_pred_raw = new TH2F("VBF_PTjj_jj_prediction","PTjj", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_MET_jj_pred_raw = new TH2F("VBF_MET_jj_prediction","MET", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_METsoft_jj_pred_raw = new TH2F("VBF_METsoft_jj_prediction","METsoft", NbinsMHTsoft, MHTsoftmin, MHTsoftmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_METsig_jj_pred_raw = new TH2F("VBF_METsig_jj_prediction","METsig", NbinsMHTsig, MHTsigmin, MHTsigmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_MHTsig_jj_pred_raw = new TH2F("VBF_MHTsig_jj_prediction","METsig", NbinsMHTsig, MHTsigmin, MHTsigmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_minDeltaPhiPTj12_jj_pred_raw = new TH2F("VBF_minDeltaPhiPTj12_jj_prediction","minDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    VBF_maxDeltaPhiPTj12_jj_pred_raw = new TH2F("VBF_maxDeltaPhiPTj12_jj_prediction","maxDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    VBF_DeltaPhiPTj3_jj_pred_raw = new TH2F("VBF_DeltaPhiPTj3_jj_prediction","DeltaPhiPTj3", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);

    VBF_dPhi_dEta_3JV_pred_raw = new TH2F("VBF_dPhi_dEta_3JV_prediction","dPhi", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    VBF_dEta_dEta_3JV_pred_raw = new TH2F("VBF_dEta_dEta_3JV_prediction","dEta", NbinsJetEta, 0., 10., Npseudo, 0.5, Npseudo + 0.5);
    VBF_Mjj_dEta_3JV_pred_raw = new TH2F("VBF_Mjj_dEta_3JV_prediction","Mjj", NbinsHT, HTmin, HTmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet1Pt_dEta_3JV_pred_raw = new TH2F("VBF_Jet1Pt_dEta_3JV_prediction","Jet1Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet2Pt_dEta_3JV_pred_raw = new TH2F("VBF_Jet2Pt_dEta_3JV_prediction","Jet2Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet3Pt_dEta_3JV_pred_raw = new TH2F("VBF_Jet3Pt_dEta_3JV_prediction","Jet3Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax/2., Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet1Eta_dEta_3JV_pred_raw = new TH2F("VBF_Jet1Eta_dEta_3JV_prediction","Jet1Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet2Eta_dEta_3JV_pred_raw = new TH2F("VBF_Jet2Eta_dEta_3JV_prediction","Jet2Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet3Eta_dEta_3JV_pred_raw = new TH2F("VBF_Jet3Eta_dEta_3JV_prediction","Jet3Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_PTjj_dEta_3JV_pred_raw = new TH2F("VBF_PTjj_dEta_3JV_prediction","PTjj", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_MET_dEta_3JV_pred_raw = new TH2F("VBF_MET_dEta_3JV_prediction","MET", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_METsoft_dEta_3JV_pred_raw = new TH2F("VBF_METsoft_dEta_3JV_prediction","METsoft", NbinsMHTsoft, MHTsoftmin, MHTsoftmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_METsig_dEta_3JV_pred_raw = new TH2F("VBF_METsig_dEta_3JV_prediction","METsig", NbinsMHTsig, MHTsigmin, MHTsigmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_MHTsig_dEta_3JV_pred_raw = new TH2F("VBF_MHTsig_dEta_3JV_prediction","MHTsig", NbinsMHTsig, MHTsigmin, MHTsigmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_minDeltaPhiPTj12_dEta_3JV_pred_raw = new TH2F("VBF_minDeltaPhiPTj12_dEta_3JV_prediction","minDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    VBF_maxDeltaPhiPTj12_dEta_3JV_pred_raw = new TH2F("VBF_maxDeltaPhiPTj12_dEta_3JV_prediction","maxDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    VBF_DeltaPhiPTj3_dEta_3JV_pred_raw = new TH2F("VBF_DeltaPhiPTj3_dEta_3JV_prediction","DeltaPhiPTj3", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);

    VBF_dPhi_dEta_3JV_dPhiPTjj_pred_raw = new TH2F("VBF_dPhi_dEta_3JV_dPhiPTjj_prediction","dPhi", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    VBF_dEta_dEta_3JV_dPhiPTjj_pred_raw = new TH2F("VBF_dEta_dEta_3JV_dPhiPTjj_prediction","dEta", NbinsJetEta, 0., 10., Npseudo, 0.5, Npseudo + 0.5);
    VBF_Mjj_dEta_3JV_dPhiPTjj_pred_raw = new TH2F("VBF_Mjj_dEta_3JV_dPhiPTjj_prediction","Mjj", NbinsHT, HTmin, HTmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet1Pt_dEta_3JV_dPhiPTjj_pred_raw = new TH2F("VBF_Jet1Pt_dEta_3JV_dPhiPTjj_prediction","Jet1Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet2Pt_dEta_3JV_dPhiPTjj_pred_raw = new TH2F("VBF_Jet2Pt_dEta_3JV_dPhiPTjj_prediction","Jet2Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet3Pt_dEta_3JV_dPhiPTjj_pred_raw = new TH2F("VBF_Jet3Pt_dEta_3JV_dPhiPTjj_prediction","Jet3Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax/2., Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet1Eta_dEta_3JV_dPhiPTjj_pred_raw = new TH2F("VBF_Jet1Eta_dEta_3JV_dPhiPTjj_prediction","Jet1Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet2Eta_dEta_3JV_dPhiPTjj_pred_raw = new TH2F("VBF_Jet2Eta_dEta_3JV_dPhiPTjj_prediction","Jet2Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_Jet3Eta_dEta_3JV_dPhiPTjj_pred_raw = new TH2F("VBF_Jet3Eta_dEta_3JV_dPhiPTjj_prediction","Jet3Eta", NbinsJetEta, JetEtamin, JetEtamax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_PTjj_dEta_3JV_dPhiPTjj_pred_raw = new TH2F("VBF_PTjj_dEta_3JV_dPhiPTjj_prediction","PTjj", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_MET_dEta_3JV_dPhiPTjj_pred_raw = new TH2F("VBF_MET_dEta_3JV_dPhiPTjj_prediction","MET", NbinsMHT, MHTmin, MHTmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_METsoft_dEta_3JV_dPhiPTjj_pred_raw = new TH2F("VBF_METsoft_dEta_3JV_dPhiPTjj_prediction","METsoft", NbinsMHTsoft, MHTsoftmin, MHTsoftmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_METsig_dEta_3JV_dPhiPTjj_pred_raw = new TH2F("VBF_METsig_dEta_3JV_dPhiPTjj_prediction","METsig", NbinsMHTsig, MHTsigmin, MHTsigmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_MHTsig_dEta_3JV_dPhiPTjj_pred_raw = new TH2F("VBF_MHTsig_dEta_3JV_dPhiPTjj_prediction","MHTsig", NbinsMHTsig, MHTsigmin, MHTsigmax, Npseudo, 0.5, Npseudo + 0.5);
    VBF_minDeltaPhiPTj12_dEta_3JV_dPhiPTjj_pred_raw = new TH2F("VBF_minDeltaPhiPTj12_dEta_3JV_dPhiPTjj_prediction","minDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    VBF_maxDeltaPhiPTj12_dEta_3JV_dPhiPTjj_pred_raw = new TH2F("VBF_maxDeltaPhiPTj12_dEta_3JV_dPhiPTjj_prediction","maxDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);
    VBF_DeltaPhiPTj3_dEta_3JV_dPhiPTjj_pred_raw = new TH2F("VBF_DeltaPhiPTj3_dEta_3JV_dPhiPTjj_prediction","DeltaPhiPTj3", NbinsJetPhi, 0., TMath::Pi(), Npseudo, 0.5, Npseudo + 0.5);

    // NJets
    NJets_baseline_withoutMET_pred_raw = new TH2F("NJets_baseline_withoutMET_pred", "NJets baseline", 15, 0, 15, Npseudo, 0.5, Npseudo + 0.5);
    NJets_baseline_pred_raw = new TH2F("NJets_baseline_pred", "NJets baseline", 15, 0, 15, Npseudo, 0.5, Npseudo + 0.5);
    NJets_baseline_withoutDeltaPhi_withoutMET_pred_raw = new TH2F("NJets_baseline_withoutDeltaPhi_withoutMET_pred", "NJets baseline", 15, 0, 15, Npseudo, 0.5, Npseudo + 0.5);
    NJets_baseline_withoutDeltaPhi_pred_raw = new TH2F("NJets_baseline_withoutDeltaPhi_pred", "NJets baseline", 15, 0, 15, Npseudo, 0.5, Npseudo + 0.5);
    NJets_presel_pred_raw = new TH2F("NJets_presel_pred_raw", "NJets presel", 15, 0, 15, Npseudo, 0.5, Npseudo + 0.5);

    // NBJets
    NBJets_baseline_withoutMET_pred_raw = new TH2F("NBJets_baseline_withoutMET_pred", "NBJets baseline", 5, 0, 5, Npseudo, 0.5, Npseudo + 0.5);
    NBJets_baseline_pred_raw = new TH2F("NBJets_baseline_pred", "NBJets baseline", 5, 0, 5, Npseudo, 0.5, Npseudo + 0.5);
    NBJets_baseline_withoutDeltaPhi_withoutMET_pred_raw = new TH2F("NBJets_baseline_withoutDeltaPhi_withoutMET_pred", "NBJets baseline", 5, 0, 5, Npseudo, 0.5, Npseudo + 0.5);
    NBJets_baseline_withoutDeltaPhi_pred_raw = new TH2F("NBJets_baseline_withoutDeltaPhi_pred", "NBJets baseline", 5, 0, 5, Npseudo, 0.5, Npseudo + 0.5);
    NBJets_presel_pred_raw = new TH2F("NBJets_presel_pred_raw", "NBJets presel", 5, 0, 5, Npseudo, 0.5, Npseudo + 0.5);

    // ------------------------------------------------------------------------------ //

    // define selection histograms
    // preselection
    HT_presel_sel = new TH1F("presel_HT_selection", "presel_HT_selection", NbinsHT, HTmin, HTmax);
    MHT_presel_sel = new TH1F("presel_MHT_selection", "presel_MHT_selection", NbinsMHT, MHTmin, MHTmax);
    MET_presel_sel = new TH1F("presel_MET_selection", "presel_MET_selection", NbinsMHT, MHTmin, MHTmax);
    Jet1Pt_presel_sel = new TH1F("presel_Jet1_Pt_selection", "presel_Jet1_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet2Pt_presel_sel = new TH1F("presel_Jet2_Pt_selection", "presel_Jet2_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet3Pt_presel_sel = new TH1F("presel_Jet3_Pt_selection", "presel_Jet3_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet1Eta_presel_sel = new TH1F("presel_Jet1_Eta_selection", "presel_Jet1_Eta", NbinsJetEta, JetEtamin, JetEtamax);
    Jet2Eta_presel_sel = new TH1F("presel_Jet2_Eta_selection", "presel_Jet2_Eta", NbinsJetEta, JetEtamin, JetEtamax);
    Jet3Eta_presel_sel = new TH1F("presel_Jet3_Eta_selection", "presel_Jet3_Eta", NbinsJetEta, JetEtamin, JetEtamax);
    DeltaPhi1_presel_sel = new TH1F("presel_DeltaPhi1_selection", "DeltaPhi1", NbinsJetPhi, 0., TMath::Pi());
    DeltaPhi2_presel_sel = new TH1F("presel_DeltaPhi2_selection", "DeltaPhi2", NbinsJetPhi, 0., TMath::Pi());
    DeltaPhi3_presel_sel = new TH1F("presel_DeltaPhi3_selection", "DeltaPhi3", NbinsJetPhi, 0., TMath::Pi());

    // preselection + delta Phi
    HT_deltaPhi_sel = new TH1F("deltaPhi_HT_selection", "deltaPhi_HT_selection", NbinsHT, HTmin, HTmax);
    MHT_deltaPhi_sel = new TH1F("deltaPhi_MHT_selection", "deltaPhi_MHT_selection", NbinsMHT, MHTmin, MHTmax);
    MET_deltaPhi_sel = new TH1F("deltaPhi_MET_selection", "deltaPhi_MET_selection", NbinsMHT, MHTmin, MHTmax);
    Jet1Pt_deltaPhi_sel = new TH1F("deltaPhi_Jet1_Pt_selection", "deltaPhi_Jet1_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet2Pt_deltaPhi_sel = new TH1F("deltaPhi_Jet2_Pt_selection", "deltaPhi_Jet2_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet3Pt_deltaPhi_sel = new TH1F("deltaPhi_Jet3_Pt_selection", "deltaPhi_Jet3_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet1Eta_deltaPhi_sel = new TH1F("deltaPhi_Jet1_Eta_selection", "deltaPhi_Jet1_Eta", NbinsJetEta, JetEtamin, JetEtamax);
    Jet2Eta_deltaPhi_sel = new TH1F("deltaPhi_Jet2_Eta_selection", "deltaPhi_Jet2_Eta", NbinsJetEta, JetEtamin, JetEtamax);
    Jet3Eta_deltaPhi_sel = new TH1F("deltaPhi_Jet3_Eta_selection", "deltaPhi_Jet3_Eta", NbinsJetEta, JetEtamin, JetEtamax);

    // HT inclusive, baseline
    MHT_JetBin1_HTinclusive_sel = new TH1F("MHT_JetBin1_HTinclusive_sel", "MHT_JetBin1_HTinclusive_sel", NbinsMHT, MHTmin, MHTmax);
    MHT_JetBin2_HTinclusive_sel = new TH1F("MHT_JetBin2_HTinclusive_sel", "MHT_JetBin2_HTinclusive_sel", NbinsMHT, MHTmin, MHTmax);
    MHT_JetBin3_HTinclusive_sel = new TH1F("MHT_JetBin3_HTinclusive_sel", "MHT_JetBin3_HTinclusive_sel", NbinsMHT, MHTmin, MHTmax);
    MHT_JetBin4_HTinclusive_sel = new TH1F("MHT_JetBin4_HTinclusive_sel", "MHT_JetBin4_HTinclusive_sel", NbinsMHT, MHTmin, MHTmax);

    MET_JetBin1_HTinclusive_sel = new TH1F("MET_JetBin1_HTinclusive_sel", "MET_JetBin1_HTinclusive_sel", NbinsMHT, MHTmin, MHTmax);
    MET_JetBin2_HTinclusive_sel = new TH1F("MET_JetBin2_HTinclusive_sel", "MET_JetBin2_HTinclusive_sel", NbinsMHT, MHTmin, MHTmax);
    MET_JetBin3_HTinclusive_sel = new TH1F("MET_JetBin3_HTinclusive_sel", "MET_JetBin3_HTinclusive_sel", NbinsMHT, MHTmin, MHTmax);
    MET_JetBin4_HTinclusive_sel = new TH1F("MET_JetBin4_HTinclusive_sel", "MET_JetBin4_HTinclusive_sel", NbinsMHT, MHTmin, MHTmax);

    // baseline
    HT_baseline_sel = new TH1F("HT_baseline_sel", "HT baseline", NbinsHT, HTmin, HTmax);
    MHT_baseline_sel = new TH1F("MHT_baseline_sel", "MHT baseline", NbinsMHT, MHTmin, MHTmax);
    MET_baseline_sel = new TH1F("MET_baseline_sel", "MET baseline", NbinsMHT, MHTmin, MHTmax);

    // baseline jet bin 1
    Jet1Pt_JetBin1_baseline_sel = new TH1F("baseline_Jet1_Pt_JetBin1_selection", "Jet1_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet2Pt_JetBin1_baseline_sel = new TH1F("baseline_Jet2_Pt_JetBin1_selection", "Jet2_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet1Eta_JetBin1_baseline_sel = new TH1F("baseline_Jet1_Eta_JetBin1_selection", "Jet1_Eta", NbinsJetEta, JetEtamin, JetEtamax);
    Jet2Eta_JetBin1_baseline_sel = new TH1F("baseline_Jet2_Eta_JetBin1_selection", "Jet2_Eta", NbinsJetEta, JetEtamin, JetEtamax);
    DeltaPhi1_JetBin1_baseline_sel = new TH1F("baseline_DeltaPhi1_JetBin1_selection", "DeltaPhi1", NbinsJetPhi, 0., TMath::Pi());
    DeltaPhi2_JetBin1_baseline_sel = new TH1F("baseline_DeltaPhi2_JetBin1_selection", "DeltaPhi2", NbinsJetPhi, 0., TMath::Pi());

    // baseline jet bin 2
    Jet1Pt_JetBin2_baseline_sel = new TH1F("baseline_Jet1_Pt_JetBin2_selection", "Jet1_Pt", NbinsJetPt,JetPtmin, JetPtmax);
    Jet2Pt_JetBin2_baseline_sel = new TH1F("baseline_Jet2_Pt_JetBin2_selection", "Jet2_Pt", NbinsJetPt,JetPtmin, JetPtmax);
    Jet3Pt_JetBin2_baseline_sel = new TH1F("baseline_Jet3_Pt_JetBin2_selection", "Jet3_Pt", NbinsJetPt,JetPtmin, JetPtmax);
    Jet1Eta_JetBin2_baseline_sel = new TH1F("baseline_Jet1_Eta_JetBin2_selection", "Jet1_Eta", NbinsJetEta,JetEtamin, JetEtamax);
    Jet2Eta_JetBin2_baseline_sel = new TH1F("baseline_Jet2_Eta_JetBin2_selection", "Jet2_Eta", NbinsJetEta,JetEtamin, JetEtamax);
    Jet3Eta_JetBin2_baseline_sel = new TH1F("baseline_Jet3_Eta_JetBin2_selection", "Jet3_Eta", NbinsJetEta,JetEtamin, JetEtamax);
    DeltaPhi1_JetBin2_baseline_sel = new TH1F("baseline_DeltaPhi1_JetBin2_selection", "DeltaPhi1", NbinsJetPhi,0., TMath::Pi());
    DeltaPhi2_JetBin2_baseline_sel = new TH1F("baseline_DeltaPhi2_JetBin2_selection", "DeltaPhi2", NbinsJetPhi,0., TMath::Pi());
    DeltaPhi3_JetBin2_baseline_sel = new TH1F("baseline_DeltaPhi3_JetBin2_selection", "DeltaPhi3", NbinsJetPhi,0., TMath::Pi());

    // baseline jet bin 3
    Jet1Pt_JetBin3_baseline_sel = new TH1F("baseline_Jet1_Pt_JetBin3_selection", "Jet1_Pt", NbinsJetPt,JetPtmin, JetPtmax);
    Jet2Pt_JetBin3_baseline_sel = new TH1F("baseline_Jet2_Pt_JetBin3_selection", "Jet2_Pt", NbinsJetPt,JetPtmin, JetPtmax);
    Jet3Pt_JetBin3_baseline_sel = new TH1F("baseline_Jet3_Pt_JetBin3_selection", "Jet3_Pt", NbinsJetPt,JetPtmin, JetPtmax);
    Jet1Eta_JetBin3_baseline_sel = new TH1F("baseline_Jet1_Eta_JetBin3_selection", "Jet1_Eta", NbinsJetEta,JetEtamin, JetEtamax);
    Jet2Eta_JetBin3_baseline_sel = new TH1F("baseline_Jet2_Eta_JetBin3_selection", "Jet2_Eta", NbinsJetEta,JetEtamin, JetEtamax);
    Jet3Eta_JetBin3_baseline_sel = new TH1F("baseline_Jet3_Eta_JetBin3_selection", "Jet3_Eta", NbinsJetEta,JetEtamin, JetEtamax);
    DeltaPhi1_JetBin3_baseline_sel = new TH1F("baseline_DeltaPhi1_JetBin3_selection", "DeltaPhi1", NbinsJetPhi,0., TMath::Pi());
    DeltaPhi2_JetBin3_baseline_sel = new TH1F("baseline_DeltaPhi2_JetBin3_selection", "DeltaPhi2", NbinsJetPhi,0., TMath::Pi());
    DeltaPhi3_JetBin3_baseline_sel = new TH1F("baseline_DeltaPhi3_JetBin3_selection", "DeltaPhi3", NbinsJetPhi,0., TMath::Pi());

    // baseline jet bin 4
    Jet1Pt_JetBin4_baseline_sel = new TH1F("baseline_Jet1_Pt_JetBin4_selection", "Jet1_Pt", NbinsJetPt,JetPtmin, JetPtmax);
    Jet2Pt_JetBin4_baseline_sel = new TH1F("baseline_Jet2_Pt_JetBin4_selection", "Jet2_Pt", NbinsJetPt,JetPtmin, JetPtmax);
    Jet3Pt_JetBin4_baseline_sel = new TH1F("baseline_Jet3_Pt_JetBin4_selection", "Jet3_Pt", NbinsJetPt,JetPtmin, JetPtmax);
    Jet1Eta_JetBin4_baseline_sel = new TH1F("baseline_Jet1_Eta_JetBin4_selection", "Jet1_Eta", NbinsJetEta,JetEtamin, JetEtamax);
    Jet2Eta_JetBin4_baseline_sel = new TH1F("baseline_Jet2_Eta_JetBin4_selection", "Jet2_Eta", NbinsJetEta,JetEtamin, JetEtamax);
    Jet3Eta_JetBin4_baseline_sel = new TH1F("baseline_Jet3_Eta_JetBin4_selection", "Jet3_Eta", NbinsJetEta,JetEtamin, JetEtamax);
    DeltaPhi1_JetBin4_baseline_sel = new TH1F("baseline_DeltaPhi1_JetBin4_selection", "DeltaPhi1", NbinsJetPhi,0., TMath::Pi());
    DeltaPhi2_JetBin4_baseline_sel = new TH1F("baseline_DeltaPhi2_JetBin4_selection", "DeltaPhi2", NbinsJetPhi,0., TMath::Pi());
    DeltaPhi3_JetBin4_baseline_sel = new TH1F("baseline_DeltaPhi3_JetBin4_selection", "DeltaPhi3", NbinsJetPhi,0., TMath::Pi());

    // baseline without delta Phi jet bin 1
    HT_JetBin1_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_HT_JetBin1_selection", "HT_selection", NbinsHT, HTmin, HTmax);
    MHT_JetBin1_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_MHT_JetBin1_selection", "MHT_selection", NbinsMHT, MHTmin, MHTmax);
    MET_JetBin1_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_MET_JetBin1_selection", "MET_selection", NbinsMHT, MHTmin, MHTmax);
    Jet1Pt_JetBin1_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet1_Pt_JetBin1_selection", "Jet1_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet2Pt_JetBin1_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet2_Pt_JetBin1_selection", "Jet2_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet1Eta_JetBin1_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet1_Eta_JetBin1_selection", "Jet1_Eta", NbinsJetEta, JetEtamin, JetEtamax);
    Jet2Eta_JetBin1_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet2_Eta_JetBin1_selection", "Jet2_Eta", NbinsJetEta, JetEtamin, JetEtamax);
    DeltaPhi1_JetBin1_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_DeltaPhi1_JetBin1_selection", "DeltaPhi1", NbinsJetPhi, 0., TMath::Pi());
    DeltaPhi2_JetBin1_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_DeltaPhi2_JetBin1_selection", "DeltaPhi2", NbinsJetPhi, 0., TMath::Pi());

    // baseline without delta Phi jet bin 2
    HT_JetBin2_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_HT_JetBin2_selection", "HT_selection", NbinsHT, HTmin, HTmax);
    MHT_JetBin2_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_MHT_JetBin2_selection", "MHT_selection", NbinsMHT, MHTmin, MHTmax);
    MET_JetBin2_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_MET_JetBin2_selection", "MET_selection", NbinsMHT, MHTmin, MHTmax);
    Jet1Pt_JetBin2_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet1_Pt_JetBin2_selection", "Jet1_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet2Pt_JetBin2_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet2_Pt_JetBin2_selection", "Jet2_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet3Pt_JetBin2_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet3_Pt_JetBin2_selection", "Jet3_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet1Eta_JetBin2_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet1_Eta_JetBin2_selection", "Jet1_Eta", NbinsJetEta, JetEtamin, JetEtamax);
    Jet2Eta_JetBin2_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet2_Eta_JetBin2_selection", "Jet2_Eta", NbinsJetEta, JetEtamin, JetEtamax);
    Jet3Eta_JetBin2_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet3_Eta_JetBin2_selection", "Jet3_Eta", NbinsJetEta, JetEtamin, JetEtamax);
    DeltaPhi1_JetBin2_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_DeltaPhi1_JetBin2_selection", "DeltaPhi1", NbinsJetPhi, 0., TMath::Pi());
    DeltaPhi2_JetBin2_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_DeltaPhi2_JetBin2_selection", "DeltaPhi2", NbinsJetPhi, 0., TMath::Pi());
    DeltaPhi3_JetBin2_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_DeltaPhi3_JetBin2_selection", "DeltaPhi3", NbinsJetPhi, 0., TMath::Pi());

    // baseline without delta Phi jet bin 3
    HT_JetBin3_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_HT_JetBin3_selection", "HT_selection", NbinsHT, HTmin, HTmax);
    MHT_JetBin3_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_MHT_JetBin3_selection", "MHT_selection", NbinsMHT, MHTmin, MHTmax);
    MET_JetBin3_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_MET_JetBin3_selection", "MET_selection", NbinsMHT, MHTmin, MHTmax);
    Jet1Pt_JetBin3_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet1_Pt_JetBin3_selection", "Jet1_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet2Pt_JetBin3_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet2_Pt_JetBin3_selection", "Jet2_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet3Pt_JetBin3_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet3_Pt_JetBin3_selection", "Jet3_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet1Eta_JetBin3_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet1_Eta_JetBin3_selection", "Jet1_Eta", NbinsJetEta, JetEtamin, JetEtamax);
    Jet2Eta_JetBin3_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet2_Eta_JetBin3_selection", "Jet2_Eta", NbinsJetEta, JetEtamin, JetEtamax);
    Jet3Eta_JetBin3_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet3_Eta_JetBin3_selection", "Jet3_Eta", NbinsJetEta, JetEtamin, JetEtamax);
    DeltaPhi1_JetBin3_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_DeltaPhi1_JetBin3_selection", "DeltaPhi1", NbinsJetPhi, 0., TMath::Pi());
    DeltaPhi2_JetBin3_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_DeltaPhi2_JetBin3_selection", "DeltaPhi2", NbinsJetPhi, 0., TMath::Pi());
    DeltaPhi3_JetBin3_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_DeltaPhi3_JetBin3_selection", "DeltaPhi3", NbinsJetPhi, 0., TMath::Pi());

    // baseline without delta Phi jet bin 4
    HT_JetBin4_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_HT_JetBin4_selection", "HT_selection", NbinsHT, HTmin, HTmax);
    MHT_JetBin4_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_MHT_JetBin4_selection", "MHT_selection", NbinsMHT, MHTmin, MHTmax);
    MET_JetBin4_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_MET_JetBin4_selection", "MET_selection", NbinsMHT, MHTmin, MHTmax);
    Jet1Pt_JetBin4_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet1_Pt_JetBin4_selection", "Jet1_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet2Pt_JetBin4_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet2_Pt_JetBin4_selection", "Jet2_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet3Pt_JetBin4_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet3_Pt_JetBin4_selection", "Jet3_Pt", NbinsJetPt, JetPtmin, JetPtmax);
    Jet1Eta_JetBin4_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet1_Eta_JetBin4_selection", "Jet1_Eta", NbinsJetEta, JetEtamin, JetEtamax);
    Jet2Eta_JetBin4_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet2_Eta_JetBin4_selection", "Jet2_Eta", NbinsJetEta, JetEtamin, JetEtamax);
    Jet3Eta_JetBin4_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_Jet3_Eta_JetBin4_selection", "Jet3_Eta", NbinsJetEta, JetEtamin, JetEtamax);
    DeltaPhi1_JetBin4_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_DeltaPhi1_JetBin4_selection", "DeltaPhi1", NbinsJetPhi, 0., TMath::Pi());
    DeltaPhi2_JetBin4_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_DeltaPhi2_JetBin4_selection", "DeltaPhi2", NbinsJetPhi, 0., TMath::Pi());
    DeltaPhi3_JetBin4_baseline_withoutDeltaPhi_sel = new TH1F("baseline_withoutDeltaPhi_DeltaPhi3_JetBin4_selection", "DeltaPhi3", NbinsJetPhi, 0., TMath::Pi());

    // Mjj
    VBF_dPhi_presel_sel = new TH1F("VBF_dPhi_presel_selection","dPhi", NbinsJetPhi, 0., TMath::Pi());
    VBF_dEta_presel_sel = new TH1F("VBF_dEta_presel_selection","dEta", NbinsJetEta, 0., 10.);
    VBF_Mjj_presel_sel = new TH1F("VBF_Mjj_presel_selection","Mjj", NbinsHT, HTmin, HTmax);
    VBF_Jet1Pt_presel_sel = new TH1F("VBF_Jet1Pt_presel_selection","Jet1Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax);
    VBF_Jet2Pt_presel_sel = new TH1F("VBF_Jet2Pt_presel_selection","Jet2Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax);
    VBF_Jet3Pt_presel_sel = new TH1F("VBF_Jet3Pt_presel_selection","Jet3Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax/2.);
    VBF_Jet1Eta_presel_sel = new TH1F("VBF_Jet1Eta_presel_selection","Jet1Eta", NbinsJetEta, JetEtamin, JetEtamax);
    VBF_Jet2Eta_presel_sel = new TH1F("VBF_Jet2Eta_presel_selection","Jet2Eta", NbinsJetEta, JetEtamin, JetEtamax);
    VBF_Jet3Eta_presel_sel = new TH1F("VBF_Jet3Eta_presel_selection","Jet3Eta", NbinsJetEta, JetEtamin, JetEtamax);
    VBF_PTjj_presel_sel = new TH1F("VBF_PTjj_presel_selection","PTjj", NbinsMHT, MHTmin, MHTmax);
    VBF_MET_presel_sel = new TH1F("VBF_MET_presel_selection","MET", NbinsMHT, MHTmin, MHTmax);
    VBF_METsoft_presel_sel = new TH1F("VBF_METsoft_presel_selection","METsoft", NbinsMHTsoft, MHTsoftmin, MHTsoftmax);
    VBF_METsig_presel_sel = new TH1F("VBF_METsig_presel_selection","METsig", NbinsMHTsig, MHTsigmin, MHTsigmax);
    VBF_MHTsig_presel_sel = new TH1F("VBF_MHTsig_presel_selection","MHTsig", NbinsMHTsig, MHTsigmin, MHTsigmax);
    VBF_minDeltaPhiPTj12_presel_sel = new TH1F("VBF_minDeltaPhiPTj12_presel_selection","minDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi());
    VBF_maxDeltaPhiPTj12_presel_sel = new TH1F("VBF_maxDeltaPhiPTj12_presel_selection","maxDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi());
    VBF_DeltaPhiPTj3_presel_sel = new TH1F("VBF_DeltaPhiPTj3_presel_selection","DeltaPhiPTj3", NbinsJetPhi, 0., TMath::Pi());

    VBF_dPhi_presel_4JV_dPhiSide_sel = new TH1F("VBF_dPhi_presel_4JV_dPhiSide_selection","dPhi", NbinsJetPhi, 0., TMath::Pi());
    VBF_dEta_presel_4JV_dPhiSide_sel = new TH1F("VBF_dEta_presel_4JV_dPhiSide_selection","dEta", NbinsJetEta, 0., 10.);
    VBF_Mjj_presel_4JV_dPhiSide_sel = new TH1F("VBF_Mjj_presel_4JV_dPhiSide_selection","Mjj", NbinsHT, HTmin, HTmax);
    VBF_Jet1Pt_presel_4JV_dPhiSide_sel = new TH1F("VBF_Jet1Pt_presel_4JV_dPhiSide_selection","Jet1Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax);
    VBF_Jet2Pt_presel_4JV_dPhiSide_sel = new TH1F("VBF_Jet2Pt_presel_4JV_dPhiSide_selection","Jet2Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax);
    VBF_Jet3Pt_presel_4JV_dPhiSide_sel = new TH1F("VBF_Jet3Pt_presel_4JV_dPhiSide_selection","Jet3Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax/2.);
    VBF_Jet1Eta_presel_4JV_dPhiSide_sel = new TH1F("VBF_Jet1Eta_presel_4JV_dPhiSide_selection","Jet1Eta", NbinsJetEta, JetEtamin, JetEtamax);
    VBF_Jet2Eta_presel_4JV_dPhiSide_sel = new TH1F("VBF_Jet2Eta_presel_4JV_dPhiSide_selection","Jet2Eta", NbinsJetEta, JetEtamin, JetEtamax);
    VBF_Jet3Eta_presel_4JV_dPhiSide_sel = new TH1F("VBF_Jet3Eta_presel_4JV_dPhiSide_selection","Jet3Eta", NbinsJetEta, JetEtamin, JetEtamax);
    VBF_PTjj_presel_4JV_dPhiSide_sel = new TH1F("VBF_PTjj_presel_4JV_dPhiSide_selection","PTjj", NbinsMHT, MHTmin, MHTmax);
    VBF_MET_presel_4JV_dPhiSide_sel = new TH1F("VBF_MET_presel_4JV_dPhiSide_selection","MET", NbinsMHT, MHTmin, MHTmax);
    VBF_METsoft_presel_4JV_dPhiSide_sel = new TH1F("VBF_METsoft_presel_4JV_dPhiSide_selection","METsoft", NbinsMHTsoft, MHTsoftmin, MHTsoftmax);
    VBF_METsig_presel_4JV_dPhiSide_sel = new TH1F("VBF_METsig_presel_4JV_dPhiSide_selection","METsig", NbinsMHTsig, MHTsigmin, MHTsigmax);
    VBF_MHTsig_presel_4JV_dPhiSide_sel = new TH1F("VBF_MHTsig_presel_4JV_dPhiSide_selection","MHTsig", NbinsMHTsig, MHTsigmin, MHTsigmax);
    VBF_minDeltaPhiPTj12_presel_4JV_dPhiSide_sel = new TH1F("VBF_minDeltaPhiPTj12_presel_4JV_dPhiSide_selection","minDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi());
    VBF_maxDeltaPhiPTj12_presel_4JV_dPhiSide_sel = new TH1F("VBF_maxDeltaPhiPTj12_presel_4JV_dPhiSide_selection","maxDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi());
    VBF_DeltaPhiPTj3_presel_4JV_dPhiSide_sel = new TH1F("VBF_DeltaPhiPTj3_presel_4JV_dPhiSide_selection","DeltaPhiPTj3", NbinsJetPhi, 0., TMath::Pi());

    VBF_dPhi_dEta_sel = new TH1F("VBF_dPhi_dEta_selection","dPhi", NbinsJetPhi, 0., TMath::Pi());
    VBF_dEta_dEta_sel = new TH1F("VBF_dEta_dEta_selection","dEta", NbinsJetEta, 0., 10.);
    VBF_Mjj_dEta_sel = new TH1F("VBF_Mjj_dEta_selection","Mjj", NbinsHT, HTmin, HTmax);
    VBF_Jet1Pt_dEta_sel = new TH1F("VBF_Jet1Pt_dEta_selection","Jet1Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax);
    VBF_Jet2Pt_dEta_sel = new TH1F("VBF_Jet2Pt_dEta_selection","Jet2Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax);
    VBF_Jet3Pt_dEta_sel = new TH1F("VBF_Jet3Pt_dEta_selection","Jet3Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax/2.);
    VBF_Jet1Eta_dEta_sel = new TH1F("VBF_Jet1Eta_dEta_selection","Jet1Eta", NbinsJetEta, JetEtamin, JetEtamax);
    VBF_Jet2Eta_dEta_sel = new TH1F("VBF_Jet2Eta_dEta_selection","Jet2Eta", NbinsJetEta, JetEtamin, JetEtamax);
    VBF_Jet3Eta_dEta_sel = new TH1F("VBF_Jet3Eta_dEta_selection","Jet3Eta", NbinsJetEta, JetEtamin, JetEtamax);
    VBF_PTjj_dEta_sel = new TH1F("VBF_PTjj_dEta_selection","PTjj", NbinsMHT, MHTmin, MHTmax);
    VBF_MET_dEta_sel = new TH1F("VBF_MET_dEta_selection","MET", NbinsMHT, MHTmin, MHTmax);
    VBF_METsoft_dEta_sel = new TH1F("VBF_METsoft_dEta_selection","METsoft", NbinsMHTsoft, MHTsoftmin, MHTsoftmax);
    VBF_METsig_dEta_sel = new TH1F("VBF_METsig_dEta_selection","METsig", NbinsMHTsig, MHTsigmin, MHTsigmax);
    VBF_MHTsig_dEta_sel = new TH1F("VBF_MHTsig_dEta_selection","MHTsig", NbinsMHTsig, MHTsigmin, MHTsigmax);
    VBF_minDeltaPhiPTj12_dEta_sel = new TH1F("VBF_minDeltaPhiPTj12_dEta_selection","minDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi());
    VBF_maxDeltaPhiPTj12_dEta_sel = new TH1F("VBF_maxDeltaPhiPTj12_dEta_selection","maxDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi());
    VBF_DeltaPhiPTj3_dEta_sel = new TH1F("VBF_DeltaPhiPTj3_dEta_selection","DeltaPhiPTj3", NbinsJetPhi, 0., TMath::Pi());

    VBF_dPhi_jj_sel = new TH1F("VBF_dPhi_jj_selection","dPhi", NbinsJetPhi, 0., TMath::Pi());
    VBF_dEta_jj_sel = new TH1F("VBF_dEta_jj_selection","dEta", NbinsJetEta, 0., 10.);
    VBF_Mjj_jj_sel = new TH1F("VBF_Mjj_jj_selection","Mjj", NbinsHT, HTmin, HTmax);
    VBF_Jet1Pt_jj_sel = new TH1F("VBF_Jet1Pt_jj_selection","Jet1Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax);
    VBF_Jet2Pt_jj_sel = new TH1F("VBF_Jet2Pt_jj_selection","Jet2Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax);
    VBF_Jet3Pt_jj_sel = new TH1F("VBF_Jet3Pt_jj_selection","Jet3Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax/2.);
    VBF_Jet1Eta_jj_sel = new TH1F("VBF_Jet1Eta_jj_selection","Jet1Eta", NbinsJetEta, JetEtamin, JetEtamax);
    VBF_Jet2Eta_jj_sel = new TH1F("VBF_Jet2Eta_jj_selection","Jet2Eta", NbinsJetEta, JetEtamin, JetEtamax);
    VBF_Jet3Eta_jj_sel = new TH1F("VBF_Jet3Eta_jj_selection","Jet3Eta", NbinsJetEta, JetEtamin, JetEtamax);
    VBF_PTjj_jj_sel = new TH1F("VBF_PTjj_jj_selection","PTjj", NbinsMHT, MHTmin, MHTmax);
    VBF_MET_jj_sel = new TH1F("VBF_MET_jj_selection","MET", NbinsMHT, MHTmin, MHTmax);
    VBF_METsoft_jj_sel = new TH1F("VBF_METsoft_jj_selection","METsoft", NbinsMHTsoft, MHTsoftmin, MHTsoftmax);
    VBF_METsig_jj_sel = new TH1F("VBF_METsig_jj_selection","METsig", NbinsMHTsig, MHTsigmin, MHTsigmax);
    VBF_MHTsig_jj_sel = new TH1F("VBF_MHTsig_jj_selection","MHTsig", NbinsMHTsig, MHTsigmin, MHTsigmax);
    VBF_minDeltaPhiPTj12_jj_sel = new TH1F("VBF_minDeltaPhiPTj12_jj_selection","minDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi());
    VBF_maxDeltaPhiPTj12_jj_sel = new TH1F("VBF_maxDeltaPhiPTj12_jj_selection","maxDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi());
    VBF_DeltaPhiPTj3_jj_sel = new TH1F("VBF_DeltaPhiPTj3_jj_selection","DeltaPhiPTj3", NbinsJetPhi, 0., TMath::Pi());

    VBF_dPhi_dEta_3JV_sel = new TH1F("VBF_dPhi_dEta_3JV_selection","dPhi", NbinsJetPhi, 0., TMath::Pi());
    VBF_dEta_dEta_3JV_sel = new TH1F("VBF_dEta_dEta_3JV_selection","dEta", NbinsJetEta, 0., 10.);
    VBF_Mjj_dEta_3JV_sel = new TH1F("VBF_Mjj_dEta_3JV_selection","Mjj", NbinsHT, HTmin, HTmax);
    VBF_Jet1Pt_dEta_3JV_sel = new TH1F("VBF_Jet1Pt_dEta_3JV_selection","Jet1Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax);
    VBF_Jet2Pt_dEta_3JV_sel = new TH1F("VBF_Jet2Pt_dEta_3JV_selection","Jet2Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax);
    VBF_Jet3Pt_dEta_3JV_sel = new TH1F("VBF_Jet3Pt_dEta_3JV_selection","Jet3Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax/2.);
    VBF_Jet1Eta_dEta_3JV_sel = new TH1F("VBF_Jet1Eta_dEta_3JV_selection","Jet1Eta", NbinsJetEta, JetEtamin, JetEtamax);
    VBF_Jet2Eta_dEta_3JV_sel = new TH1F("VBF_Jet2Eta_dEta_3JV_selection","Jet2Eta", NbinsJetEta, JetEtamin, JetEtamax);
    VBF_Jet3Eta_dEta_3JV_sel = new TH1F("VBF_Jet3Eta_dEta_3JV_selection","Jet3Eta", NbinsJetEta, JetEtamin, JetEtamax);
    VBF_PTjj_dEta_3JV_sel = new TH1F("VBF_PTjj_dEta_3JV_selection","PTjj", NbinsMHT, MHTmin, MHTmax);
    VBF_MET_dEta_3JV_sel = new TH1F("VBF_MET_dEta_3JV_selection","MET", NbinsMHT, MHTmin, MHTmax);
    VBF_METsoft_dEta_3JV_sel = new TH1F("VBF_METsoft_dEta_3JV_selection","METsoft", NbinsMHTsoft, MHTsoftmin, MHTsoftmax);
    VBF_METsig_dEta_3JV_sel = new TH1F("VBF_METsig_dEta_3JV_selection","METsig", NbinsMHTsig, MHTsigmin, MHTsigmax);
    VBF_MHTsig_dEta_3JV_sel = new TH1F("VBF_MHTsig_dEta_3JV_selection","MHTsig", NbinsMHTsig, MHTsigmin, MHTsigmax);
    VBF_minDeltaPhiPTj12_dEta_3JV_sel = new TH1F("VBF_minDeltaPhiPTj12_dEta_3JV_selection","minDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi());
    VBF_maxDeltaPhiPTj12_dEta_3JV_sel = new TH1F("VBF_maxDeltaPhiPTj12_dEta_3JV_selection","maxDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi());
    VBF_DeltaPhiPTj3_dEta_3JV_sel = new TH1F("VBF_DeltaPhiPTj3_dEta_3JV_selection","DeltaPhiPTj3", NbinsJetPhi, 0., TMath::Pi());

    VBF_dPhi_dEta_3JV_dPhiPTjj_sel = new TH1F("VBF_dPhi_dEta_3JV_dPhiPTjj_selection","dPhi", NbinsJetPhi, 0., TMath::Pi());
    VBF_dEta_dEta_3JV_dPhiPTjj_sel = new TH1F("VBF_dEta_dEta_3JV_dPhiPTjj_selection","dEta", NbinsJetEta, 0., 10.);
    VBF_Mjj_dEta_3JV_dPhiPTjj_sel = new TH1F("VBF_Mjj_dEta_3JV_dPhiPTjj_selection","Mjj", NbinsHT, HTmin, HTmax);
    VBF_Jet1Pt_dEta_3JV_dPhiPTjj_sel = new TH1F("VBF_Jet1Pt_dEta_3JV_dPhiPTjj_selection","Jet1Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax);
    VBF_Jet2Pt_dEta_3JV_dPhiPTjj_sel = new TH1F("VBF_Jet2Pt_dEta_3JV_dPhiPTjj_selection","Jet2Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax);
    VBF_Jet3Pt_dEta_3JV_dPhiPTjj_sel = new TH1F("VBF_Jet3Pt_dEta_3JV_dPhiPTjj_selection","Jet3Pt", NbinsJetPt, JetVBFPtmin, JetVBFPtmax/2.);
    VBF_Jet1Eta_dEta_3JV_dPhiPTjj_sel = new TH1F("VBF_Jet1Eta_dEta_3JV_dPhiPTjj_selection","Jet1Eta", NbinsJetEta, JetEtamin, JetEtamax);
    VBF_Jet2Eta_dEta_3JV_dPhiPTjj_sel = new TH1F("VBF_Jet2Eta_dEta_3JV_dPhiPTjj_selection","Jet2Eta", NbinsJetEta, JetEtamin, JetEtamax);
    VBF_Jet3Eta_dEta_3JV_dPhiPTjj_sel = new TH1F("VBF_Jet3Eta_dEta_3JV_dPhiPTjj_selection","Jet3Eta", NbinsJetEta, JetEtamin, JetEtamax);
    VBF_PTjj_dEta_3JV_dPhiPTjj_sel = new TH1F("VBF_PTjj_dEta_3JV_dPhiPTjj_selection","PTjj", NbinsMHT, MHTmin, MHTmax);
    VBF_MET_dEta_3JV_dPhiPTjj_sel = new TH1F("VBF_MET_dEta_3JV_dPhiPTjj_selection","MET", NbinsMHT, MHTmin, MHTmax);
    VBF_METsoft_dEta_3JV_dPhiPTjj_sel = new TH1F("VBF_METsoft_dEta_3JV_dPhiPTjj_selection","METsoft", NbinsMHTsoft, MHTsoftmin, MHTsoftmax);
    VBF_METsig_dEta_3JV_dPhiPTjj_sel = new TH1F("VBF_METsig_dEta_3JV_dPhiPTjj_selection","METsig", NbinsMHTsig, MHTsigmin, MHTsigmax);
    VBF_MHTsig_dEta_3JV_dPhiPTjj_sel = new TH1F("VBF_MHTsig_dEta_3JV_dPhiPTjj_selection","MHTsig", NbinsMHTsig, MHTsigmin, MHTsigmax);
    VBF_minDeltaPhiPTj12_dEta_3JV_dPhiPTjj_sel = new TH1F("VBF_minDeltaPhiPTj12_dEta_3JV_dPhiPTjj_selection","minDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi());
    VBF_maxDeltaPhiPTj12_dEta_3JV_dPhiPTjj_sel = new TH1F("VBF_maxDeltaPhiPTj12_dEta_3JV_dPhiPTjj_selection","maxDeltaPhiPTj12", NbinsJetPhi, 0., TMath::Pi());
    VBF_DeltaPhiPTj3_dEta_3JV_dPhiPTjj_sel = new TH1F("VBF_DeltaPhiPTj3_dEta_3JV_dPhiPTjj_selection","DeltaPhiPTj3", NbinsJetPhi, 0., TMath::Pi());

    // NJets
    NJets_baseline_withoutMET_sel = new TH1F("NJets_baseline_withoutMET_sel", "NJets baseline", 15, 0, 15);
    NJets_baseline_sel = new TH1F("NJets_baseline_sel", "NJets baseline", 15, 0, 15);
    NJets_baseline_withoutDeltaPhi_withoutMET_sel = new TH1F("NJets_baseline_withoutDeltaPhi_withoutMET_sel", "NJets baseline", 15, 0, 15);
    NJets_baseline_withoutDeltaPhi_sel = new TH1F("NJets_baseline_withoutDeltaPhi_sel", "NJets baseline", 15, 0, 15);
    NJets_presel_sel = new TH1F("NJets_presel_sel", "NJets presel", 15, 0, 15);

    // NBJets
    NBJets_baseline_withoutMET_sel = new TH1F("NBJets_baseline_withoutMET_sel", "NBJets baseline", 5, 0, 5);
    NBJets_baseline_sel = new TH1F("NBJets_baseline_sel", "NBJets baseline", 5, 0, 5);
    NBJets_baseline_withoutDeltaPhi_withoutMET_sel = new TH1F("NBJets_baseline_withoutDeltaPhi_withoutMET_sel", "NBJets baseline", 5, 0, 5);
    NBJets_baseline_withoutDeltaPhi_sel = new TH1F("NBJets_baseline_withoutDeltaPhi_sel", "NBJets baseline", 5, 0, 5);
    NBJets_presel_sel = new TH1F("NBJets_presel_sel", "NBJets presel", 5, 0, 5);

    // ------------------------------------------------------------------------------ //

    // get tree with predictions
    cout << "entries prediction tree:" << QCDPrediction.GetEntries() << endl;

    JetPt = new std::vector<Float_t>;
    JetEta = new std::vector<Float_t>;
    JetPhi = new std::vector<Float_t>;
    JetM = new std::vector<Float_t>;
    DeltaPhi = new std::vector<Float_t>;
    QCDPrediction.SetBranchAddress("Ntries",&Ntries);
    QCDPrediction.SetBranchAddress("NJets",&NJets);
    QCDPrediction.SetBranchAddress("BTags",&BTags);
    QCDPrediction.SetBranchAddress("Weight",&weight0);
    QCDPrediction.SetBranchAddress("TriggerWeight",&triggerWeight);
    QCDPrediction.SetBranchAddress("HT",&HT);
    QCDPrediction.SetBranchAddress("MHT",&MHT);
    QCDPrediction.SetBranchAddress("MHTphi",&MHTphi);
    QCDPrediction.SetBranchAddress("MET",&MET);
    QCDPrediction.SetBranchAddress("METphi",&METphi);
    QCDPrediction.SetBranchAddress("METsig",&METsig);
    QCDPrediction.SetBranchAddress("MHTsig",&MHTsig);
    QCDPrediction.SetBranchAddress("METsoft",&METsoft);
    QCDPrediction.SetBranchAddress("JetPt",&JetPt);
    QCDPrediction.SetBranchAddress("JetEta",&JetEta);
    QCDPrediction.SetBranchAddress("JetPhi",&JetPhi);
    QCDPrediction.SetBranchAddress("JetM",&JetM);
    QCDPrediction.SetBranchAddress("DeltaPhi",&DeltaPhi);

    cout << "Needed trees are loaded!" << endl;

    // loop over entries and fill prediction histos
    ULong_t nentries = QCDPrediction.GetEntries();

    for ( ULong_t i = 0 ; i<nentries ; i++) {

        if (i == 0 ) cout << "load first event!" << endl;
        QCDPrediction.GetEntry(i);
        if (i == 0 ) cout << "loaded first event successfully!" << endl;

		float tw = triggerWeight;
		tw = 1.; // for CRs
		//if (Ntries > 0 && triggerWeight < 5.e-2) tw = 0.;
        float weight = weight0 * tw;

        if( i%100000 == 0 ) std::cout << "event (prediction): " << i << " (" << 100*double(i)/nentries << "%)"<< '\n';

        double Mjj = CalcMjj();
        double MHTjj = CalcMHTjj();
        double DPhi = CalcDeltaPhi();
        double DEta = CalcDeltaEta();

        // apply some baseline cuts as used for storing in result tree
        if( HT > HTSave && ( MHTjj > MHTjjSave || MHT > MHTSave || MET > METSave ) && NJets >= NJetsSave && Mjj > MjjSave && DPhi < dPhijjSave && DEta > dEtajjSave) {

            //std::cout << "Ntries, weight, MC weight, triggerWeight: " << Ntries << ", " << weight << ", " << weight0 << ", " << triggerWeight << std::endl;
            //std::cout << "pTjj, HT: " << MHTjj <<  ", " << HT << std::endl;
            //std::cout << "MET (pt, phi): " << MET <<  ", " << METphi << std::endl;
            //std::cout << "DPhi, DEta: " << DPhi <<  ", " << DEta << std::endl;
            //std::cout << "Mjj: " << Mjj << std::endl;
            //std::cout << "1st jet (pt, eta, phi): " << JetPt->at(0) << ", " << JetEta->at(0)  << ", " << JetPhi->at(0) << std::endl;
            //std::cout << "2nd jet (pt, eta, phi): " << JetPt->at(1) << ", " << JetEta->at(1)  << ", " << JetPhi->at(1) << std::endl;
            //std::cout << "3rd jet (pt, eta, phi): " << JetPt->at(2) << ", " << JetEta->at(2)  << ", " << JetPhi->at(2) << std::endl;
            //std::cout << std::endl;

            double DPhiMET1, DPhiMET2, DPhiMET3;
            TLorentzVector vMET(0., 0., 0., 0.);
            vMET.SetPtEtaPhiM(MET, 0., METphi, 0.);
            CalcDPhiMET(DPhiMET1, DPhiMET2, DPhiMET3, vMET);

            if (VBF) {

                if (MjjJetSel() && Mjj > MjjCut && MET > METCut && DEta > DEtaLoose && DPhi < DPhiCRMax && DPhi > DPhiCRMin) {

                    if ( Soft3rd() && Veto4th() ) {

                        if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                            VBF_dPhi_presel_4JV_dPhiSide_pred_raw->Fill(DPhi, Ntries, weight);
                            VBF_dEta_presel_4JV_dPhiSide_pred_raw->Fill(DEta, Ntries, weight);
                            VBF_Mjj_presel_4JV_dPhiSide_pred_raw->Fill(Mjj, Ntries, weight);
                            VBF_Jet1Pt_presel_4JV_dPhiSide_pred_raw->Fill(JetPt->at(0), Ntries, weight);
                            VBF_Jet2Pt_presel_4JV_dPhiSide_pred_raw->Fill(JetPt->at(1), Ntries, weight);
                            VBF_Jet3Pt_presel_4JV_dPhiSide_pred_raw->Fill(JetPt->at(2), Ntries, weight);
                            VBF_Jet1Eta_presel_4JV_dPhiSide_pred_raw->Fill(JetEta->at(0), Ntries, weight);
                            VBF_Jet2Eta_presel_4JV_dPhiSide_pred_raw->Fill(JetEta->at(1), Ntries, weight);
                            VBF_Jet3Eta_presel_4JV_dPhiSide_pred_raw->Fill(JetEta->at(2), Ntries, weight);
                            VBF_PTjj_presel_4JV_dPhiSide_pred_raw->Fill(MHTjj, Ntries, weight);
                            VBF_MET_presel_4JV_dPhiSide_pred_raw->Fill(MET, Ntries, weight);
                            VBF_METsoft_presel_4JV_dPhiSide_pred_raw->Fill(METsoft, Ntries, weight);
                            VBF_METsig_presel_4JV_dPhiSide_pred_raw->Fill(METsig, Ntries, weight);
                            VBF_MHTsig_presel_4JV_dPhiSide_pred_raw->Fill(MHTsig, Ntries, weight);
                            VBF_minDeltaPhiPTj12_presel_4JV_dPhiSide_pred_raw->Fill(DPhiMET1, Ntries, weight);
                            VBF_maxDeltaPhiPTj12_presel_4JV_dPhiSide_pred_raw->Fill(DPhiMET2, Ntries, weight);
                            VBF_DeltaPhiPTj3_presel_4JV_dPhiSide_pred_raw->Fill(DPhiMET3, Ntries, weight);
                        }  else if (Ntries == -2) {
                            VBF_dPhi_presel_4JV_dPhiSide_sel->Fill(DPhi, weight);
                            VBF_dEta_presel_4JV_dPhiSide_sel->Fill(DEta, weight);
                            VBF_Mjj_presel_4JV_dPhiSide_sel->Fill(Mjj, weight);
                            VBF_Jet1Pt_presel_4JV_dPhiSide_sel->Fill(JetPt->at(0), weight);
                            VBF_Jet2Pt_presel_4JV_dPhiSide_sel->Fill(JetPt->at(1), weight);
                            VBF_Jet3Pt_presel_4JV_dPhiSide_sel->Fill(JetPt->at(2), weight);
                            VBF_Jet1Eta_presel_4JV_dPhiSide_sel->Fill(JetEta->at(0), weight);
                            VBF_Jet2Eta_presel_4JV_dPhiSide_sel->Fill(JetEta->at(1), weight);
                            VBF_Jet3Eta_presel_4JV_dPhiSide_sel->Fill(JetEta->at(2), weight);
                            VBF_PTjj_presel_4JV_dPhiSide_sel->Fill(MHTjj, weight);
                            VBF_MET_presel_4JV_dPhiSide_sel->Fill(MET, weight);
                            VBF_METsoft_presel_4JV_dPhiSide_sel->Fill(METsoft, weight);
							VBF_METsig_presel_4JV_dPhiSide_sel->Fill(METsig, weight);
                            VBF_MHTsig_presel_4JV_dPhiSide_sel->Fill(MHTsig, weight);
							VBF_minDeltaPhiPTj12_presel_4JV_dPhiSide_sel->Fill(DPhiMET1, weight);
                            VBF_maxDeltaPhiPTj12_presel_4JV_dPhiSide_sel->Fill(DPhiMET2, weight);
                            VBF_DeltaPhiPTj3_presel_4JV_dPhiSide_sel->Fill(DPhiMET3, weight);
                        }

                    }

                }

                if (MjjJetSel() && MET > METCut && DPhi < DPhiCRMax) {

                    if ( ((Soft3rd() && Veto4th()) || (Veto3rd() && Veto4th())) && DEta > DEtajj && Mjj > MjjCut) {

                        if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                            VBF_dPhi_presel_pred_raw->Fill(DPhi, Ntries, weight);
                            VBF_dEta_presel_pred_raw->Fill(DEta, Ntries, weight);
                            VBF_Mjj_presel_pred_raw->Fill(Mjj, Ntries, weight);
                            VBF_Jet1Pt_presel_pred_raw->Fill(JetPt->at(0), Ntries, weight);
                            VBF_Jet2Pt_presel_pred_raw->Fill(JetPt->at(1), Ntries, weight);
                            VBF_Jet3Pt_presel_pred_raw->Fill(JetPt->at(2), Ntries, weight);
                            VBF_Jet1Eta_presel_pred_raw->Fill(JetEta->at(0), Ntries, weight);
                            VBF_Jet2Eta_presel_pred_raw->Fill(JetEta->at(1), Ntries, weight);
                            VBF_Jet3Eta_presel_pred_raw->Fill(JetEta->at(2), Ntries, weight);
                            VBF_PTjj_presel_pred_raw->Fill(MHTjj, Ntries, weight);
                            VBF_MET_presel_pred_raw->Fill(MET, Ntries, weight);
                            VBF_METsoft_presel_pred_raw->Fill(METsoft, Ntries, weight);
                            VBF_METsig_presel_pred_raw->Fill(METsig, Ntries, weight);
                            VBF_MHTsig_presel_pred_raw->Fill(MHTsig, Ntries, weight);
                            VBF_minDeltaPhiPTj12_presel_pred_raw->Fill(DPhiMET1, Ntries, weight);
                            VBF_maxDeltaPhiPTj12_presel_pred_raw->Fill(DPhiMET2, Ntries, weight);
                            VBF_DeltaPhiPTj3_presel_pred_raw->Fill(DPhiMET3, Ntries, weight);
                        }  else if (Ntries == -2) {
                            VBF_dPhi_presel_sel->Fill(DPhi, weight);
                            VBF_dEta_presel_sel->Fill(DEta, weight);
                            VBF_Mjj_presel_sel->Fill(Mjj, weight);
                            VBF_Jet1Pt_presel_sel->Fill(JetPt->at(0), weight);
                            VBF_Jet2Pt_presel_sel->Fill(JetPt->at(1), weight);
                            VBF_Jet3Pt_presel_sel->Fill(JetPt->at(2), weight);
                            VBF_Jet1Eta_presel_sel->Fill(JetEta->at(0), weight);
                            VBF_Jet2Eta_presel_sel->Fill(JetEta->at(1), weight);
                            VBF_Jet3Eta_presel_sel->Fill(JetEta->at(2), weight);
                            VBF_PTjj_presel_sel->Fill(MHTjj, weight);
                            VBF_MET_presel_sel->Fill(MET, weight);
                            VBF_METsoft_presel_sel->Fill(METsoft, weight);
                            VBF_METsig_presel_sel->Fill(METsig, weight);
                            VBF_MHTsig_presel_sel->Fill(MHTsig, weight);
                            VBF_minDeltaPhiPTj12_presel_sel->Fill(DPhiMET1, weight);
                            VBF_maxDeltaPhiPTj12_presel_sel->Fill(DPhiMET2, weight);
                            VBF_DeltaPhiPTj3_presel_sel->Fill(DPhiMET3, weight);
                        }

                        if (DEta > DEtaTight && DPhi < DPhiSR) {

                            if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                                VBF_dPhi_dEta_pred_raw->Fill(DPhi, Ntries, weight);
                                VBF_dEta_dEta_pred_raw->Fill(DEta, Ntries, weight);
                                VBF_Mjj_dEta_pred_raw->Fill(Mjj, Ntries, weight);
                                VBF_Jet1Pt_dEta_pred_raw->Fill(JetPt->at(0), Ntries, weight);
                                VBF_Jet2Pt_dEta_pred_raw->Fill(JetPt->at(1), Ntries, weight);
                                VBF_Jet3Pt_dEta_pred_raw->Fill(JetPt->at(2), Ntries, weight);
                                VBF_Jet1Eta_dEta_pred_raw->Fill(JetEta->at(0), Ntries, weight);
                                VBF_Jet2Eta_dEta_pred_raw->Fill(JetEta->at(1), Ntries, weight);
                                VBF_Jet3Eta_dEta_pred_raw->Fill(JetEta->at(2), Ntries, weight);
                                VBF_PTjj_dEta_pred_raw->Fill(MHTjj, Ntries, weight);
                                VBF_MET_dEta_pred_raw->Fill(MET, Ntries, weight);
                                VBF_METsoft_dEta_pred_raw->Fill(METsoft, Ntries, weight);
                                VBF_METsig_dEta_pred_raw->Fill(METsig, Ntries, weight);
                                VBF_MHTsig_dEta_pred_raw->Fill(MHTsig, Ntries, weight);
                                VBF_minDeltaPhiPTj12_dEta_pred_raw->Fill(DPhiMET1, Ntries, weight);
                                VBF_maxDeltaPhiPTj12_dEta_pred_raw->Fill(DPhiMET2, Ntries, weight);
                                VBF_DeltaPhiPTj3_dEta_pred_raw->Fill(DPhiMET3, Ntries, weight);
                            }  else if (Ntries == -2) {
                                VBF_dPhi_dEta_sel->Fill(DPhi, weight);
                                VBF_dEta_dEta_sel->Fill(DEta, weight);
                                VBF_Mjj_dEta_sel->Fill(Mjj, weight);
                                VBF_Jet1Pt_dEta_sel->Fill(JetPt->at(0), weight);
                                VBF_Jet2Pt_dEta_sel->Fill(JetPt->at(1), weight);
                                VBF_Jet3Pt_dEta_sel->Fill(JetPt->at(2), weight);
                                VBF_Jet1Eta_dEta_sel->Fill(JetEta->at(0), weight);
                                VBF_Jet2Eta_dEta_sel->Fill(JetEta->at(1), weight);
                                VBF_Jet3Eta_dEta_sel->Fill(JetEta->at(2), weight);
                                VBF_PTjj_dEta_sel->Fill(MHTjj, weight);
                                VBF_MET_dEta_sel->Fill(MET, weight);
                                VBF_METsoft_dEta_sel->Fill(METsoft, weight);
                                VBF_METsig_dEta_sel->Fill(METsig, weight);
                                VBF_MHTsig_dEta_sel->Fill(MHTsig, weight);
                                VBF_minDeltaPhiPTj12_dEta_sel->Fill(DPhiMET1, weight);
                                VBF_maxDeltaPhiPTj12_dEta_sel->Fill(DPhiMET2, weight);
                                VBF_DeltaPhiPTj3_dEta_sel->Fill(DPhiMET3, weight);
                            }
                        }

                    }

                    if (DPhi < DPhiSR && DEta > DEtajj && Mjj < MjjCut && Veto3rd()) { // && ( (fabs(JetEta->at(0))<2.4 && fabs(JetEta->at(1))>2.4 ) || (fabs(JetEta->at(1))<2.4 && fabs(JetEta->at(0))>2.4 ) ) ) {

                        if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                            VBF_dPhi_jj_pred_raw->Fill(DPhi, Ntries, weight);
                            VBF_dEta_jj_pred_raw->Fill(DEta, Ntries, weight);
                            VBF_Mjj_jj_pred_raw->Fill(Mjj, Ntries, weight);
                            VBF_Jet1Pt_jj_pred_raw->Fill(JetPt->at(0), Ntries, weight);
                            VBF_Jet2Pt_jj_pred_raw->Fill(JetPt->at(1), Ntries, weight);
                            VBF_Jet3Pt_jj_pred_raw->Fill(JetPt->at(2), Ntries, weight);
                            VBF_Jet1Eta_jj_pred_raw->Fill(JetEta->at(0), Ntries, weight);
                            VBF_Jet2Eta_jj_pred_raw->Fill(JetEta->at(1), Ntries, weight);
                            VBF_Jet3Eta_jj_pred_raw->Fill(JetEta->at(2), Ntries, weight);
                            VBF_PTjj_jj_pred_raw->Fill(MHTjj, Ntries, weight);
                            VBF_MET_jj_pred_raw->Fill(MET, Ntries, weight);
                            VBF_METsoft_jj_pred_raw->Fill(METsoft, Ntries, weight);
                            VBF_METsig_jj_pred_raw->Fill(METsig, Ntries, weight);
                            VBF_MHTsig_jj_pred_raw->Fill(MHTsig, Ntries, weight);
                            VBF_minDeltaPhiPTj12_jj_pred_raw->Fill(DPhiMET1, Ntries, weight);
                            VBF_maxDeltaPhiPTj12_jj_pred_raw->Fill(DPhiMET2, Ntries, weight);
                            VBF_DeltaPhiPTj3_jj_pred_raw->Fill(DPhiMET3, Ntries, weight);
                        }  else if (Ntries == -2) {
                            VBF_dPhi_jj_sel->Fill(DPhi, weight);
                            VBF_dEta_jj_sel->Fill(DEta, weight);
                            VBF_Mjj_jj_sel->Fill(Mjj, weight);
                            VBF_Jet1Pt_jj_sel->Fill(JetPt->at(0), weight);
                            VBF_Jet2Pt_jj_sel->Fill(JetPt->at(1), weight);
                            VBF_Jet3Pt_jj_sel->Fill(JetPt->at(2), weight);
                            VBF_Jet1Eta_jj_sel->Fill(JetEta->at(0), weight);
                            VBF_Jet2Eta_jj_sel->Fill(JetEta->at(1), weight);
                            VBF_Jet3Eta_jj_sel->Fill(JetEta->at(2), weight);
                            VBF_PTjj_jj_sel->Fill(MHTjj, weight);
                            VBF_MET_jj_sel->Fill(MET, weight);
                            VBF_METsoft_jj_sel->Fill(METsoft, weight);
                            VBF_METsig_jj_sel->Fill(METsig, weight);
                            VBF_MHTsig_jj_sel->Fill(MHTsig, weight);
                            VBF_minDeltaPhiPTj12_jj_sel->Fill(DPhiMET1,weight);
                            VBF_maxDeltaPhiPTj12_jj_sel->Fill(DPhiMET2, weight);
                            VBF_DeltaPhiPTj3_jj_sel->Fill(DPhiMET3, weight);
                        }

                    }

                    if (DEta > DEtaTight && Veto3rd() && Mjj > MjjCut && DPhi < DPhiSR) {

                        //std::cout << "Ntries, weight, MC weight, triggerWeight: " << Ntries << ", " << weight << ", " << weight0 << ", " << triggerWeight << std::endl;
                        //std::cout << "pTjj, HT: " << MHTjj <<  ", " << HT << std::endl;
                        //std::cout << "MET (pt, phi): " << MET <<  ", " << METphi << std::endl;
                        //std::cout << "1st jet (pt, eta, phi): " << JetPt->at(0) << ", " << JetEta->at(0)  << ", " << JetPhi->at(0) << std::endl;
                        //std::cout << "2nd jet (pt, eta, phi): " << JetPt->at(1) << ", " << JetEta->at(1)  << ", " << JetPhi->at(1) << std::endl;
                        //std::cout << "3rd jet (pt, eta, phi): " << JetPt->at(2) << ", " << JetEta->at(2)  << ", " << JetPhi->at(2) << std::endl;
                        //std::cout << std::endl;

                        if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                            VBF_dPhi_dEta_3JV_pred_raw->Fill(DPhi, Ntries, weight);
                            VBF_dEta_dEta_3JV_pred_raw->Fill(DEta, Ntries, weight);
                            VBF_Mjj_dEta_3JV_pred_raw->Fill(Mjj, Ntries, weight);
                            VBF_Jet1Pt_dEta_3JV_pred_raw->Fill(JetPt->at(0), Ntries, weight);
                            VBF_Jet2Pt_dEta_3JV_pred_raw->Fill(JetPt->at(1), Ntries, weight);
                            VBF_Jet3Pt_dEta_3JV_pred_raw->Fill(JetPt->at(2), Ntries, weight);
                            VBF_Jet1Eta_dEta_3JV_pred_raw->Fill(JetEta->at(0), Ntries, weight);
                            VBF_Jet2Eta_dEta_3JV_pred_raw->Fill(JetEta->at(1), Ntries, weight);
                            VBF_Jet3Eta_dEta_3JV_pred_raw->Fill(JetEta->at(2), Ntries, weight);
                            VBF_PTjj_dEta_3JV_pred_raw->Fill(MHTjj, Ntries, weight);
                            VBF_MET_dEta_3JV_pred_raw->Fill(MET, Ntries, weight);
                            VBF_METsoft_dEta_3JV_pred_raw->Fill(METsoft, Ntries, weight);
                            VBF_METsig_dEta_3JV_pred_raw->Fill(METsig, Ntries, weight);
                            VBF_MHTsig_dEta_3JV_pred_raw->Fill(MHTsig, Ntries, weight);
                            VBF_minDeltaPhiPTj12_dEta_3JV_pred_raw->Fill(DPhiMET1, Ntries, weight);
                            VBF_maxDeltaPhiPTj12_dEta_3JV_pred_raw->Fill(DPhiMET2, Ntries, weight);
                            VBF_DeltaPhiPTj3_dEta_3JV_pred_raw->Fill(DPhiMET3, Ntries, weight);
                        }  else if (Ntries == -2) {
                            VBF_dPhi_dEta_3JV_sel->Fill(DPhi, (!blindSR)*weight);
                            VBF_dEta_dEta_3JV_sel->Fill(DEta, (!blindSR)*weight);
                            VBF_Mjj_dEta_3JV_sel->Fill(Mjj, (!blindSR)*weight);
                            VBF_Jet1Pt_dEta_3JV_sel->Fill(JetPt->at(0), (!blindSR)*weight);
                            VBF_Jet2Pt_dEta_3JV_sel->Fill(JetPt->at(1), (!blindSR)*weight);
                            VBF_Jet3Pt_dEta_3JV_sel->Fill(JetPt->at(2), (!blindSR)*weight);
                            VBF_Jet1Eta_dEta_3JV_sel->Fill(JetEta->at(0), (!blindSR)*weight);
                            VBF_Jet2Eta_dEta_3JV_sel->Fill(JetEta->at(1), (!blindSR)*weight);
                            VBF_Jet3Eta_dEta_3JV_sel->Fill(JetEta->at(2), (!blindSR)*weight);
                            VBF_PTjj_dEta_3JV_sel->Fill(MHTjj, (!blindSR)*weight);
                            VBF_MET_dEta_3JV_sel->Fill(MET, (!blindSR)*weight);
                            VBF_METsoft_dEta_3JV_sel->Fill(METsoft, (!blindSR)*weight);
                            VBF_METsig_dEta_3JV_sel->Fill(METsig, (!blindSR)*weight);
                            VBF_MHTsig_dEta_3JV_sel->Fill(MHTsig, (!blindSR)*weight);
                            VBF_minDeltaPhiPTj12_dEta_3JV_sel->Fill(DPhiMET1, (!blindSR)*weight);
                            VBF_maxDeltaPhiPTj12_dEta_3JV_sel->Fill(DPhiMET2, (!blindSR)*weight);
                            VBF_DeltaPhiPTj3_dEta_3JV_sel->Fill(DPhiMET3, (!blindSR)*weight);
                        }

                        if (DPhiMET1 > 1.0 && DPhiMET2 > 1.0) {

                            if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                                VBF_dPhi_dEta_3JV_dPhiPTjj_pred_raw->Fill(DPhi, Ntries, weight);
                                VBF_dEta_dEta_3JV_dPhiPTjj_pred_raw->Fill(DEta, Ntries, weight);
                                VBF_Mjj_dEta_3JV_dPhiPTjj_pred_raw->Fill(Mjj, Ntries, weight);
                                VBF_Jet1Pt_dEta_3JV_dPhiPTjj_pred_raw->Fill(JetPt->at(0), Ntries, weight);
                                VBF_Jet2Pt_dEta_3JV_dPhiPTjj_pred_raw->Fill(JetPt->at(1), Ntries, weight);
                                VBF_Jet3Pt_dEta_3JV_dPhiPTjj_pred_raw->Fill(JetPt->at(2), Ntries, weight);
                                VBF_Jet1Eta_dEta_3JV_dPhiPTjj_pred_raw->Fill(JetEta->at(0), Ntries, weight);
                                VBF_Jet2Eta_dEta_3JV_dPhiPTjj_pred_raw->Fill(JetEta->at(1), Ntries, weight);
                                VBF_Jet3Eta_dEta_3JV_dPhiPTjj_pred_raw->Fill(JetEta->at(2), Ntries, weight);
                                VBF_PTjj_dEta_3JV_dPhiPTjj_pred_raw->Fill(MHTjj, Ntries, weight);
                                VBF_MET_dEta_3JV_dPhiPTjj_pred_raw->Fill(MET, Ntries, weight);
                                VBF_METsoft_dEta_3JV_dPhiPTjj_pred_raw->Fill(METsoft, Ntries, weight);
                                VBF_METsig_dEta_3JV_dPhiPTjj_pred_raw->Fill(METsig, Ntries, weight);
                                VBF_MHTsig_dEta_3JV_dPhiPTjj_pred_raw->Fill(MHTsig, Ntries, weight);
                                VBF_minDeltaPhiPTj12_dEta_3JV_dPhiPTjj_pred_raw->Fill(DPhiMET1, Ntries, weight);
                                VBF_maxDeltaPhiPTj12_dEta_3JV_dPhiPTjj_pred_raw->Fill(DPhiMET2, Ntries, weight);
                                VBF_DeltaPhiPTj3_dEta_3JV_dPhiPTjj_pred_raw->Fill(DPhiMET3, Ntries, weight);
                            }  else if (Ntries == -2) {
                                VBF_dPhi_dEta_3JV_dPhiPTjj_sel->Fill(DPhi, (!blindSR)*weight);
                                VBF_dEta_dEta_3JV_dPhiPTjj_sel->Fill(DEta, (!blindSR)*weight);
                                VBF_Mjj_dEta_3JV_dPhiPTjj_sel->Fill(Mjj, (!blindSR)*weight);
                                VBF_Jet1Pt_dEta_3JV_dPhiPTjj_sel->Fill(JetPt->at(0), (!blindSR)*weight);
                                VBF_Jet2Pt_dEta_3JV_dPhiPTjj_sel->Fill(JetPt->at(1), (!blindSR)*weight);
                                VBF_Jet3Pt_dEta_3JV_dPhiPTjj_sel->Fill(JetPt->at(2), (!blindSR)*weight);
                                VBF_Jet1Eta_dEta_3JV_dPhiPTjj_sel->Fill(JetEta->at(0), (!blindSR)*weight);
                                VBF_Jet2Eta_dEta_3JV_dPhiPTjj_sel->Fill(JetEta->at(1), (!blindSR)*weight);
                                VBF_Jet3Eta_dEta_3JV_dPhiPTjj_sel->Fill(JetEta->at(2), (!blindSR)*weight);
                                VBF_PTjj_dEta_3JV_dPhiPTjj_sel->Fill(MHTjj, (!blindSR)*weight);
                                VBF_MET_dEta_3JV_dPhiPTjj_sel->Fill(MET, (!blindSR)*weight);
                                VBF_METsoft_dEta_3JV_dPhiPTjj_sel->Fill(METsoft, (!blindSR)*weight);
                                VBF_METsig_dEta_3JV_dPhiPTjj_sel->Fill(METsig, (!blindSR)*weight);
                                VBF_MHTsig_dEta_3JV_dPhiPTjj_sel->Fill(MHTsig, (!blindSR)*weight);
                                VBF_minDeltaPhiPTj12_dEta_3JV_dPhiPTjj_sel->Fill(DPhiMET1, (!blindSR)*weight);
                                VBF_maxDeltaPhiPTj12_dEta_3JV_dPhiPTjj_sel->Fill(DPhiMET2, (!blindSR)*weight);
                                VBF_DeltaPhiPTj3_dEta_3JV_dPhiPTjj_sel->Fill(DPhiMET3, (!blindSR)*weight);
                            }

                        }

                    }

                }

            }

            if (HTMHT) {

                // ------------------------------------------------------------- //
                // fill preselection histos
                // ------------------------------------------------------------- //

                if( NJets >=2 ) {
                    if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                        HT_presel_pred_raw->Fill(HT, Ntries, weight);
                        MHT_presel_pred_raw->Fill(MHT, Ntries, weight);
                        MET_presel_pred_raw->Fill(MET, Ntries, weight);
                        NJets_presel_pred_raw->Fill(NJets, Ntries, weight);
                        NBJets_presel_pred_raw->Fill(BTags, Ntries, weight);
                        Jet1Pt_presel_pred_raw->Fill(JetPt->at(0), Ntries, weight);
                        Jet2Pt_presel_pred_raw->Fill(JetPt->at(1), Ntries, weight);
                        Jet3Pt_presel_pred_raw->Fill(JetPt->at(2), Ntries, weight);
                        Jet1Eta_presel_pred_raw->Fill(JetEta->at(0), Ntries, weight);
                        Jet2Eta_presel_pred_raw->Fill(JetEta->at(1), Ntries, weight);
                        Jet3Eta_presel_pred_raw->Fill(JetEta->at(2), Ntries, weight);
                        DeltaPhi1_presel_pred_raw->Fill(DeltaPhi->at(0), Ntries, weight);
                        DeltaPhi2_presel_pred_raw->Fill(DeltaPhi->at(1), Ntries, weight);
                        DeltaPhi3_presel_pred_raw->Fill(DeltaPhi->at(2), Ntries, weight);
                    } else if (Ntries == -2) {
                        HT_presel_sel->Fill(HT, weight);
                        MHT_presel_sel->Fill(MHT, weight);
                        MET_presel_sel->Fill(MET, weight);
                        NJets_presel_sel->Fill(NJets, weight);
                        NBJets_presel_sel->Fill(BTags, weight);
                        Jet1Pt_presel_sel->Fill(JetPt->at(0), weight);
                        Jet2Pt_presel_sel->Fill(JetPt->at(1), weight);
                        Jet3Pt_presel_sel->Fill(JetPt->at(2), weight);
                        Jet1Eta_presel_sel->Fill(JetEta->at(0), weight);
                        Jet2Eta_presel_sel->Fill(JetEta->at(1), weight);
                        Jet3Eta_presel_sel->Fill(JetEta->at(2), weight);
                        DeltaPhi1_presel_sel->Fill(DeltaPhi->at(0), weight);
                        DeltaPhi2_presel_sel->Fill(DeltaPhi->at(1), weight);
                        DeltaPhi3_presel_sel->Fill(DeltaPhi->at(2), weight);
                    }
                }

                // ------------------------------------------------------------- //
                // baseline closure without deltaPhi
                // ------------------------------------------------------------- //
                // HT baseline cut
                if( HT > 500. ) {
                    if( NJets >= 4 ) {
                        if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                            NJets_baseline_withoutDeltaPhi_withoutMET_pred_raw->Fill(NJets, Ntries, weight);
                            NBJets_baseline_withoutDeltaPhi_withoutMET_pred_raw->Fill(BTags, Ntries, weight);
                        } else if (Ntries == -2) {
                            NJets_baseline_withoutDeltaPhi_withoutMET_sel->Fill(NJets, weight);
                            NBJets_baseline_withoutDeltaPhi_withoutMET_sel->Fill(BTags, weight);
                        }
                    }
                    if( NJets == 2 || NJets == 3 ) {
                        if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                            MHT_JetBin1_baseline_withoutDeltaPhi_pred_raw->Fill(MHT, Ntries, weight);
                            MET_JetBin1_baseline_withoutDeltaPhi_pred_raw->Fill(MET, Ntries, weight);
                        } else if (Ntries == -2) {
                            MHT_JetBin1_baseline_withoutDeltaPhi_sel->Fill(MHT, weight);
                            MET_JetBin1_baseline_withoutDeltaPhi_sel->Fill(MET, weight);
                        }
                    }
                    else if( NJets == 4 ||  NJets == 5 || NJets == 6) {
                        if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                            MHT_JetBin2_baseline_withoutDeltaPhi_pred_raw->Fill(MHT, Ntries, weight);
                            MET_JetBin2_baseline_withoutDeltaPhi_pred_raw->Fill(MET, Ntries, weight);
                        } else if (Ntries == -2) {
                            MHT_JetBin2_baseline_withoutDeltaPhi_sel->Fill(MHT, weight);
                            MET_JetBin2_baseline_withoutDeltaPhi_sel->Fill(MET, weight);
                        }
                    }
                    else if( NJets == 7  ||  NJets == 8 ) {
                        if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                            MHT_JetBin3_baseline_withoutDeltaPhi_pred_raw->Fill(MHT, Ntries, weight);
                            MET_JetBin3_baseline_withoutDeltaPhi_pred_raw->Fill(MET, Ntries, weight);
                        } else if (Ntries == -2) {
                            MHT_JetBin3_baseline_withoutDeltaPhi_sel->Fill(MHT, weight);
                            MET_JetBin3_baseline_withoutDeltaPhi_sel->Fill(MET, weight);
                        }
                    }
                    else if( NJets >= 9 ) {
                        if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                            MHT_JetBin4_baseline_withoutDeltaPhi_pred_raw->Fill(MHT, Ntries, weight);
                            MET_JetBin4_baseline_withoutDeltaPhi_pred_raw->Fill(MET, Ntries, weight);
                        } else if (Ntries == -2) {
                            MHT_JetBin4_baseline_withoutDeltaPhi_sel->Fill(MHT, weight);
                            MET_JetBin4_baseline_withoutDeltaPhi_sel->Fill(MET, weight);
                        }
                    }

                    // MHT baseline cut
                    if( MET > 200. ) {
                        if( NJets >= 4 ) {
                            if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                                NJets_baseline_withoutDeltaPhi_pred_raw->Fill(NJets, Ntries, weight);
                                NBJets_baseline_withoutDeltaPhi_pred_raw->Fill(BTags, Ntries, weight);
                            } else if (Ntries == -2) {
                                NJets_baseline_withoutDeltaPhi_sel->Fill(NJets, weight);
                                NBJets_baseline_withoutDeltaPhi_sel->Fill(BTags, weight);
                            }
                        }
                        // jet bin 1
                        if( NJets == 2 || NJets == 3) {
                            if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                                Jet1Pt_JetBin1_baseline_withoutDeltaPhi_pred_raw->Fill(JetPt->at(0), Ntries, weight);
                                Jet2Pt_JetBin1_baseline_withoutDeltaPhi_pred_raw->Fill(JetPt->at(1), Ntries, weight);
                                Jet1Eta_JetBin1_baseline_withoutDeltaPhi_pred_raw->Fill(JetEta->at(0), Ntries, weight);
                                Jet2Eta_JetBin1_baseline_withoutDeltaPhi_pred_raw->Fill(JetEta->at(1), Ntries, weight);
                                DeltaPhi1_JetBin1_baseline_withoutDeltaPhi_pred_raw->Fill(DeltaPhi->at(0), Ntries, weight);
                                DeltaPhi2_JetBin1_baseline_withoutDeltaPhi_pred_raw->Fill(DeltaPhi->at(1), Ntries, weight);
                            } else if (Ntries == -2) {
                                Jet1Pt_JetBin1_baseline_withoutDeltaPhi_sel->Fill(JetPt->at(0), weight);
                                Jet2Pt_JetBin1_baseline_withoutDeltaPhi_sel->Fill(JetPt->at(1), weight);
                                Jet1Eta_JetBin1_baseline_withoutDeltaPhi_sel->Fill(JetEta->at(0), weight);
                                Jet2Eta_JetBin1_baseline_withoutDeltaPhi_sel->Fill(JetEta->at(1), weight);
                                DeltaPhi1_JetBin1_baseline_withoutDeltaPhi_sel->Fill(DeltaPhi->at(0), weight);
                                DeltaPhi2_JetBin1_baseline_withoutDeltaPhi_sel->Fill(DeltaPhi->at(1), weight);
                            }
                        }
                        // jet bin 2
                        else if( NJets == 4 || NJets == 5 ||  NJets == 6 ) {
                            if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                                Jet1Pt_JetBin2_baseline_withoutDeltaPhi_pred_raw->Fill(JetPt->at(0), Ntries, weight);
                                Jet2Pt_JetBin2_baseline_withoutDeltaPhi_pred_raw->Fill(JetPt->at(1), Ntries, weight);
                                Jet3Pt_JetBin2_baseline_withoutDeltaPhi_pred_raw->Fill(JetPt->at(2), Ntries, weight);
                                Jet1Eta_JetBin2_baseline_withoutDeltaPhi_pred_raw->Fill(JetEta->at(0), Ntries, weight);
                                Jet2Eta_JetBin2_baseline_withoutDeltaPhi_pred_raw->Fill(JetEta->at(1), Ntries, weight);
                                Jet3Eta_JetBin2_baseline_withoutDeltaPhi_pred_raw->Fill(JetEta->at(2), Ntries, weight);
                                DeltaPhi1_JetBin2_baseline_withoutDeltaPhi_pred_raw->Fill(DeltaPhi->at(0), Ntries, weight);
                                DeltaPhi2_JetBin2_baseline_withoutDeltaPhi_pred_raw->Fill(DeltaPhi->at(1), Ntries, weight);
                                DeltaPhi3_JetBin2_baseline_withoutDeltaPhi_pred_raw->Fill(DeltaPhi->at(2), Ntries, weight);
                            } else if (Ntries == -2) {
                                Jet1Pt_JetBin2_baseline_withoutDeltaPhi_sel->Fill(JetPt->at(0), weight);
                                Jet2Pt_JetBin2_baseline_withoutDeltaPhi_sel->Fill(JetPt->at(1), weight);
                                Jet3Pt_JetBin2_baseline_withoutDeltaPhi_sel->Fill(JetPt->at(2), weight);
                                Jet1Eta_JetBin2_baseline_withoutDeltaPhi_sel->Fill(JetEta->at(0), weight);
                                Jet2Eta_JetBin2_baseline_withoutDeltaPhi_sel->Fill(JetEta->at(1), weight);
                                Jet3Eta_JetBin2_baseline_withoutDeltaPhi_sel->Fill(JetEta->at(2), weight);
                                DeltaPhi1_JetBin2_baseline_withoutDeltaPhi_sel->Fill(DeltaPhi->at(0), weight);
                                DeltaPhi2_JetBin2_baseline_withoutDeltaPhi_sel->Fill(DeltaPhi->at(1), weight);
                                DeltaPhi3_JetBin2_baseline_withoutDeltaPhi_sel->Fill(DeltaPhi->at(2), weight);
                            }
                        }
                        // jet bin 3
                        else if( NJets == 7  ||  NJets == 8 ) {
                            if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                                Jet1Pt_JetBin3_baseline_withoutDeltaPhi_pred_raw->Fill(JetPt->at(0), Ntries, weight);
                                Jet2Pt_JetBin3_baseline_withoutDeltaPhi_pred_raw->Fill(JetPt->at(1), Ntries, weight);
                                Jet3Pt_JetBin3_baseline_withoutDeltaPhi_pred_raw->Fill(JetPt->at(2), Ntries, weight);
                                Jet1Eta_JetBin3_baseline_withoutDeltaPhi_pred_raw->Fill(JetEta->at(0), Ntries, weight);
                                Jet2Eta_JetBin3_baseline_withoutDeltaPhi_pred_raw->Fill(JetEta->at(1), Ntries, weight);
                                Jet3Eta_JetBin3_baseline_withoutDeltaPhi_pred_raw->Fill(JetEta->at(2), Ntries, weight);
                                DeltaPhi1_JetBin3_baseline_withoutDeltaPhi_pred_raw->Fill(DeltaPhi->at(0), Ntries, weight);
                                DeltaPhi2_JetBin3_baseline_withoutDeltaPhi_pred_raw->Fill(DeltaPhi->at(1), Ntries, weight);
                                DeltaPhi3_JetBin3_baseline_withoutDeltaPhi_pred_raw->Fill(DeltaPhi->at(2), Ntries, weight);
                            } else if (Ntries == -2) {
                                Jet1Pt_JetBin3_baseline_withoutDeltaPhi_sel->Fill(JetPt->at(0), weight);
                                Jet2Pt_JetBin3_baseline_withoutDeltaPhi_sel->Fill(JetPt->at(1), weight);
                                Jet3Pt_JetBin3_baseline_withoutDeltaPhi_sel->Fill(JetPt->at(2), weight);
                                Jet1Eta_JetBin3_baseline_withoutDeltaPhi_sel->Fill(JetEta->at(0), weight);
                                Jet2Eta_JetBin3_baseline_withoutDeltaPhi_sel->Fill(JetEta->at(1), weight);
                                Jet3Eta_JetBin3_baseline_withoutDeltaPhi_sel->Fill(JetEta->at(2), weight);
                                DeltaPhi1_JetBin3_baseline_withoutDeltaPhi_sel->Fill(DeltaPhi->at(0), weight);
                                DeltaPhi2_JetBin3_baseline_withoutDeltaPhi_sel->Fill(DeltaPhi->at(1), weight);
                                DeltaPhi3_JetBin3_baseline_withoutDeltaPhi_sel->Fill(DeltaPhi->at(2), weight);
                            }
                        }
                        // jet bin 4
                        else if( NJets >= 9 ) {
                            if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                                Jet1Pt_JetBin4_baseline_withoutDeltaPhi_pred_raw->Fill(JetPt->at(0), Ntries, weight);
                                Jet2Pt_JetBin4_baseline_withoutDeltaPhi_pred_raw->Fill(JetPt->at(1), Ntries, weight);
                                Jet3Pt_JetBin4_baseline_withoutDeltaPhi_pred_raw->Fill(JetPt->at(2), Ntries, weight);
                                Jet1Eta_JetBin4_baseline_withoutDeltaPhi_pred_raw->Fill(JetEta->at(0), Ntries, weight);
                                Jet2Eta_JetBin4_baseline_withoutDeltaPhi_pred_raw->Fill(JetEta->at(1), Ntries, weight);
                                Jet3Eta_JetBin4_baseline_withoutDeltaPhi_pred_raw->Fill(JetEta->at(2), Ntries, weight);
                                DeltaPhi1_JetBin4_baseline_withoutDeltaPhi_pred_raw->Fill(DeltaPhi->at(0), Ntries, weight);
                                DeltaPhi2_JetBin4_baseline_withoutDeltaPhi_pred_raw->Fill(DeltaPhi->at(1), Ntries, weight);
                                DeltaPhi3_JetBin4_baseline_withoutDeltaPhi_pred_raw->Fill(DeltaPhi->at(2), Ntries, weight);
                            } else if (Ntries == -2) {
                                Jet1Pt_JetBin4_baseline_withoutDeltaPhi_sel->Fill(JetPt->at(0), weight);
                                Jet2Pt_JetBin4_baseline_withoutDeltaPhi_sel->Fill(JetPt->at(1), weight);
                                Jet3Pt_JetBin4_baseline_withoutDeltaPhi_sel->Fill(JetPt->at(2), weight);
                                Jet1Eta_JetBin4_baseline_withoutDeltaPhi_sel->Fill(JetEta->at(0), weight);
                                Jet2Eta_JetBin4_baseline_withoutDeltaPhi_sel->Fill(JetEta->at(1), weight);
                                Jet3Eta_JetBin4_baseline_withoutDeltaPhi_sel->Fill(JetEta->at(2), weight);
                                DeltaPhi1_JetBin4_baseline_withoutDeltaPhi_sel->Fill(DeltaPhi->at(0), weight);
                                DeltaPhi2_JetBin4_baseline_withoutDeltaPhi_sel->Fill(DeltaPhi->at(1), weight);
                                DeltaPhi3_JetBin4_baseline_withoutDeltaPhi_sel->Fill(DeltaPhi->at(2), weight);
                            }
                        }
                    }
                }

                // MHT baseline cut
                if( MET > 200. ) {
                    if( NJets == 2 || NJets == 3 ) {
                        if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                            HT_JetBin1_baseline_withoutDeltaPhi_pred_raw->Fill(HT, Ntries, weight);
                        } else if (Ntries == -2) {
                            HT_JetBin1_baseline_withoutDeltaPhi_sel->Fill(HT, weight);
                        }
                    }
                    else if( NJets == 4 || NJets == 5 ||  NJets == 6 ) {
                        if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                            HT_JetBin2_baseline_withoutDeltaPhi_pred_raw->Fill(HT, Ntries, weight);
                        } else if (Ntries == -2) {
                            HT_JetBin2_baseline_withoutDeltaPhi_sel->Fill(HT, weight);
                        }
                    }
                    else if(  NJets == 7  ||  NJets == 8 ) {
                        if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                            HT_JetBin3_baseline_withoutDeltaPhi_pred_raw->Fill(HT, Ntries, weight);
                        } else if (Ntries == -2) {
                            HT_JetBin3_baseline_withoutDeltaPhi_sel->Fill(HT, weight);
                        }
                    }
                    else if(  NJets >= 9 ) {
                        if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                            HT_JetBin4_baseline_withoutDeltaPhi_pred_raw->Fill(HT, Ntries, weight);
                        } else if (Ntries == -2) {
                            HT_JetBin4_baseline_withoutDeltaPhi_sel->Fill(HT, weight);
                        }
                    }
                }


                // ------------------------------------------------------------- //
                // check deltaPhi cut
                if( DeltaPhiCut() ) {

                    // ------------------------------------------------------------- //
                    // fill histos after deltaPhi cut
                    // ------------------------------------------------------------- //
                    if( NJets >= 4 ) {
                        if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                            HT_deltaPhi_pred_raw->Fill(HT, Ntries, weight);
                            MHT_deltaPhi_pred_raw->Fill(MHT, Ntries, weight);
                            MET_deltaPhi_pred_raw->Fill(MET, Ntries, weight);
                            Jet1Pt_deltaPhi_pred_raw->Fill(JetPt->at(0), Ntries, weight);
                            Jet2Pt_deltaPhi_pred_raw->Fill(JetPt->at(1), Ntries, weight);
                            Jet3Pt_deltaPhi_pred_raw->Fill(JetPt->at(2), Ntries, weight);
                            Jet1Eta_deltaPhi_pred_raw->Fill(JetEta->at(0), Ntries, weight);
                            Jet2Eta_deltaPhi_pred_raw->Fill(JetEta->at(1), Ntries, weight);
                            Jet3Eta_deltaPhi_pred_raw->Fill(JetEta->at(2), Ntries, weight);
                        } else if (Ntries == -2) {
                            HT_deltaPhi_sel->Fill(HT, weight);
                            MHT_deltaPhi_sel->Fill(MHT, weight);
                            MET_deltaPhi_sel->Fill(MET, weight);
                            Jet1Pt_deltaPhi_sel->Fill(JetPt->at(0), weight);
                            Jet2Pt_deltaPhi_sel->Fill(JetPt->at(1), weight);
                            Jet3Pt_deltaPhi_sel->Fill(JetPt->at(2), weight);
                            Jet1Eta_deltaPhi_sel->Fill(JetEta->at(0), weight);
                            Jet2Eta_deltaPhi_sel->Fill(JetEta->at(1), weight);
                            Jet3Eta_deltaPhi_sel->Fill(JetEta->at(2), weight);
                        }

                        // ------------------------------------------------------------- //
                        // fill baseline histos
                        // ------------------------------------------------------------- //
                        if( HT > 500. ) {
                            if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                                MHT_baseline_pred_raw->Fill(MHT, Ntries, weight);
                                MET_baseline_pred_raw->Fill(MET, Ntries, weight);
                                NJets_baseline_withoutMET_pred_raw->Fill(NJets, Ntries, weight);
                                NBJets_baseline_withoutMET_pred_raw->Fill(BTags, Ntries, weight);
                            } else if (Ntries == -2) {
                                MHT_baseline_sel->Fill(MHT, weight);
                                MET_baseline_sel->Fill(MET, weight);
                                NJets_baseline_withoutMET_sel->Fill(NJets, weight);
                                NBJets_baseline_withoutMET_sel->Fill(BTags, weight);
                            }
                            if( MET > 200. ) {
                                if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                                    NJets_baseline_pred_raw->Fill(NJets, Ntries, weight);
                                    NBJets_baseline_pred_raw->Fill(BTags, Ntries, weight);
                                } else if (Ntries == -2) {
                                    NJets_baseline_sel->Fill(NJets, weight);
                                    NBJets_baseline_sel->Fill(BTags, weight);
                                }
                            }
                        }
                        if( MET > 200. ) {
                            if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                                HT_baseline_pred_raw->Fill(HT, Ntries, weight);
                            } else if (Ntries == -2) {
                                HT_baseline_sel->Fill(HT, weight);
                            }
                        }
                    }

                    // ------------------------------------------------------------- //
                    //  fill different jet multiplicity bins
                    // ------------------------------------------------------------- //
                    // jet bin 1
                    if( NJets == 2 || NJets == 3) {
                        if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                            MHT_JetBin1_HTinclusive_pred_raw->Fill(MHT, Ntries, weight);
                            MET_JetBin1_HTinclusive_pred_raw->Fill(MET, Ntries, weight);
                        } else if (Ntries == -2) {
                            MHT_JetBin1_HTinclusive_sel->Fill(MHT, weight);
                            MET_JetBin1_HTinclusive_sel->Fill(MET, weight);
                        }

                        if( MET > 200. ) {
                            if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                                Jet1Pt_JetBin1_baseline_pred_raw->Fill(JetPt->at(0), Ntries, weight);
                                Jet2Pt_JetBin1_baseline_pred_raw->Fill(JetPt->at(1), Ntries, weight);
                                Jet1Eta_JetBin1_baseline_pred_raw->Fill(JetEta->at(0), Ntries, weight);
                                Jet2Eta_JetBin1_baseline_pred_raw->Fill(JetEta->at(1), Ntries, weight);
                                DeltaPhi1_JetBin1_baseline_pred_raw->Fill(DeltaPhi->at(0), Ntries, weight);
                                DeltaPhi2_JetBin1_baseline_pred_raw->Fill(DeltaPhi->at(1), Ntries, weight);
                            } else if (Ntries == -2) {
                                Jet1Pt_JetBin1_baseline_sel->Fill(JetPt->at(0), weight);
                                Jet2Pt_JetBin1_baseline_sel->Fill(JetPt->at(1), weight);
                                Jet1Eta_JetBin1_baseline_sel->Fill(JetEta->at(0), weight);
                                Jet2Eta_JetBin1_baseline_sel->Fill(JetEta->at(1), weight);
                                DeltaPhi1_JetBin1_baseline_sel->Fill(DeltaPhi->at(0), weight);
                                DeltaPhi2_JetBin1_baseline_sel->Fill(DeltaPhi->at(1), weight);
                            }
                        }

                    }
                    // jet bin 2
                    else if( NJets == 4 || NJets == 5 ||  NJets == 6 ) {
                        if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                            MHT_JetBin2_HTinclusive_pred_raw->Fill(MHT, Ntries, weight);
                            MET_JetBin2_HTinclusive_pred_raw->Fill(MET, Ntries, weight);
                        } else if (Ntries == -2) {
                            MHT_JetBin2_HTinclusive_sel->Fill(MHT, weight);
                            MET_JetBin2_HTinclusive_sel->Fill(MET, weight);
                        }

                        if( MET > 200. ) {
                            if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                                Jet1Pt_JetBin2_baseline_pred_raw->Fill(JetPt->at(0), Ntries, weight);
                                Jet2Pt_JetBin2_baseline_pred_raw->Fill(JetPt->at(1), Ntries, weight);
                                Jet3Pt_JetBin2_baseline_pred_raw->Fill(JetPt->at(2), Ntries, weight);
                                Jet1Eta_JetBin2_baseline_pred_raw->Fill(JetEta->at(0), Ntries, weight);
                                Jet2Eta_JetBin2_baseline_pred_raw->Fill(JetEta->at(1), Ntries, weight);
                                Jet3Eta_JetBin2_baseline_pred_raw->Fill(JetEta->at(2), Ntries, weight);
                                DeltaPhi1_JetBin2_baseline_pred_raw->Fill(DeltaPhi->at(0), Ntries, weight);
                                DeltaPhi2_JetBin2_baseline_pred_raw->Fill(DeltaPhi->at(1), Ntries, weight);
                                DeltaPhi3_JetBin2_baseline_pred_raw->Fill(DeltaPhi->at(2), Ntries, weight);
                            } else if (Ntries == -2) {
                                Jet1Pt_JetBin2_baseline_sel->Fill(JetPt->at(0), weight);
                                Jet2Pt_JetBin2_baseline_sel->Fill(JetPt->at(1), weight);
                                Jet3Pt_JetBin2_baseline_sel->Fill(JetPt->at(2), weight);
                                Jet1Eta_JetBin2_baseline_sel->Fill(JetEta->at(0), weight);
                                Jet2Eta_JetBin2_baseline_sel->Fill(JetEta->at(1), weight);
                                Jet3Eta_JetBin2_baseline_sel->Fill(JetEta->at(2), weight);
                                DeltaPhi1_JetBin2_baseline_sel->Fill(DeltaPhi->at(0), weight);
                                DeltaPhi2_JetBin2_baseline_sel->Fill(DeltaPhi->at(1), weight);
                                DeltaPhi3_JetBin2_baseline_sel->Fill(DeltaPhi->at(2), weight);
                            }
                        }

                    }
                    // jet bin 3
                    else if(  NJets == 7  ||  NJets == 8 ) {
                        if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                            MHT_JetBin3_HTinclusive_pred_raw->Fill(MHT, Ntries, weight);
                            MET_JetBin3_HTinclusive_pred_raw->Fill(MET, Ntries, weight);
                        } else if (Ntries == -2) {
                            MHT_JetBin3_HTinclusive_sel->Fill(MHT, weight);
                            MET_JetBin3_HTinclusive_sel->Fill(MET, weight);
                        }

                        if( MET > 200. ) {
                            if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax) {
                                Jet1Pt_JetBin3_baseline_pred_raw->Fill(JetPt->at(0), Ntries, weight);
                                Jet2Pt_JetBin3_baseline_pred_raw->Fill(JetPt->at(1), Ntries, weight);
                                Jet3Pt_JetBin3_baseline_pred_raw->Fill(JetPt->at(2), Ntries, weight);
                                Jet1Eta_JetBin3_baseline_pred_raw->Fill(JetEta->at(0), Ntries, weight);
                                Jet2Eta_JetBin3_baseline_pred_raw->Fill(JetEta->at(1), Ntries, weight);
                                Jet3Eta_JetBin3_baseline_pred_raw->Fill(JetEta->at(2), Ntries, weight);
                                DeltaPhi1_JetBin3_baseline_pred_raw->Fill(DeltaPhi->at(0), Ntries, weight);
                                DeltaPhi2_JetBin3_baseline_pred_raw->Fill(DeltaPhi->at(1), Ntries, weight);
                                DeltaPhi3_JetBin3_baseline_pred_raw->Fill(DeltaPhi->at(2), Ntries, weight);
                            } else if (Ntries == -2) {
                                Jet1Pt_JetBin3_baseline_sel->Fill(JetPt->at(0), weight);
                                Jet2Pt_JetBin3_baseline_sel->Fill(JetPt->at(1), weight);
                                Jet3Pt_JetBin3_baseline_sel->Fill(JetPt->at(2), weight);
                                Jet1Eta_JetBin3_baseline_sel->Fill(JetEta->at(0), weight);
                                Jet2Eta_JetBin3_baseline_sel->Fill(JetEta->at(1), weight);
                                Jet3Eta_JetBin3_baseline_sel->Fill(JetEta->at(2), weight);
                                DeltaPhi1_JetBin3_baseline_sel->Fill(DeltaPhi->at(0), weight);
                                DeltaPhi2_JetBin3_baseline_sel->Fill(DeltaPhi->at(1), weight);
                                DeltaPhi3_JetBin3_baseline_sel->Fill(DeltaPhi->at(2), weight);
                            }
                        }

                    }
                    // jet bin 4
                    else if(  NJets >= 9 ) {
                        if (Ntries > 0 && METsig < METSigSeedMax && MHTsig < MHTSigSeedMax && METsoft < METSoftSeedMax) {
                            MHT_JetBin4_HTinclusive_pred_raw->Fill(MHT, Ntries, weight);
                            MET_JetBin4_HTinclusive_pred_raw->Fill(MET, Ntries, weight);
                        } else if (Ntries == -2) {
                            MHT_JetBin4_HTinclusive_sel->Fill(MHT, weight);
                            MET_JetBin4_HTinclusive_sel->Fill(MET, weight);
                        }

                        if( MET > 200. ) {
                            if (Ntries > 0) {
                                Jet1Pt_JetBin4_baseline_pred_raw->Fill(JetPt->at(0), Ntries, weight);
                                Jet2Pt_JetBin4_baseline_pred_raw->Fill(JetPt->at(1), Ntries, weight);
                                Jet3Pt_JetBin4_baseline_pred_raw->Fill(JetPt->at(2), Ntries, weight);
                                Jet1Eta_JetBin4_baseline_pred_raw->Fill(JetEta->at(0), Ntries, weight);
                                Jet2Eta_JetBin4_baseline_pred_raw->Fill(JetEta->at(1), Ntries, weight);
                                Jet3Eta_JetBin4_baseline_pred_raw->Fill(JetEta->at(2), Ntries, weight);
                                DeltaPhi1_JetBin4_baseline_pred_raw->Fill(DeltaPhi->at(0), Ntries, weight);
                                DeltaPhi2_JetBin4_baseline_pred_raw->Fill(DeltaPhi->at(1), Ntries, weight);
                                DeltaPhi3_JetBin4_baseline_pred_raw->Fill(DeltaPhi->at(2), Ntries, weight);
                            } else if (Ntries == -2) {
                                Jet1Pt_JetBin4_baseline_sel->Fill(JetPt->at(0), weight);
                                Jet2Pt_JetBin4_baseline_sel->Fill(JetPt->at(1), weight);
                                Jet3Pt_JetBin4_baseline_sel->Fill(JetPt->at(2), weight);
                                Jet1Eta_JetBin4_baseline_sel->Fill(JetEta->at(0), weight);
                                Jet2Eta_JetBin4_baseline_sel->Fill(JetEta->at(1), weight);
                                Jet3Eta_JetBin4_baseline_sel->Fill(JetEta->at(2), weight);
                                DeltaPhi1_JetBin4_baseline_sel->Fill(DeltaPhi->at(0), weight);
                                DeltaPhi2_JetBin4_baseline_sel->Fill(DeltaPhi->at(1), weight);
                                DeltaPhi3_JetBin4_baseline_sel->Fill(DeltaPhi->at(2), weight);
                            }
                        }
                    }
                }
            }
        }
    }

    cout << "after filling all histograms" << endl;
    //----------------------------------------------------------//


    //----------------------------------------------------------//
    // rebin histos
    // ------------------------------------------------------------- //
    DoRebinning(HT_presel_pred_raw, HT_presel_sel , 2);
    DoRebinning(MHT_presel_pred_raw, MHT_presel_sel , -2);
    DoRebinning(MET_presel_pred_raw, MET_presel_sel , -2);
    DoRebinning(Jet1Pt_presel_pred_raw, Jet1Pt_presel_sel , 2);
    DoRebinning(Jet2Pt_presel_pred_raw, Jet2Pt_presel_sel , 2);
    DoRebinning(Jet3Pt_presel_pred_raw, Jet3Pt_presel_sel , 2);
    DoRebinning(Jet1Eta_presel_pred_raw, Jet1Eta_presel_sel , 2);
    DoRebinning(Jet2Eta_presel_pred_raw, Jet2Eta_presel_sel , 2);
    DoRebinning(Jet3Eta_presel_pred_raw, Jet3Eta_presel_sel , 2);
    DoRebinning(DeltaPhi1_presel_pred_raw, DeltaPhi1_presel_sel , 5);
    DoRebinning(DeltaPhi2_presel_pred_raw, DeltaPhi2_presel_sel , 5);
    DoRebinning(DeltaPhi3_presel_pred_raw, DeltaPhi3_presel_sel , 5);

    DoRebinning(HT_deltaPhi_pred_raw, HT_deltaPhi_sel , 2);
    DoRebinning(MHT_deltaPhi_pred_raw, MHT_deltaPhi_sel , -2);
    DoRebinning(MET_deltaPhi_pred_raw, MET_deltaPhi_sel , -2);
    DoRebinning(Jet1Pt_deltaPhi_pred_raw, Jet1Pt_deltaPhi_sel , 2);
    DoRebinning(Jet2Pt_deltaPhi_pred_raw, Jet2Pt_deltaPhi_sel , 2);
    DoRebinning(Jet3Pt_deltaPhi_pred_raw, Jet3Pt_deltaPhi_sel , 2);
    DoRebinning(Jet1Eta_deltaPhi_pred_raw, Jet1Eta_deltaPhi_sel , 2);
    DoRebinning(Jet2Eta_deltaPhi_pred_raw, Jet2Eta_deltaPhi_sel , 2);
    DoRebinning(Jet3Eta_deltaPhi_pred_raw, Jet3Eta_deltaPhi_sel , 2);

    DoRebinning(MHT_JetBin1_HTinclusive_pred_raw, MHT_JetBin1_HTinclusive_sel , -2);
    DoRebinning(MHT_JetBin2_HTinclusive_pred_raw, MHT_JetBin2_HTinclusive_sel , -2);
    DoRebinning(MHT_JetBin3_HTinclusive_pred_raw, MHT_JetBin3_HTinclusive_sel , -2);
    DoRebinning(MHT_JetBin4_HTinclusive_pred_raw, MHT_JetBin4_HTinclusive_sel , -2);

    DoRebinning(MET_JetBin1_HTinclusive_pred_raw, MET_JetBin1_HTinclusive_sel , -2);
    DoRebinning(MET_JetBin2_HTinclusive_pred_raw, MET_JetBin2_HTinclusive_sel , -2);
    DoRebinning(MET_JetBin3_HTinclusive_pred_raw, MET_JetBin3_HTinclusive_sel , -2);
    DoRebinning(MET_JetBin4_HTinclusive_pred_raw, MET_JetBin4_HTinclusive_sel , -2);

    DoRebinning(HT_baseline_pred_raw, HT_baseline_sel , 2);
    DoRebinning(MHT_baseline_pred_raw, MHT_baseline_sel , -2);
    DoRebinning(MET_baseline_pred_raw, MET_baseline_sel , -2);

    DoRebinning(Jet1Pt_JetBin1_baseline_pred_raw, Jet1Pt_JetBin1_baseline_sel, 2);
    DoRebinning(Jet2Pt_JetBin1_baseline_pred_raw, Jet2Pt_JetBin1_baseline_sel, 2);
    DoRebinning(Jet1Eta_JetBin1_baseline_pred_raw, Jet1Eta_JetBin1_baseline_sel, 2);
    DoRebinning(Jet2Eta_JetBin1_baseline_pred_raw, Jet2Eta_JetBin1_baseline_sel, 2);
    DoRebinning(DeltaPhi1_JetBin1_baseline_pred_raw, DeltaPhi1_JetBin1_baseline_sel, 5);
    DoRebinning(DeltaPhi2_JetBin1_baseline_pred_raw, DeltaPhi2_JetBin1_baseline_sel, 5);

    DoRebinning(Jet1Pt_JetBin2_baseline_pred_raw, Jet1Pt_JetBin2_baseline_sel, 2);
    DoRebinning(Jet2Pt_JetBin2_baseline_pred_raw, Jet2Pt_JetBin2_baseline_sel, 2);
    DoRebinning(Jet3Pt_JetBin2_baseline_pred_raw, Jet3Pt_JetBin2_baseline_sel, 2);
    DoRebinning(Jet1Eta_JetBin2_baseline_pred_raw, Jet1Eta_JetBin2_baseline_sel, 2);
    DoRebinning(Jet2Eta_JetBin2_baseline_pred_raw, Jet2Eta_JetBin2_baseline_sel, 2);
    DoRebinning(Jet3Eta_JetBin2_baseline_pred_raw, Jet3Eta_JetBin2_baseline_sel, 2);
    DoRebinning(DeltaPhi1_JetBin2_baseline_pred_raw, DeltaPhi1_JetBin2_baseline_sel, 5);
    DoRebinning(DeltaPhi2_JetBin2_baseline_pred_raw, DeltaPhi2_JetBin2_baseline_sel, 5);
    DoRebinning(DeltaPhi3_JetBin2_baseline_pred_raw, DeltaPhi3_JetBin2_baseline_sel, 5);

    DoRebinning(Jet1Pt_JetBin3_baseline_pred_raw, Jet1Pt_JetBin3_baseline_sel, 2);
    DoRebinning(Jet2Pt_JetBin3_baseline_pred_raw, Jet2Pt_JetBin3_baseline_sel, 2);
    DoRebinning(Jet3Pt_JetBin3_baseline_pred_raw, Jet3Pt_JetBin3_baseline_sel, 2);
    DoRebinning(Jet1Eta_JetBin3_baseline_pred_raw, Jet1Eta_JetBin3_baseline_sel, 2);
    DoRebinning(Jet2Eta_JetBin3_baseline_pred_raw, Jet2Eta_JetBin3_baseline_sel, 2);
    DoRebinning(Jet3Eta_JetBin3_baseline_pred_raw, Jet3Eta_JetBin3_baseline_sel, 2);
    DoRebinning(DeltaPhi1_JetBin3_baseline_pred_raw, DeltaPhi1_JetBin3_baseline_sel, 5);
    DoRebinning(DeltaPhi2_JetBin3_baseline_pred_raw, DeltaPhi2_JetBin3_baseline_sel, 5);
    DoRebinning(DeltaPhi3_JetBin3_baseline_pred_raw, DeltaPhi3_JetBin3_baseline_sel, 5);

    DoRebinning(Jet1Pt_JetBin4_baseline_pred_raw, Jet1Pt_JetBin4_baseline_sel, 2);
    DoRebinning(Jet2Pt_JetBin4_baseline_pred_raw, Jet2Pt_JetBin4_baseline_sel, 2);
    DoRebinning(Jet3Pt_JetBin4_baseline_pred_raw, Jet3Pt_JetBin4_baseline_sel, 2);
    DoRebinning(Jet1Eta_JetBin4_baseline_pred_raw, Jet1Eta_JetBin4_baseline_sel, 2);
    DoRebinning(Jet2Eta_JetBin4_baseline_pred_raw, Jet2Eta_JetBin4_baseline_sel, 2);
    DoRebinning(Jet3Eta_JetBin4_baseline_pred_raw, Jet3Eta_JetBin4_baseline_sel, 2);
    DoRebinning(DeltaPhi1_JetBin4_baseline_pred_raw, DeltaPhi1_JetBin4_baseline_sel, 5);
    DoRebinning(DeltaPhi2_JetBin4_baseline_pred_raw, DeltaPhi2_JetBin4_baseline_sel, 5);
    DoRebinning(DeltaPhi3_JetBin4_baseline_pred_raw, DeltaPhi3_JetBin4_baseline_sel, 5);

    DoRebinning(HT_JetBin1_baseline_withoutDeltaPhi_pred_raw, HT_JetBin1_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(MHT_JetBin1_baseline_withoutDeltaPhi_pred_raw, MHT_JetBin1_baseline_withoutDeltaPhi_sel, -2);
    DoRebinning(MET_JetBin1_baseline_withoutDeltaPhi_pred_raw, MET_JetBin1_baseline_withoutDeltaPhi_sel, -2);
    DoRebinning(Jet1Pt_JetBin1_baseline_withoutDeltaPhi_pred_raw, Jet1Pt_JetBin1_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(Jet2Pt_JetBin1_baseline_withoutDeltaPhi_pred_raw, Jet2Pt_JetBin1_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(Jet1Eta_JetBin1_baseline_withoutDeltaPhi_pred_raw, Jet1Eta_JetBin1_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(Jet2Eta_JetBin1_baseline_withoutDeltaPhi_pred_raw, Jet2Eta_JetBin1_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(DeltaPhi1_JetBin1_baseline_withoutDeltaPhi_pred_raw, DeltaPhi1_JetBin1_baseline_withoutDeltaPhi_sel, 5);
    DoRebinning(DeltaPhi2_JetBin1_baseline_withoutDeltaPhi_pred_raw, DeltaPhi2_JetBin1_baseline_withoutDeltaPhi_sel, 5);

    DoRebinning(HT_JetBin2_baseline_withoutDeltaPhi_pred_raw, HT_JetBin2_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(MHT_JetBin2_baseline_withoutDeltaPhi_pred_raw, MHT_JetBin2_baseline_withoutDeltaPhi_sel, -2);
    DoRebinning(MET_JetBin2_baseline_withoutDeltaPhi_pred_raw, MET_JetBin2_baseline_withoutDeltaPhi_sel, -2);
    DoRebinning(Jet1Pt_JetBin2_baseline_withoutDeltaPhi_pred_raw, Jet1Pt_JetBin2_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(Jet2Pt_JetBin2_baseline_withoutDeltaPhi_pred_raw, Jet2Pt_JetBin2_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(Jet3Pt_JetBin2_baseline_withoutDeltaPhi_pred_raw, Jet3Pt_JetBin2_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(Jet1Eta_JetBin2_baseline_withoutDeltaPhi_pred_raw, Jet1Eta_JetBin2_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(Jet2Eta_JetBin2_baseline_withoutDeltaPhi_pred_raw, Jet2Eta_JetBin2_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(Jet3Eta_JetBin2_baseline_withoutDeltaPhi_pred_raw, Jet3Eta_JetBin2_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(DeltaPhi1_JetBin2_baseline_withoutDeltaPhi_pred_raw, DeltaPhi1_JetBin2_baseline_withoutDeltaPhi_sel, 5);
    DoRebinning(DeltaPhi2_JetBin2_baseline_withoutDeltaPhi_pred_raw, DeltaPhi2_JetBin2_baseline_withoutDeltaPhi_sel, 5);
    DoRebinning(DeltaPhi3_JetBin2_baseline_withoutDeltaPhi_pred_raw, DeltaPhi3_JetBin2_baseline_withoutDeltaPhi_sel, 5);

    DoRebinning(HT_JetBin3_baseline_withoutDeltaPhi_pred_raw, HT_JetBin3_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(MHT_JetBin3_baseline_withoutDeltaPhi_pred_raw, MHT_JetBin3_baseline_withoutDeltaPhi_sel, -2);
    DoRebinning(MET_JetBin3_baseline_withoutDeltaPhi_pred_raw, MET_JetBin3_baseline_withoutDeltaPhi_sel, -2);
    DoRebinning(Jet1Pt_JetBin3_baseline_withoutDeltaPhi_pred_raw, Jet1Pt_JetBin3_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(Jet2Pt_JetBin3_baseline_withoutDeltaPhi_pred_raw, Jet2Pt_JetBin3_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(Jet3Pt_JetBin3_baseline_withoutDeltaPhi_pred_raw, Jet3Pt_JetBin3_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(Jet1Eta_JetBin3_baseline_withoutDeltaPhi_pred_raw, Jet1Eta_JetBin3_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(Jet2Eta_JetBin3_baseline_withoutDeltaPhi_pred_raw, Jet2Eta_JetBin3_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(Jet3Eta_JetBin3_baseline_withoutDeltaPhi_pred_raw, Jet3Eta_JetBin3_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(DeltaPhi1_JetBin3_baseline_withoutDeltaPhi_pred_raw, DeltaPhi1_JetBin3_baseline_withoutDeltaPhi_sel, 5);
    DoRebinning(DeltaPhi2_JetBin3_baseline_withoutDeltaPhi_pred_raw, DeltaPhi2_JetBin3_baseline_withoutDeltaPhi_sel, 5);
    DoRebinning(DeltaPhi3_JetBin3_baseline_withoutDeltaPhi_pred_raw, DeltaPhi3_JetBin3_baseline_withoutDeltaPhi_sel, 5);

    DoRebinning(HT_JetBin4_baseline_withoutDeltaPhi_pred_raw, HT_JetBin4_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(MHT_JetBin4_baseline_withoutDeltaPhi_pred_raw, MHT_JetBin4_baseline_withoutDeltaPhi_sel, -2);
    DoRebinning(MET_JetBin4_baseline_withoutDeltaPhi_pred_raw, MET_JetBin4_baseline_withoutDeltaPhi_sel, -2);
    DoRebinning(Jet1Pt_JetBin4_baseline_withoutDeltaPhi_pred_raw, Jet1Pt_JetBin4_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(Jet2Pt_JetBin4_baseline_withoutDeltaPhi_pred_raw, Jet2Pt_JetBin4_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(Jet3Pt_JetBin4_baseline_withoutDeltaPhi_pred_raw, Jet3Pt_JetBin4_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(Jet1Eta_JetBin4_baseline_withoutDeltaPhi_pred_raw, Jet1Eta_JetBin4_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(Jet2Eta_JetBin4_baseline_withoutDeltaPhi_pred_raw, Jet2Eta_JetBin4_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(Jet3Eta_JetBin4_baseline_withoutDeltaPhi_pred_raw, Jet3Eta_JetBin4_baseline_withoutDeltaPhi_sel, 2);
    DoRebinning(DeltaPhi1_JetBin4_baseline_withoutDeltaPhi_pred_raw, DeltaPhi1_JetBin4_baseline_withoutDeltaPhi_sel, 5);
    DoRebinning(DeltaPhi2_JetBin4_baseline_withoutDeltaPhi_pred_raw, DeltaPhi2_JetBin4_baseline_withoutDeltaPhi_sel, 5);
    DoRebinning(DeltaPhi3_JetBin4_baseline_withoutDeltaPhi_pred_raw, DeltaPhi3_JetBin4_baseline_withoutDeltaPhi_sel, 5);

    DoRebinning(VBF_dPhi_presel_pred_raw, VBF_dPhi_presel_sel, 5);
    DoRebinning(VBF_dEta_presel_pred_raw, VBF_dEta_presel_sel, 5);
    DoRebinning(VBF_Mjj_presel_pred_raw, VBF_Mjj_presel_sel, -3);
    DoRebinning(VBF_Jet1Pt_presel_pred_raw, VBF_Jet1Pt_presel_sel, 2);
    DoRebinning(VBF_Jet2Pt_presel_pred_raw, VBF_Jet2Pt_presel_sel, 2);
    DoRebinning(VBF_Jet3Pt_presel_pred_raw, VBF_Jet3Pt_presel_sel, 2);
    DoRebinning(VBF_Jet1Eta_presel_pred_raw, VBF_Jet1Eta_presel_sel, 5);
    DoRebinning(VBF_Jet2Eta_presel_pred_raw, VBF_Jet2Eta_presel_sel, 5);
    DoRebinning(VBF_Jet3Eta_presel_pred_raw, VBF_Jet3Eta_presel_sel, 5);
    DoRebinning(VBF_PTjj_presel_pred_raw, VBF_PTjj_presel_sel, -2);
    DoRebinning(VBF_MET_presel_pred_raw, VBF_MET_presel_sel, -2);
    DoRebinning(VBF_METsoft_presel_pred_raw, VBF_METsoft_presel_sel, 2);
    DoRebinning(VBF_METsig_presel_pred_raw, VBF_METsig_presel_sel, 2);
    DoRebinning(VBF_MHTsig_presel_pred_raw, VBF_MHTsig_presel_sel, 2);
    DoRebinning(VBF_minDeltaPhiPTj12_presel_pred_raw, VBF_minDeltaPhiPTj12_presel_sel, 5);
    DoRebinning(VBF_maxDeltaPhiPTj12_presel_pred_raw, VBF_maxDeltaPhiPTj12_presel_sel, 5);
    DoRebinning(VBF_DeltaPhiPTj3_presel_pred_raw, VBF_DeltaPhiPTj3_presel_sel, 5);

    DoRebinning(VBF_dPhi_presel_4JV_dPhiSide_pred_raw, VBF_dPhi_presel_4JV_dPhiSide_sel, 5);
    DoRebinning(VBF_dEta_presel_4JV_dPhiSide_pred_raw, VBF_dEta_presel_4JV_dPhiSide_sel, 5);
    DoRebinning(VBF_Mjj_presel_4JV_dPhiSide_pred_raw, VBF_Mjj_presel_4JV_dPhiSide_sel, -3);
    DoRebinning(VBF_Jet1Pt_presel_4JV_dPhiSide_pred_raw, VBF_Jet1Pt_presel_4JV_dPhiSide_sel, 2);
    DoRebinning(VBF_Jet2Pt_presel_4JV_dPhiSide_pred_raw, VBF_Jet2Pt_presel_4JV_dPhiSide_sel, 2);
    DoRebinning(VBF_Jet3Pt_presel_4JV_dPhiSide_pred_raw, VBF_Jet3Pt_presel_4JV_dPhiSide_sel, 2);
    DoRebinning(VBF_Jet1Eta_presel_4JV_dPhiSide_pred_raw, VBF_Jet1Eta_presel_4JV_dPhiSide_sel, 5);
    DoRebinning(VBF_Jet2Eta_presel_4JV_dPhiSide_pred_raw, VBF_Jet2Eta_presel_4JV_dPhiSide_sel, 5);
    DoRebinning(VBF_Jet3Eta_presel_4JV_dPhiSide_pred_raw, VBF_Jet3Eta_presel_4JV_dPhiSide_sel, 5);
    DoRebinning(VBF_PTjj_presel_4JV_dPhiSide_pred_raw, VBF_PTjj_presel_4JV_dPhiSide_sel, -2);
    DoRebinning(VBF_MET_presel_4JV_dPhiSide_pred_raw, VBF_MET_presel_4JV_dPhiSide_sel, -2);
    DoRebinning(VBF_METsoft_presel_4JV_dPhiSide_pred_raw, VBF_METsoft_presel_4JV_dPhiSide_sel, 2);
    DoRebinning(VBF_METsig_presel_4JV_dPhiSide_pred_raw, VBF_METsig_presel_4JV_dPhiSide_sel, 2);
    DoRebinning(VBF_MHTsig_presel_4JV_dPhiSide_pred_raw, VBF_MHTsig_presel_4JV_dPhiSide_sel, 2);
    DoRebinning(VBF_minDeltaPhiPTj12_presel_4JV_dPhiSide_pred_raw, VBF_minDeltaPhiPTj12_presel_4JV_dPhiSide_sel, 5);
    DoRebinning(VBF_maxDeltaPhiPTj12_presel_4JV_dPhiSide_pred_raw, VBF_maxDeltaPhiPTj12_presel_4JV_dPhiSide_sel, 5);
    DoRebinning(VBF_DeltaPhiPTj3_presel_4JV_dPhiSide_pred_raw, VBF_DeltaPhiPTj3_presel_4JV_dPhiSide_sel, 5);

    DoRebinning(VBF_dPhi_dEta_pred_raw, VBF_dPhi_dEta_sel, 5);
    DoRebinning(VBF_dEta_dEta_pred_raw, VBF_dEta_dEta_sel, 5);
    DoRebinning(VBF_Mjj_dEta_pred_raw, VBF_Mjj_dEta_sel, -3);
    DoRebinning(VBF_Jet1Pt_dEta_pred_raw, VBF_Jet1Pt_dEta_sel, 2);
    DoRebinning(VBF_Jet2Pt_dEta_pred_raw, VBF_Jet2Pt_dEta_sel, 2);
    DoRebinning(VBF_Jet3Pt_dEta_pred_raw, VBF_Jet3Pt_dEta_sel, 2);
    DoRebinning(VBF_Jet1Eta_dEta_pred_raw, VBF_Jet1Eta_dEta_sel, 5);
    DoRebinning(VBF_Jet2Eta_dEta_pred_raw, VBF_Jet2Eta_dEta_sel, 5);
    DoRebinning(VBF_Jet3Eta_dEta_pred_raw, VBF_Jet3Eta_dEta_sel, 5);
    DoRebinning(VBF_PTjj_dEta_pred_raw, VBF_PTjj_dEta_sel, -2);
    DoRebinning(VBF_MET_dEta_pred_raw, VBF_MET_dEta_sel, -2);
    DoRebinning(VBF_METsoft_dEta_pred_raw, VBF_METsoft_dEta_sel, 2);
    DoRebinning(VBF_METsig_dEta_pred_raw, VBF_METsig_dEta_sel, 2);
    DoRebinning(VBF_MHTsig_dEta_pred_raw, VBF_MHTsig_dEta_sel, 2);
    DoRebinning(VBF_minDeltaPhiPTj12_dEta_pred_raw, VBF_minDeltaPhiPTj12_dEta_sel, 5);
    DoRebinning(VBF_maxDeltaPhiPTj12_dEta_pred_raw, VBF_maxDeltaPhiPTj12_dEta_sel, 5);
    DoRebinning(VBF_DeltaPhiPTj3_dEta_pred_raw, VBF_DeltaPhiPTj3_dEta_sel, 5);

    DoRebinning(VBF_dPhi_jj_pred_raw, VBF_dPhi_jj_sel, 5);
    DoRebinning(VBF_dEta_jj_pred_raw, VBF_dEta_jj_sel, 5);
    DoRebinning(VBF_Mjj_jj_pred_raw, VBF_Mjj_jj_sel, -3);
    DoRebinning(VBF_Jet1Pt_jj_pred_raw, VBF_Jet1Pt_jj_sel, 2);
    DoRebinning(VBF_Jet2Pt_jj_pred_raw, VBF_Jet2Pt_jj_sel, 2);
    DoRebinning(VBF_Jet3Pt_jj_pred_raw, VBF_Jet3Pt_jj_sel, 2);
    DoRebinning(VBF_Jet1Eta_jj_pred_raw, VBF_Jet1Eta_jj_sel, 5);
    DoRebinning(VBF_Jet2Eta_jj_pred_raw, VBF_Jet2Eta_jj_sel, 5);
    DoRebinning(VBF_Jet3Eta_jj_pred_raw, VBF_Jet3Eta_jj_sel, 5);
    DoRebinning(VBF_PTjj_jj_pred_raw, VBF_PTjj_jj_sel, -2);
    DoRebinning(VBF_MET_jj_pred_raw, VBF_MET_jj_sel, -2);
    DoRebinning(VBF_METsoft_jj_pred_raw, VBF_METsoft_jj_sel, 2);
    DoRebinning(VBF_METsig_jj_pred_raw, VBF_METsig_jj_sel, 2);
    DoRebinning(VBF_MHTsig_jj_pred_raw, VBF_MHTsig_jj_sel, 2);
    DoRebinning(VBF_minDeltaPhiPTj12_jj_pred_raw, VBF_minDeltaPhiPTj12_jj_sel, 5);
    DoRebinning(VBF_maxDeltaPhiPTj12_jj_pred_raw, VBF_maxDeltaPhiPTj12_jj_sel, 5);
    DoRebinning(VBF_DeltaPhiPTj3_jj_pred_raw, VBF_DeltaPhiPTj3_jj_sel, 5);

    DoRebinning(VBF_dPhi_dEta_3JV_pred_raw, VBF_dPhi_dEta_3JV_sel, 5);
    DoRebinning(VBF_dEta_dEta_3JV_pred_raw, VBF_dEta_dEta_3JV_sel, 5);
    DoRebinning(VBF_Mjj_dEta_3JV_pred_raw, VBF_Mjj_dEta_3JV_sel, -3);
    DoRebinning(VBF_Jet1Pt_dEta_3JV_pred_raw, VBF_Jet1Pt_dEta_3JV_sel, 2);
    DoRebinning(VBF_Jet2Pt_dEta_3JV_pred_raw, VBF_Jet2Pt_dEta_3JV_sel, 2);
    DoRebinning(VBF_Jet3Pt_dEta_3JV_pred_raw, VBF_Jet3Pt_dEta_3JV_sel, 2);
    DoRebinning(VBF_Jet1Eta_dEta_3JV_pred_raw, VBF_Jet1Eta_dEta_3JV_sel, 5);
    DoRebinning(VBF_Jet2Eta_dEta_3JV_pred_raw, VBF_Jet2Eta_dEta_3JV_sel, 5);
    DoRebinning(VBF_Jet3Eta_dEta_3JV_pred_raw, VBF_Jet3Eta_dEta_3JV_sel, 5);
    DoRebinning(VBF_PTjj_dEta_3JV_pred_raw, VBF_PTjj_dEta_3JV_sel, -2);
    DoRebinning(VBF_MET_dEta_3JV_pred_raw, VBF_MET_dEta_3JV_sel, -2);
    DoRebinning(VBF_METsoft_dEta_3JV_pred_raw, VBF_METsoft_dEta_3JV_sel, 2);
    DoRebinning(VBF_METsig_dEta_3JV_pred_raw, VBF_METsig_dEta_3JV_sel, 2);
    DoRebinning(VBF_MHTsig_dEta_3JV_pred_raw, VBF_MHTsig_dEta_3JV_sel, 2);
    DoRebinning(VBF_minDeltaPhiPTj12_dEta_3JV_pred_raw, VBF_minDeltaPhiPTj12_dEta_3JV_sel, 5);
    DoRebinning(VBF_maxDeltaPhiPTj12_dEta_3JV_pred_raw, VBF_maxDeltaPhiPTj12_dEta_3JV_sel, 5);
    DoRebinning(VBF_DeltaPhiPTj3_dEta_3JV_pred_raw, VBF_DeltaPhiPTj3_dEta_3JV_sel, 5);

    DoRebinning(VBF_dPhi_dEta_3JV_dPhiPTjj_pred_raw, VBF_dPhi_dEta_3JV_dPhiPTjj_sel, 5);
    DoRebinning(VBF_dEta_dEta_3JV_dPhiPTjj_pred_raw, VBF_dEta_dEta_3JV_dPhiPTjj_sel, 5);
    DoRebinning(VBF_Mjj_dEta_3JV_dPhiPTjj_pred_raw, VBF_Mjj_dEta_3JV_dPhiPTjj_sel, -3);
    DoRebinning(VBF_Jet1Pt_dEta_3JV_dPhiPTjj_pred_raw, VBF_Jet1Pt_dEta_3JV_dPhiPTjj_sel, 2);
    DoRebinning(VBF_Jet2Pt_dEta_3JV_dPhiPTjj_pred_raw, VBF_Jet2Pt_dEta_3JV_dPhiPTjj_sel, 2);
    DoRebinning(VBF_Jet3Pt_dEta_3JV_dPhiPTjj_pred_raw, VBF_Jet3Pt_dEta_3JV_dPhiPTjj_sel, 2);
    DoRebinning(VBF_Jet1Eta_dEta_3JV_dPhiPTjj_pred_raw, VBF_Jet1Eta_dEta_3JV_dPhiPTjj_sel, 5);
    DoRebinning(VBF_Jet2Eta_dEta_3JV_dPhiPTjj_pred_raw, VBF_Jet2Eta_dEta_3JV_dPhiPTjj_sel, 5);
    DoRebinning(VBF_Jet3Eta_dEta_3JV_dPhiPTjj_pred_raw, VBF_Jet3Eta_dEta_3JV_dPhiPTjj_sel, 5);
    DoRebinning(VBF_PTjj_dEta_3JV_dPhiPTjj_pred_raw, VBF_PTjj_dEta_3JV_dPhiPTjj_sel, -2);
    DoRebinning(VBF_MET_dEta_3JV_dPhiPTjj_pred_raw, VBF_MET_dEta_3JV_dPhiPTjj_sel, -2);
    DoRebinning(VBF_METsoft_dEta_3JV_dPhiPTjj_pred_raw, VBF_METsoft_dEta_3JV_dPhiPTjj_sel, 2);
    DoRebinning(VBF_METsig_dEta_3JV_dPhiPTjj_pred_raw, VBF_METsig_dEta_3JV_dPhiPTjj_sel, 2);
    DoRebinning(VBF_MHTsig_dEta_3JV_dPhiPTjj_pred_raw, VBF_MHTsig_dEta_3JV_dPhiPTjj_sel, 2);
    DoRebinning(VBF_minDeltaPhiPTj12_dEta_3JV_dPhiPTjj_pred_raw, VBF_minDeltaPhiPTj12_dEta_3JV_dPhiPTjj_sel, 5);
    DoRebinning(VBF_maxDeltaPhiPTj12_dEta_3JV_dPhiPTjj_pred_raw, VBF_maxDeltaPhiPTj12_dEta_3JV_dPhiPTjj_sel, 5);
    DoRebinning(VBF_DeltaPhiPTj3_dEta_3JV_dPhiPTjj_pred_raw, VBF_DeltaPhiPTj3_dEta_3JV_dPhiPTjj_sel, 5);

    cout << "after rebinning all histograms" << endl;
    //----------------------------------------------------------//

    //----------------------------------------------------------//
    // fill prediction histos
    HT_presel_pred = CalcPrediction(HT_presel_pred_raw);
    MHT_presel_pred = CalcPrediction(MHT_presel_pred_raw);
    MET_presel_pred = CalcPrediction(MET_presel_pred_raw);
    NJets_presel_pred = CalcPrediction(NJets_presel_pred_raw);
    NBJets_presel_pred = CalcPrediction(NBJets_presel_pred_raw);
    Jet1Pt_presel_pred = CalcPrediction(Jet1Pt_presel_pred_raw);
    Jet2Pt_presel_pred = CalcPrediction(Jet2Pt_presel_pred_raw);
    Jet3Pt_presel_pred = CalcPrediction(Jet3Pt_presel_pred_raw);
    Jet1Eta_presel_pred = CalcPrediction(Jet1Eta_presel_pred_raw);
    Jet2Eta_presel_pred = CalcPrediction(Jet2Eta_presel_pred_raw);
    Jet3Eta_presel_pred = CalcPrediction(Jet3Eta_presel_pred_raw);
    DeltaPhi1_presel_pred = CalcPrediction(DeltaPhi1_presel_pred_raw);
    DeltaPhi2_presel_pred = CalcPrediction(DeltaPhi2_presel_pred_raw);
    DeltaPhi3_presel_pred = CalcPrediction(DeltaPhi3_presel_pred_raw);

    HT_deltaPhi_pred = CalcPrediction(HT_deltaPhi_pred_raw);
    MHT_deltaPhi_pred = CalcPrediction(MHT_deltaPhi_pred_raw);
    MET_deltaPhi_pred = CalcPrediction(MET_deltaPhi_pred_raw);
    Jet1Pt_deltaPhi_pred = CalcPrediction(Jet1Pt_deltaPhi_pred_raw);
    Jet2Pt_deltaPhi_pred = CalcPrediction(Jet2Pt_deltaPhi_pred_raw);
    Jet3Pt_deltaPhi_pred = CalcPrediction(Jet3Pt_deltaPhi_pred_raw);
    Jet1Eta_deltaPhi_pred = CalcPrediction(Jet1Eta_deltaPhi_pred_raw);
    Jet2Eta_deltaPhi_pred = CalcPrediction(Jet2Eta_deltaPhi_pred_raw);
    Jet3Eta_deltaPhi_pred = CalcPrediction(Jet3Eta_deltaPhi_pred_raw);

    MHT_JetBin1_HTinclusive_pred = CalcPrediction(MHT_JetBin1_HTinclusive_pred_raw);
    MHT_JetBin2_HTinclusive_pred = CalcPrediction(MHT_JetBin2_HTinclusive_pred_raw);
    MHT_JetBin3_HTinclusive_pred = CalcPrediction(MHT_JetBin3_HTinclusive_pred_raw);
    MHT_JetBin4_HTinclusive_pred = CalcPrediction(MHT_JetBin4_HTinclusive_pred_raw);

    MET_JetBin1_HTinclusive_pred = CalcPrediction(MET_JetBin1_HTinclusive_pred_raw);
    MET_JetBin2_HTinclusive_pred = CalcPrediction(MET_JetBin2_HTinclusive_pred_raw);
    MET_JetBin3_HTinclusive_pred = CalcPrediction(MET_JetBin3_HTinclusive_pred_raw);
    MET_JetBin4_HTinclusive_pred = CalcPrediction(MET_JetBin4_HTinclusive_pred_raw);

    HT_baseline_pred = CalcPrediction( HT_baseline_pred_raw);
    MHT_baseline_pred = CalcPrediction( MHT_baseline_pred_raw);
    MET_baseline_pred = CalcPrediction( MET_baseline_pred_raw);

    Jet1Pt_JetBin1_baseline_pred = CalcPrediction(Jet1Pt_JetBin1_baseline_pred_raw);
    Jet2Pt_JetBin1_baseline_pred = CalcPrediction(Jet2Pt_JetBin1_baseline_pred_raw);
    Jet1Eta_JetBin1_baseline_pred = CalcPrediction(Jet1Eta_JetBin1_baseline_pred_raw);
    Jet2Eta_JetBin1_baseline_pred = CalcPrediction(Jet2Eta_JetBin1_baseline_pred_raw);
    DeltaPhi1_JetBin1_baseline_pred = CalcPrediction(DeltaPhi1_JetBin1_baseline_pred_raw);
    DeltaPhi2_JetBin1_baseline_pred = CalcPrediction(DeltaPhi2_JetBin1_baseline_pred_raw);

    Jet1Pt_JetBin2_baseline_pred = CalcPrediction(Jet1Pt_JetBin2_baseline_pred_raw);
    Jet2Pt_JetBin2_baseline_pred = CalcPrediction(Jet2Pt_JetBin2_baseline_pred_raw);
    Jet3Pt_JetBin2_baseline_pred = CalcPrediction(Jet3Pt_JetBin2_baseline_pred_raw);
    Jet1Eta_JetBin2_baseline_pred = CalcPrediction(Jet1Eta_JetBin2_baseline_pred_raw);
    Jet2Eta_JetBin2_baseline_pred = CalcPrediction(Jet2Eta_JetBin2_baseline_pred_raw);
    Jet3Eta_JetBin2_baseline_pred = CalcPrediction(Jet3Eta_JetBin2_baseline_pred_raw);
    DeltaPhi1_JetBin2_baseline_pred = CalcPrediction(DeltaPhi1_JetBin2_baseline_pred_raw);
    DeltaPhi2_JetBin2_baseline_pred = CalcPrediction(DeltaPhi2_JetBin2_baseline_pred_raw);
    DeltaPhi3_JetBin2_baseline_pred = CalcPrediction(DeltaPhi3_JetBin2_baseline_pred_raw);

    Jet1Pt_JetBin3_baseline_pred = CalcPrediction(Jet1Pt_JetBin3_baseline_pred_raw);
    Jet2Pt_JetBin3_baseline_pred = CalcPrediction(Jet2Pt_JetBin3_baseline_pred_raw);
    Jet3Pt_JetBin3_baseline_pred = CalcPrediction(Jet3Pt_JetBin3_baseline_pred_raw);
    Jet1Eta_JetBin3_baseline_pred = CalcPrediction(Jet1Eta_JetBin3_baseline_pred_raw);
    Jet2Eta_JetBin3_baseline_pred = CalcPrediction(Jet2Eta_JetBin3_baseline_pred_raw);
    Jet3Eta_JetBin3_baseline_pred = CalcPrediction(Jet3Eta_JetBin3_baseline_pred_raw);
    DeltaPhi1_JetBin3_baseline_pred = CalcPrediction(DeltaPhi1_JetBin3_baseline_pred_raw);
    DeltaPhi2_JetBin3_baseline_pred = CalcPrediction(DeltaPhi2_JetBin3_baseline_pred_raw);
    DeltaPhi3_JetBin3_baseline_pred = CalcPrediction(DeltaPhi3_JetBin3_baseline_pred_raw);

    Jet1Pt_JetBin4_baseline_pred = CalcPrediction(Jet1Pt_JetBin4_baseline_pred_raw);
    Jet2Pt_JetBin4_baseline_pred = CalcPrediction(Jet2Pt_JetBin4_baseline_pred_raw);
    Jet3Pt_JetBin4_baseline_pred = CalcPrediction(Jet3Pt_JetBin4_baseline_pred_raw);
    Jet1Eta_JetBin4_baseline_pred = CalcPrediction(Jet1Eta_JetBin4_baseline_pred_raw);
    Jet2Eta_JetBin4_baseline_pred = CalcPrediction(Jet2Eta_JetBin4_baseline_pred_raw);
    Jet3Eta_JetBin4_baseline_pred = CalcPrediction(Jet3Eta_JetBin4_baseline_pred_raw);
    DeltaPhi1_JetBin4_baseline_pred = CalcPrediction(DeltaPhi1_JetBin4_baseline_pred_raw);
    DeltaPhi2_JetBin4_baseline_pred = CalcPrediction(DeltaPhi2_JetBin4_baseline_pred_raw);
    DeltaPhi3_JetBin4_baseline_pred = CalcPrediction(DeltaPhi3_JetBin4_baseline_pred_raw);

    HT_JetBin1_baseline_withoutDeltaPhi_pred = CalcPrediction(HT_JetBin1_baseline_withoutDeltaPhi_pred_raw);
    MHT_JetBin1_baseline_withoutDeltaPhi_pred = CalcPrediction(MHT_JetBin1_baseline_withoutDeltaPhi_pred_raw);
    MET_JetBin1_baseline_withoutDeltaPhi_pred = CalcPrediction(MET_JetBin1_baseline_withoutDeltaPhi_pred_raw);
    Jet1Pt_JetBin1_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet1Pt_JetBin1_baseline_withoutDeltaPhi_pred_raw);
    Jet2Pt_JetBin1_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet2Pt_JetBin1_baseline_withoutDeltaPhi_pred_raw);
    Jet1Eta_JetBin1_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet1Eta_JetBin1_baseline_withoutDeltaPhi_pred_raw);
    Jet2Eta_JetBin1_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet2Eta_JetBin1_baseline_withoutDeltaPhi_pred_raw);
    DeltaPhi1_JetBin1_baseline_withoutDeltaPhi_pred = CalcPrediction(DeltaPhi1_JetBin1_baseline_withoutDeltaPhi_pred_raw);
    DeltaPhi2_JetBin1_baseline_withoutDeltaPhi_pred = CalcPrediction(DeltaPhi2_JetBin1_baseline_withoutDeltaPhi_pred_raw);

    HT_JetBin2_baseline_withoutDeltaPhi_pred = CalcPrediction(HT_JetBin2_baseline_withoutDeltaPhi_pred_raw);
    MHT_JetBin2_baseline_withoutDeltaPhi_pred = CalcPrediction(MHT_JetBin2_baseline_withoutDeltaPhi_pred_raw);
    MET_JetBin2_baseline_withoutDeltaPhi_pred = CalcPrediction(MET_JetBin2_baseline_withoutDeltaPhi_pred_raw);
    Jet1Pt_JetBin2_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet1Pt_JetBin2_baseline_withoutDeltaPhi_pred_raw);
    Jet2Pt_JetBin2_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet2Pt_JetBin2_baseline_withoutDeltaPhi_pred_raw);
    Jet3Pt_JetBin2_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet3Pt_JetBin2_baseline_withoutDeltaPhi_pred_raw);
    Jet1Eta_JetBin2_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet1Eta_JetBin2_baseline_withoutDeltaPhi_pred_raw);
    Jet2Eta_JetBin2_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet2Eta_JetBin2_baseline_withoutDeltaPhi_pred_raw);
    Jet3Eta_JetBin2_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet3Eta_JetBin2_baseline_withoutDeltaPhi_pred_raw);
    DeltaPhi1_JetBin2_baseline_withoutDeltaPhi_pred = CalcPrediction(DeltaPhi1_JetBin2_baseline_withoutDeltaPhi_pred_raw);
    DeltaPhi2_JetBin2_baseline_withoutDeltaPhi_pred = CalcPrediction(DeltaPhi2_JetBin2_baseline_withoutDeltaPhi_pred_raw);
    DeltaPhi3_JetBin2_baseline_withoutDeltaPhi_pred = CalcPrediction(DeltaPhi3_JetBin2_baseline_withoutDeltaPhi_pred_raw);

    HT_JetBin3_baseline_withoutDeltaPhi_pred = CalcPrediction(HT_JetBin3_baseline_withoutDeltaPhi_pred_raw);
    MHT_JetBin3_baseline_withoutDeltaPhi_pred = CalcPrediction(MHT_JetBin3_baseline_withoutDeltaPhi_pred_raw);
    MET_JetBin3_baseline_withoutDeltaPhi_pred = CalcPrediction(MET_JetBin3_baseline_withoutDeltaPhi_pred_raw);
    Jet1Pt_JetBin3_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet1Pt_JetBin3_baseline_withoutDeltaPhi_pred_raw);
    Jet2Pt_JetBin3_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet2Pt_JetBin3_baseline_withoutDeltaPhi_pred_raw);
    Jet3Pt_JetBin3_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet3Pt_JetBin3_baseline_withoutDeltaPhi_pred_raw);
    Jet1Eta_JetBin3_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet1Eta_JetBin3_baseline_withoutDeltaPhi_pred_raw);
    Jet2Eta_JetBin3_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet2Eta_JetBin3_baseline_withoutDeltaPhi_pred_raw);
    Jet3Eta_JetBin3_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet3Eta_JetBin3_baseline_withoutDeltaPhi_pred_raw);
    DeltaPhi1_JetBin3_baseline_withoutDeltaPhi_pred = CalcPrediction(DeltaPhi1_JetBin3_baseline_withoutDeltaPhi_pred_raw);
    DeltaPhi2_JetBin3_baseline_withoutDeltaPhi_pred = CalcPrediction(DeltaPhi2_JetBin3_baseline_withoutDeltaPhi_pred_raw);
    DeltaPhi3_JetBin3_baseline_withoutDeltaPhi_pred = CalcPrediction(DeltaPhi3_JetBin3_baseline_withoutDeltaPhi_pred_raw);

    HT_JetBin4_baseline_withoutDeltaPhi_pred = CalcPrediction(HT_JetBin4_baseline_withoutDeltaPhi_pred_raw);
    MHT_JetBin4_baseline_withoutDeltaPhi_pred = CalcPrediction(MHT_JetBin4_baseline_withoutDeltaPhi_pred_raw);
    MET_JetBin4_baseline_withoutDeltaPhi_pred = CalcPrediction(MET_JetBin4_baseline_withoutDeltaPhi_pred_raw);
    Jet1Pt_JetBin4_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet1Pt_JetBin4_baseline_withoutDeltaPhi_pred_raw);
    Jet2Pt_JetBin4_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet2Pt_JetBin4_baseline_withoutDeltaPhi_pred_raw);
    Jet3Pt_JetBin4_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet3Pt_JetBin4_baseline_withoutDeltaPhi_pred_raw);
    Jet1Eta_JetBin4_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet1Eta_JetBin4_baseline_withoutDeltaPhi_pred_raw);
    Jet2Eta_JetBin4_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet2Eta_JetBin4_baseline_withoutDeltaPhi_pred_raw);
    Jet3Eta_JetBin4_baseline_withoutDeltaPhi_pred = CalcPrediction(Jet3Eta_JetBin4_baseline_withoutDeltaPhi_pred_raw);
    DeltaPhi1_JetBin4_baseline_withoutDeltaPhi_pred = CalcPrediction(DeltaPhi1_JetBin4_baseline_withoutDeltaPhi_pred_raw);
    DeltaPhi2_JetBin4_baseline_withoutDeltaPhi_pred = CalcPrediction(DeltaPhi2_JetBin4_baseline_withoutDeltaPhi_pred_raw);
    DeltaPhi3_JetBin4_baseline_withoutDeltaPhi_pred = CalcPrediction(DeltaPhi3_JetBin4_baseline_withoutDeltaPhi_pred_raw);

    NJets_baseline_withoutMET_pred = CalcPrediction( NJets_baseline_withoutMET_pred_raw);
    NJets_baseline_pred = CalcPrediction( NJets_baseline_pred_raw);
    NJets_baseline_withoutDeltaPhi_withoutMET_pred = CalcPrediction(NJets_baseline_withoutDeltaPhi_withoutMET_pred_raw);
    NJets_baseline_withoutDeltaPhi_pred = CalcPrediction(NJets_baseline_withoutDeltaPhi_pred_raw);

    NBJets_baseline_withoutMET_pred = CalcPrediction( NBJets_baseline_withoutMET_pred_raw);
    NBJets_baseline_pred = CalcPrediction( NBJets_baseline_pred_raw);
    NBJets_baseline_withoutDeltaPhi_withoutMET_pred = CalcPrediction(NBJets_baseline_withoutDeltaPhi_withoutMET_pred_raw);
    NBJets_baseline_withoutDeltaPhi_pred = CalcPrediction(NBJets_baseline_withoutDeltaPhi_pred_raw);

    VBF_dPhi_presel_pred = CalcPrediction( VBF_dPhi_presel_pred_raw);
    VBF_dEta_presel_pred = CalcPrediction( VBF_dEta_presel_pred_raw);
    VBF_Mjj_presel_pred = CalcPrediction( VBF_Mjj_presel_pred_raw);
    VBF_Jet1Pt_presel_pred = CalcPrediction( VBF_Jet1Pt_presel_pred_raw);
    VBF_Jet2Pt_presel_pred = CalcPrediction( VBF_Jet2Pt_presel_pred_raw);
    VBF_Jet3Pt_presel_pred = CalcPrediction( VBF_Jet3Pt_presel_pred_raw);
    VBF_Jet1Eta_presel_pred = CalcPrediction( VBF_Jet1Eta_presel_pred_raw);
    VBF_Jet2Eta_presel_pred = CalcPrediction( VBF_Jet2Eta_presel_pred_raw);
    VBF_Jet3Eta_presel_pred = CalcPrediction( VBF_Jet3Eta_presel_pred_raw);
    VBF_PTjj_presel_pred = CalcPrediction( VBF_PTjj_presel_pred_raw);
    VBF_MET_presel_pred = CalcPrediction( VBF_MET_presel_pred_raw);
    VBF_METsoft_presel_pred = CalcPrediction( VBF_METsoft_presel_pred_raw);
    VBF_METsig_presel_pred = CalcPrediction( VBF_METsig_presel_pred_raw);
    VBF_MHTsig_presel_pred = CalcPrediction( VBF_MHTsig_presel_pred_raw);
    VBF_minDeltaPhiPTj12_presel_pred = CalcPrediction( VBF_minDeltaPhiPTj12_presel_pred_raw);
    VBF_maxDeltaPhiPTj12_presel_pred = CalcPrediction( VBF_maxDeltaPhiPTj12_presel_pred_raw);
    VBF_DeltaPhiPTj3_presel_pred = CalcPrediction( VBF_DeltaPhiPTj3_presel_pred_raw);

    VBF_dPhi_presel_4JV_dPhiSide_pred = CalcPrediction( VBF_dPhi_presel_4JV_dPhiSide_pred_raw);
    VBF_dEta_presel_4JV_dPhiSide_pred = CalcPrediction( VBF_dEta_presel_4JV_dPhiSide_pred_raw);
    VBF_Mjj_presel_4JV_dPhiSide_pred = CalcPrediction( VBF_Mjj_presel_4JV_dPhiSide_pred_raw);
    VBF_Jet1Pt_presel_4JV_dPhiSide_pred = CalcPrediction( VBF_Jet1Pt_presel_4JV_dPhiSide_pred_raw);
    VBF_Jet2Pt_presel_4JV_dPhiSide_pred = CalcPrediction( VBF_Jet2Pt_presel_4JV_dPhiSide_pred_raw);
    VBF_Jet3Pt_presel_4JV_dPhiSide_pred = CalcPrediction( VBF_Jet3Pt_presel_4JV_dPhiSide_pred_raw);
    VBF_Jet1Eta_presel_4JV_dPhiSide_pred = CalcPrediction( VBF_Jet1Eta_presel_4JV_dPhiSide_pred_raw);
    VBF_Jet2Eta_presel_4JV_dPhiSide_pred = CalcPrediction( VBF_Jet2Eta_presel_4JV_dPhiSide_pred_raw);
    VBF_Jet3Eta_presel_4JV_dPhiSide_pred = CalcPrediction( VBF_Jet3Eta_presel_4JV_dPhiSide_pred_raw);
    VBF_PTjj_presel_4JV_dPhiSide_pred = CalcPrediction( VBF_PTjj_presel_4JV_dPhiSide_pred_raw);
    VBF_MET_presel_4JV_dPhiSide_pred = CalcPrediction( VBF_MET_presel_4JV_dPhiSide_pred_raw);
    VBF_METsoft_presel_4JV_dPhiSide_pred = CalcPrediction( VBF_METsoft_presel_4JV_dPhiSide_pred_raw);
    VBF_METsig_presel_4JV_dPhiSide_pred = CalcPrediction( VBF_METsig_presel_4JV_dPhiSide_pred_raw);
    VBF_MHTsig_presel_4JV_dPhiSide_pred = CalcPrediction( VBF_MHTsig_presel_4JV_dPhiSide_pred_raw);
    VBF_minDeltaPhiPTj12_presel_4JV_dPhiSide_pred = CalcPrediction( VBF_minDeltaPhiPTj12_presel_4JV_dPhiSide_pred_raw);
    VBF_maxDeltaPhiPTj12_presel_4JV_dPhiSide_pred = CalcPrediction( VBF_maxDeltaPhiPTj12_presel_4JV_dPhiSide_pred_raw);
    VBF_DeltaPhiPTj3_presel_4JV_dPhiSide_pred = CalcPrediction( VBF_DeltaPhiPTj3_presel_4JV_dPhiSide_pred_raw);

    VBF_dPhi_dEta_pred = CalcPrediction( VBF_dPhi_dEta_pred_raw);
    VBF_dEta_dEta_pred = CalcPrediction( VBF_dEta_dEta_pred_raw);
    VBF_Mjj_dEta_pred = CalcPrediction( VBF_Mjj_dEta_pred_raw);
    VBF_Jet1Pt_dEta_pred = CalcPrediction( VBF_Jet1Pt_dEta_pred_raw);
    VBF_Jet2Pt_dEta_pred = CalcPrediction( VBF_Jet2Pt_dEta_pred_raw);
    VBF_Jet3Pt_dEta_pred = CalcPrediction( VBF_Jet3Pt_dEta_pred_raw);
    VBF_Jet1Eta_dEta_pred = CalcPrediction( VBF_Jet1Eta_dEta_pred_raw);
    VBF_Jet2Eta_dEta_pred = CalcPrediction( VBF_Jet2Eta_dEta_pred_raw);
    VBF_Jet3Eta_dEta_pred = CalcPrediction( VBF_Jet3Eta_dEta_pred_raw);
    VBF_PTjj_dEta_pred = CalcPrediction( VBF_PTjj_dEta_pred_raw);
    VBF_MET_dEta_pred = CalcPrediction( VBF_MET_dEta_pred_raw);
    VBF_METsoft_dEta_pred = CalcPrediction( VBF_METsoft_dEta_pred_raw);
    VBF_METsig_dEta_pred = CalcPrediction( VBF_METsig_dEta_pred_raw);
    VBF_MHTsig_dEta_pred = CalcPrediction( VBF_MHTsig_dEta_pred_raw);
    VBF_minDeltaPhiPTj12_dEta_pred = CalcPrediction( VBF_minDeltaPhiPTj12_dEta_pred_raw);
    VBF_maxDeltaPhiPTj12_dEta_pred = CalcPrediction( VBF_maxDeltaPhiPTj12_dEta_pred_raw);
    VBF_DeltaPhiPTj3_dEta_pred = CalcPrediction( VBF_DeltaPhiPTj3_dEta_pred_raw);

    VBF_dPhi_jj_pred = CalcPrediction( VBF_dPhi_jj_pred_raw);
    VBF_dEta_jj_pred = CalcPrediction( VBF_dEta_jj_pred_raw);
    VBF_Mjj_jj_pred = CalcPrediction( VBF_Mjj_jj_pred_raw);
    VBF_Jet1Pt_jj_pred = CalcPrediction( VBF_Jet1Pt_jj_pred_raw);
    VBF_Jet2Pt_jj_pred = CalcPrediction( VBF_Jet2Pt_jj_pred_raw);
    VBF_Jet3Pt_jj_pred = CalcPrediction( VBF_Jet3Pt_jj_pred_raw);
    VBF_Jet1Eta_jj_pred = CalcPrediction( VBF_Jet1Eta_jj_pred_raw);
    VBF_Jet2Eta_jj_pred = CalcPrediction( VBF_Jet2Eta_jj_pred_raw);
    VBF_Jet3Eta_jj_pred = CalcPrediction( VBF_Jet3Eta_jj_pred_raw);
    VBF_PTjj_jj_pred = CalcPrediction( VBF_PTjj_jj_pred_raw);
    VBF_MET_jj_pred = CalcPrediction( VBF_MET_jj_pred_raw);
    VBF_METsoft_jj_pred = CalcPrediction( VBF_METsoft_jj_pred_raw);
    VBF_METsig_jj_pred = CalcPrediction( VBF_METsig_jj_pred_raw);
    VBF_MHTsig_jj_pred = CalcPrediction( VBF_MHTsig_jj_pred_raw);
    VBF_minDeltaPhiPTj12_jj_pred = CalcPrediction( VBF_minDeltaPhiPTj12_jj_pred_raw);
    VBF_maxDeltaPhiPTj12_jj_pred = CalcPrediction( VBF_maxDeltaPhiPTj12_jj_pred_raw);
    VBF_DeltaPhiPTj3_jj_pred = CalcPrediction( VBF_DeltaPhiPTj3_jj_pred_raw);

    VBF_dPhi_dEta_3JV_pred = CalcPrediction( VBF_dPhi_dEta_3JV_pred_raw);
    VBF_dEta_dEta_3JV_pred = CalcPrediction( VBF_dEta_dEta_3JV_pred_raw);
    VBF_Mjj_dEta_3JV_pred = CalcPrediction( VBF_Mjj_dEta_3JV_pred_raw);
    VBF_Jet1Pt_dEta_3JV_pred = CalcPrediction( VBF_Jet1Pt_dEta_3JV_pred_raw);
    VBF_Jet2Pt_dEta_3JV_pred = CalcPrediction( VBF_Jet2Pt_dEta_3JV_pred_raw);
    VBF_Jet3Pt_dEta_3JV_pred = CalcPrediction( VBF_Jet3Pt_dEta_3JV_pred_raw);
    VBF_Jet1Eta_dEta_3JV_pred = CalcPrediction( VBF_Jet1Eta_dEta_3JV_pred_raw);
    VBF_Jet2Eta_dEta_3JV_pred = CalcPrediction( VBF_Jet2Eta_dEta_3JV_pred_raw);
    VBF_Jet3Eta_dEta_3JV_pred = CalcPrediction( VBF_Jet3Eta_dEta_3JV_pred_raw);
    VBF_PTjj_dEta_3JV_pred = CalcPrediction( VBF_PTjj_dEta_3JV_pred_raw);
    VBF_MET_dEta_3JV_pred = CalcPrediction( VBF_MET_dEta_3JV_pred_raw);
    VBF_METsoft_dEta_3JV_pred = CalcPrediction( VBF_METsoft_dEta_3JV_pred_raw);
    VBF_METsig_dEta_3JV_pred = CalcPrediction( VBF_METsig_dEta_3JV_pred_raw);
    VBF_MHTsig_dEta_3JV_pred = CalcPrediction( VBF_MHTsig_dEta_3JV_pred_raw);
    VBF_minDeltaPhiPTj12_dEta_3JV_pred = CalcPrediction( VBF_minDeltaPhiPTj12_dEta_3JV_pred_raw);
    VBF_maxDeltaPhiPTj12_dEta_3JV_pred = CalcPrediction( VBF_maxDeltaPhiPTj12_dEta_3JV_pred_raw);
    VBF_DeltaPhiPTj3_dEta_3JV_pred = CalcPrediction( VBF_DeltaPhiPTj3_dEta_3JV_pred_raw);

    VBF_dPhi_dEta_3JV_dPhiPTjj_pred = CalcPrediction( VBF_dPhi_dEta_3JV_dPhiPTjj_pred_raw);
    VBF_dEta_dEta_3JV_dPhiPTjj_pred = CalcPrediction( VBF_dEta_dEta_3JV_dPhiPTjj_pred_raw);
    VBF_Mjj_dEta_3JV_dPhiPTjj_pred = CalcPrediction( VBF_Mjj_dEta_3JV_dPhiPTjj_pred_raw);
    VBF_Jet1Pt_dEta_3JV_dPhiPTjj_pred = CalcPrediction( VBF_Jet1Pt_dEta_3JV_dPhiPTjj_pred_raw);
    VBF_Jet2Pt_dEta_3JV_dPhiPTjj_pred = CalcPrediction( VBF_Jet2Pt_dEta_3JV_dPhiPTjj_pred_raw);
    VBF_Jet3Pt_dEta_3JV_dPhiPTjj_pred = CalcPrediction( VBF_Jet3Pt_dEta_3JV_dPhiPTjj_pred_raw);
    VBF_Jet1Eta_dEta_3JV_dPhiPTjj_pred = CalcPrediction( VBF_Jet1Eta_dEta_3JV_dPhiPTjj_pred_raw);
    VBF_Jet2Eta_dEta_3JV_dPhiPTjj_pred = CalcPrediction( VBF_Jet2Eta_dEta_3JV_dPhiPTjj_pred_raw);
    VBF_Jet3Eta_dEta_3JV_dPhiPTjj_pred = CalcPrediction( VBF_Jet3Eta_dEta_3JV_dPhiPTjj_pred_raw);
    VBF_PTjj_dEta_3JV_dPhiPTjj_pred = CalcPrediction( VBF_PTjj_dEta_3JV_dPhiPTjj_pred_raw);
    VBF_MET_dEta_3JV_dPhiPTjj_pred = CalcPrediction( VBF_MET_dEta_3JV_dPhiPTjj_pred_raw);
    VBF_METsoft_dEta_3JV_dPhiPTjj_pred = CalcPrediction( VBF_METsoft_dEta_3JV_dPhiPTjj_pred_raw);
    VBF_METsig_dEta_3JV_dPhiPTjj_pred = CalcPrediction( VBF_METsig_dEta_3JV_dPhiPTjj_pred_raw);
    VBF_MHTsig_dEta_3JV_dPhiPTjj_pred = CalcPrediction( VBF_MHTsig_dEta_3JV_dPhiPTjj_pred_raw);
    VBF_minDeltaPhiPTj12_dEta_3JV_dPhiPTjj_pred = CalcPrediction( VBF_minDeltaPhiPTj12_dEta_3JV_dPhiPTjj_pred_raw);
    VBF_maxDeltaPhiPTj12_dEta_3JV_dPhiPTjj_pred = CalcPrediction( VBF_maxDeltaPhiPTj12_dEta_3JV_dPhiPTjj_pred_raw);
    VBF_DeltaPhiPTj3_dEta_3JV_dPhiPTjj_pred = CalcPrediction( VBF_DeltaPhiPTj3_dEta_3JV_dPhiPTjj_pred_raw);

    cout << "after calculation of prediction" << endl;
    //----------------------------------------------------------//

    //----------------------------------------------------------//
    // write histos to file
    TFile* prediction_histos = new TFile("output_GetPrediction/prediction_histos" + postfix + ".root", "RECREATE");

    // Save histograms: prediction

    if (HTMHT) {
        HT_presel_pred->Write();
        MHT_presel_pred->Write();
        MET_presel_pred->Write();
        NJets_presel_pred->Write();
        NBJets_presel_pred->Write();
        Jet1Pt_presel_pred->Write();
        Jet2Pt_presel_pred->Write();
        Jet3Pt_presel_pred->Write();
        Jet1Eta_presel_pred->Write();
        Jet2Eta_presel_pred->Write();
        Jet3Eta_presel_pred->Write();

        HT_deltaPhi_pred->Write();
        MHT_deltaPhi_pred->Write();
        MET_deltaPhi_pred->Write();
        Jet1Pt_deltaPhi_pred->Write();
        Jet2Pt_deltaPhi_pred->Write();
        Jet3Pt_deltaPhi_pred->Write();
        Jet1Eta_deltaPhi_pred->Write();
        Jet2Eta_deltaPhi_pred->Write();
        Jet3Eta_deltaPhi_pred->Write();

        NJets_baseline_withoutMET_pred->Write();
        NJets_baseline_pred->Write();
        NJets_baseline_withoutDeltaPhi_withoutMET_pred->Write();
        NJets_baseline_withoutDeltaPhi_pred ->Write();

        NBJets_baseline_withoutMET_pred->Write();
        NBJets_baseline_pred->Write();
        NBJets_baseline_withoutDeltaPhi_withoutMET_pred->Write();
        NBJets_baseline_withoutDeltaPhi_pred ->Write();

        HT_baseline_pred->Write();
        MHT_baseline_pred->Write();
        MET_baseline_pred->Write();

        MHT_JetBin1_HTinclusive_pred->Write();
        MHT_JetBin2_HTinclusive_pred->Write();
        MHT_JetBin3_HTinclusive_pred->Write();
        MHT_JetBin4_HTinclusive_pred->Write();

        MET_JetBin1_HTinclusive_pred->Write();
        MET_JetBin2_HTinclusive_pred->Write();
        MET_JetBin3_HTinclusive_pred->Write();
        MET_JetBin4_HTinclusive_pred->Write();

        Jet1Pt_JetBin1_baseline_pred->Write();
        Jet2Pt_JetBin1_baseline_pred->Write();
        Jet1Eta_JetBin1_baseline_pred->Write();
        Jet2Eta_JetBin1_baseline_pred->Write();
        DeltaPhi1_JetBin1_baseline_pred->Write();
        DeltaPhi2_JetBin1_baseline_pred->Write();

        Jet1Pt_JetBin2_baseline_pred->Write();
        Jet2Pt_JetBin2_baseline_pred->Write();
        Jet3Pt_JetBin2_baseline_pred->Write();
        Jet1Eta_JetBin2_baseline_pred->Write();
        Jet2Eta_JetBin2_baseline_pred->Write();
        Jet3Eta_JetBin2_baseline_pred->Write();
        DeltaPhi1_JetBin2_baseline_pred->Write();
        DeltaPhi2_JetBin2_baseline_pred->Write();
        DeltaPhi3_JetBin2_baseline_pred->Write();

        Jet1Pt_JetBin3_baseline_pred->Write();
        Jet2Pt_JetBin3_baseline_pred->Write();
        Jet3Pt_JetBin3_baseline_pred->Write();
        Jet1Eta_JetBin3_baseline_pred->Write();
        Jet2Eta_JetBin3_baseline_pred->Write();
        Jet3Eta_JetBin3_baseline_pred->Write();
        DeltaPhi1_JetBin3_baseline_pred->Write();
        DeltaPhi2_JetBin3_baseline_pred->Write();
        DeltaPhi3_JetBin3_baseline_pred->Write();

        Jet1Pt_JetBin4_baseline_pred->Write();
        Jet2Pt_JetBin4_baseline_pred->Write();
        Jet3Pt_JetBin4_baseline_pred->Write();
        Jet1Eta_JetBin4_baseline_pred->Write();
        Jet2Eta_JetBin4_baseline_pred->Write();
        Jet3Eta_JetBin4_baseline_pred->Write();
        DeltaPhi1_JetBin4_baseline_pred->Write();
        DeltaPhi2_JetBin4_baseline_pred->Write();
        DeltaPhi3_JetBin4_baseline_pred->Write();

        HT_JetBin1_baseline_withoutDeltaPhi_pred->Write();
        MHT_JetBin1_baseline_withoutDeltaPhi_pred->Write();
        MET_JetBin1_baseline_withoutDeltaPhi_pred->Write();
        Jet1Pt_JetBin1_baseline_withoutDeltaPhi_pred->Write();
        Jet2Pt_JetBin1_baseline_withoutDeltaPhi_pred->Write();
        Jet1Eta_JetBin1_baseline_withoutDeltaPhi_pred->Write();
        Jet2Eta_JetBin1_baseline_withoutDeltaPhi_pred->Write();
        DeltaPhi1_JetBin1_baseline_withoutDeltaPhi_pred->Write();
        DeltaPhi2_JetBin1_baseline_withoutDeltaPhi_pred->Write();

        HT_JetBin2_baseline_withoutDeltaPhi_pred->Write();
        MHT_JetBin2_baseline_withoutDeltaPhi_pred->Write();
        MET_JetBin2_baseline_withoutDeltaPhi_pred->Write();
        Jet1Pt_JetBin2_baseline_withoutDeltaPhi_pred->Write();
        Jet2Pt_JetBin2_baseline_withoutDeltaPhi_pred->Write();
        Jet3Pt_JetBin2_baseline_withoutDeltaPhi_pred->Write();
        Jet1Eta_JetBin2_baseline_withoutDeltaPhi_pred->Write();
        Jet2Eta_JetBin2_baseline_withoutDeltaPhi_pred->Write();
        Jet3Eta_JetBin2_baseline_withoutDeltaPhi_pred->Write();
        DeltaPhi1_JetBin2_baseline_withoutDeltaPhi_pred->Write();
        DeltaPhi2_JetBin2_baseline_withoutDeltaPhi_pred->Write();
        DeltaPhi3_JetBin2_baseline_withoutDeltaPhi_pred->Write();

        HT_JetBin3_baseline_withoutDeltaPhi_pred->Write();
        MHT_JetBin3_baseline_withoutDeltaPhi_pred->Write();
        MET_JetBin3_baseline_withoutDeltaPhi_pred->Write();
        Jet1Pt_JetBin3_baseline_withoutDeltaPhi_pred->Write();
        Jet2Pt_JetBin3_baseline_withoutDeltaPhi_pred->Write();
        Jet3Pt_JetBin3_baseline_withoutDeltaPhi_pred->Write();
        Jet1Eta_JetBin3_baseline_withoutDeltaPhi_pred->Write();
        Jet2Eta_JetBin3_baseline_withoutDeltaPhi_pred->Write();
        Jet3Eta_JetBin3_baseline_withoutDeltaPhi_pred->Write();
        DeltaPhi1_JetBin3_baseline_withoutDeltaPhi_pred->Write();
        DeltaPhi2_JetBin3_baseline_withoutDeltaPhi_pred->Write();
        DeltaPhi3_JetBin3_baseline_withoutDeltaPhi_pred->Write();

        HT_JetBin4_baseline_withoutDeltaPhi_pred->Write();
        MHT_JetBin4_baseline_withoutDeltaPhi_pred->Write();
        MET_JetBin4_baseline_withoutDeltaPhi_pred->Write();
        Jet1Pt_JetBin4_baseline_withoutDeltaPhi_pred->Write();
        Jet2Pt_JetBin4_baseline_withoutDeltaPhi_pred->Write();
        Jet3Pt_JetBin4_baseline_withoutDeltaPhi_pred->Write();
        Jet1Eta_JetBin4_baseline_withoutDeltaPhi_pred->Write();
        Jet2Eta_JetBin4_baseline_withoutDeltaPhi_pred->Write();
        Jet3Eta_JetBin4_baseline_withoutDeltaPhi_pred->Write();
        DeltaPhi1_JetBin4_baseline_withoutDeltaPhi_pred->Write();
        DeltaPhi2_JetBin4_baseline_withoutDeltaPhi_pred->Write();
        DeltaPhi3_JetBin4_baseline_withoutDeltaPhi_pred->Write();

    }

    if (VBF) {

        VBF_dPhi_presel_pred->Write();
        VBF_dEta_presel_pred->Write();
        VBF_Mjj_presel_pred->Write();
        VBF_Jet1Pt_presel_pred->Write();
        VBF_Jet2Pt_presel_pred->Write();
        VBF_Jet3Pt_presel_pred->Write();
        VBF_Jet1Eta_presel_pred->Write();
        VBF_Jet2Eta_presel_pred->Write();
        VBF_Jet3Eta_presel_pred->Write();
        VBF_PTjj_presel_pred->Write();
        VBF_MET_presel_pred->Write();
        VBF_METsoft_presel_pred->Write();
        VBF_METsig_presel_pred->Write();
        VBF_MHTsig_presel_pred->Write();
        VBF_minDeltaPhiPTj12_presel_pred->Write();
        VBF_maxDeltaPhiPTj12_presel_pred->Write();
        VBF_DeltaPhiPTj3_presel_pred->Write();

        VBF_dPhi_presel_4JV_dPhiSide_pred->Write();
        VBF_dEta_presel_4JV_dPhiSide_pred->Write();
        VBF_Mjj_presel_4JV_dPhiSide_pred->Write();
        VBF_Jet1Pt_presel_4JV_dPhiSide_pred->Write();
        VBF_Jet2Pt_presel_4JV_dPhiSide_pred->Write();
        VBF_Jet3Pt_presel_4JV_dPhiSide_pred->Write();
        VBF_Jet1Eta_presel_4JV_dPhiSide_pred->Write();
        VBF_Jet2Eta_presel_4JV_dPhiSide_pred->Write();
        VBF_Jet3Eta_presel_4JV_dPhiSide_pred->Write();
        VBF_PTjj_presel_4JV_dPhiSide_pred->Write();
        VBF_MET_presel_4JV_dPhiSide_pred->Write();
        VBF_METsoft_presel_4JV_dPhiSide_pred->Write();
        VBF_METsig_presel_4JV_dPhiSide_pred->Write();
        VBF_MHTsig_presel_4JV_dPhiSide_pred->Write();
        VBF_minDeltaPhiPTj12_presel_4JV_dPhiSide_pred->Write();
        VBF_maxDeltaPhiPTj12_presel_4JV_dPhiSide_pred->Write();
        VBF_DeltaPhiPTj3_presel_4JV_dPhiSide_pred->Write();

        VBF_dPhi_dEta_pred->Write();
        VBF_dEta_dEta_pred->Write();
        VBF_Mjj_dEta_pred->Write();
        VBF_Jet1Pt_dEta_pred->Write();
        VBF_Jet2Pt_dEta_pred->Write();
        VBF_Jet3Pt_dEta_pred->Write();
        VBF_Jet1Eta_dEta_pred->Write();
        VBF_Jet2Eta_dEta_pred->Write();
        VBF_Jet3Eta_dEta_pred->Write();
        VBF_PTjj_dEta_pred->Write();
        VBF_MET_dEta_pred->Write();
        VBF_METsoft_dEta_pred->Write();
        VBF_METsig_dEta_pred->Write();
        VBF_MHTsig_dEta_pred->Write();
        VBF_minDeltaPhiPTj12_dEta_pred->Write();
        VBF_maxDeltaPhiPTj12_dEta_pred->Write();
        VBF_DeltaPhiPTj3_dEta_pred->Write();

        VBF_dPhi_jj_pred->Write();
        VBF_dEta_jj_pred->Write();
        VBF_Mjj_jj_pred->Write();
        VBF_Jet1Pt_jj_pred->Write();
        VBF_Jet2Pt_jj_pred->Write();
        VBF_Jet3Pt_jj_pred->Write();
        VBF_Jet1Eta_jj_pred->Write();
        VBF_Jet2Eta_jj_pred->Write();
        VBF_Jet3Eta_jj_pred->Write();
        VBF_PTjj_jj_pred->Write();
        VBF_MET_jj_pred->Write();
        VBF_METsoft_jj_pred->Write();
        VBF_METsig_jj_pred->Write();
        VBF_MHTsig_jj_pred->Write();
        VBF_minDeltaPhiPTj12_jj_pred->Write();
        VBF_maxDeltaPhiPTj12_jj_pred->Write();
        VBF_DeltaPhiPTj3_jj_pred->Write();

        VBF_dPhi_dEta_3JV_pred->Write();
        VBF_dEta_dEta_3JV_pred->Write();
        VBF_Mjj_dEta_3JV_pred->Write();
        VBF_Jet1Pt_dEta_3JV_pred->Write();
        VBF_Jet2Pt_dEta_3JV_pred->Write();
        VBF_Jet3Pt_dEta_3JV_pred->Write();
        VBF_Jet1Eta_dEta_3JV_pred->Write();
        VBF_Jet2Eta_dEta_3JV_pred->Write();
        VBF_Jet3Eta_dEta_3JV_pred->Write();
        VBF_PTjj_dEta_3JV_pred->Write();
        VBF_MET_dEta_3JV_pred->Write();
        VBF_METsoft_dEta_3JV_pred->Write();
        VBF_METsig_dEta_3JV_pred->Write();
        VBF_MHTsig_dEta_3JV_pred->Write();
        VBF_minDeltaPhiPTj12_dEta_3JV_pred->Write();
        VBF_maxDeltaPhiPTj12_dEta_3JV_pred->Write();
        VBF_DeltaPhiPTj3_dEta_3JV_pred->Write();

        VBF_dPhi_dEta_3JV_dPhiPTjj_pred->Write();
        VBF_dEta_dEta_3JV_dPhiPTjj_pred->Write();
        VBF_Mjj_dEta_3JV_dPhiPTjj_pred->Write();
        VBF_Jet1Pt_dEta_3JV_dPhiPTjj_pred->Write();
        VBF_Jet2Pt_dEta_3JV_dPhiPTjj_pred->Write();
        VBF_Jet3Pt_dEta_3JV_dPhiPTjj_pred->Write();
        VBF_Jet1Eta_dEta_3JV_dPhiPTjj_pred->Write();
        VBF_Jet2Eta_dEta_3JV_dPhiPTjj_pred->Write();
        VBF_Jet3Eta_dEta_3JV_dPhiPTjj_pred->Write();
        VBF_PTjj_dEta_3JV_dPhiPTjj_pred->Write();
        VBF_MET_dEta_3JV_dPhiPTjj_pred->Write();
        VBF_METsoft_dEta_3JV_dPhiPTjj_pred->Write();
        VBF_METsig_dEta_3JV_dPhiPTjj_pred->Write();
        VBF_MHTsig_dEta_3JV_dPhiPTjj_pred->Write();
        VBF_minDeltaPhiPTj12_dEta_3JV_dPhiPTjj_pred->Write();
        VBF_maxDeltaPhiPTj12_dEta_3JV_dPhiPTjj_pred->Write();
        VBF_DeltaPhiPTj3_dEta_3JV_dPhiPTjj_pred->Write();

    }

    cout << "after saving predictions" << endl;

    // Save histograms: expectation

    if (HTMHT) {

        HT_presel_sel->Write();
        MHT_presel_sel->Write();
        MET_presel_sel->Write();
        NJets_presel_sel->Write();
        NBJets_presel_sel->Write();
        Jet1Pt_presel_sel->Write();
        Jet2Pt_presel_sel->Write();
        Jet3Pt_presel_sel->Write();
        Jet1Eta_presel_sel->Write();
        Jet2Eta_presel_sel->Write();
        Jet3Eta_presel_sel->Write();
        HT_deltaPhi_sel->Write();
        MHT_deltaPhi_sel->Write();
        MET_deltaPhi_sel->Write();
        Jet1Pt_deltaPhi_sel->Write();
        Jet2Pt_deltaPhi_sel->Write();
        Jet3Pt_deltaPhi_sel->Write();
        Jet1Eta_deltaPhi_sel->Write();
        Jet2Eta_deltaPhi_sel->Write();
        Jet3Eta_deltaPhi_sel->Write();

        NJets_baseline_withoutMET_sel->Write();
        NJets_baseline_sel->Write();
        NJets_baseline_withoutDeltaPhi_withoutMET_sel->Write();
        NJets_baseline_withoutDeltaPhi_sel ->Write();

        NBJets_baseline_withoutMET_sel->Write();
        NBJets_baseline_sel->Write();
        NBJets_baseline_withoutDeltaPhi_withoutMET_sel->Write();
        NBJets_baseline_withoutDeltaPhi_sel ->Write();

        HT_baseline_sel->Write();
        MHT_baseline_sel->Write();
        MET_baseline_sel->Write();

        MHT_JetBin1_HTinclusive_sel->Write();
        MHT_JetBin2_HTinclusive_sel->Write();
        MHT_JetBin3_HTinclusive_sel->Write();
        MHT_JetBin4_HTinclusive_sel->Write();

        MET_JetBin1_HTinclusive_sel->Write();
        MET_JetBin2_HTinclusive_sel->Write();
        MET_JetBin3_HTinclusive_sel->Write();
        MET_JetBin4_HTinclusive_sel->Write();

        Jet1Pt_JetBin1_baseline_sel->Write();
        Jet2Pt_JetBin1_baseline_sel->Write();
        Jet1Eta_JetBin1_baseline_sel->Write();
        Jet2Eta_JetBin1_baseline_sel->Write();
        DeltaPhi1_JetBin1_baseline_sel->Write();
        DeltaPhi2_JetBin1_baseline_sel->Write();

        Jet1Pt_JetBin2_baseline_sel->Write();
        Jet2Pt_JetBin2_baseline_sel->Write();
        Jet3Pt_JetBin2_baseline_sel->Write();
        Jet1Eta_JetBin2_baseline_sel->Write();
        Jet2Eta_JetBin2_baseline_sel->Write();
        Jet3Eta_JetBin2_baseline_sel->Write();
        DeltaPhi1_JetBin2_baseline_sel->Write();
        DeltaPhi2_JetBin2_baseline_sel->Write();
        DeltaPhi3_JetBin2_baseline_sel->Write();

        Jet1Pt_JetBin3_baseline_sel->Write();
        Jet2Pt_JetBin3_baseline_sel->Write();
        Jet3Pt_JetBin3_baseline_sel->Write();
        Jet1Eta_JetBin3_baseline_sel->Write();
        Jet2Eta_JetBin3_baseline_sel->Write();
        Jet3Eta_JetBin3_baseline_sel->Write();
        DeltaPhi1_JetBin3_baseline_sel->Write();
        DeltaPhi2_JetBin3_baseline_sel->Write();
        DeltaPhi3_JetBin3_baseline_sel->Write();

        Jet1Pt_JetBin4_baseline_sel->Write();
        Jet2Pt_JetBin4_baseline_sel->Write();
        Jet3Pt_JetBin4_baseline_sel->Write();
        Jet1Eta_JetBin4_baseline_sel->Write();
        Jet2Eta_JetBin4_baseline_sel->Write();
        Jet3Eta_JetBin4_baseline_sel->Write();
        DeltaPhi1_JetBin4_baseline_sel->Write();
        DeltaPhi2_JetBin4_baseline_sel->Write();
        DeltaPhi3_JetBin4_baseline_sel->Write();

        HT_JetBin1_baseline_withoutDeltaPhi_sel->Write();
        MHT_JetBin1_baseline_withoutDeltaPhi_sel->Write();
        MET_JetBin1_baseline_withoutDeltaPhi_sel->Write();
        Jet1Pt_JetBin1_baseline_withoutDeltaPhi_sel->Write();
        Jet2Pt_JetBin1_baseline_withoutDeltaPhi_sel->Write();
        Jet1Eta_JetBin1_baseline_withoutDeltaPhi_sel->Write();
        Jet2Eta_JetBin1_baseline_withoutDeltaPhi_sel->Write();
        DeltaPhi1_JetBin1_baseline_withoutDeltaPhi_sel->Write();
        DeltaPhi2_JetBin1_baseline_withoutDeltaPhi_sel->Write();

        HT_JetBin2_baseline_withoutDeltaPhi_sel->Write();
        MHT_JetBin2_baseline_withoutDeltaPhi_sel->Write();
        MET_JetBin2_baseline_withoutDeltaPhi_sel->Write();
        Jet1Pt_JetBin2_baseline_withoutDeltaPhi_sel->Write();
        Jet2Pt_JetBin2_baseline_withoutDeltaPhi_sel->Write();
        Jet3Pt_JetBin2_baseline_withoutDeltaPhi_sel->Write();
        Jet1Eta_JetBin2_baseline_withoutDeltaPhi_sel->Write();
        Jet2Eta_JetBin2_baseline_withoutDeltaPhi_sel->Write();
        Jet3Eta_JetBin2_baseline_withoutDeltaPhi_sel->Write();
        DeltaPhi1_JetBin2_baseline_withoutDeltaPhi_sel->Write();
        DeltaPhi2_JetBin2_baseline_withoutDeltaPhi_sel->Write();
        DeltaPhi3_JetBin2_baseline_withoutDeltaPhi_sel->Write();

        HT_JetBin3_baseline_withoutDeltaPhi_sel->Write();
        MHT_JetBin3_baseline_withoutDeltaPhi_sel->Write();
        MET_JetBin3_baseline_withoutDeltaPhi_sel->Write();
        Jet1Pt_JetBin3_baseline_withoutDeltaPhi_sel->Write();
        Jet2Pt_JetBin3_baseline_withoutDeltaPhi_sel->Write();
        Jet3Pt_JetBin3_baseline_withoutDeltaPhi_sel->Write();
        Jet1Eta_JetBin3_baseline_withoutDeltaPhi_sel->Write();
        Jet2Eta_JetBin3_baseline_withoutDeltaPhi_sel->Write();
        Jet3Eta_JetBin3_baseline_withoutDeltaPhi_sel->Write();
        DeltaPhi1_JetBin3_baseline_withoutDeltaPhi_sel->Write();
        DeltaPhi2_JetBin3_baseline_withoutDeltaPhi_sel->Write();
        DeltaPhi3_JetBin3_baseline_withoutDeltaPhi_sel->Write();

        HT_JetBin4_baseline_withoutDeltaPhi_sel->Write();
        MHT_JetBin4_baseline_withoutDeltaPhi_sel->Write();
        MET_JetBin4_baseline_withoutDeltaPhi_sel->Write();
        Jet1Pt_JetBin4_baseline_withoutDeltaPhi_sel->Write();
        Jet2Pt_JetBin4_baseline_withoutDeltaPhi_sel->Write();
        Jet3Pt_JetBin4_baseline_withoutDeltaPhi_sel->Write();
        Jet1Eta_JetBin4_baseline_withoutDeltaPhi_sel->Write();
        Jet2Eta_JetBin4_baseline_withoutDeltaPhi_sel->Write();
        Jet3Eta_JetBin4_baseline_withoutDeltaPhi_sel->Write();
        DeltaPhi1_JetBin4_baseline_withoutDeltaPhi_sel->Write();
        DeltaPhi2_JetBin4_baseline_withoutDeltaPhi_sel->Write();
        DeltaPhi3_JetBin4_baseline_withoutDeltaPhi_sel->Write();

    }

    if (VBF) {

        VBF_dPhi_presel_sel->Write();
        VBF_dEta_presel_sel->Write();
        VBF_Mjj_presel_sel->Write();
        VBF_Jet1Pt_presel_sel->Write();
        VBF_Jet2Pt_presel_sel->Write();
        VBF_Jet3Pt_presel_sel->Write();
        VBF_Jet1Eta_presel_sel->Write();
        VBF_Jet2Eta_presel_sel->Write();
        VBF_Jet3Eta_presel_sel->Write();
        VBF_PTjj_presel_sel->Write();
        VBF_MET_presel_sel->Write();
        VBF_METsoft_presel_sel->Write();
        VBF_METsig_presel_sel->Write();
        VBF_MHTsig_presel_sel->Write();
        VBF_minDeltaPhiPTj12_presel_sel->Write();
        VBF_maxDeltaPhiPTj12_presel_sel->Write();
        VBF_DeltaPhiPTj3_presel_sel->Write();

        VBF_dPhi_presel_4JV_dPhiSide_sel->Write();
        VBF_dEta_presel_4JV_dPhiSide_sel->Write();
        VBF_Mjj_presel_4JV_dPhiSide_sel->Write();
        VBF_Jet1Pt_presel_4JV_dPhiSide_sel->Write();
        VBF_Jet2Pt_presel_4JV_dPhiSide_sel->Write();
        VBF_Jet3Pt_presel_4JV_dPhiSide_sel->Write();
        VBF_Jet1Eta_presel_4JV_dPhiSide_sel->Write();
        VBF_Jet2Eta_presel_4JV_dPhiSide_sel->Write();
        VBF_Jet3Eta_presel_4JV_dPhiSide_sel->Write();
        VBF_PTjj_presel_4JV_dPhiSide_sel->Write();
        VBF_MET_presel_4JV_dPhiSide_sel->Write();
        VBF_METsoft_presel_4JV_dPhiSide_sel->Write();
        VBF_METsig_presel_4JV_dPhiSide_sel->Write();
        VBF_MHTsig_presel_4JV_dPhiSide_sel->Write();
        VBF_minDeltaPhiPTj12_presel_4JV_dPhiSide_sel->Write();
        VBF_maxDeltaPhiPTj12_presel_4JV_dPhiSide_sel->Write();
        VBF_DeltaPhiPTj3_presel_4JV_dPhiSide_sel->Write();

        VBF_dPhi_dEta_sel->Write();
        VBF_dEta_dEta_sel->Write();
        VBF_Mjj_dEta_sel->Write();
        VBF_Jet1Pt_dEta_sel->Write();
        VBF_Jet2Pt_dEta_sel->Write();
        VBF_Jet3Pt_dEta_sel->Write();
        VBF_Jet1Eta_dEta_sel->Write();
        VBF_Jet2Eta_dEta_sel->Write();
        VBF_Jet3Eta_dEta_sel->Write();
        VBF_PTjj_dEta_sel->Write();
        VBF_MET_dEta_sel->Write();
        VBF_METsoft_dEta_sel->Write();
        VBF_METsig_dEta_sel->Write();
        VBF_MHTsig_dEta_sel->Write();
        VBF_minDeltaPhiPTj12_dEta_sel->Write();
        VBF_maxDeltaPhiPTj12_dEta_sel->Write();
        VBF_DeltaPhiPTj3_dEta_sel->Write();

        VBF_dPhi_jj_sel->Write();
        VBF_dEta_jj_sel->Write();
        VBF_Mjj_jj_sel->Write();
        VBF_Jet1Pt_jj_sel->Write();
        VBF_Jet2Pt_jj_sel->Write();
        VBF_Jet3Pt_jj_sel->Write();
        VBF_Jet1Eta_jj_sel->Write();
        VBF_Jet2Eta_jj_sel->Write();
        VBF_Jet3Eta_jj_sel->Write();
        VBF_PTjj_jj_sel->Write();
        VBF_MET_jj_sel->Write();
        VBF_METsoft_jj_sel->Write();
        VBF_METsig_jj_sel->Write();
        VBF_MHTsig_jj_sel->Write();
        VBF_minDeltaPhiPTj12_jj_sel->Write();
        VBF_maxDeltaPhiPTj12_jj_sel->Write();
        VBF_DeltaPhiPTj3_jj_sel->Write();

        VBF_dPhi_dEta_3JV_sel->Write();
        VBF_dEta_dEta_3JV_sel->Write();
        VBF_Mjj_dEta_3JV_sel->Write();
        VBF_Jet1Pt_dEta_3JV_sel->Write();
        VBF_Jet2Pt_dEta_3JV_sel->Write();
        VBF_Jet3Pt_dEta_3JV_sel->Write();
        VBF_Jet1Eta_dEta_3JV_sel->Write();
        VBF_Jet2Eta_dEta_3JV_sel->Write();
        VBF_Jet3Eta_dEta_3JV_sel->Write();
        VBF_PTjj_dEta_3JV_sel->Write();
        VBF_MET_dEta_3JV_sel->Write();
        VBF_METsoft_dEta_3JV_sel->Write();
        VBF_METsig_dEta_3JV_sel->Write();
        VBF_MHTsig_dEta_3JV_sel->Write();
        VBF_minDeltaPhiPTj12_dEta_3JV_sel->Write();
        VBF_maxDeltaPhiPTj12_dEta_3JV_sel->Write();
        VBF_DeltaPhiPTj3_dEta_3JV_sel->Write();

        VBF_dPhi_dEta_3JV_dPhiPTjj_sel->Write();
        VBF_dEta_dEta_3JV_dPhiPTjj_sel->Write();
        VBF_Mjj_dEta_3JV_dPhiPTjj_sel->Write();
        VBF_Jet1Pt_dEta_3JV_dPhiPTjj_sel->Write();
        VBF_Jet2Pt_dEta_3JV_dPhiPTjj_sel->Write();
        VBF_Jet3Pt_dEta_3JV_dPhiPTjj_sel->Write();
        VBF_Jet1Eta_dEta_3JV_dPhiPTjj_sel->Write();
        VBF_Jet2Eta_dEta_3JV_dPhiPTjj_sel->Write();
        VBF_Jet3Eta_dEta_3JV_dPhiPTjj_sel->Write();
        VBF_PTjj_dEta_3JV_dPhiPTjj_sel->Write();
        VBF_MET_dEta_3JV_dPhiPTjj_sel->Write();
        VBF_METsoft_dEta_3JV_dPhiPTjj_sel->Write();
        VBF_METsig_dEta_3JV_dPhiPTjj_sel->Write();
        VBF_MHTsig_dEta_3JV_dPhiPTjj_sel->Write();
        VBF_minDeltaPhiPTj12_dEta_3JV_dPhiPTjj_sel->Write();
        VBF_maxDeltaPhiPTj12_dEta_3JV_dPhiPTjj_sel->Write();
        VBF_DeltaPhiPTj3_dEta_3JV_dPhiPTjj_sel->Write();

    }

    cout << "after saving expectation" << endl;

    prediction_histos->Write();

    double value, error;
    cout << "Yields from R+S:" << endl << endl;

    value = VBF_Mjj_presel_4JV_dPhiSide_pred->IntegralAndError (2, 6, error);
    cout << "VR (dPhi side band): " << value << "+-" << error << endl;

    value = VBF_Mjj_jj_pred->IntegralAndError (2, 6, error);
    cout << "VR (jj, loose dEta): " << value << "+-" << error << endl;

    value = VBF_Mjj_dEta_3JV_pred->IntegralAndError (2, 6, error);
    cout << "SR: " << value << "+-" << error << endl;

    value = VBF_Mjj_dEta_3JV_pred->IntegralAndError (3, 3, error);
    cout << "SR1: " << value << "+-" << error << endl;

    value = VBF_Mjj_dEta_3JV_pred->IntegralAndError (4, 4, error);
    cout << "SR2: " << value << "+-" << error << endl;

    value = VBF_Mjj_dEta_3JV_pred->IntegralAndError (5, 6, error);
    cout << "SR3: " << value << "+-" << error << endl;

    value = VBF_Mjj_dEta_3JV_dPhiPTjj_pred->IntegralAndError (2, 6, error);
    cout << "SR (dPhi(MET,j1/2) > 1.0): " << value << "+-" << error << endl;

    value = VBF_Mjj_dEta_3JV_dPhiPTjj_pred->IntegralAndError (3, 3, error);
    cout << "SR1 (dPhi(MET,j1/2) > 1.0): " << value << "+-" << error << endl;

    value = VBF_Mjj_dEta_3JV_dPhiPTjj_pred->IntegralAndError (4, 4, error);
    cout << "SR2 (dPhi(MET,j1/2) > 1.0): " << value << "+-" << error << endl;

    value = VBF_Mjj_dEta_3JV_dPhiPTjj_pred->IntegralAndError (5, 6, error);
    cout << "SR3 (dPhi(MET,j1/2) > 1.0): " << value << "+-" << error << endl;

}


////////////////////////////////////////////////////////////////////////////////////////
bool Prediction::DeltaPhiCut()
{
    bool deltaPhiCut = true;
    if( NJets == 2 ) {
        if( DeltaPhi->at(0) < 0.5 ||
                DeltaPhi->at(1) < 0.5 ) deltaPhiCut = false;
    }
    if( NJets >= 3 ) {
        if( DeltaPhi->at(0) < 0.5 ||
                DeltaPhi->at(1) < 0.5 ||
                DeltaPhi->at(2) < 0.3 ) deltaPhiCut = false;
    }

    return deltaPhiCut;
}
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
double Prediction::CalcMHTjj() {
    TLorentzVector v1(0., 0., 0., 0.);
    v1.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetM->at(0));
    TLorentzVector v2(0., 0., 0., 0.);
    v2.SetPtEtaPhiM(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetM->at(1));
    double result = (v1+v2).Pt();
    //cout << "MHTjj: " << result << endl;

    return result;
}
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
double Prediction::CalcMjj() {
    TLorentzVector v1(0., 0., 0., 0.);
    v1.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetM->at(0));
    TLorentzVector v2(0., 0., 0., 0.);
    v2.SetPtEtaPhiM(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetM->at(1));
    double result = (v1+v2).M();
    return result;
}
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
bool Prediction::CalcDPhiMET(double& DPhiMET1, double& DPhiMET2, double& DPhiMET3, TLorentzVector& MET ) {
    TLorentzVector v1(0., 0., 0., 0.);
    v1.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetM->at(0));
    TLorentzVector v2(0., 0., 0., 0.);
    v2.SetPtEtaPhiM(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetM->at(1));
    TLorentzVector v3(0., 0., 0., 0.);
    v2.SetPtEtaPhiM(JetPt->at(2), JetEta->at(2), JetPhi->at(2), JetM->at(2));
    DPhiMET1 = fabs(MET.DeltaPhi(v1));
    DPhiMET2 = fabs(MET.DeltaPhi(v2));
    DPhiMET3 = 0.;
    if (JetPt->at(2) > 0.) fabs(MET.DeltaPhi(v3));
    return true;
}
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
double Prediction::CalcDeltaPhi() {
    TLorentzVector v1(0., 0., 0., 0.);
    v1.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetM->at(0));
    TLorentzVector v2(0., 0., 0., 0.);
    v2.SetPtEtaPhiM(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetM->at(1));
    double result = fabs(v1.DeltaPhi(v2));
    return result;
}
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
double Prediction::CalcDeltaEta() {
    TLorentzVector v1(0., 0., 0., 0.);
    v1.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetM->at(0));
    TLorentzVector v2(0., 0., 0., 0.);
    v2.SetPtEtaPhiM(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetM->at(1));
    double result = fabs(v1.Eta()-v2.Eta());
    return result;
}
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
bool Prediction::Soft3rd() {
    bool result = (JetPt->at(2) > 25. && JetPt->at(2) < 50. && fabs(JetEta->at(2)) < 4.5);
    return result;
}
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
bool Prediction::Veto3rd() {
    bool result = (JetPt->at(2) < 25. || fabs(JetEta->at(2)) > 4.5);
    return result;
}
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
bool Prediction::Veto4th() {
    bool result = (JetPt->at(3) < 25. || fabs(JetEta->at(3)) > 4.5);
    return result;
}
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
bool Prediction::MjjJetSel() {
    bool result = (JetPt->at(0) > 80. && JetPt->at(1) > 50. && fabs(JetEta->at(0)) < 4.5 && fabs(JetEta->at(1)) < 4.5 && (JetEta->at(0)*JetEta->at(1)) < 0.);
    return result;
}
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
void Prediction::DoRebinning(TH2F* prediction_raw, TH1F* selection_raw, int Nbins)
{

    TH2F* temp = (TH2F*) prediction_raw->Clone();

    //do some non-equidistant re-binning

    // MHT & MET binning
    if (Nbins == -2) {

        int nbins = 16;
        double vbins[17] = { 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 120., 150., 200., 250., 350., 500.};

        prediction_raw->GetXaxis()->Set(nbins, &vbins[0]);
        for (int j = 0; j <= prediction_raw->GetYaxis()->GetNbins() + 1; ++j) {
            for (int i = 0; i <= prediction_raw->GetXaxis()->GetNbins() + 1; ++i) {
                prediction_raw->SetBinContent(i, j, 0);
                prediction_raw->SetBinError(i, j, 0);
            }
        }

        //loop over y-axis
        for (int j = 1; j < temp->GetYaxis()->GetNbins() + 1; ++j) {
            int bin = 0;
            double sum2 = 0., content = 0.;
            for (int i = 0; i <= temp->GetXaxis()->GetNbins() + 1; ++i) {
                int this_bin = prediction_raw->GetXaxis()->FindBin(temp->GetXaxis()->GetBinCenter(i));
                if (this_bin > bin) {
                    double binWidth = prediction_raw->GetXaxis()->GetBinWidth(bin);
                    binWidth = 1.; //// for event counts
                    prediction_raw->SetBinContent(bin, j, content/binWidth);
                    prediction_raw->SetBinError(bin, j, sqrt(sum2)/binWidth);
                    bin = this_bin;
                    sum2 = content = 0.;
                }
                sum2 += pow(temp->GetBinError(i, j), 2);
                content += temp->GetBinContent(i, j);
            }

        }

        TH1F* temp2 = (TH1F*) selection_raw->Clone();
        selection_raw->GetXaxis()->Set(nbins, &vbins[0]);
        for (int i = 0; i <= selection_raw->GetXaxis()->GetNbins() + 1; ++i) {
            selection_raw->SetBinContent(i, 0);
            selection_raw->SetBinError(i, 0);
        }

        int bin = 0;
        double sum2 = 0., content = 0.;
        for (int i = 0; i <= temp2->GetXaxis()->GetNbins() + 1; ++i) {
            int this_bin = selection_raw->GetXaxis()->FindBin(temp->GetXaxis()->GetBinCenter(i));
            if (this_bin > bin) {
                double binWidth = selection_raw->GetXaxis()->GetBinWidth(bin);
                binWidth = 1.; //// for event counts
                selection_raw->SetBinContent(bin, content/binWidth);
                selection_raw->SetBinError(bin, sqrt(sum2)/binWidth);
                bin = this_bin;
                sum2 = content = 0.;
            }
            sum2 += pow(temp2->GetBinError(i), 2);
            content += temp2->GetBinContent(i);
        }

    }

    // Mjj binning
    if (Nbins == -3) {

        int nbins = 6;
        double vbins[7] = { 0., 600., 1000., 1500., 2000., 3000., 4000.};

        prediction_raw->GetXaxis()->Set(nbins, &vbins[0]);
        for (int j = 0; j <= prediction_raw->GetYaxis()->GetNbins() + 1; ++j) {
            for (int i = 0; i <= prediction_raw->GetXaxis()->GetNbins() + 1; ++i) {
                prediction_raw->SetBinContent(i, j, 0);
                prediction_raw->SetBinError(i, j, 0);
            }
        }

        //loop over y-axis
        for (int j = 1; j < temp->GetYaxis()->GetNbins() + 1; ++j) {
            int bin = 0;
            double sum2 = 0., content = 0.;
            for (int i = 0; i <= temp->GetXaxis()->GetNbins() + 1; ++i) {
                int this_bin = prediction_raw->GetXaxis()->FindBin(temp->GetXaxis()->GetBinCenter(i));
                if (this_bin > bin) {
                    double binWidth = prediction_raw->GetXaxis()->GetBinWidth(bin);
                    binWidth = 1.;
                    prediction_raw->SetBinContent(bin, j, content/binWidth);
                    prediction_raw->SetBinError(bin, j, sqrt(sum2)/binWidth);
                    bin = this_bin;
                    sum2 = content = 0.;
                }
                sum2 += pow(temp->GetBinError(i, j), 2);
                content += temp->GetBinContent(i, j);
            }

        }

        TH1F* temp2 = (TH1F*) selection_raw->Clone();
        selection_raw->GetXaxis()->Set(nbins, &vbins[0]);
        for (int i = 0; i <= selection_raw->GetXaxis()->GetNbins() + 1; ++i) {
            selection_raw->SetBinContent(i, 0);
            selection_raw->SetBinError(i, 0);
        }

        int bin = 0;
        double sum2 = 0., content = 0.;
        for (int i = 0; i <= temp2->GetXaxis()->GetNbins() + 1; ++i) {
            int this_bin = selection_raw->GetXaxis()->FindBin(temp->GetXaxis()->GetBinCenter(i));
            if (this_bin > bin) {
                double binWidth = selection_raw->GetXaxis()->GetBinWidth(bin);
                binWidth = 1.;
                selection_raw->SetBinContent(bin, content/binWidth);
                selection_raw->SetBinError(bin, sqrt(sum2)/binWidth);
                bin = this_bin;
                sum2 = content = 0.;
            }
            sum2 += pow(temp2->GetBinError(i), 2);
            content += temp2->GetBinContent(i);
        }

    }

    // standard equidistant re-binning
    if (Nbins > 0) {
        prediction_raw->Rebin2D(Nbins, 1);
        selection_raw->Rebin(Nbins);
    }

}
////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////
TH1F* Prediction::CalcPrediction(TH2F* prediction_raw) {

    TH1F* prediction = new TH1F();
    prediction = (TH1F*) prediction_raw->ProjectionX();
    prediction->Reset();
    for (int i = 0; i <= prediction_raw->GetXaxis()->GetNbins() + 1; ++i) {
        TH1F h = *((TH1F*) prediction_raw->ProjectionY("py", i, i));

        double summ = 0;
        double sumv = 0;
        int N = 0;

        //// Calculate mean
        for (int j = 1; j <= h.GetNbinsX(); ++j) {
            summ += h.GetBinContent(j);
            ++N;
        }
        double mean = summ / N;

        //// Calculated variance
        for (int j = 1; j <= h.GetNbinsX(); ++j) {
            sumv += pow(mean - h.GetBinContent(j), 2);
        }
        double variance = sqrt(sumv / N);

        prediction->SetBinContent(i, mean);
        prediction->SetBinError(i, variance);
    }

    return prediction;
}
////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////
TH1F* Prediction::GetSelectionHisto(TString type) {

    if ( type == "HT_presel") return HT_presel_sel;
    if ( type == "MHT_presel") return MHT_presel_sel;
    if ( type == "MET_presel") return MET_presel_sel;
    if ( type == "NJets_presel") return NJets_presel_sel;
    if ( type == "NBJets_presel") return NBJets_presel_sel;
    if ( type == "Jet1Pt_presel") return Jet1Pt_presel_sel;
    if ( type == "Jet2Pt_presel") return Jet2Pt_presel_sel;
    if ( type == "Jet3Pt_presel") return Jet3Pt_presel_sel;
    if ( type == "Jet1Eta_presel") return Jet1Eta_presel_sel;
    if ( type == "Jet2Eta_presel") return Jet2Eta_presel_sel;
    if ( type == "Jet3Eta_presel") return Jet3Eta_presel_sel;
    if ( type == "DeltaPhi1_presel") return DeltaPhi1_presel_sel;
    if ( type == "DeltaPhi2_presel") return DeltaPhi2_presel_sel;
    if ( type == "DeltaPhi3_presel") return DeltaPhi3_presel_sel;

    if ( type == "HT_deltaPhi") return HT_deltaPhi_sel;
    if ( type == "MHT_deltaPhi") return MHT_deltaPhi_sel;
    if ( type == "MET_deltaPhi") return MET_deltaPhi_sel;
    if ( type == "Jet1Pt_deltaPhi") return Jet1Pt_deltaPhi_sel;
    if ( type == "Jet2Pt_deltaPhi") return Jet2Pt_deltaPhi_sel;
    if ( type == "Jet3Pt_deltaPhi") return Jet3Pt_deltaPhi_sel;
    if ( type == "Jet1Eta_deltaPhi") return Jet1Eta_deltaPhi_sel;
    if ( type == "Jet2Eta_deltaPhi") return Jet2Eta_deltaPhi_sel;
    if ( type == "Jet3Eta_deltaPhi") return Jet3Eta_deltaPhi_sel;

    if ( type == "MHT_JetBin1_HTinclusive") return MHT_JetBin1_HTinclusive_sel;
    if ( type == "MHT_JetBin2_HTinclusive") return MHT_JetBin2_HTinclusive_sel;
    if ( type == "MHT_JetBin3_HTinclusive") return MHT_JetBin3_HTinclusive_sel;
    if ( type == "MHT_JetBin4_HTinclusive") return MHT_JetBin4_HTinclusive_sel;

    if ( type == "MET_JetBin1_HTinclusive") return MET_JetBin1_HTinclusive_sel;
    if ( type == "MET_JetBin2_HTinclusive") return MET_JetBin2_HTinclusive_sel;
    if ( type == "MET_JetBin3_HTinclusive") return MET_JetBin3_HTinclusive_sel;
    if ( type == "MET_JetBin4_HTinclusive") return MET_JetBin4_HTinclusive_sel;

    if ( type == "HT_baseline") return HT_baseline_sel;
    if ( type == "MHT_baseline") return MHT_baseline_sel;
    if ( type == "MET_baseline") return MET_baseline_sel;

    if( type == "Jet1Pt_JetBin1_baseline" ) return Jet1Pt_JetBin1_baseline_sel;
    if( type == "Jet2Pt_JetBin1_baseline" ) return Jet2Pt_JetBin1_baseline_sel;
    if( type == "Jet1Eta_JetBin1_baseline" ) return Jet1Eta_JetBin1_baseline_sel;
    if( type == "Jet2Eta_JetBin1_baseline" ) return Jet2Eta_JetBin1_baseline_sel;
    if( type == "DeltaPhi1_JetBin1_baseline") return DeltaPhi1_JetBin1_baseline_sel;
    if( type == "DeltaPhi2_JetBin1_baseline") return DeltaPhi2_JetBin1_baseline_sel;

    if( type == "Jet1Pt_JetBin2_baseline" ) return Jet1Pt_JetBin2_baseline_sel;
    if( type == "Jet2Pt_JetBin2_baseline" ) return Jet2Pt_JetBin2_baseline_sel;
    if( type == "Jet3Pt_JetBin2_baseline" ) return Jet3Pt_JetBin2_baseline_sel;
    if( type == "Jet1Eta_JetBin2_baseline" ) return Jet1Eta_JetBin2_baseline_sel;
    if( type == "Jet2Eta_JetBin2_baseline" ) return Jet2Eta_JetBin2_baseline_sel;
    if( type == "Jet3Eta_JetBin2_baseline" ) return Jet3Eta_JetBin2_baseline_sel;
    if( type == "DeltaPhi1_JetBin2_baseline") return DeltaPhi1_JetBin2_baseline_sel;
    if( type == "DeltaPhi2_JetBin2_baseline") return DeltaPhi2_JetBin2_baseline_sel;
    if( type == "DeltaPhi3_JetBin2_baseline") return DeltaPhi3_JetBin2_baseline_sel;

    if( type == "Jet1Pt_JetBin3_baseline" ) return Jet1Pt_JetBin3_baseline_sel;
    if( type == "Jet2Pt_JetBin3_baseline" ) return Jet2Pt_JetBin3_baseline_sel;
    if( type == "Jet3Pt_JetBin3_baseline" ) return Jet3Pt_JetBin3_baseline_sel;
    if( type == "Jet1Eta_JetBin3_baseline" ) return Jet1Eta_JetBin3_baseline_sel;
    if( type == "Jet2Eta_JetBin3_baseline" ) return Jet2Eta_JetBin3_baseline_sel;
    if( type == "Jet3Eta_JetBin3_baseline" ) return Jet3Eta_JetBin3_baseline_sel;
    if( type == "DeltaPhi1_JetBin3_baseline") return DeltaPhi1_JetBin3_baseline_sel;
    if( type == "DeltaPhi2_JetBin3_baseline") return DeltaPhi2_JetBin3_baseline_sel;
    if( type == "DeltaPhi3_JetBin3_baseline") return DeltaPhi3_JetBin3_baseline_sel;

    if( type == "Jet1Pt_JetBin4_baseline" ) return Jet1Pt_JetBin4_baseline_sel;
    if( type == "Jet2Pt_JetBin4_baseline" ) return Jet2Pt_JetBin4_baseline_sel;
    if( type == "Jet3Pt_JetBin4_baseline" ) return Jet3Pt_JetBin4_baseline_sel;
    if( type == "Jet1Eta_JetBin4_baseline" ) return Jet1Eta_JetBin4_baseline_sel;
    if( type == "Jet2Eta_JetBin4_baseline" ) return Jet2Eta_JetBin4_baseline_sel;
    if( type == "Jet3Eta_JetBin4_baseline" ) return Jet3Eta_JetBin4_baseline_sel;
    if( type == "DeltaPhi1_JetBin4_baseline") return DeltaPhi1_JetBin4_baseline_sel;
    if( type == "DeltaPhi2_JetBin4_baseline") return DeltaPhi2_JetBin4_baseline_sel;
    if( type == "DeltaPhi3_JetBin4_baseline") return DeltaPhi3_JetBin4_baseline_sel;

    if( type == "HT_JetBin1_baseline_withoutDeltaPhi") return HT_JetBin1_baseline_withoutDeltaPhi_sel;
    if( type == "MHT_JetBin1_baseline_withoutDeltaPhi") return MHT_JetBin1_baseline_withoutDeltaPhi_sel;
    if( type == "MET_JetBin1_baseline_withoutDeltaPhi") return MET_JetBin1_baseline_withoutDeltaPhi_sel;
    if( type == "Jet1Pt_JetBin1_baseline_withoutDeltaPhi" ) return Jet1Pt_JetBin1_baseline_withoutDeltaPhi_sel;
    if( type == "Jet2Pt_JetBin1_baseline_withoutDeltaPhi" ) return Jet2Pt_JetBin1_baseline_withoutDeltaPhi_sel;
    if( type == "Jet1Eta_JetBin1_baseline_withoutDeltaPhi" ) return Jet1Eta_JetBin1_baseline_withoutDeltaPhi_sel;
    if( type == "Jet2Eta_JetBin1_baseline_withoutDeltaPhi" ) return Jet2Eta_JetBin1_baseline_withoutDeltaPhi_sel;
    if( type == "DeltaPhi1_JetBin1_baseline_withoutDeltaPhi") return DeltaPhi1_JetBin1_baseline_withoutDeltaPhi_sel;
    if( type == "DeltaPhi2_JetBin1_baseline_withoutDeltaPhi") return DeltaPhi2_JetBin1_baseline_withoutDeltaPhi_sel;

    if( type == "HT_JetBin2_baseline_withoutDeltaPhi") return HT_JetBin2_baseline_withoutDeltaPhi_sel;
    if( type == "MHT_JetBin2_baseline_withoutDeltaPhi") return MHT_JetBin2_baseline_withoutDeltaPhi_sel;
    if( type == "MET_JetBin2_baseline_withoutDeltaPhi") return MET_JetBin2_baseline_withoutDeltaPhi_sel;
    if( type == "Jet1Pt_JetBin2_baseline_withoutDeltaPhi" ) return Jet1Pt_JetBin2_baseline_withoutDeltaPhi_sel;
    if( type == "Jet2Pt_JetBin2_baseline_withoutDeltaPhi" ) return Jet2Pt_JetBin2_baseline_withoutDeltaPhi_sel;
    if( type == "Jet3Pt_JetBin2_baseline_withoutDeltaPhi" ) return Jet3Pt_JetBin2_baseline_withoutDeltaPhi_sel;
    if( type == "Jet1Eta_JetBin2_baseline_withoutDeltaPhi" ) return Jet1Eta_JetBin2_baseline_withoutDeltaPhi_sel;
    if( type == "Jet2Eta_JetBin2_baseline_withoutDeltaPhi" ) return Jet2Eta_JetBin2_baseline_withoutDeltaPhi_sel;
    if( type == "Jet3Eta_JetBin2_baseline_withoutDeltaPhi" ) return Jet3Eta_JetBin2_baseline_withoutDeltaPhi_sel;
    if( type == "DeltaPhi1_JetBin2_baseline_withoutDeltaPhi") return DeltaPhi1_JetBin2_baseline_withoutDeltaPhi_sel;
    if( type == "DeltaPhi2_JetBin2_baseline_withoutDeltaPhi") return DeltaPhi2_JetBin2_baseline_withoutDeltaPhi_sel;
    if( type == "DeltaPhi3_JetBin2_baseline_withoutDeltaPhi") return DeltaPhi3_JetBin2_baseline_withoutDeltaPhi_sel;

    if( type == "HT_JetBin3_baseline_withoutDeltaPhi") return HT_JetBin3_baseline_withoutDeltaPhi_sel;
    if( type == "MHT_JetBin3_baseline_withoutDeltaPhi") return MHT_JetBin3_baseline_withoutDeltaPhi_sel;
    if( type == "MET_JetBin3_baseline_withoutDeltaPhi") return MET_JetBin3_baseline_withoutDeltaPhi_sel;
    if( type == "Jet1Pt_JetBin3_baseline_withoutDeltaPhi" ) return Jet1Pt_JetBin3_baseline_withoutDeltaPhi_sel;
    if( type == "Jet2Pt_JetBin3_baseline_withoutDeltaPhi" ) return Jet2Pt_JetBin3_baseline_withoutDeltaPhi_sel;
    if( type == "Jet3Pt_JetBin3_baseline_withoutDeltaPhi" ) return Jet3Pt_JetBin3_baseline_withoutDeltaPhi_sel;
    if( type == "Jet1Eta_JetBin3_baseline_withoutDeltaPhi" ) return Jet1Eta_JetBin3_baseline_withoutDeltaPhi_sel;
    if( type == "Jet2Eta_JetBin3_baseline_withoutDeltaPhi" ) return Jet2Eta_JetBin3_baseline_withoutDeltaPhi_sel;
    if( type == "Jet3Eta_JetBin3_baseline_withoutDeltaPhi" ) return Jet3Eta_JetBin3_baseline_withoutDeltaPhi_sel;
    if( type == "DeltaPhi1_JetBin3_baseline_withoutDeltaPhi") return DeltaPhi1_JetBin3_baseline_withoutDeltaPhi_sel;
    if( type == "DeltaPhi2_JetBin3_baseline_withoutDeltaPhi") return DeltaPhi2_JetBin3_baseline_withoutDeltaPhi_sel;
    if( type == "DeltaPhi3_JetBin3_baseline_withoutDeltaPhi") return DeltaPhi3_JetBin3_baseline_withoutDeltaPhi_sel;

    if( type == "HT_JetBin4_baseline_withoutDeltaPhi") return HT_JetBin4_baseline_withoutDeltaPhi_sel;
    if( type == "MHT_JetBin4_baseline_withoutDeltaPhi") return MHT_JetBin4_baseline_withoutDeltaPhi_sel;
    if( type == "MET_JetBin4_baseline_withoutDeltaPhi") return MET_JetBin4_baseline_withoutDeltaPhi_sel;
    if( type == "Jet1Pt_JetBin4_baseline_withoutDeltaPhi" ) return Jet1Pt_JetBin4_baseline_withoutDeltaPhi_sel;
    if( type == "Jet2Pt_JetBin4_baseline_withoutDeltaPhi" ) return Jet2Pt_JetBin4_baseline_withoutDeltaPhi_sel;
    if( type == "Jet3Pt_JetBin4_baseline_withoutDeltaPhi" ) return Jet3Pt_JetBin4_baseline_withoutDeltaPhi_sel;
    if( type == "Jet1Eta_JetBin4_baseline_withoutDeltaPhi" ) return Jet1Eta_JetBin4_baseline_withoutDeltaPhi_sel;
    if( type == "Jet2Eta_JetBin4_baseline_withoutDeltaPhi" ) return Jet2Eta_JetBin4_baseline_withoutDeltaPhi_sel;
    if( type == "Jet3Eta_JetBin4_baseline_withoutDeltaPhi" ) return Jet3Eta_JetBin4_baseline_withoutDeltaPhi_sel;
    if( type == "DeltaPhi1_JetBin4_baseline_withoutDeltaPhi") return DeltaPhi1_JetBin4_baseline_withoutDeltaPhi_sel;
    if( type == "DeltaPhi2_JetBin4_baseline_withoutDeltaPhi") return DeltaPhi2_JetBin4_baseline_withoutDeltaPhi_sel;
    if( type == "DeltaPhi3_JetBin4_baseline_withoutDeltaPhi") return DeltaPhi3_JetBin4_baseline_withoutDeltaPhi_sel;

    if ( type == "NJets_baseline_withoutMET") return NJets_baseline_withoutMET_sel;
    if ( type == "NJets_baseline") return NJets_baseline_sel;
    if ( type == "NJets_baseline_withoutDeltaPhi_withoutMET") return NJets_baseline_withoutDeltaPhi_withoutMET_sel;
    if ( type == "NJets_baseline_withoutDeltaPhi") return NJets_baseline_withoutDeltaPhi_sel;

    if ( type == "NBJets_baseline_withoutMET") return NBJets_baseline_withoutMET_sel;
    if ( type == "NBJets_baseline") return NBJets_baseline_sel;
    if ( type == "NBJets_baseline_withoutDeltaPhi_withoutMET") return NBJets_baseline_withoutDeltaPhi_withoutMET_sel;
    if ( type == "NBJets_baseline_withoutDeltaPhi") return NBJets_baseline_withoutDeltaPhi_sel;

    if ( type == "VBF_dPhi_presel") return VBF_dPhi_presel_sel;
    if ( type == "VBF_dEta_presel") return VBF_dEta_presel_sel;
    if ( type == "VBF_Mjj_presel") return VBF_Mjj_presel_sel;
    if ( type == "VBF_Jet1Pt_presel") return VBF_Jet1Pt_presel_sel;
    if ( type == "VBF_Jet2Pt_presel") return VBF_Jet2Pt_presel_sel;
    if ( type == "VBF_Jet3Pt_presel") return VBF_Jet3Pt_presel_sel;
    if ( type == "VBF_Jet1Eta_presel") return VBF_Jet1Eta_presel_sel;
    if ( type == "VBF_Jet2Eta_presel") return VBF_Jet2Eta_presel_sel;
    if ( type == "VBF_Jet3Eta_presel") return VBF_Jet3Eta_presel_sel;
    if ( type == "VBF_PTjj_presel") return VBF_PTjj_presel_sel;
    if ( type == "VBF_MET_presel") return VBF_MET_presel_sel;
    if ( type == "VBF_METsoft_presel") return VBF_METsoft_presel_sel;
    if ( type == "VBF_METsig_presel") return VBF_METsig_presel_sel;
    if ( type == "VBF_MHTsig_presel") return VBF_MHTsig_presel_sel;
    if ( type == "VBF_minDeltaPhiPTj12_presel") return VBF_minDeltaPhiPTj12_presel_sel;
    if ( type == "VBF_maxDeltaPhiPTj12_presel") return VBF_maxDeltaPhiPTj12_presel_sel;
    if ( type == "VBF_DeltaPhiPTj3_presel") return VBF_DeltaPhiPTj3_presel_sel;

    if ( type == "VBF_dPhi_presel_4JV_dPhiSide") return VBF_dPhi_presel_4JV_dPhiSide_sel;
    if ( type == "VBF_dEta_presel_4JV_dPhiSide") return VBF_dEta_presel_4JV_dPhiSide_sel;
    if ( type == "VBF_Mjj_presel_4JV_dPhiSide") return VBF_Mjj_presel_4JV_dPhiSide_sel;
    if ( type == "VBF_Jet1Pt_presel_4JV_dPhiSide") return VBF_Jet1Pt_presel_4JV_dPhiSide_sel;
    if ( type == "VBF_Jet2Pt_presel_4JV_dPhiSide") return VBF_Jet2Pt_presel_4JV_dPhiSide_sel;
    if ( type == "VBF_Jet3Pt_presel_4JV_dPhiSide") return VBF_Jet3Pt_presel_4JV_dPhiSide_sel;
    if ( type == "VBF_Jet1Eta_presel_4JV_dPhiSide") return VBF_Jet1Eta_presel_4JV_dPhiSide_sel;
    if ( type == "VBF_Jet2Eta_presel_4JV_dPhiSide") return VBF_Jet2Eta_presel_4JV_dPhiSide_sel;
    if ( type == "VBF_Jet3Eta_presel_4JV_dPhiSide") return VBF_Jet3Eta_presel_4JV_dPhiSide_sel;
    if ( type == "VBF_PTjj_presel_4JV_dPhiSide") return VBF_PTjj_presel_4JV_dPhiSide_sel;
    if ( type == "VBF_MET_presel_4JV_dPhiSide") return VBF_MET_presel_4JV_dPhiSide_sel;
    if ( type == "VBF_METsoft_presel_4JV_dPhiSide") return VBF_METsoft_presel_4JV_dPhiSide_sel;
    if ( type == "VBF_METsig_presel_4JV_dPhiSide") return VBF_METsig_presel_4JV_dPhiSide_sel;
    if ( type == "VBF_MHTsig_presel_4JV_dPhiSide") return VBF_MHTsig_presel_4JV_dPhiSide_sel;
    if ( type == "VBF_minDeltaPhiPTj12_presel_4JV_dPhiSide") return VBF_minDeltaPhiPTj12_presel_4JV_dPhiSide_sel;
    if ( type == "VBF_maxDeltaPhiPTj12_presel_4JV_dPhiSide") return VBF_maxDeltaPhiPTj12_presel_4JV_dPhiSide_sel;
    if ( type == "VBF_DeltaPhiPTj3_presel_4JV_dPhiSide") return VBF_DeltaPhiPTj3_presel_4JV_dPhiSide_sel;

    if ( type == "VBF_dPhi_dEta") return VBF_dPhi_dEta_sel;
    if ( type == "VBF_dEta_dEta") return VBF_dEta_dEta_sel;
    if ( type == "VBF_Mjj_dEta") return VBF_Mjj_dEta_sel;
    if ( type == "VBF_Jet1Pt_dEta") return VBF_Jet1Pt_dEta_sel;
    if ( type == "VBF_Jet2Pt_dEta") return VBF_Jet2Pt_dEta_sel;
    if ( type == "VBF_Jet3Pt_dEta") return VBF_Jet3Pt_dEta_sel;
    if ( type == "VBF_Jet1Eta_dEta") return VBF_Jet1Eta_dEta_sel;
    if ( type == "VBF_Jet2Eta_dEta") return VBF_Jet2Eta_dEta_sel;
    if ( type == "VBF_Jet3Eta_dEta") return VBF_Jet3Eta_dEta_sel;
    if ( type == "VBF_PTjj_dEta") return VBF_PTjj_dEta_sel;
    if ( type == "VBF_MET_dEta") return VBF_MET_dEta_sel;
    if ( type == "VBF_METsoft_dEta") return VBF_METsoft_dEta_sel;
    if ( type == "VBF_METsig_dEta") return VBF_METsig_dEta_sel;
    if ( type == "VBF_MHTsig_dEta") return VBF_MHTsig_dEta_sel;
    if ( type == "VBF_minDeltaPhiPTj12_dEta") return VBF_minDeltaPhiPTj12_dEta_sel;
    if ( type == "VBF_maxDeltaPhiPTj12_dEta") return VBF_maxDeltaPhiPTj12_dEta_sel;
    if ( type == "VBF_DeltaPhiPTj3_dEta") return VBF_DeltaPhiPTj3_dEta_sel;

    if ( type == "VBF_dPhi_jj") return VBF_dPhi_jj_sel;
    if ( type == "VBF_dEta_jj") return VBF_dEta_jj_sel;
    if ( type == "VBF_Mjj_jj") return VBF_Mjj_jj_sel;
    if ( type == "VBF_Jet1Pt_jj") return VBF_Jet1Pt_jj_sel;
    if ( type == "VBF_Jet2Pt_jj") return VBF_Jet2Pt_jj_sel;
    if ( type == "VBF_Jet3Pt_jj") return VBF_Jet3Pt_jj_sel;
    if ( type == "VBF_Jet1Eta_jj") return VBF_Jet1Eta_jj_sel;
    if ( type == "VBF_Jet2Eta_jj") return VBF_Jet2Eta_jj_sel;
    if ( type == "VBF_Jet3Eta_jj") return VBF_Jet3Eta_jj_sel;
    if ( type == "VBF_PTjj_jj") return VBF_PTjj_jj_sel;
    if ( type == "VBF_MET_jj") return VBF_MET_jj_sel;
    if ( type == "VBF_METsoft_jj") return VBF_METsoft_jj_sel;
    if ( type == "VBF_METsig_jj") return VBF_METsig_jj_sel;
    if ( type == "VBF_MHTsig_jj") return VBF_MHTsig_jj_sel;
    if ( type == "VBF_minDeltaPhiPTj12_jj") return VBF_minDeltaPhiPTj12_jj_sel;
    if ( type == "VBF_maxDeltaPhiPTj12_jj") return VBF_maxDeltaPhiPTj12_jj_sel;
    if ( type == "VBF_DeltaPhiPTj3_jj") return VBF_DeltaPhiPTj3_jj_sel;

    if ( type == "VBF_dPhi_dEta_3JV") return VBF_dPhi_dEta_3JV_sel;
    if ( type == "VBF_dEta_dEta_3JV") return VBF_dEta_dEta_3JV_sel;
    if ( type == "VBF_Mjj_dEta_3JV") return VBF_Mjj_dEta_3JV_sel;
    if ( type == "VBF_Jet1Pt_dEta_3JV") return VBF_Jet1Pt_dEta_3JV_sel;
    if ( type == "VBF_Jet2Pt_dEta_3JV") return VBF_Jet2Pt_dEta_3JV_sel;
    if ( type == "VBF_Jet3Pt_dEta_3JV") return VBF_Jet3Pt_dEta_3JV_sel;
    if ( type == "VBF_Jet1Eta_dEta_3JV") return VBF_Jet1Eta_dEta_3JV_sel;
    if ( type == "VBF_Jet2Eta_dEta_3JV") return VBF_Jet2Eta_dEta_3JV_sel;
    if ( type == "VBF_Jet3Eta_dEta_3JV") return VBF_Jet3Eta_dEta_3JV_sel;
    if ( type == "VBF_PTjj_dEta_3JV") return VBF_PTjj_dEta_3JV_sel;
    if ( type == "VBF_MET_dEta_3JV") return VBF_MET_dEta_3JV_sel;
    if ( type == "VBF_METsoft_dEta_3JV") return VBF_METsoft_dEta_3JV_sel;
    if ( type == "VBF_METsig_dEta_3JV") return VBF_METsig_dEta_3JV_sel;
    if ( type == "VBF_MHTsig_dEta_3JV") return VBF_MHTsig_dEta_3JV_sel;
    if ( type == "VBF_minDeltaPhiPTj12_dEta_3JV") return VBF_minDeltaPhiPTj12_dEta_3JV_sel;
    if ( type == "VBF_maxDeltaPhiPTj12_dEta_3JV") return VBF_maxDeltaPhiPTj12_dEta_3JV_sel;
    if ( type == "VBF_DeltaPhiPTj3_dEta_3JV") return VBF_DeltaPhiPTj3_dEta_3JV_sel;

    if ( type == "VBF_dPhi_dEta_3JV_dPhiPTjj") return VBF_dPhi_dEta_3JV_dPhiPTjj_sel;
    if ( type == "VBF_dEta_dEta_3JV_dPhiPTjj") return VBF_dEta_dEta_3JV_dPhiPTjj_sel;
    if ( type == "VBF_Mjj_dEta_3JV_dPhiPTjj") return VBF_Mjj_dEta_3JV_dPhiPTjj_sel;
    if ( type == "VBF_Jet1Pt_dEta_3JV_dPhiPTjj") return VBF_Jet1Pt_dEta_3JV_dPhiPTjj_sel;
    if ( type == "VBF_Jet2Pt_dEta_3JV_dPhiPTjj") return VBF_Jet2Pt_dEta_3JV_dPhiPTjj_sel;
    if ( type == "VBF_Jet3Pt_dEta_3JV_dPhiPTjj") return VBF_Jet3Pt_dEta_3JV_dPhiPTjj_sel;
    if ( type == "VBF_Jet1Eta_dEta_3JV_dPhiPTjj") return VBF_Jet1Eta_dEta_3JV_dPhiPTjj_sel;
    if ( type == "VBF_Jet2Eta_dEta_3JV_dPhiPTjj") return VBF_Jet2Eta_dEta_3JV_dPhiPTjj_sel;
    if ( type == "VBF_Jet3Eta_dEta_3JV_dPhiPTjj") return VBF_Jet3Eta_dEta_3JV_dPhiPTjj_sel;
    if ( type == "VBF_PTjj_dEta_3JV_dPhiPTjj") return VBF_PTjj_dEta_3JV_dPhiPTjj_sel;
    if ( type == "VBF_MET_dEta_3JV_dPhiPTjj") return VBF_MET_dEta_3JV_dPhiPTjj_sel;
    if ( type == "VBF_METsoft_dEta_3JV_dPhiPTjj") return VBF_METsoft_dEta_3JV_dPhiPTjj_sel;
    if ( type == "VBF_METsig_dEta_3JV_dPhiPTjj") return VBF_METsig_dEta_3JV_dPhiPTjj_sel;
    if ( type == "VBF_MHTsig_dEta_3JV_dPhiPTjj") return VBF_MHTsig_dEta_3JV_dPhiPTjj_sel;
    if ( type == "VBF_minDeltaPhiPTj12_dEta_3JV_dPhiPTjj") return VBF_minDeltaPhiPTj12_dEta_3JV_dPhiPTjj_sel;
    if ( type == "VBF_maxDeltaPhiPTj12_dEta_3JV_dPhiPTjj") return VBF_maxDeltaPhiPTj12_dEta_3JV_dPhiPTjj_sel;
    if ( type == "VBF_DeltaPhiPTj3_dEta_3JV_dPhiPTjj") return VBF_DeltaPhiPTj3_dEta_3JV_dPhiPTjj_sel;

    else {
        cout << "Error: No valid hist type" << endl;
        return dummy;
    }
}
////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////
TH1F* Prediction::GetPredictionHisto(TString type) {

    if ( type == "HT_presel") return HT_presel_pred;
    if ( type == "MHT_presel") return MHT_presel_pred;
    if ( type == "MET_presel") return MET_presel_pred;
    if ( type == "NJets_presel") return NJets_presel_pred;
    if ( type == "NBJets_presel") return NBJets_presel_pred;
    if ( type == "Jet1Pt_presel") return Jet1Pt_presel_pred;
    if ( type == "Jet2Pt_presel") return Jet2Pt_presel_pred;
    if ( type == "Jet3Pt_presel") return Jet3Pt_presel_pred;
    if ( type == "Jet1Eta_presel") return Jet1Eta_presel_pred;
    if ( type == "Jet2Eta_presel") return Jet2Eta_presel_pred;
    if ( type == "Jet3Eta_presel") return Jet3Eta_presel_pred;
    if ( type == "DeltaPhi1_presel") return DeltaPhi1_presel_pred;
    if ( type == "DeltaPhi2_presel") return DeltaPhi2_presel_pred;
    if ( type == "DeltaPhi3_presel") return DeltaPhi3_presel_pred;

    if ( type == "HT_deltaPhi") return HT_deltaPhi_pred;
    if ( type == "MHT_deltaPhi") return MHT_deltaPhi_pred;
    if ( type == "MET_deltaPhi") return MET_deltaPhi_pred;
    if ( type == "Jet1Pt_deltaPhi") return Jet1Pt_deltaPhi_pred;
    if ( type == "Jet2Pt_deltaPhi") return Jet2Pt_deltaPhi_pred;
    if ( type == "Jet3Pt_deltaPhi") return Jet3Pt_deltaPhi_pred;
    if ( type == "Jet1Eta_deltaPhi") return Jet1Eta_deltaPhi_pred;
    if ( type == "Jet2Eta_deltaPhi") return Jet2Eta_deltaPhi_pred;
    if ( type == "Jet3Eta_deltaPhi") return Jet3Eta_deltaPhi_pred;

    if ( type == "MHT_JetBin1_HTinclusive") return MHT_JetBin1_HTinclusive_pred;
    if ( type == "MHT_JetBin2_HTinclusive") return MHT_JetBin2_HTinclusive_pred;
    if ( type == "MHT_JetBin3_HTinclusive") return MHT_JetBin3_HTinclusive_pred;
    if ( type == "MHT_JetBin4_HTinclusive") return MHT_JetBin4_HTinclusive_pred;

    if ( type == "MET_JetBin1_HTinclusive") return MET_JetBin1_HTinclusive_pred;
    if ( type == "MET_JetBin2_HTinclusive") return MET_JetBin2_HTinclusive_pred;
    if ( type == "MET_JetBin3_HTinclusive") return MET_JetBin3_HTinclusive_pred;
    if ( type == "MET_JetBin4_HTinclusive") return MET_JetBin4_HTinclusive_pred;

    if ( type == "HT_baseline") return HT_baseline_pred;
    if ( type == "MHT_baseline") return MHT_baseline_pred;
    if ( type == "MET_baseline") return MET_baseline_pred;

    if( type == "Jet1Pt_JetBin1_baseline" ) return Jet1Pt_JetBin1_baseline_pred;
    if( type == "Jet2Pt_JetBin1_baseline" ) return Jet2Pt_JetBin1_baseline_pred;
    if( type == "Jet1Eta_JetBin1_baseline" ) return Jet1Eta_JetBin1_baseline_pred;
    if( type == "Jet2Eta_JetBin1_baseline" ) return Jet2Eta_JetBin1_baseline_pred;
    if( type == "DeltaPhi1_JetBin1_baseline") return DeltaPhi1_JetBin1_baseline_pred;
    if( type == "DeltaPhi2_JetBin1_baseline") return DeltaPhi2_JetBin1_baseline_pred;

    if( type == "Jet1Pt_JetBin2_baseline" ) return Jet1Pt_JetBin2_baseline_pred;
    if( type == "Jet2Pt_JetBin2_baseline" ) return Jet2Pt_JetBin2_baseline_pred;
    if( type == "Jet3Pt_JetBin2_baseline" ) return Jet3Pt_JetBin2_baseline_pred;
    if( type == "Jet1Eta_JetBin2_baseline" ) return Jet1Eta_JetBin2_baseline_pred;
    if( type == "Jet2Eta_JetBin2_baseline" ) return Jet2Eta_JetBin2_baseline_pred;
    if( type == "Jet3Eta_JetBin2_baseline" ) return Jet3Eta_JetBin2_baseline_pred;
    if( type == "DeltaPhi1_JetBin2_baseline") return DeltaPhi1_JetBin2_baseline_pred;
    if( type == "DeltaPhi2_JetBin2_baseline") return DeltaPhi2_JetBin2_baseline_pred;
    if( type == "DeltaPhi3_JetBin2_baseline") return DeltaPhi3_JetBin2_baseline_pred;

    if( type == "Jet1Pt_JetBin3_baseline" ) return Jet1Pt_JetBin3_baseline_pred;
    if( type == "Jet2Pt_JetBin3_baseline" ) return Jet2Pt_JetBin3_baseline_pred;
    if( type == "Jet3Pt_JetBin3_baseline" ) return Jet3Pt_JetBin3_baseline_pred;
    if( type == "Jet1Eta_JetBin3_baseline" ) return Jet1Eta_JetBin3_baseline_pred;
    if( type == "Jet2Eta_JetBin3_baseline" ) return Jet2Eta_JetBin3_baseline_pred;
    if( type == "Jet3Eta_JetBin3_baseline" ) return Jet3Eta_JetBin3_baseline_pred;
    if( type == "DeltaPhi1_JetBin3_baseline") return DeltaPhi1_JetBin3_baseline_pred;
    if( type == "DeltaPhi2_JetBin3_baseline") return DeltaPhi2_JetBin3_baseline_pred;
    if( type == "DeltaPhi3_JetBin3_baseline") return DeltaPhi3_JetBin3_baseline_pred;

    if( type == "Jet1Pt_JetBin4_baseline" ) return Jet1Pt_JetBin4_baseline_pred;
    if( type == "Jet2Pt_JetBin4_baseline" ) return Jet2Pt_JetBin4_baseline_pred;
    if( type == "Jet3Pt_JetBin4_baseline" ) return Jet3Pt_JetBin4_baseline_pred;
    if( type == "Jet1Eta_JetBin4_baseline" ) return Jet1Eta_JetBin4_baseline_pred;
    if( type == "Jet2Eta_JetBin4_baseline" ) return Jet2Eta_JetBin4_baseline_pred;
    if( type == "Jet3Eta_JetBin4_baseline" ) return Jet3Eta_JetBin4_baseline_pred;
    if( type == "DeltaPhi1_JetBin4_baseline") return DeltaPhi1_JetBin4_baseline_pred;
    if( type == "DeltaPhi2_JetBin4_baseline") return DeltaPhi2_JetBin4_baseline_pred;
    if( type == "DeltaPhi3_JetBin4_baseline") return DeltaPhi3_JetBin4_baseline_pred;

    if( type == "HT_JetBin1_baseline_withoutDeltaPhi") return HT_JetBin1_baseline_withoutDeltaPhi_pred;
    if( type == "MHT_JetBin1_baseline_withoutDeltaPhi") return MHT_JetBin1_baseline_withoutDeltaPhi_pred;
    if( type == "MET_JetBin1_baseline_withoutDeltaPhi") return MET_JetBin1_baseline_withoutDeltaPhi_pred;
    if( type == "Jet1Pt_JetBin1_baseline_withoutDeltaPhi" ) return Jet1Pt_JetBin1_baseline_withoutDeltaPhi_pred;
    if( type == "Jet2Pt_JetBin1_baseline_withoutDeltaPhi" ) return Jet2Pt_JetBin1_baseline_withoutDeltaPhi_pred;
    if( type == "Jet1Eta_JetBin1_baseline_withoutDeltaPhi" ) return Jet1Eta_JetBin1_baseline_withoutDeltaPhi_pred;
    if( type == "Jet2Eta_JetBin1_baseline_withoutDeltaPhi" ) return Jet2Eta_JetBin1_baseline_withoutDeltaPhi_pred;
    if( type == "DeltaPhi1_JetBin1_baseline_withoutDeltaPhi") return DeltaPhi1_JetBin1_baseline_withoutDeltaPhi_pred;
    if( type == "DeltaPhi2_JetBin1_baseline_withoutDeltaPhi") return DeltaPhi2_JetBin1_baseline_withoutDeltaPhi_pred;

    if( type == "HT_JetBin2_baseline_withoutDeltaPhi") return HT_JetBin2_baseline_withoutDeltaPhi_pred;
    if( type == "MHT_JetBin2_baseline_withoutDeltaPhi") return MHT_JetBin2_baseline_withoutDeltaPhi_pred;
    if( type == "MET_JetBin2_baseline_withoutDeltaPhi") return MET_JetBin2_baseline_withoutDeltaPhi_pred;
    if( type == "Jet1Pt_JetBin2_baseline_withoutDeltaPhi" ) return Jet1Pt_JetBin2_baseline_withoutDeltaPhi_pred;
    if( type == "Jet2Pt_JetBin2_baseline_withoutDeltaPhi" ) return Jet2Pt_JetBin2_baseline_withoutDeltaPhi_pred;
    if( type == "Jet3Pt_JetBin2_baseline_withoutDeltaPhi" ) return Jet3Pt_JetBin2_baseline_withoutDeltaPhi_pred;
    if( type == "Jet1Eta_JetBin2_baseline_withoutDeltaPhi" ) return Jet1Eta_JetBin2_baseline_withoutDeltaPhi_pred;
    if( type == "Jet2Eta_JetBin2_baseline_withoutDeltaPhi" ) return Jet2Eta_JetBin2_baseline_withoutDeltaPhi_pred;
    if( type == "Jet3Eta_JetBin2_baseline_withoutDeltaPhi" ) return Jet3Eta_JetBin2_baseline_withoutDeltaPhi_pred;
    if( type == "DeltaPhi1_JetBin2_baseline_withoutDeltaPhi") return DeltaPhi1_JetBin2_baseline_withoutDeltaPhi_pred;
    if( type == "DeltaPhi2_JetBin2_baseline_withoutDeltaPhi") return DeltaPhi2_JetBin2_baseline_withoutDeltaPhi_pred;
    if( type == "DeltaPhi3_JetBin2_baseline_withoutDeltaPhi") return DeltaPhi3_JetBin2_baseline_withoutDeltaPhi_pred;

    if( type == "HT_JetBin3_baseline_withoutDeltaPhi") return HT_JetBin3_baseline_withoutDeltaPhi_pred;
    if( type == "MHT_JetBin3_baseline_withoutDeltaPhi") return MHT_JetBin3_baseline_withoutDeltaPhi_pred;
    if( type == "MET_JetBin3_baseline_withoutDeltaPhi") return MET_JetBin3_baseline_withoutDeltaPhi_pred;
    if( type == "Jet1Pt_JetBin3_baseline_withoutDeltaPhi" ) return Jet1Pt_JetBin3_baseline_withoutDeltaPhi_pred;
    if( type == "Jet2Pt_JetBin3_baseline_withoutDeltaPhi" ) return Jet2Pt_JetBin3_baseline_withoutDeltaPhi_pred;
    if( type == "Jet3Pt_JetBin3_baseline_withoutDeltaPhi" ) return Jet3Pt_JetBin3_baseline_withoutDeltaPhi_pred;
    if( type == "Jet1Eta_JetBin3_baseline_withoutDeltaPhi" ) return Jet1Eta_JetBin3_baseline_withoutDeltaPhi_pred;
    if( type == "Jet2Eta_JetBin3_baseline_withoutDeltaPhi" ) return Jet2Eta_JetBin3_baseline_withoutDeltaPhi_pred;
    if( type == "Jet3Eta_JetBin3_baseline_withoutDeltaPhi" ) return Jet3Eta_JetBin3_baseline_withoutDeltaPhi_pred;
    if( type == "DeltaPhi1_JetBin3_baseline_withoutDeltaPhi") return DeltaPhi1_JetBin3_baseline_withoutDeltaPhi_pred;
    if( type == "DeltaPhi2_JetBin3_baseline_withoutDeltaPhi") return DeltaPhi2_JetBin3_baseline_withoutDeltaPhi_pred;
    if( type == "DeltaPhi3_JetBin3_baseline_withoutDeltaPhi") return DeltaPhi3_JetBin3_baseline_withoutDeltaPhi_pred;

    if( type == "HT_JetBin4_baseline_withoutDeltaPhi") return HT_JetBin4_baseline_withoutDeltaPhi_pred;
    if( type == "MHT_JetBin4_baseline_withoutDeltaPhi") return MHT_JetBin4_baseline_withoutDeltaPhi_pred;
    if( type == "MET_JetBin4_baseline_withoutDeltaPhi") return MET_JetBin4_baseline_withoutDeltaPhi_pred;
    if( type == "Jet1Pt_JetBin4_baseline_withoutDeltaPhi" ) return Jet1Pt_JetBin4_baseline_withoutDeltaPhi_pred;
    if( type == "Jet2Pt_JetBin4_baseline_withoutDeltaPhi" ) return Jet2Pt_JetBin4_baseline_withoutDeltaPhi_pred;
    if( type == "Jet3Pt_JetBin4_baseline_withoutDeltaPhi" ) return Jet3Pt_JetBin4_baseline_withoutDeltaPhi_pred;
    if( type == "Jet1Eta_JetBin4_baseline_withoutDeltaPhi" ) return Jet1Eta_JetBin4_baseline_withoutDeltaPhi_pred;
    if( type == "Jet2Eta_JetBin4_baseline_withoutDeltaPhi" ) return Jet2Eta_JetBin4_baseline_withoutDeltaPhi_pred;
    if( type == "Jet3Eta_JetBin4_baseline_withoutDeltaPhi" ) return Jet3Eta_JetBin4_baseline_withoutDeltaPhi_pred;
    if( type == "DeltaPhi1_JetBin4_baseline_withoutDeltaPhi") return DeltaPhi1_JetBin4_baseline_withoutDeltaPhi_pred;
    if( type == "DeltaPhi2_JetBin4_baseline_withoutDeltaPhi") return DeltaPhi2_JetBin4_baseline_withoutDeltaPhi_pred;
    if( type == "DeltaPhi3_JetBin4_baseline_withoutDeltaPhi") return DeltaPhi3_JetBin4_baseline_withoutDeltaPhi_pred;

    if ( type == "NJets_baseline_withoutMET") return NJets_baseline_withoutMET_pred;
    if ( type == "NJets_baseline") return NJets_baseline_pred;
    if ( type == "NJets_baseline_withoutDeltaPhi_withoutMET") return NJets_baseline_withoutDeltaPhi_withoutMET_pred;
    if ( type == "NJets_baseline_withoutDeltaPhi") return NJets_baseline_withoutDeltaPhi_pred;

    if ( type == "NBJets_baseline_withoutMET") return NBJets_baseline_withoutMET_pred;
    if ( type == "NBJets_baseline") return NBJets_baseline_pred;
    if ( type == "NBJets_baseline_withoutDeltaPhi_withoutMET") return NBJets_baseline_withoutDeltaPhi_withoutMET_pred;
    if ( type == "NBJets_baseline_withoutDeltaPhi") return NBJets_baseline_withoutDeltaPhi_pred;

    if ( type == "VBF_dPhi_presel") return VBF_dPhi_presel_pred;
    if ( type == "VBF_dEta_presel") return VBF_dEta_presel_pred;
    if ( type == "VBF_Mjj_presel") return VBF_Mjj_presel_pred;
    if ( type == "VBF_Jet1Pt_presel") return VBF_Jet1Pt_presel_pred;
    if ( type == "VBF_Jet2Pt_presel") return VBF_Jet2Pt_presel_pred;
    if ( type == "VBF_Jet3Pt_presel") return VBF_Jet3Pt_presel_pred;
    if ( type == "VBF_Jet1Eta_presel") return VBF_Jet1Eta_presel_pred;
    if ( type == "VBF_Jet2Eta_presel") return VBF_Jet2Eta_presel_pred;
    if ( type == "VBF_Jet3Eta_presel") return VBF_Jet3Eta_presel_pred;
    if ( type == "VBF_PTjj_presel") return VBF_PTjj_presel_pred;
    if ( type == "VBF_MET_presel") return VBF_MET_presel_pred;
    if ( type == "VBF_METsoft_presel") return VBF_METsoft_presel_pred;
    if ( type == "VBF_METsig_presel") return VBF_METsig_presel_pred;
    if ( type == "VBF_MHTsig_presel") return VBF_MHTsig_presel_pred;
    if ( type == "VBF_minDeltaPhiPTj12_presel") return VBF_minDeltaPhiPTj12_presel_pred;
    if ( type == "VBF_maxDeltaPhiPTj12_presel") return VBF_maxDeltaPhiPTj12_presel_pred;
    if ( type == "VBF_DeltaPhiPTj3_presel") return VBF_DeltaPhiPTj3_presel_pred;

    if ( type == "VBF_dPhi_presel_4JV_dPhiSide") return VBF_dPhi_presel_4JV_dPhiSide_pred;
    if ( type == "VBF_dEta_presel_4JV_dPhiSide") return VBF_dEta_presel_4JV_dPhiSide_pred;
    if ( type == "VBF_Mjj_presel_4JV_dPhiSide") return VBF_Mjj_presel_4JV_dPhiSide_pred;
    if ( type == "VBF_Jet1Pt_presel_4JV_dPhiSide") return VBF_Jet1Pt_presel_4JV_dPhiSide_pred;
    if ( type == "VBF_Jet2Pt_presel_4JV_dPhiSide") return VBF_Jet2Pt_presel_4JV_dPhiSide_pred;
    if ( type == "VBF_Jet3Pt_presel_4JV_dPhiSide") return VBF_Jet3Pt_presel_4JV_dPhiSide_pred;
    if ( type == "VBF_Jet1Eta_presel_4JV_dPhiSide") return VBF_Jet1Eta_presel_4JV_dPhiSide_pred;
    if ( type == "VBF_Jet2Eta_presel_4JV_dPhiSide") return VBF_Jet2Eta_presel_4JV_dPhiSide_pred;
    if ( type == "VBF_Jet3Eta_presel_4JV_dPhiSide") return VBF_Jet3Eta_presel_4JV_dPhiSide_pred;
    if ( type == "VBF_PTjj_presel_4JV_dPhiSide") return VBF_PTjj_presel_4JV_dPhiSide_pred;
    if ( type == "VBF_MET_presel_4JV_dPhiSide") return VBF_MET_presel_4JV_dPhiSide_pred;
    if ( type == "VBF_METsoft_presel_4JV_dPhiSide") return VBF_METsoft_presel_4JV_dPhiSide_pred;
    if ( type == "VBF_METsig_presel_4JV_dPhiSide") return VBF_METsig_presel_4JV_dPhiSide_pred;
    if ( type == "VBF_MHTsig_presel_4JV_dPhiSide") return VBF_MHTsig_presel_4JV_dPhiSide_pred;
    if ( type == "VBF_minDeltaPhiPTj12_presel_4JV_dPhiSide") return VBF_minDeltaPhiPTj12_presel_4JV_dPhiSide_pred;
    if ( type == "VBF_maxDeltaPhiPTj12_presel_4JV_dPhiSide") return VBF_maxDeltaPhiPTj12_presel_4JV_dPhiSide_pred;
    if ( type == "VBF_DeltaPhiPTj3_presel_4JV_dPhiSide") return VBF_DeltaPhiPTj3_presel_4JV_dPhiSide_pred;

    if ( type == "VBF_dPhi_dEta") return VBF_dPhi_dEta_pred;
    if ( type == "VBF_dEta_dEta") return VBF_dEta_dEta_pred;
    if ( type == "VBF_Mjj_dEta") return VBF_Mjj_dEta_pred;
    if ( type == "VBF_Jet1Pt_dEta") return VBF_Jet1Pt_dEta_pred;
    if ( type == "VBF_Jet2Pt_dEta") return VBF_Jet2Pt_dEta_pred;
    if ( type == "VBF_Jet3Pt_dEta") return VBF_Jet3Pt_dEta_pred;
    if ( type == "VBF_Jet1Eta_dEta") return VBF_Jet1Eta_dEta_pred;
    if ( type == "VBF_Jet2Eta_dEta") return VBF_Jet2Eta_dEta_pred;
    if ( type == "VBF_Jet3Eta_dEta") return VBF_Jet3Eta_dEta_pred;
    if ( type == "VBF_PTjj_dEta") return VBF_PTjj_dEta_pred;
    if ( type == "VBF_MET_dEta") return VBF_MET_dEta_pred;
    if ( type == "VBF_METsoft_dEta") return VBF_METsoft_dEta_pred;
    if ( type == "VBF_METsig_dEta") return VBF_METsig_dEta_pred;
    if ( type == "VBF_MHTsig_dEta") return VBF_MHTsig_dEta_pred;
    if ( type == "VBF_minDeltaPhiPTj12_dEta") return VBF_minDeltaPhiPTj12_dEta_pred;
    if ( type == "VBF_maxDeltaPhiPTj12_dEta") return VBF_maxDeltaPhiPTj12_dEta_pred;
    if ( type == "VBF_DeltaPhiPTj3_dEta") return VBF_DeltaPhiPTj3_dEta_pred;

    if ( type == "VBF_dPhi_jj") return VBF_dPhi_jj_pred;
    if ( type == "VBF_dEta_jj") return VBF_dEta_jj_pred;
    if ( type == "VBF_Mjj_jj") return VBF_Mjj_jj_pred;
    if ( type == "VBF_Jet1Pt_jj") return VBF_Jet1Pt_jj_pred;
    if ( type == "VBF_Jet2Pt_jj") return VBF_Jet2Pt_jj_pred;
    if ( type == "VBF_Jet3Pt_jj") return VBF_Jet3Pt_jj_pred;
    if ( type == "VBF_Jet1Eta_jj") return VBF_Jet1Eta_jj_pred;
    if ( type == "VBF_Jet2Eta_jj") return VBF_Jet2Eta_jj_pred;
    if ( type == "VBF_Jet3Eta_jj") return VBF_Jet3Eta_jj_pred;
    if ( type == "VBF_PTjj_jj") return VBF_PTjj_jj_pred;
    if ( type == "VBF_MET_jj") return VBF_MET_jj_pred;
    if ( type == "VBF_METsoft_jj") return VBF_METsoft_jj_pred;
    if ( type == "VBF_METsig_jj") return VBF_METsig_jj_pred;
    if ( type == "VBF_MHTsig_jj") return VBF_MHTsig_jj_pred;
    if ( type == "VBF_minDeltaPhiPTj12_jj") return VBF_minDeltaPhiPTj12_jj_pred;
    if ( type == "VBF_maxDeltaPhiPTj12_jj") return VBF_maxDeltaPhiPTj12_jj_pred;
    if ( type == "VBF_DeltaPhiPTj3_jj") return VBF_DeltaPhiPTj3_jj_pred;

    if ( type == "VBF_dPhi_dEta_3JV") return VBF_dPhi_dEta_3JV_pred;
    if ( type == "VBF_dEta_dEta_3JV") return VBF_dEta_dEta_3JV_pred;
    if ( type == "VBF_Mjj_dEta_3JV") return VBF_Mjj_dEta_3JV_pred;
    if ( type == "VBF_Jet1Pt_dEta_3JV") return VBF_Jet1Pt_dEta_3JV_pred;
    if ( type == "VBF_Jet2Pt_dEta_3JV") return VBF_Jet2Pt_dEta_3JV_pred;
    if ( type == "VBF_Jet3Pt_dEta_3JV") return VBF_Jet3Pt_dEta_3JV_pred;
    if ( type == "VBF_Jet1Eta_dEta_3JV") return VBF_Jet1Eta_dEta_3JV_pred;
    if ( type == "VBF_Jet2Eta_dEta_3JV") return VBF_Jet2Eta_dEta_3JV_pred;
    if ( type == "VBF_Jet3Eta_dEta_3JV") return VBF_Jet3Eta_dEta_3JV_pred;
    if ( type == "VBF_PTjj_dEta_3JV") return VBF_PTjj_dEta_3JV_pred;
    if ( type == "VBF_MET_dEta_3JV") return VBF_MET_dEta_3JV_pred;
    if ( type == "VBF_METsoft_dEta_3JV") return VBF_METsoft_dEta_3JV_pred;
    if ( type == "VBF_METsig_dEta_3JV") return VBF_METsig_dEta_3JV_pred;
    if ( type == "VBF_MHTsig_dEta_3JV") return VBF_MHTsig_dEta_3JV_pred;
    if ( type == "VBF_minDeltaPhiPTj12_dEta_3JV") return VBF_minDeltaPhiPTj12_dEta_3JV_pred;
    if ( type == "VBF_maxDeltaPhiPTj12_dEta_3JV") return VBF_maxDeltaPhiPTj12_dEta_3JV_pred;
    if ( type == "VBF_DeltaPhiPTj3_dEta_3JV") return VBF_DeltaPhiPTj3_dEta_3JV_pred;

    if ( type == "VBF_dPhi_dEta_3JV_dPhiPTjj") return VBF_dPhi_dEta_3JV_dPhiPTjj_pred;
    if ( type == "VBF_dEta_dEta_3JV_dPhiPTjj") return VBF_dEta_dEta_3JV_dPhiPTjj_pred;
    if ( type == "VBF_Mjj_dEta_3JV_dPhiPTjj") return VBF_Mjj_dEta_3JV_dPhiPTjj_pred;
    if ( type == "VBF_Jet1Pt_dEta_3JV_dPhiPTjj") return VBF_Jet1Pt_dEta_3JV_dPhiPTjj_pred;
    if ( type == "VBF_Jet2Pt_dEta_3JV_dPhiPTjj") return VBF_Jet2Pt_dEta_3JV_dPhiPTjj_pred;
    if ( type == "VBF_Jet3Pt_dEta_3JV_dPhiPTjj") return VBF_Jet3Pt_dEta_3JV_dPhiPTjj_pred;
    if ( type == "VBF_Jet1Eta_dEta_3JV_dPhiPTjj") return VBF_Jet1Eta_dEta_3JV_dPhiPTjj_pred;
    if ( type == "VBF_Jet2Eta_dEta_3JV_dPhiPTjj") return VBF_Jet2Eta_dEta_3JV_dPhiPTjj_pred;
    if ( type == "VBF_Jet3Eta_dEta_3JV_dPhiPTjj") return VBF_Jet3Eta_dEta_3JV_dPhiPTjj_pred;
    if ( type == "VBF_PTjj_dEta_3JV_dPhiPTjj") return VBF_PTjj_dEta_3JV_dPhiPTjj_pred;
    if ( type == "VBF_MET_dEta_3JV_dPhiPTjj") return VBF_MET_dEta_3JV_dPhiPTjj_pred;
    if ( type == "VBF_METsoft_dEta_3JV_dPhiPTjj") return VBF_METsoft_dEta_3JV_dPhiPTjj_pred;
    if ( type == "VBF_METsig_dEta_3JV_dPhiPTjj") return VBF_METsig_dEta_3JV_dPhiPTjj_pred;
    if ( type == "VBF_MHTsig_dEta_3JV_dPhiPTjj") return VBF_MHTsig_dEta_3JV_dPhiPTjj_pred;
    if ( type == "VBF_minDeltaPhiPTj12_dEta_3JV_dPhiPTjj") return VBF_minDeltaPhiPTj12_dEta_3JV_dPhiPTjj_pred;
    if ( type == "VBF_maxDeltaPhiPTj12_dEta_3JV_dPhiPTjj") return VBF_maxDeltaPhiPTj12_dEta_3JV_dPhiPTjj_pred;
    if ( type == "VBF_DeltaPhiPTj3_dEta_3JV_dPhiPTjj") return VBF_DeltaPhiPTj3_dEta_3JV_dPhiPTjj_pred;

    else {
        cout << "Error: No valid hist type" << endl;
        return dummy;
    }
}
////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////
double Prediction::GetResultValue(TH1F* histo, double MHTlow, double MHTup)
{
    double result_value;
    if( MHTlow == MHTup ) {
        result_value = histo->Integral(histo->FindBin(MHTlow),histo->GetNbinsX());
    }
    else result_value = histo->Integral(histo->FindBin(MHTlow), histo->FindBin(MHTup)-1);

    return result_value;
}
////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////
double Prediction::GetResultError(TH1F* histo, double MHTlow, double MHTup)
{
    double result_error;
    if( MHTlow == MHTup ) {
        histo->IntegralAndError(histo->FindBin(MHTlow),histo->GetNbinsX(), result_error);
    }
    else histo->IntegralAndError(histo->FindBin(MHTlow), histo->FindBin(MHTup)-1, result_error);

    return result_error;
}
////////////////////////////////////////////////////////////////////////////////////////



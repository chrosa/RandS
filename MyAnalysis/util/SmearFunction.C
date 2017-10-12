#include "SmearFunction.h"

#include <TROOT.h>
#include <TKey.h>
#include <TFile.h>
#include <TMath.h>
#include <TArray.h>

#include <TStyle.h>
#include <TCanvas.h>

#include <memory>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

SmearFunction::SmearFunction(std::string smearingfile,
                             std::string inputhistPtHF, std::string inputhistEtaHF, std::string inputhistPhiHF, std::string inputhistResMuHF,
                             std::string inputhistPtLF, std::string inputhistEtaLF, std::string inputhistPhiLF, std::string inputhistResMuLF,
                             std::vector<double> PtBinEdges, std::vector<double> EtaBinEdges, std::vector<double> ResBinEdges,
                             std::vector<double> PtBinEdges_scaling, std::vector<double> EtaBinEdges_scaling,
                             std::vector<double> AdditionalSmearing,
                             std::vector<double> LowerTailScaling, std::vector<double> UpperTailScaling,
                             double AdditionalSmearing_variation, double LowerTailScaling_variation, double UpperTailScaling_variation,
                             bool absoluteTailScaling,
                             double A0RMS, double A1RMS, double probExtreme) {


    smearingfile_ = smearingfile;
    inputhistPtHF_ = inputhistPtHF;
    inputhistEtaHF_ = inputhistEtaHF;
    inputhistPhiHF_ = inputhistPhiHF;
    inputhistPtLF_ = inputhistPtLF;
    inputhistEtaLF_ = inputhistEtaLF;
    inputhistPhiLF_ = inputhistPhiLF;
    inputhistResMuHF_ = inputhistResMuHF;
    inputhistResMuLF_ = inputhistResMuLF;
    PtBinEdges_ = PtBinEdges;
    EtaBinEdges_ = EtaBinEdges;
    ResBinEdges_ = ResBinEdges;
    PtBinEdges_scaling_ = PtBinEdges_scaling;
    EtaBinEdges_scaling_ = EtaBinEdges_scaling;
    AdditionalSmearing_ = AdditionalSmearing;
    LowerTailScaling_ = LowerTailScaling;
    UpperTailScaling_ = UpperTailScaling;
    AdditionalSmearing_variation_ = AdditionalSmearing_variation;
    LowerTailScaling_variation_ = LowerTailScaling_variation;
    UpperTailScaling_variation_ = UpperTailScaling_variation;
    absoluteTailScaling_ = absoluteTailScaling;
    A0RMS_ = A0RMS;
    A1RMS_ = A1RMS;
    probExtreme_ = probExtreme;

    // Get correct dimensions for smear functions
    ResizeSmearFunctions();

    // Get/scale/fill smear functions
    CalculateSmearFunctions();
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
TH1F* SmearFunction::getSmearFunc(int i_flav, int i_eta, int i_pt) const {
    return smearFunc_scaled.at(i_flav).at(i_eta).at(i_pt);
}

TF1* SmearFunction::getSigmaPtForRebalancing(int i_eta) const {
    return  SigmaPt.at(0).at(i_eta);
}

TF1* SmearFunction::getSigmaPtScaledForRebalancing(int i_eta) const {
    return SigmaPt_scaled.at(0).at(i_eta);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void SmearFunction::CalculateSmearFunctions() {

    //// open root file/tree and create SmearingFunction histo
    TFile *f1 = new TFile(smearingfile_.c_str(), "READ", "", 0);

    //// Fetch counts for efficiency calculation
    for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
        char hname[100];
        sprintf(hname, "h_NjetGen_tot_Eta%i", i_eta);
        cout << hname << endl;
        TH1F* gen_tot = (TH1F*) f1->FindObjectAny(hname);
        sprintf(hname, "h_NjetGen_b_Eta%i", i_eta);
        cout << hname << endl;
        TH1F* gen_b = (TH1F*) f1->FindObjectAny(hname);
        sprintf(hname, "h_NjetGen_nob_Eta%i", i_eta);
        cout << hname << endl;
        TH1F* gen_nob = (TH1F*) f1->FindObjectAny(hname);

        sprintf(hname, "h_NjetReco_tot_Eta%i", i_eta);
        cout << hname << endl;
        TH1F* reco_tot = (TH1F*) f1->FindObjectAny(hname);
        sprintf(hname, "h_NjetReco_b_Eta%i", i_eta);
        cout << hname << endl;
        TH1F* reco_b = (TH1F*) f1->FindObjectAny(hname);
        sprintf(hname, "h_NjetReco_nob_Eta%i", i_eta);
        cout << hname << endl;
        TH1F* reco_nob = (TH1F*) f1->FindObjectAny(hname);

        TH1F* eff_tot = new TH1F(*reco_tot);
        eff_tot->Divide(gen_tot);
        TH1F* eff_b = new TH1F(*reco_b);
        eff_b->Divide(gen_b);
        TH1F* eff_nob = new TH1F(*reco_nob);
        eff_nob->Divide(gen_nob);

        RecoEff_tot.push_back(eff_tot);
        RecoEff_b.push_back(eff_b);
        RecoEff_nob.push_back(eff_nob);
    }

    //// Fetch histos of muon response and merge eta bins
    for (unsigned int i_Pt = 0; i_Pt < PtBinEdges_.size() - 1; ++i_Pt) {
        //// Get the histos
        char hname[100];
        int i_eta = 0;
        //// get no heavy flavor hists
        sprintf(hname, "%s_Pt%i_Eta%i", inputhistResMuLF_.c_str(), i_Pt, i_eta);
        cout << hname << endl;
        TH2F* tmp = (TH2F*) f1->FindObjectAny(hname)->Clone();
        i_eta = 1;
        sprintf(hname, "%s_Pt%i_Eta%i", inputhistResMuLF_.c_str(), i_Pt, i_eta);
        cout << hname << endl;
        tmp->Add((TH2F*) f1->FindObjectAny(hname));
        i_eta = 2;
        sprintf(hname, "%s_Pt%i_Eta%i", inputhistResMuLF_.c_str(), i_Pt, i_eta);
        cout << hname << endl;
        tmp->Add((TH2F*) f1->FindObjectAny(hname));
        for (unsigned int i_Res = 0; i_Res < ResBinEdges_.size() - 1; ++i_Res) {
			sprintf(hname, "%s_Pt%i_Res%i", inputhistResMuLF_.c_str(), i_Pt, i_Res);
            muRes.at(0).at(i_Pt).at(i_Res) = tmp->ProjectionY(hname, i_Res+1, i_Res+1);
        }

        //// get heavy flavor hists
        i_eta = 0;
        sprintf(hname, "%s_Pt%i_Eta%i", inputhistResMuHF_.c_str(), i_Pt, i_eta);
        cout << hname << endl;
        TH2F* tmp2 = (TH2F*) f1->FindObjectAny(hname)->Clone();
        i_eta = 1;
        sprintf(hname, "%s_Pt%i_Eta%i", inputhistResMuHF_.c_str(), i_Pt, i_eta);
        cout << hname << endl;
        tmp2->Add((TH2F*) f1->FindObjectAny(hname));
        i_eta = 2;
        sprintf(hname, "%s_Pt%i_Eta%i", inputhistResMuHF_.c_str(), i_Pt, i_eta);
        cout << hname << endl;
        tmp2->Add((TH2F*) f1->FindObjectAny(hname));
        for (unsigned int i_Res = 0; i_Res < ResBinEdges_.size() - 1; ++i_Res) {
			sprintf(hname, "%s_Pt%i_Res%i", inputhistResMuHF_.c_str(), i_Pt, i_Res);
            muRes.at(1).at(i_Pt).at(i_Res) = tmp2->ProjectionY(hname, i_Res+1, i_Res+1);
        }
    }

    //// Fetch histos and fit gaussian core
    for (unsigned int i_Pt = 0; i_Pt < PtBinEdges_.size() - 1; ++i_Pt) {
        for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
            //// Get the histos
            char hname[100];
            cout << "Histos, i_Pt: " << i_Pt <<  " i_eta: " << i_eta << endl;
            //// get no heavy flavor hists
            sprintf(hname, "%s_Pt%i_Eta%i", inputhistPtLF_.c_str(), i_Pt, i_eta);
            cout << hname << endl;
            smearFunc.at(0).at(i_eta).at(i_Pt) = (TH1F*) f1->FindObjectAny(hname);
            sprintf(hname, "%s_Pt%i_Eta%i", inputhistEtaLF_.c_str(), i_Pt, i_eta);
            smearFuncEta.at(0).at(i_eta).at(i_Pt) = (TH1F*) f1->FindObjectAny(hname);
            sprintf(hname, "%s_Pt%i_Eta%i", inputhistPhiLF_.c_str(), i_Pt, i_eta);
            smearFuncPhi.at(0).at(i_eta).at(i_Pt) = (TH1F*) f1->FindObjectAny(hname);
            //// get heavy flavor hists
            sprintf(hname, "%s_Pt%i_Eta%i", inputhistPtHF_.c_str(), i_Pt, i_eta);
            smearFunc.at(1).at(i_eta).at(i_Pt) = (TH1F*) f1->FindObjectAny(hname);
            sprintf(hname, "%s_Pt%i_Eta%i", inputhistEtaHF_.c_str(), i_Pt, i_eta);
            smearFuncEta.at(1).at(i_eta).at(i_Pt) = (TH1F*) f1->FindObjectAny(hname);
            sprintf(hname, "%s_Pt%i_Eta%i", inputhistPhiHF_.c_str(), i_Pt, i_eta);
            smearFuncPhi.at(1).at(i_eta).at(i_Pt) = (TH1F*) f1->FindObjectAny(hname);

            //// Phi resolution

            for (unsigned int i_flav = 0; i_flav < 2; ++i_flav) {
                //// Get RMS
                if (smearFuncPhi.at(i_flav).at(i_eta).at(i_Pt)->GetEntries() > 100) {
                    double RMS = smearFuncPhi.at(i_flav).at(i_eta).at(i_Pt)->GetRMS();
                    SigmaPhi.at(i_flav).at(i_eta).at(i_Pt) = RMS;
                } else {
                    double RMS = 0.02;
                    SigmaPhi.at(i_flav).at(i_eta).at(i_Pt) = RMS;
                }
            }

            //// Eta resolution

            for (unsigned int i_flav = 0; i_flav < 2; ++i_flav) {
                if (smearFuncEta.at(i_flav).at(i_eta).at(i_Pt)->GetEntries() > 100) {
                    double RMS = smearFuncEta.at(i_flav).at(i_eta).at(i_Pt)->GetRMS();
                    SigmaEta.at(i_flav).at(i_eta).at(i_Pt) = RMS;
                } else {
                    double RMS = 0.03;
                    SigmaEta.at(i_flav).at(i_eta).at(i_Pt) = RMS;
                }
            }

            //// Pt response

            for (unsigned int i_flav = 0; i_flav < 2; ++i_flav) {
                if (probExtreme_ > 0) {
                    double p = probExtreme_ * smearFunc.at(i_flav).at(i_eta).at(i_Pt)->Integral();
                    smearFunc.at(i_flav).at(i_eta).at(i_Pt)->SetBinContent(1, p);
                }
                //// Get width of gaussian core

                // check if bin is meaningfull (Pt(E_bin)>>20 GeV)
                // otherwise reco jet pT threshold biases jet response to larger values
                bool BinIsOK = (PtBinEdges_.at(i_Pt)/cosh(EtaBinEdges_.at(i_eta)) > 25.);
                //bool BinIsOK = true;
                if (smearFunc.at(i_flav).at(i_eta).at(i_Pt)->GetEntries() > 100 && BinIsOK) {
                    double RMS = smearFunc.at(i_flav).at(i_eta).at(i_Pt)->GetRMS();
                    double MEAN = smearFunc.at(i_flav).at(i_eta).at(i_Pt)->GetMean();
                    TF1* fitfunction = new TF1("f", "gaus(0)", MEAN - 1 * RMS, MEAN + 1 * RMS);
                    fitfunction->SetParameters(smearFunc.at(i_flav).at(i_eta).at(i_Pt)->GetMaximum(), MEAN, RMS);
                    smearFunc.at(i_flav).at(i_eta).at(i_Pt)->Fit(fitfunction, "LLRQN");
                    double Pt = (PtBinEdges_.at(i_Pt) + PtBinEdges_.at(i_Pt + 1)) / 2;
                    double eta = (EtaBinEdges_.at(i_eta) + EtaBinEdges_.at(i_eta + 1)) / 2;
                    double f = GetAdditionalSmearing(Pt, eta);
                    SigmaPtHist.at(i_flav).at(i_eta)->SetBinContent(i_Pt + 1, std::abs(fitfunction->GetParameter(2)));
                    SigmaPtHist.at(i_flav).at(i_eta)->SetBinError(i_Pt + 1, fitfunction->GetParError(2));
                    SigmaPtHist_scaled.at(i_flav).at(i_eta)->SetBinContent(i_Pt+1, std::abs(fitfunction->GetParameter(2))*f);
                    SigmaPtHist_scaled.at(i_flav).at(i_eta)->SetBinError(i_Pt+1, fitfunction->GetParError(2) * f);

                    //// Split smearFunc in Core and Tail
                    TH1F* hResponseFit = new TH1F(*smearFunc.at(i_flav).at(i_eta).at(i_Pt));
                    hResponseFit->Reset();
                    for (int i = 0; i < hResponseFit->GetNbinsX(); ++i) {
                        hResponseFit->SetBinContent(i, fitfunction->Eval(hResponseFit->GetBinCenter(i)));
                    }

                    //// Split lower tail
                    smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt) = new TH1F(*smearFunc.at(i_flav).at(i_eta).at(i_Pt));
                    smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt)->Add(hResponseFit, -1.);
                    for (int i = 0; i <= smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt)->GetNbinsX(); ++i) {
                        double tmp_v = smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt)->GetBinContent(i);
                        double tmp_e = smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt)->GetBinError(i);
                        //// by definition a tail has positive entries
                        if (tmp_v < 0) {
                            tmp_v = 0;
                            tmp_e = 0;
                        } else {
                            //// suppress everything except for low response tail
                            double scale = 1;
                            double x = smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt)->GetBinCenter(i);
                            if (x > MEAN - 1 * RMS)
                                scale = 0.;
                            tmp_v = scale * smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt)->GetBinContent(i);
                            tmp_e = scale * smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt)->GetBinError(i);
                        }
                        smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt)->SetBinContent(i, tmp_v);
                        smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt)->SetBinError(i, tmp_e);
                    }

                    //// Split upper tail
                    smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt) = new TH1F(*smearFunc.at(i_flav).at(i_eta).at(i_Pt));
                    smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt)->Add(hResponseFit, -1.);
                    for (int i = 0; i <= smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt)->GetNbinsX(); ++i) {
                        double tmp_v = smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt)->GetBinContent(i);
                        double tmp_e = smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt)->GetBinError(i);
                        //// by definition a tail has positive entries
                        if (tmp_v < 0) {
                            tmp_v = 0;
                            tmp_e = 0;
                        } else {
                            //// suppress everything except for high response tail
                            double scale = 1;
                            double x = smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt)->GetBinCenter(i);
                            if (x < MEAN + 1 * RMS)
                                scale = 0.;
                            tmp_v = scale * smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt)->GetBinContent(i);
                            tmp_e = scale * smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt)->GetBinError(i);
                        }
                        smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt)->SetBinContent(i, tmp_v);
                        smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt)->SetBinError(i, tmp_e);
                    }

                    smearFunc_Core.at(i_flav).at(i_eta).at(i_Pt) = new TH1F(*smearFunc.at(i_flav).at(i_eta).at(i_Pt));
                    smearFunc_Core.at(i_flav).at(i_eta).at(i_Pt)->Add(smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt), -1.);
                    smearFunc_Core.at(i_flav).at(i_eta).at(i_Pt)->Add(smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt), -1.);

                } else {
                    //// Set core and tail if needed
                    smearFunc_Core.at(i_flav).at(i_eta).at(i_Pt) = new TH1F(*smearFunc.at(i_flav).at(i_eta).at(i_Pt));
                    smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt) = new TH1F(*smearFunc.at(i_flav).at(i_eta).at(i_Pt));
                    smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt)->Reset();
                    smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt) = new TH1F(*smearFunc.at(i_flav).at(i_eta).at(i_Pt));
                    smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt)->Reset();
                }
            }
        }
    }

    //// Fit scaled gaussian sigma as function of pt
    for (unsigned int i_flav = 0; i_flav < 2; ++i_flav) {
        for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
            char fname[100];
            sprintf(fname, "SigmaPtScaled_JetFlavor%i_Eta%i", i_flav, i_eta);
            bool first = false;
            int FirstBin = 1;
            int LastBin = SigmaPtHist_scaled.at(i_flav).at(i_eta)->GetNbinsX();
            for (int j = 1; j <= SigmaPtHist_scaled.at(i_flav).at(i_eta)->GetNbinsX(); ++j) {
                if (!first && SigmaPtHist_scaled.at(i_flav).at(i_eta)->GetBinContent(j) > 0.001) {
                    first = true;
                    FirstBin = j;
                }
                if (first && SigmaPtHist_scaled.at(i_flav).at(i_eta)->GetBinContent(j) < 0.001) {
                    LastBin = j - 1;
                    break;
                }
            }
            SigmaPt_scaled.at(i_flav).at(i_eta) = new TF1(fname,
                    "sqrt(sign(1,[0])*pow([0]/x,2)+pow([1],2)*pow(x,[2]-1.)+pow([3],2))",
                    SigmaPtHist_scaled.at(i_flav).at(i_eta)->GetBinCenter(FirstBin),
                    SigmaPtHist_scaled.at(i_flav).at(i_eta)->GetBinCenter(LastBin));
            SigmaPt_scaled.at(i_flav).at(i_eta)->SetParameters(1.2, 0., 0.03);
            SigmaPtHist_scaled.at(i_flav).at(i_eta)->Fit(SigmaPt_scaled.at(i_flav).at(i_eta), "LLRQ");
        }
    }

    //// Fit gaussian sigma as function of pt
    for (unsigned int i_flav = 0; i_flav < 2; ++i_flav) {
        for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
            char fname[100];
            sprintf(fname, "SigmaPt_JetFlavor%i_Eta%i", i_flav, i_eta);
            bool first = false;
            int FirstBin = 1;
            int LastBin = SigmaPtHist.at(i_flav).at(i_eta)->GetNbinsX();
            for (int j = 1; j <= SigmaPtHist.at(i_flav).at(i_eta)->GetNbinsX(); ++j) {
                if (!first && SigmaPtHist.at(i_flav).at(i_eta)->GetBinContent(j) > 0.001) {
                    first = true;
                    FirstBin = j;
                }
                if (first && SigmaPtHist.at(i_flav).at(i_eta)->GetBinContent(j) < 0.001) {
                    LastBin = j - 1;
                    break;
                }
            }
            SigmaPt.at(i_flav).at(i_eta) = new TF1(fname,"sqrt(sign(1,[0])*pow([0]/x,2)+pow([1],2)*pow(x,[2]-1.)+pow([3],2))",
                                                   SigmaPtHist.at(i_flav).at(i_eta)->GetBinCenter(FirstBin),
                                                   SigmaPtHist.at(i_flav).at(i_eta)->GetBinCenter(LastBin));
            SigmaPt.at(i_flav).at(i_eta)->SetParameters(1.2, 0., 0.03);
            SigmaPtHist.at(i_flav).at(i_eta)->Fit(SigmaPt.at(i_flav).at(i_eta), "LLRQ");
        }
    }

    //// Book and fill histograms for smeared and scaled response functions
    for (unsigned int i_flav = 0; i_flav < 2; ++i_flav) {
        for (unsigned int i_Pt = 0; i_Pt < PtBinEdges_.size() - 1; ++i_Pt) {
            for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
                char hname[100];
                sprintf(hname, "SmearedAndScaledResolution_Pt%i_Eta%i_JetFlavor%i", i_Pt, i_eta, i_flav);
                smearFunc_scaled.at(i_flav).at(i_eta).at(i_Pt) = new TH1F(hname, hname,
                        smearFunc.at(i_flav).at(i_eta).at(i_Pt)->GetNbinsX(),
                        smearFunc.at(i_flav).at(i_eta).at(i_Pt)->GetXaxis()->GetXmin(),
                        smearFunc.at(i_flav).at(i_eta).at(i_Pt)->GetXaxis()->GetXmax());

                bool BinIsOK = (PtBinEdges_.at(i_Pt)/cosh(EtaBinEdges_.at(i_eta)) > 25.);
                //bool BinIsOK = true;
                if (smearFunc.at(i_flav).at(i_eta).at(i_Pt)->GetEntries() > 100 && BinIsOK) {
                    //// fold core and tail with additional gaussian
                    TH1F smearFunc_Core_tmp(*smearFunc_Core.at(i_flav).at(i_eta).at(i_Pt));
                    TH1F smearFunc_LowerTail_tmp(*smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt));
                    TH1F smearFunc_UpperTail_tmp(*smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt));
                    double AddSmear = GetAdditionalSmearing(PtBinEdges_.at(i_Pt), EtaBinEdges_.at(i_eta));
                    if (AddSmear > 1) {
                        //// additional sigma from (1+x)*sigma = sqrt(sigma^2+add_sigma^2)
                        //// or from sigma' = sqrt((sigma'/(1+x))^2+add_sigma^2)
                        double sigma = SigmaPtHist_scaled.at(i_flav).at(i_eta)->GetBinContent(i_Pt + 1);
                        //// if no sigma was fitted use the extrapolation
                        if (sigma == 0)
                            sigma = SigmaPt_scaled.at(i_flav).at(i_eta)->Eval(SigmaPtHist_scaled.at(i_flav).at(i_eta)->GetBinCenter(i_Pt + 1));
                        double AdditionalSigma = TMath::Sqrt(1. - 1. / pow(AddSmear, 2)) * sigma;
                        smearFunc_Core_tmp.Reset();
                        smearFunc_LowerTail_tmp.Reset();
                        smearFunc_UpperTail_tmp.Reset();
                        FoldWithGaussian(*smearFunc_Core.at(i_flav).at(i_eta).at(i_Pt), smearFunc_Core_tmp, AdditionalSigma);
                        FoldWithGaussian(*smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt), smearFunc_LowerTail_tmp, AdditionalSigma);
                        FoldWithGaussian(*smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt), smearFunc_UpperTail_tmp, AdditionalSigma);
                    } else if (AddSmear < 1) {
                        smearFunc_Core_tmp.Reset();
                        StretchHisto(*smearFunc_Core.at(i_flav).at(i_eta).at(i_Pt), smearFunc_Core_tmp, AddSmear);
                    }

                    //// Scale tails
                    double LowerTailScale = GetLowerTailScaling(PtBinEdges_.at(i_Pt), EtaBinEdges_.at(i_eta));
                    double UpperTailScale = GetUpperTailScaling(PtBinEdges_.at(i_Pt), EtaBinEdges_.at(i_eta));
                    //cout << "absolute scaling factor: " << TailScale << endl;
                    if (absoluteTailScaling_) {
                        double RMS = smearFunc.at(i_flav).at(i_eta).at(i_Pt)->GetRMS();
                        //cout << "Integral from " << 1-A1RMS_*RMS << " to " << 1-A0RMS_*RMS << endl;

                        //// get integral of tails
                        int i_min = smearFunc.at(i_flav).at(i_eta).at(i_Pt)->FindBin(1 - A1RMS_ * RMS);
                        int i_max = smearFunc.at(i_flav).at(i_eta).at(i_Pt)->FindBin(1 - A0RMS_ * RMS);
                        double RLowerTail = smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt)->Integral(i_min, i_max);
                        double Rcore = smearFunc_Core.at(i_flav).at(i_eta).at(i_Pt)->Integral(i_min, i_max);
                        if (RLowerTail > 0)
                            LowerTailScale = (LowerTailScale * (RLowerTail + Rcore) - Rcore) / RLowerTail;

                        i_min = smearFunc.at(i_flav).at(i_eta).at(i_Pt)->FindBin(1 + A0RMS_ * RMS);
                        i_max = smearFunc.at(i_flav).at(i_eta).at(i_Pt)->FindBin(1 + A1RMS_ * RMS);
                        double RUpperTail = smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt)->Integral(i_min, i_max);
                        Rcore = smearFunc_Core.at(i_flav).at(i_eta).at(i_Pt)->Integral(i_min, i_max);
                        if (RUpperTail > 0)
                            UpperTailScale = (UpperTailScale * (RUpperTail + Rcore) - Rcore) / RUpperTail;

                    }
                    //cout << "tail scaling factor: " << TailScale << endl;
                    smearFunc_scaled.at(i_flav).at(i_eta).at(i_Pt)->Add(&smearFunc_Core_tmp);
                    smearFunc_scaled.at(i_flav).at(i_eta).at(i_Pt)->Add(&smearFunc_LowerTail_tmp, LowerTailScale);
                    smearFunc_scaled.at(i_flav).at(i_eta).at(i_Pt)->Add(&smearFunc_UpperTail_tmp, UpperTailScale);
                } else {
                    //// Replace Histograms with only few entries by gaussians
                    double N = 1;
                    if (smearFunc.at(i_flav).at(i_eta).at(i_Pt)->GetEntries() > 10) {
                        N = smearFunc.at(i_flav).at(i_eta).at(i_Pt)->Integral();
                    }
                    cout << "Too few entries for (i_Pt, i_eta, i_flav): " << i_Pt << ", " << i_eta << ", " << i_flav
                         << ", entries = " << smearFunc.at(i_flav).at(i_eta).at(i_Pt)->GetEntries() << endl;
                    for (int j = 1; j <= smearFunc_scaled.at(i_flav).at(i_eta).at(i_Pt)->GetNbinsX(); ++j) {
                        double pt = (PtBinEdges_.at(i_Pt) + PtBinEdges_.at(i_Pt + 1)) / 2;
                        double g = N * smearFunc_scaled.at(0).at(i_eta).at(i_Pt)->GetBinWidth(j) * TMath::Gaus(smearFunc_scaled.at(0).at(i_eta).at(i_Pt)->GetBinCenter(j), 1., SigmaPt_scaled.at(0).at(i_eta)->Eval(pt), true);
                        smearFunc_scaled.at(i_flav).at(i_eta).at(i_Pt)->SetBinContent(j, g);
                        smearFunc_Core.at(i_flav).at(i_eta).at(i_Pt)->SetBinContent(j, g);
                        smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt)->SetBinContent(j, 0.);
                        smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt)->SetBinContent(j, 0.);
                    }
                }
            }
        }
    }

    gROOT->Reset();
    gROOT->SetStyle("Plain");
    gStyle->SetStatColor(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetOptStat(0);
    gStyle->SetStatBorderSize(2);
    gStyle->SetOptTitle(1);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadBorderSize(2);
    gStyle->SetPalette(51, 0);
    gStyle->SetPadBottomMargin(0.25);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadLeftMargin(0.2);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetTitleOffset(1.2, "X");
    gStyle->SetTitleOffset(1.6, "Y");
    gStyle->SetTitleOffset(1.0, "Z");
    gStyle->SetLabelSize(0.05, "X");
    gStyle->SetLabelSize(0.05, "Y");
    gStyle->SetLabelSize(0.05, "Z");
    gStyle->SetLabelOffset(0.02, "X");
    gStyle->SetLabelOffset(0.02, "Y");
    gStyle->SetLabelOffset(0.02, "Z");
    gStyle->SetTitleSize(0.05, "X");
    gStyle->SetTitleSize(0.05, "Y");
    gStyle->SetTitleSize(0.05, "Z");
    gStyle->SetTitleColor(1);
    gStyle->SetTitleFillColor(0);
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetTitleY(0.99);
    gStyle->SetTitleX(0.15);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetLineWidth(2);
    gStyle->SetHistLineWidth(2);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetNdivisions(505, "X");
    gStyle->SetMarkerSize(0.8);
    gStyle->SetTickLength(0.03);
    gROOT->ForceStyle();
    TString psfile = "SmearFunctions";
    TCanvas *c = new TCanvas("", "", 800, 800);
    c->cd();
    c->Print(psfile + ".pdf(");
    for (unsigned int i_flav = 0; i_flav < 2; ++i_flav) {
        for (unsigned int i_Pt = 0; i_Pt < PtBinEdges_.size() - 1; ++i_Pt) {
            for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
                char cname[100];
                sprintf(cname, "c_Pt%i_Eta%i_JetFlavor%i", i_Pt, i_eta, i_flav);
                c->SetName(cname);
                c->SetLogy();
                smearFunc.at(i_flav).at(i_eta).at(i_Pt)->SetTitle(cname);
                smearFunc.at(i_flav).at(i_eta).at(i_Pt)->SetLineColor(kBlack);
                smearFunc_Core.at(i_flav).at(i_eta).at(i_Pt)->SetLineColor(kGreen);
                smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt)->SetLineColor(kRed);
                smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt)->SetLineColor(kMagenta);
                smearFunc_scaled.at(i_flav).at(i_eta).at(i_Pt)->SetLineColor(kBlue);
                smearFunc.at(i_flav).at(i_eta).at(i_Pt)->Draw("hist");
                smearFunc_Core.at(i_flav).at(i_eta).at(i_Pt)->Draw("hist same");
                smearFunc_LowerTail.at(i_flav).at(i_eta).at(i_Pt)->Draw("hist same");
                smearFunc_UpperTail.at(i_flav).at(i_eta).at(i_Pt)->Draw("hist same");
                smearFunc_scaled.at(i_flav).at(i_eta).at(i_Pt)->Draw("hist same");
                c->Print(psfile + ".pdf");
            }
        }
    }
    for (unsigned int i_flav = 0; i_flav < 2; ++i_flav) {
        for (unsigned int i_Pt = 0; i_Pt < PtBinEdges_.size() - 1; ++i_Pt) {
            for (unsigned int i_res = 0; i_res < ResBinEdges_.size() - 1; ++i_res) {
                char cname[100];
                sprintf(cname, "c_Pt%i_Res%i_JetFlavor%i", i_Pt, i_res, i_flav);
                c->SetName(cname);
                c->SetLogy();
                muRes.at(i_flav).at(i_Pt).at(i_res)->SetTitle(cname);
                muRes.at(i_flav).at(i_Pt).at(i_res)->SetLineColor(kBlack);
                muRes.at(i_flav).at(i_Pt).at(i_res)->Draw("hist");
                c->Print(psfile + ".pdf");
            }
        }
    }
    for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
        char cname[100];
        sprintf(cname, "c_Eff_Eta%i", i_eta);
        c->SetName(cname);
        RecoEff_nob.at(i_eta)->SetTitle(cname);
        RecoEff_nob.at(i_eta)->SetLineColor(kBlue);
        RecoEff_nob.at(i_eta)->Draw("histe");
        //RecoEff_b.at(i_eta)->SetTitle(cname);
        //RecoEff_b.at(i_eta)->SetLineColor(kRed);
        //RecoEff_b.at(i_eta)->Draw("hist same");
        c->Print(psfile + ".pdf");

        sprintf(cname, "c_sigmaPt_Eta%i", i_eta);
        c->SetName(cname);
        SigmaPtHist.at(0).at(i_eta)->SetTitle(cname);
        SigmaPtHist.at(0).at(i_eta)->SetLineColor(kBlue);
        SigmaPtHist.at(0).at(i_eta)->Draw("histe");
        //SigmaPtHist.at(1).at(i_eta)->SetTitle(cname);
        //SigmaPtHist.at(1).at(i_eta)->SetLineColor(kRed);
        //SigmaPtHist.at(1).at(i_eta)->Draw("hist same");
        c->Print(psfile + ".pdf");

        sprintf(cname, "c_sigmaPtscaled_Eta%i", i_eta);
        c->SetName(cname);
        SigmaPtHist_scaled.at(0).at(i_eta)->SetTitle(cname);
        SigmaPtHist_scaled.at(0).at(i_eta)->SetLineColor(kBlue);
        SigmaPtHist_scaled.at(0).at(i_eta)->Draw("histe");
        //SigmaPtHist_scaled.at(1).at(i_eta)->SetTitle(cname);
        //SigmaPtHist_scaled.at(1).at(i_eta)->SetLineColor(kRed);
        //SigmaPtHist_scaled.at(1).at(i_eta)->Draw("hist same");
        c->Print(psfile + ".pdf");
    }
    c->Print(psfile + ".pdf)");

}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double SmearFunction::GetAdditionalSmearing(const double& pt, const double& eta) {
    int i_Pt = GetIndex(pt, &PtBinEdges_scaling_);
    int i_eta = GetIndex(eta, &EtaBinEdges_scaling_);
    double result = AdditionalSmearing_.at(i_eta * (PtBinEdges_scaling_.size() - 1) + i_Pt) * AdditionalSmearing_variation_;
    return result;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double SmearFunction::GetLowerTailScaling(const double& pt, const double& eta) {
    int i_Pt = GetIndex(pt, &PtBinEdges_scaling_);
    int i_eta = GetIndex(eta, &EtaBinEdges_scaling_);
    double result = LowerTailScaling_.at(i_eta * (PtBinEdges_scaling_.size() - 1) + i_Pt) * LowerTailScaling_variation_;
    return result;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double SmearFunction::GetUpperTailScaling(const double& pt, const double& eta) {
    int i_Pt = GetIndex(pt, &PtBinEdges_scaling_);
    int i_eta = GetIndex(eta, &EtaBinEdges_scaling_);
    double result = UpperTailScaling_.at(i_eta * (PtBinEdges_scaling_.size() - 1) + i_Pt) * UpperTailScaling_variation_;
    return result;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// Fold histogram input with gaussian resolution
void SmearFunction::FoldWithGaussian(const TH1& input, TH1& output, const double& sigma) {

    double min = input.GetXaxis()->GetXmin();
    double max = input.GetXaxis()->GetXmax();
    for (int i = 0; i < input.GetNbinsX(); ++i) {
        double weight = input.GetBinContent(i);
        double mean = input.GetBinCenter(i);
        TF1 gauss("gauss", "gaus(0)", min, max);
        gauss.SetParameters(weight * 1. / sigma / sqrt(2 * TMath::Pi()), mean, sigma);
        for (int j = 0; j < output.GetNbinsX(); ++j) {
            double xmin = output.GetBinLowEdge(j);
            double xmax = output.GetBinLowEdge(j) + output.GetBinWidth(j);
            output.AddBinContent(j, gauss.Integral(xmin, xmax));
        }
    }

}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// Fold histogram input with gaussian resolution
void SmearFunction::StretchHisto(const TH1& input, TH1& output, const double& f) {

    if (input.Integral() > 0) {
        double mean = input.GetMean();
        for (int i = 0; i < 1000000; ++i) {
            double r = input.GetRandom();
            double rprime = mean + (r - mean) * f;
            output.Fill(rprime);
        }
        output.Scale(input.Integral() / output.Integral());
    }

}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
int SmearFunction::GetIndex(const double& x, const std::vector<double>* vec) {
    int index = -1;
    for (std::vector<double>::const_iterator it = vec->begin(); it != vec->end(); ++it) {
        if ((*it) > fabs(x))
            break;
        ++index;
    }
    if (index < 0)
        index = 0;
    if (index > (int) vec->size() - 2)
        index = vec->size() - 2;

    return index;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void SmearFunction::ResizeSmearFunctions() {

    cout << "EtaBinEdges: ";
    for (std::vector<double>:: iterator it = EtaBinEdges_.begin(); it != EtaBinEdges_.end(); ++it) {
        cout << *it << ", ";
    }
    cout << endl;
    cout << "PtBinEdges: ";
    for (std::vector<double>:: iterator it = PtBinEdges_.begin(); it != PtBinEdges_.end(); ++it) {
        cout << *it << ", ";
    }
    cout << endl;
    cout << "ResBinEdges: ";
    for (std::vector<double>:: iterator it = ResBinEdges_.begin(); it != ResBinEdges_.end(); ++it) {
        cout << *it << ", ";
    }
    cout << endl;

    cout << "smearFunc" << endl;
    smearFunc.resize(2); //// two bins for jet flavour
    for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc.begin(); it != smearFunc.end(); ++it) {
        it->resize(EtaBinEdges_.size() - 1);
        for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
            jt->resize(PtBinEdges_.size() - 1);
        }
    }

    cout << "smearFuncEta" << endl;
    smearFuncEta.resize(2); //// two bins for jet flavour
    for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFuncEta.begin(); it != smearFuncEta.end(); ++it) {
        it->resize(EtaBinEdges_.size() - 1);
        for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
            jt->resize(PtBinEdges_.size() - 1);
        }
    }

    cout << "smearFuncPhi" << endl;
    smearFuncPhi.resize(2); //// two bins for jet flavour
    for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFuncPhi.begin(); it != smearFuncPhi.end(); ++it) {
        it->resize(EtaBinEdges_.size() - 1);
        for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
            jt->resize(PtBinEdges_.size() - 1);
        }
    }

    cout << "smearFunc_Core" << endl;
    smearFunc_Core.resize(2); //// two bins for jet flavour
    for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_Core.begin(); it != smearFunc_Core.end(); ++it) {
        it->resize(EtaBinEdges_.size() - 1);
        for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
            jt->resize(PtBinEdges_.size() - 1);
        }
    }

    cout << "smearFunc_LowerTail" << endl;
    smearFunc_LowerTail.resize(2); //// two bins for jet flavour
    for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_LowerTail.begin(); it != smearFunc_LowerTail.end(); ++it) {
        it->resize(EtaBinEdges_.size() - 1);
        for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
            jt->resize(PtBinEdges_.size() - 1);
        }
    }

    cout << "smearFunc_UpperTail" << endl;
    smearFunc_UpperTail.resize(2); //// two bins for jet flavour
    for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_UpperTail.begin(); it != smearFunc_UpperTail.end(); ++it) {
        it->resize(EtaBinEdges_.size() - 1);
        for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
            jt->resize(PtBinEdges_.size() - 1);
        }
    }

    cout << "smearFunc_scaled" << endl;
    smearFunc_scaled.resize(2); //// two bins for jet flavour
    for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_scaled.begin(); it != smearFunc_scaled.end(); ++it) {
        it->resize(EtaBinEdges_.size() - 1);
        for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
            jt->resize(PtBinEdges_.size() - 1);
        }
    }

    cout << "SigmaEta" << endl;
    SigmaEta.resize(2);
    for (std::vector<std::vector<std::vector<double> > >::iterator it = SigmaEta.begin(); it != SigmaEta.end(); ++it) {
        it->resize(EtaBinEdges_.size() - 1);
        for (std::vector<std::vector<double> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
            jt->resize(PtBinEdges_.size() - 1);
        }
    }

    cout << "SigmaPhi" << endl;
    SigmaPhi.resize(2);
    for (std::vector<std::vector<std::vector<double> > >::iterator it = SigmaPhi.begin(); it != SigmaPhi.end(); ++it) {
        it->resize(EtaBinEdges_.size() - 1);
        for (std::vector< std::vector<double> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
            jt->resize(PtBinEdges_.size() - 1);
        }
    }

    cout << "SigmaPtHist" << endl;
    SigmaPtHist.resize(2);
    for (unsigned int i_flav = 0; i_flav < 2; ++i_flav) {
        SigmaPtHist.at(i_flav).resize(EtaBinEdges_.size() - 1);
        for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
            char hname[100];
            sprintf(hname, "SigmaPtHist_JetFlavor%i_Eta%i", i_flav, i_eta);
            SigmaPtHist.at(i_flav).at(i_eta) = new TH1F(hname, hname, PtBinEdges_.size() - 1, &(PtBinEdges_.at(0)));
            SigmaPtHist.at(i_flav).at(i_eta)->Sumw2();
        }
    }

    cout << "SigmaPtHist_scaled" << endl;
    SigmaPtHist_scaled.resize(2);
    for (unsigned int i_flav = 0; i_flav < 2; ++i_flav) {
        SigmaPtHist_scaled.at(i_flav).resize(EtaBinEdges_.size() - 1);
        for (unsigned int i_eta = 0; i_eta < EtaBinEdges_.size() - 1; ++i_eta) {
            char hname[100];
            sprintf(hname, "SigmaPtHist_scaled_JetFlavor%i_Eta%i", i_flav, i_eta);
            SigmaPtHist_scaled.at(i_flav).at(i_eta) = new TH1F(hname, hname, PtBinEdges_.size() - 1, &(PtBinEdges_.at(0)));
            SigmaPtHist_scaled.at(i_flav).at(i_eta)->Sumw2();
        }
    }

    cout << "SigmaPt" << endl;
    SigmaPt.resize(2);
    for (std::vector<std::vector<TF1*> >::iterator it = SigmaPt.begin(); it != SigmaPt.end(); ++it) {
        it->resize(EtaBinEdges_.size() - 1);
    }

    cout << "SigmaPt_scaled" << endl;
    SigmaPt_scaled.resize(2);
    for (std::vector<std::vector<TF1*> >::iterator it = SigmaPt_scaled.begin(); it != SigmaPt_scaled.end(); ++it) {
        it->resize(EtaBinEdges_.size() - 1);
    }

    cout << "MuRes" << endl;
    muRes.resize(2); //// two bins for jet flavour
    for (std::vector<std::vector<std::vector<TH1D*> > >::iterator it = muRes.begin(); it != muRes.end(); ++it) {
        it->resize(PtBinEdges_.size() - 1);
        for (std::vector<std::vector<TH1D*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
            jt->resize(ResBinEdges_.size() - 1);
        }
    }

}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
SmearFunction::~SmearFunction() {

    for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc.begin(); it != smearFunc.end(); ++it) {
        for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
            for (std::vector<TH1F*>::iterator kt = jt->begin(); kt != jt->end(); ++kt) {
                delete *kt;
            }
        }
    }
    smearFunc.clear();

    for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_Core.begin(); it != smearFunc_Core.end(); ++it) {
        for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
            for (std::vector<TH1F*>::iterator kt = jt->begin(); kt != jt->end(); ++kt) {
                delete *kt;
            }
        }
    }
    smearFunc_Core.clear();

    for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_LowerTail.begin(); it != smearFunc_LowerTail.end(); ++it) {
        for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
            for (std::vector<TH1F*>::iterator kt = jt->begin(); kt != jt->end(); ++kt) {
                delete *kt;
            }
        }
    }
    smearFunc_LowerTail.clear();

    for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_UpperTail.begin(); it != smearFunc_UpperTail.end(); ++it) {
        for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
            for (std::vector<TH1F*>::iterator kt = jt->begin(); kt != jt->end(); ++kt) {
                delete *kt;
            }
        }
    }
    smearFunc_UpperTail.clear();

    for (std::vector<std::vector<std::vector<TH1F*> > >::iterator it = smearFunc_scaled.begin(); it != smearFunc_scaled.end(); ++it) {
        for (std::vector<std::vector<TH1F*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
            for (std::vector<TH1F*>::iterator kt = jt->begin(); kt != jt->end(); ++kt) {
                delete *kt;
            }
        }
    }
    smearFunc_scaled.clear();

    for (std::vector<std::vector<std::vector<TH1D*> > >::iterator it = muRes.begin(); it != muRes.end(); ++it) {
        for (std::vector<std::vector<TH1D*> >::iterator jt = it->begin(); jt != it->end(); ++jt) {
            for (std::vector<TH1D*>::iterator kt = jt->begin(); kt != jt->end(); ++kt) {
                delete *kt;
            }
        }
    }
    muRes.clear();

}
//--------------------------------------------------------------------------



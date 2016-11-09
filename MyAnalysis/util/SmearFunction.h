#ifndef SMEAR_FUNCTION_H
#define SMEAR_FUNCTION_H

#include "TF1.h"
#include "TH1F.h"

#include <string>
#include <vector>

using namespace std;

class SmearFunction {

    public:
        SmearFunction(std::string smearingfile_,
                      std::string inputhistPtHF_, std::string inputhistEtaHF_, std::string inputhistPhiHF_,
                      std::string inputhistPtLF_, std::string inputhistEtaLF_, std::string inputhistPhiLF_,
                      std::vector<double> PtBinEdges_, std::vector<double> EtaBinEdges_, std::vector<double> PtBinEdges_scaling_, std::vector<double> EtaBinEdges_scaling_,
                      std::vector<double> AdditionalSmearing_,
                      std::vector<double> LowerTailScaling_, std::vector<double> UpperTailScaling_,
                      double AdditionalSmearing_variation_, double LowerTailScaling_variation_, double UpperTailScaling_variation_,
                      bool absoluteTailScaling_,
                      double A0RMS_, double A1RMS_, double probExtreme_);
        ~SmearFunction();

        TF1* getSigmaPtForRebalancing(int i_eta) const;
        TF1* getSigmaPtScaledForRebalancing(int i_eta) const;
        TH1F* getSmearFunc(int i_flav, int i_eta, int i_pt) const;

        std::vector<std::vector<std::vector<double> > > SigmaEta;
        std::vector<std::vector<std::vector<double> > > SigmaPhi;

        std::vector<TH1F*> RecoEff_tot;
        std::vector<TH1F*> RecoEff_b;
        std::vector<TH1F*> RecoEff_nob;


    private:
        typedef std::vector<std::string>::const_iterator StrIter;

        void FillSigmaHistsForRebalancing();
        void ResizeSmearFunctions();
        void CalculateSmearFunctions();
        double GetAdditionalSmearing(const double&, const double&);
        double GetLowerTailScaling(const double&, const double&);
        double GetUpperTailScaling(const double&, const double&);
        void FoldWithGaussian(const TH1&, TH1&, const double&);
        void StretchHisto(const TH1&, TH1&, const double&);
        int GetIndex(const double&, const std::vector<double>*);

        std::string smearingfile_;
        std::string inputhistPtHF_;
        std::string inputhistEtaHF_;
        std::string inputhistPhiHF_;
        std::string inputhistPtLF_;
        std::string inputhistEtaLF_;
        std::string inputhistPhiLF_;

        std::vector<double> PtBinEdges_;
        std::vector<double> EtaBinEdges_;
        std::vector<double> PtBinEdges_scaling_;
        std::vector<double> EtaBinEdges_scaling_;
        std::vector<double> AdditionalSmearing_;
        std::vector<double> LowerTailScaling_;
        std::vector<double> UpperTailScaling_;
        double AdditionalSmearing_variation_;
        double LowerTailScaling_variation_;
        double UpperTailScaling_variation_;
        bool absoluteTailScaling_;

        double A0RMS_;
        double A1RMS_;
        double probExtreme_;

        std::string uncertaintyName_;

        //// vectors of response functions
        std::vector<std::vector<std::vector<TH1F*> > > smearFunc;
        std::vector<std::vector<std::vector<TH1F*> > > smearFuncEta;
        std::vector<std::vector<std::vector<TH1F*> > > smearFuncPhi;
        std::vector<std::vector<std::vector<TH1F*> > > smearFunc_Core;
        std::vector<std::vector<std::vector<TH1F*> > > smearFunc_LowerTail;
        std::vector<std::vector<std::vector<TH1F*> > > smearFunc_UpperTail;
        std::vector<std::vector<std::vector<TH1F*> > > smearFunc_scaled;
        std::vector<std::vector<std::vector<TH1F*> > > smearFunc_total;

        std::vector<std::vector<TH1F*> > SigmaPtHist;
        std::vector<std::vector<TH1F*> > SigmaPtHist_scaled;

        std::vector<std::vector<TF1*> > SigmaPt;
        std::vector<std::vector<TF1*> > SigmaPt_scaled;

        std::vector<TH1F*> NReco_tot;
        std::vector<TH1F*> NReco_b;
        std::vector<TH1F*> NReco_nob;

        std::vector<TH1F*> NGen_tot;
        std::vector<TH1F*> NGen_b;
        std::vector<TH1F*> NGen_nob;

};

#endif

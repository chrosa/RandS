#ifndef MyAnalysis_RandS_H
#define MyAnalysis_RandS_H

#include "MyAnalysis/SmearFunction.h"

#include <EventLoop/Algorithm.h>
#include "TRandom3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TMatrixD.h"

#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include "JetSelectorTools/JetCleaningTool.h"
#include "SUSYTools/SUSYToolsDict.h"
#include "SUSYTools/ISUSYObjDef_xAODTool.h"


#include <memory>



class RandS : public EL::Algorithm
{
    // put your configuration variables here as public variables.
    // that way they can be set directly from CINT and python.

    public:

        int m_eventCounter; //!
        int m_numCleanEvents; //!
        GoodRunsListSelectionTool *m_grl; //!
        JetCleaningTool *m_jetCleaning; //!
   		ST::SUSYObjDef_xAOD *objTool; //!
		std::string prw_file_;

        std::vector<double> PtBinEdges_;
        std::vector<double> EtaBinEdges_;

        std::string uncertaintyName_;

        bool cleverPrescaleTreating_;
        bool controlPlots_;
        int debug_;

        bool JetEffEmulation_;
        double smearedJetPt_;
        std::vector<double> PtBinEdges_scaling_;
        std::vector<double> EtaBinEdges_scaling_;
        std::vector<double> AdditionalSmearing_;
        std::vector<double> LowerTailScaling_;
        std::vector<double> UpperTailScaling_;

        double AdditionalSmearing_variation_;
        double LowerTailScaling_variation_;
        double UpperTailScaling_variation_;
        std::string smearCollection_; // "Gen" or "Reco"

        std::string jetTag_;
        std::string genJetTag_;
        std::string electronTag_;
        std::string muonTag_;
        std::string btagTag_;
        double btagCut_;

        std::string smearingfile_;
        std::string inputhistPtHF_;
        std::string inputhistEtaHF_;
        std::string inputhistPhiHF_;
        std::string inputhistPtLF_;
        std::string inputhistEtaLF_;
        std::string inputhistPhiLF_;

        bool absoluteTailScaling_;
        double A0RMS_;
        double A1RMS_;
        double probExtreme_;

        double rebalancedJetPt_;
        std::string rebalanceMode_; // "MHTall", "MHThigh" or "MET" only for smearCollection = "Reco"
        std::string RebalanceCorrectionFile_;
        std::string genMHTprobFile_;
        std::string METsoftResolutionFile_;
        bool useRebalanceCorrectionFactors_;
        bool useCleverRebalanceCorrectionFactors_;
        bool useGenMHTprob_;
        bool useMETsoftResolution_;

        double JetsHTPt_, JetsHTEta_;
        double JetsMHTPt_, JetsMHTEta_;
        double JetsPt_, JetsEta_;
        double BJetsPt_, BJetsEta_;
        vector<double> JetDeltaMin_;

        double HTSeed;
        int NJetSeed;
        double HTSeedMin_;
        int NJetsSeedMin_;
        bool doSmearing_;

        std::string outputfile_;
        bool storeMHTtree_;
        std::string outputfileMHT_;
        int NJetsStored_;
        int Ntries_;
        int NJetsSave_;
        double HTSave_;
        double METSave_;
        
        UShort_t Njets_stored;

        // variables that don't get filled at submission time should be
        // protected from being send from the submission node to the worker
        // node (done by the //!)

        // this is a standard constructor
        RandS ();

    private:

        typedef TLorentzVector LorentzVector;
        typedef std::vector<std::string>::const_iterator StrIter;
        typedef struct myJet {
            LorentzVector momentum;
            bool btag;
        } myJet;

        TRandom3 *rand_; //!
        SmearFunction *smearFunc_; //!

		bool IsReconstructed(const double&, const double&);
        double JetResolution_Pt2(const double&, const double&);
        double JetResolution_Ptrel(const double&, const double&);
        double JetResolution_Eta(const double&, const double&);
        double JetResolution_Phi(const double&, const double&);
        double JetResolutionHist_Pt_Smear(const double&, const double&, const int&);
        int GetIndex(const double&, const std::vector<double>*);

        double calcHT(std::vector<myJet>&);
        LorentzVector calcMHT(std::vector<myJet>&, const double&, const double&);
        int calcNJets(std::vector<myJet>&);
        int calcNBJets(std::vector<myJet>&);
        bool calcMinDeltaPhi(std::vector<myJet>&, LorentzVector&);
        void calcPredictions(std::vector<myJet>&, LorentzVector&, const int&, const double&);
        void calcLeadingJetPredictions(std::vector<myJet>&, LorentzVector&);
        double GetRebalanceCorrection(double jet_pt, bool btag);

        bool RebalanceJets_KinFitter(std::vector<myJet>&, std::vector<myJet>&, LorentzVector&);
        void SmearingJets(std::vector<myJet>&, LorentzVector&, const double&);

        // these are the functions inherited from Algorithm
        virtual EL::StatusCode setupJob (EL::Job& job);
        virtual EL::StatusCode fileExecute ();
        virtual EL::StatusCode histInitialize ();
        virtual EL::StatusCode changeInput (bool firstFile);
        virtual EL::StatusCode initialize ();
        virtual EL::StatusCode execute ();
        virtual EL::StatusCode postExecute ();
        virtual EL::StatusCode finalize ();
        virtual EL::StatusCode histFinalize ();

        double weight_;

        TTree *PredictionTree; //!
        UShort_t vtxN; //!
        UShort_t Njets_pred; //!
        UShort_t BTags_pred; //!
        Int_t Ntries_pred; //!
        Float_t HT_pred; //!
        Float_t MHT_pred; //!
        Float_t MET_pred; //!
        std::vector<Float_t>  JetPt_p; //!
		std::vector<Float_t> * JetPt_pred = &JetPt_p; //!
        std::vector<Float_t>  JetEta_p; //!
		std::vector<Float_t> * JetEta_pred = &JetEta_p; //!
        std::vector<Float_t>  DeltaPhi_p; //!
		std::vector<Float_t> * DeltaPhi_pred = &DeltaPhi_p; //!
        Float_t weight; //!

        TTree *MHTTree; //!
        Float_t HTreco; //!
        Float_t MHTreco_pt; //!
        Float_t MHTreco_phi; //!
        Float_t MHTrecolow_pt; //!
        Float_t MHTrecolow_phi; //!
        Float_t MHTreb_pt; //!
        Float_t MHTreb_phi; //!
        Float_t MHTreblow_pt; //!
        Float_t MHTreblow_phi; //!
        Float_t HTgen; //!
        Float_t METgen_pt; //!
        Float_t METgen_phi; //!
        Float_t MHTgen_pt; //!
        Float_t MHTgen_phi; //!
        Float_t MHTgenreb_pt; //!
        Float_t MHTgenreb_phi; //!
        Float_t MHTtruereb_pt; //!
        Float_t MHTtruereb_phi; //!
        Float_t MET_pt; //!
        Float_t MET_phi; //!

        std::map <const myJet*, bool> genJet_btag; //!

        TH2F* h_RebCorrection_vsReco, *h_RebCorrection_vsReco_b; //!
        TH1F* h_RebCorrectionFactor, *h_RebCorrectionFactor_b; //!
        TH2F* h_2DRebCorrectionFactor, *h_2DRebCorrectionFactor_b; //!
        vector <TH1D*> h_2DRebCorrectionFactor_py, h_2DRebCorrectionFactor_b_py;  //!
        TH3F* h_MHTtrueProb, *h_MHTtrueProb_input; //!
        TH1F* h_METsoft, *h_METsoft_resPt, *h_METsoft_resPhi; //!
        vector < vector <TH1D*> > h_MHTtrueProb_pz; //!

        // this is needed to distribute the algorithm to the workers
        ClassDef(RandS, 1);
};


template<typename T>
struct GreaterByPt {
    typedef T first_argument_type;
    typedef T second_argument_type;
    bool operator()( const T & t1, const T & t2 ) const {
        return t1.momentum.Pt() > t2.momentum.Pt();
    }
};

#endif

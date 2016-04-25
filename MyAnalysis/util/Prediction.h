#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <TH1.h>
#include <TH2.h>

// Headers needed by this particular selector
#include <vector>

using namespace std;

class Prediction {

    public:
        Prediction (TChain&);
        ~Prediction();

        TH1F* GetSelectionHisto(TString type);
        TH1F* GetPredictionHisto(TString type);
        double GetResultValue(TH1F* histo, double MHTBound1, double MHTBound2);
        double GetResultError(TH1F* histo, double MHTBound1, double MHTBound2);

    private:
        //TTree* QCDPrediction;
        Int_t Ntries;
        UShort_t NJets;
        UShort_t BTags;
        Float_t weight;
        Float_t HT;
        Float_t MHT;
        std::vector<float> *JetPt;
        std::vector<float> *JetEta;
        std::vector<float> *DeltaPhi;

        // store deltaPhi cut
        vector<bool> MinDeltaPhiCut;

        // ------------------------------------------------------- //
        // raw prediction histograms
        // preselection
        TH2F* HT_presel_pred_raw;
        TH2F* MHT_presel_pred_raw;
        TH2F* NJets_presel_pred_raw;
        TH2F* NBJets_presel_pred_raw;
        TH2F* Jet1Pt_presel_pred_raw;
        TH2F* Jet2Pt_presel_pred_raw;
        TH2F* Jet3Pt_presel_pred_raw;
        TH2F* Jet1Eta_presel_pred_raw;
        TH2F* Jet2Eta_presel_pred_raw;
        TH2F* Jet3Eta_presel_pred_raw;
        TH2F* DeltaPhi1_presel_pred_raw;
        TH2F* DeltaPhi2_presel_pred_raw;
        TH2F* DeltaPhi3_presel_pred_raw;

        // preselection + delta phi cut
        TH2F* HT_deltaPhi_pred_raw;
        TH2F* MHT_deltaPhi_pred_raw;
        TH2F* Jet1Pt_deltaPhi_pred_raw;
        TH2F* Jet2Pt_deltaPhi_pred_raw;
        TH2F* Jet3Pt_deltaPhi_pred_raw;
        TH2F* Jet1Eta_deltaPhi_pred_raw;
        TH2F* Jet2Eta_deltaPhi_pred_raw;
        TH2F* Jet3Eta_deltaPhi_pred_raw;

        // NJets
        TH2F* NJets_baseline_withoutMHT_pred_raw;
        TH2F* NJets_baseline_pred_raw;
        TH2F* NJets_baseline_withoutDeltaPhi_withoutMHT_pred_raw;
        TH2F* NJets_baseline_withoutDeltaPhi_pred_raw;

        // NBJets
        TH2F* NBJets_baseline_withoutMHT_pred_raw;
        TH2F* NBJets_baseline_pred_raw;
        TH2F* NBJets_baseline_withoutDeltaPhi_withoutMHT_pred_raw;
        TH2F* NBJets_baseline_withoutDeltaPhi_pred_raw;

        // baseline
        TH2F* HT_baseline_pred_raw;
        TH2F* MHT_baseline_pred_raw;

        // baseline jet bin 1
        TH2F* Jet1Pt_JetBin1_baseline_pred_raw;
        TH2F* Jet2Pt_JetBin1_baseline_pred_raw;
        TH2F* Jet1Eta_JetBin1_baseline_pred_raw;
        TH2F* Jet2Eta_JetBin1_baseline_pred_raw;
        TH2F* DeltaPhi1_JetBin1_baseline_pred_raw;
        TH2F* DeltaPhi2_JetBin1_baseline_pred_raw;

        // baseline jet bin 2
        TH2F* Jet1Pt_JetBin2_baseline_pred_raw;
        TH2F* Jet2Pt_JetBin2_baseline_pred_raw;
        TH2F* Jet3Pt_JetBin2_baseline_pred_raw;
        TH2F* Jet1Eta_JetBin2_baseline_pred_raw;
        TH2F* Jet2Eta_JetBin2_baseline_pred_raw;
        TH2F* Jet3Eta_JetBin2_baseline_pred_raw;
        TH2F* DeltaPhi1_JetBin2_baseline_pred_raw;
        TH2F* DeltaPhi2_JetBin2_baseline_pred_raw;
        TH2F* DeltaPhi3_JetBin2_baseline_pred_raw;

        // baseline jet bin 3
        TH2F* Jet1Pt_JetBin3_baseline_pred_raw;
        TH2F* Jet2Pt_JetBin3_baseline_pred_raw;
        TH2F* Jet3Pt_JetBin3_baseline_pred_raw;
        TH2F* Jet1Eta_JetBin3_baseline_pred_raw;
        TH2F* Jet2Eta_JetBin3_baseline_pred_raw;
        TH2F* Jet3Eta_JetBin3_baseline_pred_raw;
        TH2F* DeltaPhi1_JetBin3_baseline_pred_raw;
        TH2F* DeltaPhi2_JetBin3_baseline_pred_raw;
        TH2F* DeltaPhi3_JetBin3_baseline_pred_raw;

        // baseline jet bin 4
        TH2F* Jet1Pt_JetBin4_baseline_pred_raw;
        TH2F* Jet2Pt_JetBin4_baseline_pred_raw;
        TH2F* Jet3Pt_JetBin4_baseline_pred_raw;
        TH2F* Jet1Eta_JetBin4_baseline_pred_raw;
        TH2F* Jet2Eta_JetBin4_baseline_pred_raw;
        TH2F* Jet3Eta_JetBin4_baseline_pred_raw;
        TH2F* DeltaPhi1_JetBin4_baseline_pred_raw;
        TH2F* DeltaPhi2_JetBin4_baseline_pred_raw;
        TH2F* DeltaPhi3_JetBin4_baseline_pred_raw;

        // baseline without deltaPhi jet bin 1
        TH2F* HT_JetBin1_baseline_withoutDeltaPhi_pred_raw;
        TH2F* MHT_JetBin1_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet1Pt_JetBin1_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet2Pt_JetBin1_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet3Pt_JetBin1_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet1Eta_JetBin1_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet2Eta_JetBin1_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet3Eta_JetBin1_baseline_withoutDeltaPhi_pred_raw;
        TH2F* DeltaPhi1_JetBin1_baseline_withoutDeltaPhi_pred_raw;
        TH2F* DeltaPhi2_JetBin1_baseline_withoutDeltaPhi_pred_raw;
        TH2F* DeltaPhi3_JetBin1_baseline_withoutDeltaPhi_pred_raw;

        // baseline without deltaPhi jet bin 2
        TH2F* HT_JetBin2_baseline_withoutDeltaPhi_pred_raw;
        TH2F* MHT_JetBin2_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet1Pt_JetBin2_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet2Pt_JetBin2_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet3Pt_JetBin2_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet1Eta_JetBin2_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet2Eta_JetBin2_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet3Eta_JetBin2_baseline_withoutDeltaPhi_pred_raw;
        TH2F* DeltaPhi1_JetBin2_baseline_withoutDeltaPhi_pred_raw;
        TH2F* DeltaPhi2_JetBin2_baseline_withoutDeltaPhi_pred_raw;
        TH2F* DeltaPhi3_JetBin2_baseline_withoutDeltaPhi_pred_raw;

        // baseline without deltaPhi jet bin 3
        TH2F* HT_JetBin3_baseline_withoutDeltaPhi_pred_raw;
        TH2F* MHT_JetBin3_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet1Pt_JetBin3_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet2Pt_JetBin3_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet3Pt_JetBin3_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet1Eta_JetBin3_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet2Eta_JetBin3_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet3Eta_JetBin3_baseline_withoutDeltaPhi_pred_raw;
        TH2F* DeltaPhi1_JetBin3_baseline_withoutDeltaPhi_pred_raw;
        TH2F* DeltaPhi2_JetBin3_baseline_withoutDeltaPhi_pred_raw;
        TH2F* DeltaPhi3_JetBin3_baseline_withoutDeltaPhi_pred_raw;

        // baseline without deltaPhi jet bin 4
        TH2F* HT_JetBin4_baseline_withoutDeltaPhi_pred_raw;
        TH2F* MHT_JetBin4_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet1Pt_JetBin4_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet2Pt_JetBin4_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet3Pt_JetBin4_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet1Eta_JetBin4_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet2Eta_JetBin4_baseline_withoutDeltaPhi_pred_raw;
        TH2F* Jet3Eta_JetBin4_baseline_withoutDeltaPhi_pred_raw;
        TH2F* DeltaPhi1_JetBin4_baseline_withoutDeltaPhi_pred_raw;
        TH2F* DeltaPhi2_JetBin4_baseline_withoutDeltaPhi_pred_raw;
        TH2F* DeltaPhi3_JetBin4_baseline_withoutDeltaPhi_pred_raw;

        // HT inclusive 2-3 jets
        TH2F* MHT_JetBin1_HTinclusive_pred_raw;
        // HT inclusive 4-6 jets
        TH2F* MHT_JetBin2_HTinclusive_pred_raw;
        // HT inclusive 7-8 jets
        TH2F* MHT_JetBin3_HTinclusive_pred_raw;
        // HT inclusive >=9 jets
        TH2F* MHT_JetBin4_HTinclusive_pred_raw;

        // ------------------------------------------------------- //
        // prediction histograms
        // preselection
        TH1F* HT_presel_pred;
        TH1F* MHT_presel_pred;
        TH1F* NJets_presel_pred;
        TH1F* NBJets_presel_pred;
        TH1F* Jet1Pt_presel_pred;
        TH1F* Jet2Pt_presel_pred;
        TH1F* Jet3Pt_presel_pred;
        TH1F* Jet1Eta_presel_pred;
        TH1F* Jet2Eta_presel_pred;
        TH1F* Jet3Eta_presel_pred;
        TH1F* DeltaPhi1_presel_pred;
        TH1F* DeltaPhi2_presel_pred;
        TH1F* DeltaPhi3_presel_pred;

        // preselection + delta phi
        TH1F* HT_deltaPhi_pred;
        TH1F* MHT_deltaPhi_pred;
        TH1F* Jet1Pt_deltaPhi_pred;
        TH1F* Jet2Pt_deltaPhi_pred;
        TH1F* Jet3Pt_deltaPhi_pred;
        TH1F* Jet1Eta_deltaPhi_pred;
        TH1F* Jet2Eta_deltaPhi_pred;
        TH1F* Jet3Eta_deltaPhi_pred;

        // NJets
        TH1F* NJets_baseline_withoutMHT_pred;
        TH1F* NJets_baseline_pred;
        TH1F* NJets_baseline_withoutDeltaPhi_withoutMHT_pred;
        TH1F* NJets_baseline_withoutDeltaPhi_pred;

        // NBJets
        TH1F* NBJets_baseline_withoutMHT_pred;
        TH1F* NBJets_baseline_pred;
        TH1F* NBJets_baseline_withoutDeltaPhi_withoutMHT_pred;
        TH1F* NBJets_baseline_withoutDeltaPhi_pred;

        // baseline
        TH1F* HT_baseline_pred;
        TH1F* MHT_baseline_pred;

        // baseline jet bin 1
        TH1F* Jet1Pt_JetBin1_baseline_pred;
        TH1F* Jet2Pt_JetBin1_baseline_pred;
        TH1F* Jet3Pt_JetBin1_baseline_pred;
        TH1F* Jet1Eta_JetBin1_baseline_pred;
        TH1F* Jet2Eta_JetBin1_baseline_pred;
        TH1F* Jet3Eta_JetBin1_baseline_pred;
        TH1F* DeltaPhi1_JetBin1_baseline_pred;
        TH1F* DeltaPhi2_JetBin1_baseline_pred;
        TH1F* DeltaPhi3_JetBin1_baseline_pred;

        // baseline jet bin 2
        TH1F* Jet1Pt_JetBin2_baseline_pred;
        TH1F* Jet2Pt_JetBin2_baseline_pred;
        TH1F* Jet3Pt_JetBin2_baseline_pred;
        TH1F* Jet1Eta_JetBin2_baseline_pred;
        TH1F* Jet2Eta_JetBin2_baseline_pred;
        TH1F* Jet3Eta_JetBin2_baseline_pred;
        TH1F* DeltaPhi1_JetBin2_baseline_pred;
        TH1F* DeltaPhi2_JetBin2_baseline_pred;
        TH1F* DeltaPhi3_JetBin2_baseline_pred;

        // baseline jet bin 3
        TH1F* Jet1Pt_JetBin3_baseline_pred;
        TH1F* Jet2Pt_JetBin3_baseline_pred;
        TH1F* Jet3Pt_JetBin3_baseline_pred;
        TH1F* Jet1Eta_JetBin3_baseline_pred;
        TH1F* Jet2Eta_JetBin3_baseline_pred;
        TH1F* Jet3Eta_JetBin3_baseline_pred;
        TH1F* DeltaPhi1_JetBin3_baseline_pred;
        TH1F* DeltaPhi2_JetBin3_baseline_pred;
        TH1F* DeltaPhi3_JetBin3_baseline_pred;

        // baseline jet bin 4
        TH1F* Jet1Pt_JetBin4_baseline_pred;
        TH1F* Jet2Pt_JetBin4_baseline_pred;
        TH1F* Jet3Pt_JetBin4_baseline_pred;
        TH1F* Jet1Eta_JetBin4_baseline_pred;
        TH1F* Jet2Eta_JetBin4_baseline_pred;
        TH1F* Jet3Eta_JetBin4_baseline_pred;
        TH1F* DeltaPhi1_JetBin4_baseline_pred;
        TH1F* DeltaPhi2_JetBin4_baseline_pred;
        TH1F* DeltaPhi3_JetBin4_baseline_pred;

        // baseline without delta Phi jet bin 1
        TH1F* HT_JetBin1_baseline_withoutDeltaPhi_pred;
        TH1F* MHT_JetBin1_baseline_withoutDeltaPhi_pred;
        TH1F* Jet1Pt_JetBin1_baseline_withoutDeltaPhi_pred;
        TH1F* Jet2Pt_JetBin1_baseline_withoutDeltaPhi_pred;
        TH1F* Jet3Pt_JetBin1_baseline_withoutDeltaPhi_pred;
        TH1F* Jet1Eta_JetBin1_baseline_withoutDeltaPhi_pred;
        TH1F* Jet2Eta_JetBin1_baseline_withoutDeltaPhi_pred;
        TH1F* Jet3Eta_JetBin1_baseline_withoutDeltaPhi_pred;
        TH1F* DeltaPhi1_JetBin1_baseline_withoutDeltaPhi_pred;
        TH1F* DeltaPhi2_JetBin1_baseline_withoutDeltaPhi_pred;
        TH1F* DeltaPhi3_JetBin1_baseline_withoutDeltaPhi_pred;

        // baseline without delta Phi jet bin 2
        TH1F* HT_JetBin2_baseline_withoutDeltaPhi_pred;
        TH1F* MHT_JetBin2_baseline_withoutDeltaPhi_pred;
        TH1F* Jet1Pt_JetBin2_baseline_withoutDeltaPhi_pred;
        TH1F* Jet2Pt_JetBin2_baseline_withoutDeltaPhi_pred;
        TH1F* Jet3Pt_JetBin2_baseline_withoutDeltaPhi_pred;
        TH1F* Jet1Eta_JetBin2_baseline_withoutDeltaPhi_pred;
        TH1F* Jet2Eta_JetBin2_baseline_withoutDeltaPhi_pred;
        TH1F* Jet3Eta_JetBin2_baseline_withoutDeltaPhi_pred;
        TH1F* DeltaPhi1_JetBin2_baseline_withoutDeltaPhi_pred;
        TH1F* DeltaPhi2_JetBin2_baseline_withoutDeltaPhi_pred;
        TH1F* DeltaPhi3_JetBin2_baseline_withoutDeltaPhi_pred;

        // baseline without delta Phi jet bin 3
        TH1F* HT_JetBin3_baseline_withoutDeltaPhi_pred;
        TH1F* MHT_JetBin3_baseline_withoutDeltaPhi_pred;
        TH1F* Jet1Pt_JetBin3_baseline_withoutDeltaPhi_pred;
        TH1F* Jet2Pt_JetBin3_baseline_withoutDeltaPhi_pred;
        TH1F* Jet3Pt_JetBin3_baseline_withoutDeltaPhi_pred;
        TH1F* Jet1Eta_JetBin3_baseline_withoutDeltaPhi_pred;
        TH1F* Jet2Eta_JetBin3_baseline_withoutDeltaPhi_pred;
        TH1F* Jet3Eta_JetBin3_baseline_withoutDeltaPhi_pred;
        TH1F* DeltaPhi1_JetBin3_baseline_withoutDeltaPhi_pred;
        TH1F* DeltaPhi2_JetBin3_baseline_withoutDeltaPhi_pred;
        TH1F* DeltaPhi3_JetBin3_baseline_withoutDeltaPhi_pred;

        // baseline without delta Phi jet bin 4
        TH1F* HT_JetBin4_baseline_withoutDeltaPhi_pred;
        TH1F* MHT_JetBin4_baseline_withoutDeltaPhi_pred;
        TH1F* Jet1Pt_JetBin4_baseline_withoutDeltaPhi_pred;
        TH1F* Jet2Pt_JetBin4_baseline_withoutDeltaPhi_pred;
        TH1F* Jet3Pt_JetBin4_baseline_withoutDeltaPhi_pred;
        TH1F* Jet1Eta_JetBin4_baseline_withoutDeltaPhi_pred;
        TH1F* Jet2Eta_JetBin4_baseline_withoutDeltaPhi_pred;
        TH1F* Jet3Eta_JetBin4_baseline_withoutDeltaPhi_pred;
        TH1F* DeltaPhi1_JetBin4_baseline_withoutDeltaPhi_pred;
        TH1F* DeltaPhi2_JetBin4_baseline_withoutDeltaPhi_pred;
        TH1F* DeltaPhi3_JetBin4_baseline_withoutDeltaPhi_pred;

        // HT inclusive 2-3 jets
        TH1F* MHT_JetBin1_HTinclusive_pred;
        // HT inclusive 4-6 jets
        TH1F* MHT_JetBin2_HTinclusive_pred;
        // HT inclusive 7-8 jets
        TH1F* MHT_JetBin3_HTinclusive_pred;
        // HT inclusive >=9 jets
        TH1F* MHT_JetBin4_HTinclusive_pred;
        // ------------------------------------------------------- //

        // selection histograms
        // preselection
        TH1F* HT_presel_sel;
        TH1F* MHT_presel_sel;
        TH1F* NJets_presel_sel;
        TH1F* NBJets_presel_sel;
        TH1F* Jet1Pt_presel_sel;
        TH1F* Jet2Pt_presel_sel;
        TH1F* Jet3Pt_presel_sel;
        TH1F* Jet1Eta_presel_sel;
        TH1F* Jet2Eta_presel_sel;
        TH1F* Jet3Eta_presel_sel;
        TH1F* DeltaPhi1_presel_sel;
        TH1F* DeltaPhi2_presel_sel;
        TH1F* DeltaPhi3_presel_sel;

        // preselection + delta phi
        TH1F* HT_deltaPhi_sel;
        TH1F* MHT_deltaPhi_sel;
        TH1F* Jet1Pt_deltaPhi_sel;
        TH1F* Jet2Pt_deltaPhi_sel;
        TH1F* Jet3Pt_deltaPhi_sel;
        TH1F* Jet1Eta_deltaPhi_sel;
        TH1F* Jet2Eta_deltaPhi_sel;
        TH1F* Jet3Eta_deltaPhi_sel;

        // NJets
        TH1F* NJets_baseline_withoutMHT_sel;
        TH1F* NJets_baseline_sel;
        TH1F* NJets_baseline_withoutDeltaPhi_withoutMHT_sel;
        TH1F* NJets_baseline_withoutDeltaPhi_sel;

        // NJets
        TH1F* NBJets_baseline_withoutMHT_sel;
        TH1F* NBJets_baseline_sel;
        TH1F* NBJets_baseline_withoutDeltaPhi_withoutMHT_sel;
        TH1F* NBJets_baseline_withoutDeltaPhi_sel;

        // baseline
        TH1F* HT_baseline_sel;
        TH1F* MHT_baseline_sel;

        // baseline jet bin 1
        TH1F* Jet1Pt_JetBin1_baseline_sel;
        TH1F* Jet2Pt_JetBin1_baseline_sel;
        TH1F* Jet3Pt_JetBin1_baseline_sel;
        TH1F* Jet1Eta_JetBin1_baseline_sel;
        TH1F* Jet2Eta_JetBin1_baseline_sel;
        TH1F* Jet3Eta_JetBin1_baseline_sel;
        TH1F* DeltaPhi1_JetBin1_baseline_sel;
        TH1F* DeltaPhi2_JetBin1_baseline_sel;
        TH1F* DeltaPhi3_JetBin1_baseline_sel;

        // baseline jet bin 2
        TH1F* Jet1Pt_JetBin2_baseline_sel;
        TH1F* Jet2Pt_JetBin2_baseline_sel;
        TH1F* Jet3Pt_JetBin2_baseline_sel;
        TH1F* Jet1Eta_JetBin2_baseline_sel;
        TH1F* Jet2Eta_JetBin2_baseline_sel;
        TH1F* Jet3Eta_JetBin2_baseline_sel;
        TH1F* DeltaPhi1_JetBin2_baseline_sel;
        TH1F* DeltaPhi2_JetBin2_baseline_sel;
        TH1F* DeltaPhi3_JetBin2_baseline_sel;

        // baseline jet bin 3
        TH1F* Jet1Pt_JetBin3_baseline_sel;
        TH1F* Jet2Pt_JetBin3_baseline_sel;
        TH1F* Jet3Pt_JetBin3_baseline_sel;
        TH1F* Jet1Eta_JetBin3_baseline_sel;
        TH1F* Jet2Eta_JetBin3_baseline_sel;
        TH1F* Jet3Eta_JetBin3_baseline_sel;
        TH1F* DeltaPhi1_JetBin3_baseline_sel;
        TH1F* DeltaPhi2_JetBin3_baseline_sel;
        TH1F* DeltaPhi3_JetBin3_baseline_sel;

        // baseline jet bin 4
        TH1F* Jet1Pt_JetBin4_baseline_sel;
        TH1F* Jet2Pt_JetBin4_baseline_sel;
        TH1F* Jet3Pt_JetBin4_baseline_sel;
        TH1F* Jet1Eta_JetBin4_baseline_sel;
        TH1F* Jet2Eta_JetBin4_baseline_sel;
        TH1F* Jet3Eta_JetBin4_baseline_sel;
        TH1F* DeltaPhi1_JetBin4_baseline_sel;
        TH1F* DeltaPhi2_JetBin4_baseline_sel;
        TH1F* DeltaPhi3_JetBin4_baseline_sel;

        // baseline without delta Phi jet bin 1
        TH1F* HT_JetBin1_baseline_withoutDeltaPhi_sel;
        TH1F* MHT_JetBin1_baseline_withoutDeltaPhi_sel;
        TH1F* Jet1Pt_JetBin1_baseline_withoutDeltaPhi_sel;
        TH1F* Jet2Pt_JetBin1_baseline_withoutDeltaPhi_sel;
        TH1F* Jet3Pt_JetBin1_baseline_withoutDeltaPhi_sel;
        TH1F* Jet1Eta_JetBin1_baseline_withoutDeltaPhi_sel;
        TH1F* Jet2Eta_JetBin1_baseline_withoutDeltaPhi_sel;
        TH1F* Jet3Eta_JetBin1_baseline_withoutDeltaPhi_sel;
        TH1F* DeltaPhi1_JetBin1_baseline_withoutDeltaPhi_sel;
        TH1F* DeltaPhi2_JetBin1_baseline_withoutDeltaPhi_sel;
        TH1F* DeltaPhi3_JetBin1_baseline_withoutDeltaPhi_sel;

        // baseline without delta Phi jet bin 2
        TH1F* HT_JetBin2_baseline_withoutDeltaPhi_sel;
        TH1F* MHT_JetBin2_baseline_withoutDeltaPhi_sel;
        TH1F* Jet1Pt_JetBin2_baseline_withoutDeltaPhi_sel;
        TH1F* Jet2Pt_JetBin2_baseline_withoutDeltaPhi_sel;
        TH1F* Jet3Pt_JetBin2_baseline_withoutDeltaPhi_sel;
        TH1F* Jet1Eta_JetBin2_baseline_withoutDeltaPhi_sel;
        TH1F* Jet2Eta_JetBin2_baseline_withoutDeltaPhi_sel;
        TH1F* Jet3Eta_JetBin2_baseline_withoutDeltaPhi_sel;
        TH1F* DeltaPhi1_JetBin2_baseline_withoutDeltaPhi_sel;
        TH1F* DeltaPhi2_JetBin2_baseline_withoutDeltaPhi_sel;
        TH1F* DeltaPhi3_JetBin2_baseline_withoutDeltaPhi_sel;

        // baseline without delta Phi jet bin 3
        TH1F* HT_JetBin3_baseline_withoutDeltaPhi_sel;
        TH1F* MHT_JetBin3_baseline_withoutDeltaPhi_sel;
        TH1F* Jet1Pt_JetBin3_baseline_withoutDeltaPhi_sel;
        TH1F* Jet2Pt_JetBin3_baseline_withoutDeltaPhi_sel;
        TH1F* Jet3Pt_JetBin3_baseline_withoutDeltaPhi_sel;
        TH1F* Jet1Eta_JetBin3_baseline_withoutDeltaPhi_sel;
        TH1F* Jet2Eta_JetBin3_baseline_withoutDeltaPhi_sel;
        TH1F* Jet3Eta_JetBin3_baseline_withoutDeltaPhi_sel;
        TH1F* DeltaPhi1_JetBin3_baseline_withoutDeltaPhi_sel;
        TH1F* DeltaPhi2_JetBin3_baseline_withoutDeltaPhi_sel;
        TH1F* DeltaPhi3_JetBin3_baseline_withoutDeltaPhi_sel;

        // baseline without delta Phi jet bin 4
        TH1F* HT_JetBin4_baseline_withoutDeltaPhi_sel;
        TH1F* MHT_JetBin4_baseline_withoutDeltaPhi_sel;
        TH1F* Jet1Pt_JetBin4_baseline_withoutDeltaPhi_sel;
        TH1F* Jet2Pt_JetBin4_baseline_withoutDeltaPhi_sel;
        TH1F* Jet3Pt_JetBin4_baseline_withoutDeltaPhi_sel;
        TH1F* Jet1Eta_JetBin4_baseline_withoutDeltaPhi_sel;
        TH1F* Jet2Eta_JetBin4_baseline_withoutDeltaPhi_sel;
        TH1F* Jet3Eta_JetBin4_baseline_withoutDeltaPhi_sel;
        TH1F* DeltaPhi1_JetBin4_baseline_withoutDeltaPhi_sel;
        TH1F* DeltaPhi2_JetBin4_baseline_withoutDeltaPhi_sel;
        TH1F* DeltaPhi3_JetBin4_baseline_withoutDeltaPhi_sel;

        // HT inclusive 2-3 jets
        TH1F* MHT_JetBin1_HTinclusive_sel;
        // HT inclusive 4-6 jets
        TH1F* MHT_JetBin2_HTinclusive_sel;
        // HT inclusive 7-8 jets
        TH1F* MHT_JetBin3_HTinclusive_sel;
        // HT inclusive >=9 jets
        TH1F* MHT_JetBin4_HTinclusive_sel;
        // ------------------------------------------------------- //

        // dummy histo
        TH1F* dummy;

        bool DeltaPhiCut();
        void DoRebinning(TH2F* prediction_raw, TH1F* selection_raw, int Nbins);
        TH1F* CalcPrediction(TH2F* prediction_raw);

};


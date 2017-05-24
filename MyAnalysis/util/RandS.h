//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 20 14:03:29 2016 by ROOT version 6.06/08
// from TTree EventTree/EventTree
// found on file: sample.root
//////////////////////////////////////////////////////////

#ifndef RandS_h
#define RandS_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>
#include <iostream>
#include <string>
#include <map>

#include "TRandom3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TMatrixD.h"

#include "SmearFunction.h"
#include "MyJet.h"

using namespace std;

class RandS : public TSelector {

    private :
        bool isMC_;
        
        double jvtcut_;
        double lumi_;

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
        std::string smearCollection_; // "Gen for GenSmearing" or "Reco for R+S"

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
        std::string rebalanceMode_; // "MHTall", "MHThigh" or "METsoft" only for smearCollection = "Reco"
        std::string METsoftResolutionFile_;
        std::string triggerTurnOnFile_;
        double maxCleverWeight_;
        bool useMETsoftResolution_;
        bool useTrueMETsoftForRebalance_;
        bool useTriggerTurnOn_;

        double JetsHTPt_, JetsHTEta_;
        double JetsMHTPt_, JetsMHTEta_;
        double JetsPt_, JetsEta_;
        double BJetsPt_, BJetsEta_;
        vector<double> JetDeltaMin_;
        double MjjFirstPt_;
        double MjjSecondPt_;

        double HTSeed;
        int NJetSeed;
        double MHTSeedMax_;
        double MHTSigSeedMax_;
        double HTSeedMin_;
        int NJetsSeedMin_;
        bool doSmearing_;

        std::string outputfile_;
        int NJetsStored_;
        int Ntries_;
        int NJetsSave_;
        double HTSave_;
        double METSave_;
        double MHTSave_;
        double MHTjjSave_;
        double MjjSave_;
        double dPhiSave_;
        double jet3PtSave_;

        UShort_t Njets_stored;

        std::map<UInt_t, UInt_t> ProcessedEvents;
        std::map<UInt_t, UInt_t> AvailableEvents;

        Long64_t NEvents = 0;
        Long64_t NTotEvents = 1;

    public :
        TFile *outputFile = 0;

        TTreeReader fReader; //!the tree reader
        TTree *fChain = 0; //!pointer to the analyzed TTree or TChain

        // Readers to access the data (delete the ones you do not need).
        TTreeReaderValue<Float_t> Weight = {fReader, "Weight"};
        TTreeReaderValue<UInt_t> DatasetID = {fReader, "DatasetID"};
        TTreeReaderValue<Bool_t> PrimaryVtx = {fReader, "PrimaryVtx"};
		TTreeReaderValue<std::vector<Float_t>> JetPt = {fReader, "JetPt"};
        TTreeReaderValue<std::vector<Float_t>> JetEta = {fReader, "JetEta"};
        TTreeReaderValue<std::vector<Float_t>> JetPhi = {fReader, "JetPhi"};
        TTreeReaderValue<std::vector<Float_t>> JetM = {fReader, "JetM"};
        TTreeReaderValue<std::vector<bool>> JetBtag = {fReader, "JetBtag"};
        TTreeReaderValue<std::vector<Float_t>> JetJVT = {fReader, "JetJVT"};
        TTreeReaderValue<std::vector<bool>> JetGood = {fReader, "JetGood"};
        TTreeReaderValue<std::vector<Float_t>> GenJetPt = {fReader, "GenJetPt"};
        TTreeReaderValue<std::vector<Float_t>> GenJetEta = {fReader, "GenJetEta"};
        TTreeReaderValue<std::vector<Float_t>> GenJetPhi = {fReader, "GenJetPhi"};
        TTreeReaderValue<std::vector<Float_t>> GenJetM = {fReader, "GenJetM"};
        TTreeReaderValue<std::vector<bool>> GenJetBtag = {fReader, "GenJetBtag"};
        TTreeReaderValue<std::vector<Float_t>> ElePt = {fReader, "ElePt"};
        TTreeReaderValue<std::vector<Float_t>> EleEta = {fReader, "EleEta"};
        TTreeReaderValue<std::vector<Float_t>> ElePhi = {fReader, "ElePhi"};
        TTreeReaderValue<std::vector<bool>> EleIsSignal = {fReader, "EleIsSignal"};
        TTreeReaderValue<std::vector<Float_t>> PhotonPt = {fReader, "PhotonPt"};
        TTreeReaderValue<std::vector<Float_t>> PhotonEta = {fReader, "PhotonEta"};
        TTreeReaderValue<std::vector<Float_t>> PhotonPhi = {fReader, "PhotonPhi"};
        TTreeReaderValue<std::vector<bool>> PhotonIsSignal = {fReader, "PhotonIsSignal"};
        TTreeReaderValue<std::vector<Float_t>> MuonPt = {fReader, "MuonPt"};
        TTreeReaderValue<std::vector<Float_t>> MuonEta = {fReader, "MuonEta"};
        TTreeReaderValue<std::vector<Float_t>> MuonPhi = {fReader, "MuonPhi"};
        TTreeReaderValue<std::vector<bool>> MuonIsSignal = {fReader, "MuonIsSignal"};
        TTreeReaderValue<std::vector<bool>> MuonIsBad = {fReader, "MuonIsBad"};
        TTreeReaderValue<std::vector<Float_t>> TauPt = {fReader, "TauPt"};
        TTreeReaderValue<std::vector<Float_t>> TauEta = {fReader, "TauEta"};
        TTreeReaderValue<std::vector<Float_t>> TauPhi = {fReader, "TauPhi"};
        TTreeReaderValue<std::vector<bool>> TauIsSignal = {fReader, "TauIsSignal"};
        TTreeReaderValue<Float_t> MET_pt = {fReader, "MET_pt"};
        TTreeReaderValue<Float_t> MET_phi = {fReader, "MET_phi"};
        TTreeReaderValue<Float_t> METmu_pt = {fReader, "METmu_pt"};
        TTreeReaderValue<Float_t> METmu_phi = {fReader, "METmu_phi"};
        TTreeReaderValue<Float_t> GenMET_pt = {fReader, "GenMET_pt"};
        TTreeReaderValue<Float_t> GenMET_phi = {fReader, "GenMET_phi"};
        TTreeReaderValue<Float_t> TrueMHT_pt = {fReader, "TrueMHT_pt"};
        TTreeReaderValue<Float_t> TrueMHT_phi = {fReader, "TrueMHT_phi"};

        RandS(TTree * /*tree*/ =0) { }
        virtual ~RandS() { }
        virtual Int_t   Version() const {
            return 2;
        }
        virtual void Begin(TTree *tree);
        virtual void SlaveBegin(TTree *tree);
        virtual void Init(TTree *tree);
        virtual Bool_t Notify();
        virtual Bool_t Process(Long64_t entry);
        virtual Int_t GetEntry(Long64_t entry, Int_t getall = 0) {
            return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
        }
        virtual void SetOption(const char *option) {
            fOption = option;
        }
        virtual void SetObject(TObject *obj) {
            fObject = obj;
        }
        virtual void SetInputList(TList *input) {
            fInput = input;
        }
        virtual TList *GetOutputList() const {
            return fOutput;
        }
        virtual void SlaveTerminate();
        virtual void Terminate();

        TRandom3 *rand_;
        SmearFunction *smearFunc_;

        bool IsReconstructed(const double&, const double&);
        double JetResolution_Pt2(const double&, const double&);
        double JetResolution_Ptrel(const double&, const double&);
        double JetResolution_Eta(const double&, const double&);
        double JetResolution_Phi(const double&, const double&);
        double JetResolutionHist_Pt_Smear(const double&, const double&, const int&);
        int GetIndex(const double&, const std::vector<double>*);

        double calcHT(std::vector<MyJet>&);
        TLorentzVector calcMHT(std::vector<MyJet>&, const double&, const double&);
        int calcNJets(std::vector<MyJet>&);
        int calcNBJets(std::vector<MyJet>&);
        bool calcMinDeltaPhi(std::vector<MyJet>&, TLorentzVector&);
        void calcPredictions(std::vector<MyJet>&, TLorentzVector&, const int&, const float&);
        void calcLeadingJetPredictions(std::vector<MyJet>&, TLorentzVector&);
        double calcMjj(std::vector<MyJet>&);
        double calcMHTjj(std::vector<MyJet>&);
        double calcDPhijj(std::vector<MyJet>&);
        double calcJet3Pt(std::vector<MyJet>&);
        bool calcMjjSeed(std::vector<MyJet>&, const float&, const float&, const float&);

        bool RebalanceJets_KinFitter(std::vector<MyJet>&, std::vector<MyJet>&, TLorentzVector&);
        void SmearingJets(std::vector<MyJet>&, TLorentzVector&, const float&);

        std::map <const MyJet*, bool> genJet_btag;

		TH2F* h_MHTvsHT_all;
        TH2F* h_MHTvsHT_triggered;
        TH2F* h_MHTvsHT_eff;

        TH2F* h_METsoft_Pt, *h_METsoft_Phi;
        vector <TH1D*> h_METsoft_Pt_px;
        vector <TH1D*> h_METsoft_Phi_px;

        TTree *PredictionTree;
        UShort_t vtxN;
        UShort_t Njets_pred;
        UShort_t BTags_pred;
        Int_t Ntries_pred;
        Float_t HT_pred;
        Float_t MHT_pred;
        Float_t MET_pred;
        std::vector<Float_t>  JetPt_p; //!
        std::vector<Float_t> * JetPt_pred = &JetPt_p;
        std::vector<Float_t>  JetEta_p;
        std::vector<Float_t> * JetEta_pred = &JetEta_p;
        std::vector<Float_t>  JetPhi_p;
        std::vector<Float_t> * JetPhi_pred = &JetPhi_p;
        std::vector<Float_t>  JetM_p;
        std::vector<Float_t> * JetM_pred = &JetM_p;
        std::vector<Float_t>  DeltaPhi_p;
        std::vector<Float_t> * DeltaPhi_pred = &DeltaPhi_p;
        Float_t weight;
        Float_t triggerWeight;

        //ClassDef(RandS,0);

};

#endif

#ifdef RandS_cxx
void RandS::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the reader is initialized.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    std::cout << "Init()" << std::endl;

    fReader.SetTree(tree);
    
    std::cout << "End of Init()" << std::endl;
    
}

Bool_t RandS::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

template<typename T>
struct GreaterByPt {
    typedef T first_argument_type;
    typedef T second_argument_type;
    bool operator()( const T & t1, const T & t2 ) const {
        return t1.Pt() > t2.Pt();
    }
};

#endif // #ifdef RandS_cxx

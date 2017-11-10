//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 20 14:03:29 2016 by ROOT version 6.06/08
// from TTree EventTree/EventTree
// found on file: sample.root
//////////////////////////////////////////////////////////

#ifndef MyABCDStudies_h
#define MyABCDStudies_h

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

#include <TH1F.h>
#include <TH2F.h>

class MyABCDStudies : public TSelector {

    private :
        double m_jvtcut;
        double m_lumi;

        double dPhijjMin_SR;
        double dPhijjMax_SR;
        double dPhijjMin_CR;
        double dPhijjMax_CR;

        double METMin_SR;
        double METMax_SR;
        double METMin_CR;
        double METMax_CR;

		double dEtajjMin_SR;
		double dEtajjMin_CR;

		double dPhiJet1METMin_SR;
		double dPhiJet1METMin_CR;
		
		double dPhiJet2METMin_SR;
		double dPhiJet2METMin_CR;
       
        bool isMC = false;

		//// 1D HISTOGRAMS

        TH1F* h_Jet1_Pt;
        TH1F* h_Jet2_Pt;
        TH1F* h_Jet3_Pt;
        TH1F* h_Jet1_Eta;
        TH1F* h_Jet2_Eta;
        TH1F* h_Jet3_Eta;
        TH1F* h_Jet1_Phi;
        TH1F* h_Jet2_Phi;
        TH1F* h_Jet3_Phi;
        TH1F* h_Jet1_DeltaPhi;
        TH1F* h_Jet2_DeltaPhi;
        TH1F* h_Jet3_DeltaPhi;
        
        TH1F* h_MET;
        TH1F* h_METsig;
        TH1F* h_METsoft;
		TH1F* h_DeltaPhijj;
		TH1F* h_DeltaEtajj;
		TH1F* h_Mjj;

		//// 2D HISTOGRAMS

        TH2F* h_MET_vs_dPhi_Incl;
        TH2F* h_MET_vs_dPhi_Mjj600;
        TH2F* h_MET_vs_dPhi_Mjj1000;
        TH2F* h_MET_vs_dPhi_Mjj1500;
        TH2F* h_MET_vs_dPhi_Mjj2000;

        TH2F* h_MET_vs_dEta_Incl;
        TH2F* h_MET_vs_dEta_Mjj600;
        TH2F* h_MET_vs_dEta_Mjj1000;
        TH2F* h_MET_vs_dEta_Mjj1500;
        TH2F* h_MET_vs_dEta_Mjj2000;

        std::vector<TH1F*> histos_1D;
        std::vector<TH2F*> histos_2D;

        std::map<UInt_t, UInt_t> ProcessedEvents;
        std::map<UInt_t, UInt_t> AvailableEvents;

        Long64_t NEvents = 0;
        Long64_t NTotEvents = 1;
        
        double N_METCR_dPhiCR = 0;
        double N_METCR_dPhiSR = 0;
        double N_METSR_dPhiCR = 0;
        double N_METSR_dPhiSR = 0;

        double N_METCR_dPhiCR_Mjj600 = 0;
        double N_METCR_dPhiSR_Mjj600 = 0;
        double N_METSR_dPhiCR_Mjj600 = 0;
        double N_METSR_dPhiSR_Mjj600 = 0;

        double N_METCR_dPhiCR_Mjj1000 = 0;
        double N_METCR_dPhiSR_Mjj1000 = 0;
        double N_METSR_dPhiCR_Mjj1000 = 0;
        double N_METSR_dPhiSR_Mjj1000 = 0;

        double N_METCR_dPhiCR_Mjj1500 = 0;
        double N_METCR_dPhiSR_Mjj1500 = 0;
        double N_METSR_dPhiCR_Mjj1500 = 0;
        double N_METSR_dPhiSR_Mjj1500 = 0;

        double N_METCR_dPhiCR_Mjj2000 = 0;
        double N_METCR_dPhiSR_Mjj2000 = 0;
        double N_METSR_dPhiCR_Mjj2000 = 0;
        double N_METSR_dPhiSR_Mjj2000 = 0;

        double w2_METCR_dPhiCR = 0;
        double w2_METCR_dPhiSR = 0;
        double w2_METSR_dPhiCR = 0;
        double w2_METSR_dPhiSR = 0;

        double w2_METCR_dPhiCR_Mjj600 = 0;
        double w2_METCR_dPhiSR_Mjj600 = 0;
        double w2_METSR_dPhiCR_Mjj600 = 0;
        double w2_METSR_dPhiSR_Mjj600 = 0;

        double w2_METCR_dPhiCR_Mjj1000 = 0;
        double w2_METCR_dPhiSR_Mjj1000 = 0;
        double w2_METSR_dPhiCR_Mjj1000 = 0;
        double w2_METSR_dPhiSR_Mjj1000 = 0;

        double w2_METCR_dPhiCR_Mjj1500 = 0;
        double w2_METCR_dPhiSR_Mjj1500 = 0;
        double w2_METSR_dPhiCR_Mjj1500 = 0;
        double w2_METSR_dPhiSR_Mjj1500 = 0;

        double w2_METCR_dPhiCR_Mjj2000 = 0;
        double w2_METCR_dPhiSR_Mjj2000 = 0;
        double w2_METSR_dPhiCR_Mjj2000 = 0;
        double w2_METSR_dPhiSR_Mjj2000 = 0;

        double N_METCR_dPhiCR_dEtaCR = 0;
        double N_METCR_dPhiSR_dEtaCR = 0;
        double N_METSR_dPhiCR_dEtaCR = 0;
        double N_METSR_dPhiSR_dEtaCR = 0;

        double N_METCR_dPhiCR_dEtaCR_Mjj600 = 0;
        double N_METCR_dPhiSR_dEtaCR_Mjj600 = 0;
        double N_METSR_dPhiCR_dEtaCR_Mjj600 = 0;
        double N_METSR_dPhiSR_dEtaCR_Mjj600 = 0;

        double N_METCR_dPhiCR_dEtaCR_Mjj1000 = 0;
        double N_METCR_dPhiSR_dEtaCR_Mjj1000 = 0;
        double N_METSR_dPhiCR_dEtaCR_Mjj1000 = 0;
        double N_METSR_dPhiSR_dEtaCR_Mjj1000 = 0;

        double N_METCR_dPhiCR_dEtaCR_Mjj1500 = 0;
        double N_METCR_dPhiSR_dEtaCR_Mjj1500 = 0;
        double N_METSR_dPhiCR_dEtaCR_Mjj1500 = 0;
        double N_METSR_dPhiSR_dEtaCR_Mjj1500 = 0;

        double N_METCR_dPhiCR_dEtaCR_Mjj2000 = 0;
        double N_METCR_dPhiSR_dEtaCR_Mjj2000 = 0;
        double N_METSR_dPhiCR_dEtaCR_Mjj2000 = 0;
        double N_METSR_dPhiSR_dEtaCR_Mjj2000 = 0;

        double w2_METCR_dPhiCR_dEtaCR = 0;
        double w2_METCR_dPhiSR_dEtaCR = 0;
        double w2_METSR_dPhiCR_dEtaCR = 0;
        double w2_METSR_dPhiSR_dEtaCR = 0;

        double w2_METCR_dPhiCR_dEtaCR_Mjj600 = 0;
        double w2_METCR_dPhiSR_dEtaCR_Mjj600 = 0;
        double w2_METSR_dPhiCR_dEtaCR_Mjj600 = 0;
        double w2_METSR_dPhiSR_dEtaCR_Mjj600 = 0;

        double w2_METCR_dPhiCR_dEtaCR_Mjj1000 = 0;
        double w2_METCR_dPhiSR_dEtaCR_Mjj1000 = 0;
        double w2_METSR_dPhiCR_dEtaCR_Mjj1000 = 0;
        double w2_METSR_dPhiSR_dEtaCR_Mjj1000 = 0;

        double w2_METCR_dPhiCR_dEtaCR_Mjj1500 = 0;
        double w2_METCR_dPhiSR_dEtaCR_Mjj1500 = 0;
        double w2_METSR_dPhiCR_dEtaCR_Mjj1500 = 0;
        double w2_METSR_dPhiSR_dEtaCR_Mjj1500 = 0;

        double w2_METCR_dPhiCR_dEtaCR_Mjj2000 = 0;
        double w2_METCR_dPhiSR_dEtaCR_Mjj2000 = 0;
        double w2_METSR_dPhiCR_dEtaCR_Mjj2000 = 0;
        double w2_METSR_dPhiSR_dEtaCR_Mjj2000 = 0;
        
		int N_CR = 0;
        int N_Gen50Lost = 0;
        int N_LowResGen12 = 0;
        int N_HighResReco60 = 0;
        int N_PU60tag = 0;
        int N_PU60notag = 0;
        int N_PU50notag = 0;
        int N_PU40notag = 0;
        int N_Good60tag = 0;
        int N_Good50tag = 0;
        int N_Good40tag = 0;

	public :
        TFile *outputfile = 0;

        TTreeReader fReader; //!the tree reader
        TTree *fChain = 0; //!pointer to the analyzed TTree or TChain

        // Readers to access the data (delete the ones you do not need).
        TTreeReaderValue<Float_t> Weight = {fReader, "Weight"};
        TTreeReaderValue<UInt_t> DatasetID = {fReader, "DatasetID"};
        TTreeReaderValue<UInt_t> EventNo = {fReader, "EventNo"};
        TTreeReaderValue<Bool_t> PrimaryVtx = {fReader, "PrimaryVtx"};
        TTreeReaderValue<Bool_t> xe90triggered = {fReader, "xe90triggered"};
        TTreeReaderValue<Bool_t> xe110triggered = {fReader, "xe110triggered"};
        TTreeReaderValue<std::vector<Float_t>> JetPt = {fReader, "JetPt"};
        TTreeReaderValue<std::vector<Float_t>> JetEta = {fReader, "JetEta"};
        TTreeReaderValue<std::vector<Float_t>> JetPhi = {fReader, "JetPhi"};
        TTreeReaderValue<std::vector<Float_t>> JetM = {fReader, "JetM"};
        TTreeReaderValue<std::vector<bool>> JetBtag = {fReader, "JetBtag"};
        TTreeReaderValue<std::vector<Float_t>> JetJVT = {fReader, "JetJVT"};
        TTreeReaderValue<std::vector<bool>> JetFJVT = {fReader, "JetFJVT"};
        TTreeReaderValue<std::vector<Float_t>> JetSumPtTracks = {fReader, "JetSumPtTracks"};
        TTreeReaderValue<std::vector<Float_t>> JetTrackWidth = {fReader, "JetTrackWidth"};
        TTreeReaderValue<std::vector<UShort_t>> JetNTracks = {fReader, "JetNTracks"};
        TTreeReaderValue<std::vector<bool>> JetGood = {fReader, "JetGood"};
        TTreeReaderValue<std::vector<bool>> JetPassOR = {fReader, "JetPassOR"};
        TTreeReaderValue<std::vector<Float_t>> GenJetPt = {fReader, "GenJetPt"};
        TTreeReaderValue<std::vector<Float_t>> GenJetEta = {fReader, "GenJetEta"};
        TTreeReaderValue<std::vector<Float_t>> GenJetPhi = {fReader, "GenJetPhi"};
        TTreeReaderValue<std::vector<Float_t>> GenJetM = {fReader, "GenJetM"};
        TTreeReaderValue<std::vector<bool>> GenJetBtag = {fReader, "GenJetBtag"};
        TTreeReaderValue<std::vector<Float_t>> GenJetNoNuMuPt = {fReader, "GenJetNoNuMuPt"};
        TTreeReaderValue<std::vector<Float_t>> GenJetNoNuMuEta = {fReader, "GenJetNoNuMuEta"};
        TTreeReaderValue<std::vector<Float_t>> GenJetNoNuMuPhi = {fReader, "GenJetNoNuMuPhi"};
        TTreeReaderValue<std::vector<Float_t>> GenJetNoNuMuM = {fReader, "GenJetNoNuMuM"};
        TTreeReaderValue<std::vector<bool>> GenJetNoNuMuBtag = {fReader, "GenJetNoNuMuBtag"};
        TTreeReaderValue<std::vector<Float_t>> ElePt = {fReader, "ElePt"};
        TTreeReaderValue<std::vector<Float_t>> EleEta = {fReader, "EleEta"};
        TTreeReaderValue<std::vector<Float_t>> ElePhi = {fReader, "ElePhi"};
        TTreeReaderValue<std::vector<bool>> EleIsSignal = {fReader, "EleIsSignal"};
        TTreeReaderValue<std::vector<Int_t>> EleCharge = {fReader, "EleCharge"};
        TTreeReaderValue<std::vector<bool>> ElePassOR = {fReader, "ElePassOR"};
        TTreeReaderValue<std::vector<Float_t>> PhotonPt = {fReader, "PhotonPt"};
        TTreeReaderValue<std::vector<Float_t>> PhotonEta = {fReader, "PhotonEta"};
        TTreeReaderValue<std::vector<Float_t>> PhotonPhi = {fReader, "PhotonPhi"};
        TTreeReaderValue<std::vector<bool>> PhotonIsSignal = {fReader, "PhotonIsSignal"};
        TTreeReaderValue<std::vector<bool>> PhotonPassOR = {fReader, "PhotonPassOR"};
        TTreeReaderValue<std::vector<Float_t>> MuonPt = {fReader, "MuonPt"};
        TTreeReaderValue<std::vector<Float_t>> MuonEta = {fReader, "MuonEta"};
        TTreeReaderValue<std::vector<Float_t>> MuonPhi = {fReader, "MuonPhi"};
        TTreeReaderValue<std::vector<bool>> MuonIsSignal = {fReader, "MuonIsSignal"};
        TTreeReaderValue<std::vector<bool>> MuonIsBad = {fReader, "MuonIsBad"};
        TTreeReaderValue<std::vector<Int_t>> MuonCharge = {fReader, "MuonCharge"};
        TTreeReaderValue<std::vector<bool>> MuonPassOR = {fReader, "MuonPassOR"};
        TTreeReaderValue<std::vector<Float_t>> TauPt = {fReader, "TauPt"};
        TTreeReaderValue<std::vector<Float_t>> TauEta = {fReader, "TauEta"};
        TTreeReaderValue<std::vector<Float_t>> TauPhi = {fReader, "TauPhi"};
        TTreeReaderValue<std::vector<bool>> TauIsSignal = {fReader, "TauIsSignal"};
        TTreeReaderValue<std::vector<Int_t>> TauCharge = {fReader, "TauCharge"};
        TTreeReaderValue<std::vector<bool>> TauPassOR = {fReader, "TauPassOR"};
        TTreeReaderValue<Float_t> MET_pt = {fReader, "MET_pt"};
        TTreeReaderValue<Float_t> MET_phi = {fReader, "MET_phi"};
        TTreeReaderValue<Float_t> METjet_pt = {fReader, "METjet_pt"};
        TTreeReaderValue<Float_t> METjet_phi = {fReader, "METjet_phi"};
        TTreeReaderValue<Float_t> METmu_pt = {fReader, "METmu_pt"};
        TTreeReaderValue<Float_t> METmu_phi = {fReader, "METmu_phi"};
        TTreeReaderValue<Float_t> METele_pt = {fReader, "METele_pt"};
        TTreeReaderValue<Float_t> METele_phi = {fReader, "METele_phi"};
        TTreeReaderValue<Float_t> METgamma_pt = {fReader, "METgamma_pt"};
        TTreeReaderValue<Float_t> METgamma_phi = {fReader, "METgamma_phi"};
        TTreeReaderValue<Float_t> METtrack_pt = {fReader, "METtrack_pt"};
        TTreeReaderValue<Float_t> METtrack_phi = {fReader, "METtrack_phi"};
        TTreeReaderValue<Float_t> GenMET_pt = {fReader, "GenMET_pt"};
        TTreeReaderValue<Float_t> GenMET_phi = {fReader, "GenMET_phi"};
        //TTreeReaderValue<Float_t> TrueMHT_pt = {fReader, "TrueMHT_pt"};
        //TTreeReaderValue<Float_t> TrueMHT_phi = {fReader, "TrueMHT_phi"};

        MyABCDStudies(TTree * /*tree*/ =0) { }
        virtual ~MyABCDStudies() { }
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

        //ClassDef(MyABCDStudies,0);

};

#endif

#ifdef MyABCDStudies_cxx
void MyABCDStudies::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the reader is initialized.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    std::cout << "Init()" << std::endl;

    fReader.SetTree(tree);
}

Bool_t MyABCDStudies::Notify()
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

#endif // #ifdef MyABCDStudies_cxx

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 20 14:03:29 2016 by ROOT version 6.06/08
// from TTree EventTree/EventTree
// found on file: sample.root
//////////////////////////////////////////////////////////

#ifndef MyMETStudies_h
#define MyMETStudies_h

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

class MyMETStudies : public TSelector {

    private :
        double m_MatchingCone;
        double m_jvtcut;
        double m_lumi;
        
        bool isMC = false;
        bool cutFlowStudies = false;

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
        
        TH1F* h_MHT;
        TH1F* h_MHTnoOR;
        TH1F* h_MHTnoJVT;
        TH1F* h_MHTnoJVTnoOR;
        
        TH1F* h_TruthMHT;

        TH1F* h_MET;
        TH1F* h_TruthMET;

        TH1F* h_METreplaced;

		//// 2D HISTOGRAMS

        TH2F* h_MET_vs_MHT;
        TH2F* h_MET_vs_MHTnoOR;
        TH2F* h_MET_vs_MHTnoJVT;
        TH2F* h_MET_vs_MHTnoJVTnoOR;
        
        TH2F* h_MHT_vs_MHTnoOR;
        TH2F* h_MHT_vs_MHTnoJVT;
        TH2F* h_MHT_vs_MHTnoJVTnoOR;

        TH2F* h_JVTPhi_vs_METPhi;
        TH2F* h_JVTPhi_vs_MHTPhi;
        TH2F* h_JVTPhi_vs_MHTnoORPhi;

        TH2F* h_MET_vs_METreplaced;
        TH2F* h_MHTnoOR_vs_MHTnoORreplaced;
        TH2F* h_MHTnoJVTnoOR_vs_MHTnoJVTnoORreplaced;

        TH2F* h_MET_vs_ElePt;
        TH2F* h_MET_vs_MuPt;
        TH2F* h_MET_vs_PhoPt;
        TH2F* h_MET_vs_LepPt;
        TH2F* h_MET_vs_METjet;
        TH2F* h_MET_vs_METele;
        TH2F* h_MET_vs_METmu;
        TH2F* h_MET_vs_METgamma;
        TH2F* h_MET_vs_METtrack;

        std::vector<TH1F*> histos_1D;
        std::vector<TH2F*> histos_2D;

        std::map<UInt_t, UInt_t> ProcessedEvents;
        std::map<UInt_t, UInt_t> AvailableEvents;

        Long64_t NEvents = 0;
        Long64_t NTotEvents = 1;
        
        Long64_t Ntot_VR1 = 0;
        double Ntot_w_VR1 = 0;
        Long64_t N2jets_VR1 = 0;
        double N2jets_w_VR1 = 0;
        Long64_t NjetPt_VR1 = 0;
        double NjetPt_w_VR1 = 0;
        Long64_t NTJV_VR1 = 0;
        double NTJV_w_VR1 = 0;
        Long64_t NdPhi_VR1 = 0;
        double NdPhi_w_VR1 = 0;
        Long64_t NdEta_VR1 = 0;
        double NdEta_w_VR1 = 0;
        Long64_t Nhemi_VR1 = 0;
        double Nhemi_w_VR1 = 0;
        Long64_t NpTjj_VR1 = 0;
        double NpTjj_w_VR1 = 0;
        Long64_t Nmjj_VR1 = 0;
        double Nmjj_w_VR1 = 0;
        Long64_t Nlv_VR1 = 0;
        double Nlv_w_VR1 = 0;

        Long64_t Ntot_VR1_trig = 0;
        double Ntot_w_VR1_trig = 0;
        Long64_t N2jets_VR1_trig = 0;
        double N2jets_w_VR1_trig = 0;
        Long64_t NjetPt_VR1_trig = 0;
        double NjetPt_w_VR1_trig = 0;
        Long64_t NTJV_VR1_trig = 0;
        double NTJV_w_VR1_trig = 0;
        Long64_t NdPhi_VR1_trig = 0;
        double NdPhi_w_VR1_trig = 0;
        Long64_t NdEta_VR1_trig = 0;
        double NdEta_w_VR1_trig = 0;
        Long64_t Nhemi_VR1_trig = 0;
        double Nhemi_w_VR1_trig = 0;
        Long64_t NpTjj_VR1_trig = 0;
        double NpTjj_w_VR1_trig = 0;
        Long64_t Nmjj_VR1_trig = 0;
        double Nmjj_w_VR1_trig = 0;
        Long64_t Nlv_VR1_trig = 0;
        double Nlv_w_VR1_trig = 0;

        Long64_t Ntot_VR2 = 0;
        double Ntot_w_VR2 = 0;
        Long64_t N2jets_VR2 = 0;
        double N2jets_w_VR2 = 0;
        Long64_t NjetPt_VR2 = 0;
        double NjetPt_w_VR2 = 0;
        Long64_t NTJV_VR2 = 0;
        double NTJV_w_VR2 = 0;
        Long64_t NdPhi_VR2 = 0;
        double NdPhi_w_VR2 = 0;
        Long64_t NdEta_VR2 = 0;
        double NdEta_w_VR2 = 0;
        Long64_t Nhemi_VR2 = 0;
        double Nhemi_w_VR2 = 0;
        Long64_t NpTjj_VR2 = 0;
        double NpTjj_w_VR2 = 0;
        Long64_t Nmjj_VR2 = 0;
        double Nmjj_w_VR2 = 0;
        Long64_t Nlv_VR2 = 0;
        double Nlv_w_VR2 = 0;

        Long64_t Ntot_VR2_trig = 0;
        double Ntot_w_VR2_trig = 0;
        Long64_t N2jets_VR2_trig = 0;
        double N2jets_w_VR2_trig = 0;
        Long64_t NjetPt_VR2_trig = 0;
        double NjetPt_w_VR2_trig = 0;
        Long64_t NTJV_VR2_trig = 0;
        double NTJV_w_VR2_trig = 0;
        Long64_t NdPhi_VR2_trig = 0;
        double NdPhi_w_VR2_trig = 0;
        Long64_t NdEta_VR2_trig = 0;
        double NdEta_w_VR2_trig = 0;
        Long64_t Nhemi_VR2_trig = 0;
        double Nhemi_w_VR2_trig = 0;
        Long64_t NpTjj_VR2_trig = 0;
        double NpTjj_w_VR2_trig = 0;
        Long64_t Nmjj_VR2_trig = 0;
        double Nmjj_w_VR2_trig = 0;
        Long64_t Nlv_VR2_trig = 0;
        double Nlv_w_VR2_trig = 0;

        Long64_t Ntot_VR3 = 0;
        double Ntot_w_VR3 = 0;
        Long64_t Ntrig_VR3 = 0;
        double Ntrig_w_VR3 = 0;
        Long64_t N2jets_VR3 = 0;
        double N2jets_w_VR3 = 0;
        Long64_t NjetPt_VR3 = 0;
        double NjetPt_w_VR3 = 0;
        Long64_t NTJV_VR3 = 0;
        double NTJV_w_VR3 = 0;
        Long64_t NdPhi_VR3 = 0;
        double NdPhi_w_VR3 = 0;
        Long64_t NdEta_VR3 = 0;
        double NdEta_w_VR3 = 0;
        Long64_t Nhemi_VR3 = 0;
        double Nhemi_w_VR3 = 0;
        Long64_t NpTjj_VR3 = 0;
        double NpTjj_w_VR3 = 0;
        Long64_t Nmjj_VR3 = 0;
        double Nmjj_w_VR3 = 0;
        Long64_t Nlv_VR3 = 0;
        double Nlv_w_VR3 = 0;
        
        Long64_t Ntot_VR3_trig = 0;
        double Ntot_w_VR3_trig = 0;
        Long64_t Ntrig_VR3_trig = 0;
        double Ntrig_w_VR3_trig = 0;
        Long64_t N2jets_VR3_trig = 0;
        double N2jets_w_VR3_trig = 0;
        Long64_t NjetPt_VR3_trig = 0;
        double NjetPt_w_VR3_trig = 0;
        Long64_t NTJV_VR3_trig = 0;
        double NTJV_w_VR3_trig = 0;
        Long64_t NdPhi_VR3_trig = 0;
        double NdPhi_w_VR3_trig = 0;
        Long64_t NdEta_VR3_trig = 0;
        double NdEta_w_VR3_trig = 0;
        Long64_t Nhemi_VR3_trig = 0;
        double Nhemi_w_VR3_trig = 0;
        Long64_t NpTjj_VR3_trig = 0;
        double NpTjj_w_VR3_trig = 0;
        Long64_t Nmjj_VR3_trig = 0;
        double Nmjj_w_VR3_trig = 0;
        Long64_t Nlv_VR3_trig = 0;
        double Nlv_w_VR3_trig = 0;

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

        MyMETStudies(TTree * /*tree*/ =0) { }
        virtual ~MyMETStudies() { }
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

        //ClassDef(MyMETStudies,0);

};

#endif

#ifdef MyMETStudies_cxx
void MyMETStudies::Init(TTree *tree)
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

Bool_t MyMETStudies::Notify()
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

#endif // #ifdef MyMETStudies_cxx

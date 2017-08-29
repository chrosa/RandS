//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 20 14:03:29 2016 by ROOT version 6.06/08
// from TTree EventTree/EventTree
// found on file: sample.root
//////////////////////////////////////////////////////////

#ifndef MyResolution_h
#define MyResolution_h

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

class MyResolution : public TSelector {

    private :
        double m_VetoCone;
        double m_MatchingCone;
        double m_RelGenActivityVeto;
        double m_RelRecoActivityVeto;
        double m_jvtcut;
        double m_lumi;

        std::vector<double> PtBinEdges;
        std::vector<double> EtaBinEdges;

        // Resize histo vectors
        void ResizeHistoVector(std::vector<std::vector<TH1F*> > &histoVector);

        std::vector< std::vector<TH1F*> > PtResolution_tot;
        std::vector< std::vector<TH1F*> > EtaResolution_tot;
        std::vector< std::vector<TH1F*> > PhiResolution_tot;

        std::vector< std::vector<TH1F*> > PtResolution_LF;
        std::vector< std::vector<TH1F*> > EtaResolution_LF;
        std::vector< std::vector<TH1F*> > PhiResolution_LF;

        std::vector< std::vector<TH1F*> > PtResolution_HF;
        std::vector< std::vector<TH1F*> > EtaResolution_HF;
        std::vector< std::vector<TH1F*> > PhiResolution_HF;

        std::vector< std::vector<TH1F*> > PtResolution_b;
        std::vector< std::vector<TH1F*> > EtaResolution_b;
        std::vector< std::vector<TH1F*> > PhiResolution_b;

        std::vector< std::vector<TH1F*> > PtResolution_nob;
        std::vector< std::vector<TH1F*> > EtaResolution_nob;
        std::vector< std::vector<TH1F*> > PhiResolution_nob;

        std::vector<TH1F*> NReco_tot;
        std::vector<TH1F*> NReco_b;
        std::vector<TH1F*> NReco_nob;

        std::vector<TH1F*> NGen_tot;
        std::vector<TH1F*> NGen_b;
        std::vector<TH1F*> NGen_nob;

		std::map<UInt_t, UInt_t> ProcessedEvents;
		std::map<UInt_t, UInt_t> AvailableEvents;

        std::string GetHistName(unsigned int i_pt, unsigned int i_eta, std::string s1, std::string s2);
        std::string GetHistNameEta(unsigned int i_pt, std::string s1);
        unsigned int GetPtBin(double pt);
        unsigned int GetEtaBin(double eta);
        
        Long64_t NEvents = 0;
        Long64_t NTotEvents = 1;

    public :
        TFile *outputfile = 0;

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
        TTreeReaderValue<std::vector<Float_t>> GenJetPt = {fReader, "GenJetNoNuMuPt"};
        TTreeReaderValue<std::vector<Float_t>> GenJetEta = {fReader, "GenJetNoNuMuEta"};
        TTreeReaderValue<std::vector<Float_t>> GenJetPhi = {fReader, "GenJetNoNuMuPhi"};
        TTreeReaderValue<std::vector<Float_t>> GenJetM = {fReader, "GenJetNoNuMuM"};
        TTreeReaderValue<std::vector<bool>> GenJetBtag = {fReader, "GenJetNoNuMuBtag"};
        //TTreeReaderValue<std::vector<Float_t>> GenJetPt = {fReader, "GenJetPt"};
        //TTreeReaderValue<std::vector<Float_t>> GenJetEta = {fReader, "GenJetEta"};
        //TTreeReaderValue<std::vector<Float_t>> GenJetPhi = {fReader, "GenJetPhi"};
        //TTreeReaderValue<std::vector<Float_t>> GenJetM = {fReader, "GenJetM"};
        //TTreeReaderValue<std::vector<bool>> GenJetBtag = {fReader, "GenJetBtag"};
        TTreeReaderValue<std::vector<Float_t>> ElePt = {fReader, "ElePt"};
        TTreeReaderValue<std::vector<Float_t>> EleEta = {fReader, "EleEta"};
        TTreeReaderValue<std::vector<Float_t>> ElePhi = {fReader, "ElePhi"};
        TTreeReaderValue<std::vector<bool>> EleIsSignal = {fReader, "EleIsSignal"};
        TTreeReaderValue<std::vector<Float_t>> PhotonPt = {fReader, "PhotonPt"};
        TTreeReaderValue<std::vector<Float_t>> PhotonEta = {fReader, "PhotonEta"};
        TTreeReaderValue<std::vector<Float_t>> PhotonPhi = {fReader, "PhotonPhi"};
        TTreeReaderValue<std::vector<Float_t>> MuonPt = {fReader, "MuonPt"};
        TTreeReaderValue<std::vector<Float_t>> MuonEta = {fReader, "MuonEta"};
        TTreeReaderValue<std::vector<Float_t>> MuonPhi = {fReader, "MuonPhi"};
        TTreeReaderValue<std::vector<bool>> MuonIsSignal = {fReader, "MuonIsSignal"};
        TTreeReaderValue<std::vector<bool>> MuonIsBad = {fReader, "MuonIsBad"};
        TTreeReaderValue<Float_t> MET_pt = {fReader, "MET_pt"};
        TTreeReaderValue<Float_t> MET_phi = {fReader, "MET_phi"};
        TTreeReaderValue<Float_t> GenMET_pt = {fReader, "GenMET_pt"};
        TTreeReaderValue<Float_t> GenMET_phi = {fReader, "GenMET_phi"};
        //TTreeReaderValue<Float_t> TrueMHT_pt = {fReader, "TrueMHT_pt"};
        //TTreeReaderValue<Float_t> TrueMHT_phi = {fReader, "TrueMHT_phi"};

        MyResolution(TTree * /*tree*/ =0) { }
        virtual ~MyResolution() { }
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

        //ClassDef(MyResolution,0);

};

#endif

#ifdef MyResolution_cxx
void MyResolution::Init(TTree *tree)
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

Bool_t MyResolution::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}


#endif // #ifdef MyResolution_cxx

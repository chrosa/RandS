//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 11 10:47:53 2016 by ROOT version 6.06/03
// from TTree MHTTree/MHTTree
// found on file: sample.root
//////////////////////////////////////////////////////////

#ifndef METvsMHTSelector_h
#define METvsMHTSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <TCanvas.h>
#include <TPad.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLorentzVector.h>

#include <vector>
#include <map>

#include "MyJet.h"


// Headers needed by this particular selector


class METvsMHTSelector : public TSelector {
    public :
        TTreeReader     fReader;  //!the tree reader
        TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

        // Readers to access the data (delete the ones you do not need).
        TTreeReaderValue<Float_t> Weight = {fReader, "Weight"};
        TTreeReaderValue<UInt_t> DatasetID = {fReader, "DatasetID"};
        TTreeReaderValue<Bool_t> PrimaryVtx = {fReader, "PrimaryVtx"};
        TTreeReaderValue<std::vector<Float_t>> JetPt = {fReader, "JetNoMuPt"};
        TTreeReaderValue<std::vector<Float_t>> JetEta = {fReader, "JetNoMuEta"};
        TTreeReaderValue<std::vector<Float_t>> JetPhi = {fReader, "JetNoMuPhi"};
        TTreeReaderValue<std::vector<Float_t>> JetM = {fReader, "JetNoMuM"};
        TTreeReaderValue<std::vector<bool>> JetBtag = {fReader, "JetNoMuBtag"};
        TTreeReaderValue<std::vector<Float_t>> JetJVT = {fReader, "JetNoMuJVT"};
        TTreeReaderValue<std::vector<bool>> JetGood = {fReader, "JetNoMuGood"};
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
        TTreeReaderValue<std::vector<Float_t>> PhotonPt = {fReader, "PhotonPt"};
        TTreeReaderValue<std::vector<Float_t>> PhotonEta = {fReader, "PhotonEta"};
        TTreeReaderValue<std::vector<Float_t>> PhotonPhi = {fReader, "PhotonPhi"};
        TTreeReaderValue<std::vector<Float_t>> MuonPt = {fReader, "MuonPt"};
        TTreeReaderValue<std::vector<Float_t>> MuonEta = {fReader, "MuonEta"};
        TTreeReaderValue<std::vector<Float_t>> MuonPhi = {fReader, "MuonPhi"};
        TTreeReaderValue<std::vector<bool>> MuonIsBad = {fReader, "MuonIsBad"};
        TTreeReaderValue<std::vector<bool>> MuonIsSignal = {fReader, "MuonIsSignal"};
        TTreeReaderValue<Float_t> MET_pt = {fReader, "MET_pt"};
        TTreeReaderValue<Float_t> MET_phi = {fReader, "MET_phi"};
        TTreeReaderValue<Float_t> METmu_pt = {fReader, "METmu_pt"};
        TTreeReaderValue<Float_t> METmu_phi = {fReader, "METmu_phi"};
        TTreeReaderValue<Float_t> METtrack_pt = {fReader, "METtrack_pt"};
        TTreeReaderValue<Float_t> METtrack_phi = {fReader, "METtrack_phi"};
        TTreeReaderValue<Float_t> GenMET_pt = {fReader, "GenMET_pt"};
        TTreeReaderValue<Float_t> GenMET_phi = {fReader, "GenMET_phi"};
        TTreeReaderValue<Float_t> TrueMHT_pt = {fReader, "TrueMHT_pt"};
        TTreeReaderValue<Float_t> TrueMHT_phi = {fReader, "TrueMHT_phi"};

        METvsMHTSelector(TTree * /*tree*/ =0) { }
        virtual ~METvsMHTSelector() { }
        virtual Int_t   Version() const {
            return 2;
        }
        virtual void    Begin(TTree *tree);
        virtual void    SlaveBegin(TTree *tree);
        virtual void    Init(TTree *tree);
        virtual Bool_t  Notify();
        virtual Bool_t  Process(Long64_t entry);
        virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) {
            return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
        }
        virtual void    SetOption(const char *option) {
            fOption = option;
        }
        virtual void    SetObject(TObject *obj) {
            fObject = obj;
        }
        virtual void    SetInputList(TList *input) {
            fInput = input;
        }
        virtual TList  *GetOutputList() const {
            return fOutput;
        }
        virtual void    SlaveTerminate();
        virtual void    Terminate(); 

        typedef struct MyJet {
            TLorentzVector momentum;
            bool btag;
            bool good;
            double jvt;
        } MyJet;

        ClassDef(METvsMHTSelector,0);

    private:

        TH2F *h_MHTtruerebPhiRes_vs_MHTrebMinusMET;
        TH2F *h_MHTtruerebPtRes_vs_MHTrebMinusMET;
        TH2F *h_MHTtruerebPt_vs_MHTrebMinusMET;
        TH2F *h_MHT_vs_MET;

        int NEvents = 0;
		std::map<UInt_t, UInt_t> AvailableEvents;
		double rebalancedJetPt_ = 20.;
		double jvtcut_ = 0.59;
		double lumi_ = 30000.;

};

#endif

#ifdef METvsMHTSelector_cxx
void METvsMHTSelector::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the reader is initialized.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    fReader.SetTree(tree);
}

Bool_t METvsMHTSelector::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}


#endif // #ifdef METvsMHTSelector_cxx

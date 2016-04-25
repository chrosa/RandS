#ifndef Expectation_h
#define Expectation_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include "RA2bBin.h"

// Headers needed by this particular selector
#include <string>
#include <vector>

// Header file for the classes stored in the TTree if any.

class Expectation : public TSelector {

    public :

        TTreeReader fReader;  //!the tree reader
        TTree *fChain = 0;   //!pointer to the analyzed TTree or TChain

        // Readers to access the data (delete the ones you do not need).
        TTreeReaderValue<Int_t> Ntries = {fReader, "Ntries"};
        TTreeReaderValue<UShort_t> NJets = {fReader, "NJets"};
        TTreeReaderValue<UShort_t> BTags = {fReader, "BTags"};
        TTreeReaderValue<Float_t> Weight = {fReader, "Weight"};
        TTreeReaderValue<Float_t> HT = {fReader, "HT"};
        TTreeReaderValue<Float_t> MHT = {fReader, "MHT"};
        TTreeReaderArray<float> JetPt = {fReader, "JetPt"};
        TTreeReaderArray<float> JetEta = {fReader, "JetEta"};
        TTreeReaderArray<float> DeltaPhi = {fReader, "DeltaPhi"};

        Expectation(TTree * /*tree*/ =0) : fChain(0) { }
        virtual ~Expectation() { }
        virtual Int_t Version() const {
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
        virtual bool DeltaPhiCut();

        std::vector<RA2bBin*> SB;
        TH1F* yields;
        int TotEvents;

        ClassDef(Expectation,0);
};

#endif

#ifdef Expectation_cxx
void Expectation::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    // Set branch addresses and branch pointers
    if (!tree) return;
    fReader.SetTree(tree);

}

Bool_t Expectation::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

#endif // #ifdef Expectation_cxx

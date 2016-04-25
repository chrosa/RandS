#define Expectation_cxx
// The class definition in Expectation.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("Expectation.C")
// root> T->Process("Expectation.C","some options")
// root> T->Process("Expectation.C+")
//

#include "Expectation.h"
#include "RA2bBin.h"
#include <TH1F.h>
#include <TStyle.h>
#include <vector>


void Expectation::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
   
   TotEvents = 0;
   
   TString option = GetOption();
   
   //Define search bins
   //￼Njets bins: 4-6, 7-8, 9+
   //Nb bins: 0, 1, 2, 3+
   //MHT,HT = 200,500,500,800
   //MHT,HT = 200,500,800,1200
   //MHT,HT = 200,500,1200+
   //MHT,HT = 500,750,500,1200
   //MHT,HT = 500,750,1200+
   //MHT,HT = 750+,800+
   
   RA2bBin* SB1 = new RA2bBin("Njet456-Nb0-MHT1-HT1", 4, 6, 0, 0, 200., 500., 500., 800., true);
   RA2bBin* SB2 = new RA2bBin("Njet456-Nb0-MHT1-HT2", 4, 6, 0, 0, 200., 500., 800., 1200., true);
   RA2bBin* SB3 = new RA2bBin("Njet456-Nb0-MHT1-HT3", 4, 6, 0, 0, 200., 500., 1200., 9999., true);
   RA2bBin* SB4 = new RA2bBin("Njet456-Nb0-MHT2-HT1", 4, 6, 0, 0, 500., 750., 500., 1200., true);
   RA2bBin* SB5 = new RA2bBin("Njet456-Nb0-MHT2-HT23", 4, 6, 0, 0, 500., 750., 1200., 9999., true);
   RA2bBin* SB6 = new RA2bBin("Njet456-Nb0-MHT2-HT123", 4, 6, 0, 0, 750., 9999., 800., 9999., true);
   
   RA2bBin* SB7 = new RA2bBin("Njet456-Nb1-MHT1-HT1", 4, 6, 1, 1, 200., 500., 500., 800., true);
   RA2bBin* SB8 = new RA2bBin("Njet456-Nb1-MHT1-HT2", 4, 6, 1, 1, 200., 500., 800., 1200., true);
   RA2bBin* SB9 = new RA2bBin("Njet456-Nb1-MHT1-HT3", 4, 6, 1, 1, 200., 500., 1200., 9999., true);
   RA2bBin* SB10 = new RA2bBin("Njet456-Nb1-MHT2-HT1", 4, 6, 1, 1, 500., 750., 500., 1200., true);
   RA2bBin* SB11 = new RA2bBin("Njet456-Nb1-MHT2-HT23", 4, 6, 1, 1, 500., 750., 1200., 9999., true);
   RA2bBin* SB12 = new RA2bBin("Njet456-Nb1-MHT2-HT123", 4, 6, 1, 1, 750., 9999., 800., 9999., true);
   
   RA2bBin* SB13 = new RA2bBin("Njet456-Nb2-MHT1-HT1", 4, 6, 2, 2, 200., 500., 500., 800., true);
   RA2bBin* SB14 = new RA2bBin("Njet456-Nb2-MHT1-HT2", 4, 6, 2, 2, 200., 500., 800., 1200., true);
   RA2bBin* SB15 = new RA2bBin("Njet456-Nb2-MHT1-HT3", 4, 6, 2, 2, 200., 500., 1200., 9999., true);
   RA2bBin* SB16 = new RA2bBin("Njet456-Nb2-MHT2-HT1", 4, 6, 2, 2, 500., 750., 500., 1200., true);
   RA2bBin* SB17 = new RA2bBin("Njet456-Nb2-MHT2-HT23", 4, 6, 2, 2, 500., 750., 1200., 9999., true);
   RA2bBin* SB18 = new RA2bBin("Njet456-Nb2-MHT2-HT123", 4, 6, 2, 2, 750., 9999., 800., 9999., true);
   
   RA2bBin* SB19 = new RA2bBin("Njet456-Nb3+-MHT1-HT1", 4, 6, 3, 9, 200., 500., 500., 800., true);
   RA2bBin* SB20 = new RA2bBin("Njet456-Nb3+-MHT1-HT2", 4, 6, 3, 9, 200., 500., 800., 1200., true);
   RA2bBin* SB21 = new RA2bBin("Njet456-Nb3+-MHT1-HT3", 4, 6, 3, 9, 200., 500., 1200., 9999., true);
   RA2bBin* SB22 = new RA2bBin("Njet456-Nb3+-MHT2-HT1", 4, 6, 3, 9, 500., 750., 500., 1200., true);
   RA2bBin* SB23 = new RA2bBin("Njet456-Nb3+-MHT2-HT23", 4, 6, 3, 9, 500., 750., 1200., 9999., true);
   RA2bBin* SB24 = new RA2bBin("Njet456-Nb3+-MHT2-HT123", 4, 6, 3, 9, 750., 9999., 800., 9999., true);
   
   RA2bBin* SB25 = new RA2bBin("Njet78-Nb0-MHT1-HT1", 7, 8, 0, 0, 200., 500., 500., 800., true);
   RA2bBin* SB26 = new RA2bBin("Njet78-Nb0-MHT1-HT2", 7, 8, 0, 0, 200., 500., 800., 1200., true);
   RA2bBin* SB27 = new RA2bBin("Njet78-Nb0-MHT1-HT3", 7, 8, 0, 0, 200., 500., 1200., 9999., true);
   RA2bBin* SB28 = new RA2bBin("Njet78-Nb0-MHT2-HT1", 7, 8, 0, 0, 500., 750., 500., 1200., true);
   RA2bBin* SB29 = new RA2bBin("Njet78-Nb0-MHT2-HT23", 7, 8, 0, 0, 500., 750., 1200., 9999., true);
   RA2bBin* SB30 = new RA2bBin("Njet78-Nb0-MHT2-HT123", 7, 8, 0, 0, 750., 9999., 800., 9999., true);
   
   RA2bBin* SB31 = new RA2bBin("Njet78-Nb1-MHT1-HT1", 7, 8, 1, 1, 200., 500., 500., 800., true);
   RA2bBin* SB32 = new RA2bBin("Njet78-Nb1-MHT1-HT2", 7, 8, 1, 1, 200., 500., 800., 1200., true);
   RA2bBin* SB33 = new RA2bBin("Njet78-Nb1-MHT1-HT3", 7, 8, 1, 1, 200., 500., 1200., 9999., true);
   RA2bBin* SB34 = new RA2bBin("Njet78-Nb1-MHT2-HT1", 7, 8, 1, 1, 500., 750., 500., 1200., true);
   RA2bBin* SB35 = new RA2bBin("Njet78-Nb1-MHT2-HT23", 7, 8, 1, 1, 500., 750., 1200., 9999., true);
   RA2bBin* SB36 = new RA2bBin("Njet78-Nb1-MHT2-HT123", 7, 8, 1, 1, 750., 9999., 800., 9999., true);
   
   RA2bBin* SB37 = new RA2bBin("Njet78-Nb2-MHT1-HT1", 7, 8, 2, 2, 200., 500., 500., 800., true);
   RA2bBin* SB38 = new RA2bBin("Njet78-Nb2-MHT1-HT2", 7, 8, 2, 2, 200., 500., 800., 1200., true);
   RA2bBin* SB39 = new RA2bBin("Njet78-Nb2-MHT1-HT3", 7, 8, 2, 2, 200., 500., 1200., 9999., true);
   RA2bBin* SB40 = new RA2bBin("Njet78-Nb2-MHT2-HT1", 7, 8, 2, 2, 500., 750., 500., 1200., true);
   RA2bBin* SB41 = new RA2bBin("Njet78-Nb2-MHT2-HT23", 7, 8, 2, 2, 500., 750., 1200., 9999., true);
   RA2bBin* SB42 = new RA2bBin("Njet78-Nb2-MHT2-HT123", 7, 8, 2, 2, 750., 9999., 800., 9999., true);
   
   RA2bBin* SB43 = new RA2bBin("Njet78-Nb3+-MHT1-HT1", 7, 8, 3, 9, 200., 500., 500., 800., true);
   RA2bBin* SB44 = new RA2bBin("Njet78-Nb3+-MHT1-HT2", 7, 8, 3, 9, 200., 500., 800., 1200., true);
   RA2bBin* SB45 = new RA2bBin("Njet78-Nb3+-MHT1-HT3", 7, 8, 3, 9, 200., 500., 1200., 9999., true);
   RA2bBin* SB46 = new RA2bBin("Njet78-Nb3+-MHT2-HT1", 7, 8, 3, 9, 500., 750., 500., 1200., true);
   RA2bBin* SB47 = new RA2bBin("Njet78-Nb3+-MHT2-HT23", 7, 8, 3, 9, 500., 750., 1200., 9999., true);
   RA2bBin* SB48 = new RA2bBin("Njet78-Nb3+-MHT2-HT123", 7, 8, 3, 9, 750., 9999., 800., 9999., true);
   
   RA2bBin* SB49 = new RA2bBin("Njet9+-Nb0-MHT1-HT1", 9, 99, 0, 0, 200., 500., 500., 800., true);
   RA2bBin* SB50 = new RA2bBin("Njet9+-Nb0-MHT1-HT2", 9, 99, 0, 0, 200., 500., 800., 1200., true);
   RA2bBin* SB51 = new RA2bBin("Njet9+-Nb0-MHT1-HT3", 9, 99, 0, 0, 200., 500., 1200., 9999., true);
   RA2bBin* SB52 = new RA2bBin("Njet9+-Nb0-MHT2-HT1", 9, 99, 0, 0, 500., 750., 500., 1200., true);
   RA2bBin* SB53 = new RA2bBin("Njet9+-Nb0-MHT2-HT23", 9, 99, 0, 0, 500., 750., 1200., 9999., true);
   RA2bBin* SB54 = new RA2bBin("Njet9+-Nb0-MHT2-HT123", 9, 99, 0, 0, 750., 9999., 800., 9999., true);
   
   RA2bBin* SB55 = new RA2bBin("Njet9+-Nb1-MHT1-HT1", 9, 99, 1, 1, 200., 500., 500., 800., true);
   RA2bBin* SB56 = new RA2bBin("Njet9+-Nb1-MHT1-HT2", 9, 99, 1, 1, 200., 500., 800., 1200., true);
   RA2bBin* SB57 = new RA2bBin("Njet9+-Nb1-MHT1-HT3", 9, 99, 1, 1, 200., 500., 1200., 9999., true);
   RA2bBin* SB58 = new RA2bBin("Njet9+-Nb1-MHT2-HT1", 9, 99, 1, 1, 500., 750., 500., 1200., true);
   RA2bBin* SB59 = new RA2bBin("Njet9+-Nb1-MHT2-HT23", 9, 99, 1, 1, 500., 750., 1200., 9999., true);
   RA2bBin* SB60 = new RA2bBin("Njet9+-Nb1-MHT2-HT123", 9, 99, 1, 1, 750., 9999., 800., 9999., true);
   
   RA2bBin* SB61 = new RA2bBin("Njet9+-Nb2-MHT1-HT1", 9, 99, 2, 2, 200., 500., 500., 800., true);
   RA2bBin* SB62 = new RA2bBin("Njet9+-Nb2-MHT1-HT2", 9, 99, 2, 2, 200., 500., 800., 1200., true);
   RA2bBin* SB63 = new RA2bBin("Njet9+-Nb2-MHT1-HT3", 9, 99, 2, 2, 200., 500., 1200., 9999., true);
   RA2bBin* SB64 = new RA2bBin("Njet9+-Nb2-MHT2-HT1", 9, 99, 2, 2, 500., 750., 500., 1200., true);
   RA2bBin* SB65 = new RA2bBin("Njet9+-Nb2-MHT2-HT23", 9, 99, 2, 2, 500., 750., 1200., 9999., true);
   RA2bBin* SB66 = new RA2bBin("Njet9+-Nb2-MHT2-HT123", 9, 99, 2, 2, 750., 9999., 800., 9999., true);
   
   RA2bBin* SB67 = new RA2bBin("Njet9+-Nb3+-MHT1-HT1", 9, 99, 3, 9, 200., 500., 500., 800., true);
   RA2bBin* SB68 = new RA2bBin("Njet9+-Nb3+-MHT1-HT2", 9, 99, 3, 9, 200., 500., 800., 1200., true);
   RA2bBin* SB69 = new RA2bBin("Njet9+-Nb3+-MHT1-HT3", 9, 99, 3, 9, 200., 500., 1200., 9999., true);
   RA2bBin* SB70 = new RA2bBin("Njet9+-Nb3+-MHT2-HT1", 9, 99, 3, 9, 500., 750., 500., 1200., true);
   RA2bBin* SB71 = new RA2bBin("Njet9+-Nb3+-MHT2-HT23", 9, 99, 3, 9, 500., 750., 1200., 9999., true);
   RA2bBin* SB72 = new RA2bBin("Njet9+-Nb3+-MHT2-HT123", 9, 99, 3, 9, 750., 9999., 800., 9999., true);
   
   SB.push_back(SB1);
   SB.push_back(SB2);
   SB.push_back(SB3);
   SB.push_back(SB4);
   SB.push_back(SB5);
   SB.push_back(SB6);
   SB.push_back(SB7);
   SB.push_back(SB8);
   SB.push_back(SB9);
   SB.push_back(SB10);
   SB.push_back(SB11);
   SB.push_back(SB12);
   SB.push_back(SB13);
   SB.push_back(SB14);
   SB.push_back(SB15);
   SB.push_back(SB16);
   SB.push_back(SB17);
   SB.push_back(SB18);
   SB.push_back(SB19);
   SB.push_back(SB20);
   SB.push_back(SB21);
   SB.push_back(SB22);
   SB.push_back(SB23);
   SB.push_back(SB24);
   SB.push_back(SB25);
   SB.push_back(SB26);
   SB.push_back(SB27);
   SB.push_back(SB28);
   SB.push_back(SB29);
   SB.push_back(SB30);
   SB.push_back(SB31);
   SB.push_back(SB32);
   SB.push_back(SB33);
   SB.push_back(SB34);
   SB.push_back(SB35);
   SB.push_back(SB36);
   SB.push_back(SB37);
   SB.push_back(SB38);
   SB.push_back(SB39);
   SB.push_back(SB40);
   SB.push_back(SB41);
   SB.push_back(SB42);
   SB.push_back(SB43);
   SB.push_back(SB44);
   SB.push_back(SB45);
   SB.push_back(SB46);
   SB.push_back(SB47);
   SB.push_back(SB48);
   SB.push_back(SB49);
   SB.push_back(SB50);
   SB.push_back(SB51);
   SB.push_back(SB52);
   SB.push_back(SB53);
   SB.push_back(SB54);
   SB.push_back(SB55);
   SB.push_back(SB56);
   SB.push_back(SB57);
   SB.push_back(SB58);
   SB.push_back(SB59);
   SB.push_back(SB60);
   SB.push_back(SB61);
   SB.push_back(SB62);
   SB.push_back(SB63);
   SB.push_back(SB64);
   SB.push_back(SB65);
   SB.push_back(SB66);
   SB.push_back(SB67);
   SB.push_back(SB68);
   SB.push_back(SB69);
   SB.push_back(SB70);
   SB.push_back(SB71);
   SB.push_back(SB72);
   
   yields = new TH1F("expectation", "expectation", SB.size(), 0, SB.size());
   yields->Sumw2();
   for (int ibin=1; ibin <= (int) yields->GetNbinsX(); ibin++){
      yields->GetXaxis()->SetBinLabel(ibin, SB.at(ibin-1)->BinName.c_str());
   }
   
}

void Expectation::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
   
   TString option = GetOption();
   
}

Bool_t Expectation::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either Expectation::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
   
   fReader.SetEntry(entry);
   
   if( TotEvents%100000 == 0 ) std::cout << "event (expectation): " << TotEvents << '\n';
   ++TotEvents;

   bool mdp = DeltaPhiCut();
   if( *HT > 500. && *MHT > 200. && *NJets >= 4 && mdp && *Ntries == -2) {
      
      for (int it = 0; it < (int) SB.size(); ++it){
         if (SB.at(it)->Fill((int) *NJets, (int) *BTags, (double) *MHT, (double) *HT, mdp))
         {
            yields->Fill(it, *Weight);
         }
      }
      
   }

   return kTRUE;

}

void Expectation::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
   
}

void Expectation::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   
}

bool Expectation::DeltaPhiCut()
{
   bool deltaPhiCut = true;
   if( *NJets == 2 ) {
      if( DeltaPhi[0] < 0.5 ||
         DeltaPhi[1] < 0.5 ) deltaPhiCut = false;
   }
   if( *NJets >= 3 ) {
      if( DeltaPhi[0] < 0.5 ||
         DeltaPhi[1] < 0.5 ||
         DeltaPhi[2] < 0.3 ) deltaPhiCut = false;
   }
   
   return deltaPhiCut;
}

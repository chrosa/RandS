//
//  RA2bBin.cpp
//
//
//  Created by Christian Sander on 02/06/15.
//
//

#include "RA2bBin.h"

RA2bBin::RA2bBin(std::string BinName_, int NJets_low_, int NJets_high_, int NBJets_low_, int NBJets_high_, double MHT_low_, double MHT_high_, double HT_low_, double HT_high_, bool DeltaPhiMin_) {
   
   BinName = BinName_;
   MHT_low = MHT_low_;
   MHT_high = MHT_high_;
   HT_low = HT_low_;
   HT_high = HT_high_;
   NJets_low = NJets_low_;
   NJets_high = NJets_high_;
   NBJets_low = NBJets_low_;
   NBJets_high = NBJets_high_;
   DeltaPhiMin = DeltaPhiMin_;
   
}

bool RA2bBin::Fill(int NJets, int NBJets, double MHT, double HT, bool minDeltaPhi){
   if (HT < HT_low || HT >= HT_high) return false;
   if (MHT < MHT_low || MHT >= MHT_high) return false;
   if (NJets < NJets_low || NJets > NJets_high) return false;
   if (NBJets < NBJets_low || NBJets > NBJets_high) return false;
   if (minDeltaPhi != DeltaPhiMin) return false;
   return true;
}

RA2bBin::~RA2bBin(){
   
}


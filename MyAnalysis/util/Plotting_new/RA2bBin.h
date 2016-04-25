//
//  RA2bBin.h
//
//
//  Created by Christian Sander on 02/06/15.
//
//

#ifndef ____RA2bBin__
#define ____RA2bBin__

#include <stdio.h>
#include <vector>
#include <string>
#include <numeric>
#include <cmath>
#include <iostream>

using namespace std;

class RA2bBin{
   
public:
   
   RA2bBin(std::string, int, int, int, int, double, double, double, double, bool);
   ~RA2bBin();
   
   bool Fill(int, int, double, double, bool);
   
   std::string BinName;
   
private:
   
   double MHT_low, MHT_high;
   double HT_low, HT_high;
   int NJets_low, NJets_high;
   int NBJets_low, NBJets_high;
   bool DeltaPhiMin;
   
};

#endif /* defined(____RA2bBin__) */

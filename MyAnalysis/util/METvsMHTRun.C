//
//  METvsMHTRun.c
//
//
//  Created by Christian Sander on 11/05/16.
//
//

#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TFile.h>
#include <TChain.h>
#include <TSelector.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

int METvsMHTRun()
{

    TFile *f = new TFile("sample_AOD.root");
    TTree *t = (TTree *) f->Get("EventTree");
    t->Process("METvsMHTSelector.C+");

    return 1;

}

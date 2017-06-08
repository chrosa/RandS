//
//  RunMyResolution.C
//
//
//  Created by Christian Sander on 20/10/16.
//
//

#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TFile.h>
#include <TChain.h>
#include <TSelector.h>

#include <fstream>

#include "MyResolution.h"

int main()
{

    MyResolution* res = new MyResolution();
    TChain* res_chain = new TChain("EventTree");

    // ------------------------------------------------------------------- //

    std::ifstream myfile ("filelist_mc_v2.txt");

    std::string root_file;
    if (myfile.is_open()) {
        while( myfile.good() ) {
            std::getline (myfile,root_file);
            std::cout << root_file << std::endl;
            if (root_file.length() > 0) {
                TString path = root_file;
                res_chain->Add(path);
            }
        }
        myfile.close();
    }

	Long64_t nentries = res_chain->GetEntries();
	std::cout << "Number of total events: " << nentries << std::endl; 
	//Int_t cachesize = 30*1024*1024;
	//res_chain->SetCacheSize(cachesize);
	//res_chain->AddBranchToCache("*",kTRUE); //<<< add all branches to the cache


    // ------------------------------------------------------------------- //

    res_chain->Process(res);

    // ------------------------------------------------------------------- //

    return 1;

}

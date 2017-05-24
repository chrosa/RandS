#define METvsMHTSelector_cxx
// The class definition in METvsMHTSelector.h has been generated automatically
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
// root> T->Process("METvsMHTSelector.C")
// root> T->Process("METvsMHTSelector.C","some options")
// root> T->Process("METvsMHTSelector.C+")
//


#include "METvsMHTSelector.h"
#include <TH2.h>
#include <TStyle.h>

void METvsMHTSelector::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

    h_MHTtruerebPhiRes_vs_MHTrebMinusMET = new TH2F("h_MHTtruerebPhiRes_vs_MHTrebMinusMET","h_MHTtruerebPhiRes_vs_MHTrebMinusMET",60,-3.,3.,10,0.,100.);
    h_MHTtruerebPtRes_vs_MHTrebMinusMET = new TH2F("h_MHTtruerebPtRes_vs_MHTrebMinusMET","h_MHTtruerebPtRes_vs_MHTrebMinusMET",40,0.,2.,10,0.,100.);
    h_MHTtruerebPt_vs_MHTrebMinusMET = new TH2F("h_MHTtruerebPt_vs_MHTrebMinusMET","h_MHTtruerebPt_vs_MHTrebMinusMET",40,0.,100.,10,0.,100.);
	h_MHT_vs_MET = new TH2F("h_MHT_vs_MET","h_MHT_vs_MET",50,0.,500.,50,0.,500.);

    // NTulpleMaker_v1
	//AvailableEvents[361022] = 1993647;
	//AvailableEvents[361023] = 7514494;
	//AvailableEvents[361024] = 811300;
	//AvailableEvents[361025] = 570000;
	//AvailableEvents[361026] = 1224000;

    // test
	AvailableEvents[361023] = 10000;
	AvailableEvents[361024] = 10000;
	AvailableEvents[361025] = 9800;
	AvailableEvents[361026] = 10000;

}

void METvsMHTSelector::SlaveBegin(TTree * /*tree*/)
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

}

Bool_t METvsMHTSelector::Process(Long64_t entry)
{
    // The Process() function is called for each entry in the tree (or possibly
    // keyed object in the case of PROOF) to be processed. The entry argument
    // specifies which entry in the currently loaded tree is to be processed.
    // When processing keyed objects with PROOF, the object is already loaded
    // and is available via the fObject pointer.
    //
    // This function should contain the \"body\" of the analysis. It can contain
    // simple or elaborate selection criteria, run algorithms on the data
    // of the event and typically fill histograms.
    //
    // The processing can be stopped by calling Abort().
    //
    // Use fStatus to set the return value of TTree::Process().
    //
    // The return value is currently not used.

    fReader.SetLocalEntry(entry);

    NEvents += 1;
    //if (NEvents%1000 == 0) std::cout << NEvents << " processed: " << 100*NEvents/NTotEvents << "% done" << std::endl;
    if (NEvents%1000 == 0) std::cout << NEvents << " processed!" << std::endl;

    //std::cout << "Weight: " << *Weight << std::endl;
    double eventWeight = *Weight;
    eventWeight *= lumi_ / AvailableEvents[*DatasetID];
    
    eventWeight = 1.;

    std::vector<MyJet> Jets_gen;
    std::vector<MyJet> Jets_rec;
    std::vector<TLorentzVector> recoElectrons;
    std::vector<TLorentzVector> recoMuons;
    std::vector<TLorentzVector> recoPhotons;

    int NJets = JetPt->size();
    //std::cout << "Njets: " << NJets << std::endl;

    for (int i = 0; i < NJets; ++i) {
        float pt = JetPt->at(i);
        float eta = JetEta->at(i);
        float phi = JetPhi->at(i);
        float m = JetM->at(i);
        float jvt = JetJVT->at(i);
        bool btag = JetBtag->at(i);
        bool good = JetGood->at(i);
        MyJet jet;
        jet.momentum.SetPtEtaPhiM(pt, eta, phi,m);
        jet.jvt = jvt;
        jet.btag = btag;
        jet.good = good;

		//if (jet.momentum.Pt() > 20.){
		//	if ( !jet.good ) {
        //        std::cout << "Reject event because of bad jet!" << std::endl;
        //        return 1;				
		//	}
		//}

        if (jet.momentum.Pt() > 20. && jet.momentum.Pt() < 60. && fabs(jet.momentum.Eta()) < 2.4) {
            if ( !jet.good && jet.jvt > jvtcut_ ) {
                std::cout << "Reject event because of bad central jet!" << std::endl;
                return 1;
            }
        }
        if ( jet.momentum.Pt() > 20. && jet.momentum.Pt() < 60. && fabs(jet.momentum.Eta()) >= 2.4) {
            if ( !jet.good ) {
                std::cout << "Reject event because of bad forward jet!" << std::endl;
                return 1;
            }
        }
        if (jet.momentum.Pt() > 60.) {
            if ( !jet.good ) {
                std::cout << "Reject event because of bad high pT jet!" << std::endl;
                return 1;
            }
        }

        Jets_rec.push_back(jet);

    }

    int NGenJets = GenJetPt->size();
    //std::cout << "NGenjets: " << NGenJets << std::endl;

    for (int i = 0; i < NGenJets; ++i) {
        float pt = GenJetPt->at(i);
        float eta = GenJetEta->at(i);
        float phi = GenJetPhi->at(i);
        float m = GenJetM->at(i);
        bool btag = GenJetBtag->at(i);
        MyJet genjet;
        genjet.momentum.SetPtEtaPhiM(pt, eta, phi, m);
        genjet.btag = btag;
        genjet.jvt = 1.;
        genjet.good = true;
        Jets_gen.push_back(genjet);
    }


    int NMuons = MuonPt->size();
    //std::cout << "NMuons: " << NMuons << std::endl;
    //if ( NMuons > 0 ) return 1;

    for (int i = 0; i < NMuons; ++i) {
        float pt = MuonPt->at(i);
        float eta = MuonEta->at(i);
        float phi = MuonPhi->at(i);
        if (MuonIsBad->at(i)) {
            std::cout << "Reject event because of bad muon!" << std::endl;
            return 1;
        }
        if (MuonIsSignal->at(i) && MuonIsIso->at(i)) {
            std::cout << "Reject event because of isolated muon!" << std::endl;
        }
        if (MuonIsSignal->at(i)) {
            TLorentzVector muon(0.,0.,0.,0.);
            muon.SetPtEtaPhiM(pt, eta, phi, 0.1057);
            recoMuons.push_back(muon);
        }
    }

    int NElectrons = ElePt->size();
    //std::cout << "NElectrons: " << NElectrons << std::endl;
    if ( NElectrons > 0 ) {
        std::cout << "Reject event because of isolated electon!" << std::endl;
        return 1;
    }
    for (int i = 0; i < NElectrons; ++i) {
        float pt = ElePt->at(i);
        float eta = EleEta->at(i);
        float phi = ElePhi->at(i);
        TLorentzVector electron(0.,0.,0.,0.);
        electron.SetPtEtaPhiM(pt, eta, phi, 0.000511);
        recoElectrons.push_back(electron);
    }

    TLorentzVector recoMHTreb(0., 0., 0., 0.);
    for (vector<MyJet>::iterator it = Jets_rec.begin(); it != Jets_rec.end(); ++it) {
        if (it->momentum.Pt() > rebalancedJetPt_ && (it->momentum.Pt() > 60. || it->jvt > jvtcut_ || fabs(it->momentum.Eta()) > 2.4) ) {
        //if (it->momentum.Pt() > rebalancedJetPt_ ) {
            recoMHTreb -= it->momentum;
        }
    }
    
    TLorentzVector MET;
    MET.SetPtEtaPhiM(*MET_pt, 0, *MET_phi, 0);
    TLorentzVector METmu;
    METmu.SetPtEtaPhiM(*METmu_pt, 0, *METmu_phi, 0);
    TLorentzVector METtrack;
    METtrack.SetPtEtaPhiM(*METtrack_pt, 0, *METtrack_phi, 0);
    TLorentzVector METsoft = MET - METmu - recoMHTreb;
    TLorentzVector trueMHTreb;
    trueMHTreb.SetPtEtaPhiM(*TrueMHT_pt, 0, *TrueMHT_phi, 0);
    
    h_MHTtruerebPtRes_vs_MHTrebMinusMET->Fill(trueMHTreb.Pt()/METsoft.Pt(), (METsoft).Pt(), eventWeight);
    h_MHTtruerebPhiRes_vs_MHTrebMinusMET->Fill(trueMHTreb.DeltaPhi(METsoft), (METsoft).Pt(), eventWeight);
    h_MHTtruerebPhiRes_vs_MHTrebMinusMET->Fill(-trueMHTreb.DeltaPhi(METsoft), (METsoft).Pt(), eventWeight);
    h_MHTtruerebPt_vs_MHTrebMinusMET->Fill(trueMHTreb.Pt(), (METsoft).Pt(), eventWeight);
    h_MHT_vs_MET->Fill(recoMHTreb.Pt(), MET.Pt(), eventWeight);

    return kTRUE;
}

void METvsMHTSelector::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.

}

void METvsMHTSelector::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.
    
    TFile file("METsoft_resolutions.root","RECREATE");
	
    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
    c1->cd();
    c1->SetRightMargin(0.15);
    //c1->SetLogz(1);
    h_MHTtruerebPtRes_vs_MHTrebMinusMET->SetStats(0);
	//h_MHTtruerebPtRes_vs_MHTrebMinusMET->Smooth();
    h_MHTtruerebPtRes_vs_MHTrebMinusMET->GetXaxis()->SetTitle("MET^{soft}_{true}/(MET-MHT)");
    h_MHTtruerebPtRes_vs_MHTrebMinusMET->GetYaxis()->SetTitle("MET-MHT");
    h_MHTtruerebPtRes_vs_MHTrebMinusMET->GetZaxis()->SetTitle("events");
    h_MHTtruerebPtRes_vs_MHTrebMinusMET->GetXaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPtRes_vs_MHTrebMinusMET->GetYaxis()->SetTitleOffset(1.5);
    h_MHTtruerebPtRes_vs_MHTrebMinusMET->GetZaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPtRes_vs_MHTrebMinusMET->Draw("COLZ");
    h_MHTtruerebPtRes_vs_MHTrebMinusMET->Write();
	c1->SaveAs("h_MHTtruerebPtRes_vs_MHTrebMinusMET.pdf");

    TCanvas *c2 = new TCanvas("c2", "c2", 1200, 800);
    c2->cd();
    c2->SetRightMargin(0.15);
    //c2->SetLogz(1);
    h_MHTtruerebPt_vs_MHTrebMinusMET->SetStats(0);
	//h_MHTtruerebPt_vs_MHTrebMinusMET->Smooth();
    h_MHTtruerebPt_vs_MHTrebMinusMET->GetXaxis()->SetTitle("MET^{soft}_{true}");
    h_MHTtruerebPt_vs_MHTrebMinusMET->GetYaxis()->SetTitle("MET-MHT");
    h_MHTtruerebPt_vs_MHTrebMinusMET->GetZaxis()->SetTitle("events");
    h_MHTtruerebPt_vs_MHTrebMinusMET->GetXaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPt_vs_MHTrebMinusMET->GetYaxis()->SetTitleOffset(1.5);
    h_MHTtruerebPt_vs_MHTrebMinusMET->GetZaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPt_vs_MHTrebMinusMET->Draw("COLZ");
    h_MHTtruerebPt_vs_MHTrebMinusMET->Write();
	c2->SaveAs("h_MHTtruerebPt_vs_MHTrebMinusMET.pdf");

    TCanvas *c3 = new TCanvas("c3", "c3", 1200, 800);
    c3->cd();
    c3->SetRightMargin(0.15);
    //c3->SetLogz(1);
    h_MHTtruerebPhiRes_vs_MHTrebMinusMET->SetStats(0);
    //h_MHTtruerebPhiRes_vs_MHTrebMinusMET->Smooth();
    h_MHTtruerebPhiRes_vs_MHTrebMinusMET->GetXaxis()->SetTitle("#phi(MET^{soft}_{true}-(MET-MHT))");
    h_MHTtruerebPhiRes_vs_MHTrebMinusMET->GetYaxis()->SetTitle("MET-MHT");
    h_MHTtruerebPhiRes_vs_MHTrebMinusMET->GetZaxis()->SetTitle("events");
    h_MHTtruerebPhiRes_vs_MHTrebMinusMET->GetXaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPhiRes_vs_MHTrebMinusMET->GetYaxis()->SetTitleOffset(1.5);
    h_MHTtruerebPhiRes_vs_MHTrebMinusMET->GetZaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPhiRes_vs_MHTrebMinusMET->Draw("COLZ");
    h_MHTtruerebPhiRes_vs_MHTrebMinusMET->Write();
	c3->SaveAs("h_MHTtruerebPhiRes_vs_MHTrebMinusMET.pdf");
    
    TCanvas *c4 = new TCanvas("c4", "c4", 1200, 800);
    c4->cd();
    c4->SetRightMargin(0.15);
    //c4->SetLogz(1);
    h_MHT_vs_MET->SetStats(0);
    //h_MHT_vs_MET->Smooth();
    h_MHT_vs_MET->GetXaxis()->SetTitle("MHT");
    h_MHT_vs_MET->GetYaxis()->SetTitle("MET");
    h_MHT_vs_MET->GetZaxis()->SetTitle("events");
    h_MHT_vs_MET->GetXaxis()->SetTitleOffset(1.3);
    h_MHT_vs_MET->GetYaxis()->SetTitleOffset(1.5);
    h_MHT_vs_MET->GetZaxis()->SetTitleOffset(1.3);
    h_MHT_vs_MET->Draw("COLZ");
    h_MHT_vs_MET->Write();
	c4->SaveAs("h_MHT_vs_MET.pdf");

    file.Close();

}

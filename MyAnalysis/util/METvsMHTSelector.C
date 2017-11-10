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

    h_MHTtruerebPhiRes_vs_MET = new TH2F("h_MHTtruerebPhiRes_vs_MET","h_MHTtruerebPhiRes_vs_MET",60,-3.,3.,10,0.,100.);
    h_MHTtruerebPt_vs_MET = new TH2F("h_MHTtruerebPt_vs_MET","h_MHTtruerebPt_vs_MET",50,0.,50.,10,0.,100.);
    h_MHTtruerebPhiRes_vs_MHT = new TH2F("h_MHTtruerebPhiRes_vs_MHT","h_MHTtruerebPhiRes_vs_MHT",60,-3.,3.,10,0.,100.);
    h_MHTtruerebPt_vs_MHT = new TH2F("h_MHTtruerebPt_vs_MHT","h_MHTtruerebPt_vs_MHT",50,0.,50.,10,0.,100.);
    h_MHTtruerebPhiRes_vs_MHTnoJVT = new TH2F("h_MHTtruerebPhiRes_vs_MHTnoJVT","h_MHTtruerebPhiRes_vs_MHTnoJVT",60,-3.,3.,10,0.,100.);
    h_MHTtruerebPt_vs_MHTnoJVT = new TH2F("h_MHTtruerebPt_vs_MHTnoJVT","h_MHTtruerebPt_vs_MHTnoJVT",50,0.,50.,10,0.,100.);
    h_MHTtruerebPhiRes_vs_METsoft = new TH2F("h_MHTtruerebPhiRes_vs_METsoft","h_MHTtruerebPhiRes_vs_METsoft",60,-3.,3.,10,0.,25.);
    h_MHTtruerebPt_vs_METsoft = new TH2F("h_MHTtruerebPt_vs_METsoft","h_MHTtruerebPt_vs_METsoft",50,0.,50.,10,0.,25.);
    h_MHT_vs_MET = new TH2F("h_MHT_vs_MET","h_MHT_vs_MET",100,0.,100.,100,0.,100.);
    h_MHTnoJVT_vs_MET = new TH2F("h_MHTnoJVT_vs_MET","h_MHTnoJVT_vs_MET",100,0.,100.,100,0.,100.);
    h_MHTtruerebPhiRes_vs_MHTtruereb = new TH2F("h_MHTtruerebPhiRes_vs_MHTtruereb","h_MHTtruerebPhiRes_vs_MHTtruereb",60,-3.,3.,10,0.,50.);

    // NTulpleMaker_v1
	//AvailableEvents[361022] = 1993647;
	//AvailableEvents[361023] = 7514494;
	//AvailableEvents[361024] = 811300;
	//AvailableEvents[361025] = 570000;
	//AvailableEvents[361026] = 1224000;

    // test
	AvailableEvents[361023] = 10000;
	AvailableEvents[361024] = 10000;
	AvailableEvents[361025] = 10000;
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
    //if (*DatasetID != 361023) eventWeight = 0.;

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

		if (jet.momentum.Pt() > 20.){
			if ( !jet.good ) {
                std::cout << "Reject event because of bad jet!" << std::endl;
                return 1;				
			}
		}

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

    for (int i = 0; i < NMuons; ++i) {
        float pt = MuonPt->at(i);
        float eta = MuonEta->at(i);
        float phi = MuonPhi->at(i);
        if (MuonIsBad->at(i)) {
            std::cout << "Reject event because of bad muon!" << std::endl;
            return 1;
        }
        if (MuonIsSignal->at(i)) {
            std::cout << "Reject event because of isolated muon!" << std::endl;
            return 1;
        }
        if (MuonIsSignal->at(i)) {
            TLorentzVector muon(0.,0.,0.,0.);
            muon.SetPtEtaPhiM(pt, eta, phi, 0.1057);
            recoMuons.push_back(muon);
        }
    }

    int NElectrons = ElePt->size();
    //std::cout << "NElectrons: " << NElectrons << std::endl;

    for (int i = 0; i < NElectrons; ++i) {
        float pt = ElePt->at(i);
        float eta = EleEta->at(i);
        float phi = ElePhi->at(i);
        if (EleIsSignal->at(i)) {
            std::cout << "Reject event because of isolated electron!" << std::endl;
            return 1;
        }
        TLorentzVector electron(0.,0.,0.,0.);
        electron.SetPtEtaPhiM(pt, eta, phi, 0.000511);
        recoElectrons.push_back(electron);
    }

    TLorentzVector MHT(0., 0., 0., 0.);
    TLorentzVector MHTnoJVT(0., 0., 0., 0.);
    for (vector<MyJet>::iterator it = Jets_rec.begin(); it != Jets_rec.end(); ++it) {
        if (it->momentum.Pt() > rebalancedJetPt_) {
			MHTnoJVT -= it->momentum;
			if (it->momentum.Pt() > 60. || it->jvt > jvtcut_ || fabs(it->momentum.Eta()) > 2.4) {
				MHT -= it->momentum;
			}
        }
    }
    
    TLorentzVector MET;
    MET.SetPtEtaPhiM(*MET_pt, 0, *MET_phi, 0);

    TLorentzVector METtrack;
    METtrack.SetPtEtaPhiM(*METtrack_pt, 0, *METtrack_phi, 0);

    TLorentzVector METsoft = METtrack;

    TLorentzVector trueMHTreb;
    trueMHTreb.SetPtEtaPhiM(*TrueMHT_pt, 0, *TrueMHT_phi, 0);
    
	if (MET.Pt() > 0.) {
		h_MHTtruerebPhiRes_vs_MHTtruereb->Fill(trueMHTreb.DeltaPhi(MET),trueMHTreb.Pt(),eventWeight);
		h_MHTtruerebPhiRes_vs_MET->Fill(trueMHTreb.DeltaPhi(MET),MET.Pt(),eventWeight);
		h_MHTtruerebPt_vs_MET->Fill(trueMHTreb.Pt(),MET.Pt(),eventWeight);
		h_MHTtruerebPhiRes_vs_MHT->Fill(trueMHTreb.DeltaPhi(MHT),MHT.Pt(),eventWeight);
		h_MHTtruerebPt_vs_MHT->Fill(trueMHTreb.Pt(),MHT.Pt(),eventWeight);
		h_MHTtruerebPhiRes_vs_MHTnoJVT->Fill(trueMHTreb.DeltaPhi(MHTnoJVT),MHTnoJVT.Pt(),eventWeight);
		h_MHTtruerebPt_vs_MHTnoJVT->Fill(trueMHTreb.Pt(),MHTnoJVT.Pt(),eventWeight);
		h_MHTtruerebPhiRes_vs_METsoft->Fill(trueMHTreb.DeltaPhi(METsoft),METsoft.Pt(),eventWeight);
		h_MHTtruerebPt_vs_METsoft->Fill(trueMHTreb.Pt(),METsoft.Pt(),eventWeight);
		h_MHT_vs_MET->Fill(MHT.Pt(),MET.Pt(),eventWeight);
		h_MHTnoJVT_vs_MET->Fill(MHTnoJVT.Pt(),MET.Pt(),eventWeight);
	}

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
    
    TFile file("METsoft_resolutions_new.root","RECREATE");
	
    TCanvas *c1a = new TCanvas("c1a", "c1a", 1200, 800);
    c1a->cd();
    c1a->SetRightMargin(0.15);
    //c1a->SetLogz(1);
    h_MHTtruerebPhiRes_vs_MET->SetStats(0);
	//h_MHTtruerebPhiRes_vs_MET->Smooth();
    h_MHTtruerebPhiRes_vs_MET->GetXaxis()->SetTitle("#Delta#phi(MET^{soft}_{true},MET)");
    h_MHTtruerebPhiRes_vs_MET->GetYaxis()->SetTitle("MET");
    h_MHTtruerebPhiRes_vs_MET->GetZaxis()->SetTitle("events");
    h_MHTtruerebPhiRes_vs_MET->GetXaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPhiRes_vs_MET->GetYaxis()->SetTitleOffset(1.5);
    h_MHTtruerebPhiRes_vs_MET->GetZaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPhiRes_vs_MET->Draw("COLZ");
    h_MHTtruerebPhiRes_vs_MET->Write();
	c1a->SaveAs("h_MHTtruerebPhiRes_vs_MET.pdf");

    TCanvas *c2a = new TCanvas("c2a", "c2a", 1200, 800);
    c2a->cd();
    c2a->SetRightMargin(0.15);
    //c2a->SetLogz(1);
    h_MHTtruerebPt_vs_MET->SetStats(0);
	//h_MHTtruerebPt_vs_MET->Smooth();
    h_MHTtruerebPt_vs_MET->GetXaxis()->SetTitle("MET^{soft}_{true}");
    h_MHTtruerebPt_vs_MET->GetYaxis()->SetTitle("MET");
    h_MHTtruerebPt_vs_MET->GetZaxis()->SetTitle("events");
    h_MHTtruerebPt_vs_MET->GetXaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPt_vs_MET->GetYaxis()->SetTitleOffset(1.5);
    h_MHTtruerebPt_vs_MET->GetZaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPt_vs_MET->Draw("COLZ");
    h_MHTtruerebPt_vs_MET->Write();
	c2a->SaveAs("h_MHTtruerebPt_vs_MET.pdf");

    TCanvas *c1b = new TCanvas("c1b", "c1b", 1200, 800);
    c1b->cd();
    c1b->SetRightMargin(0.15);
    //c1b->SetLogz(1);
    h_MHTtruerebPhiRes_vs_MHT->SetStats(0);
	//h_MHTtruerebPhiRes_vs_MHT->Smooth();
    h_MHTtruerebPhiRes_vs_MHT->GetXaxis()->SetTitle("#Delta#phi(MET^{soft}_{true},MHT)");
    h_MHTtruerebPhiRes_vs_MHT->GetYaxis()->SetTitle("MHT");
    h_MHTtruerebPhiRes_vs_MHT->GetZaxis()->SetTitle("events");
    h_MHTtruerebPhiRes_vs_MHT->GetXaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPhiRes_vs_MHT->GetYaxis()->SetTitleOffset(1.5);
    h_MHTtruerebPhiRes_vs_MHT->GetZaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPhiRes_vs_MHT->Draw("COLZ");
    h_MHTtruerebPhiRes_vs_MHT->Write();
	c1b->SaveAs("h_MHTtruerebPhiRes_vs_MHT.pdf");

    TCanvas *c2b = new TCanvas("c2b", "c2b", 1200, 800);
    c2b->cd();
    c2b->SetRightMargin(0.15);
    //c2b->SetLogz(1);
    h_MHTtruerebPt_vs_MHT->SetStats(0);
	//h_MHTtruerebPt_vs_MHT->Smooth();
    h_MHTtruerebPt_vs_MHT->GetXaxis()->SetTitle("MET^{soft}_{true}");
    h_MHTtruerebPt_vs_MHT->GetYaxis()->SetTitle("MHT");
    h_MHTtruerebPt_vs_MHT->GetZaxis()->SetTitle("events");
    h_MHTtruerebPt_vs_MHT->GetXaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPt_vs_MHT->GetYaxis()->SetTitleOffset(1.5);
    h_MHTtruerebPt_vs_MHT->GetZaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPt_vs_MHT->Draw("COLZ");
    h_MHTtruerebPt_vs_MHT->Write();
	c2b->SaveAs("h_MHTtruerebPt_vs_MHT.pdf");

    TCanvas *c1c = new TCanvas("c1c", "c1c", 1200, 800);
    c1c->cd();
    c1c->SetRightMargin(0.15);
    //c1c->SetLogz(1);
    h_MHTtruerebPhiRes_vs_MHTnoJVT->SetStats(0);
	//h_MHTtruerebPhiRes_vs_MHTnoJVT->Smooth();
    h_MHTtruerebPhiRes_vs_MHTnoJVT->GetXaxis()->SetTitle("#Delta#phi(MET^{soft}_{true},MHTnoJVT)");
    h_MHTtruerebPhiRes_vs_MHTnoJVT->GetYaxis()->SetTitle("MHTnoJVT");
    h_MHTtruerebPhiRes_vs_MHTnoJVT->GetZaxis()->SetTitle("events");
    h_MHTtruerebPhiRes_vs_MHTnoJVT->GetXaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPhiRes_vs_MHTnoJVT->GetYaxis()->SetTitleOffset(1.5);
    h_MHTtruerebPhiRes_vs_MHTnoJVT->GetZaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPhiRes_vs_MHTnoJVT->Draw("COLZ");
    h_MHTtruerebPhiRes_vs_MHTnoJVT->Write();
	c1c->SaveAs("h_MHTtruerebPhiRes_vs_MHTnoJVT.pdf");

    TCanvas *c2c = new TCanvas("c2c", "c2c", 1200, 800);
    c2c->cd();
    c2c->SetRightMargin(0.15);
    //c2c->SetLogz(1);
    h_MHTtruerebPt_vs_MHTnoJVT->SetStats(0);
	//h_MHTtruerebPt_vs_MHTnoJVT->Smooth();
    h_MHTtruerebPt_vs_MHTnoJVT->GetXaxis()->SetTitle("MET^{soft}_{true}");
    h_MHTtruerebPt_vs_MHTnoJVT->GetYaxis()->SetTitle("MHTnoJVT");
    h_MHTtruerebPt_vs_MHTnoJVT->GetZaxis()->SetTitle("events");
    h_MHTtruerebPt_vs_MHTnoJVT->GetXaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPt_vs_MHTnoJVT->GetYaxis()->SetTitleOffset(1.5);
    h_MHTtruerebPt_vs_MHTnoJVT->GetZaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPt_vs_MHTnoJVT->Draw("COLZ");
    h_MHTtruerebPt_vs_MHTnoJVT->Write();
	c2c->SaveAs("h_MHTtruerebPt_vs_MHTnoJVT.pdf");

    TCanvas *c1d = new TCanvas("c1d", "c1d", 1200, 800);
    c1d->cd();
    c1d->SetRightMargin(0.15);
    //c1d->SetLogz(1);
    h_MHTtruerebPhiRes_vs_METsoft->SetStats(0);
	//h_MHTtruerebPhiRes_vs_METsoft->Smooth();
    h_MHTtruerebPhiRes_vs_METsoft->GetXaxis()->SetTitle("#Delta#phi(MET^{soft}_{true},MET^{soft})");
    h_MHTtruerebPhiRes_vs_METsoft->GetYaxis()->SetTitle("MET^{soft}");
    h_MHTtruerebPhiRes_vs_METsoft->GetZaxis()->SetTitle("events");
    h_MHTtruerebPhiRes_vs_METsoft->GetXaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPhiRes_vs_METsoft->GetYaxis()->SetTitleOffset(1.5);
    h_MHTtruerebPhiRes_vs_METsoft->GetZaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPhiRes_vs_METsoft->Draw("COLZ");
    h_MHTtruerebPhiRes_vs_METsoft->Write();
	c1d->SaveAs("h_MHTtruerebPhiRes_vs_METsoft.pdf");

    TCanvas *c2d = new TCanvas("c2d", "c2d", 1200, 800);
    c2d->cd();
    c2d->SetRightMargin(0.15);
    //c2d->SetLogz(1);
    h_MHTtruerebPt_vs_METsoft->SetStats(0);
	//h_MHTtruerebPt_vs_METsoft->Smooth();
    h_MHTtruerebPt_vs_METsoft->GetXaxis()->SetTitle("MET^{soft}_{true}");
    h_MHTtruerebPt_vs_METsoft->GetYaxis()->SetTitle("MET^{soft}");
    h_MHTtruerebPt_vs_METsoft->GetZaxis()->SetTitle("events");
    h_MHTtruerebPt_vs_METsoft->GetXaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPt_vs_METsoft->GetYaxis()->SetTitleOffset(1.5);
    h_MHTtruerebPt_vs_METsoft->GetZaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPt_vs_METsoft->Draw("COLZ");
    h_MHTtruerebPt_vs_METsoft->Write();
	c2d->SaveAs("h_MHTtruerebPt_vs_METsoft.pdf");

    TCanvas *c1e = new TCanvas("c1e", "c1e", 1200, 800);
    c1e->cd();
    c1e->SetRightMargin(0.15);
    //c1e->SetLogz(1);
    h_MHTtruerebPhiRes_vs_MHTtruereb->SetStats(0);
	//h_MHTtruerebPhiRes_vs_MHTtruereb->Smooth();
    h_MHTtruerebPhiRes_vs_MHTtruereb->GetXaxis()->SetTitle("#Delta#phi(MET^{soft}_{true},MET)");
    h_MHTtruerebPhiRes_vs_MHTtruereb->GetYaxis()->SetTitle("MET^{soft}_{true}");
    h_MHTtruerebPhiRes_vs_MHTtruereb->GetZaxis()->SetTitle("events");
    h_MHTtruerebPhiRes_vs_MHTtruereb->GetXaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPhiRes_vs_MHTtruereb->GetYaxis()->SetTitleOffset(1.5);
    h_MHTtruerebPhiRes_vs_MHTtruereb->GetZaxis()->SetTitleOffset(1.3);
    h_MHTtruerebPhiRes_vs_MHTtruereb->Draw("COLZ");
    h_MHTtruerebPhiRes_vs_MHTtruereb->Write();
	c1e->SaveAs("h_MHTtruerebPhiRes_vs_MHTtruereb.pdf");

    TCanvas *c3a = new TCanvas("c3a", "c3a", 1200, 800);
    c3a->cd();
    c3a->SetRightMargin(0.15);
    //c3a->SetLogz(1);
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
	c3a->SaveAs("h_MHT_vs_MET.pdf");

    TCanvas *c3b = new TCanvas("c3b", "c3b", 1200, 800);
    c3b->cd();
    c3b->SetRightMargin(0.15);
    //c3b->SetLogz(1);
    h_MHTnoJVT_vs_MET->SetStats(0);
	//h_MHTnoJVT_vs_MET->Smooth();
    h_MHTnoJVT_vs_MET->GetXaxis()->SetTitle("MHTnoJVT");
    h_MHTnoJVT_vs_MET->GetYaxis()->SetTitle("MET");
    h_MHTnoJVT_vs_MET->GetZaxis()->SetTitle("events");
    h_MHTnoJVT_vs_MET->GetXaxis()->SetTitleOffset(1.3);
    h_MHTnoJVT_vs_MET->GetYaxis()->SetTitleOffset(1.5);
    h_MHTnoJVT_vs_MET->GetZaxis()->SetTitleOffset(1.3);
    h_MHTnoJVT_vs_MET->Draw("COLZ");
    h_MHTnoJVT_vs_MET->Write();
	c3b->SaveAs("h_MHTnoJVT_vs_MET.pdf");

    file.Close();

}

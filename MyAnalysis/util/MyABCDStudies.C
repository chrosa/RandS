#define MyABCDStudies_cxx
// The class definition in MyABCDStudies.h has been generated automatically
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
// root> T->Process("MyABCDStudies.C")
// root> T->Process("MyABCDStudies.C","some options")
// root> T->Process("MyABCDStudies.C+")
//

#include "MyJet.h"
#include "MyElectron.h"
#include "MyPhoton.h"
#include "MyMuon.h"
#include "MyTau.h"
#include "MyABCDStudies.h"

void MyABCDStudies::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

    outputfile = new TFile("ABCDStudiesOutput.root","RECREATE");

    m_jvtcut = 0.59;
    m_lumi = 32900.;

    dPhijjMin_SR = 0.;
    dPhijjMax_SR = 1.8;
    dPhijjMin_CR = 1.8;
    dPhijjMax_CR = TMath::Pi();

    METMin_SR = 150.;
    METMax_SR = 9999.;
    METMin_CR = 80.;
    METMax_CR = 150.;

    dEtajjMin_SR = 4.8;
    dEtajjMin_CR = 4.8;

    dPhiJet1METMin_SR = 0.0;
    dPhiJet1METMin_CR = 0.0;

    dPhiJet2METMin_SR = 0.0;
    dPhiJet2METMin_CR = 0.0;

    //// Book 1d histograms

    std::cout << "Booking 1D histos" << std::endl;

    h_Jet1_Pt = new TH1F("h_Jet1_Pt", "h_Jet1_Pt", 100, 0., 4000.);
    h_Jet1_Pt->Sumw2();
    histos_1D.push_back(h_Jet1_Pt);

    h_Jet2_Pt = new TH1F("h_Jet2_Pt", "h_Jet2_Pt", 100, 0., 4000.);
    h_Jet2_Pt->Sumw2();
    histos_1D.push_back(h_Jet2_Pt);

    h_Jet3_Pt = new TH1F("h_Jet3_Pt", "h_Jet3_Pt", 100, 0., 1000.);
    h_Jet3_Pt->Sumw2();
    histos_1D.push_back(h_Jet3_Pt);

    h_Jet1_Eta = new TH1F("h_Jet1_Eta", "h_Jet1_Eta", 100, -5., 5.);
    h_Jet1_Eta->Sumw2();
    histos_1D.push_back(h_Jet1_Eta);

    h_Jet2_Eta = new TH1F("h_Jet2_Eta", "h_Jet2_Eta", 100, -5., 5.);
    h_Jet2_Eta->Sumw2();
    histos_1D.push_back(h_Jet2_Eta);

    h_Jet3_Eta = new TH1F("h_Jet3_Eta", "h_Jet3_Eta", 100, -5., 5.);
    h_Jet3_Eta->Sumw2();
    histos_1D.push_back(h_Jet3_Eta);

    h_Jet1_Phi = new TH1F("h_Jet1_Phi", "h_Jet1_Phi", 100, -TMath::Pi(), TMath::Pi());
    h_Jet1_Phi->Sumw2();
    histos_1D.push_back(h_Jet1_Phi);

    h_Jet2_Phi = new TH1F("h_Jet2_Phi", "h_Jet2_Phi", 100, -TMath::Pi(), TMath::Pi());
    h_Jet2_Phi->Sumw2();
    histos_1D.push_back(h_Jet2_Phi);

    h_Jet3_Phi = new TH1F("h_Jet3_Phi", "h_Jet3_Phi", 100, -TMath::Pi(), TMath::Pi());
    h_Jet3_Phi->Sumw2();
    histos_1D.push_back(h_Jet3_Phi);

    h_Jet1_DeltaPhi = new TH1F("h_Jet1_DeltaPhi", "h_Jet1_DeltaPhi", 100, 0., TMath::Pi());
    h_Jet1_DeltaPhi->Sumw2();
    histos_1D.push_back(h_Jet1_DeltaPhi);

    h_Jet2_DeltaPhi = new TH1F("h_Jet2_DeltaPhi", "h_Jet2_DeltaPhi", 100, 0., TMath::Pi());
    h_Jet2_DeltaPhi->Sumw2();
    histos_1D.push_back(h_Jet2_DeltaPhi);

    h_Jet3_DeltaPhi = new TH1F("h_Jet3_DeltaPhi", "h_Jet3_DeltaPhi", 100, 0., TMath::Pi());
    h_Jet3_DeltaPhi->Sumw2();
    histos_1D.push_back(h_Jet3_DeltaPhi);

    h_MET = new TH1F("h_MET", "h_MET", 100, 0., 1000.);
    h_MET->Sumw2();
    histos_1D.push_back(h_MET);

    h_METsig = new TH1F("h_METsig", "h_METsig", 100, 0., 20.);
    h_METsig->Sumw2();
    histos_1D.push_back(h_METsig);

    h_METsoft = new TH1F("h_METsoft", "h_METsoft", 100, 0., 250.);
    h_METsoft->Sumw2();
    histos_1D.push_back(h_METsoft);

    h_DeltaPhijj = new TH1F("h_DeltaPhijj", "h_DeltaPhijj", 100, 0., TMath::Pi());
    h_DeltaPhijj->Sumw2();
    histos_1D.push_back(h_DeltaPhijj);

    h_DeltaEtajj = new TH1F("h_DeltaEtajj", "h_DeltaEtajj", 100, 0., 10.);
    h_DeltaEtajj->Sumw2();
    histos_1D.push_back(h_DeltaEtajj);

    h_Mjj = new TH1F("h_Mjj", "h_Mjj", 100, 0., 5000.);
    h_Mjj->Sumw2();
    histos_1D.push_back(h_Mjj);

    //// Book 2d histograms

    std::cout << "Booking 2D histos" << std::endl;

    h_MET_vs_dPhi_Incl = new TH2F("h_MET_vs_dPhi_Incl", "h_MET_vs_dPhi_Incl", 100, 0., 500., 100, 0., TMath::Pi());
    h_MET_vs_dPhi_Incl->Sumw2();
    histos_2D.push_back(h_MET_vs_dPhi_Incl);

    h_MET_vs_dPhi_Mjj600 = new TH2F("h_MET_vs_dPhi_Mjj600", "h_MET_vs_dPhi_Mjj600", 100, 0., 500., 100, 0., TMath::Pi());
    h_MET_vs_dPhi_Mjj600->Sumw2();
    histos_2D.push_back(h_MET_vs_dPhi_Mjj600);

    h_MET_vs_dPhi_Mjj1000 = new TH2F("h_MET_vs_dPhi_Mjj1000", "h_MET_vs_dPhi_Mjj1000", 100, 0., 500., 100, 0., TMath::Pi());
    h_MET_vs_dPhi_Mjj1000->Sumw2();
    histos_2D.push_back(h_MET_vs_dPhi_Mjj1000);

    h_MET_vs_dPhi_Mjj1500 = new TH2F("h_MET_vs_dPhi_Mjj1500", "h_MET_vs_dPhi_Mjj1500", 100, 0., 500., 100, 0., TMath::Pi());
    h_MET_vs_dPhi_Mjj1500->Sumw2();
    histos_2D.push_back(h_MET_vs_dPhi_Mjj1500);

    h_MET_vs_dPhi_Mjj2000 = new TH2F("h_MET_vs_dPhi_Mjj2000", "h_MET_vs_dPhi_Mjj2000", 100, 0., 500., 100, 0., TMath::Pi());
    h_MET_vs_dPhi_Mjj2000->Sumw2();
    histos_2D.push_back(h_MET_vs_dPhi_Mjj2000);

    h_MET_vs_dEta_Incl = new TH2F("h_MET_vs_dEta_Incl", "h_MET_vs_dEta_Incl", 100, 0., 500., 100, 0., 10.);
    h_MET_vs_dEta_Incl->Sumw2();
    histos_2D.push_back(h_MET_vs_dEta_Incl);

    h_MET_vs_dEta_Mjj600 = new TH2F("h_MET_vs_dEta_Mjj600", "h_MET_vs_dEta_Mjj600", 100, 0., 500., 100, 0., 10.);
    h_MET_vs_dEta_Mjj600->Sumw2();
    histos_2D.push_back(h_MET_vs_dEta_Mjj600);

    h_MET_vs_dEta_Mjj1000 = new TH2F("h_MET_vs_dEta_Mjj1000", "h_MET_vs_dEta_Mjj1000", 100, 0., 500., 100, 0., 10.);
    h_MET_vs_dEta_Mjj1000->Sumw2();
    histos_2D.push_back(h_MET_vs_dEta_Mjj1000);

    h_MET_vs_dEta_Mjj1500 = new TH2F("h_MET_vs_dEta_Mjj1500", "h_MET_vs_dEta_Mjj1500", 100, 0., 500., 100, 0., 10.);
    h_MET_vs_dEta_Mjj1500->Sumw2();
    histos_2D.push_back(h_MET_vs_dEta_Mjj1500);

    h_MET_vs_dEta_Mjj2000 = new TH2F("h_MET_vs_dEta_Mjj2000", "h_MET_vs_dEta_Mjj2000", 100, 0., 500., 100, 0., 10.);
    h_MET_vs_dEta_Mjj2000->Sumw2();
    histos_2D.push_back(h_MET_vs_dEta_Mjj2000);

    std::cout << "Booking histos: DONE" << std::endl;

    //NTotEvents = fChain->GetEntries();

    //// Not very elegant! TODO: Store this info in and read from file

    // [v1]
    AvailableEvents[361022] = 1993647;
    AvailableEvents[361023] = 7724495;
    AvailableEvents[361024] = 7890000;
    AvailableEvents[361025] = 7977600;
    AvailableEvents[361026] = 1833400;

    std::cout << "Init: DONE" << std::endl;

}

void MyABCDStudies::SlaveBegin(TTree * /*tree*/)
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

}

Bool_t MyABCDStudies::Process(Long64_t entry)
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
    if (NEvents%100000 == 0) std::cout << NEvents << " processed!" << std::endl;


    if (ProcessedEvents.find(*DatasetID) == ProcessedEvents.end()) {
        ProcessedEvents[*DatasetID] = 1;
    } else {
        ProcessedEvents[*DatasetID] += 1;
    }

    //if ( *DatasetID != 361022 ) return 0;

    if ( *PrimaryVtx == 0 ) return 0;

    if (isinf(*Weight)) return 0;
    if (isnan(*Weight)) return 0;
    
    //std::cout << "Weight: " << *Weight << std::endl;
    double eventWeight = *Weight;
    if (isMC) eventWeight *= m_lumi / AvailableEvents[*DatasetID];
    //eventWeight = 1.; // for signal samples

    std::vector<MyJet> genJets;
    std::vector<MyJet> recoJets;
    std::vector<MyElectron> recoElectrons;
    std::vector<MyMuon> recoMuons;
    std::vector<MyTau> recoTaus;
    std::vector<MyPhoton> recoPhotons;

    //// Jets ////////////////

    int NJets = JetPt->size();
    //std::cout << "Njets: " << NJets << std::endl;

    for (int i = 0; i < NJets; ++i) {
        bool OR = JetPassOR->at(i);
        float pt = JetPt->at(i);
        float eta = JetEta->at(i);
        float phi = JetPhi->at(i);
        float m = JetM->at(i);
        float jvt = JetJVT->at(i);
        bool fjvt = JetFJVT->at(i);
        float sumpt = JetSumPtTracks->at(i);
        float tw = JetTrackWidth->at(i);
        unsigned short ntracks = JetNTracks->at(i);
        bool btag = JetBtag->at(i);
        bool good = JetGood->at(i);
        MyJet jet(pt, eta, phi,m);
        jet.SetJVT(jvt);
        if (fabs(jet.Eta()) > 2.4) jet.SetJVT(1.);
        jet.SetFJVT(fjvt);
        jet.SetSumPtTracks(sumpt);
        jet.SetTrackWidth(tw);
        jet.SetNTracks(ntracks);
        jet.SetBTag(btag);
        jet.SetJetID(good);
        jet.SetPassOR(OR);

        if (jet.Pt() > 20. && jet.Pt() < 60. && fabs(jet.Eta()) < 2.4) {
            if (jet.IsBad() && jet.IsNoPU(m_jvtcut)) {
                //std::cout << "Reject event because of bad central jet!" << std::endl;
                return 1;
            }
        }
        if (jet.Pt() > 20. && jet.Pt() < 60. && fabs(jet.Eta()) >= 2.4) {
            if (jet.IsBad()) {
                //std::cout << "Reject event because of bad forward jet!" << std::endl;
                return 1;
            }
        }
        if (jet.Pt() > 60.) {
            if (jet.IsBad()) {
                //std::cout << "Reject event because of bad high pT jet!" << std::endl;
                return 1;
            }
        }
        if ( jet.IsGood() && jet.Pt() > 20. && fabs(jet.Eta()) < 4.5 ) recoJets.push_back(jet);
        //if ( jet.IsGood() && jet.Pt() > 20. && fabs(jet.Eta()) < 4.5 && (jet.IsNoPU(m_jvtcut) || jet.Pt() > 60. || fabs(jet.Eta()) > 2.4 )) recoJets.push_back(jet);

    }
    GreaterByPt<MyJet> ptComparator_;
    std::sort(recoJets.begin(), recoJets.end(), ptComparator_);

    if (isMC) {
        int NGenJets = GenJetPt->size();
        //std::cout << "NGenjets: " << NGenJets << std::endl;

        for (int i = 0; i < NGenJets; ++i) {
            float pt = GenJetPt->at(i);
            float eta = GenJetEta->at(i);
            float phi = GenJetPhi->at(i);
            float m = GenJetM->at(i);
            bool btag = GenJetBtag->at(i);
            MyJet genjet(pt, eta, phi, m);
            genjet.SetBTag(btag);
            genJets.push_back(genjet);
        }
        std::sort(genJets.begin(), genJets.end(), ptComparator_);
    }

    //// Muons ////////////////

    int NMuons = MuonPt->size();
    //std::cout << "NMuons: " << NMuons << std::endl;

    for (int i = 0; i < NMuons; ++i) {
        bool OR = MuonPassOR->at(i);
        bool signal = MuonIsSignal->at(i);
        bool bad = MuonIsBad->at(i);
        float pt = MuonPt->at(i);
        float eta = MuonEta->at(i);
        float phi = MuonPhi->at(i);
        int q = MuonCharge->at(i);

        if (bad) {
            //std::cout << "Reject event because of bad muon!" << std::endl;
            return 1;
        }

        if (signal & OR) {
            //std::cout << "Reject event because of isolated muon!" << std::endl;
            return 1;

        }

        /*
        if (OR) {
            //std::cout << "Reject event because of baseline muon!" << std::endl;
            lv = true;
            return 1;
        }
        */

        MyMuon muon(pt, eta, phi);
        muon.SetIsSignal(signal);
        muon.SetPassOR(OR);
        muon.SetIsBad(bad);
        muon.SetCharge(q);
        recoMuons.push_back(muon);
    }
    GreaterByPt<MyMuon> ptComparator2_;
    std::sort(recoMuons.begin(), recoMuons.end(), ptComparator2_);
    double PtMu = 0;
    if (recoMuons.size() > 0) PtMu = recoMuons.at(0).Pt();

    //// Electrons ////////////////

    int NElectrons = ElePt->size();
    //std::cout << "NElectrons: " << NElectrons << std::endl;

    for (int i = 0; i < NElectrons; ++i) {
        bool OR = ElePassOR->at(i);
        bool signal = EleIsSignal->at(i);
        float pt = ElePt->at(i);
        float eta = EleEta->at(i);
        float phi = ElePhi->at(i);
        int q = EleCharge->at(i);

        if (signal & OR) {
            //std::cout << "Reject event because of isolated electon!" << std::endl;
            return 1;
        }

        /*
        if (OR) {
            //std::cout << "Reject event because of baseline electon!" << std::endl;
            lv = true;
            return 1;
        }
        */

        MyElectron electron(pt, eta, phi);
        electron.SetIsSignal(signal);
        electron.SetPassOR(OR);
        electron.SetCharge(q);
        recoElectrons.push_back(electron);
    }
    GreaterByPt<MyElectron> ptComparator3_;
    std::sort(recoElectrons.begin(), recoElectrons.end(), ptComparator3_);
    double PtEle = 0;
    if (recoElectrons.size() > 0) PtEle = recoElectrons.at(0).Pt();

    //// Taus ////////////////

    int NTaus = TauPt->size();
    //std::cout << "NTaus: " << NTaus << std::endl;

    for (int i = 0; i < NTaus; ++i) {
        bool OR = TauPassOR->at(i);
        bool signal = TauIsSignal->at(i);
        float pt = TauPt->at(i);
        float eta = TauEta->at(i);
        float phi = TauPhi->at(i);
        int q = TauCharge->at(i);
        MyTau tau(pt, eta, phi);
        tau.SetIsSignal(signal);
        tau.SetPassOR(OR);
        tau.SetCharge(q);
        recoTaus.push_back(tau);
    }
    GreaterByPt<MyTau> ptComparator4_;
    std::sort(recoTaus.begin(), recoTaus.end(), ptComparator4_);

    //// Photons ////////////////

    int NPhotons = PhotonPt->size();
    //std::cout << "NPhoton: " << NPhotons << std::endl;
    for (int i = 0; i < NPhotons; ++i) {
        bool OR = PhotonPassOR->at(i);
        bool signal = PhotonIsSignal->at(i);
        float pt = PhotonPt->at(i);
        float eta = PhotonEta->at(i);
        float phi = PhotonPhi->at(i);
        MyPhoton photon(pt, eta, phi);
        photon.SetIsSignal(signal);
        photon.SetPassOR(OR);
        recoPhotons.push_back(photon);
    }
    GreaterByPt<MyPhoton> ptComparator5_;
    std::sort(recoPhotons.begin(), recoPhotons.end(), ptComparator5_);
    double PtPho = 0;
    if (recoPhotons.size() > 0) PtPho = recoPhotons.at(0).Pt();

    //// Triggers ////////////////

    bool triggered = false;
    if (*xe90triggered || *xe110triggered) triggered = true;

    //// Let's start doing something ////////////////

    MyJet* firstJet = 0;
    MyJet* secondJet = 0;
    MyJet* thirdJet = 0;
    MyJet* fourthJet = 0;
    double HT = 0;
    TLorentzVector MHT(0., 0., 0., 0.);
    for ( auto& jet : recoJets) {
        if (!(jet.IsNoPU(m_jvtcut) || jet.Pt() > 60. || fabs(jet.Eta()) > 2.4 )) continue;
        HT += jet.Pt();
        MHT -= jet;
        if (firstJet == 0) {
            firstJet = &jet;
            //std::cout << "1: " << firstJet->Pt() << ", " << firstJet->Eta() << ", " << firstJet->Phi()<< std::endl;
        } else if (secondJet == 0) {
            secondJet = &jet;
            //std::cout << "2: " << secondJet->Pt() << ", " << secondJet->Eta() << ", " << secondJet->Phi()<< std::endl;
        } else if (thirdJet == 0) {
            thirdJet = &jet;
            //std::cout << "3: " << thirdJet->Pt() << ", " << thirdJet->Eta() << ", " << thirdJet->Phi()<< std::endl;
        } else if (fourthJet == 0) {
            fourthJet = &jet;
            //std::cout << "4: " << fourthJet->Pt() << ", " << fourthJet->Eta() << ", " << fourthJet->Phi()<< std::endl;
        } else {
            break;
        }
    }

    TLorentzVector MET;
    MET.SetPtEtaPhiM(*MET_pt, 0, *MET_phi, 0);
    TLorentzVector METsoft = MHT - MET;

    double dPhijj = 0.;
    if (secondJet) dPhijj = fabs(firstJet->DeltaPhi(*secondJet));
    double dEtajj = 0.;
    if (secondJet) dEtajj = fabs(firstJet->Eta() - secondJet->Eta());
    bool opphemijj = false;
    if (secondJet) opphemijj = (firstJet->Eta() * secondJet->Eta() < 0);
    bool soft3rd = false;
    if (thirdJet && !fourthJet) {
        if (thirdJet->Pt() > 25. && thirdJet->Pt() < 50.) {
            soft3rd = true;
        }
    }
    double Mjj = 0;
    if (secondJet) Mjj = (*firstJet + *secondJet).M();

    if (secondJet && dEtajj > 3.0) {
        if (firstJet->Pt() > 80. && secondJet->Pt() > 50.) {
            h_METsig->Fill(MET.Pt()/sqrt(HT), eventWeight);
            h_METsoft->Fill(METsoft.Pt(), eventWeight);
        }
    }

    if (secondJet && !thirdJet && !fourthJet) {
        //if (opphemijj && Mjj > 600.) {

        if (firstJet->Pt() < 80. || secondJet->Pt() < 50.) return 1;

        //if (thirdJet->Pt() > 50.) return 1;

        if (MET.Pt() > 100. && dEtajj > 0.0 && dPhijj < 2.7) {

            ++N_CR;

            int i_jet = 1;
            for ( auto& jet : recoJets) {
                std::cout << i_jet << "th reco jet (pT, eta, phi, btag, jvt): " << jet.Pt() << ", " << jet.Eta() << ", " <<  jet.Phi() << ", " <<  jet.IsB() << ", " << jet.IsPU(0.59) << std::endl;
                double dRmin = 999;
                MyJet* match = 0;
                for ( auto& genjet : genJets) {
                    double dR = jet.DeltaR(genjet);
                    if (dR < dRmin){
						dRmin = dR;
						match  = &genjet;
					}
                }
                if (dRmin < 0.15) {
                    if (jet.Pt() > 40. && jet.IsPU(0.59)) ++N_Good40tag;
                    if (jet.Pt() > 50. && jet.IsPU(0.59)) ++N_Good50tag;
                    if (jet.Pt() > 60. && jet.IsPU(0.59)) ++N_Good60tag;
                }
                if (dRmin > 0.15) {
                    if (jet.Pt() > 40. && jet.IsNoPU(0.59)) ++N_PU40notag;
                    if (jet.Pt() > 50. && jet.IsNoPU(0.59)) ++N_PU50notag;
                    if (jet.Pt() > 60. && jet.IsNoPU(0.59)) ++N_PU60notag;
                    if (jet.Pt() > 60. && jet.IsPU(0.59)) ++N_PU60tag;
                }
                ++i_jet;
            }

            i_jet = 1;
            for ( auto& genjet : genJets) {
                std::cout << i_jet << "th gen jet (pT, eta, phi, hf): " << genjet.Pt() << ", " << genjet.Eta() << ", " <<  genjet.Phi() << ", " <<  genjet.IsB() << std::endl;
                double dRmin = 999;
                MyJet* match = 0;
                for ( auto& jet : recoJets) {
                    double dR = jet.DeltaR(genjet);
                    if (dR < dRmin) {
						dRmin = dR;
						match  = &jet;
					}
                }
                if (dRmin > 0.15) {
                    if (genjet.Pt() > 50) ++N_Gen50Lost;
                }
                if (dRmin < 0.15 && i_jet < 3) {
                    if ( (match->Pt()/genjet.Pt()) < 0.5) ++N_LowResGen12;
                    if (match->Pt() > 60. && (match->Pt()/genjet.Pt()) > 2.) ++N_HighResReco60;
                }
                ++i_jet;
            }

            std::cout << "MET (pt, phi): " << *MET_pt << ", " << *MET_phi << std::endl;
            std::cout << "METmu (pt, phi): " << *METmu_pt << ", " << *METmu_phi << std::endl;
            std::cout << "Mjj: " << Mjj << std::endl;
            std::cout << "dEtajj: " << dEtajj << std::endl;

        }
        //}
    }

    if (secondJet && !thirdJet && opphemijj && Mjj > 0. && MET.Pt() > 0.) {

        if (firstJet->Pt() < 80. || secondJet->Pt() < 50.) return 1;

        if (MET.Pt() > 150. && !triggered) return 1;
        if (dEtajj > 4.8) h_MET_vs_dPhi_Incl->Fill(MET.Pt(), dPhijj, eventWeight);
        if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) h_MET_vs_dEta_Incl->Fill(MET.Pt(), dEtajj, eventWeight);
        if (Mjj > 600. && Mjj < 1000.) {
            if (dEtajj > 4.8) h_MET_vs_dPhi_Mjj600->Fill(MET.Pt(), dPhijj, eventWeight);
            if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) h_MET_vs_dEta_Mjj600->Fill(MET.Pt(), dEtajj, eventWeight);
        }
        if (Mjj > 1000. && Mjj < 1500.) {
            if (dEtajj > 4.8) h_MET_vs_dPhi_Mjj1000->Fill(MET.Pt(), dPhijj, eventWeight);
            if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) h_MET_vs_dEta_Mjj1000->Fill(MET.Pt(), dEtajj, eventWeight);
        }
        if (Mjj > 1500. && Mjj < 2000.) {
            if (dEtajj > 4.8) h_MET_vs_dPhi_Mjj1500->Fill(MET.Pt(), dPhijj, eventWeight);
            if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) h_MET_vs_dEta_Mjj1500->Fill(MET.Pt(), dEtajj, eventWeight);
        }
        if (Mjj > 2000. ) {
            if (dEtajj > 4.8) h_MET_vs_dPhi_Mjj2000->Fill(MET.Pt(), dPhijj, eventWeight);
            if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) h_MET_vs_dEta_Mjj2000->Fill(MET.Pt(), dEtajj, eventWeight);
        }

        if (dEtajj > 0) {

            if (firstJet) {
                double dphi = fabs(MET.DeltaPhi(*firstJet));
                h_Jet1_Pt->Fill(firstJet->Pt(),eventWeight);
                h_Jet1_Eta->Fill(firstJet->Eta(),eventWeight);
                h_Jet1_Phi->Fill(firstJet->Phi(),eventWeight);
                h_Jet1_DeltaPhi->Fill(dphi ,eventWeight);
            }

            if (secondJet) {
                double dphi = fabs(MET.DeltaPhi(*secondJet));
                h_Jet2_Pt->Fill(secondJet->Pt(),eventWeight);
                h_Jet2_Eta->Fill(secondJet->Eta(),eventWeight);
                h_Jet2_Phi->Fill(secondJet->Phi(),eventWeight);
                h_Jet2_DeltaPhi->Fill(dphi ,eventWeight);
            }

            if (thirdJet) {
                double dphi = fabs(MET.DeltaPhi(*thirdJet));
                h_Jet3_Pt->Fill(thirdJet->Pt(),eventWeight);
                h_Jet3_Eta->Fill(thirdJet->Eta(),eventWeight);
                h_Jet3_Phi->Fill(thirdJet->Phi(),eventWeight);
                h_Jet3_DeltaPhi->Fill(dphi ,eventWeight);
            }

            if (dEtajj > 0.0 && Mjj > 0.) h_MET->Fill(MET.Pt(), eventWeight);
            if (dEtajj > 0.0 && Mjj > 0.) h_DeltaPhijj->Fill(dPhijj, eventWeight);
            if (dEtajj > 0.0 && Mjj > 0.) h_DeltaEtajj->Fill(dEtajj, eventWeight);
            if (dEtajj > 0.0 && Mjj > 0.) h_Mjj->Fill(Mjj, eventWeight);

            if (MET.Pt() > METMin_CR && MET.Pt() < METMax_CR) {
                if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                    N_METCR_dPhiCR += eventWeight;
                    w2_METCR_dPhiCR += eventWeight*eventWeight;
                }
                if (dPhijj < dPhijjMax_SR) {
                    N_METCR_dPhiSR += eventWeight;
                    w2_METCR_dPhiSR += eventWeight*eventWeight;
                }
            }
            if (MET.Pt() > METMin_SR) {
                if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                    N_METSR_dPhiCR += eventWeight*eventWeight;
                    w2_METSR_dPhiCR += eventWeight;
                }
                if (dPhijj < dPhijjMax_SR) {
                    N_METSR_dPhiSR += eventWeight;
                    w2_METSR_dPhiSR += eventWeight*eventWeight;
                }
            }

            if (Mjj > 600. && Mjj < 1000.) {
                if (MET.Pt() > METMin_CR && MET.Pt() < METMax_CR) {
                    if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                        N_METCR_dPhiCR_Mjj600 += eventWeight;
                        w2_METCR_dPhiCR_Mjj600 += eventWeight*eventWeight;
                    }
                    if (dPhijj < dPhijjMax_SR) {
                        N_METCR_dPhiSR_Mjj600 += eventWeight;
                        w2_METCR_dPhiSR_Mjj600 += eventWeight*eventWeight;
                    }
                }
                if (MET.Pt() > METMin_SR) {
                    if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                        N_METSR_dPhiCR_Mjj600 += eventWeight*eventWeight;
                        w2_METSR_dPhiCR_Mjj600 += eventWeight;
                    }
                    if (dPhijj < dPhijjMax_SR) {
                        N_METSR_dPhiSR_Mjj600 += eventWeight;
                        w2_METSR_dPhiSR_Mjj600 += eventWeight*eventWeight;
                    }
                }
            }

            if (Mjj > 1000. && Mjj < 1500.) {
                if (MET.Pt() > METMin_CR && MET.Pt() < METMax_CR) {
                    if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                        N_METCR_dPhiCR_Mjj1000 += eventWeight;
                        w2_METCR_dPhiCR_Mjj1000 += eventWeight*eventWeight;
                    }
                    if (dPhijj < dPhijjMax_SR) {
                        N_METCR_dPhiSR_Mjj1000 += eventWeight;
                        w2_METCR_dPhiSR_Mjj1000 += eventWeight*eventWeight;
                    }
                }
                if (MET.Pt() > METMin_SR) {
                    if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                        N_METSR_dPhiCR_Mjj1000 += eventWeight*eventWeight;
                        w2_METSR_dPhiCR_Mjj1000 += eventWeight;
                    }
                    if (dPhijj < dPhijjMax_SR) {
                        N_METSR_dPhiSR_Mjj1000 += eventWeight;
                        w2_METSR_dPhiSR_Mjj1000 += eventWeight*eventWeight;
                    }
                }
            }

            if (Mjj > 1500. && Mjj < 2000.) {
                if (MET.Pt() > METMin_CR && MET.Pt() < METMax_CR) {
                    if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                        N_METCR_dPhiCR_Mjj1500 += eventWeight;
                        w2_METCR_dPhiCR_Mjj1500 += eventWeight*eventWeight;
                    }
                    if (dPhijj < dPhijjMax_SR) {
                        N_METCR_dPhiSR_Mjj1500 += eventWeight;
                        w2_METCR_dPhiSR_Mjj1500 += eventWeight*eventWeight;
                    }
                }
                if (MET.Pt() > METMin_SR) {
                    if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                        N_METSR_dPhiCR_Mjj1500 += eventWeight*eventWeight;
                        w2_METSR_dPhiCR_Mjj1500 += eventWeight;
                    }
                    if (dPhijj < dPhijjMax_SR) {
                        N_METSR_dPhiSR_Mjj1500 += eventWeight;
                        w2_METSR_dPhiSR_Mjj1500 += eventWeight*eventWeight;
                    }
                }
            }

            if (Mjj > 2000. ) {
                if (MET.Pt() > METMin_CR && MET.Pt() < METMax_CR) {
                    if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                        N_METCR_dPhiCR_Mjj2000 += eventWeight;
                        w2_METCR_dPhiCR_Mjj2000 += eventWeight*eventWeight;
                    }
                    if (dPhijj < dPhijjMax_SR) {
                        N_METCR_dPhiSR_Mjj2000 += eventWeight;
                        w2_METCR_dPhiSR_Mjj2000 += eventWeight*eventWeight;
                    }
                }
                if (MET.Pt() > METMin_SR) {
                    if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                        N_METSR_dPhiCR_Mjj2000 += eventWeight*eventWeight;
                        w2_METSR_dPhiCR_Mjj2000 += eventWeight;
                    }
                    if (dPhijj < dPhijjMax_SR) {
                        N_METSR_dPhiSR_Mjj2000 += eventWeight;
                        w2_METSR_dPhiSR_Mjj2000 += eventWeight*eventWeight;
                    }
                }
            }

        } else { // if dEtajj<4.8

            if (MET.Pt() > METMin_CR && MET.Pt() < METMax_CR) {
                if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                    N_METCR_dPhiCR_dEtaCR += eventWeight;
                    w2_METCR_dPhiCR_dEtaCR += eventWeight*eventWeight;
                }
                if (dPhijj < dPhijjMax_SR) {
                    N_METCR_dPhiSR_dEtaCR += eventWeight;
                    w2_METCR_dPhiSR_dEtaCR += eventWeight*eventWeight;
                }
            }
            if (MET.Pt() > METMin_SR) {
                if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                    N_METSR_dPhiCR_dEtaCR += eventWeight*eventWeight;
                    w2_METSR_dPhiCR_dEtaCR += eventWeight;
                }
                if (dPhijj < dPhijjMax_SR) {
                    N_METSR_dPhiSR_dEtaCR += eventWeight;
                    w2_METSR_dPhiSR_dEtaCR += eventWeight*eventWeight;
                }
            }

            if (Mjj > 600. && Mjj < 1000.) {
                if (MET.Pt() > METMin_CR && MET.Pt() < METMax_CR) {
                    if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                        N_METCR_dPhiCR_dEtaCR_Mjj600 += eventWeight;
                        w2_METCR_dPhiCR_dEtaCR_Mjj600 += eventWeight*eventWeight;
                    }
                    if (dPhijj < dPhijjMax_SR) {
                        N_METCR_dPhiSR_dEtaCR_Mjj600 += eventWeight;
                        w2_METCR_dPhiSR_dEtaCR_Mjj600 += eventWeight*eventWeight;
                    }
                }
                if (MET.Pt() > METMin_SR) {
                    if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                        N_METSR_dPhiCR_dEtaCR_Mjj600 += eventWeight*eventWeight;
                        w2_METSR_dPhiCR_dEtaCR_Mjj600 += eventWeight;
                    }
                    if (dPhijj < dPhijjMax_SR) {
                        N_METSR_dPhiSR_dEtaCR_Mjj600 += eventWeight;
                        w2_METSR_dPhiSR_dEtaCR_Mjj600 += eventWeight*eventWeight;
                    }
                }
            }

            if (Mjj > 1000. && Mjj < 1500.) {
                if (MET.Pt() > METMin_CR && MET.Pt() < METMax_CR) {
                    if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                        N_METCR_dPhiCR_dEtaCR_Mjj1000 += eventWeight;
                        w2_METCR_dPhiCR_dEtaCR_Mjj1000 += eventWeight*eventWeight;
                    }
                    if (dPhijj < dPhijjMax_SR) {
                        N_METCR_dPhiSR_dEtaCR_Mjj1000 += eventWeight;
                        w2_METCR_dPhiSR_dEtaCR_Mjj1000 += eventWeight*eventWeight;
                    }
                }
                if (MET.Pt() > METMin_SR) {
                    if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                        N_METSR_dPhiCR_dEtaCR_Mjj1000 += eventWeight*eventWeight;
                        w2_METSR_dPhiCR_dEtaCR_Mjj1000 += eventWeight;
                    }
                    if (dPhijj < dPhijjMax_SR) {
                        N_METSR_dPhiSR_dEtaCR_Mjj1000 += eventWeight;
                        w2_METSR_dPhiSR_dEtaCR_Mjj1000 += eventWeight*eventWeight;
                    }
                }
            }

            if (Mjj > 1500. && Mjj < 2000.) {
                if (MET.Pt() > METMin_CR && MET.Pt() < METMax_CR) {
                    if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                        N_METCR_dPhiCR_dEtaCR_Mjj1500 += eventWeight;
                        w2_METCR_dPhiCR_dEtaCR_Mjj1500 += eventWeight*eventWeight;
                    }
                    if (dPhijj < dPhijjMax_SR) {
                        N_METCR_dPhiSR_dEtaCR_Mjj1500 += eventWeight;
                        w2_METCR_dPhiSR_dEtaCR_Mjj1500 += eventWeight*eventWeight;
                    }
                }
                if (MET.Pt() > METMin_SR) {
                    if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                        N_METSR_dPhiCR_dEtaCR_Mjj1500 += eventWeight*eventWeight;
                        w2_METSR_dPhiCR_dEtaCR_Mjj1500 += eventWeight;
                    }
                    if (dPhijj < dPhijjMax_SR) {
                        N_METSR_dPhiSR_dEtaCR_Mjj1500 += eventWeight;
                        w2_METSR_dPhiSR_dEtaCR_Mjj1500 += eventWeight*eventWeight;
                    }
                }
            }

            if (Mjj > 2000. ) {
                if (MET.Pt() > METMin_CR && MET.Pt() < METMax_CR) {
                    if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                        N_METCR_dPhiCR_dEtaCR_Mjj2000 += eventWeight;
                        w2_METCR_dPhiCR_dEtaCR_Mjj2000 += eventWeight*eventWeight;
                    }
                    if (dPhijj < dPhijjMax_SR) {
                        N_METCR_dPhiSR_dEtaCR_Mjj2000 += eventWeight;
                        w2_METCR_dPhiSR_dEtaCR_Mjj2000 += eventWeight*eventWeight;
                    }
                }
                if (MET.Pt() > METMin_SR) {
                    if (dPhijj > dPhijjMin_CR && dPhijj < dPhijjMax_CR) {
                        N_METSR_dPhiCR_dEtaCR_Mjj2000 += eventWeight*eventWeight;
                        w2_METSR_dPhiCR_dEtaCR_Mjj2000 += eventWeight;
                    }
                    if (dPhijj < dPhijjMax_SR) {
                        N_METSR_dPhiSR_dEtaCR_Mjj2000 += eventWeight;
                        w2_METSR_dPhiSR_dEtaCR_Mjj2000 += eventWeight*eventWeight;
                    }
                }
            }

        }

    }

    return kTRUE;

}

void MyABCDStudies::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.

    std::cout << "N_METCR_dPhiCR = " << N_METCR_dPhiCR << "+-" << sqrt(w2_METCR_dPhiCR) << std::endl;
    std::cout << "N_METCR_dPhiSR = " << N_METCR_dPhiSR << "+-" << sqrt(w2_METCR_dPhiSR) << std::endl;
    std::cout << "N_METSR_dPhiCR = " << N_METSR_dPhiCR << "+-" << sqrt(w2_METSR_dPhiCR) << std::endl;
    std::cout << "N_METSR_dPhiSR = " << N_METSR_dPhiSR << "+-" << sqrt(w2_METSR_dPhiSR) << std::endl;

    std::cout << std::endl;

    std::cout << "N_METCR_dPhiCR_Mjj600 = " << N_METCR_dPhiCR_Mjj600 << "+-" << sqrt(w2_METCR_dPhiCR_Mjj600) << std::endl;
    std::cout << "N_METCR_dPhiSR_Mjj600 = " << N_METCR_dPhiSR_Mjj600 << "+-" << sqrt(w2_METCR_dPhiSR_Mjj600) << std::endl;
    std::cout << "N_METSR_dPhiCR_Mjj600 = " << N_METSR_dPhiCR_Mjj600 << "+-" << sqrt(w2_METSR_dPhiCR_Mjj600) << std::endl;
    std::cout << "N_METSR_dPhiSR_Mjj600 = " << N_METSR_dPhiSR_Mjj600 << "+-" << sqrt(w2_METSR_dPhiSR_Mjj600) << std::endl;

    std::cout << std::endl;

    std::cout << "N_METCR_dPhiCR_Mjj1000 = " << N_METCR_dPhiCR_Mjj1000 << "+-" << sqrt(w2_METCR_dPhiCR_Mjj1000) << std::endl;
    std::cout << "N_METCR_dPhiSR_Mjj1000 = " << N_METCR_dPhiSR_Mjj1000 << "+-" << sqrt(w2_METCR_dPhiSR_Mjj1000) << std::endl;
    std::cout << "N_METSR_dPhiCR_Mjj1000 = " << N_METSR_dPhiCR_Mjj1000 << "+-" << sqrt(w2_METSR_dPhiCR_Mjj1000) << std::endl;
    std::cout << "N_METSR_dPhiSR_Mjj1000 = " << N_METSR_dPhiSR_Mjj1000 << "+-" << sqrt(w2_METSR_dPhiSR_Mjj1000) << std::endl;

    std::cout << std::endl;

    std::cout << "N_METCR_dPhiCR_Mjj1500 = " << N_METCR_dPhiCR_Mjj1500 << "+-" << sqrt(w2_METCR_dPhiCR_Mjj1500) << std::endl;
    std::cout << "N_METCR_dPhiSR_Mjj1500 = " << N_METCR_dPhiSR_Mjj1500 << "+-" << sqrt(w2_METCR_dPhiSR_Mjj1500) << std::endl;
    std::cout << "N_METSR_dPhiCR_Mjj1500 = " << N_METSR_dPhiCR_Mjj1500 << "+-" << sqrt(w2_METSR_dPhiCR_Mjj1500) << std::endl;
    std::cout << "N_METSR_dPhiSR_Mjj1500 = " << N_METSR_dPhiSR_Mjj1500 << "+-" << sqrt(w2_METSR_dPhiSR_Mjj1500) << std::endl;

    std::cout << std::endl;

    std::cout << "N_METCR_dPhiCR_Mjj2000 = " << N_METCR_dPhiCR_Mjj2000 << "+-" << sqrt(w2_METCR_dPhiCR_Mjj2000) << std::endl;
    std::cout << "N_METCR_dPhiSR_Mjj2000 = " << N_METCR_dPhiSR_Mjj2000 << "+-" << sqrt(w2_METCR_dPhiSR_Mjj2000) << std::endl;
    std::cout << "N_METSR_dPhiCR_Mjj2000 = " << N_METSR_dPhiCR_Mjj2000 << "+-" << sqrt(w2_METSR_dPhiCR_Mjj2000) << std::endl;
    std::cout << "N_METSR_dPhiSR_Mjj2000 = " << N_METSR_dPhiSR_Mjj2000 << "+-" << sqrt(w2_METSR_dPhiSR_Mjj2000) << std::endl;

    std::cout << std::endl;

    std::cout << "N_METCR_dPhiCR_dEtaCR = " << N_METCR_dPhiCR_dEtaCR << "+-" << sqrt(w2_METCR_dPhiCR_dEtaCR) << std::endl;
    std::cout << "N_METCR_dPhiSR_dEtaCR = " << N_METCR_dPhiSR_dEtaCR << "+-" << sqrt(w2_METCR_dPhiSR_dEtaCR) << std::endl;
    std::cout << "N_METSR_dPhiCR_dEtaCR = " << N_METSR_dPhiCR_dEtaCR << "+-" << sqrt(w2_METSR_dPhiCR_dEtaCR) << std::endl;
    std::cout << "N_METSR_dPhiSR_dEtaCR = " << N_METSR_dPhiSR_dEtaCR << "+-" << sqrt(w2_METSR_dPhiSR_dEtaCR) << std::endl;

    std::cout << std::endl;

    std::cout << "N_METCR_dPhiCR_dEtaCR_Mjj600 = " << N_METCR_dPhiCR_dEtaCR_Mjj600 << "+-" << sqrt(w2_METCR_dPhiCR_dEtaCR_Mjj600) << std::endl;
    std::cout << "N_METCR_dPhiSR_dEtaCR_Mjj600 = " << N_METCR_dPhiSR_dEtaCR_Mjj600 << "+-" << sqrt(w2_METCR_dPhiSR_dEtaCR_Mjj600) << std::endl;
    std::cout << "N_METSR_dPhiCR_dEtaCR_Mjj600 = " << N_METSR_dPhiCR_dEtaCR_Mjj600 << "+-" << sqrt(w2_METSR_dPhiCR_dEtaCR_Mjj600) << std::endl;
    std::cout << "N_METSR_dPhiSR_dEtaCR_Mjj600 = " << N_METSR_dPhiSR_dEtaCR_Mjj600 << "+-" << sqrt(w2_METSR_dPhiSR_dEtaCR_Mjj600) << std::endl;

    std::cout << std::endl;

    std::cout << "N_METCR_dPhiCR_dEtaCR_Mjj1000 = " << N_METCR_dPhiCR_dEtaCR_Mjj1000 << "+-" << sqrt(w2_METCR_dPhiCR_dEtaCR_Mjj1000) << std::endl;
    std::cout << "N_METCR_dPhiSR_dEtaCR_Mjj1000 = " << N_METCR_dPhiSR_dEtaCR_Mjj1000 << "+-" << sqrt(w2_METCR_dPhiSR_dEtaCR_Mjj1000) << std::endl;
    std::cout << "N_METSR_dPhiCR_dEtaCR_Mjj1000 = " << N_METSR_dPhiCR_dEtaCR_Mjj1000 << "+-" << sqrt(w2_METSR_dPhiCR_dEtaCR_Mjj1000) << std::endl;
    std::cout << "N_METSR_dPhiSR_dEtaCR_Mjj1000 = " << N_METSR_dPhiSR_dEtaCR_Mjj1000 << "+-" << sqrt(w2_METSR_dPhiSR_dEtaCR_Mjj1000) << std::endl;

    std::cout << std::endl;

    std::cout << "N_METCR_dPhiCR_dEtaCR_Mjj1500 = " << N_METCR_dPhiCR_dEtaCR_Mjj1500 << "+-" << sqrt(w2_METCR_dPhiCR_dEtaCR_Mjj1500) << std::endl;
    std::cout << "N_METCR_dPhiSR_dEtaCR_Mjj1500 = " << N_METCR_dPhiSR_dEtaCR_Mjj1500 << "+-" << sqrt(w2_METCR_dPhiSR_dEtaCR_Mjj1500) << std::endl;
    std::cout << "N_METSR_dPhiCR_dEtaCR_Mjj1500 = " << N_METSR_dPhiCR_dEtaCR_Mjj1500 << "+-" << sqrt(w2_METSR_dPhiCR_dEtaCR_Mjj1500) << std::endl;
    std::cout << "N_METSR_dPhiSR_dEtaCR_Mjj1500 = " << N_METSR_dPhiSR_dEtaCR_Mjj1500 << "+-" << sqrt(w2_METSR_dPhiSR_dEtaCR_Mjj1500) << std::endl;

    std::cout << std::endl;

    std::cout << "N_METCR_dPhiCR_dEtaCR_Mjj2000 = " << N_METCR_dPhiCR_dEtaCR_Mjj2000 << "+-" << sqrt(w2_METCR_dPhiCR_dEtaCR_Mjj2000) << std::endl;
    std::cout << "N_METCR_dPhiSR_dEtaCR_Mjj2000 = " << N_METCR_dPhiSR_dEtaCR_Mjj2000 << "+-" << sqrt(w2_METCR_dPhiSR_dEtaCR_Mjj2000) << std::endl;
    std::cout << "N_METSR_dPhiCR_dEtaCR_Mjj2000 = " << N_METSR_dPhiCR_dEtaCR_Mjj2000 << "+-" << sqrt(w2_METSR_dPhiCR_dEtaCR_Mjj2000) << std::endl;
    std::cout << "N_METSR_dPhiSR_dEtaCR_Mjj2000 = " << N_METSR_dPhiSR_dEtaCR_Mjj2000 << "+-" << sqrt(w2_METSR_dPhiSR_dEtaCR_Mjj2000) << std::endl;

    std::cout << std::endl;

    std::cout << "N_CR: " << N_CR << std::endl;
    std::cout << "N_Gen50Lost: " << N_Gen50Lost << std::endl;
    std::cout << "N_LowResGen12: " << N_LowResGen12 << std::endl;
    std::cout << "N_HighResReco60: " << N_HighResReco60 << std::endl;
    std::cout << "N_PU60tag: " << N_PU60tag << std::endl;
    std::cout << "N_PU60notag: " << N_PU60notag << std::endl;
    std::cout << "N_PU50notag: " << N_PU50notag << std::endl;
    std::cout << "N_PU40notag: " << N_PU40notag << std::endl;
    std::cout << "N_Good60tag: " << N_Good60tag << std::endl;
    std::cout << "N_Good50tag: " << N_Good50tag << std::endl;
    std::cout << "N_Good40tag: " << N_Good40tag << std::endl;


    for (unsigned int i = 0; i < histos_1D.size(); ++i) {
        histos_1D.at(i)->Write();
    }

    for (unsigned int i = 0; i < histos_2D.size(); ++i) {
        histos_2D.at(i)->Write();
    }

    outputfile->Close();

    for (std::map<UInt_t, UInt_t>::iterator it=ProcessedEvents.begin(); it != ProcessedEvents.end(); ++it) {
        std::cout << it->first << " " << it->second << std::endl;
    }

}

void MyABCDStudies::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

}

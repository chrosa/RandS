#define MyMETStudies_cxx
// The class definition in MyMETStudies.h has been generated automatically
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
// root> T->Process("MyMETStudies.C")
// root> T->Process("MyMETStudies.C","some options")
// root> T->Process("MyMETStudies.C+")
//

#include "MyJet.h"
#include "MyElectron.h"
#include "MyPhoton.h"
#include "MyMuon.h"
#include "MyTau.h"
#include "MyMETStudies.h"

void MyMETStudies::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

    outputfile = new TFile("METStudiesOutput_data.root","RECREATE");

    m_MatchingCone = 0.1;
    m_jvtcut = 0.59;
    m_lumi = 36100.;

    //// Book 1d histograms

    h_Jet1_Pt = new TH1F("h_Jet1_Pt", "h_Jet1_Pt", 100, 0., 1000.);
    h_Jet1_Pt->Sumw2();
    histos_1D.push_back(h_Jet1_Pt);

    h_Jet2_Pt = new TH1F("h_Jet2_Pt", "h_Jet2_Pt", 100, 0., 1000.);
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

    h_MHT = new TH1F("h_MHT", "h_MHT", 100, 0., 1000.);
    h_MHT->Sumw2();
    histos_1D.push_back(h_MHT);

    h_MHTnoOR = new TH1F("h_MHTnoOR", "h_MHTnoOR", 100, 0., 1000.);
    h_MHTnoOR->Sumw2();
    histos_1D.push_back(h_MHTnoOR);

    h_MHTnoJVT = new TH1F("h_MHTnoJVT", "h_MHTnoJVT", 100, 0., 1000.);
    h_MHTnoJVT->Sumw2();
    histos_1D.push_back(h_MHTnoJVT);

    h_MHTnoJVTnoOR = new TH1F("h_MHTnoJVTnoOR", "h_MHTnoJVTnoOR", 100, 0., 1000.);
    h_MHTnoJVTnoOR->Sumw2();
    histos_1D.push_back(h_MHTnoJVTnoOR);

    h_TruthMHT = new TH1F("h_TruthMHT", "h_TruthMHT", 100, 0., 1000.);
    h_TruthMHT->Sumw2();
    histos_1D.push_back(h_TruthMHT);

    h_MET = new TH1F("h_MET", "h_MET", 100, 0., 1000.);
    h_MET->Sumw2();
    histos_1D.push_back(h_MET);

    h_TruthMET = new TH1F("h_TruthMET", "h_TruthMET", 100, 0., 1000.);
    h_TruthMET->Sumw2();
    histos_1D.push_back(h_TruthMET);

    h_METreplaced = new TH1F("h_METreplaced", "h_METreplaced", 100, 0., 1000.);
    h_METreplaced->Sumw2();
    histos_1D.push_back(h_METreplaced);

    //// Book 2d histograms

    h_MET_vs_MHT = new TH2F("h_MET_vs_MHT", "h_MET_vs_MHT", 100, 0., 1000., 100, 0., 1000.);
    h_MET_vs_MHT->Sumw2();
    histos_2D.push_back(h_MET_vs_MHT);

    h_MET_vs_MHTnoOR = new TH2F("h_MET_vs_MHTnoOR", "h_MET_vs_MHTnoOR", 100, 0., 1000., 100, 0., 1000.);
    h_MET_vs_MHTnoOR->Sumw2();
    histos_2D.push_back(h_MET_vs_MHTnoOR);

    h_MET_vs_MHTnoJVT = new TH2F("h_MET_vs_MHTnoJVT", "h_MET_vs_MHTnoJVT", 100, 0., 1000., 100, 0., 1000.);
    h_MET_vs_MHTnoJVT->Sumw2();
    histos_2D.push_back(h_MET_vs_MHTnoJVT);

    h_MET_vs_MHTnoJVTnoOR = new TH2F("h_MET_vs_MHTnoJVTnoOR", "h_MET_vs_MHTnoJVTnoOR", 100, 0., 1000., 100, 0., 1000.);
    h_MET_vs_MHTnoJVTnoOR->Sumw2();
    histos_2D.push_back(h_MET_vs_MHTnoJVTnoOR);

    h_MHT_vs_MHTnoOR = new TH2F("h_MHT_vs_MHTnoOR", "h_MHT_vs_MHTnoOR", 100, 0., 1000., 100, 0., 1000.);
    h_MHT_vs_MHTnoOR->Sumw2();
    histos_2D.push_back(h_MHT_vs_MHTnoOR);

    h_MHT_vs_MHTnoJVT = new TH2F("h_MHT_vs_MHTnoJVT", "h_MHT_vs_MHTnoJVT", 100, 0., 1000., 100, 0., 1000.);
    h_MHT_vs_MHTnoJVT->Sumw2();
    histos_2D.push_back(h_MHT_vs_MHTnoJVT);

    h_MHT_vs_MHTnoJVTnoOR = new TH2F("h_MHT_vs_MHTnoJVTnoOR", "h_MHT_vs_MHTnoJVTnoOR", 100, 0., 1000., 100, 0., 1000.);
    h_MHT_vs_MHTnoJVTnoOR->Sumw2();
    histos_2D.push_back(h_MHT_vs_MHTnoJVTnoOR);

    h_JVTPhi_vs_METPhi = new TH2F("h_JVTPhi_vs_METPhi", "h_JVTPhi_vs_METPhi", 40, -TMath::Pi(), TMath::Pi(), 40, -TMath::Pi(), TMath::Pi());
    h_JVTPhi_vs_METPhi->Sumw2();
    histos_2D.push_back(h_JVTPhi_vs_METPhi);

    h_JVTPhi_vs_MHTPhi = new TH2F("h_JVTPhi_vs_MHTPhi", "h_JVTPhi_vs_MHTPhi", 40, -TMath::Pi(), TMath::Pi(), 40, -TMath::Pi(), TMath::Pi());
    h_JVTPhi_vs_MHTPhi->Sumw2();
    histos_2D.push_back(h_JVTPhi_vs_MHTPhi);

    h_JVTPhi_vs_MHTnoORPhi = new TH2F("h_JVTPhi_vs_MHTnoORPhi", "h_JVTPhi_vs_MHTnoORPhi", 40, -TMath::Pi(), TMath::Pi(), 40, -TMath::Pi(), TMath::Pi());
    h_JVTPhi_vs_MHTnoORPhi->Sumw2();
    histos_2D.push_back(h_JVTPhi_vs_MHTnoORPhi);

    h_MET_vs_METreplaced = new TH2F("h_MET_vs_METreplaced", "h_MET_vs_METreplaced", 100, 0., 1000., 100, 0., 1000.);
    h_MET_vs_METreplaced->Sumw2();
    histos_2D.push_back(h_MET_vs_METreplaced);

    h_MHTnoOR_vs_MHTnoORreplaced = new TH2F("h_MHTnoOR_vs_MHTnoORreplaced", "h_MHTnoOR_vs_MHTnoORreplaced", 100, 0., 1000., 100, 0., 1000.);
    h_MHTnoOR_vs_MHTnoORreplaced->Sumw2();
    histos_2D.push_back(h_MHTnoOR_vs_MHTnoORreplaced);

    h_MHTnoJVTnoOR_vs_MHTnoJVTnoORreplaced = new TH2F("h_MHTnoJVTnoOR_vs_MHTnoJVTnoORreplaced", "h_MHTnoJVTnoOR_vs_MHTnoJVTnoORreplaced", 100, 0., 1000., 100, 0., 1000.);
    h_MHTnoJVTnoOR_vs_MHTnoJVTnoORreplaced->Sumw2();
    histos_2D.push_back(h_MHTnoJVTnoOR_vs_MHTnoJVTnoORreplaced);

    h_MET_vs_ElePt = new TH2F("h_MET_vs_ElePt", "h_MET_vs_ElePt", 100, 0., 1000., 100, 0., 1000.);
    h_MET_vs_ElePt->Sumw2();
    histos_2D.push_back(h_MET_vs_ElePt);

    h_MET_vs_MuPt = new TH2F("h_MET_vs_MuPt", "h_MET_vs_MuPt", 100, 0., 1000., 100, 0., 1000.);
    h_MET_vs_MuPt->Sumw2();
    histos_2D.push_back(h_MET_vs_MuPt);

    h_MET_vs_PhoPt = new TH2F("h_MET_vs_PhoPt", "h_MET_vs_PhoPt", 100, 0., 1000., 100, 0., 1000.);
    h_MET_vs_PhoPt->Sumw2();
    histos_2D.push_back(h_MET_vs_PhoPt);

    h_MET_vs_LepPt = new TH2F("h_MET_vs_LepPt", "h_MET_vs_LepPt", 100, 0., 1000., 100, 0., 1000.);
    h_MET_vs_LepPt->Sumw2();
    histos_2D.push_back(h_MET_vs_LepPt);

    h_MET_vs_METjet = new TH2F("h_MET_vs_METjet", "h_MET_vs_METjet", 100, 0., 1000., 100, 0., 1000.);
    h_MET_vs_METjet->Sumw2();
    histos_2D.push_back(h_MET_vs_METjet);

    h_MET_vs_METele = new TH2F("h_MET_vs_METele", "h_MET_vs_METele", 100, 0., 1000., 100, 0., 1000.);
    h_MET_vs_METele->Sumw2();
    histos_2D.push_back(h_MET_vs_METele);

    h_MET_vs_METmu = new TH2F("h_MET_vs_METmu", "h_MET_vs_METmu", 100, 0., 1000., 100, 0., 1000.);
    h_MET_vs_METmu->Sumw2();
    histos_2D.push_back(h_MET_vs_METmu);

    h_MET_vs_METgamma = new TH2F("h_MET_vs_METgamma", "h_MET_vs_METgamma", 100, 0., 1000., 100, 0., 1000.);
    h_MET_vs_METgamma->Sumw2();
    histos_2D.push_back(h_MET_vs_METgamma);

    h_MET_vs_METtrack = new TH2F("h_MET_vs_METtrack", "h_MET_vs_METtrack", 100, 0., 1000., 100, 0., 1000.);
    h_MET_vs_METtrack->Sumw2();
    histos_2D.push_back(h_MET_vs_METtrack);

    //NTotEvents = fChain->GetEntries();

    //// Not very elegant! TODO: Store this info in and read from file

    // [v3]
    AvailableEvents[361022] = 1993647;
    AvailableEvents[361023] = 7884494;
    AvailableEvents[361024] = 7889800;
    AvailableEvents[361025] = 7977600;
    AvailableEvents[361026] = 1893400;

}

void MyMETStudies::SlaveBegin(TTree * /*tree*/)
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

}

Bool_t MyMETStudies::Process(Long64_t entry)
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

    //std::cout << entry << std::endl;
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
    //std::cout << "Weight: " << *Weight << std::endl;
    double eventWeight = *Weight;
    if (isMC) eventWeight *= m_lumi / AvailableEvents[*DatasetID];

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

        recoJets.push_back(jet);
        //if ( jet.IsNoPU(m_jvtcut) || jet.Pt() > 60. || fabs(jet.Eta()) > 2.4 ) recoJets.push_back(jet);

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

    bool lv = false;

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

        if (signal) {
            //std::cout << "Reject event because of isolated muon!" << std::endl;
            lv = true;
            if (!cutFlowStudies) return 1;

        }

        if (OR) {
            //std::cout << "Reject event because of baseline muon!" << std::endl;
            lv = true;
            if (!cutFlowStudies) return 1;
        }

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

        if (signal) {
            //std::cout << "Reject event because of isolated electon!" << std::endl;
            lv = true;
            if (!cutFlowStudies) return 1;
        }

        if (OR) {
            //std::cout << "Reject event because of baseline electon!" << std::endl;
            lv = true;
            if (!cutFlowStudies) return 1;
        }

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

    Ntot_VR1 += 1;
    Ntot_w_VR1 += eventWeight;
    Ntot_VR2 += 1;
    Ntot_w_VR2 += eventWeight;
    Ntot_VR3 += 1;
    Ntot_w_VR3 += eventWeight;
    if (triggered) {
        Ntot_VR1_trig += 1;
        Ntot_w_VR1_trig += eventWeight;
        Ntot_VR2_trig += 1;
        Ntot_w_VR2_trig += eventWeight;
        Ntot_VR3_trig += 1;
        Ntot_w_VR3_trig += eventWeight;
    }


    MyJet* firstJet = 0;
    MyJet* secondJet = 0;
    MyJet* thirdJet = 0;
    MyJet* fourthJet = 0;

    for ( auto& jet : recoJets) {
        if (jet.IsGood() && jet.Pt() > 25. && fabs(jet.Eta()) < 4.5 && (jet.Pt() > 60. || jet.IsNoPU(m_jvtcut) || fabs(jet.Eta()) > 2.4)) {
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
        } else {
            //std::cout << "removed: " << jet.Pt() << ", " << jet.Eta() << ", " << jet.Phi()<< std::endl;
        }
    }

    /*
    if  ( ( *DatasetID == 361022 && *EventNo == 1218277) ||
            ( *DatasetID == 361022 && *EventNo == 277914) ||
            ( *DatasetID == 361023 && *EventNo == 4886599) ||
            ( *DatasetID == 361023 && *EventNo == 7488620) ||
            ( *DatasetID == 361023 && *EventNo == 608320) ||
            ( *DatasetID == 361023 && *EventNo == 5217909) ) {
        std::cout << "Dataset/runnumber: " << *DatasetID << std::endl;
        std::cout << "*EventNo: " << *EventNo << std::endl;
        if (firstJet) std::cout << "jet 1 (pT, eta, phi): " << firstJet->Pt() << ", " << firstJet->Eta() << ", " << firstJet->Phi()<< std::endl;
        if (secondJet) std::cout << "jet 2 (pT, eta, phi): " << secondJet->Pt() << ", " << secondJet->Eta() << ", " << secondJet->Phi()<< std::endl;
        if (thirdJet) std::cout << "jet 3 (pT, eta, phi): " << thirdJet->Pt() << ", " << thirdJet->Eta() << ", " << thirdJet->Phi()<< std::endl;
        if (fourthJet) std::cout << "jet 4 (pT, eta, phi): " << fourthJet->Pt() << ", " << fourthJet->Eta() << ", " << fourthJet->Phi()<< std::endl;
    	if (secondJet) std::cout << "pTjj: " << (*firstJet + *secondJet).Pt() << std::endl;
    	if (secondJet) std::cout << "Mjj: " << (*firstJet + *secondJet).M() << std::endl;
    	std::cout << "LV: " << lv << std::endl;

    }
    * */

    if (secondJet && cutFlowStudies) {
        N2jets_VR1 += 1;
        N2jets_w_VR1 += eventWeight;
        N2jets_VR2 += 1;
        N2jets_w_VR2 += eventWeight;
        N2jets_VR3 += 1;
        N2jets_w_VR3 += eventWeight;
        if (triggered) {
            N2jets_VR1_trig += 1;
            N2jets_w_VR1_trig += eventWeight;
            N2jets_VR2_trig += 1;
            N2jets_w_VR2_trig += eventWeight;
            N2jets_VR3_trig += 1;
            N2jets_w_VR3_trig += eventWeight;
        }
        if (firstJet->Pt()>80. && secondJet->Pt() > 50.) {
            NjetPt_VR1 += 1;
            NjetPt_w_VR1 += eventWeight;
            NjetPt_VR2 += 1;
            NjetPt_w_VR2 += eventWeight;
            NjetPt_VR3 += 1;
            NjetPt_w_VR3 += eventWeight;
            if (triggered) {
                NjetPt_VR1_trig += 1;
                NjetPt_w_VR1_trig += eventWeight;
                NjetPt_VR2_trig += 1;
                NjetPt_w_VR2_trig += eventWeight;
                NjetPt_VR3_trig += 1;
                NjetPt_w_VR3_trig += eventWeight;
            }
            // VR1
            if (thirdJet && !fourthJet) {
                if (thirdJet->Pt() > 25. && thirdJet->Pt() < 50.) {
                    NTJV_VR1 += 1;
                    NTJV_w_VR1 += eventWeight;
                    if (triggered) {
                        NTJV_VR1_trig += 1;
                        NTJV_w_VR1_trig += eventWeight;
                    }
                    double dPhi = fabs(firstJet->DeltaPhi(*secondJet));
                    if (dPhi < 2.5 && dPhi > 1.8) {
                        NdPhi_VR1 += 1;
                        NdPhi_w_VR1 += eventWeight;
                        if (triggered) {
                            NdPhi_VR1_trig += 1;
                            NdPhi_w_VR1_trig += eventWeight;
                        }
                        double dEta = fabs(firstJet->Eta() - secondJet->Eta());
                        if (dEta > 3.0) {
                            NdEta_VR1 += 1;
                            NdEta_w_VR1 += eventWeight;
                            if (triggered) {
                                NdEta_VR1_trig += 1;
                                NdEta_w_VR1_trig += eventWeight;
                            }
                            if ( (firstJet->Eta() * secondJet->Eta()) < 0 ) {
                                Nhemi_VR1 += 1;
                                Nhemi_w_VR1 += eventWeight;
                                if (triggered) {
                                    Nhemi_VR1_trig += 1;
                                    Nhemi_w_VR1_trig += eventWeight;
                                }
                                if ((*firstJet + *secondJet).Pt() > 150.) {
                                    NpTjj_VR1 += 1;
                                    NpTjj_w_VR1 += eventWeight;
                                    if (triggered) {
                                        NpTjj_VR1_trig += 1;
                                        NpTjj_w_VR1_trig += eventWeight;
                                    }
                                    if ((*firstJet + *secondJet).M() > 1000.) {
                                        Nmjj_VR1 += 1;
                                        Nmjj_w_VR1 += eventWeight;
                                        if (triggered) {
                                            Nmjj_VR1_trig += 1;
                                            Nmjj_w_VR1_trig += eventWeight;
                                        }
                                        if (!lv) {
                                            Nlv_VR1 += 1;
                                            Nlv_w_VR1 += eventWeight;
                                            if (triggered) {
                                                Nlv_VR1_trig += 1;
                                                Nlv_w_VR1_trig += eventWeight;
                                            }
                                            std::cout << "Selected by VR1: " << std::endl;
                                            std::cout << "Dataset/runnumber: " << *DatasetID << std::endl;
                                            std::cout << "EventNo: " << *EventNo << std::endl;
                                            std::cout << "jet 1 (pT, eta, phi): " << firstJet->Pt() << ", " << firstJet->Eta() << ", " << firstJet->Phi() << std::endl;
                                            std::cout << "jet 2 (pT, eta, phi): " << secondJet->Pt() << ", " << secondJet->Eta() << ", " << secondJet->Phi() << std::endl;
                                            std::cout << "jet 3 (pT, eta, phi): " << thirdJet->Pt() << ", " << thirdJet->Eta() << ", " << thirdJet->Phi() << std::endl;
                                            std::cout << "pTjj: " << (*firstJet + *secondJet).Pt() << std::endl;
                                            std::cout << "Mjj: " << (*firstJet + *secondJet).M() << std::endl;
                                            std::cout << "LV: " << lv << std::endl;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else if (thirdJet) { // VR2 and VR3
                if (thirdJet->Pt() > 25. && thirdJet->Pt() < 50.) {
                    NTJV_VR2 += 1;
                    NTJV_w_VR2 += eventWeight;
                    NTJV_VR3 += 1;
                    NTJV_w_VR3 += eventWeight;
                    if (triggered) {
                        NTJV_VR2_trig += 1;
                        NTJV_w_VR2_trig += eventWeight;
                        NTJV_VR3_trig += 1;
                        NTJV_w_VR3_trig += eventWeight;
                    }
                    double dPhi = fabs(firstJet->DeltaPhi(*secondJet));
                    if (dPhi < 2.5) {
                        NdPhi_VR2 += 1;
                        NdPhi_w_VR2 += eventWeight;
                        NdPhi_VR3 += 1;
                        NdPhi_w_VR3 += eventWeight;
                        if (triggered) {
                            NdPhi_VR2_trig += 1;
                            NdPhi_w_VR2_trig += eventWeight;
                            NdPhi_VR3_trig += 1;
                            NdPhi_w_VR3_trig += eventWeight;
                        }
                        double dEta = fabs(firstJet->Eta() - secondJet->Eta());
                        if (dEta > 4.5) {
                            NdEta_VR2 += 1;
                            NdEta_w_VR2 += eventWeight;
                            if (triggered) {
                                NdEta_VR2_trig += 1;
                                NdEta_w_VR2_trig += eventWeight;
                            }
                            if ( (firstJet->Eta() * secondJet->Eta()) < 0 ) {
                                Nhemi_VR2 += 1;
                                Nhemi_w_VR2 += eventWeight;
                                if (triggered) {
                                    Nhemi_VR2_trig += 1;
                                    Nhemi_w_VR2_trig += eventWeight;
                                }
                                if ((*firstJet + *secondJet).Pt() > 150.) {
                                    NpTjj_VR2 += 1;
                                    NpTjj_w_VR2 += eventWeight;
                                    if (triggered) {
                                        NpTjj_VR2_trig += 1;
                                        NpTjj_w_VR2_trig += eventWeight;
                                    }
                                    if ((*firstJet + *secondJet).M() > 1000.) {
                                        Nmjj_VR2 += 1;
                                        Nmjj_w_VR2 += eventWeight;
                                        if (triggered) {
                                            Nmjj_VR2_trig += 1;
                                            Nmjj_w_VR2_trig += eventWeight;
                                        }
                                        if (!lv) {
                                            Nlv_VR2 += 1;
                                            Nlv_w_VR2 += eventWeight;
                                            if (triggered) {
                                                Nlv_VR2_trig += 1;
                                                Nlv_w_VR2_trig += eventWeight;
                                            }
                                            std::cout << "Selected by VR2: " << std::endl;
                                            std::cout << "Dataset/runnumber: " << *DatasetID << std::endl;
                                            std::cout << "EventNo: " << *EventNo << std::endl;
                                            std::cout << "jet 1 (pT, eta, phi): " << firstJet->Pt() << ", " << firstJet->Eta() << ", " << firstJet->Phi() << std::endl;
                                            std::cout << "jet 2 (pT, eta, phi): " << secondJet->Pt() << ", " << secondJet->Eta() << ", " << secondJet->Phi() << std::endl;
                                            std::cout << "jet 3 (pT, eta, phi): " << thirdJet->Pt() << ", " << thirdJet->Eta() << ", " << thirdJet->Phi() << std::endl;
                                            std::cout << "pTjj: " << (*firstJet + *secondJet).Pt() << std::endl;
                                            std::cout << "Mjj: " << (*firstJet + *secondJet).M() << std::endl;
                                            std::cout << "LV: " << lv << std::endl;
                                        }
                                    }
                                }
                            }
                        }

                        if (dEta > 3.0) { // VR3
                            NdEta_VR3 += 1;
                            NdEta_w_VR3 += eventWeight;
                            if (triggered) {
                                NdEta_VR3_trig += 1;
                                NdEta_w_VR3_trig += eventWeight;
                            }
                            if ( (firstJet->Eta() * secondJet->Eta()) < 0 ) {
                                Nhemi_VR3 += 1;
                                Nhemi_w_VR3 += eventWeight;
                                if (triggered) {
                                    Nhemi_VR3_trig += 1;
                                    Nhemi_w_VR3_trig += eventWeight;
                                }
                                if ((*firstJet + *secondJet).Pt() > 150.) {
                                    NpTjj_VR3 += 1;
                                    NpTjj_w_VR3 += eventWeight;
                                    if (triggered) {
                                        NpTjj_VR3_trig += 1;
                                        NpTjj_w_VR3_trig += eventWeight;
                                    }
                                    if ((*firstJet + *secondJet).M() > 1000.) {
                                        Nmjj_VR3 += 1;
                                        Nmjj_w_VR3 += eventWeight;
                                        if (triggered) {
                                            Nmjj_VR3_trig += 1;
                                            Nmjj_w_VR3_trig += eventWeight;
                                        }
                                        if (!lv) {
                                            Nlv_VR3 += 1;
                                            Nlv_w_VR3 += eventWeight;
                                            if (triggered) {
                                                Nlv_VR3_trig += 1;
                                                Nlv_w_VR3_trig += eventWeight;
                                            }
                                            std::cout << "Selected by VR3: " << std::endl;
                                            std::cout << "Dataset/runnumber: " << *DatasetID << std::endl;
                                            std::cout << "EventNo: " << *EventNo << std::endl;
                                            std::cout << "jet 1 (pT, eta, phi): " << firstJet->Pt() << ", " << firstJet->Eta() << ", " << firstJet->Phi() << std::endl;
                                            std::cout << "jet 2 (pT, eta, phi): " << secondJet->Pt() << ", " << secondJet->Eta() << ", " << secondJet->Phi() << std::endl;
                                            std::cout << "jet 3 (pT, eta, phi): " << thirdJet->Pt() << ", " << thirdJet->Eta() << ", " << thirdJet->Phi() << std::endl;
                                            std::cout << "pTjj: " << (*firstJet + *secondJet).Pt() << std::endl;
                                            std::cout << "Mjj: " << (*firstJet + *secondJet).M() << std::endl;
                                            std::cout << "LV: " << lv << std::endl;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //}

    if (isMC) {
        // Loop over all genjets in this container
        TLorentzVector TruthMHT(0.,0.,0.,0.);
        for ( auto& genjet : genJets ) {
            if (genjet.Pt() > 20. && fabs(genjet.Eta()) < 2.4) {
                TruthMHT -= genjet;
            }
        }
        h_TruthMHT->Fill(TruthMHT.Pt(), eventWeight);
	}

        TLorentzVector MET;
        MET.SetPtEtaPhiM(*MET_pt, 0, *MET_phi, 0);
        h_MET->Fill(MET.Pt(), eventWeight);

        TLorentzVector METjet;
        METjet.SetPtEtaPhiM(*METjet_pt, 0, *METjet_phi, 0);

        TLorentzVector METmu;
        METmu.SetPtEtaPhiM(*METmu_pt, 0, *METmu_phi, 0);

        TLorentzVector METele;
        METele.SetPtEtaPhiM(*METele_pt, 0, *METele_phi, 0);

        TLorentzVector METgamma;
        METgamma.SetPtEtaPhiM(*METgamma_pt, 0, *METgamma_phi, 0);

        TLorentzVector METtrack;
        METtrack.SetPtEtaPhiM(*METtrack_pt, 0, *METtrack_phi, 0);

        TLorentzVector genMET;
        genMET.SetPtEtaPhiM(*GenMET_pt, 0, *GenMET_phi, 0);
        h_TruthMET->Fill(genMET.Pt(), eventWeight);

        TLorentzVector MHT(0.,0.,0.,0.);
        TLorentzVector MHTnoOR(0.,0.,0.,0.);
        TLorentzVector MHTnoJVT(0.,0.,0.,0.);
        TLorentzVector MHTnoJVTnoOR(0.,0.,0.,0.);
        TLorentzVector JVTjets(0.,0.,0.,0.);
        TLorentzVector TruthMatchedMHT(0.,0.,0.,0.);
        TLorentzVector TruthMatchedMHTnoOR(0.,0.,0.,0.);
        TLorentzVector TruthMatchedMHTnoJVT(0.,0.,0.,0.);
        TLorentzVector TruthMatchedMHTnoJVTnoOR(0.,0.,0.,0.);

        int JetCount = 0;
        int GoodJetCount = 0;

        firstJet = 0;
        secondJet = 0;

        for ( auto& jet : recoJets) {

            bool isNoPU = false;
            if (jet.Pt() > 60. || jet.IsNoPU(m_jvtcut) || fabs(jet.Eta()) > 2.4) isNoPU = true;

            if (jet.IsGood() && isNoPU) {

                double dphi = fabs(MET.DeltaPhi(jet));
                if (JetCount == 0) {
                    h_Jet1_Pt->Fill(jet.Pt(),eventWeight);
                    h_Jet1_Eta->Fill(jet.Eta(),eventWeight);
                    h_Jet1_Phi->Fill(jet.Phi(),eventWeight);
                    h_Jet1_DeltaPhi->Fill(dphi ,eventWeight);
                }

                if (JetCount == 1) {
                    h_Jet2_Pt->Fill(jet.Pt(),eventWeight);
                    h_Jet2_Eta->Fill(jet.Eta(),eventWeight);
                    h_Jet2_Phi->Fill(jet.Phi(),eventWeight);
                    h_Jet2_DeltaPhi->Fill(dphi ,eventWeight);
                }

                if (JetCount == 2) {
                    h_Jet3_Pt->Fill(jet.Pt(),eventWeight);
                    h_Jet3_Eta->Fill(jet.Eta(),eventWeight);
                    h_Jet3_Phi->Fill(jet.Phi(),eventWeight);
                    h_Jet3_DeltaPhi->Fill(dphi ,eventWeight);
                }

            }

            if (jet.IsGood()) MHTnoJVTnoOR -= jet;
            if (jet.IsGood() && jet.PassOR() ) MHTnoJVT -= jet;
            if (jet.IsGood() && isNoPU) MHTnoOR -= jet;
            if (jet.IsGood() && jet.PassOR() && isNoPU) MHT -= jet;

            if (jet.IsGood()) {

                if (GoodJetCount == 0) firstJet = &jet;
                if (GoodJetCount == 1) secondJet = &jet;

                //// Match reco jets to truth jets
                MyJet* matchedJet = 0;
                double dRmin = 999.;
                for ( auto& genjet : genJets ) {
                    if (genjet.DeltaR(jet) < dRmin) {
                        matchedJet = &genjet;
                        dRmin = genjet.DeltaR(jet);
                    }
                }

                //// If match is found replace reco jet by truth jet
                if (dRmin < 0.15) {
                    if (jet.IsGood()) TruthMatchedMHTnoJVTnoOR -= *matchedJet;
                    if (jet.IsGood() && jet.PassOR() ) TruthMatchedMHTnoJVT -= *matchedJet;
                    if (jet.IsGood() && isNoPU) TruthMatchedMHTnoOR -= *matchedJet;
                    if (jet.IsGood() && jet.PassOR() && isNoPU) TruthMatchedMHT -= *matchedJet;
                } else {
                    if (jet.IsGood()) TruthMatchedMHTnoJVTnoOR -= jet;
                    if (jet.IsGood() && jet.PassOR() ) TruthMatchedMHTnoJVT -= jet;
                    if (jet.IsGood() && isNoPU) TruthMatchedMHTnoOR -= jet;
                    if (jet.IsGood() && jet.PassOR() && isNoPU) TruthMatchedMHT -= jet;
                }
                ++GoodJetCount;

            }

            if (jet.IsGood() && jet.Pt() > 25. && !isNoPU) {
                JVTjets += jet;
            }

            ++JetCount;
        }

        for ( auto& mu : recoMuons) {
            if (mu.PassOR()) {
                MHT-=mu;
                MHTnoJVT-=mu;
            }
            MHTnoOR-=mu;
            MHTnoJVTnoOR-=mu;
        }

        for ( auto& ele : recoElectrons) {
            if (ele.PassOR()) {
                MHT-=ele;
                MHTnoJVT-=ele;
            }
        }

		/*
        for ( auto& pho : recoPhotons) {
            if (pho.PassOR()) {
                MHT-=pho;
                MHTnoJVT-=pho;
            }
        }
        */

        if (GoodJetCount > 1) {

            //if (firstJet->Pt() > 80. && secondJet->Pt() > 50. && fabs(firstJet->Eta() - secondJet->Eta()) > 3.5 && (firstJet->Eta()*secondJet->Eta()) < 0. ) {
            //if ( MHTnoJVT.Pt() > 150 || MHTnoJVTnoOR.Pt() > 150. ) {
            //if ( fabs(MET.Pt() - MHTnoOR.Pt()) > 200) {
            //if ( METgamma.Pt() > 100 ) {
            //if ( MET.Pt() > 150 && ( (METgamma.Pt() > 200) || (METtrack.Pt() > 200) || (JVTjets.Pt() > 200) ) ){
            //if ( MET.Pt() > 150 && ( (JVTjets.Pt() > 100) ) ){
            //if ( MET.Pt() > 400 && METele.Pt() > 400.){
            //if ( MHT.Pt() < 100 && MET.Pt() > 200 && ( (METtrack.Pt() > 300) ) ) {
            //if ( MHT.Pt() > 150 && METele.Pt() < 20. && METmu.Pt() < 20. ) {
			if ( JVTjets.Pt() > 11150 ) {

                std::cout << "----------------------" << std::endl;
                std::cout << "EventNo, DatasetID/RunNo: " << *EventNo << ", " << *DatasetID << std::endl;
                std::cout << "MET            = " << MET.Pt() << ", " << MET.Phi() << std::endl;
                std::cout << "METjet         = " << METjet.Pt() << ", " << METjet.Phi() << std::endl;
                std::cout << "METmu          = " << METmu.Pt() << ", " << METmu.Phi() << std::endl;
                std::cout << "METele         = " << METele.Pt() << ", " << METele.Phi() << std::endl;
                std::cout << "METgamma       = " << METgamma.Pt() << ", " << METgamma.Phi() << std::endl;
                std::cout << "METtrack       = " << METtrack.Pt() << ", " << METtrack.Phi() << std::endl;
                std::cout << "genMET         = " << genMET.Pt() << ", " << genMET.Phi() << std::endl;
                std::cout << "MHT (j,e,mu,g)        = " << MHT.Pt() << ", " << MHT.Phi() << std::endl;
                std::cout << "MHT (j no OR)         = " << MHTnoOR.Pt() << ", " << MHTnoOR.Phi() << std::endl;
                std::cout << "MHT (j no OR, no jvt) = " << MHTnoJVTnoOR.Pt() << ", " << MHTnoJVTnoOR.Phi() << std::endl;
                for ( auto& genjet : genJets ) {
                    if (genjet.Pt() > 20.) {
                        std::cout << "genjet (pt, eta, phi): " << genjet.Pt() << ", " << genjet.Eta() << ", " << genjet.Phi()<< std::endl;
                    }
                }
                for ( auto& jet : recoJets) {
                    if (jet.IsGood() && jet.Pt() > 20. && fabs(jet.Eta()) < 4.5 && (jet.Pt() > 60. || jet.IsNoPU(m_jvtcut) || fabs(jet.Eta()) > 2.4)) {
                        std::cout << "jet (pt, eta, phi, fwJVT, passOR): " << jet.Pt() << ", " << jet.Eta() << ", " << jet.Phi() << ", " << jet.IsFJVT() << ", " << jet.PassOR() << std::endl;
                    } else {
                        std::cout << "removed jet (pt, eta, phi, fwJVT, passOR): " << jet.Pt() << ", " << jet.Eta() << ", " << jet.Phi()  << ", " << jet.IsFJVT() << ", " << jet.PassOR() << std::endl;
                    }
                }
                for ( auto& ele : recoElectrons) {
                    std::cout << "electron (pt, eta, phi, signal, passOR): " << ele.Pt() << ", " << ele.Eta() << ", " << ele.Phi() << ", " << ele.IsSignal() << ", " << ele.PassOR() << std::endl;
                }
                for ( auto& muon : recoMuons) {
                    std::cout << "muon (pt, eta, phi, signal, passOR): " << muon.Pt() << ", " << muon.Eta() << ", " << muon.Phi() << ", " << muon.IsSignal() << ", " << muon.PassOR() << std::endl;
                }
                for ( auto& tau : recoTaus) {
                    std::cout << "tau (pt, eta, phi, signal, passOR): " << tau.Pt() << ", " << tau.Eta() << ", " << tau.Phi() << ", " << tau.IsSignal() << ", " << tau.PassOR() << std::endl;
                }
                for ( auto& gamma : recoPhotons) {
                    std::cout << "photon (pt, eta, phi, signal, passOR): " << gamma.Pt() << ", " << gamma.Eta() << ", " << gamma.Phi() << ", " << gamma.IsSignal() << ", " << gamma.PassOR() << std::endl;
                }
            }

            TLorentzVector METreplaced(0.,0.,0.,0.);
            METreplaced = MET - MHT + TruthMatchedMHT;
            h_METreplaced->Fill(METreplaced.Pt(), eventWeight);

			h_MHT->Fill(MHT.Pt(), eventWeight);
            h_MHTnoOR->Fill(MHTnoOR.Pt(), eventWeight);
            h_MHTnoJVT->Fill(MHTnoJVT.Pt(), eventWeight);
            h_MHTnoJVTnoOR->Fill(MHTnoJVTnoOR.Pt(), eventWeight);

            h_MET_vs_MHT->Fill(MET.Pt(), MHT.Pt(), eventWeight);
            h_MET_vs_MHTnoOR->Fill(MET.Pt(), MHTnoOR.Pt(), eventWeight);
            h_MET_vs_MHTnoJVT->Fill(MET.Pt(), MHTnoJVT.Pt(), eventWeight);
            h_MET_vs_MHTnoJVTnoOR->Fill(MET.Pt(), MHTnoJVTnoOR.Pt(), eventWeight);
            
            h_MHT_vs_MHTnoOR->Fill(MHT.Pt(), MHTnoJVT.Pt(), eventWeight);
            h_MHT_vs_MHTnoJVT->Fill(MHT.Pt(), MHTnoJVT.Pt(), eventWeight);
            h_MHT_vs_MHTnoJVTnoOR->Fill(MHT.Pt(), MHTnoJVTnoOR.Pt(), eventWeight);
            
            h_JVTPhi_vs_METPhi->Fill(JVTjets.Phi(), MET.Phi(), eventWeight);
            h_JVTPhi_vs_MHTPhi->Fill(JVTjets.Phi(), MHT.Phi(), eventWeight);
            h_JVTPhi_vs_MHTnoORPhi->Fill(JVTjets.Phi(), MHTnoOR.Phi(), eventWeight);
            
            h_MET_vs_METreplaced->Fill(MET.Pt(), METreplaced.Pt(), eventWeight);
            h_MHTnoOR_vs_MHTnoORreplaced->Fill(MHTnoOR.Pt(), TruthMatchedMHTnoOR.Pt(), eventWeight);
            h_MHTnoJVTnoOR_vs_MHTnoJVTnoORreplaced->Fill(MHTnoJVTnoOR.Pt(), TruthMatchedMHTnoJVTnoOR.Pt(), eventWeight);

            h_MET_vs_ElePt->Fill(MET.Pt(), PtEle, eventWeight);
            h_MET_vs_MuPt->Fill(MET.Pt(), PtMu, eventWeight);
            h_MET_vs_PhoPt->Fill(MET.Pt(), PtPho, eventWeight);
            h_MET_vs_LepPt->Fill(MET.Pt(), (PtEle>PtMu ? PtEle : PtMu), eventWeight);

            h_MET_vs_METjet->Fill(MET.Pt(), METjet.Pt(), eventWeight);
            h_MET_vs_METele->Fill(MET.Pt(), METele.Pt(), eventWeight);
            h_MET_vs_METmu->Fill(MET.Pt(), METmu.Pt(), eventWeight);
            h_MET_vs_METgamma->Fill(MET.Pt(), METgamma.Pt(), eventWeight);
            h_MET_vs_METtrack->Fill(MET.Pt(), METtrack.Pt(), eventWeight);

        }

    //} // isMC

    return kTRUE;
}

void MyMETStudies::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.


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

    /*
    std::cout << "Selected by VR1 (no trigger): " << std::endl;
    std::cout << std::endl;
    std::cout << "Ntot   : " << Ntot_VR1 << std::endl;
    std::cout << "N2jets : " << N2jets_VR1 << std::endl;
    std::cout << "NjetPt : " << NjetPt_VR1 << std::endl;
    std::cout << "NTJV   : " << NTJV_VR1 << std::endl;
    std::cout << "NdPhi  : " << NdPhi_VR1 << std::endl;
    std::cout << "NdEta  : " << NdEta_VR1 << std::endl;
    std::cout << "Nhemi  : " << Nhemi_VR1 << std::endl;
    std::cout << "NpTjj  : " << NpTjj_VR1 << std::endl;
    std::cout << "Nmjj   : " << Nmjj_VR1 << std::endl;
    std::cout << "Nlv    : " << Nlv_VR1 << std::endl;
    std::cout << std::endl;
    std::cout << "Ntot_w   : " << Ntot_w_VR1 << std::endl;
    std::cout << "N2jets_w : " << N2jets_w_VR1 << std::endl;
    std::cout << "NjetPt_w : " << NjetPt_w_VR1 << std::endl;
    std::cout << "NTJV_w   : " << NTJV_w_VR1 << std::endl;
    std::cout << "NdPhi_w  : " << NdPhi_w_VR1 << std::endl;
    std::cout << "NdEta_w  : " << NdEta_w_VR1 << std::endl;
    std::cout << "Nhemi_w  : " << Nhemi_w_VR1 << std::endl;
    std::cout << "NpTjj_w  : " << NpTjj_w_VR1 << std::endl;
    std::cout << "Nmjj_w   : " << Nmjj_w_VR1 << std::endl;
    std::cout << "Nlv_w    : " << Nlv_w_VR1 << std::endl;
    std::cout << std::endl;
    std::cout << "Selected by VR1 (with trigger): " << std::endl;
    std::cout << std::endl;
    std::cout << "Ntot   : " << Ntot_VR1_trig << std::endl;
    std::cout << "N2jets : " << N2jets_VR1_trig << std::endl;
    std::cout << "NjetPt : " << NjetPt_VR1_trig << std::endl;
    std::cout << "NTJV   : " << NTJV_VR1_trig << std::endl;
    std::cout << "NdPhi  : " << NdPhi_VR1_trig << std::endl;
    std::cout << "NdEta  : " << NdEta_VR1_trig << std::endl;
    std::cout << "Nhemi  : " << Nhemi_VR1_trig << std::endl;
    std::cout << "NpTjj  : " << NpTjj_VR1_trig << std::endl;
    std::cout << "Nmjj   : " << Nmjj_VR1_trig << std::endl;
    std::cout << "Nlv    : " << Nlv_VR1_trig << std::endl;
    std::cout << std::endl;
    std::cout << "Ntot_w   : " << Ntot_w_VR1_trig << std::endl;
    std::cout << "N2jets_w : " << N2jets_w_VR1_trig << std::endl;
    std::cout << "NjetPt_w : " << NjetPt_w_VR1_trig << std::endl;
    std::cout << "NTJV_w   : " << NTJV_w_VR1_trig << std::endl;
    std::cout << "NdPhi_w  : " << NdPhi_w_VR1_trig << std::endl;
    std::cout << "NdEta_w  : " << NdEta_w_VR1_trig << std::endl;
    std::cout << "Nhemi_w  : " << Nhemi_w_VR1_trig << std::endl;
    std::cout << "NpTjj_w  : " << NpTjj_w_VR1_trig << std::endl;
    std::cout << "Nmjj_w   : " << Nmjj_w_VR1_trig << std::endl;
    std::cout << "Nlv_w    : " << Nlv_w_VR1_trig << std::endl;
    std::cout << std::endl;
    std::cout << "Selected by VR2 (no trigger): " << std::endl;
    std::cout << std::endl;
    std::cout << "Ntot   : " << Ntot_VR2 << std::endl;
    std::cout << "N2jets : " << N2jets_VR2 << std::endl;
    std::cout << "NjetPt : " << NjetPt_VR2 << std::endl;
    std::cout << "NTJV   : " << NTJV_VR2 << std::endl;
    std::cout << "NdPhi  : " << NdPhi_VR2 << std::endl;
    std::cout << "NdEta  : " << NdEta_VR2 << std::endl;
    std::cout << "Nhemi  : " << Nhemi_VR2 << std::endl;
    std::cout << "NpTjj  : " << NpTjj_VR2 << std::endl;
    std::cout << "Nmjj   : " << Nmjj_VR2 << std::endl;
    std::cout << "Nlv    : " << Nlv_VR2 << std::endl;
    std::cout << std::endl;
    std::cout << "Ntot_w   : " << Ntot_w_VR2 << std::endl;
    std::cout << "N2jets_w : " << N2jets_w_VR2 << std::endl;
    std::cout << "NjetPt_w : " << NjetPt_w_VR2 << std::endl;
    std::cout << "NTJV_w   : " << NTJV_w_VR2 << std::endl;
    std::cout << "NdPhi_w  : " << NdPhi_w_VR2 << std::endl;
    std::cout << "NdEta_w  : " << NdEta_w_VR2 << std::endl;
    std::cout << "Nhemi_w  : " << Nhemi_w_VR2 << std::endl;
    std::cout << "NpTjj_w  : " << NpTjj_w_VR2 << std::endl;
    std::cout << "Nmjj_w   : " << Nmjj_w_VR2 << std::endl;
    std::cout << "Nlv_w    : " << Nlv_w_VR2 << std::endl;
    std::cout << std::endl;
    std::cout << "Selected by VR2 (with trigger): " << std::endl;
    std::cout << std::endl;
    std::cout << "Ntot   : " << Ntot_VR2_trig << std::endl;
    std::cout << "N2jets : " << N2jets_VR2_trig << std::endl;
    std::cout << "NjetPt : " << NjetPt_VR2_trig << std::endl;
    std::cout << "NTJV   : " << NTJV_VR2_trig << std::endl;
    std::cout << "NdPhi  : " << NdPhi_VR2_trig << std::endl;
    std::cout << "NdEta  : " << NdEta_VR2_trig << std::endl;
    std::cout << "Nhemi  : " << Nhemi_VR2_trig << std::endl;
    std::cout << "NpTjj  : " << NpTjj_VR2_trig << std::endl;
    std::cout << "Nmjj   : " << Nmjj_VR2_trig << std::endl;
    std::cout << "Nlv    : " << Nlv_VR2_trig << std::endl;
    std::cout << std::endl;
    std::cout << "Ntot_w   : " << Ntot_w_VR2_trig << std::endl;
    std::cout << "N2jets_w : " << N2jets_w_VR2_trig << std::endl;
    std::cout << "NjetPt_w : " << NjetPt_w_VR2_trig << std::endl;
    std::cout << "NTJV_w   : " << NTJV_w_VR2_trig << std::endl;
    std::cout << "NdPhi_w  : " << NdPhi_w_VR2_trig << std::endl;
    std::cout << "NdEta_w  : " << NdEta_w_VR2_trig << std::endl;
    std::cout << "Nhemi_w  : " << Nhemi_w_VR2_trig << std::endl;
    std::cout << "NpTjj_w  : " << NpTjj_w_VR2_trig << std::endl;
    std::cout << "Nmjj_w   : " << Nmjj_w_VR2_trig << std::endl;
    std::cout << "Nlv_w    : " << Nlv_w_VR2_trig << std::endl;
    std::cout << std::endl;
    std::cout << "Selected by VR3 (no trigger): " << std::endl;
    std::cout << std::endl;
    std::cout << "Ntot   : " << Ntot_VR3 << std::endl;
    std::cout << "N2jets : " << N2jets_VR3 << std::endl;
    std::cout << "NjetPt : " << NjetPt_VR3 << std::endl;
    std::cout << "NTJV   : " << NTJV_VR3 << std::endl;
    std::cout << "NdPhi  : " << NdPhi_VR3 << std::endl;
    std::cout << "NdEta  : " << NdEta_VR3 << std::endl;
    std::cout << "Nhemi  : " << Nhemi_VR3 << std::endl;
    std::cout << "NpTjj  : " << NpTjj_VR3 << std::endl;
    std::cout << "Nmjj   : " << Nmjj_VR3 << std::endl;
    std::cout << "Nlv    : " << Nlv_VR3 << std::endl;
    std::cout << std::endl;
    std::cout << "Ntot_w   : " << Ntot_w_VR3 << std::endl;
    std::cout << "N2jets_w : " << N2jets_w_VR3 << std::endl;
    std::cout << "NjetPt_w : " << NjetPt_w_VR3 << std::endl;
    std::cout << "NTJV_w   : " << NTJV_w_VR3 << std::endl;
    std::cout << "NdPhi_w  : " << NdPhi_w_VR3 << std::endl;
    std::cout << "NdEta_w  : " << NdEta_w_VR3 << std::endl;
    std::cout << "Nhemi_w  : " << Nhemi_w_VR3 << std::endl;
    std::cout << "NpTjj_w  : " << NpTjj_w_VR3 << std::endl;
    std::cout << "Nmjj_w   : " << Nmjj_w_VR3 << std::endl;
    std::cout << "Nlv_w    : " << Nlv_w_VR3 << std::endl;
    std::cout << std::endl;
    std::cout << "Selected by VR3 (with trigger): " << std::endl;
    std::cout << std::endl;
    std::cout << "Ntot   : " << Ntot_VR3_trig << std::endl;
    std::cout << "N2jets : " << N2jets_VR3_trig << std::endl;
    std::cout << "NjetPt : " << NjetPt_VR3_trig << std::endl;
    std::cout << "NTJV   : " << NTJV_VR3_trig << std::endl;
    std::cout << "NdPhi  : " << NdPhi_VR3_trig << std::endl;
    std::cout << "NdEta  : " << NdEta_VR3_trig << std::endl;
    std::cout << "Nhemi  : " << Nhemi_VR3_trig << std::endl;
    std::cout << "NpTjj  : " << NpTjj_VR3_trig << std::endl;
    std::cout << "Nmjj   : " << Nmjj_VR3_trig << std::endl;
    std::cout << "Nlv    : " << Nlv_VR3_trig << std::endl;
    std::cout << std::endl;
    std::cout << "Ntot_w   : " << Ntot_w_VR3_trig << std::endl;
    std::cout << "N2jets_w : " << N2jets_w_VR3_trig << std::endl;
    std::cout << "NjetPt_w : " << NjetPt_w_VR3_trig << std::endl;
    std::cout << "NTJV_w   : " << NTJV_w_VR3_trig << std::endl;
    std::cout << "NdPhi_w  : " << NdPhi_w_VR3_trig << std::endl;
    std::cout << "NdEta_w  : " << NdEta_w_VR3_trig << std::endl;
    std::cout << "Nhemi_w  : " << Nhemi_w_VR3_trig << std::endl;
    std::cout << "NpTjj_w  : " << NpTjj_w_VR3_trig << std::endl;
    std::cout << "Nmjj_w   : " << Nmjj_w_VR3_trig << std::endl;
    std::cout << "Nlv_w    : " << Nlv_w_VR3_trig << std::endl;
    */

}

void MyMETStudies::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

}

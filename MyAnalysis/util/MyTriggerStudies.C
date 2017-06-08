#define MyTriggerStudies_cxx
// The class definition in MyTriggerStudies.h has been generated automatically
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
// root> T->Process("MyTriggerStudies.C")
// root> T->Process("MyTriggerStudies.C","some options")
// root> T->Process("MyTriggerStudies.C+")
//

#include "MyJet.h"
#include "MyTriggerStudies.h"

void MyTriggerStudies::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

    outputfile = new TFile("TriggerStudiesOutput_VBF_mc.root","RECREATE");

    m_jvtcut = 0.59;
    m_lumi = 30000.;
    isData = true;

    ////Book histograms

    h_MHT_all = new TH1F("h_MHT_all", "h_MHT_all", 100, 0., 500.);
    h_MHT_all->Sumw2();
    histos_1D.push_back(h_MHT_all);

    h_MHT_triggered = new TH1F("h_MHT_triggered", "h_MHT_triggered", 100, 0., 500.);
    h_MHT_triggered->Sumw2();
    histos_1D.push_back(h_MHT_triggered);

    h_MHTvsHT_all = new TH2F("h_MHTvsHT_all", "h_MHTvsHT_all", 100, 0., 500., 30, 0., 3000.);
    h_MHTvsHT_all->Sumw2();
    histos_2D.push_back(h_MHTvsHT_all);

    h_MHTvsHT_triggered = new TH2F("h_MHTvsHT_triggered", "h_MHTvsHT_triggered", 100, 0., 500., 30, 0., 3000.);
    h_MHTvsHT_triggered->Sumw2();
    histos_2D.push_back(h_MHTvsHT_triggered);

    //NTotEvents = fChain->GetEntries();

    //// Not very elegant! TODO: Store this info in and read from file

    // [v2]
    AvailableEvents[361022] = 499750;
    AvailableEvents[361023] = 6445895;
    AvailableEvents[361024] = 7439800;
    AvailableEvents[361025] = 7619000;
    AvailableEvents[361026] = 1683400;

}

void MyTriggerStudies::SlaveBegin(TTree * /*tree*/)
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

}

Bool_t MyTriggerStudies::Process(Long64_t entry)
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
    if (NEvents%1000 == 0) std::cout << NEvents << " processed!" << std::endl;

    if (ProcessedEvents.find(*DatasetID) == ProcessedEvents.end()) {
        ProcessedEvents[*DatasetID] = 1;
    } else {
        ProcessedEvents[*DatasetID] += 1;
    }

    //std::cout << "Weight: " << *Weight << std::endl;
    double eventWeight = *Weight;
    if (!isData) {
        eventWeight *= m_lumi / AvailableEvents[*DatasetID];
    }

    std::vector<MyJet> genJets;
    std::vector<MyJet> recoJets;
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
        MyJet jet(pt, eta, phi,m);
        jet.SetJVT(jvt);
        jet.SetBTag(btag);
        jet.SetJetID(good);

        if (jet.Pt() > 20. && jet.Pt() < 60. && fabs(jet.Eta()) < 2.4) {
            if (jet.IsBad() && jet.IsNoPU(m_jvtcut)) {
                std::cout << "Reject event because of bad central jet!" << std::endl;
                return 1;
            }
        }
        if (jet.Pt() > 20. && jet.Pt() < 60. && fabs(jet.Eta()) >= 2.4) {
            if (jet.IsBad()) {
                std::cout << "Reject event because of bad forward jet!" << std::endl;
                return 1;
            }
        }
        if (jet.Pt() > 60.) {
            if (jet.IsBad()) {
                std::cout << "Reject event because of bad high pT jet!" << std::endl;
                return 1;
            }
        }

        recoJets.push_back(jet);
        //if ( jet.IsNoPU(m_jvtcut) || jet.Pt() > 60. || fabs(jet.Eta()) > 2.4 ) recoJets.push_back(jet);

    }
    GreaterByPt<MyJet> ptComparator_;
    std::sort(recoJets.begin(), recoJets.end(), ptComparator_);

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
        if (MuonIsSignal->at(i)) {
            std::cout << "Reject event because of isolated muon!" << std::endl;
            return 1;
        }
        TLorentzVector muon(0.,0.,0.,0.);
        muon.SetPtEtaPhiM(pt, eta, phi, 0.1057);
        recoMuons.push_back(muon);
    }
    GreaterByPt<TLorentzVector> ptComparator2_;
    std::sort(recoMuons.begin(), recoMuons.end(), ptComparator2_);
    double PtMu = 0;
    if (NMuons > 0) PtMu = recoMuons.at(0).Pt();

    int NElectrons = ElePt->size();
    //std::cout << "NElectrons: " << NElectrons << std::endl;
    for (int i = 0; i < NElectrons; ++i) {
        if (EleIsSignal->at(i)) {
            std::cout << "Reject event because of isolated electron!" << std::endl;
            return 1;
        }
        float pt = ElePt->at(i);
        float eta = EleEta->at(i);
        float phi = ElePhi->at(i);
        TLorentzVector electron(0.,0.,0.,0.);
        electron.SetPtEtaPhiM(pt, eta, phi, 0.000511);
        recoElectrons.push_back(electron);
    }
    std::sort(recoElectrons.begin(), recoElectrons.end(), ptComparator2_);
    double PtEle = 0;
    if (NElectrons > 0) PtEle = recoElectrons.at(0).Pt();

    int NPhotons = PhotonPt->size();
    //std::cout << "NPhoton: " << NPhotons << std::endl;
    for (int i = 0; i < NPhotons; ++i) {
        float pt = PhotonPt->at(i);
        float eta = PhotonEta->at(i);
        float phi = PhotonPhi->at(i);
        TLorentzVector photon(0.,0.,0.,0.);
        photon.SetPtEtaPhiM(pt, eta, phi, 0.);
        recoPhotons.push_back(photon);
    }
    std::sort(recoPhotons.begin(), recoPhotons.end(), ptComparator2_);
    double PtPho = 0;
    if (NPhotons > 0) PtPho = recoPhotons.at(0).Pt();

    TLorentzVector MHT(0.,0.,0.,0.);
    TLorentzVector NoJVTMHT(0.,0.,0.,0.);
    double HT = 0;
    double NoJVTHT = 0;

    int JetCount = 0;
    int NoJVTJetCount = 0;

    for ( auto& jet : recoJets) {

        if (jet.Pt() < 30.) continue;

        if (jet.IsGood()) {
            ++NoJVTJetCount;
            NoJVTMHT -= jet;
            NoJVTHT += jet.Pt();
        }
        if (jet.IsGood() && (jet.Pt() > 60. || jet.IsNoPU(m_jvtcut) || fabs(jet.Eta()) > 2.4)) {
            ++JetCount;
            MHT-=jet;
            HT += jet.Pt();
        }

    }

    if (JetCount == 2) {
        h_MHT_all->Fill(MHT.Pt());
        h_MHTvsHT_all->Fill(MHT.Pt(),HT);
        if (*xe90triggered || *xe110triggered) {
            h_MHT_triggered->Fill(MHT.Pt());
            h_MHTvsHT_triggered->Fill(MHT.Pt(),HT);
        }
    }


    return kTRUE;
}

void MyTriggerStudies::SlaveTerminate()
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
}

void MyTriggerStudies::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

}

#define MyResolution_cxx
// The class definition in MyResolution.h has been generated automatically
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
// root> T->Process("MyResolution.C")
// root> T->Process("MyResolution.C","some options")
// root> T->Process("MyResolution.C+")
//

#include "MyJet.h"
#include "MyMuon.h"
#include "MyResolution.h"

void MyResolution::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

    outputfile = new TFile("resolutions.root","RECREATE");

    m_VetoCone = 0.8;
    m_addActivityCone = 0.6;
    m_MatchingCone = 0.1;
    m_METmuCone = 0.1;
    m_RelGenActivityVeto = 0.01;
    m_RelRecoActivityVeto = 0.05;
    //m_RelRecoActivityVeto = 999.;
    m_jvtcut = 0.59;
    m_lumi = 36100.;
    PtBinEdges.push_back(0);
    PtBinEdges.push_back(10);
    PtBinEdges.push_back(20);
    PtBinEdges.push_back(30);
    PtBinEdges.push_back(40);
    PtBinEdges.push_back(50);
    PtBinEdges.push_back(70);
    PtBinEdges.push_back(100);
    PtBinEdges.push_back(140);
    PtBinEdges.push_back(190);
    PtBinEdges.push_back(250);
    PtBinEdges.push_back(320);
    PtBinEdges.push_back(400);
    PtBinEdges.push_back(490);
    PtBinEdges.push_back(590);
    PtBinEdges.push_back(700);
    PtBinEdges.push_back(820);
    PtBinEdges.push_back(950);
    PtBinEdges.push_back(1090);
    PtBinEdges.push_back(1240);
    PtBinEdges.push_back(1400);
    PtBinEdges.push_back(1570);
    PtBinEdges.push_back(1750);
    PtBinEdges.push_back(1940);
    PtBinEdges.push_back(2140);
    PtBinEdges.push_back(2350);
    PtBinEdges.push_back(2600);
    PtBinEdges.push_back(3000);
    /*
    PtBinEdges.push_back(0);
    PtBinEdges.push_back(20);
    PtBinEdges.push_back(30);
    PtBinEdges.push_back(40);
    PtBinEdges.push_back(50);
    PtBinEdges.push_back(70);
    PtBinEdges.push_back(90);
    PtBinEdges.push_back(110);
    PtBinEdges.push_back(140);
    PtBinEdges.push_back(170);
    PtBinEdges.push_back(200);
    PtBinEdges.push_back(240);
    PtBinEdges.push_back(280);
    PtBinEdges.push_back(330);
    PtBinEdges.push_back(380);
    PtBinEdges.push_back(440);
    PtBinEdges.push_back(500);
    PtBinEdges.push_back(570);
    PtBinEdges.push_back(640);
    PtBinEdges.push_back(720);
    PtBinEdges.push_back(810);
    PtBinEdges.push_back(910);
    PtBinEdges.push_back(1020);
    PtBinEdges.push_back(1140);
    PtBinEdges.push_back(1270);
    PtBinEdges.push_back(1410);
    PtBinEdges.push_back(1660);
    PtBinEdges.push_back(1820);
    PtBinEdges.push_back(1990);
    PtBinEdges.push_back(2170);
    PtBinEdges.push_back(2360);
    PtBinEdges.push_back(2560);
    PtBinEdges.push_back(3000);
    PtBinEdges.push_back(9999);
    */

    EtaBinEdges.push_back(0.0);
    EtaBinEdges.push_back(0.7);
    EtaBinEdges.push_back(1.3);
    EtaBinEdges.push_back(1.8);
    EtaBinEdges.push_back(2.5);
    EtaBinEdges.push_back(3.2);
    EtaBinEdges.push_back(5.0);

    //// Array of histograms for jet resolutions (all jet multiplicities)
    ResizeHistoVector(PtResolution_tot);
    ResizeHistoVector(EtaResolution_tot);
    ResizeHistoVector(PhiResolution_tot);
    ResizeHisto2Vector(MuRes_tot);
    ResizeHistoVector(PtResolution_LF);
    ResizeHistoVector(EtaResolution_LF);
    ResizeHistoVector(PhiResolution_LF);
    ResizeHisto2Vector(MuRes_LF);
    ResizeHistoVector(PtResolution_HF);
    ResizeHistoVector(EtaResolution_HF);
    ResizeHistoVector(PhiResolution_HF);
    ResizeHisto2Vector(MuRes_HF);
    ResizeHistoVector(PtResolution_nob);
    ResizeHistoVector(EtaResolution_nob);
    ResizeHistoVector(PhiResolution_nob);
    ResizeHisto2Vector(MuRes_nob);
    ResizeHistoVector(PtResolution_b);
    ResizeHistoVector(EtaResolution_b);
    ResizeHistoVector(PhiResolution_b);
    ResizeHisto2Vector(MuRes_b);

    for (unsigned int i_pt = 0; i_pt < PtBinEdges.size() - 1; ++i_pt) {
        for (unsigned int i_eta = 0; i_eta < EtaBinEdges.size() - 1; ++i_eta) {

            //// Book histograms Pt response
            TH1F* h_jetRes_tot_pt = new TH1F(GetHistName(i_pt, i_eta,"tot","Pt").c_str(), "p_T^{reco}/p_T^{gen}", 300, 0.0, 3.0);
            PtResolution_tot.at(i_pt).at(i_eta) = h_jetRes_tot_pt;

            TH1F* h_jetRes_LF_pt = new TH1F(GetHistName(i_pt, i_eta,"LF","Pt").c_str(), "p_T^{reco}/p_T^{gen}", 300, 0.0, 3.0);
            PtResolution_LF.at(i_pt).at(i_eta) = h_jetRes_LF_pt;

            TH1F* h_jetRes_HF_pt = new TH1F(GetHistName(i_pt, i_eta,"HF","Pt").c_str(), "p_T^{reco}/p_T^{gen}", 300, 0.0, 3.0);
            PtResolution_HF.at(i_pt).at(i_eta) = h_jetRes_HF_pt;

            TH1F* h_jetRes_nob_pt = new TH1F(GetHistName(i_pt, i_eta,"nob","Pt").c_str(), "p_T^{reco}/p_T^{gen}", 300, 0.0, 3.0);
            PtResolution_nob.at(i_pt).at(i_eta) = h_jetRes_nob_pt;

            TH1F* h_jetRes_b_pt = new TH1F(GetHistName(i_pt, i_eta,"b","Pt").c_str(), "p_T^{reco}/p_T^{gen}", 300, 0.0, 3.0);
            PtResolution_b.at(i_pt).at(i_eta) = h_jetRes_b_pt;

            //// Book histograms Mu response
            TH2F* h_muRes_tot_pt = new TH2F(GetHistName(i_pt, i_eta,"tot","Mu").c_str(), "p_T^{mu}/p_T^{gen}", 20, 0.0, 1.0, 20, 0.0, 1.0);
            MuRes_tot.at(i_pt).at(i_eta) = h_muRes_tot_pt;

            TH2F* h_muRes_LF_pt = new TH2F(GetHistName(i_pt, i_eta,"LF","Mu").c_str(), "p_T^{mu}/p_T^{gen}", 20, 0.0, 1.0, 20, 0.0, 1.0);
            MuRes_LF.at(i_pt).at(i_eta) = h_muRes_LF_pt;

            TH2F* h_muRes_HF_pt = new TH2F(GetHistName(i_pt, i_eta,"HF","Mu").c_str(), "p_T^{mu}/p_T^{gen}", 20, 0.0, 1.0, 20, 0.0, 1.0);
            MuRes_HF.at(i_pt).at(i_eta) = h_muRes_HF_pt;

            TH2F* h_muRes_nob_pt = new TH2F(GetHistName(i_pt, i_eta,"nob","Mu").c_str(), "p_T^{mu}/p_T^{gen}", 20, 0.0, 1.0, 20, 0.0, 1.0);
            MuRes_nob.at(i_pt).at(i_eta) = h_muRes_nob_pt;

            TH2F* h_muRes_b_pt = new TH2F(GetHistName(i_pt, i_eta,"b","Mu").c_str(), "p_T^{mu}/p_T^{gen}", 20, 0.0, 1.0, 20, 0.0, 1.0);
            MuRes_b.at(i_pt).at(i_eta) = h_muRes_b_pt;

            //// Book histograms Phi resolution
            TH1F* h_jetRes_tot_phi = new TH1F(GetHistName(i_pt, i_eta,"tot","Phi").c_str(), "#Delta#phi(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            PhiResolution_tot.at(i_pt).at(i_eta) = h_jetRes_tot_phi;

            TH1F* h_jetRes_LF_phi = new TH1F(GetHistName(i_pt, i_eta,"LF","Phi").c_str(), "#Delta#phi(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            PhiResolution_LF.at(i_pt).at(i_eta) = h_jetRes_LF_phi;

            TH1F* h_jetRes_HF_phi = new TH1F(GetHistName(i_pt, i_eta,"HF","Phi").c_str(), "#Delta#phi(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            PhiResolution_HF.at(i_pt).at(i_eta) = h_jetRes_HF_phi;

            TH1F* h_jetRes_nob_phi = new TH1F(GetHistName(i_pt, i_eta,"nob","Phi").c_str(), "#Delta#phi(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            PhiResolution_nob.at(i_pt).at(i_eta) = h_jetRes_nob_phi;

            TH1F* h_jetRes_b_phi = new TH1F(GetHistName(i_pt, i_eta,"b","Phi").c_str(), "#Delta#phi(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            PhiResolution_b.at(i_pt).at(i_eta) = h_jetRes_b_phi;

            //// Book histograms Eta resolution
            TH1F* h_jetRes_tot_eta = new TH1F(GetHistName(i_pt, i_eta,"tot","Eta").c_str(), "#Delta#eta(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            EtaResolution_tot.at(i_pt).at(i_eta) = h_jetRes_tot_eta;

            TH1F* h_jetRes_LF_eta = new TH1F(GetHistName(i_pt, i_eta,"LF","Eta").c_str(), "#Delta#eta(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            EtaResolution_LF.at(i_pt).at(i_eta) = h_jetRes_LF_eta;

            TH1F* h_jetRes_HF_eta = new TH1F(GetHistName(i_pt, i_eta,"HF","Eta").c_str(), "#Delta#eta(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            EtaResolution_HF.at(i_pt).at(i_eta) = h_jetRes_HF_eta;

            TH1F* h_jetRes_nob_eta = new TH1F(GetHistName(i_pt, i_eta,"nob","Eta").c_str(), "#Delta#eta(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            EtaResolution_nob.at(i_pt).at(i_eta) = h_jetRes_nob_eta;

            TH1F* h_jetRes_b_eta = new TH1F(GetHistName(i_pt, i_eta,"b","Eta").c_str(), "#Delta#eta(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            EtaResolution_b.at(i_pt).at(i_eta) = h_jetRes_b_eta;
        }
    }

    for (unsigned int i_eta = 0; i_eta < EtaBinEdges.size() - 1; ++i_eta) {

        //// Book histograms for jet counts (needed for reconstruction efficiency)
        TH1F* h_NjetReco_tot_pt = new TH1F(GetHistNameEta(i_eta,"h_NjetReco_tot").c_str(), "N(jet^{reco})", 100, 0., 500.);
        NReco_tot.push_back(h_NjetReco_tot_pt);

        TH1F* h_NjetReco_b_pt = new TH1F(GetHistNameEta(i_eta,"h_NjetReco_b").c_str(), "N(jet^{reco})", 100, 0., 500.);
        NReco_b.push_back(h_NjetReco_b_pt);

        TH1F* h_NjetReco_nob_pt = new TH1F(GetHistNameEta(i_eta,"h_NjetReco_nob").c_str(), "N(jet^{reco})", 100, 0., 500.);
        NReco_nob.push_back(h_NjetReco_nob_pt);

        TH1F* h_NjetGen_tot_pt = new TH1F(GetHistNameEta(i_eta,"h_NjetGen_tot").c_str(), "N(jet^{gen})", 100, 0., 500.);
        NGen_tot.push_back(h_NjetGen_tot_pt);

        TH1F* h_NjetGen_b_pt = new TH1F(GetHistNameEta(i_eta,"h_NjetGen_b").c_str(), "N(jet^{gen})", 100, 0., 500.);
        NGen_b.push_back(h_NjetGen_b_pt);

        TH1F* h_NjetGen_nob_pt = new TH1F(GetHistNameEta(i_eta,"h_NjetGen_nob").c_str(), "N(jet^{gen})", 100, 0., 500.);
        NGen_nob.push_back(h_NjetGen_nob_pt);
    }

    //NTotEvents = fChain->GetEntries();

    //// Not very elegant! TODO: Store this info in and read from file

    // [v3]
    AvailableEvents[361022] = 1993647;
    AvailableEvents[361023] = 7884494;
    AvailableEvents[361024] = 7889800;
    AvailableEvents[361025] = 7977600;
    AvailableEvents[361026] = 1893400;

}

void MyResolution::SlaveBegin(TTree * /*tree*/)
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

}

Bool_t MyResolution::Process(Long64_t entry)
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

    //std::cout << "Weight: " << *Weight << std::endl;
    double eventWeight = *Weight;
    eventWeight *= m_lumi / AvailableEvents[*DatasetID];

    eventWeight = 1.;

    std::vector<MyJet> genJets;
    std::vector<MyJet> recoJets;
    std::vector<MyMuon> recoMuons;

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
            return 1;

        }

        /*
                if (OR) {
                    //std::cout << "Reject event because of baseline muon!" << std::endl;
                    lv = true;
                    if (!cutFlowStudies) return 1;
                }
        */
        MyMuon muon(pt, eta, phi);
        muon.SetIsSignal(signal);
        muon.SetPassOR(OR);
        muon.SetIsBad(bad);
        muon.SetCharge(q);
        recoMuons.push_back(muon);
    }

    int NElectrons = ElePt->size();

    for (int i = 0; i < NElectrons; ++i) {
        bool OR = ElePassOR->at(i);
        bool signal = EleIsSignal->at(i);
        float pt = ElePt->at(i);
        float eta = EleEta->at(i);
        float phi = ElePhi->at(i);
        int q = EleCharge->at(i);

        if (signal) {
            //std::cout << "Reject event because of isolated electon!" << std::endl;
            return 1;
        }

        /*
                if (OR) {
                    //std::cout << "Reject event because of baseline electon!" << std::endl;
                    lv = true;
                    if (!cutFlowStudies) return 1;
                }
        */

    }


    // Loop over all genjets in this container
    for ( auto& genjet : genJets ) {

        // neglect genjets in very forward region
        if (fabs(genjet.Eta()) > 4.5) continue;

        // neglect very soft genjets
        if (fabs(genjet.Pt()) < 10.) continue;

        // check for no additional genJet activity
        bool noGenActivity = true;
        for ( auto& genjet2 : genJets ) {
            if (genjet2 == genjet) continue;
            double dR = genjet.DeltaR(genjet2);
            if (dR < m_VetoCone && (genjet2.Pt()/genjet.Pt()) > m_RelGenActivityVeto ) {
                noGenActivity = false;
            }
        }

        //if (!noGenActivity) continue; // continue with next genJet if another genJet is closeby

        // check for additional recoJet activity
        MyJet* matchedJet = 0;
        TLorentzVector addRecoActivity(0., 0., 0., 0.);
        double dRmin_matched = 999.;
        //std::cout << "GenJet (pt, eta, phi): " << genjet.Pt() << ", " << genjet.Eta() << ", " << genjet.Phi() << std::endl;
        for ( auto& jet : recoJets) {
            double dR = jet.DeltaR(genjet);
            if (dR > m_addActivityCone) continue;
            if (dR < dRmin_matched) {
                if (matchedJet) {
                    addRecoActivity += *matchedJet;
                    //std::cout << "Add previous matched jet (pt, eta, phi): " << matchedJet->Pt() << ", " << matchedJet->Eta() << ", " << matchedJet->Phi() << std::endl;
                    dRmin_matched = dR;
                    matchedJet = &jet;
                } else {
                    dRmin_matched = dR;
                    matchedJet = &jet;
                }
            } else {
                addRecoActivity += jet;
                //std::cout << "Add Jet (pt, eta, phi): " << jet.Pt() << ", " << jet.Eta() << ", " << jet.Phi() << std::endl;
            }
        } // end for loop over jets

        TLorentzVector JetMu(0.,0.,0.,0.);
        for ( auto& mu : recoMuons) {
            if (mu.DeltaR(genjet) < 0.4) JetMu += mu;
        }

        bool noRecoActivity = true;
        if (matchedJet) {
            if ((addRecoActivity.Pt()/matchedJet->Pt())>m_RelRecoActivityVeto) noRecoActivity = false;
            //if (matchedJet->Pt() < (genjet.Pt() - 100.) && matchedJet->Pt() < genjet.Pt()/2. && addRecoActivity.Pt() < 50. ) {
            //std::cout << "GenJet (pt, eta, phi): " << genjet.Pt() << ", " << genjet.Eta() << ", " << genjet.Phi() << std::endl;
            //std::cout << "MatchedJet (pt, eta, phi): " << matchedJet->Pt() << ", " << matchedJet->Eta() << ", " << matchedJet->Phi() << std::endl;
            //std::cout << "dR: " << dRmin_matched << std::endl;
            //std::cout << "addRecoActivity (pt, eta, phi): " << addRecoActivity.Pt() << ", " << addRecoActivity.Eta() << ", " << addRecoActivity.Phi() << std::endl;
            //}
        }

        if (noGenActivity) {
            int ii_eta = GetEtaBin(genjet.Eta());
            NGen_tot.at(ii_eta)->Fill(genjet.Pt(), eventWeight);
            if (genjet.IsB()) {
                NGen_b.at(ii_eta)->Fill(genjet.Pt(), eventWeight);
            } else {
                NGen_nob.at(ii_eta)->Fill(genjet.Pt(), eventWeight);
            }
        }

        if (noGenActivity && dRmin_matched < 0.2) {
            int ii_eta = GetEtaBin(genjet.Eta());
            NReco_tot.at(ii_eta)->Fill(genjet.Pt(), eventWeight);
            if (genjet.IsB()) {
                NReco_b.at(ii_eta)->Fill(genjet.Pt(), eventWeight);
            } else {
                NReco_nob.at(ii_eta)->Fill(genjet.Pt(), eventWeight);
            }
        }

        if (dRmin_matched < m_MatchingCone && noRecoActivity && noGenActivity) {
            int i_pt = GetPtBin(genjet.E());
            int i_eta = GetEtaBin(genjet.Eta());
            PtResolution_tot.at(i_pt).at(i_eta)->Fill((*matchedJet+addRecoActivity).Pt()/genjet.Pt(), eventWeight);
            PhiResolution_tot.at(i_pt).at(i_eta)->Fill(matchedJet->DeltaPhi(genjet), eventWeight);
            EtaResolution_tot.at(i_pt).at(i_eta)->Fill(matchedJet->Eta()-genjet.Eta(), eventWeight);
            MuRes_tot.at(i_pt).at(i_eta)->Fill((*matchedJet+addRecoActivity).Pt()/genjet.Pt(), JetMu.Pt()/genjet.Pt(), eventWeight);
            if (matchedJet->IsB()) {
                PtResolution_b.at(i_pt).at(i_eta)->Fill((*matchedJet+addRecoActivity).Pt()/genjet.Pt(), eventWeight);
                PhiResolution_b.at(i_pt).at(i_eta)->Fill(matchedJet->DeltaPhi(genjet), eventWeight);
                EtaResolution_b.at(i_pt).at(i_eta)->Fill(matchedJet->Eta()-genjet.Eta(), eventWeight);
                MuRes_b.at(i_pt).at(i_eta)->Fill((*matchedJet+addRecoActivity).Pt()/genjet.Pt(), JetMu.Pt()/genjet.Pt(), eventWeight);
            } else {
                PtResolution_nob.at(i_pt).at(i_eta)->Fill((*matchedJet+addRecoActivity).Pt()/genjet.Pt(), eventWeight);
                PhiResolution_nob.at(i_pt).at(i_eta)->Fill(matchedJet->DeltaPhi(genjet), eventWeight);
                EtaResolution_nob.at(i_pt).at(i_eta)->Fill(matchedJet->Eta()-genjet.Eta(), eventWeight);
                MuRes_nob.at(i_pt).at(i_eta)->Fill((*matchedJet+addRecoActivity).Pt()/genjet.Pt(), JetMu.Pt()/genjet.Pt(), eventWeight);
            }
            if (genjet.IsB()) {
                PtResolution_HF.at(i_pt).at(i_eta)->Fill((*matchedJet+addRecoActivity).Pt()/genjet.Pt(), eventWeight);
                PhiResolution_HF.at(i_pt).at(i_eta)->Fill(matchedJet->DeltaPhi(genjet), eventWeight);
                EtaResolution_HF.at(i_pt).at(i_eta)->Fill(matchedJet->Eta()-genjet.Eta(), eventWeight);
                MuRes_HF.at(i_pt).at(i_eta)->Fill((*matchedJet+addRecoActivity).Pt()/genjet.Pt(), JetMu.Pt()/genjet.Pt(), eventWeight);
            } else {
                PtResolution_LF.at(i_pt).at(i_eta)->Fill((*matchedJet+addRecoActivity).Pt()/genjet.Pt(), eventWeight);
                PhiResolution_LF.at(i_pt).at(i_eta)->Fill(matchedJet->DeltaPhi(genjet), eventWeight);
                EtaResolution_LF.at(i_pt).at(i_eta)->Fill(matchedJet->Eta()-genjet.Eta(), eventWeight);
                MuRes_LF.at(i_pt).at(i_eta)->Fill((*matchedJet+addRecoActivity).Pt()/genjet.Pt(), JetMu.Pt()/genjet.Pt(), eventWeight);
            }
        }
    }

    return kTRUE;
}

void MyResolution::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.

    for (unsigned int i_pt = 0; i_pt < PtBinEdges.size() - 1; ++i_pt) {
        for (unsigned int i_eta = 0; i_eta < EtaBinEdges.size() - 1; ++i_eta) {
            PtResolution_tot.at(i_pt).at(i_eta)->Write();
            PtResolution_LF.at(i_pt).at(i_eta)->Write();
            PtResolution_HF.at(i_pt).at(i_eta)->Write();
            PtResolution_nob.at(i_pt).at(i_eta)->Write();
            PtResolution_b.at(i_pt).at(i_eta)->Write();

            PhiResolution_tot.at(i_pt).at(i_eta)->Write();
            PhiResolution_LF.at(i_pt).at(i_eta)->Write();
            PhiResolution_HF.at(i_pt).at(i_eta)->Write();
            PhiResolution_nob.at(i_pt).at(i_eta)->Write();
            PhiResolution_b.at(i_pt).at(i_eta)->Write();

            EtaResolution_tot.at(i_pt).at(i_eta)->Write();
            EtaResolution_LF.at(i_pt).at(i_eta)->Write();
            EtaResolution_HF.at(i_pt).at(i_eta)->Write();
            EtaResolution_nob.at(i_pt).at(i_eta)->Write();
            EtaResolution_b.at(i_pt).at(i_eta)->Write();

            MuRes_tot.at(i_pt).at(i_eta)->Write();
            MuRes_LF.at(i_pt).at(i_eta)->Write();
            MuRes_HF.at(i_pt).at(i_eta)->Write();
            MuRes_nob.at(i_pt).at(i_eta)->Write();
            MuRes_b.at(i_pt).at(i_eta)->Write();
        }
    }

    for (unsigned int i_eta = 0; i_eta < EtaBinEdges.size() - 1; ++i_eta) {
        NReco_tot.at(i_eta)->Write();
        NReco_b.at(i_eta)->Write();
        NReco_nob.at(i_eta)->Write();
        NGen_tot.at(i_eta)->Write();
        NGen_b.at(i_eta)->Write();
        NGen_nob.at(i_eta)->Write();
    }

    outputfile->Close();

    for (std::map<UInt_t, UInt_t>::iterator it=ProcessedEvents.begin(); it != ProcessedEvents.end(); ++it) {
        std::cout << it->first << " " << it->second << std::endl;
    }
}

void MyResolution::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

}

std::string MyResolution::GetHistName(unsigned int i_pt, unsigned int i_eta, std::string s1, std::string s2) {
    std::string hname = "h_"+s1+"_JetAll_Res"+s2+"_Pt";
    hname += std::to_string(i_pt);
    hname += "_Eta";
    hname += std::to_string(i_eta);
    return hname;
}



std::string MyResolution::GetHistNameEta(unsigned int i_eta, std::string s1) {
    std::string hname = s1+"_Eta";
    hname += std::to_string(i_eta);
    return hname;
}



unsigned int MyResolution::GetPtBin(double pt) {
    int i_pt = -1;
    for (std::vector<double>::const_iterator it = PtBinEdges.begin(); it != PtBinEdges.end(); ++it) {
        if ((*it) > pt)
            break;
        ++i_pt;
    }
    if (i_pt < 0)
        i_pt = 0;
    if (i_pt > (int) PtBinEdges.size() - 2)
        i_pt = (int) PtBinEdges.size() - 2;

    return i_pt;
}



unsigned int MyResolution::GetEtaBin(double eta) {
    int i_eta = -1;
    for (std::vector<double>::const_iterator it = EtaBinEdges.begin(); it != EtaBinEdges.end(); ++it) {
        if ((*it) > fabs(eta))
            break;
        ++i_eta;
    }
    if (i_eta < 0)
        i_eta = 0;
    if (i_eta > (int) EtaBinEdges.size() - 2)
        i_eta = (int) EtaBinEdges.size() - 2;
    return i_eta;
}



void MyResolution::ResizeHistoVector(std::vector<std::vector<TH1F*> > &histoVector) {

    histoVector.resize(PtBinEdges.size() - 1);
    for (std::vector<std::vector<TH1F*> >::iterator it = histoVector.begin(); it != histoVector.end(); ++it) {
        it->resize(PtBinEdges.size() - 1);
    }
}

void MyResolution::ResizeHisto2Vector(std::vector<std::vector<TH2F*> > &histoVector) {

    histoVector.resize(PtBinEdges.size() - 1);
    for (std::vector<std::vector<TH2F*> >::iterator it = histoVector.begin(); it != histoVector.end(); ++it) {
        it->resize(PtBinEdges.size() - 1);
    }
}

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
#include "MyMETStudies.h"

void MyMETStudies::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

    outputfile = new TFile("METStudiesOutput_0L.root","RECREATE");

    m_MatchingCone = 0.1;
    m_jvtcut = 0.59;
    m_lumi = 10000.;

    ////Book histograms

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

    h_MHT = new TH1F("h_MHT", "h_MHT", 100, 0., 1000.);
    h_MHT->Sumw2();
    histos_1D.push_back(h_MHT);

    h_NoJVTMHT = new TH1F("h_NoJVTMHT", "h_NoJVTMHT", 100, 0., 1000.);
    h_NoJVTMHT->Sumw2();
    histos_1D.push_back(h_NoJVTMHT);

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

    h_MHTvsNoJVTMHT = new TH2F("h_MHTvsNoJVTMHT", "h_MHTvsNoJVTMHT", 100, 0., 1000., 100, 0., 1000.);
    h_MHTvsNoJVTMHT->Sumw2();
    histos_2D.push_back(h_MHTvsNoJVTMHT);

    h_METvsMHT = new TH2F("h_METvsMHT", "h_METvsMHT", 100, 0., 1000., 100, 0., 1000.);
    h_METvsMHT->Sumw2();
    histos_2D.push_back(h_METvsMHT);

    h_METvsNoJVTMHT = new TH2F("h_METvsNoJVTMHT", "h_METvsNoJVTMHT", 100, 0., 1000., 100, 0., 1000.);
    h_METvsNoJVTMHT->Sumw2();
    histos_2D.push_back(h_METvsNoJVTMHT);

    h_METvsMETreplaced = new TH2F("h_METvsMETreplaced", "h_METvsMETreplaced", 100, 0., 1000., 100, 0., 1000.);
    h_METvsMETreplaced->Sumw2();
    histos_2D.push_back(h_METvsMETreplaced);

    h_METvsElePt = new TH2F("h_METvsElePt", "h_METvsElePt", 100, 0., 1000., 100, 0., 1000.);
    h_METvsElePt->Sumw2();
    histos_2D.push_back(h_METvsElePt);

    h_METvsMuPt = new TH2F("h_METvsMuPt", "h_METvsMuPt", 100, 0., 1000., 100, 0., 1000.);
    h_METvsMuPt->Sumw2();
    histos_2D.push_back(h_METvsMuPt);

    h_METvsPhoPt = new TH2F("h_METvsPhoPt", "h_METvsPhoPt", 100, 0., 1000., 100, 0., 1000.);
    h_METvsPhoPt->Sumw2();
    histos_2D.push_back(h_METvsPhoPt);

    h_METvsLepPt = new TH2F("h_METvsLepPt", "h_METvsLepPt", 100, 0., 1000., 100, 0., 1000.);
    h_METvsLepPt->Sumw2();
    histos_2D.push_back(h_METvsLepPt);

    h_METvsMETjet = new TH2F("h_METvsMETjet", "h_METvsMETjet", 100, 0., 1000., 100, 0., 1000.);
    h_METvsMETjet->Sumw2();
    histos_2D.push_back(h_METvsMETjet);

    h_METvsMETele = new TH2F("h_METvsMETele", "h_METvsMETele", 100, 0., 1000., 100, 0., 1000.);
    h_METvsMETele->Sumw2();
    histos_2D.push_back(h_METvsMETele);

    h_METvsMETmu = new TH2F("h_METvsMETmu", "h_METvsMETmu", 100, 0., 1000., 100, 0., 1000.);
    h_METvsMETmu->Sumw2();
    histos_2D.push_back(h_METvsMETmu);

    h_METvsMETgamma = new TH2F("h_METvsMETgamma", "h_METvsMETgamma", 100, 0., 1000., 100, 0., 1000.);
    h_METvsMETgamma->Sumw2();
    histos_2D.push_back(h_METvsMETgamma);

    h_METvsMETtrack = new TH2F("h_METvsMETtrack", "h_METvsMETtrack", 100, 0., 1000., 100, 0., 1000.);
    h_METvsMETtrack->Sumw2();
    histos_2D.push_back(h_METvsMETtrack);

    //NTotEvents = fChain->GetEntries();

    //// Not very elegant! TODO: Store this info in and read from file

    // NTulpleMaker_mc_default_v2
    AvailableEvents[361022] = 1993647;
    AvailableEvents[361023] = 7634495;
    AvailableEvents[361024] = 7979800;
    AvailableEvents[361025] = 7907600;
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
    if (NEvents%1000 == 0) std::cout << NEvents << " processed!" << std::endl;

    if (ProcessedEvents.find(*DatasetID) == ProcessedEvents.end()) {
        ProcessedEvents[*DatasetID] = 1;
    } else {
        ProcessedEvents[*DatasetID] += 1;
    }

    if ( *DatasetID == 361026 ) return 0;

    //std::cout << "Weight: " << *Weight << std::endl;
    double eventWeight = *Weight;
    eventWeight *= m_lumi / AvailableEvents[*DatasetID];

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

    bool lv = false;

    int NMuons = MuonPt->size();
    //std::cout << "NMuons: " << NMuons << std::endl;
    //if ( NMuons > 0 ) return 1;

    for (int i = 0; i < NMuons; ++i) {
        float pt = MuonPt->at(i);
        float eta = MuonEta->at(i);
        float phi = MuonPhi->at(i);
        if (MuonIsBad->at(i)) {
            //std::cout << "Reject event because of bad muon!" << std::endl;
            return 1;
        }
        //if (MuonIsSignal->at(i) && MuonIsIso->at(i)) {
        if (MuonIsSignal->at(i)) {
            //std::cout << "Reject event because of isolated muon!" << std::endl;
            lv = true;
            //return 1;

        }
        //if (MuonIsSignal->at(i)) {
        TLorentzVector muon(0.,0.,0.,0.);
        muon.SetPtEtaPhiM(pt, eta, phi, 0.1057);
        recoMuons.push_back(muon);
        //}
    }
    GreaterByPt<TLorentzVector> ptComparator2_;
    std::sort(recoMuons.begin(), recoMuons.end(), ptComparator2_);
    double PtMu = 0;
    if (NMuons > 0) PtMu = recoMuons.at(0).Pt();

    int NElectrons = ElePt->size();
    //std::cout << "NElectrons: " << NElectrons << std::endl;
    if ( NElectrons > 0 ) {
        //std::cout << "Reject event because of isolated electon!" << std::endl;
        lv = true;
        //return 1;
    }
    for (int i = 0; i < NElectrons; ++i) {
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

    Ntot += 1;
    Ntot_w += eventWeight;

	MyJet* firstJet = 0;
	MyJet* secondJet = 0;
	MyJet* thirdJet = 0;
	MyJet* fourthJet = 0;

    if (*xe90triggered || *xe110triggered) {
        Ntrig += 1;
        Ntrig_w += eventWeight;

        for ( auto& jet : recoJets) {
            if (jet.IsGood() && jet.Pt() > 25. && (jet.Pt() > 60. || jet.IsNoPU(m_jvtcut) || fabs(jet.Eta()) > 2.4)) {
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

        if (secondJet) {
            N2jets += 1;
            N2jets_w += eventWeight;
            if (firstJet->Pt()>80. && secondJet->Pt() > 50.) {
                NjetPt += 1;
                NjetPt_w += eventWeight;
                if (thirdJet && !fourthJet) {
                    if (thirdJet->Pt() > 25. && thirdJet->Pt() < 50.) {
                        NTJV += 1;
                        NTJV_w += eventWeight;
                        double dPhi = fabs(firstJet->DeltaPhi(*secondJet));
                        if (dPhi < 2.5 && dPhi > 1.8) {
                            NdPhi += 1;
                            NdPhi_w += eventWeight;
                            double dEta = fabs(firstJet->Eta() - secondJet->Eta());
                            if (dEta > 3.0) {
                                NdEta += 1;
                                NdEta_w += eventWeight;
								std::cout << "1: " << firstJet->Pt() << ", " << firstJet->Eta() << ", " << firstJet->Phi()<< std::endl;
								std::cout << "2: " << secondJet->Pt() << ", " << secondJet->Eta() << ", " << secondJet->Phi()<< std::endl;
								std::cout << "3: " << thirdJet->Pt() << ", " << thirdJet->Eta() << ", " << thirdJet->Phi()<< std::endl;
                                if ( (firstJet->Eta() * secondJet->Eta()) < 0 ) {
                                    Nhemi += 1;
                                    Nhemi_w += eventWeight;
                                    if ((*firstJet + *secondJet).Pt() > 150.) {
                                        NpTjj += 1;
                                        NpTjj_w += eventWeight;
                                        if ((*firstJet + *secondJet).M() > 1000.) {
                                            Nmjj += 1;
                                            Nmjj_w += eventWeight;
                                            if (!lv) {
                                                Nlv += 1;
                                                Nlv_w += eventWeight;
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
    }

    return 0;

    // Loop over all genjets in this container
    TLorentzVector TruthMHT(0.,0.,0.,0.);
    for ( auto& genjet : genJets ) {
        if (genjet.Pt() > 20. && fabs(genjet.Eta()) < 2.4) {
            TruthMHT -= genjet;
        }
    }
    h_TruthMHT->Fill(TruthMHT.Pt(), eventWeight);

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
    TLorentzVector NoJVTMHT(0.,0.,0.,0.);
    TLorentzVector TruthMHTmatched(0.,0.,0.,0.);

    int JetCount = 0;
    int GoodJetCount = 0;
    //MyJet* firstJet = 0;
    firstJet = 0;
    //MyJet* secondJet = 0;
    secondJet = 0;

    for ( auto& jet : recoJets) {

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

        NoJVTMHT-=jet;
        if (jet.IsGood() && (jet.Pt() > 60. || jet.IsNoPU(m_jvtcut) || fabs(jet.Eta()) > 2.4)) {
            if (GoodJetCount == 0) firstJet = &jet;
            if (GoodJetCount == 1) secondJet = &jet;
            MHT-=jet;
            MyJet* matchedJet = 0;
            double dRmin = 999.;
            for ( auto& genjet : genJets ) {
                if (genjet.DeltaR(jet) < dRmin) {
                    matchedJet = &genjet;
                    dRmin = genjet.DeltaR(jet);
                }
            }
            if (dRmin < 0.1) {
                TruthMHTmatched -= *matchedJet;
            } else {
                TruthMHTmatched -= jet;
            }
            ++GoodJetCount;
        }

        ++JetCount;
    }

    if (GoodJetCount > 1) {

        if (firstJet->Pt() > 80. && secondJet->Pt() > 50. && fabs(firstJet->Eta() - secondJet->Eta()) > 3.5 && (firstJet->Eta()*secondJet->Eta()) < 0. ) {

            TLorentzVector METreplaced(0.,0.,0.,0.);
            METreplaced = MET - MHT + TruthMHTmatched;
            h_METreplaced->Fill(METreplaced.Pt(), eventWeight);

            h_MHT->Fill(MHT.Pt(), eventWeight);
            h_NoJVTMHT->Fill(NoJVTMHT.Pt(), eventWeight);


            h_MHTvsNoJVTMHT->Fill(MHT.Pt(), NoJVTMHT.Pt(), eventWeight);

            h_METvsNoJVTMHT->Fill(MET.Pt(), NoJVTMHT.Pt(), eventWeight);

            h_METvsMHT->Fill(MET.Pt(), MHT.Pt(), eventWeight);

            h_METvsMETreplaced->Fill(MET.Pt(), METreplaced.Pt(), eventWeight);

            h_METvsElePt->Fill(MET.Pt(), PtEle, eventWeight);

            h_METvsMuPt->Fill(MET.Pt(), PtMu, eventWeight);

            h_METvsPhoPt->Fill(MET.Pt(), PtPho, eventWeight);

            h_METvsLepPt->Fill(MET.Pt(), (PtEle>PtMu ? PtEle : PtMu), eventWeight);

            h_METvsMETjet->Fill(MET.Pt(), METjet.Pt(), eventWeight);

            h_METvsMETele->Fill(MET.Pt(), METele.Pt(), eventWeight);

            h_METvsMETmu->Fill(MET.Pt(), METmu.Pt(), eventWeight);

            h_METvsMETgamma->Fill(MET.Pt(), METgamma.Pt(), eventWeight);

            h_METvsMETtrack->Fill(MET.Pt(), METtrack.Pt(), eventWeight);

        }

    }

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

    std::cout << "Ntot : " << Ntot << std::endl;
    std::cout << "Ntrig : " << Ntrig << std::endl;
    std::cout << "N2jets : " << N2jets << std::endl;
    std::cout << "NjetPt : " << NjetPt << std::endl;
    std::cout << "NTJV : " << NTJV << std::endl;
    std::cout << "NdPhi : " << NdPhi << std::endl;
    std::cout << "NdEta : " << NdEta << std::endl;
    std::cout << "Nhemi : " << Nhemi << std::endl;
    std::cout << "NpTjj : " << NpTjj << std::endl;
    std::cout << "Nmjj : " << Nmjj << std::endl;
    std::cout << "Nlv : " << Nlv << std::endl;

    std::cout << "Ntot_w : " << Ntot_w << std::endl;
    std::cout << "Ntrig_w : " << Ntrig_w << std::endl;
    std::cout << "N2jets_w : " << N2jets_w << std::endl;
    std::cout << "NjetPt_w : " << NjetPt_w << std::endl;
    std::cout << "NTJV_w : " << NTJV_w << std::endl;
    std::cout << "NdPhi_w : " << NdPhi_w << std::endl;
    std::cout << "NdEta_w : " << NdEta_w << std::endl;
    std::cout << "Nhemi_w : " << Nhemi_w << std::endl;
    std::cout << "NpTjj_w : " << NpTjj_w << std::endl;
    std::cout << "Nmjj_w : " << Nmjj_w << std::endl;
    std::cout << "Nlv_w : " << Nlv_w << std::endl;

}

void MyMETStudies::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

}

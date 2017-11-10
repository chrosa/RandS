#define RandS_cxx
// The class definition in RandS.h has been generated automatically
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
// root> T->Process("RandS.C")
// root> T->Process("RandS.C","some options")
// root> T->Process("RandS.C+")
//

#include "TKinFitter/TKinFitter.h"
#include "TKinFitter/TFitParticleEtEtaPhi.h"
#include "TKinFitter/TFitConstraintEp.h"
#include "RandS.h"

void RandS::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

    isMC_ = false;
    jvtcut_= 0.59; //// 0.59 (medium)
    lumi_ = 36100.;
    smearingfile_ = "/afs/desy.de/user/c/csander/xxl-af-cms/testarea/2.4.8/MyAnalysis/util/resolutions_GenJetMuNu_RecoNoMu_E_OR_v2.root";
    //inputhistPtHF_ = "h_tot_JetAll_ResPt";
    //inputhistEtaHF_ = "h_tot_JetAll_ResEta";
    //inputhistPhiHF_ = "h_tot_JetAll_ResPhi";
    //inputhistPtLF_ = "h_tot_JetAll_ResPt";
    //inputhistEtaLF_ = "h_tot_JetAll_ResEta";
    //inputhistPhiLF_ = "h_tot_JetAll_ResPhi";
    inputhistPtHF_ = "h_b_JetAll_ResPt";
    inputhistEtaHF_ = "h_b_JetAll_ResEta";
    inputhistPhiHF_ = "h_b_JetAll_ResPhi";
    inputhistPtLF_ = "h_nob_JetAll_ResPt";
    inputhistEtaLF_ = "h_nob_JetAll_ResEta";
    inputhistPhiLF_ = "h_nob_JetAll_ResPhi";
    //inputhistPtHF_ = "h_HF_JetAll_ResPt";
    //inputhistEtaHF_ = "h_HF_JetAll_ResEta";
    //inputhistPhiHF_ = "h_HF_JetAll_ResPhi";
    //inputhistPtLF_ = "h_LF_JetAll_ResPt";
    //inputhistEtaLF_ = "h_LF_JetAll_ResEta";
    //inputhistPhiLF_ = "h_LF_JetAll_ResPhi";
    inputhistResMuHF_ = "h_b_JetAll_ResMu";
    inputhistResMuLF_ = "h_nob_JetAll_ResMu";
    // Reminder here E is used instead pT for binning (variabel name not changed yet, maybe later)
    //PtBinEdges_ = {0,20,30,40,50,70,90,110,140,170,200,240,280,330,380,440,500,570,640,720,810,910,1020,1140,1270,1410,1660,1820,1990,2170,2360,2560,3000,9999};
    PtBinEdges_ = {0,10,20,30,40,50,70,100,140,190,250,320,400,490,590,700,820,950,1090,1240,1400,1570,1750,1940,2140,2350,2600,3000};
    //PtBinEdges_ = {0,20,40,70,140,250,400,590,820,1090,1400,1750,2140,2600};
    EtaBinEdges_ = {0.0,0.7,1.3,1.8,2.5,3.2,5.0};
    ResBinEdges_ = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
    rebalancedJetPt_ = 20.;
    rebalancedNJet_ = 99; // Number of leading jets used in rebalancing
    rebalanceMode_ = "METsoft"; // "METsoft" (should be best), "MHTall", "MHThigh"
    smearCollection_ = "Reco"; // "Gen" for truth jet smearing only; "Reco" for full R+S
    smearedJetPt_ = 15.;
    smearedJetEta_ = 5.0;
    smearedNJet_ = 99; // Number of leading jets used in smearing
    doSmearing_ = true; // only "false" for test purposes
    doJVT_ = true; // if "false" decide on used jets prior to rebalancing
    fixJVTjets_ = false; // if "true" jvt jets will not be changed in rebalancing process
    JVTeta_ = 5.; // default is 2.4
    JetEffEmulation_ = false;
    doMETmu_ = true;
    PtBinEdges_scaling_ = {0.,3000.};
    EtaBinEdges_scaling_ = {0.,5.};
    AdditionalSmearing_ = {1.0};
    LowerTailScaling_ = {1.0};
    UpperTailScaling_ = {1.0};
    AdditionalSmearing_variation_ = 1.0;
    LowerTailScaling_variation_ = 1.0;
    UpperTailScaling_variation_ = 1.0;
    absoluteTailScaling_ = false; // depends on definition of tails scail factors ("true" for M. Schroeder's definition)
    A0RMS_ = 2.5;
    A1RMS_ = 10.;
    probExtreme_ = 0.; // possibilty to emulate total jet loss
    uncertaintyName_ = "";
    useMETsoftResolution_ = true; // Smear METsoft component prior to rebalancing
    useTrueMETsoftForRebalance_ = false; // only possible on simulated events, for testing purposes (best performance)
    useTriggerTurnOn_ = true; //save event weight from trigger turn on
    METsoftResolutionFile_ = "/afs/desy.de/user/c/csander/xxl-af-cms/testarea/2.4.8/MyAnalysis/util/METsoft_resolutions_new.root";
    triggerTurnOnFile_ = "/afs/desy.de/user/c/csander/xxl-af-cms/testarea/2.4.8/MyAnalysis/util/TriggerStudiesOutput_data.root";
    controlPlots_ = false;
    debug_ = 0;
    outputfile_ = "RandS.root";
    cleverPrescaleTreating_ = true; // "true", to get better statistical  precision for high weight seed events
    maxCleverWeight_ = 5000; // the larger, the better (but also much slower), not greater than O(1000)
    HTSeedMin_ = 0.;
    HTSeedMax_ = 99999.;
    MHTSeedMax_ = 99999.;
    NJetsSeedMin_ = 2;
    NJetsSeedMax_ = 99;
    NJetsStored_ = 4;
    Ntries_ = 20;
    NJetsSave_ = 0;
    MjjSave_ = 0.;
    HTSave_ = 0;
    METSave_ = 100.;
    MHTSave_ = 9999.;
    MHTjjSave_ = 150.;
    BJetsPt_ = 25.;
    BJetsEta_ = 2.4;
    JetsPt_ = 25.;
    JetsEta_ = 5.0;
    JetsHTPt_ = 25.;
    JetsHTEta_ = 5.0;
    JetsMHTPt_ = 25.;
    JetsMHTEta_ = 5.0;
    JetDeltaMin_ = {0.5,0.5,0.3};
    MjjFirstPt_ = 80.;
    MjjSecondPt_ = 50.;
    dPhiSave_ = 2.7;
    dEtaSave_ = 2.5;
    //dEtaSave_ = 0.;
    jet3PtSave_ = 50.;

    smearFunc_ = new SmearFunction(smearingfile_,
                                   inputhistPtHF_,inputhistEtaHF_,inputhistPhiHF_,inputhistResMuHF_,
                                   inputhistPtLF_,inputhistEtaLF_,inputhistPhiLF_,inputhistResMuLF_,
                                   PtBinEdges_,EtaBinEdges_,ResBinEdges_,
                                   PtBinEdges_scaling_,EtaBinEdges_scaling_,
                                   AdditionalSmearing_,LowerTailScaling_,UpperTailScaling_,AdditionalSmearing_variation_,LowerTailScaling_variation_,UpperTailScaling_variation_,absoluteTailScaling_,
                                   A0RMS_,A1RMS_,probExtreme_
                                  );

    if( useMETsoftResolution_ ) {
        TFile *f_METsoft = new TFile(METsoftResolutionFile_.c_str(), "READ", "", 0);
        h_METsoft_Pt = (TH2F*) f_METsoft->FindObjectAny("h_MHTtruerebPt_vs_MHTnoJVT");
        h_METsoft_Pt_px.resize(h_METsoft_Pt->GetYaxis()->GetNbins());
        for (int jj = 1; jj <= h_METsoft_Pt->GetYaxis()->GetNbins(); ++jj) {
            TH1D* tmp = new TH1D(*h_METsoft_Pt->ProjectionX("px", jj, jj));
            h_METsoft_Pt_px.at(jj-1) = tmp;
        }
        h_METsoft_Phi = (TH2F*) f_METsoft->FindObjectAny("h_MHTtruerebPhiRes_vs_MHTtruereb");
        h_METsoft_Phi_px.resize(h_METsoft_Phi->GetYaxis()->GetNbins());
        for (int jj = 1; jj <= h_METsoft_Phi->GetYaxis()->GetNbins(); ++jj) {
            TH1D* tmp = new TH1D(*h_METsoft_Phi->ProjectionX("px", jj, jj));
            h_METsoft_Phi_px.at(jj-1) = tmp;
        }
    }

    if( useTriggerTurnOn_ ) {
        TFile *f_triggerTurnOn = new TFile(triggerTurnOnFile_.c_str(), "READ", "", 0);
        h_MHTvsHT_all  =  (TH2F*) f_triggerTurnOn->FindObjectAny("h_MHT2jetvsHT_all");
        h_MHTvsHT_triggered  =  (TH2F*) f_triggerTurnOn->FindObjectAny("h_MHT2jetvsHT_triggered");
        h_MHTvsHT_eff = new TH2F(*h_MHTvsHT_triggered);
        //cout << "Pointers (all, triggered, eff): " << h_MHTvsHT_all << ", " << h_MHTvsHT_triggered << ", " << h_MHTvsHT_eff << endl;
        h_MHTvsHT_eff->Divide(h_MHTvsHT_all);
    }

    // define output tree
    cout << "outputfile_: " << outputfile_ << endl;
    outputFile = new TFile(outputfile_.c_str(),"RECREATE");
    PredictionTree = new TTree("PredictionTree", "PredictionTree");
    PredictionTree->SetDirectory(outputFile);
    PredictionTree->SetAutoSave(10000000000);
    PredictionTree->SetAutoFlush(100000000);

    //cout << PredictionTree << endl;

    // set branches for output tree
    //PredictionTree->Branch("NVtx", &vtxN);
    PredictionTree->Branch("Ntries",&Ntries_pred);
    PredictionTree->Branch("NJets",&Njets_pred);
    PredictionTree->Branch("BTags",&BTags_pred);
    PredictionTree->Branch("Weight",&weight);
    PredictionTree->Branch("TriggerWeight",&triggerWeight);
    PredictionTree->Branch("HT", &HT_pred);
    PredictionTree->Branch("MHT", &MHT_pred);
    PredictionTree->Branch("MHTphi", &MHTphi_pred);
    PredictionTree->Branch("MET", &MET_pred);
    PredictionTree->Branch("METphi", &METphi_pred);
    PredictionTree->Branch("METsig", &METsig_seed);
    PredictionTree->Branch("MHTsig", &MHTsig_seed);
    PredictionTree->Branch("METsoft", &METsoft_seed);
    PredictionTree->Branch("JetPt", "std::vector<Float_t>", &JetPt_pred);
    PredictionTree->Branch("JetEta", "std::vector<Float_t>", &JetEta_pred);
    PredictionTree->Branch("JetPhi", "std::vector<Float_t>", &JetPhi_pred);
    PredictionTree->Branch("JetM", "std::vector<Float_t>", &JetM_pred);
    PredictionTree->Branch("DeltaPhi", "std::vector<Float_t>", &DeltaPhi_pred);

    //NTotEvents = fChain->GetEntries();

    // Different seed per initialization
    gRandom->SetSeed(0);
    rand_ = new TRandom3(0);
    gErrorIgnoreLevel = kError;

    //// Not very elegant! TODO: Store this info in and read from file

    // [v1]
    AvailableEvents[361022] = 1993647;
    AvailableEvents[361023] = 7724495;
    AvailableEvents[361024] = 7890000;
    AvailableEvents[361025] = 7977600;
    AvailableEvents[361026] = 1833400;

    AvailableEvents[364184] = 22856611;
    AvailableEvents[364185] = 9379984;
    AvailableEvents[364186] = 16307429;
    AvailableEvents[364187] = 14554282;
    AvailableEvents[364188] = 9845431;
    AvailableEvents[364189] = 8995506;
    AvailableEvents[364190] = 9716288;
    AvailableEvents[364191] = 7364236;
    AvailableEvents[364192] = 24232000;


    if (controlPlots_) {

        h_RebRes_genPt_eta0 = new TH2F("h_RebRes_genPt_eta0","h_RebRes_genPt_eta0", 50, 0., 500., 100, 0., 2.);
        h_RebRes_genPt_eta1 = new TH2F("h_RebRes_genPt_eta1","h_RebRes_genPt_eta1", 50, 0., 500., 100, 0., 2.);
        h_RebRes_genPt_eta2 = new TH2F("h_RebRes_genPt_eta2","h_RebRes_genPt_eta2", 50, 0., 500., 100, 0., 2.);

        h_RebRes_genPt_jet1 = new TH2F("h_RebRes_genPt_jet1","h_RebRes_genPt_jet1", 50, 0., 500., 100, 0., 2.);
        h_RebRes_genPt_jet2 = new TH2F("h_RebRes_genPt_jet2","h_RebRes_genPt_jet2", 50, 0., 500., 100, 0., 2.);
        h_RebRes_genPt_jet3 = new TH2F("h_RebRes_genPt_jet3","h_RebRes_genPt_jet3", 50, 0., 500., 100, 0., 2.);

        h_HTgenVsHTreb3 = new TH2F("h_HTgenVsHTreb3","h_HTgenVsHTreb3", 50, 0., 1000., 50, 0., 1000.);
        h_HTgenVsHTreb4 = new TH2F("h_HTgenVsHTreb4","h_HTgenVsHTreb4", 50, 0., 1000., 50, 0., 1000.);
        h_HTgenVsHTreb5 = new TH2F("h_HTgenVsHTreb5","h_HTgenVsHTreb5", 50, 0., 1000., 50, 0., 1000.);

        h_MinDphiJJgenVsMinDphiJJreb3 = new TH2F("h_MinDphiJJgenVsMinDphiJJreb3","h_MinDphiJJgenVsMinDphiJJreb3", 50, 0., TMath::Pi(), 50, 0., TMath::Pi());
        h_MinDphiJJgenVsMinDphiJJreb4 = new TH2F("h_MinDphiJJgenVsMinDphiJJreb4","h_MinDphiJJgenVsMinDphiJJreb4", 50, 0., TMath::Pi(), 50, 0., TMath::Pi());
        h_MinDphiJJgenVsMinDphiJJreb5 = new TH2F("h_MinDphiJJgenVsMinDphiJJreb5","h_MinDphiJJgenVsMinDphiJJreb5", 50, 0., TMath::Pi(), 50, 0., TMath::Pi());
    }

}

void RandS::SlaveBegin(TTree * /*tree*/)
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

}

Bool_t RandS::Process(Long64_t entry)
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

    if (isinf(*Weight)) return 0;
    //std::cout << "Weight: " << *Weight << std::endl;
    float eventWeight = *Weight;
    if (isMC_) eventWeight *= lumi_ / AvailableEvents[*DatasetID];

    std::vector<MyJet> Jets_gen;
    std::vector<MyJet> Jets_rec;
    std::vector<MyElectron> Electrons_rec;
    std::vector<MyMuon> Muons_rec;
    std::vector<MyTau> Taus_rec;
    std::vector<MyPhoton> Photons_rec;

    int NJets = JetPt->size();
    //std::cout << "Njets: " << NJets << std::endl;

    for (int i = 0; i < NJets; ++i) {
        bool OR = JetPassOR->at(i);
        float pt = JetPt->at(i);
        float eta = JetEta->at(i);
        float phi = JetPhi->at(i);
        float m = JetM->at(i);
        float jvt = JetJVT->at(i);
        if (fabs(eta) > 2.4) jvt = 1.;
        if (pt > 60.) jvt = 1.;
        if (smearCollection_ == "Gen") jvt = 1.;
        bool fjvt = JetFJVT->at(i);
        if (smearCollection_ == "Gen") fjvt = true;
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

        //// To avoid that after R+S bad jets may end up in MHT ...
        if (jet.Pt() > 20. && jet.IsBad()) {
            std::cout << "Reject event because of bad jet!" << std::endl;
            return 1;
        }

        if (jet.Pt() > 20. && jet.Pt() < 60. && fabs(jet.Eta()) < 2.4) {
            if (jet.IsBad() && jet.IsNoPU(jvtcut_)) {
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

        // Keep all good jets (important for rebalancing, and not critical since also PU events should be balanced in pT)
        if (doJVT_) {
            if (jet.IsGood()) Jets_rec.push_back(jet);
        } else { // Keep only signal jets and don't worry about jvt any further
            if (jet.IsGood() && jet.Pt() > 20. && fabs(jet.Eta()) < JetsEta_ && (jet.Pt() > 60. || jet.IsNoPU(jvtcut_) || fabs(jet.Eta()) > 2.4) ) {
                jet.SetJVT(1.);
                Jets_rec.push_back(jet);
            }
        }

    }
    GreaterByPt<MyJet> ptComparator_;
    std::sort(Jets_rec.begin(), Jets_rec.end(), ptComparator_);


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
        genjet.SetJVT(1.);
        genjet.SetJetID(true);
        Jets_gen.push_back(genjet);
    }
    std::sort(Jets_gen.begin(), Jets_gen.end(), ptComparator_);

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
        Muons_rec.push_back(muon);
    }
    GreaterByPt<MyMuon> ptComparator2_;
    std::sort(Muons_rec.begin(), Muons_rec.end(), ptComparator2_);

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
            return 1;
        }

        /*
                if (OR) {
                    //std::cout << "Reject event because of baseline electon!" << std::endl;
                    lv = true;
                    if (!cutFlowStudies) return 1;
                }
        */

        MyElectron electron(pt, eta, phi);
        electron.SetIsSignal(signal);
        electron.SetPassOR(OR);
        electron.SetCharge(q);
        Electrons_rec.push_back(electron);
    }
    GreaterByPt<MyElectron> ptComparator3_;
    std::sort(Electrons_rec.begin(), Electrons_rec.end(), ptComparator3_);

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
        Taus_rec.push_back(tau);
    }
    GreaterByPt<MyTau> ptComparator4_;
    std::sort(Taus_rec.begin(), Taus_rec.end(), ptComparator4_);

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
        Photons_rec.push_back(photon);
    }
    GreaterByPt<MyPhoton> ptComparator5_;
    std::sort(Photons_rec.begin(), Photons_rec.end(), ptComparator5_);

    //// Calculate some MET related quantities

    double recoHT = calcHT(Jets_rec, doJVT_);
    TLorentzVector recoMHT = calcMHT(Jets_rec, JetsMHTPt_, JetsMHTEta_, doJVT_);
    TLorentzVector recoMHTreb = calcMHT(Jets_rec, rebalancedJetPt_, JetsMHTEta_, false);
    TLorentzVector recoMHTall = calcMHT(Jets_rec, 0., JetsMHTEta_, false);

    TLorentzVector genMHT = calcMHT(Jets_gen, JetsMHTPt_, JetsMHTEta_, false);

    TLorentzVector MET;
    MET.SetPtEtaPhiM(*MET_pt, 0, *MET_phi, 0);

    TLorentzVector METmu;
    METmu.SetPtEtaPhiM(*METmu_pt, 0, *METmu_phi, 0);

    TLorentzVector METtrack;
    METtrack.SetPtEtaPhiM(*METtrack_pt, 0, *METtrack_phi, 0);
    METsoft_seed = METtrack.Pt();


    TLorentzVector METpho;
    METpho.SetPtEtaPhiM(*METgamma_pt, 0, *METgamma_phi, 0);

    TLorentzVector METele;
    METele.SetPtEtaPhiM(*METele_pt, 0, *METele_phi, 0);

    TLorentzVector METjet;
    METjet.SetPtEtaPhiM(*METjet_pt, 0, *METjet_phi, 0);

    //// if you want to compare to MET without muon term
    //MET = MET - METmu;

    TLorentzVector METsoft = METtrack;

    float HTSeedJVT = calcHT(Jets_rec, true);
    //cout << "HTSeedJVT: " << HTSeedJVT << endl;
    float HTSeednoJVT = calcHT(Jets_rec, false);
    //cout << "HTSeednoJVT: " << HTSeednoJVT << endl;
    float METSig = 99999;
    if (HTSeedJVT > 0) METSig = (MET.Pt()/sqrt(HTSeedJVT));
    METsig_seed = METSig;
    if (debug_ > 0) cout << "MET significance: " << METSig << endl;
    if (debug_ > 0) cout << "HTSeedJVT: " << HTSeedJVT << endl;
    float MHTSig = 99999;
    if (HTSeednoJVT > 0) MHTSig = (recoMHTall.Pt()/sqrt(HTSeednoJVT));
    MHTsig_seed = MHTSig;
    if (debug_ > 0) cout << "MHT significance: " << MHTSig << endl;
    if (debug_ > 0) cout << "HTSeednoJVT: " << HTSeednoJVT << endl;

    if (debug_ > 0) {
        cout << "---------------" << endl;
        cout << "recoMHTall (pt,phi): " << recoMHTall.Pt() << ", " << recoMHTall.Phi() << endl;
        cout << "recoMHTreb (pt,phi): " << recoMHTreb.Pt() << ", " << recoMHTreb.Phi() << endl;
        cout << "MET (pt,phi)    : " << MET.Pt() << ", " << MET.Phi() << endl;
        cout << "METjet (pt,phi)    : " << METjet.Pt() << ", " << METjet.Phi() << endl;
        cout << "METpho (pt,phi)     : " << METpho.Pt() << ", " << METpho.Phi() << endl;
        cout << "METmu (pt,phi)     : " << METmu.Pt() << ", " << METmu.Phi() << endl;
        cout << "METtrack (pt,phi)  : " << METtrack.Pt() << ", " << METtrack.Phi() << endl;
        cout << "METsoft (pt,phi)   : " << METsoft.Pt() << ", " << METsoft.Phi() << endl;
    }

    calcPredictions(Jets_rec, MET, -2, eventWeight);

    TLorentzVector genMET;
    genMET.SetPtEtaPhiM(*GenMET_pt, 0, *GenMET_phi, 0);

    TLorentzVector trueMHTreb;
    trueMHTreb.SetPtEtaPhiM(*TrueMHT_pt, 0, *TrueMHT_phi, 0);

    //// Rebalance multi jet system

    bool isRebalanced = false;

    std::vector<MyJet> Jets_reb;
    Jets_reb.reserve(Jets_rec.size());

    //// Save reco event information
    //// This is the MC expectation, which is compared to the data driven prediction in a closure test
    float mjj = 0; //calcMjj(Jets_rec);
    float MHTjj =  0; //calcMHTjj(Jets_rec);
    float dPhijj = 0; //calcDPhijj(Jets_rec);
    float dEtajj = 0; //calcDEtajj(Jets_rec);
    float jet3Pt = 0; //calcJet3Pt(Jets_rec);

    calcJJ(Jets_rec, MHTjj, mjj, dPhijj, dEtajj, jet3Pt);

    int Nreco = Njets_pred;

    if ( HT_pred > HTSave_ && ( MET_pred > METSave_ || MHT_pred > MHTSave_ || MHTjj > MHTjjSave_ ) && Njets_pred >= NJetsSave_ && mjj > MjjSave_ && dPhijj < dPhiSave_ && dEtajj > dEtaSave_ && jet3Pt < jet3PtSave_ ) {

        //std::cout << "HT_pred: " << HT_pred<< std::endl;
        //std::cout << "MET_pred: " << MET_pred<< std::endl;
        //std::cout << "MHT_pred: " << MHT_pred<< std::endl;
        //std::cout << "MHTjj: " << MHTjj<< std::endl;
        //std::cout << "Njets_pred: " << Njets_pred<< std::endl;
        //std::cout << "mjj: " << mjj<< std::endl;
        //std::cout << "dPhijj: " << dPhijj<< std::endl;
        //std::cout << "dEtajj: " << dEtajj<< std::endl;
        //std::cout << "jet3Pt: " << jet3Pt<< std::endl;

        if ( *xe90triggered || *xe110triggered) {
            triggerWeight = 1.;
        } else {
            triggerWeight = 0.;
        }
        /*
        if (useTriggerTurnOn_) {
            Int_t binx = h_MHTvsHT_eff->GetXaxis()->FindFixBin(MET_pred);
            Int_t biny = h_MHTvsHT_eff->GetYaxis()->FindFixBin(HT_pred);
            triggerWeight = h_MHTvsHT_eff->GetBinContent(binx, biny);
        } else {
            triggerWeight = 1.;
        }
        */
        PredictionTree->Fill();
    }
    //// clean variables in tree
    weight = 0.;
    Ntries_pred = 0.;
    Njets_pred = 0;
    BTags_pred = 0;
    HT_pred = 0.;
    MHT_pred = 0.;
    MHTphi_pred = 0.;
    MET_pred = 0.;
    METphi_pred = 0.;
    JetPt_pred->clear();
    JetEta_pred->clear();
    JetPhi_pred->clear();
    JetM_pred->clear();
    DeltaPhi_pred->clear();

    float deta_min = 2.4;
    float dphi_max = 2.8;
    float mjj_min = 0.;
    bool MjjSeed = calcMjjSeed(Jets_rec, deta_min, dphi_max, mjj_min);
    //cout << "MjjSeed: " << MjjSeed << endl;

    int NJetSeed = 0;
    if (fixJVTjets_) {
        NJetSeed = calcNJets(Jets_rec, true);
    } else {
        NJetSeed = calcNJets(Jets_rec, false);
    }
    //cout << "NJetSeed: " << NJetSeed << endl;

    bool IsSeed = true;
    //if ((METjet+METpho-recoMHTreb).Pt() > 10.) IsSeed = false;
    //if (METpho.Pt() > 0.) IsSeed = false;
    //if (NJetSeed > 1) IsSeed = RebPossible(Jets_rec, recoMHTall);
    //cout << "IsSeed: " << IsSeed << endl;

    //// Seed event selection (CAREFUL!!!)

    if (IsSeed && HTSeednoJVT > HTSeedMin_ && HTSeednoJVT < HTSeedMax_ && NJetSeed >= NJetsSeedMin_ && NJetSeed <= NJetsSeedMax_&& MjjSeed ) {

        //cout << "Will be rebalanced" << endl;

        if (debug_ > 99) {

            int ii = 0;
            TLorentzVector vhigh(0,0,0,0);
            TLorentzVector vlow(0,0,0,0);
            for (vector<MyJet>::iterator it = Jets_gen.begin(); it != Jets_gen.end(); ++it) {
                if (it->Pt() < 25.) vlow += (*it);
                cout << ii << "th gen jet (pt, eta, phi, flav, good, jvt): " <<
                     it->Pt() << ", " <<
                     it->Eta() << ", " <<
                     it->Phi() << ", " <<
                     it->IsB() << ", " <<
                     it->IsPU(jvtcut_) <<
                     std::endl;
                if (it->Pt() > 25.) vhigh += (*it);
                ++ii;
            }
            std::cout << "gen JetVecHigh (pt, phi): " << vhigh.Pt() << ", " << vhigh.Phi() << std::endl;
            std::cout << "gen JetVecLow  (pt, phi): " << vlow.Pt() << ", " << vlow.Phi() << std::endl;

            ii = 0;
            TLorentzVector vhigh1(0,0,0,0);
            TLorentzVector vlow1(0,0,0,0);
            for (vector<MyJet>::iterator it = Jets_rec.begin(); it != Jets_rec.end(); ++it) {
                if (it->Pt() < 25.) vlow1 += (*it);
                cout << ii << "th reco jet (pt, eta, phi, flav, good, jvt): " <<
                     it->Pt() << ", " <<
                     it->Eta() << ", " <<
                     it->Phi() << ", " <<
                     it->IsB() << ", " <<
                     it->IsPU(jvtcut_) <<
                     std::endl;
                if (it->Pt() > 25.) vhigh1 += (*it);
                ++ii;
            }

            std::cout << "reco JetVecHigh (pt, phi): " << vhigh1.Pt() << ", " << vhigh1.Phi() << std::endl;
            std::cout << "reco JetVecLow  (pt, phi): " << vlow1.Pt() << ", " << vlow1.Phi() << std::endl;

        } // end of debug

        if (smearCollection_ == "Reco") {

            double METsoftSmeared_Pt = METsoft.Pt();
            double METsoftSmeared_Phi = METsoft.Phi();
            if (useMETsoftResolution_) {
                //// Smear soft MET
                //int yBin = h_METsoft_Pt->GetYaxis()->FindBin(recoMHTreb.Pt());
                //if (yBin > h_METsoft_Pt->GetYaxis()->GetNbins()) yBin = h_METsoft_Pt->GetYaxis()->GetNbins();
                //METsoftSmeared_Pt = h_METsoft_Pt_px.at(yBin-1)->GetRandom();

                int yBin = h_METsoft_Phi->GetYaxis()->FindBin(METsoftSmeared_Pt);
                if (yBin > h_METsoft_Phi->GetYaxis()->GetNbins()) yBin = h_METsoft_Phi->GetYaxis()->GetNbins();
                double dPhi = h_METsoft_Phi_px.at(yBin-1)->GetRandom();
                METsoftSmeared_Phi = MET.Phi() + dPhi;
            }

            TLorentzVector vMETsoft_smeared(0,0,0,0);
            vMETsoft_smeared.SetPtEtaPhiM(METsoftSmeared_Pt, 0., METsoftSmeared_Phi, 0.);
            if (debug_ > 0) cout << "METsoft smeared (pt,phi)   : " << vMETsoft_smeared.Pt() << ", " << vMETsoft_smeared.Phi() << endl;

            if (!useTrueMETsoftForRebalance_) {
                isRebalanced = RebalanceJets_KinFitter(Jets_rec, Jets_gen, Jets_reb, vMETsoft_smeared);
            } else {
                isRebalanced = RebalanceJets_KinFitter(Jets_rec, Jets_gen, Jets_reb, trueMHTreb);
            }
            if (!isRebalanced) {
                cout << "Bad event: Not possible to rebalance!" << endl;
                return kTRUE;
            }

            // sort rebalanced jets
            std::sort(Jets_reb.begin(), Jets_reb.end(), ptComparator_);

            if (controlPlots_) {

                int rr = 0;
                double HTreb3 = 0.;
                double HTreb4 = 0.;
                double HTreb5 = 0.;
                MyJet* Reb1 = 0;
                MyJet* Reb2 = 0;
                MyJet* Reb3 = 0;
                for (vector<MyJet>::iterator rt = Jets_reb.begin(); rt != Jets_reb.end(); ++rt) {
                    ++rr;
                    if (rr==1) Reb1 = &(*rt);
                    if (rr==2) Reb2 = &(*rt);
                    if (rr==3) Reb3 = &(*rt);
                    if (rr<3) HTreb3 += rt->Pt();
                    if (rr<4) HTreb4 += rt->Pt();
                    if (rr<5) HTreb5 += rt->Pt();
                }

                int gg = 0;
                double HTgen3 = 0.;
                double HTgen4 = 0.;
                double HTgen5 = 0.;
                MyJet* Gen1 = 0;
                MyJet* Gen2 = 0;
                MyJet* Gen3 = 0;
                for (vector<MyJet>::iterator gt = Jets_gen.begin(); gt != Jets_gen.end(); ++gt) {
                    ++gg;
                    if (gg==1) Gen1 = &(*gt);
                    if (gg==2) Gen2 = &(*gt);
                    if (gg==3) Gen3 = &(*gt);
                    if (gg<3) HTgen3 += gt->Pt();
                    if (gg<4) HTgen4 += gt->Pt();
                    if (gg<5) HTgen5 += gt->Pt();
                }
                h_HTgenVsHTreb3-> Fill(HTgen3, HTreb3);
                h_HTgenVsHTreb4-> Fill(HTgen4, HTreb4);
                h_HTgenVsHTreb5-> Fill(HTgen5, HTreb5);

                double dPhiGenMin = 999;
                double dPhiRebMin = 999;
                double dPhiGen12 = 999;
                double dPhiReb12 = 999;
                double dPhiGen13 = 999;
                double dPhiReb13 = 999;
                double dPhiGen23 = 999;
                double dPhiReb23 = 999;
                if (Reb1 && Reb2 && Gen1 && Gen2) {
                    if ((Gen1->Pt()>80. && Gen2->Pt() > 50) && (Reb1->Pt()>80. && Reb2->Pt() > 50)) {
                        dPhiReb12 = fabs(Reb1->DeltaPhi(*Reb2));
                        dPhiGen12 = fabs(Gen1->DeltaPhi(*Gen2));
                    }
                    if (Gen3 && Reb3) {
                        if ((Gen1->Pt()>80. && Gen3->Pt() > 50) && (Reb1->Pt()>80. && Reb3->Pt() > 50)) {
                            dPhiReb13 = fabs(Reb1->DeltaPhi(*Reb3));
                            dPhiGen13 = fabs(Gen1->DeltaPhi(*Gen3));
                        }
                        if ((Gen2->Pt()>80. && Gen3->Pt() > 50) && (Reb2->Pt()>80. && Reb3->Pt() > 50)) {
                            dPhiReb23 = fabs(Reb2->DeltaPhi(*Reb3));
                            dPhiGen23 = fabs(Gen2->DeltaPhi(*Gen3));
                        }
                    }
                }
                dPhiGenMin = dPhiGen12;
                if (dPhiGen13 < dPhiGenMin) dPhiGenMin = dPhiGen13;
                if (dPhiGen23 < dPhiGenMin) dPhiGenMin = dPhiGen23;

                dPhiRebMin = dPhiReb12;
                if (dPhiReb13 < dPhiRebMin) dPhiRebMin = dPhiReb13;
                if (dPhiReb23 < dPhiRebMin) dPhiRebMin = dPhiReb23;

                if (Nreco == 3) h_MinDphiJJgenVsMinDphiJJreb3->Fill(dPhiGenMin, dPhiRebMin);
                if (Nreco == 4) h_MinDphiJJgenVsMinDphiJJreb4->Fill(dPhiGenMin, dPhiRebMin);
                if (Nreco == 5) h_MinDphiJJgenVsMinDphiJJreb5->Fill(dPhiGenMin, dPhiRebMin);

                if (Gen3 && Reb3) {
                    h_RebRes_genPt_jet1-> Fill(Gen1->Pt(), Reb1->Pt()/Gen1->Pt());
                    h_RebRes_genPt_jet2-> Fill(Gen2->Pt(), Reb2->Pt()/Gen2->Pt());
                    h_RebRes_genPt_jet3-> Fill(Gen3->Pt(), Reb3->Pt()/Gen3->Pt());
                }

            } //end of control plots

        } else {

            isRebalanced = true;
            Jets_reb = Jets_gen; // for GenJet smearing no rebalancing is needed
            METsoft.SetPtEtaPhiM(0., 0., 0., 0.);

        }

        if (isRebalanced) {

            //// comment in, if you want to save truth event information
            //triggerWeight = 1.;
            //calcPredictions(Jets_gen, genMET, -1, eventWeight);
            //PredictionTree->Fill();
            //// clean variables in tree
            //weight = 0.;
            //Ntries_pred = 0.;
            //Njets_pred = 0;
            //BTags_pred = 0;
            //HT_pred = 0.;
            //MHT_pred = 0.;
            //MET_pred = 0.;
            //JetPt_pred->clear();
            //JetEta_pred->clear();
            //JetPhi_pred->clear();
            //JetM_pred->clear();
            //DeltaPhi_pred->clear();


            //// comment in, if you want to save rebalanced event information
            //TLorentzVector MET_reb(0., 0., 0., 0.);
            //calcPredictions(Jets_reb, MET_reb, 0, eventWeight);
            //PredictionTree->Fill();
            //// clean variables in tree
            //weight = 0.;
            //Ntries_pred = 0.;
            //Njets_pred = 0;
            //BTags_pred = 0;
            //HT_pred = 0.;
            //MHT_pred = 0.;
            //MET_pred = 0.;
            //JetPt_pred->clear();
            //JetEta_pred->clear();
            //JetPhi_pred->clear();
            //JetM_pred->clear();
            //DeltaPhi_pred->clear();

            //// Smear rebalanced multi jet system
            if (doSmearing_) {
                SmearingJets(Jets_reb, METsoft, eventWeight);
            }

        }

    }

    return kTRUE;
}

void RandS::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.

    if (controlPlots_) {

        h_RebRes_genPt_eta0->Write();
        h_RebRes_genPt_eta1->Write();
        h_RebRes_genPt_eta2->Write();

        h_RebRes_genPt_jet1->Write();
        h_RebRes_genPt_jet2->Write();
        h_RebRes_genPt_jet3->Write();

        h_HTgenVsHTreb3->Write();
        h_HTgenVsHTreb4->Write();
        h_HTgenVsHTreb5->Write();

        h_MinDphiJJgenVsMinDphiJJreb3->Write();
        h_MinDphiJJgenVsMinDphiJJreb4->Write();
        h_MinDphiJJgenVsMinDphiJJreb5->Write();

    }

    outputFile->Write();

    outputFile->Close();

    for (std::map<UInt_t, UInt_t>::iterator it=ProcessedEvents.begin(); it != ProcessedEvents.end(); ++it) {
        std::cout << it->first << " " << it->second << std::endl;
    }
}

void RandS::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

}

//--------------------------------------------------------------------------
int RandS::GetIndex(const double& x, const std::vector<double>* vec) {
    int index = -1;
    for (std::vector<double>::const_iterator it = vec->begin(); it != vec->end(); ++it) {
        if ((*it) > fabs(x))
            break;
        ++index;
    }
    if (index < 0)
        index = 0;
    if (index > (int) vec->size() - 2)
        index = vec->size() - 2;

    return index;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// dice if gen (or rebalanced jet) is reconstructed
bool RandS::IsReconstructed(const double& pt, const double& eta) {
    int i_Eta = GetIndex(eta, &EtaBinEdges_);
    int i_bin = smearFunc_->RecoEff_nob.at(i_Eta)->GetXaxis()->FindBin(pt);
    double eff = 1.;
    if (pt < 500.) {
        eff = smearFunc_->RecoEff_nob.at(i_Eta)->GetBinContent(i_bin);
    }
    double random = rand_->Rndm();
    //cout << pt << ", " << eff << ", " << (random < eff) << endl;
    return (random < eff);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// pt resolution for KinFitter
double RandS::JetResolution_Pt2(const double& pt, const double& eta) {
    double ptNotE = pt/cosh(eta);
    int i_eta = GetIndex(eta, &EtaBinEdges_);
    return pow(ptNotE, 2) * pow(smearFunc_->getSigmaPtForRebalancing(i_eta)->Eval(pt), 2);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// relative pt resolution for KinFitter
double RandS::JetResolution_Ptrel(const double& pt, const double& eta) {
    int i_eta = GetIndex(eta, &EtaBinEdges_);
    return smearFunc_->getSigmaPtScaledForRebalancing(i_eta)->Eval(pt);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// eta resolution for KinFitter
double RandS::JetResolution_Eta(const double& pt, const double& eta) {
    int i_eta = GetIndex(eta, &EtaBinEdges_);
    int i_Pt = GetIndex(pt, &PtBinEdges_);
    return smearFunc_->SigmaEta.at(0).at(i_eta).at(i_Pt);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// phi resolution for KinFitter
double RandS::JetResolution_Phi(const double& pt, const double& eta) {
    int i_eta = GetIndex(eta, &EtaBinEdges_);
    int i_Pt = GetIndex(pt, &PtBinEdges_);
    return smearFunc_->SigmaPhi.at(0).at(i_eta).at(i_Pt);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// pt resolution for smearing
double RandS::JetResolutionHist_Pt_Smear(const double& pt, const double& eta, const int& i_flav) {
    int i_Pt = GetIndex(pt, &PtBinEdges_);
    int i_eta = GetIndex(eta, &EtaBinEdges_);
    double res = 1.0;
    res = smearFunc_->getSmearFunc(i_flav, i_eta, i_Pt)->GetRandom();
    return res;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// randomized muonic response of a jet
double RandS::MuResponse(const double& pt, const double& res, const int& i_flav) {
    int i_Pt = GetIndex(pt, &PtBinEdges_);
    int i_res = GetIndex(res, &ResBinEdges_);
    double mures = 0.;
    if (smearFunc_->muRes.at(i_flav).at(i_Pt).at(i_res)->GetEntries() > 0) {
        mures = smearFunc_->muRes.at(i_flav).at(i_Pt).at(i_res)->GetRandom();
        if (mures > (1.-res)) mures = 1.-res;
    }
    return mures;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
int RandS::calcNJets(std::vector<MyJet>& Jets, const bool& dojvt) {
    int NJets = 0;
    for (vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        if (it->IsGood() &&  it->Pt() > JetsPt_ && fabs(it->Eta())< JetsEta_ && ( !dojvt || (it->Pt() > 60. || it->Pt() < 20. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > JVTeta_)) )  {
            ++NJets;
        }
    }
    return NJets;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
int RandS::calcNBJets(std::vector<MyJet>& Jets, const bool& dojvt) {
    int NBJets = 0;
    for (vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        if (it->IsGood() && it->Pt() > BJetsPt_ && fabs(it->Eta()) < BJetsEta_ && it->IsB() && ( !dojvt || (it->Pt() > 60. || it->Pt() < 20. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > JVTeta_)) ) {
            ++NBJets;
        }
    }
    return NBJets;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double RandS::calcHT(std::vector<MyJet>& Jets, const bool& dojvt) {
    double HT = 0;
    for (vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        if (it->IsGood() && it->Pt() > JetsHTPt_ && fabs(it->Eta()) < JetsHTEta_ && ( !dojvt || (it->Pt() > 60. || it->Pt() < 20. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > JVTeta_)) ) {
            HT += it->Pt();
        }
    }
    return HT;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double RandS::calcMjj(std::vector<MyJet>& Jets) {
    double Mjj = 0;
    bool first = false;
    bool second = false;
    MyJet* firstJet = 0;
    MyJet* secondJet = 0;
    for (vector<MyJet>::iterator it = Jets.begin(); ( it != Jets.end() && !second ); ++it) {
        if (it->IsGood() && ( !doJVT_ || (it->Pt() > 60. || it->Pt() < 20. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > JVTeta_))) {
            if (!first && it->Pt() > MjjFirstPt_) {
                firstJet = &(*it);
                first = true;
            } else if (first && !second && it->Pt() > MjjSecondPt_) {
                secondJet = &(*it);
                second = true;
                break;
            }
        }
    }
    if (second) {
        Mjj = ((*firstJet) + (*secondJet)).M();
    }
    return Mjj;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double RandS::calcMHTjj(std::vector<MyJet>& Jets) {
    double MHTjj = 0;
    bool first = false;
    bool second = false;
    MyJet* firstJet = 0;
    MyJet* secondJet = 0;
    for (vector<MyJet>::iterator it = Jets.begin(); ( it != Jets.end() && !second ); ++it) {
        if (it->IsGood() && ( !doJVT_ || (it->Pt() > 60. || it->Pt() < 20. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > JVTeta_))) {
            if (!first && it->Pt() > MjjFirstPt_) {
                firstJet = &(*it);
                first = true;
            } else if (first && !second && it->Pt() > MjjSecondPt_) {
                secondJet = &(*it);
                second = true;
                break;
            }
        }
    }
    if (second) {
        MHTjj = ((*firstJet) + (*secondJet)).Pt();
    }
    return MHTjj;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double RandS::calcDPhijj(std::vector<MyJet>& Jets) {
    double dPhi = TMath::Pi();
    bool first = false;
    bool second = false;
    MyJet* firstJet = 0;
    MyJet* secondJet = 0;
    for (vector<MyJet>::iterator it = Jets.begin(); ( it != Jets.end() && !second ); ++it) {
        if (it->IsGood() && ( !doJVT_ || (it->Pt() > 60. || it->Pt() < 20. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > JVTeta_))) {
            if (!first && it->Pt() > MjjFirstPt_) {
                firstJet = &(*it);
                first = true;
            } else if (first && !second && it->Pt() > MjjSecondPt_) {
                secondJet = &(*it);
                second = true;
                break;
            }
        }
    }
    if (second) {
        dPhi = fabs((*firstJet).DeltaPhi(*secondJet));
    }
    return dPhi;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double RandS::calcDEtajj(std::vector<MyJet>& Jets) {
    double dEta = 0.;
    bool first = false;
    bool second = false;
    MyJet* firstJet = 0;
    MyJet* secondJet = 0;
    for (vector<MyJet>::iterator it = Jets.begin(); ( it != Jets.end() && !second ); ++it) {
        if (it->IsGood() && ( !doJVT_ || (it->Pt() > 60. || it->Pt() < 20. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > JVTeta_))) {
            if (!first && it->Pt() > MjjFirstPt_) {
                firstJet = &(*it);
                first = true;
            } else if (first && !second && it->Pt() > MjjSecondPt_) {
                secondJet = &(*it);
                second = true;
                break;
            }
        }
    }
    if (second) {
        dEta = fabs((*firstJet).Eta()-(*secondJet).Eta());
    }
    return dEta;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double RandS::calcJet3Pt(std::vector<MyJet>& Jets) {
    double Pt3 = 0;
    bool first = false;
    bool second = false;
    bool third = false;
    MyJet* firstJet = 0;
    MyJet* secondJet = 0;
    MyJet* thirdJet = 0;
    for (vector<MyJet>::iterator it = Jets.begin(); ( it != Jets.end() && !third); ++it) {
        if (it->IsGood() && ( !doJVT_ || (it->Pt() > 60. || it->Pt() < 20. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > JVTeta_))) {
            if (!first) {
                firstJet = &(*it);
                first = true;
            } else if (first && !second) {
                secondJet = &(*it);
                second = true;
            } else if (first && second && !third) {
                thirdJet = &(*it);
                third = true;
                break;
            }
        }
    }
    if (third) {
        Pt3 = (*thirdJet).Pt();
    }
    return Pt3;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
bool RandS::calcJJ(std::vector<MyJet>& Jets, float& MHTjj, float& Mjj, float& dPhijj, float& dEtajj, float& pT3) {
    MHTjj = 0;
    Mjj = 0;
    dPhijj = TMath::Pi();
    dEtajj = 0;
    pT3 = 0;
    bool first = false;
    bool second = false;
    bool third = false;
    MyJet* firstJet = 0;
    MyJet* secondJet = 0;
    MyJet* thirdJet = 0;
    for (vector<MyJet>::iterator it = Jets.begin(); ( it != Jets.end() && !third); ++it) {
        if (it->IsGood() && ( !doJVT_ || (it->Pt() > 60. || it->Pt() < 20. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > JVTeta_))) {
            if (!first && it->Pt() > MjjFirstPt_) {
                firstJet = &(*it);
                first = true;
            } else if (first && !second && it->Pt() > MjjSecondPt_) {
                secondJet = &(*it);
                second = true;
            } else if (first && second && !third) {
                thirdJet = &(*it);
                third = true;
                break;
            }
        }
    }
    if (second) {
        Mjj = ((*firstJet) + (*secondJet)).M();
        MHTjj = ((*firstJet) + (*secondJet)).Pt();
        dPhijj = fabs((*firstJet).DeltaPhi(*secondJet));
        dEtajj = fabs((*firstJet).Eta()-(*secondJet).Eta());
    }
    if (third) {
        pT3 = (*thirdJet).Pt();
    }
    return true;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
bool RandS::calcMjjSeed(std::vector<MyJet>& Jets, const float& dEta_min, const float& dPhi_max, const float& Mjj_min) {

    float Mjj_max = 0;
    for (vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        //// no jvt jet removal for seed event selection
        if (it->Pt() > 30. && it->IsGood()) {
            for (vector<MyJet>::iterator jt = it; jt != Jets.end(); ++jt) {
                if ( it == jt ) continue;
                if (jt->Pt() > rebalancedJetPt_ && jt->IsGood()) {
                    if ( fabs(it->Eta()-jt->Eta()) < dEta_min ) continue;
                    if ( fabs(it->DeltaPhi(*jt)) > dPhi_max ) continue;
                    float Mjj = ((*it) + (*jt)).M();
                    if ( Mjj > Mjj_min ) return true;
                }
            }
        }
    }
    return false;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
bool RandS::calcMinDeltaPhi(std::vector<MyJet>& Jets, TLorentzVector& MHT) {
    bool result = true;
    unsigned int i = 0;
    for (vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        if (it->IsGood() && it->Pt() > JetsPt_ && fabs(it->Eta()) < JetsEta_ && ( !doJVT_ ||  (it->Pt() > 60. || it->Pt() < 20. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > JVTeta_)) ) {
            if (i < JetDeltaMin_.size()) {
                if (fabs(it->DeltaPhi(MHT)) < JetDeltaMin_.at(i))
                    result = false;
                ++i;
            } else {
                break;
            }
        }
    }
    return result;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
TLorentzVector RandS::calcMHT(std::vector<MyJet>& Jets, const double& ptmin, const double& etamax, const bool& dojvt) {
    TLorentzVector MHT(0, 0, 0, 0);
    for (vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        if (it->IsGood() && it->Pt() > ptmin && fabs(it->Eta()) < etamax) {
            if (!dojvt || (it->Pt() > 60. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > JVTeta_)) {
                MHT -= *it;
            }
        }
    }
    return MHT;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// check if rebalancing is possible
bool RandS::RebPossible(std::vector<MyJet> &Jets_rec, TLorentzVector& vMET) {
    bool result = false;

    for (vector<MyJet>::iterator it = Jets_rec.begin(); it != Jets_rec.end(); ++it) {

        if ( it->IsBad() || it->Pt() < rebalancedJetPt_ || it->Pt() < vMET.Pt()/3.) continue;

        double dPhi1 = fabs(it->DeltaPhi(vMET));

        for (vector<MyJet>::iterator jt = Jets_rec.begin(); jt != Jets_rec.end(); ++jt) {

            if ( jt->IsBad() || jt->Pt() < rebalancedJetPt_ || jt->Pt() < vMET.Pt()/3. || (*it) == (*jt) ) continue;

            double dPhi2 = fabs(jt->DeltaPhi(vMET));

            if ( (dPhi1 + dPhi2) > TMath::Pi() ) return true;

        }

    }

    return result;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// rebalance the events using a kinematic fit and transverse momentum balance
bool RandS::RebalanceJets_KinFitter(std::vector<MyJet> &Jets_rec,   std::vector<MyJet> &Jets_gen, std::vector<MyJet> &Jets_reb, TLorentzVector& vMETsoft) {

    bool result = true;

    TKinFitter* myFit = new TKinFitter();

    std::vector<TLorentzVector*> lvec_m;

    std::vector<TMatrixD*> covMat_m;

    std::vector<TFitParticleEtEtaPhi*> fitted;
    std::vector<TFitParticleEtEtaPhi*> measured;
    std::map<int, MyJet*> JetMap;
    double MHTx_low = 0;
    double MHTy_low = 0;

    //// Fill measured particles to vector
    int i = 0;
    for (vector<MyJet>::iterator it = Jets_rec.begin(); it != Jets_rec.end(); ++it) {

        if ( it->IsBad() ) continue;

        if ( it->Pt() < rebalancedJetPt_ || i > rebalancedNJet_) {

            if (rebalanceMode_ == "MHTall") {
                MHTx_low -= it->Px();
                MHTy_low -= it->Py();
                MyJet rebalancedJet = (*it);
                Jets_reb.push_back(rebalancedJet);
            }

            if (rebalanceMode_ == "METsoft" && i > rebalancedNJet_) {
                MyJet rebalancedJet = (*it);
                Jets_reb.push_back(rebalancedJet);
                MHTx_low -= it->Px();
                MHTy_low -= it->Py();
            }

        } else {

            if (fixJVTjets_ && it->IsPU(jvtcut_)) {

                MyJet rebalancedJet = (*it);
                Jets_reb.push_back(rebalancedJet);
                MHTx_low -= it->Px();
                MHTy_low -= it->Py();

            } else {

                JetMap[i] = &(*it);

                // The particles before fitting
                double tmppx, tmppy, tmppz, tmpe;
                tmppx = it->Px();
                tmppy = it->Py();
                tmppz = it->Pz();
                tmpe = it->Energy();

                TLorentzVector* lv = new TLorentzVector(tmppx, tmppy, tmppz, tmpe);
                lvec_m.push_back(lv);
                TMatrixD* cM = new TMatrixD(3, 3);
                (*cM)(0, 0) = JetResolution_Pt2(it->E(), it->Eta());
                (*cM)(1, 1) = pow(JetResolution_Eta(it->E(), it->Eta()), 2);
                (*cM)(2, 2) = pow(JetResolution_Phi(it->E(), it->Eta()), 2);
                covMat_m.push_back(cM);
                char name[10];
                sprintf(name, "jet%i", i);
                TFitParticleEtEtaPhi* jet1 = new TFitParticleEtEtaPhi(name, name, lvec_m.back(), covMat_m.back());
                measured.push_back(jet1);
                TFitParticleEtEtaPhi* jet2 = new TFitParticleEtEtaPhi(name, name, lvec_m.back(), covMat_m.back());
                fitted.push_back(jet2);
                myFit->addMeasParticle(fitted.back());
                ++i;
            }
        }
    }

    //// Add momentum constraints
    double MET_constraint_x = 0.;
    double MET_constraint_y = 0.;

    if (rebalanceMode_ == "MHTall") {
        //// rebalance MHT of all jets
        MET_constraint_x = MHTx_low;
        MET_constraint_y = MHTy_low;
    } else if (rebalanceMode_ == "MHThigh") {
        //// rebalance MHT of fitted jets
        MET_constraint_x = 0.;
        MET_constraint_y = 0.;
    } else if (rebalanceMode_ == "METsoft") {
        //// rebalance MHT of fitted jets to soft MET
        MET_constraint_x = vMETsoft.Px() + MHTx_low;
        MET_constraint_y = vMETsoft.Py() + MHTy_low;
    } else {
        //// default: rebalance MHT of fitted jets
        MET_constraint_x = 0.;
        MET_constraint_y = 0.;
    }

    TFitConstraintEp* momentumConstr1 = new TFitConstraintEp("px", "px", 0, TFitConstraintEp::pX, MET_constraint_x);
    TFitConstraintEp* momentumConstr2 = new TFitConstraintEp("py", "py", 0, TFitConstraintEp::pY, MET_constraint_y);
    for (unsigned int i = 0; i < fitted.size(); ++i) {
        momentumConstr1->addParticle(fitted.at(i));
        momentumConstr2->addParticle(fitted.at(i));
    }
    myFit->addConstraint(momentumConstr1);
    myFit->addConstraint(momentumConstr2);

    //// Set fit parameters
    myFit->setVerbosity(0);
    myFit->setMaxNbIter(100);
    myFit->setMaxF(0.01 * 2);
    myFit->setMaxDeltaS(1.e-3);
    myFit->fit();

    int status = myFit->getStatus();

    //double chi2 = 0;
    //double prob = 0;
    if (status == 0) {
        //chi2 = myFit->getS();
        //int dof = myFit->getNDF();
        //prob = TMath::Prob(chi2, dof);
    } else {
        //chi2 = 99999;
        //prob = 0;
        result = false;
    }
    //cout << "status, chi2, prob: " << status << ", " << chi2 << ", " << prob << endl;

    //// Get the output of KinFitter
    for (unsigned int i = 0; i < measured.size(); ++i) {
        // create new rebalanced Jet
        MyJet rebalancedJet = *JetMap[i];
        rebalancedJet.SetPtEtaPhiM(fitted.at(i)->getCurr4Vec()->Pt(), fitted.at(i)->getCurr4Vec()->Eta(), fitted.at(i)->getCurr4Vec()->Phi(), fitted.at(i)->getCurr4Vec()->M());
        Jets_reb.push_back(rebalancedJet);

        /*
        if (controlPlots_) {
            int ii = 0;
            for (vector<MyJet>::iterator gt = Jets_gen.begin(); gt != Jets_gen.end(); ++gt) {
                ++ii;
                if (gt->DeltaR(rebalancedJet) < 0.1) {
                    if (ii > 3) continue;
                    if (fabs(gt->Eta())<EtaBinEdges_.at(1)) {
                        h_RebRes_genPt_eta0 -> Fill(gt->Pt(), rebalancedJet.Pt() / gt->Pt());
                    } else if (fabs(gt->Eta())<EtaBinEdges_.at(2)) {
                        h_RebRes_genPt_eta1-> Fill(gt->Pt(), rebalancedJet.Pt() / gt->Pt());
                    } else {
                        h_RebRes_genPt_eta2-> Fill(gt->Pt(), rebalancedJet.Pt() / gt->Pt());
                    }
                    if (ii==1) {
                        h_RebRes_genPt_jet1 -> Fill(gt->Pt(), rebalancedJet.Pt() / gt->Pt());
                    } else if (ii==2) {
                        h_RebRes_genPt_jet2-> Fill(gt->Pt(), rebalancedJet.Pt() / gt->Pt());
                    } else {
                        h_RebRes_genPt_jet3-> Fill(gt->Pt(), rebalancedJet.Pt() / gt->Pt());
                    }
                }
            }
        } // end of control plots
        */
    }

    delete myFit;
    for (unsigned int i = 0; i < measured.size(); ++i) {
        delete lvec_m.at(i);
        delete covMat_m.at(i);
        delete measured.at(i);
        delete fitted.at(i);
    }
    delete momentumConstr1;
    delete momentumConstr2;

    return result;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void RandS::SmearingJets(std::vector<MyJet> &Jets, TLorentzVector& vMETsoft, const float& w) {

    std::vector<MyJet> Jets_smeared;
    Jets_smeared.reserve(Jets.size());

    unsigned long Ntries2 = 1;
    double eventWeight = w;
    if (cleverPrescaleTreating_ == true && eventWeight > 1.) {
        //std::cout << "eventWeight: " << eventWeight << std::endl;
        Ntries2 = (unsigned long) w;
        if (Ntries2 > maxCleverWeight_) Ntries2 = maxCleverWeight_;
        eventWeight = w / double(Ntries2);
        //std::cout << "new eventWeight, Ntries2: " << eventWeight << ", " << Ntries2 << std::endl;
    }


    for (int i = 1; i <= Ntries_; ++i) {
        for (unsigned long j = 1; j <= Ntries2; ++j) {
            Jets_smeared.clear();
            int i_jet = 0;
            TLorentzVector METmu_sim(0.,0.,0.,0.);
            for (std::vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
                int i_flav = 0;
                if (it->IsB()) {
                    i_flav = 1;
                }
                if (IsReconstructed(it->Pt(), it->Eta()) || !JetEffEmulation_) {
                    if (it->Pt() > smearedJetPt_ && fabs(it->Eta()) < smearedJetEta_ && i_jet < smearedNJet_ && !(fixJVTjets_ && it->IsPU(jvtcut_))) {
                        double scale = JetResolutionHist_Pt_Smear(it->E(), it->Eta(), i_flav);
                        if (scale <= 0) scale = 1.;
                        if (doMETmu_ && fabs(it->Eta()) < 2.4) {
                            double mu_scale = MuResponse(it->E(), scale, i_flav);
                            METmu_sim -= (*it)*mu_scale;
                        }
                        double newE = it->Energy() * scale;
                        double newMass = it->M() * scale;
                        double newEta = it->Eta();
                        //double newEta = rand_->Gaus(it->Eta(), JetResolution_Eta(it->E(), it->Eta()));
                        double newPhi = it->Phi();
                        //double newPhi = rand_->Gaus(it->Phi(), JetResolution_Phi(it->E(), it->Eta()));
                        double newPt = sqrt(newE*newE-newMass*newMass)/cosh(newEta);

                        MyJet smearedJet;
                        smearedJet.SetPtEtaPhiM(newPt, newEta, newPhi, newMass);
                        smearedJet.SetBTag(it->IsB());
                        smearedJet.SetJVT(it->GetJVT());
                        smearedJet.SetJetID(true);
                        Jets_smeared.push_back(smearedJet);
                        ++i_jet;
                    } else {
                        MyJet smearedJet = (*it);
                        Jets_smeared.push_back(smearedJet);
                    }
                }
            }
            GreaterByPt<MyJet> ptComparator_;
            std::sort(Jets_smeared.begin(), Jets_smeared.end(), ptComparator_);

            TLorentzVector vMETpred = calcMHT(Jets_smeared, 0., JetsMHTEta_, doJVT_ );
            TLorentzVector vMHTnoJVTpred = calcMHT(Jets_smeared, 0., JetsMHTEta_, false );
            float HTnoJVTpred = calcHT(Jets_smeared, false );

            if (rebalanceMode_ == "METsoft") vMETpred += vMETsoft;
            vMETpred += METmu_sim;
            /*
            if (vMETpred.Pt() > 100. && (vMETpred-METmu_sim).Pt() < 100.) {
            	if ( vMETpred.Pt() > (vMETpred-METmu_sim).Pt() ) std::cout << "larger:  MET, METsoft, METmu, METnomu: " << vMETpred.Pt() << ", " <<  vMETsoft.Pt() << ", " << METmu_sim.Pt() << ", " << (vMETpred-METmu_sim).Pt() << std::endl;
            	if ( vMETpred.Pt() < (vMETpred-METmu_sim).Pt() ) std::cout << "smaller: MET, METsoft, METmu, METnomu: " << vMETpred.Pt() << ", " <<  vMETsoft.Pt() << ", " << METmu_sim.Pt() << ", " << (vMETpred-METmu_sim).Pt() << std::endl;
            }
            if (vMETpred.Pt() < 100. && (vMETpred-METmu_sim).Pt() > 100.) {
            	if ( vMETpred.Pt() > (vMETpred-METmu_sim).Pt() ) std::cout << "larger:  MET, METsoft, METmu, METnomu: " << vMETpred.Pt() << ", " <<  vMETsoft.Pt() << ", " << METmu_sim.Pt() << ", " << (vMETpred-METmu_sim).Pt() << std::endl;
            	if ( vMETpred.Pt() < (vMETpred-METmu_sim).Pt() ) std::cout << "smaller: MET, METsoft, METmu, METnomu: " << vMETpred.Pt() << ", " <<  vMETsoft.Pt() << ", " << METmu_sim.Pt() << ", " << (vMETpred-METmu_sim).Pt() << std::endl;
            }
            */
            calcPredictions(Jets_smeared, vMETpred, i, eventWeight);

            float mjj = 0; //calcMjj(Jets_smeared);
            float MHTjj = 0; // calcMHTjj(Jets_smeared);
            float dPhijj = 0; // calcDPhijj(Jets_smeared);
            float dEtajj = 0; // calcDEtajj(Jets_smeared);
            float jet3Pt = 0; // calcJet3Pt(Jets_smeared);
            calcJJ(Jets_smeared, MHTjj, mjj, dPhijj, dEtajj, jet3Pt);

            if( HT_pred > HTSave_ && ( MET_pred > METSave_ || MHT_pred > MHTSave_ || MHTjj > MHTjjSave_ ) && Njets_pred >= NJetsSave_ && mjj > MjjSave_ && dPhijj < dPhiSave_ && dEtajj > dEtaSave_ && jet3Pt < jet3PtSave_) {

                if (debug_ && MET_pred > METSave_ && mjj > MjjSave_ && dPhijj < 1.8 && jet3Pt < 25. ) {
                    cout << "ALARM !!!! " << endl;
                    int ii = 0;
                    TLorentzVector vhigh(0,0,0,0);
                    TLorentzVector vlow(0,0,0,0);
                    for (vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
                        if (it->Pt() < 25.) vlow += (*it);
                        cout << ii << "th rebalanced jet (pt, eta, phi, flav, jvt): " <<
                             it->Pt() << ", " <<
                             it->Eta() << ", " <<
                             it->Phi() << ", " <<
                             it->IsB() << ", " <<
                             it->IsPU(jvtcut_) <<
                             std::endl;
                        if (it->Pt() > 25.) vhigh += (*it);
                        ++ii;
                    }
                    std::cout << "rebalanced JetVecHigh (pt, phi): " << vhigh.Pt() << ", " << vhigh.Phi() << std::endl;
                    std::cout << "rebalanced JetVecLow  (pt, phi): " << vlow.Pt() << ", " << vlow.Phi() << std::endl;

                    ii = 0;
                    TLorentzVector vhigh1(0,0,0,0);
                    TLorentzVector vlow1(0,0,0,0);
                    for (vector<MyJet>::iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
                        if (it->Pt() < 25.) vlow1 += (*it);
                        cout << ii << "th smeared jet (pt, eta, phi, flav, jvt): " <<
                             it->Pt() << ", " <<
                             it->Eta() << ", " <<
                             it->Phi() << ", " <<
                             it->IsB() << ", " <<
                             it->IsPU(jvtcut_) <<
                             std::endl;
                        if (it->Pt() > 25.) vhigh1 += (*it);
                        ++ii;
                    }
                    std::cout << "smeared JetVecHigh (pt, phi): " << vhigh1.Pt() << ", " << vhigh1.Phi() << std::endl;
                    std::cout << "smeared JetVecLow  (pt, phi): " << vlow1.Pt() << ", " << vlow1.Phi() << std::endl;
                    cout << "MET_pred, MHTjj: " << MET_pred << ", " << MHTjj << endl;
                    std::cout << std::endl;
                } //end of debug

                if (useTriggerTurnOn_) {
                    Int_t binx = h_MHTvsHT_eff->GetXaxis()->FindFixBin(vMHTnoJVTpred.Pt());
                    Int_t biny = h_MHTvsHT_eff->GetYaxis()->FindFixBin(HTnoJVTpred);
                    triggerWeight = h_MHTvsHT_eff->GetBinContent(binx, biny);
                    //if (MET_pred > 200) triggerWeight = 1.;
                } else {
                    triggerWeight = 1.;
                }

                PredictionTree->Fill();
            }

            // clean variables in tree
            weight = 0.;
            Ntries_pred = 0.;
            Njets_pred = 0;
            BTags_pred = 0;
            HT_pred = 0.;
            MHT_pred = 0.;
            MET_pred = 0.;
            JetPt_pred->clear();
            JetEta_pred->clear();
            JetPhi_pred->clear();
            JetM_pred->clear();
            DeltaPhi_pred->clear();
        }
    }

    return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void RandS::calcPredictions(std::vector<MyJet>& Jets, TLorentzVector& vMET, const int& i, const float& w) {

    int NJets = calcNJets(Jets, doJVT_ );
    int NBJets = calcNBJets(Jets, doJVT_ );
    double HT = calcHT(Jets, doJVT_ );
    TLorentzVector vMHT = calcMHT(Jets, JetsMHTPt_, JetsMHTEta_, doJVT_ );
    double MHT = vMHT.Pt();
    double MET = vMET.Pt();
    double MHTphi = vMHT.Phi();
    double METphi = vMET.Phi();

    weight = w;
    Ntries_pred = i;
    Njets_pred = NJets;
    BTags_pred = NBJets;
    HT_pred = HT;
    MHT_pred = MHT;
    MHTphi_pred = MHTphi;
    MET_pred = MET;
    METphi_pred = METphi;
    calcLeadingJetPredictions(Jets, vMET);

    return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void RandS::calcLeadingJetPredictions(std::vector<MyJet>& Jets, TLorentzVector& vMHT) {
    int NJets = 0;
    for (vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        //// Fill all leading jets to output ntuple, and not only those passing jet counting thresholds
        if (it->IsGood() && fabs(it->Eta()) < JetsEta_ && ( !doJVT_ || (it->Pt() > 60. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > JVTeta_)) ) {
            ++NJets;

            if( NJets <= NJetsStored_ ) {
                JetPt_pred->push_back(it->Pt());
                JetEta_pred->push_back(it->Eta());
                JetPhi_pred->push_back(it->Phi());
                JetM_pred->push_back(it->M());
                double dphi = fabs(it->DeltaPhi(vMHT));
                DeltaPhi_pred->push_back(dphi);
            } else {
                return;
            }
        }
    }
    while ( NJets < NJetsStored_ ) {
        ++NJets;
        JetPt_pred->push_back(0.);
        JetEta_pred->push_back(0.);
        JetPhi_pred->push_back(0.);
        JetM_pred->push_back(0.);
        DeltaPhi_pred->push_back(999.);

    }

    return;
}
//--------------------------------------------------------------------------

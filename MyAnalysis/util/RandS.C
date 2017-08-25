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

    isMC_ = true;
    jvtcut_= 0.59; //// 0.59 (medium)
    lumi_ = 32900.;
    smearingfile_ = "/afs/desy.de/user/c/csander/xxl-af-cms/testarea/2.4.8/MyAnalysis/util/resolutions_GenJetMuNu_RecoNoMu_E_OR_v1.root";
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
    // Reminder here E is used instead pT for binning (variabel name not changed yet, maybe later)
    PtBinEdges_ = {0,10,20,30,40,50,70,100,140,190,250,320,400,490,590,700,820,950,1090,1240,1400,1570,1750,1940,2140,2350,2600,3000};
    EtaBinEdges_ = {0.0,0.7,1.3,1.8,2.5,3.2,5.0};
    rebalancedJetPt_ = 20.;
    rebalanceMode_ = "MHTall"; // "METsoft" (should be best), "MHTall", "MHThigh"
    smearCollection_ = "Reco"; // "Gen" for truth jet smearing only; "Reco" for full R+S
    smearedJetPt_ = 0.;
    doSmearing_ = true; // only "false" for test purposes
    JetEffEmulation_ = false;
    PtBinEdges_scaling_ = {0.,3000.};
    EtaBinEdges_scaling_ = {0.,5.};
    AdditionalSmearing_ = {1.0};
    LowerTailScaling_ = {1.0};
    UpperTailScaling_ = {1.0};
    AdditionalSmearing_variation_ = 1.0;
    LowerTailScaling_variation_ = 1.0;
    UpperTailScaling_variation_ = 1.0;
    absoluteTailScaling_ = false; // depends on definition of tails scail factors ("true" for M. Schroeder definition)
    A0RMS_ = 2.5;
    A1RMS_ = 10.;
    probExtreme_ = 0.; // possibilty to emulate total jet loss
    uncertaintyName_ = "";
    useMETsoftResolution_ = false; // Smear METsoft component prior to rebalancing
    useTrueMETsoftForRebalance_ = false; // only possible on simulated events, for testing purposes (best performance)
    useTriggerTurnOn_ = true; //save event weight from trigger turn on
    METsoftResolutionFile_ = "/afs/desy.de/user/c/csander/xxl-af-cms/testarea/2.4.8/MyAnalysis/util/METsoft_resolutions.root";
    triggerTurnOnFile_ = "/afs/desy.de/user/c/csander/xxl-af-cms/testarea/2.4.8/MyAnalysis/util/TriggerStudiesOutput_data.root";
    controlPlots_ = false;
    debug_ = 0;
    outputfile_ = "RandS.root";
    cleverPrescaleTreating_ = true; // "true", to get better statistical  precision for high weight seed events
    maxCleverWeight_ = 200; // the larger, the better (but also the slower), not greater than O(1000)
    HTSeedMin_ = 0.;
    MHTSeedMax_ = 99999.;
    MHTSigSeedMax_ = 9999.;
    NJetsSeedMin_ = 0;
    NJetsStored_ = 4;
    Ntries_ = 20;
    NJetsSave_ = 0;
    MjjSave_ = 600.;
    HTSave_ = 0;
    METSave_ = 100.;
    MHTSave_ = 9999.;
    MHTjjSave_ = 9999.;
    BJetsPt_ = 40.;
    BJetsEta_ = 2.4;
    JetsPt_ = 40.;
    JetsEta_ = 2.4;
    JetsHTPt_ = 25.;
    JetsHTEta_ = 4.5;
    JetsMHTPt_ = 25.;
    JetsMHTEta_ = 4.5;
    JetDeltaMin_ = {0.5,0.5,0.3};
    MjjFirstPt_ = 80.;
    MjjSecondPt_ = 50.;
    dPhiSave_ = 2.7;
    jet3PtSave_ = 50.;

    smearFunc_ = new SmearFunction(smearingfile_,
                                   inputhistPtHF_,inputhistEtaHF_,inputhistPhiHF_,
                                   inputhistPtLF_,inputhistEtaLF_,inputhistPhiLF_,
                                   PtBinEdges_,EtaBinEdges_,
                                   PtBinEdges_scaling_,EtaBinEdges_scaling_,
                                   AdditionalSmearing_,LowerTailScaling_,UpperTailScaling_,AdditionalSmearing_variation_,LowerTailScaling_variation_,UpperTailScaling_variation_,absoluteTailScaling_,
                                   A0RMS_,A1RMS_,probExtreme_
                                  );

    if( useMETsoftResolution_ ) {
        TFile *f_METsoft = new TFile(METsoftResolutionFile_.c_str(), "READ", "", 0);
        h_METsoft_Pt = (TH2F*) f_METsoft->FindObjectAny("h_MHTtruerebPt_vs_MHTrebMinusMET");
        h_METsoft_Pt_px.resize(h_METsoft_Pt->GetYaxis()->GetNbins());
        for (int jj = 1; jj <= h_METsoft_Pt->GetYaxis()->GetNbins(); ++jj) {
            TH1D* tmp = new TH1D(*h_METsoft_Pt->ProjectionX("px", jj, jj));
            h_METsoft_Pt_px.at(jj-1) = tmp;
        }
        h_METsoft_Phi = (TH2F*) f_METsoft->FindObjectAny("h_MHTtruerebPhiRes_vs_MHTrebMinusMET");
        h_METsoft_Phi_px.resize(h_METsoft_Phi->GetYaxis()->GetNbins());
        for (int jj = 1; jj <= h_METsoft_Phi->GetYaxis()->GetNbins(); ++jj) {
            TH1D* tmp = new TH1D(*h_METsoft_Phi->ProjectionX("px", jj, jj));
            h_METsoft_Phi_px.at(jj-1) = tmp;
        }
    }

    if( useTriggerTurnOn_ ) {
        TFile *f_triggerTurnOn = new TFile(triggerTurnOnFile_.c_str(), "READ", "", 0);
        h_MHTvsHT_all  =  (TH2F*) f_triggerTurnOn->FindObjectAny("h_MET2jetvsHT_all");
        h_MHTvsHT_triggered  =  (TH2F*) f_triggerTurnOn->FindObjectAny("h_MET2jetvsHT_triggered");
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
    PredictionTree->Branch("JetPt", "std::vector<Float_t>", &JetPt_pred);
    PredictionTree->Branch("JetEta", "std::vector<Float_t>", &JetEta_pred);
    PredictionTree->Branch("JetPhi", "std::vector<Float_t>", &JetPhi_pred);
    PredictionTree->Branch("JetM", "std::vector<Float_t>", &JetM_pred);
    PredictionTree->Branch("DeltaPhi", "std::vector<Float_t>", &DeltaPhi_pred);

    //NTotEvents = fChain->GetEntries();

    // Different seed per initialization
    gRandom->SetSeed(0);
    rand_ = new TRandom3(0);

    //// Not very elegant! TODO: Store this info in and read from file

    // [v1]
    AvailableEvents[361022] = 1993647;
    AvailableEvents[361023] = 7724495;
    AvailableEvents[361024] = 7890000;
    AvailableEvents[361025] = 7977600;
    AvailableEvents[361026] = 1833400;

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

        // Keep all jets (important for rebalancing, and not critical since also PU events should be balanced in pT)
        if (jet.IsGood()) Jets_rec.push_back(jet);

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

    double recoHT = calcHT(Jets_rec);
    TLorentzVector recoMHT = calcMHT(Jets_rec, JetsMHTPt_, JetsMHTEta_, true);
    TLorentzVector recoMHTreb = calcMHT(Jets_rec, rebalancedJetPt_, 99., false);

    TLorentzVector genMHT = calcMHT(Jets_gen, JetsMHTPt_, JetsMHTEta_, false);

    TLorentzVector MET;
    MET.SetPtEtaPhiM(*MET_pt, 0, *MET_phi, 0);

    TLorentzVector METmu;
    METmu.SetPtEtaPhiM(*METmu_pt, 0, *METmu_phi, 0);

    //// if you want to compare to MET without muon term
    //MET = MET - METmu;

    TLorentzVector METsoft = MET - recoMHTreb;

    calcPredictions(Jets_rec, MET, -2, eventWeight);

    TLorentzVector genMET;
    genMET.SetPtEtaPhiM(*GenMET_pt, 0, *GenMET_phi, 0);

    TLorentzVector trueMHTreb;
    trueMHTreb.SetPtEtaPhiM(*TrueMHT_pt, 0, *TrueMHT_phi, 0);

    //// Seed event selection (CAREFUL!!!)
    //if (recoMHT.Pt() > MHTSeedMax_ || recoMHT.Pt()/sqrt(recoHT) > MHTSigSeedMax_ ) {
    //std::cout << "Reject event because of high MHT or MHT significance!" << std::endl;
    //std::cout << recoMHT.Pt() << ", " << recoMHT.Pt()/sqrt(recoHT) << std::endl;
    //return 1;
    //}

    //// Rebalance multi jet system

    bool isRebalanced = false;

    std::vector<MyJet> Jets_reb;
    Jets_reb.reserve(Jets_rec.size());

    //// Save reco event information
    //// This is the MC expectation, which is compared to the data driven prediction in a closure test
    double mjj = calcMjj(Jets_rec);
    double MHTjj = calcMHTjj(Jets_rec);
    double dPhijj = calcDPhijj(Jets_rec);
    double jet3Pt = calcJet3Pt(Jets_rec);

    if ( HT_pred > HTSave_ && ( MET_pred > METSave_ || MHT_pred > MHTSave_ || MHTjj > MHTjjSave_ ) && Njets_pred >= NJetsSave_ && mjj > MjjSave_ && dPhijj < dPhiSave_ && jet3Pt < jet3PtSave_) {
        //std::cout << "HT_pred: " << HT_pred<< std::endl;
        //std::cout << "MET_pred: " << MET_pred<< std::endl;
        //std::cout << "MHT_pred: " << MHT_pred<< std::endl;
        //std::cout << "MHTjj: " << MHTjj<< std::endl;
        //std::cout << "Njets_pred: " << Njets_pred<< std::endl;
        //std::cout << "mjj: " << mjj<< std::endl;
        //std::cout << "dPhijj: " << dPhijj<< std::endl;
        //std::cout << "jet3Pt: " << jet3Pt<< std::endl;
        if (useTriggerTurnOn_) {
            Int_t binx = h_MHTvsHT_eff->GetXaxis()->FindFixBin(MET_pred);
            Int_t biny = h_MHTvsHT_eff->GetYaxis()->FindFixBin(HT_pred);
            triggerWeight = h_MHTvsHT_eff->GetBinContent(binx, biny);
        } else {
            triggerWeight = 1.;
        }
        PredictionTree->Fill();
    }
    //// clean variables in tree
    //weight = 0.;
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

    float deta_min = 2.9;
    float dphi_max = 2.8;
    float mjj_min = 400.;
    bool MjjSeed = calcMjjSeed(Jets_rec, deta_min, dphi_max, mjj_min);

    if (Jets_rec.size() >= 2 && HTSeed > HTSeedMin_ && NJetSeed >= NJetsSeedMin_ && MjjSeed) {

        if (smearCollection_ == "Reco") {

            if (!useTrueMETsoftForRebalance_) {
                isRebalanced = RebalanceJets_KinFitter(Jets_rec, Jets_reb, METsoft);
            } else {
                isRebalanced = RebalanceJets_KinFitter(Jets_rec, Jets_reb, trueMHTreb);
            }
            if (!isRebalanced) {
                cout << "Bad event: Not possible to rebalance!" << endl;
                return kTRUE;
            }

            // sort rebalanced jets
            std::sort(Jets_reb.begin(), Jets_reb.end(), ptComparator_);

        } else {

            isRebalanced = true;
            Jets_reb = Jets_gen; // for GenJet smearing no rebalancing is needed
            METsoft.SetPtEtaPhiM(0., 0., 0., 0.);

        }

        if (isRebalanced) {

            //// comment in, if you want to save truth event information
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
int RandS::calcNJets(std::vector<MyJet>& Jets) {
    int NJets = 0;
    for (vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        if (it->IsGood() &&  it->Pt() > JetsPt_ && std::abs(it->Eta())< JetsEta_ && (it->Pt() > 60. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > 2.4) )  {
            ++NJets;
        }
    }
    return NJets;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
int RandS::calcNBJets(std::vector<MyJet>& Jets) {
    int NBJets = 0;
    for (vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        if (it->IsGood() && it->Pt() > BJetsPt_ && std::abs(it->Eta()) < BJetsEta_ && it->IsB() && (it->Pt() > 60. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > 2.4) ) {
            ++NBJets;
        }
    }
    return NBJets;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double RandS::calcHT(std::vector<MyJet>& Jets) {
    double HT = 0;
    for (vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        if (it->IsGood() && it->Pt() > JetsHTPt_ && std::abs(it->Eta()) < JetsHTEta_ && (it->Pt() > 60. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > 2.4) ) {
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
        if (it->IsGood() && (it->Pt() > 60. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > 2.4)) {
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
    /*
    for (vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        if (it->IsPU(jvtcut_) && it->Pt() < 60. && fabs(it->Eta()) < 2.4) {
            if (first)
                if (fabs(firstJet->DeltaPhi(*it)) > 2.5 && it->Pt()>(firstJet->Pt()/3.)) return 0;
            if (second)
                if (fabs(secondJet->DeltaPhi(*it)) > 2.5 && it->Pt()>(secondJet->Pt()/3.)) return 0;
        }
    } */
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
        if (it->IsGood() && (it->Pt() > 60. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > 2.4)) {
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
        if (it->IsGood() && (it->Pt() > 60. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > 2.4)) {
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
bool RandS::calcMjjSeed(std::vector<MyJet>& Jets, const float& dEta_min, const float& dPhi_max, const float& Mjj_min) {

    float Mjj_max = 0;
    for (vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        //// no jvt jet removal for seed event selection
        //if (it->Pt() > 30. && it->IsGood() && (it->Pt() > 60. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > 2.4)) {
        if (it->Pt() > 30. && it->IsGood()) {
            for (vector<MyJet>::iterator jt = it; jt != Jets.end(); ++jt) {
                if ( it == jt ) continue;
                //if (jt->Pt() > 30. && jt->IsGood() && (jt->Pt() > 60. || jt->IsNoPU(jvtcut_) || fabs(jt->Eta()) > 2.4)) {
                if (jt->Pt() > 30. && jt->IsGood()) {
                    if ( fabs(it->Eta()-jt->Eta()) < dEta_min ) continue;
                    if ( fabs(it->DeltaPhi(*jt)) > dPhi_max ) continue;
                    float Mjj = ((*it) + (*jt)).M();
                    //if ( Mjj > Mjj_max ) Mjj_max = Mjj;
                    if ( Mjj > Mjj_min ) return true;
                }
            }
        }
    }
    //return (Mjj_max > Mjj_min);
    return false;
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
        if (it->IsGood() && (it->Pt() > 60. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > 2.4)) {
            if (!first) {
                //cout << "1st: " << it->Pt() << endl;
                firstJet = &(*it);
                first = true;
            } else if (first && !second) {
                //cout << "2nd: " << it->Pt() << endl;
                secondJet = &(*it);
                second = true;
            } else if (first && second && !third) {
                //cout << "3rd: " << it->Pt() << endl;
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
bool RandS::calcMinDeltaPhi(std::vector<MyJet>& Jets, TLorentzVector& MHT) {
    bool result = true;
    unsigned int i = 0;
    for (vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        if (it->IsGood() && it->Pt() > JetsPt_ && std::abs(it->Eta()) < JetsEta_ && (it->Pt() > 60. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > 2.4) ) {
            if (i < JetDeltaMin_.size()) {
                if (std::abs(it->DeltaPhi(MHT)) < JetDeltaMin_.at(i))
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
        if (it->IsGood() && it->Pt() > ptmin && std::abs(it->Eta()) < etamax) {
            if (!dojvt || (it->Pt() > 60. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > 2.4)) {
                MHT -= *it;
            }
        }
    }
    return MHT;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// rebalance the events using a kinematic fit and transverse momentum balance
bool RandS::RebalanceJets_KinFitter(std::vector<MyJet> &Jets_rec, std::vector<MyJet> &Jets_reb, TLorentzVector& vMETsoft) {

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

        if ( it->Pt() < rebalancedJetPt_ ) {

            if (rebalanceMode_ == "MHTall") {
                MHTx_low -= it->Px();
                MHTy_low -= it->Py();
                MyJet rebalancedJet = (*it);
                Jets_reb.push_back(rebalancedJet);
            }

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

        double METsoft_Pt = 0;
        double METsoft_Phi = 0;
        if (!useMETsoftResolution_) {
            METsoft_Pt = vMETsoft.Pt();
            METsoft_Phi = vMETsoft.Phi();
        } else {
            //// Smear soft MET
            int yBin = h_METsoft_Pt->GetYaxis()->FindBin(vMETsoft.Pt());
            if (yBin > h_METsoft_Pt->GetYaxis()->GetNbins()) yBin = h_METsoft_Pt->GetYaxis()->GetNbins();
            METsoft_Pt = h_METsoft_Pt_px.at(yBin-1)->GetRandom();
            double dPhi = h_METsoft_Phi_px.at(yBin-1)->GetRandom();
            METsoft_Phi = vMETsoft.Phi() + dPhi;
        }

        TLorentzVector vMETsoft_smeared(0,0,0,0);
        vMETsoft_smeared.SetPtEtaPhiM(METsoft_Pt, 0., METsoft_Phi, 0.);

        MET_constraint_x = vMETsoft_smeared.Px();
        MET_constraint_y = vMETsoft_smeared.Py();
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
        Ntries2 = (unsigned long) weight;
        if (Ntries2 > maxCleverWeight_) Ntries2 = maxCleverWeight_;
        eventWeight = w / double(Ntries2);
    }


    for (int i = 1; i <= Ntries_; ++i) {
        for (unsigned long j = 1; j <= Ntries2; ++j) {
            Jets_smeared.clear();
            int i_jet = 0;
            for (std::vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
                int i_flav = 0;
                if (it->IsB()) {
                    i_flav = 1;
                }
                if (IsReconstructed(it->Pt(), it->Eta()) || !JetEffEmulation_) {
                    if (it->Pt() > smearedJetPt_) {
                        double scale = JetResolutionHist_Pt_Smear(it->E(), it->Eta(), i_flav);
                        double newE = it->Energy() * scale;
                        double newMass = it->M() * scale;
                        //double newEta = it->Eta();
                        //double newPhi = it->Phi();
                        double newEta = rand_->Gaus(it->Eta(), JetResolution_Eta(it->E(), it->Eta()));
                        double newPhi = rand_->Gaus(it->Phi(), JetResolution_Phi(it->E(), it->Eta()));
                        double newPt = sqrt(newE*newE-newMass*newMass)/cosh(newEta);
                        MyJet smearedJet;
                        smearedJet.SetPtEtaPhiM(newPt, newEta, newPhi, newMass);
                        smearedJet.SetBTag(it->IsB());
                        smearedJet.SetJVT(it->GetJVT());
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

            //TLorentzVector vMETpred = calcMHT(Jets_smeared, 0., 5., true) + vMETsoft;
            TLorentzVector vMETpred = calcMHT(Jets_smeared, 0., 5., true);
            calcPredictions(Jets_smeared, vMETpred, i, eventWeight);

            double mjj = calcMjj(Jets_smeared);
            double MHTjj = calcMHTjj(Jets_smeared);
            double dPhijj = calcDPhijj(Jets_smeared);
            double jet3Pt = calcJet3Pt(Jets_smeared);

            if( HT_pred > HTSave_ && ( MET_pred > METSave_ || MHT_pred > MHTSave_ || MHTjj > MHTjjSave_ ) && Njets_pred >= NJetsSave_ && mjj > MjjSave_ && dPhijj < dPhiSave_ && jet3Pt < jet3PtSave_) {
                if (useTriggerTurnOn_) {
                    Int_t binx = h_MHTvsHT_eff->GetXaxis()->FindFixBin(MET_pred);
                    Int_t biny = h_MHTvsHT_eff->GetYaxis()->FindFixBin(HT_pred);
                    triggerWeight = h_MHTvsHT_eff->GetBinContent(binx, biny);
                } else {
                    triggerWeight = 1.;
                }
                //cout << "i, j, MHTjj, MET, mjj, dPhijj, jet3Pt: " << i << ", " << j << ", " << MHTjj << ", " << MET_pred << mjj << ", " << dPhijj << ", " << jet3Pt << endl;
                PredictionTree->Fill();
            }

            // clean variables in tree
            //weight = 0.;
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

    int NJets = calcNJets(Jets);
    int NBJets = calcNBJets(Jets);
    double HT = calcHT(Jets);
    TLorentzVector vMHT = calcMHT(Jets, JetsMHTPt_, JetsMHTEta_, true);
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

    if (i == -2) {
        NJetSeed = Njets_pred;
        HTSeed = HT_pred;
    }

    return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void RandS::calcLeadingJetPredictions(std::vector<MyJet>& Jets, TLorentzVector& vMHT) {
    int NJets = 0;
    for (vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        //// Fill all leading jets to output ntuple, and not only those passing jet counting thresholds
        if (it->IsGood() && fabs(it->Eta()) < JetsMHTEta_ && (it->Pt() > 60. || it->IsNoPU(jvtcut_) || fabs(it->Eta()) > 2.4) ) {
            ++NJets;

            if( NJets <= NJetsStored_ ) {
                JetPt_pred->push_back(it->Pt());
                JetEta_pred->push_back(it->Eta());
                JetPhi_pred->push_back(it->Phi());
                JetM_pred->push_back(it->M());
                double dphi = std::abs(it->DeltaPhi(vMHT));
                DeltaPhi_pred->push_back(dphi);
            } else {
                break;
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

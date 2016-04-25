#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <MyAnalysis/RandS.h>

#include <KinFitter/TKinFitter.h>
#include <KinFitter/TFitParticleEtEtaPhi.h>
#include <KinFitter/TFitConstraintEp.h>


// Infrastructure include(s):
#include <AsgTools/MessageCheck.h>

#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/Message.h"

// EDM includes:
#include "xAODEventInfo/EventInfo.h"
#include "xAODJet/JetContainer.h"
#include "xAODBTagging/BTagging.h"
#include "xAODTruth/TruthEvent.h"

#include <TSystem.h>

/// Helper macro for checking xAOD::TReturnCode return values
#ifndef EL_RETURN_CHECK
#define EL_RETURN_CHECK( CONTEXT, EXP )                     \
	   do {                                                     \
	      if( ! EXP.isSuccess() ) {                             \
	         Error( CONTEXT,                                    \
	                XAOD_MESSAGE( "Failed to execute: %s" ),    \
	                #EXP );                                     \
	         return EL::StatusCode::FAILURE;                    \
	      }                                                     \
	   } while( false )
#endif

// this is needed to distribute the algorithm to the workers
ClassImp(RandS)



RandS :: RandS ()
{
    jetTag_ = "AntiKt4EMTopoJets";
    genJetTag_ ="AntiKt4TruthJets";
    btagTag_ = "MV2c20";
    btagCut_ = -0.7887;
    electronTag_ = "";
    muonTag_ = "";
    smearingfile_ = "$ROOTCOREBIN/data/MyAnalysis/resolutions_v3.root";
    //inputhistPtHF_ = "h_b_JetAll_ResPt";
    //inputhistEtaHF_ = "h_b_JetAll_ResEta";
    //inputhistPhiHF_ = "h_b_JetAll_ResPhi";
    //inputhistPtLF_ = "h_nob_JetAll_ResPt";
    //inputhistEtaLF_ = "h_nob_JetAll_ResEta";
    //inputhistPhiLF_ = "h_nob_JetAll_ResPhi";
    inputhistPtHF_ = "h_HF_JetAll_ResPt";
    inputhistEtaHF_ = "h_HF_JetAll_ResEta";
    inputhistPhiHF_ = "h_HF_JetAll_ResPhi";
    inputhistPtLF_ = "h_LF_JetAll_ResPt";
    inputhistEtaLF_ = "h_LF_JetAll_ResEta";
    inputhistPhiLF_ = "h_LF_JetAll_ResPhi";
    PtBinEdges_ = {0,20,30,50,80,120,170,230,300,380,470,570,680,800,1000,1300,1700,2200,2800,3500,4300,5200,6500};
    EtaBinEdges_ = {0.0,0.7,1.3,1.8,3.2,5.0};
    rebalancedJetPt_ = 15000.;
    rebalanceMode_ = "MHTall";
    smearCollection_ = "Reco";
    smearedJetPt_ = 0.;
    PtBinEdges_scaling_ = {0.,7000.};
    EtaBinEdges_scaling_ = {0.,5.};
    AdditionalSmearing_ = {1.0};
    LowerTailScaling_ = {1.0};
    UpperTailScaling_ = {1.0};
    AdditionalSmearing_variation_ = 1.0;
    LowerTailScaling_variation_ = 1.0;
    UpperTailScaling_variation_ = 1.0;
    absoluteTailScaling_ = true;
    A0RMS_ = 2.5;
    A1RMS_ = 10.;
    probExtreme_ = 0.;
    uncertaintyName_ = "";
    useRebalanceCorrectionFactors_ = false;
    useCleverRebalanceCorrectionFactors_ = false;
    RebalanceCorrectionFile_ = "RebalanceCorrection.root";
    controlPlots_ = false;
    outputfile_ = "RandS.root";
    cleverPrescaleTreating_ = true;
    HTSeedMin_ = 350.;
    NJetsSeedMin_ = 2;
    NJetsStored_ = 3;
    Ntries_ = 10;
    NJetsSave_ = 2;
    HTSave_ = 500;
    MHTSave_ = 0;
    BJetsPt_ = 40000.;
    BJetsEta_ = 2.5;
    JetsPt_ = 40000.;
    JetsEta_ = 2.5;
    JetsHTPt_ = 40000.;
    JetsHTEta_ = 2.5;
    JetsMHTPt_ = 30000.;
    JetsMHTEta_ = 5.0;
    JetDeltaMin_ = {0.5,0.5,0.3};
    // Here you put any code for the base initialization of variables,
    // e.g. initialize all pointers to 0.  Note that you should only put
    // the most basic initialization here, since this method will be
    // called on both the submission and the worker node.  Most of your
    // initialization code will go into histInitialize() and
    // initialize().

    unsigned int needed_dim = (PtBinEdges_scaling_.size() - 1) * (EtaBinEdges_scaling_.size() - 1);
    if (AdditionalSmearing_.size() != needed_dim) {
        cout << "AdditionalSmearing has not correct dimension" << endl;
    }
    if (LowerTailScaling_.size() != needed_dim) {
        cout << "LowerTailScaling has not correct dimension" << endl;
    }
    if (UpperTailScaling_.size() != needed_dim) {
        cout << "UpperTailScaling has not correct dimension" << endl;
    }
}



EL::StatusCode RandS :: setupJob (EL::Job& job)
{
    // Here you put code that sets up the job on the submission object
    // so that it is ready to work with your algorithm, e.g. you can
    // request the D3PDReader service or add output files.  Any code you
    // put here could instead also go into the submission script.  The
    // sole advantage of putting it here is that it gets automatically
    // activated/deactivated when you add/remove the algorithm from your
    // job, which may or may not be of value to you.

    // let's initialize the algorithm to use the xAODRootAccess package
    job.useXAOD ();
    EL_RETURN_CHECK( "setupJob()", xAOD::Init() ); // call before opening first file

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode RandS :: histInitialize ()
{
    // Here you do everything that needs to be done at the very
    // beginning on each worker node, e.g. create histograms and output
    // trees.  This method gets called before any input files are
    // connected.

    // get object of class SmearFunction
    // I assume there is a better way of passing config parameters, but for now it will be like this

    cout << "Before initialization of SmearFunc" << endl;

    smearFunc_ = new SmearFunction(smearingfile_,
                                   inputhistPtHF_,inputhistEtaHF_,inputhistPhiHF_,
                                   inputhistPtLF_,inputhistEtaLF_,inputhistPhiLF_,
                                   PtBinEdges_,EtaBinEdges_,
                                   PtBinEdges_scaling_,EtaBinEdges_scaling_,
                                   AdditionalSmearing_,LowerTailScaling_,UpperTailScaling_,AdditionalSmearing_variation_,LowerTailScaling_variation_,UpperTailScaling_variation_,absoluteTailScaling_,
                                   A0RMS_,A1RMS_,probExtreme_
                                  );

    cout << "After initialization of SmearFunc" << endl;

    //// get rebalance correction histo
    if( useRebalanceCorrectionFactors_ ) {
        TFile *f_rebCorr = new TFile(RebalanceCorrectionFile_.c_str(), "READ", "", 0);
        h_RebCorrectionFactor = (TH1F*) f_rebCorr->FindObjectAny("RebCorrection_vsReco_px");
        h_RebCorrectionFactor_b = (TH1F*) f_rebCorr->FindObjectAny("RebCorrection_vsReco_b_px");
        h_2DRebCorrectionFactor = (TH2F*) f_rebCorr->FindObjectAny("RebCorrection_vsReco");
        h_2DRebCorrectionFactor_b = (TH2F*) f_rebCorr->FindObjectAny("RebCorrection_vsReco_b");
        //// Do projections for each x-bin
        for (int ii = 1; ii <= h_2DRebCorrectionFactor_b->GetXaxis()->GetNbins(); ++ii) {
            TH1D* tmp_py = new TH1D(*h_2DRebCorrectionFactor_b->ProjectionY("py", ii, ii));
            h_2DRebCorrectionFactor_b_py.push_back(tmp_py);
        }
        for (int ii = 1; ii <= h_2DRebCorrectionFactor->GetXaxis()->GetNbins(); ++ii) {
            TH1D* tmp_py = new TH1D(*h_2DRebCorrectionFactor->ProjectionY("py", ii, ii));
            h_2DRebCorrectionFactor_py.push_back(tmp_py);
        }
    }

    // define output tree
    cout << "outputfile_: " << outputfile_ << endl;
    TFile *outputFile = wk()->getOutputFile(outputfile_);
    PredictionTree = new TTree("PredictionTree", "PredictionTree");
    PredictionTree->SetDirectory(outputFile);

    cout << PredictionTree << endl;
    PredictionTree->SetAutoSave(10000000000);
    PredictionTree->SetAutoFlush(1000000);

    // set branches for output tree
    //PredictionTree->Branch("NVtx", &vtxN);
    PredictionTree->Branch("Ntries",&Ntries_pred);
    PredictionTree->Branch("NJets",&Njets_pred);
    PredictionTree->Branch("BTags",&BTags_pred);
    PredictionTree->Branch("Weight",&weight);
    PredictionTree->Branch("HT", &HT_pred);
    PredictionTree->Branch("MHT", &MHT_pred);
    PredictionTree->Branch("JetPt", "std::vector<Float_t>", &JetPt_pred);
    PredictionTree->Branch("JetEta", "std::vector<Float_t>", &JetEta_pred);
    PredictionTree->Branch("DeltaPhi", "std::vector<Float_t>", &DeltaPhi_pred);

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode RandS :: fileExecute ()
{
    // Here you do everything that needs to be done exactly once for every
    // single file, e.g. collect a list of all lumi-blocks processed
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode RandS :: changeInput (bool firstFile)
{
    // Here you do everything you need to do when we change input files,
    // e.g. resetting branch addresses on trees.  If you are using
    // D3PDReader or a similar service this method is not needed.
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode RandS :: initialize ()
{
    // Here you do everything that you need to do after the first input
    // file has been connected and before the first event is processed,
    // e.g. create additional histograms based on which variables are
    // available in the input files.  You can also create all of your
    // histograms and trees in here, but be aware that this method
    // doesn't get called if no events are processed.  So any objects
    // you create here won't be available in the output if you have no
    // input events.

    // Different seed per initialization
    gRandom->SetSeed(0);
    rand_ = new TRandom3(0);

    // count number of events
    m_eventCounter = 0;
    m_numCleanEvents = 0;

    xAOD::TEvent* event = wk()->xaodEvent();

    // tell ANA_CHECK to return an EL::StatusCode
    ANA_CHECK_SET_TYPE (EL::StatusCode);

    // make and initialize the QuickAna tool
    std::unique_ptr<ana::QuickAna> myQuickAna (new ana::QuickAna ("quickana"));
    myQuickAna->setConfig (*this);
    quickAna = std::move (myQuickAna);
    ANA_CHECK (quickAna->initialize());

    // GRL
    m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    const char* grlFilePath = "$ROOTCOREBIN/data/MyAnalysis/data15_13TeV.periodAllYear_DetStatus-v62-pro18_DQDefects-00-01-02_PHYS_StandardGRL_All_Good.xml";
    const char* fullGRLFilePath = gSystem->ExpandPathName (grlFilePath);
    std::vector<std::string> vecStringGRL;
    vecStringGRL.push_back(fullGRLFilePath);
    ANA_CHECK(m_grl->setProperty( "GoodRunsListVec", vecStringGRL));
    EL_RETURN_CHECK("initialize()",m_grl->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
    EL_RETURN_CHECK("initialize()",m_grl->initialize());

    // as a check, let's see the number of events in our xAOD
    Info("initialize()", "Number of events = %lli", event->getEntries() ); // print long long int

    // initialize and configure the jet cleaning tool
    m_jetCleaning = new JetCleaningTool("JetCleaning");
    m_jetCleaning->msg().setLevel( MSG::DEBUG );
    EL_RETURN_CHECK("initialize()",m_jetCleaning->setProperty( "CutLevel", "LooseBad"));
    EL_RETURN_CHECK("initialize()",m_jetCleaning->setProperty("DoUgly", false));
    EL_RETURN_CHECK("initialize()",m_jetCleaning->initialize());

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode RandS :: execute ()
{

    // Some debug level for detailed root output
    //gDebug = 2;

    // Here you do everything that needs to be done on every single
    // events, e.g. read input variables, apply cuts, and fill
    // histograms and trees.  This is where most of your actual analysis
    // code will go.

    ANA_CHECK_SET_TYPE (EL::StatusCode);
    ANA_CHECK (quickAna->process ());

    // Here you do everything that needs to be done on every single
    // events, e.g. read input variables, apply cuts, and fill
    // histograms and trees.  This is where most of your actual analysis
    // code will go.

    xAOD::TEvent* event = wk()->xaodEvent();

    // print every 100 events, so we know where we are:
    if( (m_eventCounter % 1000) ==0 ) Info("execute()", "Event number = %i", m_eventCounter );
    m_eventCounter++;

    //----------------------------
    // Event information
    //---------------------------
    const xAOD::EventInfo* eventInfo = 0;
    EL_RETURN_CHECK("execute",event->retrieve( eventInfo, "EventInfo"));

    // check if the event is data or MC
    // (many tools are applied either to data or MC)
    bool isMC = false;
    // check if the event is MC
    int datasetID = 0;

    double eventWeight = 1;
    if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ) {
        isMC = true; // can do something with this later
        //  extra event-level information you might need:
        datasetID =  eventInfo->mcChannelNumber();
        //const std::vector< float > weights = eventInfo->mcEventWeights();
        //if( weights.size() > 0 ) eventWeight = weights[0];
        eventWeight = quickAna->eventWeight();
        //std::cout << "datasetID = " << datasetID << std::endl;
        double lumi = 10000;
        if (datasetID == 426133) {
            double XS = 13075;
            int events = 1968000;
            double weight = XS * lumi / events;
            eventWeight *= weight;
        }
        if (datasetID == 426134) {
            double XS = 96.076;
            int events = 1984000;
            double weight = XS * lumi / events;
            eventWeight *= weight;
        }
        if (datasetID == 426135) {
            double XS = 2.725;
            int events = 1961500;
            double weight = XS * lumi / events;
            eventWeight *= weight;
        }
        if (datasetID == 426136) {
            double XS = 0.20862;
            int events = 1967800;
            double weight = XS * lumi / events;
            eventWeight *= weight;
        }
        //std::cout << "eventWeight = " << eventWeight << std::endl;
    }

    // if data check if event passes GRL
    if(!isMC) { // it's data!
        if(!m_grl->passRunLB(*eventInfo)) {
            return EL::StatusCode::SUCCESS; // go to next event
        }
    } // end if not MC

    //------------------------------------------------------------
    // Apply event cleaning to remove events due to
    // problematic regions of the detector
    // or incomplete events.
    // Apply to data.
    //------------------------------------------------------------
    // reject event if:
    if(!isMC) {
        if(   (eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::SCT)==xAOD::EventInfo::Error ) || (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) ) )
        {
            return EL::StatusCode::SUCCESS; // go to the next event
        } // end if event flags check
    } // end if the event is data
    m_numCleanEvents++;

    // get jet container of interest
    const xAOD::JetContainer* jets = 0;
    EL_RETURN_CHECK("execute()", event->retrieve( jets, jetTag_ ));
    const xAOD::JetContainer* genjets = 0;
    EL_RETURN_CHECK("execute()", event->retrieve( genjets, genJetTag_ ));
    const xAOD::TruthParticleContainer* genparticles = 0;
    EL_RETURN_CHECK("execute()", event->retrieve( genparticles, "TruthParticles"));

    // loop over the emjets and reject events with at least one bad jet
    xAOD::JetContainer::const_iterator jet_itr = jets->begin();
    xAOD::JetContainer::const_iterator jet_end = jets->end();
    for( ; jet_itr != jet_end; ++jet_itr ) {
        //cout << (*jet_itr)->pt() << endl;
        if( !m_jetCleaning->accept( **jet_itr ) && (*jet_itr)->pt() > 20.) {
            Info("execute()", "Reject event because of a bad jet");
            return EL::StatusCode::SUCCESS;
        }
    }

    std::vector<myJet> Jets_gen;
    for ( const auto* genjet : *genjets ) {
        myJet newJet;
        TLorentzVector numuActivity(0., 0., 0., 0.);
        for ( const auto* it : *genparticles ) {
            if (abs(it->pdgId()) == 13 || abs(it->pdgId()) == 12 || abs(it->pdgId()) == 14 || abs(it->pdgId()) == 16 ) {
                double dR = genjet->p4().DeltaR(it->p4());
                if (dR < 0.4) {
                    numuActivity += it->p4();
                }
            }
        }
        newJet.momentum = genjet->p4() + numuActivity;
        newJet.btag = false;
        if (abs(genjet->getAttribute<int>("PartonTruthLabelID")) == 5) newJet.btag = true;
        Jets_gen.push_back(newJet);
    }

    std::vector<myJet> Jets_rec;
    for ( auto jet : *quickAna->jets()) {
        const xAOD::BTagging * myBTag = jet->btagging();
        double mv2val=-999.;
        myBTag->MVx_discriminant(btagTag_, mv2val);
        myJet newJet;
        newJet.momentum = jet->p4();
        newJet.btag = false;
        if (mv2val > btagCut_) newJet.btag = true;
        Jets_rec.push_back(newJet);
    }

    std::vector<myJet> Jets_reb;
    Jets_reb.reserve(Jets_rec.size());

    GreaterByPt<myJet> ptComparator_;

    //
    // Rebalance multi jet system
    //
    bool isRebalanced = false;

    FillPredictions(Jets_rec, -2, eventWeight);

    if (HTSeed > HTSeedMin_ && NJetSeed >= NJetsSeedMin_) {

        //cout << "HTSeed = " << HTSeed << endl;
        //cout << "NJetSeed = " << NJetSeed << endl;

        PredictionTree->Fill();

        if (smearCollection_ == "Reco") {
            if (Jets_rec.size() > 2) {
                isRebalanced = RebalanceJets_KinFitter(Jets_rec, Jets_reb);
            } else {
                cout << "Bad event: Not possible to rebalance becaue of too few jets!" << endl;
            }

            if (!isRebalanced && Jets_rec.size() > 2) {
                cout << "Bad event: Not possible to rebalance!" << endl;
                weight_ = 0;
            }

            // sort rebalanced jets
            std::sort(Jets_reb.begin(), Jets_reb.end(), ptComparator_);
        } else {
            isRebalanced = true;
            Jets_reb = Jets_gen; // for GenJet smearing no rebalancing is needed
        }

        if (isRebalanced) {
            FillPredictions(Jets_gen, -1, eventWeight);
            //cout << "HTgen = " << HT_pred << endl;
            //cout << "NJetgen = " << Njets_pred << endl;
            PredictionTree->Fill();

            FillPredictions(Jets_reb, 0, eventWeight);
            //cout << "HTreb = " << HT_pred << endl;
            //cout << "NJetreb = " << Njets_pred << endl;
            PredictionTree->Fill();
        }

        //
        // Smear rebalanced multi jet system
        //
        if (doSmearing_) {
            SmearingJets(Jets_reb, eventWeight);
        }

    }

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode RandS :: postExecute ()
{
    // Here you do everything that needs to be done after the main event
    // processing.  This is typically very rare, particularly in user
    // code.  It is mainly used in implementing the NTupleSvc.
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode RandS :: finalize ()
{
    // This method is the mirror image of initialize(), meaning it gets
    // called after the last event has been processed on the worker node
    // and allows you to finish up any objects you created in
    // initialize() before they are written to disk.  This is actually
    // fairly rare, since this happens separately for each worker node.
    // Most of the time you want to do your post-processing on the
    // submission node after all your histogram outputs have been
    // merged.  This is different from histFinalize() in that it only
    // gets called on worker nodes that processed input events.

    xAOD::TEvent* event = wk()->xaodEvent();

    Info("finalize()", "Number of clean events = %i", m_numCleanEvents);

    // GRL
    if (m_grl) {
        delete m_grl;
        m_grl = 0;
    }

    if( m_jetCleaning ) {
        delete m_jetCleaning;
        m_jetCleaning = 0;
    }
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode RandS :: histFinalize ()
{
    // This method is the mirror image of histInitialize(), meaning it
    // gets called after the last event has been processed on the worker
    // node and allows you to finish up any objects you created in
    // histInitialize() before they are written to disk.  This is
    // actually fairly rare, since this happens separately for each
    // worker node.  Most of the time you want to do your
    // post-processing on the submission node after all your histogram
    // outputs have been merged.  This is different from finalize() in
    // that it gets called on all worker nodes regardless of whether
    // they processed input events.
    return EL::StatusCode::SUCCESS;
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
// pt resolution for KinFitter
double RandS::JetResolution_Pt2(const double& pt, const double& eta) {
    int i_eta = GetIndex(eta, &EtaBinEdges_);
    return pow(pt, 2) * pow(smearFunc_->getSigmaPtForRebalancing(i_eta)->Eval(pt/1000.), 2);
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
double RandS::JetResolution_Eta(const double& pt, const double& eta, const int& i_flav) {
    int i_eta = GetIndex(eta, &EtaBinEdges_);
    int i_Pt = GetIndex(pt, &PtBinEdges_);
    return smearFunc_->SigmaEta.at(i_flav).at(i_eta).at(i_Pt);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// phi resolution for KinFitter
double RandS::JetResolution_Phi(const double& pt, const double& eta, const int& i_flav) {
    int i_eta = GetIndex(eta, &EtaBinEdges_);
    int i_Pt = GetIndex(pt, &PtBinEdges_);
    return smearFunc_->SigmaPhi.at(i_flav).at(i_eta).at(i_Pt);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// pt resolution for smearing
double RandS::JetResolutionHist_Pt_Smear(const double& pt, const double& eta, const int& i_flav) {
    int i_Pt = GetIndex(pt, &PtBinEdges_);
    //if (i_Pt == 0) i_Pt = 1; // first pt bin is biased by jet selection
    int i_eta = GetIndex(eta, &EtaBinEdges_);
    double res = 1.0;
    res = smearFunc_->getSmearFunc(i_flav, i_eta, i_Pt)->GetRandom();
    //cout << pt << "->" << i_Pt << endl;
    return res;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
int RandS::calcNJets(std::vector<myJet>& Jets) {
    int NJets = 0;
    for (vector<myJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        if (it->momentum.Pt() > JetsPt_ && std::abs(it->momentum.Eta()) < JetsEta_) {
            ++NJets;
        }
    }
    return NJets;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
int RandS::calcNBJets(std::vector<myJet>& Jets) {
    int NBJets = 0;
    for (vector<myJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        if (it->momentum.Pt() > BJetsPt_ && std::abs(it->momentum.Eta()) < BJetsEta_ && it->btag) {
            ++NBJets;
        }
    }
    return NBJets;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double RandS::calcHT(std::vector<myJet>& Jets) {
    double HT = 0;
    for (vector<myJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        if (it->momentum.Pt() > JetsHTPt_ && std::abs(it->momentum.Eta()) < JetsHTEta_) {
            HT += it->momentum.Pt();
        }
    }
    return HT;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
bool RandS::calcMinDeltaPhi(std::vector<myJet>& Jets, LorentzVector& MHT) {
    bool result = true;
    unsigned int i = 0;
    for (vector<myJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        if (it->momentum.Pt() > JetsPt_ && std::abs(it->momentum.Eta()) < JetsEta_) {
            if (i < JetDeltaMin_.size()) {
                if (std::abs(it->momentum.DeltaPhi(MHT)) < JetDeltaMin_.at(i))
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
TLorentzVector RandS::calcMHT(std::vector<myJet>& Jets) {
    LorentzVector MHT(0, 0, 0, 0);
    for (vector<myJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        if (it->momentum.Pt() > JetsMHTPt_ && std::abs(it->momentum.Eta()) < JetsMHTEta_) {
            //cout << "(pt, eta, phi)" << it->momentum.Pt() << "," << it->momentum.Eta() << "," << it->momentum.Phi() << endl;
            MHT -= it->momentum;
        }
    }
    return MHT;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// rebalance the events using a kinematic fit and transverse momentum balance
bool RandS::RebalanceJets_KinFitter(std::vector<myJet> &Jets_rec, std::vector<myJet> &Jets_reb) {

    bool result = true;

    TKinFitter* myFit = new TKinFitter();

    std::vector<TLorentzVector*> lvec_m;

    std::vector<TMatrixD*> covMat_m;

    std::vector<TFitParticleEtEtaPhi*> fitted;
    std::vector<TFitParticleEtEtaPhi*> measured;
    std::map<int, myJet*> JetMap;
    double MHTx_low = 0;
    double MHTy_low = 0;

    //// Fill measured particles to vector
    int i = 0;
    for (vector<myJet>::iterator it = Jets_rec.begin(); it != Jets_rec.end(); ++it) {

        if (it->momentum.Pt() < rebalancedJetPt_) {

            if (rebalanceMode_ == "MHTall") {
                MHTx_low -= it->momentum.Px();
                MHTy_low -= it->momentum.Py();
                myJet rebalancedJet = (*it);
                Jets_reb.push_back(rebalancedJet);
            }

        } else {

            JetMap[i] = &(*it);

            // The particles before fitting
            double tmppx, tmppy, tmppz, tmpe;

            if( useRebalanceCorrectionFactors_ ) {
                bool btag = (it->btag);
                tmppx = it->momentum.Px()/GetRebalanceCorrection( it->momentum.Pt()/1000., btag );
                tmppy = it->momentum.Py()/GetRebalanceCorrection( it->momentum.Pt()/1000., btag );
                tmppz = it->momentum.Pz()/GetRebalanceCorrection( it->momentum.Pt()/1000., btag );
                tmpe = it->momentum.Energy()/GetRebalanceCorrection( it->momentum.Pt()/1000., btag );
            }
            else {
                tmppx = it->momentum.Px();
                tmppy = it->momentum.Py();
                tmppz = it->momentum.Pz();
                tmpe = it->momentum.Energy();
            }

            TLorentzVector* lv = new TLorentzVector(tmppx, tmppy, tmppz, tmpe);
            lvec_m.push_back(lv);
            TMatrixD* cM = new TMatrixD(3, 3);
            (*cM)(0, 0) = JetResolution_Pt2(it->momentum.Pt(), it->momentum.Eta());
            (*cM)(1, 1) = pow(JetResolution_Eta(it->momentum.Pt()/1000, it->momentum.Eta(), 0), 2);
            (*cM)(2, 2) = pow(JetResolution_Phi(it->momentum.Pt()/1000, it->momentum.Eta(), 0), 2);
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

    double chi2 = 0;
    double prob = 0;
    if (status == 0) {
        chi2 = myFit->getS();
        int dof = myFit->getNDF();
        prob = TMath::Prob(chi2, dof);
    } else {
        chi2 = 99999;
        prob = 0;
        result = false;
    }
    //cout << "status, chi2, prob: " << status << ", " << chi2 << ", " << prob << endl;
    if (controlPlots_) {
        //h_fitProb->Fill(prob)
    };

    //// Get the output of KinFitter
    for (unsigned int i = 0; i < measured.size(); ++i) {
        // create new rebalanced Jet
        TLorentzVector newP4(fitted.at(i)->getCurr4Vec()->Px(), fitted.at(i)->getCurr4Vec()->Py(), fitted.at(i)->getCurr4Vec()->Pz(), fitted.at(i)->getCurr4Vec()->E());
        myJet rebalancedJet = *JetMap[i];
        rebalancedJet.momentum = (LorentzVector) newP4;
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
double RandS::GetRebalanceCorrection(double jet_pt, bool btag)
{
    if( jet_pt > 1000. ) jet_pt = 999.;

    if ( jet_pt > rebalancedJetPt_ ) {

        int i_bin = h_RebCorrectionFactor->FindBin(jet_pt);

        if (!useCleverRebalanceCorrectionFactors_) {
            double result = 0;
            if (btag) {
                result = h_RebCorrectionFactor_b->GetBinContent(i_bin);
            } else {
                result = h_RebCorrectionFactor->GetBinContent(i_bin);
            }
            if (result < 0.01) result = 1.;
            return result;
        } else {
            double result = 0;
            if (btag) {
                result = h_2DRebCorrectionFactor_b_py.at(i_bin-1)->GetRandom();;
            } else {
                result = h_2DRebCorrectionFactor_py.at(i_bin-1)->GetRandom();;
            }
            if (result < 0.01) result = 1.;
            return result;
        }
    }

    else return 1.;

}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void RandS::SmearingJets(std::vector<myJet> &Jets, const double& w) {

    double dPx = 0;
    double dPy = 0;

    std::vector<myJet> Jets_smeared;
    Jets_smeared.reserve(Jets.size());

    for (int i = 1; i <= Ntries_; ++i) {
        //cout << "Ntries: " << i << endl;
        int Ntries2 = 1;
        weight = w;
        if (cleverPrescaleTreating_ == true && weight > 1) {
            Ntries2 = (int) weight;
            if (Ntries2 > 100) Ntries2 = 100;
            weight = w / Ntries2;
        }
        for (int j = 1; j <= Ntries2; ++j) {
            //cout << "Ntries2: " << j << endl;
            Jets_smeared.clear();
            int i_jet = 0;
            for (std::vector<myJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
                if (it->momentum.Pt() > smearedJetPt_) {
                    int i_flav = 0;
                    if (it->btag) {
                        i_flav = 1;
                    }
                    //cout << "Original pT, eta, phi: " << it->momentum.Pt() << "," << it->momentum.Eta() << "," << it->momentum.Phi() << endl;
                    double scale = JetResolutionHist_Pt_Smear(it->momentum.Pt()/1000., it->momentum.Eta(), i_flav);
                    double newE = it->momentum.Energy() * scale;
                    double newMass = it->momentum.M() * scale;
                    //double newEta = it->momentum.Eta();
                    //double newPhi = it->momentum.Phi();
                    double newEta = rand_->Gaus(it->momentum.Eta(), JetResolution_Eta(it->momentum.Pt()/1000., it->momentum.Eta(), i_flav));
                    double newPhi = rand_->Gaus(it->momentum.Phi(), JetResolution_Phi(it->momentum.Pt()/1000., it->momentum.Eta(), i_flav));
                    double newPt = sqrt(newE*newE-newMass*newMass)/cosh(newEta);
                    //cout << "New pT, eta, phi: " << newPt << "," << newEta << "," << newPhi << endl;
                    LorentzVector newP4(0.,0.,0.,0.);
                    newP4.SetPtEtaPhiM(newPt, newEta, newPhi, it->momentum.M());
                    myJet smearedJet;
                    smearedJet.momentum = newP4;
                    smearedJet.btag = it->btag;
                    Jets_smeared.push_back(smearedJet);
                    dPx -= newP4.Px() - it->momentum.Px();
                    dPy -= newP4.Py() - it->momentum.Py();
                    ++i_jet;
                } else {
                    myJet smearedJet = (*it);
                    Jets_smeared.push_back(smearedJet);
                }
            }
            GreaterByPt<myJet> ptComparator_;
            std::sort(Jets_smeared.begin(), Jets_smeared.end(), ptComparator_);

            FillPredictions(Jets_smeared, i, weight);

            if( HT_pred > HTSave_ && MHT_pred > MHTSave_ && Njets_pred >= NJetsSave_) {
                //cout << "HTpred = " << HT_pred << endl;
                PredictionTree->Fill();
            }

            // clean variables in tree
            //weight = 0.;
            Ntries_pred = 0.;
            Njets_pred = 0;
            BTags_pred = 0;
            HT_pred = 0.;
            MHT_pred = 0.;
            JetPt_pred->clear();
            JetEta_pred->clear();
            DeltaPhi_pred->clear();
        }
    }

    return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void RandS::FillPredictions(std::vector<myJet>& Jets, const int& i, const double& w) {

    int NJets = calcNJets(Jets);
    int NBJets = calcNBJets(Jets);
    double HT = calcHT(Jets);
    LorentzVector vMHT = calcMHT(Jets);
    double MHT = vMHT.Pt();

    weight = w;
    Ntries_pred = i;
    Njets_pred = NJets;
    BTags_pred = NBJets;
    HT_pred = HT/1000.;
    MHT_pred = MHT/1000.;
    FillLeadingJetPredictions(Jets, vMHT);

    if (i == -2) {
        NJetSeed = Njets_pred;
        HTSeed = HT_pred;
    }

    //cout << "weight, Ntries_pred, Njets_pred, BTags_pred, HT_pred, MHT_pred: " << weight << "," << Ntries_pred << "," << Njets_pred << "," << BTags_pred << "," << HT_pred << "," << MHT_pred << endl;

    return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void RandS::FillLeadingJetPredictions(std::vector<myJet>& Jets, LorentzVector& vMHT) {
    int NJets = 0;
    for (vector<myJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        if (it->momentum.Pt() > JetsPt_ && std::abs(it->momentum.Eta()) < JetsEta_) {
            ++NJets;

            if( NJets <= NJetsStored_ ) {
                JetPt_pred->push_back(it->momentum.Pt()/1000.);
                JetEta_pred->push_back(it->momentum.Eta());
                double dphi = std::abs(it->momentum.DeltaPhi(vMHT));
                DeltaPhi_pred->push_back(dphi);
            }
        }
    }
    return;
}
//--------------------------------------------------------------------------

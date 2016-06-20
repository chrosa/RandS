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
#include "xAODJet/JetAuxContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/ElectronAuxContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODEgamma/PhotonAuxContainer.h"
#include "xAODMuon/Muon.h"
#include "xAODBTagging/BTagging.h"
#include "xAODTruth/TruthEvent.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"

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
    jetTag_ = "AntiKt4EMTopoJets"; // not used since moving to SUSYTOOLS
    genJetTag_ ="AntiKt4TruthJets"; // not used since moving to SUSYTOOLS
    btagTag_ = "MV2c20"; // not used since moving to SUSYTOOLS
    btagCut_ = -0.7887; // not used since moving to SUSYTOOLS
    electronTag_ = ""; // not used since moving to SUSYTOOLS
    muonTag_ = ""; // not used since moving to SUSYTOOLS
    smearingfile_ = "$ROOTCOREBIN/data/MyAnalysis/resolutions_E_withNuMu_v1.root";
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
    //// Reminder here E is used instead pT for binning (variabel name not changed yet, maybe later)
    PtBinEdges_ = {0,20,30,50,80,120,170,230,300,380,470,570,680,800,1000,1300,1700,2200,2800,3500,4300,5200,6500};
    EtaBinEdges_ = {0.0,0.7,1.3,1.8,3.2,5.0};
    rebalancedJetPt_ = 20000.;
    rebalanceMode_ = "METsoft";
    smearCollection_ = "Reco";
    smearedJetPt_ = 0.;
    JetEffEmulation_ = true;
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
    useGenMHTprob_ = false;
    genMHTprobFile_ = "$ROOTCOREBIN/data/MyAnalysis/genMHTprob_noNuMu_v1.root";
    useMETsoftResolution_ = true;
    METsoftResolutionFile_ = "$ROOTCOREBIN/data/MyAnalysis/METsoft_resolutions.root";
    controlPlots_ = false;
    outputfile_ = "RandS.root";
    outputfileMHT_ = "MHT.root";
    storeMHTtree_ = false;
    cleverPrescaleTreating_ = true;
    HTSeedMin_ = 350.;
    NJetsSeedMin_ = 2;
    NJetsStored_ = 3;
    Ntries_ = 10;
    NJetsSave_ = 3;
    HTSave_ = 500;
    METSave_ = 0;
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

    if( useGenMHTprob_ ) {
        TFile *f_GenMHTprob = new TFile(genMHTprobFile_.c_str(), "READ", "", 0);
        h_MHTtrueProb_input = (TH3F*) f_GenMHTprob->FindObjectAny("h_MHTtrueProb");
        //// Do projections for each x-bin and y-bin (doing an average of neighbouring HT bins)
        h_MHTtrueProb_pz.resize(h_MHTtrueProb_input->GetXaxis()->GetNbins());
        for (int ii = 1; ii <= h_MHTtrueProb_input->GetXaxis()->GetNbins(); ++ii) {
            h_MHTtrueProb_pz.at(ii-1).resize(h_MHTtrueProb_input->GetYaxis()->GetNbins());
            for (int jj = 1; jj <= h_MHTtrueProb_input->GetYaxis()->GetNbins(); ++jj) {
                cout << "ii, jj: " << ii << ", " << jj << endl;
                TH1D* tmp = new TH1D(*h_MHTtrueProb_input->ProjectionZ("pz", ii-1, ii+1, jj, jj));
                h_MHTtrueProb_pz.at(ii-1).at(jj-1) = tmp;
            }
        }
    }

    if( useMETsoftResolution_ ) {
        TFile *f_METsoft = new TFile(METsoftResolutionFile_.c_str(), "READ", "", 0);
        h_METsoft_resPt = (TH1F*) f_METsoft->FindObjectAny("h_METsoft_resPt");
        h_METsoft_resPhi = (TH1F*) f_METsoft->FindObjectAny("h_METsoft_resPhi");
        h_METsoft = (TH1F*) f_METsoft->FindObjectAny("h_METsoft");
	}
    
    h_MHTtrueProb = new TH3F("h_MHTtrueProb","h_MHTtrueProb", 100, 0., 10000000, 5, -0.5, 4.5, 100, -100000., 400000.);
    wk()->addOutput(h_MHTtrueProb);

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
    PredictionTree->Branch("MET", &MET_pred);
    PredictionTree->Branch("JetPt", "std::vector<Float_t>", &JetPt_pred);
    PredictionTree->Branch("JetEta", "std::vector<Float_t>", &JetEta_pred);
    PredictionTree->Branch("DeltaPhi", "std::vector<Float_t>", &DeltaPhi_pred);

    // define output tree
    cout << "outputfileMHT_: " << outputfileMHT_ << endl;
    //TFile *outputFileMHT = wk()->getOutputFile(outputfileMHT_);
    MHTTree = new TTree("MHTTree", "MHTTree");
    MHTTree->SetDirectory(outputFile);

    cout << MHTTree << endl;
    MHTTree->SetAutoSave(10000000000);
    MHTTree->SetAutoFlush(1000000);

    // set branches for output tree
    MHTTree->Branch("HTreco",&HTreco);
    MHTTree->Branch("MHTreco_pt",&MHTreco_pt);
    MHTTree->Branch("MHTreco_phi",&MHTreco_phi);
    MHTTree->Branch("MHTrecolow_pt",&MHTrecolow_pt);
    MHTTree->Branch("MHTrecolow_phi",&MHTrecolow_phi);
    MHTTree->Branch("MHTreb_pt",&MHTreb_pt);
    MHTTree->Branch("MHTreb_phi",&MHTreb_phi);
    MHTTree->Branch("MHTreblow_pt",&MHTreblow_pt);
    MHTTree->Branch("MHTreblow_phi",&MHTreblow_phi);
    MHTTree->Branch("HTgen", &HTgen);
    MHTTree->Branch("METgen_pt", &METgen_pt);
    MHTTree->Branch("METgen_phi", &METgen_phi);
    MHTTree->Branch("MHTgen_pt", &MHTgen_pt);
    MHTTree->Branch("MHTgen_phi", &MHTgen_phi);
    MHTTree->Branch("MHTgenreb_pt", &MHTgenreb_pt);
    MHTTree->Branch("MHTgenreb_phi", &MHTgenreb_phi);
    MHTTree->Branch("MHTtruereb_pt", &MHTtruereb_pt);
    MHTTree->Branch("MHTtruereb_phi", &MHTtruereb_phi);
    MHTTree->Branch("MET_pt",&MET_pt);
    MHTTree->Branch("MET_phi",&MET_phi);
    MHTTree->Branch("Weight",&weight);

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


    ST::ISUSYObjDef_xAODTool::DataSource datasource = ST::ISUSYObjDef_xAODTool::FullSim;

    //const xAOD::EventInfo* eventInfo = 0;
    //EL_RETURN_CHECK("initialize",event->retrieve( eventInfo, "EventInfo"));
    //if (eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ) {
    //    datasource = 1;
    //}

    objTool = new ST::SUSYObjDef_xAOD("SUSYObjDef_xAOD");

    prw_file_ = "DUMMY";
    std::vector<std::string> prw_conf;
    if (prw_file_ == "DUMMY") {
        prw_conf.push_back("/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/PileupReweighting/mc15ab_defaults.NotRecommended.prw.root");
        prw_conf.push_back("/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root");
    }
    else {
        prw_conf.push_back(prw_file_);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    // Configure the SUSYObjDef instance
    EL_RETURN_CHECK("initialize", objTool->setProperty("PRWConfigFiles", prw_conf) );
    EL_RETURN_CHECK("initialize", objTool->setProperty("DataSource", datasource) ) ;
    EL_RETURN_CHECK("initialize", objTool->setProperty("ConfigFile", "SUSYTools/SUSYTools_Default.conf") );

    //Manually setting additional properties will override what's in the configuration file
    EL_RETURN_CHECK("initialize", objTool->setProperty("EleId", "TightLLH") );
    EL_RETURN_CHECK("initialize", objTool->setProperty("EleBaselineId", "LooseLLH") );

    EL_RETURN_CHECK("initialize", objTool->initialize() );

    // GRL
    m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    const char* grlFilePath = "$ROOTCOREBIN/data/MyAnalysis/data15_13TeV.periodAllYear_DetStatus-v62-pro18_DQDefects-00-01-02_PHYS_StandardGRL_All_Good.xml";
    const char* fullGRLFilePath = gSystem->ExpandPathName (grlFilePath);
    std::vector<std::string> vecStringGRL;
    vecStringGRL.push_back(fullGRLFilePath);
    EL_RETURN_CHECK("initialize()",m_grl->setProperty( "GoodRunsListVec", vecStringGRL));
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

    // Here you do everything that needs to be done on every single
    // events, e.g. read input variables, apply cuts, and fill
    // histograms and trees.  This is where most of your actual analysis
    // code will go.

    // Some debug level for detailed root output
    //gDebug = 2;

    xAOD::TEvent* event = wk()->xaodEvent();

    // print every 100 events, so we know where we are:
    if( (m_eventCounter % 100) ==0 ) Info("execute()", "Event number = %i", m_eventCounter );
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
        //extra event-level information you might need:
        datasetID =  eventInfo->mcChannelNumber();
        const std::vector< float > weights = eventInfo->mcEventWeights();
        if( weights.size() > 0 ) eventWeight = weights[0];
        //std::cout << "eventWeight = " << eventWeight << std::endl;
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
        if (datasetID == 361022) {
            double XS = 2433200;
            int events = 1992912;
            double eff = 0.00033264;
            double weight = XS * lumi * eff / events;
            eventWeight *= weight;
        }
        if (datasetID == 361023) {
            double XS = 26454;
            int events = 7882487;
            double eff = 0.00031953;
            double weight = XS * lumi * eff / events;
            eventWeight *= weight;
        }
        if (datasetID == 361024) {
            double XS = 254.63;
            int events = 7979799;
            double eff = 0.00053009;
            double weight = XS * lumi * eff / events;
            eventWeight *= weight;
        }
        if (datasetID == 361025) {
            double XS = 4.5535;
            int events = 7981600;
            double eff = 0.00092325;
            double weight = XS * lumi * eff / events;
            eventWeight *= weight;
        }
        if (datasetID == 361026) {
            double XS = 0.25753;
            int events = 1893400;
            double eff = 0.00094016;
            double weight = XS * lumi * eff / events;
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

    // Electrons
    //xAOD::ElectronContainer* electrons_nominal(0);
    //xAOD::ShallowAuxContainer* electrons_nominal_aux(0);
    //EL_RETURN_CHECK("execute()", objTool->GetElectrons(electrons_nominal, electrons_nominal_aux) );

    // Photons
    //xAOD::PhotonContainer* photons_nominal(0);
    //xAOD::ShallowAuxContainer* photons_nominal_aux(0);
    //EL_RETURN_CHECK("execute()", objTool->GetPhotons(photons_nominal,photons_nominal_aux) );

    // Muons
    //xAOD::MuonContainer* muons_nominal(0);
    //xAOD::ShallowAuxContainer* muons_nominal_aux(0);
    //EL_RETURN_CHECK("execute()", objTool->GetMuons(muons_nominal, muons_nominal_aux) );

    // Jets
    xAOD::JetContainer* jets_nominal(0);
    xAOD::ShallowAuxContainer* jets_nominal_aux(0);
    EL_RETURN_CHECK("execute()", objTool->GetJets(jets_nominal, jets_nominal_aux) );

    // Taus
    //xAOD::TauJetContainer* taus_nominal(0);
    //xAOD::ShallowAuxContainer* taus_nominal_aux(0);
    //EL_RETURN_CHECK("execute()", objTool->GetTaus(taus_nominal,taus_nominal_aux) );

    // MET
    xAOD::MissingETContainer* mettst_nominal = new xAOD::MissingETContainer;
    xAOD::MissingETAuxContainer* mettst_nominal_aux = new xAOD::MissingETAuxContainer;
    mettst_nominal->setStore(mettst_nominal_aux);
    mettst_nominal->reserve(10);

    //Call Overlap Removal
    //EL_RETURN_CHECK("execute()", objTool->OverlapRemoval(electrons_nominal, muons_nominal, jets) );
    //EL_RETURN_CHECK("execute()", objTool->GetMET(*mettst_nominal,jets,electrons_nominal,muons_nominal,0,0,true) );
    EL_RETURN_CHECK("execute()", objTool->GetMET(*mettst_nominal,jets_nominal,0,0,0,0,true) );

    // get jet container of interest
    //const xAOD::JetContainer* jets = 0;
    //EL_RETURN_CHECK("execute()", event->retrieve( jets, "AntiKt4EMTopoJets" ));
    const xAOD::JetContainer* genjets = 0;
    EL_RETURN_CHECK("execute()", event->retrieve( genjets, "AntiKt4TruthJets"));
    const xAOD::TruthParticleContainer* genparticles = 0;
    EL_RETURN_CHECK("execute()", event->retrieve( genparticles, "TruthParticles"));

    // loop over the emjets and reject events with at least one bad jet
    xAOD::JetContainer::const_iterator jet_itr = jets_nominal->begin();
    xAOD::JetContainer::const_iterator jet_end = jets_nominal->end();
    for( ; jet_itr != jet_end; ++jet_itr ) {
        if( !m_jetCleaning->accept( **jet_itr ) && (*jet_itr)->pt() > 20000.) {
            Info("execute()", "Reject event because of a bad jet");
            return EL::StatusCode::SUCCESS;
        }
    }

	// Fill truth jets to simplified jet container
	// add back neutrinos and muons to truth jets
    std::vector<myJet> Jets_gen;
    for ( const auto& genjet : *genjets ) {
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

    GreaterByPt<myJet> ptComparator_;
    std::sort(Jets_gen.begin(), Jets_gen.end(), ptComparator_);

	// Fill reco jets to simplified jet container
    std::vector<myJet> Jets_rec;
    for ( const auto& jet : *jets_nominal) {
        myJet newJet;
        newJet.momentum = jet->p4();
        newJet.btag = false;
        if (objTool->IsBJet(*jet)) newJet.btag = true;
        Jets_rec.push_back(newJet);
    }
    std::sort(Jets_rec.begin(), Jets_rec.end(), ptComparator_);

    std::vector<myJet> Jets_reb;
    Jets_reb.reserve(Jets_rec.size());

    LorentzVector vrecoMHT = calcMHT(Jets_rec, JetsMHTPt_, JetsMHTEta_);
    LorentzVector vrecoMHTreb = calcMHT(Jets_rec, rebalancedJetPt_, 5.);
    LorentzVector vgenMHT = calcMHT(Jets_gen, JetsMHTPt_, JetsMHTEta_);
    
    LorentzVector MET(0,0,0,0);
    MET.SetPxPyPzE((*mettst_nominal)["Final"]->mpx(), (*mettst_nominal)["Final"]->mpx(), 0, (*mettst_nominal)["Final"]->met() );
    LorentzVector METsoft = MET - vrecoMHTreb;
    calcPredictions(Jets_rec, MET, -2, eventWeight);

    //// Fill 3D plots for true MHT hypothesis used in rebalancing
    //if (controlPlots_) {
    //    vgenMHT.SetPtEtaPhiE(vgenMHT.Pt(), 0., vgenMHT.Phi(), vgenMHT.Pt());
    //    vrecoMHT.SetPtEtaPhiE(vrecoMHT.Pt(), 0., vrecoMHT.Phi(), vrecoMHT.Pt());
    //    double genMHTp = (vgenMHT.Vect() * vrecoMHT.Vect()) / sqrt(vrecoMHT.Vect() * vrecoMHT.Vect());
    //    h_MHTtrueProb->Fill(HT_pred*1000, BTags_pred, genMHTp, eventWeight);
    //    if (vrecoMHT.Pt() > 200000.) {
    //        cout << "HT (reco): " << HT_pred << endl;
    //        cout << "BTags (reco): " << BTags_pred << endl;
    //        cout << "MHT, phi (reco): " << vrecoMHT.Pt() << ", " << vrecoMHT.Phi() << endl;
    //        cout << "MHT, phi (gen): " << vgenMHT.Pt() << ", " << vgenMHT.Phi() << endl;
    //        cout << "genMHTp, factor: " <<  genMHTp << ", " << genMHTp << endl;
    //    }
    //}

    // gen MET from all status==1 gen particles
    // true MHT from all status==1 gen particles matched (dR < 0.4) to jets above rebalancing threshold
    LorentzVector vgenMET(0, 0, 0, 0);
    LorentzVector vtrueMHTreb(0, 0, 0, 0);
    for ( const auto* it : *genparticles ) {
        if (it->pt() < 1. || it->status() != 1 || it->hasDecayVtx()) {
            //cout << "Pid, status, pt, eta: " << it->pdgId() << ", " << it->status() << ", " << it->pt() << ", " << it->eta() << endl;
            continue;
        }
        vgenMET -= it->p4();
        //cout << "MET (pt, phi): " << vgenMET.Pt() << ", " << vgenMET.Phi() << endl;
        //cout << "Pid, status, pt, eta, phi: " << it->pdgId() << ", " << it->status() << ", " << it->pt() << ", " << it->eta() << ", " << it->phi() << endl;
        bool particleAdded = false;
        for (vector<myJet>::iterator jt = Jets_rec.begin(); (jt != Jets_rec.end() && !particleAdded) ; ++jt) {
            if (jt->momentum.Pt() > rebalancedJetPt_ ) {
                double dR = jt->momentum.DeltaR(it->p4());
                if (dR < 0.4) {
                    vtrueMHTreb -= it->p4();
                    particleAdded = true;
                }
            }
        }
    }

    // reco HT
    HTreco = calcHT(Jets_rec);

    // reco MHT
    MHTreco_pt = vrecoMHT.Pt();
    MHTreco_phi = vrecoMHT.Phi();
	LorentzVector vrecoMHTall = calcMHT(Jets_rec, 0., 5.);
	LorentzVector vrecoMHTlow = vrecoMHT - vrecoMHTall;
    MHTrecolow_pt = vrecoMHTlow.Pt();
    MHTrecolow_phi = vrecoMHTlow.Phi();

    // reco MHT with rebalancing threshold
    LorentzVector vrecoMHTreblow = vrecoMHTreb - vrecoMHTall;
    MHTreb_pt = vrecoMHTreb.Pt();
    MHTreb_phi = vrecoMHTreb.Phi();
    MHTreblow_pt = vrecoMHTreblow.Pt();
    MHTreblow_phi = vrecoMHTreblow.Phi();

    // gen HT
    HTgen = calcHT(Jets_gen);

    // gen MHT
    MHTgen_pt = vgenMHT.Pt();
    MHTgen_phi = vgenMHT.Phi();

    // gen MHT with rebalancing threshold
    LorentzVector vgenMHTreb = calcMHT(Jets_gen, rebalancedJetPt_, JetsMHTEta_);
    MHTgenreb_pt = vgenMHTreb.Pt();
    MHTgenreb_phi = vgenMHTreb.Phi();

    METgen_pt = vgenMET.Pt();
    METgen_phi = vgenMET.Phi();
    MHTtruereb_pt = vtrueMHTreb.Pt();
    MHTtruereb_phi = vtrueMHTreb.Phi();

    // MET
    MET_pt = (*mettst_nominal)["Final"]->met();
    MET_phi = (*mettst_nominal)["Final"]->phi();

    if (storeMHTtree_) {
        MHTTree->Fill();
    }

    //
    // Rebalance multi jet system
    //
    bool isRebalanced = false;

    if (HTSeed > HTSeedMin_ && NJetSeed >= NJetsSeedMin_) {

        //cout << "HTSeed = " << HTSeed << endl;
        //cout << "NJetSeed = " << NJetSeed << endl;

        PredictionTree->Fill();

        if (smearCollection_ == "Reco") {
            if (Jets_rec.size() > 2) {
                isRebalanced = RebalanceJets_KinFitter(Jets_rec, Jets_reb, METsoft);
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
            calcPredictions(Jets_gen, vgenMET, -1, eventWeight);
            //cout << "HTgen = " << HT_pred << endl;
            //cout << "NJetgen = " << Njets_pred << endl;
            PredictionTree->Fill();

            LorentzVector MET_reb(0,0,0,0);
            calcPredictions(Jets_reb, MET_reb, 0, eventWeight);
            //cout << "HTreb = " << HT_pred << endl;
            //cout << "NJetreb = " << Njets_pred << endl;
            PredictionTree->Fill();
        }

        //
        // Smear rebalanced multi jet system
        //
        if (doSmearing_) {
            SmearingJets(Jets_reb, METsoft, eventWeight);
        }

    }

    // The containers created by the shallow copy are owned by you. Remember to delete them
    delete jets_nominal; // not these, we put them in the store!
    delete jets_nominal_aux;
    // delete taus_nominal;
    // delete taus_nominal_aux;
    // delete muons_nominal;
    // delete muons_nominal_aux;
    // delete electrons_nominal;
    // delete electrons_nominal_aux;
    // delete photons_nominal;
    // delete photons_nominal_aux;
    // delete metcst_nominal;
    // delete metcst_nominal_aux;
    delete mettst_nominal;
    delete mettst_nominal_aux;

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

    //xAOD::TEvent* event = wk()->xaodEvent();

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
// dice if gen (or rebalanced jet) is reconstructed
bool RandS::IsReconstructed(const double& pt, const double& eta) {
    int i_Eta = GetIndex(eta, &EtaBinEdges_);
    int i_bin = smearFunc_->RecoEff_b.at(i_Eta)->GetXaxis()->FindBin(pt);
    double eff = 1.;
    if (pt < 40.) {
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
TLorentzVector RandS::calcMHT(std::vector<myJet>& Jets, const double& ptmin, const double& etamax) {
    LorentzVector MHT(0, 0, 0, 0);
    for (vector<myJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
        if (it->momentum.Pt() > ptmin && std::abs(it->momentum.Eta()) < etamax) {
            //cout << "(pt, eta, phi)" << it->momentum.Pt() << "," << it->momentum.Eta() << "," << it->momentum.Phi() << endl;
            MHT -= it->momentum;
        }
    }
    return MHT;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// rebalance the events using a kinematic fit and transverse momentum balance
bool RandS::RebalanceJets_KinFitter(std::vector<myJet> &Jets_rec, std::vector<myJet> &Jets_reb, LorentzVector& vMETsoft) {

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

            if (rebalanceMode_ == "MHTall" && !useGenMHTprob_) {
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
            (*cM)(0, 0) = JetResolution_Pt2(it->momentum.E(), it->momentum.Eta());
            (*cM)(1, 1) = pow(JetResolution_Eta(it->momentum.E()/1000, it->momentum.Eta()), 2);
            (*cM)(2, 2) = pow(JetResolution_Phi(it->momentum.E()/1000, it->momentum.Eta()), 2);
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

    if (useGenMHTprob_) {
		
        //// Now get from simulation the expected true MHT for a given reco HT and reco MHT
        //
        //LorentzVector vrecoMHT = calcMHT(Jets_rec, JetsMHTPt_, JetsMHTEta_);
        //double recoHT = calcHT(Jets_rec);
        //int NBJets = calcNBJets(Jets_rec);
        //// Do projections for the given HT and MHT bin (including neighbouring bins)
        //int ii = h_MHTtrueProb_input->GetXaxis()->FindBin(recoHT);
        //if (ii > h_MHTtrueProb_input->GetXaxis()->GetNbins() - 1) ii = h_MHTtrueProb_input->GetXaxis()->GetNbins() - 1;
        //int jj = h_MHTtrueProb_input->GetYaxis()->FindBin((double) NBJets);
        //if (jj > h_MHTtrueProb_input->GetYaxis()->GetNbins() - 1) jj = h_MHTtrueProb_input->GetYaxis()->GetNbins() - 1;
        //double genMHTp = h_MHTtrueProb_pz.at(ii-1).at(jj-1)->GetRandom();
        //LorentzVector vgenMHT = vrecoMHT;
        //double px = vrecoMHT.Px() * genMHTp / vrecoMHT.Pt();
        //double py = vrecoMHT.Py() * genMHTp / vrecoMHT.Pt();
        //vgenMHT.SetPx(px);
        //vgenMHT.SetPy(py);
        //if (vrecoMHT.Pt() > 200000.) {
        //    cout << "vrecoMHT, phi (reco): " << vrecoMHT.Pt() << ", " << vrecoMHT.Phi() << endl;
        //    cout << "vrecoMHT, phi (gen): " << vgenMHT.Pt() << ", " << vgenMHT.Phi() << endl;
        //    cout << "genMHTp: " << genMHTp << endl;
        //}
        //MET_constraint_x = vgenMHT.Px();
        //MET_constraint_y = vgenMHT.Py();
        
    } else {
		
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

			double METsoft_Pt = vMETsoft.Pt();
            double METsoft_Phi =  vMETsoft.Phi();

			//// Smear soft MET
            //double METsoft_Pt = h_METsoft->GetRandom();
            //double METsoft_Phi = vMETsoft.Phi() + h_METsoft_resPhi->GetRandom();

            LorentzVector vMETsoft_smeared(0,0,0,0);
            vMETsoft_smeared.SetPtEtaPhiE(METsoft_Pt, 0., METsoft_Phi, METsoft_Pt);

            MET_constraint_x = vMETsoft_smeared.Px();
            MET_constraint_y = vMETsoft_smeared.Py();
        } else {
            //// default: rebalance MHT of fitted jets
            MET_constraint_x = 0.;
            MET_constraint_y = 0.;
        }
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
void RandS::SmearingJets(std::vector<myJet> &Jets, LorentzVector& vMETsoft, const double& w) {

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
                int i_flav = 0;
                if (it->btag) {
                    i_flav = 1;
                }
                if (IsReconstructed(it->momentum.Pt(), it->momentum.Eta()) || !JetEffEmulation_) {
                    if (it->momentum.Pt() > smearedJetPt_) {
                        //cout << "Original pT, eta, phi: " << it->momentum.Pt() << "," << it->momentum.Eta() << "," << it->momentum.Phi() << endl;
                        double scale = JetResolutionHist_Pt_Smear(it->momentum.E()/1000., it->momentum.Eta(), i_flav);
                        double newE = it->momentum.Energy() * scale;
                        double newMass = it->momentum.M() * scale;
                        //double newEta = it->momentum.Eta();
                        //double newPhi = it->momentum.Phi();
                        double newEta = rand_->Gaus(it->momentum.Eta(), JetResolution_Eta(it->momentum.E()/1000., it->momentum.Eta()));
                        double newPhi = rand_->Gaus(it->momentum.Phi(), JetResolution_Phi(it->momentum.E()/1000., it->momentum.Eta()));
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
            }
            GreaterByPt<myJet> ptComparator_;
            std::sort(Jets_smeared.begin(), Jets_smeared.end(), ptComparator_);

            //cout << "Reco METsoft (pt, phi): " << vMETsoft.Pt() << ", " << vMETsoft.Phi() << endl;
            double METsoft_Pt = vMETsoft.Pt();
            double METsoft_Phi =  vMETsoft.Phi();
  
            //// Smear soft MET part
            //double METsoft_Pt = h_METsoft->GetRandom();
            //double METsoft_Phi = vMETsoft.Phi() + h_METsoft_resPhi->GetRandom();
  
            LorentzVector vMETsoft_smeared(0,0,0,0);
            vMETsoft_smeared.SetPtEtaPhiE(METsoft_Pt, 0., METsoft_Phi, METsoft_Pt);
  
            //cout << "Smeared METsoft (pt, phi): " << vMETsoft_smeared.Pt() << ", " << vMETsoft_smeared.Phi() << endl;
            LorentzVector vMETpred = calcMHT(Jets_smeared, 0., 5.) + vMETsoft_smeared;
            calcPredictions(Jets_smeared, vMETpred, i, weight);

            if( HT_pred > HTSave_ && vMETpred.Pt() > METSave_ && Njets_pred >= NJetsSave_) {
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
            MET_pred = 0.;
            JetPt_pred->clear();
            JetEta_pred->clear();
            DeltaPhi_pred->clear();
        }
    }

    return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void RandS::calcPredictions(std::vector<myJet>& Jets, LorentzVector& vMET, const int& i, const double& w) {

    int NJets = calcNJets(Jets);
    int NBJets = calcNBJets(Jets);
    double HT = calcHT(Jets);
    LorentzVector vMHT = calcMHT(Jets, JetsMHTPt_, JetsMHTEta_);
    double MHT = vMHT.Pt();
    double MET = vMET.Pt();

    weight = w;
    Ntries_pred = i;
    Njets_pred = NJets;
    BTags_pred = NBJets;
    HT_pred = HT/1000.;
    MHT_pred = MHT/1000.;
    MET_pred = MET/1000.;
    calcLeadingJetPredictions(Jets, vMET);

    if (i == -2) {
        NJetSeed = Njets_pred;
        HTSeed = HT_pred;
    }

    //cout << "weight, Ntries_pred, Njets_pred, BTags_pred, HT_pred, MHT_pred: " << weight << "," << Ntries_pred << "," << Njets_pred << "," << BTags_pred << "," << HT_pred << "," << MHT_pred << endl;

    return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void RandS::calcLeadingJetPredictions(std::vector<myJet>& Jets, LorentzVector& vMHT) {
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

#include <MyAnalysis/MyResolution.h>

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
#include "xAODMuon/Muon.h"

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>

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


std::vector<std::string> getTokens(TString line, TString delim) {
    std::vector<std::string> vtokens;
    TObjArray* tokens = TString(line).Tokenize(delim); //delimiters
    if(tokens->GetEntriesFast()) {
        TIter iString(tokens);
        TObjString* os=0;
        while ((os=(TObjString*)iString())) {
            vtokens.push_back( os->GetString().Data() );
        }
    }
    delete tokens;
    return vtokens;
}

// this is needed to distribute the algorithm to the workers
ClassImp(MyResolution)

MyResolution :: MyResolution ():
    m_VetoCone(1.2),
    m_AddActivityCone(0.8),
    m_MatchingCone(0.1),
    m_RelGenActivityVeto(0.01),
    m_RelRecoActivityVeto(0.05)
{
    // Here you put any code for the base initialization of variables,
    // e.g. initialize all pointers to 0.  Note that you should only put
    // the most basic initialization here, since this method will be
    // called on both the submission and the worker node.  Most of your
    // initialization code will go into histInitialize() and
    // initialize().

    //
    // Property declaration
    //

    PtBinEdges.push_back(0);
    PtBinEdges.push_back(20);
    PtBinEdges.push_back(30);
    PtBinEdges.push_back(50);
    PtBinEdges.push_back(80);
    PtBinEdges.push_back(120);
    PtBinEdges.push_back(170);
    PtBinEdges.push_back(230);
    PtBinEdges.push_back(300);
    PtBinEdges.push_back(380);
    PtBinEdges.push_back(470);
    PtBinEdges.push_back(570);
    PtBinEdges.push_back(680);
    PtBinEdges.push_back(800);
    PtBinEdges.push_back(1000);
    PtBinEdges.push_back(1300);
    PtBinEdges.push_back(1700);
    PtBinEdges.push_back(2200);
    PtBinEdges.push_back(2800);
    PtBinEdges.push_back(3500);
    PtBinEdges.push_back(4300);
    PtBinEdges.push_back(5200);
    PtBinEdges.push_back(6500);

    EtaBinEdges.push_back(0.0);
    EtaBinEdges.push_back(0.7);
    EtaBinEdges.push_back(1.3);
    EtaBinEdges.push_back(1.8);
    EtaBinEdges.push_back(3.2);
    EtaBinEdges.push_back(5.0);

    //// Array of histograms for jet resolutions (all jet multiplicities)
    ResizeHistoVector(PtResolution_tot);
    ResizeHistoVector(EtaResolution_tot);
    ResizeHistoVector(PhiResolution_tot);
    ResizeHistoVector(PtResolution_LF);
    ResizeHistoVector(EtaResolution_LF);
    ResizeHistoVector(PhiResolution_LF);
    ResizeHistoVector(PtResolution_HF);
    ResizeHistoVector(EtaResolution_HF);
    ResizeHistoVector(PhiResolution_HF);
    ResizeHistoVector(PtResolution_nob);
    ResizeHistoVector(EtaResolution_nob);
    ResizeHistoVector(PhiResolution_nob);
    ResizeHistoVector(PtResolution_b);
    ResizeHistoVector(EtaResolution_b);
    ResizeHistoVector(PhiResolution_b);
}



EL::StatusCode MyResolution :: setupJob (EL::Job& job)
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



EL::StatusCode MyResolution :: histInitialize ()
{
    // Here you do everything that needs to be done at the very
    // beginning on each worker node, e.g. create histograms and output
    // trees.  This method gets called before any input files are
    // connected.

    for (unsigned int i_pt = 0; i_pt < PtBinEdges.size() - 1; ++i_pt) {
        for (unsigned int i_eta = 0; i_eta < EtaBinEdges.size() - 1; ++i_eta) {

            //// Book histograms Pt response
            TH1F* h_jetRes_tot_pt = new TH1F(GetHistName(i_pt, i_eta,"tot","Pt").c_str(), "p_T^{reco}/p_T^{gen}", 300, 0.0, 3.0);
            wk()->addOutput(h_jetRes_tot_pt);
            PtResolution_tot.at(i_pt).at(i_eta) = h_jetRes_tot_pt;

            TH1F* h_jetRes_LF_pt = new TH1F(GetHistName(i_pt, i_eta,"LF","Pt").c_str(), "p_T^{reco}/p_T^{gen}", 300, 0.0, 3.0);
            wk()->addOutput(h_jetRes_LF_pt);
            PtResolution_LF.at(i_pt).at(i_eta) = h_jetRes_LF_pt;

            TH1F* h_jetRes_HF_pt = new TH1F(GetHistName(i_pt, i_eta,"HF","Pt").c_str(), "p_T^{reco}/p_T^{gen}", 300, 0.0, 3.0);
            wk()->addOutput(h_jetRes_HF_pt);
            PtResolution_HF.at(i_pt).at(i_eta) = h_jetRes_HF_pt;

            TH1F* h_jetRes_nob_pt = new TH1F(GetHistName(i_pt, i_eta,"nob","Pt").c_str(), "p_T^{reco}/p_T^{gen}", 300, 0.0, 3.0);
            wk()->addOutput(h_jetRes_nob_pt);
            PtResolution_nob.at(i_pt).at(i_eta) = h_jetRes_nob_pt;

            TH1F* h_jetRes_b_pt = new TH1F(GetHistName(i_pt, i_eta,"b","Pt").c_str(), "p_T^{reco}/p_T^{gen}", 300, 0.0, 3.0);
            wk()->addOutput(h_jetRes_b_pt);
            PtResolution_b.at(i_pt).at(i_eta) = h_jetRes_b_pt;

            //// Book histograms Phi resolution
            TH1F* h_jetRes_tot_phi = new TH1F(GetHistName(i_pt, i_eta,"tot","Phi").c_str(), "#Delta#phi(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            wk()->addOutput(h_jetRes_tot_phi);
            PhiResolution_tot.at(i_pt).at(i_eta) = h_jetRes_tot_phi;

            TH1F* h_jetRes_LF_phi = new TH1F(GetHistName(i_pt, i_eta,"LF","Phi").c_str(), "#Delta#phi(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            wk()->addOutput(h_jetRes_LF_phi);
            PhiResolution_LF.at(i_pt).at(i_eta) = h_jetRes_LF_phi;

            TH1F* h_jetRes_HF_phi = new TH1F(GetHistName(i_pt, i_eta,"HF","Phi").c_str(), "#Delta#phi(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            wk()->addOutput(h_jetRes_HF_phi);
            PhiResolution_HF.at(i_pt).at(i_eta) = h_jetRes_HF_phi;

            TH1F* h_jetRes_nob_phi = new TH1F(GetHistName(i_pt, i_eta,"nob","Phi").c_str(), "#Delta#phi(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            wk()->addOutput(h_jetRes_nob_phi);
            PhiResolution_nob.at(i_pt).at(i_eta) = h_jetRes_nob_phi;

            TH1F* h_jetRes_b_phi = new TH1F(GetHistName(i_pt, i_eta,"b","Phi").c_str(), "#Delta#phi(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            wk()->addOutput(h_jetRes_b_phi);
            PhiResolution_b.at(i_pt).at(i_eta) = h_jetRes_b_phi;

            //// Book histograms Eta resolution
            TH1F* h_jetRes_tot_eta = new TH1F(GetHistName(i_pt, i_eta,"tot","Eta").c_str(), "#Delta#eta(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            wk()->addOutput(h_jetRes_tot_eta);
            EtaResolution_tot.at(i_pt).at(i_eta) = h_jetRes_tot_eta;

            TH1F* h_jetRes_LF_eta = new TH1F(GetHistName(i_pt, i_eta,"LF","Eta").c_str(), "#Delta#eta(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            wk()->addOutput(h_jetRes_LF_eta);
            EtaResolution_LF.at(i_pt).at(i_eta) = h_jetRes_LF_eta;

            TH1F* h_jetRes_HF_eta = new TH1F(GetHistName(i_pt, i_eta,"HF","Eta").c_str(), "#Delta#eta(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            wk()->addOutput(h_jetRes_HF_eta);
            EtaResolution_HF.at(i_pt).at(i_eta) = h_jetRes_HF_eta;

            TH1F* h_jetRes_nob_eta = new TH1F(GetHistName(i_pt, i_eta,"nob","Eta").c_str(), "#Delta#eta(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            wk()->addOutput(h_jetRes_nob_eta);
            EtaResolution_nob.at(i_pt).at(i_eta) = h_jetRes_nob_eta;

            TH1F* h_jetRes_b_eta = new TH1F(GetHistName(i_pt, i_eta,"b","Eta").c_str(), "#Delta#eta(jet^{reco},jet^{gen})", 80, -0.2, 0.2);
            wk()->addOutput(h_jetRes_b_eta);
            EtaResolution_b.at(i_pt).at(i_eta) = h_jetRes_b_eta;
        }
    }

    for (unsigned int i_eta = 0; i_eta < EtaBinEdges.size() - 1; ++i_eta) {

        //// Book histograms for jet counts (needed for reconstruction efficiency)
        TH1F* h_NjetReco_tot_pt = new TH1F(GetHistNameEta(i_eta,"h_NjetReco_tot").c_str(), "N(jet^{reco})", 100, 0., 500000.);
        wk()->addOutput(h_NjetReco_tot_pt);
        NReco_tot.push_back(h_NjetReco_tot_pt);

        TH1F* h_NjetReco_b_pt = new TH1F(GetHistNameEta(i_eta,"h_NjetReco_b").c_str(), "N(jet^{reco})", 100, 0., 500000.);
        wk()->addOutput(h_NjetReco_b_pt);
        NReco_b.push_back(h_NjetReco_b_pt);

        TH1F* h_NjetReco_nob_pt = new TH1F(GetHistNameEta(i_eta,"h_NjetReco_nob").c_str(), "N(jet^{reco})", 100, 0., 500000.);
        wk()->addOutput(h_NjetReco_nob_pt);
        NReco_nob.push_back(h_NjetReco_nob_pt);

        TH1F* h_NjetGen_tot_pt = new TH1F(GetHistNameEta(i_eta,"h_NjetGen_tot").c_str(), "N(jet^{gen})", 100, 0., 500000.);
        wk()->addOutput(h_NjetGen_tot_pt);
        NGen_tot.push_back(h_NjetGen_tot_pt);

        TH1F* h_NjetGen_b_pt = new TH1F(GetHistNameEta(i_eta,"h_NjetGen_b").c_str(), "N(jet^{gen})", 100, 0., 500000.);
        wk()->addOutput(h_NjetGen_b_pt);
        NGen_b.push_back(h_NjetGen_b_pt);

        TH1F* h_NjetGen_nob_pt = new TH1F(GetHistNameEta(i_eta,"h_NjetGen_nob").c_str(), "N(jet^{gen})", 100, 0., 500000.);
        wk()->addOutput(h_NjetGen_nob_pt);
        NGen_nob.push_back(h_NjetGen_nob_pt);
    }

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyResolution :: fileExecute ()
{
    // Here you do everything that needs to be done exactly once for every
    // single file, e.g. collect a list of all lumi-blocks processed
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyResolution :: changeInput (bool firstFile)
{
    // Here you do everything you need to do when we change input files,
    // e.g. resetting branch addresses on trees.  If you are using
    // D3PDReader or a similar service this method is not needed.
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyResolution :: initialize ()
{
    // Here you do everything that you need to do after the first input
    // file has been connected and before the first event is processed,
    // e.g. create additional histograms based on which variables are
    // available in the input files.  You can also create all of your
    // histograms and trees in here, but be aware that this method
    // doesn't get called if no events are processed.  So any objects
    // you create here won't be available in the output if you have no
    // input events.

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
        prw_conf = getTokens(prw_file_,",");
        //prw_conf.push_back(prw_file_);
    }
    EL_RETURN_CHECK("initialize", objTool->setProperty("PRWConfigFiles", prw_conf) );

    std::vector<std::string> prw_lumicalc;
    ilumicalc_file_ = "DUMMY";
    if (ilumicalc_file_ == "DUMMY") {
        prw_lumicalc.push_back("/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/SUSYTools/ilumicalc_histograms_None_276262-284154_IBLOFF.root");
        prw_lumicalc.push_back("/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/SUSYTools/ilumicalc_histograms_None_297730-299243.root");

    } else {
        prw_lumicalc = getTokens(ilumicalc_file_, ",");
        //prw_lumicalc.push_back(ilumicalc_file_);
    }
    EL_RETURN_CHECK("initialize", objTool->setProperty("PRWLumiCalcFiles", prw_lumicalc) );

    ///////////////////////////////////////////////////////////////////////////////////////////
    // Configure the SUSYObjDef instance
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



EL::StatusCode MyResolution :: execute ()
{
    // Here you do everything that needs to be done on every single
    // events, e.g. read input variables, apply cuts, and fill
    // histograms and trees.  This is where most of your actual analysis
    // code will go.

    xAOD::TEvent* event = wk()->xaodEvent();

    // print every 100 events, so we know where we are:
    if( (m_eventCounter % 100) ==0 ) Info("execute()", "Event number = %i", m_eventCounter );
    m_eventCounter++;

    //----------------------------
    // Event information
    //---------------------------
    const xAOD::EventInfo* eventInfo = 0;
    EL_RETURN_CHECK("execute", event->retrieve( eventInfo, "EventInfo"));

    EL_RETURN_CHECK("execute", objTool->ApplyPRWTool() );


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
    xAOD::MuonContainer* muons_nominal(0);
    xAOD::ShallowAuxContainer* muons_nominal_aux(0);
    EL_RETURN_CHECK("execute()", objTool->GetMuons(muons_nominal, muons_nominal_aux) );

    // Jets
    xAOD::JetContainer* jets(0);
    xAOD::ShallowAuxContainer* jets_aux(0);
    EL_RETURN_CHECK("execute()", objTool->GetJets(jets, jets_aux) );

    // Taus
    //xAOD::TauJetContainer* taus_nominal(0);
    //xAOD::ShallowAuxContainer* taus_nominal_aux(0);
    //EL_RETURN_CHECK("execute()", objTool->GetTaus(taus_nominal,taus_nominal_aux) );

    // MET
    //xAOD::MissingETContainer* mettst_nominal = new xAOD::MissingETContainer;
    //xAOD::MissingETAuxContainer* mettst_nominal_aux = new xAOD::MissingETAuxContainer;
    //mettst_nominal->setStore(mettst_nominal_aux);
    //mettst_nominal->reserve(10);


    // get jet container of interest
    //const xAOD::JetContainer* jets = 0;
    //EL_RETURN_CHECK("execute()", event->retrieve( jets, "AntiKt4EMTopoJets" ));
    const xAOD::JetContainer* genjets = 0;
    EL_RETURN_CHECK("execute()", event->retrieve( genjets, "AntiKt4TruthJets"));
    const xAOD::TruthParticleContainer* genparticles = 0;
    EL_RETURN_CHECK("execute()", event->retrieve( genparticles, "TruthParticles"));

    // loop over the emjets and reject events with at least one bad jet
    xAOD::JetContainer::const_iterator jet_itr = jets->begin();
    xAOD::JetContainer::const_iterator jet_end = jets->end();
    for( ; jet_itr != jet_end; ++jet_itr ) {
        if( !m_jetCleaning->accept( **jet_itr ) && (*jet_itr)->pt() > 20000.) {
            Info("execute()", "Reject event because of a bad jet");
            return EL::StatusCode::SUCCESS;
        }
    }

    // Print out number of jets
    //Info("execute()", "  number of jets = %lu", jets->size());

    // Loop over all jets in this container
    for ( const auto& genjet : *genjets ) {

        TLorentzVector numuActivity(0., 0., 0., 0.);

        for ( const auto& it : *genparticles ) {
            if (abs(it->pdgId()) == 13 || abs(it->pdgId()) == 12 || abs(it->pdgId()) == 14 || abs(it->pdgId()) == 16 ) {
                if (it->pt() > 0.0 && it->status() == 1 && !(it->hasDecayVtx())) {
                    //std::cout << "id, status, pt, eta, phi: " << it->pdgId() << ", " << it->status() << ", " << it->pt() << ", "<< it->eta()<< ", " << it->phi()<< std::endl;
                    double dR = genjet->p4().DeltaR(it->p4());
                    if (dR < 0.4) {
                        numuActivity += it->p4();
                    }
                }
            }
        }
        //if (numuActivity.Pt() > 0) std::cout << "Old, new pt: " <<  genjet->pt() << ", " << (genjet->p4() + numuActivity).Pt() << std::endl;

        // remove jets in very forward region
        if (fabs(genjet->eta()) > 4.5) continue;

        // check for no additional genJet activity
        bool noGenActivity = true;
        for ( const auto& genjet2 : *genjets ) {
            if (genjet2 == genjet) continue;
            double dR = genjet->p4().DeltaR(genjet2->p4());
            if (dR < m_VetoCone && genjet2->pt()/genjet->pt() > m_RelGenActivityVeto ) {
                noGenActivity = false;
            }
        }

        //if (!noGenActivity) continue; // continue with next genJet if another genJet is closeby

        // check for additional recoJet activity
        const xAOD::Jet* matchedJet = 0; //create a new jet object
        const xAOD::Jet* nextJet = 0; //create a new jet object
        TLorentzVector addRecoActivity(0., 0., 0., 0.);
        double dRmin_matched = 999.;
        double dRmin_next = 999.;
        for ( const auto& jet : *jets) {
            double dR = jet->p4().DeltaR(genjet->p4());
            if (dR < dRmin_matched && dR < m_VetoCone) {
                if (dRmin_matched < m_VetoCone ) {
                    dRmin_next = dRmin_matched;
                    nextJet = matchedJet;
                    if (dR < m_AddActivityCone) addRecoActivity += nextJet->p4();
                }
                dRmin_matched = dR;
                matchedJet = jet;
            } else if (dR < dRmin_next && dR < m_VetoCone) {
                dRmin_next = dR;
                nextJet = jet;
                if (dR < m_AddActivityCone) addRecoActivity += nextJet->p4();
            }
        } // end for loop over jets

        bool noRecoActivity = true;
        TLorentzVector MuonContribution(0., 0., 0., 0.);
        if (dRmin_matched < m_MatchingCone) {
            if (dRmin_next < m_VetoCone) {
                if (nextJet->pt()/matchedJet->pt()>m_RelRecoActivityVeto) noRecoActivity = false;
            }
            //for ( const auto& muon : *muons_nominal) {
            //    if (matchedJet->p4().DeltaR(muon->p4()) < 0.4) {
            //        MuonContribution += muon->p4();
            //    }
            //}
        }

        if (noGenActivity) {
            int ii_eta = GetEtaBin(genjet->eta());
            NGen_tot.at(ii_eta)->Fill(genjet->pt(), eventWeight);
            //if (objTool->IsTruthBJet(*genjet)) {
            if (abs(genjet->getAttribute<int>("PartonTruthLabelID")) == 5) {
                NGen_b.at(ii_eta)->Fill(genjet->pt(), eventWeight);
            } else {
                NGen_nob.at(ii_eta)->Fill(genjet->pt(), eventWeight);
            }
        }

        if (noGenActivity && dRmin_matched < 0.2) {
            int ii_eta = GetEtaBin(genjet->eta());
            NReco_tot.at(ii_eta)->Fill(genjet->pt(), eventWeight);
            //if (objTool->IsTruthBJet(*genjet)) {
            if (abs(genjet->getAttribute<int>("PartonTruthLabelID")) == 5) {
                NReco_b.at(ii_eta)->Fill(genjet->pt(), eventWeight);
            } else {
                NReco_nob.at(ii_eta)->Fill(genjet->pt(), eventWeight);
            }
        }

        if (dRmin_matched < m_MatchingCone && noRecoActivity && noGenActivity) {
            //if (dRmin_matched < m_MatchingCone && noGenActivity) {
            int i_pt = GetPtBin((genjet->p4()+numuActivity).E()/1000.); // Reminder here E is used instead pT for binning (variabel name not changed yet, maybe later)
            int i_eta = GetEtaBin(genjet->eta());
            PtResolution_tot.at(i_pt).at(i_eta)->Fill((matchedJet->p4()+addRecoActivity+MuonContribution).Pt()/(genjet->p4()+numuActivity).Pt(), eventWeight);
            PhiResolution_tot.at(i_pt).at(i_eta)->Fill(matchedJet->p4().DeltaPhi(genjet->p4()), eventWeight);
            EtaResolution_tot.at(i_pt).at(i_eta)->Fill(matchedJet->eta()-genjet->eta(), eventWeight);
            if (objTool->IsBJet(*matchedJet)) {
                PtResolution_b.at(i_pt).at(i_eta)->Fill((matchedJet->p4()+addRecoActivity+MuonContribution).Pt()/(genjet->p4()+numuActivity).Pt(), eventWeight);
                PhiResolution_b.at(i_pt).at(i_eta)->Fill(matchedJet->p4().DeltaPhi(genjet->p4()), eventWeight);
                EtaResolution_b.at(i_pt).at(i_eta)->Fill(matchedJet->eta()-genjet->eta(), eventWeight);
            } else {
                PtResolution_nob.at(i_pt).at(i_eta)->Fill((matchedJet->p4()+addRecoActivity+MuonContribution).Pt()/(genjet->p4()+numuActivity).Pt(), eventWeight);
                PhiResolution_nob.at(i_pt).at(i_eta)->Fill(matchedJet->p4().DeltaPhi(genjet->p4()), eventWeight);
                EtaResolution_nob.at(i_pt).at(i_eta)->Fill(matchedJet->eta()-genjet->eta(), eventWeight);
            }
            //if (objTool->IsTruthBJet(*genjet)) {
            if (abs(genjet->getAttribute<int>("PartonTruthLabelID")) == 5) {
                PtResolution_HF.at(i_pt).at(i_eta)->Fill((matchedJet->p4()+addRecoActivity+MuonContribution).Pt()/(genjet->p4()+numuActivity).Pt(), eventWeight);
                PhiResolution_HF.at(i_pt).at(i_eta)->Fill(matchedJet->p4().DeltaPhi(genjet->p4()), eventWeight);
                EtaResolution_HF.at(i_pt).at(i_eta)->Fill(matchedJet->eta()-genjet->eta(), eventWeight);
            } else {
                PtResolution_LF.at(i_pt).at(i_eta)->Fill((matchedJet->p4()+addRecoActivity+MuonContribution).Pt()/(genjet->p4()+numuActivity).Pt(), eventWeight);
                PhiResolution_LF.at(i_pt).at(i_eta)->Fill(matchedJet->p4().DeltaPhi(genjet->p4()), eventWeight);
                EtaResolution_LF.at(i_pt).at(i_eta)->Fill(matchedJet->eta()-genjet->eta(), eventWeight);
            }
        }

        //delete matchedJet;
        //delete nextJet;

    } // end for loop over genjets

    // The containers created by the shallow copy are owned by you. Remember to delete them
    delete jets; // not these, we put them in the store!
    delete jets_aux;
    // delete taus_nominal;
    // delete taus_nominal_aux;
    delete muons_nominal;
    delete muons_nominal_aux;
    // delete electrons_nominal;
    // delete electrons_nominal_aux;
    // delete photons_nominal;
    // delete photons_nominal_aux;
    // delete metcst_nominal;
    // delete metcst_nominal_aux;
    // delete mettst_nominal;
    // delete mettst_nominal_aux;

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyResolution :: postExecute ()
{
    // Here you do everything that needs to be done after the main event
    // processing.  This is typically very rare, particularly in user
    // code.  It is mainly used in implementing the NTupleSvc.
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyResolution :: finalize ()
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



EL::StatusCode MyResolution :: histFinalize ()
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

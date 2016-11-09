#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <MyAnalysis/NtupleMaker.h>

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
#define EL_RETURN_CHECK( CONTEXT, EXP )                         \
	   do {                                                     \
	      if( ! EXP.isSuccess() ) {                             \
	         Error( CONTEXT,                                    \
	                XAOD_MESSAGE( "Failed to execute: %s" ),    \
	                #EXP );                                     \
	         return EL::StatusCode::FAILURE;                    \
	      }                                                     \
	   } while( false )
#endif


std::vector<std::string> getTokens2(TString line, TString delim) {
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
ClassImp(NtupleMaker)

NtupleMaker :: NtupleMaker ()
{
    debug_ = 0;
    outputfile_ = "NtupleMaker.root";
    calculateTrueMET_ = false;
    rebalancedJetPt_ = 20000.;
    AddMuToJets_ = true;

    // Here you put any code for the base initialization of variables,
    // e.g. initialize all pointers to 0.  Note that you should only put
    // the most basic initialization here, since this method will be
    // called on both the submission and the worker node.  Most of your
    // initialization code will go into histInitialize() and
    // initialize().

}



EL::StatusCode NtupleMaker :: setupJob (EL::Job& job)
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



EL::StatusCode NtupleMaker :: histInitialize ()
{
    // Here you do everything that needs to be done at the very
    // beginning on each worker node, e.g. create histograms and output
    // trees.  This method gets called before any input files are
    // connected.

    // get object of class SmearFunction
    // I assume there is a better way of passing config parameters, but for now it will be like this

    // define output tree
    std::cout << "outputfile_: " << outputfile_ << std::endl;
    TFile *outputFile = wk()->getOutputFile(outputfile_);
    EventTree = new TTree("EventTree", "EventTree");
    EventTree->SetDirectory(outputFile);

    std::cout << EventTree << std::endl;
    EventTree->SetAutoSave(10000000000);
    EventTree->SetAutoFlush(100000000);

    // set branches for output tree

    //EventTree->Branch("NVtx", &NVtx_);
    EventTree->Branch("Weight",&weight_);
    EventTree->Branch("DatasetID",&dsid_);
    EventTree->Branch("PrimaryVtx",&pvtx_);

    EventTree->Branch("JetPt", "std::vector<Float_t>", &JetPt_n);
    EventTree->Branch("JetEta", "std::vector<Float_t>", &JetEta_n);
    EventTree->Branch("JetPhi", "std::vector<Float_t>", &JetPhi_n);
    EventTree->Branch("JetM", "std::vector<Float_t>", &JetM_n);
    EventTree->Branch("JetBtag", "std::vector<bool>", &JetBtag_n);
    EventTree->Branch("JetJVT", "std::vector<Float_t>", &JetJVT_n);
    EventTree->Branch("JetGood", "std::vector<bool>", &JetGood_n);

    EventTree->Branch("JetNoMuPt", "std::vector<Float_t>", &JetPt_n);
    EventTree->Branch("JetNoMuEta", "std::vector<Float_t>", &JetEta_n);
    EventTree->Branch("JetNoMuPhi", "std::vector<Float_t>", &JetPhi_n);
    EventTree->Branch("JetNoMuM", "std::vector<Float_t>", &JetM_n);
    EventTree->Branch("JetNoMuBtag", "std::vector<bool>", &JetBtag_n);
    EventTree->Branch("JetNoMuJVT", "std::vector<Float_t>", &JetJVT_n);
    EventTree->Branch("JetNoMuGood", "std::vector<bool>", &JetGood_n);

    EventTree->Branch("GenJetPt","std::vector<Float_t>", &GenJetPt_n);
    EventTree->Branch("GenJetEta","std::vector<Float_t>", &GenJetEta_n);
    EventTree->Branch("GenJetPhi","std::vector<Float_t>", &GenJetPhi_n);
    EventTree->Branch("GenJetM","std::vector<Float_t>", &GenJetM_n);
    EventTree->Branch("GenJetBtag","std::vector<bool>", &GenJetBtag_n);

    //EventTree->Branch("GenJetNoNuPt","std::vector<Float_t>", &GenJetNoNuPt_n);
    //EventTree->Branch("GenJetNoNuEta","std::vector<Float_t>", &GenJetNoNuEta_n);
    //EventTree->Branch("GenJetNoNuPhi","std::vector<Float_t>", &GenJetNoNuPhi_n);
    //EventTree->Branch("GenJetNoNuM","std::vector<Float_t>", &GenJetNoNuM_n);
    //EventTree->Branch("GenJetNoNuBtag","std::vector<bool>", &GenJetNoNuBtag_n);

    EventTree->Branch("GenJetNoNuMuPt","std::vector<Float_t>", &GenJetNoNuMuPt_n);
    EventTree->Branch("GenJetNoNuMuEta","std::vector<Float_t>", &GenJetNoNuMuEta_n);
    EventTree->Branch("GenJetNoNuMuPhi","std::vector<Float_t>", &GenJetNoNuMuPhi_n);
    EventTree->Branch("GenJetNoNuMuM","std::vector<Float_t>", &GenJetNoNuMuM_n);
    EventTree->Branch("GenJetNoNuMuBtag","std::vector<bool>", &GenJetNoNuMuBtag_n);

    EventTree->Branch("ElePt","std::vector<Float_t>", &ElePt_n);
    EventTree->Branch("EleEta","std::vector<Float_t>", &EleEta_n);
    EventTree->Branch("ElePhi","std::vector<Float_t>", &ElePhi_n);

    EventTree->Branch("PhotonPt","std::vector<Float_t>", &PhotonPt_n);
    EventTree->Branch("PhotonEta","std::vector<Float_t>", &PhotonEta_n);
    EventTree->Branch("PhotonPhi","std::vector<Float_t>", &PhotonPhi_n);

    EventTree->Branch("MuonPt","std::vector<Float_t>", &MuonPt_n);
    EventTree->Branch("MuonEta","std::vector<Float_t>", &MuonEta_n);
    EventTree->Branch("MuonPhi","std::vector<Float_t>", &MuonPhi_n);
    EventTree->Branch("MuonIsBad","std::vector<bool>", &MuonIsBad_n);
    EventTree->Branch("MuonIsSignal","std::vector<bool>", &MuonIsSignal_n);
    EventTree->Branch("MuonIsIso","std::vector<bool>", &MuonIsIso_n);


    EventTree->Branch("MET_pt", &MET_pt_);
    EventTree->Branch("MET_phi", &MET_phi_);
    EventTree->Branch("METjet_pt", &METjet_pt_);
    EventTree->Branch("METjet_phi", &METjet_phi_);
    EventTree->Branch("METmu_pt", &METmu_pt_);
    EventTree->Branch("METmu_phi", &METmu_phi_);
    EventTree->Branch("METele_pt", &METele_pt_);
    EventTree->Branch("METele_phi", &METele_phi_);
    EventTree->Branch("METgamma_pt", &METgamma_pt_);
    EventTree->Branch("METgamma_phi", &METgamma_phi_);
    EventTree->Branch("METtrack_pt", &METtrack_pt_);
    EventTree->Branch("METtrack_phi", &METtrack_phi_);
    EventTree->Branch("GenMET_pt", &GenMET_pt_);
    EventTree->Branch("GenMET_phi", &GenMET_phi_);
    EventTree->Branch("TrueMHT_pt", &TrueMHT_pt_);
    EventTree->Branch("TrueMHT_phi", &TrueMHT_phi_);


    return EL::StatusCode::SUCCESS;
}



EL::StatusCode NtupleMaker :: fileExecute ()
{
    // Here you do everything that needs to be done exactly once for every
    // single file, e.g. collect a list of all lumi-blocks processed
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode NtupleMaker :: changeInput (bool firstFile)
{
    // Here you do everything you need to do when we change input files,
    // e.g. resetting branch addresses on trees.  If you are using
    // D3PDReader or a similar service this method is not needed.
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode NtupleMaker :: initialize ()
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


    isMC = false;
    const xAOD::EventInfo* eventInfo = 0;
    EL_RETURN_CHECK("initialize",event->retrieve( eventInfo, "EventInfo"));
    if (eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ) {
        isMC = true;
    }
    ST::ISUSYObjDef_xAODTool::DataSource datasource = ( isMC ? ST::ISUSYObjDef_xAODTool::FullSim : ST::ISUSYObjDef_xAODTool::Data );

    objTool = new ST::SUSYObjDef_xAOD("SUSYObjDef_xAOD");

    prw_file_ = "DUMMY";
    std::vector<std::string> prw_conf;
    if (prw_file_ == "DUMMY") {
        prw_conf.push_back("/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/PileupReweighting/mc15ab_defaults.NotRecommended.prw.root");
        prw_conf.push_back("/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root");
    }
    else {
        prw_conf = getTokens2(prw_file_,",");
        //prw_conf.push_back(prw_file_);
    }
    EL_RETURN_CHECK("initialize", objTool->setProperty("PRWConfigFiles", prw_conf) );

    std::vector<std::string> prw_lumicalc;
    ilumicalc_file_ = "DUMMY";
    if (ilumicalc_file_ == "DUMMY") {
        prw_lumicalc.push_back("/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/SUSYTools/ilumicalc_histograms_None_276262-284154_IBLOFF.root");
        prw_lumicalc.push_back("/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/SUSYTools/ilumicalc_histograms_None_297730-299243.root");

    } else {
        prw_lumicalc = getTokens2(ilumicalc_file_, ",");
        //prw_lumicalc.push_back(ilumicalc_file_);
    }
    EL_RETURN_CHECK("initialize", objTool->setProperty("PRWLumiCalcFiles", prw_lumicalc) );

    ///////////////////////////////////////////////////////////////////////////////////////////
    // Configure the SUSYObjDef instance
    EL_RETURN_CHECK("initialize", objTool->setProperty("DataSource", datasource) ) ;
    EL_RETURN_CHECK("initialize", objTool->setProperty("ConfigFile", "SUSYTools/SUSYTools_Default.conf") );

    //Manually setting additional properties will override what's in the configuration file
    //EL_RETURN_CHECK("initialize", objTool->setProperty("EleId", "TightLLH") );
    //EL_RETURN_CHECK("initialize", objTool->setProperty("EleBaselineId", "LooseLLH") );

    EL_RETURN_CHECK("initialize", objTool->initialize() );

    // SUSYTools cross section data base
    my_XsecDB = new SUSY::CrossSectionDB(gSystem->ExpandPathName("$ROOTCOREBIN/data/SUSYTools/mc15_13TeV/"));

    // GRL
    m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    const char* grlFilePath = "$ROOTCOREBIN/data/MyAnalysis/data16_13TeV.periodAllYear_DetStatus-v83-pro20-15_DQDefects-00-02-04_PHYS_StandardGRL_All_Good_25ns.xml";
    const char* fullGRLFilePath = gSystem->ExpandPathName (grlFilePath);
    std::vector<std::string> vecStringGRL;
    vecStringGRL.push_back(fullGRLFilePath);
    EL_RETURN_CHECK("initialize()",m_grl->setProperty( "GoodRunsListVec", vecStringGRL));
    EL_RETURN_CHECK("initialize()",m_grl->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
    EL_RETURN_CHECK("initialize()",m_grl->initialize());

    // Initialize and configure trigger tools
    m_trigConfigTool = new TrigConf::xAODConfigTool("xAODConfigTool"); // gives us access to the meta-data
    EL_RETURN_CHECK("initialize",  m_trigConfigTool->initialize() );
    ToolHandle< TrigConf::ITrigConfigTool > trigConfigHandle( m_trigConfigTool );
    m_trigDecisionTool = new Trig::TrigDecisionTool("TrigDecisionTool");
    EL_RETURN_CHECK("initialize", m_trigDecisionTool->setProperty( "ConfigTool", trigConfigHandle ) ); // connect the TrigDecisionTool to the ConfigTool
    EL_RETURN_CHECK("initialize", m_trigDecisionTool->setProperty( "TrigDecisionKey", "xTrigDecision" ) );
    EL_RETURN_CHECK("initialize", m_trigDecisionTool->initialize() );

    // as a check, let's see the number of events in our xAOD
    Info("initialize()", "Number of events = %lli", event->getEntries() ); // print long long int

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode NtupleMaker :: execute ()
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
    EL_RETURN_CHECK("execute",event->retrieve( eventInfo, "EventInfo"));

    // check if the event is data or MC
    // (many tools are applied either to data or MC)
    bool isMC = false;
    // check if the event is MC
    double eventWeight = 1;

    if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ) {

        EL_RETURN_CHECK("execute", objTool->ApplyPRWTool() );
        isMC = true; // can do something with this later
        //extra event-level information you might need:
        dsid_ =  eventInfo->mcChannelNumber();
        const std::vector< float > weights = eventInfo->mcEventWeights();
        if( weights.size() > 0 ) eventWeight = weights[0];
        double XS =  my_XsecDB->rawxsect(dsid_);
        double eff = my_XsecDB->efficiency(dsid_);
        double k = my_XsecDB->kfactor(dsid_);
        double weight = XS * k * eff;
        eventWeight *= weight;
        //std::cout << "name: " << my_XsecDB->name(datasetID) << std::endl;
        //std::cout << "XStimesEff: " << my_XsecDB->xsectTimesEff(datasetID) << std::endl;
        //std::cout << "XS: " << my_XsecDB->rawxsect(datasetID) << std::endl;
        //std::cout << "kFactor: " << my_XsecDB->kfactor(datasetID) << std::endl;
        //std::cout << "efficiency: " << my_XsecDB->efficiency(datasetID) << std::endl;

        //std::cout << "eventWeight = " << eventWeight << std::endl;
    }
    weight_ = eventWeight;

    if(!isMC) { // it's data!

		// No dataset ID for data
		dsid_ = eventInfo->runNumber();;

        // check if event passes GRL
        if(!m_grl->passRunLB(*eventInfo)) {
            return EL::StatusCode::SUCCESS;
        }

        //// Get list of single jet triggers

        auto chainGroup = m_trigDecisionTool->getChainGroup("HLT_j.*");
        std::map<int, std::string> PresenTriggers;
        for (auto &trig : chainGroup->getListOfTriggers()) {
            std::string thisTrig = trig;
            if (thisTrig.length() < 9) {
                //Info( "execute()", thisTrig.c_str() );
                //// Extract threshold
                std::string str_thres = thisTrig.substr(5);
                //std::cout << str_thres << " to " << atoi(str_thres.c_str()) << std::endl;
                PresenTriggers[atoi(str_thres.c_str())] = thisTrig;
            }
        }
        std::vector<std::string> single_jet_triggers;
        std::vector<double> single_jet_theshold;
        for (std::map<int, std::string>::iterator it = PresenTriggers.begin(); it != PresenTriggers.end(); it++) {
            single_jet_triggers.push_back(it->second);
            single_jet_theshold.push_back(int(it->first));
        }

        //// calculate the prescale weight from single jet triggers

        std::vector<bool > triggersPass;
        for(const std::string &sjt : single_jet_triggers) {
            //// use susytools or your favourite tool to test if each single jet trigger has fired
            triggersPass.push_back(objTool->IsTrigPassed(sjt));
        }

        //// retrieve the HLT jets
        const xAOD::JetContainer * hlt_jet = 0;
        TString mc_name="HLT_xAOD__JetContainer_a4tcemsubjesFS";
        if( ! event->retrieve( hlt_jet, mc_name.Data()).isSuccess() )
            Error("execute()", Form("failed to retrieve %s", mc_name.Data()));

        double HLTjetPt = 0;
        if (hlt_jet->size() > 0) {
            HLTjetPt = hlt_jet->at(0)->pt()/1000.;
            //Info( "execute()", "Leading jet pT: %.1f ", HLTjetPt);
        }

        int i_highest_fire = -1;
        for (int i  = triggersPass.size()-1; i > 0; i--) {
            //std::cout << "fire i: " << i << std::endl;
            if (triggersPass.at(i)) {
                i_highest_fire = i;
                break;
            }
        }

        int i_highest_threshold = -1;
        for (int i  = single_jet_theshold.size()-1; i > 0; i--) {
            //std::cout << "threshold i: " << i << std::endl;
            if (single_jet_theshold.at(i) < HLTjetPt) {
                i_highest_threshold = i;
                break;
            }
        }

        //// if another (higher) trigger could have fired -> event weight = 0;
        if (i_highest_fire == -1) {
            eventWeight = 0.;
            Info( "execute()", "No relevant single jet trigger has fired" );
            return EL::StatusCode::SUCCESS;
        } else if ( i_highest_fire < i_highest_threshold) {
            eventWeight = 0.;
            //Info( "execute()", "%30s chain passed, but %.1f threshold could be passed", single_jet_triggers.at(i_highest_fire).c_str(), single_jet_theshold.at(i_highest_threshold));
            return EL::StatusCode::SUCCESS;
        } else {
            //// else event weight = prescale
            //// examine the HLT_jxxx.* chains, see if they passed/failed and their total prescale
            auto chainGroup = m_trigDecisionTool->getChainGroup(single_jet_triggers.at(i_highest_fire));
            std::map<std::string,int> triggerCounts;
            for (auto &trig : chainGroup->getListOfTriggers()) {
                auto cg = m_trigDecisionTool->getChainGroup(trig);
                std::string thisTrig = trig;
                //Info( "execute()", "%30s chain passed(1)/failed(0): %d ; total chain prescale (L1*HLT): %.1f", thisTrig.c_str(), cg->isPassed(), cg->getPrescale() );
                eventWeight = cg->getPrescale();
            }
        }
        
        weight_ = eventWeight;

    } // end if not MC

    pvtx_ = true;
    if (objTool->GetPrimVtx() == nullptr) {
        Info("execute()" , "No PV found for this event! ...");
        pvtx_ = false;
    }

    //------------------------------------------------------------
    // Apply event cleaning to remove events due to
    // problematic regions of the detector
    // or incomplete events.
    // Apply to data.
    //------------------------------------------------------------
    // reject event if:
    if(!isMC) {
        if( (eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::SCT)==xAOD::EventInfo::Error ) || (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) ) )
        {
            return EL::StatusCode::SUCCESS; // go to the next event
        } // end if event flags check
    } // end if the event is data
    m_numCleanEvents++;

    // Electrons
    xAOD::ElectronContainer* electrons_nominal(0);
    xAOD::ShallowAuxContainer* electrons_nominal_aux(0);
    EL_RETURN_CHECK("execute()", objTool->GetElectrons(electrons_nominal, electrons_nominal_aux) );

    // Photons
    xAOD::PhotonContainer* photons_nominal(0);
    xAOD::ShallowAuxContainer* photons_nominal_aux(0);
    EL_RETURN_CHECK("execute()", objTool->GetPhotons(photons_nominal,photons_nominal_aux) );

    // Muons
    xAOD::MuonContainer* muons_nominal(0);
    xAOD::ShallowAuxContainer* muons_nominal_aux(0);
    EL_RETURN_CHECK("execute()", objTool->GetMuons(muons_nominal, muons_nominal_aux) );

    // Jets
    xAOD::JetContainer* jets_nominal(0);
    xAOD::ShallowAuxContainer* jets_nominal_aux(0);
    EL_RETURN_CHECK("execute()", objTool->GetJets(jets_nominal, jets_nominal_aux) );

    // MET
    xAOD::MissingETContainer* mettst_nominal = new xAOD::MissingETContainer;
    xAOD::MissingETAuxContainer* mettst_nominal_aux = new xAOD::MissingETAuxContainer;
    mettst_nominal->setStore(mettst_nominal_aux);
    mettst_nominal->reserve(10);

    // Call Overlap Removal
    //EL_RETURN_CHECK("execute()", objTool->OverlapRemoval(electrons_nominal, muons_nominal, jets_nominal) );

    //// muon vector
    //std::vector<myMuon> Muons;
    float mu_ptcut = 10000.;
    float mu_d0sigcut = 3.;
    float mu_z0cut = 0.5;
    float mu_etacut = 2.7;
    for (const auto& mu : *muons_nominal) {
        if (!(objTool->IsBaselineMuon(*mu))) continue;
        myMuon newMuon;
        newMuon.momentum = mu->p4();
        newMuon.isSignal = false;
        if (objTool->IsNoIsoSignalMuon(*mu, mu_ptcut, mu_d0sigcut, mu_z0cut, mu_etacut)) {
            newMuon.isSignal = true;
        }
        newMuon.isIso = false;
        if (objTool->IsSignalMuon(*mu, mu_ptcut, mu_d0sigcut, mu_z0cut, mu_etacut)) {
            newMuon.isIso = true;
        }
        newMuon.isBad = false;
        if (objTool->IsBadMuon(*mu, 0.2)) {
            newMuon.isBad = true;
        }
        if (newMuon.momentum.Pt() > mu_ptcut && fabs (newMuon.momentum.Eta()) < mu_etacut) {
            MuonPt_n->push_back(newMuon.momentum.Pt()/1000.);
            MuonEta_n->push_back(newMuon.momentum.Eta());
            MuonPhi_n->push_back(newMuon.momentum.Phi());
            MuonIsIso_n->push_back(newMuon.isIso);
            MuonIsSignal_n->push_back(newMuon.isSignal);
            MuonIsBad_n->push_back(newMuon.isBad);
        }
    }


    //// electron vector
    //std::vector<myElectron> Electrons;
    float el_etcut = 10000;
    float el_d0sigcut = 5.;
    float el_z0cut = 0.5;
    float el_etacut = 2.47;
    for (const auto& ele : *electrons_nominal) {
        if (objTool->IsSignalElectron(*ele, el_etcut, el_d0sigcut, el_z0cut, el_etacut)) {
            myElectron newElectron;
            newElectron.momentum = ele->p4();
            ElePt_n->push_back(newElectron.momentum.Pt()/1000.);
            EleEta_n->push_back(newElectron.momentum.Eta());
            ElePhi_n->push_back(newElectron.momentum.Phi());
        }
    }

    //// photon vector
    //std::vector<myPhoton> Photons;
    float pho_ptcut = 20000.;
    float pho_etacut = 2.5;
    for (const auto& pho : *photons_nominal) {
        if (objTool->IsSignalPhoton(*pho, pho_ptcut, pho_etacut)) {
            myPhoton newPhoton;
            newPhoton.momentum = pho->p4();
            PhotonPt_n->push_back(newPhoton.momentum.Pt()/1000.);
            PhotonEta_n->push_back(newPhoton.momentum.Eta());
            PhotonPhi_n->push_back(newPhoton.momentum.Phi());
        }
    }


    // Calculate MET
    //EL_RETURN_CHECK("execute()", objTool->GetMET(*mettst_nominal,jets_nominal,electrons_nominal,muons_nominal,photons_nominal) );
    EL_RETURN_CHECK("execute()", objTool->GetMET(*mettst_nominal,jets_nominal,0,muons_nominal) );

    std::vector<myJet> GenJets;
    //std::vector<myJet> GenJetsNoNu;
    std::vector<myJet> GenJetsNoNuMu;

    // get truth container of interest
    const xAOD::JetContainer* genjets = 0;
    if (isMC) EL_RETURN_CHECK("execute()", event->retrieve( genjets, "AntiKt4TruthJets"));
    const xAOD::TruthParticleContainer* genparticles = 0;
    if (isMC) EL_RETURN_CHECK("execute()", event->retrieve( genparticles, "TruthParticles"));

    if (isMC) {

        // Fill truth jets to simplified jet container
        // add back neutrinos and muons to truth jets
        // CAREFUL: for now particles can be assigned to more than one truth jet
        for ( const auto& genjet : *genjets ) {
            myJet newGenJet;
            //myJet newGenJetNoNu;
            myJet newGenJetNoNuMu;
            TLorentzVector nuActivity(0., 0., 0., 0.);
            TLorentzVector muActivity(0., 0., 0., 0.);
            for ( const auto* it : *genparticles ) {
                if (it->pt() > 0.0 && it->status() == 1 && !(it->hasDecayVtx())) {
                    if ( abs(it->pdgId()) == 12 || abs(it->pdgId()) == 14 || abs(it->pdgId()) == 16 ) {
                        double dR = genjet->p4().DeltaR(it->p4());
                        if (dR < 0.4) {
                            nuActivity += it->p4();
                        }
                    }
                    if ( abs(it->pdgId()) == 13 ) {
                        double dR = genjet->p4().DeltaR(it->p4());
                        if (dR < 0.4) {
                            muActivity += it->p4();
                        }
                    }
                }
            }
            newGenJet.momentum = genjet->p4() + nuActivity  + muActivity;
            //newGenJetNoNu.momentum = genjet->p4() + muActivity;
            newGenJetNoNuMu.momentum = genjet->p4();
            newGenJet.btag = false;
            //newGenJetNoNu.btag = false;
            newGenJetNoNuMu.btag = false;
            if (abs(genjet->getAttribute<int>("PartonTruthLabelID")) == 5) {
                newGenJet.btag = true;
                //newGenJetNoNu.btag = true;
                newGenJetNoNuMu.btag = true;
            }
            GenJets.push_back(newGenJet);
            //GenJetsNoNu.push_back(newGenJetNoNu);
            GenJetsNoNuMu.push_back(newGenJetNoNuMu);

            GenJetPt_n->push_back(newGenJet.momentum.Pt()/1000.);
            GenJetEta_n->push_back(newGenJet.momentum.Eta());
            GenJetPhi_n->push_back(newGenJet.momentum.Phi());
            GenJetM_n->push_back(newGenJet.momentum.M()/1000.);
            GenJetBtag_n->push_back(newGenJet.btag);

            GenJetNoNuMuPt_n->push_back(newGenJetNoNuMu.momentum.Pt()/1000.);
            GenJetNoNuMuEta_n->push_back(newGenJetNoNuMu.momentum.Eta());
            GenJetNoNuMuPhi_n->push_back(newGenJetNoNuMu.momentum.Phi());
            GenJetNoNuMuM_n->push_back(newGenJetNoNuMu.momentum.M()/1000.);
            GenJetNoNuMuBtag_n->push_back(newGenJetNoNuMu.btag);

        }

    }

    GreaterByPt2<myJet> ptComparator_;
    std::sort(GenJets.begin(), GenJets.end(), ptComparator_);
    //std::sort(GenJetsNoNu.begin(), GenJetsNoNu.end(), ptComparator_);
    std::sort(GenJetsNoNuMu.begin(), GenJetsNoNuMu.end(), ptComparator_);

    // Fill reco jets to simplified jet container
    std::vector<myJet> JetsNoMu_rec;
    for ( const auto& jet : *jets_nominal) {
        myJet newJet;
        newJet.momentum = jet->p4();
        newJet.btag = false;
        if (objTool->IsBJet(*jet)) newJet.btag = true;
        newJet.jvt = jet->auxdata<float>("Jvt");
        newJet.good = true;
        if (objTool->IsBadJet(*jet)) newJet.good = false;
        JetsNoMu_rec.push_back(newJet);
    }

    //// Modify reco jets momentum by reconstruced muons in jet cone
    std::vector<myJet> Jets_rec = JetsNoMu_rec;
    if (AddMuToJets_) {
        for ( const auto& muon : *muons_nominal) {
            float mu_ptcut = 5000.;
            float mu_d0sigcut = 3.;
            float mu_z0cut = 0.5;
            float mu_etacut = 2.7;
            //// Add only high quality muons (without isolation requirement)
            if (!(objTool->IsNoIsoSignalMuon(*muon, mu_ptcut, mu_d0sigcut, mu_z0cut, mu_etacut)) || objTool->IsBadMuon(*muon, 0.2)) continue;
            bool muonAdded = false;
            for (std::vector<myJet>::iterator it = Jets_rec.begin(); (it != Jets_rec.end() && !muonAdded); ++it) {
                if (muon->p4().DeltaR(it->momentum) < 0.4) {
                    //std::cout << "vorher (pT, eta): " << it->momentum.Pt() << ", " << it->momentum.Eta() << std::endl;
                    it->momentum += muon->p4();
                    //std::cout << "nacher (pT, eta): " << it->momentum.Pt() << ", " << it->momentum.Eta() << std::endl;
                    muonAdded = true;
                }
            }
        }
    }

    LorentzVector MHTtotal(0,0,0,0);
    for (std::vector<myJet>::iterator it = Jets_rec.begin(); it != Jets_rec.end(); ++it) {
        JetPt_n->push_back(it->momentum.Pt()/1000.);
        JetEta_n->push_back(it->momentum.Eta());
        JetPhi_n->push_back(it->momentum.Phi());
        JetM_n->push_back(it->momentum.M()/1000.);
        JetBtag_n->push_back(it->btag);
        JetJVT_n->push_back(it->jvt);
        JetGood_n->push_back(it->good);
        if (it->good && ( it->momentum.Pt()>60000. || fabs(it->momentum.Eta()) > 2.4 || it->jvt > 0.59 ) ) MHTtotal -= it->momentum;
    }

    LorentzVector MHTNoMutotal(0,0,0,0);
    for (std::vector<myJet>::iterator it = JetsNoMu_rec.begin(); it != JetsNoMu_rec.end(); ++it) {
        JetNoMuPt_n->push_back(it->momentum.Pt()/1000.);
        JetNoMuEta_n->push_back(it->momentum.Eta());
        JetNoMuPhi_n->push_back(it->momentum.Phi());
        JetNoMuM_n->push_back(it->momentum.M()/1000.);
        JetNoMuBtag_n->push_back(it->btag);
        JetNoMuJVT_n->push_back(it->jvt);
        JetNoMuGood_n->push_back(it->good);
        if (it->good && ( it->momentum.Pt()>60000. || fabs(it->momentum.Eta()) > 2.4 || it->jvt > 0.59 ) ) MHTNoMutotal -= it->momentum;
    }

    LorentzVector METtotal(0,0,0,0);
    METtotal.SetPxPyPzE((*mettst_nominal)["Final"]->mpx(), (*mettst_nominal)["Final"]->mpy(), 0, (*mettst_nominal)["Final"]->met() );
    //std::cout << "METtotal (pt, phi): " << METtotal.Pt() << ", " << METtotal.Phi() << std::endl;
    MET_pt_ = METtotal.Pt() / 1000.;
    MET_phi_ = METtotal.Phi();

    LorentzVector METjet(0,0,0,0);
    METjet.SetPxPyPzE((*mettst_nominal)["RefJet"]->mpx(), (*mettst_nominal)["RefJet"]->mpy(), 0, (*mettst_nominal)["RefJet"]->met() );
    //std::cout << "METjet (pt, phi): " << METjet.Pt() << ", " << METjet.Phi() << std::endl;
    METjet_pt_ = METjet.Pt() / 1000.;
    METjet_phi_ = METjet.Phi();

    LorentzVector METmuon(0,0,0,0);
    METmuon.SetPxPyPzE((*mettst_nominal)["Muons"]->mpx(), (*mettst_nominal)["Muons"]->mpy(), 0, (*mettst_nominal)["Muons"]->met() );
    //std::cout << "METmuon (pt, phi): " << METmuon.Pt() << ", " << METmuon.Phi() << std::endl;
    METmu_pt_ = METmuon.Pt() / 1000.;
    METmu_phi_ = METmuon.Phi();

    LorentzVector METele(0,0,0,0);
    //METele.SetPxPyPzE((*mettst_nominal)["RefEle"]->mpx(), (*mettst_nominal)["RefEle"]->mpy(), 0, (*mettst_nominal)["RefEle"]->met() );
    //std::cout << "METele (pt, phi): " << METele.Pt() << ", " << METele.Phi() << std::endl;
    METele_pt_ = METele.Pt() / 1000.;
    METele_phi_ = METele.Phi();

    LorentzVector METgamma(0,0,0,0);
    //METgamma.SetPxPyPzE((*mettst_nominal)["RefGamma"]->mpx(), (*mettst_nominal)["RefGamma"]->mpy(), 0, (*mettst_nominal)["RefGamma"]->met() );
    //std::cout << "METgamma (pt, phi): " << METgamma.Pt() << ", " << METgamma.Phi() << std::endl;
    METgamma_pt_ = METgamma.Pt() / 1000.;
    METgamma_phi_ = METgamma.Phi();

    LorentzVector METtrack(0,0,0,0);
    METtrack.SetPxPyPzE((*mettst_nominal)[3]->mpx(), (*mettst_nominal)[3]->mpy(), 0, (*mettst_nominal)[3]->met() );
    //std::cout << "METtrack (pt, phi): " << METtrack.Pt() << ", " << METtrack.Phi() << std::endl;
    METtrack_pt_ = METtrack.Pt() / 1000.;
    METtrack_phi_ = METtrack.Phi();

    // gen MET from all status==1 gen particles
    // true MHT from all status==1 gen particles matched (dR < 0.4) to jets above rebalancing threshold
    LorentzVector vgenMET(0, 0, 0, 0);
    LorentzVector vtrueMHTreb(0, 0, 0, 0);

    if (calculateTrueMET_ && isMC) {
        for ( const auto* it : *genparticles ) {
            if (it->pt() > 0. && it->status() == 1 && !(it->hasDecayVtx())) {
                if (abs(it->pdgId()) == 12 || abs(it->pdgId()) == 14 || abs(it->pdgId()) == 16 )
                    vgenMET += it->p4();
                //std::cout << "Pid, status, pt, eta, phi: " << it->pdgId() << ", " << it->status() << ", " << it->pt() << ", " << it->eta() << ", " << it->phi() << std::endl;
                bool particleAdded = false;
                for (std::vector<myJet>::iterator jt = Jets_rec.begin(); (jt != Jets_rec.end() && !particleAdded) ; ++jt) {
                    if (jt->momentum.Pt() > rebalancedJetPt_ && jt->good) {
                        double dR = jt->momentum.DeltaR(it->p4());
                        if (dR < 0.4) {
                            vtrueMHTreb += it->p4();
                            particleAdded = true;
                        }
                    }
                }
            }
        }
    }
    GenMET_pt_ = vgenMET.Pt() / 1000.;
    GenMET_phi_ = vgenMET.Phi();
    TrueMHT_pt_ = vtrueMHTreb.Pt() / 1000.;
    TrueMHT_phi_ = vtrueMHTreb.Phi();


    //// what is going on with soft and muonic MET term
    if (( MET_pt_>100 || MHTtotal.Pt() > 100000) && debug_ >= 1 && isMC) {
        std::cout << "MHTtotal (pt, phi): " << MHTtotal.Pt() << ", " << MHTtotal.Phi() << std::endl;
        std::cout << "MHTNoMutotal (pt, phi): " << MHTNoMutotal.Pt() << ", " << MHTNoMutotal.Phi() << std::endl;
        std::cout << "METgamma (pt, phi): " << METgamma.Pt() << ", " << METgamma.Phi() << std::endl;
        std::cout << "METele (pt, phi): " << METele.Pt() << ", " << METele.Phi() << std::endl;
        std::cout << "METtrack (pt, phi): " << METtrack.Pt() << ", " << METtrack.Phi() << std::endl;
        std::cout << "METmuon (pt, phi): " << METmuon.Pt() << ", " << METmuon.Phi() << std::endl;
        std::cout << "METjet (pt, phi): " << METjet.Pt() << ", " << METjet.Phi() << std::endl;
        std::cout << "METtotal (pt, phi): " << METtotal.Pt()<< ", " << METtotal.Phi() << std::endl;
        int index = 0;
        for ( const auto& jet : *jets_nominal) {
            if (jet->pt() > 10000.) {
                std::cout << "Jet (" << index << ") (pt, eta, phi): " << jet->pt() << ", " << jet->eta() << ", " << jet->phi() << std::endl;
                ++index;
            }
        }
        index = 0;
        for (std::vector<myJet>::iterator it = Jets_rec.begin(); it != Jets_rec.end(); ++it) {
            if (it->momentum.Pt() > 10.) {
                std::cout << "MyRecoJet (" << index << ") (pt, eta, phi): " << it->momentum.Pt() << ", " << it->momentum.Eta() << ", " << it->momentum.Phi() << std::endl;
                ++index;
            }
        }
        index = 0;
        for ( const auto& genjet : *genjets ) {
            if (genjet->pt() > 10000.) {
                std::cout << "GenJet (" << index << ") (pt, eta, phi): " << genjet->pt() << ", " << genjet->eta() << ", " << genjet->phi() << std::endl;
                ++index;
            }
        }
        index = 0;
        for (std::vector<myJet>::iterator it = GenJets.begin(); it != GenJets.end(); ++it) {
            if (it->momentum.Pt() > 10.) {
                std::cout << "MyGenJet (" << index << ") (pt, eta, phi): " << it->momentum.Pt() << ", " << it->momentum.Eta() << ", " << it->momentum.Phi() << std::endl;
                ++index;
            }
        }
        index = 0;
        for ( const auto& photon : *photons_nominal) {
            if (photon->pt() > 5000.) {
                std::cout << "Photon (" << index << ") (pt, eta, phi): " << photon->pt() << ", " << photon->eta() << ", " << photon->phi() << std::endl;
                ++index;
            }
        }
        index = 0;
        for ( const auto& muons : *muons_nominal) {
            if (muons->pt() > 5000.) {
                std::cout << "Muon (" << index << ") (pt, eta, phi): " << muons->pt() << ", " << muons->eta() << ", " << muons->phi() << std::endl;
                ++index;
            }
        }
        index = 0;
        for ( const auto& electrons : *electrons_nominal) {
            if (electrons->pt() > 5000.) {
                std::cout << "Electrons (" << index << ") (pt, eta, phi): " << electrons->pt() << ", " << electrons->eta() << ", " << electrons->phi() << std::endl;
                ++index;
            }
        }
    }

    EventTree->Fill();

    GenJetPt_n->clear();
    GenJetEta_n->clear();
    GenJetPhi_n->clear();
    GenJetM_n->clear();
    GenJetBtag_n->clear();

    GenJetNoNuMuPt_n->clear();
    GenJetNoNuMuEta_n->clear();
    GenJetNoNuMuPhi_n->clear();
    GenJetNoNuMuM_n->clear();
    GenJetNoNuMuBtag_n->clear();

    JetPt_n->clear();
    JetEta_n->clear();
    JetPhi_n->clear();
    JetM_n->clear();
    JetBtag_n->clear();
    JetJVT_n->clear();
    JetGood_n->clear();

    JetNoMuPt_n->clear();
    JetNoMuEta_n->clear();
    JetNoMuPhi_n->clear();
    JetNoMuM_n->clear();
    JetNoMuBtag_n->clear();
    JetNoMuJVT_n->clear();
    JetNoMuGood_n->clear();

    ElePt_n->clear();
    EleEta_n->clear();
    ElePhi_n->clear();

    MuonPt_n->clear();
    MuonEta_n->clear();
    MuonPhi_n->clear();
    MuonIsBad_n->clear();
    MuonIsIso_n->clear();
    MuonIsSignal_n->clear();

    PhotonPt_n->clear();
    PhotonEta_n->clear();
    PhotonPhi_n->clear();

    // The containers created by the shallow copy are owned by you. Remember to delete them
    delete jets_nominal; // not these, we put them in the store!
    delete jets_nominal_aux;
    delete muons_nominal;
    delete muons_nominal_aux;
    delete electrons_nominal;
    delete electrons_nominal_aux;
    delete photons_nominal;
    delete photons_nominal_aux;
    delete mettst_nominal;
    delete mettst_nominal_aux;

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode NtupleMaker :: postExecute ()
{
    // Here you do everything that needs to be done after the main event
    // processing.  This is typically very rare, particularly in user
    // code.  It is mainly used in implementing the NTupleSvc.
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode NtupleMaker :: finalize ()
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

    // close our output file:
    //TFile *file = wk()->getOutputFile(outputfile_);
    //file->Close();

    Info("finalize()", "Number of clean events = %i", m_numCleanEvents);

    // GRL
    if (m_grl) {
        delete m_grl;
        m_grl = 0;
    }

    // XS data base
    delete my_XsecDB;

    // SUSYTools
    delete objTool;

    // cleaning up trigger tools
    if( m_trigConfigTool ) {
        delete m_trigConfigTool;
        m_trigConfigTool = 0;
    }
    if( m_trigDecisionTool ) {
        delete m_trigDecisionTool;
        m_trigDecisionTool = 0;
    }

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode NtupleMaker :: histFinalize ()
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

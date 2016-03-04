#include <MyAnalysis/MyResolution.h>

// Infrastructure include(s):
#include <AsgTools/MessageCheck.h>

#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/Message.h"

// EDM includes:
#include "xAODEventInfo/EventInfo.h"
#include "xAODJet/JetContainer.h"

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>

#include <TSystem.h>

/// Helper macro for checking xAOD::TReturnCode return values
#define EL_RETURN_CHECK( CONTEXT, EXP )                     \
	   do {                                                     \
	      if( ! EXP.isSuccess() ) {                             \
	         Error( CONTEXT,                                    \
	                XAOD_MESSAGE( "Failed to execute: %s" ),    \
	                #EXP );                                     \
	         return EL::StatusCode::FAILURE;                    \
	      }                                                     \
	   } while( false )


// this is needed to distribute the algorithm to the workers
ClassImp(MyResolution)

MyResolution :: MyResolution ()
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
    m_VetoCone=0.8;
    m_MatchingCone=0.1;
    m_RelGenActivityVeto=0.01;
    m_RelRecoActivityVeto=0.05;

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
    EtaBinEdges.push_back(0.5);
    EtaBinEdges.push_back(1.0);
    EtaBinEdges.push_back(1.5);
    EtaBinEdges.push_back(2.0);
    EtaBinEdges.push_back(2.5);
    EtaBinEdges.push_back(3.5);
    EtaBinEdges.push_back(5.0);

    //// Array of histograms for jet resolutions (all jet multiplicities)
    ResizeHistoVector(PtResponse);
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

            //// Book histograms Pt resolution
            TH1F* h_jetRes = new TH1F(GetHistName(i_pt, i_eta).c_str(), "p_T^{reco}/p_T^{gen}", 250, 0.0, 5.0);
            wk()->addOutput(h_jetRes);
            PtResponse.at(i_pt).at(i_eta) = h_jetRes;

        }
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

    // tell ANA_CHECK to return an EL::StatusCode
    ANA_CHECK_SET_TYPE (EL::StatusCode);

    // make and initialize the QuickAna tool
    std::unique_ptr<ana::QuickAna> myQuickAna (new ana::QuickAna ("quickana"));
    myQuickAna->setConfig (*this);
    quickAna = std::move (myQuickAna);
    ANA_CHECK (quickAna->initialize());

    // GRL
    m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    const char* GRLFilePath = "~/sonas/GRL/data15_13TeV.periodAllYear_DetStatus-v73-pro19-08_DQDefects-00-01-02_PHYS_StandardGRL_All_Good_25ns.xml";
    const char* fullGRLFilePath = gSystem->ExpandPathName (GRLFilePath);
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
    ANA_CHECK_SET_TYPE (EL::StatusCode);
    ANA_CHECK (quickAna->process ());

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
    int datasetID = 0;
    double eventWeight = 1;
    if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ) {
        isMC = true; // can do something with this later
        //  extra event-level information you might need:
        datasetID =  eventInfo->mcChannelNumber();
        eventWeight = eventInfo->mcEventWeight() ;
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

    //for (auto jet_itr : *quickAna->jets()) {
    //    if (jet_itr->auxdata<ana::SelectType> ("ana_select"))
    //        std::cout << jet_itr->pt() << std::endl;
    //}

    // get jet container of interest
    const xAOD::JetContainer* jets = 0;
    EL_RETURN_CHECK("execute()", event->retrieve( jets, "AntiKt4EMTopoJets" ));
    const xAOD::JetContainer* genjets = 0;
    EL_RETURN_CHECK("execute()", event->retrieve( genjets, "AntiKt4TruthJets"));

    // loop over the emjets and reject events with at least one bad jet
    xAOD::JetContainer::const_iterator jet_itr = jets->begin();
    xAOD::JetContainer::const_iterator jet_end = jets->end();
    for( ; jet_itr != jet_end; ++jet_itr ) {
        if( !m_jetCleaning->accept( **jet_itr ) && (*jet_itr)->pt() > 20.) {
            Info("execute()", "Reject event because of a bad jet");
            return EL::StatusCode::SUCCESS;
        }
    }

    // Print out number of jets
    //Info("execute()", "  number of jets = %lu", jets->size());

    // Loop over all jets in this container
    for ( const auto* genjet : *genjets ) {

        // check for no additional genJet activity
        bool noGenActivity = true;
        for ( const auto* genjet2 : *genjets ) {
            if (genjet2 == genjet) continue;
            double dR = genjet->p4().DeltaR(genjet2->p4());
            if (dR < m_VetoCone && genjet2->pt()/genjet->pt() > m_RelGenActivityVeto) {
                noGenActivity = false;
            }
        }
        //if (!noGenActivity) continue; // continue with next genJet if another genJet is closeby

        // check for additional recoJet activity
        const xAOD::Jet* matchedJet = new xAOD::Jet; //create a new jet object
        const xAOD::Jet* nextJet = new xAOD::Jet; //create a new jet object
        TLorentzVector addRecoActivity(0., 0., 0., 0.);
        double dRmin_matched = 999.;
        double dRmin_next = 999.;
        for ( auto jet : *quickAna->jets()) {
        //for ( const auto* jet : *jets ) {
            double dR = jet->p4().DeltaR(genjet->p4());
            if (dR < dRmin_matched && dR < m_VetoCone) {
                if (dRmin_matched < m_VetoCone ) {
                    dRmin_next = dRmin_matched;
                    nextJet = matchedJet;
                    addRecoActivity += nextJet->p4();
                }
                dRmin_matched = dR;
                matchedJet = jet;
            } else if (dR < dRmin_next && dR < m_VetoCone) {
                dRmin_next = dR;
                nextJet = jet;
                addRecoActivity += nextJet->p4();
            }
        } // end for loop over jets

        bool noRecoActivity = true;
        if (dRmin_matched < m_MatchingCone) {
            if (dRmin_next < m_VetoCone) {
                if (nextJet->pt()/matchedJet->pt()>m_RelRecoActivityVeto) noRecoActivity = false;
            }
        }

        if (dRmin_matched < m_MatchingCone && noRecoActivity && noGenActivity) {
            int i_pt = GetPtBin(genjet->pt()/1000.);
            int i_eta = GetEtaBin(genjet->eta());
            PtResponse.at(i_pt).at(i_eta)->Fill((matchedJet->p4()+addRecoActivity).Pt()/genjet->pt(), eventWeight);
        }

        //delete matchedJet; //why does this not work???

    } // end for loop over genjets

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



std::string MyResolution::GetHistName(unsigned int i_pt, unsigned int i_eta) {
    std::string hname = "h_tot_JetAll_ResponsePt_Pt";
    hname += std::to_string(i_pt);
    hname += "_Eta";
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


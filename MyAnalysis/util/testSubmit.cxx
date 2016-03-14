#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ScanDir.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoopGrid/PrunDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <TSystem.h>


#include <MyAnalysis/MyResolution.h>

int main( int argc, char* argv[] ) {

    // Take the submit directory from the input if provided:
    std::string submitDir = "submitDir";
    if( argc > 1 ) submitDir = argv[ 1 ];

    // Set up the job for xAOD access:
    xAOD::Init().ignore();

    // Construct the samples to run on:
    SH::SampleHandler sh;

    // use SampleHandler to scan all of the subdirectories of a directory for particular MC single file:
    //const char* inputFilePath = gSystem->ExpandPathName ("~/sonas/xAOD/mc15_13TeV.426136.Sherpa_CT10_jets_JZ6.merge.DAOD_SUSY1.e4355_s2608_r6869_r6282_p2470_tid07480677_00/");
    //SH::ScanDir().filePattern("DAOD_SUSY1.07480677._000001.pool.root.1").scan(sh, inputFilePath);

	SH::scanDQ2 (sh, "mc15_13TeV.426136.Sherpa_CT10_jets_JZ6.merge.DAOD_SUSY1.e4355_s2608_r6869_r6282_p2470_tid07480677_00");

    // Set the name of the input TTree. It's always "CollectionTree"
    // for xAOD files.
    sh.setMetaString( "nc_tree", "CollectionTree" );

    // Print what we found:
    sh.print();

    // Create an EventLoop job:
    EL::Job job;
    job.sampleHandler( sh );
    job.options()->setDouble (EL::Job::optMaxEvents, -1);

    // Add our analysis to the job:
    MyResolution* alg = new MyResolution();
    //alg->electronDef = "default";
    //alg->muonDef     = "default";
    //alg->tauDef      = "default";
    //alg->photonDef   = "default";
    alg->jetDef      = "default";
    alg->jetKine      = "pt>0";
    //alg->metDef      = "default";
    //alg->orDef       = "default";

    job.algsAdd( alg );

    // Run the job using the local/direct driver:
    //EL::DirectDriver driver;
	EL::PrunDriver driver;
	driver.options()->setString("nc_outputSampleName", "user.csander.test8.%in:name[2]%.%in:name[6]%");
    driver.submit( job, submitDir );

    return 0;
}

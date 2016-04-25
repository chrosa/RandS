void testRun( const std::string& submitDir ) {

    // Set up the job for xAOD access:
    xAOD::Init().ignore();

    // Construct the samples to run on:
    SH::SampleHandler sh;

    // use SampleHandler to scan all of the subdirectories of a directory for particular MC single file:
    //const char* inputFilePath = gSystem->ExpandPathName ("~/sonas/xAOD/");
    //SH::ScanDir().filePattern("DAOD_TOPQ1.07458545._000054.pool.root.1").scan(sh, inputFilePath);
    const char* inputFilePath = gSystem->ExpandPathName ("~/sonas/xAOD/mc15_13TeV.426134.Sherpa_CT10_jets_JZ4.merge.DAOD_SUSY1.e4355_s2608_r6869_r6282_p2470_tid07480766_00/");
    SH::ScanDir().filePattern("DAOD_SUSY1.07480766._000001.pool.root.1").scan(sh, inputFilePath);

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
    alg->jetKine     = "pt>0";
    //alg->metDef      = "default";
    //alg->orDef       = "default";

    job.algsAdd( alg );

    // Run the job using the local/direct driver:
    EL::DirectDriver driver;
    driver.submit( job, submitDir );

    return 0;
}

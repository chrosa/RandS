void RandSSubmit( const std::string& submitDir ) {

    // Set up the job for xAOD access:
    xAOD::Init().ignore();

    // Construct the samples to run on:
    SH::SampleHandler sh;

    //SH::addGrid (sh, "mc15_13TeV.426131.Sherpa_CT10_jets_JZ1.merge.DAOD_SUSY1.e4355_s2608_r6869_r6282_p2470_tid07480924_00");
    //SH::addGrid (sh, "mc15_13TeV.426132.Sherpa_CT10_jets_JZ2.merge.DAOD_SUSY1.e4355_s2608_r6869_r6282_p2470_tid07480940_00");
    SH::addGrid (sh, "mc15_13TeV.426133.Sherpa_CT10_jets_JZ3.merge.DAOD_SUSY1.e4355_s2608_r6869_r6282_p2470_tid07480935_00");
    SH::addGrid (sh, "mc15_13TeV.426134.Sherpa_CT10_jets_JZ4.merge.DAOD_SUSY1.e4355_s2608_r6869_r6282_p2470_tid07480766_00");
    SH::addGrid (sh, "mc15_13TeV.426135.Sherpa_CT10_jets_JZ5.merge.DAOD_SUSY1.e4355_s2608_r6869_r6282_p2470_tid07480768_00");
    SH::addGrid (sh, "mc15_13TeV.426136.Sherpa_CT10_jets_JZ6.merge.DAOD_SUSY1.e4355_s2608_r6869_r6282_p2470_tid07480677_00");

    //if you want to run over all datasets of one class using wildcards
    //SH::scanDQ2 (sh, "mc15_13TeV.*.Sherpa_CT10_jets_JZ*.merge.DAOD_SUSY1.e4355_s2608_r6869_r6282_p2470_*");

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
    RandS* alg = new RandS();
    //alg->electronDef = "default";
    //alg->muonDef     = "default";
    //alg->tauDef      = "default";
    //alg->photonDef   = "default";
    alg->jetDef      = "default";
    alg->jetKine     = "pt>0";
    //alg->metDef      = "default";
    //alg->orDef       = "default";

    // define an output and an ntuple associated to that output
	EL::OutputStream output("RandS");
	job.outputAdd (output);
	EL::NTupleSvc *ntuple = new EL::NTupleSvc("RandS");
	job.algsAdd(ntuple);

    job.algsAdd( alg );
    alg->outputfile_ = "RandS"; // give the name of the output to our algorithm
    
    // Run the job using the local/direct driver:
    //EL::DirectDriver driver;
    EL::PrunDriver driver;
    driver.options()->setString("nc_outputSampleName", "user.csander.RandS_RS_v1.%in:name[2]%.%in:name[6]%");
    driver.submit( job, submitDir );

}

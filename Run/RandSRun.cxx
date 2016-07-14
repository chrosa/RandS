void RandSRun( const std::string& submitDir ) {

    // Set up the job for xAOD access:
    xAOD::Init().ignore();

    // Construct the samples to run on:
    SH::SampleHandler sh;

	//SH::readFileList (sh, "sample", "filelist_input_pythia_test.txt");
	SH::readFileList (sh, "sample", "filelist_input_AOD.txt");

    // Set the name of the input TTree. It's always "CollectionTree"
    // for xAOD files.
    sh.setMetaString( "nc_tree", "CollectionTree" );

    // Print what we found:
    sh.print();

    // Create an EventLoop job:
    EL::Job job;
    job.sampleHandler( sh );
    job.options()->setDouble (EL::Job::optMaxEvents, -1);
    job.options()->setString (EL::Job::optXaodAccessMode, EL::Job::optXaodAccessMode_athena);

    // Add our analysis to the job:
    RandS* alg = new RandS();

    // define an output and an ntuple associated to that output
    EL::OutputStream output("RandS");
    job.outputAdd (output);
    EL::NTupleSvc *ntuple = new EL::NTupleSvc("RandS");
    job.algsAdd(ntuple);

    job.algsAdd( alg );
    alg->outputfile_ = "RandS"; // give the name of the output to our algorithm

    // Run the job using the local/direct driver:
    EL::DirectDriver driver;
    driver.submit( job, submitDir );

    return 0;
}

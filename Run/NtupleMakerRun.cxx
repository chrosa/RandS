void NtupleMakerRun( const std::string& submitDir ) {

    // Set up the job for xAOD access:
    xAOD::Init().ignore();

    // Construct the samples to run on:
    SH::SampleHandler sh;

	//SH::readFileList (sh, "sample", "filelist_input_AOD.txt");
	//SH::readFileList (sh, "sample", "filelist_input_Wtau.txt");
	//SH::readFileList (sh, "sample", "filelist_input_data2015.txt");
	SH::readFileList (sh, "sample", "filelist_input_data2016.txt");

    // Set the name of the input TTree. It's always "CollectionTree"
    // for xAOD files.
    sh.setMetaString( "nc_tree", "CollectionTree" );

    // Print what we found:
    sh.print();

    // Create an EventLoop job:
    EL::Job job;
    job.sampleHandler( sh );
    job.options()->setDouble (EL::Job::optMaxEvents, 10000);
    job.options()->setString (EL::Job::optXaodAccessMode, EL::Job::optXaodAccessMode_athena);

    // Add our analysis to the job:
    NtupleMaker* alg = new NtupleMaker();

    // define an output and an ntuple associated to that output
    EL::OutputStream output("NtupleMaker");
    job.outputAdd (output);
    EL::NTupleSvc *ntuple = new EL::NTupleSvc("NtupleMaker");
    job.algsAdd(ntuple);

    job.algsAdd( alg );
    alg->outputfile_ = "NtupleMaker"; // give the name of the output to our algorithm

    // Run the job using the local/direct driver:
    EL::DirectDriver driver;
    driver.submit( job, submitDir );

    return 0;
}

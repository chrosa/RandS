void NtupleMakerSubmit( const std::string& submitDir ) {

    // Set up the job for xAOD access:
    xAOD::Init().ignore();

    // Construct the samples to run on:
    SH::SampleHandler sh;

    SH::addGrid (sh, "data16_13TeV.periodA.physics_Main.PhysCont.DAOD_SUSY11.grp16_v01_p2667");
    SH::addGrid (sh, "data16_13TeV.periodB.physics_Main.PhysCont.DAOD_SUSY11.grp16_v01_p2667");
    SH::addGrid (sh, "data16_13TeV.periodC.physics_Main.PhysCont.DAOD_SUSY11.grp16_v01_p2689");
    SH::addGrid (sh, "data16_13TeV.periodD.physics_Main.PhysCont.DAOD_SUSY11.grp16_v01_p2689");
    SH::addGrid (sh, "data16_13TeV.periodE.physics_Main.PhysCont.DAOD_SUSY11.grp16_v01_p2689");
    SH::addGrid (sh, "data16_13TeV.periodF.physics_Main.PhysCont.DAOD_SUSY11.grp16_v01_p2689");
    SH::addGrid (sh, "data16_13TeV.periodG.physics_Main.PhysCont.DAOD_SUSY11.grp16_v01_p2769");
    SH::addGrid (sh, "data16_13TeV.periodI.physics_Main.PhysCont.DAOD_SUSY11.grp16_v01_p2769");


    // Set the name of the input TTree. It's always "CollectionTree"
    // for xAOD files.
    sh.setMetaString( "nc_tree", "CollectionTree" );

    // Print what we found:
    sh.print();

    // Create an EventLoop job:
    EL::Job job;
    job.sampleHandler( sh );
    job.options()->setString (EL::Job::optXaodAccessMode, EL::Job::optXaodAccessMode_athena);
    job.options()->setDouble (EL::Job::optMaxEvents, -1);
	job.options()->setDouble (EL::Job::optGridMergeOutput, 0);

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
	EL::PrunDriver driver;
    driver.options()->setString("nc_outputSampleName", "user.csander.NtupleMaker_v5.%in:name[2]%.%in:name[6]%");
    driver.submitOnly( job, submitDir );


}

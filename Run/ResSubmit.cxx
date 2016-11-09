void ResSubmit( const std::string& submitDir ) {

    // Set up the job for xAOD access:
    xAOD::Init().ignore();

    // Construct the samples to run on:
    SH::SampleHandler sh;

    SH::addGrid (sh, "mc15_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.merge.DAOD_SUSY11.e3668_s2576_s2132_r7725_r7676_p2666_tid08818957_00");
    SH::addGrid (sh, "mc15_13TeV.361023.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3W.merge.DAOD_SUSY11.e3668_s2576_s2132_r7725_r7676_p2666_tid08821580_00");
    SH::addGrid (sh, "mc15_13TeV.361023.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3W.merge.DAOD_SUSY11.e3668_s2576_s2132_r7725_r7676_p2666_tid08821586_00");
    SH::addGrid (sh, "mc15_13TeV.361024.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4W.merge.DAOD_SUSY11.e3668_s2576_s2132_r7725_r7676_p2666_tid08818756_00");
    SH::addGrid (sh, "mc15_13TeV.361024.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ4W.merge.DAOD_SUSY11.e3668_s2576_s2132_r7725_r7676_p2666_tid08818760_00");
    SH::addGrid (sh, "mc15_13TeV.361025.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5W.merge.DAOD_SUSY11.e3668_s2576_s2132_r7725_r7676_p2666_tid08818562_00");
    SH::addGrid (sh, "mc15_13TeV.361025.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ5W.merge.DAOD_SUSY11.e3668_s2576_s2132_r7725_r7676_p2666_tid08818565_00");
    SH::addGrid (sh, "mc15_13TeV.361026.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ6W.merge.DAOD_SUSY11.e3569_s2608_s2183_r7725_r7676_p2666_tid08821040_00");

    // Set the name of the input TTree. It's always "CollectionTree" for xAOD files.
    sh.setMetaString( "nc_tree", "CollectionTree");

    // Print what we found:
    sh.print();

    // Create an EventLoop job:
    EL::Job job;
    job.sampleHandler( sh );
    job.options()->setDouble (EL::Job::optMaxEvents, -1);

    // Add our analysis to the job:
    MyResolution* alg = new MyResolution();
    job.algsAdd( alg );

    // Run the job using the local/direct driver:
	EL::PrunDriver driver;
    driver.options()->setString("nc_outputSampleName", "user.csander.Resolution_FineEta_v1.%in:name[2]%.%in:name[6]%");
    driver.submit( job, submitDir );
}

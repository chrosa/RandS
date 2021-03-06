void NtupleMakerNafSubmit( const std::string& submitDir ) {

    // Set up the job for xAOD access:
    xAOD::Init().ignore();

    // Construct the samples to run on:
    SH::SampleHandler sh;

	SH::readFileList (sh, "sample", "filelist_input_Wtau.txt");
	//SH::readFileList (sh, "sample", "filelist_input_pythia.txt");
	//SH::readFileList (sh, "sample", "filelist_input_data2015.txt");
	//SH::readFileList (sh, "sample", "filelist_input_data2016.txt");
    
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
    job.options()->setDouble(EL::Job::optFilesPerWorker, 3); // 10 for data // 3 for MC
	job.options()->setString(EL::Job::optSubmitFlags, "-o /nfs/dust/atlas/user/csander/logs -e /nfs/dust/atlas/user/csander/logs -S /bin/bash -l h_rt=24:00:00 -l h_vmem=4000M -l distro=sld6");


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
    EL::SoGEDriver driver;
    driver.shellInit = "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase; export DQ2_LOCAL_SITE_ID=DESY-HH_SCRATCHDISK; source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh; cd /afs/desy.de/user/c/csander/xxl-af-cms/testarea/2.4.8; lsetup rcsetup; cd -";
    driver.options()->setString("nc_outputSampleName", "user.csander.NtupleMaker_v1.%in:name[2]%.%in:name[6]%");
    driver.submitOnly( job, submitDir );

}

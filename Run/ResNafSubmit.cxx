void ResNafSubmit( const std::string& submitDir ) {

    // Set up the job for xAOD access:
    xAOD::Init().ignore();

    // Construct the samples to run on:
    SH::SampleHandler sh;

	SH::readFileList (sh, "sample", "filelist_input_pythia.txt");

    // Set the name of the input TTree. It's always "CollectionTree" for xAOD files.
    sh.setMetaString( "nc_tree", "CollectionTree");

    // Print what we found:
    sh.print();

    // Create an EventLoop job:
    EL::Job job;
    job.sampleHandler( sh );
    job.options()->setDouble(EL::Job::optMaxEvents, -1);
    job.options()->setDouble(EL::Job::optFilesPerWorker, 5);
    job.options()->setString(EL::Job::optSubmitFlags, "-o /nfs/dust/atlas/user/csander/logs -e /nfs/dust/atlas/user/csander/logs -S /bin/bash -l h_cpu=01:00:00 -l h_vmem=3000M -l distro=sld6");

    // Add our analysis to the job:
    MyResolution* alg = new MyResolution();
    job.algsAdd( alg );

    // Run the job using the local/direct driver:
    EL::SoGEDriver driver;
    driver.shellInit = "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase; export DQ2_LOCAL_SITE_ID=DESY-HH_SCRATCHDISK; source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh; cd /afs/desy.de/user/c/csander/xxl-af-cms/testarea/2.4.8; lsetup rcsetup; cd -";
    driver.options()->setString("nc_outputSampleName", "user.csander.MCRes_pythia_v5.%in:name[2]%.%in:name[6]%");
    driver.submit( job, submitDir );

}

void nafSubmit( const std::string& submitDir ) {

    // Set up the job for xAOD access:
    xAOD::Init().ignore();

    // Construct the samples to run on:
    SH::SampleHandler sh;

	SH::readFileList (sh, "sample", "filelist_input_pythia.txt");

    //if you want to run over all datasets of one class using wildcards
    //SH::scanDQ2 (sh, "mc15_13TeV.*.Sherpa_CT10_jets_JZ*.merge.DAOD_SUSY1.e4355_s2608_r6869_r6282_p2470_*");

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
    EL::SoGEDriver driver;
    driver.shellInit = "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase; export DQ2_LOCAL_SITE_ID=DESY-HH_SCRATCHDISK; source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh; cd /afs/desy.de/user/c/csander/xxl-af-cms/testarea/2.3.41; lsetup rcsetup; cd -";
    driver.options()->setString("nc_outputSampleName", "user.csander.MCRes_pythia_v3.%in:name[2]%.%in:name[6]%");
    driver.submit( job, submitDir );

}

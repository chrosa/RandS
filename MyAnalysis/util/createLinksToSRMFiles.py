#!/usr/bin/env python2
import os, sys, subprocess, argparse

CHRISTIAN_DEFAULT_PATH = "/afs/desy.de/user/c/cluedtke/NEWLINKS"
KILIAN_DEFAULT_PATH = "/afs/desy.de/user/k/krosbach/NEWLINKS"
PHILIPP_DEFAULT_PATH = "/afs/desy.de/user/p/pmogg/storage/public/CxAOD/NEWLINKS"
MICHAEL_DEFAULT_PATH = "/afs/desy.de/user/m/mischube/NEWLINKS"
CSANDER_DEFAULT_PATH = "/afs/desy.de/user/c/csander/sonas/NEWLINKS"

parser = argparse.ArgumentParser(description='Look up the path to each file of a sample and create a symbolic link to it.\nIMPORTANT NOTE: You need a valid grid token and rucio set up for this script, this is not checked by the script.')
parser.add_argument('samples', type=str,
                    help='Pattern which will be used by rucio to get a list of datasamples')
parser.add_argument('--user', '-u', type=str,
                    help='Grid username. In case your local username is different from the grid username, you need to specify this.')
parser.add_argument('--output_dir', '-o',
    help='Path to the output directory where the symbolic links will be stored. Will be created in case it doesn\'t exist yet. (default: hard-coded, user-specific)')
args = parser.parse_args()
pattern = args.samples
outputdir = args.output_dir
username = args.user
# Check if rucio was set
if not os.getenv("RUCIO_ACCOUNT"):
  print "Set up rucio before using this script by typing \"lsetup rucio\". Exiting."
  exit()
  #subprocess.call("localSetupRucioClients", shell=True)
import rucio.client

ruc = rucio.client.Client()

if not username: username = ruc.whoami()["account"] 
# Assume, that we use this script only on samples, which were made by the user with name username
#scopename = "user." + username
#scopename = "data16_13TeV"
scopename = "mc15_13TeV"
print "Search for samples with pattern \"" + pattern + "\" in scope \"" + scopename + "\""

#####################
# Start of the script
#####################
samples = []
for did in ruc.list_dids(scope=scopename, filters={'name':pattern}, type='container', long=False):
  samples.append(did)
if not samples:
  sys.exit("There were no samples found for the pattern \"" + pattern + "\" in scope \"" + scopename + "\". (Maybe change the scopename by specifying the user?) Exiting.")

# Make a list of sample using rucio
#searchagain = True
#samples = []
#tryDiffPattern = False
#while searchagain:
#    if not pattern or tryDiffPattern:
#        pattern = raw_input("Please enter the name of the dataset: ")
#    samplename = pattern
#    for did in ruc.list_dids(scope=scopename, filters={'name':samplename}, type='container', long=False):
#        samples.append(did)
#    if not samples:
#        print samples
#        user_dec = ""
#        while user_dec!="y" and user_dec!="n":
#            user_dec = raw_input("There were no samples found for the pattern " + pattern + ", do you want to try a different one? [y/n] ")
#        if user_dec=="n": exit()
#        searchagain = True
#        tryDiffPattern = True
#    else: searchagain = False

# For each sample, make a list of corresponding files
# For each file, check whether the file can be found on the DESY-HH_LOCALGROUPDISK
# If a file is not there, print a message and exit
print "Checking if all datasets were properly replicated to DESY-HH_LOCALGROUPDISK. This might take a few minutes."
samples_files = []
bad_samples = []
i_sample = 0
n_samples = str(len(samples))
for sample in samples: 
    i_sample += 1
    print "Checking sample " + str(i_sample) + " of " + n_samples + ": " + sample
    files = []
    n_bad_files = 0
    n_tot_files = 0
    all_files_ok = True
    for f in ruc.list_replicas([{'scope':scopename, 'name':sample}]):
        rse_path = f['rses'].get('DESY-HH_LOCALGROUPDISK')
        n_tot_files += 1
        if not rse_path:
            print "WARNING: Some files of dataset " + sample + " cannot be found on DESY-HH_LOCALGROUPDISK. Maybe the sample is not finished on grid or replication is stuck/failed."
            all_files_ok = False
            n_bad_files += 1
            continue
        srm_path_unformatted = rse_path[0]
        if not 'srm://dcache-se-atlas.desy.de' in srm_path_unformatted or not '/pnfs/' in srm_path_unformatted:
            print "WARNING: rucio output not formatted as expected. Cannot properly extract the logical file location from string:\n" + srm_path_unformatted + ". Skipping file."
            all_files_ok = False
            n_bad_files += 1
            continue
        srm_path = srm_path_unformatted[srm_path_unformatted.find('/pnfs/'):]
        files.append(srm_path)
    if not all_files_ok: bad_samples.append([sample, n_tot_files, n_bad_files])
    samples_files.append([sample, files])

if not bad_samples: print "All samples are OK! Proceeding to create links."
else:
    print "These sample cannot be found at all or only incomplete at DESY-HH_LOCALGROUPDISK:"
    for bs in bad_samples: 
        print "\t" + bs[0] + ": missing " + str(bs[2]) + " out of " + str(bs[1]) + " files."
    user_dec = ""
    while user_dec!="y" and user_dec!="n":
        user_dec = raw_input("Do you want to proceed to create links? [y/n] ")
    if user_dec=="n": exit()

# Check if the directory exists. If not, create it.
if not outputdir:
    if username=="cluedtke":
        outputdir = CHRISTIAN_DEFAULT_PATH
    elif username=="krosbach":
        outputdir = KILIAN_DEFAULT_PATH
    elif username=="pmogg":
        outputdir = PHILIPP_DEFAULT_PATH
    elif username=="mischube":
        outputdir = MICHAEL_DEFAULT_PATH
    elif username=="csander":
        outputdir = CSANDER_DEFAULT_PATH
    else:
        outputdir = raw_input("Please enter the path to the folder which shall hold the links: ")

if not outputdir.endswith('/'): outputdir+='/'

if os.path.isdir(outputdir):
    print "Write links into " + outputdir
else:
    os.makedirs(outputdir)
    print "Created " + outputdir + ". Write symbolic links here."

for sf in samples_files:
    sample = sf[0]
    dirname = sample.split(".root")[0] # drop .root extension if present
    if os.path.isdir(os.path.join(outputdir, dirname)):
        print "Output directory", outputdir+dirname, "exists, skipping."
        continue

    os.makedirs(os.path.join(outputdir, dirname))
    for logical_file_path in sf[1]:
        linkname = logical_file_path.split("/")[-1]
        os.symlink(logical_file_path, os.path.join(outputdir, dirname, linkname))

print "Done creating links. Check if everything went as expected e.g.: ls " + outputdir
print "If you're pleased with the result, proceed with the organizeSamples.py script."

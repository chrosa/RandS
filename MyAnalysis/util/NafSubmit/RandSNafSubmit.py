#!/usr/bin/python
from os import path, listdir, system
from ROOT import TTree, TFile, gDirectory

def submit_job(index, jobDir):
	system("sed 's/jobid/"+str(index)+"/' run.template > tmp.sh")
	system("sed 's/jobdir/"+str(jobDir)+"/' tmp.sh > run.sh")
	system("mkdir /nfs/dust/atlas/user/csander/RandS/Output/"+jobDir+"/"+str(index))
	system("cp filelist_mc.txt /nfs/dust/atlas/user/csander/RandS/Output/"+jobDir+"/"+str(index))
	system("cp run.sh /nfs/dust/atlas/user/csander/RandS/Output/"+jobDir+"/"+str(index))
	system("cd /nfs/dust/atlas/user/csander/RandS/Output/"+jobDir+"/"+str(index))
	system("qsub -cwd -o /nfs/dust/atlas/user/csander/RandS/Output/"+jobDir+"/"+str(index)+" -e /nfs/dust/atlas/user/csander/RandS/Output/"+jobDir+"/"+str(index)+" -S /bin/bash -l h_rt=24:00:00 -l h_vmem=4000M -l distro=sld6 run.sh");
	system("cd /afs/desy.de/user/c/csander/xxl-af-cms/testarea/2.4.8/MyAnalysis/util/NafSubmit")

#filelist = open("../filelist_mc_Wtau_v1.txt")
#filelist = open("../filelist_mc_v1.txt")
filelist = open("../filelist_data_v2.txt")
line = "init"
jobDir = "MyTest_data_METsoftPhiSmear_noAngSmear_smear15_N20_v3"
#500k for sim
#2.5M for data
#2.5M for Wtau
Nevts_max = 500000

system("mkdir /nfs/dust/atlas/user/csander/RandS/Output/"+jobDir)
# loop over ntuple files
i = 0
Nevts_job = 0
submit = True
while line:
	line = filelist.readline()
	print line
	if len(line) != 0:
		if submit:
			fout = open("filelist_mc.txt","w")
		fout.write(line)
		# open the input file
		fin = TFile(line.strip('\n\r'))
		tree = gDirectory.Get("EventTree")
		Nevts_file = tree.GetEntries()
		print "Entries in file: ", Nevts_file
		fin.Close()
		Nevts_job += Nevts_file
		if Nevts_job > Nevts_max and Nevts_file > 0:
			print "Time to submit, NEvts: ", Nevts_job
			i = i+1
			fout.close()
			submit_job(i, jobDir)
			Nevts_job = 0
			submit = True
		else:
			print "Not ready to submit, NEvts: ", Nevts_job
			submit = False
	else:
		if not submit:
			i = i+1
			fout.close()
			submit_job(i, jobDir)

# end of loop over events

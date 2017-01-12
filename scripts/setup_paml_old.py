## Author:       Zane Goodwin
## Contact:      gzane01@gmail.com
## Date Written: 1/5/2015
## Last Edited:  11/9/2015

import sys
from os import path
from subprocess import Popen, PIPE


## Runs a process through the bash terminal
## Parameters:
#		process_string: A string containing the commmand text
## Return value:
#		None
def runProcess(process_string):

	cmd = process_string.split(" ")
	process = Popen(cmd, stdout=PIPE, stderr=PIPE)
	stdout, stderr = process.communicate()
	return

args=sys.argv

ctlFile = args[1]
seqFile = args[2]
treeFile = args[3]
outFile = args[4]
outDir = args[5]
hyp = args[6]

# Find the gene name in the working directory
geneName = path.basename(seqFile).split("_")[0]
newCtlFile = geneName + "_" + hyp + ".ctl"

outDir = outDir + "/" + geneName
outCtl = outDir + "/" + newCtlFile

# Create directory to store the PAML runs, move the tree and seq files into 
#  that directory
runProcess("mkdir " + outDir)
runProcess("cp " + treeFile + " " + outDir)
runProcess("cp " + seqFile + " " + outDir)

# Rename the sequence file as its base name
seqFile = path.basename(seqFile)


# Write the new control file
o = open(outCtl, 'w')
with open(ctlFile, 'r') as f:
	for line in f:

		l=line.strip()
		fields = l.split(' ')
		
		# varname is the name of a variable in the control file
		varname = fields[0]

		if varname == "seqfile":
			fields[2] = seqFile

		elif varname == "treefile":
			fields[2] = treeFile

		elif varname == "outfile":
			fields[2] = outFile

		newLine = " ".join(fields)
		
		o.write(newLine + "\n")

o.close()
"""

Date Written: 1/5/2015
Last Edited:  1/12/2017

@author: Zane Goodwin
"""


def setup_parser():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description='''For setting up folders to run PAML jobs.''')
	### Required positional arguments
	parser.add_argument('ctlFile',help="Control file", type=str)
	parser.add_argument('seqFileDir',help="Directory of alignment files in Phylip format", type=str)
	parser.add_argument('treeFile',help="Unweighted Tree file in newick format", type=str)
	parser.add_argument('model', help="Name of folder in which to store model results", type=str)
    
	return parser


def edit_ctl_file(ctl_file_name, new_ctl_file_name, seqFile, treeFile, outFile):
	"""Replaces the default entries in the control file with those on the command line.
	Args:
		ctl_file_name: Path of the PAML control file (str)
		new_ctl_file-name: Desired name of the edited PAML control file (str)
		seqFile: Path of the current sequence file (str)
		treeFile: Path of the tree file (str)
		outFile: Path of the output file

	Returns: Nothing.
	"""
	o = open(new_ctl_file_name, 'w')
	with open(ctl_file_name, 'r') as f:
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
	return


def runProcess(process_string):
	cmd = process_string.split(" ")
	process = Popen(cmd, stdout=PIPE, stderr=PIPE)
	stdout, stderr = process.communicate()
	return


def main():
	global parser
	global gene_list

	# Get parser arguments
	parser = setup_parser()
	args   = parser.parse_args()

	model      = args.model
	treeFile   = args.treeFile
	seqFileDir = args.seqFileDir

	# Set up model directory and get the list of genes
	runProcess("mkdir " + model)
	gene_list = listdir(seqFileDir)

	for g in gene_list:
		
		# make a new folder for each gene to store the PAML run
		geneName = path.basename(g).split("_")[0]
		geneDir = "./%s/%s" % (model, geneName)
		runProcess("mkdir " + geneDir)

		# Generate new paml control file 
		newCtlFile = "%s/%s_%s.ctl" % (geneDir, geneName, model)
		outFile = "%s-%s.mlc" % (geneName, model)
		edit_ctl_file(args.ctlFile, newCtlFile,
						g, path.basename(treeFile), outFile)

		# Move tree and alignment file into the seq dir
		runProcess("cp " + treeFile + " " + geneDir)
		runProcess("cp " + seqFileDir + "/" + g + " " + geneDir)


if __name__ == "__main__":
	import sys, argparse
	from os import path, listdir
	from subprocess import Popen
	from subprocess import PIPE
	main()

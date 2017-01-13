"""

Date Written: 1/12/2017
Last Edited:  1/12/2017

@author: Zane Goodwin
"""


def setup_parser():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description='''Runs PAML jobs.''')
	### Required positional arguments
	parser.add_argument('model', help="Name of model to run in PAML", type=str)
	parser.add_argument('paml_executable',help="Location of the desired PAML script (usually codeml)", type=str)
	return parser


def runProcess(process_string):
	cmd = process_string.split(" ")
	process = Popen(cmd, stdout=PIPE, stderr=PIPE)
	stdout, stderr = process.communicate()
	return [stdout, stderr]


def writeOutput(stdout, stderr, gene):
	with open(gene + ".stdout", 'w') as o:
		for line in stdout:
			o.write(str(line))
	o.close()
	
	with open(gene + ".stderr", 'w') as e:
		for line in stderr:
			e.write(str(line))
	e.close()
	return


def main():
	global paml_executable
	global gene_list
	global parser

	# Parse command line arguments
	parser = setup_parser()
	args = parser.parse_args()
	model = args.model
	paml_executable = args.paml_executable

	# PAML can only print all of its output files if you run it from the same
	# directory as the seqeunce file. Here, we change into the directory for a 
	# given gene, run paml, then change to the directory of the next gene.

	chdir(model)
	for g in glob.glob('*'):
		chdir(g)
		ctlFile = glob.glob('./*.ctl')[0]
		print("Now running PAML on " + g + "...")
		[out, err] = runProcess(paml_executable + " " + ctlFile)
		writeOutput(out,err,g)
		chdir('../')
	chdir('../')


if __name__ == "__main__":
	import sys, argparse, glob
	from os import path, listdir, chdir
	from subprocess import Popen
	from subprocess import PIPE
	main()

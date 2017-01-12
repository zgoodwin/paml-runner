"""

Date Written: 1/12/2017
Last Edited:  1/12/2017

@author: Zane Goodwin
"""


def runProcess(process_string):
	cmd = process_string.split(" ")
	process = Popen(cmd, stdout=PIPE, stderr=PIPE)
	stdout, stderr = process.communicate()
	return [stdout, stderr]


def writeOutput(stdout, stderr, gene):
	with open(gene + ".stdout", 'w') as o:
		for line in stdout:
			o.write(line)
	o.close()
	
	with open(gene + ".stderr", 'w') as e:
		for line in stderr:
			e.write(line)
	e.close()
	return


def main():
	global paml_executable
	global gene_list

	args = sys.argv
	model = args[1]
	paml_executable = args[2]

	chdir(model)
	for g in listdir('./'):
		chdir(g)
		ctlFile = glob.glob('./*.ctl')[0]
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

"""

Date Written: 1/12/2017
Last Edited:  1/12/2017

@author: Zane Goodwin
"""

def setup_parser():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description='''For parsing dN and dS matrices from PAML and for calculating chi2 p-values for each gene''')
	parser.add_argument("altHypDir", help="Path to the alternative hypothesis directory", type=str)
	parser.add_argument("nullHypDir", help="Path to the null hypothesis directory", type=str)
	parser.add_argument("outfile", help="Name of text file to store estimated rate values", type=str)
	# parser.add_argument("foregroundClade", help="Comma-separated list of species to include in the foreground clade rate calculations", type=str)
	# parser.add_argument("backgroundClade", help="Comma-separated list of species to include in the background clade rate calculations", type=str)
	return parser


def parseDistanceMatrix(m):
	"""Reads the distance matrix 
	Args:
		m: The name of the dN or dS file produced by PAML

	Returns:
		matrix: a dictionary representation of the dN or dS matrix 
	"""

	lineCount = 0
	matrix = {}
	with open(m, 'r') as f:
		for line in f:
			if lineCount != 0:
				l = line.strip().split()
				l.append('0.000')		# Add a zero for the diagonal
				species = l[0]			# Get the species name
				values = l[1:(len(l))]	# Get values attached to the species
				newValues = []
				for v in values:
					newValues.append(float(v))
				matrix[species] = newValues
			lineCount += 1
	return matrix


def replaceZeros(nums, value):
	newNums = []
	for i in nums:
		if i == 0.0:
			newNums.append(value)
		else:
			newNums.append(i)
	newNums = np.array(newNums)
	return newNums


def divideMatrix(matrixA, matrixB):
	"""Element-wise division of matrixA by matrixB.
	
	Args:
		matrixA: A list of numpy arrays representing the lower half of a  
			diagonal matrix
		matrixB: Same as matrixA
	Returns:
		quotientMatrix: A matrix containing the element-wise quotients of
			matrixA and matrixB
	"""

	# Check that matrixA and matrixB have the same dimensions
	assert len(matrixA.keys()) == len(matrixB.keys()), "Matrix sizes not equal!"

	quotientMatrix = {}
	pseudoDenominator = 999.0
	for k in matrixA.keys():
		# Convert the rows to numpy array objects
		rowMatrixA = np.array(matrixA[k], dtype=np.float32)
		rowMatrixB = np.array(matrixB[k], dtype=np.float32)
		# Replace zeros in matrixB to avoid divide-by-zero errors
		nonZeroRowMatrixB = replaceZeros(rowMatrixB, pseudoDenominator)
		quotientRow = rowMatrixA/nonZeroRowMatrixB
		quotientMatrix[k] = quotientRow
	return quotientMatrix


def calculateAverageRates(m):
	total = 0.0
	numElements = 0.0
	averageRate = 0.0
	for k in m.keys():
		total += sum(m[k])
		numElements += float(len(m[k]))
	averageRate = (total/numElements)
	return averageRate


def findLikelihood(mlcFile):
	"""Locates the likelihood value in the PAML mlc file
	"""

	lnLikelihood = 0.0
	with open(mlcFile,'r') as f:
		for line in f:
			if re.match('^lnL', line):
				lnLikelihood = line.split()[4]
	return float(lnLikelihood)


def calculateChi2(null, alt):
	chi2 = (-2.0) * (null - alt)
	return float(chi2)


def calculatePval(chi2_stat, df):
	pval = 2 * (stats.chi2.pdf(chi2_stat, df))
	return float(pval)

# This is supposed to return a pandas data frame with the LRT results.
# NOTE: You need to add error handling here because sometimes PAML does
#		not produce output files!! 
def parse_rates(altHypDir, nullHypDir):
	"""Parses the dN, dS, dN/dS, chi2 and p-values from PAML output.

	Reads an mlc output file from paml to get likelihood values from the null  
	and alternate models. Obtains dN and dS values from the "2NG.dN" and  
	"2NG.dS" output files from PAML and calculates the dN/dS ratios, plus their
	associated chi2 statistics and p-values.

	Args:
		altHypDir: Path to the alternative hypothesis folder
		nullHypDir: Path to the null hypothesis folder
	Returns:
		results_frame: a pandas data frame containing all of the statistics
	"""

	dnList = []
	dsList = []
	dndsList = []
	geneList = []
	altLikelihoods = []
	nullLikelihoods = []
	chi2_vals = []
	pvals = []

	altDirs = glob.glob(altHypDir + '*')
	nullDirs = glob.glob(nullHypDir + '*')

	# Loop through genes in the alternate hypothesis folder
	for gene in altDirs:
		dnFile = gene + "/2NG.dN"
		dsFile = gene + "/2NG.dS"
		resultFile = glob.glob(gene + "/*.mlc")[0]
		altLike = findLikelihood(resultFile)

		if path.isfile(dnFile):
			dnMatrix = parseDistanceMatrix(dnFile)
			dsMatrix = parseDistanceMatrix(dsFile)
			dndsMatrix = divideMatrix(dnMatrix, dsMatrix)
			
			try:
				avgdN = calculateAverageRates(dnMatrix)
				avgdS = calculateAverageRates(dsMatrix)
				avgdNdS = calculateAverageRates(dndsMatrix)
				# fg = foregroundClade.split(",")
				# bg = backgroundClade.split(",")
				# calculateAverageCladeRates(dndsMatrix, fg, bg)
			except ZeroDivisionError:
				print("The test failed on " + gene + "... check mlc file")
				avgdN = 0.0
				avgdS = 0.0
				avgdNdS = 0.0

			geneList.append(gene.split("/")[1])
			dnList.append(avgdN)
			dsList.append(avgdS)
			dndsList.append(avgdNdS)
			altLikelihoods.append(altLike)

	# Loop through genes in the null hypothesis folder
	for gene in nullDirs:
		resultFile = glob.glob(gene + "/*.mlc")[0]
		nullLike = findLikelihood(resultFile)
		nullLikelihoods.append(nullLike)

	# Calculate likelihoods and p-values
	for i in range(0,len(geneList)):
		n = float(nullLikelihoods[i])
		a = float(altLikelihoods[i])
		chi2 = calculateChi2(n,a)
		pval = calculatePval(chi2, 2)
		chi2_vals.append(chi2)
		pvals.append(pval)

	# Put the data into a pandas data frame
	results = {'genes' : geneList,
			   'dN' : dnList,
			   'dS' : dsList,
			   'dN_dS' : dndsList,
			   'alt_likelihood' : altLikelihoods,
			   'null_likelihood' : nullLikelihoods,
			   'chi2_stat' : chi2_vals,
			   'pval' : pvals}
	results_frame = pd.DataFrame.from_dict(results)
	return results_frame


def main():
	parser = setup_parser()
	args = parser.parse_args()
	altHypDir = args.altHypDir
	nullHypDir = args.nullHypDir
	outfile = args.outfile
	# foregroundClade = args.foregroundClade
	# backgroundClade = args.backgroundClade
	rates = parse_rates(altHypDir, nullHypDir)
	rates.to_csv(outfile, index=False)


if __name__ == "__main__":
	import sys, argparse, glob, re
	import pandas as pd
	import numpy as np
	from scipy import stats
	from os import path, listdir, chdir
	main()

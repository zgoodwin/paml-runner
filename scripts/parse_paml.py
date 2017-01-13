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
	lineCount = 0
	matrix = {}
	with open(m, 'r') as f:
		for line in f:
			if lineCount != 0:
				l = line.strip().split()
				l.append('0.000')
				species = l[0]
				values = l[1:(len(l))]
				newValues = []
				for v in values:
					newValues.append(float(v))
				matrix[species] = newValues
			lineCount += 1
	return matrix


def replaceZeros(nums):
	newNums = []
	for i in nums:
		if i == 0.0:
			newNums.append(999.99)
		else:
			newNums.append(i)
	# Recast the new numbers array as a numerical matrix
	newNums = np.array(newNums)
	return newNums


def divideMatrix(matrixA, matrixB):
	quotientMatrix = {}
	for k in matrixA.keys():
		# Convert the rows to numpy array objects
		rowMatrixA = np.array(matrixA[k], dtype=np.float32)
		rowMatrixB = np.array(matrixB[k], dtype=np.float32)
		# Replace zeros in matrixB to avoid divide-by-zero errors
		nonZeroRowMatrixB = replaceZeros(rowMatrixB)
		quotientRow = rowMatrixA/nonZeroRowMatrixB
		quotientMatrix[k] = quotientRow
	return quotientMatrix


def symmetrize(m):
	keys = []
	max_length = len(m.keys())
	final = zeros(max_length, dtype=float32)
	symDict = {}
	for sorted_key in sorted(m, key=lambda k: len(m[k])):
		keys.append(sorted_key)
		numbers = m[sorted_key]
		n = numbers
		if len(numbers) < max_length:
			nZeros = max_length - len(numbers)
			n = append(numbers, zeros(nZeros, dtype=float32))

		final = vstack([final, n])

	# Now, remove the top row
	final = final[1:(max_length+1)]
	
	# Symmetrize by transposing the lower diagonal, then adding
	sym = (final + final.T)

	# convert the numpy array into a dictionary with species names
	for i in range(0,len(keys)):
		symDict[keys[i]] = sym[i]

	return symDict


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
	print(str(chi2_stat))
	pval = 2 * (stats.chi2.pdf(chi2_stat, df))
	return float(pval)

# This is supposed to return a pandas data frame with the LRT results.
# NOTE: You need to add error handling here because sometimes PAML does
#		not produce output files!! 
def parse_rates(altHypDir, nullHypDir):

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

	for d in altDirs:
		# genePath = "".join([altHypDir,d]) + "/"
		
		pamlFiles = glob.glob(d + '/*')
		# dnFile = "".join([genePath,"2NG.dN"])
		# dsFile = "".join([genePath,"2NG.dS"])

		dnFile = d + "/2NG.dN"
		dsFile = d + "/2NG.dS"

		# resultFile = "".join([genePath,geneName + "-H1.mlc"])
		resultFile = glob.glob(d + "/*.mlc")[0]
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
				print("The test failed on " + d)
				
				avgdN = 0.0
				avgdS = 0.0
				avgdNdS = 0.0

			geneList.append(d.split("/")[1])
			dnList.append(avgdN)
			dsList.append(avgdS)
			dndsList.append(avgdNdS)
			altLikelihoods.append(altLike)

	for d in nullDirs:
		resultFile = glob.glob(d + "/*.mlc")[0]
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

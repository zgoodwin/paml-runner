"""

Date Written: 1/12/2017
Last Edited:  1/12/2017

@author: Zane Goodwin
"""

def setup_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument("altHypDir", help="Path to the alternative hypothesis directory", type=str)
	parser.add_argument("nullHypDir", help="Path to the null hypothesis directory", type=str)
	parser.add_argument("outfile", help="Name of text file to store estimated rate values", type=str)
	parser.add_argument("pvalOutfile", help="Name of text file to store LRT likelihoods and pvalues", type=str)
	parser.add_argument("foregroundClade", help="Comma-separated list of species to include in the foreground clade rate calculations", type=str)
	parser.add_argument("backgroundClade", help="Comma-separated list of species to include in the background clade rate calculations", type=str)
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
	newNums = array(newNums)
	return newNums


def divideMatrix(matrixA, matrixB):
	quotientMatrix = {}
	for k in matrixA.keys():
		# Convert the rows to numpy array objects
		rowMatrixA = array(matrixA[k], dtype=float32)
		rowMatrixB = array(matrixB[k], dtype=float32)
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
		numElements += double(len(m[k]))
	averageRate = (total/numElements)
	return averageRate


def findLikelihood(mlcFile):
	lnLikeilihood = 0.0 # The natural log of the likelihood value
	with open(mlcFile,'r') as f:
		for line in f:
			if re.match('^lnL', line):
				lnLikelihood = line.split()[4]
	return float(lnLikelihood)



def parse_rates(altHypDir, nullHypDir):

	dnList = []
	dsList = []
	dndsList = []
	geneList = []
	altLikelihoods = []
	nullLikelihoods = []

	altDirs = os.listdir(altHypDir)
	nullDirs = os.listdir(nullHypDir)

	for d in altDirs:
		geneName = d

		genePath = "".join([altHypDir,d]) + "/"
		pamlFiles = os.listdir(genePath)

		dnFile = "".join([genePath,"2NG.dN"])
		dsFile = "".join([genePath,"2NG.dS"])
		resultFile = "".join([genePath,geneName + "-H1.mlc"])
		altLike = findLikelihood(resultFile)

		if os.path.isfile(dnFile):

			dnMatrix = parseDistanceMatrix(dnFile)
			dsMatrix = parseDistanceMatrix(dsFile)
			dndsMatrix = divideMatrix(dnMatrix, dsMatrix)
			
			try:
				avgdN = calculateAverageRates(dnMatrix)
				avgdS = calculateAverageRates(dsMatrix)
				avgdNdS = calculateAverageRates(dndsMatrix)

				fg = foregroundClade.split(",")
				bg = backgroundClade.split(",")

				calculateAverageCladeRates(dndsMatrix, fg, bg)

				
			except ZeroDivisionError:
				print "The test failed on " + geneName
				
				avgdN = 0.0
				avgdS = 0.0
				avgdNdS = 0.0

			geneList.append(geneName)
			dnList.append(avgdN)
			dsList.append(avgdS)
			dndsList.append(avgdNdS)
			altLikelihoods.append(altLike)

	for d in nullDirs:
		geneName = d
		genePath = "".join([nullHypDir,d]) + "/"
		pamlFiles = os.listdir(genePath)
		resultFile = "".join([genePath,geneName + "-H0.mlc"])
		nullLike = findLikelihood(resultFile)
		nullLikelihoods.append(nullLike)

	

	# Put the data into a pandas data frame
	results = {'genes' : geneList,
			   'dN' : dnList,
			   'dS' : dsList,
			   'dN_dS' : dndsList,
			   'alt_likelihood' : altLikelihoods,
			   'null_likelihood' : nullLikelihoods}
	results_frame = pd.DataFrame.from_dict(results)
	return results_frame

def main():
	parser = setup_parser()
	args = parser.parse_args()

	altHypDir = args.altHypDir
	nullHypDir = args.nullHypDir
	outfile = args.outfile
	pvalOutfile = args.pvalOutfile
	foregroundClade = args.foregroundClade
	backgroundClade = args.backgroundClade

	rates = parse_rates(altHypDir, nullHypDir)



if __name__ == "__main__":
	import sys, argparse, glob, re
	import pandas as pd
	import numpy
	from scipy import stats
	from os import path, listdir, chdir
	main()

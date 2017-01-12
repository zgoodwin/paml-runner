import sys
import os
import re
import argparse
from numpy import *
# import pandas as pd
from scipy import stats
# import matplotlib.pyplot as plt
# import matplotlib.ticker as mtick

font = {'size' : 12}
# plt.rc('font', **font)

##### FUNCTIONS #####

## Read text file containing dN and dS values
## Parameters:
#		m: Path to the dN or dS matrix file 
## Return value:
#		matrix: a dictionary containing the dN and dS value
# 
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


## Helper function for divideMatrix to replace zeros in a matrix with high numbers.
## Parameters:
#		nums: A numpy array of float32 values
## Return value:
#		newNums: A numpy array of float32 values with zeros replaced with 999.99
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


## Elementwise division of one matrix (matrixA) by another matrix (matrixB)
## Parameters:
#		matrixA: an NxM matrix of float values
#		matrixB: an NxM matrix of float values (same dimensions as matrixA)
## Return value:
#		An NxM matrix representing the result of matrixA/matrixB
# 
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

## Calculate average rates for a matrix
## Parameters:
#		m: A numpy matrix of dN or dS value
## Return values:
#		averageRate: a numpy matrix containing the average dN/dS rate for a gene,
#		across all species.
def calculateAverageRates(m):

	total = 0.0
	numElements = 0.0
	averageRate = 0.0
	for k in m.keys():

		total += sum(m[k])
		numElements += double(len(m[k]))

	averageRate = (total/numElements)

	return averageRate

## Symmetrize the lower half of a numpy matrix
## Parameters:
#		m: A dictionary with dN/dS rate values
## Return value:
#		symDict: a new dictionary object with the entire square matrix
def symmetrize(m):
	keys = []
	max_length = len(m.keys())
	final = zeros(max_length, dtype=float32)
	symDict = {}
	# secondRow = vstack([firstRow, ones(max_length, dtype=float32)])
	
	for sorted_key in sorted(m, key=lambda k: len(m[k])):
		keys.append(sorted_key)
		numbers = m[sorted_key]

		n = numbers
		if len(numbers) < max_length:
			nZeros = max_length - len(numbers)
			n = append(numbers, zeros(nZeros, dtype=float32))

		final = vstack([final, n])
	# 	# print array_str(n, max_line_width=2000)

	# Now, remove the top row
	final = final[1:(max_length+1)]
	
	# Symmetrize by transposing the lower diagonal, then adding
	sym = (final + final.T)
	
	# convert the numpy array into a dictionary with species names
	for i in range(0,len(keys)):
		symDict[keys[i]] = sym[i]

	return symDict

## Calculate average rates for a matrix, given a foreground clade and a background
#  clade
## Parameters:
#		m: a matrix of rates
#		groupA: a python list of species belonging to the foreground clade
#		groupB: a python list of species belonging to the background clade
#
## Return value:
#		rateDict: A dictionary containing the average rates for the foreground 
#		and the background 
#
# def calculateAverageCladeRates(m, groupA, groupB):
# 	rateDict = {"FG":0.0,"BG":0.0}

# 	sym_rate_matrix = symmetrize(m)

# 	aGroupRates = [m[i] for i in groupA]
# 	bGroupRates = [m[j] for j in groupB]

	

# 	sumGroupA = 0.0
# 	sumGroupB = 0.0


# 	# rateDict["FG"] = 
# 	# rateDict["BG"] = 
	
# 	return rateDict

## Locate the likelihood value in PAML's mlc file
## Parameters: 
#		mlcFile: Path to PAML's MLC, which contains the results of a PAML run.
## Return value:
#		lnLikelihood: The natural log of the likelihood value from PAML
def findLikelihood(mlcFile):

	lnLikeilihood = 0.0 # The natural log of the likelihood value
	with open(mlcFile,'r') as f:
		for line in f:
			if re.match('^lnL', line):
				lnLikelihood = line.split()[4]

	return float(lnLikelihood)

##### MAIN #####

# args = sys.argv
# altHypDir = args[1]
# nullHypDir = args[2]
# outfile = args[3]
# pvalOutfile = args[4]

parser = argparse.ArgumentParser()
parser.add_argument("altHypDir", help="Path to the alternative hypothesis directory", type=str)
parser.add_argument("nullHypDir", help="Path to the null hypothesis directory", type=str)
parser.add_argument("outfile", help="Name of text file to store estimated rate values", type=str)
parser.add_argument("pvalOutfile", help="Name of text file to store LRT likelihoods and pvalues", type=str)
parser.add_argument("foregroundClade", help="Comma-separated list of species to include in the foreground clade rate calculations", type=str)
parser.add_argument("backgroundClade", help="Comma-separated list of species to include in the background clade rate calculations", type=str)

args = parser.parse_args()
altHypDir = args.altHypDir
nullHypDir = args.nullHypDir
outfile = args.outfile
pvalOutfile = args.pvalOutfile
foregroundClade = args.foregroundClade
backgroundClade = args.backgroundClade


dirs = os.listdir(altHypDir)

dnList = []
dsList = []
dndsList = []
geneList = []
altLikelihoods = []
for d in dirs:
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
		
## Collect the directories containing the PAML runs for each gene		
dirs = os.listdir(nullHypDir)

## Collect the likelihood values for each gene
nullLikelihoods = []
for d in dirs:
	geneName = d
	genePath = "".join([nullHypDir,d]) + "/"
	pamlFiles = os.listdir(genePath)
	resultFile = "".join([genePath,geneName + "-H0.mlc"])
	nullLike = findLikelihood(resultFile)
	nullLikelihoods.append(nullLike)

## Draw the bar plot that plots dN, dS and dN/dS for each gene
# See : http://stackoverflow.com/questions/21810823/clustered-barchart-in-matplotlib
# for info on how to draw the bar plot
		
ind = arange(len(dndsList))
width = 0.35

# fig, ax = plt.subplots()

geneDict = {}
for i in range(len(geneList)):
	gene = geneList[i]
	dN = dnList[i]
	dS = dsList[i]
	dnds = dndsList[i]

	rates = [dN, dS, dnds]

	geneDict.update({ gene : rates })

## Write the p-values, alternative likelihood, null likelihood and 
pvalDict = {}
phandle = open(pvalOutfile, 'w')
for j in range(len(geneList)):
	gene = geneList[j]
	alt = altLikelihoods[j]
	null = nullLikelihoods[j]

	# Get the difference between the null and alternative likelihood values
	dist = (-2.0)*(null - alt)
	
	# Look up the p-value using the distance between the null and alternative likelihoods
	pval = 2 * (stats.chi2.pdf(dist, 2))

	# Reformat the likelihoods and p-values to have 8 significant figures
	alt =  "%.8f" % alt
	null =  "%.8f" % null
	dist =  "%.8f" % dist
	pval = "%.8f" % pval

	phandle.write("%s\t%s\t%s\t%s\t%s\n" % (gene, alt, null, dist, pval))

phandle.close()

## Set parameters for the dN/dS graph
initial_gap = 0.8
start = initial_gap
width = 1.0
gap = 0.9
colors="rgb"


## Write the dN, dS and dN/dS values for each gene to a TSV file,
#  and set up the labels for the dN/dS graph
out = open(outfile, 'w')
out.write("gene\tdN\tdS\tdN/dS\n")
for g in geneDict:
	
	rateList = ('dN', 'dS', 'dN/dS')
	size = len(geneDict[g])
	ind = linspace(start, start + width, size+1)[:-1]
	w = (ind[1]-ind[0])
	start = start + width + gap
	# plt.bar(ind, geneDict[g], w, color=list(colors[:size]), alpha = 0.3)

	numlist = ["%.8f" % n for n in geneDict[g]]
	out.write(g + "\t" + "\t".join(numlist) + "\n")

out.close()
	

## Draw the dN/dS graph
# tick_loc = ((arange(len(geneDict)) * (width + gap) + initial_gap + width/2)-0.5)
# ax.set_xticklabels([g for g in geneDict])
# ax.xaxis.set_major_locator(mtick.FixedLocator(tick_loc))
# plt.xticks(rotation=45)
# plt.ylabel('Rate Value')
# plt.xlabel('Gene')
# plt.title('dN, dS, and dN/dS in the EDC')
# plt.savefig("test.pdf", format='pdf')





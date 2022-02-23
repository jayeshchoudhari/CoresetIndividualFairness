import numpy as np
import scipy.stats as stats
import random
import sys



if __name__ == "__main__":

	if len(sys.argv) < 5: 
		print("1) inputFileName, 2) outputFileName, 3) num UAR Points, 4) num PowerLaw Points")
		sys.exit()
		
	inputFileName = sys.argv[1]
	outputFileName = sys.argv[2]
	numUARPointsToSample = int(sys.argv[3])
	finalNumOfPoints = int(sys.argv[4])

	powerLawAlpha = 1.5


	fin = open(inputFileName)

	allDataLines = []
	allUARLines = []

	for line in fin:
		if line[0] != "#":
			line = line.strip()
			allDataLines.append(line)

	fin.close()

	sampleIds = np.random.choice(range(len(allDataLines)), numUARPointsToSample)

	for eachId in sampleIds:
		allUARLines.append(allDataLines[eachId])


	print("Got UAR points...Now sampling with Power Law....")

	powerLawVals = [(x+1)**(-1.5) for x in range(numUARPointsToSample)]
	powerLawDist = [x/sum(powerLawVals) for x in powerLawVals]


	newSampleIds = np.random.choice(range(len(allUARLines)), finalNumOfPoints, p=powerLawDist)

	fout = open(outputFileName, 'w')

	for eachId in newSampleIds:
		fout.write(allUARLines[eachId] + "\n")

	fout.close()
	print("Done writing to file --- ", outputFileName, len(newSampleIds))
import numpy as np
import sys

if len(sys.argv) != 3:
	print("Error (Not enough input parameters):")
	print("1: FileName")
	print("2: Num of Points to sample")
	exit(1)


datasetFileName = sys.argv[1]
numPointsToSample = int(sys.argv[2])

fin = open(datasetFileName)

allDataLines = []

for line in fin:
	if line[0] != "#":
		line = line.strip()
		allDataLines.append(line)

fin.close()

sampleIds = np.random.choice(allDataLines, numPointsToSample)
## sampleIds will be data points (string)


print("Got ", len(sampleIds), " point ids... Storing to file...")

lastdotIndex = datasetFileName.rfind('.')

foutFileName = datasetFileName[:lastdotIndex] + "_" + str(numPointsToSample) + ".txt"

fout = open(foutFileName, "w")

# print(sampleIds)

for eachId in sampleIds:
	lineToWrite = eachId
	fout.write(lineToWrite + "\n")

fout.close()
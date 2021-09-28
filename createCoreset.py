import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
from kmedianInitialCenters import _kmeans_plusplus
import math
from sklearn.cluster import KMeans
from collections import defaultdict, Counter
from time import perf_counter
import json
import sys


if len(sys.argv) != 4:
	print("Error (Not enough input parameters):")
	print("1: Dataset Name (make sure that a folder named 'datasetNameCoreset' exists inside 'datasets' folder)")
	print("2: FileName")
	print("3: k(num Centers)")
	exit(1)


datasetName = sys.argv[1]
datasetFileName = sys.argv[2]
k = int(sys.argv[3])

datasetFolder = "datasets/" + datasetName + "Coreset"

percentageSizes = [0.5, 1.0, 2.0, 5.0, 7.0, 10.0, 15.0, 20.0]
numCopies = 5
coresetAlgos = ['kmedian', 'indFair', 'uni']
listOfPoints =  []

def localEuclideanDist(point_a, point_b):
	if len(point_a) != len(point_b):
		print("dim mismatch")
		exit(1)

	sqDist = 0
	for i in range(len(point_a)):
		sqDist +=  (point_a[i] - point_b[i])**2

	return math.sqrt(sqDist)


def getNearestCenter(biCriteriaCenters, point):
	# print(type(point), np.array(point).shape)
	# print(type(biCriteriaCenters[0]), np.array(biCriteriaCenters[0]).shape)

	minDist = localEuclideanDist(point, biCriteriaCenters[0])
	# print(minDist)

	centerId = 0

	for i in range(1, len(biCriteriaCenters)):
		idist = localEuclideanDist(point, biCriteriaCenters[i])
		if minDist > idist:
			minDist = idist
			centerId = i

	return [minDist, centerId]


def KMedianCost(biCriteriaCenters):
	totalCost = 0
	clusterAvgWiseCost = [0]*k
	clusterWiseNumPoints = [0]*k
	minDistOfAllPoints = [0]*len(listOfPoints)
	centerIdsOfAllPoints = [] 

	for i in range(len(listOfPoints)):
		iCost, iCenter = getNearestCenter(biCriteriaCenters, listOfPoints[i])
		totalCost += iCost 
		clusterAvgWiseCost[iCenter] += iCost
		clusterWiseNumPoints[iCenter] += 1
		minDistOfAllPoints[i] = iCost
		centerIdsOfAllPoints.append(iCenter)

	print(clusterAvgWiseCost)
	print(clusterWiseNumPoints)

	# average Cost...
	for i in range(len(clusterAvgWiseCost)):
		clusterAvgWiseCost[i] = clusterAvgWiseCost[i]/clusterWiseNumPoints[i]


	avgTotalCost = totalCost / len(listOfPoints)

	return [avgTotalCost, clusterAvgWiseCost, minDistOfAllPoints, centerIdsOfAllPoints]


def computeKMedianSensitivity(biCriteriaCenters):
	p = 1
	alphaVal = 2**(p+3) * (math.log(k) + 2)
	avgTotalCost, clusterWiseAvgCost, minDistOfAllPoints, centerIdsOfAllPoints = KMedianCost(biCriteriaCenters)

	# gives the size of each cluster
	eachClusterSize = Counter(centerIdsOfAllPoints)

	sensAllPoints = [0]*len(listOfPoints)

	for i in range(len(listOfPoints)):

		firstTerm = (alphaVal * (2**p) * minDistOfAllPoints[i])/(2 * avgTotalCost)
		secondTerm = (alphaVal * (4**p) * clusterWiseAvgCost[centerIdsOfAllPoints[i]]) / (4 * avgTotalCost)
		thirdTerm = ((4**p) * len(listOfPoints))/ eachClusterSize[centerIdsOfAllPoints[i]]

		sensAllPoints[i] = firstTerm +  secondTerm + thirdTerm

	return sensAllPoints


def IFFairCoresetSensitivity(kmedianSensAllPoints, epsilon=0.1):
	ifSensAllPoints = []
	epsVal = epsilon
	numPoints = len(listOfPoints)

	constTerm = 16 * k * math.log(numPoints) / (epsVal * numPoints)

	for i in range(len(listOfPoints)):
		tempSens_i = kmedianSensAllPoints[i]
		ifSensAllPoints.append(tempSens_i + constTerm)

	return ifSensAllPoints


def getProbDist(weightVector):
	normalization = sum(weightVector)
	probDist = [(x / normalization) for x in weightVector]
	return probDist


def writeCoresetToFile(fileName, numpyArray, weightVector):
	f = open(fileName, 'w')

	for i in range(len(numpyArray)):
		point_i_str = " ".join([str(x) for x in numpyArray[i]]) + " " + str(weightVector[i][0]) + "\n"
		f.write(point_i_str)

	f.close()


def writeToJsonFile(dictObj, fileName):
	jsonObj = json.dumps(dictObj)
	f = open(fileName, 'w')
	f.write(jsonObj)
	f.close()


f = open(datasetFileName, 'r')

## main part....
for line in f:
	flds = line.strip().split()
	point = [float(x) for x in flds]
	listOfPoints.append(point)

f.close()

n = len(listOfPoints)
dim = len(listOfPoints[0])

listOfPoints = np.array(listOfPoints)
# pairWiseDistances = computePairwiseDistances();

# Start the stopwatch / counter
kmSensComp_start = perf_counter()
centers, ind = _kmeans_plusplus(listOfPoints, 5, 1)
# print(centers)

kmedian = KMeans(n_clusters=k, random_state=0, max_iter=1, n_init = 1).fit(listOfPoints)
# print(kmedian.cluster_centers_)
biCriteriaCenters = kmedian.cluster_centers_

## Sensitivity Scores
kmedianSensAllPoints = computeKMedianSensitivity(biCriteriaCenters)
kmSensComp_end = perf_counter()

kmSensCompTime = kmSensComp_end - kmSensComp_start


ifSensComp_start = perf_counter()
indFairSensAllPoints = IFFairCoresetSensitivity(kmedianSensAllPoints)
ifSensComp_end = perf_counter()

ifSensCompTime = kmSensCompTime + (ifSensComp_end - ifSensComp_start)


## Sampling Distributions
kmedianSensProbDist = getProbDist(kmedianSensAllPoints)
indFairSensProbDist = getProbDist(indFairSensAllPoints)
uniProbDist = [1/n]*n


timeForEachMethod = defaultdict(lambda: defaultdict(float))

# fileName convention -- dataset_k_coreset_Method_size_copy

allPointsIds = range(n)

for eachPerc in percentageSizes:

	for copyId in range(numCopies):

		numPointsToSample = int((eachPerc * n)/100)

		kmSampleStart = perf_counter()
		kmedianSampleIds = np.random.choice(allPointsIds, numPointsToSample, p=kmedianSensProbDist)
		kmSampleEnd = perf_counter()

		kmCoresetTime = kmSensCompTime + (kmSampleEnd - kmSampleStart)

		ifSampleStart = perf_counter()
		indFairSampleIds = np.random.choice(allPointsIds, numPointsToSample, p=indFairSensProbDist)
		ifSampleEnd = perf_counter()
		
		ifCoresetTime = ifSensCompTime + (ifSampleEnd - ifSampleStart)


		uniSampleStart = perf_counter()
		uniSampleIds = np.random.choice(allPointsIds, numPointsToSample, p=uniProbDist)
		uniSampleEnd = perf_counter()

		uniCoresetTime = (uniSampleEnd - uniSampleStart)

		if copyId == numCopies - 1:
			timeForEachMethod[eachPerc]['KMedian'] = kmCoresetTime
			timeForEachMethod[eachPerc]['IndFair'] = ifCoresetTime
			timeForEachMethod[eachPerc]['Uni'] = uniCoresetTime

		kmSamplePoints = np.empty([numPointsToSample, dim])
		ifSamplePoints = np.empty([numPointsToSample, dim])
		uniSamplePoints = np.empty([numPointsToSample, dim])

		kmWeights = np.empty([numPointsToSample, 1])
		ifWeights = np.empty([numPointsToSample, 1])
		uniWeights = np.empty([numPointsToSample, 1])


		for i in range(numPointsToSample):
			kmSamplePoints[i] = listOfPoints[kmedianSampleIds[i]]
			kmWeights[i] = 1/(numPointsToSample * kmedianSensProbDist[kmedianSampleIds[i]])

			ifSamplePoints[i] = listOfPoints[indFairSampleIds[i]]
			ifWeights[i] = 1/(numPointsToSample * indFairSensProbDist[indFairSampleIds[i]])
			
			uniSamplePoints[i] = listOfPoints[uniSampleIds[i]]
			uniWeights[i] = 1/(numPointsToSample * uniProbDist[uniSampleIds[i]])

		print("Got samples for each coreset method")
		print("writing to file")

		kmFileName = datasetFolder + "/" + datasetName + "_" + str(k) + "_" + "coreset" + "_" + "KMedian" + "_" + str(eachPerc) + "_" + str(copyId) + ".txt"
		ifFileName = datasetFolder + "/" + datasetName + "_" + str(k) + "_" + "coreset" + "_" + "IndFair" + "_" + str(eachPerc) + "_" + str(copyId) + ".txt"
		uniFileName = datasetFolder + "/" + datasetName + "_" + str(k) + "_" + "coreset" + "_" + "Uni" + "_" + str(eachPerc) + "_" + str(copyId) + ".txt"


		writeCoresetToFile(kmFileName, kmSamplePoints, kmWeights)
		writeCoresetToFile(ifFileName, ifSamplePoints, ifWeights)
		writeCoresetToFile(uniFileName, uniSamplePoints, uniWeights)

		print("Done creating coreset for perc = ", eachPerc, " withCopyId = ", copyId)


writeToJsonFile(timeForEachMethod, datasetFolder + "/" + datasetName + "_" + str(k) + "_coresetCreationTimeForEachMethod.json")

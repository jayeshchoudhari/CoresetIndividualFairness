from __future__ import division
import sys
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
import pathlib



allDatasetNames = ['Bank', 'Diabetes', 'Census'] 

fileOfDName = {'Bank':'./outputForPlots/Bank_AllOutputForPlots.txt', 'Diabetes':'./outputForPlots/Diabetes_AllOutputForPlots.txt', 'Census':'./outputForPlots/Census_AllOutputForPlots.txt'}


kValues = [5, 10, 15]
algos = ['Uni', 'KMedian', 'IndFair']


barWidth = 0.25

	



# costDataForDname = defaultdict(list)
# fairDataForDname = defaultdict(list)
# timeDataForDname = defaultdict(list)

avgCostValuesAlgoKVal = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
stdCostValuesAlgoKVal = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

avgFairnessValuesAlgoKVal = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
stdFairnessValuesAlgoKVal = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

avgFairnessOnFulldataAlgoKVal = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
stdFairnessOnFulldataAlgoKVal = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
ogFairOnFullDataList = defaultdict(lambda: defaultdict(list))

timeValuesAlgoKVal = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
timeOnFullDataList = defaultdict(lambda: defaultdict(list))


for eachDName in allDatasetNames:

	fname = fileOfDName[eachDName]

	fin = open(fname)

	headerLine = fin.readline()

	# coresetTechnique coresetSize k avgRelErrOfCostOnFullData stdDevRelErrOfCostOnFullData avgDiffOfFairnessOnFullData stdDevDiffOfFairnessOnFullData avgTimeOnCoreset timeOnFullData

	for line in fin:
		flds = line.strip().split()
		
		algoName = flds[0].strip()
		corSize = float(flds[1].strip())
		kval = int(flds[2].strip())

		avgRelErrOfCostOnFullData = float(flds[3].strip())
		stdDevRelErrOfCostOnFullData = float(flds[4].strip())
		
		avgDiffOfFairnessOnFullData = float(flds[5].strip())
		stdDevDiffOfFairnessOnFullData = float(flds[6].strip())

		avgFairnessOnFulldata = float(flds[7].strip())
		stdFairnessOnFulldata = float(flds[8].strip())
		ogFairness = float(flds[9].strip())

		avgTimeOnCoreset = float(flds[10].strip()) * 10**(-6)
		timeOnFullData = float(flds[11].strip()) * 10**(-6)
		
		avgCostValuesAlgoKVal[kval][eachDName][algoName].append(avgRelErrOfCostOnFullData)
		stdCostValuesAlgoKVal[kval][eachDName][algoName].append(stdDevRelErrOfCostOnFullData)

		avgFairnessValuesAlgoKVal[kval][eachDName][algoName].append(avgDiffOfFairnessOnFullData)
		stdFairnessValuesAlgoKVal[kval][eachDName][algoName].append(stdDevDiffOfFairnessOnFullData)

		avgFairnessOnFulldataAlgoKVal[kval][eachDName][algoName].append(avgFairnessOnFulldata)
		stdFairnessOnFulldataAlgoKVal[kval][eachDName][algoName].append(stdFairnessOnFulldata)
		ogFairOnFullDataList[kval][eachDName].append(ogFairness)
		
		timeValuesAlgoKVal[kval][eachDName][algoName].append(avgTimeOnCoreset)
		timeOnFullDataList[kval][eachDName] = timeOnFullData

	fin.close()




print("Plotting....")

print(timeValuesAlgoKVal)


# '''
#### Err in Cost...
numRows = 3
fig = plt.figure(figsize=(35, 11))
fig.subplots_adjust(hspace=0.35, wspace=0.15)
plt.subplots_adjust(left=0.07, right=0.93)

# did = 1
rowCounter = 0
for ind in (range(len(kValues))):
# for kval in kValues:
	kval = kValues[ind]

	avgCostErrVals = defaultdict(lambda: defaultdict(list))
	stdCostErrVals = defaultdict(lambda: defaultdict(list))
	avgFairErrVals = defaultdict(lambda: defaultdict(list))
	stdFairErrVals = defaultdict(lambda: defaultdict(list))

	avgFairValsOnFullData = defaultdict(lambda: defaultdict(list))
	stdFairValsOnFullData = defaultdict(lambda: defaultdict(list))


	for eachDName in allDatasetNames:
		for algoName in algos:

			avgCostErrVals[eachDName][algoName] = avgCostValuesAlgoKVal[kval][eachDName][algoName]
			stdCostErrVals[eachDName][algoName] = stdCostValuesAlgoKVal[kval][eachDName][algoName]

			avgFairErrVals[eachDName][algoName] = avgFairnessValuesAlgoKVal[kval][eachDName][algoName]
			stdFairErrVals[eachDName][algoName] = stdFairnessValuesAlgoKVal[kval][eachDName][algoName]

			avgFairValsOnFullData[eachDName][algoName] = avgFairnessOnFulldataAlgoKVal[kval][eachDName][algoName]
			stdFairValsOnFullData[eachDName][algoName] = stdFairnessOnFulldataAlgoKVal[kval][eachDName][algoName]


	print(avgCostErrVals)

	did = 1
	for eachDName in allDatasetNames:

		if eachDName == 'Bank':
			coresetSizes = [2.0, 5.0, 7.0, 10.0]
		else:
			coresetSizes = [0.5, 1.0, 2.0, 5.0]


		i = 0
		print("plotting -- ", did + (rowCounter*len(allDatasetNames)))
		ax = fig.add_subplot(numRows, len(allDatasetNames), did + (rowCounter*len(allDatasetNames)))

		for algoName in algos:
			# Cost Error Plot
			r1 = [(x + i*barWidth) for x in range(len(coresetSizes))]
			ax.bar(r1, avgCostErrVals[eachDName][algoName], yerr=stdCostErrVals[eachDName][algoName], width=barWidth, edgecolor='white', label=algoName + 'Cor', capsize=3)
			i = i + 1
			
			# ax.tick_params(axis="x", labelsize=20)
			ax.tick_params(axis="y", labelsize=20)

			print(eachDName, " ", algoName, " --- " , avgCostErrVals[eachDName][algoName])


		# ax.legend(fontsize=15)
		ax.set_xticks([r + barWidth for r in range(len(coresetSizes))])
		ax.set_xticklabels(coresetSizes, fontsize=20)
		ax.tick_params(axis="x", labelsize=20)


		titleStr = eachDName + "(k=" + str(kval) + ")"
		ax.set_title(titleStr, fontsize = 22)

		if did == 1:
			ax.set_ylabel('% Relative Error', fontsize=20)

		if ind == len(kValues) - 1:
			ax.set_xlabel('% Coreset Size', fontsize=20)

		did += 1

		

	rowCounter += 1

newHandles, newLabels = ax.get_legend_handles_labels()
fig.legend(newHandles, newLabels, loc= 'upper center', ncol=3, fontsize=20)
plt.show()

# '''

'''
#### Avg Fairness...
numRows = 3
fig = plt.figure(figsize=(35, 11))
fig.subplots_adjust(hspace=0.35, wspace=0.15)
plt.subplots_adjust(left=0.07, right=0.93)

# did = 1
rowCounter = 0
for ind in (range(len(kValues))):
# for kval in kValues:
	kval = kValues[ind]

	avgCostErrVals = defaultdict(lambda: defaultdict(list))
	stdCostErrVals = defaultdict(lambda: defaultdict(list))
	avgFairErrVals = defaultdict(lambda: defaultdict(list))
	stdFairErrVals = defaultdict(lambda: defaultdict(list))

	avgFairValsOnFullData = defaultdict(lambda: defaultdict(list))
	stdFairValsOnFullData = defaultdict(lambda: defaultdict(list))


	for eachDName in allDatasetNames:
		for algoName in algos:

			avgCostErrVals[eachDName][algoName] = avgCostValuesAlgoKVal[kval][eachDName][algoName]
			stdCostErrVals[eachDName][algoName] = stdCostValuesAlgoKVal[kval][eachDName][algoName]

			avgFairErrVals[eachDName][algoName] = avgFairnessValuesAlgoKVal[kval][eachDName][algoName]
			stdFairErrVals[eachDName][algoName] = stdFairnessValuesAlgoKVal[kval][eachDName][algoName]

			avgFairValsOnFullData[eachDName][algoName] = avgFairnessOnFulldataAlgoKVal[kval][eachDName][algoName]
			stdFairValsOnFullData[eachDName][algoName] = stdFairnessOnFulldataAlgoKVal[kval][eachDName][algoName]


	# print(avgCostErrVals)

	did = 1
	for eachDName in allDatasetNames:

		if eachDName == 'Bank':
			coresetSizes = [2.0, 5.0, 7.0, 10.0]
		else:
			coresetSizes = [0.5, 1.0, 2.0, 5.0]


		i = 0
		print("plotting -- ", did + (rowCounter*len(allDatasetNames)))
		ax = fig.add_subplot(numRows, len(allDatasetNames), did + (rowCounter*len(allDatasetNames)))

		for algoName in ['IndFair']:
			# Cost Error Plot
			r1 = [(x + i*barWidth) for x in range(len(coresetSizes))]
			# ax.bar(r1, avgCostErrVals[eachDName][algoName], yerr=stdCostErrVals[eachDName][algoName], width=barWidth, edgecolor='white', label=algoName, capsize=3)

			# ax.plot(range(len(coresetSizes)), avgFairValsOnFullData[eachDName][algoName], yerr=stdFairValsOnFullData[eachDName][algoName], label=algoName)
			ax.plot(range(len(coresetSizes)), avgFairValsOnFullData[eachDName][algoName], 'g-*', label=algoName, markersize=15)
			ax.plot(range(len(coresetSizes)), [ogFairOnFullDataList[kval][eachDName] for i in range(len(r1))], 'ko-', label='Orig. Fairness', markersize=15)

			if did == 1 and rowCounter == 0:
				finalLabels = [algoName + 'Cor' + ' Fairness', 'MV Fairness']

			i = i + 1
			
			ax.tick_params(axis="y", labelsize=20)
			ax.tick_params(axis="x", labelsize=20)
			ax.set_xticks([r for r in range(len(coresetSizes))])
			ax.set_xticklabels(coresetSizes, fontsize=20)

			# print(eachDName, " ", algoName, " --- " , avgCostErrVals[eachDName][algoName])

			if did == 1 and rowCounter == 0:
				newHandles, newLabels = ax.get_legend_handles_labels()
				print("Getting handles.....", newLabels)


		# ax.legend(fontsize=15)
		# ax.set_xticks([r + barWidth for r in range(len(coresetSizes))])
		# ax.set_xticklabels(coresetSizes, fontsize=20)
		# ax.tick_params(axis="x", labelsize=20)


		titleStr = eachDName + "(k=" + str(kval) + ")"
		ax.set_title(titleStr, fontsize = 22)

		if did == 1:
			ax.set_ylabel('Avg. Max. Fairness', fontsize=20)

		if ind == len(kValues) - 1:
			ax.set_xlabel('% Coreset Size', fontsize=20)

		did += 1

		

	rowCounter += 1

fig.legend(newHandles, finalLabels, loc= 'upper center', ncol=3, fontsize=20)
plt.show()
'''


'''
#### Avg Time...
numRows = 3
fig = plt.figure(figsize=(35, 11))
fig.subplots_adjust(hspace=0.35, wspace=0.15)
plt.subplots_adjust(left=0.07, right=0.93)

# did = 1
rowCounter = 0
for ind in (range(len(kValues))):
# for kval in kValues:
	kval = kValues[ind]

	avgCostErrVals = defaultdict(lambda: defaultdict(list))
	stdCostErrVals = defaultdict(lambda: defaultdict(list))
	avgFairErrVals = defaultdict(lambda: defaultdict(list))
	stdFairErrVals = defaultdict(lambda: defaultdict(list))

	avgFairValsOnFullData = defaultdict(lambda: defaultdict(list))
	stdFairValsOnFullData = defaultdict(lambda: defaultdict(list))

	avgTimeValsOnFullData = defaultdict(lambda: defaultdict(list))
	timeOnFullData = defaultdict(list)

	for eachDName in allDatasetNames:
		for algoName in algos:

			avgTimeValsOnFullData[eachDName][algoName] = timeValuesAlgoKVal[kval][eachDName][algoName]
			timeOnFullData[eachDName] = timeOnFullDataList[kval][eachDName]


	print(avgTimeValsOnFullData)

	did = 1
	for eachDName in allDatasetNames:

		if eachDName == 'Bank':
			coresetSizes = [2.0, 5.0, 7.0, 10.0]
		else:
			coresetSizes = [0.5, 1.0, 2.0, 5.0]


		i = 0
		print("plotting -- ", did + (rowCounter*len(allDatasetNames)))
		ax = fig.add_subplot(numRows, len(allDatasetNames), did + (rowCounter*len(allDatasetNames)))

		for algoName in ['IndFair']:

			print('indfair:', avgTimeValsOnFullData[eachDName][algoName])
			print('og:', timeOnFullData[eachDName])

			ax.plot(range(len(coresetSizes)), avgTimeValsOnFullData[eachDName][algoName], 'g-*', label=algoName + 'Cor', markersize=15)
			ax.plot(range(len(coresetSizes)), [timeOnFullData[eachDName] for i in range(len(coresetSizes))], 'ko-', label='Orig. Fairness', markersize=15)

			if did == 1 and rowCounter == 0:
				finalLabels = [algoName + 'Cor' + ' Time', 'MV Time on full dataset']

			i = i + 1
			
			ax.tick_params(axis="y", labelsize=20)
			ax.tick_params(axis="x", labelsize=20)
			ax.set_xticks([r for r in range(len(coresetSizes))])
			ax.set_xticklabels(coresetSizes, fontsize=20)

			# print(eachDName, " ", algoName, " --- " , avgCostErrVals[eachDName][algoName])

			if did == 1 and rowCounter == 0:
				newHandles, newLabels = ax.get_legend_handles_labels()
				print("Getting handles.....", newLabels)


		# ax.legend(fontsize=15)
		# ax.set_xticks([r + barWidth for r in range(len(coresetSizes))])
		# ax.set_xticklabels(coresetSizes, fontsize=20)
		# ax.tick_params(axis="x", labelsize=20)


		titleStr = eachDName + "(k=" + str(kval) + ")"
		ax.set_title(titleStr, fontsize = 22)

		if did == 1:
			ax.set_ylabel('Time (secs)', fontsize=20)

		if ind == len(kValues) - 1:
			ax.set_xlabel('% Coreset Size', fontsize=20)

		did += 1

		

	rowCounter += 1

fig.legend(newHandles, finalLabels, loc= 'upper center', ncol=3, fontsize=20)

plt.show()
'''
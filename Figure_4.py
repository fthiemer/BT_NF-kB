from Bielefeld_dataset_basic_python_dicts import *
from scipy.stats import multinomial
from itertools import permutations

def openReadableCSV(filename):#muss für jede Iteration neu kreiert werden
    readableCSV=csv.reader(open(str(filename),newline=''),delimiter='\t')
    return readableCSV
def openWriteableCSV(filename):
    writeableCSV=csv.writer(open(str(filename),'w+',newline=''),delimiter='\t')
    return writeableCSV
	
#create summarized Dicts which only hold BS numbers of Promoter/Gene
promoterDict = get_prom_BS_sum_dict()
geneDict = {}
for key in promoterDict:
	geneDict[extractEPDGeneName(key)] = [[]]
#put gene info in first entry of the list. Promoter BS Lists in the second entry
for key in promoterDict:
	if len(geneDict[extractEPDGeneName(key)]) < 2:
		geneDict[extractEPDGeneName(key)].insert(0, promoterDict[key]["geneInfo"])
	geneDict[extractEPDGeneName(key)][1].append(promoterDict[key]["BS_data"])
summarizedGeneDict = {}
for key in geneDict:
	summarizedGeneDict[key] = {"RELA":0, "RELB":0, "REL":0, "all":0}
for key in geneDict:
	for promoterBSList in geneDict[key][1]: #for every promoter in the gene
		for BS in promoterBSList: #for every BS in the promoter
			summarizedGeneDict[key]["all"] += 1
			summarizedGeneDict[key][BS["CRE"]] += 1
summarizedPromoterDict = {}
for key in promoterDict:
	summarizedPromoterDict[key] = {"RELA":0, "RELB":0, "REL":0, "all":0}
for key in promoterDict:
	for BS in promoterDict[key]["BS_data"]:
		summarizedPromoterDict[key]["all"] += 1
		summarizedPromoterDict[key][BS["CRE"]] += 1

#print data for 2 a) (old version mit ooen)
printDict2a = {"A":0, "B":0, "C":0, "AB":0, "AC":0, "BC":0, "ABC":0}
for promoter in summarizedPromoterDict:
	classificationString = ""
	for CRE in summarizedPromoterDict[promoter]:
		if summarizedPromoterDict[promoter][CRE] and CRE != "all":
			if CRE[-1] == "L":
				classificationString += "C"
			else:
				classificationString += CRE[-1]
	printDict2a[classificationString] += 1
print("2 a)")
for key in printDict2a:
	print(key, " ", str(printDict2a[key]))
#print data for 2 b)
printDict2b = {"A":0, "B":0, "C":0, "AB":0, "AC":0, "BC":0, "ABC":0}
for gene in summarizedGeneDict:
	classificationString = ""
	for CRE in summarizedGeneDict[gene]:
		if summarizedGeneDict[gene][CRE] and CRE != "all":
			if CRE[-1] == "L":
				classificationString += "C"
			else:
				classificationString += CRE[-1]
	printDict2b[classificationString] += 1
print("2 b)")
for key in printDict2b:
	print(key, " ", str(printDict2b[key]))

#2 c) and d)
def pmf(BSList):
	return multinomial.pmf(BSList, n = sum(BSList), p = len(BSList)*[1/len(BSList)])
def get_predominance_string(numOfBSDict, pValue):
	#returns string with predominance
	returnPieces = ["A", "B", "C"]
	tmpBSList = [numOfBSDict["RELA"], numOfBSDict["RELB"], numOfBSDict["REL"]]
	if pmf(tmpBSList) > pValue:
		return "ABC"
	else:
		for index, CRE in enumerate(returnPieces):
			indexList = [0,1,2]
			indexList.remove(index)
			if tmpBSList.index(max(tmpBSList)) == index: #if one of the max BS site values is the current one
				if(pmf([tmpBSList[index], tmpBSList[indexList[0]]]) <= pValue and 
				pmf([tmpBSList[indexList[1]], tmpBSList[indexList[0]]]) <= pValue and
				pmf([tmpBSList[index], tmpBSList[indexList[1]]]) > pValue):
					return "".join(sorted(CRE + returnPieces[indexList[1]]))
				elif( pmf([tmpBSList[index], tmpBSList[indexList[1]]]) <= pValue and 
				pmf([tmpBSList[indexList[0]], tmpBSList[indexList[1]]]) <= pValue and
				pmf([tmpBSList[index], tmpBSList[indexList[0]]]) > pValue):
					return "".join(sorted(CRE + returnPieces[indexList[0]]))
				elif tmpBSList[index] > tmpBSList[indexList[0]] and tmpBSList[index] > tmpBSList[indexList[1]]:
					if(pmf([tmpBSList[index], tmpBSList[indexList[0]]]) <= pValue and 
					pmf([tmpBSList[index], tmpBSList[indexList[0]]]) <= pValue):
						return CRE
					else:
						return "ABC"
				elif (pmf([tmpBSList[index], tmpBSList[indexList[0]]]) > pValue and 
					pmf([tmpBSList[index], tmpBSList[indexList[0]]]) > pValue):
					return "ABC"
				else:
					print("special case which was not considered")
					input()

print("2 c)")
printDict2c = {"A":0, "B":0, "C":0, "AB":0, "AC":0, "BC":0, "ABC":0}
for key in summarizedPromoterDict:
	printDict2c[get_predominance_string(summarizedPromoterDict[key], 0.05)] += 1
for key in printDict2c:
	print(key, " ", str(printDict2c[key]))
	
print("2 d)")
printDict2d = {"A":0, "B":0, "C":0, "AB":0, "AC":0, "BC":0, "ABC":0}
for key in summarizedGeneDict:
	printDict2d[get_predominance_string(summarizedGeneDict[key], 0.05)] += 1
for key in printDict2d:
	print(key, " ", str(printDict2d[key]))
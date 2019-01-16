from Bielefeld_dataset_basic_python_dicts import *

# create folder to store file in
if not os.path.exists(os.path.join(".", "data")):
    os.mkdir(os.path.join(".", "data"))

with open(os.path.join(".", "data", "summary_for_t-tests.txt"), "w+", newline = "") as f:
    tTestFile = csv.writer(f, delimiter = "\t")
    geneDict = get_gene_dict_with_QS()
    tTestDict = {}
    for QS in geneDict["IGF2"]["QSs"]:
        tTestDict[str(QS) + "_N"] = []
        tTestDict[str(QS) + "_Y"] = []
    for QS in geneDict["IGF2"]["avgHighestX"]:
        tTestDict[str(QS) + "_N"] = []
        tTestDict[str(QS) + "_Y"] = []
    tTestDict["num_of_promoters_N"] = []
    tTestDict["num_of_promoters_Y"] = []
    tTestFile.writerow([i for i in tTestDict])
    for gene in geneDict:
        for QS in geneDict[gene]["QSs"]:
            tTestDict[str(QS) + "_" + str(geneDict[gene]["geneInfo"]["exp_val"][0])].append(geneDict[gene]["QSs"][QS])
        for QS in geneDict[gene]["avgHighestX"]:
            tTestDict[str(QS) + "_" + str(geneDict[gene]["geneInfo"]["exp_val"][0])].append(geneDict[gene]["avgHighestX"][QS])
        tTestDict["num_of_promoters_" +str(geneDict[gene]["geneInfo"]["exp_val"][0])].append(int(geneDict[gene]["geneInfo"]["num_of_promoters"]))
    emptyColumns = 0
    num_of_categories = len(tTestDict.keys())
    print("num of categories are: " + str(num_of_categories))
    while emptyColumns != num_of_categories:
        emptyColumns = 0
        rowToWrite = []
        for subdict in [tTestDict[i] for i in tTestDict.keys()]:
            if subdict == []:
                emptyColumns += 1
                rowToWrite.append("")
            else:
                rowToWrite.append(subdict.pop())
        if emptyColumns != num_of_categories:
            tTestFile.writerow(rowToWrite)
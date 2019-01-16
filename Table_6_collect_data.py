#gathers info about the Target Genes predicted by QS_avg_10_highest_BS
from Bielefeld_dataset_basic_python_dicts import *

geneDict = get_gene_dict_with_best_QS()
def sort_key(item):
    return item[2]

#returns list with highest 10 BS in a gene
def hi10(gene, geneDict = geneDict):
    return [(str(geneDict[gene]["BS_data"][i]["CRE"])) for i in range(10)]
avg10List = sorted([
            (gene,
            ["RelA: " + str(hi10(gene).count("RELA")),"RelB: " + str(hi10(gene).count("RELB")),
            "c-Rel: " + str(hi10(gene).count("REL"))],
            geneDict[gene]["avgHighestX"]["avg10"], 
            geneDict[gene]["geneInfo"]["exp_val"],
            geneDict[gene]["BS_data"][0]["Chromosome"]) for gene in geneDict.keys()]
             ,reverse = True, key = sort_key)

def sort_keyW(item):
    return item[1]

Wang1List = sorted([
            (gene,
            geneDict[gene]["best_QS"]["Wang1"], 
            geneDict[gene]["geneInfo"]["exp_val"],
            geneDict[gene]["BS_data"][0]["Chromosome"])
             for gene in geneDict.keys()
             ], reverse = True, key = sort_keyW)

print("QS_avg_10_highest_BS data:")
print(avg10List[:10])
print("QS_Wang data:")
print(Wang1List[:11])

for X in [6995]:
    print("top " +str(X))
    avg10ExpVal = 0
    nosumavg10 = 0
    for i, tuple in enumerate(avg10List):
        # if tuple[2]>=405 and tuple[3] == "Yes":
            # avg10ExpVal += 1
        # if tuple[2]>=405 and tuple[3] == "No":
            # nosumavg10 += 1
        if i < X and tuple[3] == "Yes":
            avg10ExpVal += 1
        if i < X and tuple[3] == "No":
            nosumavg10 += 1

    print(avg10ExpVal)
    print(nosumavg10)

    Wang1ExpVal = 0
    nosumwang = 0
    for i, tuple in enumerate(Wang1List):
        # if tuple[1]> 0 and tuple[2] == "Yes":
            # Wang1ExpVal += 1
        # elif(tuple[1]> 0 and tuple[2] == "No"):
            # nosumwang += 1
        if i < X and tuple[2] == "Yes":
            Wang1ExpVal += 1
        elif(i < X and tuple[2] == "No"):
            nosumwang += 1

    print(Wang1ExpVal)
    print(nosumwang)

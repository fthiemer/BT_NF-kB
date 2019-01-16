# creates Figure 8 and prints used data to the terminal
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os 
from Bielefeld_dataset_basic_python_dicts import *

# create folder to store Figure in
if not os.path.exists(os.path.join(".", "Figures")):
    print("creating Folder for histogram")
    os.mkdir(os.path.join(".", "Figures"))

geneDict = get_gene_dict_with_best_QS()
#open file and assign gene name in a set and chromosome of not before saved genes in a list
BSFile = openReadableCSV(os.path.join("source_files", "all_data.txt"))
usedGeneNames = set()
chromosomeList = []
for row in BSFile:
    if BSFile.line_num != 1:
        if extractEPDGeneName(row[0]) not in usedGeneNames:
            usedGeneNames.add(extractEPDGeneName(row[0]))
            chromosomeList.append(row[3])

# create chrLabelList
chrLabelList = []
for i in range(1, 23):
    chrLabelList.append("chr"+str(i))
for i in ["chrX", "chrY"]:
    chrLabelList.append(i)
# create dict to map indexes of QS lists to the chromosomes
chrIndexDict = {}
for i, c in enumerate(chrLabelList):
    chrIndexDict[c] = i
avgAvg10Score = []
avgWang1Score = []
wang1Frequencies, avg10Frequencies = ([[] for i in range(24)], [[] for i in range(24)])
for gene in geneDict:
    chrIndex = chrIndexDict[geneDict[gene]["BS_data"][0]["Chromosome"]]
    avgAvg10Score.append(geneDict[gene]["avgHighestX"]["avg10"])
    avgWang1Score.append(geneDict[gene]["best_QS"]["Wang1"])
    avg10Frequencies[chrIndex].append(geneDict[gene]["avgHighestX"]["avg10"])
    wang1Frequencies[chrIndex].append(geneDict[gene]["best_QS"]["Wang1"])
    
avgAvg10Score = mean(avgAvg10Score)
print("average QS_avg_10_highest_BS Score over all chromosomes")
print(avgAvg10Score)
print("average QS_Wang1 over all chromosomes")
avgWang1Score = mean(avgWang1Score)
print(avgWang1Score)



for i in range(24):
    avg10Frequencies[i] = list(filter(lambda x: x >= avgAvg10Score, avg10Frequencies[i]))
    avg10Frequencies[i] = len(avg10Frequencies[i])
    wang1Frequencies[i] = list(filter(lambda x: x > 0, wang1Frequencies[i]))
    wang1Frequencies[i] = len(wang1Frequencies[i])
print("Number of genes with QS_avg_10_highest_BS higher than its average on all chromosomes")
print(sum(avg10Frequencies))

# create bar chart for fig 6
fig, ax = plt.subplots()
ax.set(xlabel = "Chromosome", ylabel = "Frequency of Genes")
# set xticks
plt.xticks(np.arange(0.5, 48, 2), chrLabelList, rotation = 60, size = "xx-small")
frequencies = []
exp_valFrequencies = []
for i in chrLabelList:
    frequencies.append(chromosomeList.count(i))
    exp_valFrequencies.append(0)

#print information for table 5 and description of Figure 6
print(frequencies)
print("Frequencies of QS_avg_10_highest_BS per chromosome")
print(avg10Frequencies)
print("number of genes with QS_avg_10_highest_BS over mean")
print(sum(avg10Frequencies))
print("mean of QS_avg_10_highest_BS predicted frequencies over all chromosomes")
print(mean(avg10Frequencies))
print("frequencies of target genes predicted by QS_Wang1 (value above mean) per chromosome")
print(wang1Frequencies)
print("absolute number of target genes predicted by QS_Wang1 on all chromosomes")
print(sum(wang1Frequencies))
print("average frequency of QS_Wang1 predicted target genes per chromosome.")
print(mean(wang1Frequencies))
print("frequencies_avg10/frequencies")
relAvg10 = [round(avg10Frequencies[i]/frequencies[i], 3) for i in range(24)]
print(relAvg10)
print("average relative frequency of QS_avg_10_highest_BS predicted target genes")
print(mean(relAvg10))

relWang1 = [round(wang1Frequencies[i]/frequencies[i], 3) for i in range(24)]
print("frequencies_QS_Wang1_predicted/frequency_on_chr")
print(relWang1)
print("average of frequencies_QS_Wang1_predicted/frequency_on_chr")
print(mean(relWang1))


for gene in geneDict:
    if geneDict[gene]["geneInfo"]["exp_val"] == "Yes":
        exp_valFrequencies[chrLabelList.index(geneDict[gene]["BS_data"][0]["Chromosome"])] += 1
print("chromosomal distribution of verified genes")
print(exp_valFrequencies)
print("proportion of verified genes of genes on chr in bi dataset")
print([round(exp_valFrequencies[i]/frequencies[i], 3) for i in range(24)])
plt.bar(np.arange(0.5, 48, 2), frequencies, width = 2, label = "Target genes per chromosome", 
color = "#202020cc", edgecolor = "#202020cc")
plt.bar(np.arange(0, 48, 2), wang1Frequencies, width = 1, label = "Genes with QS_Wang_1 higher than 0.",
color = "#cc17177d", edgecolor = "#cc17176d")
plt.bar(np.arange(1, 48, 2), avg10Frequencies, width = 1, label = "Genes with QS_avg_10_highest_BS higher than mean.",
color = "#17cccc7d", edgecolor = "#17cccc6d")
plt.bar(np.arange(0.5, 48, 2), exp_valFrequencies, width = 2, label = "Verified target genes.",
color = "#888888cc", edgecolor = "#404040cc")
plt.legend()
# plt.show()
fig.savefig(os.path.join(".", "Figures","BT_figure8.png"))
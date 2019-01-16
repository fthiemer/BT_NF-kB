from Bielefeld_dataset_basic_python_dicts import *

# create folder to store Figure in
if not os.path.exists(os.path.join(".", "data")):
    os.mkdir(os.path.join(".", "data"))

# creates datafile for Table 4
with open(os.path.join(".", "data", "QS_avgHighestX_comparison_motifspecific.txt"), "w+", newline = "") as f:
    tryoutFile = csv.writer(f, delimiter = "\t")
    geneDict = get_gene_dict_with_QS()
    header = ["gene", "promoters", "Ensembl_ID", "verified?", "BSs_RelA", "BSs_RelB", "BSs_c-Rel", "BSs_total", "QS_avg_10_highest_BS(all BSs)"]
    header.extend([motif + "_avg_"+str(i) for motif in ["RelA", "RelB", "Rel"] for i in range(2,16)])
    tryoutFile.writerow(header)
    for gene in geneDict:
        rowToWrite = [gene, 
                      geneDict[gene]["geneInfo"]["num_of_promoters"],
                      geneDict[gene]["geneInfo"]["ENSEMBL_ID"],
                      geneDict[gene]["geneInfo"]["exp_val"]
                      ]
        for motifType in ["RELA", "RELB", "REL", "all"]:   
            rowToWrite.append(geneDict[gene]["QSs"]["BS_sum_all_motifs"][motifType])
        try:
            rowToWrite.append(geneDict[gene]["avgHighestXmsp"]["avg10"]["all"])
        except:
            print(geneDict[gene]["avgHighestXmsp"])
            input()
        for motifType in ["RELA", "RELB", "REL"]:
            for number in range(2, 16):
                rowToWrite.append(geneDict[gene]["avgHighestXmsp"]["avg" + str(number)][motifType])
        tryoutFile.writerow(rowToWrite)
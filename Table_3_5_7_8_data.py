from Bielefeld_dataset_basic_python_dicts import *

# used for the binary logistic regression displayed in Table 3 and 5
with open(os.path.join(".", "important_source_files", "Bachelor_Thesis_Frederick_Thiemer_dataset", "All_infos_gene_summary.txt"), "w+", newline = "") as f:
	thresholdFile = csv.writer(f, delimiter = "\t")
	thresholdDict = {}
	geneDict = get_gene_dict_with_QS()
	for gene in geneDict.keys():
		thresholdDict[gene] = [gene, geneDict[gene]["geneInfo"]["ENSEMBL_ID"], 
								geneDict[gene]["geneInfo"]["exp_val"],
								geneDict[gene]["BS_data"][0]["Chromosome"]
								]
	QSList = list(geneDict["IGF2"]["QSs"].keys())
	header = ["geneName", "EnsemblID", "exp_val", "Chromosome"]
	header.append("QS_avg_10_highest_BS")
	header.extend(["QS_" + str(i) for i in QSList])
	for gene in geneDict.keys():
			thresholdDict[gene].append(geneDict[gene]["avgHighestX"]["avg10"])
			for QS in [geneDict[gene]["QSs"][j] for j in geneDict[gene]["QSs"].keys()]:
				thresholdDict[gene].append(QS)
	#add QSs without threshold
	header.append("QS_avg_10_highest_BS")
	header.extend(["QS_" + str(i)) for i in QSList])
	geneDict = get_gene_dict_with_QS()
	for gene in geneDict.keys():
		thresholdDict[gene].append(geneDict[gene]["avgHighestX"]["avg10"])
		for QS in [geneDict[gene]["QSs"][j] for j in geneDict[gene]["QSs"].keys()]:
			thresholdDict[gene].append(QS)
	# add QSs with threshold
	for threshold in range(240, 551, 10):
		header.append(str(threshold)+ "_QS_avg_10_highest_BS")
		header.extend([(str(threshold) + "_QS_" + str(i)) for i in QSList])
		geneDict = get_gene_dict_with_QS(threshold)
		for gene in geneDict.keys():
			thresholdDict[gene].append(geneDict[gene]["avgHighestX"]["avg10"])
			for QS in [geneDict[gene]["QSs"][j] for j in geneDict[gene]["QSs"].keys()]:
				thresholdDict[gene].append(QS)
	thresholdFile.writerow(header)
	for geneQSList in [thresholdDict[gene] for gene in thresholdDict.keys()]:
		thresholdFile.writerow(geneQSList)

# used for Table 7 and 8, includes QSs with optimal threshold applied as displayed in Table 5
with open(os.path.join(".", "important_source_files", "Bachelor_Thesis_Frederick_Thiemer_dataset", "gene_top_QSs.txt"), "w+", newline = "") as f:
	thresholdFile = csv.writer(f, delimiter = "\t")
	thresholdDict = {}
	geneDict = get_gene_dict_with_best_QS()
	for gene in geneDict.keys():
		thresholdDict[gene] = [gene, geneDict[gene]["geneInfo"]["ENSEMBL_ID"], 
								geneDict[gene]["geneInfo"]["exp_val"],
								geneDict[gene]["BS_data"][0]["Chromosome"],
								geneDict[gene]["geneInfo"]["num_of_promoters"]
								]
	header = ["geneName", "EnsemblID", "exp_val", "Chromosome", "QS_num_of_promoters"]
	thresholdList = ["", 370,310,350,450,480,330,480,"",""]
	QSList = list(geneDict["IGF2"]["best_QS"].keys()) 
	QSList[QSList.index("avg10")] = "avg_10_highest_BS"
	header.extend(["QS_" + str(i) + "_" + str(thresholdList[threshold]) for threshold, i in enumerate(QSList)])
	for gene in geneDict.keys():
		for QS in [geneDict[gene]["best_QS"][j] for j in geneDict[gene]["best_QS"].keys()]:
			thresholdDict[gene].append(QS)
	thresholdFile.writerow(header)
	for geneQSList in [thresholdDict[gene] for gene in thresholdDict.keys()]:
		thresholdFile.writerow(geneQSList)

# 
with open(os.path.join(".", "important_source_files", "Bachelor_Thesis_Frederick_Thiemer_dataset", "gene_top_QSs_motif_separated.txt"), "w+", newline = "") as f:
	thresholdFile = csv.writer(f, delimiter = "\t")
	thresholdDict = {}
	geneDict = get_motif_specific_gene_dict_with_best_QS()
	for gene in geneDict.keys():
		thresholdDict[gene] = [gene, geneDict[gene]["geneInfo"]["ENSEMBL_ID"], 
								geneDict[gene]["geneInfo"]["exp_val"],
								geneDict[gene]["BS_data"][0]["Chromosome"],
								geneDict[gene]["geneInfo"]["num_of_promoters"]
								]
	header = ["geneName", "EnsemblID", "exp_val", "Chromosome", "QS_num_of_promoters"]
	thresholdList = ["", 370,310,350,450,480,330,480,"",""]
	QSList = list(geneDict[list(geneDict.keys())[0]]["best_QS"].keys())
	QSList[QSList.index("avg10")] = "avg_10_highest_BS"
	QSList = ["QS_" + str(i) + "_" + str(thresholdList[threshold]) for threshold, i in enumerate(QSList)]
	motifSpecificQSList = [str(i) + "_" + str(x) for i in QSList for x in ["RELA", "RELB", "REL", "all"]]
	for i in ["QS_Wang1_480_RELB", "QS_Wang1_480_REL", "QS_Wang1_480_all", "QS_Wang5__RELB", "QS_Wang5__REL", "QS_Wang5__all"]:
		del motifSpecificQSList[motifSpecificQSList.index(i)]
	header.extend(motifSpecificQSList)
	for gene in geneDict.keys():
		for QS in list(geneDict[gene]["best_QS"].keys()):
			for key in ["RELA", "RELB", "REL", "all"]:
				if str(QS)[:4] == "Wang":
					if key == "RELA":
						thresholdDict[gene].append(geneDict[gene]["best_QS"][QS])
				else:
					thresholdDict[gene].append(geneDict[gene]["best_QS"][QS][key])
	thresholdFile.writerow(header)
	for geneQSList in [thresholdDict[gene] for gene in thresholdDict.keys()]:
		thresholdFile.writerow(geneQSList)
		
# Used for Top Ten BSs column in Table 8
with open(os.path.join(".", "ABC_and_final_motif_specific_QSs.txt"), "w+", newline = "") as f:
    finalFile = csv.writer(f, delimiter = "\t")
    geneDict = get_gene_dict_with_QS()
	header = ["gene", "promoters", "Ensembl_ID", "Chromosome", "verified?", "BSs_RelA", "BSs_RelB", "BSs_c-Rel", "BSs_total", "QS_avg_10_highest_BS(all BSs)","Top_10_BSs_RelA", "Top_10_BSs_RelB", "Top_10_BSs_c-Rel"]
    header.extend(["QS_" + motif[0] + "_avg_" + str(motif[1]) for motif in [("RelA",4), ("RelB", 10), ("c-Rel", 5)]])
    finalFile.writerow(header)
    for gene in geneDict:
        rowToWrite = [gene, 
                      geneDict[gene]["geneInfo"]["num_of_promoters"],
                      geneDict[gene]["geneInfo"]["ENSEMBL_ID"],
                      geneDict[gene]["BS_data"][0]["Chromosome"],
                      geneDict[gene]["geneInfo"]["exp_val"]
                      ]
        for motifType in ["RELA", "RELB", "REL", "all"]:
            rowToWrite.append(geneDict[gene]["QSs"]["BS_sum_all_motifs"][motifType])
        rowToWrite.append(geneDict[gene]["avgHighestXmsp"]["avg10"]["all"])
        abcCountList = [geneDict[gene]["BS_data"][i]["CRE"] for i in range(10)]
        rowToWrite.extend([abcCountList.count(i) for i in ["RELA", "RELB", "REL"]])
        for motifType,number in [("RELA",4), ("RELB", 10), ("REL", 5)]:
            rowToWrite.append(geneDict[gene]["avgHighestXmsp"]["avg" + str(number)][motifType])
        finalFile.writerow(rowToWrite)
import matplotlib.pyplot as plt
import numpy as np
import os
from Bielefeld_dataset_basic_python_dicts import *
from scipy import mean

plt.rc("font", size = 9)
def create_QS_dict(floaty = True):
	QS_Dict = {}
	QS_Dict["QS_0"] = {"Yes": [], "No": []}
	QS_Dict["QS_1"] = {"Yes": [], "No": []}
	QS_Dict["QS_2"] = {"Yes": [], "No": []}
	QS_Dict["QS_3"] = {"Yes": [], "No": []}
	QS_Dict["QS_Wang"] = {"Yes": [], "No": []}
	QS_Dict["QS_avgScore"] = {"Yes": [], "No": []}
	geneDict = get_BB8_gene_BS_sum_dict()
	#calculate expected scores for QS_0
	expScoreDict = {"RELA": [], "RELB": [], "REL": [], "all": []}
	expLenDict = {"RELA": [], "RELB": [], "REL": [], "all": []}
	for gene in geneDict:
		expLenDict["all"].append(int(len(geneDict[gene]["BS_data"])))
		BS_counter = {"RELA":0, "RELB":0, "REL":0}
		for i, BS in enumerate(geneDict[gene]["BS_data"]):
			BS_counter[BS["CRE"]] += 1
			expScoreDict["all"].append(int(BS["Score"]))
			expScoreDict[BS["CRE"]].append(int(BS["Score"]))
		for countCategory in BS_counter:
			expLenDict[countCategory].append(BS_counter[countCategory])
	for key in expScoreDict:
		expScoreDict[key] = mean(expScoreDict[key])
		expLenDict[key] = mean(expLenDict[key])
	#calculate QSs
	for gene in geneDict:
		#create Dicts
		RELA_cW_sum = 0	
		RELA_cW_sur_fit_sum = 0
		sumDict = {"RELA": 0, "RELB": 0, "REL": 0, "all":0}
		avgScoreDict = {"RELA": [], "RELB": [], "REL": [], "all": []}
		avgTSSDict = {"RELA": [], "RELB": [], "REL": [], "all": []}
		QS_0 = {"RELA": [], "RELB": [], "REL": [], "all": []}
		QS_1 = {"RELA": [], "RELB": [], "REL": [], "all": []}
		QS_2 = {"RELA": [], "RELB": [], "REL": [], "all": []}
		QS_3 = {"RELA": [], "RELB": [], "REL": [], "all": []}
		avgHighest15 = {"RELA": [], "RELB": [], "REL": [], "all": []}
		avgHighest10 = {"RELA": [], "RELB": [], "REL": [], "all": []}
		#gather data
		for i, BS in enumerate(geneDict[gene]["BS_data"]):
			sumDict["all"] += 1
			sumDict[BS["CRE"]] += 1
			avgScoreDict["all"].append(int(BS["Score"]))
			avgScoreDict[BS["CRE"]].append(int(BS["Score"]))
			avgTSSDict["all"].append(int(BS["dist_TSS"]))
			avgTSSDict[BS["CRE"]].append(int(BS["dist_TSS"]))
			RELA_cW_sum += int(BS["RELA_central_BS_W"]) # werden in vorherigen Files zusammengefasst, je nach File das get_BB8_gene_BS_sum_dict benutzt
			RELA_cW_sur_fit_sum += int(BS["RELA_cW_sur_fit"]) # 24 nimmt nur RelA binding sites, 25 alle, die die Voraussetzung erfüllt
		if floaty:
			# calculate QSs
			for key in avgScoreDict:
				currentScoreList = sorted(avgScoreDict[key], reverse = True)
				if currentScoreList == []:
					QS_0[key] = 0
					QS_1[key] = 0
					QS_2[key] = 0
					QS_3[key] = 0
					avgHighest15[key] = 0
					avgHighest10[key] = 0
				else:
					QS_0[key] = round(((mean(currentScoreList) * len(currentScoreList))/(expScoreDict[key] * expLenDict[key])), 3)
					QS_1[key] = 0
					QS_2[key] = 0
					for rank, score in enumerate(currentScoreList):
						QS_1[key] += score/(rank + 1)
						QS_2[key] += score**(1/(rank + 1))
					QS_1[key] = round(QS_1[key], 3)
					QS_2[key] = round(QS_2[key], 3)
					QS_3[key] = int(currentScoreList[0])
					avgHighest15[key] = round(mean(currentScoreList[:15]), 3)
					avgHighest10[key] = round(mean(currentScoreList[:10]), 3)
		exp_val = geneDict[gene]["geneInfo"]["exp_val"]
		QS_Dict["QS_0"][exp_val].append(QS_0["all"])
		QS_Dict["QS_1"][exp_val].append(QS_1["all"])
		QS_Dict["QS_2"][exp_val].append(QS_2["all"])
		QS_Dict["QS_3"][exp_val].append(QS_3["all"])
		QS_Dict["QS_Wang"][exp_val].append(RELA_cW_sur_fit_sum)
		QS_Dict["QS_avgScore"][exp_val].append(round(mean(avgScoreDict["all"]), 3))
	return QS_Dict

def bi_vs_all(title, hist_filename, dataset_bi, dataset_others):#datasets as list
	fig, ax = plt.subplots()
	ax.set( xlim = (min(dataset_bi) - min(dataset_bi)/10, max(dataset_bi) + max(dataset_bi)/10), xlabel = QS,
	ylabel = "Relative Frequency", title = title )
	
	plt.hist(np.array(dataset_others), bins = 50, density = True,
	histtype = "stepfilled", alpha = 0.5, color = "#ff7f0e", label ="Verified target genes")
	
	plt.hist(np.array(dataset_bi), bins = 50, density = True,
	histtype = "step", alpha = 1, color = "#1f77b4", label = "Computationally identified\ntarget genes (present study)")
	

	plt.legend()
	#plt.show()
	fig.savefig(os.path.join(".", "Figures",str(hist_filename)+".png"), dpi = 600)

QS_Dict = create_QS_dict()
for QS in QS_Dict:
	bi_vs_all("Comparison based on " + str(QS), str("Figure_5_" + str(QS) + "_comparison_bi_vs_exp_val"), QS_Dict[QS]["Yes"], QS_Dict[QS]["No"])
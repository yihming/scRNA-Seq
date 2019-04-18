#!/usr/bin/env python

import os
from sys import argv, exit

import numpy as np
import pandas as pd
import xlsxwriter
import scCloud

def generate_de_results(input_file, attr_label, attr_cond, cond_order, output_file):
	adata = scCloud.tools.read_input(input_file, mode = 'a')
	results = scCloud.misc.mwu_test(adata, attr_label, attr_cond, cond_order, n_jobs = 7)

	cell_types = []
	num_de_genes = []

	writer = pd.ExcelWriter(output_file, engine='xlsxwriter')
	for label in adata.obs[attr_label].cat.categories:
		df = results[label]
		passed = (df['mwu_qval'] <= 0.05).values

		cell_types.append(label)
		num_de_genes.append(passed.sum())

		idx_up = (df['log_fc'] > 0).values
		idx_down = (df['log_fc'] < 0).values
		df_res = df.loc[passed & idx_up].copy()
		if df_res.shape[0] > 0:
			df_res.sort_values(by = "mwu_qval", ascending = True, inplace = True)
		df_res.to_excel(writer, sheet_name = '{}, up'.format(label))
		df_res = df.loc[passed & idx_down].copy()
		if df_res.shape[0] > 0:
			df_res.sort_values(by = "mwu_qval", ascending = True, inplace = True)
		df_res.to_excel(writer, sheet_name = '{}, down'.format(label))
	writer.save()

	return cell_types, num_de_genes


if __name__ == "__main__":
	if len(argv) != 2:
		print("python generate_de_results.py output_directory")
		exit(-1)

	cell_types, num_de_genes = generate_de_results(argv[1] + 'experiment1_stonly_lowremove.h5ad', 'coarse', 'Sample', ['human_st', 'human_control'], argv[1] + 'table_s4.de.xlsx')
	if os.system("scCloud de_analysis --labels coarse -p 7 --mwu --subset Sample:human_control " + argv[1] + "experiment1_stonly_lowremove.h5ad " + argv[1] + "experiment1_stonly_lowremove_control.de.xlsx"):
		exit(-1)

	adata = scCloud.tools.read_input(argv[1] + "experiment1_stonly_lowremove_control.de.h5ad")
	num_control_de = []
	for cell_type in cell_types:
		num_control_de.append((adata.var["mwu_qval_" + cell_type] <= 0.05).values.sum())
	percentages = []
	for i in range(len(cell_types)):
		percentages.append("{:.2%}".format(num_de_genes[i] * 1.0 / num_control_de[i]))

	df = pd.DataFrame({1: num_de_genes, 2: num_control_de, 3: percentages}, index = cell_types)
	df.index.name = "Cell type"
	df.columns = ['Number of DE genes between hashed and control', 'Number of DE genes between control nuclei of the same cell type  and all other control nuclei', 'Percentage of column 2 over column 3']
	df.to_csv(argv[1] + "table_s4.txt", sep = '\t')

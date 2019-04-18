import scCloud
import os, re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

input_dir = '/software/inputs/'
output_dir = '/software/outputs/'

def compare_results(file_list):
	# Load ground truth.
	ground_res = scCloud.tools.read_input('{output}experiment1_human_st_demux.h5ad'.format(output = output_dir), mode = 'a')
	# Load names.
	names = np.loadtxt('{output}experiment1_stonly_barcodes.txt'.format(output = output_dir), dtype = 'str', delimiter = '\n')
	# Iterate over experiment output files.
	corr_arr = []
	for file_name in file_list:
		exp_res = scCloud.tools.read_input(output_dir + file_name, mode = 'a')
		n_total = exp_res.obs.loc[names, 'assignment'].size
		percent = (ground_res.obs.loc[names, 'assignment'].astype('str') == exp_res.obs.loc[names, 'assignment'].astype('str')).sum() * 1.0 / n_total
		corr_arr.append(percent)
	return corr_arr

def plot_result(df_list):
	df = pd.concat(df_list).reset_index()
	#sns.set(style = "whitegrid")
	ax = sns.barplot(x = 'experiment', y = 'consistency rate', data = df)
	ax.set(ylabel = 'Consistency', xlabel = '')
	ax.set_ylim(0.9, 1.01)

	vals = ax.get_yticks()
	ax.set_yticklabels(['{0:.0%}'.format(x) for x in vals])
	
	ax.get_figure().savefig("{output}figure_s3a.pdf".format(output = output_dir))
	plt.close()

def generate_consistency_dataframe(experiment_name, consist_arr):
	n_obs = len(consist_arr)
	df = pd.DataFrame({'experiment': [experiment_name] * n_obs, 'consistency rate': consist_arr})
	return df
        
def random_state_experiment(low, high, size):
	# Generate random seed array.
	file_list = []

	seed_arr = np.array([], dtype = int)
	for i in range(size):
		num = np.random.randint(low, high, 1)[0]
		while (seed_arr == num).sum() > 0:
			num = np.random.randint(low, high, 1)[0]
		seed_arr = np.append(seed_arr, num)

	# Execute with random states.
	for seed in seed_arr:
		os.system("scCloud demuxEM -p 8 --hash-type nuclei-hashing --min-num-genes 200 --generate-diagnostic-plots --random-state {state} {input}experiment1_human_st_ADT.csv {input}experiment1_human_st_raw_10x.h5 {output}random{state}".format(state = seed, input = input_dir, output = output_dir))

	# Compare Results with original analysis.
	file_list = ["random{seed}_demux.h5ad".format(seed = seed) for seed in seed_arr]


	consist_arr = compare_results(file_list)

	df = generate_consistency_dataframe('k-means', consist_arr)
	return df


def hierarchical_clustering_experiment():
	# Four different linkage methods
	link_list = ['ward', 'complete', 'average', 'single']

	# Execute with linkage methods.
	for link in link_list:
		os.system("scCloud demuxEM -p 8 --hash-type nuclei-hashing --min-num-genes 200 --hier-linkage {linkage} --generate-diagnostic-plots {input}experiment1_human_st_ADT.csv {input}experiment1_human_st_raw_10x.h5 {output}linkage_{linkage}".format(linkage = link, input = input_dir, output = output_dir))


	# Compare Results with original analysis.
	file_list = ["linkage_{linkage}_demux.h5ad".format(linkage = link) for link in link_list]
	consist_arr = compare_results(file_list)

	df = generate_consistency_dataframe('Hierarchical clustering', consist_arr)
	return df

"""
For threshold in [T_min, T_max].
"""
def hard_threshold_experiment(T_min, T_max, step):
	
	# Execute with thresholds.
	T_arr = np.arange(T_min, T_max + step, step)
	for t in T_arr:
		os.system("scCloud demuxEM -p 8 --hash-type nuclei-hashing --min-num-genes 200 --hard-threshold {threshold} --generate-diagnostic-plots {input}experiment1_human_st_ADT.csv {input}experiment1_human_st_raw_10x.h5 {output}threshold_{threshold}".format(threshold = t, input = input_dir, output = output_dir))
		        

	# Compare Results with orginal analysis.
	file_list = ["threshold_{thresh}_demux.h5ad".format(thresh = t) for t in T_arr]
	consist_arr = compare_results(file_list)

	# Print Consistency Rates for numeric summary.
	df_summary = pd.DataFrame({'threshold': T_arr, 'consistency rate': consist_arr}).sort_values(by = ['threshold'])

	df = generate_consistency_dataframe('Threshold', consist_arr)
	return df

def prior_on_samples_experiment(P_min, P_max, step):
	
	# Execute with priors.
	P_arr = np.arange(P_min, P_max + step, step)
	for p in P_arr:
		os.system("scCloud demuxEM -p 8 --hash-type nuclei-hashing --min-num-genes 200 --prior-on-samples {prior} --generate-diagnostic-plots {input}experiment1_human_st_ADT.csv {input}experiment1_human_st_raw_10x.h5 {output}prior_{prior}".format(prior = p, input = input_dir, output = output_dir))

	# Compare Results with original analysis.
	file_list = ["prior_{prior}_demux.h5ad".format(prior = p) for p in P_arr]
	consist_arr = compare_results(file_list)

	# Print Consistency Rates for numeric summary.
	df_summary = pd.DataFrame({'prior': P_arr, 'consistency rate': consist_arr}).sort_values(by = ['prior'])
	
	df = pd.DataFrame({'prior on samples': P_arr, 'consistency rate': consist_arr}).sort_values(by = ['prior on samples'])
	ax = sns.lineplot(x = 'prior on samples', y = 'consistency rate', data = df)
	ax = sns.scatterplot(x = 'prior on samples', y = 'consistency rate', data = df)
	ax.set(ylim = (0.79, 1.01), xlabel = 'Prior on Samples', ylabel = 'Consistency')


	vals = ax.get_yticks()
	ax.set_yticklabels(['{0:.0%}'.format(x) for x in vals])

	ax.get_figure().savefig("{output}figure_s3b.pdf".format(output = output_dir))
	plt.close()

def main():
	df_list = []

	df_random_state = random_state_experiment(0, 2**32, 10)
	df_list.append(df_random_state)

	df_hierarchical = hierarchical_clustering_experiment()
	df_list.append(df_hierarchical)

	df_threshold = hard_threshold_experiment(10, 100, 10)
	df_list.append(df_threshold)

	plot_result(df_list)

	df_prior = prior_on_samples_experiment(0.1, 1.0, 0.1)
        

if __name__ == '__main__':
	main()

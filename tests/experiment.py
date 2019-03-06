import scCloud
import os, re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def compare_results(file_list):
	# Load ground truth.
	ground_res = scCloud.tools.read_input('pilot4_human_st_demux.h5ad', mode = 'a')
	# Load names.
	names = np.loadtxt('pilot4_st_valid.txt', dtype = 'str', delimiter = '\n')
	# Iterate over experiment output files.
	corr_arr = []
	for file_name in file_list:
		exp_res = scCloud.tools.read_input(file_name, mode = 'a')
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
	
	ax.get_figure().savefig("summary_consistency.png")
	plt.close()

def generate_consistency_dataframe(experiment_name, consist_arr):
	n_obs = len(consist_arr)
	df = pd.DataFrame({'experiment': [experiment_name] * n_obs, 'consistency rate': consist_arr})
	return df
        
def random_state_experiment(low, high, size, start_over = True):
	# Generate random seed array.
	file_list = []

	if start_over:
		seed_arr = np.array([], dtype = int)
		for i in range(size):
			num = np.random.randint(low, high, 1)[0]
			while (seed_arr == num).sum() > 0:
				num = np.random.randint(low, high, 1)[0]
			seed_arr = np.append(seed_arr, num)

		# Execute with random states.
		for seed in seed_arr:
			os.system("scCloud demuxEM -p 8 --hash-type nuclei-hashing --min-num-genes 200 --generate-diagnostic-plots --random-state {state} pilot4_human_st_ADT.csv pilot4_human_st_raw_10x.h5 random{state}".format(state = seed))

		# Compare Results with original analysis.
		file_list = ["random{seed}_demux.h5ad".format(seed = seed) for seed in seed_arr]
	else:
		file_list = [f for f in os.listdir('.') if re.match(r'random[0-9]+_demux\.h5ad', f)]


	consist_arr = compare_results(file_list)
	print(consist_arr)

	df = generate_consistency_dataframe('KMeans', consist_arr)
	return df


def hierarchical_clustering_experiment(start_over = True):
	# Four different linkage methods
	link_list = ['ward', 'complete', 'average', 'single']

	# Execute with linkage methods.
	if start_over:
		for link in link_list:
			os.system("scCloud demuxEM -p 8 --hash-type nuclei-hashing --min-num-genes 200 --hier-linkage {linkage} --generate-diagnostic-plots pilot4_human_st_ADT.csv pilot4_human_st_raw_10x.h5 linkage_{linkage}".format(linkage = link))

		# Plot histograms of separated droplets under different linkage methods.
		for link in link_list:
			adt = scCloud.tools.read_input("linkage_{link}_ADTs.h5ad".format(link = link), mode = 'a')
			scCloud.demuxEM.plot_adt_hist(adt, 'hto_type', "{link}_hist".format(link = link))

	# Compare Results with original analysis.
	file_list = ["linkage_{linkage}_demux.h5ad".format(linkage = link) for link in link_list]
	consist_arr = compare_results(file_list)

	# Print consistency rates for numeric summary.
	dict_summary = dict([tuple(elem) for elem in np.vstack([link_list, consist_arr]).T])
	print(dict_summary)

	df = generate_consistency_dataframe('Hierarchical clustering', consist_arr)
	return df

"""
For threshold in [T_min, T_max].
"""
def hard_threshold_experiment(T_min, T_max, step, start_over = True):
	T_arr = []
	# Execute with thresholds.
	if start_over:
		T_arr = np.arange(T_min, T_max + step, step)
		for t in T_arr:
			os.system("scCloud demuxEM -p 8 --hash-type nuclei-hashing --min-num-genes 200 --hard-threshold {threshold} --generate-diagnostic-plots pilot4_human_st_ADT.csv pilot4_human_st_raw_10x.h5 threshold_{threshold}".format(threshold = t))
		        
		# Plot histograms of separated droplets under different thresholds.
		for t in T_arr:
			adt = scCloud.tools.read_input("threshold_{threshold}_ADTs.h5ad".format(threshold = t), mode = 'a')
			scCloud.demuxEM.plot_adt_hist(adt, 'hto_type', "threshold_{threshold}_hist".format(threshold = t))

	# Compare Results with orginal analysis.
	file_list = [f for f in os.listdir('.') if re.match(r'threshold_[0-9]+_demux\.h5ad', f)]
	if not start_over:
		T_arr = np.array([float(re.findall(r'\d+', f)[0]) for f in file_list])
	consist_arr = compare_results(file_list)

	# Print Consistency Rates for numeric summary.
	df_summary = pd.DataFrame({'threshold': T_arr, 'consistency rate': consist_arr}).sort_values(by = ['threshold'])
	print(df_summary['consistency rate'].values.tolist())

	df = generate_consistency_dataframe('Threshold', consist_arr)
	return df

def prior_on_samples_experiment(P_min, P_max, step, start_over = True):
	P_arr = []
	# Execute with priors.
	if start_over:
		P_arr = np.arange(P_min, P_max + step, step)
		for p in P_arr:
			os.system("scCloud demuxEM -p 8 --hash-type nuclei-hashing --min-num-genes 200 --prior-on-samples {prior} --generate-diagnostic-plots pilot4_human_st_ADT.csv pilot4_human_st_raw_10x.h5 prior_{prior}".format(prior = p))

		# Plot histograms of separated droplets under different priors.
		for p in P_arr:
			adt = scCloud.tools.read_input("prior_{prior}_ADTs.h5ad".format(prior = p), mode = 'a')
			scCloud.demuxEM.plot_adt_hist(adt, 'hto_type', "prior_{prior}_hist".format(prior = p))

	# Compare Results with original analysis.
	file_list = [f for f in os.listdir('.') if re.match(r'prior_[0-9]+(\.[0-9]+)?_demux\.h5ad', f)]
	if not start_over:
		P_arr = np.array([float(re.findall(r'[0-9]+[\.[0-9]+]?', f)[0]) for f in file_list])

	consist_arr = compare_results(file_list)

	# Print Consistency Rates for numeric summary.
	df_summary = pd.DataFrame({'prior': P_arr, 'consistency rate': consist_arr}).sort_values(by = ['prior'])
	print(df_summary['consistency rate'].values.tolist())

	df = pd.DataFrame({'prior on samples': P_arr, 'consistency rate': consist_arr}).sort_values(by = ['prior on samples'])
	ax = sns.lineplot(x = 'prior on samples', y = 'consistency rate', data = df)
	ax = sns.scatterplot(x = 'prior on samples', y = 'consistency rate', data = df)
	ax.set(ylim = (0.79, 1.01), xlabel = 'Prior on Samples', ylabel = 'Consistency')


	vals = ax.get_yticks()
	ax.set_yticklabels(['{0:.0%}'.format(x) for x in vals])

	ax.get_figure().savefig("prior_on_samples_plot.pdf")
	plt.close()

def main():
	df_list = []

	df_random_state = random_state_experiment(0, 2**32, 10, start_over = True)
	df_list.append(df_random_state)

	df_hierarchical = hierarchical_clustering_experiment(start_over = True)
	df_list.append(df_hierarchical)

	df_threshold = hard_threshold_experiment(10, 100, 10, start_over = True)
	df_list.append(df_threshold)

	plot_result(df_list)

	df_prior = prior_on_samples_experiment(0.1, 1.0, 0.1, start_over = True)
        

if __name__ == '__main__':
	main()

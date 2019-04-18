import time
import re
import numpy as np 
import pandas as pd
import anndata
import xlsxwriter
from collections import Counter

from scipy.sparse import issparse
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.utils.sparsefuncs import mean_variance_axis
from sklearn.utils.extmath import randomized_svd



transfer_gene_name = [(358, 'ENSG00000268991', 'FAM231C.2'), (921, 'ENSG00000278139', 'AL358075.4'), (2207, 'ENSG00000232995', 'RGS5.2'), (5847, 'ENSG00000282827', 'AC134772.2'), (5938, 'ENSG00000271858', 'CYB561D2.2'), (6087, 'ENSG00000241572', 'PRICKLE2-AS1.2'), (7213, 'ENSG00000249428', 'CFAP99.2'), (9596, 'ENSG00000280987', 'MATR3.2'), (9605, 'ENSG00000279686', 'AC142391.1'), (10277, 'ENSG00000282913', 'BLOC1S5.2'), (10867, 'ENSG00000124593', 'AL365205.1'), (11619, 'ENSG00000268592', 'RAET1E-AS1.2'), (13877, 'ENSG00000231963', 'AL662864.1'), (16117, 'ENSG00000225655', 'BX255923.1'), (16938, 'ENSG00000282955', 'RABL6.2'), (17241, 'ENSG00000265264', 'TIMM10B.2'), (18626, 'ENSG00000282682', 'C11orf71.2'), (18984, 'ENSG00000282883', 'AKR1C3.2'), (19226, 'ENSG00000150076', 'CCDC7.2'), (19346, 'ENSG00000264404', 'BX547991.1'), (21184, 'ENSG00000282031', 'TMBIM4.2'), (21230, 'ENSG00000257815', 'LINC01481.2'), (22033, 'ENSG00000228741', 'SPATA13.2'), (22037, 'ENSG00000281899', 'AL359736.3'), (22654, 'ENSG00000274827', 'LINC01297.2'), (23662, 'ENSG00000273259', 'AL049839.2'), (24019, 'ENSG00000211974', 'AC245369.1'), (26919, 'ENSG00000279257', 'C17orf100.2'), (26962, 'ENSG00000187838', 'PLSCR3'), (27137, 'ENSG00000255104', 'AC005324.4'), (27884, 'ENSG00000263715', 'LINC02210-CRHR1'), (28407, 'ENSG00000281844', 'FBF1.2'), (30440, 'ENSG00000283027', 'CAPS.2'), (32648, 'ENSG00000235271', 'LINC01422.2')]

def update_var_names(data):
	start = time.time()

	prefix = re.compile('^(' + '|'.join(data.uns['genome'].split(',')) + ')_+')
	if prefix.match(data.var_names[0]):
		data.var['gene_ids'] = [prefix.sub('', x) for x in data.var['gene_ids']]
		data.var_names = pd.Index([prefix.sub('', x) for x in data.var_names])

	gsyms = data.var_names.values
	
	if data.uns['genome'] == 'GRCh38':
		for pos, gid, gsym in transfer_gene_name:
			assert data.var.iloc[pos, 0] == gid
			gsyms[pos] = gsym
	else:	
		dup_ids = Counter()
		for i in range(gsyms.size):
			idn = dup_ids[gsyms[i]]
			dup_ids[gsyms[i]] += 1
			if idn > 0:
				gsyms[i] = gsyms[i] + ".{}".format(idn)
	
	data.var_names = pd.Index(gsyms)

	end = time.time()
	print("update_var_names is finished. Time spent = {:.2f}s.".format(end - start))



def filter_data(data, output_filt = None, plot_filt = None, plot_filt_figsize = None, mito_prefix = 'MT-', min_genes = 500, max_genes = 6000, min_umis = 100, max_umis = 600000, percent_mito = 0.1, percent_cells = 0.0005, min_genes_on_raw = 100):
	start = time.time()

	data.obs['n_genes'] = data.X.getnnz(axis = 1)
	if data.obs['n_genes'].min() == 0: # if raw.h5
		data._inplace_subset_obs(data.obs['n_genes'].values >= min_genes_on_raw)
	data.obs['n_counts'] = data.X.sum(axis = 1).A1

	mito_prefixes = mito_prefix.split(',')

	def startswith(name):
		for prefix in mito_prefixes:
			if name.startswith(prefix):
				return True
		return False

	mito_genes = data.var_names.map(startswith).values.nonzero()[0]
	data.obs['percent_mito'] = data.X[:, mito_genes].sum(axis=1).A1 / np.maximum(data.obs['n_counts'].values, 1.0)

	if output_filt is not None:
		writer = pd.ExcelWriter(output_filt + '.filt.xlsx', engine='xlsxwriter')
		gb1 = data.obs.groupby('Channel')
		df_before = gb1.median()
		df_before = df_before.assign(total = gb1.size())
		df_before.rename(columns = {'n_genes' : 'median_n_genes_before', 'n_counts' : 'median_n_umis_before', 'percent_mito' : 'median_percent_mito_before'}, inplace = True)

	if plot_filt is not None:
		df_plot_before = data.obs[['Channel', 'n_genes', 'n_counts', 'percent_mito']].copy()
		df_plot_before.reset_index(drop = True, inplace = True)
		df_plot_before['status'] = 'original'

	# Filter cells		
	obs_index = np.logical_and.reduce((data.obs['n_genes'] >= min_genes, 
									   data.obs['n_genes'] < max_genes,
									   data.obs['n_counts'] >= min_umis,
									   data.obs['n_counts'] < max_umis,
									   data.obs['percent_mito'] < percent_mito))
	data._inplace_subset_obs(obs_index)

	if output_filt is not None:
		gb2 = data.obs.groupby('Channel')
		df_after = gb2.median()
		df_after = df_after.assign(kept = gb2.size())
		df_after.rename(columns = {'n_genes' : 'median_n_genes', 'n_counts' : 'median_n_umis', 'percent_mito' : 'median_percent_mito'}, inplace = True)
		df = pd.concat((df_before, df_after), axis = 1, sort = False)
		df.fillna(0, inplace = True)
		df['kept'] = df['kept'].astype(int)
		df['filt'] = df['total'] - df['kept']
		df = df[['kept', 'median_n_genes', 'median_n_umis', 'median_percent_mito', 'filt', 'total', 'median_n_genes_before', 'median_n_umis_before', 'median_percent_mito_before']]
		df.sort_values('kept', inplace = True)
		df.to_excel(writer, sheet_name = "Cell filtration stats")

	if plot_filt is not None:
		df_plot_after = data.obs[['Channel', 'n_genes', 'n_counts', 'percent_mito']].copy()
		df_plot_after.reset_index(drop = True, inplace = True)
		df_plot_after['status'] = 'filtered'
		df_plot = pd.concat((df_plot_before, df_plot_after), axis = 0)
		from scCloud.plotting import plot_qc_violin
		figsize = None
		if plot_filt_figsize is not None:
			width, height = plot_filt_figsize.split(',')
			figsize = (int(width), int(height))
		plot_qc_violin(df_plot, 'count', plot_filt + '.filt.UMI.pdf', xattr = 'Channel', hue = 'status', xlabel = 'Channel', split = True, linewidth = 0, figsize = figsize)
		plot_qc_violin(df_plot, 'gene', plot_filt + '.filt.gene.pdf', xattr = 'Channel', hue = 'status', xlabel = 'Channel', split = True, linewidth = 0, figsize = figsize)
		plot_qc_violin(df_plot, 'mito', plot_filt + '.filt.mito.pdf', xattr = 'Channel', hue = 'status', xlabel = 'Channel', split = True, linewidth = 0, figsize = figsize)
		print("Filtration plots are generated.")

	# Filter genes
	data.var['n_cells'] = data.X.getnnz(axis = 0)
	data.var['percent_cells'] = data.var['n_cells'] / data.shape[0]
	data.var['robust'] = data.var['percent_cells'] >= percent_cells

	if output_filt is not None:
		idx = data.var['robust'] == False
		df = pd.DataFrame({'n_cells': data.var.loc[idx, 'n_cells'], 'percent_cells': data.var.loc[idx, 'percent_cells']})
		df.index.name = 'gene'
		df.sort_values('n_cells', ascending = False, inplace = True)
		df.to_excel(writer, sheet_name = "Gene filtration stats")
		writer.save()
		print("Filtration results are written.")

	var_index = (data.var['n_cells'] > 0).values
	data._inplace_subset_var(var_index)
	print("After filteration, {nc} cells and {ng} genes are kept. Among {ng} genes, {nrb} genes are robust.".format(nc = data.shape[0], ng = data.shape[1], nrb = data.var['robust'].sum()))

	end = time.time()
	print("filter_data is finished. Time spent = {:.2f}s.".format(end - start))


def filter_cells_cite_seq(data, max_cells):
	data.obs['n_counts'] = data.X.sum(axis = 1).A1
	obs_index = np.zeros(data.shape[0], dtype = bool)
	obs_index[np.argsort(data.obs['n_counts'].values)[::-1][:max_cells]] = True
	data._inplace_subset_obs(obs_index)
	data.var['robust'] = True
	print("After filteration, {nc} cells are kept, with the minimum nUMI = {numi}.".format(nc = max_cells, numi = data.obs['n_counts'].min()))

def log_norm(data, norm_count):
	""" Normalization and then take log """
	start = time.time()

	assert issparse(data.X)
	mat = data.X[:, data.var['robust'].values]
	scale = norm_count / mat.sum(axis = 1).A1
	data.X.data *= np.repeat(scale, np.diff(data.X.indptr))
	data.X = data.X.log1p()

	end = time.time()
	print("Normalization is finished. Time spent = {:.2f}s.".format(end - start))

def run_pca(data, standardize = True, max_value = 10, nPC = 50, random_state = 0):
	start = time.time()
	if issparse(data.X):
		data.X = data.X.toarray()

	if standardize:
		scaler = StandardScaler(copy = False)
		scaler.fit_transform(data.X)

	if max_value is not None:
		data.X[data.X > max_value] = max_value

	pca = PCA(n_components = nPC, random_state = random_state)
	X_pca = pca.fit_transform(data.X)	
	data.obsm['X_pca'] = X_pca
	data.varm['PCs'] = pca.components_.T
	data.uns['pca'] = {}
	data.uns['pca']['variance'] = pca.explained_variance_
	data.uns['pca']['variance_ratio'] = pca.explained_variance_ratio_
	end = time.time()
	print("PCA is done. Time spent = {:.2f}s.".format(end - start))

def run_rpca(data, scale = False, max_value = 10.0, nPC = 50, random_state = 0):
	""" smooth outliers, then no center/scale data """
	start = time.time()

	# Smooth out outliers
	means, variances = mean_variance_axis(data.X, axis = 0)
	stds = np.sqrt(variances * (data.X.shape[0] / (data.X.shape[0] - 1))) # make it unbiased
	assert (stds == 0.0).sum() == 0

	data_new = (data.X.data - means[data.X.indices]) / stds[data.X.indices]
	outliers = data_new > max_value
	data.X.data[outliers] = max_value * stds[data.X.indices[outliers]] + means[data.X.indices[outliers]]

	if scale:
		data.X.data /= stds[data.X.indices]

	U, S, VT = randomized_svd(data.X, n_components = nPC, random_state = random_state)
	data.obsm['X_rpca'] = U * S

	end = time.time()
	print("RPCA is done. Time spent = {:.2f}s.".format(end - start))

def parse_subset_selections(subset_selections):	
	subsets_dict = {}
	for subset_str in subset_selections:
		print(subset_str)
		attr, value_str = subset_str.split(':')
		if attr in subsets_dict:
			subsets_dict[attr].extend(value_str.split(','))
		else:
			subsets_dict[attr] = value_str.split(',')
	return subsets_dict

def get_anndata_for_subclustering(data, subset_selections):	
	obs_index = np.full_like(data.obs[data.obs.columns.values[0]], True)
	subsets_dict = parse_subset_selections(subset_selections)
	for key, value in subsets_dict.items():
		print(key, 'corresponds to', subsets_dict[key])
		obs_index = obs_index & np.isin(data.obs[key], value)	
	data = data[obs_index, :]
	obs_dict = {"obs_names" : data.obs_names.values}
	for attr in data.obs.columns:
		if attr != "pseudotime":
			if attr.find("_labels") < 0:
				obs_dict[attr] = data.obs[attr].values
			else:
				obs_dict['parent_' + attr] = data.obs[attr].values

	var_dict = {"var_names" : data.var_names.values, "gene_ids": data.var["gene_ids"].values, "robust": data.var["robust"]}

	newdata = anndata.AnnData(X = data.X, obs = obs_dict, var = var_dict)
	if "Channels" in data.uns:
		newdata.uns["Channels"] = data.uns["Channels"]
	if "Groups" in data.uns:
		newdata.uns["Groups"] = data.uns["Groups"]
	if "means" in data.varm.keys():
		newdata.varm["means"] = data.varm["means"]
	if "stds" in data.varm.keys():
		newdata.varm["stds"] = data.varm["stds"]

	print("{0} cells are selected.".format(newdata.shape[0]))
	
	return newdata

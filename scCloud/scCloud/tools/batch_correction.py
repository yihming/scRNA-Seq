import time
import numpy as np
import numpy.ma as ma
import pandas as pd
from scipy.sparse import issparse
from scipy.stats.mstats import gmean
from collections import defaultdict

def set_group_attribute(data, attribute_string):
	if attribute_string is None:
		data.obs['Group'] = 'one_group'
	elif attribute_string.find('=') >= 0:
		attr, value_str = attribute_string.split('=')
		assert attr in data.obs.columns
		values = value_str.split(';')
		data.obs['Group'] = '0'
		for group_id, value in enumerate(values):
			vals = value.split(',')
			idx = np.isin(data.obs[attr], vals)
			data.obs.loc[idx, 'Group'] = str(group_id + 1)
	elif attribute_string.find('+') >= 0:
		attrs = attribute_string.split('+')
		assert np.isin(attrs, data.obs.columns).sum() == len(attrs)
		data.obs['Group'] = data.obs[attrs].apply(lambda x: '+'.join(x), axis = 1)
	else:
		assert attribute_string in data.obs.columns
		data.obs['Group'] = data.obs[attribute_string]

def estimate_adjustment_matrices(data):
	start = time.time()

	channels = data.obs['Channel'].unique()
	if channels.size <= 1:
		print("Warning: data only contains 1 channel. Batch correction disabled!")
		return None

	means = np.zeros((data.shape[1], channels.size))
	stds = np.zeros((data.shape[1], channels.size))

	group_dict = defaultdict(list)
	is_sparse = issparse(data.X)

	for i, channel in enumerate(channels):
		idx = np.isin(data.obs['Channel'], channel)
		mat = data.X[idx]

		if is_sparse:
			if mat.shape[0] == 1:
				means[:, i] = mat.toarray()[0]
				stds[:, i] = np.zeros(mat.shape[1])
			else:
				means[:, i] = mat.mean(axis = 0).A1
				m2 = mat.power(2).sum(axis = 0).A1
				stds[:, i] = ((m2 - mat.shape[0] * (means[:, i] ** 2)) / (mat.shape[0] - 1.0)) ** 0.5
		else:
			means[:, i] = np.mean(mat, axis = 0)
			stds[:, i] = np.std(mat, axis = 0)
			
		group = data.obs['Group'][idx.nonzero()[0][0]]
		group_dict[group].append(i)

	stds[stds < 1e-6] = 1e-12 # avoid gmean warning	

	groups = data.obs['Group'].unique()
	for group in groups:
		gm = np.mean(means[:, group_dict[group]], axis = 1)
		gs = ma.getdata(gmean(ma.masked_less(stds[:, group_dict[group]], 1e-6), axis = 1))

		for i in group_dict[group]:
			outliers = stds[:, i] < 1e-6
			normals = np.logical_not(outliers)
			stds[outliers, i] = 1.0
			stds[normals, i] = gs[normals] / stds[normals, i]
			means[:, i] = gm - stds[:, i] * means[:, i]

	means[abs(means) < 1e-6] = 0.0
	stds[abs(stds - 1.0) < 1e-6] = 1.0

	data.uns['Channels'] = channels
	data.uns['Groups'] = groups
	data.varm['means'] = means
	data.varm['stds'] = stds

	end = time.time()
	print("batch_correction.estimate_adjustment_matrices is done. Time spent = {:.2f}s.".format(end - start))

def filter_genes_dispersion(data, consider_batch, min_disp=0.5, max_disp=None, min_mean=0.0125, max_mean=7):
	start = time.time()

	X = data.X[:,data.var['robust'].values].expm1()

	mean = X.mean(axis = 0).A1
	var = np.zeros(X.shape[1])

	if consider_batch and 'Channels' not in data.uns:
		print("Warning: Batch correction parameters are not calculated. Switch to not considering batch for variable gene selection.")
		consider_batch = False
	
	if consider_batch:
		for channel in data.uns['Channels']:
			mat = X[np.isin(data.obs['Channel'], channel)]
			if mat.shape[0] == 0:
				continue
			m1 = mat.mean(axis = 0).A1
			m2 = mat.power(2).sum(axis = 0).A1
			var += m2 - mat.shape[0] * (m1 ** 2)

		if data.uns['Groups'].size > 1:
			for group in data.uns['Groups']:
				mat = X[np.isin(data.obs['Group'], group)]
				if mat.shape[0] == 0:
					continue
				var += mat.shape[0] * (mat.mean(axis = 0).A1 - mean) ** 2
		var /= X.shape[0] - 1 
	else:
		m2 = X.power(2).sum(axis = 0).A1
		var = (m2 - X.shape[0] * (mean ** 2)) / (X.shape[0] - 1)

	# now actually compute the dispersion
	mean[mean == 0] = 1e-12  # set entries equal to zero to small value
	dispersion = var / mean

	dispersion[dispersion == 0] = np.nan
	dispersion = np.log(dispersion)
	mean = np.log1p(mean)
	# all of the following quantities are "per-gene" here

	df = pd.DataFrame({'mean' : mean, 'dispersion' : dispersion, 'mean_bin' : pd.cut(mean, bins=20)})
	disp_grouped = df.groupby('mean_bin')['dispersion']
	disp_mean_bin = disp_grouped.mean()
	disp_std_bin = disp_grouped.std(ddof=1)
	df['dispersion_norm'] = (df['dispersion'].values  # use values here as index differs
							 - disp_mean_bin.loc[df['mean_bin']].values) \
							 / disp_std_bin.loc[df['mean_bin']].values
	dispersion_norm = df['dispersion_norm'].values.astype('float32')
	max_disp = np.inf if max_disp is None else max_disp
	dispersion_norm[np.isnan(dispersion_norm)] = 0  # similar to Seurat
	gene_subset = np.logical_and.reduce((mean > min_mean, mean < max_mean,
										 dispersion_norm > min_disp,
										 dispersion_norm < max_disp))

	end = time.time()
	print("batch_correction.filter_genes_dispersion done. Time spent = {:.2f}s.".format(end - start))
	print("{0} genes are selected.".format(gene_subset.sum()))

	return np.rec.fromarrays((gene_subset,
							  df['mean'].values,
							  df['dispersion'].values,
							  df['dispersion_norm'].values.astype('float32', copy=False)),
							  dtype=[('gene_subset', bool),
									 ('means', 'float32'),
									 ('dispersions', 'float32'),
									 ('dispersions_norm', 'float32')])

def collect_variable_gene_matrix(data, gene_subset):
	variable_gene_index = np.zeros(data.shape[1], dtype = bool)
	variable_gene_index[data.var['robust'].values] = gene_subset
	data.var['selected'] = variable_gene_index
	data_c = data[:, variable_gene_index].copy()
	return data_c

def correct_batch_effects(data):
	start = time.time()

	if 'Channels' not in data.uns:
		print("Warning: Batch correction parameters are not calculated. Batch correction skipped!")
		return None

	if issparse(data.X):
		for i, channel in enumerate(data.uns['Channels']):
			idx = np.isin(data.obs['Channel'], channel)
			if idx.sum() == 0:
				continue
			idx = np.repeat(idx, np.diff(data.X.indptr))
			data.X.data[idx] = data.X.data[idx] * data.varm['stds'][:,i][data.X.indices[idx]] + data.varm['means'][:, i][data.X.indices[idx]]
		data.X.data[data.X.data < 0.0] = 0.0
	else:		
		m = data.shape[1]
		for i, channel in enumerate(data.uns['Channels']):
			idx = np.isin(data.obs['Channel'], channel)
			if idx.sum() == 0:
				continue
			data.X[idx] = data.X[idx] * np.reshape(data.varm['stds'][:, i], newshape = (1, m)) + np.reshape(data.varm['means'][:, i], newshape = (1, m))
		data.X[data.X < 0.0] = 0.0

	end = time.time()
	print("batch_correction.correct_batch_effects done. Time spent = {:.2f}s".format(end - start))

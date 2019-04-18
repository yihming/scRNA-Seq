import scCloud
import os, sys
import tables

from shutil import rmtree
from collections import Counter
from scipy.sparse import csr_matrix
from termcolor import cprint
from scCloud.demuxEM.down_sampling import plot_down_sampling

import numpy as np
import pandas as pd
import scCloud.demuxEM.plot as sdp

import s3_experiment as s3exp

input_dir = '/software/inputs/'
output_dir = '/software/outputs/'


def preprocess():
	if os.path.exists(output_dir):
		rmtree(output_dir)
	os.mkdir(output_dir)

def generate_figure1():
	cprint("Now generating Figure 1's...", "yellow")	

	if os.system("scCloud demuxEM -p 8 --hash-type nuclei-hashing --min-num-genes 200 --generate-diagnostic-plots {input}experiment1_human_st_ADT.csv {input}experiment1_human_st_raw_10x.h5 {output}experiment1_human_st".format(input = input_dir, output = output_dir)):
		sys.exit(1)
	if os.system("scCloud demuxEM -p 8 --hash-type nuclei-hashing --min-num-genes 200 --generate-diagnostic-plots {input}experiment1_human_pbs_ADT.csv {input}experiment1_human_pbs_raw_10x.h5 {output}experiment1_human_pbs".format(input = input_dir, output = output_dir)):
		sys.exit(1)
	if os.system("scCloud aggregate_matrix --restriction Sample:human_st,human_control --attributes Sample --minimum-number-of-genes 100 {input}experiment1_count_matrix.csv {output}experiment1_stonly".format(input = input_dir, output = output_dir)):
		sys.exit(1)
	if os.system("scCloud cluster -p 8 --min-genes 200 --run-louvain --run-tsne {output}experiment1_stonly_10x.h5 {output}experiment1_stonly".format(output = output_dir)):
		sys.exit(1)
	if os.system("scCloud de_analysis -p 8 --fisher {output}experiment1_stonly.h5ad {output}experiment1_stonly.de.xlsx".format(output = output_dir)):
		sys.exit(1)
	if os.system("scCloud annotate_cluster --json-file human_brain {output}experiment1_stonly.h5ad {output}experiment1_stonly.anno.txt".format(output = output_dir)):
		sys.exit(1)

	if os.system('scCloud annotate_cluster --annotation "annotation:1. Oligodendrocyte;2. Glutamatergic neuron;3. Astrocyte;4. GABAergic neuron;5. Glutamatergic neuron;6. Glutamatergic neuron;7. Low quality;8. Endothelial;9. Microglia;10. OPC" {output}experiment1_stonly.h5ad'.format(output = output_dir)):
		sys.exit(1)
	if os.system('scCloud annotate_cluster --annotation "coarse:Oligodendrocyte;Glutamatergic neuron;Astrocyte;GABAergic neuron;Glutamatergic neuron;Glutamatergic neuron;Low quality;Endothelial;Microglia;OPC" {output}experiment1_stonly.h5ad'.format(output = output_dir)):
		sys.exit(1)

	adata = scCloud.tools.read_input('{output}experiment1_stonly.h5ad'.format(output = output_dir), mode = 'a')
	idx = np.isin(adata.obs['louvain_labels'], '7')
	bdata = adata[~idx,:].copy()
	filter_result = scCloud.tools.filter_genes_dispersion(bdata, False)
	bdata_c = scCloud.tools.collect_variable_gene_matrix(bdata, filter_result.gene_subset)
	scCloud.tools.run_pca(bdata_c)
	bdata.obsm['X_pca'] = bdata_c.obsm['X_pca']
	scCloud.tools.run_diffmap(bdata, 'X_pca')
	scCloud.tools.run_louvain(bdata)
	scCloud.tools.run_tsne(bdata, 'X_pca', n_jobs = 8)
	scCloud.tools.run_umap(bdata, 'X_pca')
	bdata.write('{output}experiment1_stonly_lowremove.h5ad'.format(output = output_dir))

	if os.system("scCloud de_analysis -p 8 --fisher {output}experiment1_stonly_lowremove.h5ad {output}experiment1_stonly_lowremove.de.xlsx".format(output = output_dir)):
		sys.exit(1)
	if os.system("scCloud annotate_cluster --json-file human_brain {output}experiment1_stonly_lowremove.h5ad {output}experiment1_stonly_lowremove.anno.txt".format(output = output_dir)):
		sys.exit(1)
	if os.system('scCloud annotate_cluster --annotation "annotation:1. Oligodendrocyte;2. Glutamatergic neuron;3. Astrocyte;4. GABAergic neuron;5. Glutamatergic neuron;6. Glutamatergic neuron;7. Oligodendrocyte;8. Endothelial cell;9. Microglia;10. OPC" {output}experiment1_stonly_lowremove.h5ad'.format(output = output_dir)):
		sys.exit(1)
	if os.system('scCloud annotate_cluster --annotation "coarse:Oligodendrocyte;Glutamatergic neuron;Astrocyte;GABAergic neuron;Glutamatergic neuron;Glutamatergic neuron;Oligodendrocyte;Endothelial cell;Microglia;OPC" {output}experiment1_stonly_lowremove.h5ad'.format(output = output_dir)):
		sys.exit(1)

	adata = scCloud.tools.read_input('{output}experiment1_stonly_lowremove.h5ad'.format(output = output_dir), mode = 'a')
	idx = adata.obs['demux_type'] == 'singlet'
	adata.obs['Donor'] = ''
	adata.obs.loc[idx, 'Donor'] = adata.obs.loc[idx, 'assignment'].astype(str)
	adata.write('{output}experiment1_stonly_lowremove.h5ad'.format(output = output_dir))
	barcodes = np.vectorize(lambda x: x[9:])(adata.obs_names[adata.obs['Sample'] == 'human_st'])
	np.savetxt('{output}experiment1_stonly_barcodes.txt'.format(output = output_dir), barcodes, fmt = "%s", delimiter = '\n')
	barcodes_demuxlet = [x + '-1' for x in barcodes]
	np.savetxt('{output}experiment1_stonly_barcodes_for_demuxlet.txt'.format(output = output_dir), barcodes_demuxlet, fmt = "%s", delimiter = '\n')

	if os.system("scCloud plot scatter --attributes coarse --wspace 0.6 {output}experiment1_stonly_lowremove.h5ad {output}figure_1b.pdf".format(output = output_dir)):
		sys.exit(1)
	cprint("Figure 1b is generated.", "green")
	if os.system("scCloud plot scatter --attributes Sample --alpha 0.5 --wspace 0.6 {output}experiment1_stonly_lowremove.h5ad {output}figure_1c.pdf".format(output = output_dir)):
		sys.exit(1)
	cprint("Figure 1c is generated.", "green")
	if os.system("scCloud plot composition --cluster-labels coarse --attribute Sample --style normalized --not-stacked --bottom 0.5 --wspace 0.4 {output}experiment1_stonly_lowremove.h5ad {output}figure_1d.pdf".format(output = output_dir)):
		sys.exit(1)
	cprint("Figure 1d is generated.", "green")
	
	scCloud.plotting.plot_qc_violin(adata, 'gene', '{output}figure_1e.pdf'.format(output = output_dir), xattr = 'coarse', hue = 'Sample', xlabel = 'Cell type', xtick_font = 5, split = True, figsize = (7, 5), linewidth = 0.5)
	cprint("Figure 1e is generated.", "green")
	
	if os.system("scCloud plot scatter --restriction Donor:S1HuF,S2HuM,S3HuF,S4HuM,S5HuF,S6HuM,S7HuF,S8HuM --show-background --attributes Donor --alpha 0.5 --wspace 0.6 {output}experiment1_stonly_lowremove.h5ad {output}figure_1f.pdf".format(output = output_dir)):
		sys.exit(1)
	cprint("Figure 1f is generated.", "green")
	if os.system("scCloud plot composition --restriction Donor:S1HuF,S2HuM,S3HuF,S4HuM,S5HuF,S6HuM,S7HuF,S8HuM --cluster-labels coarse --attribute Donor --style normalized --not-stacked --bottom 0.5 --wspace 0.4 {output}experiment1_stonly_lowremove.h5ad {output}figure_1g.pdf".format(output = output_dir)):
		sys.exit(1)
	cprint("Figure 1g is generated.", "green")

	cprint("All Figure 1's are generated!", "yellow")

def update_var_names(data):
	gsyms = data.var_names.values
	dup_ids = Counter()
	for i in range(gsyms.size):
		idn = dup_ids[gsyms[i]]
		dup_ids[gsyms[i]] += 1
		if idn > 0:
			gsyms[i] = gsyms[i] + ".{}".format(idn)	
	data.var_names = pd.Index(gsyms)

def apply_filter(data, mito_prefixes):
	if np.unique(data.var_names).size != data.var_names.size:
		update_var_names(data)

	data.obs['n_genes'] = data.X.getnnz(axis = 1)
	data.obs['n_counts'] = data.X.sum(axis = 1).A1

	def startswith(name):
		for prefix in mito_prefixes:
			if name.startswith(prefix):
				return True
		return False

	mito_genes = [name for name in data.var_names if startswith(name)]
	data.obs['percent_mito'] = data[:, mito_genes].X.sum(axis=1).A1 / np.maximum(data.obs['n_counts'].values, 1.0)

	obs_index = np.logical_and.reduce((data.obs['n_genes'] >= 200, data.obs['n_genes'] < 6000, data.obs['percent_mito'] < 0.1))
	data._inplace_subset_obs(obs_index)

asc2bin = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}

def barcode_asc2bin(barcode):
	int_code = 0
	for x in barcode:
		int_code = (int_code << 2) + asc2bin[x]
	return int_code

def load_human_mouse_mix_data(demux_h5ad_file, molecule_info_file):
	data = scCloud.tools.read_input(demux_h5ad_file, mode = 'a')
	M = data.shape[1]
	barcode_dict = {barcode_asc2bin(x) : i for i, x in enumerate(data.obs_names)}
	with tables.open_file(molecule_info_file) as h5_in:
		reads = h5_in.get_node('/reads').read()
		idx = reads > 0
		reads = reads[idx]
		barcodes = h5_in.get_node('/barcode').read()[idx]
		genes = h5_in.get_node('/gene').read()[idx]
	mat = np.zeros(data.shape, dtype = int)
	for i in range(reads.size):
		barcode_id = barcode_dict.get(barcodes[i], None)
		if barcode_id is not None:
			assert genes[i] < M
			if reads[i] > 1:
				mat[barcode_id][genes[i]] += 1
	data.X = csr_matrix(mat)
	apply_filter(data, ['mm10_premrna___mt-', 'GRCh38_premrna_MT-'])

	mm10_genes = [name for name in data.var_names if name.startswith('mm10_premrna')]
	grch38_genes = [name for name in data.var_names if name.startswith('GRCh38_premrna')]
	data.obs['mm10_n_genes'] = (data[:, mm10_genes].X > 0).sum(axis = 1).A1
	data.obs['grch38_n_genes'] = (data[:, grch38_genes].X > 0).sum(axis = 1).A1
	data.obs['mm10_n_counts'] = data[:, mm10_genes].X.sum(axis = 1).A1
	data.obs['grch38_n_counts'] = data[:, grch38_genes].X.sum(axis = 1).A1

	return data

def filter_csv(input_file, barcodes, output_file):
	with open(input_file) as fin, open(output_file, 'w') as fout:
		cols = pd.Index(next(fin).strip().split(',')[1:])
		idx = cols.isin(barcodes)
		fout.write(',{}\n'.format(','.join(cols[idx])))
		for line in fin:
			fields = line.strip().split(',')
			samples = np.array(fields[1:])
			fout.write('{},{}\n'.format(fields[0], ','.join(samples[idx])))

def generate_figure2():
	cprint("Now generating Figure 2's...", "yellow")

	if os.system("scCloud demuxEM -p 8 --hash-type nuclei-hashing --min-num-genes 200 --generate-diagnostic-plots {input}experiment2_mouse_pbs_ADT.csv {input}experiment2_mouse_pbs_raw_10x.h5 {output}experiment2_mouse_pbs".format(input = input_dir, output = output_dir)):
		sys.exit(1)

	adata = scCloud.tools.read_input('{output}experiment2_mouse_pbs_demux.h5ad'.format(output = output_dir), mode = 'a')
	scCloud.tools.update_var_names(adata)
	scCloud.tools.filter_data(adata, mito_prefix = 'mt-', min_genes = 200)
	scCloud.tools.log_norm(adata, 1e5)
	adata.write('{output}experiment2_mouse_pbs_demux_filter.h5ad'.format(output = output_dir))
	scCloud.demuxEM.plot_violin(adata, {'gene' : 'Xist'}, '{output}Figure_2b.pdf'.format(output = output_dir), title = 'Xist')
	cprint("Figure 2b is generated.", "green")

	if os.system("scCloud demuxEM -p 8 --hash-type nuclei-hashing --min-num-genes 200 --generate-diagnostic-plots {input}experiment3_human_mouse_pbs_ADT.csv {input}experiment3_human_mouse_pbs_raw_10x.h5 {output}experiment3_human_mouse_pbs".format(input = input_dir, output = output_dir)):
		sys.exit(1)

	adata = load_human_mouse_mix_data('{output}experiment3_human_mouse_pbs_demux.h5ad'.format(output = output_dir), '{input}experiment3_human_mouse_pbs_molecule_info.h5'.format(input = input_dir))
	adata.write('{output}experiment3_human_mouse_pbs_robust.h5ad'.format(output = output_dir))
	sdp.plot_human_vs_mouse(adata, '{output}figure_2c.pdf'.format(output = output_dir), 'count')
	cprint("Figure 2c is generated.", "green")

	if os.system("mv {output}experiment3_human_mouse_pbs.background_probabilities.bar.pdf {output}figure_2d.pdf".format(output = output_dir)):
		sys.exit(1)
	cprint("Figure 2d is generated.", "green")

	adata = scCloud.tools.read_input('{output}experiment1_human_st_demux.h5ad'.format(output = output_dir), mode = 'a')
	barcodes = np.loadtxt('{output}experiment1_stonly_barcodes.txt'.format(output = output_dir), dtype = 'str', delimiter = '\n')
	bdata = adata[barcodes,].copy()
	demuxlet = pd.read_csv('{input}experiment1_demuxlet.best'.format(input = input_dir), header = 0, index_col = 0, sep = '\t')
	demuxlet.index = demuxlet.index.map(lambda x: x[:-2])
	demuxlet['demux_type'] = 'doublet'
	idx_amb = demuxlet['BEST'].map(lambda x: x.startswith('AMB'))
	demuxlet.loc[idx_amb, 'demux_type'] = 'unknown'
	idx_sng = (demuxlet['PRB.DBL'] <= 0.99) & (~idx_amb)
	demuxlet.loc[idx_sng, 'demux_type'] = 'singlet'
	demuxlet['assignment'] = ''
	demuxlet.loc[idx_sng, 'assignment'] = demuxlet.loc[idx_sng, 'SNG.1ST']
	idx_dbt = np.isin(demuxlet['demux_type'], 'doublet')
	demuxlet.loc[idx_dbt, 'assignment'] = [x + ',' + y for x, y in zip(demuxlet.loc[idx_dbt, 'DBL.1ST'].values, demuxlet.loc[idx_dbt, 'DBL.2ND'].values)]
	bdata.obs[['demuxlet_demux_type', 'demuxlet_assignment']] = demuxlet.loc[bdata.obs_names, ['demux_type', 'assignment']]
	bdata.write('{output}experiment1_stonly_dds.h5ad'.format(output = output_dir)) # dds: demuxEM, demuxlet, and seurat

	demuxEM = bdata.obs['demux_type'].astype(str).values
	idx_sng = demuxEM == 'singlet'
	demuxEM[idx_sng] = bdata.obs.loc[idx_sng, 'assignment'].astype(str).values

	demuxlet = bdata.obs['demuxlet_demux_type'].astype(str).values
	idx_sng = demuxlet == 'singlet'
	demuxlet[idx_sng] = bdata.obs.loc[idx_sng, 'demuxlet_assignment'].astype(str).values
	
	sdp.plot_heatmap(demuxEM, demuxlet, '{output}figure_2e.pdf'.format(output = output_dir), xlabel = 'Demuxlet', ylabel = 'DemuxEM')
	cprint("Figure 2e is generated.", "green")


	filter_csv('{input}experiment1_human_st_ADT.seurat.csv'.format(input = input_dir), bdata.obs_names, '{output}experiment1_human_st_ADT.seurat.selected.csv'.format(output = output_dir))
	if os.system("Rscript ./helper.R"):
		sys.exit(1)

	df_seurat = pd.read_csv('{output}experiment1_human_st_demux.seurat.txt'.format(output = output_dir), header = None, index_col = 0, sep = '\t')
	df_seurat.columns = ['demux_type', '1st', '2nd']
	df_seurat['demux_type'] = df_seurat['demux_type'].apply(lambda x: x.lower() if x != 'Negative' else 'unknown')
	df_seurat['assignment'] = [row['1st'][:5] + ',' + row['2nd'][:5] if row['demux_type'] == 'doublet' else (row['1st'][:5] if row['demux_type'] == 'singlet' else '') for name, row in df_seurat.iterrows()]

	bdata.obs[['seurat_demux_type', 'seurat_assignment']] = df_seurat.loc[bdata.obs_names, ['demux_type', 'assignment']]
	bdata.write('{output}experiment1_stonly_dds.h5ad'.format(output = output_dir))

	seurat = bdata.obs['seurat_demux_type'].astype(str).values
	idx_sng = seurat == 'singlet'
	seurat[idx_sng] = bdata.obs.loc[idx_sng, 'seurat_assignment'].astype(str).values

	sdp.plot_heatmap(seurat, demuxlet, '{output}figure_2f.pdf'.format(output = output_dir), xlabel = 'Demuxlet', ylabel = 'Seurat')
	cprint("Figure 2f is generated.", "green")

	if os.system("scCloud demuxEM -p 8 --hash-type nuclei-hashing --min-num-genes 200 --generate-diagnostic-plots {input}experiment4_human_st_500_ADT.csv {input}experiment4_human_st_500_raw_10x.h5 {output}experiment4_human_st_500".format(input = input_dir, output = output_dir)):
		sys.exit(1)
	if os.system("scCloud demuxEM -p 8 --hash-type nuclei-hashing --min-num-genes 200 --generate-diagnostic-plots {input}experiment4_human_st_1500_ADT.csv {input}experiment4_human_st_1500_raw_10x.h5 {output}experiment4_human_st_1500".format(input = input_dir, output = output_dir)):
		sys.exit(1)
	if os.system("scCloud demuxEM -p 8 --hash-type nuclei-hashing --min-num-genes 200 --generate-diagnostic-plots {input}experiment4_human_st_3000_ADT.csv {input}experiment4_human_st_3000_raw_10x.h5 {output}experiment4_human_st_3000".format(input = input_dir, output = output_dir)):
		sys.exit(1)
	if os.system("scCloud demuxEM -p 8 --hash-type nuclei-hashing --min-num-genes 200 --generate-diagnostic-plots {input}experiment4_human_st_4500_ADT.csv {input}experiment4_human_st_4500_raw_10x.h5 {output}experiment4_human_st_4500".format(input = input_dir, output = output_dir)):
		sys.exit(1)

	if os.system("scCloud aggregate_matrix --attributes Sample --select-only-singlets --minimum-number-of-genes 200 {input}experiment4_count_matrix.csv {output}experiment4_human_st".format(input = input_dir, output = output_dir)):
		sys.exit(1)
	if os.system("scCloud cluster -p 8 --min-genes 200 --run-louvain --run-tsne {output}experiment4_human_st_10x.h5 {output}experiment4_human_st".format(output = output_dir)):
		sys.exit(1)
	if os.system("scCloud de_analysis -p 8 --fisher {output}experiment4_human_st.h5ad {output}experiment4_human_st.de.xlsx".format(output = output_dir)):
		sys.exit(1)
	if os.system("scCloud annotate_cluster --json-file human_brain {output}experiment4_human_st.h5ad {output}experiment4_human_st.anno.txt".format(output = output_dir)):
		sys.exit(1)

	adata = scCloud.tools.read_input('{output}experiment4_human_st.h5ad'.format(output = output_dir), mode = 'a')
	idx = np.isin(adata.obs['louvain_labels'], '5')
	bdata = adata[~idx,:].copy()
	filter_result = scCloud.tools.filter_genes_dispersion(bdata, False)
	bdata_c = scCloud.tools.collect_variable_gene_matrix(bdata, filter_result.gene_subset)
	scCloud.tools.run_pca(bdata_c)
	bdata.obsm['X_pca'] = bdata_c.obsm['X_pca']
	scCloud.tools.run_diffmap(bdata, 'X_pca')
	scCloud.tools.run_louvain(bdata)
	scCloud.tools.run_tsne(bdata, 'X_pca', n_jobs = 8)
	bdata.write('{output}experiment4_human_st_lowremove.h5ad'.format(output = output_dir))

	if os.system("scCloud de_analysis -p 8 --fisher {output}experiment4_human_st_lowremove.h5ad {output}experiment4_human_st_lowremove.de.xlsx".format(output = output_dir)):
		sys.exit(1)
	if os.system("scCloud annotate_cluster --json-file human_brain {output}experiment4_human_st_lowremove.h5ad {output}experiment4_human_st_lowremove.anno.txt".format(output = output_dir)):
		sys.exit(1)
	if os.system('scCloud annotate_cluster --annotation "annotation:1. Oligodendrocyte;2. Astrocyte;3. Glutamatergic neuron;4. Glutamatergic neuron;5. Glutamatergic neuron;6. Glutamatergic neuron;7. GABAergic neuron;8. Glutamatergic neuron;9. Glutamatergic neuron;10. GABAergic neuron;11. Astrocyte;12. OPC;13. Glutamatergic neuron;14. Endothelial cell;15. Microglia;16. Glutamatergic neuron;17. GABAergic neuron;18. Oligodendrocyte;19. Glutamatergic neuron" {output}experiment4_human_st_lowremove.h5ad'.format(output = output_dir)):
		sys.exit(1)
	if os.system('scCloud annotate_cluster --annotation "coarse:Oligodendrocyte;Astrocyte;Glutamatergic neuron;Glutamatergic neuron;Glutamatergic neuron;Glutamatergic neuron;GABAergic neuron;Glutamatergic neuron;Glutamatergic neuron;GABAergic neuron;Astrocyte;OPC;Glutamatergic neuron;Endothelial cell;Microglia;Glutamatergic neuron;GABAergic neuron;Oligodendrocyte;Glutamatergic neuron" {output}experiment4_human_st_lowremove.h5ad'.format(output = output_dir)):
		sys.exit(1)
	if os.system("scCloud plot scatter --attributes coarse --wspace 0.6 {output}experiment4_human_st_lowremove.h5ad {output}figure_2g.pdf".format(output = output_dir)):
		sys.exit(1)
	cprint("Figure 2g is generated.", "green")

	adata = scCloud.tools.read_input('{output}experiment4_human_st_lowremove.h5ad'.format(output = output_dir), mode = 'a')
	scCloud.plotting.plot_qc_violin(adata, 'gene', '{output}figure_2h.pdf'.format(output = output_dir), xattr = 'coarse', hue = 'Sample', xlabel = 'Cell type', xtick_font = 5, linewidth = 0, figsize = (7, 4))
	cprint("Figure 2h is generated.", "green")

	if os.system("scCloud plot scatter --attributes Sample --alpha 0.5 --wspace 0.6 {output}experiment4_human_st_lowremove.h5ad {output}figure_2i.pdf".format(output = output_dir)):
		sys.exit(1)
	cprint("Figure 2i is generated.", "green")
	if os.system("scCloud plot composition --cluster-labels coarse --attribute Sample --style normalized --not-stacked --bottom 0.5 --wspace 0.4 {output}experiment4_human_st_lowremove.h5ad {output}figure_2j.pdf".format(output = output_dir)):
		sys.exit(1)
	cprint("Figure 2j is generated.", "green")

	cprint("All Figure 2's are generated!", "yellow")


def generate_figure_s1():
	cprint("Now generating Figure S1's...", "yellow")

	if os.system("scCloud aggregate_matrix --attributes Sample --minimum-number-of-genes 100 {input}experiment1_count_matrix.csv {output}experiment1_all".format(input = input_dir, output = output_dir)):
		sys.exit(1)
	if os.system("scCloud cluster -p 8 --min-genes 200 --run-louvain --run-tsne {output}experiment1_all_10x.h5 {output}experiment1_all".format(output = output_dir)):
		sys.exit(1)
	if os.system("scCloud de_analysis -p 8 --fisher {output}experiment1_all.h5ad {output}experiment1_all.de.xlsx".format(output = output_dir)):
		sys.exit(1)
	if os.system("scCloud annotate_cluster --json-file human_brain {output}experiment1_all.h5ad {output}experiment1_all.anno.txt".format(output = output_dir)):
		sys.exit(1)
	if os.system('scCloud annotate_cluster --annotation "annotation:1. Glutamatergic neuron;2. Glutamatergic neuron;3. Glutamatergic neuron;4. Astrocyte;5. Oligodendrocyte;6. Oligodendrocyte;7. GABAergic neuron;8. GABAergic neuron;9. Low quality;10. Endothelial;11. Microglia;12. OPC" {output}experiment1_all.h5ad'.format(output = output_dir)):
		sys.exit(1)

	adata = scCloud.tools.read_input('{output}experiment1_all.h5ad'.format(output = output_dir), mode = 'a')
	idx = np.isin(adata.obs['annotation'], '9. Low quality')
	bdata = adata[~idx,:].copy()
	filter_result = scCloud.tools.filter_genes_dispersion(bdata, False)
	bdata_c = scCloud.tools.collect_variable_gene_matrix(bdata, filter_result.gene_subset)
	scCloud.tools.run_pca(bdata_c)
	bdata.obsm['X_pca'] = bdata_c.obsm['X_pca']
	scCloud.tools.run_diffmap(bdata, 'X_pca')
	scCloud.tools.run_louvain(bdata)
	scCloud.tools.run_tsne(bdata, 'X_pca', n_jobs = 8)
	scCloud.tools.run_umap(bdata, 'X_pca')
	bdata.write('{output}experiment1_all_lowremove.h5ad'.format(output = output_dir))

	if os.system("scCloud de_analysis -p 8 --fisher {output}experiment1_all_lowremove.h5ad {output}experiment1_all_lowremove.de.xlsx".format(output = output_dir)):
		sys.exit(1)
	if os.system("scCloud annotate_cluster --json-file human_brain {output}experiment1_all_lowremove.h5ad {output}experiment1_all_lowremove.anno.txt".format(output = output_dir)):
		sys.exit(1)
	if os.system('scCloud annotate_cluster --annotation "annotation:1. Glutamatergic neuron;2. Astrocyte;3. Glutamatergic neuron;4. Oligodendrocyte;5. Oligodendrocyte;6. Glutamatergic neuron;7. Glutamatergic neuron;8. GABAergic neuron;9. GABAergic neuron;10. Microglia;11. OPC;12. Endothelial cell" {output}experiment1_all_lowremove.h5ad'.format(output = output_dir)):
		sys.exit(1)
	if os.system('scCloud annotate_cluster --annotation "coarse:Glutamatergic neuron;Astrocyte;Glutamatergic neuron;Oligodendrocyte;Oligodendrocyte;Glutamatergic neuron;Glutamatergic neuron;GABAergic neuron;GABAergic neuron;Microglia;OPC;Endothelial cell" {output}experiment1_all_lowremove.h5ad'.format(output = output_dir)):
		sys.exit(1)
	if os.system("scCloud plot scatter --attributes coarse,Sample --alpha 1.0,0.5 --wspace 0.6 {output}experiment1_all_lowremove.h5ad {output}figure_s1a.pdf".format(output = output_dir)):
		sys.exit(1)
	cprint("Figure S1a is generated.", "green")

	adata = scCloud.tools.read_input('{output}experiment1_all_lowremove.h5ad'.format(output = output_dir), mode = 'a')
	bdata = adata[adata.obs['Sample'] != 'human_pbs', ].copy()
	scCloud.plotting.plot_qc_violin(bdata, 'gene', '{output}figure_s1b.pdf'.format(output = output_dir), xattr = 'coarse', hue = 'Sample', xlabel = 'Cell type', xtick_font = 5, split = True, figsize = (7, 5), linewidth = 0.5)
	cprint("Figure S1b is generated.", "green")
	bdata = adata[adata.obs['Sample'] != 'human_st', ].copy()
	scCloud.plotting.plot_qc_violin(bdata, 'gene', '{output}figure_s1c.pdf'.format(output = output_dir), xattr = 'coarse', hue = 'Sample', xlabel = 'Cell type', xtick_font = 5, split = True, figsize = (7, 5), linewidth = 0.5)
	cprint("Figure S1c is generated.", "green")

	cprint("All Figure S1's are generated!", "yellow")

def generate_figure_s2():
	cprint("Now generating Figure S2's...", "yellow")

	data_500 = scCloud.tools.read_input('{output}experiment4_human_st_500_demux.h5ad'.format(output = output_dir), mode = 'a')
	gt = scCloud.tools.read_input('{output}experiment4_human_st_500_demux_10x.h5'.format(output = output_dir))
	mito_genes = [name for name in data_500.var_names if name.startswith('MT-')]
	data_500.obs['percent_mito'] = gt[data_500.obs_names,][:,mito_genes].X.sum(axis=1).A1 / np.maximum(data_500.obs['n_counts'].values, 1.0)
	idx = (data_500.obs['n_genes'] < 6000) & (data_500.obs['percent_mito'] < 0.1)
	data_500 = data_500[idx,].copy()
	data_500.write('{output}experiment4_human_st_500_demux_filter.h5ad'.format(output = output_dir))
	scCloud.demuxEM.plot_rna_hist(data_500, '{output}figure_s2a.pdf'.format(output = output_dir))
	cprint("Figure S2a is generated.", "green")

	data_1500 = scCloud.tools.read_input('{output}experiment4_human_st_1500_demux.h5ad'.format(output = output_dir), mode = 'a')
	gt = scCloud.tools.read_input('{output}experiment4_human_st_1500_demux_10x.h5'.format(output = output_dir))
	data_1500.obs['percent_mito'] = gt[data_1500.obs_names,][:,mito_genes].X.sum(axis=1).A1 / np.maximum(data_1500.obs['n_counts'].values, 1.0)
	idx = (data_1500.obs['n_genes'] < 6000) & (data_1500.obs['percent_mito'] < 0.1)
	data_1500 = data_1500[idx,].copy()
	data_1500.write('{output}experiment4_human_st_1500_demux_filter.h5ad'.format(output = output_dir))
	scCloud.demuxEM.plot_rna_hist(data_1500, '{output}figure_s2b.pdf'.format(output = output_dir))
	cprint("Figure S2b is generated.", "green")

	data_3000 = scCloud.tools.read_input('{output}experiment4_human_st_3000_demux.h5ad'.format(output = output_dir), mode = 'a')
	gt = scCloud.tools.read_input('{output}experiment4_human_st_3000_demux_10x.h5'.format(output = output_dir))
	data_3000.obs['percent_mito'] = gt[data_3000.obs_names,][:,mito_genes].X.sum(axis=1).A1 / np.maximum(data_3000.obs['n_counts'].values, 1.0)
	idx = (data_3000.obs['n_genes'] < 6000) & (data_3000.obs['percent_mito'] < 0.1)
	data_3000 = data_3000[idx,].copy()
	data_3000.write('{output}experiment4_human_st_3000_demux_filter.h5ad'.format(output = output_dir))
	scCloud.demuxEM.plot_rna_hist(data_3000, '{output}figure_s2c.pdf'.format(output = output_dir))
	cprint("Figure S2c is generated.", "green")

	data_4500 = scCloud.tools.read_input('{output}experiment4_human_st_4500_demux.h5ad'.format(output = output_dir), mode = 'a')
	gt = scCloud.tools.read_input('{output}experiment4_human_st_4500_demux_10x.h5'.format(output = output_dir))
	data_4500.obs['percent_mito'] = gt[data_4500.obs_names,][:,mito_genes].X.sum(axis=1).A1 / np.maximum(data_4500.obs['n_counts'].values, 1.0)
	idx = (data_4500.obs['n_genes'] < 6000) & (data_4500.obs['percent_mito'] < 0.1)
	data_4500 = data_4500[idx,].copy()
	data_4500.write('{output}experiment4_human_st_4500_demux_filter.h5ad'.format(output = output_dir))
	scCloud.demuxEM.plot_rna_hist(data_4500, '{output}figure_s2d.pdf'.format(output = output_dir))   
	cprint("Figure S2d is generated.", "green")

	cprint("All Figure S2's are generated!", "yellow")

def generate_figure_s3():
	cprint("Now generating Figure S3's...", "yellow")

	s3exp.main()

	barcodes = np.loadtxt('{output}experiment1_stonly_barcodes.txt'.format(output = output_dir), dtype = 'str', delimiter = '\n')
	adata = scCloud.tools.read_input('{output}experiment1_human_st_demux.h5ad'.format(output = output_dir), mode = 'a')
	adata[barcodes].write('{output}experiment1_human_st_demux_fig1b.h5ad'.format(output = output_dir))
	plot_down_sampling('{output}experiment1_human_st_demux_fig1b.h5ad'.format(output = output_dir), '{output}experiment1_human_st_ADTs.h5ad'.format(output = output_dir), '{output}figure_s3c.pdf'.format(output = output_dir), n_threads = 8)
	cprint("figure S3c is generated.", "green")

	cprint("All Figure S3's are generated!", "yellow")

def main():
	preprocess()

	# Figure 1
	generate_figure1()

	# Figure 2
	generate_figure2()

	# Figure S1
	generate_figure_s1()

	# Figure S2
	generate_figure_s2()

	# Figure S3
	generate_figure_s3()


if __name__ == "__main__":
	main()
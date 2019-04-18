import scCloud
import os, sys
import numpy as np
import pandas as pd
from termcolor import cprint

input_dir = '/software/inputs/'
output_dir = '/software/outputs/'

def gen_table_s1():
	cprint("Now generating Table S1...", "yellow")
	exp1_data = scCloud.tools.read_input('{out}experiment1_all.h5ad'.format(out = output_dir))

	nuclei_counts = exp1_data.obs['Sample'].value_counts()
	mean_info = exp1_data.obs.groupby(['Sample']).mean()
	median_info = exp1_data.obs.groupby(['Sample']).median()

	f = open("{out}table_s1.txt".format(out = output_dir), "w")
	f.write("Experiment\tStaining Buffer\tNuclei loading concentration per ul\tNumber of nuclei passing QC (no low-quality cluster removal)\tMean # of UMI\tMedian # of UMI\tMean # of gene\tMedian # of gene\n")
	f.write("1\tST-SB\t500\t{nuclei}\t{mean_umi:.0f}\t{median_umi:.0f}\t{mean_gene:.0f}\t{median_gene:.0f}\n".format(nuclei = nuclei_counts['human_st'], mean_umi = mean_info['n_counts']['human_st'], median_umi = median_info['n_counts']['human_st'], mean_gene = mean_info['n_genes']['human_st'], median_gene = median_info['n_genes']['human_st']))
	f.write("1\tPBS-SB\t500\t{nuclei}\t{mean_umi:.0f}\t{median_umi:.0f}\t{mean_gene:.0f}\t{median_gene:.0f}\n".format(nuclei = nuclei_counts['human_pbs'], mean_umi = mean_info['n_counts']['human_pbs'], median_umi = median_info['n_counts']['human_pbs'], mean_gene = mean_info['n_genes']['human_pbs'], median_gene = median_info['n_genes']['human_pbs']))
	f.write("1\tNA-control\t500\t{nuclei}\t{mean_umi:.0f}\t{median_umi:.0f}\t{mean_gene:.0f}\t{median_gene:.0f}\n".format(nuclei = nuclei_counts['human_control'], mean_umi = mean_info['n_counts']['human_control'], median_umi = median_info['n_counts']['human_control'], mean_gene = mean_info['n_genes']['human_control'], median_gene = median_info['n_genes']['human_control'] + 0.5))

	exp2_data = scCloud.tools.read_input('{out}experiment2_mouse_pbs_demux_filter.h5ad'.format(out = output_dir))

	f.write("2\tPBS-SB\t500\t{nuclei}\t{mean_umi:.0f}\t{median_umi:.0f}\t{mean_gene:.0f}\t{median_gene:.0f}\n".format(nuclei = exp2_data.shape[0], mean_umi = exp2_data.obs['n_counts'].mean(), median_umi = int(exp2_data.obs['n_counts'].median() + 0.5), mean_gene = exp2_data.obs['n_genes'].mean(), median_gene = exp2_data.obs['n_genes'].median()))

	exp3_data = scCloud.tools.read_input('{out}experiment3_human_mouse_pbs_robust.h5ad'.format(out = output_dir))
	f.write("3\tPBS-SB\t500\t{nuclei}\t{mean_umi:.0f}\t{median_umi:.0f}\t{mean_gene:.0f}\t{median_gene:.0f}\n".format(nuclei = exp3_data.shape[0], mean_umi = exp3_data.obs['n_counts'].mean(), median_umi = exp3_data.obs['n_counts'].median(), mean_gene = exp3_data.obs['n_genes'].mean(), median_gene = exp3_data.obs['n_genes'].median()))

	exp4_g1_data = scCloud.tools.read_input('{out}experiment4_human_st_500_demux_filter.h5ad'.format(out = output_dir))
	f.write("4\tST-SB\t500\t{nuclei}\t{mean_umi:.0f}\t{median_umi:.0f}\t{mean_gene:.0f}\t{median_gene:.0f}\n".format(nuclei = exp4_g1_data.shape[0], mean_umi = exp4_g1_data.obs['n_counts'].mean(), median_umi = int(exp4_g1_data.obs['n_counts'].median() + 0.5), mean_gene = exp4_g1_data.obs['n_genes'].mean(), median_gene = int(exp4_g1_data.obs['n_genes'].median() + 0.5)))

	exp4_g2_data = scCloud.tools.read_input('{out}experiment4_human_st_1500_demux_filter.h5ad'.format(out = output_dir))
	f.write("4\tST-SB\t1500\t{nuclei}\t{mean_umi:.0f}\t{median_umi:.0f}\t{mean_gene:.0f}\t{median_gene:.0f}\n".format(nuclei = exp4_g2_data.shape[0], mean_umi = exp4_g2_data.obs['n_counts'].mean(), median_umi = exp4_g2_data.obs['n_counts'].median(), mean_gene = exp4_g2_data.obs['n_genes'].mean(), median_gene = exp4_g2_data.obs['n_genes'].median()))

	exp4_g3_data = scCloud.tools.read_input('{out}experiment4_human_st_3000_demux_filter.h5ad'.format(out = output_dir))
	f.write("4\tST-SB\t3000\t{nuclei}\t{mean_umi:.0f}\t{median_umi:.0f}\t{mean_gene:.0f}\t{median_gene:.0f}\n".format(nuclei = exp4_g3_data.shape[0], mean_umi = exp4_g3_data.obs['n_counts'].mean(), median_umi = exp4_g3_data.obs['n_counts'].median(), mean_gene = exp4_g3_data.obs['n_genes'].mean(), median_gene = exp4_g3_data.obs['n_genes'].median()))
	
	exp4_g4_data = scCloud.tools.read_input('{out}experiment4_human_st_4500_demux_filter.h5ad'.format(out = output_dir))
	f.write("4\tST-SB\t4500\t{nuclei}\t{mean_umi:.0f}\t{median_umi:.0f}\t{mean_gene:.0f}\t{median_gene:.0f}\n".format(nuclei = exp4_g4_data.shape[0], mean_umi = exp4_g4_data.obs['n_counts'].mean(), median_umi = exp4_g4_data.obs['n_counts'].median(), mean_gene = exp4_g4_data.obs['n_genes'].mean(), median_gene = exp4_g4_data.obs['n_genes'].median()))


	f.close()

def gen_table_s2():
	cprint("Now generating Table S2...", "yellow")
	adata = scCloud.tools.read_input('{out}experiment1_stonly_dds.h5ad'.format(out = output_dir))
	demuxEM_counts = adata.obs['demux_type'].value_counts()
	demuxlet_counts = adata.obs['demuxlet_demux_type'].value_counts()
	seurat_counts = adata.obs['seurat_demux_type'].value_counts()
	n_rows = demuxEM_counts.shape[0]

	f = open("{out}table_s2.txt".format(out = output_dir), "w")

	f.write("Method\tSinglet\tMultiple\tUnknown\tTotal\tMultiple rate = Multiple / Total * 100%\n")
	f.write("demuxEM\t{singlet_count}\t{doublet_count}\t{unknown_count}\t{total}\t{rate:.1%}\n".format(singlet_count = demuxEM_counts['singlet'], doublet_count = demuxEM_counts['doublet'], unknown_count = demuxEM_counts['unknown'], total = np.sum(demuxEM_counts), rate = demuxEM_counts['doublet'] / np.sum(demuxEM_counts)))
	f.write("demuxlet\t{singlet_count}\t{doublet_count}\t{unknown_count}\t{total}\t{rate:.1%}\n".format(singlet_count = demuxlet_counts['singlet'], doublet_count = demuxlet_counts['doublet'], unknown_count = demuxlet_counts['unknown'], total = np.sum(demuxlet_counts), rate = demuxlet_counts['doublet'] / np.sum(demuxlet_counts)))
	f.write("Seurat\t{singlet_count}\t{doublet_count}\t{unknown_count}\t{total}\t{rate:.1%}\n".format(singlet_count = seurat_counts['singlet'], doublet_count = seurat_counts['doublet'], unknown_count = seurat_counts['unknown'], total = np.sum(seurat_counts), rate = seurat_counts['doublet'] / np.sum(seurat_counts)))


	f.close()

def gen_table_s3():
	cprint("Now generating Table S3...", "yellow")
	data_500 = scCloud.tools.read_input('{output}experiment4_human_st_500_demux_filter.h5ad'.format(output = output_dir), mode = 'a')
	g1_counts = data_500.obs['demux_type'].value_counts()


	data_1500 = scCloud.tools.read_input('{output}experiment4_human_st_1500_demux_filter.h5ad'.format(output = output_dir), mode = 'a')
	g2_counts = data_1500.obs['demux_type'].value_counts()
	
	data_3000 = scCloud.tools.read_input('{output}experiment4_human_st_3000_demux_filter.h5ad'.format(output = output_dir), mode = 'a')
	g3_counts = data_3000.obs['demux_type'].value_counts()

	data_4500 = scCloud.tools.read_input('{output}experiment4_human_st_4500_demux_filter.h5ad'.format(output = output_dir), mode = 'a')
	g4_counts = data_4500.obs['demux_type'].value_counts()


	f = open("{out}table_s3.txt".format(out = output_dir), "w")

	f.write("\t Nuclei loading concentration\n")
	f.write("Nuclei Type\t500nuc/ul\t1500nuc/ul\t3000nuc/ul\t4500nuc/ul\n")
	f.write("Singlet\t{g1}\t{g2}\t{g3}\t{g4}\n".format(g1 = g1_counts['singlet'], g2 = g2_counts['singlet'], g3 = g3_counts['singlet'], g4 = g4_counts['singlet']))
	f.write("Multiplet\t{g1}\t{g2}\t{g3}\t{g4}\n".format(g1 = g1_counts['doublet'], g2 = g2_counts['doublet'], g3 = g3_counts['doublet'], g4 = g4_counts['doublet']))
	f.write("Unassigned\t{g1}\t{g2}\t{g3}\t{g4}\n".format(g1 = g1_counts['unknown'], g2 = g2_counts['unknown'], g3 = g3_counts['unknown'], g4 = g4_counts['unknown']))
	f.write("Total number of nuclei\t{g1}\t{g2}\t{g3}\t{g4}\n".format(g1 = np.sum(g1_counts), g2 = np.sum(g2_counts), g3 = np.sum(g3_counts), g4 = np.sum(g4_counts)))

	f.close()

def main():
	gen_table_s1()
	gen_table_s2()
	gen_table_s3()

	cprint("Now generating Table S4...", "yellow")
	if os.system("python generate_de_results.py {out}".format(out = output_dir)):
		sys.exit(1)

if __name__ == "__main__":
	main()
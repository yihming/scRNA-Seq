API
===

``scrtools`` can also be used as a python package. Import ``scrtools`` by::

	import scrtools

Tools:
------

**Aggregate channel-specific count matrices**

.. autosummary::
	:toctree: .

	tools.aggregate_10x_matrices

**Preprocess**

.. autosummary::
	:toctree: .

	tools.read_input
	tools.update_var_names
	tools.filter_data
	tools.log_norm
	tools.run_pca
	tools.run_rpca
	tools.get_anndata_for_subclustering

**Batch correction**

.. autosummary::
	:toctree: .

	tools.set_group_attribute
	tools.estimate_adjustment_matrices
	tools.filter_genes_dispersion
	tools.collect_variable_gene_matrix
	tools.correct_batch_effects

**Diffusion map**

.. autosummary::
	:toctree: .

	tools.run_diffmap
	tools.run_pseudotime_calculation

**Cluster algorithms**

.. autosummary::
	:toctree: .

	tools.run_louvain
	tools.run_hdbscan
	tools.run_kmeans
	tools.run_approximated_louvain

**Visualization algorithms**

.. autosummary::
	:toctree: .

	tools.run_tsne
	tools.run_fitsne
	tools.run_umap
	tools.run_force_directed_layout

**Differential expression analysis**

.. autosummary::
	:toctree: .

	tools.run_de_analysis

Annotate clusters:
------------------

.. autosummary::
	:toctree: .

	annotate_cluster.annotate_clusters

Plotting:
---------

**Static plots**

.. autosummary::
	:toctree: .

	plotting.plot_composition
	plotting.plot_scatter
	plotting.plot_scatter_groups
	plotting.plot_scatter_genes
	plotting.plot_scatter_gene_groups
	plotting.plot_heatmap

**Interactive plots**

.. autosummary::
	:toctree: .

	plotting.scatter
	plotting.scatter_real
	plotting.scatter3d
	plotting.scatter3d_real


Miscellaneous:
--------------

.. autosummary::
	:toctree: .

	misc.search_genes
	misc.search_de_genes
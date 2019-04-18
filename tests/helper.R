library(Seurat)

input_dir = "/software/inputs/"
output_dir = "/software/outputs/"

seurat_adt = read.csv(paste(output_dir, "experiment1_human_st_ADT.seurat.selected.csv", sep=""), sep = ",", header = TRUE, row.names = 1)
seurat_adt = seurat_adt[1:8,]

dummy_rna_file = paste(input_dir, "experiment1_human_st_raw_GRCh38_premrna", sep="")
dummy_rna = Read10X(data.dir = dummy_rna_file)

joint_overload = intersect(colnames(dummy_rna),colnames(seurat_adt))
seurat_rna = dummy_rna[,joint_overload]
seurat_hto = as.matrix(seurat_adt[,joint_overload])


seurat_obj = CreateSeuratObject(raw.data = seurat_rna)
seurat_obj = NormalizeData(seurat_obj, display.progress = FALSE)
seurat_obj = FindVariableGenes(seurat_obj, do.plot = F, display.progress = FALSE)
seurat_obj = ScaleData(seurat_obj, genes.use = seurat_obj@var.genes, display.progress = FALSE)

seurat_obj = SetAssayData(seurat_obj, assay.type = "HTO", slot = "raw.data", new.data = seurat_hto)
seurat_obj = NormalizeData(seurat_obj, assay.type = "HTO", normalization.method = "genesCLR", display.progress = FALSE)
seurat_obj = HTODemux(seurat_obj, assay.type = "HTO", positive_quantile = 0.99, print.output = FALSE)

seurat_obj = SetAllIdent(seurat_obj,id = "hto_classification_global")
demux_type = as.data.frame(seurat_obj@ident)
seurat_obj = SetAllIdent(seurat_obj,id = "hash_maxID")
max_id = as.data.frame(seurat_obj@ident)
seurat_obj = SetAllIdent(seurat_obj,id = "hash_secondID")
second_id = as.data.frame(seurat_obj@ident)
df = cbind(demux_type, max_id, second_id)
write.table(df, file = paste(output_dir, "experiment1_human_st_demux.seurat.txt", sep=""), quote = FALSE, sep = "\t", col.names = FALSE)


PDAC_agg <- AggregateExpression(obj2,
                                group.by =c( "orig.ident","seurat_clusters","SCT_snn_res.0.5_annotation","fov"),
                                assays = "SCT",
                                slot = "counts",
                                return.seurat = T)
View(PDAC_agg@meta.data)

### View and Extract Data

View(PDAC_agg@assays$SCT$data)
c_names = colnames(PDAC_agg@assays$SCT$data)

###Parse Column Names and Create Metadata

m1 = do.call(rbind, strsplit(c_names, "_"))
print(head(m1))
meta_text = cbind(c_names, m1)
colnames(meta_text) = c ("Sample", "Cluster", "CellType","fov")
head(meta_text)

### Create Sample_CellType fov Column

Sample_celltype_fov = paste(meta_text[, "Cluster"], meta_text[,"CellType"],meta_text[,"fov"], sep="_")
meta_text = cbind(meta_text, Sample_celltype_fov);
colnames(meta_text)
head(meta_text)

### Define Output Folder and File Paths
output_folder <- "/home/rstudio/"
Project_Name <- "PDAC_PseudoBulk"

### Write Data to Files
Matrix_file = paste(output_folder, Project_Name, "Matrix.txt", sep="");
write.table(PDAC_agg@assays$SCT$data, file = Matrix_file, sep="\t", row.names = TRUE,col.names=NA)
head(PDAC_agg@assays$SCT$data)

meta_table = paste(output_folder, Project_Name, "_meta.txt", sep="");
write.table(meta_text, file = meta_table, sep="\t", row.names=F)

### Write out the matrix ###
outputmatrix = matrix_file
outputmeta = meta_file

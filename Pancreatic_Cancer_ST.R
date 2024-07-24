
library(tidyverse)
library(Seurat)
library(data.table)

exp = fread("/home/rstudio/SMI_images/T1-b/T1-b_exprMat_file.csv") %>%
  mutate(row_name = paste0("pan_",fov,"_",cell_ID))%>%
  column_to_rownames("row_name")%>%
  select(-fov,-cell_ID) %>%
  t()

metadata = fread("/home/rstudio/SMI_images/T1-b/T1-b_metadata_file.csv")%>%
  mutate(row_name = paste0("pan_",fov,"_",cell_ID))%>%
  column_to_rownames("row_name") 

obj2 = CreateSeuratObject(exp, project = "pan", assay = "RNA", min.cells = 1, min.features = 20, meta.data = metadata)
View(obj2@meta.data)

unique(metadata$fov)

fov1 = metadata %>%
  filter(fov == 1)

fov1_clusters = obj2@meta.data %>%
  filter(fov == 1)%>%
  ggplot()+geom_point(aes(x = CenterX_local_px, y = CenterY_local_px, color = seurat_clusters))+
  coord_fixed()
fov1 %>%
  ggplot()+
  geom_point(aes(x = CenterX_global_px, y = CenterY_global_px, color = seurat_clusters))+
  coord_fixed()

fov2 = metadata %>%
  filter(fov == 2)
         
fov2 %>%
  ggplot()+
  geom_point(aes(x = CenterX_local_px, y = CenterY_local_px))+
  coord_fixed()

obj2 = SCTransform(obj2, assay = "RNA", new.assay.name = "SCT")
obj2 = RunPCA(obj2, assay = "SCT", reduction.name = "PCA", npcs = 50)
DimPlot(obj2, reduction = "PCA")

obj2 = RunUMAP(obj2, reduction = "PCA", reduction.name = "UMAP",dims = 1:30, repulsion.strength = 10)
DimPlot(obj2, reduction = "UMAP",label = T)

# Finding clusters --------------------------------------------------------
obj2 = FindNeighbors(obj2, reduction = "PCA", dims = 1:30)
obj2 = FindClusters(obj2, algorithm = 1,resolution = 0.5)


# View the clusters on UMAP
DimPlot(obj2, reduction = "UMAP", group.by = "SCT_snn_res.0.5",label = T)

#ImageDimPlot(obj2, fov = "1" ,axes = T,cols = "glasby")

#Find allmarkers
All.markers = FindAllMarkers(obj2, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
All.markers %>%
  group_by(cluster)%>%
  dplyr::filter()
All.markers %>%
  group_by(cluster)%>%
  dplyr::filter(avg_log2FC > 0.25)%>%
  slice_head(n = 10) %>%
  ungroup()-> top10.markers
DoHeatmap(obj2, features = top10.markers$gene) + NoLegend()

fwrite(All.markers,"/home/rstudio/pdac_markers_2024_07.csv")
library(data.table)

saveRDS(obj2,"/home/rstudio/PDAC_Seurat.rds")

# Rename cell clusters

#Idents(obj2) = "SCT_snn_res.0.5"

obj2 = RenameIdents(obj2, '0'= "Stroma", '1' = "Stroma", '2'="Early Macrophage Ag" , '3'="PDAC Cancer",'4'="Stroma", '5'="PDAC Cancer", '6'="Tcells(CD45,CD3)",'7'="Stroma",'8'="Stroma",'9'="Bcells(CD45,CD3)",'10'="Stroma",'11'= "Late Macrophage Ag",'12'="moDC(CD45,CD68,CD11b)",'13'="Stroma",'14'="Bcell/Mix",'15'= "Monocytes",'16'="Stroma/B/Tcell" )
obj2@meta.data$SCT_snn_res.0.5_annotation = Idents(obj2)

DimPlot(obj2, reduction = "UMAP", label = T)

#dotplot

plot = DotPlot(object = obj2, features = c("IL1R1","IL1RAP","IL32","IL2RG","IL1R2","IL1RL1","IL18","IL7R","IL13RA1","IL15RA","IL6ST","IL11","IL10RA","IL20RA","IL16"), assay="SCT")
plot + coord_flip()

plot = DotPlot(object = obj2, features = c("TNFRSF12A","TNFSF10","TNFRSF21","TNFRSF1A","IFNGR1","IFNGR2","TGFBR2","TGFB3","TGFB1","CSF1R","CSF2RA","CSF3R","LTBR","LTB","FLT1"), assay="SCT")


plot = DotPlot(object = obj2, features = c("CXCL12","CXCR4","CXCL16","CXCL14","CX3CL1","CCL5","CCL4","CCL19"), assay="SCT")


plot + coord_flip()


metafile = cbind(rownames(obj2@meta.data), obj2@meta.data)
metafile = cbind(rownames(Samples), obj2@meta.data)

### Define Output Folder and File Paths
output_folder <- "/home/rstudio/"
Project_Name <- "Pancreatic_cancer"
meta_file = paste(output_folder, Project_Name, "_meta.txt", sep="");
write.table(metafile, file = meta_file, sep="\t",row.names = F)




matrix_file = paste(output_folder, Project_Name, "exp.txt", sep="");

### Write Data to Files
write.table( obj2@assays$SCT$data, file = matrix_file, sep="\t", row.names = TRUE,col.names=NA)

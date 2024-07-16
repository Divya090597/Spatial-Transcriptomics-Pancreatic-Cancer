library(tidyverse)
library(Seurat)

exp = fread("/home/rstudio/SMI_images/T1-b/T1-b_exprMat_file.csv") %>%
  mutate(row_name = paste0("pan_",fov,"_",cell_ID))%>%
  column_to_rownames("row_name")%>%
  select(-fov,-cell_ID) %>%
  t()

metadata = fread("SMI_images/T1-b/T1-b_metadata_file.csv")%>%
  mutate(row_name = paste0("pan_",fov,"_",cell_ID))%>%
  column_to_rownames("row_name") 

obj2 = CreateSeuratObject(exp, project = "pan", assay = "RNA", min.cells = 1, min.features = 20, meta.data = metadata)
View(obj2@meta.data)

unique(metadata$fov)

s1_fov1 = metadata %>%
  filter(fov == 1)

obj2@meta.data %>%
  filter(fov == 1)%>%
  ggplot()+geom_point(aes(x = CenterX_local_px, y = CenterY_local_px, color = seurat_clusters))+
  coord_fixed()
s1_fov1 %>%
  ggplot()+
  geom_point(aes(x = CenterX_global_px, y = CenterY_global_px))+
  coord_fixed()

s2_fov1 = metadata %>%
  filter(fov == 1,
         cell_ID == 2)
s2_fov1 %>%
  ggplot()+
  geom_point(aes(x = CenterX_local_px, y = CenterY_local_px))+
  coord_fixed()

obj2 = SCTransform(obj2, assay = "RNA", new.assay.name = "SCT")
obj2 = RunPCA(obj2, assay = "SCT", reduction.name = "PCA", npcs = 50)
DimPlot(obj2, reduction = "PCA")

obj2 = RunUMAP(obj2, reduction = "PCA", reduction.name = "UMAP",dims = 1:30, repulsion.strength = 10)
DimPlot(obj2, reduction = "UMAP")
# Finding clusters --------------------------------------------------------
obj2 = FindNeighbors(obj2, reduction = "PCA", dims = 1:30)
obj2 = FindClusters(obj2, algorithm = 1,resolution = 0.3)


# View the clusters on UMAP
DimPlot(obj2, reduction = "UMAP", group.by = "SCT_snn_res.0.3")
ImageDimPlot(obj2, fov = "1" ,axes = T,cols = "glasby")

obj2@meta.data%>%
  filter()

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

fwrite(All.markers,"/home/rstudio/Pan_gene.markers.csv")
library(data.table)

saveRDS(obj2,"/home/rstudio/PDAC_Seurat.rds")
 
# Rename cell clusters
obj2 = RenameIdents(obj2, '0'="Stellate",'1'="myeloid",'2'="Ductalcells", '3'="Bcells",'4'="Ductal", '5'="Endothelial cells", '6'="Tcells",'7'="Stellate",'8'="Plasma",'9'="Dendritic_cells",'10'="Mast_cells",'11'= "B cells",'12'="Plasma",'13'="Fibroblasts",'14'="B cells",'15'= "Monocytes",'16'="Stellate" )
DimPlot(obj2, reduction = "UMAP", label = T)

View(obj2@meta.data)

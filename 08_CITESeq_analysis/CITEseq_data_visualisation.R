# CRYO PAPER
# CID44971-2 CHUNKS
# 
# 01 SET UP ------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(stringr)

# library(ggplot2)
# library(cowplot)
# library(reshape2)
# library(stringr)
# library(data.table)

# 02 LOAD DATA ---------------------------------------------------------------

seurat_10X <- readRDS("/share/ScratchGeneral/sunwu/projects/cryopreservation_paper/seurat_v3/07_CITESeq/seurat_CID44971CHUNKS2/Output/Rdata/03_seurat_object_processed.RData")

# 03 DIRECTORIES -------------------------------------------------------------

setwd("/share/ScratchGeneral/sunwu/projects/cryopreservation_paper/seurat_v3/07_CITESeq/")
dir.create("01_annotation_CID44971CHUNKS2")
setwd("01_annotation_CID44971CHUNKS2")


# 04 PLOTS --------------------------------------------

dir.create("01_dimplots_featureplots")

temp_png_function <-
  function(x) {
    png(
      file = (x), 
      width = 6, 
      height = 6, 
      res = 300, 
      units = 'in'
    )
  }

temp_dimplot <- DimPlot(
  object = seurat_10X,
  pt.size = 0.1,
  reduction = "UMAPB",
  group.by = "PC_B_res.1.2", 
  label = T
)
temp_png_function("01_dimplots_featureplots/01_dimplot_UMAPB_res.1.2.png")
print(temp_dimplot)
dev.off()


temp_png_function <-
  function(x) {
    png(
      file = (x), 
      width = 12, 
      height = 6, 
      res = 300, 
      units = 'in'
    )
  }

temp_vlnPlot <- VlnPlot(object = seurat_10X,
                        features = c("nCount_RNA",
                                     "nFeature_RNA",
                                     "percent.mito"),
                        pt.size = 0,
                        group.by = "PC_B_res.1.2",
                        log=T
                        )

temp_png_function("01_dimplots_featureplots/02_metrics.png")
print(temp_vlnPlot)
dev.off()

temp_png_function <-
  function(x) {
    png(
      file = (x), 
      width = 14, 
      height = 14, 
      res = 300, 
      units = 'in'
    )
  }

temp_featureplot <- FeaturePlot(object = seurat_10X,
                                features = c("ACTB",
                                             "PTPRC",
                                             "CD3D",
                                             "CD68",
                                             "MKI67",
                                             "MS4A1",
                                             "EPCAM",
                                             "KRT14",
                                             "PECAM1", 
                                             "COL1A1",
                                             "PDGFRB",
                                             "MCAM"),
                                pt.size = 0.05, 
                                reduction = "UMAPB"
                                )
temp_png_function("01_dimplots_featureplots/03_featureplot.png")
print(temp_featureplot)
dev.off()

# DEG
Idents(object = seurat_10X) <- "PC_B_res.1.2"
temp_cluster_allmarkers <- FindAllMarkers(
  only.pos = T,
  object = seurat_10X,
  min.pct = 0.25, 
  logfc.threshold = 0.25,
  min.diff.pct = 0.1, 
  test.use = 'MAST', 
  print.bar = T
) 

temp_cluster_allmarkers <- arrange(temp_cluster_allmarkers,
                                   (cluster),
                                   desc(avg_logFC))      

write.csv(temp_cluster_allmarkers,
          paste0("01_dimplots_featureplots/04_FindAllMarkers.csv"))
 

temp_png_function <-
  function(x) {
    png(
      file = (x), 
      width = 12, 
      height = 14, 
      res = 300, 
      units = 'in'
    )
  }

temp_genes_for_heatmap <- 
  (temp_cluster_allmarkers %>% group_by(cluster) %>% top_n(10,
                                                           avg_logFC))

seurat_10X@meta.data[,"PC_B_res.1.2"] <- 
  factor(seurat_10X@meta.data[,"PC_B_res.1.2"],
         levels=str_sort(unique(seurat_10X@meta.data[,"PC_B_res.1.2"]),numeric = T))


seurat_10X <- ScaleData(object = seurat_10X,
                        verbose = FALSE,
                        features = as.vector(temp_genes_for_heatmap$gene))

temp_png_function(paste0("01_dimplots_featureplots/05_DoHeatMap.png"))
print(DoHeatmap(
  object = seurat_10X,
  features = temp_genes_for_heatmap$gene,
  group.by = "PC_B_res.1.2",
  raster=F
))
dev.off()




# 05 FILTER CLUSTERS ---------------------------------------------------------

Idents(seurat_10X) <- "PC_B_res.1.2"
temp_ids <- unique(seurat_10X@meta.data$PC_B_res.1.2)

seurat_10X <- subset(seurat_10X,
                     idents = temp_ids[!temp_ids %in% c(1,2,10)])

# 06 ANNOTATE ----------------------------------------------------------------

Idents(seurat_10X) <- temp_colname
seurat_10X <- RenameIdents(seurat_10X,
                           `0` = "Cancer/Epithelial 1 - ELF5+",
                           `3` = "CD4+ T-cells",
                           `4` = "Cancer/Epithelial 2 - ESR1+",
                           `5` = "CD8+ T-cells",
                           `6` = "Myoepithelial",
                           `7` = "Endothelial",
                           `8` = "Cancer/Epithelial 3 - LTF+",
                           `9` = "Endothelial",
                           `11` = "SMCs",
                           `12` = "Monocyte/Macrophage",
                           `13` = "CAFs",
                           `14` = "Cancer/Epithelial 4 - Cycling",
                           `15` = "Cancer/Epithelial 5")

seurat_10X@meta.data$celltype_major <- seurat_10X@active.ident

temp_png_function <-
  function(x) {
    png(
      file = (x), 
      width = 10, 
      height = 10, 
      res = 300, 
      units = 'in'
    )
  }

temp_dimplot <- DimPlot(
  object = seurat_10X,
  pt.size = 0.1,
  reduction = "UMAPB",
  group.by = "celltype_major", 
  label = T,
  repel = T
)
temp_png_function("01_dimplots_featureplots/06_annotated.png")
print(temp_dimplot)
dev.off()



seurat_10X_newUMAP <-
  RunUMAP(seurat_10X,
          dims =  1:50,
          reduction.key = paste0("UMAPNEW",50,"_"),
          reduction.name = paste0("UMAPNEW",50),
          verbose = F
  )

temp_dimplot <- DimPlot(
  object = seurat_10X_newUMAP,
  pt.size = 0.1,
  reduction = "UMAPNEW50",
  group.by = "celltype_major", 
  label = T,
  repel = T
)
temp_png_function("01_dimplots_featureplots/07_rerunUMAP__annotated.png")
print(temp_dimplot)
dev.off()



# 07 SAVERDS FOR NENAD -----------------------------------------------------------------

dir.create("02_Rdata")
saveRDS(seurat_10X_newUMAP,
        "02_Rdata/04_seurat_object_processed_filtered_annotated.RData")




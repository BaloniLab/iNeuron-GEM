source('~/.Rprofile')

library(Seurat)
library(dplyr)
library(patchwork)
library(SeuratWrappers)
library(Matrix)
library(hdf5r)
library(presto)

setwd("/scratch/negishi/jiang817/Neuron_GEM/Data_collection/ROSMAP_sample_merged")
folder_rosmap <- "/scratch/negishi/jiang817/Neuron_GEM/Data_collection/ROSMAP_sample_merged"

dimension_table <- read.csv('/depot/pbaloni/data/Lab_members/Boyu_Jiang/Neuron_GEM/Extract_neuron_data/ROSMAP_DATA/Step_2_findmarkers/dimsions_rosmap.csv')
samples <- unlist(as.list(dimension_table['Sample_id']))
dimensions <- unlist(as.list(dimension_table['Dim']))
names(dimensions) = NULL

for (i in 1:127){
  print(i)
  sample_id <- samples[i]
  rdsfile_path <- paste0(folder_rosmap, '/', sample_id, '/', sample_id, '_ALRA.rds' )
  print(rdsfile_path)
  Object<-readRDS(rdsfile_path)
  dimension = dimensions[i]
  Object <- FindNeighbors(Object, dims = 1:dimension)
  for (res in c(0.5,0.6,0.7)) {
    Object=FindClusters(Object, resolution = res, algorithm = 1)
    Object <- RunUMAP(Object, dims = 1:dimension)
    print(paste0('-----DimPlot-----', res))
    
    DimPlot_path = paste0(folder_rosmap, '/', sample_id, '/', sample_id,'_', res, '_DimPlot.png' )
    png(DimPlot_path, width = 800, height = 600)
    D_plot <- DimPlot(Object, reduction = "umap")
    print(D_plot)
    dev.off()
    Object.markers <- FindAllMarkers(Object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    marker_path <- paste0(folder_rosmap, '/', sample_id, '/', sample_id,'_', res, '_markers.csv' )
    write.csv(Object.markers, marker_path)
    
    remove_file1 =  paste0(folder_rosmap, '/', sample_id, '/', sample_id,'_', '0.8', '_DimPlot.png' )
    remove_file2 =  paste0(folder_rosmap, '/', sample_id, '/', sample_id,'_', '0.8', '_markers.csv' )
    file.remove(remove_file1)
    file.remove(remove_file2)
    }

}
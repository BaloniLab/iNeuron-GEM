source('~/.Rprofile')
library(Seurat)
library(SeuratData)
library(dplyr)
library(patchwork)
library(SeuratWrappers)
library(Matrix)
library(hdf5r)
library(presto)
library(clustree)
library(stringr)
library(BPCells)
library(future)
#plan("multisession", workers = 4)
#options(future.globals.maxSize = 30000 * 1024^2)

## extract sample folder names
setwd("/scratch/negishi/jiang817/Neuron_GEM/Data_collection/ROSMAP_sample_merged")
folder_rosmap <- "/scratch/negishi/jiang817/Neuron_GEM/Data_collection/ROSMAP_sample_merged"
batch <- c()
files <- list.files(folder_rosmap)
for (i in files){
  batch_group = strsplit(i, '[.]')[[1]][1]
  batch = c(batch, batch_group )
}
individual_batch = unique(batch)
individual_batch

for(i in 1:127){
  print(i)
  sampleID <- individual_batch[i]
  print(sampleID)
  sample_path = paste0(folder_rosmap, '/', sampleID)
  data = Read10X(sampleID)
  Object = CreateSeuratObject(counts = data)
  Object[["percent.mt"]] <- PercentageFeatureSet(Object, pattern = "^MT-")
  Object <- subset(Object, subset = nFeature_RNA > 1000  & nFeature_RNA < 50000 & percent.mt < 25)
  Object <- NormalizeData(Object, normalization.method = "LogNormalize", scale.factor = 10000)
  Object <- RunALRA(Object)
  Object <- FindVariableFeatures(Object, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(Object)
  Object<- ScaleData(Object, features = all.genes)
  Object <- RunPCA(Object, features = VariableFeatures(object = Object))
  Object <- JackStraw(Object, num.replicate = 100)
  Object <- ScoreJackStraw(Object, dims = 1:20)
  
  print('-----JackStrawPlot-----')
  JackStrawPlot_path = paste0(sample_path, '/', sampleID, '_JackstrawPlot.png' )
  png(JackStrawPlot_path, width = 800, height = 600)
  J_plot <- JackStrawPlot(Object, dims = 1:15)
  print(J_plot)
  dev.off()
  print('-----Elbowplot---')
  ElbowPlot_path = paste0(sample_path, '/', sampleID, '_ElbowPlot.png')
  png(ElbowPlot_path, width = 800, height = 600)
  E_plot <- ElbowPlot(Object)
  print(E_plot)
  dev.off()
  rds_path <- paste0(sample_path, '/', sampleID, '_ALRA.rds')
  saveRDS(Object, rds_path)

}




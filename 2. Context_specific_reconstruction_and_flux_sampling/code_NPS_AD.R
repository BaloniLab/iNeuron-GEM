library(spatstat.explore)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(anndata)
library(SeuratObject)
library(biomaRt)
#library(schard)
#library(sceasy)
#install.packages("reticulate", dependencies = TRUE, INSTALL_opts = '--no-lock')
library(reticulate)
reticulate::use_condaenv('/depot/pbaloni/data/Lab_members/Boyu_Jiang/Software/.conda/envs/scrnaseq', required = TRUE)
reticulate::py_config()
library(anndata)

#options(Seurat.object.assay.version = "v3")
FOLDER = "~/BaloniLab_Depot/data/Lab_members/Boyu_Jiang/Neuron_GEM/Case_Studies_manuscript/NPS_AD/Data/Individuals/"
files = list.files(FOLDER)

# convert gene id
adata <-  read_h5ad("~/BaloniLab_Depot/data/Lab_members/Boyu_Jiang/Neuron_GEM/Case_Studies_manuscript/NPS_AD/Data/Individuals/AMPAD_HBCC_0000000002.h5ad")
mart <- useMart("ensembl","hsapiens_gene_ensembl")
list <- listFilters(mart)
symbol_id = adata$var_names
symbol2entrez <- getBM(attributes=c('hgnc_symbol','entrezgene_id'), filters='hgnc_symbol', values=symbol_id, mart=mart)
symbol2entrez <- as.data.frame(symbol2entrez)

infeasible_sample = list()
sample_name = c()
n_EN_list = c()
mean_EN_genes_list = c()
mean_EN_counts_list = c()

n_IN_list = c()
mean_IN_genes_list = c()
mean_IN_counts_list = c()

for (i in 1:1494){
  filename = files[i]
  print(i)
  print(filename)
  samplePATH <- paste0(FOLDER, filename)
  adata = read_h5ad(samplePATH)
  count = adata$raw$X
  colnames(count) <- adata$var_names
  rownames(count) <- gsub("_", "-", adata$obs_names)
  count = t(count)
  meta = adata$obs
  colnames(meta) = gsub("_", "-", colnames(meta))
  obj = CreateSeuratObject(counts = count, meta.data = meta)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- RunALRA(obj)
  obj_seurat = CreateSeuratObject(counts = obj@assays$alra$data, meta.data = meta)
  
  cell = list(unique(unique(obj$class)))
  neuon_exist= all(c("IN", "EN") %in% cell[[1]])
  
  if (neuon_exist){
    EX_neuron <- subset(obj, subset = class == "EN")
    IN_neuron <- subset(obj, subset = class == "IN")
    
    
    
    num_EX <- Dim(EX_neuron)[2]
    num_IN <- Dim(IN_neuron)[2]
  
    check_neuron_num <- num_EX > 2 & num_IN > 2
    
    # extract info
    
    EN_neuron_seurat <- subset(obj_seurat, subset = class == "EN")
    IN_neuron_seurat <- subset(obj_seurat, subset = class == "IN")
    
    EN_n_cells = Dim(EN_neuron_seurat)[2]
    IN_n_cells = Dim(IN_neuron_seurat)[2]
    
    EN_mean_genes = mean(EN_neuron_seurat@meta.data$nFeature_RNA)
    IN_mean_genes = mean(IN_neuron_seurat@meta.data$nFeature_RNA)
    
    EN_mean_counts = mean(EN_neuron_seurat@meta.data$nCount_RNA)
    IN_mean_counts= mean(IN_neuron_seurat@meta.data$nCount_RNA)
    
    sample_name = c(sample_name, filename)
    
    n_EN_list = c(n_EN_list, EN_n_cells)
    mean_EN_genes_list = c(mean_EN_genes_list, EN_mean_genes)
    mean_EN_counts_list = c(mean_EN_counts_list, EN_mean_counts)
    
    n_IN_list = c(n_IN_list, IN_n_cells)
    mean_IN_genes_list = c(mean_IN_genes_list, IN_mean_genes)
    mean_IN_counts_list = c(mean_IN_counts_list, IN_mean_counts)
    
    
    if (check_neuron_num){
      genes_EX <- rownames(EX_neuron)
      profile_EX <- rowSums(EX_neuron@assays$alra$data)
      EX_CPM <- tibble::tibble(genes=genes_EX, EX_neuron = profile_EX)
      EX_CPM <- EX_CPM[[2]]*10^6 / sum(EX_CPM[[2]])
      EX_CPM <- as.data.frame(EX_CPM)
      
      
      genes_IN <- rownames(IN_neuron)
      profile_IN <- rowSums(IN_neuron@assays$alra$data)
      IN_CPM <- tibble::tibble(genes=genes_IN, IN_neuron = profile_IN)
      IN_CPM <- IN_CPM[[2]]*10^6 / sum(IN_CPM[[2]])
      IN_CPM <- as.data.frame(IN_CPM)
      
      
      EX_CPM$symbol <- rownames(EX_CPM)
      EX_CPM_save <- merge(EX_CPM, symbol2entrez,
                           by.x = "symbol",
                           by.y = "hgnc_symbol", 
                           all.x = TRUE,
                           sort = F, 
                           remove.duplicates = T)
      EX_CPM_save$BIGG <- paste0(EX_CPM_save$entrezgene_id, "_AT1")
      
      IN_CPM$symbol <- rownames(IN_CPM)
      IN_CPM_save <- merge(IN_CPM, symbol2entrez,
                           by.x = "symbol",
                           by.y = "hgnc_symbol", 
                           all.x = TRUE,
                           sort = F, 
                           remove.duplicates = T)
      IN_CPM_save$BIGG <- paste0(IN_CPM_save$entrezgene_id, "_AT1")
      
      
      
      
      
      save_folder = '~/BaloniLab_Depot/data/Lab_members/Boyu_Jiang/Neuron_GEM/Case_Studies_manuscript/NPS_AD/Data/expressionSet_CPM/'
      
      samplename = strsplit(filename, '.h5ad')[[1]]
      
      save_path_EX_CPM <- paste0(save_folder, "/", "EN/", samplename, "_EN_CPM_ALRA.tsv")
      save_path_IN_CPM <- paste0(save_folder, "/", "IN/", samplename, "_IN_CPM_ALRA.tsv")
      
      #write.csv(EX_CPM_save, save_path_EX_CPM)
      #write.csv(IN_CPM_save, save_path_IN_CPM)
      
    }else{
      infeasible_sample <- c(infeasible_sample, list(filename))
    }
    
  } else{
    infeasible_sample <- c(infeasible_sample, list(filename))
  }
  
  
}

## Whoel data set 
#sceasy::convertFormat("~/BaloniLab_Depot/data/Lab_members/Boyu_Jiang/Neuron_GEM/Case_Studies_manuscript/NPS_AD/Data/psychAD_snRNAseq_rawCounts.h5ad", 
#                      from="anndata", 
#                      to="seurat",
#                      outFile='~/BaloniLab_Depot/data/Lab_members/Boyu_Jiang/Neuron_GEM/Case_Studies_manuscript/NPS_AD/Data/psychAD_snRNAseq_rawCounts.rds')

for (i in 1:1494){
  filename = files[i]
  print(i)
  print(filename)
  samplePATH <- paste0(FOLDER, filename)
  adata = read_h5ad(samplePATH)
  count = adata$raw$X
  colnames(count) <- adata$var_names
  rownames(count) <- gsub("_", "-", adata$obs_names)
  count = t(count)
  meta = adata$obs
  colnames(meta) = gsub("_", "-", colnames(meta))
  obj = CreateSeuratObject(counts = count, meta.data = meta)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- RunALRA(obj)
  
  cell = list(unique(unique(obj$class)))
  neuon_exist= all(c("IN", "EN") %in% cell[[1]])
  
  if (neuon_exist){
    EX_neuron <- subset(obj, subset = class == "EN")
    IN_neuron <- subset(obj, subset = class == "IN")
    
    num_EX <- Dim(EX_neuron)[2]
    num_IN <- Dim(IN_neuron)[2]
    check_neuron_num <- num_EX > 2 & num_IN > 2
    
    if (check_neuron_num){
      genes_EX <- rownames(EX_neuron)
      profile_EX <- rowSums(EX_neuron@assays$alra$data)
      EX_CPM <- tibble::tibble(genes=genes_EX, EX_neuron = profile_EX)
      EX_CPM <- EX_CPM[[2]]*10^6 / sum(EX_CPM[[2]])
      EX_CPM <- as.data.frame(EX_CPM)
      
      
      genes_IN <- rownames(IN_neuron)
      profile_IN <- rowSums(IN_neuron@assays$alra$data)
      IN_CPM <- tibble::tibble(genes=genes_IN, IN_neuron = profile_IN)
      IN_CPM <- IN_CPM[[2]]*10^6 / sum(IN_CPM[[2]])
      IN_CPM <- as.data.frame(IN_CPM)
      
      
      EX_CPM$symbol <- rownames(EX_CPM)
      EX_CPM_save <- merge(EX_CPM, symbol2entrez,
                           by.x = "symbol",
                           by.y = "hgnc_symbol", 
                           all.x = TRUE,
                           sort = F, 
                           remove.duplicates = T)
      EX_CPM_save$BIGG <- paste0(EX_CPM_save$entrezgene_id, "_AT1")
      
      IN_CPM$symbol <- rownames(IN_CPM)
      IN_CPM_save <- merge(IN_CPM, symbol2entrez,
                           by.x = "symbol",
                           by.y = "hgnc_symbol", 
                           all.x = TRUE,
                           sort = F, 
                           remove.duplicates = T)
      IN_CPM_save$BIGG <- paste0(IN_CPM_save$entrezgene_id, "_AT1")
      
      
      
      
      
      save_folder = '~/BaloniLab_Depot/data/Lab_members/Boyu_Jiang/Neuron_GEM/Case_Studies_manuscript/NPS_AD/Data/expressionSet_CPM/'
      
      samplename = strsplit(filename, '.h5ad')[[1]]
      
      save_path_EX_CPM <- paste0(save_folder, "/", "EN/", samplename, "_EN_CPM_ALRA.tsv")
      save_path_IN_CPM <- paste0(save_folder, "/", "IN/", samplename, "_IN_CPM_ALRA.tsv")
      
      write.csv(EX_CPM_save, save_path_EX_CPM)
      write.csv(IN_CPM_save, save_path_IN_CPM)
      
    }else{
      infeasible_sample <- c(infeasible_sample, list(filename))
    }
    
  } else{
    infeasible_sample <- c(infeasible_sample, list(filename))
  }
  
  
}


## Calaculate the thresold of 75 percentile of average gene expression values of protein-coding genes across all samples
df_IN$row_avg <- rowMeans(df_IN[, -1], na.rm = TRUE)
df_EN$row_avg <- rowMeans(df_EN[, -1], na.rm = TRUE)



ineuron_genes <- read.csv('~/BaloniLab_Depot/data/Lab_members/Boyu_Jiang/Neuron_GEM/Case_Studies_manuscript/NPS_AD/Data/expressionSet_CPM/ineuron_gene_id.csv')
df_IN_ienruon <- df_IN[df_IN$symbol %in% ineuron_genes$gene_id, ]
df_EN_ienruon <- df_EN[df_EN$symbol %in% ineuron_genes$gene_id, ]

IN_75percentile = quantile(df_IN_ienruon$row_avg,0.75)
# IN_75percentile = 44.97852
# IN_75percentile_ineuron = 72.5
# EN_75percentile = quantile(df_EN_ienruon$row_avg,0.75)
# EN_75percentile = 45.6044 
# EN_75percentile_ineuron = 69.1 


IN_data_df <- data.frame(
  IN_n_cells = n_IN_list,
  IN_mean_genes = mean_IN_genes_list,
  IN_mean_counts = mean_IN_counts_list
)

write.csv(IN_data_df, '/depot/pbaloni/data/Lab_members/Boyu_Jiang/Manuscripts/Neuron_GEM/Revision/IN_NPS_AD_summary_alra.csv')

EN_data_df <- data.frame(
  EN_n_cells = n_EN_list,
  EN_mean_genes = mean_EN_genes_list,
  EN_mean_counts = mean_EN_counts_list
)

write.csv(IN_data_df, '/depot/pbaloni/data/Lab_members/Boyu_Jiang/Manuscripts/Neuron_GEM/Revision/EN_NPS_AD_summary_alra.csv')




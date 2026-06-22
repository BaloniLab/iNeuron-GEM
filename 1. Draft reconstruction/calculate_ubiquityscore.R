library(stringr)  

library(data.table)
filtered_cell_meta = read.csv('/depot/pbaloni/data/Lab_members/Boyu_Jiang/Neuron_GEM/Extract_neuron_data/ROSMAP_DATA/Step_1_imputation/filtered_cell_meta.csv')
individual_sample_id = list.files('/scratch/negishi/jiang817/Neuron_GEM/Data_collection/ROSMAP_sample_splited')
mergeed_sample_id = list.files('/depot/pbaloni/data/Lab_members/Boyu_Jiang/Neuron_GEM/Extract_neuron_data/ROSMAP_DATA/Step_3_annoteCells/Extracted_Neuron_counts')


for (i in 112 : 127){
  print(i)
  print(mergeed_sample_id[i])
  merged_counts_path = paste0('/depot/pbaloni/data/Lab_members/Boyu_Jiang/Neuron_GEM/Extract_neuron_data/ROSMAP_DATA/Step_3_annoteCells/Extracted_Neuron_counts/', mergeed_sample_id[i])
  merged_counts = fread(merged_counts_path, header = TRUE)
  batch_id = strsplit(mergeed_sample_id[i], '_Neuron_Counts.csv')[[1]][1]
  batch_id = str_replace_all(batch_id, '-', "_")
  batch_meta = filtered_cell_meta[filtered_cell_meta$libraryBatch == batch_id, ]
  batch_samples = unique(batch_meta$sampleID)
  print(batch_samples)
  for (sample_id in batch_samples){
    print(sample_id)
    batch_one_sample <- batch_meta[batch_meta$sampleID == sample_id, ]
    cell_barcode_one_sample <- batch_one_sample$cellBarcode
    if(grepl("[.]", names(merged_counts)[2])){
      cell_barcode_one_sample <- gsub("-", ".", cell_barcode_one_sample, fixed=TRUE)
    }
    selected_columns = names(merged_counts)[names(merged_counts) %in% cell_barcode_one_sample]
    count_one_sample = merged_counts[, ..selected_columns]
    save_onesample_counts_path = paste0('/scratch/negishi/jiang817/Neuron_GEM/Data_collection/ROSMAP_extracted_fileterd_Neuron_counts/', sample_id, '.csv')
    write.csv(count_one_sample, save_onesample_counts_path)
  }
}





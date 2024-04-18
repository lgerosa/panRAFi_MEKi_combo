#combination into gDR format

library(gDR)
library(gDRimport)

#set working directory
cwd <- "/gstore/home/gerosal/projects/collab/CellBasedDataImports/gDRinternalData/2022_10_05_384w_belva_cobi_synergy_NEW"
setwd(cwd)

#set directory names
data_dir='data'
results_dir='results'
figures_dir='figures'
clines_dir ='clines_annotations'
drugs_dir ='drug_annotations'

### IMPORT THE NRAS DATA

#import data
manifest <- './data/manifest_10052022_384w_belva_cobi_synergy.xlsx'
treatment <- list.files(file.path(cwd, data_dir), "trt_", full.names = TRUE)
raw_data <- './data/Raw_data_10052022_384w_belva_cobi_synergy.xlsx'    

data_imported_1 <- gDR::import_data(manifest, treatment, raw_data, instrument = "Tecan")

#remove statutosporin data
idx <- data_imported_1$DrugName != "Staurosporine" 
data_imported_1 <- data_imported_1[idx,]

#change reference division time for IPC-298
idx <- data_imported_1$CellLineName == "IPC-298"
data_imported_1[idx, c('ReferenceDivisionTime')] <- 60


### IMPORT THE A375 DATA

#import data
manifest <- './data/manifest_07172023_A375_belva_synergy.xlsx'
treatment <- list.files(file.path(cwd, data_dir), "trt_", full.names = TRUE)
raw_data <- './data/07172023_A375_belva_synergy.xlsx'    

data_imported_2 <- gDR::import_data(manifest, treatment, raw_data, instrument = "Tecan")

#remove statutosporin data
idx <- data_imported_2$DrugName != "Staurosporine" 
data_imported_2 <- data_imported_2[idx,]

#change the doubling rate of A375 to the one measured by Fran with incucyte
data_imported_2$ReferenceDivisionTime <- 45

### MERGE

data_imported <- rbind(data_imported_1 , data_imported_2)


### SAVE DATA

MAE <- gDRcore::runDrugResponseProcessingPipeline(data_imported, nested_confounders = "Plate")
saveRDS(MAE, "MAE_20221005_belva_cobi_NEW.RDS")



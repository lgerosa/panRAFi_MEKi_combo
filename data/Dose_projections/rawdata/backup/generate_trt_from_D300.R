#this script imports the KRAS185 dataset with GDC-6036 and GDC-1971 drug dose
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

#generate treatment files from D300 file 
D300 <- './data/D300_10052022_384w_belva_cobi_synergy.tdd'
metadata <- './data/D300_Gnumbers_10052022_384w_belva_cobi_synergy.xlsx'
#note 1: D300 file was manualy changed to remove "<Mode> Concentration </Mode> 
#in the .tdd definition because it was not dealt with by the import_D300 function
#note 2: trt_P1-2 and trt_P5-P8 deleted manually after generation
gDRimport::import_D300(D300, metadata , './data')
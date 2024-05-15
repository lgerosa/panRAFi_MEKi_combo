# Genentech Inc
# by Luca Gerosa and Andrew Goetz
#
# Obtains average PK values 
library(ggplot2)
library(ggnewscale)
library(gridExtra)
library(dplyr)
library(hash)
library(viridis)
library(gDRutils)
library(ggbeeswarm)

#set working directory
cwd <- "/Users/andrewgoetz/Documents/Luca_Projects/panRAFi_MEKi_combo"
setwd(cwd)

#load utility functions
source(paste(cwd,'/scripts/utility.R', sep=''))

#set directory names
data_dir = 'data'
results_dir = 'results'
figures_dir='figures'


### LOAD DATA ###

### load combo data from higher resolution matrixes

#load combo
#combo <- readRDS(file.path(cwd, data_dir, 'Drug_screen', 'Drug_screen_Belva_Cobi_combo_metrics.RDS'))
combo <- readRDS(file.path(cwd, data_dir, 'Dose_projections', 'Zoomed_in_Belva_Cobi_combo_metrics.RDS'))
#load smooth matrix
smooth <- combo$SmoothMatrix
#keep only Belvarafenib and Cobimetinib
smooth <- smooth[smooth$DrugName == 'panRAFi_Belvarafenib' & smooth$DrugName_2 =='MEKi_Cobimetinib',]
#prepare annotations of mutations
anno <- readRDS(file.path(cwd, data_dir, 'Drug_screen', 'Drug_screen_cell_line_annotations.RDS'))
#add annotations
smooth <-  merge(smooth, anno , by='CellLineName',  all.x = TRUE)
#smooth <- filter(smooth,CellLineName == "IPC-298")
# Full time range
time_plot_range <- c(528,528+96) 

# loads belva PK data and filters for desired time range
belva_pk <- data.frame(read.csv(file.path(cwd, data_dir, 'PK_variability', 'belva_sim_pk_data_pruned.csv')))
belva_pk <- dplyr::select(belva_pk,ID,time,IPRED,TRT01P)
belva_pk <- rename(belva_pk, Dose_ID = TRT01P)
belva_pk <- filter(belva_pk,time <= max(time_plot_range), time >= min(time_plot_range))


# loads cobi PK data and filters for desired time range
cobi_pk <- data.frame(read.csv(file.path(cwd, data_dir, 'PK_variability', 'cobi_sim_pk_data_pruned.csv')))
cobi_pk <- dplyr::select(cobi_pk,ID,time,IPREDnormal,TRT01P)
cobi_pk <- rename(cobi_pk, Dose_ID_2 = TRT01P)
#cobi_pk$Dose_ID_2 <- gsub("Cobi","Cobimetinib",cobi_pk$Dose_ID_2)
cobi_pk <- filter(cobi_pk,time <= max(time_plot_range), time >= min(time_plot_range))

# specify which doses to use
belva_Doses_To_Use <- c("Belvarafenib 50mg QD","Belvarafenib 100mg BID" ,"Belvarafenib 200mg BID","Belvarafenib 400mg BID")
cobi_Doses_To_Use <- c("Cobi 20mg QOD" , "Cobi 20mg QD" ,"Cobi 40mg QD","Cobi 60mg QD")
all_dose_combo <- expand.grid(belva_Doses_To_Use, cobi_Doses_To_Use,stringsAsFactors = FALSE)
belva_Doses_To_Use <- all_dose_combo$Var1
cobi_Doses_To_Use <- all_dose_combo$Var2

# specify which cell lines to use
cellLinesToUse <- c("A-375", "IPC-298")
#cellLinesToUse <- c("MEL-JUSO",  "SK-MEL-2",  "SK-MEL-30")

# Finds time for full cycle for each dose combo
#full_cycle_time <- c()
#for (i in 1:length(belva_Doses_To_Use)){
#  full_cycle_time <- c(full_cycle_time,max(get_cycle_length(belva_Doses_To_Use[i]),get_cycle_length(cobi_Doses_To_Use[i]),na.rm=TRUE))
#}
full_cycle_time <- rep(48,length(cobi_Doses_To_Use))
dose_count <- length(full_cycle_time)

#set % FBS in media
used_FBS_perc = 10

#decide which metrics to use
gtf <- list()
choise_gtf <- 1
if (choise_gtf == 0){
  gtf$long <- 'RelativeViability'
  gtf$short <- 'RV'
} else{
  gtf$long <- 'GRvalue' 
  gtf$short <- 'GR'
}



#go through each cell line and drug combination
assay_ID <- c("SmoothMatrix")
title_ID <- c(gtf$long)
field_ID <- c(gtf$long)

colors_fields <- list()
if (gtf$short =='RV'){
  colors_fields[[1]] <- viridis(51)
}else if (gtf$short =='GR') {
  colors_fields[[1]] <- c(rev(rocket(60)[10:60]), viridis(50))
}

mixmax_fields <- c(1,2,2) #1 for GR, Rel, 2 for HSA and Bliss, method to calculate min and max


#calculate free conc
fu_human_Cobi <- 0.052
#Cobi MW: 531.3 g/mol
Cobi_MW <- 531.3

#calculate free conc
fu_human_Belva <- 0.00258
#Belva MW: 478.93 g/mol
Belva_MW <- 478.93


# converts from ng/mL to umol/L
belva_pk$free_Concentration <- (belva_pk$IPRED/Belva_MW)*fu_human_Belva # (ng/mL / g/mol) = 10^-6 mol/L = uMol
cobi_pk$free_Concentration_2 <- (cobi_pk$IPREDnormal/Cobi_MW)*fu_human_Cobi # (ng/mL / g/mol) = 10^-6 mol/L = uMol

# adds drug names
belva_pk$DrugName <- "panRAFi_Belvarafenib"
cobi_pk$DrugName_2 <- "MEKi_Cobimetinib"

#flatten
field <- sprintf('%s_x', gtf$short)
groups <- c("normalization_type")
wide_cols <- c('x')
smooth <- gDRutils::flatten(smooth, groups = groups, wide_cols = wide_cols)

# grabs cell lines which will be used
smooth <- filter(smooth,CellLineName %in% cellLinesToUse)


# Any concentrations above in vitro levels are set to max in vitro level.
# This wont significantly impact results assuming the in vitro value
# was able to come sufficiently close to Emax
# Combines PK data into single dataframe
dose_response <- merge(belva_pk, cobi_pk, by=c("ID","time"))
dose_response <- dose_response[order( dose_response$ID, dose_response$time ),]
dose_response_filtered <- filter(dose_response,dose_response$Dose_ID_2 %in% cobi_Doses_To_Use & dose_response$Dose_ID %in% belva_Doses_To_Use)
dose_response_filtered <- dose_response_filtered[order( dose_response_filtered$ID, dose_response_filtered$time ),]

# specify which doses to use
belva_Doses_To_Use <- c("Belvarafenib 50mg QD","Belvarafenib 100mg BID" ,"Belvarafenib 200mg BID" ,"Belvarafenib 400mg BID")
cobi_Doses_To_Use <- c("Cobi 20mg QOD" , "Cobi 20mg QD" ,"Cobi 40mg QD","Cobi 60mg QD")

doses <- c(belva_Doses_To_Use,cobi_Doses_To_Use)
DDoses <- data.frame(DrugName = character(), dd = double(), dd_min = double(), dd_max = double(), title = character(),Dose_ID = character())
DDoses_var <- data.frame(DrugName = character(), intrinsic_var = double(), extrinsic_var = double(),  title = character(),Dose_ID = character())
for (i in 1:length(belva_Doses_To_Use)){
  belva_pk_filtered <- filter(belva_pk,Dose_ID == belva_Doses_To_Use[i])
  dose_mean <- mean(belva_pk_filtered$free_Concentration)
  low_q <- quantile(belva_pk_filtered$free_Concentration)
  patient_IDs <- unique(belva_pk_filtered$ID)
  patient_variances <- c()
  patient_means <- c()
  for (ID_val in patient_IDs){
    patient_response <- filter(belva_pk_filtered, ID == ID_val)
    patient_variances <- c(patient_variances,var(patient_response$free_Concentration)*(length(patient_response$free_Concentration)-1)/length(patient_response$free_Concentration))
    patient_means <- c(patient_means,mean(patient_response$free_Concentration))
  }
  pop_var <- var(belva_pk_filtered$free_Concentration)*(length(belva_pk_filtered$free_Concentration)-1)/length(belva_pk_filtered$free_Concentration)
  mean_var <- var(patient_means)*(length(patient_means)-1)/length(patient_means) 
  print(mean(patient_variances)/pop_var)
  print(mean_var/pop_var)
  print(mean(patient_variances) + var(patient_means)*(length(patient_means)-1)/length(patient_means) - var(belva_pk_filtered$free_Concentration)*(length(belva_pk_filtered$free_Concentration)-1)/length(belva_pk_filtered$free_Concentration))
  dose_sd <- sd(belva_pk_filtered$free_Concentration)
  title_part <- strsplit(belva_Doses_To_Use[i]," ")
  title <- paste(title_part[[1]][2],"_",title_part[[1]][3],sep="")
  DrugName <- belva_pk_filtered$DrugName[1]
  Dose_ID <- paste(DrugName,"[",title,"]",sep = "")
  DDoses[i,] = c(DrugName,dose_mean,dose_mean-dose_sd,dose_mean+dose_sd,title,Dose_ID)
} 

for (i in 1:length(cobi_Doses_To_Use)){
  
  cobi_pk_filtered <- filter(cobi_pk,Dose_ID_2 == cobi_Doses_To_Use[i])
  dose_mean <- mean(cobi_pk_filtered$free_Concentration)
  dose_sd <- sd(cobi_pk_filtered$free_Concentration)
  
  patient_IDs <- unique(cobi_pk_filtered$ID)
  patient_variances <- c()
  patient_means <- c()
  for (ID_val in patient_IDs){
    patient_response <- filter(cobi_pk_filtered, ID == ID_val)
    patient_variances <- c(patient_variances,var(patient_response$free_Concentration)*(length(patient_response$free_Concentration)-1)/length(patient_response$free_Concentration))
    patient_means <- c(patient_means,mean(patient_response$free_Concentration))
  }
  pop_var <- var(cobi_pk_filtered$free_Concentration)*(length(cobi_pk_filtered$free_Concentration)-1)/length(cobi_pk_filtered$free_Concentration)
  mean_var <- var(patient_means)*(length(patient_means)-1)/length(patient_means) 
  print(mean(patient_variances)/pop_var)
  print(mean_var/pop_var)
  
  title_part <- strsplit(cobi_Doses_To_Use[i]," ")
  title <- paste(title_part[[1]][2],"_",title_part[[1]][3],sep="")
  DrugName <- cobi_pk_filtered$DrugName[1]
  Dose_ID <- paste(DrugName,"[",title,"]",sep = "")
  DDoses[i+length(belva_Doses_To_Use),] = c(DrugName,dose_mean,dose_mean-dose_sd,dose_mean+dose_sd,title,Dose_ID)
} 


file_res <- 'DDoses.csv'
write.csv( DDoses, file.path(cwd, data_dir, 'Dose_projections' ,file_res))

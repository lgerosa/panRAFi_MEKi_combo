# Genentech Inc
# by Luca Gerosa and Andrew Goetz
#
# Project individual patient data onto in-vitro dose responses
library(ggplot2)
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
source(paste(cwd,'/scripts/Andrew_utility.R', sep=''))

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

# Full time range
time_plot_range <- c(528,528+96) 

# loads belva PK data
belva_pk <- data.frame(read.csv(file.path(cwd, data_dir, 'PK_variability', 'belva_sim_pk_data_pruned.csv')))
belva_pk <- select(belva_pk,ID,time,IPRED,TRT01P)
belva_pk <- rename(belva_pk, Dose_ID = TRT01P)
belva_pk <- filter(belva_pk,time <= max(time_plot_range), time >= min(time_plot_range))

#belva_pk <- filter(belva_pk, time < )
# loads cobi PK data
cobi_pk <- data.frame(read.csv(file.path(cwd, data_dir, 'PK_variability', 'cobi_sim_pk_data_pruned.csv')))
cobi_pk <- select(cobi_pk,ID,time,IPREDnormal,TRT01P)
cobi_pk <- rename(cobi_pk, Dose_ID_2 = TRT01P)
cobi_pk <- filter(cobi_pk,time <= max(time_plot_range), time >= min(time_plot_range))

# specify which doses to use
unique(belva_pk$Dose_ID)
unique(cobi_pk$Dose_ID_2)
belva_Doses_To_Use <- c("Belvarafenib 400mg BID", "Belvarafenib 400mg BID","Belvarafenib 100mg BID","Belvarafenib 100mg BID","Belvarafenib 50mg QD","Belvarafenib 100mg BID")
cobi_Doses_To_Use <- c("Cobi 20mg QOD",  "Cobi 20mg QD","Cobi 20mg QOD","Cobi 20mg QD", "Cobi 40mg QD","Cobi 40mg QD")

# Finds time for full cycle for each dose combo
full_cycle_time <- c()
for (i in 1:length(belva_Doses_To_Use)){
  full_cycle_time <- c(full_cycle_time,max(get_cycle_length(belva_Doses_To_Use[i]),get_cycle_length(cobi_Doses_To_Use[i])))
}

dose_count <- length(full_cycle_time)

# specify which cell lines to use
#cellLinesToUse <- c("A-375",     "IPC-298")
cellLinesToUse <- c("MEL-JUSO",  "SK-MEL-2",  "SK-MEL-30")


# specify how many patients to add to violin plots
violin_patient_count <- 100
violin_patients <- 1:violin_patient_count

# specify which patients to plot on heatmaps, 12
patient_colors <- c("magenta","red","orange","yellow", "black","cyan","grey","brown","purple","maroon","magenta","red","orange","yellow", "black","cyan","grey","brown","purple","maroon","magenta","red","orange","yellow", "black","cyan","grey","brown","purple","maroon","magenta","red","orange","yellow", "black","cyan","grey","brown","purple","maroon","magenta","red","orange","yellow", "black","cyan","grey","brown","purple","maroon","magenta","red","orange","yellow", "black","cyan","grey","brown","purple","maroon")
num_patients_to_plot <- 10



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

# Filters doses for specified conditions
cobi_pk <- filter(cobi_pk,cobi_pk$Dose_ID_2 %in% cobi_Doses_To_Use)
belva_pk <- filter(belva_pk,belva_pk$Dose_ID %in% belva_Doses_To_Use)

#convert Concentrations to Concentrations_free using FuDrugs
FuDrugs <- data.frame(read.csv(file.path(cwd, data_dir, 'Dose_projections', 'FuDrugs.csv')))
FuDrugs <- FuDrugs[FuDrugs$FBS_perc==used_FBS_perc, ]

##### Define function to add free drug concentrations to a dataset #####
add_free_Concentrations <- function(dt_data, fu_drugs){
  #add fields for free Concentrations
  dt_data$free_Concentration <- NaN
  dt_data$free_Concentration <- as.numeric(dt_data$free_Concentration)
  dt_data$free_Concentration_2 <- NaN
  dt_data$free_Concentration_2 <- as.numeric(dt_data$free_Concentration_2)
  #add free drug concentrations
  for (i in 1:length(fu_drugs$DrugName)){
    #for Drug
    idx <- dt_data$DrugName == fu_drugs[i, c('DrugName')]
    dt_data[ idx ,c('free_Concentration')] <-  dt_data[idx ,c('Concentration')] * fu_drugs[i, c('fu_FBS')]
    #for Drug_2
    idx <- dt_data$DrugName_2 == fu_drugs[i, c('DrugName')]
    dt_data[ idx ,c('free_Concentration_2')] <-  dt_data[idx ,c('Concentration_2')] * fu_drugs[i, c('fu_FBS')]
  }
  return(dt_data)
}
#add free drug calculations to combo dataset
assay_ID <- c("SmoothMatrix", "HSAExcess", "BlissExcess")
for (i in 1:length(assay_ID)){
  combo[[assay_ID[i]]] <- add_free_Concentrations(combo[[assay_ID[i]]], FuDrugs)
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

#add free concentration vectors  
smooth <- add_free_Concentrations(smooth, FuDrugs)

#belva_pk <- filter(belva_pk,Dose_ID %in% )


### Prompt dataset information ###


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

# Combines PK data into single dataframe
dose_response <- merge(belva_pk, cobi_pk, by=c("ID","time"))

#single_time_DR <- filter(dose_response,time == 714)

dose_response <- dose_response[order( dose_response$ID, dose_response$time ),]

#flatten
field <- sprintf('%s_x', gtf$short)
groups <- c("normalization_type")
wide_cols <- c('x')
smooth <- gDRutils::flatten(smooth, groups = groups, wide_cols = wide_cols)

smooth <- filter(smooth,CellLineName %in% cellLinesToUse)



#go through each cell line and drug combination
assay_ID <- c("SmoothMatrix", "HSAExcess", "BlissExcess")
title_ID <- c(gtf$long, "HSA Excess", "Bliss Excess")
field_ID <- c(gtf$long, gtf$long, gtf$long)
colors_fields <- list()
if (gtf$short =='RV'){
  colors_fields[[1]] <- viridis(51)
}else if (gtf$short =='GR') {
  colors_fields[[1]] <- c(rev(rocket(60)[10:60]), viridis(50))
}

colors_fields[[2]] <- colorRampPalette(c("royalblue3", "royalblue1", "grey95" , "grey95" , "firebrick1", "firebrick3"))(51)
colors_fields[[3]] <- colorRampPalette(c("royalblue3", "royalblue1", "grey95" , "grey95" , "firebrick1", "firebrick3"))(51)

projected_doses <- getProjectedPKEffects(smooth,dose_response,gtf)
uIDs <- unique(projected_doses$ID)
# IDs without values lying outside of interpolation range
vIDs <- c()
for (i in 1:length(uIDs)){
  t<- filter(projected_doses,ID == uIDs[i])
  if (sum(is.na(t$Combo)) == 0){
    vIDs <- c(vIDs,uIDs[i])
  }
}
projected_doses <- filter(projected_doses,ID %in% vIDs)
#filter dt_dd
dt_dd_filt <- projected_doses
uclines_drugs <- unique(dplyr::select(dt_dd_filt, c('CellLineName', 'DrugName', 'DrugName_2')))




patient_ID_to_plot <- hash()
for (k in 1:dose_count){
  projected_doses_filtered <- filter(projected_doses,projected_doses$CellLineName==uclines_drugs[1, c('CellLineName')],projected_doses$Dose_ID==belva_Doses_To_Use[k],projected_doses$Dose_ID_2==cobi_Doses_To_Use[k])
  
  plot_able_IDs <- c()
  
  for (ID_val in unique(projected_doses_filtered$ID)){
    min_invitro_1 <- min(filter(projected_doses_filtered,projected_doses_filtered$ID == ID_val)$free_Concentration)
    max_invitro_1 <- max(filter(projected_doses_filtered,projected_doses_filtered$ID == ID_val)$free_Concentration)
    min_invitro_2 <- min(filter(projected_doses_filtered,projected_doses_filtered$ID == ID_val)$free_Concentration_2)
    max_invitro_2 <- max(filter(projected_doses_filtered,projected_doses_filtered$ID == ID_val)$free_Concentration_2)
    if ((min_invitro_1>min(ud1)) & (max_invitro_1<max(ud1))&(min_invitro_2>min(ud2)) & (max_invitro_2<max(ud2))){
      plot_able_IDs <- c(plot_able_IDs,ID_val)
    }
  }
  
  
  
  drug_1_aucs <- c()
  drug_2_aucs <- c()
  for (ID_val in plot_able_IDs){
    single_patient_dose <- mean(filter(projected_doses_filtered,projected_doses_filtered$ID == ID_val)$free_Concentration)
    drug_1_aucs <- c(drug_1_aucs,single_patient_dose)
    single_patient_dose_2 <- mean(filter(projected_doses_filtered,projected_doses_filtered$ID == ID_val)$free_Concentration_2)
    drug_2_aucs <- c(drug_2_aucs,single_patient_dose_2)
  }
  summed_z_score <- (drug_1_aucs - mean(drug_1_aucs))/sd(drug_1_aucs) + (drug_2_aucs - mean(drug_2_aucs))/sd(drug_2_aucs)
  patients_to_map <- plot_able_IDs[order(summed_z_score)[round(seq(1, length(summed_z_score), length.out = num_patients_to_plot+2)[2:(num_patients_to_plot+1)])]]

  #patients_to_map <- c(low_drug_1_ID,high_drug_1_ID,low_drug_2_ID,high_drug_2_ID)
  #patients_to_map <- c(LD1_LD2,HD1_LD2,HD1_HD2,LD1_HD2)
  patient_ID_to_plot[[toString(k)]] <- patients_to_map
  }

rownames(uclines_drugs) <- NULL
#for each cell line and drug combination
p_matv <- list()
l_matv <- list()
t_matv <- list()
#for each cell line and drug combo
for (i in 1:nrow(uclines_drugs)){
  #extract Rel.vial or GR, HSA and Bliss from QCS_combo
  for (j in 1:length(assay_ID)){
    assay <- assay_ID[j]
    field <- field_ID[j]
    matv <- combo[[assay]]
    field <- gtf$long
    matv <- combo[[assay]]
    idx <- ((matv$CellLineName==uclines_drugs[i, c('CellLineName')]) &
              (matv$DrugName==uclines_drugs[i, c('DrugName')]) &
              (matv$DrugName_2==uclines_drugs[i, c('DrugName_2')]))
    matv <- matv[idx, ]
    if ('normalization_type' %in% colnames(matv)){
      groups <- c("normalization_type")
      wide_cols <- c('x')
      matv <- gDRutils::flatten(matv, groups = groups, wide_cols = wide_cols)
    }
    
    if (nrow(matv)>1) {
      #use free concentrations to plot
      matv <- add_free_Concentrations(matv, FuDrugs)
      matv$Concentration <- matv$free_Concentration
      matv$Concentration_2 <- matv$free_Concentration_2
      matv <- createLogConcentrations(matv)
      if (mixmax_fields[j]==1) {
        mine <- min(c(0.0, min(na.omit(matv[,..field]))))
        maxe <- max(c(1.1, max(na.omit(matv[,..field])))) 
        cdotdose <- 'black'
      } else if (mixmax_fields[j]==2) {
        tope <- max(abs(c(min(na.omit(matv[,..field])), 
                          max(na.omit(matv[,..field])))))
        mine <- min(c(-0.5,-tope))
        maxe <- max(c(0.5, tope))
        cdotdose <- 'black'
      }
      limits <- c(mine,maxe)
      
      #calculate width and height of each tile for geom_tile
      matv$x <- matv$logConcentration
      matv$y <- matv$logConcentration_2
      ux <- unique(matv$x)
      uy <- unique(matv$y)
      diffux <- diff(ux)/2
      diffuy <- diff(uy)/2
      #calculate width left and right (wl, wr) and height up and down (hu, hd)
      matv$wl=plyr::mapvalues(matv$x, ux, c(diffux[1],diffux))
      matv$wr=plyr::mapvalues(matv$x, ux, c(diffux,diffux[length(diffux)]))
      matv$hu=plyr::mapvalues(matv$y, uy, c(diffuy,diffuy[length(diffuy)]))
      matv$hd=plyr::mapvalues(matv$y, uy, c(diffuy[1],diffuy))
      #calculate xmin, xmax, ymin, ymax for geom_rect
      matv$xmin <- (matv$x - matv$wl)
      matv$xmax <- (matv$x + matv$wr)
      matv$ymin <- (matv$y - matv$hd)
      matv$ymax <- (matv$y + matv$hu)
      colors_matv <- colors_fields[[j]]
      p_matv[[length(p_matv)+1]] <- ggplot()  +
        geom_rect(data=matv, 
                  aes_string(xmin = 'xmin', xmax = 'xmax', 
                             ymin = 'ymin', ymax = 'ymax', 
                             fill = field)) +
        labs(title=sprintf('%s', 
                           uclines_drugs[i, c('CellLineName')])) +
        xlab(paste(unique(matv$DrugName),'uM'))+
        ylab(paste(unique(matv$DrugName_2),'uM')) +
        scale_fill_gradientn(colours = colors_matv, limits=limits) +
        ggpubr::theme_pubr() +
        labs(fill = title_ID[j]) +
        theme(text = element_text(size=7),
              axis.text = element_text(size = 7),
              axis.text.x = element_text(angle = 0, hjust = 0.5),
              plot.title = element_text(hjust = 0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "right"
        )
      for (k in 1:dose_count){
        ind_heatmap <- length(l_matv)+1
        l_matv[[ind_heatmap]] <- p_matv[[length(p_matv)]]
        t_matv[[ind_heatmap]] <- assay
        projected_doses_filtered <- filter(projected_doses,projected_doses$CellLineName==uclines_drugs[i, c('CellLineName')],projected_doses$Dose_ID==belva_Doses_To_Use[k],projected_doses$Dose_ID_2==cobi_Doses_To_Use[k])
        
        #patients_to_map <- c(low_drug_1_ID,high_drug_1_ID,low_drug_2_ID,high_drug_2_ID)
        patients_to_map <- patient_ID_to_plot[[toString(k)]]
        for (patient_ID_ind in 1:length(patients_to_map)){
          projected_curve <- filter(projected_doses_filtered,ID == patients_to_map[patient_ID_ind])
          l_matv[[ind_heatmap]] <- l_matv[[ind_heatmap]]+geom_path(data = projected_curve,aes(x=log10(free_Concentration), y=log10(free_Concentration_2)), linewidth = 1,col = patient_colors[patient_ID_ind+(k-1)*length(patients_to_map)])+ggtitle(sprintf("%s\n%s\n%s",uclines_drugs[i, c('CellLineName')],belva_Doses_To_Use[k],cobi_Doses_To_Use[k]))
        }
      }
    }
  }
}
t_matv
drug_1_aucs
match(max(drug_1_aucs[drug_1_aucs<quantile(drug_1_aucs)[2]]),drug_1_aucs)
match(max(drug_1_aucs[drug_1_aucs<(mean(drug_1_aucs)-var(drug_1_aucs)^.5/1.5) & drug_2_aucs<(mean(drug_2_aucs)-var(drug_2_aucs)^.5/1.5) ]),drug_1_aucs)
match(max(drug_1_aucs[drug_1_aucs>(mean(drug_1_aucs)+var(drug_1_aucs)^.5/1.5) & drug_2_aucs<(mean(drug_2_aucs)-var(drug_2_aucs)^.5/1.5) ]),drug_1_aucs)

match(min(drug_1_aucs[drug_1_aucs>(mean(drug_1_aucs)+var(drug_1_aucs)^.5/3) & drug_2_aucs>(mean(drug_2_aucs)+var(drug_2_aucs)^.5/3) ]),drug_1_aucs)
unique(projected_doses_filtered$ID)[match(min(drug_1_aucs[drug_1_aucs>(mean(drug_1_aucs)+var(drug_1_aucs)^.5/3) & drug_2_aucs>(mean(drug_2_aucs)+var(drug_2_aucs)^.5/3) ]),drug_1_aucs)]


min(drug_1_aucs[drug_1_aucs>(mean(drug_1_aucs)+var(drug_1_aucs)^.5/3) & drug_2_aucs>(mean(drug_2_aucs)+var(drug_2_aucs)^.5/3) ])
match(max(drug_1_aucs[drug_1_aucs>(mean(drug_1_aucs)+var(drug_1_aucs)^.5/3) & drug_2_aucs>(mean(drug_2_aucs)+var(drug_2_aucs)^.5/3) ]),drug_1_aucs)
match(max(drug_1_aucs[drug_1_aucs<(mean(drug_1_aucs)-var(drug_1_aucs)^.5/3) & drug_2_aucs>(mean(drug_2_aucs)+var(drug_2_aucs)^.5/3) ]),drug_1_aucs)

projected_doses_filtered_time <- filter(projected_doses_filtered,time == min(time))
drug_1_aucs <- c()
drug_2_aucs <- c()
for (ID_val in unique(projected_doses_filtered$ID)){
  single_patient_dose <- mean(filter(projected_doses_filtered,projected_doses_filtered$ID == ID_val)$free_Concentration)
  single_patient_dose_2 <- mean(filter(projected_doses_filtered,projected_doses_filtered$ID == ID_val)$free_Concentration_2)
  drug_1_aucs <- c(drug_1_aucs,single_patient_dose)
  drug_2_aucs <- c(drug_2_aucs,single_patient_dose_2)
}
quantile(drug_1_aucs)[2]
quantile(drug_1_aucs)[4]
unique(projected_doses_filtered$ID)[match(quantile(drug_1_aucs)[2],drug_1_aucs)]
unique(projected_doses_filtered$ID)[match(quantile(drug_1_aucs)[4],drug_1_aucs)]

unique(projected_doses_filtered$ID)[match(quantile(drug_2_aucs)[2],drug_2_aucs)]
unique(projected_doses_filtered$ID)[match(quantile(drug_2_aucs)[4],drug_2_aucs)]
drug_1_aucs[53]
projected_doses_filtered_time$free_Concentration_2
length(unique(projected_doses_filtered$ID))
#save heatmaps with concentration plots
#file_res <- sprintf('Heatmaps_DA_all_doses_in_one_%s.pdf', gtf$short)
file_res <- sprintf('Single_PK_traj_projection_%s_%d_FBS_%d.pdf', gtf$short,length(cellLinesToUse),used_FBS_perc)
ncol <- dose_count
nrow <- length(l_matv)/dose_count
pdf(file.path(cwd, figures_dir, 'PK_variability' ,file_res), width= 5.5 * ncol , height=5 * nrow  )  
print(grid.arrange(grobs = l_matv, ncol=dose_count))
dev.off()







dt_dd_filt <- projected_doses
uclines_drugs <- unique(dplyr::select(dt_dd_filt, c('CellLineName', 'DrugName', 'DrugName_2')))
rownames(uclines_drugs) <- NULL
#for each cell line and drug combination
cplot <- list()
#for each cell line and drug combo
for (j in 1:dose_count){
  if (j== 1){
    max_effect <- max(filter(projected_doses,time > 528, time <= 528+full_cycle_time[j])$Combo)
    min_effect <- min(filter(projected_doses,time > 528, time <= 528+full_cycle_time[j])$Combo)
  } else{
  max_effect <- max(max_effect,filter(projected_doses,time > 528, time <= 528+full_cycle_time[j])$Combo)
  min_effect <- min(min_effect,filter(projected_doses,time > 528, time <= 528+full_cycle_time[j])$Combo)
  }
}
for (j in 1:dose_count){
  
  projected_doses_filtered <- filter(projected_doses,projected_doses$Dose_ID==belva_Doses_To_Use[j],projected_doses$Dose_ID_2==cobi_Doses_To_Use[j],time > 528, time <= 528+full_cycle_time[j])
  patients_to_map <- unique(projected_doses_filtered$ID)[violin_patients]
  projected_doses_filtered <- filter(projected_doses_filtered,projected_doses_filtered$ID %in% patients_to_map)
  

  cplot[[length(cplot)+1]] <- ggplot(data=projected_doses_filtered, 
                                     aes(x = CellLineName, y = Combo,color = time))+geom_quasirandom()+ggtitle(sprintf("%s\n%s",belva_Doses_To_Use[j],cobi_Doses_To_Use[j])) + ylim(min(min_effect,0),max_effect+.01) 
}
#+ scale_color_viridis(discrete = FALSE)
length(unique(projected_doses_filtered$ID))


length(unique(projected_doses_filtered$ID))
#save heatmaps with concentration plots
#file_res <- sprintf('Heatmaps_DA_all_doses_in_one_%s.pdf', gtf$short)
file_res <- sprintf('Violin_Plot_%s_%d_FBS_%d.pdf', gtf$short,length(cellLinesToUse),used_FBS_perc)
ncol <- length(vplot)
nrow <- 1
pdf(file.path(cwd, figures_dir, 'PK_variability' ,file_res), width= 10 * ncol , height=5 * nrow  )  
print(grid.arrange(grobs = cplot, ncol=length(cplot)))
dev.off()








projected_doses$BlissExcess <- projected_doses$Bliss - projected_doses$Combo




dt_dd_filt <- projected_doses
uclines_drugs <- unique(dplyr::select(dt_dd_filt, c('CellLineName', 'DrugName', 'DrugName_2')))
rownames(uclines_drugs) <- NULL
#for each cell line and drug combination
bplot <- list()
#for each cell line and drug combo
for (j in 1:dose_count){
  if (j== 1){
    max_effect <- max(filter(projected_doses,time > 528, time <= 528+full_cycle_time[j])$BlissExcess)
    min_effect <- min(filter(projected_doses,time > 528, time <= 528+full_cycle_time[j])$BlissExcess)
  } else{
    max_effect <- max(max_effect,filter(projected_doses,time > 528, time <= 528+full_cycle_time[j])$BlissExcess)
    min_effect <- min(min_effect,filter(projected_doses,time > 528, time <= 528+full_cycle_time[j])$BlissExcess)
  }
}
for (j in 1:dose_count){
  
  projected_doses_filtered <- filter(projected_doses,projected_doses$Dose_ID==belva_Doses_To_Use[j],projected_doses$Dose_ID_2==cobi_Doses_To_Use[j],time > 528, time <= 528+full_cycle_time[j])
  patients_to_map <- unique(projected_doses_filtered$ID)[violin_patients]
  projected_doses_filtered <- filter(projected_doses_filtered,projected_doses_filtered$ID %in% patients_to_map)
  
  
  bplot[[length(bplot)+1]] <- ggplot(data=projected_doses_filtered, 
                                     aes(x = CellLineName, y = BlissExcess,color = time))+geom_quasirandom()+ggtitle(sprintf("%s\n%s",belva_Doses_To_Use[j],cobi_Doses_To_Use[j])) + ylim(min(min_effect,0),max_effect+.01) 
}
#+ scale_color_viridis(discrete = FALSE)
length(unique(projected_doses_filtered$ID))


length(unique(projected_doses_filtered$ID))
#save heatmaps with concentration plots
#file_res <- sprintf('Heatmaps_DA_all_doses_in_one_%s.pdf', gtf$short)
file_res <- sprintf('Violin_Plot_Bliss_Excess_%s_%d_FBS_%d.pdf', gtf$short,length(cellLinesToUse),used_FBS_perc)
ncol <- length(vplot)
nrow <- 1
pdf(file.path(cwd, figures_dir, 'PK_variability' ,file_res), width= 10 * ncol , height=5 * nrow  )  
print(grid.arrange(grobs = bplot, ncol=length(bplot)))
dev.off()














dt_dd_filt <- projected_doses
uclines_drugs <- unique(dplyr::select(dt_dd_filt, c('CellLineName', 'DrugName', 'DrugName_2')))
rownames(uclines_drugs) <- NULL
#for each cell line and drug combination
pkplot <- list()
#for each cell line and drug combo

for (j in 1:dose_count){
  projected_doses_filtered <- filter(projected_doses,projected_doses$CellLineName==uclines_drugs[i, c('CellLineName')],projected_doses$Dose_ID==belva_Doses_To_Use[j],projected_doses$Dose_ID_2==cobi_Doses_To_Use[j])
  patients_to_map <- patient_ID_to_plot[[toString(j)]]
  all_dose_values <- log10(select(filter(projected_doses_filtered,projected_doses_filtered$ID %in% patients_to_map),free_Concentration,free_Concentration_2))
  if (j == 1){
  projected_doses_min <- min(all_dose_values,na.rm=TRUE)
  projected_doses_max <- max(all_dose_values,na.rm=TRUE)
  } else{
    projected_doses_min <- min(all_dose_values,projected_doses_min,na.rm=TRUE)
    projected_doses_max <- max(all_dose_values,projected_doses_max,na.rm=TRUE)
  }
}

for (j in 1:dose_count){
  projected_doses_filtered <- filter(projected_doses,projected_doses$CellLineName==uclines_drugs[i, c('CellLineName')],projected_doses$Dose_ID==belva_Doses_To_Use[j],projected_doses$Dose_ID_2==cobi_Doses_To_Use[j])
  patients_to_map <- patient_ID_to_plot[[toString(j)]]
  for (patient_ID_ind in 1:length(patients_to_map)){
    projected_curve <- filter(projected_doses_filtered,ID == patients_to_map[patient_ID_ind])
    pkplot[[patient_ID_ind+(j-1)*length(patients_to_map)]]  <- ggplot() + geom_path(data = projected_curve,aes(x=time, y=log10(free_Concentration_2)), linewidth = 1,linetype="dashed")+ggtitle(sprintf("Patient ID: %d\n%s\n%s",patients_to_map[patient_ID_ind],belva_Doses_To_Use[j],cobi_Doses_To_Use[j])) + theme(
      panel.background = element_rect(fill = "white", colour = "white",
                                      linewidth = 2, linetype = "solid"),
      panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                      colour = "gray85") + theme_classic()
    )
    pkplot[[patient_ID_ind+(j-1)*length(patients_to_map)]]  <- pkplot[[patient_ID_ind+(j-1)*length(patients_to_map)]] + geom_path(data = projected_curve,aes(x=time, y=log10(free_Concentration)), linewidth = 1)+ggtitle(sprintf("Patient ID: %d\n%s\n%s",patients_to_map[patient_ID_ind],belva_Doses_To_Use[j],cobi_Doses_To_Use[j]))+ylim(projected_doses_min,projected_doses_max) + theme_classic()
  #,col = patient_colors[patient_ID_ind+(j-1)*length(patients_to_map)]
    }
}


length(unique(projected_doses_filtered$ID))


length(unique(projected_doses_filtered$ID))
#save heatmaps with concentration plots
#file_res <- sprintf('Heatmaps_DA_all_doses_in_one_%s.pdf', gtf$short)
file_res <- sprintf('PK_trajectories_%d.pdf',length(cellLinesToUse))
ncol <- 4
nrow <- length(pkplot)/4
pdf(file.path(cwd, figures_dir, 'PK_variability' ,file_res), width= 5 * ncol , height=5 * nrow  )  
print(grid.arrange(grobs = pkplot, ncol=4))
dev.off()




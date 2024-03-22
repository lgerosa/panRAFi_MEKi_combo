# Genentech Inc
# by Luca Gerosa and Andrew Goetz
#
#analyse the RAF/MEK screen:
#1. Study dose range of combination of RAF/MEK

library(ggplot2)
library(gridExtra)
library(dplyr)
library(viridis)

#set working directory
cwd <- "/gstore/home/gerosal/projects/work/panRAFi_MEKi_combo"
setwd(cwd)

#load utility functions
source(paste(cwd,'/scripts/utility.R', sep=''))

#set directory names
data_dir = 'data'
results_dir = 'results'
figures_dir='figures'


### LOAD DATA ###

#load combo
combo <- readRDS(file.path(cwd, data_dir, 'Drug_screen', 'Drug_screen_Belva_Cobi_combo_metrics.RDS'))
#load smooth matrix
smooth <- combo$SmoothMatrix
#keep only Belvarafenib and Cobimetinib
smooth <- smooth[smooth$DrugName == 'panRAFi_Belvarafenib' & smooth$DrugName_2 =='MEKi_Cobimetinib',]

#prepare annotations of mutations
anno <- readRDS(file.path(cwd, data_dir, 'Drug_screen', 'Drug_screen_cell_line_annotations.RDS'))
#add annotations
smooth <-  merge(smooth, anno , by='CellLineName')

### Prompt dataset information ###
message('Cell lines: ', length(unique(smooth$CellLineName)))
message('Drug combos: ', nrow(unique(dplyr::select(smooth, c('DrugName','DrugName_2')))))

#decide which metrics to use
gtf <- list()
choise_gtf <- 0
if (choise_gtf == 0){
  gtf$long <- 'RelativeViability'
  gtf$short <- 'RV'
} else{
  gtf$long <- 'GRvalue' 
  gtf$short <- 'GR'
}

#### 1. PERFORM DOSE PROJECTION ON IN VITRO CELL LINE DATA 
## Inputs:
## DDose -> dataframe that defines SA dose regime for treatments
## CGroups -> list that define cell line groups to be assessed
## DAsses -> dataframe that specifies for which SA or Combo we should asses viability values

#set % FBS in media
used_FBS_perc = 10

#load the DDoses definitions and generate unique IDs
DDoses <- data.frame(read.csv(file.path(cwd, data_dir, 'Dose_projections', 'DDoses.csv')))
DDoses$Dose_ID <- sprintf('%s[%s]', DDoses$DrugName, DDoses$title)

#create CGroups (cell line groups) list
CGroups <- list()
fgroup <- 'CellLineName'

#NRAS mutant
idx <- smooth[['NRAS_mut']]=='yes' & smooth[['BRAF_mut']]=='no'
CGroups[['NRAS_mut']] <-  unique(smooth[idx,..fgroup])

#BRAF V600E selected
#idx <- smooth$CellLineName %in% c('A-375', 'WM−266−4', 'SK-MEL-28')
#CGroups[['BRAF_mut']] <-  unique(smooth[idx,..fgroup])

#WT sensitive
#idx <- smooth$CellLineName %in% c('OVCA 420','HT-115')
#CGroups[['WT']] <-  unique(smooth[idx,..fgroup])


#load isobologram data and define quantities to use
clids_drugs <- unique(select(smooth, c('CellLineName','DrugName', 'DrugName_2')))
ord_cols <- c('DrugName', 'DrugName','CellLineName')
clids_drugs <- data.table::setorderv(clids_drugs, ord_cols)

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
##### END FUNCTION #####

#add free concentration vectors  
smooth <- add_free_Concentrations(smooth, FuDrugs)

#merge doses 
DDoses_temp <- DDoses
DDoses_temp$Dose_ID_2 <- DDoses_temp$Dose_ID
DDoses_temp$DrugName_2 <- DDoses_temp$DrugName
clids_drugs <- merge(clids_drugs, dplyr::select(DDoses_temp, c('DrugName','Dose_ID')), by='DrugName',  allow.cartesian=TRUE)
clids_drugs <- merge(clids_drugs, dplyr::select(DDoses_temp, c('DrugName_2','Dose_ID_2')), by='DrugName_2',  allow.cartesian=TRUE)

#create struture to save results 
col_names_dd <- c('CellLineName',
                  'DrugName', 'DrugName_2', 
                  'Dose_ID', 'Dose_ID_2', 'Metric',
                  'DD', 'DD_min', 'DD_max',
                  'DD_2', 'DD_min_2', 'DD_max_2',
                  'SA', 'SA_min', 'SA_max',
                  'SA_2', 'SA_min_2', 'SA_max_2',
                  'Combo', 'Combo_min', 'Combo_max',
                  'HSA', 'HSA_min', 'HSA_max',
                  'Bliss', 'Bliss_min', 'Bliss_max')
dt_dd <- data.frame(matrix(ncol = length(col_names_dd), nrow = 0))
colnames(dt_dd) <- col_names_dd

#do calculations only for cells that are at least in one group
allCellLines <- c()
allNames <- names(CGroups)
for (i in 1:length(allNames)){
  allCellLines <- c(allCellLines, CGroups[[allNames[i]]][[1]])
}
allCellLines <- unique(allCellLines)
idx <- clids_drugs$CellLineName %in% allCellLines
clids_drugs <- clids_drugs[idx,]

#extract SA or combo metrics for each drug and cell line at defined concentration range
for (i in 1:nrow(clids_drugs)) {
  #extract cell lines and drugs
  clid <- clids_drugs[['clid']][i]
  cln <- clids_drugs[['CellLineName']][i]
  drug <- clids_drugs[['Gnumber']][i]
  drug_2 <- clids_drugs[['Gnumber_2']][i]
  dnp <- clids_drugs[['DrugName']][i]
  dnp_2 <- clids_drugs[['DrugName_2']][i]   
  dose_id <- clids_drugs[['Dose_ID']][i]
  dose_id_2 <- clids_drugs[['Dose_ID_2']][i]
  
  #use smooth matrix to get measured values
  idx <- (smooth$clid==clid) &
    (smooth$Gnumber==drug) &
    (smooth$Gnumber_2==drug_2)
  dt_smooth_sub <- smooth[idx,]
  #use interpolation to extract growth metrics at SA and combo doses
  smooth_matrix <- reshape2::acast(dt_smooth_sub, free_Concentration ~ free_Concentration_2, value.var = gtf$long)
  idx <- (dt_smooth_sub$free_Concentration==0)
  x <- sort(dt_smooth_sub$free_Concentration_2[idx])
  idx <- (dt_smooth_sub$free_Concentration_2==0)
  y <- sort(dt_smooth_sub$free_Concentration[idx])
  
  #extract values to interpolate
  idx <- (DDoses$Dose_ID==dose_id)
  dd <- DDoses[idx, c('dd', 'dd_min', 'dd_max')]
  dd <- c(dd$dd, dd$dd_min, dd$dd_max)
  idx <- (DDoses$Dose_ID==dose_id_2)
  dd_2 <- DDoses[idx, c('dd', 'dd_min', 'dd_max')]
  dd_2 <- c(dd_2$dd, dd_2$dd_min, dd_2$dd_max)
  #interpolate at SA and combo areas
  iv <- pracma::interp2(x  = x, 
                        y  = y, 
                        Z = smooth_matrix, 
                        xp = c(dd_2, c(0,0,0), dd_2), 
                        yp = c(c(0,0,0), dd, dd)
  )
  
  #HSA direct calculation
  iv_hsa <- c()
  if((iv[1]<iv[4]) | is.na(iv[4])){
    iv_hsa[1] <- iv[1]
    iv_hsa[2] <- iv[2]
    iv_hsa[3] <- iv[3]
  }else{
    iv_hsa[1] <- iv[4]
    iv_hsa[2] <- iv[5]
    iv_hsa[3] <- iv[6]
  }
  
  #calculate Bliss using gDR
  sa2 <- dt_smooth_sub[dt_smooth_sub[['free_Concentration']] == 0,]
  sa2 <- as.data.frame(dplyr::select(sa2, c('free_Concentration', 'free_Concentration_2', gtf$long)))
  sa1 <- dt_smooth_sub[dt_smooth_sub[['free_Concentration_2']] == 0,]
  sa1 <- as.data.frame(dplyr::select(sa1, c('free_Concentration', 'free_Concentration_2', gtf$long)))
  dt_bliss <- gDRcore::calculate_Bliss(sa1, 'free_Concentration', sa2, 'free_Concentration_2', gtf$long)
  bliss_matrix <- reshape2::acast(dt_bliss, free_Concentration ~ free_Concentration_2, value.var = 'metric')
  y_bliss <- sort(sa1$free_Concentration)
  x_bliss <- sort(sa2$free_Concentration_2)
  iv_bliss <- pracma::interp2(x  = x_bliss, 
                              y  = y_bliss, 
                              Z = bliss_matrix, 
                              xp = dd_2, 
                              yp = dd
  )
  
  #save values 
  dd_temp <- c(clid, cln, 
               drug, drug_2, 
               dnp, dnp_2,
               dose_id, dose_id_2, gtf$long,
               dd[1], dd[2], dd[3],
               dd_2[1], dd_2[2], dd_2[3],
               iv[4], iv[5], iv[6],
               iv[1], iv[2], iv[3],
               iv[7], iv[8], iv[9], 
               iv_hsa[1], iv_hsa[2], iv_hsa[3],
               iv_bliss[1], iv_bliss[2], iv_bliss[3]
  )
  dt_dd_temp <- data.frame(matrix(ncol = length(col_names_dd), nrow = 1, dd_temp))
  colnames(dt_dd_temp) <- col_names_dd
  #add to full dataset
  dt_dd <- rbind(dt_dd, dt_dd_temp)
}  
#convert numbers back froms string
dt_dd[, c(10:30)] <- sapply(dt_dd[, c(10:30)], as.numeric)

## Plot dose regime for each group (each a  pdf page)
cgroups_ids <- names(CGroups) 
pall <- list()
for (j in 1:length(cgroups_ids)){
  celllinename_sel <- CGroups[[cgroups_ids[j]]]$CellLineName
  p <- list()
  for  (i in 1:nrow(DDoses)) {
    DDn <- DDoses[i, c('DrugName')]
    #get single agent response of each cell line
    idx1 <- (smooth$DrugName == DDn) & 
      (smooth$free_Concentration_2==0) &
      (smooth$CellLineName %in% celllinename_sel)
    idx2 <- (smooth$DrugName_2 == DDn) & 
      (smooth$free_Concentration==0) & 
      (smooth$CellLineName %in% celllinename_sel)
    
    #extract SA
    dt_SA <- smooth[idx1,]
    dt_SA_2 <- smooth[idx2,]
    #swap Drug and Drug_2 if 
    if (nrow(dt_SA)<1)  {
      dt_SA <-dt_SA_2
      dt_SA$free_Concentration <- dt_SA$free_Concentration_2
      dt_SA$Gnumber <- dt_SA$Gnumber_2
      dt_SA$DrugName <- dt_SA$DrugName_2
    }
    #get the dose and min and max used for calculations 
    dd_ave <-DDoses[i, c('dd')]
    dd_min <- DDoses[i, c('dd_min')]
    dd_max <- DDoses[i, c('dd_max')]
    #take the mean across cell lines that had different smoothing
    dt_SA <- dplyr::select(dt_SA, c('clid', 'CellLineName',
                                    'Gnumber',
                                    'free_Concentration', 'DrugName', gtf$long))
    dt_SA<- dt_SA %>%
      group_by(clid,CellLineName,Gnumber, free_Concentration, DrugName) %>%
      summarise(across(.data[[gtf$long]], mean), .groups = 'drop')
    #plot the single agent resopnse with dose
    p[[length(p)+1]] <- ggplot(data = dt_SA, aes_string(x = 'free_Concentration', y = gtf$long, color='CellLineName'))+  
      geom_line(alpha=0.8, size = 0.7) +
      geom_hline(yintercept=0, color='royalblue1') +
      geom_vline(xintercept=dd_ave, color='firebrick2', size=1.5, linetype='solid') +
      geom_vline(xintercept=dd_min, color='firebrick2', size=1, alpha=0.5, linetype='dashed') +
      geom_vline(xintercept=dd_max, color='firebrick2', size=1, alpha=0.5, linetype='dashed') + 
      scale_x_continuous(trans = 'log10') +
      scale_y_continuous(lim=c(NA, 1.2)) +
      xlab(paste(unique(dt_SA$DrugName),'uM'))+
      scale_colour_grey(start = 0.0, end = 0.8, na.value = "red", aesthetics = "colour") +
      labs(title=paste(cgroups_ids[j],'\n',DDoses[i, c('Dose_ID')],'N=',length(unique(dt_SA$clid)))) +
      ggpubr::theme_pubr() +
      theme(text = element_text(size=6),
            axis.text = element_text(size = 6),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none"
      )
  }
  pall[[cgroups_ids[j]]] <- p
}

#save dose range plot
file_res <- sprintf('Dose_range_DA_%s_%s.pdf', gtf$short, dataset_name)
ncol <- length(pall[[cgroups_ids[1]]])
pdf(file.path(cwd, figures_dir, fig_sub_dir ,file_res), width=ncol * 4, height= 4)  
for (j in 1:length(cgroups_ids)){
  p_temp <- pall[[cgroups_ids[j]]]
  ncol <- length(pall[[cgroups_ids[j]]])
  print(grid.arrange(grobs = p_temp, ncol=ncol))
}  
dev.off() 

## Plot barplot, violin plot and waterfall for each subgroup (each a pdf page)

clines_hl <- c()
cgroups_ids <- names(CGroups) 
drug_combos <- unique(select(dt_dd, c('Gnumber', 'Gnumber_2', 'DD','DD_2','Dose_ID', 'Dose_ID_2')))
ord_cols <- c('Gnumber', 'Gnumber_2', 'DD','DD_2','Dose_ID', 'Dose_ID_2')
drug_combos <- data.table::setorderv(drug_combos, ord_cols)
pall_wf <- list() #this is for waterfall plots
pall_br <- list() #this is for barplot
pall_vi <- list() #this is for violinplot
for (j in 1:length(cgroups_ids)){
  celllinename_sel <- CGroups[[cgroups_ids[j]]]$CellLineName
  p_wf <- list()
  p_br <- list()
  p_vi <- list()
  for  (i in 1:nrow(drug_combos)) {
    #extract drug pairs
    drug <- drug_combos[['Gnumber']][i] 
    drug_2 <- drug_combos[['Gnumber_2']][i]
    DD <- drug_combos[['DD']][i] 
    DD_2 <- drug_combos[['DD_2']][i]
    Dose_ID <- drug_combos[['Dose_ID']][i] 
    Dose_ID_2 <- drug_combos[['Dose_ID_2']][i]
    #extract corresponding smooth data
    idx <- (dt_dd$Dose_ID==Dose_ID) &
      (dt_dd$Dose_ID_2==Dose_ID_2) &
      (dt_dd$CellLineName %in% celllinename_sel)
    dt_dd_sub <- dt_dd[idx,]
    
    #remove entries with NAN values
    dt_dd_sub <- dt_dd_sub %>% tidyr::drop_na(SA,SA_2,Combo)
    
    #create a ordering for SA, SA_2, Co_meas, Co_sim
    dt_dd_sub <- data.table::setorderv(dt_dd_sub, c('SA'))
    dt_dd_sub$SA_ord <- nrow(dt_dd_sub):1 
    dt_dd_sub$SA_ord <-  dt_dd_sub$SA_ord/max(dt_dd_sub$SA_ord) * 100
    dt_dd_sub <- data.table::setorderv(dt_dd_sub, c('SA_2'))
    dt_dd_sub$SA_2_ord <- nrow(dt_dd_sub):1 
    dt_dd_sub$SA_2_ord <-  dt_dd_sub$SA_2_ord/max(dt_dd_sub$SA_2_ord) * 100
    dt_dd_sub <- data.table::setorderv(dt_dd_sub, c('Combo'))
    dt_dd_sub$Combo_ord <- nrow(dt_dd_sub):1 
    dt_dd_sub$Combo_ord <-  dt_dd_sub$Combo_ord/max(dt_dd_sub$Combo_ord) * 100
    dt_dd_sub <- data.table::setorderv(dt_dd_sub, c('HSA'))
    dt_dd_sub$HSA_ord <- nrow(dt_dd_sub):1 
    dt_dd_sub$HSA_ord <-  dt_dd_sub$HSA_ord/max(dt_dd_sub$HSA_ord) * 100 
    dt_dd_sub <- data.table::setorderv(dt_dd_sub, c('Bliss'))
    dt_dd_sub$Bliss_ord <- nrow(dt_dd_sub):1 
    dt_dd_sub$Bliss_ord <-  dt_dd_sub$Bliss_ord/max(dt_dd_sub$Bliss_ord) * 100 
    
    #define colors and legend names
    lbl_x <- unique(dt_dd_sub$DrugName)
    color_x <- 'firebrick2'
    lbl_y <- unique(dt_dd_sub$DrugName_2)
    color_y <- 'royalblue2'
    lbl_combo <- 'Combo'
    color_combo <-  'darkolivegreen4'
    lbl_hsa <- 'HSA'
    color_hsa <-  'magenta3' 
    lbl_bliss <- 'Bliss'
    color_bliss <- 'orange3'
    values_col <- c(color_x, color_y, color_combo, color_hsa, color_bliss)
    names(values_col) <- c(lbl_x, lbl_y, lbl_combo, lbl_hsa, lbl_bliss)
    #define ylimits
    ymin <- NA
    ymax <- max(1.2, 
                max(na.omit(dt_dd_sub$SA)), 
                max(na.omit(dt_dd_sub$SA_2)), 
                max(na.omit(dt_dd_sub$Combo)), 
                max(na.omit(dt_dd_sub$HSA)),
                max(na.omit(dt_dd_sub$Bliss)))
    #do waterfall plot
    p_wf[[length(p_wf)+1]] <- ggplot() +
      geom_line(data=dt_dd_sub, 
                aes(x = SA_ord, y = SA, color=lbl_x), color=color_x)+ 
      geom_text(data=dt_dd_sub,
                aes(x = SA_ord, y = SA, color=lbl_x, 
                    label=ifelse(CellLineName %in% clines_hl,as.character(CellLineName),'')),
                color=color_x, size=3)+
      geom_errorbar(data=dt_dd_sub,
                    aes(x = SA_ord, ymin=SA_min, ymax=SA_max), width=.1, color=color_x, alpha=0.3) +
      geom_line(data=dt_dd_sub, 
                aes(x = SA_2_ord, y = SA_2, color=lbl_y), color=color_y) +
      geom_text(data=dt_dd_sub, 
                aes(x = SA_2_ord, y = SA_2, color=lbl_y, 
                    label=ifelse(CellLineName %in% clines_hl,as.character(CellLineName),'')),
                color=color_y, size=3) +
      geom_errorbar(data=dt_dd_sub,
                    aes(x = SA_2_ord, ymin=SA_min_2, ymax=SA_max_2), width=.1, color=color_y , alpha=0.3) +
      geom_line(data=dt_dd_sub, 
                aes(x = Combo_ord, y = Combo, color=lbl_combo)) +
      geom_text(data=dt_dd_sub, 
                aes(x = Combo_ord, y = Combo, color=lbl_combo,
                    label=ifelse(CellLineName %in% clines_hl,as.character(CellLineName),'')),
                color=color_combo, size=3)+
      geom_errorbar(data=dt_dd_sub,
                    aes(x = Combo_ord, ymin=Combo_min, ymax=Combo_max), width=.1, color=color_combo , alpha=0.3) +
      geom_line(data=dt_dd_sub, 
                aes(x = HSA_ord, y = HSA, color=lbl_hsa), color=color_hsa, linetype='dashed') +
      geom_errorbar(data=dt_dd_sub,
                    aes(x = HSA_ord, ymin=HSA_min, ymax=HSA_max), width=.1, color=color_combo , alpha=0.3) +
      geom_line(data=dt_dd_sub, 
                aes(x = Bliss_ord, y = Bliss, color=lbl_bliss), color=color_bliss, linetype='dashed') +
      geom_errorbar(data=dt_dd_sub,
                    aes(x = Bliss_ord, ymin=Bliss_min, ymax=Bliss_max), width=.1, color=color_bliss , alpha=0.3) +
      geom_hline(yintercept=0, color='gray70') +
      #scale_x_continuous(breaks=seq(0,100,25) ,labels=c('0','25','50','75','100')) +
      ylim(c(ymin,ymax)) +
      #labs(title=sprintf('%s\n%s (%guM) x %s (%guM)',cgroups_ids[j],unique(dt_dd_sub$DrugName),DD,
      #                   unique( dt_dd_sub$DrugName_2),DD_2)) +
      labs(title=sprintf('%s x\n %s',unique(dt_dd_sub$Dose_ID), unique(dt_dd_sub$Dose_ID_2)))+
      xlab(paste('% Cell Lines ','(N=',nrow(dt_dd_sub),')',sep='')) +
      ylab(unique(dt_dd_sub$Metric)) +
      scale_color_manual(values = values_col) +
      labs(color='') +
      guides(color=guide_legend(nrow=2, byrow=TRUE)) +
      ggpubr::theme_pubr() +
      theme(text = element_text(size=7),
            legend.text=element_text(size=7),
            #axis.text.x=element_blank(),
            #axis.ticks.x=element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            #legend.position = "right",
            legend.margin = margin(0, 0, 0, 0),
            plot.margin = margin(.5, .5, .5, 2, "cm")
      )
    
    #melt dataset for plotting bar format
    dt_dd_sub_melted <- reshape2::melt(dt_dd_sub, id.vars = c("clid", "CellLineName", "DrugName","DrugName_2"),
                                       measure.vars = c("SA", 
                                                        "SA_2", 
                                                        "Combo", 
                                                        "HSA",
                                                        "Bliss") 
    )
    
    #do bar plot
    values_col_bar <- c(color_x, color_y, color_combo, color_hsa, color_bliss)
    names(values_col_bar) <- c("SA", "SA_2", "Combo", "HSA","Bliss")
    p_br[[length(p_br)+1]] <- ggplot(dt_dd_sub_melted, aes(x=reorder(CellLineName, dplyr::desc(-value)), y=value, fill=variable)) + 
      geom_bar(stat="identity", position="dodge") + 
      scale_fill_manual(values = values_col_bar) +
      guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
      ylab(unique(dt_dd_sub$Metric)) +
      geom_hline(yintercept = 1.0, color='gray50') +
      geom_hline(yintercept = 0.1, color='black', linetype = "dashed") + #90%viability reduction
      #labs(title=sprintf('%s\n%s (%guM) x %s (%guM)',cgroups_ids[j],unique(dt_dd_sub_melted$DrugName),DD,
      #                   unique(dt_dd_sub_melted$DrugName_2),DD_2)) +
      labs(title=sprintf('%s x\n %s',unique(dt_dd_sub$Dose_ID), unique(dt_dd_sub$Dose_ID_2)))+
      coord_flip() +
      ggpubr::theme_pubr() +
      theme(text = element_text(size=7),
            legend.text=element_text(size=7),
            #axis.text.x=element_blank(),
            #axis.ticks.x=element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.title = element_blank(),
            #legend.position = "right",
            legend.margin = margin(0, 0, 0, 0),
            plot.margin = margin(.5, .5, .5, 2, "cm")
      )
    
    #do violin plot
    p_vi[[length(p_vi)+1]] <- ggplot(dt_dd_sub_melted, aes(x=reorder(variable, dplyr::desc(variable)), y=value, color=variable)) + 
      geom_violin(alpha=.5)+
      #geom_boxplot(alpha=.5, outlier.shape = NA) +
      geom_quasirandom(size=1.5) +
      geom_vline(xintercept=0, color='gray80') +
      #labs(title=sprintf('%s\n%s (%guM) x %s (%guM)',cgroups_ids[j],unique(dt_dd_sub_melted$DrugName),DD,
      #                   unique(dt_dd_sub_melted$DrugName_2),DD_2)) +
      labs(title=sprintf('%s x\n %s',unique(dt_dd_sub$Dose_ID), unique(dt_dd_sub$Dose_ID_2)))+
      ggpubr::theme_pubr() +
      theme(text = element_text(size=6),
            axis.text = element_text(size = 6),
            axis.text.x = element_text(angle = 0, hjust = 1),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            #legend.position = "none",
            #plot.margin = margin(0.2, 2, 0.2, 0.2, "cm")
      ) +
      coord_flip()
    
  }
  pall_wf[[cgroups_ids[j]]] <- p_wf
  pall_br[[cgroups_ids[j]]] <- p_br
  pall_vi[[cgroups_ids[j]]] <- p_vi
}

#save bar plots
file_res <- sprintf('Barplot_DA_%s_%s.pdf', gtf$short, dataset_name)
ncol <- length(pall_br[[cgroups_ids[1]]])
pdf(file.path(cwd, figures_dir, fig_sub_dir ,file_res), width=5 * ncol, height=5 )  
for (j in 1:length(cgroups_ids)){
  p_temp <- pall_br[[cgroups_ids[j]]]
  ncol <- length(p_temp)
  print(grid.arrange(grobs = p_temp, ncol=ncol))
}  
dev.off() 

#save violin plots
file_res <- sprintf('Violin_DA_%s_%s.pdf', gtf$short, dataset_name)
ncol <- length(pall_vi[[cgroups_ids[1]]])
pdf(file.path(cwd, figures_dir, fig_sub_dir ,file_res), width=3 * ncol, height=5 )  
for (j in 1:length(cgroups_ids)){
  p_temp <- pall_vi[[cgroups_ids[j]]]
  ncol <- length(p_temp)
  print(grid.arrange(grobs = p_temp, ncol=ncol))
}  
dev.off() 


#save waterfall plots
file_res <- sprintf('Waterfall_DA_%s_%s.pdf', gtf$short, dataset_name)
ncol <- length(pall_wf[[cgroups_ids[1]]])
pdf(file.path(cwd, figures_dir, fig_sub_dir ,file_res), width=5 * ncol, height=5 )  
for (j in 1:length(cgroups_ids)){
  p_temp <- pall_wf[[cgroups_ids[j]]]
  ncol <- length(p_temp)
  print(grid.arrange(grobs = p_temp, ncol=ncol))
}  
dev.off() 


## Plot individual cells heatmaps of 1) Rel Vial or GR, 2) HSA, 3) Bliss

#go through each cell line and drug combination
assay_ID <- c("SmoothMatrix", "HSAExcess", "BlissExcess")
title_ID <- c(gtf$long, "HSA Excess", "Bliss Excess")
field_ID <- c(gtf$long, "excess", "excess")
colors_fields <- list()
if (gtf$short =='RV'){
  colors_fields[[1]] <- viridis(51)
}else if (gtf$short =='GR') {
  colors_fields[[1]] <- c(rev(rocket(60)[10:60]), viridis(50))
}

colors_fields[[2]] <- colorRampPalette(c("royalblue3", "royalblue1", "grey95" , "grey95" , "firebrick1", "firebrick3"))(51)
colors_fields[[3]] <- colorRampPalette(c("royalblue3", "royalblue1", "grey95" , "grey95" , "firebrick1", "firebrick3"))(51)
mixmax_fields <- c(1,2,2) #1 for GR, Rel, 2 for HSA and Bliss, method to calculate min and max


#add free drug calculations to combo dataset
for (i in 1:length(assay_ID)){
  dt_QCS_combo[[assay_ID[i]]] <- add_free_Concentrations(dt_QCS_combo[[assay_ID[i]]], FuDrugs)
}



#for each cell line and drug combination
p_smooth <- list()
for (i in 1:nrow(dt_dd)){
  #extract Rel.vial or GR, HSA and Bliss from QCS_combo
  for (j in 1:length(assay_ID)){
    assay <- assay_ID[j]
    field <- field_ID[j]
    smooth <- dt_QCS_combo[[assay]]
    idx <- ((smooth$CellLineName==dt_dd[i, c('CellLineName')]) &
              (smooth$Gnumber==dt_dd[i, c('Gnumber')]) &
              (smooth$Gnumber_2==dt_dd[i, c('Gnumber_2')]))
    smooth <- smooth[idx, ]
    if ('normalization_type' %in% colnames(smooth)){
      idx <- (smooth$normalization_type==gtf$short)
      smooth <- smooth[idx, ]
    }
    if (nrow(smooth)>1) {
      #use free concentrations to plot
      smooth$Concentration <- smooth$free_Concentration
      smooth$Concentration_2 <- smooth$free_Concentration_2
      smooth <- createLogConcentrations(smooth)
      if (mixmax_fields[j]==1) {
        mine <- min(c(0.0, min(na.omit(smooth[,..field]))))
        maxe <- max(c(1.1, max(na.omit(smooth[,..field])))) 
        cdotdose <- 'black'
      } else if (mixmax_fields[j]==2) {
        tope <- max(abs(c(min(na.omit(smooth[,..field])), 
                          max(na.omit(smooth[,..field])))))
        mine <- min(c(-0,3,-tope))
        maxe <- max(c(0.3, tope))
        cdotdose <- 'black'
      }
      limits <- c(mine,maxe)
      #calculate width and height of each tile for geom_tile
      smooth$x <- smooth$logConcentration
      smooth$y <- smooth$logConcentration_2
      ux <- unique(smooth$x)
      uy <- unique(smooth$y)
      diffux <- diff(ux)/2
      diffuy <- diff(uy)/2
      #calculate width left and right (wl, wr) and height up and down (hu, hd)
      smooth$wl=plyr::mapvalues(smooth$x, ux, c(diffux[1],diffux))
      smooth$wr=plyr::mapvalues(smooth$x, ux, c(diffux,diffux[length(diffux)]))
      smooth$hu=plyr::mapvalues(smooth$y, uy, c(diffuy,diffuy[length(diffuy)]))
      smooth$hd=plyr::mapvalues(smooth$y, uy, c(diffuy[1],diffuy))
      #calculate xmin, xmax, ymin, ymax for geom_rect
      smooth$xmin <- (smooth$x - smooth$wl)
      smooth$xmax <- (smooth$x + smooth$wr)
      smooth$ymin <- (smooth$y - smooth$hd)
      smooth$ymax <- (smooth$y + smooth$hu)
      colors_smooth <- colors_fields[[j]]
      p_smooth[[length(p_smooth)+1]] <- ggplot()  +
        geom_rect(data=smooth, 
                  aes_string(xmin = 'xmin', xmax = 'xmax', 
                             ymin = 'ymin', ymax = 'ymax', 
                             fill = field)) +
        labs(title=sprintf('%s\n%s x %s', 
                           dt_dd[i, c('CellLineName')],
                           dt_dd[i, c('Dose_ID')],
                           dt_dd[i, c('Dose_ID_2')])) +
        xlab(paste(unique(smooth$DrugName),'uM'))+
        ylab(paste(unique(smooth$DrugName_2),'uM')) +
        scale_fill_gradientn(colours = colors_smooth, limits=limits) +
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
      
      #plot the dd concentrations on x and y axis 
      dt_dd_sub <- dt_dd[i,]
      p_smooth[[length(p_smooth)]] <- p_smooth[[length(p_smooth)]] +
        geom_point(data=dt_dd_sub,
                   aes(x = log10(DD), y = log10(DD_2)), colour = cdotdose, shape= 10, size = 5)
    }
    
  }
}

#save heatmaps with concentration plots
file_res <- sprintf('Heatmaps_DA_%s_%s.pdf', gtf$short, dataset_name)
nrow <- length(p_smooth)/3
pdf(file.path(cwd, figures_dir, fig_sub_dir ,file_res), width= 15 , height=5 * nrow  )  
print(grid.arrange(grobs = p_smooth, ncol=3))
dev.off()

#### Plot individual cells heatmaps with all the doses in the same plot #####

#filter dt_dd
dt_dd_filt <- dt_dd
#idx <- (dt_dd_filt$Dose_ID != 'SHP2i_GDC1971[SHP2i_50mg.csv]')
# idx <- ((dt_dd_filt$Dose_ID == 'SHP2i_GDC1971[SHP2i_20mg.csv]') &
#        (dt_dd_filt$Dose_ID_2 == 'KRASG12Ci_GDC6036[KRASi_400mg.csv]')) |
#        ((dt_dd_filt$Dose_ID == 'SHP2i_GDC1971[SHP2i_60mg.csv]') &
#        (dt_dd_filt$Dose_ID_2 == 'KRASG12Ci_GDC6036[KRASi_400mg.csv]')) |
#        ((dt_dd_filt$Dose_ID == 'SHP2i_GDC1971[SHP2i_20mg.csv]') &
#        (dt_dd_filt$Dose_ID_2 == 'KRASG12Ci_GDC6036[KRASi_100mg.csv]')) |
#        ((dt_dd_filt$Dose_ID == 'SHP2i_GDC1971[SHP2i_60mg.csv]') &
#         (dt_dd_filt$Dose_ID_2 == 'KRASG12Ci_GDC6036[KRASi_100mg.csv]')) 
#dt_dd_filt <- dt_dd_filt[idx,]

uclid_drugs <- unique(dplyr::select(dt_dd_filt, c('CellLineName', 'Gnumber', 'Gnumber_2')))
rownames(uclid_drugs) <- NULL
#for each cell line and drug combination
p_smooth <- list()
#for each cell line and drug combo
for (i in 1:nrow(uclid_drugs)){
  #extract Rel.vial or GR, HSA and Bliss from QCS_combo
  for (j in 1:length(assay_ID)){
    assay <- assay_ID[j]
    field <- field_ID[j]
    smooth <- dt_QCS_combo[[assay]]
    idx <- ((smooth$CellLineName==uclid_drugs[i, c('CellLineName')]) &
              (smooth$Gnumber==uclid_drugs[i, c('Gnumber')]) &
              (smooth$Gnumber_2==uclid_drugs[i, c('Gnumber_2')]))
    smooth <- smooth[idx, ]
    if ('normalization_type' %in% colnames(smooth)){
      idx <- (smooth$normalization_type==gtf$short)
      smooth <- smooth[idx, ]
    }
    
    if (nrow(smooth)>1) {
      #use free concentrations to plot
      smooth$Concentration <- smooth$free_Concentration
      smooth$Concentration_2 <- smooth$free_Concentration_2
      smooth <- createLogConcentrations(smooth)
      if (mixmax_fields[j]==1) {
        mine <- min(c(0.0, min(na.omit(smooth[,..field]))))
        maxe <- max(c(1.1, max(na.omit(smooth[,..field])))) 
        cdotdose <- 'black'
      } else if (mixmax_fields[j]==2) {
        tope <- max(abs(c(min(na.omit(smooth[,..field])), 
                          max(na.omit(smooth[,..field])))))
        mine <- min(c(-0,3,-tope))
        maxe <- max(c(0.3, tope))
        cdotdose <- 'black'
      }
      limits <- c(mine,maxe)
      #calculate width and height of each tile for geom_tile
      smooth$x <- smooth$logConcentration
      smooth$y <- smooth$logConcentration_2
      ux <- unique(smooth$x)
      uy <- unique(smooth$y)
      diffux <- diff(ux)/2
      diffuy <- diff(uy)/2
      #calculate width left and right (wl, wr) and height up and down (hu, hd)
      smooth$wl=plyr::mapvalues(smooth$x, ux, c(diffux[1],diffux))
      smooth$wr=plyr::mapvalues(smooth$x, ux, c(diffux,diffux[length(diffux)]))
      smooth$hu=plyr::mapvalues(smooth$y, uy, c(diffuy,diffuy[length(diffuy)]))
      smooth$hd=plyr::mapvalues(smooth$y, uy, c(diffuy[1],diffuy))
      #calculate xmin, xmax, ymin, ymax for geom_rect
      smooth$xmin <- (smooth$x - smooth$wl)
      smooth$xmax <- (smooth$x + smooth$wr)
      smooth$ymin <- (smooth$y - smooth$hd)
      smooth$ymax <- (smooth$y + smooth$hu)
      colors_smooth <- colors_fields[[j]]
      p_smooth[[length(p_smooth)+1]] <- ggplot()  +
        geom_rect(data=smooth, 
                  aes_string(xmin = 'xmin', xmax = 'xmax', 
                             ymin = 'ymin', ymax = 'ymax', 
                             fill = field)) +
        labs(title=sprintf('%s', 
                           uclid_drugs[i, c('CellLineName')])) +
        xlab(paste(unique(smooth$DrugName),'uM'))+
        ylab(paste(unique(smooth$DrugName_2),'uM')) +
        scale_fill_gradientn(colours = colors_smooth, limits=limits) +
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
      
      #if requested, plot IC50 isobolograms
      # if (IC50bols[j] == 1){
      #   ic_sel=c("0.5")
      #   dt_isobolograms <- dt_QCS_combo[['isobolograms']]
      #   idx <- ((dt_isobolograms$CellLineName==uclid_drugs[i, c('CellLineName')]) &
      #             (dt_isobolograms$Gnumber==uclid_drugs[i, c('Gnumber')]) &
      #             (dt_isobolograms$Gnumber_2==uclid_drugs[i, c('Gnumber_2')]))
      #   dt_isobolograms <- dt_isobolograms[idx, ]
      #   dt_isobolograms <- dt_isobolograms[dt_isobolograms$normalization_type==gtf$short]
      #   dt_isobolograms <- dt_isobolograms[dt_isobolograms$iso_level %in% ic_sel]
      #   if (nrow(dt_isobolograms)>1) {
      #     ##### ATTENTION !!!!! ### CONVERSION DONE MANUALLY BECAUSE NEED NEW CODE TO CONVERT
      #     #### POS_REF AND POS_ in ISOBOLOGRAM FOR FREE DRUG
      #     dt_isobolograms$pos_x <- dt_isobolograms$pos_x + (log10(0.06))
      #     dt_isobolograms$pos_y <- dt_isobolograms$pos_y + (log10(0.11))
      #     dt_isobolograms$pos_x_ref <- dt_isobolograms$pos_x_ref + (log10(0.06))
      #     dt_isobolograms$pos_y_ref <- dt_isobolograms$pos_y_ref + (log10(0.11))
      #     #prepare isobologram for proper plotting with x_pos and x_pos_ref
      #     dt_isob_proper <- dplyr::select(dt_isobolograms, -c('pos_x_ref','pos_y_ref'))
      #     dt_isob_proper_2 <- dplyr::select(dt_isobolograms, -c('pos_x','pos_y'))
      #     dt_isob_proper_2 <- dplyr::rename(dt_isob_proper_2, pos_x = pos_x_ref , pos_y = pos_y_ref)
      #     dt_isob_proper$iso_source <- 'measured'
      #     dt_isob_proper_2$iso_source <- 'expected'
      #     dt_isob_proper <- rbind(dt_isob_proper, dt_isob_proper_2)
      #     colors_iso <- colorRampPalette(c("red", "purple"))(length(unique(dt_isob_proper$iso_level)))
      #     dt_isob_proper$iso_level <- factor(dt_isob_proper$iso_level)
      #     dt_isob_proper$iso_source <- factor(dt_isob_proper$iso_source)
      #     #add isobologram to plot
      #     p_smooth[[length(p_smooth)]] <- p_smooth[[length(p_smooth)]] + 
      #       geom_line(data=dt_isob_proper, 
      #                 aes(x = pos_y, y = pos_x, color = iso_level, linetype= iso_source) ,size = 2) +
      #       scale_colour_manual(values = colors_iso) +
      #       scale_x_continuous(expand = c(0, 0)) +
      #       scale_y_continuous(expand = c(0, 0)) 
      #   }
      # }
      #plot all the concentrations that were tested
      dt_dd_sub <- dt_dd_filt
      idx <- ((dt_dd_sub$CellLineName==uclid_drugs[i, c('CellLineName')]) &
                (dt_dd_sub$Gnumber==uclid_drugs[i, c('Gnumber')]) &
                (dt_dd_sub$Gnumber_2==uclid_drugs[i, c('Gnumber_2')]))
      dt_dd_sub <- dt_dd_sub[idx, ]
      for (j in 1:nrow(dt_dd_sub))
        #plot the dd concentrations on x and y axis 
        p_smooth[[length(p_smooth)]] <- p_smooth[[length(p_smooth)]] +
        geom_point(data=dt_dd_sub[j,],
                   aes(x = log10(DD), y = log10(DD_2)), colour = cdotdose, shape= 10, size = 3, stroke = 1)
      
      
    }
  }
}

#save heatmaps with concentration plots
file_res <- sprintf('Heatmaps_DA_all_doses_in_one_%s_%s.pdf', gtf$short, dataset_name)
ncol <- length(assay_ID)
nrow <- length(p_smooth)/ncol
pdf(file.path(cwd, figures_dir, fig_sub_dir ,file_res), width= 5.5 * ncol , height=5 * nrow  )  
print(grid.arrange(grobs = p_smooth, ncol=length(assay_ID)))
dev.off()



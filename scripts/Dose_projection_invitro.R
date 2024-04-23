# Genentech Inc
# by Luca Gerosa and Andrew Goetz
#
#analyse the RAF/MEK screen:
#1. Study dose range of combination of RAF/MEK in cell lines

library(ggplot2)
library(gridExtra)
library(dplyr)
library(viridis)
library(gDRutils)
library(ggbeeswarm)

#set working directory
#cwd <- "/gstore/home/gerosal/projects/work/panRAFi_MEKi_combo"
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

### Prompt dataset information ###
message('Cell lines: ', length(unique(smooth$CellLineName)))
message('Drug combos: ', nrow(unique(dplyr::select(smooth, c('DrugName','DrugName_2')))))

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

plot_DDoses <- DDoses
plot_DDoses$title <- gsub("_"," ",DDoses$title)
#plot_DDoses <- arrange(plot_DDoses,plot_DDoses$loc)
d_temp <- c()
for (drug_ind in 1:length(unique(plot_DDoses$DrugName))){
  d_temp[[drug_ind]] <- ggplot(data = filter(plot_DDoses,plot_DDoses$DrugName==unique(plot_DDoses$DrugName)[drug_ind]),aes(y=title,x = dd)) + geom_col() + geom_errorbar(aes(y = title,xmin = dd_min,xmax = dd_max), width = 0.5) + xlab(expression(paste("Plasma Free Drug ", mu, "M"))) + ylab(strsplit(unique(plot_DDoses$DrugName)[drug_ind],"_")[[1]][2])+
    scale_y_discrete(limits=filter(plot_DDoses,plot_DDoses$DrugName==unique(plot_DDoses$DrugName)[drug_ind])$title)+ theme_classic()
}
strsplit(unique(plot_DDoses$DrugName)[drug_ind],"_")[[1]][2]


file_res <- 'DDose.pdf'
ncol <- length(unique(DDoses$DrugName))
pdf(file.path(cwd, figures_dir, 'Dose_projections' ,file_res), width=ncol * 4, height= 4)  
print(grid.arrange(grobs = d_temp, ncol=ncol))

dev.off() 




#create CGroups (cell line groups) list
CGroups <- list()
fgroup <- 'CellLineName'

#NRAS mutant
idx <- smooth[['NRAS_mut']]=='yes' & smooth[['BRAF_mut']]=='no'
CGroups[['NRAS_mut']] <-  unique(smooth[idx,..fgroup])

#BRAF V600E selected
idx <- smooth$CellLineName %in% c('A-375')
#idx <- smooth$CellLineName %in% c('A-375', 'WM−266−4', 'SK-MEL-28')
CGroups[['BRAF_mut']] <-  unique(smooth[idx,..fgroup])

#WT sensitive
#idx <- smooth$CellLineName %in% c('OVCA 420','HT-115')
#CGroups[['WT']] <-  unique(smooth[idx,..fgroup])

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

#flatten
field <- sprintf('%s_x', gtf$short)
groups <- c("normalization_type")
wide_cols <- c('x')
smooth <- gDRutils::flatten(smooth, groups = groups, wide_cols = wide_cols)


# Converts dose scheme into format with each row corresponding to specific PK condition
dose_conditions <- flatten_Doses(smooth,DDoses)

# Projects free drug concentrations derived from PK range (min, mean, max) to in-vitro response
projected_doses <- getProjectedPKEffects(smooth,dose_conditions,gtf)
#dt_dd <- getProjectedPKEffects_dep(smooth, CGroups, DDoses,gtf)


# re-formats output from pk projection to give PK min/max/avg different columns
projected_doses_dd <- projected_doses %>% 
  filter(PK_subcondition == 'dd') %>%
  rename(DD = free_Concentration, DD_2 = free_Concentration_2) %>% 
  dplyr::select(-PK_subcondition)
projected_doses_dd_min <- projected_doses %>% 
  filter(PK_subcondition == 'dd_min') %>% 
  rename(SA_min = SA, SA_min_2 = SA_2, DD_min = free_Concentration, DD_min_2 = free_Concentration_2, Combo_min = Combo,HSA_min = HSA,Bliss_min = Bliss) %>%
  dplyr::select(-PK_subcondition)
projected_doses_dd_max <- projected_doses %>%
  filter(PK_subcondition == 'dd_max') %>%
  rename(SA_max = SA, SA_max_2 = SA_2, DD_max = free_Concentration, DD_max_2 = free_Concentration_2, Combo_max = Combo,HSA_max = HSA,Bliss_max = Bliss) %>%
  dplyr::select(-PK_subcondition)
projected_doses_new <- merge(projected_doses_dd,projected_doses_dd_min)
dt_dd <- merge(projected_doses_new,projected_doses_dd_max)


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
      dt_SA$DrugName <- dt_SA$DrugName_2
    }
    #get the dose and min and max used for calculations 
    dd_ave <-DDoses[i, c('dd')]
    dd_min <- DDoses[i, c('dd_min')]
    dd_max <- DDoses[i, c('dd_max')]
    #take the mean across cell lines that had different smoothing
    dt_SA <- dplyr::select(dt_SA, c('CellLineName',
                                    'DrugName',
                                    'free_Concentration', gtf$long))
    dt_SA<- dt_SA %>%
      group_by(CellLineName,DrugName, free_Concentration) %>%
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
      labs(title=paste(cgroups_ids[j],'\n',DDoses[i, c('Dose_ID')],'N=',length(unique(dt_SA$CellLineName)))) +
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
#file_res <- sprintf('Dose_range_DA_%s.pdf', gtf$short)
file_res <- sprintf('Dose_range_DA_%s_FBS_%d.pdf', gtf$short,used_FBS_perc)
ncol <- length(pall[[cgroups_ids[1]]])
pdf(file.path(cwd, figures_dir, 'Dose_projections' ,file_res), width=ncol * 4, height= 4)  
for (j in 1:length(cgroups_ids)){
  p_temp <- pall[[cgroups_ids[j]]]
  ncol <- length(pall[[cgroups_ids[j]]])
  print(grid.arrange(grobs = p_temp, ncol=ncol))
}  
dev.off() 

## Plot barplot, violin plot and waterfall for each subgroup (each a pdf page)

clines_hl <- c()
cgroups_ids <- names(CGroups) 
drug_combos <- unique(dplyr::select(dt_dd, c('DrugName', 'DrugName_2', 'DD','DD_2','Dose_ID', 'Dose_ID_2')))
ord_cols <- c('DrugName', 'DrugName_2', 'DD','DD_2','Dose_ID', 'Dose_ID_2')
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
    drug <- drug_combos[['DrugName']][i] 
    drug_2 <- drug_combos[['DrugName_2']][i]
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
    dt_dd_sub_melted <- reshape2::melt(dt_dd_sub, id.vars = c("CellLineName", "DrugName","DrugName_2"),
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
#file_res <- sprintf('Barplot_DA_%s.pdf', gtf$short)
file_res <- sprintf('Barplot_DA_%s_FBS_%d.pdf', gtf$short,used_FBS_perc)
ncol <- length(pall_br[[cgroups_ids[1]]])
pdf(file.path(cwd, figures_dir, 'Dose_projections' ,file_res), width=5 * ncol, height=5 )  
for (j in 1:length(cgroups_ids)){
  p_temp <- pall_br[[cgroups_ids[j]]]
  ncol <- length(p_temp)
  print(grid.arrange(grobs = p_temp, ncol=ncol))
}  
dev.off() 

#save violin plots
#file_res <- sprintf('Violin_DA_%s.pdf', gtf$short)
file_res <- sprintf('Violin_DA_%s_%d.pdf', gtf$short,used_FBS_perc)
ncol <- length(pall_vi[[cgroups_ids[1]]])
pdf(file.path(cwd, figures_dir, 'Dose_projections',file_res), width=3 * ncol, height=5 )  
for (j in 1:length(cgroups_ids)){
  p_temp <- pall_vi[[cgroups_ids[j]]]
  ncol <- length(p_temp)
  print(grid.arrange(grobs = p_temp, ncol=ncol))
}  
dev.off() 


#save waterfall plots
#file_res <- sprintf('Waterfall_DA_%s.pdf', gtf$short)
file_res <- sprintf('Waterfall_DA_%s_%d.pdf', gtf$short,used_FBS_perc)
ncol <- length(pall_wf[[cgroups_ids[1]]])
pdf(file.path(cwd, figures_dir, 'Dose_projections' ,file_res), width=5 * ncol, height=5 )  
for (j in 1:length(cgroups_ids)){
  p_temp <- pall_wf[[cgroups_ids[j]]]
  ncol <- length(p_temp)
  print(grid.arrange(grobs = p_temp, ncol=ncol))
}  
dev.off() 

## Generate in vivo prediction heatmaps
filter(combo$BlissExcess,CellLineName == uclines_drugs$CellLineName[1],DrugName == uclines_drugs$DrugName[1],DrugName_2 == uclines_drugs$DrugName_2[1],normalization_type == gtf$short)
pall_hm <- list() # List for Bliss and Combo plots
for (i in 1:nrow(uclines_drugs)){
  cell_line_dt_dd <- dt_dd %>%
    filter(dt_dd$CellLineName == uclines_drugs$CellLineName[i])
  
  combo_matrix <- reshape2::acast(cell_line_dt_dd, DD ~ DD_2, value.var = 'Combo')

  bliss_matrix <- reshape2::acast(cell_line_dt_dd, DD ~ DD_2, value.var = 'Bliss')
  
  # formats row and column names
  DD_convert <- unique(dplyr::select(cell_line_dt_dd,Dose_ID,DD))
  rownames(combo_matrix) <- plyr::mapvalues(rownames(combo_matrix), from=DD_convert$DD, to=DD_convert$Dose_ID)
  rownames(combo_matrix) <- gsub("\\]","",gsub("_"," ",unlist(strsplit(rownames(combo_matrix),"\\["))[c(FALSE, TRUE)]))
  DD_2_convert <- unique(dplyr::select(cell_line_dt_dd,Dose_ID_2,DD_2))
  colnames(combo_matrix) <- plyr::mapvalues(colnames(combo_matrix), from=DD_2_convert$DD_2, to=DD_2_convert$Dose_ID_2)
  colnames(combo_matrix) <- gsub("\\]","",gsub("_"," ",unlist(strsplit(colnames(combo_matrix),"\\["))[c(FALSE, TRUE)]))
  
  # formats row and column names
  DD_convert <- unique(dplyr::select(cell_line_dt_dd,Dose_ID,DD))
  rownames(bliss_matrix) <- plyr::mapvalues(rownames(bliss_matrix), from=DD_convert$DD, to=DD_convert$Dose_ID)
  rownames(bliss_matrix) <- gsub("\\]","",gsub("_"," ",unlist(strsplit(rownames(bliss_matrix),"\\["))[c(FALSE, TRUE)]))
  DD_2_convert <- unique(dplyr::select(cell_line_dt_dd,Dose_ID_2,DD_2))
  colnames(bliss_matrix) <- plyr::mapvalues(colnames(bliss_matrix), from=DD_2_convert$DD_2, to=DD_2_convert$Dose_ID_2)
  colnames(bliss_matrix) <- gsub("\\]","",gsub("_"," ",unlist(strsplit(colnames(bliss_matrix),"\\["))[c(FALSE, TRUE)]))
  
  # Gets drug names for x and y titles
  #drug_1 <- strsplit(cell_line_dt_dd$DrugName[[1]],"_")[[1]][2]
  drug_1 <- cell_line_dt_dd$DrugName[[1]]
  #drug_2 <- strsplit(cell_line_dt_dd$DrugName_2[[1]],"_")[[1]][2]
  drug_2 <- cell_line_dt_dd$DrugName_2[[1]]
  
  # Generates bliss heatmaps
  bliss_long <- reshape2::melt(bliss_matrix-combo_matrix)
  
  bliss_smooth <- filter(combo$BlissExcess,CellLineName == uclines_drugs$CellLineName[i],DrugName == uclines_drugs$DrugName[i],DrugName_2 == uclines_drugs$DrugName_2[i])

  flatten
  field <- sprintf('%s_x', gtf$short)
  groups <- c("normalization_type")
  wide_cols <- c("excess")
  bliss_smooth <- gDRutils::flatten(bliss_smooth, groups = groups, wide_cols = wide_cols)
  
  tope <- max(abs(c(min(bliss_smooth[[paste(gtf$short,"_excess",sep="")]],na.rm = TRUE), 
                          max(bliss_smooth[[paste(gtf$short,"_excess",sep="")]],na.rm = TRUE))))
  max_val <- max(tope,.5)
  min_val <- min(-tope,-.5)
  
  pall_hm[[length(pall_hm)+1]] <- ggplot(bliss_long, aes(x = Var1, y = Var2))+
    geom_tile(aes(fill=value))+xlab(drug_1) +
    geom_text(aes(label = round(value, 2)))+ylab(drug_2) +
    ggtitle(sprintf("%s\nBliss Excess",uclines_drugs$CellLineName[i])) + 
    scale_fill_gradientn( colours=c("royalblue3", "royalblue1", "grey95","grey95", "firebrick1", "firebrick3"),
                          limits=range(min_val,max_val),
                          name = "Bliss Excess"
    )+theme(plot.title = element_text(hjust = 0.5))+ 
    theme(panel.background = element_blank())+ 
    theme(axis.line = element_line(colour = "black"))
  #scale_fill_viridis(discrete=FALSE,limits=range(-.5,1))
  # Generates combo heatmaps
  combo_long <- reshape2::melt(combo_matrix)
  
  smooth_ph <- filter(smooth,CellLineName == uclines_drugs$CellLineName[i],DrugName == uclines_drugs$DrugName[i],DrugName_2 == uclines_drugs$DrugName_2[i])

  
  if (gtf$long=='GRvalue'){
    colors <-  c(rev(rocket(60)[10:60]), viridis(50))
    max_val <- max(1.1,smooth_ph[[gtf$long]],na.rm = TRUE)
    min_val <- min(-0.5,smooth_ph[[gtf$long]],na.rm = TRUE) 
    #max_val <- max(1.1,combo_matrix,na.rm = TRUE)
    #min_val <- min(-0.5,combo_matrix,na.rm = TRUE) 
  }else{
    colors <- viridis(50)
    max_val <- max(1.1,smooth_ph[[gtf$long]],na.rm = TRUE)
    min_val <- min(0,smooth_ph[[gtf$long]],na.rm = TRUE) 
  }
  limits = range(min_val,max_val)

  pall_hm[[length(pall_hm)+1]] <- ggplot(combo_long, aes(x = Var1, y = Var2))+
    geom_tile(aes(fill=value))+xlab(drug_1) +
    geom_text(aes(label = round(value, 2)))+ylab(drug_2) +
    ggtitle(sprintf("%s\nCombo",uclines_drugs$CellLineName[i])) + 
    scale_fill_gradientn(colours=colors,
                          limits=limits,
                          name = "Combo"
    )+theme(plot.title = element_text(hjust = 0.5))+ 
    theme(panel.background = element_blank())+ 
    theme(axis.line = element_line(colour = "black"))
}

#scale_fill_gradient2(
#  low = "royalblue2", 
#  mid = "white", 
#  high = "firebrick2", 
#  midpoint = midpoint_val,
#  limits=range(min_val,1),
#  name = "Combo"
#)



#save in dose projection heatmaps
file_res <- sprintf('Dose_Projection_Heatmaps_%s_FBS_%d.pdf', gtf$short,used_FBS_perc)
nrow <- length(pall_hm) %/% 2
pdf(file.path(cwd, figures_dir, 'Dose_projections' ,file_res), width=6 * 2, height=5*nrow )  
print(grid.arrange(grobs = pall_hm, ncol=2))
dev.off() 

## Plot individual cells heatmaps of 1) Rel Vial or GR, 2) HSA, 3) Bliss

#go through each cell line and drug combination
assay_ID <- c("SmoothMatrix", "HSAExcess", "BlissExcess")
title_ID <- c(gtf$long, "HSA Excess", "Bliss Excess")
field_ID <- c(gtf$long, paste(gtf$short,"_excess", sep=""), paste(gtf$short,"_excess", sep=""))
wide_ID <- c("x","excess","excess")
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
  combo[[assay_ID[i]]] <- add_free_Concentrations(combo[[assay_ID[i]]], FuDrugs)
}



#for each cell line and drug combination
p_matv <- list()
for (i in 1:nrow(dt_dd)){
  #extract Rel.vial or GR, HSA and Bliss from QCS_combo
  for (j in 1:length(assay_ID)){
    assay <- assay_ID[j]
    field <- field_ID[j]
    matv <- combo[[assay]]
    idx <- ((matv$CellLineName==dt_dd[i, c('CellLineName')]) &
              (matv$DrugName==dt_dd[i, c('DrugName')]) &
              (matv$DrugName_2==dt_dd[i, c('DrugName_2')]))
    matv <- matv[idx, ]
    if ('normalization_type' %in% colnames(matv)){
      groups <- c("normalization_type")
      wide_cols <- c(wide_ID[j])
      matv <- gDRutils::flatten(matv, groups = groups, wide_cols = wide_cols)
    }
    if (nrow(matv)>1) {
      #use free concentrations to plot
      matv$Concentration <- matv$free_Concentration
      matv$Concentration_2 <- matv$free_Concentration_2
      matv <- createLogConcentrations(matv)
      if (mixmax_fields[j]==1) {
        if (gtf$short =='RV'){
          mine <- min(c(0.0, min(na.omit(matv[,..field]))))
          maxe <- max(c(1.1, max(na.omit(matv[,..field])))) 
        }else if (gtf$short =='GR') {
          mine <- min(c(-0.5, min(na.omit(matv[,..field]))))
          maxe <- max(c(1.1, max(na.omit(matv[,..field])))) 
        }
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
        labs(title=sprintf('%s\n%s x %s', 
                           dt_dd[i, c('CellLineName')],
                           dt_dd[i, c('Dose_ID')],
                           dt_dd[i, c('Dose_ID_2')])) +
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
      
      #plot the dd concentrations on x and y axis 
      dt_dd_sub <- dt_dd[i,]
      p_matv[[length(p_matv)]] <- p_matv[[length(p_matv)]] +
        geom_point(data=dt_dd_sub,
                   aes(x = log10(DD), y = log10(DD_2)), colour = cdotdose, shape= 10, size = 5)
    }
    
  }
}

#save heatmaps with concentration plots
#file_res <- sprintf('Heatmaps_DA_%s.pdf', gtf$short)
file_res <- sprintf('Heatmaps_DA_%s_FBS_%d.pdf', gtf$short,used_FBS_perc)
nrow <- length(p_matv)/3
pdf(file.path(cwd, figures_dir, 'Dose_projections' ,file_res), width= 15 , height=5 * nrow  )  
print(grid.arrange(grobs = p_matv, ncol=3))
dev.off()

#### Plot individual cells heatmaps with all the doses in the same plot #####

#filter dt_dd
dt_dd_filt <- dt_dd
uclines_drugs <- unique(dplyr::select(dt_dd_filt, c('CellLineName', 'DrugName', 'DrugName_2')))
rownames(uclines_drugs) <- NULL
#for each cell line and drug combination
p_matv <- list()
#for each cell line and drug combo
for (i in 1:nrow(uclines_drugs)){
  #extract Rel.vial or GR, HSA and Bliss from QCS_combo
  for (j in 1:length(assay_ID)){
    assay <- assay_ID[j]
    field <- field_ID[j]
    matv <- combo[[assay]]
    idx <- ((matv$CellLineName==uclines_drugs[i, c('CellLineName')]) &
              (matv$DrugName==uclines_drugs[i, c('DrugName')]) &
              (matv$DrugName_2==uclines_drugs[i, c('DrugName_2')]))
    matv <- matv[idx, ]
    if ('normalization_type' %in% colnames(matv)){
      groups <- c("normalization_type")
      wide_cols <- c(wide_ID[j])
      matv <- gDRutils::flatten(matv, groups = groups, wide_cols = wide_cols)
    }
    
    if (nrow(matv)>1) {
      #use free concentrations to plot
      matv$Concentration <- matv$free_Concentration
      matv$Concentration_2 <- matv$free_Concentration_2
      matv <- createLogConcentrations(matv)
      if (mixmax_fields[j]==1) {
        if (gtf$short =='RV'){
          mine <- min(c(0.0, min(na.omit(matv[,..field]))))
          maxe <- max(c(1.1, max(na.omit(matv[,..field])))) 
        }else if (gtf$short =='GR') {
          mine <- min(c(-0.5, min(na.omit(matv[,..field]))))
          maxe <- max(c(1.1, max(na.omit(matv[,..field])))) 
        }
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
      
      #plot all the concentrations that were tested
      dt_dd_sub <- dt_dd_filt
      idx <- ((dt_dd_sub$CellLineName==uclines_drugs[i, c('CellLineName')]) &
                (dt_dd_sub$DrugName==uclines_drugs[i, c('DrugName')]) &
                (dt_dd_sub$DrugName_2==uclines_drugs[i, c('DrugName_2')]))
      dt_dd_sub <- dt_dd_sub[idx, ]
      for (j in 1:nrow(dt_dd_sub))
        #plot the dd concentrations on x and y axis 
        p_matv[[length(p_matv)]] <- p_matv[[length(p_matv)]] +
        geom_point(data=dt_dd_sub[j,],
                   aes(x = log10(DD), y = log10(DD_2)), colour = cdotdose, shape= 10, size = 3, stroke = 1)
      
      
    }
  }
}

#save heatmaps with concentration plots
#file_res <- sprintf('Heatmaps_DA_all_doses_in_one_%s.pdf', gtf$short)
file_res <- sprintf('Heatmaps_DA_all_doses_in_one_%s_FBS_%d.pdf', gtf$short,used_FBS_perc)
ncol <- length(assay_ID)
nrow <- length(p_matv)/ncol
pdf(file.path(cwd, figures_dir, 'Dose_projections' ,file_res), width= 5.5 * ncol , height=5 * nrow  )  
print(grid.arrange(grobs = p_matv, ncol=length(assay_ID)))
dev.off()








#### Plot individual cells heatmaps #####

#filter dt_dd
dt_dd_filt <- dt_dd
uclines_drugs <- unique(dplyr::select(dt_dd_filt, c('CellLineName', 'DrugName', 'DrugName_2')))
rownames(uclines_drugs) <- NULL
#for each cell line and drug combination
p_matv <- list()
#for each cell line and drug combo
for (i in 1:nrow(uclines_drugs)){
  #extract Rel.vial or GR, HSA and Bliss from QCS_combo
  for (j in 1:length(assay_ID)){
    assay <- assay_ID[j]
    field <- field_ID[j]
    matv <- combo[[assay]]
    idx <- ((matv$CellLineName==uclines_drugs[i, c('CellLineName')]) &
              (matv$DrugName==uclines_drugs[i, c('DrugName')]) &
              (matv$DrugName_2==uclines_drugs[i, c('DrugName_2')]))
    matv <- matv[idx, ]
    if ('normalization_type' %in% colnames(matv)){
      groups <- c("normalization_type")
      wide_cols <- c(wide_ID[j])
      matv <- gDRutils::flatten(matv, groups = groups, wide_cols = wide_cols)
    }
    
    if (nrow(matv)>1) {
      #use free concentrations to plot
      matv$Concentration <- matv$free_Concentration
      matv$Concentration_2 <- matv$free_Concentration_2
      matv <- createLogConcentrations(matv)
      if (mixmax_fields[j]==1) {
        if (gtf$short =='RV'){
          mine <- min(c(0.0, min(na.omit(matv[,..field]))))
          maxe <- max(c(1.1, max(na.omit(matv[,..field])))) 
        }else if (gtf$short =='GR') {
          mine <- min(c(-0.5, min(na.omit(matv[,..field]))))
          maxe <- max(c(1.1, max(na.omit(matv[,..field])))) 
        }
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
    }
  }
}

#save heatmaps with concentration plots
#file_res <- sprintf('Heatmaps_DA_all_doses_in_one_%s.pdf', gtf$short)
file_res <- sprintf('Heatmaps_Blank_%s_FBS_%d.pdf', gtf$short,used_FBS_perc)
ncol <- length(assay_ID)
nrow <- length(p_matv)/ncol
pdf(file.path(cwd, figures_dir, 'Dose_projections' ,file_res), width= 5.5 * ncol , height=5 * nrow  )  
print(grid.arrange(grobs = p_matv, ncol=length(assay_ID)))
dev.off()

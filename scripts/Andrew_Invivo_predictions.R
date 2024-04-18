# Genentech Inc
# by Luca Gerosa and Andrew Goetz
#
#analyse the RAF/MEK screen:
#1. Predict xenograft tumor growth of IPC-298 cells from in-vitro responses

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
source(paste(cwd,'/scripts/Andrew_utility.R', sep=''))

#set directory names
data_dir = 'data'
results_dir = 'results'
figures_dir='figures'

#set directory names
data_dir = 'data'
results_dir = 'results'
figures_dir='figures'
clines_dir ='clines_annotations'
drugs_dir ='drug_annotations'
dda_dir = 'DDA'

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
choise_gtf <- 0
if (choise_gtf == 0){
  gtf$long <- 'RelativeViability'
  gtf$short <- 'RV'
} else{
  gtf$long <- 'GRvalue' 
  gtf$short <- 'GR'
}

#load the DDoses definitions

DDoses <- data.frame(read.csv(file.path(cwd, data_dir, 'Invivo_predictions', 'DDoses_mouse.csv')))
#generate a unique identifier
DDoses$Dose_ID <- sprintf('%s[%s]', DDoses$DrugName, DDoses$title)
#create CGroups (cell line groups) list
CGroups <- list()
fgroup <- 'CellLineName'

#IPC298
CGroups[['IPC298']] <-  'IPC-298'

#convert Concentrations to Concentrations_free using FuDrugs
used_FBS_perc=10
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


#### Plot individual cells heatmaps with all the doses in the same plot #####

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
mixmax_fields <- c(1,2,2) #1 for GR, Rel, 2 for HSA and Bliss, method to calculate min and max


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
        mine <- min(c(-0,3,-tope))
        maxe <- max(c(0.3, tope))
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
file_res <- sprintf('Andrew_Heatmaps_DA_all_doses_in_one_%s_FBS_%d.pdf', gtf$short,used_FBS_perc)
ncol <- length(assay_ID)
nrow <- length(p_matv)/ncol
pdf(file.path(cwd, figures_dir, 'Invivo_predictions' ,file_res), width= 5.5 * ncol , height=5 * nrow)
print(grid.arrange(grobs = p_matv, ncol=length(assay_ID)))
dev.off()



### calculate and plot time course response with different baseline growth rate

#create DDoses (doses) dataframe
CDynam = data.frame(Sim_ID=character(), 
                    CellLineName= character(), 
                    Dose_ID = character(), 
                    Dose_ID_2 = character(), 
                    growth_rate= numeric())
cln_sel=c("IPC-298")
gr_baseline=data.frame(CellLineName= cln_sel, 
                       Growth_rate= log(2)/c(18)
)
for (i in 1:length(cln_sel)) {
  #dv_time=unique(dt_QCS_all_sa$ReferenceDivisionTime[dt_QCS_all_sa$CellLineName==cln_sel[i]])/24
  idx = gr_baseline$CellLineName==cln_sel[i]
  dv_time=gr_baseline[idx, c('Growth_rate')]
  #doses single agents
  CDynam[nrow(CDynam)+1,] <- c(paste('1',cln_sel[i],sep='_'), cln_sel[i], 'panRAFi_Belvarafenib[15mgkg_QD]', '', dv_time)
  CDynam[nrow(CDynam)+1,] <- c(paste('1',cln_sel[i],sep='_'), cln_sel[i], 'panRAFi_Belvarafenib[30mgkg_QD]', '', dv_time)
  CDynam[nrow(CDynam)+1,] <- c(paste('1',cln_sel[i],sep='_'), cln_sel[i], '', 'MEKi_Cobimetinib[5mgkg_QD]', dv_time)
  CDynam[nrow(CDynam)+1,] <- c(paste('1',cln_sel[i],sep='_'), cln_sel[i], 'panRAFi_Belvarafenib[15mgkg_QD]', 'MEKi_Cobimetinib[5mgkg_QD]', dv_time)
  CDynam[nrow(CDynam)+1,] <- c(paste('1',cln_sel[i],sep='_'), cln_sel[i], 'panRAFi_Belvarafenib[30mgkg_QD]', 'MEKi_Cobimetinib[5mgkg_QD]', dv_time)
}
CDynam[, c(5)] <- sapply(CDynam[, c(5)], as.numeric)


#set color and order for this particular plotting 
conds_u <- c("Untreated", 
             " x\nMEKi_Cobimetinib[5mgkg_QD]", 
             "panRAFi_Belvarafenib[15mgkg_QD] x\n", 
             "panRAFi_Belvarafenib[30mgkg_QD] x\n", 
             "panRAFi_Belvarafenib[15mgkg_QD] x\nMEKi_Cobimetinib[5mgkg_QD]",
             "panRAFi_Belvarafenib[30mgkg_QD] x\nMEKi_Cobimetinib[5mgkg_QD]")
levels_bar <- rev(conds_u)

colors_conds <- c(rgb(175/255, 206/255, 225/255),
                  rgb(60/255, 120/255, 176/255),
                  rgb(187/255, 221/255, 147/255),
                  rgb(85/255, 158/255, 64/255),
                  rgb(237/255, 159/255, 155/255),
                  rgb(208/255, 53/255, 43/255))
names(colors_conds) <- conds_u

#simulate plot of time course growth
t_tot<-35
Sim_ID <- unique(CDynam$Sim_ID)
strt_vol <- 200 
p_bars <- list()
p_lines <- list()
for (i in 1:length(Sim_ID)){
  idx <- CDynam$Sim_ID==Sim_ID[i]
  CDynam_sub <- CDynam[idx,]
  gr_base <- unique(CDynam_sub$growth_rate)
  cln <-  unique(CDynam_sub$CellLineName)
  #get growth rate for each cell line and drug combo
  growth_rate <- c(1.0) #initialize untreated control to 1.0
  for (j in 1:nrow(CDynam_sub)){
    Dose_ID <- CDynam_sub[j, c('Dose_ID')]
    Dose_ID_2 <- CDynam_sub[j, c('Dose_ID_2')]
    #find the drug dose in dt_dd
    if(Dose_ID_2==''){
      idx <- (dt_dd$CellLineName==cln) &
        (dt_dd$Dose_ID==Dose_ID)
      dt_dd_sub <- dt_dd[idx,]
      growth_rate[j+1] <- mean(dt_dd_sub$SA)
    }else if(Dose_ID==''){
      idx <- (dt_dd$CellLineName==cln) &
        (dt_dd$Dose_ID_2==Dose_ID_2)
      dt_dd_sub <- dt_dd[idx,]
      growth_rate[j+1] <- mean(dt_dd_sub$SA_2)
    }else{
      idx <- (dt_dd$CellLineName==cln) &
        (dt_dd$Dose_ID==Dose_ID) &
        (dt_dd$Dose_ID_2==Dose_ID_2) 
      dt_dd_sub <- dt_dd[idx,]
      growth_rate[j+1] <-dt_dd_sub$Combo
    }
  }
  #calculate growth rate lines
  growth_rate <- gr_base * growth_rate
  df_growth_rate <- data.frame(growth_rate)
  g_names <- c('Untreated', paste(CDynam_sub$Dose_ID,CDynam_sub$Dose_ID_2, sep=' x\n'))
  #rownames(gr_drug) <- g_names
  df_growth_rate$condition <- g_names
  #plot
  p_bars[[length(p_bars)+1]] <- ggplot(data=df_growth_rate, 
                                       aes(y=growth_rate, x=factor(condition, levels=levels_bar), fill=condition)) +
    geom_bar(stat="identity") +
    labs(title=cln)+
    ggpubr::theme_pubr() +
    guides(fill=guide_legend(nrow=4, byrow=TRUE)) +
    theme(text = element_text(size=10),
          axis.text = element_text(size = 10),
          #axis.text.x = element_text(angle = 0, hjust = 1),
          #axis.ticks.x=element_blank(),
          #axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank()
          #legend.position = "none",
          #plot.margin = margin(0.2, 2, 0.2, 0.2, "cm")
    ) + 
    scale_fill_manual(values=colors_conds) +
    coord_flip()
  
  #calculate volume changes
  sims_celllinename <- list()
  sims_cond <- list()
  sims_time <- list()
  sims_vol <- list()
  for (z in 1:nrow(df_growth_rate)){
    sims_celllinename[length(sims_celllinename)+1]<-cln 
    sims_celllinename[length(sims_celllinename)+1]<- cln 
    sims_cond[length(sims_cond)+1]<-df_growth_rate[z,c('condition')]
    sims_cond[length(sims_cond)+1]<-df_growth_rate[z,c('condition')]
    sims_time[length(sims_time)+1]<- 0
    sims_time[length(sims_time)+1]<- t_tot
    sims_vol[length(sims_vol)+1]<- strt_vol
    gr <- df_growth_rate[z,c('growth_rate')]
    sims_vol[[length(sims_vol)+1]]<- strt_vol * exp(gr*t_tot)
  }
  dt_sim <-  data.frame(CellLineName= unlist(sims_celllinename),
                        Condition= unlist(sims_cond),
                        Time = unlist(sims_time),
                        Volume= unlist(sims_vol))
  
  #plot time course simulations 
  p_lines[[length(p_lines)+1]] <- ggplot(data=dt_sim, 
                                         aes(x=Time, y=Volume, color=Condition)) +    
    geom_line(size=2) +
    labs(title=cln)+
    scale_y_continuous(trans='log2') +
    #ggpubr::theme_pubr() +
    guides(color=guide_legend(nrow=4, byrow=TRUE)) +
    theme(
      panel.background = element_rect(fill = "white", colour = "white",
                                      linewidth = 2, linetype = "solid"),
      panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                      colour = "gray85"),
      legend.position = "top",
      #panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
      #                                colour = "gray85")
    ) +
    scale_color_manual(values=colors_conds) 
    # theme(text = element_text(size=10),
    #       axis.text = element_text(size = 10),
    #       axis.text.x = element_text(angle = 0, hjust = 1),
    #       plot.title = element_text(hjust = 0.5),
    #       #panel.grid.major = element_blank(),
    #       #panel.grid.minor = element_blank(),
    #       legend.title = element_blank()
    #       #legend.position = "none",
    #       #plot.margin = margin(0.2, 2, 0.2, 0.2, "cm")
    # )
}

#plot growth rates
#file_res <- sprintf('Sim_Growth_Rates_DA_%s.pdf', gtf$short)
file_res <- sprintf('Andrew_Sim_Growth_Rates_DA_%s_FBS_%d.pdf', gtf$short,used_FBS_perc)
pdf(file.path(cwd, figures_dir, 'Invivo_predictions' ,file_res), width=10, height=10)  
print(grid.arrange(grobs = p_bars, ncol=length(cln_sel)))
dev.off() 

#plot simulations
#file_res <- sprintf('Sim_Time_Course_DA_%s.pdf', gtf$short)
file_res <- sprintf('Andrew_Sim_Time_Course_DA_%s_FBS_%d.pdf', gtf$short,used_FBS_perc)
pdf(file.path(cwd, figures_dir, 'Invivo_predictions' ,file_res), width=8, height=10)  
print(grid.arrange(grobs = p_lines, ncol=length(cln_sel)))
dev.off() 


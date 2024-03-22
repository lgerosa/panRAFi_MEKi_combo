# Genentech Inc
# by Luca Gerosa and Andrew Goetz
#
#analyse the RAF/MEK screen:
#1. Plot SA response metrics (IC50, AUC, Emax)
#2. Plot drug combo response metrics (HSA, Bliss)
#3. Plot individual cell line heatmaps 

library(ggplot2)
library(gridExtra)
library(dplyr)
library(viridis)
library(gDRutils)

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

sa <- readRDS(file.path(cwd, data_dir, 'Drug_screen', 'Drug_screen_Belva_Cobi_sa_metrics.RDS'))
combo <- readRDS(file.path(cwd, data_dir, 'Drug_screen', 'Drug_screen_Belva_Cobi_combo_metrics.RDS'))
#prepare annotations of mutations
anno <- readRDS(file.path(cwd, data_dir, 'Drug_screen', 'Drug_screen_cell_line_annotations.RDS'))
#prepare annotations for plotting
anno_hm <- anno
rownames(anno_hm) <-anno_hm$CellLineName
anno_hm <- dplyr::select(anno_hm, -c('CellLineName'))
#define colors
col_anno <- list()
col_anno$NRAS_mut <- c(yes='black', no='white')
col_anno$BRAF_mut <- c(yes='black', no='white')

### Prompt dataset information ###
message('Cell lines: ', length(unique(sa$CellLineName)))
message('Drugs: ', length(unique(sa$DrugName)))
message('Drug combos: ', nrow(unique(dplyr::select(combo$RawTreated, c('DrugName','DrugName_2')))))


### GENERATE SINGLE-AGENT METRIC FIGURES ###

#decide which metrics to plot
gtf <- list()
choise_gtf <- 1
if (choise_gtf == 0){
  gtf$long <- 'RelativeViability'
  gtf$short <- 'RV'
} else{
  gtf$long <- 'GRvalue' 
  gtf$short <- 'GR'
}


#quantity used as quantification of drug effects 
aqm <- c('%s_gDR_xc50', '%s_gDR_x_max', '%s_gDR_x_mean')
for (i in 1:length(aqm)){
  aqm[i] <- sprintf(aqm[i], gtf$short)
}
#name plot for quantification of drug effects
aqmlab <-c('IC50_uM', 'E_max', 'AUC')
#scale of y axis
qmfunc <- c(log10, identity, identity)
p <- list()
for (i in 1:length(aqm)) {
  
  #order according to mutations 
  sa_mut <- merge(sa, anno , by='CellLineName')
  sa_mut <- data.table::setorderv(sa_mut, c('BRAF_mut', "NRAS_mut",  aqm[i]))
  
  #extract metrics to plot in matrix format
  sa_mat <- data.table::dcast(sa_mut, 
                                  factor(CellLineName, levels =unique(CellLineName)) ~  factor(DrugNamePlot), 
                                  value.var = aqm[i])
  #order drugs
  sa_mat <- sa_mat[, c('CellLineName', 'MEKi_Cobimetinib', 'panRAFi_Belvarafenib', 'BRAFi_Vemurafenib'), drop=FALSE ]
  
  #make it a data.frame or pheatmap messes up annotations
  sa_mat <- as.data.frame(sa_mat)
  rownames(sa_mat)<- sa_mat$CellLineName
  sa_mat <- dplyr::select(sa_mat, -c('CellLineName'))
  #remove columns and rows that are all NAN
  sa_mat <- sa_mat[, !apply(is.na(sa_mat), 2, all)]
  sa_mat <- sa_mat[!apply(is.na(sa_mat), 1, all), ]
  #convert to right scale
  sa_mat[] <- lapply(sa_mat, function(x) qmfunc[[i]](x))
  
  #flip 
  t_sa_mat <- data.table::transpose(sa_mat)
  # get row and colnames in order
  colnames(t_sa_mat) <- rownames(sa_mat)
  rownames(t_sa_mat) <- colnames(sa_mat)
  
  #plot heatmap
  breaks <- seq(from=min(na.omit(sa_mat)), to=1.0, length.out=50)
  hmcol <- rev(colorRampPalette(c("firebrick2", "white" ))(51))
  p[[i]]<- pheatmap::pheatmap(t_sa_mat, scale="none", display_numbers = TRUE, fontsize_number=4, number_color = "black",
                              color=rev(hmcol), breaks= breaks, na_col = "white", angle_col = 90, fontsize=6,
                              treeheight_row = 30, treeheight_col = 30,
                              cluster_rows = F,
                              cluster_cols = F,
                              annotation_col= anno_hm,
                              annotation_colors= col_anno,
                              main= paste(aqmlab[i], gtf$long)
  )
}  

#save heatmaps
file_res <- sprintf('sa_metrics_heatmaps_%s.pdf', gtf$short)
pdf(file.path(cwd, figures_dir, 'Drug_screen' ,file_res), width=10, height= 2)  
for (i in 1:length(p)){
  print(p[[i]])
  if (i < length(p))
    grid::grid.newpage()
}  
dev.off() 

### GENERATE COMBO METRIC FIGURES ###

#quantity used as quantification of drug effects 
aqm <- c('HSAScore', 'BlissScore')
p<- list()
#for each aggregate quantity
for (i in 1:length(aqm)) {
  combo_met <- combo[[aqm[i]]]
  
  #filter for the wanted growth measure
  field <- sprintf('%s_x', gtf$short)
  groups <- c("normalization_type")
  wide_cols <- c('x')
  combo_met <- gDRutils::flatten(combo_met, groups = groups, wide_cols = wide_cols)
  
  #combo_met <- combo_met[combo_met$normalization_type == gtf$short, ]
  
  #order according to mutations 
  combo_met <- merge(combo_met, anno , by='CellLineName')
  combo_met <- data.table::setorderv(combo_met, c('BRAF_mut', "NRAS_mut", gtf$long))
   
  #create field with both drugs
  combo_met$Drugs_combo_name <- paste(combo_met$DrugNamePlot, combo_met$DrugNamePlot_2, sep=' x ')
  #dcast to matrix format the HSA score
  Combo_heatmap <- data.table::dcast(combo_met, 
                                     factor(CellLineName, levels =unique(CellLineName)) ~  factor(Drugs_combo_name), 
                                     value.var = gtf$long)
  Combo_heatmap <- as.data.frame(Combo_heatmap)
  # make clids rownames
  rownames(Combo_heatmap)<- Combo_heatmap$CellLineName
  Combo_heatmap <- select(Combo_heatmap, select = -c('CellLineName') )
  # remove cell lines or drugs with all NA
  Combo_heatmap <- Combo_heatmap[, !apply(is.na(Combo_heatmap), 2, all)]
  Combo_heatmap <- Combo_heatmap[!apply(is.na(Combo_heatmap), 1, all), ]
  
  #define colors and breaks
  breaks <- seq(from=-0.7, to=0.7, length.out=50)
  hmcol <- rev(colorRampPalette(c("royalblue2", "royalblue1", "grey95" , "grey95" , "firebrick1", "firebrick2"))(51))
  
  #transpose
  t_Combo_heatmap <- data.table::transpose(Combo_heatmap)
  # get row and colnames in order
  colnames(t_Combo_heatmap) <- rownames(Combo_heatmap)
  rownames(t_Combo_heatmap) <- colnames(Combo_heatmap)
  #heatmap 
  p[[i]] <- pheatmap::pheatmap(t_Combo_heatmap, scale="none", display_numbers = TRUE, fontsize_number=10, number_color = "black",
                               na_col = "white", angle_col = 45, fontsize=10, breaks=breaks, color=rev(hmcol),
                               treeheight_row = 30, treeheight_col = 30,
                               #labels_col=labels_row$CellLineName,
                               annotation_col=anno_hm, 
                               annotation_colors=col_anno,
                               show_rownames=T,
                               main= aqm[i],
                               cluster_rows = FALSE,
                               cluster_cols = FALSE
  )
  #dev.off()
}

#save heatmaps
file_res <- sprintf('combo_metrics_heatmap_%s.pdf', gtf$short)
pdf(file.path(cwd, figures_dir, 'Drug_screen' ,file_res), width=25, height= 3.5)  
for (i in 1:length(p)){
  print(p[[i]])
  if (i < length(p))
    grid::grid.newpage()
}  
dev.off() 


### GENERATE CELL LINE SPECIFIC DOSE RESPONSE HEATMAPS  ###

#plot selected lines
cell_lines_sel <- c( 
                    #'A-375', # BRAF V600E
                    'WM−266−4', # BRAF V600E
                    'SK-MEL-28',# BRAF V600E
                    'MEL-JUSO', #, # NRAS Q61
                    'SK-MEL-2' # NRAS Q61
)

#plot all lines
cell_lines_sel <- unique(combo[['Averaged']]$CellLineName)

#define assays to plot
#plot_ID <- c("SmoothMatrix", "SmoothMatrix_swapped", "isobolograms", "HSAExcess", "BlissExcess", "isobologramsIC50") 
plot_ID <- c("SmoothMatrix", "BlissExcess") 

nplot_ID <- length(plot_ID)
#select combinations of cell lines and drugs to plot (sort by cell line)
clines_drugs <- unique(select( combo[['Averaged']], c('CellLineName', 'DrugName', 'DrugName_2')))
clines_drugs <- clines_drugs[ clines_drugs$DrugName=='panRAFi_Belvarafenib' &
                              clines_drugs$DrugName_2=='MEKi_Cobimetinib', ]

ord_cols <- c('DrugName', 'DrugName_2','CellLineName')
clines_drugs <- data.table::setorderv(clines_drugs, ord_cols)
#filter cell lines
idx <- clines_drugs$CellLineName %in% cell_lines_sel
clines_drugs <- clines_drugs[idx,]

#extract cell lines and drugs to plot
clines <- cell_lines_sel
drug <- clines_drugs$DrugName
drug_2 <- clines_drugs$DrugName_2



#plots
p <- list()
#for each cell line and drug combo combination
for (i in 1:length(clines)){
  #extract combo data for cell lines 
  combo_sub <- list()
  for (j in 1:length(plot_ID)){
    aID <- plot_ID[[j]]
    idx <- (combo[[aID]]$CellLineName==clines[i]) &
      (combo[[aID]]$DrugName==drug[i]) &
      (combo[[aID]]$DrugName_2==drug_2[i])
    combo_sub[[aID]] <- combo[[aID]][idx]
  }  
  

  #use a coloring schema and min value different between GR and RV
  if (gtf$long=='GRvalue'){
    colors <-  c(rev(rocket(60)[10:60]), viridis(50))
    minF <- -1.0
    maxF <- 1.1
  }else{
    colors <- viridis(50)
    minF <- 0.0
    maxF <- 1.1
  }
  
  #plot matrix
  field <- sprintf('%s', gtf$long)
  dt_smooth <- combo_sub[['SmoothMatrix']]
  groups <- c("normalization_type")
  wide_cols <- c('x')
  dt_smooth <- gDRutils::flatten(dt_smooth, groups = groups, wide_cols = wide_cols)
  mine <- min(c(minF, min(na.omit(dt_smooth[,..field]))))
  maxe <- max(c(maxF, max(na.omit(dt_smooth[,..field])))) 
  limits <- c(mine,maxe)
  p[[length(p)+1]]<- plotHeatMapCombo(dt_smooth, field, limits, colors)

  #plot Bliss Excess
  field <- sprintf('%s_excess', gtf$short)
  dt_smooth <- combo_sub[['BlissExcess']]
  groups <- c("normalization_type")
  wide_cols <- c('excess')
  dt_smooth <- gDRutils::flatten(dt_smooth, groups = groups, wide_cols = wide_cols)
  colors <- colorRampPalette(c("royalblue3", "royalblue1", "grey95" , "grey95" , "firebrick1", "firebrick3"))(51)
  mine <- min(c(-0.7, min(na.omit(dt_smooth[,..field]))))
  maxe <- max(c(0.7, max(na.omit(dt_smooth[,..field])))) 
  limits <- c(mine,maxe)
  p[[length(p)+1]] <- plotHeatMapCombo(dt_smooth, field, limits, colors) 
}

#plot figure
ncol <- length(plot_ID)
nrow <- length(clines)
file_res <- sprintf('combo_averaged_bliss_heatmap_%s.pdf', gtf$short)
pdf(file.path(cwd, figures_dir, 'Drug_screen' ,file_res), width= 4.5 * ncol, height= 3 * nrow)  
print(grid.arrange(grobs = p, ncol=ncol))
dev.off()


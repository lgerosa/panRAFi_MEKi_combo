#analyse the RAF/MEK screen:
#1. Plot SA response metrics with annotated mutations and lineage of cell lines
#2. Plot drug combo response metrics with annotated mutations and lineage of cell lines
#3. Study dose range of combination of RAF/MEK

library(ggplot2)
library(gridExtra)
library(dplyr)

#set working directory
cwd <- "/gstore/home/gerosal/projects/work/panRAFi_MEKi_combo"
setwd(cwd)

#load utility functions
#source(paste(cwd,'/scripts/utility_scripts.R', sep=''))

#set directory names
data_dir = 'data'
results_dir = 'results'
figures_dir='figures'


### Load drug screen data ###
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

### Generate figures ###

#decide which metrics to plot
gtf <- list()
gtf$long <- 'RelativeViability' # 'GRvalue')
gtf$short <- 'RV' #'GR'

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
  #add annotation mutations for sorting
  sa_mut <- merge(sa, anno , by='CellLineName')
  sa_mut <- data.table::setorderv(sa_mut, c('BRAF_mut', "NRAS_mut",  aqm[i]))

  #extract metrics to plot in matrix format
  sa_mat <- data.table::dcast(sa_mut, 
                                  factor(CellLineName, levels =unique(CellLineName)) ~  factor(DrugNamePlot), 
                                  value.var = aqm[i])
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
  #dev.off()
}  

#save heatmaps
file_res <- sprintf('sa_metrics_Heatmaps_%s.pdf', gtf$short)
pdf(file.path(cwd, figures_dir, 'Drug_screen' ,file_res), width=10, height= 3)  
for (i in 1:length(p)){
  print(p[[i]])
  if (i < length(p))
    grid::grid.newpage()
}  
dev.off() 



### THE CODE NEEDS TO BE UDPDATE FROM HERE ONE BY LUCA

#### LOAD ANNOTATION INFORMATION FOR CELL LINES
#add mutation information 

#add zygosity information 
CL_ann_sa_zygo <- unique(dplyr::select(dt_QCS_sa_metrics, c('CellLineName', 'clid')))
CL_ann_sa_zygo <- merge(CL_ann_sa_zygo, CL_anno$CL_zygo, by= 'clid', all.x=T)
CL_ann_sa_zygo[is.na(CL_ann_sa_zygo)] <- "WT"
#addd tissue information
CL_ann_sa_lin <- unique(dplyr::select(dt_QCS_sa_metrics, c('CellLineName', 'clid', 'Tissue')))
#merge together
CL_ann_sa <- unique(dplyr::select(dt_QCS_sa_metrics, c('clid', 'CellLineName')))
CL_ann_sa <- merge(CL_ann_sa, CL_ann_sa_mut, by= c('CellLineName', 'clid'), all.x=T)
#CL_ann_sa <- merge(CL_ann_sa, CL_ann_sa_zygo, by= 'clid', all.x=T)
CL_ann_sa <- merge(CL_ann_sa, CL_ann_sa_lin, by= c('CellLineName', 'clid'), all.x=T)
#create color list for annotation
#col_ann_sa <- c(CL_anno$col_anno_CL_mut, CL_anno$col_anno_CL_zygo)
col_ann_sa <- c(CL_anno$col_anno_CL_mut)


#### LOAD DRUG ANNNOTATIONS 
drug_anno_manual <- data.frame(read.csv(file.path(cwd, data_dir, drugs_dir, 'drug_annotation_manual.csv')))
drug_anno_manual$DrugNamePlot <- paste(drug_anno_manual$drug_type, drug_anno_manual$DrugName,  sep='_')
drug_mapping_manual <- mapDrugNamePlot(dt_QCS_sa_averaged, drug_anno_manual)
#sa agent
dt_QCS_all_sa_averaged <- addDrugNamePlot(dt_QCS_sa_averaged, drug_mapping_manual)
dt_QCS_sa_metrics <- addDrugNamePlot(dt_QCS_sa_metrics, drug_mapping_manual)

#combos
assay_ID <- names(dt_QCS_combo)
for (i in 1:length(assay_ID)){
  drug_mapping_manual <- mapDrugNamePlot(dt_QCS_combo[[assay_ID[i]]], drug_anno_manual)
  dt_QCS_combo[[assay_ID[i]]] <- addDrugNamePlot(dt_QCS_combo[[assay_ID[i]]], drug_mapping_manual)
}

#save tables with the right columns for plotting for paper
remove_columns_sa <- c('clid', 'parental_identifier', 'Gnumber', 'rId', 'cId')
sa_metrics <- dplyr::select(dt_QCS_sa_metrics, -all_of(remove_columns_sa))
saveRDS(sa_metrics, file.path(cwd, 'results_paper', 'Drug_screen_Belva_Cobi_sa_metrics.RDS'))

remove_columns_combo <- c('clid', 'parental_identifier', 'Gnumber', 'Gnumber_2', 'rId', 'cId', 'record_id', 'barcode')
combo_metrics <- dt_QCS_combo
assay_ID <- names(combo_metrics)
for (i in 1:length(assay_ID)){
  combo_metrics[[assay_ID[i]]] <- dplyr::select(combo_metrics[[assay_ID[i]]], -all_of(remove_columns_sa))
}
saveRDS(combo_metrics, file.path(cwd, 'results_paper', 'Drug_screen_Belva_Cobi_combo_metrics.RDS'))


#####  DECIDE METRICS TO USE 

gr_mtr <- 'RelativeViability' #'GRvalue' # # # #
gtf <- getGrowthTypeFormat(gr_mtr)

#replace undefned or overestimaed IC50 with max concentration used
xc50_str <- sprintf('%s_gDR_xc50', gtf$short)
inf_xc50 <-is.infinite(dt_QCS_sa_metrics[[xc50_str]])
dt_QCS_sa_metrics[inf_xc50, c(xc50_str)] <- 10^dt_QCS_sa_metrics[inf_xc50, c('maxlog10Concentration')]
over_xc50 <- dt_QCS_sa_metrics[[xc50_str]] > 10^dt_QCS_sa_metrics$maxlog10Concentration
dt_QCS_sa_metrics[over_xc50, c(xc50_str)] <- 10^dt_QCS_sa_metrics[over_xc50, c('maxlog10Concentration')]


#### 0) PLOT CENSUS OF MUTATIONS AND LINEAGES

# CL_ann_sa_melt <- reshape2::melt(CL_ann_sa, id.vars = c('clid','CellLineName'))
# p <- ggplot( data=CL_ann_sa_melt, aes(x=variable, fill=value)) + geom_bar(color='black') #+  geom_text(size = 3, position = position_stack(vjust = 0.5))
# file_res <- sprintf('Census_Mutations_Lineages_%s.pdf', dataset_name)
# pdf(file.path(cwd, figures_dir, dataset_name ,file_res), width=6, height=6)  
# print(p)
# dev.off()

##### 1) PLOT SA HEAMTAPS  AUC, IC50, Emax ##### 

#save cell annotations
CL_ann_sa_ready$CellLineName <- rownames(CL_ann_sa_ready)
saveRDS(CL_ann_sa_ready, file.path(cwd, 'results_paper', 'Drug_screen_cell_line_annotations.RDS'))

#prepare cell annotations
CL_ann_sa_ready <- CL_ann_sa
CL_ann_sa_ready <- dplyr::select(CL_ann_sa_ready, -c("CellLineName","clid"))
rownames(CL_ann_sa_ready)<- CL_ann_sa$CellLineName

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
  #extract metrics to plot in matrix format
  SA_mat_ori <- data.table::dcast(dt_QCS_sa_metrics, 
                                      factor(CellLineName, levels =unique(CellLineName)) ~  factor(DrugNamePlot), 
                                      value.var = aqm[i])
  #make it a data.frame or pheatmap messes up annotations
  SA_mat_cvd <- as.data.frame(SA_mat_ori)
  rownames(SA_mat_cvd)<- SA_mat_cvd$CellLineName
  SA_mat_cvd <- dplyr::select(SA_mat_cvd, -c('CellLineName'))
  #remove columns and rows that are all NAN
  SA_mat_cvd <- SA_mat_cvd[, !apply(is.na(SA_mat_cvd), 2, all)]
  SA_mat_cvd <- SA_mat_cvd[!apply(is.na(SA_mat_cvd), 1, all), ]
  #convert to right scale
  SA_mat_cvd[] <- lapply(SA_mat_cvd, function(x) qmfunc[[i]](x))
  
  #flip 
  t_SA_mat_cvd <- data.table::transpose(SA_mat_cvd)
  # get row and colnames in order
  colnames(t_SA_mat_cvd) <- rownames(SA_mat_cvd)
  rownames(t_SA_mat_cvd) <- colnames(SA_mat_cvd)
  
  #plot heatmap
  breaks <- seq(from=min(na.omit(SA_mat_cvd)), to=1.0, length.out=50)
  hmcol <- rev(colorRampPalette(c("firebrick2", "white" ))(51))
  #file_res <- sprintf('%s_%s_%s_Heatmap_SA.pdf',aqmlab[i], gtf$short, dataset_name)
  #pdf(file.path(cwd, figures_dir, dataset_name ,file_res), width=15, height=4)  
  p[[i]]<- pheatmap::pheatmap(t_SA_mat_cvd, scale="none", display_numbers = TRUE, fontsize_number=8, number_color = "black",
                             color=rev(hmcol), breaks= breaks, na_col = "white", angle_col = 90, fontsize=6,
                             treeheight_row = 30, treeheight_col = 30,
                             annotation_col= CL_ann_sa_ready,
                             annotation_colors=col_ann_sa,
                             main= paste(aqmlab[i], dataset_name)
  )
  #dev.off()
}  

#save heatmaps
file_res <- sprintf('SA_Metrics_Heatmaps_%s_%s.pdf', gtf$short, dataset_name)
pdf(file.path(cwd, figures_dir, 'gCSI' ,file_res), width=15, height= 3.5)  
for (i in 1:length(p)){
  print(p[[i]])
  if (i < length(p))
     grid::grid.newpage()
}  
dev.off() 


##### 1.1) PLOT BAR PLOTS AUC, IC50, Emax for each drug  with mutation status as color##### 


# drug_sel <- unique(dt_QCS_sa_metrics$Gnumber)
# p<- list()
# for (i in 1:length(drug_sel)){
#   #select the drug
#   idx <- dt_QCS_sa_metrics$Gnumber== drug_sel[i] 
#   dt_QCS_sa_metrics_sel <- dt_QCS_sa_metrics[idx,]
#   #add the mutational information
#   CL_ann_sa_sel <- dplyr::select(CL_ann_sa, c('clid', 'BRAF_mut', 'KRAS_mut', 'NRAS_mut','RASRAF_wt'))
#   CL_ann_sa_sel$mutations <- apply(CL_ann_sa_sel,1,function(x) paste(names(x)[which(x=="yes")], collapse='_'))
#   idx <- CL_ann_sa_sel$mutations %in% c('BRAF_mut_KRAS_mut','BRAF_mut_NRAS_mut')
#   CL_ann_sa_sel[idx, c('mutations')] <- 'BRAF_RAS_mut'
#   
#   dt_QCS_sa_metrics_sel <- merge(dt_QCS_sa_metrics_sel, CL_ann_sa_sel, by='clid', x.all=TRUE)
#   col_sel <- colnames(CL_ann_sa_sel)
#   idx <- col_sel[CL_ann_sa_sel[, ..col_sel] == 'yes']
#   
#   #plot 
#   p[[i]] <- ggplot(dt_QCS_sa_metrics_sel, aes(x = reorder(clid, RV_gDR_xc50 ), y = RV_gDR_xc50, fill=mutations)) +
#     geom_bar(stat = "identity") +
#     xlab(paste('Cell Lines ','(N=',nrow(dt_QCS_sa_metrics_sel),')',sep='')) +
#     ylab(unique(paste(dt_QCS_sa_metrics_sel$DrugNamePlot, 'IC50 uM'))) +
#     scale_fill_manual("", values = c('KRAS_mut' = "orange", 'NRAS_mut' = "firebrick", 'RASRAF_wt' = "gray", 'BRAF_mut' = 'blue', 'BRAF_mut_RAS_mut' = 'black')) +
#     scale_y_continuous(trans='log10') +
#     ggpubr::theme_pubr() +
#     theme(text = element_text(size=15),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank(),
#           plot.title = element_text(hjust = 0.5),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           #legend.position = "none"
#           plot.margin = margin(.5, .5, .5, 2, "cm")
#     )
# }
# 
# 
# #save dose range plot
# file_res <- sprintf('SA_Metrics_BarPlots_%s_%s.pdf', gtf$short, dataset_name)
# pdf(file.path(cwd, figures_dir, dataset_name ,file_res), width=8, height= 5)  
# for (i in 1:length(p)){
#   print(p[[i]])
# }  
# dev.off() 


##### 2) PLOT Combo HSA and Bliss score ####

#prepare cell annotations
CL_ann_combo_ready <- CL_ann_sa
CL_ann_combo_ready <- dplyr::select(CL_ann_combo_ready, -c("CellLineName","clid"))
rownames(CL_ann_combo_ready)<- CL_ann_sa$CellLineName

#quantity used as quantification of drug effects 
aqm <- c('HSAScore', 'BlissScore')
p<- list()
#for each aggregate quantity
for (i in 1:length(aqm)) {
  Combo_QCS <- dt_QCS_combo[[aqm[i]]]
  
  #flatten columns
  field <- sprintf('%s_x', gtf$short)
  groups <- c("normalization_type")
  wide_cols <- c('x')
  Combo_QCS <- gDRutils::flatten(Combo_QCS, groups = groups, wide_cols = wide_cols)
  

  #create field with both drugs
  Combo_QCS$Drugs_combo_name <- paste(Combo_QCS$DrugNamePlot, Combo_QCS$DrugNamePlot_2, sep=' x ')
  #dcast to matrix format the HSA score
  Combo_heatmap <- data.table::dcast(Combo_QCS, 
                                     factor(CellLineName, levels =unique(CellLineName)) ~  factor(Drugs_combo_name), 
                                     value.var = gtf$long)
  Combo_heatmap <- as.data.frame(Combo_heatmap)
  # make clids rownames
  rownames(Combo_heatmap)<- Combo_heatmap$CellLineName
  Combo_heatmap <- select(Combo_heatmap, select = -c('CellLineName') )
  # remove cell lines or drugs with all NA
  #Combo_heatmap <- Combo_heatmap[, !apply(is.na(Combo_heatmap), 2, all)]
  #Combo_heatmap <- Combo_heatmap[!apply(is.na(Combo_heatmap), 1, all), ]
  
  #create row annotations and labels
  #labels_row <- unique(select(Combo_QCS, c('clid','CellLineName')))
  #rownames(labels_row)<- labels_row$clid
  #labels_row <- select(labels_row, select = -c('clid') )
  breaks <- seq(from=-0.7, to=0.7, length.out=50)
  hmcol <- rev(colorRampPalette(c("royalblue2", "royalblue1", "grey95" , "grey95" , "firebrick1", "firebrick2"))(51))
  
  #transpose
  t_Combo_heatmap <- data.table::transpose(Combo_heatmap)
  # get row and colnames in order
  colnames(t_Combo_heatmap) <- rownames(Combo_heatmap)
  rownames(t_Combo_heatmap) <- colnames(Combo_heatmap)
  #heatmap 
  #file_res <- sprintf('%s_%s_%s.pdf', aqm[i], gtf$short, dataset_name)
  #pdf(file.path(cwd, figures_dir , dataset_name, file_res), width=25, height=4)  
  p[[i]] <- pheatmap::pheatmap(t_Combo_heatmap, scale="none", display_numbers = TRUE, fontsize_number=10, number_color = "black",
                     na_col = "white", angle_col = 45, fontsize=10, breaks=breaks, color=rev(hmcol),
                     treeheight_row = 30, treeheight_col = 30,
                     #labels_col=labels_row$CellLineName,
                     annotation_col=CL_ann_combo_ready, 
                     annotation_colors=col_ann_sa,
                     show_rownames=T,
                     main= aqm[i],
                     cluster_rows = FALSE
  )
  #dev.off()
}

#save heatmaps
file_res <- sprintf('Combo_Metrics_Heatmaps_%s_%s.pdf', gtf$short, dataset_name)
pdf(file.path(cwd, figures_dir,  'gCSI' ,file_res), width=25, height= 3.5)  
for (i in 1:length(p)){
  print(p[[i]])
  if (i < length(p))
    grid::grid.newpage()
}  
dev.off() 



# ##### 2.1) PLOT HSA and Bliss score for each drug with mutation status as color ##### 
# 
# aqm <- c('BlissScore')
# drug_combo_sel <- unique(dt_QCS_combo[[aqm]]$rId)
# p<- list()
# for (i in 1:length(drug_combo_sel)){
#   #select the drug
#   idx <- dt_QCS_combo[[aqm]]$rId== drug_combo_sel[i] 
#   dt_QCS_combo_sel <- dt_QCS_combo[[aqm]][idx,]
#   #add the mutational information
#   CL_ann_sa_sel <- dplyr::select(CL_ann_sa, c('clid', 'BRAF_mut', 'KRAS_mut', 'NRAS_mut','RASRAF_wt'))
#   CL_ann_sa_sel$mutations <- apply(CL_ann_sa_sel,1,function(x) paste(names(x)[which(x=="yes")], collapse='_'))
#   idx <- CL_ann_sa_sel$mutations %in% c('BRAF_mut_KRAS_mut','BRAF_mut_NRAS_mut')
#   CL_ann_sa_sel[idx, c('mutations')] <- 'BRAF_RAS_mut'
#   
#   dt_QCS_combo_sel <- merge(dt_QCS_combo_sel, CL_ann_sa_sel, by='clid', x.all=TRUE)
# 
#   #plot 
#   p[[i]] <- ggplot(dt_QCS_combo_sel, aes(x = reorder(clid, RV ), y = RV, fill=mutations)) +
#     geom_bar(stat = "identity") +
#     xlab(paste('Cell Lines ','(N=',nrow(dt_QCS_combo_sel),')',sep='')) +
#     ylab(unique(paste(dt_QCS_combo_sel$DrugNamePlot, dt_QCS_combo_sel$DrugNamePlot_2, 'Bliss'))) +
#     scale_fill_manual("", values = c('KRAS_mut' = "orange", 'NRAS_mut' = "firebrick", 'RASRAF_wt' = "gray", 'BRAF_mut' = 'blue', 'BRAF_mut_RAS_mut' = 'black')) +
#     #scale_y_continuous(trans='log10') +
#     ggpubr::theme_pubr() +
#     theme(text = element_text(size=15),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank(),
#           plot.title = element_text(hjust = 0.5),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           #legend.position = "none"
#           plot.margin = margin(.5, .5, .5, 2, "cm")
#     )
# }
# 
# #save dose range plot
# file_res <- sprintf('Combo_Metrics_BarPlots_%s_%s.pdf', gtf$short, dataset_name)
# pdf(file.path(cwd, figures_dir, dataset_name ,file_res), width=8, height= 5)  
# for (i in 1:length(p)){
#   print(p[[i]])
# }  
# dev.off() 

##### 3) PLOT ISOBOLOGRAM, HSA AND BLISS SCORE FOR SELECTED LINES, DRUGS and COMBOS

cell_lines_sel <- c('A-375', # BRAF V600E
                    #'SK-MEL-28',# BRAF V600E
                    'MEL-JUSO' #, # NRAS Q61
                    #'SK-MEL-2' # NRAS Q61
                    )
#combo_sel <- list(c('panRAFi_Belvarafenib','MEKi_Cobimetinib')) $TO DO

#define assays to plot
#assay_ID <- c("SmoothMatrix", "SmoothMatrix_swapped", "isobolograms", "HSAExcess", "BlissExcess", "isobologramsIC50") 
assay_ID <- c("isobolograms", "BlissExcess") 

nassay_ID <- length(assay_ID)
#extract values
QCS_values <-  dt_QCS_combo[['Averaged']]
QCS_metrics <- dt_QCS_combo[['Metrics']]
#select combinations of cell lines and drugs to plot (sort by cell line)
clids_drugs <- unique(select(QCS_values, c('clid','CellLineName','Gnumber', 'Gnumber_2', 'DrugNamePlot', 'DrugNamePlot_2')))
ord_cols <- c('Gnumber', 'Gnumber_2','clid')
clids_drugs <- data.table::setorderv(clids_drugs, ord_cols)
#filter cell lines
idx <- clids_drugs$CellLineName %in% cell_lines_sel
clids_drugs <- clids_drugs[idx,]
#filter drug combos (TO DO)

p <- plotComboAgentFit(SE_drug_combo, clids_drugs, gtf, assay_ID, drug_mapping_manual)
#generate pdf
file_res <- sprintf('%s_%s_selected_combos_individual_cell_lines.pdf', gtf$short, dataset_name)
pdffilepath <- file.path(cwd, figures_dir , dataset_name, file_res)
#rearrange the plot with cell lines on the columns, 6 different plots on rows, 1 page per drug combo
nclines <-length(unique(clids_drugs$clid))
ndrugcombos <-nrow(unique(dplyr::select(clids_drugs, c('Gnumber','Gnumber_2')))) 
nplots <- nassay_ID
#generate new indexes
new_indexes <- c()
for (di in 1:ndrugcombos) {
  for (pi in 1:nplots) {
    for (ci in 1:nclines) {
      new_indexes <- c(new_indexes, 
                       (((ci-1) * (nplots)+1) + (pi-1)) +
                         ((di-1) * (nclines * nplots)))
      
    }
  }
}
#reorder plots to have lines on columns and one page per drug combo
p_ord <- p[new_indexes]
#print plot
plotWithMultiplePages(p_ord, pdffilepath, nplots, nclines, 3.5, 5 )


## PLOT SINGLE AGENT IC50 VS SYNERGY BLISS SCORE FOR BELVA and COBI ACROSS MUTATIONAL STATUS

#IC50 Belva
IC50_Belva <- dt_QCS_sa_metrics[ dt_QCS_sa_metrics$DrugNamePlot == 'panRAFi_Belvarafenib', 
                                 c('clid', 'Tissue','CellLineName', 'RV_gDR_xc50')]
IC50_Belva <- dplyr::rename(IC50_Belva, panRAFi_Belvarafenib_IC50_uM = RV_gDR_xc50 )

#IC50 Cobi
IC50_Cobi <- dt_QCS_sa_metrics[ dt_QCS_sa_metrics$DrugNamePlot == 'MEKi_Cobimetinib', 
                                c('clid', 'RV_gDR_xc50')]
IC50_Cobi <- dplyr::rename(IC50_Cobi,  MEKi_Cobimetinib_IC50_uM = RV_gDR_xc50)
#merge
IC50_Bliss <- merge(IC50_Belva, IC50_Cobi, by='clid')
#add bliss data
Bliss_QCS <- dt_QCS_combo$BlissScore
#filter to keep only IC50 for Belva and Cobi
idx <- (Bliss_QCS$DrugNamePlot == 'panRAFi_Belvarafenib') & 
       (Bliss_QCS$DrugNamePlot_2 == 'MEKi_Cobimetinib')
Bliss_QCS <- Bliss_QCS[idx, c('clid', 'RV')]  
Bliss_QCS <- dplyr::rename(Bliss_QCS,  Bliss_Score = RV)
#merge Bliss
IC50_Bliss <- merge(IC50_Bliss, Bliss_QCS, by='clid')
#add NRAS Q61 status
mut_status <- dplyr::select(CL_ann_sa_mut, c('clid','NRAS_61Q/X' , 'BRAF_600V/EK', 'BRAF_class2', 'BRAF_class3'))
colnames(mut_status) <- c('clid','NRAS_Q61X', 'BRAF_V600', 'BRAF_class2', 'BRAF_class3')
mut_status$mut_status <- 'other'
idx <- mut_status$BRAF_V600 == 'yes'
mut_status[idx, c('mut_status')] <- 'BRAF_V600'
idx <- mut_status$BRAF_class2 == 'yes'
mut_status[idx, c('mut_status')] <- 'BRAF_class_2_3'
idx <- mut_status$BRAF_class3 == 'yes'
mut_status[idx, c('mut_status')] <- 'BRAF_class_2_3'
idx <- mut_status$NRAS_Q61X == 'yes'
mut_status[idx, c('mut_status')] <- 'NRAS_Q61'

IC50_Bliss <- merge(IC50_Bliss, mut_status, by='clid')
IC50_Bliss$NRAS_Q61X <- as.factor(IC50_Bliss$NRAS_Q61X)
IC50_Bliss$clabel <- IC50_Bliss$CellLineName
idx <- IC50_Bliss$NRAS_Q61X == 'no' | IC50_Bliss$Tissue != 'Skin'
IC50_Bliss[idx, c('clabel')] <- NaN

#order the cell lines with Q61 melanoma at bottom to be on top in the plot
IC50_Bliss <- IC50_Bliss %>% arrange(clabel) %>% arrange(desc(row_number()))

co_Cobi_lw <- 0.03
co_Cobi_hg <- 0.15

#plot IC50 and Bliss
p <- ggplot(IC50_Bliss, aes(panRAFi_Belvarafenib_IC50_uM, MEKi_Cobimetinib_IC50_uM, 
                       size=Bliss_Score, fill=mut_status, label=clabel)) +
  geom_point(shape = 21, colour = "black") +
  scale_size(range = c(0, 10)) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  geom_text(size=4, nudge_x=0, nudge_y=0.1)+ 
  geom_hline(yintercept=co_Cobi_lw, color="gray50", 
             linetype="dashed", size=0.2)+
  geom_hline(yintercept=co_Cobi_hg, color="gray50", 
             linetype="dashed", size=0.2)+
  ggpubr::theme_pubr() 
file_res <- sprintf('IC50_vs_Comb_Belva_Cobi_%s_%s_.pdf', gtf$short, dataset_name)
pdf(file.path(cwd, figures_dir, dataset_name ,file_res), width=8, height=6)  
print(p)
dev.off()

#extract the names of cell lines in the same sector as NRAS Q61 melanoma and plot their 
#isobolograms

idx <- (IC50_Bliss$MEKi_Cobimetinib_IC50_uM > co_Cobi_lw) &
       (IC50_Bliss$MEKi_Cobimetinib_IC50_uM < co_Cobi_hg)
IC50_Bliss[idx,]
 
  
  


### 4) PERFORM COMBINATION ANALYSIS  AND PLOTS ##### 

#set name of the plotting directory
fig_sub_dir <- dataset_name

##MODIFICATION NOTE:: load manual cell lines quantification from Fran and annotate them with drugs names

used_FBS_perc_select = 5
use_Fran_data = FALSE

if (use_Fran_data){

#load combo
SE_name <- '/gstore/home/gerosal/projects/collab/CellBasedDataImports/gDRinternalData/2022_10_05_384w_belva_cobi_synergy/MAE_20221005_belva_cobi/MAE_20221005_belva_cobi_SE_drug_combo.RDS'
SE_drug_combo <- readRDS(SE_name)
dt_QCS_combo <- list()
assay_ID <- c("RawTreated", "Normalized", "Averaged", "SmoothMatrix", "isobolograms", "HSAExcess", "BlissExcess", "HSAScore", "BlissScore")
for (i in 1:length(assay_ID)){
  dt_QCS_combo[[assay_ID[i]]] <- gDRutils::convert_se_assay_to_dt(SE_drug_combo, assay_name = assay_ID[i])
}
#combos
assay_ID <- names(dt_QCS_combo)
for (i in 1:length(assay_ID)){
  drug_mapping_manual <- mapDrugNamePlot(dt_QCS_combo[[assay_ID[i]]], drug_anno_manual)
  dt_QCS_combo[[assay_ID[i]]] <- addDrugNamePlot(dt_QCS_combo[[assay_ID[i]]], drug_mapping_manual)
}

#### LOAD ANNOTATION INFORMATION FOR CELL LINES
CL_anno <- readRDS(file.path(cwd, data_dir, clines_dir, 'CL_anno.RDS'))
#add mutation information
CL_ann_sa_mut <- unique(dplyr::select(dt_QCS_combo$Averaged, c('clid')))
CL_ann_sa_mut <- merge(CL_ann_sa_mut, CL_anno$CL_mut, by= 'clid', all.x=T)
CL_ann_sa_mut[is.na(CL_ann_sa_mut)] <- "no"
#add lineage information
CL_ann_sa_zygo <- unique(dplyr::select(dt_QCS_combo$Averaged, c('clid')))
CL_ann_sa_zygo <- merge(CL_ann_sa_zygo, CL_anno$CL_zygo, by= 'clid', all.x=T)
CL_ann_sa_zygo[is.na(CL_ann_sa_zygo)] <- "WT"
#addd tissue information
CL_ann_sa_lin <- unique(dplyr::select(dt_QCS_combo$Averaged, c('clid','Tissue')))
#merge together
CL_ann_sa <- unique(dplyr::select(dt_QCS_combo$Averaged, c('clid', 'CellLineName')))
CL_ann_sa <- merge(CL_ann_sa, CL_ann_sa_mut, by= 'clid', all.x=T)
CL_ann_sa <- merge(CL_ann_sa, CL_ann_sa_zygo, by= 'clid', all.x=T)
CL_ann_sa <- merge(CL_ann_sa, CL_ann_sa_lin, by= 'clid', all.x=T)
#create color list for annotation
col_ann_sa <- c(CL_anno$col_anno_CL_mut, CL_anno$col_anno_CL_zygo)

fig_sub_dir <- "QCS-34656_DS000012937_Fran_data"

used_FBS_perc_select = 10

}

### END MODIFICATION TO USE FRANS DATA


#### Dose-Dependent Analysis (DDA):

## Inputs:
## DDose -> dataframe that defines SA dose regime for treatments
## CGroups -> list that define cell line groups to be assessed
## DAsses -> dataframe that specifies for which SA or Combo we should asses viability values

#load the DDoses definitions
DDoses <- data.frame(read.csv(file.path(cwd, data_dir, dda_dir, 'DDoses.csv')))
#generate a unique identifier
#DDoses$Dose_ID <- sprintf('%s(%guM[%g-%g])', DDoses$DrugNamePlot, DDoses$dd, DDoses$dd_min, DDoses$dd_max)
DDoses$Dose_ID <- sprintf('%s[%s]', DDoses$DrugNamePlot, DDoses$title)

#create CGroups (cell line groups) list
CGroups <- list()
fgroup <- 'CellLineName'
CL_ann_group <- CL_ann_sa

#NRAS Q61
idx <- ((CL_ann_sa[['NRAS_61Q/X']]=='yes') &
       (CL_ann_sa[['Tissue']]=='Skin')) 
CGroups[['Skin_NRAS_Q61X']] <-  unique(CL_ann_sa[idx,..fgroup])

#BRAF V600E selected
idx <- CL_ann_sa$CellLineName %in% c('A-375','SK-MEL-1', 'SK-MEL-28')
CGroups[['Skin_BRAF_V600E']] <-  unique(CL_ann_sa[idx,..fgroup])

#WT sensitive
#idx <- CL_ann_sa$CellLineName %in% c('OVCA 420','HT-115')
#CGroups[['WT_sensitive']] <-  unique(CL_ann_sa[idx,..fgroup])

#load isobologram data and define quantities to use
dt_smooth <- dt_QCS_combo$SmoothMatrix
clids_drugs <- unique(select(dt_smooth, c('clid','CellLineName','Gnumber', 'Gnumber_2', 'DrugNamePlot', 'DrugNamePlot_2')))
ord_cols <- c('Gnumber', 'Gnumber_2','clid')
clids_drugs <- data.table::setorderv(clids_drugs, ord_cols)

#convert Concentrations to Concentrations_free using FuDrugs
FuDrugs <- data.frame(read.csv(file.path(cwd, data_dir, dda_dir, 'FuDrugs.csv')))
used_FBS_perc = used_FBS_perc_select
idx <- (FuDrugs$FBS_perc==used_FBS_perc)
FuDrugs <- FuDrugs[idx, ]
#check there aren't duplicates
unFuDrugs <- unique(FuDrugs$DrugNamePlot)
if (!(length(unFuDrugs)==length(FuDrugs$DrugNamePlot))){
  stop('Duplicate entries in unFuDrugs, cannot proceed.')
}

##### Define function to add free drug concentrations to a dataset #####
add_free_Concentrations <- function(dt_data, fu_drugs){
  #add fields for free Concentrations
  dt_data$free_Concentration <- NaN
  dt_data$free_Concentration <- as.numeric(dt_data$free_Concentration)
  dt_data$free_Concentration_2 <- NaN
  dt_data$free_Concentration_2 <- as.numeric(dt_data$free_Concentration_2)
  #add free drug concentrations
  for (i in 1:length(fu_drugs$DrugNamePlot)){
    #for Drug
    idx <- dt_data$DrugNamePlot == fu_drugs[i, c('DrugNamePlot')]
    dt_data[ idx ,c('free_Concentration')] <-  dt_data[idx ,c('Concentration')] * fu_drugs[i, c('fu_FBS')]
    #for Drug_2
    idx <- dt_data$DrugNamePlot_2 == fu_drugs[i, c('DrugNamePlot')]
    dt_data[ idx ,c('free_Concentration_2')] <-  dt_data[idx ,c('Concentration_2')] * fu_drugs[i, c('fu_FBS')]
  }
  return(dt_data)
}
##### END FUNCTION #####

#add free concentration vectors  
dt_smooth <- add_free_Concentrations(dt_smooth, FuDrugs)

#merge doses 
DDoses_temp <- DDoses
DDoses_temp$Dose_ID_2 <- DDoses_temp$Dose_ID
DDoses_temp$DrugNamePlot_2 <- DDoses_temp$DrugNamePlot
clids_drugs <- merge(clids_drugs, dplyr::select(DDoses_temp, c('DrugNamePlot','Dose_ID')), by='DrugNamePlot',  allow.cartesian=TRUE)
clids_drugs <- merge(clids_drugs, dplyr::select(DDoses_temp, c('DrugNamePlot_2','Dose_ID_2')), by='DrugNamePlot_2',  allow.cartesian=TRUE)

#create struture to save results 
col_names_dd <- c('clid', 'CellLineName',
                  'Gnumber','Gnumber_2',
                  'DrugNamePlot', 'DrugNamePlot_2', 
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
  dnp <- clids_drugs[['DrugNamePlot']][i]
  dnp_2 <- clids_drugs[['DrugNamePlot_2']][i]   
  dose_id <- clids_drugs[['Dose_ID']][i]
  dose_id_2 <- clids_drugs[['Dose_ID_2']][i]
  
  #use smooth matrix to get measured values
  idx <- (dt_smooth$clid==clid) &
    (dt_smooth$Gnumber==drug) &
    (dt_smooth$Gnumber_2==drug_2)
  dt_smooth_sub <- dt_smooth[idx,]
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
    DDn <- DDoses[i, c('DrugNamePlot')]
    #get single agent response of each cell line
    idx1 <- (dt_smooth$DrugNamePlot == DDn) & 
      (dt_smooth$free_Concentration_2==0) &
      (dt_smooth$CellLineName %in% celllinename_sel)
    idx2 <- (dt_smooth$DrugNamePlot_2 == DDn) & 
      (dt_smooth$free_Concentration==0) & 
      (dt_smooth$CellLineName %in% celllinename_sel)
    
    #extract SA
    dt_SA <- dt_smooth[idx1,]
    dt_SA_2 <- dt_smooth[idx2,]
    #swap Drug and Drug_2 if 
    if (nrow(dt_SA)<1)  {
      dt_SA <-dt_SA_2
      dt_SA$free_Concentration <- dt_SA$free_Concentration_2
      dt_SA$Gnumber <- dt_SA$Gnumber_2
      dt_SA$DrugNamePlot <- dt_SA$DrugNamePlot_2
    }
    #get the dose and min and max used for calculations 
    dd_ave <-DDoses[i, c('dd')]
    dd_min <- DDoses[i, c('dd_min')]
    dd_max <- DDoses[i, c('dd_max')]
    #take the mean across cell lines that had different smoothing
    dt_SA <- dplyr::select(dt_SA, c('clid', 'CellLineName',
                                    'Gnumber',
                                    'free_Concentration', 'DrugNamePlot', gtf$long))
    dt_SA<- dt_SA %>%
      group_by(clid,CellLineName,Gnumber, free_Concentration, DrugNamePlot) %>%
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
      xlab(paste(unique(dt_SA$DrugNamePlot),'uM'))+
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
    lbl_x <- unique(dt_dd_sub$DrugNamePlot)
    color_x <- 'firebrick2'
    lbl_y <- unique(dt_dd_sub$DrugNamePlot_2)
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
      #labs(title=sprintf('%s\n%s (%guM) x %s (%guM)',cgroups_ids[j],unique(dt_dd_sub$DrugNamePlot),DD,
      #                   unique( dt_dd_sub$DrugNamePlot_2),DD_2)) +
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
    dt_dd_sub_melted <- reshape2::melt(dt_dd_sub, id.vars = c("clid", "CellLineName", "DrugNamePlot","DrugNamePlot_2"),
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
      #labs(title=sprintf('%s\n%s (%guM) x %s (%guM)',cgroups_ids[j],unique(dt_dd_sub_melted$DrugNamePlot),DD,
      #                   unique(dt_dd_sub_melted$DrugNamePlot_2),DD_2)) +
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
      #labs(title=sprintf('%s\n%s (%guM) x %s (%guM)',cgroups_ids[j],unique(dt_dd_sub_melted$DrugNamePlot),DD,
      #                   unique(dt_dd_sub_melted$DrugNamePlot_2),DD_2)) +
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
    dt_smooth <- dt_QCS_combo[[assay]]
    idx <- ((dt_smooth$CellLineName==dt_dd[i, c('CellLineName')]) &
            (dt_smooth$Gnumber==dt_dd[i, c('Gnumber')]) &
            (dt_smooth$Gnumber_2==dt_dd[i, c('Gnumber_2')]))
    dt_smooth <- dt_smooth[idx, ]
    if ('normalization_type' %in% colnames(dt_smooth)){
      idx <- (dt_smooth$normalization_type==gtf$short)
      dt_smooth <- dt_smooth[idx, ]
    }
    if (nrow(dt_smooth)>1) {
      #use free concentrations to plot
      dt_smooth$Concentration <- dt_smooth$free_Concentration
      dt_smooth$Concentration_2 <- dt_smooth$free_Concentration_2
      dt_smooth <- createLogConcentrations(dt_smooth)
      if (mixmax_fields[j]==1) {
        mine <- min(c(0.0, min(na.omit(dt_smooth[,..field]))))
        maxe <- max(c(1.1, max(na.omit(dt_smooth[,..field])))) 
        cdotdose <- 'black'
      } else if (mixmax_fields[j]==2) {
        tope <- max(abs(c(min(na.omit(dt_smooth[,..field])), 
                 max(na.omit(dt_smooth[,..field])))))
        mine <- min(c(-0,3,-tope))
        maxe <- max(c(0.3, tope))
        cdotdose <- 'black'
      }
      limits <- c(mine,maxe)
      #calculate width and height of each tile for geom_tile
      dt_smooth$x <- dt_smooth$logConcentration
      dt_smooth$y <- dt_smooth$logConcentration_2
      ux <- unique(dt_smooth$x)
      uy <- unique(dt_smooth$y)
      diffux <- diff(ux)/2
      diffuy <- diff(uy)/2
      #calculate width left and right (wl, wr) and height up and down (hu, hd)
      dt_smooth$wl=plyr::mapvalues(dt_smooth$x, ux, c(diffux[1],diffux))
      dt_smooth$wr=plyr::mapvalues(dt_smooth$x, ux, c(diffux,diffux[length(diffux)]))
      dt_smooth$hu=plyr::mapvalues(dt_smooth$y, uy, c(diffuy,diffuy[length(diffuy)]))
      dt_smooth$hd=plyr::mapvalues(dt_smooth$y, uy, c(diffuy[1],diffuy))
      #calculate xmin, xmax, ymin, ymax for geom_rect
      dt_smooth$xmin <- (dt_smooth$x - dt_smooth$wl)
      dt_smooth$xmax <- (dt_smooth$x + dt_smooth$wr)
      dt_smooth$ymin <- (dt_smooth$y - dt_smooth$hd)
      dt_smooth$ymax <- (dt_smooth$y + dt_smooth$hu)
      colors_smooth <- colors_fields[[j]]
      p_smooth[[length(p_smooth)+1]] <- ggplot()  +
        geom_rect(data=dt_smooth, 
                  aes_string(xmin = 'xmin', xmax = 'xmax', 
                             ymin = 'ymin', ymax = 'ymax', 
                             fill = field)) +
        labs(title=sprintf('%s\n%s x %s', 
                           dt_dd[i, c('CellLineName')],
                           dt_dd[i, c('Dose_ID')],
                           dt_dd[i, c('Dose_ID_2')])) +
        xlab(paste(unique(dt_smooth$DrugNamePlot),'uM'))+
        ylab(paste(unique(dt_smooth$DrugNamePlot_2),'uM')) +
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
    dt_smooth <- dt_QCS_combo[[assay]]
    idx <- ((dt_smooth$CellLineName==uclid_drugs[i, c('CellLineName')]) &
              (dt_smooth$Gnumber==uclid_drugs[i, c('Gnumber')]) &
              (dt_smooth$Gnumber_2==uclid_drugs[i, c('Gnumber_2')]))
    dt_smooth <- dt_smooth[idx, ]
    if ('normalization_type' %in% colnames(dt_smooth)){
      idx <- (dt_smooth$normalization_type==gtf$short)
      dt_smooth <- dt_smooth[idx, ]
    }
    
    if (nrow(dt_smooth)>1) {
      #use free concentrations to plot
      dt_smooth$Concentration <- dt_smooth$free_Concentration
      dt_smooth$Concentration_2 <- dt_smooth$free_Concentration_2
      dt_smooth <- createLogConcentrations(dt_smooth)
      if (mixmax_fields[j]==1) {
        mine <- min(c(0.0, min(na.omit(dt_smooth[,..field]))))
        maxe <- max(c(1.1, max(na.omit(dt_smooth[,..field])))) 
        cdotdose <- 'black'
      } else if (mixmax_fields[j]==2) {
        tope <- max(abs(c(min(na.omit(dt_smooth[,..field])), 
                          max(na.omit(dt_smooth[,..field])))))
        mine <- min(c(-0,3,-tope))
        maxe <- max(c(0.3, tope))
        cdotdose <- 'black'
      }
      limits <- c(mine,maxe)
      #calculate width and height of each tile for geom_tile
      dt_smooth$x <- dt_smooth$logConcentration
      dt_smooth$y <- dt_smooth$logConcentration_2
      ux <- unique(dt_smooth$x)
      uy <- unique(dt_smooth$y)
      diffux <- diff(ux)/2
      diffuy <- diff(uy)/2
      #calculate width left and right (wl, wr) and height up and down (hu, hd)
      dt_smooth$wl=plyr::mapvalues(dt_smooth$x, ux, c(diffux[1],diffux))
      dt_smooth$wr=plyr::mapvalues(dt_smooth$x, ux, c(diffux,diffux[length(diffux)]))
      dt_smooth$hu=plyr::mapvalues(dt_smooth$y, uy, c(diffuy,diffuy[length(diffuy)]))
      dt_smooth$hd=plyr::mapvalues(dt_smooth$y, uy, c(diffuy[1],diffuy))
      #calculate xmin, xmax, ymin, ymax for geom_rect
      dt_smooth$xmin <- (dt_smooth$x - dt_smooth$wl)
      dt_smooth$xmax <- (dt_smooth$x + dt_smooth$wr)
      dt_smooth$ymin <- (dt_smooth$y - dt_smooth$hd)
      dt_smooth$ymax <- (dt_smooth$y + dt_smooth$hu)
      colors_smooth <- colors_fields[[j]]
      p_smooth[[length(p_smooth)+1]] <- ggplot()  +
        geom_rect(data=dt_smooth, 
                  aes_string(xmin = 'xmin', xmax = 'xmax', 
                             ymin = 'ymin', ymax = 'ymax', 
                             fill = field)) +
        labs(title=sprintf('%s', 
                           uclid_drugs[i, c('CellLineName')])) +
        xlab(paste(unique(dt_smooth$DrugNamePlot),'uM'))+
        ylab(paste(unique(dt_smooth$DrugNamePlot_2),'uM')) +
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





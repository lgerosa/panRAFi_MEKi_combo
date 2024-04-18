library(SummarizedExperiment)
library(MultiAssayExperiment)
library(ggplot2)
library(gridExtra)
library(dplyr)

#set working directory
cwd <- "/gstore/home/gerosal/projects/work/Belva_Cobi_combo_2023_07"
setwd(cwd)



#Load the xenograft data

studyIDs <- list(
  #c('1232', 'Disco'),  #GDC-6036 TEAD
  c('55182', 'Divos'), # IPC298
  c('49622', 'Divos')  # A375 
) 


study <- list()
for (i in 1:length(studyIDs)){
  #retrieve MAE
  ID <- studyIDs[[i]][1]
  Platform <- studyIDs[[i]][2]
  studyTemp <- readRDS(file.path(cwd, data_dir, invivo_dir, paste('MAE_',Platform,'_',ID, '.RDS', sep='')))
  study[[ID]] <- studyTemp
}


#plot each study
assayPlotName <- 'TUMOR_VOLUME'
p <- list()
for (i in 1:length(study)){
  
  ID <- names(study)[i]
  name_SE <- names(study[[i]])
  
  #extract data
  studyMetaData <- metadata(study[[i]])
  treatments <- studyMetaData$treatments
  groups <- studyMetaData$groups
  studyMasterData <- study[[i]][[name_SE]]
  assaysStudy <- assays(studyMasterData)
  metaDataStudy <- metadata(studyMasterData)
  rowDataStudy <- as.data.frame(rowData(studyMasterData))
  colDataStudy <-  as.data.frame(colData(studyMasterData))
  assayPlot <- assay(studyMasterData, assayPlotName)
  
  #convert some fields to right type
  rowDataStudy$numTimeValue <- as.numeric(rowDataStudy$timeValue)

  
  #melt the assay data
  assayPlotMelt <- as.data.frame(reshape2::melt(assayPlot))
  colnames(assayPlotMelt) <- c('rowid', 'colid', 'numVal')
  assayPlotMelt$colid <- as.character(assayPlotMelt$colid)
  
  #change name of columns because there is a bug in merge later on otherwise..
  rowDataStudy$rowid <- rowDataStudy$id
  colDataStudy$colid <- colDataStudy$id
  
  #annotate with row and column information
  assayPlotMelt <- merge(assayPlotMelt, colDataStudy, by='colid', all.x=TRUE)
  assayPlotMelt <- merge(assayPlotMelt, rowDataStudy, by='rowid', all.x=TRUE)
  #annotate with group information
  groups$groups.id <- groups$id
  groups <- dplyr::select(groups, -c('id'))
  assayPlotMelt <- merge(assayPlotMelt, groups, by='groups.id', all.x=TRUE)
   
  
  ### PLOT INDIVIDUAL MOUSE TUMOR VOLUMES AND AGGREGATES
  
  #plot each mouse as individual
  p[[i]] <- list()
  p[[i]][[1]] <- ggplot(data=assayPlotMelt, aes(x=numTimeValue, y=numVal, color=animalCode)) +
    geom_line() +
    scale_y_continuous(trans = 'log10') +
    ggtitle(sprintf('%s: Individual mice', ID)) +
    xlab('Time (days)') +
    ylab(assayPlotName) +
    theme_minimal() #+
  #theme(legend.position="top")
  
  
  #plot each mouse as individual but with group coloring
  p[[i]][[2]] <- ggplot(data=assayPlotMelt, aes(x=numTimeValue, y=numVal, text=animalCode , color=groupName)) +
    geom_line() +
    scale_y_continuous(trans = 'log10') +
    ggtitle(sprintf('%s: Individual mice', ID)) +
    xlab('Time (days)') +
    ylab(assayPlotName) +
    theme_minimal() #+ 
  #theme(legend.position="top")
  
  assayPlotMelt$groupIdFactor <- factor(assayPlotMelt$groups.id)
  aggr_meas <- assayPlotMelt %>%
    group_by(groupName, groups.id, numTimeValue) %>%
    summarise(mean_numVal = mean(numVal, na.rm=TRUE),
              sd_numVal = sd(numVal, na.rm=TRUE))
  
  #plot means
  p[[i]][[3]] <- ggplot(data=aggr_meas, aes(x=numTimeValue, y=mean_numVal, color=groupName)) +
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin=mean_numVal-sd_numVal, ymax=mean_numVal+sd_numVal), 
                  width=.2,  
                  position=position_dodge(0.05)) +
    #geom_ribbon(aes(x=numTimeValue, y=mean_numVal, ymax=mean_numVal+sd_numVal, ymin=mean_numVal-sd_numVal, color=groupName), 
    #            alpha=0.5) +
    scale_y_continuous(trans = 'log10') +
    ggtitle(sprintf('%s %s: Mice group average', assayPlotMelt$modelName, ID)) +
    xlab('Time (days)') +
    ylab(assayPlotName) +
    theme_minimal() 
  
}

#save figures in pdf
pdffilepath <- file.path(cwd, data_dir, invivo_dir, 'Divos_Tumor_Volumes.pdf')
pdf(pdffilepath, width= 8  , height= 20)
for (i in 1:length(p)){
  print(grid.arrange(grobs = p[[i]], ncol=1))
}
dev.off()

library(codetools)
library(dplyr)
library(modules)

### FUNCTION plotHeatMapCombo ###

#this function plots heatmaps for combos
plotHeatMapCombo <- function(QCS_values, value, limits, colors) {
  #factorize concentrations to have equally spaced heatmaps
  QCS_values$Concentration <- factor(QCS_values$Concentration)
  QCS_values$Concentration_2 <- factor(QCS_values$Concentration_2)
  p<- ggplot(QCS_values, aes_string(x = 'Concentration', y = 'Concentration_2', fill = value))+
    geom_tile() +
    labs(title=unique(QCS_values$CellLineName)) +
    xlab(paste(unique(QCS_values$DrugNamePlot),'(uM)'))+
    ylab(paste(unique(QCS_values$DrugNamePlot_2),'(uM)')) +
    scale_fill_gradientn(colours = colors, limits=limits) +
    ggpubr::theme_pubr() +
    theme(text = element_text(size=7),
          axis.text = element_text(size = 7),
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right"
    ) 
  return(p)
}


### FUNCTION ploDoseResponseCombo ###

ploDoseResponseCombo <- function(QCS_values, drug_x, drug_f, conc_x, conc_f, value, limits) {
  #transform second concentration in factors
  QCS_values[[conc_f]] <- factor(QCS_values[[conc_f]])
  #define colors
  ncolors <- length(unique(QCS_values[[conc_f]])) 
  colors <- viridis(ncolors)
  p <-ggplot(QCS_values, aes_string(x = conc_x, y=value)) +
    geom_hline(yintercept=0, color='gray80') +
    geom_point( aes_string(color = conc_f)) + 
    geom_line( aes_string(color = conc_f)) + 
    scale_y_continuous(limits = limits) +
    scale_x_continuous(trans = 'log10') +
    scale_colour_manual(values=colors) +
    labs(title=unique(QCS_values[['CellLineName']])) +
    xlab(paste(unique(QCS_values[[drug_x]]),'(uM)'))+
    ylab(value) +
    labs(color=paste(unique(QCS_values[[drug_f]]),'(uM)')) +
    ggpubr::theme_pubr() +
    theme(text = element_text(size=7),
          axis.text = element_text(size = 7),
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right"
    )
  return(p)
}

### FUNCTION createLogConcentrations ### 

#this function create a log Concentration field for combos to facilitate plotting
createLogConcentrations <- function(QCS_values){
  #convert concentrations to log
  QCS_values$logConcentration <- log10(QCS_values$Concentration)
  QCS_values$logConcentration_2 <- log10(QCS_values$Concentration_2)
  #replace the -Inf value coming from the 0 dose with one step less in the dose dilution 
  idx_inf<- (QCS_values$Concentration == 0)
  doses<- sort(unique(QCS_values$logConcentration))
  zero_value <- doses[2]  + (doses[2] - doses[3])
  QCS_values[idx_inf]$logConcentration <- zero_value
  #do the same for Conc_2 
  idx_inf<- (QCS_values$Concentration_2 == 0)
  doses<- sort(unique(QCS_values$logConcentration_2))
  zero_value <- doses[2]  + (doses[2] - doses[3])
  QCS_values[idx_inf]$logConcentration_2 <- zero_value
  return(QCS_values)
}



plotComboAgentFit <- function(combo, clines_drugs, gtf, plot_ID) {
  
  #define allowd plots
  plot_ID_allowed <-c("RawTreated", 
                      "Normalized", 
                      "Averaged", 
                      "SmoothMatrix", 
                      "isobolograms",
                      "isobologramsIC50",
                      "HSAExcess", 
                      "BlissExcess")
  plot_ID_allowed <- c(plot_ID_allowed, 
                       paste(plot_ID_allowed, 'swapped', sep='_'))
  
  #test plot_ID are valid
  if (!(all(plot_ID %in% plot_ID_allowed))) {
    stop('plotComboAgentFit, plot ID not valid.')     
  }
  
  #extract assays to extract and eventually plot
  assay_ID <- c("RawTreated", 
                "Normalized", 
                "Averaged", 
                "isobolograms",
                "excess", 
                "scores", 
                'HSAExcess', 
                'BlissExcess', 
                'SmoothMatrix', 
                'HSAScore', 
                'BlissScore') 
  
  #extract cell lines and drugs to plot
  clines <- clines_drugs$clines
  drug <- clines_drugs$Gnumber
  drug_2 <- clines_drugs$Gnumber_2
  #vector of ggplots
  p <- list()
  #for each cell line and drug combo combination
  for (i in 1:length(clines)){
    #print(i)
    #extract the relevant clines and drug combination
    combo_sub <- list()
    for (j in 1:length(assay_ID)){
      aID <- assay_ID[[j]]
      idx <- (combo[[aID]]$clines==clines[i]) &
        (combo[[aID]]$Gnumber==drug[i]) &
        (combo[[aID]]$Gnumber_2==drug_2[i])
      combo_sub[[aID]] <- combo[[aID]][idx]
    }  
    
    #plot dose response raw treated (swapping the two drugs if needed)
    field <- 'ReadoutValue'
    dt_rawtreated <- combo_sub[['RawTreated']]
    mine <- min(c(0.0, min(na.omit(dt_rawtreated[,..field]))))
    maxe <- max(c(1.1, max(na.omit(dt_rawtreated[,..field]))))
    limits <- c(mine,maxe)
    if ( 'RawTreated' %in% plot_ID) {
      #Drug on x axis and Drug_2 as factors
      p[[length(p)+1]] <- ploDoseResponseCombo(dt_rawtreated,
                                               'DrugNamePlot','DrugNamePlot_2',
                                               'Concentration','Concentration_2',
                                               field, limits)
    }
    if ( 'RawTreated_swapped' %in% plot_ID) {
      #Drug on x axis and Drug_2 as factors
      p[[length(p)+1]] <-  ploDoseResponseCombo(dt_rawtreated,
                                                'DrugNamePlot_2','DrugNamePlot',
                                                'Concentration_2','Concentration',
                                                field, limits) 
    }
    
    #plot dose response averaged (swapping the two drugs if needed)
    field <- sprintf('%s', gtf$long)
    dt_averaged <- combo_sub[['Averaged']]
    groups <- c("normalization_type")
    wide_cols <- c('x', 'x_std')
    dt_averaged <- gDRutils::flatten(dt_averaged, groups = groups, wide_cols = wide_cols)
    
    
    mine <- min(c(0.0, min(na.omit(dt_averaged[,..field]))))
    maxe <- max(c(1.1, max(na.omit(dt_averaged[,..field]))))
    limits <- c(mine,maxe)
    if ('Averaged' %in% plot_ID) {
      #Drug on x axis and Drug_2 as factors
      p[[length(p)+1]] <- ploDoseResponseCombo(dt_averaged,
                                               'DrugNamePlot','DrugNamePlot_2',
                                               'Concentration','Concentration_2',
                                               field, limits)
    }
    if ( 'Averaged_swapped' %in% plot_ID) {
      #Drug on x axis and Drug_2 as factors
      p[[length(p)+1]] <-  ploDoseResponseCombo(dt_averaged,
                                                'DrugNamePlot_2','DrugNamePlot',
                                                'Concentration_2','Concentration',
                                                field, limits)
    }
    
    #plot dose response smoothed (swapping the two drugs)
    field <- sprintf('%s', gtf$long)
    dt_smooth <- combo_sub[['SmoothMatrix']]
    groups <- c("normalization_type")
    wide_cols <- c('x')
    dt_smooth <- gDRutils::flatten(dt_smooth, groups = groups, wide_cols = wide_cols)
    
    
    mine <- min(c(0.0, min(na.omit(dt_smooth[,..field]))))
    maxe <- max(c(1.1, max(na.omit(dt_smooth[,..field])))) 
    limits <- c(mine,maxe)
    if ( 'SmoothMatrix' %in% plot_ID) {
      #Drug on x axis and Drug_2 as factors
      p[[length(p)+1]] <- ploDoseResponseCombo(dt_smooth, 
                                               'DrugNamePlot','DrugNamePlot_2',
                                               'Concentration','Concentration_2',
                                               field, limits)
    }
    if ( 'SmoothMatrix_swapped' %in% plot_ID) {
      #Drug on x axis and Drug_2 as factors
      p[[length(p)+1]] <- ploDoseResponseCombo(dt_smooth, 
                                               'DrugNamePlot_2','DrugNamePlot',
                                               'Concentration_2','Concentration',
                                               field, limits) 
    }
    
    #colors for RV and GR
    field <- gtf$long
    #use a coloring schema and min value different between GR and RV
    if (field=='GRvalue'){
      colors <-  c(rev(rocket(60)[10:60]), viridis(50))
      minF <- -1.0
      maxF <- 1.1
    }else{
      colors <- viridis(50)
      minF <- 0.0
      maxF <- 1.1
    }
    
    field <- sprintf('%s', gtf$long)
    dt_smooth <- combo_sub[['SmoothMatrix']]
    groups <- c("normalization_type")
    wide_cols <- c('x')
    dt_smooth <- gDRutils::flatten(dt_smooth, groups = groups, wide_cols = wide_cols)
    
    mine <- min(c(minF, min(na.omit(dt_smooth[,..field]))))
    maxe <- max(c(maxF, max(na.omit(dt_smooth[,..field])))) 
    limits <- c(mine,maxe)
    
    if ( 'isobolograms' %in% plot_ID) {
      p[[length(p)+1]]<- plotHeatMapCombo(dt_smooth, field, limits, colors)
    }
    
    if ( 'isobolograms_swapped' %in% plot_ID) {
      p[[length(p)+1]]<- plotHeatMapCombo(dt_smooth, field, limits, colors) + 
        coord_flip()
    }
    
    
    #plot HSAExcess
    colors <- colorRampPalette(c("royalblue3", "royalblue1", "grey95" , "grey95" , "firebrick1", "firebrick3"))(51)
    dt_HSAExcess <- combo_sub[['HSAExcess']]
    dt_HSAExcess <- dt_HSAExcess[dt_HSAExcess$normalization_type==gtf$short]
    mine <- min(c(-0.5, min(na.omit(dt_HSAExcess$x))))
    maxe <- max(c(0.5, max(na.omit(dt_HSAExcess$x)))) 
    limits <- c(mine,maxe)
    idx <- ((dt_HSAExcess$Concentration==0) | (dt_HSAExcess$Concentration_2==0))
    dt_HSAExcess[idx]$x <- 0
    #plot HSA Excess
    if ( 'HSAExcess' %in% plot_ID) {
      p[[length(p)+1]]<- plotHeatMapCombo(dt_HSAExcess, 'x', limits, colors) +
        labs(fill="HSA excess")
    }
    #plot HSA Excess swapped
    if ( 'HSAExcess_swapped' %in% plot_ID) {
      p[[length(p)+1]]<- plotHeatMapCombo(dt_HSAExcess, 'x', limits, colors) +
        labs(fill="HSA excess") +
        coord_flip()
    }
    
    #plot BlissExcess
    dt_BlissExcess <- combo_sub[['BlissExcess']]
    dt_BlissExcess <- dt_BlissExcess[dt_BlissExcess$normalization_type==gtf$short]
    mine <- min(c(-0.5, min(na.omit(dt_BlissExcess$x))))
    maxe <- max(c(0.5, max(na.omit(dt_BlissExcess$x)))) 
    limits <- c(mine,maxe)
    idx <- ((dt_BlissExcess$Concentration==0) | (dt_BlissExcess$Concentration_2==0))
    dt_BlissExcess[idx]$x <- 0
    #plot Bliss Excess
    if ( 'BlissExcess' %in% plot_ID) {
      p[[length(p)+1]]<- plotHeatMapCombo(dt_BlissExcess, 'x', limits, colors) +
        labs(fill="Bliss excess") 
    }
    #plot Bliss Excess swapped
    if ( 'BlissExcess_swapped' %in% plot_ID) {
      p[[length(p)+1]]<- plotHeatMapCombo(dt_BlissExcess, 'x', limits, colors) +
        labs(fill="Bliss excess") +
        coord_flip()
    }
    
    
    
    #Plot Isobologram: first plot smooth data with complicate spacing
    #get min and max from isoblogram and also log scale values in smooth
    
    
    #not we haven't implemented the swapped one yet
    if ('isobologramsIC50' %in% plot_ID) {
      
      colors <- viridis(50)
      field <- sprintf('%s', gtf$long)
      dt_smooth <- combo_sub[['SmoothMatrix']]
      groups <- c("normalization_type")
      wide_cols <- c('x')
      dt_smooth <- gDRutils::flatten(dt_smooth, groups = groups, wide_cols = wide_cols)
      if (nrow(dt_smooth)>1) {
        dt_smooth <- createLogConcentrations(dt_smooth)
        mine <- min(c(0.0, min(na.omit(dt_smooth[,..field]))))
        maxe <- max(c(1.1, max(na.omit(dt_smooth[,..field])))) 
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
        colors_smooth <- viridis(50)
        p_smooth <- ggplot()  +
          geom_rect(data=dt_smooth, 
                    aes_string(xmin = 'xmin', xmax = 'xmax', 
                               ymin = 'ymin', ymax = 'ymax', 
                               fill = field)) +
          labs(title=unique(dt_smooth$CellLineName)) +
          xlab(paste(unique(dt_smooth$DrugNamePlot),'uM'))+
          ylab(paste(unique(dt_smooth$DrugNamePlot_2),'uM')) +
          scale_fill_gradientn(colours = colors_smooth, limits=limits) +
          ggpubr::theme_pubr() +
          theme(text = element_text(size=7),
                axis.text = element_text(size = 7),
                axis.text.x = element_text(angle = 0, hjust = 0.5),
                plot.title = element_text(hjust = 0.5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.position = "right"
          )
        
        #plot isobologram
        ic_sel=c("0.5")
        dt_isobolograms <- combo_sub[['isobolograms']]
        dt_isobolograms <- dt_isobolograms[dt_isobolograms$normalization_type==gtf$short]
        dt_isobolograms <- dt_isobolograms[dt_isobolograms$iso_level %in% ic_sel]
        if (nrow(dt_isobolograms)>1) {
          #prepare isobologram for proper plotting with x_pos and x_pos_ref
          dt_isob_proper <- dplyr::select(dt_isobolograms, -c('pos_x_ref','pos_y_ref'))
          dt_isob_proper_2 <- dplyr::select(dt_isobolograms, -c('pos_x','pos_y'))
          dt_isob_proper_2 <- dplyr::rename(dt_isob_proper_2, pos_x = pos_x_ref , pos_y = pos_y_ref)
          dt_isob_proper$iso_source <- 'measured'
          dt_isob_proper_2$iso_source <- 'expected'
          dt_isob_proper <- rbind(dt_isob_proper, dt_isob_proper_2)
          colors_iso <- colorRampPalette(c("red", "purple"))(length(unique(dt_isob_proper$iso_level)))
          dt_isob_proper$iso_level <- factor(dt_isob_proper$iso_level)
          dt_isob_proper$iso_source <- factor(dt_isob_proper$iso_source)
          #add isobologram to plot
          p_smooth <- p_smooth + 
            geom_line(data=dt_isob_proper, 
                      aes(x = pos_y, y = pos_x, color = iso_level, linetype= iso_source) ,size = 2) +
            scale_colour_manual(values = colors_iso) +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0)) 
        }
        #add plot
        p[[length(p)+1]]<- p_smooth 
      }else{
        p[[length(p)+1]]<- ggplot()
      }
    }
    
    if ('isobologramsIC50_swapped' %in% plot_ID) {  
      stop('isobologramsIC50_swapped not implemented yet')
    }
    
  }
  
  return(p)
}

getProjectedPKEffects_dep <- function(smooth_loc,CGroups_loc,DDoses_loc,gtf_loc){
  # Projects free drug concentrations derived from PK range (min, mean, max) to in-vitro response
  # Previous version which is more suited to min/max/avg PK, but less suited for timecourse
  
  #load isobologram data and define quantities to use
  clines_drugs <- unique(dplyr::select(smooth_loc, c('CellLineName','DrugName', 'DrugName_2')))
  ord_cols <- c('DrugName', 'DrugName','CellLineName')
  clines_drugs <- data.table::setorderv(clines_drugs, ord_cols)
    
  #merge doses 
  DDoses_temp <- DDoses_loc
  DDoses_temp$Dose_ID_2 <- DDoses_temp$Dose_ID
  DDoses_temp$DrugName_2 <- DDoses_temp$DrugName
  clines_drugs <- merge(clines_drugs, dplyr::select(DDoses_temp, c('DrugName','Dose_ID')), by='DrugName',  allow.cartesian=TRUE)
  clines_drugs <- merge(clines_drugs, dplyr::select(DDoses_temp, c('DrugName_2','Dose_ID_2')), by='DrugName_2',  allow.cartesian=TRUE)  
    
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
  allNames <- names(CGroups_loc)
  for (i in 1:length(allNames)){
    allCellLines <- c(allCellLines, CGroups_loc[[allNames[i]]][[1]])
  }
  allCellLines <- unique(allCellLines)
  idx <- clines_drugs$CellLineName %in% allCellLines
  clines_drugs <- clines_drugs[idx,]
  
  #extract SA or combo metrics for each drug and cell line at defined concentration range
  for (i in 1:nrow(clines_drugs)) {
    #extract cell lines and drugs
    cln <- clines_drugs[['CellLineName']][i]
    drug <- clines_drugs[['DrugName']][i]
    drug_2 <- clines_drugs[['DrugName_2']][i]
    dose_id <- clines_drugs[['Dose_ID']][i]
    dose_id_2 <- clines_drugs[['Dose_ID_2']][i]
    
    #use smooth_loc matrix to get measured values
    idx <- (smooth_loc$CellLineName==cln) &
      (smooth_loc$DrugName==drug) &
      (smooth_loc$DrugName_2==drug_2) 
    dt_smooth_loc_sub <- smooth_loc[idx,]
    #use interpolation to extract growth metrics at SA and combo doses
    smooth_loc_matrix <- reshape2::acast(dt_smooth_loc_sub, free_Concentration ~ free_Concentration_2, value.var = gtf_loc$long)
    idx <- (dt_smooth_loc_sub$free_Concentration==0)
    x <- sort(dt_smooth_loc_sub$free_Concentration_2[idx])
    idx <- (dt_smooth_loc_sub$free_Concentration_2==0)
    y <- sort(dt_smooth_loc_sub$free_Concentration[idx])
    
    #extract values to interpolate
    idx <- (DDoses_loc$Dose_ID==dose_id)
    dd <- DDoses_loc[idx, c('dd', 'dd_min', 'dd_max')]
    dd <- c(dd$dd, dd$dd_min, dd$dd_max)
    idx <- (DDoses_loc$Dose_ID==dose_id_2)
    dd_2 <- DDoses_loc[idx, c('dd', 'dd_min', 'dd_max')]
    dd_2 <- c(dd_2$dd, dd_2$dd_min, dd_2$dd_max)
    #interpolate at SA and combo areas
    iv <- pracma::interp2(x  = x, 
                          y  = y, 
                          Z = smooth_loc_matrix, 
                          xp = c(dd_2, c(0,0,0), dd_2), 
                          yp = c(c(0,0,0), dd, dd)
    )
    
    #HSA direct calculation
    # Andrew - There are cases where the mean HSA and min/max HSA
    # result from different agents
    iv_hsa <- c()
    iv_hsa[1] <- min(c(iv[1],iv[4]),na.rm = TRUE)
    iv_hsa[2] <- min(c(iv[2],iv[5]),na.rm = TRUE)
    iv_hsa[3] <- min(c(iv[3],iv[6]),na.rm = TRUE)
    #if((iv[1]<iv[4]) | is.na(iv[4])){
    #  iv_hsa[1] <- iv[1]
    #  iv_hsa[2] <- iv[2]
    #  iv_hsa[3] <- iv[3]
    #}else{
    #  iv_hsa[1] <- iv[4]
    #  iv_hsa[2] <- iv[5]
    #  iv_hsa[3] <- iv[6]
    #}
    
    #calculate Bliss using gDR
    sa2 <- dt_smooth_loc_sub[dt_smooth_loc_sub[['free_Concentration']] == 0,]
    sa2 <- dplyr::select(sa2, c('free_Concentration', 'free_Concentration_2', gtf_loc$long))
    colnames(sa2) <- c("free_Concentration",   "free_Concentration_2", "x")
    sa1 <- dt_smooth_loc_sub[dt_smooth_loc_sub[['free_Concentration_2']] == 0,]
    sa1 <- dplyr::select(sa1, c('free_Concentration', 'free_Concentration_2', gtf_loc$long))
    colnames(sa1) <- c("free_Concentration",   "free_Concentration_2", "x")
    dt_bliss <- gDRcore::calculate_Bliss(sa1, 'free_Concentration', sa2, 'free_Concentration_2', 'x')
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
    dd_temp <- c(cln, 
                 drug, drug_2, 
                 dose_id, dose_id_2, gtf_loc$long,
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
  dt_dd[, c(7:27)] <- sapply(dt_dd[, c(7:27)], as.numeric)
  return(dt_dd)
}

getProjectedPKEffects <- function(smooth_loc,dose_conditions_loc,gtf_loc){
  # Projects free drug concentrations derived from PK range (min, mean, max) to in-vitro response
  
  
  
  unique_conditions <- unique(select(smooth_loc,DrugName,DrugName_2,CellLineName))
  
  for (i in 1:nrow(unique_conditions)){
    do <- filter(smooth_loc,DrugName == unique_conditions$DrugName[i], DrugName_2 == unique_conditions$DrugName_2[i], CellLineName == unique_conditions$CellLineName[i])
    
    dc <- filter(dose_conditions_loc,DrugName == unique_conditions$DrugName[i], DrugName_2 == unique_conditions$DrugName_2[i])
    
    dc$CellLineName <- unique_conditions$CellLineName[i]
    #rep(list(x), 2)
    
    do_matrix <- reshape2::acast(do, free_Concentration ~ free_Concentration_2, value.var = gtf_loc$long)
    idx <- (do$free_Concentration==0)
    x <- sort(do$free_Concentration_2[idx])
    idy <- (do$free_Concentration_2==0)
    y <- sort(do$free_Concentration[idy])
    
    #extract values to interpolate
    
    #interpolate at SA and combo areas
    SA_2 <- pracma::interp2(x  = x, 
                            y  = y, 
                            Z = do_matrix, 
                            xp = c(dc$free_Concentration_2), 
                            yp = c(numeric(length(dc$free_Concentration)))
    )
    SA <- pracma::interp2(x  = x, 
                          y  = y, 
                          Z = do_matrix, 
                          xp = c(numeric(length(dc$free_Concentration_2))), 
                          yp = c(dc$free_Concentration)
    )
    Combo <- pracma::interp2(x  = x, 
                             y  = y, 
                             Z = do_matrix, 
                             xp = c(dc$free_Concentration_2), 
                             yp = c(dc$free_Concentration)
    )
    
    # Calculate HSA
    HSA <- pmin(SA_2,SA,na.rm = TRUE)
    
    # Calculate Bliss scores
    sa2 <- do[do[['free_Concentration']] == 0,]
    sa2 <- dplyr::select(sa2, c('free_Concentration', 'free_Concentration_2', gtf_loc$long))
    colnames(sa2) <- c("free_Concentration",   "free_Concentration_2", "x")
    sa1 <- do[do[['free_Concentration_2']] == 0,]
    sa1 <- dplyr::select(sa1, c('free_Concentration', 'free_Concentration_2', gtf_loc$long))
    colnames(sa1) <- c("free_Concentration",   "free_Concentration_2", "x")
    dt_bliss <- gDRcore::calculate_Bliss(sa1, 'free_Concentration', sa2, 'free_Concentration_2', 'x')
    bliss_matrix <- reshape2::acast(dt_bliss, free_Concentration ~ free_Concentration_2, value.var = 'metric')
    y_bliss <- sort(sa1$free_Concentration)
    x_bliss <- sort(sa2$free_Concentration_2)
    iv_bliss <- pracma::interp2(x  = x_bliss, 
                                y  = y_bliss, 
                                Z = bliss_matrix, 
                                xp = dc$free_Concentration_2, 
                                yp = dc$free_Concentration
    )
    
    dc$SA <- SA
    dc$SA_2 <- SA_2
    dc$Combo <- Combo
    dc$HSA <- HSA
    dc$Bliss <- iv_bliss
    
    if (i == 1){
      projected_doses <- dc
    }
    else{
      projected_doses <- bind_rows(projected_doses,dc)
    }
    
  }
  return(projected_doses)
  
}

flatten_Doses <- function(smooth_loc,DDoses_loc){
  # Converts dose scheme into format with each row corresponding to specific PK condition
  clines_drugs <- unique(dplyr::select(smooth_loc, c('DrugName', 'DrugName_2')))
  ord_cols <- c('DrugName', 'DrugName')
  clines_drugs <- data.table::setorderv(clines_drugs, ord_cols)
  
  
  DDoses_loc$Dose_ID_2 <- DDoses_loc$Dose_ID
  DDoses_loc$DrugName_2 <- DDoses_loc$DrugName
  clines_drugs <- merge(clines_drugs, dplyr::select(DDoses_loc, c('DrugName','Dose_ID','dd','dd_min','dd_max')), by='DrugName',  allow.cartesian=TRUE)
  clines_drugs <- merge(clines_drugs, dplyr::select(DDoses_loc, c('DrugName_2','Dose_ID_2','dd','dd_min','dd_max')), by='DrugName_2',  allow.cartesian=TRUE) 
  
  
  cline_drugs_1 <- dplyr::select(clines_drugs,-dd.x,-dd_min.x,-dd_max.x)
  cline_drugs_1 <- rename(cline_drugs_1, dd = dd.y,dd_min = dd_min.y,dd_max = dd_max.y)
  cline_drugs_1 <- reshape2::melt(cline_drugs_1,measure.vars = c("dd", 
                                                                 "dd_min", 
                                                                 "dd_max") ,value.name = "free_Concentration_2",variable.name = "PK_subcondition")
  
  cline_drugs_2 <- dplyr::select(clines_drugs,-dd.y,-dd_min.y,-dd_max.y)
  cline_drugs_2 <- rename(cline_drugs_2, dd = dd.x,dd_min = dd_min.x,dd_max = dd_max.x)
  cline_drugs_2 <- reshape2::melt(cline_drugs_2,measure.vars = c("dd", 
                                                                 "dd_min", 
                                                                 "dd_max"),value.name = "free_Concentration",variable.name = "PK_subcondition") 
  dose_conditions <- merge(cline_drugs_1,cline_drugs_2)
  return(dose_conditions)
}

get_cycle_length <- function(Dose_name){
  result = switch(   
    strsplit(Dose_name," ")[[1]][3],   
    "BID"= 12,   
    "QD"= 24,
    "QOD"= 48
  ) 
  return(result)
}


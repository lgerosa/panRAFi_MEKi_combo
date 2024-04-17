### FUNCTION plotHeatMapCombo ###

#this function plots heatmaps for combos
plotHeatMapCombo <- function(QCS_values, value, limits, colors) {
  #factorize concentrations to have equally spaced heatmaps
  QCS_values$Concentration <- factor(QCS_values$Concentration)
  QCS_values$Concentration_2 <- factor(QCS_values$Concentration_2)
  p<- ggplot(QCS_values, aes_string(x = 'Concentration', y = 'Concentration_2', fill = value))+
    geom_tile() +
    labs(title=unique(QCS_values$CellLineName)) +
    xlab(paste(unique(QCS_values$DrugName),'(uM)'))+
    ylab(paste(unique(QCS_values$DrugName_2),'(uM)')) +
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
                                               'DrugName','DrugName_2',
                                               'Concentration','Concentration_2',
                                               field, limits)
    }
    if ( 'RawTreated_swapped' %in% plot_ID) {
      #Drug on x axis and Drug_2 as factors
      p[[length(p)+1]] <-  ploDoseResponseCombo(dt_rawtreated,
                                                'DrugName_2','DrugName',
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
                                               'DrugName','DrugName_2',
                                               'Concentration','Concentration_2',
                                               field, limits)
    }
    if ( 'Averaged_swapped' %in% plot_ID) {
      #Drug on x axis and Drug_2 as factors
      p[[length(p)+1]] <-  ploDoseResponseCombo(dt_averaged,
                                                'DrugName_2','DrugName',
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
                                               'DrugName','DrugName_2',
                                               'Concentration','Concentration_2',
                                               field, limits)
    }
    if ( 'SmoothMatrix_swapped' %in% plot_ID) {
      #Drug on x axis and Drug_2 as factors
      p[[length(p)+1]] <- ploDoseResponseCombo(dt_smooth, 
                                               'DrugName_2','DrugName',
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
          xlab(paste(unique(dt_smooth$DrugName),'uM'))+
          ylab(paste(unique(dt_smooth$DrugName_2),'uM')) +
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
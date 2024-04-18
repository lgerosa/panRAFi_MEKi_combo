# Genentech Inc
# by Luca Gerosa and Andrew Goetz
#
# Prunes PK data
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

# Full time range
time_plot_range <- c(528,528+96) 

# loads belva PK data
belva_pk <- data.frame(read.csv(file.path(cwd, data_dir, 'PK_variability', 'belva_sim_pk_data.csv')))
unique(belva_pk$TRT01P)
#belva_pk <- select(belva_pk,ID,time,IPRED,TRT01P)
#belva_pk <- rename(belva_pk, Dose_ID = TRT01P)
belva_pk <- filter(belva_pk, time == round(time))
belva_pk <- filter(belva_pk,time <= max(time_plot_range), time >= min(time_plot_range))

#belva_pk <- filter(belva_pk, time < )
# loads cobi PK data
cobi_pk <- data.frame(read.csv(file.path(cwd, data_dir, 'PK_variability', 'cobi_sim_pk_data.csv')))
#cobi_pk <- select(cobi_pk,ID,time,IPREDnormal,TRT01P)
#cobi_pk <- rename(cobi_pk, Dose_ID_2 = TRT01P)
cobi_pk$time <- cobi_pk$time*24
cobi_pk <- filter(cobi_pk,time <= max(time_plot_range), time >= min(time_plot_range))


write.csv(belva_pk, file.path(cwd, data_dir, 'PK_variability', 'belva_sim_pk_data_pruned.csv'))
write.csv(cobi_pk, file.path(cwd, data_dir, 'PK_variability', 'cobi_sim_pk_data_pruned.csv'))



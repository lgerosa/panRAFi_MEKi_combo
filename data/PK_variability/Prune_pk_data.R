# Genentech Inc
# by Luca Gerosa and Andrew Goetz
#
# Prunes PK data from large csv file to one more manageable, formats both
# drugs to have similar form.
library(dplyr)

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

# Full time range
time_plot_range <- c(528,528+96) 

# loads belva PK data
belva_pk <- data.frame(read.csv(file.path(cwd, data_dir, 'PK_variability', 'belva_sim_pk_data.csv')))
belva_pk <- filter(belva_pk, time == round(time))
belva_pk <- filter(belva_pk,time <= max(time_plot_range), time >= min(time_plot_range))

# loads cobi PK data
cobi_pk <- data.frame(read.csv(file.path(cwd, data_dir, 'PK_variability', 'cobi_sim_pk_data.csv')))
cobi_pk$time <- cobi_pk$time*24
cobi_pk <- filter(cobi_pk,time <= max(time_plot_range), time >= min(time_plot_range))

# saves csv files
write.csv(belva_pk, file.path(cwd, data_dir, 'PK_variability', 'belva_sim_pk_data_pruned.csv'))
write.csv(cobi_pk, file.path(cwd, data_dir, 'PK_variability', 'cobi_sim_pk_data_pruned.csv'))



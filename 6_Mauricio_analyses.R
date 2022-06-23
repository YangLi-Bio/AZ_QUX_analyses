#########################################################
#                                                       #
#         Integration of IFP data from Mauricio         #
#                                                       #
#########################################################


# Parameters
IFP.dir <- "/fs/ess/scratch/PCON0022/liyang/AstraZeneca/QUX/Mauricio/IPF_Single_Cell_15JUN22/"
out.dir <- "R_files/"


# Load data
library(dplyr)
input.list <- list.dirs(path = IFP.dir) %>% grep("\\/outs$", ., value = T)
input.list

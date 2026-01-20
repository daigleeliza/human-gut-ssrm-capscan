#######################################
# SSRM project

# load libraries
# configure data directories
# source base functions, if needed
#######################################
library()


#--------------------------------------------
# define repository paths
#--------------------------------------------
configDir = paste0(here::here(), "/0-Cluster/config/")
treeInfoDir = paste0(here::here(), "/0-Cluster/workflow/out/treeInfo/")
abundanceDir = paste0(here::here(), "/0-Cluster/workflow/out/dsrAB_CapsuleStool_Abundances/")

dataDir = paste0(here::here(), "/data/")
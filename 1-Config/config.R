#######################################
# SSRM project

# load libraries
# configure data directories
# source base functions, if needed
#######################################
library(readr)
library(janitor)
library(arsenal)
library(renv)
library(here)
library(readxl)
# library(haven)
library(plyr)
library(dplyr)
# library(ggcorrplot)
library(tidyr)
library(reshape2)
library(ggplot2)
library(patchwork)
library(colorspace)
library(ggpubr)
library(grid)
library(gridExtra)
library(stringr)
library(purrr)
library(viridis)
library(tidyverse)
library(factoextra)
#library(pheatmap)
#library(ComplexHeatmap)
library(ape)
library(phytools)
library(treeio)
library(ggtree)
library(tidytree)
library(RColorBrewer)
library(patchwork)

# fix conflicts between packages
library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::count)
conflicts_prefer(dplyr::slice)
conflicts_prefer(wrapr::coalesce)
conflicts_prefer(dplyr::combine)
conflicts_prefer(tibble::has_name)
conflicts_prefer(dplyr::lag)
conflicts_prefer(purrr::modify)
conflicts_prefer(wrapr::pack)
conflicts_prefer(wrapr::unpack)
conflicts_prefer(wrapr::view)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::summarize)
conflicts_prefer(dplyr::summarise)
conflicts_prefer(plyr::mutate)
conflicts_prefer(dplyr::left_join)

#--------------------------------------------
# define repository paths
#--------------------------------------------
configDir = paste0(here::here(), "/0-Cluster/config/")
treeInfoDir = paste0(here::here(), "/0-Cluster/workflow/out/treeInfo/")
raxmlDir = paste0(here::here(), "/0-Cluster/workflow/out/raxmlOutput")
abundanceDir = paste0(here::here(), "/0-Cluster/workflow/out/dsrAB_CapsuleStool_Abundances/")

dataDir = paste0(here::here(), "/data/")
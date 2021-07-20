require(tidyverse)
require(foreach)
require(cowplot)
require(RColorBrewer)
require(seqinr)
library(ape)


source('initial_processing.r')
source('call_mutations.r')

source('determine_first_mutation.r')
source('process_kaplan_meier_curves.r')
source('CalculatingTargetSize_MutRate.R')
source('test_first_mutation.r')
source('plot_figure_2.r')



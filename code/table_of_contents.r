#Import required packages
require('foreach')
require('tidyverse')
require('condSURV')
require('lmtest')
require('cowplot')
require('ggtext')
require('emdbook')
require('doParallel')


##############################################
############# Produce Figure 2 ###############
##############################################

#There are three major components that go into producing Figure 2:

#(A) A survival analysis of the data in Rocheleau 2018, Clin Microbiol Infect
# https://pubmed.ncbi.nlm.nih.gov/28652115/

## (A.1) First, we use the condSURV package to fit conditional survival estimates 
## (i.e., the probably of no resistance in a given year). 
## Note - this data is not *publicly available*. The running of this script must be
## coordinated with the BC Centre for HIV Excellence 
## However, the output of this script is given in $filename, so you can still plot 
## Fig 2 in its entirety. 

source('conditional_survival_analysis_BCdata.r')

##(A. 2) Second, we fit exponential models to the data
source('fit_exponential_decay.r')

#(B) A breakdown of the number of mutations present in 
#     initially drug naive, triple-drug treated, virologically-rebounding individuals    
#(C) and the mutation present among individuals with one DRM

##For both of these analyses (B and C), we initially run a (slow) script to process
##  the raw data and call mutations
##Once these have been run once, they should print intermediate data formats that
##  don't require these to be re-run

source('initial_processing.r')
source('call_mutations.r')

##Then we determine the first mutation (and the total number of mutations)
## - this creates intermediate data structures that will be plot in Figure 2B and C

source('determine_first_mutation.r')

##We also perform statistical tests to test the significance of the NNRTI first 
## (on NNRTI-based therapies) and 3TC first (on PI-based therapies) hypotheses
## To do this, we first need to calculate mutation rates:

source('CalculatingTargetSize_MutRate.R')
source('test_first_mutation.r')

##Figure 2: Once all of these pieces are in place, we plot everything

source('plot_figure_2.r')



##############################################
############# Produce Figure 5 ###############
##############################################

#This second section has the spatially- and temporally- varying models of Figure 5

#This code was originally run in parallel on a cluster


#These have drug specific parameter values that are called in both the spatial and temporal
#models
source('shared_functions.r')

#Temporal heterogeneity model
source('time_runner.r')

#Spatial heterogeneity model
source('space_runner.r')

#Now, we can plot these pieces together into Figure 5
source('plot_figure_5.r')

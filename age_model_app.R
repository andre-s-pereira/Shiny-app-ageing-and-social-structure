age_model <- function(pold, nkin, nsims) {

outputdata <- data.frame(n=rep(50, nsims), n.old.individuals = NA, n.kin = NA, perc.old.individuals= NA, diameter= NA, transitivity= NA, mean.degree= NA) # storage for output
  
for(cursim in 1:nsims){ # for each simulation (cursim = current simulation)
    
## set the model
library(igraph)
library(magrittr)
  
gsize = 50                                                # number of individuals (nodes)
n_kin_groups <- nkin                                       # number of individuals in each kin group. 1 = all individuals are non-kin (i.e., 1 individual per kin group = 50 kin groups). 50= all individuals are kin (i.e., 50 individuals per kin group = 1 kin group)
kingroupraw <- rep(1:n_kin_groups,each=(50/n_kin_groups))
typenames <- c('ook', 'oon', 'oyk', 'oyn', 'yyk', 'yyn')  # names of dyad types. Must fit with code and typeparams. o: Old; y: Young; k: kin; n: nonkin
typeparams <- c(0.33, 0.02, 0.37, 0.05, 0.27, 0.08)       # parameters from data for each dyad type. Must be in the same order as the typenames (the format could be changed for security)
nold <- round(pold*gsize/100)

## draw and store number of group members that are old ####

outputdata[cursim, 'n.old.individuals'] <- nold  # store it in the output
outputdata[cursim, 'n.kin'] <- nkin
outputdata[cursim, 'perc.old.individuals'] <- pold # store the percentage of old individuals in the output
    
## create node data ####
    
agegroup <- c(rep(1,nold), rep(0, gsize-nold))  # vector with agegroup for each individual (old = 1, young = 0)
kingroup <- sample(kingroupraw)             # vector with kingroup for each individual (we must randomise one of the vectors to get independence between age and kingroup, so We randomise this one). 
    
## change the format to dyad-based ####
    
# two vectors of individual names that, when matched, give all the unique dyads once
inames <- rep(1:(gsize-1),(gsize-1):1)           # names of all i individuals (the first individual of each dyad)
jnames <- sequence((gsize-1):1, from = 2:gsize)  # names of all j individuals (the second individual of each dyad)
    
# vectors with the agegroup and kingroup for all i and j individuals, respectively
agegroupalli <- agegroup[inames]
agegroupallj <- agegroup[jnames]
kingroupalli <- kingroup[inames]
kingroupallj <- kingroup[jnames]
    
    
## get dyadtypes ####
    
ndyad <- length(inames)     # number of dyads
dyadtype <- rep(NA, ndyad)  # storage for dyadtype
    
# find and store the type of dyad for each dyad. The type is determined by their ages and kinship.
dyadtype[agegroupalli == 1 & agegroupallj == 1 & kingroupalli == kingroupallj] <- 'ook'  # old, old, kin
dyadtype[agegroupalli == 1 & agegroupallj == 1 & kingroupalli != kingroupallj] <- 'oon'  # old, old, nonkin
dyadtype[agegroupalli == 1 & agegroupallj == 0 & kingroupalli == kingroupallj] <- 'oyk'  # old, young, kin
dyadtype[agegroupalli == 1 & agegroupallj == 0 & kingroupalli != kingroupallj] <- 'oyn'  # old, young, nonkin
dyadtype[agegroupalli == 0 & agegroupallj == 0 & kingroupalli == kingroupallj] <- 'yyk'  # young, young, kin
dyadtype[agegroupalli == 0 & agegroupallj == 0 & kingroupalli != kingroupallj] <- 'yyn'  # young, young, nonkin
    
    
## get links ####
    
link <- rep(NA, ndyad)  # storage for links for all dyads (will contain 1 for connected dyads and 0 for unconnected dyads)
    
for (typenum in seq_along(typenames)){     # for each dyadtype
  curtypename <- typenames[typenum]          # name of current dyadtype
  curtypeparam <- typeparams[typenum]        # parameter for current dyadtype
  ncurtype <- sum(dyadtype == curtypename)   # the number of dyads of the current type
  link[dyadtype == curtypename] <- rbinom(n = ncurtype, size = 1, prob = curtypeparam) #(rtruncnorm(n=1,a=0,b=1,mean=curtypeparam, sd=0.07)))  # draw and store links for the current dyadtype
}
    
    
## create igraph network ####
    
dyaddata <- data.frame(inames, jnames, agegroupalli, agegroupallj, kingroupalli, kingroupallj, dyadtype, link) # dataframe with dyad data, for all dyads
edgelistdata <- dyaddata[link==1,]                                              # dataframe with edgelist in the two first columns, link attributes in the rest (i.e. dyads without links are not included here)
nodedata <- data.frame(nodenames = 1:gsize, agegroup, kingroup)                     # dataframe with node data (not sure whether this is needed but maybe igraph uses it to see nodes that have no links)
inet <-graph_from_data_frame(d=edgelistdata, vertices=nodedata, directed=FALSE) # get network in igraph format
    
## calculate and store network metrics ####
  
curdegrees <- degree(inet) # degrees for all individuals
    
# global network metrics
outputdata[cursim, 'diameter'] <- diameter(inet, unconnected = T) # ANDRE's COMMENT: CHECK THIS ONE, -> check the unconnected=T
outputdata[cursim, 'transitivity'] <- transitivity(inet, type = "global") # ANDRE's COMMENT: CHECK THIS ONE
outputdata[cursim, 'mean.degree'] <- mean(curdegrees)
} ## end of for loop
save(outputdata, file="age_simulation_dataset.Rdata")    

#### Save igraph object for last simulated network

dyaddata <- data.frame(inames, jnames, agegroupalli, agegroupallj, kingroupalli, kingroupallj, dyadtype, link) # dataframe with dyad data, for all dyads
edgelistdata <- dyaddata[link==1,]                                              # dataframe with edgelist in the two first columns, link attributes in the rest (i.e. dyads without links are not included here)
nodedata <- data.frame(nodenames = 1:gsize, agegroup, kingroup)                     # dataframe with node data (not sure whether this is needed but maybe igraph uses it to see nodes that have no links)
nodedata %<>% mutate(agegroup = as.factor(agegroup))
nodedata %<>% mutate(kingroup = as.factor(kingroup))
# dataframe with node data (not sure whether this is needed but maybe igraph uses it to see nodes that have no links)
inet <-graph_from_data_frame(d=edgelistdata, vertices=nodedata, directed=FALSE) %>% 
  as_tbl_graph	
save(inet,file="age_network_plot.Rdata")
save(nodedata,file="nodedata.Rdata")

} ## end of simulation

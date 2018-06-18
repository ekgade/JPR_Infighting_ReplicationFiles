#JPR Replication Files - Fratricide in Rebel Movements: 
#A Network Analysis of Syrian Militant Infighting
# Emily Kalah Gade - ekgade@uw.edu
# 15 May 2018

rm(list = ls())
set.seed(32314)
options(digits = 3)

library(statnet)
library(network)
library(tidyverse)
library(amen)
library(tile)
require(RColorBrewer)
col <- brewer.pal(9, "Set1")
# Main document replication (all suplementary analysis below)



################################
# AME Regression Analysis - Table 5 #
################################

# Loading DV matricies for AME regression analysis for 30 groups 
DV_30<-read.csv("Infighting30_2014_6Sep2017.csv", stringsAsFactors = F, header = T) 

# Making a var with the names of groups in order, to use to match the order of all other variables 
# (This is left as a column because R adds "X" before numeric indicators when it imports, and we have some groups that are abreviated by number)

colnames<-DV_30$X # Taken from row names in this first dataset
rownames(DV_30)<-colnames # Make the list rownames for this dataframe
DV_30<-DV_30[,-1] # Subtract duplicate row
colnames(DV_30)<-colnames # Also rename colnames (as some will be listed as "X1st" rather than "1st")


## Setting up dyadic level and node level controls 

ideol<-read.csv("IdeologyVars_JPR.csv", header = T, stringsAsFactors = F)

# Select coveariate data for only 30 infighting groups we are currently analyzing  
ideol2<-filter(ideol, ideol$GroupCode %in% colnames) 
                  
ideolVars<-ideol2[,-1] # Get rid of long form names 

# Cut down to only vars of intrest for this analysis 
ideolVars<-ideolVars[,-10] # Drop power broker score - alternative measure of power
ideolVars<-ideolVars[,-6:-7] # Drop location and state sponsor, will add in later in matrix form
ideolVars$GroupSize <- as.numeric(gsub(",","",ideolVars$GroupSize)) # Remove commas
ideolVars<-ideolVars[,-7] # Remove group size "full"
idv<-ideolVars[match(colnames, ideolVars$GroupCode),]
idv<-na.omit(idv) # AME won't take NAs (that's why we removed state sponsor etc for now)

## now state sponsor node var (just does a group have a sponsor, Y/N)

Nodespons<-as.data.frame(cbind(ideol2$GroupCode, ideol2$State.Sponsors), stringsAsFactors = F) # Create two column dataframe of state sponsors
Nodespons$V2[Nodespons$V2==""]<-NA # Make sure all blanks are read as NA
Nodespons$V2<- ifelse(is.na(Nodespons$V2), 0, 1) # Replace groups that are NA with Zero (no sponsor) or 1 (sponsor)
names(Nodespons)<-c("GroupCode", "spons")
NodeS<-Nodespons[match(colnames, Nodespons$GroupCode),]
a<-array(NodeS[,2]) # Convert to array 

## Add a node level var for ISIL

NodeS$V3[NodeS$GroupCode=="ISIL"]<-"1" # Give ISIL a 1
NodeS$V3<-as.numeric(ifelse(is.na(NodeS$V3), 0, NodeS$V3)) # Give all non-ISIL groups  zero
ISIL<-array(NodeS[,3]) # Convert to array

## Creating Xn and Xd for AME analysis 
b<-array(idv[,5]) # Average ideology, could use a compoent of ideology instead
c<-array(idv[,6]) # Size
a_a<-c/1000 # Putting them on similar scale for rough interp. 

Xn<-as.array(cbind(b, a_a, a, ISIL)) # AME takes arrays
dimnames(Xn)[[1]]<-colnames # Labling rows (groups), already ordered/matched above
dimnames(Xn)[[2]]<-c("averageId", "size", "spons", "ISIL") # Lable columns


## Creating Xn and Xd for AME analysis  WITHOUT ISIS

Xn_b<-as.array(cbind(b, a_a, a)) # AME takes arrays
dimnames(Xn_b)[[1]]<-colnames # Labling rows (groups), already ordered/matched above
dimnames(Xn_b)[[2]]<-c("averageId", "size", "spons") # Lable columns


## Constructing dyadic vars
ad<-as.data.frame(Xn) # Generate a dataframe from your array
ideol_diff<- outer(ad$averageId, ad$averageId, "-") # Make a matrix where each square represents the difference in ideology scores of a given dyad
ideol_diff<-abs(ideol_diff) # Has to be postive - so take the absolute value 

power_diff<- outer(ad$size, ad$size, "-") # Replicates the above for power
power_diff<-abs(power_diff)

# Load State Sponsor and Location Matrcies 

loc<-read.csv("LocationDiff_Infighting_6June2018.csv", stringsAsFactors = F, header = T) # Load shared location matrix 

x<-loc$X # Deal with row and column names as above
loc<-loc[,-1]
colnames(loc)<-x
rownames(loc)<-x
dim(loc) # Check dimentions 
loc<-as.matrix(loc) # Convert to matrix format
dimnames(loc)[[1]]<-x
dimnames(loc)[[2]]<-x

locFinal<-loc[colnames,] # Make order of "loc" match order of ideology and everything else
locFinal<-locFinal[,colnames] # And match on the other axis... 

## state sponorship

spons<-read.csv("stateSponsonership_InfightingGroups_6June2018.csv", stringsAsFactors = F, header = T)
x<-spons$X 
spons<-spons[,-1]
rownames(spons)<-x
colnames(spons)<-x
spons<-as.matrix(spons)
dimnames(spons)[[1]]<-x
dimnames(spons)[[2]]<-x
sponsFinal<-spons[colnames,]
sponsFinal<-sponsFinal[,colnames]

# Create an empty array to hold your matrcies (dyad vars) for AME
Xd<- array(dim=c(nrow(DV_30), nrow(DV_30), 4)) # Must be same dimentions as number of groups (rows and columns) and have enough "slices" to fit all dyad vars  
Xd[,,1]<- ideol_diff # Average ideology dyad
Xd[,,2]<-power_diff  # Size dyad
Xd[,,3]<-as.matrix(locFinal) # location shared dyad
Xd[,,4]<-as.matrix(sponsFinal) # sponsorship shared dyad 

dimnames(Xd)[[1]]<-colnames
dimnames(Xd)[[2]]<-colnames
dimnames(Xd)[[3]]<-c("ideol_diff", "powerdiff", "loc", "spons")

# Modeling 
diag(DV_30)<-0 #AMEN won't take NAs
Yreal<-as.matrix(DV_30) # For raw count 
Yrealsqr<-sqrt(Yreal) # For square root - best practice as AME doesn't take a count model, 
#this should better approximate a normal distirubiton

ybin<-Yreal
ybin[ybin>0]<-1


### for sqrt transformed 
fit_infighting_Bivariate_power_sqrt<-ame(Yrealsqr, power_diff, R=1, model="nrm",
                                    symmetric=TRUE,burn=10000,nscan=10000,odens=10)
fit_infighting_Bivariate_ideol_sqrt<-ame(Yrealsqr, ideol_diff, R=1, model="nrm",
                                    symmetric=TRUE,burn=10000,nscan=10000,odens=10)
fit_infighting_Bivariate_spons_sqrt<-ame(Yrealsqr, sponsFinal,  R=1, model="nrm",
                                    symmetric=TRUE,burn=10000,nscan=10000,odens=10)

# sqrt - no ISIS
fit_nodesdyads_noIsis_sqrt<-ame(Yrealsqr, Xd, Xn_b, R=1, model="nrm",
                                     symmetric=TRUE,burn=10000,nscan=10000,odens=10)

# sqrt - WITH ISIS
fit_nodesdyads_Isis_sqrt<-ame(Yrealsqr, Xd, Xn, R=1, model="nrm",
                                symmetric=TRUE,burn=10000,nscan=10000,odens=10)

## binary 
fit_nodesdyads_Isis_bin<-ame(ybin, Xd, Xn, R=1, model="bin",
                              symmetric=TRUE,burn=10000,nscan=10000,odens=10)
# ordinal
fit_nodesdyads_Isis_ord<-ame(Yreal, Xd, Xn, R=1, model="ord", symmetric=TRUE,burn=10000,nscan=10000,odens=10)

summary(fit_infighting_Bivariate_power_sqrt)
summary(fit_infighting_Bivariate_ideol_sqrt)
summary(fit_infighting_Bivariate_spons_sqrt)
summary(fit_nodesdyads_noIsis_sqrt)
summary(fit_nodesdyads_Isis_sqrt)

# To see AME dignostic plots: plot(name of fit object) - these are in Suplementary Materials

summary(fit_infighting_Bivariate_ideol)
summary(fit_infighting_Bivariate_power)
summary(fit_infighting_Bivariate_spons)
summary(fit_nodesdyads_Isis)
summary(fit_nodesdyads_noIsis)

################################
# AME Analysis for 44 Groups - Robsutness Test - Suplemenary Materials #
################################

# Load matrix for 44 groups
DV_44<-read.csv("fullfighting_Matrix_all44_8June2018.csv", stringsAsFactors = F, header = T)

colnames2<-DV_44$X # Taken from row names in this first dataset
rownames(DV_44)<-colnames2 # Make the list rownames for this dataframe
DV_44<-DV_44[,-1] # Subtract duplicate row
colnames(DV_44)<-colnames2 # Also rename colnames (as some will be listed as "X1st" rather than "1st")

ideolVars2<-ideol[,-1] # Get rid of long form names 

# Cut down to only vars of intrest for this analysis 
ideolVars2<-ideolVars2[,-10] # Drop power broker score - alternative measure of power
ideolVars2<-ideolVars2[,-6:-7] # Drop location and state sponsor, will add in later in matrix form
ideolVars2$GroupSize <- as.numeric(gsub(",","",ideolVars2$GroupSize)) # Remove commas
ideolVars2<-ideolVars2[,-7] # Remove group size "full"
idv2<-ideolVars2[match(colnames2, ideolVars2$GroupCode),]
idv2<-na.omit(idv2) # AME won't take NAs 


## Creating Xn and Xd for AME analysis 
b2<-array(idv2[,5]) # Average ideology, could use a compoent of ideology instead
c2<-array(idv2[,6]) # Size
a_a2<-c2/1000 # Putting them on similar scale for rough interp. 

Xn2<-as.array(cbind(b2, a_a2)) # AME takes arrays
dimnames(Xn2)[[1]]<-colnames2 # Labling rows (groups), already ordered/matched above
dimnames(Xn2)[[2]]<-c("averageId", "size") # Lable columns

## Constructing dyadic vars
ad2<-as.data.frame(Xn2) # Generate a dataframe from your array
ideol_diff2<- outer(ad2$averageId, ad2$averageId, "-") # Make a matrix where each square represents the difference in ideology scores of a given dyad
ideol_diff2<-abs(ideol_diff2) # Has to be postive - so take the absolute value 

power_diff2<- outer(ad2$size, ad2$size, "-") # Replicates the above for power
power_diff2<-abs(power_diff2)

# Create an empty array to hold your matrcies (dyad vars) for AME
Xd2<- array(dim=c(nrow(DV_44), nrow(DV_44), 2)) # Must be same dimentions as number of groups (rows and columns) and have enough "slices" to fit all dyad vars  
Xd2[,,1]<- ideol_diff2 # Average ideology dyad
Xd2[,,2]<-power_diff2  # Size dyad

dimnames(Xd2)[[1]]<-colnames2
dimnames(Xd2)[[2]]<-colnames2
dimnames(Xd2)[[3]]<-c("ideol_diff", "powerdiff")


# Modeling 
diag(DV_44)<-0 #AMEN won't take NAs
Yreal2<-as.matrix(DV_44) # For raw count 
Yrealsqr2<-sqrt(Yreal2) # For square root - best practice as AME doesn't take a count model, 
#this should better approximate a normal distirubiton

### for sqrt transformed 
fit_infighting_Bivariate_power_sqrt_44<-ame(Yrealsqr2, power_diff2, R=1, model="nrm", #substitute Yreal2 for non-squareroot transformed
                                         symmetric=TRUE,burn=10000,nscan=10000,odens=10)
fit_infighting_Bivariate_ideol_sqrt_44<-ame(Yrealsqr2, ideol_diff2, R=1, model="nrm",
                                         symmetric=TRUE,burn=10000,nscan=10000,odens=10)

# sqrt - no ISIS
fit_nodesdyads_noIsis_sqrt_44<-ame(Yrealsqr2, Xd2, Xn2, R=1, model="nrm",
                                symmetric=TRUE,burn=10000,nscan=10000,odens=10)

summary(fit_infighting_Bivariate_power_sqrt_44)
summary(fit_infighting_Bivariate_ideol_sqrt_44)
summary(fit_nodesdyads_noIsis_sqrt_44)

# count DV
fit_nodesdyads_noIsis_count_44<-ame(Yreal2, Xd2, Xn2, R=1, model="nrm",
                                   symmetric=TRUE,burn=10000,nscan=10000,odens=10)

# To see AME dignostic plots: plot(name of fit object)


################################
# ERGM Analysis for 30 and 44 Groups - Suplementary Materials #
################################
detach("package:igraph", unload=TRUE)
library(statnet)
library(ergm)
library(ergm.count)

# For 30 groups 

net <- as.network(x = DV_30,  # The network object
                  directed = FALSE, # Specify whether the network is directed
                  loops = FALSE, # do we allow self ties (should not allow them)
                  matrix.type = "adjacency" # the type of input
)

network.vertex.names(net)

## set attributes for nodes 
set.edge.value(net, "Ideol_dif", ideol_diff) 
set.edge.value(net, "power_dif", power_diff)
set.edge.value(net, "shared_spon", sponsFinal)
set.vertex.attribute(net, "size", as.numeric(a_a)) # grab these from above 
set.vertex.attribute(net, "ideol", as.numeric(b))
set.vertex.attribute(net, "spons", as.numeric(a))

m30<-ergm(net ~ edges + 
           absdiff("ideol") +
           absdiff("size"), 
         directed = F)


m30_2<-ergm(net ~ edges+ 
           nodecov("ideol") + 
           nodecov("size") + 
           absdiff("ideol") +
           absdiff("size"),
         directed = F)

par(mfrow = c(2, 2))
plot(gof(m30_2))



# For 44 groups 

net2 <- as.network(x = DV_44,  # The network object
                  directed = FALSE, # Specify whether the network is directed
                  loops = FALSE, # do we allow self ties (should not allow them)
                  matrix.type = "adjacency" # the type of input
)

network.vertex.names(net2)

## set attributes for nodes 
set.edge.value(net2, "Ideol_dif", ideol_diff2) 
set.edge.value(net2, "power_dif", power_diff2)
set.vertex.attribute(net2, "size", as.numeric(a_a2)) 
set.vertex.attribute(net2, "ideol", as.numeric(b2))

m44<-ergm(net2 ~ edges + 
           absdiff("ideol") +
           absdiff("size"), 
         directed = F)


m44_2<-ergm(net2 ~ edges+ 
           nodecov("ideol") + 
           nodecov("size") + 
           absdiff("ideol") +
           absdiff("size"),
         directed = F)

par(mfrow = c(2, 2))
plot(gof(m5))



summary(m30)
summary(m30_2)
summary(m44)
summary(m44_2)




#' This script runs SaTScan on the simulated data, 
#' to evaluate how well it detects clusters of infection risk.
#' @author Amelie Dolfi
#' @Date September 22 2023
#' @Updated : Feb 06 2026
#' @note To use this code, you need to download the Satscan software
#' You also need to adjust the work directory and the name for the parameter file

library(ggplot2)
library(ggforce) # for plotting circles with ggplot
library(RColorBrewer)
library(terra)
library(rsatscan)
library(beepr)



## Function to convert back from the element number, find the place in the matrix
#' From LAUREN WHITE 2018
Rmatrix<-function (vector_loc, nrows){
  col <- floor(vector_loc/nrows) + 1
  row <- vector_loc - (col - 1) * nrows + 1
  return(list(col = col, row = row))
}


# Get the parameters
# Here, I run this on OK. If needed to run on CS or SK, change the parameter and saving folder later on
parameters <- read.csv('Data_input/Parameters_OK.csv') # I run on OK here, same thing for Sk and CS


### Check options at https://rdrr.io/cran/rsatscan/src/R/zzz.R

# an analysis should start by resetting the parameter file for "good style"; this just resets everything to default
invisible(ss.options(reset=T))

# once reset, change parameter options as desired
ss.options(list(CaseFile="/test_cases_r.cas")) # where to find the case file

# various analysis choices (need to specify about not using dates/times otherwise SaTScan gets confused)
# NOTE: Reading a bit more about the exponential versus normal models used in SaTScan, I think you actually want the NORMAL model
ss.options(list(CoordinatesFile="/Coordinates_r.geo", 
                CoordinatesType=0, AnalysisType=1, ModelType=5, ScanAreas=3, 
                PrecisionCaseTimes=0, TimeAggregationUnits=0))
# some options for results
ss.options(list(ResultsFile="/test_r_results.txt", 
                OutputShapefiles="n", MostLikelyClusterEachCentroidASCII="y"))
# can play with max population to evaluate; absolute max is 50% - large values may predispose you to "clusters" that aren't really what we're looking for (see results of the toy analysis below)
ss.options(list(MaxSpatialSizeInPopulationAtRisk=30, MaxSpatialSizeInMaxCirclePopulationFile=30)) 

head(ss.options(), 12)

# write out the parameter file
write.ss.prm(location = "/", filename = "test_r_params", matchout = F) # FYI: help file recommends "matchout = T" 
# I'm not sure why "matchout = T" is preferred - the helpful gives a cryptic "this will likely result in undesirable performance"
# but if "matchout = F," the parameter and output files get the same name, which may be annoying for automation purposes. 
# So...disclaimer/word of caution about that.



##########
nsim = 20 # I have 20 simulations by parameters Index

## Base Habitat DF to get the clusters 
All_hab <- data.frame(y = rep(1:33, 33),
                      x = rep(1:33, each = 33),
                      vec = seq(1, 33*33))


Satscan_all <- data.frame()

for(i in (1:nrow(parameters))){

  para <- parameters[i,]
  
  Clusters <- data.frame()
  
  for(s in (1:nsim)){
    
    if(file.exists('/HabRisk_',i,'_Sim',s,'_Para',para$para,'.csv')){
      Probs <- read.csv('/HabRisk_',i,'_Sim',s,'_Para',para$para,'.csv')
      
      for(st in seq(300,max(Probs$step), by = 600)){
        
        Probs_st <- subset(Probs, step == st )
        
        # Only enter if there are at least 1% of the environment that is infected 
        if(nrow(Probs_st) > 0.01 * (33*33)){
          
          Probs_st$no.cases <- 1 # assign each cell as 1 "case" (SaTScan requires "number of cases per location")
          
          Probs_st <- merge(Probs_st, All_hab, by ='vec', all = T)
          
          Probs_st[is.na(Probs_st$P_inf), 'P_inf'] <- 0
          Probs_st[is.na(Probs_st$no.cases), 'no.cases'] <- 0
          
          names(Probs_st)[c(1,5)] <- c('id','attribute')
          
          # case file includes location id, number of cases at that location (always 1), and the attribute value (infection probability) at that location
          case.file <- Probs_st[,c("id", "no.cases",'attribute')] 
          # coordinate file includes location id and the x and y coordinates for the location
          coord.file <- Probs_st[,c("id", "x", "y")]
          
          # can adjust location and filename as needed, when saving:
          write.cas(case.file, location= "/", filename = paste0("test_cases_r"))
          write.geo(coord.file, location= "/", filename = paste0("Coordinates_r"))
          
          
          #### RUN SATSCAN AND EXTRACT RESULTS ####
          
          ## now run SaTScan; note, with the 10,000 cells here, this takes a few minutes
          test.run <- satscan(prmlocation = "/",
                              prmfilename = "test_r_params",
                              sslocation = "C:/Program Files/SaTScan", # you may need to change this, depending on where SaTScan is saved on your computer
                              ssbatchfilename = "SaTScanBatch64",
                              cleanup = F,
                              verbose = T)
          
          # SatScan outputs a text file with the results, which we can read in:
          f <- readLines("/test_r_results.txt", warn = F)
          
          cr <- grep("Coordinates / radius", f, value = T) # locations and radii of detected clusters
          # coord <- sub('\\).*','',sub('.*\\(', '', cr)) # Coordinates of clusters
          coord <- sub('.*\\(', '', cr)
          coord <- sub('\\).*','',coord)
          num <- grep("Number of cases", f, value = T) # number of cells in cluster
          mi <- grep("Mean inside", f, value = T) # mean inside cluster
          mo <- grep("Mean outside", f, value = T)  # mean outside cluster
          llr <- grep("Log likelihood ratio", f, value = T) # log-likelihood ratio
          pv <- sub('.*: ', '', grep("P-value", f, value = T))  # p-value
          
          ## compile into dataframe
          
          if(length(cr) != 0){
            
            if(length(cr) == length(pv)){
              p.val =  as.numeric(pv)
            }else{p.val = as.numeric(pv[-length(pv)])}
            
            sat.out <- data.frame(clust.id = 1:length(cr),
                                  radius = as.numeric(sub('.*/ ', '', cr)),
                                  x = sub(',.*', '', coord),
                                  y =  sub('.*,', '', coord),
                                  num.cells = as.numeric(sub('.*: ', '', num)),
                                  mean.inside = as.numeric(sub('.*: ', '', mi)),
                                  mean.outside = as.numeric(sub('.*: ', '', mo)),
                                  log.like.ratio = as.numeric(sub('.*: ', '', llr)),
                                  p.value = p.val,
                                  Para = rep(para$para, length(cr)),
                                  Index = rep(i, length(cr)),
                                  Sim = rep(s,length(cr)),
                                  Step = rep(st,length(cr)))
            
          Clusters <- rbind(Clusters, sat.out)
        
          }
        }
      }
    }
  }
  saveRDS(Clusters, paste0('/Satscan_',i,'_Para_',para$para,'.RData'))
  
  Satscan_all <- rbind(Satscan_all, Clusters)
}

# You can save all the data together
saveRDS(Satscan_all, '/Satscan_ALL.RData')

t2 = Sys.time()

t2 - t1

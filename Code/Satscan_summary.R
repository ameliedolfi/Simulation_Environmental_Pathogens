#' @title Analysis of Rsatscan results
#' @author Amelie Dolfi
#' @date  June 1st 2026
#' 

#rm(list = ls())

lapply(c('ggplot2','ggpubr', 'gganimate', 'ggforce',  ### Plotting
         'tidyverse', 'dplyr','tidyr',  ## Anslysing / transforming data
         'sf' ## for the spatial analysis
), library, character.only=TRUE)


### Obtain the number of unique clusters, their maximum, minimum, mean, median risk/area values (mean.inside output from Satscan)
## ant the number of steps with those clusters recorded
# The function arguments are 
#' @param Satscan The dataframe with all the Satscan information

Analysis_satscan <- function(Satscan){
  
  Unique_clusters <- c()
  overlap <- c()
  
  ## If there are more than one recorded cluster, look at the overlaps
  if(nrow(Satscan)!=1){
    
    theta = runif(250, max = 2*pi)
    
    # Create the circle polygon
    Circle1 <- tibble(theta = theta,
                      x = as.numeric(Satscan[1,'x']) + as.numeric(Satscan[1,'radius']) * sin(theta),
                      y = as.numeric(Satscan[1,'y']) + as.numeric(Satscan[1,'radius']) * cos(theta),
                      atan2 = atan2(x, y))
    
    
    Circle <- Circle1 %>% 
      arrange(Circle1) %>% 
      st_as_sf(coords = c("x", "y")) %>% 
      st_combine() %>% 
      st_cast("POLYGON")  %>%
      st_sf()
    
    Circle$Name <- 1
    Circle$Circle_area <- st_area(Circle)
    
    ## For each subsequent clusters, make their circle polygon and calculate the percentage of area that overlap
    for(i in (2:nrow(Satscan))){
      
      Circlei <- tibble(theta = theta,
                        x = Satscan[i,'x'] + Satscan[i,'radius'] * sin(theta),
                        y = Satscan[i,'y'] + Satscan[i,'radius'] * cos(theta),
                        atan2 = atan2(x, y))
      
      Circlei <- Circlei %>% 
        arrange(Circlei) %>% 
        st_as_sf(coords = c("x", "y")) %>% 
        st_combine() %>% 
        st_cast("POLYGON") %>%
        st_sf()
      
      Circlei$Name <- i
      Circlei$Circle_area <- st_area(Circlei)
      
      Circle <- bind_rows(Circle, Circlei)
      
      # Calculate overlap areas from cluster i with all previously calculated clusters
      intersect_pct <- st_intersection(Circle[seq(1,i-1),], Circle[i,] ) %>% 
        mutate(intersect_area = st_area(.)) %>%   # create new column with shape area
        dplyr::select(Name, intersect_area) %>%   # only select columns needed to merge
        st_drop_geometry()  # drop geometry as we don't need it
      
      # Merge by circle name
      Circle_overlap <- merge(Circle[seq(1,i-1),], intersect_pct, by = "Name", all.x = TRUE) %>% 
        mutate(coverage = as.numeric(intersect_area/Circle_area))
      
      overlap <- rbind(overlap,
                       data.frame(Circle1 = seq(1,i-1),
                                  Circle2 = rep(i, length(seq(1,i-1))),
                                  Overlap = Circle_overlap$coverage))
      
    }
    
    Satscan$Lap <- rep(NA)
    Satscan$ID <- seq(1,nrow(Satscan))
    Association <- data.frame(Cluster = seq(1:nrow(Satscan)),
                              Overlap = seq(1:nrow(Satscan)))
    
    ## Then, ensure we get all the unique clusters 
    for(n in (1:(nrow(Circle)-1))){
      for(d in overlap[which(overlap$Overlap > 0.5 & overlap$Circle1 == n),'Circle2']){
        Association[Association$Cluster ==d , 'Overlap'] <- min(n, Association[Association$Cluster == Association[Association$Cluster ==d , 'Overlap'] , 'Overlap'])
      }
    }
    
    ## Finally, obtain the different values desired per unique clusters
    for(i in unique(Association$Overlap)){
      Unique_clusters <- rbind(Unique_clusters, data.frame(Simulation = unique(Satscan$Index),
                                                           Parameter = unique(Satscan$Sim),
                                                           Clust.id = i,
                                                           RA_max = max(Satscan[which(Association$Overlap==i),'mean.inside']), 
                                                           RA_min = min(Satscan[which(Association$Overlap==i),'mean.inside']),
                                                           RA_mean = mean(Satscan[which(Association$Overlap==i),'mean.inside']),
                                                           RA_median = median(Satscan[which(Association$Overlap==i),'mean.inside']),
                                                           RA_numbers = length(Satscan[which(Association$Overlap==i),'mean.inside'])))
      
    }
    
    ## If there was only 1 cluster, record it as unique
  }else{
    
    Unique_clusters <- rbind(Unique_clusters, data.frame(Simulation = unique(Satscan$Index),
                                                         Parameter = unique(Satscan$Sim),
                                                         Clust.id = 1,
                                                         RA_max = Satscan[1,'mean.inside'], 
                                                         RA_min = Satscan[1,'mean.inside'],
                                                         RA_mean = Satscan[1,'mean.inside'],
                                                         RA_median = Satscan[1,'mean.inside'],
                                                         RA_numbers = 1))
    
  }
  return(Unique_clusters)
}


# Load parameters
parameter_OK = read.csv('Data_input/Parameters_OK.csv')
parameter_CS = read.csv('Data_input/Parameters_CS.csv')
parameter_SK = read.csv('Data_input/Parameters_SK.csv')


# Create empty df
All_Satscans_OK <- data.frame()
All_Satscans_CS <- data.frame()
All_Satscans_SK <- data.frame()

# Ensure the name correspond between the rep and simulations
nrep = 5
parameter_num <- rep(seq(1:500), each = nrep)

## Start with CS
Unique_clusters_CS <- data.frame()

## Run the first for loop if doing so on the entire results. Here I use parameter 249 as an example
#for(s in (1:2500)){
for(s in (1241:1245)){
  print(s)
  Satscan_i <- tryCatch({
    read.csv(paste0('Satscan_output_249/Satscan_',s,'_Para_',parameter_num[s],'_CS.csv'))
  }, error = function(err){NA})
  
  
  if(!is.null(nrow(Satscan_i))){
    All_Satscans_CS <- rbind(All_Satscans_CS, Satscan_i)
    
    for(t in unique(Satscan_i$Sim)){
      Satscan_t <- Satscan_i[Satscan_i$Sim == t,]
      Sat_analysis <- Analysis_satscan(Satscan = Satscan_t)
      Unique_clusters_CS <- rbind(Unique_clusters_CS, Sat_analysis)
    }
  }
}

saveRDS(Unique_clusters_CS, 'Satscan_output_249/Unique_clusters_CS.RData')
saveRDS(All_Satscans_CS, 'Satscan_output_249/All_Satscans_CS.RData')


## Now OK 

Unique_clusters_OK <- data.frame()

#for(s in (1:2500)){
for(s in (1241:1245)){
  print(s)
  Satscan_i <- tryCatch({
    read.csv('Satscan_output_249/Satscan_',s,'_Para_',parameter_num[s],'_OK.csv')
  }, error = function(err){NA})
  
  if(!is.null(nrow(Satscan_i))){
    All_Satscans_OK <- rbind(All_Satscans_OK, Satscan_i)
    
    for(t in unique(Satscan_i$Sim)){
      Satscan_t <- Satscan_i[Satscan_i$Sim == t,]
      Sat_analysis <- Analysis_satscan(Satscan = Satscan_t)
      Unique_clusters_OK <- rbind(Unique_clusters_OK, Sat_analysis)
    }
  }
}

saveRDS(Unique_clusters_OK, 'Satscan_output_249/Unique_clusters_OK.RData')
saveRDS(All_Satscans_OK, 'Satscan_output_249/All_Satscans_OK.RData')


## Now SK
Unique_clusters_SK <- data.frame()

t3 <- Sys.time()
for(s in (1:2500)){
  print(s)
  Satscan_i <- tryCatch({
    read.csv(paste0('Satscan_output_249/Satscan_',s,'_Para_',parameter_num[s],'_SK.csv'))
  }, error = function(err){NA})
  
  if(!is.null(nrow(Satscan_i))){
    All_Satscans_SK <- rbind(All_Satscans_SK, Satscan_i)
    
    for(t in unique(Satscan_i$Sim)){
      Satscan_t <- Satscan_i[Satscan_i$Sim == t,]
      Sat_analysis <- Analysis_satscan(Satscan = Satscan_t)
      Unique_clusters_SK <- rbind(Unique_clusters_SK, Sat_analysis)
    }
  }
}

saveRDS(Unique_clusters_SK, 'Satscan_output_249/Unique_clusters_SK.RData')
saveRDS(All_Satscans_SK, 'Satscan_output_249/All_Satscans_SK.RData')


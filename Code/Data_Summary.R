#' @author Amelie Dolfi
#' @Date April 24 2024
#' @Updated October 22 2024
#' @title  Code summarising the model outputs
#' 

### You can test this code on the results for parameter 249, for OK, SK and CS

#rm(list=ls())
##Load packages

lapply(c('ggplot2','ggpubr', 'gganimate', 'ggforce',  ### Plotting
         'tidyverse', 'dplyr','tidyr', 'omnibus',  ## Anslysing / transforming data
         'party', 'randomForest',  'pdp', ## Random Forest, Variable importance, PDF
         'ape',  # Moran
         'sf', ## for the spatial analysis
         'beepr'), library, character.only=TRUE)



## Functions

## transform the environment matrix 
#' From LAUREN WHITE PAPER

Cmatrix<-function (row, col, nrows) { #y,x
  vector_loc <- row + (col - 1) * nrows
  return(vector_loc)
}

## Convert back : from the element number, find, the place in the matrix
#' From LAUREN WHITE PAPER
Rmatrix<-function (vector_loc, nrows){
  col <- floor(vector_loc/nrows) + 1
  row <- vector_loc - (col - 1) * nrows + 1
  return(list(col = col, row = row))
}


### Get all the outputs summarize into one dataset called Numbers
#' This dataset as one row per simulation and calculates
#' The last step
#' The maximum number of infected individuals at one step 
#' The total number of infected individuals individuals at the end of the simulation
#' The number of Susceptible individuals at the end of the simulation
#' The maximum number of reservoir at one step
#' The number of unique environmental cells that were infected at least once
#' The last step at which an infection occured
#' The maximum number of continuous steps where no individuals infected
#' The maximum number of continuous steps where at least one individual is infected
#' The step at which 50% of the population got infected
#' The total number of simulation ran for this parameter run 
#' The number of sucessful simulation for this parameter run

# The arguments to enter are the 
#' @param  parameter parameter dataset, 
#' @param wd the work directory where the outputs are 
#' @param nrep The number of Code repetition of the same paramater
#' @param out any simulation number that is inexistant beacuse it did not work

Infection_visual <- function(parameters, wd, wd_output, nrep, nstep, out=c(), 
                             scenario, last_name = ''){
  Numbers <- c() # initialize the output dataframe
  
  for(p in seq(parameters[1, 'para'], max(parameters$para))){ # For every parameter sets
    print(p)
    para <- parameters[parameters$para==p,][1,] # take the parameter values to open the file
    
    rep <- 0 # Ensure each repetition with the same parameter values have a different repetition number at the end (for example I launched 2 25 simulation run so this ensure that in my final df I have the simulation number from 1 to 50 and not twice 1 to 25)
    
    for(s in seq(para$para*nrep-(nrep-1), para$para*nrep)){ # ensure the simulations with the same parameter values are put together
      if(!s %in% out){
        
        dsum <- read.csv(paste0('Data_output/Sum_',s,'_Para',p,last_name,'.csv'))
        Num_sim = dsum[nrow(dsum),6]
        Total_sim = dsum[nrow(dsum),8]
        dsum <- dsum[-nrow(dsum),]
        
        dsum <- mutate(dsum, Infected = ifelse( Num_I - lag(Num_I) > Num_I_dead - lag(Num_I_dead), TRUE, FALSE))
        
        
        if(Num_sim !=0){
          for(n in (1:Num_sim)){
            
            Test_habrisk <- read.csv(paste0('Data_output/HabRisk_',s,'_Sim',n,'_Para',p,last_name,'.csv'))[-1,]
            
            Inds <- read.csv(paste0('Data_output/Inds_',s,'_Para',p,last_name,'.csv'))[-1,]
            
            rep <- rep +1
            ls <- nrow(subset(dsum, Index == s & sim == n))
            MaxinfT <- max(subset(dsum, Index == s & sim == n)$Num_I, na.rm = T)
            Infection <- subset(dsum, Index == s & sim == n  & step == ls)$Num_I_dead + subset(dsum, Index == s & sim == n & step == ls)$Num_I
            Max_reservoir <- max(subset(dsum, Index == s & sim == n)$Num_res)
            Final_S_population <- subset(dsum, Index == s & sim == n  & step == ls)$Num_S
            Number_cell_infected <- length(unique(Test_habrisk$vec))
            last_inf <- max(subset(dsum, Index == s & sim == n & Num_I !=0)$step)
            Max_period_no_inf <- longRun(subset(dsum, Index == s & sim == n & step < last_inf)$Infected[-1], val = FALSE)
            Max_period_inf <- longRun(subset(dsum, Index == s & sim == n)$Infected[-1] , val = TRUE)
            Step50 = min(subset(dsum, Index == s & sim == n & Num_S <= para$N/2)$step)
            
            Crossed_res_I <- mean(subset(Inds, sim == n & status == 'I' & inf_step != -1)$res_crossed)
            Crossed_res_S <- mean(subset(Inds, sim == n & status == 'S')$res_crossed)
            
            
            Numbers <- rbind(Numbers , data.frame(para = para$para, 
                                                  rep = rep, 
                                                  Sim = s,
                                                  last_step = ls, 
                                                  Max_infection = MaxinfT, 
                                                  Total_infection= Infection,
                                                  Final_S_population = Final_S_population, 
                                                  Max_reservoir_step = Max_reservoir,
                                                  Number_cell_infected = Number_cell_infected,
                                                  last_inf = last_inf, 
                                                  Max_period_no_inf = Max_period_no_inf,
                                                  Max_period_inf = Max_period_inf,
                                                  Step50 = ifelse(Step50== Inf, NA, Step50),
                                                  Crossed_res_I = Crossed_res_I,
                                                  Crossed_res_S = Crossed_res_S,
                                                  Total_sim = Total_sim,
                                                  Total_successful =  Num_sim))
          }
        }
      }
    }
  }
  
  write.csv(Numbers, paste0('Data_output/Numbers', scenario,'.csv'), row.names = FALSE)
  
  
  
  Numbers <- merge(Numbers, unique(parameter_OK[,-11]), by = 'para', all.x=T, all.y =F)
  
  Numbers$Disease_die <- ifelse(Numbers$last_step != nstep & Numbers$Final_S_population != 0, 1,0) # number of simulation with disease die out
  Numbers$Pop_die <- ifelse(Numbers$Final_S_population == 0, 1,0) # number of simulation with population die out
  Numbers$Max_step_Prevalence <- Numbers$Max_infection / Numbers$N # Maximum step prevalence reached in the simulation
  Numbers$Prop_Infected_Final <- Numbers$Total_infection / Numbers$N # Proportion of the population that got infected
  
  ## Number of simulations
  Number_simulations <- unique(subset(Numbers, select = c(para, Sim,Total_sim, Total_successful)))
  
  Number_sim <-  Number_simulations %>% group_by(para) %>%
    summarise(Total_sim = sum(as.numeric(Total_sim)),
              Total_success = sum(Total_successful))
  
  ## Get the average values for each variables
  
  Numbers_means <- Numbers %>% group_by(para) %>%
    summarize(Percent_simulation_Pop_die = (sum(Pop_die)),
              Percent_simulation_Disease_die = (sum(Disease_die)),
              Mean_step_Prevalence = mean(Max_step_Prevalence),
              Mean_Prop_Total_Infected = mean(Prop_Infected_Final),
              Mean_Number_cell_infected = mean(Number_cell_infected),
              Mean_Continuous_Period_Infected = mean(Max_period_inf, na.rm=T),
              Mean_Continuous_Period_No_Infected = mean(Max_period_no_inf, na.rm=T),
              Mean_Step50 = mean(Step50, na.rm= T),
              Mean_Crossed_res_I = mean(Crossed_res_I, na.rm = T),
              Mean_Crossed_res_S = mean(Crossed_res_S, na.rm=T))
  
  
  Numbers_means <- merge(Numbers_means, Number_sim, by ='para', all.x = T)
  
  Numbers_means$Percent_simulation_Pop_die <- Numbers_means$Percent_simulation_Pop_die / Numbers_means$Total_success 
  Numbers_means$Percent_simulation_Disease_die <- Numbers_means$Percent_simulation_Disease_die / Numbers_means$Total_success 
  Numbers_means$Prop_Sucess <- Numbers_means$Total_success / as.numeric(Numbers_means$Total_sim)
  Numbers_means <- merge(Numbers_means, unique(parameters[,-11]), by = 'para')
  
  write.csv(Numbers, paste0('Data_output/Numbers_final_',scenario,'.csv'), row.names = FALSE)
  write.csv(Numbers_means, paste0('Data_output/Numbers_means_',scenario,'.csv'), row.names = FALSE)
  
  return(list(Numbers, Numbers_means))
  
}



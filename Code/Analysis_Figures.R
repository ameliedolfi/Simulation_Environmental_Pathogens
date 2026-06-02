## Code for plots and results 
#' @Date June 1st 2026
#' @title  Analysis and Figures for 
#' "Where curse of the Pharaoh meets landscape of disgust: host behavior and pathogen release mode jointly shape disease dynamics and spatial distribution'

# rm(list = ls())

lapply(c('ggplot2','ggpubr', 'ggforce','ggnewscale',  ### Plotting
         'tidyverse', 'dplyr','tidyr',  ## Analysing / transforming data
         'sf', ## for the spatial analysis
         'rstatix',
         'terra',
         'randomForest','pdp'), library, character.only=TRUE)

## Theme for all the plots
plot_theme <- theme_pubclean()+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 21),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 21),
        title = element_text(size = 22),
        legend.key.size = unit(1.5,'cm'))


#####################################################################
######################## Statistical Analysis #######################
#####################################################################


###" A few functions first
# Function to obtain the parameter landscape info
#' @name Param_land
#' @description
#' This function adds the land as a numeric values, and it adds the info of p
#' and H for each Land type
#' @param df data to add the info on

Param_land <- function(df){
  
  df1 <- df %>% 
    mutate(Landi = case_when(
      Land == 'Land1' ~ 1,
      Land == 'Land2' ~ 2,
      Land == 'Land3' ~ 3,
      Land == 'Land4' ~ 4),
      p = case_when(
        Land == 'Land1' ~ 0.3,
        Land == 'Land2' ~ 0.3,
        Land == 'Land3' ~ 0.7,
        Land == 'Land4' ~ 0.7),
      H = case_when(
        Land == 'Land1' ~ 0.2,
        Land == 'Land2' ~ 0.8,
        Land == 'Land3' ~ 0.2,
        Land == 'Land4' ~ 0.8))
  
  return(df1)
}


## Function to obtain the random forest models
#' @name Compute_rf
#' @description
#' This function computes the random forests and gives out the Variable Importance Analysis info
#' @param df data to run the random forest on
#' @param variable variable to test
#' @param rhoD add ' + rhoD' to the random forest function if death events occurs, otherwise keeps empty ''
#' @param rhoS add ' + rhoS' to the random forest function if shed events occurs, otherwise keeps empty ''
#' @param nt number of trees to run the random forests on. Fixed at 1000 unless stated otherwise
Compute_rf <- function(df, Variable, rhoD ='', rhoS = '' , nt = 1000) {
  
  rf_model <- randomForest(formula = as.formula(paste0(Variable, ' ~ N + p + H + B + B1 + gamma + mu + per ', rhoD, rhoS)), 
                           data=df, ntree=nt, keep.forest=TRUE, 
                           importance=TRUE)
  
  X <- df %>% select(c('N', 'p', 'H', 'B', 'B1', 'gamma', 'mu', 'per', 'rhoD', 'rhoS') )
  y <- df[,Variable]
  
  return(rf_model)
}

#' @name VIA_plot
#' @description
#' This function makes the Variable Importance Analysis plot
#' @param rf_OK random forest results for the OK scenario
#' @param rf_CS random forest results for the CS scenario
#' @param rf_SK random forest results for thtitle to write on the plot
#' @param name_VIA name of the VIA for saving purposes

VIA_Plot <- function(rf_OK, rf_CS, rf_SK, title_VIA, name_VIA){
  
  # Get all the parameters names
  VIA <- data.frame(Parameters = attributes(rf_SK$importance)$dimnames[[1]])
  
  # Had the importance value of the corresponding variable
  VIA$All_OK <- c(rf_OK$importance[seq(1,9),1], NA)
  VIA$All_CS <- c(rf_CS$importance[seq(1,8),1], NA, rf_CS$importance[9,1])
  VIA$All_SK <- rf_SK$importance[,1]
  
  # pivot this for the plot
  VIA <- pivot_longer(VIA, cols = c(2,3,4), names_to = 'Type_All',values_to = 'Importance_all')
  
  VIA <- VIA %>% mutate(Type = case_when(
    Type_All == 'All_OK' ~ 'Obligate-Killer',
    Type_All == 'All_CS' ~ 'Continuous-Shedder',
    Type_All == 'All_SK' ~ 'Shedder-Killer'
  ))
  
  # Visual mark of the 3 best variable per scenario
  Order_parameters <- c("p", "H", "N", "B", "B1",  "per", "gamma", "mu", "rhoD", "rhoS")
  
  VIA <- VIA %>% mutate(Parameters = factor(Parameters, levels = Order_parameters)) %>%
    arrange(Parameters)
  
  
  VIA_plot <- ggplot(VIA, aes(x = Parameters, y =Importance_all, fill = as.factor(Parameters)))+
    facet_wrap(~Type)+
    geom_bar(position = position_dodge2(), stat = 'identity', col ='black')+
    ylab('Variable Importance')+
    scale_fill_manual('Parameter type', 
                      values = c('gold','goldenrod',
                                 'midnightblue','blue2','turquoise','darkblue','purple','deeppink','orchid','magenta'),
                      labels = c('p' = 'p',
                                 'H' = 'H',
                                 'per' = 'per',
                                 'N' = 'N',
                                 'gamma' = expression(gamma),
                                 'rhoD' = expression(rho['D']),
                                 'rhoS' = expression(rho['S']),
                                 'B' =  expression(beta),
                                 'B1' = expression(beta[1]),
                                 'mu' = expression(mu))) +
    scale_x_discrete(labels = c('p' = 'p',
                                'H' = 'H',
                                'per' = 'per',
                                'N' = 'N',
                                'gamma' = expression(gamma),
                                'rhoD' = expression(rho['D']),
                                'rhoS' = expression(rho['S']),
                                'B' =  expression(beta),
                                'B1' = expression(beta[1]),
                                'mu' = expression(mu)))+
    theme_pubclean()+
    plot_theme + 
    theme(strip.text = element_text(size = 20))+
    ggtitle(title_VIA)
  
  return(VIA_plot)
}


#' @name Pdp_function
#' @description
#' This function makes the Partial dependence plots (pdp) depending on the random forests results
#' and the variable of interest
#' @param rf_OK random forest results for the OK scenario
#' @param rf_CS random forest results for the CS scenario
#' @param rf_SK random forest results for thtitle to write on the plot
#' @param sub.data_OK data used for training the pdp on OK
#' @param sub.data_Cs data used for training the pdp on CS
#' @param sub.data_SK data used for training the pdp on SK
#' @param var variable of interest
#' @param title plot title
#' @param xlab label for the variable

Pdp_function <- function(rf_OK,rf_CS,rf_SK, 
                         sub.data_OK, sub.data_CS, sub.data_SK, 
                         var, title, xlab, ylab){
  
  if(var != 'rhoD' & var != 'rhoS'){
    Plot_pdp_OK <- rf_OK %>% 
      partial(pred.var = var , train =  sub.data_OK) %>%
      mutate(Scenario = 'OK')
    Plot_pdp_CS <- rf_CS %>% 
      partial(pred.var = var , train =  sub.data_CS) %>%
      mutate(Scenario = 'CS')
    Plot_pdp_SK <- rf_SK %>% 
      partial(pred.var = var , train =  sub.data_SK) %>%
      mutate(Scenario = 'SK')
    
    plot_PDP <- ggplot()+
      geom_line(data = Plot_pdp_OK, aes(x= Plot_pdp_SK[,1], y = yhat, color = 'OK')) +
      geom_smooth(data = Plot_pdp_OK, aes(x= Plot_pdp_SK[,1], y = yhat, color = 'OK', fill = 'OK'), alpha = 0.2) +
      geom_line(data = Plot_pdp_CS, aes(x= Plot_pdp_SK[,1], y = yhat, color = 'CS')) +
      geom_smooth(data = Plot_pdp_CS, aes(x= Plot_pdp_SK[,1], y = yhat, color = 'CS', fill = 'CS'), alpha = 0.2) +
      geom_line(data = Plot_pdp_SK, aes(x= Plot_pdp_SK[,1], y = yhat, color = 'SK')) +
      geom_smooth(data = Plot_pdp_SK, aes(x= Plot_pdp_SK[,1], y = yhat, color = 'SK', fill = 'SK'), alpha = 0.2) +
      scale_color_manual('Scenario', values= c('red','blue','forestgreen'))+
      scale_fill_manual('Scenario', values= c('red','blue','forestgreen'))+
      xlab(xlab)+
      ylab(ylab)+
      
      ggtitle(title) +
      plot_theme
  }
  
  if(var == 'rhoD' ){
    Plot_pdp_OK <- rf_OK %>% 
      partial(pred.var = var , train =  sub.data_OK) %>%
      mutate(Scenario = 'OK')
    Plot_pdp_SK <- rf_SK %>% 
      partial(pred.var = var , train =  sub.data_SK) %>%
      mutate(Scenario = 'SK')
    
    plot_PDP <- ggplot()+
      geom_line(data = Plot_pdp_OK, aes(x= Plot_pdp_SK[,1], y = yhat, color = 'OK')) +
      geom_smooth(data = Plot_pdp_OK, aes(x= Plot_pdp_SK[,1], y = yhat, color = 'OK', fill = 'OK'), alpha = 0.2) +
      geom_line(data = Plot_pdp_SK, aes(x= Plot_pdp_SK[,1], y = yhat, color = 'SK')) +
      geom_smooth(data = Plot_pdp_SK, aes(x= Plot_pdp_SK[,1], y = yhat, color = 'SK', fill = 'SK'), alpha = 0.2) +
      scale_color_manual('Scenario', values= c('blue','forestgreen'))+
      scale_fill_manual('Scenario', values= c('blue','forestgreen'))+
      
      xlab(xlab)+
      ylab(ylab)+
      
      ggtitle(title) +
      plot_theme
  }
  
  if(var == 'rhoS'){
    Plot_pdp_CS <- rf_CS %>% 
      partial(pred.var = var , train =  sub.data_CS) %>%
      mutate(Scenario = 'CS')
    Plot_pdp_SK <- rf_SK %>% 
      partial(pred.var = var , train =  sub.data_SK) %>%
      mutate(Scenario = 'SK')
    
    plot_PDP <- ggplot()+
      geom_line(data = Plot_pdp_CS, aes(x= Plot_pdp_SK[,1], y = yhat, color = 'CS')) +
      geom_smooth(data = Plot_pdp_CS, aes(x= Plot_pdp_SK[,1], y = yhat, color = 'CS', fill = 'CS'), alpha = 0.2) +
      geom_line(data = Plot_pdp_SK, aes(x= Plot_pdp_SK[,1], y = yhat, color = 'SK')) +
      geom_smooth(data = Plot_pdp_SK, aes(x= Plot_pdp_SK[,1], y = yhat, color = 'SK', fill = 'SK'), alpha = 0.2) +
      scale_color_manual('Scenario', values= c('red','forestgreen'))+
      scale_fill_manual('Scenario', values= c('red','forestgreen'))+
      xlab(xlab)+
      ylab(ylab)+
      
      ggtitle(title) +
      plot_theme
  }
  
  return(plot_PDP)
}

#' @name Pdp_function_all
#' @description
#' This function makes the Partial dependence plots (pdp) for all parameters
#' @param rf_OK random forest results for the OK scenario
#' @param rf_CS random forest results for the CS scenario
#' @param rf_SK random forest results for thtitle to write on the plot
#' @param sub.data_OK data used for training the pdp on OK
#' @param sub.data_Cs data used for training the pdp on CS
#' @param sub.data_SK data used for training the pdp on SK

Pdp_function_all <- function(rf_OK,rf_CS,rf_SK, 
                             sub.data_OK, sub.data_CS, sub.data_SK, ylab){
  

  
  Plot_b1 <- Pdp_function(rf_OK, rf_CS, rf_SK,
                          sub.data_OK, sub.data_CS, sub.data_SK,
                          var = 'B1',
                          title = expression(paste('a. PDP for ', beta[1] )),
                          xlab = expression(paste(beta[1])), ylab)
  
  Plot_b <- Pdp_function(rf_OK, rf_CS, rf_SK,
                         sub.data_OK, sub.data_CS, sub.data_SK,
                         var = 'B',
                         title = expression(paste('b. PDP for ', beta )),
                         xlab = expression(paste(beta)), ylab)
  
  Plot_per <- Pdp_function(rf_OK, rf_CS, rf_SK,
                           sub.data_OK, sub.data_CS, sub.data_SK,
                           var = 'per',
                           title = expression(paste('c. PDP for per' )),
                           xlab = 'per', ylab)
  
  Plot_N <- Pdp_function(rf_OK, rf_CS, rf_SK,
                         sub.data_OK, sub.data_CS, sub.data_SK,
                         var = 'N',
                         title = expression(paste('d. PDP for N' )),
                         xlab = 'N', ylab)
  
  Plot_gamma <- Pdp_function(rf_OK, rf_CS, rf_SK,
                             sub.data_OK, sub.data_CS, sub.data_SK,
                             var = 'gamma',
                             title = expression(paste('e. PDP for ', gamma )),
                             xlab = expression(paste(gamma)), ylab)
  
  Plot_mu <- Pdp_function(rf_OK, rf_CS, rf_SK,
                          sub.data_OK, sub.data_CS, sub.data_SK,
                          var = 'mu',
                          title = expression(paste('f. PDP for ', mu )),
                          xlab = expression(paste(mu)), ylab)
  
  Plot_rhoD <- Pdp_function(rf_OK, rf_CS, rf_SK,
                            sub.data_OK, sub.data_CS, sub.data_SK,
                            var = 'rhoD',
                            title = expression(paste('g. PDP for ', rho['D'] )),
                            xlab = expression(paste(rho['D'])), ylab)
  
  Plot_rhoS <- Pdp_function(rf_OK, rf_CS, rf_SK,
                            sub.data_OK, sub.data_CS, sub.data_SK,
                            var = 'rhoS',
                            title = expression(paste('h. PDP for ', rho['S'] )),
                            xlab = expression(paste(rho['S'])), ylab)
  
  Plot_p <- Pdp_function(rf_OK, rf_CS, rf_SK,
                         sub.data_OK, sub.data_CS, sub.data_SK,
                         var = 'p',
                         title = expression(paste('i. PDP for p')),
                         xlab = 'p', ylab)
  Plot_h <- Pdp_function(rf_OK, rf_CS, rf_SK,
                         sub.data_OK, sub.data_CS, sub.data_SK,
                         var = 'H',
                         title = expression(paste('j. PDP for H')),
                         xlab = 'H', ylab)
  
  plot_PDP_all <-  ggarrange(Plot_b1, Plot_b, Plot_per, Plot_N, Plot_gamma, 
                             Plot_mu, Plot_rhoD, Plot_rhoS,Plot_p, Plot_h, 
                             nrow = 4, ncol = 3, common.legend = TRUE)
  
  
  return(plot_PDP_all)
}


### now the data
## Get the datasets

# Parameter file
para_OK = read.csv('Data_input/Parameters_OK.csv')
para_CS = read.csv('Data_input/Parameters_CS.csv')
para_SK = read.csv('Data_input/Parameters_SK.csv')

## Results summarized at the parameter set level (500 line by file)
Res_para_CS <- read.csv('Data_output/Numbers_means_CS.csv') %>%
  mutate(Scenario = 'CS')
Res_para_SK <- read.csv('Data_output/Numbers_means_SK.csv') %>%
  mutate(Scenario = 'SK')
Res_para_OK <- read.csv('Data_output/Numbers_means_OK.csv') %>%
  mutate(Scenario = 'OK')

Res_para_OK <- Param_land(df = Res_para_OK)
Res_para_SK <- Param_land(df = Res_para_SK)
Res_para_CS <- Param_land(df = Res_para_CS)

Res_para <- bind_rows(Res_para_SK, Res_para_OK, Res_para_CS)


## Results at the simulation level (ALL simulations results here, almost 50,000lines by line) 
Res_sim_CS <- read.csv('Data_output/Numbers_final_CS.csv') %>%
  mutate(Scenario = 'CS')
Res_sim_SK <- read.csv('Data_output/Numbers_final_SK.csv') %>%
  mutate(Scenario = 'SK')
Res_sim_OK <- read.csv('Data_output/Numbers_final_OK.csv') %>%
  mutate(Scenario = 'OK')

Res_sim_OK <- Param_land(df = Res_sim_OK)
Res_sim_SK <- Param_land(df = Res_sim_SK)
Res_sim_CS <- Param_land(df = Res_sim_CS)

Res_sim <- bind_rows(Res_sim_SK, Res_sim_OK, Res_sim_CS) %>%
  mutate(Gamma_plot = case_when(
    gamma < 0.85 ~'Short',
    TRUE ~ 'Long')
  )

# Data for spatial analysis
Unique_clusters_CS <- readRDS('Data_output/Unique_clusters_CS.RData')
Unique_clusters_SK <- readRDS('Data_output/Unique_clusters_SK.RData')
Unique_clusters_OK <- readRDS('Data_output/Unique_clusters_OK.RData')

Num_unique <- bind_rows(Unique_clusters_CS,Unique_clusters_OK, Unique_clusters_SK) %>%
  group_by(para2, Sim, Scenario) %>%
  summarise(Unique_clust = n(),
            .groups = 'drop') %>% 
  left_join(para_SK  , by = 'para2')
Num_unique <- as.data.frame(Param_land(Num_unique))

Unique_clusters_CS_p <- Unique_clusters_CS %>% 
  left_join(para_CS, by = 'para2')
Unique_clusters_CS_p <- Param_land(Unique_clusters_CS_p)

Unique_clusters_SK_p <- Unique_clusters_SK %>% 
  left_join(para_SK, by = 'para2' )
Unique_clusters_SK_p <- Param_land(Unique_clusters_SK_p)

Unique_clusters_OK_p <- Unique_clusters_OK %>% 
  left_join(para_OK, by = 'para2') 
Unique_clusters_OK_p <- Param_land(Unique_clusters_OK_p)
 
Unique_cluster_p <- bind_rows(Unique_clusters_CS_p,Unique_clusters_OK_p, Unique_clusters_SK_p)

#########################################################################
###################### Statistical analysis #############################
#########################################################################

#################################################
## Successfull outbreaks (+5 infected hosts)  ##
################################################
nrow(Res_para[which(Res_para$Total_success < 100 ),])
## 28 simulations failed out of 1500 (i.e. have less than 100 successful simulations)
Res_para[which(Res_para$Total_success < 100), c(1, 26)]
# 6 from OK
# 19 from CS
# 3 from SK
## Proportion of success ## 
summary(Res_para[, 'Prop_Sucess']) ; sd(Res_para[, 'Prop_Sucess'], na.rm = T)

## Number of encounters before infection ##
summary(Res_para[Res_para$Scenario == 'OK', 'Prop_Sucess']) ; sd(Res_para[Res_para$Scenario == 'OK', 'Prop_Sucess'], na.rm = T)
# Number of encounters , average of 68 % (+- 27%)
summary(Res_para[Res_para$Scenario == 'CS', 'Prop_Sucess']) ; sd(Res_para[Res_para$Scenario == 'CS', 'Prop_Sucess'], na.rm = T)
# Number of encounters , average of 47% (+- 29%)
summary(Res_para[Res_para$Scenario == 'SK', 'Prop_Sucess']) ; sd(Res_para[Res_para$Scenario == 'SK', 'Prop_Sucess'], na.rm = T)
# Number of encounters , average of 76% (+- 24)

# Wilcoxon
Signif_continuous_PropSuccess <- Res_para %>%
  pairwise_wilcox_test(Prop_Sucess ~ Scenario, p.adjust.method = "bonferroni")
# All significantly different from each others

## SK > OK and CS

set.seed(1234)
res_PropSucess_OK <- Compute_rf(df = Res_para_OK, Variable = 'Prop_Sucess', rhoD = ' + rhoD')
res_PropSucess_CS <- Compute_rf(df = Res_para_CS, Variable = 'Prop_Sucess', rhoS = ' + rhoS')
res_PropSucess_SK <- Compute_rf(df = Res_para_SK, Variable = 'Prop_Sucess', rhoD = ' + rhoD',  rhoS = ' + rhoS')

SI2a <- VIA_Plot(rf_OK = res_PropSucess_OK,
                 rf_CS = res_PropSucess_CS,
                 rf_SK = res_PropSucess_SK,
                 title_VIA = 'a. Proportion of simulations with a successful outbreak',
                 name_VIA = 'Prop_Sucess')

####################################
## Depletion of susceptible hosts ##
####################################
summary(Res_para[,'Percent_simulation_Pop_die']); sd(Res_para[,'Percent_simulation_Pop_die'])
# Mean of 9.6 % of simulation
summary(Res_para[Res_para$Scenario == 'OK','Percent_simulation_Pop_die']); sd(Res_para[Res_para$Scenario == 'OK','Percent_simulation_Pop_die'])
# Mean of 7.5 %
summary(Res_para[Res_para$Scenario == 'CS','Percent_simulation_Pop_die']); sd(Res_para[Res_para$Scenario == 'CS','Percent_simulation_Pop_die'])
# Mean of 6.0 %
summary(Res_para[Res_para$Scenario == 'SK','Percent_simulation_Pop_die']); sd(Res_para[Res_para$Scenario == 'SK','Percent_simulation_Pop_die'])
# Mean of 15.2 %

# Steps when no susceptible individuals left
summary(Res_sim[Res_sim$Pop_die == 1 & Res_sim$Scenario == 'OK',]$last_step); sd(Res_sim[Res_sim$Pop_die == 1 & Res_sim$Scenario == 'OK',]$last_step)
# Mean 4444 +- 2278 / Max 9983
summary(Res_sim[Res_sim$Pop_die == 1 & Res_sim$Scenario == 'CS',]$last_step); sd(Res_sim[Res_sim$Pop_die == 1 & Res_sim$Scenario == 'CS',]$last_step)
# Mean = 4611 +_ 2431 / Max = 9994
summary(Res_sim[Res_sim$Pop_die == 1 & Res_sim$Scenario == 'SK',]$last_step); sd(Res_sim[Res_sim$Pop_die == 1 & Res_sim$Scenario == 'SK',]$last_step)
# Mean = 3654 +_ 2208 / Max = 10000

# Wilcoxon test
Signif_Pop_Die <- Res_para %>%
  pairwise_wilcox_test(Percent_simulation_Pop_die ~ Scenario, p.adjust.method = "bonferroni")
# SK significantly different from OK and CS
# OK and CS not significantly different

set.seed(1234)
res_Perc_PopDie_OK <- Compute_rf(df = Res_para_OK, Variable = 'Percent_simulation_Pop_die', rhoD = ' + rhoD')
res_Perc_PopDie_CS <- Compute_rf(df = Res_para_CS, Variable = 'Percent_simulation_Pop_die', rhoS = ' + rhoS')
res_Perc_PopDie_SK <- Compute_rf(df = Res_para_SK, Variable = 'Percent_simulation_Pop_die', rhoD = ' + rhoD',  rhoS = ' + rhoS')

SI2b <- VIA_Plot(rf_OK = res_Perc_PopDie_OK,
                 rf_CS = res_Perc_PopDie_CS,
                 rf_SK = res_Perc_PopDie_SK,
                 title_VIA = 'b. Proportion of simulation with susceptible host depletion',
                 name_VIA = 'Perc_PopDie')

##########################
### Pathogen faded-out ##
#########################
summary(Res_para[,'Percent_simulation_Disease_die'])
# Mean of 36.4 %
summary(Res_para[Res_para$Scenario == 'OK','Percent_simulation_Disease_die']); sd(Res_para[Res_para$Scenario == 'OK','Percent_simulation_Disease_die'])
# Mean of 30.4 %
summary(Res_para[Res_para$Scenario == 'CS','Percent_simulation_Disease_die']); sd(Res_para[Res_para$Scenario == 'CS','Percent_simulation_Disease_die'])
# Mean of 48.4 %
summary(Res_para[Res_para$Scenario == 'SK','Percent_simulation_Disease_die']); sd(Res_para[Res_para$Scenario == 'SK','Percent_simulation_Disease_die'])
# mean of 30.3 %

summary(Res_sim[Res_sim$Disease_die == 1 & Res_sim$Scenario == 'OK',]$last_step); sd(Res_sim[Res_sim$Disease_die == 1 & Res_sim$Scenario == 'OK',]$last_step)
# Mean = 6817 +- 1928 / Max = 9999
summary(Res_sim[Res_sim$Disease_die == 1 & Res_sim$Scenario == 'CS',]$last_step); sd(Res_sim[Res_sim$Disease_die == 1 & Res_sim$Scenario == 'CS',]$last_step)
# Mean = 5643 +- 2262 / Max = 9999
summary(Res_sim[Res_sim$Disease_die == 1 & Res_sim$Scenario == 'SK',]$last_step); sd(Res_sim[Res_sim$Disease_die == 1 & Res_sim$Scenario == 'SK',]$last_step)
# mean = 7047 +- 1816 / Max = 9999

# Wilcoxon test
Signif_Disease_die <- Res_para %>%
  pairwise_wilcox_test(Percent_simulation_Disease_die ~ Scenario, p.adjust.method = "bonferroni")
# CS significantly different from OK and SK
# OK and SK not significantly different

set.seed(1234)
res_Perc_DiseaseDie_OK <- Compute_rf(df = Res_para_OK, Variable = 'Percent_simulation_Disease_die', rhoD = ' + rhoD')
res_Perc_DiseaseDie_CS <- Compute_rf(df = Res_para_CS, Variable = 'Percent_simulation_Disease_die', rhoS = ' + rhoS')
res_Perc_DiseaseDie_SK <- Compute_rf(df = Res_para_SK, Variable = 'Percent_simulation_Disease_die', rhoD = ' + rhoD',  rhoS = ' + rhoS')

SI2c <- VIA_Plot(rf_OK = res_Perc_DiseaseDie_OK,
                 rf_CS = res_Perc_DiseaseDie_CS,
                 rf_SK = res_Perc_DiseaseDie_SK,
                 title_VIA = 'c. Proportion of simulation with pathogen fade-out',
                 name_VIA = 'Perc_DzDie')


########################################
## Proportion of infected individuals ##
########################################
summary(Res_sim[Res_sim$Scenario == 'OK', 'Prop_Infected_Final']) ; sd(Res_sim[Res_sim$Scenario == 'OK', 'Prop_Infected_Final'])
# for OK scenario, the proportion of population infected at the end is 60% (+- 37%)
summary(Res_sim[Res_sim$Scenario == 'CS', 'Prop_Infected_Final']) ; sd(Res_sim[Res_sim$Scenario == 'CS', 'Prop_Infected_Final'])
# for CS scenario, the proportion of population infected at the end is 47% (+- 40%)
summary(Res_sim[Res_sim$Scenario == 'SK', 'Prop_Infected_Final']) ; sd(Res_sim[Res_sim$Scenario == 'SK', 'Prop_Infected_Final'])
# for SK scenario, the proportion of population infected at the end is 69% (+- 34%)

# Wilcoxon test
Signif_Prop_Inf <- Res_sim %>%
  pairwise_wilcox_test(Prop_Infected_Final ~ Scenario, p.adjust.method = "bonferroni")
# All significantly different from each others

set.seed(1234)
res_Perc_PropInd_OK <- Compute_rf(df = Res_sim_OK, Variable = 'Prop_Infected_Final', rhoD = ' + rhoD')
res_Perc_PropInd_CS <- Compute_rf(df = Res_sim_CS, Variable = 'Prop_Infected_Final', rhoS = ' + rhoS')
res_Perc_PropInd_SK <- Compute_rf(df = Res_sim_SK, Variable = 'Prop_Infected_Final', rhoD = ' + rhoD',  rhoS = ' + rhoS')

SI2d <- VIA_Plot(rf_OK = res_Perc_PropInd_OK,
                 rf_CS = res_Perc_PropInd_CS,
                 rf_SK = res_Perc_PropInd_SK,
                 title_VIA = 'd. Proportion of individuals infected per simulation',
                 name_VIA = 'Perc_PropInD')

###############################
## Period without infection ##
##############################
summary(Res_sim$Max_period_no_inf) ; sd(Res_sim$ Max_period_no_inf, na.rm = T)
# Mean of 1719
summary(Res_sim[Res_sim$Scenario == 'OK', 'Max_period_no_inf']) ; sd(Res_sim[Res_sim$Scenario == 'OK', 'Max_period_no_inf'], na.rm = T)
# Max number of continuous steps without infections, average of 1736 (+- 1209; max 9706)
summary(Res_sim[Res_sim$Scenario == 'CS', 'Max_period_no_inf']) ; sd(Res_sim[Res_sim$Scenario == 'CS', 'Max_period_no_inf'], na.rm = T)
# Max number of continuous steps without infections, average of 1659 (+- 1059; max 9760)
summary(Res_sim[Res_sim$Scenario == 'SK', 'Max_period_no_inf']) ; sd(Res_sim[Res_sim$Scenario == 'SK', 'Max_period_no_inf'], na.rm = T)
# Max number of continuous steps without infections, average of 1760 (+- 1274; max 9826)

# Wilcoxon test
Signif_continuous_Noinf_step <- Res_sim %>%
  pairwise_wilcox_test(Max_period_no_inf ~ Scenario, p.adjust.method = "bonferroni")
# CS significantly different from OK and SK
# OK and SK not significantly different from each others
#  CS < SK and OK

set.seed(1234)
res_Perc_MaxPeirod_OK <- Compute_rf(df = Res_sim_OK, Variable = 'Max_period_no_inf', rhoD = ' + rhoD')
res_Perc_MaxPeirod_CS <- Compute_rf(df = Res_sim_CS, Variable = 'Max_period_no_inf', rhoS = ' + rhoS')
res_Perc_MaxPeirod_SK <- Compute_rf(df = Res_sim_SK, Variable = 'Max_period_no_inf', rhoD = ' + rhoD',  rhoS = ' + rhoS')

SI2e <- VIA_Plot(rf_OK = res_Perc_MaxPeirod_OK,
                 rf_CS = res_Perc_MaxPeirod_CS,
                 rf_SK = res_Perc_MaxPeirod_SK,
                 title_VIA = 'e. Maximum period between consecutive infection',
                 name_VIA = 'Max_period_no_inf')


###########################################
## Number of encounters before infection ##
###########################################
summary(Res_para[Res_para$Scenario == 'OK', 'Mean_Crossed_res_I']) ; sd(Res_para[Res_para$Scenario == 'OK', 'Mean_Crossed_res_I'], na.rm = T)
# Number of encounters , average of 62 (+- 51; max 483, min = 1.6)
summary(Res_para[Res_para$Scenario == 'CS', 'Mean_Crossed_res_I']) ; sd(Res_para[Res_para$Scenario == 'CS', 'Mean_Crossed_res_I'], na.rm = T)
# Number of encounters , average of 131 (+- 102; max 932, min = 3.4)
summary(Res_para[Res_para$Scenario == 'SK', 'Mean_Crossed_res_I']) ; sd(Res_para[Res_para$Scenario == 'SK', 'Mean_Crossed_res_I'], na.rm = T)
# Number of encounters , average of 87 (+- 67; max 468, min = 3)

# Wilcoxon
Signif_continuous_Noinf_step <- Res_para %>%
  pairwise_wilcox_test(Mean_Crossed_res_I ~ Scenario, p.adjust.method = "bonferroni")
# All significantly different from each others
# OK < SK < CS

set.seed(1234)
res_MeanInteract_OK <- Compute_rf(df = Res_para_OK, Variable = 'Mean_Crossed_res_I', rhoD = ' + rhoD')
res_MeanInteract_CS <- Compute_rf(df = Res_para_CS, Variable = 'Mean_Crossed_res_I', rhoS = ' + rhoS')
res_MeanInteract_SK <- Compute_rf(df = Res_para_SK, Variable = 'Mean_Crossed_res_I', rhoD = ' + rhoD',  rhoS = ' + rhoS')

SI2f <- VIA_Plot(rf_OK = res_MeanInteract_OK,
                 rf_CS = res_MeanInteract_CS,
                 rf_SK = res_MeanInteract_SK,
                 title_VIA = 'f. Mean number of reservoir interaction before infection',
                 name_VIA = 'Mean_Crossed_res_I')


##################################################
## Proportion of infectious environmental cells ##
##################################################
summary(Res_para$Mean_Number_cell_infected)/ 1089 ; sd(Res_para$Mean_Number_cell_infected) / 1089
# Mean of 17 %
summary(Res_para[Res_para$Scenario == 'OK', 'Mean_Number_cell_infected']) / 1089 ; sd(Res_para[Res_para$Scenario == 'OK', 'Mean_Number_cell_infected']) / 1089
# for OK scenario, the average % of environment infected at some point is 10% (+- 7%; max = 43%)
summary(Res_para[Res_para$Scenario == 'CS', 'Mean_Number_cell_infected']) / 1089 ; sd(Res_para[Res_para$Scenario == 'CS', 'Mean_Number_cell_infected']) / 1089
# for CS scenario, the average % of environment infected at some point is 16% (+- 17%; max = 99.4%)
summary(Res_para[Res_para$Scenario == 'SK', 'Mean_Number_cell_infected']) / 1089 ; sd(Res_para[Res_para$Scenario == 'SK', 'Mean_Number_cell_infected']) / 1089
# for SK scenario, the average % of environment infected at some point is 25% (+- 20%; max = 99.6%)

# Wilcoxon
Signif_Prop_cell <- Res_para %>%
  pairwise_wilcox_test(Mean_Number_cell_infected ~ Scenario, p.adjust.method = "bonferroni")
# All significantly different fro; each others 
## OK < CS < SK

set.seed(1234)
res_Perc_PropCell_OK <- Compute_rf(df = Res_para_OK, Variable = 'Mean_Number_cell_infected', rhoD = ' + rhoD')
res_Perc_PropCell_CS <- Compute_rf(df = Res_para_CS, Variable = 'Mean_Number_cell_infected', rhoS = ' + rhoS')
res_Perc_PropCell_SK <- Compute_rf(df = Res_para_SK, Variable = 'Mean_Number_cell_infected', rhoD = ' + rhoD',  rhoS = ' + rhoS')

SI2g <- VIA_Plot(rf_OK = res_Perc_PropCell_OK,
                 rf_CS = res_Perc_PropCell_CS,
                 rf_SK = res_Perc_PropCell_SK,
                 title_VIA = 'g. Proportion of cumulative environmental pathogen reservoirs',
                 name_VIA = 'Perc_PropCell')

#######################
## Mean cluster peak ##
#######################
summary(Unique_cluster_p$RA_mean) ; sd(Unique_cluster_p$ RA_mean, na.rm = T)
summary(Unique_cluster_p[Unique_cluster_p$Scenario == 'OK', 'RA_mean']) ; sd(Unique_cluster_p[Unique_cluster_p$Scenario == 'OK', 'RA_mean'], na.rm = T)
# Infection risk peak mean : 0.034 +_ 0.032 / max = 0.51
summary(Unique_cluster_p[Unique_cluster_p$Scenario == 'CS', 'RA_mean']) ; sd(Unique_cluster_p[Unique_cluster_p$Scenario == 'CS', 'RA_mean'], na.rm = T)
# Infection risk peak mean : 0.0097 +_ 0.011 / Max = 0.33
summary(Unique_cluster_p[Unique_cluster_p$Scenario == 'SK', 'RA_mean']) ; sd(Unique_cluster_p[Unique_cluster_p$Scenario == 'SK', 'RA_mean'], na.rm = T)
# Infection risk peak mean : 0.027 +_ 0.031 / Max = 0.69

# Wilcoxon
Signif_mean_RA <- Unique_cluster_p %>%
  pairwise_wilcox_test(RA_mean ~ Scenario, p.adjust.method = "bonferroni")
# ALL SIGNIFCANT : CS < SK < OK
# CS < SK < OK

set.seed(1234)
res_MeanCluster_OK <- Compute_rf(df = Unique_cluster_p[Unique_cluster_p$Scenario == 'OK',], Variable = 'RA_mean', rhoD = ' + rhoD', nt = 250)
res_MeanCluster_CS <- Compute_rf(df = Unique_cluster_p[Unique_cluster_p$Scenario == 'CS',], Variable = 'RA_mean', rhoS = ' + rhoS', nt = 250)
res_MeanCluster_SK <- Compute_rf(df = Unique_cluster_p[Unique_cluster_p$Scenario == 'SK',], Variable = 'RA_mean', rhoD = ' + rhoD',  rhoS = ' + rhoS', nt = 250)

SI2h <- VIA_Plot(rf_OK = res_MeanCluster_OK,
                 rf_CS = res_MeanCluster_CS,
                 rf_SK = res_MeanCluster_SK,
                 title_VIA = 'h. Mean cluster infectivity peak',
                 name_VIA = 'RA_mean')


######################################
## Number of reservoirs in clusters ##
######################################

summary(Unique_cluster_p$Size_mean) ; sd(Unique_cluster_p$ Size_mean, na.rm = T)
summary(Unique_cluster_p[Unique_cluster_p$Scenario == 'OK', 'Size_mean']) ; sd(Unique_cluster_p[Unique_cluster_p$Scenario == 'OK', 'Size_mean'], na.rm = T)
# Mean number of reservoir in a cluster :  Median :9 / Mean sd 15 +_ 14.7 / Max = 127
summary(Unique_cluster_p[Unique_cluster_p$Scenario == 'CS', 'Size_mean']) ; sd(Unique_cluster_p[Unique_cluster_p$Scenario == 'CS', 'Size_mean'], na.rm = T)
# Mean number of reservoir in a cluster :  Median : 23 / Mean sd 32 +_ 32 / Max = 323
summary(Unique_cluster_p[Unique_cluster_p$Scenario == 'SK', 'Size_mean']) ; sd(Unique_cluster_p[Unique_cluster_p$Scenario == 'SK', 'Size_mean'], na.rm = T)
# Mean number of reservoir in a cluster : Median : 29 / Mean sd 42 +_ 44 / Max = 326

# Wilcoxon
Signif_size_mean <- Unique_cluster_p %>%
  pairwise_wilcox_test(Size_mean ~ Scenario, p.adjust.method = "bonferroni")
# ALL SIGNIFICANT : OK < CS < SK

set.seed(1234)
res_Size_mean_OK <- Compute_rf(df = Unique_cluster_p[Unique_cluster_p$Scenario == 'OK',], Variable = 'Size_mean', rhoD = ' + rhoD', nt = 250)
res_Size_mean_CS <- Compute_rf(df = Unique_cluster_p[Unique_cluster_p$Scenario == 'CS',], Variable = 'Size_mean', rhoS = ' + rhoS', nt = 250)
res_Size_mean_SK <- Compute_rf(df = Unique_cluster_p[Unique_cluster_p$Scenario == 'SK',], Variable = 'Size_mean', rhoD = ' + rhoD',  rhoS = ' + rhoS', nt = 250)

SI2i <- VIA_Plot(rf_OK = res_Size_mean_OK,
                 rf_CS = res_Size_mean_CS,
                 rf_SK = res_Size_mean_SK,
                 title_VIA = 'i. Mean number of reservoirs within a cluster',
                 name_VIA = 'Size_mean')

##############################
## Number of unique cluster ##
##############################

summary(Num_unique$Unique_clust) ; sd(Num_unique$ Unique_clust, na.rm = T)
summary(Num_unique[Num_unique$Scenario == 'OK',]$Unique_clust) ; sd(Num_unique[Num_unique$Scenario == 'OK',]$Unique_clust, na.rm = T)
# Number of unique clusters : 1.71 +_ 0.93 / Max = 8
summary(Num_unique[Num_unique$Scenario == 'CS',]$Unique_clust) ; sd(Num_unique[Num_unique$Scenario == 'CS',]$Unique_clust, na.rm = T)
# Number of unique clusters  : 4.39 +_ 2.78 / Max = 23
summary(Num_unique[Num_unique$Scenario == 'SK',]$Unique_clust) ; sd(Num_unique[Num_unique$Scenario == 'SK',]$Unique_clust, na.rm = T)
# Number of unique clusters  : 2.88 +_ 1.76 / Max = 15

# Wilcoxon
Signif_Unique_clust <- Num_unique %>%
  as.data.frame() %>%
  pairwise_wilcox_test(Unique_clust ~ Scenario, p.adjust.method = "bonferroni")

## ALL SIGNIFICANT : SK < OK < CS

set.seed(1234)
res_Unique_clust_OK <- Compute_rf(df = Num_unique[Num_unique$Scenario == 'OK',], Variable = 'Unique_clust', rhoD = ' + rhoD')
res_Unique_clust_CS <- Compute_rf(df = Num_unique[Num_unique$Scenario == 'CS',], Variable = 'Unique_clust', rhoS = ' + rhoS')
res_Unique_clust_SK <- Compute_rf(df = Num_unique[Num_unique$Scenario == 'SK',], Variable = 'Unique_clust', rhoD = ' + rhoD',  rhoS = ' + rhoS')


SI2j <- VIA_Plot(rf_OK = res_Unique_clust_OK,
                 rf_CS = res_Unique_clust_CS,
                 rf_SK = res_Unique_clust_SK,
                 title_VIA = 'j. Number of unique clusters created over a simulation',
                 name_VIA = 'Unique_clust')



#####################################################################
############################## Figures ##############################
#####################################################################

############################
######## MS figures ########
############################


##############
## Figure 1 ##
##############
# Background landscape 

# Get the Landscape as a df

#' @name Plot_land
#' @details
#' This function allows to transform the csv Land data into a usable 
#' data frame
#' @param Land_p Name of the file 
#' @param wd workdirectory to find the file
Plot_Land <- function(Land_p){
  Land <-  rast(as.matrix(read.csv(paste0('Data_input/',Land_p,'.csv')))) # convert our matrix to a raster
  
  # Convert raster → data frame
  Land <- as.data.frame(Land, xy = TRUE)
  
  names(Land)[3] <- 'Resource'
  Land$Resource <- ifelse(Land$Resource == 0, 'Low','High')
  Land$xloc <- Land$x + 0.5
  Land$yloc <- Land$y + 0.5
  
  return(Land)
}

# Get the landscape df
Land1 <- Plot_Land('Land1')
Land2 <- Plot_Land('Land2')
Land3 <- Plot_Land('Land3')
Land4 <- Plot_Land('Land4')


## Make background plot for figure 1 :

pLand1 <- ggplot(Land1)+
  geom_tile(aes(x=x, y=y, fill = Resource)) +
  scale_fill_manual(name = 'Resource', values = c('Low'="navajowhite2", 'High'= "forestgreen"))+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16),
        panel.background = element_blank())


## plot with the individuals on it
set.seed(1234) # Put a seed for reproducibility
random_ind <- data.frame(x = sample(0.5:32.5, size = 200, replace = TRUE), y = sample(0.5:32.5,200, replace = TRUE),
                         Statut = c('I',rep('S',199)))

pLand1 +
  geom_point(data = random_ind, aes(x=x, y=y, color = Statut), size = 4)+
  scale_color_manual(values = c('red','blue'))

ggsave('Figures/Land1.png', width = 9, height = 8)


# Plot with the reservoirs / clusters
# I take one of the simulation to plot the reservoirs
Cluster <- readRDS('/Data_output/Satscan_6_Para_2_OK.RData')
Cluster$x <- as.integer(Cluster$x)
Cluster$y <- as.integer(Cluster$y)
Cluster <- Cluster[Cluster$p.value < 0.05 & Cluster$Sim == 19 & Cluster$Step == 5100,]

Hab <- read.csv('/Data_output/HabRisk_6_Sim19_Para2.csv')
Hab$xloc <- as.integer(Hab$xloc)
Hab$yloc <- as.integer(Hab$yloc)
Hab$P_inf <- as.numeric(Hab$P_inf)

Hab_s <- Hab[Hab$step == 5100,]

Land_pinf <- merge(Hab_s, Land1, by = c('xloc','yloc'), all = T)
Land_pinf$P_inf <- ifelse(is.na(Land_pinf$P_inf),0, round(as.numeric(Land_pinf$P_inf),5))

pLand1+
  geom_point(data=Land_pinf[Land_pinf$P_inf!=0,], aes(x=x,y=y, col = P_inf), size = 4)+
  scale_color_gradient('Infectivity', low = 'grey50', high = 'grey10', limits = c(0,  round(max(summary(Land_pinf$P_inf)),2)))+
  geom_circle(data = Cluster, aes(x0 = x, y0 = y, r = radius), color = "red") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave('Figures/Land1_pinf.png', width = 9, height = 8)


# Plot the disease dynaimc:

SK_dynamic <- readRDS('/Data_output/All_SK_time_dynamic_Dec2025.RData')
SK_dynamic <- merge(SK_dynamic, unique(para_SK[,-11]), by = 'para')


One_ex <- SK_dynamic %>%
  filter(para == 1) %>%
  group_by(step) %>%
  mutate(ymin = ifelse(mean_p - sd_p <0, 0, mean_p - sd_p),
         ymax = ifelse(mean_p + sd_p > 100, 100, mean_p + sd_p),
         yminS = ifelse(mean_pS - sd_pS <0, 0, mean_pS - sd_pS),
         ymaxS = ifelse(mean_pS + sd_pS > 100, 100, mean_pS + sd_pS))

ggplot(One_ex, aes(x=step, y = mean_p))+
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.1, lty = 2, linewidth = 0.4, fill = 'red') + 
  geom_ribbon(aes(ymin = yminS, ymax = ymaxS), alpha = 0.1, lty = 2, linewidth = 0.4, fill = 'blue') + 
  geom_line(linewidth = 0.7, aes(color = 'Infected'))+
  geom_line(aes(y = mean_pS, color = 'Susceptible'), linewidth = 0.7)+
  scale_color_manual('Status', values = c('Infected' = 'red', 'Susceptible'='blue'))+
  scale_y_continuous('% individuals', breaks = seq(0,100, by = 25), labels = seq(0,100, by = 25), limits = c(0,100))+
  theme_pubclean()+
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))

ggsave('Figures/F1_DD.png', width = 7, height = 5)


##############
## Figure 2 ##
##############
# Heatmap for the percentage of simulations with disease fade out (a) 
# and population die out (b)
# depending on parameter mu (low: 0.0002-0.00085, med: 0.00085-0.0015, high: 0.0015-0.002)


Res_para$mu_plot <- factor(ifelse(Res_para$mu <= 0.00085, 'low', 
                                    ifelse(Res_para$mu > 0.0015, 'high', 'medium')),
                             levels = c('low','medium','high'))


#a.
df_popdieout <- Res_para %>%
  group_by(Scenario, mu_plot) %>%
  summarise(value = mean(Percent_simulation_Pop_die, na.rm = TRUE),
            .groups ='drop')

Plot2a <-ggplot(df_popdieout, aes(x=Scenario, y = mu_plot, fill=value)) +
  geom_tile(col = 'black')+
  ylab('Pathogen decay rate')+
  scale_fill_gradient2('Proportion of simulation with\n susceptible host depletion', low = 'turquoise3', mid = 'yellow', high = 'red', midpoint = 0.5, limit = c(0,1))+
  ggtitle('a. ')+
  plot_theme

Plot2a

#b.
df_dzfadeout <- Res_para %>%
  group_by(Scenario, mu_plot) %>%
  summarise(value = mean(Percent_simulation_Disease_die, na.rm = TRUE))


Plot2b <- ggplot(df_dzfadeout, aes(x=Scenario, y = mu_plot, fill=value)) +
  geom_tile(col = 'black')+
  ylab('Pathogen decay rate')+
  scale_fill_gradient2('Proportion of simulation\nwith disease fade-out', low = 'turquoise3', mid = 'yellow', high = 'red', midpoint = 0.5, limit = c(0,1))+
  plot_theme +
  ggtitle('b.')

Plot2b

ggarrange(Plot2a, Plot2b, nrow = 1)

ggsave('/Figures/Figure2_HeatMap.png', width = 18, height = 10)


##############
## Figure 3 ##
##############

OK_dynamic <- readRDS('D:/Amelie/These/MS_Simulation/Code et analysis/GitHubVersion//Data_output/All_OK_time_dynamic_Dec2025.RData')
OK_dynamic <- merge(OK_dynamic, unique(para_OK[,-11]), by = 'para')

CS_dynamic <- readRDS('/Data_output/All_CS_time_dynamic_Dec2025.RData')
CS_dynamic <- merge(CS_dynamic, unique(para_CS[,-11]), by = 'para')


One_mean <- bind_rows(SK_dynamic, OK_dynamic, CS_dynamic) %>%
  group_by(Scenario, step) %>%
  summarise( N_scenar = n(),
             Mean_pI = mean(mean_p),
             sd_pI = sd(mean_p, na.rm = T),
             Mean_pS = mean(mean_pS),
             sd_pS = sd(mean_pS, na.rm = T)) %>%
  mutate(ymin = ifelse(Mean_pI - sd_pI <0, 0, Mean_pI - sd_pI),
         ymax = ifelse(Mean_pI + sd_pI > 100, 100, Mean_pI + sd_pI),
         yminS = ifelse(Mean_pS - sd_pS <0, 0, Mean_pS - sd_pS),
         ymaxS = ifelse(Mean_pS + sd_pS > 100, 100, Mean_pS + sd_pS))

Plot_3a <- ggplot(One_mean, aes(x=step, y = Mean_pI, col = Scenario))+
  geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = Scenario), alpha = 0.1, lty = 2, linewidth = 0.4) + 
  geom_line(linewidth = 1.2)+
  xlab('Time step')+
  scale_y_continuous('Proportion of cumulative infections (%)', breaks = seq(0,100, by = 25), labels = seq(0,100, by = 25), limits = c(0,100))+
  scale_color_manual('Release mode', values = c('#0072B2','#E69F00','#CC79A7'))+
  scale_fill_manual('Release mode', values = c('#0072B2','#E69F00','#CC79A7'))+
  plot_theme 

# Plot_3b <- ggplot(data2_final, aes(x = Scenario, y = Prop_Infected_Final * 100)) +
#      geom_violin(aes(fill = Scenario)) +
#      geom_boxplot(width = 0.05)+
#   scale_color_manual(values = c('black','black'))+
#   scale_fill_manual('Release mode', values = c('#0072B2','#E69F00','#CC79A7'))+
#   ylab('Proportion of cumulative infections (%)')+
#   guides(color = 'none')+
#   plot_theme +
#   theme(axis.title.x = element_blank())
# Plot_3b

ggsave('/Figures/Figure3_CumulativeI.png', width = 18, height = 10)


##############
## Figure 4 ##
##############

Plot_4a <- ggplot(Res_sim, aes(x = Scenario, y = Max_period_no_inf )) +
  geom_violin(aes(fill = Scenario), position = position_dodge(0.6))+
  geom_boxplot(width = 0.05, outlier.shape = NA, position = position_dodge(0.6))+
  scale_color_manual(values = c('black','black'))+
  scale_fill_manual('Release mode', values = c('#0072B2','#E69F00','#CC79A7'))+
  guides(color = 'none')+
  scale_y_log10(breaks = c(10,25,50,100,250,500,1000,2500,5000,7500,10000))+
  ylab('Number of time steps between consecutive infections')+
  plot_theme+
  theme(axis.title.x = element_blank())+
  ggtitle('a.')
Plot_4a

Plot_4b <- ggplot(Res_sim, aes(x = Scenario, y = Crossed_res_I)) +
  geom_violin(aes(fill = Scenario), position = position_dodge(0.6))+
  geom_boxplot(width = 0.05, outlier.shape = NA, position = position_dodge(0.6))+
  scale_color_manual(values = c('black','black'))+
  scale_fill_manual('Release mode', values = c('#0072B2','#E69F00','#CC79A7'))+
  guides(color = 'none')+
  scale_y_log10(breaks = c(1,5,10,50,100,250,500,1000,1500,3000))+
  ylab('Number of pathogen reservoir encounter before host infection')+
  plot_theme+
  theme(axis.title.x = element_blank())+
  ggtitle('b.')

Plot_4b

ggarrange(Plot_4a, Plot_4b)
ggsave('Figures/Figure4_Violin.png', width = 20, height = 12)



##############
## Figure 5 ##
##############

Plot_5 <- ggplot(Res_sim, aes(x = Scenario, y = (Number_cell_infected/1089) * 100)) +
  geom_violin(aes(fill = Gamma_plot), position = position_dodge(0.6))+
  scale_fill_manual('Infection duration', values = c('yellow','purple'))+
  geom_boxplot(width = 0.05, outlier.shape = NA, position = position_dodge(0.6), aes(col = Gamma_plot))+
  scale_color_manual(values = c('black','black'))+
  guides(color = 'none')+
  ylab('Cumulative propotion of pathogen reservoir cells (%)')+
  xlab('Release mode')+
  plot_theme

Plot_5

ggsave('Figure/Figure5_Violin_Env.png', width = 20, height = 12)


##############
## Figure 6 ##
##############
## A few functions to facilitate the plot making

#' @name Snapshot_Pinf
#' @description
#' This function allows to capture a snapshot of the pathogen spatial distribution
#' and to see the clusters created
#' @param Sce Pathogen Scenario (OK, CS or SK)
#' @param Para Parameter set indice
#' @param Para2 Indice of the parameter used for running the model
#' @param Sim Simulation number
#' @param Land Which landscape to use as background
#' @param st Step for the snapshot
#' @param Satscan_all file containing all the unique clusters
Snapshot_Pinf <- function(Sce, Para, Para2, Sim, Land, st, Satscan_all){
  
  Hab <- read.csv(paste0('/Model_output_249/HabRisk_',Para2,'_Sim',Sim,'_Para',Para,'_',Sce,'.csv'))

  Hab$xloc <- as.integer(Hab$xloc)
  Hab$yloc <- as.integer(Hab$yloc)
  Hab$P_inf <- as.numeric(Hab$P_inf)
  
  Cluster0 <- Satscan_all[Satscan_all$Scenario == Sce &
                            Satscan_all$para2 == Para2 &
                            Satscan_all$Sim == Sim,]
  
  Cluster <- Cluster0[Cluster0$First_step_appear <= st & 
                        Cluster0$Last_step_appear >= st,]
  
  Hab <- Hab[Hab$step == st,]
  
  
  Land_pinf <- merge(Hab[,c("xloc","yloc","step","P_inf")], 
                     Land[,c("Resource", "xloc", "yloc")] , by = c('xloc','yloc'), all = T)
  Land_pinf$P_inf <- ifelse(is.na(Land_pinf$P_inf),0, round(as.numeric(Land_pinf$P_inf),5))
  
  return(list(Land_pinf, Cluster, Cluster0))
}

#' @name Plot_snapshot
#' @description
#' This function allows to plot the snapshot obtained via the function Snapshot_Pinf
#' @param data data obtain using Snapshot_Pinf function
#' @param lim limits for the infectivity legend
#' @param Sce Pathogen scenario (OK, SK or CS) 
Plot_snapshot <- function(data, lim = '',Sce){
  
  d1 <- data[[1]]
  d2 <- data[[2]]
  
  if(lim[1] == ''){
    lim = c(min(d1$P_inf), max(d1$P_inf) +0.001)
  }
  
  Plot_Land <- ggplot(data=d1)+
    geom_raster(aes(x = xloc, y = yloc, fill = Resource)) +
    scale_fill_manual(name = "Resource", values = c('Low'= "navajowhite2", 'High'= "forestgreen"), na.value = NA) +  
    new_scale_fill()+
    geom_circle(data = d2, aes(x0 = as.numeric(mean.x), y0 = as.numeric(mean.y), r = as.numeric(round(Radius_mean, 0.0001)), fill = RA_mean), color = "black", alpha = 0.4) +
    geom_point(data=d1[d1$P_inf!=0,], aes(x=xloc, y=yloc, col = P_inf), size = 2) +
    scale_color_gradient('Infectivity', low = 'blue', high = 'red', limits = lim)+
    scale_fill_gradient('Infectivity', low = 'blue', high = 'red', limits = lim)+
    
    plot_theme +
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank())+
    ggtitle(Sce)
  plot(Plot_Land)
}

st = 3300
para_p <- para_CS[para_CS$para == 249,][1,] ## Low medium gamma

OK_p <- Snapshot_Pinf(Sce = 'OK', Para = para_p$para, Para2 = para_p$para2, Sim = 17, st = st, 
                      Land = get(para_p$Land), Satscan_all = Unique_cluster_p)
P_OK_lowg <- Plot_snapshot(data = OK_p, Sce= 'OK', lim = c(0, 0.05))

SK_p <- Snapshot_Pinf(Sce = 'SK', Para = para_p$para, Para2 = para_p$para2+2, Sim = 17, st = st, 
                      Land = get(para_p$Land), Satscan_all = Unique_cluster_p)
P_SK_lowg <- Plot_snapshot(data = SK_p, Sce = 'SK', lim = c(0, 0.05))

CS_p <- Snapshot_Pinf(Sce = 'CS', Para = para_p$para, Para2 = para_p$para2, Sim = 10, st = st, 
                      Land = get(para_p$Land), Satscan_all = Unique_cluster_p)
P_CS_lowg <- Plot_snapshot(data = CS_p, Sce = 'CS', lim = c(0, 0.05))


ggarrange(P_OK_lowg, P_CS_lowg, P_SK_lowg, common.legend = T, ncol = 3, nrow = 1)

ggsave('/Figures/Figure6_Spatial.png', width = 20, height = 13)


############################
######## SI figures ########
############################

###############
## Figure S1 ##
###############

pLand1 <- ggplot(Land1)+
  geom_raster(aes(x=x, y=y, fill = Resource)) +
  scale_fill_manual(name = 'Resource', values = c('Low'="navajowhite2", 'High'= "forestgreen"))+
  plot_theme+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16),
        panel.background = element_blank())+
  ggtitle('Land 1')


pLand2 <- ggplot(Land2)+
  geom_raster(aes(x=x, y=y, fill = Resource)) +
  scale_fill_manual(name = 'Resource', values = c('Low'="navajowhite2", 'High'= "forestgreen"))+
  plot_theme +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16),
        panel.background = element_blank())+
  ggtitle('Land 2')


pLand3 <- ggplot(Land3)+
  geom_raster(aes(x=x, y=y, fill = Resource)) +
  scale_fill_manual(name = 'Resource', values = c('Low'="navajowhite2", 'High'= "forestgreen"))+
  plot_theme +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16),
        panel.background = element_blank())+
  ggtitle('Land 3')


pLand4 <- ggplot(Land4)+
  geom_raster(aes(x=x, y=y, fill = Resource)) +
  scale_fill_manual(name = 'Resource', values = c('Low'="navajowhite2", 'High'= "forestgreen"))+
  plot_theme +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16),
        panel.background = element_blank())+
  ggtitle('Land 4')


ggarrange(pLand1, pLand2, pLand3, pLand4, common.legend = T)
ggsave('Figures/FigureSI1.png')

###############
## Figure S2 ##
###############

VIA_plot_SI2a <- ggarrange(SI2a, SI2b, SI2c,  
                          ncol = 1, nrow = 3, common.legend = TRUE)

ggsave(plot = VIA_plot_SI2a, 'Figures/FigureSI2a.png', width = 12, height = 18)


VIA_plot_SI2b <- ggarrange(SI2d, SI2e, SI2f,
                          ncol = 1, nrow = 3, common.legend = TRUE)

ggsave(plot = VIA_plot_SI2b, 'Figures/FigureSI2b.png', width = 12, height = 18)


VIA_plot_SI2c <- ggarrange(SI2g, SI2h, SI2i, SI2j, 
                         ncol = 1, nrow = 4, common.legend = TRUE)

ggsave(plot = VIA_plot_SI2c, 'Figures/FigureSI2c.png', width = 14, height = 22)


###############
## Figure S3 ##
###############

pdp_PropSucess <- Pdp_function_all(rf_OK = res_PropSucess_OK,
                                   rf_CS = res_PropSucess_CS,
                                   rf_SK = res_PropSucess_SK,
                                   sub.data_OK = Res_para_OK, 
                                   sub.data_CS = Res_para_CS, 
                                   sub.data_SK = Res_para_SK,
                                   ylab = 'Successful outbreaks')
ggsave(plot = pdp_PropSucess, 'Figures/FigureSI3.png' , width = 15, height = 15)

###############
## Figure S4 ##
###############
pdp_Perc_PopDie <- Pdp_function_all(rf_OK = res_Perc_PopDie_OK,
                                    rf_CS = res_Perc_PopDie_CS,
                                    rf_SK = res_Perc_PopDie_SK, 
                                    sub.data_OK = Res_para_OK, 
                                    sub.data_CS = Res_para_CS, 
                                    sub.data_SK = Res_para_SK,
                                    ylab = 'Host depletion')
ggsave(plot = pdp_Perc_PopDie, 'Figures/FigureSI4.png' , width = 15, height = 15)


###############
## Figure S5 ##
###############
pdp_Perc_DiseaseDie <- Pdp_function_all(rf_OK = res_Perc_DiseaseDie_OK,
                                        rf_CS = res_Perc_DiseaseDie_CS,
                                        rf_SK = res_Perc_DiseaseDie_SK, 
                                        sub.data_OK = Res_para_OK, 
                                        sub.data_CS = Res_para_CS, 
                                        sub.data_SK = Res_para_SK,
                                        ylab = 'Disease fade-out')
ggsave(plot = pdp_Perc_DiseaseDie, 'Figures/FigureSI5.png' , width = 15, height = 15)

###############
## Figure S6 ##
###############
pdp_PropInd <- Pdp_function_all(rf_OK = res_Perc_PropInd_OK,
                                rf_CS = res_Perc_PropInd_CS,
                                rf_SK = res_Perc_PropInd_SK, 
                                sub.data_OK = Res_sim_OK, 
                                sub.data_CS = Res_sim_CS, 
                                sub.data_SK = Res_sim_SK,
                                ylab = 'Infected individuals')
ggsave(plot = pdp_PropInd, 'Figures/FigureSI6.png' , width = 15, height = 15)

###############
## Figure S7 ##
###############
pdp_MaxPeirod <- Pdp_function_all(rf_OK = res_Perc_MaxPeirod_OK,
                                  rf_CS = res_Perc_MaxPeirod_CS,
                                  rf_SK = res_Perc_MaxPeirod_SK,
                                  sub.data_OK = Res_sim_OK, 
                                  sub.data_CS = Res_sim_CS, 
                                  sub.data_SK = Res_sim_SK,
                                  ylab = 'Period between infections')
ggsave(plot = pdp_MaxPeirod, 'Figures/FigureSI7.png', width = 15, height = 15)

###############
## Figure S3 ##
############8##
pdp_MeanInteract <- Pdp_function_all(rf_OK = res_MeanInteract_OK,
                                     rf_CS = res_MeanInteract_CS,
                                     rf_SK = res_MeanInteract_SK,
                                     sub.data_OK = Res_para_OK, 
                                     sub.data_CS = Res_para_CS, 
                                     sub.data_SK = Res_para_SK,
                                     ylab = '# encounters')
ggsave(plot = pdp_MeanInteract, 'Figures/FigureSI8.png' , width = 15, height = 15)

###############
## Figure S9 ##
###############
pdp_PropCell <- Pdp_function_all(rf_OK = res_Perc_PropCell_OK,
                                 rf_CS = res_Perc_PropCell_CS,
                                 rf_SK = res_Perc_PropCell_SK, 
                                 sub.data_OK = Res_para_OK, 
                                 sub.data_CS = Res_para_CS, 
                                 sub.data_SK = Res_para_SK,
                                 ylab = '# pathogen resevoirs')
ggsave(plot = pdp_PropCell, 'Figures/FigureSI9.png' , width = 15, height = 15)

###############
## Figure S10 ##
###############
pdp_MeanCluster <- Pdp_function_all(rf_OK = res_MeanCluster_OK,
                                    rf_CS = res_MeanCluster_CS,
                                    rf_SK = res_MeanCluster_SK,
                                    sub.data_OK = Unique_cluster_p[Unique_cluster_p$Scenario == 'OK',], 
                                    sub.data_CS = Unique_cluster_p[Unique_cluster_p$Scenario == 'CS',], 
                                    sub.data_SK = Unique_cluster_p[Unique_cluster_p$Scenario == 'SK',],
                                    ylab = 'Mean cluster peak')

ggsave(plot = pdp_MeanCluster, 'Figures/FigureSI10.png' , width = 15, height = 15)

###############
## Figure S11 ##
###############
pdp_Size_mean <- Pdp_function_all(rf_OK = res_Size_mean_OK,
                                  rf_CS = res_Size_mean_CS,
                                  rf_SK = res_Size_mean_SK,
                                  sub.data_OK = Unique_cluster_p[Unique_cluster_p$Scenario == 'OK',], 
                                  sub.data_CS = Unique_cluster_p[Unique_cluster_p$Scenario == 'CS',], 
                                  sub.data_SK = Unique_cluster_p[Unique_cluster_p$Scenario == 'SK',],
                                  ylab = '# reservoirs per cluster')
ggsave(plot = pdp_Size_mean, 'Figures/FigureSI11.png' , width = 15, height = 15)

###############
## Figure S12 ##
###############
pdp_Unique_clust <- Pdp_function_all(rf_OK = res_Unique_clust_OK,
                                     rf_CS = res_Unique_clust_CS,
                                     rf_SK = res_Unique_clust_SK,
                                     sub.data_OK = Num_unique[Num_unique$Scenario == 'OK',], 
                                     sub.data_CS = Num_unique[Num_unique$Scenario == 'CS',], 
                                     sub.data_SK = Num_unique[Num_unique$Scenario == 'SK',],
                                     ylab = '# unique clusters')
ggsave(plot = pdp_Unique_clust, 'Figures/FigureSI12.png' , width = 15, height = 15)

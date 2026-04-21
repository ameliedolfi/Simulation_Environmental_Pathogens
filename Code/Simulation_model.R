### Simulation model ###

### Final code to use on the server

#' @author Amélie Dolfi
#' @Date September 22 2023
#' @Updated : Feb 06 2026
#' @description This model aims to determine how individuals movement decisions 
#' (pertaining to environmental landscape and reservoir presence), disease dynamics 
#' (disease releasing type (Continuous-shedding, Obligate-Killer, Killer-Shedder), 
#' infectivity of individual, persistance in the environment)
#' shapes the dynamic of environmnetally transmitted pathogens
#' This model is a theoretical simulation 
#' The parameter definitions are based on real life diseases but do not reflet their dynamics per se.
#' 

Simulation_loop <- function(nsim, max_step, Shed, Death,
                            k, Land, 
                            n.initial.ind, n.initial.res, 
                            rho, gamma, mu, xi, 
                            beta, beta1, Thalf, torus, neighbor,
                            Num, Para ){
  
  #' parameter description :
  
  #' @param nsim- Number of simulation
  #' @param max_step- Number of steps per simulation
  #' @param k- Number of landscape cells : 2^k+1 / I actually created fixed landscape to run my simulations on
  #' @param Land # Landscape to run the simulation on >> Fixed landscape (Landscapes are binary 0= low-quality / 1 = high-quality)
  #' @param n.initial.ind- The number of individuals in the model at the initialization
  #' @param n.initial.res- The number of individuals that will get infected at initialization 
  #' @param rho- Initial Probability of infection risk for the reservoir (rhoD = at death; rhoS = at shedding)
  #' @param gamma- Number of steps between infection and death or Recovery (recovered individuals are removed from the simulation > Assuming they are immune)
  #' @param xi- Probability of contact with the pathogen (xi0 = in low-quality habitat and xi1 = in high-quality habitat)
  #'  @param beta- Parameter weight on the habitat (the higher weight = the higher the attraction for high-quality habitat)
  #' @param beta1- Parameter weight on the carcass effect > 0 = no effect, <0 repulsion; >0 attraction
  #' @param mu - Pathogen decay rate
  
  
  ### Packages to load
  library(DescTools)
  library(ape)  
  library(dplyr)
  library(extraDistr)
  library(tidyr)
  
  
  ################################# ######### ######################################
  ################################# FUNCTIONS ######################################
  ################################# ######### ######################################
  
  ######## Most functions comes from Lauren White code (White et al., 2018), with adaptation for some
  
  # A convenience function to use with apply when sampling values from rows within a matrix.
  #' From LAUREN WHITE PAPER
  sampfun<-function(x){
    sample(1:length(x),1,prob=x) 
  }
  
  
  #'Convenience function to calculate xy coordinates that correspond to numeric length of matrix
  #' From LAUREN WHITE PAPER
  #'assigns numbers vertically per default settings in R
  #'column-x // row-y
  #'@param dim: dimension (1D) of landscape
  longxy<-function(dim){
    ly <- rep(1:dim, dim)
    lx <- rep(1:dim, each=dim)
    lxy <- data.frame(lx,ly)
    lxy
  }
  
  
  #' Calculate the product of 2 matrix
  #' From LAUREN WHITE PAPER
  repmat<-function (x, m = 1, n = 1) {
    if (length(x) * n * m > 1e+07) {
      return(kronecker(matrix(1, n, m), x))
    }
    else {
      mx = dim(x)[1]
      nx = dim(x)[2]
      return(matrix(t(matrix(x, mx, nx * n)), mx * m, nx * 
                      n, byrow = T))
    }
  }
  
  
  ## Convert the row/col number into the element number in the habitat/ matrix
  #' From LAUREN WHITE PAPER
  Cmatrix<-function (row, col, nrows) { #y,x
    vector_loc <- row + (col - 1) * nrows
    return(vector_loc)
  }
  
  ## Convert back : from the element number, find the place in the matrix
  #' From LAUREN WHITE PAPER
  Rmatrix <- function (vector_loc, nrows){
    col <- floor(vector_loc/nrows) + 1
    row <- vector_loc - (col - 1) * nrows + 1
    return(list(row = row, col = col))
  }
  
  
  #' Function to randomly generate initial individuals across a landscape
  #' Adapted from LAUREN WHITE PAPER
  #' 
  #' @param n.initial.ind- the number of desired individuals to start with
  #' @param n.initial.res- The number of infected individuals at initialization
  #' @param dim- the length (1D) of the landscape
  #' 
  #' Individual data frame is composed of
  #' ID of individuals, xloc, ylox and vec as coordinates in the landscape
  #' Status : either S for susceptible or I for infected
  #' Res_ID for the ID of the reservoir responsible for infection
  #' Res_crossed for the number of pathogen reservoir they interacted with
  #' Step_dead for the time step at which they died
  
  
  make.inds <- function(n.initial.ind, n.initial.res, dim, st_created){
    
    ID <- 1:n.initial.ind # ID of individual
    xloc <- sample(seq(1, dim), n.initial.ind, replace = T) # Random x location to start
    yloc <- sample(seq(1, dim), n.initial.ind, replace = T) # Random y location to start
    vec <- Cmatrix(row = yloc, col = xloc, dim) # Vector value of the location
    status <- rep("S", times=n.initial.ind) # All individuals are Susceptible
    inf_step <- rep(NA, times = n.initial.ind) # Step when infection occurred
    res_ID <- rep(NA, times= n.initial.ind) # ID of the reservoir that led to the infection
    res_crossed <-  rep(0, times = n.initial.ind) # Number of reservoir the individual engaged with
    step_created <- rep(st_created, times = n.initial.ind) # Step where the individual was created (not used here but useful if you have individual dynamic)
    step_dead <- rep(NA, times = n.initial.ind) # Step where the individual died or became immune to the disease >> Step where individual is removed from the environmnet
    
    # make the dataframe
    inds <- data.frame(ID = ID, xloc=xloc, yloc=yloc, vec=vec, status=status,
                       inf_step = inf_step, res_ID = res_ID, res_crossed = res_crossed,
                       step_created = step_created, step_dead = step_dead, 
                       stringsAsFactors=FALSE) 
    
    # Randomly infect n.initial.res individuals
    
    Id_infected <- sample(inds$ID, n.initial.res, replace = FALSE)
    inds[inds$ID %in% Id_infected, 'status'] <- rep('I', length(Id_infected)) # Status of infected individuals becomes infected
    inds[inds$ID %in% Id_infected, 'res_ID'] <- rep(-1, length(Id_infected)) # -1 to say those individuals where the infected at initialization
    inds[inds$ID %in% Id_infected, 'inf_step'] <- rep(-1, length(Id_infected)) #  -1 to say those individuals where the infected at initialization
    
    return(inds)
  }
  
  
  #' Convert distances in x,y (row, col on a matrix) to distances and then weight according to a function.
  #'  From LAUREN WHITE PAPER
  #'  
  #' @param celldist matrix. Should represent a two-column matrix of change in row and column ranks, respectively.
  #' @param min.dist numeric. The distance to use for celldist==0 (i.e., for the cell of origin).
  #' @param dist.weight.fun function. A function to use for calculating the distance weight (e.g., dgamma)
  #' @param ... Other parameters to pass to dist.weight.fun.
  #' @author James Forester
  #' @export
  #' 
  get.dist.weight <-function(celldist, min.dist=0.5, dist.weight.fun=NULL,...){
    ##calculate distance from central cell to others in the neighborhood
    celldist=sqrt(rowSums(celldist^2)) #sqrt((x-0)^2+ (y-0)^2)
    celldist[celldist==0]=min.dist
    
    if(is.function(dist.weight.fun)){
      ##calculate weight, assuming a uniform bearing distribution
      cell.dist.weight=dist.weight.fun(celldist,...)/(2*pi*celldist)
    }else{
      ##default to inverse distance weight
      cell.dist.weight=1/(2*pi*celldist^2)
    }
    
    return(cell.dist.weight)
  }
  
  
  #' Calculate coordinates of neighboring cells
  #' Given a neighborhood size and starting coordinates, calculate the coordinates of all cells in the neighborhoods around each starting location.
  #' from LAUREN WHITE PAPER
  #' 
  #' @param loc matrix. Two columns representing matrix coordinates of current locations.
  #' @param mapdim vector. The number of rows and columns in the map of interest.
  #' @param rowcol.delta matrix. The relative coordinates of the neighborhood of interest (identical for all locations).
  #' @param n.offset integer. Used to define the neighborhood if needed (i.e., the number of cells on any side of the starting cell to include in the neighborhood)
  #' @param torus logical. If TRUE will wrap the map matrix into a torus to remove edge effects.
  #' @param na.val numeric. Value to assign locations that fall off of map (if torus == FALSE).
  #' @param ... Other parameters to pass through.
  #' @author James Forester
  #' @modified May 12, 2017 by Lauren White; switched order of row/columns in Cmatrix conversion
  #' @export
  #' get.neighbors(loc=loc,mapdim=mapdim,rowcol.delta=rowcol.delta,torus=TRUE, na.val=0)
  
  get.neighbors <- function(loc, mapdim, rowcol.delta=NULL, n.offset, torus, na.val=NA, ...){
    
    if(is.null(rowcol.delta)){
      rowcol.delta<- expand.grid(-n.offset:n.offset,-n.offset:n.offset)
    }
    
    neighborhood<-nrow(rowcol.delta) #e.g. an offset of 1, yields 9 neighbors
    outmat<-matrix(NA,nrow(loc),neighborhood) #rows= number of inds./col.= size of neighborhood
    
    for(i in 1:neighborhood){
      newloc=cbind(loc[,1]+rowcol.delta[i,1], loc[,2]+rowcol.delta[i,2])
      
      if(torus){
        ##Wrap map to torus
        newloc[newloc[,1]<=0,1]=mapdim[1]+newloc[newloc[,1]<=0,1]
        newloc[newloc[,2]<=0,2]=mapdim[2]+newloc[newloc[,2]<=0,2]
        newloc[newloc[,1]>mapdim[1],1] = newloc[newloc[,1]>mapdim[1],1]-mapdim[1]
        newloc[newloc[,2]>mapdim[2],2] = newloc[newloc[,2]>mapdim[2],2]-mapdim[2]
      }else{
        newloc[newloc[,1]<=0,1]=na.val
        newloc[newloc[,2]<=0,2]=na.val
        newloc[newloc[,1]>mapdim[1],1] = na.val
        newloc[newloc[,2]>mapdim[2],2] = na.val
      }
      
      outmat[,i]=Cmatrix(newloc[,1],newloc[,2],nrows=mapdim[1]) #row, col, nrows
    }
    
    return(outmat)
  }
  
  
  
  #' Movement function : Includes perception
  #' Adapted from LAUREN WHITE PAPER
  #' Individuals move one cell at a time in a Moore neighborhood
  #' They can perceive the direct cell only or up to 3 cells ahead
  #' Movement follows a Resource Selection Function based on high-habitat quality and Reservoir
  #' @param loc- A vector recording the x and y loc of each individuals
  #' @param map.mat- Matrix of environmental weight for the RSF
  
  sample.mod <- function(loc,map.mat, torus, n.offset){
    
    rowcol_step <- expand.grid(row = -n.offset:n.offset, col = -n.offset:n.offset)
    ##calculate the map locations of the neighborhood surrounding each current animal location
    mapdim <- dim(map.mat)
    
    # Get the cells within the perception range + their values
    neighbors <- get.neighbors(loc=loc, mapdim=mapdim, rowcol.delta=rowcol_step, torus = torus, na.val=NA)
    neighbor.vals = list()
    neighbor.vals <- matrix(map.mat[as.vector(neighbors)],nrow=nrow(neighbors), ncol=ncol(neighbors))
    
    ##calculate the selection function based on map values within the neighborhood
    rsf.vals <- exp(neighbor.vals)
    
    cell.dist.weight <- get.dist.weight(rowcol_step, dist.weight.fun=NULL)
    pmat <- rsf.vals * repmat(matrix(cell.dist.weight,1,length(cell.dist.weight)),nrow(rsf.vals),1)
    Final_Proba<-pmat/rowSums(pmat, na.rm=TRUE)
    Final_Proba[is.na(Final_Proba)]=0
    
    # We select the final cell according to the final probabilities
    which.cell <- apply(Final_Proba,1,sampfun) #sample(1:length(x),1,prob=x) 
    newloc <- cbind(row = loc[,1]+rowcol_step[which.cell,1], col = loc[,2]+rowcol_step[which.cell,2])
    
    ## Record new location, wrapping on torus if needed
    ##Wrap map to torus
    if(any(newloc<=0)){
      newloc[newloc[,1]<=0,1]=mapdim[1]+newloc[newloc[,1]<=0,1]
      newloc[newloc[,2]<=0,2]=mapdim[2]+newloc[newloc[,2]<=0,2]
    }
    if(any(newloc>mapdim[1])){ #assuming symmetrical landscape
      newloc[newloc[,1]>mapdim[1],1] = newloc[newloc[,1]>mapdim[1],1]-mapdim[1]
      newloc[newloc[,2]>mapdim[2],2] = newloc[newloc[,2]>mapdim[2],2]-mapdim[2]
    }
    return(newloc)
  }
  
  
  #' INFECTION
  #' From Amélie Dolfi
  #' Function to calculate infection : Look at the probability of contact and the risk of infection given contact
  #' @param inds the dataset of individuals
  #' @param reservoir the dataset of reservoirs
  #' @param P_contact Probability that the individual graze in the cell, depending on the habitat P_inter_0. P_inter_1
  #' @param hab- The habitat resource-quality layer
  
  Infection <- function(inds, reservoir, P_contact, hab){
    
    ## For each individuals on a cell with a reservoir
    if (any(inds$vec %in% reservoir$vec)){
      
      Risk_cells <- which(reservoir$vec %in% inds$vec)  ## Which cells infection might happen = cells in which both a reservoir and individual are present
      Ind_risk <- data.frame(ID = inds[inds$vec %in% reservoir[Risk_cells,'vec'],'ID'], ## Get the ID and vector location of the individuals at risk
                             vec = inds[inds$vec %in% reservoir[Risk_cells,'vec'],'vec'])
      Ind_risk <- merge(Ind_risk, subset(reservoir, select = c(vec, type)), by = 'vec', all.x = T, all.y = F) # Get the type of reservoir (ie: Death or Shed)
      
      Hab_type <- hab[inds[inds$ID %in% Ind_risk$ID, 'vec']] # Extract the cell habitat type (low or high quality)
      Proba <- ifelse(Hab_type == 0, P_contact[1], P_contact[2]) ## Contact rate differ per habitat quality
      Interact <- ifelse(runif(length(Proba)) > Proba, 'N','Y') # Draw the probability of interaction depending on the Probability of contact
      
      if (length(Interact == 'Y')!=0){ # if contact occurs: Draw the probability of infection
        Ind_Interact <- Ind_risk[which(Interact=='Y'),] # Which individuals interact with the reservoir (ie. individuals at risk of infection)
        
        P_inf <- reservoir[reservoir$vec %in% inds[inds$vec %in% Ind_Interact$vec, 'vec'], 'P_inf']  # Probability of infection for each individual (obtain from the reservoir value)
        Ind_I <- which(runif(length(P_inf)) < P_inf ) ## Get which individuals are infected
        
        # Modify the attributes of infected agents
        if(length(Ind_I) != 0 ){
          inds[inds$ID %in% Ind_Interact[Ind_I, 'ID'], 'status'] <- 'I' ## Change their status to 'I'
          inds[inds$ID %in% Ind_Interact[Ind_I, 'ID'], 'inf_step'] <- step ## Record infection step
          
          
          for (j in Ind_I){
            inds[inds$ID == Ind_Interact[j,'ID'], 'res_ID']  <- reservoir[reservoir$vec == Ind_Interact[j,'vec'],'ID'] # add which reservoir led to the infection
            reservoir[reservoir$vec == Ind_Interact[j,'vec'], 'ind_inf'] <- reservoir[reservoir$vec == Ind_Interact[j,'vec'], 'ind_inf']  + 1 # Add 1 infection for each reservoir that led to infection
          }
        }
        
        inds[inds$ID %in% Ind_Interact[,'ID'], 'res_crossed']  <-  inds[inds$ID %in% Ind_Interact[,'ID'], 'res_crossed']+ 1 # Increase the number of reservoir an individual interacted with 
        
        
        ## Add the number of individuals interacting with the reservoir
        for (k in (Ind_Interact$ID)){
          reservoir[reservoir$vec %in% inds[inds$ID == k,'vec'], 'ind_crossed'] <- reservoir[reservoir$vec %in% inds[inds$ID == k,'vec'], 'ind_crossed'] +1 # Add 1 interaction to each reservoir where individuals were interacting with 
          
        }
      }
    }
    return(list(inds, reservoir))
  }
  
  
  #' Release of pathogen into the environment
  #' From Amélie Dolfi
  #' @param release_type the type of release: Shed or Death
  #' @param Release_inds The infected individuals dataset
  #' @param res the reservoir dataset
  #' @param P_infection the intial pathogen load
  
  Release_pathogen <- function(release_type, Release_inds, res, P_infection){
    
    ## Obtain the value of the infection probability depending on the release type
    if(release_type == 'Shed'){ Release_inds$P_inf <- P_infection[2]
    }else{ Release_inds$P_inf <- P_infection[1]}
    
    ## If multiple die in same cell, aggregate the probability     
    if(length(unique(Release_inds$vec))!= nrow(Release_inds)){ 
      Release_inds <- aggregate(P_inf ~  vec + xloc + yloc, Release_inds, sum)
    }
    
    ## If a reservoir site already exist, then aggregate the P_inf          
    if(any(Release_inds$vec %in% res$vec)){ 
      ID_in_res <- Release_inds[which(Release_inds$vec %in% res$vec), ]
      res[res$vec %in% ID_in_res$vec, 'P_inf'] <- res[res$vec %in% ID_in_res$vec, 'P_inf'] + ID_in_res$P_inf
      Release_inds <- Release_inds[-which(Release_inds$vec %in% res$vec),]
    }
    
    # Update the reservoir data
    
    if(nrow(Release_inds) !=0){
      
      res <- rbind(res, cbind(ID = seq(ID_res+1, ID_res+nrow(Release_inds)), 
                              xloc =Release_inds[,'xloc'], 
                              yloc = Release_inds[,'yloc'], 
                              vec= Release_inds['vec'], 
                              type = rep(release_type, nrow(Release_inds)),
                              status = rep('I', nrow(Release_inds)), 
                              P_inf = Release_inds[, 'P_inf'],
                              ind_crossed = rep(0, nrow(Release_inds)),
                              ind_inf = rep(0, nrow(Release_inds)),
                              step_created = rep(step, nrow(Release_inds)),
                              step_dead = rep(NA, nrow(Release_inds))))
    }
    
    # Ensure the probability never goes above 1
    res$P_inf <- ifelse(res$P_inf >= 1 , 1, res$P_inf)
    ID_res <- ID_res + nrow(Release_inds) ## Get the ID number to start at
    
    return(list(res, ID_res))
  }
  
  
  
  
  ################################# ##### ######################################
  ################################# MODEL ######################################
  ################################# ##### ######################################
  
  # File name to save the outputs
  
  # filename <- paste('D-',Death,'-S-',Shed, '-k-',k,'_',Land,
  #                   '-NS', n.initial.ind,  '-NI',n.initial.res, 
  #                   '-xi', xi[[1]],'-', xi[[2]],
  #                   '-rho', rho[1],'-',rho[2], 
  #                   '-B', beta,'-B1_', beta1,
  #                   '-gamma',gamma, '-mu', mu, '-neigh',neighbor,'.csv',sep= '')
  
  
  # Calculating the minimum value of Probability of infection of reservoirs before it is removed from the environment 
  Pinf_minimum <- rho[2] * (1 - mu) ** ( (-log(2)) / log((1 - mu))*Thalf)
  
  # Summary output
  summary_tot <- data.frame(Scenario = character(),
                            Index = numeric(),
                            Para = numeric(),
                            sim = numeric(),
                            step= numeric(),
                            xi0 = numeric(),
                            xi1 = numeric(),
                            mu = numeric(),
                            Land = character(),
                            beta=numeric(),
                            beta1=numeric(),
                            Gamma = numeric(),
                            RhoD = numeric(),
                            RhoS = numeric(),
                            Per = numeric(),
                            Num_S = numeric(),
                            Num_I=numeric(),
                            Num_I_dead = numeric(),
                            Num_res = numeric(),
                            Num_res_dead = numeric())
  # Individual output
  Inds_tot <- data.frame(ID = numeric(),
                        xloc = numeric(),
                        yloc = numeric(),
                        vec = numeric(),
                        status = character(),
                        inf_step = numeric(),
                        res_ID = numeric(),
                        res_crossed = numeric(),
                        step_created = numeric(),
                        step_dead = numeric(),
                        sim = numeric())
  # Dead individuals output
  Ind_Dead_tot <- Inds_tot  
  #Reservoir output
  Res_tot <- data.frame(ID = numeric(),
                        xloc = numeric(),
                        yloc = numeric(),
                        vec = numeric(),
                        type = character(),
                        status = character(),
                        P_inf = numeric(),
                        ind_crossed = numeric(),
                        ind_inf = numeric(),
                        step_dead = numeric(),
                        step_created = numeric(),
                        sim = numeric())
  #Dead reservoir output
  Res_Dead_tot <- Res_tot
  
  # I run the model to obtain nsim simulation that are successfull
  Failed_simulation <- 0
  Successfull_Simulation <- 0
  sim <- 0

  while(Successfull_Simulation != nsim){
    print(Successfull_Simulation)
    sim <- Successfull_Simulation +1
    print(paste0('Failed ', Failed_simulation))
    print(paste0('Sucess ', Successfull_Simulation))
    
    ### Create the habitat
    hab1 <- as.matrix(read.csv(paste0(Land, '.csv'))) ## 
    #hab1 <- as.matrix(read.csv(paste0('C:/Users/Amelie/Documents/PhD/Chap2_Simulation/BACKON_SEPT2024/DATA/',Land, '.csv')))
    lxy <- longxy(nrow(hab1))
    vec <- Cmatrix(row = lxy$ly, col = lxy$lx, nrow(hab1))
    dim = nrow(hab1)
    
    
    hab_risk <- data.frame(vec = numeric(),
                           xloc = numeric(),
                           yloc = numeric(),
                           step = numeric(),
                           P_inf = numeric())
    
    ## Create the initiate individuals
    
    inds <- make.inds(n.initial.ind =n.initial.ind, n.initial.res = n.initial.res, dim=nrow(hab1), st_created = 0)
    #ggplot()+ geom_spatraster(data = rast(hab1)$lyr.1) + geom_point(data=inds, aes(x=xloc-0.5, y=yloc-0.5), col = 'red') 
    
    ## Create the reservoir, dead individuals and dead reservoir output df
    reservoir <- data.frame(ID = numeric(), xloc=numeric(), yloc=numeric(), vec=numeric(), 
                            type = character(), status = character(), P_inf = numeric(),
                            ind_crossed = numeric(), ind_inf= numeric(),
                            step_dead = numeric(), step_created = numeric()) 
    
    inds_dead <- inds[0,]
    reservoir_dead <- reservoir[0,]
    


    ##### SAVE INITIAL DATA
    ID_res <- 0
    
    summary_step <- data.frame(Scenario = paste0('Death_', Death,'_Shed_',Shed),
                               Index = Num,
                               Para = Para,
                               sim = sim,
                               step= 1, 
                               xi0 = xi[1],
                               xi1 = xi[2],
                               mu = mu,
                               Land=Land,  
                               beta=beta, 
                               beta1=beta1, 
                               Gamma = gamma,
                               RhoD = rho[1], 
                               RhoS = rho[2], 
                               Per = neighbor,
                               Num_S = n.initial.ind - n.initial.res,
                               Num_I= n.initial.res , 
                               Num_I_dead = 0,
                               Num_res = 0,
                               Num_res_dead = 0)
    
    
    
    
    for (step in (2: max_step)){
      
      #### SUBMODEL 1 : MOVEMENT
      ##Calculate  RSF map
      ## beta is for the resource impact
      ## beta1 is for carcass impact
      ## Environmental map risk
      
      # I save the habitat risk map (ie the P_inf of every cells) every 150 steps
      if(step %% 150 == 0 || step == max_step){
        hab_risk <- rbind(hab_risk,
                          cbind(vec = reservoir[,'vec'],
                                xloc = reservoir[,'xloc'],
                                yloc = reservoir[,'yloc'],
                                step = rep(step, nrow(reservoir)),
                                P_inf = reservoir[,'P_inf']))
        
      }
      
      ### Environmental reservoir map
      habcarc <- matrix(0,2^k+1, 2^k+1)
      habcarc[reservoir[,'vec']] <- 1
      
      # Get the linear predictor for the RSF 
      linear.pred <- (beta*hab1) + (beta1*habcarc)
      
      ##Calculate the new locations based on the neighborhood around each of previous locations
      newloc <- as.data.frame(sample.mod(loc=as.data.frame(Rmatrix(inds$vec-1,nrow(hab1))), n.offset = neighbor, map.mat=linear.pred, torus = torus))
      inds$xloc <- newloc$col # col
      inds$yloc <- newloc$row # row
      inds$vec <-Cmatrix(row = newloc$row, col= newloc$col, nrow(hab1)) #row, col, nrow
      
      
      ##### SUBMODEL 2-3 : INTERACTION WITH RESERVOIR AND INFECTION
      
      ### Infection risk : If hosts are in the same cell as a reservoir, 
      ###   infection may occur depending on the probability of interaction and
      ###   the intensity of the interaction
      
      
      Inf1 <- Infection(inds[inds$status == 'S',], reservoir,
                        P_contact = xi,  hab=hab1)
      
      inds <- rbind(Inf1[[1]], inds[inds$status == 'I',])
      reservoir <- Inf1[[2]]
      
      
      #### SUBMODEL 4 : PATHOGEN RELEASE: EITHER AT DEATH OR SHED and remove the individuals not infected anymore
      
      Inf_inds <- inds[inds$status=='I',]
      # Draw which infected individuals will die/ be removed from the simulation
      Max_inf_ind <- Inf_inds[which(runif(nrow(Inf_inds)) >= gamma), ]
      
      ## If individuals shed they do so at every steps (even the one they will be removed/died)
      if (Shed == TRUE){
        
        Shed_inds <- subset(inds, status == 'I')
        
        if (nrow(Shed_inds) != 0){
          Release_shed <- Release_pathogen(release_type = 'Shed', 
                                           Release_inds = Shed_inds, res = reservoir, P_infection = rho)
          
          reservoir <- Release_shed[[1]]
          ID_res <- Release_shed[[2]]
          
        }
      }
      
      # If any individual die or is removed 
      if (nrow(Max_inf_ind) != 0){
        if (Death == TRUE ){
          
          Release_death <- Release_pathogen(release_type = 'Carcass',  Release_inds = Max_inf_ind,
                                            res = reservoir, P_infection = rho)
          
          reservoir <- Release_death[[1]]
          ID_res <- Release_death[[2]]
          
        }
        
        ## Remove the deads >> All individuals die of disease, but not all of them create an infectious carcass       
        inds[inds$ID %in% Max_inf_ind$ID, 'step_dead'] <- step
        inds_dead <- rbind(inds_dead, inds[inds$ID %in% Max_inf_ind$ID,])
        inds <- inds[!inds$ID %in% Max_inf_ind$ID,]
        
      }
      
      
      ### SUBMODEL 5 : PATHOGEN CONCENTRATION DECREASE : 
      ## When Pathogen concentration is below P_inf minimum, take out for the environment
      
      
      reservoir[reservoir$step_created != step, 'P_inf'] <- reservoir[reservoir$step_created != step, 'P_inf'] * (1-mu)   
      reservoir[reservoir$step_created != step, 'status'] <- ifelse(reservoir[reservoir$step_created != step, 'P_inf'] < Pinf_minimum, 'D', reservoir[reservoir$step_created != step, 'status'])
      reservoir[reservoir$step_created != step, 'step_dead'] <- ifelse(reservoir[reservoir$step_created != step, 'status'] == 'D', step, NA)
      
      reservoir_dead <- rbind(reservoir_dead, reservoir[reservoir$status == 'D',])
      reservoir <- reservoir[reservoir$status != 'D',]

      
      ### Save the outputs
      
      summary_step <-rbind(summary_step,
                           data.frame(Scenario = paste('Death_', Death,'_Shed_',Shed, sep ='_'),
                                      Index = Num,
                                      Para = Para,
                                      sim = sim,
                                      step= step, 
                                      xi0 = xi[1],
                                      xi1 = xi[2],
                                      mu = mu,
                                      Land=Land,  
                                      beta=beta, 
                                      beta1=beta1, 
                                      Gamma = gamma,
                                      RhoD = rho[1], 
                                      RhoS = rho[2], 
                                      Per = neighbor,
                                      Num_S =  nrow(inds[inds$status=='S',]),
                                      Num_I= nrow(inds[inds$status=='I',]) , 
                                      Num_I_dead =  nrow(inds_dead),
                                      Num_res = nrow(reservoir),
                                      Num_res_dead =  nrow(reservoir_dead)))
      
  
      if(nrow(inds) == 0){ # If no individuals left. stop running
        print(paste0('no indiviuals left after ',step,' steps'))
        break
      }
      
      if(nrow(reservoir) == 0 & sum(inds$status=='I') == 0){ #If number of individuals infected and reservoir is 0, stop running
        print(paste0('no reservoir no infectious left after ',step,' steps'))
        break
      }
    }
    
    if(summary_step[nrow(summary_step), 'Num_I_dead'] < 5){
      Failed_simulation <- Failed_simulation + 1
    }else{
      Successfull_Simulation <- Successfull_Simulation  +1
      #Add the summary to the overall total
      summary_tot <- rbind(summary_tot, summary_step)
      
      All_inds <- rbind(inds, inds_dead)
      All_inds$sim <- sim
      
      Inds_tot <- rbind(Inds_tot, All_inds)
      
      
      
      #Save the infectious habitat
      write.csv(hab_risk, file=paste0('HabRisk_',Num,'_Sim',sim,'_Para',Para,'.csv'), row.names = F)
      
    }
    
    if(Failed_simulation > nsim * 100){
      break
    }
    
    print('DONE')
  }

  summary_tot <- rbind(summary_tot,
                       data.frame(Scenario = 'Finished',
                                  Index = ifelse(Successfull_Simulation != nsim, F, T),
                                  Para = 'Simulation_failed',
                                  sim = Failed_simulation,
                                  step= 'Simulation_Sucess',
                                  xi0 = Successfull_Simulation, 
                                  xi1 = 'Total_simulations',
                                  mu = Successfull_Simulation + Failed_simulation,
                                  Land = NA,
                                  beta= NA,
                                  beta1= NA,
                                  Gamma = NA,
                                  RhoD = NA,
                                  RhoS = NA,
                                  Per = NA,
                                  Num_S = NA,
                                  Num_I= NA,
                                  Num_I_dead = NA,
                                  Num_res = NA,
                                  Num_res_dead = NA))
  
  write.csv(summary_tot, file=paste0('Sum_',Num,'_Para',Para,'.csv'), row.names = F)
  write.csv(Inds_tot, file=paste0('Inds_',Num,'_Para',Para,'.csv'), row.names = F)
  
}



## OK run
parameter_OK <- read.csv("Data_input/Parameters_OK.csv", sep=',')

for (i in nrow(parameter_OK)){
  
  para <- parameter_OK[index,]
  set.seed(i)
  Simulation_loop(nsim = 20, max_step = 10000, Shed = para$Shed, Death = para$Death,
                  k = 5 , Land = para$Land ,
                  n.initial.ind = para$N, n.initial.res= 1,
                  rho = c(para$rhoD, para$rhoS), gamma = para$gamma, mu = para$mu,
                  xi = c(0.25, 0.5),
                  beta = para$B, beta1 = para$B1, Thalf = 2, torus = FALSE, neighbor = para$per,
                  Num = index, Para = para$para)
}


## CS run

parameter_CS <- read.csv("Data_input/Parameters_CS.csv", sep=',')


for (i in nrow(parameter_CS)){
  
  para <- parameter_CS[index,]
  set.seed(i)
  Simulation_loop(nsim = 20, max_step = 10000, Shed = para$Shed, Death = para$Death,
                  k = 5 , Land = para$Land ,
                  n.initial.ind = para$N, n.initial.res= 1,
                  rho = c(para$rhoD, para$rhoS), gamma = para$gamma, mu = para$mu,
                  xi = c(0.25, 0.5),
                  beta = para$B, beta1 = para$B1, Thalf = 2, torus = FALSE, neighbor = para$per,
                  Num = index, Para = para$para)
}


## SK run

parameter_SK <- read.csv("Data_input/Parameters_SK.csv", sep=',')


for (i in nrow(parameter_SK)){
  
  para <- parameter_SK[index,]
  set.seed(i)
  Simulation_loop(nsim = 20, max_step = 10000, Shed = para$Shed, Death = para$Death,
                  k = 5 , Land = para$Land ,
                  n.initial.ind = para$N, n.initial.res= 1,
                  rho = c(para$rhoD, para$rhoS), gamma = para$gamma, mu = para$mu,
                  xi = c(0.25, 0.5),
                  beta = para$B, beta1 = para$B1, Thalf = 2, torus = FALSE, neighbor = para$per,
                  Num = index, Para = para$para)
}
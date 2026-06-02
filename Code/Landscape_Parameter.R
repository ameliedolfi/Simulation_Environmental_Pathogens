#' @title Parameter tests and landscape test
#' @author Amelie Dolfi
#' @date  June 1st 2026


#rm(list=ls())

library(tidyverse)
library(terra)
library(lhs)

# Parameter range determined for the MS:

set.seed(123)

# Define number of samples and parameters
n_samples <- 500  # Number of parameter sets
n_parameters <- 8 # Number of varying parameters

# Generate Optimized Latin Hypercube Sample (better spread)
lhs_sample <- maximinLHS(n_samples, n_parameters)  

# Adjust LHS to better cover extreme values
lhs_sample <- lhs_sample  # Forces more sampling near boundaries

# Define parameter ranges
N <- c(300, 800)        # Number of individuals
rhoD <- c(0.0025, 0.025) # Infective cells released
mu <- c(0.0002, 0.002)  # Pathogen decay rate
gamma <- c(0.7, 0.95)   # Infectivity time
B <- c(0, 6)            # Attraction toward good habitat
B1 <- c(-1,1)          # Attraction/Avoidance of infectious cells
Land <- c(1, 4)         # Landscape type
per <- c(1, 3)          # Perception (movement range)

# Convert LHS continuous values to parameter values
Test_LHS_random <- data.frame(
  N = round(((N[2] - N[1]) * lhs_sample[,1]) + N[1], 0),
  B = round(((B[2] - B[1]) * lhs_sample[,2]) + B[1], 1),
  B1 = round(((B1[2] - B1[1]) * lhs_sample[,3]) + B1[1], 1),
  gamma = round(((gamma[2] - gamma[1]) * lhs_sample[,4]) + gamma[1], 2),
  mu = round(((mu[2] - mu[1]) * lhs_sample[,5]) + mu[1], 5),
  rhoD = round(((rhoD[2] - rhoD[1]) * lhs_sample[,6]) + rhoD[1], 5),
  per = cut(lhs_sample[,7], breaks=3, labels=c(1, 2, 3)),  # Stratify `per` into 3 equal groups
  Land = cut(lhs_sample[,8], breaks=4, labels=c(1, 2, 3, 4))  # Stratify `Land` into 4 equal groups
)

# Convert factors to integers
Test_LHS_random$per <- as.integer(as.character(Test_LHS_random$per))
Test_LHS_random$Land <- as.integer(as.character(Test_LHS_random$Land))


### Let's check the representation of each parameters
ggplot(data = Test_LHS_random, aes(x=N))+
  geom_density()
ggplot(data = Test_LHS_random, aes(x=rhoD))+
  geom_density()
ggplot(data = Test_LHS_random, aes(x=mu))+
  geom_density()
ggplot(data = Test_LHS_random, aes(x=gamma))+
  geom_density()
ggplot(data = Test_LHS_random, aes(x=per))+
  geom_density()
ggplot(data = Test_LHS_random, aes(x=Land))+
  geom_density()
ggplot(data = Test_LHS_random, aes(x=B))+
  geom_density()
ggplot(data = Test_LHS_random, aes(x=B1))+
  geom_density()

table(Test_LHS_random$Land)
table(Test_LHS_random$per)
table(Test_LHS_random$B)
table(Test_LHS_random$B1)
table(Test_LHS_random$N)

# Change the name back to Land
Test_LHS_random$Land <- paste0('Land',Test_LHS_random$Land)
# We assume that the number of pathogen released at sheding is 10 times lower than realeased at death
Test_LHS_random$rhoS <- Test_LHS_random$rhoD /10 

Test_LHS_random$para <- 1:500 # Get the number of parameters
# Replicating each row 5 times as I will run 20 simulations per rows on the server (makes it easier and faster)
Test_LHS_random <- Test_LHS_random[rep(1:nrow(Test_LHS_random), each = 5),]
Test_LHS_random$para2 <- 1:nrow(Test_LHS_random) # Get the number of lines, to make it easier to navigate the results from the server

parameters_OK <- Test_LHS_random
parameters_OK$Death = TRUE
parameters_OK$Shed = FALSE

write.csv(parameters_OK, 'Data_input/Parameters_OK.csv', row.names = FALSE)

parameters_CS <- Test_LHS_random
parameters_CS$Death = FALSE
parameters_CS$Shed = TRUE
write.csv(parameters_CS, 'Data_input/Parameters_CS.csv', row.names = FALSE)


parameters_SK <- Test_LHS_random
parameters_SK$Death = TRUE
parameters_SK$Shed = TRUE
write.csv(parameters_SK, 'Data_input/Parameters_SK.csv', row.names = FALSE)



### Make the landscape:

#' Create neutral landscape maps
#' From LAUREN WHITE PAPER
#' Use standard methods to generate fractal maps. Binary and continuous surfaces may be produced.
#' 
#' @param k integer. The extent of the map (2^k+1)^2 pixels
#' @param h numeric. Level of aggregation in the map.
#' @param p numeric (0,1). The proportion of map in habitat=1
#' @param binary logical. If TRUE, a 0/1 categorical landscape is produced.
#' @author Shannon Pittman, James Forester, modified by Lauren White
#' @export
#' @example examples/neutral.landscape_example.R
#' 
fracland_mod <- function(k, h, p, binary = TRUE) {
  
  ## Function for creating neutral landscapes Shannon Pittman University of Minnesota May, 2013 k = the extent of the map (2^k+1)^2 pixels h =
  ## how clumped the map should be (ranging from ?? to ??) -- weird behavior at higher values p = proportion of map in habitat 1 binary =
  ## plotflag == if TRUE will plot a filled contour version of the matrix
  
  ## function call: testmap=land(6,1,.5,FALSE,TRUE)
  A <- 2^k + 1  # Scalar-determines length of landscape matrix
  
  if(p > 0.95){
    return(matrix(1, A, A) )
  }
  else{
    
    #Right now, as written (1-p) represents the amount of habitat listed as "1"
    B <- matrix(0, A, A)  # Creates landscape matrix
    
    B[1, 1] <- 0
    B[1, A] <- 0
    B[A, 1] <- 0
    B[A, A] <- 0
    
    iter <- 1
    for (iter in 1:k) {
      scalef <- (0.5 + (1 - h)/2)^(iter)
      
      d <- 2^(k - iter)
      
      # ALL SQUARE STEPS#
      for (i in seq(d + 1, A - d, 2 * d)) {
        for (j in seq(d + 1, A - d, 2 * d)) {
          B[i, j] <- mean(c(B[i - d, j - d], B[i - d, j + d], B[i + d, j - d], B[i + d, j + d])) + scalef * rnorm(n = 1)
        }
      }
      
      
      # OUTSIDE DIAMOND STEP#
      for (j in seq(d + 1, A - d, 2 * d)) {
        B[1, j] <- mean(c(B[1, j - d], B[1, j + d], B[1 + d, j])) + scalef * rnorm(n = 1)
        B[A, j] <- mean(c(B[A, j - d], B[A, j + d], B[A - d, j])) + scalef * rnorm(n = 1)
      }
      
      for (i in seq(d + 1, A - d, 2 * d)) {
        B[i, 1] <- mean(c(B[i - d, 1], B[i + d, 1], B[i, 1 + d])) + scalef * rnorm(n = 1)
        B[i, A] <- mean(c(B[i - d, A], B[i + d, A], B[i, A - d])) + scalef * rnorm(n = 1)
      }
      
      # INSIDE DIAMOND STEP#
      if (2 * d + 1 <= A - 2 * d) {
        for (i in seq(d + 1, A - d, 2 * d)) {
          for (j in seq(2 * d + 1, A - 2 * d, 2 * d)) {
            B[i, j] <- mean(c(B[i - d, j], B[i + d, j], B[i, j - d], B[i, j + d])) + scalef * rnorm(n = 1)
          }
        }
        
        
        for (i in seq(2 * d + 1, A - 2 * d, 2 * d)) {
          for (j in seq(d + 1, A - d, 2 * d)) {
            B[i, j] <- mean(c(B[i - d, j], B[i + d, j], B[i, j - d], B[i, j + d])) + scalef * rnorm(n = 1)
          }
        }
        
      }
      
      iter <- iter + 1
    }
    
    if (binary == T) {
      R <- sort(B)
      PosR <- (1 - p) * length(R)  #larger values become habitat, designated as 1
      pval <- R[PosR]
      T1 <- which(B > pval)
      T2 <- which(B <= pval)
      B[T1] <- 1  #habitat is 1
      B[T2] <- 0
    } 
    return(B)
  }
}

#set.seed(12)
set.seed(3)
Land1 <- fracland_mod(k = 5, h = 0.2, p = 0.3)
Land2 <- fracland_mod(k = 5, h = 0.8, p = 0.3)
Land3 <- fracland_mod(k = 5, h = 0.2, p = 0.7)
set.seed(123)
Land4 <- fracland_mod(k = 5, h = 0.8, p = 0.7)

par(mfrow=c(2,2))
plot(rast(Land1))
plot(rast(Land2))
plot(rast(Land3))
plot(rast(Land4))


write.csv(Land1, 'Data_input/Land1.csv', row.names = FALSE)
write.csv(Land2, 'Data_input/Land2.csv', row.names = FALSE)
write.csv(Land3, 'Data_input/Land3.csv', row.names = FALSE)
write.csv(Land4, 'Data_input/Land4.csv', row.names = FALSE)



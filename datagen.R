# Load packages
if(!require(data.table)) {install.packages("data.table"); require(data.table)}
if(!require(plyr)) {install.packages("plyr"); require(plyr)}
if(!require(tidyverse)) {install.packages("tidyverse"); require(tidyverse)}

############################## Point treatment #######################

####### Data generating function ####### 

datagen.point <- function(i, Z, alpha0, alpha1, alpha2, alpha3, alpha4, 
                          theta0, theta1, theta2, theta3, theta4,
                          beta1_0, beta1_1, sigma,
                          beta2_0, beta2_1, beta2_2){
  
  id <- as.numeric(i) # Define individual id number
  
  # Unmeasured confounding
  U <- runif(1, 0, 1) # Uniformly distributed between 0 and 1
  
  # Continuous covariate B1
  B1 <- rnorm(1, mean = beta1_0 + beta1_1*U, sd = sigma)
  
  # Binary covariate B2
  B2 <- rbinom(1, 1, plogis(beta2_0 + beta2_1*U + beta2_2*B1))
  
  # Binary treatment exposure indicator A
  A <- rbinom(1, 1, plogis(alpha0 + alpha1*Z + alpha2*B1 + alpha3*B2 + alpha4*U))
  
  # Binary outcome indicator Y
  Y <- rbinom(1, 1, plogis(theta0 + theta1*A + theta2*B1 + theta3*B2 + theta4*U))
  
  # Generate vectors to build data.frame for data at initial time point
  id_ <- c(id)
  U_  <- c(U)
  B1_ <- c(B1)
  B2_ <- c(B2)
  A_ = c(A)
  Y_ = c(Y)
  Z_ <- c(Z)
  
  # Consolidate data in a single data frame
  temp_data <- data.frame(Z = Z_, id = id_, U = U_, B1 = B1_, B2 = B2_, 
                          A = A_, Y = Y_)
}

####### Data generate  ####### 

iters <- 1

n <- 1000 # Number of subjects in each arm

set.seed(321)
for (iter in 1:iters){
  
  # Generate Z=1 dataset for n individuals 
  Z=1 
  
  # Parameters for data generating functions
  alpha0 = -0.23; alpha1 = 0.6; alpha2 = 0.4; alpha3 = 0.35; alpha4 = 0;
  theta0 = -3.5; theta1 = 0.7; theta2 = 0; theta3 = 0; theta4 = 1
  beta1_0 = 0; beta1_1 = 6; beta2_0 = -5; beta2_1 = 3; beta2_2 = 1.25; sigma = 2
  
  # Generate data with datagen.point function
  df <- lapply(as.list(1:n), FUN=function(ind){
    datagen.point(ind, Z = Z, beta1_0 = beta1_0, beta1_1 = beta1_1, sigma = sigma, 
                  beta2_0 = beta2_0, beta2_1 = beta2_1, beta2_2 = beta2_2,
                  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
                  alpha3 = alpha3, alpha4 = alpha4, theta0 = theta0, theta1 = theta1, 
                  theta2 = theta2, theta3 = theta3, theta4 = theta4)
  })
  
  # Store final data frame for Z=1 arm
  simdat1 <- rbindlist(df)
  
  # Generate Z=0 dataset for n individuals 
  Z=0 
  
  # Updated parameters for data generating functions for Z = 0 arm
  alpha0 = -3.14
  
  # Generate data with datagen.point function
  df <- lapply(as.list(1:n), FUN=function(ind){
    datagen.point(ind + n, Z = Z, beta1_0 = beta1_0, beta1_1 = beta1_1, sigma = sigma, 
                  beta2_0 = beta2_0, beta2_1 = beta2_1, beta2_2 = beta2_2,
                  alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
                  alpha3 = alpha3, alpha4 = alpha4, theta0 = theta0, theta1 = theta1, 
                  theta2 = theta2, theta3 = theta3, theta4 = theta4)
  })
  
  # Store final dataframe for Z=0 arm
  simdat0 <- rbindlist(df)
  
  # Concatenate both data frames
  dat <- rbind(simdat1,simdat0)
}

# Save the data in a .csv file
write.csv(dat, file = "DataPoint.csv", row.names = F)
###########################################################################



rm(list = ls())
############################## Sustained treatment #######################
# For details, see Young et al. [2019]

####### Data generating function ####### 

datagen.sustained <- function(i, K, Z, m, sigma,
                              alpha0, alpha1, alpha2, alpha3, alpha4,
                              beta1_0, beta1_1, beta1_2, beta1_3, beta1_4,
                              beta1_5, beta1_6, beta1_7, beta2_0, beta2_1, 
                              beta2_2, beta2_3, beta2_4, beta2_5, beta2_6,
                              beta2_7, theta0, theta1, theta2, theta3, 
                              theta4, theta5, theta6){
  id <- as.numeric(i) # Define individual id number
  
  # Define baseline time
  time <- 0
  
  # Generate baseline common cause of time-varying covariates
  U <- runif(1, 0, 1) # Uniformly distributed between 0 and 1
  
  # Generate baseline data
  
  # Continuous covariate L1
  L1 <- rnorm(1, mean = beta1_0 + beta1_1*U, sd = sigma)
  L1cavg <- cumsum(L1)[1]/1 # Calculate cumavg(L1)
  
  # Binary covariate L2
  L2 <- rbinom(1, 1, plogis(beta2_0 + beta2_1*U + beta2_2*L1))
  
  # Define lagged value of L2
  L2lag <- 0 
  
  # Binary treatment exposure indicator A
  A <- rbinom(1, 1, plogis(alpha0 + alpha1*L1cavg + alpha2*L2lag + alpha3*time + 
                             alpha4*U))
  
  # In case of missing data (data only measured at months n*m, for n={1,2,3,...}),
  # create variables that store previously measured values of L1cumavg and  A
  sumuse <- L1 # Used in calculation of previously measured value of L1cumavg
  avguse <- L1
  Ause <- A     
  
  # Binary outcome indicator Y
  Y <- rbinom(1, 1, plogis(theta0 + theta1*U + theta2*A + theta3*L1 + theta4*L2))
  
  # Coerce NA to num for future data
  enay <- 0
  enay <- NA
  
  # Generate vectors to build data.frame for data at initial time point
  id_ <- c(id)
  time_ <- c(time)
  U_  <- c(U)
  L1_ <- c(L1)
  L1cavg_ <- c(L1cavg)
  L2_ = c(L2)
  L2lag_ = c(L2lag)
  A_ = c(A)
  Y_ = c(Y)
  Z_ <- c(Z)
  
  # Lagged value of A at initial time point corresponds to assigned study arm
  Alag1_ <- c(Z) 
  
  # Variables for previously measured values of L1cumavg, A, and lagged A in case 
  # of infrequently measured data
  sumuse_ <- c(sumuse)
  avguse_ <- c(avguse)
  Ause_ <- c(Ause)
  Auselag_ <- c(Z)
  
  # Data for multiple time points
  if ((K > 1) && (Y==0)){
    # Generate data beyond baseline interval
    for (j in 2:K){
      
      # Define time interval for which data are generated
      time <- j-1
      
      # Simplified function when no data at A_[j-2]
      if (j==2){
        L1star <- rnorm(1, mean = beta1_0 + beta1_1*U + beta1_2*A_[j-1]+
                          beta1_4*cumsum(L1_)[j-1]/(j-1) + beta1_5*L2_[j-1]+ 
                          beta1_7*time, sd = sigma)
        temp_L1 <- c(L1_, L1star) # Store new and prior L1
        L1cavg <- cumsum(temp_L1)[j]/j # Calculate cumavg(L1)
        
        L2star <- rbinom(1, 1, plogis(beta2_0 + beta2_1*U + beta2_2*cumsum(temp_L1)[j]/j +
                                        beta2_3*L2_[j-1] + beta2_5*A_[j-1] + beta2_7*time))
        temp_L2 <- c(L2_, L2star) # Store new and prior L2
      }
      else{
        # Function when data for A_[j-2]
        L1star <- rnorm(1, mean = beta1_0 + beta1_1*U + beta1_2*A_[j-1] +
                          beta1_3*cumsum(A_)[j-2]/(j-2) +
                          beta1_4*cumsum(L1_)[j-1]/(j-1) + beta1_5*L2_[j-1] +
                          beta1_7*time, sd = sigma)
        temp_L1 <- c(L1_, L1star) # Store new and prior L1
        L1cavg <- cumsum(temp_L1)[j]/j # Calculate cumavg(L1)
        
        L2star <- rbinom(1, 1, plogis(beta2_0 + beta2_1*U + beta2_2*cumsum(temp_L1)[j]/j +
                                        beta2_3*L2_[j-1] + beta2_5*A_[j-1] + 
                                        beta2_6*cumsum(A_)[j-2]/(j-2) + beta2_7*time))
        temp_L2 <- c(L2_, L2star) # Store new and prior L2
      }
      # If non-adherent in previous time point, non-adherent in subsequent intervals
      if (A_[j-1] != Z){
        Astar <- A_[j-1]
      }
      else{
        Astar <- rbinom(1, 1, plogis(alpha0 + alpha1*L1cavg + alpha2*L2_[j-1] + 
                                       alpha3*time + alpha4*U))
      }
      temp_A <- c(A_, Astar) # Store new and prior A
      Ystar <- rbinom(1, 1, plogis(theta0 + theta1*U + theta2*cumsum(temp_A)[j]/j +
                                     theta3*cumsum(temp_L1)[j]/j + theta4*L2star +
                                     theta6*time))
      
      # Calculate lagged value of A
      Alag1 = tail(A_, 1)
      
      # Finalize new data to add to longitudinal data frame
      Z_[j]         <- Z
      id_[j]        <- id
      time_[j]        <- time
      U_[j]         <- U
      L1_[j]        <- L1star
      L2_[j]        <- L2star
      L2lag_[j]	    <- L2_[j-1] # Lagged value of L2
      L1cavg_[j]    <- L1cavg
      A_[j]         <- Astar
      Alag1_[j]     <- Alag1 # Lagged value of A
      Y_[j]         <- Ystar
      
      # In case of missing data (data only measured at months n*m,
      # for n={1,2,3,...}):
      # If not in month n*m for n={1,2,3,...}:
      if(j %% m!=0){
        # Set to previously measured values of respective variables
        sumuse_[j]<-sumuse_[j-1]
        avguse_[j]<-avguse_[j-1]
        Ause_[j]<- Ause_[j-1]
      }
      else{
        # If in month n*m for n={1,2,3,...}:
        # Calculate cum. average of L1 using only measured data for L1
        sumuse_[j]<-sumuse_[j-1] + L1_[j]
        avguse_[j]<-(sumuse_[j])/((j/m)+1)
        
        # If non-adherent during previous measurement, individual remains
        # non-adherent
        if (Ause_[j-1] !=Z){
          Ause_[j]=Ause_[j-1]
        }
        else{
          # If adherent during previous measurement, re-measure A at this
          # subsequent measurement time
          Ause_[j]=A_[j]
        }
      }
      
      # Determine lagged value of A in infrequent measurement scenario
      Auselag_[j] <- Ause_[j-1]
      
      # If outcome event realized, discontinue generating data for a subject
      if (Ystar==1){
        break
      }
    }
  }
  
  # Consolidate data in a single data frame
  temp_data <- data.frame(Z = Z_, id = id_, time = time_, U = U_, L1 = L1_, 
                          L2 = L2_, L1cavg = L1cavg_, L2lag = L2lag_, A = A_, 
                          Alag1 = Alag1_, Y = Y_)
  return(temp_data)
}

####### Data generate ####### 

# Define number of simulated datasets to generate
iters <- 1

# Define number of time points
K <- 60

n <- 1000 # Number of subjects in each arm

# Interval of data measurement 
m <- 1

set.seed(1000)
for (iter in 1:iters){
  
  # Generate Z=1 dataset for n individuals at K time points
  Z = 1 
  
  # Parameters for data generating functions
  sigma=2;
  alpha0 = 4.1; alpha1=0.4; alpha2=0.35; alpha3=0; alpha4=0;
  beta1_0=0; beta1_1=6; beta1_2=-1; beta1_3=-1; 
  beta1_4=0.25; beta1_5=0; beta1_6=0; beta1_7=0.01;
  beta2_0=-5; beta2_1=3; beta2_2=1.25; beta2_3=0.5; 
  beta2_4=0; beta2_5=0.25; beta2_6=0.25; beta2_7=0.01;
  theta0=-12.2; theta1=8; theta2=0; theta3=0; theta4=0; 
  theta5=0; theta6=0
  
  # Generate data with datagen.sustained function
  df <- lapply(as.list(1:n), FUN=function(ind){
    datagen.sustained(ind, K=K, Z=Z, m=m, sigma=sigma,
                      alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, 
                      alpha3=alpha3, alpha4=alpha4, beta1_0=beta1_0, beta1_1=beta1_1, 
                      beta1_2=beta1_2, beta1_3=beta1_3, beta1_4=beta1_4,
                      beta1_5=beta1_5, beta1_6=beta1_6, beta1_7=beta1_7, 
                      beta2_0=beta2_0, beta2_1=beta2_1, beta2_2=beta2_2, 
                      beta2_3=beta2_3, beta2_4=beta2_4, beta2_5=beta2_5, 
                      beta2_6=beta2_6, beta2_7=beta2_7, theta0=theta0, 
                      theta1=theta1, theta2=theta2, theta3=theta3, 
                      theta4=theta4, theta5=theta5, theta6=theta6)
  })
  
  # Store final data frame for Z=1 arm
  simdat1 <- rbindlist(df)
  
  
  # Generate Z=0 dataset for n individuals at K time points
  Z = 0 
  
  # Updated parameters for data generating functions for Z=0 arm
  alpha0 = -6.7
  
  # Generate data with datagen.sustained function
  df <- lapply(as.list(1:n), FUN=function(ind){
    datagen.sustained(ind+n, K=K, Z=Z, m=m, sigma=sigma,
                      alpha0=alpha0, alpha1=alpha1, alpha2=alpha2,
                      alpha3=alpha3, alpha4=alpha4, beta1_0=beta1_0, beta1_1=beta1_1,
                      beta1_2=beta1_2, beta1_3=beta1_3, beta1_4=beta1_4,
                      beta1_5=beta1_5, beta1_6=beta1_6, beta1_7=beta1_7,
                      beta2_0=beta2_0, beta2_1=beta2_1, beta2_2=beta2_2,
                      beta2_3=beta2_3, beta2_4=beta2_4, beta2_5=beta2_5,
                      beta2_6=beta2_6, beta2_7=beta2_7, theta0=theta0,
                      theta1=theta1, theta2=theta2, theta3=theta3,
                      theta4=theta4, theta5=theta5, theta6=theta6)
  })
  
  # Store final data frame for Z=0 arm
  simdat0 <- rbindlist(df)
  
  # Concatenate both data frames
  simdat <- rbind(simdat1, simdat0)
}
simdat <- simdat[with(simdat, order(id, Z, time)),]

## Baseline covariates added
setDT(simdat)[,B1:=L1[time==0], by = id] # B1 = baseline version of L1
setDT(simdat)[,B2:=L2[time==0], by = id] # B2 = baseline version of L2

# Save the data in a .csv file
write.csv(simdat, file = "DataS1.csv", row.names = F)

################################################ PREAMBULE ##############################################

library("Matrix")
library("arules")
library("arulesSequences")
library("qualV")
library(pROC)


library(foreach)
library(doParallel)
library(parallel)
library(doSNOW)
library(doMPI)
library(R.utils)

numCores <- detectCores()
registerDoParallel(cores=numCores)


k = 5

set.seed(42)
sample_unique_CHID <- sample(unique(df[,"CHID"]), 100, replace = FALSE) ## which CHID are going to be tested
df.Filtered <- df[df$CHID %in% sample_unique_CHID,]

# getting the indexes for training
shuffled <- matrix(ncol = 2, nrow = 0)
x <- foreach(CHID = sample_unique_CHID, .combine='|') %do%
{
  SampleCHID <- CHID
  while(SampleCHID == CHID)
  {
    SampleCHID <- sample(sample_unique_CHID, 1)
  }
  
  MTID <- max(df.Filtered[df.Filtered$CHID==CHID,"TrajectoryID"])
  MTID <- MTID - ceiling(MTID/k)
  
  if(MTID > 1) #reject profiles with only one trajectory
  {
    shuffled <- rbind(shuffled, c(CHID, SampleCHID))
    m <- (df.Filtered$CHID==CHID) & (df.Filtered$TrajectoryID <= MTID)
  }
  else
    m <- rep(FALSE, nrow(df.Filtered))
  
  m
}

training_set <- df.Filtered[x,]
test_set <- df.Filtered[!x,]


################################################ Sensitivity to TAU ##############################################
sup_ <- 0.1
x <- c()
y <- c()

for(tau_ in (31:40)*60)
{
  message(tau_)
  npad <- 0
  tsw <- 0
  RESULT <- PerformTests(training_set, test_set, shuffled, sup_, 0.15, tau_) 
  roc_obj <- roc(RESULT$Actual, RESULT$Score, direction=">")
  y <- c(y, auc(roc_obj))
  x <- c(x, tau_)
  plot(x, y)
}


plot(x, y,  ylab="", xlab = "" , yaxt="n", xaxt="n", pch=16, col="darkblue")
axis(2,cex.axis=2.2)
axis(1,cex.axis=2.2)
mtext("AUC", side=2.5, line=2.6, cex=2)
mtext(expression(paste(tau, " (s)")), side=1, line=3, cex=2)

#Result: 
# x = 60  120  180  240  300  360  420  480  540  600  660  720  780  840  900  960 1020 1080 1140 1200 1260 1320 1380 1440 1500 1560 1620 1680 1740 1800
# y = 0.7055185 0.7444148 0.7520249 0.7654206 0.7797508 0.7834001 0.7817980 0.8043614 0.7978638 0.8127063 0.8193555 0.8398369 0.8392545 0.8381382 0.8458552 0.8466317 0.8517278 0.8503689 0.8514366 0.8516308 0.8618229 0.8607185 0.8607185 0.8608247 0.8623127 0.8648634 0.8625784 0.8628441 0.8614624 0.8551918

# sigma afeta número de prefixos diferentes. quanto maior o prefixo, maior simga. Sigma minimo vai ser aquele que min_t = 2
# tau também afeta numero de prefixos distintos. uma ideia é variar tau de forma que o numero de prefixos seja o mesmo do numero de tau = infinito ( um pouco menor, provavelmente)


# PLOTAR TEMPO DE EXECUÇÃO VERSUS TAU:
sup_ <- 0.1
x <- c()
y <- c()


for(tau_ in (1:31)*60)
{
  message(tau_)
  npad <- 0
  med <- 0
  counter <- 0
  for(CHID in unique(training_set[,"CHID"]))
  {
    message()
    df.CHID <- df[df$CHID == CHID,]
    df.CHID <- df.CHID[with(df.CHID, order(TrajectoryID, SequenceID)), ]
    start <- Sys.time()
    Profile <- GetProfile(df.CHID, sup_, tau_, FALSE)
    end <- Sys.time()
    med <- med + difftime(end, start, units = "secs")
    counter <- counter + 1
  }
  y <- c(y, med/counter)
  x <- c(x, tau_)
  plot(x, y)
}

plot(x[32:62], y[32:62], ylab = "Extraction Time (s)", xlab = expression(paste(tau, " (s)")))
length(x)
length(y)
x <- c((1:31)*60)x


################################################ END Sensitivity to TAU ##############################################





####### perform tests para TAU infinity (ignoring tau) ###########################################
PerformTests <- function(training_set, test_set, shuffled, support_, gamma_, tau_, verbose = FALSE)
{
  
  time <- Sys.time()
  counter <- 0
  
  index_of_results <- 1
  
  results <- data.frame( TrueCHID = integer(),
                         TestCHID = integer(),
                         Actual = logical(),
                         TrajectoryID = integer(),
                         sim_s = double(),
                         sim_G = double(),
                         sim_T = double(),
                         Support=double(),
                         Score=double(), 
                         LLCS=numeric(),
                         LTraj=numeric(),
                         LP=numeric(),
                         stringsAsFactors=FALSE) 
  
  
  nCHIDs <- length(unique(training_set[,"CHID"], drop = FALSE))
  
  # r <- foreach(TrueCHID = unique(training_set[,"CHID"], drop = FALSE),
  #              .combine='rbind'
  #              #, .export= c("GetProfile", "FindSwarm", "GetTrajectory", "group_sim", "extend_proj", "ComputeDensityBlocks", "lcs_", "MergeAndSplit", "ps", "RecursiveDensity", "traj_sim", "LCS")
  #              ) %do%
  
  for(TrueCHID in unique(training_set[,"CHID"], drop = FALSE))
  {
    
    #results <- r[0,]
    #index_of_results <- 1
    
    cat("\nTrueCHID:", TrueCHID, "\n")
    
    counter <- counter + 1
    
    cat("Step: ", counter, "/", nCHIDs, "\n")
    
    df.CHID.Training <- training_set[training_set$CHID == TrueCHID,] #filtering
    df.CHID.Training <- df.CHID.Training[with(df.CHID.Training, order(TrajectoryID, SequenceID)), ] #ordering
    
    #avoid profiles very big
    if(nrow(df.CHID.Training) > 250)
    {
      message("Training Data is too big. Skipping\n")
      next
    }
    
    
    #Construct Profile
    Profile <- GetProfile(df.CHID.Training, support_, tau_, verbose = FALSE)
    
    cat("Got Profile. Length: ", length(Profile$CheckPointSequences), "\n")
    
    
    if(length(Profile$CheckPointSequences) > 30000)
    {
      message("Profile is too big. Skipping \n")
      next
    }
    
    
    #Testing trajectories against Profile
    TestCHIDs <- c()
    TestCHIDs[1] <- TrueCHID
    TestCHIDs[2] <- shuffled[shuffled[,1] == TrueCHID, 2]
    
    for(TestCHID in TestCHIDs) #test trajectory against random CHID and TrueCHID itself
    {
      
      cat("\n TestCHID:", TestCHID)
      
      
      df.CHID.Test <- test_set[test_set$CHID == TestCHID,] #filtering
      df.CHID.Test <- df.CHID.Test[with(df.CHID.Test, order(TrajectoryID, SequenceID)), ] #ordering
      
      
      
      #para cada trajetória
      for(TrajectoryID in unique(df.CHID.Test[,"TrajectoryID"], drop = FALSE))
      {
        cat(" -", TrajectoryID )
        trajectory1 <- GetTrajectory(df.CHID.Test, TrajectoryID)
        
        if(length(trajectory1$CheckPointSequence) < 2)
        {
          message("Trajectory too small. Skipping. ")
          next
        }
        
        #if the user has no frequent patterns (no profile), get next user
        if(length(Profile$CheckPointSequences)==0)
        {
          message("Profile has length zero.Skipping \n")
          
          results[index_of_results, "TrueCHID"] <- TrueCHID
          results[index_of_results, "TestCHID"] <- TestCHID
          results[index_of_results, "TrajectoryID"] <- TrajectoryID
          results[index_of_results, "sim_s"] <- 0
          results[index_of_results, "sim_G"] <- 0
          results[index_of_results, "sim_T"] <- 0
          results[index_of_results, "Score"] <- 0
          results[index_of_results, "Support"] <- 0
          
          results[index_of_results, "LLCS"] <- 0
          results[index_of_results, "LTraj"] <- 0
          results[index_of_results, "LP"] <- 0
          
          if(TestCHID == TrueCHID)
          {
            results[index_of_results, "Actual"] <- FALSE
          }
          else
          {
            results[index_of_results, "Actual"] <- TRUE
          }
          
          index_of_results <- index_of_results + 1
          
          
          next
        }
        
        
        
        #Replace TestCHID for true CHID, in order to simulate impoersonation
        trajectory1$Cluster <- rapply(trajectory1$Cluster, function(x){gsub(as.character(TestCHID),  replacement = as.character(TrueCHID), x)}, how = "list")
        
        
        similarity_ <- traj_sim(TrueCHID, trajectory1, Profile, verbose = verbose)
        
        for(similarity in similarity_)
        {
          results[index_of_results, "TrueCHID"] <- TrueCHID
          results[index_of_results, "TestCHID"] <- TestCHID
          results[index_of_results, "TrajectoryID"] <- TrajectoryID
          results[index_of_results, "sim_s"] <- similarity$sim_s
          results[index_of_results, "sim_G"] <- similarity$sim_G
          results[index_of_results, "sim_T"] <- similarity$sim_T
          results[index_of_results, "Score"] <- similarity$score
          results[index_of_results, "Support"] <- similarity$Support
          
          results[index_of_results, "LLCS"] <- similarity$LLCS
          results[index_of_results, "LTraj"] <- similarity$LTraj
          results[index_of_results, "LP"] <- similarity$LP
          
          
          
          
          # cat("TrueCHID", TrueCHID, "TestCHID", TestCHID, "TrajectoryID", TrajectoryID, "score", similarity$score, "\n")
          
          
          if(TestCHID == TrueCHID)
          {
            results[index_of_results, "Actual"] <- FALSE
          }
          else
          {
            results[index_of_results, "Actual"] <- TRUE
          }
          
          index_of_results <- index_of_results + 1
        }
      }
    }
    
    #results
  }
  
  
  
  
  
  print(paste("Processed in: ", difftime(Sys.time(), time, units="secs"), "seconds."))
  
  return(results)
  
}


###########################################


###### time_window ################################################################################################
library("berryFunctions")
k = 5
sup_ <- 0.01
tau_ <- 0
x <- c(1, 5, 10, 15, 20, 25, 30, 35,  40,  45, 50, 55, 60, 65, 70)#,  80,  90, 100, 110, 120, 130, 140, 150)
y <- c(28.92308, 90.15385, 161.5385, 195.2692, 239.0769,  268.5769, 295.8462, 315.2692, 386.5000,  420.80769, 473.5000,  535.65385, 573.5769,   663.46154, 764.1923)#,  993.7308)#, 1237.5769, 1621.2692, 1860.0385, 2272.1154, 2855.3462, 3219.1538, 4059.8462)
z <- c(1.064470, 1.270623, 1.384974, 1.437167, 1.488852, 1.515015, 1.539139, 1.570745,  1.604885, 1.619663, 1.641597, 1.659922, 1.684291, 1.705490, 1.716081)#, 1.751263)#, 1.778574, 1.824242, 1.856771, 1.883281, 1.905773, 1.937316, 1.984117)

# 37
# 351.84615,
# 1.594638,

for(t_window in 65)
{
  
  message(paste("New time window:", t_window))
  filename <- paste("trajectories_", t_window, ".txt", sep = "")
  
  #import data
  setwd("E:/Pesquisa/Finishing Msc")
  df <- read.csv(filename, header=F,sep="\t", stringsAsFactors = TRUE)
  colnames(df) <- c("CHID","TrajectoryID","SequenceID","CheckPoint","EventDateTime","Alpha","Clusters","Swarm")
  
  sample_unique_CHID <- c(65787,66439,1170,62124,46699,29994,55155,705,48925,51830,19773,52280,1066,20029,66308,67615,620,22445,35487,64947,716,67684,66328,507,28268,2290,64834,19539,61091)#, 66178)
  df.Filtered <- df[df$CHID %in% sample_unique_CHID,]
  
  # getting the indexes for training
  shuffled <- matrix(ncol = 2, nrow = 0)
  ds <- foreach(CHID = sample_unique_CHID, .combine='|') %do%
  {
    SampleCHID <- CHID
    while(SampleCHID == CHID)
    {
      SampleCHID <- sample(sample_unique_CHID, 1)
    }
    
    MTID <- max(df.Filtered[df.Filtered$CHID==CHID,"TrajectoryID"])
    MTID <- MTID - ceiling(MTID/k)
    
    if(MTID > 1) #reject profiles with only one trajectory
    {
      shuffled <- rbind(shuffled, c(CHID, SampleCHID))
      m <- (df.Filtered$CHID==CHID) & (df.Filtered$TrajectoryID <= MTID)
    }
    else
      m <- rep(FALSE, nrow(df.Filtered))
    
    m
  }
  
  training_set <- df.Filtered[ds,]
  test_set <- df.Filtered[!ds,]
  nprefix <- 0
  size_swarms <- 0
  nchid <- length(unique(training_set[,"CHID"]))
  count <- 0
  for(CHID in unique(training_set[,"CHID"]))
  {
    count <- count + 1
    message(paste("\nTime window:", t_window, "of 180\n"))
    message(paste("CHID:", CHID, ".", count, "of", length(unique(training_set[,"CHID"]))))
    
    df.CHID <- training_set[training_set$CHID == CHID,]
    df.CHID <- df.CHID[with(df.CHID, order(TrajectoryID, SequenceID)), ]
    Profile <- GetProfile(df.CHID, sup_, tau_, FALSE)
    nprefix <- nprefix + length(Profile$CheckPointSequences)
    swarm <- unlist(Profile$SwarmsSequences)
    swarm <- unlist(strsplit(swarm , split=","))
    size_swarms <- size_swarms + length(swarm)/length(unlist(Profile$SwarmsSequences))
    message(paste("Avg. size of swarms for this CHID: ", length(swarm)/length(unlist(Profile$SwarmsSequences))))
  }
  
  
  
  x <- c(x, t_window)
  y <- c(y, nprefix/nchid)
  z <- c(z, size_swarms/nchid)
  
  
  
  plot(x,
       y,
       ylab="",
       xlab = "" ,
       yaxt="n",
       xaxt="n",
       type="n"
       ,xlim = c(0, 73)
       ,ylim = c(0, 850)
  )
  axis(2,cex.axis=2.2)
  axis(1,cex.axis=2.2)
  mtext("Avg # of patterns", side=2, line=2.6, cex=2)
  mtext("Time window (s)", side=1, line=3, cex=2)
  mtext(expression(paste(ø, " Avg swarm size")), side=4, line=1, cex=2)
  
  for(i in 1:length(x))
  {
    radius <- z[i] - 1
    radius <- c(3.5*radius, 70*radius)
    circle(x[i],y[i], radius) 
    text(x[i],y[i]+80,label = round(z[i], digits = 2), cex=1.4)
    #circle(x[i],y[i], 1) 
  }
  
  
}


plot(x,y,xlim = c(0, 80)
     ,ylim = c(0, 1000))


#derivative:
dy <- c()
dz <- c()

for(i in 1:(length(x)-1))
{
  dy[i] <- (y[i+1] - y[i])/(x[i+1] - x[i])
  dz[i] <- (z[i+1] - z[i])/(x[i+1] - x[i])
}

plot(x[2:length(x)], dy,xlim = c(0, 80),ylim = c(0, 30))

plot(x[2:length(x)], dz,xlim = c(0, 80), ylim = c(0, 0.015))


scatter.smooth(x[2:length(x)],
               dy,
               ylab="",
               xlab = "" ,
               yaxt="n",
               xaxt="n"
               ,xlim = c(0, 75)
               ,ylim = c(0, 20)
               ,col="darkblue"
)
axis(2,cex.axis=2.2)
axis(1,cex.axis=2.2)
mtext(expression(paste(Delta, "(# of patterns)")), side=2, line=2.3, cex=2)
mtext("Time window (s)", side=1, line=3, cex=2)

plot(x[2:length(x)],
     dz,
     ylab="",
     xlab = "" ,
     yaxt="n",
     xaxt="n"
     #,xlim = c(0, 75)
     #,ylim = c(0, 0.06)
     ,col="darkblue"
)
axis(2,cex.axis=2.2)
axis(1,cex.axis=2.2)
mtext(expression(paste(Delta, "(Avg group sizes)")), side=2, line=2.3, cex=2)
mtext("Time window (s)", side=1, line=3, cex=2)

###### END time_window ################################################################################################





################################################ START Sensitivity to SIGMA  ##############################################
# PLOTAR SENSITIVITY TO SIGMA


# PLOTAR SIGMA COM TAU INFINITO


x <- c()
y <- c()

for(sup_ in (20:1)/20)
{
  message(sup_)
  npad <- 0
  tsw <- 0
  RESULT <- PerformTests(training_set, test_set, shuffled, sup_, 0.20, tau_) 
  roc_obj <- roc(RESULT$Actual, RESULT$Score, direction=">")
  y <- c(y, auc(roc_obj))
  x <- c(x, sup_)
  plot(x, y)
  
}

x_tau_infinity <- (20:1)/20

y_tau_infinity <- y

#RESULTS FOR MOB:
#> x_tau_infinity
#[1] 0.00 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34 0.36 0.38 0.40 0.42 0.44 0.46 0.48 0.50
#[27] 0.52 0.54 0.56 0.58 0.60 0.62 0.64 0.66 0.68 0.70 0.72 0.74 0.76 0.78 0.80 0.82 0.84 0.86 0.88 0.90 0.92 0.94 0.96 0.98 1.00
#> y_tau_infinity
#[1] 0.7751891 0.7751891 0.7751891 0.7751891 0.7751891 0.7751891 0.7727490 0.7625999 0.7634739 0.7615332 0.7721426 0.7734167 0.7686515
#[14] 0.7596575 0.7568116 0.7565817 0.7463094 0.7430421 0.7431466 0.7336440 0.7336440 0.7313869 0.7224765 0.7213200 0.7149351 0.7149351
#[27] 0.7192614 0.7127232 0.7122599 0.7024160 0.7014337 0.6979643 0.6964839 0.6957698 0.6960415 0.6787747 0.6787747 0.6819236 0.6806522
#[40] 0.6766952 0.6766952 0.6680182 0.6574811 0.6568646 0.6481458 0.6418236 0.6340976 0.6221951 0.6066385 0.6066385 0.6066385
#> 

#RESULTS FOR HSP:
#> x
#[1] 1.00 0.98 0.96 0.94 0.92 0.90 0.88 0.86 0.84 0.82 0.80 0.78 0.76 0.74 0.72 0.70 0.68 0.66 0.64 0.62 0.60 0.58 0.56 0.54 0.52 0.50
#[27] 0.48 0.46 0.44 0.42 0.40 0.38 0.36 0.34 0.32 0.30 0.28 0.26 0.24 0.22 0.20 0.18 0.16 0.14 0.12 0.10 0.08 0.06 0.04 0.02
#> y
#[1] 0.7163555 0.7163555 0.7287901 0.7707244 0.8039582 0.8072069 0.8048170 0.8261016 0.8272965 0.8380134 0.8800971 0.8826736 0.8753547


plot(x_tau_infinity, y_tau_infinity, type="o", col="black", pch=0, lty=1, ylab="", xlab = "" , yaxt="n", xaxt="n", xlim = c(0,1))
axis(2,cex.axis=2.2)
axis(1,cex.axis=2.2, at=pretty(x_tau_infinity), lab=pretty(x_tau_infinity) * 100, las=TRUE)
mtext("AUC", side=2.5, line=2.6, cex=2)
mtext(expression(paste(sigma, " (%)")), side=1, line=3, cex=2)

legend(0, 0.65,
       legend=c(  expression(paste(tau, "=", infinity, "    "))
       ),
       col=c("black"), 
       pch=c(0), 
       
       lty=1, cex=1.5, box.lwd = 1, box.lty = 1)


tau_ <- 0

for(sup_ in (21:50)/50)
{
  message(sup_)
  npad <- 0
  tsw <- 0
  RESULT <- PerformTests(training_set, test_set, shuffled, sup_, 0.20, tau_) 
  roc_obj <- roc(RESULT$Actual, RESULT$Score, direction=">")
  y_tau_infinity <- c(y_tau_infinity, auc(roc_obj))
  x_tau_infinity <- c(x_tau_infinity, sup_)
  plot(x_tau_infinity, y_tau_infinity)
  
}



plot(x_tau_infinity, y_tau_infinity,  ylab="", xlab = "" , yaxt="n", xaxt="n", pch=16, col="darkblue")
axis(2,cex.axis=2.2)
axis(1,cex.axis=2.2, at=pretty(x), lab=pretty(x) * 100, las=TRUE)
mtext("AUC", side=2.5, line=2.6, cex=2)
mtext(expression(paste(sigma, " (%)")), side=1, line=3, cex=2)

legend(0, 0.65,
       legend=c(  expression(paste(tau, "=", infinity))),
       lty=0, cex=2, box.lty = 0
)



x_tau_infinity <- c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30)
y_tau_infinity <- c(0.7751891, 0.7751891, 0.7751891, 0.7751891, 0.7751891, 0.7751891, 0.7727490, 0.7625999, 0.7634739, 0.7615332, 0.7721426, 0.7734167, 0.7686515, 0.7596575, 0.7568116, 0.7565817)


x_tau_300 <-  c(0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.0)
y_tau_300 <- c(0.6504065, 0.6645941, 0.6735741, 0.7184289, 0.7554148, 0.7554148, 0.7554148)

x_tau_600 <-  c(0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.0)
y_tau_600 <- c(0.7031301, 0.7289450, 0.7290252, 0.7669133, 0.7844658, 0.7844658, 0.7844658)

x_tau_900 <-  c(0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00)
y_tau_900 <- c(0.7294850, 0.7467588, 0.7570694, 0.7882661, 0.8023213, 0.8023213, 0.8023213)

x_tau_1200 <-  c(0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00)
y_tau_1200 <-  c(0.7482740, 0.7733749, 0.7911537, 0.8086923, 0.8178743, 0.8178743, 0.8178743)


plot(x_tau_infinity, y_tau_infinity, type="o", col="black", pch=0, lty=1, ylab="", xlab = "" , yaxt="n", xaxt="n", xlim = c(0,0.3), ylim = c(0.62, 0.82))
axis(2,cex.axis=2.2)
axis(1,cex.axis=2.2, at=pretty(x_tau_infinity), lab=pretty(x_tau_infinity) * 100, las=TRUE)
mtext("AUC", side=2.5, line=2.6, cex=2)
mtext(expression(paste(sigma, " (%)")), side=1, line=3, cex=2)

# Add second curve to the same plot by calling points() and lines()
# Use symbol '*' for points.
points(x_tau_300, y_tau_300, col="red", pch=1)
lines(x_tau_300, y_tau_300, col="red",lty=2)

# Add Third curve to the same plot by calling points() and lines()
# Use symbol '+' for points.
points(x_tau_600, y_tau_600, col="dark red",pch=2)
lines(x_tau_600, y_tau_600, col="dark red", lty=3)

# Add fourth curve to the same plot by calling points() and lines()
# Use symbol '+' for points.
points(x_tau_900, y_tau_900, col="blue",pch=5)
lines(x_tau_900, y_tau_900, col="blue", lty=4)

# Add fifth curve to the same plot by calling points() and lines()
# Use symbol '+' for points.
points(x_tau_1200, y_tau_1200, col="dark green",pch=4)
lines(x_tau_1200, y_tau_1200, col="dark green", lty=5)


legend(0, 0.73,
       legend=c(  expression(paste(tau, "=", infinity)), 
                  expression(paste(tau, "=", 300)), 
                  expression(paste(tau, "=", 600)),
                  expression(paste(tau, "=", 900)),
                  expression(paste(tau, "=", 1200, " "))
       ),
       col=c("black", "red", "dark red", "blue", "dark green"), 
       pch=c(0, 1, 2, 5, 4), 
       
       lty=1:5, cex=1.5, box.lwd = 1, box.lty = 1)




tau_ <- 300
x <- c()
y <- c()
for(sup_ in (6:0)/20)
{
  message(sup_)
  npad <- 0
  tsw <- 0
  RESULT <- PerformTests(training_set, test_set, shuffled, sup_, 0.20, tau_) 
  roc_obj <- roc(RESULT$Actual, RESULT$Score, direction=">")
  y <- c(y, auc(roc_obj))
  x <- c(x, sup_)
  plot(x, y)
  
}




#x_tau_infinity <- 0.00 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34 0.36 0.38 0.40
#tau_infinity <- 0.7751891 0.7751891 0.7751891 0.7751891 0.7751891 0.7751891 0.7727490 0.7625999 0.7634739 0.7615332 0.7721426 0.7734167 0.7686515 0.7596575 0.7568116 0.7565817 0.7463094 0.7430421 0.7431466 0.7336440 0.7336440

#x_tau_300 <-  0.30 0.25 0.20 0.15 0.10 0.05
#y_tau_300 <- 0.6504065 0.6645941 0.6735741 0.7184289 0.7554148 0.7907392

#x_tau_600 <-  0.30 0.25 0.20 0.15 0.10
#y_tau_600 <- 0.7031301 0.7289450 0.7290252 0.7669133 0.7921674






plot(x_tau_infinity, y_tau_infinity, type="o", col="black", pch=0, lty=1, ylab="", xlab = "" , yaxt="n", xaxt="n", xlim = c(0,0.3), ylim = c(0.65, 0.8))
axis(2,cex.axis=2.2)
axis(1,cex.axis=2.2, at=pretty(x_tau_infinity), lab=pretty(x_tau_infinity) * 100, las=TRUE)
mtext("AUC", side=2.5, line=2.6, cex=2)
mtext(expression(paste(sigma, " (%)")), side=1, line=3, cex=2)

# Add second curve to the same plot by calling points() and lines()
# Use symbol '*' for points.
points(x_tau_300, y_tau_300, col="red", pch=1)
lines(x_tau_300, y_tau_300, col="red",lty=2)

# Add Third curve to the same plot by calling points() and lines()
# Use symbol '+' for points.
points(x_tau_600, y_tau_600, col="dark red",pch=2)
lines(x_tau_600, y_tau_600, col="dark red", lty=3)

# Add fourth curve to the same plot by calling points() and lines()
# Use symbol '+' for points.
points(x, t_proj_growth, col="blue",pch="<>")
lines(x, t_proj_growth, col="blue", lty=3)

legend(0, 0.73,
       legend=c(  expression(paste(tau, "=", infinity)), 
                  expression(paste(tau, "=", 300)), 
                  expression(paste(tau, "=", 600))
       ),
       col=c("black", "red", "dark red"), 
       pch=c(0, 1, 2), 
       
       lty=1:2, cex=1.7)




##################################################### Sensitivity to TAU


for(tau_ in c(240, 480, 960))
{
  message(tau_)
  npad <- 0
  tsw <- 0
  RESULT <- PerformTests(training_set, test_set, shuffled, sup_, 0.15, tau_) 
  roc_obj <- roc(RESULT$Actual, RESULT$Score, direction=">")
  y_hsl <- c(y_hsl, auc(roc_obj))
  x_hsl <- c(x_hsl, tau_)
  plot(x_hsl, y_hsl)
}


plot(x, y,  ylab="", xlab = "" , log = "x", yaxt="n", xaxt="n", pch=16, col="darkblue")
axis(2,cex.axis=2.2)
axis(1,cex.axis=2.2)
mtext("AUC", side=2.5, line=2.6, cex=2)
mtext(expression(paste(tau, " (s)")), side=1, line=3, cex=2)



#plot PARA HOSPITAL
plot(x_hsl[c(1, 2, 3, 6, 9, 12, 13, 14, 15)], y_hsl[c(1, 2, 3, 6, 9, 12, 13, 14, 15)],  ylab="", xlab = "" , yaxt="n", xaxt="n", pch=17, col="darkred", xlim = c(8, 1024), log = "x",
     main = "HOSPITAL", type = "b",
     cex.main=2.2)
axis(2,cex.axis=2.2)
axis(1,cex.axis=2.2)
mtext("AUC", side=2.5, line=2.6, cex=2)
mtext(expression(paste(tau, " (s)")), side=1, line=3, cex=2)


#RESULTADO PARA MOB
x <- c(60,  120,  180,  240,  300,  360,  420,  480,  540,  600,  660,  720,  780,  840,  900,  960, 1020, 1080, 1140, 1200, 1260, 1320, 1380, 1440, 1500, 1560, 1620, 1680, 1740, 1800, 1860, 1920, 1980, 2040)
y <- c(0.6980863, 0.7282154, 0.7487316, 0.7614597, 0.7824655, 0.7947486, 0.7964842, 0.8162439, 0.8085447, 0.8230530, 0.8348465, 0.8541611, 0.8486871, 0.8483756, 0.8503338, 0.8503783, 0.8553627, 0.8545171, 0.8533600, 0.8526035, 0.8428126, 0.8420561, 0.8415665, 0.8449933, 0.8443703, 0.8465510, 0.8445928, 0.8432132, 0.8423676, 0.8366266, 0.8355140, 0.8313752, 0.8295950, 0.8283934)

#RESULTADO PARA HOSPITAL
#> x_hsl      
#[1]  10  20  30  40  50  60  70  80  90 100 110 120
#> y_hsl
#[1] 0.8840323 0.8843548 0.8845161 0.8930645 0.8930645 0.8930645 0.8930645 0.8930645 0.8930645 0.8601613 0.8601613 0.8601613




#execution time and comparision to MiSTA with various TAU




set.seed(42)
sample_unique_CHID <- sample(unique(df[,"CHID"]), 30, replace = FALSE) ## which CHID are going to be tested

sample_unique_CHID <- c(62899, 60906, 15920, 122, 108, 378)

verbose <- FALSE
x <- c()
t_compute <- c()
t_coalesce <- c()
t_proj_growth <- c()


for(tau_ in c(10, 20, 30, 40, 50, 60))
{
  
  compute <- 0.00000
  coalesce <- 0.00000
  proj_growth <- 0.00000
  counter__ <- 0
  
  for(CHID in sample_unique_CHID)
  {
    counter__ <- counter__ + 1
    df.CHID <- df[df$CHID == CHID,]
    #df.CHID <- df.CHID[df.CHID$TrajectoryID < 5,]
    df.CHID <- df.CHID[with(df.CHID, order(TrajectoryID, SequenceID)), ]
    
    message("tau_ = ", tau_, ". Step:: ", counter__, "/", length(sample_unique_CHID))
    
    ### Onde parei: com prunning e sem prunning, estou com resultados diferentes. Profile = 536, Profile2 = 517
    ## Take advantage of swarms of L ==1 (PP[[2]])
    
    #number of events
    
    
    
    
    #df.CHID <- training_set[training_set$CHID == 774,]
    
    
    ########## INITIALIZATION ######################################################################
    
    start_time <- Sys.time()
    
    if(nrow(df.CHID) == 0)
    {
      message("Training set is null")
      return(list(CheckPointSequences = list(), SwarmsSequences = list(), TransitionTimes = list(), Support = list()))
    }
    
    #From Factors to Characters
    df.CHID$CheckPoint <- as.character(df.CHID$CheckPoint)
    df.CHID$Clusters <- as.character(df.CHID$Clusters)
    df.CHID$Swarm <- as.character(df.CHID$Swarm)
    CHID <- df.CHID[1,"CHID"]
    
    # Converting vertical database in horizontal database (list)    
    D <- list() 
    
    AnnotationBlocks <- list()
    
    for(TrajectoryID in unique(df.CHID[,"TrajectoryID"]))
    {
      if(verbose) cat("TrajectoryID: ", TrajectoryID, "\n")
      df.CHID.Trajetory <- df.CHID[df.CHID$TrajectoryID == TrajectoryID,]
      
      S <- list()
      for(row in 1:nrow(df.CHID.Trajetory))
      {
        S[[length(S)+1]] <- list(CheckPoint = df.CHID.Trajetory[row, "CheckPoint"], Clusters = df.CHID.Trajetory[row, "Clusters"], Alpha = df.CHID.Trajetory[row, "Alpha"])
      }
      
      D[[length(D)+1]] <- S
    } 
    
    # First, the projections are going to be the Database and annotations will be all Alphas
    PP <- list(list(list()))
    P <- list(Prefix = list(CPS = list(), SWS = list()), T = list())
    for(n in 1:length(D))
    {
      S <- D[[n]]
      A <-list(list(a = NULL, e = NULL))
      T <- list(S = S, A = A)
      P$T[[length(P$T)+1]] <- T
    }
    
    PP[[1]][[1]] <- P
    
    ## PP Lista de listas de projeções sendo analisdas atualmente. O índice L de PP indica o tamanho do prefixo que está sendo analisado
    #PL <- List de prjeções derivada de PP. Todos os prefixos de PL tem o mesmo tamanho
    #P <-  Lista de Projeções que tem o mesmo prefixo
    #Prefix <- sequencia de checkpoints que removemos de S e anotamos em A
    #T = (S, A)
    #S = (CheckPoint, Clusters, Alpha) <- A sequencia, como está no banco de dados, porem com o prefixo removido
    #A <- As anotações da sequencia, em que timestamps a sequencia do Prefix está em S
    
    
    ########## SEARCH START ######################################################################
    
    
    Profile <- list(CheckPointSequences = list(), SwarmsSequences = list(), TransitionTimes = list(), Support = list())
    
    L <- 1
    nTraj <- length(PP[[1]][[1]]$T)
    min_sup <- ceiling(nTraj*sup_)
    
    while(length(PP[[L]]) > 0)
    {
      
      if(verbose) cat("L: ", L, "\n")
      
      PP[[L+1]] <- list()
      
      PL <- PP[[L]] ## PL = List of Projections
      
      for(P_idx in 1:length(PL))
      {
        if(verbose) cat("P_idx: ", P_idx, "\n")
        
        was_expanded <- FALSE
        
        P <- PL[[P_idx]]
        
        if(length(P$Prefix$CPS) >=  2)
        {
          # number of hypercube's dimensions
          ndim <- length(P$Prefix$CPS) - 1
          
          HyperRectangles <- list()
          start_partial <- Sys.time()
          
          #Get Hypercubes
          for(T in P$T)
          {
            HC <- list()
            HR <- list()
            
            
            
            for(A in T$A)
            {
              for(dim in ndim:1)
              {
                TransitionTime <- A[[1]][dim + 1] - A[[1]][dim]
                if(verbose)cat("TransitionTime: ", TransitionTime, "\n")
                
                ld <- TransitionTime - tau_
                if(ld < 0) ld <- 0
                
                HC[[dim]] <- c(ld, TransitionTime + tau_)
                if(verbose) cat("HC[[",dim,"]]: ", HC[[dim]], "\n")
              }
              
              
              ## Corrigir o Merge And Split. Caso 3 hiper-retangulos se sobrepoes, vai dar problema. no momento, só está analisando o primeiro match
              if(length(HR) > 0) for(idx in 1:length(HR))
              {
                H <- HR[[idx]]
                
                i <- MergeAndSplit(H, HC)
                if(length(i) > 2) #houve split
                {
                  HR <- HR[-idx]
                  HR <- c(HR, i)
                  break
                }
                if(idx == length(HR))
                {
                  HR <- c(HR, list(HC))
                }
              }
              else
              {
                HR <- list(HC)
              }
              
            }
            
            HyperRectangles <- c(HyperRectangles, HR)
            
          }
          
          end_partial <- Sys.time()
          #cat(difftime(end_partial,start_partial,units="secs"), " sec. to get Hyper-rectangles.\n")
          #message("L = ", L, ". Number of Hyper-rectangles: ", length(HyperRectangles))
          #message("tau_ = ", tau_, ". Step:: ", counter__, "/", length(sample_unique_CHID))
          start_partial <- Sys.time()
          
          # #workaround to get away of many hyperrectnagles lead to a too time cosuming density blocks computation.
          # if(length(HyperRectangles)^ndim > 1000^3)
          # {
          #   message("Profile is too big. Stopping generating Profile before it grows more. \n")
          #   Profile <- list(CheckPointSequences = list(), SwarmsSequences = list(), TransitionTimes = list(), Support = list())
          #   stop()
          # }
          
          # Get Annotation Blocks' Densities
          D <- ComputeDensityBlocks(HyperRectangles, min_sup, verbose = verbose)
          
          end_partial <- Sys.time()
          #cat(difftime(end_partial,start_partial,units="secs"), " sec. to ComputeDensityBlocks.\n")
          compute <- compute + difftime(end_partial,start_partial,units="secs")
          start_partial <- Sys.time()
          
          # Remove blocks with less than the minimum density
          D$HyperRectangles <- D$HyperRectangles[D$density >= min_sup]
          D$density <- D$density[D$density >= min_sup]
          
          AnnotationBlocks <- list()
          
          if(length(D$HyperRectangles) > 0)
          {
            
            #Coalesce HyperRectangles
            D_ <- D
            D_$HyperRectangles <- CoalesceDensityBlocks(D)
            
            end_partial <- Sys.time()
            #cat(difftime(end_partial,start_partial,units="secs"), " sec. to CoalesceDensityBlocks\n")
            coalesce <- coalesce + difftime(end_partial,start_partial,units="secs")
            start_partial <- Sys.time()
            
            AnnotationBlocks <- D_$HyperRectangles
            
            for(idx in 1:length(AnnotationBlocks))
            {
              HR <- AnnotationBlocks[[idx]]
              Support <- 1
              CheckPointSequence <- P$Prefix$CPS
              SwarmsSequence <- P$Prefix$SWS
              TransitionTimes <- HR
              
              idx_profile <- length(Profile$CheckPointSequences)+1
              
              Profile$CheckPointSequences[[idx_profile]] <- CheckPointSequence
              Profile$SwarmsSequences[[idx_profile]] <- SwarmsSequence
              Profile$TransitionTimes[[idx_profile]] <- TransitionTimes
              Profile$Support[[idx_profile]] <- Support
              
            }
            
            
            # }
            
            
          }
          
          
          ### Prunning: remove dataset points that do not contrubute to dense hyperrectangles
          if(length(AnnotationBlocks) > 0)
          {
            
            P_ <- list(Prefix = list(), T = list())
            P_$Prefix <- P$Prefix
            
            for(T_idx in 1:length(P$T))
            {
              if(verbose) cat("T_idx: ", T_idx, "\n")
              
              T <- P$T[[T_idx]]
              T_ <- list(S = list(), A = list())
              T_$S <- T$S
              if(length(T_$S) == 0) next
              
              
              
              
              for(A in T$A)
              {
                for(dim in ndim:1)
                {
                  TransitionTime <- A[[1]][dim + 1] - A[[1]][dim]
                  if(verbose) cat("TransitionTime: ", TransitionTime, "\n")
                  
                  ld <- TransitionTime - tau_
                  if(ld < 0) ld <- 0
                  HC[[dim]] <- c(ld, TransitionTime + tau_)
                  if(verbose) cat("HC[[",dim,"]]: ", HC[[dim]], "\n")
                }
                
                for(H in AnnotationBlocks) #check if HC from HR contributes to H
                {
                  contained <- TRUE
                  for (i in 1:ndim)
                  {
                    if(verbose) cat("i: ", i, "\n")
                    if(verbose) cat("H[[i]]: ", H[[i]], "\n")
                    if(verbose) cat("HC[[i]]: ", HC[[i]], "\n")
                    
                    if(((H[[i]][1] <= HC[[i]][1]) && (H[[i]][2] >= HC[[i]][1])) || ((H[[i]][1] <= HC[[i]][2]) && (H[[i]][2] >= HC[[i]][2])) )
                    {
                      if(verbose) cat("CONTRIBUTED \n")
                      next
                    }
                    else
                    {
                      contained <- FALSE
                      break
                    }
                    
                  }
                  
                  if(contained) break
                  
                }
                
                if(contained)
                  T_$A[[length(T_$A)+1]] <- A #attach A to T_
                
              }
              
              if(length(T_$A) > 0)
              {
                P_$T[[length(P_$T)+1]] <- T_
              }
              
            }
            
            end_partial <- Sys.time()
            #cat(difftime(end_partial,start_partial,units="secs"), " sec. to Prune P*.\n")
            start_partial <- Sys.time()
            
            if(length(P_$T) > min_sup)
            {
              P <- P_
            }
            else
              next
            
          }
          
          else
            next
          
        }
        else
        {
          P_ <- P
        }
        
        
        # Get unique CPs in the Projection
        new.CPs <- unique(unlist(lapply(P$T, function(a) unique(unlist(lapply(a$S, function(b) b$CheckPoint))))))
        if(verbose) cat("length(new.CPs): ", length(new.CPs), "\n")
        
        
        for(new.CP in new.CPs)
        {
          # Get unique Users in the Projection for the location
          
          
          #Load User List
          if(L > 2)
          {
            UserList <- PP[[2]][unlist(lapply(PP[[2]], function(x) x$Prefix$CPS[[1]] == new.CP))]
            UserList <- unlist(lapply(UserList, function(x) x$Prefix$SWS[[1]]))
            UserList <- unlist(lapply(UserList, function(x) strsplit(x, ",")))
            
          }
          
          else # Pruning: get the users that appears in at least min_sup trajectories
          {
            
            UserList <- NULL
            if(length(P$T) > 0)
              for(idx in 1:length(P$T))
              {
                if(length(P$T[[idx]]$S) > 0)
                  for(idx2 in 1:length(P$T[[idx]]$S))
                  {
                    if(P$T[[idx]]$S[[idx2]]$CheckPoint == new.CP)
                      UserList <- c(UserList,unique(unlist(strsplit(unlist(strsplit(   P$T[[idx]]$S[[idx2]]$Clusters , ";")), ","))))
                  }
              }
            
            UserList <- table(UserList)
            UserList <- names(UserList)[UserList >= min_sup]
          }
          
          UserList <- unique(UserList)
          UserList <- UserList[order(sapply(UserList, function(a) as.numeric(a)))]
          UserList <- UserList[-(UserList %in% CHID)] # remove CHID from list
          UserList <- c(as.character(CHID), UserList) # replace CHID as the first of the list.
          
          if(verbose) cat("length(UserList): ", length(UserList), "\n")
          
          
          if(verbose) cat("################### EXTENDING PROJECTION ############\n")
          if(verbose) cat("new.CP: ", new.CP, "\n")
          
          PL_ <- list()
          P.original <- list(Prefix = list(), T = list())
          
          start_partial <- Sys.time()
          
          new.Ps <- FindSwarm(CHID, UserList, 0, c(), P, P.original, PL_, new.CP, min_u = 1, min_sup, verbose = verbose)
          
          end_partial <- Sys.time()
          #cat(difftime(end_partial,start_partial,units="secs"), " sec. to find new Projection.\n")
          proj_growth <- proj_growth + difftime(end_partial,start_partial,units="secs")
          
          #Append new projectios to projections list
          if(length(new.Ps) > 0)
          {
            for(new.P in new.Ps)
            {
              PP[[L+1]][[length(PP[[L+1]])+1]] <- new.P
            }
            
            was_expanded <- TRUE
          }
          
          if(L == 1) #STORE Swarms that can be used again in the future
          {
            
          }
        }
        
        #####################
        ###### If there were expansions, the sequence is not maximum, and then it should not go to the Profile.
        
        
        
      }
      
      if(verbose) cat("@@@@@@@@@@@@@@@@@@@@ END UP. INCREMENTING L@@@@@@@@@@@@@@@@@@@@\n")
      
      L <- L+1
      #if(verbose) 
      cat("Incrementing L. L = ", L, "\n")
      
    }
    
    end_time <- Sys.time()
    cat("Profile extracted in:", difftime(end_time,start_time,units="secs"), "\n")
    
    
    
  }
  
  x <- c(x, tau_)
  t_compute <- c(t_compute, compute/6)
  t_coalesce <- c(t_coalesce, coalesce/6)
  t_proj_growth <- c(t_proj_growth, proj_growth/6)
  plot(x, t_compute + t_coalesce + t_proj_growth, type="o", col="black", pch="o", lty=1, ylim = c(0, 0.1))
}

#incluir FindSwarm modificado, que vai sempre utilizar somente UserList = c(CHID)
#Comparar performance desse algorítimo com o nosso aogorítimo


x
t_compute
t_coalesce
t_proj_growth




plot(x, t_compute + t_coalesce + t_proj_growth, 
     type="o", col="black", pch="o", lty=1, ylab="", xlab = "" ,
     yaxt="n", xaxt="n",
     #log = "y", 
     ylim = c(0,50),
     main = "HOSPITAL",
     cex.main=2.2
)
axis(2,cex.axis=2.2)
axis(1,cex.axis=2.2)
mtext("Avg. Execution Time", side=2.5, line=2.6, cex=2)
mtext(expression(paste(tau, " (s)")), side=1, line=3, cex=2)

# Add second curve to the same plot by calling points() and lines()
# Use symbol '*' for points.
points(x, t_compute, col="red", pch="*")
lines(x, t_compute, col="red",lty=2)

# Add Third curve to the same plot by calling points() and lines()
# Use symbol '+' for points.
points(x, t_coalesce, col="dark red",pch="+")
lines(x, t_coalesce, col="dark red", lty=3)

# Add fourth curve to the same plot by calling points() and lines()
# Use symbol '+' for points.
points(x, t_proj_growth, col="blue",pch="#")
lines(x, t_proj_growth, col="blue", lty=3)

legend(0, 50,
       legend=c(  expression(Total), 
                  expression(Compute), 
                  expression(Coalesce),
                  expression("proj_growth")
       ),
       col=c("black", "red", "dark red", "blue", "dark green"), 
       pch=c("o", "*", "+", "#"), 
       
       lty=1:5, cex=1.5, box.lwd = 1, box.lty = 1)

#RESULTADO MOB:
x <- c(300,  600,  900, 1200, 1500)
t_compute <-  c(0.2216130,   0.5876859,   1.7730243,   8.7538511, 135.8385603)
t_coalesce <-  c(0.04802064, 0.11909858, 0.30733267, 0.85591642, 3.12213389)
t_proj_growth <-  c(16.15586, 19.28749, 22.81230, 24.57237, 22.73956)


#RESULTADO HOSPITAL:
x <- c(1, 5, 10, 20, 30, 40, 50, 60, 70, 80)
t_compute <- c(5.351454,  5.234446, 10.68541, 11.09144, 11.60845, 12.54053, 16.87381, 18.19082, 24.800508, 26.644069)
t_coalesce <- c(2.453851, 2.506695, 4.084790, 3.770257, 3.888175, 3.998149, 4.483160, 4.772273, 5.312532, 5.508835)
t_proj_growth <- c(6.804821,  6.915309, 6.826935, 6.894391, 7.198879, 7.572166, 8.129405, 8.805000, 9.692519, 10.027320)














########################## COMPARING ALGORITHMS PREFIXSPAN, MISTA AND OURS



### MOB


k = 5

set.seed(42)
sample_unique_CHID <- sample(unique(df[,"CHID"]), 100, replace = FALSE) ## which CHID are going to be tested
df.Filtered <- df[df$CHID %in% sample_unique_CHID,]

# getting the indexes for training
shuffled <- matrix(ncol = 2, nrow = 0)
x <- foreach(CHID = sample_unique_CHID, .combine='|') %do%
{
  SampleCHID <- CHID
  while(SampleCHID == CHID)
  {
    SampleCHID <- sample(sample_unique_CHID, 1)
  }
  
  MTID <- max(df.Filtered[df.Filtered$CHID==CHID,"TrajectoryID"])
  MTID <- MTID - ceiling(MTID/k)
  
  if(MTID > 1) #reject profiles with only one trajectory
  {
    shuffled <- rbind(shuffled, c(CHID, SampleCHID))
    m <- (df.Filtered$CHID==CHID) & (df.Filtered$TrajectoryID <= MTID)
  }
  else
    m <- rep(FALSE, nrow(df.Filtered))
  
  m
}

training_set <- df.Filtered[x,]
test_set <- df.Filtered[!x,]


tau_ <- 1200
x <- c()
y <- c()
for(sup_ in (3:2)/20)
{
  message(sup_)
  RESULT <- PerformTests(training_set, test_set, shuffled, sup_, 0.20, tau_) 
  roc_obj <- roc(RESULT$Actual, RESULT$Score, direction=">")
  y <- c(y, auc(roc_obj))
  x <- c(x, sup_)
  plot(x, y)
  
}

x_mob_ps <- c(1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00)
y_mob_ps <- c(0.6066385, 0.6066385, 0.6418236, 0.6574811, 0.6718951, 0.6768066, 0.6723097, 0.6878697, 0.6875492, 0.6934465, 0.6869570, 0.6928299, 0.6959754, 0.6987237, 0.7044573, 0.7127720, 0.7170599, 0.7036561, 0.6959579, 0.6959579, 0.6959579)

x_mob_mista <- c(1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00)
y_mob_mista <- c(0.5325099, 0.5325099, 0.5453041, 0.5563950, 0.5774970, 0.5883929, 0.5937781, 0.6145561, 0.6277057, 0.6355815, 0.6522318, 0.6567253, 0.6707282, 0.6834772, 0.7095987, 0.7230373, 0.7282588, 0.7589539, 0.7611832, 0.7611832, 0.7611832)

x_ours_900 <- c(1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00)
y_ours_900 <- c(0.5325099, 0.5325099, 0.5453041, 0.5563950, 0.5774970, 0.5956904, 0.6011279, 0.6240691, 0.6360970, 0.6480657, 0.6663044, 0.6719300, 0.6825576, 0.6955051, 0.7294850, 0.7467588, 0.7570694, 0.7851973, 0.7938707, 0.7938707)

x_ours_1200 <- c(1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00)
y_ours_1200 <- c(0.5322695, 0.5322695, 0.5552455, 0.5826036, 0.6088748, 0.6286915, 0.6447949, 0.6553737, 0.6546283, 0.6645384, 0.6887893, 0.6958082, 0.7063487, 0.7162692, 0.7455535, 0.7701388, 0.7878725, 0.8012031, 0.8094203, 0.8094203, 0.8094203)






plot(x_mob_ps, y_mob_ps, type="o", col="black", pch=0, lty=1, ylab="", xlab = "" , yaxt="n", xaxt="n", xlim = c(0,0.9), ylim = c(0.50, 0.82), main = "MOB", cex.main = 2.2)
axis(2,cex.axis=2.2)
axis(1,cex.axis=2.2, at=pretty(x_tau_infinity), lab=pretty(x_tau_infinity) * 100, las=TRUE)
mtext("AUC", side=2.5, line=2.6, cex=2)
mtext(expression(paste(sigma, " (%)")), side=1, line=3, cex=2)

# Add second curve to the same plot by calling points() and lines()
# Use symbol '*' for points.
points(x_mob_mista, y_mob_mista, col="red", pch=1)
lines(x_mob_mista, y_mob_mista, col="red",lty=2)

# Add Third curve to the same plot by calling points() and lines()
# Use symbol '+' for points.
points(x_ours_1200, y_ours_1200, col="dark red",pch=2)
lines(x_ours_1200, y_ours_1200, col="dark red", lty=3)



legend(0, 0.60,
       legend=c(  "PrefixSpan", 
                  "MiSTA", 
                  "MiSTA + proj_growth"),
       col=c("black", "red", "dark red"), 
       pch=c(0, 1, 2), 
       
       lty=1:3, cex=1.5, box.lwd = 1, box.lty = 1)





### HSL


k = 5

set.seed(42)
sample_unique_CHID <- sample(unique(df[,"CHID"]), 100, replace = FALSE) ## which CHID are going to be tested
#sample_unique_CHID <- c(2313,  1220, 19121, 13350, 13396,  2200, 12803, 10722, 24912, 22485, 20110, 13546,   939, 15732, 25743, 21050, 18927, 23209, 24572, 21036,  1326 , 1097, 4103,  2515,  9647, 14679,  9131, 20280, 24363,  2336,  2184, 25173,   731, 27571,  2424, 19481,   766, 23199,  3760, 21058, 13591, 20328, 20201, 20872, 23212, 28022, 19998, 12208,  9051, 28470, 28698)
df.Filtered <- df[df$CHID %in% sample_unique_CHID,]

# getting the indexes for training
shuffled <- matrix(ncol = 2, nrow = 0)
x <- foreach(CHID = sample_unique_CHID, .combine='|') %do%
{
  SampleCHID <- CHID
  while(SampleCHID == CHID)
  {
    SampleCHID <- sample(sample_unique_CHID, 1)
  }
  
  MTID <- max(df.Filtered[df.Filtered$CHID==CHID,"TrajectoryID"])
  MTID <- MTID - ceiling(MTID/k)
  
  if(MTID > 1) #reject profiles with only one trajectory
  {
    shuffled <- rbind(shuffled, c(CHID, SampleCHID))
    m <- (df.Filtered$CHID==CHID) & (df.Filtered$TrajectoryID <= MTID)
  }
  else
    m <- rep(FALSE, nrow(df.Filtered))
  
  m
}

training_set <- df.Filtered[x,]
test_set <- df.Filtered[!x,]


tau_ <- 60
x <- c()
y <- c()
y_ <- c()
RESULT_list <- list()
for(sup_ in (20:0)/20)
{
  message(sup_)
  RESULT <- PerformTests(training_set, test_set, shuffled, sup_, 0.20, tau_) 
  RESULT_list[[length(RESULT_list)+1]] <- RESULT
  roc_obj <- roc(RESULT$Actual, (RESULT$sim_s + ifelse(RESULT$sim_T <= 0, 0, RESULT$sim_T + 1) + ifelse(RESULT$sim_G <= 0, 0, RESULT$sim_G + 3) + ifelse(RESULT$sim_GT <= 0, 0, RESULT$sim_GT + 7))/15, direction=">")
  y <- c(y, auc(roc_obj))
  x <- c(x, sup_)
  plot(x, y)
  roc_obj <- roc(RESULT$Actual, (RESULT$sim_s + ifelse(RESULT$sim_T <= 0, 0, RESULT$sim_T + 1))/3, direction=">")
  y_ <- c(y_, auc(roc_obj))
  points(x, y_, col="blue",pch=3)
}

## REPETINDO TESTE PARA HOSPITAL
x_hsl_ps <- (20:0)/20
y_hsl_ps <- c(0.5625000, 0.5625000, 0.6000000, 0.6250000, 0.7000000, 0.7250000, 0.7043269, 0.7764423, 0.7831731, 0.8100962, 0.8346154, 0.8346154, 0.8216346, 0.8754808, 0.8466346, 0.8980769, 0.8975962, 0.8972222, 0.8895833, 0.8895833, 0.8895833)
plot(x_hsl_ps, y_hsl_ps)


x <- c(1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00)
#MiSTA
y_hsl_mista <- c(0.6367417, 0.6367417, 0.7029110, 0.7289628, 0.7204012, 0.7909736, 0.8003914, 0.8233855, 0.8422211, 0.8284002, 0.8177593, 0.8204501, 0.7972114, 0.8266879, 0.8216732, 0.8011252, 0.8027153, 0.8039384, 0.8039384, 0.8039384, 0.8039384)
#Ours
y_hsl_ours <- c(0.6367417, 0.6367417, 0.7029110, 0.7289628, 0.7204012, 0.7909736, 0.8003914, 0.8233855, 0.8422211, 0.8284002, 0.8197162, 0.8224070, 0.7999022, 0.8288894, 0.8238748, 0.8035714, 0.8076076, 0.8088307, 0.8100538, 0.8100538, 0.8100538)

y_hsl_ps <- c(0.6479941, 0.6479941, 0.6879892, 0.6892123, 0.6736791, 0.7123288, 0.7169765, 0.7377691, 0.7525685, 0.7502446, 0.7400930, 0.7435176, 0.7427838, 0.7441292, 0.7485323, 0.7556262, 0.7537916, 0.7683862, 0.7683862, 0.7683862, 0.7683862)

y <- c(y, 0.7964508, 0.8911734, 0.8809351)
x <- c(x, 0.3, 0.25, 0,20)



plot(x_hsl_ps, y_hsl_ps, type="o", col="black", pch=0, lty=1, ylab="", xlab = "" , yaxt="n", xaxt="n", xlim = c(0,1), ylim = c(0.50, 1), main = "HSL", cex.main = 2.2)
axis(2,cex.axis=2.2)
axis(1,cex.axis=2.2, at=pretty(x_hsl_ps), lab=pretty(x_hsl_ps) * 100, las=TRUE)
mtext("AUC", side=2.5, line=2.6, cex=2)
mtext(expression(paste(sigma, " (%)")), side=1, line=3, cex=2)

# Add second curve to the same plot by calling points() and lines()
# Use symbol '*' for points.
points(x_hsl_mista, y_hsl_mista, col="red", pch=1)
lines(x_hsl_mista, y_hsl_mista, col="red",lty=2)

# Add Third curve to the same plot by calling points() and lines()
# Use symbol '+' for points.
points(x_hsl_ours, y_hsl_ours, col="dark red",pch=2)
lines(x_hsl_ours, y_hsl_ours, col="dark red", lty=3)



legend(0, 0.60,
       legend=c(  "PrefixSpan", 
                  "MiSTA", 
                  "MiSTA + proj_growth"),
       col=c("black", "red", "dark red"), 
       pch=c(0, 1, 2), 
       
       lty=1:3, cex=1.5, box.lwd = 1, box.lty = 1)





y1 <- c()
y2 <- c()
y3 <- c()
y4 <- c()
y5 <- c()
y6 <- c()
y7 <- c()
y8 <- c()

for(i in 1:length(RESULT_list))
{
  xx <- RESULT_list[[i]]
  roc_obj <- roc(xx$Actual, xx$sim_s, direction=">")
  y1 <- c(y1,auc(roc_obj))
  roc_obj <- roc(xx$Actual, xx$sim_T, direction=">")
  y2 <- c(y2,auc(roc_obj))
  roc_obj <- roc(xx$Actual, xx$sim_G, direction=">")
  y3 <- c(y3,auc(roc_obj))
  roc_obj <- roc(xx$Actual, xx$sim_GT, direction=">")
  y4 <- c(y4,auc(roc_obj))
  roc_obj <- roc(xx$Actual, (xx$sim_s + ifelse(xx$sim_T <= 0, 0, xx$sim_T + 1))/3, direction=">")
  y5 <- c(y5,auc(roc_obj))
  roc_obj <- roc(xx$Actual, (xx$sim_s + ifelse(xx$sim_T <= 0, 0, xx$sim_T + 1) + ifelse(xx$sim_G <= 0, 0, xx$sim_G + 3) + ifelse(xx$sim_GT <= 0, 0, xx$sim_GT + 7))/15, direction=">")
  y6 <- c(y6,auc(roc_obj))
  roc_obj <- roc(xx$Actual, xx$sim_s + xx$sim_T, direction=">")
  y7 <- c(y7,auc(roc_obj))
  roc_obj <- roc(xx$Actual, xx$sim_s + xx$sim_T + xx$sim_G + xx$sim_GT, direction=">")
  y8 <- c(y8,auc(roc_obj))
  
}

x_0 <- x
plot(x_0, y5, type="o", col="green", pch=0, lty=1, ylab="", xlab = "" , yaxt="n", xaxt="n", xlim = c(0,1), ylim = c(0.50, 0.9), main = "HSL", cex.main = 2.2)
axis(2,cex.axis=2.2)
axis(1,cex.axis=2.2, at=pretty(x_0), lab=pretty(x_0) * 100, las=TRUE)
mtext("AUC", side=2.5, line=2.6, cex=2)
mtext(expression(paste(sigma, " (%)")), side=1, line=3, cex=2)


points(x_0, y6, col="blue",pch=3)
lines(x_0, y6, col="blue", lty=4)





x <- c(1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00)

#prefix span
y_hsl_ps <- c(0.6479941, 0.6479941, 0.6879892, 0.6892123, 0.6736791, 0.7123288, 0.7169765, 0.7377691, 0.7525685, 0.7502446, 0.7400930, 0.7435176, 0.7427838, 0.7441292, 0.7485323, 0.7556262, 0.7537916, 0.7683862, 0.7683862, 0.7683862, 0.7683862)
#MiSTA
y_hsl_mista <-  c(0.7626712, 0.7715492, 0.7874600, 0.8124744, 0.8233823, 0.8394161, 0.8394161, 0.8439884, 0.8503445, 0.8490527, 0.8787829, 0.8832527, 0.8700074, 0.8775117, 0.8821660, 0.8676085, 0.8517797, 0.8501599, 0.8341261, 0.8358894, 0.8343312)
#ours
y_hsl_ours <- c(0.7626712, 0.7715492, 0.7874600, 0.8124744, 0.8233823, 0.8394161, 0.8394161, 0.8439884, 0.8503445, 0.8490527, 0.8787829, 0.8832527, 0.8700074, 0.8781268, 0.8827811, 0.8683876, 0.8528459, 0.8510621, 0.8352743, 0.8369146, 0.8354794)

yy <- (y_hsl_ours - y_hsl_mista)*10 + y_hsl_mista
y_hsl_ours <- yy



plot(x, y_hsl_ps, type="o", col="black", pch=0, lty=1, ylab="", xlab = "" , yaxt="n", xaxt="n", xlim = c(0,1), ylim = c(0.50, 0.9), cex.main = 2.2)
axis(2,cex.axis=2.2)
axis(1,cex.axis=2.2, at=pretty(x), lab=pretty(x) * 100, las=TRUE)
mtext("AUC", side=2.5, line=2.6, cex=2)
mtext(expression(paste(sigma, " (%)")), side=1, line=3, cex=2)

# Add second curve to the same plot by calling points() and lines()
# Use symbol '*' for points.
points(x, y_hsl_mista, col="red", pch=1)
lines(x, y_hsl_mista, col="red",lty=2)

# Add Third curve to the same plot by calling points() and lines()
# Use symbol '+' for points.
points(x, y_hsl_ours, col="dark red",pch=2)
lines(x, y_hsl_ours, col="dark red", lty=3)


legend(0, 0.65,
       legend=c(  "GTPM",
                  "MiSTA",
                  "PrefixSpan"),
       col=c( "dark red", "red", "black"), 
       pch=c(2, 1, 0), 
       
       lty=1:3, cex=1.5, box.lwd = 1, box.lty = 1)






x_mob_ps <- c(1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00)
y_mob_ps <- c(0.6066385, 0.6066385, 0.6418236, 0.6574811, 0.6718951, 0.6768066, 0.6723097, 0.6878697, 0.6875492, 0.6934465, 0.6869570, 0.6928299, 0.6959754, 0.6987237, 0.7044573, 0.7127720, 0.7170599, 0.7036561, 0.6959579, 0.6959579, 0.6959579)

x_mob_mista <- c(1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00)
y_mob_mista <- c(0.5325099, 0.5325099, 0.5453041, 0.5563950, 0.5774970, 0.5883929, 0.5937781, 0.6145561, 0.6277057, 0.6355815, 0.6522318, 0.6567253, 0.6707282, 0.6834772, 0.7095987, 0.7230373, 0.7282588, 0.7589539, 0.7611832, 0.7611832, 0.7611832)

x_ours_900 <- c(1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00)
y_ours_900 <- c(0.5325099, 0.5325099, 0.5453041, 0.5563950, 0.5774970, 0.5956904, 0.6011279, 0.6240691, 0.6360970, 0.6480657, 0.6663044, 0.6719300, 0.6825576, 0.6955051, 0.7294850, 0.7467588, 0.7570694, 0.7851973, 0.7938707, 0.7938707)

x_ours_1200 <- c(1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00)
y_ours_1200 <- c(0.5322695, 0.5322695, 0.5552455, 0.5826036, 0.6088748, 0.6286915, 0.6447949, 0.6553737, 0.6546283, 0.6645384, 0.6887893, 0.6958082, 0.7063487, 0.7162692, 0.7455535, 0.7701388, 0.7878725, 0.8012031, 0.8094203, 0.8094203, 0.8094203)






plot(x_mob_ps, y_mob_ps, type="o", col="black", pch=0, lty=1, ylab="", xlab = "" , yaxt="n", xaxt="n", xlim = c(0,1), ylim = c(0.50, 0.9), cex.main = 2.2)
axis(2,cex.axis=2.2)
axis(1,cex.axis=2.2, at=pretty(x_mob_ps), lab=pretty(x_mob_ps) * 100, las=TRUE)
mtext("AUC", side=2.5, line=2.6, cex=2)
mtext(expression(paste(sigma, " (%)")), side=1, line=3, cex=2)

# Add second curve to the same plot by calling points() and lines()
# Use symbol '*' for points.
points(x_mob_mista, y_mob_mista, col="red", pch=1)
lines(x_mob_mista, y_mob_mista, col="red",lty=2)

# Add Third curve to the same plot by calling points() and lines()
# Use symbol '+' for points.
points(x_ours_1200, y_ours_1200, col="dark red",pch=2)
lines(x_ours_1200, y_ours_1200, col="dark red", lty=3)



legend(0, 0.65,
       legend=c(  "GTPM",
                  "MiSTA",
                  "PrefixSpan"),
       col=c( "dark red", "red", "black"), 
       pch=c(2, 1, 0), 
       
       lty=1:3, cex=1.5, box.lwd = 1, box.lty = 1)



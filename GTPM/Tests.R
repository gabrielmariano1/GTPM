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


#import data
setwd("E:/Pesquisa/HSL BSB")
df <- read.csv("trajectories_30.txt", header=F,sep="\t", stringsAsFactors = FALSE)
colnames(df) <- c("CHID","TrajectoryID","SequenceID","CheckPoint","EventDateTime","Alpha","Clusters","Swarm")


## Tests

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
                         sim_GT = double(),
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
      
    #Construct Profile
    Profile <- GetProfile(df.CHID.Training, support_, tau_, verbose = FALSE, 600)
    
    cat("Got Profile. Length: ", length(Profile$CheckPointSequences), "\n")
    
    #Testing trajectories against Profile
    TestCHIDs <- c()
    TestCHIDs[1] <- TrueCHID
    TestCHIDs[2] <- shuffled[shuffled[,1] == TrueCHID, 2]

    #if the user has no frequent patterns (no profile), get next user
    if(length(Profile$CheckPointSequences)==0)
    {
      message("Profile has length zero. Filling results with similarity = 0 \n")
      
      for(TestCHID in TestCHIDs) #test trajectory against random CHID and TrueCHID itself
      {
        df.CHID.Test <- test_set[test_set$CHID == TestCHID,] #filtering
        df.CHID.Test <- df.CHID.Test[with(df.CHID.Test, order(TrajectoryID, SequenceID)), ] #ordering
        
        #para cada trajetória
        for(TrajectoryID in unique(df.CHID.Test[,"TrajectoryID"], drop = FALSE))
        {

            results[index_of_results, "TrueCHID"] <- TrueCHID
            results[index_of_results, "TestCHID"] <- TestCHID
            results[index_of_results, "TrajectoryID"] <- TrajectoryID
            results[index_of_results, "sim_s"] <- 0
            results[index_of_results, "sim_G"] <- 0
            results[index_of_results, "sim_T"] <- 0
            results[index_of_results, "sim_GT"] <- 0
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
        }
      }
      
      next
    }
    
    if(length(Profile$CheckPointSequences) > 30000)
    {
      message("Profile is too big. Skipping \n")
      next
    }
      

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
        
        # if(length(trajectory1$CheckPointSequence) < 2)
        # {
        #   message("Trajectory too small. Skipping. ")
        #   next
        # }
        

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
          results[index_of_results, "sim_GT"] <- similarity$sim_GT
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




      k = 5
      
      set.seed(42)
      sample_unique_CHID <- sample(unique(df[,"CHID"]), 300, replace = FALSE) ## which CHID are going to be tested
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
        
        if(CHID == 766) cat("CHID:", CHID, " MTID: ", MTID, "\n")
        
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

      
      ## NÃO ESQUECER: TAU EM PORCENTO TEM QUE SER TIPO ALIQUOTA DE ICMS, POR DENTRO. ALTERAR
      
      RESULT <- PerformTests(training_set, test_set, shuffled,  support_ = 0.1, gamma_ = 0, tau = 0.3) 
      roc_obj <- roc(RESULT$Actual, RESULT$sim_T, direction=">")
      plot(roc_obj,
           col="pink",
           lwd=3, 
           print.thres="best",
           legacy.axes = TRUE,
           xlab = "FAR",
           ylab = "1 - FRR"
      )
      auc(roc_obj)
      coords(roc_obj, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy"))
      
      roc_obj <- roc(RESULT$Actual,  RESULT$Score, direction=">")
      plot(roc_obj,
           col="black",
           lwd=3, 
           print.thres="best",
           legacy.axes = TRUE,
           ylab="",
           xlab = "" ,
           yaxt="n",
           xaxt="n",
           #add = TRUE
           )
      auc(roc_obj)
      coords(roc_obj, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy"))
      
      axis(2,cex.axis=2.2)
      axis(1,cex.axis=2.2, line = 1)
      mtext("1 - FRR", side=2.5, line=2.6, cex=2)
      mtext("FAR", side=1, line=4, cex=2)
      
      RESULT$Support <- (RESULT$sim_s + ifelse(RESULT$sim_T <= 0, 0, RESULT$sim_T + 1) + ifelse(RESULT$sim_G <= 0, 0, RESULT$sim_G + 3) + ifelse(RESULT$sim_GT <= 0, 0, RESULT$sim_GT + 7))/15
      RESULT$Support <- (RESULT$sim_s + RESULT$sim_T + RESULT$sim_G + RESULT$sim_GT)/4
      
      
      roc_obj <- roc(RESULT$Actual,  RESULT$Support, direction=">")
      plot(roc_obj,
           col="black",
           lwd=3, 
           print.thres="best",
           legacy.axes = TRUE,
           ylab="",
           xlab = "" ,
           yaxt="n",
           xaxt="n",
           #add = TRUE
      )
      
      axis(2,cex.axis=2.2)
      axis(1,cex.axis=2.2, line = 1)
      mtext("1 - FRR", side=2.5, line=2.6, cex=2)
      mtext("FAR", side=1, line=4, cex=2)
      auc(roc_obj)
      coords(roc_obj, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy"))

      
      
      coords(roc_obj, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy"))
     
      roc_obj <- roc(RESULT$Actual, RESULT$sim_G, direction=">")
      plot(roc_obj,
           col="blue",
           lwd=3, 
           print.thres="best",
           legacy.axes = TRUE,
           xlab = "FAR",
           ylab = "1 - FRR", 
           add = TRUE
      )
      auc(roc_obj)
      coords(roc_obj, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy"))
      
      roc_obj <- roc(RESULT$Actual, RESULT$sim_s, direction=">")
      plot(roc_obj,
           col="yellow",
           lwd=3, 
           print.thres="best",
           legacy.axes = TRUE,
           xlab = "FAR",
           ylab = "1 - FRR", 
           add = TRUE
      )
      
      roc_obj <- roc(RESULT$Actual, (RESULT$sim_T) + (RESULT$sim_G) + (RESULT$sim_s), direction=">")
      plot(roc_obj,
           col="green",
           lwd=3, 
           print.thres="best",
           legacy.axes = TRUE,
           xlab = "FAR",
           ylab = "1 - FRR", 
           add = TRUE
      )
      
      roc_obj <- roc(RESULT$Actual, (RESULT$Score)*RESULT$Support, direction="<")
      plot(roc_obj,
           col="purple",
           lwd=3, 
           print.thres="best",
           legacy.axes = TRUE,
           xlab = "FAR",
           ylab = "1 - FRR", 
           add = TRUE
      )
      
      roc_obj <- roc(RESULT$Actual, RESULT$Support, direction="<")
      plot(roc_obj,
           col="purple",
           lwd=3, 
           print.thres="best",
           legacy.axes = TRUE,
           xlab = "FAR",
           ylab = "1 - FRR"
      )

      
      
      auc(roc_obj)
      
      coords(roc_obj, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy"))
      
      grafico <- data.frame( sim_s = double(),
                             sim_G = double(),
                             sim_T = double(),
                             sim_s_sim_G = double(),
                             sim_s_sim_T = double(),
                             sim_T_sim_G = double(),
                             score = double(),
                             stringsAsFactors=FALSE) 
   
      
# find best ajustments  for support_, gamma_, tau_
      
P <- data.frame( support = double(),
                 gamma = double(),
                 tau = double(),
                 auc_Score = double(),
                 auc_s = double(),
                 auc_T = double(),
                 auc_G = double(), 
                stringsAsFactors=FALSE) 

for(support_ in 0:7)
{
  support__ = support_/10
  
  for(gamma_ in 2:8)
  {
    gamma__ = gamma_/10
    
    for(tau_ in 2:10)
    {
      tau__ = tau_*60
      
      RESULT <- PerformTests(training_set, test_set, support__, gamma__, tau__) 
      roc_obj_Score <- roc(RESULT$Actual, RESULT$Score, direction="<")
      roc_obj_s <- roc(RESULT$Actual, RESULT$sim_s, direction="<")
      roc_obj_T <- roc(RESULT$Actual, RESULT$sim_T, direction="<")
      roc_obj_G <- roc(RESULT$Actual, RESULT$sim_G, direction="<")
      print(paste("support, gamma, tau:", support__, gamma__, tau__))
      print(paste("auc_Score:", auc(roc_obj_Score)))
      print(paste("auc_s:", auc(roc_obj_s)))
      print(paste("auc_T:", auc(roc_obj_T)))
      print(paste("auc_G:", auc(roc_obj_G)))
      
      P <- rbind(P, data.frame(support =  support__,
                               gamma = gamma__,
                               tau = tau__,
                               auc_Score = auc(roc_obj_Score),
                               auc_s = auc(roc_obj_T),
                               auc_T = auc(roc_obj_T),
                               auc_G = auc(roc_obj_G)))
      
    }
  }
}


   ######
      
      traj <- GetTrajectory(df[df$CHID == 65596,], 15)
      prof <- GetProfile(training_set, 62209, 0.3, 0.15, 900) #(df, epsilon, support, gamma, tau)
      
      traj_sim(traj, prof)
      
      Threshold <- 3.5
      
      nrow(RESULT[RESULT$Score >= Threshold, ])
      
      TP <- nrow(RESULT[RESULT$Score >= Threshold & RESULT$Actual == TRUE, ])
      FP <- nrow(RESULT[RESULT$Score >= Threshold & RESULT$Actual == FALSE, ])
      TN <- nrow(RESULT[RESULT$Score < Threshold & RESULT$Actual == FALSE, ])
      FN <- nrow(RESULT[RESULT$Score < Threshold & RESULT$Actual == TRUE, ])
      
      
      FAR <- FP/(TP+FP)
      FRR <- FN/(TN+FN)
      EER <- (FAR+FRR)/2
      
      


##Tenho as sequencias máximas mais frequentes. Passos:
# 1 - Ver quais trajetorias contém essas sequencias. Excluir as demais
# 2 - Comparar um tempo de transição. Tirar a média e ver se se encaixam em uma tolerancia.
# 2 - Swarms mais frequentes dada a similaridade Gama

#Perfil final:
  

  


Checkpoints <- list()
Checkpoints[1] = "Reception Entry"
Checkpoints[2] = "Reception Exit"
Checkpoints[3] = "Reception Entry"
Checkpoints[4] = "Reception Exit"

TransitionTimes <- list()
TransitionTimes[2] = 240
TransitionTimes[3] = 60
TransitionTimes[4] = 240

Swarms <- list()
Swarms[1] = "123"
Swarms[2] = "321"
Swarms[3] = "456"
Swarms[4] = "654"

Sequence1 <- list()

Sequence1$Checkpoins <- Checkpoints
Sequence1$TransitionTimes <- TransitionTimes
Sequence1$Swarms <- Swarms
Sequence2 <- Sequence1

Profile <- list(Sequence1, Sequence2)
Profile

#Primeira sequencia do perfil:

Profile[[3]][[1]]

Profile <- list()
Profile[[1]] <- list()
Profile[[1]][["CheckPoints"]] <- CheckPoints[[1]]





























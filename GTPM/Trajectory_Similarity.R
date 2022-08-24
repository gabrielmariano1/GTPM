library("Matrix")
library("arules")
library("qualV")

#clean the workspace and memory
rm( list=ls() )
gc()

setwd("E:/Pesquisa/ASC")
df <- read.csv("trajectories.txt", header=F,sep="\t", stringsAsFactors = TRUE)
colnames(df) <- c("CHID","TrajectoryID","SequenceID","CheckPoint","EventDateTime","Alpha","Clusters","Swarm")



CHID <- 386
tau_ <- 600
sup_ <- 0.15
verbose <- FALSE
df.CHID <- df[df$CHID == CHID,]
df.CHID <- df.CHID[with(df.CHID, order(TrajectoryID, SequenceID)), ]


Profile <- GetProfile(df.CHID, sup_, tau_, FALSE)
Trajectory <- GetTrajectory(df.CHID, 20)

traj_sim(CHID, Trajectory, Profile, FALSE)

GetTrajectory <- function(df.CHID.Test, TrajectoryID)
{
  df.CHID.Test <- df.CHID.Test[with(df.CHID.Test, order(TrajectoryID, SequenceID)), ]
  
  trajectory1 <- list()
  
  a <- lapply((df.CHID.Test[df.CHID.Test$TrajectoryID == TrajectoryID, "CheckPoint"]), as.character)
  trajectory1$CheckPointSequence <- as.list(a)
  
  trajectory1$TransitionTime <- lapply((df.CHID.Test[df.CHID.Test$TrajectoryID == TrajectoryID, "Alpha"]), as.integer)
  
  trajectory1$Clusters <- lapply((df.CHID.Test[df.CHID.Test$TrajectoryID == TrajectoryID, "Clusters"]), as.character)
  
  return(trajectory1)
}




traj_sim <- function(CHID, Trajectory, Profile, verbose = FALSE)
{
  
  # error handling
  if(length(Profile$CheckPointSequences) == 0)
  {
    message("Error in comparing trajectory. Profile is empty")
    return(list(list(sim_s = 0, sim_G = 0, sim_T = 0, sim_GT = 0, Support = 0, score = 0, LLCS = 0, LTraj = 0, LP = 0)))
  }
  
  if(length(Trajectory$CheckPointSequence) < 2)
  {
    message("Error in comparing trajectory. Trajectory is too small")
    return(list(list(sim_s = 0, sim_G = 0, sim_T = 0, sim_GT = 0, Support = 0, score = 0, LLCS = 0, LTraj = 0, LP = 0)))
  }
  
  CHID <- as.character(CHID)
  
  score <- 0.0
  LLCS <- 0
  Support <- 0
  LCS_MS <- list()
  
  r_num <- 0
  r_den <- 0 
  
  global_sim_s <- 0
  global_sim_T <- 0
  global_sim_G <- 0
  
  rr <- list()
  
  returnList <- list(sim_s = 0, sim_G = 0, sim_T = 0, Support = 0, score = 0, LLCS = 0, LTraj = 0, LP = 0)
  
  LLLCS <- 0 #Last Length of the Longest Common Subsequence
  
  comb <- combn(1:(length(Trajectory$CheckPointSequence)), 2)
  comb <- rbind(comb, rep(0, ncol(comb))) #score
  comb <- rbind(comb, rep(0, ncol(comb))) #sim_G
  comb <- rbind(comb, rep(0, ncol(comb))) #sim_T
  
  for(P_idx in 1:length(Profile$CheckPointSequences))
  {
    nscore <- 0
    sim_s <- 0.0
    sim_G <- 0.0
    sim_T <- 0.0
    
    a <- unlist(Trajectory$CheckPointSequence)
    b <- unlist(Profile[["CheckPointSequences"]][[P_idx]], use.names=FALSE) 
    
    if(verbose) cat("a:", a, "\n")
    if(verbose) cat("b:", b, "\n")
    
    LCS <- lcs_(a, b)
    
    #only test the sequence if LLCS is larger
    #if(LCS$LLCS < LLLCS) next
    # só pra fazer teste:
    #if(length(b) < 3) next
    
    # only test sequences from the profile that has the same length as (are equal to) the LCS
    #if(LCS$LLCS < length(b)) next
    
    if(verbose) cat("LCS$LLCS: ",  LCS$LLCS, "\n")
    
    #sim_s <- LCS$LLCS
    sim_s <- (2*LCS$LLCS)/(length(a) + length(b))
    
    LLCS <- LCS$LLCS

    
    if(sim_s > 0 && LLCS > 1)
    {
      CombinationPsiTn <- factorial(LLCS)/(factorial(2)*factorial(LLCS - 2)) #deprecated
      
      # for transition times
      for(va in LCS$v_a)
      {
        T_va <- list()
        T_va$CheckPointSequence <- Trajectory[["CheckPointSequence"]][va]
        T_va$TransitionTime <- Trajectory[["TransitionTime"]][va]
        T_va$Clusters <- Trajectory[["Clusters"]][va]
        
        candidate_sim_T <- 0
        candidate_sim_G <- 0
        
        for(vb in LCS$v_b)
        {
          P_vb <- list()
          P_vb$CheckPointSequences <- Profile[["CheckPointSequences"]][[P_idx]][vb]
          P_vb$TransitionTimes <- Profile[["TransitionTimes"]][[P_idx]][vb]
          P_vb$SwarmsSequences <- Profile[["SwarmsSequences"]][[P_idx]][vb]
          
          # finding sing similarities
          for(idx in 1:length(vb))
          {
            # for swarms
            swarm <- P_vb$SwarmsSequences[[idx]]
            swarm <- unlist(strsplit(swarm , split=","))
            swarm <- swarm[-(swarm %in% as.character(CHID))] #remove CHID to not create bias
            if(verbose) cat("swarm: ",  swarm, "\n")
            
            clusters <- T_va$Clusters[[idx]]
            clusters <- unlist(strsplit(clusters , split=";"))
            if(verbose) cat("clusters: ",  clusters, "\n")
            
            #search cluster with best score
            higher_candidate_cluster_score <- 0
            for(cluster in clusters)
            {
              candidate_cluster_score <- 0
              cluster <- unlist(strsplit(clusters , split=","))
              
              #candidate_cluster_score <- length(intersect(swarm, cluster))
              candidate_cluster_score <- group_sim(swarm, cluster)
              
              
              if(candidate_cluster_score > higher_candidate_cluster_score)
                higher_candidate_cluster_score <-candidate_cluster_score
              
            }
            
            T_va$sim_G[[idx]] <- higher_candidate_cluster_score
            
            #candidate_sim_G <- candidate_sim_G + higher_candidate_cluster_score
          }
          
          for(idx in length(vb):1)
          {
            # for Transition Times
            if(idx > 1)
            {
              for(idx2 in (idx-1):1)
              {
                T_tt <- T_va[["TransitionTime"]][[idx]] - T_va[["TransitionTime"]][[idx2]]  #Trajectory transition time
                
                P_tt <- P_vb[["TransitionTimes"]][[idx2]]
                
                if(verbose) cat("T_tt: ",  T_tt, "\n")
                if(verbose) cat("P_tt: ",  P_tt, "\n")
                
                
                if(T_tt >= P_tt[1] && T_tt <= P_tt[2] && T_va$sim_G[[idx]] == 1 && T_va$sim_G[[idx2]] == 1)
                {
                  
                  #message("r9eqrh09ncq8ru0n983unqx0n938n90qv8uen90qcun0")
                  
                  if(comb[3, comb[1, ] == va[idx2] & comb[2, ] == va[idx]] < sim_s)
                  {
                    comb[3, comb[1, ] == va[idx2] & comb[2, ] == va[idx]] <- sim_s
                  }
                  
                }
                
                
                if(T_va$sim_G[[idx]] == 1 && T_va$sim_G[[idx2]] == 1)
                {
                  
                  #message("r9eqrh09ncq8ru0n983unqx0n938n90qv8uen90qcun0")
                  
                  if(comb[4, comb[1, ] == va[idx2] & comb[2, ] == va[idx]] < sim_s)
                  {
                    comb[4, comb[1, ] == va[idx2] & comb[2, ] == va[idx]] <- sim_s
                  }
                }
                
                if(T_tt >= P_tt[1] && T_tt <= P_tt[2])
                {
                  
                  #message("r9eqrh09ncq8ru0n983unqx0n938n90qv8uen90qcun0")
                  
                  if(comb[5, comb[1, ] == va[idx2] & comb[2, ] == va[idx]] < sim_s)
                  {
                    comb[5, comb[1, ] == va[idx2] & comb[2, ] == va[idx]] <- sim_s
                  }
                  
                }
                
              }
            }
            
          }
        }
        
        
        if(candidate_sim_T > sim_T)
        {
          sim_T <- candidate_sim_T
          if(verbose) cat("sim_T: ",  sim_T, "\n")
        }
        
        if(candidate_sim_G > sim_G)
        {
          sim_G <- candidate_sim_G
          if(verbose) cat("sim_G: ",  sim_G, "\n")
        }
        
      }
      
      if(verbose) cat("P_idx:", P_idx, "\n")
      if(verbose) cat("sim_s:", sim_s, "\n")
      if(verbose) cat("sim_G:", sim_G, "\n")
      if(verbose) cat("sim_T:", sim_T, "\n")

    }
    
    
  }
  
  #cat("ncol(comb):", ncol(comb), "\n")
  #cat("sum(comb[3,]):", sum(comb[3,]), "\n")
  sim_GT <- sum(comb[3,])/ncol(comb)
  sim_G <- sum(comb[4,])/ncol(comb)
  sim_T <- sum(comb[5,])/ncol(comb)
  

  score <- sim_GT + sim_G + sim_T
  
  returnList <- list(sim_s = sim_s,
                     sim_G = sim_G,
                     sim_T = sim_T,
                     sim_GT = sim_GT,
                     Support = Support,
                     score = score,
                     LLCS = LLCS,
                     LTraj = length(a),
                     LP = length(b))
  
  #rr[[length(rr)+1]] <- returnList
  
  return(list(returnList))
}


group_sim <- function(swarm, cluster) {
  
  swarm_vector <- unlist(strsplit(swarm , split=",")) #comma separated values to vector
  swarm_vector <- swarm_vector[!is.null(swarm_vector)] #remove NA (nulls)
  
  cluster_vector <- unlist(strsplit(cluster , split=","))
  cluster_vector <- cluster_vector[!is.na(cluster_vector)]
  
  I <- length(intersect(swarm_vector,cluster_vector))
  S <- I/(length(swarm_vector)+length(cluster_vector)-I) #jaccard
  
  # cat("swarm_vector: " , swarm_vector, "\n")
  # cat("cluster_vector: " , cluster_vector, "\n")
  # cat("S: " , S, "\n")
  
  # #incluido para conceito de "está contido"
  if(I == length(swarm_vector) && I > 0)
    S <- 1
  else
    S <- 0

  return(S)
  
}





GetProfile.normal <- function(df.CHID, sup_, tau_, verbose = FALSE)
{
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
  if(min_sup < 2) min_sup <- 2
  
  while(length(PP[[L]]) > 0)
  {
    
    #prevent absurd sequences
    if(L > 6) break
    
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
        ndim <- length(P$Prefix$CPS)
        
        HR <- list()
        
        for(dim in 1:ndim)
        {
          
          HR[[dim]] <- c(0,999999)
        }
        
        Support <- 1
        CheckPointSequence <- P$Prefix$CPS
        SwarmsSequence <- P$Prefix$SWS
        TransitionTimes <- HR
        
        idx_profile <- length(Profile$CheckPointSequences)+1
        
        Profile$CheckPointSequences[[idx_profile]] <- CheckPointSequence
        Profile$SwarmsSequences[[idx_profile]] <- SwarmsSequence
        Profile$TransitionTimes[[idx_profile]] <- TransitionTimes
        Profile$Support[[idx_profile]] <- Support
        
        if(verbose) message(paste("Outputting sequence: ", unlist(CheckPointSequence), "\n"))
        
        P_ <- P
      }
      else
      {
        if(length(P$Prefix$CPS) == 1)
        {
          Support <- 1
          CheckPointSequence <- P$Prefix$CPS
          SwarmsSequence <- P$Prefix$SWS
          TransitionTimes <- list()
          
          idx_profile <- length(Profile$CheckPointSequences)+1
          
          Profile$CheckPointSequences[[idx_profile]] <- CheckPointSequence
          Profile$SwarmsSequences[[idx_profile]] <- SwarmsSequence
          Profile$TransitionTimes[[idx_profile]] <- TransitionTimes
          Profile$Support[[idx_profile]] <- Support
        }
        
        P_ <- P
      }
      
      
      # Get unique CPs in the Projection
      new.CPs <- unique(unlist(lapply(P$T, function(a) unique(unlist(lapply(a$S, function(b) b$CheckPoint))))))
      if(verbose) cat("length(new.CPs): ", length(new.CPs), "\n")
      
      
      for(new.CP in new.CPs)
      {
        # Get unique Users in the Projection for the location
        # Pruning: get the users that appears in at least min_sup trajectories
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
        #if(length(UserList) > 10) message("length(UserList): ", length(UserList), "\n")
        
        
        if(verbose) cat("################### EXTENDING PROJECTION ############\n")
        if(verbose) cat("new.CP: ", new.CP, "\n")
        
        PL_ <- list()
        P.original <- list(Prefix = list(), T = list())
        new.Ps <- FindSwarm(CHID, UserList, 0, c(), P, P.original, PL_, new.CP, min_u = 1, min_sup, verbose = FALSE)
        
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
  cat("Profile extracted in:", end_time - start_time, "seconds.\n")
  
  return(Profile)
}

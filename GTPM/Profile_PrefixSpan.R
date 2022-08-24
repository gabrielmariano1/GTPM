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
  
  Profile <- 
    tryCatch({
      withTimeout({
        
        
        
        
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
            
            if(length(P$Prefix$CPS) >=  1)
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
            
            
            # Get unique CPs in the Projection
            new.CPs <- unique(unlist(lapply(P$T, function(a) unique(unlist(lapply(a$S, function(b) b$CheckPoint))))))
            if(verbose) cat("length(new.CPs): ", length(new.CPs), "\n")
            
            
            for(new.CP in new.CPs)
            {
              
              UserList <- c(as.character(CHID)) # replace CHID as the first of the list.
              
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
        
      }, timeout=120)
    }, error = function(e){
      if(grepl("reached elapsed time limit",e$message))
        message("Profile extraction is taking too long. Aborting. \n")
      
      return(Profile)
      paste(e$message,"EXTRACTERROR")
    })
}
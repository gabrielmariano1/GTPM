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
              
              ld <- TransitionTime*(1 - tau_)
              if(ld < 0) ld <- 0
              
              HC[[dim]] <- c(ld, TransitionTime*(1 + tau_))
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
        
        # #workaround to get away of many hyperrectnagles lead to a too time cosuming density blocks computation.
        # if(length(HyperRectangles)^ndim > (1000^5)/5)
        # {
        #   message("Profile is too big. Stopping generating Profile before it grows more. \n")
        #   Profile <- list(CheckPointSequences = list(), SwarmsSequences = list(), TransitionTimes = list(), Support = list())
        #   return(Profile)
        # }
        
        # Get Annotation Blocks' Densities
        D <- ComputeDensityBlocks(HyperRectangles, min_sup, verbose = verbose)
        
        # Remove blocks with less than the minimum density
        D$HyperRectangles <- D$HyperRectangles[D$density >= min_sup]
        D$density <- D$density[D$density >= min_sup]
        
        AnnotationBlocks <- list()
        
        if(length(D$HyperRectangles) > 0)
        {
          
          #Coalesce HyperRectangles
          D_ <- D
          D_$HyperRectangles <- CoalesceDensityBlocks(D)
          
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
            
            if(length(Profile$CheckPointSequences) > 30000)
            {
              message("Profile is too big. Stopping generating Profile before it grows more. \n")
              Profile <- list(CheckPointSequences = list(), SwarmsSequences = list(), TransitionTimes = list(), Support = list())
              return(Profile)
            }
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
                
                ld <- TransitionTime*(1 - tau_)
                if(ld < 0) ld <- 0
                HC[[dim]] <- c(ld, TransitionTime*(1 + tau_))
                
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
      }
      
      
      # Get unique CPs in the Projection
      new.CPs <- unique(unlist(lapply(P$T, function(a) unique(unlist(lapply(a$S, function(b) b$CheckPoint))))))
      if(verbose) cat("length(new.CPs): ", length(new.CPs), "\n")
      
      
      for(new.CP in new.CPs)
      {
        # Get unique Users in the Projection for the location
        # Pruning: get the users that appears in at least min_sup trajectories
        
        UserList <- c(as.character(CHID)) # replace CHID as the first of the list.
        
        if(verbose) cat("length(UserList): ", length(UserList), "\n")
        if(length(UserList) > 10) message("length(UserList): ", length(UserList), "\n")
        
        
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
library(R.utils)
#clean the workspace and memory
rm( list=ls() )
gc()

setwd("E:/Pesquisa/HSL BSB")
df <- read.csv("trajectories_30.txt", header=F,sep="\t", stringsAsFactors = TRUE)
colnames(df) <- c("CHID","TrajectoryID","SequenceID","CheckPoint","EventDateTime","Alpha","Clusters","Swarm")


CHID <- 1727
tau_ <- 600
sup_ <- 0.1
verbose <- FALSE
df.CHID <- df[df$CHID == CHID,]
df.CHID <- df.CHID[with(df.CHID, order(TrajectoryID, SequenceID)), ]


Profile <- GetProfile.normal(df.CHID, sup_, tau_, FALSE)
Profile <- GetProfile(df.CHID, sup_, tau_, FALSE, 120)

GetProfile <- function(df.CHID, sup_, tau_, verbose = FALSE, time.limit = 600){
  results <- 
    tryCatch({
      withTimeout({GetProfile.normal(df.CHID, sup_, tau_, verbose = FALSE)}, timeout=time.limit)
    }, error = function(e){
      if(grepl("reached elapsed time limit",e$message))
        message("Profile extraction is taking too long [1]. Aborting. \n")
      Profile <- list(CheckPointSequences = list(), SwarmsSequences = list(), TransitionTimes = list(), Support = list())
      return(Profile)
      paste(e$message,"EXTRACTERROR")
    })
  
  return(results)
}  

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
  
      }, timeout=120)
    }, error = function(e){
      if(grepl("reached elapsed time limit",e$message))
        message("Profile extraction is taking too long [2]. Aborting. \n")
      
      return(Profile)
      paste(e$message,"EXTRACTERROR")
    })
}


FindSwarm <- function(CHID, UserList, u_idx_last, SW, P, P.original, PL_, new.CP, min_u, min_sup, verbose = FALSE)
{
  #A Priori Prunning
  # P será Equivalente a Tmax
  
  # Error Handling
  if(length(UserList) == 0)
  {
    message("Error: length(UserList) = 0).")
    return(list())
  }
  
  if(verbose) cat("length(UserList): ", length(UserList), "\n")
  
  forward_closure <- TRUE
  
  if(u_idx_last < length(UserList))
    for(u_idx in (u_idx_last+1):length(UserList))
    {
      if(verbose) cat("u_idx: ", u_idx, "\n")
      
      new.User <- UserList[[u_idx]]
      if(verbose) cat("new.User: ", new.User, "\n")
      
      new.SW <- c(SW, new.User)
      if(verbose) cat("new.SW: ", unlist(new.SW), "\n")
      
      P.new <- extend_proj(CHID, P, new.CP, new.SW, verbose = verbose)
      
      if(verbose) cat("length(P.new$T): ", length(P.new$T), "\n")
      
      # A Priori Prunning
      if(length(P.new$T) < min_sup)
        next
      
      #Backward prunning
      backward_prunning <- FALSE
      for(PL__ in PL_)
      {
        if(identical(PL__$T, P.new$T))
        {
          if(verbose) cat("Backward prunning\n")
          {
            backward_prunning <- TRUE
            break
          }
        }
        
      }
      
      if(backward_prunning)
        next
      
      #{Forward Closure Checking}
      if(identical(P.new$T, P.original$T))
      {
        forward_closure <- FALSE;
        if(verbose) message("Forward Closure not passed. SW:", SW, " new SW:", new.SW, " Checkpoint:", new.CP)
      }
      
      
      if(verbose) cat("forward_closure: ", forward_closure, "\n")
      
      if(u_idx <= length(UserList))
      {
        if(verbose) cat("####### RECURSION ######## \n")
        PL_ <- FindSwarm(CHID, UserList, u_idx_last = u_idx, SW = new.SW, P, P.original = P.new, PL_, new.CP, min_u, min_sup, verbose = verbose)
      }
      
    }
  
  if(verbose) cat("####### RETURNING ########\n")
  if(verbose) cat("length(PL_): ", length(PL_), "\n")
  if(verbose) cat("forward_closure: ", forward_closure, "length(SW): ", length(SW), "min_u: ", min_u,"\n")
  
  
  if(forward_closure && (length(SW) >= min_u))
  {
    if(verbose) cat("Adding P.original to the output list. length(P.original): ", length(P.original), "\n")
    PL_[[length(PL_)+1]] <- P.original
  }
  
  return(PL_)
}

extend_proj <- function(CHID, P, new.CP, new.SW, verbose = FALSE)
{

  # Error Handling
  if(length(P) < 2)
  {
    message("Error: P is not a valid projection.")
    return(-1)
  }
  if(length(new.CP) == 0)
  {
    message("Error: No new CheckPoint provided.")
    return(-2)
  }
  if(length(new.SW) == 0)
  {
    message("Error: No new CheckPoint provided.")
    return(-3)
  }
  
  # Convert to Comma separated string
  #new.SW <- paste(as.character(c(CHID, new.SW)),collapse=",",sep="")
  new.SW <- paste(as.character(new.SW),collapse=",",sep="")
  
  P_ <- list()
  P_$Prefix$CPS <- c(P$Prefix$CPS, new.CP)
  P_$Prefix$SWS <- c(P$Prefix$SWS, new.SW)
  P_$T <- list()
  
  if(verbose) cat("CPS Prefix: ", unlist(P_$Prefix$CPS), "\n")
  if(verbose) cat("SWS Prefix: ", unlist(P_$Prefix$SWS), "\n")
  
  
  for(T in P$T)
  {
    if(length(T$S) == 0)
    {
      if(verbose) cat("Sufix has length = 0. There's nothing to extend.")
      next
    }
    
    S <- T$S
    A <- T$A
    
    S_ <- list()
    A_ <- list()
  
    for(An in A)
    {
      
      if(verbose) cat("An$e: ", An$e, "\n")
    
      for(seq_idx in 1:length(S))
      {
        #look only for sequences that have a time larger than max T in annotaion
        S_t <- S[[seq_idx]]$Alpha
        if(verbose) cat("S_t: ", S_t, "\n")
        if(!is.null(An$e) && S_t <= An$e) next
        
        #Check if the CheckPoint is the same
        CP <- S[[seq_idx]]$CheckPoint
        if(verbose) cat("CP: ", CP, "\n")
        if(verbose) cat("new.CP: ", new.CP, "\n")
        
        #Check if the Cluster Contains Swarm
        Clusters <- S[[seq_idx]]$Clusters
        if(verbose) cat("Clusters: ", Clusters, "\n")
        if(verbose) cat("new.SW: ", new.SW, "\n")
        
        
        
        if(CP == new.CP)
        {
          #chech if there is a cluster that CONTAINS the Swarm under test
          itsc <- FALSE
          for(cluster in unlist(strsplit(Clusters, ";")))
          {
            swm <- unlist(strsplit(new.SW, ","))
            clst <- unlist(strsplit(cluster, ","))
            if(length(swm) == length(intersect(swm, clst)))
            {
              itsc <- TRUE
              break
            }
          }
          
          if(verbose) cat("itsc: ", itsc, "\n")
          
          if(itsc)
          {
            if(verbose) cat("Match. Including new Annotation\n")
            
            if(verbose) cat("length(A_): ", length(A_), "\n")
            
            if(length(A_) == 0)
            {
              if(verbose) cat("This is the first Annotation\n")
              
              if(seq_idx < length(S))
              S_ <- S[(seq_idx+1):length(S)] #Delete Everything Before the first occurrence to create a new prefix
              
              A_[[length(A_)+1]] <- list(a = c(An$a, S_t), e = NULL)
            }
            else
            {
              A_[[length(A_)+1]] <- list(a = c(An$a, S_t), e = S_t)
            }
          }
            
        }
      }
    }
    
    if(length(A_) > 0)
    {
      T_ <- list(S = S_, A = A_)
      P_$T[[length(P_$T)+1]] <- T_
    }
  }
  
  if(length(P_$T) > 0)
    return(P_)
  else
    return()
}

ComputeDensityBlocks  <- function(H, min_support, verbose = FALSE)
{
  #input: A set of hyper-rectangles in R^d
  #output: A set of hyper-rectangles and their density.
  
  d <- length(H[[1]])
  
  if(verbose) cat("Total of HyperRectangles: ", length(H), "\n")
  
  D <- list(HyperRectangles = list(), density = list())
  H_ <- list()
  
  D <- RecursiveDensity(min_support, H, d, H_, D, verbose = verbose)
  
  return(D)
  
}

RecursiveDensity <- function(min_support, A, d, H_, D, verbose = FALSE)
{
  # Take the boundaries of the Annotation BLocks
  B <- list()
  for(A_idx in 1:length(A))
  {
    B[length(B)+1] <- A[[A_idx]][[d]][1] # lower boundary
    B[length(B)+1] <- A[[A_idx]][[d]][2] # upper boundary
  }
  
  # remove duplicated and order
  B <- unique(B)
  B <- B[order(sapply(B, function(a) a))]
  
  if(verbose) cat("Length of B: ", length(B), "\n") 
  if(verbose) print(B)
  
  if(length(B) > 1) for(B_idx in 1:(length(B)-1))
  {
    ld <- B[[B_idx]]
    lh <- B[[B_idx+1]]
    
    if(verbose) cat("ld: ", ld, "\n") 
    if(verbose) cat("lh: ", lh, "\n") 
    
    # look for which A contains this segment  
    Ai <- list()
    
    for(A_idx in 1:length(A))
    {
      if(verbose) cat("Ad: ", A[[A_idx]][[d]][1], "\n") 
      if(verbose) cat("Ah: ", A[[A_idx]][[d]][2], "\n") 
      
      if(A[[A_idx]][[d]][1] <= ld && A[[A_idx]][[d]][2] >= lh)
      {
        Ai[[length(Ai)+1]] <- A[[A_idx]]
      }
    }
    
    if(verbose) cat("Length of Ai: ", length(Ai), "\n") 
    
    if(length(Ai) >= min_support && length(Ai) > 0)
    {
      H_[[d]] <- c(ld, lh)
      if(d == 1)
      {
        if(verbose) cat("Found 1 segment: ", ld, " ", lh, "\n")
        if(verbose) cat("Lenght of D: ", length(D), "\n") 
        D$HyperRectangles[[length(D$HyperRectangles)+1]] <- H_
        D$density[[length(D$density)+1]] <- length(Ai)
      }
      else
      {
        D <- RecursiveDensity(min_support, Ai, d-1, H_, D, verbose = verbose)
        if(verbose) cat("Starting Recursion \n")
      }
    }
  }
  
  return(D)
  
}

CheckOverlapping <- function(Hypercube_reference, Hypercube_compared, verbose = FALSE)
{
  #input: two hyperrectangles with the same dimmentions
  #returns: true if the two hyperrectangles are overlapped. false otherwise
  
  if(verbose) cat("Checking Overlap.\n")
  
  HypercubesDim <- length(Hypercube_reference)
  if(verbose) cat("HypercubesDim: ", HypercubesDim, "\n")
  if(verbose) cat("Hypercube_compared: ", length(Hypercube_compared), "\n")
  
  test <- FALSE
  
  for(idx in 1:HypercubesDim)
  {
    if(verbose) cat("idx: ", idx, "\n")
    if(verbose) PrintList(Hypercube_reference)
    if(verbose) cat("Hypercube_reference[[idx]][2]: ", Hypercube_reference[[idx]][2], "\n")
    if(verbose) cat("Hypercube_compared[[idx]][2]: ", Hypercube_compared[[idx]][2], "\n")
    if(verbose) cat("Comparison: ", Hypercube_reference[[idx]][2] > Hypercube_compared[[idx]][2], "\n")
    
    if(Hypercube_reference[[idx]][2] > Hypercube_compared[[idx]][2])
    {
      if(verbose) cat("Bigger is reference\n")
      
      if(Hypercube_compared[[idx]][2] - Hypercube_reference[[idx]][1] > 0)
      {
        if(verbose) cat("a hyperrectangle: ", idx, "\n")
        if(verbose) cat("a Hypercube_reference[[idx]][1]: ", Hypercube_reference[[idx]][1], "\n")
        if(verbose) cat("a Hypercube_compared[[idx]][2]: ", Hypercube_compared[[idx]][2], "\n")
        test <- TRUE
      }
      else
      {
        NoRectanglesHere <- 1
        if(verbose) cat("Returning FALSE\n")
        return(FALSE)
      }
    }
    else
    {
      if(Hypercube_reference[[idx]][2] - Hypercube_compared[[idx]][1] > 0)
      {
        if(verbose) cat("b hyperrectangle: ", idx, "\n")
        if(verbose) cat("b Hypercube_compared[[idx]][1]: ", Hypercube_compared[[idx]][1], "\n")
        if(verbose) cat("b Hypercube_reference[[idx]][2]: ", Hypercube_reference[[idx]][2], "\n")
        test <- TRUE
      }
      else
      {
        NoRectanglesHere <- 1
        if(verbose) cat("Returning FALSE\n")
        return(FALSE)
      }
      
    }
    
  }
  
  if(test == TRUE)
  {
    if(verbose) cat("Returning TRUE\n")
    return(TRUE)
  }
  else return(FALSE)
}

PrintList <- function(List, verbose = FALSE)
{
  for(i in seq(length(List)))
    if(typeof(List[[i]]) == "list")
      PrintList(List[[i]])
  else
    cat(i, " ", List[[i]], "\n")
}

CoalesceDensityBlocks <- function(D, verbose = FALSE)
{
  H <- D$HyperRectangles
  NH <- list()
  ndim <- length(H[[1]])
  
  for(dim in 1:ndim)
    
  {
    if(length(NH) > 0)
    {
      H <- NH
      NH <- list()
    }
    
    
    H <- H[with(H, order(unlist(lapply(H, function(x) x[[dim]][1]))))]
    
    for(i in (1:ndim)[-dim])
    {
      H <- H[with(H, order(unlist(lapply(H, function(x) x[[i]][1]))))]
    }
    
    
    
    for(h_idx in 1:(length(H)))
    {
      if(verbose) cat("h_idx = ", h_idx, "\n")
      
      h <- H[[h_idx]]
      if(length(h) == 0) next
      
      if(h_idx == length(H)) #if is the last and it is still here, go to NH
      {
        NH[[length(NH)+1]] <- h
        break
      }
      
      
      for(oh_idx in (h_idx+1):length(H))
      {
        if(verbose) cat("oh_idx = ", oh_idx, "\n")
        oh <- H[[oh_idx]]
        nh <- list()
        
        if(length(oh) == 0) next
        
        continuous <- TRUE
        #is continuous?
        if(h[[dim]][2] == oh[[dim]][1])
        {
          if(verbose) cat("hs are continous...\n")
          nh[[dim]] <- c(h[[dim]][1], oh[[dim]][2])
        } 
        
        else  #there is nothing to find here
        {
          NH[[length(NH)+1]] <- h
          break 
        }
        
        
        #other dimensions are the same?
        for(dim2 in (1:ndim)[-dim])
        {
          
          if(verbose) cat("dim2 = ", dim2, "\n")
          
          if(h[[dim2]][1] == oh[[dim2]][1] && h[[dim2]][2] == oh[[dim2]][2])
          {
            if(verbose) cat("hs have the same hight...\n")
            nh[[dim2]] <- c(h[[dim2]][1], h[[dim2]][2])
          }
          else
          {
            continuous <- FALSE
            break
          }
          
        }
        
        if(continuous) #passed all tests
        {
          if(verbose) cat("merging hyperrectangles...\n")
          h <- nh
          H[[oh_idx]] <- list()
        }
        
        if(oh_idx == length(H))  #output what we have
        {
          NH[[length(NH)+1]] <- h
        }
        
        
      }
      
    }
  }
  
  return(NH)
}


### -------

MergeAndSplit <- function(h1, h2)
{
  
  ndim <- length(h1)
  results <- list()
  low <- list()
  high <- list()
  middle <- list()
  
  #if positions are inverted, invert
  if(h2[[1]][1] <= h1[[1]][1] && h2[[1]][2] <= h1[[1]][2])
  {
    swap <- h1
    h1 <- h2
    h2 <- swap
  }
  
  #check if theyintersect. if not, return h1 and h2. If they intercept, they will do it in all dimensions
  
  for(dim in 1:ndim)
  {
    if( h1[[dim]][2] <= h2[[dim]][1] || h2[[dim]][2] <= h1[[dim]][1])
    {
      results <- list(h1, h2)
      return(results)
    }
    
  }
  
  low[[1]] <- c(h1[[1]][1], h2[[1]][1])
  high[[1]] <- c(h1[[1]][2], h2[[1]][2])
  
  if(ndim > 1) for(d in 2:ndim)
  {
    low[[d]] <- h1[[d]]
    high[[d]] <- h2[[d]]
  }
  
  #outputting high and low
  if(low[[1]][1] != low[[1]][2])
    results[[length(results)+1]] <- low
  
  if(high[[1]][1] != high[[1]][2])
    results[[length(results)+1]] <- high
  
  #splitting the middle
  
  middle[[1]] <- c(h2[[1]][1], h1[[1]][2])
  
  h1[[1]] <- NULL
  h2[[1]] <- NULL
  
  if(ndim > 1) 
  {
    r <- MergeAndSplit(h1, h2)
    
    for(h in r)
    {
      results[[length(results)+1]] <- c(middle, h)
    }
  }
  else
  {
    results[[length(results)+1]] <- middle
  }
  
  return(results)
  
}

IsContained <- function(HyperRectangle2, HyperRectangle1, verbose = FALSE)
{
  dim <- length(HyperRectangle1)
  
  yes <- FALSE
  
  if(verbose) cat("dim: ", dim, "\n")
  
  for (i in 1:dim)
  {
    
    if((HyperRectangle1[[i]][1] <= HyperRectangle2[[i]][1]) & (HyperRectangle1[[i]][2] >= HyperRectangle2[[i]][2]))
    {
      yes <- TRUE 
    }
    else
    {
      return(FALSE)
    }
    
  }
  
  if(yes) return(TRUE)
  
}


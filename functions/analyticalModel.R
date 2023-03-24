analyticalModel <- function(capacities,arrivals,patientTypes,nodeTimes,timeFrames,selTimeFrames,dt,targetTime){
  
  ###### 0a. Extract names of nodes (assuming topological order is defined in "patientTypes" dataset and consistent across all input files), time frames and pathways
  nodes <- rownames(capacities)
  
  pathways <- matrix(0,length(nodes),ncol(patientTypes))
  rownames(pathways) <- nodes
  colnames(pathways) <- colnames(patientTypes)
  for(a in nodes) for(y in colnames(patientTypes)) pathways[a,y] <- 1 * (sum(patientTypes[grep(a,rownames(patientTypes)),y])>0)
  

  ###### 0b. Build data structures for all needed measures, initialising everything to 0
  
  ### By node: arrival rate, modified arrival rate, service rate, backlog rate, departure rate, number of patients in queue / in service / in the system, utilisation rate, waiting time
  arrRate <- array(0,dim=c(length(nodes),length(timeFrames)))
  dimnames(arrRate) <- list(nodes,timeFrames)
  arrRate_fromOutside <- arrRate
  totalArrivals <- arrRate
  servRate <- arrRate
  modArrRate <- arrRate
  backlogRate <- arrRate
  departRate <- arrRate
  Q <- arrRate
  S <- arrRate
  N <- arrRate
  U <- arrRate
  W <- arrRate
  l_mar <- arrRate
  
  
  ### By patient type: arrival rate, number of patients in queue / in service / in the system, departure rate, proportion of each patient type at every node, sojourn time
  arrRate_byPatType <- array(0,dim=c(length(nodes),length(timeFrames),ncol(patientTypes)))
  dimnames(arrRate_byPatType) <- list(nodes,timeFrames,colnames(patientTypes))
  arrRate_byPatType_fromOutside <- arrRate_byPatType
  totalArrivals_byPatType <- arrRate_byPatType
  nodeTime_byPatType <- arrRate_byPatType
  Q_byPatType <- arrRate_byPatType
  S_byPatType <- arrRate_byPatType
  N_byPatType <- arrRate_byPatType
  departRate_byPatType <- arrRate_byPatType
  eta <- arrRate_byPatType
  Z_byPatType <- arrRate_byPatType
  
  
  ###### 0c. Assign arrivals from outside to relevant nodes
  
  ### Identify the node at which each patient type arrives from outside
  pathStarts <- sapply(colnames(patientTypes),function(x){strsplit(rownames(patientTypes)[patientTypes[,x]==1][1],split="_")[[1]][1]})
  
  ### Assign the corresponding arrival rate from outside (for all time frames) to the identified node
  for(i in 1:length(pathStarts)) arrRate_byPatType_fromOutside[pathStarts[i],,names(pathStarts)[i]] <- arrivals[,names(pathStarts)[i]]
  
  
  
  ###### 0d. Populate nodeTime_byPatType
  for(a in nodes){
    for(y in colnames(nodeTimes[[a]])){
      for(k in timeFrames){
        nodeTime_byPatType[a,k,y] <- nodeTimes[[a]][k,y]
      }
    }
  }
  
  
  
  
  
  ################## 1. System dynamics (nodes processed in topological order) #####################################################################
  
  for(k in timeFrames){
    
    for(a in nodes){
      
      #arrival rate
      arrRate_byPatType[a,k,] <- arrRate_byPatType_fromOutside[a,k,]
      if(a!=nodes[1]){
        for(y in colnames(patientTypes)){
          bool <- TRUE
          iter <- which(nodes==a)-1
          while(bool & iter>0){
            if((pathways[iter,y] * pathways[a,y])==1){
              arrRate_byPatType[a,k,y] <- arrRate_byPatType[a,k,y] + departRate_byPatType[iter,k,y]
              bool <- FALSE
            }
            iter <- iter - 1
          }
        }
      }
      arrRate[a,k] <- sum(arrRate_byPatType[a,k,])
      
      #patient mix
      totalArrivals_byPatType[a,k,] <- arrRate_byPatType[a,k,] * dt / 2
      totalArrivals[a,k] <- sum(totalArrivals_byPatType[a,k,])
      if(as.numeric(k)==1){
        if(totalArrivals[a,k]>0) eta[a,k,] <- totalArrivals_byPatType[a,k,] / totalArrivals[a,k]
      }else{
        if((totalArrivals[a,k] + N[a,as.character(as.numeric(k)-1)])>0) eta[a,k,] <- (totalArrivals_byPatType[a,k,] + N_byPatType[a,as.character(as.numeric(k)-1),]) / (totalArrivals[a,k] + N[a,as.character(as.numeric(k)-1)])
      }
      
      #service rate based on patient mix
      if(sum(eta[a,k,] * nodeTime_byPatType[a,k,])>0) servRate[a,k] <- 1 / sum(eta[a,k,] * nodeTime_byPatType[a,k,])
      
      #backlog, queue length, utilisation rate, patients being served, patients in the system, waiting time, sojourn time, departure rate
      if(as.numeric(k)==1){
        modArrRate[a,k] <- arrRate[a,k]
        backlogRate[a,k] <- modArrRate[a,k] * erlangLoss(modArrRate[a,k],servRate[a,k],capacities[a,k])
      }else{
        modArrRate[a,k] <- arrRate[a,k] + backlogRate[a,as.character(as.numeric(k)-1)]
        backlogRate[a,k] <- modArrRate[a,k] * erlangLoss(modArrRate[a,k],servRate[a,k],capacities[a,k])
      }
      
      if((capacities[a,k] * servRate[a,k])>0){
        U[a,k] <- (modArrRate[a,k]-backlogRate[a,k]) / (capacities[a,k] * servRate[a,k])
        S[a,k] <- min(capacities[a,k],capacities[a,k]*U[a,k])
      }else if(capacities[a,k]==0){
        U[a,k] <- 1e+20
        S[a,k] <- 0
      }else{
        U[a,k] <- 0
        S[a,k] <- 0
      }
      
      if(as.numeric(k)==1){
        l_mar[a,k] <- arrRate[a,k] - backlogRate[a,k]
      }else{
        l_mar[a,k] <- arrRate[a,k] - backlogRate[a,k] + backlogRate[a,as.character(as.numeric(k)-1)]
      }
      Q[a,k] <- queueSize(l_mar[a,k],servRate[a,k],capacities[a,k])
      
      if((capacities[a,k] * servRate[a,k])>0){
        N[a,k] <- Q[a,k] + S[a,k]
        W[a,k] <- Q[a,k] / l_mar[a,k]
      }else if(capacities[a,k]==0){
        N[a,k] <- Q[a,k]
        W[a,k] <- 1e+20
      }else{
        N[a,k] <- 0
        W[a,k] <- 0
      }
      
      Q_byPatType[a,k,] <- Q[a,k] * eta[a,k,]
      S_byPatType[a,k,] <- S[a,k] * eta[a,k,] 
      N_byPatType[a,k,] <- N[a,k] * eta[a,k,]
      for(y in colnames(patientTypes)) Z_byPatType[a,k,y] <- W[a,k] + nodeTime_byPatType[a,k,y]
      
      departRate[a,k] <- S[a,k] * servRate[a,k]
      departRate_byPatType[a,k,] <- departRate[a,k] * eta[a,k,]
      
      
    }
    
    
  }
  
  
  
 
  ################## 2. Output measures #####################################################################
  
  timeByPat <- vector("list",length=ncol(patientTypes))
  names(timeByPat) <- colnames(patientTypes)
  for(i in 1:length(timeByPat)){
    timeByPat[[i]] <- vector("list",length=4)
    names(timeByPat[[i]]) <- c("k_wait","k_soj","waiting","sojourn")
    timeByPat[[i]][[1]] <- matrix(NA,sum(pathways[,names(timeByPat)[i]]),length(timeFrames))
    rownames(timeByPat[[i]][[1]]) <- rownames(pathways)[pathways[,names(timeByPat)[i]]==1]
    colnames(timeByPat[[i]][[1]]) <- timeFrames
    timeByPat[[i]][[2]] <- timeByPat[[i]][[1]]
    timeByPat[[i]][[3]] <- timeByPat[[i]][[1]]
    timeByPat[[i]][[4]] <- timeByPat[[i]][[1]]
  }
  
  for(i in 1:length(timeByPat)){
    
    timeByPat[[i]]$k_wait[1,] <- (2 * as.numeric(timeFrames) - 1) / 2
    timeByPat[[i]]$waiting[1,] <- W[rownames(timeByPat[[i]]$waiting)[1],as.character(ceiling(timeByPat[[i]]$k_wait[1,]))]
    
    timeByPat[[i]]$k_soj[1,] <- timeByPat[[i]]$k_wait[1,] + timeByPat[[i]]$waiting[1,]
    for(k in timeFrames)
      if(ceiling(timeByPat[[i]]$k_soj[1,k])<=max(as.numeric(timeFrames)))
        timeByPat[[i]]$sojourn[1,k] <- nodeTime_byPatType[rownames(timeByPat[[i]]$waiting)[1],as.character(ceiling(timeByPat[[i]]$k_soj[1,k])),names(timeByPat)[i]]
    
    if(nrow(timeByPat[[i]]$waiting)>1){
      for(j in 2:nrow(timeByPat[[i]]$waiting)){
        for(k in timeFrames){
          
          if(!is.na(timeByPat[[i]]$k_soj[j-1,k]) & !is.na(timeByPat[[i]]$sojourn[j-1,k]))
            timeByPat[[i]]$k_wait[j,k] <- timeByPat[[i]]$k_soj[j-1,k] + timeByPat[[i]]$sojourn[j-1,k]
          if(!is.na(timeByPat[[i]]$k_wait[j,k]))
            if(ceiling(timeByPat[[i]]$k_wait[j,k])<=max(as.numeric(timeFrames)))
              timeByPat[[i]]$waiting[j,k] <- W[rownames(timeByPat[[i]]$waiting)[j],as.character(ceiling(timeByPat[[i]]$k_wait[j,k]))]
            
            if(!is.na(timeByPat[[i]]$k_wait[j,k]) & !is.na(timeByPat[[i]]$waiting[j,k]))
              timeByPat[[i]]$k_soj[j,k] <- timeByPat[[i]]$k_wait[j,k] + timeByPat[[i]]$waiting[j,k]
            
            if(!is.na(timeByPat[[i]]$k_soj[j,k]))
              if(ceiling(timeByPat[[i]]$k_soj[j,k])<=max(as.numeric(timeFrames)))
                timeByPat[[i]]$sojourn[j,k] <- nodeTime_byPatType[rownames(timeByPat[[i]]$waiting)[j],as.character(ceiling(timeByPat[[i]]$k_soj[j,k])),names(timeByPat)[i]]
              
        }
      }
    }
    
  }
 
  ######### Expected sojourn time
  expectedSojournTime_byArrTime_byPatType <- do.call(cbind.data.frame,lapply(timeByPat,function(x){apply(cbind(t(x$waiting),t(x$sojourn)),1,sum)}))
  rownames(expectedSojournTime_byArrTime_byPatType) <- timeFrames
  expectedSojournTime_byArrTime <- apply(expectedSojournTime_byArrTime_byPatType * arrivals[rownames(expectedSojournTime_byArrTime_byPatType),] / apply(data.frame(arrivals[rownames(expectedSojournTime_byArrTime_byPatType),]),1,sum),1,sum)
  
  temp <- data.frame(na.omit(expectedSojournTime_byArrTime_byPatType[selTimeFrames,]))
  rownames(temp) <- selTimeFrames[selTimeFrames%in%rownames(na.omit(expectedSojournTime_byArrTime_byPatType))]
  expectedSojournTime <- sum(temp * arrivals[rownames(temp),] / sum(arrivals[rownames(temp),]))

  
  ######### Proportion of patients spending more than a target time period in the system
  propWithinTarget_byArrTime_byPatType <- do.call(cbind.data.frame,lapply(timeByPat,function(x){apply(cbind(t(x$waiting),t(x$sojourn))[timeFrames%in%rownames(na.omit(expectedSojournTime_byArrTime_byPatType)),],1,function(y){phypoexp(targetTime,rate=1/y[y>0])})}))
  propWithinTarget_byArrTime_byPatType[propWithinTarget_byArrTime_byPatType==Inf] <- 1
  rownames(propWithinTarget_byArrTime_byPatType) <- rownames(expectedSojournTime_byArrTime_byPatType)[timeFrames%in%rownames(na.omit(expectedSojournTime_byArrTime_byPatType))]
  propWithinTarget_byArrTime <- apply(propWithinTarget_byArrTime_byPatType * arrivals[rownames(propWithinTarget_byArrTime_byPatType),] / apply(data.frame(arrivals[rownames(propWithinTarget_byArrTime_byPatType),]),1,sum),1,sum)
  
  temp <- data.frame(na.omit(propWithinTarget_byArrTime_byPatType[selTimeFrames,]))
  rownames(temp) <- selTimeFrames[selTimeFrames%in%rownames(na.omit(propWithinTarget_byArrTime_byPatType))]
  propWithinTarget <- sum(temp * arrivals[rownames(temp),] / sum(arrivals[rownames(temp),]))
  
  
  ######### Patients in the system
  patInSyst <- round(sum(expectedSojournTime_byArrTime[selTimeFrames] * apply(arrivals[selTimeFrames,],1,sum) * apply(arrivals[selTimeFrames,],1,sum)) / sum(arrivals[selTimeFrames,]),0)
  patInSyst_byHour <- round(expectedSojournTime_byArrTime[selTimeFrames] * apply(arrivals[selTimeFrames,],1,sum),0)
  
  
  
############################### The output measures in the following section depend on specific assumptions, i.e.: ################
############################### - time to initial assessment is the time until a patient is seen at triage #########################
############################### - time to treatment is the time until treatment for a patient starts in the last ED area visited ###
############################### - time to treatment is specifically computed for a subset of ED areas #################################
  
  ######### Expected time to initial assessment
  expectedTimeToInitialAssessment_byArrTime <- do.call(cbind.data.frame,lapply(timeByPat[grep("_w",names(timeByPat))[1]],function(x){
    if("TRIAGE"%in%rownames(x$waiting)) x$waiting["TRIAGE",]
  }))
  rownames(expectedTimeToInitialAssessment_byArrTime) <- timeFrames
  colnames(expectedTimeToInitialAssessment_byArrTime) <- "timeToInitialAssessment"
  
  temp <- data.frame(na.omit(expectedTimeToInitialAssessment_byArrTime[selTimeFrames,]))
  rownames(temp) <- selTimeFrames[selTimeFrames%in%rownames(na.omit(expectedTimeToInitialAssessment_byArrTime))]
  expectedTimeToInitialAssessment <- sum(temp * apply(arrivals[rownames(temp),grep("_w",colnames(arrivals))],1,sum) / sum(arrivals[rownames(temp),grep("_w",colnames(arrivals))]))


  ######### Expected time to treatment
  expectedTimeToTreatment_byArrTime_byPatType <- do.call(cbind.data.frame,lapply(timeByPat,function(x){
    if(nrow(x$sojourn)==1){
      apply(t(x$waiting),1,sum)
    }else if(nrow(x$sojourn)==2){
      apply(cbind(t(x$waiting),t(x$sojourn)[,1]),1,sum)
    }
  }))
  rownames(expectedTimeToTreatment_byArrTime_byPatType) <- timeFrames
  
  expectedTimeToTreatment_byArrTime <- apply(expectedTimeToTreatment_byArrTime_byPatType * arrivals[rownames(expectedTimeToTreatment_byArrTime_byPatType),] / apply(data.frame(arrivals[rownames(expectedTimeToTreatment_byArrTime_byPatType),]),1,sum),1,sum)
  
  temp <- data.frame(na.omit(expectedTimeToTreatment_byArrTime_byPatType[selTimeFrames,]))
  rownames(temp) <- selTimeFrames[selTimeFrames%in%rownames(na.omit(expectedTimeToTreatment_byArrTime_byPatType))]
  expectedTimeToTreatment <- sum(temp * arrivals[rownames(temp),] / sum(arrivals[rownames(temp),]))
  
  
  ######### Expected time to treatment - GP
  expectedTimeToTreatment_byArrTime_byPatType_GP <- do.call(cbind.data.frame,lapply(timeByPat[grep("GP",names(timeByPat))],function(x){
    if(nrow(x$sojourn)==1){
      apply(t(x$waiting),1,sum)
    }else if(nrow(x$sojourn)==2){
      apply(cbind(t(x$waiting),t(x$sojourn)[,1]),1,sum)
    }
  }))
  rownames(expectedTimeToTreatment_byArrTime_byPatType_GP) <- timeFrames
  
  expectedTimeToTreatment_byArrTime_GP <- apply(expectedTimeToTreatment_byArrTime_byPatType_GP * arrivals[rownames(expectedTimeToTreatment_byArrTime_byPatType_GP),colnames(expectedTimeToTreatment_byArrTime_byPatType_GP)] / apply(data.frame(arrivals[rownames(expectedTimeToTreatment_byArrTime_byPatType_GP),colnames(expectedTimeToTreatment_byArrTime_byPatType_GP)]),1,sum),1,sum)
  
  temp <- data.frame(na.omit(expectedTimeToTreatment_byArrTime_byPatType_GP[selTimeFrames,]))
  rownames(temp) <- selTimeFrames[selTimeFrames%in%rownames(na.omit(expectedTimeToTreatment_byArrTime_byPatType_GP))]
  expectedTimeToTreatment_GP <- sum(temp * arrivals[rownames(temp),colnames(temp)] / sum(arrivals[rownames(temp),colnames(temp)]))
  
  
  ######### Expected time to treatment - UTC
  expectedTimeToTreatment_byArrTime_byPatType_UTC <- do.call(cbind.data.frame,lapply(timeByPat[grep("UTC",names(timeByPat))],function(x){
    if(nrow(x$sojourn)==1){
      apply(t(x$waiting),1,sum)
    }else if(nrow(x$sojourn)==2){
      apply(cbind(t(x$waiting),t(x$sojourn)[,1]),1,sum)
    }
  }))
  rownames(expectedTimeToTreatment_byArrTime_byPatType_UTC) <- timeFrames
  
  expectedTimeToTreatment_byArrTime_UTC <- apply(expectedTimeToTreatment_byArrTime_byPatType_UTC * arrivals[rownames(expectedTimeToTreatment_byArrTime_byPatType_UTC),colnames(expectedTimeToTreatment_byArrTime_byPatType_UTC)] / apply(data.frame(arrivals[rownames(expectedTimeToTreatment_byArrTime_byPatType_UTC),colnames(expectedTimeToTreatment_byArrTime_byPatType_UTC)]),1,sum),1,sum)
  
  temp <- data.frame(na.omit(expectedTimeToTreatment_byArrTime_byPatType_UTC[selTimeFrames,]))
  rownames(temp) <- selTimeFrames[selTimeFrames%in%rownames(na.omit(expectedTimeToTreatment_byArrTime_byPatType_UTC))]
  expectedTimeToTreatment_UTC <- sum(temp * arrivals[rownames(temp),colnames(temp)] / sum(arrivals[rownames(temp),colnames(temp)]))
  
  
  
  ######### Expected time to treatment - MAJORS
  expectedTimeToTreatment_byArrTime_byPatType_MAJORS <- do.call(cbind.data.frame,lapply(timeByPat[grep("MAJORS",names(timeByPat))],function(x){
    if(nrow(x$sojourn)==1){
      apply(t(x$waiting),1,sum)
    }else if(nrow(x$sojourn)==2){
      apply(cbind(t(x$waiting),t(x$sojourn)[,1]),1,sum)
    }
  }))
  rownames(expectedTimeToTreatment_byArrTime_byPatType_MAJORS) <- timeFrames
  
  expectedTimeToTreatment_byArrTime_MAJORS <- apply(expectedTimeToTreatment_byArrTime_byPatType_MAJORS * arrivals[rownames(expectedTimeToTreatment_byArrTime_byPatType_MAJORS),colnames(expectedTimeToTreatment_byArrTime_byPatType_MAJORS)] / apply(data.frame(arrivals[rownames(expectedTimeToTreatment_byArrTime_byPatType_MAJORS),colnames(expectedTimeToTreatment_byArrTime_byPatType_MAJORS)]),1,sum),1,sum)
  
  temp <- data.frame(na.omit(expectedTimeToTreatment_byArrTime_byPatType_MAJORS[selTimeFrames,]))
  rownames(temp) <- selTimeFrames[selTimeFrames%in%rownames(na.omit(expectedTimeToTreatment_byArrTime_byPatType_MAJORS))]
  expectedTimeToTreatment_MAJORS <- sum(temp * arrivals[rownames(temp),colnames(temp)] / sum(arrivals[rownames(temp),colnames(temp)]))
  
  
  
  ######### Expected time to treatment - RESUS
  expectedTimeToTreatment_byArrTime_byPatType_RESUS <- do.call(cbind.data.frame,lapply(timeByPat[grep("RESUS",names(timeByPat))],function(x){
    if(nrow(x$sojourn)==1){
      apply(t(x$waiting),1,sum)
    }else if(nrow(x$sojourn)==2){
      apply(cbind(t(x$waiting),t(x$sojourn)[,1]),1,sum)
    }
  }))
  rownames(expectedTimeToTreatment_byArrTime_byPatType_RESUS) <- timeFrames
  
  expectedTimeToTreatment_byArrTime_RESUS <- apply(expectedTimeToTreatment_byArrTime_byPatType_RESUS * arrivals[rownames(expectedTimeToTreatment_byArrTime_byPatType_RESUS),colnames(expectedTimeToTreatment_byArrTime_byPatType_RESUS)] / apply(data.frame(arrivals[rownames(expectedTimeToTreatment_byArrTime_byPatType_RESUS),colnames(expectedTimeToTreatment_byArrTime_byPatType_RESUS)]),1,sum),1,sum)
  
  temp <- data.frame(na.omit(expectedTimeToTreatment_byArrTime_byPatType_RESUS[selTimeFrames,]))
  rownames(temp) <- selTimeFrames[selTimeFrames%in%rownames(na.omit(expectedTimeToTreatment_byArrTime_byPatType_RESUS))]
  expectedTimeToTreatment_RESUS <- sum(temp * arrivals[rownames(temp),colnames(temp)] / sum(arrivals[rownames(temp),colnames(temp)]))
  
#########################################################################################################################################
#########################################################################################################################################
#########################################################################################################################################
  
  

  return(list(
    propWithinTarget=propWithinTarget,
    expectedSojournTime=expectedSojournTime,
    propWithinTarget_byArrTime=propWithinTarget_byArrTime[selTimeFrames],
    expectedSojournTime_byArrTime=expectedSojournTime_byArrTime[selTimeFrames],
    propWithinTarget_byArrTime_byPatType=propWithinTarget_byArrTime_byPatType[selTimeFrames,],
    expectedSojournTime_byArrTime_byPatType=expectedSojournTime_byArrTime_byPatType[selTimeFrames,],
    expectedTimeToInitialAssessment=expectedTimeToInitialAssessment,
    expectedTimeToInitialAssessment_byArrTime=expectedTimeToInitialAssessment_byArrTime[selTimeFrames,],
    expectedTimeToTreatment=expectedTimeToTreatment,
    expectedTimeToTreatment_GP=expectedTimeToTreatment_GP,
    expectedTimeToTreatment_UTC=expectedTimeToTreatment_UTC,
    expectedTimeToTreatment_MAJORS=expectedTimeToTreatment_MAJORS,
    expectedTimeToTreatment_RESUS=expectedTimeToTreatment_RESUS,
    expectedTimeToTreatment_byArrTime_byPatType=expectedTimeToTreatment_byArrTime_byPatType[selTimeFrames,],
    expectedTimeToTreatment_byArrTime=expectedTimeToTreatment_byArrTime[selTimeFrames],
    expectedTimeToTreatment_byArrTime_GP=expectedTimeToTreatment_byArrTime_GP[selTimeFrames],
    expectedTimeToTreatment_byArrTime_UTC=expectedTimeToTreatment_byArrTime_UTC[selTimeFrames],
    expectedTimeToTreatment_byArrTime_MAJORS=expectedTimeToTreatment_byArrTime_MAJORS[selTimeFrames],
    expectedTimeToTreatment_byArrTime_RESUS=expectedTimeToTreatment_byArrTime_RESUS[selTimeFrames],
    patInSyst = patInSyst,
    patInSyst_byHour = patInSyst_byHour
  ))
  
  
  
  
  
  
}





	










library(sdprisk)

#load functions
source("functions/erlangLoss.R")
source("functions/queueSize.R")
source("functions/analyticalModel.R")


#we set the reference time frame to 1 (hour) and define 72 equidistant time frames (3 days)
dt <- 1
timeFrames <- as.character(seq(1,72,dt))

#each row is an ED area; each column is a time frame; capacity represents the number of patients that can be in that area in a given time frame
capacities <- as.matrix(t(read.table("parameters/capacities.txt",header=TRUE,row.names=1,dec=".",sep="\t")))
capacities <- matrix(rep(capacities,3),nrow=nrow(capacities),byrow=FALSE,dimnames=list(rownames(capacities),timeFrames))

#each row is an ED area; each column is a patient type (any name allowed); 1 if a given patient type goes through a given area, 0 otherwise
patientTypes <- as.matrix(read.table("parameters/patientTypes.txt",header=TRUE,row.names=1,dec=".",sep="\t"))

#each row is a time frame; each column is a patient type; the number of patients per time frame (hour) arriving at the ED
arrivals <- as.matrix(read.table("parameters/arrivals.txt",header=TRUE,row.names=1,dec=".",sep="\t"))
arrivals <- matrix(rep(t(arrivals),3),ncol=ncol(arrivals),byrow=TRUE,dimnames=list(timeFrames,colnames(arrivals)))

#list of data frames (one for each ED area); each data frame reports the treatment time for a given patient type in that area in each time frame
#any additional ED area should be added manually below
nodeTimes <- lapply(rownames(capacities),function(n){
  tmp <- as.matrix(read.table(paste("parameters/nodeTimes_",n,".txt",sep=''),header=TRUE,row.names=1,dec=".",sep="\t"))
  matrix(rep(t(tmp),3),ncol=ncol(tmp),byrow=TRUE,dimnames=list(timeFrames,colnames(tmp)))
})
names(nodeTimes) <- rownames(capacities)

#we choose to consider system performance from the 25th to 48th hour
selTimeFrames <- timeFrames[floor(length(timeFrames)/3+1):floor(length(timeFrames)*2/3)]

#output will be the proportion of patients leaving the ED within this target time (i.e. 4-hour target)
targetTime <- 4

#### 1. Running the model once, with above parameters
modRes_basic <- analyticalModel(capacities,arrivals,patientTypes,nodeTimes,timeFrames,selTimeFrames,dt,targetTime)

#### 2. Running the model for a set of scenarios

#modifying arrival rates
scenArrivals <- read.table("scenarios/scenArrivals.txt",header=TRUE,row.names=1,dec=".",sep="\t")
modRes_scenArrivals <- vector("list",length=nrow(scenArrivals)+1)
names(modRes_scenArrivals) <- c("basic",rownames(scenArrivals))

modRes_scenArrivals$basic <- analyticalModel(capacities,arrivals,patientTypes,nodeTimes,timeFrames,selTimeFrames,dt,targetTime)
for(i in 2:length(modRes_scenArrivals)){
  tempArr <- arrivals %*% diag(1+scenArrivals[i-1,])
  rownames(tempArr) <- rownames(arrivals)
  colnames(tempArr) <- colnames(arrivals)
  modRes_scenArrivals[[i]] <- analyticalModel(capacities,tempArr,patientTypes,nodeTimes,timeFrames,selTimeFrames,dt,targetTime)
  rm(tempArr)
}

#modifying capacities
scenCapacities <- read.table("scenarios/scenCapacities.txt",header=TRUE,row.names=1,dec=".",sep="\t")
modRes_scenCapacities <- vector("list",length=nrow(scenCapacities)+1)
names(modRes_scenCapacities) <- c("basic",rownames(scenCapacities))

modRes_scenCapacities$basic <- analyticalModel(capacities,arrivals,patientTypes,nodeTimes,timeFrames,selTimeFrames,dt,targetTime)
for(i in 2:length(modRes_scenCapacities)){
  tempCap <- capacities
  for(j in 1:nrow(capacities)){
    tempCap[j,] <- rep(scenCapacities[i-1,rownames(capacities)[j]],ncol(capacities))
  }
  modRes_scenCapacities[[i]] <- analyticalModel(tempCap,arrivals,patientTypes,nodeTimes,timeFrames,selTimeFrames,dt,targetTime)
  rm(tempCap)
}

#modifying node times
scenNodeTimes <- read.table("scenarios/scenNodeTimes.txt",header=TRUE,row.names=1,dec=".",sep="\t")
modRes_scenNodeTimes <- vector("list",length=nrow(scenNodeTimes)+1)
names(modRes_scenNodeTimes) <- c("basic",rownames(scenNodeTimes))

modRes_scenNodeTimes$basic <- analyticalModel(capacities,arrivals,patientTypes,nodeTimes,timeFrames,selTimeFrames,dt,targetTime)
for(i in 2:length(modRes_scenNodeTimes)){
  tempNodeTimes <- nodeTimes
  for(j in 1:length(tempNodeTimes)){
    tempNodeTimes[[j]] <- tempNodeTimes[[j]] * (1+scenNodeTimes[i-1,names(tempNodeTimes)[j]])
  }
  modRes_scenNodeTimes[[i]] <- analyticalModel(capacities,arrivals,patientTypes,tempNodeTimes,timeFrames,selTimeFrames,dt,targetTime)
  rm(tempNodeTimes)
}


### 3. Result summaries

#performance (4-hour target)
sapply(modRes_scenArrivals,function(x)x$propWithinTarget)
sapply(modRes_scenCapacities,function(x)x$propWithinTarget)
sapply(modRes_scenNodeTimes,function(x)x$propWithinTarget)

#expected sojourn time
sapply(modRes_scenArrivals,function(x)x$expectedSojournTime)
sapply(modRes_scenCapacities,function(x)x$expectedSojournTime)
sapply(modRes_scenNodeTimes,function(x)x$expectedSojournTime)

#expected time to treatment
sapply(modRes_scenArrivals,function(x)x$expectedTimeToTreatment)
sapply(modRes_scenCapacities,function(x)x$expectedTimeToTreatment)
sapply(modRes_scenNodeTimes,function(x)x$expectedTimeToTreatment)

#expected number of patients in the system
sapply(modRes_scenArrivals,function(x)x$patInSyst)
sapply(modRes_scenCapacities,function(x)x$patInSyst)
sapply(modRes_scenNodeTimes,function(x)x$patInSyst)

#expected number of patients in the system by hour
sapply(modRes_scenArrivals,function(x)x$patInSyst_byHour)
sapply(modRes_scenCapacities,function(x)x$patInSyst_byHour)
sapply(modRes_scenNodeTimes,function(x)x$patInSyst_byHour)


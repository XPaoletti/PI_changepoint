##################################################################
###          FUNCTIONS TO USE IN MAKEDATA AND LOGLIK           ###
Sigma.is <- function(upars){  
  sig0 <- (upars[1])^2
}

L.is <- function(p=nt,gampars){
  gam0 <- gampars[1]
  mat <- matrix(0,nrow=p,ncol=1)
  mat[,1] <- gam0
  mat
}


#############################################################
#####                DATA GENERATION                    #####

makedata<- function(meanpars,varpars, gampars, nt, n1, dose, d, nvisit, pl){
  # Generates both survival and longitudinal data before starting the trial
  # For each subject it generates all possible outcomes for all doses and all treatment cycles, for the survival outcome, and all possible outcomes for all doses and all time visits for the longitudinal data 
  
  
  by <- meanpars[1:3]  # parameters for the longitudinal outcome
  bt<- meanpars[4:6]   # for the probit outcome
  sig.y <- varpars[1]  # residual variance
  upars <- varpars[-1] # random effect variance 
  
  
  Sig.u <- Sigma.is(upars) 
  L <- L.is(nt,gampars) 
  
  
  # Simulates from multivariate normal distribution the random effects.
  u <- rmnorm(n1, mean=rep(0,1), varcov=Sig.u) 
  Lu <- L%*%t(u)  
  Lu <- t(Lu)
  
  
  # Data generation for the probit model
  surv<- p <- array(0, dim=c(d,nt,n1)) 
  for (i in 1:d){
    for (j in 1:nt){
      p[i,j,] <- pnorm(bt[1] + bt[2]*(j-1)+ bt[3]*dose[i] + Lu[,j])  
      surv[i,j,] <- rbinom(n=n1, size=1, prob=p[i,j,])  
      
    }
  }
  
  # Data generation for the longitudinal model
  # If pl=1 then there is no plateau and activity is decreasing or increasing with dose
  # if pl>1 then there is a plateau at the specific dose level
  # After the plateau beginning, all doses are assumed to be the same
  visit <- c(0.1, 0.4, 0.8, 1.2, 1.4, 1.8, 2.1, 2.4, 2.8, 3.1, 3.4, 3.8, 4.1, 4.4, 4.8, 5.1, 5.4, 5.8)
  
  
  if(pl==1){
    y<- mu.y <- array(0, dim=c(d,nvisit,n1))
    u<- as.array(u, dim=c(1,1,n1))
    for (i in 1:d){
      for (j in 1:nvisit){
        mu.y[i,j,] <- (by[1] + by[2]*(visit[j])^2 + by[3]*visit[j]*dose[i] + u*visit[j])
        y[i,j,] <- rnorm(length(mu.y[i,j,]), mu.y[i,j,], sd=sig.y)
      }
    }
  }else{
    dose1<- ifelse((dose<dose[pl]),dose, dose[pl] ) 
    
    y<- mu.y <- array(0, dim=c(d,nvisit,n1))
    u<- as.array(u, dim=c(1,1,n1))
    for (i in 1:d){
      for (j in 1:nvisit){
        mu.y[i,j,] <- (by[1] + by[2]*(visit[j])^2 + by[3]*visit[j]*dose1[i]  + u*visit[j]) 
        y[i,j,] <- rnorm(length(mu.y[i,j,]), mu.y[i,j,], sd=sig.y)
      }
    }
  }  
  
  out<- list(surv=surv, rawlong=y)  
  
}




#############################################################
#####       SURV AND LONG DATA 2+2 BEGINNING            #####

databegin <- function(met, surv, nvisit, dose, inter_miss){
  # Extracts the longitudinal data from the initial simulated data after the first cohort of subjects has been enrolled
  # Longitudinal data are selected from the simulated data based on the dose administered to each subject 
  # Censoring is applied to the longitudinal data
  # Second censoring due to lack of activity is applied to the longitudinal data
  # Finally the same censoring is applied to the survival data  
  
  
  met1 <- surv  # survival data after 2+2
  nsurv <- nrow(met1)
  met2 <- met$rawlong  # simulated longitudinal data  
  
  
  visittime <- matrix(NA, nrow = (nvisit*nsurv), ncol=1)
  for(i in 1:nsurv){
    x <- rep(NA, nvisit)
    if(met1$time[i]==1){
      x <- c(c(0.1, 0.4, 0.8), rep(NA,15))
    } else if(met1$time[i]==2){
      x <- c(c(0.1, 0.4, 0.8, 1.1, 1.4, 1.8), rep(NA, 12))
    } else if(met1$time[i]==3){
      x <- c(c(0.1, 0.4, 0.8, 1.1, 1.4, 1.8, 2.1, 2.4, 2.8), rep(NA, 9))
    } else if(met1$time[i]==4){
      x <- c(c(0.1, 0.4, 0.8, 1.1, 1.4, 1.8, 2.1, 2.4, 2.8, 3.1, 3.4, 3.8), rep(NA, 6))
    } else if(met1$time[i]==5){
      x <- c(c(0.1, 0.4, 0.8, 1.1, 1.4, 1.8, 2.1, 2.4, 2.8, 3.1, 3.4, 3.8, 4.1, 4.4, 4.8), rep(NA, 3))
    }else{
      x <- c(0.1, 0.4, 0.8, 1.1, 1.4, 1.8, 2.1, 2.4, 2.8, 3.1, 3.4, 3.8, 4.1, 4.4, 4.8, 5.1, 5.4, 5.8)
    }
    x <- as.matrix(x)
    visittime <- rbind(visittime, x)
  }  
  visittime <- visittime[(nvisit*nsurv+1):(nvisit*nsurv+nvisit*nsurv),1] 
  reps <- rep(nvisit, nsurv)
  id <- rep(met1$id, reps) 
  tempdose <- rep(met1$dose, reps) 
  
  y <- matrix(NA, nrow = (nvisit*nsurv), ncol=1)  
  for(i in 1:nsurv){
    x <- rep(NA, nvisit)
    indexdose <- which(dose==met1$dose[i])
    x <- met2[indexdose,,i]
    x <- as.matrix(x)
    y <- rbind(y, x)
    
  }
  y <- y[(nvisit*nsurv+1):(nvisit*nsurv+nvisit*nsurv),1]  
  
  tempdat <- data.frame(id, y, visittime, dose=tempdose)
  tempdat <- na.omit(tempdat)  
  
  
  
  listid <- split.data.frame(tempdat, as.factor(tempdat$id))
  
  findmax <- function(y){
    
    chose <- c(0, rbinom((nrow(y)-1), 1, inter_miss)) 
    for(i in 1:length(chose)){
      if(chose[i]==1){
        y[i,2] <- NA
      }
    }
    y <- na.omit(y)
    
    # Censoring due to lack of efficacy
    upplim <- 35
    maxeff <- min(y[2]) 
    pos_max=which(y[2]==maxeff)  
    
    if(pos_max!=dim(y)[1]){ 
      
      var_a_test=c(rep(NA, pos_max), seq(from=pos_max+1,to=dim(y)[1]))  
      if(maxeff<= upplim){
        pos=which(y[var_a_test,2]>= 2*upplim) 
      }else{
        pos=which(y[var_a_test,2]>= 2*maxeff) 
      }
      round <- ceiling(y[pos[1],3]) 
      if(any(which(y[,3]>=round))){
        sup <- which(y[,3]>=round) 
        y[sup,2]=NA 
        y <- na.omit(y)
      }
    }    
    
    maxtime <- ceiling(max(y[,3]))
    out <- list(y=y, maxtime=maxtime)    
  }
  
  
  res <- lapply(listid, findmax)
  
  maxvisit <- rep(NA, nsurv)
  for(i in 1:nsurv){
    maxvisit[i] <- res[[i]][2]  
  }  
  maxvisit <- unlist(maxvisit)
  maxvisit <- as.matrix(maxvisit)
  
  for(i in 1:nsurv){
    if((met1$time[i]>maxvisit[i]) & (met1$status[i]==0)){
      met1$time[i] <- maxvisit[i] 
      met1$status[i] <- 1
    } else{
      met1$time[i] <- maxvisit[i] 
    }
  }  
  
  
  dat1 <- data.frame(y.id=NA, y.y=NA, y.visittime=NA, y.dose=NA)  
  tempdat <-  data.frame(res[[1]][1][1])
  dat1 <- rbind(dat1,tempdat)  
  for (i in 2:nsurv) {
    tempdat <-  data.frame(res[[i]][1][1])
    dat1 <- rbind(dat1[],tempdat)
  }
  
  dat1 <- na.omit(dat1)
  names(dat1)[names(dat1)=="y.id"] <- "id"
  names(dat1)[names(dat1)=="y.y"] <- "y"
  names(dat1)[names(dat1)=="y.visittime"] <- "visittime"
  names(dat1)[names(dat1)=="y.dose"] <- "dose"
  
  out <- list(dat=dat1, surv=met1, rawlong=met2)
}





#############################################################
#####            UP AND DOWN SCHEME FOR 2+2             #####

updown <- function(dose, i , y, test2, d){
  # Applies the 2+2 escalation scheme
  # If no DLT is observed add 2 more subjects to the next dose level
  # If 1 DLT is observed add 2 more subjects to the same dose level
  # If 2 DLTs are observed add 2 more subjects to the previous dose level
  # The escalation procedure depends on the first treatment cycle only  
  # In case of 2 DLTs in the first cycle, the next 2 subjects will be assigned to the same (first) dose level 
  # In of case no DLTs at the last dose level the next 2 subjects will be assigned to the same (last) dose level 
  doselevel<- which(dose==y[i-1])
  
  if(doselevel==1 & (test2[doselevel,1,i-2]==(0)) & (test2[doselevel,1,i-1]==(0))){ 
    (y[i]<-dose[1]) & (y[i+1]<-dose[1]) 
  }else if (doselevel==d & (test2[d,1,i-2]==(1)) & (test2[d,1,i-1]==(1))){ 
    (y[i]<-dose[d]) & (y[i+1]<-dose[d]) 
  } else if ((test2[doselevel,1,i-2]==(1)) & (test2[doselevel,1,i-1]==(1))){ 
    (y[i]<- dose[doselevel+1]) & (y[i+1]<- dose[doselevel+1])
  }else if ((test2[doselevel,1,i-2]==(0)) & (test2[doselevel,1,i-1]==(1)) | 
            (test2[doselevel,1,i-1]==(0)) & (test2[doselevel,1,i-2]==(1))){ 
    (y[i]<- dose[doselevel]) & (y[i+1]<- dose[doselevel]) 
  }else if ((test2[doselevel,1,i-2]==(0)) & (test2[doselevel,1,i-1]==(0))) { 
    (y[i]<-dose[doselevel-1]) & (y[i+1]<- dose[doselevel-1])
  }
  out<- y
} 



#############################################################
#####     SURV AND LONG DATA DURING 2+2 CRM AND JM      #####

data2plcrmjm <- function(cool, survdata, longdata, nvisit, dose, inter_miss){
  # Same function as databegin for all subjects entering the trial during the 2+2, CRM and JM
  
  met1 <- survdata
  nsurv <- nrow(met1)
  dat1 <- longdata  
  met2 <- cool$rawlong  # simulated longitudinal
  
  x <- rep(NA, nvisit)
  i<- nsurv
  if(met1$time[i]==1){
    x <- c(c(0.1, 0.4, 0.8), rep(NA,15))
  } else if(met1$time[i]==2){
    x <- c(c(0.1, 0.4, 0.8, 1.1, 1.4, 1.8), rep(NA, 12))
  } else if(met1$time[i]==3){
    x <- c(c(0.1, 0.4, 0.8, 1.1, 1.4, 1.8, 2.1, 2.4, 2.8), rep(NA, 9))
  } else if(met1$time[i]==4){
    x <- c(c(0.1, 0.4, 0.8, 1.1, 1.4, 1.8, 2.1, 2.4, 2.8, 3.1, 3.4, 3.8), rep(NA, 6))
  } else if(met1$time[i]==5){
    x <- c(c(0.1, 0.4, 0.8, 1.1, 1.4, 1.8, 2.1, 2.4, 2.8, 3.1, 3.4, 3.8, 4.1, 4.4, 4.8), rep(NA, 3))
  }else{
    x <- c(0.1, 0.4, 0.8, 1.1, 1.4, 1.8, 2.1, 2.4, 2.8, 3.1, 3.4, 3.8, 4.1, 4.4, 4.8, 5.1, 5.4, 5.8)
  }
  visittime <- as.matrix(x)
  
  
  id <- as.matrix(rep(nsurv, nvisit)) 
  tempdose <- as.matrix(rep(met1$dose[nsurv], nvisit)) 
  indexdose <- which(dose==met1$dose[nsurv])
  y <- as.matrix(met2[indexdose,,nsurv])
  
  
  
  tempdat <- data.frame(id, y, visittime, dose=tempdose)
  tempdat <- na.omit(tempdat)  
  
  # Intermittent censoring is applied to all visits apart from the first one, so as to have at least one measurement per subject
  chose <- c(0, rbinom((nrow(tempdat)-1), 1, inter_miss))
  for(i in 1:length(chose)){
    if(chose[i]==1){
      tempdat[i,2] <- NA
    }
  }
  tempdat <- na.omit(tempdat)
  upplim <- 35
  maxeff <- min(tempdat$y) 
  pos_max=which(tempdat$y==maxeff) 
  
  if(pos_max!=dim(tempdat)[1]){ 
    var_a_test=c(rep(NA, pos_max), seq(from=pos_max+1,to=dim(tempdat)[1]))  
    if(maxeff<= upplim){
      pos=which(tempdat[var_a_test,2]>= 2*upplim) 
    }else{
      pos=which(tempdat[var_a_test,2]>= 2*maxeff) 
    }
    
    round <- ceiling(tempdat[pos[1],3]) 
    if(any(which(tempdat[,3]>=round))){
      sup <- which(tempdat[,3]>=round) 
      tempdat[sup,2]=NA 
      tempdat <- na.omit(tempdat)
    }
  } 
  maxtime <- ceiling(max(tempdat[,3]))
  
  if((met1$time[nsurv]>maxtime) & (met1$status[nsurv]==0)){
    met1$time[nsurv] <- maxtime 
    met1$status[nsurv] <- 1
  } else{
    met1$time[nsurv] <- maxtime     
  }  
  dat1 <- rbind(dat1,tempdat)
  
  out <- list(dat=dat1, surv=met1, long=met2)
}



#######################################################################
#####      SEQUENTIAL ENTRY PER CYCLE FOR SUBJECTS IN 2+2         #####

sequentialentry0 <- function(test, test1, i){
  # Analyzes subject data during the 2+2 design
  # The first time a cohort is analyzed it takes the data from their first treatment cycle 
  # When a new cohort is enrolled for the previous cohort we use the data till the second treatment cycle (if there is one), whereas for the new cohort we use the data from their first cycle
  
  vb <- vbt <- vbc <- vbk <- vbl <-  vb1 <- vbt1 <- vbc1 <- vbk1 <- vbl1 <- NULL
  mb <- mbt <- mbc <- mbk <- mbl <-  mb1 <- mbt1 <- mbc1 <- mbk1 <- mbl1 <- NULL
  
  
  if(isTRUE((test$id[i]==test$counter[i]))){ 
    if((test$status[i]==0 & test$time[i]==3) | (test$status[i]==0 & test$time[i]==2) |(test$status[i]==1 & test$time[i]==3) |
       (test$status[i]==1 & test$time[i]==2) | (test$status[i]==1 & test$time[i]==1) | (test$status[i]==0 & test$time[i]==4) |
       (test$status[i]==0 & test$time[i]==5) | (test$status[i]==0 & test$time[i]==6)| (test$status[i]==1 & test$time[i]==4) | 
       (test$status[i]==1 & test$time[i]==5) | (test$status[i]==1 & test$time[i]==6) ){
      vb<- c(i,i,1,1,test$dose[i])
      mb <- test1[which(test1$id==i & test1$visittime<1),]
      
    }else{
      vb <- c(i,i,0,1,test$dose[i])
      mb <- test1[which(test1$id==i & test1$visittime<1),]
    }
  }
  
  if(isTRUE((test$id[i+1]-test$counter[i])==1)){ 
    if((test$status[i+1]==0 & test$time[i+1]==3) | (test$status[i+1]==0 & test$time[i+1]==2) |(test$status[i+1]==1 & test$time[i+1]==3) |
       (test$status[i+1]==1 & test$time[i+1]==2) | (test$status[i+1]==1 & test$time[i+1]==1) | (test$status[i+1]==0 & test$time[i+1]==4) |
       (test$status[i+1]==0 & test$time[i+1]==5) | (test$status[i+1]==0 & test$time[i+1]==6)| (test$status[i+1]==1 & test$time[i+1]==4) | 
       (test$status[i+1]==1 & test$time[i+1]==5) | (test$status[i+1]==1 & test$time[i+1]==6) ){
      vb1<- c(i+1,i+1,1,1,test$dose[i+1])
      mb1 <- test1[which(test1$id==(i+1) & test1$visittime<1),]
    }else{
      vb1 <- c(i+1,i+1,0,1,test$dose[i+1])
      mb1 <- test1[which(test1$id==(i+1) & test1$visittime<1),]
    }
  }
  
  
  if(isTRUE((test$counter[i]-test$id[i-1])==1)){
    if((test$status[i-1]==0 & test$time[i-1]==3) | (test$status[i-1]==0 & test$time[i-1]==4) |
       (test$status[i-1]==0 & test$time[i-1]==5) | (test$status[i-1]==0 & test$time[i-1]==6) ){
      vbt<- c(i-1,i-1,1,2,test$dose[i-1])
      mbt <- test1[which(test1$id==(i-1) & test1$visittime<2),]
    }else if ((test$status[i-1]==0 & test$time[i-1]==2)){
      vbt <- c(i-1,i-1,0,2,test$dose[i-1])
      mbt <- test1[which(test1$id==(i-1) & test1$visittime<2),]
    }else if((test$status[i-1]==0 & test$time[i-1]==1)){
      vbt <- c(i-1,i-1,0,1,test$dose[i-1])
      mbt <- test1[which(test1$id==(i-1) & test1$visittime<1),]
    }else if((test$status[i-1]==1 & test$time[i-1]==2) | (test$status[i-1]==1 & test$time[i-1]==3) |
             (test$status[i-1]==1 & test$time[i-1]==4) | (test$status[i-1]==1 & test$time[i-1]==5) |
             (test$status[i-1]==1 & test$time[i-1]==6)){
      vbt <- c(i-1,i-1,1,2,test$dose[i-1])
      mbt <- test1[which(test1$id==(i-1) & test1$visittime<2),]
    }else{
      vbt <- c(i-1,i-1,1,1,test$dose[i-1])
      mbt <- test1[which(test1$id==(i-1) & test1$visittime<1),]
    }
  }
  
  if(isTRUE((test$counter[i]-test$id[i-2])==2)){
    if((test$status[i-2]==0 & test$time[i-2]==3) | (test$status[i-2]==0 & test$time[i-2]==4) |
       (test$status[i-2]==0 & test$time[i-2]==5) | (test$status[i-2]==0 & test$time[i-2]==6) ){
      vbt1 <- c(i-2,i-2,1,2,test$dose[i-2])
      mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<2),]
    }else if ((test$status[i-2]==0 & test$time[i-2]==2)){
      vbt1 <- c(i-2,i-2,0,2,test$dose[i-2])
      mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<2),]
    }else if((test$status[i-2]==0 & test$time[i-2]==1)){
      vbt1 <- c(i-2,i-2,0,1,test$dose[i-2])
      mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<1),]
    }else if((test$status[i-2]==1 & test$time[i-2]==2) | (test$status[i-2]==1 & test$time[i-2]==3) |
             (test$status[i-2]==1 & test$time[i-2]==4) | (test$status[i-2]==1 & test$time[i-2]==5) |
             (test$status[i-2]==1 & test$time[i-2]==6)){
      vbt1 <- c(i-2,i-2,1,2,test$dose[i-2])
      mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<2),]
    }else{
      vbt1 <- c(i-2,i-2,1,1,test$dose[i-2])
      mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<1),]
    }
  }   
  
  if(isTRUE((test$counter[i]-test$id[i-3])==3)){
    if((test$status[i-3]==0 & test$time[i-3]==3) ){
      vbc<- c(i-3,i-3,0,3,test$dose[i-3])
      mbc <- test1[which(test1$id==(i-3) & test1$visittime<3),]
    }else if ((test$status[i-3]==0 & test$time[i-3]==2)){
      vbc <- c(i-3,i-3,0,2,test$dose[i-3])
      mbc <- test1[which(test1$id==(i-3) & test1$visittime<2),]
    }else if((test$status[i-3]==0 & test$time[i-3]==1)){
      vbc <- c(i-3,i-3,0,1,test$dose[i-3])
      mbc <- test1[which(test1$id==(i-3) & test1$visittime<1),]
    }else if((test$status[i-3]==0 & test$time[i-3]==4) | (test$status[i-3]==0 & test$time[i-3]==5) |
             (test$status[i-3]==0 & test$time[i-3]==6) | (test$status[i-3]==1 & test$time[i-3]==3) |
             (test$status[i-3]==1 & test$time[i-3]==4) | (test$status[i-3]==1 & test$time[i-3]==5) |
             (test$status[i-3]==1 & test$time[i-3]==6)){
      vbc <- c(i-3,i-3,1,3,test$dose[i-3])
      mbc <- test1[which(test1$id==(i-3) & test1$visittime<3),]
    }else if((test$status[i-3]==1 & test$time[i-3]==2)){
      vbc <- c(i-3,i-3,1,2,test$dose[i-3])
      mbc <- test1[which(test1$id==(i-3) & test1$visittime<2),]
    }else{
      vbc <- c(i-3,i-3,1,1,test$dose[i-3])
      mbc <- test1[which(test1$id==(i-3) & test1$visittime<1),]
    }
  }
  
  if(isTRUE((test$counter[i]-test$id[i-4])==4)){
    if((test$status[i-4]==0 & test$time[i-4]==3) ){
      vbc1<- c(i-4,i-4,0,3,test$dose[i-4])
      mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<3),]
    }else if ((test$status[i-4]==0 & test$time[i-4]==2)){
      vbc1 <- c(i-4,i-4,0,2,test$dose[i-4])
      mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<2),]
    }else if((test$status[i-4]==0 & test$time[i-4]==1)){
      vbc1 <- c(i-4,i-4,0,1,test$dose[i-4])
      mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<1),]
    }else if((test$status[i-4]==0 & test$time[i-4]==4) | (test$status[i-4]==0 & test$time[i-4]==5) |
             (test$status[i-4]==0 & test$time[i-4]==6) | (test$status[i-4]==1 & test$time[i-4]==3) |
             (test$status[i-4]==1 & test$time[i-4]==4) | (test$status[i-4]==1 & test$time[i-4]==5) |
             (test$status[i-4]==1 & test$time[i-4]==6)){
      vbc1 <- c(i-4,i-4,1,3,test$dose[i-4])
      mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<3),]
    }else if((test$status[i-4]==1 & test$time[i-4]==2)){
      vbc1 <- c(i-4,i-4,1,2,test$dose[i-4])
      mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<2),]
    }else{
      vbc1 <- c(i-4,i-4,1,1,test$dose[i-4])
      mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<1),]
    }
  }
  
  if(isTRUE((test$counter[i]-test$id[i-5])==5)){
    if((test$status[i-5]==0 & test$time[i-5]==4) ){
      vbk <- c(i-5,i-5,0,4,test$dose[i-5])
      mbk <- test1[which(test1$id==(i-5) & test1$visittime<4),]
    }else if((test$status[i-5]==0 & test$time[i-5]==3) ){
      vbk <- c(i-5,i-5,0,3,test$dose[i-5])
      mbk <- test1[which(test1$id==(i-5) & test1$visittime<3),]
    }else if ((test$status[i-5]==0 & test$time[i-5]==2)){
      vbk <- c(i-5,i-5,0,2,test$dose[i-5])
      mbk <- test1[which(test1$id==(i-5) & test1$visittime<2),]
    }else if((test$status[i-5]==0 & test$time[i-5]==1)){
      vbk <- c(i-5,i-5,0,1,test$dose[i-5])
      mbk <- test1[which(test1$id==(i-5) & test1$visittime<1),]
    }else if((test$status[i-5]==0 & test$time[i-5]==5) |
             (test$status[i-5]==0 & test$time[i-5]==6) |
             (test$status[i-5]==1 & test$time[i-5]==4) | (test$status[i-5]==1 & test$time[i-5]==5) |
             (test$status[i-5]==1 & test$time[i-5]==6)){
      vbk <- c(i-5,i-5,1,4,test$dose[i-5])
      mbk <- test1[which(test1$id==(i-5) & test1$visittime<4),]
    }else if((test$status[i-5]==1 & test$time[i-5]==2)){
      vbk <- c(i-5,i-5,1,2,test$dose[i-5])
      mbk <- test1[which(test1$id==(i-5) & test1$visittime<2),]
    }else if((test$status[i-5]==1 & test$time[i-5]==3)){
      vbk <- c(i-5,i-5,1,3,test$dose[i-5])
      mbk <- test1[which(test1$id==(i-5) & test1$visittime<3),]
    }else{
      vbk <- c(i-5,i-5,1,1,test$dose[i-5])
      mbk <- test1[which(test1$id==(i-5) & test1$visittime<1),]
    }
  }
  
  if(isTRUE((test$counter[i]-test$id[i-6])==6)){
    if((test$status[i-6]==0 & test$time[i-6]==4) ){
      vbk1 <- c(i-6,i-6,0,4,test$dose[i-6])
      mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<4),]
    }else if((test$status[i-6]==0 & test$time[i-6]==3) ){
      vbk1 <- c(i-6,i-6,0,3,test$dose[i-6])
      mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<3),]
    }else if ((test$status[i-6]==0 & test$time[i-6]==2)){
      vbk1 <- c(i-6,i-6,0,2,test$dose[i-6])
      mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<2),]
    }else if((test$status[i-6]==0 & test$time[i-6]==1)){
      vbk1 <- c(i-6,i-6,0,1,test$dose[i-6])
      mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<1),]
    }else if((test$status[i-6]==0 & test$time[i-6]==5) |
             (test$status[i-6]==0 & test$time[i-6]==6) |
             (test$status[i-6]==1 & test$time[i-6]==4) | (test$status[i-6]==1 & test$time[i-6]==5) |
             (test$status[i-6]==1 & test$time[i-6]==6)){
      vbk1 <- c(i-6,i-6,1,4,test$dose[i-6])
      mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<4),]
    }else if((test$status[i-6]==1 & test$time[i-6]==2)){
      vbk1 <- c(i-6,i-6,1,2,test$dose[i-6])
      mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<2),]
    }else if((test$status[i-6]==1 & test$time[i-6]==3)){
      vbk1 <- c(i-6,i-6,1,3,test$dose[i-6])
      mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<3),]
    }else{
      vbk1 <- c(i-6,i-6,1,1,test$dose[i-6])
      mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<1),]
    }
  }
  
  if(isTRUE((test$counter[i]-test$id[i-7])==7)){
    if((test$status[i-7]==0 & test$time[i-7]==5) ){
      vbl <- c(i-7,i-7,0,5,test$dose[i-7])
      mbl <- test1[which(test1$id==(i-7) & test1$visittime<5),]
    }else if((test$status[i-7]==0 & test$time[i-7]==4) ){
      vbl <- c(i-7,i-7,0,4,test$dose[i-7])
      mbl <- test1[which(test1$id==(i-7) & test1$visittime<4),]
    }else if((test$status[i-7]==0 & test$time[i-7]==3) ){
      vbl <- c(i-7,i-7,0,3,test$dose[i-7])
      mbl <- test1[which(test1$id==(i-7) & test1$visittime<3),]
    }else if ((test$status[i-7]==0 & test$time[i-7]==2)){
      vbl <- c(i-7,i-7,0,2,test$dose[i-7])
      mbl <- test1[which(test1$id==(i-7) & test1$visittime<2),]
    }else if((test$status[i-7]==0 & test$time[i-7]==1)){
      vbl <- c(i-7,i-7,0,1,test$dose[i-7])
      mbl <- test1[which(test1$id==(i-7) & test1$visittime<1),]
    }else if((test$status[i-7]==0 & test$time[i-7]==6) | (test$status[i-7]==1 & test$time[i-7]==5) |
             (test$status[i-7]==1 & test$time[i-7]==6)){
      vbl <- c(i-7,i-7,1,5,test$dose[i-7])
      mbl <- test1[which(test1$id==(i-7) & test1$visittime<5),]
    }else if((test$status[i-7]==1 & test$time[i-7]==2)){
      vbl <- c(i-7,i-7,1,2,test$dose[i-7])
      mbl <- test1[which(test1$id==(i-7) & test1$visittime<2),]
    }else if((test$status[i-7]==1 & test$time[i-7]==3)){
      vbl <- c(i-7,i-7,1,3,test$dose[i-7])
      mbl <- test1[which(test1$id==(i-7) & test1$visittime<3),]
    }else if((test$status[i-7]==1 & test$time[i-7]==4)){
      vbl <- c(i-7,i-7,1,4,test$dose[i-7])
      mbl <- test1[which(test1$id==(i-7) & test1$visittime<4),]
    }else{
      vbl <- c(i-7,i-7,1,1,test$dose[i-7])
      mbl <- test1[which(test1$id==(i-7) & test1$visittime<1),]
    }
  }  
  
  if(isTRUE((test$counter[i]-test$id[i-8])==8)){
    if((test$status[i-8]==0 & test$time[i-8]==5) ){
      vbl1 <- c(i-8,i-8,0,5,test$dose[i-8])
      mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<5),]
    }else if((test$status[i-8]==0 & test$time[i-8]==4) ){
      vbl1 <- c(i-8,i-8,0,4,test$dose[i-8])
      mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<4),]
    }else if((test$status[i-8]==0 & test$time[i-8]==3) ){
      vbl1 <- c(i-8,i-8,0,3,test$dose[i-8])
      mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<3),]
    }else if ((test$status[i-8]==0 & test$time[i-8]==2)){
      vbl1 <- c(i-8,i-8,0,2,test$dose[i-8])
      mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<2),]
    }else if((test$status[i-8]==0 & test$time[i-8]==1)){
      vbl1 <- c(i-8,i-8,0,1,test$dose[i-8])
      mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<1),]
    }else if((test$status[i-8]==0 & test$time[i-8]==6) | (test$status[i-8]==1 & test$time[i-8]==5) |
             (test$status[i-8]==1 & test$time[i-8]==6)){
      vbl1 <- c(i-8,i-8,1,5,test$dose[i-8])
      mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<5),]
    }else if((test$status[i-8]==1 & test$time[i-8]==2)){
      vbl1 <- c(i-8,i-8,1,2,test$dose[i-8])
      mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<2),]
    }else if((test$status[i-8]==1 & test$time[i-8]==3)){
      vbl1 <- c(i-8,i-8,1,3,test$dose[i-8])
      mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<3),]
    }else if((test$status[i-8]==1 & test$time[i-8]==4)){
      vbl1 <- c(i-8,i-8,1,4,test$dose[i-8])
      mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<4),]
    }else{
      vbl1 <- c(i-8,i-8,1,1,test$dose[i-8])
      mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<1),]
    }
  }  
  
  
  ver <- c(!is.null(vb[1]), !is.null(vb1[1]), !is.null(vbt[1]), !is.null(vbt1[1]), !is.null(vbc[1]), !is.null(vbc1[1]), !is.null(vbk[1]), 
           !is.null(vbk1[1]), !is.null(vbl[1]),!is.null(vbl1[1]))
  
  
  testtemp <- test[,]
  testtemp <- test[ which(test$id<=((i+1)-(sum(ver, na.rm=TRUE)))),] 
  
  testtemp <- rbind(testtemp, vbl1) 
  testtemp <- rbind(testtemp, vbl) 
  testtemp <- rbind(testtemp, vbk1) 
  testtemp <- rbind(testtemp, vbk) 
  testtemp <- rbind(testtemp, vbc1) 
  testtemp <- rbind(testtemp, vbc) 
  testtemp <- rbind(testtemp, vbt1) 
  testtemp <- rbind(testtemp, vbt) 
  testtemp <- rbind(testtemp, vb) 
  testtemp <- rbind(testtemp, vb1)
  
  colnames(testtemp)[1] <- "id"
  colnames(testtemp)[2] <- "counter"
  colnames(testtemp)[3] <- "status"
  colnames(testtemp)[4] <- "time"
  colnames(testtemp)[5] <- "dose"
  
  
  test1temp <- test1[,]
  test1temp <- test1[which(test1$id<=((i+1)-(sum(ver, na.rm=TRUE)))),]
  test1temp <- rbind(test1temp, mbl1) 
  test1temp <- rbind(test1temp, mbl) 
  test1temp <- rbind(test1temp, mbk1) 
  test1temp <- rbind(test1temp, mbk) 
  test1temp <- rbind(test1temp, mbc1) 
  test1temp <- rbind(test1temp, mbc) 
  test1temp <- rbind(test1temp, mbt1) 
  test1temp <- rbind(test1temp, mbt) 
  test1temp <- rbind(test1temp, mb) 
  test1temp <- rbind(test1temp, mb1) 
  
  colnames(test1temp)[1] <- "id"
  colnames(test1temp)[2] <- "y"
  colnames(test1temp)[3] <- "visittime"
  colnames(test1temp)[4] <- "dose"
  
  out <- list(survtemp=testtemp, longtemp=test1temp)
}




#############################################################
#####                      2+2                          #####

escalate<- function(meanpars,varpars, gampars, nt, n1, d, nvisit, dose, censor, pl, stopping_rule, inter_miss){
  # This function is for the 2+2 design
  # Assigns sequentially cohorts of subjects to doses, till the 2+2 is finished
  # A stopping rule for excessive toxicity is applied on the first 15 subjects and for the first treatment cycle only
  # When we observe two DLTs and two non-DLTs in two different cycles, the 2+2 design will stop to proceed to the joint modeling or the CRM

  toxic <- NULL # it takes value 0 if trial does not stop for excessive toxicity and 1 otherwise
  x<-rep(NA,2) # empty vector for the dose level that will be assigned to every person
  test <- makedata(meanpars, varpars, gampars, nt, n1, dose, d, nvisit, pl)   
  data <- test$surv
  
  id<-c(1:2) 
  counter <-c(1:2) 
  x[1:2]<- dose[1]  # the first 2 subjects start by dose level 1
  c1 <-c2 <- c3 <- c4 <- c5 <- c6 <- rep(NA,2) 
  for (i in 1:2){
    c1[i]<-data[which(dose==x[i]),1,i]
    c2[i]<-data[which(dose==x[i]),2,i]
    c3[i]<-data[which(dose==x[i]),3,i]
    c4[i]<-data[which(dose==x[i]),4,i]
    c5[i]<-data[which(dose==x[i]),5,i]
    c6[i]<-data[which(dose==x[i]),6,i]
  }
  
  dtime <- function(x){                
    if(sum(x)==nt) out <- nt
    else out <- min(which(x==0))
    out
  }
  surv <- matrix(c(c1,c2,c3, c4, c5, c6),nrow=2,ncol=nt)
  time <- rep(0,2)  
  
  time <- apply(surv,1,dtime) 
  status <- apply(surv,1,min)  
  if(censor){					
    censortime <- runif(2,0,nt)
    censortime <- ceiling(censortime)
    status <- ifelse(censortime < time,1,status)
    time <- apply(cbind(time,censortime),1,min)    
  }
  time <- as.integer(time) 
  
  survdata <- data.frame(id=id,counter=counter,status=status,time=time, dose=x)
  survdata <- databegin(test, survdata, nvisit, dose, inter_miss)$surv
  longdata <- databegin(test, survdata, nvisit, dose, inter_miss)$dat  
  
  # Adds 2 subjects per iteration and assigns the dose
  for (i in seq(3,n1,2)){ 
    x <- updown(dose,i,x,data,d)
    c1<-data[which(dose==x[i]),1,i]
    c11<-data[which(dose==x[i+1]),1,(i+1)]
    
    c2<-data[which(dose==x[i]),2,i]
    c22<-data[which(dose==x[i+1]),2,(i+1)]
    
    c3 <-data[which(dose==x[i]),3,i]
    c33 <-data[which(dose==x[i+1]),3,(i+1)]
    
    c4 <-data[which(dose==x[i]),4,i]
    c44 <-data[which(dose==x[i+1]),4,(i+1)]
    
    c5 <-data[which(dose==x[i]),5,i]
    c55 <-data[which(dose==x[i+1]),5,(i+1)]
    
    c6 <-data[which(dose==x[i]),6,i]
    c66 <-data[which(dose==x[i+1]),6,(i+1)]
    
    dtime <- function(x){                
      if(sum(x)==nt) out <- nt
      else out <- min(which(x==0))
      out
    }
    surv <- matrix(c(c1,c2,c3, c4, c5, c6),nrow=1,ncol=nt)
    surv1 <- matrix(c(c11, c22, c33, c44, c55, c66),nrow=1,ncol=nt)
    
    time <- rep(0,1) 
    time1 <- rep(0,1)
    time <- apply(surv,1,dtime) 
    time1 <- apply(surv1,1,dtime) 
    status <- apply(surv,1,min) 
    status1 <- apply(surv1,1,min)
    
    # Uniform censoring
    if(censor){					
      censortime <- runif(1,0,nt)
      censortime <- ceiling(censortime)
      status <- ifelse(censortime < time,1,status)
      time <- apply(cbind(time,censortime),1,min)    
    }
    
    # Uniform censoring
    if(censor){					
      censortime <- runif(1,0,nt)
      censortime <- ceiling(censortime)
      status1 <- ifelse(censortime < time1,1,status1)
      time1 <- apply(cbind(time1,censortime),1,min)    
    }
    
    
    addin <- c(i, i, status, time, x[i]) 
    addin1 <- c((i+1), (i+1), status1, time1, x[i+1])
    survdata = rbind(survdata, addin)
    survdata <- data2plcrmjm(test, survdata, longdata, nvisit, dose, inter_miss)$surv
    longdata <- data2plcrmjm(test, survdata, longdata, nvisit, dose, inter_miss)$dat
    survdata = rbind(survdata, addin1)
    survdata <- data2plcrmjm(test, survdata, longdata, nvisit, dose, inter_miss)$surv
    longdata <- data2plcrmjm(test, survdata, longdata, nvisit, dose, inter_miss)$dat
    survtemp <- sequentialentry0(survdata, longdata, i)$survtemp
    longtemp <- sequentialentry0(survdata, longdata, i)$longtemp
    # stopping rule for excessive toxicity
    count_DLTs_i <- NROW(which((survtemp[1:i, 3]==0) & (survtemp[1:i, 4]==1) & (survtemp[1:i, 5]==1.2)))
    count_DLTs_ii <- NROW(which((survtemp$status==0) & (survtemp$time==1) & (survtemp$dose==1.2)))
    if((count_DLTs_i >= stopping_rule[i]) |  (count_DLTs_ii >= stopping_rule[i+1])){
      toxic <- 1 
      break
    }
    
    
    toxic <- 0
    # Summary of DLTs and non DLTs per cycle
    perc <- table(survdata$status, survdata$time)
    if(nrow(perc)==1){
      next
    }
    
    
    chose0 <- as.vector(perc[1,])
    chose1 <- as.vector(perc[2,])
    
    for(j in 1:length(chose0)){
      if(chose0[j]>0){
        chose0[j] <- 1 
      }
    }
    
    
    for(j in 1:length(chose1)){
      if(chose1[j]>0){
        chose1[j] <- 1 
      }
    }
    
    
    if(sum(chose0)>1 & sum(chose1)>1){
      perctemp <- table(survtemp$status, survtemp$time)
      if(nrow(perctemp)==1){
        next
      }
      
      
      chose2 <- as.vector(perctemp[1,])
      chose3 <- as.vector(perctemp[2,])
      
      
      for(j in 1:length(chose2)){
        if(chose2[j]>0){
          chose2[j] <- 1 
        }
      }
      
      
      for(j in 1:length(chose3)){
        if(chose3[j]>0){
          chose3[j] <- 1 
        }
      }
      
      if(sum(chose2)>1 & sum(chose3)>1){
        break
      }
    }
  }
  
  # If too many subjects are treated with the 2+2 stop and recreate data. No point to continue with model based.
  if(nrow(survdata)>=(round((5/6)*(n1)))){   
    return(escalate(meanpars, varpars, gampars, nt, n1, d, nvisit, dose, censor, pl, stopping_rule, inter_miss))
  }
  out <- list(survdata=survdata,  survtemp=survtemp, rawsurv=data, rawlong=test$rawlong, longdata=longdata, longtemp=longtemp, toxic= toxic)
  
}



#############################################################
#####          CALCULATES THE CRM LOGLIKELIHOOD         #####

loglikcrm<- function(pars1, ses=NULL, surv, pmnormtol=1e-3, printpars=F, fix_param){
  # Calculates the log likelihood of the crm model with two parameters estimated (intercept and time)
  # pmnormtol specifies the accuracy to which multivariate normal CDF's should be calculated
  # printpars to print parameters for all calls to loglikcrm evaluated by nlminb
  
  if(!is.null(ses))	pars1 <- pars1*ses
  if(printpars) print(pars)
  n <- dim(surv)[1]
  nt <- max(surv$time) # max number of time intervals
  pyt <- dim(surv)[2] - 1 
  bt <- pars1[1:(pyt-2)] # point estimates of the survival data (intercept and time)
  
  status <- surv$status
  time <- surv$time  # max survival time for every subject
  X.t <- surv[,-(1:4)] # creates a data frame with only the dose from the survival data
  
  cov.const <- bt[1] + as.matrix(X.t) %*% as.matrix(fix_param) # dose slope is fixed 
  betaX <- matrix(cov.const,nrow=nt,ncol=n,byrow=T) + bt[2]*matrix(0:(nt-1),nrow=nt,ncol=n)
  q.Phi <- betaX
  
  # calculate the multivariate normal CDFs and put results in Phi,
  # Phi[,i] is the probability of surviving to the start of the ith interval
  Phi <- matrix(0,nrow=n,ncol=(nt+1)) 
  # First column of Phi is the probability of surviving to time 0
  Phi[,1] <- 1
  I <- array(diag(nt),dim=c(nt,nt,n)) 
  
  for(i in 1:nt){
    Inew <- I[1:i,1:i,]
    
    Sig.Phi <- Inew
    mu.Phi <- rep(0,i)
    
    # only calculate the components of Phi that are actually going to be used
    # (i.e. time for censored individuals and time and (time - 1) for deaths)
    ids1 <- c(which((time==i)&status==0),which((time==(i+1))&status==0))
    ids2 <- which((time==i)&status==1)
    ids <- c(ids1,ids2)
    
    p.Phi <- rep(0,n)
    if(i==1){
      p.Phi[ids] <- pnorm(q.Phi[1,ids],mean=0,sd=sqrt(Sig.Phi[ids]))
    } else {
      for(k in ids){
        p.Phi[k] <- pmnorm(q.Phi[1:i,k],mu.Phi,Sig.Phi[,,k],abseps=pmnormtol)
      }
    }
    Phi[,(i+1)] <- p.Phi
  }
  
  # Psi[,i] is the probability of dying in time interval i
  Psi <- matrix(nrow=n,ncol=nt)
  for(i in 1:nt){
    Psi[,i] <- Phi[,i] - Phi[,(i+1)]
  }
  
  # calculate all the terms in the likelihood
  term1 <- 0
  for(i in 1:nt){
    ids <- which((time==i)&status==0)
    term1 <- term1 + sum(log(Psi[ids,i]))
  }
  term2 <- 0
  for(i in 1:nt){
    ids <- which((time==i)&status==1)
    term2 <- term2 + sum(log(Phi[ids,(i+1)]))
  }
  
  loglik <- term1 + term2
  # Return minus loglik - optimisation routines will minimise this
  if (is.finite(loglik)) {
    return(-loglik)
  } else {
    return(1e+09)
  }
  
  
}




#######################################################################
#####    CALCULATES DOSE ASSINGED TO EVERY SUBJECT IN THE CRM     #####


nextdose <- function(estim, nt, n1, d, dose, target, fix_param){
  # It calculates the dose for the next subject 
  # First it estimates the cumulative probability of toxicity per dose level 
  # Then it calculates the distance between the target toxicity level and the cumulative probability of toxicity at each dose level 
  # This function calculates the MTD not the OD
  
  # estim are the parameters estimated from the CRM
  
  
  tox <- rep(NA,d)
  
  for(i in 1:d){
    tox[i] <- ((1-pnorm(estim[1] + fix_param*dose[i]))+ (pnorm(estim[1] + fix_param*dose[i]))*(1-pnorm(estim[1] + fix_param*dose[i] +estim[2]))+ 
                 (pnorm(estim[1] + fix_param*dose[i]))*(pnorm(estim[1] + fix_param*dose[i]+estim[2]))*(1-pnorm(estim[1] + fix_param*dose[i]+ estim[2]*2)) +
                 (pnorm(estim[1] + fix_param*dose[i])*(pnorm(estim[1] + fix_param*dose[i]+estim[2]))*(pnorm(estim[1] + fix_param*dose[i]+estim[2]*2))*(1-pnorm(estim[1] + fix_param*dose[i]+estim[2]*3)))+
                 (pnorm(estim[1] + fix_param*dose[i])*(pnorm(estim[1] + fix_param*dose[i]+estim[2]))*(pnorm(estim[1] + fix_param*dose[i]+estim[2]*2))*
                    (pnorm(estim[1] + fix_param*dose[i]+estim[2]*3))*(1-pnorm(estim[1] + fix_param*dose[i]+estim[2]*4)))+
                 (pnorm(estim[1] + fix_param*dose[i])*(pnorm(estim[1] + fix_param*dose[i]+estim[2]))*(pnorm(estim[1] + fix_param*dose[i]+estim[2]*2))*
                    (pnorm(estim[1] + fix_param*dose[i]+estim[2]*3))*(pnorm(estim[1] + fix_param*dose[i]+estim[2]*4))*(1-pnorm(estim[1] + fix_param*dose[i]+estim[2]*5))))
  }
  
  dist <-(abs(tox-target)) 
  out <- which.min(dist)
  
} 



#######################################################################
#####   SEQUENTIAL ENTRY PER CYCLE FOR SUBJECTS IN CRM AND JM     #####

sequentialentry <- function(test, test1, i, count){
  # Similar function to sequentialentry0
  # Analyzes subject data during the CRM and JM 
  vb <- vbt <- vbt1 <- vbc <- vbc1 <- vbk <- vbk1 <- vbl <- vbl1 <- NULL # for survival objects
  mb <- mbt <- mbt1 <- mbc <- mbc1 <- mbk <- mbk1 <- mbl <- mbl1 <- NULL # for longitudinal objects
  
  if(isTRUE(i-count==1)){
    
    if(isTRUE(test$id[i]==test$counter[i])){ 
      if((test$status[i]==0 & test$time[i]==3) | (test$status[i]==0 & test$time[i]==2) |(test$status[i]==1 & test$time[i]==3) |
         (test$status[i]==1 & test$time[i]==2) | (test$status[i]==1 & test$time[i]==1) | (test$status[i]==0 & test$time[i]==4) |
         (test$status[i]==0 & test$time[i]==5) | (test$status[i]==0 & test$time[i]==6)| (test$status[i]==1 & test$time[i]==4) | 
         (test$status[i]==1 & test$time[i]==5) | (test$status[i]==1 & test$time[i]==6) ){
        vb<- c(i,i,1,1,test$dose[i])
        mb <- test1[which(test1$id==i & test1$visittime<1),]
      }else{
        vb <- c(i,i,0,1,test$dose[i])
        mb <- test1[which(test1$id==i & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-1])==1)){
      if((test$status[i-1]==0 & test$time[i-1]==3) | (test$status[i-1]==0 & test$time[i-1]==4) |
         (test$status[i-1]==0 & test$time[i-1]==5) | (test$status[i-1]==0 & test$time[i-1]==6) ){
        vbt<- c(i-1,i-1,1,2,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<2),]
      }else if ((test$status[i-1]==0 & test$time[i-1]==2)){
        vbt <- c(i-1,i-1,0,2,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<2),]
      }else if((test$status[i-1]==0 & test$time[i-1]==1)){
        vbt <- c(i-1,i-1,0,1,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<1),]
      }else if((test$status[i-1]==1 & test$time[i-1]==2) | (test$status[i-1]==1 & test$time[i-1]==3) |
               (test$status[i-1]==1 & test$time[i-1]==4) | (test$status[i-1]==1 & test$time[i-1]==5) |
               (test$status[i-1]==1 & test$time[i-1]==6)){
        vbt <- c(i-1,i-1,1,2,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<2),]
      }else{
        vbt <- c(i-1,i-1,1,1,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-2])==2)){
      if((test$status[i-2]==0 & test$time[i-2]==3) | (test$status[i-2]==0 & test$time[i-2]==4) |
         (test$status[i-2]==0 & test$time[i-2]==5) | (test$status[i-2]==0 & test$time[i-2]==6) ){
        vbt1 <- c(i-2,i-2,1,2,test$dose[i-2])
        mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<2),]
      }else if ((test$status[i-2]==0 & test$time[i-2]==2)){
        vbt1 <- c(i-2,i-2,0,2,test$dose[i-2])
        mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<2),]
      }else if((test$status[i-2]==0 & test$time[i-2]==1)){
        vbt1 <- c(i-2,i-2,0,1,test$dose[i-2])
        mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<1),]
      }else if((test$status[i-2]==1 & test$time[i-2]==2) | (test$status[i-2]==1 & test$time[i-2]==3) |
               (test$status[i-2]==1 & test$time[i-2]==4) | (test$status[i-2]==1 & test$time[i-2]==5) |
               (test$status[i-2]==1 & test$time[i-2]==6)){
        vbt1 <- c(i-2,i-2,1,2,test$dose[i-2])
        mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<2),]
      }else{
        vbt1 <- c(i-2,i-2,1,1,test$dose[i-2])
        mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<1),]
      }
    } 
    
    if(isTRUE((test$counter[i]-test$id[i-3])==3)){
      if((test$status[i-3]==0 & test$time[i-3]==3) ){
        vbc<- c(i-3,i-3,0,3,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<3),]
      }else if ((test$status[i-3]==0 & test$time[i-3]==2)){
        vbc <- c(i-3,i-3,0,2,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<2),]
      }else if((test$status[i-3]==0 & test$time[i-3]==1)){
        vbc <- c(i-3,i-3,0,1,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<1),]
      }else if((test$status[i-3]==0 & test$time[i-3]==4) | (test$status[i-3]==0 & test$time[i-3]==5) |
               (test$status[i-3]==0 & test$time[i-3]==6) | (test$status[i-3]==1 & test$time[i-3]==3) |
               (test$status[i-3]==1 & test$time[i-3]==4) | (test$status[i-3]==1 & test$time[i-3]==5) |
               (test$status[i-3]==1 & test$time[i-3]==6)){
        vbc <- c(i-3,i-3,1,3,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<3),]
      }else if((test$status[i-3]==1 & test$time[i-3]==2)){
        vbc <- c(i-3,i-3,1,2,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<2),]
      }else{
        vbc <- c(i-3,i-3,1,1,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-4])==4)){
      if((test$status[i-4]==0 & test$time[i-4]==3) ){
        vbc1<- c(i-4,i-4,0,3,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<3),]
      }else if ((test$status[i-4]==0 & test$time[i-4]==2)){
        vbc1 <- c(i-4,i-4,0,2,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<2),]
      }else if((test$status[i-4]==0 & test$time[i-4]==1)){
        vbc1 <- c(i-4,i-4,0,1,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<1),]
      }else if((test$status[i-4]==0 & test$time[i-4]==4) | (test$status[i-4]==0 & test$time[i-4]==5) |
               (test$status[i-4]==0 & test$time[i-4]==6) | (test$status[i-4]==1 & test$time[i-4]==3) |
               (test$status[i-4]==1 & test$time[i-4]==4) | (test$status[i-4]==1 & test$time[i-4]==5) |
               (test$status[i-4]==1 & test$time[i-4]==6)){
        vbc1 <- c(i-4,i-4,1,3,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<3),]
      }else if((test$status[i-4]==1 & test$time[i-4]==2)){
        vbc1 <- c(i-4,i-4,1,2,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<2),]
      }else{
        vbc1 <- c(i-4,i-4,1,1,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-5])==5)){
      if((test$status[i-5]==0 & test$time[i-5]==4) ){
        vbk <- c(i-5,i-5,0,4,test$dose[i-5])
        mbk <- test1[which(test1$id==(i-5) & test1$visittime<4),]
      }else if((test$status[i-5]==0 & test$time[i-5]==3) ){
        vbk <- c(i-5,i-5,0,3,test$dose[i-5])
        mbk <- test1[which(test1$id==(i-5) & test1$visittime<3),]
      }else if ((test$status[i-5]==0 & test$time[i-5]==2)){
        vbk <- c(i-5,i-5,0,2,test$dose[i-5])
        mbk <- test1[which(test1$id==(i-5) & test1$visittime<2),]
      }else if((test$status[i-5]==0 & test$time[i-5]==1)){
        vbk <- c(i-5,i-5,0,1,test$dose[i-5])
        mbk <- test1[which(test1$id==(i-5) & test1$visittime<1),]
      }else if((test$status[i-5]==0 & test$time[i-5]==5) |
               (test$status[i-5]==0 & test$time[i-5]==6) |
               (test$status[i-5]==1 & test$time[i-5]==4) | (test$status[i-5]==1 & test$time[i-5]==5) |
               (test$status[i-5]==1 & test$time[i-5]==6)){
        vbk <- c(i-5,i-5,1,4,test$dose[i-5])
        mbk <- test1[which(test1$id==(i-5) & test1$visittime<4),]
      }else if((test$status[i-5]==1 & test$time[i-5]==2)){
        vbk <- c(i-5,i-5,1,2,test$dose[i-5])
        mbk <- test1[which(test1$id==(i-5) & test1$visittime<2),]
      }else if((test$status[i-5]==1 & test$time[i-5]==3)){
        vbk <- c(i-5,i-5,1,3,test$dose[i-5])
        mbk <- test1[which(test1$id==(i-5) & test1$visittime<3),]
      }else{
        vbk <- c(i-5,i-5,1,1,test$dose[i-5])
        mbk <- test1[which(test1$id==(i-5) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-6])==6)){
      if((test$status[i-6]==0 & test$time[i-6]==4) ){
        vbk1 <- c(i-6,i-6,0,4,test$dose[i-6])
        mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<4),]
      }else if((test$status[i-6]==0 & test$time[i-6]==3) ){
        vbk1 <- c(i-6,i-6,0,3,test$dose[i-6])
        mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<3),]
      }else if ((test$status[i-6]==0 & test$time[i-6]==2)){
        vbk1 <- c(i-6,i-6,0,2,test$dose[i-6])
        mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<2),]
      }else if((test$status[i-6]==0 & test$time[i-6]==1)){
        vbk1 <- c(i-6,i-6,0,1,test$dose[i-6])
        mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<1),]
      }else if((test$status[i-6]==0 & test$time[i-6]==5) |
               (test$status[i-6]==0 & test$time[i-6]==6) |
               (test$status[i-6]==1 & test$time[i-6]==4) | (test$status[i-6]==1 & test$time[i-6]==5) |
               (test$status[i-6]==1 & test$time[i-6]==6)){
        vbk1 <- c(i-6,i-6,1,4,test$dose[i-6])
        mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<4),]
      }else if((test$status[i-6]==1 & test$time[i-6]==2)){
        vbk1 <- c(i-6,i-6,1,2,test$dose[i-6])
        mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<2),]
      }else if((test$status[i-6]==1 & test$time[i-6]==3)){
        vbk1 <- c(i-6,i-6,1,3,test$dose[i-6])
        mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<3),]
      }else{
        vbk1 <- c(i-6,i-6,1,1,test$dose[i-6])
        mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-7])==7)){
      if((test$status[i-7]==0 & test$time[i-7]==5) ){
        vbl <- c(i-7,i-7,0,5,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<5),]
      }else if((test$status[i-7]==0 & test$time[i-7]==4) ){
        vbl <- c(i-7,i-7,0,4,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<4),]
      }else if((test$status[i-7]==0 & test$time[i-7]==3) ){
        vbl <- c(i-7,i-7,0,3,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<3),]
      }else if ((test$status[i-7]==0 & test$time[i-7]==2)){
        vbl <- c(i-7,i-7,0,2,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<2),]
      }else if((test$status[i-7]==0 & test$time[i-7]==1)){
        vbl <- c(i-7,i-7,0,1,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<1),]
      }else if((test$status[i-7]==0 & test$time[i-7]==6) | (test$status[i-7]==1 & test$time[i-7]==5) |
               (test$status[i-7]==1 & test$time[i-7]==6)){
        vbl <- c(i-7,i-7,1,5,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<5),]
      }else if((test$status[i-7]==1 & test$time[i-7]==2)){
        vbl <- c(i-7,i-7,1,2,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<2),]
      }else if((test$status[i-7]==1 & test$time[i-7]==3)){
        vbl <- c(i-7,i-7,1,3,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<3),]
      }else if((test$status[i-7]==1 & test$time[i-7]==4)){
        vbl <- c(i-7,i-7,1,4,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<4),]
      }else{
        vbl <- c(i-7,i-7,1,1,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<1),]
      }
    }  
    
    if(isTRUE((test$counter[i]-test$id[i-8])==8)){
      if((test$status[i-8]==0 & test$time[i-8]==5) ){
        vbl1 <- c(i-8,i-8,0,5,test$dose[i-8])
        mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<5),]
      }else if((test$status[i-8]==0 & test$time[i-8]==4) ){
        vbl1 <- c(i-8,i-8,0,4,test$dose[i-8])
        mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<4),]
      }else if((test$status[i-8]==0 & test$time[i-8]==3) ){
        vbl1 <- c(i-8,i-8,0,3,test$dose[i-8])
        mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<3),]
      }else if ((test$status[i-8]==0 & test$time[i-8]==2)){
        vbl1 <- c(i-8,i-8,0,2,test$dose[i-8])
        mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<2),]
      }else if((test$status[i-8]==0 & test$time[i-8]==1)){
        vbl1 <- c(i-8,i-8,0,1,test$dose[i-8])
        mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<1),]
      }else if((test$status[i-8]==0 & test$time[i-8]==6) | (test$status[i-8]==1 & test$time[i-8]==5) |
               (test$status[i-8]==1 & test$time[i-8]==6)){
        vbl1 <- c(i-8,i-8,1,5,test$dose[i-8])
        mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<5),]
      }else if((test$status[i-8]==1 & test$time[i-8]==2)){
        vbl1 <- c(i-8,i-8,1,2,test$dose[i-8])
        mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<2),]
      }else if((test$status[i-8]==1 & test$time[i-8]==3)){
        vbl1 <- c(i-8,i-8,1,3,test$dose[i-8])
        mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<3),]
      }else if((test$status[i-8]==1 & test$time[i-8]==4)){
        vbl1 <- c(i-8,i-8,1,4,test$dose[i-8])
        mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<4),]
      }else{
        vbl1 <- c(i-8,i-8,1,1,test$dose[i-8])
        mbl1 <- test1[which(test1$id==(i-8) & test1$visittime<1),]
      }
    }  
  }
  
  if(isTRUE(i-count==2)){
    
    if(isTRUE(test$id[i]==test$counter[i])){ 
      if((test$status[i]==0 & test$time[i]==3) | (test$status[i]==0 & test$time[i]==2) |(test$status[i]==1 & test$time[i]==3) |
         (test$status[i]==1 & test$time[i]==2) | (test$status[i]==1 & test$time[i]==1) | (test$status[i]==0 & test$time[i]==4) |
         (test$status[i]==0 & test$time[i]==5) | (test$status[i]==0 & test$time[i]==6)| (test$status[i]==1 & test$time[i]==4) | 
         (test$status[i]==1 & test$time[i]==5) | (test$status[i]==1 & test$time[i]==6) ){
        vb<- c(i,i,1,1,test$dose[i])
        mb <- test1[which(test1$id==i & test1$visittime<1),]
      }else{
        vb <- c(i,i,0,1,test$dose[i])
        mb <- test1[which(test1$id==i & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-1])==1)){
      if((test$status[i-1]==0 & test$time[i-1]==3) | (test$status[i-1]==0 & test$time[i-1]==4) |
         (test$status[i-1]==0 & test$time[i-1]==5) | (test$status[i-1]==0 & test$time[i-1]==6) ){
        vbt<- c(i-1,i-1,1,2,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<2),]
      }else if ((test$status[i-1]==0 & test$time[i-1]==2)){
        vbt <- c(i-1,i-1,0,2,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<2),]
      }else if((test$status[i-1]==0 & test$time[i-1]==1)){
        vbt <- c(i-1,i-1,0,1,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<1),]
      }else if((test$status[i-1]==1 & test$time[i-1]==2) | (test$status[i-1]==1 & test$time[i-1]==3) |
               (test$status[i-1]==1 & test$time[i-1]==4) | (test$status[i-1]==1 & test$time[i-1]==5) |
               (test$status[i-1]==1 & test$time[i-1]==6)){
        vbt <- c(i-1,i-1,1,2,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<2),]
      }else{
        vbt <- c(i-1,i-1,1,1,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-2])==2)){
      if((test$status[i-2]==0 & test$time[i-2]==3) ){
        vbt1<- c(i-2,i-2,0,3,test$dose[i-2])
        mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<3),]
      }else if ((test$status[i-2]==0 & test$time[i-2]==2)){
        vbt1 <- c(i-2,i-2,0,2,test$dose[i-2])
        mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<2),]
      }else if((test$status[i-2]==0 & test$time[i-2]==1)){
        vbt1 <- c(i-2,i-2,0,1,test$dose[i-2])
        mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<1),]
      }else if((test$status[i-2]==0 & test$time[i-2]==4) | (test$status[i-2]==0 & test$time[i-2]==5) |
               (test$status[i-2]==0 & test$time[i-2]==6) | (test$status[i-2]==1 & test$time[i-2]==3) |
               (test$status[i-2]==1 & test$time[i-2]==4) | (test$status[i-2]==1 & test$time[i-2]==5) |
               (test$status[i-2]==1 & test$time[i-2]==6)){
        vbt1 <- c(i-2,i-2,1,3,test$dose[i-2])
        mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<3),]
      }else if((test$status[i-2]==1 & test$time[i-2]==2)){
        vbt1 <- c(i-2,i-2,1,2,test$dose[i-2])
        mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<2),]
      }else{
        vbt1 <- c(i-2,i-2,1,1,test$dose[i-2])
        mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<1),]
      }
    }  
    
    if(isTRUE((test$counter[i]-test$id[i-3])==3)){
      if((test$status[i-3]==0 & test$time[i-3]==3) ){
        vbc<- c(i-3,i-3,0,3,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<3),]
      }else if ((test$status[i-3]==0 & test$time[i-3]==2)){
        vbc <- c(i-3,i-3,0,2,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<2),]
      }else if((test$status[i-3]==0 & test$time[i-3]==1)){
        vbc <- c(i-3,i-3,0,1,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<1),]
      }else if((test$status[i-3]==0 & test$time[i-3]==4) | (test$status[i-3]==0 & test$time[i-3]==5) |
               (test$status[i-3]==0 & test$time[i-3]==6) | (test$status[i-3]==1 & test$time[i-3]==3) |
               (test$status[i-3]==1 & test$time[i-3]==4) | (test$status[i-3]==1 & test$time[i-3]==5) |
               (test$status[i-3]==1 & test$time[i-3]==6)){
        vbc <- c(i-3,i-3,1,3,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<3),]
      }else if((test$status[i-3]==1 & test$time[i-3]==2)){
        vbc <- c(i-3,i-3,1,2,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<2),]
      }else{
        vbc <- c(i-3,i-3,1,1,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-4])==4)){
      if((test$status[i-4]==0 & test$time[i-4]==4) ){
        vbc1 <- c(i-4,i-4,0,4,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<4),]
      }else if((test$status[i-4]==0 & test$time[i-4]==3) ){
        vbc1 <- c(i-4,i-4,0,3,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<3),]
      }else if ((test$status[i-4]==0 & test$time[i-4]==2)){
        vbc1 <- c(i-4,i-4,0,2,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<2),]
      }else if((test$status[i-4]==0 & test$time[i-4]==1)){
        vbc1 <- c(i-4,i-4,0,1,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<1),]
      }else if((test$status[i-4]==0 & test$time[i-4]==5) |
               (test$status[i-4]==0 & test$time[i-4]==6) |
               (test$status[i-4]==1 & test$time[i-4]==4) | (test$status[i-4]==1 & test$time[i-4]==5) |
               (test$status[i-4]==1 & test$time[i-4]==6)){
        vbc1 <- c(i-4,i-4,1,4,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<4),]
      }else if((test$status[i-4]==1 & test$time[i-4]==2)){
        vbc1 <- c(i-4,i-4,1,2,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<2),]
      }else if((test$status[i-4]==1 & test$time[i-4]==3)){
        vbc1 <- c(i-4,i-4,1,3,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<3),]
      }else{
        vbc1 <- c(i-4,i-4,1,1,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-5])==5)){
      if((test$status[i-5]==0 & test$time[i-5]==4) ){
        vbk <- c(i-5,i-5,0,4,test$dose[i-5])
        mbk <- test1[which(test1$id==(i-5) & test1$visittime<4),]
      }else if((test$status[i-5]==0 & test$time[i-5]==3) ){
        vbk <- c(i-5,i-5,0,3,test$dose[i-5])
        mbk <- test1[which(test1$id==(i-5) & test1$visittime<3),]
      }else if ((test$status[i-5]==0 & test$time[i-5]==2)){
        vbk <- c(i-5,i-5,0,2,test$dose[i-5])
        mbk <- test1[which(test1$id==(i-5) & test1$visittime<2),]
      }else if((test$status[i-5]==0 & test$time[i-5]==1)){
        vbk <- c(i-5,i-5,0,1,test$dose[i-5])
        mbk <- test1[which(test1$id==(i-5) & test1$visittime<1),]
      }else if((test$status[i-5]==0 & test$time[i-5]==5) |
               (test$status[i-5]==0 & test$time[i-5]==6) |
               (test$status[i-5]==1 & test$time[i-5]==4) | (test$status[i-5]==1 & test$time[i-5]==5) |
               (test$status[i-5]==1 & test$time[i-5]==6)){
        vbk <- c(i-5,i-5,1,4,test$dose[i-5])
        mbk <- test1[which(test1$id==(i-5) & test1$visittime<4),]
      }else if((test$status[i-5]==1 & test$time[i-5]==2)){
        vbk <- c(i-5,i-5,1,2,test$dose[i-5])
        mbk <- test1[which(test1$id==(i-5) & test1$visittime<2),]
      }else if((test$status[i-5]==1 & test$time[i-5]==3)){
        vbk <- c(i-5,i-5,1,3,test$dose[i-5])
        mbk <- test1[which(test1$id==(i-5) & test1$visittime<3),]
      }else{
        vbk <- c(i-5,i-5,1,1,test$dose[i-5])
        mbk <- test1[which(test1$id==(i-5) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-6])==6)){
      if((test$status[i-6]==0 & test$time[i-6]==5) ){
        vbk1 <- c(i-6,i-6,0,5,test$dose[i-6])
        mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<5),]
      }else if((test$status[i-6]==0 & test$time[i-6]==4) ){
        vbk1 <- c(i-6,i-6,0,4,test$dose[i-6])
        mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<4),]
      }else if((test$status[i-6]==0 & test$time[i-6]==3) ){
        vbk1 <- c(i-6,i-6,0,3,test$dose[i-6])
        mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<3),]
      }else if ((test$status[i-6]==0 & test$time[i-6]==2)){
        vbk1 <- c(i-6,i-6,0,2,test$dose[i-6])
        mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<2),]
      }else if((test$status[i-6]==0 & test$time[i-6]==1)){
        vbk1 <- c(i-6,i-6,0,1,test$dose[i-6])
        mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<1),]
      }else if((test$status[i-6]==0 & test$time[i-6]==6) | (test$status[i-6]==1 & test$time[i-6]==5) |
               (test$status[i-6]==1 & test$time[i-6]==6)){
        vbk1 <- c(i-6,i-6,1,5,test$dose[i-6])
        mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<5),]
      }else if((test$status[i-6]==1 & test$time[i-6]==2)){
        vbk1 <- c(i-6,i-6,1,2,test$dose[i-6])
        mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<2),]
      }else if((test$status[i-6]==1 & test$time[i-6]==3)){
        vbk1 <- c(i-6,i-6,1,3,test$dose[i-6])
        mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<3),]
      }else if((test$status[i-6]==1 & test$time[i-6]==4)){
        vbk1 <- c(i-6,i-6,1,4,test$dose[i-6])
        mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<4),]
      }else{
        vbk1 <- c(i-6,i-6,1,1,test$dose[i-6])
        mbk1 <- test1[which(test1$id==(i-6) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-7])==7)){
      if((test$status[i-7]==0 & test$time[i-7]==5) ){
        vbl <- c(i-7,i-7,0,5,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<5),]
      }else if((test$status[i-7]==0 & test$time[i-7]==4) ){
        vbl <- c(i-7,i-7,0,4,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<4),]
      }else if((test$status[i-7]==0 & test$time[i-7]==3) ){
        vbl <- c(i-7,i-7,0,3,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<3),]
      }else if ((test$status[i-7]==0 & test$time[i-7]==2)){
        vbl <- c(i-7,i-7,0,2,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<2),]
      }else if((test$status[i-7]==0 & test$time[i-7]==1)){
        vbl <- c(i-7,i-7,0,1,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<1),]
      }else if((test$status[i-7]==0 & test$time[i-7]==6) | (test$status[i-7]==1 & test$time[i-7]==5) |
               (test$status[i-7]==1 & test$time[i-7]==6)){
        vbl <- c(i-7,i-7,1,5,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<5),]
      }else if((test$status[i-7]==1 & test$time[i-7]==2)){
        vbl <- c(i-7,i-7,1,2,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<2),]
      }else if((test$status[i-7]==1 & test$time[i-7]==3)){
        vbl <- c(i-7,i-7,1,3,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<3),]
      }else if((test$status[i-7]==1 & test$time[i-7]==4)){
        vbl <- c(i-7,i-7,1,4,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<4),]
      }else{
        vbl <- c(i-7,i-7,1,1,test$dose[i-7])
        mbl <- test1[which(test1$id==(i-7) & test1$visittime<1),]
      }
    } 
  }
  
  if(isTRUE(i-count==3)){
    if(isTRUE(test$id[i]==test$counter[i])){ 
      if((test$status[i]==0 & test$time[i]==3) | (test$status[i]==0 & test$time[i]==2) |(test$status[i]==1 & test$time[i]==3) |
         (test$status[i]==1 & test$time[i]==2) | (test$status[i]==1 & test$time[i]==1) | (test$status[i]==0 & test$time[i]==4) |
         (test$status[i]==0 & test$time[i]==5) | (test$status[i]==0 & test$time[i]==6)| (test$status[i]==1 & test$time[i]==4) | 
         (test$status[i]==1 & test$time[i]==5) | (test$status[i]==1 & test$time[i]==6) ){
        vb<- c(i,i,1,1,test$dose[i])
        mb <- test1[which(test1$id==i & test1$visittime<1),]
      }else{
        vb <- c(i,i,0,1,test$dose[i])
        mb <- test1[which(test1$id==i & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-1])==1)){
      if((test$status[i-1]==0 & test$time[i-1]==3) | (test$status[i-1]==0 & test$time[i-1]==4) |
         (test$status[i-1]==0 & test$time[i-1]==5) | (test$status[i-1]==0 & test$time[i-1]==6) ){
        vbt<- c(i-1,i-1,1,2,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<2),]
      }else if ((test$status[i-1]==0 & test$time[i-1]==2)){
        vbt <- c(i-1,i-1,0,2,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<2),]
      }else if((test$status[i-1]==0 & test$time[i-1]==1)){
        vbt <- c(i-1,i-1,0,1,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<1),]
      }else if((test$status[i-1]==1 & test$time[i-1]==2) | (test$status[i-1]==1 & test$time[i-1]==3) |
               (test$status[i-1]==1 & test$time[i-1]==4) | (test$status[i-1]==1 & test$time[i-1]==5) |
               (test$status[i-1]==1 & test$time[i-1]==6)){
        vbt <- c(i-1,i-1,1,2,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<2),]
      }else{
        vbt <- c(i-1,i-1,1,1,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<1),]
      }
    } 
    
    if(isTRUE((test$counter[i]-test$id[i-2])==2)){
      if((test$status[i-2]==0 & test$time[i-2]==3) ){
        vbt1<- c(i-2,i-2,0,3,test$dose[i-2])
        mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<3),]
      }else if ((test$status[i-2]==0 & test$time[i-2]==2)){
        vbt1 <- c(i-2,i-2,0,2,test$dose[i-2])
        mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<2),]
      }else if((test$status[i-2]==0 & test$time[i-2]==1)){
        vbt1 <- c(i-2,i-2,0,1,test$dose[i-2])
        mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<1),]
      }else if((test$status[i-2]==0 & test$time[i-2]==4) | (test$status[i-2]==0 & test$time[i-2]==5) |
               (test$status[i-2]==0 & test$time[i-2]==6) | (test$status[i-2]==1 & test$time[i-2]==3) |
               (test$status[i-2]==1 & test$time[i-2]==4) | (test$status[i-2]==1 & test$time[i-2]==5) |
               (test$status[i-2]==1 & test$time[i-2]==6)){
        vbt1 <- c(i-2,i-2,1,3,test$dose[i-2])
        mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<3),]
      }else if((test$status[i-2]==1 & test$time[i-2]==2)){
        vbt1 <- c(i-2,i-2,1,2,test$dose[i-2])
        mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<2),]
      }else{
        vbt1 <- c(i-2,i-2,1,1,test$dose[i-2])
        mbt1 <- test1[which(test1$id==(i-2) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-3])==3)){
      if((test$status[i-3]==0 & test$time[i-3]==4) ){
        vbc <- c(i-3,i-3,0,4,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<4),]
      }else if((test$status[i-3]==0 & test$time[i-3]==3) ){
        vbc <- c(i-3,i-3,0,3,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<3),]
      }else if ((test$status[i-3]==0 & test$time[i-3]==2)){
        vbc <- c(i-3,i-3,0,2,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<2),]
      }else if((test$status[i-3]==0 & test$time[i-3]==1)){
        vbc <- c(i-3,i-3,0,1,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<1),]
      }else if((test$status[i-3]==0 & test$time[i-3]==5) |
               (test$status[i-3]==0 & test$time[i-3]==6) |
               (test$status[i-3]==1 & test$time[i-3]==4) | (test$status[i-3]==1 & test$time[i-3]==5) |
               (test$status[i-3]==1 & test$time[i-3]==6)){
        vbc <- c(i-3,i-3,1,4,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<4),]
      }else if((test$status[i-3]==1 & test$time[i-3]==2)){
        vbc <- c(i-3,i-3,1,2,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<2),]
      }else if((test$status[i-3]==1 & test$time[i-3]==3)){
        vbc <- c(i-3,i-3,1,3,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<3),]
      }else{
        vbc <- c(i-3,i-3,1,1,test$dose[i-3])
        mbc <- test1[which(test1$id==(i-3) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-4])==4)){
      if((test$status[i-4]==0 & test$time[i-4]==4) ){
        vbc1 <- c(i-4,i-4,0,4,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<4),]
      }else if((test$status[i-4]==0 & test$time[i-4]==3) ){
        vbc1 <- c(i-4,i-4,0,3,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<3),]
      }else if ((test$status[i-4]==0 & test$time[i-4]==2)){
        vbc1 <- c(i-4,i-4,0,2,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<2),]
      }else if((test$status[i-4]==0 & test$time[i-4]==1)){
        vbc1 <- c(i-4,i-4,0,1,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<1),]
      }else if((test$status[i-4]==0 & test$time[i-4]==5) |
               (test$status[i-4]==0 & test$time[i-4]==6) |
               (test$status[i-4]==1 & test$time[i-4]==4) | (test$status[i-4]==1 & test$time[i-4]==5) |
               (test$status[i-4]==1 & test$time[i-4]==6)){
        vbc1 <- c(i-4,i-4,1,4,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<4),]
      }else if((test$status[i-4]==1 & test$time[i-4]==2)){
        vbc1 <- c(i-4,i-4,1,2,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<2),]
      }else if((test$status[i-4]==1 & test$time[i-4]==3)){
        vbc1 <- c(i-4,i-4,1,3,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<3),]
      }else{
        vbc1 <- c(i-4,i-4,1,1,test$dose[i-4])
        mbc1 <- test1[which(test1$id==(i-4) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-5])==5)){
      if((test$status[i-5]==0 & test$time[i-5]==5) ){
        vbk1 <- c(i-5,i-5,0,5,test$dose[i-5])
        mbk1 <- test1[which(test1$id==(i-5) & test1$visittime<5),]
      }else if((test$status[i-5]==0 & test$time[i-5]==4) ){
        vbk1 <- c(i-5,i-5,0,4,test$dose[i-5])
        mbk1 <- test1[which(test1$id==(i-5) & test1$visittime<4),]
      }else if((test$status[i-5]==0 & test$time[i-5]==3) ){
        vbk1 <- c(i-5,i-5,0,3,test$dose[i-5])
        mbk1 <- test1[which(test1$id==(i-5) & test1$visittime<3),]
      }else if ((test$status[i-5]==0 & test$time[i-5]==2)){
        vbk1 <- c(i-5,i-5,0,2,test$dose[i-5])
        mbk1 <- test1[which(test1$id==(i-5) & test1$visittime<2),]
      }else if((test$status[i-5]==0 & test$time[i-5]==1)){
        vbk1 <- c(i-5,i-5,0,1,test$dose[i-5])
        mbk1 <- test1[which(test1$id==(i-5) & test1$visittime<1),]
      }else if((test$status[i-5]==0 & test$time[i-5]==6) | (test$status[i-5]==1 & test$time[i-5]==5) |
               (test$status[i-5]==1 & test$time[i-5]==6)){
        vbk1 <- c(i-5,i-5,1,5,test$dose[i-5])
        mbk1 <- test1[which(test1$id==(i-5) & test1$visittime<5),]
      }else if((test$status[i-5]==1 & test$time[i-5]==2)){
        vbk1 <- c(i-5,i-5,1,2,test$dose[i-5])
        mbk1 <- test1[which(test1$id==(i-5) & test1$visittime<2),]
      }else if((test$status[i-5]==1 & test$time[i-5]==3)){
        vbk1 <- c(i-5,i-5,1,3,test$dose[i-5])
        mbk1 <- test1[which(test1$id==(i-5) & test1$visittime<3),]
      }else if((test$status[i-5]==1 & test$time[i-5]==4)){
        vbk1 <- c(i-5,i-5,1,4,test$dose[i-5])
        mbk1 <- test1[which(test1$id==(i-5) & test1$visittime<4),]
      }else{
        vbk1 <- c(i-5,i-5,1,1,test$dose[i-5])
        mbk1 <- test1[which(test1$id==(i-5) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-6])==6)){
      if((test$status[i-6]==0 & test$time[i-6]==5) ){
        vbl1 <- c(i-6,i-6,0,5,test$dose[i-6])
        mbl1 <- test1[which(test1$id==(i-6) & test1$visittime<5),]
      }else if((test$status[i-6]==0 & test$time[i-6]==4) ){
        vbl1 <- c(i-6,i-6,0,4,test$dose[i-6])
        mbl1 <- test1[which(test1$id==(i-6) & test1$visittime<4),]
      }else if((test$status[i-6]==0 & test$time[i-6]==3) ){
        vbl1 <- c(i-6,i-6,0,3,test$dose[i-6])
        mbl1 <- test1[which(test1$id==(i-6) & test1$visittime<3),]
      }else if ((test$status[i-6]==0 & test$time[i-6]==2)){
        vbl1 <- c(i-6,i-6,0,2,test$dose[i-6])
        mbl1 <- test1[which(test1$id==(i-6) & test1$visittime<2),]
      }else if((test$status[i-6]==0 & test$time[i-6]==1)){
        vbl1 <- c(i-6,i-6,0,1,test$dose[i-6])
        mbl1 <- test1[which(test1$id==(i-6) & test1$visittime<1),]
      }else if((test$status[i-6]==0 & test$time[i-6]==6) | (test$status[i-6]==1 & test$time[i-6]==5) |
               (test$status[i-6]==1 & test$time[i-6]==6)){
        vbl1 <- c(i-6,i-6,1,5,test$dose[i-6])
        mbl1 <- test1[which(test1$id==(i-6) & test1$visittime<5),]
      }else if((test$status[i-6]==1 & test$time[i-6]==2)){
        vbl1 <- c(i-6,i-6,1,2,test$dose[i-6])
        mbl1 <- test1[which(test1$id==(i-6) & test1$visittime<2),]
      }else if((test$status[i-6]==1 & test$time[i-6]==3)){
        vbl1 <- c(i-6,i-6,1,3,test$dose[i-6])
        mbl1 <- test1[which(test1$id==(i-6) & test1$visittime<3),]
      }else if((test$status[i-6]==1 & test$time[i-6]==4)){
        vbl1 <- c(i-6,i-6,1,4,test$dose[i-6])
        mbl1 <- test1[which(test1$id==(i-6) & test1$visittime<4),]
      }else{
        vbl1 <- c(i-6,i-6,1,1,test$dose[i-6])
        mbl1 <- test1[which(test1$id==(i-6) & test1$visittime<1),]
      }
    }
  }
  
  if(isTRUE(i-count==4)){
    if(isTRUE(test$id[i]==test$counter[i])){ 
      if((test$status[i]==0 & test$time[i]==3) | (test$status[i]==0 & test$time[i]==2) |(test$status[i]==1 & test$time[i]==3) |
         (test$status[i]==1 & test$time[i]==2) | (test$status[i]==1 & test$time[i]==1) | (test$status[i]==0 & test$time[i]==4) |
         (test$status[i]==0 & test$time[i]==5) | (test$status[i]==0 & test$time[i]==6)| (test$status[i]==1 & test$time[i]==4) | 
         (test$status[i]==1 & test$time[i]==5) | (test$status[i]==1 & test$time[i]==6) ){
        vb<- c(i,i,1,1,test$dose[i])
        mb <- test1[which(test1$id==i & test1$visittime<1),]
      }else{
        vb <- c(i,i,0,1,test$dose[i])
        mb <- test1[which(test1$id==i & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-1])==1)){
      if((test$status[i-1]==0 & test$time[i-1]==3) | (test$status[i-1]==0 & test$time[i-1]==4) |
         (test$status[i-1]==0 & test$time[i-1]==5) | (test$status[i-1]==0 & test$time[i-1]==6) ){
        vbt<- c(i-1,i-1,1,2,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<2),]
      }else if ((test$status[i-1]==0 & test$time[i-1]==2)){
        vbt <- c(i-1,i-1,0,2,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<2),]
      }else if((test$status[i-1]==0 & test$time[i-1]==1)){
        vbt <- c(i-1,i-1,0,1,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<1),]
      }else if((test$status[i-1]==1 & test$time[i-1]==2) | (test$status[i-1]==1 & test$time[i-1]==3) |
               (test$status[i-1]==1 & test$time[i-1]==4) | (test$status[i-1]==1 & test$time[i-1]==5) |
               (test$status[i-1]==1 & test$time[i-1]==6)){
        vbt <- c(i-1,i-1,1,2,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<2),]
      }else{
        vbt <- c(i-1,i-1,1,1,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<1),]
      }
    } 
    
    if(isTRUE((test$counter[i]-test$id[i-2])==2)){
      if((test$status[i-2]==0 & test$time[i-2]==3) ){
        vbc<- c(i-2,i-2,0,3,test$dose[i-2])
        mbc <- test1[which(test1$id==(i-2) & test1$visittime<3),]
      }else if ((test$status[i-2]==0 & test$time[i-2]==2)){
        vbc <- c(i-2,i-2,0,2,test$dose[i-2])
        mbc <- test1[which(test1$id==(i-2) & test1$visittime<2),]
      }else if((test$status[i-2]==0 & test$time[i-2]==1)){
        vbc <- c(i-2,i-2,0,1,test$dose[i-2])
        mbc <- test1[which(test1$id==(i-2) & test1$visittime<1),]
      }else if((test$status[i-2]==0 & test$time[i-2]==4) | (test$status[i-2]==0 & test$time[i-2]==5) |
               (test$status[i-2]==0 & test$time[i-2]==6) | (test$status[i-2]==1 & test$time[i-2]==3) |
               (test$status[i-2]==1 & test$time[i-2]==4) | (test$status[i-2]==1 & test$time[i-2]==5) |
               (test$status[i-2]==1 & test$time[i-2]==6)){
        vbc <- c(i-2,i-2,1,3,test$dose[i-2])
        mbc <- test1[which(test1$id==(i-2) & test1$visittime<3),]
      }else if((test$status[i-2]==1 & test$time[i-2]==2)){
        vbc <- c(i-2,i-2,1,2,test$dose[i-2])
        mbc <- test1[which(test1$id==(i-2) & test1$visittime<2),]
      }else{
        vbc <- c(i-2,i-2,1,1,test$dose[i-2])
        mbc <- test1[which(test1$id==(i-2) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-3])==3)){
      if((test$status[i-3]==0 & test$time[i-3]==4) ){
        vbk <- c(i-3,i-3,0,4,test$dose[i-3])
        mbk <- test1[which(test1$id==(i-3) & test1$visittime<4),]
      }else if((test$status[i-3]==0 & test$time[i-3]==3) ){
        vbk <- c(i-3,i-3,0,3,test$dose[i-3])
        mbk <- test1[which(test1$id==(i-3) & test1$visittime<3),]
      }else if ((test$status[i-3]==0 & test$time[i-3]==2)){
        vbk <- c(i-3,i-3,0,2,test$dose[i-3])
        mbk <- test1[which(test1$id==(i-3) & test1$visittime<2),]
      }else if((test$status[i-3]==0 & test$time[i-3]==1)){
        vbk <- c(i-3,i-3,0,1,test$dose[i-3])
        mbk <- test1[which(test1$id==(i-3) & test1$visittime<1),]
      }else if((test$status[i-3]==0 & test$time[i-3]==5) |
               (test$status[i-3]==0 & test$time[i-3]==6) |
               (test$status[i-3]==1 & test$time[i-3]==4) | (test$status[i-3]==1 & test$time[i-3]==5) |
               (test$status[i-3]==1 & test$time[i-3]==6)){
        vbk <- c(i-3,i-3,1,4,test$dose[i-3])
        mbk <- test1[which(test1$id==(i-3) & test1$visittime<4),]
      }else if((test$status[i-3]==1 & test$time[i-3]==2)){
        vbk <- c(i-3,i-3,1,2,test$dose[i-3])
        mbk <- test1[which(test1$id==(i-3) & test1$visittime<2),]
      }else if((test$status[i-3]==1 & test$time[i-3]==3)){
        vbk <- c(i-3,i-3,1,3,test$dose[i-3])
        mbk <- test1[which(test1$id==(i-3) & test1$visittime<3),]
      }else{
        vbk <- c(i-3,i-3,1,1,test$dose[i-3])
        mbk <- test1[which(test1$id==(i-3) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-4])==4)){
      if((test$status[i-4]==0 & test$time[i-4]==5) ){
        vbk1 <- c(i-4,i-4,0,5,test$dose[i-4])
        mbk1 <- test1[which(test1$id==(i-4) & test1$visittime<5),]
      }else if((test$status[i-4]==0 & test$time[i-4]==4) ){
        vbk1 <- c(i-4,i-4,0,4,test$dose[i-4])
        mbk1 <- test1[which(test1$id==(i-4) & test1$visittime<4),]
      }else if((test$status[i-4]==0 & test$time[i-4]==3) ){
        vbk1 <- c(i-4,i-4,0,3,test$dose[i-4])
        mbk1 <- test1[which(test1$id==(i-4) & test1$visittime<3),]
      }else if ((test$status[i-4]==0 & test$time[i-4]==2)){
        vbk1 <- c(i-4,i-4,0,2,test$dose[i-4])
        mbk1 <- test1[which(test1$id==(i-4) & test1$visittime<2),]
      }else if((test$status[i-4]==0 & test$time[i-4]==1)){
        vbk1 <- c(i-4,i-4,0,1,test$dose[i-4])
        mbk1 <- test1[which(test1$id==(i-4) & test1$visittime<1),]
      }else if((test$status[i-4]==0 & test$time[i-4]==6) | (test$status[i-4]==1 & test$time[i-4]==5) |
               (test$status[i-4]==1 & test$time[i-4]==6)){
        vbk1 <- c(i-4,i-4,1,5,test$dose[i-4])
        mbk1 <- test1[which(test1$id==(i-4) & test1$visittime<5),]
      }else if((test$status[i-4]==1 & test$time[i-4]==2)){
        vbk1 <- c(i-4,i-4,1,2,test$dose[i-4])
        mbk1 <- test1[which(test1$id==(i-4) & test1$visittime<2),]
      }else if((test$status[i-4]==1 & test$time[i-4]==3)){
        vbk1 <- c(i-4,i-4,1,3,test$dose[i-4])
        mbk1 <- test1[which(test1$id==(i-4) & test1$visittime<3),]
      }else if((test$status[i-4]==1 & test$time[i-4]==4)){
        vbk1 <- c(i-4,i-4,1,4,test$dose[i-4])
        mbk1 <- test1[which(test1$id==(i-4) & test1$visittime<4),]
      }else{
        vbk1 <- c(i-4,i-4,1,1,test$dose[i-4])
        mbk1 <- test1[which(test1$id==(i-4) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-5])==5)){
      if((test$status[i-5]==0 & test$time[i-5]==5) ){
        vbl1 <- c(i-5,i-5,0,5,test$dose[i-5])
        mbl1 <- test1[which(test1$id==(i-5) & test1$visittime<5),]
      }else if((test$status[i-5]==0 & test$time[i-5]==4) ){
        vbl1 <- c(i-5,i-5,0,4,test$dose[i-5])
        mbl1 <- test1[which(test1$id==(i-5) & test1$visittime<4),]
      }else if((test$status[i-5]==0 & test$time[i-5]==3) ){
        vbl1 <- c(i-5,i-5,0,3,test$dose[i-5])
        mbl1 <- test1[which(test1$id==(i-5) & test1$visittime<3),]
      }else if ((test$status[i-5]==0 & test$time[i-5]==2)){
        vbl1 <- c(i-5,i-5,0,2,test$dose[i-5])
        mbl1 <- test1[which(test1$id==(i-5) & test1$visittime<2),]
      }else if((test$status[i-5]==0 & test$time[i-5]==1)){
        vbl1 <- c(i-5,i-5,0,1,test$dose[i-5])
        mbl1 <- test1[which(test1$id==(i-5) & test1$visittime<1),]
      }else if((test$status[i-5]==0 & test$time[i-5]==6) | (test$status[i-5]==1 & test$time[i-5]==5) |
               (test$status[i-5]==1 & test$time[i-5]==6)){
        vbl1 <- c(i-5,i-5,1,5,test$dose[i-5])
        mbl1 <- test1[which(test1$id==(i-5) & test1$visittime<5),]
      }else if((test$status[i-5]==1 & test$time[i-5]==2)){
        vbl1 <- c(i-5,i-5,1,2,test$dose[i-5])
        mbl1 <- test1[which(test1$id==(i-5) & test1$visittime<2),]
      }else if((test$status[i-5]==1 & test$time[i-5]==3)){
        vbl1 <- c(i-5,i-5,1,3,test$dose[i-5])
        mbl1 <- test1[which(test1$id==(i-5) & test1$visittime<3),]
      }else if((test$status[i-5]==1 & test$time[i-5]==4)){
        vbl1 <- c(i-5,i-5,1,4,test$dose[i-5])
        mbl1 <- test1[which(test1$id==(i-5) & test1$visittime<4),]
      }else{
        vbl1 <- c(i-5,i-5,1,1,test$dose[i-5])
        mbl1 <- test1[which(test1$id==(i-5) & test1$visittime<1),]
      }
    }
  }
  
  if(isTRUE(i-count>=5)){
    if(isTRUE(test$id[i]==test$counter[i])){ 
      if((test$status[i]==0 & test$time[i]==3) | (test$status[i]==0 & test$time[i]==2) |(test$status[i]==1 & test$time[i]==3) |
         (test$status[i]==1 & test$time[i]==2) | (test$status[i]==1 & test$time[i]==1) | (test$status[i]==0 & test$time[i]==4) |
         (test$status[i]==0 & test$time[i]==5) | (test$status[i]==0 & test$time[i]==6)| (test$status[i]==1 & test$time[i]==4) | 
         (test$status[i]==1 & test$time[i]==5) | (test$status[i]==1 & test$time[i]==6) ){
        vb<- c(i,i,1,1,test$dose[i])
        mb <- test1[which(test1$id==i & test1$visittime<1),]
      }else{
        vb <- c(i,i,0,1,test$dose[i])
        mb <- test1[which(test1$id==i & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-1])==1)){
      if((test$status[i-1]==0 & test$time[i-1]==3) | (test$status[i-1]==0 & test$time[i-1]==4) |
         (test$status[i-1]==0 & test$time[i-1]==5) | (test$status[i-1]==0 & test$time[i-1]==6) ){
        vbt<- c(i-1,i-1,1,2,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<2),]
      }else if ((test$status[i-1]==0 & test$time[i-1]==2)){
        vbt <- c(i-1,i-1,0,2,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<2),]
      }else if((test$status[i-1]==0 & test$time[i-1]==1)){
        vbt <- c(i-1,i-1,0,1,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<1),]
      }else if((test$status[i-1]==1 & test$time[i-1]==2) | (test$status[i-1]==1 & test$time[i-1]==3) |
               (test$status[i-1]==1 & test$time[i-1]==4) | (test$status[i-1]==1 & test$time[i-1]==5) |
               (test$status[i-1]==1 & test$time[i-1]==6)){
        vbt <- c(i-1,i-1,1,2,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<2),]
      }else{
        vbt <- c(i-1,i-1,1,1,test$dose[i-1])
        mbt <- test1[which(test1$id==(i-1) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-2])==2)){
      if((test$status[i-2]==0 & test$time[i-2]==3) ){
        vbc<- c(i-2,i-2,0,3,test$dose[i-2])
        mbc <- test1[which(test1$id==(i-2) & test1$visittime<3),]
      }else if ((test$status[i-2]==0 & test$time[i-2]==2)){
        vbc <- c(i-2,i-2,0,2,test$dose[i-2])
        mbc <- test1[which(test1$id==(i-2) & test1$visittime<2),]
      }else if((test$status[i-2]==0 & test$time[i-2]==1)){
        vbc <- c(i-2,i-2,0,1,test$dose[i-2])
        mbc <- test1[which(test1$id==(i-2) & test1$visittime<1),]
      }else if((test$status[i-2]==0 & test$time[i-2]==4) | (test$status[i-2]==0 & test$time[i-2]==5) |
               (test$status[i-2]==0 & test$time[i-2]==6) | (test$status[i-2]==1 & test$time[i-2]==3) |
               (test$status[i-2]==1 & test$time[i-2]==4) | (test$status[i-2]==1 & test$time[i-2]==5) |
               (test$status[i-2]==1 & test$time[i-2]==6)){
        vbc <- c(i-2,i-2,1,3,test$dose[i-2])
        mbc <- test1[which(test1$id==(i-2) & test1$visittime<3),]
      }else if((test$status[i-2]==1 & test$time[i-2]==2)){
        vbc <- c(i-2,i-2,1,2,test$dose[i-2])
        mbc <- test1[which(test1$id==(i-2) & test1$visittime<2),]
      }else{
        vbc <- c(i-2,i-2,1,1,test$dose[i-2])
        mbc <- test1[which(test1$id==(i-2) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-3])==3)){
      if((test$status[i-3]==0 & test$time[i-3]==4) ){
        vbk <- c(i-3,i-3,0,4,test$dose[i-3])
        mbk <- test1[which(test1$id==(i-3) & test1$visittime<4),]
      }else if((test$status[i-3]==0 & test$time[i-3]==3) ){
        vbk <- c(i-3,i-3,0,3,test$dose[i-3])
        mbk <- test1[which(test1$id==(i-3) & test1$visittime<3),]
      }else if ((test$status[i-3]==0 & test$time[i-3]==2)){
        vbk <- c(i-3,i-3,0,2,test$dose[i-3])
        mbk <- test1[which(test1$id==(i-3) & test1$visittime<2),]
      }else if((test$status[i-3]==0 & test$time[i-3]==1)){
        vbk <- c(i-3,i-3,0,1,test$dose[i-3])
        mbk <- test1[which(test1$id==(i-3) & test1$visittime<1),]
      }else if((test$status[i-3]==0 & test$time[i-3]==5) |
               (test$status[i-3]==0 & test$time[i-3]==6) |
               (test$status[i-3]==1 & test$time[i-3]==4) | (test$status[i-3]==1 & test$time[i-3]==5) |
               (test$status[i-3]==1 & test$time[i-3]==6)){
        vbk <- c(i-3,i-3,1,4,test$dose[i-3])
        mbk <- test1[which(test1$id==(i-3) & test1$visittime<4),]
      }else if((test$status[i-3]==1 & test$time[i-3]==2)){
        vbk <- c(i-3,i-3,1,2,test$dose[i-3])
        mbk <- test1[which(test1$id==(i-3) & test1$visittime<2),]
      }else if((test$status[i-3]==1 & test$time[i-3]==3)){
        vbk <- c(i-3,i-3,1,3,test$dose[i-3])
        mbk <- test1[which(test1$id==(i-3) & test1$visittime<3),]
      }else{
        vbk <- c(i-3,i-3,1,1,test$dose[i-3])
        mbk <- test1[which(test1$id==(i-3) & test1$visittime<1),]
      }
    }
    
    if(isTRUE((test$counter[i]-test$id[i-4])==4)){
      if((test$status[i-4]==0 & test$time[i-4]==5) ){
        vbl <- c(i-4,i-4,0,5,test$dose[i-4])
        mbl <- test1[which(test1$id==(i-4) & test1$visittime<5),]
      }else if((test$status[i-4]==0 & test$time[i-4]==4) ){
        vbl <- c(i-4,i-4,0,4,test$dose[i-4])
        mbl <- test1[which(test1$id==(i-4) & test1$visittime<4),]
      }else if((test$status[i-4]==0 & test$time[i-4]==3) ){
        vbl <- c(i-4,i-4,0,3,test$dose[i-4])
        mbl <- test1[which(test1$id==(i-4) & test1$visittime<3),]
      }else if ((test$status[i-4]==0 & test$time[i-4]==2)){
        vbl <- c(i-4,i-4,0,2,test$dose[i-4])
        mbl <- test1[which(test1$id==(i-4) & test1$visittime<2),]
      }else if((test$status[i-4]==0 & test$time[i-4]==1)){
        vbl <- c(i-4,i-4,0,1,test$dose[i-4])
        mbl <- test1[which(test1$id==(i-4) & test1$visittime<1),]
      }else if((test$status[i-4]==0 & test$time[i-4]==6) | (test$status[i-4]==1 & test$time[i-4]==5) |
               (test$status[i-4]==1 & test$time[i-4]==6)){
        vbl <- c(i-4,i-4,1,5,test$dose[i-4])
        mbl <- test1[which(test1$id==(i-4) & test1$visittime<5),]
      }else if((test$status[i-4]==1 & test$time[i-4]==2)){
        vbl <- c(i-4,i-4,1,2,test$dose[i-4])
        mbl <- test1[which(test1$id==(i-4) & test1$visittime<2),]
      }else if((test$status[i-4]==1 & test$time[i-4]==3)){
        vbl <- c(i-4,i-4,1,3,test$dose[i-4])
        mbl <- test1[which(test1$id==(i-4) & test1$visittime<3),]
      }else if((test$status[i-4]==1 & test$time[i-4]==4)){
        vbl <- c(i-4,i-4,1,4,test$dose[i-4])
        mbl <- test1[which(test1$id==(i-4) & test1$visittime<4),]
      }else{
        vbl <- c(i-4,i-4,1,1,test$dose[i-4])
        mbl <- test1[which(test1$id==(i-4) & test1$visittime<1),]
      }
    }  
  }
  
  ver <- c(!is.null(vb[1]), !is.null(vbt[1]), !is.null(vbt1[1]), !is.null(vbc[1]), !is.null(vbc1[1]), !is.null(vbk[1]), 
           !is.null(vbk1[1]), !is.null(vbl[1]),!is.null(vbl1[1]))
  
  testtemp <- test[,]
  testtemp <- test[which(test$id<=(i-(sum(ver, na.rm=TRUE)))),] 
  testtemp <- rbind(testtemp, vbl1) 
  testtemp <- rbind(testtemp, vbl) 
  testtemp <- rbind(testtemp, vbk1) 
  testtemp <- rbind(testtemp, vbk) 
  testtemp <- rbind(testtemp, vbc1) 
  testtemp <- rbind(testtemp, vbc) 
  testtemp <- rbind(testtemp, vbt1) 
  testtemp <- rbind(testtemp, vbt) 
  testtemp <- rbind(testtemp, vb) 
  
  colnames(testtemp)[1] <- "id"
  colnames(testtemp)[2] <- "counter"
  colnames(testtemp)[3] <- "status"
  colnames(testtemp)[4] <- "time"
  colnames(testtemp)[5] <- "dose"
  
  
  test1temp <- test1[,]
  test1temp <- test1[which(test1$id<=(i-(sum(ver, na.rm=TRUE)))),] 
  test1temp <- rbind(test1temp, mbl1) 
  test1temp <- rbind(test1temp, mbl) 
  test1temp <- rbind(test1temp, mbk1) 
  test1temp <- rbind(test1temp, mbk) 
  test1temp <- rbind(test1temp, mbc1) 
  test1temp <- rbind(test1temp, mbc) 
  test1temp <- rbind(test1temp, mbt1) 
  test1temp <- rbind(test1temp, mbt) 
  test1temp <- rbind(test1temp, mb) 
  
  out <- list(survtemp=testtemp, longtemp=test1temp) 
}




#############################################################
#####      CALCULATES THE 2+2 - CRM PROCESS             #####

crm <- function(meanpars, varpars, gampars, nt, n1, d, dose, nvisit, maxiter, target, censor, pl, stopping_rule, inter_miss, fix_param){
  # Runs the 2+2 and CRM design
  # Stop_moment takes value 0 if design stops at 2+2, 1 for CRM and 2 if it continues till the end 
  # If stop_moment is 0 or 1 then toxic should be 1
  
  temp <- escalate(meanpars, varpars, gampars, nt, n1, d, nvisit, dose, censor, pl, stopping_rule, inter_miss) 
  survdata <- temp$survdata
  longdata <- temp$longdata
  survtemp <- temp$survtemp
  longtemp <- temp$longtemp
  count1<- nrow(survdata)
  toxic <- temp$toxic
  lastestimacrm <- 0 
  lastsecrm <- 0
  stop_moment <- 0 
  if(toxic!=1){
    
    myprobit <- glm(status ~ time + offset(fix_param*dose), family=binomial(link="probit"), control = list(maxit = 100), data=survdata) 
    inits <- c(myprobit$coefficients[1], myprobit$coefficients[2])
    fit <- nlminb(start = inits, objective=loglikcrm, surv=survtemp, pmnormtol=1e-3, fix_param=fix_param, control=list(iter.max=maxiter, eval.max=maxiter))
    fithess <- hessian(func =loglikcrm, x=fit$par, surv=survtemp, pmnormtol=1e-3, fix_param=fix_param)
    foo <- function(X) {
      RES <- try(solve(X,tol=1e-08))
      if (class(RES) == 'try-error'){
        return(rep(NA, nrow(X)))
      } else {
        return(sqrt(diag(RES)))
      }
    }
    se <- foo(fithess)
    estim <- fit$par
    
    # Sequential allocation of subjects to a dose level
    for(i in (count1+1):n1){
      
      temp1 <- nextdose(estim, nt, n1, d, dose, target, fix_param) 
      
      # Condition for not skipping dose levels
      previousdose <- which(dose==survdata$dose[i-1]) 
      if(temp1>(previousdose + 1)){ 
        temp1 <- previousdose + 1
      }else if(temp1<(previousdose - 1)){
        temp1 <- previousdose - 1
      }else{
        temp1 <- temp1
      }
      
      temp2 <- temp$rawsurv                        
      temp3 <- temp2[temp1, , i]                   
      
      if(sum(temp3)==nt){   
        tm <- nt
      } else{
        tm <- min(which(temp3==0))
      } 
      status <- min(temp3)	
      
      if(censor){ 
        chose <- rbinom(1, 1, 0.2)
        censortime <- runif(1,0,nt)
        censortime <- ceiling(censortime)
        status <- ifelse((censortime < tm) & (chose==1),1,status)
        tm <- ifelse((censortime < tm) & (chose==1),censortime,tm)
      }
      
      temp1 <- dose[temp1]                
      addin <- c(i, i, status, tm, temp1)  
      survdata = rbind(survdata,addin)        
      survdata <- data2plcrmjm(temp, survdata, longdata, nvisit, dose, inter_miss)$surv 
      longdata <- data2plcrmjm(temp, survdata, longdata, nvisit, dose, inter_miss)$dat  
      survtemp <- sequentialentry(survdata, longdata, i, count1)$survtemp
      longtemp <- sequentialentry(survdata, longdata, i, count1)$longtemp
      # Counts the DLTs of the first dose and of the first cycle
      count_DLTs <- NROW(which((survtemp$status==0) & (survtemp$time==1) & (survtemp$dose==1.2)))
      if(count_DLTs >= stopping_rule[i]){
        stop_moment <- 1 # stop during crm
        toxic <- 1 
        se <- se
        estim <- estim  
        break
      }
      fit <- nlminb(start = estim, objective=loglikcrm, surv=survtemp, pmnormtol=1e-3, fix_param=fix_param, control=list(iter.max=maxiter, eval.max=maxiter))
      fithess <- hessian(func=loglikcrm, x=fit$par, surv=survtemp, pmnormtol=1e-3, fix_param=fix_param)
      foo <- function(X) {
        RES <- try(solve(X,tol=1e-08))
        if (class(RES) == 'try-error'){
          return(rep(NA, nrow(X)))
        } else {
          return(sqrt(diag(RES)))
        }
      }
      se <- foo(fithess)
      estim <- fit$par  
      lastestimacrm <- estim
      lastsecrm <- se
      toxic <- 0
      stop_moment <- 2 # no stopping for excessive toxicity
    }
  } 
  out<- list(survtemp= survtemp, survdata=survdata, longtemp=longtemp, longdata=longdata, rawsurv=temp$rawsurv, rawlong=temp$rawlong, 
             count1=count1, lastestimacrm=lastestimacrm, lastsecrm=lastsecrm, toxic=toxic, stop_moment=stop_moment)
  
}




#############################################################
#####         CALCULATES THE JM LOGLIKELIHOOD           #####

loglik <- function(pars, ses=NULL, dat, surv, dose, pmnormtol=1e-3, printpars=F, ds, fix_param){
  # Calculates the log likelihood of the JM 
  # pmnormtol specifies the accuracy to which multivariate normal CDF's should be calculated
  # printpars to print parameters for all calls to loglik evaluated by nlminb
  
  n <- dim(surv)[1]
  ny <- dim(dat)[1] 
  nt <- max(surv$time) # max number of time intervals
  
  
  py <- dim(dat)[2]  - 2 
  pt <- dim(surv)[2] - 2 
  pyt <- py+pt
  by <- pars[1:(py+1)]  # point estimates of the longitudinal data
  bt <- pars[(py+2):pyt] # point estimates of the survival data
  
  sig.y <- exp(pars[(pyt+1)]) # residual variance
  
  
  upars <- c(exp(pars[(pyt+2)])) 
  gampars <- pars[(pyt+3)]  
  
  
  Sig.u <- Sigma.is(upars)  
  L <- L.is(p=nt, gampars)   
  nu <- 1   
  
  
  id <- dat$id
  y <- dat$y
  time.y <- dat$visittime 
  if(!is.integer(time.y)){ 
    timeinterval.y <- as.integer(floor(time.y)+1)
  } else {
    timeinterval.y <- time.y	
  }
  X.y <- dat[,-(1:3)] 
  status <- surv$status
  time <- surv$time  
  X.t <- surv[,-(1:4)] 
  
  # calculate the matrix (vector here) R = A^T A / \v^2
  
  R <- array(0,dim=c(1,1,n))
  R[1,1,] <- tapply(time.y^2,factor(id,levels=surv$id),sum) 
  R <- ifelse(is.na(R),0,R)
  R <- R / sig.y^2  
  
  
  # Calculate V = H^{-1} = (A^T A / \nu^2 + \Sigma^{-1})^{-1}
  
  Sig.u.inv <- 1/Sig.u  
  getV <- function(M){solve(M+Sig.u.inv)} 
  V <- array(0,dim=c(1,1,n))
  V[1,1,] <- aaply(R,3,getV) 
  
  
  if(ds!=1){
    X.y1<- ifelse((X.y<dose[ds]),X.y, dose[ds] ) # ds is the candidate dose where the function begins to plateau
    
    # remove fixed effects from y intercept, time and time dose interaction
    yran <- y - by[1] - by[2]*(time.y)^2 - by[3]*time.y*as.matrix(X.y1) 
  }
  
  if(ds==1){
    X.y1<- X.y 
    
    yran <- y - by[1] - by[2]*(time.y)^2 - by[3]*time.y*as.matrix(X.y1) 
  }
  
  # calculate r = A^T (y - X *beta) / \v^2
  r <- matrix(nrow=n,ncol=1)
  r[,1] <- tapply(yran*time.y,factor(id,levels=surv$id),sum) /sig.y^2 
  r <- ifelse(is.na(r),0,r)
  
  
  
  # calculate q.Phi, arguments of multivariate normal CDFs
  Vr <- matrix(nrow=nu,ncol=n)
  for(j in 1:nu){
    Vr[j,] <- V[j,,] * t(r) 
  }
  cov.const <- bt[1] + as.matrix(X.t) %*% as.matrix(fix_param) 
  betaX <- matrix(cov.const,nrow=nt,ncol=n,byrow=T) + bt[2]*matrix(0:(nt-1),nrow=nt,ncol=n)
  q.Phi <- betaX + L%*% Vr
  
  # calculate the multivariate normal CDFs and put results in Phi,
  # Phi[,i] is the probability of surviving to the start of the ith interval
  Phi <- matrix(0,nrow=n,ncol=(nt+1))  
  # First column of Phi is the probability of surviving to time 0
  Phi[,1] <- 1
  I <- array(diag(nt),dim=c(nt,nt,n)) 
  for(i in 1:nt){
    Inew <- I[1:i,1:i,]
    
    Lnew <- L[1:i,1]
    if(i==1){
      tmp <- tensor(as.matrix(Lnew),V,1,1) 
      LVL <- tensor(tmp,as.matrix(Lnew),1,1) 
    } else {
      tmp <- tensor(as.matrix(Lnew),V,2,1)
      LVL <- tensor(as.matrix(Lnew),tmp,2,2)
    }
    
    
    Sig.Phi <- Inew + LVL
    mu.Phi <- rep(0,i)
    
    # only calculate the components of Phi that are actually going to be used
    # (i.e. time for censored individuals and time and (time - 1) for deaths)
    ids1 <- c(which((time==i)&status==0),which((time==(i+1))&status==0))
    ids2 <- which((time==i)&status==1)
    ids <- c(ids1,ids2)
    
    p.Phi <- rep(0,n)
    if(i==1){
      p.Phi[ids] <- pnorm(q.Phi[1,ids],mean=0,sd=sqrt(Sig.Phi[ids]))
    } else {
      for(k in ids){
        p.Phi[k] <- pmnorm(q.Phi[1:i,k],mu.Phi,Sig.Phi[,,k],abseps=pmnormtol)
      }
    }
    Phi[,(i+1)] <- p.Phi
  }
  
  # Psi[,i] is the probability of dying in time interval i
  Psi <- matrix(nrow=n,ncol=nt)
  for(i in 1:nt){
    Psi[,i] <- Phi[,i] - Phi[,(i+1)]
  }
  
  # calculate all the terms in the likelihood
  term1 <- 0
  for(i in 1:nt){
    ids <- which((time==i)&status==0)
    term1 <- term1 + sum(log(Psi[ids,i]))
  }
  
  term2 <- 0
  for(i in 1:nt){
    ids <- which((time==i)&status==1)
    term2 <- term2 + sum(log(Phi[ids,(i+1)]))
  }
  term3 <- - sum(yran^2)/(2*sig.y^2)
  tmp <- matrix(nrow=nu,ncol=nu)
  for(i in 1:nu){
    for(j in 1:nu){
      tmp[i,j] <- sum(r[,i] * V[i,j,] * r[,j])
    }	
  }
  term4 <- 0.5 * sum(tmp)
  term5 <- -ny*log(sig.y)
  term6 <- - 0.5*n*log(Sig.u)
  term7 <- 0.5*sum(apply(V,3,log))  
  
  loglik <- term1 + term2 + term3 + term4 + term5 + term6 + term7
  # Return minus loglik - optimisation routines will minimise this
  if (is.finite(loglik)) {
    return(-loglik)
  } else {
    return(1e+09)
  }
  
}





#######################################################################
#####                    MIN OBS TO FIT THE JM                    #####

checkconv <- function(dose, count1, n_min_jm, n1, kl, maxiter, se_max_jm, fix_param){
  # This function checks after how many observations the algorithm can proceed to the JM
  
  count2 <- NULL
  # Safety constraint. JM cannot start before a certain number of subjects is enrolled (n_min_jm)
  if(count1< n_min_jm){
    count2 <- n_min_jm
  }else{
    count2 <- count1
  }
  
  for(i in count2:round((5/6)*n1)){
    print(i)
    survdata <- kl$survdata[1:i,]
    longdata <- kl$longdata
    longdata <- kl$longdata[which(longdata$id<=i),]
    
    
    if(i==count1){
      survtemp <- sequentialentry0(survdata, longdata, (i-1))$survtemp
      longtemp <- sequentialentry0(survdata, longdata, (i-1))$longtemp
    }else{
      survtemp <- sequentialentry(survdata, longdata, i, count1)$survtemp
      longtemp <- sequentialentry(survdata, longdata, i, count1)$longtemp
    }
    
    count2 <- nrow(survdata)
    
    # Initial values for the probit model
    myprobit <- glm(status ~ time + offset(fix_param*dose), family=binomial(link="probit"), control = list(maxit = 100), data=survtemp) 
    probinits <- c(myprobit$coefficients[1], myprobit$coefficients[2])
    
    
    count_doses <- NROW(unique(survdata$dose)) 
    if(count_doses==1){
      aic <- NA
      fit1 <- matrix(NA, nrow = 1, ncol = 8)
      lmedat=try(lme(y~  I(visittime*visittime) + I(dose*visittime), 
                     random=~ visittime| id, control =list(maxiter=maxiter, msMaxIter=maxiter, msMaxEval=maxiter, opt="optim"), data=longtemp))
      if(class(lmedat) == 'try-error'){
        next 
      }
      coef <- fixef(lmedat)
      sigma <- lmedat$sigma  
      re= try(intervals(lmedat)$reStruct$id[2,2]) 
      if(class(re) == 'try-error'){
        re <- runif(1,1,10) 
      }
      # Initial values for the gamma parameter
      gamma_inits <- round(runif(1,0,1),2)
      inits <- c(coef[1:3], probinits,  log(sigma), log(re), gamma_inits)
      fit <- nlminb(start = inits, objective=loglik, dat=longtemp, surv=survtemp, dose=dose, pmnormtol=1e-3, ds=count_doses, fix_param=fix_param, control=list(iter.max=maxiter, eval.max=maxiter))
      fithess <- hessian(func=loglik, x=fit$par, dat=longtemp, surv=survtemp, dose=dose, pmnormtol=1e-3, ds=count_doses, fix_param=fix_param)
      RES <- try(solve(fithess,tol=1e-08))
      # Constraint to proceed to JM
      if(class(RES) != 'try-error'){ 
        if(sum(is.nan(sqrt(diag(RES))))==0 & max(sqrt(diag(RES))) < se_max_jm ){
          aic <- (2*length(fit$par))-(2*(-fit$objective))
          fit1[1,] <- fit$par
        }
      }
    }
    
    if(count_doses>=2){
      aic <- rep(NA, (count_doses-1))
      fit1 <- matrix(NA, nrow = (count_doses-1), ncol = 8)
      for(j in 1:(count_doses-1)){
        if(j==1){
          lmedat=try(lme(y~  I(visittime*visittime) + I(dose*visittime), 
                         random=~ visittime| id, control =list(maxiter=maxiter, msMaxIter=maxiter, msMaxEval=maxiter, opt="optim"), data=longtemp))
          if(class(lmedat) == 'try-error'){
            next 
          }
          
        }else{
          dose_temp1 <- ifelse(longtemp$dose<dose[j], longtemp$dose, dose[j])
          long_temp <- longtemp
          long_temp$dose_temp1 <- dose_temp1
          
          
          lmedat=try(lme(y~  I(visittime*visittime) + I(dose_temp1*visittime), 
                         random=~ visittime| id, control =list(maxiter=maxiter, msMaxIter=maxiter, msMaxEval=maxiter, opt="optim"), data=long_temp))
          if(class(lmedat) == 'try-error'){
            next 
          }
        }
        coef <- fixef(lmedat)
        sigma <- lmedat$sigma  
        re= try(intervals(lmedat)$reStruct$id[2,2]) 
        if(class(re) == 'try-error'){
          re <- runif(1,1,10) 
        }
        # Initial values for the gamma parameter
        gamma_inits <- round(runif(1,0,1),2)
        inits <- c(coef[1:3], probinits, log(sigma), log(re), gamma_inits)
        fit <- nlminb(start = inits, objective=loglik, dat=longtemp, surv=survtemp, dose=dose, pmnormtol=1e-3, ds=j, fix_param=fix_param, control=list(iter.max=maxiter, eval.max=maxiter))
        fithess <- hessian(func=loglik, x=fit$par, dat=longtemp, surv=survtemp, dose=dose, pmnormtol=1e-3, ds=j, fix_param=fix_param)
        RES <- try(solve(fithess,tol=1e-08))
        if(class(RES) != 'try-error'){ 
          if(sum(is.nan(sqrt(diag(RES))))==0 & max(sqrt(diag(RES))) < se_max_jm ){
            aic[j] <- (2*length(fit$par))-(2*(-fit$objective))
            fit1[j,] <- fit$par
          }
        }
      }
    }
    
    
    if(any(!is.na(aic))==TRUE){
      break
    }
    
  }  
  
  
  if(exists('fit1')==TRUE){
    festim <- fit1
  }else{
    festim <- NULL
    aic <- NULL
  }
  
  out <- list(count2=count2, festim=festim, aic=aic, survdata=survdata, survtemp=survtemp, longdata=longdata, longtemp=longtemp)  
}



#######################################################################
#####    CALCULATES DOSE ASSINGED TO EVERY SUBJECT ON THE JM     #####

nextdose1 <- function(estima, d, dose, target, nvisit, plateau, meandiff, fix_param){
  # It calculates the dose for the next subject 
  # First it estimates the cumulative probability of toxicity per dose level 
  # Then it calculates the distance between the target toxicity level and the cumulative probability of toxicity at each dose level
  # MTD is the dose that minimizes this distance
  # For each dose below and equal to the MTD it calculates the biomarker activity
  # Select the OD based on the mean difference 
  
  # estima are the parameters estimated from the joint modeling
  
  tox <- rep(NA,d)
  
  for(i in 1:d){
    tox[i] <- ((1-pnorm(estima[4] + fix_param*dose[i]))+ (pnorm(estima[4] + fix_param*dose[i]))*(1-pnorm(estima[4] + fix_param*dose[i] +estima[5]))+ 
                 (pnorm(estima[4] + fix_param*dose[i]))*(pnorm(estima[4] + fix_param*dose[i]+estima[5]))*(1-pnorm(estima[4] + fix_param*dose[i]+ estima[5]*2)) +
                 (pnorm(estima[4] + fix_param*dose[i])*(pnorm(estima[4] + fix_param*dose[i]+estima[5]))*(pnorm(estima[4] + fix_param*dose[i]+estima[5]*2))*(1-pnorm(estima[4] + fix_param*dose[i]+estima[5]*3)))+
                 (pnorm(estima[4] + fix_param*dose[i])*(pnorm(estima[4] + fix_param*dose[i]+estima[5]))*(pnorm(estima[4] + fix_param*dose[i]+estima[5]*2))*
                    (pnorm(estima[4] + fix_param*dose[i]+estima[5]*3))*(1-pnorm(estima[4] + fix_param*dose[i]+estima[5]*4)))+
                 (pnorm(estima[4] + fix_param*dose[i])*(pnorm(estima[4] + fix_param*dose[i]+estima[5]))*(pnorm(estima[4] + fix_param*dose[i]+estima[5]*2))*
                    (pnorm(estima[4] + fix_param*dose[i]+estima[5]*3))*(pnorm(estima[4] + fix_param*dose[i]+estima[5]*4))*(1-pnorm(estima[4] + fix_param*dose[i]+estima[5]*5))))
  }
  dist <-(abs(tox-target)) 
  MTD <- which.min(dist)
  
  visit <- c(0.1, 0.4, 0.8, 1.2, 1.4, 1.8, 2.1, 2.4, 2.8, 3.1, 3.4, 3.8, 4.1, 4.4, 4.8, 5.1, 5.4, 5.8)
  
  
  if(plateau==1){
    eff <- matrix(NA, nrow = MTD, ncol = nvisit)
    for(i in 1:nvisit){
      for(j in 1:MTD){
        eff[j,i] <- estima[1]+ estima[2]*visit[i]^2 + estima[3]*visit[i]*dose[j]
      }
    }
    
  }else{
    dose1<- ifelse((dose<dose[plateau]),dose, dose[plateau]) 
    eff <- matrix(NA, nrow = MTD, ncol = nvisit)
    for(i in 1:nvisit){
      for(j in 1:MTD){
        eff[j,i] <- estima[1]+ estima[2]*visit[i]^2 + estima[3]*visit[i]*dose1[j]
      }
    }
    
  }
  
  maxeff <- rep(NA, MTD)
  maxeff <- apply(eff, 1, min)
  
  
  dist1 <- (abs(maxeff[MTD]- maxeff))
  
  if(max(dist1) <= meandiff){
    OD <- which.min(dose)
  }else{
    OD <- min(which((dist1)>=0 & (dist1)<= meandiff)) 
  }
  
  
  out <- list(MTD=MTD, OD=OD)
  
} 




jmproc <- function(meanpars, varpars, gampars, nt, n1, d ,dose, nvisit, maxiter, target, censor, pl, meandiff, stopping_rule, inter_miss, n_min_jm, se_max_jm, fix_param, seed){
  # This function calls all of the above functions 
  # Runs the entire algorithm from data simulaition to JM and provides the results
  
  # meanpars sets the parameters for the longitudinal and the survival models (in the order longitudinal: intercept, time, dose survival: intercept, cycle, dose)
  # varpars sets the variance parameters (in the order residual variance and random effect variance)
  # gampars sets the gamma parameter 
  # nt sets the number of treatment cycles
  # n1 sets the number of subjects in the trial
  # d sets the total number of doses to be tested
  # dose sets the vector of the dose levels
  # nvisit sets the maximum number of repeated measurements
  # maxiter sets the maximum number of iterations in nlminb
  # target sets the target toxicity level
  # censor sets the uniform censoring (TRUE/FALSE)
  # pl indicates the dose where the plateau begins. In case of no plateau set pl=1
  # meandiff sets the maximum accepted difference so that one or more doses can be considered equivalent to the MTD in terms of activity
  # stopping_rule sets the sequential boundaries to test for excessive toxicity
  # inter_miss sets the probability of intermittent missing visits (per visit)
  # n_min_jm sets the minimum acceptable number of subjects before applying the jm
  # se_max_jm sets the maximum acceptable standard error for parameter estimation in the jm
  # fix_param sets the dose slope of the survival model to a fixed value
  # seed sets the seed to run the simulation study
  
  set.seed(seed)
  # If trial stops for excessive toxicity then function stops here and saves the objects below
  met <- crm(meanpars, varpars, gampars, nt, n1, d, dose, nvisit, maxiter, target, censor, pl, stopping_rule, inter_miss, fix_param) 
  count1 <- met$count1 
  toxic <- met$toxic 
  survdata <- met$survdata
  survtemp <- met$survtemp
  longdata <- met$longdata
  longtemp <- met$longtemp
  stop_moment <- met$stop_moment 
  estima <- 0
  estimase <- 0
  count2 <- 0
  lastdose <- survdata$dose[count1]
  lastMTD <- survdata$dose[count1]
  plateau_overall <- 0
  
  # If trial does not stop for excessive toxicity it proceeds below
  if(toxic!=1){
    # It tests whether it can proceed to the JM
    check1 <- checkconv(dose, count1, n_min_jm, n1, met, maxiter, se_max_jm, fix_param)
    count2 <- check1$count2 
    # If trial proceeds with the CRM exclusively then it saves the objects below
    if(count2==round((5/6)*n1)){ 
      survdata <- met$survdata
      survtemp <- met$survtemp
      longdata <- met$longdata
      longtemp <- met$longtemp
      estima <- met$lastestimacrm
      estimase <- met$lastsecrm
      lastMTD <- survdata$dose[n1]
      lastdose <- survdata$dose[n1]
      plateau_overall <- 0
      count2 <- n1
      # Proceed to the JM
    }else{
      plateau <- which.min(check1$aic) # finds which is the best model 
      estima <- check1$festim[plateau,] # selects the parameters only from the above model
      survdata <- check1$survdata
      longdata <- check1$longdata 
      
      plateau_overall <- rep(NA, n1)
      # Sequential allocation of subjects to a dose level
      for(i in (count2 + 1):n1){
        print(i)
        count_doses <- NROW(unique(survdata$dose)) 
        tempor1 <- nextdose1(estima, d, dose, target, nvisit, plateau, meandiff, fix_param)$OD 
        lastMTD <- nextdose1(estima, d, dose, target, nvisit, plateau, meandiff, fix_param)$MTD
        lastMTD <- dose[lastMTD]

        # Condition for not skipping dose levels
        previousdose <- which(dose==survdata$dose[i-1]) 
        if(tempor1>(previousdose + 1)){ 
          tempor1 <- previousdose + 1
        }else if(tempor1<(previousdose - 1)){
          tempor1 <- previousdose - 1
        }else{
          tempor1 <- tempor1
        }
        
        tempor2 <- met$rawsurv                            
        tempor3 <- tempor2[tempor1, , i]                   
        
        if(sum(tempor3)==nt){      
          tim <- nt
        } else{
          tim <- min(which(tempor3==0))
        } 
        status <- min(tempor3)	
        
        # Uniform censoring
        if(censor){ 
          chose <- rbinom(1, 1, 0.2)
          censortime <- runif(1,0,nt)
          censortime <- ceiling(censortime)
          status <- ifelse((censortime < tim) & (chose==1),1,status)
          tim <- ifelse((censortime < tim) & (chose==1),censortime,tim)
        }
        
        tempor4 <- dose[tempor1]               
        addin <- c(i, i, status, tim, tempor4)  
        survdata = rbind(survdata,addin)        
        survdata <- data2plcrmjm(met, survdata, longdata, nvisit, dose, inter_miss)$surv
        longdata <- data2plcrmjm(met, survdata, longdata, nvisit, dose, inter_miss)$dat   
        survtemp <- sequentialentry(survdata, longdata, i, count1)$survtemp
        longtemp <- sequentialentry(survdata, longdata, i, count1)$longtemp
        
        # Take initial values to run nlminb
        myprobit <- glm(status ~ time + offset(fix_param*dose), family=binomial(link="probit"), control = list(maxit = 100), data=survtemp) 
        probinits <- c(myprobit$coefficients[1], myprobit$coefficients[2])
        
        
        count_doses <- NROW(unique(survdata$dose)) 
        
        
        if(count_doses==1){
          aic <- NA
          fit1 <- matrix(NA, nrow = count_doses, ncol = 8)
          estimase <- matrix(NA, nrow = count_doses, ncol = 8)
          lmedat=try(lme(y~  I(visittime*visittime) + I(dose*visittime), 
                         random=~ visittime| id, control =list(maxiter=maxiter, msMaxIter=maxiter, msMaxEval=maxiter, opt="optim"), data=longtemp))
          if(class(lmedat) == 'try-error'){
            coef <- estima[1:3]
            sigma <- exp(estima[6])
            re <-  exp(estima[7])
          }else{
            coef <- fixef(lmedat)
            sigma <- lmedat$sigma  
            re= try(intervals(lmedat)$reStruct$id[2,2]) 
            if(class(re) == 'try-error'){
              re <- runif(1,1,10) 
            }
          }
          
          # Initial values for the gamma parameter
          gamma_inits <- round(runif(1,0,1),2)
          inits <- c(coef[1:3], probinits,  log(sigma), log(re), gamma_inits)
          fit <- nlminb(start = inits, objective=loglik, dat=longtemp, surv=survtemp, dose=dose, pmnormtol=1e-3, ds=count_doses, fix_param=fix_param, control=list(iter.max=maxiter, eval.max=maxiter))
          fithess <- hessian(func=loglik, x=fit$par, dat=longtemp, surv=survtemp, dose=dose, pmnormtol=1e-3, ds=count_doses, fix_param=fix_param)
          foo <- function(X) {
            RES <- try(solve(X,tol=1e-08))
            if (class(RES) == 'try-error'){
              return(rep(NA, nrow(X)))
            } else {
              return(sqrt(diag(RES)))
            }
          }
          se<-foo(fithess)
          estimase[count_doses,]<-c(se[1:5], (exp(se[6:7])*se[6:7]), se[8])
          aic <- (2*length(fit$par))-(2*(-fit$objective))
          fit1[count_doses,] <- c(fit$par[1:5], fit$par[6:8])
          plateau <- 1
          estima <- fit1
        }
        
        
        if(count_doses>=2){
          aic <- rep(NA, (count_doses-1))
          fit1 <- matrix(NA, nrow = (count_doses-1), ncol = 8)
          estimase <- matrix(NA, nrow = (count_doses-1), ncol = 8)
          
          for(j in 1:(count_doses-1)){
            if(j==1){
              lmedat=try(lme(y~  I(visittime*visittime) + I(dose*visittime), 
                             random=~ visittime| id, control =list(maxiter=maxiter, msMaxIter=maxiter, msMaxEval=maxiter, opt="optim"), data=longtemp))
              if(class(lmedat) == 'try-error'){
                coef <- estima[1:3]
                sigma <- exp(estima[6])
                re <-  exp(estima[7])
              }else{
                coef <- fixef(lmedat)
                sigma <- lmedat$sigma  
                re= try(intervals(lmedat)$reStruct$id[2,2]) 
                if(class(re) == 'try-error'){
                  re <- runif(1,1,10) 
                }
              }
              
            }else{
              dose_temp1 <- ifelse(longtemp$dose<dose[j], longtemp$dose, dose[j])
              long_temp <- longtemp
              long_temp$dose_temp1 <- dose_temp1
              
              lmedat=try(lme(y~  I(visittime*visittime) + I(dose_temp1*visittime), 
                             random=~ visittime| id, control =list(maxiter=maxiter, msMaxIter=maxiter, msMaxEval=maxiter, opt="optim"), data=long_temp))
              if(class(lmedat) == 'try-error'){
                if(is.na(estima[6])){
                  estima[6] <- runif(1,-90,-10)
                }
                coef <- c(estima[1:3])
                sigma <- exp(estima[6])
                re <-  exp(estima[7])
              }else{
                coef <- fixef(lmedat)
                sigma <- lmedat$sigma  
                re= try(intervals(lmedat)$reStruct$id[2,2]) 
                if(class(re) == 'try-error'){
                  re <- runif(1,1,10) 
                }
              }
            }     
            # Initial values for the gamma parameter
            gamma_inits <- round(runif(1,0,1),2)
            inits <- c(coef[1:3], probinits, log(sigma), log(re), gamma_inits)
            fit <- nlminb(start = inits, objective=loglik, dat=longtemp, surv=survtemp, dose=dose, pmnormtol=1e-3, ds=j, fix_param=fix_param, control=list(iter.max=maxiter, eval.max=maxiter))
            fithess <- hessian(func=loglik, x=fit$par, dat=longtemp, surv=survtemp, dose=dose, pmnormtol=1e-3, ds=j, fix_param=fix_param)
            foo <- function(X) {
              RES <- try(solve(X,tol=1e-08))
              if (class(RES) == 'try-error'){
                return(rep(NA, nrow(X)))
              } else {
                return(sqrt(diag(RES)))
              }
            }
            se<-foo(fithess)
            estimase[j,]<-c(se[1:5], (exp(se[6:7])*se[6:7]), se[8])
            aic[j] <- (2*length(fit$par))-(2*(-fit$objective))
            fit1[j,] <- fit$par
          }
          
          plateau <- which.min(aic) # finds which is the best model -> where plateau begins
          estima <- fit1[plateau,]
        }
        
        
        plateau_overall[i] <- plateau
      }
      lastdose <- survdata$dose[n1]
    }
  }  
 
  cr <- rep(NA, d) 
  for(i in 1:d){
    cr[i] <- NROW(which(survdata$dose==dose[i]))
  }
  
  
  doseassignment <-  c(cr, lastdose) 
  lastMTD <- lastMTD
  plateau_output <- plateau_overall
  lastestimase <- estimase
  jminitiation <-  count2
  crminitiation <- count1
  survivaldata <-  survdata
  survivaltemp <- survtemp
  longitudinaldata <-  longdata
  longitudinaltemp <- longtemp
  lastestimations <-  estima
  toxic <- toxic
  stop_moment <- stop_moment
  
 
  
  # SAVED ELEMENTS AT THE END OF THE TRIAL
  # doseassignment: number of subjects assigned at each dose level and final recommendation
  # lastMTD: the MTD at the end of the trial
  # plateau_output: indicates the dose that is the beginning of the plateau at the final recommendation
  # jminitiation: the number of subjects enrolled in the trial when proceeding to the JM
  # crminitiation: the number of subjects enrolled in the trial when proceeding to the crm
  # survivaldata: the survival data at the end of the trial
  # longitudinaldata: the longitudinal data at the end of the trial
  # survivaltemp: the survival data at the end of the trial in real time. For instance for the last subject only the information of the first cycle is listed
  # longitudinaltemp: the longitudinal data at the end of the trial in real time
  # lastestimations: the model estimations used to select the dose for the last subject
  # lastestimase: SEs of the above estimations
  # toxic: indicates whether the trial stopped for excessive toxicity or not
  # stop_moment: indicates at what point the trial stopped for excessive toxicity (2+2 or CRM or not at all) 
  out <- return(list(doseassignment=doseassignment, lastMTD=lastMTD, plateau_output=plateau_output, jminitiation=jminitiation, crminitiation=crminitiation, survivaldata=survivaldata, longitudinaldata=longitudinaldata, 
                     survivaltemp=survivaltemp, longitudinaltemp=longitudinaltemp, lastestimations=lastestimations, lastestimase=lastestimase, toxic=toxic, stop_moment=stop_moment))
}











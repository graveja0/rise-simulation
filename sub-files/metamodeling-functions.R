#Function to determine the number of ticks in ggplot2 plots
number_ticks <- function(n) {function(limits) pretty(limits, n)} #function to determine limits of axis' ticks 
number_ticks_x0 <- function(n,change) {function(limits) c(pretty(limits, n),change)} #function to determine limits of axis' ticks 
fmt <- function(){function(x) format(x,nsmall = 2,scientific = FALSE)}


Get_ICER <- function(result) {
  result %>% arrange(desc(dQALY)) %>% 
    mutate(ICER = (dCOST-lead(dCOST))/(dQALY-lead(dQALY))) %>% 
    mutate(dominated.strong = as.integer(lag(ICER)<0)) -> temp
  
  dominated.strong <- temp %>% select(strategy,dominated.strong)
  dominated.strong[is.na(dominated.strong)] = 0
  temp %>% filter(dominated.strong!=1|is.na(dominated.strong)) %>% 
    mutate(ICER = (dCOST-lead(dCOST))/(dQALY-lead(dQALY))) %>% 
    mutate(dominated.extended = as.integer(ICER>lag(ICER))) -> temp2 
  
  dominated.extended = temp2 %>% select(strategy,dominated.extended)
  dominated.extended[is.na(dominated.extended)] = 0
  
  temp2 %>% filter(dominated.extended!=1|is.na(dominated.extended)) %>% 
    mutate(ICER = (dCOST-lead(dCOST))/(dQALY-lead(dQALY))) %>% 
    select(strategy,cost=dCOST,effectiveness=dQALY,ICER=ICER) %>% select(strategy,ICER) -> out1
  
  result %>% full_join(out1,"strategy") %>% full_join(dominated.strong,"strategy") %>% full_join(dominated.extended,"strategy") %>% arrange(desc(dQALY))
}



#Define offset as a new axis transformation. 
offset_trans <- function(offset=0) {
  trans_new(paste0("offset-", format(offset)), function(x) x-offset, function(x) x+offset)
}

OneWaySA<-function(Strategies,Parms,Outcomes,parm,range){
  #Extract parameter column number in Parms matrix
  x<-which(colnames(Parms)==parm)
  dep<-length(Strategies) #Number of dependent variables, i.e., strategies outcomes
  indep<-ncol(Parms) #Number of independent variables, i.e., parameters
  Sim <- data.frame(Outcomes,Parms)
  #Determine range of of the parameer to be plotted
  if (!missing("range")){ #If user defines a range
    vector<-seq(range[1],range[2],length.out=400)
  }
  else{ #Default range given by the domanin of the parameter's sample
    #vector to define 400 samples between the 2.5th and 97.5th percentiles
    y = seq(2.5,97.5,length=400) 
    j = round(y*(length(Parms[,x])/100)) #indexing vector;j=round(y*n/100) where n is the size of vector of interest
    vector<-sort(Parms[j,x])    
  }
  
  #Generate a formula by pasting column names for both dependent and independent variables. Imposes a 1 level interaction
  f <- as.formula(paste('cbind(',paste(colnames(Sim)[1:dep],collapse=','), ') ~ (','poly(',parm,',2)+' ,paste(colnames(Parms)[-x], collapse='+'),')'))
  #Run Multiple Multivariate Regression (MMR) Metamodel
  Oway.mlm = lm(f,data=Sim)
  
  #Generate matrix to use for prediction 
  Sim.fit<-matrix(rep(colMeans(Parms)),nrow=length(vector),ncol=ncol(Parms), byrow=T)
  Sim.fit[,x]<-vector
  Sim.fit<-data.frame(Sim.fit) #Transform to data frame, the format required for predict
  colnames(Sim.fit)<-colnames(Parms) #Name data frame's columns with parameters' names
  
  #Predict Outcomes using MMMR Metamodel fit
  plotdata = data.frame(predict(Oway.mlm, newdata = Sim.fit))
  colnames(plotdata) <- Strategies #Name the predicted outcomes columns with strategies names

  #Reshape dataframe for ggplot
  plotdata = stack(plotdata, select=Strategies) #
  plotdata = cbind(Sim.fit, plotdata) #Append parameter's dataframe to predicted outcomes dataframe
  
  #A simple trick to define my variables in my functions environment
  plotdata$parm<-plotdata[,parm];
  
  txtsize<-12 #Text size for the graphs
  ggplot(data = plotdata, aes(x = parm, y = values, color = ind)) +
    geom_line() +
    ggtitle("One-way sensitivity analysis \n Net Health Benefit") + 
    xlab(parm) +
    ylab("E[NHB]") +
    scale_colour_hue("Strategy", l=50) +
    scale_x_continuous(breaks=number_ticks(6)) + #Adjust for number of ticks in x axis
    scale_y_continuous(breaks=number_ticks(6)) +
    theme_bw() +
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
}

TwoWaySA<-function(Strategies,Parms,Outcomes,parm1,parm2,range1,range2){
  #Extract parameter column number in Parms matrix
  x1<-which(colnames(Parms)==parm1)
  x2<-which(colnames(Parms)==parm2)
  dep<-length(Strategies) #Number of dependent variables, i.e., strategies
  indep<-ncol(Parms) #Number of independent variables, i.e., parameters
  
  Sim <- data.frame(Outcomes,Parms)
  
  if (ncol(Parms)==2) {
    Parms$constant = 1
    Sim$constant = 1
  }
  
  #Determine range of of the parameer to be plotted
  if (!missing("range1")&!missing("range2")){ #If user defines a range
    vector1<-seq(from=range1[1],to=range1[2],length.out=301)
    vector2<-seq(from=range2[1],to=range2[2],length.out=301)
  } else if (!missing("range1")&missing("range2")){ #Default range given by the domanin of the parameter's sample
    #vector to define 400 samples between the 2.5th and 97.5th percentiles
    vector1<-seq(from=range1[1],to=range1[2],length.out=301)
    y2 = seq(2.5,97.5,length.out=301)
    j2 = round(y2*(length(Parms[,x2])/100)) #indexing vector;j=round(y*n/100) where n is the size of vector of interest
    vector2<-sort(Parms[j2,x2])    
  } else if (missing("range1")&!missing("range2")){ #Default range given by the domanin of the parameter's sample
    #vector to define 400 samples between the 2.5th and 97.5th percentiles
    vector2<-seq(from=range2[1],to=range2[2],length.out=301)
    y1 = seq(2.5,97.5,length.out=301)
    j1 = round(y1*(length(Parms[,x1])/100)) #indexing vector;j=round(y*n/100) where n is the size of vector of interest
    vector1<-sort(Parms[j1,x1])    
  } else{
    y1 = seq(2.5,97.5,length.out=301) 
    y2 = seq(2.5,97.5,length.out=301) 
    j1 = round(y1*(length(Parms[,x1])/100)) #indexing vector;j=round(y*n/100) where n is the size of vector of interest
    j2 = round(y2*(length(Parms[,x2])/100))
    vector1<-sort(Parms[j1,x1])
    vector2<-sort(Parms[j2,x2])
  }
  #Generate a formula by pasting column names for both dependent and independent variables
  f <- as.formula(paste('cbind(',paste(colnames(Sim)[1:dep],collapse=','), 
                        ') ~ (','poly(',parm1,',2)+','poly(',parm2,',8)+' ,
                        paste(colnames(Parms)[c(-x1,-x2)], collapse='+'),')'))
  #Run Multiple Multivariate Regression (MMR) Metamodel
  Tway.mlm = lm(f,data=Sim)
  
  TWSA <- expand.grid(parm1=vector1,parm2=vector2)
  
  #Generate matrix to use for prediction 
  Sim.fit<-matrix(rep(colMeans(Parms)),nrow=nrow(TWSA),ncol=ncol(Parms), byrow=T)
  Sim.fit[,x1]<-TWSA[,1]
  Sim.fit[,x2]<-TWSA[,2]
  Sim.fit<-data.frame(Sim.fit) #Transform to data frame, the format required for predict
  colnames(Sim.fit)<-colnames(Parms) #Name data frame's columns with parameters' names
  
  #Predict Outcomes using MMMR Metamodel fit
  Sim.TW = data.frame(predict(Tway.mlm, newdata = Sim.fit))
  #Find optimal strategy in terms of maximum Outcome
  Optimal <- max.col(Sim.TW)
  #Get Outcome of Optimal strategy
  OptimalOut<-apply(Sim.TW,1,max)
  
  plotdata = Sim.fit #Append parameter's dataframe to predicted outcomes dataframe
  
  #A simple trick to define my variables in my functions environment
  plotdata$parm1<-plotdata[,parm1];
  plotdata$parm2<-plotdata[,parm2];
  
  plotdata$Strategy<-factor(Optimal,labels=Strategies[as.numeric(names(table(Optimal)))])
  plotdata$value<-OptimalOut
  
  txtsize<-12
  ggplot(plotdata, aes(x=parm1,y=parm2))+ 
    geom_tile(aes(fill=Strategy)) +
    theme_bw() +
    ggtitle(expression(atop("Two-way sensitivity analysis", 
                            atop("Net Health Benefit")))) + 
    scale_fill_discrete("Strategy: ", l=50)+
    xlab(parm1)+
    ylab(parm2)+
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
}

TornadoOpt <-function(Parms,Outcomes){
  # Grouped Bar Plot
  # Determine the overall optimal strategy
  opt<-which.max(colMeans(Outcomes))
  # calculate min and max vectors of the parameters (e.g., lower 2.5% and 97.5%)
  X <- as.matrix(Parms)
  y <- as.matrix(Outcomes[,opt])
  ymean <- mean(y)
  n <- nrow(Parms)
  nParams <- ncol(Parms)
  #paramNames <- Names[seq(8,7+nParams)]
  paramNames <- colnames(Parms)
  Parms.sorted <- apply(Parms,2,sort,decreasing=F)#Sort in increasing order each column of Parms
  lb <- 2.5
  ub <- 97.5 
  Xmean <- rep(1,nParams) %*% t(colMeans(X))
  XMin <- Xmean
  XMax <- Xmean
  paramMin <- as.vector(Parms.sorted[round(lb*n/100),])
  paramMax <- as.vector(Parms.sorted[round(ub*n/100),])
  paramNames2 <- paste(paramNames, "[", round(paramMin,2), ",", round(paramMax,2), "]")
  
  diag(XMin) <- paramMin
  diag(XMax) <- paramMax
  
  XMin <- cbind(1, XMin)
  XMax <- cbind(1, XMax)
  
  X <- cbind(1,X)
  B <- solve(t(X) %*% X) %*% t(X) %*% y
  yMin <- XMin %*% B - ymean
  yMax <- XMax %*% B - ymean
  ySize <- abs(yMax - yMin) 
  
  rankY<- order(ySize)
  xmin <- min(c(yMin, yMax)) + ymean
  xmax <- max(c(yMin, yMax)) + ymean
  
  Tor <- data.frame(
    Parameter=c(paramNames2[rankY],paramNames2[rankY]),  
    Level=c(rep("Low",nParams),rep("High",nParams)),
    value=ymean+c(yMin[rankY],yMax[rankY]),
    sort=seq(1,nParams)
  )
  #re-order the levels in the order of appearance in the data.frame
  Tor$Parameter2 <- factor(Tor$Parameter, as.character(Tor$Parameter))
  #Define offset as a new axis transformation. Source: http://blog.ggplot2.org/post/25938265813/defining-a-new-transformation-for-ggplot2-scales  
  offset_trans <- function(offset=0) {
    trans_new(paste0("offset-", format(offset)), function(x) x-offset, function(x) x+offset)
  }
  #Plot the Tornado diagram.
  txtsize<-12
  ggplot(Tor[Tor$Level=="Low",], aes(x=Parameter2,y=value, fill=Level)) +
    geom_bar(stat="identity") +
    ggtitle("Tornado Diagram")+
    scale_fill_discrete("Parameter Level: ", l=50)+
    scale_y_continuous(name="Net Benefit",trans=offset_trans(offset=ymean)) +
    scale_x_discrete(name="Parameter") +
    geom_bar(data=Tor[Tor$Level=="High",], aes(x=Parameter2,y=value, fill=Level), stat="identity") +
    geom_hline(yintercept = ymean, linetype = "dotted", size=0.5) +
    theme_bw()+
    theme(legend.position="bottom",legend.title=element_text(size = txtsize,angle = 0, hjust = 1),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize),
          axis.ticks.y = element_blank())+
    coord_flip()  
}

TornadoAll <-function(Strategies,Parms,Outcomes){
  opt<-which.max(colMeans(Outcomes))
  # calculate min and max vectors of the parameters (e.g., lower 2.5% and 97.5%)
  X <- as.matrix(Parms)
  y <- as.matrix(Outcomes[,opt])
  Y <- as.matrix(Outcomes)
  ymean <- mean(y)
  n <- nrow(Parms)
  nParams <- ncol(Parms)
  #paramNames <- Names[seq(8,7+nParams)]
  paramNames <- colnames(Parms)
  Parms.sorted <- apply(Parms,2,sort,decreasing=F)#Sort in increasing order each column of Parms
  lb <- 2.5
  ub <- 97.5 
  Xmean <- rep(1,nParams) %*% t(colMeans(X))
  XMin <- Xmean
  XMax <- Xmean
  paramMin <- as.vector(Parms.sorted[round(lb*n/100),])
  paramMax <- as.vector(Parms.sorted[round(ub*n/100),])
  
  diag(XMin) <- paramMin
  diag(XMax) <- paramMax
  
  XMin <- cbind(1, XMin)
  XMax <- cbind(1, XMax)
  
  X <- cbind(1,X)
  B <- solve(t(X) %*% X) %*% t(X) %*% y
  #install.packages("matrixStats")
  library(matrixStats)
  bigBeta <- solve(t(X) %*% X) %*% t(X) %*% Y
  yMin <- rowMaxs(XMin %*% bigBeta - ymean)
  yMax <- rowMaxs(XMax %*% bigBeta - ymean)
  ySize <- abs(yMax - yMin) 
  
  rankY<- order(ySize)
  xmin <- min(c(yMin, yMax)) + ymean
  xmax <- max(c(yMin, yMax)) + ymean
  
  paramNames2 <- paste(paramNames, "[", round(paramMin,2), ",", round(paramMax,2), "]")
  
  strategyNames<-Strategies
  strategyColors <- c("red","darkgreen","blue")
  
  ## Polygon graphs:
  nRect <- 0
  x1Rect <- NULL
  x2Rect <- NULL
  ylevel <- NULL
  colRect <- NULL
  
  for (p in 1:nParams){
    xMean <- colMeans(X)
    xStart = paramMin[rankY[p]]
    xEnd = paramMax[rankY[p]]
    xStep = (xEnd-xStart)/1000
    for (x in seq(xStart,xEnd, by = xStep)){
      #for each point determine which one is the optimal strategy
      xMean[rankY[p] + 1] <- x 
      yOutcomes <- xMean %*% bigBeta
      yOptOutcomes <- max(yOutcomes)
      yOpt <- which.max(yOutcomes)
      if (x == xStart){
        yOptOld <- yOpt
        y1 <- yOptOutcomes
      }
      #if yOpt changes, then plot a rectangle for that region
      if (yOpt != yOptOld | x == xEnd){
        nRect <- nRect + 1
        x1Rect[nRect] <- y1
        x2Rect[nRect] <- yOptOutcomes
        ylevel[nRect] <- p
        colRect[nRect] <- strategyColors[yOptOld]
        yOptOld <- yOpt
        y1 <- yOptOutcomes
      }
    }
  }
  
  txtsize <- 12
  d=data.frame(x1=x2Rect, x2=x1Rect, y1=ylevel-0.4, y2=ylevel+0.4, t=colRect, r = ylevel)
  ggplot(d, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = t)) +
    ggtitle("Tornado Diagram") + 
    xlab("Expected NHB") +
    ylab("Parameters") + 
    geom_rect()+
    theme_bw() + 
    scale_y_continuous(limits = c(0.5, nParams + 0.5),breaks=seq(1:ncol(Parms)), labels=paramNames2[rankY]) +
    #scale_y_discrete(breaks=seq(1:8), labels=paramNames2[rankY]) + 
    scale_fill_discrete(name="Optimal\nStrategy",
                        #breaks=c("ctrl", "trt1", "trt2"),
                        labels=strategyNames,
                        l=50) + 
    geom_vline(xintercept=ymean, linetype="dotted") + 
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
}


###copy TornadoAll and add threshold values as second y-axis
TornadoAll2 <-function(Strategies,Parms,Outcomes,Highlight){ #add parameter name list to highlight labels
  opt<-which.max(colMeans(Outcomes))
  # calculate min and max vectors of the parameters (e.g., lower 2.5% and 97.5%)
  X <- as.matrix(Parms)
  y <- as.matrix(Outcomes[,opt])
  Y <- as.matrix(Outcomes)
  ymean <- mean(y)
  n <- nrow(Parms)
  nParams <- ncol(Parms)
  #paramNames <- Names[seq(8,7+nParams)]
  paramNames <- colnames(Parms)
  Parms.sorted <- apply(Parms,2,sort,decreasing=F)#Sort in increasing order each column of Parms
  lb <- 2.5
  ub <- 97.5 
  Xmean <- rep(1,nParams) %*% t(colMeans(X))
  XMin <- Xmean
  XMax <- Xmean
  paramMin <- as.vector(Parms.sorted[round(lb*n/100),])
  paramMax <- as.vector(Parms.sorted[round(ub*n/100),])
  
  diag(XMin) <- paramMin
  diag(XMax) <- paramMax
  
  XMin <- cbind(1, XMin)
  XMax <- cbind(1, XMax)
  
  X <- cbind(1,X)
  B <- solve(t(X) %*% X) %*% t(X) %*% y
  #install.packages("matrixStats")
  library(matrixStats)
  bigBeta <- solve(t(X) %*% X) %*% t(X) %*% Y
  yMin <- rowMaxs(XMin %*% bigBeta - ymean)
  yMax <- rowMaxs(XMax %*% bigBeta - ymean)
  ySize <- abs(yMax - yMin) 
  
  rankY<- order(ySize)
  xmin <- min(c(yMin, yMax)) + ymean
  xmax <- max(c(yMin, yMax)) + ymean
  
  paramNames2 <- paste(paramNames, "[", round(paramMin,2), ",", round(paramMax,2), "]")
  
  strategyNames<-Strategies
  strategyColors <- c("red","darkgreen","blue")
  
  ## Polygon graphs:
  nRect <- 0
  x1Rect <- NULL
  x2Rect <- NULL
  ylevel <- NULL
  colRect <- NULL
  namRect <- NULL #add for validation
  paraRect <- NULL #add for validation
  
  for (p in 1:nParams){
    xMean <- colMeans(X)
    xStart = paramMin[rankY[p]]
    xEnd = paramMax[rankY[p]]
    xStep = (xEnd-xStart)/1000
    #dd <- data.frame(NULL) #used for validation
    for (x in seq(xStart,xEnd, by = xStep)){
      #for each point determine which one is the optimal strategy
      xMean[rankY[p] + 1] <- x 
      yOutcomes <- xMean %*% bigBeta
      yOptOutcomes <- max(yOutcomes)
      yOpt <- which.max(yOutcomes)
      if (x == xStart){
        yOptOld <- yOpt
        y1 <- yOptOutcomes
      }
      #if yOpt changes, then plot a rectangle for that region
      if (yOpt != yOptOld | x == seq(xStart,xEnd, by = xStep)[1001]){ #x==xEnd fails to catch 1 end point
        nRect <- nRect + 1
        x1Rect[nRect] <- y1
        x2Rect[nRect] <- yOptOutcomes
        ylevel[nRect] <- p
        colRect[nRect] <- strategyColors[yOptOld]
        namRect[nRect] <- strategyNames[yOptOld] #add for validation
        paraRect[nRect] <- paramNames[rankY[p]] #add for validation
        yOptOld <- yOpt
        y1 <- yOptOutcomes
      }
      
      ## used for validation
      # temp <- data.frame(x=x,yOpt=yOpt,ct=nRect) %>% cbind(data.frame(yOutcomes))
      # if(length(dd)==0) {dd <- temp } else {dd <- rbind(dd,temp)}
      
    }
  }
  
  ###calculate threshold
  testx <- colMeans(X)
  testb <- bigBeta
  threshold <- c()
  for (i in 1:nParams){
    #grab metamodel estimations and calculate all intersections
    m1 <- testb[(i+1),]
    m2 <- colSums((testx*testb)[-(i+1),])
    xx <- (-1)*sapply(combn(m2,2,simplify = F),FUN=diff)/sapply(combn(m1,2,simplify = F),FUN=diff)
    if(any(xx<=paramMax[i] & xx>=paramMin[i])) { #2.5%~97.5%
      temp <- xx[xx<=paramMax[i] & xx>=paramMin[i]] #must be within parameter range
      for(j in seq(temp)) {
        tt <- temp[j]*m1+m2 
        if(length(which(tt==max(tt)))==1) {temp <- temp[-j]} #drop this intersection if a single strategy dominates all
      }
      threshold[i] <- paste(round(temp,digits=2)) 
    } else {threshold[i] <- "NA"}
  }
  
  txtsize <- 12
  d=data.frame(x1=x2Rect, x2=x1Rect, y1=ylevel-0.4, y2=ylevel+0.4, t=colRect, n=namRect, para=paraRect, r = ylevel)
  #customize how to highlight parameters (color & font)
  couleur <- ifelse(paramNames[rankY] %in% Highlight,"red","black")
  ecriture <- ifelse(paramNames[rankY] %in% Highlight,"bold","plain")
  
  ggplot(d, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill=n)) +
    ggtitle("Torando Diagram") + 
    xlab("Expected NHB") +
    ylab("Parameters") + 
    geom_rect()+
    scale_fill_manual(
      name="Optimal\nStrategy",
      labels =strategyNames,
      values =strategyColors) +
    theme_bw() + 
    scale_y_continuous(limits = c(0.5, nParams + 0.5),breaks=seq(1:ncol(Parms)), labels=paramNames2[rankY],
                       ###add threshold label
                       sec.axis = sec_axis(~.,breaks=derive(),name="Threshold",labels = threshold[rankY])) +
    #scale_y_discrete(breaks=seq(1:8), labels=paramNames2[rankY]) + 
    # scale_fill_discrete(name="Optimal\nStrategy",
    #                     #breaks=c("ctrl", "trt1", "trt2"),
    #                     labels=rev(strategyNames),
    #                     l=50) + 
    geom_vline(xintercept=ymean, linetype="dotted") + 
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize,color=couleur,face=ecriture), 
          #if not want to highlight right y axis, add axis.text.y.right specification
          axis.text.x = element_text(size=txtsize)) 
  
}



PlaneCE<-function(Strategies,Outcomes){
  ndep<-length(Strategies)*2 #Determine number of outcomes for all starteges, i.e., cost and effectiveness
  ind_c<- grep("COST",colnames(Outcomes)) #seq(1,(ndep-1),by=2) #Index to extract the costs from matrix Outcomes
  ind_e<-grep("QALY",colnames(Outcomes)) #seq(2,ndep,by=2) #Index to extract the effectiveness from matrix Outcomes
  Cost<-reshape2::melt(Outcomes[,ind_c],variable_name="Strategy") %>% rename(Strategy=variable)
  levels(Cost$Strategy)<- gsub("dCOST_","",levels(Cost$Strategy))
  Eff<-reshape2::melt(Outcomes[,ind_e])   %>% rename(Strategy=variable)
  levels(Eff$Strategy)<- gsub("dQALY_","",levels(Eff$Strategy))
  CE<-cbind(Cost,Eff[,2])
  colnames(CE)<-c("Strategy","Cost","Effectiveness")
  
  #Dataframe with means of strategies.
  Means <- CE %>% group_by(Strategy) %>% summarise(
                 N = length(Cost),
                 Cost.mean = mean(Cost),
                 Eff.mean = mean(Effectiveness))
  
  #Define ggplot object
  txtsize<-12
  ggplot(Means, aes(x = Eff.mean, y = Cost.mean, color=Strategy)) + 
    geom_point(size=4, aes(shape=Strategy)) +
    ggtitle("Cost-Effectiveness Plane") +
    scale_colour_discrete(l=50) +  # Use a slightly darker palette than normal
    scale_y_continuous(labels = dollar)+
    scale_x_continuous(breaks=number_ticks(6))+
    xlab("Effectiveness")+
    ylab("Cost")+
    theme_bw() +
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
}

ScatterCE<-function(Strategies,Outcomes){
  ndep<-length(Strategies)*2 #Determine number of outcomes for all starteges, i.e., cost and effectiveness
  ind_c<- grep("COST",colnames(Outcomes)) #seq(1,(ndep-1),by=2) #Index to extract the costs from matrix Outcomes
  ind_e<-grep("QALY",colnames(Outcomes)) #seq(2,ndep,by=2) #Index to extract the effectiveness from matrix Outcomes
  Cost<-reshape2::melt(Outcomes[,ind_c],variable_name="Strategy") %>% rename(Strategy=variable)
  levels(Cost$Strategy)<- gsub("dCOST_","",levels(Cost$Strategy))
  Eff<-reshape2::melt(Outcomes[,ind_e])   %>% rename(Strategy=variable)
  levels(Eff$Strategy)<- gsub("dQALY_","",levels(Eff$Strategy))
  CE<-cbind(Cost,Eff[,2])
  colnames(CE)<-c("Strategy","Cost","Effectiveness")
  
  # Ellipses code
  df_ell <- data.frame() #create an empty dataframe
  # for each level in df$groups 
  for(g in levels(CE$Strategy)){
    # create 100 points per variable around the mean of each group
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(CE[CE$Strategy==g,], 
                                                     ellipse(cor(Effectiveness, Cost), 
                                                             scale=c(sd(Effectiveness),sd(Cost)), 
                                                             centre=c(mean(Effectiveness),mean(Cost)))
    )),group=g))
  }
  Means <- CE %>% group_by(Strategy) %>% summarise(
                 N = length(Cost),
                 Cost.mean = mean(Cost),
                 Eff.mean = mean(Effectiveness))  
  #Define ggplot object
  txtsize<-12
  ggplot(CE, aes(x=Effectiveness, y=Cost, color=Strategy)) + 
    geom_point(size=0.7) +
    geom_point(data = Means,aes(x = Eff.mean, y = Cost.mean, shape=Strategy),size=8,fill="white") +
    geom_text(data = Means,aes(x = Eff.mean, y = Cost.mean, label=c(seq(nrow(Means)))),size=5,colour="gray",alpha=1) +
    geom_path(data=df_ell, aes(x=x, y=y,colour=group), size=1, linetype=2, alpha=1) + # draw ellipse lines
    ggtitle("Cost-Effectiveness Scatterplot") +
    scale_colour_discrete(l=50) +  # Use a slightly darker palette than normal
    scale_y_continuous(labels = dollar)+
    scale_x_continuous(breaks=number_ticks(6))+
    theme_bw() +
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
}

CEAC<-function(range,Strategies,Outcomes){
  # Outcomes must be ordered in a way that for each strategy the cost must appear first then the effectiveness
  lambda<-seq(range[1],range[2], length.out=range[3])
  NHB <- array(0, dim=c(dim(Outcomes)[1],length(Strategies))) # Matrix to store NHB for each strategy
  colnames(NHB)<-Strategies
  CEA<-array(0,dim=c(length(lambda),length(Strategies)))
  
  for(l in 1:length(lambda)){
    for(i in 1:dim(NHB)[2]){
      NHB[,i] <-  Outcomes[,(2*i)]-Outcomes[,(2*i-1)]/lambda[l] # Effectiveness minus Costs
    }
    Max.NHB <- max.col(NHB)
    Optimal <- array(0, dim=c(dim(Outcomes)[1],length(Strategies))) # Matrix with dummy variables indicating if strategy is optimal
    for (j in 1:dim(Optimal)[1]){
      k<-Max.NHB[j]
      Optimal[j,k]<-1
    }
    CEA[l,]<-colMeans(Optimal)
  }
  CEA<-data.frame(cbind(lambda,CEA))
  colnames(CEA)<-c("Lambda", Strategies)
  
  CEAC<-reshape2::melt(CEA, id.vars = "Lambda") 
  
  txtsize<-12
  ggplot(data = CEAC, aes(x = Lambda, y = value, color = variable)) +
    geom_point() +
    geom_line() +
    ggtitle("Cost-Effectiveness Acceptability Curves") + 
    scale_colour_hue("Strategies: ",l=50) +
    scale_x_continuous(breaks=number_ticks(6))+
    xlab(expression("Willingness to Pay "(lambda))) +
    ylab("Pr Cost-Effective") +
    theme_bw() +
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
}



survSplit<-function(data, cut, end,event,start,id=NULL,
                    zero=0,episode=NULL){

  cut<-sort(cut)
  ntimes <- length(cut)
  n <- nrow(data)
  p <- ncol(data)
  
  newdata <- lapply(data,rep,ntimes+1)

  endtime <- rep(c(cut, Inf) ,each=n)

  eventtime<-newdata[[end]]

  if( start %in% names(data))
    starttime<-data[[start]]
  else
    starttime<-rep(zero,length.out=n)

  starttime<-c(starttime, pmax(starttime, rep(cut,each=n)))
  
  epi<-rep(0:ntimes,each=n)
  
  status <- ifelse( eventtime <= endtime & eventtime>starttime,
                   newdata[[event]], 0)
  endtime<- pmin(endtime,eventtime)

  drop<-starttime>=endtime
  
  newdata<-do.call("data.frame",newdata)
  newdata[,start]<-starttime
  newdata[,end]<-endtime
  newdata[,event]<-status
  if (!is.null(id))
    newdata[,id]<-rep(rownames(data),ntimes+1)
  if (!is.null(episode))
    newdata[,episode]<-epi
  
  newdata<-newdata[!drop,]

  newdata

}

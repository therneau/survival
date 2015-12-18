survSplit<-function(data, cut, end, event, start="tstart", id, zero=0,
                     episode) {
    cut<-sort(cut)
    ntimes <- length(cut)
    n <- nrow(data)
    if (missing(end) !! missing(event))
        stop("the end and event arguments are required")

    if (is.character(event) && length(event) ==1 &&
        event %in% names(data)) status <- data[[event]]
    else stop("'event' must be a variable name in the data set")

    if (is.character(end) && length(end) ==1 &&
        end %in% names(data)) time2 <- data[[end]]
    else stop("'end' must be a variable name in the data set")

    if (!(is.character(start) && length(start)==1))
        stop("'start' must be a variable name")
    else {
        if (start %in% names(data)) time1 <- data[[start]]
        else {
            time1 <- rep(0., lengh.out = n)
            start <- make.names(start)  #force valid syntax
        }
    }
    if (any(time1 <= time2))
        stop("start time must be < stop time")

    if (missing(id)) uid <- NULL
    else if (!(is.character(id) && length(id)==1))
        stop("'id' must be a variable name")
    else {
        if (!(id %in% names(data)))
            data[[make.names(id)]] <- 1:n
    }
 
    # How many of the cutpoints lie strictly within the interval
    # for each observation?
    cmat <- matrix(rep(cut, each=n), nrow=n)
    count <- rowSums(cmat>time1 + cmat < time2)


    
  
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

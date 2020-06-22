print.survcheck <- function(x, ...){
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
    }
    temp <- x$n
    names(temp) <- c("Unique identifiers", "Observations", "Transitions")
    print(temp)
    if(!is.null(x$na.action)){
        cat(length(x$na.action), "observations removed due to missing","\n")
    } 
    
    cat("\nTransitions table:\n")
    print(x$transitions)
    cat('\n')

    ## how many of each state does each id have?
    cat("Number of subjects with 0, 1, ... transitions to each state:\n")
    print(x$events)
    cat("\n")
    
    if(x$flag["overlap"]>0) {
        cat("Overlap check: ", 
            length(x$overlap$id),
            ifelse(length(x$overlap$id)==1," id (", " ids ("),
            length(x$overlap$row),
            " rows)\n", sep="")
    } 
    if(x$flag["gap"]>0) {
        cat("Gap check: ", 
            length(x$gap$id),
            ifelse(length(x$gap$id)==1," id (", " ids ("),
            length(x$gap$row),
            " rows)\n", sep="")
    } 
    if(x$flag["teleport"]>0) {
        cat("Teleport check: ", 
            length(x$teleport$id),
            ifelse(length(x$teleport$id)==1," id (", " ids ("),
            length(x$teleport$row),
            " rows)\n", sep="")
    } 
    if(x$flag["jump"] >0){
        cat("Jump check: ", 
            length(x$jump$id),
            ifelse(length(x$jump$id)==1," id (", " ids ("),
            length(x$jump$row),
            " rows)\n", sep="")
    } 
    if(x$flag["duplicate"] >0) {
        cat("Duplicate time check: ", 
            length(x$duplicate$id),
            ifelse(length(x$duplicate$id)==1," id (", " ids ("),
            length(x$duplicate$row),
            " rows)\n", sep="")
    } 
}


summary.survcheck <- function(object, max.show=5, ...){
    if (!is.null(cl <- object$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
    }

    temp <- object$n
    names(temp) <- c("Unique identifiers", "Observations", "Transitions")
    print(temp)

    if(!is.null(object$na.action)){
        cat(length(object$na.action), "observations removed due to missing","\n")
    }
    
    ## change between states
    cat("Transitions table:\n")
    print(object$transitions)
    cat('\n')

    ## how many of each state does each id have?
    cat("Number of subjects with 1, 2, ... copies of each state:\n")
    print(object$events)
    cat("\n")
    
    if(object$flag["overlap"]>0) {
        cat("Overlap: ", 
            length(object$overlap$id),
            ifelse(length(object$overlap$id)==1," id (", " ids ("),
            length(object$overlap$row),
            " rows)\n", sep="")
        odat <- data.frame(id=object$id, y= object$Y)[object$overlap$row,]
        nshow <- min(max.show, length(object$overlap$id))
        print(odat[odat[,1] %in% object$overlap$id[1:nshow],])
        cat('\n')
    } 
    if(object$flag["gap"]>0) {
        cat("Gap: ", 
            length(object$gap$id),
            ifelse(length(object$gap$id)==1," id (", " ids ("),
            length(object$gap$row),
            " rows)\n", sep="")
        gdat <- data.frame(id=object$id, y=object$Y)[object$gap$row,]
        nshow <- min(max.show, length(object$gap$id))
        print(gdat[gdat[,1] %in% object$gap$id[1:nshow],])
        cat('\n')
    } 

    if(object$flag["teleport"]>0) {
        cat("Teleport: ", 
            length(object$teleport$id),
            ifelse(length(object$teleport$id)==1," id (", " ids ("),
            length(object$teleport$row),
            " rows)\n", sep="")
        tdat <- data.frame(id=object$id,object$Y)[object$teleport$row,]
        nshow <- min(max.show, length(object$teleport$id))
        print(tdat[tdat[,1] %in% object$teleport$id[1:nshow],])
        cat('\n')
    } 
 
   if(object$flag["jump"] >0){
        cat("Jump: ", 
            length(object$jump$id),
            ifelse(length(object$jump$id)==1," id (", " ids ("),
            length(object$jump$row),
            " rows)\n", sep="")
        jdat <- data.frame(id=object$id, y= object$Y)[object$jump$row,]
        nshow <- min(max.show, length(object$jump$id))
        print(jdat[jdat[,1] %in% object$jump$id[1:nshow],])
        cat('\n')
    }  

   if(object$flag["duplicate"] >0){
        cat("Duplicate times: ", 
            length(object$duplicate$id),
            ifelse(length(object$duplicate$id)==1," id (", " ids ("),
            length(object$duplicate$row),
            " rows)\n", sep="")
        jdat <- data.frame(id=object$id, y= object$Y)[object$duplicate$row,]
        nshow <- min(max.show, length(object$duplicate$id))
        print(jdat[jdat[,1] %in% object$duplicate$id[1:nshow],])
        cat('\n')
    }  
}

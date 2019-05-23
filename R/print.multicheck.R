print.multicheck <- function(x, ...){
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
    }
    #browser()
    if(!is.null(x$na.action)){
        cat(length(x$na.action), "observations with missing values:","\n")
    }
    cat(sum(x$statecount[,1]),"subjects available for analysis","\n")
    
    ## how many of each state does each id have?
    cat("Number of subjects with 1, 2, ... copies of each state:\n")
    print(x$statecount[-1,-ncol(x$statecount),drop=FALSE])
    cat("\n")
    
    ## change between states
    if(nrow(x$transitions)>1){
        cat("# subjects moving between states:\n")
        print(x$transitions)
        cat('\n')
    }
    
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
    if(x$flag["jump"]){
        cat("Jump check: ", 
            length(x$jump$id),
            ifelse(length(x$jump$id)==1," id (", " ids ("),
            length(x$jump$row),
            " rows)\n", sep="")
    } 
}


summary.multicheck <- function(object, max.show=5, ...){
    if (!is.null(cl <- object$call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
    }

    if(!is.null(object$na.action)){
        cat(length(object$na.action), "observations with missing values:","\n")
    }
    cat(sum(object$statecount[,1]),"subjects available for analysis","\n")
    
    ## how many of each state does each id have?
    cat("Number of subjects with 1, 2, ... copies of each state:\n")
    print(object$statecount[-1,-ncol(object$statecount),drop=FALSE])
    cat("\n")
    
    ## change between states
    if(nrow(object$transitions)>1){
        cat("# subjects moving between states:\n")
        print(object$transitions)
        cat('\n')
    }
    if(object$flag["overlap"]>0) {
        cat("Overlap check: ", 
            length(object$overlap$id),
            ifelse(length(object$overlap$id)==1," id (", " ids ("),
            length(object$overlap$row),
            " rows)\n", sep="")
        o.mat <- cbind(id=object$id,object$Y)[object$overlap$row,]
        max.show <- min(max.show, length(object$overlap$id))
        print(o.mat[o.mat[,1] %in% object$overlap$id[1:max.show],])
        cat('\n')
    } 
    if(object$flag["gap"]>0) {
        cat("Gap check: ", 
            length(object$gap$id),
            ifelse(length(object$gap$id)==1," id (", " ids ("),
            length(object$gap$row),
            " rows)\n", sep="")
        g.mat <- cbind(id=object$id,object$Y)[object$gap$row,]
        max.show <- min(max.show, length(object$gap$id))
        print(g.mat[g.mat[,1] %in% object$gap$id[1:max.show],])
        cat('\n')
    } 
    if(object$flag["teleport"]>0) {
        cat("Teleport check: ", 
            length(object$teleport$id),
            ifelse(length(object$teleport$id)==1," id (", " ids ("),
            length(object$teleport$row),
            " rows)\n", sep="")
        t.mat <- cbind(id=object$id,object$Y)[object$teleport$row,]
        max.show <- min(max.show, length(object$teleport$id))
        print(t.mat[t.mat[,1] %in% object$teleport$id[1:max.show],])
        cat('\n')
    } 
    if(object$flag["jump"]){
        cat("Jump check: ", 
            length(object$jump$id),
            ifelse(length(object$jump$id)==1," id (", " ids ("),
            length(object$jump$row),
            " rows)\n", sep="")
        j.mat <- cbind(id=object$id,object$Y)[object$jump$row,]
        max.show <- min(max.show, length(object$jump$id))
        print(j.mat[j.mat[,1] %in% object$jump$id[1:max.show],])
        cat('\n')
    }  
}

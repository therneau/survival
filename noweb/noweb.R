# See the inst/noweb directory in the coxme package for the original source.

nwread <- function(file) {
    program <- tab.to.blank(readLines(file))
    codestart <- which(grepl("^ *<<[^>]*>>=", program))
    textstart <- which(grepl("^@$|^@ ", program))
    program <- nwkillat(program)
    if (min(codestart, textstart) > 1) textstart <- c(1, textstart)
    
    temp <- rbind( c(codestart, textstart),
                   c(rep(2, length(codestart)), rep(1, length(textstart))))
    temp <- temp[,order(temp[1,])]
    endline <- c(temp[1,-1] -1, length(program))

    output <- vector("list", ncol(temp))  #number of chunks
    oname <- rep("", ncol(temp))
    for (i in 1:ncol(temp)) {
        if (temp[2,i]==1) { # text
            blankline <- sub("^@ *","", program[temp[1,i]])
            if (blankline=="" || substring(blankline,1,1)=="%") {
                # The line is blank
                if (temp[1,i]==endline[i])
                    text <- vector("character",0)  #Nothing there!
                else text <- program[(temp[1,i]+1):endline[i]]
                attr(text, "blankline") <- TRUE
                }
            else {
                text <- blankline
                if (temp[1,i] < endline[i])
                    text <- c(text, program[(temp[1,i]+1):endline[i]])
                attr(text, "blankline") <- FALSE
                }
            class(text) <- "nwtext"
            output[[i]] <- text
            }
        
        else {  #code
            cname <-  sub(">>=.*$", "", sub("<<", "", program[temp[1,i]]))
            if (temp[1,i] == endline[i]) code <- vector("character", 0)
            else code <- program[(temp[1,i]+1):endline[i]]
            oname[i] <- cname
            output[[i]] <- nwparse(code)
            }
        }
    
    names(output) <- oname

    for (i in which(oname!= "")) {  # all the code chunks
        if (any(is.na(match(output[[i]]$xref, oname)))) {
            indx <- which(is.na(match(output[[i]]$xref, oname)))
            stop(paste("Code referenced but not found:",
                       paste((output[[i]]$xref)[indx], collapse=", ")))
            }
        }
    
    temp <- nwloop(output)
    if (length(temp)) 
        stop(paste("Code structure has circular references: ",
                   paste(temp, collapse=" --> ")))

    class(output) <- "noweb"
    output
    }
nwparse <- function(lines) {
    # Look for references to other code
    indx <- which(grepl("<<[^>]*>>", lines))
    if (length(indx)) {
        xref <- sub(">>.*$", "", sub("[ \t]*<<", "", lines[indx]))
        indent <- sub("<<.*", "", lines[indx])
        out <- list(lines=lines, xref=xref, indent=indent, xindex=indx)
        }
    else out <- list(lines=lines, xref=NULL)
    
    class(out) <- "nwcode"
    out
    }
nwloop <- function(code) {   
    xref <- lapply(code, function(x) 
                   if (class(x)=="nwcode") unique(x$xref) else NULL)

    nwchase <- function(chain) {
        xtemp <- xref[[chain[1]]]  #routines called by the head of the chain
        if (length(xtemp) ==0) return(NULL)
        
        for (i in 1:length(xtemp)) {
            if (!is.na(match(xtemp[i], chain))) return(c(rev(chain), xtemp[i]))
            temp <- nwchase(c(xtemp[i], chain))
            if (!is.null(temp)) return(temp)
            }
        NULL
        }
    

    cnames <- names(code)
    temp <- lapply(cnames[cnames!=""], nwchase)
    templen <- sapply(temp,length)
    if (any(templen) > 0) 
        temp[[min(which(templen==min(templen[templen>0])))]]
    else NULL
    }
nwkillat <- function(program) {
    suspectlines <- which(grepl("<<[~>]>>", program))

    # This is slower than Hades, but I expect to see only 0-3
    #   lines in the suspectlines set
    for (i in suspectlines) {
        line <- strsplit(program[i], split='')
        nl <- length(line)
        state <- 0  # 0=ordirnay text, 1=inside a [[
        cstate <- integer(nl)  # 1= an ampersand to be removed
        for (j in 1:nl) {
            if (state==0) {
                if ((i+1 <nl) && line[i]=='@' && line[i+1]=='<' && 
                    line[i+2]=='<') cstate[i]=1
                if (i<nl && line[i]=='[' && line[i+1]=='[') state <-1
                }
            else {
                if (i<nl && line[i]==']' && line[i+1]==']') state<-0
                }
            }
        if (any(cstate)) program[i] <- paste(line[cstate==0], collapse='')
        }
    program
    }
notangle <- function(noweb, target='*', file.) {
    cname <- names(noweb)
    indx <- match(target, cname)
    if (is.na(indx)) {
        if (missing(target) && any(cname != '')) 
            target <- (cname[cname!=''])[1]
        else stop(paste("Code chunk", target, "not found in source file"))
        }
    program <- nwextract(noweb, target, prefix="")
    class(program) <-"notangle"
    program
    }

print.notangle <- function(x, file., ...) {
    if (missing(file.)) cat(x, sep='\n')
    else cat(x, sep='\n', file=file.)
    invisible(x)
    }
nwextract<- function(code, target, prefix="") {
    mycode <- code[names(code)==target]
    if (is.null(mycode)) stop("Program logic flaw 1")
    
    for (chunk in 1:length(mycode)) {
        ctemp <- mycode[[chunk]]
        if (length(ctemp$xref) ==0) temp <- ctemp$lines
        else {
            inclusions <- length(ctemp$xref)
            temp <- vector("list", 2*inclusions +1)
            for (i in 1:length(ctemp$xref))
                temp[[2*i]] <- nwextract(code, ctemp$xref[i], ctemp$indent[i])
            start <- c(1, ctemp$xindex+1) #start and end of non-inclusions
            end   <- c(ctemp$xindex-1, length(ctemp$lines))
            for (i in 1:length(start)) 
                if (start[i]<=end[i]) 
                    temp[[2*i -1]] <- ctemp$lines[start[i]:end[i]]
            temp <- unlist(temp)
            }
        mycode[[chunk]] <- ifelse(temp=="", "", paste(prefix, temp, sep=''))
        }
    as.vector(unlist(mycode))   #kill any names added to the vector
    }

tab.to.blank <- function(x, tabstop=8) {
    blanks <- rep(" ", tabstop)
    for (i in (tabstop-1):1) blanks[i] <- paste(blanks[i +0:1], collapse='')

    temp <- strsplit(x, '')
    linefix <- function(x) {
        n <- length(x)
        if (n==0) ""
        else {
            since.last.tab <- 1:n - cummax(ifelse(x=='\t', 1:n, 0))
            newx <- ifelse(x=='\t', blanks[1+ since.last.tab%%tabstop], x)
            paste(newx, collapse='')
            }
        }
    unlist(lapply(temp, linefix))
    }

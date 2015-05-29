##' gtTable extracts all the targets and their p-values from a GT object.
##'
##' Details follow.
##' 
##' @title Extract gtTable
##' @param object gt-object
##' @return gt-table
##' @author Maarten van Iterson, Sander Bervoets
##' @export
gtTable <- function(object)
{
  test <- function(set) object@functions$test(set, calculateP=TRUE)

  leaves <- t(sapply(1:size(object), function(i) test(i)))
            
  ##calculate zscores  
  zsc <- (leaves[,"S"]  - leaves[,"ES"]) / leaves[,"sdS"]
            
  ##association
  positive <- object@functions$positive()
	
  tbl <- data.frame(Pvalue=leaves[,"p"], Association=positive, Weights = weights(object), zscores=zsc)
  
  rownames(tbl) <- object@functions$cov.names(1:size(object))	
  tbl
}

##' toTable generates table from gtrun result list
##'
##' Details follow
##' @title toTable
##' @param results gtrun results 
##' @param method see p.adjust methods
##' @param alpha singificance level
##' @param level local or global
##' @param org either Hs or Mm are implemented yet
##' @return data.frame
##' @author Maarten van Iterson
##' @export
toTable <- function(results, method="BH", alpha=0.05, level=c("global", "local"), org=c("Hs", "Mm", "None"))
  {
    if(org == "Hs")
      require(org.Hs.eg.db)
    else  if(org == "Mm")
      require(org.Mm.eg.db)
    
    mirs <- results$mirs
    targets <- results$targets
    table <- data.frame()
    cnames <- c("mRNA", "miRNA", "Pvalue")
    level <- match.arg(level)
    
    if(length(method)==2)
      {
        method.global <- method[1]
        method.local <- method[2]
      }
    else
      method.global <- method.local <- method

    if(length(alpha) == 2)
      {
        alpha.global <- alpha[1]
        alpha.local <- alpha[2]
      }
    else
      alpha.global <- alpha.local <- alpha
    
    if(!is.data.frame(mirs))
      {
        cnames <- c("miRNA", "mRNA", "Pvalue")
        tmp <- mirs
        mirs <- targets
        targets <- tmp
      }
    mirs <- mirs[p.adjust(mirs$Pvalue, method=method.global) < alpha.global,]
    for(row in rownames(mirs))
      {
        trgts <- targets[[row]]
        if(level=="local")
          {
            if(sum(p.adjust(trgts$Pvalue, method=method.local) < alpha.local) > 0)
              trgts <- trgts[p.adjust(trgts$Pvalue, method=method.local) < alpha.local, ]              
          }
        if(nrow(trgts) == 0)
          next
        
        rows <- cbind(rownames(trgts),
                      row,
                      mirs[rownames(mirs) == row, "Pvalue"],
                      trgts
                      )
        table <- rbind(table, rows)           
      }
    
    colnames(table)[1:3] <- cnames
    rownames(table) <- 1:nrow(table)
    if(org == "Hs")
      table$Symbol <- unlist(mget(as.character(table$mRNA), org.Hs.egSYMBOL))
    else if( org == "Mm")
      table$Symbol <- unlist(mget(as.character(table$mRNA), org.Mm.egSYMBOL))
    table
  }

##' rungt is a wrapper for the gt() function
##'
##' Details follow.
##' 
##' @title Run the global test on the list of mirs
##' @param mirs list of mirs
##' @param X mRNA expression
##' @param Y miRNA expression
##' @param path path to database
##' @param dbName database name
##' @param tables prediction databases
##' @param numOverlapping number of at least overlapping targets between databases
##' @param top number of significant targets returned -1 is all 
##' @return list of microRNAs and targets
##' @author Maarten van Iterson, Sander Bervoets
##' @export
rungt <- function (mirs, X, Y, path, dbName, tables, numOverlapping, top = -1) 
{
    targets <- list()
    micrornas <- matrix(NA, nrow = length(mirs), ncol = 2)
    for (i in 1:length(mirs)) {
        ovl <- overlap(path, dbName, tables = tables, mirs[i], numOverlapping, full = FALSE)
        if (is.null(ovl)) 
            next
        sX <- X[which(rownames(X) %in% rownames(ovl)), , drop = FALSE]
        if(dim(sX)[1] == 0) {
                message("Skipping... No targets in data for this microRNA!")
                next
            }
        sy <- Y[which(rownames(Y) %in% mirs[i]), ]
        obj <- gt(sy, t(sX), directional = 1e-06)        
        tbl <- gtTable(obj)
        tbl <- tbl[order(tbl$Association, tbl$Pvalue), ]
        if (top == -1) 
            targets[[mirs[i]]] <- tbl
        else targets[[mirs[i]]] <- tbl[1:max(c(top, nrow(tbl))),]
        micrornas[i, 1] <- p.value(obj)
        micrornas[i, 2] <- size(obj)
    }
    colnames(micrornas) <- c("Pvalue", "targets")
    rownames(micrornas) <- mirs
    micrornas <- data.frame(na.omit(micrornas))
    micrornas <- micrornas[order(micrornas$Pvalue), ]
    list(mirs = micrornas, targets = targets)
}

##' reversed version rungt is a wrapper for the gt() function
##'
##' Details follow.
##' 
##' @title Run the global test on the list of mirs
##' @param targets list of targets
##' @param Y miRNA expression
##' @param X mRNA expression
##' @param A data.frame containing predicted microRNA target pairs
##' @return list of microRNAs and targets
##' @author Maarten van Iterson
##' @export
runrgt <- function (targets, Y, X, A) 
{
    mirs <- list()
    Targets <- matrix(NA, nrow = length(targets), ncol = 2)
    
    for (i in 1:length(targets)) {

        miRNAs <- subset(A, mRNA == targets[i])$miRNA
        
        sx <- X[which(rownames(X) == targets[i]),]
        sY <- Y[which(rownames(Y) %in% miRNAs), , drop = FALSE]

        if(nrow(sY) == 0)
          next
               
        obj <- gt(sx, t(sY), directional = 1e-06)
        
        tbl <- gtTable(obj)
        tbl <- tbl[order(tbl$Association, tbl$Pvalue), ]
        mirs[[targets[i]]] <- tbl
        Targets[i, 1] <- p.value(obj)
        Targets[i, 2] <- size(obj)
    }
    colnames(Targets) <- c("Pvalue", "mirs")
    rownames(Targets) <- targets
    Targets <- data.frame(na.omit(Targets))
    Targets <- Targets[order(Targets$Pvalue), ]
    list(targets = Targets, mirs = mirs)
}


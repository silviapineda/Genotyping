#' Plotting a tree extracted from a random forest
#' 
#' @param rforest A \code{randomForest} object
#' @param tr A \code{tree} object. Either rforest or tr must be input
#' @param k The index of the tree to be plotted
#' @param depth The depth of the tree to be plotted
#' @param main The title to put on the graph
#' @param ... Additional parameters to be passed to \code{text.tree}
#' @export
#' @examples
#' library(randomForest)
#' rforest <- randomForest(Species~., data=iris, ntree=20)
#' plot.getTree(rforest, k=3, depth=4)
plot.getTree <- function(rforest=NULL,tr=NULL,k=1, depth=0,main=NULL, ...){
  require(randomForest)
  if(is.null(rforest) && is.null(tr))stop('One of a random forest object or a tree object must be input')
  if(!is.null(rforest)){
    gTree <- getTree(rforest, k=k, labelVar=TRUE)
    x <- as.tree(gTree, rforest)
  } else {
    x <- tr
  }
  if(depth>0){
    x <- snip.depth(x,depth)
  }
  plot(x, type='uniform')
  text(x,split=FALSE,...)
  labelBG(x)
  labelYN(x)
  title(main=main)
}
Contact GitHub API Training Shop Blog About
© 2017 GitHub, Inc. Terms Privacy Security Status Help





#' Representative trees from ensembles
#' 
#' This package implements the concept of representative trees from ensemble
#' tree-based learners introduced by Banerjee, et al, (2012). Representative trees
#' are, in some sense, trees in the ensemble which are on average the "closest" to 
#' all the other trees in the ensemble. Several trees can be representative of the 
#' ensemble in general. This package currently implements the d2 metric of tree closeness
#' (close in prediction) from Banerjee, et al.
#' 
#' @name reprtree-package
#' @references M. Banerjee, Y. Ding and A-M Noone (2012) "Identifying representative
#'     trees from ensembles". Statistics in Medicine, 31(15):1601-1616.
#' @import randomForest tree
#' @docType package
#' @name reprtree
#' @examples
#' library(reprtree)
#' rforest <- randomForest(Species~., data=iris)
#' reptree <- ReprTree(rforest, iris, metric='d2')
#' plot(reptree, index=1)
NULL


#' Identifying and extracting representative trees from a random forest 
#' 
#' This function takes a random forest object and data to run predictions on
#' and identifies representative trees based on the d0, d1 or d2 metric defined
#' in Banerjee, et al (2012). Currently only the d2 metric is implemented, using
#' either a euclidean distance (for numeric predictions) or a mismatch distance 
#' (for categorical predictions). The average distance D(T) of each tree in the 
#' set of trees in computed, and trees with the lowest D(T) value are extracted
#' and formatted to be compatible with the \code{\link{tree}} class. Trees can then
#' be visualized using \code{plot.tree} and \code{text.tree}, or using a custom
#' plot function (to be defined)
#' 
#' @param rforest A randomForest object
#' @param newdata The data on which predictions will be computed
#' @param metric The metric to be used to evaluate distance between trees. Currently
#'    only the d2 metric is implemented
#' @return A list object containing representations of the representative trees
#'    conformable with the \code{tree} class. Names of the list give the indices
#'    of the representative trees in the set of trees. 
#' @import randomForest tree
#' @export
#' @references M. Banerjee, Y. Ding and A-M Noone (2012) "Identifying representative
#'     trees from ensembles". Statistics in Medicine, 31(15):1601-1616.
#' @examples
#' library(randomForest)
#' library(tree)
#' rforest <- randomForest(Species~., data=iris)
#' reptree <- ReprTree(rforest, iris, metric='d2')
ReprTree <- function(rforest, newdata, metric='d2'){
  if(metric!='d2') stop('invalid metric!')
  require(randomForest)
  print('Constructing distance matrix...')
  preds <- predict2(rforest, newdata=newdata, predict.all=T)
  preds.indiv <- preds$individual
  d <- dist.fn(t(preds.indiv), method=ifelse(rforest$type=='classification',
                                             'mismatch',
                                             'euclidean'))
  print('Finding representative trees...')
  D <- colMeans(d)
  index <- which(D==min(D))
  trees <- lapply(as.list(index), function(i) getTree(rforest, i, labelVar=TRUE))
  names(trees) <- as.character(index)
  trees <- lapply(trees, as.tree, rforest)
  out <- list(trees=trees,D = D)
  class(out) <- c('reprtree','list')
  return(out)
}



# Convert the result of a getTree call to a format compatible with tree
# 
# This function takes the results of a \code{randomForest::getTree} call and 
# converts the results to a form compatible with \code{tree}
# @param gTree The results of a call to \code{getTree}
# @param rforest The randomForest object 
# @return An object of class \code{tree}, which has a \code{frame} and sufficient
#     attributes to enable plotting
as.tree <- function(gTree,rforest,max.depth=3){
  if(is.numeric(gTree[,'split var'])) stop("labelVar=T required")
  bl <- matrix("", nrow=nrow(gTree), ncol=3)
  for(row in 1:nrow(gTree)){
    if(row==1){
      bl[row, 1:2] <- c('10','11')
      next
    }
    if(gTree[row,1]>0){
      bl[row,1:2] <- paste0(bl[which(gTree[,1:2]==row,arr.ind=T)], c('0','1'))
    } else {
      bl[row,3] <- bl[which(gTree[,1:2]==row, arr.ind=T)]
    }
  }
  bl <- data.frame(bl, stringsAsFactors=F); names(bl) <- c('left','right','terminal')
  fr <- list()
  fr$var <- as.character(gTree[,"split var"])
  fr$var[is.na(fr$var)] <- '<leaf>'
  fr$n <- fr$dev <- rep(0,length(fr$var))
  fr$yval <- gTree[,'prediction']
  
  # Need to work out split points based on classes of the splitting vars
  classes <- attributes(rforest$terms)$dataClasses
  blah <- data.frame(var=fr$var, splits=as.character(gTree[,'split point']), 
                     classes=classes[fr$var], stringsAsFactors=F)
  index <- which(blah$classes=='factor' & !is.na(blah$classes))
  blah$splits[index] <- sapply(blah$splits[index], factor.repr)  
  
  
  splits <- cbind(
    cutleft=paste0(ifelse(blah$classes=='factor' & !is.na(blah$classes),': ','<'),
                   blah$splits), 
    cutright=paste0(ifelse(blah$classes=='factor' & !is.na(blah$classes),
                           ': ','>'),
                    blah$splits))
  splits[fr$var=='<leaf>',] <- ""
  
  fr <- as.data.frame(fr, stringsAsFactors=F)
  fr$splits <- splits
  x <- ifelse(fr$var=='<leaf>', bl[,3], gsub('.{1}$', '', bl[,1]))
  if(rforest$type=='classification'){
    fr$yprob = matrix(1/length(rforest$classes),nrow=nrow(fr), ncol=length(rforest$classes))
  }
  row.names(fr) <- strtoi(x,2)
  fr <- fr[order(x),]
  
  newtr <- list()
  newtr$frame=fr
  attr(newtr,'xlevels') <- rforest$forest$xlevels
  if(rforest$type=='classification') attr(newtr,'ylevels') <- rforest$classes
  class(newtr) <- 'tree'
  return(newtr)
}

# Compute a distance matrix between rows of a data matrix
# 
# This function takes a matrix or a data.frame, and computes the distance
# between the between the rows of the data matrix. It extends the function 
# \code{\link{dist}} by adding a new metric defined by the proportion of 
# mismatches between two vectors.
# 
# @param x a numeric matrix or data frame
# @param method the distance measure to be used. This must be one of "mismatch",
#      "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#      Any unambiguous substring can be given
# @param ... additional arguments to be passed to the \code{dist} function
# @keywords dist
dist.fn <- function(x, method='mismatch',...){
  METHODS <- c("euclidean", "maximum", "manhattan", "canberra", 
               "binary", "minkowski", "mismatch")
  method <- pmatch(method, METHODS)
  if(is.na(method)) stop("invalid distance method")
  if(METHODS[method] !="mismatch"){
    z <- as.matrix(dist(x, method=METHODS[method], ...))
  } else {
    z = matrix(0, nrow=nrow(x), ncol=nrow(x))
    for(k in 1:(nrow(x)-1)){
      for (l in (k+1):nrow(x)){
        z[k,l] <- mean(x[k,]!=x[l,])
        z[l,k] <- z[k,l]
      }}
  }
  dimnames(z)  <- list(dimnames(x)[[1]],dimnames(x)[[1]])
  return(z)
}

# Convert integers to binary representation
# 
# @param x integer to be converted
# @param reverse Should the ordering be reversed
int2bin <- function(x, reverse=F){
  y <- intToBits(x)
  yy <- paste(sapply(strsplit(paste(rev(y)),""),`[[`,2),collapse="")
  out <- gsub('^[0]+','',yy)
  if(reverse){
    bl <- rev(unlist(strsplit(out,'')))
    out <- paste(bl, collapse='')
  }
  return(out)
}

# Represent factor splits using letters
# 
# @param x character representation of integer in "split point"
factor.repr <- function(x){
  x <- int2bin(as.integer(x), reverse=T)
  n <- nchar(x)
  paste(letters[1:n][unlist(strsplit(x,''))=='1'],collapse='')
}

# Alternative version of the predict function
# 
# @param object Object from which predictions will be derived
predict2 <- function(object, ...){
  UseMethod('predict2')
}

# Alternative version of predict.randomForest, from randomForest package
# 
# @description
# This function adapts the predict.randomForest function to accomodate the situation
# where factors can have NA as a legitimate level
# 
predict2.randomForest <- function (object, newdata, type = "response", norm.votes = TRUE, 
                                   predict.all = FALSE, proximity = FALSE, nodes = FALSE, cutoff,
                                   ...) 
{
  if (!inherits(object, "randomForest")) 
    stop("object not of class randomForest")
  if (is.null(object$forest)) 
    stop("No forest component in the object")
  out.type <- charmatch(tolower(type), c("response", "prob", 
                                         "vote", "class"))
  if (is.na(out.type)) 
    stop("type must be one of 'response', 'prob', 'vote'")
  if (out.type == 4) 
    out.type <- 1
  if (out.type != 1 && object$type == "regression") 
    stop("'prob' or 'vote' not meaningful for regression")
  if (out.type == 2) 
    norm.votes <- TRUE
  if (missing(newdata)) {
    p <- if (!is.null(object$na.action)) {
      napredict(object$na.action, object$predicted)
    }
    else {
      object$predicted
    }
    if (object$type == "regression") 
      return(p)
    if (proximity & is.null(object$proximity)) 
      warning("cannot return proximity without new data if random forest object does not already have proximity")
    if (out.type == 1) {
      if (proximity) {
        return(list(pred = p, proximity = object$proximity))
      }
      else return(p)
    }
    v <- object$votes
    if (!is.null(object$na.action)) 
      v <- napredict(object$na.action, v)
    if (norm.votes) {
      t1 <- t(apply(v, 1, function(x) {
        x/sum(x)
      }))
      class(t1) <- c(class(t1), "votes")
      if (proximity) 
        return(list(pred = t1, proximity = object$proximity))
      else return(t1)
    }
    else {
      if (proximity) 
        return(list(pred = v, proximity = object$proximity))
      else return(v)
    }
  }
  if (missing(cutoff)) {
    cutoff <- object$forest$cutoff
  }
  else {
    if (sum(cutoff) > 1 || sum(cutoff) < 0 || !all(cutoff > 
                                                   0) || length(cutoff) != length(object$classes)) {
      stop("Incorrect cutoff specified.")
    }
    if (!is.null(names(cutoff))) {
      if (!all(names(cutoff) %in% object$classes)) {
        stop("Wrong name(s) for cutoff")
      }
      cutoff <- cutoff[object$classes]
    }
  }
  if (object$type == "unsupervised") 
    stop("Can't predict unsupervised forest.")
  if (inherits(object, "randomForest.formula")) {
    newdata <- as.data.frame(newdata)
    rn <- row.names(newdata)
    Terms <- delete.response(object$terms)
    x <- model.frame(Terms, newdata, na.action = na.omit)
    keep <- match(row.names(x), rn)
  }
  else {
    if (is.null(dim(newdata))) 
      dim(newdata) <- c(1, length(newdata))
    x <- newdata
    if (nrow(x) == 0) 
      stop("newdata has 0 rows")
    if (any(is.na(x))) 
      stop("missing values in newdata")
    keep <- 1:nrow(x)
    rn <- rownames(x)
    if (is.null(rn)) 
      rn <- keep
  }
  vname <- if (is.null(dim(object$importance))) {
    names(object$importance)
  }
  else {
    rownames(object$importance)
  }
  if (is.null(colnames(x))) {
    if (ncol(x) != length(vname)) {
      stop("number of variables in newdata does not match that in the training data")
    }
  }
  else {
    if (any(!vname %in% colnames(x))) 
      stop("variables in the training data missing in newdata")
    x <- x[, vname, drop = FALSE]
  }
  if (is.data.frame(x)) {
    isFactor <- function(x) is.factor(x) & !is.ordered(x)
    xfactor <- which(sapply(x, isFactor))
    if (length(xfactor) > 0 && "xlevels" %in% names(object$forest)) {
      for (i in xfactor) {
        if (any(!levels(x[[i]]) %in% object$forest$xlevels[[i]])) 
          stop("New factor levels not present in the training data")
        excl <- NA
        if(any(is.na(levels(x[[i]])))) excl <- NULL
        x[[i]] <- factor(x[[i]], levels = levels(x[[i]])[match(levels(x[[i]]), 
                                                               object$forest$xlevels[[i]])],
                         exclude=excl)
      }
    }
    cat.new <- sapply(x, function(x) if (is.factor(x) && 
                                         !is.ordered(x)) 
      length(levels(x))
      else 1)
    if (!all(object$forest$ncat == cat.new)) 
      stop("Type of predictors in new data do not match that of the training data.")
  }
  mdim <- ncol(x)
  ntest <- nrow(x)
  ntree <- object$forest$ntree
  maxcat <- max(object$forest$ncat)
  nclass <- object$forest$nclass
  nrnodes <- object$forest$nrnodes
  op <- options(warn = -1)
  on.exit(options(op))
  x <- t(data.matrix(x))
  if (predict.all) {
    treepred <- if (object$type == "regression") {
      matrix(double(ntest * ntree), ncol = ntree)
    }
    else {
      matrix(integer(ntest * ntree), ncol = ntree)
    }
  }
  else {
    treepred <- numeric(ntest)
  }
  proxmatrix <- if (proximity) 
    matrix(0, ntest, ntest)
  else numeric(1)
  nodexts <- if (nodes) 
    integer(ntest * ntree)
  else integer(ntest)
  if (object$type == "regression") {
    if (!is.null(object$forest$treemap)) {
      object$forest$leftDaughter <- object$forest$treemap[, 
                                                          1, , drop = FALSE]
      object$forest$rightDaughter <- object$forest$treemap[, 
                                                           2, , drop = FALSE]
      object$forest$treemap <- NULL
    }
    keepIndex <- "ypred"
    if (predict.all) 
      keepIndex <- c(keepIndex, "treepred")
    if (proximity) 
      keepIndex <- c(keepIndex, "proximity")
    if (nodes) 
      keepIndex <- c(keepIndex, "nodexts")
    if (!is.integer(object$forest$leftDaughter)) 
      storage.mode(object$forest$leftDaughter) <- "integer"
    if (!is.integer(object$forest$rightDaughter)) 
      storage.mode(object$forest$rightDaughter) <- "integer"
    if (!is.integer(object$forest$nodestatus)) 
      storage.mode(object$forest$nodestatus) <- "integer"
    if (!is.double(object$forest$xbestsplit)) 
      storage.mode(object$forest$xbestsplit) <- "double"
    if (!is.double(object$forest$nodepred)) 
      storage.mode(object$forest$nodepred) <- "double"
    if (!is.integer(object$forest$bestvar)) 
      storage.mode(object$forest$bestvar) <- "integer"
    if (!is.integer(object$forest$ndbigtree)) 
      storage.mode(object$forest$ndbigtree) <- "integer"
    if (!is.integer(object$forest$ncat)) 
      storage.mode(object$forest$ncat) <- "integer"
    ans <- .C("regForest", as.double(x), ypred = double(ntest), 
              as.integer(mdim), as.integer(ntest), as.integer(ntree), 
              object$forest$leftDaughter, object$forest$rightDaughter, 
              object$forest$nodestatus, nrnodes, object$forest$xbestsplit, 
              object$forest$nodepred, object$forest$bestvar, object$forest$ndbigtree, 
              object$forest$ncat, as.integer(maxcat), as.integer(predict.all), 
              treepred = as.double(treepred), as.integer(proximity), 
              proximity = as.double(proxmatrix), nodes = as.integer(nodes), 
              nodexts = as.integer(nodexts), DUP = FALSE, PACKAGE = "randomForest")[keepIndex]
    yhat <- rep(NA, length(rn))
    names(yhat) <- rn
    if (!is.null(object$coefs)) {
      yhat[keep] <- object$coefs[1] + object$coefs[2] * 
        ans$ypred
    }
    else {
      yhat[keep] <- ans$ypred
    }
    if (predict.all) {
      treepred <- matrix(NA, length(rn), ntree, dimnames = list(rn, 
                                                                NULL))
      treepred[keep, ] <- ans$treepred
    }
    if (!proximity) {
      res <- if (predict.all) 
        list(aggregate = yhat, individual = treepred)
      else yhat
    }
    else {
      res <- list(predicted = yhat, proximity = structure(ans$proximity, 
                                                          dim = c(ntest, ntest), dimnames = list(rn, rn)))
    }
    if (nodes) {
      attr(res, "nodes") <- matrix(ans$nodexts, ntest, 
                                   ntree, dimnames = list(rn[keep], 1:ntree))
    }
  }
  else {
    countts <- matrix(0, ntest, nclass)
    t1 <- .C("classForest", mdim = as.integer(mdim), ntest = as.integer(ntest), 
             nclass = as.integer(object$forest$nclass), maxcat = as.integer(maxcat), 
             nrnodes = as.integer(nrnodes), jbt = as.integer(ntree), 
             xts = as.double(x), xbestsplit = as.double(object$forest$xbestsplit), 
             pid = object$forest$pid, cutoff = as.double(cutoff), 
             countts = as.double(countts), treemap = as.integer(aperm(object$forest$treemap, 
                                                                      c(2, 1, 3))), nodestatus = as.integer(object$forest$nodestatus), 
             cat = as.integer(object$forest$ncat), nodepred = as.integer(object$forest$nodepred), 
             treepred = as.integer(treepred), jet = as.integer(numeric(ntest)), 
             bestvar = as.integer(object$forest$bestvar), nodexts = as.integer(nodexts), 
             ndbigtree = as.integer(object$forest$ndbigtree), 
             predict.all = as.integer(predict.all), prox = as.integer(proximity), 
             proxmatrix = as.double(proxmatrix), nodes = as.integer(nodes), 
             DUP = FALSE, PACKAGE = "randomForest")
    if (out.type > 1) {
      out.class.votes <- t(matrix(t1$countts, nrow = nclass, 
                                  ncol = ntest))
      if (norm.votes) 
        out.class.votes <- sweep(out.class.votes, 1, 
                                 rowSums(out.class.votes), "/")
      z <- matrix(NA, length(rn), nclass, dimnames = list(rn, 
                                                          object$classes))
      z[keep, ] <- out.class.votes
      class(z) <- c(class(z), "votes")
      res <- z
    }
    else {
      out.class <- factor(rep(NA, length(rn)), levels = 1:length(object$classes), 
                          labels = object$classes)
      out.class[keep] <- object$classes[t1$jet]
      names(out.class)[keep] <- rn[keep]
      res <- out.class
    }
    if (predict.all) {
      treepred <- matrix(object$classes[t1$treepred], nrow = length(keep), 
                         dimnames = list(rn[keep], NULL))
      res <- list(aggregate = res, individual = treepred)
    }
    if (proximity) 
      res <- list(predicted = res, proximity = structure(t1$proxmatrix, 
                                                         dim = c(ntest, ntest), dimnames = list(rn[keep], 
                                                                                                rn[keep])))
    if (nodes) 
      attr(res, "nodes") <- matrix(t1$nodexts, ntest, ntree, 
                                   dimnames = list(rn[keep], 1:ntree))
  }
  res
}


#' Plotting representative trees
#' 
#' This function creates either a single plot or a panel of plots for 
#' visualizing representative trees from an ensemble
#' 
#' @param reptree An object of class \code{reprtree)}
#' @param index The index of the reprtree object you want to plot
#' @param all (logical) Do you want to create a panel of plots?
#' @param depth The maximum depth of the tree to be plotted. \code{depth=0} means the full tree.
#' @param adj 
#' @param main Title of plot (default is NULL)
#' @param ... additional arguments to pass to text.tree. In particular, suppress node labels using \code{label=NULL}
#' @export
#' @S3method plot reprtree
#' @section Details:
#' This plot function takes a \code{reprtree} object, and then plots a 
#' single representative tree or a sequence of representative trees (using \code{all=T}).
#' 
#' If only one tree needs to be visualized, the index of the reprtree object to
#' be visualized can be provided.
plot.reprtree <- function(reptree, all=F, index = ifelse(all,NULL, 1), depth=0,
                          main=NULL,adj = 0.5,  ...){
  require(plotrix)
  require(tree)
  if(!is(reptree,'reprtree')) stop('Wrong class!')
  reptree <- reptree[['trees']]
  n <- length(reptree)
  if(all){
    for(i in 1:n){
      plot.getTree(tr=reptree[[i]], depth=depth, main=main, adj=0.5, ...)
      #       tr <- reptree[[i]]
      #       if(depth>0){
      #         tr <- snip.depth(reptree[[i]], depth)
      #       } 
      #       plot(tr, type='uniform')
      #       text(tr,adj=adj,cex=0.7, split=F,...)
      #       labelBG(tr)
      #       labelYN(tr)
      #       title(main=main)
    }
  } else {
    plot.getTree(tr=reptree[[index]],depth=depth, main=main, adj=0.5,...)
    #     tr <- reptree[[index]]
    #     if(depth>0){
    #       tr <- snip.depth(tr, depth)
    #     }
    #     plot(tr, type='uniform') 
    #     text(tr,adj=adj,split=F, cex=0.7, digits=2, ...)
    #     labelBG(tr)
    #     labelYN(tr)
    #     title(main=main)
  }
}

labelBG <- function(tr){
  require(plotrix)
  charht <- par('cxy')[2L]
  xy <- tree:::treeco(tr,uniform=TRUE)
  nodes <- as.integer(row.names(tr$frame))
  left.child <- match(2*nodes, nodes)
  rows <- tree:::labels.tree(tr)[left.child]
  rows <- gsub('NA,','',rows)
  ind <- !is.na(left.child)
  boxed.labels(xy$x[ind],xy$y[ind]+0.5*charht, rows[ind] , border=F, bg='white',
               cex=0.8, xpad=0.5, ypad=1)
}
labelYN <- function(tr){
  charht <- par('cxy')[2L]
  xy <- tree:::treeco(tr, uniform=T)
  nodes <- as.integer(row.names(tr$frame))
  left.child <- match(2*nodes, nodes)
  ind <- !is.na(left.child)
  text(xy$x[ind]-0.1, xy$y[ind]-0.2*charht, '<< Y',cex=0.6, adj=1)
  text(xy$x[ind]+0.1, xy$y[ind]-0.2*charht, 'N >>', cex=0.6, adj=0)
}

plot.tree <- function (x, y = NULL, type = c("proportional", "uniform"), ...) 
{
  if (inherits(x, "singlenode")) 
    stop("cannot plot singlenode tree")
  if (!inherits(x, "tree")) 
    stop("not legitimate tree")
  type <- match.arg(type)
  uniform <- type == "uniform"
  dev <- dev.cur()
  if (dev == 1L) 
    dev <- 2L
  #assign(paste0("device", dev), uniform, envir = tree_env)
  invisible(treepl(tree:::treeco(x, uniform), node = as.integer(row.names(x$frame)), 
                   ...))
}


treepl <- function (xy, node, erase = FALSE, ...) 
{
  # Modified from tree:::treepl
  x <- xy$x
  y <- xy$y
  parent <- match((node%/%2L), node)
  sibling <- match(ifelse(node%%2L, node - 1L, node + 1L), 
                   node)
  xx <- rbind(x, x, x[sibling], x[sibling], NA)
  yy <- rbind(y, y[parent], y[parent], y[sibling], NA)
  if (any(erase)) {
    lines(c(xx[, erase]), c(yy[, erase]), col = par("bg"))
    return(x = x[!erase], y = y[!erase])
  }
  plot(range(x), c(min(y)-0.5, max(y)+0.5), type = "n", axes = FALSE, xlab = "", 
       ylab = "")
  text(x[1L], y[1L], "|", ...)
  lines(c(xx[, -1L]), c(yy[, -1L]), ...)
  list(x = x, y = y)
}

#' Annotate a Tree Plot
#' 
#' @S3method text tree
#' @description Add text to a tree plot. Modification of \code{tree:::text.tree} to add uniform type of tree.
#' 
#' @param x an object of class "tree"
#' @param splits logical. If \code{TRUE} the splits are labelled
#' @param label The name of column in the \code{frame} component of \code{x}, to be used to label the nodes. Can be \code{NULL} to suppress node-labelling
#' @param all logical. By default, only the leaves are labelled, but if true interior nodes are also labelled
#' @param pretty the manipulation used for split labels infolving attributes. See Details.
#' @param digits significant digits for numerical labels
#' @param adj,xpd,... graphical parameters such as \code{cex} and \code{font}
#' @param uniform logical. Is \code{plot(..., type='uniform')} used to create the tree plot?
#' 
#' @details
#' If pretty = 0 then the level names of a factor split attributes are used unchanged. 
#' If pretty = NULL, the levels are presented by a, b, ... z, 0 ... 5. 
#' If pretty is a positive integer, abbreviate is applied to the labels with that value for its argument minlength.
#' 
#' If the lettering is vertical (par srt = 90) and adj is not supplied it is adjusted appropriately.
#' 
#' @author Abhijit Dasgupta, modifying original code by B.D. Ripley
#' @seealso \link{\code{plot.tree}}
text.tree <- function (x, splits = TRUE, label = "yval", all = FALSE, pretty = NULL, 
                       digits = getOption("digits") - 3, adj = par("adj"), xpd = TRUE, uniform=T,
                       ...) 
{
  oldxpd <- par(xpd = xpd)
  on.exit(par(oldxpd))
  if (inherits(x, "singlenode")) 
    stop("cannot plot singlenode tree")
  if (!inherits(x, "tree")) 
    stop("not legitimate tree")
  frame <- x$frame
  column <- names(frame)
  if (!is.null(ylevels <- attr(x, "ylevels"))) 
    column <- c(column, ylevels)
  if (!is.null(label) && is.na(match(label, column))) 
    stop("label must be a column label of the frame component of the tree")
  charht <- par("cxy")[2L]
  if (!is.null(srt <- list(...)$srt) && srt == 90) {
    if (missing(adj)) 
      adj <- 0
    ladj <- 1 - adj
  }
  else ladj <- adj
  xy <- tree:::treeco(x, uniform=uniform)
  if (splits) {
    node <- as.integer(row.names(frame))
    left.child <- match(2 * node, node)
    rows <- tree:::labels.tree(x, pretty = pretty)[left.child]
    ind <- !is.na(rows)
    text(xy$x[ind], xy$y[ind] + 0.5 * charht, rows[ind], 
         adj = adj, ...)
  }
  if (!is.null(label)) {
    leaves <- if (all) 
      rep(TRUE, nrow(frame))
    else frame$var == "<leaf>"
    if (label == "yval" & !is.null(ylevels)) 
      stat <- as.character(frame$yval[leaves])
    else if (!is.null(ylevels) && !is.na(lev <- match(label, 
                                                      ylevels))) 
      stat <- format(signif(frame$yprob[leaves, lev], digits = digits))
    else stat <- format(signif(frame[leaves, label], digits = digits))
    if (!is.null(dim(stat)) && dim(stat)[2L] > 1) {
      if (length(dimnames(stat)[[2L]])) 
        stat[1L, ] <- paste(sep = ":", dimnames(stat)[[2L]], 
                            stat[1L, ])
      stat <- do.call("paste", c(list(sep = "\n"), split(stat, 
                                                         col(stat))))
    }
    text(xy$x[leaves], xy$y[leaves] - 0.5 * charht, labels = stat, 
         adj = ladj, ...)
  }
  invisible()
}

#' Snipping a tree to a particular depth
#' 
#' @param tr The tree to be snipped
#' @param depth The depth at which to snip the tree
#' @return A tree object of the specified depth
snip.depth <- function(tr, depth){
  require(tree)
  nodes <- sapply(as.integer(row.names(tr$frame)), int2bin)
  nodes.to.snip <- strtoi(nodes[nchar(nodes)==depth & tr$frame$var !='<leaf>'],2)
  return(snip.tree(tr, nodes.to.snip))
}



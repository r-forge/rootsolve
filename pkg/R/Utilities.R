### ============================================================================
### first some common functions
### ============================================================================

## =============================================================================
## Set the mfrow parameters and whether to "ask" for opening a new device
## =============================================================================

setplotpar <- function(nmdots, dots, nv, ask) {
    if (!any(match(nmdots, c("mfrow", "mfcol"), nomatch = 0))) {
      nc <- min(ceiling(sqrt(nv)),3)
      nr <- min(ceiling(nv/nc),3)
      mfrow <- c(nr, nc)
    }
    else if ("mfcol" %in% nmdots)
        mfrow <- rev(dots$mfcol)
    else mfrow <- dots$mfrow

    if (! is.null(mfrow)) {
      mf <- par(mfrow=mfrow)
    }

   ## interactively wait if there are remaining figures
    if (is.null(ask))
      ask <- prod(par("mfrow")) < length(which) && dev.interactive()

    return(ask)
}

## =============================================================================
## find a variable
## =============================================================================

selectstvar <- function (which,var) {

    if (!is.numeric(which)) {
        ln <- length(which)
        Select <- which(var %in% which)
        if (length(Select) != ln)
            stop("not all variables in 'which' are in 'x'")
    }
    else {
        Select <- which   # "Select now refers to the column number
        if (max(Select) > length(var))
            stop("index in 'which' too large")
        if (min(Select) < 1)
            stop("index in 'which' should be > 0")
    }
  return(Select)
}

### ============================================================================
### S3 methods
### ============================================================================

plot.steady1D <- function (x, which = NULL, grid = NULL, xyswap=FALSE, 
  ask = NULL, ...) {

# if x is vector, check if there is more than one species...  
    X <- x$y
    if (is.vector(X)) {
      nspec <- attributes(x)$nspec
      if (length(X)%%nspec != 0) 
        stop("length of 'x' should be a multiple of 'nspec' if x is a vector")
      x <- matrix(nc = nspec, data = X)
    } else x <- X   # only state variables
   
      
    if (is.null(which)) which <- 1:ncol(x)
    var <- colnames(x)
    if(is.null(var)) var <- 1:ncol(x)
    which <- selectstvar(which,var)
    
    if (is.null(grid)) 
       grid <- 1:nrow(x)
    if (length(grid) != nrow(x)) 
      stop("length of grid (x-axis) should be = number of rows in 'x$y'")
    
    np <- length(which)

    dots <- list(...)
    nmdots <- names(dots)

    # number of figures in a row and 
    # interactively wait if there are remaining figures
   
    ask <- setplotpar(nmdots, dots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    
    Main <- is.null(dots$main)

    labs <- (is.null(dots$xlab) && is.null(dots$ylab))
    xxlab <- if (is.null(dots$xlab))  "x"  else dots$xlab
    yylab <- if (is.null(dots$ylab))  ""   else dots$ylab 

    ## allow individual xlab and ylab (vectorized)
    xxlab <- rep(xxlab, length.out = np)
    yylab <- rep(yylab, length.out = np)

    for (i in which) {
        if (Main)
            dots$main <- colnames(x)[i]
            
        if (! xyswap) {            
          dots$xlab <- xxlab[i]
          dots$ylab <- yylab[i]
          do.call("plot", c(alist(grid, x[, i]), dots))
        } else {
          if (is.null(labs)){
            dots$ylab <- xxlab[i]
            dots$xlab <- yylab[i]
          } else {
            dots$xlab <- xxlab[i]
            dots$ylab <- yylab[i]
            dots$ylim <- rev(range(grid))    # y-axis reversed
          }
          do.call("plot", c(alist(x[, i], grid), dots))
        } 
    }
}


### ============================================================================

plot.steady2D <- function (x, which = NULL, image= TRUE, ask = NULL, ...) {

# if x is vector, check if there is more than one species...  
    X <- x$y
    out <- list()
    if (is.vector(X)) {
      nspec <- attributes(x)$nspec
      dimens <- attributes(x)$dimens
      if (length(X) - nspec*prod(dimens) != 0) 
        stop("length of 'x' should be = 'nspec' * prod(dimens) if x is a vector")
      x <- matrix(nc = nspec, data = X)
      
      for ( i in 1:nspec){
        istart <- (i-1)*prod(dimens) 
        out[[i]] <- matrix(nr=dimens[1], nc=dimens[2], data =
          X[(istart+1):(istart+prod(dimens))])
      }
    } else x <- X   # only state variables
   
      
    if (is.null(which)) which <- 1:ncol(x)
    var <- colnames(x)
    if(is.null(var)) var <- 1:ncol(x)
    which <- selectstvar(which,var)
    
    np <- length(which)

    dots <- list(...)
    nmdots <- names(dots)

    # number of figures in a row and 
    # interactively wait if there are remaining figures
   
    ask <- setplotpar(nmdots, dots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    
    Main <- is.null(dots$main)

    labs <- (is.null(dots$xlab) && is.null(dots$ylab))
    xxlab <- if (is.null(dots$xlab))  "x"  else dots$xlab
    yylab <- if (is.null(dots$ylab))  ""   else dots$ylab 
    if (image) {
    dots$col <- if (is.null(dots$col)) 
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100) else dots$col
    } else
    dots$color <- if (is.null(dots$color)) 
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")) else dots$color
    
    ## allow individual xlab and ylab (vectorized)
    xxlab <- rep(xxlab, length.out = np)
    yylab <- rep(yylab, length.out = np)

    for (i in which) {
        if (Main)
            dots$main <- colnames(x)[i]
            
        dots$xlab <- xxlab[i]
        dots$ylab <- yylab[i]
        if(image)
          do.call("image", c(alist(out[[i]]), dots)) else
          do.call("filled.contour", c(alist(out[[i]]), dots)) 
          
    }
}


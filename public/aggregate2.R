aggregate2 <-function(x,y,by,FUN,...,simplify = TRUE, drop = TRUE){
  if(!is.data.frame(x))
    x <- as.data.frame(x)
  if(!is.data.frame(y))
    y <- as.data.frame(y)
  FUN <- match.fun(FUN)
  if (NROW(x) == 0L|NROW(y) == 0L|NCOL(x) == 0L|NCOL(y) == 0L) 
    stop("no rows to aggregate")
  if (!is.list(by)) 
    stop("'by' must be a list")
  if (is.null(names(by)) && length(by)) 
    names(by) <- paste0("Group.", seq_along(by))
  else {
    nam <- names(by)
    ind <- which(!nzchar(nam))
    names(by)[ind] <- paste0("Group.", ind)
  }
  if (any(lengths(by) != NROW(x)) | any(lengths(by) != NROW(y))) 
    stop("arguments must have same length")
  yy <- as.data.frame(by, stringsAsFactors = FALSE)
  keep <- complete.cases(by)
  yy <- yy[keep, , drop = FALSE]
  x <- x[keep, , drop = FALSE]
  y <- y[keep, , drop = FALSE]
  nrx <- NROW(x)
  nry <- NROW(y)
  ident <- function(x) {
    y <- as.factor(x)
    l <- length(levels(y))
    s <- as.character(seq_len(l))
    n <- nchar(s)
    levels(y) <- paste0(strrep("0", n[l] - n), s)
    y
  }
  grp <- lapply(yy, ident)
  multi.yy <- !drop && ncol(yy)
  if (multi.yy) {
    lev <- lapply(grp, levels)
    yy <- as.list(yy)
    for (i in seq_along(yy)) {
      z <- yy[[i]][match(lev[[i]], grp[[i]])]
      if (is.factor(z) && any(keep <- is.na(z))) 
        z[keep] <- levels(z)[keep]
      yy[[i]] <- z
    }
    eGrid <- function(L) expand.grid(L, KEEP.OUT.ATTRS = FALSE, 
                                     stringsAsFactors = FALSE)
    yy <- eGrid(yy)
  }
  grp <- if (ncol(yy)) {
    names(grp) <- NULL
    do.call(paste, c(rev(grp), list(sep = ".")))
  }
  else integer(nrx)
  if (multi.yy) {
    lev <- as.list(eGrid(lev))
    names(lev) <- NULL
    lev <- do.call(paste, c(rev(lev), list(sep = ".")))
  }
  else yy <- yy[match(sort(unique(grp)), grp, 0L), , drop = FALSE]
  
  
  z <- lapply2(x,y, function(e,f) {
    ans <- lapply2(X = unname(split(e, grp)),Y = unname(split(f, grp)), FUN = FUN, ...)
    if (simplify && length(len <- unique(lengths(ans))) == 
        1L) {
      if (len == 1L) {
        cl <- lapply(ans, oldClass)
        cl1 <- cl[[1L]]
        ans <- if (!is.null(cl1) && all(vapply(cl, identical, 
                                               NA, y = cl1))) 
          do.call(c, ans)
        else unlist(ans, recursive = FALSE, use.names = FALSE)
      }
      else if (len > 1L) 
        ans <- matrix(unlist(ans, recursive = FALSE, 
                             use.names = FALSE), ncol = len, byrow = TRUE, 
                      dimnames = if (!is.null(nms <- names(ans[[1L]]))) 
                        list(NULL, nms))
    }
    ans
  })
  len <- length(yy)
  if (multi.yy){
    keep <- match(lev, sort(unique(grp)))
    for (i in seq_along(z)) yy[[len + i]] <- if (is.matrix(z[[i]])) 
      z[[i]][keep, , drop = FALSE]
    else z[[i]][keep]
  }
  else for (i in seq_along(z)) yy[[len + i]] <- z[[i]]
  names(yy) <- c(names(by), names(x))
  row.names(yy) <- NULL
  yy
}
  
lapply2 <- function(X,Y,FUN,...){
  FUN <- match.fun(FUN)
  if (!is.vector(X) || is.object(X)) 
    X <- as.list(X)
  if (!is.vector(Y) || is.object(Y)) 
    Y <- as.list(Y)
  if(!length(X)==length(Y))
    stop("arguments must have same length")
  Z <- list()
  for (i in seq_along(X)) {
    Z[[i]] <- FUN(X[[i]],Y[[i]],...)
  }
  Z
}

"ade4toR" <-
function (fictab, ficcolnames = NULL, ficrownames = NULL) 
{
    if (!file.exists(fictab)) 
        stop(paste("file", fictab, "not found"))
    if (!is.null(ficrownames) && !file.exists(ficrownames)) 
        stop(paste("file", ficrownames, "not found"))
    if (!is.null(ficcolnames) && !file.exists(ficcolnames)) 
        stop(paste("file", ficcolnames, "not found"))
    w <- read.table(fictab, h = FALSE)
    nl <- nrow(w)
    nc <- ncol(w)
    if (!is.null(ficcolnames)) 
        provicol <- as.character((read.table(ficcolnames, h = FALSE))$V1)
    else provicol <- as.character(1:nc)
    if ((length(provicol)) != nc) {
        stop(paste("Non convenient row number in file", ficcolnames, 
            "- Expected:", nc, "- Input:", length(provicol)))
    }
    if (is.null(ficcolnames)) 
        names(w) <- paste("v", provicol, sep = "")
    else names(w) <- provicol
    if (!is.null(ficrownames)) 
        provirow <- as.character((read.table(ficrownames, h = FALSE))$V1)
    else provirow <- as.character(1:nl)
    if ((length(provirow)) != nl) {
        stop(paste("Non convenient row number in file", ficrownames, 
            "- Expected:", nl, "- Input:", length(provirow)))
    }
    row.names(w) <- provirow
    return(w)
}

"bicenter.wt" <-
function (X, row.wt = rep(1, nrow(X)), col.wt = rep(1, ncol(X))) 
{
    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    if (length(row.wt) != n) 
        stop("length of row.wt must equal the number of rows in x")
    if (any(row.wt < 0) || (sr <- sum(row.wt)) == 0) 
        stop("weights must be non-negative and not all zero")
    row.wt <- row.wt/sr
    if (length(col.wt) != p) 
        stop("length of col.wt must equal the number of columns in x")
    if (any(col.wt < 0) || (st <- sum(col.wt)) == 0) 
        stop("weights must be non-negative and not all zero")
    col.wt <- col.wt/st
    row.mean <- apply(row.wt * X, 2, sum)
    col.mean <- apply(col.wt * t(X), 2, sum)
    col.mean <- col.mean - sum(row.mean * col.wt)
    X <- sweep(X, 2, row.mean)
    X <- t(sweep(t(X), 2, col.mean))
    return(X)
}

"prep.fuzzy.var" <-
function (df, col.blocks, row.w = rep(1, nrow(df))) 
{
    if (!is.data.frame(df)) 
        stop("data.frame expected")
    if (!is.null(row.w)) {
        if (length(row.w) != nrow(df)) 
            stop("non convenient dimension")
    }
    if (sum(col.blocks) != ncol(df)) {
        stop("non convenient data in col.blocks")
    }
    if (is.null(row.w)) 
        row.w <- rep(1, nrow(df))/nrow(df)
    row.w <- row.w/sum(row.w)
    if (is.null(names(col.blocks))) {
        names(col.blocks) <- paste("FV", as.character(1:length(col.blocks)), 
            sep = "")
    }
    f1 <- function(x) {
        a <- sum(x)
        if (is.na(a)) 
            return(rep(0, length(x)))
        if (a == 0) 
            return(rep(0, length(x)))
        return(x/a)
    }
    k2 <- 0
    col.w <- rep(1, ncol(df))
    for (k in 1:(length(col.blocks))) {
        k1 <- k2 + 1
        k2 <- k2 + col.blocks[k]
        X <- df[, k1:k2]
        X <- t(apply(X, 1, f1))
        X.marge <- apply(X, 1, sum)
        X.marge <- X.marge * row.w
        X.marge <- X.marge/sum(X.marge)
        X.mean <- apply(X * X.marge, 2, sum)
        nr <- sum(X.marge == 0)
        if (nr > 0) {
            nc <- col.blocks[k]
            X[X.marge == 0, ] <- rep(X.mean, rep(nr, nc))
            cat(nr, "missing data found in block", k, "\n")
        }
        df[, k1:k2] <- X
        col.w[k1:k2] <- X.mean
    }
    attr(df, "col.blocks") <- col.blocks
    attr(df, "row.w") <- row.w
    attr(df, "col.freq") <- col.w
    col.num <- factor(rep((1:length(col.blocks)), col.blocks))
    attr(df, "col.num") <- col.num
    return(df)
}
"Rtoade4" <-
function (x) 
{
    if (!is.data.frame(x)) 
        stop("x is not a data.frame")
    nombase <- deparse(substitute(x))
    if (all(unlist(lapply(x, is.factor)))) {
        z <- matrix(0, nrow(x), ncol(x))
        for (j in 1:(ncol(x))) {
            toto <- x[, j]
            z[, j] <- unlist(lapply(toto, function(x, fac) which(x == 
                levels(fac)), fac = toto))
        }
        nomfic <- paste(nombase, ".txt", sep = "")
        write.table(z, file = nomfic, quote = FALSE, sep = "    ", 
            eol = "\n", na = "-999", row.names = FALSE, col.names = FALSE, 
            qmethod = c("escape", "double"))
        cat("File creation", nomfic, "\n")
        if (!is.null(attr(x, "names"))) {
            y <- attr(x, "names")
            nomfic <- paste(nombase, "_var_lab.txt", sep = "")
            write.table(y, file = nomfic, quote = FALSE, sep = "    ", 
                eol = "\n", na = "999", row.names = FALSE, col.names = FALSE, 
                qmethod = c("escape", "double"))
            cat("File creation", nomfic, "\n")
        }
        nommoda <- NULL
        for (j in 1:(ncol(x))) {
            toto <- x[, j]
            nommoda <- c(nommoda, levels(toto))
        }
        nomfic <- paste(nombase, "_moda_lab.txt", sep = "")
        write.table(y, file = nomfic, quote = FALSE, sep = "    ", 
            eol = "\n", na = "999", row.names = FALSE, col.names = FALSE, 
            qmethod = c("escape", "double"))
        cat("File creation", nomfic, "\n")
        return(invisible())
    }
    nomfic <- paste(nombase, ".txt", sep = "")
    write.table(x, file = nomfic, quote = FALSE, sep = "    ", eol = "\n", 
        na = "-999", row.names = FALSE, col.names = FALSE, qmethod = c("escape", 
            "double"))
    cat("File creation", nomfic, "\n")
    if (!is.null(attr(x, "names"))) {
        y <- attr(x, "names")
        nomfic <- paste(nombase, "_col_lab.txt", sep = "")
        write.table(y, file = nomfic, quote = FALSE, sep = "    ", 
            eol = "\n", na = "-999", row.names = FALSE, col.names = FALSE, 
            qmethod = c("escape", "double"))
        cat("File creation", nomfic, "\n")
    }
    if (!is.null(attr(x, "row.names"))) {
        y <- attr(x, "row.names")
        nomfic <- paste(nombase, "_row_lab.txt", sep = "")
        write.table(y, file = nomfic, quote = FALSE, sep = "    ", 
            eol = "\n", na = "-999", row.names = FALSE, col.names = FALSE, 
            qmethod = c("escape", "double"))
        cat("File creation", nomfic, "\n")
    }
    if (!is.null(attr(x, "col.blocks"))) {
        y <- as.vector(attr(x, "col.blocks"))
        nomfic <- paste(nombase, "_col_bloc.txt", sep = "")
        write.table(y, file = nomfic, quote = FALSE, sep = "    ", 
            eol = "\n", na = "-999", row.names = FALSE, col.names = FALSE, 
            qmethod = c("escape", "double"))
        cat("File creation", nomfic, "\n")
        y <- names(attr(x, "col.blocks"))
        nomfic <- paste(nombase, "_col_bloc_lab.txt", sep = "")
        write.table(y, file = nomfic, quote = FALSE, sep = "    ", 
            eol = "\n", na = "-999", row.names = FALSE, col.names = FALSE, 
            qmethod = c("escape", "double"))
        cat("File creation", nomfic, "\n")
    }
}
"scalewt" <-
function (X, wt = rep(1, nrow(X)), center = TRUE, scale = TRUE) 
{
    X <- as.matrix(X)
    n <- nrow(X)
    if (length(wt) != n) 
        stop("length of wt must equal the number of rows in x")
    if (any(wt < 0) || (s <- sum(wt)) == 0) 
        stop("weights must be non-negative and not all zero")
    wt <- wt/s
    center <- if (center) 
        apply(wt * X, 2, sum)
    else 0
    X <- sweep(X, 2, center)
    norm <- apply(X * X * wt, 2, sum)
    norm[norm <= 1e-07 * max(norm)] <- 1
    if (scale) 
        X <- sweep(X, 2, sqrt(norm), "/")
    return(X)
}
"unique.df" <-
function (x) 
{
    x <- data.frame(x)
    lig <- nrow(x)
    col <- ncol(x)
    w <- unlist(x[1])
    for (j in 2:col) {
        w <- paste(w, x[, j], sep = "")
    }
    w <- factor(w, unique(w))
    levels(w) <- 1:length(unique(w))
    select <- match(1:length(w), w)[1:nlevels(w)]
    x <- x[select, ]
    attr(x, "factor") <- w
    attr(x, "len.class") <- as.vector(table(w))
    return(x)
}

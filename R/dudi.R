"as.dudi" <-
function (df, col.w, row.w, scannf, nf, call, type, tol = 1e-07, 
    full = FALSE) 
{
    if (!is.data.frame(df)) 
        stop("data.frame expected")
    lig <- nrow(df)
    col <- ncol(df)
    if (length(col.w) != col) 
        stop("Non convenient col weights")
    if (length(row.w) != lig) 
        stop("Non convenient row weights")
    if (any(col.w) < 0) 
        stop("col weight < 0")
    if (any(row.w) < 0) 
        stop("row weight < 0")
    if (full) 
        scannf <- FALSE
    res <- list(tab = df, cw = col.w, lw = row.w)
    df <- as.matrix(df)
    df <- df * sqrt(row.w)
    df <- sweep(df, 2, sqrt(col.w), "*")
    svd1 <- svd(df)
    eig <- svd1$d^2
    rank <- sum((eig/eig[1]) > tol)
    if (scannf) {
        barplot(eig[1:rank])
        cat("Select the number of axes: ")
        nf <- as.integer(readLines(n = 1))
    }
    if (nf <= 0) 
        nf <- 2
    if (nf > rank) 
        nf <- rank
    if (full) 
        nf <- rank
    res$eig <- eig[1:rank]
    res$rank <- rank
    res$nf <- nf
    col.w[which(col.w == 0)] <- 1
    col.w <- 1/sqrt(col.w)
    auxi <- data.frame(svd1$v[, 1:nf] * col.w)
    names(auxi) <- paste("CS", (1:nf), sep = "")
    row.names(auxi) <- names(res$tab)
    res$c1 <- auxi
    row.w[which(row.w == 0)] <- 1
    row.w <- 1/sqrt(row.w)
    auxi <- data.frame(svd1$u[, 1:nf] * row.w)
    names(auxi) <- paste("RS", (1:nf), sep = "")
    row.names(auxi) <- row.names(res$tab)
    res$l1 <- auxi
    w <- matrix(svd1$d[1:nf], col, nf, byr = TRUE)
    auxi <- data.frame(as.matrix(res$c1) * w)
    names(auxi) <- paste("Comp", (1:nf), sep = "")
    row.names(auxi) <- names(res$tab)
    res$co <- auxi
    w <- matrix(svd1$d[1:nf], lig, nf, byr = TRUE)
    auxi <- data.frame(as.matrix(res$l1) * w)
    names(auxi) <- paste("Axis", (1:nf), sep = "")
    row.names(auxi) <- row.names(res$tab)
    res$li <- auxi
    res$call <- call
    class(res) <- c(type, "dudi")
    return(res)
}


"dudi.dec" <-
function (df, eff, scannf = TRUE, nf = 2) 
{
    if (!is.data.frame(df)) 
        stop("data.frame expected")
    lig <- nrow(df)
    col <- ncol(df)
    if (any(df < 0)) 
        stop("negative entries in table")
    if ((N <- sum(df)) == 0) 
        stop("all frequencies are zero")
    if (length(eff) != lig) 
        stop("non convenient dimension")
    if (any(eff) <= 0) 
        stop("non convenient vector eff")
    rtot <- sum(eff)
    row.w <- eff/rtot
    col.w <- apply(df, 2, sum)
    col.w <- col.w/rtot
    df <- sweep(df, 1, eff, "/")
    df <- sweep(df, 2, col.w, "/") - 1
    if (any(is.na(df))) {
        fun1 <- function(x) {
            if (is.na(x)) 
                return(0)
            else return(x)
        }
        df <- apply(df, c(1, 2), fun1)
        df <- data.frame(df)
    }
    X <- as.dudi(df, col.w, row.w, scannf = scannf, nf = nf, 
        call = match.call(), type = "dec")
    X$R <- rtot
    return(X)
}
"dudi.fca" <-
function (df, scannf = TRUE, nf = 2) 
{
    if (!is.data.frame(df)) 
        stop("data.frame expected")
    if (is.null(attr(df, "col.blocks"))) 
        stop("attribute 'col.blocks' expected for df")
    if (is.null(attr(df, "row.w"))) 
        stop("attribute 'row.w' expected for df")
    bloc <- attr(df, "col.blocks")
    row.w <- attr(df, "row.w")
    indica <- attr(df, "col.num")
    nvar <- length(bloc)
    col.w <- apply(df * row.w, 2, sum)
    df <- sweep(df, 2, col.w, "/") - 1
    col.w <- col.w/length(bloc)
    X <- as.dudi(df, col.w, row.w, scannf = scannf, nf = nf, 
        call = match.call(), type = "fca")
    rcor <- matrix(0, nvar, X$nf)
    rcor <- row(rcor) + 0 + (0+1i) * col(rcor)
    floc <- function(x) {
        i <- Re(x)
        j <- Im(x)
        if (i == 1) 
            k1 <- 0
        else k1 <- cumsum(bloc)[i - 1]
        k2 <- k1 + bloc[i]
        k1 <- k1 + 1
        z <- X$co[k1:k2, j]
        poicla <- X$cw[k1:k2] * nvar
        return(sum(poicla * z * z))
    }
    rcor <- apply(rcor, c(1, 2), floc)
    rcor <- data.frame(rcor)
    row.names(rcor) <- names(bloc)
    names(rcor) <- names(X$l1)
    X$cr <- rcor
    X$blo <- bloc
    X$indica <- indica
    return(X)
}
"dudi.mix" <-
function (df, add.square = FALSE, scannf = TRUE, nf = 2) 
{
    if (!is.data.frame(df)) 
        stop("data.frame expected")
    row.w <- rep(1, nrow(df))/nrow(df)
    acm.util <- function(cl) {
        n <- length(cl)
        cl <- as.factor(cl)
        x <- matrix(0, n, length(levels(cl)))
        x[(1:n) + n * (unclass(cl) - 1)] <- 1
        dimnames(x) <- list(names(cl), as.character(levels(cl)))
        data.frame(x)
    }
    f1 <- function(v) {
        moy <- sum(v)/length(v)
        v <- v - moy
        et <- sqrt(sum(v * v)/length(v))
        return(v/et)
    }
    df <- data.frame(df)
    nc <- ncol(df)
    nl <- nrow(df)
    if (any(is.na(df))) 
        stop("na entries in table")
    index <- rep("", nc)
    for (j in 1:nc) {
        w1 <- "q"
        if (is.factor(df[, j])) 
            w1 <- "f"
        if (is.ordered(df[, j])) 
            w1 <- "o"
        index[j] <- w1
    }
    res <- matrix(0, nl, 1)
    provinames <- "0"
    col.w <- NULL
    col.assign <- NULL
    k <- 0
    for (j in 1:nc) {
        if (index[j] == "q") {
            if (!add.square) {
                res <- cbind(res, f1(df[, j]))
                provinames <- c(provinames, names(df)[j])
                col.w <- c(col.w, 1)
                k <- k + 1
                col.assign <- c(col.assign, k)
            }
            else {
                w <- df[, j]
                deg.poly <- 2
                w <- sqrt(nl - 1) * poly(w, deg.poly)
                cha <- paste(names(df)[j], c(".L", ".Q"), sep = "")
                res <- cbind(res, as.matrix(w))
                provinames <- c(provinames, cha)
                col.w <- c(col.w, rep(1, deg.poly))
                k <- k + 1
                col.assign <- c(col.assign, rep(k, deg.poly))
            }
        }
        else if (index[j] == "o") {
            w <- as.numeric(df[, j])
            deg.poly <- min(nlevels(df[, j]) - 1, 2)
            w <- sqrt(nl - 1) * poly(w, deg.poly)
            if (deg.poly == 1) 
                cha <- names(df)[j]
            else cha <- paste(names(df)[j], c(".L", ".Q"), sep = "")
            res <- cbind(res, as.matrix(w))
            provinames <- c(provinames, cha)
            col.w <- c(col.w, rep(1, deg.poly))
            k <- k + 1
            col.assign <- c(col.assign, rep(k, deg.poly))
        }
        else if (index[j] == "f") {
            w <- acm.util(factor(df[, j]))
            cha <- paste(substr(names(df)[j], 1, 5), ".", names(w), 
                sep = "")
            col.w.provi <- drop(row.w %*% as.matrix(w))
            w <- t(t(w)/col.w.provi) - 1
            col.w <- c(col.w, col.w.provi)
            res <- cbind(res, w)
            provinames <- c(provinames, cha)
            k <- k + 1
            col.assign <- c(col.assign, rep(k, length(cha)))
        }
    }
    res <- data.frame(res)
    names(res) <- make.names(provinames, unique = TRUE)
    res <- res[, -1]
    names(col.w) <- provinames[-1]
    X <- as.dudi(res, col.w, row.w, scannf = scannf, nf = nf, 
        call = match.call(), type = "mix")
    X$assign <- factor(col.assign)
    X$index <- factor(index)
    rcor <- matrix(0, nc, X$nf)
    rcor <- row(rcor) + 0 + (0+1i) * col(rcor)
    floc <- function(x) {
        i <- Re(x)
        j <- Im(x)
        if (index[i] == "q") {
            if (sum(col.assign == i)) {
                w <- X$l1[, j] * X$lw * X$tab[, col.assign == 
                  i]
                return(sum(w)^2)
            }
            else {
                w <- X$lw * X$l1[, j]
                w <- X$tab[, col.assign == i] * w
                w <- apply(w, 2, sum)
                return(sum(w^2))
            }
        }
        else if (index[i] == "o") {
            w <- X$lw * X$l1[, j]
            w <- X$tab[, col.assign == i] * w
            w <- apply(w, 2, sum)
            return(sum(w^2))
        }
        else if (index[i] == "f") {
            x <- X$l1[, j] * X$lw
            qual <- df[, i]
            poicla <- unlist(tapply(X$lw, qual, sum))
            z <- unlist(tapply(x, qual, sum))/poicla
            return(sum(poicla * z * z))
        }
        else return(NA)
    }
    rcor <- apply(rcor, c(1, 2), floc)
    rcor <- data.frame(rcor)
    row.names(rcor) <- names(df)
    names(rcor) <- names(X$l1)
    X$cr <- rcor
    X
}
"dudi.nsc" <-
function (df, scannf = TRUE, nf = 2) 
{
    df <- data.frame(df)
    lig <- nrow(df)
    col <- ncol(df)
    if (any(df < 0)) 
        stop("negative entries in table")
    if ((N <- sum(df)) == 0) 
        stop("all frequencies are zero")
    row.w <- apply(df, 1, sum)/N
    col.w <- apply(df, 2, sum)/N
    df <- t(apply(df, 1, function(x) if (sum(x) == 0) 
        col.w
    else x/sum(x)))
    df <- sweep(df, 2, col.w)
    df <- data.frame(col * df)
    X <- as.dudi(df, rep(1, col)/col, row.w, scannf = scannf, 
        nf = nf, call = match.call(), type = "nsc")
    X$N <- N
    return(X)
}
"dudi.pca" <-
function (df, row.w = rep(1, nrow(df))/nrow(df), col.w = rep(1, 
    ncol(df)), center = TRUE, scale = TRUE, scannf = TRUE, nf = 2) 
{
    df <- data.frame(df)
    nc <- ncol(df)
    if (any(is.na(df))) 
        stop("na entries in table")
    f1 <- function(v) sum(v * row.w)/sum(row.w)
    f2 <- function(v) sqrt(sum(v * v * row.w)/sum(row.w))
    if (is.logical(center)) {
        if (center) {
            center <- apply(df, 2, f1)
            df <- sweep(df, 2, center)
        }
        else center <- rep(0, nc)
    }
    else if (is.numeric(center) && (length(center) == nc)) 
        df <- sweep(df, 2, center)
    else stop("Non convenient selection for center")
    if (scale) {
        norm <- apply(df, 2, f2)
        norm[norm < 1e-08] <- 1
        df <- sweep(df, 2, norm, "/")
    }
    else norm <- rep(1, nc)
    X <- as.dudi(df, col.w, row.w, scannf = scannf, nf = nf, 
        call = match.call(), type = "pca")
    X$cent <- center
    X$norm <- norm
    X
}


"dudi.pco" <-
function (d, row.w = "uniform", scannf = TRUE, nf = 2, full = FALSE, 
    tol = 1e-07) 
{
    if (!inherits(d, "dist")) 
        stop("Distance matrix expected")
    if (full) 
        scannf <- FALSE
    distmat <- dist2mat(d)
    n <- ncol(distmat)
    rownames <- attr(d, "Labels")
    if (any(is.na(d))) 
        stop("missing value in d")
    if (is.null(rownames)) 
        rownames <- as.character(1:n)
    if (row.w == "uniform") {
        row.w <- rep(1, n)
    }
    else {
        if (length(row.w) != n) 
            stop("Non convenient length(row.w)")
        if (any(row.w < 0)) 
            stop("Non convenient row.w (p<0)")
        if (any(row.w == 0)) 
            stop("Non convenient row.w (p=0)")
    }
    row.w <- row.w/sum(row.w)
    delta <- -0.5 * bicenter.wt(distmat * distmat, row.wt = row.w, 
        col.wt = row.w)
    wsqrt <- sqrt(row.w)
    delta <- delta * wsqrt
    delta <- t(t(delta) * wsqrt)
    eig <- La.eigen(delta, symmetric = TRUE)
    lambda <- eig$values
    w0 <- lambda[n]/lambda[1]
    if (w0 < -tol) 
        warning("Non euclidean distance")
    r <- sum(lambda > (lambda[1] * tol))
    if (scannf) {
        barplot(lambda)
        cat("Select the number of axes: ")
        nf <- as.integer(readLines(n = 1))
    }
    if (nf <= 0) 
        nf <- 2
    if (nf > r) 
        nf <- r
    if (full) 
        nf <- r
    res <- list()
    res$eig <- lambda[1:r]
    res$rank <- r
    res$nf <- nf
    res$cw <- rep(1, r)
    w <- t(t(eig$vectors[, 1:r]) * sqrt(lambda[1:r]))/wsqrt
    w <- data.frame(w)
    names(w) <- paste("A", 1:r, sep = "")
    row.names(w) <- rownames
    res$tab <- w
    res$li <- w[, 1:nf]
    w <- t(t(eig$vectors[, 1:nf])/wsqrt)
    w <- data.frame(w)
    names(w) <- paste("RS", 1:nf, sep = "")
    row.names(w) <- rownames
    res$l1 <- w
    w <- data.frame(diag(1, r))
    names(w) <- paste("CS", (1:r), sep = "")
    row.names(w) <- names(res$tab)
    res$c1 <- w[, 1:nf]
    w <- data.frame(matrix(0, r, nf))
    w[1:nf, 1:nf] <- diag(sqrt(lambda[1:nf]))
    names(w) <- paste("Comp", (1:nf), sep = "")
    row.names(w) <- names(res$tab)
    res$co <- w
    res$lw <- row.w
    res$call <- match.call()
    class(res) <- c("pco", "dudi")
    return(res)
}

"inertia.dudi" <-
function (dudi, row.inertia = FALSE, col.inertia = FALSE) 
{
    if (!inherits(dudi, "dudi")) 
        stop("Object of class 'dudi' expected")
    app <- function(x) {
        if (is.na(x)) 
            return(x)
        if (is.infinite(x)) 
            return(NA)
        if ((ceiling(x) - x) > (x - floor(x))) 
            return(floor(x))
        else return(ceiling(x))
    }
    nf <- dudi$nf
    inertia <- dudi$eig
    cum <- cumsum(inertia)
    ratio <- cum/sum(inertia)
    TOT <- cbind.data.frame(inertia, cum, ratio)
    listing <- list(TOT = TOT)
    if (row.inertia) {
        w <- dudi$tab * sqrt(dudi$lw)
        w <- sweep(w, 2, sqrt(dudi$cw), "*")
        w <- w * w
        con.tra <- apply(w, 1, sum)/sum(w)
        w <- dudi$li * dudi$li * dudi$lw
        w <- sweep(w, 2, dudi$eig[1:nf], "/")
        listing$row.abs <- apply(10000 * w, c(1, 2), app)
        w <- dudi$tab
        w <- sweep(w, 2, sqrt(dudi$cw), "*")
        d2 <- apply(w * w, 1, sum)
        w <- dudi$li * dudi$li
        w <- sweep(w, 1, d2, "/")
        w <- w * sign(dudi$li)
        names(w) <- names(dudi$li)
        w <- cbind.data.frame(w, con.tra)
        listing$row.rel <- apply(10000 * w, c(1, 2), app)
        w <- dudi$li * dudi$li
        w <- sweep(w, 1, d2, "/")
        w <- data.frame(t(apply(w, 1, cumsum)))
        names(w) <- names(dudi$li)
        remain <- 1 - w[, ncol(w)]
        w <- cbind.data.frame(w, remain)
        listing$row.cum <- apply(10000 * w, c(1, 2), app)
    }
    if (col.inertia) {
        w <- dudi$tab * sqrt(dudi$lw)
        w <- sweep(w, 2, sqrt(dudi$cw), "*")
        w <- w * w
        con.tra <- apply(w, 2, sum)/sum(w)
        w <- dudi$co * dudi$co * dudi$cw
        w <- sweep(w, 2, dudi$eig[1:nf], "/")
        listing$col.abs <- apply(10000 * w, c(1, 2), app)
        w <- dudi$tab
        w <- sweep(w, 1, sqrt(dudi$lw), "*")
        d2 <- apply(w * w, 2, sum)
        w <- dudi$co * dudi$co
        w <- sweep(w, 1, d2, "/")
        w <- w * sign(dudi$co)
        names(w) <- names(dudi$co)
        w <- cbind.data.frame(w, con.tra)
        listing$col.rel <- apply(10000 * w, c(1, 2), app)
        w <- dudi$co * dudi$co
        w <- sweep(w, 1, d2, "/")
        w <- data.frame(t(apply(w, 1, cumsum)))
        names(w) <- names(dudi$co)
        remain <- 1 - w[, ncol(w)]
        w <- cbind.data.frame(w, remain)
        listing$col.cum <- apply(10000 * w, c(1, 2), app)
    }
    return(listing)
}
"is.dudi" <-
function (x) 
{
    inherits(x, "dudi")
}

"reconst" <-
function (dudi, ...) 
{
    UseMethod("reconst")
}
"reconst.pca" <-
function (dudi, nf = 1, ...) 
{
    if (!inherits(dudi, "dudi")) 
        stop("Object of class 'dudi' expected")
    if (nf > dudi$nf) 
        stop(paste(nf, "factors need >", dudi$nf, "factors available\n"))
    if (!inherits(dudi, "pca")) 
        stop("Object of class 'dudi' expected")
    cent <- dudi$cent
    norm <- dudi$norm
    n <- nrow(dudi$tab)
    p <- ncol(dudi$tab)
    res <- matrix(0, n, p)
    for (i in 1:nf) {
        xli <- dudi$li[, i]
        yc1 <- dudi$c1[, i]
        res <- res + matrix(xli, n, 1) %*% matrix(yc1, 1, p)
    }
    res <- t(apply(res, 1, function(x) x * norm))
    res <- t(apply(res, 1, function(x) x + cent))
    res <- data.frame(res)
    names(res) <- names(dudi$tab)
    row.names(res) <- row.names(dudi$tab)
    return(res)
}
"redo.dudi" <-
function (dudi, newnf = 2) 
{
    if (!inherits(dudi, "dudi")) 
        stop("Object of class 'dudi' expected")
    appel <- as.list(dudi$call)
    if (appel[[1]] == "t.dudi") {
        dudiold <- eval(appel[[2]], sys.frame(0))
        appel <- as.list(dudiold$call)
        appel$nf <- newnf
        appel$scannf <- FALSE
        dudinew <- eval(as.call(appel), sys.frame(0))
        return(t.dudi(dudinew))
    }
    appel$nf <- newnf
    appel$scannf <- FALSE
    eval(as.call(appel), sys.frame(0))
}
"scatter.dudi" <-
function (x, xax = 1, yax = 2, clab.row = 0.5, clab.col = 1, 
    permute = FALSE, posieig = "top", sub = NULL, ...) 
{
    if (!inherits(x, "dudi")) 
        stop("Object of class 'dudi' expected")
    opar <- par(mar = par("mar"))
    on.exit(par(opar))
    coolig <- x$li[, c(xax, yax)]
    coocol <- x$c1[, c(xax, yax)]
    if (permute) {
        coolig <- x$co[, c(xax, yax)]
        coocol <- x$l1[, c(xax, yax)]
    }
    s.label(coolig, clab = clab.row)
    born <- par("usr")
    k1 <- min(coocol[, 1])/born[1]
    k2 <- max(coocol[, 1])/born[2]
    k3 <- min(coocol[, 2])/born[3]
    k4 <- max(coocol[, 2])/born[4]
    k <- c(k1, k2, k3, k4)
    coocol <- 0.9 * coocol/max(k)
    s.arrow(coocol, clab = clab.col, add.p = TRUE, sub = sub, 
        possub = "bottomright")
    add.scatter.eig(x$eig, x$nf, xax, yax, posi = posieig, ratio = 1/4)
}
"scatter.fca" <-
function (x, xax = 1, yax = 2, clab.moda = 1, labels = names(x$tab), 
    sub = NULL, csub = 2, ...) 
{
    opar <- par(mfrow = par("mfrow"))
    on.exit(par(opar))
    if ((xax == yax) || (x$nf == 1)) 
        stop("Unidimensional plot (xax=yax) not yet implemented")
    par(mfrow = n2mfrow(length(x$blo)))
    oritab <- eval(as.list(x$call)[[2]], sys.frame(0))
    indica <- factor(rep(names(x$blo), x$blo))
    for (j in levels(indica)) s.distri(x$l1, oritab[, which(indica == 
        j)], clab = clab.moda, sub = as.character(j), cell = 0, 
        csta = 0.5, csub = csub, label = labels[which(indica == 
            j)])
}

"scatter.pco" <-
function (x, xax = 1, yax = 2, clab.row = 1, posieig = "top", 
    sub = NULL, csub = 2, ...) 
{
    if (!inherits(x, "pco")) 
        stop("Object of class 'pco' expected")
    opar <- par(mar = par("mar"))
    on.exit(par(opar))
    coolig <- x$li[, c(xax, yax)]
    s.label(coolig, clab = clab.row)
    add.scatter.eig(x$eig, x$nf, xax, yax, posi = posieig, ratio = 1/4)
}

"score.mix" <-
function (x, xax = 1, csub = 2, mfrow = NULL, which.var = NULL, ...) 
{
    if (!inherits(x, "mix")) 
        stop("For 'mix' object")
    if (x$nf == 1) 
        xax <- yax <- 1
    lm.pcaiv <- function(x, df, weights, use) {
        if (!inherits(df, "data.frame")) 
            stop("data.frame expected")
        reponse.generic <- x
        begin <- "reponse.generic ~ "
        fmla <- as.formula(paste(begin, paste(names(df), collapse = "+")))
        df <- cbind.data.frame(reponse.generic, df)
        lm0 <- lm(fmla, data = df, weights = weights)
        if (use == 0) 
            return(predict(lm0))
        else if (use == 1) 
            return(residuals(lm0))
        else if (use == -1) 
            return(lm0)
        else stop("Non convenient use")
    }
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    oritab <- eval(as.list(x$call)[[2]], sys.frame(0))
    nvar <- length(x$index)
    if (is.null(which.var)) 
        which.var <- (1:nvar)
    index <- as.character(x$index)
    if (is.null(mfrow)) 
        par(mfrow = n2mfrow(length(which.var)))
    if (prod(par("mfrow")) < length(which.var)) 
        par(ask = TRUE)
    sub <- names(oritab)
    par(mar = c(2.6, 2.6, 1.1, 1.1))
    score <- x$l1[, xax]
    for (i in which.var) {
        type.var <- index[i]
        col.var <- which(x$assign == i)
        if (type.var == "q") {
            if (length(col.var) == 1) {
                y <- x$tab[, col.var]
                plot(score, y, type = "n")
                points(score, y, pch = 20)
                abline(lm(y ~ score), lwd = 2)
            }
            else {
                y <- x$tab[, col.var]
                plot(score, y[, 1], type = "n")
                points(score, y[, 1], pch = 20)
                score.est <- lm.pcaiv(score, y, w = rep(1, nrow(y))/nrow(y), 
                  use = 0)
                ord0 <- order(y[, 1])
                lines(score.est[ord0], y[, 1][ord0], lwd = 2)
            }
        }
        else if (type.var == "f") {
            y <- oritab[, i]
            moy <- unlist(tapply(score, y, mean))
            plot(score, score, type = "n")
            h <- (max(score) - min(score))/40
            abline(h = moy)
            segments(score, moy[y] - h, score, moy[y] + h)
            abline(0, 1)
            scatterutil.eti(moy, moy, label = as.character(levels(y)), 
                clab = 1)
        }
        else if (type.var == "o") {
            y <- x$tab[, col.var]
            plot(score, y[, 1], type = "n")
            points(score, y[, 1], pch = 20)
            score.est <- lm.pcaiv(score, y, w = rep(1, nrow(y))/nrow(y), 
                use = 0)
            ord0 <- order(y[, 1])
            lines(score.est[ord0], y[, 1][ord0])
        }
        scatterutil.sub(sub[i], csub, "topleft")
    }
}
"score.pca" <-
function (x, xax = 1, which.var = NULL, mfrow = NULL, csub = 2, 
    sub = names(x$tab), abline = TRUE, ...) 
{
    if (!inherits(x, "pca")) 
        stop("Object of class 'pca' expected")
    if (x$nf == 1) 
        xax <- 1
    if ((xax < 1) || (xax > x$nf)) 
        stop("non convenient axe number")
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    oritab <- eval(as.list(x$call)[[2]], sys.frame(0))
    nvar <- ncol(oritab)
    if (is.null(which.var)) 
        which.var <- (1:nvar)
    if (is.null(mfrow)) 
        mfrow <- n2mfrow(length(which.var))
    par(mfrow = mfrow)
    if (prod(par("mfrow")) < length(which.var)) 
        par(ask = TRUE)
    par(mar = c(2.6, 2.6, 1.1, 1.1))
    score <- x$l1[, xax]
    for (i in which.var) {
        y <- oritab[, i]
        plot(score, y, type = "n")
        points(score, y, pch = 20)
        if (abline) 
            abline(lm(y ~ score))
        scatterutil.sub(sub[i], csub = csub, "topleft")
    }
}
"supcol" <-
function (x, ...) 
UseMethod("supcol")
"supcol.coa" <-
function (x, Xsup, ...) 
{
    Xsup <- data.frame(Xsup)
    if (!inherits(x, "dudi")) 
        stop("Object of class 'dudi' expected")
    if (!inherits(x, "coa")) 
        stop("Object of class 'coa' expected")
    if (!inherits(Xsup, "data.frame")) 
        stop("Xsup is not a data.frame")
    if (nrow(Xsup) != nrow(x$tab)) 
        stop("non convenient row numbers")
    cwsup <- apply(Xsup, 2, sum)
    cwsup[cwsup == 0] <- 1
    Xsup <- sweep(Xsup, 2, cwsup, "/")
    coosup <- t(as.matrix(Xsup)) %*% as.matrix(x$l1)
    coosup <- data.frame(coosup, row.names = names(Xsup))
    names(coosup) <- names(x$co)
    return(coosup)
}
"supcol.default" <-
function (x, Xsup, ...) 
{
    Xsup <- data.frame(Xsup)
    if (!inherits(x, "dudi")) 
        stop("Object of class 'dudi' expected")
    if (!inherits(Xsup, "data.frame")) 
        stop("Xsup is not a data.frame")
    if (nrow(Xsup) != nrow(x$tab)) 
        stop("non convenient row numbers")
    coosup <- t(as.matrix(Xsup)) %*% (as.matrix(x$l1) * x$lw)
    coosup <- data.frame(coosup, row.names = names(Xsup))
    names(coosup) <- names(x$co)
    return(coosup)
}
"suprow" <-
function (x, ...) 
UseMethod("suprow")
"suprow.coa" <-
function (x, Xsup, ...) 
{
    Xsup <- data.frame(Xsup)
    if (!inherits(x, "dudi")) 
        stop("Object of class 'dudi' expected")
    if (!inherits(x, "coa")) 
        stop("Object of class 'coa' expected")
    if (!inherits(Xsup, "data.frame")) 
        stop("Xsup is not a data.frame")
    if (ncol(Xsup) != ncol(x$tab)) 
        stop("non convenient col numbers")
    lwsup <- apply(Xsup, 1, sum)
    lwsup[lwsup == 0] <- 1
    Xsup <- sweep(Xsup, 1, lwsup, "/")
    coosup <- as.matrix(Xsup) %*% as.matrix(x$c1)
    coosup <- data.frame(coosup, row.names = row.names(Xsup))
    names(coosup) <- names(x$li)
    return(coosup)
}
"suprow.default" <-
function (x, Xsup, ...) 
{
    Xsup <- data.frame(Xsup)
    if (!inherits(x, "dudi")) 
        stop("Object of class 'dudi' expected")
    if (!inherits(Xsup, "data.frame")) 
        stop("Xsup is not a data.frame")
    if (ncol(Xsup) != ncol(x$tab)) 
        stop("non convenient col numbers")
    coosup <- as.matrix(Xsup) %*% t(t(as.matrix(x$c1)) * x$cw)
    coosup <- data.frame(coosup, row.names = row.names(Xsup))
    names(coosup) <- names(x$li)
    return(coosup)
}
"suprow.pca" <-
function (x, Xsup, ...) 
{
    Xsup <- data.frame(Xsup)
    if (!inherits(x, "dudi")) 
        stop("Object of class 'dudi' expected")
    if (!inherits(x, "pca")) 
        stop("Object of class 'pca' expected")
    if (!inherits(Xsup, "data.frame")) 
        stop("Xsup is not a data.frame")
    if (ncol(Xsup) != ncol(x$tab)) 
        stop("non convenient col numbers")
    f1 <- function(w) (w - x$cent)/x$norm
    Xsup <- t(apply(Xsup, 1, f1))
    coosup <- as.matrix(Xsup) %*% as.matrix(x$c1)
    coosup <- data.frame(coosup, row.names = row.names(Xsup))
    names(coosup) <- names(x$li)
    return(coosup)
}
"t.dudi" <-
function (x) 
{
    if (!inherits(x, "dudi")) 
        stop("Object of class 'dudi' expected")
    res <- list()
    res$tab <- data.frame(t(x$tab))
    res$cw <- x$lw
    res$lw <- x$cw
    res$eig <- x$eig
    res$rank <- x$rank
    res$nf <- x$nf
    res$c1 <- x$l1
    res$l1 <- x$c1
    res$co <- x$li
    res$li <- x$co
    res$call <- match.call()
    class(res) <- c("transpo", "dudi")
    return(res)
}

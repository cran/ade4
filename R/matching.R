"between" <-
function (dudi, fac, scannf = TRUE, nf = 2) 
{
    if (!inherits(dudi, "dudi")) 
        stop("Object of class dudi expected")
    if (!is.factor(fac)) 
        stop("factor expected")
    lig <- nrow(dudi$tab)
    col <- ncol(dudi$tab)
    if (length(fac) != lig) 
        stop("Non convenient dimension")
    cla.w <- tapply(dudi$lw, fac, sum)
    mean.w <- function(x, w, fac, cla.w) {
        z <- x * w
        z <- tapply(z, fac, sum)/cla.w
        return(z)
    }
    tabmoy <- apply(dudi$tab, 2, mean.w, w = dudi$lw, fac = fac, 
        cla.w = cla.w)
    tabmoy <- data.frame(tabmoy)
    row.names(tabmoy) <- levels(fac)
    names(tabmoy) <- names(dudi$tab)
    X <- as.dudi(tabmoy, dudi$cw, as.vector(cla.w), scannf = scannf, 
        nf = nf, call = match.call(), type = "bet")
    X$ratio <- sum(X$eig)/sum(dudi$eig)
    U <- as.matrix(X$c1) * unlist(X$cw)
    U <- data.frame(as.matrix(dudi$tab) %*% U)
    row.names(U) <- row.names(dudi$tab)
    names(U) <- names(X$c1)
    X$ls <- U
    U <- as.matrix(X$c1) * unlist(X$cw)
    U <- data.frame(t(as.matrix(dudi$c1)) %*% U)
    row.names(U) <- names(dudi$li)
    names(U) <- names(X$li)
    X$as <- U
    class(X) <- c("between", "dudi")
    return(X)
}
"cca" <-
function (sitspe, sitenv, scannf = TRUE, nf = 2) 
{
    sitenv <- data.frame(sitenv)
    if (!inherits(sitspe, "data.frame")) 
        stop("data.frame expected")
    if (!inherits(sitenv, "data.frame")) 
        stop("data.frame expected")
    coa1 <- dudi.coa(sitspe, scannf = FALSE, nf = 8)
    x <- pcaiv(coa1, sitenv, scannf = scannf, nf = nf)
    class(x) <- c("cca", "pcaiv", "dudi")
    x$call <- match.call()
    return(x)
}
"coinertia" <-
function (dudiX, dudiY, scannf = TRUE, nf = 2) 
{
    normalise.w <- function(X, w) {
        f2 <- function(v) sqrt(sum(v * v * w)/sum(w))
        norm <- apply(X, 2, f2)
        X <- sweep(X, 2, norm, "/")
        return(X)
    }
    if (!inherits(dudiX, "dudi")) 
        stop("Object of class dudi expected")
    lig1 <- nrow(dudiX$tab)
    col1 <- ncol(dudiX$tab)
    if (!inherits(dudiY, "dudi")) 
        stop("Object of class dudi expected")
    lig2 <- nrow(dudiY$tab)
    col2 <- ncol(dudiY$tab)
    if (lig1 != lig2) 
        stop("Non equal row numbers")
    if (any((dudiX$lw - dudiY$lw)^2 > 1e-07)) 
        stop("Non equal row weights")
    tabcoiner <- t(as.matrix(dudiY$tab)) %*% (as.matrix(dudiX$tab) * 
        dudiX$lw)
    tabcoiner <- data.frame(tabcoiner)
    names(tabcoiner) <- names(dudiX$tab)
    row.names(tabcoiner) <- names(dudiY$tab)
    if (nf > dudiX$nf) 
        nf <- dudiX$nf
    if (nf > dudiY$nf) 
        nf <- dudiY$nf
    coi <- as.dudi(tabcoiner, dudiX$cw, dudiY$cw, scannf = scannf, 
        nf = nf, call = match.call(), type = "coinertia")
    U <- as.matrix(coi$c1) * unlist(coi$cw)
    U <- data.frame(as.matrix(dudiX$tab) %*% U)
    row.names(U) <- row.names(dudiX$tab)
    names(U) <- paste("AxcX", (1:coi$nf), sep = "")
    coi$lX <- U
    U <- normalise.w(U, dudiX$lw)
    names(U) <- paste("NorS", (1:coi$nf), sep = "")
    coi$mX <- U
    U <- as.matrix(coi$l1) * unlist(coi$lw)
    U <- data.frame(as.matrix(dudiY$tab) %*% U)
    row.names(U) <- row.names(dudiY$tab)
    names(U) <- paste("AxcY", (1:coi$nf), sep = "")
    coi$lY <- U
    U <- normalise.w(U, dudiY$lw)
    names(U) <- paste("NorS", (1:coi$nf), sep = "")
    coi$mY <- U
    U <- as.matrix(coi$c1) * unlist(coi$cw)
    U <- data.frame(t(as.matrix(dudiX$c1)) %*% U)
    row.names(U) <- paste("Ax", (1:dudiX$nf), sep = "")
    names(U) <- paste("AxcX", (1:coi$nf), sep = "")
    coi$aX <- U
    U <- as.matrix(coi$l1) * unlist(coi$lw)
    U <- data.frame(t(as.matrix(dudiY$c1)) %*% U)
    row.names(U) <- paste("Ax", (1:dudiY$nf), sep = "")
    names(U) <- paste("AxcY", (1:coi$nf), sep = "")
    coi$aY <- U
    RV <- sum(coi$eig)/sqrt(sum(dudiX$eig^2))/sqrt(sum(dudiY$eig^2))
    coi$RV <- RV
    return(coi)
}
"discrimin" <-
function (dudi, fac, scannf = TRUE, nf = 2) 
{
    if (!inherits(dudi, "dudi")) 
        stop("Object of class dudi expected")
    if (!is.factor(fac)) 
        stop("factor expected")
    lig <- nrow(dudi$tab)
    col <- ncol(dudi$tab)
    if (length(fac) != lig) 
        stop("Non convenient dimension")
    rank <- dudi$rank
    dudi <- redo.dudi(dudi, rank)
    deminorm <- as.matrix(dudi$c1) * dudi$cw
    deminorm <- t(t(deminorm)/sqrt(dudi$eig))
    cla.w <- tapply(dudi$lw, fac, sum)
    mean.w <- function(x) {
        z <- x * dudi$lw
        z <- tapply(z, fac, sum)/cla.w
        return(z)
    }
    tabmoy <- apply(dudi$l1, 2, mean.w)
    tabmoy <- data.frame(tabmoy)
    row.names(tabmoy) <- levels(fac)
    cla.w <- cla.w/sum(cla.w)
    X <- as.dudi(tabmoy, rep(1, rank), as.vector(cla.w), scannf = scannf, 
        nf = nf, call = match.call(), type = "dis")
    res <- list()
    res$eig <- X$eig
    res$nf <- X$nf
    res$fa <- deminorm %*% as.matrix(X$c1)
    res$li <- as.matrix(dudi$tab) %*% res$fa
    w <- scalewt(dudi$tab, dudi$lw)
    res$va <- t(as.matrix(w)) %*% (res$li * dudi$lw)
    res$cp <- t(as.matrix(dudi$l1)) %*% (dudi$lw * res$li)
    res$fa <- data.frame(res$fa)
    row.names(res$fa) <- names(dudi$tab)
    names(res$fa) <- paste("DS", 1:X$nf, sep = "")
    res$li <- data.frame(res$li)
    row.names(res$li) <- row.names(dudi$tab)
    names(res$li) <- names(res$fa)
    w <- apply(res$li, 2, mean.w)
    res$gc <- data.frame(w)
    row.names(res$gc) <- as.character(levels(fac))
    names(res$gc) <- names(res$fa)
    res$cp <- data.frame(res$cp)
    row.names(res$cp) <- names(dudi$l1)
    names(res$cp) <- names(res$fa)
    res$call <- match.call()
    class(res) <- "discrimin"
    return(res)
}
"discrimin.coa" <-
function (df, fac, scannf = TRUE, nf = 2) 
{
    if (!is.factor(fac)) 
        stop("factor expected")
    lig <- nrow(df)
    if (length(fac) != lig) 
        stop("Non convenient dimension")
    dudi.coarp <- function(df) {
        if (!is.data.frame(df)) 
            stop("data.frame expected")
        lig <- nrow(df)
        col <- ncol(df)
        if (any(df < 0)) 
            stop("negative entries in table")
        if ((N <- sum(df)) == 0) 
            stop("all frequencies are zero")
        df <- df/N
        row.w <- apply(df, 1, sum)
        col.w <- apply(df, 2, sum)
        if (any(col.w == 0)) 
            stop("null column found in data")
        df <- df/row.w
        df <- sweep(df, 2, col.w)
        X <- as.dudi(df, 1/col.w, row.w, scannf = FALSE, nf = 2, 
            call = match.call(), type = "coarp", full = TRUE)
        X$N <- N
        class(X) <- "dudi"
        return(X)
    }
    dudi <- dudi.coarp(df)
    rank <- dudi$rank
    deminorm <- as.matrix(dudi$c1) * dudi$cw
    deminorm <- t(t(deminorm)/sqrt(dudi$eig))
    cla.w <- as.vector(tapply(dudi$lw, fac, sum))
    mean.w <- function(x) {
        z <- x * dudi$lw
        z <- tapply(z, fac, sum)/cla.w
        return(z)
    }
    tabmoy <- apply(dudi$l1, 2, mean.w)
    tabmoy <- data.frame(tabmoy)
    row.names(tabmoy) <- levels(fac)
    X <- as.dudi(tabmoy, rep(1, rank), cla.w, scannf = scannf, 
        nf = nf, call = match.call(), type = "dis")
    res <- list(eig = X$eig)
    res$nf <- X$nf
    res$fa <- deminorm %*% as.matrix(X$c1)
    res$li <- as.matrix(dudi$tab) %*% res$fa
    w <- scalewt(dudi$tab, dudi$lw)
    res$va <- t(as.matrix(w)) %*% (res$li * dudi$lw)
    res$cp <- t(as.matrix(dudi$l1)) %*% (dudi$lw * res$li)
    res$fa <- data.frame(res$fa)
    row.names(res$fa) <- names(dudi$tab)
    names(res$fa) <- paste("DS", 1:X$nf, sep = "")
    res$li <- data.frame(res$li)
    row.names(res$li) <- row.names(dudi$tab)
    names(res$li) <- names(res$fa)
    w <- apply(res$li, 2, mean.w)
    res$gc <- data.frame(w)
    row.names(res$gc) <- as.character(levels(fac))
    names(res$gc) <- names(res$fa)
    res$va <- data.frame(res$va)
    row.names(res$va) <- names(dudi$tab)
    names(res$va) <- names(res$fa)
    res$cp <- data.frame(res$cp)
    row.names(res$cp) <- names(dudi$l1)
    names(res$cp) <- names(res$fa)
    res$call <- match.call()
    class(res) <- c("coadisc", "discrimin")
    return(res)
}
"niche" <-
function (dudiX, Y, scannf = TRUE, nf = 2) 
{
    if (!inherits(dudiX, "dudi")) 
        stop("Object of class dudi expected")
    lig1 <- nrow(dudiX$tab)
    col1 <- ncol(dudiX$tab)
    if (!is.data.frame(Y)) 
        stop("Y is not a data.frame")
    lig2 <- nrow(Y)
    col2 <- ncol(Y)
    if (lig1 != lig2) 
        stop("Non equal row numbers")
    w1 <- apply(Y, 2, sum)
    if (any(w1 <= 0)) 
        stop(paste("Column sum <=0 in Y"))
    Y <- sweep(Y, 2, w1, "/")
    w1 <- w1/sum(w1)
    tabcoiner <- t(as.matrix(Y)) %*% (as.matrix(dudiX$tab))
    tabcoiner <- data.frame(tabcoiner)
    names(tabcoiner) <- names(dudiX$tab)
    row.names(tabcoiner) <- names(Y)
    if (nf > dudiX$nf) 
        nf <- dudiX$nf
    nic <- as.dudi(tabcoiner, dudiX$cw, w1, scannf = scannf, 
        nf = nf, call = match.call(), type = "niche")
    U <- as.matrix(nic$c1) * unlist(nic$cw)
    U <- data.frame(as.matrix(dudiX$tab) %*% U)
    row.names(U) <- row.names(dudiX$tab)
    names(U) <- names(nic$c1)
    nic$ls <- U
    U <- as.matrix(nic$c1) * unlist(nic$cw)
    U <- data.frame(t(as.matrix(dudiX$c1)) %*% U)
    row.names(U) <- names(dudiX$li)
    names(U) <- names(nic$li)
    nic$as <- U
    return(nic)
}
"pcaiv" <-
function (dudi, df, scannf = TRUE, nf = 2) 
{
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
    if (!inherits(dudi, "dudi")) 
        stop("dudi is not a 'dudi' object")
    df <- data.frame(df)
    if (!inherits(df, "data.frame")) 
        stop("df is not a 'data.frame'")
    if (nrow(df) != length(dudi$lw)) 
        stop("Non convenient dimensions")
    weights <- dudi$lw
    isfactor <- unlist(lapply(as.list(df), is.factor))
    for (i in 1:ncol(df)) {
        if (!isfactor[i]) 
            df[, i] <- scalewt(df[, i], weights)
    }
    tab <- data.frame(apply(dudi$tab, 2, lm.pcaiv, df = df, use = 0, 
        weights = dudi$lw))
    X <- as.dudi(tab, dudi$cw, dudi$lw, scannf = scannf, nf = nf, 
        call = match.call(), type = "pcaiv")
    X$X <- df
    X$Y <- dudi$tab
    U <- as.matrix(X$c1) * unlist(X$cw)
    U <- as.matrix(dudi$tab) %*% U
    U <- data.frame(U)
    row.names(U) <- row.names(dudi$tab)
    names(U) <- names(X$li)
    X$ls <- U
    sumry <- array("", c(X$nf, 7), list(rep("", X$nf), c("iner", 
        "inercum", "inerC", "inercumC", "ratio", "R2", "lambda")))
    sumry[, 1] <- signif(dudi$eig[1:X$nf], dig = 3)
    sumry[, 2] <- signif(cumsum(dudi$eig[1:X$nf]), dig = 3)
    varpro <- apply(U, 2, function(x) sum(x * x * dudi$lw))
    sumry[, 3] <- signif(varpro, dig = 3)
    sumry[, 4] <- signif(cumsum(varpro), dig = 3)
    sumry[, 5] <- signif(cumsum(varpro)/cumsum(dudi$eig[1:X$nf]), 
        dig = 3)
    sumry[, 6] <- signif(X$eig[1:X$nf]/varpro, dig = 3)
    sumry[, 7] <- signif(X$eig[1:X$nf], dig = 3)
    class(sumry) <- "table"
    X$param <- sumry
    U <- as.matrix(X$c1) * unlist(X$cw)
    U <- data.frame(t(as.matrix(dudi$c1)) %*% U)
    row.names(U) <- names(dudi$li)
    names(U) <- names(X$li)
    X$as <- U
    w <- apply(X$ls, 2, function(x) coefficients(lm.pcaiv(x, 
        df, weights, -1)))
    w <- data.frame(w)
    names(w) <- names(X$l1)
    X$fa <- w
    fmla <- as.formula(paste("~ ", paste(names(df), collapse = "+")))
    w <- scalewt(model.matrix(fmla, data = df), weights) * weights
    w <- t(w) %*% as.matrix(X$l1)
    w <- data.frame(w)
    X$cor <- w
    return(X)
}
"pcaivortho" <-
function (dudi, df, scannf = TRUE, nf = 2) 
{
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
    if (!inherits(dudi, "dudi")) 
        stop("dudi is not a 'dudi' object")
    df <- data.frame(df)
    if (!inherits(df, "data.frame")) 
        stop("df is not a 'data.frame'")
    if (nrow(df) != length(dudi$lw)) 
        stop("Non convenient dimensions")
    weights <- dudi$lw
    isfactor <- unlist(lapply(as.list(df), is.factor))
    for (i in 1:ncol(df)) {
        if (!isfactor[i]) 
            df[, i] <- scalewt(df[, i], weights)
    }
    tab <- data.frame(apply(dudi$tab, 2, lm.pcaiv, df = df, use = 1, 
        weights = dudi$lw))
    X <- as.dudi(tab, dudi$cw, dudi$lw, scannf = scannf, nf = nf, 
        call = match.call(), type = "pcaivortho")
    X$X <- df
    X$Y <- dudi$tab
    U <- as.matrix(X$c1) * unlist(X$cw)
    U <- as.matrix(dudi$tab) %*% U
    U <- data.frame(U)
    row.names(U) <- row.names(dudi$tab)
    names(U) <- names(X$li)
    X$ls <- U
    sumry <- array("", c(X$nf, 7), list(rep("", X$nf), c("iner", 
        "inercum", "inerC", "inercumC", "ratio", "R2", "lambda")))
    sumry[, 1] <- signif(dudi$eig[1:X$nf], dig = 3)
    sumry[, 2] <- signif(cumsum(dudi$eig[1:X$nf]), dig = 3)
    varpro <- apply(U, 2, function(x) sum(x * x * dudi$lw))
    sumry[, 3] <- signif(varpro, dig = 3)
    sumry[, 4] <- signif(cumsum(varpro), dig = 3)
    sumry[, 5] <- signif(cumsum(varpro)/cumsum(dudi$eig[1:X$nf]), 
        dig = 3)
    sumry[, 6] <- signif(X$eig[1:X$nf]/varpro, dig = 3)
    sumry[, 7] <- signif(X$eig[1:X$nf], dig = 3)
    class(sumry) <- "table"
    X$param <- sumry
    U <- as.matrix(X$c1) * unlist(X$cw)
    U <- data.frame(t(as.matrix(dudi$c1)) %*% U)
    row.names(U) <- names(dudi$li)
    names(U) <- names(X$li)
    X$as <- U
    return(X)
}
"procuste" <-
function (df1, df2, scale = TRUE, nf = 4, tol = 1e-07) 
{
    df1 <- data.frame(df1)
    df2 <- data.frame(df2)
    if (!is.data.frame(df1)) 
        stop("data.frame expected")
    if (!is.data.frame(df2)) 
        stop("data.frame expected")
    l1 <- nrow(df1)
    if (nrow(df2) != l1) 
        stop("Row numbers are different")
    if (any(row.names(df2) != row.names(df1))) 
        stop("row names are different")
    c1 <- ncol(df1)
    c2 <- ncol(df2)
    X <- scale(df1, scale = FALSE)
    Y <- scale(df2, scale = FALSE)
    var1 <- apply(X, 2, function(x) sum(x^2))
    var2 <- apply(Y, 2, function(x) sum(x^2))
    tra1 <- sum(var1)
    tra2 <- sum(var2)
    if (scale) {
        X <- X/sqrt(tra1)
        Y <- Y/sqrt(tra2)
    }
    PS <- t(X) %*% Y
    svd1 <- svd(PS)
    rank <- sum((svd1$d/svd1$d[1]) > tol)
    if (nf > rank) 
        nf <- rank
    u <- svd1$u[, 1:nf]
    v <- svd1$v[, 1:nf]
    scor1 <- X %*% u
    scor2 <- Y %*% v
    rot1 <- X %*% u %*% t(v)
    rot2 <- Y %*% v %*% t(u)
    res <- list()
    X <- data.frame(X)
    row.names(X) <- row.names(df1)
    names(X) <- names(df1)
    Y <- data.frame(Y)
    row.names(Y) <- row.names(df2)
    names(Y) <- names(df2)
    res$d <- svd1$d
    res$rank <- rank
    res$nfact <- nf
    u <- data.frame(u)
    row.names(u) <- names(df1)
    names(u) <- paste("ax", 1:nf, sep = "")
    v <- data.frame(v)
    row.names(v) <- names(df2)
    names(v) <- paste("ax", 1:nf, sep = "")
    scor1 <- data.frame(scor1)
    row.names(scor1) <- row.names(df1)
    names(scor1) <- paste("ax", 1:nf, sep = "")
    scor2 <- data.frame(scor2)
    row.names(scor2) <- row.names(df1)
    names(scor2) <- paste("ax", 1:nf, sep = "")
    if ((nf == c1) & (nf == c2)) {
        rot1 <- data.frame(rot1)
        row.names(rot1) <- row.names(df1)
        names(rot1) <- names(df2)
        rot2 <- data.frame(rot2)
        row.names(rot2) <- row.names(df1)
        names(rot2) <- names(df1)
        res$rot1 <- rot1
        res$rot2 <- rot2
    }
    res$tab1 <- X
    res$tab2 <- Y
    res$load1 <- u
    res$load2 <- v
    res$scor1 <- scor1
    res$scor2 <- scor2
    res$call <- match.call()
    class(res) <- "procuste"
    return(res)
}
"within" <-
function (dudi, fac, scannf = TRUE, nf = 2) 
{
    if (!inherits(dudi, "dudi")) 
        stop("Object of class dudi expected")
    if (!is.factor(fac)) 
        stop("factor expected")
    lig <- nrow(dudi$tab)
    col <- ncol(dudi$tab)
    if (length(fac) != lig) 
        stop("Non convenient dimension")
    cla.w <- tapply(dudi$lw, fac, sum)
    mean.w <- function(x, w, fac, cla.w) {
        z <- x * w
        z <- tapply(z, fac, sum)/cla.w
        return(z)
    }
    tabmoy <- apply(dudi$tab, 2, mean.w, w = dudi$lw, fac = fac, 
        cla.w = cla.w)
    tabw <- unlist(tapply(dudi$lw, fac, sum))
    tabw <- tabw/sum(tabw)
    tabwit <- dudi$tab - tabmoy[fac, ]
    X <- as.dudi(tabwit, dudi$cw, dudi$lw, scannf = scannf, nf = nf, 
        call = match.call(), type = "wit")
    X$ratio <- sum(X$eig)/sum(dudi$eig)
    U <- as.matrix(X$c1) * unlist(X$cw)
    U <- data.frame(as.matrix(dudi$tab) %*% U)
    row.names(U) <- row.names(dudi$tab)
    names(U) <- names(X$li)
    X$ls <- U
    U <- as.matrix(X$c1) * unlist(X$cw)
    U <- data.frame(t(as.matrix(dudi$c1)) %*% U)
    row.names(U) <- names(dudi$li)
    names(U) <- names(X$li)
    X$as <- U
    X$tabw <- tabw
    X$fac <- fac
    class(X) <- c("within", "dudi")
    return(X)
}
"within.pca" <-
function (df, fac, scaling = c("partial", "total"), scannf = TRUE, 
    nf = 2) 
{
    if (!inherits(df, "data.frame")) 
        stop("Object of class 'data.frame' expected")
    if (!is.factor(fac)) 
        stop("factor expected")
    lig <- nrow(df)
    col <- ncol(df)
    if (length(fac) != lig) 
        stop("Non convenient dimension")
    cla.w <- tapply(rep(1, length(fac)), fac, sum)
    df <- data.frame(scalewt(df))
    mean.w <- function(x) tapply(x, fac, sum)/cla.w
    tabmoy <- apply(df, 2, mean.w)
    tabw <- cla.w
    tabw <- tabw/sum(tabw)
    tabwit <- df
    tabwit <- tabwit - tabmoy[fac, ]
    if (scaling == "total") {
        tabwit <- scalewt(tabwit, center = FALSE, scale = TRUE)
    }
    else if (scaling == "partial") {
        for (j in levels(fac)) {
            w <- tabwit[fac == j, ]
            w <- scalewt(w)
            tabwit[fac == j, ] <- w
        }
    }
    else stop("unknown scaling value")
    tabwit <- data.frame(tabwit)
    for (i in 1:nrow(df)) {
        df[i, ] <- tabwit[i, ] + tabmoy[fac[i], ]
    }
    dudi <- as.dudi(df, row.w = rep(1, nrow(df))/nrow(df), col.w = rep(1, 
        ncol(df)), scannf = FALSE, nf = 4, call = match.call(), 
        type = "tmp")
    X <- as.dudi(tabwit, row.w = rep(1, nrow(df))/nrow(df), col.w = rep(1, 
        ncol(df)), scannf = scannf, nf = nf, call = match.call(), 
        type = "wit")
    X$ratio <- sum(X$eig)/sum(dudi$eig)
    U <- as.matrix(X$c1) * unlist(X$cw)
    U <- data.frame(as.matrix(dudi$tab) %*% U)
    row.names(U) <- row.names(dudi$tab)
    names(U) <- names(X$c1)
    X$ls <- U
    U <- as.matrix(X$c1) * unlist(X$cw)
    U <- data.frame(t(as.matrix(dudi$c1)) %*% U)
    row.names(U) <- names(dudi$li)
    names(U) <- names(X$li)
    X$as <- U
    X$tabw <- tabw
    X$fac <- fac
    class(X) <- c("within", "dudi")
    return(X)
}

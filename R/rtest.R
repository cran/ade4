"as.rtest" <-
function (sim, obs, call = match.call()) 
{
    res <- list(sim = sim, obs = obs)
    res$rep <- length(sim)
    res$pvalue <- (sum(sim >= obs) + 1)/(length(sim) + 1)
    res$call <- call
    class(res) <- "rtest"
    return(res)
}
"plot.rtest" <-
function (x, nclass = 10, coeff = 1, ...) 
{
    if (!inherits(x, "rtest")) 
        stop("Non convenient data")
    obs <- x$obs
    sim <- x$sim
    r0 <- c(sim, obs)
    h0 <- hist(sim, plot = FALSE, nclass = nclass, xlim = xlim0)
    y0 <- max(h0$counts)
    l0 <- max(sim) - min(sim)
    w0 <- l0/(log(length(sim), base = 2) + 1)
    w0 <- w0 * coeff
    xlim0 <- range(r0) + c(-w0, w0)
    hist(sim, plot = TRUE, nclass = nclass, xlim = xlim0, col = grey(0.8), 
        ...)
    lines(c(obs, obs), c(y0/2, 0))
    points(obs, y0/2, pch = 18, cex = 2)
    invisible()
}
"print.rtest" <-
function (x, ...) 
{
    if (!inherits(x, "rtest")) 
        stop("Non convenient data")
    cat("Monte-Carlo test\n")
    cat("Observation:", x$obs, "\n")
    cat("Call: ")
    print(x$call)
    cat("Based on", x$rep, "replicates\n")
    cat("Simulated p-value:", x$pvalue, "\n")
}
"rtest" <-
function (xtest, ...) 
{
    UseMethod("rtest")
}
"rtest.between" <-
function (xtest, nrepet = 99, ...) 
{
    if (!inherits(xtest, "dudi")) 
        stop("Object of class dudi expected")
    if (!inherits(xtest, "between")) 
        stop("Type 'between' expected")
    appel <- as.list(xtest$call)
    dudi1 <- eval(appel$dudi, sys.frame(0))
    fac <- eval(appel$fac, sys.frame(0))
    X <- dudi1$tab
    X.lw <- dudi1$lw
    X.lw <- X.lw/sum(X.lw)
    X.cw <- sqrt(dudi1$cw)
    X <- t(t(X) * X.cw)
    inertot <- sum(dudi1$eig)
    inerinter <- function(perm = TRUE) {
        if (perm) 
            sel <- sample(nrow(X))
        else sel <- 1:nrow(X)
        Y <- X[sel, ]
        Y.lw <- X.lw[sel]
        cla.w <- tapply(Y.lw, fac, sum)
        Y1 <- Y * Y.lw
        Y <- apply(Y * Y.lw, 2, function(x) tapply(x, fac, sum)/cla.w)
        inerb <- sum(apply(Y, 2, function(x) sum(x * x * cla.w)))
        return(inerb/inertot)
    }
    obs <- inerinter(FALSE)
    sim <- unlist(lapply(1:nrepet, inerinter))
    return(as.rtest(sim, obs))
}
"rtest.discrimin" <-
function (xtest, nrepet = 99, ...) 
{
    if (!inherits(xtest, "discrimin")) 
        stop("'discrimin' object expected")
    appel <- as.list(xtest$call)
    dudi <- eval(appel$dudi, sys.frame(0))
    fac <- eval(appel$fac, sys.frame(0))
    lig <- nrow(dudi$tab)
    if (length(fac) != lig) 
        stop("Non convenient dimension")
    rank <- dudi$rank
    dudi <- redo.dudi(dudi, rank)
    dudi.lw <- dudi$lw
    dudi <- dudi$l1
    between.var <- function(x, w, group, group.w) {
        z <- x * w
        z <- tapply(z, group, sum)/group.w
        return(sum(z * z * group.w))
    }
    inertia.ratio <- function(perm = TRUE) {
        if (perm) {
            sigma <- sample(lig)
            Y <- dudi[sigma, ]
            Y.w <- dudi.lw[sigma]
        }
        else {
            Y <- dudi
            Y.w <- dudi.lw
        }
        cla.w <- tapply(Y.w, fac, sum)
        ww <- apply(Y, 2, between.var, w = Y.w, group = fac, 
            group.w = cla.w)
        return(sum(ww)/rank)
    }
    obs <- inertia.ratio(perm = FALSE)
    sim <- unlist(lapply(1:nrepet, inertia.ratio))
    return(as.rtest(sim, obs))
}
"mantel.rtest" <-
function (m1, m2, nrepet = 99) 
{
    if (!inherits(m1, "dist")) 
        stop("Object of class 'dist' expected")
    if (!inherits(m2, "dist")) 
        stop("Object of class 'dist' expected")
    n <- attr(m1, "Size")
    if (n != attr(m2, "Size")) 
        stop("Non convenient dimension")
    permutedist <- function(m) {
        permutevec <- function(v, perm) return(v[perm])
        m <- dist2mat(m)
        n <- ncol(m)
        w0 <- sample(n)
        mperm <- apply(m, 1, permutevec, perm = w0)
        mperm <- t(mperm)
        mperm <- apply(mperm, 2, permutevec, perm = w0)
        return(mat2dist(t(mperm)))
    }
    mantelnoneuclid <- function(m1, m2, nrepet) {
        obs <- cor(unclass(m1), unclass(m2))
        if (nrepet == 0) 
            return(obs)
        perm <- matrix(0, nrow = nrepet, ncol = 1)
        perm <- apply(perm, 1, function(x) cor(unclass(m1), unclass(permutedist(m2))))
        w <- as.rtest(obs = obs, sim = perm, , call = match.call())
        return(w)
    }
    if (is.euclid(m1) & is.euclid(m2)) {
        tab1 <- pcoscaled(m1)
        obs <- cor(dist.quant(tab1, 1), m2)
        if (nrepet == 0) 
            return(obs)
        perm <- rep(0, nrepet)
        perm <- unlist(lapply(perm, function(x) cor(dist(tab1[sample(n), 
            ]), m2)))
        w <- as.rtest(obs = obs, sim = perm, call = match.call())
        return(w)
    }
    w <- mantelnoneuclid(m1, m2, nrepet = nrepet)
    return(w)
}
"procuste.rtest" <-
function (df1, df2, nrepet = 99) 
{
    if (!is.data.frame(df1)) 
        stop("data.frame expected")
    if (!is.data.frame(df2)) 
        stop("data.frame expected")
    l1 <- nrow(df1)
    if (nrow(df2) != l1) 
        stop("Row numbers are different")
    if (any(row.names(df2) != row.names(df1))) 
        stop("row names are different")
    X <- scale(df1, scale = FALSE)
    Y <- scale(df2, scale = FALSE)
    var1 <- apply(X, 2, function(x) sum(x^2))
    var2 <- apply(Y, 2, function(x) sum(x^2))
    tra1 <- sum(var1)
    tra2 <- sum(var2)
    X <- X/sqrt(tra1)
    Y <- Y/sqrt(tra2)
    obs <- sum(svd(t(X) %*% Y)$d)
    if (nrepet == 0) 
        return(obs)
    perm <- matrix(0, nrow = nrepet, ncol = 1)
    perm <- apply(perm, 1, function(x) sum(svd(t(X) %*% Y[sample(l1), 
        ])$d))
    w <- as.rtest(obs = obs, sim = perm, call = match.call())
    return(w)
}

"RV.rtest" <-
function (df1, df2, nrepet = 99) 
{
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
    X <- X/(sum(svd(X)$d^4)^0.25)
    Y <- Y/(sum(svd(Y)$d^4)^0.25)
    obs <- sum(svd(t(X) %*% Y)$d^2)
    if (nrepet == 0) 
        return(obs)
    perm <- matrix(0, nrow = nrepet, ncol = 1)
    perm <- apply(perm, 1, function(x) sum(svd(t(X) %*% Y[sample(l1), 
        ])$d^2))
    w <- as.rtest(obs = obs, sim = perm, call = match.call())
    return(w)
}

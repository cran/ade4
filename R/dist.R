"cailliez" <-
function (distmat, print = FALSE) 
{
    if (is.euclid(distmat)) {
        warning("Euclidean distance found : no correction need")
        return(distmat)
    }
    distmat <- dist2mat(distmat)
    size <- ncol(distmat)
    m1 <- matrix(0, size, size)
    m1 <- rbind(m1, -diag(size))
    m2 <- -bicenter.wt(distmat * distmat)
    m2 <- rbind(m2, 2 * bicenter.wt(distmat))
    m1 <- cbind(m1, m2)
    lambda <- La.eigen(m1, only = TRUE)$values
    c <- max(Re(lambda)[Im(lambda) < 1e-08])
    if (print) 
        cat(paste("Cailliez constant =", round(c, dig = 5), "\n"))
    distmat <- mat2dist(distmat + c)
    attr(distmat, "call") <- match.call()
    attr(distmat, "method") <- "Cailliez"
    return(distmat)
}

"dist.binary" <-
function (df, method = NULL, diag = FALSE, upper = FALSE) 
{
    METHODS <- c("JACCARD", "SOCKAL & MICHENER", "SOCKAL & SNEATH S5", 
        "ROGERS & TANIMOTO", "CZEKANOWSKI", "S9", "OCHIAI", "SOKAL & SNEATH S13", 
        "Phi of PEARSON", "GOWER & LEGENDRE S2")
    if (!inherits(df, "data.frame")) 
        stop("df is not a data.frame")
    if (any(df < 0)) 
        stop("non negative value expected in df")
    d.names <- row.names(df)
    nlig <- nrow(df)
    df <- as.matrix(1 * (df > 0))
    if (is.null(method)) {
        cat("1 = JACCARD index (1901) S3 coefficient of GOWER & LEGENDRE\n")
        cat("s1 = a/(a+b+c) --> d = sqrt(1 - s)\n")
        cat("2 = SOCKAL & MICHENER index (1958) S4 coefficient of GOWER & LEGENDRE \n")
        cat("s2 = (a+d)/(a+b+c+d) --> d = sqrt(1 - s)\n")
        cat("3 = SOCKAL & SNEATH(1963) S5 coefficient of GOWER & LEGENDRE\n")
        cat("s3 = a/(a+2(b+c)) --> d = sqrt(1 - s)\n")
        cat("4 = ROGERS & TANIMOTO (1960) S5 coefficient of GOWER & LEGENDRE\n")
        cat("s4 = (a+d)/(a+2(b+c)+d) --> d = sqrt(1 - s)\n")
        cat("5 = CZEKANOWSKI (1913) or SORENSEN (1948)\n")
        cat("s5 = 2*a/(2*a+b+c) --> d = sqrt(1 - s)\n")
        cat("6 = S9 index of GOWER & LEGENDRE (1986)\n")
        cat("s6 = (a-(b+c)+d)/(a+b+c+d) --> d = sqrt(1 - s)\n")
        cat("7 = OCHIAI (1957) S12 coefficient of GOWER & LEGENDRE\n")
        cat("s7 = a/sqrt((a+b)(a+c)) --> d = sqrt(1 - s)\n")
        cat("8 = SOKAL & SNEATH (1963) S13 coefficient of GOWER & LEGENDRE\n")
        cat("s8 = ad/sqrt((a+b)(a+c)(d+b)(d+c)) --> d = sqrt(1 - s)\n")
        cat("9 = Phi of PEARSON = S14 coefficient of GOWER & LEGENDRE\n")
        cat("s9 = ad-bc)/sqrt((a+b)(a+c)(b+d)(d+c)) --> d = sqrt(1 - s)\n")
        cat("10 = S2 coefficient of GOWER & LEGENDRE\n")
        cat("s10 =  a/(a+b+c+d) --> d = sqrt(1 - s) and unit self-similarity\n")
        cat("Selec an integer (1-10): ")
        method <- as.integer(readLines(n = 1))
    }
    df <- as.matrix(df)
    a <- df %*% t(df)
    b <- df %*% (1 - t(df))
    c <- (1 - df) %*% t(df)
    d <- nlig - a - b - c
    if (method == 1) {
        d <- a/(a + b + c)
    }
    else if (method == 2) {
        d <- (a + d)/(a + b + c + d)
    }
    else if (method == 3) {
        d <- a/(a + 2 * (b + c))
    }
    else if (method == 4) {
        d <- (a + d)/(a + 2 * (b + c) + d)
    }
    else if (method == 5) {
        d <- (a + d)/(a + 2 * (b + c) + d)
    }
    else if (method == 6) {
        d <- 2 * a/(2 * a + b + c)
    }
    else if (method == 7) {
        d <- (a - (b + c) + d)/(a + b + c + d)
    }
    else if (method == 8) {
        d <- a * d/sqrt((a + b) * (a + c) * (d + b) * (d + c))
    }
    else if (method == 9) {
        d <- (a * d - b * c)/sqrt((a + b) * (a + c) * (b + d) * 
            (d + c))
    }
    else if (method == 10) {
        d <- a/(a + b + c + d)
        diag(d) <- 1
    }
    else stop("Non convenient method")
    d <- sqrt(1 - d)
    d <- mat2dist(d)
    attr(d, "Size") <- nlig
    attr(d, "Labels") <- d.names
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- METHODS[method]
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    return(d)
}
"dist.dudi" <-
function (dudi, amongrow = TRUE) 
{
    if (!inherits(dudi, "dudi")) 
        stop("Object of class 'dudi' expected")
    nr <- nrow(dudi$tab)
    nc <- ncol(dudi$tab)
    lw <- dudi$lw
    cw <- dudi$cw
    if (amongrow) {
        x <- t(t(dudi$tab) * sqrt(dudi$cw))
        x <- x %*% t(x)
        y <- diag(x)
        x <- (-2) * x + y
        x <- t(t(x) + y)
        x <- (x + t(x))/2
        diag(x) <- 0
        x <- mat2dist(sqrt(x))
        attr(x, "Labels") <- row.names(dudi$tab)
        attr(x, "method") <- "DUDI"
        return(x)
    }
    else {
        x <- as.matrix(dudi$tab) * sqrt(dudi$lw)
        x <- t(x) %*% x
        y <- diag(x)
        x <- (-2) * x + y
        x <- t(t(x) + y)
        x <- (x + t(x))/2
        diag(x) <- 0
        x <- mat2dist(sqrt(x))
        attr(x, "Labels") <- names(dudi$tab)
        attr(x, "method") <- "DUDI"
        return(x)
    }
}
"dist.neig" <-
function (neig) 
{
    if (!inherits(neig, "neig")) 
        stop("Object of class 'neig' expected")
    res <- neig.util.LtoG(neig)
    n <- nrow(res)
    auxi1 <- res
    auxi2 <- res
    for (itour in 2:n) {
        auxi2 <- auxi2 %*% auxi1
        auxi2[res != 0] <- 0
        diag(auxi2) <- 0
        auxi2 <- (auxi2 > 0) * itour
        if (sum(auxi2) == 0) 
            break
        res <- res + auxi2
    }
    return(mat2dist(res))
}
"dist.prop" <-
function (df, method = NULL, diag = FALSE, upper = FALSE) 
{
    METHODS <- c("d1 Manly", "Overlap index Manly", "Rogers 1972", 
        "Nei 1972", "Edwards 1971")
    if (!inherits(df, "data.frame")) 
        stop("df is not a data.frame")
    if (any(df < 0)) 
        stop("non negative value expected in df")
    dfs <- apply(df, 1, sum)
    if (any(dfs == 0)) 
        stop("row with all zero value")
    df <- df/dfs
    if (is.null(method)) {
        cat("1 = d1 Manly\n")
        cat("d1 = Sum|p(i)-q(i)|/2\n")
        cat("2 = Overlap index Manly\n")
        cat("d2=1-Sum(p(i)q(i))/sqrt(Sum(p(i)^2)/sqrt(Sum(q(i)^2)\n")
        cat("3 = Rogers 1972 (one locus)\n")
        cat("d3=sqrt(0.5*Sum(p(i)-q(i)^2))\n")
        cat("4 = Nei 1972 (one locus)\n")
        cat("d4=-ln(Sum(p(i)q(i)/sqrt(Sum(p(i)^2)/sqrt(Sum(q(i)^2))\n")
        cat("5 = Edwards 1971 (one locus)\n")
        cat("d5= sqrt (1 - (Sum(sqrt(p(i)q(i))))\n")
        cat("Selec an integer (1-5): ")
        method <- as.integer(readLines(n = 1))
    }
    nlig <- nrow(df)
    d <- matrix(0, nlig, nlig)
    d.names <- row.names(df)
    df <- as.matrix(df)
    fun1 <- function(x) {
        p <- df[x[1], ]
        q <- df[x[2], ]
        w <- sum(abs(p - q))/2
        return(w)
    }
    fun2 <- function(x) {
        p <- df[x[1], ]
        q <- df[x[2], ]
        w <- 1 - sum(p * q)/sqrt(sum(p * p))/sqrt(sum(q * q))
        return(w)
    }
    fun3 <- function(x) {
        p <- df[x[1], ]
        q <- df[x[2], ]
        w <- sqrt(0.5 * sum((p - q)^2))
        return(w)
    }
    fun4 <- function(x) {
        p <- df[x[1], ]
        q <- df[x[2], ]
        if (sum(p * q) == 0) 
            stop("sum(p*q)==0 -> non convenient data")
        w <- -log(sum(p * q)/sqrt(sum(p * p))/sqrt(sum(q * q)))
        return(w)
    }
    fun5 <- function(x) {
        p <- df[x[1], ]
        q <- df[x[2], ]
        w <- sqrt(1 - sum(sqrt(p * q)))
        return(w)
    }
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    if (method == 1) 
        d <- unlist(apply(index, 1, fun1))
    else if (method == 2) 
        d <- unlist(apply(index, 1, fun2))
    else if (method == 3) 
        d <- unlist(apply(index, 1, fun3))
    else if (method == 4) 
        d <- unlist(apply(index, 1, fun4))
    else if (method == 5) 
        d <- unlist(apply(index, 1, fun5))
    else stop("Non convenient method")
    attr(d, "Size") <- nlig
    attr(d, "Labels") <- d.names
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- METHODS[method]
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    return(d)
}
"dist.quant" <-
function (df, method = NULL, diag = FALSE, upper = FALSE, tol = 1e-07) 
{
    METHODS <- c("Canonical", "Joreskog", "Mahalanobis")
    df <- data.frame(df)
    if (!inherits(df, "data.frame")) 
        stop("df is not a data.frame")
    if (is.null(method)) {
        cat("1 = Canonical\n")
        cat("d1 = ||x-y|| A=Identity\n")
        cat("2 = Joreskog\n")
        cat("d2=d2 = ||x-y|| A=1/diag(cov)\n")
        cat("3 = Mahalanobis\n")
        cat("d3 = ||x-y|| A=inv(cov)\n")
        cat("Selec an integer (1-3): ")
        method <- as.integer(readLines(n = 1))
    }
    nlig <- nrow(df)
    d <- matrix(0, nlig, nlig)
    d.names <- row.names(df)
    fun1 <- function(x) {
        sqrt(sum((df[x[1], ] - df[x[2], ])^2))
    }
    df <- as.matrix(df)
    index <- cbind(col(d)[col(d) < row(d)], row(d)[col(d) < row(d)])
    if (method == 1) {
        d <- unlist(apply(index, 1, fun1))
    }
    else if (method == 2) {
        dfcov <- cov(df) * (nlig - 1)/nlig
        jor <- diag(dfcov)
        jor[jor == 0] <- 1
        jor <- 1/sqrt(jor)
        df <- t(t(df) * jor)
        d <- unlist(apply(index, 1, fun1))
    }
    else if (method == 3) {
        dfcov <- cov(df) * (nlig - 1)/nlig
        maha <- La.eigen(dfcov, sym = TRUE)
        maha.r <- sum(maha$values > (maha$values[1] * tol))
        maha.e <- 1/sqrt(maha$values[1:maha.r])
        maha.v <- maha$vectors[, 1:maha.r]
        maha.v <- t(t(maha.v) * maha.e)
        df <- df %*% maha.v
        d <- unlist(apply(index, 1, fun1))
    }
    else stop("Non convenient method")
    attr(d, "Size") <- nlig
    attr(d, "Labels") <- d.names
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- METHODS[method]
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    return(d)
}
"dist2mat" <-
function (x) 
{
    size <- attr(x, "Size")
    df <- matrix(0, size, size)
    df[row(df) > col(df)] <- x
    df <- df + t(df)
    labels <- attr(x, "Labels")
    dimnames(df) <- if (is.null(labels)) 
        list(1:size, 1:size)
    else list(labels, labels)
    df
}


"is.euclid" <-
function (distmat, plot = FALSE, print = FALSE, tol = 1e-07) 
{
    if (!inherits(distmat, "dist")) 
        stop("Object of class 'dist' expected")
    distmat <- dist2mat(distmat)
    n <- ncol(distmat)
    delta <- -0.5 * bicenter.wt(distmat * distmat)
    lambda <- La.eigen(delta, symmetric = TRUE, only = TRUE)$values
    w0 <- lambda[n]/lambda[1]
    if (plot) 
        barplot(lambda)
    if (print) 
        print(lambda)
    return((w0 > -tol))
}

"lingoes" <-
function (distmat, print = FALSE) 
{
    if (is.euclid(distmat)) {
        warning("Euclidean distance found : no correction need")
        return(distmat)
    }
    distmat <- dist2mat(distmat)
    n <- ncol(distmat)
    delta <- -0.5 * bicenter.wt(distmat * distmat)
    lambda <- eigen(delta, sym = TRUE)$values
    lder <- lambda[ncol(distmat)]
    distmat <- sqrt(distmat * distmat + 2 * abs(lder))
    if (print) 
        cat("Lingoes constant =", round(abs(lder), dig = 6), 
            "\n")
    distmat <- mat2dist(distmat)
    attr(distmat, "call") <- match.call()
    attr(distmat, "method") <- "Lingoes"
    return(distmat)
}
"mat2dist" <-
function (m, diag = FALSE, upper = FALSE) 
{
    m <- as.matrix(m)
    retval <- m[row(m) > col(m)]
    attributes(retval) <- NULL
    attr(retval, "Labels") <- as.character(1:nrow(m))
    if (!is.null(rownames(m))) 
        attr(retval, "Labels") <- rownames(m)
    else if (!is.null(colnames(m))) 
        attr(retval, "Labels") <- colnames(m)
    attr(retval, "Size") <- nrow(m)
    attr(retval, "Diag") <- diag
    attr(retval, "Upper") <- upper
    attr(retval, "call") <- match.call()
    class(retval) <- "dist"
    retval
}
"pcoscaled" <-
function (distmat, tol = 1e-07) 
{
    if (!inherits(distmat, "dist")) 
        stop("Object of class 'dist' expected")
    if (!is.euclid(distmat)) 
        stop("Euclidean distance expected")
    lab <- attr(distmat, "Labels")
    distmat <- dist2mat(distmat)
    n <- ncol(distmat)
    if (is.null(lab)) 
        lab <- as.character(1:n)
    delta <- -0.5 * bicenter.wt(distmat * distmat)
    eig <- eigen(delta, symmetric = TRUE)
    w0 <- eig$values[n]/eig$values[1]
    if ((w0 < -tol)) 
        stop("Euclidean distance matrix expected")
    ncomp <- sum(eig$values > (eig$values[1] * tol))
    x <- eig$vectors[, 1:ncomp]
    variances <- eig$values[1:ncomp]
    x <- t(apply(x, 1, "*", sqrt(variances)))
    inertot <- sum(variances)
    x <- x/sqrt(inertot)
    x <- data.frame(x)
    names(x) <- paste("C", 1:ncomp, sep = "")
    row.names(x) <- lab
    return(x)
}
"quasieuclid" <-
function (distmat) 
{
    if (is.euclid(distmat)) {
        warning("Euclidean distance found : no correction need")
        return(distmat)
    }
    distmat <- dist2mat(distmat)
    n <- ncol(distmat)
    delta <- -0.5 * bicenter.wt(distmat * distmat)
    eig <- La.eigen(delta, sym = TRUE)
    ncompo <- sum(eig$value > 0)
    tabnew <- eig$vectors[, 1:ncompo] * rep(sqrt(eig$values[1:ncompo]), 
        rep(n, ncompo))
    distmat <- dist.quant(tabnew, 1)
    attr(distmat, "call") <- match.call()
    return(distmat)
}

"print.between" <-
function (x, ...) 
{
    if (!inherits(x, "between")) 
        stop("to be used with 'between' object")
    cat("Between analysis\n")
    cat("call: ")
    print(x$call)
    cat("class: ")
    cat(class(x), "\n")
    cat("\n$nf (axis saved) :", x$nf)
    cat("\n$rank: ", x$rank)
    cat("\n$ratio: ", x$ratio)
    cat("\n\neigen values: ")
    l0 <- length(x$eig)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n\n")
    else cat("\n\n")
    sumry <- array("", c(3, 4), list(1:3, c("vector", "length", 
        "mode", "content")))
    sumry[1, ] <- c("$eig", length(x$eig), mode(x$eig), "eigen values")
    sumry[2, ] <- c("$lw", length(x$lw), mode(x$lw), "group weigths")
    sumry[3, ] <- c("$cw", length(x$cw), mode(x$cw), "col weigths")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(7, 4), list(1:7, c("data.frame", "nrow", 
        "ncol", "content")))
    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "array class-variables")
    sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "class coordinates")
    sumry[3, ] <- c("$l1", nrow(x$l1), ncol(x$l1), "class normed scores")
    sumry[4, ] <- c("$co", nrow(x$co), ncol(x$co), "column coordinates")
    sumry[5, ] <- c("$c1", nrow(x$c1), ncol(x$c1), "column normed scores")
    sumry[6, ] <- c("$ls", nrow(x$ls), ncol(x$ls), "row coordinates")
    sumry[7, ] <- c("$as", nrow(x$as), ncol(x$as), "inertia axis onto between axis")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
}
"print.coinertia" <-
function (x, ...) 
{
    if (!inherits(x, "coinertia")) 
        stop("to be used with 'coinertia' object")
    cat("Coinertia analysis\n")
    cat("call: ")
    print(x$call)
    cat("class: ")
    cat(class(x), "\n")
    cat("\n$rank (rank)     :", x$rank)
    cat("\n$nf (axis saved) :", x$nf)
    cat("\n$RV (RV coeff)   :", x$RV)
    cat("\n\neigen values: ")
    l0 <- length(x$eig)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n\n")
    else cat("\n\n")
    sumry <- array("", c(3, 4), list(1:3, c("vector", "length", 
        "mode", "content")))
    sumry[1, ] <- c("$eig", length(x$eig), mode(x$eig), "eigen values")
    sumry[2, ] <- c("$lw", length(x$lw), mode(x$lw), "row weigths (crossed array)")
    sumry[3, ] <- c("$cw", length(x$cw), mode(x$cw), "col weigths (crossed array)")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(11, 4), list(1:11, c("data.frame", "nrow", 
        "ncol", "content")))
    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "crossed array (CA)")
    sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "Y col = CA row: coordinates")
    sumry[3, ] <- c("$l1", nrow(x$l1), ncol(x$l1), "Y col = CA row: normed scores")
    sumry[4, ] <- c("$co", nrow(x$co), ncol(x$co), "X col = CA column: coordinates")
    sumry[5, ] <- c("$c1", nrow(x$c1), ncol(x$c1), "X col = CA column: normed scores")
    sumry[6, ] <- c("$lX", nrow(x$lX), ncol(x$lX), "row coordinates (X)")
    sumry[7, ] <- c("$mX", nrow(x$mX), ncol(x$mX), "normed row scores (X)")
    sumry[8, ] <- c("$lY", nrow(x$lY), ncol(x$lY), "row coordinates (Y)")
    sumry[9, ] <- c("$mY", nrow(x$mY), ncol(x$mY), "normed row scores (Y)")
    sumry[10, ] <- c("$aX", nrow(x$aX), ncol(x$aX), "axis onto co-inertia axis (X)")
    sumry[11, ] <- c("$aY", nrow(x$aX), ncol(x$aX), "axis onto co-inertia axis (Y)")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
}
"print.discrimin" <-
function (x, ...) 
{
    if (!inherits(x, "discrimin")) 
        stop("to be used with 'discrimin' object")
    cat("Discriminant analysis\n")
    cat("call: ")
    print(x$call)
    cat("class: ")
    cat(class(x), "\n")
    cat("\n$nf (axis saved) :", x$nf)
    cat("\n\neigen values: ")
    l0 <- length(x$eig)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n\n")
    else cat("\n\n")
    sumry <- array("", c(5, 4), list(1:5, c("data.frame", "nrow", 
        "ncol", "content")))
    sumry[1, ] <- c("$fa", nrow(x$fa), ncol(x$fa), "loadings / canonical weights")
    sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "canonical scores")
    sumry[3, ] <- c("$va", nrow(x$va), ncol(x$va), "cos(variables, canonical scores)")
    sumry[4, ] <- c("$cp", nrow(x$cp), ncol(x$cp), "cos(components, canonical scores)")
    sumry[5, ] <- c("$gc", nrow(x$gc), ncol(x$gc), "class scores")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
}
"print.dudi" <-
function (x, ...) 
{
    cat("Duality diagramm\n")
    cat("class: ")
    cat(class(x))
    cat("\n$call: ")
    print(x$call)
    cat("\n$nf:", x$nf, "axis-components saved")
    cat("\n$rank: ")
    cat(x$rank)
    cat("\neigen values: ")
    l0 <- length(x$eig)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n")
    else cat("\n")
    sumry <- array("", c(3, 4), list(1:3, c("vector", "length", 
        "mode", "content")))
    sumry[1, ] <- c("$cw", length(x$cw), mode(x$cw), "column weights")
    sumry[2, ] <- c("$lw", length(x$lw), mode(x$lw), "row weights")
    sumry[3, ] <- c("$eig", length(x$eig), mode(x$eig), "eigen values")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(5, 4), list(1:5, c("data.frame", "nrow", 
        "ncol", "content")))
    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "modified array")
    sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "row coordinates")
    sumry[3, ] <- c("$l1", nrow(x$l1), ncol(x$l1), "row normed scores")
    sumry[4, ] <- c("$co", nrow(x$co), ncol(x$co), "column coordinates")
    sumry[5, ] <- c("$c1", nrow(x$c1), ncol(x$c1), "column normed scores")
    class(sumry) <- "table"
    print(sumry)
    cat("other elements: ")
    if (length(names(x)) > 11) 
        cat(names(x)[12:(length(x))], "\n")
    else cat("NULL\n")
}

"print.neig" <-
function (x, ...) 
{
    deg <- attr(x, "degrees")
    n <- length(deg)
    labels <- names(deg)
    df <- neig.util.LtoG(x)
    for (i in 1:n) {
        w <- c(".", "1")[df[i, 1:i] + 1]
        cat(labels[i], " ", w, "\n", sep = "")
    }
    invisible(df)
}
"print.niche" <-
function (x, ...) 
{
    if (!inherits(x, "niche")) 
        stop("to be used with 'niche' object")
    cat("Niche analysis\n")
    cat("call: ")
    print(x$call)
    cat("class: ")
    cat(class(x), "\n")
    cat("\n$rank (rank)     :", x$rank)
    cat("\n$nf (axis saved) :", x$nf)
    cat("\n$RV (RV coeff)   :", x$RV)
    cat("\n\neigen values: ")
    l0 <- length(x$eig)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n\n")
    else cat("\n\n")
    sumry <- array("", c(3, 4), list(1:3, c("vector", "length", 
        "mode", "content")))
    sumry[1, ] <- c("$eig", length(x$eig), mode(x$eig), "eigen values")
    sumry[2, ] <- c("$lw", length(x$lw), mode(x$lw), "row weigths (crossed array)")
    sumry[3, ] <- c("$cw", length(x$cw), mode(x$cw), "col weigths (crossed array)")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(7, 4), list(1:7, c("data.frame", "nrow", 
        "ncol", "content")))
    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "crossed array (averaging species/sites)")
    sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "species coordinates")
    sumry[3, ] <- c("$l1", nrow(x$l1), ncol(x$l1), "species normed scores")
    sumry[4, ] <- c("$co", nrow(x$co), ncol(x$co), "variables coordinates")
    sumry[5, ] <- c("$c1", nrow(x$c1), ncol(x$c1), "variables normed scores")
    sumry[6, ] <- c("$ls", nrow(x$ls), ncol(x$ls), "sites coordinates")
    sumry[7, ] <- c("$as", nrow(x$as), ncol(x$as), "axis upon niche axis")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
}
"print.pcaiv" <-
function (x, ...) 
{
    if (!inherits(x, "pcaiv")) 
        stop("to be used with 'pcaiv' object")
    cat("Principal Component Analysis with Instrumental Variables\n")
    cat("call: ")
    print(x$call)
    cat("class: ")
    cat(class(x), "\n")
    cat("\n$rank (rank)     :", x$rank)
    cat("\n$nf (axis saved) :", x$nf)
    cat("\n\neigen values: ")
    l0 <- length(x$eig)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n\n")
    else cat("\n\n")
    sumry <- array("", c(3, 4), list(rep("", 3), c("vector", 
        "length", "mode", "content")))
    sumry[1, ] <- c("$eig", length(x$eig), mode(x$eig), "eigen values")
    sumry[2, ] <- c("$lw", length(x$lw), mode(x$lw), "row weigths (from dudi)")
    sumry[3, ] <- c("$cw", length(x$cw), mode(x$cw), "col weigths (from dudi)")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(3, 4), list(rep("", 3), c("data.frame", 
        "nrow", "ncol", "content")))
    sumry[1, ] <- c("$Y", nrow(x$Y), ncol(x$Y), "Dependant variables")
    sumry[2, ] <- c("$X", nrow(x$X), ncol(x$X), "Explanatory variables")
    sumry[3, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "modified array (projected variables)")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(4, 4), list(rep("", 4), c("data.frame", 
        "nrow", "ncol", "content")))
    sumry[1, ] <- c("$c1", nrow(x$c1), ncol(x$c1), "PPA Pseudo Principal Axes")
    sumry[2, ] <- c("$as", nrow(x$as), ncol(x$as), "Principal axis of dudi$tab on PAP")
    sumry[3, ] <- c("$ls", nrow(x$ls), ncol(x$ls), "projection of lines of dudi$tab on PPA")
    sumry[4, ] <- c("$li", nrow(x$li), ncol(x$li), "$ls predicted by X")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(4, 4), list(rep("", 4), c("data.frame", 
        "nrow", "ncol", "content")))
    sumry[1, ] <- c("$fa", nrow(x$fa), ncol(x$fa), "Loadings (CPC as linear combinations of X")
    sumry[2, ] <- c("$l1", nrow(x$l1), ncol(x$l1), "CPC Constraint Principal Components")
    sumry[3, ] <- c("$co", nrow(x$co), ncol(x$co), "inner product CPC - Y")
    sumry[4, ] <- c("$cor", nrow(x$cor), ncol(x$cor), "correlation CPC - X")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    print(x$param)
    cat("\n")
}
"print.procuste" <-
function (x, ...) 
{
    cat("Procustes rotation\n")
    cat("call: ")
    print(x$call)
    cat(paste("class:", class(x)))
    cat(paste("\nrank:", x$rank))
    cat(paste("\naxis number:", x$nfact))
    cat("\nSingular value decomposition: ")
    l0 <- length(x$d)
    cat(signif(x$d, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n")
    else cat("\n")
    cat("tab1   data.frame  ", nrow(x$tab1), "  ", ncol(x$tab1), 
        "   scaled array 1\n")
    cat("tab2   data.frame  ", nrow(x$tab2), "  ", ncol(x$tab2), 
        "   scaled array 2\n")
    cat("scor1  data.frame  ", nrow(x$scor1), " ", ncol(x$scor1), 
        "   row coordinates 1\n")
    cat("scor2  data.frame  ", nrow(x$scor2), " ", ncol(x$scor2), 
        "   row coordinates 2\n")
    cat("load1  data.frame  ", nrow(x$load1), " ", ncol(x$load1), 
        "   loadings 1\n")
    cat("load2  data.frame  ", nrow(x$load2), " ", ncol(x$load2), 
        "   loadings 2\n")
    if (length(names(x)) > 12) {
        cat("other elements: ")
        cat(names(x)[11:(length(x))], "\n")
    }
}

"print.within" <-
function (x, ...) 
{
    if (!inherits(x, "within")) 
        stop("to be used with 'within' object")
    cat("Within analysis\n")
    cat("call: ")
    print(x$call)
    cat("class: ")
    cat(class(x), "\n")
    cat("\n$nf (axis saved) :", x$nf)
    cat("\n$rank: ", x$rank)
    cat("\n$ratio: ", x$ratio)
    cat("\n\neigen values: ")
    l0 <- length(x$eig)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n\n")
    else cat("\n\n")
    sumry <- array("", c(5, 4), list(1:5, c("vector", "length", 
        "mode", "content")))
    sumry[1, ] <- c("$eig", length(x$eig), mode(x$eig), "eigen values")
    sumry[2, ] <- c("$lw", length(x$lw), mode(x$lw), "row weigths")
    sumry[3, ] <- c("$cw", length(x$cw), mode(x$cw), "col weigths")
    sumry[4, ] <- c("$tabw", length(x$tabw), mode(x$tabw), "table weigths")
    sumry[5, ] <- c("$fac", length(x$fac), mode(x$fac), "factor for grouping")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(7, 4), list(1:7, c("data.frame", "nrow", 
        "ncol", "content")))
    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "array class-variables")
    sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "row coordinates")
    sumry[3, ] <- c("$l1", nrow(x$l1), ncol(x$l1), "row normed scores")
    sumry[4, ] <- c("$co", nrow(x$co), ncol(x$co), "column coordinates")
    sumry[5, ] <- c("$c1", nrow(x$c1), ncol(x$c1), "column normed scores")
    sumry[6, ] <- c("$ls", nrow(x$ls), ncol(x$ls), "supplementary row coordinates")
    sumry[7, ] <- c("$as", nrow(x$as), ncol(x$as), "inertia axis onto within axis")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
}

"summary.coinertia" <-
function (object, ...) 
{
    if (!inherits(object, "coinertia")) 
        stop("to be used with 'coinertia' object")
    appel <- as.list(object$call)
    dudiX <- eval(appel$dudiX, sys.frame(0))
    dudiY <- eval(appel$dudiY, sys.frame(0))
    norm.w <- function(X, w) {
        f2 <- function(v) sqrt(sum(v * v * w)/sum(w))
        norm <- apply(X, 2, f2)
        return(norm)
    }
    util <- function(n) {
        x <- "1"
        for (i in 2:n) x[i] <- paste(x[i - 1], i, sep = "")
        return(x)
    }
    eig <- object$eig[1:object$nf]
    covar <- sqrt(eig)
    sdX <- norm.w(object$lX, dudiX$lw)
    sdY <- norm.w(object$lY, dudiX$lw)
    corr <- covar/sdX/sdY
    U <- cbind.data.frame(eig, covar, sdX, sdY, corr)
    row.names(U) <- as.character(1:object$nf)
    cat("\nEigenvalues decomposition:\n")
    print(U)
    cat("\nInertia & coinertia X:\n")
    inertia <- cumsum(sdX^2)
    max <- cumsum(dudiX$eig[1:object$nf])
    ratio <- inertia/max
    U <- cbind.data.frame(inertia, max, ratio)
    row.names(U) <- util(object$nf)
    print(U)
    cat("\nInertia & coinertia Y:\n")
    inertia <- cumsum(sdY^2)
    max <- cumsum(dudiY$eig[1:object$nf])
    ratio <- inertia/max
    U <- cbind.data.frame(inertia, max, ratio)
    row.names(U) <- util(object$nf)
    print(U)
    RV <- sum(object$eig)/sqrt(sum(dudiX$eig^2))/sqrt(sum(dudiY$eig^2))
    cat("\nRV:\n", RV, "\n")
}
"summary.dist" <-
function (object, ...) 
{
    if (!inherits(object, "dist")) 
        stop("For use on the class 'dist'")
    cat("Class: ")
    cat(class(object), "\n")
    cat("Distance matrix by lower triangle : d21, d22, ..., d2n, d32, ...\n")
    cat("Size:", attr(object, "Size"), "\n")
    cat("Labels:", attr(object, "Labels"), "\n")
    cat("call: ")
    print(attr(object, "call"))
    cat("method:", attr(object, "method"), "\n")
    cat("Euclidean matrix (Gower 1966):", is.euclid(object), "\n")
}

"summary.witwit" <-
function (object, ...) 
{
    if (!inherits(object, "witwit")) 
        stop("For 'witwit' object")
    cat("Internal correspondence analysis\n")
    cat("class: ")
    cat(class(object))
    cat("\n$call: ")
    print(object$call)
    cat(object$nf, "axis-components saved")
    cat("\neigen values: ")
    l0 <- length(object$eig)
    cat(signif(object$eig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n")
    else cat("\n")
    cat("\n")
    cat("Eigen value decomposition among row blocks\n")
    nf <- object$nf
    nrb <- nrow(object$lbvar)
    aa <- as.matrix(object$lbvar)
    sumry <- array("", c(nrb + 1, nf + 1), list(c(row.names(object$lbvar), 
        "mean"), c(names(object$lbvar), "weights")))
    sumry[(1:nrb), (1:nf)] <- round(aa, dig = 4)
    sumry[(1:nrb), (nf + 1)] <- round(object$lbw, dig = 4)
    sumry[(nrb + 1), (1:nf)] <- round(object$eig[1:nf], dig = 4)
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(nrb + 1, nf), list(c(row.names(object$lbvar), 
        "sum"), names(object$lbvar)))
    aa <- object$lbvar * object$lbw
    aa <- 1000 * t(t(aa)/object$eig[1:nf])
    sumry[(1:nrb), (1:nf)] <- round(aa, dig = 0)
    sumry[(nrb + 1), (1:nf)] <- rep(1000, nf)
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    cat("Eigen value decomposition among column blocks\n")
    nrb <- nrow(object$cbvar)
    aa <- as.matrix(object$cbvar)
    sumry <- array("", c(nrb + 1, nf + 1), list(c(row.names(object$cbvar), 
        "mean"), c(names(object$cbvar), "weights")))
    sumry[(1:nrb), (1:nf)] <- round(aa, dig = 4)
    sumry[(1:nrb), (nf + 1)] <- round(object$cbw, dig = 4)
    sumry[(nrb + 1), (1:nf)] <- round(object$eig[1:nf], dig = 4)
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(nrb + 1, nf), list(c(row.names(object$cbvar), 
        "sum"), names(object$cbvar)))
    aa <- object$cbvar * object$cbw
    aa <- 1000 * t(t(aa)/object$eig[1:nf])
    sumry[(1:nrb), (1:nf)] <- round(aa, dig = 0)
    sumry[(nrb + 1), (1:nf)] <- rep(1000, nf)
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
}

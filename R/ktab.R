"[.ktab" <-
function (object, selection) 
{
    blocks <- object$blo
    nblo <- length(blocks)
    if (is.logical(selection)) 
        selection <- which(selection)
    if (any(selection > nblo)) 
        stop("Non convenient selection")
    indica <- as.factor(rep(1:nblo, blocks))
    res <- unclass(object)[selection]
    cw <- object$cw
    cw <- split(cw, indica)
    cw <- unlist(cw[selection])
    res$cw <- cw
    res$lw <- object$lw
    nr <- length(res$lw)
    blocks <- unlist(lapply(res, function(x) ncol(x)))
    nblo <- length(blocks)
    res$blo <- blocks
    ktab.util.addfactor(res) <- list(blocks, length(res$lw))
    res$call <- match.call()
    class(res) <- "ktab"
    return(res)
}

"c.ktab" <-
function (...) 
{
    x <- list(...)
    n <- length(x)
    if (any(lapply(x, class) != "ktab")) 
        stop("arguments imply object without 'ktab' class")
    nr <- unlist(lapply(x, function(x) nrow(x[[1]])))
    if (length(unique(nr)) != 1) 
        stop("arguments imply object with non constant row numbers")
    lw <- x[[1]]$lw
    nr <- length(lw)
    noms <- row.names(x[[1]][[1]])
    res <- NULL
    cw <- NULL
    blocks <- NULL
    for (i in 1:n) {
        if (any(x[[i]]$lw != lw)) 
            stop("arguments imply object with non constant row weights")
        if (any(row.names(x[[i]][[1]]) != noms)) 
            stop("arguments imply object with non constant row.names")
        blo.i <- x[[i]]$blo
        nblo.i <- length(blo.i)
        res <- c(res, unclass(x[[i]])[1:nblo.i])
        cw <- c(cw, x[[i]]$cw)
        blocks <- c(blocks, blo.i)
    }
    names(res) <- make.names(names(res), TRUE)
    res$lw <- lw
    res$cw <- cw
    res$blo <- blocks
    nblo <- length(blocks)
    ktab.util.addfactor(res) <- list(blocks, length(lw))
    res$call <- match.call()
    class(res) <- "ktab"
    return(res)
}

"col.names" <-
function (x) 
UseMethod("col.names")
"col.names.ktab" <-
function (x) 
{
    if (!inherits(x, "ktab")) 
        stop("to be used with 'ktab' object")
    ntab <- length(x$blo)
    cha <- unlist(lapply(1:ntab, function(y) attr(x[[y]], "names")))
    return(cha)
}
"col.names<-" <-
function (x, value) 
UseMethod("col.names<-")
"col.names<-.ktab" <-
function (x, value) 
{
    if (!inherits(x, "ktab")) 
        stop("to be used with 'ktab' object")
    ntab <- length(x$blo)
    old <- unlist(lapply(1:ntab, function(y) attr(x[[y]], "names")))
    if (!is.null(old) && length(value) != length(old)) 
        stop("invalid col.names length")
    value <- as.character(value)
    indica <- as.factor(rep(1:ntab, x$blo))
    for (i in 1:ntab) {
        if (any(duplicated(value[indica == i]))) 
            stop("duplicate col.names are not allowed in the same array")
        attr(x[[i]], "names") <- value[indica == i]
    }
    x
}

"foucart" <-
function (X, scannf = TRUE, nf = 2) 
{
    if (!is.list(X)) 
        stop("X is not a list")
    nblo <- length(X)
    if (!all(unlist(lapply(X, is.data.frame)))) 
        stop("a component of X is not a data.frame")
    blocks <- unlist(lapply(X, ncol))
    if (length(unique(blocks)) != 1) 
        stop("non equal col numbers among array")
    blocks <- unlist(lapply(X, nrow))
    if (length(unique(blocks)) != 1) 
        stop("non equal row numbers among array")
    r.n <- row.names(X[[1]])
    for (i in 1:nblo) {
        r.new <- row.names(X[[i]])
        if (any(r.new != r.n)) 
            stop("non equal row.names among array")
    }
    unique.col.names <- names(X[[1]])
    for (i in 1:nblo) {
        c.new <- names(X[[i]])
        if (any(c.new != unique.col.names)) 
            stop
        ("non equal col.names among array")
    }
    for (i in 1:nblo) {
        if (any(X[[i]] < 0)) 
            stop(paste("negative entries in data.frame", i))
        if (sum(X[[i]]) <= 0) 
            stop(paste("Non convenient sum in data.frame", i))
    }
    X <- ktab.list.df(X)
    auxinames <- ktab.util.names(X)
    blocks <- X$blo
    nblo <- length(blocks)
    tnames <- tab.names(X)
    tabm <- X[[1]]/sum(X[[1]])
    for (k in 2:nblo) tabm <- tabm + X[[k]]/sum(X[[k]])
    tabm <- tabm/nblo
    row.names(tabm) <- row.names(X)
    names(tabm) <- unique.col.names
    fouc <- dudi.coa(tabm, scannf = scannf, nf = nf)
    fouc$call <- match.call()
    class(fouc) <- c("foucart", "coa", "dudi")
    cooli <- suprow(fouc, X[[1]])
    for (k in 2:nblo) {
        cooli <- rbind(cooli, suprow(fouc, X[[k]]))
    }
    row.names(cooli) <- auxinames$row
    fouc$Tli <- cooli
    cooco <- supcol(fouc, X[[1]])
    for (k in 2:nblo) {
        cooco <- rbind(cooco, supcol(fouc, X[[k]]))
    }
    row.names(cooco) <- auxinames$col
    fouc$Tco <- cooco
    fouc$TL <- X$TL
    fouc$TC <- X$TC
    fouc$blocks <- blocks
    fouc$tab.names <- tnames
    fouc$call <- match.call()
    return(fouc)
}

"is.ktab" <-
function (x) 
inherits(x, "ktab")
"kplot" <-
function (object, ...) 
{
    UseMethod("kplot")
}

"kplot.foucart" <-
function (object, xax = 1, yax = 2, mfrow = NULL, which.tab = 1:length(object$blo), 
    clab.r = 1, clab.c = 1.25, csub = 2, possub = "bottomright", ...) 
{
    if (!inherits(object, "foucart")) 
        stop("Object of type 'foucart' expected")
    opar <- par(ask = par("ask"), mfrow = par("mfrow"), mar = par("mar"))
    on.exit(par(opar))
    if (is.null(mfrow)) 
        mfrow <- n2mfrow(length(which.tab))
    par(mfrow = mfrow)
    nblo <- length(object$blo)
    if (length(which.tab) > prod(mfrow)) 
        par(ask = TRUE)
    rank.fac <- factor(rep(1:nblo, object$rank))
    nf <- ncol(object$li)
    coolig <- object$Tli[, c(xax, yax)]
    coocol <- object$Tco[, c(xax, yax)]
    names(coocol) <- names(coolig)
    cootot <- rbind.data.frame(coocol, coolig)
    if (clab.r > 0) 
        cpoi <- 0
    else cpoi <- 2
    for (ianal in which.tab) {
        coolig <- object$Tli[object$TL[, 1] == ianal, c(xax, yax)]
        coocol <- object$Tco[object$TC[, 1] == ianal, c(xax, yax)]
        s.label(cootot, clab = 0, cpoi = 0, sub = object$tab.names[ianal], 
            csub = csub, possub = possub)
        s.label(coolig, clab = clab.r, cpoi = cpoi, add.p = TRUE)
        s.label(coocol, clab = clab.c, add.p = TRUE)
    }
}

"kplot.mcoa" <-
function (object, xax = 1, yax = 2, which.tab = 1:nrow(object$cov2), 
    mfrow = NULL, option = c("points", "axis", "columns"), clab = 1, 
    cpoint = 2, csub = 2, possub = "bottomright", ...) 
{
    if (!inherits(object, "mcoa")) 
        stop("Object of type 'mcoa' expected")
    opar <- par(ask = par("ask"), mfrow = par("mfrow"), mar = par("mar"))
    on.exit(par(opar))
    if (option == "points") {
        if (is.null(mfrow)) 
            mfrow <- n2mfrow(length(which.tab) + 1)
        par(mfrow = mfrow)
        if (length(which.tab) > prod(mfrow) - 1) 
            par(ask = TRUE)
        par(mar = c(0.1, 0.1, 0.1, 0.1))
        coo1 <- object$SynVar[, c(xax, yax)]
        cootot <- object$Tl1[, c(xax, yax)]
        names(cootot) <- names(coo1)
        coofull <- coo1
        for (i in which.tab) coofull <- rbind.data.frame(coofull, 
            cootot[object$TL[, 1] == i, ])
        s.label(coo1, clab = clab, sub = "Reference", possub = "bottomright", 
            csub = csub)
        for (ianal in which.tab) {
            scatterutil.base(coofull, 1, 2, xlim = NULL, ylim = NULL, 
                grid = TRUE, addaxes = TRUE, cgrid = 1, include.origin = TRUE, 
                origin = c(0, 0), sub = row.names(object$cov2)[ianal], 
                csub = csub, possub = possub, pixmap = NULL, 
                contour = NULL, area = NULL, add.plot = FALSE)
            coo2 <- cootot[object$TL[, 1] == ianal, 1:2]
            s.match(coo1, coo2, clab = 0, add.p = TRUE)
            s.label(coo1, clab = 0, cpoi = cpoint, add.p = TRUE)
        }
        return(invisible())
    }
    if (is.null(mfrow)) 
        mfrow <- n2mfrow(length(which.tab))
    par(mfrow = mfrow)
    if (option == "axis") {
        if (length(which.tab) > prod(mfrow)) 
            par(ask = TRUE)
        for (ianal in which.tab) {
            coo2 <- object$Tax[object$T4[, 1] == ianal, c(xax, yax)]
            row.names(coo2) <- as.character(1:4)
            s.corcircle(coo2, clab = clab, sub = row.names(object$cov2)[ianal], 
                csub = csub, possub = possub)
        }
        return(invisible())
    }
    if (option == "columns") {
        if (length(which.tab) > prod(mfrow)) 
            par(ask = TRUE)
        for (ianal in which.tab) {
            coo2 <- object$Tco[object$TC[, 1] == ianal, c(xax, yax)]
            s.arrow(coo2, clab = clab, sub = row.names(object$cov2)[ianal], 
                csub = csub, possub = possub)
        }
        return(invisible())
    }
}
"kplot.mfa" <-
function (object, xax = 1, yax = 2, mfrow = NULL, which.tab = 1:length(object$blo), 
    row.names = FALSE, col.names = TRUE, traject = FALSE, permute.row.col = FALSE, 
    clab = 1, csub = 2, possub = "bottomright", ...) 
{
    if (!inherits(object, "mfa")) 
        stop("Object of type 'mfa' expected")
    opar <- par(ask = par("ask"), mfrow = par("mfrow"), mar = par("mar"))
    on.exit(par(opar))
    if (is.null(mfrow)) 
        mfrow <- n2mfrow(length(which.tab))
    par(mfrow = mfrow)
    nbloc <- length(object$blo)
    if (length(which.tab) > prod(mfrow)) 
        par(ask = TRUE)
    rank.fac <- factor(rep(1:nbloc, object$rank))
    nf <- ncol(object$li)
    for (ianal in which.tab) {
        coolig <- object$lisup[object$TL[, 1] == ianal, c(xax, yax)]
        coocol <- object$co[object$TC[, 1] == ianal, c(xax, yax)]
        if (permute.row.col) {
            auxi <- coolig
            coolig <- coocol
            coocol <- auxi
        }
        cl <- clab * row.names
        if (cl > 0) 
            cpoi <- 0
        else cpoi <- 2
        s.label(coolig, clab = cl, cpoi = cpoi)
        if (traject) 
            s.traject(coolig, clab = 0, add.p = TRUE)
        born <- par("usr")
        k1 <- min(coocol[, 1])/born[1]
        k2 <- max(coocol[, 1])/born[2]
        k3 <- min(coocol[, 2])/born[3]
        k4 <- max(coocol[, 2])/born[4]
        k <- c(k1, k2, k3, k4)
        coocol <- 0.7 * coocol/max(k)
        s.arrow(coocol, clab = clab * col.names, add.p = TRUE, 
            sub = object$tab.names[ianal], possub = possub, csub = csub)
    }
}
"kplot.pta" <-
function (object, xax = 1, yax = 2, which.tab = 1:nrow(object$RV), 
    mfrow = NULL, which.graph = 1:4, clab = 1, cpoint = 2, csub = 2, 
    possub = "bottomright", ask = par("ask"), ...) 
{
    if (!inherits(object, "pta")) 
        stop("Object of type 'pta' expected")
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    show <- rep(FALSE, 4)
    if (!is.numeric(which.graph) || any(which.graph < 1) || any(which.graph > 
        4)) 
        stop("`which' must be in 1:4")
    show[which.graph] <- TRUE
    if (is.null(mfrow)) {
        mfcol <- c(length(which.tab), length(which.graph))
        par(mfcol = mfcol)
    }
    else par(mfrow = mfrow)
    par(ask = ask)
    if (show[1]) {
        for (ianal in which.tab) {
            coo2 <- object$Tax[object$T4[, 1] == ianal, c(xax, yax)]
            row.names(coo2) <- as.character(1:4)
            s.corcircle(coo2, clab = clab, cgrid = 0, sub = row.names(object$RV)[ianal], 
                csub = csub, possub = possub)
        }
    }
    if (show[2]) {
        par(mar = c(0.1, 0.1, 0.1, 0.1))
        coo1 <- object$li[, c(xax, yax)]
        cootot <- object$Tli[, c(xax, yax)]
        names(cootot) <- names(coo1)
        coofull <- coo1
        for (i in which.tab) coofull <- rbind.data.frame(coofull, 
            cootot[object$TL[, 1] == i, ])
        for (ianal in which.tab) {
            scatterutil.base(coofull, 1, 2, xlim = NULL, ylim = NULL, 
                grid = TRUE, addaxes = TRUE, cgrid = 1, include.origin = TRUE, 
                origin = c(0, 0), sub = row.names(object$RV)[ianal], 
                csub = csub, possub = possub, pixmap = NULL, 
                contour = NULL, area = NULL, add.plot = FALSE)
            coo2 <- cootot[object$TL[, 1] == ianal, 1:2]
            s.label(coo2, add.p = TRUE, clab = clab, label = row.names(object$Tli)[object$TL[, 
                1] == ianal])
        }
    }
    if (show[3]) {
        par(mar = c(0.1, 0.1, 0.1, 0.1))
        coo1 <- object$co[, c(xax, yax)]
        cootot <- object$Tco[, c(xax, yax)]
        names(cootot) <- names(coo1)
        coofull <- coo1
        for (i in which.tab) coofull <- rbind.data.frame(coofull, 
            cootot[object$TC[, 1] == i, ])
        for (ianal in which.tab) {
            scatterutil.base(coofull, 1, 2, xlim = NULL, ylim = NULL, 
                grid = TRUE, addaxes = TRUE, cgrid = 1, include.origin = TRUE, 
                origin = c(0, 0), sub = row.names(object$RV)[ianal], 
                csub = csub, possub = possub, pixmap = NULL, 
                contour = NULL, area = NULL, add.plot = FALSE)
            coo2 <- object$Tco[object$TC[, 1] == ianal, c(xax, yax)]
            s.arrow(coo2, add.p = TRUE, clab = clab, sub = row.names(object$RV)[ianal], 
                csub = csub, possub = possub)
        }
    }
    if (show[4]) {
        for (ianal in which.tab) {
            coo2 <- object$Tcomp[object$T4[, 1] == ianal, c(xax, yax)]
            row.names(coo2) <- as.character(1:4)
            s.corcircle(coo2, clab = clab, cgrid = 0, sub = row.names(object$RV)[ianal], 
                csub = csub, possub = possub)
        }
    }
}
"kplot.sepan" <-
function (object, xax = 1, yax = 2, which.tab = 1:length(object$blo), 
    mfrow = NULL, permute.row.col = FALSE, clab.row = 1, clab.col = 1.25, 
    traject.row = FALSE, csub = 2, possub = "bottomright", show.eigen.value = TRUE, ...) 
{
    if (!inherits(object, "sepan")) 
        stop("Object of type 'sepan' expected")
    opar <- par(ask = par("ask"), mfrow = par("mfrow"), mar = par("mar"))
    on.exit(par(opar))
    nbloc <- length(object$blo)
    if (is.null(mfrow)) 
        mfrow <- n2mfrow(length(which.tab))
    par(mfrow = mfrow)
    if (length(which.tab) > prod(mfrow)) 
        par(ask = TRUE)
    rank.fac <- factor(rep(1:nbloc, object$rank))
    nf <- ncol(object$Li)
    neig <- max(object$rank)
    appel <- as.list(object$call)
    X <- eval(appel$X, sys.frame(0))
    names.li <- row.names(X[[1]])
    for (ianal in which.tab) {
        coolig <- object$Li[object$TL[, 1] == ianal, c(xax, yax)]
        row.names(coolig) <- names.li
        coocol <- object$Co[object$TC[, 1] == ianal, c(xax, yax)]
        row.names(coocol) <- names(X[[ianal]])
        if (permute.row.col) {
            auxi <- coolig
            coolig <- coocol
            coocol <- auxi
        }
        if (clab.row > 0) 
            cpoi <- 0
        else cpoi <- 2
        if (!traject.row) 
            s.label(coolig, clab = clab.row, cpoi = cpoi)
        else s.traject(coolig, clab = 0, cpoi = 2)
        born <- par("usr")
        k1 <- min(coocol[, 1])/born[1]
        k2 <- max(coocol[, 1])/born[2]
        k3 <- min(coocol[, 2])/born[3]
        k4 <- max(coocol[, 2])/born[4]
        k <- c(k1, k2, k3, k4)
        coocol <- 0.7 * coocol/max(k)
        s.arrow(coocol, clab = clab.col, add.p = TRUE, sub = object$tab.names[ianal], 
            csub = csub, possub = possub)
        w <- object$Eig[rank.fac == ianal]
        if (length(w) < neig) 
            w <- c(w, rep(0, neig - length(w)))
        if (show.eigen.value) 
            add.scatter.eig(w, nf, xax, yax, posi = c("bottom", 
                "top"), ratio = 1/4)
    }
}
"kplot.sepan.coa" <-
function (object, xax = 1, yax = 2, which.tab = 1:length(object$blo), 
    mfrow = NULL, permute.row.col = FALSE, clab.row = 1, clab.col = 1.25, 
    csub = 2, possub = "bottomright", show.eigen.value = TRUE, 
    poseig = c("bottom", "top"), ...) 
{
    if (!inherits(object, "sepan")) 
        stop("Object of type 'sepan' expected")
    opar <- par(ask = par("ask"), mfrow = par("mfrow"), mar = par("mar"))
    on.exit(par(opar))
    nbloc <- length(object$blo)
    if (is.null(mfrow)) 
        mfrow <- n2mfrow(length(which.tab))
    par(mfrow = mfrow)
    if (length(which.tab) > prod(mfrow)) 
        par(ask = TRUE)
    rank.fac <- factor(rep(1:nbloc, object$rank))
    nf <- ncol(object$Li)
    neig <- max(object$rank)
    appel <- as.list(object$call)
    X <- eval(appel$X, sys.frame(0))
    names.li <- row.names(X[[1]])
    for (ianal in which.tab) {
        coocol <- object$C1[object$TC[, 1] == ianal, c(xax, yax)]
        row.names(coocol) <- names(X[[ianal]])
        coolig <- object$Li[object$TL[, 1] == ianal, c(xax, yax)]
        row.names(coolig) <- names.li
        if (permute.row.col) {
            auxi <- coolig
            coolig <- coocol
            coocol <- auxi
        }
        if (clab.col > 0) 
            cpoi <- 0
        else cpoi <- 3
        s.label(coocol, clab = 0, cpoi = 0, sub = object$tab.names[ianal], 
            csub = csub, possub = possub)
        s.label(coocol, clab = clab.col, cpoi = cpoi, add.p = TRUE)
        s.label(coolig, clab = clab.row, add.p = TRUE)
        if (permute.row.col) {
            auxi <- coolig
            coolig <- coocol
            coocol <- auxi
        }
        w <- object$Eig[rank.fac == ianal]
        if (length(w) < neig) 
            w <- c(w, rep(0, neig - length(w)))
        if (show.eigen.value) 
            add.scatter.eig(w, nf, xax, yax, posi = poseig, ratio = 1/4)
    }
}
"kplot.statis" <-
function (object, xax = 1, yax = 2, mfrow = NULL, which.tab = 1:length(object$tab.names), 
    clab = 1.5, cpoi = 2, traject = FALSE, arrow = TRUE, class = NULL, 
    unique.scale = FALSE, csub = 2, possub = "bottomright", ...) 
{
    if (!inherits(object, "statis")) 
        stop("Object of type 'statis' expected")
    opar <- par(ask = par("ask"), mfrow = par("mfrow"), mar = par("mar"))
    on.exit(par(opar))
    if (is.null(mfrow)) 
        mfrow <- n2mfrow(length(which.tab))
    par(mfrow = mfrow)
    if (length(which.tab) > prod(mfrow)) 
        par(ask = TRUE)
    nbloc <- length(object$RV.tabw)
    nf <- ncol(object$C.Co)
    if (xax > nf) 
        stop("Non convenient xax")
    if (yax > nf) 
        stop("Non convenient yax")
    cootot <- object$C.Co[, c(xax, yax)]
    label <- TRUE
    if (!is.null(class)) {
        class <- factor(class)
        if (length(class) != length(object$TC[, 1])) 
            class <- NULL
        else label <- FALSE
    }
    for (ianal in which.tab) {
        coocol <- cootot[object$TC[, 1] == ianal, ]
        if (unique.scale) 
            s.label(cootot, clab = 0, cpoi = 0, sub = object$tab.names[ianal], 
                possub = possub, csub = csub)
        else s.label(coocol, clab = 0, cpoi = 0, sub = object$tab.names[ianal], 
            possub = possub, csub = csub)
        if (arrow) {
            s.arrow(coocol, clab = clab, add.p = TRUE)
            label <- FALSE
        }
        if (label) 
            s.label(coocol, clab = clab, cpoi = cpoi, add.p = TRUE)
        if (traject) 
            s.traject(coocol, clab = 0, add.p = TRUE)
        if (!is.null(class)) {
            f1 <- as.factor(class[object$TC[, 1] == ianal])
            s.class(coocol, f1, clab = clab, cpoi = 2, 
                pch = 20, axese = FALSE, cell = 0, add.plot = TRUE)
        }
    }
}
"ktab.data.frame" <-
function (df, blocks, rownames = NULL, colnames = NULL, tabnames = NULL, 
    w.row = rep(1, nrow(df))/nrow(df), w.col = rep(1, ncol(df))) 
{
    if (!inherits(df, "data.frame")) 
        stop("object 'data.frame' expected")
    nblo <- length(blocks)
    if (sum(blocks) != ncol(df)) 
        stop("Non convenient 'blocks' parameter")
    if (is.null(rownames)) 
        rownames <- row.names(df)
    else if (length(rownames) != length(row.names(df))) 
        stop("Non convenient rownames length")
    if (is.null(colnames)) 
        colnames <- names(df)
    else if (length(colnames) != length(names(df))) 
        stop("Non convenient colnames length")
    if (is.null(names(blocks))) 
        tn <- paste("Ana", 1:nblo, sep = "")
    else tn <- names(blocks)
    if (is.null(tabnames)) 
        tabnames <- tn
    else if (length(tabnames) != length(tn)) 
        stop("Non convenient tabnames length")
    for (x in c("lw", "cw", "blo", "TL", "TC", "T4")) tabnames[tabnames == 
        x] <- paste(x, "*", sep = "")
    indica <- as.factor(rep(1:nblo, blocks))
    res <- list()
    for (i in 1:nblo) {
        res[[i]] <- df[, indica == i]
    }
    names(blocks) <- tabnames
    res$lw <- w.row
    res$cw <- w.col
    res$blo <- blocks
    ktab.util.addfactor(res) <- list(blocks, length(res$lw))
    res$call <- match.call()
    class(res) <- "ktab"
    row.names(res) <- rownames
    col.names(res) <- colnames
    tab.names(res) <- tabnames
    return(res)
}
"ktab.list.df" <-
function (obj, rownames = NULL, colnames = NULL, tabnames = NULL, 
    w.row = rep(1, nrow(obj[[1]])), w.col = lapply(obj, function(x) rep(1/ncol(x), 
        ncol(x)))) 
{
    obj <- as.list(obj)
    if (any(unlist(lapply(obj, function(x) !inherits(x, "data.frame"))))) 
        stop("list of 'data.frame' object expected")
    nblo <- length(obj)
    res <- list()
    nlig <- nrow(obj[[1]])
    blocks <- unlist(lapply(obj, function(x) ncol(x)))
    cn <- unlist(lapply(obj, names))
    if (is.null(rownames)) 
        rownames <- row.names(obj[[1]])
    else if (length(rownames) != length(row.names(obj[[1]]))) 
        stop("Non convenient rownames length")
    if (is.null(colnames)) 
        colnames <- cn
    else if (length(colnames) != length(cn)) 
        stop("Non convenient colnames length")
    if (is.null(names(obj))) 
        tn <- paste("Ana", 1:nblo, sep = "")
    else tn <- names(obj)
    if (is.null(tabnames)) 
        tabnames <- tn
    else if (length(tabnames) != length(tn)) 
        stop("Non convenient tabnames length")
    if (nlig != length(w.row)) 
        stop("Non convenient length for w.row")
    n1 <- unlist(lapply(w.col, length))
    n2 <- unlist(lapply(obj, ncol))
    if (any(n1 != n2)) 
        stop("Non convenient length in  w.col")
    for (i in 1:nblo) {
        res[[i]] <- obj[[i]]
    }
    lw <- w.row
    cw <- unlist(w.col)
    names(cw) <- NULL
    names(blocks) <- tabnames
    res$blo <- blocks
    res$lw <- lw
    res$cw <- cw
    ktab.util.addfactor(res) <- list(blocks, length(lw))
    res$call <- match.call()
    class(res) <- "ktab"
    row.names(res) <- rownames
    col.names(res) <- colnames
    tab.names(res) <- tabnames
    return(res)
}
"ktab.list.dudi" <-
function (obj, rownames = NULL, colnames = NULL, tabnames = NULL) 
{
    obj <- as.list(obj)
    if (any(unlist(lapply(obj, function(x) !inherits(x, "dudi"))))) 
        stop("list of object 'dudi' expected")
    nblo <- length(obj)
    res <- list()
    lw <- obj[[1]]$lw
    cw <- NULL
    nlig <- nrow(obj[[1]]$tab)
    blocks <- unlist(lapply(obj, function(x) ncol(x$tab)))
    for (i in 1:nblo) {
        if (any(obj[[i]]$lw != lw)) 
            stop("Non equal row weights among arrays")
        res[[i]] <- obj[[i]]$tab
        cw <- c(cw, obj[[i]]$cw)
    }
    cn <- unlist(lapply(obj, function(x) names(x$tab)))
    if (is.null(rownames)) 
        rownames <- row.names(obj[[1]]$tab)
    else if (length(rownames) != length(row.names(obj[[1]]$tab))) 
        stop("Non convenient rownames length")
    if (is.null(colnames)) 
        colnames <- cn
    else if (length(colnames) != length(cn)) 
        stop("Non convenient colnames length")
    if (is.null(names(obj))) 
        tn <- paste("Ana", 1:nblo, sep = "")
    else tn <- names(obj)
    if (is.null(tabnames)) 
        tabnames <- tn
    else if (length(tabnames) != length(tn)) 
        stop("Non convenient tabnames length")
    names(blocks) <- tabnames
    res$blo <- blocks
    res$lw <- lw
    res$cw <- cw
    ktab.util.addfactor(res) <- list(blocks, length(lw))
    res$call <- match.call()
    class(res) <- "ktab"
    row.names(res) <- rownames
    col.names(res) <- colnames
    tab.names(res) <- tabnames
    return(res)
}
"ktab.util.addfactor<-" <-
function (x, value) 
{
    blocks <- value[[1]]
    nlig <- value[[2]]
    nblo <- length(blocks)
    w <- cbind.data.frame(gl(nblo, nlig), factor(rep(1:nlig, 
        nblo)))
    names(w) <- c("T", "L")
    x$TL <- w
    w <- NULL
    for (i in 1:nblo) w <- c(w, 1:blocks[i])
    w <- cbind.data.frame(factor(rep(1:nblo, blocks)), factor(w))
    names(w) <- c("T", "C")
    x$TC <- w
    w <- cbind.data.frame(gl(nblo, 4), factor(rep(1:4, nblo)))
    names(w) <- c("T", "4")
    x$T4 <- w
    x
}
"ktab.util.names" <-
function (x) 
{
    w <- row.names(x)
    w1 <- paste(w, as.character(x$TL[, 1]), sep = ".")
    w <- col.names(x)
    if (any(dups <- duplicated(w))) 
        w <- paste(w, as.character(x$TC[, 1]), sep = ".")
    w2 <- w
    w <- tab.names(x)
    l0 <- length(w)
    w3 <- paste(rep(w, rep(4, l0)), as.character(1:4), sep = ".")
    return(list(row = w1, col = w2, tab = w3))
}
"ktab.within" <-
function (dudiwit, rownames = NULL, colnames = NULL, tabnames = NULL) 
{
    if (!inherits(dudiwit, "within")) 
        stop("Result from within expected for dudiwit")
    fac <- dudiwit$fac
    res <- list()
    nblo <- nlevels(fac)
    res <- list()
    blocks <- rep(0, nblo)
    if (is.null(rownames)) 
        rownames <- names(dudiwit$tab)
    else if (length(rownames) != length(names(dudiwit$tab))) 
        stop("Non convenient rownames length")
    if (is.null(colnames)) 
        colnames <- row.names(dudiwit$tab)
    else if (length(colnames) != length(row.names(dudiwit$tab))) 
        stop("Non convenient colnames length")
    if (is.null(tabnames)) 
        tabnames <- as.character(unique(fac))
    else if (length(tabnames) != length(as.character(unique(fac)))) 
        stop("Non convenient tabnames length")
    nlig <- ncol(dudiwit$tab)
    cw <- NULL
    for (i in 1:nblo) {
        k <- unique(fac)[i]
        w1 <- dudiwit$lw[fac == k]
        w1 <- w1/sum(w1)
        cw <- c(cw, w1)
        res[[i]] <- data.frame(t(dudiwit$tab[fac == k, ]))
        blocks[i] <- ncol(res[[i]])
    }
    names(blocks) <- tabnames
    res$lw <- dudiwit$cw
    res$cw <- cw
    res$blo <- blocks
    ktab.util.addfactor(res) <- list(blocks, length(res$lw))
    res$call <- match.call()
    res$tabw <- dudiwit$tabw
    class(res) <- "ktab"
    row.names(res) <- rownames
    col.names(res) <- colnames
    tab.names(res) <- tabnames
    return(res)
}

"mcoa" <-
function (X, option = c("inertia", "lambda1", "uniform", "internal"), 
    scannf = TRUE, nf = 3, tol = 1e-07) 
{
    if (!inherits(X, "ktab")) 
        stop("object 'ktab' expected")
    if (option == "internal") {
        if (is.null(X$tabw)) {
            warning("Internal weights not found: uniform weigths are used")
            option <- "uniform"
        }
    }
    lw <- X$lw
    nlig <- length(lw)
    cw <- X$cw
    ncol <- length(cw)
    blo <- X$blo
    nbloc <- length(X$blo)
    indicablo <- X$TC[, 1]
    Xsepan <- sepan(X, nf = 4)
    rank.fac <- factor(rep(1:nbloc, Xsepan$rank))
    tabw <- NULL
    auxinames <- ktab.util.names(X)
    if (option == "lambda1") {
        for (i in 1:nbloc) tabw <- c(tabw, 1/Xsepan$Eig[rank.fac == 
            i][1])
    }
    else if (option == "inertia") {
        for (i in 1:nbloc) tabw <- c(tabw, 1/sum(Xsepan$Eig[rank.fac == 
            i]))
    }
    else if (option == "uniform") {
        tabw <- rep(1, nbloc)
    }
    else if (option == "internal") 
        tabw <- X$tabw
    else stop("Unknown option")
    for (i in 1:nbloc) X[[i]] <- X[[i]] * sqrt(tabw[i])
    Xsepan <- sepan(X, nf = 4)
    normaliserparbloc <- function(scorcol) {
        for (i in 1:nbloc) {
            w1 <- scorcol[indicablo == i]
            w2 <- sqrt(sum(w1 * w1))
            if (w2 > tol) 
                w1 <- w1/w2
            scorcol[indicablo == i] <- w1
        }
        return(scorcol)
    }
    recalculer <- function(tab, scorcol) {
        for (k in 1:nbloc) {
            soustabk <- tab[, indicablo == k]
            uk <- scorcol[indicablo == k]
            wk <- cw[indicablo == k]
            soustabk.hat <- t(apply(soustabk, 1, function(x) sum(x * 
                uk) * uk))
            soustabk <- soustabk - soustabk.hat
            tab[, indicablo == k] <- soustabk
        }
        return(tab)
    }
    tab <- as.matrix(X[[1]])
    for (i in 2:nbloc) {
        tab <- cbind(tab, X[[i]])
    }
    names(tab) <- auxinames$col
    tab <- tab * sqrt(lw)
    tab <- t(t(tab) * sqrt(cw))
    compogene <- list()
    uknorme <- list()
    valsing <- NULL
    nfprovi <- min(c(20, nlig, ncol))
    for (i in 1:nfprovi) {
        af <- svd(tab)
        w <- af$u[, 1]
        w <- w/sqrt(lw)
        compogene[[i]] <- w
        w <- af$v[, 1]
        w <- normaliserparbloc(w)
        tab <- recalculer(tab, w)
        w <- w/sqrt(cw)
        uknorme[[i]] <- w
        w <- af$d[1]
        valsing <- c(valsing, w)
    }
    pseudoeig <- valsing^2
    if (scannf) {
        barplot(pseudoeig)
        cat("Select the number of axes: ")
        nf <- as.integer(readLines(n = 1))
    }
    if (nf <= 0) 
        nf <- 2
    acom <- list()
    acom$pseudoeig <- pseudoeig
    w <- matrix(0, nbloc, nf)
    for (i in 1:nbloc) {
        w1 <- Xsepan$Eig[rank.fac == i]
        r0 <- Xsepan$rank[i]
        if (r0 > nf) 
            r0 <- nf
        w[i, 1:r0] <- w1[1:r0]
    }
    w <- data.frame(w)
    row.names(w) <- Xsepan$tab.names
    names(w) <- paste("lam", 1:nf, sep = "")
    acom$lambda <- w
    w <- matrix(0, nlig, nf)
    for (j in 1:nf) w[, j] <- compogene[[j]]
    w <- data.frame(w)
    names(w) <- paste("SynVar", 1:nf, sep = "")
    row.names(w) <- row.names(X)
    acom$SynVar <- w
    w <- matrix(0, ncol, nf)
    for (j in 1:nf) w[, j] <- uknorme[[j]]
    w <- data.frame(w)
    names(w) <- paste("Axis", 1:nf, sep = "")
    row.names(w) <- auxinames$col
    acom$axis <- w
    w <- matrix(0, nlig * nbloc, nf)
    covar <- matrix(0, nbloc, nf)
    i1 <- 0
    i2 <- 0
    for (k in 1:nbloc) {
        i1 <- i2 + 1
        i2 <- i2 + nlig
        urk <- as.matrix(acom$axis[indicablo == k, ])
        tab <- as.matrix(X[[k]])
        urk <- urk * cw[indicablo == k]
        urk <- tab %*% urk
        w[i1:i2, ] <- urk
        urk <- urk * acom$SynVar * lw
        covar[k, ] <- apply(urk, 2, sum)
    }
    w <- data.frame(w, row.names = auxinames$row)
    names(w) <- paste("Axis", 1:nf, sep = "")
    acom$Tli <- w
    covar <- data.frame(covar)
    row.names(covar) <- tab.names(X)
    names(covar) <- paste("cov2", 1:nf, sep = "")
    acom$cov2 <- covar^2
    w <- matrix(0, nlig * nbloc, nf)
    i1 <- 0
    i2 <- 0
    for (k in 1:nbloc) {
        i1 <- i2 + 1
        i2 <- i2 + nlig
        tab <- acom$Tli[i1:i2, ]
        tab <- scalewt(tab, wt = lw, center = FALSE, scale = TRUE)
        w[i1:i2, ] <- tab
    }
    w <- data.frame(w, row.names = auxinames$row)
    names(w) <- paste("Axis", 1:nf, sep = "")
    acom$Tl1 <- w
    w <- matrix(0, ncol, nf)
    i1 <- 0
    i2 <- 0
    for (k in 1:nbloc) {
        i1 <- i2 + 1
        i2 <- i2 + ncol(X[[k]])
        urk <- as.matrix(acom$SynVar)
        tab <- as.matrix(X[[k]])
        urk <- urk * lw
        w[i1:i2, ] <- t(tab) %*% urk
    }
    w <- data.frame(w, row.names = auxinames$col)
    names(w) <- paste("SV", 1:nf, sep = "")
    acom$Tco <- w
    var.names <- NULL
    w <- matrix(0, nbloc * 4, nf)
    i1 <- 0
    i2 <- 0
    for (k in 1:nbloc) {
        i1 <- i2 + 1
        i2 <- i2 + 4
        urk <- as.matrix(acom$axis[indicablo == k, ])
        tab <- as.matrix(Xsepan$C1[indicablo == k, ])
        urk <- urk * cw[indicablo == k]
        tab <- t(tab) %*% urk
        for (i in 1:min(nf, 4)) {
            if (tab[i, i] < 0) {
                for (j in 1:nf) tab[i, j] <- -tab[i, j]
            }
        }
        w[i1:i2, ] <- tab
        var.names <- c(var.names, paste(Xsepan$tab.names[k], 
            ".a", 1:4, sep = ""))
    }
    w <- data.frame(w, row.names = auxinames$tab)
    names(w) <- paste("Axis", 1:nf, sep = "")
    acom$Tax <- w
    acom$nf <- nf
    acom$TL <- X$TL
    acom$TC <- X$TC
    acom$T4 <- X$T4
    class(acom) <- "mcoa"
    acom$call <- match.call()
    return(acom)
}
"mfa" <-
function (X, option = c("lambda1", "inertia", "uniform", "internal"), 
    scannf = TRUE, nf = 3) 
{
    if (!inherits(X, "ktab")) 
        stop("object 'ktab' expected")
    if (option == "internal") {
        if (is.null(X$tabw)) {
            warning("Internal weights not found: uniform weigths are used")
            option <- "uniform"
        }
    }
    lw <- X$lw
    cw <- X$cw
    sepan <- sepan(X, nf = 4)
    nbloc <- length(sepan$blo)
    indicablo <- factor(rep(1:nbloc, sepan$blo))
    rank.fac <- factor(rep(1:nbloc, sepan$rank))
    ncw <- NULL
    tab.names <- names(X)[1:nbloc]
    auxinames <- ktab.util.names(X)
    if (option == "lambda1") {
        for (i in 1:nbloc) {
            ncw <- c(ncw, rep(1/sepan$Eig[rank.fac == i][1], 
                sepan$blo[i]))
        }
    }
    else if (option == "inertia") {
        for (i in 1:nbloc) {
            ncw <- c(ncw, rep(1/sum(sepan$Eig[rank.fac == i]), 
                sepan$blo[i]))
        }
    }
    else if (option == "uniform") 
        ncw <- rep(1, sum(sepan$blo))
    else if (option == "internal") 
        ncw <- rep(X$tabw, sepan$blo)
    else stop("unknown option")
    ncw <- cw * ncw
    tab <- X[[1]]
    for (i in 2:nbloc) {
        tab <- cbind.data.frame(tab, X[[i]])
    }
    names(tab) <- auxinames$col
    anaco <- as.dudi(tab, col.w = ncw, row.w = lw, nf = nf, scannf = scannf, 
        call = match.call(), type = "mfa")
    nf <- anaco$nf
    afm <- list()
    afm$tab.names <- names(X)[1:nbloc]
    afm$blo <- X$blo
    afm$TL <- X$TL
    afm$TC <- X$TC
    afm$T4 <- X$T4
    afm$tab <- anaco$tab
    afm$eig <- anaco$eig
    afm$rank <- anaco$rank
    afm$li <- anaco$li
    afm$l1 <- anaco$l1
    afm$nf <- anaco$nf
    afm$lw <- anaco$lw
    afm$cw <- anaco$cw
    afm$co <- anaco$co
    afm$c1 <- anaco$c1
    projiner <- function(xk, qk, d, z) {
        w7 <- t(as.matrix(xk) * d) %*% as.matrix(z)
        iner <- apply(w7 * w7 * qk, 2, sum)
        return(iner)
    }
    link <- matrix(0, nbloc, nf)
    for (k in 1:nbloc) {
        xk <- X[[k]]
        q <- ncw[indicablo == k]
        link[k, ] <- projiner(xk, q, lw, anaco$l1)
    }
    link <- as.data.frame(link)
    names(link) <- paste("Comp", 1:nf, sep = "")
    row.names(link) <- tab.names
    afm$link <- link
    w <- matrix(0, nbloc * 4, nf)
    i1 <- 0
    i2 <- 0
    matl1 <- as.matrix(afm$l1)
    for (k in 1:nbloc) {
        i1 <- i2 + 1
        i2 <- i2 + 4
        tab <- as.matrix(sepan$L1[sepan$TL[, 1] == k, ])
        if (ncol(tab) > 4) 
            tab <- tab[, 1:4]
        if (ncol(tab) < 4) 
            tab <- cbind(tab, matrix(0, nrow(tab), 4 - ncol(tab)))
        tab <- t(tab * lw) %*% matl1
        for (i in 1:min(nf, 4)) {
            if (tab[i, i] < 0) {
                for (j in 1:nf) tab[i, j] <- -tab[i, j]
            }
        }
        w[i1:i2, ] <- tab
    }
    w <- data.frame(w)
    names(w) <- paste("Comp", 1:nf, sep = "")
    row.names(w) <- auxinames$tab
    afm$T4comp <- w
    w <- matrix(0, nrow(sepan$TL), ncol = nf)
    i1 <- 0
    i2 <- 0
    for (k in 1:nbloc) {
        i1 <- i2 + 1
        i2 <- i2 + length(lw)
        qk <- ncw[indicablo == k]
        xk <- as.matrix(X[[k]])
        w[i1:i2, ] <- (xk %*% (qk * t(xk))) %*% (matl1 * lw)
    }
    w <- data.frame(w)
    row.names(w) <- auxinames$row
    names(w) <- paste("Fac", 1:nf, sep = "")
    afm$lisup <- w
    afm$tabw <- X$tabw
    afm$call <- match.call()
    class(afm) <- c("mfa", "list")
    return(afm)
}
"plot.foucart" <-
function (x, xax = 1, yax = 2, clab = 1, csub = 2, possub = "bottomright", ...) 
{
    if (!inherits(x, "foucart")) 
        stop("Object of type 'foucart' expected")
    opar <- par(ask = par("ask"), mfrow = par("mfrow"), mar = par("mar"))
    on.exit(par(opar))
    par(mfrow = c(2, 2))
    nf <- ncol(x$li)
    cootot <- x$li[, c(xax, yax)]
    auxi <- x$li[, c(xax, yax)]
    names(auxi) <- names(cootot)
    cootot <- rbind.data.frame(cootot, auxi)
    auxi <- x$Tli[, c(xax, yax)]
    names(auxi) <- names(cootot)
    cootot <- rbind.data.frame(cootot, auxi)
    auxi <- x$Tco[, c(xax, yax)]
    names(auxi) <- names(cootot)
    cootot <- rbind.data.frame(cootot, auxi)
    s.label(cootot, clab = 0, cpoi = 0, sub = "Rows (Base)", 
        csub = csub, possub = possub)
    s.label(x$li, xax, yax, clab = clab, add.p = TRUE)
    s.label(cootot, clab = 0, cpoi = 0, sub = "Columns (Base)", 
        csub = csub, possub = possub)
    s.label(x$co, xax, yax, clab = clab, add.p = TRUE)
    s.label(cootot, clab = 0, cpoi = 0, sub = "Rows", csub = csub, 
        possub = possub)
    s.class(x$Tli, x$TL[, 2], xax = xax, yax = yax, 
        axesell = FALSE, clab = clab, add.p = TRUE)
    s.label(cootot, clab = 0, cpoi = 0, sub = "Columns", 
        csub = csub, possub = possub)
    s.class(x$Tco, x$TC[, 2], xax = xax, yax = yax, 
        axesell = FALSE, clab = clab, add.p = TRUE)
}
"plot.mcoa" <-
function (x, xax = 1, yax = 2, eig.bottom = TRUE, ...) 
{
    if (!inherits(x, "mcoa")) 
        stop("Object of type 'mcoa' expected")
    nf <- x$nf
    if (xax > nf) 
        stop("Non convenient xax")
    if (yax > nf) 
        stop("Non convenient yax")
    opar <- par(mar = par("mar"), mfrow = par("mfrow"), xpd = par("xpd"))
    on.exit(par(opar))
    par(mfrow = c(2, 2))
    coolig <- x$SynVar[, c(xax, yax)]
    for (k in 2:nrow(x$cov2)) {
        coolig <- rbind.data.frame(coolig, x$SynVar[, c(xax, 
            yax)])
    }
    names(coolig) <- names(x$Tl1)[c(xax, yax)]
    row.names(coolig) <- row.names(x$Tl1)
    s.match(x$Tl1[, c(xax, yax)], coolig, clab = 0, 
        sub = "Row projection", csub = 1.5, edge = FALSE)
    s.label(x$SynVar[, c(xax, yax)], add.plot = TRUE)
    coocol <- x$Tco[, c(xax, yax)]
    s.arrow(coocol, sub = "Col projection", csub = 1.5)
    valpr <- function(x) {
        opar <- par(mar = par("mar"))
        on.exit(par(opar))
        born <- par("usr")
        w <- x$pseudoeig
        col <- rep(grey(1), length(w))
        col[1:nf] <- grey(0.8)
        col[c(xax, yax)] <- grey(0)
        l0 <- length(w)
        xx <- seq(born[1], born[1] + (born[2] - born[1]) * l0/60, 
            le = l0 + 1)
        w <- w/max(w)
        w <- w * (born[4] - born[3])/4
        par(mar = c(0.1, 0.1, 0.1, 0.1))
        if (eig.bottom) 
            m3 <- born[3]
        else m3 <- born[4] - w[1]
        w <- m3 + w
        rect(xx[1], m3, xx[l0 + 1], w[1], col = grey(1))
        for (i in 1:l0) rect(xx[i], m3, xx[i + 1], w[i], col = col[i])
    }
    s.corcircle(x$Tax[x$T4[, 2] == 1, ], full = FALSE, 
        sub = "First axis projection", possub = "topright", csub = 1.5)
    valpr(x)
    plot(x$cov2[, c(xax, yax)])
    scatterutil.grid(0)
    title(main = "Pseudo-eigen values")
    par(xpd = TRUE)
    scatterutil.eti(x$cov2[, xax], x$cov2[, yax], label = row.names(x$cov2), 
        clabel = 1)
}
"plot.mfa" <-
function (x, xax = 1, yax = 2, option.plot = 1:4, ...) 
{
    if (!inherits(x, "mfa")) 
        stop("Object of type 'mfa' expected")
    nf <- x$nf
    if (xax > nf) 
        stop("Non convenient xax")
    if (yax > nf) 
        stop("Non convenient yax")
    opar <- par(mar = par("mar"), mfrow = par("mfrow"), xpd = par("xpd"))
    on.exit(par(opar))
    mfrow <- n2mfrow(length(option.plot))
    par(mfrow = mfrow)
    for (j in option.plot) {
        if (j == 1) {
            coolig <- x$lisup[, c(xax, yax)]
            s.class(coolig, fac = as.factor(x$TL[, 2]), 
                label = row.names(x$li), cell = 0, sub = "Row projection", 
                csub = 1.5)
            add.scatter.eig(x$eig, x$nf, xax, yax, posi = "top", 
                ratio = 1/5)
        }
        if (j == 2) {
            coocol <- x$co[, c(xax, yax)]
            s.arrow(coocol, sub = "Col projection", csub = 1.5)
            add.scatter.eig(x$eig, x$nf, xax, yax, posi = "top", 
                ratio = 1/5)
        }
        if (j == 3) {
            s.corcircle(x$T4comp[x$T4[, 2] == 1, ], 
                full = FALSE, sub = "Component projection", possub = "topright", 
                csub = 1.5)
            add.scatter.eig(x$eig, x$nf, xax, yax, posi = "bottom", 
                ratio = 1/5)
        }
        if (j == 4) {
            plot(x$link[, c(xax, yax)])
            scatterutil.grid(0)
            title(main = "Link")
            par(xpd = TRUE)
            scatterutil.eti(x$link[, xax], x$link[, yax], 
                label = row.names(x$link), clabel = 1)
        }
        if (j == 5) {
            scatterutil.eigen(x$eig, wsel = 1:x$nf, sub = "Eigen values", 
                csub = 2, possub = "topright")
        }
    }
}

"plot.pta" <-
function (x, xax = 1, yax = 2, option = 1:4, ...) 
{
    if (!inherits(x, "pta")) 
        stop("Object of type 'pta' expected")
    nf <- x$nf
    if (xax > nf) 
        stop("Non convenient xax")
    if (yax > nf) 
        stop("Non convenient yax")
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    mfrow <- n2mfrow(length(option))
    par(mfrow = mfrow)
    for (j in option) {
        if (j == 1) {
            coolig <- x$RV.coo[, c(1, 2)]
            s.corcircle(coolig, label = x$tab.names, 
                cgrid = 0, sub = "Interstructure", csub = 1.5, 
                possub = "topleft", full = TRUE)
            l0 <- length(x$RV.eig)
            add.scatter.eig(x$RV.eig, l0, 1, 2, posi = "bottom", 
                ratio = 1/4)
        }
        if (j == 2) {
            coolig <- x$li[, c(xax, yax)]
            s.label(coolig, sub = "Compromise", csub = 1.5, 
                possub = "topleft", )
            add.scatter.eig(x$eig, x$nf, xax, yax, posi = "bottom", 
                ratio = 1/4)
        }
        if (j == 3) {
            cooco <- x$co[, c(xax, yax)]
            s.arrow(cooco, sub = "Compromise", csub = 1.5, 
                possub = "topleft")
        }
        if (j == 4) {
            plot(x$tabw, x$cos2, xlab = "Tables weights", 
                ylab = "Cos 2")
            scatterutil.grid(0)
            title(main = "Typological value")
            par(xpd = TRUE)
            scatterutil.eti(x$tabw, x$cos2, label = x$tab.names, 
                clabel = 1)
        }
    }
}
"plot.sepan" <-
function (x, mfrow = NULL, csub = 2, ...) 
{
    if (!inherits(x, "sepan")) 
        stop("Object of type 'sepan' expected")
    opar <- par(ask = par("ask"), mfrow = par("mfrow"), mar = par("mar"))
    on.exit(par(opar))
    par(mar = c(0.6, 2.6, 0.6, 0.6))
    nbloc <- length(x$blo)
    if (is.null(mfrow)) 
        mfrow <- n2mfrow(nbloc)
    par(mfrow = mfrow)
    if (nbloc > prod(mfrow)) 
        par(ask = TRUE)
    rank.fac <- factor(rep(1:nbloc, x$rank))
    nf <- ncol(x$Li)
    neig <- max(x$rank)
    appel <- as.list(x$call)
    X <- eval(appel$X, sys.frame(0))
    maxeig <- max(x$Eig)
    for (ianal in 1:nbloc) {
        w <- x$Eig[rank.fac == ianal]
        scatterutil.eigen(w, xmax = neig, ymax = maxeig, wsel = 1:nf, 
            sub = x$tab.names[ianal], csub = csub, possub = "topright")
    }
}
"plot.statis" <-
function (x, xax = 1, yax = 2, option = 1:4, ...) 
{
    if (!inherits(x, "statis")) 
        stop("Object of type 'statis' expected")
    nf <- x$C.nf
    if (xax > nf) 
        stop("Non convenient xax")
    if (yax > nf) 
        stop("Non convenient yax")
    opar <- par(mar = par("mar"), mfrow = par("mfrow"), xpd = par("xpd"))
    on.exit(par(opar))
    mfrow <- n2mfrow(length(option))
    par(mfrow = mfrow)
    for (j in option) {
        if (j == 1) {
            coolig <- x$RV.coo[, c(1, 2)]
            s.corcircle(coolig, label = x$tab.names, 
                cgrid = 0, sub = "Interstructure", csub = 1.5, 
                possub = "topleft", full = TRUE)
            l0 <- length(x$RV.eig)
            add.scatter.eig(x$RV.eig, l0, 1, 2, posi = "bottom", 
                ratio = 1/4)
        }
        if (j == 2) {
            coolig <- x$C.li[, c(xax, yax)]
            s.label(coolig, sub = "Compromise", csub = 1.5, 
                possub = "topleft", )
            add.scatter.eig(x$C.eig, x$C.nf, xax, yax, 
                posi = "bottom", ratio = 1/4)
        }
        if (j == 4) {
            cooax <- x$C.T4[x$T4[, 2] == 1, ]
            s.corcircle(cooax, xax, yax, full = TRUE, sub = "Component projection", 
                possub = "topright", csub = 1.5)
            add.scatter.eig(x$C.eig, x$C.nf, xax, yax, 
                posi = "bottom", ratio = 1/5)
        }
        if (j == 3) {
            plot(x$RV.tabw, x$cos2, xlab = "Tables weights", 
                ylab = "Cos 2")
            scatterutil.grid(0)
            title(main = "Typological value")
            par(xpd = TRUE)
            scatterutil.eti(x$RV.tabw, x$cos2, label = x$tab.names, 
                clabel = 1)
        }
    }
}

"print.foucart" <-
function (x, ...) 
{
    cat("Foucart's  COA\n")
    cat("class: ")
    cat(class(x))
    cat("\n$call: ")
    print(x$call)
    cat("table  number:", length(x$blo), "\n")
    cat("\n$nf:", x$nf, "axis-components saved")
    cat("\n$rank: ")
    cat(x$rank)
    cat("\neigen values: ")
    l0 <- length(x$eig)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n")
    else cat("\n")
    cat("blo	vector		", length(x$blo), "		blocks\n")
    sumry <- array("", c(3, 4), list(rep("", 3), c("vector", 
        "length", "mode", "content")))
    sumry[1, ] <- c("$cw", length(x$cw), mode(x$cw), "column weights")
    sumry[2, ] <- c("$lw", length(x$lw), mode(x$lw), "row weights")
    sumry[3, ] <- c("$eig", length(x$eig), mode(x$eig), "eigen values")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(5, 4), list(rep("", 5), c("data.frame", 
        "nrow", "ncol", "content")))
    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "modified array")
    sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "row coordinates")
    sumry[3, ] <- c("$l1", nrow(x$l1), ncol(x$l1), "row normed scores")
    sumry[4, ] <- c("$co", nrow(x$co), ncol(x$co), "column coordinates")
    sumry[5, ] <- c("$c1", nrow(x$c1), ncol(x$c1), "column normed scores")
    class(sumry) <- "table"
    print(sumry)
    cat("\n     **** Intrastructure ****\n\n")
    sumry <- array("", c(4, 4), list(rep("", 4), c("data.frame", 
        "nrow", "ncol", "content")))
    sumry[1, ] <- c("$Tli", nrow(x$Tli), ncol(x$Tli), "row coordinates (each table)")
    sumry[2, ] <- c("$Tco", nrow(x$Tco), ncol(x$Tco), "col coordinates (each table)")
    sumry[3, ] <- c("$TL", nrow(x$TL), ncol(x$TL), "factors for Tli")
    sumry[4, ] <- c("$TC", nrow(x$TC), ncol(x$TC), "factors for Tco")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
}
"print.ktab" <-
function (x, ...) 
{
    if (!inherits(x, "ktab")) 
        stop("to be used with 'ktab' object")
    cat("class:", class(x), "\n")
    ntab <- length(x$blo)
    cat("\ntab number:	", ntab, "\n")
    sumry <- array("", c(ntab, 3), list(1:ntab, c("data.frame", 
        "nrow", "ncol")))
    for (i in 1:ntab) {
        sumry[i, ] <- c(names(x)[i], nrow(x[[i]]), ncol(x[[i]]))
    }
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(4, 4), list((ntab + 1):(ntab + 4), c("vector", 
        "length", "mode", "content")))
    sumry[1, ] <- c("$lw", length(x$lw), mode(x$lw), "row weigths")
    sumry[2, ] <- c("$cw", length(x$cw), mode(x$cw), "column weights")
    sumry[3, ] <- c("$blo", length(x$blo), mode(x$blo), "column numbers")
    sumry[4, ] <- c("$tabw", length(x$tabw), mode(x$tabw), "array weights")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(3, 4), list((ntab + 5):(ntab + 7), c("data.frame", 
        "nrow", "ncol", "content")))
    sumry[1, ] <- c("$TL", nrow(x$TL), ncol(x$TL), "Factors Table number Line number")
    sumry[2, ] <- c("$TC", nrow(x$TC), ncol(x$TC), "Factors Table number Col number")
    sumry[3, ] <- c("$T4", nrow(x$T4), ncol(x$T4), "Factors Table number 1234")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    cat((ntab + 8), "$call: ")
    print(x$call)
    cat("\n")
    cat("names :\n")
    for (i in 1:ntab) {
        cat(names(x)[i], ":", names(x[[i]]), "\n")
    }
    cat("\n")
    indica <- as.factor(rep(1:ntab, x$blo))
    w <- split(x$cw, indica)
    cat("Col weigths :\n")
    for (i in 1:ntab) {
        cat(names(x)[i], ":", w[[i]], "\n")
    }
    cat("\n")
    cat("Row weigths :\n")
    cat(x$lw)
    cat("\n")
}
"print.mcoa" <-
function (x, ...) 
{
    if (!inherits(x, "mcoa")) 
        stop("non convenient data")
    cat("Multiple Co-inertia Analysis\n")
    cat(paste("list of class", class(x)))
    l0 <- length(x$pseudoeig)
    cat("\n\n$pseudoeig:", l0, "pseudo eigen values\n")
    cat(signif(x$pseudoeig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n")
    else cat("\n")
    cat("\n$call: ")
    print(x$call)
    cat("\n$nf:", x$nf, "axis saved\n\n")
    sumry <- array("", c(11, 4), list(1:11, c("data.frame", "nrow", 
        "ncol", "content")))
    sumry[1, ] <- c("$SynVar", nrow(x$SynVar), ncol(x$SynVar), 
        "synthetic scores")
    sumry[2, ] <- c("$axis", nrow(x$axis), ncol(x$axis), 
        "co-inertia axis")
    sumry[3, ] <- c("$Tli", nrow(x$Tli), ncol(x$Tli), "co-inertia coordinates")
    sumry[4, ] <- c("$Tl1", nrow(x$Tl1), ncol(x$Tl1), "co-inertia normed scores")
    sumry[5, ] <- c("$Tax", nrow(x$Tax), ncol(x$Tax), "inertia axes onto co-inertia axis")
    sumry[6, ] <- c("$Tco", nrow(x$Tco), ncol(x$Tco), "columns onto synthetic scores")
    sumry[7, ] <- c("$TL", nrow(x$TL), ncol(x$TL), "factors for Tli Tl1")
    sumry[8, ] <- c("$TC", nrow(x$TC), ncol(x$TC), "factors for Tco")
    sumry[9, ] <- c("$T4", nrow(x$T4), ncol(x$T4), "factors for Tax")
    sumry[10, ] <- c("$lambda", nrow(x$lambda), ncol(x$lambda), 
        "eigen values (separate analysis)")
    sumry[11, ] <- c("$cov2", nrow(x$cov2), ncol(x$cov2), 
        "pseudo eigen values (synthetic analysis)")
    class(sumry) <- "table"
    print(sumry)
    cat("other elements: ")
    if (length(names(x)) > 14) 
        cat(names(x)[15:(length(x))], "\n")
    else cat("NULL\n")
}
"print.mfa" <-
function (x, ...) 
{
    if (!inherits(x, "mfa")) 
        stop("non convenient data")
    cat("Multiple Factorial Analysis\n")
    cat(paste("list of class", class(mfa)))
    cat("\n$call: ")
    print(x$call)
    cat("$nf:", x$nf, "axis-components saved\n\n")
    sumry <- array("", c(6, 4), list(1:6, c("vector", "length", 
        "mode", "content")))
    sumry[1, ] <- c("$tab.names", length(x$tab.names), mode(x$tab.names), 
        "tab names")
    sumry[2, ] <- c("$blo", length(x$blo), mode(x$blo), "column number")
    sumry[3, ] <- c("$rank", length(x$rank), mode(x$rank), 
        "tab rank")
    sumry[4, ] <- c("$eig", length(x$eig), mode(x$eig), "eigen values")
    sumry[5, ] <- c("$lw", length(x$lw), mode(x$lw), "row weights")
    sumry[6, ] <- c("$tabw", length(x$tabw), mode(x$tabw), 
        "array weights")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(11, 4), list(1:11, c("data.frame", "nrow", 
        "ncol", "content")))
    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "modified array")
    sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "row coordinates")
    sumry[3, ] <- c("$l1", nrow(x$l1), ncol(x$l1), "row normed scores")
    sumry[4, ] <- c("$co", nrow(x$co), ncol(x$co), "column coordinates")
    sumry[5, ] <- c("$c1", nrow(x$c1), ncol(x$c1), "column normed scores")
    sumry[6, ] <- c("$lisup", nrow(x$lisup), ncol(x$lisup), 
        "row coordinates from each table")
    sumry[7, ] <- c("$TL", nrow(x$TL), ncol(x$TL), "factors for li l1")
    sumry[8, ] <- c("$TC", nrow(x$TC), ncol(x$TC), "factors for co c1")
    sumry[9, ] <- c("$T4", nrow(x$T4), ncol(x$T4), "factors for T4comp")
    sumry[10, ] <- c("$T4comp", nrow(x$T4comp), ncol(x$T4comp), 
        "component projection")
    sumry[11, ] <- c("$link", nrow(x$link), ncol(x$link), 
        "link array-total")
    class(sumry) <- "table"
    print(sumry)
    cat("other elements: ")
    if (length(names(mcoa)) > 19) 
        cat(names(mcoa)[20:(length(mcoa))], "\n")
    else cat("NULL\n")
}

"print.pta" <-
function (x, ...) 
{
    cat("Partial Triadic Analysis\n")
    cat("class:")
    cat(class(x), "\n")
    cat("table number:", length(x$blo), "\n")
    cat("row number:", length(x$lw), "	column number:", length(x$cw), 
        "\n")
    cat("\n     **** Interstructure ****\n")
    cat("\neigen values: ")
    l0 <- length(x$RV.eig)
    cat(signif(x$RV.eig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n")
    else cat("\n")
    cat(" $RV		matrix		", nrow(x$RV), "	", ncol(x$RV), "	RV coefficients\n")
    cat(" $RV.eig	vector		", length(x$RV.eig), "		eigenvalues\n")
    cat(" $RV.coo	data.frame	", nrow(x$RV.coo), "	", ncol(x$RV.coo), 
        "	array scores\n")
    cat(" $tab.names	vector		", length(x$tab.names), "		array names\n")
    cat("\n      **** Compromise ****\n")
    cat("\neigen values: ")
    l0 <- length(x$eig)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n")
    else cat("\n")
    cat("\n $nf:", x$nf, "axis-components saved")
    cat("\n $rank: ")
    cat(x$rank, "\n\n")
    sumry <- array("", c(5, 4), list(rep("", 5), c("vector", 
        "length", "mode", "content")))
    sumry[1, ] <- c("$tabw", length(x$tabw), mode(x$tabw), "array weights")
    sumry[2, ] <- c("$cw", length(x$cw), mode(x$cw), "column weights")
    sumry[3, ] <- c("$lw", length(x$lw), mode(x$lw), "row weights")
    sumry[4, ] <- c("$eig", length(x$eig), mode(x$eig), "eigen values")
    sumry[5, ] <- c("$cos2", length(x$cos2), mode(x$cos2), "cosine^2 between compromise and arrays")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(5, 4), list(rep("", 5), c("data.frame", 
        "nrow", "ncol", "content")))
    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "modified array")
    sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "row coordinates")
    sumry[3, ] <- c("$l1", nrow(x$l1), ncol(x$l1), "row normed scores")
    sumry[4, ] <- c("$co", nrow(x$co), ncol(x$co), "column coordinates")
    sumry[5, ] <- c("$c1", nrow(x$c1), ncol(x$c1), "column normed scores")
    class(sumry) <- "table"
    print(sumry)
    cat("\n     **** Intrastructure ****\n\n")
    sumry <- array("", c(7, 4), list(rep("", 7), c("data.frame", 
        "nrow", "ncol", "content")))
    sumry[1, ] <- c("$Tli", nrow(x$Tli), ncol(x$Tli), "row coordinates (each table)")
    sumry[2, ] <- c("$Tco", nrow(x$Tco), ncol(x$Tco), "col coordinates (each table)")
    sumry[3, ] <- c("$Tcomp", nrow(x$Tcomp), ncol(x$Tcomp), "principal components (each table)")
    sumry[4, ] <- c("$Tax", nrow(x$Tax), ncol(x$Tax), "principal axis (each table)")
    sumry[5, ] <- c("$TL", nrow(x$TL), ncol(x$TL), "factors for Tli")
    sumry[6, ] <- c("$TC", nrow(x$TC), ncol(x$TC), "factors for Tco")
    sumry[7, ] <- c("$T4", nrow(x$T4), ncol(x$T4), "factors for Tax Tcomp")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
}
"print.sepan" <-
function (x, ...) 
{
    if (!inherits(x, "sepan")) 
        stop("to be used with 'sepan' object")
    cat("class:", class(x), "\n")
    cat("$call: ")
    print(x$call)
    sumry <- array("", c(4, 4), list(1:4, c("vector", "length", 
        "mode", "content")))
    sumry[1, ] <- c("$tab.names", length(x$tab.names), mode(x$tab.names), 
        "tab names")
    sumry[2, ] <- c("$blo", length(x$blo), mode(x$blo), "column number")
    sumry[3, ] <- c("$rank", length(x$rank), mode(x$rank), "tab rank")
    sumry[4, ] <- c("$Eig", length(x$Eig), mode(x$Eig), "All the eigen values")
    class(sumry) <- "table"
    print(sumry)
    sumry <- array("", c(6, 4), list(1:6, c("data.frame", "nrow", 
        "ncol", "content")))
    sumry[1, ] <- c("$Li", nrow(x$Li), ncol(x$Li), "row coordinates")
    sumry[2, ] <- c("$L1", nrow(x$L1), ncol(x$L1), "row normed scores")
    sumry[3, ] <- c("$Co", nrow(x$Co), ncol(x$Co), "column coordinates")
    sumry[4, ] <- c("$C1", nrow(x$C1), ncol(x$C1), "column normed coordinates")
    sumry[5, ] <- c("$TL", nrow(x$TL), ncol(x$TL), "factors for Li L1")
    sumry[6, ] <- c("$TC", nrow(x$TC), ncol(x$TC), "factors for Co C1")
    class(sumry) <- "table"
    print(sumry)
}
"print.statis" <-
function (x, ...) 
{
    cat("STATIS Analysis\n")
    cat("class:")
    cat(class(x), "\n")
    cat("table number:", length(x$RV.tabw), "\n")
    cat("row number:", nrow(x$C.li), "	total column number:", 
        nrow(x$C.Co), "\n")
    cat("\n     **** Interstructure ****\n")
    cat("\neigen values: ")
    l0 <- length(x$RV.eig)
    cat(signif(x$RV.eig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n")
    else cat("\n")
    cat(" $RV		matrix		", nrow(x$RV), "	", ncol(x$RV), "	RV coefficients\n")
    cat(" $RV.eig	vector		", length(x$RV.eig), "		eigenvalues\n")
    cat(" $RV.coo	data.frame	", nrow(x$RV.coo), "	", ncol(x$RV.coo), 
        "	array scores\n")
    cat(" $tab.names	vector		", length(x$tab.names), "		array names\n")
    cat(" $RV.tabw	vector		", length(x$RV.tabw), "		array weigths\n")
    cat("\nRV coefficient\n")
    w <- x$RV
    w[row(w) < col(w)] <- NA
    print(w, na = "")
    cat("\n      **** Compromise ****\n")
    cat("\neigen values: ")
    l0 <- length(x$C.eig)
    cat(signif(x$C.eig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n")
    else cat("\n")
    cat("\n $nf:", x$C.nf, "axis-components saved")
    cat("\n $rank: ")
    cat(x$C.rank, "\n")
    sumry <- array("", c(6, 4), list(rep("", 6), c("data.frame", 
        "nrow", "ncol", "content")))
    sumry[1, ] <- c("$C.li", nrow(x$C.li), ncol(x$C.li), "row coordinates")
    sumry[2, ] <- c("$C.Co", nrow(x$C.Co), ncol(x$C.Co), "column coordinates")
    sumry[3, ] <- c("$T4", nrow(x$T4), ncol(x$T4), "principal vectors (each table)")
    sumry[4, ] <- c("$TL", nrow(x$TL), ncol(x$TL), "factors (not used)")
    sumry[5, ] <- c("$TC", nrow(x$TC), ncol(x$TC), "factors for Co")
    sumry[6, ] <- c("$T4", nrow(x$T4), ncol(x$T4), "factors for T4")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
}

"pta" <-
function (X, scannf = TRUE, nf = 2) 
{
    if (!inherits(X, "ktab")) 
        stop("object 'ktab' expected")
    auxinames <- ktab.util.names(X)
    sepa <- sepan(X, nf = 4)
    blocks <- X$blo
    nblo <- length(blocks)
    tnames <- tab.names(X)
    lw <- X$lw
    lwsqrt <- sqrt(X$lw)
    nl <- length(lw)
    r.n <- row.names(X[[1]])
    for (i in 1:nblo) {
        r.new <- row.names(X[[i]])
        if (any(r.new != r.n)) 
            stop("non equal row.names among array")
    }
    if (length(unique(blocks)) != 1) 
        stop("non equal col numbers among array")
    unique.col.names <- names(X[[1]])
    for (i in 1:nblo) {
        c.new <- names(X[[i]])
        if (any(c.new != unique.col.names)) 
            stop("non equal col.names among array")
    }
    indica <- as.factor(rep(1:nblo, blocks))
    w <- split(X$cw, indica)
    cw <- w[[1]]
    for (i in 1:nblo) {
        col.w.new <- w[[i]]
        if (any(cw != col.w.new)) 
            stop("non equal column weights among array")
    }
    cwsqrt <- sqrt(cw)
    nc <- length(cw)
    atp <- list()
    for (i in 1:nblo) {
        w <- as.matrix(X[[i]]) * lwsqrt
        w <- t(t(w) * cwsqrt)
        atp[[i]] <- w
    }
    atp <- matrix(unlist(atp), nl * nc, nblo)
    RV <- t(atp) %*% atp
    ak <- sqrt(diag(RV))
    RV <- sweep(RV, 1, ak, "/")
    RV <- sweep(RV, 2, ak, "/")
    dimnames(RV) <- list(tnames, tnames)
    atp <- list()
    inter <- eigen(as.matrix(RV))
    if (any(inter$vectors[, 1] < 0)) 
        inter$vectors[, 1] <- -inter$vectors[, 1]
    is <- inter$vectors[, (1:min(c(nblo, 4)))]
    tabw <- as.vector(is[, 1])
    is <- t(t(is) * sqrt(inter$values[1:ncol(is)]))
    is <- as.data.frame(is)
    row.names(is) <- tnames
    names(is) <- paste("IS", 1:ncol(is), sep = "")
    atp$RV <- RV
    atp$RV.eig <- inter$values
    atp$RV.coo <- is
    atp$tabw <- tabw
    tab <- X[[1]] * tabw[1]
    for (i in 2:nblo) {
        tab <- tab + X[[i]] * tabw[i]
    }
    tab <- as.data.frame(tab, row.names = row.names(X))
    names(tab) <- unique.col.names
    comp <- as.dudi(tab, col.w = cw, row.w = lw, nf = nf, scannf = scannf, 
        call = match.call(), type = "pta")
    atp$rank <- comp$rank
    nf <- atp$nf <- comp$nf
    atp$tab <- comp$tab
    atp$lw <- comp$lw
    atp$cw <- comp$cw
    atp$eig <- comp$eig
    atp$li <- comp$li
    atp$co <- comp$co
    atp$l1 <- comp$li
    atp$c1 <- comp$co
    w1 <- matrix(0, nblo * 4, nf)
    w2 <- matrix(0, nblo * 4, nf)
    i1 <- 0
    i2 <- 0
    for (k in 1:nblo) {
        i1 <- i2 + 1
        i2 <- i2 + 4
        tab1 <- as.matrix(sepa$L1[X$TL[, 1] == k, ])
        tab1 <- t(tab1 * lw) %*% as.matrix(comp$l1)
        tab2 <- as.matrix(sepa$C1[X$TC[, 1] == k, ])
        tab2 <- (t(tab2) * cw) %*% as.matrix(comp$c1)
        for (i in 1:min(nf, 4)) {
            if (tab2[i, i] < 0) {
                for (j in 1:nf) tab2[i, j] <- -tab2[i, j]
            }
            if (tab1[i, i] < 0) {
                for (j in 1:nf) tab1[i, j] <- -tab1[i, j]
            }
        }
        w1[i1:i2, ] <- tab1
        w2[i1:i2, ] <- tab2
    }
    w1 <- data.frame(w1, row.names = auxinames$tab)
    w2 <- data.frame(w2, row.names = auxinames$tab)
    names(w2) <- names(w1) <- paste("C", 1:nf, sep = "")
    atp$Tcomp <- w1
    atp$Tax <- w2
    tab <- as.matrix(X[[1]])
    w <- as.matrix(comp$c1)
    cooli <- t(t(tab) * cw) %*% w
    for (k in 2:nblo) {
        tab <- as.matrix(X[[k]])
        cooliauxi <- t(t(tab) * cw) %*% w
        cooli <- rbind(cooli, cooliauxi)
    }
    cooli <- data.frame(cooli, row.names = auxinames$row)
    atp$Tli <- cooli
    tab <- as.matrix(X[[1]])
    w <- as.matrix(comp$l1) * lw
    cooco <- t(tab) %*% w
    for (k in 2:nblo) {
        tab <- as.matrix(X[[k]])
        coocoauxi <- t(tab) %*% w
        cooco <- rbind(cooco, coocoauxi)
    }
    cooco <- data.frame(cooco, row.names = auxinames$col)
    atp$Tco <- cooco
    normcompro <- sum(atp$eig)
    indica <- as.factor(rep(1:nblo, sepa$rank))
    w <- split(sepa$Eig, indica)
    normtab <- unlist(lapply(w, sum))
    covv <- rep(0, nblo)
    w1 <- atp$tab * lwsqrt
    w1 <- t(t(w1) * cwsqrt)
    for (k in 1:nblo) {
        wk <- X[[k]] * lwsqrt
        wk <- t(t(wk) * cwsqrt)
        covv[k] <- sum(w1 * wk)
    }
    atp$cos2 <- covv/sqrt(normcompro)/sqrt(normtab)
    atp$TL <- X$TL
    atp$TC <- X$TC
    atp$T4 <- X$T4
    atp$blo <- X$blo
    atp$tab.names <- tnames
    atp$call <- match.call()
    class(atp) <- c("pta", "dudi")
    return(atp)
}

"row.names.ktab" <-
function (x) 
{
    if (!inherits(x, "ktab")) 
        stop("to be used with 'ktab' object")
    ntab <- length(x$blo)
    cha <- attr(x[[1]], "row.names")
    for (i in 1:ntab) {
        if (any(attr(x[[i]], "row.names") != cha)) 
            warnings(paste("array", i, "and array 1 have different row.names"))
    }
    return(cha)
}
"row.names<-.ktab" <-
function (x, value) 
{
    if (!inherits(x, "ktab")) 
        stop("to be used with 'ktab' object")
    ntab <- length(x$blo)
    old <- attr(x[[1]], "row.names")
    if (!is.null(old) && length(value) != length(old)) 
        stop("invalid row.names length")
    value <- as.character(value)
    if (any(duplicated(value))) 
        stop("duplicate row.names are not allowed")
    for (i in 1:ntab) {
        attr(x[[i]], "row.names") <- value
    }
    x
}

"sepan" <-
function (X, nf = 2) 
{
    if (!inherits(X, "ktab")) 
        stop("object 'ktab' expected")
    complete.dudi <- function(dudi, nf1, nf2) {
        pcolzero <- nf2 - nf1 + 1
        w <- data.frame(matrix(0, nrow(dudi$li), pcolzero))
        names(w) <- paste("Axis", (nf1:nf2), sep = "")
        dudi$li <- cbind.data.frame(dudi$li, w)
        w <- data.frame(matrix(0, nrow(dudi$li), pcolzero))
        names(w) <- paste("RS", (nf1:nf2), sep = "")
        dudi$l1 <- cbind.data.frame(dudi$l1, w)
        w <- data.frame(matrix(0, nrow(dudi$co), pcolzero))
        names(w) <- paste("Comp", (nf1:nf2), sep = "")
        wco <- data.frame(matrix(0, nrow(dudi$co), pcolzero))
        dudi$co <- cbind.data.frame(dudi$co, w)
        w <- data.frame(matrix(0, nrow(dudi$co), pcolzero))
        names(w) <- paste("CS", (nf1:nf2), sep = "")
        dudi$c1 <- cbind.data.frame(dudi$c1, w)
        return(dudi)
    }
    lw <- X$lw
    cw <- X$cw
    blo <- X$blo
    ntab <- length(blo)
    auxinames <- ktab.util.names(X)
    tab <- as.data.frame(X[[1]])
    j1 <- 1
    j2 <- as.numeric(blo[1])
    auxi <- as.dudi(tab, col.w = cw[j1:j2], row.w = lw, nf = nf, 
        scannf = FALSE, call = match.call(), type = "sepan")
    if (auxi$nf < nf) 
        auxi <- complete.dudi(auxi, auxi$nf + 1, nf)
    Eig <- auxi$eig
    Co <- auxi$co
    Li <- auxi$li
    C1 <- auxi$c1
    L1 <- auxi$l1
    row.names(Li) <- paste(row.names(Li), j1, sep = ".")
    row.names(L1) <- paste(row.names(L1), j1, sep = ".")
    row.names(Co) <- paste(row.names(Co), j1, sep = ".")
    row.names(C1) <- paste(row.names(C1), j1, sep = ".")
    rank <- auxi$rank
    for (i in 2:ntab) {
        j1 <- j2 + 1
        j2 <- j2 + as.numeric(blo[i])
        tab <- as.data.frame(X[[i]])
        auxi <- as.dudi(tab, col.w = cw[j1:j2], row.w = lw, nf = nf, 
            scannf = FALSE, call = match.call(), type = "sepan")
        Eig <- c(Eig, auxi$eig)
        row.names(auxi$li) <- paste(row.names(auxi$li), i, sep = ".")
        row.names(auxi$l1) <- paste(row.names(auxi$l1), i, sep = ".")
        row.names(auxi$co) <- paste(row.names(auxi$co), i, sep = ".")
        row.names(auxi$c1) <- paste(row.names(auxi$c1), i, sep = ".")
        if (auxi$nf < nf) 
            auxi <- complete.dudi(auxi, auxi$nf + 1, nf)
        Co <- rbind.data.frame(Co, auxi$co)
        Li <- rbind.data.frame(Li, auxi$li)
        C1 <- rbind.data.frame(C1, auxi$c1)
        L1 <- rbind.data.frame(L1, auxi$l1)
        rank <- c(rank, auxi$rank)
    }
    res <- list()
    res$Li <- Li
    res$L1 <- L1
    res$Co <- Co
    res$C1 <- C1
    res$Eig <- Eig
    res$TL <- X$TL
    res$TC <- X$TC
    res$blo <- blo
    res$rank <- rank
    res$tab.names <- names(X)[1:ntab]
    res$call <- match.call()
    class(res) <- c("sepan", "list")
    return(res)
}
"statis" <-
function (X, scannf = TRUE, nf = 3, tol = 1e-07) 
{
    if (!inherits(X, "ktab")) 
        stop("object 'ktab' expected")
    lw <- X$lw
    nlig <- length(lw)
    cw <- X$cw
    ncol <- length(cw)
    blo <- X$blo
    ntab <- length(X$blo)
    indicablo <- X$TC[, 1]
    tab.names <- tab.names(X)
    auxinames <- ktab.util.names(X)
    statis <- list()
    sep <- list()
    lwsqrt <- sqrt(lw)
    for (k in 1:ntab) {
        ak <- sqrt(cw[indicablo == k])
        wk <- as.matrix(X[[k]]) * lwsqrt
        wk <- t(t(wk) * ak)
        wk <- wk %*% t(wk)
        sep[[k]] <- wk
    }
    sep <- matrix(unlist(sep), nlig * nlig, ntab)
    RV <- t(sep) %*% sep
    ak <- sqrt(diag(RV))
    RV <- sweep(RV, 1, ak, "/")
    RV <- sweep(RV, 2, ak, "/")
    dimnames(RV) <- list(tab.names, tab.names)
    statis$RV <- RV
    eig1 <- La.eigen(RV, sym = TRUE)
    statis$RV.eig <- eig1$values
    if (any(eig1$vectors[, 1] < 0)) 
        eig1$vectors[, 1] <- -eig1$vectors[, 1]
    tabw <- eig1$vectors[, 1]
    statis$RV.tabw <- tabw
    w <- t(t(eig1$vectors) * sqrt(eig1$values))
    w <- as.data.frame(w)
    row.names(w) <- tab.names
    names(w) <- paste("S", 1:ncol(w), sep = "")
    statis$RV.coo <- w[, 1:min(4, ncol(w))]
    sep <- t(t(sep)/ak)
    C.ro <- apply(t(sep) * tabw, 2, sum)
    C.ro <- matrix(unlist(C.ro), nlig, nlig)
    eig1 <- La.eigen(C.ro, sym = TRUE)
    eig <- eig1$values
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
    statis$C.eig <- eig[1:rank]
    statis$C.nf <- nf
    statis$C.rank <- rank
    wref <- eig1$vectors[, 1:nf]
    wref <- wref/lwsqrt
    w <- data.frame(t(t(wref) * sqrt(eig[1:nf])))
    row.names(w) <- row.names(X)
    names(w) <- paste("C", 1:nf, sep = "")
    statis$C.li <- w
    w <- as.matrix(X[[1]])
    for (k in 2:ntab) {
        w <- cbind(w, as.matrix(X[[k]]))
    }
    w <- w * lw
    w <- t(w) %*% wref
    w <- data.frame(w, row.names = auxinames$col)
    names(w) <- paste("C", 1:nf, sep = "")
    statis$C.Co <- w
    sepanL1 <- sepan(X, nf = 4)$L1
    w <- matrix(0, ntab * 4, nf)
    i1 <- 0
    i2 <- 0
    for (k in 1:ntab) {
        i1 <- i2 + 1
        i2 <- i2 + 4
        tab <- as.matrix(sepanL1[X$TL[, 1] == k, ])
        tab <- t(tab * lw) %*% wref
        for (i in 1:min(nf, 4)) {
            if (tab[i, i] < 0) {
                for (j in 1:nf) tab[i, j] <- -tab[i, j]
            }
        }
        w[i1:i2, ] <- tab
    }
    w <- data.frame(w, row.names = auxinames$tab)
    names(w) <- paste("C", 1:nf, sep = "")
    statis$C.T4 <- w
    w <- as.matrix(statis$C.li) * lwsqrt
    w <- w %*% t(w)
    w <- w/sqrt(sum(w * w))
    w <- as.vector(unlist(w))
    sep <- sep * unlist(w)
    w <- apply(sep, 2, sum)
    statis$cos2 <- w
    statis$tab.names <- tab.names
    statis$TL <- X$TL
    statis$TC <- X$TC
    statis$T4 <- X$T4
    class(statis) <- "statis"
    return(statis)
}

"summary.mcoa" <-
function (object, ...) 
{
    if (!inherits(object, "mcoa")) 
        stop("non convenient data")
    cat("Multiple Co-inertia Analysis\n")
    appel <- as.list(object$call)
    X <- eval(appel$X, sys.frame(0))
    lw <- sqrt(X$lw)
    nlig <- length(lw)
    cw <- X$cw
    ncol <- length(cw)
    blo <- X$blo
    nbloc <- length(X$blo)
    indicablo <- X$TC[, 1]
    nf <- object$nf
    Xsepan <- sepan(X, nf)
    rank.fac <- factor(rep(1:nbloc, Xsepan$rank))
    for (i in 1:nbloc) {
        cat("Array n°", i, names(X)[[i]], "Rows", nrow(X[[i]]), 
            "Cols", ncol(X[[i]]), "\n")
        eigval <- unlist(object$lambda[i, ])
        eigval <- zapsmall(eigval)
        eigvalplus <- zapsmall(cumsum(eigval))
        w <- object$Tli[object$TL[, 1] == i, ]
        w <- w * lw
        varproj <- zapsmall(apply(w * w, 2, sum))
        varprojplus <- zapsmall(cumsum(varproj))
        w1 <- object$SynVar
        w1 <- w1 * lw
        cos2 <- apply(w * w1, 2, sum)
        cos2 <- cos2^2/varproj
        cos2[is.infinite(cos2)] <- NA
        cos2 <- zapsmall(cos2)
        sumry <- array("", c(nf, 6), list(1:nf, c("Iner", "Iner+", 
            "Var", "Var+", "cos2", "cov2")))
        sumry[, 1] <- round(eigval, dig = 3)
        sumry[, 2] <- round(eigvalplus, dig = 3)
        sumry[, 3] <- round(varproj, dig = 3)
        sumry[, 4] <- round(varprojplus, dig = 3)
        sumry[, 5] <- round(cos2, dig = 3)
        sumry[, 6] <- round(object$cov2[i, ], dig = 3)
        class(sumry) <- "table"
        print(sumry)
        cat("\n")
    }
}
"summary.mfa" <-
function (object, ...) 
{
    if (!inherits(object, "mfa")) 
        stop("non convenient data")
    cat("Multiple Factorial Analysis\n")
    cat("rows:", nrow(object$tab), "columns:", ncol(object$tab))
    l0 <- length(object$eig)
    cat("\n\n$eig:", l0, "eigen values\n")
    cat(signif(object$eig, 4)[1:(min(5, l0))])
    if (l0 > 5) 
        cat(" ...\n")
    else cat("\n")
}
"summary.neig" <-
function (object, ...) 
{
    cat("Neigbourhood undirected graph\n")
    deg <- attr(object, "degrees")
    size <- length(deg)
    cat("Vertices:", size, "\n")
    cat("Degrees:", deg, "\n")
    m <- sum(deg)/2
    cat("Edges (pairs of vertices):", m, "\n")
}
"summary.sepan" <-
function (object, ...) 
{
    if (!inherits(object, "sepan")) 
        stop("to be used with 'sepan' object")
    cat("Separate Analyses of a 'ktab' object\n")
    x1 <- object$tab.names
    ntab <- length(x1)
    nf <- min(max(object$rank), 4)
    indica <- factor(rep(1:length(object$blo), object$rank))
    nrow <- nlevels(object$TL[, 2])
    sumry <- array("", c(ntab, 9), list(1:ntab, c("names", "nrow", 
        "ncol", "rank", "lambda1", "lambda2", "lambda3", "lambda4", 
        "")))
    for (k in 1:ntab) {
        eig <- zapsmall(object$Eig[indica == k], dig = 4)
        l0 <- min(length(eig), 4)
        sumry[k, 4 + (1:l0)] <- round(eig[1:l0], dig = 3)
        if (length(eig) > 4) 
            sumry[k, 9] <- "..."
    }
    sumry[, 1] <- x1
    sumry[, 2] <- rep(nrow, ntab)
    sumry[, 3] <- object$blo
    sumry[, 4] <- object$rank
    class(sumry) <- "table"
    print(sumry)
}
"t.ktab" <-
function (x) 
{
    if (!inherits(x, "ktab")) 
        stop("object 'ktab' expected")
    blocks <- x$blo
    nblo <- length(blocks)
    res <- x
    r.n <- row.names(x[[1]])
    for (i in 1:nblo) {
        r.new <- row.names(x[[i]])
        if (any(r.new != r.n)) 
            stop("non equal row.names among array")
    }
    if (length(unique(blocks)) != 1) 
        stop("non equal col numbers among array")
    c.n <- names(x[[1]])
    for (i in 1:nblo) {
        c.new <- names(x[[i]])
        if (any(c.new != c.n)) 
            stop("non equal col.names among array")
    }
    nr <- blocks[1]
    new.row.names <- names(x[[1]])
    indica <- as.factor(rep(1:nblo, blocks))
    w <- split(x$cw, indica)
    col.w <- w[[1]]
    for (i in 1:nblo) {
        col.w.new <- w[[i]]
        if (any(col.w != col.w.new)) 
            stop("non equal column weights among array")
    }
    for (j in 1:nblo) {
        w <- x[[j]]
        w <- data.frame(t(w))
        row.names(w) <- new.row.names
        res[[j]] <- w
        blocks[j] <- ncol(w)
    }
    res$lw <- col.w
    res$cw <- rep(x$lw, nblo)
    res$blo <- blocks
    ktab.util.addfactor(res) <- list(blocks, length(res$lw))
    res$call <- match.call()
    class(res) <- "ktab"
    return(res)
}

"tab.names" <-
function (x) 
UseMethod("tab.names")
"tab.names.ktab" <-
function (x) 
{
    if (!inherits(x, "ktab")) 
        stop("to be used with 'ktab' object")
    ntab <- length(x$blo)
    cha <- names(x)[1:ntab]
    return(cha)
}
"tab.names<-" <-
function (x, value) 
UseMethod("tab.names<-")
"tab.names<-.ktab" <-
function (x, value) 
{
    if (!inherits(x, "ktab")) 
        stop("to be used with 'ktab' object")
    ntab <- length(x$blo)
    old <- tab.names(x)[1:ntab]
    if (!is.null(old) && length(value) != length(old)) 
        stop("invalid tab.names length")
    value <- as.character(value)
    if (any(duplicated(value))) 
        stop("duplicate tab.names are not allowed")
    names(x)[1:ntab] <- value
    x
}


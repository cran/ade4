"dudi.coa" <-
function (df, scannf = TRUE, nf = 2) 
{
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
    df <- df/row.w
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
        call = match.call(), type = "coa")
    X$N <- N
    return(X)
}

"score.coa" <-
function (x, xax = 1, dotchart = FALSE, clab.r = 1, clab.c = 1, 
    csub = 1, cpoi = 1.5, cet = 1.5, ...) 
{
    if (!inherits(x, "coa")) 
        stop("Object of class 'coa' expected")
    if (x$nf == 1) 
        xax <- 1
    if ((xax < 1) || (xax > x$nf)) 
        stop("non convenient axe number")
    "dudi.coa.dotchart" <- function(dudi, numfac, clab) {
        if (!inherits(dudi, "coa")) 
            stop("Object of class 'coa' expected")
        sli <- dudi$li[, numfac]
        sco <- dudi$co[, numfac]
        oli <- order(sli)
        oco <- order(sco)
        a <- c(sli[oli], sco[oco])
        gr <- as.factor(rep(c("Rows", "Columns"), c(length(sli), 
            length(sco))))
        lab <- c(row.names(dudi$li)[oli], row.names(dudi$co)[oco])
        if (clab > 0) 
            labels <- lab
        else labels <- NULL
        dotchart(a, labels = labels, groups = gr, pch = 20)
    }
    if (dotchart) {
        clab <- clab.r * clab.c
        dudi.coa.dotchart(x, xax, clab)
        return(invisible())
    }
    def.par <- par(mar = par("mar"))
    on.exit(par(def.par))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    sco.distri.class.2g <- function(score, fac1, fac2, weight, 
        labels1 = as.character(levels(fac1)), labels2 = as.character(levels(fac2)), 
        clab1, clab2, cpoi, cet) {
        n <- length(score)
        nvar1 <- nlevels(fac1)
        nvar2 <- nlevels(fac2)
        nvar <- nvar1 + nvar2
        ymin <- scoreutil.base(y = score, xlim = NULL, grid = TRUE, 
            cgrid = 0.75, include.origin = TRUE, origin = 0, 
            sub = NULL, csub = 0)
        ymax <- par("usr")[4]
        ylabel <- strheight("A", cex = par("cex") * max(1, clab1, 
            clab2)) * 1.4
        xmin <- par("usr")[1]
        xmax <- par("usr")[2]
        xaxp <- par("xaxp")
        nline <- xaxp[3] + 1
        v0 <- seq(xaxp[1], xaxp[2], le = nline)
        segments(v0, rep(ymin, nline), v0, rep(ymax, nline), 
            col = gray(0.5), lty = 1)
        rect(xmin, ymin, xmax, ymax)
        sum.col1 <- unlist(tapply(weight, fac1, sum))
        sum.col2 <- unlist(tapply(weight, fac2, sum))
        sum.col1[sum.col1 == 0] <- 1
        sum.col2[sum.col2 == 0] <- 1
        weight1 <- weight/sum.col1[fac1]
        weight2 <- weight/sum.col2[fac2]
        y.distri1 <- tapply(score * weight1, fac1, sum)
        y.distri1 <- rank(y.distri1)
        y.distri2 <- tapply(score * weight2, fac2, sum)
        y.distri2 <- rank(y.distri2) + nvar1 + 2
        y.distri <- c(y.distri1, y.distri2)
        ylabel <- strheight("A", cex = par("cex") * max(1, clab1, 
            clab2)) * 1.4
        y.distri1 <- (y.distri1 - min(y.distri))/(max(y.distri) - 
            min(y.distri))
        y.distri1 <- ymin + ylabel + (ymax - ymin - 2 * ylabel) * 
            y.distri1
        y.distri2 <- (y.distri2 - min(y.distri))/(max(y.distri) - 
            min(y.distri))
        y.distri2 <- ymin + ylabel + (ymax - ymin - 2 * ylabel) * 
            y.distri2
        for (i in 1:nvar1) {
            w <- weight1[fac1 == levels(fac1)[i]]
            y0 <- y.distri1[i]
            score0 <- score[fac1 == levels(fac1)[i]]
            x.moy <- sum(w * score0)
            x.et <- sqrt(sum(w * (score0 - x.moy)^2))
            x1 <- x.moy - cet * x.et
            x2 <- x.moy + cet * x.et
            etiagauche <- TRUE
            if ((x1 - xmin) < (xmax - x2)) 
                etiagauche <- FALSE
            segments(x1, y0, x2, y0)
            if (clab1 > 0) {
                cha <- labels1[i]
                cex0 <- par("cex") * clab1
                xh <- strwidth(cha, cex = cex0)
                xh <- xh + strwidth("x", cex = cex0)
                yh <- strheight(cha, cex = cex0) * 5/6
                if (etiagauche) 
                  x0 <- x1 - xh/2
                else x0 <- x2 + xh/2
                rect(x0 - xh/2, y0 - yh, x0 + xh/2, y0 + yh, 
                  col = "white", border = 1)
                text(x0, y0, cha, cex = cex0)
            }
            points(x.moy, y0, pch = 20, cex = par("cex") * cpoi)
        }
        for (i in 1:nvar2) {
            w <- weight2[fac2 == levels(fac2)[i]]
            y0 <- y.distri2[i]
            score0 <- score[fac2 == levels(fac2)[i]]
            x.moy <- sum(w * score0)
            x.et <- sqrt(sum(w * (score0 - x.moy)^2))
            x1 <- x.moy - cet * x.et
            x2 <- x.moy + cet * x.et
            etiagauche <- TRUE
            if ((x1 - xmin) < (xmax - x2)) 
                etiagauche <- FALSE
            segments(x1, y0, x2, y0)
            if (clab2 > 0) {
                cha <- labels2[i]
                cex0 <- par("cex") * clab2
                xh <- strwidth(cha, cex = cex0)
                xh <- xh + strwidth("x", cex = cex0)
                yh <- strheight(cha, cex = cex0) * 5/6
                if (etiagauche) 
                  x0 <- x1 - xh/2
                else x0 <- x2 + xh/2
                rect(x0 - xh/2, y0 - yh, x0 + xh/2, y0 + yh, 
                  col = "white", border = 1)
                text(x0, y0, cha, cex = cex0)
            }
            points(x.moy, y0, pch = 20, cex = par("cex") * cpoi)
        }
    }
    nl <- nrow(x$l1)
    nc <- nrow(x$c1)
    if (inherits(x, "witwit")) {
        y <- eval(as.list(x$call)[[2]], sys.frame(0))
        oritab <- eval(as.list(y$call)[[2]], sys.frame(0))
    }
    else oritab <- eval(as.list(x$call)[[2]], sys.frame(0))
    l.names <- row.names(oritab)
    c.names <- names(oritab)
    oritab <- as.matrix(oritab)
    a <- x$co[col(oritab), xax]
    a <- a + x$li[row(oritab), xax]
    a <- a/sqrt(2 * x$eig[xax] * (1 + sqrt(x$eig[xax])))
    a <- a[oritab > 0]
    aco <- col(oritab)[oritab > 0]
    aco <- factor(aco)
    levels(aco) <- c.names
    ali <- row(oritab)[oritab > 0]
    ali <- factor(ali)
    levels(ali) <- l.names
    aw <- oritab[oritab > 0]/sum(oritab)
    sco.distri.class.2g(a, aco, ali, aw, clab1 = clab.c, clab2 = clab.r, 
        cpoi = cpoi, cet = cet)
    scatterutil.sub("Rows", csub = csub, possub = "topleft")
    scatterutil.sub("Columns", csub = csub, possub = "bottomright")
}

"witwit.coa" <-
function (dudi, row.blocks, col.blocks, scannf = TRUE, nf = 2) 
{
    if (!inherits(dudi, "coa")) 
        stop("Object of class coa expected")
    lig <- nrow(dudi$tab)
    col <- ncol(dudi$tab)
    row.fac <- rep(1:length(row.blocks), row.blocks)
    col.fac <- rep(1:length(col.blocks), col.blocks)
    if (length(col.fac) != col) 
        stop("Non convenient col.fac")
    if (length(row.fac) != lig) 
        stop("Non convenient row.fac")
    tabinit <- as.matrix(eval(as.list(dudi$call)$df, sys.frame(0)))
    tabinit <- tabinit/sum(tabinit)
    wrmat <- rowsum(tabinit, row.fac, reorder = FALSE)[row.fac, 
        ]
    wrvec <- tapply(dudi$lw, row.fac, sum)[row.fac]
    wrvec <- dudi$lw/wrvec
    wrmat <- wrmat * wrvec
    wcmat <- rowsum(t(tabinit), col.fac, reorder = FALSE)[col.fac, 
        ]
    wcvec <- tapply(dudi$cw, col.fac, sum)[col.fac]
    wcvec <- dudi$cw/wcvec
    wcmat <- t(wcmat * wcvec)
    wcmat <- wrmat + wcmat
    wrmat <- rowsum(tabinit, row.fac, reorder = FALSE)
    wrmat <- t(rowsum(t(wrmat), col.fac, reorder = FALSE))
    wrmat <- wrmat[row.fac, col.fac]
    wrmat <- wrmat * wrvec
    wrmat <- t(t(wrmat) * wcvec)
    tabinit <- tabinit - wcmat + wrmat
    tabinit <- tabinit/dudi$lw
    tabinit <- t(t(tabinit)/dudi$cw)
    tabinit <- data.frame(tabinit + wrmat)
    ww <- as.dudi(tabinit, dudi$cw, dudi$lw, scannf = scannf, 
        nf = nf, call = match.call(), type = "witwit")
    class(ww) <- c("witwit", "coa", "dudi")
    wr <- ww$li * ww$li * wrvec
    wr <- rowsum(as.matrix(wr), row.fac, reorder = FALSE)
    cha <- names(row.blocks)
    if (is.null(cha)) 
        cha <- as.character(1:length(row.blocks))
    wr <- data.frame(wr)
    names(wr) <- names(ww$li)
    row.names(wr) <- cha
    ww$lbvar <- wr
    ww$lbw <- tapply(dudi$lw, row.fac, sum)
    wr <- ww$co * ww$co * wcvec
    wr <- rowsum(as.matrix(wr), col.fac, reorder = FALSE)
    cha <- names(col.blocks)
    if (is.null(cha)) 
        cha <- as.character(1:length(col.blocks))
    wr <- data.frame(wr)
    names(wr) <- names(ww$co)
    row.names(wr) <- cha
    ww$cbvar <- wr
    ww$cbw <- tapply(dudi$cw, col.fac, sum)
    return(ww)
}

"score" <-
function (x, ...) 
UseMethod("score")

"sco.boxplot" <-
function (score, df, labels = names(df), clabel = 1, xlim = NULL, 
    grid = TRUE, cgrid = 0.75, include.origin = TRUE, origin = 0, sub = NULL, 
    csub = 1) 
{
    if (!is.vector(score)) 
        stop("vector expected for score")
    if (!is.numeric(score)) 
        stop("numeric expected for score")
    if (!is.data.frame(df)) 
        stop("data.frame expected for df")
    if (!all(unlist(lapply(df, is.factor)))) 
        stop("All variables must be factors")
    n <- length(score)
    if ((nrow(df) != n)) 
        stop("Non convenient match")
    n <- length(score)
    nvar <- ncol(df)
    nlev <- unlist(lapply(df, nlevels))
    opar <- par(mar = par("mar"))
    on.exit(par(opar))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    ymin <- scoreutil.base(y = score, xlim = xlim, grid = grid, 
        cgrid = cgrid, include.origin = include.origin, origin = origin, 
        sub = sub, csub = csub)
    n1 <- sum(nlev)
    ymax <- par("usr")[4]
    ylabel <- strheight("A", cex = par("cex") * max(1, clabel)) * 
        1.4
    yunit <- (ymax - ymin - nvar * ylabel)/n1
    y1 <- ymin + ylabel
    xmin <- par("usr")[1]
    xmax <- par("usr")[2]
    xaxp <- par("xaxp")
    nline <- xaxp[3] + 1
    v0 <- seq(xaxp[1], xaxp[2], le = nline)
    for (i in 1:nvar) {
        y2 <- y1 + nlev[i] * yunit
        rect(xmin, y1, xmax, y2)
        if (clabel > 0) {
            text((xmin + xmax)/2, y1 - ylabel/2, labels[i], cex = par("cex") * 
                clabel)
        }
        param <- tapply(score, df[, i], function(x) quantile(x, 
            seq(0, 1, by = 0.25)))
        moy <- tapply(score, df[, i], mean)
        nbox <- length(param)
        namebox <- names(param)
        pp <- ppoints(n = (nbox + 2), a = 1)
        pp <- pp[2:(nbox + 1)]
        ypp <- y1 + (y2 - y1) * pp
        hbar <- (y2 - y1)/nbox/4
        if (grid) {
            segments(v0, rep(y1, nline), v0, rep(y2, nline), 
                col = gray(0.5), lty = 1)
        }
        for (j in 1:nbox) {
            stat <- unlist(param[j])
            amin <- stat[1]
            aq1 <- stat[2]
            amed <- stat[3]
            aq2 <- stat[4]
            amax <- stat[5]
            rect(aq1, ypp[j] - hbar, aq2, ypp[j] + hbar, col = "white")
            segments(amed, ypp[j] - hbar, amed, ypp[j] + hbar, 
                lwd = 2)
            segments(amin, ypp[j], aq1, ypp[j])
            segments(amax, ypp[j], aq2, ypp[j])
            segments(amin, ypp[j] - hbar, amin, ypp[j] + hbar)
            segments(amax, ypp[j] - hbar, amax, ypp[j] + hbar)
            points(moy[j], ypp[j], pch = 20)
            if (clabel > 0) {
                text(amax, ypp[j], namebox[j], pos = 4, cex = par("cex") * 
                  clabel * 0.8, offset = 0.2)
            }
        }
        y1 <- y2 + ylabel
    }
    invisible()
}
"sco.quant" <-
function (score, df, fac = NULL, clabel = 1, abline = FALSE, 
    sub = names(df), csub = 2, possub = "topleft") 
{
    if (!is.vector(score)) 
        stop("vector expected for score")
    if (!is.numeric(score)) 
        stop("numeric expected for score")
    if (!is.data.frame(df)) 
        stop("data.frame expected for df")
    if (nrow(df) != length(score)) 
        stop("Not convenient dimensions")
    if (!is.null(fac)) {
        fac <- factor(fac)
        if (length(fac) != length(score)) 
            stop("Not convenient dimensions")
    }
    opar <- par(mar = par("mar"), mfrow = par("mfrow"))
    on.exit(par(opar))
    par(mar = c(2.6, 2.6, 1.1, 1.1))
    nfig <- ncol(df)
    par(mfrow = n2mfrow(nfig))
    for (i in 1:nfig) {
        plot(score, df[, i], type = "n")
        if (!is.null(fac)) {
            s.class(cbind.data.frame(score, df[, i]), fac, 
                axesell = FALSE, add.plot = TRUE, clab = clabel)
        }
        else points(score, df[, i])
        if (abline) {
            abline(lm(df[, i] ~ score))
        }
        scatterutil.sub(sub[i], csub, possub)
    }
}
"scoreutil.base" <-
function (y, xlim, grid, cgrid, include.origin, origin, sub, 
    csub) 
{
    if (is.null(xlim)) {
        x1 <- y
        if (include.origin) 
            x1 <- c(x1, origin)
        x1 <- c(x1 - diff(range(x1)/10), x1 + diff(range(x1))/10)
        xlim <- range(x1)
    }
    ylim <- c(0, 1)
    plot.default(0, 0, type = "n", xlab = "", ylab = "", xaxt = "n", 
        yaxt = "n", xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i", 
        frame.plot = FALSE)
    href <- max(3, 2 * cgrid, 2 * csub)
    href <- strheight("A", cex = par("cex") * href)
    if (grid) {
        xaxp <- par("xaxp")
        nline <- xaxp[3] + 1
        v0 <- seq(xaxp[1], xaxp[2], le = nline)
        segments(v0, rep(par("usr")[3], nline), v0, rep(par("usr")[3] + 
            href, nline), col = gray(0.5), lty = 1)
        segments(0, par("usr")[3], 0, par("usr")[3] + href, col = 1, 
            lwd = 3)
        if (cgrid > 0) {
            a <- (xaxp[2] - xaxp[1])/xaxp[3]
            cha <- paste("d = ", a, sep = "")
            cex0 <- par("cex") * cgrid
            xh <- strwidth(cha, cex = cex0)
            yh <- strheight(cha, cex = cex0) + strheight(" ", 
                cex = cex0)/2
            x0 <- strwidth("  ", cex = cex0)
            y0 <- strheight(" ", cex = cex0)/2
            x1 <- par("usr")[1]
            y1 <- par("usr")[3]
            rect(x1 + x0, y1 + y0, x1 + xh + x0, y1 + yh + y0, 
                col = "white", border = 0)
            text(x1 + xh/2 + x0/2, y1 + yh/2 + y0/2, cha, cex = cex0)
        }
    }
    y1 <- rep(par("usr")[3] + href/2, length(y))
    y2 <- rep(par("usr")[3] + href, length(y))
    segments(y, y1, y, y2)
    if (csub > 0) {
        cha <- as.character(sub)
        if (all(c(length(cha) > 0, !is.null(cha), !is.na(cha), 
            cha != ""))) {
            cex0 <- par("cex") * csub
            xh <- strwidth(cha, cex = cex0)
            yh <- strheight(cha, cex = cex0)
            x0 <- strwidth(" ", cex = cex0)
            y0 <- strheight(" ", cex = cex0)
            x1 <- par("usr")[2]
            y1 <- par("usr")[3]
            rect(x1 - x0 - xh, y1, x1, y1 + yh + y0, col = "white", 
                border = 0)
            text(x1 - xh/2 - x0/2, y1 + yh/2 + y0/2, cha, cex = cex0)
        }
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[3] + 
        href)
    return(par("usr")[3] + href)
}
"sco.distri" <-
function (score, df, y.rank = TRUE, csize = 1, labels = names(df), 
    clabel = 1, xlim = NULL, grid = TRUE, cgrid = 0.75, include.origin = TRUE, 
    origin = 0, sub = NULL, csub = 1) 
{
    if (!is.vector(score)) 
        stop("vector expected for score")
    if (!is.numeric(score)) 
        stop("numeric expected for score")
    if (!is.data.frame(df)) 
        stop("data.frame expected for df")
    if (any(df < 0)) 
        stop("data >=0 expected in df")
    n <- length(score)
    if ((nrow(df) != n)) 
        stop("Non convenient match")
    n <- length(score)
    nvar <- ncol(df)
    opar <- par(mar = par("mar"))
    on.exit(par(opar))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    ymin <- scoreutil.base(y = score, xlim = xlim, grid = grid, 
        cgrid = cgrid, include.origin = include.origin, origin = origin, 
        sub = sub, csub = csub)
    ymax <- par("usr")[4]
    ylabel <- strheight("A", cex = par("cex") * max(1, clabel)) * 
        1.4
    xmin <- par("usr")[1]
    xmax <- par("usr")[2]
    xaxp <- par("xaxp")
    nline <- xaxp[3] + 1
    v0 <- seq(xaxp[1], xaxp[2], le = nline)
    if (grid) {
        segments(v0, rep(ymin, nline), v0, rep(ymax, nline), 
            col = gray(0.5), lty = 1)
    }
    rect(xmin, ymin, xmax, ymax)
    sum.col <- apply(df, 2, sum)
    df <- df[, sum.col > 0]
    labels <- labels[sum.col > 0]
    nvar <- ncol(df)
    sum.col <- apply(df, 2, sum)
    df <- sweep(df, 2, sum.col, "/")
    y.distri <- (nvar:1)
    if (y.rank) {
        y.distri <- drop(score %*% as.matrix(df))
        y.distri <- rank(y.distri)
    }
    ylabel <- strheight("A", cex = par("cex") * max(1, clabel)) * 
        1.4
    y.distri <- (y.distri - min(y.distri))/(max(y.distri) - min(y.distri))
    y.distri <- ymin + ylabel + (ymax - ymin - 2 * ylabel) * 
        y.distri
    for (i in 1:nvar) {
        w <- df[, i]
        y0 <- y.distri[i]
        x.moy <- sum(w * score)
        x.et <- sqrt(sum(w * (score - x.moy)^2))
        x1 <- x.moy - x.et * csize
        x2 <- x.moy + x.et * csize
        etiagauche <- TRUE
        if ((x1 - xmin) < (xmax - x2)) 
            etiagauche <- FALSE
        segments(x1, y0, x2, y0)
        if (clabel > 0) {
            cha <- labels[i]
            cex0 <- par("cex") * clabel
            xh <- strwidth(cha, cex = cex0)
            xh <- xh + strwidth("x", cex = cex0)
            yh <- strheight(cha, cex = cex0) * 5/6
            if (etiagauche) 
                x0 <- x1 - xh/2
            else x0 <- x2 + xh/2
            text(x0, y0, cha, cex = cex0)
        }
        points(x.moy, y0, pch = 20, cex = par("cex") * 2)
    }
    invisible()
}

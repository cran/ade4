"add.scatter.eig" <-
function (w, nf, xax, yax, posi = c("bottom", "top", "none"), 
    ratio = 1/4) 
{
    if (posi == "none") 
        return(invisible())
    born <- par("usr")
    opar <- par(mar = par("mar"))
    on.exit(par(opar))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    neig <- min(length(w), 15)
    w <- w[1:neig]
    col <- rep(grey(1), length(w))
    col[1:nf] <- grey(0.8)
    col[c(xax, yax)] <- grey(0)
    x <- seq(born[1], born[1] + (born[2] - born[1]) * ratio, 
        le = neig + 1)
    w <- w/max(w)
    w <- w * (born[4] - born[3]) * ratio
    if (posi == "bottom") 
        m3 <- born[3]
    else m3 <- born[4] - w[1]
    w <- m3 + w
    rect(x[1], m3, x[neig + 1], w[1], col = grey(1))
    for (i in 1:neig) {
        rect(x[i], m3, x[i + 1], w[i], col = col[i])
    }
}
"scatter" <-
function (x, ...) 
UseMethod("scatter")
"scatter.acm" <-
function (x, xax = 1, yax = 2, csub = 2, possub = "topleft", ...) 
{
    if (!inherits(x, "acm")) 
        stop("For 'acm' object")
    if (x$nf == 1) {
        score.(x, 1)
        return(invisible())
    }
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    oritab <- eval(as.list(x$call)[[2]], sys.frame(0))
    nvar <- ncol(oritab)
    par(mfrow = n2mfrow(nvar))
    for (i in 1:(nvar)) s.class(x$li, oritab[, i], clab = 1.5, 
        sub = names(oritab)[i], csub = csub, possub = possub, 
        cgrid = 0, csta = 0)
}
"s.arrow" <-
function (dfxy, xax = 1, yax = 2, label = row.names(dfxy), clabel = 1, 
    pch = 20, cpoint = 0, edge = TRUE, origin = c(0, 0), xlim = NULL, 
    ylim = NULL, grid = TRUE, addaxes = TRUE, cgrid = 1, sub = "", 
    csub = 1.25, possub = "bottomleft", pixmap = NULL, contour = NULL, 
    area = NULL, add.plot = FALSE) 
{
    arrow1 <- function(x0, y0, x1, y1, len = 0.1, ang = 15, lty = 1, 
        edge) {
        d0 <- sqrt((x0 - x1)^2 + (y0 - y1)^2)
        if (d0 < 1e-07) 
            return(invisible())
        segments(x0, y0, x1, y1, lty = lty)
        h <- strheight("A", cex = par("cex"))
        if (d0 > 2 * h) {
            x0 <- x1 - h * (x1 - x0)/d0
            y0 <- y1 - h * (y1 - y0)/d0
            if (edge) 
                arrows(x0, y0, x1, y1, ang = ang, len = len, 
                  lty = 1)
        }
    }
    dfxy <- data.frame(dfxy)
    opar <- par(mar = par("mar"))
    on.exit(par(opar))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax, 
        xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, 
        cgrid = cgrid, include.origin = TRUE, origin = origin, 
        sub = sub, csub = csub, possub = possub, pixmap = pixmap, 
        contour = contour, area = area, add.plot = add.plot)
    if (grid & !add.plot) 
        scatterutil.grid(cgrid)
    if (addaxes & !add.plot) 
        abline(h = 0, v = 0, lty = 1)
    if (cpoint > 0) 
        points(coo$x, coo$y, pch = pch, cex = par("cex") * cpoint)
    for (i in 1:(length(coo$x))) arrow1(origin[1], origin[2], 
        coo$x[i], coo$y[i], edge = edge)
    if (clabel > 0) 
        scatterutil.eti.circ(coo$x, coo$y, label, clabel)
    if (csub > 0) 
        scatterutil.sub(sub, csub, possub)
    box()
}
"s.chull" <-
function (dfxy, fac, xax = 1, yax = 2, optchull = c(0.25, 0.5, 
    0.75, 1), label = levels(fac), clabel = 1, cpoint = 0, xlim = NULL, 
    ylim = NULL, grid = TRUE, addaxes = TRUE, origin = c(0, 0), 
    include.origin = TRUE, sub = "", csub = 1, possub = "bottomleft", 
    cgrid = 1, pixmap = NULL, contour = NULL, area = NULL, add.plot = FALSE) 
{
    dfxy <- data.frame(dfxy)
    opar <- par(mar = par("mar"))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    on.exit(par(opar))
    coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax, 
        xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, 
        cgrid = cgrid, include.origin = include.origin, origin = origin, 
        sub = sub, csub = csub, possub = possub, pixmap = pixmap, 
        contour = contour, area = area, add.plot = add.plot)
    scatterutil.chull(coo$x, coo$y, fac, optchull = optchull)
    if (cpoint > 0) 
        points(coo$x, coo$y, pch = 20, cex = par("cex") * cpoint)
    if (clabel > 0) {
        coox <- tapply(coo$x, fac, mean)
        cooy <- tapply(coo$y, fac, mean)
        scatterutil.eti(coox, cooy, label, clabel)
    }
    box()
}
"s.class" <-
function (dfxy, fac, wt = rep(1, length(fac)), xax = 1, yax = 2, 
    cstar = 1, cellipse = 1.5, axesell = TRUE, label = levels(fac), 
    clabel = 1, cpoint = 1, pch = 20, xlim = NULL, ylim = NULL, 
    grid = TRUE, addaxes = TRUE, origin = c(0, 0), include.origin = TRUE, 
    sub = "", csub = 1, possub = "bottomleft", cgrid = 1, pixmap = NULL, 
    contour = NULL, area = NULL, add.plot = FALSE) 
{
    f1 <- function(cl) {
        n <- length(cl)
        cl <- as.factor(cl)
        x <- matrix(0, n, length(levels(cl)))
        x[(1:n) + n * (unclass(cl) - 1)] <- 1
        dimnames(x) <- list(names(cl), levels(cl))
        data.frame(x)
    }
    opar <- par(mar = par("mar"))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    on.exit(par(opar))
    dfxy <- data.frame(dfxy)
    if (!is.data.frame(dfxy)) 
        stop("Non convenient selection for dfxy")
    if (any(is.na(dfxy))) 
        stop("NA non implemented")
    if (!is.factor(fac)) 
        stop("factor expected for fac")
    dfdistri <- f1(fac) * wt
    w1 <- unlist(lapply(dfdistri, sum))
    dfdistri <- t(t(dfdistri)/w1)
    coox <- as.matrix(t(dfdistri)) %*% dfxy[, xax]
    cooy <- as.matrix(t(dfdistri)) %*% dfxy[, yax]
    if (nrow(dfxy) != nrow(dfdistri)) 
        stop(paste("Non equal row numbers", nrow(dfxy), nrow(dfdistri)))
    coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax, 
        xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, 
        cgrid = cgrid, include.origin = include.origin, origin = origin, 
        sub = sub, csub = csub, possub = possub, pixmap = pixmap, 
        contour = contour, area = area, add.plot = add.plot)
    if (cpoint > 0) 
        points(coo$x, coo$y, pch = pch, cex = par("cex") * cpoint)
    if (cstar > 0) 
        for (i in 1:ncol(dfdistri)) {
            scatterutil.star(coo$x, coo$y, dfdistri[, i], cstar = cstar)
        }
    if (cellipse > 0) 
        for (i in 1:ncol(dfdistri)) {
            scatterutil.ellipse(coo$x, coo$y, dfdistri[, i], 
                cellipse = cellipse, axesell = axesell)
        }
    if (clabel > 0) 
        scatterutil.eti(coox, cooy, label, clabel)
    box()
}
"scatter.coa" <-
function (x, xax = 1, yax = 2, method = 1:3, clab.row = 0.75, 
    clab.col = 1.25, posieig = "top", sub = NULL, csub = 2, ...) 
{
    if (!inherits(x, "dudi")) 
        stop("Object of class 'dudi' expected")
    if (!inherits(x, "coa")) 
        stop("Object of class 'coa' expected")
    nf <- x$nf
    if ((xax > nf) || (xax < 1) || (yax > nf) || (yax < 1) || 
        (xax == yax)) 
        stop("Non convenient selection")
    if (method == 1) {
        coolig <- x$li[, c(xax, yax)]
        coocol <- x$co[, c(xax, yax)]
        names(coocol) <- names(coolig)
        s.label(rbind.data.frame(coolig, coocol), clab = 0, 
            cpoi = 0, sub = sub, csub = csub)
        s.label(coolig, xax, yax, clab = clab.row, add.p = TRUE)
        s.label(coocol, xax, yax, clab = clab.col, add.p = TRUE)
    }
    else if (method == 2) {
        coocol <- x$c1[, c(xax, yax)]
        coolig <- x$li[, c(xax, yax)]
        s.label(coocol, clab = clab.col, sub = sub, csub = csub)
        s.label(coolig, clab = clab.row, add.plot = TRUE)
    }
    else if (method == 3) {
        coolig <- x$l1[, c(xax, yax)]
        coocol <- x$co[, c(xax, yax)]
        s.label(coolig, clab = clab.col, sub = sub, csub = csub)
        s.label(coocol, clab = clab.row, add.plot = TRUE)
    }
    else stop("Unknown method")
    add.scatter.eig(x$eig, x$nf, xax, yax, posi = posieig, ratio = 1/4)
}
"s.corcircle" <-
function (dfxy, xax = 1, yax = 2, label = row.names(df), clabel = 1, 
    grid = TRUE, sub = "", csub = 1, possub = "bottomleft", cgrid = 0, 
    fullcircle = TRUE, box = FALSE, add.plot = FALSE) 
{
    arrow1 <- function(x0, y0, x1, y1, len = 0.1, ang = 15, lty = 1, 
        edge) {
        d0 <- sqrt((x0 - x1)^2 + (y0 - y1)^2)
        if (d0 < 1e-07) 
            return(invisible())
        segments(x0, y0, x1, y1, lty = lty)
        h <- strheight("A", cex = par("cex"))
        if (d0 > 2 * h) {
            x0 <- x1 - h * (x1 - x0)/d0
            y0 <- y1 - h * (y1 - y0)/d0
            if (edge) 
                arrows(x0, y0, x1, y1, ang = ang, len = len, 
                  lty = 1)
        }
    }
    scatterutil.circ <- function(cgrid, h) {
        cc <- seq(from = -1, to = 1, by = h)
        col <- "lightgray"
        lty <- 1
        for (i in 1:(length(cc))) {
            x <- cc[i]
            a1 <- sqrt(1 - x * x)
            a2 <- (-a1)
            segments(x, a1, x, a2, col = col)
            segments(a1, x, a2, x, col = col)
        }
        symbols(0, 0, circles = 1, inches = FALSE, add = TRUE)
        segments(-1, 0, 1, 0)
        segments(0, -1, 0, 1)
        if (cgrid <= 0) 
            return(invisible())
        cha <- paste("d = ", h, sep = "")
        cex0 <- par("cex") * cgrid
        xh <- strwidth(cha, cex = cex0)
        yh <- strheight(cha, cex = cex0) + strheight(" ", cex = cex0)/2
        x0 <- strwidth(" ", cex = cex0)
        y0 <- strheight(" ", cex = cex0)/2
        x1 <- par("usr")[2]
        y1 <- par("usr")[4]
        rect(x1 - x0, y1 - y0, x1 - xh - x0, y1 - yh - y0, col = "white", 
            border = 0)
        text(x1 - xh/2 - x0/2, y1 - yh/2 - y0/2, cha, cex = cex0)
    }
    df <- data.frame(dfxy)
    if (!is.data.frame(df)) 
        stop("Non convenient selection for df")
    if ((xax < 1) || (xax > ncol(df))) 
        stop("Non convenient selection for xax")
    if ((yax < 1) || (yax > ncol(df))) 
        stop("Non convenient selection for yax")
    if (!is.null(label)) 
        showpoint <- FALSE
    x <- df[, xax]
    y <- df[, yax]
    if (add.plot) {
        for (i in 1:length(x)) arrow1(0, 0, x[i], y[i], len = 0.1, 
            ang = 15, edge = TRUE)
        if (clabel > 0) 
            scatterutil.eti.circ(x, y, label, clabel)
        return(invisible())
    }
    opar <- par(mar = par("mar"))
    on.exit(par(opar))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    x1 <- x
    y1 <- y
    x1 <- c(x1, -0.01, +0.01)
    y1 <- c(y1, -0.01, +0.01)
    if (fullcircle) {
        x1 <- c(x1, -1, 1)
        y1 <- c(y1, -1, 1)
    }
    x1 <- c(x1 - diff(range(x1)/20), x1 + diff(range(x1))/20)
    y1 <- c(y1 - diff(range(y1)/20), y1 + diff(range(y1))/20)
    plot(x1, y1, type = "n", ylab = "", asp = 1, xaxt = "n", 
        yaxt = "n", frame.plot = FALSE)
    if (grid) 
        scatterutil.circ(cgrid = cgrid, h = 0.2)
    for (i in 1:length(x)) arrow1(0, 0, x[i], y[i], len = 0.1, 
        ang = 15, edge = TRUE)
    if (clabel > 0) 
        scatterutil.eti.circ(x, y, label, clabel)
    if (csub > 0) 
        scatterutil.sub(sub, csub, possub)
    if (box) 
        box()
}
"s.distri" <-
function (dfxy, dfdistri, xax = 1, yax = 2, cstar = 1, cellipse = 1.5, 
    axesell = TRUE, label = names(dfdistri), clabel = 0, cpoint = 1, 
    pch = 20, xlim = NULL, ylim = NULL, grid = TRUE, addaxes = TRUE, 
    origin = c(0, 0), include.origin = TRUE, sub = "", csub = 1, 
    possub = "bottomleft", cgrid = 1, pixmap = NULL, contour = NULL, 
    area = NULL, add.plot = FALSE) 
{
    opar <- par(mar = par("mar"))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    on.exit(par(opar))
    dfxy <- data.frame(dfxy)
    dfdistri <- data.frame(dfdistri)
    if (!is.data.frame(dfxy)) 
        stop("Non convenient selection for dfxy")
    if (!is.data.frame(dfdistri)) 
        stop("Non convenient selection for dfdistri")
    if (any(dfdistri < 0)) 
        stop("Non convenient selection for dfdistri")
    if (nrow(dfxy) != nrow(dfdistri)) 
        stop("Non equal row numbers")
    if (any(is.na(dfxy))) 
        stop("NA non implemented")
    w1 <- unlist(lapply(dfdistri, sum))
    label <- label
    dfdistri <- t(t(dfdistri)/w1)
    coox <- as.matrix(t(dfdistri)) %*% as.matrix(dfxy[, xax])
    cooy <- as.matrix(t(dfdistri)) %*% as.matrix(dfxy[, yax])
    coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax, 
        xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, 
        cgrid = cgrid, include.origin = include.origin, origin = origin, 
        sub = sub, csub = csub, possub = possub, pixmap = pixmap, 
        contour = contour, area = area, add.plot = add.plot)
    if (cpoint > 0) 
        points(coo$x, coo$y, pch = pch, cex = par("cex") * cpoint)
    if (cstar > 0) 
        for (i in 1:ncol(dfdistri)) {
            scatterutil.star(coo$x, coo$y, dfdistri[, i], cstar = cstar)
        }
    if (cellipse > 0) 
        for (i in 1:ncol(dfdistri)) {
            scatterutil.ellipse(coo$x, coo$y, dfdistri[, i], 
                cellipse = cellipse, axesell = axesell)
        }
    if (clabel > 0) 
        scatterutil.eti(unlist(coox), unlist(cooy), label, clabel)
    box()
}

"s.label" <-
function (dfxy, xax = 1, yax = 2, label = row.names(dfxy), clabel = 1, 
    pch = 20, cpoint = if (clabel == 0) 1 else 0, neig = NULL, 
    cneig = 2, xlim = NULL, ylim = NULL, grid = TRUE, addaxes = TRUE, 
    cgrid = 1, include.origin = TRUE, origin = c(0, 0), sub = "", 
    csub = 1.25, possub = "bottomleft", pixmap = NULL, contour = NULL, 
    area = NULL, add.plot = FALSE) 
{
    dfxy <- data.frame(dfxy)
    opar <- par(mar = par("mar"))
    on.exit(par(opar))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax, 
        xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, 
        cgrid = cgrid, include.origin = include.origin, origin = origin, 
        sub = sub, csub = csub, possub = possub, pixmap = pixmap, 
        contour = contour, area = area, add.plot = add.plot)
    if (!is.null(neig)) {
        if (is.null(class(neig))) 
            neig <- NULL
        if (class(neig) != "neig") 
            neig <- NULL
        deg <- attr(neig, "degrees")
        if ((length(deg)) != (length(coo$x))) 
            neig <- NULL
    }
    if (!is.null(neig)) {
        fun <- function(x, coo) {
            segments(coo$x[x[1]], coo$y[x[1]], coo$x[x[2]], coo$y[x[2]], 
                lwd = par("lwd") * cneig)
        }
        apply(unclass(neig), 1, fun, coo = coo)
    }
    if (clabel > 0) 
        scatterutil.eti(coo$x, coo$y, label, clabel)
    if (cpoint > 0) 
        points(coo$x, coo$y, pch = pch, cex = par("cex") * cpoint)
    box()
}
"s.match" <-
function (df1xy, df2xy, xax = 1, yax = 2, pch = 20, cpoint = 1, 
    label = row.names(df1xy), clabel = 1, edge = TRUE, xlim = NULL, 
    ylim = NULL, grid = TRUE, addaxes = TRUE, cgrid = 1, include.origin = TRUE, 
    origin = c(0, 0), sub = "", csub = 1.25, possub = "bottomleft", 
    pixmap = NULL, contour = NULL, area = NULL, add.plot = FALSE) 
{
    arrow1 <- function(x0, y0, x1, y1, len = 0.1, ang = 15, lty = 1, 
        edge) {
        d0 <- sqrt((x0 - x1)^2 + (y0 - y1)^2)
        if (d0 < 1e-07) 
            return(invisible())
        segments(x0, y0, x1, y1, lty = lty)
        h <- strheight("A", cex = par("cex"))
        if (d0 > 2 * h) {
            x0 <- x1 - h * (x1 - x0)/d0
            y0 <- y1 - h * (y1 - y0)/d0
            if (edge) 
                arrows(x0, y0, x1, y1, ang = ang, len = len, 
                  lty = 1)
        }
    }
    df1xy <- data.frame(df1xy)
    df2xy <- data.frame(df2xy)
    if (!is.data.frame(df1xy)) 
        stop("Non convenient selection for df1xy")
    if (!is.data.frame(df2xy)) 
        stop("Non convenient selection for df2xy")
    if (any(is.na(df1xy))) 
        stop("NA non implemented")
    if (any(is.na(df2xy))) 
        stop("NA non implemented")
    n <- nrow(df1xy)
    if (n != nrow(df2xy)) 
        stop("Non equal row numbers")
    opar <- par(mar = par("mar"))
    on.exit(par(opar))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    coo <- scatterutil.base(dfxy = rbind.data.frame(df1xy, df2xy), 
        xax = xax, yax = yax, xlim = xlim, ylim = ylim, grid = grid, 
        addaxes = addaxes, cgrid = cgrid, include.origin = include.origin, 
        origin = origin, sub = sub, csub = csub, possub = possub, 
        pixmap = pixmap, contour = contour, area = area, add.plot = add.plot)
    for (i in 1:n) {
        arrow1(coo$x[i], coo$y[i], coo$x[i + n], coo$y[i + n], 
            lty = 1, edge = edge)
    }
    if (cpoint > 0) 
        points(coo$x[1:n], coo$y[1:n], pch = pch, cex = par("cex") * 
            cpoint)
    if (clabel > 0) {
        a <- (coo$x[1:n] + coo$x[(n + 1):(2 * n)])/2
        b <- (coo$y[1:n] + coo$y[(n + 1):(2 * n)])/2
        scatterutil.eti(a, b, label, clabel)
    }
    box()
}
"s.traject" <-
function (dfxy, fac = factor(rep(1, nrow(dfxy))), ord = (1:length(fac)), 
    xax = 1, yax = 2, label = levels(fac), clabel = 1, cpoint = 1, 
    pch = 20, xlim = NULL, ylim = NULL, grid = TRUE, addaxes = TRUE, 
    edge = TRUE, origin = c(0, 0), include.origin = TRUE, sub = "", 
    csub = 1, possub = "bottomleft", cgrid = 1, pixmap = NULL, 
    contour = NULL, area = NULL, add.plot = FALSE) 
{
    opar <- par(mar = par("mar"))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    on.exit(par(opar))
    dfxy <- data.frame(dfxy)
    if (!is.data.frame(dfxy)) 
        stop("Non convenient selection for dfxy")
    if (any(is.na(dfxy))) 
        stop("NA non implemented")
    if (!is.factor(fac)) 
        stop("factor expected for fac")
    if (length(fac) != nrow(dfxy)) 
        stop("Non convenient length (fac)")
    if (length(ord) != nrow(dfxy)) 
        stop("Non convenient length (ord)")
    coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax, 
        xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, 
        cgrid = cgrid, include.origin = include.origin, origin = origin, 
        sub = sub, csub = csub, possub = possub, pixmap = pixmap, 
        contour = contour, area = area, add.plot = add.plot)
    arrow1 <- function(x0, y0, x1, y1, len = 0.15, ang = 15, 
        lty = 1, edge) {
        d0 <- sqrt((x0 - x1)^2 + (y0 - y1)^2)
        if (d0 < 1e-07) 
            return(invisible())
        segments(x0, y0, x1, y1, lty = lty)
        h <- strheight("A", cex = par("cex"))
        x0 <- x1 - h * (x1 - x0)/d0
        y0 <- y1 - h * (y1 - y0)/d0
        if (edge) 
            arrows(x0, y0, x1, y1, ang = 15, len = 0.1, lty = 1)
    }
    trajec <- function(X, cpoint, clabel, label) {
        if (nrow(X) == 1) 
            return(as.numeric(X[1, ]))
        x <- X$x
        y <- X$y
        ord <- order(X$ord)
        fac <- as.numeric(X$fac)
        dmax <- 0
        xmax <- 0
        ymax <- 0
        for (i in 1:(length(x) - 1)) {
            x0 <- x[ord[i]]
            y0 <- y[ord[i]]
            x1 <- x[ord[i + 1]]
            y1 <- y[ord[i + 1]]
            arrow1(x0, y0, x1, y1, lty = fac, edge = edge)
            if (cpoint > 0) 
                points(x0, y0, pch = 14 + fac, cex = par("cex") * 
                  cpoint)
            d0 <- sqrt((origin[1] - (x0 + x1)/2)^2 + (origin[2] - 
                (y0 + y1)/2)^2)
            if (d0 > dmax) {
                xmax <- (x0 + x1)/2
                ymax <- (y0 + y1)/2
                dmax <- d0
            }
        }
        if (cpoint > 0) 
            points(x[ord[length(x)]], y[ord[length(x)]], pch = 14 + 
                fac, cex = par("cex") * cpoint)
        return(c(xmax, ymax))
    }
    provi <- cbind.data.frame(x = coo$x, y = coo$y, fac = fac, 
        ord = ord)
    provi <- split(provi, fac)
    w <- lapply(provi, trajec, cpoint = cpoint, clabel = clabel, 
        label = label)
    w <- t(data.frame(w))
    if (clabel > 0) 
        scatterutil.eti(w[, 1], w[, 2], label, clabel)
    box()
}
"scatterutil.base" <-
function (dfxy, xax, yax, xlim, ylim, grid, addaxes, cgrid, include.origin, 
    origin, sub, csub, possub, pixmap, contour, area, add.plot) 
{
    df <- data.frame(dfxy)
    if (!is.data.frame(df)) 
        stop("Non convenient selection for df")
    if ((xax < 1) || (xax > ncol(df))) 
        stop("Non convenient selection for xax")
    if ((yax < 1) || (yax > ncol(df))) 
        stop("Non convenient selection for yax")
    x <- df[, xax]
    y <- df[, yax]
    if (is.null(xlim)) {
        x1 <- x
        if (include.origin) 
            x1 <- c(x1, origin[1])
        x1 <- c(x1 - diff(range(x1)/10), x1 + diff(range(x1))/10)
        xlim <- range(x1)
    }
    if (is.null(ylim)) {
        y1 <- y
        if (include.origin) 
            y1 <- c(y1, origin[2])
        y1 <- c(y1 - diff(range(y1)/10), y1 + diff(range(y1))/10)
        ylim <- range(y1)
    }
    if (!is.null(pixmap)) {
        if (is.null(class(pixmap))) 
            pixmap <- NULL
        if (is.na(charmatch("pixmap", class(pixmap)))) 
            pixmap <- NULL
    }
    if (!is.null(pixmap)) {
        dimobj <- attr(pixmap, "dim")
    }
    if (!is.null(contour)) {
        if (!is.data.frame(contour)) 
            contour <- NULL
        if (ncol(contour) != 4) 
            contour <- NULL
    }
    if (!is.null(area)) {
        if (!is.data.frame(area)) 
            area <- NULL
        if (!is.factor(area[, 1])) 
            area <- NULL
        if (ncol(area) < 3) 
            area <- NULL
    }
    if (add.plot) 
        return(list(x = x, y = y))
    plot.default(0, 0, type = "n", asp = 1, xlab = "", ylab = "", 
        xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim, xaxs = "i", 
        yaxs = "i", frame.plot = FALSE)
    if (!is.null(pixmap)) {
        plot(pixmap, add = TRUE)
    }
    if (!is.null(contour)) {
        apply(contour, 1, function(x) segments(x[1], x[2], x[3], 
            x[4], lwd = 1))
    }
    if (!is.null(area)) {
        nlev <- nlevels(area[, 1])
        x1 <- area[, 2]
        x2 <- area[, 3]
        for (i in 1:nlev) {
            lev <- levels(area[, 1])[i]
            a1 <- x1[area[, 1] == lev]
            a2 <- x2[area[, 1] == lev]
            polygon(a1, a2)
        }
    }
    if (grid) 
        scatterutil.grid(cgrid)
    if (addaxes) 
        abline(h = 0, v = 0, lty = 1)
    if (csub > 0) 
        scatterutil.sub(sub, csub, possub)
    return(list(x = x, y = y))
}
"scatterutil.chull" <-
function (x, y, fac, optchull = c(0.25, 0.5, 0.75, 1)) 
{
    if (!is.factor(fac)) 
        return(invisible())
    if (length(x) != length(fac)) 
        return(invisible())
    if (length(y) != length(fac)) 
        return(invisible())
    for (i in 1:nlevels(fac)) {
        x1 <- x[fac == levels(fac)[i]]
        y1 <- y[fac == levels(fac)[i]]
        long <- length(x1)
        longinit <- long
        cref <- 1
        repeat {
            if (long < 3) 
                break
            if (cref == 0) 
                break
            num <- chull(x1, y1)
            x2 <- x1[num]
            y2 <- y1[num]
            taux <- long/longinit
            if ((taux <= cref) & (cref == 1)) {
                cref <- 0.75
                if (any(optchull == 1)) 
                  polygon(x2, y2, lty = 1)
            }
            if ((taux <= cref) & (cref == 0.75)) {
                if (any(optchull == 0.75)) 
                  polygon(x2, y2, lty = 5)
                cref <- 0.5
            }
            if ((taux <= cref) & (cref == 0.5)) {
                if (any(optchull == 0.5)) 
                  polygon(x2, y2, lty = 3)
                cref <- 0.25
            }
            if ((taux <= cref) & (cref == 0.25)) {
                if (any(optchull == 0.25)) 
                  polygon(x2, y2, lty = 2)
                cref <- 0
            }
            x1 <- x1[-num]
            y1 <- y1[-num]
            long <- length(x1)
        }
    }
}
"scatterutil.eigen" <-
function (w, xmax = length(w), ymax = max(w), wsel = 1, sub = "Eigen values", 
    csub = 2, possub = "topright") 
{
    opar <- par(mar = par("mar"))
    on.exit(par(opar))
    par(mar = c(0.8, 2.8, 0.8, 0.8))
    if (length(w) < xmax) 
        w <- c(w, rep(0, xmax - length(w)))
    col.w <- rep(grey(0.8), length(w))
    col.w[wsel] <- grey(0)
    barplot(w, col = col.w, ylim = c(0, ymax))
    scatterutil.sub(cha = sub, csub = csub, possub = possub)
}
"scatterutil.ellipse" <-
function (x, y, z, cellipse, axesell) 
{
    util.ellipse <- function(mx, my, vx, cxy, vy, coeff) {
        lig <- 100
        epsi <- 1e-05
        x <- 0
        y <- 0
        if (vx < 0) 
            vx <- 0
        if (vy < 0) 
            vy <- 0
        if (vx == 0 && vy == 0) 
            return(NULL)
        delta <- (vx - vy) * (vx - vy) + 4 * cxy * cxy
        delta <- sqrt(delta)
        l1 <- (vx + vy + delta)/2
        l2 <- vx + vy - l1
        if (l1 < 0) 
            l1 <- 0
        if (l2 < 0) 
            l2 <- 0
        l1 <- sqrt(l1)
        l2 <- sqrt(l2)
        test <- 0
        if (vx == 0) {
            a0 <- 0
            b0 <- 1
            test <- 1
        }
        if ((vy == 0) && (test == 0)) {
            a0 <- 1
            b0 <- 0
            test <- 1
        }
        if (((abs(cxy)) < epsi) && (test == 0)) {
            a0 <- 1
            b0 <- 0
            test <- 1
        }
        if (test == 0) {
            a0 <- 1
            b0 <- (l1 * l1 - vx)/cxy
            norm <- sqrt(a0 * a0 + b0 * b0)
            a0 <- a0/norm
            b0 <- b0/norm
        }
        a1 <- 2 * pi/lig
        c11 <- coeff * a0 * l1
        c12 <- (-coeff) * b0 * l2
        c21 <- coeff * b0 * l1
        c22 <- coeff * a0 * l2
        angle <- 0
        for (i in 1:lig) {
            cosinus <- cos(angle)
            sinus <- sin(angle)
            x[i] <- mx + c11 * cosinus + c12 * sinus
            y[i] <- my + c21 * cosinus + c22 * sinus
            angle <- angle + a1
        }
        return(list(x = x, y = y, seg1 = c(mx + c11, my + c21, 
            mx - c11, my - c21), seg2 = c(mx + c12, my + c22, 
            mx - c12, my - c22)))
    }
    z <- z/sum(z)
    m1 <- sum(x * z)
    m2 <- sum(y * z)
    v1 <- sum((x - m1) * (x - m1) * z)
    v2 <- sum((y - m2) * (y - m2) * z)
    cxy <- sum((x - m1) * (y - m2) * z)
    ell <- util.ellipse(m1, m2, v1, cxy, v2, cellipse)
    if (is.null(ell)) 
        return(invisible())
    polygon(ell$x, ell$y)
    if (axesell) 
        segments(ell$seg1[1], ell$seg1[2], ell$seg1[3], ell$seg1[4], 
            lty = 2)
    if (axesell) 
        segments(ell$seg2[1], ell$seg2[2], ell$seg2[3], ell$seg2[4], 
            lty = 2)
}
"scatterutil.eti" <-
function (x, y, label, clabel) 
{
    if (length(label) == 0) 
        return(invisible())
    if (is.null(label)) 
        return(invisible())
    if (label == "") 
        return(invisible())
    for (i in 1:(length(x))) {
        cha <- as.character(label[i])
        cha <- paste(" ", cha, " ", sep = "")
        cex0 <- par("cex") * clabel
        xh <- strwidth(cha, cex = cex0)
        yh <- strheight(cha, cex = cex0) * 5/3
        x1 <- x[i]
        y1 <- y[i]
        rect(x1 - xh/2, y1 - yh/2, x1 + xh/2, y1 + yh/2, col = "white", 
            border = 1)
        text(x1, y1, cha, cex = cex0)
    }
}
"scatterutil.eti.circ" <-
function (x, y, label, clabel) 
{
    if (is.null(label)) 
        return(invisible())
    if (is.na(label)) 
        return(invisible())
    if (label == "") 
        return(invisible())
    for (i in 1:(length(x))) {
        cha <- as.character(label[i])
        cha <- paste(" ", cha, " ", sep = "")
        cex0 <- par("cex") * clabel
        xh <- strwidth(cha, cex = cex0)
        yh <- strheight(cha, cex = cex0) * 5/6
        if ((x[i] > y[i]) & (x[i] > -y[i])) {
            x1 <- x[i] + xh/2
            y1 <- y[i]
        }
        else if ((x[i] > y[i]) & (x[i] <= (-y[i]))) {
            x1 <- x[i]
            y1 <- y[i] - yh
        }
        else if ((x[i] <= y[i]) & (x[i] <= (-y[i]))) {
            x1 <- x[i] - xh/2
            y1 <- y[i]
        }
        else if ((x[i] <= y[i]) & (x[i] > (-y[i]))) {
            x1 <- x[i]
            y1 <- y[i] + yh
        }
        rect(x1 - xh/2, y1 - yh, x1 + xh/2, y1 + yh, col = "white", 
            border = 1)
        text(x1, y1, cha, cex = cex0)
    }
}
"scatterutil.grid" <-
function (cgrid) 
{
    col <- "lightgray"
    lty <- 1
    xaxp <- par("xaxp")
    ax <- (xaxp[2] - xaxp[1])/xaxp[3]
    yaxp <- par("yaxp")
    ay <- (yaxp[2] - yaxp[1])/yaxp[3]
    a <- min(ax, ay)
    v0 <- seq(xaxp[1], xaxp[2], by = a)
    h0 <- seq(yaxp[1], yaxp[2], by = a)
    abline(v = v0, col = col, lty = lty)
    abline(h = h0, col = col, lty = lty)
    if (cgrid <= 0) 
        return(invisible())
    cha <- paste(" d = ", a, " ", sep = "")
    cex0 <- par("cex") * cgrid
    xh <- strwidth(cha, cex = cex0)
    yh <- strheight(cha, cex = cex0) * 5/3
    x1 <- par("usr")[2]
    y1 <- par("usr")[4]
    rect(x1 - xh, y1 - yh, x1 + xh, y1 + yh, col = "white", border = 0)
    text(x1 - xh/2, y1 - yh/2, cha, cex = cex0)
}
"scatterutil.legend.bw.square" <-
function (br0, sq0, sig0, clegend) 
{
    br0 <- round(br0, dig = 6)
    cha <- as.character(br0[1])
    for (i in (2:(length(br0)))) cha <- paste(cha, br0[i], sep = " ")
    cex0 <- par("cex") * clegend
    yh <- max(c(strheight(cha, cex = cex0), sq0))
    h <- strheight(cha, cex = cex0)
    y0 <- par("usr")[3] + yh/2 + h/2
    ltot <- strwidth(cha, cex = cex0) + sum(sq0) + h
    rect(par("usr")[1] + h/4, y0 - yh/2 - h/4, par("usr")[1] + 
        ltot + h/4, y0 + yh/2 + h/4, col = "white")
    x0 <- par("usr")[1] + h/2
    for (i in (1:(length(sq0)))) {
        cha <- br0[i]
        cha <- paste(" ", cha, sep = "")
        xh <- strwidth(cha, cex = cex0)
        text(x0 + xh/2, y0, cha, cex = cex0)
        z0 <- sq0[i]
        x0 <- x0 + xh + z0/2
        if (sig0[i] >= 0) 
            symbols(x0, y0, squares = z0, bg = "black", fg = "white", 
                add = TRUE, inch = FALSE)
        else symbols(x0, y0, squares = z0, bg = "white", fg = "black", 
            add = TRUE, inch = FALSE)
        x0 <- x0 + z0/2
    }
    invisible()
}
"scatterutil.legend.square.grey" <-
function (br0, valgris, h, clegend) 
{
    if (clegend <= 0) 
        return(invisible())
    br0 <- round(br0, dig = 6)
    nborn <- length(br0)
    cex0 <- par("cex") * clegend
    x0 <- par("usr")[1] + h
    x1 <- x0
    for (i in (2:(nborn))) {
        x1 <- x1 + h
        cha <- br0[i]
        cha <- paste(cha, "]", sep = "")
        xh <- strwidth(cha, cex = cex0)
        if (i == (nborn)) 
            break
        x1 <- x1 + xh + h
    }
    yh <- max(strheight(paste(br0), cex = cex0), h)
    y0 <- par("usr")[3] + yh/2 + h/2
    rect(par("usr")[1] + h/4, y0 - yh/2 - h/4, x1 - h/4, y0 + 
        yh/2 + h/4, col = "white")
    x0 <- par("usr")[1] + h
    for (i in (2:(nborn))) {
        symbols(x0, y0, squares = h, bg = gray(valgris[i - 1]), 
            add = TRUE, inch = FALSE)
        x0 <- x0 + h
        cha <- br0[i]
        if (cha < 1e-05) 
            cha <- round(cha, dig = 3)
        cha <- paste(cha, "]", sep = "")
        xh <- strwidth(cha, cex = cex0)
        if (i == (nborn)) 
            break
        text(x0 + xh/2, y0, cha, cex = cex0)
        x0 <- x0 + xh + h
    }
    invisible()
}
"scatterutil.legendgris" <-
function (w, nclasslegend, clegend) 
{
    l0 <- as.integer(nclasslegend)
    if (l0 == 0) 
        return(invisible())
    if (l0 == 1) 
        l0 <- 2
    if (l0 > 10) 
        l0 <- 10
    h0 <- 1/(l0 + 1)
    mid0 <- seq(h0/2, 1 - h0/2, le = l0 + 1)
    qq <- quantile(w, seq(0, 1, le = l0 + 1))
    w0 <- as.numeric(cut(w, br = qq, inc = TRUE))
    w0 <- seq(0, 1, le = l0)[w0]
    opar <- par(new = par("new"), mar = par("mar"), usr = par("usr"))
    on.exit(par(opar))
    par(new = TRUE)
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    plot(0, 0, type = "n", xlab = "", ylab = "", xaxt = "n", 
        yaxt = "n", xlim = c(0, 2), ylim = c(0, 1.5))
    rect(rep(0, l0), seq(h0/2, by = h0, le = l0), rep(h0, l0), 
        seq(3 * h0/2, by = h0, le = l0), col = gray(seq(1, 0, 
            le = l0)))
    text(rep(h0, 9), mid0, as.character(signif(qq, dig = 2)), 
        pos = 4, cex = par("cex") * clegend)
    box(col = "white")
}
"scatterutil.scaling" <-
function (refold, refnew, xyold) 
{
    refold <- as.matrix(data.frame(refold))
    refnew <- as.matrix(data.frame(refnew))
    meanold <- apply(refold, 2, mean)
    meannew <- apply(refnew, 2, mean)
    refold0 <- sweep(refold, 2, meanold)
    refnew0 <- sweep(refnew, 2, meannew)
    sold <- sqrt(sum(refold0^2))
    snew <- sqrt(sum(refnew0^2))
    xyold <- sweep(xyold, 2, meanold)
    xyold <- t(t(xyold)/sold)
    xynew <- t(t(xyold) * snew)
    xynew <- sweep(xynew, 2, meannew, "+")
    xynew <- data.frame(xynew)
    names(xynew) <- names(xyold)
    row.names(xynew) <- row.names(xyold)
    return(xynew)
}
"scatterutil.star" <-
function (x, y, z, cstar) 
{
    z <- z/sum(z)
    x1 <- sum(x * z)
    y1 <- sum(y * z)
    for (i in which(z > 0)) {
        hx <- cstar * (x[i] - x1)
        hy <- cstar * (y[i] - y1)
        segments(x1, y1, x1 + hx, y1 + hy)
    }
}
"scatterutil.sub" <-
function (cha, csub, possub = "bottomleft") 
{
    cha <- as.character(cha)
    if (length(cha) == 0) 
        return(invisible())
    if (is.null(cha)) 
        return(invisible())
    if (is.na(cha)) 
        return(invisible())
    if (cha == "") 
        return(invisible())
    if (csub == 0) 
        return(invisible())
    cex0 <- par("cex") * csub
    cha <- paste(" ", cha, " ", sep = "")
    xh <- strwidth(cha, cex = cex0)
    yh <- strheight(cha, cex = cex0) * 5/3
    if (possub == "bottomleft") {
        x1 <- par("usr")[1]
        y1 <- par("usr")[3]
        rect(x1, y1, x1 + xh, y1 + yh, col = "white", border = 0)
        text(x1 + xh/2, y1 + yh/2, cha, cex = cex0)
    }
    else if (possub == "topleft") {
        x1 <- par("usr")[1]
        y1 <- par("usr")[4]
        rect(x1, y1, x1 + xh, y1 - yh, col = "white", border = 0)
        text(x1 + xh/2, y1 - yh/2, cha, cex = cex0)
    }
    else if (possub == "bottomright") {
        x1 <- par("usr")[2]
        y1 <- par("usr")[3]
        rect(x1, y1, x1 - xh, y1 + yh, col = "white", border = 0)
        text(x1 - xh/2, y1 + yh/2, cha, cex = cex0)
    }
    else if (possub == "topright") {
        x1 <- par("usr")[2]
        y1 <- par("usr")[4]
        rect(x1, y1, x1 - xh, y1 - yh, col = "white", border = 0)
        text(x1 - xh/2, y1 - yh/2, cha, cex = cex0)
    }
}
"s.value" <-
function (dfxy, z, xax = 1, yax = 2, method = c("squaresize", 
    "greylevel", "circle"), csize = 1, cpoint = 0, pch = 20, 
    clegend = 0.75, neig = NULL, cneig = 1, xlim = NULL, ylim = NULL, 
    grid = TRUE, addaxes = TRUE, cgrid = 0.75, include.origin = TRUE, 
    origin = c(0, 0), sub = "", csub = 1, possub = "topleft", 
    pixmap = NULL, contour = NULL, area = NULL, add.plot = FALSE) 
{
    dfxy <- data.frame(dfxy)
    if (length(z) != nrow(dfxy)) 
        stop(paste("Non equal row numbers", nrow(dfxy), length(z)))
    opar <- par(mar = par("mar"))
    on.exit(par(opar))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax, 
        xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes, 
        cgrid = cgrid, include.origin = include.origin, origin = origin, 
        sub = sub, csub = csub, possub = possub, pixmap = pixmap, 
        contour = contour, area = area, add.plot = add.plot)
    if (!is.null(neig)) {
        if (is.null(class(neig))) 
            neig <- NULL
        if (class(neig) != "neig") 
            neig <- NULL
        deg <- attr(neig, "degrees")
        if ((length(deg)) != (length(coo$x))) 
            neig <- NULL
    }
    if (!is.null(neig)) {
        fun <- function(x, coo) {
            segments(coo$x[x[1]], coo$y[x[1]], coo$x[x[2]], coo$y[x[2]], 
                lwd = par("lwd") * cneig)
        }
        apply(unclass(neig), 1, fun, coo = coo)
    }
    if (method == "greylevel") {
        br0 <- pretty(z, 6)
        nborn <- length(br0)
        coeff <- diff(range(coo$x))/15
        numclass <- cut.default(z, br0, include = TRUE, lab = FALSE)
        valgris <- seq(1, 0, le = (nborn - 1))
        h <- csize * coeff
        for (i in 1:(nrow(dfxy))) {
            symbols(coo$x[i], coo$y[i], squares = h, bg = gray(valgris[numclass[i]]), 
                add = TRUE, inch = FALSE)
        }
        scatterutil.legend.square.grey(br0, valgris, h/2, clegend)
        if (cpoint > 0) 
            points(coo$x, coo$y, pch = pch, cex = par("cex") * 
                cpoint)
    }
    else if (method == "squaresize") {
        coeff <- diff(range(coo$x))/15
        sq <- sqrt(abs(z))
        w1 <- max(sq)
        sq <- csize * coeff * sq/w1
        for (i in 1:(nrow(dfxy))) {
            if (sign(z[i]) >= 0) {
                symbols(coo$x[i], coo$y[i], squares = sq[i], 
                  bg = 1, fg = 0, add = TRUE, inch = FALSE)
            }
            else {
                symbols(coo$x[i], coo$y[i], squares = sq[i], 
                  bg = "white", fg = 1, add = TRUE, inch = FALSE)
            }
        }
        br0 <- pretty(z, 4)
        l0 <- length(br0)
        br0 <- (br0[1:(l0 - 1)] + br0[2:l0])/2
        sq0 <- sqrt(abs(br0))
        sq0 <- csize * coeff * sq0/w1
        sig0 <- sign(br0)
        if (clegend > 0) 
            scatterutil.legend.bw.square(br0, sq0, sig0, clegend)
        if (cpoint > 0) 
            points(coo$x, coo$y, pch = pch, cex = par("cex") * 
                cpoint)
    }
    else if (method == "circles") {
        print("not yet implemented")
    }
    box()
}

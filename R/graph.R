######################### TRIANGLE ######################################
"add.position.triangle" <-
function (d) 
{
    opar <- par(new = par("new"), mar = par("mar"))
    on.exit(par(opar))
    par(new = TRUE)
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    w <- matrix(0, 3, 3)
    w[1, 1] <- d$mini[1]
    w[1, 2] <- d$mini[2]
    w[1, 3] <- d$maxi[3]
    w[2, 1] <- d$maxi[1]
    w[2, 2] <- d$mini[2]
    w[2, 3] <- d$mini[3]
    w[3, 1] <- d$mini[1]
    w[3, 2] <- d$maxi[2]
    w[3, 3] <- d$mini[3]
    A <- triangle.posipoint(c(0, 0, 1), c(0, 0, 0), c(1, 1, 1))
    B <- triangle.posipoint(c(1, 0, 0), c(0, 0, 0), c(1, 1, 1))
    C <- triangle.posipoint(c(0, 1, 0), c(0, 0, 0), c(1, 1, 1))
    a <- triangle.posipoint(w[1, ], c(0, 0, 0), c(1, 1, 1))
    b <- triangle.posipoint(w[2, ], c(0, 0, 0), c(1, 1, 1))
    c <- triangle.posipoint(w[3, ], c(0, 0, 0), c(1, 1, 1))
    plot(0, 0, type = "n", xlim = c(-0.71, 4 - 0.71), ylim = c(-4 + 
        0.85, 0.85), xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
        asp = 1, frame.plot = FALSE)
    polygon(c(A[1], B[1], C[1]), c(A[2], B[2], C[2]))
    polygon(c(a[1], b[1], c[1]), c(a[2], b[2], c[2]), col = grey(0.75))
}

"triangle.biplot" <-
function (ta1, ta2, label = as.character(1:nrow(ta1)), draw.line = TRUE, 
    show.position = TRUE, scale = TRUE) 
{
    seg <- function(a, b, col = 1) {
        segments(a[1], a[2], b[1], b[2], col = col)
    }
    nam <- names(ta1)
    ta1 <- t(apply(ta1, 1, function(x) x/sum(x)))
    ta2 <- t(apply(ta2, 1, function(x) x/sum(x)))
    d <- triangle.param(rbind(ta1, ta2), scale = scale)
    opar <- par(mar = par("mar"))
    on.exit(par(opar))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    A <- d$A
    B <- d$B
    C <- d$C
    xy <- d$xy
    mini <- d$mini
    maxi <- d$maxi
    plot(0, 0, type = "n", xlim = c(-0.8, 0.8), ylim = c(-0.6, 
        1), xlab = "", ylab = "", xaxt = "n", yaxt = "n", asp = 1, 
        frame.plot = FALSE)
    seg(A, B)
    seg(B, C)
    seg(C, A)
    text(C[1], C[2], labels = paste(mini[1]), pos = 2)
    text(C[1], C[2], labels = paste(maxi[3]), pos = 4)
    text((A + C)[1]/2, (A + C)[2]/2, labels = nam[1], cex = 1.5, 
        pos = 2)
    text(A[1], A[2], labels = paste(maxi[1]), pos = 2)
    text(A[1], A[2], labels = paste(mini[2]), pos = 1)
    text((A + B)[1]/2, (A + B)[2]/2, labels = nam[2], cex = 1.5, 
        pos = 1)
    text(B[1], B[2], labels = paste(maxi[2]), pos = 1)
    text(B[1], B[2], labels = paste(mini[3]), pos = 4)
    text((B + C)[1]/2, (B + C)[2]/2, labels = nam[3], cex = 1.5, 
        pos = 4)
    if (draw.line) {
        nlg <- 10 * (maxi[1] - mini[1])
        for (i in (1:(nlg - 1))) {
            x1 <- A + (i/nlg) * (B - A)
            x2 <- C + (i/nlg) * (B - C)
            seg(x1, x2, col = "lightgrey")
            x1 <- A + (i/nlg) * (B - A)
            x2 <- A + (i/nlg) * (C - A)
            seg(x1, x2, col = "lightgrey")
            x1 <- C + (i/nlg) * (A - C)
            x2 <- C + (i/nlg) * (B - C)
            seg(x1, x2, col = "lightgrey")
        }
    }
    nl <- nrow(ta1)
    for (i in (1:nl)) {
        arrows(xy[i, 1], xy[i, 2], xy[i + nl, 1], xy[i + nl, 
            2], le = 0.1, ang = 15)
    }
    points(xy[1:nrow(ta1), ])
    text(xy[1:nrow(ta1), ], label, pos = 4)
    if (show.position) 
        add.position.triangle(d)
}
"triangle.param" <-
function (ta, scale = TRUE, min3 = NULL, max3 = NULL) 
{
    if (ncol(ta) != 3) 
        stop("Non convenient data")
    if (min(ta) < 0) 
        stop("Non convenient data")
    if ((!is.null(min3)) & (!is.null(max3))) 
        scale <- TRUE
    cal <- matrix(0, 9, 3)
    tb <- t(apply(ta, 1, function(x) x/sum(x)))
    mini <- apply(tb, 2, min)
    maxi <- apply(tb, 2, max)
    mini <- (floor(mini/0.1))/10
    maxi <- (floor(maxi/0.1) + 1)/10
    if (!is.null(min3)) 
        mini <- min3
    if (!is.null(max3)) 
        maxi <- min3
    ampli <- maxi - mini
    amplim <- max(ampli)
    for (j in 1:3) {
        k <- amplim - ampli[j]
        while (k > 0) {
            if ((k > 0) & (maxi[j] < 1)) {
                maxi[j] <- maxi[j] + 0.1
                k <- k - 1
            }
            if ((k > 0) & (mini[j] > 0)) {
                mini[j] <- mini[j] - 0.1
                k <- k - 1
            }
        }
    }
    cal[1, 1] <- mini[1]
    cal[1, 2] <- mini[2]
    cal[1, 3] <- 1 - cal[1, 1] - cal[1, 2]
    cal[2, 1] <- mini[1]
    cal[2, 2] <- maxi[2]
    cal[2, 3] <- 1 - cal[2, 1] - cal[2, 2]
    cal[3, 1] <- maxi[1]
    cal[3, 2] <- mini[2]
    cal[3, 3] <- 1 - cal[3, 1] - cal[3, 2]
    cal[4, 1] <- mini[1]
    cal[4, 3] <- mini[3]
    cal[4, 2] <- 1 - cal[4, 1] - cal[4, 3]
    cal[5, 1] <- mini[1]
    cal[5, 3] <- maxi[3]
    cal[5, 2] <- 1 - cal[5, 1] - cal[5, 3]
    cal[6, 1] <- maxi[1]
    cal[6, 3] <- mini[3]
    cal[6, 2] <- 1 - cal[6, 1] - cal[6, 3]
    cal[7, 2] <- mini[2]
    cal[7, 3] <- mini[3]
    cal[7, 1] <- 1 - cal[7, 2] - cal[7, 3]
    cal[8, 2] <- mini[2]
    cal[8, 3] <- maxi[3]
    cal[8, 1] <- 1 - cal[8, 2] - cal[8, 3]
    cal[9, 2] <- maxi[2]
    cal[9, 3] <- mini[3]
    cal[9, 1] <- 1 - cal[9, 2] - cal[9, 3]
    mini <- apply(cal, 2, min)
    mini <- round(mini, dig = 4)
    maxi <- apply(cal, 2, max)
    maxi <- round(maxi, dig = 4)
    ampli <- maxi - mini
    if (!scale) {
        mini <- c(0, 0, 0)
        maxi <- c(1, 1, 1)
    }
    A <- c(-1/sqrt(2), -1/sqrt(6))
    B <- c(1/sqrt(2), -1/sqrt(6))
    C <- c(0, 2/sqrt(6))
    xy <- t(apply(tb, 1, FUN = triangle.posipoint, mini = mini, 
        maxi = maxi))
    return(list(A = A, B = B, C = C, xy = xy, mini = mini, maxi = maxi))
}
"triangle.plot" <-
function (ta, label = as.character(1:nrow(ta)), clabel = 0, cpoint = 1, 
    draw.line = TRUE, addaxes = FALSE, addmean = FALSE, labeltriangle = TRUE, 
    sub = "", csub = 0, possub = "topright", show.position = TRUE, 
    scale = TRUE, min3 = NULL, max3 = NULL) 
{
    seg <- function(a, b, col = par("col")) {
        segments(a[1], a[2], b[1], b[2], col = col)
    }
    nam <- names(ta)
    ta <- t(apply(ta, 1, function(x) x/sum(x)))
    d <- triangle.param(ta, scale = scale, min3 = min3, max3 = max3)
    opar <- par(mar = par("mar"))
    on.exit(par(opar))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    A <- d$A
    B <- d$B
    C <- d$C
    xy <- d$xy
    mini <- d$mini
    maxi <- d$maxi
    plot(0, 0, type = "n", xlim = c(-0.8, 0.8), ylim = c(-0.6, 
        1), xlab = "", ylab = "", xaxt = "n", yaxt = "n", asp = 1, 
        frame.plot = FALSE)
    seg(A, B)
    seg(B, C)
    seg(C, A)
    text(C[1], C[2], labels = paste(mini[1]), pos = 2)
    text(C[1], C[2], labels = paste(maxi[3]), pos = 4)
    if (labeltriangle) 
        text((A + C)[1]/2, (A + C)[2]/2, labels = nam[1], cex = 1.5, 
            pos = 2)
    text(A[1], A[2], labels = paste(maxi[1]), pos = 2)
    text(A[1], A[2], labels = paste(mini[2]), pos = 1)
    if (labeltriangle) 
        text((A + B)[1]/2, (A + B)[2]/2, labels = nam[2], cex = 1.5, 
            pos = 1)
    text(B[1], B[2], labels = paste(maxi[2]), pos = 1)
    text(B[1], B[2], labels = paste(mini[3]), pos = 4)
    if (labeltriangle) 
        text((B + C)[1]/2, (B + C)[2]/2, labels = nam[3], cex = 1.5, 
            pos = 4)
    if (draw.line) {
        nlg <- 10 * (maxi[1] - mini[1])
        for (i in 1:(nlg - 1)) {
            x1 <- A + (i/nlg) * (B - A)
            x2 <- C + (i/nlg) * (B - C)
            seg(x1, x2, col = "lightgrey")
            x1 <- A + (i/nlg) * (B - A)
            x2 <- A + (i/nlg) * (C - A)
            seg(x1, x2, col = "lightgrey")
            x1 <- C + (i/nlg) * (A - C)
            x2 <- C + (i/nlg) * (B - C)
            seg(x1, x2, col = "lightgrey")
        }
    }
    if (cpoint > 0) 
        points(xy, pch = 20, cex = par("cex") * cpoint)
    if (clabel > 0) 
        scatterutil.eti(xy[, 1], xy[, 2], label, clabel)
    if (addaxes) {
        pr0 <- dudi.pca(ta, scale = FALSE, scann = FALSE)$c1
        w1 <- triangle.posipoint(apply(ta, 2, mean), mini, maxi)
        points(w1[1], w1[2], pch = 16, cex = 2)
        a1 <- pr0[, 1]
        x1 <- a1[1] * A + a1[2] * B + a1[3] * C
        seg(w1 - x1, w1 + x1)
        a1 <- pr0[, 2]
        x1 <- a1[1] * A + a1[2] * B + a1[3] * C
        seg(w1 - x1, w1 + x1)
    }
    if (addmean) {
        m <- apply(ta, 2, mean)
        w1 <- triangle.posipoint(m, mini, maxi)
        points(w1[1], w1[2], pch = 16, cex = 2)
        w2 <- triangle.posipoint(c(m[1], mini[2], 1 - m[1] - 
            mini[2]), mini, maxi)
        w3 <- triangle.posipoint(c(1 - m[2] - mini[3], m[2], 
            mini[3]), mini, maxi)
        w4 <- triangle.posipoint(c(mini[1], 1 - m[3] - mini[1], 
            m[3]), mini, maxi)
        points(w2[1], w2[2], pch = 20, cex = 2)
        points(w3[1], w3[2], pch = 20, cex = 2)
        points(w4[1], w4[2], pch = 20, cex = 2)
        seg(w1, w2)
        seg(w1, w3)
        seg(w1, w4)
        text(w2[1], w2[2], labels = as.character(round(m[1], 
            dig = 3)), cex = 1.5, pos = 2)
        text(w3[1], w3[2], labels = as.character(round(m[2], 
            dig = 3)), cex = 1.5, pos = 1)
        text(w4[1], w4[2], labels = as.character(round(m[3], 
            dig = 3)), cex = 1.5, pos = 4)
    }
    if (csub > 0) 
        scatterutil.sub(sub, csub, possub)
    if (show.position) 
        add.position.triangle(d)
}
"triangle.posipoint" <-
function (x, mini, maxi) 
{
    x <- (x - mini)/(maxi - mini)
    x <- x/sum(x)
    x1 <- (x[2] - x[1])/sqrt(2)
    y1 <- (2 * x[3] - x[2] - x[1])/sqrt(6)
    return(c(x1, y1))
}
############################### AREA ###################################
"area.plot" <-
function (x, values = NULL, graph = NULL, lwdgraph = 2, nclasslegend = 8, 
    clegend = 0.75, sub = "", csub = 1, possub = "topleft", cpoint = 0, 
    label = NULL, clabel = 0, ...) 
{
    x.area <- x
    opar <- par(mar = par("mar"), new = par("new"))
    on.exit(par(opar))
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    if (!is.factor(x.area[, 1])) 
        stop("Factor expected in x.area[1,]")
    fac <- x.area[, 1]
    lev.poly <- unique(fac)
    nlev <- nlevels(lev.poly)
    label.poly <- as.character(unique(x.area[, 1]))
    x1 <- x.area[, 2]
    x2 <- x.area[, 3]
    r1 <- range(x1)
    r2 <- range(x2)
    plot(r1, r2, type = "n", asp = 1, xlab = "", ylab = "", xaxt = "n", 
        yaxt = "n", frame.plot = FALSE)
    if (!is.null(values)) {
        if (!is.vector(values)) 
            values <- as.vector(values)
        if (length(values) != nlev) 
            values <- rep(values, le = nlev)
        br0 <- pretty(values, 6)
        nborn <- length(br0)
        h <- diff(range(x1))/20
        numclass <- cut.default(values, br0, include = TRUE, 
            lab = FALSE, right = TRUE)
        valgris <- seq(1, 0, le = (nborn - 1))
    }
    if (!is.null(graph)) {
        if (class(graph) != "neig") 
            stop("graph need an object of class 'ng'")
    }
    if (cpoint != 0) 
        points(x1, x2, pch = 20, cex = par("cex") * cpoint)
    for (i in 1:nlev) {
        a1 <- x1[fac == lev.poly[i]]
        a2 <- x2[fac == lev.poly[i]]
        if (!is.null(values)) 
            polygon(a1, a2, col = grey(valgris[numclass[i]]))
        else polygon(a1, a2)
    }
    w <- area.util.xy(x.area)
    if (!is.null(graph)) {
        for (i in 1:nrow(graph)) {
            segments(w$x[graph[i, 1]], w$y[graph[i, 1]], w$x[graph[i, 
                2]], w$y[graph[i, 2]], lwd = lwdgraph)
        }
    }
    if (clabel > 0) {
        if (is.null(label)) 
            label <- row.names(w)
        scatterutil.eti(w$x, w$y, label, clabel = clabel)
    }
    scatterutil.sub(sub, csub, possub)
    if (!is.null(values)) 
        scatterutil.legend.square.grey(br0, valgris, h, clegend)
}
"area.util.contour" <-
function (area) 
{
    poly <- area[, 1]
    x <- area[, 2]
    y <- area[, 3]
    res <- NULL
    f1 <- function(x) {
        if (x[1] > x[3]) {
            s <- x[1]
            x[1] <- x[3]
            x[3] <- s
            s <- x[2]
            x[2] <- x[4]
            x[4] <- s
        }
        if (x[1] == x[3]) {
            if (x[2] > x[4]) {
                s <- x[2]
                x[2] <- x[4]
                x[4] <- s
            }
        }
        return(paste(x[1], x[2], x[3], x[4], sep = "A"))
    }
    for (i in 1:(nlevels(poly))) {
        xx <- x[poly == levels(poly)[i]]
        yy <- y[poly == levels(poly)[i]]
        n0 <- length(xx)
        xx <- c(xx, xx[1])
        yy <- c(yy, yy[1])
        z <- cbind(xx[1:n0], yy[1:n0], xx[2:(n0 + 1)], yy[2:(n0 + 
            1)])
        z <- apply(z, 1, f1)
        res <- c(res, z)
    }
    res <- res[table(res)[res] < 2]
    res <- unlist(lapply(res, function(x) as.numeric(unlist(strsplit(x, 
        "A")))))
    res <- matrix(res, ncol = 4, byr = TRUE)
    res <- data.frame(res)
    names(res) <- c("x1", "y1", "x2", "y2")
    return(res)
}
"area.util.xy" <-
function (area) 
{
    fac <- area[, 1]
    lev.poly <- unique(fac)
    npoly <- length(lev.poly)
    x <- rep(0, npoly)
    y <- rep(0, npoly)
    for (i in 1:npoly) {
        lev <- lev.poly[i]
        a1 <- area[fac == lev, 2]
        a2 <- area[fac == lev, 3]
        x[i] <- mean(a1)
        y[i] <- mean(a2)
    }
    cbind.data.frame(x = x, y = y, row.names = as.character(lev.poly))
}
"area2poly" <-
function (area) 
{
    if (!is.factor(area[, 1])) 
        stop("Factor expected in area[,1]")
    fac <- area[, 1]
    lev.poly <- unique(fac)
    nlev <- nlevels(lev.poly)
    label.poly <- as.character(lev.poly)
    x1 <- area[, 2]
    x2 <- area[, 3]
    res <- list()
    for (i in 1:nlev) {
        a1 <- x1[fac == lev.poly[i]]
        a2 <- x2[fac == lev.poly[i]]
        res <- c(res, list(as.matrix(cbind(a1, a2))))
    }
    r0 <- matrix(0, nlev, 4)
    r0[, 1] <- tapply(x1, fac, min)
    r0[, 2] <- tapply(x2, fac, min)
    r0[, 3] <- tapply(x1, fac, max)
    r0[, 4] <- tapply(x2, fac, max)
    class(res) <- "polylist"
    attr(res, "region.id") <- label.poly
    attr(res, "region.rect") <- r0
    return(res)
}
"poly2area" <-
function (polys) 
{
    if (!inherits(polys, "polylist")) 
        stop("Non convenient data")
    if (!is.null(attr(polys, "region.id"))) 
        reg.names <- attr(polys, "region.id")
    else reg.names <- paste("R", 1:length(polys), sep = "")
    area <- data.frame(polys[[1]])
    area <- cbind(rep(reg.names[1], nrow(area)), area)
    names(area) <- c("reg", "x", "y")
    for (i in 2:length(polys)) {
        provi <- data.frame(polys[[i]])
        provi <- cbind(rep(reg.names[i], nrow(provi)), provi)
        names(provi) <- c("reg", "x", "y")
        area <- rbind.data.frame(area, provi)
    }
    area$reg <- factor(area$reg)
    return(area)
}
################################# TABLE ###############################################
"table.cont" <-
function (df, x = 1:ncol(df), y = 1:nrow(df), row.labels = row.names(df), 
    col.labels = names(df), clabel.row = 1, clabel.col = 1, abmean.x = FALSE, 
    abline.x = FALSE, abmean.y = FALSE, abline.y = FALSE, csize = 1, clegend = 0, 
    grid = TRUE) 
{
    opar <- par(mai = par("mai"), srt = par("srt"))
    on.exit(par(opar))
    if (any(df < 0)) 
        stop("Non negative values expected")
    N <- sum(df)
    df <- df/sum(df)
    table.prepare(x = x, y = y, row.labels = row.labels, col.labels = col.labels, 
        clabel.row = clabel.row, clabel.col = clabel.col, grid = grid, 
        pos = "leftbottom")
    xtot <- x[col(as.matrix(df))]
    ytot <- y[row(as.matrix(df))]
    coeff <- diff(range(x))/15
    z <- unlist(df)
    sq <- sqrt(abs(z))
    w1 <- max(sq)
    sq <- csize * coeff * sq/w1
    for (i in 1:length(z)) symbols(xtot[i], ytot[i], squares = sq[i], 
        bg = "white", fg = 1, add = TRUE, inch = FALSE)
    f1 <- function(x) {
        w1 <- weighted.mean(val, x)
        val <- (val - w1)^2
        w2 <- sqrt(weighted.mean(val, x))
        return(c(w1, w2))
    }
    if (abmean.x) {
        val <- y
        w <- t(apply(df, 2, f1))
        points(x, w[, 1], pch = 20, cex = 2)
        segments(x, w[, 1] - w[, 2], x, w[, 1] + w[, 2])
    }
    if (abmean.y) {
        val <- x
        w <- t(apply(df, 1, f1))
        points(w[, 1], y, pch = 20, cex = 2)
        segments(w[, 1] - w[, 2], y, w[, 1] + w[, 2], y)
    }
    df <- as.matrix(df)
    x <- x[col(df)]
    y <- y[row(df)]
    df <- as.vector(df)
    if (abline.x) {
        abline(lm(y ~ x, wei = df))
    }
    if (abline.y) {
        w <- coefficients(lm(x ~ y, wei = df))
        if (w[2] == 0) 
            abline(h = w[1])
        else abline(c(-w[1]/w[2], 1/w[2]))
    }
    br0 <- pretty(z, 4)
    l0 <- length(br0)
    br0 <- (br0[1:(l0 - 1)] + br0[2:l0])/2
    sq0 <- sqrt(abs(br0))
    sq0 <- csize * coeff * sq0/w1
    sig0 <- sign(br0)
    if (clegend > 0) 
        scatterutil.legend.bw.square(br0, sq0, sig0, clegend)
}
"table.dist" <-
function (d, x = 1:(attr(d, "Size")), labels = as.character(x), 
    clabel = 1, csize = 1, grid = TRUE) 
{
    opar <- par(mai = par("mai"), srt = par("srt"))
    on.exit(par(opar))
    if (!inherits(d, "dist")) 
        stop("object of class 'dist expected")
    table.prepare(x, x, labels, labels, clabel, clabel, grid, 
        "leftbottom")
    n <- attr(d, "Size")
    d <- dist2mat(d)
    xtot <- x[col(d)]
    ytot <- x[row(d)]
    coeff <- diff(range(x))/n
    z <- as.vector(d)
    sq <- sqrt(z * pi)
    w1 <- max(sq)
    sq <- csize * coeff * sq/w1
    symbols(xtot, ytot, circles = sq, fg = 1, bg = grey(0.8), 
        add = TRUE, inch = FALSE)
}
"table.paint" <-
function (df, x = 1:ncol(df), y = nrow(df):1, row.labels = row.names(df), 
    col.labels = names(df), clabel.row = 1, clabel.col = 1, csize = 1, 
    clegend = 1) 
{
    x <- rank(x)
    y <- rank(y)
    opar <- par(mai = par("mai"), srt = par("srt"))
    on.exit(par(opar))
    table.prepare(x = x, y = y, row.labels = row.labels, col.labels = col.labels, 
        clabel.row = clabel.row, clabel.col = clabel.col, grid = FALSE, 
        pos = "paint")
    xtot <- x[col(as.matrix(df))]
    ytot <- y[row(as.matrix(df))]
    xdelta <- (max(x) - min(x))/(length(x) - 1)/2
    ydelta <- (max(y) - min(y))/(length(y) - 1)/2
    coeff <- diff(range(xtot))/15
    z <- unlist(df)
    br0 <- pretty(z, 6)
    nborn <- length(br0)
    coeff <- diff(range(x))/15
    numclass <- cut.default(z, br0, include = TRUE, lab = FALSE)
    valgris <- seq(1, 0, le = (nborn - 1))
    h <- csize * coeff
    rect(xtot - xdelta, ytot - ydelta, xtot + xdelta, ytot + 
        ydelta, col = gray(valgris[numclass]))
    if (clegend > 0) 
        scatterutil.legend.square.grey(br0, valgris, h/2, clegend)
}
"table.prepare" <-
function (x, y, row.labels, col.labels, clabel.row, clabel.col, 
    grid, pos) 
{
    cexrow <- par("cex") * clabel.row
    cexcol <- par("cex") * clabel.col
    wx <- range(x)
    wy <- range(y)
    maxx <- max(x)
    maxy <- max(y)
    minx <- min(x)
    miny <- min(y)
    dx <- diff(wx)/(length(x))
    dy <- diff(wy)/(length(y))
    if (cexrow > 0) {
        ncar <- max(nchar(paste(" ", row.labels, " ", sep = "")))
        strx <- par("cin")[1] * ncar * cexrow/2 + 0.1
    }
    else strx <- 0.1
    if (cexcol > 0) {
        ncar <- max(nchar(paste(" ", col.labels, " ", sep = "")))
        stry <- par("cin")[1] * ncar * cexcol/2 + 0.1
    }
    else stry <- 0.1
    if (pos == "righttop") {
        par(mai = c(0.1, 0.1, stry, strx))
        xlim <- wx + c(-dx, 2 * dx)
        ylim <- wy + c(-2 * dy, 2 * dy)
        plot.default(0, 0, type = "n", xlab = "", ylab = "", 
            xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim, 
            xaxs = "i", yaxs = "i", frame.plot = FALSE)
        if (cexrow > 0) {
            for (i in 1:length(y)) {
                ynew <- seq(miny, maxy, le = length(y))
                ynew <- ynew[rank(y)]
                text(maxx + 2 * dx, ynew[i], row.labels[i], adj = 0, 
                  cex = cexrow, xpd = NA)
                segments(maxx + 2 * dx, ynew[i], maxx + dx, y[i])
            }
        }
        if (cexcol > 0) {
            par(srt = 90)
            for (i in 1:length(x)) {
                xnew <- seq(minx, maxx, le = length(x))
                xnew <- xnew[rank(x)]
                text(xnew[i], maxy + 2 * dy, col.labels[i], adj = 0, 
                  cex = cexcol, xpd = NA)
                segments(xnew[i], maxy + 2 * dy, x[i], maxy + 
                  dy)
            }
            par(srt = 0)
        }
        if (grid) {
            col <- "lightgray"
            for (i in 1:length(y)) segments(maxx + dx, y[i], 
                minx - dx, y[i], col = col)
            for (i in 1:length(x)) segments(x[i], miny - dy, 
                x[i], maxy + dy, col = col)
        }
        rect(minx - dx, miny - dy, maxx + dx, maxy + dy)
        return(invisible())
    }
    if (pos == "phylog") {
        par(mai = c(0.1, 0.1, stry, strx))
        xlim <- wx + c(-dx, 2 * dx)
        ylim <- wy + c(-dy, 2 * dy)
        plot.default(0, 0, type = "n", xlab = "", ylab = "", 
            xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim, 
            xaxs = "i", yaxs = "i", frame.plot = FALSE)
        if (cexrow > 0) {
            for (i in 1:length(y)) {
                ynew <- seq(miny, maxy, le = length(y))
                ynew <- ynew[rank(y)]
                text(maxx + 2 * dx, ynew[i], row.labels[i], adj = 0, 
                  cex = cexrow, xpd = NA)
                segments(maxx + 2 * dx, ynew[i], maxx + dx, y[i])
            }
        }
        if (cexcol > 0) {
            par(srt = 90)
            xnew <- x[2:length(x)]
            x <- xnew
            for (i in 1:length(x)) {
                text(xnew[i], maxy + 2 * dy, col.labels[i], adj = 0, 
                  cex = cexcol, xpd = NA)
                segments(xnew[i], maxy + 2 * dy, x[i], maxy + 
                  dy)
            }
            par(srt = 0)
        }
        minx <- min(x)
        if (grid) {
            col <- "lightgray"
            for (i in 1:length(y)) segments(maxx + dx, y[i], 
                minx - dx, y[i], col = col)
            for (i in 1:length(x)) segments(x[i], miny - dy, 
                x[i], maxy + dy, col = col)
        }
        rect(minx - dx, miny - dy, maxx + dx, maxy + dy)
        rect(-dx, miny - dy, minx - dx, maxy + dy)
        return(c(0, minx - dx))
    }
    if (pos == "leftbottom") {
        par(mai = c(stry, strx, 0.05, 0.05))
        xlim <- wx + c(-2 * dx, dx)
        ylim <- wy + c(-2 * dy, dy)
        plot.default(0, 0, type = "n", xlab = "", ylab = "", 
            xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim, 
            xaxs = "i", yaxs = "i", frame.plot = FALSE)
        if (cexrow > 0) {
            for (i in 1:length(y)) {
                ynew <- seq(miny, maxy, le = length(y))
                ynew <- ynew[rank(y)]
                w9 <- strwidth(row.labels[i], cex = cexrow)
                text(minx - w9 - 2 * dx, ynew[i], row.labels[i], 
                  adj = 0, cex = cexrow, xpd = NA)
                segments(minx - 2 * dx, ynew[i], minx - dx, y[i])
            }
        }
        if (cexcol > 0) {
            par(srt = -90)
            for (i in 1:length(x)) {
                xnew <- seq(minx, maxx, le = length(x))
                xnew <- xnew[rank(x)]
                text(xnew[i], miny - 2 * dy, col.labels[i], adj = 0, 
                  cex = cexcol, xpd = NA)
                segments(xnew[i], miny - 2 * dy, x[i], miny - 
                  dy)
            }
            par(srt = 0)
        }
        if (grid) {
            col <- "lightgray"
            for (i in 1:length(y)) segments(maxx + 2 * dx, y[i], 
                minx - dx, y[i], col = col)
            for (i in 1:length(x)) segments(x[i], miny - 2 * 
                dy, x[i], maxy + dy, col = col)
        }
        rect(minx - dx, miny - dy, maxx + dx, maxy + dy)
        return(invisible())
    }
    if (pos == "paint") {
        par(mai = c(0.2, strx, stry, 0.1))
        xlim <- wx + c(-dx, dx)
        ylim <- wy + c(-dy, dy)
        plot.default(0, 0, type = "n", xlab = "", ylab = "", 
            xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim, 
            xaxs = "i", yaxs = "i", frame.plot = FALSE)
        if (cexrow > 0) {
            ynew <- seq(miny, maxy, le = length(y))
            ynew <- ynew[rank(y)]
            w9 <- strwidth(row.labels, cex = cexrow)
            text(minx - w9 - 3 * dx/4, ynew, row.labels, adj = 0, 
                cex = cexrow, xpd = NA)
        }
        if (cexcol > 0) {
            xnew <- seq(minx, maxx, le = length(x))
            xnew <- xnew[rank(x)]
            par(srt = 90)
            text(xnew, maxy + 3 * dy/4, col.labels, adj = 0, 
                cex = cexcol, xpd = NA)
            par(srt = 0)
        }
        return(invisible())
    }
}
"table.value" <-
function (df, x = 1:ncol(df), y = nrow(df):1, row.labels = row.names(df), 
    col.labels = names(df), clabel.row = 1, clabel.col = 1, csize = 1, 
    clegend = 1, grid = TRUE) 
{
    opar <- par(mai = par("mai"), srt = par("srt"))
    on.exit(par(opar))
    table.prepare(x = x, y = y, row.labels = row.labels, col.labels = col.labels, 
        clabel.row = clabel.row, clabel.col = clabel.col, grid = grid, 
        pos = "righttop")
    xtot <- x[col(as.matrix(df))]
    ytot <- y[row(as.matrix(df))]
    coeff <- diff(range(xtot))/15
    z <- unlist(df)
    sq <- sqrt(abs(z))
    w1 <- max(sq)
    sq <- csize * coeff * sq/w1
    for (i in 1:length(z)) {
        if (sign(z[i]) >= 0) {
            symbols(xtot[i], ytot[i], squares = sq[i], bg = 1, 
                fg = 0, add = TRUE, inch = FALSE)
        }
        else {
            symbols(xtot[i], ytot[i], squares = sq[i], bg = "white", 
                fg = 1, add = TRUE, inch = FALSE)
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
}

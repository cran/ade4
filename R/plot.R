"plot.between" <-
function (x, xax = 1, yax = 2, ...) 
{
    bet <- x
    if (!inherits(bet, "between")) 
        stop("Use only with 'between' objects")
    if ((bet$nf == 1) || (xax == yax)) {
        appel <- as.list(bet$call)
        dudi <- eval(appel$dudi, sys.frame(0))
        fac <- eval(appel$fac, sys.frame(0))
        lig <- nrow(dudi$tab)
        if (length(fac) != lig) 
            stop("Non convenient dimension")
        sco.quant(bet$ls[, 1], dudi$tab, fac = fac)
        return(invisible())
    }
    if (xax > bet$nf) 
        stop("Non convenient xax")
    if (yax > bet$nf) 
        stop("Non convenient yax")
    fac <- eval(as.list(bet$call)$fac, sys.frame(0))
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    nf <- layout(matrix(c(1, 2, 3, 4, 4, 5, 4, 4, 6), 3, 3), 
        respect = TRUE)
    par(mar = c(0.2, 0.2, 0.2, 0.2))
    s.arrow(bet$c1, xax = xax, yax = yax, sub = "Canonical weights", 
        csub = 2, clab = 1.25)
    s.arrow(bet$co, xax = xax, yax = yax, sub = "Variables", 
        csub = 2, cgrid = 0, clab = 1.25)
    scatterutil.eigen(bet$eig, wsel = c(xax, yax))
    s.class(bet$ls, fac, xax = xax, yax = yax, sub = "Scores and classes", 
        csub = 2, clab = 1.25)
    s.corcircle(bet$as, xax = xax, yax = yax, sub = "Inertia axes", 
        csub = 2, cgrid = 0, clab = 1.25)
    s.label(bet$li, xax = xax, yax = yax, sub = "Classes", 
        csub = 2, clab = 1.25)
}
"plot.coinertia" <-
function (x, xax = 1, yax = 2, ...) 
{
    if (!inherits(x, "coinertia")) 
        stop("Use only with 'coinertia' objects")
    if (x$nf == 1) {
        warnings("One axis only : not yet implemented")
        return(invisible())
    }
    if (xax > x$nf) 
        stop("Non convenient xax")
    if (yax > x$nf) 
        stop("Non convenient yax")
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    nf <- layout(matrix(c(1, 2, 3, 4, 4, 5, 4, 4, 6), 3, 3), 
        respect = TRUE)
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    s.corcircle(x$aX, xax, yax, sub = "X axes", csub = 2, 
        clab = 1.25)
    s.corcircle(x$aY, xax, yax, sub = "Y axes", csub = 2, 
        clab = 1.25)
    scatterutil.eigen(x$eig, wsel = c(xax, yax))
    s.match(x$mX, x$mY, xax, yax, clab = 1.5)
    s.arrow(x$l1, xax = xax, yax = yax, sub = "Y Canonical weights", 
        csub = 2, clab = 1.25)
    s.arrow(x$c1, xax = xax, yax = yax, sub = "X Canonical weights", 
        csub = 2, clab = 1.25)
}
"plot.discrimin" <-
function (x, xax = 1, yax = 2, ...) 
{
    if (!inherits(x, "discrimin")) 
        stop("Use only with 'discrimin' objects")
    if ((x$nf == 1) || (xax == yax)) {
        if (inherits(x, "coadisc")) {
            appel <- as.list(x$call)
            df <- eval(appel$df, sys.frame(0))
            fac <- eval(appel$fac, sys.frame(0))
            lig <- nrow(df)
            if (length(fac) != lig) 
                stop("Non convenient dimension")
            lig.w <- apply(df, 1, sum)
            lig.w <- lig.w/sum(lig.w)
            cla.w <- as.vector(tapply(lig.w, fac, sum))
            mean.w <- function(x) {
                z <- x * lig.w
                z <- tapply(z, fac, sum)/cla.w
                return(z)
            }
            w <- apply(df, 2, mean.w)
            w <- data.frame(t(w))
            sco.distri(x$fa[, xax], w, clabel = 1, xlim = NULL, 
                grid = TRUE, cgrid = 1, include.origin = TRUE, origin = 0, 
                sub = NULL, csub = 1)
            return(invisible())
        }
        appel <- as.list(x$call)
        dudi <- eval(appel$dudi, sys.frame(0))
        fac <- eval(appel$fac, sys.frame(0))
        lig <- nrow(dudi$tab)
        if (length(fac) != lig) 
            stop("Non convenient dimension")
        sco.quant(x$li[, 1], dudi$tab, fac = fac)
        return(invisible())
    }
    if (xax > x$nf) 
        stop("Non convenient xax")
    if (yax > x$nf) 
        stop("Non convenient yax")
    fac <- eval(as.list(x$call)$fac, sys.frame(0))
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    nf <- layout(matrix(c(1, 2, 3, 4, 4, 5, 4, 4, 6), 3, 3), 
        respect = TRUE)
    par(mar = c(0.2, 0.2, 0.2, 0.2))
    s.arrow(x$fa, xax = xax, yax = yax, sub = "Canonical weights", 
        csub = 2, clab = 1.25)
    s.corcircle(x$va, xax = xax, yax = yax, sub = "Cos(variates,canonical variates)", 
        csub = 2, cgrid = 0, clab = 1.25)
    scatterutil.eigen(x$eig, wsel = c(xax, yax))
    s.class(x$li, fac, xax = xax, yax = yax, sub = "Scores and classes", 
        csub = 2, clab = 1.5)
    s.corcircle(x$cp, xax = xax, yax = yax, sub = "Cos(components,canonical variates)", 
        csub = 2, cgrid = 0, clab = 1.25)
    s.label(x$gc, xax = xax, yax = yax, sub = "Class scores", 
        csub = 2, clab = 1.25)
}

"plot.niche" <-
function (x, xax = 1, yax = 2, ...) 
{
    if (!inherits(x, "niche")) 
        stop("Use only with 'niche' objects")
    if (x$nf == 1) {
        warnings("One axis only : not yet implemented")
        return(invisible())
    }
    if (xax > x$nf) 
        stop("Non convenient xax")
    if (yax > x$nf) 
        stop("Non convenient yax")
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    nf <- layout(matrix(c(1, 2, 3, 4, 4, 5, 4, 4, 6), 3, 3), 
        respect = TRUE)
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    s.corcircle(x$as, xax, yax, sub = "Axis", csub = 2, 
        clab = 1.25)
    s.arrow(x$c1, xax, yax, sub = "Variables", csub = 2, 
        clab = 1.25)
    scatterutil.eigen(x$eig, wsel = c(xax, yax))
    s.label(x$ls, xax, yax, clab = 0, cpo = 2, sub = "Samples and Species", 
        csub = 2)
    s.label(x$li, xax, yax, clab = 1.5, add.p = TRUE)
    s.label(x$ls, xax, yax, clab = 1.25, sub = "Samples", 
        csub = 2)
    s.distri(x$ls, eval(as.list(x$call)[[3]], sys.frame(0)), 
        cstar = 0, axesell = FALSE, cell = 1, sub = "Niches", csub = 2)
}
"plot.pcaiv" <-
function (x, xax = 1, yax = 2, ...) 
{
    if (!inherits(x, "pcaiv")) 
        stop("Use only with 'pcaiv' objects")
    if (x$nf == 1) {
        warnings("One axis only : not yet implemented")
        return(invisible())
    }
    if (xax > x$nf) 
        stop("Non convenient xax")
    if (yax > x$nf) 
        stop("Non convenient yax")
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    nf <- layout(matrix(c(1, 2, 3, 4, 4, 5, 4, 4, 6), 3, 3), 
        respect = TRUE)
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    s.arrow(x$fa, xax, yax, sub = "Loadings", csub = 2, 
        clab = 1.25)
    s.arrow(na.omit(x$cor), xax = xax, yax = yax, sub = "Correlation", 
        csub = 2, clab = 1.25)
    s.corcircle(x$as, xax, yax, sub = "Inertia axes", csub = 2)
    s.match(x$li, x$ls, xax, yax, clab = 1.5, sub = "Scores and predictions", 
        csub = 2)
    if (inherits(x, "cca")) 
        s.label(x$co, xax, yax, clab = 0, cpoi = 3, add.p = TRUE)
    if (inherits(x, "cca")) 
        s.label(x$co, xax, yax, clab = 1.25, sub = "Species", 
            csub = 2)
    else s.arrow(x$c1, xax = xax, yax = yax, sub = "Variables", 
        csub = 2, clab = 1.25)
    scatterutil.eigen(x$eig, wsel = c(xax, yax))
}
"plot.procuste" <-
function (x, xax = 1, yax = 2, ...) 
{
    if (!inherits(x, "procuste")) 
        stop("Use only with 'procuste' objects")
    if (x$nf == 1) {
        warnings("One axis only : not yet implemented")
        return(invisible())
    }
    if (xax > x$nf) 
        stop("Non convenient xax")
    if (yax > x$nf) 
        stop("Non convenient yax")
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    nf <- layout(matrix(c(1, 2, 3, 4, 4, 5, 4, 4, 6), 3, 3), 
        respect = TRUE)
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    s.arrow(x$load1, xax, yax, sub = "Loadings 1", csub = 2, 
        clab = 1.25)
    s.arrow(x$load2, xax, yax, sub = "Loadings 2", csub = 2, 
        clab = 1.25)
    scatterutil.eigen(x$d^2, wsel = c(xax, yax))
    s.match(x$scor1, x$scor2, xax, yax, clab = 1.5, sub = "Common projection", 
        csub = 2)
    s.label(x$scor1, xax = xax, yax = yax, sub = "Array 1", 
        csub = 2, clab = 1.25)
    s.label(x$scor2, xax = xax, yax = yax, sub = "Array 2", 
        csub = 2, clab = 1.25)
}

"plot.within" <-
function (x, xax = 1, yax = 2, ...) 
{
    if (!inherits(x, "within")) 
        stop("Use only with 'within' objects")
    if ((x$nf == 1) || (xax == yax)) {
        return(invisible())
    }
    if (xax > x$nf) 
        stop("Non convenient xax")
    if (yax > x$nf) 
        stop("Non convenient yax")
    fac <- x$fac
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    nf <- layout(matrix(c(1, 2, 3, 4, 4, 5, 4, 4, 6), 3, 3), 
        respect = TRUE)
    par(mar = c(0.2, 0.2, 0.2, 0.2))
    s.arrow(x$c1, xax = xax, yax = yax, sub = "Canonical weights", 
        csub = 2, clab = 1.25)
    s.arrow(x$co, xax = xax, yax = yax, sub = "Variables", 
        csub = 2, clab = 1.25)
    scatterutil.eigen(x$eig, wsel = c(xax, yax))
    s.class(x$ls, fac, xax = xax, yax = yax, sub = "Scores and classes", 
        csub = 2, clab = 1.5, cpoi = 2)
    s.corcircle(x$as, xax = xax, yax = yax, sub = "Inertia axes", 
        csub = 2, cgrid = 0, clab = 1.25)
    s.class(x$li, fac, xax = xax, yax = yax, axesell = FALSE, 
        clab = 0, cstar = 0, sub = "Common centring", csub = 2)
}

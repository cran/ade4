"randtest" <-
function (xtest, ...) 
{
    UseMethod("randtest")
}

"as.randtest" <-
function (sim, obs, call = match.call()) 
{
    res <- list(sim = sim, obs = obs)
    res$rep <- length(sim)
    res$pvalue <- (sum(sim >= obs) + 1)/(length(sim) + 1)
    res$call <- call
    class(res) <- "randtest"
    return(res)
}

"print.randtest" <-
function (x, ...) 
{
    if (!inherits(x, "randtest")) 
        stop("Non convenient data")
    cat("Monte-Carlo test\n")
    cat("Observation:", x$obs, "\n")
    cat("Call: ")
    print(x$call)
    cat("Based on", x$rep, "replicates\n")
    cat("Simulated p-value:", x$pvalue, "\n")
}

"plot.randtest" <-
function (x, nclass = 10, coeff = 1, ...) 
{
    if (!inherits(x, "randtest")) 
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

"randtest.between" <-
function(xtest, nrepet=999, ...) {
	nrepet<-nrepet+1
	if (!inherits(xtest,"dudi"))
		stop("Object of class dudi expected")
	if (!inherits(xtest,"between"))
		stop ("Type 'between' expected")
	appel<-as.list(xtest$call)
	dudi1<-eval(appel$dudi,sys.frame(0))
	fac<-eval(appel$fac,sys.frame(0))
	X<-dudi1$tab
	X.lw<-dudi1$lw
	X.lw<-X.lw/sum(X.lw)
	X.cw<-sqrt(dudi1$cw)
	X<-t(t(X)*X.cw)
	inertot<-sum(dudi1$eig)
#	isim<-testinter(nrepet, X.lw, X.cw, length(unique(fac)), fac, X, nrow(X), ncol(X))/inertot
	isim<-testinter(nrepet, dudi1$lw, dudi1$cw, length(unique(fac)), fac, dudi1$tab, nrow(X), ncol(X))/inertot
	obs<-isim[1]
	return(as.randtest(isim[-1],obs))
}

testinter <- function(npermut, pl, pc, moda, indica, tab, l1, c1)
	.C("testinter",
		as.integer(npermut),
		as.double(pl),
		as.integer(length(pl)),
		as.double(pc),
		as.integer(length(pc)),
		as.integer(moda),
		as.double(indica),
		as.integer(length(indica)),
		as.double(t(tab)),
		as.integer(l1),
		as.integer(c1),
		inersim = double(npermut))$inersim

"procuste.randtest" <-
function(df1, df2, nrepet=999) {
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
	lig<-nrow(X)
	c1<-ncol(X)
	c2<-ncol(Y)
	isim<-testprocuste(nrepet, lig, c1, c2, as.matrix(X), as.matrix(Y))
	obs<-isim[1]
	return(as.randtest(isim[-1],obs))
}

testprocuste <- function(npermut, lig, c1, c2, tab1, tab2)
	.C("testprocuste",
		as.integer(npermut),
		as.integer(lig),
		as.integer(c1),
		as.integer(c2),
		as.double(t(tab1)),
		as.double(t(tab2)),
		inersim = double(npermut))$inersim

"mantel.randtest" <-
function(m1, m2, nrepet=999) {
    if (!inherits(m1, "dist")) 
        stop("Object of class 'dist' expected")
    if (!inherits(m2, "dist")) 
        stop("Object of class 'dist' expected")
    n <- attr(m1, "Size")
    if (n != attr(m2, "Size")) 
        stop("Non convenient dimension")
    m1 <- dist2mat(m1)
    m2 <- dist2mat(m2)
    col <- ncol(m1)
	isim<-testmantel(nrepet, col, as.matrix(m1), as.matrix(m2))
	obs<-isim[1]
	return(as.randtest(isim[-1],obs))
}

testmantel <- function(npermut, col, tab1, tab2)
	.C("testmantel",
		as.integer(npermut),
		as.integer(col),
		as.double(t(tab1)),
		as.double(t(tab2)),
		inersim = double(npermut))$inersim

"randtest.discrimin" <-
function(xtest, nrepet=999, ...) {
	if (!inherits(xtest, "discrimin"))
		stop("'discrimin' object expected")
	appel<-as.list(xtest$call)
	dudi<-eval(appel$dudi,sys.frame(0))
	fac<-eval(appel$fac,sys.frame(0))
	lig<-nrow(dudi$tab)
	if (length(fac)!=lig) stop ("Non convenient dimension")
	rank<-dudi$rank
	dudi<-redo.dudi(dudi,rank)
	# dudi.lw<-dudi$lw
	# dudi<-dudi$l1
	X<-dudi$l1
	X.lw<-dudi$lw
	# dudi et dudi.lw sont soumis a la permutation
	# fac reste fixe

	isim<-testdiscrimin(nrepet, rank, X.lw, length(unique(fac)), fac, X, nrow(X), ncol(X))
	obs<-isim[1]
	return(as.randtest(isim[-1],obs))
}

testdiscrimin <- function(npermut, rank, pl, moda, indica, tab, l1, c1)
	.C("testdiscrimin",
		as.integer(npermut),
		as.double(rank),
		as.double(pl),
		as.integer(length(pl)),
		as.integer(moda),
		as.double(indica),
		as.integer(length(indica)),
		as.double(t(tab)),
		as.integer(l1),
		as.integer(c1),
		inersim = double(npermut))$inersim

"randtest.coinertia" <-
function(xtest, nrepet=999, fixed=0, ...) {
	nrepet<-nrepet+1
	if (!inherits(xtest,"dudi"))
		stop("Object of class dudi expected")
	if (!inherits(xtest,"coinertia"))
		stop("Object of class 'coinertia' expected")
	appel<-as.list(xtest$call)
	dudiX<-eval(appel$dudiX,sys.frame(0))
	dudiY<-eval(appel$dudiY,sys.frame(0))
	X<-dudiX$tab
	X.cw<-dudiX$cw
	X.lw<-dudiX$lw
	appelX<-as.list(dudiX$call)
	apx<-appelX$df
	Xinit<-eval(appelX$df,sys.frame(0))
	if (appelX[[1]] == "dudi.pca") {		
		if (is.null(appelX$scale)) appelX$scale<-TRUE
		if (appelX$scale=="TRUE") appelX$scale<-TRUE
		if (appelX$scale=="FALSE") appelX$scale<-FALSE
		if (is.null(appelX$center)) appelX$center<-TRUE
		if (appelX$center=="TRUE") appelX$center<-TRUE
		if (appelX$center=="FALSE") appelX$center<-FALSE
		if (appelX$center == FALSE && appelX$scale == FALSE) typX<-"nc"
		if (appelX$center == FALSE && appelX$scale == TRUE) typX<-"cs"
		if (appelX$center == TRUE  && appelX$scale == FALSE) typX<-"cp"
		if (appelX$center == TRUE  && appelX$scale == TRUE) typX<-"cn"
	} else if (appelX[[1]] == "dudi.coa") {
		typX<-"fc"
	} else if (appelX[[1]] == "dudi.fca") {
		typX<-"fc"
	} else if (appelX[[1]] == "dudi.mca") {
		typX<-"cm"
	}
	Y<-dudiY$tab
	Y.cw<-dudiY$cw
	Y.lw<-dudiY$lw
	appelY<-as.list(dudiY$call)
	apy<-appelY$df
	Yinit<-eval(appelY$df,sys.frame(0))
	if (appelY[[1]] == "dudi.pca") {		
		if (is.null(appelY$scale)) appelY$scale<-TRUE
		if (appelY$scale=="TRUE") appelY$scale<-TRUE
		if (appelY$scale=="FALSE") appelY$scale<-FALSE
		if (is.null(appelY$center)) appelY$center<-TRUE
		if (appelY$center=="TRUE") appelY$center<-TRUE
		if (appelY$center=="FALSE") appelY$center<-FALSE
		if (appelY$center == FALSE && appelY$scale == FALSE) typY<-"nc"
		if (appelY$center == FALSE && appelY$scale == TRUE) typY<-"cs"
		if (appelY$center == TRUE  && appelY$scale == FALSE) typY<-"cp"
		if (appelY$center == TRUE  && appelY$scale == TRUE) typY<-"cn"
	} else if (appelY[[1]] == "dudi.coa") {
		typY<-"fc"
	} else if (appelY[[1]] == "dudi.fca") {
		typY<-"fc"
	} else if (appelY[[1]] == "dudi.mca") {
		typY<-"cm"
	}
	if (all(X.lw==Y.lw)) {
		if ( all(X.lw==rep(1/nrow(X), nrow(X))) ) {
			isim<-testertrace(nrepet, X.cw, Y.cw, X, Y, nrow(X), ncol(X), ncol(Y))
		} else {
			if (fixed==0) {
				cat("Warning: non uniform weight. The results from simulations\n")
				cat("are not valid if weights are computed from analysed data.\n")
				isim<-testertracenu(nrepet, X.cw, Y.cw, X.lw, X, Y, nrow(X), ncol(X), ncol(Y), Xinit, Yinit, typX, typY)
			} else if (fixed==1) {
				cat("Warning: non uniform weight. The results from permutations\n")
				cat("are valid only if the row weights come from the fixed table.\n")
				cat("The fixed table is table X : ")
				print(apx)
				isim<-testertracenubis(nrepet, X.cw, Y.cw, X.lw, X, Y, nrow(X), ncol(X), ncol(Y), Xinit, Yinit, typX, typY, fixed)
			} else if (fixed==2) {
				cat("Warning: non uniform weight. The results from permutations\n")
				cat("are valid only if the row weights come from the fixed table.\n")
				cat("The fixed table is table Y : ")
				print(apy)
				isim<-testertracenubis(nrepet, X.cw, Y.cw, X.lw, X, Y, nrow(X), ncol(X), ncol(Y), Xinit, Yinit, typX, typY, fixed)
			} else if (fixed==3) stop ("Error : fixed must be =< 2")
		}
		# On calcule le RV a partir de la coinertie
		isim<-isim/sqrt(sum(dudiX$eig^2))/sqrt(sum(dudiY$eig^2))
		obs<-isim[1]
		return(as.randtest(isim[-1],obs))
	} else {
		stop ("Equal row weights expected")
	}
}

testertrace <- function(npermut, pc1, pc2, tab1, tab2, l1, c1, c2)
	.C("testertrace",
		as.integer(npermut),
		as.double(pc1),
		as.integer(length(pc1)),
		as.double(pc2),
		as.integer(length(pc2)),
		as.double(t(tab1)),
		as.integer(l1),
		as.integer(c1),
		as.double(t(tab2)),
		as.integer(l1),
		as.integer(c2),
		inersim = double(npermut))$inersim

testertracenu <- function(npermut, pc1, pc2, pl, tab1, tab2, l1, c1, c2, Xinit, Yinit, typX, typY)
	.C("testertracenu",
		as.integer(npermut),
		as.double(pc1),
		as.integer(length(pc1)),
		as.double(pc2),
		as.integer(length(pc2)),
		as.double(pl),
		as.integer(length(pl)),
		as.double(t(tab1)),
		as.integer(l1),
		as.integer(c1),
		as.double(t(tab2)),
		as.integer(l1),
		as.integer(c2),
		as.double(t(Xinit)),
		as.double(t(Yinit)),
		as.character(typX),
		as.character(typY),
		inersim = double(npermut))$inersim

testertracenubis <- function(npermut, pc1, pc2, pl, tab1, tab2, l1, c1, c2, Xinit, Yinit, typX, typY, fixed)
	.C("testertracenubis",
		as.integer(npermut),
		as.double(pc1),
		as.integer(length(pc1)),
		as.double(pc2),
		as.integer(length(pc2)),
		as.double(pl),
		as.integer(length(pl)),
		as.double(t(tab1)),
		as.integer(l1),
		as.integer(c1),
		as.double(t(tab2)),
		as.integer(l1),
		as.integer(c2),
		as.double(t(Xinit)),
		as.double(t(Yinit)),
		as.character(typX),
		as.character(typY),
		as.integer(fixed),
		inersim = double(npermut))$inersim


.First.lib <- function(lib, pkg) {
  library.dynam("ade4", pkg, lib)
}

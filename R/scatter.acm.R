"scatter.acm" <- function (x, xax = 1, yax = 2, csub = 2, possub = "topleft", ...) {
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
    # modif lundi, d�cembre 16, 2002 at 16:48 
    # suite � message d'Alain Guerreau  
    for (i in 1:(nvar)) s.class(x$li, oritab[, i], xax=xax, yax=yax, clab = 1.5, 
        sub = names(oritab)[i], csub = csub, possub = possub, 
        cgrid = 0, csta = 0)
}

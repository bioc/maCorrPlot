CorrSample = function(x, np, seed, rp, ndx)
#
# Name: CorrSample
# Desc: samples np random pairs of genes from the expression matrix x and 
#       computes their means, variances, and correlation
# Auth: Alexander.Ploner@ki.se 200206
#
# Chng:
#
{
    ## Find the random pairs
    if (!missing(rp)) {
        if (!missing(np) | !missing(seed)) {
            warning("rp specified, seed and np are ignored")
        }
    } else {
        if (!missing(seed)) {
            set.seed(seed)
        } else {
            seed = NULL
        }
        rp = RandPairs(1:nrow(x), np)
    }
    
    ## We delete unspecified genes
    if (!missing(ndx)) {
        x[!ndx] = NA
    }
    
    ## Worker function - not very efficient
    func = function(ndx) {
            g1 = x[ndx[1],]
            g2 = x[ndx[2],]
            c(cor(g1, g2, use="pairwise"),sd(g1, na.rm=TRUE), 
              sd(g2, na.rm=TRUE), mean(g1, na.rm=TRUE), mean(g2, na.rm=TRUE))
    }
    
    ret = t(apply(rp, 1, func))
    ret = data.frame(ret[,1], ret[,2]*ret[,3], rowMeans(ret[,4:5], na.rm=TRUE), 
                     ret[,2:5], rp)
    colnames(ret) = c("Correlation","StdDev","Mean","sd1","sd2","m1","m2",
                      "ndx1","ndx2")
    class(ret) = c("corr.sample","data.frame")
    ret
}

RandPairs = function(probes, number)
#
# Name: RandPairs
# Desc: given a vector of either integers or characters, selects random pairs
# Auth: Alexander.Ploner@ki.se 200206
#
# Chng:
#
{
    ret = matrix(probes[1], nrow=number, ncol=2)
    for (i in 1:number) {
        ret[i,] = sample(probes, 2)
    }
    ret
}    

plot.corr.sample = function(x, ..., cond, groups, grid=TRUE, refline=TRUE, 
                            xlog=TRUE, scatter=FALSE, curve=FALSE, ci=TRUE, 
                            nint=10, alpha=0.95, length=0.1, xlab="Standard Deviation")
#
# Name: plot.corr.sample
# Desc: plotting method for output from CorrSample
# Auth: Alexander.Ploner@ki.se 200206
#
# Chng: 270206 added ylim, xlab, ylab
#       060306 multiple groups and conds use the same set of factor levels
#               across corr.sample objects - better auto.key
#              check for cond that averages across objects   
#
{
    require(lattice)
    
    # Separate correlations and plotting stuff
    args = list(x, ...)
    ndx = sapply(args, inherits, "corr.sample")
    corr.arg = args[ndx]  # Not empty
    latt.arg = args[!ndx] # May be empty
    
    # Set up the common arguments
    common.arg = list(grid=grid, refline=refline, scatter=scatter, curve=curve,
                      ci=ci, nint=nint, length=length, xlab=xlab, alpha=alpha)
    common.arg = c(common.arg, latt.arg)
    # Translate the convenience parameters
    common.arg$scales$x$log = xlog
    
    # One or several correlations?
    nc = length(corr.arg) 
    if (nc == 1) {
        args = list(formula=Correlation~StdDev, data=corr.arg[[1]])
        if (missing(groups)) {
            args$panel = panel.corr.sample
        } else {
            if (length(groups)!=nrow(x)) {
                stop("length group vector does not match x")
            }
            args$panel = panel.superpose
            args$panel.groups = panel.corr.sample
            args$groups = groups
        }
    } else {
        x = data.frame(Correlation =  unlist(lapply(corr.arg, function(x) x$Corr)),
                       StdDev = unlist(lapply(corr.arg, function(x) x$StdDev)))
        if (missing(cond)) {
            cond = factor(1:nc)
        }
        if (!is.list(cond)) {
            cond = list(cond)
        }
        ncond = length(cond)
        if (ncond > 2) {
            warning("only the first two conditions are processed")
            cond = cond[1:2]
            ncond = 2
        }
        formula.char = "Correlation~StdDev"
        formula.sep  = c("|","*")
        cond.check   = NULL
        for (i in 1:ncond) {
            cc = cond[[i]]
            if (length(cc)!=nc) {
                stop("length of condition",i, "must match number of arguments")
            }
            cc = rep(cc, sapply(corr.arg, nrow))            
            x = cbind(x, cc)
            namu = paste("cond",i,sep="")
            colnames(x)[ncol(x)] = namu
            formula.char = paste(formula.char, namu, sep=formula.sep[i])
            cond.check   = paste(cond.check, cond[[i]], sep="/")
        }
        # Now check whether we have (combined) as many conditions as objects
        if (length(unique(cond.check))!=nc)
            warning("Correlations are averaged across objects - check cond")
        args = list(formula=as.formula(formula.char), data=x)
        if (missing(groups)) {
            args$panel = panel.corr.sample
        } else { # we require a list of group vectors (eg a data.frame)
            if (length(groups)!=nc) {
                stop("length groups must match number of data sets")
            }
            for (i in 1:nc) {
                if (length(groups[[i]])!=nrow(corr.arg[[i]])) {
                    stop("length group vector",i," does not match x")
                }
            }
            # Build the combined vector of factors
            groups = lapply(groups, factor)
            # Try to keep some order among levels
            levls  = NULL
            for (i in 1:nc)
                levls = union(levls, levels(groups[[i]]))
            groups = unlist(lapply(groups, as.character))
            groups = factor(groups, levels=levls)
            args$groups = groups
            args$panel = panel.superpose
            args$panel.groups = panel.corr.sample
        }        
    }
    
    # Combine & call; if we're really smooth, we check for overlaps
    args = c(args, common.arg)
    ret = do.call("xyplot", args)
    ret
}

panel.corr.sample = function(x,y, grid=TRUE, refline=TRUE, xlog=TRUE, scatter=FALSE, 
                             curve=FALSE, ci=TRUE, nint=10, alpha=0.95, 
                             length=0.1, col.line, col.symbol, ...)
#
# Name: panel.corr.sample
# Desc: panel function for plot.corr.sample
# Auth: Alexander.Ploner@ki.se 200206
#
# Chng: 060306 added check for empty panel (sigh)
#
{
    if (missing(col.line)) 
        col.line = trellis.par.get("plot.line")$col
    if (missing(col.symbol)) 
        col.symbol = gray(0.4)    
    if (scatter) {
        panel.xyplot(x,y, col.symbol=col.symbol, ...)
    }        
    if (grid) {
        panel.grid(h=-1, v=-1)
    }
    if (refline) {
        panel.abline(h=0, lty=2, col.line="black",lwd=0.5)
    }
    # We skip here if we have no data 
    if (length(x) <= 0) return()
    # Ok, still with us        
    xx = cbind(y,x)
    cc = CutCI(xx, number=nint, alpha=alpha)
    panel.xyplot(cc$x, cc$y, type="o", col.line=col.line, ...)
    if (ci) {
        panel.arrows(cc$x, cc$yci[,1], cc$x, cc$yci[,2], code=3,  
                     angle=90, length=length, col=col.line)
    }
    if (curve) {
        xx = seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length=101)
        xt = if (xlog) 10^(-x) else 1/x
        fit = lm(y~xt)
        yy = predict(fit, list(xt=if (xlog) 10^(-xx) else 1/xx))
        panel.xyplot(xx, yy, type="l",lty=1, lwd=2, col="red")
    }   
}

CIrho = function(rho, alpha=0.95)
#
# Name: CIrho
# Desc: computes a confidence interval for a set of correlation coefficients
{
    rho = rho[!is.na(rho)]
    nn = length(rho)
    qq = qnorm((1+alpha)/2)    
    mn = mean(rho)
    se = sd(rho)/sqrt(nn)
    CL=mn-qq*se
    CU=mn+qq*se

    c(CL=CL, CU=CU)
}    

## Combine the Cut-function with CI
CutCI = function(dat, number=10, func=mean, alpha=0.95)
{
    x = dat[,2]
    y = dat[,1]
    int = co.intervals(x, number=number, overlap=0)
    cutt = c(int[1,1], int[,2])
    xgrp = cut(x, cutt)
    ysum = tapply(y, xgrp, func, na.rm=TRUE)
    yci  = matrix(unlist(tapply(y, xgrp, CIrho, alpha=alpha)), ncol=2, byrow=TRUE)
    colnames(yci) = c("CL","CU")
    rownames(yci) = names(ysum)
    xmid = (int[,1]+int[,2])/2
    invisible(list(x=xmid, y=ysum, yci=yci))
}        

\name{boomSim}
\alias{boomSim}
main title here 
\title{
main title here 
}
\description{
main title here 
}
\usage{
boom(obs = obs, d = d, X = X, riskcat = riskcat, risklab = risklab, deltaCPE = deltaCPE, FinalNum = FinalNum, minatrisk = minatrisk, minevent = minevent)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{
%%     ~~Describe \code{obs} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{riskcat}{
%%     ~~Describe \code{riskcat} here~~
}
  \item{risklab}{
%%     ~~Describe \code{risklab} here~~
}
  \item{deltaCPE}{
%%     ~~Describe \code{deltaCPE} here~~
}
  \item{FinalNum}{
%%     ~~Describe \code{FinalNum} here~~
}
  \item{minatrisk}{
%%     ~~Describe \code{minatrisk} here~~
}
  \item{minevent}{
%%     ~~Describe \code{minevent} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (obs = obs, d = d, X = X, riskcat = riskcat, risklab = risklab, 
    deltaCPE = deltaCPE, FinalNum = FinalNum, minatrisk = minatrisk, 
    minevent = minevent) 
{
    X <- as.matrix(X)
    p <- dim(X)[2]
    i <- order(obs, 1 - d)
    obs <- obs[i]
    d <- d[i]
    fst <- !duplicated(obs)
    X <- X[i, , drop = F]
    riskcat <- riskcat[i]
    N <- length(obs)
    ux <- NULL
    jx <- NULL
    for (j in 1:p) {
        u <- sort(unique(X[, j]))
        if (length(u) > 1) {
            ux <- c(ux, u[-1] - diff(u)/2)
            jx <- c(jx, rep(j, length(u) - 1))
        }
    }
    if (any(table(riskcat) < minatrisk)) 
        print(paste("Starting partition size less than ", minatrisk, 
            sep = ""))
    if (any(table(riskcat[d == 1]) < minevent)) 
        print(paste("Starting partition event size less than ", 
            minevent, sep = ""))
    stopifnot(any(table(riskcat) < minatrisk) == FALSE)
    stopifnot(any(table(riskcat[d == 1]) < minevent) == FALSE)
    starthere <- .Call("findmaxCPE", as.numeric(obs), as.integer(d), 
        1 * riskcat, as.numeric(unique(c(riskcat))), as.numeric(rep(0, 
            N)), as.numeric(rep(1, N)), as.integer(rep(0, N)), 
        as.double(1e-09), as.double(.Machine$double.eps^0.75))
    cpecat <- starthere[[1]][1]
    riskcatORG <- riskcat
    flag <- 0
    node <- 1
    expandgroup <- list(list(riskcat = riskcat, catlabels = risklab, 
        splitGrp = NULL, splitVar = NULL, splitVal = NULL, cpecat = cpecat))
    currenthouse <- NULL
    for (gg in 1:length(unique(riskcatORG))) {
        cols <- sort(unique(riskcatORG))[gg]
        currenthouse <- c(currenthouse, list(list(InitRisk = cols, 
            Xinc = rep(999, dim(X)[2]))))
    }
    while (flag < 1) {
        currentrisk <- expandgroup[[node]]$riskcat
        listsplit <- NULL
        for (ii in 1:dim(X)[2]) {
            listsplit <- c(listsplit, any(apply(table(currentrisk, 
                t(t(X[, jx]) > ux)[, ii]) > minatrisk, 1, all) & 
                apply(table(currentrisk[d == 1], t(t(X[, jx]) > 
                  ux)[d == 1, ii]) > minevent, 1, all)))
        }
        if (!any(listsplit)) 
            flag <- 1
        if (!any(listsplit)) 
            next
        comball <- .Call("combinedEXPANDMAX", as.integer(currentrisk), 
            as.integer(unique(currentrisk)), as.numeric(X), as.numeric(ux), 
            as.integer(jx), as.integer(minatrisk), as.integer(minevent), 
            as.numeric(obs), as.integer(d), as.numeric(rep(0, 
                N)), as.numeric(rep(1, N)), as.integer(rep(0, 
                N)), as.double(1e-09), as.double(.Machine$double.eps^0.75))
        if (comball$metrics[1] == 999999) 
            flag <- 1
        if (comball$metrics[1] == 999999) 
            next
        if (comball$metrics[1] == 0) 
            flag <- 1
        if (comball$metrics[1] == 0) 
            next
        Grphit <- comball$metrics[2]
        jxhit <- comball$metrics[3]
        Valhit <- comball$metrics[4]
        Varhit <- paste("X", jxhit, sep = "")
        if ((comball$metrics[1] - expandgroup[[node]]$cpecat) < 
            deltaCPE) 
            flag <- 1
        if ((comball$metrics[1] - expandgroup[[node]]$cpecat) < 
            deltaCPE) 
            next
        newhouse <- NULL
        if (Grphit == min(currentrisk)) {
            newlabels <- c(paste(expandgroup[[node]]$catlabels[Grphit], 
                "&", Varhit, "<", Valhit), paste(expandgroup[[node]]$catlabels[Grphit], 
                "&", Varhit, ">", Valhit), expandgroup[[node]]$catlabels[(Grphit + 
                1):(max(currentrisk))])
            lowgroup <- highgroup <- currenthouse[[Grphit]]$Xinc
            lowgroup[jxhit] <- 0
            highgroup[jxhit] <- 1
            newhouse <- c(list(list(InitRisk = currenthouse[[Grphit]]$InitRisk, 
                Xinc = lowgroup)), list(list(InitRisk = currenthouse[[Grphit]]$InitRisk, 
                Xinc = highgroup)), currenthouse[(Grphit + 1):(max(currentrisk))])
        }
        else {
            if (Grphit == max(currentrisk)) {
                newlabels <- c(expandgroup[[node]]$catlabels[1:(Grphit - 
                  1)], paste(expandgroup[[node]]$catlabels[Grphit], 
                  "&", Varhit, "<", Valhit), paste(expandgroup[[node]]$catlabels[Grphit], 
                  "&", Varhit, ">", Valhit))
                lowgroup <- highgroup <- currenthouse[[Grphit]]$Xinc
                lowgroup[jxhit] <- 0
                highgroup[jxhit] <- 1
                newhouse <- c(currenthouse[1:(Grphit - 1)], list(list(InitRisk = currenthouse[[Grphit]]$InitRisk, 
                  Xinc = lowgroup)), list(list(InitRisk = currenthouse[[Grphit]]$InitRisk, 
                  Xinc = highgroup)))
            }
            else {
                newlabels <- c(expandgroup[[node]]$catlabels[1:(Grphit - 
                  1)], paste(expandgroup[[node]]$catlabels[Grphit], 
                  "&", Varhit, "<", Valhit), paste(expandgroup[[node]]$catlabels[Grphit], 
                  "&", Varhit, ">", Valhit), expandgroup[[node]]$catlabels[(Grphit + 
                  1):(max(currentrisk))])
                lowgroup <- highgroup <- currenthouse[[Grphit]]$Xinc
                lowgroup[jxhit] <- 0
                highgroup[jxhit] <- 1
                newhouse <- c(currenthouse[1:(Grphit - 1)], list(list(InitRisk = currenthouse[[Grphit]]$InitRisk, 
                  Xinc = lowgroup)), list(list(InitRisk = currenthouse[[Grphit]]$InitRisk, 
                  Xinc = highgroup)), currenthouse[(Grphit + 
                  1):(max(currentrisk))])
            }
        }
        currenthouse <- newhouse
        expandgroup <- c(expandgroup, list(list(riskcat = comball$newriskgrp, 
            catlabels = newlabels, splitGrp = Grphit, splitVar = Varhit, 
            splitVal = Valhit, cpecat = comball$metrics[1])))
        node <- node + 1
    }
    lastnode <- length(expandgroup)
    mingroupnum <- FinalNum
    riskcat <- expandgroup[[lastnode]]$riskcat
    risklab <- expandgroup[[lastnode]]$catlabels
    reducegroup <- list(list(riskcat = riskcat, catlabels = risklab, 
        combGrp = NULL))
    reducehouse <- currenthouse
    node <- 1
    flag <- 0
    if (length(unique(riskcat)) <= FinalNum) 
        flag = 1
    if (length(unique(riskcat)) < FinalNum) 
        print(paste("Growing did not reach ", FinalNum, " categories", 
            sep = ""))
    if (length(unique(riskcat)) == FinalNum) 
        print(paste("Growing only reached ", FinalNum, " categories", 
            sep = ""))
    while (flag < 1 & length(expandgroup) > 1) {
        currentrisk <- reducegroup[[node]]$riskcat
        mergefind <- .Call("combinedMERGEMAX", as.integer(currentrisk), 
            as.integer(unique(sort(currentrisk))), as.numeric(obs), 
            as.integer(d), as.numeric(rep(0, N)), as.numeric(rep(1, 
                N)), as.integer(rep(0, N)), as.double(1e-09), 
            as.double(.Machine$double.eps^0.75))
        reducegroup <- c(reducegroup, list(list(riskcat = mergefind$newriskgrp, 
            catlabels = NA)))
        newhouse <- NULL
        removedgrp <- c(mergefind$metrics[2], mergefind$metrics[3])
        removedlist <- list(list(InitRisk = c(reducehouse[[removedgrp[1]]]$InitRisk, 
            reducehouse[[removedgrp[2]]]$InitRisk), Xinc = cbind(reducehouse[[removedgrp[1]]]$Xinc, 
            reducehouse[[removedgrp[2]]]$Xinc)))
        if (removedgrp[1] == min(currentrisk)) {
            newhouse <- c(removedlist, reducehouse[-c(removedgrp)])
        }
        else {
            if (removedgrp[1] == (max(currentrisk) - 1)) {
                newhouse <- c(reducehouse[-c(removedgrp)], removedlist)
            }
            else {
                newhouse <- c(reducehouse[1:(removedgrp[1] - 
                  1)], removedlist, reducehouse[((removedgrp[1] + 
                  1):max(currentrisk))[((removedgrp[1] + 1):max(currentrisk)) != 
                  removedgrp[2]]])
            }
        }
        reducehouse <- newhouse
        if (max(mergefind$newriskgrp) == mingroupnum) 
            flag <- 1
        node <- node + 1
    }
    newgroup <- reducegroup[[length(reducegroup)]]$riskcat
    newlabels <- reducegroup[[length(reducegroup)]]$catlabels
    finalhouse <- NULL
    for (bb in 1:length(reducehouse)) {
        Xmatrix <- reducehouse[[bb]]$Xinc
        orggrps <- reducehouse[[bb]]$InitRisk
        if (is.null(dim(Xmatrix)[2])) 
            finalhouse <- c(finalhouse, list(reducehouse[[bb]]))
        if (is.null(dim(Xmatrix)[2])) 
            next
        if (dim(Xmatrix)[2] == 1) 
            finalhouse <- c(finalhouse, list(reducehouse[[bb]]))
        if (dim(Xmatrix)[2] == 1) 
            next
        newmerge <- 1
        while (newmerge == 1) {
            newmerge = 0
            if (is.null(dim(Xmatrix)[2])) 
                next
            if (dim(Xmatrix)[2] == 1) 
                next
            for (i in 1:(dim(Xmatrix)[2] - 1)) {
                for (j in (i + 1):dim(Xmatrix)[2]) {
                  if (newmerge == 1) 
                    next
                  test <- Xmatrix[, i] != Xmatrix[, j]
                  possmerge <- (sum(test) < 2) & all(Xmatrix[test, 
                    i] != 999) & all(Xmatrix[test, j] != 999) & 
                    (orggrps[i] == orggrps[j])
                  if (possmerge) {
                    Xmatrix[Xmatrix[, i] != Xmatrix[, j], i] <- 999
                    Xmatrix <- Xmatrix[, -j]
                    orggrps <- orggrps[-j]
                    newmerge <- 1
                  }
                }
            }
        }
        finalhouse <- c(finalhouse, list(list(InitRisk = orggrps, 
            Xinc = Xmatrix)))
    }
    resultslabels <- finalhouse
    return(resultslabels)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

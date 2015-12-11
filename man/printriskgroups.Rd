\name{printriskgroups}
\alias{printriskgroups}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main title here 
}
\description{
Main title here 
}
\usage{
printriskgroups(finalhouse)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{finalhouse}{
%%     ~~Describe \code{finalhouse} here~~
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
function (finalhouse) 
{
    curlabmatrix <- NULL
    for (i in 1:length(finalhouse)) {
        for (j in 1:length(finalhouse[[i]]$InitRisk)) {
            curlab <- paste("New Risk Group ", i, "=", risklab[finalhouse[[i]]$InitRisk[j]], 
                sep = "")
            if (length(finalhouse[[i]]$InitRisk) == 1) {
                if (any(finalhouse[[i]]$Xinc == 0)) {
                  lablow <- paste("X", (1:length(finalhouse[[i]]$Xinc))[finalhouse[[i]]$Xinc == 
                    0], "=0", sep = "")
                  for (kk in 1:sum(finalhouse[[i]]$Xinc == 0)) {
                    curlab <- paste(curlab, " And ", lablow[kk], 
                      sep = "")
                  }
                }
                if (any(finalhouse[[i]]$Xinc == 1)) {
                  lablow <- paste("X", (1:length(finalhouse[[i]]$Xinc))[finalhouse[[i]]$Xinc == 
                    1], "=1", sep = "")
                  for (kk in 1:sum(finalhouse[[i]]$Xinc == 1)) {
                    curlab <- paste(curlab, " And ", lablow[kk], 
                      sep = "")
                  }
                }
            }
            else {
                if (any(finalhouse[[i]]$Xinc[, j] == 0)) {
                  lablow <- paste("X", (1:length(finalhouse[[i]]$Xinc[, 
                    j]))[finalhouse[[i]]$Xinc[, j] == 0], "=0", 
                    sep = "")
                  for (kk in 1:sum(finalhouse[[i]]$Xinc[, j] == 
                    0)) {
                    curlab <- paste(curlab, " And ", lablow[kk], 
                      sep = "")
                  }
                }
                if (any(finalhouse[[i]]$Xinc[, j] == 1)) {
                  lablow <- paste("X", (1:length(finalhouse[[i]]$Xinc[, 
                    j]))[finalhouse[[i]]$Xinc[, j] == 1], "=1", 
                    sep = "")
                  for (kk in 1:sum(finalhouse[[i]]$Xinc[, j] == 
                    1)) {
                    curlab <- paste(curlab, " And ", lablow[kk], 
                      sep = "")
                  }
                }
            }
            curlabmatrix <- rbind(curlabmatrix, curlab)
        }
    }
    print(curlabmatrix, rownames = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

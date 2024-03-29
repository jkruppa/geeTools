% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gee.R
\name{gee_var_select}
\alias{gee_var_select}
\title{Small wrapper function of the geesmv package}
\usage{
gee_var_select(formula, idvar, data, varianceEstimator = NULL)
}
\arguments{
\item{formula}{The formula of the analysis}

\item{idvar}{The id cluster variable}

\item{data}{The data set}

\item{varianceEstimator}{Choose one of the following list:
\describe{
\item{\emph{FayGrauband2001}}{}
\item{\emph{GoshoSatoTakeuchi2014}}{}
\item{\emph{KauermannCarroll2001}}{}
\item{\emph{MorelBokossaNeerchal2003}}{}
\item{\emph{ManclDeRouen2001}}{}
\item{\emph{Mackinnon1985}}{}
\item{\emph{Pan2001}}{}
\item{\emph{WangLong2011}}{}
##' }}
}
\value{
\describe{
\item{\emph{type}}{Which sandwich variance estimator was chosen?}
\item{\emph{df.residual}}{Degree of freedom of the residuals}
\item{\emph{beta}}{Effect of the treatments}
\item{\emph{vbeta}}{Variance/Covariance matrix of beta} 
}
}
\description{
Small wrapper function of the geesmv package
}
\details{
Small wrapper function of the geesmv package
}
\examples{
require(multcomp)
geeFit <- gee_var_select(formula = "resp ~ 0 + treat",
                         idvar = "litter",
                         data = litterGeeDf,
                         varianceEstimator = "WangLong2011")

summary(geeFit)

contrMat_n <- setNames(rep(1, length(names(geeFit$beta))),
                       names(geeFit$beta))

ggGee <- glht(parm(coef = geeFit$beta, vcov = geeFit$vbeta), 
              linfct = contrMat(contrMat_n, type = "Dunnett"))
ggGee$df <- geeFit$df.residual

summary(ggGee)
confint(ggGee)

}
\author{
Jochen Kruppa
}


##' Small wrapper function of the geesmv package
##'
##' Small wrapper function of the geesmv package
##' @title Small wrapper function of the geesmv package
##' @param formula The formula of the analysis
##' @param idvar The id cluster variable
##' @param data The data set
##' @param varianceEstimator
##' Choose one of the following list:
##' \describe{
##' \item{\emph{FayGrauband2001}}{}
##' \item{\emph{GoshoSatoTakeuchi2014}}{}
##' \item{\emph{KauermannCarroll2001}}{}
##' \item{\emph{MorelBokossaNeerchal2003}}{}
##' \item{\emph{ManclDeRouen2001}}{}
##' \item{\emph{Mackinnon1985}}{}
##' \item{\emph{Pan2001}}{}
##' \item{\emph{WangLong2011}}{}
##' }
##' @return
##' \describe{
##' \item{\emph{type}}{Which sandwich variance estimator was chosen?}
##' \item{\emph{df.residual}}{Degree of freedom of the residuals}
##' \item{\emph{beta}}{Effect of the treatments}
##' \item{\emph{vbeta}}{Variance/Covariance matrix of beta} 
##' }
##' @author Jochen Kruppa
##' @export
##' @examples
##' require(multcomp)
##' geeFit <- gee_var_select(formula = "resp ~ 0 + treat",
##'                          idvar = "litter",
##'                          data = litterGeeDf,
##'                          varianceEstimator = "WangLong2011")
##' 
##' summary(geeFit)
##' 
##' contrMat_n <- setNames(rep(1, length(names(geeFit$beta))),
##'                        names(geeFit$beta))
##' 
##' ggGee <- glht(parm(coef = geeFit$beta, vcov = geeFit$vbeta), 
##'               linfct = contrMat(contrMat_n, type = "Dunnett"))
##' ggGee$df <- geeFit$df.residual
##' 
##' summary(ggGee)
##' confint(ggGee)
##' 
gee_var_select <- function(formula, idvar, data, varianceEstimator = NULL){
  require(geesmv)
  require(geepack)
  if(any(is.na(data))) {
    stop("No missing values are allowed!")
  }
  varEstList <- list("FayGrauband2001" = "fg",
                     "GoshoSatoTakeuchi2014" = "gst",
                     "KauermannCarroll2001" = "kc",
                     "LiangZeger1986" = "lz",
                     "MorelBokossaNeerchal2003" = "mbn",
                     "ManclDeRouen2001" = "md",
                     "Mackinnon1985" = "mk",
                     "Pan2001" = "pan",
                     "WangLong2011" = "wl")
  if(is.null(varianceEstimator)){
    message("Specifiy a sandwich variance estimator!\n Choose one of the following: ",
            paste(varEstList, collapse = ", "), ".\n Specified by:")
    print(varEstList)
  }
  fit <- geeglm(as.formula(formula), id = data[[idvar]], data,
                  family = poisson(), corstr="exch")
  if(varianceEstimator == "default"){
    estimates <- list(type = "default.gee.estimation",
                      df.residual = df.residual(fit),
                      beta = coef(fit),
                      vbeta = summary(fit)$cov.scaled)
  } else {
    stringFit <- paste0("GEE.var.", varEstList[varianceEstimator],
                        "(", formula, ", id = '", idvar, "', data, family = poisson)")
    capture.output(suppressMessages(vcovEst <- eval(parse(text = stringFit))))
    estimates <- list(type = paste0(varEstList[varianceEstimator],".gee.estimation_",
                                    varianceEstimator),
                      df.residual = df.residual(fit),
                      beta = coef(fit),
                      vbeta = diag(vcovEst$cov.beta))        
  }
  return(estimates)    
}

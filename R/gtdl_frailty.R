
#' GTDL regression model with covariate in the alpha parameter
#' @param a0 Intercept
#' @param a1 Regression coefficient
#' @param lambda scale>0
#' @param beta Regression coefficient
#'
#' @usage dens(t,x,param)
#' @usage surv(t,x,param)
#' @usage hazard(t,x,param)
#' @usage long_term(t,x,param)
#'
#'
#' @details  Survival function  S(t|x)=[\frac{1+exp(alpha t+x \beta)}{1+exp(x \beta)}]
#'
#' @references Calsavara, V. F., Milani, E. A., Bertolli, E., & Tomazella, V. (2020). Long-term frailty modeling using a non-proportional hazards model: Application with a melanoma dataset. Statistical Methods in Medical Research, 29(8), 2100-2118.
#'
#' @examples library(GTDLfrailty)
#' @examples  t <- seq(0,20,by = 0.1)
#' @examples  x <- 1
#' @examples  a0 <- -1.00
#' @examples  a1 <- -0.05
#' @examples  lambda <- 2
#' @examples  beta <- 1
#' @examples  param <- c(a0,a1,lambda,beta)
#' @examples  y1 <- dens(t,x,param)
#' @examples  y2 <- surv(t,x,param)
#' @examples  y3 <- hazard(t,x,param)
#' @examples  tt <- as.matrix(cbind(t,t,t))
#' @examples  yy <- as.matrix(cbind(y1,y2,y3))
#' @examples  matplot(tt,yy,type="l",xlab="Time",ylab="",ylim=c(0,1),lty = 1:3,col=1:3,lwd=2)
#' @export




dens=function(t,x,param){
  a=param[1]+param[2]*x

  f=surv(t,x,param)*hazard(t,x,param)
  return(f)

}


surv=function(t,x,param){
  a=param[1]+param[2]*x

  surv=( (1+exp(a*t+x*param[4]) )/( 1+exp(x*param[4]) ) )^(-param[3]/a)
  return(surv)

}


hazard=function(t,x,param){
  a=param[1]+param[2]*x

  hazard=param[3]*(exp(a*t+x*param[4]))/( 1+exp(x*param[4]) )
  return(hazard)

}



long_term=function(x,param){
  a=param[1]+param[2]*x

  cured=ifelse(a<0, ( 1+exp(x*param[4]) )^(param[3]/a),0)
  return(cured)

}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gtdl_frailty.R
\name{dens}
\alias{dens}
\title{GTDL regression model with covariate in the alpha parameter}
\usage{
dens(t,x,param)

surv(t,x,param)

hazard(t,x,param)

long_term(t,x,param)
}
\arguments{
\item{a0}{Intercept}

\item{a1}{Regression coefficient}

\item{lambda}{scale>0}

\item{beta}{Regression coefficient}
}
\description{
GTDL regression model with covariate in the alpha parameter
}
\details{
Survival function  S(t|x)=[\frac{1+exp(alpha t+x \beta)}{1+exp(x \beta)}]
}
\examples{
library(GTDLfrailty)
 t <- seq(0,20,by = 0.1)
 x <- 1
 a0 <- -1.00
 a1 <- -0.05
 lambda <- 2
 beta <- 1
 param <- c(a0,a1,lambda,beta)
 y1 <- dens(t,x,param)
 y2 <- surv(t,x,param)
 y3 <- hazard(t,x,param)
 tt <- as.matrix(cbind(t,t,t))
 yy <- as.matrix(cbind(y1,y2,y3))
 matplot(tt,yy,type="l",xlab="Time",ylab="",ylim=c(0,1),lty = 1:3,col=1:3,lwd=2)
}
\references{
Calsavara, V. F., Milani, E. A., Bertolli, E., & Tomazella, V. (2020). Long-term frailty modeling using a non-proportional hazards model: Application with a melanoma dataset. Statistical Methods in Medical Research, 29(8), 2100-2118.
}

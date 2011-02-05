#Reproduce Table 20.1 from
#big Wooldridge

require(stats4)
read.recid <- function(){
  return(read.table("recid.raw",col.names=c("black",
                                  "alcohol",
                                  "drugs",
                                  "super",
                                  "married",
                                  "felon",
                                  "workprg",
                                  "property",
                                  "person",
                                  "priors",
                                  "educ",
                                  "rules",
                                  "age",
                                  "tserved",
                                  "follow",
                                  "durat",
                                  "cens",
                                  "ldurat")))
}

depvars <- c("workprg","priors","tserved","felon","alcohol",
             "drugs","black","married","educ","age","constant")

ll.factory <- function(){
  recid <- read.recid()
  recid$constant <- 1
  y <- recid$durat
  X <- as.matrix(recid[,depvars])
  censored <- recid$cens==1
  N <- length(depvars)

  log.weibull.density <- Vectorize(function(t,xbeta,alpha){
    g <- exp(xbeta)
    xbeta+log(alpha*t^(alpha-1))-g*t^alpha
  },c("t","xbeta"))
 

  log.weibull.survival <- Vectorize(function(t,xbeta,alpha){
    g <- exp(xbeta)
    -g*t^alpha
  },c("t","xbeta"))
  
  negloglik <- function(alpha,workprg,priors,tserved,felon,alcohol,
                        drugs,black,married,educ,age,constant){
    beta <- c(workprg,priors,tserved,felon,alcohol,
              drugs,black,married,educ,age,constant)
    X.beta <- X%*%beta

    -sum(log.weibull.density(y[!censored],X.beta[!censored],alpha))                 -sum(log.weibull.survival(y[censored],X.beta[censored],alpha))
  }

  
  negloglik
}

main <- function(){
  negloglik <- ll.factory()
  start <- c(1,rep(1e-4,length(depvars)))
  names(start) <- c("alpha",depvars)
  mle(negloglik,as.list(start))
}

    

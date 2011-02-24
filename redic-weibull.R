#Reproduce Table 20.1 from
#big Wooldridge

require(stats4)
dyn.load(paste("recid_c" ,.Platform$dynlib.ext,sep=""))

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
  censored <- recid$cens==1
  ycens <- recid$durat[censored]
  yuncens <- recid$durat[!censored]
  Xcens <- as.matrix(recid[censored,depvars])
  Xuncens <- as.matrix(recid[!censored,depvars])
  N <- length(depvars)

  handle<-.Call("jmoy_loaddata",N,
                as.double(ycens),as.double(Xcens),length(ycens),
                as.double(yuncens),as.double(Xuncens),length(yuncens))
  
  negloglik <- function(workprg,priors,tserved,felon,alcohol,
                        drugs,black,married,educ,age,constant,alpha){
    beta <- c(workprg,priors,tserved,felon,alcohol,
              drugs,black,married,educ,age,constant)
    .Call("jmoy_negloglik",handle,as.double(beta),as.double(alpha))
 }

  gradient <- function(pars){
    .Call("jmoy_gradient",handle,as.double(pars[1:N]),as.double(pars[N+1]))
 }


  list(nll=negloglik,gr=gradient)
}

main <- function(){
  fns <- ll.factory()
  truestart <- c(0.091,0.089,0.014,-0.299,0.447,0.281,0.454,-0.152,
             -0.023,-0.0037,-3.402,0.806)
  start <- c(rep(0,length(depvars)),1)
  names(start) <- c(depvars,"alpha")
  mle(fns$nll,as.list(start),gr=fns$gr,
           control=list(trace=1))
}

    

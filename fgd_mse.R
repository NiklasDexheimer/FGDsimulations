x=c("RColorBrewer","lattice","plotly","gridExtra","colorspace","Matrix","purrr","SLOPE",
    "sigmoid","Rcpp","grDevices","sfsmisc")
lapply(x, require, character.only = TRUE)
#set seed for reproducability
set.seed(0)
#set dimension, number of iterations and amount of training data
d=30
n=1000000
it=10
errsgd=rep(0,n*it)
errsfgd=rep(0,n*it)
errmfgd=rep(0,n*it)
errnfgd=rep(0,n*it)
errnsfgd=rep(0,n*it)
theta_star=runif(d,-10,10)
for(m in (1:it)){
  #generate training data
  x=runif(d*n,-sqrt(3),sqrt(3))
  x=matrix(x,nrow=n)
  ep=rnorm(n)
  y=x%*%theta_star+ep
  theta_sgd=rep(0,d)
  #compute SGD estimator
  for(i in (1:n)){
    theta_sgd=theta_sgd+(1/(d+i))*(y[i]-as.numeric(x[i,]%*%theta_sgd))*x[i,]
    errsgd[(m-1)*n+i]=sum((theta_sgd-theta_star)^2)
  }
  theta_sfgd=rep(0,d)
  #compute FGD(1) estimator
  for(i in (1:n)){
    xi=rnorm(d)
    theta_sfgd=theta_sfgd+(1/(d^2+i))*(y[i]-as.numeric(x[i,]%*%theta_sfgd))*
      as.numeric(x[i,]%*%xi)*xi
    errsfgd[(m-1)*n+i]=sum((theta_sfgd-theta_star)^2)
  }
  theta_mfgd=rep(0,d)
  #Compute FGD(d) estimator
  for(i in (1:n)){
    for(k in (0:(d-1))){
      xi=rnorm(d)
      theta_mfgd=theta_mfgd+(1/(d*(d+i)))*(y[i]-as.numeric(x[i,]%*%theta_mfgd))*
        as.numeric(x[i,]%*%xi)*xi
    }
    errmfgd[(m-1)*n+i]=sum((theta_mfgd-theta_star)^2)
  }
  #compute aFGD(d) estimator
  theta_nfgd=rep(0,d)
  for(i in (1:n)){
    for(k in (0:(d-1))){
      xi=rnorm(d)
      theta_nfgd=theta_nfgd+(sqrt(pi/2)/(d*(d+i)))*(y[i]-as.numeric(x[i,]%*%theta_nfgd))*sqrt(sum(x[i,]^2))*
        sign(as.numeric(x[i,]%*%xi))*xi
    }
    errnfgd[(m-1)*n+i]=sum((theta_nfgd-theta_star)^2)
  }
    theta_nsfgd=rep(0,d)
    #compute aFGD(1) estimator
    for(i in (1:n)){
        xi=rnorm(d)
        theta_nsfgd=theta_nsfgd+(sqrt(pi/2)/(d^2+i))*(y[i]-as.numeric(x[i,]%*%theta_nsfgd))*sqrt(sum(x[i,]^2))*
          sign(as.numeric(x[i,]%*%xi))*xi
      errnsfgd[(m-1)*n+i]=sum((theta_nsfgd-theta_star)^2)
  }
}
#generate plots
plot((1:n),errsgd[1:n],type="l",xlab="Steps (log-scale)",ylab="MSE (log-scale)",lwd=0.15,xaxs="i",log="xy",axes="false")
for(i in (2:it)){
  lines((1:n),errsgd[((i-1)*n+1):(i*n)],lwd=0.15)
}
for(i in (1:it)){
  lines((1:n),errsfgd[((i-1)*n+1):(i*n)],col="blue",lwd=0.15)
}
for(i in (1:it)){
  lines((1:n),errmfgd[((i-1)*n+1):(i*n)],col="red",lwd=0.15)
}
for(i in (1:it)){
  lines((1:n),errnfgd[((i-1)*n+1):(i*n)],col="green",lwd=0.15)
}
for(i in (1:it)){
  lines((1:n),errnsfgd[((i-1)*n+1):(i*n)],col="orange",lwd=0.15)
}
eaxis(side=1,at.small=FALSE)
eaxis(side=2,at.small=FALSE,at=c(10^(-6),10^(-5),10^(-3),10^(-1),10,1000))
curve(d^2/x,lty=2,lwd=2,col="purple",add=TRUE,from=1)
curve(d/x,lty=3,lwd=2,col="purple",add=TRUE,from=1)
legend(1,10^(-0.5),cex=0.8,legend=c("SGD","FGD(1)","FGD(d)","aFGD(1)","aFGD(d)","d^2/n","d/n"),col=c("black","blue","red","orange","green","purple","purple"),lty=c(1,1,1,1,1,2,3), lwd=1.5)





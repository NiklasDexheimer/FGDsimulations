x=c("RColorBrewer","lattice","plotly","gridExtra","colorspace","Matrix","purrr","SLOPE","sigmoid","Rcpp","grDevices","sfsmisc")
lapply(x, require, character.only = TRUE)
#set seed for reproducability
set.seed(1)
#compute a logarithmic grid for the dimensions
d=exp(seq(log(10),log(200),length.out=10))
d=round(d)
#amount of training data
n=50000
#number of iterations
it=10
errsgd=rep(0,length(d)*it)
errsfgd=rep(0,length(d)*it)
errmfgd=rep(0,length(d)*it)
for(j in (1:length(d))){
  theta_star=runif(d[j],-10,10)
  theta_star=theta_star/sqrt(sum(theta_star^2))
  for(m in (1:it)){
    #generate training data
x=runif(d[j]*n,-sqrt(3),sqrt(3))
x=matrix(x,nrow=n)
ep=rnorm(n)
y=x%*%theta_star+ep
theta_sgd=rep(0,d[j])
#compute SGD estimator
for(i in (1:n)){
  theta_sgd=theta_sgd+(1/(d[j]+i))*(y[i]-as.numeric(x[i,]%*%theta_sgd))*x[i,]
}
theta_sfgd=rep(0,d[j])
#compute FGD(1) estimator
for(i in (1:n)){
  xi=rnorm(d[j])
  theta_sfgd=theta_sfgd+(1/(d[j]^2+i))*(y[i]-as.numeric(x[i,]%*%theta_sfgd))*
    as.numeric(x[i,]%*%xi)*xi
}
errsgd[(j-1)*it+m]=mean((theta_sgd-theta_star)^2)
errsfgd[(j-1)*it+m]=mean((theta_sfgd-theta_star)^2)
theta_mfgd=rep(0,d[j])
#compute FGD(d) estimator
for(i in (1:n)){
  for(k in (0:(d[j]-1))){
    xi=rnorm(d[j])
  theta_mfgd=theta_mfgd+(1/(d[j]*(d[j]+i)))*(y[i]-as.numeric(x[i,]%*%theta_mfgd))*
    as.numeric(x[i,]%*%xi)*xi
  }
}
errmfgd[(j-1)*it+m]=mean((theta_mfgd-theta_star)^2)
}
}
err=matrix(0,nrow=3,ncol=length(d))
sds=matrix(0,nrow=3,ncol=length(d))
#compute average error and standard deviations
for(i in (1:length(d))){
  err[1,i]=mean(errsgd[(1+(i-1)*it):(i*it)])
  err[2,i]=mean(errsfgd[(1+(i-1)*it):(i*it)])
  err[3,i]=mean(errmfgd[(1+(i-1)*it):(i*it)])
}
for(i in (1:length(d))){
  sds[1,i]=sd(errsgd[(1+(i-1)*it):(i*it)])
  sds[2,i]=sd(errsfgd[(1+(i-1)*it):(i*it)])
  sds[3,i]=sd(errmfgd[(1+(i-1)*it):(i*it)])
}
#generate plots
logsds=(sds/err)/log(10)
tred=adjustcolor("red",alpha.f=0.2)
tblue=adjustcolor("blue",alpha.f=0.2)
tblack=adjustcolor("black",alpha.f=0.2)
par(mar=c(5,5,3,3))
plot(d,err[2,],type="l",ylim=range(err),xlab="Dimension (log-scale)",
     ylab="",xaxs="i",lwd=1.5,axes="false",log="xy")
polygon(c(d,rev(d)),c(err[1,],rev(exp((log(err[1,])+logsds[1,])))),col=tblack,border= NA )
polygon(c(d,rev(d)),c(err[1,],rev(exp(log(err[1,])-logsds[1,]))),col=tblack,border= NA )
polygon(c(d,rev(d)),c(err[3,],rev(exp(log(err[3,])+logsds[3,]))),col=tred,border= NA )
polygon(c(d,rev(d)),c(err[3,],rev(exp(log(err[3,])-logsds[3,]))),col=tred,border= NA )
polygon(c(d,rev(d)),c(err[2,],rev(exp(log(err[2,])+logsds[2,]))),col=tblue,border= NA )
polygon(c(d,rev(d)),c(err[2,],rev(exp(log(err[2,])-logsds[2,]))),col=tblue,border= NA )
lines(d,err[1,],col="black",lwd=1.5)
lines(d,err[3,],col="red",lwd=1.5)
lines(d,err[2,],col="blue",lwd=1.5)
lines(d,d/50000,col="darkgreen",lwd=1.5,lty=2)
eaxis(side=1,at.small=FALSE,use.expr=FALSE)
eaxis(side=2,at.small=FALSE,use.expr=TRUE,at=c(10^-5,10^-4,10^-3,10^-2))
mtext(side=2,text="MSE/Dimension (log-scale)",line=3)


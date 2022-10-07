library(ggplot2)
library(tidyjson)
library(readr)
library(bain)
library(dplyr)
library(xlsx)
options(scipen =999)
setwd("C:/Users/svill/Desktop/mg/bayesian/ass")
df <- read.csv("df.csv")
df<-df[-c(1)]
N<-(length(df$scale))
Y<-df$scale
X1<-df$prl
X2<-df$lgl
X3<-df$plc
#################### DEFINE INFORMATIVE PRIORS ####################
#B0
mu00<-5.41177
tau00<-0.24158 
#B1
mu10<-0.06598
tau10<-0.03578  
#B2
mu20<- -0.19384 
tau20<-0.03637 
#B3
mu30<-0.08784 
tau30<-0.03812 
#S2
alfa0<-.01
beta0<-.01
# DEFINE NUMBER OF ITERATIONS AND STORAGE OF RESULTS
n.iterations <- 1000000
results1 <- matrix(0,n.iterations,5)
results2 <- matrix(0,n.iterations,5)
colnames(results1)<-c("b0", "b1", "b2","b3", "var")
colnames(results2)<-c("b0", "b1", "b2","b3", "var")
# DEFINE INITIAL VALUES FOR ITERATIONS OF TWO CHAINS
b0=-3
b1=5
b2=-8
b3=.6
var=1
b0_2=3
b1_2=-5
b2_2=8
b3_2=-.6
var1=2
results1[1,] <- c(b0, b1, b2,b3, var)
results2[1,] <- c(b0_2, b1_2, b2_2,b3_2, var1)

# STORAGE OF RESULTADOS
resultados1 <- matrix(0,n.iterations,5)
resultados2 <- matrix(0,n.iterations,5)
colnames(resultados1)<-c("b0", "b1", "b2","b3", "var")
colnames(resultados2)<-c("b0", "b1", "b2","b3", "var")
# INITIAL VALUES ALREADY DEFINED
resultados1[1,] <- c(b0, b1, b2,b3, var)
resultados2[1,] <- c(b0_2, b1_2, b2_2,b3_2, var1)
#################### GIBBS SAMPLER WITH MH STEP FOR B1-INFORMATIVE PRIORS ####################
set.seed(6060714)
for(iteration in 2:n.iterations){
  mu01_1<-((sum(Y-b1*X1-b2*X2)/var) + (mu00/tau00)) / (N/var + 1/tau00)
  mu01_2<-((sum(Y-b1_2*X1-b2_2*X2)/var1) + (mu00/tau00)) / (N/var1 + 1/tau00)
  tau01_1<-1 / (N/var + 1/tau00) 
  tau01_2<-1 / (N/var1 + 1/tau00)
  
  b0   <-rnorm(1, mu01_1, sqrt(tau01_1))
  b0_2 <-rnorm(1, mu01_2, sqrt(tau01_2))
  
  betaStar<-rnorm(1,b1,1) ;betaStar
  ref<-runif(1,0,1) ;ref
  yStar<-sum(X1*(Y-b0-b2*X2)) ;yStar
  mnStar1<-(yStar/var+ mu10/tau10)/(sum(X1^2)/var+1/tau10)
  sgmStar1<-1/sqrt((sum(X1^2)/var+1/tau10))
  pc<-dnorm(b1,b1,1)
  pn<-dnorm(betaStar,b1,1)
  tc<-dnorm(b1,mnStar1,sgmStar1)
  tn<-dnorm(betaStar,mnStar1,sgmStar1)
  prob<-tn/tc*pc/pn
  if (ref<prob)(b1<-betaStar)
  
  betaStar2<-rnorm(1,b1_2,1) ;betaStar2
  ref2<-runif(1,0,1) ;ref2
  yStar2<-sum(X1*(Y-b0_2-b2_2*X2)) ;yStar2
  mnStar12<-(yStar2/var1+ mu10/tau10)/(sum(X1^2)/var1+1/tau10)
  sgmStar12<-1/sqrt((sum(X1^2)/var1+1/tau10))
  pc2<-dnorm(b1_2,b1_2,1)
  pn2<-dnorm(betaStar2,b1_2,1)
  tc2<-dnorm(b1_2,mnStar12,sgmStar12)
  tn2<-dnorm(betaStar2,mnStar12,sgmStar12)
  prob2<-tn2/tc2*pc2/pn2
  if (ref2<prob2)(b1_2<-betaStar2)
  
  mu21      <- ((sum(X2*(Y-b0-b1*X1))/var) + (mu20/tau20)) / (sum(X2^2)/var + 1/tau20)
  mu21_2    <- ((sum(X2*(Y-b0_2-b1_2*X1))/var1) + (mu20/tau20)) / (sum(X2^2)/var1 + 1/tau20)
  tau21     <- 1 / (sum(X2^2)/var + 1/tau20)
  tau21_2   <- 1 / (sum(X2^2/var1 + 1/tau20))
  
  b2         <- rnorm(1, mu21, sqrt(tau21))
  b2_2       <- rnorm(1,mu21_2,sqrt(tau21_2))
  
  mu31<-  ((sum(X3*(Y-b0-b1*X1-b2*X2))/var)+(mu30/tau30))/(sum(X3^2)/var +1/tau30)
  mu31_2<-((sum(X3*(Y-b0_2-b1_2*X1-b2_2*X2))/var1)+(mu30/tau30))/(sum(X3^2)/var1 +1/tau30)
  tau31<- 1/(sum(X3^2)/var + 1/tau30)
  tau31_2<- 1/(sum(X3^2)/var1 + 1/tau30)
  
  b3<-rnorm(1,mu31,sqrt(tau31))
  b3_2<-rnorm(1,mu31_2,sqrt(tau31_2))
  
  alfa1  <- N/2 + alfa0                                   ########################
  beta1  <-   sum((Y-(b0+b1*X1+b2*X2+b3*X3))^2)/2 + beta0 ########################
  beta1_2  <- sum((Y-(b0_2+b1_2*X1+b2_2*X2+b3_2*X3))^2)/2 + beta0 ################
  
  var <- 1/rgamma(1, alfa1, beta1)
  var1 <- 1/rgamma(1, alfa1, beta1_2)
  
  results1[iteration,]<- c(b0, b1, b2,b3, var)
  results2[iteration,]<-c(b0_2,b1_2,b2_2,b3_2,var1)
}

results.1<-as.data.frame(results1)
results.1$g<-1
results.1$n<-1:1000000
results.2<-as.data.frame(results2)
results.2$g<-2
results.2$n<-1:1000000
chains<-bind_rows(results.1,results.2)
#################### DEFINE UNINFORMATIVE PRIORS ####################
#B0
mu00<-0
tau00<-1000
#B1
mu10<-0
tau10<-1000
#B2
mu20<-0
tau20<-1000
#B3
mu30<-0
tau30<-1000
#S2
alfa0<-.01
beta0<-.01
#################### GIBBS SAMPLER WITH MH STEP FOR B1 - UNINFORMATIVE PRIORS ####################
set.seed(18573)
for(iteration in 2:n.iterations){
  mu01_1<-((sum(Y-b1*X1-b2*X2)/var) + (mu00/tau00)) / (N/var + 1/tau00)
  mu01_2<-((sum(Y-b1_2*X1-b2_2*X2)/var1) + (mu00/tau00)) / (N/var1 + 1/tau00)
  tau01_1<-1 / (N/var + 1/tau00) #sqrt???
  tau01_2<-1 / (N/var1 + 1/tau00)
  
  b0   <-rnorm(1, mu01_1, sqrt(tau01_1))
  b0_2 <-rnorm(1, mu01_2, sqrt(tau01_2))
  
  betaStar<-rnorm(1,b1,1) ;betaStar
  ref<-runif(1,0,1) ;ref
  yStar<-sum(X1*(Y-b0-b2*X2)) ;yStar
  mnStar1<-(yStar/var+ mu10/tau10)/(sum(X1^2)/var+1/tau10)
  sgmStar1<-1/sqrt((sum(X1^2)/var+1/tau10))
  pc<-dnorm(b1,b1,1)
  pn<-dnorm(betaStar,b1,1)
  tc<-dnorm(b1,mnStar1,sgmStar1)
  tn<-dnorm(betaStar,mnStar1,sgmStar1)
  prob<-tn/tc*pc/pn
  if (ref<prob)(b1<-betaStar)
  
  betaStar2<-rnorm(1,b1_2,1) ;betaStar2
  ref2<-runif(1,0,1) ;ref2
  yStar2<-sum(X1*(Y-b0_2-b2_2*X2)) ;yStar2
  mnStar12<-(yStar2/var1+ mu10/tau10)/(sum(X1^2)/var1+1/tau10)
  sgmStar12<-1/sqrt((sum(X1^2)/var1+1/tau10))
  pc2<-dnorm(b1_2,b1_2,1)
  pn2<-dnorm(betaStar2,b1_2,1)
  tc2<-dnorm(b1_2,mnStar12,sgmStar12)
  tn2<-dnorm(betaStar2,mnStar12,sgmStar12)
  prob2<-tn2/tc2*pc2/pn2
  if (ref2<prob2)(b1_2<-betaStar2)
  
  mu21      <- ((sum(X2*(Y-b0-b1*X1))/var) + (mu20/tau20)) / (sum(X2^2)/var + 1/tau20)
  mu21_2    <- ((sum(X2*(Y-b0_2-b1_2*X1))/var1) + (mu20/tau20)) / (sum(X2^2)/var1 + 1/tau20)
  tau21     <- 1 / (sum(X2^2)/var + 1/tau20)
  tau21_2   <- 1 / (sum(X2^2/var1 + 1/tau20))
  
  b2         <- rnorm(1, mu21, sqrt(tau21))
  b2_2       <- rnorm(1,mu21_2,sqrt(tau21_2))
  
  mu31<-  ((sum(X3*(Y-b0-b1*X1-b2*X2))/var)+(mu30/tau30))/(sum(X3^2)/var +1/tau30)
  mu31_2<-((sum(X3*(Y-b0_2-b1_2*X1-b2_2*X2))/var1)+(mu30/tau30))/(sum(X3^2)/var1 +1/tau30)
  tau31<- 1/(sum(X3^2)/var + 1/tau30)
  tau31_2<- 1/(sum(X3^2)/var1 + 1/tau30)
  
  b3<-rnorm(1,mu31,sqrt(tau31))
  b3_2<-rnorm(1,mu31_2,sqrt(tau31_2))
  
  alfa1  <- N/2 + alfa0                                   ########################
  beta1  <-   sum((Y-(b0+b1*X1+b2*X2+b3*X3))^2)/2 + beta0 ########################
  beta1_2  <- sum((Y-(b0_2+b1_2*X1+b2_2*X2+b3_2*X3))^2)/2 + beta0 ################
  
  var <- 1/rgamma(1, alfa1, beta1)
  var1 <- 1/rgamma(1, alfa1, beta1_2)
  
  resultados1[iteration,]<- c(b0, b1, b2,b3, var)
  resultados2[iteration,]<-c(b0_2,b1_2,b2_2,b3_2,var1)
}

resultados.1<-as.data.frame(resultados1)
resultados.1$g<-1
resultados.1$n<-1:10000
resultados.2<-as.data.frame(resultados2)
resultados.2$g<-2
resultados.2$n<-1:10000
chains_uni<-bind_rows(resultados.1,resultados.2)

#################### BURN-IN PERIOD: TRACE PLOT- INFORMATIVE PRIORS ####################
cnv.b0<-ggplot(chains,aes(x=n,y=b0,group=g,color=factor(g)))+
  geom_line(size=1.5,alpha=.1)+geom_point(size=.5) ;cnv.b0
cnv.b1<-ggplot(chains,aes(x=n,y=b1,group=g,color=factor(g)))+
  geom_line(size=.5,alpha=.7)+geom_point(size=.5) ;cnv.b1 
cnv.b2<-ggplot(chains,aes(x=n,y=b2,group=g,color=factor(g)))+
  geom_line(size=1.5,alpha=.1)+geom_point(size=.5) ;cnv.b2 
cnv.b3<-ggplot(chains,aes(x=n,y=b3,group=g,color=factor(g)))+
  geom_line(size=1.5,alpha=.1)+geom_point(size=.5) ;cnv.b3 
cnv.va<-ggplot(chains,aes(x=n,y=var,group=g,color=factor(g)))+
  geom_line(size=1.5,alpha=.1)+geom_point(size=.5) ;cnv.va 

chains<-chains[-c(1000001:1900000),] 
chains<-chains[-c(1:900000),]

cnv.b0<-ggplot(chains,aes(x=n,y=b0,group=g,color=factor(g)))+
  geom_line(size=1.5,alpha=.1)+geom_point(size=.5) ;cnv.b0
cnv.b1<-ggplot(chains,aes(x=n,y=b1,group=g,color=factor(g)))+
  geom_line(size=.5,alpha=.7)+geom_point(size=.5) ;cnv.b1 
cnv.b2<-ggplot(chains,aes(x=n,y=b2,group=g,color=factor(g)))+
  geom_line(size=1.5,alpha=.1)+geom_point(size=.5) ;cnv.b2 
cnv.b3<-ggplot(chains,aes(x=n,y=b3,group=g,color=factor(g)))+
  geom_line(size=1.5,alpha=.1)+geom_point(size=.5) ;cnv.b3 
cnv.va<-ggplot(chains,aes(x=n,y=var,group=g,color=factor(g)))+
  geom_line(size=1.5,alpha=.1)+geom_point(size=.5) ;cnv.va 
#################### BURN-IN PERIOD: TRACE PLOT - UNINFORMATIVE PRIORS####################
  cnv.b0_uni<-ggplot(chains_uni,aes(x=n,y=b0,group=g,color=factor(g)))+
    geom_line(size=1.5,alpha=.1)+geom_point(size=.5) ;cnv.b0_uni #-3000
  cnv.b1_uni<-ggplot(chains_uni,aes(x=n,y=b1,group=g,color=factor(g)))+
    geom_line(size=.5,alpha=.7)+geom_point(size=.5) ;cnv.b1_uni
  cnv.b2_uni<-ggplot(chains_uni,aes(x=n,y=b2,group=g,color=factor(g)))+
    geom_line(size=1.5,alpha=.1)+geom_point(size=.5) ;cnv.b2_uni
  cnv.b3_uni<-ggplot(chains_uni,aes(x=n,y=b3,group=g,color=factor(g)))+
    geom_line(size=1.5,alpha=.1)+geom_point(size=.5) ;cnv.b3_uni 
  cnv.va_uni<-ggplot(chains_uni,aes(x=n,y=var,group=g,color=factor(g)))+
    geom_line(size=1.5,alpha=.1)+geom_point(size=.5) ;cnv.va_uni 
  
  chains_uni<-chains_uni[-c(1000001:1900000),] 
  chains_uni<-chains_uni[-c(1:900000),]
  
  cnv.b0_uni<-ggplot(chains_uni,aes(x=n,y=b0,group=g,color=factor(g)))+
    geom_line(size=1.5,alpha=.1)+geom_point(size=.5) ;cnv.b0_uni
  cnv.b1_uni<-ggplot(chains_uni,aes(x=n,y=b1,group=g,color=factor(g)))+
    geom_line(size=.5,alpha=.7)+geom_point(size=.5) ;cnv.b1_uni 
  cnv.b2_uni<-ggplot(chains_uni,aes(x=n,y=b2,group=g,color=factor(g)))+
    geom_line(size=1.5,alpha=.1)+geom_point(size=.5) ;cnv.b2_uni 
  cnv.b3_uni<-ggplot(chains_uni,aes(x=n,y=b3,group=g,color=factor(g)))+
    geom_line(size=1.5,alpha=.1)+geom_point(size=.5) ;cnv.b3_uni 
  cnv.va_uni<-ggplot(chains_uni,aes(x=n,y=var,group=g,color=factor(g)))+
    geom_line(size=1.5,alpha=.1)+geom_point(size=.5) ;cnv.va_uni
#################### COnVERGENCE: Autocorrelation plots ####################
acf(chains$b0)
acf(chains$b1)
acf(chains$b2)
acf(chains$b3)
acf(chains$var)

acf(chains_uni$b0)
acf(chains_uni$b1)
acf(chains_uni$b2)
acf(chains_uni$b3)
acf(chains_uni$var)
################### CONVERGENCE: MONTE CARLO ERROR - INFORMATIVE PRIORS ####################
#SD of estimates divided by the square root of the number of iterations
chains<-chains[-c(6:7)]
mce <- matrix(0,5,1)
sd1<-sqrt(200000)
mce[,1] <- apply(chains,2,sd)
mc_b0<-(mce[1,1]/sd1)
mc_b1<-(mce[2,1]/sd1)
mc_b2<-(mce[3,1]/sd1)
mc_b3<-(mce[4,1]/sd1)
mc_var<-(mce[5,1]/sd1)
mc<-rbind(mc_b0,mc_b1,mc_b2,mc_b3,mc_var) 
mce<-cbind(mce,mc)
colnames(mce) <- c("Sd","Monte Carlo error")
rownames(mce) <- c(c("b0", "b1", "b2","b3", "var")) 
mce
#################### CONVERGENCE: MONTE CARLO ERROR - UNINFORMATIVE PRIORS####################
#SD of estimates divided by the square root of the number of iterations
chains_uni<-chains_uni[-c(6:7)]
mce_uni <- matrix(0,5,1)
mce_uni [,1] <- apply(chains_uni,2,sd)
mc_b0_uni <-(mce_uni[1,1]/sd1)
mc_b1_uni <-(mce_uni[2,1]/sd1)
mc_b2_uni <-(mce_uni[3,1]/sd1)
mc_b3_uni <-(mce_uni[4,1]/sd1)
mc_var_uni <-(mce_uni[5,1]/sd1)
mc_uni <-rbind(mc_b0_uni,mc_b1_uni,mc_b2_uni,mc_b3_uni,mc_var_uni ) 
mce_uni <-cbind(mce_uni ,mc_uni )
colnames(mce_uni ) <- c("Sd","MC error")
rownames(mce_uni ) <- c(c("b0", "b1", "b2","b3", "var")) 
  mce_uni
#################### POSTERIOR PREDICTIVE P VALUE FOR NORMALITY ASSUMPTION -INFORMATIVE PRIORS ####################
Y.est<-matrix(0,1524,200000) #Posterior predictive distribution 
for(i in 1:200000){
  Yi<-chains[i,1]+chains[i,2]*X1+chains[i,3]*X2+chains[i,4]*X3
  Y.est[,i]<-Yi}
#Compute residuals for observed data
res.obs<-matrix(0,1524,200000)
for(i in 1:200000){
  r.obs<-Y-chains[i,1]-chains[i,2]*X1-chains[i,3]*X2-chains[i,4]*X3
  res.obs[,i]<-r.obs}
#Compute residuals for simulated data
res.sim<-matrix(0,1524,200000)
for(i in 1:200000){
  r.sim<-Y.est[i]-chains[i,1]-chains[i,2]*X1-chains[i,3]*X2-chains[i,4]*X3
  res.sim[,i]<-r.sim}
#Skewness of the observed data
skw.obs<-matrix(0,1,200000)
for (i in 1:200000){
  s.obs<-abs(3*(mean(res.obs[,i])-median(res.obs[,i])/sd(res.obs[,i])))
  skw.obs[,i]<-s.obs}
#Skewness of the simulated data
skw.sim<-matrix(0,1,200000)
for (i in 1:200000){
  s.sim<-abs(3*(mean(res.sim[,i])-median(res.sim[,i])/sd(res.sim[,i])))
  skw.sim[,i]<-s.sim}
#Compute how many times skewness is larger in simulated data than in observed data
test.stat<-matrix(0,1,200000)
for(i in 1:200000){
  if(skw.sim[i]>skw.obs[i])
  {    (test.stat[i]<-1);  }
  else{    (test.stat[i]<-0); }}

table(test.stat) # 0= 8036 1=191964 
mean(test.stat)  # 0.95892
#################### POSTERIOR PREDICTIVE P VALUE FOR NORMALITY ASSUMPTION - UNINFORMATIVE PRIORS ####################
Y.est.unin<-matrix(0,1524,200000) #Posterior predictive distribution 
for(i in 1:200000){
  Yi.uni<-chains_uni[i,1]+chains_uni[i,2]*X1+chains_uni[i,3]*X2+chains_uni[i,4]*X3
  Y.est.unin[,i]<-Yi.uni}
#Compute residuals for observed data
res.obs.uni<-matrix(0,1524,200000)
for(i in 1:200000){
  r.obs.uni<-Y-chains_uni[i,1]-chains_uni[i,2]*X1-chains_uni[i,3]*X2-chains_uni[i,4]*X3
  res.obs.uni[,i]<-r.obs.uni}
#Compute residuals for simulated data
res.sim.uni<-matrix(0,1524,200000)
for(i in 1:200000){
  r.sim.uni<-Y.est.unin[i]-chains_uni[i,1]-chains_uni[i,2]*X1-chains_uni[i,3]*X2-chains_uni[i,4]*X3
  res.sim.uni[,i]<-r.sim.uni}
#Skweness of the observed data
skw.obs.uni<-matrix(0,1,200000)
for (i in 1:200000){
  s.obs.uni<-abs(3*(mean(res.obs.uni[,i])-median(res.obs.uni[,i])/sd(res.obs.uni[,i])))
  skw.obs.uni[,i]<-s.obs.uni}
#Skweness of the simulated data
skw.sim.uni<-matrix(0,1,200000)
for (i in 1:200000){
  s.sim.uni<-abs(3*(mean(res.sim.uni[,i])-median(res.sim.uni[,i])/sd(res.sim.uni[,i])))
  skw.sim.uni[,i]<-s.sim.uni}
#Compute how many times skweness is larger in simulated data than in observed data
test.stat.uni<-matrix(0,1,200000)
for(i in 1:200000){
  if(skw.sim.uni[i]>skw.obs.uni[i])
  { (test.stat.uni[i]<-1);  }
  else{    (test.stat.uni[i]<-0); }}

table(test.stat.uni) # 0= 8529 1=191741 
mean(test.stat.uni)  # 0.958705
#################### RESULTS #################### 
chains<-chains[-c(6:7)]
output_inf<-matrix(0,5,4)
rownames(output_inf)<-c("b0", "b1", "b2","b3", "var")
colnames(output_inf)<-c("Mean","SD","Q025","Q975")
output_inf[,1]<-apply(chains,2,mean)
output_inf[,2]<-apply(chains,2,sd)
output_inf[,3]<-apply(chains,2,quantile,0.025)
output_inf[,4]<-apply(chains,2,quantile,0.975)
output_inf

chains_uni<-chains_uni[-c(6:7)]
output_unin<-matrix(0,5,4)
rownames(output_unin)<-c("b0", "b1", "b2","b3", "var")
colnames(output_unin)<-c("Mean","SD","Q025","Q975")
output_unin[,1]<-apply(chains_uni,2,mean)
output_unin[,2]<-apply(chains_uni,2,sd)
output_unin[,3]<-apply(chains_uni,2,quantile,0.025)
output_unin[,4]<-apply(chains_uni,2,quantile,0.975)
output_unin

#################### DIC- informative ###########################
##DHAT
y_data<-df$scale  #outcome from the data
x<-cbind(chains$b0,chains$b1,chains$b2,chains$b3) #sampled beta coefficients
post_mean<-(colMeans(x)) #mean of sampled coefficientes

post_var<-mean(chains$var) #mean of the sampled variances



dhat<- -2*sum(dnorm(y_data,mean=post_mean,sd=sqrt(post_var)),log=TRUE)
dhat

##DBAR
sample_lik<-rep(NA,nrow(chains))
#Estimate likelihood for all samples
for (i in 1:nrow(chains)){  
  sample_lik[i]<- -2*sum(dnorm(y_data,mean=(chains$b0[i]+chains$b1[i]*X1+
                        chains$b2[i]*X2+chains$b3[i]*X3),sd=sqrt(chains$var[i]), 
                        log=TRUE))
}
dbar<-mean(sample_lik)
dbar
## DIC
pd<-dbar-dhat
dic_inf<-dhat+2*pd ;dic_inf

#################### DIC- uninformative ###########################
x2<-cbind(chains_uni$b0,chains_uni$b1,chains_uni$b2,chains_uni$b3) #sampled beta coefficients
post_mean2<-(colMeans(x2)) #mean
post_var2<-mean(chains_uni$var) #mean of the sampled variance


dhat2<- -2*sum(dnorm(y_data,mean=post_mean2,sd=sqrt(post_var2)),log=TRUE)

dhat2

##DBAR
sample_lik2<-rep(NA,nrow(chains_uni))
#Estimate likelihood for all samples
for (i in 1:nrow(chains_uni)){  
  sample_lik2[i]<- -2*sum(dnorm(y_data,mean=(chains_uni$b0[i]+
                                              chains_uni$b1[i]*X1+
                                              chains_uni$b2[i]*X2+
                                              chains_uni$b3[i]*X3),
                               sd=sqrt(chains_uni$var[i]), log=TRUE))
}
dbar2<-mean(sample_lik2)
dbar2
## DIC
pd2<-dbar2-dhat2
dic_uninf<-dhat2+2*pd2 ;dic_uninf

#################### BAYESIAN HYPOTHESES EVALUATION ####################
m1<-lm(scale~prl+lgl+plc,data = df)
set.seed(43838)
#set.seed(78874) #This seed was also run to check stability
bf<-bain(m1," prl=0 & lgl=0 & plc=0;
              prl>0 & lgl>0 & plc>0; 
              prl>0 & lgl>0 & plc<0;
              prl>0 & lgl<0 & plc>0;
              prl>0 & lgl<0 & plc<0;
              prl<0 & lgl<0 & plc<0;
              prl<0 & lgl<0 & plc>0;
              prl<0 & lgl>0 & plc>0;
              prl>0 & lgl<0 & plc<0;
              prl<0 & lgl>0 & plc<0",fraction=1,standardize=F) ;bf



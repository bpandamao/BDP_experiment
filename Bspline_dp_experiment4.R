library(mgcv);library(mvtnorm);library(truncnorm);library(GIGrvg);library(invgamma);library(splines)
##experiment 3
XS=runif(1000,0,1)
X=matrix(XS,100,10)
bb=c(1,1,1,0,0,0,0,0,0,0)/sqrt(3)
ww=X%*%bb
Y=sin((ww-sqrt(3)/2-1.645/sqrt(12))*pi/(3.29/sqrt(12)))+0.1*rnorm(100)

##experiment 4
sigma2=matrix(rep(0,25),5,5)
for (i in 1:5){
  for (j in 1:5){
    ep=abs(j-i)
    sigma2[i,j]=0.5^ep
  }
}
X=rmvnorm(100,rep(0,5),sqrt(sigma2))
bb=c(2,0,-1,0,2)/3
ww=X%*%bb
Y=cos(ww)+rnorm(100)


##### define the function of DP to get the weights Gj
weights2=function(K,L,V,Z){
  p = rep(NA, L)
  p[1] = V[1]
  for (i in 2:L) {
    p[i] = prod(1 - V[1:(i - 1)]) * V[i]
  }
  p = c(1 - sum(p), p)  
  
  Zbound = matrix(0,K,2)
  for (i in 1:K) {
    Zbound[i,] = c((i-1)/K,i/K)
  }
  Zg = rep(0,L+1)
  # weights of each beta density
  G = rep(0, K) 
  for (j in 1:L+1){
    for (i in 1:nrow(Zbound)) {
      if ((Z[j]>Zbound[i,1]) & (Z[j]<=Zbound[i,2])){
        G[i]=G[i]+p[j]
        Zg[j]=i
        break
      }
    }
  }
  out=list(G=G,Zg=Zg,P=p)
  return(out)
}

##### data structure
n=nrow(X)
p=ncol(X)
S=10000
L=max(20,n^(1/3))
Kmax=50
M=1
a=0.01
b=0.01
K=10
d=3
nu0=500
s0=0.05
lambda=10
dk=4
sigma=0.5
baseb=1.15

##### The initialization

# Initialise G
G=rep(NA,K);
# Initialise V,Z
V = rep(1/2,L)
Z = runif(L + 1)
U = rep(1/2,L)
C = runif(L + 1)

#obtain G(j)
DP=weights2(K,L,V,Z)
G=DP$G

#Initialise knots
DP2=weights2(K-d,L,U,C)
knot_diff=c(0,DP2$G)
in_knots=rep(0,(K-d+1))
in_knots=cumsum(knot_diff)

#Initialise r
r=rbinom(p,1,1/2)
if(all(r==0)){r[sample(1:p,1)]=1}


#Initialise beta
mcmc_r_beta_ratio=rexp(1,3)
mcmc_beta_r_rho=0.2
beta=rep(0,p)
beta[which(r==1)]=rnorm(length(which(r==1)),0,1) #unit norm
beta=beta/c(sqrt(crossprod(beta)));
if(beta[which(r==1)[1]]<0){beta=-beta} #beta1>0

#Initialise kc
kc=runif(n,0,3) 

##### Define Parameter space
KK=Lambda=Bspline_dens=Raa=Sigma=NULL
R=Beta=matrix(nr=S,nc=p)
VV=matrix(nr=S,nc=L)
ZZ=matrix(nr=S,nc=L+1)
UU=matrix(nr=S,nc=L)
CC=matrix(nr=S,nc=L+1)
KC=matrix(nr=S,nc=n)
GG=matrix(nr=S,nc=Kmax)
Rho=R_beta_ratio=NULL

##### the updating process

for(s in 1:S){
  ### the K
  Lk=max(4,round(K-dk));Uk=min(Kmax,ceiling(K+dk));K.s=sample(Lk:Uk,1);
  if(K.s!=K){
    G.s=rep(0,K.s)
    DP.s=weights2(K.s,L,V,Z)
    G.s= DP.s$G
    ### knots
    DP2.s=weights2(K.s-d,L,U,C)
   knot_diff.s=c(0,DP2.s$G)
   in_knots.s=rep(0,(K.s-d+1))
   in_knots.s=cumsum(knot_diff.s)
   Bspline_dens.s=matrix(nr=K.s,nc=n)
   Bspline_dens=matrix(nr=K,nc=n)
   XB=X%*%beta
   s1=s2=0
   bspline_dens=bspline_dens.s=NULL
    for(i in 1:n){
     w=exp(XB[i])/(1+exp(XB[i]))
     bspline_dens=dbspline(w,in_knots)
     Bspline_dens[,i]=bspline_dens
     bspline_dens.s=dbspline(w,in_knots.s)
     Bspline_dens.s[,i]=bspline_dens.s
     s1=s1+dnorm(Y[i],logb(t(G.s)%*%bspline_dens.s,baseb),sqrt(8*kc[i]*sigma),log=TRUE)
     s2=s2+dnorm(Y[i],logb(t(G)%*%bspline_dens,baseb),sqrt(8*kc[i]*sigma),log=TRUE)
    }
    lnr=s1-s2+(K.s-K)*log(lambda)+log(factorial(K))-log(factorial(K.s))
    if(min(exp(lnr),1)>runif(1)){
      K=K.s
      G=NULL
      G=G.s
      DP=DP2=NULL
      DP=DP.s
      DP2=DP2.s
      in_knots=NULL
      in_knots=in_knots.s
      Bspline_dens=matrix(nr=K.s,nc=n)
      Bspline_dens=Bspline_dens.s}
    if(exp(lnr)>0.65) dk=min(4,dk*2)
    else if(exp(lnr)<0.15) dk=max(2,dk/2)
  } ### K.s!=K
  KK[s]=K
  
  ### the lambda
  lambda=rgamma(1,a+K,b+1)
  Lambda[s]=lambda
  
  
  ### the r
  for(j in sample(1:p,p)){
    q=rbeta(1,1+r[j],2-r[j]) 
    rj.s=rbinom(1,1,q) #0,1
    r.s=NULL;r.s=r
    r.s[j]=rj.s
    #### if Beta_dens is NULL
    if(length(bspline_dens)==0){ 
      Bspline_dens=matrix(nr=K,nc=n)
      XB=X%*%beta
      for(i in 1:n){
        w=exp(XB[i])/(1+exp(XB[i]))
        bspline_dens=dbspline(w,in_knots)
        Bspline_dens[,i]=bspline_dens
      }
    }
    
    if(rj.s!=r[j] & sum(r.s)>=1){
      q.s=rbeta(1,1+rj.s,2-rj.s)
      beta.s=NULL
      if(sum(r.s)==1){
        beta.s=rep(0,p)
        beta.s[which(r.s==1)]=1
      }
      else{
        beta.s=beta
        wns=which(Beta[1:(s-1),j]!=0)
        if(length(wns)<=1){
          beta.s[j]=rj.s*rnorm(1,0,1)
        }
        else{
          beta.s[j]=rj.s*rnorm(1,mean(Beta[wns,j]),sd(Beta[wns,j]))
        }
        beta.s=beta.s/c(sqrt(crossprod(beta.s)));
        if(beta.s[which(r.s==1)[1]]<0){beta.s=-beta.s};
      }### sum(r.s)>1
      nr=sum(r);nr.s=sum(r.s)
      s1=s2=0;Bspline_dens.s=matrix(nr=K,nc=n)
      XB.s=X%*%beta.s
      bspline_dens.s=NULL
      for(i in 1:n){
 		    w.s=exp(XB.s[i])/(1+exp(XB.s[i]))
        bspline_dens.s=dbspline(w.s,in_knots)
        Bspline_dens.s[,i]=bspline_dens.s
 			  s1=s1+dnorm(Y[i],logb(t(G)%*%bspline_dens.s,baseb),sqrt(8*kc[i]*sigma),log=TRUE)
 			  s2=s2+dnorm(Y[i],logb(t(G)%*%Bspline_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
 		}   
      lnr=s1-s2+(rj.s)*log(q.s/(1-q.s))-(r[j])*log(q/(1-q))+log(1-q.s)-log(1-q)+(nr-nr.s)*log(pi)/2+log(gamma(nr.s/2))-log(gamma(nr/2))
      if(min(exp(lnr),1)>runif(1)){
        r[j]=rj.s
        beta=NULL;beta=beta.s
        Bspline_dens=matrix(nr=K,nc=n);Bspline_dens=Bspline_dens.s
      }
    }### rj.s!=r[j] & sum(r.s)>=1
  } ### j ends
  R[s,]=r
  
  
  ### the beta
  wb=which(r==1)
  if(length(wb)>1){
    beta.s=rep(0,p)
    mcmc_beta_r_rho=max(0,mcmc_beta_r_rho+(0.234-mcmc_r_beta_ratio)/sqrt(s))
    Rho[s]=mcmc_beta_r_rho
    beta.s[wb]=as.vector(rmvnorm(1,sqrt(2)*mcmc_beta_r_rho*beta[wb],diag(rep(length(wb)))))
    beta.s=beta.s/sqrt(c(crossprod(beta.s[wb])))
    if(beta.s[which(beta.s!=0)[1]]<0){beta.s=-beta.s}
    Bspline_dens.s=matrix(nr=K,nc=n)
    XB.s=X%*%beta.s
    s1=s2=0
    bspline_dens.s=NULL
 		for(i in 1:n){
 		  w.s=exp(XB.s[i])/(1+exp(XB.s[i]))
      bspline_dens.s=dbspline(w.s,in_knots)
      Bspline_dens.s[,i]=bspline_dens.s
 			s1=s1+dnorm(Y[i],logb(t(G)%*%bspline_dens.s,baseb),sqrt(8*kc[i]*sigma),log=TRUE)
 			s2=s2+dnorm(Y[i],logb(t(G)%*%Bspline_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
 		}   
    lnr=s1-s2+(t(beta.s-sqrt(2)*mcmc_beta_r_rho*beta)%*%diag(rep(p))%*%(beta.s-sqrt(2)*mcmc_beta_r_rho*beta)-t(beta-sqrt(2)*mcmc_beta_r_rho*beta.s)%*%diag(rep(p))%*%(beta-sqrt(2)*mcmc_beta_r_rho*beta.s))/2
    mcmc_r_beta_ratio=min(exp(lnr),1);R_beta_ratio[s]=mcmc_r_beta_ratio
    if(mcmc_r_beta_ratio>runif(1)){
      beta=NULL;beta=beta.s;
      Bspline_dens=matrix(nr=K,nc=n);Bspline_dens=Bspline_dens.s
    }
  }##### length(wb)>1
  Beta[s,]=beta
  
  
  ### the V
  for(l in sample(1:L,L)){
    dv=l/(l+2*sqrt(n));Lv=max(0,V[l]-dv);Uv=min(1,V[l]+dv)
    Vl.s=runif(1,Lv,Uv);V.s=V;V.s[l]=Vl.s;G.s=NULL;
    DP.s=weights2(K,L,V.s,Z)
    G.s= DP.s$G
    s1=s2=0
    for(i in 1:n){
      s1=s1+dnorm(Y[i],logb(t(G.s)%*%Bspline_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
      s2=s2+dnorm(Y[i],logb(t(G)%*%Bspline_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
    }
    lnr=s1-s2+(M-1)*(log(1-Vl.s)-log(1-V[l]))
    if(min(exp(lnr),1)>runif(1)){V[l]=Vl.s;G=NULL;G=G.s;DP=NULL;DP=DP.s}
  }## l ends
  VV[s,]=V
  
  ### the Z
  for(l in sample(1:(L+1),L+1)){
    dz=l/(l+2*sqrt(n));Lz=max(0,Z[l]-dz);Uz=min(1,Z[l]+dz)
    Zl.s=runif(1,Lz,Uz);G.s=NULL;G.s=G;Zg=DP$Zg;Zg.s=Zg
    for (i in 1:K) {
      if ((Zl.s>(i-1)/K) & (Zl.s<=i/K)){
        Zg.s[l]=i
        break
      }
    }
    if(Zg.s[l]!=Zg[l]){
      P=DP$P
      G.s[Zg[l]]=G.s[Zg[l]]-P[l]
      G.s[Zg.s[l]]=G.s[Zg.s[l]]+P[l]
      s1=s2=0
      for (i in 1:n){
        s1=s1+dnorm(Y[i],logb(t(G.s)%*%Bspline_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
        s2=s2+dnorm(Y[i],logb(t(G)%*%Bspline_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
      }
      lnr=s1-s2
      if(min(exp(lnr),1)>runif(1)){Z[l]=Zl.s;G=G.s;DP=NULL;DP=list(G=G,Zg=Zg.s,P=P)}
    }
  }
  ZZ[s,]=Z;GG[s,1:K]=G
  
   ### the U
 for(l in sample(1:L,L)){
 	du=l/(l+2*sqrt(n));Lu=max(0,U[l]-du);Uu=min(1,U[l]+du)
 	Ul.s=runif(1,Lu,Uu);U.s=U;U.s[l]=Ul.s;
  ### knots.s
  DP2.s=weights2(K-d,L,U.s,C)
   knot_diff.s=c(0,DP2.s$G)
   in_knots.s=rep(0,(K-d+1))
   in_knots.s=cumsum(knot_diff.s)
   Bspline_dens.s=matrix(nr=K,nc=n)
   XB=X%*%beta
   s1=s2=0
   bspline_dens=bspline_dens.s=NULL
   for(i in 1:n){
     w=exp(XB[i])/(1+exp(XB[i]))
     bspline_dens.s=dbspline(w,in_knots.s)
     Bspline_dens.s[,i]=bspline_dens.s
     s1=s1+dnorm(Y[i],logb(t(G)%*%bspline_dens.s,baseb),sqrt(8*kc[i]*sigma),log=TRUE)
     s2=s2+dnorm(Y[i],logb(t(G)%*%Bspline_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
    }
 	lnr=s1-s2+(M-1)*(log(1-Ul.s)-log(1-U[l]))
 	if(min(exp(lnr),1)>runif(1)){
     U[l]=Ul.s
     in_knots=NULL
     in_knots=in_knots.s
     DP2=NULL;DP2=DP2.s
     Bspline_dens=matrix(nr=K,nc=n);Bspline_dens=Bspline_dens.s}
 }## l ends
 UU[s,]=U

 ### the C
 for(l in sample(1:(L+1),L+1)){
    dc=l/(l+2*sqrt(n));Lc=max(0,C[l]-dc);Uc=min(1,C[l]+dc)
    Cl.s=runif(1,Lc,Uc);in_knots.s=NULL;in_knots.s=in_knots
    Cg=DP2$Zg;Cg.s=Cg
    for (i in 1:(K-d)) {
      if ((Cl.s>(i-1)/(K-d)) & (Cl.s<=i/(K-d))){
        Cg.s[l]=i
        break
      }
    }
    if(Cg.s[l]!=Cg[l]){
      P2=DP2$P
      G2.s=DP2$G
      G2.s[Cg[l]]=G2.s[Cg[l]]-P2[l]
      G2.s[Cg.s[l]]=G2.s[Cg.s[l]]+P2[l]
    ### knots.s
    knot_diff.s=c(0,G2.s)
    in_knots.s=rep(0,(K-d+1))
    in_knots.s=cumsum(knot_diff.s)
    Bspline_dens.s=matrix(nr=K,nc=n)
    XB=X%*%beta
    s1=s2=0
    bspline_dens=bspline_dens.s=NULL

    for(i in 1:n){
     w=exp(XB[i])/(1+exp(XB[i]))
     bspline_dens.s=dbspline(w,in_knots.s)
     Bspline_dens.s[,i]=bspline_dens.s
     s1=s1+dnorm(Y[i],logb(t(G)%*%bspline_dens.s,baseb),sqrt(8*kc[i]*sigma),log=TRUE)
     s2=s2+dnorm(Y[i],logb(t(G)%*%Bspline_dens[,i],baseb),sqrt(8*kc[i]*sigma),log=TRUE)
    }
 		lnr=s1-s2
    if(min(exp(lnr),1)>runif(1)){
      C[l]=Cl.s
      in_knots=NULL;in_knots=in_knots.s
      DP2=NULL;DP2=list(G=G2.s,Zg=Cg.s,P=P2)
      Bspline_dens=matrix(nr=K,nc=n);Bspline_dens=Bspline_dens.s}
    }
  }
   CC[s,]=C

  ### the sigma
  mu=t(Y)-logb(t(G)%*%Bspline_dens,baseb)
  sy=mu%*%diag(1/kc)%*%t(mu)
  sigma=rinvgamma(1,(nu0+3*n)/2,(nu0*s0+2*sum(kc)+sy/8)/2)
  Sigma[s]=sigma
  
  ### the kc
  eta2=2/sigma
  for (i in 1:n){
    eta1=((Y[i]-logb(t(G)%*%Bspline_dens[,i],baseb))^2)/(8*sigma)
    kc[i]=rgig(1,1/2,sqrt(eta1),sqrt(eta2))
  } 
  KC[s,]=kc
  
  print(paste(s,Sys.time()," "))
}

###the posterior estimate
R.pred=Beta.pred=G.pred=NULL

### the K.pred
### the function to get mode
getmode = function(v) {
  uniqv = unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
K.pred=getmode(KK);
### the R.pred
wke=which(KK==K.pred)
lwke=length(wke)
R0=R[wke,]
Beta0=Beta[wke,]
R.uniqv=unique(R0)
lR=nrow(R.uniqv)
R.count=rep(0,lR)
R.cg=rep(0,lwke)
for (i in 1: lwke){
  for (j in 1: lR){
    if (isTRUE(all.equal(R0[i,],R.uniqv[j,]))){
      R.count[j]=R.count[j]+1
      R.cg[i]=j
    }
  }
  
} 
R.pred=R.uniqv[which.max(R.count),]

### the beta.pred
beta0.order=which(R.cg==which.max(R.count))
Beta00=matrix(nr=max(R.count),nc=p)
Beta00=Beta0[beta0.order,]
beta=apply(Beta00,2,mean,na.rm=TRUE)
beta=beta/c(sqrt(crossprod(beta)))
if(beta[which(R.pred==1)[1]]<0){beta=-beta}
Beta.est=beta;


Beta.est

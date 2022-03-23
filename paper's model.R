Postsample=function(X,Y,S=10000,L=max(20,n^(1/3)),Kmax=50,t=0.5,M=1,nu0=500,s0=0.05,a=0.01,b=0.01,K=8,lambda=10,dk=2,sigma=0.5){
  
  n=nrow(X);p=ncol(X) ### n: sample size; p: dimension of predictor
  
  ########### Initialise
  rz=P=G=NULL;
  V=rep(1/2,L);Z=runif(L+1); # Initialise V,Z
  h=1;L.s=1:(L+1);L.ex=NULL;
  while(h<=K){
    for(j in L.s){
      if(Z[j]>(h-1)/K & Z[j]<=h/K){rz[j]=h;L.ex=c(L.ex,j)}
    }
    if(length(L.ex)>0){L.s=(1:(L+1))[-L.ex]};if(length(L.s)<1){break}
    h=h+1
  }
  P[2]=V[1];l=3;
  while(l<=(L+1)){P[l]=0.5^(l-1);l=l+1};P[1]=1-sum(P[-1]); #Initialise P
  h=1;while(h<=K){G[h]=sum(P[which(rz==h)]);h=h+1}; #Initialise G
  
  r=rbinom(p,1,1/2);if(all(r==0)){r[sample(1:p,1)]=1}; #Initialise r
  r_beta=rexp(1,3);rho=0.2; #Initialise beta
  beta=rep(0,p);beta[which(r==1)]=rnorm(length(which(r==1)),0,1)
  beta=beta/c(sqrt(crossprod(beta)));
  if(beta[which(r==1)[1]]<0){beta=-beta}
  kc=runif(n,0,3) #Initialise \xi
  
  tt=t*(1-t)
  
  ########## Define Parameter space
  KK=Lambda=beta_dens=Raa=Sigma=NULL
  R=Beta=matrix(nr=S,nc=p);VV=matrix(nr=S,nc=L);ZZ=matrix(nr=S,nc=L+1);KC=matrix(nr=S,nc=n);GG=matrix(nr=S,nc=Kmax);Rho=R_beta=NULL;
  
  
  
  #################### S iterations for posterior update steps
  ########   '*.s': proposal candidate of parameter '*'
  s=1;
  while(s<=S){
    
    ###### Update K
    Lk=max(3,round(K-dk));Uk=min(Kmax,ceiling(K+dk));K.s=sample(Lk:Uk,1);
    if(K.s!=K){
      G.s=rep(0,K.s);rz.s=NULL;
      l=1;
      while(l<=(L+1)){
        j=1;
        while(j<=K.s){
          if(Z[l]>(j-1)/K.s & Z[l]<=j/K.s){
            rz.s[l]=j;G.s[j]=G.s[j]+P[l];break
          }
          j=j+1
        }
        l=l+1
      }
      
      Beta_dens.s=matrix(nr=K.s,nc=n);Beta_dens=matrix(nr=K,nc=n)
      s1=s2=0;i=1;
      while(i<=n){
        beta_dens=beta_dens.s=NULL;j=1;w=X[i,]%*%beta;
        while(j<=K.s){
          beta_dens.s[j]=(exp(w)/(1+exp(w)))^(j-1)*(1/(1+exp(w)))^(K.s-j)/beta(j,K.s-j+1)
          j=j+1
        }
        Beta_dens.s[,i]=beta_dens.s;j=1
        while(j<=K){
          beta_dens[j]=(exp(w)/(1+exp(w)))^(j-1)*(1/(1+exp(w)))^(K-j)/beta(j,K-j+1)
          j=j+1
        }
        Beta_dens[,i]=beta_dens
        s1=s1+dnorm(Y[i],t(G.s)%*%beta_dens.s+(1-2*t)*kc[i]/tt,sqrt(2*kc[i]*sigma/tt),log=TRUE)
        s2=s2+dnorm(Y[i],t(G)%*%beta_dens+(1-2*t)*kc[i]/tt,sqrt(2*kc[i]*sigma/tt),log=TRUE)
        i=i+1
      }
      lnr=s1-s2+(K.s-K)*log(lambda)+log(factorial(K))-log(factorial(K.s));
      if(min(exp(lnr),1)>runif(1)){K=K.s;G=rz=NULL;G=G.s;rz=rz.s;Beta_dens=matrix(nr=K.s,nc=n);Beta_dens=Beta_dens.s}
      if(exp(lnr)>0.65) dk=min(4,dk*2) else if(exp(lnr)<0.15) dk=max(2,dk/2)
    } ### K.s!=K
    KK[s]=K;##GG[s,1:K]=G
    
    ######## Update lambda
    lambda=rgamma(1,a+K,b+1)
    Lambda[s]=lambda
    
    ######## Update r
    for(j in sample(1:p,p)){
      pp=rbeta(1,1+r[j],2-r[j]);rj.s=rbinom(1,1,pp);r.s=NULL;r.s=r;r.s[j]=rj.s
      
      if(length(beta_dens)==0){ ####Beta_dens is NULL
        Beta_dens=matrix(nr=K,nc=n);i=1
        while(i<=n){
          w=X[i,]%*%beta;
          h=1;while(h<=K){
            beta_dens[h]=(exp(w)/(1+exp(w)))^(h-1)*(1/(1+exp(w)))^(K-h)/beta(h,K-h+1)
            h=h+1
          }
          Beta_dens[,i]=beta_dens;i=i+1
        }###i
      }##### Beta_dens is NULL
      
      if(rj.s!=r[j] & sum(r.s)>=1){
        p.s=rbeta(1,1+rj.s,2-rj.s);beta.s=NULL;
        if(sum(r.s)==1){beta.s=rep(0,p);beta.s[which(r.s==1)]=1}else{
          beta.s=beta;wns=which(Beta[1:(s-1),j]!=0);
          if(length(wns)<=1){beta.s[j]=rj.s*rnorm(1,0,1)}
          else{beta.s[j]=rj.s*rnorm(1,mean(Beta[wns,j]),sd(Beta[wns,j]))}
          beta.s=beta.s/c(sqrt(crossprod(beta.s)));
          if(beta.s[which(r.s==1)[1]]<0){beta.s=-beta.s};
          
        }### sum(r.s)>1
        nr=sum(r);nr.s=sum(r.s);s1=s2=0;i=1;Beta_dens.s=matrix(nr=K,nc=n)
        while(i<=n){
          w.s=X[i,]%*%beta.s;beta_dens.s=NULL;h=1
          while(h<=K){
            beta_dens.s[h]=(exp(w.s)/(1+exp(w.s)))^(h-1)*(1/(1+exp(w.s)))^(K-h)/beta(h,K-h+1)
            h=h+1
          }
          s1=s1+dnorm(Y[i],t(G)%*%beta_dens.s+(1-2*t)*kc[i]/tt,sqrt(2*kc[i]*sigma/tt),log=TRUE)
          s2=s2+dnorm(Y[i],t(G)%*%Beta_dens[,i]+(1-2*t)*kc[i]/tt,sqrt(2*kc[i]*sigma/tt),log=TRUE)
          Beta_dens.s[,i]=beta_dens.s;i=i+1
        }   
        lnr=s1-s2+r[j]*log(p.s/pp)+(1-r[j])*log((1-p.s)/(1-pp))+(nr-nr.s)*log(pi)/2+log(gamma(nr.s/2))-log(gamma(nr/2))
        if(min(exp(lnr),1)>runif(1)){r[j]=rj.s;beta=NULL;beta=beta.s;Beta_dens=matrix(nr=K,nc=n);Beta_dens=Beta_dens.s}
      }### rj.s!=r[j] & sum(r.s)>=1
    }### j ends
    R[s,]=r
    
    ######### Update beta
    wb=which(r==1)
    if(length(wb)>1){
      beta.s=rep(0,p);rho=max(0,rho+(0.234-r_beta)/sqrt(s));Rho[s]=rho
      beta.s[wb]=as.vector(rmvnorm(1,sqrt(2)*rho*beta[wb],diag(rep(length(wb)))))
      beta.s=beta.s/sqrt(as.numeric(crossprod(beta.s[wb])));
      if(beta.s[which(beta.s!=0)[1]]<0){beta.s=-beta.s}
      Beta_dens.s=matrix(nr=K,nc=n);s1=s2=0;i=1
      while(i<=n){
        w.s=X[i,]%*%beta.s;beta_dens.s=NULL;j=1
        while(j<=K){beta_dens.s[j]=dbeta(exp(w.s)/(1+exp(w.s)),j,K-j+1);j=j+1}
        Beta_dens.s[,i]=beta_dens.s
        s1=s1+dnorm(Y[i],t(G)%*%beta_dens.s+(1-2*t)*kc[i]/tt,sqrt(2*kc[i]*sigma/tt),log=TRUE)
        s2=s2+dnorm(Y[i],t(G)%*%Beta_dens[,i]+(1-2*t)*kc[i]/tt,sqrt(2*kc[i]*sigma/tt),log=TRUE)
        i=i+1
      }
      lnr=s1-s2+(t(beta.s-sqrt(2)*rho*beta)%*%solve(1*diag(rep(p)))%*%(beta.s-sqrt(2)*rho*beta)-t(beta-sqrt(2)*rho*beta.s)%*%solve(1*diag(rep(p)))%*%(beta-sqrt(2)*rho*beta.s))/2
      r_beta=exp(lnr);
      if(min(exp(lnr),1)>runif(1)){beta=NULL;beta=beta.s;Beta_dens=matrix(nr=K,nc=n);Beta_dens=Beta_dens.s}
    }##### length(wb)>1
    Beta[s,]=beta
    
    ############ Update V
    for(l in sample(1:L,L)){
      dv=l/(l+2*sqrt(n));Lv=max(0,V[l]-dv);Uv=min(1,V[l]+dv)
      Vl.s=runif(1,Lv,Uv);V.s=V;V.s[l]=Vl.s;P.s=P;G.s=NULL;
      if(l==1){P.s[l+1]=Vl.s}else{
        P.s[l+1]=P.s[l+1]/V[l]*Vl.s
      }#l>1
      if(l<L){
        j=l+2
        while(j<=(L+1)){
          P.s[j]=P.s[j-1]/V.s[j-2]*(1-V.s[j-2])*V.s[j-1];j=j+1
        }
      }
      P.s[1]=1-sum(P.s[-1])
      j=1;while(j<=K){G.s[j]=sum(P.s[which(rz==j)]);j=j+1}
      s1=s2=0;i=1
      while(i<=n){
        s1=s1+dnorm(Y[i],t(G.s)%*%Beta_dens[,i]+(1-2*t)*kc[i]/tt,sqrt(2*kc[i]*sigma/tt),log=TRUE)
        s2=s2+dnorm(Y[i],t(G)%*%Beta_dens[,i]+(1-2*t)*kc[i]/tt,sqrt(2*kc[i]*sigma/tt),log=TRUE)
        i=i+1
      }
      lnr=s1-s2+(M-1)*(log(1-Vl.s)-log(1-V[l]))
      if(min(exp(lnr),1)>runif(1)){V[l]=Vl.s;P=G=NULL;P=P.s;G=G.s}
    }## l ends
    VV[s,]=V
    
    ###################  Update Z
    for(l in sample(1:(L+1),L+1)){
      dz=l/(l+2*sqrt(n));Lz=max(0,Z[l]-dz);Uz=min(1,Z[l]+dz)
      Zl.s=runif(1,Lz,Uz);G.s=NULL;G.s=G;h=1
      while(h<=K){
        if(Zl.s>(h-1)/K & Zl.s<=h/K){j.s=h;break}
        h=h+1
      }
      if(j.s!=rz[l]){
        G.s[rz[l]]=G.s[rz[l]]-P[l];G.s[j.s]=G.s[j.s]+P[l]
        s1=s2=0;i=1
        while(i<=n){
          s1=s1+dnorm(Y[i],t(G.s)%*%Beta_dens[,i]+(1-2*t)*kc[i]/tt,sqrt(2*kc[i]*sigma/tt),log=TRUE)
          s2=s2+dnorm(Y[i],t(G)%*%Beta_dens[,i]+(1-2*t)*kc[i]/tt,sqrt(2*kc[i]*sigma/tt),log=TRUE)
          i=i+1
        }
        lnr=s1-s2
        if(min(exp(lnr),1)>runif(1)){Z[l]=Zl.s;G[rz[l]]=G.s[rz[l]];G[j.s]=G.s[j.s];rz[l]=j.s}
      }## j.s!=rz[l]
    }## l ends
    ZZ[s,]=Z;GG[s,1:K]=G
    
    ############# Update sigma
    sy=(t(Y)-t(G)%*%Beta_dens-(1-2*t)*kc/tt)%*%diag(1/kc)%*%t(t(Y)-t(G)%*%Beta_dens-(1-2*t)*kc/tt)
    sigma=rinvgamma(1,(nu0+3*n)/2,(nu0*s0+2*sum(kc)+tt*sy/2)/2)
    Sigma[s]=sigma
    
    ############# Update kc
    eta=(2+(1-2*t)^2/(2*tt))/sigma;i=1
    while(i<=n){
      lmdd=(tt*(Y[i]-t(G)%*%Beta_dens[,i])^2)/(2*sigma)
      kc[i]=rgig(1,1/2,lmdd,eta);
      i=i+1
    } 
    KC[s,]=kc
    
    print(paste(s,Sys.time()," "))
    
    s=s+1
  }################ S ends
  
  
  ############### Obtain posterior estimate
  #### '*.est': estimate of parameter '*'
  R.est=Beta.est=G.est=NULL
  k=1;count=NULL;while(k<=Kmax){count[k]=length(which(KK==k));k=k+1};K1=which.max(count);K.est=K1;
  
  wke=which(KK==K1);lwke=length(wke);
  K=K1;lp=2^p;R00=matrix(nr=min(lp,lwke),nc=p);R0=R[wke,];Beta0=Beta[wke,];lR=dim(R0)[1];R00[1,]=R0[1,];count=rep(0,nrow(R00));
  count[1]=1;record_r0=1;record_r00=matrix(nr=min(lp,lwke),nc=lwke);record_r00[1,1]=1;l=2
  while(l<=lR){
    j=1
    while(j<=record_r0){
      if(all(R0[l,]==R00[j,])){count[j]=count[j]+1;record_r00[j,count[j]]=l;break}
      else{j=j+1}
    }
    if(j==(record_r0+1)){record_r0=record_r0+1;R00[record_r0,]=R0[l,];count[record_r0]=1;record_r00[record_r0,1]=l}
    l=l+1
  }
  mc=which.max(count);r00=R00[mc,];r=r00;R.est=r
  
  Beta00=matrix(nr=max(count)-1,nc=p);Beta00=Beta0[record_r00[mc,1:(max(count)-1)],];
  beta=apply(Beta00,2,mean,na.rm=TRUE);beta=beta/c(sqrt(crossprod(beta)));if(beta[which(r00==1)[1]]<0){beta=-beta}
  Beta.est=beta;
  
  if(K1==1){G.est1=1}else{G.est=apply(GG[wke,1:K1],2,mean,na.rm=TRUE)}
  kc=apply(KC[wke[floor(lwke/3):lwke],],2,mean)
  sigma.est=mean(Sigma[wke[floor(lwke/3):lwke]]);
  Lambda.est=mean(Lambda)
  
  RL=list(K.est=K.est,R.est=R.est,Beta.est=Beta.est,G.est=G.est,sigma.est=sigma.est)
  
  RL
  
  
}

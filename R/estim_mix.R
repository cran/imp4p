
####################################
#
# Estimating the mixture model
#
####################################


estim.mix=function(tab, tab.imp, conditions, x.step.mod=100, x.step.pi=100, nb.rei=50, method=1, gridsize=100){

  #require(Iso);
  min.x=min(tab.imp,na.rm=TRUE);
  max.x=max(tab.imp,na.rm=TRUE);

  x.min=min.x-(max.x-min.x)/2;
  x.max=max.x+(max.x-min.x)/2;

  #Function to minimize : Weibull assumption
  fr <- function(pi,x,pi_est,w,Ft) {
    p <- pi[1];
    b <- pi[2];
    c <- pi[3];
    f=p+(1-p)*exp(-b*(x)^c)/(1-Ft);
    return(sum(w*(pi_est-f)^2))
  }

  #Number of conditions
  nb_cond=length(levels(conditions));

  nb_rep=rep(0,nb_cond);
  pi_abs=rep(0,nb_cond);
  pi_na=rep(0,length(tab[1,]));
  pi_mcar=rep(0,length(tab[1,]));
  sta=rep(0,length(tab[1,]));
  dat_nbNA=matrix(0,length(tab[,1]),nb_cond);
  abs=matrix(0,x.step.pi,length(tab[1,]));
  absi=abs;
  PI_INIT=matrix(0,x.step.pi,length(tab[1,]));
  VAR_PI_INIT=matrix(0,x.step.pi,length(tab[1,]));
  pi.trend=matrix(0,x.step.pi,length(tab[1,]));
  pi.trend.pava=matrix(0,x.step.pi,length(tab[1,]));
  pi.init.pava=matrix(0,x.step.pi,length(tab[1,]));
  FTOT.i=matrix(0,x.step.pi,length(tab[1,]));
  FNA.i=matrix(0,x.step.pi,length(tab[1,]));

  q.n=matrix(0,100,length(tab[1,]));
  q.c=matrix(0,100,length(tab[1,]));
  s.q.n=rep(0,length(tab[1,]));
  s.q.c=rep(0,length(tab[1,]));

  ABSC=matrix(0,x.step.mod,length(tab[1,]));
  FTOT=matrix(0,x.step.mod,length(tab[1,]));
  FNA=matrix(0,x.step.mod,length(tab[1,]));
  FOBS=matrix(0,x.step.mod,length(tab[1,]));
  FMNAR=matrix(0,x.step.mod,length(tab[1,]));

  k=1;
  j=1;
  for (i in 1:nb_cond){
    #Number of replicates in each condition
    nb_rep[i]=sum((conditions==levels(conditions)[i]));
    #Number of missing values on each row in the condition
    dat_nbNA[,i]=apply(tab[,(k:(k+nb_rep[i]-1))],1, function(x) sum(is.na(x)));
    #Percentage of rows without observed values in the condition
    pi_abs[i]=sum(dat_nbNA[,i]==nb_rep[i])/sum(dat_nbNA[,i]>0);
    #Deletion of rows without data in the condition
    liste_retire=which(dat_nbNA[,i]==nb_rep[i]);
    liste_garde=which(dat_nbNA[,i]!=nb_rep[i]);
    #Data after deletion
    tab2=as.matrix(tab[liste_garde,])
    ############
    #Roughly imputed data after deletion
    tab_imp2=as.matrix(tab.imp[liste_garde,]);

    ############

    while (j<(k+nb_rep[i])){
      #Percentage of missing values in sample j
      nna=sum(is.na(tab2[,j]));
      pi_na[j]=nna/length(tab2[,j]);

      ############
      #Rough estimation of the distribution of missing values in sample j
      v_na=tab_imp2[is.na(tab2[,j]),j];
      F_na=ecdf(v_na);

      ############
      #Distribution of observed values in sample j
      F_obs=ecdf(tab2[which(!is.na(tab2[,j])),j]);

      ############
      #Determination of the interval on which is estimated pi^MLE(x)
      #F_na(x)>0 et F_obs(x)>0
      xmin=max(sort(v_na,na.last=T)[1],sort(tab2[which(!is.na(tab2[,j])),j],na.last=T)[1]);
      #F_na(x)<1 et F_obs(x)<1
      xmax=min(-sort(-v_na,na.last=T)[2],-sort(-tab2[which(!is.na(tab2[,j])),j],na.last=T)[2]);
      absi[,j]=seq(xmin,xmax,length.out=x.step.pi)
      pi_init=NULL;
      for (x in absi[,j]){
        #Initial estimator of the proportion (1-F_na(x))/(1-F(x))
        if (((pi_na[j]*F_na(x)+(1-pi_na[j])*F_obs(x))!=1)){
            pi_init_x=(1-F_na(x))/(1-(pi_na[j]*F_na(x)+(1-pi_na[j])*F_obs(x)));
        }else{
            pi_init_x=(1-F_na(x))/(1-(1-1e-10));
        }
        pi_init=c(pi_init,pi_init_x);
      }
      kx=1
      while (pi_init[kx]>=1){
        kx=kx+1;
      }
      xmin=absi[kx,j];
      kx=1;
      while ((kx<length(pi_init))&(pi_init[kx]>=(max(pi_init,na.rm=TRUE)-min(pi_init,na.rm=TRUE))/2)){
        kx=kx+1;
      }
      xmean=absi[kx,j];
      if(xmin>(xmean-(xmax-xmean))){
        xmin=(xmean-(xmax-xmean));
        absi[,j]=seq(xmin,xmax,length.out=x.step.pi)
        pi_init=NULL;
        for (x in absi[,j]){
          #Initial estimator of the proportion (1-F_na(x))/(1-F(x))
          pi_init_x=(1-F_na(x))/(1-(pi_na[j]*F_na(x)+(1-pi_na[j])*F_obs(x)))
          pi_init=c(pi_init,pi_init_x)
        }
      }
      abs[,j]=seq(xmin,xmax,length.out=x.step.pi)

      ############
      #Estimation of pi^MLE(x)
      pi_init=NULL;
      F_tot1=NULL;
      F_na1=NULL;
      Fobs=NULL;

      varasy=NULL;
      del=NULL;
      mm=NULL;
      gg=NULL;
      hh=NULL;
      kk=NULL;
      for (x in abs[,j]){
        ############
        #Initial estimator of the proportion (1-F_na(x))/(1-F(x))
        pi_init_x=(1-F_na(x))/(1-(pi_na[j]*F_na(x)+(1-pi_na[j])*F_obs(x)));
        pi_init=c(pi_init,pi_init_x);
        pi_init[pi_init>=1]=1;
        if (pi_init_x>=1){pi_init_x=1;}
        ############
        #Estimation of the cdfs
        F_na1=c(F_na1,F_na(x))
        Fobs=c(Fobs,F_obs(x))
        F_tot1=c(F_tot1,pi_na[j]*F_na(x)+(1-pi_na[j])*F_obs(x))

        ############
        #Estimation of the asymptotic variance of the estimator
        del=(pi_na[j]-1)*pi_init_x/(pi_na[j]*pi_init_x-1)

        p=(1-F_obs(x))
        if(p==1){p=1-1e-10;}

        mm=(1-del*(1-F_obs(x)))/(pi_init_x*(p*(pi_na[j]-1)-pi_na[j])+1)^2
        a=pi_na[j]

        I22=del*p*(1-del*p)*((1/pi_init_x+p*(1-a)/(a*(p-1)*pi_init_x-p*pi_init_x+1))^2)/(1-a*pi_init_x)^2
        I11=(1/a-1)/(p*(1-p))+del*p*(1-del*p)*((1/p-pi_init_x*(1-a)/(a*(p-1)*pi_init_x-p*pi_init_x+1))^2)
        I12=(1-pi_na[j])*mm
        I=matrix(c(I11,I12,I12,I22),2,2)
        v=solve(I)[2,2]
        varasy=c(varasy,((1-pi_na[j])/pi_na[j])*v)

      }
      varasy=varasy/length(tab2[which(!is.na(tab2[,j])),j])
      varasy[varasy<(0.005^2)]=(0.005^2);
      PI_INIT[,j]=pi_init;
      VAR_PI_INIT[,j]=varasy;

      #Function to minimize
      #The minimization is performed nb.rei times to try to find a global minimum
      pim=matrix(0,nb.rei,4)
      for (nbit in 1:nb.rei){
        init=c(runif(1,0,0.5),runif(1,1,50),runif(1,1,10))
        nbtest=1;
        #while (!inherits(try(optim(init, fr, gr=NULL, x=(abs[,j]-xmin) , pi_est=pi_init, lower=rep(0,0,3), upper=c(1,Inf,Inf), method="L-BFGS-B", w=1/varasy, Ft=F_tot1), TRUE), "try-error")==FALSE){
        while ((!inherits(try(optim(init, fr, gr=NULL, x=(abs[,j]-xmin) , pi_est=pi_init, lower=c(0.01,0.0005,2), upper=c(0.5,0.01,5), method="L-BFGS-B", w=1/varasy, Ft=F_tot1), TRUE), "try-error")==FALSE)&(nbtest<nb.rei)){

            init=c(runif(1,0,1),runif(1,1,50),runif(1,1,10));
            nbtest=nbtest+1;
        }
        #re=optim(init, fr, gr=NULL, x=(abs[,j]-xmin) , pi_est=pi_init, lower=rep(0,0,3), upper=c(1,Inf,Inf), method="L-BFGS-B", w=1/varasy, Ft=F_tot1);
        re=optim(init, fr, gr=NULL, x=(abs[,j]-xmin) , pi_est=pi_init, lower=c(0.01,0.0005,2), upper=c(0.5,0.01,5), method="L-BFGS-B", w=1/varasy, Ft=F_tot1);
        pim[nbit,1]=re$par[1]
        pim[nbit,2]=re$par[2]
        pim[nbit,3]=re$par[3]
        pim[nbit,4]=re$value
      }

      #Final estimation of the proportion of MCAR values with different methods
      pi_trend=rep(0,length(abs[,j]))
      pi_trend_pava=rep(0,length(abs[,j]))
      pi_init_pava=rep(0,length(abs[,j]))
      if (method==1){
        pi_mcar[j]=pim[which.min(pim[,4]),1];
      }
      if (method==2){
        pi_trend=pim[which.min(pim[,4]),1]+(1-pim[which.min(pim[,4]),1])*exp(-pim[which.min(pim[,4]),2]*((abs[,j]-xmin))^pim[which.min(pim[,4]),3])/(1-F_tot1)
        pi_mcar[j]=pi_trend[length(pi_trend)];
      }
      if (method==3){
        pi_trend=pim[which.min(pim[,4]),1]+(1-pim[which.min(pim[,4]),1])*exp(-pim[which.min(pim[,4]),2]*((abs[,j]-xmin))^pim[which.min(pim[,4]),3])/(1-F_tot1)
        pi_trend_pava=pava(pi_trend,decreasing=T);
        pi_mcar[j]=pi_trend_pava[length(pi_trend_pava)];
      }
      if (method==4){
        pi_init_pava=pava(pi_init,decreasing=T);
        pi_mcar[j]=pi_init_pava[length(pi_init_pava)];
      }
      if (method==5){
        rr=hist(pi_init,breaks=20,plot=F);
        ri=diff(rr$counts);
        ri=ri<0;
        kh=1;
        while (ri[kh]==F){kh=kh+1;}
        pi_mcar[j]=rr$mids[kh];
      }
      if (method==6){
        Freq=diff(c(0, F_na1))
        gridsize=300
        distance=NULL
        for(kl in 1:gridsize){
          a=kl/gridsize ## Assumes a value of the mixing proportion
          F.hat=(F_na1-(1-a)*F_tot1)/a ## Computes the naive estimator
          F.is=pava(F.hat,Freq,decreasing=FALSE) ## Computes the isotonic estimator
          F.is[which(F.is<=0)]=0
          F.is[which(F.is>=1)]=1
          distance=c(distance,a*sqrt(t((F.hat[!is.nan(F.is)]-F.is[!is.nan(F.is)])^2)%*%Freq[!is.nan(F.is)]));
        }
        dder=diff(distance[2:length(distance)])-diff(distance[1:(length(distance)-1)]);
        pi_mcar[j] <- 1-(which.max(dder)+1)/gridsize;
      }
      pi_mcar[j] = max(min(pi_mcar[j],1-1e-5),1e-5);
      pi.trend[,j]=pi_trend;
      pi.trend.pava[,j]=pi_trend_pava;
      pi.init.pava[,j]=pi_init_pava;

      ############
      #Find the normal distribution of complete values with
      #Regression between sufficiently high quantiles of observed values and normal quantiles
      h=1;
      x=abs[h,j];
      Fmn=0;
      while ((Fmn[length(Fmn)]<0.9)&(h<length(abs[,j]))){
        Fmn=c(Fmn,(1-pi_na[j]*pi_mcar[j])*F_na(x)/(1-pi_mcar[j])-(1-pi_na[j])*pi_mcar[j]*F_obs(x)/(1-pi_mcar[j]));
        h=h+1;
        x=abs[h,j];
      }
      sta[j]=F_obs(x);

      gamma=pi_na[j]*(1-pi_mcar[j])/(1-pi_mcar[j]*pi_na[j]);
      upper.q=0.99
      interv=gamma+(1-gamma)*seq(sta[j],upper.q,length.out=100);

      q.n[,j]=qnorm(gamma+(1-gamma)*seq(0.01,upper.q,length.out=100), mean = 0, sd = 1);
      q.c[,j]=quantile(tab2[,j], probs = seq(0.01,upper.q,length.out=100), na.rm = T);
      s.q.n[j]=qnorm(gamma+(1-gamma)*sta[j], mean = 0, sd = 1);
      s.q.c[j]=quantile(tab2[,j], probs = sta[j], na.rm = T);

      q.normal = qnorm(interv, mean = 0, sd = 1);
      q.curr.sample = quantile(tab2[,j], probs = seq(sta[j],upper.q,length.out=100), na.rm = T);
      temp.QR = lm(q.curr.sample ~ q.normal);
      #mean and variance of the estimated normal density of complete values
      m = temp.QR$coefficients[1];
      v = (as.numeric(temp.QR$coefficients[2]))^2;

      ############
      #Final estimation of the mixture model on the interval precised in input

      ABSC[,j]=seq(x.min,x.max,length.out=x.step.mod);

      F_tot=pnorm(ABSC[,j],mean=m,sd=sqrt(v));
      theta=pi_mcar[j]*(1-pi_na[j])/(1-pi_mcar[j]*pi_na[j]);
      F_mna=(F_tot/pi_na[j]-((1-pi_na[j])/pi_na[j]+theta)*F_obs(ABSC[,j]))/(1-theta);
      F_mnar=pava(F_mna,decreasing=F);
      F_mnar[F_mnar>1]=1;
      F_mnar[F_mnar<0]=0;

      #Reestimation of the cdfs
      kl=1
      F_na2=NULL
      F_tot2=NULL
      Fobs=NULL
      for (x in ABSC[,j]){
        aa=((1-pi_mcar[j])*F_mnar[kl])/(1-pi_mcar[j]*pi_na[j]);
        bb=(pi_mcar[j]*(1-pi_na[j])*F_obs(x))/(1-pi_mcar[j]*pi_na[j]);
        F_na2=c(F_na2,aa+bb);
        Fobs=c(Fobs,F_obs(x));
        F_tot2=c(F_tot2,pi_na[j]*(aa+bb)+(1-pi_na[j])*F_obs(x));
        kl=kl+1;
      }

      FTOT.i[,j]=F_tot1;
      FNA.i[,j]=F_na1;
      FTOT[,j]=F_tot2;
      FNA[,j]=F_na2;
      FOBS[,j]=Fobs;
      FMNAR[,j]=F_mnar;

      j=j+1;
    }

    k=k+nb_rep[i];
  }

  return(list(abs.pi=abs,pi.init=PI_INIT,var.pi.init=VAR_PI_INIT,
              abs.mod=ABSC[,1],pi.na=pi_na,F.na=FNA,F.tot=FTOT,F.obs=FOBS,F.mnar=FMNAR,pi.mcar=pi_mcar,
              pi.m=pim,
              F.na.i=FNA.i,F.tot.i=FTOT.i,pi.trend=pi.trend,pi.trend.pava=pi.trend.pava,pi.init.pava=pi.init.pava,
              eta=sta,q.n=q.n,q.c=q.c,s.q.n=s.q.n,s.q.c=s.q.c
  ));
}

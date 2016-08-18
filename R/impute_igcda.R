
#impute the data to get a gaussian distribution in each column

impute.igcda <-function (tab, tab.imp, conditions, q=0.95){
  
  dataSet.imputed=tab;
  nb_cond=length(levels(conditions));
  nb_rep=rep(0,nb_cond);
  k=1;
  j=1;
  for (i in 1:nb_cond){
    #Number of replicates in each condition
    nb_rep[i]=sum((conditions==levels(conditions)[i]));
    dataSet.mvs=tab[,(k:(k+nb_rep[i]-1))];

    mim=apply(dataSet.mvs,1,min,na.rm=TRUE);
    mam=apply(dataSet.mvs,1,max,na.rm=TRUE);
    rmammim=mam-mim;
    rangem=quantile(rmammim,probs=q);
    mini=min(mim,na.rm=T);
    
    while (j<(k+nb_rep[i])){
      curr.sample=tab[, j];
      
      if (sum(is.na(curr.sample))!=0){
        Fobs=ecdf(curr.sample);
        sta=Fobs(quantile(tab.imp[which(is.na(curr.sample)),j],0.9));
        #
        gamma=sum(is.na(curr.sample))/length(curr.sample);
        interv=gamma+(1-gamma)*seq(sta,0.999,length.out=100);
        q.normal = qnorm(interv, mean = 0, sd = 1);
        q.curr.sample = quantile(curr.sample, probs = seq(sta,0.999,length.out=100), na.rm = T);
        temp.QR = lm(q.curr.sample ~ q.normal);
        #mean and variance of the estimated normal density of complete values
        m = temp.QR$coefficients[1];
        v = (as.numeric(temp.QR$coefficients[2]))^2;
        
        #Computation of missing value density
        dx=seq(max(c(0,mini-rangem)),max(curr.sample,na.rm=T),length.out=1e3);
        Fn=pnorm(dx,mean=m,sd=sqrt(v));
        Fna=(Fn-(1-gamma)*Fobs(dx))/gamma;
        Fna[Fna>1]=1;
        Fna=pava(Fna);
        #approximate probability to get a value in each interval of dx
        pna=diff(Fna);
        gen.sample=NULL;
        for (ll in 1:length(pna)){
          gen.sample=c(gen.sample,runif(floor(pna[ll]*max(1e4,sum(is.na(curr.sample))+1)),dx[ll],dx[ll+1]));
        }
        
        curr.sample.imputed=curr.sample;
        gss=gen.sample[floor(runif(sum(is.na(curr.sample)),1,length(gen.sample)+1))];
        curr.sample.imputed[which(is.na(curr.sample))][order(mim[which(is.na(curr.sample))])]=sort(gss);
        
        dataSet.imputed[,j] = curr.sample.imputed;
        
      }
      
      j=j+1;
    }
    
    k=k+nb_rep[i];
  }
  
  return(dataSet.imputed)
}






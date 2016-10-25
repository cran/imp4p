
#impute missing values of each column of tab 
#with a small value.
#The small value is randomly draws from a uniform distribution
#between the min of observed values in the column 
#and the q% quantile of these observed values. 

impute.pa=function(tab, conditions, q.min=0, q.norm=3, eps=2){
  
  tab_imp=tab;
  
  qu=apply(tab_imp,2,quantile,na.rm=T,q.min);
  
  nb_cond=length(levels(conditions));
  nb_rep=rep(0,nb_cond);
  k=1;
  j=1;
  for (i in 1:nb_cond){
    
    nb_rep[i]=sum((conditions==levels(conditions)[i]));
    sde=apply(tab_imp[,(k:(k+nb_rep[i]-1))],1,sd,na.rm=T);
    
    while (j<(k+nb_rep[i])){
      tab_imp[which(is.na(tab_imp[,j])),j]=runif(n=sum(is.na(tab_imp[,j])),min=qu[j]-eps-q.norm*median(sde,na.rm=T),max=qu[j]-eps);
      j=j+1;
    }
    
    k=k+nb_rep[i];
  }
  
  return(tab_imp);
}


 
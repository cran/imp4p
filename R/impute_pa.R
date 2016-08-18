
#impute missing values of each column of tab 
#with a small value.
#The small value is randomly draws from a uniform distribution
#between the min of observed values in the column 
#and the q% quantile of these observed values. 

impute.pa=function(tab, q=0.01){
  
  tab_imp=tab;
  
  qu=apply(tab_imp,2,quantile,na.rm=T,q);
  mi=apply(tab_imp,2,min,na.rm=T);
  
  for (j in 1:ncol(tab)){
      tab_imp[which(is.na(tab_imp[,j])),j]=runif(n=sum(is.na(tab_imp[,j])),min=mi[j],max=qu[j]);
  }
  
  return(tab_imp);
}
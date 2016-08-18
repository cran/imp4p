
#return a tab with at least one observed value in each condition.
#liste= liste of deleted rows in the original tab

delete.na.rows=function(tab,tab.c,conditions,list.MCAR){
  nb_cond=length(levels(conditions));
  nb_rep=rep(0,nb_cond);
  k=1;
  liste=NULL;
  for (i in 1:nb_cond){
      nb_rep[i]=sum((conditions==levels(conditions)[i]));
      nb.NA<-apply(tab[,(k:(k+nb_rep[i]-1))],1,function(x) sum(is.na(x)));
      liste=c(liste,which(nb.NA==nb_rep[i]));
      k=k+nb_rep[i];
  }
  
  tab.r=tab[-liste,];
  tab.comp.r=tab.c[-liste,];
  pi.na=apply(tab.r,2,function(x) sum(is.na(x))/length(x));
  nb.mcar=apply(list.MCAR,2,function(x) sum(!x%in%liste));
  pi.mcar=nb.mcar/(pi.na*nrow(tab.r));
  
  return(list(tab.mod=tab.r,tab.comp=tab.comp.r,list.delete=liste,pi.na=pi.na,pi.mcar=pi.mcar))
}



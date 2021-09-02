
impute.rand=function(tab,conditions){

  tab_imp=as.matrix(tab);

  new_tab=tab
  new_tab.imp=tab_imp
  new_conditions=NULL
  index=NULL
  k=1
  for (j in 1:length(levels(conditions))){
    index=c(index,which(conditions==levels(conditions)[j]))
    nb_rep=sum((conditions==levels(conditions)[j]));
    new_tab[,(k:(k+nb_rep-1))]=tab[,which(conditions==levels(conditions)[j])]
    new_tab.imp[,(k:(k+nb_rep-1))]=tab_imp[,which(conditions==levels(conditions)[j])]
    new_conditions=c(new_conditions,conditions[which(conditions==levels(conditions)[j])])
    k=k+nb_rep
  }

  tab=new_tab
  tab_imp=new_tab.imp
  conditions=new_conditions
  conditions=factor(as.character(conditions),levels=as.character(unique(conditions)));

  nb_cond=length(levels(conditions));
  nb_rep=rep(0,nb_cond);
  k=1;

  for (n in 1:nb_cond){

    nb_rep[n]=sum((conditions==levels(conditions)[n]));

    md=suppressWarnings(apply(new_tab[,(k:(k+nb_rep[n]-1))],1,mean,na.rm=T));
    asd=suppressWarnings(apply(new_tab[,(k:(k+nb_rep[n]-1))],1,sd,na.rm=T));
    masd=quantile(asd,0.25,na.rm=T);

    for (i in 1:nrow(tab_imp)){
      tab_imp[i,(k:(k+nb_rep[n]-1))][which(is.na(tab_imp[i,(k:(k+nb_rep[n]-1))]))]=md[i]+rnorm(sum(is.na(tab_imp[i,(k:(k+nb_rep[n]-1))])),0,masd);
    }

    k=k+nb_rep[n];
  }

  tab_imp[is.nan(tab_imp)]=NA;
  tab_imp[,index]=tab_imp
  colnames(tab_imp)=colnames(tab)

  return(tab_imp);
}






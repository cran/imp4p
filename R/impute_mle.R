
impute.mle=function (tab, conditions) {

  tab_imp=as.matrix(tab);

  nb_cond=length(levels(conditions));
  nb_rep=rep(0,nb_cond);
  k=1;

  for (n in 1:nb_cond){
    nb_rep[n]=sum((conditions==levels(conditions)[n]));
    xincomplete=as.matrix(tab[,(k:(k+nb_rep[n]-1))]);
    nbna=apply(xincomplete,1,function(x){sum(is.na(x));})
    xincomplete1=xincomplete[which(nbna!=nb_rep[n]),];
    s <- prelim.norm(xincomplete1);
    thetahat <- em.norm(s, showits = FALSE);
    rngseed(1234567);
    xcomplete1 <- imp.norm(s, thetahat, xincomplete1);
    tab_imp[which(nbna!=nb_rep[n]),(k:(k+nb_rep[n]-1))]=xcomplete1;
    k=k+nb_rep[n];
  }

  return(tab_imp)
}

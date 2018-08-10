
impute.mle=function (tab, conditions) {
  
  tab_imp=as.matrix(tab);
  
  nb_cond=length(levels(conditions));
  nb_rep=rep(0,nb_cond);
  k=1;
  
  for (n in 1:nb_cond){
    nb_rep[n]=sum((conditions==levels(conditions)[n]));
    xincomplete=as.matrix(tab[,(k:(k+nb_rep[n]-1))]);
    s <- prelim.norm(xincomplete);
    thetahat <- em.norm(s, showits = FALSE);
    rngseed(1234567);
    xcomplete <- imp.norm(s, thetahat, xincomplete);
    tab_imp[,(k:(k+nb_rep[n]-1))]=xcomplete;
    k=k+nb_rep[n];
  }
  
  return(tab_imp)
}
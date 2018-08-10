
impute.mi=function(tab, conditions, repbio=NULL, reptech=NULL, nb.iter=3, nknn=15, selec=1000, siz=900, weight=1, ind.comp=1, progress.bar=TRUE, x.step.mod=300, x.step.pi=300, nb.rei=100, method=1, gridsize=300, q=0.95, q.min=0, q.norm=3, eps=2, methodi="slsa"){

  if (is.null(repbio)){repbio=as.factor(1:length(conditions));}
  if (is.null(reptech)){reptech=as.factor(1:length(conditions));}

  if (progress.bar==TRUE){cat(paste("\n 1/ Initial imputation under the MCAR assumption... \n  "));}

  #Imputation of missing values with the slsa algorithm
  dat.slsa=impute.slsa(tab=tab, conditions=conditions, repbio=repbio, reptech=reptech, nknn=nknn, selec=selec, weight=weight, ind.comp=ind.comp, progress.bar=progress.bar);
  if (methodi == "mle"){
     dat.slsa=impute.mle(tab=tab, conditions=conditions);
  }

  if (progress.bar==TRUE){cat(paste("\n\n 2/ Estimation of the mixture model in each sample... \n  "));}

  #Estimation of the mixture model
  res=estim.mix(tab=tab, tab.imp=dat.slsa, conditions=conditions, x.step.mod=x.step.mod, x.step.pi=x.step.pi, nb.rei=nb.rei, method=method, gridsize=gridsize);

  if (progress.bar==TRUE){cat(paste("\n 3/ Estimation of the probabilities each missing value is MCAR... \n  "));}

  #Computing probabilities to be MCAR
  born=estim.bound(tab=tab,conditions=conditions,q=q);
  proba=prob.mcar.tab(born$tab.lower,born$tab.upper,res);

  if (progress.bar==TRUE){cat(paste("\n 4/ Multiple imputation strategy with mi.mix... \n  "));}

  #Multiple imputation strategy
  data.mi=mi.mix(tab=tab, tab.imp=dat.slsa, prob.MCAR=proba, conditions=conditions, repbio=repbio, reptech=reptech, nb.iter=nb.iter, nknn=nknn, weight=weight, selec=selec, siz=siz, ind.comp=ind.comp, methodi=methodi, q=q, progress.bar=progress.bar);

  if (progress.bar==TRUE){cat(paste("\n\n 5/ Imputation of rows with only missing values in a condition with impute.pa... \n  "));}

  data.final=impute.pa(tab=data.mi, conditions=conditions, q.min=q.min, q.norm=q.norm, eps=eps);

  return(data.final)
}








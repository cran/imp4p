
#Function to perform igcda algorithm for selected MNAR values and slsa algorithm or mle algorithm for MCAR values

impute.mix <-function (tab, prob.MCAR, conditions, repbio=NULL, reptech=NULL, method="slsa", nknn=15, weight=1, selec="all", ind.comp=1, progress.bar=TRUE, q=0.95){
  
  if (is.null(repbio)){repbio=as.factor(1:length(conditions));}
  if (is.null(reptech)){reptech=as.factor(1:length(conditions));}
  
  tab.mvs=tab;
    
  #Random draw of MCAR and MNAR values
  l.MCAR=matrix(0,nrow(tab),ncol(tab));
  l.MCAR[which(prob.MCAR>0.5)]=1;
  l.MCAR[which(!is.na(tab))]=0;

  #Impute MCAR values
  if (method=="slsa"){
      tab.imp=impute.slsa(tab=tab.mvs, conditions=conditions, repbio=repbio, reptech=reptech, nknn=nknn, selec=min(selec,nrow(tab.mvs)), weight=weight, ind.comp=ind.comp, progress.bar=progress.bar);
  }else{
      tab.imp=impute.mle(tab=tab.mvs, conditions=conditions);
  }
  tab.mvs[which(l.MCAR==1)]=tab.imp[which(l.MCAR==1)];
    
  #Impute MNAR values 
  tab.mvs.imp=impute.igcda(tab=tab.mvs, tab.imp=tab.imp, conditions=conditions, q=q);

  return(tab.mvs.imp);
}




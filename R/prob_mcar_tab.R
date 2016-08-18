
#Compute a matrix of probabilities to be MCAR

prob.mcar.tab=function(tab.l,tab.u,res){
  nri=nrow(tab.l);
  nci=ncol(tab.l);
  prob=matrix(0,nri,nci);
  for (j in 1:nci){
    prob[,j]=prob.mcar(b.l=tab.l[,j], b.u=tab.u[,j], absc=res$abs.mod, pi.mcar=res$pi.mcar[j], F.tot=res$F.tot[,j], F.na=res$F.na[,j]);
  }
  return(prob)
}
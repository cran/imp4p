miss.mnar.process=function(abs,pi_mcar,F_mnar,F_na){
  prob.miss.mnar=NULL
  for (x in 1:(length(abs)-1)){
       prob.miss.mnar=c(prob.miss.mnar,(1-pi_mcar)*(F_mnar[x+1]-F_mnar[x])/(F_na[x+1]-F_na[x]));
  }
  return(list(abs=abs[1:(length(abs)-1)],p=prob.miss.mnar))
}

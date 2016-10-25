
#Compute a vector of probabilities to be MCAR

prob.mcar=function(b.l,b.u,absc,pi.mcar,F.tot,F.na){
  prob=rep(0,length(b.l));
  for (i in 1:length(b.l)){
    
      if (!is.na(b.l[i])){  
      if (!is.na(b.u[i])){
      if (!is.nan(b.l[i])){
      if (!is.nan(b.u[i])){  
      if (b.l[i]!=Inf){  
      if (b.u[i]!=Inf){

        
        i.binf=which.min(abs(absc-b.l[i]));
        i.bsup=which.min(abs(absc-b.u[i]));
        if ( ((F.na[i.bsup]-F.na[i.binf])!=0) ){
          prob[i]=pi.mcar*(F.tot[i.bsup]-F.tot[i.binf])/(F.na[i.bsup]-F.na[i.binf]);
        }else{
          while((F.na[i.bsup]==F.na[i.binf])&(i.binf>1)){
              i.binf=i.binf-1;
          }
          if (i.binf==1){
            while((F.na[i.bsup]==F.na[i.binf])){i.bsup=i.bsup+1;}
          }
          prob[i]=pi.mcar*(F.tot[i.bsup]-F.tot[i.binf])/(F.na[i.bsup]-F.na[i.binf]);
        }
        
      }
      }
      }
      }
      }
      }

  }
  
  prob[prob>1]=1;
  prob[prob<0]=0;
  
  return(prob)
}
 
 
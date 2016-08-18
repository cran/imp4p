#Simulation of data sets
#with a hierarchical design

sim.data=function(nb.pept=2000,nb.miss=600,pi.mcar=0.2,para=10,nb.cond=2,nb.repbio=3,nb.sample=5,m.c=25,sd.c=2,sd.rb=0.5,sd.r=0.2){
  
  tab=matrix(0,nb.pept,nb.sample*nb.cond*nb.repbio);
  pmnar=tab;
  moyenne.cond=matrix(0,nb.pept,nb.cond);
  moyenne.repbio=matrix(0,nb.pept,nb.repbio);
  
  #Generating the averages and the probabilities to be MNAR
  for (i in 0:(nb.cond-1)){
    moyenne.cond[,i+1]=rnorm(nb.pept,mean=m.c,sd=sd.c);
    for (rb in 0:(nb.repbio-1)){
      moyenne.repbio[,rb+1]=rnorm(nb.pept,mean=0,sd=sd.rb);
      mimo=min(moyenne.cond[,i+1]+moyenne.repbio[,rb+1]);
      mamo=max(moyenne.cond[,i+1]+moyenne.repbio[,rb+1]);
      for (j in (nb.sample*(rb+i*nb.repbio)):(nb.sample*(rb+i*nb.repbio)+nb.sample-1)){
        for (k in 1:nb.pept){
           pmnar[k,j+1]=dbeta((moyenne.cond[k,i+1]+moyenne.repbio[k,rb+1]-mimo)/(mamo-mimo),1,para)/para;
        }
      }
    }
  }
  
  
  tab.comp=tab;

  #Drawing to generate MCAR values
  nbMCAR=nb.miss*pi.mcar
  nbMNAR=nb.miss*(1-pi.mcar)
  listeMCAR=matrix(0,nbMCAR,nb.sample*nb.cond*nb.repbio)
  
  if (pi.mcar!=0){
    for (i in 0:(nb.cond-1)){
      
      for (rb in 0:(nb.repbio-1)){
        
        for (j in (nb.sample*(rb+i*nb.repbio)):(nb.sample*(rb+i*nb.repbio)+nb.sample-1)){
          liste=sample(1:nb.pept,size=nbMCAR);
          tab[liste,j+1]=NA;
          listeMCAR[,j+1]=liste;
        }
        
      }
    }
  }

  #Drawing to find MNAR values
  if (pi.mcar!=1){
    for (i in 0:(nb.cond-1)){
      
      for (rb in 0:(nb.repbio-1)){
        
        for (j in (nb.sample*(rb+i*nb.repbio)):(nb.sample*(rb+i*nb.repbio)+nb.sample-1)){
          
          d=1:nb.pept
          sa=d[-c(which(is.na(tab[,j+1])))];
          if (length(sa)>=nbMNAR){
              liste=sample(sa,size=nbMNAR,prob=pmnar[-c(which(is.na(tab[,j+1]))),j+1]);
          }else{warning(paste("\n The proportion of generated MCAR values is superior to the chosen proportion in sample ",j+1,".\n"));
              liste=sa;
          }

          tab[liste,j+1]=NA;  

        }
        
      }
      
    }
    
    #Generating the ground truth of the data
    for (i in 0:(nb.cond-1)){
        
        for (rb in 0:(nb.repbio-1)){
          
          for (j in (nb.sample*(rb+i*nb.repbio)+1):(nb.sample*(rb+i*nb.repbio)+nb.sample)){            
              
              #MNAR values
              listeMNAR=which(is.na(tab[-listeMCAR[,j],j]));
              qq=qnorm(pmnar[-listeMCAR[,j],j][listeMNAR],mean=moyenne.cond[-listeMCAR[,j],i+1][listeMNAR]+moyenne.repbio[-listeMCAR[,j],rb+1][listeMNAR],sd=sd.r);
              tab.comp[-listeMCAR[,j],j][listeMNAR]=rtruncnorm(length(listeMNAR),mean=moyenne.cond[-listeMCAR[,j],i+1][listeMNAR]+moyenne.repbio[-listeMCAR[,j],rb+1][listeMNAR], sd=sd.r, b=qq);   
              
              #Observed values
              qq=qnorm(pmnar[-listeMCAR[,j],j][-listeMNAR],mean=moyenne.cond[-listeMCAR[,j],i+1][-listeMNAR]+moyenne.repbio[-listeMCAR[,j],rb+1][-listeMNAR],sd=sd.r);
              tab.comp[-listeMCAR[,j],j][-listeMNAR]=rtruncnorm(length(tab.comp[-listeMCAR[,1],1][-listeMNAR]),mean=moyenne.cond[-listeMCAR[,j],i+1][-listeMNAR]+moyenne.repbio[-listeMCAR[,j],rb+1][-listeMNAR], sd=sd.r, a=qq);   
              
              #MCAR values
              tab.comp[listeMCAR[,j],j]=rnorm(length(tab.comp[listeMCAR[,1],1]),mean=moyenne.cond[listeMCAR[,j],i+1]+moyenne.repbio[listeMCAR[,j],rb+1], sd=sd.r);   
          }
          
        }
    
    }
    tab.r=tab.comp;
    tab.comp[which(is.na(tab))]=NA;
  }
  
  return(list(dat.obs=tab.comp, dat.comp=tab.r, list.MCAR=listeMCAR, conditions=gen.cond(nb.cond,nb.repbio*nb.sample), repbio=gen.cond(nb.cond*nb.repbio,nb.sample)))
}





#library(Rcpp)
#sourceCpp("C:/Users/Quentin/Documents/article_imp4p/article/correlation_mod2.cpp")
#sourceCpp("C:/Users/Quentin/Documents/article_imp4p/article/fast_apply_nb_na.cpp")
#sourceCpp("C:/Users/Quentin/Documents/article_imp4p/article/fast_apply_nb_not_na.cpp")
#sourceCpp("C:/Users/Quentin/Documents/article_imp4p/article/fast_apply_sd_na_rm_T.cpp")
#sourceCpp("C:/Users/Quentin/Documents/article_imp4p/article/fast_apply_sum_na_rm_T.cpp")
#sourceCpp("C:/Users/Quentin/Documents/article_imp4p/article/fast_sim.cpp")
#Function to impute data sets with the SLSA algorithm

impute.slsa=function(tab, conditions, repbio=NULL, reptech=NULL, nknn=15, selec="all", weight=1, ind.comp=1, progress.bar=TRUE){

  if (is.null(repbio)){repbio=as.factor(1:length(conditions));}
  if (is.null(reptech)){reptech=as.factor(1:length(conditions));}

  tab_imp=as.matrix(tab);

  nb_cond=length(levels(conditions));
  nb_rep=rep(0,nb_cond);
  k=1;

  for (n in 1:nb_cond){

    if (progress.bar==TRUE){
      cat(paste("\n Imputation in condition ",n,"...\n In progress: "));
    }

    nb_rep[n]=sum((conditions==levels(conditions)[n]));
    repb=repbio[(k:(k+nb_rep[n]-1))];
    rept=reptech[(k:(k+nb_rep[n]-1))];

    xincomplete1=tab_imp[,(k:(k+nb_rep[n]-1))];

    asd=fast_apply_sd_na_rm_T(xincomplete1,1);#apply(xincomplete1,1,sd,na.rm=T);

    #Imputation function for each row
    imputeHLS=function(roww){
      #Progress bar
      if (progress.bar==TRUE){
        if (roww[1]%/%(floor(0.1*nrow(xincomplete1)))==roww[1]/(floor(0.1*nrow(xincomplete1)))){cat(c(0.1*(roww[1]%/%(floor(0.1*nrow(xincomplete1))))*100,"% - "));}
      }
      #print(roww[1])
      roww=roww[2:length(roww)];
      naroww=is.na(roww);
      # if at least one missing value else no imputation of the row
      if (sum(naroww)!=0){
        #if at least one observed value else no imputation of the row
        if (sum(naroww)!=length(roww)){

          row.exp=which(naroww);
          prot=roww[-row.exp];

          if (is.numeric(selec)){
            xincomplete=xincomplete1[base::sample.int(n=nrow(xincomplete1),size=selec), ];
          }else{
            xincomplete=xincomplete1;
          }
          cand_x=xincomplete[, -row.exp, drop=F];

          #similarity measures:
          #pairwise correlation if at least three side-by-side observed values
          #euclidean distance between side-by-side observed values if length(prot)<3
          sim=fast_sim(prot,cand_x);
          sim=abs(sim);
          osim=order(sim, decreasing=T);
          #row.cand is the set of values from which linear models are fitted
          #row.impcand is the set of values allowing to predict responses from fitted linear models
          #At least nknn observed values are needed in each column of row.cand and row.impcand
          Kmin=0;
          kl=1;
          ind.comp.temp=ind.comp
          while (Kmin<(nknn)){
            kl=kl+1;
            row.idx=osim[2:(kl+1)];
            row.r=sim[row.idx];

            row.cand=cand_x[row.idx, , drop=F];
            indcand=fast_apply_nb_na(row.cand,1);
            indic=which(indcand!=length(row.cand[1,]));
            if (ind.comp.temp==1){indic=which(indcand==0);}
            row.idx=row.idx[indic];
            row.r=row.r[indic];

            row.cand=row.cand[indic,];

            row.impcand=xincomplete[row.idx, row.exp, drop=F];

            if (is.matrix(row.cand)){
              nborc=fast_apply_nb_not_na(row.cand,2);
            }else{nborc=sum(!is.na(row.cand));}
            nboric=fast_apply_nb_not_na(row.impcand,2);
            Kmin=min(c(nborc,nboric));
            if (ind.comp.temp==1){
              if (kl>=700){
                    kl=1;ind.comp.temp=0;
              }
              #if (ind.comp==0){Kmin=nknn;}
            }
          }

          #Let's go to fit linear models and responses
          #1) creating design matrices
          rrepb=as.factor(as.numeric(repb[-row.exp]));
          rrept=as.factor(as.numeric(rept[-row.exp]));
          rrepb1=as.factor(as.numeric(repb[row.exp]));
          rrept1=as.factor(as.numeric(rept[row.exp]));
          lrb=levels(as.factor(as.numeric(repb)));
          lrt=levels(as.factor(as.numeric(rept)));
          llrb=length(lrb);
          llrt=length(lrt);
          mm=rep(1,length(prot));
          mm1=rep(1,length(row.exp));
          if (llrb>1){
            for (ki in 2:llrb){
              mm=cbind(mm,1*(rrepb==lrb[ki]));
              mm1=cbind(mm1,1*(rrepb1==lrb[ki]));
            }
          }
          if (llrt>1){
            for (ki in 2:llrt){
              mm=cbind(mm,1*(rrept==lrt[ki]));
              mm1=cbind(mm1,1*(rrept1==lrt[ki]));
            }
          }

          if (nrow(mm)>1){
            mm[,c(FALSE,(apply(mm[,2:length(mm[1,])],2,sum)<=1))]=rep(0,length(mm[,1]));
          }else{
            mm[2:length(mm)][mm[2:length(mm)]<=1]=0;}

          #2) compute matrix of responses (y)
          nri=nrow(row.impcand);
          nci=ncol(row.impcand);
          y=matrix(0, nrow=nri, ncol=nci);

          for(i in 1:nri) {
            if (is.matrix(row.cand)){
              mx=cbind(mm,row.cand[i,]);
              ll=!is.na(row.cand[i,]);
              mx=mx[ll,];
              if (is.matrix(mx)){
                lg=lm.fit(x=mx,y=prot[ll])$coefficients;
                lg[is.na(lg)]=0;
              }else{lg=rep(0,length(mx));
              lg[length(mx)]=prot[ll]/mx[length(mx)];
              }
            }else{lg=rep(0,length(mm)+1);
            lg[length(lg)]=prot/row.cand[i];
            }

            mx1=cbind(mm1,row.impcand[i,]);
            y[i, ]=mx1%*%lg;
          }

          #delete outliers
          for (j in 1:ncol(y)){
            qyj=quantile(y[,j],na.rm=T);
            bi=qyj[3]-0.5*(qyj[4]-qyj[2]);
            bs=qyj[3]+0.5*(qyj[4]-qyj[2]);
            y[y[,j]<bi,j]=NA;
            y[y[,j]>bs,j]=NA;
          }

          #compute the weight function
          if (is.numeric(weight)==TRUE){w=row.r^weight;}
          if (weight=="o"){w=(row.r**2/(1-row.r**2+0.000001))**2;}

          #final imputation by weighting the observed responses
          roww[row.exp]<-apply(y, 2, function(x){xx=x[!is.na(x)];ww=w[!is.na(x)];ww[1:min(c(nknn,length(ww)))]=ww[1:min(c(nknn,length(ww)))]/sum(ww[1:min(c(nknn,length(ww)))]);sum(ww[1:min(c(nknn,length(ww)))]*xx[1:min(c(nknn,length(ww)))]);})

          #to limit imputations in a reasonable neighborhood
          l.ol=sum(roww[row.exp]<(mean(roww[-row.exp])-2*quantile(asd,probs=0.9,na.rm=T)));
          l.ou=sum(roww[row.exp]>(mean(roww[-row.exp])+2*quantile(asd,probs=0.9,na.rm=T)));
          if (l.ol>0){roww[row.exp][(roww[row.exp]<(mean(roww[-row.exp])-2*quantile(asd,probs=0.9,na.rm=T)))]=(mean(roww[-row.exp])-rnorm(length(roww[row.exp][(roww[row.exp]<(mean(roww[-row.exp])-2*quantile(asd,probs=0.9,na.rm=T)))]),mean=2,sd=0.005)*quantile(asd,probs=0.9,na.rm=T));};
          if (l.ou>0){roww[row.exp][(roww[row.exp]>(mean(roww[-row.exp])+2*quantile(asd,probs=0.9,na.rm=T)))]=(mean(roww[-row.exp])+rnorm(length(roww[row.exp][(roww[row.exp]>(mean(roww[-row.exp])+2*quantile(asd,probs=0.9,na.rm=T)))]),mean=2,sd=0.005)*quantile(asd,probs=0.9,na.rm=T));};
        }
      }

      return(roww);
    }

    xmiss = t(apply(cbind(seq_len(nrow(xincomplete1)),xincomplete1), 1, imputeHLS));

    tab_imp[,(k:(k+nb_rep[n]-1))] = xmiss;

    k=k+nb_rep[n];
  }

  return (tab_imp)
}

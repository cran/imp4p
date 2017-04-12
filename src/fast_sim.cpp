
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

inline NumericVector fast_si( NumericVector& prot, NumericMatrix& mat) {
  
  int nr = mat.nrow();
  
  int n = prot.size();
  LogicalVector out_prot(n);
  LogicalVector out_mat(n);
  
  double m_p=0;
  double m_m=0;
  double k=0;
  double v_p=0;
  double v_m=0;
  double cov_mp=0;
  double d_m=0;
  NumericVector prot2(n);
  NumericVector corr(nr);
  for (int j = 0; j < nr; ++j) {
    prot2=mat(j,_);
    m_p=0;
    m_m=0;
    k=0;
    d_m=0;
    //selection des observations pairwise
    for (int i = 0; i < n; ++i) {
      //prot2[i]=mat(j,i);
      out_prot[i] = NumericVector::is_na(prot[i]);
      out_mat[i] = NumericVector::is_na(prot2[i]);
      if ( out_prot[i] == FALSE ){
        if ( out_mat[i] == FALSE ){
          m_p+=prot[i];
          m_m+=prot2[i];
          d_m+=(prot[i]-prot2[i])*(prot[i]-prot2[i]);
          ++k;
        }
      }
    }
    //Plus de 2 valeurs observées de chaque côtés -> corrélation simple
    if ( k>2 ){
        m_p=m_p/k;
        m_m=m_m/k;
        v_p=0;
        v_m=0;
        cov_mp=0;
        for (int i = 0; i < n; ++i) {
            out_prot[i] = NumericVector::is_na(prot[i]);
            out_mat[i] = NumericVector::is_na(prot2[i]);
            if ( out_prot[i] == FALSE ){
              if ( out_mat[i] == FALSE ){
                v_p += (prot[i]-m_p)*(prot[i]-m_p);
                v_m += (prot2[i]-m_m)*(prot2[i]-m_m);
                cov_mp += (prot[i]-m_p)*(prot2[i]-m_m);
              }
            }
        }
        if ((v_m != 0)){
          if ((v_p != 0)){
              corr[j]=(cov_mp)/( sqrt(v_m)*sqrt(v_p) );
          }
          else{corr[j]=NA_REAL;}
        }else{corr[j]=NA_REAL;}
    }
    //Moins de 3 valeurs observées de chaque côté
    if ( k<3 ){
        if ( k>0 ){
            //Si prot a moins de trois valeurs : distance euclidienne entre prot et prot2
            if ( n<3 ){
                corr[j]=1/(1+sqrt(d_m));
            }else{
              //sinon si prot a plus de deux valeurs et moins de 3 observations "pairwise" avec prot2: 
              //on met NA
              corr[j]=NA_REAL;                
            }
        }else{corr[j]=NA_REAL;}
    }
  }
  
  return corr;
}

// [[Rcpp::export]]
NumericVector fast_sim(NumericVector prot, NumericMatrix mat) {
  return fast_si(prot, mat);
}

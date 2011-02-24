/*
(c)2011 Jyotirmoy Bhattacharya
Licence: GPL
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

typedef struct {
  int nvars;

  double *c_y;
  double *c_x;
  int c_nobs;

  double *uc_y;
  double *uc_x;
  int  uc_nobs;
} Data;

void *check_alloc(size_t size)
{
  void *p=malloc(size);
  if (p==NULL){
    fprintf(stderr,"Memory allocation failed\n");
    exit(1);
  }
  return p;
}

double *vecdup(double *src,int len)
{
  double *dst=(double *)check_alloc(sizeof(double)*len);
  memcpy(dst,src,sizeof(double)*len);
  return dst;
}

void data_finalizer(SEXP handle)
{
  Data *data=R_ExternalPtrAddr(handle);
  if (data==NULL)
    return;

  free(data->c_y);
  free(data->c_x);
  free(data->uc_y);
  free(data->uc_x);
  free(data);

  R_ClearExternalPtr(handle);
}

SEXP jmoy_loaddata(SEXP nvars,
		   SEXP c_y,SEXP c_x,SEXP c_nobs,
		   SEXP uc_y,SEXP uc_x,SEXP uc_nobs)
{
  Data *data;
  SEXP handle;
  
  data=(Data *)check_alloc(sizeof(Data));
  data->nvars=INTEGER(nvars)[0];

  data->c_nobs=INTEGER(c_nobs)[0];
  data->c_y=vecdup(REAL(c_y),data->c_nobs);
  data->c_x=vecdup(REAL(c_x),data->c_nobs*data->nvars);

  data->uc_nobs=INTEGER(uc_nobs)[0];
  data->uc_y=vecdup(REAL(uc_y),data->uc_nobs);
  data->uc_x=vecdup(REAL(uc_x),data->uc_nobs*data->nvars);

  handle=R_MakeExternalPtr(data,R_NilValue,R_NilValue);
  PROTECT(handle);
  R_RegisterCFinalizer(handle,data_finalizer);
  UNPROTECT(1); /*handle*/
  return handle;
}

SEXP jmoy_negloglik(SEXP handle,SEXP sbeta,SEXP salpha){
  Data *data;
  SEXP ans;

  double *beta;
  double alpha;
  double loglik=0.0;

  int i,j;
  
  data=R_ExternalPtrAddr(handle);
  beta=REAL(sbeta);
  alpha=REAL(salpha)[0];

  /*Uncensored observations*/
  
  for(i=0;i<data->uc_nobs;i++){
    double xbeta=0.0;
    double t=data->uc_y[i];

    for(j=0;j<data->nvars;j++)
      xbeta+=data->uc_x[i+j*data->uc_nobs]*beta[j];

    loglik+=xbeta+log(alpha*pow(t,alpha-1))-exp(xbeta)*pow(t,alpha);
  }

  /*Censored observations*/
  
  for(i=0;i<data->c_nobs;i++){
    double xbeta=0.0;
    double t=data->c_y[i];

    for(j=0;j<data->nvars;j++)
      xbeta+=data->c_x[i+j*data->c_nobs]*beta[j];

    loglik+=-exp(xbeta)*pow(t,alpha);
  }

  ans=allocVector(REALSXP,1);
  REAL(ans)[0]=-loglik;

  return ans;
}

SEXP jmoy_gradient(SEXP handle,SEXP sbeta,SEXP salpha){
  Data *data;
  SEXP sgrad;

  double *beta;
  double alpha;
  double *grad;

  int i,j;
  
  data=R_ExternalPtrAddr(handle);
  beta=REAL(sbeta);
  alpha=REAL(salpha)[0];
  PROTECT(sgrad=allocVector(REALSXP,data->nvars+1));
  grad=REAL(sgrad);
  for(i=0;i<data->nvars+1;i++)
    grad[i]=0;

  /*Uncensored observations*/
  
  for(i=0;i<data->uc_nobs;i++){
    double xbeta=0.0;
    double t=data->uc_y[i];
    double t_alpha=pow(t,alpha);
    double exp_xbeta;

    double grad_xbeta;

    for(j=0;j<data->nvars;j++)
      xbeta+=data->uc_x[i+j*data->uc_nobs]*beta[j];
    exp_xbeta=exp(xbeta);

    grad_xbeta=1.0-exp_xbeta*t_alpha;
    for (j=0;j<data->nvars;j++)
      grad[j]-=data->uc_x[i+j*data->uc_nobs]*grad_xbeta;
    
    grad[data->nvars]-=
      1.0/alpha
      +(1.0-exp_xbeta*t_alpha)*log(t);
  }

  /*Censored observations*/
  
  for(i=0;i<data->c_nobs;i++){
    double xbeta=0.0;
    double t=data->c_y[i];
    double t_alpha=pow(t,alpha);
    double exp_xbeta;

    double grad_xbeta;

    for(j=0;j<data->nvars;j++)
      xbeta+=data->c_x[i+j*data->c_nobs]*beta[j];
    exp_xbeta=exp(xbeta);

    grad_xbeta=-exp_xbeta*t_alpha;
    for (j=0;j<data->nvars;j++)
      grad[j]-=data->c_x[i+j*data->c_nobs]*grad_xbeta;
    grad[data->nvars]-=-exp_xbeta*t_alpha*log(t);
  }

  UNPROTECT(1); /*sgrad*/

  return sgrad;
}
  
R_CallMethodDef callMethods[]={
  {"jmoy_loaddata",(DL_FUNC) &jmoy_loaddata,7},
  {"jmoy_negloglik",(DL_FUNC) &jmoy_negloglik,3},
  {"jmoy_gradient",(DL_FUNC) &jmoy_gradient,3},
  {NULL,NULL,0}
};

void R_init_recid_c(DllInfo *info)
{
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
}

void R_unload_mylib(DllInfo *info)
{
  ;/* Release resources. */
}

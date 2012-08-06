#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

/* Funcao que calcula a densidade acumulada para L, variavel latente considerando uma distribuicao normal*/

void acumN(int *y, double *L, int *esc, int *max, int *min, int *tame, int *n, double *vero, double *m, double *tau, double *ve, double *verossim)
 {
  int register i,j;
  double ident;
  GetRNGstate();
  *verossim = 0.0;
  for(i = 0; i < *n; i++) {
    ident = runif(0.0,1.0);
    if(y[i] == *max) {
      L[i]  = qnorm((pnorm(tau[*tame-2],m[i],*ve,1,0) + ((1-pnorm(tau[*tame-2],m[i],*ve,1,0))*ident)), m[i], sqrt(*ve),1,0);
      vero[i] = 1 - pnorm(-m[i]/sqrt(*ve),0,1,1,0);
    }
    else if (y[i] == *min) {
      L[i]    = qnorm((pnorm(tau[0],m[i],*ve,1,0))*ident, m[i], sqrt(*ve),1,0);
      vero[i] = pnorm(((tau[0] - m[i])/sqrt(*ve)),0,1,1,0);
    }
    else {
      for(j = 1; j < (*tame-1); j++) {
        if(y[i] == esc[j]) {
          L[i]    = qnorm((pnorm(tau[(j-1)], m[i], *ve,1,0) + (pnorm(tau[j], m[i], *ve,1,0) - pnorm(tau[(j-1)],m[i],*ve,1,0))*ident),m[i],sqrt(*ve),1,0);
          vero[i] = pnorm(((tau[j]-m[i])/sqrt(*ve)),0,1,1,0)-pnorm(((tau[(j-1)]-m[i])/sqrt(*ve)),0,1,1,0);
      }
     }
    }
    *verossim += log(vero[i]);
   }
   PutRNGstate();
  }




/* Funcao que calcula a densidade acumulada para L, variavel latente considerando uma distribuicao t-Student */


void acumt(int *y, double *L, int *esc, int *max, int *min, int *tame, int *n, double *vero, double *m, double *tau, double *ve,
          double *verossim)
 {
  int register i,j;
  double ident;
  *verossim = 0.0;
  GetRNGstate();
  for(i = 0; i < *n; i++) {
    ident = runif(0.0,1.0);
    if(y[i] == *max) {
      L[i]  = qnorm((pnorm(tau[*tame-2],m[i],ve[i],1,0) + ((1-pnorm(tau[*tame-2],m[i],ve[i],1,0))*ident)), m[i], ve[i],1,0);
      vero[i] = 1 - pnorm(-m[i]/ve[i],0,1,1,0);
    }
    else if (y[i] == *min) {
      L[i]    = qnorm((pnorm(tau[0],m[i],ve[i],1,0))*ident, m[i], ve[i],1,0);
      vero[i] = pnorm(((tau[0] - m[i])/ve[i]),0,1,1,0);
    }
    else {
      for(j = 1; j < (*tame-1); j++) {
        if(y[i] == esc[j]) {
          L[i]    = qnorm((pnorm(tau[(j-1)], m[i], ve[i],1,0) + (pnorm(tau[j], m[i], ve[i],1,0) - pnorm(tau[(j-1)],m[i],ve[i],1,0))*ident),m[i],ve[i],1,0);
          vero[i] = pnorm(((tau[j]-m[i])/ve[i]),0,1,1,0)-pnorm(((tau[(j-1)]-m[i])/ve[i]),0,1,1,0);
      }
     }
    }
    *verossim += log(vero[i]);
   }
   PutRNGstate();
  }

void calcpN(double *tau, double *taunovo, double *dpc, double *mp, int *nl)
{
  int register i,j;
  int id;
  id = *nl;
  for(i = 0; i < 4; i++){
    for(j = 0; j < *nl; j++){
      if(i == 0) mp[i*id + j] = pnorm(((tau[j+2] - tau[j+1])/(*dpc)), 0, 1, 1,0);
      else if(i == 1) mp[i*id + j] = pnorm(((taunovo[j] - tau[j+1])/(*dpc)), 0, 1, 1, 0);
      else if(i == 2) mp[i*id + j] = pnorm(((taunovo[j+2] - taunovo[j+1])/(*dpc)), 0, 1, 1, 0);
      else {
        mp[i*id + j] = pnorm(((tau[j] - taunovo[j+1])/(*dpc)), 0, 1, 1, 0);
      }
    }
  }
}

void calcpt(double *tau, double *taunovo, double *dpc, double *mp, int *nl)
{
  int register i,j;
  int id;
  id = *nl;
  for(i = 0; i < 4; i++){
    for(j = 0; j < *nl; j++){
      if(i == 0) mp[i*id + j] = pnorm(((tau[j+2] - tau[j+1])/(*dpc)), 0, 1, 1,0);
      else if(i == 1) mp[i*id + j] = pnorm(((taunovo[j] - tau[j+1])/(*dpc)), 0, 1, 1, 0);
      else if(i == 2) mp[i*id + j] = pnorm(((taunovo[j+2] - taunovo[j+1])/(*dpc)), 0, 1, 1, 0);
      else {
        mp[i*id + j] = pnorm(((tau[j] - taunovo[j+1])/(*dpc)), 0, 1, 1, 0);
      }
    }
  }
}

void taunN(double *taunovo, double *tau, double *dp, int *k)
{
 int register i;
 double area,area1,area2;
 GetRNGstate();
 for(i = 1; i < (*k-1); i++){
   area1 = pnorm(taunovo[i-1], tau[i], *dp, 1, 0);
   area2 = pnorm(tau[i+1], tau[i], *dp, 1, 0) - pnorm(taunovo[i-1], tau[i], *dp, 1, 0);
   area  = area1 + runif(0.0,1.0)*area2;
   taunovo[i] = qnorm(area,tau[i],*dp,1,0);
 }
 PutRNGstate();
}


void taunt(double *taunovo, double *tau, double *dpc, int *tame)
{
 int register i;
 double area,area1,area2;
 GetRNGstate();
 for(i = 1; i < *tame; i++){
   area1 = pnorm(taunovo[i-1], tau[i], *dpc, 1, 0);
   area2 = pnorm(tau[i+1], tau[i], *dpc, 1, 0) - pnorm(taunovo[i-1], tau[i], *dpc, 1, 0);
   area  = area1 + runif(0.0,1.0)*area2;
   taunovo[i] = qnorm(area,tau[i],*dpc,1,0);
 }
 PutRNGstate();
}



void probtauN(int *y, double *ftau, double *ftaunovo, int *esc, int *max, int *min, int *tame, int *n, double *m, double *tau, double *taunovo)
 {
  int register i,j;
  double ident;
  for(i = 0; i < *n; i++) {
    if(y[i] == *max) {
      ftau[i]     = 1 - pnorm(tau[(*tame)-2] - m[i], 0, 1, 1, 0);
      ftaunovo[i] = 1 - pnorm(taunovo[(*tame)-2] - m[i], 0, 1, 1, 0);
    }
    else if (y[i] == *min) {
      ftau[i]     = pnorm(tau[0] - m[i], 0, 1, 1, 0);
      ftaunovo[i] = pnorm(taunovo[0] - m[i], 0, 1, 1, 0);
    }
    else {
      for(j = 1; j < (*tame-1); j++) {
        if(y[i] == esc[j]) {
	  ftau[i]  = pnorm((tau[j] - m[i]), 0, 1, 1, 0) - pnorm((tau[j-1] - m[i]), 0, 1, 1, 0);
	  ftaunovo[i] = pnorm((taunovo[j] - m[i]), 0, 1, 1, 0) - pnorm((taunovo[j-1] - m[i]), 0, 1, 1, 0);
     }
    }
   }
  }
}



void probtaut(int *y, double *ftau, double *ftaunovo, int *esc, int *max, int *min, int *tame, int *n, double *m, double *tau, double *taunovo, double *dp)
 {
  int register i,j;
  double ident;
  for(i = 0; i < *n; i++) {
    if(y[i] == *max) {
      ftau[i]     = 1 - pnorm(((tau[*tame-2] - m[i])/dp[i]), 0, 1, 1, 0);
      ftaunovo[i] = 1 - pnorm(((taunovo[*tame-2] - m[i])/dp[i]), 0, 1, 1, 0);
    }
    else if (y[i] == *min) {
      ftau[i]     = pnorm(((tau[0] - m[i])/dp[i]), 0, 1, 1, 0);
      ftaunovo[i] = pnorm(((taunovo[0] - m[i])/dp[i]), 0, 1, 1, 0);
    }
    else {
      for(j = 1; j < (*tame-1); j++) {
        if(y[i] == esc[j]) {
	  ftau[i]  = pnorm(((tau[j] - m[i])/dp[i]), 0, 1, 1, 0) - pnorm(((tau[j-1] - m[i])/dp[i]), 0, 1, 1, 0);
	  ftaunovo[i] = pnorm(((taunovo[j] - m[i])/dp[i]), 0, 1, 1, 0) - pnorm(((taunovo[j-1] - m[i])/dp[i]), 0, 1, 1, 0);
    }
   }
  }
 }
}  

  

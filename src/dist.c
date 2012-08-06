/* Esta funcao contem as distribuicoes a serem utilizadas na funcao rigtf
   e na funcao principal */


#include <R.h>
#include <Rmath.h>

void F77_SUB(fseedi)(void){
  GetRNGstate();
}

void F77_SUB(fseedo)(void){
  PutRNGstate();
}


void F77_SUB(runif2)(double *px){
  *px = runif(0.00,1.00);
}

void F77_SUB(rgamma2)(double *px, double *shape, double *scale){
  *px = rgamma(*shape,*scale);
}



! Amostra lambda (GS elemento a elemento)
!
!
   SUBROUTINE gslambda (lambda, nu, L, W, theta,n,m)
   IMPLICIT NONE
   INTEGER :: i,j,n,m
   REAL(8) :: nu, sp, rate
   REAL(8) :: suml, scum, resp
   REAL(8), DIMENSION(n) :: L, lambda
   REAL(8), DIMENSION(m) :: theta
   REAL(8), DIMENSION(n,m) :: W
!
!
   sp = (nu+1)/2
!
   CALL fseedi()
   DO i = 1, n
      suml = 0.0
      DO j = 1, m
         suml = (W(i,j)*theta(j))+suml
      ENDDO
      scum = L(i)-suml
      rate = 1/((nu + scum*scum)/2)
      CALL rgamma2(resp, sp, rate)
      lambda(i) = resp
   ENDDO
   CALL fseedo()
   END SUBROUTINE gslambda

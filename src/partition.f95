! Funcao que particiona a matriz A para estimar os componentes de variancia
! A = matriz A de dimensao m x m
! Re = vetor dos efeitos aleatorios estimados
! Vu = Vetor com a variancia dos efeitos aleatorios
! ru e su sao constantes
! tc = Comprimento do vetor Vu
! il e ic = inicio e fim das matrizes quadradas particionadas a partir de A
!
!
!
!
 SUBROUTINE part(A,m,Re,Vu,ru,su,tc,il,ic)
 IMPLICIT NONE
 INTEGER :: i,j,k
 INTEGER :: tc,m
 INTEGER, DIMENSION(tc) :: il, ic
 REAL(8) :: ru, su, c1, c2, resp
 REAL(8), DIMENSION(tc) :: Vu, sumc
 REAL(8), DIMENSION(m) :: Re, result
 REAL(8), DIMENSION(m,m) :: A
!
!   
!
   CALL FSEEDI()
   DO i = 1, tc
      sumc(i) = 0.0
      DO j = il(i), ic(i)
         result(j) = 0.0
         DO k = il(i), ic(i)
            result(j) = A(j,k)*Re(k) + result(j)
         ENDDO
         sumc(i) = result(j)*Re(j) + sumc(i)
      ENDDO
      c1 = (ic(i)+2*ru)/2
      c2 = 1/((sumc(i)+2*su)/2)
      CALL rgamma2(resp,c1,c2)
      Vu(i) = 1/resp
   ENDDO
   CALL FSEEDO()
END SUBROUTINE part



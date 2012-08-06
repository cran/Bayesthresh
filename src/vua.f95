        SUBROUTINE vua(vu, A, ic, il, m, n, t)
        IMPLICIT NONE
        INTEGER :: i,j,m,n,t
        INTEGER, DIMENSION(t) :: ic,il
	REAL(8), DIMENSION(t) :: vu
	REAL(8), DIMENSION(m,n) :: A
!
!
        DO i = 1, t
           DO j = il(i), ic(i)
              A(:,j) = A(:,j)*vu(i)
           ENDDO
        ENDDO
        END SUBROUTINE vua


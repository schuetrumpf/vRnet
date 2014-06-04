! Simple example of use of HSL_MA48
PROGRAM MAIN
   USE HSL_ZD11_DOUBLE
   USE HSL_MA48_DOUBLE
   IMPLICIT NONE
   TYPE(ZD11_TYPE) MATRIX
   TYPE(MA48_CONTROL) CONTROL
   TYPE(MA48_AINFO) AINFO
   TYPE(MA48_FINFO) FINFO
   TYPE(MA48_SINFO) SINFO
   TYPE(MA48_FACTORS) FACTORS

   DOUBLE PRECISION, ALLOCATABLE :: B(:),X(:)
   DOUBLE PRECISION RES(2),ERR
   INTEGER I,INFO,FAST,M,N,NE

! Read matrix order and number of entries

      READ (5,*) M,N,NE
      MATRIX%M = M
      MATRIX%N = N
      MATRIX%NE = NE

! Allocate arrays of appropriate sizes
      ALLOCATE(MATRIX%VAL(NE), MATRIX%ROW(NE), MATRIX%COL(NE))
      ALLOCATE(B(N),X(N))

! Read matrix and right-hand side
      READ (5,*) (MATRIX%ROW(I),MATRIX%COL(I),MATRIX%VAL(I),I=1,NE)
      READ (5,*) B

! Initialize the structures
      CALL MA48_INITIALIZE(FACTORS,CONTROL)

! Analyse and factorize

      CALL MA48_ANALYSE(MATRIX,FACTORS,CONTROL,AINFO,FINFO)
      IF(AINFO%FLAG .ne. 0) THEN
         WRITE(6,'(A,I2)') &
            ' Failure of MA48_ANALYSE with AINFO%FLAG=', AINFO%FLAG
         STOP
      END IF

! Solve without iterative refinement
      CALL MA48_SOLVE(MATRIX,FACTORS,B,X,CONTROL,SINFO)
      IF(SINFO%FLAG == 0) WRITE(6,'(A,/,(6ES11.3))')  &
         'Solution of first set of equations without refinement is',X

! Read new matrix and right-hand side
      READ (5,*) (MATRIX%VAL(I),I=1,NE)
      READ (5,*) B

! Fast factorize
      CALL MA48_FACTORIZE(MATRIX,FACTORS,CONTROL,FINFO,FAST)
      IF(FINFO%FLAG .ne. 0) THEN
         WRITE(6,'(A,I2)') &
            ' Failure of MA48_FACTORIZE with FINFO%FLAG=', FINFO%FLAG
         STOP
      END IF

! Solve with iterative refinement
      CALL MA48_SOLVE(MATRIX,FACTORS,B,X,CONTROL,SINFO,RESID=RES,ERROR=ERR)
      IF(SINFO%FLAG == 0)THEN
         WRITE(6,'(//A,/,(6ES11.3))')       &
              'Solution of second system with refinement is',X
         WRITE(6,'(A/(6ES11.3))') 'Scaled residual is',RES
         WRITE(6,'(A/(6ES11.3))') 'Estimated error is',ERR
      ENDIF

! Clean up
      DEALLOCATE(MATRIX%VAL, MATRIX%ROW, MATRIX%COL)
      CALL MA48_FINALIZE(FACTORS,CONTROL,INFO)

END PROGRAM MAIN

* COPYRIGHT (c) 1993 Council for the Central Laboratory
*                    of the Research Councils

C Original date 29 Jan 2001
C 29 January 2001. Modified from MC49 to be threadsafe.

C 12th July 2004 Version 1.0.0. Version numbering added.
C 28 February 2008. Version 1.0.1. Comments flowed to column 72.
C 21 September 2009. Version 1.0.2. Minor change to documentation.

      SUBROUTINE MC59AD(ICNTL,NC,NR,NE,IRN,LJCN,JCN,LA,A,LIP,IP,
     &                  LIW,IW,INFO)
C
C To sort the sparsity pattern of a matrix to an ordering by columns.
C There is an option for ordering the entries within each column by
C increasing row indices and an option for checking the user-supplied
C matrix entries for indices which are out-of-range or duplicated.
C
C ICNTL:  INTEGER array of length 10. Intent(IN). Used to specify
C         control parameters for the subroutine.
C ICNTL(1): indicates whether the user-supplied matrix entries are to
C           be checked for duplicates, and out-of-range indices.
C           Note  simple checks are always performed.
C           ICNTL(1) = 0, data checking performed.
C           Otherwise, no data checking.
C ICNTL(2): indicates the ordering requested.
C           ICNTL(2) = 0, input is by rows and columns in arbitrary
C           order and the output is sorted by columns.
C           ICNTL(2) = 1, the output is also row ordered
C           within each column.
C           ICNTL(2) = 2, the input is already ordered by
C           columns and is to be row ordered within each column.
C           Values outside the range 0 to 2 are flagged as an error.
C ICNTL(3): indicates whether matrix entries are also being ordered.
C           ICNTL(3) = 0, matrix entries are ordered.
C           Otherwise, only the sparsity pattern is ordered
C           and the array A is not accessed by the routine.
C ICNTL(4): the unit number of the device to
C           which error messages are sent. Error messages
C           can be suppressed by setting ICNTL(4) < 0.
C ICNTL(5): the unit number of the device to
C           which warning messages are sent. Warning
C           messages can be suppressed by setting ICNTL(5) < 0.
C ICNTL(6)  indicates whether matrix symmetric. If unsymmetric, ICNTL(6)
C           must be set to 0.
C           If ICNTL(6) = -1 or 1, symmetric and only the lower
C           triangular part of the reordered matrix is returned.
C           If ICNTL(6) = -2 or 2, Hermitian and only the lower
C           triangular part of the reordered matrix is returned.
C           If error checks are performed (ICNTL(1) = 0)
C           and ICNTL(6)> 1 or 2, the values of duplicate
C           entries are added together; if ICNTL(6) < -1 or -2, the
C           value of the first occurrence of the entry is used.
C ICNTL(7) to ICNTL(10) are not currently accessed by the routine.
C
C NC:      INTEGER variable. Intent(IN). Must be set by the user
C          to the number of columns in the matrix.
C NR:      INTEGER variable. Intent(IN). Must be set by the user
C          to the number of rows in the matrix.
C NE:      INTEGER variable. Intent(IN). Must be set by the user
C          to the number of entries in the matrix.
C IRN: INTEGER array of length NE. Intent (INOUT). Must be set by the
C            user to hold the row indices of the entries in the matrix.
C          If ICNTL(2).NE.2, the entries may be in any order.
C          If ICNTL(2).EQ.2, the entries in column J must be in
C            positions IP(J) to IP(J+1)-1 of IRN. On exit, the row
C            indices are reordered so that the entries of a single
C            column are contiguous with column J preceding column J+1, J
C            = 1, 2, ..., NC-1, with no space between columns.
C          If ICNTL(2).EQ.0, the order within each column is arbitrary;
C            if ICNTL(2) = 1 or 2, the order within each column is by
C            increasing row indices.
C LJCN:    INTEGER variable. Intent(IN). Defines length array
C JCN:     INTEGER array of length LJCN. Intent (INOUT).
C          If ICNTL(2) = 0 or 1, JCN(K) must be set by the user
C          to the column index of the entry
C          whose row index is held in IRN(K), K = 1, 2, ..., NE.
C          On exit, the contents of this array  will have been altered.
C          If ICNTL(2) = 2, the array is not accessed.
C LA:      INTEGER variable. Intent(IN). Defines length of array
C          A.
C A:       is a REAL (DOUBLE PRECISION in the D version, INTEGER in
C          the I version, COMPLEX in the C version,
C          or COMPLEX"*"16 in the Z version) array of length LA.
C          Intent(INOUT).
C          If ICNTL(3).EQ.0, A(K) must be set by the user to
C          hold the value of the entry with row index IRN(K),
C          K = 1, 2, ..., NE. On exit, the array will have been
C          permuted in the same way as the array IRN.
C          If ICNTL(3).NE.0, the array is not accessed.
C LIP:     INTEGER variable. Intent(IN). Defines length of array
C          IP.
C IP:      INTEGER array of length LIP. Intent(INOUT). IP
C          need only be set by the user if ICNTL(2) = 2.
C          In this case, IP(J) holds the position in
C          the array IRN of the first entry in column J, J = 1, 2,
C          ..., NC, and IP(NC+1) is one greater than the number of
C          entries in the matrix.
C          In all cases, the array IP will have this meaning on exit
C          from the subroutine and is altered when ICNTL(2) = 2 only
C          when ICNTL(1) =  0 and there are out-of-range
C          indices or duplicates.
C LIW:     INTEGER variable. Intent(IN). Defines length of array
C          IW.
C IW:      INTEGER array of length LIW. Intent(OUT). Used by the
C          routine as workspace.
C INFO:    INTEGER array of length 10.  Intent(OUT). On exit,
C          a negative value of INFO(1) is used to signal a fatal
C          error in the input data, a positive value of INFO(1)
C          indicates that a warning has been issued, and a
C          zero value is used to indicate a successful call.
C          In cases of error, further information is held in INFO(2).
C          For warnings, further information is
C          provided in INFO(3) to INFO(6).  INFO(7) to INFO(10) are not
C          currently used and are set to zero.
C          Possible nonzero values of INFO(1):
C         -1 -  The restriction ICNTL(2) = 0, 1, or 2 violated.
C               Value of ICNTL(2) is given by INFO(2).
C         -2 -  NC.LE.0. Value of NC is given by INFO(2).
C         -3 -  Error in NR. Value of NR is given by INFO(2).
C         -4 -  NE.LE.0. Value of NE is given by INFO(2).
C         -5 -  LJCN too small. Min. value of LJCN is given by INFO(2).
C         -6 -  LA too small. Min. value of LA is given by INFO(2).
C         -7 -  LIW too small. Value of LIW is given by INFO(2).
C         -8 -  LIP too small. Value of LIP is given by INFO(2).
C         -9 -  The entries of IP not monotonic increasing.
C        -10 -  For each I, IRN(I) or JCN(I) out-of-range.
C        -11 -  ICNTL(6) is out of range.
C         +1 -  One or more duplicated entries. One copy of
C               each such entry is kept and, if ICNTL(3) = 0 and
C               ICNTL(6).GE.0, the values of these entries are
C               added together. If  ICNTL(3) = 0 and ICNTL(6).LT.0,
C               the value of the first occurrence of the entry is used.
C               Initially INFO(3) is set to zero. If an entry appears
C               k times, INFO(3) is incremented by k-1 and INFO(6)
C               is set to the revised number of entries in the
C               matrix.
C         +2 - One or more of the entries in IRN out-of-range. These
C               entries are removed by the routine.`INFO(4) is set to
C               the number of entries which were out-of-range and
C               INFO(6) is set to the revised number of entries in the
C               matrix.
C         +4 - One or more of the entries in JCN out-of-range. These
C               entries are removed by the routine. INFO(5) is set to
C               the number of entries which were out-of-range and
C               INFO(6) is set to the revised number of entries in the
C               matrix. Positive values of INFO(1) are summed so that
C               the user can identify all warnings.
C
C     .. Scalar Arguments ..
      INTEGER LA,LIP,LIW,LJCN,NC,NE,NR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER ICNTL(10),IP(LIP),INFO(10),IRN(NE),IW(LIW),JCN(LJCN)
C     ..
C     .. Local Scalars ..
      INTEGER I,ICNTL1,ICNTL2,ICNTL3,ICNTL6,LAA
      INTEGER IDUP,IOUT,IUP,JOUT,LP,MP,KNE,PART
      LOGICAL LCHECK
C     ..
C     .. External Subroutines ..
      EXTERNAL MC59BD,MC59CD,MC59DD,MC59ED,MC59FD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX
C     ..
C     .. Executable Statements ..

C Initialise
      DO 10 I = 1,10
         INFO(I) = 0
   10 CONTINUE

      ICNTL1 = ICNTL(1)
      ICNTL2 = ICNTL(2)
      ICNTL3 = ICNTL(3)
      ICNTL6 = ICNTL(6)
      LCHECK = (ICNTL1.EQ.0)
C Streams for errors/warnings
      LP = ICNTL(4)
      MP = ICNTL(5)

C  Check the input data
      IF (ICNTL2.GT.2 .OR. ICNTL2.LT.0) THEN
         INFO(1) = -1
         INFO(2) = ICNTL2
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9010) ICNTL2
         END IF
         GO TO 70
      END IF

      IF (ICNTL6.GT.2 .OR. ICNTL6.LT.-2) THEN
         INFO(1) = -11
         INFO(2) = ICNTL6
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9150) ICNTL6
         END IF
         GO TO 70
      END IF
C For real matrices, symmetric = Hermitian so only
C have to distinguish between unsymmetric (ICNTL6 = 0) and
C symmetric (ICNTL6.ne.0)

      IF (NC.LT.1) THEN
        INFO(1) = -2
        INFO(2) = NC
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9020) NC
        END IF
        GO TO 70
      END IF

      IF (NR.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9030) NR
        END IF
        GO TO 70
      END IF

      IF (ICNTL6.NE.0 .AND. NR.NE.NC) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9035) NC,NR
        END IF
        GO TO 70
      END IF

      IF (NE.LT.1) THEN
        INFO(1) = -4
        INFO(2) = NE
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9040) NE
        END IF
        GO TO 70
      END IF

      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.1) THEN
        IF (LJCN.LT.NE) THEN
          INFO(1) = -5
          INFO(2) = NE
        END IF
      ELSE
        IF (LJCN.LT.1) THEN
          INFO(1) = -5
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-5) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9050) LJCN,INFO(2)
         END IF
         GO TO 70
      END IF

      IF (ICNTL3.EQ.0) THEN
        IF (LA.LT.NE) THEN
          INFO(1) = -6
          INFO(2) = NE
        END IF
      ELSE
        IF (LA.LT.1) THEN
          INFO(1) = -6
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-6) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9060) LA,INFO(2)
         END IF
         GO TO 70
      END IF

      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.2) THEN
        IF (LIP.LT.NC+1) THEN
          INFO(1) = -7
          INFO(2) = NC+1
        END IF
      ELSE IF (LIP.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -7
        INFO(2) = MAX(NR,NC)+1
      END IF
      IF (INFO(1).EQ.-7) THEN
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9065) LIP,INFO(2)
        END IF
        GO TO 70
      END IF

C Check workspace sufficient
      IF (LIW.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -8
        INFO(2) = MAX(NR,NC)+1
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9070) LIW,INFO(2)
        END IF
        GO TO 70
      END IF

      LAA = NE
      IF (ICNTL3.NE.0) LAA = 1
C Initialise counts of number of out-of-range entries and duplicates
      IOUT = 0
      JOUT = 0
      IDUP = 0
      IUP = 0

C PART is used by MC59BD to indicate if upper or lower or
C all of matrix is required.
C PART =  0 : unsymmetric case, whole matrix wanted
C PART =  1 : symmetric case, lower triangular part of matrix wanted
C PART = -1 : symmetric case, upper triangular part of matrix wanted
      PART = 0
      IF (ICNTL6.NE.0) PART = 1

      IF (ICNTL2.EQ.0) THEN

C Order directly by columns
C On exit from MC59BD, KNE holds number of entries in matrix
C after removal of out-of-range entries. If no data checking, KNE = NE.
        CALL MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LAA,A,IP,IW,
     +              IOUT,JOUT,KNE)
C Return if ALL entries out-of-range.
        IF (KNE.EQ.0) GO TO 50

C Check for duplicates
        IF (LCHECK) CALL MC59ED(NC,NR,NE,IRN,LIP,IP,LAA,A,IW,IDUP,
     &                          KNE,ICNTL6)

      ELSE IF (ICNTL2.EQ.1) THEN

C First order by rows.
C Interchanged roles of IRN and JCN, so set PART = -1
C if matrix is symmetric case
        IF (ICNTL6.NE.0) PART = -1
        CALL MC59BD(LCHECK,PART,NR,NC,NE,JCN,IRN,LAA,A,IW,IP,
     +              JOUT,IOUT,KNE)
C Return if ALL entries out-of-range.
        IF (KNE.EQ.0) GO TO 50

C At this point, JCN and IW hold column indices and row pointers
C Optionally, check for duplicates.
        IF (LCHECK) CALL MC59ED(NR,NC,NE,JCN,NR+1,IW,LAA,A,IP,
     &                          IDUP,KNE,ICNTL6)

C Now order by columns and by rows within each column
        CALL MC59CD(NC,NR,KNE,IRN,JCN,LAA,A,IP,IW)

      ELSE IF (ICNTL2.EQ.2) THEN
C Input is using IP, IRN.
C Optionally check for duplicates and remove out-of-range entries
        IF (LCHECK) THEN
          CALL MC59FD(NC,NR,NE,IRN,NC+1,IP,LAA,A,LIW,IW,IDUP,
     +                IOUT,IUP,KNE,ICNTL6,INFO)
C Return if IP not monotonic.
          IF (INFO(1).EQ.-9) GO TO 40
C Return if ALL entries out-of-range.
          IF (KNE.EQ.0) GO TO 50
        ELSE
           KNE = NE
        END IF

C  Order by rows within each column
        CALL MC59DD(NC,KNE,IRN,IP,LAA,A)

      END IF

      INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(6) = KNE
      INFO(7) = IUP

C Set warning flag if out-of-range /duplicates found
      IF (IDUP.GT.0) INFO(1) = INFO(1) + 1
      IF (IOUT.GT.0) INFO(1) = INFO(1) + 2
      IF (JOUT.GT.0) INFO(1) = INFO(1) + 4
      IF (INFO(1).GT.0 .AND. MP.GT.0) THEN
        WRITE (MP,FMT=9080) INFO(1)
        IF (IOUT.GT.0) WRITE (MP,FMT=9090) IOUT
        IF (JOUT.GT.0) WRITE (MP,FMT=9110) JOUT
        IF (IDUP.GT.0) WRITE (MP,FMT=9100) IDUP
        IF (IUP.GT.0)  WRITE (MP,FMT=9130) IUP
      END IF
      GO TO 70

   40 INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(7) = IUP
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9140)
      END IF
      GO TO 70

   50 INFO(1) = -10
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(2) = IOUT + JOUT
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9120)
      END IF
   70 RETURN

 9000 FORMAT (/,' *** Error return from MC59AD *** INFO(1) = ',I3)
 9010 FORMAT (1X,'ICNTL(2) = ',I2,' is out of range')
 9020 FORMAT (1X,'NC = ',I6,' is out of range')
 9030 FORMAT (1X,'NR = ',I6,' is out of range')
 9035 FORMAT (1X,'Symmetric case. NC = ',I6,' but NR = ',I6)
 9040 FORMAT (1X,'NE = ',I10,' is out of range')
 9050 FORMAT (1X,'Increase LJCN from ',I10,' to at least ',I10)
 9060 FORMAT (1X,'Increase LA from ',I10,' to at least ',I10)
 9065 FORMAT (1X,'Increase LIP from ',I8,' to at least ',I10)
 9070 FORMAT (1X,'Increase LIW from ',I8,' to at least ',I10)
 9080 FORMAT (/,' *** Warning message from MC59AD *** INFO(1) = ',I3)
 9090 FORMAT (1X,I8,' entries in IRN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9100 FORMAT (1X,I8,' duplicate entries were supplied by the user')
 9110 FORMAT (1X,I8,' entries in JCN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9120 FORMAT (1X,'All entries out of range')
 9130 FORMAT (1X,I8,' of these entries were in the upper triangular ',
     +       /,'       part of matrix')
 9140 FORMAT (1X,'Entries in IP are not monotonic increasing')
 9150 FORMAT (1X,'ICNTL(6) = ',I2,' is out of range')
      END
C***********************************************************************
      SUBROUTINE MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LA,A,IP,IW,IOUT,
     +                  JOUT,KNE)
C
C   To sort a sparse matrix from arbitrary order to
C   column order, unordered within each column. Optionally
C   checks for out-of-range entries in IRN,JCN.
C
C LCHECK - logical variable. Intent(IN). If true, check
C          for out-of-range indices.
C PART -   integer variable. Intent(IN)
C PART =  0 : unsymmetric case, whole matrix wanted
C PART =  1 : symmetric case, lower triangular part of matrix wanted
C             (ie IRN(K) .ge. JCN(K) on exit)
C PART = -1 : symmetric case, upper triangular part of matrix wanted
C             (ie IRN(K) .le. JCN(K) on exit)
C   NC - integer variable. Intent(IN)
C      - on entry must be set to the number of columns in the matrix
C   NR - integer variable. Intent(IN)
C      - on entry must be set to the number of rows in the matrix
C   NE - integer variable. Intent(IN)
C      - on entry, must be set to the number of nonzeros in the matrix
C  IRN - integer array of length NE. Intent(INOUT)
C      - on entry set to contain the row indices of the nonzeros
C        in arbitrary order.
C      - on exit, the entries in IRN are reordered so that the row
C        indices for column 1 precede those for column 2 and so on,
C        but the order within columns is arbitrary.
C  JCN - integer array of length NE. Intent(INOUT)
C      - on entry set to contain the column indices of the nonzeros
C      - JCN(K) must be the column index of
C        the entry in IRN(K)
C      - on exit, JCN(K) is the column index for the entry with
C        row index IRN(K) (K=1,...,NE).
C  LA  - integer variable which defines the length of the array A.
C        Intent(IN)
C   A  - real (double precision/complex/complex*16) array of length LA
C        Intent(INOUT)
C      - if LA > 1, the array must be of length NE, and A(K)
C        must be set to the value of the entry in (IRN(K), JCN(K));
C        on exit A is reordered in the same way as IRN
C      - if LA = 1, the array is not accessed
C  IP  - integer array of length NC+1. Intent(INOUT)
C      - not set on entry
C      - on exit, IP(J) contains the position in IRN (and A) of the
C        first entry in column J (J=1,...,NC)
C      - IP(NC+1) is set to NE+1
C  IW  - integer array of length NC+1.  Intent(INOUT)
C      - the array is used as workspace
C      - on exit IW(I) = IP(I) (so IW(I) points to the beginning
C        of column I).
C IOUT - integer variable. Intent(OUT). On exit, holds number
C        of entries in IRN found to be out-of-range
C JOUT - integer variable. Intent(OUT). On exit, holds number
C        of entries in JCN found to be out-of-range
C  KNE - integer variable. Intent(OUT). On exit, holds number
C        of entries in matrix after removal of out-of-range entries.
C        If no data checking, KNE = NE.

C     .. Scalar Arguments ..
      INTEGER LA,NC,NE,NR,IOUT,JOUT,KNE,PART
      LOGICAL LCHECK
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NC+1),JCN(NE)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,JCE,JCEP,K,L,LOC
C     ..
C     .. Executable Statements ..

C Initialise IW
      DO 10 J = 1,NC + 1
        IW(J) = 0
   10 CONTINUE

      KNE = 0
      IOUT = 0
      JOUT = 0
C Count the number of entries in each column and store in IW.
C We also allow checks for out-of-range indices
      IF (LCHECK) THEN
C Check data.
C Treat case of pattern only separately.
        IF (LA.GT.1) THEN
          IF (PART.EQ.0) THEN
C Unsymmetric
            DO 20 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
C IRN out-of-range. Is JCN also out-of-range?
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                A(KNE) = A(K)
                IW(J) = IW(J) + 1
              END IF
   20       CONTINUE
          ELSE IF (PART.EQ.1) THEN
C Symmetric, lower triangle
            DO 21 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
C IRN out-of-range. Is JCN also out-of-range?
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
C Lower triangle ... swap if necessary
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   21       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
C Symmetric, upper triangle
            DO 22 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
C IRN out-of-range. Is JCN also out-of-range?
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
C Upper triangle ... swap if necessary
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   22       CONTINUE
          END IF
        ELSE
C Pattern only
          IF (PART.EQ.0) THEN
            DO 25 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                IW(J) = IW(J) + 1
              END IF
   25       CONTINUE
          ELSE IF (PART.EQ.1) THEN
            DO 26 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
C Lower triangle ... swap if necessary
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   26       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
            DO 27 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
C Upper triangle ... swap if necessary
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   27       CONTINUE
          END IF
        END IF
C Return if ALL entries out-of-range.
        IF (KNE.EQ.0) GO TO 130

      ELSE

C No checks
        KNE = NE
        IF (PART.EQ.0) THEN
          DO 30 K = 1,NE
            J = JCN(K)
            IW(J) = IW(J) + 1
   30     CONTINUE
        ELSE IF (PART.EQ.1) THEN
          DO 35 K = 1,NE
            I = IRN(K)
            J = JCN(K)
C Lower triangle ... swap if necessary
            IF (I.LT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   35     CONTINUE
        ELSE IF (PART.EQ.-1) THEN
          DO 36 K = 1,NE
            I = IRN(K)
            J = JCN(K)
C Upper triangle ... swap if necessary
            IF (I.GT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   36     CONTINUE
        END IF
      END IF

C KNE is now the number of nonzero entries in matrix.

C Put into IP and IW the positions where each column
C would begin in a compressed collection with the columns
C in natural order.

      IP(1) = 1
      DO 37 J = 2,NC + 1
        IP(J) = IW(J-1) + IP(J-1)
        IW(J-1) = IP(J-1)
   37 CONTINUE

C Reorder the elements into column order.
C Fill in each column from the front, and as a new entry is placed
C in column K increase the pointer IW(K) by one.

      IF (LA.EQ.1) THEN
C Pattern only
        DO 70 L = 1,NC
          DO 60 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            DO 40 J = 1,NE
              IF (JCE.EQ.L) GO TO 50
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
   40       CONTINUE
   50       JCN(K) = JCE
            IRN(K) = ICE
   60     CONTINUE
   70   CONTINUE
      ELSE

        DO 120 L = 1,NC
          DO 110 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            ACE = A(K)
            DO 90 J = 1,NE
              IF (JCE.EQ.L) GO TO 100
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
              ACEP = A(LOC)
              A(LOC) = ACE
              ACE = ACEP
   90       CONTINUE
  100       JCN(K) = JCE
            IRN(K) = ICE
            A(K) = ACE
  110     CONTINUE
  120   CONTINUE
      END IF

  130 CONTINUE

      RETURN
      END
C
C**********************************************************
      SUBROUTINE MC59CD(NC,NR,NE,IRN,JCN,LA,A,IP,IW)
C
C   To sort a sparse matrix stored by rows,
C   unordered within each row, to ordering by columns, with
C   ordering by rows within each column.
C
C   NC - integer variable. Intent(IN)
C      - on entry must be set to the number of columns in the matrix
C   NR - integer variable. Intent(IN)
C      - on entry must be set to the number of rows in the matrix
C  NE - integer variable. Intent(IN)
C      - on entry, must be set to the number of nonzeros in the matrix
C  IRN - integer array of length NE. Intent(OUT).
C      - not set on entry.
C      - on exit,  IRN holds row indices with the row
C        indices for column 1 preceding those for column 2 and so on,
C        with ordering by rows within each column.
C  JCN - integer array of length NE. Intent(INOUT)
C      - on entry set to contain the column indices of the nonzeros
C        with indices for column 1 preceding those for column 2
C        and so on, with the order within columns is arbitrary.
C      - on exit, contents destroyed.
C  LA  - integer variable which defines the length of the array A.
C        Intent(IN)
C   A  - real (double precision/complex/complex*16) array of length LA
C        Intent(INOUT)
C      - if LA > 1, the array must be of length NE, and A(K)
C        must be set to the value of the entry in JCN(K);
C        on exit A, A(K) holds the value of the entry in IRN(K).
C      - if LA = 1, the array is not accessed
C  IP  - integer array of length NC+1. Intent(INOUT)
C      - not set on entry
C      - on exit, IP(J) contains the position in IRN (and A) of the
C        first entry in column J (J=1,...,NC)
C      - IP(NC+1) is set to NE+1
C  IW  - integer array of length NR+1.  Intent(IN)
C      - on entry, must be set on entry so that IW(J) points to the
C        position in JCN of the first entry in row J, J=1,...,NR, and
C        IW(NR+1) must be set to NE+1
C
C     .. Scalar Arguments ..
      INTEGER LA,NC,NE,NR
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NR+1),JCN(NE)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,J1,J2,K,L,LOC,LOCP
C     ..
C     .. Executable Statements ..

C  Count the number of entries in each column

      DO 10 J = 1,NC
        IP(J) = 0
   10 CONTINUE

      IF (LA.GT.1) THEN

        DO 20 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
          IRN(K) = JCN(K)
   20   CONTINUE
        IP(NC+1) = NE + 1

C  Set IP so that IP(I) points to the first entry in column I+1

        IP(1) = IP(1) + 1
        DO 30 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
   30   CONTINUE

        DO 50 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 40 J = J1,J2
            K = IRN(J)
            L = IP(K) - 1
            JCN(J) = L
            IRN(J) = I
            IP(K) = L
   40     CONTINUE
   50   CONTINUE
        IP(NC+1) = NE + 1
        DO 70 J = 1,NE
          LOC = JCN(J)
          IF (LOC.EQ.0) GO TO 70
          ICE = IRN(J)
          ACE = A(J)
          JCN(J) = 0
          DO 60 K = 1,NE
            LOCP = JCN(LOC)
            ICEP = IRN(LOC)
            ACEP = A(LOC)
            JCN(LOC) = 0
            IRN(LOC) = ICE
            A(LOC) = ACE
            IF (LOCP.EQ.0) GO TO 70
            ICE = ICEP
            ACE = ACEP
            LOC = LOCP
   60     CONTINUE
   70   CONTINUE
      ELSE

C Pattern only

C  Count the number of entries in each column

        DO 90 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
   90   CONTINUE
        IP(NC+1) = NE + 1

C  Set IP so that IP(I) points to the first entry in column I+1

        IP(1) = IP(1) + 1
        DO 100 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
  100   CONTINUE

        DO 120 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 110 J = J1,J2
            K = JCN(J)
            L = IP(K) - 1
            IRN(L) = I
            IP(K) = L
  110     CONTINUE
  120   CONTINUE

      END IF

      RETURN
      END

C**********************************************************

      SUBROUTINE MC59DD(NC,NE,IRN,IP,LA,A)
C
C To sort from arbitrary order within each column to order
C by increasing row index. Note: this is taken from MC20B/BD.
C
C   NC - integer variable. Intent(IN)
C      - on entry must be set to the number of columns in the matrix
C   NE - integer variable. Intent(IN)
C      - on entry, must be set to the number of nonzeros in the matrix
C  IRN - integer array of length NE. Intent(INOUT)
C      - on entry set to contain the row indices of the nonzeros
C        ordered so that the row
C        indices for column 1 precede those for column 2 and so on,
C        but the order within columns is arbitrary.
C        On exit, the order within each column is by increasing
C        row indices.
C   LA - integer variable which defines the length of the array A.
C        Intent(IN)
C    A - real (double precision/complex/complex*16) array of length LA
C        Intent(INOUT)
C      - if LA > 1, the array must be of length NE, and A(K)
C        must be set to the value of the entry in IRN(K);
C        on exit A is reordered in the same way as IRN
C      - if LA = 1, the array is not accessed
C  IP  - integer array of length NC. Intent(IN)
C      - on entry, IP(J) contains the position in IRN (and A) of the
C        first entry in column J (J=1,...,NC)
C     . .
C     .. Scalar Arguments ..
      INTEGER LA,NC,NE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(NC)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ACE
      INTEGER ICE,IK,J,JJ,K,KDUMMY,KLO,KMAX,KOR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Executable Statements ..

C Jump if pattern only.
      IF (LA.GT.1) THEN
        KMAX = NE
        DO 50 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 40
          KOR = KMAX
          DO 30 KDUMMY = KLO,KMAX
C Items KOR, KOR+1, .... ,KMAX are in order
            ACE = A(KOR-1)
            ICE = IRN(KOR-1)
            DO 10 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 20
              IRN(K-1) = IK
              A(K-1) = A(K)
   10       CONTINUE
            K = KMAX + 1
   20       IRN(K-1) = ICE
            A(K-1) = ACE
            KOR = KOR - 1
   30     CONTINUE
C Next column
   40     KMAX = KLO - 2
   50   CONTINUE
      ELSE

C Pattern only.
        KMAX = NE
        DO 150 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 140
          KOR = KMAX
          DO 130 KDUMMY = KLO,KMAX
C Items KOR, KOR+1, .... ,KMAX are in order
            ICE = IRN(KOR-1)
            DO 110 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 120
              IRN(K-1) = IK
  110       CONTINUE
            K = KMAX + 1
  120       IRN(K-1) = ICE
            KOR = KOR - 1
  130     CONTINUE
C Next column
  140     KMAX = KLO - 2
  150   CONTINUE
      END IF
      END
C***********************************************************************

      SUBROUTINE MC59ED(NC,NR,NE,IRN,LIP,IP,LA,A,IW,IDUP,KNE,ICNTL6)

C Checks IRN for duplicate entries.
C On exit, IDUP holds number of duplicates found and KNE is number
C of entries in matrix after removal of duplicates
C     . .
C     .. Scalar Arguments ..
      INTEGER ICNTL6,IDUP,KNE,LIP,LA,NC,NR,NE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(LIP),IW(NR)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,K,KSTART,KSTOP,NZJ

      IDUP = 0
      KNE = 0
C Initialise IW
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE

      KSTART = IP(1)
      IF (LA.GT.1) THEN
C Matrix entries considered
        NZJ = 0
        DO 30 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
C We have a duplicate in column J
              IDUP = IDUP + 1
C If requested, sum duplicates
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE

      ELSE

C Pattern only
        DO 50 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 40 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
C  We have a duplicate in column J
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF

      RETURN
      END
C***********************************************************************

      SUBROUTINE MC59FD(NC,NR,NE,IRN,LIP,IP,LA,A,LIW,IW,IDUP,IOUT,
     +                  IUP,KNE,ICNTL6,INFO)

C Checks IRN for duplicate and out-of-range entries.
C For symmetric matrix, also checks NO entries lie in upper triangle.
C Also checks IP is monotonic.
C On exit:
C IDUP holds number of duplicates found
C IOUT holds number of out-of-range entries
C For symmetric matrix, IUP holds number of entries in upper
C triangular part.
C KNE holds number of entries in matrix after removal of
C out-of-range and duplicate entries.
C Note: this is similar to MC59ED except it also checks IP is
C monotonic and removes out-of-range entries in IRN.
C     . .
C     .. Scalar Arguments ..
      INTEGER ICNTL6,IDUP,IOUT,IUP,KNE,LA,LIP,LIW,NC,NR,NE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(LIP),IW(LIW),INFO(2)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,K,KSTART,KSTOP,NZJ,LOWER

      IDUP = 0
      IOUT = 0
      IUP = 0
      KNE = 0
C Initialise IW
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE

      KSTART = IP(1)
      LOWER = 1
      IF (LA.GT.1) THEN
        NZJ = 0
        DO 30 J = 1,NC
C In symmetric case, entries out-of-range if they lie
C in upper triangular part.
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
C Check for out-of-range
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
C In symmetric case, check if entry is out-of-range because
C it lies in upper triangular part.
              IF (ICNTL6.NE.0 .AND. I.LT.J) IUP = IUP + 1
            ELSE IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
C  We have a duplicate in column J
              IDUP = IDUP + 1
C If requested, sum duplicates
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE

      ELSE

C Pattern only
        DO 50 J = 1,NC
C In symmetric case, entries out-of-range if lie
C in upper triangular part.
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO  40 K = KSTART,KSTOP - 1
            I = IRN(K)
C Check for out-of-range
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
              IF (ICNTL6.NE.0 .AND. I.GT.1) IUP = IUP + 1
            ELSE IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
C  We have a duplicate in column J
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF

      RETURN
      END

C COPYRIGHT (c) 1993 Council for the Central Laboratory
*                    of the Research Councils
C Original date 20 May 1993
C 24 September 1993 Some IVDEP comments added for speed on the Cray,
C     some minor bugs fixed, default changed to BLAS 3 with block size
C     32.
C 6 December 1993. Minor bug fixed re threshold test for pivots.
C 4/10/95. IQ in MA50BD made assumed size
C 1/11/95. IQ in MA50CD made assumed size
C 1/11/95. DTRSV not called for zero-sized array.
C 14/2/96. NP initialized to 0.
C 13/11/97 INFO(4) and INFO(6) in MA50AD made to reflect the situation
C          at the point of failure in the case of insufficient storage.
C 17/3/98  In MA50AD, copy a row forward if there is space at its front,
C          rather than put new entry at front. Makes the result
C          repeatable as LA is altered.

C 12th July 2004 Version 1.0.0. Version numbering added.
C 29 November 2006 Version 2.0.0. Sizes of CNTL, ICNTL, INFO, RINFO
C          increased and ICNTL(8) added.

      SUBROUTINE MA50AD(M,N,NE,LA,A,IRN,JCN,IQ,CNTL,ICNTL,IP,NP,JFIRST,
     +                  LENR,LASTR,NEXTR,IW,IFIRST,LENC,LASTC,NEXTC,
     +                  INFO,RINFO)

C MA50A/AD chooses a pivot sequence using a Markowitz criterion with
C     threshold pivoting.

C If  the user requires a more convenient data interface then the MA48
C     package should be used. The MA48 subroutines call the MA50
C     subroutines after checking the user's input data and optionally
C     permute the matrix to block triangular form.

      INTEGER M,N,NE,LA
      DOUBLE PRECISION A(LA)
      DOUBLE PRECISION CNTL(10)
      INTEGER IRN(LA),JCN(LA),IQ(N)
      INTEGER ICNTL(20),IP(M),NP,JFIRST(M),LENR(M),LASTR(M),NEXTR(M),
     +        IW(M),IFIRST(N),LENC(N),LASTC(N),NEXTC(N),INFO(15)
      DOUBLE PRECISION RINFO(10)

C M is an integer variable that must be set to the number of rows.
C      It is not altered by the subroutine.
C N is an integer variable that must be set to the number of columns.
C      It is not altered by the subroutine.
C NE is an integer variable that must be set to the number of entries
C      in the input matrix. It is not altered by the subroutine.
C LA is an integer variable that must be set to the size of A, IRN, and
C      JCN. It is not altered by the subroutine.
C A is an array that holds the input matrix on entry and is used as
C      workspace.
C IRN  is an integer array.  Entries 1 to NE must be set to the
C      row indices of the corresponding entries in A.  IRN is used
C      as workspace and holds the row indices of the reduced matrix.
C JCN  is an integer array that need not be set by the user. It is
C      used to hold the column indices of entries in the reduced
C      matrix.
C IQ is an integer array of length N. On entry, it holds pointers
C      to column starts. During execution, IQ(j) holds the position of
C      the start of column j of the reduced matrix or -IQ(j) holds the
C      column index in the permuted matrix of column j. On exit, IQ(j)
C      holds the index of the column that is in position j of the
C      permuted matrix.
C CNTL must be set by the user as follows and is not altered.
C     CNTL(1)  Full matrix processing will be used if the density of
C       the reduced matrix is MIN(CNTL(1),1.0) or more.
C     CNTL(2) determines the balance between pivoting for sparsity and
C       for stability, values near zero emphasizing sparsity and values
C       near one emphasizing stability. Each pivot must have absolute
C       value at least CNTL(2) times the greatest absolute value in the
C       same column of the reduced matrix.
C     CNTL(3) If this is set to a positive value, any entry of the
C       reduced matrix whose modulus is less than CNTL(3) will be
C       dropped.
C     CNTL(4)  Any entry of the reduced matrix whose modulus is less
C       than or equal to CNTL(4) will be regarded as zero from the
C        point of view of rank.
C ICNTL must be set by the user as follows and is not altered.
C     ICNTL(1)  must be set to the stream number for error messages.
C       A value less than 1 suppresses output.
C     ICNTL(2) must be set to the stream number for diagnostic output.
C       A value less than 1 suppresses output.
C     ICNTL(3) must be set to control the amount of output:
C       0 None.
C       1 Error messages only.
C       2 Error and warning messages.
C       3 As 2, plus scalar parameters and a few entries of array
C         parameters on entry and exit.
C       4 As 3, plus all parameters on entry and exit.
C     ICNTL(4) If set to a positive value, the pivot search is limited
C       to ICNTL(4) columns (Zlatev strategy). This may result in
C       different fill-in and execution time. If ICNTL(4) is positive,
C       the workspace arrays LASTR and NEXTR are not referenced.
C     ICNTL(5) The block size to be used for full-matrix processing.
C     ICNTL(6) The last ICNTL(6) columns of A must be the last
C       ICNTL(6) columns of the permuted matrix. A value outside the
C       range 1 to N-1 is treated as zero.
C     ICNTL(7) If given the value 1, pivots are limited to
C       the main diagonal, which may lead to a premature switch to full
C       processing if no suitable diagonal entries are available.
C       If given the value 2, IFIRST must be set so that IFIRST(i) is
C       the column in position i of the permuted matrix and IP must
C       be set so that IP(i) < IP(j) if row i is recommended to
C       precede row j in the pivot sequence.
C IP is an integer array of length M that need not be set on entry
C      unless ICNTL(7)=2 (see ICNTL(7) for details of this case).
C      During execution, IP(i) holds the position of the start of row i
C      of the reduced matrix or -IP(i) holds the row index in the
C      permuted matrix of row i. Before exit, IP(i) is made positive.
C NP is an integer variable. It need not be set on entry. On exit,
C     it will be set to the number of columns to be processed in
C     packed storage.
C JFIRST is an integer workarray of length M. JFIRST(i) is the
C      first column of the reduced matrix to have i entries or is
C      zero if no column has i entries.
C LENR is an integer workarray of length M that is used to hold the
C      numbers of entries in the rows of the reduced matrix.
C LASTR is an integer workarray of length M, used only if ICNTL(4) = 0.
C      For rows in the reduced matrix, LASTR(i) indicates the previous
C      row to i with the same number of entries. LASTR(i) is zero if
C      no such row exists.
C NEXTR is an integer workarray of length M, used only if ICNTL(4) = 0
C      or ICNTL(7)=2. If ICNTL(4)=0, for rows in the reduced matrix,
C      NEXTR(i) indicates the next row to i with the same number of
C      entries; and if row i is the last in the chain, NEXTR is
C      equal to zero. If ICNTL(7)=2, NEXTR is a copy of the value of
C      IP on entry.
C IW is an integer array of length M used as workspace and is used to
C     assist the detection of duplicate entries and the sparse SAXPY
C     operations. It is reset to zero each time round the main loop.
C IFIRST is an integer array of length N, used only if ICNTL(4) = 0
C      or ICNTL(7)=2. If ICNTL(4) = 0, it is a workarray; IFIRST(i)
C      points to the first row of the reduced matrix to have i entries
C      or is zero if no row has i entries. If ICNTL(7)=2, IFIRST
C      must be set on entry (see ICNTL(7) for details of this case).
C LENC is an integer workarray of length N that is used to hold
C      the numbers of entries in the columns of the reduced matrix.
C LASTC is an integer workarray of length N.  For columns in the reduced
C      matrix, LASTC(j) indicates the previous column to j with the same
C      number of entries.  If column j is the first in the chain,
C      LASTC(j) is equal to zero.
C NEXTC is an integer workarray of length N.  For columns in the reduced
C      matrix, NEXTC(j) indicates the next column to j with the same
C      number of entries.  If column j is the last column in the chain,
C      NEXTC(j) is zero.
C INFO need not be set on entry. On exit, it holds the following:
C    INFO(1):
C       0  Successful entry.
C      -1  M < 1 or N < 1.
C      -2  NE < 1.
C      -3  Insufficient space.
C      -4  Duplicated entries.
C      -5  Faulty column permutation in IFIRST when ICNTL(7)=2.
C      -6  ICNTL(4) not equal to 1 when ICNTL(7)=2.
C      +1  Rank deficient.
C      +2  Premature switch to full processing because of failure to
C          find a stable diagonal pivot (ICNTL(7)>=1 case only).
C      +3  Both of these warnings.
C    INFO(2) Number of compresses of the arrays.
C    INFO(3) Minimum LA recommended to analyse matrix.
C    INFO(4) Minimum LFACT required to factorize matrix.
C    INFO(5) Upper bound on the rank of the matrix.
C    INFO(6) Number of entries dropped from the data structure.
C    INFO(7) Number of rows processed in full storage.
C RINFO need not be set on entry. On exit, RINFO(1) holds the number of
C    floating-point operations needed for the factorization.

      INTEGER IDAMAX
      EXTERNAL IDAMAX
      EXTERNAL MA50DD
      INTRINSIC ABS,MAX,MIN

      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0D0,ONE=1.0D0)

      DOUBLE PRECISION ALEN,AMULT,ANEW,ASW,AU,COST,CPIV
      INTEGER DISPC,DISPR,EYE,I,IDROP,IDUMMY,IEND,IFILL,IFIR,II,IJ,
     +        IJPOS,IOP,IPIV,IPOS,ISRCH,I1,I2,J,JBEG,JEND,JJ,JLAST,
     +        JMORE,JNEW,JPIV,JPOS,J1,J2,L,LC,LEN,LENPIV,LP,LR
      DOUBLE PRECISION MAXENT
      INTEGER MINC,MORD,MP,MSRCH,NC,NDROP,NEFACT,NEPR,NERED,NE1,NORD,
     +        NORD1,NR,NULLC,NULLI,NULLJ,NULLR,PIVBEG,PIVCOL,PIVEND,
     +        PIVOT
      DOUBLE PRECISION PIVR,PIVRAT,U

C ALEN Real(LEN-1).
C AMULT Temporary variable used to store current multiplier.
C ANEW Temporary variable used to store value of fill-in.
C ASW Temporary variable used when swopping two real quantities.
C AU Temporary variable used in threshold test.
C COST Markowitz cost of current potential pivot.
C CPIV Markowitz cost of best pivot so far found.
C DISPC is the first free location in the column file.
C DISPR is the first free location in the row file.
C EYE Running relative position when processing pivot row.
C I Temporary variable holding row number. Also used as index in DO
C     loops used in initialization of arrays.
C IDROP Temporary variable used to accumulate number of entries dropped.
C IDUMMY DO index not referenced in the loop.
C IEND Position of end of pivot row.
C IFILL is the fill-in to the non-pivot column.
C IFIR Temporary variable holding first entry in chain.
C II Running position for current column.
C IJ Temporary variable holding row/column index.
C IJPOS Position of current pivot in A/IRN.
C IOP holds a running count of the number of rows with entries in both
C     the pivot and the non-pivot column.
C IPIV Row of the pivot.
C IPOS Temporary variable holding position in column file.
C ISRCH Temporary variable holding number of columns searched for pivot.
C I1 Position of the start of the current column.
C I2 Position of the end of the current column.
C J Temporary variable holding column number.
C JBEG Position of beginning of non-pivot column.
C JEND Position of end of non-pivot column.
C JJ Running position for current row.
C JLAST Last column acceptable as pivot.
C JMORE Temporary variable holding number of locations still needed
C     for fill-in in non-pivot column.
C JNEW Position of end of changed non-pivot column.
C JPIV Column of the pivot.
C JPOS Temporary variable holding position in row file.
C J1 Position of the start of the current row.
C J2 Position of the end of the current row.
C L Loop index.
C LC Temporary variable holding previous column in sequence.
C LEN Length of column or row.
C LENPIV Length of pivot column.
C LP Unit for error messages.
C LR Temporary variable holding previous row in sequence.
C MAXENT Temporary variable used to hold value of largest entry in
C    column.
C MINC Minimum number of entries of any row or column of the reduced
C     matrix, or in any column if ICNTL(4) > 0.
C MORD Number of rows ordered, excluding null rows.
C MP Unit for diagnostic messages.
C MSRCH Number of columns to be searched.
C NC Temporary variable holding next column in sequence.
C NDROP Number of entries dropped because of being in a column all of
C   whose entries are smaller than the pivot threshold.
C NEFACT Number of entries in factors.
C NEPR Number of entries in pivot row, excluding the pivot.
C NERED Number of entries in reduced matrix.
C NE1 Temporary variable used to hold number of entries in row/column
C     and to hold temporarily value of MINC.
C NORD Number of columns ordered, excluding null columns beyond JLAST.
C NORD1 Value of NORD at start of step.
C NR Temporary variable holding next row in sequence.
C NULLC Number of structurally zero columns found before any entries
C     dropped for being smaller than CNTL(3).
C NULLR Number of structurally zero rows found before any entries
C     dropped for being smaller than CNTL(3).
C NULLI Number of zero rows found.
C NULLJ Number of zero columns found beyond column JLAST.
C PIVBEG Position of beginning of pivot column.
C PIVCOL Temporary variable holding position in pivot column.
C PIVEND Position of end of pivot column.
C PIVOT Current step in Gaussian elimination.
C PIVR ratio of current pivot candidate to largest in its column.
C PIVRAT ratio of best pivot candidate to largest in its column.
C U Used to hold local copy of CNTL(2), changed if necessary so that it
C    is in range.

      LP = ICNTL(1)
      IF (ICNTL(3).LE.0) LP = 0
      MP = ICNTL(2)
      IF (ICNTL(3).LE.1) MP = 0
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = NE
      INFO(4) = NE
      INFO(5) = 0
      INFO(6) = 0
      INFO(7) = 0
      RINFO(1) = ZERO
      NP = 0

C Make some simple checks
      IF (M.LT.1 .OR. N.LT.1) GO TO 690
      IF (NE.LT.1) GO TO 700
      IF (LA.LT.NE) THEN
         INFO(3) = NE
         GO TO 710
      END IF

C Initial printing
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(/2(A,I6),A,I8,A,I8/A,1P,4E10.2/A,7I4)')
     +     ' Entering MA50AD with M =',M,' N =',N,' NE =',NE,' LA =',LA,
     +     ' CNTL =',(CNTL(I),I=1,4),' ICNTL =',(ICNTL(I),I=1,7)
         IF (N.EQ.1 .OR. ICNTL(3).GT.3) THEN
            DO 10 J = 1,N - 1
               IF (IQ(J).LT.IQ(J+1)) WRITE (MP,
     +             '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',J,
     +             (A(II),IRN(II),II=IQ(J),IQ(J+1)-1)
   10       CONTINUE
            IF (IQ(N).LE.NE) WRITE (MP,'(A,I5,(T13,3(1P,E12.4,I5)))')
     +          ' Column',N, (A(II),IRN(II),II=IQ(N),NE)
         ELSE
            IF (IQ(1).LT.IQ(2)) WRITE (MP,
     +          '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',1,
     +          (A(II),IRN(II),II=IQ(1),IQ(2)-1)
         END IF
         IF (ICNTL(7).EQ.2) THEN
            WRITE (MP,'(A,(T10,10(I7)))') ' IP = ',IP
            WRITE (MP,'(A,(T10,10(I7)))') ' IFIRST = ',IFIRST
         END IF
      END IF

C Initialization of counts etc.
      MINC = 1
      NERED = NE
      U = MIN(CNTL(2),ONE)
      U = MAX(U,ZERO)
      MSRCH = ICNTL(4)
      IF (MSRCH.EQ.0) MSRCH = N
      JLAST = N - ICNTL(6)
      IF (JLAST.LT.1 .OR. JLAST.GT.N) JLAST = N
      NULLI = 0
      NULLJ = 0
      MORD = 0
      NORD = 0
      NDROP = 0
      NEFACT = 0
      DO 20 I = 1,N - 1
         LENC(I) = IQ(I+1) - IQ(I)
   20 CONTINUE
      LENC(N) = NE + 1 - IQ(N)

      IF (CNTL(3).GT.ZERO) THEN
C Drop small entries
         NERED = 0
         DO 40 J = 1,N
            I = IQ(J)
            IQ(J) = NERED + 1
            DO 30 II = I,I + LENC(J) - 1
               IF (ABS(A(II)).GE.CNTL(3)) THEN
                  NERED = NERED + 1
                  A(NERED) = A(II)
                  IRN(NERED) = IRN(II)
               ELSE
                  INFO(6) = INFO(6) + 1
               END IF
   30       CONTINUE
            LENC(J) = NERED + 1 - IQ(J)
   40    CONTINUE
      END IF

      IF (ICNTL(7).EQ.2) THEN
C Column order specified - copy the row ordering array
         DO 50 I = 1,M
            NEXTR(I) = IP(I)
   50    CONTINUE
C Check ICNTL(4)
         IF (ICNTL(4).NE.1) GO TO 740
      END IF

      DISPR = NERED + 1
      DISPC = NERED + 1
C
C Set up row oriented storage.
      DO 60 I = 1,M
         IW(I) = 0
         LENR(I) = 0
         JFIRST(I) = 0
   60 CONTINUE
C Calculate row counts.
      DO 70 II = 1,NERED
         I = IRN(II)
         LENR(I) = LENR(I) + 1
   70 CONTINUE
C Set up row pointers so that IP(i) points to position after end
C     of row i in row file.
      IP(1) = LENR(1) + 1
      DO 80 I = 2,M
         IP(I) = IP(I-1) + LENR(I)
   80 CONTINUE
C Generate row file.
      DO 100 J = 1,N
         I = IQ(J)
         DO 90 II = I,I + LENC(J) - 1
            I = IRN(II)
C Check for duplicate entry.
            IF (IW(I).EQ.J) GO TO 720
            IW(I) = J
            IPOS = IP(I) - 1
            JCN(IPOS) = J
            IP(I) = IPOS
   90    CONTINUE
  100 CONTINUE
      DO 110 I = 1,M
         IW(I) = 0
  110 CONTINUE

C Check for zero rows and (unless ICNTL(4) > 0), compute chains of rows
C    with equal numbers of entries.
      IF (ICNTL(4).LE.0) THEN
         DO 120 I = 1,N
            IFIRST(I) = 0
  120    CONTINUE
         DO 130 I = M,1,-1
            NE1 = LENR(I)
            IF (NE1.GT.0) THEN
               IFIR = IFIRST(NE1)
               IFIRST(NE1) = I
               LASTR(I) = 0
               NEXTR(I) = IFIR
               IF (IFIR.GT.0) LASTR(IFIR) = I
            ELSE
               IP(I) = -M + NULLI
               NULLI = NULLI + 1
            END IF
  130    CONTINUE
      ELSE
         DO 140 I = M,1,-1
            NE1 = LENR(I)
            IF (NE1.EQ.0) THEN
               IP(I) = -M + NULLI
               NULLI = NULLI + 1
            END IF
  140    CONTINUE
      END IF
C Check for zero columns and compute chains of columns with equal
C   numbers of entries.
      DO 150 J = N,1,-1
         NE1 = LENC(J)
         IF (NE1.EQ.0) THEN
            IF (ICNTL(7).NE.2) THEN
               IF (J.LE.JLAST) THEN
                  NORD = NORD + 1
                  IQ(J) = -NORD
                  IF (NORD.EQ.JLAST) THEN
C We have ordered the first N - ICNTL(6) columns.
                     NORD = NORD + NULLJ
                     JLAST = N
                     NULLJ = 0
                  END IF
               ELSE
                  NULLJ = NULLJ + 1
                  IQ(J) = - (JLAST+NULLJ)
               END IF
               LASTC(J) = 0
               NEXTC(J) = 0
            END IF
         ELSE
            IFIR = JFIRST(NE1)
            JFIRST(NE1) = J
            NEXTC(J) = IFIR
            LASTC(J) = 0
            IF (IFIR.GT.0) LASTC(IFIR) = J
         END IF
  150 CONTINUE
      IF (INFO(6).EQ.0) THEN
         NULLC = NORD + NULLJ
         NULLR = NULLI
      END IF

C
C **********************************************
C ****    Start of main elimination loop    ****
C **********************************************
      DO 630 PIVOT = 1,N
C Check to see if reduced matrix should be considered as full.
         IF (NERED.GE. (MIN(CNTL(1),ONE)*(N-NORD))*
     +       (M-MORD)) GO TO 640

         IF (ICNTL(7).EQ.2) THEN
C Column order specified - choose the pivot within the column
            IPIV = 0
            J = IFIRST(PIVOT)
            IF (J.LT.1 .OR. J.GT.N) GO TO 730
            IF (IQ(J).LT.0) GO TO 730
            LEN = LENC(J)
            IF (LEN.LE.0) GO TO 320
            ALEN = LEN - 1
            I1 = IQ(J)
            I2 = I1 + LEN - 1
C Find largest entry in column
            II = IDAMAX(LEN,A(I1),1)
            MAXENT = ABS(A(I1+II-1))
C Is every entry in the column below the pivot threshold?
            IF (MAXENT.LE.CNTL(4)) GO TO 320
            AU = MAX(MAXENT*U,CNTL(4))
C Scan column for pivot
            DO 160 II = I1,I2
               IF (ABS(A(II)).LT.AU) GO TO 160
C Candidate satisfies threshold criterion.
               I = IRN(II)
               IF (IPIV.NE.0) THEN
                  IF (NEXTR(I).GE.NEXTR(IPIV)) GO TO 160
               END IF
               CPIV = ALEN*(LENR(I)-1)
               IJPOS = II
               IPIV = I
               JPIV = J
  160       CONTINUE
            GO TO 330
         END IF

C Find the least number of entries in a row or column (column only if
C   the Zlatev strategy is in use)
         LEN = MINC
         DO 170 MINC = LEN,M - MORD
            IF (JFIRST(MINC).NE.0) GO TO 180
            IF (ICNTL(4).LE.0) THEN
               IF (IFIRST(MINC).NE.0) GO TO 180
            END IF
  170    CONTINUE

C Find the next pivot or a column whose entries are all very small.
C CPIV is the Markowitz cost of the best pivot so far and PIVRAT is the
C      ratio of its absolute value to that of the largest entry in its
C      column.
  180    CPIV = M
         CPIV = CPIV*N
         PIVRAT = ZERO
C Examine columns/rows in order of ascending count.
         ISRCH = 0
         DO 300 LEN = MINC,M - MORD
            ALEN = LEN - 1
C Jump if Markowitz count cannot be bettered.
            IF (CPIV.LE.ALEN**2 .AND. ICNTL(4).LE.0) GO TO 310
            IJ = JFIRST(LEN)
C Scan columns with LEN entries.
            DO 220 IDUMMY = 1,N
C If no more columns with LEN entries, exit loop.
               IF (IJ.LE.0) GO TO 230
               J = IJ
               IJ = NEXTC(J)
               IF (J.GT.JLAST) GO TO 220
C Column J is now examined.
C First calculate multiplier threshold level.
               MAXENT = ZERO
               I1 = IQ(J)
               I2 = I1 + LEN - 1
               II = IDAMAX(LEN,A(I1),1)
               MAXENT = ABS(A(I1+II-1))
C Exit loop if every entry in the column is below the pivot threshold.
               IF (MAXENT.LE.CNTL(4)) GO TO 320
               AU = MAX(MAXENT*U,CNTL(4))
C If diagonal pivoting requested, look for diagonal entry.
               IF (ICNTL(7).EQ.1) THEN
                  DO 190 II = I1,I2
                     IF (IRN(II).EQ.J) GO TO 200
  190             CONTINUE
                  GO TO 220
  200             I1 = II
                  I2 = II
               END IF
C Scan column for possible pivots
               DO 210 II = I1,I2
                  IF (ABS(A(II)).LT.AU) GO TO 210
C Candidate satisfies threshold criterion.
                  I = IRN(II)
                  COST = ALEN*(LENR(I)-1)
                  IF (COST.GT.CPIV) GO TO 210
                  PIVR = ABS(A(II))/MAXENT
                  IF (COST.EQ.CPIV) THEN
                     IF (PIVR.LE.PIVRAT) GO TO 210
                  END IF
C Best pivot so far is found.
                  CPIV = COST
                  IJPOS = II
                  IPIV = I
                  JPIV = J
                  IF (CPIV.LE.ALEN**2 .AND. ICNTL(4).LE.0) GO TO 330
                  PIVRAT = PIVR
  210          CONTINUE
C Increment number of columns searched.
               ISRCH = ISRCH + 1
C Jump if we have searched the number of columns stipulated and found a
C   pivot.
               IF (ISRCH.GE.MSRCH) THEN
                  IF (PIVRAT.GT.ZERO) GO TO 330
               END IF
  220       CONTINUE
C
C Rows with LEN entries now examined.
  230       IF (ICNTL(4).GT.0) GO TO 300
            IF (CPIV.LE.ALEN*(ALEN+1)) GO TO 310
            IF (LEN.GT.N-NORD) GO TO 300
            IJ = IFIRST(LEN)
            DO 290 IDUMMY = 1,M
               IF (IJ.EQ.0) GO TO 300
               I = IJ
               IJ = NEXTR(IJ)
               J1 = IP(I)
               J2 = J1 + LEN - 1
C If diagonal pivoting requested, look for diagonal entry.
               IF (ICNTL(7).EQ.1) THEN
                  DO 240 JJ = J1,J2
                     IF (JCN(JJ).EQ.I) GO TO 250
  240             CONTINUE
                  GO TO 290
  250             J1 = JJ
                  J2 = JJ
               END IF
C Scan row I.
               DO 280 JJ = J1,J2
                  J = JCN(JJ)
                  IF (J.GT.JLAST) GO TO 280
                  COST = ALEN*(LENC(J)-1)
                  IF (COST.GE.CPIV) GO TO 280
C Pivot has best Markowitz count so far. Now check its suitability
C     on numerical grounds by examining other entries in its column.
                  I1 = IQ(J)
                  I2 = I1 + LENC(J) - 1
                  II = IDAMAX(LENC(J),A(I1),1)
                  MAXENT = ABS(A(I1+II-1))
                  DO 260 II = I1,I2 - 1
                     IF (IRN(II).EQ.I) GO TO 270
  260             CONTINUE
  270             JPOS = II
C Exit loop if every entry in the column is below the pivot threshold.
                  IF (MAXENT.LE.CNTL(4)) GO TO 320
                  IF (ABS(A(JPOS)).LT.MAXENT*U) GO TO 280
C Candidate satisfies threshold criterion.
                  CPIV = COST
                  IPIV = I
                  JPIV = J
                  IJPOS = JPOS
                  PIVRAT = ABS(A(JPOS))/MAXENT
                  IF (CPIV.LE.ALEN*(ALEN+1)) GO TO 330
  280          CONTINUE

  290       CONTINUE
C
  300    CONTINUE
  310    IF (PIVRAT.GT.ZERO) GO TO 330
C No pivot found. Switch to full matrix processing.
         INFO(1) = INFO(1) + 2
         IF (MP.GT.0) WRITE (MP,'(A/A)')
     +       ' Warning message from MA50AD: no suitable diagonal pivot',
     +       ' found, so switched to full matrix processing.'
         GO TO 640

C Every entry in the column is below the pivot threshold.
  320    IPIV = 0
         JPIV = J

C The pivot has now been found in position (IPIV,JPIV) in location
C     IJPOS in column file or all entries of column JPIV are very small
C     (IPIV=0).
C Update row and column ordering arrays to correspond with removal
C     of the active part of the matrix. Also update NEFACT.
  330    NEFACT = NEFACT + LENC(JPIV)
         PIVBEG = IQ(JPIV)
         PIVEND = PIVBEG + LENC(JPIV) - 1
         NORD = NORD + 1
         NORD1 = NORD
         IF (NORD.EQ.JLAST) THEN
C We have ordered the first N - ICNTL(6) columns.
            NORD = NORD + NULLJ
            JLAST = N
            NULLJ = 0
         END IF
         IF (ICNTL(4).LE.0) THEN
C Remove active rows from their row ordering chains.
            DO 340 II = PIVBEG,PIVEND
               I = IRN(II)
               LR = LASTR(I)
               NR = NEXTR(I)
               IF (NR.NE.0) LASTR(NR) = LR
               IF (LR.EQ.0) THEN
                  NE1 = LENR(I)
                  IFIRST(NE1) = NR
               ELSE
                  NEXTR(LR) = NR
               END IF
  340       CONTINUE
         END IF
         IF (IPIV.GT.0) THEN
C NEPR is number of entries in strictly U part of pivot row.
            NEPR = LENR(IPIV) - 1
            NEFACT = NEFACT + NEPR
            RINFO(1) = RINFO(1) + CPIV*2 + LENR(IPIV)
            J1 = IP(IPIV)
C Remove active columns from their column ordering chains.
            DO 350 JJ = J1,J1 + NEPR
               J = JCN(JJ)
               LC = LASTC(J)
               NC = NEXTC(J)
               IF (NC.NE.0) LASTC(NC) = LC
               IF (LC.EQ.0) THEN
                  NE1 = LENC(J)
                  JFIRST(NE1) = NC
               ELSE
                  NEXTC(LC) = NC
               END IF
  350       CONTINUE
C Move pivot to beginning of pivot column.
            IF (PIVBEG.NE.IJPOS) THEN
               ASW = A(PIVBEG)
               A(PIVBEG) = A(IJPOS)
               A(IJPOS) = ASW
               IRN(IJPOS) = IRN(PIVBEG)
               IRN(PIVBEG) = IPIV
            END IF
         ELSE
            NEPR = 0
            NE1 = LENC(JPIV)
            IF (CNTL(3).GT.ZERO) NDROP = NDROP + NE1
            IF (NE1.GT.0) THEN
C Remove column of small entries from its column ordering chain.
               LC = LASTC(JPIV)
               NC = NEXTC(JPIV)
               IF (NC.NE.0) LASTC(NC) = LC
               IF (LC.EQ.0) THEN
                  JFIRST(NE1) = NC
               ELSE
                  NEXTC(LC) = NC
               END IF
            END IF
         END IF
C
C Set up IW array so that IW(i) holds the relative position of row i
C    entry from beginning of pivot column.
         DO 360 II = PIVBEG + 1,PIVEND
            I = IRN(II)
            IW(I) = II - PIVBEG
  360    CONTINUE
C LENPIV is length of strictly L part of pivot column.
         LENPIV = PIVEND - PIVBEG
C
C Remove pivot column (including pivot) from row oriented file.
         DO 390 II = PIVBEG,PIVEND
            I = IRN(II)
            LENR(I) = LENR(I) - 1
            J1 = IP(I)
C J2 is last position in old row.
            J2 = J1 + LENR(I)
            DO 370 JJ = J1,J2 - 1
               IF (JCN(JJ).EQ.JPIV) GO TO 380
  370       CONTINUE
  380       JCN(JJ) = JCN(J2)
            JCN(J2) = 0
  390    CONTINUE

C For each active column, add the appropriate multiple of the pivot
C     column to it.
C We loop on the number of entries in the pivot row since the position
C     of this row may change because of compresses.
         DO 600 EYE = 1,NEPR
            J = JCN(IP(IPIV)+EYE-1)
C Search column J for entry to be eliminated, calculate multiplier,
C     and remove it from column file.
C  IDROP is the number of nonzero entries dropped from column J
C        because these fall beneath tolerance level.
            IDROP = 0
            JBEG = IQ(J)
            JEND = JBEG + LENC(J) - 1
            DO 400 II = JBEG,JEND - 1
               IF (IRN(II).EQ.IPIV) GO TO 410
  400       CONTINUE
  410       AMULT = -A(II)/A(IQ(JPIV))
            A(II) = A(JEND)
            IRN(II) = IRN(JEND)
            LENC(J) = LENC(J) - 1
            IRN(JEND) = 0
            JEND = JEND - 1
C Jump if pivot column is a singleton.
            IF (LENPIV.EQ.0) GO TO 600
C Now perform necessary operations on rest of non-pivot column J.
            IOP = 0
C Innermost loop.
CDIR$ IVDEP
            DO 420 II = JBEG,JEND
               I = IRN(II)
               IF (IW(I).GT.0) THEN
C Row i is involved in the pivot column.
                  IOP = IOP + 1
                  PIVCOL = IQ(JPIV) + IW(I)
C Flag IW(I) to show that the operation has been done.
                  IW(I) = -IW(I)
                  A(II) = A(II) + AMULT*A(PIVCOL)
               END IF
  420       CONTINUE

            IF (CNTL(3).GT.ZERO) THEN
C  Run through non-pivot column compressing column so that entries less
C      than CNTL(3) are not stored. All entries less than CNTL(3) are
C      also removed from the row structure.
               JNEW = JBEG
               DO 450 II = JBEG,JEND
                  IF (ABS(A(II)).GE.CNTL(3)) THEN
                     A(JNEW) = A(II)
                     IRN(JNEW) = IRN(II)
                     JNEW = JNEW + 1
                  ELSE
C  Remove non-zero entry from row structure.
                     I = IRN(II)
                     J1 = IP(I)
                     J2 = J1 + LENR(I) - 1
                     DO 430 JJ = J1,J2 - 1
                        IF (JCN(JJ).EQ.J) GO TO 440
  430                CONTINUE
  440                JCN(JJ) = JCN(J2)
                     JCN(J2) = 0
                     LENR(I) = LENR(I) - 1
                  END IF
  450          CONTINUE
               DO 460 II = JNEW,JEND
                  IRN(II) = 0
  460          CONTINUE
               IDROP = JEND + 1 - JNEW
               JEND = JNEW - 1
               LENC(J) = LENC(J) - IDROP
               NERED = NERED - IDROP
               INFO(6) = INFO(6) + IDROP
            END IF

C IFILL is fill-in left to do to non-pivot column J.
            IFILL = LENPIV - IOP
            NERED = NERED + IFILL
            INFO(3) = MAX(INFO(3),NERED+LENC(J))

C Treat no-fill case
            IF (IFILL.EQ.0) THEN
CDIR$ IVDEP
               DO 470 II = PIVBEG + 1,PIVEND
                  I = IRN(II)
                  IW(I) = -IW(I)
  470          CONTINUE
               GO TO 600
            END IF

C See if there is room for fill-in at end of the column.
            DO 480 IPOS = JEND + 1,MIN(JEND+IFILL,DISPC-1)
               IF (IRN(IPOS).NE.0) GO TO 490
  480       CONTINUE
            IF (IPOS.EQ.JEND+IFILL+1) GO TO 540
            IF (JEND+IFILL+1.LE.LA+1) THEN
               DISPC = JEND + IFILL + 1
               GO TO 540
            END IF
            IPOS = LA
            DISPC = LA + 1
C JMORE more spaces for fill-in are required.
  490       JMORE = JEND + IFILL - IPOS + 1
C We now look in front of the column to see if there is space for
C     the rest of the fill-in.
            DO 500 IPOS = JBEG - 1,MAX(JBEG-JMORE,1),-1
               IF (IRN(IPOS).NE.0) GO TO 510
  500       CONTINUE
            IPOS = IPOS + 1
            IF (IPOS.EQ.JBEG-JMORE) GO TO 520
C Column must be moved to the beginning of available storage.
  510       IF (DISPC+LENC(J)+IFILL.GT.LA+1) THEN
               INFO(2) = INFO(2) + 1
               CALL MA50DD(LA,A,IRN,IQ,N,DISPC,.TRUE.)
               JBEG = IQ(J)
               JEND = JBEG + LENC(J) - 1
               PIVBEG = IQ(JPIV)
               PIVEND = PIVBEG + LENC(JPIV) - 1
               IF (DISPC+LENC(J)+IFILL.GT.LA+1) GO TO 705
            END IF
            IPOS = DISPC
            DISPC = DISPC + LENC(J) + IFILL
C Move non-pivot column J.
  520       IQ(J) = IPOS
            DO 530 II = JBEG,JEND
               A(IPOS) = A(II)
               IRN(IPOS) = IRN(II)
               IPOS = IPOS + 1
               IRN(II) = 0
  530       CONTINUE
            JBEG = IQ(J)
            JEND = IPOS - 1
C Innermost fill-in loop which also resets IW.
C We know at this stage that there are IFILL positions free after JEND.
  540       IDROP = 0
            DO 580 II = PIVBEG + 1,PIVEND
               I = IRN(II)
               INFO(3) = MAX(INFO(3),NERED+LENR(I)+1)
               IF (IW(I).LT.0) THEN
                  IW(I) = -IW(I)
                  GO TO 580
               END IF
               ANEW = AMULT*A(II)
               IF (ABS(ANEW).LT.CNTL(3)) THEN
                  IDROP = IDROP + 1
               ELSE
                  JEND = JEND + 1
                  A(JEND) = ANEW
                  IRN(JEND) = I

C Put new entry in row file.
                  IEND = IP(I) + LENR(I)
                  IF (IEND.LT.DISPR) THEN
                     IF (JCN(IEND).EQ.0) GO TO 560
                  ELSE
                     IF (DISPR.LE.LA) THEN
                        DISPR = DISPR + 1
                        GO TO 560
                     END IF
                  END IF
                  IF (IP(I).GT.1) THEN
                     IF (JCN(IP(I)-1).EQ.0) THEN
C Copy row forward
                        IEND = IEND - 1
                        DO 545 JJ = IP(I),IEND
                           JCN(JJ-1) = JCN(JJ)
  545                   CONTINUE
                        IP(I) = IP(I) - 1
                        GO TO 560
                     END IF
                  END IF
                  IF (DISPR+LENR(I).GT.LA) THEN
C Compress.
                     INFO(2) = INFO(2) + 1
                     CALL MA50DD(LA,A,JCN,IP,M,DISPR,.FALSE.)
                     IF (DISPR+LENR(I).GT.LA) GO TO 705
                  END IF
C Copy row to first free position.
                  J1 = IP(I)
                  J2 = IP(I) + LENR(I) - 1
                  IP(I) = DISPR
                  DO 550 JJ = J1,J2
                     JCN(DISPR) = JCN(JJ)
                     JCN(JJ) = 0
                     DISPR = DISPR + 1
  550             CONTINUE
                  IEND = DISPR
                  DISPR = IEND + 1
  560             JCN(IEND) = J
                  LENR(I) = LENR(I) + 1
C End of adjustment to row file.
               END IF
  580       CONTINUE
            INFO(6) = INFO(6) + IDROP
            NERED = NERED - IDROP
            DO 590 II = 1,IDROP
               IRN(JEND+II) = 0
  590       CONTINUE
            LENC(J) = LENC(J) + IFILL - IDROP
C End of scan of pivot row.
  600    CONTINUE


C Remove pivot row from row oriented storage and update column
C     ordering arrays.  Remember that pivot row no longer includes
C     pivot.
         DO 610 EYE = 1,NEPR
            JJ = IP(IPIV) + EYE - 1
            J = JCN(JJ)
            JCN(JJ) = 0
            NE1 = LENC(J)
            LASTC(J) = 0
            IF (NE1.GT.0) THEN
               IFIR = JFIRST(NE1)
               JFIRST(NE1) = J
               NEXTC(J) = IFIR
               IF (IFIR.NE.0) LASTC(IFIR) = J
               MINC = MIN(MINC,NE1)
            ELSE IF (ICNTL(7).NE.2) THEN
               IF (INFO(6).EQ.0) NULLC = NULLC + 1
               IF (J.LE.JLAST) THEN
                  NORD = NORD + 1
                  IQ(J) = -NORD
                  IF (NORD.EQ.JLAST) THEN
C We have ordered the first N - ICNTL(6) columns.
                     NORD = NORD + NULLJ
                     JLAST = N
                     NULLJ = 0
                  END IF
               ELSE
                  NULLJ = NULLJ + 1
                  IQ(J) = - (JLAST+NULLJ)
               END IF
            END IF
  610    CONTINUE
         NERED = NERED - NEPR

C Restore IW and remove pivot column from column file.
C    Record the row permutation in IP(IPIV) and the column
C    permutation in IQ(JPIV), flagging them negative so that they
C    are not confused with real pointers in compress routine.
         IF (IPIV.NE.0) THEN
            LENR(IPIV) = 0
            IW(IPIV) = 0
            IRN(PIVBEG) = 0
            MORD = MORD + 1
            PIVBEG = PIVBEG + 1
            IP(IPIV) = -MORD
         END IF
         NERED = NERED - LENPIV - 1
         DO 620 II = PIVBEG,PIVEND
            I = IRN(II)
            IW(I) = 0
            IRN(II) = 0
            NE1 = LENR(I)
            IF (NE1.EQ.0) THEN
               IF (INFO(6).EQ.0) NULLR = NULLR + 1
               IP(I) = -M + NULLI
               NULLI = NULLI + 1
            ELSE IF (ICNTL(4).LE.0) THEN
C Adjust row ordering arrays.
               IFIR = IFIRST(NE1)
               LASTR(I) = 0
               NEXTR(I) = IFIR
               IFIRST(NE1) = I
               IF (IFIR.NE.0) LASTR(IFIR) = I
               MINC = MIN(MINC,NE1)
            END IF
  620    CONTINUE
         IQ(JPIV) = -NORD1
  630 CONTINUE
C We may drop through this loop with NULLI nonzero.

C ********************************************
C ****    End of main elimination loop    ****
C ********************************************

C Complete the permutation vectors
  640 INFO(5) = MORD + MIN(M-MORD-NULLI,N-NORD-NULLJ)
      DO 650 L = 1,MIN(M-MORD,N-NORD)
         RINFO(1) = RINFO(1) + M - MORD - L + 1 +
     +                   REAL(M-MORD-L)*(N-NORD-L)*2
  650 CONTINUE
      NP = NORD
      INFO(4) = 2 + NEFACT + M*2 + MAX(N-NORD+M-MORD,
     +          (N-NORD)*(M-MORD))
      INFO(6) = INFO(6) + NDROP
      INFO(7) = M - MORD
      DO 660 L = 1,M
         IF (IP(L).LT.0) THEN
            IP(L) = -IP(L)
         ELSE
            MORD = MORD + 1
            IP(L) = MORD
         END IF
  660 CONTINUE
      DO 670 L = 1,N
         IF (IQ(L).LT.0) THEN
            LASTC(L) = -IQ(L)
         ELSE
            IF (NORD.EQ.JLAST) NORD = NORD + NULLJ
            NORD = NORD + 1
            LASTC(L) = NORD
         END IF
  670 CONTINUE
C Store the inverse permutation
      DO 680 L = 1,N
         IQ(LASTC(L)) = L
  680 CONTINUE

C Test for rank deficiency
      IF (INFO(5).LT.MIN(M,N)) INFO(1) = INFO(1) + 1

      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(A,I6,A,F12.1/A,7I8)') ' Leaving MA50AD with NP =',
     +     NP,' RINFO(1) =',RINFO(1),' INFO =',(INFO(I),I=1,7)
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,'(A,(T6,10(I7)))') ' IP = ',IP
            WRITE (MP,'(A,(T6,10(I7)))') ' IQ = ',IQ
         END IF
      END IF

      GO TO 750

C Error conditions.
  690 INFO(1) = -1
      IF (LP.GT.0) WRITE (LP,'(/A/(2(A,I8)))')
     +    ' **** Error return from MA50AD ****',' M =',M,' N =',N
      GO TO 750
  700 INFO(1) = -2
      IF (LP.GT.0) WRITE (LP,'(/A/(A,I10))')
     +    ' **** Error return from MA50AD ****',' NE =',NE
      GO TO 750
  705 INFO(4) =  NEFACT + NERED
      INFO(6) = INFO(6) + NDROP
  710 INFO(1) = -3
      IF (LP.GT.0) WRITE (LP,'(/A/A,I9,A,I9)')
     +    ' **** Error return from MA50AD ****',
     +    ' LA  must be increased from',LA,' to at least',INFO(3)
      GO TO 750
  720 INFO(1) = -4
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))')
     +    ' **** Error return from MA50AD ****',' Entry in row',I,
     +    ' and column',J,' duplicated'
      GO TO 750
  730 INFO(1) = -5
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))')
     +    ' **** Error return from MA50AD ****',' Fault in component ',
     +    PIVOT,' of column permutation given in IFIRST'
      GO TO 750
  740 INFO(1) = -6
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))')
     +    ' **** Error return from MA50AD ****',' ICNTL(4) = ',ICNTL(4),
     +    ' when ICNTL(6) = 2'
  750 END


      SUBROUTINE MA50BD(M,N,NE,JOB,AA,IRNA,IPTRA,CNTL,ICNTL,IP,IQ,NP,
     +                  LFACT,FACT,IRNF,IPTRL,IPTRU,W,IW,INFO,RINFO)
C MA50B/BD factorizes the matrix in AA/IRNA/IPTRA as P L U Q where
C     P and Q are permutations, L is lower triangular, and U is unit
C     upper triangular. The prior information that it uses depends on
C     the value of the parameter JOB.
C
      INTEGER M,N,NE,JOB
      DOUBLE PRECISION AA(NE)
      INTEGER IRNA(NE),IPTRA(N)
      DOUBLE PRECISION CNTL(10)
      INTEGER ICNTL(20),IP(M),IQ(*),NP,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER IRNF(LFACT),IPTRL(N),IPTRU(N)
      DOUBLE PRECISION W(M)
      INTEGER IW(M+2*N),INFO(15)
      DOUBLE PRECISION RINFO(10)
C
C M is an integer variable that must be set to the number of rows.
C      It is not altered by the subroutine.
C N is an integer variable that must be set to the number of columns.
C      It is not altered by the subroutine.
C NE is an integer variable that must be set to the number of entries
C      in the input matrix.  It is not altered by the subroutine.
C JOB is an integer variable that must be set to the value 1, 2, or 3.
C     If JOB is equal to 1 and any of the first NP recommended pivots
C      fails to satisfy the threshold pivot tolerance, the row is
C      interchanged with the earliest row in the recommended sequence
C      that does satisfy the tolerance. Normal row interchanges are
C      performed in the last N-NP columns.
C     If JOB is equal to 2, then M, N, NE, IRNA, IPTRA, IP, IQ,
C      LFACT, NP, IRNF, IPTRL, and IPTRU must be unchanged since a
C      JOB=1 entry for the same matrix pattern and no interchanges are
C      performed among the first NP pivots; if ICNTL(6) > 0, the first
C      N-ICNTL(6) columns of AA must also be unchanged.
C     If JOB is equal to 3, ICNTL(6) must be in the range 1 to N-1.
C      The effect is as for JOB=2 except that interchanges are
C      performed.
C     JOB is not altered by the subroutine.
C AA is an array that holds the entries of the matrix and
C      is not altered.
C IRNA is an integer array of length NE that must be set to hold the
C      row indices of the corresponding entries in AA. It is not
C      altered.
C IPTRA is an integer array that holds the positions of the starts of
C      the columns of AA. It is not altered by the subroutine.
C CNTL  must be set by the user as follows and is not altered.
C     CNTL(2) determines the balance between pivoting for sparsity and
C       for stability, values near zero emphasizing sparsity and values
C       near one emphasizing stability.
C     CNTL(3) If this is set to a positive value, any entry whose
C       modulus is less than CNTL(3) will be dropped from the factors.
C       The factorization will then require less storage but will be
C       inaccurate.
C     CNTL(4)  Any entry of the reduced matrix whose modulus is less
C       than or equal to CNTL(4) will be regarded as zero from the
C        point of view of rank.
C ICNTL must be set by the user as follows and is not altered.
C     ICNTL(1)  must be set to the stream number for error messages.
C       A value less than 1 suppresses output.
C     ICNTL(2) must be set to the stream number for diagnostic output.
C       A value less than 1 suppresses output.
C     ICNTL(3) must be set to control the amount of output:
C       0 None.
C       1 Error messages only.
C       2 Error and warning messages.
C       3 As 2, plus scalar parameters and a few entries of array
C         parameters on entry and exit.
C       4 As 2, plus all parameters on entry and exit.
C     ICNTL(5) The block size to be used for full-matrix processing.
C       If <=0, the BLAS1 version is used.
C       If =1, the BLAS2 version is used.
C     ICNTL(6) If N > ICNTL(6) > 0, only the columns of A that
C       correspond to the last ICNTL(6) columns of the permuted matrix
C       may change prior to an entry with JOB > 1.
C     ICNTL(8) If this has value 1, there is no attempt to compute a
C       recommended value for LFACT if it is too small.
C IP is an integer array. If JOB=1, it must be set so that IP(I) < IP(J)
C      if row I is recommended to precede row J in the pivot sequence.
C      If JOB>1, it need not be set. If JOB=1 or JOB=3, IP(I) is set
C      to -K when row I is chosen for pivot K and IP is eventually
C      reset to recommend the chosen pivot sequence to a subsequent
C      JOB=1 entry. If JOB=2, IP is not be referenced.
C IQ is an integer array that must be set so that either IQ(J) is the
C      column in position J in the pivot sequence, J=1,2,...,N,
C      or IQ(1)=0 and the columns are taken in natural order.
C      It is not altered by the subroutine.
C NP is an integer variable that holds the number of columns to be
C      processed in packed storage. It is not altered by the subroutine.
C LFACT is an integer variable set to the size of FACT and IRNF.
C      It is not altered by the subroutine.
C FACT is an array that need not be set on a JOB=1 entry and must be
C      unchanged since the previous entry if JOB>1. On return, FACT(1)
C      holds the value of CNTL(3) used, FACT(2) will holds the value
C      of CNTL(4) used, FACT(3:IPTRL(N)) holds the packed part of L/U
C      by columns, and the full part of L/U is held by columns
C      immediately afterwards. U has unit diagonal entries, which are
C      not stored. In each column of the packed part, the entries of
C      U precede the entries of L; also the diagonal entries of L
C      head each column of L and are reciprocated.
C IRNF is an integer array of length LFACT that need not be set on
C      a JOB=1 entry and must be unchanged since the previous entry
C      if JOB>1. On exit, IRNF(1) holds the number of dropped entries,
C      IRNF(2) holds the number of rows MF in full storage,
C      IRNF(3:IPTRL(N)) holds the row numbers of the packed part
C      of L/U, IRNF(IPTRL(N)+1:IPTRL(N)+MF) holds the row indices
C      of the full part of L/U, and IRNF(IPTRL(N)+MF+I), I=1,2,..,N-NP
C      holds the vector IPIV output by MA50GD.
C      If JOB=2, IRNF will be unaltered.
C IPTRL is an integer array that need not be set on a JOB=1 entry and
C     must be unchanged since the previous entry if JOB>1.
C     For J = 1,..., NP, IPTRL(J) holds the position in
C     FACT and IRNF of the end of column J of L.
C     For J = NP+1,..., N, IPTRL(J) is equal to IPTRU(J).
C IPTRU is an integer array that need not be set on a JOB=1 entry and
C     must be unchanged since the previous entry if JOB>1.
C     For J = 1,..., N, IPTRU(J) holds the position in
C     FACT and IRNF of the end of the packed part of column J of U.
C W is an array of length M used as workspace for holding
C      the expanded form of a sparse vector.
C IW is an integer array of length M+2*N used as workspace.
C INFO need not be set on entry. On exit, it holds the following:
C    INFO(1) A negative value will indicate an error return and a
C       positive value a warning. Possible nonzero values are:
C      -1  M < 1 or N < 1.
C      -2  NE < 0.
C      -3  Insufficient space.
C      -4  There are duplicated entries.
C      -5  JOB < 1, 3 when ICNTL(6)=0, or > 3.
C      -6  JOB = 2, but entries were dropped in the corresponding JOB=1
C          entry.
C      -7  NP < 0 or NP > N.
C     -(7+K) Pivot too small in column K when JOB=2.
C      +1  Rank deficient.
C    INFO(4) Minimum storage required to factorize matrix
C            if INFO(1) >= 0. Recommended value for LFACT
C            if ICNTL(8) = 0 and INFO(1) = -3.
C    INFO(5) Computed rank of the matrix.
C    INFO(6) Number of entries dropped from the data structure.
C    INFO(7) Number of rows processed in full storage.
C RINFO need not be set on entry. On exit, RINFO(1) holds the number of
C    floating-point operations performed.

      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0D0,ONE=1.0D0)
      DOUBLE PRECISION AMULT,ASW
      INTEGER BEGCOL
      LOGICAL DROP
      INTEGER ENDCOL,EYE,EYE1,I,IA1,IA2,IF1,IF2,II,IL1,IL2,IPIV,IQPIV,
     +        IU1,IU2,ISW,J,JDUMMY,JJ,JLAST,K,LP
      DOUBLE PRECISION MAXENT
      INTEGER MF,MORD,MP,NEU,NF,NULLC
      DOUBLE PRECISION PIVLIM
      INTEGER RANK
      DOUBLE PRECISION U
C AMULT Temporary variable used to store current multiplier.
C ASW Temporary variable used when swopping two real quantities.
C BEGCOL is pointer to beginning of section of column when pruning.
C DROP True if any entries dropped from current column.
C ENDCOL is pointer to end of section of column when pruning.
C EYE Running position for current column.
C EYE1 Position of the start of second current column.
C I Temporary variable holding row number. Also used as index in DO
C     loops used in initialization of arrays.
C IA1 Position of the start of the current column in AA.
C IA2 Position of the end of the current column in AA.
C IF1 Position of the start of the full submatrix.
C IF2 Position of the end of the full submatrix.
C II Running position for current column.
C IL1 Position of the first entry of the current column of L.
C IL2 Position of the last entry of the current column of L.
C IPIV Position of the pivot in FACT and IRNF.
C IQPIV Recommended position of the pivot row in the pivot sequence.
C IU1 Position of the start of current column of U.
C IU2 Position of the end of the current column of U.
C ISW Temporary variable used when swopping two integer quantities.
C J Temporary variable holding column number.
C JDUMMY DO index not referenced in the loop.
C JJ Running position for current column.
C JLAST The lesser of NP and the last column of A for which no new
C     factorization operations are needed.
C K Temporary variable holding the current pivot step in the elimination
C LP Unit for error messages.
C MAXENT Temporary variable used to hold value of largest entry in
C    column.
C MF Number of rows in full block.
C MORD Number of rows ordered.
C MP Unit for diagnostic messages.
C NEU Number of entries omitted from U and the full block in order to
C    calculate INFO(4) (0 unless INFO(1)=-3).
C NF Number of columns in full block.
C NULLC Number of columns found null before dropping any elements.
C PIVLIM Limit on pivot size.
C RANK Value returned by MA50E/ED or MA50F/FD
C U Used to hold local copy of CNTL(2), changed if necessary so that it
C    is in range.
C
      EXTERNAL MA50ED,MA50FD,MA50GD
      INTRINSIC ABS,MAX,MIN
C LAPACK subroutine for triangular factorization.

      INFO(1) = 0
      INFO(4) = 0
      INFO(5) = 0
      INFO(6) = 0
      INFO(7) = 0
      RINFO(1) = ZERO
      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (ICNTL(3).LE.0) LP = 0
      IF (ICNTL(3).LE.1) MP = 0

C Check input values
      IF (M.LT.1 .OR. N.LT.1) THEN
         INFO(1) = -1
         IF (LP.GT.0) WRITE (LP,'(/A/A,I8,A,I8)')
     +       ' **** Error return from MA50BD ****',' M =',M,' N =',N
         GO TO 550
      END IF
      IF (NE.LE.0) THEN
         INFO(1) = -2
         IF (LP.GT.0) WRITE (LP,'(/A/A,I6)')
     +       ' **** Error return from MA50BD ****',' NE =',NE
         GO TO 550
      END IF
      IF (NP.LT.0 .OR. NP.GT.N) THEN
         INFO(1) = -7
         IF (LP.GT.0) WRITE (LP,'(/A/A,I8,A,I8)')
     +       ' **** Error return from MA50BD ****',' NP =',NP,' N =',N
         GO TO 550
      END IF
      IF (LFACT.LT.MAX(M,NE+2)) THEN
         INFO(4) = MAX(M,NE+2)
         GO TO 520
      END IF
      IF (JOB.EQ.1) THEN
      ELSE IF (JOB.EQ.2 .OR. JOB.EQ.3) THEN
         IF (IRNF(1).NE.0) THEN
            INFO(1) = -6
            IF (LP.GT.0) WRITE (LP,'(/A/A,I1,A)')
     +          ' **** Error return from MA50BD ***',' Call with JOB=',
     +          JOB,' follows JOB=1 call in which entries were dropped'
            GO TO 550
         END IF
      ELSE
         INFO(1) = -5
         IF (LP.GT.0) WRITE (LP,'(/A/A,I2)')
     +       ' **** Error return from MA50BD ****',' JOB =',JOB
         GO TO 550
      END IF

C Print input data
      IF (MP.GT.0) THEN
         IF (ICNTL(3).GT.2) WRITE (MP,
     +       '(/2(A,I6),A,I8,A,I3/A,I8,A,I7/A,1P,4E10.2/A,7I8)')
     +       ' Entering MA50BD with M =',M,' N =',N,' NE =',NE,' JOB =',
     +       JOB,' LFACT =',LFACT,' NP =',NP,' CNTL =',(CNTL(I),I=1,4),
     +       ' ICNTL =',(ICNTL(I),I=1,7)
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,'(A,(T6,10(I7)))') ' IP = ',IP
            IF (IQ(1).GT.0) THEN
               WRITE (MP,'(A,(T6,10(I7)))') ' IQ = ',(IQ(J),J=1,N)
            ELSE
               WRITE (MP,'(A,(T6,I7))') ' IQ = ',IQ(1)
            END IF
            DO 10 J = 1,N - 1
               IF (IPTRA(J).LT.IPTRA(J+1)) WRITE (MP,
     +             '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',J,
     +             (AA(II),IRNA(II),II=IPTRA(J),IPTRA(J+1)-1)
   10       CONTINUE
            IF (IPTRA(N).LE.NE) WRITE (MP,
     +          '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',N,
     +          (AA(II),IRNA(II),II=IPTRA(N),NE)
         END IF
      END IF

C Initializations.
      JLAST = 0
      NULLC = 0
      IF (JOB.GT.1 .AND. ICNTL(6).GT.0 .AND.
     +    ICNTL(6).LT.N) JLAST = MIN(NP,N-ICNTL(6))

      U = MIN(CNTL(2),ONE)
      U = MAX(U,ZERO)
      DO 20 I = 1,M
         IW(I+N) = 0
         W(I) = ZERO
   20 CONTINUE
      MORD = 0
      IF1 = LFACT + 1
      IF2 = 0
      NF = N - NP
      MF = 0
      IL2 = 2
      IF (JLAST.GT.0) IL2 = IPTRL(JLAST)
      NEU = 0

C Jump if JOB is equal to 2.
      IF (JOB.EQ.2) GO TO 370

      IF (JOB.EQ.3) THEN
C Reconstruct IP and set MORD
         DO 30 J = 1,NP
            IA1 = IPTRU(J) + 1
            IF (IA1.GT.IPTRL(J)) GO TO 30
            IF (J.LE.JLAST) THEN
               MORD = MORD + 1
               IP(IRNF(IA1)) = -J
            ELSE
               IP(IRNF(IA1)) = J
            END IF
   30    CONTINUE
         MF = IRNF(2)
         IA1 = IPTRL(N)
         DO 40 J = 1,MF
            IP(IRNF(IA1+J)) = NP + J
   40    CONTINUE
      END IF

C Store copies of column ends ready for pruning
      DO 50 K = 1,JLAST
         IW(M+N+K) = IPTRL(K)
   50 CONTINUE

C Each pass through this main loop processes column K.
      DO 310 K = JLAST + 1,N
         DROP = .FALSE.
         IF (K.EQ.NP+1) THEN
C Set up data structure for full part.
            MF = M - MORD
            IF1 = LFACT + 1 - MF
            II = 0
            DO 60 I = 1,M
               IF (IP(I).GT.0) THEN
                  IW(I+N) = N
                  IRNF(IF1+II) = I
                  II = II + 1
                  IP(I) = NP + II
               END IF
   60       CONTINUE
            IF1 = LFACT + 1 - MAX(MF*NF,MF+NF)
            IF2 = IF1 - 1 + MF*MAX(0,JLAST-NP)
         END IF
         J = K
         IF (IQ(1).GT.0) J = IQ(K)
         IA1 = IPTRA(J)
         IA2 = NE
         IF (J.NE.N) IA2 = IPTRA(J+1) - 1
         IU1 = IL2 + 1
         IU2 = IU1 - 1
         IL1 = IF1 - 1 + IA1 - IA2
         IL2 = IL1 - 1
         INFO(4) = MAX(INFO(4),NEU+LFACT-IL1+IU2+M+1)
         IF (IL1-IU2.LE.M) THEN
            IF (INFO(1).NE.-3) THEN
               INFO(1) = -3
               IF (ICNTL(8).NE.0) GO TO 480
C Get rid of U info.
               NEU = IL2 + LFACT + 1 - MF - IF1
               IF1 = LFACT + 1 - MF
               IF2 = IF1 - 1
               IL2 = 0
               EYE = 0
               DO 80 J = 1,MIN(K-1,NP)
                  IU2 = IPTRU(J)
                  IPTRU(J) = EYE
                  IL2 = IPTRL(J)
                  NEU = NEU + IU2 - IL2
                  DO 70 II = IU2 + 1,IL2
                     EYE = EYE + 1
                     IRNF(EYE) = IRNF(II)
                     FACT(EYE) = FACT(II)
   70             CONTINUE
                  IPTRL(J) = EYE
                  IW(M+N+J) = EYE
   80          CONTINUE
               IU1 = EYE + 1
               IU2 = EYE
               IL1 = IF1 - 1 + IA1 - IA2
               IL2 = IL1 - 1
            END IF
C Quit if LFACT is much too small
            IF (IL1-IU2.LE.M) GO TO 480
         END IF
C Load column K of AA into full vector W and into the back of IRNF.
C Check for duplicates.
         EYE = IL1
         DO 90 II = IA1,IA2
            I = IRNA(II)
            IF (IW(I+N).EQ.-1) GO TO 540
            IW(I+N) = -1
            W(I) = AA(II)
            IRNF(EYE) = I
            EYE = EYE + 1
   90    CONTINUE
C Depth first search to find topological order for triangular solve
C     and structure of column K of L/U
C IW(J) is used to hold a pointer to next entry in column J
C     during the depth-first search at stage K, J = 1,..., N.
C IW(I+N) is set to K when row I has been processed, and to N for rows
C     of the full part once column NP has been passed. It is also
C     used for backtracking, a negative value being used to point to the
C     previous row in the chain.
C IW(M+N+I) is set to the position in FACT and IRNF of the end of the
C     active part of the column after pruning.  It is initially set to
C     IPTRL(I) and is flagged negative when column has been pruned.
C Set IPTRL temporarily for column K so that special code is
C     not required to process this column.
         IPTRL(K) = EYE - 1
         IW(M+N+K) = EYE - 1
C IW(K) is set to beginning of original column K.
         IW(K) = IL1
         J = K
C The outer loop of the depth-first search is executed once for column
C      K and twice for each entry in the upper-triangular part of column
C      K (once to initiate a search in the corresponding column and
C      once when the search in the column is finished).
         DO 120 JDUMMY = 1,2*K
C Look through column J of L (or column K of A). All the entries
C     are entries of the filled-in column K. Store new entries of the
C     lower triangle and continue until reaching an entry of the upper
C     triangle.
            DO 100 II = IW(J),ABS(IW(M+N+J))
               I = IRNF(II)
C Jump if index I already encountered in column K or is in full part.
               IF (IW(I+N).GE.K) GO TO 100
               IF (IP(I).LE.0) GO TO 110
C Entry is in lower triangle. Flag it and store it in L.
               IW(I+N) = K
               IL1 = IL1 - 1
               IRNF(IL1) = I
  100       CONTINUE
            IF (J.EQ.K) GO TO 130
C Flag J, put its row index into U, and backtrack
            IU2 = IU2 + 1
            I = IRNF(IPTRU(J)+1)
            IRNF(IU2) = I
            J = -IW(I+N)
            IW(I+N) = K
            GO TO 120
C Entry in upper triangle.  Move search to corresponding column.
  110       IW(I+N) = -J
            IW(J) = II + 1
            J = -IP(I)
            IW(J) = IPTRU(J) + 2
  120    CONTINUE
C Run through column K of U in the lexicographical order that was just
C     constructed, performing elimination operations.
  130    DO 150 II = IU2,IU1,-1
            I = IRNF(II)
            J = -IP(I)
C Add multiple of column J of L to column K
            EYE1 = IPTRU(J) + 1
            IF (ABS(W(I)).LT.CNTL(3)) GO TO 150
            AMULT = -W(I)*FACT(EYE1)
C Note we are storing negative multipliers
            W(I) = AMULT
            DO 140 EYE = EYE1 + 1,IPTRL(J)
               I = IRNF(EYE)
               W(I) = W(I) + AMULT*FACT(EYE)
  140       CONTINUE
            RINFO(1) = RINFO(1) + ONE + 2*(IPTRL(J)-EYE1)
  150    CONTINUE

C Unload reals of column of U and set pointer
         IF (CNTL(3).GT.ZERO) THEN
            EYE = IU1
            DO 160 II = IU1,IU2
               I = IRNF(II)
               IF (ABS(W(I)).LT.CNTL(3)) THEN
                  INFO(6) = INFO(6) + 1
               ELSE
                  IRNF(EYE) = -IP(I)
                  FACT(EYE) = W(I)
                  EYE = EYE + 1
               END IF
               W(I) = ZERO
  160       CONTINUE
            IU2 = EYE - 1
         ELSE
            DO 170 II = IU1,IU2
               I = IRNF(II)
               IRNF(II) = -IP(I)
               FACT(II) = W(I)
               W(I) = ZERO
  170       CONTINUE
         END IF
         IF (INFO(1).EQ.-3) THEN
            NEU = NEU + IU2 - IU1 + 1
            IU2 = IU1 - 1
         END IF
         IPTRU(K) = IU2
         IF (K.LE.NP) THEN
C Find the largest entry in the column and drop any small entries
            MAXENT = ZERO
            IF (CNTL(3).GT.ZERO) THEN
               EYE = IL1
               DO 180 II = IL1,IL2
                  I = IRNF(II)
                  IF (ABS(W(I)).LT.CNTL(3)) THEN
                     INFO(6) = INFO(6) + 1
                     W(I) = ZERO
                     DROP = .TRUE.
                  ELSE
                     IRNF(EYE) = I
                     EYE = EYE + 1
                     MAXENT = MAX(ABS(W(I)),MAXENT)
                  END IF
  180          CONTINUE
               IL2 = EYE - 1
            ELSE
               DO 190 II = IL1,IL2
                  MAXENT = MAX(ABS(W(IRNF(II))),MAXENT)
  190          CONTINUE
            END IF
C Unload column of L, performing pivoting and moving indexing
C      information.
            PIVLIM = U*MAXENT
            EYE = IU2
            IQPIV = M + N
            IF (IL1.GT.IL2) NULLC = NULLC + 1
            DO 200 II = IL1,IL2
               I = IRNF(II)
               EYE = EYE + 1
               IRNF(EYE) = I
               FACT(EYE) = W(I)
               W(I) = ZERO
C Find position of pivot
               IF (ABS(FACT(EYE)).GE.PIVLIM) THEN
                  IF (ABS(FACT(EYE)).GT.CNTL(4)) THEN
                     IF (IP(I).LT.IQPIV) THEN
                        IQPIV = IP(I)
                        IPIV = EYE
                     END IF
                  END IF
               END IF
  200       CONTINUE
            IL1 = IU2 + 1
            IL2 = EYE
            IF (IL1.LE.IL2) THEN
C Column is not null
               IF (IQPIV.EQ.M+N) THEN
C All entries in the column are too small to be pivotal. Drop them all.
                  IF (CNTL(3).GT.ZERO) INFO(6) = INFO(6) + EYE - IU2
                  IL2 = IU2
               ELSE
                  IF (IL1.NE.IPIV) THEN
C Move pivot to front of L
                     ASW = FACT(IPIV)
                     FACT(IPIV) = FACT(IL1)
                     FACT(IL1) = ASW
                     ISW = IRNF(IL1)
                     IRNF(IL1) = IRNF(IPIV)
                     IRNF(IPIV) = ISW
                  END IF
C Reciprocate pivot
                  INFO(5) = INFO(5) + 1
                  FACT(IL1) = ONE/FACT(IL1)
                  RINFO(1) = RINFO(1) + ONE
C Record pivot row
                  MORD = MORD + 1
                  IP(IRNF(IL1)) = -K
               END IF
            END IF
         ELSE
C Treat column as full
            IL2 = IPTRU(K)
CDIR$ IVDEP
            DO 210 II = LFACT - MF + 1,LFACT
               I = IRNF(II)
               IF2 = IF2 + 1
               FACT(IF2) = W(I)
               W(I) = ZERO
  210       CONTINUE
            IF (INFO(1).EQ.-3) IF2 = IF2 - MF
         END IF
         IW(M+N+K) = IL2
         IPTRL(K) = IL2
         IF (DROP) GO TO 310
C Scan columns involved in update of column K and remove trailing block.
         DO 300 II = IU1,IU2
            I = IRNF(II)
C Jump if column already pruned.
            IF (IW(M+N+I).LT.0) GO TO 300
            BEGCOL = IPTRU(I) + 2
            ENDCOL = IPTRL(I)
C Scan column to see if there is an entry in the current pivot row.
            IF (K.LE.NP) THEN
               DO 220 JJ = BEGCOL,ENDCOL
                  IF (IP(IRNF(JJ)).EQ.-K) GO TO 230
  220          CONTINUE
               GO TO 300
            END IF
C Sort the entries so that those in rows already pivoted (negative IP
C    values) precede the rest.
  230       DO 280 JDUMMY = BEGCOL,ENDCOL
               JJ = BEGCOL
               DO 240 BEGCOL = JJ,ENDCOL
                  IF (IP(IRNF(BEGCOL)).GT.0) GO TO 250
  240          CONTINUE
               GO TO 290
  250          JJ = ENDCOL
               DO 260 ENDCOL = JJ,BEGCOL,-1
                  IF (IP(IRNF(ENDCOL)).LT.0) GO TO 270
  260          CONTINUE
               GO TO 290
  270          ASW = FACT(BEGCOL)
               FACT(BEGCOL) = FACT(ENDCOL)
               FACT(ENDCOL) = ASW
               J = IRNF(BEGCOL)
               IRNF(BEGCOL) = IRNF(ENDCOL)
               IRNF(ENDCOL) = J
               BEGCOL = BEGCOL + 1
               ENDCOL = ENDCOL - 1
  280       CONTINUE
  290       IW(M+N+I) = -ENDCOL
  300    CONTINUE
  310 CONTINUE
      IF (N.EQ.NP) THEN
C Set up data structure for the (null) full part.
         MF = M - MORD
         IF1 = LFACT + 1 - MF
         II = 0
         DO 320 I = 1,M
            IF (IP(I).GT.0) THEN
               IW(I+N) = N
               IRNF(IF1+II) = I
               II = II + 1
               IP(I) = NP + II
            END IF
  320    CONTINUE
         IF1 = LFACT + 1 - MAX(MF*NF,MF+NF)
         IF2 = IF1 - 1 + MF*MAX(0,JLAST-NP)
      END IF
      IF (INFO(5).EQ.MIN(M,N)) THEN
C Restore sign of IP
         DO 330 I = 1,M
            IP(I) = ABS(IP(I))
  330    CONTINUE
      ELSE
C Complete IP
         MORD = NP
         DO 340 I = 1,M
            IF (IP(I).LT.0) THEN
               IP(I) = -IP(I)
            ELSE
               MORD = MORD + 1
               IP(I) = MORD
            END IF
  340    CONTINUE
      END IF
      IRNF(1) = INFO(6)
      IRNF(2) = MF
      INFO(7) = MF
      FACT(1) = CNTL(3)
      FACT(2) = CNTL(4)
      IF (INFO(1).EQ.-3) GO TO 520
C Move full part forward
      IF2 = IF2 - MF*NF
      DO 350 II = 1,MF*NF
         FACT(IL2+II) = FACT(IF1-1+II)
  350 CONTINUE
      DO 360 II = 1,MF
         IRNF(IL2+II) = IRNF(LFACT-MF+II)
  360 CONTINUE
      IF1 = IL2 + 1
      GO TO 440
C
C Fast factor (JOB = 2)
C Each pass through this main loop processes column K.
  370 MF = IRNF(2)
      IF1 = IPTRL(N) + 1
      IF2 = IF1 - 1
      DO 430 K = JLAST + 1,N
         J = K
         IF (IQ(1).GT.0) J = IQ(K)
         IA1 = IPTRA(J)
         IA2 = NE
         IF (J.NE.N) IA2 = IPTRA(J+1) - 1
         IU1 = IL2 + 1
         IU2 = IPTRU(K)
         IL1 = IU2 + 1
         IL2 = IPTRL(K)
C Load column K of A into full vector W
         DO 380 II = IA1,IA2
            W(IRNA(II)) = AA(II)
  380    CONTINUE
C Run through column K of U in lexicographical order, performing
C      elimination operations.
         DO 400 II = IU2,IU1,-1
            J = IRNF(II)
            I = IRNF(IPTRU(J)+1)
C Add multiple of column J of L to column K
            EYE1 = IPTRU(J) + 1
            AMULT = -W(I)*FACT(EYE1)
C Note we are storing negative multipliers
            FACT(II) = AMULT
            W(I) = ZERO
            DO 390 EYE = EYE1 + 1,IPTRL(J)
               I = IRNF(EYE)
               W(I) = W(I) + AMULT*FACT(EYE)
  390       CONTINUE
            RINFO(1) = RINFO(1) + ONE + 2*(IPTRL(J)-EYE1)
  400    CONTINUE
         IF (K.LE.NP) THEN
            IF (IL1.LE.IL2) THEN
C Load column of L.
CDIR$ IVDEP
               DO 410 II = IL1,IL2
                  I = IRNF(II)
                  FACT(II) = W(I)
                  W(I) = ZERO
  410          CONTINUE
C Test pivot. Note that this is the only numerical test when JOB = 2.
               IF (ABS(FACT(IL1)).LE.CNTL(4)) THEN
                  GO TO 530
               ELSE
C Reciprocate pivot
                  INFO(5) = INFO(5) + 1
                  FACT(IL1) = ONE/FACT(IL1)
                  RINFO(1) = RINFO(1) + ONE
               END IF
            END IF
         ELSE
C Treat column as full
            DO 420 II = IF1,IF1 + MF - 1
               I = IRNF(II)
               IF2 = IF2 + 1
               FACT(IF2) = W(I)
               W(I) = ZERO
  420       CONTINUE
         END IF
  430 CONTINUE
      INFO(4) = MAX(IF1+MF+NF-1,IF2)

  440 IF (MF.GT.0 .AND. NF.GT.0) THEN
C Factorize full block
         IF (ICNTL(5).GT.1) CALL MA50GD(MF,NF,FACT(IF1),MF,ICNTL(5),
     +                                  CNTL(4),IRNF(IF1+MF),RANK)
         IF (ICNTL(5).EQ.1) CALL MA50FD(MF,NF,FACT(IF1),MF,CNTL(4),
     +                                  IRNF(IF1+MF),RANK)
         IF (ICNTL(5).LE.0) CALL MA50ED(MF,NF,FACT(IF1),MF,CNTL(4),
     +                                  IRNF(IF1+MF),RANK)
         INFO(5) = INFO(5) + RANK
         DO 450 I = 1,MIN(MF,NF)
            RINFO(1) = RINFO(1) + MF - I + 1 + REAL(MF-I)*(NF-I)*2
  450    CONTINUE
      END IF
      IF (INFO(5).LT.MIN(M,N)) INFO(1) = 1
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(A,I6,A,F12.1/A,I3,A,4I8)')
     +     ' Leaving MA50BD with IRNF(2) =',IRNF(2),
     +     ' RINFO(1) =',RINFO(1),
     +     ' INFO(1) =',INFO(1),' INFO(4:7) =', (INFO(J),J=4,7)
         IF (ICNTL(3).GT.3) THEN
            IF (JOB.NE.2) WRITE (MP,'(A,(T6,10(I7)))') ' IP = ',IP
            DO 460 J = 1,N
               IF (J.GT.1) THEN
                  IF (IPTRL(J-1).LT.IPTRU(J)) WRITE (MP,
     +                '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,
     +                ' of U', (FACT(II),IRNF(II),II=IPTRL(J-1)+1,
     +                IPTRU(J))
               END IF
               IF (IPTRU(J).LT.IPTRL(J)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of L',
     +              (FACT(II),IRNF(II),II=IPTRU(J)+1,IPTRL(J))
  460       CONTINUE
            WRITE (MP,'(A)') ' Full part'
            WRITE (MP,'((6I12))') (IRNF(IF1+MF+J),J=0,NF-1)
            DO 470 I = 0,MF - 1
               WRITE (MP,'(I4,1P,6E12.4:/(4X,1P,6E12.4))') IRNF(IF1+I),
     +            (FACT(IF1+I+J*MF),J=0,NF-1)
  470       CONTINUE
         END IF
      END IF
      GO TO 550

C Error conditions
C LFACT is much too small or ICNTL(8)/=0. Patch up IP and quit.
  480 DO 490 I = 1,M
         IW(I) = 0
  490 CONTINUE
      DO 500 I = 1,M
         IF (IP(I).GT.0) THEN
            IW(IP(I)) = I
         ELSE
            IP(I) = -IP(I)
         END IF
  500 CONTINUE
      DO 510 I = 1,M
         IF (IW(I).GT.0) THEN
            IP(IW(I)) = K
            K = K + 1
         END IF
  510 CONTINUE
  520 INFO(1) = -3
      IF (LP.GT.0) THEN
         WRITE (LP,'(/A)')' **** Error return from MA50BD **** '
         IF (ICNTL(8).EQ.0) THEN
            WRITE (LP,'(A,I7,A,I7)')' LFACT must be increased from',
     +         LFACT,' to at least',INFO(4)
         ELSE
            WRITE (LP,'(A,I7)')' LFACT must be increased from',LFACT
         END IF
      END IF
      GO TO 550
  530 INFO(1) = - (7+K)
      IF (LP.GT.0) WRITE (LP,'(/A/A,I6,A)')
     +    ' **** Error return from MA50BD **** ',
     +    ' Small pivot found in column',K,' of the permuted matrix.'
      GO TO 550
  540 INFO(1) = -4
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))')
     +    ' **** Error return from MA50BD ****',' Entry in row',I,
     +    ' and column',J,' duplicated'
  550 END

      SUBROUTINE MA50CD(M,N,ICNTL,IQ,NP,TRANS,LFACT,FACT,IRNF,IPTRL,
     +                  IPTRU,B,X,W,INFO)
C MA50C/CD uses the factorization produced by
C     MA50B/BD to solve A x = b or (A trans) x = b.
C
      INTEGER M,N,ICNTL(20),IQ(*),NP
      LOGICAL TRANS
      INTEGER LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER IRNF(LFACT),IPTRL(N),IPTRU(N)
      DOUBLE PRECISION B(*),X(*),W(*)
      INTEGER INFO(15)
C
C M  is an integer variable set to the number of rows.
C     It is not altered by the subroutine.
C N  is an integer variable set to the number of columns.
C     It is not altered by the subroutine.
C ICNTL must be set by the user as follows and is not altered.
C     ICNTL(1)  must be set to the stream number for error messages.
C       A value less than 1 suppresses output.
C     ICNTL(2) must be set to the stream number for diagnostic output.
C       A value less than 1 suppresses output.
C     ICNTL(3) must be set to control the amount of output:
C       0 None.
C       1 Error messages only.
C       2 Error and warning messages.
C       3 As 2, plus scalar parameters and a few entries of array
C         parameters on entry and exit.
C       4 As 2, plus all parameters on entry and exit.
C     ICNTL(5) must be set to control the level of BLAS used:
C       0 Level 1 BLAS.
C      >0 Level 2 BLAS.
C IQ is an integer array holding the permutation Q.
C     It is not altered by the subroutine.
C NP is an integer variable that must be unchanged since calling
C     MA50B/BD. It holds the number of rows and columns in packed
C     storage. It is not altered by the subroutine.
C TRANS a logical variable thatmust be set to .TRUE. if (A trans)x = b
C     is to be solved and to .FALSE. if A x = b is to be solved.
C     TRANS is not altered by the subroutine.
C LFACT is an integer variable set to the size of FACT and IRNF.
C     It is not altered by the subroutine.
C FACT is an array that must be unchanged since calling MA50B/BD. It
C     holds the packed part of L/U by columns, and the full part of L/U
C     by columns. U has unit diagonal entries, which are not stored, and
C     the signs of the off-diagonal entries are inverted.  In the packed
C     part, the entries of U precede the entries of L; also the diagonal
C     entries of L head each column of L and are reciprocated.
C     FACT is not altered by the subroutine.
C IRNF is an integer array that must be unchanged since calling
C     MA50B/BD. It holds the row numbers of the packed part of L/U, and
C     the row numbers of the full part of L/U.
C     It is not altered by the subroutine.
C IPTRL is an integer array that must be unchanged since calling
C     MA50B/BD. For J = 1,..., NP, IPTRL(J) holds the position in
C     FACT and IRNF of the end of column J of L.
C     For J = NP+1,..., N, IPTRL(J) is equal to IPTRU(J).
C     It is not altered by the subroutine.
C IPTRU is an integer array that must be unchanged since calling
C     MA50B/BD. For J = 1,..., N, IPTRU(J) holds the position in
C     FACT and IRNF of the end of the packed part of column J of U.
C     It is not altered by the subroutine.
C B is an array that must be set to the vector b.
C     It is not altered.
C X is an array that need not be set on entry. On return, it holds the
C    solution x.
C W is a work array of length max(M,N).
C INFO need not be set on entry. On exit, it holds the following:
C    INFO(1) A nonzero value will indicate an error return. Possible
C      nonzero values are:
C      -1  M < 1 or N < 1

      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0D0)
      INTEGER I,II,IA1,IF1,J,LP,MF,MP,NF
      DOUBLE PRECISION PROD
C I Temporary variable holding row number.
C II Position of the current entry in IRNF.
C IA1 Position of the start of the current row or column.
C IF1 Position of the start of the full part of U.
C J Temporary variable holding column number.
C LP Unit for error messages.
C MF Number of rows held in full format.
C MP Unit for diagnostic messages.
C NF Number of columns held in full format.
C PROD Temporary variable used to accumulate inner products.

      EXTERNAL MA50HD

      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (ICNTL(3).LE.0) LP = 0
      IF (ICNTL(3).LE.1) MP = 0

C Make some simple checks
      IF (M.LT.1 .OR. N.LT.1) GO TO 250

      IF (MP.GT.0 .AND. ICNTL(3).GT.2) WRITE (MP,
     +    '(/2(A,I6),A,I4,A,L2)') ' Entering MA50CD with M=',M,' N =',N,
     +    ' NP =',NP,' TRANS =',TRANS
      IF1 = IPTRL(N) + 1
      MF = IRNF(2)
      NF = N - NP
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) WRITE (MP,
     +    '(A,I5,A,I5)') ' Size of full submatrix',MF,' by',NF
      IF (MP.GT.0 .AND. ICNTL(3).GT.3) THEN
         DO 10 J = 1,N
            IF (J.GT.1) THEN
               IF (IPTRL(J-1).LT.IPTRU(J)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of U',
     +              (FACT(II),IRNF(II),II=IPTRL(J-1)+1,IPTRU(J))
            END IF
            IF (IPTRU(J).LT.IPTRL(J)) WRITE (MP,
     +          '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of L',
     +          (FACT(II),IRNF(II),II=IPTRU(J)+1,IPTRL(J))
   10    CONTINUE
         WRITE (MP,'(A)') ' Full part'
         WRITE (MP,'((6I12))') (IRNF(IF1+MF+J),J=0,NF-1)
         DO 20 I = 0,MF - 1
            WRITE (MP,'(I4,1P,6E12.4:/(4X,1P,6E12.4))') IRNF(IF1+I),
     +        (FACT(IF1+I+J*MF),J=0,NF-1)
   20    CONTINUE
      END IF

      IF (TRANS) THEN
         IF (MP.GT.0 .AND. ICNTL(3).GT.3) WRITE (MP,
     +       '(A4,5F10.4:/(4X,5F10.4))') ' B =', (B(I),I=1,N)
         IF (IQ(1).GT.0) THEN
            DO 30 I = 1,N
               W(I) = B(IQ(I))
   30       CONTINUE
         ELSE
            DO 40 I = 1,N
               W(I) = B(I)
   40       CONTINUE
         END IF
         DO 50 I = 1,M
            X(I) = ZERO
   50    CONTINUE
C Forward substitution through packed part of (U trans).
         DO 70 I = 2,N
            PROD = ZERO
            DO 60 II = IPTRL(I-1) + 1,IPTRU(I)
               PROD = PROD + FACT(II)*W(IRNF(II))
   60       CONTINUE
            W(I) = W(I) + PROD
   70    CONTINUE
C Backsubstitute through the full part of (PL) trans.
         DO 80 I = 1,NF
            X(I) = W(NP+I)
   80    CONTINUE
         IF (MF.GT.0 .AND. NF.GT.0) THEN
            CALL MA50HD(TRANS,MF,NF,FACT(IF1),MF,IRNF(IF1+MF),X,
     +                  ICNTL(5))
         ELSE
            DO 90 I = 1,MF
               X(I) = ZERO
   90       CONTINUE
         END IF
         DO 100 I = MF,1,-1
            J = IRNF(IF1+I-1)
            IF (J.NE.I) X(J) = X(I)
  100    CONTINUE
C Backsubstitute through the packed part of (PL) trans.
         DO 120 I = NP,1,-1
            IA1 = IPTRU(I) + 1
            IF (IA1.GT.IPTRL(I)) GO TO 120
            PROD = ZERO
            DO 110 II = IA1 + 1,IPTRL(I)
               PROD = PROD + FACT(II)*X(IRNF(II))
  110       CONTINUE
            X(IRNF(IA1)) = (W(I)-PROD)*FACT(IA1)
  120    CONTINUE
         IF (MP.GT.0 .AND. ICNTL(3).GT.3) WRITE (MP,
     +       '(A/(4X,5F10.4))') ' Leaving MA50CD with X =', (X(I),I=1,M)
C
      ELSE
         IF (MP.GT.0 .AND. ICNTL(3).GT.3) WRITE (MP,
     +       '(A4,5F10.4:/(4X,5F10.4))') ' B =', (B(I),I=1,M)
C Forward substitution through the packed part of PL
         DO 130 I = 1,M
            W(I) = B(I)
  130    CONTINUE
         DO 150 I = 1,NP
            IA1 = IPTRU(I) + 1
            IF (IA1.LE.IPTRL(I)) THEN
               X(I) = W(IRNF(IA1))*FACT(IA1)
               IF (X(I).NE.ZERO) THEN
CDIR$ IVDEP
                  DO 140 II = IA1 + 1,IPTRL(I)
                     W(IRNF(II)) = W(IRNF(II)) - FACT(II)*X(I)
  140             CONTINUE
               END IF
            END IF
  150    CONTINUE
C Forward substitution through the full part of PL
         IF (MF.GT.0 .AND. NF.GT.0) THEN
            DO 160 I = 1,MF
               W(I) = W(IRNF(IF1+I-1))
  160       CONTINUE
            CALL MA50HD(TRANS,MF,NF,FACT(IF1),MF,IRNF(IF1+MF),W,
     +                  ICNTL(5))
            DO 170 I = 1,NF
               X(NP+I) = W(I)
  170       CONTINUE
         ELSE
            DO 180 I = 1,NF
               X(NP+I) = ZERO
  180       CONTINUE
         END IF
C Back substitution through the packed part of U
         DO 200 J = N,MAX(2,NP+1),-1
            PROD = X(J)
CDIR$ IVDEP
            DO 190 II = IPTRL(J-1) + 1,IPTRU(J)
               X(IRNF(II)) = X(IRNF(II)) + FACT(II)*PROD
  190       CONTINUE
  200    CONTINUE
         DO 220 J = NP,2,-1
            IA1 = IPTRU(J)
            IF (IA1.GE.IPTRL(J)) THEN
               X(J) = ZERO
            ELSE
               PROD = X(J)
CDIR$ IVDEP
               DO 210 II = IPTRL(J-1) + 1,IA1
                  X(IRNF(II)) = X(IRNF(II)) + FACT(II)*PROD
  210          CONTINUE
            END IF
  220    CONTINUE
         IF (NP.GE.1 .AND. IPTRU(1).GE.IPTRL(1)) X(1) = ZERO
         IF (IQ(1).GT.0) THEN
C         Permute X
            DO 230 I = 1,N
               W(I) = X(I)
  230       CONTINUE
            DO 240 I = 1,N
               X(IQ(I)) = W(I)
  240       CONTINUE
         END IF
         IF (MP.GT.0 .AND. ICNTL(3).GT.3) WRITE (MP,
     +       '(A/(4X,5F10.4))') ' Leaving MA50CD with X =', (X(I),I=1,N)
      END IF
      RETURN
C Error condition.
  250 INFO(1) = -1
      IF (LP.GT.0) WRITE (LP,'(/A/2(A,I8))')
     +    ' **** Error return from MA50CD ****',' M =',M,' N =',N
      END

      SUBROUTINE MA50DD(LA,A,IND,IPTR,N,DISP,REALS)
C This subroutine performs garbage collection on the arrays A and IND.
C DISP is the position in arrays A/IND immediately after the data
C     to be compressed.
C     On exit, DISP equals the position of the first entry
C     after the compressed part of A/IND.
C
      INTEGER LA,N,DISP
      DOUBLE PRECISION A(LA)
      INTEGER IPTR(N)
      LOGICAL REALS
      INTEGER IND(LA)
C Local variables.
      INTEGER J,K,KN
C Set the first entry in each row(column) to the negative of the
C     row(column) and hold the column(row) index in the row(column)
C     pointer.  This enables the start of each row(column) to be
C     recognized in a subsequent scan.
      DO 10 J = 1,N
         K = IPTR(J)
         IF (K.GT.0) THEN
            IPTR(J) = IND(K)
            IND(K) = -J
         END IF
   10 CONTINUE
      KN = 0
C Go through arrays compressing to the front so that there are no
C     zeros held in positions 1 to DISP-1 of IND.
C     Reset first entry of each row(column) and the pointer array IPTR.
      DO 20 K = 1,DISP - 1
         IF (IND(K).EQ.0) GO TO 20
         KN = KN + 1
         IF (REALS) A(KN) = A(K)
         IF (IND(K).LE.0) THEN
C First entry of row(column) has been located.
            J = -IND(K)
            IND(K) = IPTR(J)
            IPTR(J) = KN
         END IF
         IND(KN) = IND(K)
   20 CONTINUE
      DISP = KN + 1
      END


      SUBROUTINE MA50ED(M,N,A,LDA,PIVTOL,IPIV,RANK)
**
      INTEGER LDA,M,N,RANK
      DOUBLE PRECISION PIVTOL

      INTEGER IPIV(N)
      DOUBLE PRECISION A(LDA,N)

*
*  Purpose
*  =======
*
*  MA50ED computes an LU factorization of a general m-by-n matrix A.

*  The factorization has the form
*     A = P * L * U * Q
*  where P is a permutation matrix of order m, L is lower triangular
*  of order m with unit diagonal elements, U is upper trapezoidal of
*  order m * n, and Q is a permutation matrix of order n.
*
*  Row interchanges are used to ensure that the entries of L do not
*  exceed 1 in absolute value. Column interchanges are used to
*  ensure that the first r diagonal entries of U exceed PIVTOL in
*  absolute value. If r < m, the last (m-r) rows of U are zero.

*  This is the Level 1 BLAS version.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 1.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 1.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U*Q; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  PIVTOL  (input) DOUBLE PRECISION
*          The pivot tolerance. Any entry with absolute value less
*          than or equal to PIVTOL is regarded as unsuitable to be a
*          pivot.
**
*  IPIV    (output) INTEGER array, dimension (N)
*          The permutations; for 1 <= i <= RANK, row i of the
*          matrix was interchanged with row IPIV(i); for
*          RANK + 1 <= j <= N, column j of the
*          matrix was interchanged with column -IPIV(j).
*
*  RANK    (output) INTEGER
*          The computed rank of the matrix.
*
*  =====================================================================
*
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)

      INTEGER I,J,JP,K
      LOGICAL PIVOT
* I   Row index.
* J   Current column.
* JP  Pivot position.
* K   Main loop index.
* PIVOT True if there is a pivot in current column.

      INTEGER IDAMAX
      EXTERNAL IDAMAX

      EXTERNAL DAXPY,DSCAL,DSWAP
      INTRINSIC ABS

*
      J = 1
      DO 30 K = 1,N

*        Update elements in column J.
         DO 10 I = 1,J - 1
            IF (M.GT.I) CALL DAXPY(M-I,-A(I,J),A(I+1,I),1,A(I+1,J),1)
   10    CONTINUE

*        Find pivot.
         IF (J.LE.M) THEN
            JP = J - 1 + IDAMAX(M-J+1,A(J,J),1)
            IPIV(J) = JP
            PIVOT = ABS(A(JP,J)) .GT. PIVTOL
         ELSE
            PIVOT = .FALSE.
         END IF
         IF (PIVOT) THEN

*           Apply row interchange to columns 1:N+J-K.
            IF (JP.NE.J) CALL DSWAP(N+J-K,A(J,1),LDA,A(JP,1),LDA)
*
*           Compute elements J+1:M of J-th column.
            IF (J.LT.M) CALL DSCAL(M-J,ONE/A(J,J),A(J+1,J),1)

*           Update J
            J = J + 1
*
         ELSE
*
            DO 20 I = J,M
               A(I,J) = ZERO
   20       CONTINUE
*           Apply column interchange and record it.
            IF (K.LT.N) CALL DSWAP(M,A(1,J),1,A(1,N-K+J),1)
            IPIV(N-K+J) = -J
*
         END IF
*
   30 CONTINUE

      RANK = J - 1
*
*     End of MA50ED
*
      END


      SUBROUTINE MA50FD(M,N,A,LDA,PIVTOL,IPIV,RANK)
*
*  -- This is a variant of the LAPACK routine DGETF2 --
*
      INTEGER LDA,M,N,RANK
      DOUBLE PRECISION PIVTOL

      INTEGER IPIV(N)
      DOUBLE PRECISION A(LDA,N)

*
*  Purpose
*  =======
*
*  MA50FD computes an LU factorization of a general m-by-n matrix A.

*  The factorization has the form
*     A = P * L * U * Q
*  where P is a permutation matrix of order m, L is lower triangular
*  of order m with unit diagonal elements, U is upper trapezoidal of
*  order m * n, and Q is a permutation matrix of order n.
*
*  Row interchanges are used to ensure that the entries of L do not
*  exceed 1 in absolute value. Column interchanges are used to
*  ensure that the first r diagonal entries of U exceed PIVTOL in
*  absolute value. If r < m, the last (m-r) rows of U are zero.

*  This is the Level 2 BLAS version.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 1.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 1.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U*Q; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  PIVTOL  (input) DOUBLE PRECISION
*          The pivot tolerance. Any entry with absolute value less
*          than or equal to PIVTOL is regarded as unsuitable to be a
*          pivot.
**
*  IPIV    (output) INTEGER array, dimension (N)
*          The permutations; for 1 <= i <= RANK, row i of the
*          matrix was interchanged with row IPIV(i); for
*          RANK + 1 <= j <= N, column j of the
*          matrix was interchanged with column -IPIV(j).
*
*  RANK    (output) INTEGER
*          The computed rank of the matrix.
*
*  =====================================================================
*
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)

      INTEGER I,J,JP,K
      LOGICAL PIVOT
* I   Row index.
* J   Current column.
* JP  Pivot position.
* K   Main loop index.
* PIVOT True if there is a pivot in current column.

      INTEGER IDAMAX
      EXTERNAL IDAMAX

      EXTERNAL DGEMV,DSCAL,DSWAP

      INTRINSIC ABS

*
      J = 1
      DO 20 K = 1,N

         IF (J.LE.M) THEN
*           Update diagonal and subdiagonal elements in column J.
            CALL DGEMV('No transpose',M-J+1,J-1,-ONE,A(J,1),LDA,A(1,J),
     +                 1,ONE,A(J,J),1)
*          Find pivot.
            JP = J - 1 + IDAMAX(M-J+1,A(J,J),1)
            IPIV(J) = JP
            PIVOT = ABS(A(JP,J)) .GT. PIVTOL
         ELSE
            PIVOT = .FALSE.
         END IF

         IF (PIVOT) THEN

*           Apply row interchange to columns 1:N+J-K.
            IF (JP.NE.J) CALL DSWAP(N+J-K,A(J,1),LDA,A(JP,1),LDA)
*
*           Compute elements J+1:M of J-th column.
            IF (J.LT.M) CALL DSCAL(M-J,ONE/A(J,J),A(J+1,J),1)

            IF (J.LT.N) THEN
*             Compute block row of U.
               CALL DGEMV('Transpose',J-1,N-J,-ONE,A(1,J+1),LDA,A(J,1),
     +                    LDA,ONE,A(J,J+1),LDA)
            END IF

*           Update J
            J = J + 1
*
         ELSE
*
            DO 10 I = J,M
               A(I,J) = ZERO
   10       CONTINUE
*           Apply column interchange and record it.
            IF (K.LT.N) CALL DSWAP(M,A(1,J),1,A(1,N-K+J),1)
            IPIV(N-K+J) = -J
*
         END IF
*
   20 CONTINUE

      RANK = J - 1
*
*     End of MA50FD
*
      END


      SUBROUTINE MA50GD(M,N,A,LDA,NB,PIVTOL,IPIV,RANK)
*
*  -- This is a variant of the LAPACK routine DGETRF --
*
      INTEGER LDA,M,N,NB,RANK
      DOUBLE PRECISION PIVTOL

      INTEGER IPIV(N)
      DOUBLE PRECISION A(LDA,N)

*
*  Purpose
*  =======
*
*  MA50GD computes an LU factorization of a general m-by-n matrix A.
*
*  The factorization has the form
*     A = P * L * U * Q
*  where P is a permutation matrix of order m, L is lower triangular
*  of order m with unit diagonal elements, U is upper trapezoidal of
*  order m * n, and Q is a permutation matrix of order n.
*
*  Row interchanges are used to ensure that the entries of L do not
*  exceed 1 in absolute value. Column interchanges are used to
*  ensure that the first r diagonal entries of U exceed PIVTOL in
*  absolute value. If r < m, the last (m-r) rows of U are zero.
*
*  This is the Level 3 BLAS version.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 1.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 1.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U*Q; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  NB      (input) INTEGER
*          The block size for BLAS3 processing.
*
*  PIVTOL  (input) DOUBLE PRECISION
*          The pivot tolerance. Any entry with absolute value less
*          than or equal to PIVTOL is regarded as unsuitable to be a
*          pivot.
**
*  IPIV    (output) INTEGER array, dimension (N)
*          The permutations; for 1 <= i <= RANK, row i of the
*          matrix was interchanged with row IPIV(i); for
*          RANK + 1 <= j <= N, column j of the
*          matrix was interchanged with column -IPIV(j).
*
*  RANK    (output) INTEGER
*          The computed rank of the matrix.
*
*  =====================================================================
*
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)

      INTEGER I,J,JJ,JP,J1,J2,K
      LOGICAL PIVOT
      DOUBLE PRECISION TEMP

* I   DO index for applying permutations.
* J   Current column.
* JJ  Column in which swaps occur.
* JP  Pivot position.
* J1  Column at start of current block.
* J2  Column at end of current block.
* K   Main loop index.
* PIVOT True if there is a pivot in current column.
* TEMP Temporary variable for swaps.

      EXTERNAL DGEMM,DGEMV,DSWAP,DSCAL,DTRSM,DTRSV

      INTEGER IDAMAX
      EXTERNAL IDAMAX

      INTRINSIC ABS,MIN

*
      J = 1
      J1 = 1
      J2 = MIN(N,NB)
      DO 70 K = 1,N

         IF (J.LE.M) THEN

*          Update diagonal and subdiagonal elements in column J.
            CALL DGEMV('No transpose',M-J+1,J-J1,-ONE,A(J,J1),LDA,
     +                 A(J1,J),1,ONE,A(J,J),1)

*          Find pivot.
            JP = J - 1 + IDAMAX(M-J+1,A(J,J),1)
            IPIV(J) = JP
            PIVOT = ABS(A(JP,J)) .GT. PIVTOL
         ELSE
            PIVOT = .FALSE.
         END IF

         IF (PIVOT) THEN

*           Apply row interchange to columns J1:J2
            IF (JP.NE.J) CALL DSWAP(J2-J1+1,A(J,J1),LDA,A(JP,J1),LDA)
*
*           Compute elements J+1:M of J-th column.
            IF (J.LT.M) CALL DSCAL(M-J,ONE/A(J,J),A(J+1,J),1)

            IF (J+1.LE.J2) THEN
*             Compute row of U within current block
               CALL DGEMV('Transpose',J-J1,J2-J,-ONE,A(J1,J+1),LDA,
     +                    A(J,J1),LDA,ONE,A(J,J+1),LDA)
            END IF

*           Update J
            J = J + 1
*
         ELSE

            DO 10 I = J,M
               A(I,J) = ZERO
   10       CONTINUE
*
*           Record column interchange and revise J2 if necessary
            IPIV(N-K+J) = -J
*           Apply column interchange.
            IF (K.NE.N) CALL DSWAP(M,A(1,J),1,A(1,N-K+J),1)
            IF (N-K+J.GT.J2) THEN
*              Apply operations to new column.
               DO 20 I = J1,J - 1
                  JP = IPIV(I)
                  TEMP = A(I,J)
                  A(I,J) = A(JP,J)
                  A(JP,J) = TEMP
   20          CONTINUE
               IF(J.GT.J1) CALL DTRSV('Lower','No transpose','Unit',
     +                                J-J1,A(J1,J1),LDA,A(J1,J),1)
            ELSE
               J2 = J2 - 1
            END IF
*
         END IF

         IF (J.GT.J2) THEN
*           Apply permutations to columns outside the block
            DO 40 JJ = 1,J1 - 1
               DO 30 I = J1,J2
                  JP = IPIV(I)
                  TEMP = A(I,JJ)
                  A(I,JJ) = A(JP,JJ)
                  A(JP,JJ) = TEMP
   30          CONTINUE
   40       CONTINUE
            DO 60 JJ = J2 + 1,N - K + J - 1
               DO 50 I = J1,J2
                  JP = IPIV(I)
                  TEMP = A(I,JJ)
                  A(I,JJ) = A(JP,JJ)
                  A(JP,JJ) = TEMP
   50          CONTINUE
   60       CONTINUE

            IF (K.NE.N) THEN
*              Update the Schur complement
               CALL DTRSM('Left','Lower','No transpose','Unit',J2-J1+1,
     +                    N-K,ONE,A(J1,J1),LDA,A(J1,J2+1),LDA)
               IF (M.GT.J2) CALL DGEMM('No transpose','No transpose',
     +                                 M-J2,N-K,J2-J1+1,-ONE,A(J2+1,J1),
     +                                 LDA,A(J1,J2+1),LDA,ONE,
     +                                 A(J2+1,J2+1),LDA)
            END IF

            J1 = J2 + 1
            J2 = MIN(J2+NB,N-K+J-1)

         END IF

*
   70 CONTINUE
      RANK = J - 1
*
*     End of MA50GD
*
      END

      SUBROUTINE MA50HD(TRANS,M,N,A,LDA,IPIV,B,ICNTL5)
*
*  -- This is a variant of the LAPACK routine DGETRS --
*     It handles the singular or rectangular case.
*
      LOGICAL TRANS
      INTEGER LDA,M,N,ICNTL5
      INTEGER IPIV(N)
      DOUBLE PRECISION A(LDA,N),B(*)

*
*  Purpose
*  =======
*
*  Solve a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general m by n matrix A using the LU factorization computed
*  by MA50DE, MA50FD, or MA50GD.
*
*  Arguments
*  =========
*
*  TRANS   (input) LOGICAL
*          Specifies the form of the system of equations.
*          = .FALSE. :  A * X = B  (No transpose)
*          = .TRUE.  :  A'* X = B  (Transpose)
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 1.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 1.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by MA50GD.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The permutations; for 1 <= i <= RANK, row i of the
*          matrix was interchanged with row IPIV(i); for
*          RANK + 1 <= j <= N, column j of the
*          matrix was interchanged with column -IPIV(j).

*  B       (input/output) DOUBLE PRECISION array, size max(M,N)
*          On entry, the right hand side vectors B for the system of
*          linear equations.
*          On exit, the solution vectors, X.
*
*  ICNTL5  (input) INTEGER
*          0 for BLAS1 or >0 for BLAS2
*
*  =====================================================================
*
      INTEGER I,K,RANK
C I    Temporary variable.
C K    Temporary variable.
C RANK Rank of matrix.

      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0D0)
      DOUBLE PRECISION TEMP
      INTRINSIC MIN
      EXTERNAL DAXPY,DDOT,DTRSV
      DOUBLE PRECISION DDOT

*   Find the rank
      RANK = 0
      DO 10 RANK = MIN(M,N),1,-1
         IF (IPIV(RANK).GT.0) GO TO 20
   10 CONTINUE

   20 IF (.NOT.TRANS) THEN
*
*        Solve A * X = B.
*
*        Apply row interchanges to the right hand side.
         DO 30 I = 1,RANK
            K = IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
   30    CONTINUE
*
*        Solve L*X = B, overwriting B with X.
         IF (ICNTL5.GT.0) THEN
            IF (RANK.GT.0) CALL DTRSV('L','NoTrans','Unit',RANK,A,LDA,B,
     +                                1)
         ELSE
            DO 40 K = 1,RANK - 1
               IF (B(K).NE.ZERO) CALL DAXPY(RANK-K,-B(K),A(K+1,K),1,
     +                                B(K+1),1)
   40       CONTINUE
         END IF

*        Solve U*X = B, overwriting B with X.
         IF (ICNTL5.GT.0) THEN
            IF (RANK.GT.0) CALL DTRSV('U','NoTrans','NonUnit',RANK,A,
     +                                LDA,B,1)
         ELSE
            DO 50 K = RANK,2,-1
               IF (B(K).NE.ZERO) THEN
                  B(K) = B(K)/A(K,K)
                  CALL DAXPY(K-1,-B(K),A(1,K),1,B(1),1)
               END IF
   50       CONTINUE
            IF (RANK.GT.0) B(1) = B(1)/A(1,1)
         END IF

*        Set singular part to zero
         DO 60 K = RANK + 1,N
            B(K) = ZERO
   60    CONTINUE
*
*        Apply column interchanges to the right hand side.
         DO 70 I = RANK + 1,N
            K = -IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
   70    CONTINUE

      ELSE
*
*        Solve A' * X = B.

*        Apply column interchanges to the right hand side.
         DO 80 I = N,RANK + 1,-1
            K = -IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
   80    CONTINUE
*
*        Solve U'*X = B, overwriting B with X.
*
         IF (ICNTL5.GT.0) THEN
            IF (RANK.GT.0) CALL DTRSV('U','Trans','NonUnit',RANK,A,LDA,
     +                                B,1)
         ELSE
            IF (RANK.GT.0) B(1) = B(1)/A(1,1)
            DO 90 I = 2,RANK
               TEMP = B(I) - DDOT(I-1,A(1,I),1,B(1),1)
               B(I) = TEMP/A(I,I)
   90       CONTINUE
         END IF

*        Solve L'*X = B, overwriting B with X.
         IF (ICNTL5.GT.0) THEN
            IF (RANK.GT.0) CALL DTRSV('L','Trans','Unit',RANK,A,LDA,B,1)
         ELSE
            DO 100 I = RANK - 1,1,-1
               B(I) = B(I) - DDOT(RANK-I,A(I+1,I),1,B(I+1),1)
  100       CONTINUE
         END IF

*        Set singular part to zero
         DO 110 I = RANK + 1,M
            B(I) = ZERO
  110    CONTINUE
*
*        Apply row interchanges to the solution vectors.
         DO 120 I = RANK,1,-1
            K = IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
  120    CONTINUE
      END IF

      END

      SUBROUTINE MA50ID(CNTL,ICNTL)
C Set default values for the control arrays.

      DOUBLE PRECISION CNTL(10)
      INTEGER I,ICNTL(20)

      CNTL(1) = 0.5D0
      CNTL(2) = 0.1D0
      DO 10 I = 3,10
         CNTL(I) = 0.0D0
   10 CONTINUE

      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = 1
      ICNTL(4) = 3
      ICNTL(5) = 32
      DO 20 I = 6,20
        ICNTL(I) = 0
   20 CONTINUE

      END
* COPYRIGHT (c) 1976 AEA Technology
* Original date 21 Jan 1993
C       Toolpack tool decs employed.
C	Double version of MC13D (name change only)
C 10 August 2001 DOs terminated with CONTINUE
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC13DD(N,ICN,LICN,IP,LENR,IOR,IB,NUM,IW)
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUM
C     ..
C     .. Array Arguments ..
      INTEGER IB(N),ICN(LICN),IOR(N),IP(N),IW(N,3),LENR(N)
C     ..
C     .. External Subroutines ..
      EXTERNAL MC13ED
C     ..
C     .. Executable Statements ..
      CALL MC13ED(N,ICN,LICN,IP,LENR,IOR,IB,NUM,IW(1,1),IW(1,2),IW(1,3))
      RETURN

      END
      SUBROUTINE MC13ED(N,ICN,LICN,IP,LENR,ARP,IB,NUM,LOWL,NUMB,PREV)
C
C ARP(I) IS ONE LESS THAN THE NUMBER OF UNSEARCHED EDGES LEAVING
C     NODE I.  AT THE END OF THE ALGORITHM IT IS SET TO A
C     PERMUTATION WHICH PUTS THE MATRIX IN BLOCK LOWER
C     TRIANGULAR FORM.
C IB(I) IS THE POSITION IN THE ORDERING OF THE START OF THE ITH
C     BLOCK.  IB(N+1-I) HOLDS THE NODE NUMBER OF THE ITH NODE
C     ON THE STACK.
C LOWL(I) IS THE SMALLEST STACK POSITION OF ANY NODE TO WHICH A PATH
C     FROM NODE I HAS BEEN FOUND.  IT IS SET TO N+1 WHEN NODE I
C     IS REMOVED FROM THE STACK.
C NUMB(I) IS THE POSITION OF NODE I IN THE STACK IF IT IS ON
C     IT, IS THE PERMUTED ORDER OF NODE I FOR THOSE NODES
C     WHOSE FINAL POSITION HAS BEEN FOUND AND IS OTHERWISE ZERO.
C PREV(I) IS THE NODE AT THE END OF THE PATH WHEN NODE I WAS
C     PLACED ON THE STACK.
C
C
C   ICNT IS THE NUMBER OF NODES WHOSE POSITIONS IN FINAL ORDERING HAVE
C     BEEN FOUND.
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUM
C     ..
C     .. Array Arguments ..
      INTEGER ARP(N),IB(N),ICN(LICN),IP(N),LENR(N),LOWL(N),NUMB(N),
     +        PREV(N)
C     ..
C     .. Local Scalars ..
      INTEGER DUMMY,I,I1,I2,ICNT,II,ISN,IST,IST1,IV,IW,J,K,LCNT,NNM1,STP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Executable Statements ..
      ICNT = 0
C NUM IS THE NUMBER OF BLOCKS THAT HAVE BEEN FOUND.
      NUM = 0
      NNM1 = N + N - 1
C
C INITIALIZATION OF ARRAYS.
      DO 20 J = 1,N
        NUMB(J) = 0
        ARP(J) = LENR(J) - 1
   20 CONTINUE
C
C
      DO 120 ISN = 1,N
C LOOK FOR A STARTING NODE
        IF (NUMB(ISN).NE.0) GO TO 120
        IV = ISN
C IST IS THE NUMBER OF NODES ON THE STACK ... IT IS THE STACK POINTER.
        IST = 1
C PUT NODE IV AT BEGINNING OF STACK.
        LOWL(IV) = 1
        NUMB(IV) = 1
        IB(N) = IV
C
C THE BODY OF THIS LOOP PUTS A NEW NODE ON THE STACK OR BACKTRACKS.
        DO 110 DUMMY = 1,NNM1
          I1 = ARP(IV)
C HAVE ALL EDGES LEAVING NODE IV BEEN SEARCHED.
          IF (I1.LT.0) GO TO 60
          I2 = IP(IV) + LENR(IV) - 1
          I1 = I2 - I1
C
C LOOK AT EDGES LEAVING NODE IV UNTIL ONE ENTERS A NEW NODE OR
C     ALL EDGES ARE EXHAUSTED.
          DO 50 II = I1,I2
            IW = ICN(II)
C HAS NODE IW BEEN ON STACK ALREADY.
            IF (NUMB(IW).EQ.0) GO TO 100
C UPDATE VALUE OF LOWL(IV) IF NECESSARY.
            LOWL(IV) = MIN(LOWL(IV),LOWL(IW))
   50     CONTINUE
C
C THERE ARE NO MORE EDGES LEAVING NODE IV.
          ARP(IV) = -1
C IS NODE IV THE ROOT OF A BLOCK.
   60     IF (LOWL(IV).LT.NUMB(IV)) GO TO 90
C
C ORDER NODES IN A BLOCK.
          NUM = NUM + 1
          IST1 = N + 1 - IST
          LCNT = ICNT + 1
C PEEL BLOCK OFF THE TOP OF THE STACK STARTING AT THE TOP AND
C     WORKING DOWN TO THE ROOT OF THE BLOCK.
          DO 70 STP = IST1,N
            IW = IB(STP)
            LOWL(IW) = N + 1
            ICNT = ICNT + 1
            NUMB(IW) = ICNT
            IF (IW.EQ.IV) GO TO 80
   70     CONTINUE
   80     IST = N - STP
          IB(NUM) = LCNT
C ARE THERE ANY NODES LEFT ON THE STACK.
          IF (IST.NE.0) GO TO 90
C HAVE ALL THE NODES BEEN ORDERED.
          IF (ICNT.LT.N) GO TO 120
          GO TO 130
C
C BACKTRACK TO PREVIOUS NODE ON PATH.
   90     IW = IV
          IV = PREV(IV)
C UPDATE VALUE OF LOWL(IV) IF NECESSARY.
          LOWL(IV) = MIN(LOWL(IV),LOWL(IW))
          GO TO 110
C
C PUT NEW NODE ON THE STACK.
  100     ARP(IV) = I2 - II - 1
          PREV(IW) = IV
          IV = IW
          IST = IST + 1
          LOWL(IV) = IST
          NUMB(IV) = IST
          K = N + 1 - IST
          IB(K) = IV
  110   CONTINUE
C
  120 CONTINUE
C
C
C PUT PERMUTATION IN THE REQUIRED FORM.
  130 DO 140 I = 1,N
        II = NUMB(I)
        ARP(II) = I
  140 CONTINUE
      RETURN

      END
* COPYRIGHT (c) 1977 AEA Technology
* Original date 8 Oct 1992
C######8/10/92 Toolpack tool decs employed.
C######8/10/92 D version created by name change only.
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC21AD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUMNZ
C     ..
C     .. Array Arguments ..
      INTEGER ICN(LICN),IP(N),IPERM(N),IW(N,4),LENR(N)
C     ..
C     .. External Subroutines ..
      EXTERNAL MC21BD
C     ..
C     .. Executable Statements ..
      CALL MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2),
     +            IW(1,3),IW(1,4))
      RETURN
C
      END
      SUBROUTINE MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,OUT)
C   PR(I) IS THE PREVIOUS ROW TO I IN THE DEPTH FIRST SEARCH.
C IT IS USED AS A WORK ARRAY IN THE SORTING ALGORITHM.
C   ELEMENTS (IPERM(I),I) I=1, ... N  ARE NON-ZERO AT THE END OF THE
C ALGORITHM UNLESS N ASSIGNMENTS HAVE NOT BEEN MADE.  IN WHICH CASE
C (IPERM(I),I) WILL BE ZERO FOR N-NUMNZ ENTRIES.
C   CV(I) IS THE MOST RECENT ROW EXTENSION AT WHICH COLUMN I
C WAS VISITED.
C   ARP(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C WHICH HAVE NOT BEEN SCANNED WHEN LOOKING FOR A CHEAP ASSIGNMENT.
C   OUT(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C WHICH HAVE NOT BEEN SCANNED DURING ONE PASS THROUGH THE MAIN LOOP.
C
C   INITIALIZATION OF ARRAYS.
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUMNZ
C     ..
C     .. Array Arguments ..
      INTEGER ARP(N),CV(N),ICN(LICN),IP(N),IPERM(N),LENR(N),OUT(N),PR(N)
C     ..
C     .. Local Scalars ..
      INTEGER I,II,IN1,IN2,IOUTK,J,J1,JORD,K,KK
C     ..
C     .. Executable Statements ..
      DO 10 I = 1,N
        ARP(I) = LENR(I) - 1
        CV(I) = 0
        IPERM(I) = 0
   10 CONTINUE
      NUMNZ = 0
C
C
C   MAIN LOOP.
C   EACH PASS ROUND THIS LOOP EITHER RESULTS IN A NEW ASSIGNMENT
C OR GIVES A ROW WITH NO ASSIGNMENT.
      DO 100 JORD = 1,N
        J = JORD
        PR(J) = -1
        DO 70 K = 1,JORD
C LOOK FOR A CHEAP ASSIGNMENT
          IN1 = ARP(J)
          IF (IN1.LT.0) GO TO 30
          IN2 = IP(J) + LENR(J) - 1
          IN1 = IN2 - IN1
          DO 20 II = IN1,IN2
            I = ICN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
C   NO CHEAP ASSIGNMENT IN ROW.
          ARP(J) = -1
C   BEGIN LOOKING FOR ASSIGNMENT CHAIN STARTING WITH ROW J.
   30     CONTINUE
          OUT(J) = LENR(J) - 1
C INNER LOOP.  EXTENDS CHAIN BY ONE OR BACKTRACKS.
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENR(J) - 1
            IN1 = IN2 - IN1
C FORWARD SCAN.
            DO 40 II = IN1,IN2
              I = ICN(II)
              IF (CV(I).EQ.JORD) GO TO 40
C   COLUMN I HAS NOT YET BEEN ACCESSED DURING THIS PASS.
              J1 = J
              J = IPERM(I)
              CV(I) = JORD
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
C
   40       CONTINUE
C
C   BACKTRACKING STEP.
   50       CONTINUE
            J = PR(J)
            IF (J.EQ.-1) GO TO 100
   60     CONTINUE
C
   70   CONTINUE
C
C   NEW ASSIGNMENT IS MADE.
   80   CONTINUE
        IPERM(I) = J
        ARP(J) = IN2 - II - 1
        NUMNZ = NUMNZ + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 100
          II = IP(J) + LENR(J) - OUT(J) - 2
          I = ICN(II)
          IPERM(I) = J
   90   CONTINUE
C
  100 CONTINUE
C
C   IF MATRIX IS STRUCTURALLY SINGULAR, WE NOW COMPLETE THE
C PERMUTATION IPERM.
      IF (NUMNZ.EQ.N) RETURN
      DO 110 I = 1,N
        ARP(I) = 0
  110 CONTINUE
      K = 0
      DO 130 I = 1,N
        IF (IPERM(I).NE.0) GO TO 120
        K = K + 1
        OUT(K) = I
        GO TO 130
C
  120   CONTINUE
        J = IPERM(I)
        ARP(J) = I
  130 CONTINUE
      K = 0
      DO 140 I = 1,N
        IF (ARP(I).NE.0) GO TO 140
        K = K + 1
        IOUTK = OUT(K)
        IPERM(IOUTK) = I
  140 CONTINUE
      RETURN
C
      END
* COPYRIGHT (c) 1993 AEA Technology and
* Council for the Central Laboratory of the Research Councils
C Original date March 1993
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC29AD(M,N,NE,A,IRN,ICN,R,C,W,LP,IFAIL)
      INTEGER M,N,NE
      DOUBLE PRECISION A(NE)
      INTEGER IRN(NE),ICN(NE)
      DOUBLE PRECISION R(M),C(N),W(M*2+N*3)
      INTEGER LP,IFAIL
C M is an integer variable that must be set to the number of rows.
C      It is not altered by the subroutine.
C N is an integer variable that must be set to the number of columns.
C      It is not altered by the subroutine.
C NE is an integer variable that must be set to the number of entries.
C      It is not altered by the subroutine.
C A is an array that holds the values of the entries.
C IRN  is an integer array that must be set to the row indices of the
C      entries. It is not altered by the subroutine.
C ICN  is an integer array that must be set to the column indices of the
C      entries. It is not altered by the subroutine.
C R is an array that need not be set on entry. On return, it holds the
C      logarithms of the row scaling factors.
C C is an array that need not be set on entry. On return, it holds the
C      logarithms of the column scaling factors.
C W is a workarray.
C      W(1:M)  holds row non-zero counts (diagonal matrix M).
C      W(M+1:M+N) holds column non-zero counts (diagonal matrix N).
C      W(M+N+J) holds the logarithm of the column I scaling
C         factor during the iteration, J=1,2,...,N.
C      W(M+N*2+J) holds the 2-iteration change in the logarithm
C         of the column J scaling factor, J=1,2,...,N.
C      W(M+N*3+I) is used to save the average logarithm of
C          the entries of row I, I=1,2,...,M.
C LP must be set to the unit number for messages.
C      It is not altered by the subroutine.
C IFAIL need not be set by the user. On return it has one of the
C     following values:
C     0 successful entry.
C     -1 M < 1 or N < 1.
C     -2 NE < 1.

      INTRINSIC LOG,ABS,MIN

C Constants
      INTEGER MAXIT
      PARAMETER (MAXIT=100)
      DOUBLE PRECISION ONE,SMIN,ZERO
      PARAMETER (ONE=1D0,SMIN=0.1,ZERO=0D0)
C MAXIT is the maximal permitted number of iterations.
C SMIN is used in a convergence test on (residual norm)**2

C Local variables
      INTEGER I,I1,I2,I3,I4,I5,ITER,J,K
      DOUBLE PRECISION E,E1,EM,Q,Q1,QM,S,S1,SM,U,V

C Check M, N and NE.
      IFAIL = 0
      IF (M.LT.1 .OR. N.LT.1) THEN
         IFAIL = -1
         GO TO 220
      ELSE IF (NE.LE.0) THEN
         IFAIL = -2
         GO TO 220
      END IF

C     Partition W
      I1 = 0
      I2 = M
      I3 = M + N
      I4 = M + N*2
      I5 = M + N*3

C     Initialise for accumulation of sums and products.
      DO 10 I = 1,M
         R(I) = ZERO
         W(I1+I) = ZERO
   10 CONTINUE
      DO 20 J = 1,N
         C(J) = ZERO
         W(I2+J) = ZERO
         W(I3+J) = ZERO
         W(I4+J) = ZERO
   20 CONTINUE

C     Count non-zeros in the rows, and compute rhs vectors.
      DO 30 K = 1,NE
         U = ABS(A(K))
         IF (U.EQ.ZERO) GO TO 30
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 30
         U = LOG(U)
         W(I1+I) = W(I1+I) + ONE
         W(I2+J) = W(I2+J) + ONE
         R(I) = R(I) + U
         W(I3+J) = W(I3+J) + U
   30 CONTINUE
C
C     Divide rhs by diag matrices.
      DO 40 I = 1,M
         IF (W(I1+I).EQ.ZERO) W(I1+I) = ONE
         R(I) = R(I)/W(I1+I)
C     Save R(I) for use at end.
         W(I5+I) = R(I)
   40 CONTINUE
      DO 50 J = 1,N
         IF (W(I2+J).EQ.ZERO) W(I2+J) = ONE
         W(I3+J) = W(I3+J)/W(I2+J)
   50 CONTINUE
      SM = SMIN*NE

C     Sweep to compute initial residual vector
      DO 60 K = 1,NE
         IF (A(K).EQ.ZERO) GO TO 60
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 60
         R(I) = R(I) - W(I3+J)/W(I1+I)
   60 CONTINUE
C
C     Initialise iteration
      E = ZERO
      Q = ONE
      S = ZERO
      DO 70 I = 1,M
         S = S + W(I1+I)*R(I)**2
   70 CONTINUE
      IF (S.LE.SM) GO TO 160

C     Iteration loop
      DO 150 ITER = 1,MAXIT
C    Sweep through matrix to update residual vector
         DO 80 K = 1,NE
            IF (A(K).EQ.ZERO) GO TO 80
            J = ICN(K)
            I = IRN(K)
            IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 80
            C(J) = C(J) + R(I)
   80    CONTINUE
         S1 = S
         S = ZERO
         DO 90 J = 1,N
            V = -C(J)/Q
            C(J) = V/W(I2+J)
            S = S + V*C(J)
   90    CONTINUE
         E1 = E
         E = Q*S/S1
         Q = ONE - E
C      write(*,'(a,i3,a,f12.4)')' Iteration',ITER,' S =',S
         IF (S.LE.SM) E = ZERO
C     Update residual.
         DO 100 I = 1,M
            R(I) = R(I)*E*W(I1+I)
  100    CONTINUE
         IF (S.LE.SM) GO TO 180
         EM = E*E1
C    Sweep through matrix to update residual vector
         DO 110 K = 1,NE
            IF (A(K).EQ.ZERO) GO TO 110
            I = IRN(K)
            J = ICN(K)
            IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 110
            R(I) = R(I) + C(J)
  110    CONTINUE
         S1 = S
         S = ZERO
         DO 120 I = 1,M
            V = -R(I)/Q
            R(I) = V/W(I1+I)
            S = S + V*R(I)
  120    CONTINUE
         E1 = E
         E = Q*S/S1
         Q1 = Q
         Q = ONE - E
C     Special fixup for last iteration.
         IF (S.LE.SM) Q = ONE
C     Update col. scaling powers
         QM = Q*Q1
         DO 130 J = 1,N
            W(I4+J) = (EM*W(I4+J)+C(J))/QM
            W(I3+J) = W(I3+J) + W(I4+J)
  130    CONTINUE
C      write(*,'(a,i3,a,f12.4)')' Iteration',ITER,' S =',S
         IF (S.LE.SM) GO TO 160
C     UPDATE RESIDUAL.
         DO 140 J = 1,N
            C(J) = C(J)*E*W(I2+J)
  140    CONTINUE
  150 CONTINUE
  160 DO 170 I = 1,M
         R(I) = R(I)*W(I1+I)
  170 CONTINUE
C
C     Sweep through matrix to prepare to get row scaling powers
  180 DO 190 K = 1,NE
         IF (A(K).EQ.ZERO) GO TO 190
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 190
         R(I) = R(I) + W(I3+J)
  190 CONTINUE
C
C     Final conversion to output values.
      DO 200 I = 1,M
         R(I) = R(I)/W(I1+I) - W(I5+I)
  200 CONTINUE
      DO 210 J = 1,N
         C(J) = -W(I3+J)
  210 CONTINUE
      RETURN

C Error returns
  220 IF (LP.GT.0) WRITE (LP,'(/A/A,I3)')
     +    ' **** Error return from MC29AD ****',' IFAIL =',IFAIL

      END
* COPYRIGHT (c) 1988 AEA Technology and
* Council for the Central Laboratory of the Research Councils
C Original date 14 June 2001
C  June 2001: threadsafe version of MC41
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC71AD(N,KASE,X,EST,W,IW,KEEP)
C
C      MC71A/AD ESTIMATES THE 1-NORM OF A SQUARE MATRIX A.
C      REVERSE COMMUNICATION IS USED FOR EVALUATING
C      MATRIX-VECTOR PRODUCTS.
C
C
C         N       INTEGER
C                 THE ORDER OF THE MATRIX.  N .GE. 1.
C
C         KASE    INTEGER
C                 SET INITIALLY TO ZERO . IF N .LE. 0 SET TO -1
C                 ON INTERMEDIATE RETURN
C                 = 1 OR 2.
C                ON FINAL RETURN
C                 =  0  ,IF SUCCESS
C                 = -1  ,IF N .LE.0
C
C         X       DOUBLE PRECISION ARRAY OF DIMENSION (N)
C                 IF 1-NORM IS REQUIRED
C                 MUST BE OVERWRITTEN BY
C
C                      A*X,             IF KASE=1,
C                      TRANSPOSE(A)*X,  IF KASE=2,
C
C                 AND MC71 MUST BE RE-CALLED, WITH ALL THE OTHER
C                 PARAMETERS UNCHANGED.
C                 IF INFINITY-NORM IS REQUIRED
C                 MUST BE OVERWRITTEN BY
C
C                      TRANSPOSE(A)*X,  IF KASE=1,
C                      A*X,             IF KASE=2,
C
C                 AND MC71 MUST BE RE-CALLED, WITH ALL THE OTHER
C                 PARAMETERS UNCHANGED.
C
C         EST     DOUBLE PRECISION
C                 CONTAINS AN ESTIMATE (A LOWER BOUND) FOR NORM(A).
C
C         W       DOUBLE PRECISION ARRAY OF DIMENSION (N)
C                 = A*V,   WHERE  EST = NORM(W)/NORM(V)
C                          (V  IS NOT RETURNED).
C         IW      INTEGER(N) USED AS WORKSPACE.
C
C         KEEP    INTEGER ARRAY LENGTH 5 USED TO PRESERVE PRIVATE
C                 DATA, JUMP, ITER, J AND JLAST BETWEEN CALLS,
C                 KEEP(5) IS SPARE.
C
C      REFERENCE
C      N.J. HIGHAM (1987) FORTRAN CODES FOR ESTIMATING
C      THE ONE-NORM OF A
C      REAL OR COMPLEX MATRIX, WITH APPLICATIONS
C      TO CONDITION  ESTIMATION, NUMERICAL ANALYSIS REPORT NO. 135,
C      UNIVERSITY OF MANCHESTER, MANCHESTER M13 9PL, ENGLAND.
C
C      SUBROUTINES AND FUNCTIONS
C
C
C      INTERNAL VARIABLES
C
C     .. Parameters ..
      INTEGER ITMAX
      PARAMETER (ITMAX=5)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION EST
      INTEGER KASE,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION W(*),X(*)
      INTEGER IW(*),KEEP(5)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALTSGN,TEMP
      INTEGER I,ITER,J,JLAST,JUMP
C     ..
C     .. External Functions ..
      INTEGER IDAMAX
      EXTERNAL IDAMAX
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SIGN,NINT,DBLE
C     ..
C     .. Executable Statements ..
C
      IF (N.LE.0) THEN
        KASE = -1
        RETURN

      END IF

      IF (KASE.EQ.0) THEN
        DO 10 I = 1,N
          X(I) = ONE/DBLE(N)
   10   CONTINUE
        KASE = 1
        JUMP = 1
        KEEP(1) = JUMP
        KEEP(2) = 0
        KEEP(3) = 0
        KEEP(4) = 0
        RETURN

      END IF
C
      JUMP  = KEEP(1)
      ITER  = KEEP(2)
      J     = KEEP(3)
      JLAST = KEEP(4)
C
      GO TO (100,200,300,400,500) JUMP
C
C      ................ ENTRY   (JUMP = 1)
C
  100 CONTINUE
      IF (N.EQ.1) THEN
        W(1) = X(1)
        EST = ABS(W(1))
C         ... QUIT
        GO TO 510

      END IF
C
      DO 110 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  110 CONTINUE
      KASE = 2
      JUMP = 2
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 2)
C
  200 CONTINUE
      J = IDAMAX(N,X,1)
      ITER = 2
C
C      MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
C
  220 CONTINUE
      DO 230 I = 1,N
        X(I) = ZERO
  230 CONTINUE
      X(J) = ONE
      KASE = 1
      JUMP = 3
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 3)
C
  300 CONTINUE
C
C      COPY X INTO W
C
      DO 310 I = 1,N
        W(I) = X(I)
  310 CONTINUE
      DO 320 I = 1,N
        IF (NINT(SIGN(ONE,X(I))).NE.IW(I)) GO TO 330
  320 CONTINUE
C
C      REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 410
C
  330 CONTINUE
      DO 340 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  340 CONTINUE
      KASE = 2
      JUMP = 4
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 4)
C
  400 CONTINUE
      JLAST = J
      J = IDAMAX(N,X,1)
      IF ((ABS(X(JLAST)).NE.ABS(X(J))) .AND. (ITER.LT.ITMAX)) THEN
        ITER = ITER + 1
        GO TO 220

      END IF
C
C      ITERATION COMPLETE.  FINAL STAGE.
C
  410 CONTINUE
      EST = ZERO
      DO 420 I = 1,N
        EST = EST + ABS(W(I))
  420 CONTINUE
C
      ALTSGN = ONE
      DO 430 I = 1,N
        X(I) = ALTSGN* (ONE+DBLE(I-1)/DBLE(N-1))
        ALTSGN = -ALTSGN
  430 CONTINUE
      KASE = 1
      JUMP = 5
      GO TO 1010
C
C      ................ ENTRY   (JUMP = 5)
C
  500 CONTINUE
      TEMP = ZERO
      DO 520 I = 1,N
        TEMP = TEMP + ABS(X(I))
  520 CONTINUE
      TEMP = 2.0*TEMP/DBLE(3*N)
      IF (TEMP.GT.EST) THEN
C
C      COPY X INTO W
C
        DO 530 I = 1,N
          W(I) = X(I)
  530   CONTINUE
        EST = TEMP
      END IF
C
  510 KASE = 0
C
 1010 CONTINUE
      KEEP(1) = JUMP
      KEEP(2) = ITER
      KEEP(3) = J
      KEEP(4) = JLAST
      RETURN
C
      END

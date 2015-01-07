#if !defined (PRG_DALTON)
C FILE    : printpkg.F
C
C...   Dalton, Release DALTON2013, pdpack/printpkg.F
C...
C...   These routines are in the public domain and can be
C...   used freely in other programs.
C...
C
C  /* Deck output */
      SUBROUTINE OUTPUT (AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,
     *                   NCTL,LUPRI)
C.......................................................................
C Revised 15-Dec-1983 by Hans Jorgen Aa. Jensen.
C         16-Jun-1986 hjaaj ( removed Hollerith )
C
C OUTPUT PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
C AND COLUMNS.  THE INPUT IS AS FOLLOWS;
C
C        AMATRX(',').........MATRIX TO BE OUTPUT
C
C        ROWLOW..............ROW NUMBER AT WHICH OUTPUT IS TO BEGIN
C
C        ROWHI...............ROW NUMBER AT WHICH OUTPUT IS TO END
C
C        COLLOW..............COLUMN NUMBER AT WHICH OUTPUT IS TO BEGIN
C
C        COLHI...............COLUMN NUMBER AT WHICH OUTPUT IS TO END
C
C        ROWDIM..............ROW DIMENSION OF AMATRX(',')
C
C        COLDIM..............COLUMN DIMENSION OF AMATRX(',')
C
C        NCTL................CARRIAGE CONTROL FLAG; 1 FOR SINGLE SPACE
C                                                   2 FOR DOUBLE SPACE
C                                                   3 FOR TRIPLE SPACE
C                            hjaaj: negative for 132 col width
C
C THE PARAMETERS THAT FOLLOW MATRIX ARE ALL OF TYPE INTEGER.  THE
C PROGRAM IS SET UP TO HANDLE 5 COLUMNS/PAGE WITH A 1P,5D24.15 FORMAT
C FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED,
C CHANGE FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER
C OF COLUMNS.
C
C AUTHOR;  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
C          FLORIDA, GAINESVILLE
C REVISED; FEBRUARY 26, 1971
C
C.......................................................................
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER   ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,BEGIN,KCOL
      DIMENSION AMATRX(ROWDIM,COLDIM)
      CHARACTER*1 ASA(3), BLANK, CTL
      CHARACTER   PFMT*20, COLUMN*8
      LOGICAL     IS_NAN
      PARAMETER (ZERO=0.D00, KCOLP=5, KCOLN=8)
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA COLUMN/'Column  '/, BLANK/' '/, ASA/' ', '0', '-'/
C
      IF (ROWHI.LT.ROWLOW) GO TO 3
      IF (COLHI.LT.COLLOW) GO TO 3
C
      AMAX = ZERO
      N_NAN = 0
      DO 10 J = COLLOW,COLHI
         DO 10 I = ROWLOW,ROWHI
            IF ( IS_NAN(AMATRX(I,J),AMATRX(I,J)) ) THEN
               N_NAN = N_NAN + 1
            ELSE
               AMAX = MAX( AMAX, ABS(AMATRX(I,J)) )
            END IF
   10 CONTINUE
      IF (N_NAN .GT. 0) WRITE (LUPRI,'(/T6,A,I10,A)')
     &   'WARNING: matrix contains',N_NAN,' NaN.'
      IF (AMAX .EQ. ZERO) THEN
         WRITE (LUPRI,'(/T6,A)') 'Zero matrix.'
         GO TO 3
      END IF
      IF (FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX) THEN
C        use F output format
         PFMT = '(A1,I7,2X,8F15.8)'
         thrpri = 0.5D-8
      ELSE
C        use 1PD output format
         PFMT = '(A1,I7,2X,1P,8D15.6)'
         thrpri = 1.0D-8*AMAX
      END IF
C
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE.3).AND.(MCTL.GT.0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
C
      LAST = MIN(COLHI,COLLOW+KCOL-1)
      DO 2 BEGIN = COLLOW,COLHI,KCOL
         WRITE (LUPRI,1000) (COLUMN,I,I = BEGIN,LAST)
         DO 1 K = ROWLOW,ROWHI
            DO 4 I = BEGIN,LAST
               GO TO 5
               IF (abs(AMATRX(K,I)).gt.thrpri) GO TO 5
    4       CONTINUE
         GO TO 1
    5       WRITE (LUPRI,PFMT) CTL,K,(AMATRX(K,I), I = BEGIN,LAST)
    1    CONTINUE
    2 LAST = MIN(LAST+KCOL,COLHI)
    3 WRITE(LUPRI,'(A)') '    ==== End of matrix output ===='
      RETURN
 1000 FORMAT (/10X,8(5X,A6,I4))
      END
C === Deck is_nan ===============
      LOGICAL FUNCTION IS_NAN(XA,XB)
C
C     May 2010, Hans Joergen Aa. Jensen
C     Purpose: IS_NAN(X,X) is true iff X is NAN
C
      REAL*8 XA, XB
      IS_NAN = XA .NE. XB
      RETURN
      END
#endif

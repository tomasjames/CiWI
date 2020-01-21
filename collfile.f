      SUBROUTINE ORDER(NG,NTD,Y0,RAD,TIME)
C
C------------------------------------------------------------------------------
C  Subroutine to list out the formation and destruction reactions for each 
C  species in order of importance.
C    IOTYP=0: Selected species (11) only
C    IOTYP=1: All gas-phase species
C    IOTYP=2: All species
C------------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER(NSP=100,NCONS=2,NRM=1500,NSEL=16)
C
      CHARACTER*8 R1(NRM),R2(NRM),P1(NRM),P2(NRM),P3(NRM),P4(NRM),
     *            SPECI(NSP),SPSEL(NSEL),SPEC
      REAL*8 G(NRM),R(NRM),X(NSP),Y0(NSP),XS(NCONS),TOT(NCONS),FLAG(10)
      INTEGER INDXJ(NRM),LAB(10) 
C
      COMMON/BLK1/SPECI
      COMMON/BLK3/XS,TOT,R,D
      COMMON/ORDER1/R1,R2,P1,P2,P3,P4
      COMMON/ORDER2/INDXJ,IOTYP
C
      DATA  SPSEL/'N2H+','N2','NH3','ELECTR','CO','C','C+','HCO+',
     *      'CS','O','H2O','GNH3','GN2','GCO','GH2O','GO2'/
C      DATA (SPSEL(I),I=1,33)/'H2+','H3+','CH','CH+','CH2','CH2+','CH3+',
C     *     'CN','N2','N2H+','N+','NH','NH+','NH2','NH2+','NH3+','NH4+',
C     *     'NH3','H2NC+','NO','HCO+','H3O+','OH','O2','CO','C2H+','H+',
C     *     'HE+','C+','S+','C2S','HC2S+','ELECTR'/
      DATA LOR/2/
C
      WRITE(LOR,4) RAD,TIME
 4    FORMAT(1X,'Analysis at ',1PE8.2,' pc. and ',1PE8.2,' years',/)
C
      DO 112 K1=1,NRM
         G(K1) = R(K1)
 112  CONTINUE
      DO 113 K1=1,NCONS
         X(K1) = XS(K1)*D
 113  CONTINUE
      DO 114 K1=1,NTD
         X(NCONS+K1) = Y0(K1)*D
 114  CONTINUE
      NTOT=NTD+NCONS
C
      DO 1 J=1,NRM
        IF(INDXJ(J).EQ.9999) GO TO 2
        DO 14 LI=1,NTOT
           IF(R1(J).EQ.SPECI(LI)) G(J)=G(J)*X(LI)
           IF(R2(J).EQ.SPECI(LI)) G(J)=G(J)*X(LI)
 14     CONTINUE
        IF(R1(J).EQ.'G') G(J)=G(J)*D
        IF(R2(J).EQ.'G') G(J)=G(J)*D
 1    CONTINUE
 2    L=J-1
C
C------------------------------------------------------------------------------
C
      DO 90 I=1,L-1
        DO 91 K=I+1,L
          J1=INDXJ(I)
          J2=INDXJ(K)
          IF(G(J1).LT.G(J2)) THEN
             INDXJ(I)=J2
             INDXJ(K)=J1
          END IF
 91     CONTINUE
 90   CONTINUE
C
C----Select species
C
      IF(IOTYP.EQ.0) LMAX=NSEL
      IF(IOTYP.EQ.1) LMAX=NCONS+NG
      IF(IOTYP.EQ.2) LMAX=NTOT
C
      DO 13 K=1,LMAX
        IF(IOTYP.EQ.0) THEN
           SPEC=SPSEL(K)
           DO 18 ISEL=1,NTOT
 18            IF(SPECI(ISEL).EQ.SPEC) GO TO 19
 19        CONTINUE
        ELSE 
           SPEC=SPECI(K)
           ISEL=K
        END IF
        FRATE = 0.0
        DRATE = 0.0
        K2=1
        DO 16 K1=1,L
          M=INDXJ(K1)
          IF((R1(M).EQ.SPEC).OR.(R2(M).EQ.SPEC).OR.(P1(M).EQ.SPEC).OR.
     *      (P2(M).EQ.SPEC).OR.(P3(M).EQ.SPEC).OR.(P4(M).EQ.SPEC)) THEN
          FLAG(K2) = 1.0
          IF((R1(M).EQ.SPEC).OR.(R2(M).EQ.SPEC)) FLAG(K2)= -1.0
          IF(FLAG(K2).EQ.1.0) FRATE  = FRATE+G(M)
          IF(FLAG(K2).EQ.-1.0) DRATE = DRATE-G(M)
          LAB(K2)=M
          K2=K2+1
        END IF
        IF(K2.GT.10) GO TO 17
 16   CONTINUE
 17   WRITE(LOR,9)
      WRITE(LOR,11) SPEC,FRATE,-DRATE,X(ISEL)
      WRITE(LOR,10)
      DO 12 I=1,K2-1
      J=LAB(I)
      IR=0
      IF((FLAG(I).EQ.1.0).AND.(FRATE.NE.0.0)) IR=NINT(G(J)*100.0/FRATE)
      IF((FLAG(I).EQ.-1.0).AND.(DRATE.NE.0.0)) IR=NINT(G(J)*100.0/DRATE)
      IF(IR.NE.0) WRITE(LOR,15) J,R1(J),R2(J),P1(J),P2(J),P3(J),P4(J),IR
 12   CONTINUE
 13   CONTINUE
C
      WRITE(LOR,20)
 20   FORMAT(1X,78('-'),/)
 15   FORMAT(1X,I4,1X,5(1A8,1X),A8,2X,I4,'%')
 9    FORMAT(//,1X,5X,12('*'))
 10   FORMAT(1X,5X,12('*'),/)
 11   FORMAT(1X,5X,'*',1X,1A8,1X,'*',2X,'Frate=',1PE9.3,2X,'Drate=',
     *       1PE9.3,2X,'Abundance=',1PE9.3)
      RETURN
      END
C
C******************************************************************************
C******************************************************************************
C******************************************************************************
C
      SUBROUTINE DLSODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1                  ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
      EXTERNAL F, JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      DOUBLE PRECISION Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
C***BEGIN PROLOGUE  DLSODE
C***PURPOSE  Livermore Solver for Ordinary Differential Equations.
C            DLSODE solves the initial-value problem for stiff or
C            nonstiff systems of first-order ODE's,
C               dy/dt = f(t,y),   or, in component form,
C               dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(N)),  i=1,...,N.
C This version: 12/11/2003.
C

C  Declare externals.
      EXTERNAL DPREPJ, DSOLSY
      DOUBLE PRECISION DUMACH, DVNORM
C
C  Declare all other variables.
      INTEGER INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, I1, I2, IFLAG, IMXER, KGO, LF0,
     1   LENIW, LENRW, LENWM, ML, MORD, MU, MXHNL0, MXSTP0
      DOUBLE PRECISION ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI,
     1   TCRIT, TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0
      DIMENSION MORD(2)
      LOGICAL IHIT
      CHARACTER*80 MSG
      SAVE MORD, MXSTP0, MXHNL0
C-----------------------------------------------------------------------
      COMMON /DLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C
      DATA  MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
C-----------------------------------------------------------------------
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .EQ. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
 20   IF (NEQ(1) .LE. 0) GO TO 604
      IF (ISTATE .EQ. 1) GO TO 25
      IF (NEQ(1) .GT. N) GO TO 605
 25   N = NEQ(1)
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      METH = MF/10
      MITER = MF - 10*METH
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LT. 0 .OR. MITER .GT. 5) GO TO 608
      IF (MITER .LE. 3) GO TO 30
      ML = IWORK(1)
      MU = IWORK(2)
      IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609
      IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610
 30   CONTINUE
C Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .EQ. 1) H0 = 0.0D0
      HMXI = 0.0D0
      HMIN = 0.0D0
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = MIN(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .NE. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. 0.0D0) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. 0.0D0) GO TO 615
      HMXI = 0.0D0
      IF (HMAX .GT. 0.0D0) HMXI = 1.0D0/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. 0.0D0) GO TO 616
C-----------------------------------------------------------------------
 60   LYH = 21
      IF (ISTATE .EQ. 1) NYH = N
      LWM = LYH + (MAXORD + 1)*NYH
      IF (MITER .EQ. 0) LENWM = 0
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) LENWM = N*N + 2
      IF (MITER .EQ. 3) LENWM = N + 2
      IF (MITER .GE. 4) LENWM = (2*ML + MU + 1)*N + 2
      LEWT = LWM + LENWM
      LSAVF = LEWT + N
      LACOR = LSAVF + N
      LENRW = LACOR + N - 1
      IWORK(17) = LENRW
      LIWM = 1
      LENIW = 20 + N
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) LENIW = 20
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
C Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 70 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. 0.0D0) GO TO 619
        IF (ATOLI .LT. 0.0D0) GO TO 620
 70     CONTINUE
      IF (ISTATE .EQ. 1) GO TO 100
C If ISTATE = 3, set flag to signal parameter changes to DSTODE. -------
      JSTART = -1
      IF (NQ .LE. MAXORD) GO TO 90
C MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
      DO 80 I = 1,N
 80     RWORK(I+LSAVF-1) = RWORK(I+LWM-1)
C Reload WM(1) = RWORK(LWM), since LWM may have changed. ---------------
 90   IF (MITER .GT. 0) RWORK(LWM) = SQRT(UROUND)
      IF (N .EQ. NYH) GO TO 200
C NEQ was reduced.  Zero part of YH to avoid undefined references. -----
      I1 = LYH + L*NYH
      I2 = LYH + (MAXORD + 1)*NYH - 1
      IF (I1 .GT. I2) GO TO 200
      DO 95 I = I1,I2
 95     RWORK(I) = 0.0D0
      GO TO 200
C-----------------------------------------------------------------------
 100  UROUND = DUMACH()
      TN = T
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0D0) GO TO 625
      IF (H0 .NE. 0.0D0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0D0)
     1   H0 = TCRIT - T
 110  JSTART = 0
      IF (MITER .GT. 0) RWORK(LWM) = SQRT(UROUND)
      NHNIL = 0
      NST = 0
      NJE = 0
      NSLAST = 0
      HU = 0.0D0
      NQU = 0
      CCMAX = 0.3D0
      MAXCOR = 3
      MSBP = 20
      MXNCF = 10
C Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH
      CALL F (NEQ, T, Y, RWORK(LF0))
      NFE = 1
C Load the initial value vector in YH. ---------------------------------
      DO 115 I = 1,N
 115    RWORK(I+LYH-1) = Y(I)
C Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      NQ = 1
      H = 1.0D0
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 120 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 621
 120    RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
C-----------------------------------------------------------------------
      IF (H0 .NE. 0.0D0) GO TO 180
      TDIST = ABS(TOUT - T)
      W0 = MAX(ABS(T),ABS(TOUT))
      IF (TDIST .LT. 2.0D0*UROUND*W0) GO TO 622
      TOL = RTOL(1)
      IF (ITOL .LE. 2) GO TO 140
      DO 130 I = 1,N
 130    TOL = MAX(TOL,RTOL(I))
 140  IF (TOL .GT. 0.0D0) GO TO 160
      ATOLI = ATOL(1)
      DO 150 I = 1,N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        AYI = ABS(Y(I))
        IF (AYI .NE. 0.0D0) TOL = MAX(TOL,ATOLI/AYI)
 150    CONTINUE
 160  TOL = MAX(TOL,100.0D0*UROUND)
      TOL = MIN(TOL,0.001D0)
      SUM = DVNORM (N, RWORK(LF0), RWORK(LEWT))
      SUM = 1.0D0/(TOL*W0*W0) + TOL*SUM**2
      H0 = 1.0D0/SQRT(SUM)
      H0 = MIN(H0,TDIST)
      H0 = SIGN(H0,TOUT-T)
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. 1.0D0) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      DO 190 I = 1,N
 190    RWORK(I+LF0-1) = H0*RWORK(I+LF0-1)
      GO TO 270
C-----------------------------------------------------------------------
 200  NSLAST = NST
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(1.0D0 + 100.0D0*UROUND)
      IF ((TP - TOUT)*H .GT. 0.0D0) GO TO 623
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. 0.0D0) GO TO 625
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 245
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND)
      IF (ISTATE .EQ. 2) JSTART = -2
C-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 510
 260    RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*DVNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. 1.0D0) GO TO 280
      TOLSF = TOLSF*2.0D0
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'DLSODE-  Warning..internal T (=R1) and H (=R2) are'
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWD (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      (H = step size). Solver will continue anyway'
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'DLSODE-  Above warning has been issued I1 times.  '
      CALL XERRWD (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      It will not be issued again for this problem'
      CALL XERRWD (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)
 290  CONTINUE
C-----------------------------------------------------------------------
      CALL DSTODE (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   RWORK(LSAVF), RWORK(LACOR), RWORK(LWM), IWORK(LIWM),
     2   F, JAC, DPREPJ, DSOLSY)
      KGO = 1 - KFLAG
      GO TO (300, 530, 540), KGO
C-----------------------------------------------------------------------
 300  INIT = 1
      GO TO (310, 400, 330, 340, 350), ITASK
C ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 310  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. 0.0D0) GO TO 400
      GO TO 250
C ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 345
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND)
      JSTART = -2
      GO TO 250
C ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
C-----------------------------------------------------------------------
 400  DO 410 I = 1,N
 410    Y(I) = RWORK(I+LYH-1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      RETURN
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'DLSODE-  At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWD (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWD (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0D0)
      ISTATE = -1
      GO TO 580
C EWT(I) .LE. 0.0 for some I (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'DLSODE-  At T (=R1), EWT(I1) has become R2 .LE. 0.'
      CALL XERRWD (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 580
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'DLSODE-  At T (=R1), too much accuracy requested  '
      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      for precision of machine..  see TOLSF (=R2) '
      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 580
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'DLSODE-  At T(=R1) and step size H(=R2), the error'
      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      test failed repeatedly or with ABS(H) = HMIN'
      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 560
C KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'DLSODE-  At T (=R1) and step size H (=R2), the    '
      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      or with ABS(H) = HMIN   '
      CALL XERRWD (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -5
C Compute IMXER if relevant. -------------------------------------------
 560  BIG = 0.0D0
      IMXER = 1
      DO 570 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 570
        BIG = SIZE
        IMXER = I
 570    CONTINUE
      IWORK(16) = IMXER
C Set Y vector, T, and optional outputs. -------------------------------
 580  DO 590 I = 1,N
 590    Y(I) = RWORK(I+LYH-1)
      T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      RETURN
C-----------------------------------------------------------------------
 601  MSG = 'DLSODE-  ISTATE (=I1) illegal '
      CALL XERRWD (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0D0, 0.0D0)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'DLSODE-  ITASK (=I1) illegal  '
      CALL XERRWD (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 603  MSG = 'DLSODE-  ISTATE .GT. 1 but DLSODE not initialized '
      CALL XERRWD (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 604  MSG = 'DLSODE-  NEQ (=I1) .LT. 1     '
      CALL XERRWD (MSG, 30, 4, 0, 1, NEQ(1), 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 605  MSG = 'DLSODE-  ISTATE = 3 and NEQ increased (I1 to I2)  '
      CALL XERRWD (MSG, 50, 5, 0, 2, N, NEQ(1), 0, 0.0D0, 0.0D0)
      GO TO 700
 606  MSG = 'DLSODE-  ITOL (=I1) illegal   '
      CALL XERRWD (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 607  MSG = 'DLSODE-  IOPT (=I1) illegal   '
      CALL XERRWD (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 608  MSG = 'DLSODE-  MF (=I1) illegal     '
      CALL XERRWD (MSG, 30, 8, 0, 1, MF, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 609  MSG = 'DLSODE-  ML (=I1) illegal.. .LT.0 or .GE.NEQ (=I2)'
      CALL XERRWD (MSG, 50, 9, 0, 2, ML, NEQ(1), 0, 0.0D0, 0.0D0)
      GO TO 700
 610  MSG = 'DLSODE-  MU (=I1) illegal.. .LT.0 or .GE.NEQ (=I2)'
      CALL XERRWD (MSG, 50, 10, 0, 2, MU, NEQ(1), 0, 0.0D0, 0.0D0)
      GO TO 700
 611  MSG = 'DLSODE-  MAXORD (=I1) .LT. 0  '
      CALL XERRWD (MSG, 30, 11, 0, 1, MAXORD, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 612  MSG = 'DLSODE-  MXSTEP (=I1) .LT. 0  '
      CALL XERRWD (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 613  MSG = 'DLSODE-  MXHNIL (=I1) .LT. 0  '
      CALL XERRWD (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 614  MSG = 'DLSODE-  TOUT (=R1) behind T (=R2)      '
      CALL XERRWD (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)
      MSG = '      Integration direction is given by H0 (=R1)  '
      CALL XERRWD (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0D0)
      GO TO 700
 615  MSG = 'DLSODE-  HMAX (=R1) .LT. 0.0  '
      CALL XERRWD (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0D0)
      GO TO 700
 616  MSG = 'DLSODE-  HMIN (=R1) .LT. 0.0  '
      CALL XERRWD (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0D0)
      GO TO 700
 617  CONTINUE
      MSG='DLSODE-  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWD (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)
      GO TO 700
 618  CONTINUE
      MSG='DLSODE-  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWD (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0D0, 0.0D0)
      GO TO 700
 619  MSG = 'DLSODE-  RTOL(I1) is R1 .LT. 0.0        '
      CALL XERRWD (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0D0)
      GO TO 700
 620  MSG = 'DLSODE-  ATOL(I1) is R1 .LT. 0.0        '
      CALL XERRWD (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0D0)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'DLSODE-  EWT(I1) is R1 .LE. 0.0         '
      CALL XERRWD (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0D0)
      GO TO 700
 622  CONTINUE
      MSG='DLSODE-  TOUT (=R1) too close to T(=R2) to start integration'
      CALL XERRWD (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  CONTINUE
      MSG='DLSODE-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWD (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  CONTINUE
      MSG='DLSODE-  ITASK = 4 OR 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWD (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  CONTINUE
      MSG='DLSODE-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWD (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'DLSODE-  At start of problem, too much accuracy   '
      CALL XERRWD (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      requested for precision of machine..  See TOLSF (=R1) '
      CALL XERRWD (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0D0)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG = 'DLSODE-  Trouble in DINTDY.  ITASK = I1, TOUT = R1'
      CALL XERRWD (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0D0)
C
 700  ISTATE = -3
      RETURN
C
 800  MSG = 'DLSODE-  Run aborted.. apparent infinite loop     '
      CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0D0, 0.0D0)
      RETURN
C----------------------- END OF SUBROUTINE DLSODE ----------------------
      END
C
*DECK DUMACH
      DOUBLE PRECISION FUNCTION DUMACH ()
C***PURPOSE  Compute the unit roundoff of the machine.
      DOUBLE PRECISION U, COMP
      U = 1.0D0
 10   U = U*0.5D0
      CALL DUMSUM(1.0D0, U, COMP)
      IF (COMP .NE. 1.0D0) GO TO 10
      DUMACH = U*2.0D0
      RETURN
C----------------------- End of Function DUMACH ------------------------
      END
      SUBROUTINE DUMSUM(A,B,C)
C     Routine to force normal storing of A + B, for DUMACH.
      DOUBLE PRECISION A, B, C
      C = A + B
      RETURN
      END
*DECK DCFODE
      SUBROUTINE DCFODE (METH, ELCO, TESCO)
C***PURPOSE  Set ODE integrator coefficients.
      INTEGER METH
      INTEGER I, IB, NQ, NQM1, NQP1
      DOUBLE PRECISION ELCO, TESCO
      DOUBLE PRECISION AGAMQ, FNQ, FNQM1, PC, PINT, RAGQ,
     1   RQFAC, RQ1FAC, TSIGN, XPIN
      DIMENSION ELCO(13,12), TESCO(3,12)
      DIMENSION PC(12)
C
      GO TO (100, 200), METH
C
 100  ELCO(1,1) = 1.0D0
      ELCO(2,1) = 1.0D0
      TESCO(1,1) = 0.0D0
      TESCO(2,1) = 2.0D0
      TESCO(1,2) = 1.0D0
      TESCO(3,12) = 0.0D0
      PC(1) = 1.0D0
      RQFAC = 1.0D0
      DO 140 NQ = 2,12
C-----------------------------------------------------------------------
        RQ1FAC = RQFAC
        RQFAC = RQFAC/NQ
        NQM1 = NQ - 1
        FNQM1 = NQM1
        NQP1 = NQ + 1
C Form coefficients of p(x)*(x+nq-1). ----------------------------------
        PC(NQ) = 0.0D0
        DO 110 IB = 1,NQM1
          I = NQP1 - IB
 110      PC(I) = PC(I-1) + FNQM1*PC(I)
        PC(1) = FNQM1*PC(1)
C Compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
        PINT = PC(1)
        XPIN = PC(1)/2.0D0
        TSIGN = 1.0D0
        DO 120 I = 2,NQ
          TSIGN = -TSIGN
          PINT = PINT + TSIGN*PC(I)/I
 120      XPIN = XPIN + TSIGN*PC(I)/(I+1)
C Store coefficients in ELCO and TESCO. --------------------------------
        ELCO(1,NQ) = PINT*RQ1FAC
        ELCO(2,NQ) = 1.0D0
        DO 130 I = 2,NQ
 130      ELCO(I+1,NQ) = RQ1FAC*PC(I)/I
        AGAMQ = RQFAC*XPIN
        RAGQ = 1.0D0/AGAMQ
        TESCO(2,NQ) = RAGQ
        IF (NQ .LT. 12) TESCO(1,NQP1) = RAGQ*RQFAC/NQP1
        TESCO(3,NQM1) = RAGQ
 140    CONTINUE
      RETURN
C
 200  PC(1) = 1.0D0
      RQ1FAC = 1.0D0
      DO 230 NQ = 1,5
C-----------------------------------------------------------------------
        FNQ = NQ
        NQP1 = NQ + 1
C Form coefficients of p(x)*(x+nq). ------------------------------------
        PC(NQP1) = 0.0D0
        DO 210 IB = 1,NQ
          I = NQ + 2 - IB
 210      PC(I) = PC(I-1) + FNQ*PC(I)
        PC(1) = FNQ*PC(1)
C Store coefficients in ELCO and TESCO. --------------------------------
        DO 220 I = 1,NQP1
 220      ELCO(I,NQ) = PC(I)/PC(2)
        ELCO(2,NQ) = 1.0D0
        TESCO(1,NQ) = RQ1FAC
        TESCO(2,NQ) = NQP1/ELCO(1,NQ)
        TESCO(3,NQ) = (NQ+2)/ELCO(1,NQ)
        RQ1FAC = RQ1FAC/FNQ
 230    CONTINUE
      RETURN
C----------------------- END OF SUBROUTINE DCFODE ----------------------
      END
*DECK DINTDY
      SUBROUTINE DINTDY (T, K, YH, NYH, DKY, IFLAG)
C***PURPOSE  Interpolate solution derivatives.
      INTEGER K, NYH, IFLAG
      DOUBLE PRECISION T, YH, DKY
      DIMENSION YH(NYH,*), DKY(*)
      INTEGER IOWND, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      DOUBLE PRECISION ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      COMMON /DLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   IOWND(6), IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, IC, J, JB, JB2, JJ, JJ1, JP1
      DOUBLE PRECISION C, R, S, TP
      CHARACTER*80 MSG
C
      IFLAG = 0
      IF (K .LT. 0 .OR. K .GT. NQ) GO TO 80
      TP = TN - HU -  100.0D0*UROUND*(TN + HU)
      IF ((T-TP)*(T-TN) .GT. 0.0D0) GO TO 90
C
      S = (T - TN)/H
      IC = 1
      IF (K .EQ. 0) GO TO 15
      JJ1 = L - K
      DO 10 JJ = JJ1,NQ
 10     IC = IC*JJ
 15   C = IC
      DO 20 I = 1,N
 20     DKY(I) = C*YH(I,L)
      IF (K .EQ. NQ) GO TO 55
      JB2 = NQ - K
      DO 50 JB = 1,JB2
        J = NQ - JB
        JP1 = J + 1
        IC = 1
        IF (K .EQ. 0) GO TO 35
        JJ1 = JP1 - K
        DO 30 JJ = JJ1,J
 30       IC = IC*JJ
 35     C = IC
        DO 40 I = 1,N
 40       DKY(I) = C*YH(I,JP1) + S*DKY(I)
 50     CONTINUE
      IF (K .EQ. 0) RETURN
 55   R = H**(-K)
      DO 60 I = 1,N
 60     DKY(I) = R*DKY(I)
      RETURN
C
 80   MSG = 'DINTDY-  K (=I1) illegal      '
      CALL XERRWD (MSG, 30, 51, 0, 1, K, 0, 0, 0.0D0, 0.0D0)
      IFLAG = -1
      RETURN
 90   MSG = 'DINTDY-  T (=R1) illegal      '
      CALL XERRWD (MSG, 30, 52, 0, 0, 0, 0, 1, T, 0.0D0)
      MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
      CALL XERRWD (MSG, 60, 52, 0, 0, 0, 0, 2, TP, TN)
      IFLAG = -2
      RETURN
C----------------------- END OF SUBROUTINE DINTDY ----------------------
      END
*DECK DPREPJ
      SUBROUTINE DPREPJ (NEQ, Y, YH, NYH, EWT, FTEM, SAVF, WM, IWM,
     1   F, JAC)
C***PURPOSE  Compute and process Newton iteration matrix.
      EXTERNAL F, JAC
      INTEGER NEQ, NYH, IWM
      DOUBLE PRECISION Y, YH, EWT, FTEM, SAVF, WM
      DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), FTEM(*), SAVF(*),
     1   WM(*), IWM(*)
      INTEGER IOWND, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      DOUBLE PRECISION ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      COMMON /DLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   IOWND(6), IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, I1, I2, IER, II, J, J1, JJ, LENP,
     1   MBA, MBAND, MEB1, MEBAND, ML, ML3, MU, NP1
      DOUBLE PRECISION CON, DI, FAC, HL0, R, R0, SRUR, YI, YJ, YJJ,
     1   DVNORM
C
      NJE = NJE + 1
      IERPJ = 0
      JCUR = 1
      HL0 = H*EL0
      GO TO (100, 200, 300, 400, 500), MITER
C If MITER = 1, call JAC and multiply by scalar. -----------------------
 100  LENP = N*N
      DO 110 I = 1,LENP
 110    WM(I+2) = 0.0D0
      CALL JAC (NEQ, TN, Y, 0, 0, WM(3), N)
      CON = -HL0
      DO 120 I = 1,LENP
 120    WM(I+2) = WM(I+2)*CON
      GO TO 240
C If MITER = 2, make N calls to F to approximate J. --------------------
 200  FAC = DVNORM (N, SAVF, EWT)
      R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0
      SRUR = WM(1)
      J1 = 2
      DO 230 J = 1,N
        YJ = Y(J)
        R = MAX(SRUR*ABS(YJ),R0/EWT(J))
        Y(J) = Y(J) + R
        FAC = -HL0/R
        CALL F (NEQ, TN, Y, FTEM)
        DO 220 I = 1,N
 220      WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
        Y(J) = YJ
        J1 = J1 + N
 230    CONTINUE
      NFE = NFE + N
C Add identity matrix. -------------------------------------------------
 240  J = 3
      NP1 = N + 1
      DO 250 I = 1,N
        WM(J) = WM(J) + 1.0D0
 250    J = J + NP1
C Do LU decomposition on P. --------------------------------------------
      CALL DGEFA (WM(3), N, N, IWM(21), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
C If MITER = 3, construct a diagonal approximation to J and P. ---------
 300  WM(2) = HL0
      R = EL0*0.1D0
      DO 310 I = 1,N
 310    Y(I) = Y(I) + R*(H*SAVF(I) - YH(I,2))
      CALL F (NEQ, TN, Y, WM(3))
      NFE = NFE + 1
      DO 320 I = 1,N
        R0 = H*SAVF(I) - YH(I,2)
        DI = 0.1D0*R0 - H*(WM(I+2) - SAVF(I))
        WM(I+2) = 1.0D0
        IF (ABS(R0) .LT. UROUND/EWT(I)) GO TO 320
        IF (ABS(DI) .EQ. 0.0D0) GO TO 330
        WM(I+2) = 0.1D0*R0/DI
 320    CONTINUE
      RETURN
 330  IERPJ = 1
      RETURN
C If MITER = 4, call JAC and multiply by scalar. -----------------------
 400  ML = IWM(1)
      MU = IWM(2)
      ML3 = ML + 3
      MBAND = ML + MU + 1
      MEBAND = MBAND + ML
      LENP = MEBAND*N
      DO 410 I = 1,LENP
 410    WM(I+2) = 0.0D0
      CALL JAC (NEQ, TN, Y, ML, MU, WM(ML3), MEBAND)
      CON = -HL0
      DO 420 I = 1,LENP
 420    WM(I+2) = WM(I+2)*CON
      GO TO 570
C If MITER = 5, make MBAND calls to F to approximate J. ----------------
 500  ML = IWM(1)
      MU = IWM(2)
      MBAND = ML + MU + 1
      MBA = MIN(MBAND,N)
      MEBAND = MBAND + ML
      MEB1 = MEBAND - 1
      SRUR = WM(1)
      FAC = DVNORM (N, SAVF, EWT)
      R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0
      DO 560 J = 1,MBA
        DO 530 I = J,N,MBAND
          YI = Y(I)
          R = MAX(SRUR*ABS(YI),R0/EWT(I))
 530      Y(I) = Y(I) + R
        CALL F (NEQ, TN, Y, FTEM)
        DO 550 JJ = J,N,MBAND
          Y(JJ) = YH(JJ,1)
          YJJ = Y(JJ)
          R = MAX(SRUR*ABS(YJJ),R0/EWT(JJ))
          FAC = -HL0/R
          I1 = MAX(JJ-MU,1)
          I2 = MIN(JJ+ML,N)
          II = JJ*MEB1 - ML + 2
          DO 540 I = I1,I2
 540        WM(II+I) = (FTEM(I) - SAVF(I))*FAC
 550      CONTINUE
 560    CONTINUE
      NFE = NFE + MBA
C Add identity matrix. -------------------------------------------------
 570  II = MBAND + 2
      DO 580 I = 1,N
        WM(II) = WM(II) + 1.0D0
 580    II = II + MEBAND
C Do LU decomposition of P. --------------------------------------------
      CALL DGBFA (WM(3), MEBAND, N, ML, MU, IWM(21), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
C----------------------- END OF SUBROUTINE DPREPJ ----------------------
      END
*DECK DSOLSY
      SUBROUTINE DSOLSY (WM, IWM, X, TEM)
C***PURPOSE  ODEPACK linear system solver.
      INTEGER IWM
      DOUBLE PRECISION WM, X, TEM
      DIMENSION WM(*), IWM(*), X(*), TEM(*)
      INTEGER IOWND, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      DOUBLE PRECISION ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      COMMON /DLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   IOWND(6), IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, MEBAND, ML, MU
      DOUBLE PRECISION DI, HL0, PHL0, R
C
      IERSL = 0
      GO TO (100, 100, 300, 400, 400), MITER
 100  CALL DGESL (WM(3), N, N, IWM(21), X, 0)
      RETURN
C
 300  PHL0 = WM(2)
      HL0 = H*EL0
      WM(2) = HL0
      IF (HL0 .EQ. PHL0) GO TO 330
      R = HL0/PHL0
      DO 320 I = 1,N
        DI = 1.0D0 - R*(1.0D0 - 1.0D0/WM(I+2))
        IF (ABS(DI) .EQ. 0.0D0) GO TO 390
 320    WM(I+2) = 1.0D0/DI
 330  DO 340 I = 1,N
 340    X(I) = WM(I+2)*X(I)
      RETURN
 390  IERSL = 1
      RETURN
C
 400  ML = IWM(1)
      MU = IWM(2)
      MEBAND = 2*ML + MU + 1
      CALL DGBSL (WM(3), MEBAND, N, ML, MU, IWM(21), X, 0)
      RETURN
C----------------------- END OF SUBROUTINE DSOLSY ----------------------
      END
*DECK DSRCOM
      SUBROUTINE DSRCOM (RSAV, ISAV, JOB)
C***PURPOSE  Save/restore ODEPACK COMMON blocks.
      INTEGER ISAV, JOB
      INTEGER ILS
      INTEGER I, LENILS, LENRLS
      DOUBLE PRECISION RSAV,   RLS
      DIMENSION RSAV(*), ISAV(*)
      SAVE LENRLS, LENILS
      COMMON /DLS001/ RLS(218), ILS(37)
      DATA LENRLS/218/, LENILS/37/
C
      IF (JOB .EQ. 2) GO TO 100
C
      DO 10 I = 1,LENRLS
 10     RSAV(I) = RLS(I)
      DO 20 I = 1,LENILS
 20     ISAV(I) = ILS(I)
      RETURN
C
 100  CONTINUE
      DO 110 I = 1,LENRLS
 110     RLS(I) = RSAV(I)
      DO 120 I = 1,LENILS
 120     ILS(I) = ISAV(I)
      RETURN
C----------------------- END OF SUBROUTINE DSRCOM ----------------------
      END
*DECK DSTODE
      SUBROUTINE DSTODE (NEQ, Y, YH, NYH, YH1, EWT, SAVF, ACOR,
     1   WM, IWM, F, JAC, PJAC, SLVS)
C***PURPOSE  Performs one step of an ODEPACK integration.
      EXTERNAL F, JAC, PJAC, SLVS
      INTEGER NEQ, NYH, IWM
      DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, ACOR, WM
      DIMENSION NEQ(*), Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*),
     1   ACOR(*), WM(*), IWM(*)
      INTEGER IOWND, IALTH, IPUP, LMAX, MEO, NQNYH, NSLP,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, I1, IREDO, IRET, J, JB, M, NCF, NEWQ
      DOUBLE PRECISION CONIT, CRATE, EL, ELCO, HOLD, RMAX, TESCO,
     2   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION DCON, DDN, DEL, DELP, DSM, DUP, EXDN, EXSM, EXUP,
     1   R, RH, RHDN, RHSM, RHUP, TOLD, DVNORM
      COMMON /DLS001/ CONIT, CRATE, EL(13), ELCO(13,12),
     1   HOLD, RMAX, TESCO(3,12),
     2   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     3   IOWND(6), IALTH, IPUP, LMAX, MEO, NQNYH, NSLP,
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C
      KFLAG = 0
      TOLD = TN
      NCF = 0
      IERPJ = 0
      IERSL = 0
      JCUR = 0
      ICF = 0
      DELP = 0.0D0
      IF (JSTART .GT. 0) GO TO 200
      IF (JSTART .EQ. -1) GO TO 100
      IF (JSTART .EQ. -2) GO TO 160
C-----------------------------------------------------------------------
      LMAX = MAXORD + 1
      NQ = 1
      L = 2
      IALTH = 2
      RMAX = 10000.0D0
      RC = 0.0D0
      EL0 = 1.0D0
      CRATE = 0.7D0
      HOLD = H
      MEO = METH
      NSLP = 0
      IPUP = MITER
      IRET = 3
      GO TO 140
C-----------------------------------------------------------------------
 100  IPUP = MITER
      LMAX = MAXORD + 1
      IF (IALTH .EQ. 1) IALTH = 2
      IF (METH .EQ. MEO) GO TO 110
      CALL DCFODE (METH, ELCO, TESCO)
      MEO = METH
      IF (NQ .GT. MAXORD) GO TO 120
      IALTH = L
      IRET = 1
      GO TO 150
 110  IF (NQ .LE. MAXORD) GO TO 160
 120  NQ = MAXORD
      L = LMAX
      DO 125 I = 1,L
 125    EL(I) = ELCO(I,NQ)
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5D0/(NQ+2)
      DDN = DVNORM (N, SAVF, EWT)/TESCO(1,L)
      EXDN = 1.0D0/L
      RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
      RH = MIN(RHDN,1.0D0)
      IREDO = 3
      IF (H .EQ. HOLD) GO TO 170
      RH = MIN(RH,ABS(H/HOLD))
      H = HOLD
      GO TO 175
C-----------------------------------------------------------------------
 140  CALL DCFODE (METH, ELCO, TESCO)
 150  DO 155 I = 1,L
 155    EL(I) = ELCO(I,NQ)
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5D0/(NQ+2)
      GO TO (160, 170, 200), IRET
C-----------------------------------------------------------------------
 160  IF (H .EQ. HOLD) GO TO 200
      RH = H/HOLD
      H = HOLD
      IREDO = 3
      GO TO 175
 170  RH = MAX(RH,HMIN/ABS(H))
 175  RH = MIN(RH,RMAX)
      RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
      R = 1.0D0
      DO 180 J = 2,L
        R = R*RH
        DO 180 I = 1,N
 180      YH(I,J) = YH(I,J)*R
      H = H*RH
      RC = RC*RH
      IALTH = L
      IF (IREDO .EQ. 0) GO TO 690
C-----------------------------------------------------------------------
 200  IF (ABS(RC-1.0D0) .GT. CCMAX) IPUP = MITER
      IF (NST .GE. NSLP+MSBP) IPUP = MITER
      TN = TN + H
      I1 = NQNYH + 1
      DO 215 JB = 1,NQ
        I1 = I1 - NYH
Cdir$ ivdep
        DO 210 I = I1,NQNYH
 210      YH1(I) = YH1(I) + YH1(I+NYH)
 215    CONTINUE
C-----------------------------------------------------------------------
 220  M = 0
      DO 230 I = 1,N
 230    Y(I) = YH(I,1)
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      IF (IPUP .LE. 0) GO TO 250
C-----------------------------------------------------------------------
      CALL PJAC (NEQ, Y, YH, NYH, EWT, ACOR, SAVF, WM, IWM, F, JAC)
      IPUP = 0
      RC = 1.0D0
      NSLP = NST
      CRATE = 0.7D0
      IF (IERPJ .NE. 0) GO TO 430
 250  DO 260 I = 1,N
 260    ACOR(I) = 0.0D0
 270  IF (MITER .NE. 0) GO TO 350
C-----------------------------------------------------------------------
      DO 290 I = 1,N
        SAVF(I) = H*SAVF(I) - YH(I,2)
 290    Y(I) = SAVF(I) - ACOR(I)
      DEL = DVNORM (N, Y, EWT)
      DO 300 I = 1,N
        Y(I) = YH(I,1) + EL(1)*SAVF(I)
 300    ACOR(I) = SAVF(I)
      GO TO 400
C-----------------------------------------------------------------------
 350  DO 360 I = 1,N
 360    Y(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
      CALL SLVS (WM, IWM, Y, SAVF)
      IF (IERSL .LT. 0) GO TO 430
      IF (IERSL .GT. 0) GO TO 410
      DEL = DVNORM (N, Y, EWT)
      DO 380 I = 1,N
        ACOR(I) = ACOR(I) + Y(I)
 380    Y(I) = YH(I,1) + EL(1)*ACOR(I)
C-----------------------------------------------------------------------
 400  IF (M .NE. 0) CRATE = MAX(0.2D0*CRATE,DEL/DELP)
      DCON = DEL*MIN(1.0D0,1.5D0*CRATE)/(TESCO(2,NQ)*CONIT)
      IF (DCON .LE. 1.0D0) GO TO 450
      M = M + 1
      IF (M .EQ. MAXCOR) GO TO 410
      IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP) GO TO 410
      DELP = DEL
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      GO TO 270
C-----------------------------------------------------------------------
 410  IF (MITER .EQ. 0 .OR. JCUR .EQ. 1) GO TO 430
      ICF = 1
      IPUP = MITER
      GO TO 220
 430  ICF = 2
      NCF = NCF + 1
      RMAX = 2.0D0
      TN = TOLD
      I1 = NQNYH + 1
      DO 445 JB = 1,NQ
        I1 = I1 - NYH
Cdir$ ivdep
        DO 440 I = I1,NQNYH
 440      YH1(I) = YH1(I) - YH1(I+NYH)
 445    CONTINUE
      IF (IERPJ .LT. 0 .OR. IERSL .LT. 0) GO TO 680
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 670
      IF (NCF .EQ. MXNCF) GO TO 670
      RH = 0.25D0
      IPUP = MITER
      IREDO = 1
      GO TO 170
C-----------------------------------------------------------------------
 450  JCUR = 0
      IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
      IF (M .GT. 0) DSM = DVNORM (N, ACOR, EWT)/TESCO(2,NQ)
      IF (DSM .GT. 1.0D0) GO TO 500
C-----------------------------------------------------------------------
      KFLAG = 0
      IREDO = 0
      NST = NST + 1
      HU = H
      NQU = NQ
      DO 470 J = 1,L
        DO 470 I = 1,N
 470      YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
      IALTH = IALTH - 1
      IF (IALTH .EQ. 0) GO TO 520
      IF (IALTH .GT. 1) GO TO 700
      IF (L .EQ. LMAX) GO TO 700
      DO 490 I = 1,N
 490    YH(I,LMAX) = ACOR(I)
      GO TO 700
C-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1
      TN = TOLD
      I1 = NQNYH + 1
      DO 515 JB = 1,NQ
        I1 = I1 - NYH
Cdir$ ivdep
        DO 510 I = I1,NQNYH
 510      YH1(I) = YH1(I) - YH1(I+NYH)
 515    CONTINUE
      RMAX = 2.0D0
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 660
      IF (KFLAG .LE. -3) GO TO 640
      IREDO = 2
      RHUP = 0.0D0
      GO TO 540
C-----------------------------------------------------------------------
 520  RHUP = 0.0D0
      IF (L .EQ. LMAX) GO TO 540
      DO 530 I = 1,N
 530    SAVF(I) = ACOR(I) - YH(I,LMAX)
      DUP = DVNORM (N, SAVF, EWT)/TESCO(3,NQ)
      EXUP = 1.0D0/(L+1)
      RHUP = 1.0D0/(1.4D0*DUP**EXUP + 0.0000014D0)
 540  EXSM = 1.0D0/L
      RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      RHDN = 0.0D0
      IF (NQ .EQ. 1) GO TO 560
      DDN = DVNORM (N, YH(1,L), EWT)/TESCO(1,NQ)
      EXDN = 1.0D0/NQ
      RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
 560  IF (RHSM .GE. RHUP) GO TO 570
      IF (RHUP .GT. RHDN) GO TO 590
      GO TO 580
 570  IF (RHSM .LT. RHDN) GO TO 580
      NEWQ = NQ
      RH = RHSM
      GO TO 620
 580  NEWQ = NQ - 1
      RH = RHDN
      IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
      GO TO 620
 590  NEWQ = L
      RH = RHUP
      IF (RH .LT. 1.1D0) GO TO 610
      R = EL(L)/L
      DO 600 I = 1,N
 600    YH(I,NEWQ+1) = ACOR(I)*R
      GO TO 630
 610  IALTH = 3
      GO TO 700
 620  IF ((KFLAG .EQ. 0) .AND. (RH .LT. 1.1D0)) GO TO 610
      IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
C-----------------------------------------------------------------------
      IF (NEWQ .EQ. NQ) GO TO 170
 630  NQ = NEWQ
      L = NQ + 1
      IRET = 2
      GO TO 150
C-----------------------------------------------------------------------
 640  IF (KFLAG .EQ. -10) GO TO 660
      RH = 0.1D0
      RH = MAX(HMIN/ABS(H),RH)
      H = H*RH
      DO 645 I = 1,N
 645    Y(I) = YH(I,1)
      CALL F (NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      DO 650 I = 1,N
 650    YH(I,2) = H*SAVF(I)
      IPUP = MITER
      IALTH = 5
      IF (NQ .EQ. 1) GO TO 200
      NQ = 1
      L = 2
      IRET = 3
      GO TO 150
C-----------------------------------------------------------------------
 660  KFLAG = -1
      GO TO 720
 670  KFLAG = -2
      GO TO 720
 680  KFLAG = -3
      GO TO 720
 690  RMAX = 10.0D0
 700  R = 1.0D0/TESCO(2,NQU)
      DO 710 I = 1,N
 710    ACOR(I) = ACOR(I)*R
 720  HOLD = H
      JSTART = 1
      RETURN
C----------------------- END OF SUBROUTINE DSTODE ----------------------
      END
*DECK DEWSET
      SUBROUTINE DEWSET (N, ITOL, RTOL, ATOL, YCUR, EWT)
C***PURPOSE  Set error weight vector.
      INTEGER N, ITOL
      INTEGER I
      DOUBLE PRECISION RTOL, ATOL, YCUR, EWT
      DIMENSION RTOL(*), ATOL(*), YCUR(N), EWT(N)
C
      GO TO (10, 20, 30, 40), ITOL
 10   CONTINUE
      DO 15 I = 1,N
 15     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(1)
      RETURN
 20   CONTINUE
      DO 25 I = 1,N
 25     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(I)
      RETURN
 30   CONTINUE
      DO 35 I = 1,N
 35     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(1)
      RETURN
 40   CONTINUE
      DO 45 I = 1,N
 45     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I)
      RETURN
C----------------------- END OF SUBROUTINE DEWSET ----------------------
      END
*DECK DVNORM
      DOUBLE PRECISION FUNCTION DVNORM (N, V, W)
C***PURPOSE  Weighted root-mean-square vector norm.
      INTEGER N,   I
      DOUBLE PRECISION V, W,   SUM
      DIMENSION V(N), W(N)
C
      SUM = 0.0D0
      DO 10 I = 1,N
 10     SUM = SUM + (V(I)*W(I))**2
      DVNORM = SQRT(SUM/N)
      RETURN
C----------------------- END OF FUNCTION DVNORM ------------------------
      END
C
*DECK DGEFA
      SUBROUTINE DGEFA (A, LDA, N, IPVT, INFO)
C***PURPOSE  Factor a matrix using Gaussian elimination.
      INTEGER LDA,N,IPVT(*),INFO
      DOUBLE PRECISION A(LDA,*)
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
*DECK DGESL
      SUBROUTINE DGESL (A, LDA, N, IPVT, B, JOB)
C***PURPOSE  Solve the real system A*X=B or TRANS(A)*X=B using the
C            factors computed by DGECO or DGEFA.
      INTEGER LDA,N,IPVT(*),JOB
      DOUBLE PRECISION A(LDA,*),B(*)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,NM1
C
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
*DECK DGBFA
      SUBROUTINE DGBFA (ABD, LDA, N, ML, MU, IPVT, INFO)
C***PURPOSE  Factor a band matrix using Gaussian elimination.
      INTEGER LDA,N,ML,MU,IPVT(*),INFO
      DOUBLE PRECISION ABD(LDA,*)
C
      DOUBLE PRECISION T
      INTEGER I,IDAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
C
C
      M = ML + MU + 1
      INFO = 0
C
C     ZERO INITIAL FILL-IN COLUMNS
C
      J0 = MU + 2
      J1 = MIN(N,M) - 1
      IF (J1 .LT. J0) GO TO 30
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0D0
   10    CONTINUE
   20 CONTINUE
   30 CONTINUE
      JZ = J1
      JU = 0
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 130
      DO 120 K = 1, NM1
         KP1 = K + 1
C
C        ZERO NEXT FILL-IN COLUMN
C
         JZ = JZ + 1
         IF (JZ .GT. N) GO TO 50
         IF (ML .LT. 1) GO TO 50
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0D0
   40       CONTINUE
   50    CONTINUE
C
C        FIND L = PIVOT INDEX
C
         LM = MIN(ML,N-K)
         L = IDAMAX(LM+1,ABD(M,K),1) + M - 1
         IPVT(K) = L + K - M
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (ABD(L,K) .EQ. 0.0D0) GO TO 100
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. M) GO TO 60
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
   60       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/ABD(M,K)
            CALL DSCAL(LM,T,ABD(M+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            JU = MIN(MAX(JU,MU+IPVT(K)),N)
            MM = M
            IF (JU .LT. KP1) GO TO 90
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .EQ. MM) GO TO 70
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
   70          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
   80       CONTINUE
   90       CONTINUE
         GO TO 110
  100    CONTINUE
            INFO = K
  110    CONTINUE
  120 CONTINUE
  130 CONTINUE
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
*DECK DGBSL
      SUBROUTINE DGBSL (ABD, LDA, N, ML, MU, IPVT, B, JOB)
C***PURPOSE  Solve the real band system A*X=B or TRANS(A)*X=B using
C            the factors computed by DGBCO or DGBFA.
      INTEGER LDA,N,ML,MU,IPVT(*),JOB
      DOUBLE PRECISION ABD(LDA,*),B(*)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,LA,LB,LM,M,NM1
C
      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE L*Y = B
C
         IF (ML .EQ. 0) GO TO 30
         IF (NM1 .LT. 1) GO TO 30
            DO 20 K = 1, NM1
               LM = MIN(ML,N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .EQ. K) GO TO 10
                  B(L) = B(K)
                  B(K) = T
   10          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
   20       CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/ABD(M,K)
            LM = MIN(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            LM = MIN(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = DDOT(LM,ABD(LA,K),1,B(LB),1)
            B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (ML .EQ. 0) GO TO 90
         IF (NM1 .LT. 1) GO TO 90
            DO 80 KB = 1, NM1
               K = N - KB
               LM = MIN(ML,N-K)
               B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1)
               L = IPVT(K)
               IF (L .EQ. K) GO TO 70
                  T = B(L)
                  B(L) = B(K)
                  B(K) = T
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
*DECK DAXPY
      SUBROUTINE DAXPY (N, DA, DX, INCX, DY, INCY)
C***PURPOSE  Compute a constant times a vector plus a vector.
      DOUBLE PRECISION DX(*), DY(*), DA
C
      IF (N.LE.0 .OR. DA.EQ.0.0D0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C     Clean-up loop so remaining vector length is a multiple of 4.
C
   20 M = MOD(N,4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N .LT. 4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I+1) = DY(I+1) + DA*DX(I+1)
        DY(I+2) = DY(I+2) + DA*DX(I+2)
        DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DA*DX(I) + DY(I)
   70 CONTINUE
      RETURN
      END
*DECK DDOT
      DOUBLE PRECISION FUNCTION DDOT (N, DX, INCX, DY, INCY)
C***PURPOSE  Compute the inner product of two vectors.
      DOUBLE PRECISION DX(*), DY(*)
C
      DDOT = 0.0D0
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
      DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)*DY(I+2) +
     1              DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DDOT = DDOT + DX(I)*DY(I)
   70 CONTINUE
      RETURN
      END
*DECK DSCAL
      SUBROUTINE DSCAL (N, DA, DX, INCX)
C***PURPOSE  Multiply a vector by a constant.
      DOUBLE PRECISION DA, DX(*)
      INTEGER I, INCX, IX, M, MP1, N
C
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DX(IX) = DA*DX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GOTO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I+1) = DA*DX(I+1)
        DX(I+2) = DA*DX(I+2)
        DX(I+3) = DA*DX(I+3)
        DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END
*DECK IDAMAX
      INTEGER FUNCTION IDAMAX (N, DX, INCX)
C***PURPOSE  Find the smallest index of that component of a vector
C            having the maximum magnitude.

      DOUBLE PRECISION DX(*), DMAX, XMAG
      INTEGER I, INCX, IX, N
C
      IDAMAX = 0
      IF (N .LE. 0) RETURN
      IDAMAX = 1
      IF (N .EQ. 1) RETURN
C
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increments not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DMAX = ABS(DX(IX))
      IX = IX + INCX
      DO 10 I = 2,N
        XMAG = ABS(DX(IX))
        IF (XMAG .GT. DMAX) THEN
          IDAMAX = I
          DMAX = XMAG
        ENDIF
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increments equal to 1.
C
   20 DMAX = ABS(DX(1))
      DO 30 I = 2,N
        XMAG = ABS(DX(I))
        IF (XMAG .GT. DMAX) THEN
          IDAMAX = I
          DMAX = XMAG
        ENDIF
   30 CONTINUE
      RETURN
      END
*DECK XERRWD
      SUBROUTINE XERRWD (MSG, NMES, NERR, LEVEL, NI, I1, I2, NR, R1, R2)
C***PURPOSE  Write error message with values.

      DOUBLE PRECISION R1, R2
      INTEGER NMES, NERR, LEVEL, NI, I1, I2, NR
      CHARACTER*(*) MSG
C
C  Declare local variables.
C
      INTEGER LUNIT, IXSAV, MESFLG
C
C  Get logical unit number and message print flag.
C
      LUNIT = IXSAV (1, 0, .FALSE.)
      MESFLG = IXSAV (2, 0, .FALSE.)
      IF (MESFLG .EQ. 0) GO TO 100
C
C  Write the message.
C
      WRITE (LUNIT,10)  MSG
 10   FORMAT(1X,A)
      IF (NI .EQ. 1) WRITE (LUNIT, 20) I1
 20   FORMAT(6X,'In above message,  I1 =',I10)
      IF (NI .EQ. 2) WRITE (LUNIT, 30) I1,I2
 30   FORMAT(6X,'In above message,  I1 =',I10,3X,'I2 =',I10)
      IF (NR .EQ. 1) WRITE (LUNIT, 40) R1
 40   FORMAT(6X,'In above message,  R1 =',D21.13)
      IF (NR .EQ. 2) WRITE (LUNIT, 50) R1,R2
 50   FORMAT(6X,'In above,  R1 =',D21.13,3X,'R2 =',D21.13)
C
C  Abort the run if LEVEL = 2.
C
 100  IF (LEVEL .NE. 2) RETURN
      STOP
C----------------------- End of Subroutine XERRWD ----------------------
      END
*DECK XSETF
      SUBROUTINE XSETF (MFLAG)
C***PURPOSE  Reset the error print control flag.
      INTEGER MFLAG, JUNK, IXSAV
C
      IF (MFLAG .EQ. 0 .OR. MFLAG .EQ. 1) JUNK = IXSAV (2,MFLAG,.TRUE.)
      RETURN
C----------------------- End of Subroutine XSETF -----------------------
      END
*DECK XSETUN
      SUBROUTINE XSETUN (LUN)
C***PURPOSE  Reset the logical unit number for error messages.
      INTEGER LUN, JUNK, IXSAV
C
      IF (LUN .GT. 0) JUNK = IXSAV (1,LUN,.TRUE.)
      RETURN
C----------------------- End of Subroutine XSETUN ----------------------
      END
*DECK IXSAV
      INTEGER FUNCTION IXSAV (IPAR, IVALUE, ISET)
C***PURPOSE  Save and recall error message control parameters.

      LOGICAL ISET
      INTEGER IPAR, IVALUE
C-----------------------------------------------------------------------
      INTEGER IUMACH, LUNIT, MESFLG
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this routine.
C-----------------------------------------------------------------------
      SAVE LUNIT, MESFLG
      DATA LUNIT/-1/, MESFLG/1/
C
      IF (IPAR .EQ. 1) THEN
        IF (LUNIT .EQ. -1) LUNIT = IUMACH()
        IXSAV = LUNIT
        IF (ISET) LUNIT = IVALUE
        ENDIF
C
      IF (IPAR .EQ. 2) THEN
        IXSAV = MESFLG
        IF (ISET) MESFLG = IVALUE
        ENDIF
C
      RETURN
C----------------------- End of Function IXSAV -------------------------
      END
*DECK IUMACH
      INTEGER FUNCTION IUMACH()
C***PURPOSE  Provide standard output unit number.
      IUMACH = 6
C
      RETURN
C----------------------- End of Function IUMACH ------------------------
      END
C
C******************************************************************************
C******************************************************************************
C******************************************************************************
C
      SUBROUTINE DIFFUN(N,T,Y,YDOT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*8 K
C
      PARAMETER(NRM=1500,NSP=100,NCONS=2)
C
      DIMENSION K(NRM),Y(N),YDOT(N),YD(20),X(NCONS),TOTAL(NCONS)
      COMMON /BLK3/X,TOTAL,K,D
      COMMON /BLK9/FCO2,COCOV,FOX,PHOH,H2SHL,COSHL,FFRZ,FDES
C
      CALL ADJUST(N,T,Y,YDOT)
C
      X( 1) = TOTAL( 1)-0.5*(Y(1  )+Y(2  )+Y(3  )+2*Y(4  )+3*Y(
     *        5  )+Y(13 )+Y(14 )+2*Y(15 )+2*Y(16 )+3*Y(17 )+3*Y(18 )
     *        +4*Y(19 )+4*Y(20 )+5*Y(21 )+Y(24 )+Y(25 )+2*Y(26 )+2*Y
     *        (27 )+3*Y(28 )+3*Y(29 )+4*Y(30 )+Y(33 )+Y(38 )+Y(39 )+
     *        2*Y(40 )+2*Y(41 )+3*Y(42 )+Y(45 )+Y(46 )+2*Y(47 )+2*Y(
     *        48 )+Y(51 )+Y(54 )+Y(55 )+Y(56 )+Y(59 )+Y(60 )+Y(61 )+
     *        Y(61 )+2*Y(62 )+Y(63 )+2*Y(64 )+Y(67 )+Y(68 )+2*Y(69 )
     *        +2*Y(70 )+3*Y(71 )+Y(76 )+Y(77 )+Y(78 )+2*Y(79 )+4*Y(8
     *        0 )+3*Y(82 )+2*Y(85 )+2*Y(87 )+Y(89 )+Y(90 )+Y(91 )+2*
     *        Y(92 )+Y(94 ))                                        
      X( 2) = TOTAL( 2)+(Y(2  )-Y(3  )+Y(4  )+Y(5  )+Y(7  )+Y(9
     *          )-Y(10 )+Y(12 )+Y(14 )+Y(16 )+Y(18 )+Y(20 )+Y(21 )+Y
     *        (23 )+Y(25 )+Y(27 )+Y(29 )+Y(30 )+Y(32 )+Y(33 )+Y(35 )
     *        +Y(37 )+Y(39 )+Y(41 )+Y(42 )+Y(44 )+Y(46 )+Y(48 )+Y(50
     *         )+Y(51 )+Y(53 )+Y(55 )+Y(58 )+Y(60 )+Y(61 )+Y(62 )+Y(
     *        63 )+Y(64 )+Y(66 )+Y(68 )+Y(70 )+Y(71 )+Y(73 )+Y(75 )+
     *        Y(76 )+Y(78 )+Y(79 ))                                 
C
      YD( 1) = 
     *        +K(  10)*Y(2  )*Y(34 )*D  +K(  23)*Y(15 )*Y(36 )*D  +K
     *        (  23)*Y(15 )*Y(36 )*D  +K(  33)*Y(17 )*Y(34 )*D  +K( 
     *         64)*Y(11 )*Y(38 )*D  +K(  71)*Y(7  )*X(1 )*D  +K(  72
     *        )*Y(13 )*X(1 )*D  +K(  74)*Y(20 )*X(1 )*D  +K(  75)*Y(
     *        23 )*X(1 )*D  +K(  76)*Y(26 )*X(1 )*D  +K(  77)*Y(34 )
     *        *X(1 )*D  +K(  78)*Y(38 )*X(1 )*D  +K(  79)*Y(52 )*X(1
     *         )*D  +K(  80)*Y(70 )*X(1 )*D  +K(  82)*Y(8  )*Y(26 )*
     *        D  +K(  83)*Y(8  )*Y(26 )*D  +K(  86)*Y(26 )*Y(34 )*D 
     *         +K(  88)*Y(29 )*X(1 )*D  +K( 103)*Y(1  )*Y(56 )*D  +K
     *        ( 109)*Y(2  )*Y(3  )*D  +K( 109)*Y(2  )*Y(3  )*D  +K( 
     *        110)*Y(2  )*Y(10 )*D  +K( 111)*Y(2  )*Y(13 )*D  +K( 11
     *        3)*Y(2  )*Y(15 )*D  +K( 114)*Y(2  )*Y(17 )*D  +K( 116)
     *        *Y(2  )*Y(19 )*D  +K( 117)*Y(2  )*Y(24 )*D  +K( 118)*Y
     *        (2  )*Y(26 )*D  +K( 119)*Y(2  )*Y(28 )*D  +K( 120)*Y(2
     *          )*Y(36 )*D  +K( 121)*Y(2  )*Y(38 )*D  +K( 122)*Y(2  
     *        )*Y(40 )*D  +K( 125)*Y(2  )*Y(45 )*D  +K( 126)*Y(2  )*
     *        Y(47 )*D  +K( 128)*Y(2  )*Y(47 )*D  +K( 130)*Y(2  )*Y(
     *        54 )*D  +K( 132)*Y(2  )*Y(57 )*D                      
      YD( 2) = YD( 1)
     *        +K( 134)*Y(2  )*Y(65 )*D  +K( 136)*Y(2  )*Y(67 )*D  +K
     *        ( 138)*Y(2  )*Y(69 )*D  +K( 139)*Y(2  )*Y(69 )*D  +K( 
     *        140)*Y(2  )*Y(72 )*D  +K( 141)*Y(2  )*Y(74 )*D  +K( 14
     *        4)*Y(3  )*Y(4  )*D  +K( 146)*Y(3  )*Y(7  )*D  +K( 148)
     *        *Y(3  )*Y(9  )*D  +K( 154)*Y(3  )*Y(23 )*D  +K( 159)*Y
     *        (3  )*Y(35 )*D  +K( 161)*Y(3  )*Y(42 )*D  +K( 163)*Y(3
     *          )*Y(44 )*D  +K( 167)*Y(3  )*Y(66 )*D  +K( 168)*Y(4  
     *        )*Y(8  )*D  +K( 169)*Y(4  )*Y(11 )*D  +K( 171)*Y(4  )*
     *        Y(13 )*D  +K( 173)*Y(4  )*Y(15 )*D  +K( 175)*Y(4  )*Y(
     *        19 )*D  +K( 176)*Y(4  )*Y(19 )*D  +K( 178)*Y(4  )*Y(22
     *         )*D  +K( 179)*Y(4  )*Y(24 )*D  +K( 183)*Y(4  )*Y(31 )
     *        *D  +K( 184)*Y(4  )*Y(34 )*D  +K( 186)*Y(4  )*Y(38 )*D
     *          +K( 188)*Y(4  )*Y(40 )*D  +K( 192)*Y(4  )*Y(47 )*D  
     *        +K( 194)*Y(4  )*Y(49 )*D  +K( 195)*Y(4  )*Y(52 )*D  +K
     *        ( 198)*Y(4  )*Y(57 )*D  +K( 201)*Y(4  )*Y(69 )*D  +K( 
     *        203)*Y(4  )*X(1 )*D  +K( 215)*Y(5  )*Y(34 )*D  +K( 218
     *        )*Y(5  )*Y(43 )*D  +K( 235)*Y(7  )*Y(13 )*D  +K( 238)*
     *        Y(7  )*Y(15 )*D  +K( 240)*Y(7  )*Y(19 )*D             
      YD( 3) = YD( 2)
     *        +K( 243)*Y(7  )*Y(19 )*D  +K( 245)*Y(7  )*Y(24 )*D  +K
     *        ( 247)*Y(7  )*Y(26 )*D  +K( 249)*Y(7  )*Y(28 )*D  +K( 
     *        255)*Y(7  )*Y(38 )*D  +K( 257)*Y(7  )*Y(40 )*D  +K( 26
     *        0)*Y(7  )*Y(45 )*D  +K( 263)*Y(7  )*Y(47 )*D  +K( 274)
     *        *Y(7  )*Y(54 )*D  +K( 275)*Y(7  )*Y(54 )*D  +K( 277)*Y
     *        (7  )*Y(56 )*D  +K( 278)*Y(7  )*Y(56 )*D  +K( 281)*Y(7
     *          )*Y(59 )*D  +K( 283)*Y(7  )*Y(67 )*D  +K( 285)*Y(7  
     *        )*Y(69 )*D  +K( 292)*Y(7  )*Y(77 )*D  +K( 295)*Y(8  )*
     *        Y(24 )*D  +K( 304)*Y(8  )*Y(38 )*D  +K( 317)*Y(8  )*Y(
     *        67 )*D  +K( 318)*Y(8  )*Y(68 )*D  +K( 319)*Y(8  )*Y(70
     *         )*D  +K( 323)*Y(9  )*Y(24 )*D  +K( 324)*Y(9  )*Y(26 )
     *        *D  +K( 326)*Y(9  )*Y(28 )*D  +K( 331)*Y(9  )*Y(38 )*D
     *          +K( 332)*Y(9  )*Y(40 )*D  +K( 343)*Y(9  )*Y(67 )*D  
     *        +K( 344)*Y(9  )*Y(69 )*D  +K( 374)*Y(13 )*Y(22 )*D  +K
     *        ( 375)*Y(13 )*Y(23 )*D  +K( 383)*Y(13 )*Y(34 )*D  +K( 
     *        385)*Y(13 )*Y(35 )*D  +K( 408)*Y(13 )*Y(65 )*D  +K( 40
     *        9)*Y(13 )*Y(66 )*D  +K( 411)*Y(14 )*Y(22 )*D  +K( 417)
     *        *Y(14 )*Y(34 )*D  +K( 423)*Y(14 )*Y(40 )*D            
      YD( 4) = YD( 3)
     *        +K( 435)*Y(14 )*Y(65 )*D  +K( 448)*Y(15 )*Y(34 )*D  +K
     *        ( 448)*Y(15 )*Y(34 )*D  +K( 450)*Y(15 )*Y(34 )*D  +K( 
     *        454)*Y(15 )*Y(38 )*D  +K( 471)*Y(15 )*Y(65 )*D  +K( 47
     *        3)*Y(15 )*Y(66 )*D  +K( 475)*Y(16 )*Y(34 )*D  +K( 481)
     *        *Y(16 )*Y(65 )*D  +K( 483)*Y(16 )*Y(69 )*D  +K( 484)*Y
     *        (17 )*Y(34 )*D  +K( 489)*Y(17 )*Y(66 )*D  +K( 491)*Y(1
     *        8 )*Y(34 )*D  +K( 505)*Y(19 )*Y(32 )*D  +K( 515)*Y(19 
     *        )*Y(66 )*D  +K( 535)*Y(15 )*Y(22 )*D  +K( 536)*Y(15 )*
     *        Y(22 )*D  +K( 537)*Y(16 )*Y(22 )*D  +K( 538)*Y(17 )*Y(
     *        22 )*D  +K( 538)*Y(17 )*Y(22 )*D  +K( 540)*Y(18 )*Y(22
     *         )*D  +K( 541)*Y(22 )*Y(24 )*D  +K( 542)*Y(22 )*Y(25 )
     *        *D  +K( 543)*Y(22 )*Y(27 )*D  +K( 546)*Y(22 )*Y(38 )*D
     *          +K( 547)*Y(22 )*Y(39 )*D  +K( 549)*Y(22 )*Y(41 )*D  
     *        +K( 557)*Y(19 )*Y(23 )*D  +K( 558)*Y(19 )*Y(23 )*D  +K
     *        ( 559)*Y(19 )*Y(23 )*D  +K( 559)*Y(19 )*Y(23 )*D  +K( 
     *        561)*Y(23 )*Y(24 )*D  +K( 585)*Y(23 )*Y(69 )*D  +K( 59
     *        1)*Y(24 )*Y(24 )*D  +K( 591)*Y(24 )*Y(24 )*D  +K( 598)
     *        *Y(24 )*Y(34 )*D  +K( 599)*Y(24 )*Y(35 )*D            
      YD( 5) = YD( 4)
     *        +K( 602)*Y(24 )*Y(38 )*D  +K( 609)*Y(24 )*Y(57 )*D  +K
     *        ( 655)*Y(26 )*Y(57 )*D  +K( 719)*Y(32 )*Y(47 )*D  +K( 
     *        725)*Y(32 )*Y(69 )*D  +K( 733)*Y(27 )*Y(34 )*D  +K( 73
     *        9)*Y(34 )*Y(38 )*D  +K( 740)*Y(34 )*Y(39 )*D  +K( 742)
     *        *Y(34 )*Y(45 )*D  +K( 744)*Y(34 )*Y(47 )*D  +K( 761)*Y
     *        (35 )*Y(38 )*D  +K( 786)*Y(38 )*Y(46 )*D  +K( 826)*Y(4
     *        0 )*Y(53 )*D  +K( 848)*Y(42 )*Y(43 )*D  +K( 859)*Y(43 
     *        )*Y(61 )*D  +K( 860)*Y(43 )*Y(61 )*D  +K( 882)*Y(37 )*
     *        Y(47 )*D  +K( 949)*Y(12 )*X(1 )*D  +K( 950)*Y(14 )*X(1
     *         )*D  +K( 951)*Y(16 )*X(1 )*D  +K( 954)*Y(25 )*X(1 )*D
     *          +K( 955)*Y(27 )*X(1 )*D  +K( 956)*Y(32 )*X(1 )*D  +K
     *        ( 957)*Y(35 )*X(1 )*D  +K( 958)*Y(39 )*X(1 )*D  +K( 95
     *        9)*Y(41 )*X(1 )*D  +K( 960)*Y(50 )*X(1 )*D  +K( 961)*Y
     *        (53 )*X(1 )*D  +K( 962)*Y(55 )*X(1 )*D  +K( 965)*Y(73 
     *        )*X(1 )*D  +K( 966)*Y(75 )*X(1 )*D  +K( 968)*Y(2  )*X(
     *        2 )*D  +K( 969)*Y(4  )*X(2 )*D  +K( 969)*Y(4  )*X(2 )*
     *        D  +K( 970)*Y(5  )*X(2 )*D  +K( 970)*Y(5  )*X(2 )*D  +
     *        K( 970)*Y(5  )*X(2 )*D  +K( 971)*Y(5  )*X(2 )*D       
      YD( 6) = YD( 5)
     *        +K( 976)*Y(14 )*X(2 )*D  +K( 977)*Y(16 )*X(2 )*D  +K( 
     *        978)*Y(16 )*X(2 )*D  +K( 978)*Y(16 )*X(2 )*D  +K( 980)
     *        *Y(18 )*X(2 )*D  +K( 982)*Y(18 )*X(2 )*D  +K( 982)*Y(1
     *        8 )*X(2 )*D  +K( 984)*Y(20 )*X(2 )*D  +K( 985)*Y(20 )*
     *        X(2 )*D  +K( 985)*Y(20 )*X(2 )*D  +K( 987)*Y(21 )*X(2 
     *        )*D  +K( 988)*Y(21 )*X(2 )*D  +K( 988)*Y(21 )*X(2 )*D 
     *         +K( 989)*Y(21 )*X(2 )*D  +K( 992)*Y(25 )*X(2 )*D  +K(
     *         993)*Y(27 )*X(2 )*D  +K( 993)*Y(27 )*X(2 )*D  +K( 994
     *        )*Y(27 )*X(2 )*D  +K( 995)*Y(29 )*X(2 )*D  +K( 995)*Y(
     *        29 )*X(2 )*D  +K( 996)*Y(29 )*X(2 )*D  +K( 997)*Y(30 )
     *        *X(2 )*D  +K( 997)*Y(30 )*X(2 )*D  +K( 999)*Y(30 )*X(2
     *         )*D  +K(1001)*Y(33 )*X(2 )*D  +K(1005)*Y(39 )*X(2 )*D
     *          +K(1006)*Y(41 )*X(2 )*D  +K(1006)*Y(41 )*X(2 )*D  +K
     *        (1008)*Y(41 )*X(2 )*D  +K(1009)*Y(42 )*X(2 )*D  +K(101
     *        0)*Y(42 )*X(2 )*D  +K(1010)*Y(42 )*X(2 )*D  +K(1012)*Y
     *        (42 )*X(2 )*D  +K(1014)*Y(46 )*X(2 )*D  +K(1015)*Y(48 
     *        )*X(2 )*D  +K(1015)*Y(48 )*X(2 )*D  +K(1016)*Y(48 )*X(
     *        2 )*D  +K(1020)*Y(51 )*X(2 )*D                        
      YD( 7) = YD( 6)
     *        +K(1021)*Y(51 )*X(2 )*D  +K(1023)*Y(55 )*X(2 )*D  +K(1
     *        025)*Y(60 )*X(2 )*D  +K(1026)*Y(61 )*X(2 )*D  +K(1026)
     *        *Y(61 )*X(2 )*D  +K(1027)*Y(61 )*X(2 )*D  +K(1028)*Y(6
     *        1 )*X(2 )*D  +K(1030)*Y(62 )*X(2 )*D  +K(1033)*Y(64 )*
     *        X(2 )*D  +K(1035)*Y(68 )*X(2 )*D  +K(1036)*Y(70 )*X(2 
     *        )*D  +K(1037)*Y(70 )*X(2 )*D  +K(1037)*Y(70 )*X(2 )*D 
     *         +K(1039)*Y(71 )*X(2 )*D  +K(1039)*Y(71 )*X(2 )*D  +K(
     *        1040)*Y(71 )*X(2 )*D  +K(1044)*Y(76 )*X(2 )*D  +K(1045
     *        )*Y(78 )*X(2 )*D  +K(1047)*Y(79 )*X(2 )*D  +K(1047)*Y(
     *        79 )*X(2 )*D  +K(1048)*Y(79 )*X(2 )*D  +K(1055)*X(1 ) 
     *        +K(1056)*X(1 ) +K(1056)*X(1 ) +K(1059)*X(1 ) +K(1059)*
     *        X(1 ) +K(1060)*Y(3  ) +K(1061)*Y(4  ) +K(1062)*Y(5  ) 
     *        +K(1069)*Y(13 ) +K(1071)*Y(14 ) +K(1073)*Y(15 ) +K(107
     *        4)*Y(16 ) +K(1078)*Y(17 ) +K(1080)*Y(18 ) +K(1082)*Y(1
     *        9 ) +K(1084)*Y(19 ) +K(1086)*Y(20 ) +K(1088)*Y(24 ) +K
     *        (1089)*Y(26 ) +K(1093)*Y(28 ) +K(1099)*Y(38 ) +K(1100)
     *        *Y(39 ) +K(1101)*Y(40 ) +K(1103)*Y(41 ) +K(1105)*Y(45 
     *        ) +K(1107)*Y(46 ) +K(1109)*Y(47 )                     
      YD( 8) = YD( 7)
     *        +K(1109)*Y(47 ) +K(1111)*Y(47 ) +K(1114)*Y(54 ) +K(111
     *        5)*Y(56 ) +K(1118)*Y(59 ) +K(1120)*Y(67 ) +K(1121)*Y(6
     *        8 ) +K(1122)*Y(69 ) +K(1163)*Y(2  )*D +K(1164)*Y(4  )*
     *        D +K(1164)*Y(4  )*D +K(1165)*Y(5  )*D 
     *          +K(1173)*Y(14 )*D +K(1174)*Y(16 )*D +K(1174)*Y(16 )*
     *        D +K(1176)*Y(20 )*D +K(1177)*Y(21 )*D +K(1177)*Y(21 )*
     *        D +K(1178)*Y(25 )*D +K(1179)*Y(27 )*D +K(1179)*Y(27 )*
     *        D +K(1180)*Y(29 )*D +K(1181)*Y(30 )*D +K(1185)*Y(39 )*
     *        D +K(1186)*Y(41 )*D +K(1186)*Y(41 )*D +K(1187)*Y(42 )*
     *        D +K(1187)*Y(42 )*D +K(1188)*Y(46 )*D +K(1189)*Y(48 )*
     *        D +K(1189)*Y(48 )*D +K(1191)*Y(51 )*D +K(1193)*Y(55 )*
     *        D +K(1195)*Y(60 )*D +K(1196)*Y(61 )*D +K(1196)*Y(61 )*
     *        D +K(1197)*Y(62 )*D +K(1199)*Y(64 )*D +K(1200)*Y(68 )*
     *        D +K(1201)*Y(70 )*D +K(1202)*Y(71 )*D +K(1205)*Y(76 )*
     *        D +K(1207)*Y(79 )*D +K(1207)*Y(79 )*D +(Y(1  )*(-K(   
     *        1)*Y(9  )*D -K(   2)*Y(13 )*D -K(   3)*Y(24 )*D -K(   
     *        4)*Y(26 )*D -K(   5)*Y(38 )*D -K(   6)*Y(38 )*D -K(   
     *        7)*Y(47 )*D -K(   8)*Y(59 )*D ))                      
      YDOT(  1) = +YD( 8)
     *        +(Y(  1)*(
     *        -K(   9)*Y(69 )*D -K(  81)*Y(35 )*D -K(  89)*Y(4  )*D 
     *        -K(  90)*Y(7  )*D -K(  91)*Y(8  )*D -K(  92)*Y(10 )*D 
     *        -K(  93)*Y(12 )*D -K(  94)*Y(14 )*D -K(  95)*Y(15 )*D 
     *        -K(  96)*Y(20 )*D -K(  97)*Y(21 )*D -K(  98)*Y(34 )*D 
     *        -K(  99)*Y(45 )*D -K( 100)*Y(50 )*D -K( 101)*Y(53 )*D 
     *        -K( 102)*Y(55 )*D -K( 103)*Y(56 )*D -K( 104)*Y(67 )*D 
     *        -K( 105)*Y(68 )*D -K( 106)*Y(70 )*D -K( 107)*Y(71 )*D 
     *        -K( 108)*Y(2  )*D -K( 143)*Y(3  )*D -K( 967)*X(2 )*D  
     *        -K(1049) -K(1130)*D ))                             
      YD( 1) = 
     *        +K(  71)*Y(7  )*X(1 )*D  +K(  81)*Y(1  )*Y(35 )*D  +K(
     *          89)*Y(1  )*Y(4  )*D  +K(  90)*Y(1  )*Y(7  )*D  +K(  
     *        93)*Y(1  )*Y(12 )*D  +K( 101)*Y(1  )*Y(53 )*D  +K( 102
     *        )*Y(1  )*Y(55 )*D  +K( 131)*Y(2  )*Y(56 )*D  +K( 242)*
     *        Y(7  )*Y(19 )*D  +K( 256)*Y(7  )*Y(40 )*D  +K( 282)*Y(
     *        7  )*Y(59 )*D  +K( 291)*Y(7  )*Y(77 )*D  +K(1049)*Y(1 
     *         ) +K(1055)*X(1 ) +K(1057)*X(1 ) +K(1061)*Y(4  ) +K(10
     *        63)*Y(5  ) +K(1070)*Y(14 ) +K(1075)*Y(16 ) +(Y(2  )*(-
     *        K(  10)*Y(34 )*D -K( 108)*Y(1  )*D -K( 109)*Y(3  )*D -
     *        K( 110)*Y(10 )*D -K( 111)*Y(13 )*D -K( 112)*Y(15 )*D -
     *        K( 113)*Y(15 )*D -K( 114)*Y(17 )*D -K( 115)*Y(19 )*D -
     *        K( 116)*Y(19 )*D -K( 117)*Y(24 )*D -K( 118)*Y(26 )*D -
     *        K( 119)*Y(28 )*D -K( 120)*Y(36 )*D -K( 121)*Y(38 )*D -
     *        K( 122)*Y(40 )*D -K( 123)*Y(45 )*D -K( 124)*Y(45 )*D -
     *        K( 125)*Y(45 )*D -K( 126)*Y(47 )*D -K( 127)*Y(47 )*D -
     *        K( 128)*Y(47 )*D -K( 129)*Y(49 )*D -K( 130)*Y(54 )*D -
     *        K( 131)*Y(56 )*D -K( 132)*Y(57 )*D -K( 133)*Y(59 )*D -
     *        K( 134)*Y(65 )*D -K( 135)*Y(67 )*D ))                 
      YDOT(  2) = +YD( 1)
     *        +(Y(  2)*(
     *        -K( 136)*Y(67 )*D -K( 137)*Y(69 )*D -K( 138)*Y(69 )*D 
     *        -K( 139)*Y(69 )*D -K( 140)*Y(72 )*D -K( 141)*Y(74 )*D 
     *        -K( 142)*Y(77 )*D -K( 968)*X(2 )*D  -K(1163)*D ))     
      YDOT(  3) =      +K( 967)*Y(1  )*X(2 )*D  +K(1057)*X(1 ) 
     *        +(Y(3  )*(-K( 109)*Y(2  )*D -K( 143)*Y(1  )*D -K( 144)
     *        *Y(4  )*D -K( 145)*Y(5  )*D -K( 146)*Y(7  )*D -K( 147)
     *        *Y(8  )*D -K( 148)*Y(9  )*D -K( 149)*Y(11 )*D -K( 150)
     *        *Y(13 )*D -K( 151)*Y(15 )*D -K( 152)*Y(17 )*D -K( 153)
     *        *Y(22 )*D -K( 154)*Y(23 )*D -K( 155)*Y(24 )*D -K( 156)
     *        *Y(26 )*D -K( 157)*Y(30 )*D -K( 158)*Y(34 )*D -K( 159)
     *        *Y(35 )*D -K( 160)*Y(38 )*D -K( 161)*Y(42 )*D -K( 162)
     *        *Y(42 )*D -K( 163)*Y(44 )*D -K( 164)*Y(45 )*D -K( 165)
     *        *Y(46 )*D -K( 166)*Y(52 )*D -K( 167)*Y(66 )*D -K(1060)
     *         ))                                                   
      YDOT(  4) =      +K( 108)*Y(1  )*Y(2  )*D  +K( 124)*Y(2  
     *        )*Y(45 )*D  +K( 945)*Y(7  )*X(1 )*D  +K(1058)*X(1 ) +K
     *        (1062)*Y(5  ) +(Y(4  )*(-K(  89)*Y(1  )*D -K( 144)*Y(3
     *          )*D -K( 168)*Y(8  )*D -K( 169)*Y(11 )*D -K( 170)*Y(1
     *        1 )*D -K( 171)*Y(13 )*D -K( 172)*Y(13 )*D -K( 173)*Y(1
     *        5 )*D -K( 174)*Y(15 )*D -K( 175)*Y(19 )*D -K( 176)*Y(1
     *        9 )*D -K( 177)*Y(19 )*D -K( 178)*Y(22 )*D -K( 179)*Y(2
     *        4 )*D -K( 180)*Y(24 )*D -K( 181)*Y(26 )*D -K( 182)*Y(2
     *        8 )*D -K( 183)*Y(31 )*D -K( 184)*Y(34 )*D -K( 185)*Y(3
     *        6 )*D -K( 186)*Y(38 )*D -K( 187)*Y(38 )*D -K( 188)*Y(4
     *        0 )*D -K( 189)*Y(40 )*D -K( 190)*Y(45 )*D -K( 191)*Y(4
     *        5 )*D -K( 192)*Y(47 )*D -K( 193)*Y(47 )*D -K( 194)*Y(4
     *        9 )*D -K( 195)*Y(52 )*D -K( 196)*Y(52 )*D -K( 197)*Y(5
     *        4 )*D -K( 198)*Y(57 )*D -K( 199)*Y(57 )*D -K( 200)*Y(6
     *        9 )*D -K( 201)*Y(69 )*D -K( 202)*Y(69 )*D -K( 203)*X(1
     *         )*D  -K( 969)*X(2 )*D  -K(1061) -K(1164)*D ))        
      YDOT(  5) =      +K( 190)*Y(4  )*Y(45 )*D  +K( 203)*Y(4  
     *        )*X(1 )*D  +K( 953)*Y(25 )*X(1 )*D  +(Y(5  )*(-K( 145)
     *        *Y(3  )*D -K( 204)*Y(8  )*D -K( 205)*Y(11 )*D -K( 206)
     *        *Y(13 )*D -K( 207)*Y(15 )*D -K( 208)*Y(17 )*D -K( 209)
     *        *Y(19 )*D -K( 210)*Y(24 )*D -K( 211)*Y(26 )*D -K( 212)
     *        *Y(28 )*D -K( 213)*Y(31 )*D -K( 214)*Y(34 )*D -K( 215)
     *        *Y(34 )*D -K( 216)*Y(38 )*D -K( 217)*Y(40 )*D -K( 218)
     *        *Y(43 )*D -K( 219)*Y(45 )*D -K( 220)*Y(49 )*D -K( 221)
     *        *Y(52 )*D -K( 222)*Y(54 )*D -K( 223)*Y(56 )*D -K( 224)
     *        *Y(57 )*D -K( 225)*Y(59 )*D -K( 226)*Y(65 )*D -K( 227)
     *        *Y(67 )*D -K( 228)*Y(69 )*D -K( 229)*Y(72 )*D -K( 230)
     *        *Y(74 )*D -K( 231)*Y(77 )*D -K( 970)*X(2 )*D  -K( 971)
     *        *X(2 )*D  -K(1062) -K(1063) -K(1165)*D ))             
      YD( 1) = 
     *        +K(  71)*Y(7  )*X(1 )*D  +K(  90)*Y(1  )*Y(7  )*D  +K(
     *         146)*Y(3  )*Y(7  )*D  +K( 232)*Y(7  )*Y(8  )*D  +K( 2
     *        33)*Y(7  )*Y(10 )*D  +K( 234)*Y(7  )*Y(11 )*D  +K( 235
     *        )*Y(7  )*Y(13 )*D  +K( 236)*Y(7  )*Y(13 )*D  +K( 237)*
     *        Y(7  )*Y(15 )*D  +K( 238)*Y(7  )*Y(15 )*D  +K( 239)*Y(
     *        7  )*Y(17 )*D  +K( 240)*Y(7  )*Y(19 )*D  +K( 241)*Y(7 
     *         )*Y(19 )*D  +K( 242)*Y(7  )*Y(19 )*D  +K( 243)*Y(7  )
     *        *Y(19 )*D  +K( 244)*Y(7  )*Y(19 )*D  +K( 245)*Y(7  )*Y
     *        (24 )*D  +K( 246)*Y(7  )*Y(26 )*D  +K( 247)*Y(7  )*Y(2
     *        6 )*D  +K( 248)*Y(7  )*Y(28 )*D  +K( 249)*Y(7  )*Y(28 
     *        )*D  +K( 250)*Y(7  )*Y(28 )*D  +K( 251)*Y(7  )*Y(31 )*
     *        D  +K( 252)*Y(7  )*Y(31 )*D  +K( 253)*Y(7  )*Y(36 )*D 
     *         +K( 254)*Y(7  )*Y(36 )*D  +K( 255)*Y(7  )*Y(38 )*D  +
     *        K( 256)*Y(7  )*Y(40 )*D  +K( 257)*Y(7  )*Y(40 )*D  +K(
     *         258)*Y(7  )*Y(40 )*D  +K( 259)*Y(7  )*Y(45 )*D  +K( 2
     *        60)*Y(7  )*Y(45 )*D  +K( 261)*Y(7  )*Y(47 )*D  +K( 262
     *        )*Y(7  )*Y(47 )*D  +K( 263)*Y(7  )*Y(47 )*D  +K( 264)*
     *        Y(7  )*Y(47 )*D  +K( 265)*Y(7  )*Y(49 )*D             
      YDOT(  6) = +YD( 1)
     *        +K( 266)*Y(7  )*Y(49 )*D  +K( 267)*Y(7  )*Y(49 )*D  +K
     *        ( 268)*Y(7  )*Y(49 )*D  +K( 269)*Y(7  )*Y(49 )*D  +K( 
     *        270)*Y(7  )*Y(52 )*D  +K( 271)*Y(7  )*Y(52 )*D  +K( 27
     *        2)*Y(7  )*Y(54 )*D  +K( 273)*Y(7  )*Y(54 )*D  +K( 274)
     *        *Y(7  )*Y(54 )*D  +K( 275)*Y(7  )*Y(54 )*D  +K( 276)*Y
     *        (7  )*Y(56 )*D  +K( 277)*Y(7  )*Y(56 )*D  +K( 278)*Y(7
     *          )*Y(56 )*D  +K( 279)*Y(7  )*Y(57 )*D  +K( 280)*Y(7  
     *        )*Y(57 )*D  +K( 281)*Y(7  )*Y(59 )*D  +K( 282)*Y(7  )*
     *        Y(59 )*D  +K( 283)*Y(7  )*Y(67 )*D  +K( 284)*Y(7  )*Y(
     *        69 )*D  +K( 285)*Y(7  )*Y(69 )*D  +K( 286)*Y(7  )*Y(69
     *         )*D  +K( 287)*Y(7  )*Y(72 )*D  +K( 288)*Y(7  )*Y(72 )
     *        *D  +K( 289)*Y(7  )*Y(74 )*D  +K( 290)*Y(7  )*Y(74 )*D
     *          +K( 291)*Y(7  )*Y(77 )*D  +K( 292)*Y(7  )*Y(77 )*D  
     *        +K( 945)*Y(7  )*X(1 )*D  +K( 972)*Y(7  )*X(2 )*D  +K(1
     *        166)*Y(7  )*D +(Y(6  )*(-K(1050) ))                   
      YD( 1) = 
     *        +K(1050)*Y(6  ) +(Y(7  )*(-K(  71)*X(1 )*D  -K(  90)*Y
     *        (1  )*D -K( 146)*Y(3  )*D -K( 232)*Y(8  )*D -K( 233)*Y
     *        (10 )*D -K( 234)*Y(11 )*D -K( 235)*Y(13 )*D -K( 236)*Y
     *        (13 )*D -K( 237)*Y(15 )*D -K( 238)*Y(15 )*D -K( 239)*Y
     *        (17 )*D -K( 240)*Y(19 )*D -K( 241)*Y(19 )*D -K( 242)*Y
     *        (19 )*D -K( 243)*Y(19 )*D -K( 244)*Y(19 )*D -K( 245)*Y
     *        (24 )*D -K( 246)*Y(26 )*D -K( 247)*Y(26 )*D -K( 248)*Y
     *        (28 )*D -K( 249)*Y(28 )*D -K( 250)*Y(28 )*D -K( 251)*Y
     *        (31 )*D -K( 252)*Y(31 )*D -K( 253)*Y(36 )*D -K( 254)*Y
     *        (36 )*D -K( 255)*Y(38 )*D -K( 256)*Y(40 )*D -K( 257)*Y
     *        (40 )*D -K( 258)*Y(40 )*D -K( 259)*Y(45 )*D -K( 260)*Y
     *        (45 )*D -K( 261)*Y(47 )*D -K( 262)*Y(47 )*D -K( 263)*Y
     *        (47 )*D -K( 264)*Y(47 )*D -K( 265)*Y(49 )*D -K( 266)*Y
     *        (49 )*D -K( 267)*Y(49 )*D -K( 268)*Y(49 )*D -K( 269)*Y
     *        (49 )*D -K( 270)*Y(52 )*D -K( 271)*Y(52 )*D -K( 272)*Y
     *        (54 )*D -K( 273)*Y(54 )*D -K( 274)*Y(54 )*D -K( 275)*Y
     *        (54 )*D -K( 276)*Y(56 )*D -K( 277)*Y(56 )*D -K( 278)*Y
     *        (56 )*D -K( 279)*Y(57 )*D ))                          
      YDOT(  7) = +YD( 1)
     *        +(Y(  7)*(
     *        -K( 280)*Y(57 )*D -K( 281)*Y(59 )*D -K( 282)*Y(59 )*D 
     *        -K( 283)*Y(67 )*D -K( 284)*Y(69 )*D -K( 285)*Y(69 )*D 
     *        -K( 286)*Y(69 )*D -K( 287)*Y(72 )*D -K( 288)*Y(72 )*D 
     *        -K( 289)*Y(74 )*D -K( 290)*Y(74 )*D -K( 291)*Y(77 )*D 
     *        -K( 292)*Y(77 )*D -K( 945)*X(1 )*D  -K( 972)*X(2 )*D  
     *        -K(1166)*D ))                                         
      YD( 1) = 
     *        +K(   2)*Y(1  )*Y(13 )*D  +K(  16)*Y(13 )*Y(22 )*D  +K
     *        (  17)*Y(13 )*Y(34 )*D  +K(  22)*Y(13 )*Y(65 )*D  +K( 
     *        110)*Y(2  )*Y(10 )*D  +K( 148)*Y(3  )*Y(9  )*D  +K( 23
     *        3)*Y(7  )*Y(10 )*D  +K( 267)*Y(7  )*Y(49 )*D  +K( 271)
     *        *Y(7  )*Y(52 )*D  +K( 276)*Y(7  )*Y(56 )*D  +K( 287)*Y
     *        (7  )*Y(72 )*D  +K( 290)*Y(7  )*Y(74 )*D  +K( 320)*Y(9
     *          )*Y(10 )*D  +K( 320)*Y(9  )*Y(10 )*D  +K( 321)*Y(9  
     *        )*Y(13 )*D  +K( 322)*Y(9  )*Y(15 )*D  +K( 327)*Y(9  )*
     *        Y(28 )*D  +K( 333)*Y(9  )*Y(43 )*D  +K( 335)*Y(9  )*Y(
     *        45 )*D  +K( 338)*Y(9  )*Y(47 )*D  +K( 340)*Y(9  )*Y(57
     *         )*D  +K( 341)*Y(9  )*Y(65 )*D  +K( 345)*Y(9  )*Y(69 )
     *        *D  +K( 346)*Y(9  )*Y(74 )*D  +K( 348)*Y(10 )*Y(23 )*D
     *          +K( 351)*Y(10 )*Y(35 )*D  +K( 355)*Y(10 )*Y(44 )*D  
     *        +K( 357)*Y(10 )*Y(66 )*D  +K( 371)*Y(12 )*Y(13 )*D  +K
     *        ( 380)*Y(13 )*Y(29 )*D  +K( 414)*Y(14 )*Y(28 )*D  +K( 
     *        422)*Y(14 )*Y(40 )*D  +K( 431)*Y(14 )*Y(54 )*D  +K( 43
     *        2)*Y(14 )*Y(56 )*D  +K( 434)*Y(14 )*Y(65 )*D  +K( 437)
     *        *Y(14 )*Y(69 )*D  +K( 551)*Y(22 )*Y(52 )*D            
      YD( 2) = YD( 1)
     *        +K( 552)*Y(22 )*Y(53 )*D  +K( 554)*Y(11 )*Y(23 )*D  +K
     *        ( 769)*Y(35 )*Y(52 )*D  +K( 974)*Y(9  )*X(2 )*D  +K( 9
     *        75)*Y(12 )*X(2 )*D  +K( 976)*Y(14 )*X(2 )*D  +K( 978)*
     *        Y(16 )*X(2 )*D  +K( 979)*Y(16 )*X(2 )*D  +K(1022)*Y(53
     *         )*X(2 )*D  +K(1041)*Y(73 )*X(2 )*D  +K(1042)*Y(75 )*X
     *        (2 )*D  +K(1065)*Y(10 ) +K(1066)*Y(11 ) +K(1069)*Y(13 
     *        ) +K(1070)*Y(14 ) +K(1113)*Y(52 ) +K(1125)*Y(72 ) +K(1
     *        127)*Y(73 ) +K(1128)*Y(74 ) +K(1167)*Y(9  )*D +K(1172)
     *        *Y(12 )*D +K(1173)*Y(14 )*D +K(1174)*Y(16 )*D +K(1192)
     *        *Y(53 )*D +K(1203)*Y(73 )*D +(Y(8  )*(-K(  11)*Y(22 )*
     *        D -K(  12)*Y(24 )*D -K(  13)*Y(66 )*D -K(  82)*Y(26 )*
     *        D -K(  83)*Y(26 )*D -K(  91)*Y(1  )*D -K( 147)*Y(3  )*
     *        D -K( 168)*Y(4  )*D -K( 204)*Y(5  )*D -K( 232)*Y(7  )*
     *        D -K( 293)*Y(12 )*D -K( 294)*Y(21 )*D -K( 295)*Y(24 )*
     *        D -K( 296)*Y(25 )*D -K( 297)*Y(30 )*D -K( 298)*Y(32 )*
     *        D -K( 299)*Y(33 )*D -K( 300)*Y(34 )*D -K( 301)*Y(36 )*
     *        D -K( 302)*Y(37 )*D -K( 303)*Y(37 )*D -K( 304)*Y(38 )*
     *        D -K( 305)*Y(39 )*D -K( 306)*Y(41 )*D ))              
      YDOT(  8) = +YD( 2)
     *        +(Y(  8)*(
     *        -K( 307)*Y(42 )*D -K( 308)*Y(45 )*D -K( 309)*Y(46 )*D 
     *        -K( 310)*Y(51 )*D -K( 311)*Y(53 )*D -K( 312)*Y(55 )*D 
     *        -K( 313)*Y(57 )*D -K( 314)*Y(57 )*D -K( 315)*Y(60 )*D 
     *        -K( 316)*Y(65 )*D -K( 317)*Y(67 )*D -K( 318)*Y(68 )*D 
     *        -K( 319)*Y(70 )*D -K( 946)*X(1 )*D  -K( 973)*X(2 )*D  
     *        -K(1051) -K(1064) -K(1131)*D )) +K(1204)*Y(75)*D                      
      YD( 1) = 
     *        +K(  94)*Y(1  )*Y(14 )*D  +K( 232)*Y(7  )*Y(8  )*D  +K
     *        ( 234)*Y(7  )*Y(11 )*D  +K( 235)*Y(7  )*Y(13 )*D  +K( 
     *        237)*Y(7  )*Y(15 )*D  +K( 268)*Y(7  )*Y(49 )*D  +K( 27
     *        0)*Y(7  )*Y(52 )*D  +K( 274)*Y(7  )*Y(54 )*D  +K( 277)
     *        *Y(7  )*Y(56 )*D  +K( 288)*Y(7  )*Y(72 )*D  +K( 289)*Y
     *        (7  )*Y(74 )*D  +K( 293)*Y(8  )*Y(12 )*D  +K( 298)*Y(8
     *          )*Y(32 )*D  +K( 303)*Y(8  )*Y(37 )*D  +K( 311)*Y(8  
     *        )*Y(53 )*D  +K(1051)*Y(8  ) +K(1064)*Y(8  ) +K(1067)*Y
     *        (12 ) +K(1071)*Y(14 ) +K(1076)*Y(16 ) +(Y(9  )*(-K(   
     *        1)*Y(1  )*D -K(  14)*Y(22 )*D -K( 148)*Y(3  )*D -K( 32
     *        0)*Y(10 )*D -K( 321)*Y(13 )*D -K( 322)*Y(15 )*D -K( 32
     *        3)*Y(24 )*D -K( 324)*Y(26 )*D -K( 325)*Y(28 )*D -K( 32
     *        6)*Y(28 )*D -K( 327)*Y(28 )*D -K( 328)*Y(34 )*D -K( 32
     *        9)*Y(36 )*D -K( 330)*Y(36 )*D -K( 331)*Y(38 )*D -K( 33
     *        2)*Y(40 )*D -K( 333)*Y(43 )*D -K( 334)*Y(45 )*D -K( 33
     *        5)*Y(45 )*D -K( 336)*Y(47 )*D -K( 337)*Y(47 )*D -K( 33
     *        8)*Y(47 )*D -K( 339)*Y(49 )*D -K( 340)*Y(57 )*D -K( 34
     *        1)*Y(65 )*D -K( 342)*Y(65 )*D ))                      
      YDOT(  9) = +YD( 1)
     *        +(Y(  9)*(
     *        -K( 343)*Y(67 )*D -K( 344)*Y(69 )*D -K( 345)*Y(69 )*D 
     *        -K( 346)*Y(74 )*D -K( 947)*X(1 )*D  -K( 974)*X(2 )*D  
     *        -K(1167)*D -K(1208)*D ))                              
      YDOT( 10) =      +K( 973)*Y(8  )*X(2 )*D  +(Y(10 )*(-K(  
     *        92)*Y(1  )*D -K( 110)*Y(2  )*D -K( 233)*Y(7  )*D -K( 3
     *        20)*Y(9  )*D -K( 347)*Y(22 )*D -K( 348)*Y(23 )*D -K( 3
     *        49)*Y(24 )*D -K( 350)*Y(34 )*D -K( 351)*Y(35 )*D -K( 3
     *        52)*Y(36 )*D -K( 353)*Y(38 )*D -K( 354)*Y(40 )*D -K( 3
     *        55)*Y(44 )*D -K( 356)*Y(49 )*D -K( 357)*Y(66 )*D -K( 9
     *        48)*X(1 )*D  -K(1065) -K(1162)*D ))                   
      YD( 1) = 
     *        +K(  19)*Y(13 )*Y(45 )*D  +K(  21)*Y(13 )*Y(49 )*D  +K
     *        (  27)*Y(15 )*Y(36 )*D  +K(  33)*Y(17 )*Y(34 )*D  +K( 
     *         42)*Y(22 )*Y(45 )*D  +K(  43)*Y(22 )*Y(49 )*D  +K(  5
     *        9)*Y(34 )*Y(52 )*D  +K(  63)*Y(34 )*Y(72 )*D  +K(  87)
     *        *Y(36 )*Y(52 )*D  +K(  93)*Y(1  )*Y(12 )*D  +K(  99)*Y
     *        (1  )*Y(45 )*D  +K( 124)*Y(2  )*Y(45 )*D  +K( 165)*Y(3
     *          )*Y(46 )*D  +K( 190)*Y(4  )*Y(45 )*D  +K( 265)*Y(7  
     *        )*Y(49 )*D  +K( 293)*Y(8  )*Y(12 )*D  +K( 300)*Y(8  )*
     *        Y(34 )*D  +K( 301)*Y(8  )*Y(36 )*D  +K( 304)*Y(8  )*Y(
     *        38 )*D  +K( 308)*Y(8  )*Y(45 )*D  +K( 309)*Y(8  )*Y(46
     *         )*D  +K( 314)*Y(8  )*Y(57 )*D  +K( 330)*Y(9  )*Y(36 )
     *        *D  +K( 334)*Y(9  )*Y(45 )*D  +K( 336)*Y(9  )*Y(47 )*D
     *          +K( 339)*Y(9  )*Y(49 )*D  +K( 350)*Y(10 )*Y(34 )*D  
     *        +K( 356)*Y(10 )*Y(49 )*D  +K( 356)*Y(10 )*Y(49 )*D  +K
     *        ( 362)*Y(12 )*Y(36 )*D  +K( 363)*Y(12 )*Y(45 )*D  +K( 
     *        365)*Y(12 )*Y(47 )*D  +K( 366)*Y(12 )*Y(49 )*D  +K( 36
     *        7)*Y(12 )*Y(57 )*D  +K( 368)*Y(12 )*Y(65 )*D  +K( 370)
     *        *Y(12 )*Y(69 )*D  +K( 372)*Y(12 )*Y(13 )*D            
      YD( 2) = YD( 1)
     *        +K( 383)*Y(13 )*Y(34 )*D  +K( 387)*Y(13 )*Y(36 )*D  +K
     *        ( 395)*Y(13 )*Y(46 )*D  +K( 400)*Y(13 )*Y(57 )*D  +K( 
     *        426)*Y(14 )*Y(45 )*D  +K( 428)*Y(14 )*Y(47 )*D  +K( 43
     *        0)*Y(14 )*Y(49 )*D  +K( 440)*Y(12 )*Y(15 )*D  +K( 448)
     *        *Y(15 )*Y(34 )*D  +K( 449)*Y(15 )*Y(34 )*D  +K( 460)*Y
     *        (15 )*Y(45 )*D  +K( 461)*Y(15 )*Y(46 )*D  +K( 477)*Y(1
     *        6 )*Y(45 )*D  +K( 479)*Y(16 )*Y(49 )*D  +K( 487)*Y(17 
     *        )*Y(45 )*D  +K( 495)*Y(18 )*Y(45 )*D  +K( 502)*Y(12 )*
     *        Y(19 )*D  +K( 572)*Y(23 )*Y(45 )*D  +K( 588)*Y(12 )*Y(
     *        24 )*D  +K( 606)*Y(24 )*Y(46 )*D  +K( 630)*Y(25 )*Y(49
     *         )*D  +K( 640)*Y(12 )*Y(26 )*D  +K( 651)*Y(26 )*Y(46 )
     *        *D  +K( 682)*Y(12 )*Y(28 )*D  +K( 690)*Y(28 )*Y(46 )*D
     *          +K( 710)*Y(29 )*Y(45 )*D  +K( 717)*Y(32 )*Y(45 )*D  
     *        +K( 729)*Y(12 )*Y(34 )*D  +K( 743)*Y(34 )*Y(45 )*D  +K
     *        ( 744)*Y(34 )*Y(47 )*D  +K( 745)*Y(34 )*Y(50 )*D  +K( 
     *        764)*Y(35 )*Y(45 )*D  +K( 768)*Y(35 )*Y(49 )*D  +K( 77
     *        9)*Y(12 )*Y(38 )*D  +K( 785)*Y(38 )*Y(45 )*D  +K( 787)
     *        *Y(38 )*Y(46 )*D  +K( 799)*Y(39 )*Y(45 )*D            
      YD( 3) = YD( 2)
     *        +K( 814)*Y(12 )*Y(40 )*D  +K( 818)*Y(40 )*Y(46 )*D  +K
     *        ( 836)*Y(41 )*Y(45 )*D  +K( 868)*Y(45 )*Y(45 )*D  +K( 
     *        868)*Y(45 )*Y(45 )*D  +K( 869)*Y(45 )*Y(45 )*D  +K( 87
     *        1)*Y(45 )*Y(57 )*D  +K( 873)*Y(45 )*Y(66 )*D  +K( 876)
     *        *Y(45 )*Y(46 )*D  +K( 877)*Y(46 )*Y(65 )*D  +K( 878)*Y
     *        (46 )*Y(67 )*D  +K( 879)*Y(46 )*Y(69 )*D  +K( 880)*Y(4
     *        6 )*Y(72 )*D  +K( 881)*Y(46 )*Y(74 )*D  +K( 885)*Y(47 
     *        )*Y(66 )*D  +K( 889)*Y(45 )*Y(52 )*D  +K( 890)*Y(52 )*
     *        Y(57 )*D  +K( 894)*Y(36 )*Y(53 )*D  +K( 896)*Y(45 )*Y(
     *        53 )*D  +K( 904)*Y(12 )*Y(54 )*D  +K( 907)*Y(46 )*Y(54
     *         )*D  +K( 914)*Y(45 )*Y(55 )*D  +K( 923)*Y(46 )*Y(56 )
     *        *D  +K(1014)*Y(46 )*X(2 )*D  +K(1015)*Y(48 )*X(2 )*D  
     *        +K(1018)*Y(50 )*X(2 )*D  +K(1019)*Y(51 )*X(2 )*D  +K(1
     *        021)*Y(51 )*X(2 )*D  +K(1031)*Y(63 )*X(2 )*D  +K(1105)
     *        *Y(45 ) +K(1109)*Y(47 ) +K(1110)*Y(47 ) +K(1112)*Y(49 
     *        ) +K(1188)*Y(46 )*D +K(1189)*Y(48 )*D +K(1190)*Y(50 )*
     *        D +K(1191)*Y(51 )*D +K(1198)*Y(63 )*D +K(1211)*Y(81 ) 
     *        +(Y(11 )*(-K(  64)*Y(38 )*D ))                        
      YDOT( 11) = +YD( 3)
     *        +(Y( 11)*(
     *        -K( 149)*Y(3  )*D -K( 169)*Y(4  )*D -K( 170)*Y(4  )*D 
     *        -K( 205)*Y(5  )*D -K( 234)*Y(7  )*D -K( 358)*Y(32 )*D 
     *        -K( 359)*Y(33 )*D -K( 360)*Y(51 )*D -K( 361)*Y(60 )*D 
     *        -K( 517)*Y(20 )*D -K( 527)*Y(21 )*D -K( 554)*Y(23 )*D 
     *        -K( 555)*Y(23 )*D -K( 611)*Y(25 )*D -K( 791)*Y(39 )*D 
     *        -K( 832)*Y(41 )*D -K( 893)*Y(53 )*D -K( 912)*Y(55 )*D 
     *        -K(1052) -K(1066) -K(1132)*D -K(1133)*D ))            
      YD( 1) = 
     *        +K( 123)*Y(2  )*Y(45 )*D  +K( 126)*Y(2  )*Y(47 )*D  +K
     *        ( 170)*Y(4  )*Y(11 )*D  +K( 260)*Y(7  )*Y(45 )*D  +K( 
     *        262)*Y(7  )*Y(47 )*D  +K( 266)*Y(7  )*Y(49 )*D  +K( 30
     *        2)*Y(8  )*Y(37 )*D  +K( 328)*Y(9  )*Y(34 )*D  +K( 329)
     *        *Y(9  )*Y(36 )*D  +K( 331)*Y(9  )*Y(38 )*D  +K( 339)*Y
     *        (9  )*Y(49 )*D  +K( 358)*Y(11 )*Y(32 )*D  +K( 385)*Y(1
     *        3 )*Y(35 )*D  +K( 417)*Y(14 )*Y(34 )*D  +K( 418)*Y(14 
     *        )*Y(36 )*D  +K( 421)*Y(14 )*Y(38 )*D  +K( 555)*Y(11 )*
     *        Y(23 )*D  +K( 577)*Y(23 )*Y(49 )*D  +K( 752)*Y(34 )*Y(
     *        73 )*D  +K( 893)*Y(11 )*Y(53 )*D  +K(1052)*Y(11 ) +K(1
     *        107)*Y(46 ) +(Y(12 )*(-K(  93)*Y(1  )*D -K( 293)*Y(8  
     *        )*D -K( 362)*Y(36 )*D -K( 363)*Y(45 )*D -K( 364)*Y(47 
     *        )*D -K( 365)*Y(47 )*D -K( 366)*Y(49 )*D -K( 367)*Y(57 
     *        )*D -K( 368)*Y(65 )*D -K( 369)*Y(69 )*D -K( 370)*Y(69 
     *        )*D -K( 371)*Y(13 )*D -K( 372)*Y(13 )*D -K( 439)*Y(15 
     *        )*D -K( 440)*Y(15 )*D -K( 501)*Y(19 )*D -K( 502)*Y(19 
     *        )*D -K( 587)*Y(24 )*D -K( 588)*Y(24 )*D -K( 639)*Y(26 
     *        )*D -K( 640)*Y(26 )*D -K( 681)*Y(28 )*D ))            
      YDOT( 12) = +YD( 1)
     *        +(Y( 12)*(
     *        -K( 682)*Y(28 )*D -K( 729)*Y(34 )*D -K( 778)*Y(38 )*D 
     *        -K( 779)*Y(38 )*D -K( 813)*Y(40 )*D -K( 814)*Y(40 )*D 
     *        -K( 904)*Y(54 )*D -K( 949)*X(1 )*D  -K( 975)*X(2 )*D  
     *        -K(1067) -K(1172)*D ))                                
      YD( 1) = 
     *        +K(  12)*Y(8  )*Y(24 )*D  +K(  29)*Y(15 )*Y(38 )*D  +K
     *        (  31)*Y(15 )*Y(52 )*D  +K(  91)*Y(1  )*Y(8  )*D  +K( 
     *         92)*Y(1  )*Y(10 )*D  +K(  95)*Y(1  )*Y(15 )*D  +K( 14
     *        7)*Y(3  )*Y(8  )*D  +K( 273)*Y(7  )*Y(54 )*D  +K( 308)
     *        *Y(8  )*Y(45 )*D  +K( 337)*Y(9  )*Y(47 )*D  +K( 416)*Y
     *        (14 )*Y(28 )*D  +K( 425)*Y(14 )*Y(43 )*D  +K( 427)*Y(1
     *        4 )*Y(45 )*D  +K( 433)*Y(14 )*Y(57 )*D  +K( 436)*Y(14 
     *        )*Y(65 )*D  +K( 439)*Y(12 )*Y(15 )*D  +K( 474)*Y(16 )*
     *        Y(28 )*D  +K( 482)*Y(16 )*Y(69 )*D  +K( 770)*Y(35 )*Y(
     *        54 )*D  +K( 977)*Y(16 )*X(2 )*D  +K( 981)*Y(18 )*X(2 )
     *        *D  +K( 982)*Y(18 )*X(2 )*D  +K( 990)*Y(21 )*X(2 )*D  
     *        +K(1043)*Y(76 )*X(2 )*D  +K(1046)*Y(78 )*X(2 )*D  +K(1
     *        073)*Y(15 ) +K(1075)*Y(16 ) +K(1079)*Y(17 ) +K(1084)*Y
     *        (19 ) +K(1175)*Y(18 )*D +K(1206)*Y(78 )*D +(Y(13 )*(-K
     *        (   2)*Y(1  )*D -K(  15)*Y(19 )*D -K(  16)*Y(22 )*D -K
     *        (  17)*Y(34 )*D -K(  18)*Y(36 )*D -K(  19)*Y(45 )*D -K
     *        (  20)*Y(47 )*D -K(  21)*Y(49 )*D -K(  22)*Y(65 )*D -K
     *        (  72)*X(1 )*D  -K(  73)*X(1 )*D  ))                  
      YDOT( 13) = +YD( 1)
     *        +(Y( 13)*(
     *        -K( 111)*Y(2  )*D -K( 150)*Y(3  )*D -K( 171)*Y(4  )*D 
     *        -K( 172)*Y(4  )*D -K( 206)*Y(5  )*D -K( 235)*Y(7  )*D 
     *        -K( 236)*Y(7  )*D -K( 321)*Y(9  )*D -K( 371)*Y(12 )*D 
     *        -K( 372)*Y(12 )*D -K( 373)*Y(21 )*D -K( 374)*Y(22 )*D 
     *        -K( 375)*Y(23 )*D -K( 376)*Y(23 )*D -K( 377)*Y(25 )*D 
     *        -K( 378)*Y(27 )*D -K( 379)*Y(27 )*D -K( 380)*Y(29 )*D 
     *        -K( 381)*Y(32 )*D -K( 382)*Y(33 )*D -K( 383)*Y(34 )*D 
     *        -K( 384)*Y(34 )*D -K( 385)*Y(35 )*D -K( 386)*Y(35 )*D 
     *        -K( 387)*Y(36 )*D -K( 388)*Y(37 )*D -K( 389)*Y(37 )*D 
     *        -K( 390)*Y(39 )*D -K( 391)*Y(39 )*D -K( 392)*Y(41 )*D 
     *        -K( 393)*Y(41 )*D -K( 394)*Y(42 )*D -K( 395)*Y(46 )*D 
     *        -K( 396)*Y(48 )*D -K( 397)*Y(48 )*D -K( 398)*Y(53 )*D 
     *        -K( 399)*Y(55 )*D -K( 400)*Y(57 )*D -K( 401)*Y(57 )*D 
     *        -K( 402)*Y(57 )*D -K( 403)*Y(57 )*D -K( 404)*Y(59 )*D 
     *        -K( 405)*Y(60 )*D -K( 406)*Y(61 )*D -K( 407)*Y(61 )*D 
     *        -K( 408)*Y(65 )*D -K( 409)*Y(66 )*D -K( 410)*Y(68 )*D 
     *        -K(1068) -K(1069) -K(1134)*D ))                       
      YD( 1) = 
     *        +K(   1)*Y(1  )*Y(9  )*D  +K( 111)*Y(2  )*Y(13 )*D  +K
     *        ( 112)*Y(2  )*Y(15 )*D  +K( 168)*Y(4  )*Y(8  )*D  +K( 
     *        172)*Y(4  )*Y(13 )*D  +K( 204)*Y(5  )*Y(8  )*D  +K( 23
     *        6)*Y(7  )*Y(13 )*D  +K( 238)*Y(7  )*Y(15 )*D  +K( 239)
     *        *Y(7  )*Y(17 )*D  +K( 240)*Y(7  )*Y(19 )*D  +K( 259)*Y
     *        (7  )*Y(45 )*D  +K( 272)*Y(7  )*Y(54 )*D  +K( 294)*Y(8
     *          )*Y(21 )*D  +K( 296)*Y(8  )*Y(25 )*D  +K( 299)*Y(8  
     *        )*Y(33 )*D  +K( 305)*Y(8  )*Y(39 )*D  +K( 306)*Y(8  )*
     *        Y(41 )*D  +K( 309)*Y(8  )*Y(46 )*D  +K( 310)*Y(8  )*Y(
     *        51 )*D  +K( 312)*Y(8  )*Y(55 )*D  +K( 315)*Y(8  )*Y(60
     *         )*D  +K( 321)*Y(9  )*Y(13 )*D  +K( 334)*Y(9  )*Y(45 )
     *        *D  +K( 372)*Y(12 )*Y(13 )*D  +K( 376)*Y(13 )*Y(23 )*D
     *          +K( 379)*Y(13 )*Y(27 )*D  +K( 381)*Y(13 )*Y(32 )*D  
     *        +K( 386)*Y(13 )*Y(35 )*D  +K( 389)*Y(13 )*Y(37 )*D  +K
     *        ( 391)*Y(13 )*Y(39 )*D  +K( 393)*Y(13 )*Y(41 )*D  +K( 
     *        397)*Y(13 )*Y(48 )*D  +K( 398)*Y(13 )*Y(53 )*D  +K(106
     *        8)*Y(13 ) +K(1074)*Y(16 ) +K(1081)*Y(18 ) +(Y(14 )*(-K
     *        (  94)*Y(1  )*D -K( 411)*Y(22 )*D ))                  
      YDOT( 14) = +YD( 1)
     *        +(Y( 14)*(
     *        -K( 412)*Y(24 )*D -K( 413)*Y(26 )*D -K( 414)*Y(28 )*D 
     *        -K( 415)*Y(28 )*D -K( 416)*Y(28 )*D -K( 417)*Y(34 )*D 
     *        -K( 418)*Y(36 )*D -K( 419)*Y(36 )*D -K( 420)*Y(36 )*D 
     *        -K( 421)*Y(38 )*D -K( 422)*Y(40 )*D -K( 423)*Y(40 )*D 
     *        -K( 424)*Y(40 )*D -K( 425)*Y(43 )*D -K( 426)*Y(45 )*D 
     *        -K( 427)*Y(45 )*D -K( 428)*Y(47 )*D -K( 429)*Y(47 )*D 
     *        -K( 430)*Y(49 )*D -K( 431)*Y(54 )*D -K( 432)*Y(56 )*D 
     *        -K( 433)*Y(57 )*D -K( 434)*Y(65 )*D -K( 435)*Y(65 )*D 
     *        -K( 436)*Y(65 )*D -K( 437)*Y(69 )*D -K( 438)*Y(69 )*D 
     *        -K( 950)*X(1 )*D  -K( 976)*X(2 )*D  -K(1070) -K(1071) 
     *        -K(1173)*D ))                                         
      YD( 1) = 
     *        +K(  15)*Y(13 )*Y(19 )*D  +K(  19)*Y(13 )*Y(45 )*D  +K
     *        (  20)*Y(13 )*Y(47 )*D  +K(  35)*Y(17 )*Y(38 )*D  +K( 
     *         37)*Y(17 )*Y(52 )*D  +K(  72)*Y(13 )*X(1 )*D  +K(  84
     *        )*Y(17 )*Y(26 )*D  +K( 150)*Y(3  )*Y(13 )*D  +K( 404)*
     *        Y(13 )*Y(59 )*D  +K( 429)*Y(14 )*Y(47 )*D  +K( 480)*Y(
     *        16 )*Y(57 )*D  +K( 490)*Y(18 )*Y(28 )*D  +K( 508)*Y(19
     *         )*Y(39 )*D  +K( 575)*Y(23 )*Y(47 )*D  +K( 731)*Y(21 )
     *        *Y(34 )*D  +K( 946)*Y(8  )*X(1 )*D  +K( 948)*Y(10 )*X(
     *        1 )*D  +K( 980)*Y(18 )*X(2 )*D  +K( 985)*Y(20 )*X(2 )*
     *        D  +K( 989)*Y(21 )*X(2 )*D  +K(1078)*Y(17 ) +K(1083)*Y
     *        (19 ) +(Y(15 )*(-K(  23)*Y(36 )*D -K(  24)*Y(36 )*D -K
     *        (  25)*Y(36 )*D -K(  26)*Y(36 )*D -K(  27)*Y(36 )*D -K
     *        (  28)*Y(38 )*D -K(  29)*Y(38 )*D -K(  30)*Y(47 )*D -K
     *        (  31)*Y(52 )*D -K(  32)*Y(57 )*D -K(  95)*Y(1  )*D -K
     *        ( 112)*Y(2  )*D -K( 113)*Y(2  )*D -K( 151)*Y(3  )*D -K
     *        ( 173)*Y(4  )*D -K( 174)*Y(4  )*D -K( 207)*Y(5  )*D -K
     *        ( 237)*Y(7  )*D -K( 238)*Y(7  )*D -K( 322)*Y(9  )*D -K
     *        ( 439)*Y(12 )*D -K( 440)*Y(12 )*D ))                  
      YDOT( 15) = +YD( 1)
     *        +(Y( 15)*(
     *        -K( 441)*Y(21 )*D -K( 442)*Y(25 )*D -K( 443)*Y(27 )*D 
     *        -K( 444)*Y(27 )*D -K( 445)*Y(29 )*D -K( 446)*Y(32 )*D 
     *        -K( 447)*Y(33 )*D -K( 448)*Y(34 )*D -K( 449)*Y(34 )*D 
     *        -K( 450)*Y(34 )*D -K( 451)*Y(35 )*D -K( 452)*Y(37 )*D 
     *        -K( 453)*Y(37 )*D -K( 454)*Y(38 )*D -K( 455)*Y(39 )*D 
     *        -K( 456)*Y(39 )*D -K( 457)*Y(41 )*D -K( 458)*Y(41 )*D 
     *        -K( 459)*Y(42 )*D -K( 460)*Y(45 )*D -K( 461)*Y(46 )*D 
     *        -K( 462)*Y(48 )*D -K( 463)*Y(48 )*D -K( 464)*Y(53 )*D 
     *        -K( 465)*Y(55 )*D -K( 466)*Y(57 )*D -K( 467)*Y(59 )*D 
     *        -K( 468)*Y(60 )*D -K( 469)*Y(61 )*D -K( 470)*Y(61 )*D 
     *        -K( 471)*Y(65 )*D -K( 472)*Y(65 )*D -K( 473)*Y(66 )*D 
     *        -K( 535)*Y(22 )*D -K( 536)*Y(22 )*D -K( 556)*Y(23 )*D 
     *        -K(1072) -K(1073) -K(1135)*D ))                       
      YD( 1) = 
     *        +K( 113)*Y(2  )*Y(15 )*D  +K( 171)*Y(4  )*Y(13 )*D  +K
     *        ( 174)*Y(4  )*Y(15 )*D  +K( 206)*Y(5  )*Y(13 )*D  +K( 
     *        241)*Y(7  )*Y(19 )*D  +K( 261)*Y(7  )*Y(47 )*D  +K( 32
     *        2)*Y(9  )*Y(15 )*D  +K( 336)*Y(9  )*Y(47 )*D  +K( 373)
     *        *Y(13 )*Y(21 )*D  +K( 377)*Y(13 )*Y(25 )*D  +K( 378)*Y
     *        (13 )*Y(27 )*D  +K( 382)*Y(13 )*Y(33 )*D  +K( 390)*Y(1
     *        3 )*Y(39 )*D  +K( 392)*Y(13 )*Y(41 )*D  +K( 394)*Y(13 
     *        )*Y(42 )*D  +K( 395)*Y(13 )*Y(46 )*D  +K( 396)*Y(13 )*
     *        Y(48 )*D  +K( 399)*Y(13 )*Y(55 )*D  +K( 405)*Y(13 )*Y(
     *        60 )*D  +K( 406)*Y(13 )*Y(61 )*D  +K( 407)*Y(13 )*Y(61
     *         )*D  +K( 410)*Y(13 )*Y(68 )*D  +K( 426)*Y(14 )*Y(45 )
     *        *D  +K( 440)*Y(12 )*Y(15 )*D  +K( 444)*Y(15 )*Y(27 )*D
     *          +K( 446)*Y(15 )*Y(32 )*D  +K( 451)*Y(15 )*Y(35 )*D  
     *        +K( 453)*Y(15 )*Y(37 )*D  +K( 456)*Y(15 )*Y(39 )*D  +K
     *        ( 458)*Y(15 )*Y(41 )*D  +K( 463)*Y(15 )*Y(48 )*D  +K( 
     *        464)*Y(15 )*Y(53 )*D  +K( 504)*Y(19 )*Y(32 )*D  +K( 55
     *        6)*Y(15 )*Y(23 )*D  +K( 947)*Y(9  )*X(1 )*D  +K( 950)*
     *        Y(14 )*X(1 )*D  +K(1072)*Y(15 )                       
      YDOT( 16) = +YD( 1)
     *        +K(1080)*Y(18 ) +K(1085)*Y(20 ) +(Y(16 )*(-K( 474)*Y(2
     *        8 )*D -K( 475)*Y(34 )*D -K( 476)*Y(36 )*D -K( 477)*Y(4
     *        5 )*D -K( 478)*Y(47 )*D -K( 479)*Y(49 )*D -K( 480)*Y(5
     *        7 )*D -K( 481)*Y(65 )*D -K( 482)*Y(69 )*D -K( 483)*Y(6
     *        9 )*D -K( 537)*Y(22 )*D -K( 951)*X(1 )*D  -K( 977)*X(2
     *         )*D  -K( 978)*X(2 )*D  -K( 979)*X(2 )*D  -K(1074) -K(
     *        1075) -K(1076) -K(1174)*D ))                          
      YD( 1) = 
     *        +K(  15)*Y(13 )*Y(19 )*D  +K(  28)*Y(15 )*Y(38 )*D  +K
     *        (  30)*Y(15 )*Y(47 )*D  +K(  39)*Y(19 )*Y(38 )*D  +K( 
     *         40)*Y(19 )*Y(52 )*D  +K(  55)*Y(19 )*Y(34 )*D  +K(  7
     *        3)*Y(13 )*X(1 )*D  +K( 151)*Y(3  )*Y(15 )*D  +K( 242)*
     *        Y(7  )*Y(19 )*D  +K( 460)*Y(15 )*Y(45 )*D  +K( 467)*Y(
     *        15 )*Y(59 )*D  +K( 478)*Y(16 )*Y(47 )*D  +K( 494)*Y(18
     *         )*Y(43 )*D  +K( 496)*Y(18 )*Y(45 )*D  +K( 498)*Y(18 )
     *        *Y(57 )*D  +K( 501)*Y(12 )*Y(19 )*D  +K( 503)*Y(19 )*Y
     *        (29 )*D  +K( 509)*Y(19 )*Y(41 )*D  +K( 510)*Y(19 )*Y(5
     *        0 )*D  +K( 513)*Y(19 )*Y(55 )*D  +K( 516)*Y(19 )*Y(73 
     *        )*D  +K( 517)*Y(11 )*Y(20 )*D  +K( 518)*Y(19 )*Y(20 )*
     *        D  +K( 519)*Y(20 )*Y(28 )*D  +K( 522)*Y(20 )*Y(40 )*D 
     *         +K( 524)*Y(20 )*Y(49 )*D  +K( 525)*Y(20 )*Y(69 )*D  +
     *        K( 983)*Y(18 )*X(2 )*D  +K( 984)*Y(20 )*X(2 )*D  +K( 9
     *        86)*Y(21 )*X(2 )*D  +K( 988)*Y(21 )*X(2 )*D  +K(1082)*
     *        Y(19 ) +K(1176)*Y(20 )*D +K(1177)*Y(21 )*D +(Y(17 )*(-
     *        K(  33)*Y(34 )*D -K(  34)*Y(38 )*D -K(  35)*Y(38 )*D -
     *        K(  36)*Y(47 )*D -K(  37)*Y(52 )*D ))                 
      YDOT( 17) = +YD( 1)
     *        +(Y( 17)*(
     *        -K(  38)*Y(69 )*D -K(  84)*Y(26 )*D -K( 114)*Y(2  )*D 
     *        -K( 152)*Y(3  )*D -K( 208)*Y(5  )*D -K( 239)*Y(7  )*D 
     *        -K( 484)*Y(34 )*D -K( 485)*Y(36 )*D -K( 486)*Y(38 )*D 
     *        -K( 487)*Y(45 )*D -K( 488)*Y(59 )*D -K( 489)*Y(66 )*D 
     *        -K( 538)*Y(22 )*D -K( 539)*Y(22 )*D -K(1077) -K(1078) 
     *        -K(1079) -K(1136)*D ))                                
      YD( 1) = 
     *        +K(  96)*Y(1  )*Y(20 )*D  +K( 114)*Y(2  )*Y(17 )*D  +K
     *        ( 115)*Y(2  )*Y(19 )*D  +K( 173)*Y(4  )*Y(15 )*D  +K( 
     *        176)*Y(4  )*Y(19 )*D  +K( 207)*Y(5  )*Y(15 )*D  +K( 24
     *        3)*Y(7  )*Y(19 )*D  +K( 428)*Y(14 )*Y(47 )*D  +K( 441)
     *        *Y(15 )*Y(21 )*D  +K( 442)*Y(15 )*Y(25 )*D  +K( 443)*Y
     *        (15 )*Y(27 )*D  +K( 445)*Y(15 )*Y(29 )*D  +K( 447)*Y(1
     *        5 )*Y(33 )*D  +K( 455)*Y(15 )*Y(39 )*D  +K( 457)*Y(15 
     *        )*Y(41 )*D  +K( 459)*Y(15 )*Y(42 )*D  +K( 461)*Y(15 )*
     *        Y(46 )*D  +K( 462)*Y(15 )*Y(48 )*D  +K( 465)*Y(15 )*Y(
     *        55 )*D  +K( 468)*Y(15 )*Y(60 )*D  +K( 469)*Y(15 )*Y(61
     *         )*D  +K( 470)*Y(15 )*Y(61 )*D  +K( 477)*Y(16 )*Y(45 )
     *        *D  +K( 505)*Y(19 )*Y(32 )*D  +K( 557)*Y(19 )*Y(23 )*D
     *          +K( 730)*Y(20 )*Y(34 )*D  +K( 755)*Y(19 )*Y(35 )*D  
     *        +K( 951)*Y(16 )*X(1 )*D  +K(1077)*Y(17 ) +K(1086)*Y(20
     *         ) +(Y(18 )*(-K( 490)*Y(28 )*D -K( 491)*Y(34 )*D -K( 4
     *        92)*Y(34 )*D -K( 493)*Y(38 )*D -K( 494)*Y(43 )*D -K( 4
     *        95)*Y(45 )*D -K( 496)*Y(45 )*D -K( 497)*Y(47 )*D -K( 4
     *        98)*Y(57 )*D -K( 499)*Y(65 )*D ))                     
      YDOT( 18) = +YD( 1)
     *        +(Y( 18)*(
     *        -K( 500)*Y(67 )*D -K( 540)*Y(22 )*D -K( 589)*Y(24 )*D 
     *        -K( 952)*X(1 )*D  -K( 980)*X(2 )*D  -K( 981)*X(2 )*D  
     *        -K( 982)*X(2 )*D  -K( 983)*X(2 )*D  -K(1080) -K(1081) 
     *        -K(1175)*D ))                                         
      YD( 1) = 
     *        +K(  34)*Y(17 )*Y(38 )*D  +K(  36)*Y(17 )*Y(47 )*D  +K
     *        (  38)*Y(17 )*Y(69 )*D  +K( 152)*Y(3  )*Y(17 )*D  +K( 
     *        294)*Y(8  )*Y(21 )*D  +K( 373)*Y(13 )*Y(21 )*D  +K( 44
     *        1)*Y(15 )*Y(21 )*D  +K( 487)*Y(17 )*Y(45 )*D  +K( 488)
     *        *Y(17 )*Y(59 )*D  +K( 497)*Y(18 )*Y(47 )*D  +K( 520)*Y
     *        (20 )*Y(28 )*D  +K( 521)*Y(20 )*Y(36 )*D  +K( 523)*Y(2
     *        0 )*Y(47 )*D  +K( 526)*Y(20 )*Y(69 )*D  +K( 527)*Y(11 
     *        )*Y(21 )*D  +K( 528)*Y(21 )*Y(40 )*D  +K( 529)*Y(21 )*
     *        Y(45 )*D  +K( 530)*Y(21 )*Y(49 )*D  +K( 531)*Y(21 )*Y(
     *        54 )*D  +K( 532)*Y(21 )*Y(56 )*D  +K( 533)*Y(21 )*Y(65
     *         )*D  +K( 534)*Y(21 )*Y(69 )*D  +K( 590)*Y(21 )*Y(24 )
     *        *D  +K( 641)*Y(21 )*Y(26 )*D  +K( 683)*Y(21 )*Y(28 )*D
     *          +K( 780)*Y(21 )*Y(38 )*D  +K( 987)*Y(21 )*X(2 )*D  +
     *        K(1212)*Y(80 ) +(Y(19 )*(-K(  15)*Y(13 )*D -K(  39)*Y(
     *        38 )*D -K(  40)*Y(52 )*D -K(  55)*Y(34 )*D -K( 115)*Y(
     *        2  )*D -K( 116)*Y(2  )*D -K( 175)*Y(4  )*D -K( 176)*Y(
     *        4  )*D -K( 177)*Y(4  )*D -K( 209)*Y(5  )*D -K( 240)*Y(
     *        7  )*D -K( 241)*Y(7  )*D ))                           
      YDOT( 19) = +YD( 1)
     *        +(Y( 19)*(
     *        -K( 242)*Y(7  )*D -K( 243)*Y(7  )*D -K( 244)*Y(7  )*D 
     *        -K( 501)*Y(12 )*D -K( 502)*Y(12 )*D -K( 503)*Y(29 )*D 
     *        -K( 504)*Y(32 )*D -K( 505)*Y(32 )*D -K( 506)*Y(33 )*D 
     *        -K( 507)*Y(39 )*D -K( 508)*Y(39 )*D -K( 509)*Y(41 )*D 
     *        -K( 510)*Y(50 )*D -K( 511)*Y(50 )*D -K( 512)*Y(51 )*D 
     *        -K( 513)*Y(55 )*D -K( 514)*Y(60 )*D -K( 515)*Y(66 )*D 
     *        -K( 516)*Y(73 )*D -K( 518)*Y(20 )*D -K( 557)*Y(23 )*D 
     *        -K( 558)*Y(23 )*D -K( 559)*Y(23 )*D -K( 560)*Y(23 )*D 
     *        -K( 755)*Y(35 )*D -K( 756)*Y(35 )*D -K(1082) -K(1083) 
     *        -K(1084) -K(1137)*D ))                                
      YDOT( 20) =      +K(  97)*Y(1  )*Y(21 )*D  +K( 116)*Y(2  
     *        )*Y(19 )*D  +K( 177)*Y(4  )*Y(19 )*D  +K( 208)*Y(5  )*
     *        Y(17 )*D  +K( 244)*Y(7  )*Y(19 )*D  +K( 495)*Y(18 )*Y(
     *        45 )*D  +K( 502)*Y(12 )*Y(19 )*D  +K( 511)*Y(19 )*Y(50
     *         )*D  +K( 560)*Y(19 )*Y(23 )*D  +K( 756)*Y(19 )*Y(35 )
     *        *D  +(Y(20 )*(-K(  74)*X(1 )*D  -K(  96)*Y(1  )*D -K( 
     *        517)*Y(11 )*D -K( 518)*Y(19 )*D -K( 519)*Y(28 )*D -K( 
     *        520)*Y(28 )*D -K( 521)*Y(36 )*D -K( 522)*Y(40 )*D -K( 
     *        523)*Y(47 )*D -K( 524)*Y(49 )*D -K( 525)*Y(69 )*D -K( 
     *        526)*Y(69 )*D -K( 730)*Y(34 )*D -K( 984)*X(2 )*D  -K( 
     *        985)*X(2 )*D  -K(1085) -K(1086) -K(1176)*D ))         
      YDOT( 21) =      +K(  74)*Y(20 )*X(1 )*D  +K( 175)*Y(4  )
     *        *Y(19 )*D  +K( 209)*Y(5  )*Y(19 )*D  +K( 506)*Y(19 )*Y
     *        (33 )*D  +K( 507)*Y(19 )*Y(39 )*D  +K( 512)*Y(19 )*Y(5
     *        1 )*D  +K( 514)*Y(19 )*Y(60 )*D  +K( 518)*Y(19 )*Y(20 
     *        )*D  +K( 952)*Y(18 )*X(1 )*D  +(Y(21 )*(-K(  97)*Y(1  
     *        )*D -K( 294)*Y(8  )*D -K( 373)*Y(13 )*D -K( 441)*Y(15 
     *        )*D -K( 527)*Y(11 )*D -K( 528)*Y(40 )*D -K( 529)*Y(45 
     *        )*D -K( 530)*Y(49 )*D -K( 531)*Y(54 )*D -K( 532)*Y(56 
     *        )*D -K( 533)*Y(65 )*D -K( 534)*Y(69 )*D -K( 590)*Y(24 
     *        )*D -K( 641)*Y(26 )*D -K( 683)*Y(28 )*D -K( 731)*Y(34 
     *        )*D -K( 780)*Y(38 )*D -K( 986)*X(2 )*D  -K( 987)*X(2 )
     *        *D  -K( 988)*X(2 )*D  -K( 989)*X(2 )*D  -K( 990)*X(2 )
     *        *D  -K(1177)*D ))                                     
      YD( 1) = 
     *        +K(   3)*Y(1  )*Y(24 )*D  +K(  12)*Y(8  )*Y(24 )*D  +K
     *        (  32)*Y(15 )*Y(57 )*D  +K(  50)*Y(24 )*Y(52 )*D  +K( 
     *         52)*Y(24 )*Y(65 )*D  +K(  59)*Y(34 )*Y(52 )*D  +K( 15
     *        4)*Y(3  )*Y(23 )*D  +K( 251)*Y(7  )*Y(31 )*D  +K( 270)
     *        *Y(7  )*Y(52 )*D  +K( 272)*Y(7  )*Y(54 )*D  +K( 274)*Y
     *        (7  )*Y(54 )*D  +K( 277)*Y(7  )*Y(56 )*D  +K( 279)*Y(7
     *          )*Y(57 )*D  +K( 296)*Y(8  )*Y(25 )*D  +K( 314)*Y(8  
     *        )*Y(57 )*D  +K( 348)*Y(10 )*Y(23 )*D  +K( 376)*Y(13 )*
     *        Y(23 )*D  +K( 377)*Y(13 )*Y(25 )*D  +K( 403)*Y(13 )*Y(
     *        57 )*D  +K( 442)*Y(15 )*Y(25 )*D  +K( 555)*Y(11 )*Y(23
     *         )*D  +K( 556)*Y(15 )*Y(23 )*D  +K( 557)*Y(19 )*Y(23 )
     *        *D  +K( 560)*Y(19 )*Y(23 )*D  +K( 562)*Y(23 )*Y(24 )*D
     *          +K( 563)*Y(23 )*Y(26 )*D  +K( 566)*Y(23 )*Y(28 )*D  
     *        +K( 569)*Y(23 )*Y(36 )*D  +K( 570)*Y(23 )*Y(38 )*D  +K
     *        ( 571)*Y(23 )*Y(40 )*D  +K( 573)*Y(23 )*Y(45 )*D  +K( 
     *        576)*Y(23 )*Y(47 )*D  +K( 578)*Y(23 )*Y(49 )*D  +K( 57
     *        9)*Y(23 )*Y(52 )*D  +K( 580)*Y(23 )*Y(54 )*D  +K( 582)
     *        *Y(23 )*Y(57 )*D  +K( 586)*Y(23 )*Y(69 )*D            
      YD( 2) = YD( 1)
     *        +K( 587)*Y(12 )*Y(24 )*D  +K( 593)*Y(24 )*Y(27 )*D  +K
     *        ( 594)*Y(24 )*Y(29 )*D  +K( 597)*Y(24 )*Y(34 )*D  +K( 
     *        603)*Y(24 )*Y(38 )*D  +K( 605)*Y(24 )*Y(41 )*D  +K( 61
     *        1)*Y(11 )*Y(25 )*D  +K( 612)*Y(24 )*Y(25 )*D  +K( 613)
     *        *Y(25 )*Y(26 )*D  +K( 614)*Y(25 )*Y(28 )*D  +K( 616)*Y
     *        (25 )*Y(31 )*D  +K( 617)*Y(25 )*Y(34 )*D  +K( 620)*Y(2
     *        5 )*Y(38 )*D  +K( 622)*Y(25 )*Y(40 )*D  +K( 626)*Y(25 
     *        )*Y(45 )*D  +K( 631)*Y(25 )*Y(49 )*D  +K( 632)*Y(25 )*
     *        Y(52 )*D  +K( 633)*Y(25 )*Y(54 )*D  +K( 634)*Y(25 )*Y(
     *        56 )*D  +K( 637)*Y(25 )*Y(65 )*D  +K( 670)*Y(27 )*Y(52
     *         )*D  +K( 735)*Y(32 )*Y(34 )*D  +K( 759)*Y(31 )*Y(35 )
     *        *D  +K( 771)*Y(35 )*Y(54 )*D  +K( 953)*Y(25 )*X(1 )*D 
     *         +K( 991)*Y(23 )*X(2 )*D  +K( 992)*Y(25 )*X(2 )*D  +K(
     *         993)*Y(27 )*X(2 )*D  +K(1000)*Y(32 )*X(2 )*D  +K(1000
     *        )*Y(32 )*X(2 )*D  +K(1002)*Y(33 )*X(2 )*D  +K(1022)*Y(
     *        53 )*X(2 )*D  +K(1024)*Y(58 )*X(2 )*D  +K(1088)*Y(24 )
     *         +K(1094)*Y(31 ) +K(1094)*Y(31 ) +K(1113)*Y(52 ) +K(11
     *        16)*Y(57 ) +K(1168)*Y(23 )*D                          
      YDOT( 22) = +YD( 2)
     *        +K(1178)*Y(25 )*D +K(1179)*Y(27 )*D +K(1182)*Y(32 )*D 
     *        +K(1182)*Y(32 )*D +K(1183)*Y(33 )*D +K(1192)*Y(53 )*D 
     *        +K(1194)*Y(58 )*D +(Y(22 )*(-K(  11)*Y(8  )*D -K(  14)
     *        *Y(9  )*D -K(  16)*Y(13 )*D -K(  41)*Y(36 )*D -K(  42)
     *        *Y(45 )*D -K(  43)*Y(49 )*D -K(  44)*Y(57 )*D -K(  45)
     *        *Y(59 )*D -K(  46)*Y(72 )*D -K(  47)*Y(23 )*D -K( 153)
     *        *Y(3  )*D -K( 178)*Y(4  )*D -K( 347)*Y(10 )*D -K( 374)
     *        *Y(13 )*D -K( 411)*Y(14 )*D -K( 535)*Y(15 )*D -K( 536)
     *        *Y(15 )*D -K( 537)*Y(16 )*D -K( 538)*Y(17 )*D -K( 539)
     *        *Y(17 )*D -K( 540)*Y(18 )*D -K( 541)*Y(24 )*D -K( 542)
     *        *Y(25 )*D -K( 543)*Y(27 )*D -K( 544)*Y(32 )*D -K( 545)
     *        *Y(37 )*D -K( 546)*Y(38 )*D -K( 547)*Y(39 )*D -K( 548)
     *        *Y(41 )*D -K( 549)*Y(41 )*D -K( 550)*Y(45 )*D -K( 551)
     *        *Y(52 )*D -K( 552)*Y(53 )*D -K( 553)*Y(77 )*D -K(1053)
     *         -K(1138)*D ))                                        
      YDOT( 23) =      +K( 245)*Y(7  )*Y(24 )*D  +K( 246)*Y(7  
     *        )*Y(26 )*D  +K( 251)*Y(7  )*Y(31 )*D  +K( 271)*Y(7  )*
     *        Y(52 )*D  +K( 273)*Y(7  )*Y(54 )*D  +K( 280)*Y(7  )*Y(
     *        57 )*D  +K( 544)*Y(22 )*Y(32 )*D  +K(1053)*Y(22 ) +(Y(
     *        23 )*(-K(  47)*Y(22 )*D -K(  75)*X(1 )*D  -K( 154)*Y(3
     *          )*D -K( 348)*Y(10 )*D -K( 375)*Y(13 )*D -K( 376)*Y(1
     *        3 )*D -K( 554)*Y(11 )*D -K( 555)*Y(11 )*D -K( 556)*Y(1
     *        5 )*D -K( 557)*Y(19 )*D -K( 558)*Y(19 )*D -K( 559)*Y(1
     *        9 )*D -K( 560)*Y(19 )*D -K( 561)*Y(24 )*D -K( 562)*Y(2
     *        4 )*D -K( 563)*Y(26 )*D -K( 564)*Y(28 )*D -K( 565)*Y(2
     *        8 )*D -K( 566)*Y(28 )*D -K( 567)*Y(36 )*D -K( 568)*Y(3
     *        6 )*D -K( 569)*Y(36 )*D -K( 570)*Y(38 )*D -K( 571)*Y(4
     *        0 )*D -K( 572)*Y(45 )*D -K( 573)*Y(45 )*D -K( 574)*Y(4
     *        7 )*D -K( 575)*Y(47 )*D -K( 576)*Y(47 )*D -K( 577)*Y(4
     *        9 )*D -K( 578)*Y(49 )*D -K( 579)*Y(52 )*D -K( 580)*Y(5
     *        4 )*D -K( 581)*Y(57 )*D -K( 582)*Y(57 )*D -K( 583)*Y(6
     *        9 )*D -K( 584)*Y(69 )*D -K( 585)*Y(69 )*D -K( 586)*Y(6
     *        9 )*D -K( 991)*X(2 )*D  -K(1168)*D -K(1209)*D ))      
      YD( 1) = 
     *        +K(   4)*Y(1  )*Y(26 )*D  +K(  16)*Y(13 )*Y(22 )*D  +K
     *        (  42)*Y(22 )*Y(45 )*D  +K(  45)*Y(22 )*Y(59 )*D  +K( 
     *         56)*Y(26 )*Y(34 )*D  +K(  60)*Y(34 )*Y(59 )*D  +K( 15
     *        3)*Y(3  )*Y(22 )*D  +K( 378)*Y(13 )*Y(27 )*D  +K( 400)
     *        *Y(13 )*Y(57 )*D  +K( 443)*Y(15 )*Y(27 )*D  +K( 564)*Y
     *        (23 )*Y(28 )*D  +K( 574)*Y(23 )*Y(47 )*D  +K( 583)*Y(2
     *        3 )*Y(69 )*D  +K( 585)*Y(23 )*Y(69 )*D  +K( 615)*Y(25 
     *        )*Y(28 )*D  +K( 619)*Y(25 )*Y(36 )*D  +K( 625)*Y(25 )*
     *        Y(40 )*D  +K( 628)*Y(25 )*Y(47 )*D  +K( 636)*Y(25 )*Y(
     *        57 )*D  +K( 638)*Y(25 )*Y(65 )*D  +K( 639)*Y(12 )*Y(26
     *         )*D  +K( 642)*Y(26 )*Y(29 )*D  +K( 659)*Y(26 )*Y(27 )
     *        *D  +K( 660)*Y(27 )*Y(28 )*D  +K( 665)*Y(27 )*Y(40 )*D
     *          +K( 667)*Y(27 )*Y(45 )*D  +K( 671)*Y(27 )*Y(54 )*D  
     *        +K( 672)*Y(27 )*Y(56 )*D  +K( 674)*Y(27 )*Y(65 )*D  +K
     *        ( 679)*Y(27 )*Y(69 )*D  +K( 825)*Y(40 )*Y(53 )*D  +K( 
     *        994)*Y(27 )*X(2 )*D  +K( 995)*Y(29 )*X(2 )*D  +K(1002)
     *        *Y(33 )*X(2 )*D  +K(1031)*Y(63 )*X(2 )*D  +K(1089)*Y(2
     *        6 ) +K(1092)*Y(28 ) +K(1183)*Y(33 )*D                 
      YDOT( 24) = +YD( 1)
     *        +K(1198)*Y(63 )*D +(Y(24 )*(-K(   3)*Y(1  )*D -K(  12)
     *        *Y(8  )*D -K(  48)*Y(36 )*D -K(  49)*Y(36 )*D -K(  50)
     *        *Y(52 )*D -K(  51)*Y(57 )*D -K(  52)*Y(65 )*D -K( 117)
     *        *Y(2  )*D -K( 155)*Y(3  )*D -K( 179)*Y(4  )*D -K( 180)
     *        *Y(4  )*D -K( 210)*Y(5  )*D -K( 245)*Y(7  )*D -K( 295)
     *        *Y(8  )*D -K( 323)*Y(9  )*D -K( 349)*Y(10 )*D -K( 412)
     *        *Y(14 )*D -K( 541)*Y(22 )*D -K( 561)*Y(23 )*D -K( 562)
     *        *Y(23 )*D -K( 587)*Y(12 )*D -K( 588)*Y(12 )*D -K( 589)
     *        *Y(18 )*D -K( 590)*Y(21 )*D -2*K( 591)*Y(24 )*D -2*K( 
     *        592)*Y(24 )*D -K( 593)*Y(27 )*D -K( 594)*Y(29 )*D -K( 
     *        595)*Y(32 )*D -K( 596)*Y(33 )*D -K( 597)*Y(34 )*D -K( 
     *        598)*Y(34 )*D -K( 599)*Y(35 )*D -K( 600)*Y(35 )*D -K( 
     *        601)*Y(37 )*D -K( 602)*Y(38 )*D -K( 603)*Y(38 )*D -K( 
     *        604)*Y(39 )*D -K( 605)*Y(41 )*D -K( 606)*Y(46 )*D -K( 
     *        607)*Y(53 )*D -K( 608)*Y(55 )*D -K( 609)*Y(57 )*D -K( 
     *        610)*Y(60 )*D -K( 612)*Y(25 )*D -K(1087) -K(1088) -K(1
     *        139)*D ))                                             
      YD( 1) = 
     *        +K(  75)*Y(23 )*X(1 )*D  +K( 117)*Y(2  )*Y(24 )*D  +K(
     *         178)*Y(4  )*Y(22 )*D  +K( 180)*Y(4  )*Y(24 )*D  +K( 2
     *        47)*Y(7  )*Y(26 )*D  +K( 248)*Y(7  )*Y(28 )*D  +K( 276
     *        )*Y(7  )*Y(56 )*D  +K( 562)*Y(23 )*Y(24 )*D  +K( 572)*
     *        Y(23 )*Y(45 )*D  +K( 584)*Y(23 )*Y(69 )*D  +K( 588)*Y(
     *        12 )*Y(24 )*D  +K( 595)*Y(24 )*Y(32 )*D  +K( 600)*Y(24
     *         )*Y(35 )*D  +K( 607)*Y(24 )*Y(53 )*D  +K(1087)*Y(24 )
     *         +(Y(25 )*(-K( 296)*Y(8  )*D -K( 377)*Y(13 )*D -K( 442
     *        )*Y(15 )*D -K( 542)*Y(22 )*D -K( 611)*Y(11 )*D -K( 612
     *        )*Y(24 )*D -K( 613)*Y(26 )*D -K( 614)*Y(28 )*D -K( 615
     *        )*Y(28 )*D -K( 616)*Y(31 )*D -K( 617)*Y(34 )*D -K( 618
     *        )*Y(36 )*D -K( 619)*Y(36 )*D -K( 620)*Y(38 )*D -K( 621
     *        )*Y(40 )*D -K( 622)*Y(40 )*D -K( 623)*Y(40 )*D -K( 624
     *        )*Y(40 )*D -K( 625)*Y(40 )*D -K( 626)*Y(45 )*D -K( 627
     *        )*Y(47 )*D -K( 628)*Y(47 )*D -K( 629)*Y(49 )*D -K( 630
     *        )*Y(49 )*D -K( 631)*Y(49 )*D -K( 632)*Y(52 )*D -K( 633
     *        )*Y(54 )*D -K( 634)*Y(56 )*D -K( 635)*Y(57 )*D -K( 636
     *        )*Y(57 )*D -K( 637)*Y(65 )*D ))                       
      YDOT( 25) = +YD( 1)
     *        +(Y( 25)*(
     *        -K( 638)*Y(65 )*D -K( 953)*X(1 )*D  -K( 954)*X(1 )*D  
     *        -K( 992)*X(2 )*D  -K(1178)*D ))                       
      YD( 1) = 
     *        +K(  57)*Y(28 )*Y(34 )*D  +K(  65)*Y(28 )*Y(38 )*D  +K
     *        (  85)*Y(28 )*Y(52 )*D  +K( 155)*Y(3  )*Y(24 )*D  +K( 
     *        379)*Y(13 )*Y(27 )*D  +K( 444)*Y(15 )*Y(27 )*D  +K( 44
     *        5)*Y(15 )*Y(29 )*D  +K( 627)*Y(25 )*Y(47 )*D  +K( 661)
     *        *Y(27 )*Y(28 )*D  +K( 668)*Y(27 )*Y(45 )*D  +K( 673)*Y
     *        (27 )*Y(57 )*D  +K( 675)*Y(27 )*Y(65 )*D  +K( 680)*Y(2
     *        7 )*Y(69 )*D  +K( 681)*Y(12 )*Y(28 )*D  +K( 695)*Y(28 
     *        )*Y(55 )*D  +K( 707)*Y(28 )*Y(29 )*D  +K( 996)*Y(29 )*
     *        X(2 )*D  +K( 997)*Y(30 )*X(2 )*D  +K( 998)*Y(30 )*X(2 
     *        )*D  +K(1093)*Y(28 ) +K(1180)*Y(29 )*D +(Y(26 )*(-K(  
     *         4)*Y(1  )*D -K(  53)*Y(38 )*D -K(  54)*Y(57 )*D -K(  
     *        56)*Y(34 )*D -K(  76)*X(1 )*D  -K(  82)*Y(8  )*D -K(  
     *        83)*Y(8  )*D -K(  84)*Y(17 )*D -K(  86)*Y(34 )*D -K( 1
     *        18)*Y(2  )*D -K( 156)*Y(3  )*D -K( 181)*Y(4  )*D -K( 2
     *        11)*Y(5  )*D -K( 246)*Y(7  )*D -K( 247)*Y(7  )*D -K( 3
     *        24)*Y(9  )*D -K( 413)*Y(14 )*D -K( 563)*Y(23 )*D -K( 6
     *        13)*Y(25 )*D -K( 639)*Y(12 )*D -K( 640)*Y(12 )*D -K( 6
     *        41)*Y(21 )*D -K( 642)*Y(29 )*D ))                     
      YDOT( 26) = +YD( 1)
     *        +(Y( 26)*(
     *        -K( 643)*Y(32 )*D -K( 644)*Y(33 )*D -K( 645)*Y(37 )*D 
     *        -K( 646)*Y(39 )*D -K( 647)*Y(39 )*D -K( 648)*Y(41 )*D 
     *        -K( 649)*Y(41 )*D -K( 650)*Y(42 )*D -K( 651)*Y(46 )*D 
     *        -K( 652)*Y(48 )*D -K( 653)*Y(53 )*D -K( 654)*Y(55 )*D 
     *        -K( 655)*Y(57 )*D -K( 656)*Y(60 )*D -K( 657)*Y(61 )*D 
     *        -K( 658)*Y(61 )*D -K( 659)*Y(27 )*D -K( 732)*Y(34 )*D 
     *        -K( 757)*Y(35 )*D -K(1089) -K(1090) -K(1140)*D ))     
      YD( 1) = 
     *        +K( 118)*Y(2  )*Y(26 )*D  +K( 179)*Y(4  )*Y(24 )*D  +K
     *        ( 181)*Y(4  )*Y(26 )*D  +K( 210)*Y(5  )*Y(24 )*D  +K( 
     *        249)*Y(7  )*Y(28 )*D  +K( 563)*Y(23 )*Y(26 )*D  +K( 56
     *        4)*Y(23 )*Y(28 )*D  +K( 590)*Y(21 )*Y(24 )*D  +K( 596)
     *        *Y(24 )*Y(33 )*D  +K( 604)*Y(24 )*Y(39 )*D  +K( 606)*Y
     *        (24 )*Y(46 )*D  +K( 608)*Y(24 )*Y(55 )*D  +K( 610)*Y(2
     *        4 )*Y(60 )*D  +K( 612)*Y(24 )*Y(25 )*D  +K( 624)*Y(25 
     *        )*Y(40 )*D  +K( 640)*Y(12 )*Y(26 )*D  +K( 643)*Y(26 )*
     *        Y(32 )*D  +K( 645)*Y(26 )*Y(37 )*D  +K( 647)*Y(26 )*Y(
     *        39 )*D  +K( 649)*Y(26 )*Y(41 )*D  +K( 653)*Y(26 )*Y(53
     *         )*D  +K( 757)*Y(26 )*Y(35 )*D  +K( 954)*Y(25 )*X(1 )*
     *        D  +K(1090)*Y(26 ) +(Y(27 )*(-K( 378)*Y(13 )*D -K( 379
     *        )*Y(13 )*D -K( 443)*Y(15 )*D -K( 444)*Y(15 )*D -K( 543
     *        )*Y(22 )*D -K( 593)*Y(24 )*D -K( 659)*Y(26 )*D -K( 660
     *        )*Y(28 )*D -K( 661)*Y(28 )*D -K( 662)*Y(36 )*D -K( 663
     *        )*Y(36 )*D -K( 664)*Y(40 )*D -K( 665)*Y(40 )*D -K( 666
     *        )*Y(40 )*D -K( 667)*Y(45 )*D -K( 668)*Y(45 )*D -K( 669
     *        )*Y(47 )*D -K( 670)*Y(52 )*D ))                       
      YDOT( 27) = +YD( 1)
     *        +(Y( 27)*(
     *        -K( 671)*Y(54 )*D -K( 672)*Y(56 )*D -K( 673)*Y(57 )*D 
     *        -K( 674)*Y(65 )*D -K( 675)*Y(65 )*D -K( 676)*Y(69 )*D 
     *        -K( 677)*Y(69 )*D -K( 678)*Y(69 )*D -K( 679)*Y(69 )*D 
     *        -K( 680)*Y(69 )*D -K( 733)*Y(34 )*D -K( 955)*X(1 )*D  
     *        -K( 993)*X(2 )*D  -K( 994)*X(2 )*D  -K(1179)*D ))     
      YD( 1) = 
     *        +K(  53)*Y(26 )*Y(38 )*D  +K(  76)*Y(26 )*X(1 )*D  +K(
     *          84)*Y(17 )*Y(26 )*D  +K( 156)*Y(3  )*Y(26 )*D  +K( 1
     *        57)*Y(3  )*Y(30 )*D  +K( 677)*Y(27 )*Y(69 )*D  +K( 709
     *        )*Y(29 )*Y(43 )*D  +K( 711)*Y(29 )*Y(45 )*D  +K( 713)*
     *        Y(29 )*Y(57 )*D  +K( 999)*Y(30 )*X(2 )*D  +K(1181)*Y(3
     *        0 )*D +K(1213)*Y(82 ) +(Y(28 )*(-K(  57)*Y(34 )*D -K( 
     *         65)*Y(38 )*D -K(  85)*Y(52 )*D -K( 119)*Y(2  )*D -K( 
     *        182)*Y(4  )*D -K( 212)*Y(5  )*D -K( 248)*Y(7  )*D -K( 
     *        249)*Y(7  )*D -K( 250)*Y(7  )*D -K( 325)*Y(9  )*D -K( 
     *        326)*Y(9  )*D -K( 327)*Y(9  )*D -K( 414)*Y(14 )*D -K( 
     *        415)*Y(14 )*D -K( 416)*Y(14 )*D -K( 474)*Y(16 )*D -K( 
     *        490)*Y(18 )*D -K( 519)*Y(20 )*D -K( 520)*Y(20 )*D -K( 
     *        564)*Y(23 )*D -K( 565)*Y(23 )*D -K( 566)*Y(23 )*D -K( 
     *        614)*Y(25 )*D -K( 615)*Y(25 )*D -K( 660)*Y(27 )*D -K( 
     *        661)*Y(27 )*D -K( 681)*Y(12 )*D -K( 682)*Y(12 )*D -K( 
     *        683)*Y(21 )*D -K( 684)*Y(32 )*D -K( 685)*Y(33 )*D -K( 
     *        686)*Y(37 )*D -K( 687)*Y(41 )*D -K( 688)*Y(41 )*D -K( 
     *        689)*Y(42 )*D -K( 690)*Y(46 )*D ))                    
      YDOT( 28) = +YD( 1)
     *        +(Y( 28)*(
     *        -K( 691)*Y(48 )*D -K( 692)*Y(48 )*D -K( 693)*Y(50 )*D 
     *        -K( 694)*Y(51 )*D -K( 695)*Y(55 )*D -K( 696)*Y(55 )*D 
     *        -K( 697)*Y(60 )*D -K( 698)*Y(61 )*D -K( 699)*Y(61 )*D 
     *        -K( 700)*Y(66 )*D -K( 701)*Y(68 )*D -K( 702)*Y(68 )*D 
     *        -K( 703)*Y(70 )*D -K( 704)*Y(70 )*D -K( 705)*Y(71 )*D 
     *        -K( 706)*Y(78 )*D -K( 707)*Y(29 )*D -K( 758)*Y(35 )*D 
     *        -K( 792)*Y(39 )*D -K( 793)*Y(39 )*D -K(1091) -K(1092) 
     *        -K(1093) -K(1141)*D ))                                
      YD( 1) = 
     *        +K( 119)*Y(2  )*Y(28 )*D  +K( 182)*Y(4  )*Y(28 )*D  +K
     *        ( 211)*Y(5  )*Y(26 )*D  +K( 250)*Y(7  )*Y(28 )*D  +K( 
     *        327)*Y(9  )*Y(28 )*D  +K( 416)*Y(14 )*Y(28 )*D  +K( 52
     *        0)*Y(20 )*Y(28 )*D  +K( 566)*Y(23 )*Y(28 )*D  +K( 593)
     *        *Y(24 )*Y(27 )*D  +K( 613)*Y(25 )*Y(26 )*D  +K( 615)*Y
     *        (25 )*Y(28 )*D  +K( 621)*Y(25 )*Y(40 )*D  +K( 641)*Y(2
     *        1 )*Y(26 )*D  +K( 644)*Y(26 )*Y(33 )*D  +K( 646)*Y(26 
     *        )*Y(39 )*D  +K( 648)*Y(26 )*Y(41 )*D  +K( 650)*Y(26 )*
     *        Y(42 )*D  +K( 651)*Y(26 )*Y(46 )*D  +K( 652)*Y(26 )*Y(
     *        48 )*D  +K( 654)*Y(26 )*Y(55 )*D  +K( 656)*Y(26 )*Y(60
     *         )*D  +K( 657)*Y(26 )*Y(61 )*D  +K( 658)*Y(26 )*Y(61 )
     *        *D  +K( 659)*Y(26 )*Y(27 )*D  +K( 661)*Y(27 )*Y(28 )*D
     *          +K( 664)*Y(27 )*Y(40 )*D  +K( 669)*Y(27 )*Y(47 )*D  
     *        +K( 678)*Y(27 )*Y(69 )*D  +K( 682)*Y(12 )*Y(28 )*D  +K
     *        ( 684)*Y(28 )*Y(32 )*D  +K( 686)*Y(28 )*Y(37 )*D  +K( 
     *        688)*Y(28 )*Y(41 )*D  +K( 692)*Y(28 )*Y(48 )*D  +K( 69
     *        3)*Y(28 )*Y(50 )*D  +K( 696)*Y(28 )*Y(55 )*D  +K( 700)
     *        *Y(28 )*Y(66 )*D  +K( 702)*Y(28 )*Y(68 )*D            
      YDOT( 29) = +YD( 1)
     *        +K( 704)*Y(28 )*Y(70 )*D  +K( 758)*Y(28 )*Y(35 )*D  +K
     *        ( 793)*Y(28 )*Y(39 )*D  +K( 955)*Y(27 )*X(1 )*D  +K(10
     *        91)*Y(28 ) +(Y(29 )*(-K(  88)*X(1 )*D  -K( 380)*Y(13 )
     *        *D -K( 445)*Y(15 )*D -K( 503)*Y(19 )*D -K( 594)*Y(24 )
     *        *D -K( 642)*Y(26 )*D -K( 707)*Y(28 )*D -K( 708)*Y(40 )
     *        *D -K( 709)*Y(43 )*D -K( 710)*Y(45 )*D -K( 711)*Y(45 )
     *        *D -K( 712)*Y(47 )*D -K( 713)*Y(57 )*D -K( 714)*Y(69 )
     *        *D -K( 734)*Y(34 )*D -K( 781)*Y(38 )*D -K( 995)*X(2 )*
     *        D  -K( 996)*X(2 )*D  -K(1180)*D ))                    
      YD( 1) = 
     *        +K(  88)*Y(29 )*X(1 )*D  +K( 212)*Y(5  )*Y(28 )*D  +K(
     *         380)*Y(13 )*Y(29 )*D  +K( 414)*Y(14 )*Y(28 )*D  +K( 4
     *        74)*Y(16 )*Y(28 )*D  +K( 490)*Y(18 )*Y(28 )*D  +K( 503
     *        )*Y(19 )*Y(29 )*D  +K( 519)*Y(20 )*Y(28 )*D  +K( 594)*
     *        Y(24 )*Y(29 )*D  +K( 614)*Y(25 )*Y(28 )*D  +K( 642)*Y(
     *        26 )*Y(29 )*D  +K( 660)*Y(27 )*Y(28 )*D  +K( 666)*Y(27
     *         )*Y(40 )*D  +K( 676)*Y(27 )*Y(69 )*D  +K( 683)*Y(21 )
     *        *Y(28 )*D  +K( 685)*Y(28 )*Y(33 )*D  +K( 687)*Y(28 )*Y
     *        (41 )*D  +K( 689)*Y(28 )*Y(42 )*D  +K( 690)*Y(28 )*Y(4
     *        6 )*D  +K( 691)*Y(28 )*Y(48 )*D  +K( 694)*Y(28 )*Y(51 
     *        )*D  +K( 697)*Y(28 )*Y(60 )*D  +K( 698)*Y(28 )*Y(61 )*
     *        D  +K( 699)*Y(28 )*Y(61 )*D  +K( 701)*Y(28 )*Y(68 )*D 
     *         +K( 703)*Y(28 )*Y(70 )*D  +K( 705)*Y(28 )*Y(71 )*D  +
     *        K( 706)*Y(28 )*Y(78 )*D  +K( 707)*Y(28 )*Y(29 )*D  +K(
     *         708)*Y(29 )*Y(40 )*D  +K( 710)*Y(29 )*Y(45 )*D  +K( 7
     *        12)*Y(29 )*Y(47 )*D  +K( 714)*Y(29 )*Y(69 )*D  +K( 781
     *        )*Y(29 )*Y(38 )*D  +K( 792)*Y(28 )*Y(39 )*D  +(Y(30 )*
     *        (-K( 157)*Y(3  )*D -K( 297)*Y(8  )*D ))               
      YDOT( 30) = +YD( 1)
     *        +(Y( 30)*(
     *        -K( 997)*X(2 )*D  -K( 998)*X(2 )*D  -K( 999)*X(2 )*D  
     *        -K(1181)*D ))                                         
      YD( 1) = 
     *        +K(  44)*Y(22 )*Y(57 )*D  +K(  51)*Y(24 )*Y(57 )*D  +K
     *        (  54)*Y(26 )*Y(57 )*D  +K( 298)*Y(8  )*Y(32 )*D  +K( 
     *        299)*Y(8  )*Y(33 )*D  +K( 358)*Y(11 )*Y(32 )*D  +K( 35
     *        9)*Y(11 )*Y(33 )*D  +K( 381)*Y(13 )*Y(32 )*D  +K( 382)
     *        *Y(13 )*Y(33 )*D  +K( 446)*Y(15 )*Y(32 )*D  +K( 447)*Y
     *        (15 )*Y(33 )*D  +K( 504)*Y(19 )*Y(32 )*D  +K( 505)*Y(1
     *        9 )*Y(32 )*D  +K( 506)*Y(19 )*Y(33 )*D  +K( 541)*Y(22 
     *        )*Y(24 )*D  +K( 544)*Y(22 )*Y(32 )*D  +K( 551)*Y(22 )*
     *        Y(52 )*D  +K( 591)*Y(24 )*Y(24 )*D  +K( 592)*Y(24 )*Y(
     *        24 )*D  +K( 595)*Y(24 )*Y(32 )*D  +K( 596)*Y(24 )*Y(33
     *         )*D  +K( 609)*Y(24 )*Y(57 )*D  +K( 643)*Y(26 )*Y(32 )
     *        *D  +K( 644)*Y(26 )*Y(33 )*D  +K( 655)*Y(26 )*Y(57 )*D
     *          +K( 684)*Y(28 )*Y(32 )*D  +K( 685)*Y(28 )*Y(33 )*D  
     *        +K( 716)*Y(32 )*Y(36 )*D  +K( 718)*Y(32 )*Y(45 )*D  +K
     *        ( 719)*Y(32 )*Y(47 )*D  +K( 720)*Y(32 )*Y(47 )*D  +K( 
     *        721)*Y(32 )*Y(49 )*D  +K( 722)*Y(32 )*Y(57 )*D  +K( 72
     *        3)*Y(32 )*Y(65 )*D  +K( 724)*Y(32 )*Y(69 )*D  +K( 725)
     *        *Y(32 )*Y(69 )*D  +K( 726)*Y(32 )*Y(69 )*D            
      YDOT( 31) = +YD( 1)
     *        +K( 727)*Y(33 )*Y(49 )*D  +K( 728)*Y(33 )*Y(65 )*D  +K
     *        ( 736)*Y(32 )*Y(34 )*D  +K( 737)*Y(33 )*Y(34 )*D  +K( 
     *        782)*Y(32 )*Y(38 )*D  +K( 783)*Y(33 )*Y(38 )*D  +K( 81
     *        6)*Y(32 )*Y(40 )*D  +K( 817)*Y(33 )*Y(40 )*D  +K( 854)
     *        *Y(32 )*Y(43 )*D  +K( 865)*Y(33 )*Y(45 )*D  +K( 888)*Y
     *        (32 )*Y(52 )*D  +K( 890)*Y(52 )*Y(57 )*D  +K( 905)*Y(3
     *        2 )*Y(54 )*D  +K( 906)*Y(33 )*Y(54 )*D  +K( 922)*Y(33 
     *        )*Y(56 )*D  +K(1001)*Y(33 )*X(2 )*D  +K(1214)*Y(83 ) +
     *        (Y(31 )*(-K( 183)*Y(4  )*D -K( 213)*Y(5  )*D -K( 251)*
     *        Y(7  )*D -K( 252)*Y(7  )*D -K( 616)*Y(25 )*D -K( 715)*
     *        Y(60 )*D -K( 759)*Y(35 )*D -K( 794)*Y(39 )*D -K(1094) 
     *        -K(1142)*D ))                                         
      YDOT( 32) =      +K(  47)*Y(22 )*Y(23 )*D  +K( 252)*Y(7  
     *        )*Y(31 )*D  +K( 542)*Y(22 )*Y(25 )*D  +K( 552)*Y(22 )*
     *        Y(53 )*D  +K( 561)*Y(23 )*Y(24 )*D  +K( 581)*Y(23 )*Y(
     *        57 )*D  +(Y(32 )*(-K( 298)*Y(8  )*D -K( 358)*Y(11 )*D 
     *        -K( 381)*Y(13 )*D -K( 446)*Y(15 )*D -K( 504)*Y(19 )*D 
     *        -K( 505)*Y(19 )*D -K( 544)*Y(22 )*D -K( 595)*Y(24 )*D 
     *        -K( 643)*Y(26 )*D -K( 684)*Y(28 )*D -K( 716)*Y(36 )*D 
     *        -K( 717)*Y(45 )*D -K( 718)*Y(45 )*D -K( 719)*Y(47 )*D 
     *        -K( 720)*Y(47 )*D -K( 721)*Y(49 )*D -K( 722)*Y(57 )*D 
     *        -K( 723)*Y(65 )*D -K( 724)*Y(69 )*D -K( 725)*Y(69 )*D 
     *        -K( 726)*Y(69 )*D -K( 735)*Y(34 )*D -K( 736)*Y(34 )*D 
     *        -K( 782)*Y(38 )*D -K( 815)*Y(40 )*D -K( 816)*Y(40 )*D 
     *        -K( 854)*Y(43 )*D -K( 888)*Y(52 )*D -K( 905)*Y(54 )*D 
     *        -K( 956)*X(1 )*D  -K(1000)*X(2 )*D  -K(1182)*D ))     
      YDOT( 33) =      +K( 183)*Y(4  )*Y(31 )*D  +K( 213)*Y(5  
     *        )*Y(31 )*D  +K( 543)*Y(22 )*Y(27 )*D  +K( 565)*Y(23 )*
     *        Y(28 )*D  +K( 616)*Y(25 )*Y(31 )*D  +K( 635)*Y(25 )*Y(
     *        57 )*D  +K( 715)*Y(31 )*Y(60 )*D  +K( 717)*Y(32 )*Y(45
     *         )*D  +K( 794)*Y(31 )*Y(39 )*D  +K( 815)*Y(32 )*Y(40 )
     *        *D  +K( 956)*Y(32 )*X(1 )*D  +(Y(33 )*(-K( 299)*Y(8  )
     *        *D -K( 359)*Y(11 )*D -K( 382)*Y(13 )*D -K( 447)*Y(15 )
     *        *D -K( 506)*Y(19 )*D -K( 596)*Y(24 )*D -K( 644)*Y(26 )
     *        *D -K( 685)*Y(28 )*D -K( 727)*Y(49 )*D -K( 728)*Y(65 )
     *        *D -K( 737)*Y(34 )*D -K( 783)*Y(38 )*D -K( 817)*Y(40 )
     *        *D -K( 865)*Y(45 )*D -K( 906)*Y(54 )*D -K( 922)*Y(56 )
     *        *D -K(1001)*X(2 )*D  -K(1002)*X(2 )*D  -K(1183)*D ))  
      YD( 1) = 
     *        +K(   5)*Y(1  )*Y(38 )*D  +K(  18)*Y(13 )*Y(36 )*D  +K
     *        (  26)*Y(15 )*Y(36 )*D  +K(  28)*Y(15 )*Y(38 )*D  +K( 
     *         34)*Y(17 )*Y(38 )*D  +K(  41)*Y(22 )*Y(36 )*D  +K(  4
     *        4)*Y(22 )*Y(57 )*D  +K(  49)*Y(24 )*Y(36 )*D  +K(  53)
     *        *Y(26 )*Y(38 )*D  +K(  66)*Y(38 )*Y(38 )*D  +K(  67)*Y
     *        (38 )*Y(52 )*D  +K(  81)*Y(1  )*Y(35 )*D  +K( 100)*Y(1
     *          )*Y(50 )*D  +K( 129)*Y(2  )*Y(49 )*D  +K( 159)*Y(3  
     *        )*Y(35 )*D  +K( 234)*Y(7  )*Y(11 )*D  +K( 253)*Y(7  )*
     *        Y(36 )*D  +K( 259)*Y(7  )*Y(45 )*D  +K( 261)*Y(7  )*Y(
     *        47 )*D  +K( 266)*Y(7  )*Y(49 )*D  +K( 280)*Y(7  )*Y(57
     *         )*D  +K( 301)*Y(8  )*Y(36 )*D  +K( 302)*Y(8  )*Y(37 )
     *        *D  +K( 305)*Y(8  )*Y(39 )*D  +K( 313)*Y(8  )*Y(57 )*D
     *          +K( 329)*Y(9  )*Y(36 )*D  +K( 351)*Y(10 )*Y(35 )*D  
     *        +K( 386)*Y(13 )*Y(35 )*D  +K( 388)*Y(13 )*Y(37 )*D  +K
     *        ( 390)*Y(13 )*Y(39 )*D  +K( 401)*Y(13 )*Y(57 )*D  +K( 
     *        419)*Y(14 )*Y(36 )*D  +K( 451)*Y(15 )*Y(35 )*D  +K( 45
     *        2)*Y(15 )*Y(37 )*D  +K( 455)*Y(15 )*Y(39 )*D  +K( 507)
     *        *Y(19 )*Y(39 )*D  +K( 545)*Y(22 )*Y(37 )*D            
      YD( 2) = YD( 1)
     *        +K( 550)*Y(22 )*Y(45 )*D  +K( 568)*Y(23 )*Y(36 )*D  +K
     *        ( 581)*Y(23 )*Y(57 )*D  +K( 600)*Y(24 )*Y(35 )*D  +K( 
     *        601)*Y(24 )*Y(37 )*D  +K( 604)*Y(24 )*Y(39 )*D  +K( 60
     *        9)*Y(24 )*Y(57 )*D  +K( 621)*Y(25 )*Y(40 )*D  +K( 635)
     *        *Y(25 )*Y(57 )*D  +K( 646)*Y(26 )*Y(39 )*D  +K( 663)*Y
     *        (27 )*Y(36 )*D  +K( 666)*Y(27 )*Y(40 )*D  +K( 756)*Y(1
     *        9 )*Y(35 )*D  +K( 757)*Y(26 )*Y(35 )*D  +K( 758)*Y(28 
     *        )*Y(35 )*D  +K( 760)*Y(35 )*Y(36 )*D  +K( 762)*Y(35 )*
     *        Y(38 )*D  +K( 763)*Y(35 )*Y(40 )*D  +K( 765)*Y(35 )*Y(
     *        45 )*D  +K( 767)*Y(35 )*Y(47 )*D  +K( 774)*Y(35 )*Y(69
     *         )*D  +K( 778)*Y(12 )*Y(38 )*D  +K( 781)*Y(29 )*Y(38 )
     *        *D  +K( 784)*Y(38 )*Y(41 )*D  +K( 791)*Y(11 )*Y(39 )*D
     *          +K( 792)*Y(28 )*Y(39 )*D  +K( 794)*Y(31 )*Y(39 )*D  
     *        +K( 796)*Y(38 )*Y(39 )*D  +K( 797)*Y(39 )*Y(40 )*D  +K
     *        ( 800)*Y(39 )*Y(45 )*D  +K( 803)*Y(39 )*Y(49 )*D  +K( 
     *        804)*Y(39 )*Y(52 )*D  +K( 805)*Y(39 )*Y(54 )*D  +K( 80
     *        6)*Y(39 )*Y(56 )*D  +K( 807)*Y(39 )*Y(57 )*D  +K( 809)
     *        *Y(39 )*Y(65 )*D  +K( 811)*Y(39 )*Y(69 )*D            
      YD( 3) = YD( 2)
     *        +K( 824)*Y(40 )*Y(53 )*D  +K( 975)*Y(12 )*X(2 )*D  +K(
     *        1003)*Y(35 )*X(2 )*D  +K(1004)*Y(37 )*X(2 )*D  +K(1004
     *        )*Y(37 )*X(2 )*D  +K(1005)*Y(39 )*X(2 )*D  +K(1006)*Y(
     *        41 )*X(2 )*D  +K(1007)*Y(41 )*X(2 )*D  +K(1009)*Y(42 )
     *        *X(2 )*D  +K(1018)*Y(50 )*X(2 )*D  +K(1021)*Y(51 )*X(2
     *         )*D  +K(1024)*Y(58 )*X(2 )*D  +K(1066)*Y(11 ) +K(1067
     *        )*Y(12 ) +K(1095)*Y(36 ) +K(1095)*Y(36 ) +K(1097)*Y(37
     *         ) +K(1099)*Y(38 ) +K(1112)*Y(49 ) +K(1116)*Y(57 ) +K(
     *        1172)*Y(12 )*D +K(1184)*Y(37 )*D +K(1184)*Y(37 )*D +K(
     *        1185)*Y(39 )*D +K(1186)*Y(41 )*D +K(1190)*Y(50 )*D +K(
     *        1191)*Y(51 )*D +K(1194)*Y(58 )*D +(Y(34 )*(-K(  10)*Y(
     *        2  )*D -K(  17)*Y(13 )*D -K(  33)*Y(17 )*D -K(  55)*Y(
     *        19 )*D -K(  56)*Y(26 )*D -K(  57)*Y(28 )*D -K(  58)*Y(
     *        47 )*D -K(  59)*Y(52 )*D -K(  60)*Y(59 )*D -K(  61)*Y(
     *        67 )*D -K(  62)*Y(69 )*D -K(  63)*Y(72 )*D -K(  77)*X(
     *        1 )*D  -K(  86)*Y(26 )*D -K(  98)*Y(1  )*D -K( 158)*Y(
     *        3  )*D -K( 184)*Y(4  )*D -K( 214)*Y(5  )*D -K( 215)*Y(
     *        5  )*D -K( 300)*Y(8  )*D -K( 328)*Y(9  )*D ))         
      YDOT( 34) = +YD( 3)+K(1169)*Y(35 )*D
     *        +(Y( 34)*(
     *        -K( 350)*Y(10 )*D -K( 383)*Y(13 )*D -K( 384)*Y(13 )*D 
     *        -K( 417)*Y(14 )*D -K( 448)*Y(15 )*D -K( 449)*Y(15 )*D 
     *        -K( 450)*Y(15 )*D -K( 475)*Y(16 )*D -K( 484)*Y(17 )*D 
     *        -K( 491)*Y(18 )*D -K( 492)*Y(18 )*D -K( 597)*Y(24 )*D 
     *        -K( 598)*Y(24 )*D -K( 617)*Y(25 )*D -K( 729)*Y(12 )*D 
     *        -K( 730)*Y(20 )*D -K( 731)*Y(21 )*D -K( 732)*Y(26 )*D 
     *        -K( 733)*Y(27 )*D -K( 734)*Y(29 )*D -K( 735)*Y(32 )*D 
     *        -K( 736)*Y(32 )*D -K( 737)*Y(33 )*D -2*K( 738)*Y(34 )*
     *        D -K( 739)*Y(38 )*D -K( 740)*Y(39 )*D -K( 741)*Y(41 )*
     *        D -K( 742)*Y(45 )*D -K( 743)*Y(45 )*D -K( 744)*Y(47 )*
     *        D -K( 745)*Y(50 )*D -K( 746)*Y(50 )*D -K( 747)*Y(51 )*
     *        D -K( 748)*Y(53 )*D -K( 749)*Y(59 )*D -K( 750)*Y(68 )*
     *        D -K( 751)*Y(70 )*D -K( 752)*Y(73 )*D -K( 753)*Y(77 )*
     *        D -K( 754)*Y(78 )*D -K(1054) -K(1143)*D ))            
      YDOT( 35) =      +K(  10)*Y(2  )*Y(34 )*D  +K( 253)*Y(7  
     *        )*Y(36 )*D  +K( 255)*Y(7  )*Y(38 )*D  +K( 265)*Y(7  )*
     *        Y(49 )*D  +K( 279)*Y(7  )*Y(57 )*D  +K( 330)*Y(9  )*Y(
     *        36 )*D  +K( 420)*Y(14 )*Y(36 )*D  +K( 567)*Y(23 )*Y(36
     *         )*D  +K( 729)*Y(12 )*Y(34 )*D  +K( 736)*Y(32 )*Y(34 )
     *        *D  +K( 746)*Y(34 )*Y(50 )*D  +K( 748)*Y(34 )*Y(53 )*D
     *          +K(1054)*Y(34 ) +K(1097)*Y(37 ) +K(1100)*Y(39 ) +(Y(
     *        35 )*(-K(  81)*Y(1  )*D -K( 159)*Y(3  )*D -K( 351)*Y(1
     *        0 )*D -K( 385)*Y(13 )*D -K( 386)*Y(13 )*D -K( 451)*Y(1
     *        5 )*D -K( 599)*Y(24 )*D -K( 600)*Y(24 )*D -K( 755)*Y(1
     *        9 )*D -K( 756)*Y(19 )*D -K( 757)*Y(26 )*D -K( 758)*Y(2
     *        8 )*D -K( 759)*Y(31 )*D -K( 760)*Y(36 )*D -K( 761)*Y(3
     *        8 )*D -K( 762)*Y(38 )*D -K( 763)*Y(40 )*D -K( 764)*Y(4
     *        5 )*D -K( 765)*Y(45 )*D -K( 766)*Y(47 )*D -K( 767)*Y(4
     *        7 )*D -K( 768)*Y(49 )*D -K( 769)*Y(52 )*D -K( 770)*Y(5
     *        4 )*D -K( 771)*Y(54 )*D -K( 772)*Y(69 )*D -K( 773)*Y(6
     *        9 )*D -K( 774)*Y(69 )*D -K( 957)*X(1 )*D  -K(1003)*X(2
     *         )*D  -K(1169)*D -K(1210)*D ))                        
      YD( 1) = 
     *        +K(  60)*Y(34 )*Y(59 )*D  +K( 268)*Y(7  )*Y(49 )*D  +K
     *        ( 303)*Y(8  )*Y(37 )*D  +K( 389)*Y(13 )*Y(37 )*D  +K( 
     *        453)*Y(15 )*Y(37 )*D  +K( 645)*Y(26 )*Y(37 )*D  +K( 68
     *        6)*Y(28 )*Y(37 )*D  +K( 738)*Y(34 )*Y(34 )*D  +K( 739)
     *        *Y(34 )*Y(38 )*D  +K( 747)*Y(34 )*Y(51 )*D  +K( 776)*Y
     *        (37 )*Y(65 )*D  +K( 777)*Y(37 )*Y(69 )*D  +K( 855)*Y(3
     *        7 )*Y(43 )*D  +K( 867)*Y(37 )*Y(45 )*D  +K( 882)*Y(37 
     *        )*Y(47 )*D  +K( 883)*Y(37 )*Y(47 )*D  +K( 928)*Y(37 )*
     *        Y(57 )*D  +K(1215)*Y(84 ) +(Y(36 )*(-K(  18)*Y(13 )*D 
     *        -K(  23)*Y(15 )*D -K(  24)*Y(15 )*D -K(  25)*Y(15 )*D 
     *        -K(  26)*Y(15 )*D -K(  27)*Y(15 )*D -K(  41)*Y(22 )*D 
     *        -K(  48)*Y(24 )*D -K(  49)*Y(24 )*D -K(  87)*Y(52 )*D 
     *        -K( 120)*Y(2  )*D -K( 185)*Y(4  )*D -K( 253)*Y(7  )*D 
     *        -K( 254)*Y(7  )*D -K( 301)*Y(8  )*D -K( 329)*Y(9  )*D 
     *        -K( 330)*Y(9  )*D -K( 352)*Y(10 )*D -K( 362)*Y(12 )*D 
     *        -K( 387)*Y(13 )*D -K( 418)*Y(14 )*D -K( 419)*Y(14 )*D 
     *        -K( 420)*Y(14 )*D -K( 476)*Y(16 )*D -K( 485)*Y(17 )*D 
     *        -K( 521)*Y(20 )*D -K( 567)*Y(23 )*D ))                
      YDOT( 36) = +YD( 1)
     *        +(Y( 36)*(
     *        -K( 568)*Y(23 )*D -K( 569)*Y(23 )*D -K( 618)*Y(25 )*D 
     *        -K( 619)*Y(25 )*D -K( 662)*Y(27 )*D -K( 663)*Y(27 )*D 
     *        -K( 716)*Y(32 )*D -K( 760)*Y(35 )*D -K( 775)*Y(50 )*D 
     *        -K( 795)*Y(39 )*D -K( 833)*Y(41 )*D -K( 866)*Y(45 )*D 
     *        -K( 894)*Y(53 )*D -K( 895)*Y(53 )*D -K( 913)*Y(55 )*D 
     *        -K(1095) -K(1096) -K(1144)*D ))                       
      YDOT( 37) =      +K( 120)*Y(2  )*Y(36 )*D  +K( 185)*Y(4  
     *        )*Y(36 )*D  +K( 254)*Y(7  )*Y(36 )*D  +K( 267)*Y(7  )*
     *        Y(49 )*D  +K( 362)*Y(12 )*Y(36 )*D  +K( 521)*Y(20 )*Y(
     *        36 )*D  +K( 569)*Y(23 )*Y(36 )*D  +K( 619)*Y(25 )*Y(36
     *         )*D  +K( 716)*Y(32 )*Y(36 )*D  +K( 740)*Y(34 )*Y(39 )
     *        *D  +K( 741)*Y(34 )*Y(41 )*D  +K( 745)*Y(34 )*Y(50 )*D
     *          +K( 760)*Y(35 )*Y(36 )*D  +K( 761)*Y(35 )*Y(38 )*D  
     *        +K( 768)*Y(35 )*Y(49 )*D  +K( 775)*Y(36 )*Y(50 )*D  +K
     *        ( 795)*Y(36 )*Y(39 )*D  +K( 833)*Y(36 )*Y(41 )*D  +K( 
     *        895)*Y(36 )*Y(53 )*D  +K( 913)*Y(36 )*Y(55 )*D  +K(109
     *        6)*Y(36 ) +(Y(37 )*(-K( 302)*Y(8  )*D -K( 303)*Y(8  )*
     *        D -K( 388)*Y(13 )*D -K( 389)*Y(13 )*D -K( 452)*Y(15 )*
     *        D -K( 453)*Y(15 )*D -K( 545)*Y(22 )*D -K( 601)*Y(24 )*
     *        D -K( 645)*Y(26 )*D -K( 686)*Y(28 )*D -K( 776)*Y(65 )*
     *        D -K( 777)*Y(69 )*D -K( 855)*Y(43 )*D -K( 867)*Y(45 )*
     *        D -K( 882)*Y(47 )*D -K( 883)*Y(47 )*D -K( 928)*Y(57 )*
     *        D -K(1004)*X(2 )*D  -K(1097) -K(1184)*D ))            
      YD( 1) = 
     *        +K(  17)*Y(13 )*Y(34 )*D  +K(  25)*Y(15 )*Y(36 )*D  +K
     *        (  48)*Y(24 )*Y(36 )*D  +K(  51)*Y(24 )*Y(57 )*D  +K( 
     *         55)*Y(19 )*Y(34 )*D  +K(  56)*Y(26 )*Y(34 )*D  +K(  5
     *        7)*Y(28 )*Y(34 )*D  +K(  58)*Y(34 )*Y(47 )*D  +K(  61)
     *        *Y(34 )*Y(67 )*D  +K(  62)*Y(34 )*Y(69 )*D  +K(  77)*Y
     *        (34 )*X(1 )*D  +K(  98)*Y(1  )*Y(34 )*D  +K( 158)*Y(3 
     *         )*Y(34 )*D  +K( 161)*Y(3  )*Y(42 )*D  +K( 256)*Y(7  )
     *        *Y(40 )*D  +K( 306)*Y(8  )*Y(41 )*D  +K( 387)*Y(13 )*Y
     *        (36 )*D  +K( 391)*Y(13 )*Y(39 )*D  +K( 392)*Y(13 )*Y(4
     *        1 )*D  +K( 402)*Y(13 )*Y(57 )*D  +K( 418)*Y(14 )*Y(36 
     *        )*D  +K( 456)*Y(15 )*Y(39 )*D  +K( 457)*Y(15 )*Y(41 )*
     *        D  +K( 466)*Y(15 )*Y(57 )*D  +K( 476)*Y(16 )*Y(36 )*D 
     *         +K( 597)*Y(24 )*Y(34 )*D  +K( 618)*Y(25 )*Y(36 )*D  +
     *        K( 624)*Y(25 )*Y(40 )*D  +K( 647)*Y(26 )*Y(39 )*D  +K(
     *         648)*Y(26 )*Y(41 )*D  +K( 655)*Y(26 )*Y(57 )*D  +K( 6
     *        62)*Y(27 )*Y(36 )*D  +K( 664)*Y(27 )*Y(40 )*D  +K( 687
     *        )*Y(28 )*Y(41 )*D  +K( 708)*Y(29 )*Y(40 )*D  +K( 730)*
     *        Y(20 )*Y(34 )*D  +K( 743)*Y(34 )*Y(45 )*D             
      YD( 2) = YD( 1)
     *        +K( 744)*Y(34 )*Y(47 )*D  +K( 749)*Y(34 )*Y(59 )*D  +K
     *        ( 750)*Y(34 )*Y(68 )*D  +K( 751)*Y(34 )*Y(70 )*D  +K( 
     *        753)*Y(34 )*Y(77 )*D  +K( 755)*Y(19 )*Y(35 )*D  +K( 76
     *        6)*Y(35 )*Y(47 )*D  +K( 772)*Y(35 )*Y(69 )*D  +K( 793)
     *        *Y(28 )*Y(39 )*D  +K( 795)*Y(36 )*Y(39 )*D  +K( 798)*Y
     *        (39 )*Y(40 )*D  +K( 801)*Y(39 )*Y(45 )*D  +K( 802)*Y(3
     *        9 )*Y(47 )*D  +K( 808)*Y(39 )*Y(57 )*D  +K( 810)*Y(39 
     *        )*Y(65 )*D  +K( 812)*Y(39 )*Y(69 )*D  +K( 813)*Y(12 )*
     *        Y(40 )*D  +K( 815)*Y(32 )*Y(40 )*D  +K( 820)*Y(40 )*Y(
     *        50 )*D  +K( 823)*Y(40 )*Y(53 )*D  +K( 832)*Y(11 )*Y(41
     *         )*D  +K( 834)*Y(40 )*Y(41 )*D  +K( 837)*Y(41 )*Y(45 )
     *        *D  +K( 840)*Y(41 )*Y(54 )*D  +K( 841)*Y(41 )*Y(56 )*D
     *          +K( 843)*Y(41 )*Y(65 )*D  +K( 846)*Y(41 )*Y(69 )*D  
     *        +K( 866)*Y(36 )*Y(45 )*D  +K(1008)*Y(41 )*X(2 )*D  +K(
     *        1010)*Y(42 )*X(2 )*D  +K(1011)*Y(42 )*X(2 )*D  +K(1019
     *        )*Y(51 )*X(2 )*D  +K(1101)*Y(40 ) +K(1187)*Y(42 )*D +(
     *        Y(38 )*(-K(   5)*Y(1  )*D -K(   6)*Y(1  )*D -K(  28)*Y
     *        (15 )*D -K(  29)*Y(15 )*D ))                          
      YDOT( 38) = +YD( 2)
     *        +(Y( 38)*(
     *        -K(  34)*Y(17 )*D -K(  35)*Y(17 )*D -K(  39)*Y(19 )*D 
     *        -K(  53)*Y(26 )*D -K(  64)*Y(11 )*D -K(  65)*Y(28 )*D 
     *        -2*K(  66)*Y(38 )*D -K(  67)*Y(52 )*D -K(  68)*Y(69 )*
     *        D -K(  78)*X(1 )*D  -K( 121)*Y(2  )*D -K( 160)*Y(3  )*
     *        D -K( 186)*Y(4  )*D -K( 187)*Y(4  )*D -K( 216)*Y(5  )*
     *        D -K( 255)*Y(7  )*D -K( 304)*Y(8  )*D -K( 331)*Y(9  )*
     *        D -K( 353)*Y(10 )*D -K( 421)*Y(14 )*D -K( 454)*Y(15 )*
     *        D -K( 486)*Y(17 )*D -K( 493)*Y(18 )*D -K( 546)*Y(22 )*
     *        D -K( 570)*Y(23 )*D -K( 602)*Y(24 )*D -K( 603)*Y(24 )*
     *        D -K( 620)*Y(25 )*D -K( 739)*Y(34 )*D -K( 761)*Y(35 )*
     *        D -K( 762)*Y(35 )*D -K( 778)*Y(12 )*D -K( 779)*Y(12 )*
     *        D -K( 780)*Y(21 )*D -K( 781)*Y(29 )*D -K( 782)*Y(32 )*
     *        D -K( 783)*Y(33 )*D -K( 784)*Y(41 )*D -K( 785)*Y(45 )*
     *        D -K( 786)*Y(46 )*D -K( 787)*Y(46 )*D -K( 788)*Y(53 )*
     *        D -K( 789)*Y(55 )*D -K( 790)*Y(60 )*D -K( 796)*Y(39 )*
     *        D -K(1098) -K(1099) -K(1145)*D ))                     
      YD( 1) = 
     *        +K( 121)*Y(2  )*Y(38 )*D  +K( 184)*Y(4  )*Y(34 )*D  +K
     *        ( 187)*Y(4  )*Y(38 )*D  +K( 214)*Y(5  )*Y(34 )*D  +K( 
     *        257)*Y(7  )*Y(40 )*D  +K( 570)*Y(23 )*Y(38 )*D  +K( 61
     *        7)*Y(25 )*Y(34 )*D  +K( 737)*Y(33 )*Y(34 )*D  +K( 762)
     *        *Y(35 )*Y(38 )*D  +K( 764)*Y(35 )*Y(45 )*D  +K( 779)*Y
     *        (12 )*Y(38 )*D  +K( 782)*Y(32 )*Y(38 )*D  +K( 788)*Y(3
     *        8 )*Y(53 )*D  +K( 957)*Y(35 )*X(1 )*D  +K(1098)*Y(38 )
     *         +K(1103)*Y(41 ) +(Y(39 )*(-K( 305)*Y(8  )*D -K( 390)*
     *        Y(13 )*D -K( 391)*Y(13 )*D -K( 455)*Y(15 )*D -K( 456)*
     *        Y(15 )*D -K( 507)*Y(19 )*D -K( 508)*Y(19 )*D -K( 547)*
     *        Y(22 )*D -K( 604)*Y(24 )*D -K( 646)*Y(26 )*D -K( 647)*
     *        Y(26 )*D -K( 740)*Y(34 )*D -K( 791)*Y(11 )*D -K( 792)*
     *        Y(28 )*D -K( 793)*Y(28 )*D -K( 794)*Y(31 )*D -K( 795)*
     *        Y(36 )*D -K( 796)*Y(38 )*D -K( 797)*Y(40 )*D -K( 798)*
     *        Y(40 )*D -K( 799)*Y(45 )*D -K( 800)*Y(45 )*D -K( 801)*
     *        Y(45 )*D -K( 802)*Y(47 )*D -K( 803)*Y(49 )*D -K( 804)*
     *        Y(52 )*D -K( 805)*Y(54 )*D -K( 806)*Y(56 )*D -K( 807)*
     *        Y(57 )*D -K( 808)*Y(57 )*D ))                         
      YDOT( 39) = +YD( 1)
     *        +(Y( 39)*(
     *        -K( 809)*Y(65 )*D -K( 810)*Y(65 )*D -K( 811)*Y(69 )*D 
     *        -K( 812)*Y(69 )*D -K( 958)*X(1 )*D  -K(1005)*X(2 )*D  
     *        -K(1100) -K(1185)*D ))                                
      YD( 1) = 
     *        +K(   6)*Y(1  )*Y(38 )*D  +K(  27)*Y(15 )*Y(36 )*D  +K
     *        (  29)*Y(15 )*Y(38 )*D  +K(  35)*Y(17 )*Y(38 )*D  +K( 
     *         39)*Y(19 )*Y(38 )*D  +K(  54)*Y(26 )*Y(57 )*D  +K(  6
     *        5)*Y(28 )*Y(38 )*D  +K(  66)*Y(38 )*Y(38 )*D  +K(  68)
     *        *Y(38 )*Y(69 )*D  +K(  78)*Y(38 )*X(1 )*D  +K( 160)*Y(
     *        3  )*Y(38 )*D  +K( 162)*Y(3  )*Y(42 )*D  +K( 393)*Y(13
     *         )*Y(41 )*D  +K( 394)*Y(13 )*Y(42 )*D  +K( 458)*Y(15 )
     *        *Y(41 )*D  +K( 459)*Y(15 )*Y(42 )*D  +K( 485)*Y(17 )*Y
     *        (36 )*D  +K( 603)*Y(24 )*Y(38 )*D  +K( 649)*Y(26 )*Y(4
     *        1 )*D  +K( 650)*Y(26 )*Y(42 )*D  +K( 688)*Y(28 )*Y(41 
     *        )*D  +K( 689)*Y(28 )*Y(42 )*D  +K( 773)*Y(35 )*Y(69 )*
     *        D  +K( 785)*Y(38 )*Y(45 )*D  +K( 833)*Y(36 )*Y(41 )*D 
     *         +K( 835)*Y(41 )*Y(43 )*D  +K( 838)*Y(41 )*Y(45 )*D  +
     *        K( 839)*Y(41 )*Y(47 )*D  +K( 842)*Y(41 )*Y(57 )*D  +K(
     *         844)*Y(41 )*Y(65 )*D  +K( 847)*Y(41 )*Y(69 )*D  +K( 8
     *        48)*Y(42 )*Y(43 )*D  +K( 849)*Y(42 )*Y(54 )*D  +K( 850
     *        )*Y(42 )*Y(56 )*D  +K( 851)*Y(42 )*Y(69 )*D  +K( 852)*
     *        Y(42 )*Y(72 )*D  +K( 853)*Y(42 )*Y(74 )*D             
      YDOT( 40) = +YD( 1)
     *        +K(1012)*Y(42 )*X(2 )*D  +K(1216)*Y(85 ) +(Y(40 )*(-K(
     *         122)*Y(2  )*D -K( 188)*Y(4  )*D -K( 189)*Y(4  )*D -K(
     *         217)*Y(5  )*D -K( 256)*Y(7  )*D -K( 257)*Y(7  )*D -K(
     *         258)*Y(7  )*D -K( 332)*Y(9  )*D -K( 354)*Y(10 )*D -K(
     *         422)*Y(14 )*D -K( 423)*Y(14 )*D -K( 424)*Y(14 )*D -K(
     *         522)*Y(20 )*D -K( 528)*Y(21 )*D -K( 571)*Y(23 )*D -K(
     *         621)*Y(25 )*D -K( 622)*Y(25 )*D -K( 623)*Y(25 )*D -K(
     *         624)*Y(25 )*D -K( 625)*Y(25 )*D -K( 664)*Y(27 )*D -K(
     *         665)*Y(27 )*D -K( 666)*Y(27 )*D -K( 708)*Y(29 )*D -K(
     *         763)*Y(35 )*D -K( 797)*Y(39 )*D -K( 798)*Y(39 )*D -K(
     *         813)*Y(12 )*D -K( 814)*Y(12 )*D -K( 815)*Y(32 )*D -K(
     *         816)*Y(32 )*D -K( 817)*Y(33 )*D -K( 818)*Y(46 )*D -K(
     *         819)*Y(48 )*D -K( 820)*Y(50 )*D -K( 821)*Y(50 )*D -K(
     *         822)*Y(51 )*D -K( 823)*Y(53 )*D -K( 824)*Y(53 )*D -K(
     *         825)*Y(53 )*D -K( 826)*Y(53 )*D -K( 827)*Y(55 )*D -K(
     *         828)*Y(55 )*D -K( 829)*Y(60 )*D -K( 830)*Y(68 )*D -K(
     *         831)*Y(70 )*D -K( 834)*Y(41 )*D -K(1101) -K(1102) -K(
     *        1146)*D ))                                            
      YD( 1) = 
     *        +K( 122)*Y(2  )*Y(40 )*D  +K( 186)*Y(4  )*Y(38 )*D  +K
     *        ( 189)*Y(4  )*Y(40 )*D  +K( 215)*Y(5  )*Y(34 )*D  +K( 
     *        216)*Y(5  )*Y(38 )*D  +K( 258)*Y(7  )*Y(40 )*D  +K( 57
     *        1)*Y(23 )*Y(40 )*D  +K( 620)*Y(25 )*Y(38 )*D  +K( 625)
     *        *Y(25 )*Y(40 )*D  +K( 763)*Y(35 )*Y(40 )*D  +K( 780)*Y
     *        (21 )*Y(38 )*D  +K( 783)*Y(33 )*Y(38 )*D  +K( 787)*Y(3
     *        8 )*Y(46 )*D  +K( 789)*Y(38 )*Y(55 )*D  +K( 790)*Y(38 
     *        )*Y(60 )*D  +K( 796)*Y(38 )*Y(39 )*D  +K( 798)*Y(39 )*
     *        Y(40 )*D  +K( 799)*Y(39 )*Y(45 )*D  +K( 814)*Y(12 )*Y(
     *        40 )*D  +K( 816)*Y(32 )*Y(40 )*D  +K( 821)*Y(40 )*Y(50
     *         )*D  +K( 828)*Y(40 )*Y(55 )*D  +K( 958)*Y(39 )*X(1 )*
     *        D  +K(1102)*Y(40 ) +(Y(41 )*(-K( 306)*Y(8  )*D -K( 392
     *        )*Y(13 )*D -K( 393)*Y(13 )*D -K( 457)*Y(15 )*D -K( 458
     *        )*Y(15 )*D -K( 509)*Y(19 )*D -K( 548)*Y(22 )*D -K( 549
     *        )*Y(22 )*D -K( 605)*Y(24 )*D -K( 648)*Y(26 )*D -K( 649
     *        )*Y(26 )*D -K( 687)*Y(28 )*D -K( 688)*Y(28 )*D -K( 741
     *        )*Y(34 )*D -K( 784)*Y(38 )*D -K( 832)*Y(11 )*D -K( 833
     *        )*Y(36 )*D -K( 834)*Y(40 )*D ))                       
      YDOT( 41) = +YD( 1)
     *        +(Y( 41)*(
     *        -K( 835)*Y(43 )*D -K( 836)*Y(45 )*D -K( 837)*Y(45 )*D 
     *        -K( 838)*Y(45 )*D -K( 839)*Y(47 )*D -K( 840)*Y(54 )*D 
     *        -K( 841)*Y(56 )*D -K( 842)*Y(57 )*D -K( 843)*Y(65 )*D 
     *        -K( 844)*Y(65 )*D -K( 845)*Y(69 )*D -K( 846)*Y(69 )*D 
     *        -K( 847)*Y(69 )*D -K( 959)*X(1 )*D  -K(1006)*X(2 )*D  
     *        -K(1007)*X(2 )*D  -K(1008)*X(2 )*D  -K(1103) -K(1186)*
     *        D ))                                                  
      YDOT( 42) =      +K( 188)*Y(4  )*Y(40 )*D  +K( 217)*Y(5  
     *        )*Y(40 )*D  +K( 422)*Y(14 )*Y(40 )*D  +K( 508)*Y(19 )*
     *        Y(39 )*D  +K( 509)*Y(19 )*Y(41 )*D  +K( 522)*Y(20 )*Y(
     *        40 )*D  +K( 528)*Y(21 )*Y(40 )*D  +K( 605)*Y(24 )*Y(41
     *         )*D  +K( 622)*Y(25 )*Y(40 )*D  +K( 665)*Y(27 )*Y(40 )
     *        *D  +K( 731)*Y(21 )*Y(34 )*D  +K( 784)*Y(38 )*Y(41 )*D
     *          +K( 797)*Y(39 )*Y(40 )*D  +K( 817)*Y(33 )*Y(40 )*D  
     *        +K( 818)*Y(40 )*Y(46 )*D  +K( 819)*Y(40 )*Y(48 )*D  +K
     *        ( 822)*Y(40 )*Y(51 )*D  +K( 827)*Y(40 )*Y(55 )*D  +K( 
     *        829)*Y(40 )*Y(60 )*D  +K( 830)*Y(40 )*Y(68 )*D  +K( 83
     *        1)*Y(40 )*Y(70 )*D  +K( 834)*Y(40 )*Y(41 )*D  +K( 836)
     *        *Y(41 )*Y(45 )*D  +K( 845)*Y(41 )*Y(69 )*D  +K( 959)*Y
     *        (41 )*X(1 )*D  +(Y(42 )*(-K( 161)*Y(3  )*D -K( 162)*Y(
     *        3  )*D -K( 307)*Y(8  )*D -K( 394)*Y(13 )*D -K( 459)*Y(
     *        15 )*D -K( 650)*Y(26 )*D -K( 689)*Y(28 )*D -K( 848)*Y(
     *        43 )*D -K( 849)*Y(54 )*D -K( 850)*Y(56 )*D -K( 851)*Y(
     *        69 )*D -K( 852)*Y(72 )*D -K( 853)*Y(74 )*D -K(1009)*X(
     *        2 )*D  -K(1010)*X(2 )*D  -K(1011)*X(2 )*D  -K(1012)*X(
     *        2 )*D  -K(1187)*D ))                                  
      YDOT( 43) =      +K( 163)*Y(3  )*Y(44 )*D  +K( 355)*Y(10 
     *        )*Y(44 )*D  +K(1013)*Y(44 )*X(2 )*D  +K(1170)*Y(44 )*D
     *         +K(1217)*Y(86 ) +(Y(43 )*(-K( 218)*Y(5  )*D -K( 333)*
     *        Y(9  )*D -K( 425)*Y(14 )*D -K( 494)*Y(18 )*D -K( 709)*
     *        Y(29 )*D -K( 835)*Y(41 )*D -K( 848)*Y(42 )*D -K( 854)*
     *        Y(32 )*D -K( 855)*Y(37 )*D -K( 856)*Y(46 )*D -K( 857)*
     *        Y(48 )*D -K( 858)*Y(58 )*D -K( 859)*Y(61 )*D -K( 860)*
     *        Y(61 )*D -K( 861)*Y(66 )*D -K( 862)*Y(68 )*D -K( 863)*
     *        Y(70 )*D -K( 864)*Y(73 )*D -K(1104) -K(1147)*D ))     
      YDOT( 44) =      +K( 218)*Y(5  )*Y(43 )*D  +K( 333)*Y(9  
     *        )*Y(43 )*D  +K( 425)*Y(14 )*Y(43 )*D  +K( 494)*Y(18 )*
     *        Y(43 )*D  +K( 709)*Y(29 )*Y(43 )*D  +K( 835)*Y(41 )*Y(
     *        43 )*D  +K( 848)*Y(42 )*Y(43 )*D  +K( 854)*Y(32 )*Y(43
     *         )*D  +K( 855)*Y(37 )*Y(43 )*D  +K( 856)*Y(43 )*Y(46 )
     *        *D  +K( 857)*Y(43 )*Y(48 )*D  +K( 858)*Y(43 )*Y(58 )*D
     *          +K( 859)*Y(43 )*Y(61 )*D  +K( 860)*Y(43 )*Y(61 )*D  
     *        +K( 861)*Y(43 )*Y(66 )*D  +K( 862)*Y(43 )*Y(68 )*D  +K
     *        ( 863)*Y(43 )*Y(70 )*D  +K( 864)*Y(43 )*Y(73 )*D  +K(1
     *        104)*Y(43 ) +(Y(44 )*(-K( 163)*Y(3  )*D -K( 355)*Y(10 
     *        )*D -K(1013)*X(2 )*D  -K(1170)*D ))                   
      YD( 1) = 
     *        +K(   7)*Y(1  )*Y(47 )*D  +K(  18)*Y(13 )*Y(36 )*D  +K
     *        (  20)*Y(13 )*Y(47 )*D  +K(  21)*Y(13 )*Y(49 )*D  +K( 
     *         25)*Y(15 )*Y(36 )*D  +K(  30)*Y(15 )*Y(47 )*D  +K(  3
     *        6)*Y(17 )*Y(47 )*D  +K(  58)*Y(34 )*Y(47 )*D  +K(  70)
     *        *Y(47 )*Y(52 )*D  +K( 149)*Y(3  )*Y(11 )*D  +K( 353)*Y
     *        (10 )*Y(38 )*D  +K( 364)*Y(12 )*Y(47 )*D  +K( 396)*Y(1
     *        3 )*Y(48 )*D  +K( 403)*Y(13 )*Y(57 )*D  +K( 420)*Y(14 
     *        )*Y(36 )*D  +K( 450)*Y(15 )*Y(34 )*D  +K( 462)*Y(15 )*
     *        Y(48 )*D  +K( 485)*Y(17 )*Y(36 )*D  +K( 629)*Y(25 )*Y(
     *        49 )*D  +K( 652)*Y(26 )*Y(48 )*D  +K( 669)*Y(27 )*Y(47
     *         )*D  +K( 691)*Y(28 )*Y(48 )*D  +K( 712)*Y(29 )*Y(47 )
     *        *D  +K( 819)*Y(40 )*Y(48 )*D  +K( 856)*Y(43 )*Y(46 )*D
     *          +K( 886)*Y(48 )*Y(65 )*D  +K( 908)*Y(48 )*Y(54 )*D  
     *        +K( 924)*Y(48 )*Y(56 )*D  +K(1016)*Y(48 )*X(2 )*D  +(Y
     *        (45 )*(-K(  19)*Y(13 )*D -K(  42)*Y(22 )*D -K(  69)*Y(
     *        59 )*D -K(  99)*Y(1  )*D -K( 123)*Y(2  )*D -K( 124)*Y(
     *        2  )*D -K( 125)*Y(2  )*D -K( 164)*Y(3  )*D -K( 190)*Y(
     *        4  )*D -K( 191)*Y(4  )*D ))                           
      YDOT( 45) = +YD( 1)
     *        +(Y( 45)*(
     *        -K( 219)*Y(5  )*D -K( 259)*Y(7  )*D -K( 260)*Y(7  )*D 
     *        -K( 308)*Y(8  )*D -K( 334)*Y(9  )*D -K( 335)*Y(9  )*D 
     *        -K( 363)*Y(12 )*D -K( 426)*Y(14 )*D -K( 427)*Y(14 )*D 
     *        -K( 460)*Y(15 )*D -K( 477)*Y(16 )*D -K( 487)*Y(17 )*D 
     *        -K( 495)*Y(18 )*D -K( 496)*Y(18 )*D -K( 529)*Y(21 )*D 
     *        -K( 550)*Y(22 )*D -K( 572)*Y(23 )*D -K( 573)*Y(23 )*D 
     *        -K( 626)*Y(25 )*D -K( 667)*Y(27 )*D -K( 668)*Y(27 )*D 
     *        -K( 710)*Y(29 )*D -K( 711)*Y(29 )*D -K( 717)*Y(32 )*D 
     *        -K( 718)*Y(32 )*D -K( 742)*Y(34 )*D -K( 743)*Y(34 )*D 
     *        -K( 764)*Y(35 )*D -K( 765)*Y(35 )*D -K( 785)*Y(38 )*D 
     *        -K( 799)*Y(39 )*D -K( 800)*Y(39 )*D -K( 801)*Y(39 )*D 
     *        -K( 836)*Y(41 )*D -K( 837)*Y(41 )*D -K( 838)*Y(41 )*D 
     *        -K( 865)*Y(33 )*D -K( 866)*Y(36 )*D -K( 867)*Y(37 )*D 
     *        -2*K( 868)*Y(45 )*D -2*K( 869)*Y(45 )*D -K( 870)*Y(48 
     *        )*D -K( 871)*Y(57 )*D -K( 872)*Y(60 )*D -K( 873)*Y(66 
     *        )*D -K( 874)*Y(66 )*D -K( 875)*Y(70 )*D -K( 876)*Y(46 
     *        )*D -K( 889)*Y(52 )*D -K( 896)*Y(53 )*D -K( 897)*Y(53 
     *        )*D -K( 914)*Y(55 )*D -K( 915)*Y(55 )*D -K(1105) -K(11
     *        06) -K(1148)*D ))                                     
      YD( 1) = 
     *        +K( 100)*Y(1  )*Y(50 )*D  +K( 125)*Y(2  )*Y(45 )*D  +K
     *        ( 127)*Y(2  )*Y(47 )*D  +K( 129)*Y(2  )*Y(49 )*D  +K( 
     *        169)*Y(4  )*Y(11 )*D  +K( 191)*Y(4  )*Y(45 )*D  +K( 19
     *        2)*Y(4  )*Y(47 )*D  +K( 205)*Y(5  )*Y(11 )*D  +K( 263)
     *        *Y(7  )*Y(47 )*D  +K( 307)*Y(8  )*Y(42 )*D  +K( 332)*Y
     *        (9  )*Y(40 )*D  +K( 335)*Y(9  )*Y(45 )*D  +K( 337)*Y(9
     *          )*Y(47 )*D  +K( 359)*Y(11 )*Y(33 )*D  +K( 360)*Y(11 
     *        )*Y(51 )*D  +K( 361)*Y(11 )*Y(60 )*D  +K( 363)*Y(12 )*
     *        Y(45 )*D  +K( 364)*Y(12 )*Y(47 )*D  +K( 369)*Y(12 )*Y(
     *        69 )*D  +K( 371)*Y(12 )*Y(13 )*D  +K( 384)*Y(13 )*Y(34
     *         )*D  +K( 388)*Y(13 )*Y(37 )*D  +K( 419)*Y(14 )*Y(36 )
     *        *D  +K( 424)*Y(14 )*Y(40 )*D  +K( 427)*Y(14 )*Y(45 )*D
     *          +K( 429)*Y(14 )*Y(47 )*D  +K( 430)*Y(14 )*Y(49 )*D  
     *        +K( 439)*Y(12 )*Y(15 )*D  +K( 475)*Y(16 )*Y(34 )*D  +K
     *        ( 476)*Y(16 )*Y(36 )*D  +K( 478)*Y(16 )*Y(47 )*D  +K( 
     *        492)*Y(18 )*Y(34 )*D  +K( 496)*Y(18 )*Y(45 )*D  +K( 49
     *        7)*Y(18 )*Y(47 )*D  +K( 501)*Y(12 )*Y(19 )*D  +K( 517)
     *        *Y(11 )*Y(20 )*D  +K( 527)*Y(11 )*Y(21 )*D            
      YD( 2) = YD( 1)
     *        +K( 573)*Y(23 )*Y(45 )*D  +K( 574)*Y(23 )*Y(47 )*D  +K
     *        ( 587)*Y(12 )*Y(24 )*D  +K( 611)*Y(11 )*Y(25 )*D  +K( 
     *        627)*Y(25 )*Y(47 )*D  +K( 639)*Y(12 )*Y(26 )*D  +K( 66
     *        8)*Y(27 )*Y(45 )*D  +K( 681)*Y(12 )*Y(28 )*D  +K( 711)
     *        *Y(29 )*Y(45 )*D  +K( 718)*Y(32 )*Y(45 )*D  +K( 719)*Y
     *        (32 )*Y(47 )*D  +K( 747)*Y(34 )*Y(51 )*D  +K( 754)*Y(3
     *        4 )*Y(78 )*D  +K( 765)*Y(35 )*Y(45 )*D  +K( 766)*Y(35 
     *        )*Y(47 )*D  +K( 771)*Y(35 )*Y(54 )*D  +K( 778)*Y(12 )*
     *        Y(38 )*D  +K( 791)*Y(11 )*Y(39 )*D  +K( 801)*Y(39 )*Y(
     *        45 )*D  +K( 813)*Y(12 )*Y(40 )*D  +K( 825)*Y(40 )*Y(53
     *         )*D  +K( 832)*Y(11 )*Y(41 )*D  +K( 838)*Y(41 )*Y(45 )
     *        *D  +K( 867)*Y(37 )*Y(45 )*D  +K( 870)*Y(45 )*Y(48 )*D
     *          +K( 874)*Y(45 )*Y(66 )*D  +K( 875)*Y(45 )*Y(70 )*D  
     *        +K( 882)*Y(37 )*Y(47 )*D  +K( 884)*Y(47 )*Y(66 )*D  +K
     *        ( 897)*Y(45 )*Y(53 )*D  +K( 898)*Y(47 )*Y(53 )*D  +K( 
     *        912)*Y(11 )*Y(55 )*D  +K( 949)*Y(12 )*X(1 )*D  +K(1106
     *        )*Y(45 ) +K(1111)*Y(47 ) +(Y(46 )*(-K( 165)*Y(3  )*D -
     *        K( 309)*Y(8  )*D -K( 395)*Y(13 )*D ))                 
      YDOT( 46) = +YD( 2)
     *        +(Y( 46)*(
     *        -K( 461)*Y(15 )*D -K( 606)*Y(24 )*D -K( 651)*Y(26 )*D 
     *        -K( 690)*Y(28 )*D -K( 786)*Y(38 )*D -K( 787)*Y(38 )*D 
     *        -K( 818)*Y(40 )*D -K( 856)*Y(43 )*D -K( 876)*Y(45 )*D 
     *        -K( 877)*Y(65 )*D -K( 878)*Y(67 )*D -K( 879)*Y(69 )*D 
     *        -K( 880)*Y(72 )*D -K( 881)*Y(74 )*D -K( 907)*Y(54 )*D 
     *        -K( 923)*Y(56 )*D -K(1014)*X(2 )*D  -K(1107) -K(1188)*
     *        D ))                                                  
      YD( 1) = 
     *        +K(  26)*Y(15 )*Y(36 )*D  +K(  32)*Y(15 )*Y(57 )*D  +K
     *        (  69)*Y(45 )*Y(59 )*D  +K( 164)*Y(3  )*Y(45 )*D  +K( 
     *        354)*Y(10 )*Y(40 )*D  +K( 397)*Y(13 )*Y(48 )*D  +K( 45
     *        4)*Y(15 )*Y(38 )*D  +K( 463)*Y(15 )*Y(48 )*D  +K( 484)
     *        *Y(17 )*Y(34 )*D  +K( 486)*Y(17 )*Y(38 )*D  +K( 692)*Y
     *        (28 )*Y(48 )*D  +K( 857)*Y(43 )*Y(48 )*D  +K( 869)*Y(4
     *        5 )*Y(45 )*D  +K( 870)*Y(45 )*Y(48 )*D  +K( 887)*Y(48 
     *        )*Y(65 )*D  +K( 929)*Y(48 )*Y(57 )*D  +K(1017)*Y(48 )*
     *        X(2 )*D  +K(1218)*Y(87 ) +(Y(47 )*(-K(   7)*Y(1  )*D -
     *        K(  20)*Y(13 )*D -K(  30)*Y(15 )*D -K(  36)*Y(17 )*D -
     *        K(  58)*Y(34 )*D -K(  70)*Y(52 )*D -K( 126)*Y(2  )*D -
     *        K( 127)*Y(2  )*D -K( 128)*Y(2  )*D -K( 192)*Y(4  )*D -
     *        K( 193)*Y(4  )*D -K( 261)*Y(7  )*D -K( 262)*Y(7  )*D -
     *        K( 263)*Y(7  )*D -K( 264)*Y(7  )*D -K( 336)*Y(9  )*D -
     *        K( 337)*Y(9  )*D -K( 338)*Y(9  )*D -K( 364)*Y(12 )*D -
     *        K( 365)*Y(12 )*D -K( 428)*Y(14 )*D -K( 429)*Y(14 )*D -
     *        K( 478)*Y(16 )*D -K( 497)*Y(18 )*D -K( 523)*Y(20 )*D -
     *        K( 574)*Y(23 )*D -K( 575)*Y(23 )*D ))                 
      YDOT( 47) = +YD( 1)
     *        +(Y( 47)*(
     *        -K( 576)*Y(23 )*D -K( 627)*Y(25 )*D -K( 628)*Y(25 )*D 
     *        -K( 669)*Y(27 )*D -K( 712)*Y(29 )*D -K( 719)*Y(32 )*D 
     *        -K( 720)*Y(32 )*D -K( 744)*Y(34 )*D -K( 766)*Y(35 )*D 
     *        -K( 767)*Y(35 )*D -K( 802)*Y(39 )*D -K( 839)*Y(41 )*D 
     *        -K( 882)*Y(37 )*D -K( 883)*Y(37 )*D -K( 884)*Y(66 )*D 
     *        -K( 885)*Y(66 )*D -K( 898)*Y(53 )*D -K( 899)*Y(53 )*D 
     *        -K(1108) -K(1109) -K(1110) -K(1111) -K(1149)*D ))     
      YD( 1) = 
     *        +K( 128)*Y(2  )*Y(47 )*D  +K( 193)*Y(4  )*Y(47 )*D  +K
     *        ( 219)*Y(5  )*Y(45 )*D  +K( 264)*Y(7  )*Y(47 )*D  +K( 
     *        338)*Y(9  )*Y(47 )*D  +K( 365)*Y(12 )*Y(47 )*D  +K( 42
     *        3)*Y(14 )*Y(40 )*D  +K( 452)*Y(15 )*Y(37 )*D  +K( 479)
     *        *Y(16 )*Y(49 )*D  +K( 491)*Y(18 )*Y(34 )*D  +K( 493)*Y
     *        (18 )*Y(38 )*D  +K( 523)*Y(20 )*Y(47 )*D  +K( 529)*Y(2
     *        1 )*Y(45 )*D  +K( 576)*Y(23 )*Y(47 )*D  +K( 626)*Y(25 
     *        )*Y(45 )*D  +K( 628)*Y(25 )*Y(47 )*D  +K( 667)*Y(27 )*
     *        Y(45 )*D  +K( 720)*Y(32 )*Y(47 )*D  +K( 767)*Y(35 )*Y(
     *        47 )*D  +K( 800)*Y(39 )*Y(45 )*D  +K( 802)*Y(39 )*Y(47
     *         )*D  +K( 837)*Y(41 )*Y(45 )*D  +K( 839)*Y(41 )*Y(47 )
     *        *D  +K( 865)*Y(33 )*Y(45 )*D  +K( 872)*Y(45 )*Y(60 )*D
     *          +K( 876)*Y(45 )*Y(46 )*D  +K( 883)*Y(37 )*Y(47 )*D  
     *        +K( 899)*Y(47 )*Y(53 )*D  +K( 915)*Y(45 )*Y(55 )*D  +K
     *        (1108)*Y(47 ) +(Y(48 )*(-K( 396)*Y(13 )*D -K( 397)*Y(1
     *        3 )*D -K( 462)*Y(15 )*D -K( 463)*Y(15 )*D -K( 652)*Y(2
     *        6 )*D -K( 691)*Y(28 )*D -K( 692)*Y(28 )*D -K( 819)*Y(4
     *        0 )*D -K( 857)*Y(43 )*D -K( 870)*Y(45 )*D ))          
      YDOT( 48) = +YD( 1)
     *        +(Y( 48)*(
     *        -K( 886)*Y(65 )*D -K( 887)*Y(65 )*D -K( 908)*Y(54 )*D 
     *        -K( 924)*Y(56 )*D -K( 929)*Y(57 )*D -K(1015)*X(2 )*D  
     *        -K(1016)*X(2 )*D  -K(1017)*X(2 )*D  -K(1189)*D ))     
      YD( 1) = 
     *        +K(  23)*Y(15 )*Y(36 )*D  +K(  24)*Y(15 )*Y(36 )*D  +K
     *        (  64)*Y(11 )*Y(38 )*D  +K( 310)*Y(8  )*Y(51 )*D  +K( 
     *        352)*Y(10 )*Y(36 )*D  +K( 360)*Y(11 )*Y(51 )*D  +K( 51
     *        1)*Y(19 )*Y(50 )*D  +K( 512)*Y(19 )*Y(51 )*D  +K( 693)
     *        *Y(28 )*Y(50 )*D  +K( 694)*Y(28 )*Y(51 )*D  +K( 742)*Y
     *        (34 )*Y(45 )*D  +K( 746)*Y(34 )*Y(50 )*D  +K( 775)*Y(3
     *        6 )*Y(50 )*D  +K( 821)*Y(40 )*Y(50 )*D  +K( 822)*Y(40 
     *        )*Y(51 )*D  +K( 866)*Y(36 )*Y(45 )*D  +K( 930)*Y(50 )*
     *        Y(57 )*D  +K( 943)*Y(50 )*Y(69 )*D  +K(1020)*Y(51 )*X(
     *        2 )*D  +K(1219)*Y(88 ) +(Y(49 )*(-K(  21)*Y(13 )*D -K(
     *          43)*Y(22 )*D -K( 129)*Y(2  )*D -K( 194)*Y(4  )*D -K(
     *         220)*Y(5  )*D -K( 265)*Y(7  )*D -K( 266)*Y(7  )*D -K(
     *         267)*Y(7  )*D -K( 268)*Y(7  )*D -K( 269)*Y(7  )*D -K(
     *         339)*Y(9  )*D -K( 356)*Y(10 )*D -K( 366)*Y(12 )*D -K(
     *         430)*Y(14 )*D -K( 479)*Y(16 )*D -K( 524)*Y(20 )*D -K(
     *         530)*Y(21 )*D -K( 577)*Y(23 )*D -K( 578)*Y(23 )*D -K(
     *         629)*Y(25 )*D -K( 630)*Y(25 )*D -K( 631)*Y(25 )*D -K(
     *         721)*Y(32 )*D -K( 727)*Y(33 )*D ))                   
      YDOT( 49) = +YD( 1)
     *        +(Y( 49)*(
     *        -K( 768)*Y(35 )*D -K( 803)*Y(39 )*D -K( 900)*Y(53 )*D 
     *        -K( 916)*Y(55 )*D -K( 935)*Y(60 )*D -K(1112) -K(1150)*
     *        D ))                                                  
      YDOT( 50) =      +K( 269)*Y(7  )*Y(49 )*D  +K( 366)*Y(12 
     *        )*Y(49 )*D  +K( 578)*Y(23 )*Y(49 )*D  +K( 721)*Y(32 )*
     *        Y(49 )*D  +K( 900)*Y(49 )*Y(53 )*D  +(Y(50 )*(-K( 100)
     *        *Y(1  )*D -K( 510)*Y(19 )*D -K( 511)*Y(19 )*D -K( 693)
     *        *Y(28 )*D -K( 745)*Y(34 )*D -K( 746)*Y(34 )*D -K( 775)
     *        *Y(36 )*D -K( 820)*Y(40 )*D -K( 821)*Y(40 )*D -K( 930)
     *        *Y(57 )*D -K( 943)*Y(69 )*D -K( 960)*X(1 )*D  -K(1018)
     *        *X(2 )*D  -K(1190)*D ))                               
      YDOT( 51) =      +K( 194)*Y(4  )*Y(49 )*D  +K( 220)*Y(5  
     *        )*Y(49 )*D  +K( 510)*Y(19 )*Y(50 )*D  +K( 524)*Y(20 )*
     *        Y(49 )*D  +K( 530)*Y(21 )*Y(49 )*D  +K( 631)*Y(25 )*Y(
     *        49 )*D  +K( 727)*Y(33 )*Y(49 )*D  +K( 786)*Y(38 )*Y(46
     *         )*D  +K( 803)*Y(39 )*Y(49 )*D  +K( 820)*Y(40 )*Y(50 )
     *        *D  +K( 916)*Y(49 )*Y(55 )*D  +K( 935)*Y(49 )*Y(60 )*D
     *          +K( 960)*Y(50 )*X(1 )*D  +(Y(51 )*(-K( 310)*Y(8  )*D
     *         -K( 360)*Y(11 )*D -K( 512)*Y(19 )*D -K( 694)*Y(28 )*D
     *         -K( 747)*Y(34 )*D -K( 822)*Y(40 )*D -K(1019)*X(2 )*D 
     *         -K(1020)*X(2 )*D  -K(1021)*X(2 )*D  -K(1191)*D ))    
      YD( 1) = 
     *        +K(  11)*Y(8  )*Y(22 )*D  +K(  46)*Y(22 )*Y(72 )*D  +K
     *        ( 101)*Y(1  )*Y(53 )*D  +K( 295)*Y(8  )*Y(24 )*D  +K( 
     *        311)*Y(8  )*Y(53 )*D  +K( 312)*Y(8  )*Y(55 )*D  +K( 31
     *        3)*Y(8  )*Y(57 )*D  +K( 347)*Y(10 )*Y(22 )*D  +K( 374)
     *        *Y(13 )*Y(22 )*D  +K( 398)*Y(13 )*Y(53 )*D  +K( 399)*Y
     *        (13 )*Y(55 )*D  +K( 402)*Y(13 )*Y(57 )*D  +K( 464)*Y(1
     *        5 )*Y(53 )*D  +K( 465)*Y(15 )*Y(55 )*D  +K( 607)*Y(24 
     *        )*Y(53 )*D  +K( 608)*Y(24 )*Y(55 )*D  +K( 653)*Y(26 )*
     *        Y(53 )*D  +K( 654)*Y(26 )*Y(55 )*D  +K( 748)*Y(34 )*Y(
     *        53 )*D  +K( 788)*Y(38 )*Y(53 )*D  +K( 789)*Y(38 )*Y(55
     *         )*D  +K( 827)*Y(40 )*Y(55 )*D  +K( 893)*Y(11 )*Y(53 )
     *        *D  +K( 895)*Y(36 )*Y(53 )*D  +K( 897)*Y(45 )*Y(53 )*D
     *          +K( 899)*Y(47 )*Y(53 )*D  +K( 900)*Y(49 )*Y(53 )*D  
     *        +K( 901)*Y(53 )*Y(54 )*D  +K( 902)*Y(53 )*Y(57 )*D  +K
     *        ( 903)*Y(53 )*Y(65 )*D  +K( 912)*Y(11 )*Y(55 )*D  +K( 
     *        915)*Y(45 )*Y(55 )*D  +K( 916)*Y(49 )*Y(55 )*D  +K( 91
     *        7)*Y(54 )*Y(55 )*D  +K( 918)*Y(55 )*Y(56 )*D  +K( 920)
     *        *Y(55 )*Y(65 )*D  +K(1023)*Y(55 )*X(2 )*D             
      YDOT( 52) = +YD( 1)
     *        +K(1026)*Y(61 )*X(2 )*D  +K(1029)*Y(62 )*X(2 )*D  +K(1
     *        114)*Y(54 ) +K(1115)*Y(56 ) +K(1193)*Y(55 )*D +K(1196)
     *        *Y(61 )*D +(Y(52 )*(-K(  31)*Y(15 )*D -K(  37)*Y(17 )*
     *        D -K(  40)*Y(19 )*D -K(  50)*Y(24 )*D -K(  59)*Y(34 )*
     *        D -K(  67)*Y(38 )*D -K(  70)*Y(47 )*D -K(  79)*X(1 )*D
     *          -K(  85)*Y(28 )*D -K(  87)*Y(36 )*D -K( 166)*Y(3  )*
     *        D -K( 195)*Y(4  )*D -K( 196)*Y(4  )*D -K( 221)*Y(5  )*
     *        D -K( 270)*Y(7  )*D -K( 271)*Y(7  )*D -K( 551)*Y(22 )*
     *        D -K( 579)*Y(23 )*D -K( 632)*Y(25 )*D -K( 670)*Y(27 )*
     *        D -K( 769)*Y(35 )*D -K( 804)*Y(39 )*D -K( 888)*Y(32 )*
     *        D -K( 889)*Y(45 )*D -K( 890)*Y(57 )*D -K( 891)*Y(59 )*
     *        D -K( 892)*Y(60 )*D -K(1113) -K(1151)*D ))            
      YDOT( 53) =      +K(  14)*Y(9  )*Y(22 )*D  +K( 196)*Y(4  
     *        )*Y(52 )*D  +K( 275)*Y(7  )*Y(54 )*D  +K( 278)*Y(7  )*
     *        Y(56 )*D  +K( 323)*Y(9  )*Y(24 )*D  +K( 375)*Y(13 )*Y(
     *        23 )*D  +K( 411)*Y(14 )*Y(22 )*D  +K( 412)*Y(14 )*Y(24
     *         )*D  +K( 579)*Y(23 )*Y(52 )*D  +K( 888)*Y(32 )*Y(52 )
     *        *D  +(Y(53 )*(-K( 101)*Y(1  )*D -K( 311)*Y(8  )*D -K( 
     *        398)*Y(13 )*D -K( 464)*Y(15 )*D -K( 552)*Y(22 )*D -K( 
     *        607)*Y(24 )*D -K( 653)*Y(26 )*D -K( 748)*Y(34 )*D -K( 
     *        788)*Y(38 )*D -K( 823)*Y(40 )*D -K( 824)*Y(40 )*D -K( 
     *        825)*Y(40 )*D -K( 826)*Y(40 )*D -K( 893)*Y(11 )*D -K( 
     *        894)*Y(36 )*D -K( 895)*Y(36 )*D -K( 896)*Y(45 )*D -K( 
     *        897)*Y(45 )*D -K( 898)*Y(47 )*D -K( 899)*Y(47 )*D -K( 
     *        900)*Y(49 )*D -K( 901)*Y(54 )*D -K( 902)*Y(57 )*D -K( 
     *        903)*Y(65 )*D -K( 961)*X(1 )*D  -K(1022)*X(2 )*D  -K(1
     *        192)*D ))                                             
      YD( 1) = 
     *        +K(  31)*Y(15 )*Y(52 )*D  +K(  37)*Y(17 )*Y(52 )*D  +K
     *        (  40)*Y(19 )*Y(52 )*D  +K(  50)*Y(24 )*Y(52 )*D  +K( 
     *         67)*Y(38 )*Y(52 )*D  +K(  70)*Y(47 )*Y(52 )*D  +K(  7
     *        9)*Y(52 )*X(1 )*D  +K(  82)*Y(8  )*Y(26 )*D  +K(  85)*
     *        Y(28 )*Y(52 )*D  +K( 102)*Y(1  )*Y(55 )*D  +K( 103)*Y(
     *        1  )*Y(56 )*D  +K( 131)*Y(2  )*Y(56 )*D  +K( 166)*Y(3 
     *         )*Y(52 )*D  +K( 349)*Y(10 )*Y(24 )*D  +K( 401)*Y(13 )
     *        *Y(57 )*D  +K( 406)*Y(13 )*Y(61 )*D  +K( 466)*Y(15 )*Y
     *        (57 )*D  +K( 469)*Y(15 )*Y(61 )*D  +K( 535)*Y(15 )*Y(2
     *        2 )*D  +K( 538)*Y(17 )*Y(22 )*D  +K( 539)*Y(17 )*Y(22 
     *        )*D  +K( 550)*Y(22 )*Y(45 )*D  +K( 553)*Y(22 )*Y(77 )*
     *        D  +K( 657)*Y(26 )*Y(61 )*D  +K( 696)*Y(28 )*Y(55 )*D 
     *         +K( 698)*Y(28 )*Y(61 )*D  +K( 828)*Y(40 )*Y(55 )*D  +
     *        K( 859)*Y(43 )*Y(61 )*D  +K( 889)*Y(45 )*Y(52 )*D  +K(
     *         891)*Y(52 )*Y(59 )*D  +K( 898)*Y(47 )*Y(53 )*D  +K( 9
     *        13)*Y(36 )*Y(55 )*D  +K( 919)*Y(55 )*Y(57 )*D  +K( 921
     *        )*Y(55 )*Y(65 )*D  +K( 937)*Y(61 )*Y(69 )*D  +K(1027)*
     *        Y(61 )*X(2 )*D  +K(1220)*Y(89 )                       
      YDOT( 54) = +YD( 1)
     *        +(Y(54 )*(-K( 130)*Y(2  )*D -K( 197)*Y(4  )*D -K( 222)
     *        *Y(5  )*D -K( 272)*Y(7  )*D -K( 273)*Y(7  )*D -K( 274)
     *        *Y(7  )*D -K( 275)*Y(7  )*D -K( 431)*Y(14 )*D -K( 531)
     *        *Y(21 )*D -K( 580)*Y(23 )*D -K( 633)*Y(25 )*D -K( 671)
     *        *Y(27 )*D -K( 770)*Y(35 )*D -K( 771)*Y(35 )*D -K( 805)
     *        *Y(39 )*D -K( 840)*Y(41 )*D -K( 849)*Y(42 )*D -K( 901)
     *        *Y(53 )*D -K( 904)*Y(12 )*D -K( 905)*Y(32 )*D -K( 906)
     *        *Y(33 )*D -K( 907)*Y(46 )*D -K( 908)*Y(48 )*D -K( 909)
     *        *Y(60 )*D -K( 910)*Y(68 )*D -K( 911)*Y(71 )*D -K( 917)
     *        *Y(55 )*D -K(1114) -K(1152)*D ))                      
      YDOT( 55) =      +K( 130)*Y(2  )*Y(54 )*D  +K( 195)*Y(4  
     *        )*Y(52 )*D  +K( 197)*Y(4  )*Y(54 )*D  +K( 221)*Y(5  )*
     *        Y(52 )*D  +K( 324)*Y(9  )*Y(26 )*D  +K( 325)*Y(9  )*Y(
     *        28 )*D  +K( 413)*Y(14 )*Y(26 )*D  +K( 537)*Y(16 )*Y(22
     *         )*D  +K( 558)*Y(19 )*Y(23 )*D  +K( 580)*Y(23 )*Y(54 )
     *        *D  +K( 632)*Y(25 )*Y(52 )*D  +K( 804)*Y(39 )*Y(52 )*D
     *          +K( 823)*Y(40 )*Y(53 )*D  +K( 892)*Y(52 )*Y(60 )*D  
     *        +K( 896)*Y(45 )*Y(53 )*D  +K( 901)*Y(53 )*Y(54 )*D  +K
     *        ( 904)*Y(12 )*Y(54 )*D  +K( 905)*Y(32 )*Y(54 )*D  +K( 
     *        961)*Y(53 )*X(1 )*D  +(Y(55 )*(-K( 102)*Y(1  )*D -K( 3
     *        12)*Y(8  )*D -K( 399)*Y(13 )*D -K( 465)*Y(15 )*D -K( 5
     *        13)*Y(19 )*D -K( 608)*Y(24 )*D -K( 654)*Y(26 )*D -K( 6
     *        95)*Y(28 )*D -K( 696)*Y(28 )*D -K( 789)*Y(38 )*D -K( 8
     *        27)*Y(40 )*D -K( 828)*Y(40 )*D -K( 912)*Y(11 )*D -K( 9
     *        13)*Y(36 )*D -K( 914)*Y(45 )*D -K( 915)*Y(45 )*D -K( 9
     *        16)*Y(49 )*D -K( 917)*Y(54 )*D -K( 918)*Y(56 )*D -K( 9
     *        19)*Y(57 )*D -K( 920)*Y(65 )*D -K( 921)*Y(65 )*D -K( 9
     *        62)*X(1 )*D  -K(1023)*X(2 )*D  -K(1193)*D ))          
      YDOT( 56) =      +K(  83)*Y(8  )*Y(26 )*D  +K( 407)*Y(13 
     *        )*Y(61 )*D  +K( 470)*Y(15 )*Y(61 )*D  +K( 536)*Y(15 )*
     *        Y(22 )*D  +K( 658)*Y(26 )*Y(61 )*D  +K( 699)*Y(28 )*Y(
     *        61 )*D  +K( 860)*Y(43 )*Y(61 )*D  +K( 938)*Y(61 )*Y(69
     *         )*D  +K(1028)*Y(61 )*X(2 )*D  +K(1030)*Y(62 )*X(2 )*D
     *          +K(1197)*Y(62 )*D +K(1221)*Y(90 ) +(Y(56 )*(-K( 103)
     *        *Y(1  )*D -K( 131)*Y(2  )*D -K( 223)*Y(5  )*D -K( 276)
     *        *Y(7  )*D -K( 277)*Y(7  )*D -K( 278)*Y(7  )*D -K( 432)
     *        *Y(14 )*D -K( 532)*Y(21 )*D -K( 634)*Y(25 )*D -K( 672)
     *        *Y(27 )*D -K( 806)*Y(39 )*D -K( 841)*Y(41 )*D -K( 850)
     *        *Y(42 )*D -K( 918)*Y(55 )*D -K( 922)*Y(33 )*D -K( 923)
     *        *Y(46 )*D -K( 924)*Y(48 )*D -K( 925)*Y(60 )*D -K( 926)
     *        *Y(68 )*D -K( 927)*Y(71 )*D -K(1115) -K(1153)*D ))    
      YD( 1) = 
     *        +K(   8)*Y(1  )*Y(59 )*D  +K(  41)*Y(22 )*Y(36 )*D  +K
     *        (  43)*Y(22 )*Y(49 )*D  +K(  45)*Y(22 )*Y(59 )*D  +K( 
     *         48)*Y(24 )*Y(36 )*D  +K(  69)*Y(45 )*Y(59 )*D  +K(  8
     *        7)*Y(36 )*Y(52 )*D  +K( 282)*Y(7  )*Y(59 )*D  +K( 315)
     *        *Y(8  )*Y(60 )*D  +K( 361)*Y(11 )*Y(60 )*D  +K( 404)*Y
     *        (13 )*Y(59 )*D  +K( 405)*Y(13 )*Y(60 )*D  +K( 467)*Y(1
     *        5 )*Y(59 )*D  +K( 468)*Y(15 )*Y(60 )*D  +K( 488)*Y(17 
     *        )*Y(59 )*D  +K( 514)*Y(19 )*Y(60 )*D  +K( 546)*Y(22 )*
     *        Y(38 )*D  +K( 567)*Y(23 )*Y(36 )*D  +K( 577)*Y(23 )*Y(
     *        49 )*D  +K( 598)*Y(24 )*Y(34 )*D  +K( 610)*Y(24 )*Y(60
     *         )*D  +K( 656)*Y(26 )*Y(60 )*D  +K( 697)*Y(28 )*Y(60 )
     *        *D  +K( 715)*Y(31 )*Y(60 )*D  +K( 732)*Y(26 )*Y(34 )*D
     *          +K( 749)*Y(34 )*Y(59 )*D  +K( 790)*Y(38 )*Y(60 )*D  
     *        +K( 829)*Y(40 )*Y(60 )*D  +K( 858)*Y(43 )*Y(58 )*D  +K
     *        ( 872)*Y(45 )*Y(60 )*D  +K( 891)*Y(52 )*Y(59 )*D  +K( 
     *        892)*Y(52 )*Y(60 )*D  +K( 909)*Y(54 )*Y(60 )*D  +K( 92
     *        5)*Y(56 )*Y(60 )*D  +K( 935)*Y(49 )*Y(60 )*D  +K( 936)
     *        *Y(60 )*Y(65 )*D  +K(1025)*Y(60 )*X(2 )*D             
      YDOT( 57) = +YD( 1)
     *        +K(1032)*Y(64 )*X(2 )*D  +K(1118)*Y(59 ) +K(1195)*Y(60
     *         )*D +(Y(57 )*(-K(  32)*Y(15 )*D -K(  44)*Y(22 )*D -K(
     *          51)*Y(24 )*D -K(  54)*Y(26 )*D -K( 132)*Y(2  )*D -K(
     *         198)*Y(4  )*D -K( 199)*Y(4  )*D -K( 224)*Y(5  )*D -K(
     *         279)*Y(7  )*D -K( 280)*Y(7  )*D -K( 313)*Y(8  )*D -K(
     *         314)*Y(8  )*D -K( 340)*Y(9  )*D -K( 367)*Y(12 )*D -K(
     *         400)*Y(13 )*D -K( 401)*Y(13 )*D -K( 402)*Y(13 )*D -K(
     *         403)*Y(13 )*D -K( 433)*Y(14 )*D -K( 466)*Y(15 )*D -K(
     *         480)*Y(16 )*D -K( 498)*Y(18 )*D -K( 581)*Y(23 )*D -K(
     *         582)*Y(23 )*D -K( 609)*Y(24 )*D -K( 635)*Y(25 )*D -K(
     *         636)*Y(25 )*D -K( 655)*Y(26 )*D -K( 673)*Y(27 )*D -K(
     *         713)*Y(29 )*D -K( 722)*Y(32 )*D -K( 807)*Y(39 )*D -K(
     *         808)*Y(39 )*D -K( 842)*Y(41 )*D -K( 871)*Y(45 )*D -K(
     *         890)*Y(52 )*D -K( 902)*Y(53 )*D -K( 919)*Y(55 )*D -K(
     *         928)*Y(37 )*D -K( 929)*Y(48 )*D -K( 930)*Y(50 )*D -K(
     *         931)*Y(60 )*D -K( 932)*Y(66 )*D -K( 933)*Y(68 )*D -K(
     *         934)*Y(70 )*D -K(1116) -K(1117) -K(1154)*D ))        
      YD( 1) = 
     *        +K( 132)*Y(2  )*Y(57 )*D  +K( 133)*Y(2  )*Y(59 )*D  +K
     *        ( 199)*Y(4  )*Y(57 )*D  +K( 281)*Y(7  )*Y(59 )*D  +K( 
     *        340)*Y(9  )*Y(57 )*D  +K( 367)*Y(12 )*Y(57 )*D  +K( 43
     *        3)*Y(14 )*Y(57 )*D  +K( 480)*Y(16 )*Y(57 )*D  +K( 498)
     *        *Y(18 )*Y(57 )*D  +K( 545)*Y(22 )*Y(37 )*D  +K( 547)*Y
     *        (22 )*Y(39 )*D  +K( 548)*Y(22 )*Y(41 )*D  +K( 554)*Y(1
     *        1 )*Y(23 )*D  +K( 568)*Y(23 )*Y(36 )*D  +K( 575)*Y(23 
     *        )*Y(47 )*D  +K( 582)*Y(23 )*Y(57 )*D  +K( 599)*Y(24 )*
     *        Y(35 )*D  +K( 618)*Y(25 )*Y(36 )*D  +K( 629)*Y(25 )*Y(
     *        49 )*D  +K( 636)*Y(25 )*Y(57 )*D  +K( 673)*Y(27 )*Y(57
     *         )*D  +K( 713)*Y(29 )*Y(57 )*D  +K( 722)*Y(32 )*Y(57 )
     *        *D  +K( 735)*Y(32 )*Y(34 )*D  +K( 759)*Y(31 )*Y(35 )*D
     *          +K( 769)*Y(35 )*Y(52 )*D  +K( 770)*Y(35 )*Y(54 )*D  
     *        +K( 808)*Y(39 )*Y(57 )*D  +K( 842)*Y(41 )*Y(57 )*D  +K
     *        ( 894)*Y(36 )*Y(53 )*D  +K( 902)*Y(53 )*Y(57 )*D  +K( 
     *        919)*Y(55 )*Y(57 )*D  +K( 928)*Y(37 )*Y(57 )*D  +K( 92
     *        9)*Y(48 )*Y(57 )*D  +K( 930)*Y(50 )*Y(57 )*D  +K( 931)
     *        *Y(57 )*Y(60 )*D  +K( 932)*Y(57 )*Y(66 )*D            
      YDOT( 58) = +YD( 1)
     *        +K( 933)*Y(57 )*Y(68 )*D  +K( 934)*Y(57 )*Y(70 )*D  +K
     *        (1117)*Y(57 ) +(Y(58 )*(-K( 858)*Y(43 )*D -K(1024)*X(2
     *         )*D  -K(1194)*D ))                                   
      YDOT( 59) =      +K(  49)*Y(24 )*Y(36 )*D  +K(  86)*Y(26 
     *        )*Y(34 )*D  +K( 602)*Y(24 )*Y(38 )*D  +K( 871)*Y(45 )*
     *        Y(57 )*D  +K( 931)*Y(57 )*Y(60 )*D  +K(1033)*Y(64 )*X(
     *        2 )*D  +K(1199)*Y(64 )*D +K(1222)*Y(91 ) +(Y(59 )*(-K(
     *           8)*Y(1  )*D -K(  45)*Y(22 )*D -K(  60)*Y(34 )*D -K(
     *          69)*Y(45 )*D -K( 133)*Y(2  )*D -K( 225)*Y(5  )*D -K(
     *         281)*Y(7  )*D -K( 282)*Y(7  )*D -K( 404)*Y(13 )*D -K(
     *         467)*Y(15 )*D -K( 488)*Y(17 )*D -K( 749)*Y(34 )*D -K(
     *         891)*Y(52 )*D -K(1118) -K(1155)*D ))                 
      YDOT( 60) =      +K( 198)*Y(4  )*Y(57 )*D  +K( 224)*Y(5  
     *        )*Y(57 )*D  +K( 549)*Y(22 )*Y(41 )*D  +K( 601)*Y(24 )*
     *        Y(37 )*D  +K( 623)*Y(25 )*Y(40 )*D  +K( 630)*Y(25 )*Y(
     *        49 )*D  +K( 662)*Y(27 )*Y(36 )*D  +K( 733)*Y(27 )*Y(34
     *         )*D  +K( 734)*Y(29 )*Y(34 )*D  +K( 807)*Y(39 )*Y(57 )
     *        *D  +(Y(60 )*(-K( 315)*Y(8  )*D -K( 361)*Y(11 )*D -K( 
     *        405)*Y(13 )*D -K( 468)*Y(15 )*D -K( 514)*Y(19 )*D -K( 
     *        610)*Y(24 )*D -K( 656)*Y(26 )*D -K( 697)*Y(28 )*D -K( 
     *        715)*Y(31 )*D -K( 790)*Y(38 )*D -K( 829)*Y(40 )*D -K( 
     *        872)*Y(45 )*D -K( 892)*Y(52 )*D -K( 909)*Y(54 )*D -K( 
     *        925)*Y(56 )*D -K( 931)*Y(57 )*D -K( 935)*Y(49 )*D -K( 
     *        936)*Y(65 )*D -K(1025)*X(2 )*D  -K(1195)*D ))         
      YD( 1) = 
     *        +K( 222)*Y(5  )*Y(54 )*D  +K( 223)*Y(5  )*Y(56 )*D  +K
     *        ( 431)*Y(14 )*Y(54 )*D  +K( 432)*Y(14 )*Y(56 )*D  +K( 
     *        513)*Y(19 )*Y(55 )*D  +K( 531)*Y(21 )*Y(54 )*D  +K( 53
     *        2)*Y(21 )*Y(56 )*D  +K( 559)*Y(19 )*Y(23 )*D  +K( 589)
     *        *Y(18 )*Y(24 )*D  +K( 633)*Y(25 )*Y(54 )*D  +K( 634)*Y
     *        (25 )*Y(56 )*D  +K( 671)*Y(27 )*Y(54 )*D  +K( 672)*Y(2
     *        7 )*Y(56 )*D  +K( 695)*Y(28 )*Y(55 )*D  +K( 805)*Y(39 
     *        )*Y(54 )*D  +K( 806)*Y(39 )*Y(56 )*D  +K( 840)*Y(41 )*
     *        Y(54 )*D  +K( 841)*Y(41 )*Y(56 )*D  +K( 849)*Y(42 )*Y(
     *        54 )*D  +K( 850)*Y(42 )*Y(56 )*D  +K( 906)*Y(33 )*Y(54
     *         )*D  +K( 907)*Y(46 )*Y(54 )*D  +K( 908)*Y(48 )*Y(54 )
     *        *D  +K( 909)*Y(54 )*Y(60 )*D  +K( 910)*Y(54 )*Y(68 )*D
     *          +K( 911)*Y(54 )*Y(71 )*D  +K( 914)*Y(45 )*Y(55 )*D  
     *        +K( 917)*Y(54 )*Y(55 )*D  +K( 918)*Y(55 )*Y(56 )*D  +K
     *        ( 922)*Y(33 )*Y(56 )*D  +K( 923)*Y(46 )*Y(56 )*D  +K( 
     *        924)*Y(48 )*Y(56 )*D  +K( 925)*Y(56 )*Y(60 )*D  +K( 92
     *        6)*Y(56 )*Y(68 )*D  +K( 927)*Y(56 )*Y(71 )*D  +K( 962)
     *        *Y(55 )*X(1 )*D                                       
      YDOT( 61) = +YD( 1)
     *        +(Y(61 )*(-K( 406)*Y(13 )*D -K( 407)*Y(13 )*D -K( 469)
     *        *Y(15 )*D -K( 470)*Y(15 )*D -K( 657)*Y(26 )*D -K( 658)
     *        *Y(26 )*D -K( 698)*Y(28 )*D -K( 699)*Y(28 )*D -K( 859)
     *        *Y(43 )*D -K( 860)*Y(43 )*D -K( 937)*Y(69 )*D -K( 938)
     *        *Y(69 )*D -K(1026)*X(2 )*D  -K(1027)*X(2 )*D  -K(1028)
     *        *X(2 )*D  -K(1196)*D ))                               
      YDOT( 62) =      +K( 297)*Y(8  )*Y(30 )*D  +K( 326)*Y(9  
     *        )*Y(28 )*D  +K( 415)*Y(14 )*Y(28 )*D  +K( 540)*Y(18 )*
     *        Y(22 )*D  +K( 670)*Y(27 )*Y(52 )*D  +K( 824)*Y(40 )*Y(
     *        53 )*D  +(Y(62 )*(-K(1029)*X(2 )*D  -K(1030)*X(2 )*D  
     *        -K(1197)*D ))                                         
      YDOT( 63) =      +K( 826)*Y(40 )*Y(53 )*D  +(Y(63 )*(-K(1
     *        031)*X(2 )*D  -K(1198)*D ))                           
      YDOT( 64) =      +K( 225)*Y(5  )*Y(59 )*D  +K( 663)*Y(27 
     *        )*Y(36 )*D  +(Y(64 )*(-K(1032)*X(2 )*D  -K(1033)*X(2 )
     *        *D  -K(1199)*D ))                                     
      YD( 1) = 
     *        +K(  46)*Y(22 )*Y(72 )*D  +K(  61)*Y(34 )*Y(67 )*D  +K
     *        (  63)*Y(34 )*Y(72 )*D  +K( 104)*Y(1  )*Y(67 )*D  +K( 
     *        167)*Y(3  )*Y(66 )*D  +K( 288)*Y(7  )*Y(72 )*D  +K( 35
     *        7)*Y(10 )*Y(66 )*D  +K( 410)*Y(13 )*Y(68 )*D  +K( 553)
     *        *Y(22 )*Y(77 )*D  +K( 676)*Y(27 )*Y(69 )*D  +K( 700)*Y
     *        (28 )*Y(66 )*D  +K( 701)*Y(28 )*Y(68 )*D  +K( 752)*Y(3
     *        4 )*Y(73 )*D  +K( 754)*Y(34 )*Y(78 )*D  +K( 830)*Y(40 
     *        )*Y(68 )*D  +K( 861)*Y(43 )*Y(66 )*D  +K( 874)*Y(45 )*
     *        Y(66 )*D  +K( 910)*Y(54 )*Y(68 )*D  +K( 926)*Y(56 )*Y(
     *        68 )*D  +K( 932)*Y(57 )*Y(66 )*D  +K( 941)*Y(67 )*Y(67
     *         )*D  +K( 942)*Y(68 )*Y(69 )*D  +K(1034)*Y(66 )*X(2 )*
     *        D  +K(1035)*Y(68 )*X(2 )*D  +K(1037)*Y(70 )*X(2 )*D  +
     *        K(1041)*Y(73 )*X(2 )*D  +K(1046)*Y(78 )*X(2 )*D  +K(11
     *        20)*Y(67 ) +K(1123)*Y(69 ) +K(1125)*Y(72 ) +K(1171)*Y(
     *        66 )*D +K(1200)*Y(68 )*D +K(1203)*Y(73 )*D 
     *        +K(1206)*Y(78 )*D +(Y(65 )*(-K(  22)*Y(13 )*D -
     *        K(  52)*Y(24 )*D -K( 134)*Y(2  )*D -K( 226)*Y(5  )*D -
     *        K( 316)*Y(8  )*D -K( 341)*Y(9  )*D ))                 
      YDOT( 65) = +YD( 1)
     *        +(Y( 65)*(
     *        -K( 342)*Y(9  )*D -K( 368)*Y(12 )*D -K( 408)*Y(13 )*D 
     *        -K( 434)*Y(14 )*D -K( 435)*Y(14 )*D -K( 436)*Y(14 )*D 
     *        -K( 471)*Y(15 )*D -K( 472)*Y(15 )*D -K( 481)*Y(16 )*D 
     *        -K( 499)*Y(18 )*D -K( 533)*Y(21 )*D -K( 637)*Y(25 )*D 
     *        -K( 638)*Y(25 )*D -K( 674)*Y(27 )*D -K( 675)*Y(27 )*D 
     *        -K( 723)*Y(32 )*D -K( 728)*Y(33 )*D -K( 776)*Y(37 )*D 
     *        -K( 809)*Y(39 )*D -K( 810)*Y(39 )*D -K( 843)*Y(41 )*D 
     *        -K( 844)*Y(41 )*D -K( 877)*Y(46 )*D -K( 886)*Y(48 )*D 
     *        -K( 887)*Y(48 )*D -K( 903)*Y(53 )*D -K( 920)*Y(55 )*D 
     *        -K( 921)*Y(55 )*D -K( 936)*Y(60 )*D -K( 939)*Y(68 )*D 
     *        -K( 940)*Y(70 )*D -K(1119) -K(1156)*D ))              
      YD( 1) = 
     *        +K( 105)*Y(1  )*Y(68 )*D  +K( 134)*Y(2  )*Y(65 )*D  +K
     *        ( 135)*Y(2  )*Y(67 )*D  +K( 138)*Y(2  )*Y(69 )*D  +K( 
     *        200)*Y(4  )*Y(69 )*D  +K( 283)*Y(7  )*Y(67 )*D  +K( 28
     *        4)*Y(7  )*Y(69 )*D  +K( 287)*Y(7  )*Y(72 )*D  +K( 341)
     *        *Y(9  )*Y(65 )*D  +K( 368)*Y(12 )*Y(65 )*D  +K( 436)*Y
     *        (14 )*Y(65 )*D  +K( 585)*Y(23 )*Y(69 )*D  +K( 638)*Y(2
     *        5 )*Y(65 )*D  +K( 675)*Y(27 )*Y(65 )*D  +K( 723)*Y(32 
     *        )*Y(65 )*D  +K( 724)*Y(32 )*Y(69 )*D  +K( 750)*Y(34 )*
     *        Y(68 )*D  +K( 773)*Y(35 )*Y(69 )*D  +K( 776)*Y(37 )*Y(
     *        65 )*D  +K( 810)*Y(39 )*Y(65 )*D  +K( 844)*Y(41 )*Y(65
     *         )*D  +K( 887)*Y(48 )*Y(65 )*D  +K( 903)*Y(53 )*Y(65 )
     *        *D  +K( 921)*Y(55 )*Y(65 )*D  +K( 939)*Y(65 )*Y(68 )*D
     *          +K( 940)*Y(65 )*Y(70 )*D  +K(1119)*Y(65 ) +K(1121)*Y
     *        (68 ) +K(1127)*Y(73 ) +(Y(66 )*(-K(  13)*Y(8  )*D -K( 
     *        167)*Y(3  )*D -K( 357)*Y(10 )*D -K( 409)*Y(13 )*D -K( 
     *        473)*Y(15 )*D -K( 489)*Y(17 )*D -K( 515)*Y(19 )*D -K( 
     *        700)*Y(28 )*D -K( 861)*Y(43 )*D -K( 873)*Y(45 )*D -K( 
     *        874)*Y(45 )*D -K( 884)*Y(47 )*D ))                    
      YDOT( 66) = +YD( 1)
     *        +(Y( 66)*(
     *        -K( 885)*Y(47 )*D -K( 932)*Y(57 )*D -K( 963)*X(1 )*D  
     *        -K(1034)*X(2 )*D  -K(1171)*D ))                       
      YDOT( 67) =      +K(   9)*Y(1  )*Y(69 )*D  +K(  22)*Y(13 
     *        )*Y(65 )*D  +K(  38)*Y(17 )*Y(69 )*D  +K(  52)*Y(24 )*
     *        Y(65 )*D  +K(  62)*Y(34 )*Y(69 )*D  +K(  68)*Y(38 )*Y(
     *        69 )*D  +K( 369)*Y(12 )*Y(69 )*D  +K( 584)*Y(23 )*Y(69
     *         )*D  +K( 678)*Y(27 )*Y(69 )*D  +K( 702)*Y(28 )*Y(68 )
     *        *D  +K( 703)*Y(28 )*Y(70 )*D  +K( 714)*Y(29 )*Y(69 )*D
     *          +K( 831)*Y(40 )*Y(70 )*D  +K( 845)*Y(41 )*Y(69 )*D  
     *        +K( 862)*Y(43 )*Y(68 )*D  +K( 884)*Y(47 )*Y(66 )*D  +K
     *        ( 933)*Y(57 )*Y(68 )*D  +K( 939)*Y(65 )*Y(68 )*D  +K( 
     *        944)*Y(69 )*Y(70 )*D  +K(1036)*Y(70 )*X(2 )*D  +K(1039
     *        )*Y(71 )*X(2 )*D  +K(1122)*Y(69 ) +K(1201)*Y(70 )*D +(
     *        Y(67 )*(-K(  61)*Y(34 )*D -K( 104)*Y(1  )*D -K( 135)*Y
     *        (2  )*D -K( 136)*Y(2  )*D -K( 227)*Y(5  )*D -K( 283)*Y
     *        (7  )*D -K( 317)*Y(8  )*D -K( 343)*Y(9  )*D -K( 500)*Y
     *        (18 )*D -K( 878)*Y(46 )*D -2*K( 941)*Y(67 )*D -K(1120)
     *         -K(1157)*D ))                                        
      YDOT( 68) =      +K( 106)*Y(1  )*Y(70 )*D  +K( 136)*Y(2  
     *        )*Y(67 )*D  +K( 137)*Y(2  )*Y(69 )*D  +K( 201)*Y(4  )*
     *        Y(69 )*D  +K( 226)*Y(5  )*Y(65 )*D  +K( 285)*Y(7  )*Y(
     *        69 )*D  +K( 434)*Y(14 )*Y(65 )*D  +K( 533)*Y(21 )*Y(65
     *         )*D  +K( 583)*Y(23 )*Y(69 )*D  +K( 637)*Y(25 )*Y(65 )
     *        *D  +K( 674)*Y(27 )*Y(65 )*D  +K( 677)*Y(27 )*Y(69 )*D
     *          +K( 725)*Y(32 )*Y(69 )*D  +K( 728)*Y(33 )*Y(65 )*D  
     *        +K( 751)*Y(34 )*Y(70 )*D  +K( 772)*Y(35 )*Y(69 )*D  +K
     *        ( 809)*Y(39 )*Y(65 )*D  +K( 843)*Y(41 )*Y(65 )*D  +K( 
     *        873)*Y(45 )*Y(66 )*D  +K( 877)*Y(46 )*Y(65 )*D  +K( 88
     *        6)*Y(48 )*Y(65 )*D  +K( 920)*Y(55 )*Y(65 )*D  +K( 936)
     *        *Y(60 )*Y(65 )*D  +(Y(68 )*(-K( 105)*Y(1  )*D -K( 318)
     *        *Y(8  )*D -K( 410)*Y(13 )*D -K( 701)*Y(28 )*D -K( 702)
     *        *Y(28 )*D -K( 750)*Y(34 )*D -K( 830)*Y(40 )*D -K( 862)
     *        *Y(43 )*D -K( 910)*Y(54 )*D -K( 926)*Y(56 )*D -K( 933)
     *        *Y(57 )*D -K( 939)*Y(65 )*D -K( 942)*Y(69 )*D -K( 964)
     *        *X(1 )*D  -K(1035)*X(2 )*D  -K(1121) -K(1200)*D ))    
      YD( 1) = 
     *        +K( 704)*Y(28 )*Y(70 )*D  +K( 705)*Y(28 )*Y(71 )*D  +K
     *        ( 863)*Y(43 )*Y(70 )*D  +K( 875)*Y(45 )*Y(70 )*D  +K( 
     *        911)*Y(54 )*Y(71 )*D  +K( 927)*Y(56 )*Y(71 )*D  +K( 93
     *        4)*Y(57 )*Y(70 )*D  +K( 940)*Y(65 )*Y(70 )*D  +K( 941)
     *        *Y(67 )*Y(67 )*D  +K(1038)*Y(70 )*X(2 )*D  +K(1040)*Y(
     *        71 )*X(2 )*D  +K(1202)*Y(71 )*D +K(1223)*Y(92 ) +(Y(69
     *         )*(-K(   9)*Y(1  )*D -K(  38)*Y(17 )*D -K(  62)*Y(34 
     *        )*D -K(  68)*Y(38 )*D -K( 137)*Y(2  )*D -K( 138)*Y(2  
     *        )*D -K( 139)*Y(2  )*D -K( 200)*Y(4  )*D -K( 201)*Y(4  
     *        )*D -K( 202)*Y(4  )*D -K( 228)*Y(5  )*D -K( 284)*Y(7  
     *        )*D -K( 285)*Y(7  )*D -K( 286)*Y(7  )*D -K( 344)*Y(9  
     *        )*D -K( 345)*Y(9  )*D -K( 369)*Y(12 )*D -K( 370)*Y(12 
     *        )*D -K( 437)*Y(14 )*D -K( 438)*Y(14 )*D -K( 482)*Y(16 
     *        )*D -K( 483)*Y(16 )*D -K( 525)*Y(20 )*D -K( 526)*Y(20 
     *        )*D -K( 534)*Y(21 )*D -K( 583)*Y(23 )*D -K( 584)*Y(23 
     *        )*D -K( 585)*Y(23 )*D -K( 586)*Y(23 )*D -K( 676)*Y(27 
     *        )*D -K( 677)*Y(27 )*D -K( 678)*Y(27 )*D -K( 679)*Y(27 
     *        )*D -K( 680)*Y(27 )*D -K( 714)*Y(29 )*D ))            
      YDOT( 69) = +YD( 1)
     *        +(Y( 69)*(
     *        -K( 724)*Y(32 )*D -K( 725)*Y(32 )*D -K( 726)*Y(32 )*D 
     *        -K( 772)*Y(35 )*D -K( 773)*Y(35 )*D -K( 774)*Y(35 )*D 
     *        -K( 777)*Y(37 )*D -K( 811)*Y(39 )*D -K( 812)*Y(39 )*D 
     *        -K( 845)*Y(41 )*D -K( 846)*Y(41 )*D -K( 847)*Y(41 )*D 
     *        -K( 851)*Y(42 )*D -K( 879)*Y(46 )*D -K( 937)*Y(61 )*D 
     *        -K( 938)*Y(61 )*D -K( 942)*Y(68 )*D -K( 943)*Y(50 )*D 
     *        -K( 944)*Y(70 )*D -K(1122) -K(1123) -K(1124) -K(1158)*
     *        D ))                                                  
      YDOT( 70) =      +K( 107)*Y(1  )*Y(71 )*D  +K( 139)*Y(2  
     *        )*Y(69 )*D  +K( 202)*Y(4  )*Y(69 )*D  +K( 227)*Y(5  )*
     *        Y(67 )*D  +K( 286)*Y(7  )*Y(69 )*D  +K( 345)*Y(9  )*Y(
     *        69 )*D  +K( 370)*Y(12 )*Y(69 )*D  +K( 526)*Y(20 )*Y(69
     *         )*D  +K( 586)*Y(23 )*Y(69 )*D  +K( 680)*Y(27 )*Y(69 )
     *        *D  +K( 726)*Y(32 )*Y(69 )*D  +K( 774)*Y(35 )*Y(69 )*D
     *          +K( 777)*Y(37 )*Y(69 )*D  +K( 812)*Y(39 )*Y(69 )*D  
     *        +K( 847)*Y(41 )*Y(69 )*D  +K( 878)*Y(46 )*Y(67 )*D  +K
     *        ( 885)*Y(47 )*Y(66 )*D  +K( 943)*Y(50 )*Y(69 )*D  +K( 
     *        963)*Y(66 )*X(1 )*D  +K(1124)*Y(69 ) +(Y(70 )*(-K(  80
     *        )*X(1 )*D  -K( 106)*Y(1  )*D -K( 319)*Y(8  )*D -K( 703
     *        )*Y(28 )*D -K( 704)*Y(28 )*D -K( 751)*Y(34 )*D -K( 831
     *        )*Y(40 )*D -K( 863)*Y(43 )*D -K( 875)*Y(45 )*D -K( 934
     *        )*Y(57 )*D -K( 940)*Y(65 )*D -K( 944)*Y(69 )*D -K(1036
     *        )*X(2 )*D  -K(1037)*X(2 )*D  -K(1038)*X(2 )*D  -K(1201
     *        )*D ))                                                
      YDOT( 71) =      +K(  80)*Y(70 )*X(1 )*D  +K( 228)*Y(5  )
     *        *Y(69 )*D  +K( 437)*Y(14 )*Y(69 )*D  +K( 482)*Y(16 )*Y
     *        (69 )*D  +K( 525)*Y(20 )*Y(69 )*D  +K( 534)*Y(21 )*Y(6
     *        9 )*D  +K( 679)*Y(27 )*Y(69 )*D  +K( 811)*Y(39 )*Y(69 
     *        )*D  +K( 846)*Y(41 )*Y(69 )*D  +K( 851)*Y(42 )*Y(69 )*
     *        D  +K( 879)*Y(46 )*Y(69 )*D  +K( 937)*Y(61 )*Y(69 )*D 
     *         +K( 938)*Y(61 )*Y(69 )*D  +K( 942)*Y(68 )*Y(69 )*D  +
     *        K( 944)*Y(69 )*Y(70 )*D  +K( 964)*Y(68 )*X(1 )*D  +(Y(
     *        71 )*(-K( 107)*Y(1  )*D -K( 705)*Y(28 )*D -K( 911)*Y(5
     *        4 )*D -K( 927)*Y(56 )*D -K(1039)*X(2 )*D  -K(1040)*X(2
     *         )*D  -K(1202)*D ))                                   
      YDOT( 72) =      +K( 289)*Y(7  )*Y(74 )*D  +K( 291)*Y(7  
     *        )*Y(77 )*D  +K( 316)*Y(8  )*Y(65 )*D  +K( 317)*Y(8  )*
     *        Y(67 )*D  +K( 408)*Y(13 )*Y(65 )*D  +K( 472)*Y(15 )*Y(
     *        65 )*D  +K( 706)*Y(28 )*Y(78 )*D  +K( 753)*Y(34 )*Y(77
     *         )*D  +K( 864)*Y(43 )*Y(73 )*D  +K(1042)*Y(75 )*X(2 )*
     *        D  +K(1043)*Y(76 )*X(2 )*D  +K(1045)*Y(78 )*X(2 )*D  +
     *        K(1047)*Y(79 )*X(2 )*D  +K(1128)*Y(74 ) +K(1207)*Y(79 
     *        )*D +(Y(72 )*(-K(  46)*Y(22 )*D -K(  63)*Y(34 )*D -K( 
     *        140)*Y(2  )*D -K( 229)*Y(5  )*D -K( 287)*Y(7  )*D -K( 
     *        288)*Y(7  )*D -K( 852)*Y(42 )*D -K( 880)*Y(46 )*D -K(1
     *        125) -K(1126) -K(1159)*D )) +K(1204)*Y(75)*D                         
      YDOT( 73) =      +K(  13)*Y(8  )*Y(66 )*D  +K( 140)*Y(2  
     *        )*Y(72 )*D  +K( 142)*Y(2  )*Y(77 )*D  +K( 290)*Y(7  )*
     *        Y(74 )*D  +K( 292)*Y(7  )*Y(77 )*D  +K( 318)*Y(8  )*Y(
     *        68 )*D  +K( 342)*Y(9  )*Y(65 )*D  +K( 343)*Y(9  )*Y(67
     *         )*D  +K( 409)*Y(13 )*Y(66 )*D  +K( 435)*Y(14 )*Y(65 )
     *        *D  +K(1126)*Y(72 ) +(Y(73 )*(-K( 516)*Y(19 )*D -K( 75
     *        2)*Y(34 )*D -K( 864)*Y(43 )*D -K( 965)*X(1 )*D  -K(104
     *        1)*X(2 )*D  -K(1127) -K(1203)*D ))                    
      YDOT( 74) =      +K(1044)*Y(76 )*X(2 )*D  +K(1205)*Y(76 )
     *        *D +K(1224)*Y(93 ) +(Y(74 )*(-K( 141)*Y(2  )*D -K( 230
     *        )*Y(5  )*D -K( 289)*Y(7  )*D -K( 290)*Y(7  )*D -K( 346
     *        )*Y(9  )*D -K( 853)*Y(42 )*D -K( 881)*Y(46 )*D -K(1128
     *        ) -K(1160)*D ))                                       
      YDOT( 75) =      +K( 141)*Y(2  )*Y(74 )*D  +K( 346)*Y(9  
     *        )*Y(74 )*D  +(Y(75 )*(-K( 966)*X(1 )*D  -K(1042)*X(2 )
     *        *D  -K(1204)*D ))                                     
      YDOT( 76) =      +K( 230)*Y(5  )*Y(74 )*D  +K( 853)*Y(42 
     *        )*Y(74 )*D  +K( 881)*Y(46 )*Y(74 )*D  +K( 966)*Y(75 )*
     *        X(1 )*D  +(Y(76 )*(-K(1043)*X(2 )*D  -K(1044)*X(2 )*D 
     *         -K(1205)*D ))                                        
      YDOT( 77) =      +K( 471)*Y(15 )*Y(65 )*D  +K(1048)*Y(79 
     *        )*X(2 )*D  +K(1225)*Y(94 ) +(Y(77 )*(-K( 142)*Y(2  )*D
     *         -K( 231)*Y(5  )*D -K( 291)*Y(7  )*D -K( 292)*Y(7  )*D
     *         -K( 553)*Y(22 )*D -K( 753)*Y(34 )*D -K(1129) -K(1161)
     *        *D ))                                                 
      YDOT( 78) =      +K( 229)*Y(5  )*Y(72 )*D  +K( 319)*Y(8  
     *        )*Y(70 )*D  +K( 344)*Y(9  )*Y(69 )*D  +K( 438)*Y(14 )*
     *        Y(69 )*D  +K( 473)*Y(15 )*Y(66 )*D  +K( 481)*Y(16 )*Y(
     *        65 )*D  +K( 483)*Y(16 )*Y(69 )*D  +K( 499)*Y(18 )*Y(65
     *         )*D  +K( 515)*Y(19 )*Y(66 )*D  +K( 516)*Y(19 )*Y(73 )
     *        *D  +K( 852)*Y(42 )*Y(72 )*D  +K( 880)*Y(46 )*Y(72 )*D
     *          +K( 965)*Y(73 )*X(1 )*D  +K(1129)*Y(77 ) +(Y(78 )*(-
     *        K( 706)*Y(28 )*D -K( 754)*Y(34 )*D -K(1045)*X(2 )*D  -
     *        K(1046)*X(2 )*D  -K(1206)*D ))                        
      YDOT( 79) =      +K( 231)*Y(5  )*Y(77 )*D  +K( 489)*Y(17 
     *        )*Y(66 )*D  +K( 500)*Y(18 )*Y(67 )*D  +(Y(79 )*(-K(104
     *        7)*X(2 )*D  -K(1048)*X(2 )*D  -K(1207)*D ))           
      YDOT( 80) =      +K(1131)*Y(8  )*D +K(1134)*Y(13 )*D +K(1
     *        135)*Y(15 )*D +K(1136)*Y(17 )*D +K(1137)*Y(19 )*D +K(1
     *        162)*Y(10 )*D +K(1208)*Y(9  )*D +(Y(80 )*(-K(1212) )) 
      YDOT( 81) =      +K(1132)*Y(11 )*D +(Y(81 )*(-K(1211) )) 
      YDOT( 82) =      +K(1138)*Y(22 )*D +K(1139)*Y(24 )*D +K(1
     *        140)*Y(26 )*D +K(1141)*Y(28 )*D +K(1209)*Y(23 )*D +(Y(
     *        82 )*(-K(1213) ))                                     
      YDOT( 83) =      +K(1142)*Y(31 )*D +(Y(83 )*(-K(1214) )) 
      YDOT( 84) =      +K(1144)*Y(36 )*D +(Y(84 )*(-K(1215) ))
      YDOT( 85) =      +K(1143)*Y(34 )*D +K(1145)*Y(38 )*D +K(1
     *        146)*Y(40)*D+K(1210)*Y(35)*D +(Y(85 )*(-K(1216)-PHOH)) 
      YDOT( 86) =      +K(1147)*Y(43 )*D +(Y(86 )*(-K(1217) )) 
      YDOT( 87) =      +K(1148)*Y(45 )*D +K(1149)*Y(47 )*D +(Y(
     *        87 )*(-K(1218) ))                                     
      YDOT( 88) =      +K(1133)*Y(11 )*D +K(1150)*Y(49 )*D +(Y(
     *        88 )*(-K(1219) ))                                     
      YDOT( 89) =      +K(1151)*Y(52 )*D +K(1152)*Y(54 )*D +(Y(
     *        89 )*(-K(1220) ))                                     
      YDOT( 90) =      +K(1153)*Y(56 )*D +(Y(90 )*(-K(1221) )) 
      YDOT( 91) =      +K(1154)*Y(57 )*D +K(1155)*Y(59 )*D +(Y(
     *        91 )*(-K(1222) ))                                     
      YDOT( 92) =      +K(1156)*Y(65 )*D +K(1157)*Y(67 )*D +K(1
     *        158)*Y(69 )*D +(Y(92 )*(-K(1223) ))                   
      YDOT( 93) =      +K(1160)*Y(74 )*D +(Y(93 )*(-K(1224) )) 
      YDOT( 94) =      +K(1159)*Y(72 )*D +K(1161)*Y(77 )*D +(Y(
     *        94 )*(-K(1225) ))                                     
C--include conversion of O,O+,OH to CO2 on dust grains
      OFR=K(1231)*Y(34)*D
      OPLFR=K(1232)*Y(35)*D
      OHFR=K(1233)*Y(38)*D
      YDOT(88)=YDOT(88)+(OFR+OHFR+OPLFR)
      YDOT(81)=YDOT(81)-(OFR+OHFR+OPLFR)
C--include freeze-out/desorption of extra species: GO and GOH
      YDOT(95)=K(1226)*Y(34)*D+K(1227)*Y(35)*D-K(1229)*Y(95)
      YDOT(96)=K(1228)*Y(38)*D-K(1230)*Y(96)
      YDOT(34)=YDOT(34)-K(1226)*Y(34)*D+K(1229)*Y(95)-OFR
      YDOT(35)=YDOT(35)-K(1227)*Y(35)*D-OPLFR
      YDOT(38)=YDOT(38)-K(1228)*Y(38)*D+K(1230)*Y(96)-OHFR
C--include dissociative desorption: GH2O -> OH + H
      YDOT(1)=YDOT(1)+PHOH*Y(85)
      YDOT(38)=YDOT(38)+PHOH*Y(85)
C
      RETURN
      END

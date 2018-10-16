      PROGRAM  WARDLE
C
C------------------------------------------------------------------------------
C
C  Single point calculation; model of time-dep. chemistry in a 
C  static, dust-poor cloud
C
C  Based on L1544/L1544.f [15th August 2013]
C
C  H,He,C,N,O,S,Na chemistry (gas-phase, freeze-out and desorption)
C     Only H2 and Electrons are conserved.
C
C  NG     is the number of time-dependent (gas-phase) species
C  NCONS  is the number of conserved species
C  NS     is the number of time-dependent (solid-state) species
C
C  NSP    is >= the total number of species
C  NRM    is >= the number of reactions
C  NOP    is >= the number of output times
C
C------------------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER(NSP=100,NEL=6,NCONS=2,NRM=1500,NOP=2228)
C
      CHARACTER*8  RE1(NRM),RE2(NRM),P1(NRM),P2(NRM),P3(NRM),P4(NRM),
     *             SPECI(NSP),SPX(NEL)
      CHARACTER*25 UNIT2,UNIT3,UNIT7,UNIT8,UNITL,DGFILE,TENFILE,RHONFILE
      REAL*8  Y0(NSP),FRAC(NEL),DPLT(NEL),TOTAL(NCONS),X(NCONS),
     * K(NRM),TRAT(NRM,4),PRAT(NRM,4),DRAT(NRM,4),ERAT(NRM,4),
     * B(NSP,NOP),TAGE(NOP),T_INTERP(NOP),DG(NOP),TEN(NOP),
     * RHON(NOP),NON(NOP),DUMMY(NOP)
      INTEGER INDR(NRM),ISPX(NEL)
C
C** For LSODE
      REAL*8  RWORK(11000)
      INTEGER IWORK(120)
      CHARACTER*5 JAC
      EXTERNAL DIFFUN
C**
      COMMON/BLK1/SPECI
      COMMON/BLK2/B,TAGE
      COMMON/BLK3/X,TOTAL,K,D
      COMMON/BLK4/TRAT,PRAT,DRAT,ERAT,ITR,IPR,IDR,IER,IH2
      COMMON/BLK5/IFC,ICOV,ISURF
      COMMON/BLK6/DEN0,TEMP0,AV0
      COMMON/BLK8/GRAD,SAPH,YLD,YLD2,FSITES,G0,F0,FCR,YH,CRAT,
     *            STICK0,STICKP,STICKN
      COMMON/BLK9/FCO2,COCOV,FOX,PHOH,H2SHL,COSHL,FFRZ,FDES
      COMMON/BLK10/T0,TOUT,DG0,DGOUT,TEN0,TENOUT,NON0,NONOUT
      COMMON/ORDER1/RE1,RE2,P1,P2,P3,P4
      COMMON/ORDER2/INDR,IOTYP
C
      DATA  NG,NS/79,17/
C-- FRAC, DPLT are the elemental abundances and depletions factors for -----
C   species (SPX) whose index number in the species set is given by ISPX
C   It is assumed that conserved species are (1) H2 and (2) electrons.
C   IHI AND ICO are the index numbers for H and CO respectively.
C---------------------------------------------------------------------------
      DATA  FRAC/0.1,3.75E-4,1.15E-4,6.74E-4,1.62E-5,3.5E-5/
C--Hollenbach et al. (2009) values
C      DATA  DPLT/1.0,0.373,1.0,0.475,1.73,0.086/
C--Asplund, Grevesse & Sauval values (2005)
      DATA  DPLT/1.0,0.68,0.53,0.68,0.85,1.82/
C      DATA  DPLT/1.0,0.5,1.0,1.0,0.01,0.086/
      DATA  SPX/'HE','C','N','O','S','NA'/
      DATA  ISPX/6,8,22,34,65,43/
      DATA  IHI,ICO/1,11/
C
      DATA PC,YEAR/3.0856D18,3.1557D7/
      DATA BOLTK,AMU,PI/1.381D-23,1.67D-27,3.1415927/
C
C** For LSODE
      DATA ITOL,MF/1,22/
      JAC='DUMMY'
      ITASK=1
      IOPT=0
      LRW=11000
      LIW=120
      RTOL=1.0E-6
      ATOL=1.0E-16
C**
C
C--NTD is the total number of time-dependent chemical species
      NTD=NG+NS
C--N is the total number of time-dependent species
      N=NTD
C--NTOT is the total number of chemical species
      NTOT=NTD+NCONS
C---------------------------- INPUT PARAMETERS --------------------------------
C----I/O files
      LOUT=3
      LRUN=6
      UNIT2='ordout.d'
      UNIT3='output.d'
      UNITL='report.d'
      UNIT7='coll_specs.d'
      UNIT8='coll_rates.d'
      DGFILE='interp/dg_interp.dat'
      TENFILE='interp/ten_interp.dat'
      RHONFILE='interp/rhon_interp.dat'

C------------------------------------------------------------------------------
C
C     Read the data
      OPEN(20,FILE=DGFILE,STATUS='OLD')
      OPEN(21,FILE=TENFILE,STATUS='OLD')
      OPEN(22,FILE=RHONFILE,STATUS='OLD')

      DO i=1,NOP
            READ(20,*) T_INTERP(i),DG(i)
            READ(21,*) DUMMY(i),TEN(i)
            READ(22,*) DUMMY(i),RHON(i)

C           Convert mass density to number density
            NON(i)=RHON(i)/(2*1.67D-24)
      END DO
C
C------------------------------------------------------------------------------
C--ORDER  is called at time points: IT1,IT2,IT3,IT4.
C  IOTYP=0,1,2 (selected/gas-phase/full o/p)
      IT1=2
      IT2=15
      IT3=37
      IT4=53
      IOTYP=1
C--IFULL=0,1 for reduced/full output of abundances 
      IFULL=0
C--IPLOT=0,1 for reduced/full plotting
      IPLOT=1
C--Density, temperature and extinction
C      DEN0=1.0D6
C      TEMP0=10.0
C      AV0=5.0
      DEN0=NON(1)
      TEMP0=TEN(1)
      AV0=5.0
C--Fraction of dust present
C      FDUST=0.01
      FDUST=DG(1)
C--Metallicity fraction
      FRM=100.0
C
      XFRAC=1.0D-20
      HFRAC=1.0/DEN0
C
C      TMIN=100.0*YEAR
C      TMAX=1.0E7*YEAR
C      TINC=10.0**(0.1)
      TMAX=T_INTERP(NOP)
C--Photochemistry parameters
      CRAT=1.3E-16
      DEFF=200.0
      ALBEDO=0.5
C--H2 and CO shielding factors
      H2SHL=0.0
      COSHL=0.0
C
C--Freeze-out parameters: GRAD = grain radius (in microns), SAPH = grain 
C  surface area per H-nucleon (in cm2).
C  STICK0,STICKP,STICKN are the sticking efficiencies for neutrals,+ve,-ve ions
      GRAD=0.0083
      SAPH=8.0E-21*FDUST
      STICK0=1.0
      STICKP=1.0
      STICKN=1.0
C--Surface chemistry parameters: 
C  FCO2 is the proportion of O,O+,OH (+CO ice) that converts to CO2 
      FCO2=0.1
C  FOX is the fraction of the remainder O,O+,OH that converts to H2O
      FOX=1.0
C--Desorption parameters: 
C  YLD = (direct) photodesorption yield, 
C  YLD2= (CR-induced) photodesorption yield,
C      YLD/YLD2 doubled for GH2O->OH+H, reduced by 10 for GO->O
C  G0 = flux factor (multiples of Draine), 
C  F0 = I/S photon flux, 
C  FCR = CR-induced photon flux (for zeta_0),
C  SBIND = surface density of binding sites, 
C  CRDEF = default value of the CR-heating desorption rate, if not known,
C  YH = desorption yield due to H2 formation,
C  ICOV = 1/0 switches on/off H2O desorption protection by CO/CO2/CH4.
C
      YLD=1.0E-3
      YLD2=1.0E-3
      G0=1.0
      F0=1.0E8
      FCR=4875.0
      SBIND=1.0E15
      CRDEF=1.0E-17
      YH=0.0
      ICOV=0
C--FFRZ,FDES are flags (on:1.0, off:0.00) for freeze-out & desorption
C  ISURF controls the surface chemistry (1:on, 0:off)
      FFRZ=1.0
      FDES=1.0
      ISURF=1
C-------------------- END OF INPUT PARAMETERS ---------------------------------
      OPEN(LOUT,FILE=UNIT3,STATUS='NEW')
C--Write out the input parameters to the log file
      WRITE(LOUT,100)
 100  FORMAT(1X,'INPUT PARAMETERS:',/)
      DO 101 I=1,NEL
 101    WRITE(LOUT,102) SPX(I),FRAC(I),DPLT(I)
 102  FORMAT(1X,A7,'/H = ',1PE8.2,'*',1PE8.2)
      WRITE(LOUT,108) DEN0,TEMP0,AV0
 108  FORMAT(/,1X,'Density = ',1PE8.2,' cm-3',/,
     *   1X,'Temperature = ',0PF5.1,' K',/,
     *   1X,'Extinction = ',1PE8.2,' magnitudes',/)
      WRITE(LOUT,103) ALBEDO,CRAT,DEFF,H2SHL,COSHL
 103  FORMAT(1X,'Albedo = ',1PE8.2,/,
     *   1X,'Cosmic ray ionization rate = ',1PE8.2,' s-1',/,
     *   1X,'Default CRP-induced photolysis rate = ',1PE8.2,/,
     *   1X,'H2 self-shielding factor = ',0PF5.3,/,
     *   1X,'CO self-shielding factor = ',0PF5.3,/)
      WRITE(LOUT,104) FDUST,FRM,GRAD,SAPH,STICK0,STICKP,STICKN,
     *                FFRZ,FDES,ISURF
 104  FORMAT(1X,'Fraction of dust present = ',0PF5.3,/,
     *   1X,'Metallicity fraction = ',0PF6.2,/,
     *   1X,'Grain radius = ',1PE8.2,' microns',/,
     *   1X,'Eff. dust surface area per H = ',1PE8.2,' cm-2',/,
     *   1X,'Sticking coefficient (neutrals) = ',1PE8.2,/,
     *   1X,'Sticking coefficient (positive ions) = ',1PE8.2,/,
     *   1X,'Sticking coefficient (negative ions) = ',1PE8.2,/,
     *   1X,'Freeze-out efficiency = ',0PF5.3,/,
     *   1X,'Desorption efficiency = ',0PF5.3,/,
     *   1X,'Surface chemistry flag = ',I1,/)
      WRITE(LOUT,105) FCO2,FOX,ICOV,YLD,YLD2,G0,F0,FCR,SBIND,YH
 105  FORMAT(1X,'G(O/OH)+GCO->GCO2 conversion efficiency = ',1PE8.2,/,
     *   1X,'G(O/OH)[rem.]->GH2O conversion efficiency = ',1PE8.2,/,
     *   1X,'H2O desorp. shielding by CO/CO2/CH4 (1,0;on,off) = ',I1,/,
     *   1X,'(Direct) photodesorption yield = ',1PE8.2,/,
     *   1X,'(CR-induced) photodesorption yield = ',1PE8.2,/,
     *   1X,'ISRF/Draine = ',1PE8.2,/,
     *   1X,'ISRF (flux) = ',1PE8.2,' photons cm-2 s-1',/,
     *   1X,'CR-induced flux (norm.) = ',1PE8.2,' photons cm-2 s-1',/,
     *   1X,'Binding site surface density = ',1PE8.2,' cm-2',/,
     *   1X,'H2 formation desorption yield = ',1PE8.2,/)
      WRITE(LOUT,106) IFULL
 106  FORMAT(1X,'Full/Partial results output flag = ',I1,/)
      WRITE(LOUT,107)
 107  FORMAT(1X,75('-'))
C
C------------------------------------------------------------------------------
C
C  Read the species file
C
      SPECI(1)='H2'
      SPECI(2)='ELECTR'
      OPEN(7,FILE=UNIT7,STATUS='OLD')
      READ(7,1)(INDXS,SPECI(J+2),J=1,NTD)
 1    FORMAT(5(2X,I3,2X,A8,1X))
      CLOSE(7)
C
C------------------------------------------------------------------------------
C
C  Read reaction data
C
      OPEN(8,FILE=UNIT8,STATUS='OLD')
      ITR=0
      IPR=0
      IDR=0
      IER=0
C
      ALBFAC=1.0/(1.0-ALBEDO)
      FSITES=SAPH*SBIND
C
      J=0
 10   J=J+1
      READ(8,3) INDR(J),RE1(J),RE2(J),P1(J),P2(J),P3(J),P4(J),
     *          GAMMA,ALPHA,BETA
 3    FORMAT(I4,4(1X,A8),2(1X,A4),1X,1PE8.2,1X,0PF5.2,1X,F8.1)
      INDXJ=INDR(J)
C
C  store reaction data
C
      IF(INDXJ.EQ.9999) GO TO 11
      IF(RE2(J).EQ.'CRP') THEN
          K(INDXJ)=GAMMA*CRAT
      ELSE IF((RE2(J).EQ.'PHOTON').AND.(P2(J).NE.'G')) THEN
          IPR=IPR+1
          IF(BETA.EQ.0.0) BETA=DEFF
          PRAT(IPR,1)=INDXJ
          PRAT(IPR,2)=GAMMA
          PRAT(IPR,3)=ALPHA
          PRAT(IPR,4)=BETA*ALBFAC*CRAT 
      ELSE IF(RE2(J).EQ.'G') THEN
          IDR=IDR+1
          DRAT(IDR,1)=INDXJ
          DRAT(IDR,2)=GAMMA
          DRAT(IDR,3)=ALPHA
          DRAT(IDR,4)=BETA
C--Identify reaction number for: H + grain -> H2
          IF(RE1(J).EQ.'H') IH2=INDXJ
      ELSE IF((RE2(J).EQ.'PHOTON').AND.(P2(J).EQ.'G')) THEN
          IER=IER+1
          ERAT(IER,1)=INDXJ
          IF(GAMMA.EQ.0.0) GAMMA=CRDEF
          ERAT(IER,2)=GAMMA
          ERAT(IER,3)=ALPHA
          ERAT(IER,4)=BETA
      ELSE IF((ALPHA.NE.0.0).OR.(BETA.NE.0.0)) THEN
          ITR=ITR+1
          TRAT(ITR,1)=INDXJ
          TRAT(ITR,2)=GAMMA
          TRAT(ITR,3)=ALPHA
          TRAT(ITR,4)=BETA
      ELSE
          K(INDXJ)=GAMMA
      END IF
      GO TO 10
 11   NREAC=J-1
      CLOSE(8)
C
C******************************************************************************
      IF(LRUN.NE.6) OPEN(LRUN,FILE=UNITL,STATUS='NEW')
      OPEN(2,FILE=UNIT2,STATUS='NEW') 
C
C** Total abundances for H2 and electrons
      TOTAL(1)=0.5
      TOTAL(2)=0.0
C--Set gas-phase abundances to XFRAC at start of calculations
      DO 25 I=1,NG
 25   Y0(I)=XFRAC
C--except H (set to HFRAC)
      Y0(IHI)=HFRAC
C--and ice mantle components, returned to the gas-phase
      XWATER=1.0D-4*FRM
C--H2O
      Y0(40)=1.0*XWATER
C--CO2
      Y0(49)=0.2*XWATER
C--CO
      Y0(11)=0.15*XWATER
C--CH4
      Y0(19)=0.04*XWATER
C--NH3
      Y0(28)=0.01*XWATER
C--H2CO
      Y0(47)=0.03*XWATER
C--CH3OH
C     Y0(??)=0.03*XWATER
C--He
      Y0(6)=0.1
C--------------------------------------------------------
C--Set solid state abundances to zero
      DO 27 I=NG+1,NTD
 27   Y0(I)=0.0
C
C--Set the density
      D=DEN0
C
      ISTATE=1
      IFC=0
      T0=T_INTERP(1)
      DG0=DG(1)
C       DG0=FDUST
      TEN0=TEN(1)
      NON0=NON(1)
C      TOUT=TMIN/TINC
C---------------------------------INTEGRATION LOOP-------------------------------
      IRUN=2
C 23      TOUT=TOUT*TINC
 23      TOUT=T_INTERP(IRUN)
C          WRITE(*,*) "TOUT=",TOUT
         DGOUT=DG(IRUN)         
C          WRITE(*,*) "DGOUT=",DGOUT 
         TENOUT=TEN(IRUN)
C          WRITE(*,*) "TENOUT=",TENOUT
         NONOUT=NON(IRUN)
C          WRITE(*,*) "NONOUT=",NONOUT
C
         CALL DLSODE(DIFFUN,N,Y0,T0,TOUT,ITOL,RTOL,ATOL,ITASK,
     *               ISTATE,IOPT,RWORK,LRW,IWORK,LIW,JAC,MF)
         IF(ISTATE.NE.2) WRITE(LRUN,201) ISTATE
C--Store the results
         IRUN=IRUN+1
         TAGE(IRUN)=TOUT/YEAR
         B(1,IRUN)=X(1)
         B(2,IRUN)=X(2)
         DO 45 I=1,NTD
            B(I+2,IRUN)=Y0(I)
 45      CONTINUE
         IF(ISTATE.NE.2) go to 46 
C---Call order if appropriate
         IF((IRUN.EQ.IT1).OR.(IRUN.EQ.IT2).OR.
     *      (IRUN.EQ.IT3).OR.(IRUN.EQ.IT4)) THEN
            RORD=0.0
            TORD=TOUT/YEAR
            CALL ORDER(NG,NTD,Y0,RORD,TORD)
         END IF
C
         IF(TOUT.LT.TMAX) GO TO 23
C
 46   CLOSE(2)
      CALL RESULT(IRUN,NTOT,IFULL,LOUT)
C      IF(IPLOT.EQ.1) CALL PLOT(IRUN,NTOT)
C      IF(IPLOT.EQ.0) CALL PLOT1(IRUN,NTOT)
      WRITE(LRUN,'('' Run completed.'')')
      CLOSE(LOUT)
C------------------------------------------------------------------------------
 201  FORMAT(1X,'ERROR: Integration failure. ISTATE = ',I2)
      END
C
C***************** FOR USE BY INTEGRATOR WHEN CALLING DIFFUN ******************
C
      SUBROUTINE ADJUST(N,T,Y,YDOT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NCONS=2,NRM=1500)
C
      REAL*8 K(NRM),Y(N),YDOT(N),X(NCONS),TOTAL(NCONS),
     *       TRAT(NRM,4),PRAT(NRM,4),DRAT(NRM,4),ERAT(NRM,4)
      COMMON/BLK3/X,TOTAL,K,D
      COMMON/BLK4/TRAT,PRAT,DRAT,ERAT,ITR,IPR,IDR,IER,IH2
      COMMON/BLK5/IFC,ICOV,ISURF
      COMMON/BLK6/DEN0,TEMP0,AV0
      COMMON/BLK8/GRAD,SAPH,YLD,YLD2,FSITES,G0,F0,FCR,YH,CRAT,
     *            STICK0,STICKP,STICKN
      COMMON/BLK9/FCO2,COCOV,FOX,PHOH,H2SHL,COSHL,FFRZ,FDES
      COMMON/BLK10/T0,TOUT,DG0,DGOUT,TEN0,TENOUT,NON0,NONOUT
C
      AV=AV0
      IF(IFC.EQ.0) THEN
C
C================================= D:G RATIO ==============================
C
C         WRITE(*,*) "DGOUT =",DGOUT
         SAPH=DGOUT*8.0E-21
         SBIND=1.0E15
         FSITES=SAPH*SBIND
C
C================================== DENSITY ===============================
C
         D=NONOUT
C
C================================ PHOTORATES ==============================
C
C--recalculate the photorates
         DO 5 I=1,IPR
             INDXJ=PRAT(I,1)
             K(INDXJ)=G0*PRAT(I,2)*DEXP(-AV*PRAT(I,3))+PRAT(I,4)
 5       CONTINUE
C--H2 and CO self-shielding factors
         K(1059)=H2SHL*K(1059)
         K(1066)=COSHL*K(1066)
C========= TEMPERATURE-DEPENDENT RATES & GAS-GRAIN INTERACTIONS ==============
C
C         TGAS=TEMP0
         TGAS=TENOUT
C Re-calculate temperature-dependent gas-phase rates
         TINV=1.0/TGAS
         T300=TGAS/300.0
         DO 1 J=1,ITR
            INDXJ=TRAT(J,1)
            K(INDXJ)=TRAT(J,2)*(T300**TRAT(J,3))*EXP(-TRAT(J,4)*TINV)
 1       CONTINUE
C Re-calculate freeze-out rate coefficients
         PREFAC=3.637E3*SAPH*DSQRT(TGAS)
         POSFAC=1.0+(16.711/(GRAD*TGAS))
         NEGFAC=DEXP(-16.711/(GRAD*TGAS))
C----H2 formation (Buch and Zhang(1991), ApJ 379,647)
         IF(TGAS.LE.300.D0) THEN
            STICKH=((TGAS/102.D0)+1.0)**(-2.0)
            GRATH=PREFAC*STICKH
         ELSE
            GRATH=0.0
         END IF
         PREFAC=FFRZ*PREFAC
         GRAT0=PREFAC*STICK0
         GRATP=PREFAC*POSFAC*STICKP
         GRATN=PREFAC*NEGFAC*STICKN
C
         DO 2 I=1,IDR
            IND=DRAT(I,1)
            GAMMA=DRAT(I,2)
            ALPHA=DRAT(I,3)
            BETA=DRAT(I,4)
C
            K(IND)=0.0
            IF(ALPHA.EQ.0.0) K(IND)=GRAT0/DSQRT(BETA)
            IF(ALPHA.EQ.-1.0) K(IND)=GRATN/DSQRT(BETA)
            IF(ALPHA.EQ.1.0) K(IND)=GRATP/DSQRT(BETA)
C--Apply freeze-out rules (to switch off rates)
C           IF((AV.LT.AVCRIT).AND.(GAMMA.EQ.4.0)) K(IND)=0.0
C           IF((AV.GT.AVCRIT).AND.(GAMMA.EQ.3.0)) K(IND)=0.0
            IF(GAMMA.EQ.3.0) K(IND)=0.0
 2       CONTINUE
C--H2 formation always on
         K(IH2)=GRATH
      END IF
      IFC=0
C
C  Fraction of binding sites occupied by GCO 
      COCOV=Y(81)/FSITES
      COCOV=DMIN1(COCOV,1.D0)
C  Fraction of binding sites occupied by GCO2 
      CO2COV=Y(88)/FSITES
      CO2COV=DMIN1(CO2COV,1.D0)
C  Fraction of binding sites occupied by GCH4 
      CH4COV=Y(80)/FSITES
      CH4COV=DMIN1(CH4COV,1.D0)
C--switch off direct CO->CO2 conversion on grains
      K(1133)=0.0
C--surface chemistry control:
      IF(ISURF.EQ.1) THEN
C--conversion of O,O+ and OH+ to CO2
         K(1231)=FCO2*COCOV*K(1231)
         K(1232)=FCO2*COCOV*K(1232)
         K(1233)=FCO2*COCOV*K(1233)
C--conversion of O,O+ and OH+ to H2O
         K(1143)=FOX*(1.0-(FCO2*COCOV))*K(1143)
         K(1145)=FOX*(1.0-(FCO2*COCOV))*K(1145)
         K(1210)=FOX*(1.0-(FCO2*COCOV))*K(1210)
C--conversion of O,O+ and OH+ to GO,GOH
         K(1226)=(1.0-FOX)*(1.0-(FCO2*COCOV))*K(1226)
         K(1228)=(1.0-FOX)*(1.0-(FCO2*COCOV))*K(1228)
         K(1227)=(1.0-FOX)*(1.0-(FCO2*COCOV))*K(1227)
      ELSE
         DO 99 ISC=1226,1233
            K(ISC)=0.0
 99      CONTINUE
      END IF
C
C  Desorption/evaporation reactions
C  Total surface coverage & net accretion rate (excluding negative terms)
      TOTS=0.0
      ACCR=0.0
      DO 7 I=80,96
         TOTS=TOTS+DMAX1(Y(I),0.D0)
         ACCR=ACCR+DMAX1(YDOT(I),0.D0)
 7    CONTINUE
      COV=TOTS/FSITES
C  H2 formation rate (=1/2 H-atom freeze-out rate)
      H2FORM=0.5*K(IH2)*Y(1)*D
C
      DO 6 I=1,IER
         IND=ERAT(I,1)
         GAMMA=ERAT(I,2)
         ALPHA=ERAT(I,3)
         BETA=ERAT(I,4)
C--cr desorption
C--enhance H2O cr desorption rate
C        IF(IND.EQ.1216) GAMMA=GAMMA*1.0
         CRDES=GAMMA*(CRAT/1.3E-17)
C--direct photodesorption + cr-induced photodesorption
         PDIR=SAPH*YLD*G0*F0*DEXP(-1.8*AV)
         PINDIR=SAPH*YLD2*FCR*(CRAT/1.3E-17)
         PHDES=PDIR+PINDIR
C--photodesorption yield for GO->O reduced by factor of 10
         IF(IND.EQ.1229) PHDES=0.1*PHDES
C--H2 formation-induced desorption
         H2DES=YH*H2FORM
C--correct for total surface coverage, when necessary
         IF(COV.LT.1.0) THEN
            PHDES=PHDES/FSITES
            H2DES=H2DES/FSITES
         ELSE
            IDSP=IDINT(ALPHA)
C--approximation 1: fractional surface coverage = x(i)/x(tot)
            PHDES=PHDES/TOTS
            H2DES=H2DES/TOTS
         END IF
C--inhibit surface desorption yields for H2O if CO or CO2 coverage is high
         IF((IND.EQ.1216).AND.(ICOV.EQ.1)) THEN
            CAP=DMIN1((COCOV+CO2COV+CH4COV),1.D0)
            PHDES=PHDES*(1.0-CAP)
C            PHDES=0.0
C            H2DES=H2DES*(1.0-CAP)
         END IF        
C--yield for GH2O->OH+H set to photodesorption yield for GH2O->H2O times 2
         IF(IND.EQ.1216) PHOH=2.0*FDES*PHDES
         K(IND)=FDES*(CRDES+PHDES+H2DES)
 6    CONTINUE
C
      RETURN
      END
C
C==============================================================================
C
C-- Subroutine to write out the results ---------------------------------------
      SUBROUTINE RESULT(IRUN,NMAX,IFULL,LOUT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NSP=100,NOP=2228)
C
      DIMENSION B(NSP,NOP),BL(NSP,NOP),TAGE(NOP)
      CHARACTER*8 SPECI(NSP)
      COMMON/BLK1/SPECI
      COMMON/BLK2/B,TAGE
C--List of selected species numbers (for reduced output)
C  [HCO+,H2CO,CS,N2H+,NH3,C,OH,ELECTR] +2 for conserved species
      DATA J1,J2,J3,J4,J5,J6,J7,J8/48,49,74,35,30,10,40,2/
C
C---Put the B and TAGE arrays into appropriate forms
      DO 1 I=1,IRUN
         TAGE(I)=DLOG10(TAGE(I))
C--Normalise the abundances to be relative to H2 [X(1)]
C        XH2=1.0/B(1,I)
         DO 2 J=1,NMAX
C           BTEM=B(J,I)*XH2
            BTEM=B(J,I)
            IF(BTEM.LE.0.0) BTEM=1.0
            BLOG=DLOG10(BTEM)
            BL(J,I)=DMAX1(BLOG,-99.99D0)
 2       CONTINUE
 1    CONTINUE
C
      IF(IFULL.EQ.1) THEN
C---Full output of all abundances
         IST=1
         IFN=8
 4       WRITE(LOUT,104) (SPECI(J),J=IST,IFN)
         WRITE(LOUT,101)
         DO 5 I=1,IRUN
            WRITE(LOUT,105) TAGE(I),(BL(I2,I),I2=IST,IFN)
 5       CONTINUE
         WRITE(LOUT,103)
         IST=IST+8
         IFN=IFN+8
         IF(IFN.GT.NMAX) IFN=NMAX
         IF(IST.LE.NMAX) GO TO 4
C---Reduced output
      ELSE IF(IFULL.EQ.0) THEN
       WRITE(LOUT,104) SPECI(J1),SPECI(J2),SPECI(J3),SPECI(J4),
     *   SPECI(J5),SPECI(J6),SPECI(J7),SPECI(J8)
       WRITE(LOUT,101)
       DO 6 I=1,IRUN
        WRITE(LOUT,105) TAGE(I),BL(J1,I),BL(J2,I),
     *   BL(J3,I),BL(J4,I),BL(J5,I),BL(J6,I),BL(J7,I),BL(J8,I)
 6     CONTINUE
       WRITE(LOUT,103)
      END IF
C
 101  FORMAT(1X,78('-'))
 103  FORMAT(//)
 104  FORMAT(1X,'Log(t/yr)',2X,8(A8))
 105  FORMAT(4X,0PF5.3,8(2X,0PF6.2))
C
      RETURN
      END
C
C******************************************************************************

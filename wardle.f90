PROGRAM  WARDLE
!
!------------------------------------------------------------------------------
!
! Single point calculation; model of time-dep. chemistry in a 
! static, dust-poor cloud
!
! Based on L1544/L1544.f [15th August 2013]
!
! H,He,C,N,O,S,Na chemistry (gas-phase, freeze-out and desorption)
!    Only H2 and Electrons are conserved.
!
! NG     is the number of time-dependent (gas-phase) species
! NCONS  is the number of conserved species
! NS     is the number of time-dependent (solid-state) species
!
! NSP    is >= the total number of species
! NRM    is >= the number of reactions
! NOP    is >= the number of output times
!
!------------------------------------------------------------------------------
!
      USE dvode_f90_m
      use collapse
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
      PARAMETER(NSP=214,NEL=8,NRM=2345,NOP=4000,NG=163,NS=51)
!
      CHARACTER*25  RE1(NRM),RE2(NRM),P1(NRM),P2(NRM),P3(NRM),P4(NRM), &
      &             SPECI(NSP),SPX(NEL)
      CHARACTER*25 UNIT2,UNIT3,UNIT4,UNIT7,UNIT8,UNITL,DGFILE,TENFILE,RHONFILE, &
      &             RHOEFILE, RHOIFILE
      DOUBLE PRECISION  Y0(NSP),FRAC(NEL),DPLT(NEL),ATOL(NSP), &
      & RATE(NRM),TRAT(NRM,4),PRAT(NRM,4),DRAT(NRM,4),ERAT(NRM,4),DEURAT(NRM,4), &
      & DESRAT(NRM,4),DESH2RAT(NRM,4),CRPHAT(NRM,4),B(NSP,NOP),TAGE(NOP), &
      & T_INTERP(NOP),DG(NOP),TEN(NOP),RHON(NOP),RHOE(NOP),RHOI(NOP),NON(NOP), &
      & NOE(NOP),NOI(NOP),NONOUT,DUMMY(NOP)
      DOUBLE PRECISION IONSUM
      INTEGER INDR(NRM),ISPX(NEL),SSTATE(NS),N
      INTEGER GRINHIB,STATIC
!
!** For LSODE
      CHARACTER*5 JAC
      EXTERNAL DIFFUN
!**
      COMMON/BLK1/SPECI
      COMMON/BLK2/B,TAGE
      COMMON/BLK3/RATE,D
      COMMON/BLK4/TRAT,PRAT,CRPHAT,DRAT,DEURAT,DESRAT,ITR,IPR,IDR,IDEU,IDES,IDESH2,ICR,IH2
      COMMON/BLK5/ICOV,ISURF
      COMMON/BLK6/DEN0,TEMP0,AV0
      COMMON/BLK8/GRAD,SAPH,YLD,YLD2,FSITES,G0,F0,FCR,YH,CRAT, &
     &            STICK0,STICKP,STICKN,ALBFAC,CRDEF
      COMMON/BLK9/FCO2,COCOV,FOX,PHOH,H2SHL,COSHL,FFRZ,FDES
      COMMON/BLK10/T0,TOUT,DG0,DGOUT,TEN0,TENOUT,TDUST,NONOUT
      COMMON/BLK12/COVLIM,TDES
      COMMON/ORDER1/RE1,RE2,P1,P2,P3,P4
      COMMON/ORDER2/INDR,IOTYP

      ! For DVODE
      TYPE(VODE_OPTS) :: OPTIONS
!
      ! DATA  NG,NS/165,49/
!-- FRAC, DPLT are the elemental abundances and depletions factors for -----
!  species (SPX) whose index number in the species set is given by ISPX
!  It is assumed that conserved species are (1) H2 and (2) electrons.
!  IHI AND ICO are the index numbers for H and CO respectively.
!---------------------------------------------------------------------------
      DATA  FRAC/0.1,3.6E-4,1.13E-4,6.72E-4,1.47E-5,2.68E-5,3.24E-5,3.16E-7/
!--Hollenbach et al. (2009) values
!     DATA  DPLT/1.0,0.373,1.0,0.475,1.73,0.086/
!--Asplund, Grevesse & Sauval values (2005)
      DATA  DPLT/1.0,0.68,0.53,0.68,0.85,1.0,1.0,1.0/
!     DATA  DPLT/1.0,0.5,1.0,1.0,0.01,0.086/
      DATA  SPX/'HE','C','N','O','S','MG','SI','CL'/
      DATA  ISPX/7,10,16,27,103,42,71,117/
      DATA  IHI,ICO/1,65/
!
      DATA PC,YEAR/3.0856D18,3.1557D7/
      DATA BOLTK,AMU,PI/1.381D-23,1.67D-27,3.1415927/

      DATA COVLIM,TDES/0.01,120/
!
!**
!
!--NTD is the total number of time-dependent chemical species
      NTD=NG+NS
!--N is the total number of time-dependent species
      N=NTD
!--NTOT is the total number of chemical species
      NTOT=NTD
!---------------------------- INPUT PARAMETERS --------------------------------
!----I/O files
      LOUT=3
      STOUT=4
      LRUN=6
      UNIT2='ordout.d'
      UNIT3='output.d'
      UNITL='report.d'
      UNIT7='coll_specs.d'
      UNIT8='coll_rates.d'
      DGFILE='interp/dg_interp.dat'
      TENFILE='interp/ten_interp.dat'
      RHONFILE='interp/rhon_interp.dat'
      RHOEFILE='interp/rhoe_interp.dat'
      RHOIFILE='interp/rhoi_interp.dat'
!
!----Set GRINHIB to default value (i.e. do not inhibit desorption)
      GRINHIB=0

!----Set STATIC to default value (i.e. static evolution)
      STATIC=0
!
!------------------------------------------------------------------------------
!
!    Read the data
      OPEN(20,FILE=DGFILE,STATUS='OLD')
      OPEN(21,FILE=TENFILE,STATUS='OLD')
      OPEN(22,FILE=RHONFILE,STATUS='OLD')
      OPEN(23,FILE=RHOEFILE,STATUS='OLD')
      OPEN(24,FILE=RHOIFILE,STATUS='OLD')

      DO i=1,NOP
            READ(20,*) T_INTERP(i),DG(i)
            READ(21,*) DUMMY(i),TEN(i)
            READ(22,*) DUMMY(i),RHON(i)
            READ(23,*) DUMMY(i),RHOE(i)
            READ(24,*) DUMMY(i),RHOI(i)

!          Convert mass density to number density
!          Note: this is the density of all of the neutral
!          species (alongside electron and ions)
!          n = rho/(mu*mh)
           NON(i)=RHON(i)/(2*1.67D-24)
           NOE(i)=RHOE(i)/(9.109D-28)
           NOI(i)=RHOI(i)/(2*1.67D-24)
      END DO
!
!------------------------------------------------------------------------------
!--IFULL=0,1 for reduced/full output of abundances 
      IFULL=1
!--Density, temperature and extinction
      DEN0=NON(1)
      TEMP0=TEN(1)
      AV0=2.0
!--Fraction of dust present
      FDUST=DG(1)
!--Metallicity fraction
      FRM=100.0
!
      XFRAC=1.0D-06
      HFRAC=1.0/DEN0
!
!--Maximum time
      TMAX=T_INTERP(NOP)
!--Photochemistry parameters
      CRAT=1.3E-17
      DEFF=200.0
      ALBEDO=0.5
!--H2 and CO shielding factors
      ! 1.0 is no shielding; 0.0 is total shielding
      H2SHL=1.0
      COSHL=1.0
!
!--Freeze-out parameters: GRAD = grain radius (in microns), SAPH = grain 
! surface area per H-nucleon (in cm2).
! STICK0,STICKP,STICKN are the sticking efficiencies for neutrals,+ve,-ve ions
      GRAD=0.0083
      SAPH=8.0E-21*FDUST
      STICK0=1.0
      STICKP=1.0
      STICKN=1.0
!--Surface chemistry parameters: 
! FCO2 is the proportion of O,O+,OH (+CO ice) that converts to CO2 
    !   FCO2=0.1
      FCO2=0.0
! FOX is the fraction of the remainder O,O+,OH that converts to H2O
    !   FOX=1.0
      FOX=1.0
!--Desorption parameters: 
! YLD = (direct) photodesorption yield, 
! YLD2= (CR-induced) photodesorption yield,
!     YLD/YLD2 doubled for GH2O->OH+H, reduced by 10 for GO->O
! G0 = flux factor (multiples of Draine), 
! F0 = I/S photon flux, 
! FCR = CR-induced photon flux (for zeta_0),
! SBIND = surface density of binding sites, 
! CRDEF = default value of the CR-heating desorption rate, if not known,
! YH = desorption yield due to H2 formation,
! ICOV = 1/0 switches on/off H2O desorption protection by CO/CO2/CH4.
!
      YLD=1.0E-3
      YLD2=1.0E-3
      G0=1.0
      F0=1.0E8
      FCR=4875.0
      SBIND=1.0E15
      CRDEF=1.0E-17
      YH=0.0
      ICOV=0
!--FFRZ,FDES are flags (on:1.0, off:0.00) for freeze-out & desorption
! ISURF controls the surface chemistry (1:on, 0:off)
      FFRZ=1.0
      FDES=1.0
      ISURF=1
!-------------------- END OF INPUT PARAMETERS ---------------------------------
      OPEN(LOUT,FILE=UNIT3,STATUS='UNKNOWN')
!--Write out the input parameters to the log file
      WRITE(LOUT,100)
 100  FORMAT(1X,'INPUT PARAMETERS:',/)
      DO 101 I=1,NEL
 101    WRITE(LOUT,102) SPX(I),FRAC(I),DPLT(I)
 102  FORMAT(1X,A7,'/H = ',1PE8.2,'*',1PE8.2)
      WRITE(LOUT,108) DEN0,TEMP0,AV0
 108  FORMAT(/,1X,'Density = ',1PE8.2,' cm-3',/, &
     &   1X,'Temperature = ',0PF5.1,' K',/, &
     &   1X,'Extinction = ',1PE8.2,' magnitudes',/)
      WRITE(LOUT,103) ALBEDO,CRAT,DEFF,H2SHL,COSHL
 103  FORMAT(1X,'Albedo = ',1PE8.2,/, &
     &   1X,'Cosmic ray ionization rate = ',1PE8.2,' s-1',/, &
     &   1X,'Default CRP-induced photolysis rate = ',1PE8.2,/, &
     &   1X,'H2 self-shielding factor = ',0PF5.3,/, &
     &   1X,'CO self-shielding factor = ',0PF5.3,/)
      WRITE(LOUT,104) FDUST,FRM,GRAD,SAPH,STICK0,STICKP,STICKN, &
     &                FFRZ,FDES,ISURF
 104  FORMAT(1X,'Fraction of dust present = ',0PF5.3,/, &
     &   1X,'Metallicity fraction = ',0PF6.2,/, &
     &   1X,'Grain radius = ',1PE8.2,' microns',/, &
     &   1X,'Eff. dust surface area per H = ',1PE8.2,' cm-2',/, &
     &   1X,'Sticking coefficient (neutrals) = ',1PE8.2,/, &
     &   1X,'Sticking coefficient (positive ions) = ',1PE8.2,/, &
     &   1X,'Sticking coefficient (negative ions) = ',1PE8.2,/, &
     &   1X,'Freeze-out efficiency = ',0PF5.3,/, &
     &   1X,'Desorption efficiency = ',0PF5.3,/, &
     &   1X,'Surface chemistry flag = ',I1,/)
      WRITE(LOUT,105) FCO2,FOX,ICOV,YLD,YLD2,G0,F0,FCR,SBIND,YH
 105  FORMAT(1X,'G(O/OH)+GCO->GCO2 conversion efficiency = ',1PE8.2,/, &
     &   1X,'G(O/OH)[rem.]->GH2O conversion efficiency = ',1PE8.2,/, &
     &   1X,'H2O desorp. shielding by CO/CO2/CH4 (1,0;on,off) = ',I1,/, &
     &   1X,'(Direct) photodesorption yield = ',1PE8.2,/, &
     &   1X,'(CR-induced) photodesorption yield = ',1PE8.2,/, &
     &   1X,'ISRF/Draine = ',1PE8.2,/, &
     &   1X,'ISRF (flux) = ',1PE8.2,' photons cm-2 s-1',/, &
     &   1X,'CR-induced flux (norm.) = ',1PE8.2,' photons cm-2 s-1',/, &
     &   1X,'Binding site surface density = ',1PE8.2,' cm-2',/, &
     &   1X,'H2 formation desorption yield = ',1PE8.2,/)
      WRITE(LOUT,106) IFULL
 106  FORMAT(1X,'Full/Partial results output flag = ',I1,/)
      WRITE(LOUT,107)
 107  FORMAT(1X,75('-'))
!
!------------------------------------------------------------------------------
!
! Read the species file
      OPEN(7,FILE=UNIT7,STATUS='OLD')
      READ(7,1)(INDXS,SPECI(J),J=1,NTD)
      CLOSE(7)
 1    FORMAT((3X,I3,3X,A8,1X))
!
!------------------------------------------------------------------------------
!
! Read reaction data
      OPEN(8,FILE=UNIT8,STATUS='OLD')
      ITR=0
      IPR=0
      IDR=0
      IER=0
      ICR=0
!
      ALBFAC=1.0/(1.0-ALBEDO)
      FSITES=SAPH*SBIND
!
      J=0
 10   J=J+1
      READ(8,3) INDR(J),RE1(J),RE2(J),P1(J),P2(J),P3(J),P4(J), &
     &          ALPHA,BETA,GAMA
  3   FORMAT(I4,4(1X,A8),2(1X,A4),1X,1PE10.2,1X,0PF6.2,1X,F8.1)
      INDXJ=INDR(J)
!
! store reaction data
      BETA = 0
      IF (INDXJ.EQ.9999) GO TO 11
      IF (RE2(J).EQ.'CRP') THEN ! CRP is cosmic ray proton (i.e. H-nucleus)
          RATE(INDXJ)=ALPHA
      ELSE IF (RE2(J).EQ.'CRPHOT') THEN ! CRPHOT is cosmic ray induced photon emission
          ICR=ICR+1
      !     IF(BETA.EQ.0.0) BETA=DEFF
          CRPHAT(ICR,1)=INDXJ
          CRPHAT(ICR,2)=ALPHA
          CRPHAT(ICR,3)=BETA
          CRPHAT(ICR,4)=GAMA
      ELSE IF (RE2(J) .EQ. 'PHOTON') THEN ! PHOTON is for standard photoreactions
          IPR=IPR+1
          PRAT(IPR,1)=INDXJ
          PRAT(IPR,2)=ALPHA
          PRAT(IPR,3)=BETA
          PRAT(IPR,4)=GAMA
      ELSE IF (RE2(J) .EQ. 'DES') THEN ! Desorption reactions
          IER=IER+1
          ERAT(IER,1)=INDXJ
          ERAT(IER,2)=ALPHA
          ERAT(IER,3)=BETA
          ERAT(IER,4)=GAMA
      ELSE IF (RE2(J).EQ.'FREEZE') THEN ! Freeze out
          IDR=IDR+1
          DRAT(IDR,1)=INDXJ
          DRAT(IDR,2)=ALPHA
          DRAT(IDR,3)=BETA
          DRAT(IDR,4)=GAMA
          !--Identify reaction number for: H + FREEZE -> H2
          IF (RE1(J).EQ.'H') IH2=INDXJ
      ELSE
          ! These are the collisional dissociation reactions
          ITR=ITR+1
          TRAT(ITR,1)=INDXJ
          TRAT(ITR,2)=ALPHA
          TRAT(ITR,3)=BETA
          TRAT(ITR,4)=GAMA
      END IF
      GO TO 10
 11   NREAC=J-1
      CLOSE(8)
!
!******************************************************************************
      IF(LRUN.NE.6) OPEN(LRUN,FILE=UNITL,STATUS='NEW')
      OPEN(2,FILE=UNIT2,STATUS='UNKNOWN') 
!
      IONSUM = 0.0
!** Total abundances for H
      fh = 0.0
!--- Set abundances
      N = 1
      INDX = 0
      DO N = 1, NTOT
        ! Set everything to 0.0 initially
        Y0(N) = 0.0

        !-- Set elemental abundances where required
        DO M=1,NEL
          IF (N.EQ.ISPX(M)) THEN
                ! Set the elementral abundance to the
                ! fractional abundance * depletion factor
                Y0(N) = FRAC(M)*DPLT(M)
          END IF
        END DO

        !-- Tweak further abundances
        IF (SPECI(N) == "H") THEN
            Y0(N) = (0.5*(1.0e0-fh)) 
        ELSE IF (SPECI(N) == "H2") THEN
            Y0(N) = 0.5*(0.5*(1.0e0-fh))
        ELSE IF (SPECI(N) == "SI+") THEN
            Y0(72) = FRAC(7)*DPLT(7)
            Y0(71) = 0.0 ! This sets the Y0(Si) = 0
        ELSE IF (SPECI(N) == "CL+") THEN
            Y0(118) = FRAC(8)*DPLT(8)
            Y0(117) = 0.0
        END IF

        !-- Finds all ions and sums their abundance
        DO J=1,25
            IF ((SPECI(N)(J:J) .EQ. "+")) THEN
                IONSUM = IONSUM + Y0(N)
                WRITE(*,*) "IONSUM=",IONSUM
            END IF
        END DO
        Y0(214) = IONSUM ! Set the initial abundance of ions and electrons to be equal
      END DO
! 
!--------------------------------------------------------
!
!--Set the density
      D=DEN0
      AV = ((1.25D-15*1.6D21)/DEN0)
      T0 = 0.0
      DG0=DG(1)
      TEN0=TEN(1)
!
      ISTATE=1
      ITERP=0

!** For LSODE
      DATA ITOL,MF/1,10/
      JAC='DUMMY'
      ITASK=1
      IOPT=1
      IOTYP=0
      RTOL=1.0E-4
      ATOL=1.0E-14*Y0
      WHERE(ATOL<1E-30) ATOL=1E-30 ! to a minimum degree
      MXSTEP=1.0E6
      OPTIONS = SET_OPTS(METHOD_FLAG=MF, ABSERR_VECTOR=ATOL, RELERR=RTOL, USER_SUPPLIED_JACOBIAN=.FALSE., MXSTEP=MXSTEP)

!---------------------------------INTEGRATION LOOP-------------------------------
      IRUN=0
      STATIC=1
      TOUT = 1.0d-8
 23      IRUN=IRUN+1
         ! Set the physical conditions to variables for the integrator and
         ! adjust subroutine
         IF ((TOUT/YEAR .LT. 1e6) .AND. (STATIC .EQ. 1)) THEN
            IRUN = 0   
            ! TOUT = TOUT*10
            IF (TOUT/YEAR .gt. 1.0d6) THEN !code in years for readability, TOUT in s
                  TOUT=(TOUT/YEAR+1.0d4)*YEAR
            ELSE  IF (TOUT/YEAR .gt. 1.0d5) THEN
                  TOUT=(TOUT/YEAR+1.0d4)*YEAR
            ELSE IF (TOUT/YEAR .gt. 1.0d4) THEN
                  TOUT=(TOUT/YEAR+1000.0)*YEAR
            ELSE IF (TOUT/YEAR .gt. 1000) THEN
                  TOUT=(TOUT/YEAR+100.0)*YEAR
            ELSE IF (TOUT/YEAR .gt. 0.0) THEN
                  TOUT=(TOUT/YEAR*10.0)*YEAR
            ELSE
                  TOUT=3.16d7*10.d-8
            ENDIF
            DGOUT = DG(1)
            TENOUT = TEN(1)
            TDUST = TEN(1)
            NONOUT = NON(1)
            ISTATE=1
            AV = ((1.25D-15*1.6D21)/NONOUT)
            WRITE(*,*) "Y0=",Y0
            WRITE(*,*) "TOUT=",TOUT/YEAR, irun
            !WRITE(*,*) "TENOUT=",TENOUT
         ELSE
            ! Otherwise use data
            TOUT = T_INTERP(IRUN)
            DGOUT = DG(IRUN)
            TENOUT = TEN(IRUN)
            TDUST = TEN(IRUN)
            NONOUT = NON(IRUN)
            AV = ((1.25D-15*1.6D21)/NONOUT)

            ! Reset the correct parameters to ensure integrator can work
            IF (IRUN .EQ. 1) THEN
                T0 = 0.0
                ISTATE = 1
                STATIC = 0

                ! Data is in the form of the neutral and ion densities, i.e.
                ! the sum of the abundances of all of the neutrals and all of the
                ! ions. These lines scale the abundances so as the sum of all neutrals
                ! and sum of all ions is consistent with the datafiles.
                N = 1
                INDX = 0
                DO N = 1, NTOT
                    IF (SPECI(N)(1:1) .NE. "#") THEN
                        DO M=1,25
                            ! Scale the ion abundances
                            IF ((SPECI(N)(M:M) .EQ. "+")) THEN
                                Y0(N) = (Y0(N)/IONSUM) * NOI(IRUN)
                            ELSE
                                Y0(N) = 0.0 ! From UCLCHEM's chemistry.f90
                            END IF
                        END DO
                    END IF
                END DO
            END IF
         END IF

        !  WRITE(*,*) "TOUT=",TOUT
        !  WRITE(*,*) "TENOUT=",TENOUT
        !  WRITE(*,*) "NONOUT=",NONOUT

         CALL DVODE_F90(DIFFUN, NTD, Y0, T0, TOUT, ITASK, ISTATE, OPTIONS)
         IF (ISTATE.EQ.-1) CALL RESULT(IRUN,NTOT,IFULL,LOUT)
         IF(ISTATE.NE.2) WRITE(LRUN,201) ISTATE
!--Store the results
         TAGE(IRUN)=TOUT/YEAR
         DO 45 I=1,NTD
            B(I,IRUN)=Y0(I)
 45      CONTINUE
         IF(ISTATE.NE.2) go to 46 
!---Call order if appropriate
         IF ((TOUT/YEAR .gt. 4.6e2) .AND. (TOUT/YEAR .lt. 6e2)) THEN
            RORD=0.0
            TORD=TOUT/YEAR
            CALL ORDER(NG,NTD,Y0,RORD,TORD)
         END IF
!
         IF(STATIC .EQ. 1) THEN
            IRUN = 0
            GO TO 23
         END IF
         IF(TOUT.LT.TMAX) GO TO 23
!
 46   CLOSE(2)
      CALL RESULT(IRUN,NTOT,IFULL,LOUT)
!     IF(IPLOT.EQ.1) CALL PLOT(IRUN,NTOT)
!     IF(IPLOT.EQ.0) CALL PLOT1(IRUN,NTOT)
      WRITE(LRUN,'('' Run completed.'')')
      CLOSE(LOUT)
!------------------------------------------------------------------------------
 201  FORMAT(1X,'ERROR: Integration failure. ISTATE = ',I2)
      END PROGRAM WARDLE
!
!***************** FOR USE BY INTEGRATOR WHEN CALLING DIFFUN ******************
!
       SUBROUTINE ADJUST(N,Y,YDOT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NRM=2345,NG=163,NS=51)
!
      DOUBLE PRECISION RATE(NRM),Y(N),YDOT(N), &
     &       TRAT(NRM,4),PRAT(NRM,4),DRAT(NRM,4),DEURAT(NRM,4),DESRAT(NRM,4),  &
     &       CRPHAT(NRM,4)
      INTEGER SSTATE(NS)
      DOUBLE PRECISION NONOUT,MANTLE
      INTEGER GRINHIB
      COMMON/BLK3/RATE,D
      COMMON/BLK4/TRAT,PRAT,CRPHAT,DRAT,DEURAT,DESRAT,ITR,IPR,IDR,IDEU,IDES,IDESH2,ICR,IH2
      COMMON/BLK5/ICOV,ISURF
      COMMON/BLK6/DEN0,TEMP0,AV0
      COMMON/BLK8/GRAD,SAPH,YLD,YLD2,FSITES,G0,F0,FCR,YH,CRAT, &
     &            STICK0,STICKP,STICKN,ALBFAC,CRDEF
      COMMON/BLK9/FCO2,COCOV,FOX,PHOH,H2SHL,COSHL,FFRZ,FDES
      COMMON/BLK10/T0,TOUT,DG0,DGOUT,TEN0,TENOUT,TDUST,NONOUT
      COMMON/BLK11/SSTATE
      COMMON/BLK12/COVLIM,TDES
!
!
!================================= D:G RATIO ==============================
!
      SAPH=DGOUT*8.0E-21 ! Grain surface area per H nucleon (in cm2)
      SBIND=1.0E15
      FSITES=SAPH*SBIND ! Number of binding sites per unit volume
!
!================================== DENSITY ===============================
!
      D=NONOUT
      AV=((1.25D-15*1.6D21)/D)
!
!========================= GRAIN CHEMISTRY INHIBITION =====================
!
      IF ((COV .LT. COVLIM) .AND. (TDUST .GT. TDES)) THEN
          GRINHIB = 1
      ELSE IF ((COV .GT. COVLIM) .AND. (TDUST .LT. TDES)) THEN
          GRINHIB = 0
      END IF
!
!================================ PHOTORATES ==============================
!
!--recalculate the photorates
      DO 5 I=1,IPR
            INDXJ=PRAT(I,1)
            RATE(INDXJ)=G0*PRAT(I,2)*DEXP(-AV*PRAT(I,4))
5     CONTINUE
!--H2 and CO self-shielding factors
      RATE(359)=H2SHL*RATE(359) ! H2 is dissociated by CRP
      RATE(360)=H2SHL*RATE(360) 
      RATE(361)=H2SHL*RATE(361) 

      RATE(2001)=COSHL*RATE(2001)
!========= TEMPERATURE-DEPENDENT RATES & GAS-GRAIN INTERACTIONS ==============
      TGAS=TENOUT
!Re-calculate temperature-dependent gas-phase rates
      TINV=1.0/TGAS
      T300=TGAS/300.0
      DO 1 J=1,ITR
          INDXJ=TRAT(J,1)
          RATE(INDXJ)=TRAT(J,2)*(T300**TRAT(J,3))*EXP(-TRAT(J,4)*TINV)
1     CONTINUE
!Calculate the cosmic ray induced photon rates
      DO 3 J=1,ICR
          INDXJ=CRPHAT(J,1)
          RATE(INDXJ)=CRPHAT(J,2)*CRPHAT(J,4)*(T300)**CRPHAT(J,3)
3     CONTINUE
!Re-calculate freeze-out rate coefficients
      PREFAC=3.637E3*SAPH*DSQRT(TGAS)
      POSFAC=1.0+(16.711/(GRAD*TGAS))
      NEGFAC=EXP(-16.711/(GRAD*TGAS))

!----H2 formation (Buch and Zhang(1991), ApJ 379,647)
      IF(TGAS .LE. 300) THEN
          STICKH=((TGAS/102)+1.0)**(-2.0) ! Sticking probability
          GRATH=PREFAC*STICKH
      ELSE
          GRATH=0.0
      END IF

!----Freeze out reactions (FFRZ=0 for no freeze-out)
      PREFAC=FFRZ*PREFAC
      GRAT0=PREFAC*STICK0
      GRATP=PREFAC*POSFAC*STICKP
      GRATN=PREFAC*NEGFAC*STICKN
!
      DO 2 I=1,IDR
      IND=DRAT(I,1)
      ALPHA=DRAT(I,2)
      BETA=DRAT(I,3)
      GAMA=DRAT(I,4)
!
      RATE(IND)=0.0
      IF(BETA.EQ.0.0) THEN 
            RATE(IND)=GRAT0/DSQRT(GAMA)
      ELSE IF (BETA.EQ.-1.0) THEN 
            RATE(IND)=GRATN/DSQRT(GAMA)
      ELSE IF(BETA.EQ.1.0) THEN
            RATE(IND)=GRATP/DSQRT(GAMA)
      ELSE IF (BETA.EQ.2.0) THEN ! This is for electron freeze-out (not really freeze-out, just grain-e- attraction)
            CION=1.0+16.71d-4/(GRAD*TENOUT)
            RATE(IND)=4.57d4*GAMA*SAPH*CION
      END IF
!--Apply freeze-out rules (to switch off rates)
      IF (GRINHIB .EQ. 1) RATE(IND)=0.0 
2 CONTINUE
!--H2 formation is dependent on GRATH (see above)
      RATE(IH2)=GRATH
      
! Fraction of binding sites occupied by #CO 
      COCOV=Y(61)/FSITES
      COCOV=DMIN1(COCOV,1.D0)
! Fraction of binding sites occupied by #CO2 
      CO2COV=Y(148)/FSITES
      CO2COV=DMIN1(CO2COV,1.D0)
! Fraction of binding sites occupied by #CH4 
      CH4COV=Y(22)/FSITES
      CH4COV=DMIN1(CH4COV,1.D0)
!--surface chemistry control:
      IF(ISURF.EQ.1) THEN
!--conversion of CO2 and HCO2+ to #CO2
         RATE(2207)=FCO2*COCOV*RATE(2207)
      !    RATE(2239)=FCO2*COCOV*RATE(2239)
!--conversion of OH, H2O, O+, OH+, H2O+, H3O+ and O to #H2O
         RATE(2203)=FOX*(1.0-(FCO2*COCOV))*RATE(2203)
         RATE(2206)=FOX*(1.0-(FCO2*COCOV))*RATE(2206)
         RATE(2220)=FOX*(1.0-(FCO2*COCOV))*RATE(2219)
         RATE(2226)=FOX*(1.0-(FCO2*COCOV))*RATE(2225)
         RATE(2232)=FOX*(1.0-(FCO2*COCOV))*RATE(2225)
         RATE(2239)=FOX*(1.0-(FCO2*COCOV))*RATE(2238)
         RATE(2245)=FOX*(1.0-(FCO2*COCOV))*RATE(2244)
      ELSE
         ! This disables the surface chemistry if ISURF != 1
         DO 99 ISC=1,NS
            INDX = SSTATE(ISC)
            RATE(INDX)=0.0
 99      CONTINUE
      END IF
!
! Desorption/evaporation reactions
! Total surface coverage & net accretion rate (excluding negative terms)
      TOTS=0.0
      ACCR=0.0
      DO 7 I=1,NS
         TOTS=TOTS+Y(SSTATE(I))
         ACCR=ACCR+YDOT(SSTATE(I))
 7    CONTINUE
      COV=TOTS/FSITES
      ! WRITE(*,*) "==============="
      ! WRITE(*,*) "TOTS=",TOTS
      ! WRITE(*,*) "FSITES=",FSITES
      ! WRITE(*,*) "COV=",COV
      ! WRITE(*,*) "==============="
 
! H2 formation rate (=1/2 H-atom freeze-out rate)
      H2FORM=0.5*RATE(IH2)*Y(1)*D
      ! H2FORM=1.0d-17*DSQRT(TGAS)

      DO 8 I=1,IER
!--Thermal desorption
          IND=ERAT(I,1)
          ALPHA=ERAT(I,2)
          ADMASS=ERAT(I,3)
          TBIND=ERAT(I,4)

          EFAC=-TBIND/TDUST
          EFAC=DMAX1(EFAC,-100.D0)
          EFAC=DMIN1(EFAC,-40.D0)
          FREQ0=4092.0*DSQRT(SBIND*TBIND/ADMASS)
          THDES=FREQ0*DEXP(EFAC)*SAPH*SBIND

!--cr desorption
          CRDES=ALPHA*CRAT*ALBFAC

!--cr-induced photodesorption
          PINDIR=SAPH*YLD2*FCR*CRDES
          PHDES=PINDIR

 ! Turn off C ionisation
      IF (IND.EQ.356) RATE(INDXJ) = 0.0
      IF (IND.EQ.379) RATE(INDXJ) = 0.0
      IF (IND.EQ.1973) RATE(INDXJ) = 0.0

!--photodesorption yield for #O2->O2 reduced by factor of 10
         IF(IND.EQ.2313) PHDES=0.1*PHDES
      !    IF(IND.EQ.2330) PHDES=0.1*PHDES
!--H2 formation-induced desorption
          H2DES=YH*H2FORM
!--correct for total surface coverage, when necessary
          IF(COV.LT.1.0) THEN
              PHDES=PHDES/FSITES
              H2DES=H2DES/FSITES
          ELSE
              IDSP=IDINT(BETA)
!--approximation 1: fractional surface coverage = x(i)/x(tot)
              PHDES=PHDES/TOTS
              H2DES=H2DES/TOTS
          END IF
!--inhibit surface desorption yields for H2O if CO or CO2 coverage is high
         IF((IND.EQ.2298.OR.IND.EQ.2299).AND.(ICOV.EQ.1)) THEN
            CAP=DMIN1((COCOV+CO2COV+CH4COV),1.D0)
            PHDES=PHDES*(1.0-CAP)
          PHDES=0.0
          H2DES=H2DES*(1.0-CAP)
        END IF        
!--yield for GH2O->OH+H set to photodesorption yield for GH2O->H2O times 2
         IF(IND.EQ.401.OR.IND.EQ.2015) PHOH=2.0*FDES*PHDES
!--Determine the rate for all desorption processes
         IF (GRINHIB .EQ. 0) THEN 
            RATE(IND)=FDES*(THDES+CRDES+PHDES+H2DES)
         ELSE
            RATE(IND)=0.0
         END IF
8 CONTINUE
!
      RETURN
      END SUBROUTINE ADJUST
!
!---------------------- For use with DLSODE--------------------------------------
!
      SUBROUTINE DIFFUN(N,T,Y,YDOT)

        USE collapse
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DOUBLE PRECISION PROD,LOSS
        DOUBLE PRECISION RATE
!
        PARAMETER(NRM=2345,NSP=214,NG=163,NS=51)
!
        INTEGER SSTATE(NS)
        DIMENSION RATE(NRM),Y(N),YDOT(N)
        COMMON /BLK3/RATE,D
        COMMON /BLK9/FCO2,COCOV,FOX,PHOH,H2SHL,COSHL,FFRZ,FDES
        COMMON /BLK11/SSTATE
!
      !   IF (STATIC .eq. 1) YDOT(NEQ)=densdot(Y(NEQ))

        CALL ADJUST(N,Y,YDOT)
        INCLUDE 'odes.f90'

        !H2 formation should occur at both steps - however note that here there is no 
        !temperature dependence. y(nh) is hydrogen fractional abundance.
        nh=1 
        nh2=3
        YDOT(nh)  = ydot(nh) - 2.0*( h2form*y(nh)*D - h2dis*y(nh2) )
        !                             h2 formation - h2-photodissociation
        YDOT(nh2) = ydot(nh2) + h2form*y(nh)*D - h2dis*y(nh2)
        !                       h2 formation  - h2-photodissociation

        RETURN
      END SUBROUTINE DIFFUN
!
!-- Subroutine to write out the results ---------------------------------------
      SUBROUTINE RESULT(IRUN,NMAX,IFULL,LOUT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NSP=214,NOP=4000)
!
      DIMENSION B(NSP,NOP),BL(NSP,NOP),TAGE(NOP)
      CHARACTER*25 SPECI(NSP)
      COMMON/BLK1/SPECI
      COMMON/BLK2/B,TAGE
!--List of selected species numbers (for reduced output)
! [HCO,HCO+,H2CO,CH2,H2,O,#H2CO,E-] 
      DATA J1,J2,J3,J4,J5,J6,J7,J8/75,76,84,14,3,27,83,214/
!
!---Put the B and TAGE arrays into appropriate forms
      DO 1 I=1,IRUN
         IF (TAGE(I).EQ.0.0) THEN
            TAGE(I) = 0.0
         ELSE
            TAGE(I)=LOG10(TAGE(I))
         END IF
!--Normalise the abundances to be relative to H2 [X(1)]
      !    XH2=1.0/B(1,I)
         DO 2 J=1,NMAX
            ! BTEM=B(J,I)*XH2
            BTEM=B(J,I)
            IF(BTEM.LE.0.0) BTEM=1.0
            BLOG=DLOG10(BTEM)
            BL(J,I)=DMAX1(BLOG,-99.99D0)
 2       CONTINUE
 1    CONTINUE
!
      IF(IFULL.EQ.1) THEN
!---Full output of all abundances
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
!---Reduced output
      ELSE IF(IFULL.EQ.0) THEN
       WRITE(LOUT,104) SPECI(J1),SPECI(J2),SPECI(J3),SPECI(J4), &
     &   SPECI(J5),SPECI(J6),SPECI(J7),SPECI(J8)
       WRITE(LOUT,101)
       DO 6 I=1,IRUN
        WRITE(LOUT,105) TAGE(I),BL(J1,I),BL(J2,I), &
     &   BL(J3,I),BL(J4,I),BL(J5,I),BL(J6,I),BL(J7,I),BL(J8,I)
 6     CONTINUE
       WRITE(LOUT,103)
      END IF
!
 101  FORMAT(1X,78('-'))
 103  FORMAT(//)
 104  FORMAT(1X,'Log(t/yr)',2X,8(A8))
 105  FORMAT(4X,0PF5.3,8(2X,0PF6.2))
!
      RETURN
      END SUBROUTINE RESULT
!
!******************************************************************************
!
      SUBROUTINE ORDER(NG,NTD,Y0,RAD,TIME)
!
!------------------------------------------------------------------------------
!  Subroutine to list out the formation and destruction reactions for each 
!  species in order of importance.
!    IOTYP=0: Selected species (11) only
!    IOTYP=1: All gas-phase species
!    IOTYP=2: All species
!------------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!
      PARAMETER(NSP=214,NRM=2345,NSEL=7)
!
      CHARACTER*25 R1(NRM),R2(NRM),P1(NRM),P2(NRM),P3(NRM),P4(NRM), &
                 SPECI(NSP),SPSEL(NSEL),SPEC
      REAL*8 G(NRM),R(NRM),X(NSP),Y0(NSP),FLAG(10)
      INTEGER INDXJ(NRM),LAB(10) 
!
      COMMON/BLK1/SPECI
      COMMON/BLK3/R,D
      COMMON/ORDER1/R1,R2,P1,P2,P3,P4
      COMMON/ORDER2/INDXJ,IOTYP
!
      ! DATA  SPSEL/'N2H+','N2','NH3','ELECTR','CO','C','C+','HCO+', &
      !      'CS','O','H2O','#NH3','#N2','#CO','#H2O','#O2'/
!      DATA (SPSEL(I),I=1,33)/'H2+','H3+','CH','CH+','CH2','CH2+','CH3+',
!     *     'CN','N2','N2H+','N+','NH','NH+','NH2','NH2+','NH3+','NH4+',
!     *     'NH3','H2NC+','NO','HCO+','H3O+','OH','O2','CO','C2H+','H+',
!     *     'HE+','C+','S+','C2S','HC2S+','ELECTR'/
      DATA SPSEL/'HCO','HCO+','H2CO','CH2','H2','O','E-'/
      DATA LOR/2/
!
      WRITE(LOR,4) RAD,TIME
 4    FORMAT(1X,'Analysis at ',1PE8.2,' pc. and ',1PE8.2,' years',/)
!
      DO 112 K1=1,NRM
         G(K1) = R(K1)
 112  CONTINUE
!       DO 113 K1=1,NCONS
!          X(K1) = XS(K1)*D
!  113  CONTINUE
      DO 114 K1=1,NTD
         X(K1) = Y0(K1)*D
 114  CONTINUE
      NTOT=NTD
!
      DO 1 J=1,NRM
        IF(INDXJ(J).EQ.9999) GO TO 2
        DO 14 LI=1,NTOT
            IF(R1(J).EQ.SPECI(LI)) G(J)=G(J)*X(LI)
            IF(R2(J).EQ.SPECI(LI)) G(J)=G(J)*X(LI)
 14     CONTINUE
        IF(R1(J).EQ.'FREEZE') G(J)=G(J)*D
        IF(R2(J).EQ.'FREEZE') G(J)=G(J)*D
 1    CONTINUE
 2    L=J-1
!
!------------------------------------------------------------------------------
!
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
!
!----Select species
!
      IF(IOTYP.EQ.0) LMAX=NSEL
      IF(IOTYP.EQ.1) LMAX=NCONS+NG
      IF(IOTYP.EQ.2) LMAX=NTOT
!
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
          IF((R1(M).EQ.SPEC).OR.(R2(M).EQ.SPEC).OR.(P1(M).EQ.SPEC).OR. &
           (P2(M).EQ.SPEC).OR.(P3(M).EQ.SPEC).OR.(P4(M).EQ.SPEC)) THEN
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
!
      WRITE(LOR,20)
 20   FORMAT(1X,78('-'),/)
 15   FORMAT(1X,I4,1X,5(1A8,1X),A8,2X,I4,'%')
 9    FORMAT(//,1X,5X,12('*'))
 10   FORMAT(1X,5X,12('*'),/)
 11   FORMAT(1X,5X,'*',1X,1A8,1X,'*',2X,'Frate=',1PE9.3,2X,'Drate=', &
            1PE9.3,2X,'Abundance=',1PE9.3)
      RETURN
      END
!
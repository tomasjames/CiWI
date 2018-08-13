 PROGRAM CIWI
! 
!------------------------------------------------------------------------------
!
! CIWI: Chemistry In the Wardle Instability
!
!------------------------------------------------------------------------------
      USE interp
      USE params

      IMPLICIT NONE ! Ensures that implicit types are ignored
      CHARACTER*11 :: DGFILE,RELUFILE,RHOEFILE,RHOGFILE,RHOIFILE, &
            RHONFILE,TEEFILE,TEIFILE,TENFILE,TRACKFILE,ZGFILE
      INTEGER :: i,j,k,numfiles,numelem
      REAL :: a
      INTEGER, PARAMETER :: filelength=2331,interpLength=3000
      DOUBLE PRECISION :: T1,T2
      DOUBLE PRECISION :: T(filelength),DG(filelength),RELU(filelength), &
            RHOE(filelength),RHOG(filelength),RHOI(filelength),RHON(filelength), &
            TEE(filelength),TEI(filelength),TEN(filelength), X(filelength), &
            Y(filelength),ZG(filelength)

      DOUBLE PRECISION, PARAMETER :: L_BOUND=0, U_BOUND=2.33e+10
      DOUBLE PRECISION, DIMENSION(filelength) :: T_INTERP,DG_INTERP,RELU_INTERP, &
            RHOE_INTERP,RHOG_INTERP,RHOI_INTERP,RHON_INTERP,TEE_INTERP,TEI_INTERP, &
            TEN_INTERP,ZG_INTERP

      DOUBLE PRECISION :: DG_POINT

!---------------------------- INPUT PARAMETERS --------------------------------
      ! I/O files
      DGFILE='dg.xq'
      RELUFILE='relu.xq'
      RHOEFILE='rhoe.xq'
      RHOGFILE='rhog.xq'
      RHOIFILE='rhoi.xq'
      RHONFILE='rhon.xq'
      TEEFILE='tee.xq'
      TEIFILE='tei.xq'
      TENFILE='ten.xq'
      TRACKFILE='track.xq'
      ZGFILE='zg.xq'
!-------------------------- END OF INPUT PARAMETERS ---------------------------

!-------------------------------- FILE ADMIN ----------------------------------
      PRINT *, 'Opening datafiles...'
      OPEN(1,FILE=DGFILE,STATUS='OLD',ACTION='READ')
      OPEN(2,FILE=RELUFILE,STATUS='OLD',ACTION='READ')
      OPEN(3,FILE=RHOEFILE,STATUS='OLD',ACTION='READ')
      OPEN(4,FILE=RHOGFILE,STATUS='OLD',ACTION='READ')
      OPEN(5,FILE=RHOIFILE,STATUS='OLD',ACTION='READ')
      OPEN(7,FILE=RHONFILE,STATUS='OLD',ACTION='READ')
      OPEN(8,FILE=TEEFILE,STATUS='OLD',ACTION='READ')
      OPEN(9,FILE=TEIFILE,STATUS='OLD',ACTION='READ')
      OPEN(10,FILE=TENFILE,STATUS='OLD',ACTION='READ')
      OPEN(11,FILE=TRACKFILE,STATUS='OLD',ACTION='READ')
      OPEN(12,FILE=ZGFILE,STATUS='OLD',ACTION='READ')

      PRINT *, 'Opening files to write interpolated data§'
      OPEN(13,FILE='dg_interp.dat',ACTION='WRITE')
      OPEN(14,FILE='test.dat',ACTION='WRITE')

      PRINT *, 'Reading in data...'

      ! This line loops through all filelength lines in all files
      ! and reads each line, assigning the information in
      ! each column (there are 2 columns) to the ith element in 
      ! the declared arrays
      ! Any int (such as a) used after this point is a dummy
      ! variable - likely because each file contains time, so
      ! multiple copies of the same data are not required.
      DO i=1,filelength
            READ(1,FMT="(2(E12.6,1X))") T(i),DG(i)
            READ(2,FMT="(2(E12.6,1X))") a,RELU(i)
            READ(3,FMT="(2(E12.6,1X))") a,RHOE(i)
            READ(4,FMT="(2(E12.6,1X))") a,RHOG(i)
            READ(5,FMT="(2(E12.6,1X))") a,RHOI(i)
            READ(7,FMT="(2(E12.6,1X))") a,RHON(i)
            READ(8,FMT="(2(E12.6,1X))") a,TEE(i)
            READ(9,FMT="(2(E12.6,1X))") a,TEI(i)
            READ(10,FMT="(2(E12.6,1X))") a,TEN(i)
            READ(11,FMT="(2(E12.6,1X))") X(i),Y(i)
            READ(12,FMT="(2(E12.6,1X))") a,ZG(i)
      END DO

      PRINT *, 'Data read in complete.'
      PRINT *, 'Moving to interpolation...'

!--------------------------------- INTERPOLATE ----------------------------------
      ! Need to create an array of times to interpolate at
      DO i=1,filelength
            IF (i .eq. 1) THEN
                  T_INTERP(i) = L_BOUND
            ELSE IF (i .eq. filelength) THEN
                  T_INTERP(i) = U_BOUND
            ELSE 
                  T_INTERP(i) = T_INTERP(i-1) + (U_BOUND - L_BOUND)/(filelength)
            END IF
      END DO

      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=DG,Y_INTERP=DG_INTERP, &
            DATAFILE='dg_interp.dat',interpLength=interpLength, &
            filelength=filelength)
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=RELU,Y_INTERP=RELU_INTERP, &
            DATAFILE='relu_interp.dat',interpLength=interpLength, &
            filelength=filelength)
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=RHOE,Y_INTERP=RHOE_INTERP, &
            DATAFILE='rhoe_interp.dat',interpLength=interpLength, &
            filelength=filelength)
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=RHOG,Y_INTERP=RHOG_INTERP, &
            DATAFILE='rhog_interp.dat',interpLength=interpLength, &
            filelength=filelength)
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=RHOI,Y_INTERP=RHOI_INTERP, &
            DATAFILE='rhoi_interp.dat',interpLength=interpLength, &
            filelength=filelength)
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=RHON,Y_INTERP=RHON_INTERP, &
            DATAFILE='rhon_interp.dat',interpLength=interpLength, &
            filelength=filelength)
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=TEE,Y_INTERP=TEE_INTERP, &
            DATAFILE='tee_interp.dat',interpLength=interpLength, &
            filelength=filelength)
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=TEI,Y_INTERP=TEI_INTERP, &
            DATAFILE='tei_interp.dat',interpLength=interpLength, &
            filelength=filelength)
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=TEN,Y_INTERP=TEN_INTERP, &
            DATAFILE='ten_interp.dat',interpLength=interpLength, &
            filelength=filelength)
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=ZG,Y_INTERP=ZG_INTERP, &
            DATAFILE='zg_interp.dat',interpLength=interpLength, &
            filelength=filelength)

!----------------------------------- CALL ------------------------------------

      CALL phys(XDATA=T_INTERP,YDATA=DG_INTERP,XDATAPOINT=T_INTERP(3), &
            filelength=2331,YDATAPOINT=DG_POINT)
      CALL phys(XDATA=T_INTERP,YDATA=TEI_INTERP,XDATAPOINT=T_INTERP(3), &
            filelength=2331,YDATAPOINT=DG_POINT)
 END PROGRAM
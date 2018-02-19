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
      CHARACTER(5) :: DGFILE,ZGFILE
      CHARACTER(6) :: TEEFILE,TEIFILE,TENFILE
      CHARACTER(7) :: RELUFILE,RHOEFILE,RHOGFILE,RHOIFILE, &
            RHONFILE
      CHARACTER(8) :: TRACKFILE
      INTEGER :: i,j,k,numfiles,numelem
      REAL :: a
      INTEGER, PARAMETER :: filelength=2331,interpLength=3000
      DOUBLE PRECISION :: T1,T2
      DOUBLE PRECISION, DIMENSION(filelength) :: T,DG,RELU,RHOE,RHOG, &
            RHOI,RHON,TEE,TEI,TEN,X,Y,ZG

      DOUBLE PRECISION, PARAMETER :: L_BOUND=0, U_BOUND=2.330015E10
      DOUBLE PRECISION, DIMENSION(interpLength) :: T_INTERP,DG_INTERP,RELU_INTERP, &
            RHOE_INTERP,RHOG_INTERP,RHOI_INTERP,RHON_INTERP,TEE_INTERP,TEI_INTERP, &
            TEN_INTERP,ZG_INTERP

      DOUBLE PRECISION :: DG_RETURN,RELU_RETURN,RHOE_RETURN,RHOG_RETURN,RHOI_RETURN, &
            RHON_RETURN,TEE_RETURN,TEI_RETURN,TEN_RETURN,ZG_RETURN

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
      DO i=1,interpLength
            IF (i == 1) THEN
                  T_INTERP(i) = L_BOUND
            ELSE IF (i == interpLength) THEN
                  T_INTERP(i) = U_BOUND
            ELSE 
                  T_INTERP(i) = T_INTERP(i-1) + (U_BOUND - L_BOUND)/(interpLength)
            END IF
      END DO

      ! Call the function for each parameter to be interpolated
      ! NOTE: function will return an array of Y_INTERP values
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=DG, &
            DATAFILE='dg_interp.dat',interpLength=interpLength, &
            filelength=filelength,Y_INTERP=DG_INTERP)
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=RELU, &
            DATAFILE='relu_interp.dat',interpLength=interpLength, &
            filelength=filelength,Y_INTERP=RELU_INTERP)
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=RHOE, &
            DATAFILE='rhoe_interp.dat',interpLength=interpLength, &
            filelength=filelength,Y_INTERP=RHOE_INTERP)
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=RHOG, &
            DATAFILE='rhog_interp.dat',interpLength=interpLength, &
            filelength=filelength,Y_INTERP=RHOG_INTERP)
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=RHOI, &
            DATAFILE='rhoi_interp.dat',interpLength=interpLength, &
            filelength=filelength,Y_INTERP=RHOI_INTERP)
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=RHON, &
            DATAFILE='rhon_interp.dat',interpLength=interpLength, &
            filelength=filelength,Y_INTERP=RHON_INTERP)
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=TEE, &
            DATAFILE='tee_interp.dat',interpLength=interpLength, &
            filelength=filelength,Y_INTERP=TEE_INTERP)
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=TEI, &
            DATAFILE='tei_interp.dat',interpLength=interpLength, &
            filelength=filelength,Y_INTERP=TEI_INTERP)
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=TEN, &
            DATAFILE='ten_interp.dat',interpLength=interpLength, &
            filelength=filelength,Y_INTERP=TEN_INTERP)
      CALL linear_interp(T_INTERP=T_INTERP,T=T,Y=ZG, &
            DATAFILE='zg_interp.dat',interpLength=interpLength, &
            filelength=filelength,Y_INTERP=ZG_INTERP)

!----------------------------------- COUPLE ------------------------------------
      PRINT *, 'DG_INTERP: ', DG_INTERP(1)
      CALL PHYS(T0=T_INTERP(1),T=T_INTERP,DG_INTERP=DG_INTERP,RELU_INTERP=RELU_INTERP, &
            RHOE_INTERP=RHOE_INTERP,RHOG_INTERP=RHOG_INTERP,RHOI_INTERP=RHOI_INTERP, & 
            RHON_INTERP=RHON_INTERP,TEE_INTERP=TEE_INTERP,TEI_INTERP=TEI_INTERP, &
            TEN_INTERP=TEN_INTERP,ZG_INTERP=ZG_INTERP,DG_RETURN=DG_RETURN, &
            RELU_RETURN=RELU_RETURN,RHOE_RETURN=RHOE_RETURN,RHOG_RETURN=RHOG_RETURN, &
            RHOI_RETURN=RHOI_RETURN,RHON_RETURN=RHON_RETURN,TEE_RETURN=TEE_RETURN, &
            TEI_RETURN=TEI_RETURN,TEN_RETURN=TEN_RETURN,ZG_RETURN=ZG_RETURN, &
            interpLength=interpLength)
      PRINT *, 'DG_RETURN: ', DG_RETURN
      PRINT *, 'RELU_RETURN: ', RELU_RETURN

 END PROGRAM
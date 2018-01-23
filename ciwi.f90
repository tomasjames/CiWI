 PROGRAM CIWI
! 
!------------------------------------------------------------------------------
!
! CIWI: Chemistry In the Wardle Instability
!
!------------------------------------------------------------------------------
      USE interpolate

      IMPLICIT NONE ! Ensures that implicit types are ignored
      CHARACTER*11 :: DGFILE,RELUFILE,RHOEFILE,RHOGFILE,RHOIFILE, &
            RHONFILE,TEEFILE,TEIFILE,TENFILE,TRACKFILE,ZGFILE
      INTEGER :: i,numfiles
      REAL :: a
      INTEGER, PARAMETER :: filelength=2331
      DOUBLE PRECISION :: T(filelength),DG(filelength),RELU(filelength), &
            RHOE(filelength),RHOG(filelength),RHOI(filelength),RHON(filelength), &
            TEE(filelength),TEI(filelength),TEN(filelength), X(filelength), &
            Y(filelength),ZG(filelength)

      INTEGER ( KIND=4 ), PARAMETER :: DATA_NUM = 11
      INTEGER ( KIND=4 ) INTERP
      INTEGER( KIND=4 ), PARAMETER :: INTERP_NUM = 200
      REAL( KIND=8 ), ALLOCATABLE, DIMENSION ( : ) :: X_INTERP
      REAL( KIND=8 ), ALLOCATABLE, DIMENSION ( :, : ) :: Y_INTERP
      REAL( KIND=8 ), ALLOCATABLE, DIMENSION ( : ) :: Y_VALUE
      REAL( KIND=8 ), PARAMETER :: X_MAX=0
      REAL( KIND=8 ), PARAMETER :: X_MIN=3e10

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

!------------------------------- DATA READ IN ---------------------------------
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
!--------------------------------- INTERPOLATE ----------------------------------
      PRINT *, 'Moving to interpolation...'
      ALLOCATE (X_INTERP(INTERP_NUM))
      ALLOCATE (Y_INTERP(1,INTERP_NUM))

      PRINT *, 'Calling interp_linear...'
      CALL interp_linear (1, DATA_NUM, T, DG, INTERP_NUM, &
            X_INTERP, Y_INTERP)

      ALLOCATE(Y_VALUE(1:INTERP_NUM))

      DO INTERP=1, INTERP_NUM
            ! WRITE ( *, '(2x,f10.4,2x,g14.6,2x,g14.6,2x,g10.2)' ) &
            !       X_INTERP(INTERP), Y_INTERP(1,INTERP), Y_VALUE(INTERP), &
            !       Y_INTERP(1,INTERP) - Y_VALUE(INTERP)
            PRINT *, 'X_INTERP(INTERP): ', X_INTERP(INTERP)
            PRINT *, 'Y_INTERP(1,INTERP): ', Y_INTERP(1,INTERP)
            PRINT *, 'Y_VALUE(INTERP): ', Y_VALUE(INTERP)
            PRINT *, 'Y_INTERP(1,INTERP) - Y_VALUE(INTERP)', Y_INTERP(1,INTERP) - Y_VALUE(INTERP)
      END DO
 END
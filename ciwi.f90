 PROGRAM CIWI
! 
!------------------------------------------------------------------------------
!
! CIWI: Chemistry In the Wardle Instability
!
!------------------------------------------------------------------------------
      IMPLICIT NONE ! Ensures that implicit types are ignored
      CHARACTER*11 :: DGFILE,RELUFILE,RHOEFILE,RHOGFILE,RHOIFILE, &
            RHONFILE,TEEFILE,TEIFILE,TENFILE,TRACKFILE,ZGFILE
      INTEGER :: i,a
      INTEGER, PARAMETER :: filelength=2331
      DOUBLE PRECISION :: T(filelength),DG(filelength),RELU(filelength), &
            RHOE(filelength),RHOG(filelength),RHOI(filelength),RHON(filelength), &
            TEE(filelength),TEI(filelength),TEN(filelength), X(filelength), &
            Y(filelength),ZG(filelength)
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
      OPEN(1,FILE=DGFILE,STATUS='OLD',ACTION='READ')
      OPEN(2,FILE=RELUFILE,STATUS='OLD',ACTION='READ')
      OPEN(3,FILE=RHOEFILE,STATUS='OLD',ACTION='READ')
      OPEN(4,FILE=RHOGFILE,STATUS='OLD',ACTION='READ')
      OPEN(5,FILE=RHOIFILE,STATUS='OLD',ACTION='READ')
      OPEN(6,FILE=RHONFILE,STATUS='OLD',ACTION='READ')
      OPEN(7,FILE=TEEFILE,STATUS='OLD',ACTION='READ')
      OPEN(8,FILE=TEIFILE,STATUS='OLD',ACTION='READ')
      OPEN(9,FILE=TENFILE,STATUS='OLD',ACTION='READ')
      OPEN(10,FILE=TRACKFILE,STATUS='OLD',ACTION='READ')
      OPEN(11,FILE=ZGFILE,STATUS='OLD',ACTION='READ')

      ! This line loops through all 2331 lines in all files
      ! and reads each line, assigning the information in
      ! each column (there are 2 columns) to the ith element in 
      ! the declared arrays
      ! Any int (such as a) used after this point is a dummy
      ! variable - likely because each file contains time, so
      ! multiple copies of the same data are not required.
      DO i=1,2331
            READ(1,FMT="(2(E12.6,1X))") T(i),DG(i)
            READ(2,FMT="(2(E12.6,1X))") a,RELU(i)
            READ(3,FMT="(2(E12.6,1X))") a,RHOE(i)
            READ(4,FMT="(2(E12.6,1X))") a,RHOG(i)
            READ(5,FMT="(2(E12.6,1X))") a,RHOI(i)
            READ(6,FMT="(2(E12.6,1X))") a,RHON(i)
            READ(7,FMT="(2(E12.6,1X))") a,TEE(i)
            READ(8,FMT="(2(E12.6,1X))") a,TEI(i)
            READ(9,FMT="(2(E12.6,1X))") a,TEN(i)
            READ(10,FMT="(2(E12.6,1X))") X(i),Y(i)
            READ(11,FMT="(2(E12.6,1X))") a,ZG(i)
      END DO
      
      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)
      CLOSE(5)
      CLOSE(6)
      CLOSE(7)
      CLOSE(8)
      CLOSE(9)
      CLOSE(10)
      CLOSE(11)
!----------------------------------------------------------------------------
 END
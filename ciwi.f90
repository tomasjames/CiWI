 PROGRAM CIWI
! 
!------------------------------------------------------------------------------
!
! CIWI: Chemistry In the Wardle Instability
!
!------------------------------------------------------------------------------
      CHARACTER*11 :: DGFILE,RELUFILE,RHOEFILE,RHOGFILE,RHOIFILE
      INTEGER :: i
      DOUBLE PRECISION :: T(2331),DG(2331)
!---------------------------- INPUT PARAMETERS --------------------------------
!----I/O files
      DGFILE='dg.xq'
      RELUFILE='relu.xq'
      RHOEFILE='rhoe.xq'
      RHOGFILE='rhog.xq'
      RHOIFILE='rhoi.xq'
!     RHONFILE='rhon.xq'
!     TEEFILE='tee.xq'
!     TEIFILE='tei.xq'
!     TENFILE='ten.xq'
!     TRACKFILE='track.xq'
!     ZGFILE='zq.xq'
!-------------------------- END OF INPUT PARAMETERS ---------------------------

!------------------------------- DATA READ IN ---------------------------------
      OPEN(1,FILE=DGFILE,STATUS='OLD',ACTION='READ')
      DO i=1,2331
            READ(1,FMT="(2(E12.6,1X))") T(i),DG(i)
      END DO
      CLOSE(1)
!----------------------------------------------------------------------------
 END
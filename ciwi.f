      PROGRAM CIWI
C 
C------------------------------------------------------------------------------
C
C CIWI: Chemistry In the Wardle Instability
C
C------------------------------------------------------------------------------
C 
      CHARACTER*11 DGFILE,RELUFILE,RHOEFILE,RHOGFILE,RHOIFILE,RHONFILE,TEEFILE,TEIFILE,TENFILE,TRACKFILE,ZGFILE
C 
C---------------------------- INPUT PARAMETERS --------------------------------
C----I/O files
      DGFILE='dg.xq'
      RELUFILE='relu.xq'
      RHOEFILE='rhoe.xq'
      RHOGFILE='rhog.xq'
      RHOIFILE='rhoi.xq'
      RHONFILE='rhon.xq'
C       TEEFILE='tee.xq'
C       TEIFILE='tei.xq'
C       TENFILE='ten.xq'
C       TRACKFILE='track.xq'
C       ZGFILE='zq.xq'
C-------------------------- END OF INPUT PARAMETERS ---------------------------
C 
C----------------------------- PRINT STATEMENTS -------------------------------
C 
      PRINT *, DGFILE
C 
      END
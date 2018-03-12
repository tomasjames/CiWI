MODULE interp
    IMPLICIT NONE

    CONTAINS
        SUBROUTINE linear_interp (T_INTERP,T,Y,DATAFILE,interpLength,filelength)

            INTEGER :: j,k
            DOUBLE PRECISION :: T1,T2
            DOUBLE PRECISION :: Y_INTERP(interpLength)

            DOUBLE PRECISION, INTENT(IN) :: T_INTERP(interpLength)
            DOUBLE PRECISION, INTENT(IN) :: T(filelength),Y(filelength)
            CHARACTER(*), INTENT(IN) :: DATAFILE
            INTEGER, INTENT(IN) :: interpLength,filelength

            OPEN(13,FILE=DATAFILE,ACTION='WRITE')

            DO k=1,interpLength
                DO j=1,filelength
                      ! Assess to determine whether point sits in specified range 
                      IF (T_INTERP(k) > T(j) .AND. T_INTERP(k) < T(j+1)) THEN
                            T1 = T(j)
                            T2 = T(j+1)

                            Y_INTERP(k) = Y(j) + (T_INTERP(k)-T1)*((Y(j+1)-Y(j))/(T2-T1))
                            WRITE(13,*) T_INTERP(k),Y_INTERP(k)
                      ELSE IF (T_INTERP(k) < T(1)) THEN
                            T1 = T_INTERP(k)
                            T2 = T(1)

                            Y_INTERP(k) = Y(j) + (T_INTERP(k)-T1)*((Y(j+1)-Y(j))/(T2-T1))
                            WRITE(13,*) T_INTERP(k),Y_INTERP(k)
                      ELSE IF (T_INTERP(k) > T(filelength)) THEN
                            T1 = T(filelength)
                            T2 = T_INTERP(k)

                            Y_INTERP(k) = Y(j) + (T_INTERP(k)-T1)*((Y(j+1)-Y(j))/(T2-T1))
                            WRITE(13,*) T_INTERP(k),Y_INTERP(k)
                      END IF
                END DO
          END DO

        END SUBROUTINE

END MODULE interp
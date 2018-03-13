MODULE params

    IMPLICIT NONE

    CONTAINS
        SUBROUTINE phys (XDATA,YDATA,XDATAPOINT,filelength,YDATAPOINT)

            DOUBLE PRECISION, INTENT(IN) :: XDATA(filelength),YDATA(filelength)
            DOUBLE PRECISION, INTENT(IN) :: XDATAPOINT
            DOUBLE PRECISION, INTENT(OUT) :: YDATAPOINT

            INTEGER, INTENT(IN) :: filelength

            INTEGER :: i

            DO i=1,filelength
                IF (XDATA(i) == XDATAPOINT) THEN
                    YDATAPOINT = YDATA(i)
                END IF
            END DO

        END SUBROUTINE phys

END MODULE params

MODULE params
    IMPLICIT NONE

    CONTAINS
        SUBROUTINE PHYS(T0,T,DG_INTERP,RELU_INTERP,RHOE_INTERP,RHOG_INTERP,RHOI_INTERP, & 
            RHON_INTERP,TEE_INTERP,TEI_INTERP,TEN_INTERP,ZG_INTERP, &
            DG_RETURN,RELU_RETURN,RHOE_RETURN,RHOG_RETURN,RHOI_RETURN, &
            RHON_RETURN,TEE_RETURN,TEI_RETURN,TEN_RETURN,ZG_RETURN,interpLength)

            DOUBLE PRECISION :: T0
            INTEGER :: interpLength,j

            DOUBLE PRECISION, DIMENSION(interpLength) :: T,DG_INTERP,RELU_INTERP,RHOE_INTERP, &
                RHOG_INTERP,RHOI_INTERP,RHON_INTERP,TEE_INTERP,TEI_INTERP, & 
                TEN_INTERP,ZG_INTERP
            DOUBLE PRECISION :: DG_RETURN,RELU_RETURN,RHOE_RETURN,RHOG_RETURN, &
                RHOI_RETURN,RHON_RETURN,TEE_RETURN,TEI_RETURN,TEN_RETURN,ZG_RETURN

        !----------------------------- FIND RELEVANT INFO -------------------------------

        DO j=1,interpLength
            IF (T0==T(j)) THEN
                DG_RETURN = DG_INTERP(j)
                RELU_RETURN = RELU_INTERP(j)
                RHOE_RETURN = RHOE_INTERP(j)
                RHOG_RETURN = RHOG_INTERP(j)
                RHOI_RETURN = RHOI_INTERP(j)
                RHON_RETURN = RHON_INTERP(j)
                TEE_RETURN = TEE_INTERP(j)
                TEI_RETURN = TEI_INTERP(j)
                TEN_RETURN = TEN_INTERP(j)
                ZG_RETURN = ZG_INTERP(j)
                RETURN
            END IF
        END DO

        END SUBROUTINE
END MODULE params
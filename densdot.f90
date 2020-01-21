MODULE collapse
    implicit none
    double precision :: initialDens,finalDens,bc
    CONTAINS
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Returns the time derivative of the density.                                     !
        ! Analytical function taken from Rawlings et al. 1992                             !
        ! Called from chemistry.f90, density integrated with abundances so this gives ydot!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        pure FUNCTION densdot(density)
            DOUBLE PRECISION, INTENT(IN) :: density
            DOUBLE PRECISION :: densdot
            !Rawlings et al. 1992 freefall collapse. With factor bc for B-field etc
            IF (density .lt. finalDens) THEN
                densdot=bc*(density**4./initialDens)**0.33*&
                &(8.4d-30*initialDens*((density/initialDens)**0.33-1.))**0.5
            ELSE
                densdot=0.0
            ENDIF
        END FUNCTION densdot
END MODULE
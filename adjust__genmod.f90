        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar  5 16:45:36 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ADJUST__genmod
          INTERFACE 
            SUBROUTINE ADJUST(N,T,Y,YDOT)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: T
              REAL(KIND=8) :: Y(N)
              REAL(KIND=8) :: YDOT(N)
            END SUBROUTINE ADJUST
          END INTERFACE 
        END MODULE ADJUST__genmod

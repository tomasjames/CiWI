        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar  5 16:45:36 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DIFFUN__genmod
          INTERFACE 
            SUBROUTINE DIFFUN(N,T,Y,YDOT)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: T
              REAL(KIND=8) :: Y(N)
              REAL(KIND=8) :: YDOT(N)
            END SUBROUTINE DIFFUN
          END INTERFACE 
        END MODULE DIFFUN__genmod
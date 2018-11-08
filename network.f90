    PROGRAM NETWORK
    ! ------------------------------------------------------------------------------
    !
    ! Code to read in network and format it correctly for wardle.f read 
    !
    ! ------------------------------------------------------------------------------
    
    IMPLICIT NONE

    INTEGER reaction,species,reactionFileLength,speciesFileLength,body, &
        & index
    INTEGER NG,NS,EFI
    DATA    NG,NS/0,0/
    PARAMETER(reactionFileLength=2318,speciesFileLength=214)
    CHARACTER*40 reactionFileInput,reactionFileOutput,speciesFileInput, &
                & speciesFileOutput
    CHARACTER*8, DIMENSION(3,reactionFileLength) :: REACT
    CHARACTER*8, DIMENSION(4,reactionFileLength) :: PROD
    CHARACTER*8, DIMENSION(speciesFileLength) :: SPEC
    REAL, DIMENSION(reactionFileLength) :: ALPHA,BETA,GAMMA
    REAL :: DUMMY

    ! Define the end file int
    EFI = 9999

    ! Define the file to be opened
    reactionFileInput="rates/outputFiles/reactions.csv"
    speciesFileInput="rates/outputFiles/species.csv"
    reactionFileOutput="coll_rates.d"
    speciesFileoutput="coll_specs.d"

    ! Open the reactions file
    OPEN(UNIT=1,FILE=reactionFileInput,ACTION="READ")
    OPEN(UNIT=3,FILE=speciesFileInput,ACTION="READ")

    ! Read in the data
    DO reaction=1,reactionFileLength
        READ(1,*) REACT(1,reaction), REACT(2,reaction), REACT(3,reaction), &
                & PROD(1,reaction), PROD(2,reaction), PROD(3,reaction), &
                & PROD(4,reaction), ALPHA(reaction), & 
                & BETA(reaction), GAMMA(reaction)
        DO body=2,4
            IF (PROD(body,reaction) .EQ. "NAN") THEN
                PROD(body,reaction) = "        "
            END IF
        END DO
    END DO

    ! Read in the species
    DO species=1,speciesFileLength
        READ(3,*) SPEC(species), DUMMY, DUMMY
        IF (SPEC(species)(1:1) .NE. "G") THEN
            NG=NG+1
        ELSE IF (SPEC(species)(1:1) .EQ. "#") THEN
            NS=NS+1
        END IF
    END DO

    CLOSE(1)
    CLOSE(3)

    ! Open the rate file to write out to
    OPEN(UNIT=2,FILE=reactionFileOutput,ACTION="WRITE")
    OPEN(UNIT=4,FILE=speciesFileOutput,ACTION="WRITE")

    DO reaction=1,reactionFileLength
        WRITE(2,FMT="(I4,4(1X,A8),2(1X,A4),1X,1PE10.2,1X,0PF6.2,1X,F8.1)") reaction, &
                & REACT(1,reaction), REACT(2,reaction), &
                & PROD(1,reaction), PROD(2,reaction), &
                & PROD(3,reaction), PROD(4,reaction), ALPHA(reaction), & 
                & BETA(reaction), GAMMA(reaction)
    END DO
    WRITE(2,FMT="(I4)") EFI

    DO species=1,speciesFileLength
        ! WRITE(4,FMT="(A2,I4,A2,A9,T1)",ADVANCE="NO") "Y(", species, ")=", ADJUSTL(SPEC(species))
        ! IF (MOD(species,7) .EQ. 0) THEN
        !     WRITE(4,*) CHAR(13)
        ! END IF
        WRITE(4,FMT="(A2,I4,A2,A9)") "Y(", species, ")=", ADJUSTL(SPEC(species))
    END DO

    CLOSE(2)
    CLOSE(4)

    END PROGRAM NETWORK
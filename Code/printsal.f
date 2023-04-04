      SUBROUTINE PRINTSAL(IUSAL)
C     PRINTS SALINITY FROM SELECTED LOCATIONS

      INCLUDE 'comdeck'  

      PARAMETER(NPALM=12)
      INTEGER IUSAL
      INTEGER ICOORD(NPALM),JCOORD(NPALM)

      DATA ICOORD/65,65,65,64,63,62,62,61,60,59,58,57/
      DATA JCOORD/75,74,73,73,73,73,72,72,72,72,72,72/

C---------------------------------------------------------------
C     THE SELECTED SALINITY LOCATIONS ARE:
C     1) THE SURFACE LAYER 
C     2) THE BOTTOM LAYER
C     3) THE PALM RIVER X-Z PROFILES
C     4) LOCATIONS IN MCKAY BAY AND EAST BAY, X-Z AND ISOLATED POINTS
C---------------------------------------------------------------

      WRITE(IUSAL,*) 'TIME (days) FROM MIDNIGHT JAN 1, 1990'
      WRITE(IUSAL,*) TIME

C     ITEMS #1-2
      WRITE(IUSAL,*) 'SURFACE LAYER'
      WRITE(IUSAL,*) ((S(I,J,1), I=1,IM), J=1,JM)
      WRITE(IUSAL,*) 'BOTTOM LAYER'
      WRITE(IUSAL,*) ((S(I,J,KB), I=1,IM), J=1,JM)

C     ITEM #3 
      WRITE(IUSAL,*) '# OF POINTS'
      WRITE(IUSAL,*) NPALM
      DO 50 N=1,NPALM
      DO 50 K=1,KB
        I=ICOORD(N)
        J=JCOORD(N)
        WRITE(IUSAL,100) I,J,K,S(I,J,K)
 50   CONTINUE

 100  format(3(i3,2x),f10.6)


      END






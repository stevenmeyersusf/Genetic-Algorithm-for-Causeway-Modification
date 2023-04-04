      SUBROUTINE SALPOINTS(IUSAL4,OPEN4)
C     THIS ROUTINE OUTPUTS SALINITY VALUES FOR SELECT POINTS

      include 'comdeck'

      PARAMETER(NPOINTS=9)
      INTEGER IX(NPOINTS),IY(NPOINTS),IUSAL4,OPEN4

      DATA IX /65,62,59,57,54,56,54,57,55/
      DATA IY /74,73,72,74,72,71,69,69,66/

c     write headers at top of file
      IF (OPEN4.EQ.0) THEN 
         write(6,*) 'starting sal_points' 
         WRITE(IUSAL4,*) 'NUMBER OF VERTICAL LEVELS:'
         WRITE(IUSAL4,*) KB
         WRITE(IUSAL4,*) 'NUMBER OF POINTS:'
         WRITE(IUSAL4,*) NPOINTS
         DO M=1,NPOINTS
            WRITE(IUSAL4,*) IX(M),IY(M)
         ENDDO
         OPEN4=1
      ENDIF

C     EVERY 20 TIME STEPS (20 MINUTES)
      IF (MOD(INT,20).EQ.0) THEN   
         WRITE(IUSAL4,*) TIME
         DO 100 K=1,KB
            DO M=1,NPOINTS
               WRITE(IUSAL4,*) S(IX(M),IY(M),K)
            ENDDO
 100     CONTINUE
      ENDIF

      END

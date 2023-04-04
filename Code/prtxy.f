      SUBROUTINE PRTXY(FD,FMASK,IM,JM,KT,IARRAY,SCALE,IUNIT)
C     VERSION(11/18/90)
C
      DIMENSION FD(IM,JM,KT),IARRAY(IM,JM),FMASK(IM,JM)
      DIMENSION KP(2)
      DATA KP/1,10/
C
C-----------------------------------------------------------------------
C        THIS SUBROUTINE WRITES HORIZONTAL LAYERS OF A 3-D FIELD
C-----------------------------------------------------------------------
C
      DO 10 KM=1,2
      K=KP(KM)
      WRITE(IUNIT,20) K
 20   FORMAT(3X,/7H LAYER ,I2/)
      CALL PRINT(FD(1,1,K),FMASK,IM,JM,17,IARRAY,SCALE,IUNIT)
 10   CONTINUE
C
      RETURN
      END
      

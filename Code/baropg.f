      SUBROUTINE BAROPG(DRHOX,DRHOY,TRNU,TRNV)
C     VERSION(09/30/91)
      INCLUDE 'comdeck'
      DIMENSION DRHOX(IM,JM,KB),DRHOY(IM,JM,KB),TRNU(IM,JM),TRNV(IM,JM)
C
C     DO 50 K=1,KB
C     DO 50 J=1,JM
C     DO 50 I=1,IM
C50   RHO(I,J,K)=RHO(I,J,K)-RMEAN(I,J,K)
C
C-------- X COMPONENT OF BAROCLINIC PRESSURE GRADIENT ------------------
      DO 300 J=2,JMM1
      DO 300 I=2,IMM1
 300  DRHOX(I,J,1)=.25*GRAV*DZ(1)*(DT(I,J)+DT(I-1,J))
     .    *(RHO(I,J,1)-RHO(I-1,J,1))
      DO 310 K=2,KBM1
      DO 310 J=2,JMM1
      DO 310 I=2,IMM1
 310  DRHOX(I,J,K)=DRHOX(I,J,K-1)
     .      +GRAV*.125*(DZ(K)+DZ(K-1))*(DT(I,J)+DT(I-1,J))
     .      *(RHO(I,J,K)-RHO(I-1,J,K)+RHO(I,J,K-1)-RHO(I-1,J,K-1))
     .      +GRAV*.5*Z(K)*(DT(I,J)-DT(I-1,J))
     .      *(RHO(I,J,K)+RHO(I-1,J,K)-RHO(I,J,K-1)-RHO(I-1,J,K-1))
C
      DO 360 K=1,KBM1
      DO 360 J=2,JMM1
      DO 360 I=2,IMM1
 360  DRHOX(I,J,K)=.25*(DT(I,J)+DT(I-1,J))*DRHOX(I,J,K)*DUM(I,J)
     .     *(H2(I,J)+H2(I-1,J))
C
C-------- Y COMPONENT OF BAROCLINIC PRESSURE GRADIENT ------------------
      DO 500 J=3,JMM1
      DO 500 I=2,IMM1
 500  DRHOY(I,J,1)=.25*GRAV*DZ(1)*(DT(I,J)+DT(I,J-1))*
     .                 (RHO(I,J,1)-RHO(I,J-1,1))
      DO 510 K=2,KBM1
      DO 510 J=3,JMM1
      DO 510 I=2,IMM1
 510  DRHOY(I,J,K)=DRHOY(I,J,K-1)
     .        +.125*GRAV*(DZ(K)+DZ(K-1))*(DT(I,J)+DT(I,J-1))
     .        *(RHO(I,J,K)-RHO(I,J-1,K)+RHO(I,J,K-1)-RHO(I,J-1,K-1))
     .        +.5*GRAV*Z(K)*(DT(I,J)-DT(I,J-1))
     .        *(RHO(I,J,K)+RHO(I,J-1,K)-RHO(I,J,K-1)-RHO(I,J-1,K-1))
C
      DO 560 K=1,KBM1
      DO 560 J=3,JMM1
      DO 560 I=2,IMM1
 560  DRHOY(I,J,K)=.25*(DT(I,J)+DT(I,J-1))*DRHOY(I,J,K)*DVM(I,J)
     .       *(H1(I,J)+H1(I,J-1))
      DO 580 K=1,KB
      DO 580 J=1,JM
      DO 580 I=1,IM
      DRHOX(I,J,K)=RAMP*DRHOX(I,J,K)*DUM(I,J)
 580  DRHOY(I,J,K)=RAMP*DRHOY(I,J,K)*DVM(I,J)
C
C-------- VERTICALLY INTEGRATE PRESSURE GRADIENT -----------------------
      DO 590 J=1,JM
      DO 590 I=1,IM
      TRNU(I,J)=0.0
 590  TRNV(I,J)=0.0
      DO 570 K=1,KBM1
      DO 570 J=1,JM
      DO 570 I=1,IM
      TRNU(I,J)=TRNU(I,J)+DRHOX(I,J,K)*DZ(K)
 570  TRNV(I,J)=TRNV(I,J)+DRHOY(I,J,K)*DZ(K)
C
C     DO 600 K=1,KB
C     DO 600 J=1,JM
C     DO 600 I=1,IM
C600  RHO(I,J,K)=RHO(I,J,K)+RMEAN(I,J,K)
C
      RETURN
      END
      

      SUBROUTINE PROFT(F,WFSURF,DT2)
C     VERSION(03/02/90)
      INCLUDE 'comdeck'
C
      REAL*8 AF(IM,JM,KB),CF(IM,JM,KB),VHF(IM,JM,KB),VHPF(IM,JM,KB)
      REAL*8 FF(IM,JM,KB)
C
      DIMENSION F(IM,JM,KB),WFSURF(IM,JM),DH(IM,JM)
      EQUIVALENCE (A,AF),(PROD,CF),(VH,VHF),(AF,FF)
      EQUIVALENCE (TPS,DH)
C
      UMOLPR=UMOL
C
C-----------------------------------------------------------------------
C
C        THE FOLLOWING SECTION SOLVES THE EQUATION
C         DT2*(KH*F')'-F=-FB
C
C-----------------------------------------------------------------------
C
      DO 90 J=2,JMM1
      DO 90 I=2,IMM1
      DH(I,J)=H(I,J)+ETF(I,J)
  90  CONTINUE
C
      DO 100 K=2,KBM1
      DO 100 J=2,JMM1
      DO 100 I=2,IMM1
      AF(I,J,K-1)=-DT2*(KH(I,J,K)+UMOLPR)/(DZ(K-1)*DZZ(K-1)*DH(I,J)
     .     *DH(I,J))
      CF(I,J,K)=-DT2*(KH(I,J,K)+UMOLPR)/(DZ(K)*DZZ(K-1)*DH(I,J)
     .     *DH(I,J))
 100  CONTINUE
C
C-------- SURFACE BOUNDARY CONDITIONS - WFSURF -------------------------
      DO 110 J=2,JMM1
      DO 110 I=2,IMM1
      VHF(I,J,1)=AF(I,J,1)/(AF(I,J,1)-1.)
      VHPF(I,J,1)=-DT2*WFSURF(I,J)/(-DZ(1)*DH(I,J))-F(I,J,1)
 110  VHPF(I,J,1)=VHPF(I,J,1)/(AF(I,J,1)-1.)
C
C-------- SURFACE BOUNDARY CONDITIONS - FSURF --------------------------
C     DO 110 J=2,JMM1
C     DO 110 I=2,IMM1
C     VHF(I,J,1)=0.
C110  VHPF(I,J,1)=FSURF(I,J)
C
      DO 101 K=2,KBM2
      DO 101 J=2,JMM1
      DO 101 I=2,IMM1
      VHPF(I,J,K)=1./(AF(I,J,K)+CF(I,J,K)*(1.-VHF(I,J,K-1))-1.)
      VHF(I,J,K)=AF(I,J,K)*VHPF(I,J,K)
      VHPF(I,J,K)=(CF(I,J,K)*VHPF(I,J,K-1)-DBLE(F(I,J,K)))*VHPF(I,J,K)
 101  CONTINUE
C
      DO 130 K=1,KB
      DO 130 J=1,JM
      DO 130 I=1,IM
 130  FF(I,J,K)=F(I,J,K)
C
      DO 102 J=2,JMM1
      DO 102 I=2,IMM1
 102  FF(I,J,KBM1)=((CF(I,J,KBM1)*VHPF(I,J,KBM2)-FF(I,J,KBM1))
     .   /(CF(I,J,KBM1)*(1.-VHF(I,J,KBM2))-1.))
C
      DO 105 K=2,KBM1
      KI=KB-K
      DO 105 J=2,JMM1
      DO 105 I=2,IMM1
      FF(I,J,KI)=(VHF(I,J,KI)*FF(I,J,KI+1)+VHPF(I,J,KI))
 105  CONTINUE
C
      DO 140 K=1,KB
      DO 140 J=1,JM
      DO 140 I=1,IM
 140  F(I,J,K)=FF(I,J,K)
C
      RETURN
      END
      

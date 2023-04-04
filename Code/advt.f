      SUBROUTINE ADVT(FB,F,FMEAN,DTI2,FF,FDIFF)
C     VERSION(09/29/91)
      INCLUDE 'comdeck'  
C
C-----------------------------------------------------------------------
C     THIS SUBROUTINE INTEGRATES CONSERVATIVE CONSTITUENT EQUATIONS
C-----------------------------------------------------------------------
C
      DIMENSION FB(IM,JM,KB),F(IM,JM,KB),FF(IM,JM,KB),FMEAN(IM,JM,KB)
      DIMENSION FDIFF(DBCM)
      DIMENSION XFLUX(IM,JM,KB),YFLUX(IM,JM,KB)
      DIMENSION AAMAX(IM,JM,KB),AAMAY(IM,JM,KB)
      EQUIVALENCE (XFLUX,A),(YFLUX,C)
      EQUIVALENCE (AAMAX,VH),(AAMAY,VHP)
C
      DO 10 J=1,JM
      DO 10 I=1,IM
      F (I,J,KB)=F (I,J,KBM1)
   10 FB(I,J,KB)=FB(I,J,KBM1)
C
C-------- DO ADVECTION FLUXES ------------------------------------------
      DO 20 K=1,KBM1
      DO 20 J=2,JM
      DO 20 I=2,IM
      XFLUX(I,J,K)=0.5*(F(I-1,J,K)+F(I,J,K)) * XMFL3D(I,J,K)
   20 YFLUX(I,J,K)=0.5*(F(I,J-1,K)+F(I,J,K)) * YMFL3D(I,J,K)
C
C-------- ADD DIFFUSIVE FLUXES -----------------------------------------
C     DO 30 K=1,KB
C     DO 30 J=1,JM
C     DO 30 I=1,IM
C  30 FB(I,J,K)=FB(I,J,K)-FMEAN(I,J,K)
      DO 40 K=1,KBM1
      DO 40 J=2,JM
      DO 40 I=2,IM
      AAMAX(I,J,K)=.5*(AAM(I,J,K)+AAM(I-1,J,K))
      AAMAY(I,J,K)=.5*(AAM(I,J,K)+AAM(I,J-1,K))
   40 CONTINUE
      DO 50 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      IF(JD.EQ.JC) THEN
            IF(IC.GT.ID) THEN
              DO 60 K=1,KBM1
   60         AAMAX(IC,JC,K)=0.0
            ELSE
              DO 70 K=1,KBM1
   70         AAMAX(ID,JD,K)=0.0
            ENDIF
      ELSE
            IF(JC.GT.JD) THEN
              DO 80 K=1,KBM1
   80         AAMAY(IC,JC,K)=0.0
            ELSE
              DO 90 K=1,KBM1
   90         AAMAY(ID,JD,K)=0.0
            ENDIF
      ENDIF
   50 CONTINUE
      DO 100 K=1,KBM1
      DO 100 J=2,JM
      DO 100 I=2,IM
      XFLUX(I,J,K)=XFLUX(I,J,K)
     .    -AAMAX(I,J,K)/HPRNU*(H(I,J)+H(I-1,J))
     .    *(FB(I,J,K)-FB(I-1,J,K))*DUM(I,J)/(H1(I,J)+H1(I-1,J))
     .    *0.5*(H2(I,J)+H2(I-1,J))
      YFLUX(I,J,K)=YFLUX(I,J,K)
     .    -AAMAY(I,J,K)/HPRNU*(H(I,J)+H(I,J-1))
     .    *(FB(I,J,K)-FB(I,J-1,K))*DVM(I,J)/(H2(I,J)+H2(I,J-1))
     .    *0.5*(H1(I,J)+H1(I,J-1))
  100 CONTINUE
C     DO 110 K=1,KB
C     DO 110 J=1,JM
C     DO 110 I=1,IM
C 110 FB(I,J,K)=FB(I,J,K)+FMEAN(I,J,K)
C
C-------- DO VERTICAL ADVECTION ----------------------------------------
      DO 120 J=2,JMM1
      DO 120 I=2,IMM1
  120 FF(I,J,1)=-0.5*DZR(1)*(F(I,J,1)+F(I,J,2))*W(I,J,2)*ART(I,J)
      DO 130 K=2,KBM1
      DO 130 J=2,JMM1
      DO 130 I=2,IMM1
  130 FF(I,J,K)=0.5*DZR(K)*((F(I,J,K-1)+F(I,J,K))*W(I,J,K)
     .                    -(F(I,J,K)+F(I,J,K+1))*W(I,J,K+1))*ART(I,J)
C
C-------- ADD NET HORIZONTAL FLUXES & THEN STEP FORWARD IN TIME --------
      DO 140 K=1,KBM1
      DO 140 J=2,JMM1
      DO 140 I=2,IMM1
      FF(I,J,K)=FF(I,J,K)
     .              +XFLUX(I+1,J,K)-XFLUX(I,J,K)
     .              +YFLUX(I,J+1,K)-YFLUX(I,J,K)
  140 FF(I,J,K)=(FB(I,J,K)*(H(I,J)+ETB(I,J))*ART(I,J)-DTI2*FF(I,J,K))
     .                    /((H(I,J)+ETF(I,J))*ART(I,J))
C
C-----------------------------------------------------------------------
C         IMPOSE TRACER FLUX BOUNDARY CONDITIONS
C         REVISED 7-23-98 by MSV TO CORRECT
C         if QDIFF(N) is pos. use FDIFF(N)
C         if QDIFF(N) is neg. use FF(ID,JD,K)
C-----------------------------------------------------------------------
C
      DO 150 N=1,NUMDBC
      ID=IDD(N)
      JD=JDD(N)
      DO 150 K=1,KBM1
      if(QDIFF(N).ge.0.0) then
      FF(ID,JD,K)=FF(ID,JD,K)+DTI2*FDIFF(N)*QDIFF(N)*RAMP*
     .   VDDIST(N,K)/100./((H(ID,JD)+ETF(ID,JD))*ART(ID,JD)
     .  *DZ(K))
      else
C+++ case added for withdrawals by MSV 7-23-98
      FF(ID,JD,K)=FF(ID,JD,K)+DTI2*FF(ID,JD,K)*QDIFF(N)*RAMP*
     .   VDDIST(N,K)/100./((H(ID,JD)+ETF(ID,JD))*ART(ID,JD)
     .  *DZ(K))
C+++
      endif
  150 CONTINUE
C
      RETURN
      END
      

      SUBROUTINE VERTVL(DTI2)
C     VERSION(03/02/90)
      INCLUDE 'comdeck'
C
      DIMENSION XFLUX(IM,JM,KB),YFLUX(IM,JM,KB),ZMFL3D(IM,JM,KB)
      EQUIVALENCE (XFLUX,A),(YFLUX,C),(ZMFL3D,VHP)
C
C-------- CALCULATE NEW VERTICAL VELOCITY ------------------------------
      DO 100 J=1,JM
      DO 100 I=1,IM
 100  W(I,J,1)=0.0  
C
      DO 101 K=1,KB
      DO 101 J=1,JM
      DO 101 I=1,IM
 101  ZMFL3D(I,J,K)=0.0
C
C-----------------------------------------------------------------------
C         IMPOSE MASS FLUX BOUNDARY CONDITIONS
C-----------------------------------------------------------------------
      DO 120 N=1,NUMDBC
      ID=IDD(N)
      JD=JDD(N)
      DO 120 K=1,KBM1
 120  ZMFL3D(ID,JD,K)=QDIFF(N)*RAMP*(VDDIST(N,K)/100.0)/DZ(K)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C         IMPOSE MASS FLUX BOUNDARY CONDITIONS FOR PREC. AND EVAP.
C         IN SURFACE LAYER K=1
C         MSV 7-22-98
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO 122 K=1 
      DO 122 J=2,JMM1
      DO 122 I=2,IMM1
      ZMFL3D(I,J,K)=ZMFL3D(I,J,K)
     &  +(QPREC-QEVAP)*ART(I,J)*RAMP*1.0/DZ(K)
 122  CONTINUE 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DO 110 K=1,KBM1
      DO 110 J=2,JMM1
      DO 110 I=2,IMM1
 110  W(I,J,K+1)=W(I,J,K)
     .    +DZ(K)*((-ZMFL3D(I,J,K)+XMFL3D(I+1,J,K)-XMFL3D(I,J,K)
     .            +YMFL3D(I,J+1,K)-YMFL3D(I,J,K))/ART(I,J)
     .                        +(ETF(I,J)-ETB(I,J))/DTI2 )
C
C     DO 121 N=1,NUMDBC
C     ID=IDD(N)
C     JD=JDD(N)
C     DO 121 K=1,KB
C     print *, id,jd,k,xmfl3d(id,jd,k),ymfl3d(id,jd,k)
C     print *, k,xmfl3d(id+1,jd,k),ymfl3d(id,jd+1,k)
C     print *, k,xmfl3d(id+1,jd,k)-xmfl3d(id,jd,k)
C     print *, k,ymfl3d(id,jd+1,k)-ymfl3d(id,jd,k)
C     print *, k,etf(id,jd), etb(id,jd), zmfl3d(id,jd,k)           
C 121 print *, h1(id,jd)*h2(id,jd),dti2,dz(k),w(id,jd,k)
C
      DO 130 J=1,JM
      DO 130 I=1,IM
 130  W(I,J,KB)=0.0  
C
      RETURN
      END
      

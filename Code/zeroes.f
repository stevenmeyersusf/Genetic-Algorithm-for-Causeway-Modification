      SUBROUTINE ZEROES(ADVUA,ADVVA,ADVUU,ADVVV,DRHOX,DRHOY,TRNU,TRNV)
C     VERSION(03/02/90)
      INCLUDE 'comdeck'
C
      DIMENSION ADVUA(IM,JM),ADVVA(IM,JM),ADVUU(IM,JM),ADVVV(IM,JM)
      DIMENSION DRHOX(IM,JM,KB),DRHOY(IM,JM,KB),TRNU(IM,JM),TRNV(IM,JM)
C
      DO 250 K=1,KB
      DO 250 J=1,JM
      DO 250 I=1,IM
      T(I,J,K)=0.0
      TB(I,J,K)=0.0
      S(I,J,K)=0.0
      SB(I,J,K)=0.0
      CCC(I,J,K)=0.0
      CCCB(I,J,K)=0.0
      RHO(I,J,K)=0.0
      RMEAN(I,J,K)=0.0
      TMEAN(I,J,K)=0.0
      SMEAN(I,J,K)=0.0
      CCCMEAN(I,J,K)=0.0
      XMFL3D(I,J,K)=0.0
      YMFL3D(I,J,K)=0.0
      A(I,J,K)=0.0
      C(I,J,K)=0.0
      VH(I,J,K)=0.0
      VHP(I,J,K)=0.0
      PROD(I,J,K)=0.0
      DTEF(I,J,K)=0.0
      U(I,J,K)=0.0
      UB(I,J,K)=0.0
      V(I,J,K)=0.0
      VB(I,J,K)=0.0
      UF(I,J,K)=0.0
      VF(I,J,K)=0.0
      CCCF(I,J,K)=0.0
      W(I,J,K)=0.0
      DRHOX(I,J,K)=0.0
      DRHOY(I,J,K)=0.0
      Q2B(I,J,K)=0.0
      Q2(I,J,K)=0.0
      Q2LB(I,J,K)=0.0
      Q2L(I,J,K)=0.0
      L(I,J,K)=0.0
      KH(I,J,K)=0.0
      KM(I,J,K)=0.0
      KQ(I,J,K)=0.0
  250 CONTINUE
      DO 640 J=1,JM
      DO 640 I=1,IM
      UAF(I,J)=0.0
      UA(I,J)=0.0
      UAB(I,J)=0.0
      VAF(I,J)=0.0
      VA(I,J)=0.0
      VAB(I,J)=0.0
      ELF(I,J)=0.0
      EL(I,J)=0.0
      ELB(I,J)=0.0
      TPS(I,J)=0.0
      ETF(I,J)=0.0
      ET (I,J)=0.0
      ETB(I,J)=0.0
      UTF(I,J)=0.0
      UTB(I,J)=0.0
      VTF(I,J)=0.0
      VTB(I,J)=0.0
      EGF(I,J)=0.0
      EGB(I,J)=0.0
      WTSURF(I,J)=0.0
      WSSURF(I,J)=0.0
      WCCCSURF(I,J)=0.0
      WUSURF(I,J)=0.0
      WVSURF(I,J)=0.0
      WUBOT(I,J)=0.0
      WVBOT(I,J)=0.0
      TRNU(I,J)=0.0
      TRNV(I,J)=0.0
      CURV42D(I,J)=0.0
      ADVUA(I,J)=0.0
      ADVVA(I,J)=0.0
      ADVUU(I,J)=0.0
      ADVVV(I,J)=0.0
      FLUXUA(I,J)=0.0
      FLUXVA(I,J)=0.0
  640 CONTINUE
C
      DO 100 K=1,KB
      DO 100 J=1,JM
      DO 100 I=1,IM
      ARCU (I,J,K)=0.0
      ARCV (I,J,K)=0.0
      ARCUX(I,J,K)=0.0
      ARCVX(I,J,K)=0.0
      ARCS (I,J,K)=0.0
      ARCT (I,J,K)=0.0
      ARCCCC (I,J,K)=0.0
      ARCW (I,J,K)=0.0
      ARCKH(I,J,K)=0.0
 100  CONTINUE
      DO 110 J=1,JM
      DO 110 I=1,IM
 110  ARCET (I,J)=0.0
      DO 120 K=1,KB
      DO 120 J=1,JM
      DO 120 I=1,IM
      ARCTU (I,J,K)=0.0
      ARCTV (I,J,K)=0.0
      ARCTW(I,J,K)=0.0
      ARCTAAM(I,J,K)=0.0
      ARCTKH(I,J,K)=0.0
 120  CONTINUE
      DO 130 J=1,JM
      DO 130 I=1,IM
      ARCTES(I,J)=0.0
      ARCTED(I,J)=0.0
 130  CONTINUE
C
      DO 200 N=1,EPTS
      ESAVE(N)=0.0
 200  CONTINUE
      DO 210 N=1,VPTS
      DZSAVE(N)=0.0
      DO 210 K=1,KB
      UZSAVE(N,K)=0.0
      VZSAVE(N,K)=0.0
      SZSAVE(N,K)=0.0
      TZSAVE(N,K)=0.0
      CCCZSAVE(N,K)=0.0
 210  CONTINUE
      DO 220 K=1,KB
      DO 220 N=1,FPTS
      CCFLUX(N,K)=0.0
 220  CONTINUE
      ESUM=0.0
      TKE =0.0
      APE =0.0
      TSUM=0.0
      SSUM=0.0
      CCCSUM=0.0
      VART=0.0
      VARS=0.0
      VARCCC=0.0
C
C ------- SET CONSPLT AND CONSTSR TO TRUE FOR FIRST RUN THROUGH --------
      CONSPLT=.TRUE.
      CONSTSR=.TRUE.
C
      RETURN
      END

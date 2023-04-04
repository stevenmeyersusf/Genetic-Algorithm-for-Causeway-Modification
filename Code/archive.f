      SUBROUTINE ARCHIVE(IDX)
C     VERSION(11/18/90)
C
      INCLUDE 'comdeck'
CDCB modified 10/1/00 to include new point files for
c instantanious vel, el, and s files
CMSV++++ new character for tampa bay project to use to open files
CMSV++++ for individual tseries
      character*7 efile
      CHARACTER*87 catdatapath
C+++++++
C
      GO TO (10,20,30), IDX
C
C-------- TIME AVERAGING FOR 3-D AESOP COMPUTATIONS --------------------
 10   CONTINUE
      DO 8301 JTRAC=1,JTM
      IF(INT.LT.ITRAC(JTRAC,1).OR.INT.GT.ITRAC(JTRAC,2)) GO TO 8301
      DO 501 K=1,KB
      DO 501 J=2,JMM1
      DO 501 I=2,IM
      ARCTU(I,J,K)=ARCTU(I,J,K)+U(I,J,K)*DEIT*(DT(I,J)+DT(I-1,J))/2.
 501  CONTINUE
      DO 502 K=1,KB
      DO 502 J=2,JM
      DO 502 I=2,IMM1
      ARCTV(I,J,K)=ARCTV(I,J,K)+V(I,J,K)*DEIT*(DT(I,J)+DT(I,J-1))/2.
 502  CONTINUE
      DO 503 K=1,KB
      DO 503 J=2,JMM1
      DO 503 I=2,IMM1
      ARCTW(I,J,K)=ARCTW(I,J,K)+W(I,J,K)*DEIT
      ARCTAAM(I,J,K)=ARCTAAM(I,J,K)+AAM(I,J,K)/HPRNU*DEIT
      ARCTKH(I,J,K)=ARCTKH(I,J,K)+KH(I,J,K)*DEIT
 503  CONTINUE
      DO 504 J=1,JM
      DO 504 I=1,IM
      ARCTED(I,J)=ARCTED(I,J)+(ETF(I,J)-ETB(I,J))/(2.*DTI)*DEIT
 504  CONTINUE
C
      IF (INT.EQ.ITRAC(JTRAC,1)) THEN
       DO 523 J=1,JM
       DO 523 I=1,IM
523    ARCTES(I,J)=(ETB(I,J)+ET(I,J))/2.
      ENDIF
C
      IF(INT.NE.ITRAC(JTRAC,2)) GO TO 561
C
      TMIDDLE=TIME-(0.5*DTI*DAYI/DEIT)
      WRITE(IUFLW)TMIDDLE,ARCTU,ARCTV,ARCTW,ARCTAAM,ARCTKH,ARCTES,ARCTED
C
      DO 571 K=1,KB
      DO 571 J=1,JM
      DO 571 I=1,IM
      ARCTU (I,J,K)=0.0
      ARCTV (I,J,K)=0.0
      ARCTW(I,J,K)=0.0
      ARCTAAM(I,J,K)=0.0
      ARCTKH(I,J,K)=0.0
 571  CONTINUE
      DO 572 J=1,JM
      DO 572 I=1,IM
      ARCTES(I,J)=0.0
      ARCTED(I,J)=0.0
 572  CONTINUE
 561  CONTINUE
C
 8301 CONTINUE
      RETURN

C**************************************************************
C**************************************************************

C-------- TIME AVERAGING FOR RESIDENCE TIME CONCENTRATIONS  COMPUTATIONS
 30   CONTINUE
c      IF (iccctime .ne. 0) then
c
c      DO 601 K=1,KB
c      DO 601 J=2,JMM1
c      DO 601 I=2,IMM1
c      ARCTCCC(I,J,K)=ARCTCCC(I,J,K)+CCC(I,J,K)*DEITCCC
c 601  CONTINUE
cC
c      IF(MOD(INT,iccctime).EQ.0) then
cC
c      TMIDDLE=TIME-(0.5*DTI*DAYI/DEITCCC)
cCDCB set to write after spinup *************
cc      if (thour .gt. 6550.0) then
cc      WRITE(IUCCCRES,101) TMIDDLE
cc      WRITE(IUCCCRES,101) (((ARCTCCC (I,J,K),I=1,IM),J=1,JM),K=1,KB)
cc      endif
cc      WRITE(IUFLW) TMIDDLE
cc      WRITE(IUFLW) ARCTCCC
cC
cC      print *,tmiddle,arctccc(35,35,5),int,iccctime
c      DO 671 K=1,KB
c      DO 671 J=1,JM
c      DO 671 I=1,IM
c      ARCTCCC (I,J,K)=0.0
c 671  CONTINUE
c      ENDIF
cC
c      ENDIF
      RETURN
C***************************************************************
C**************************************************************



C-----------------------------------------------------------------------
C-------- COMPUTATIONAL HISTORY WRITES AND ACCUMULATION ----------------
 20   CONTINUE
      DO 8300 JHIST=1,JHM
      IF(INT.LT.IHIST(JHIST,1).OR.INT.GT.IHIST(JHIST,2)) GO TO 8300
      IF(TOR.EQ.'PROGNOSTIC' .OR. TOR.EQ.'DIAGNOSTIC') THEN
       DO 500 K=1,KB
       DO 500 J=2,JMM1
       DO 500 I=2,IMM1
       ARCU (I,J,K)=ARCU (I,J,K)+U (I,J,K)*DEI
       ARCV (I,J,K)=ARCV (I,J,K)+V (I,J,K)*DEI
       ARCUX(I,J,K)=ARCUX(I,J,K)+U (I,J,K)*DEI*
     2    0.5*(H(I,J)+ET(I,J)+H(I-1,J)+ET(I-1,J))*DZ(K)
       ARCVX(I,J,K)=ARCVX(I,J,K)+V (I,J,K)*DEI*
     2    0.5*(H(I,J)+ET(I,J)+H(I,J-1)+ET(I,J-1))*DZ(K)
       ARCS (I,J,K)=ARCS (I,J,K)+S (I,J,K)*DEI
       ARCT (I,J,K)=ARCT (I,J,K)+T (I,J,K)*DEI
       ARCCCC (I,J,K)=ARCCCC (I,J,K)+CCC (I,J,K)*DEI
       ARCW (I,J,K)=ARCW (I,J,K)+WR(I,J,K)*DEI
       ARCKH(I,J,K)=ARCKH(I,J,K)+KH(I,J,K)*DEI
 500   CONTINUE
       DO 520 J=1,JM
       DO 520 I=1,IM
 520   ARCET (I,J)=ARCET (I,J)+ET (I,J)*DEI
      ELSE
       DO 530 J=1,JM
       DO 530 I=1,IM
       ARCUX(I,J,1)=ARCUX(I,J,1)+UAF(I,J)*DEI*
     2    0.5*(H(I,J)+EL(I,J)+H(I-1,J)+EL(I-1,J))
       ARCVX(I,J,1)=ARCVX(I,J,1)+VAF(I,J)*DEI*
     2    0.5*(H(I,J)+EL(I,J)+H(I,J-1)+EL(I,J-1))
       ARCU (I,J,1)=ARCU (I,J,1)+UAF(I,J)*DEI 
       ARCV (I,J,1)=ARCV (I,J,1)+VAF(I,J)*DEI 
 530   ARCET (I,J)=ARCET (I,J)+EL (I,J)*DEI
      END IF
C
      IF(INT.NE.IHIST(JHIST,2)) GO TO 560
C
C-------- WRITE CONSTANTS FIRST TIME THROUGH ---------------------------
      IF (CONSPLT) THEN
       WRITE(IUPLT,103) DTI,GRAV,UMOL,TOR
       WRITE(IUPLT,102) NUMEBC
       WRITE(IUPLT,102) (IETA(I),JETA(I),ICON(I),JCON(I),I=1,NUMEBC)
       WRITE(IUPLT,102) NUMQBC
       WRITE(IUPLT,102) (IQC(I),JQC(I),I=1,NUMQBC)
       WRITE(IUPLT,101) ((H(I,J),I=1,IM),J=1,JM)
       WRITE(IUPLT,101) ((H1(I,J),I=1,IM),J=1,JM)
       WRITE(IUPLT,101) ((H2(I,J),I=1,IM),J=1,JM)
       WRITE(IUPLT,101) ((ANG(I,J),I=1,IM),J=1,JM)
       WRITE(IUPLT,101) ((DUM(I,J),I=1,IM),J=1,JM)
       WRITE(IUPLT,101) ((DVM(I,J),I=1,IM),J=1,JM)
       WRITE(IUPLT,101) ((FSM(I,J),I=1,IM),J=1,JM)
       CONSPLT=.FALSE.
      ENDIF
C
      TMIDDLE=TIME-(0.5*DTI*DAYI/DEI)
      WRITE(IUPLT,101) TMIDDLE
      WRITE(IUPLT,101) ((ARCET(I,J),I=1,IM),J=1,JM)
C
      IF(TOR.EQ.'BAROTROPIC') THEN
       WRITE(IUPLT,101) ((ARCU (I,J,1),I=1,IM),J=1,JM)
       WRITE(IUPLT,101) ((ARCV (I,J,1),I=1,IM),J=1,JM)
       WRITE(IUPLT,101) ((ARCUX(I,J,1),I=1,IM),J=1,JM)
       WRITE(IUPLT,101) ((ARCVX(I,J,1),I=1,IM),J=1,JM)
      ELSE
       WRITE(IUPLT,101) (Z(K),K=1,KB)
       WRITE(IUPLT,101) (ZZ(K),K=1,KB)
       WRITE(IUPLT,101) (DZ(K),K=1,KB)
       WRITE(IUPLT,101) (((ARCU (I,J,K),I=1,IM),J=1,JM),K=1,KB)
       WRITE(IUPLT,101) (((ARCV (I,J,K),I=1,IM),J=1,JM),K=1,KB)
       WRITE(IUPLT,101) (((ARCUX(I,J,K),I=1,IM),J=1,JM),K=1,KB)
       WRITE(IUPLT,101) (((ARCVX(I,J,K),I=1,IM),J=1,JM),K=1,KB)
       WRITE(IUPLT,101) (((ARCT (I,J,K),I=1,IM),J=1,JM),K=1,KB)
       WRITE(IUPLT,101) (((ARCS (I,J,K),I=1,IM),J=1,JM),K=1,KB)
       WRITE(IUPLT,101) (((ARCCCC (I,J,K),I=1,IM),J=1,JM),K=1,KB)
CCCT       Write(IUPLT,101) (((ARCCCCT (I,J,K),I=1,IM),J=1,JM),K=1,KB)
       WRITE(IUPLT,101) (((ARCW (I,J,K),I=1,IM),J=1,JM),K=1,KB)
       WRITE(IUPLT,101) (((ARCKH(I,J,K),I=1,IM),J=1,JM),K=1,KB)
      END IF
C
 101  FORMAT(8E12.5)
 102  FORMAT(8I10)
 103  FORMAT(3E12.5,A10)
C
      DO 570 K=1,KB
      DO 570 J=1,JM
      DO 570 I=1,IM
      ARCU (I,J,K)=0.0
      ARCV (I,J,K)=0.0
      ARCUX(I,J,K)=0.0
      ARCVX(I,J,K)=0.0
      ARCS (I,J,K)=0.0
      ARCCCC (I,J,K)=0.0
      ARCT (I,J,K)=0.0
      ARCW (I,J,K)=0.0
      ARCKH(I,J,K)=0.0
 570  CONTINUE
      DO 590 J=1,JM
      DO 590 I=1,IM
 590  ARCET (I,J)=0.0
 560  CONTINUE
C
 8300 CONTINUE


C
C-----------------------------------------------------------------------
C------- TIME SERIES WRITES AND ACCUMULATIONS --------------------------
      IF(TOR.EQ.'BAROTROPIC') THEN
       DO 140 N=1,EPTS
       II=INXIE(N)
       JJ=INXJE(N)
       ESAVE(N)=ESAVE(N)+EL(II,JJ)*SKILLI
140    CONTINUE
       DO 150 N=1,VPTS
       II=INXIV(N)
       JJ=INXJV(N)
       UZSAVE(N,1)=UZSAVE(N,1)+.5*(UA(II,JJ)+UA(II+1,JJ))*SKILLI
       VZSAVE(N,1)=VZSAVE(N,1)+.5*(VA(II,JJ)+VA(II,JJ+1))*SKILLI
150    CONTINUE
C
C-------- COMPUTE CROSS SECTIONAL FLUXES -------------------------------
       DO 180 N=1,FPTS
       IS=ISFLX(N)
       JS=JSFLX(N)
       IF(DIRFLX(N).EQ.'IDIR') THEN
        IE=IS+NFLXE(N)-1
        DO 181 I=IS,IE
 181    CCFLUX(N,1)=CCFLUX(N,1)+VAF(I,JS)*SKILLI*0.25*
     2     (H(I,JS)+EL(I,JS)+H(I,JS-1)+EL(I,JS-1))*(H1(I,JS)+H1(I,JS-1))
       ELSE
C        DIRFLX(N).EQ.'JDIR'
        JE=JS+NFLXE(N)-1
        DO 182 J=JS,JE
 182    CCFLUX(N,1)=CCFLUX(N,1)+UAF(IS,J)*SKILLI*0.25*
     2     (H(IS,J)+EL(IS,J)+H(IS-1,J)+EL(IS-1,J))*(H2(IS,J)+H2(IS-1,J))
       ENDIF
 180   CONTINUE
C
C-------- TOR = PROGNOSTIC OR DIAGNOSTIC -------------------------------
      ELSE
       DO 160 N=1,EPTS
       II=INXIE(N)
       JJ=INXJE(N)
cdcb old pointfile time series output
cC+++++++++++++++++++++++++++++++++++++++++++++++++
cCMSV++++ new code to open individual elev. tseries files
c       if(ifpointfile.eq.1) then
c        efile='e'
c        write(efile(2:4),194) INXIE(N)
c        write(efile(5:7),194) INXJE(N)
c        if(ifopenede.eq.0) then
cC++++ find out how long a string datapath is (non blanks)
c        do 158 ic=1,80
c        if(datapath(ic:ic).eq.' ') then
c        icc=ic-1
c        goto 159
c        endif
c  158   continue
c  159   continue
c        print *,'datapath',datapath(1:ic)
cC       catdatapath=datapath(1:ic)//efile
cC       print *,catdatapath
c        open
c     &  (19+n,form='formatted',file=datapath(1:icc)//efile,
c     &  access='append')
c        endif
cC++++now write to the individual tseries file
c        write(19+N,196) efile, time, ET(II,JJ)
c       endif
cCMSV++++++++++++++++++++++++++++++++++++++++++++++
cC+++++++++++++++++++++++++++++++++++++++++++++++++
cC++++ next two lines are from original code

       ESAVE(N)=ESAVE(N)+ET(II,JJ)*SKILLI
160    CONTINUE

cc       ifopenede=1

       DO 170 N=1,VPTS
       II=INXIV(N)
       JJ=INXJV(N)
cC+++++++++++++++++++++++++++++++++++++++++++++++++
cCMSV++++ new code to open individual tseries files
cCMSV++++ set to 1 to tell e tseries are opened
c       if(ifpointfile.eq.1) then
cC      print *,'in vtseries'
cCMSV++++
cC++++ use same character (temp) for v files as for u files
c        efile='v'
c        write(efile(2:4),194) INXIV(N)
c        write(efile(5:7),194) INXJV(N)
cC       print *,'efile',efile
cC       print *,'ifopenedv',ifopenedv
c        if(ifopenedv.eq.0) then
cC++++ find out how long a string datapath is (non blanks)
c        do 168 ic=1,80
c        if(datapath(ic:ic).eq.' ') then
c        icc=ic-1
c        goto 169
c        endif
c  168   continue
c  169   continue
cC       print *,datapath
c        open
c     &  (29+n,form='formatted',file=datapath(1:icc)//efile,
c     &  access='append')
c        endif
cC++++now write to the individual tseries file
c        write(29+N,197) efile, time, ang(ii,jj), h(ii,jj),ET(II,JJ)
c        write(29+N,198) (((u(i,j,k),i=ii,ii),j=jj,jj),k=1,kb)
c        write(29+N,198) (((v(i,j,k),i=ii,ii),j=jj,jj),k=1,kb)
c       endif
cCMSV++++++++++++++++++++++++++++++++++++++++++++++
cC+++++++++++++++++++++++++++++++++++++++++++++++++

       DZSAVE(N)=DZSAVE(N)+DT(II,JJ)*SKILLI
       DO 170 K=1,KB
Cdcb new arrays for centered u,v ts pointfiles

       ucenter(N,K)= .5*(U(II,JJ,K)+U(II+1,JJ,K))
       vcenter(N,K)= .5*(V(II,JJ,K)+V(II,JJ+1,K))

       UZSAVE(N,K)=UZSAVE(N,K)+.5*(U(II,JJ,K)+U(II+1,JJ,K))*SKILLI
       VZSAVE(N,K)=VZSAVE(N,K)+.5*(V(II,JJ,K)+V(II,JJ+1,K))*SKILLI
       SZSAVE(N,K)=SZSAVE(N,K)+S(II,JJ,K)*SKILLI
       TZSAVE(N,K)=TZSAVE(N,K)+T(II,JJ,K)*SKILLI
       CCCZSAVE(N,K)=CCCZSAVE(N,K)+T(II,JJ,K)*SKILLI
 170   CONTINUE
cCMSV++++
cC++++++ set to 1 to tell v tseries are opened
c       ifopenedv=1
cCMSV++++
C
C-------- COMPUTE CROSS SECTIONAL FLUXES -------------------------------
       DO 250 N=1,FPTS
       IS=ISFLX(N)
       JS=JSFLX(N)
       IF(DIRFLX(N).EQ.'IDIR') THEN
        IE=IS+NFLXE(N)-1
        DO 251 K=1,KBM1
        DO 251 I=IS,IE
 251    CCFLUX(N,K)=CCFLUX(N,K)+V(I,JS,K)*DZ(K)*SKILLI*0.25*
     2     (H(I,JS)+ET(I,JS)+H(I,JS-1)+ET(I,JS-1))*(H1(I,JS)+H1(I,JS-1))
       ELSE
C        DIRFLX(N).EQ.'JDIR'
        JE=JS+NFLXE(N)-1
        DO 252 K=1,KBM1
        DO 252 J=JS,JE
 252    CCFLUX(N,K)=CCFLUX(N,K)+U(IS,J,K)*DZ(K)*SKILLI*0.25*
     2     (H(IS,J)+ET(IS,J)+H(IS-1,J)+ET(IS-1,J))*(H2(IS,J)+H2(IS-1,J))
       ENDIF
 250   CONTINUE
      ENDIF
C
      IF(TOR.EQ.'BAROTROPIC') THEN
       DO 100 J=2,JMM1
       DO 100 I=2,IMM1
       ESUM=ESUM+EL(I,J)*H1(I,J)*H2(I,J)*FSM(I,J)/AREA*SKILLI
       APE=APE+.5*GRAV*EL(I,J)*EL(I,J)*FSM(I,J)*H1(I,J)*H2(I,J)/
     2     AREA*SKILLI
       TKE=TKE+0.125*D(I,J)*((UA(I,J)+UA(I+1,J))**2+
     2     (VA(I,J)+VA(I,J+1))**2)*FSM(I,J)*H1(I,J)*H2(I,J)/
     3     AREA*SKILLI
 100   CONTINUE
      ELSE
       DO 210 K=1,KBM1
       DO 210 J=2,JMM1
       DO 210 I=2,IMM1
       TSUM=TSUM+T(I,J,K)*DT(I,J)*(Z(K)-Z(K+1))*
     2      H1(I,J)*H2(I,J)*FSM(I,J)/VOLUME*SKILLI
       SSUM=SSUM+S(I,J,K)*DT(I,J)*(Z(K)-Z(K+1))*
     2      H1(I,J)*H2(I,J)*FSM(I,J)/VOLUME*SKILLI
       CCCSUM=CCCSUM+CCC(I,J,K)*DT(I,J)*(Z(K)-Z(K+1))*
     2      H1(I,J)*H2(I,J)*FSM(I,J)/VOLUME*SKILLI
210    CONTINUE
       DO 220 K=1,KBM1
       DO 220 J=2,JMM1
       DO 220 I=2,IMM1
       VART=VART+(T(I,J,K)-TB(I,J,K))**2*H1(I,J)*H2(I,J)*FSM(I,J)
       VARS=VARS+(S(I,J,K)-SB(I,J,K))**2*H1(I,J)*H2(I,J)*FSM(I,J)
       VARCCC=VARCCC+(CCC(I,J,K)-CCCB(I,J,K))**2*H1(I,J)*
     & H2(I,J)*FSM(I,J)
220    CONTINUE
      END IF
C
      IF(ISKILL.EQ.0.OR.MOD(INT,ISKILL).NE.0) GO TO 8250
C
C-------- WRITE CONSTANTS THE FIRST TIME THROUGH -----------------------
      IF (CONSTSR) THEN
       WRITE(IUTSR,191) TOR
       WRITE(IUTSR,192) EPTS
       WRITE(IUTSR,192) (INXIE(N),INXJE(N),N=1,EPTS)
       WRITE(IUTSR,192) VPTS
       WRITE(IUTSR,192) (INXIV(N),INXJV(N),N=1,VPTS)
       WRITE(IUTSR,192) FPTS
       WRITE(IUTSR,193) (ISFLX(N),JSFLX(N),DIRFLX(N),NFLXE(N),N=1,FPTS)
       CONSTSR=.FALSE.
      ENDIF
C
      TMIDDLE=TIME-(.5*DTI*DAYI/SKILLI)
      WRITE(IUTSR,101) TMIDDLE
C
      IF(TOR.EQ.'BAROTROPIC') THEN
       WRITE(IUTSR,190)  (ESAVE(N),N=1,EPTS)
       DO 195 N=1,VPTS
 195   WRITE(IUTSR,190)  UZSAVE(N,1),VZSAVE(N,1)
       WRITE(IUTSR,190)  (CCFLUX(N,1),N=1,FPTS)
       WRITE(IUTSR,190)  ESUM,TKE,APE
      ELSE
c       WRITE(IUTSR,190)  (ESAVE(N),N=1,EPTS)
c       WRITE(IUTSR,190)  (DZSAVE(N),N=1,VPTS)
       DO 200 K=1,KBM1
       DO 200 N=1,VPTS
c      WRITE(IUTSR,190)  UZSAVE(N,K),VZSAVE(N,K),SZSAVE(N,K),TZSAVE(N,K)
c     & ,CCCZSAVE(N,K)
 200      continue
c       WRITE(IUTSR,190)  ((CCFLUX(N,K),N=1,FPTS),K=1,KBM1)
C
       VART=SQRT(VART/(AREA*FLOAT(KBM1))*SKILLI)/DTI
       VARS=SQRT(VARS/(AREA*FLOAT(KBM1))*SKILLI)/DTI
c       WRITE(IUTSR,190) TSUM,SSUM,VART,VARS,VARCCC


CMSV++++++++++++++++++++++++++++++++++++++++++++++
CMSV++++ new code to open individual elev. tseries files
       DO 1140 N=1,EPTS
       II=INXIE(N)
       JJ=INXJE(N)
       if(ifpointfile.eq.1) then
        efile='e'
        write(efile(2:4),194) INXIE(N)
        write(efile(5:7),194) INXJE(N)

        if(ifopenede.eq.0) then
C++++ find out how long a string datapath is (non blanks)
        icc=index(datapath,' ')
        icc=icc-1
cc        print *,'datapath',datapath(1:icc)//efile
        open
     &(20+n,
     &form='formatted',file=datapath(1:icc)//efile,
     &  access='append')
        endif

C++++now write to the individual tseries file
c        write(20+N,196) efile, time ,ET(II,JJ)
c        write(*,196) efile(2:7), tmiddle ,ESAVE(N)
        write(20+N,196) efile(2:7), tmiddle ,ESAVE(N)        
       endif
1140   continue
       ifopenede=1
CMSV----------------------------------------------

CMSV++++++++++++++++++++++++++++++++++++++++++++++
CMSV++++ new code to open individual tseries files
CMSV++++ set to 1 to tell e tseries are opened
       DO 1170 N=1,VPTS
       II=INXIV(N)
       JJ=INXJV(N)
       if(ifpointfile.eq.1) then
C      print *,'in vtseries'
CMSV++++
C++++ use same character (temp) for v files as for u files
        efile='v'
        write(efile(2:4),194) INXIV(N)
        write(efile(5:7),194) INXJV(N)
C       print *,'efile',efile
C       print *,'ifopenedv',ifopenedv
        if(ifopenedv.eq.0) then
C++++ find out how long a string datapath is (non blanks)
        icc=index(datapath,' ')
        icc=icc-1

cc        print *,datapath(1:icc)//efile
        open
     &  (40+n,form='formatted',file=datapath(1:icc)//efile,
     &  access='append')
        endif
C++++now write to the individual tseries file
CMSV++++ use instant values for depth and et
CMSV++++ use NEW cell center instantaneous velocities
C234567
c      write(*,199)
c     &efile(2:7), tmiddle, ang(II,JJ), DZsave(N),Esave(N)
      write(40+N,199)
c     &efile, time, ang(II,JJ), DT(II,JJ),ET(II,JJ)
     &efile(2:7), tmiddle, ang(II,JJ), DZsave(N),Esave(N),
     &(uzsave(N,k),k=1,kb-1),
     &(vzsave(N,k),k=1,kb-1),
     &(Szsave(N,k),k=1,kb-1)
c      write(40+N,198) (ucenter(N,k),k=1,kb-1)
c      write(40+N,198) (uzsave(N,k),k=1,kb-1)
c      write(40+N,198) (vcenter(N,k),k=1,kb-1)
c      write(40+N,198) (vzsave(N,k),k=1,kb-1)
c      write(40+N,198) (S(II,JJ,k),k=1,kb-1)
c      write(40+N,198) (Szsave(N,k),k=1,kb-1)
c      write(40+N,198) (T(II,JJ,k),k=1,kb-1)
c      write(40+N,198) (Tzsave(N,k),k=1,kb-1)
      endif
1170  continue
C++++++ set to 1 to tell v tseries are opened
       ifopenedv=1
CMSV------------------------------------------------

      END IF
C
 190  FORMAT(8(1PE10.3))
 191  FORMAT(A10)
 192  FORMAT(8I10)
 193  FORMAT(2I5,1X,A4,I5)
C
CMSV++++
C+++++++ new formats for tampa bay individual tseries files
 194  FORMAT(i3.3) 
c 196  FORMAT(a7,2(1x,f9.4))
 196  FORMAT(a7,3(1x,f9.4))
 197  FORMAT(a7,4(1x,f9.4))
 198  FORMAT(11(f7.3))
 199  FORMAT(a7,3(1x,f9.4),11(f7.3),11(f7.3),11(f7.3))
CMSV++++
      DO 230 N=1,EPTS
 230  ESAVE(N)=0.0
      DO 235 N=1,VPTS
 235  DZSAVE(N)=0.0
      DO 240 K=1,KB
      DO 240 N=1,VPTS
      UZSAVE(N,K)=0.0
      VZSAVE(N,K)=0.0
      SZSAVE(N,K)=0.0
      CCCZSAVE(N,K)=0.0
      TZSAVE(N,K)=0.0
 240  CONTINUE
      DO 245 K=1,KB
      DO 245 N=1,FPTS
 245  CCFLUX(N,K)=0.0
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
 8250 CONTINUE
C
      RETURN
      END

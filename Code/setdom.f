      SUBROUTINE SETDOM(ADVUA,ADVVA,ADVUU,ADVVV,DRHOX,DRHOY,TRNU,TRNV)
C     VERSION(09/21/91) was the Boris Galperin domain/grid shift version
C     VERSION(04/12/99) MSV undid the domain shift
C
      INCLUDE 'comdeck'
C
C-----------------------------------------------------------------------
C     INITIAL CONDITION PREPARATION PROGRAM FOR A GENERAL CIRCULATION
C        MODEL IN ORTHOGONAL CURVILINEAR COORDINATES
C-----------------------------------------------------------------------
C
      DIMENSION ADVUA(IM,JM),ADVVA(IM,JM),ADVUU(IM,JM),ADVVV(IM,JM)
      DIMENSION DRHOX(IM,JM,KB),DRHOY(IM,JM,KB),TRNU(IM,JM),TRNV(IM,JM)
      DIMENSION IVAR(IM,JM),PRT(IM,KB)
      EQUIVALENCE (IVAR,TPS),(PRT,A)
      DIMENSION COM(80)
C 
C-------- ALL UNITS IN M.K.S. SYSTEM -----------------------------------
      ISTART=1
      IEND=1
      INT=0 
C
      WRITE(IUPRT,10)            
   10 FORMAT(/'...... MODEL STARTING UP FROM INITAL CONDITIONS .......')
C
C-------- ESTABLISH DEPTH ARRAY ----------------------------------------
      READ(IUGRD,20) (COM(I),I=1,80)
      WRITE(IUPRT,30) (COM(I),I=1,80)
c      WRITE(6,30) (COM(I),I=1,80)
   30 FORMAT(/1X,80A1/)
   20 FORMAT(80A1)
C
C-------- INITIALIZE SIGMA LEVELS --------------------------------------
      READ(IUGRD,20) (COM(I),I=1,80)
      WRITE(IUPRT,30) (COM(I),I=1,80)
c      WRITE(6,30) (COM(I),I=1,80)
C
      READ(IUGRD,40) IKB
   40 FORMAT(I5)
      WRITE(IUPRT,50) IKB
c      WRITE(6,50) IKB
   50 FORMAT(' KB = ',I5,/)
      IF(IKB.NE.KB) THEN
c      WRITE(6,60) IKB,KB
       STOP
      ENDIF
   60 FORMAT(//' NUMBER OF SIGMA LEVELS IN MODEL_GRID',I5,' (IKB)'/
     .         '           NOT EQUAL TO'/
     .         ' NUMBER OF SIGMA LEVELS IN COMDECK   ',I5,' (KB)'/
     .         ' PLEASE CORRECT THIS PROBLEM AND TRY AGAIN'//)
C
      DO 70 K=1,IKB
      READ(IUGRD,80) Z(K)
   70 CONTINUE
      WRITE(IUPRT,80) (Z(K),K=1,IKB)
   80 FORMAT(8F10.5)
C
      DO 90 K=1,KBM1
      DZ(K)=Z(K)-Z(K+1)
      DZR(K)=1.0/DZ(K)
      ZZ(K)=0.5*(Z(K)+Z(K+1))
   90 CONTINUE
      DO 100 K=1,KBM2
      DZZ(K)=ZZ(K)-ZZ(K+1)
  100 CONTINUE
      DZZ(KBM1)=0.0
      DZ(KB)=0.0
C
C-------- DEFINE THE METRICS OF THE COORDINATE TRANSFORMATION ----------
      READ(IUGRD,20) (COM(I),I=1,80)
      WRITE(IUPRT,30) (COM(I),I=1,80)
c      WRITE(6,30) (COM(I),I=1,80)
C
C     READ(IUGRD,110)   IIX, IJY
CMSV revised format for reading in IIX,IJY
      READ(IUGRD,111)   IIX, IJY
c      WRITE(*,110) IIX, IJY
      WRITE(IUPRT,120) IIX, IJY
  110 FORMAT(2I5,5F10.2)
CMSV revised format (added 111) read for model_grid file extra information (x,y,xco,yco)
  111 format(2i5,6f10.2,4f15.6)
  120 FORMAT(' IIX = ',I5,/' IJY = ',I5)
C
CMSV 4-12-99 PER ORIGINAL CODE+++++++++++++++++++++++++++++++++++++++++++++
      IF(IIX.NE.IM) THEN
c       WRITE(6,8) IIX,IM
       STOP
      ENDIF
    8 FORMAT (//'     MODEL_GRID I-INDEX',I5,' (IIX)',/
     .          '        DOES NOT EQUAL'/
     .          '     COMDECK    I-INDEX',I5,' (IM)'/
     .          ' PLEASE CORRECT THIS PROBLEM AND TRY AGAIN'//)
      IF(IJY.NE.JM) THEN
c       WRITE(6,9) IJY,JM
       STOP
      ENDIF
    9 FORMAT (//'     MODEL_GRID J-INDEX',I5,' (IJY)',/
     .          '        DOES NOT EQUAL'/
     .          '     COMDECK    J-INDEX',I5,' (JM)'/
     .          ' PLEASE CORRECT THIS PROBLEM AND TRY AGAIN'//)
CMSV END 4-12-99-----------------------------------------------------------
C..........SET PARAMETERS TO INITIAL VALUES..........
      DO 131 J=1,JM
      DO 131 I=1,IM
      H1(I,J)   = 0.0
      H2(I,J)   = 0.0
      H(I,J)    = 0.0
      ANG(I,J)  = 0.0
      COR(I,J) = 0.0
      dummy=0.0
  131 CONTINUE
C
      NUMB=IIX*IJY
C     print *,numb
      DO 140 N=1,NUMB

C 140 READ(IUGRD,110,END=150)
c      write(6,*) I,J,H1(I,J),H2(I,J),H(I,J)
      READ(IUGRD,111,END=150)
     .     I,J,H1(I,J),H2(I,J),H(I,J),ANG(I,J),COR(I,J)
     .,dummy,dummy,dummy,dummy,dummy
c      if ((H(I,J).gt.0).and.(j.le.5)) then
c         H(I,J) = 2.5
c      endif
c      write(6,*) i,j
c      ij = i
c      i=j
c      j=ij
 140  continue
  150 CONTINUE

c----------------------------------------------------
c     cut passes through causeways
c----------------------------------------------------
      call cutpasses()

c      write(6,*) 'H 54 27 ',h(54,27)
c  check for particle trapping points   May 2014 S.D. Meyers
c      do 11 i=1,im
c      do 11 j=1,jm
cc     western trap
c         if ((fsm(i,j).eq.1).and.(fsm(i-1,j).eq.0).and. 
c     &        (fsm(i,j+1).eq.0).and.(fsm(i,j-1).eq.0)) utrap(i,j)=1
cc     eastern trap
c         if ((fsm(i,j).eq.1).and.(fsm(i+1,j).eq.0).and. 
c     &        (fsm(i,j+1).eq.0).and.(fsm(i,j-1).eq.0)) utrap(i,j)=2
cc     southern trap
c         if ((fsm(i,j).eq.1).and.(fsm(i,j-1).eq.0).and. 
c     &        (fsm(i+1,j).eq.0).and.(fsm(i-1,j).eq.0)) vtrap(i,j)=1
cc     northern trap
c         if ((fsm(i,j).eq.1).and.(fsm(i,j+1).eq.0).and. 
c     &        (fsm(i+1,j).eq.0).and.(fsm(i-1,j).eq.0)) vtrap(i,j)=2
c 11   continue


CMSV++++++ now apply the depth scaling 8-16-97
c      DO 221 I=1,IM
c      DO 221 J=1,JM
c      H(I,J)=H(I,J)*depthlb
c      if (ang(i,j).gt.180) ang(i,j) = ang(i,j)-360  ! attempt to fix transport issue at cell boundar
c  221 CONTINUE


C
C-------- DEFINE MASK FOR FREE SURFACE HEIGHT = FSM --------------------
C-------- DEFINE MASK FOR (U,V) VELOCITY = (DUM,DVM) -------------------
C-------- ADD 1 FOOT (APPROX. 0.35 M) TO ACCOUNT FOR MEAN LOW WATER ----
C++++++ MSV revised for egmont mlw to msl using 1.13'/3.2808'/m=.34443

c      write(6,*) 'Using flat bottom option'
      DO 230 J=1,JM
      DO 230 I=1,IM
      IF((ibath.lt.3).and.(H(I,J).GT.0.0)) H(I,J)=H(I,J)+0.34443 ! MLLW-> MSL
      IF ((ibath.eq.3).and.(H(I,J).GT.0.0)) then    ! NOAA grid but still walled around coastline
c         if (h(i,j).lt.hdry) h(i,j) = hdry
         if (h(i,j).lt.0.644) h(i,j) = hmin
         H(I,J)=H(I,J)+doffset 
      ENDIF
c      IF(H(I,J).GT.0.0) H(I,J)=3.  ! flat bottom option
  230 CONTINUE


      DO 240 J=1,JM
      DO 240 I=1,IM
      FSM(I,J)=1.0
c      FSM0(I,J)=1.0
      WD(I,J) = FSM(I,J)
      DUM(I,J)=1.0
      DVM(I,J)=1.0
      IF(H(I,J).LE.0.0) THEN
      H(I,J)=0.0
      FSM(I,J)=0.0
c      FSM0(I,J)=0.0
      DUM(I,J)=0.0
      DVM(I,J)=0.0
      ENDIF
c      if (j.lt.80) write(6,*) 'i j fsm',i,j,fsm(i,j)
  240 H(I,J)=H(I,J)+1.0E-03
C
      DO 250 J=1,JMM1
      DO 250 I=1,IM
      IF(FSM(I,J).EQ.0.0.AND.FSM(I,J+1).NE.0.0) DVM(I,J+1)=0.0
  250 CONTINUE
      DO 260 J=1,JM
      DO 260 I=1,IMM1
      IF(FSM(I,J).EQ.0.0.AND.FSM(I+1,J).NE.0.0) DUM(I+1,J)=0.0
  260 CONTINUE
C
      DO 270 J=1,JM
      DO 270 I=1,IM
      H1(I,J)=H1(I,J)+1.E-10
      H2(I,J)=H2(I,J)+1.E-10
      ANG(I,J)=ANG(I,J)*2.*3.14159265/360.
      COR(I,J)=2.*7.292E-5*SIN(COR(I,J)*2.*3.14159/360.)*FSM(I,J)
  270 CONTINUE
C
      DO 280 J=1,JM
      DO 280 I=1,IM
      D(I,J)=H(I,J)
  280 DT(I,J)=H(I,J)
C
      DO 290 J=2,JMM1
      DO 290 I=2,IMM1
      ART(I,J)=H1(I,J)*H2(I,J)
      ARU(I,J)=.25E0*(H1(I,J)+H1(I-1,J))*(H2(I,J)+H2(I-1,J))
      ARV(I,J)=.25E0*(H1(I,J)+H1(I,J-1))*(H2(I,J)+H2(I,J-1))
  290 CONTINUE
      DO 300 J=1,JM
      ARU(IM,J)=ARU(IMM1,J)
  300 ARU(1,J)=ARU(2,J)
      DO 320 I=1,IM
      ARV(I,JM)=ARV(I,JMM1)
  320 ARV(I,1)=ARV(I,2)
C
      DO 330 J=1,JM
      DO 330 I=1,IM
  330 TPS(I,J)=SQRT(H1(I,J)**2+H2(I,J)**2)
C 
      CALL MAXMIN(TPS,FSM,IM,JM,FMAX,FMIN)
C 
      IF(TOR.NE.'BAROTROPIC') THEN
      DO 350 K=1,KBM1
      DO 350 J=1,JM
      DO 350 I=1,IM
      AAM(I,J,K)=HORCON*TPS(I,J)/FMIN*FSM(I,J)
  350 CONTINUE
      DO 375 J=1,JM
      DO 375 I=1,IM
 375  AAM2D(I,J)=0.0
      DO 380 K=1,KBM1
      DO 380 J=1,JM
      DO 380 I=1,IM
      AAM2D(I,J)=AAM2D(I,J)+AAM(I,J,K)*DZ(K)
 380  CONTINUE
C
      ELSE
C
      DO 355 J=1,JM
      DO 355 I=1,IM
      AAM2D(I,J)=HORCON*TPS(I,J)/FMIN*FSM(I,J)
  355 CONTINUE
C 
      ENDIF          
C
      DO 360 J=1,JM
      DO 360 I=1,IM
      KM(I,J,1)=0.0
      KM(I,J,KB)=0.0
      KH(I,J,1)=0.0
      KH(I,J,KB)=0.0
      KQ(I,J,1)=0.0
      KQ(I,J,KB)=0.0
      L (I,J,1)=0.0
      L (I,J,KB)=0.0
      Q2(I,J,1)=0.0
      Q2(I,J,KB)=0.0
      Q2B(I,J,1)=0.0
      Q2B(I,J,KB)=0.0
      Q2L(I,J,1)=0.0
      Q2L(I,J,KB)=0.0
      Q2LB(I,J,1)=0.0
      Q2LB(I,J,KB)=0.0
  360 CONTINUE
C
      DO 370 K=2,KBM1
      DO 370 J=1,JM
      DO 370 I=1,IM
      Q2B(I,J,K)=1.E-5*FSM(I,J)
      Q2(I,J,K)=Q2B(I,J,K)
      Q2LB(I,J,K)=Q2B(I,J,K)
      Q2L(I,J,K)=Q2B(I,J,K)
      L(I,J,K)=1.0*FSM(I,J)
      KM(I,J,K)=UMOL*FSM(I,J)
      KQ(I,J,K)=UMOL*FSM(I,J)
      KH(I,J,K)=UMOL/VPRNU*FSM(I,J)
  370 CONTINUE
C
      IF(VERTMIX.EQ.'CONSTANT  ') UMOL=0.0
C
c      write(6,*) 'Hb 54 27 ',h(54,27)

      RETURN
      END 
     

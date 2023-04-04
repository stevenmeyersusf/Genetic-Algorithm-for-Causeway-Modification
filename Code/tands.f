      SUBROUTINE TANDS(infile)
C     VERSION(03/02/90)
C     VERSION(04/12/99) MSV removed domain specific hardwiring of intitial fields
C
Cdcb added CCC to denote a passive concentration 9/10/99
Cdcb added read for tands_new.dat from epc for top runs
      INCLUDE 'comdeck'
C
      DIMENSION COM(80)
      DIMENSION ZM(KBM1),TI(KB),SI(KB),CCCI(KB)
      DIMENSION TA(KSL),SA(KSL),TS(IM,JM,KSL),SS(IM,JM,KSL)
      DIMENSION CCCA(KSL),CCCS(IM,JM,KSL)

      integer infile  ! flag to use input file
C
      IF(KSL.GT.3*KB) THEN
      WRITE(IUPRT,10)
 10   FORMAT(' PROBLEM:  KSL > 3*KB  FIX EQUIVALENCE IN tands.f ')
      STOP
      END IF
C
      DO 120 K=1,KB
      DO 120 J=1,JM
      DO 120 I=1,IM
      T(I,J,K)=25.0*FSM(I,J)
c      TS(I,J,K) = 13.*FSM(I,J)  ! jan 2001 initialization
c      S(I,J,K)=34.  !2006 cold !31.7-j*0.1  ! 31.7 for 2004 cold start
      S(I,J,K)=10.  !2006 cold !31.7-j*0.1  ! 31.7 for 2004 cold start
c      TS(I,J,K)=26.0*FSM(I,J)
c      S(I,J,K)=31.0-j*0.1
c      SS(I,J,K)=31.0-j*0.1
c      CCCS(I,J,K)=0.0*FSM(I,J)
CMSV INITIALIZE SALINITY AS NULL
C     TS(I,J,K)=00.0*FSM(I,J)
C     SS(I,J,K)=00.0*FSM(I,J)
C     CCCS(I,J,K)=00.0*FSM(I,J)
  120 CONTINUE

C
      DO 128 K=1,KSL
      DO 128 J=1,JM
      DO 128 I=1,IM
      TS(I,J,K)=25.0*FSM(I,J)
      SS(I,J,K)=10.0*FSM(I,J) !-j*1.0/4.0  ! 35 uniform for test
      CCCS(I,J,K)=0.0*FSM(I,J)
 128  continue

c      write(6,*) 'initial boundary ',t(17,3,1),s(17,3,1)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (infile.eq.1) then
C++++++ READ  A 3-D TANDS.DAT FILE
C234567
      iutands=18
c        write(6,*) 'INITIALIZING S with sal_jan2001.dat'
 11    FORMAT(80A1)
        open (iutands,form='formatted',file=
     &     'sal_jan2001.dat',status='old')
        do m=1,5  ! skip header
           READ(IUTANDS,11)  (COM(I),I=1,80)  
        enddo
  121   continue
        read(iutands,103,end=122) i,j,k,S(i,j,k)
        ss(i,j,k) = s(i,j,k)
c        write(6,*) i,j,k
        goto 121
  122   continue
        close(iutands)

      CALL DENS
C
      DO 279 K=1,KBM1
      DO 279 J=1,JM
      DO 279 I=1,IM
      TMEAN(I,J,K)=TS(I,J,K)*FSM(I,J)
      SMEAN(I,J,K)=SS(I,J,K)*FSM(I,J)
      RMEAN(I,J,K)=RHO(I,J,K)*FSM(I,J)
  279 CONTINUE

      goto 123
      endif
C---------------------------------------------------------------------------



C------
C------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++ for this option TS = 25. set above
C+++++++                 SS  @sl1 is read in then sl2,3 are set equal to
c        iutands=18
  103   format (3i8,f8.2)
C--------------------------------------------------------------------------
C
c        open (42,file='initial.tands',form='formatted')
c        do 121 n=1,(im*jm)
c        read(42,103,end=121) i,j,(TS(i,j,k),k=1,ksl),
c     &(SS(i,j,k),k=1,ksl),(CCCS(i,j,k),k=1,ksl)
c  121   continue
c  103   format (2i5,40f5.0)
C

c        write(6,*) 'get model average depth-profiles'
      DO 220 K=1,KSL
      TA(K)=0.0
      SA(K)=0.0
      CCCA(K)=0.0
      COUNT=0.0
      DO 210 J=1,JM
      DO 210 I=1,IM
      IF (-H(I,J).GT.DPTHSL(K)) GO TO 210
      COUNT=COUNT+FSM(I,J)
      TA(K)=TA(K)+TS(I,J,K)*FSM(I,J)
      SA(K)=SA(K)+SS(I,J,K)*FSM(I,J)
      CCCA(K)=CCCA(K)+CCCS(I,J,K)*FSM(I,J)
  210 CONTINUE
      IF (COUNT.LT.1.0) GO TO 220
      TA(K)=TA(K)/COUNT
      SA(K)=SA(K)/COUNT
      CCCA(K)=CCCA(K)/COUNT
  220 CONTINUE
      DO 250 J=1,JM
      DO 250 I=1,IM
      IF (FSM(I,J).EQ.0.0) GO TO 250
      DO 230 K=1,KBM1
      ZM(K)=ZZ(K)*H(I,J)
  230 CONTINUE
      CALL SINTER (DPTHSL,TA,ZM,TI,KSL,KBM1)
      CALL SINTER (DPTHSL,SA,ZM,SI,KSL,KBM1)
      CALL SINTER (DPTHSL,CCCA,ZM,CCCI,KSL,KBM1)
      DO 240 K=1,KBM1
      T(I,J,K)=TI(K)
      S(I,J,K)=SI(K)
      CCC(I,J,K)=CCCI(K)
  240 CONTINUE
  250 CONTINUE
C
      CALL DENS
C
      DO 280 K=1,KBM1
      DO 280 J=1,JM
      DO 280 I=1,IM
      TMEAN(I,J,K)=T(I,J,K)*FSM(I,J)
      SMEAN(I,J,K)=S(I,J,K)*FSM(I,J)
      CCCMEAN(I,J,K)=CCC(I,J,K)*FSM(I,J)
      RMEAN(I,J,K)=RHO(I,J,K)*FSM(I,J)
  280 CONTINUE
C
      DO 350 J=1,JM
      DO 350 I=1,IM
      IF (FSM(I,J).EQ.0.0) GO TO 350
      DO 320 K=1,KSL
      TA(K)=TS(I,J,K)
      SA(K)=SS(I,J,K)
      CCCA(K)=CCCS(I,J,K)
  320 CONTINUE
      DO 330 K=1,KBM1
      ZM(K)=ZZ(K)*H(I,J)
  330 CONTINUE
      CALL SINTER (DPTHSL,TA,ZM,TI,KSL,KBM1)
      CALL SINTER (DPTHSL,SA,ZM,SI,KSL,KBM1)
      CALL SINTER (DPTHSL,CCCA,ZM,CCCI,KSL,KBM1)
      DO 340 K=1,KBM1
      T(I,J,K)=TI(K)*FSM(I,J)
      S(I,J,K)=SI(K)*FSM(I,J)
      CCC(I,J,K)=CCCI(K)*FSM(I,J)
  340 CONTINUE
  350 CONTINUE


C****** jump to here if read in a full 3-d t and s
 123  continue
C******


      CALL DENS
C
      DO 430 K=1,KBM1
      DO 430 J=1,JM
      DO 430 I=1,IM
C++MSV
      T(I,J,K)=T(I,J,K)*FSM(I,J)
      S(I,J,K)=S(I,J,K)*FSM(I,J)

C--MSV
      TB(I,J,K)=T(I,J,K)
      SB(I,J,K)=S(I,J,K)
      CCCB(I,J,K)=CCC(I,J,K)
  430 CONTINUE
C
      DO 440 K=1,KB
      DO 440 J=1,JM
      DO 440 I=1,IM
      A(I,J,K)=0.0
      C(I,J,K)=0.0
      VH(I,J,K)=0.0
      VHP(I,J,K)=0.0
      PROD(I,J,K)=0.0
      DTEF(I,J,K)=0.0
  440 CONTINUE
C

c      write(6,*) 'final skyway ',t(25,35,1),s(25,35,1)

      RETURN
      END
      


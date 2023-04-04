c23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE BCDATA2
C     VERSION(09/21/91)
      INCLUDE 'comdeck'
C
c     altered by S.D. Meyers 6/2005 for conventional OBC
c     use with bcond2.f

      DIMENSION COM(80)
      DIMENSION ZM(KBM1),T1(KSL),T2(KBM1),S1(KSL),S2(KBM1)
      DIMENSION CCC1(KSL),CCC2(KBM1)
C
c      write(6,*) 'in bcdata2 ',iwds

C-------- READ IN STANDARD LEVELS HERE ---------------------------------
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
c      WRITE(6,12) (COM(I),I=1,80)
 11   FORMAT(80A1)
 12   FORMAT(/1X,80A1/)
C
      READ(IURUN,24) IKSL
      WRITE(IUPRT,41) IKSL
c      WRITE(6,41) IKSL
   41 FORMAT(' KSL = ',I5,/)
      IF(IKSL.NE.KSL) THEN
c       WRITE(6,42) IKSL, KSL
       STOP
      ENDIF
   42 FORMAT(//' NUMBER OF STANDARD LEVELS IN RUN_DATA',I5,' (IKSL)'/
     .         '           DO NOT EQUAL'/
     .         ' NUMBER OF STANDARD LEVELS IN COMDECK ',I5,' (KSL)'/
     .         ' PLEASE CORRECT THIS PROBLEM AND TRY AGAIN'//)
C
      READ(IURUN,77) (DPTHSL(K),K=1,KSL)
      WRITE(IUPRT,77) (DPTHSL(K),K=1,KSL)
c      write(6,*) 'DEPTHS'
c      WRITE(6,77) (DPTHSL(K),K=1,KSL)
C
C-------- READ IN BOUNDARY CONDITIONS HERE -----------------------------
C-------- ELEVATION BOUNDARY -------------------------------------------
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
c      write(6,*) 'elevation descriptor line'
c      WRITE(6,12) (COM(I),I=1,80)
      READ(IURUN,24) NUMEBC,OPTEBC
      WRITE(IUPRT,24) NUMEBC,OPTEBC
 24   FORMAT(I5,1X,A20)
      IF(NUMEBC.EQ.0) GO TO 126
      IF(OPTEBC(1:5).NE.'DATA ') THEN
      DO 110 N=1,NUMEBC
c      write(6,*) N,NUMEBC
      READ(IURUN,9) IETA(N),JETA(N),ICON(N),JCON(N),EMEAN(N)
      WRITE(IUPRT,9) IETA(N),JETA(N),ICON(N),JCON(N),EMEAN(N)
c      WRITE(6,9) IETA(N),JETA(N),ICON(N),JCON(N),EMEAN(N)
      READ(IURUN,7) (AMP(N,I),I=1,6)
      WRITE(IUPRT,7) (AMP(N,I),I=1,6)
c      WRITE(6,7) (AMP(N,I),I=1,6)
      READ(IURUN,7) (PHASE(N,I),I=1,6)   ! minutes
      WRITE(IUPRT,7) (PHASE(N,I),I=1,6)
c      WRITE(6,7) (PHASE(N,I),I=1,6)
      DO 105 I=1,6
 105  PHASE(N,I)=60.*PHASE(N,I)          ! sec
 110  CONTINUE
 8    FORMAT(2I5,4F10.5)
 9    FORMAT(4I5,1F10.5)
 7    FORMAT(6F10.5)
      ELSE
      READ(IURUN,76) (IETA(N),JETA(N),ICON(N),JCON(N),N=1,NUMEBC)
      WRITE(IUPRT,76) (IETA(N),JETA(N),ICON(N),JCON(N),N=1,NUMEBC)
76    FORMAT(16I5)
C++++ MSV changed 7-30-98100000 to 200000 to allow 5 year runs at tides on 15 minute inc
c      DO 122 I=1,100000
      DO 122 I=1,200000
      READ(IURUN,77,ERR=125) TIME
c      write(6,*) 'el t=',time
      WRITE(IUPRT,77) TIME     
      READ(IURUN,77) (EBDRY(N),N=1,NUMEBC)
c       WRITE(6,77) (EBDRY(N),N=1,NUMEBC)           

c      do n=1,numebc
c         ebdry(n) = ebdry(n)+1  ! add in for sea level rise elevation +1m
c      enddo

      WRITE(IUPRT,77) (EBDRY(N),N=1,NUMEBC)           
      WRITE(IUT90,78) TIME
 122  WRITE(IUT90,78) (EBDRY(N),N=1,NUMEBC)
 125  BACKSPACE IURUN

      END IF

c    fix velocity masks (error in setdom). Use with bcond2.f
c     adapted from setdom.f L150
c     5/31/2005
      do n=1,numebc
         ie = ieta(n)
         je = jeta(n)
         ic = icon(n)
         jc = jcon(n)
cc correct southern OBC
         IF(DVM(IE,JE).EQ.0.0.AND.ie.eq.ic) THEN
            DVM(IE,JE)=1.0
            DUM(IE,JE)=0.0
         endif
cc correct western OBC
         IF(DUM(IE,JE).EQ.0.0.AND.je.eq.jc) then
            DVM(IE,JE)=0.0
            DUM(IE,JE)=1.0
         endif
      enddo





c        write(6,*) 'reading boundary TS
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
c      write(6,*) 'temp/sal header '
c      WRITE(6,12) (COM(I),I=1,80)

      DO 435 I=1,200000
      READ(IURUN,77,ERR=127) TIME
c      if (i.lt.10)WRITE(6,*) 'ts time ',TIME
      WRITE(IUPRT,77) TIME     
      WRITE(IUT94,78) TIME
      DO 239 N=1,NUMEBC
c        write(6,*) i,n,time,numebc
        READ(IURUN,79) (TBDRYSL(N,K),K=1,KSL)
        READ(IURUN,79) (SBDRYSL(N,K),K=1,KSL)
        READ(IURUN,79) (CCCBDRYSL(N,K),K=1,KSL)  ! cut for nowcast restart
c      if (i.lt.10)WRITE(6,*) 'ts time ',TIME,(SBDRYSL(N,K),K=1,KSL)

c     wind only or
c     set T=25 everywhere to eliminate T gradients
        do 111 itz=1,ksl
           cccbdrysl(n,itz)=0.
           tbdrysl(n,itz)=25.
c           sbdrysl(n,itz)=35.
 111    enddo

        WRITE(IUPRT,79) (TBDRYSL(N,K),K=1,KSL)       
        WRITE(IUPRT,79) (SBDRYSL(N,K),K=1,KSL)
        WRITE(IUPRT,79) (CCCBDRYSL(N,K),K=1,KSL)     
c        WRITE(6,79) (TBDRYSL(N,K),K=1,KSL)       
c        WRITE(6,79) (SBDRYSL(N,K),K=1,KSL)
c        WRITE(6,79) (CCCBDRYSL(N,K),K=1,KSL)     
        
        II=IETA(N)
        JJ=JETA(N)
        DO 243 K=1,KBM1
243       ZM(K)=ZZ(K)*H(II,JJ)
        DO  241 K=1,KSL
          T1(K)=TBDRYSL(N,K)
          CCC1(K)=CCCBDRYSL(N,K)
241       S1(K)=SBDRYSL(N,K)
      CALL SINTER(DPTHSL,T1,ZM,T2,KSL,KBM1)
      CALL SINTER(DPTHSL,S1,ZM,S2,KSL,KBM1)
      CALL SINTER(DPTHSL,CCC1,ZM,CCC2,KSL,KBM1)
        DO 246 K=1,KBM1
          TBDRY(N,K)=T2(K)
          CCCBDRY(N,K)=CCC2(K)
246       SBDRY(N,K)=S2(K)
        WRITE(IUT94,78) (TBDRY(N,K),K=1,KBM1)
        WRITE(IUT94,78) (SBDRY(N,K),K=1,KBM1)
        WRITE(IUT94,78) (CCCBDRY(N,K),K=1,KBM1)
        
239   CONTINUE
435   CONTINUE
127   BACKSPACE IURUN
   79 FORMAT(16F5.1)
 126  CONTINUE
C
c      write(6,*) 'He 54 27 ',h(54,27)
C-------- RIVER/DAM AND ONSHORE INTAKE/OUTFALL DISCHARGE BOUNDARY ------
      READ(IURUN,11)  (COM(I),I=1,80)
c      write(6,11) (com(i),i=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,8) NUMQBC
      WRITE(IUPRT,8) NUMQBC
c      WRITE(6,*) 'numqbc=',NUMQBC
      IF(NUMQBC.EQ.0) GO TO 136
      DO 321 N=1,NUMQBC
c      write(6,*) n
      READ(IURUN,5) IQD(N),JQD(N),IQC(N),JQC(N),(VQDIST(N,K),K=1,KBM1)
c      write(6,*) iqd(n),jqd(n),iqc(n),jqc(n)
321   WRITE(IUPRT,6) IQD(N),JQD(N),IQC(N),JQC(N),(VQDIST(N,K),K=1,KBM1)
5     FORMAT(4I5,20F3.0)
6     FORMAT(4I5,/,20F5.0)
      DO 130 I=1,200000
      READ (IURUN,77,ERR=135) TIME
      WRITE(IUPRT,77) TIME
c      WRITE(6,*) 'riv time=',TIME
      READ (IURUN,77) (QDIS(N),N=1,NUMQBC)
 
 3022 continue
c ................................................................      
      READ (IURUN,77) (TDIS(N),N=1,NUMQBC)
      READ (IURUN,77) (SDIS(N),N=1,NUMQBC)
      READ (IURUN,77) (CCCDIS(N),N=1,NUMQBC) ! cut for nowcast restart
      WRITE(IUPRT,77) (QDIS(N),N=1,NUMQBC)
      WRITE(IUPRT,77) (TDIS(N),N=1,NUMQBC)
      WRITE(IUPRT,77) (SDIS(N),N=1,NUMQBC)
      WRITE(IUPRT,77) (CCCDIS(N),N=1,NUMQBC)
      WRITE(IUT91,78) TIME
      WRITE(IUT91,78) (QDIS(N),N=1,NUMQBC)
      WRITE(IUT91,78) (TDIS(N),N=1,NUMQBC)
      WRITE(IUT91,78) (SDIS(N),N=1,NUMQBC)
 130  WRITE(IUT91,78) (CCCDIS(N),N=1,NUMQBC)
 135  BACKSPACE IURUN
1005  format(30f10.4)         
  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DO 140 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      IF(JD.EQ.JC) THEN
            H2(IC,JC)=H2(ID,JD)
            IF(IC.GT.ID) THEN
              DUM(IC,JC)=1.0
            ELSE
              DUM(ID,JD)=1.0
            ENDIF
      ELSE
            H1(IC,JC)=H1(ID,JD)
            IF(JC.GT.JD) THEN
              DVM(IC,JC)=1.0
            ELSE
              DVM(ID,JD)=1.0
            ENDIF
      ENDIF
            H(IC,JC)=H (ID,JD)
            FSM(IC,JC)=FSM(ID,JD)
  140    CONTINUE
C
      IF(HORZMIX.EQ.'CONSTANT  ') THEN
      DO 220 N=1,NUMQBC
      IC=IQC(N)
      JC=JQC(N)
      AAM2D(IC,JC)=0.0
      DO 220 K=1,KBM1
      AAM(IC,JC,K)=0.0
  220 CONTINUE
      ENDIF
C
 136  CONTINUE
C

C---- OFFSHORE INTAKE/OUTFALL(DIFFUSER) BOUNDARY -------------------
c      write(6,*) 'runoff'
      READ(IURUN,11)  (COM(I),I=1,80)
c      WRITE(6,12) (COM(I),I=1,80)
      READ(IURUN,8) NUMDBC
c      WRITE(6,*) 'NUM RUNOFF=',NUMDBC
      IF(NUMDBC.EQ.0) GO TO 156
      DO 322 N=1,NUMDBC
      READ(IURUN,51) IDD(N),JDD(N),(VDDIST(N,K),K=1,KBM1)
322   WRITE(IUPRT,61) IDD(N),JDD(N),(VDDIST(N,K),K=1,KBM1)
51    FORMAT(2I5,20F3.0)
61    FORMAT(2I5,/,20F5.0)
      DO 151 I=1,200000
      READ(IURUN,77,ERR=155) TIME
c      if (i.lt.10)WRITE(6,*) 'runoff t=',TIME,numdbc
      WRITE(IUPRT,77) TIME
      READ(IURUN,77) (QDIFF(N),N=1,NUMDBC)
      READ(IURUN,77) (TDIFF(N),N=1,NUMDBC)
      READ(IURUN,77) (SDIFF(N),N=1,NUMDBC)
      READ(IURUN,77) (CCCDIFF(N),N=1,NUMDBC)  ! cut for nowcast restart

c     for wind only test
c      do nn=1,numdbc  
c         qdiff(nn) = 0.
c      enddo

c      do 555 n=1,numdbc
c      if (qdiff(n).eq.25) then 
c         WRITE(6,*) 'QDIFF ERROR!!!!!*****',time
c      endif
c 555  continue
      WRITE(IUPRT,77) (QDIFF(N),N=1,NUMDBC)
      WRITE(IUPRT,77) (TDIFF(N),N=1,NUMDBC)
      WRITE(IUPRT,77) (SDIFF(N),N=1,NUMDBC)
      WRITE(IUPRT,77) (CCCDIFF(N),N=1,NUMDBC)
      WRITE(IUT92,78) TIME
      WRITE(IUT92,78) (QDIFF(N),N=1,NUMDBC)
      WRITE(IUT92,78) (TDIFF(N),N=1,NUMDBC)
      WRITE(IUT92,78) (SDIFF(N),N=1,NUMDBC)
      WRITE(IUT92,78) (CCCDIFF(N),N=1,NUMDBC)
 151  continue
 155  BACKSPACE IURUN
 156  CONTINUE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 

C-------- METEOROLOGICAL BOUNDARY CONDITION ----------------------------
C-------- PRECIPITATION (m/year) & EVAPORATION (m/year) ----------------
C-------- WIND FROM DIRECTION WIND IS BLOWING --------------------------
      READ(IURUN,11)  (COM(I),I=1,80)
c      WRITE(6,12) (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      DO 161 I=1,200000
      READ(IURUN,77,END=165) TIME
      WRITE(IUPRT,77) TIME
      READ(IURUN,77) QPREC,QEVAP,WDS,WDD,HFLUX

c      qprec = 0
c      qevap = 0
c      if ((i.eq.1).and.(iwds.eq.0))write(6,*) '====> WINDS SET=0 <====='
      wds = wds*iwds

      WRITE(IUPRT,77) QPREC,QEVAP,WDS,WDD,HFLUX
c      if (time.lt.79050) then 
c      if (i.lt.10)WRITE(6,*) 'met time ',TIME,wds,iwds
c         write(6,*) wds
c      endif
 77   FORMAT(8F10.2)  ! original
c 77   FORMAT(8F11.4)  ! 10.4 curl and 11.4 from charney
C+++++ REVISED BY MSV 7-23-98
C+++++ to read in p and e per month
C+++++ then convert to per sec same as rivers and diff
cdcb revised 10/1/00 to read in pre in cm/day and convert to m/s
C     QPREC=QPREC/(86400.*365.)
C     QEVAP=QEVAP/(86400.*365.)
c      QPREC=QPREC/(720.*60.*60.)
c      QEVAP=QEVAP/(720.*60.*60.)
      QPREC=QPREC/(100.*24.*60.*60.)
      QEVAP=QEVAP/(100.*24.*60.*60.)


C-------- WIND STRESS --------------------------------------------------
      WDD=180.+WDD
      WDD=AMOD(WDD,360.)
      VWIND=WDS*COS(6.28319*WDD/360.)
      UWIND=WDS*SIN(6.28319*WDD/360.)
      CD=1.2E-3
      IF(WDS.GE.11.) CD=(0.49+0.065*WDS)*1.E-3
      IF(WDS.GE.25.) CD=(0.49+0.065*25.)*1.E-3
      TX=1.2*CD*UWIND*WDS
      TY=1.2*CD*VWIND*WDS
      WRITE(IUT93,78) TIME
c      if ((time.gt.91400).and.(time.lt.91415)) 
c     &     WRITE(6,*) time,QPREC,QEVAP,TX,TY,HFLUX
161   WRITE(IUT93,78) QPREC,QEVAP,TX,TY,HFLUX
165   CONTINUE
 78   FORMAT(8E14.7)
C 

      RETURN
      END
      


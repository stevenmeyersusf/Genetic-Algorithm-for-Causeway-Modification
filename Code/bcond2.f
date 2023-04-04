      SUBROUTINE BCOND2(IDX,DTI2)
C     VERSION(11/22/90)
c     altered to allow flux through OBC
c     S.D. Meyers 5/2005
c     use with bcdata2.f

      INCLUDE 'comdeck'
C
      REAL YHOURS
C===============>    S2     M2     N2     K1     P1     O1    
      DATA PERIOD /43200.,44712.,45570.,86164.,86637.,92950./
      DATA PI2/6.283185307/
C
C-----------------------------------------------------------------------
C        SPECIFICATION OF OPEN BOUNDARY CONDITIONS
C   NOTE THAT AT J=2 U CALCULATION EXCLUDES ADVECTION AND DIFFUSION
C-----------------------------------------------------------------------
C        I-1           I         I+1
C
C       U(JM)=      EL(JM)      U(JM)
C
C                   V(JM)
C
C     *U(JMM1)*   *EL(JMM1)*   *U(JMM1)*
C
C                  *V(JMM1)*
C
C     *U(JMM2)*   *EL(JMM2)*   *U(JMM2)*
C                       "                                         V(IM)
C                       "
C                       "                *U(IMM1)* EL(IMM1) U(M)  EL(IM)
C                       "                              = = = = = = =
C                       "                             V(IMM1)     V(IM)
C                       "
C       *U(3)*       *EL(3)*   *U(3)*                      BC ON EL & T
C
C                    *V(3)*
C                 "
C       *U(2)*    "  *EL(2)*   *U(2)*
C                 "         "         "
C                  =  V(2)  "         "  BC ON V
C                           "         "
C          U(1)       EL(1)=    U(1) =   BC ON T
C-----------------------------------------------------------------------
C         IDX IDENTIFIES WHICH VARIABLES ARE CONSIDERED
C              1=SURFACE ELEVATION  
C              2=EXTERNAL MODE U,V 
C              3=INTERNAL MODE U,V 
C              4=TEMP,SAL,and CCC (passive tracer)
C              5=W VELOCITY
C              6= KM,KH,Q2,Q2L,L 
C              7=SURFACE FORCING AND TEMPORAL CYCLING
C-----------------------------------------------------------------------
C
      GO TO (10,20,30,40,50,60,70), IDX
C
 10   CONTINUE
C
C-------- EXTERNAL ELEVATION BOUNDARY CCONDITIONS ----------------------
      DO 100 N=1,NUMEBC  
      II=IETA(N)
      JJ=JETA(N)
      IF(OPTEBC(1:5).NE.'DATA ') THEN
      FORCE=0.0
CMSV  SET THE TIME IN HOURS FOR THE TIDE
CMSV  YHOURS=(31.+28.+31.+30.+31.0+31.0+17.0)*24.*60.*60.
      YHOURS=(0.0)*24.*60.*60.
CMSV
      DO 110 I=1,6  
 110  FORCE=AMP(N,I)*COS(PI2/PERIOD(I)*
     .(DTI*FLOAT(INT)+YHOURS-PHASE(N,I)))
     .     +FORCE
      ELF(II,JJ)=(EMEAN(N)+FORCE)*RAMP
      ELSE
      ELF(II,JJ)=EBDRY(N)*RAMP
      ENDIF
 100  CONTINUE
C
      DO 120 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      ELF(IC,JC)=ELF(ID,JD)
120   CONTINUE
C
      DO 130 J=1,JM
      DO 130 I=1,IM
 130  ELF(I,J)=ELF(I,J)*FSM(I,J)
      RETURN
C
C-------- EXTERNAL VELOCITY BOUNDARY CONDITIONS-------------------------
  20  CONTINUE
      DO 200 N=1,NUMQBC  
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)  
      FRESH=QDIS(N)/(H(IC,JC)+ELF(IC,JC))*RAMP
      IF(JD.EQ.JC) THEN
         IF(ID.LT.IC) THEN
            UAF(IC,JC)=-FRESH/H2(IC,JC)
         ELSE
            UAF(ID,JD)=+FRESH/H2(ID,JD)
         ENDIF
      ELSE
         IF(JD.LT.JC) THEN
            VAF(IC,JC)=-FRESH/H1(IC,JC)
         ELSE
            VAF(ID,JD)=+FRESH/H1(ID,JD)
         ENDIF
      ENDIF
 200  CONTINUE
C
      DO 210 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
c     5/31/2005      use with bcdata2.f
      IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN  ! eastern 
        UAF(IE,JE)=UAF(IE-1,JE)
c        VAF(IE,JE)=VAF(IE-1,JE)
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN  ! western
        UAF(IE,JE)=UAF(IE+1,JE)
c        VAF(IE,JE)=VAF(IE+1,JE)
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN  ! northern
        VAF(IE,JE)=VAF(IE,JE-1)
c        UAF(IE,JE)=UAF(IE,JE-1)
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN  ! southern
        VAF(IE,JE)=VAF(IE,JE+1)
c        UAF(IE,JE)=UAF(IE,JE+1)
      END IF
 210  CONTINUE
      DO 135 J=1,JM
      DO 135 I=1,IM
      UAF(I,J)=UAF(I,J)*DUM(I,J)
 135  VAF(I,J)=VAF(I,J)*DVM(I,J)
      RETURN
C
C-------- INTERNAL VELOCITY BOUNDARY CONDITIONS ------------------------
  30  CONTINUE
      DO 320 N=1,NUMQBC  
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      DO 320 K=1,KBM1
      FRESH=QDIS(N)/(H(IC,JC)+ELF(IC,JC))*RAMP*VQDIST(N,K)/100.
      IF(JD.EQ.JC) THEN
         IF(ID.LT.IC) THEN
            UF(IC,JC,K)=-FRESH/(H2(IC,JC)*DZ(K))
         ELSE
            UF(ID,JD,K)=+FRESH/(H2(ID,JD)*DZ(K))
         ENDIF
      ELSE
         IF(JD.LT.JC) THEN
            VF(IC,JC,K)=-FRESH/(H1(IC,JC)*DZ(K))
         ELSE
            VF(ID,JD,K)=+FRESH/(H1(ID,JD)*DZ(K))
         ENDIF
      ENDIF
  320 CONTINUE
C
      DO 321 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
c     5/31/2005  use with bcdata2
      DO 321 K=1,KBM1
      IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN  ! eastern
        UF(IE,JE,K)=UF(IE-1,JE,K)
c        VF(IE,JE,K)=VF(IE-1,JE,K)
        WVBOT(IE,JE)=WVBOT(IE-1,JE)
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN  ! western
        UF(IE,JE,K)=UF(IE+1,JE,K)
c        VF(IE,JE,K)=VF(IE+1,JE,K)
        WVBOT(IE,JE)=WVBOT(IE+1,JE)
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN  ! northern
        VF(IE,JE,K)=VF(IE,JE-1,K)
c        UF(IE,JE,K)=UF(IE,JE-1,K)
        WUBOT(IE,JE)=WUBOT(IE,JE-1)
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN  ! southern
        VF(IE,JE,K)=VF(IE,JE+1,K)
c        UF(IE,JE,K)=UF(IE,JE+1,K)
        WUBOT(IE,JE)=WUBOT(IE,JE+1)
      END IF
  321 CONTINUE
C 
      DO 360 K=1,KBM1
      DO 360 I=1,IM
      DO 360 J=1,JM
      UF(I,J,K)=UF(I,J,K)*DUM(I,J)
      VF(I,J,K)=VF(I,J,K)*DVM(I,J)
 360  CONTINUE
      RETURN
C
C-------- TEMPERATURE & SALINITY BOUNDARY CONDITIONS -------------------
  40  CONTINUE
      DO 235 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      DO 235 K=1,KBM1
      IF (VQDIST(N,K).NE.0.0 .AND. QDIS(N).GT.0.0) THEN
      UF(IC,JC,K)=TDIS(N)
      VF(IC,JC,K)=SDIS(N)

      CCCF(IC,JC,K)=CCCDIS(N)

      ELSE
      UF(IC,JC,K)=UF(ID,JD,K)
      VF(IC,JC,K)=VF(ID,JD,K)

      CCCF(IC,JC,K)=CCCF(ID,JD,K)

      END IF
 235  CONTINUE
      DO 236 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
C
      DO 170 K=1,KBM1
C

      TBDY=TBDRY(N,K)
      SBDY=SBDRY(N,K)

      CCCBDY=CCCBDRY(N,K)

C
      IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
C
C-------- EAST SIDE ----------------------------------------------------
      VEL=U(IE,JE,K)
      IF (VEL.LE.0.0) THEN
      CPH=0.0
      TLAGI=TLAG
      ELSE
      TLAGI=0.0 
      CPH=VEL*DTI*2.0/(H1(IE,JE)+H1(IE-1,JE))
      END IF
      UF(IE,JE,K)=TB(IE,JE,K)+CPH*(2.0*T(IE-1,JE,K)-TB(IE,JE,K)) 
     1 +(TBDY-TB(IE,JE,K))*TLAGI*DTI2 
      UF(IE,JE,K)=UF(IE,JE,K)/(1.0+CPH)
      VF(IE,JE,K)=SB(IE,JE,K)+CPH*(2.0*S(IE-1,JE,K)-SB(IE,JE,K)) 
     1 +(SBDY-SB(IE,JE,K))*TLAGI*DTI2 
      VF(IE,JE,K)=VF(IE,JE,K)/(1.0+CPH)

      CCCF(IE,JE,K)=CCCB(IE,JE,K)+CPH*(2.0*CCC(IE-1,JE,K)-CCCB(IE,JE,K)) 
     1 +(CCCBDY-CCCB(IE,JE,K))*TLAGI*DTI2 
      CCCF(IE,JE,K)=CCCF(IE,JE,K)/(1.0+CPH)

      GO TO 170
C
      ELSE IF (FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
C
C-------- WEST SIDE ----------------------------------------------------
      VEL=U(IE+1,JE,K)
      IF (VEL.GE.0.0) THEN
      CPH=0.0
      TLAGI=TLAG
      ELSE
      CPH=VEL*DTI*2.0/(H1(IE,JE)+H1(IE+1,JE))
      TLAGI=0.0 
      END IF
      UF(IE,JE,K)=TB(IE,JE,K)-CPH*(2.0*T(IE+1,JE,K)-TB(IE,JE,K)) 
     1 +(TBDY-TB(IE,JE,K))*TLAGI*DTI2 
      UF(IE,JE,K)=UF(IE,JE,K)/(1.0-CPH)
      VF(IE,JE,K)=SB(IE,JE,K)-CPH*(2.0*S(IE+1,JE,K)-SB(IE,JE,K)) 
     1 +(SBDY-SB(IE,JE,K))*TLAGI*DTI2 
      VF(IE,JE,K)=VF(IE,JE,K)/(1.0-CPH)

      CCCF(IE,JE,K)=CCCB(IE,JE,K)-CPH*(2.0*CCC(IE+1,JE,K)-CCCB(IE,JE,K)) 
     1 +(CCCBDY-CCCB(IE,JE,K))*TLAGI*DTI2 
      CCCF(IE,JE,K)=CCCF(IE,JE,K)/(1.0-CPH)

      GO TO 170
      ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
C
C-------- NORTH SIDE ---------------------------------------------------
      VEL=V(IE,JE,K)
      IF (VEL.LE.0.0) THEN
      CPH=0.0
      TLAGI=TLAG
      ELSE
      CPH=VEL*DTI*2.0/(H2(IE,JE)+H2(IE,JE-1))
      TLAGI=0.0 
      END IF
      UF(IE,JE,K)=TB(IE,JE,K)+CPH*(2.0*T(IE,JE-1,K)-TB(IE,JE,K)) 
     1 +(TBDY-TB(IE,JE,K))*TLAGI*DTI2 
      UF(IE,JE,K)=UF(IE,JE,K)/(1.0+CPH)
      VF(IE,JE,K)=SB(IE,JE,K)+CPH*(2.0*S(IE,JE-1,K)-SB(IE,JE,K)) 
     1 +(SBDY-SB(IE,JE,K))*TLAGI*DTI2 
      VF(IE,JE,K)=VF(IE,JE,K)/(1.0+CPH)

      CCCF(IE,JE,K)=CCCB(IE,JE,K)+CPH*(2.0*CCC(IE,JE-1,K)-CCCB(IE,JE,K)) 
     1 +(CCCBDY-CCCB(IE,JE,K))*TLAGI*DTI2 
      CCCF(IE,JE,K)=CCCF(IE,JE,K)/(1.0+CPH)

      GO TO 170
      ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
C
C-------- SOUTH SIDE ---------------------------------------------------
      VEL=V(IE,JE+1,K)
      IF (VEL.GE.0.0) THEN
      CPH=0.0
      TLAGI=TLAG
      ELSE
      CPH=VEL*DTI*2.0/(H2(IE,JE)+H2(IE,JE+1))
      TLAGI=0.0 
      END IF
      UF(IE,JE,K)=TB(IE,JE,K)-CPH*(2.0*T(IE,JE+1,K)-TB(IE,JE,K)) 
     1 +(TBDY-TB(IE,JE,K))*TLAGI*DTI2 
      UF(IE,JE,K)=UF(IE,JE,K)/(1.0-CPH)
      VF(IE,JE,K)=SB(IE,JE,K)-CPH*(2.0*S(IE,JE+1,K)-SB(IE,JE,K)) 
     1 +(SBDY-SB(IE,JE,K))*TLAGI*DTI2 
      VF(IE,JE,K)=VF(IE,JE,K)/(1.0-CPH)

      CCCF(IE,JE,K)=CCCB(IE,JE,K)-CPH*(2.0*CCC(IE,JE+1,K)-CCCB(IE,JE,K)) 
     1 +(CCCBDY-CCCB(IE,JE,K))*TLAGI*DTI2 
      CCCF(IE,JE,K)=CCCF(IE,JE,K)/(1.0-CPH)

C
C-------- DONE ---------------------------------------------------------
      END IF
  170 CONTINUE
  236 CONTINUE
      DO 240 K=1,KBM1
      DO 240 I=1,IM
      DO 240 J=1,JM
      UF(I,J,K)=UF(I,J,K)*FSM(I,J)
      VF(I,J,K)=VF(I,J,K)*FSM(I,J)

      CCCF(I,J,K)=CCCF(I,J,K)*FSM(I,J)

 240  CONTINUE
      RETURN
C
C-------- VERTICAL VELOCITY BOUNDARY CONDITIONS ------------------------
  50  CONTINUE
      DO 245 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      DO 245 K=1,KBM1
      W(IC,JC,K)=W(ID,JD,K)
  245 CONTINUE
      DO 246 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
      DO 246 K=1,KBM1
      IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        W(IE,JE,K)=W(IE-1,JE,K)
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        W(IE,JE,K)=W(IE+1,JE,K)
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
        W(IE,JE,K)=W(IE,JE-1,K)
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
        W(IE,JE,K)=W(IE,JE+1,K)
      END IF
  246 CONTINUE
      DO 250 K=1,KBM1
      DO 250 J=1,JM
      DO 250 I=1,IM
      W(I,J,K)=W(I,J,K)*FSM(I,J)
 250  CONTINUE
      RETURN
C
C-------- Q2 & Q2L BOUNDARY CONDITIONS ---------------------------------
  60  CONTINUE
      DO 255 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      DO 255 K=1,KBM1
      UF(IC,JC,K)=UF(ID,JD,K)
      VF(IC,JC,K)=VF(ID,JD,K)
      L(IC,JC,K)=L(ID,JD,K)
      KM(IC,JC,K)=KM(ID,JD,K)
      KH(IC,JC,K)=KH(ID,JD,K)
      KQ(IC,JC,K)=KQ(ID,JD,K)
  255 CONTINUE
      DO 256 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
      DO 256 K=1,KBM1
      IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        UF(IE,JE,K)=UF(IE-1,JE,K)
        VF(IE,JE,K)=VF(IE-1,JE,K)
        L(IE,JE,K)=L(IE-1,JE,K)
        KM(IE,JE,K)=KM(IE-1,JE,K)
        KH(IE,JE,K)=KH(IE-1,JE,K)
        KQ(IE,JE,K)=KQ(IE-1,JE,K)
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        UF(IE,JE,K)=UF(IE+1,JE,K)
        VF(IE,JE,K)=VF(IE+1,JE,K)
        L(IE,JE,K)=L(IE+1,JE,K)
        KM(IE,JE,K)=KM(IE+1,JE,K)
        KH(IE,JE,K)=KH(IE+1,JE,K)
        KQ(IE,JE,K)=KQ(IE+1,JE,K)
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
        UF(IE,JE,K)=UF(IE,JE-1,K)
        VF(IE,JE,K)=VF(IE,JE-1,K)
        L(IE,JE,K)=L(IE,JE-1,K)
        KM(IE,JE,K)=KM(IE,JE-1,K)
        KH(IE,JE,K)=KH(IE,JE-1,K)
        KQ(IE,JE,K)=KQ(IE,JE-1,K)
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
        UF(IE,JE,K)=UF(IE,JE+1,K)
        VF(IE,JE,K)=VF(IE,JE+1,K)
        L(IE,JE,K)=L(IE,JE+1,K)
        KM(IE,JE,K)=KM(IE,JE+1,K)
        KH(IE,JE,K)=KH(IE,JE+1,K)
        KQ(IE,JE,K)=KQ(IE,JE+1,K)
      END IF
  256 CONTINUE
      DO 300 K=1,KB
      DO 300 J=1,JM
      DO 300 I=1,IM
      UF(I,J,K)=UF(I,J,K)*FSM(I,J)                                      
      VF(I,J,K)=VF(I,J,K)*FSM(I,J)
      L (I,J,K)=L (I,J,K)*FSM(I,J)
      KM(I,J,K)=KM(I,J,K)*FSM(I,J)
      KH(I,J,K)=KH(I,J,K)*FSM(I,J)
      KQ(I,J,K)=KQ(I,J,K)*FSM(I,J)
 300  CONTINUE
      RETURN
C
   70 CONTINUE
 4000 FORMAT(8E14.7)
C 
      IF(NUMEBC.EQ.0) GO TO 4045
      IF(OPTEBC(1:5).EQ.'DATA ') THEN
      IF(THOUR.LT.T2E) GO TO 4020
C 
      T1E=T2E 
      DO 4030 N=1,NUMEBC
 4030 DEBDRY(N,1)=DEBDRY(N,2) 
      READ (IUT90,4000,END=4010) T2E
      READ (IUT90,4000         ) (DEBDRY(N,2),N=1,NUMEBC)
C
 4020 CONTINUE
      FACT=(THOUR-T1E)/(T2E-T1E)
      DO 4040 N=1,NUMEBC
 4040 EBDRY(N)=DEBDRY(N,1)+FACT*(DEBDRY(N,2)-DEBDRY(N,1)) 
      END IF
 4045 CONTINUE
C 
      IF (NUMEBC.EQ.0) GOTO 4135
      IF(THOUR.LT.T2TS) GOTO 4111
C
      T1TS=T2TS
      DO 4121 N=1,NUMEBC
      DO 4122 K=1,KBM1
        DTBDRY(N,K,1)=DTBDRY(N,K,2)
        DSBDRY(N,K,1)=DSBDRY(N,K,2)
        DCCCBDRY(N,K,1)=DCCCBDRY(N,K,2)
4122   CONTINUE
4121   CONTINUE
      READ (IUT94,4000,END=4011) T2TS
      DO 4136 N=1,NUMEBC
      READ (IUT94,4000) (DTBDRY(N,K,2),K=1,KBM1)
      READ (IUT94,4000) (DSBDRY(N,K,2),K=1,KBM1)
      READ (IUT94,4000) (DCCCBDRY(N,K,2),K=1,KBM1)
4136  CONTINUE
4111  CONTINUE
      FACT=(THOUR-T1TS)/(T2TS-T1TS)
      DO 4042 N=1,NUMEBC
      DO 4042 K=1,KBM1
      TBDRY(N,K)=DTBDRY(N,K,1)+FACT*(DTBDRY(N,K,2)-DTBDRY(N,K,1)) 
      SBDRY(N,K)=DSBDRY(N,K,1)+FACT*(DSBDRY(N,K,2)-DSBDRY(N,K,1)) 
      CCCBDRY(N,K)=DCCCBDRY(N,K,1)+FACT*
     &(DCCCBDRY(N,K,2)-DCCCBDRY(N,K,1))
 4042 continue
4135  CONTINUE
C
      IF(NUMQBC.EQ.0) GO TO 4075
      IF(THOUR.LT.T2Q) GO TO 4050
C 
      T1Q=T2Q 
      DO 4060 N=1,NUMQBC
      DQDIS(N,1)=DQDIS(N,2) 
      DTDIS(N,1)=DTDIS(N,2) 
      DSDIS(N,1)=DSDIS(N,2) 
 4060 DCCCDIS(N,1)=DCCCDIS(N,2) 
C
      READ (IUT91,4000,END=4012) T2Q
      READ (IUT91,4000         ) (DQDIS(N,2),N=1,NUMQBC)
      READ (IUT91,4000         ) (DTDIS(N,2),N=1,NUMQBC)
      READ (IUT91,4000         ) (DSDIS(N,2),N=1,NUMQBC)
      READ (IUT91,4000         ) (DCCCDIS(N,2),N=1,NUMQBC)
C
 4050 CONTINUE
      FACT=(THOUR-T1Q)/(T2Q-T1Q)
      DO 4070 N=1,NUMQBC
      QDIS(N)=DQDIS(N,1)+FACT*(DQDIS(N,2)-DQDIS(N,1)) 
      TDIS(N)=DTDIS(N,1)+FACT*(DTDIS(N,2)-DTDIS(N,1)) 
      SDIS(N)=DSDIS(N,1)+FACT*(DSDIS(N,2)-DSDIS(N,1)) 
 4070 CCCDIS(N)=DCCCDIS(N,1)+FACT*(DCCCDIS(N,2)-DCCCDIS(N,1)) 
C
 4075 CONTINUE
C 
c      write(6,*) 'bcond numdbc=',numdbc
      IF(NUMDBC.EQ.0) GO TO 4105
      IF(THOUR.LT.T2D) GO TO 4080
C 
      T1D=T2D
      DO 4090 N=1,NUMDBC
      DQDIFF(N,1)=DQDIFF(N,2) 
      DTDIFF(N,1)=DTDIFF(N,2) 
      DSDIFF(N,1)=DSDIFF(N,2) 
 4090 DCCCDIFF(N,1)=DCCCDIFF(N,2) 
C
      READ (IUT92,4000,END=4014) T2D
      READ (IUT92,4000         ) (DQDIFF(N,2),N=1,NUMDBC)
      READ (IUT92,4000         ) (DTDIFF(N,2),N=1,NUMDBC)
      READ (IUT92,4000         ) (DSDIFF(N,2),N=1,NUMDBC)
      READ (IUT92,4000         ) (DCCCDIFF(N,2),N=1,NUMDBC)
C
 4080 CONTINUE
      FACT=(THOUR-T1D)/(T2D-T1D)
      DO 4100 N=1,NUMDBC
      QDIFF(N)=DQDIFF(N,1)+FACT*(DQDIFF(N,2)-DQDIFF(N,1)) 
      TDIFF(N)=DTDIFF(N,1)+FACT*(DTDIFF(N,2)-DTDIFF(N,1)) 
      SDIFF(N)=DSDIFF(N,1)+FACT*(DSDIFF(N,2)-DSDIFF(N,1)) 
 4100 CCCDIFF(N)=DCCCDIFF(N,1)+FACT*(DCCCDIFF(N,2)-DCCCDIFF(N,1)) 
C
 4105 CONTINUE
C 
      IF(THOUR.LT.T2M) GO TO 4110
C 
      T1M=T2M
      DQPREC(1)=DQPREC(2) 
      DQEVAP(1)=DQEVAP(2) 
      DTX   (1)=DTX   (2) 
      DTY   (1)=DTY   (2) 
      DHFLUX(1)=DHFLUX(2)
      READ (IUT93,4000,END=4016) T2M
      READ (IUT93,4000         ) DQPREC(2),DQEVAP(2),DTX(2),DTY(2),
     .                           DHFLUX(2)
C
 4110 CONTINUE
      FACT=(THOUR-T1M)/(T2M-T1M)
      QPREC=DQPREC(1)+FACT*(DQPREC(2)-DQPREC(1)) 
      QEVAP=DQEVAP(1)+FACT*(DQEVAP(2)-DQEVAP(1)) 
      TX   =DTX   (1)+FACT*(DTX   (2)-DTX   (1)) 
      TY   =DTY   (1)+FACT*(DTY   (2)-DTY   (1)) 
      HFLUX=DHFLUX(1)+FACT*(DHFLUX(2)-DHFLUX(1))
C
      DO 4120 J=1,JM
      DO 4120 I=1,IM
      WUSURF(I,J)=-1.E-3*(+TX*COS(ANG(I,J))+TY*SIN(ANG(I,J)))*FSM(I,J)
 4120 WVSURF(I,J)=-1.E-3*(-TX*SIN(ANG(I,J))+TY*COS(ANG(I,J)))*FSM(I,J)
C
      RETURN
C
 4010 WRITE(IUPRT,4130) THOUR
      WRITE(6,4130) THOUR
 4130 FORMAT(//' THE MODEL HAS RUN OUT OF ELEVATION DATA AT TIME '
     .    ,F10.4,' Hours'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GO TO 4140
C
4011  WRITE(IUPRT,4132) THOUR
      WRITE(6,4132) THOUR
4132  FORMAT(//' THE MODEL HAS RUN OUT OF TEMPERATURE-SALINITY DATA AT
     . TIME ' ,F10.4,' HOURS'/,'       REVISE INPUT DECK AND
     . RESUBMIT',//)
      GO TO 4140
C
 4012 WRITE(IUPRT,4013) THOUR
      WRITE(6,4013) THOUR
 4013 FORMAT(//' THE MODEL HAS RUN OUT OF DISCHARGE DATA AT TIME '
     .    ,F10.4,' Hours'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GO TO 4140

 4014 WRITE(6,4015) THOUR
      WRITE(IUPRT,4015) THOUR
 4015 FORMAT(//' THE MODEL HAS RUN OUT OF DIFFUSER DATA AT TIME '
     .    ,F10.4,' Hours'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GO TO 4140

 4016 WRITE(6,4017) THOUR
      WRITE(IUPRT,4017) THOUR
 4017 FORMAT(//' THE MODEL HAS RUN OUT OF METEORLOGICAL DATA AT TIME '
     .    ,F10.4,' Hours'/,'       REVISE INPUT DECK AND RESUBMIT '//)
C
 4140 CONTINUE
      CLOSE (IUT90)
      CLOSE (IUT91)
      CLOSE (IUT92)
      CLOSE (IUT93)
      CLOSE (IUT94)
      CALL SYSTEM('rm gcm_temp9*')
C
      STOP
      END
      

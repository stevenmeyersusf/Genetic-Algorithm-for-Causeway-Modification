      SUBROUTINE PRINTS(DRHOX,DRHOY,TRNU,TRNV)
C     VERSION(11/18/90)
      INCLUDE 'comdeck'
C
      DIMENSION DRHOX(IM,JM,KB),DRHOY(IM,JM,KB),TRNU(IM,JM),TRNV(IM,JM)
      DIMENSION IVAR(IM,JM),PRT(IM,KB)
      EQUIVALENCE (IVAR,C),(PRT,A)
C
      IF(INT.EQ.0) THEN
C
C----PRINT INITIAL FIELDS -----------------------------------------------
       WRITE(IUPRT,900)
       CALL PRINT(H,FSM,IM,JM,15,IVAR,10.,IUPRT)
       WRITE(IUPRT,913)
       CALL PRINT(FSM,FSM,IM,JM,15,IVAR,1.0,IUPRT)
       WRITE(IUPRT,914)
       CALL PRINT(DUM,DUM,IM,JM,15,IVAR,1.0,IUPRT)
       WRITE(IUPRT,915)
       CALL PRINT(DVM,DVM,IM,JM,15,IVAR,1.0,IUPRT)
       WRITE(IUPRT,901)
       CALL PRINT(H1,FSM,IM,JM,15,IVAR,1.,IUPRT)
       WRITE(IUPRT,902)
       CALL PRINT(H2,FSM,IM,JM,15,IVAR,1.,IUPRT)
       WRITE(IUPRT,903)
       CALL PRINT(ANG,FSM,IM,JM,15,IVAR,1.E2,IUPRT)
       WRITE(IUPRT,904)
       CALL PRINT(COR,FSM,IM,JM,15,IVAR,1.E8,IUPRT)
       WRITE(IUPRT,905)
       CALL PRINT(CBC,FSM,IM,JM,15,IVAR,1.E4,IUPRT)
       WRITE(IUPRT,953)
       CALL PRTXY(S ,FSM,IM,JM,KBM1,IVAR,100.,IUPRT)
       WRITE(IUPRT,954)
       CALL PRTXY(T ,FSM,IM,JM,KBM1,IVAR,100.,IUPRT)

       WRITE(IUPRT,987)
       CALL PRTXY(CCC ,FSM,IM,JM,KBM1,IVAR,100.,IUPRT)
C
C-------- COMPUTE TIME STEP RESTRICTIONS -------------------------------
      DO 500 J=1,JM
      DO 500 I=1,IM
 500  A(I,J,1)=.5/SQRT(GRAV*H(I,J)*(1./H1(I,J)**2
     2  +1./H2(I,J)**2))*FSM(I,J)
C
       WRITE(IUPRT,502)
       CALL PRINT(A(1,1,1),FSM,IM,JM,15,IVAR,10.,IUPRT)
C
C------- SPLIT INTO TWO DIRECTIONS FOR CHANNEL MODEL -------------------
      DO 510 J=1,JM
      DO 510 I=1,IM
 510  A(I,J,1)=.5/SQRT(GRAV*H(I,J))*H1(I,J)*FSM(I,J)
C
       WRITE(IUPRT,512)
       CALL PRINT(A(1,1,1),FSM,IM,JM,17,IVAR,1.,IUPRT)
C
      DO 520 J=1,JM
      DO 520 I=1,IM
 520  A(I,J,1)=.5/SQRT(GRAV*H(I,J))*H2(I,J)*FSM(I,J)
C
       WRITE(IUPRT,522)
       CALL PRINT(A(1,1,1),FSM,IM,JM,17,IVAR,1.,IUPRT)
C
C-------- INTERNAL MODE ------------------------------------------------
      DO 530 J=1,JM
      DO 530 I=1,IM
 530  A(I,J,1)=.5/(1.75*SQRT(1./H1(I,J)**2
     2  +1./H2(I,J)**2))*FSM(I,J)
C
       WRITE(IUPRT,532)
       CALL PRINT(A(1,1,1),FSM,IM,JM,15,IVAR,10.,IUPRT)
C
      DO 540 J=1,JM
      DO 540 I=1,IM
 540  A(I,J,1)=1./1.75*H1(I,J)*FSM(I,J)
C
       WRITE(IUPRT,542)
       CALL PRINT(A(1,1,1),FSM,IM,JM,17,IVAR,1.,IUPRT)
C
      DO 550 J=1,JM
      DO 550 I=1,IM
 550  A(I,J,1)=1./1.75*H2(I,J)*FSM(I,J)
C
       WRITE(IUPRT,552)
       CALL PRINT(A(1,1,1),FSM,IM,JM,17,IVAR,1.,IUPRT)
C
C-----------------------------------------------------------------------
      ELSE
C
       WRITE(IUPRT,910) TIME 
       CALL PRINT(EL,FSM,IM,JM,17,IVAR,1000.,IUPRT)
       WRITE(IUPRT,920) TIME 
       CALL PRINT(UA,DUM,IM,JM,17,IVAR,10000.,IUPRT)
       WRITE(IUPRT,930) TIME 
       CALL PRINT(VA,DVM,IM,JM,17,IVAR,10000.,IUPRT)
c     WRITE(IUPRT,975) TIME 
c     CALL PRINT(WUBOT,DUM,IM,JM,17,IVAR,1.E7,IUPRT) 
c     WRITE(IUPRT,980) TIME 
c     CALL PRINT(WVBOT,DVM,IM,JM,17,IVAR,1.E7,IUPRT) 
C
      IF(TOR.EQ.'BAROTROPIC') RETURN
C     WRITE(IUPRT,933) TIME 
C     CALL PRINT(TRNU,DUM,IM,JM,17,IVAR,.01,IUPRT)
C     WRITE(IUPRT,934) TIME 
C     CALL PRINT(TRNV,DVM,IM,JM,17,IVAR,.01,IUPRT)
       WRITE(IUPRT,935) TIME 
       CALL PRTXY(U ,DUM,IM,JM,KBM1,IVAR,10000.,IUPRT) 
       CALL SLICEXZ(U,DUM,IM,JM,KB, 4,PRT,100.,IUPRT) 
       CALL SLICEXZ(U,DUM,IM,JM,KB,19,PRT,100.,IUPRT) 
       CALL SLICEXZ(U,DUM,IM,JM,KB,30,PRT,100.,IUPRT) 
       WRITE(IUPRT,940) TIME 
      CALL PRTXY(V ,DVM,IM,JM,KBM1,IVAR,10000.,IUPRT) 
      CALL SLICEYZ(V,DVM,IM,JM,KB, 3,PRT,100.,IUPRT) 
      CALL SLICEYZ(V,DVM,IM,JM,KB, 5,PRT,100.,IUPRT) 
      CALL SLICEYZ(V,DVM,IM,JM,KB,14,PRT,100.,IUPRT) 
      CALL SLICEYZ(V,DVM,IM,JM,KB,19,PRT,100.,IUPRT) 
      WRITE(IUPRT,945) TIME
C     CALL PRTXY(WR,FSM,IM,JM,KB,IVAR,1.E4,IUPRT) 
      CALL SLICEXZ(WR,FSM,IM,JM,KB, 4,PRT,1.E4,IUPRT) 
C
      IF(TOR.EQ.'PROGNOSTIC') THEN
      WRITE(IUPRT,951) TIME 
      CALL PRTXY(S ,FSM,IM,JM,KBM1,IVAR,100.,IUPRT) 
      CALL SLICEXZ(S,FSM,IM,JM,KB,4,PRT,1.,IUPRT) 
      CALL SLICEXZ(S,FSM,IM,JM,KB,19,PRT,1.,IUPRT) 
      CALL SLICEXZ(S,FSM,IM,JM,KB,30,PRT,1.,IUPRT) 
      CALL SLICEYZ(S,FSM,IM,JM,KB, 3,PRT,1.,IUPRT) 
      CALL SLICEYZ(S,FSM,IM,JM,KB, 5,PRT,1.,IUPRT) 
      CALL SLICEYZ(S,FSM,IM,JM,KB,19,PRT,1.,IUPRT) 
      WRITE(IUPRT,950) TIME 
      WRITE(IUPRT,950) TIME 
C     CALL PRTXY(T ,FSM,IM,JM,KBM1,IVAR,100.,IUPRT) 
c     CALL SLICEYZ(T,FSM,IM,JM,KB,25,PRT,1.,IUPRT) 
c     CALL SLICEYZ(T,FSM,IM,JM,KB,28,PRT,1.,IUPRT) 
      CALL SLICEXZ(T,FSM,IM,JM,KB,4,PRT,1.,IUPRT) 
      WRITE(IUPRT,952) TIME 
C     CALL PRTXY(RHO,FSM,IM,JM,KBM1,IVAR,1.E5,IUPRT)
      CALL SLICEYZ(RHO,FSM,IM,JM,KB,4,PRT,1.E3,IUPRT)
      END IF
C
      WRITE(IUPRT,985) TIME 
c     CALL PRTXY(AAM,FSM,IM,JM,KBM1,IVAR,100.,IUPRT)
      CALL SLICEXZ(AAM,FSM,IM,JM,KB, 4,PRT,1.,IUPRT)
      WRITE(IUPRT,955) TIME 
C     CALL PRTXY(Q2,FSM,IM,JM,KB,IVAR,1.E6,IUPRT)
      CALL SLICEXZ(Q2,FSM,IM,JM,KB, 4,PRT,1.E4,IUPRT)
C     CALL SLICEXZ(Q2,FSM,IM,JM,KB, 9,PRT,1.E4,IUPRT)
C     WRITE(IUPRT,960) TIME 
C     CALL PRTXY(L  ,FSM,IM,JM,KB,IVAR,100.,IUPRT)
C     CALL SLICEXZ(L,FSM,IM,JM,KB, 4,PRT,1.,IUPRT)
C     CALL SLICEXZ(L,FSM,IM,JM,KB, 9,PRT,1.,IUPRT)
      WRITE(IUPRT,965) TIME 
C     CALL PRTXY(KM ,FSM,IM,JM,KB,IVAR,1.E4,IUPRT)
      CALL SLICExZ(KM,FSM,IM,JM,KB, 4,PRT,1.E4,IUPRT) 
C     WRITE(IUPRT,970) TIME 
C     CALL PRTXY(KH ,FSM,IM,JM,KB,IVAR,1.E4,IUPRT)
C     CALL SLICEXZ(KH,FSM,IM,JM,KB, 4,PRT,1.E4,IUPRT) 
C     CALL SLICEXZ(KH,FSM,IM,JM,KB, 9,PRT,1.E4,IUPRT) 
C     WRITE(IUPRT,971) TIME 
C     CALL PRTXY(DRHOX,DUM,IM,JM,KB,IVAR,1.,IUPRT)
C     CALL SLICEXZ(DRHOX,DUM,IM,JM,KB, 4,PRT,1.,IUPRT)
C     CALL SLICEXZ(DRHOX,DUM,IM,JM,KB, 9,PRT,1.,IUPRT)
C     WRITE(IUPRT,972)TIME
C     CALL PRTXY(DRHOY,DVM,IM,JM,KB,IVAR,1.,IUPRT)
C     CALL SLICEXZ(DRHOY,DVM,IM,JM,KB,4,PRT,.01,IUPRT)
C     CALL SLICEXZ(DRHOY,DVM,IM,JM,KB,31,PRT,.01,IUPRT)
C 
      END IF
C 
 900  FORMAT(///24H .. DEPTH (m)...........,/)
 901  FORMAT(///24H ..... METRIC H1 (m)....,/)
 902  FORMAT(///24H ..... METRIC H2 (m)....,/)
 903  FORMAT(///'  ..... ANGLE (radians) ....',/)
 904  FORMAT(///' .... CORIOLIS PARAMETER (1/sec)....',/)
 905  FORMAT(///' .... BOTTOM DRAG COEFFICIENT ....',/)
 913  FORMAT(///24H . ELEVATION MASK.......,/)
 914  FORMAT(///24H . U  VEL MASK..........,/)
 915  FORMAT(///24H . V  VEL MASK..........,/)
C
 502  FORMAT(//'MAXIMUM TIME STEP(SEC) EXTERNAL MODE '//)
 512  FORMAT(//'MAXIMUM TIME STEP(SEC) EXTERNAL MODE XI 1 DIRECTION'//)
 522  FORMAT(//'MAXIMUM TIME STEP(SEC) EXTERNAL MODE XI 2 DIRECTION'//)
 532  FORMAT(//'MAXIMUM TIME STEP(SEC) INTERNAL MODE '//)
 542  FORMAT(//'MAXIMUM TIME STEP(SEC) INTERNAL MODE XI 1 DIRECTION'//)
 552  FORMAT(//'MAXIMUM TIME STEP(SEC) INTERNAL MODE XI 2 DIRECTION'//)
 910  FORMAT(///' SURFACE ELEVATION (M) AT ',F9.4,15H DAYS AFTER T=0,/)
 920  FORMAT(///' AVERAGED U VELOCITY (M/S)',F9.4,15H DAYS AFTER T=0,/)
 930  FORMAT(///' AVERAGED V VELOCITY (M/S)',F9.4,15H DAYS AFTER T=0,/)
 931  FORMAT(//' TOTAL X-INTEGRALS AT ',F9.4,' DAYS AFTER T=0',//)
 932  FORMAT(//' TOTAL Y-INTEGRALS AT ',F9.4,' DAYS AFTER T=0',//)
 933  FORMAT(//' U-INTEGRALS AT',F9.4,15H DAYS AFTER T=0,//)
 934  FORMAT(//' V-INTEGRALS AT',F9.4,15H DAYS AFTER T=0,//)
 935  FORMAT(//' TOTAL U VELOCITY(M/S)',F9.4,15H DAYS AFTER T=0,//) 
 940  FORMAT(//' TOTAL V VELOCITY(M/S)',F9.4,15H DAYS AFTER T=0,//)
 945  FORMAT(//' TOTAL W VELOCITY(M/S)',F9.4,15H DAYS AFTER T=0,//)
 950  FORMAT(//' TEMPERATURE (C)'   ,F9.4,15H DAYS AFTER T=0,//)
 954  FORMAT(//' TEMPERATURE (C)  INITIALIZATION ',//)
 951  FORMAT(//' SALINITY(PPT)',F9.4,15H DAYS AFTER T=0,//)
 953  FORMAT(//' SALINITY(PPT)  INITIALIZATION ',//)
 952  FORMAT(//' DENSITY - 1. (GM/CM**3)' ,F9.4,15H DAYS AFTER T=0,//)
 955  FORMAT(//' TURBULENT K.E.',F9.4,15H DAYS AFTER T=0,//)
 960  FORMAT(//' MIXING LENGTH(M)',F9.4,15H DAYS AFTER T=0,//)
 965  FORMAT(//' MIXING KM(M**2/S)',F9.4,15H DAYS AFTER T=0,//)
 970  FORMAT(//' MIXING KH(M**2/S)',F9.4,15H DAYS AFTER T=0,//)
 971  FORMAT(//22H BAROCLN. PRES. DRHOX ,F9.4,15H DAYS AFTER T=0,//)
 972  FORMAT(//22H BAROCLN. PRES. DRHOY ,F9.4,15H DAYS AFTER T=0,//)
 975  FORMAT(//' U BOTTOM STRESS((M/S)**2)',F9.4,15H DAYS AFTER T=0,//)
 980  FORMAT(//' V BOTTOM STRESS((M/S)**2)',F9.4,15H DAYS AFTER T=0,//)
 985  FORMAT(//' HORIZONTAL MIXING        ',F9.4,15H DAYS AFTER T=0,//)
 987  FORMAT(//' Concentration(percent)  INITIALIZATION ',//)
C
      RETURN
      END
      
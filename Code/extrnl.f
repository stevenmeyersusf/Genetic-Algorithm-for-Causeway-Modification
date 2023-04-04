      SUBROUTINE EXTRNL(ADVUA,ADVVA,TRNU,TRNV,DTE2,BFRIC,DTI2)
C     VERSION(06/07/90)
      INCLUDE 'comdeck'
C
C-----------------------------------------------------------------------
C     CALCULATE FREE SURFACE HEIGHT
C-----------------------------------------------------------------------
C
      DIMENSION ADVUA(IM,JM),ADVVA(IM,JM),TRNU(IM,JM),TRNV(IM,JM)
C
      dmin = 100.

      DO 400 J=2,JMM1
      DO 400 I=2,IM
 400  FLUXUA(I,J)=.25*(D(I,J)+D(I-1,J))*(H2(I,J)+H2(I-1,J))*UA(I,J)
      DO 405 J=2,JM
      DO 405 I=2,IMM1
 405  FLUXVA(I,J)=.25*(D(I,J)+D(I,J-1))*(H1(I,J)+H1(I,J-1))*VA(I,J)
C
      i0 = 35
      j0 = 50
      DO 410 J=2,JMM1
      DO 410 I=2,IMM1
      ELF(I,J)=ELB(I,J)-DTE2* 
     . (FLUXUA(I+1,J)-FLUXUA(I,J)+FLUXVA(I,J+1)-FLUXVA(I,J))/ART(I,J)  
 410  CONTINUE

C-----------------------------------------------------------------------
C         IMPOSE MASS FLUX BOUNDARY CONDITIONS
C-----------------------------------------------------------------------
C
      DO 406 N=1,NUMDBC
      ID=IDD(N)
      JD=JDD(N)
      ELF(ID,JD)=ELF(ID,JD)+ DTE2*QDIFF(N)*RAMP /ART(ID,JD)
  406 CONTINUE
c      if (int.ge.1576997) write(6,*) 'ext',art(i0,j0),elf(i0,j0)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C         IMPOSE MASS FLUX BOUNDARY CONDITIONS FROM PREC. AND EVAP.
C         MSV 7-22-98
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DO 407 J=2,JMM1
      DO 407 I=2,IMM1
      ELF(I,J)=ELF(I,J)+ DTE2*(QPREC-QEVAP)*RAMP
  407 CONTINUE

c      if (int.ge.1576997) write(6,*) 'ext',elf(i0,j0),qprec,qevap
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL BCOND2(1,DTI2)
      CALL ADVAVE(ADVUA,ADVVA)

      DO 95 J=1,JM
      DO 95 I=1,IM
 95   CURV42D(I,J)=0.25*COR(I,J)
C
C-----------------------------------------------------------------------
C        COMPUTE BOTTOM FRICTION AND CURVATURE TERMS 
C            IF RUN IS BAROTROPIC MODE ONLY
C-----------------------------------------------------------------------
C
      IF(TOR.NE.'BAROTROPIC') GO TO 5000
C
      DO 60 J=2,JMM1
      DO 60 I=2,IMM1
      CURV42D(I,J)= CURV42D(I,J)
     .  +.125*((VA(I,J+1)+VA(I,J))*
     .     ((H2(I+1,J)*FSM(I+1,J)+H2(I,J)*FSM(I,J))/
     .          (FSM(I+1,J)+FSM(I,J)+1.E-30)-
     .      (H2(I,J)*FSM(I,J)+H2(I-1,J)*FSM(I-1,J))/
     .          (FSM(I,J)+FSM(I-1,J)+1.E-30)) 
     .  -(UA(I+1,J)+UA(I,J))*
     .     ((H1(I,J+1)*FSM(I,J+1)+H1(I,J)*FSM(I,J))/
     .          (FSM(I,J+1)+FSM(I,J)+1.E-30)-
     .      (H1(I,J)*FSM(I,J)+H1(I,J-1)*FSM(I,J-1))/
     .          (FSM(I,J)+FSM(I,J-1)+1.E-30)) )
     .                /(H1(I,J)*H2(I,J))
 60   CONTINUE
C
      DO 100 J=2,JMM1
      DO 100 I=2,IMM1
      WUBOT(I,J)=-BFRIC
     .     *SQRT(UAB(I,J)**2+(.25*(VAB(I,J)
     .     +VAB(I,J+1)+VAB(I-1,J)+VAB(I-1,J+1)))**2)*UAB(I,J)
 100  CONTINUE
      DO 102 J=2,JMM1
      DO 102 I=2,IMM1
      WVBOT(I,J)=-BFRIC
     .     *SQRT((.25*(UAB(I,J)+UAB(I+1,J)
     .     +UAB(I,J-1)+UAB(I+1,J-1)))**2+VAB(I,J)**2)*VAB(I,J)
 102  CONTINUE
 5000 CONTINUE
C
C-----------------------------------------------------------------------
C     CALCULATE U COMPONENT OF VELOCITY
C-----------------------------------------------------------------------
C
      DO 420 J=2,JMM1
      DO 420 I=3,IMM1
 420  UAF(I,J)=ADVUA(I,J)
     .    -ARU(I,J)*(  CURV42D(I,J)*D(I,J)*(VA(I,J+1)+VA(I,J))
     .              +CURV42D(I-1,J)*D(I-1,J)*(VA(I-1,J+1)+VA(I-1,J)) )
     .         +.25*GRAV*(H2(I,J)+H2(I-1,J))*(D(I,J)+D(I-1,J))
     .             *( (1.-2.*ALPHA)*(EL(I,J)-EL(I-1,J))
     .            +ALPHA*(ELB(I,J)-ELB(I-1,J)+ELF(I,J)-ELF(I-1,J)) )
     .           +RAMP*TRNU(I,J)
     .      -ARU(I,J)*(-.5*(WUSURF(I,J)+WUSURF(I-1,J))+WUBOT(I,J)  )
      DO 425 J=2,JMM1
      DO 425 I=3,IMM1
      UAF(I,J)=
     .         ((H(I,J)+ELB(I,J)+H(I-1,J)+ELB(I-1,J))*ARU(I,J)*UAB(I,J)
     .                -4.*DTE*UAF(I,J))
     .        /((H(I,J)+ELF(I,J)+H(I-1,J)+ELF(I-1,J))*ARU(I,J))

 425  continue
      
C
C-----------------------------------------------------------------------
C     CALCULATE V COMPONENT OF VELOCITY
C-----------------------------------------------------------------------
C
      DO 430 J=3,JMM1
      DO 430 I=2,IMM1
 430  VAF(I,J)=ADVVA(I,J)  
     .    +ARV(I,J)*(  CURV42D(I,J)*D(I,J)*(UA(I+1,J)+UA(I,J))
     .                +CURV42D(I,J-1)*D(I,J-1)*(UA(I+1,J-1)+UA(I,J-1)) )
     .         +.25*GRAV*(H1(I,J)+H1(I,J-1))*(D(I,J)+D(I,J-1))
     .                   *( (1.-2.*ALPHA)*(EL(I,J)-EL(I,J-1))
     .            +ALPHA*(ELB(I,J)-ELB(I,J-1)+ELF(I,J)-ELF(I,J-1)) )
     .                   +RAMP*TRNV(I,J)
     .    + ARV(I,J)*( .5*(WVSURF(I,J)+WVSURF(I,J-1))-WVBOT(I,J)   )
      DO 435 J=3,JMM1
      DO 435 I=2,IMM1
      VAF(I,J)=
     .        ((H(I,J)+ELB(I,J)+H(I,J-1)+ELB(I,J-1))*VAB(I,J)*ARV(I,J)
     .              -4.*DTE*VAF(I,J))
     .       /((H(I,J)+ELF(I,J)+H(I,J-1)+ELF(I,J-1))*ARV(I,J))
C  

 435  continue

      CALL BCOND2(2,DTI2)
C

c    2012 wet-dry implementation
      if (ido_wd.gt.0) then 
         ndry = 0
         DO 408 J=2,JMm1          ! CHECK WET-DRY 2012
         DO 408 I=2,IMm1
            dij = h(i,j)+elf(i,j)
            if ((fsm(i,j).eq.1).and.(dij.le.hdry)) then   
               elf(i,j) = -h(i,j)+hdry-0.0000001 ! extra is for roundoff
               wd(i,j) = 0           
               ndry = ndry+1
c               write(6,*) 'DRY at ',i,j
               if (uaf(i,j).lt.0)  uaf(i,j)=0
               if (uaf(i+1,j).gt.0) uaf(i+1,j)=0
               if (vaf(i,j).lt.0) vaf(i,j)=0
               if (vaf(i,j+1).gt.0) vaf(i,j+1)=0
            endif       
 408     continue
cc         if (ndry.gt.0) write(6,*)'ndry=',ndry,int*dti/3600,' hr'
      endif                     !ido




      RETURN
      END
      

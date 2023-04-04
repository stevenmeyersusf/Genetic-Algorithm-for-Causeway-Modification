      SUBROUTINE ADVAVE(ADVUA,ADVVA)
C     VERSION(03/02/90)
      INCLUDE 'comdeck'
C
C-----------------------------------------------------------------------
C      THIS SUBROUTINE CALCULATES ADVECTION & DIFFUSION
C-----------------------------------------------------------------------
C
      DIMENSION ADVUA(IM,JM),ADVVA(IM,JM)
C
C-------- CALCULATE U-ADVECTION & DIFFUSION ----------------------------
C
C-------- ADVECTIVE FLUXES ---------------------------------------------
      DO 10 J=1,JM
      DO 10 I=1,IM
      FLUXUA(I,J)=0.0
      FLUXVA(I,J)=0.0
      ADVUA(I,J)=0.0
 10   TPS(I,J)=0.0
      IF (ADVECT.EQ.'NON-LINEAR') THEN
C
C-------- ADVECTIVE FLUXES ---------------------------------------------
      DO 300 J=2,JMM1
      DO 300 I=2,IMM1
 300  FLUXUA(I,J)=.125*((D(I+1,J)+D(I,J))*UA(I+1,J)
     .                 +(D(I,J)+D(I-1,J))*UA(I,J))
     .                  *(UA(I+1,J)+UA(I,J))
      DO 400 J=2,JM
      DO 400 I=2,IMM1
 400  FLUXVA(I,J)=.125*((D(I,J)+D(I,J-1))*VA(I,J)
     .                 +(D(I-1,J)+D(I-1,J-1))*VA(I-1,J))
     .                 *(UA(I,J)+UA(I,J-1))
      END IF
C
C-------- ADD VISCOUS FLUXES -------------------------------------------
      DO 460 J=2,JMM1
      DO 460 I=2,IMM1
 460  FLUXUA(I,J)=FLUXUA(I,J)
     .         -D(I,J)*2.*AAM2D(I,J)*(UAB(I+1,J)-UAB(I,J))/H1(I,J)
      DO 470 J=2,JM
      DO 470 I=2,IMM1
      TPS(I,J)=.25*(D(I,J)+D(I-1,J)+D(I,J-1)+D(I-1,J-1))
     .            *(AAM2D(I,J)+AAM2D(I,J-1)+AAM2D(I-1,J)+AAM2D(I-1,J-1))
     .                *((UAB(I,J)-UAB(I,J-1))
     .                /(H2(I,J)+H2(I-1,J)+H2(I,J-1)+H2(I-1,J-1))
     .                 +(VAB(I,J)-VAB(I-1,J))
     .               /(H1(I,J)+H1(I-1,J)+H1(I,J-1)+H1(I-1,J-1)) )
      FLUXUA(I,J)=FLUXUA(I,J)*H2(I,J)
      FLUXVA(I,J)=(FLUXVA(I,J)-TPS(I,J))
     .            *.25*(H1(I,J)+H1(I-1,J)+H1(I,J-1)+H1(I-1,J-1))

c      if (((wd(i,j).eq.0).and.(UA(i,j).lt.0)).or. 
c     &     ((wd(i-1,j).eq.0).and.(UA(i,j).gt.0))) then 
c         FLUXUA(i,j)=0.
c         write(6,*) 'fluxua fixing ',i,j
c      endif

c      if (((wd(i,j).eq.0).and.(VA(i,j).lt.0)).or.
c     &     ((wd(i,j-1).eq.0).and.(VA(i,j).gt.0)))then 
c         fLUXVA(i,j)=0.
c         write(6,*) 'fluxva fixing ',i,j
c      endif

 470  CONTINUE
C
      DO 480 J=2,JMM1
      DO 480 I=3,IMM1
 480  ADVUA(I,J)=FLUXUA(I,J)-FLUXUA(I-1,J)
     .           +FLUXVA(I,J+1)-FLUXVA(I,J)
C
C-------- CALCULATE V-ADVECTION & DIFFUSION ----------------------------
      DO 20 J=1,JM
      DO 20 I=1,IM
      FLUXUA(I,J)=0.0
      FLUXVA(I,J)=0.0
 20   ADVVA(I,J)=0.0
      IF (ADVECT.EQ.'NON-LINEAR') THEN
C
C-------- ADVECTIVE FLUXES ---------------------------------------------
      DO 700 J=3,JMM1
      DO 700 I=2,IM
 700  FLUXUA(I,J)=.125*((D(I,J)+D(I-1,J))*UA(I,J)
     .         +(D(I,J-1)+D(I-1,J-1))*UA(I,J-1))*
     .                        (VA(I-1,J)+VA(I,J))
      DO 800 J=2,JMM1
      DO 800 I=2,IMM1
 800  FLUXVA(I,J)=.125*((D(I,J+1)+D(I,J))
     .       *VA(I,J+1)+(D(I,J)+D(I,J-1))*VA(I,J))
     .      *(VA(I,J+1)+VA(I,J))
      END IF
C
C-------- ADD VISCOUS FLUXES -------------------------------------------
      DO 860 J=2,JMM1
      DO 860 I=2,IMM1
 860  FLUXVA(I,J)=FLUXVA(I,J)
     .        -D(I,J)*2.*AAM2D(I,J)*(VAB(I,J+1)-VAB(I,J))/H2(I,J)
      DO 870 J=2,JMM1
      DO 870 I=2,IM
      FLUXVA(I,J)=FLUXVA(I,J)*H1(I,J)
      FLUXUA(I,J)=(FLUXUA(I,J)-TPS(I,J))
     .             *.25*(H2(I,J)+H2(I-1,J)+H2(I,J-1)+H2(I-1,J-1))
 870  continue

C
      DO 880 J=3,JMM1
      DO 880 I=2,IMM1
      ADVVA(I,J)=FLUXUA(I+1,J)-FLUXUA(I,J)
     .          +FLUXVA(I,J)-FLUXVA(I,J-1)
c      if ((wd(i,j).eq.0).and.(UA(i,j).lt.0)) FLUXUA(i,j)=0.
c      if ((wd(i-1,j).eq.0).and.(UA(i,j).gt.0)) FLUXUA(i,j)=0.

c      if ((wd(i,j).eq.0).and.(VA(i,j).lt.0)) FLUXVA(i,j)=0.
c      if ((wd(i,j-1).eq.0).and.(VA(i,j).gt.0)) FLUXVA(i,j)=0.


 880  continue

C
      DO 310 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
      IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        ADVUA(IE,JE)=0.0
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        ADVUA(IE+1,JE)=0.0
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
        ADVVA(IE,JE)=0.0
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
        ADVVA(IE,JE+1)=0.0
      ENDIF
 310  CONTINUE
C
      RETURN
      END

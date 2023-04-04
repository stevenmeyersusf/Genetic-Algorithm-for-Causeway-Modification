      SUBROUTINE DENS
C     VERSION(03/02/90)
      INCLUDE 'comdeck'
C
C-----------------------------------------------------------------------
C         THIS FUNCTION COMPUTES DENSITY-1.0
C-----------------------------------------------------------------------
C
      REAL*8 RHOF(IM,JM,KB),TF(IM,JM,KB),SF(IM,JM,KB)
      EQUIVALENCE (A,RHOF),(VH,TF),(UF,SF)
C
      ROSS=0.0
      IF (ROSS.EQ.1.0) THEN
      DO 10 K=1,KB
      DO 10 J=1,JM
      DO 10 I=1,IM
   10 RHO(I,J,K)=0.0
      RETURN
      END IF
      DO 1 K=1,KBM1
      DO 1 J=1,JM
      DO 1 I=1,IM
      TF(I,J,K)=T(I,J,K)
c      TF(I,J,K)=25.0  ! no active T in density. 
      SF(I,J,K)=S(I,J,K)
      RHOF(I,J,K)=
     . SF(I,J,K)*SF(I,J,K)*SF(I,J,K)*6.76786136E-6-SF(I,J,K)*SF(I,J,K)*
     . 4.8249614E-4+SF(I,J,K)*8.14876577E-1-0.22584586E0
C
      RHOF(I,J,K)=RHOF(I,J,K)*
     . (TF(I,J,K)*TF(I,J,K)*TF(I,J,K)*1.667E-8-TF(I,J,K)*TF(I,J,K)
     . *8.164E-7+TF(I,J,K)*1.803E-5)
C
      RHOF(I,J,K)=RHOF(I,J,K)+1.- 
     . TF(I,J,K)*TF(I,J,K)*TF(I,J,K)*1.0843E-6+TF(I,J,K)*TF(I,J,K)
     . *9.8185E-5-TF(I,J,K)*4.786E-3
C
      RHOF(I,J,K)=RHOF(I,J,K)*
     . (SF(I,J,K)*SF(I,J,K)*SF(I,J,K)*6.76786136E-6-SF(I,J,K)*SF(I,J,K)*
     . 4.8249614E-4+SF(I,J,K)*8.14876577E-1+3.895414E-2)
C
      RHOF(I,J,K)=RHOF(I,J,K)-
     . (TF(I,J,K)-3.98)**2 * (TF(I,J,K)+283.)/(503.57*(TF(I,J,K)+67.26))
    1 CONTINUE
C
      DO 2 K=1,KBM1
      DO 2 J=1,JM
      DO 2 I=1,IM
      RHO(I,J,K)= RHOF(I,J,K)     *1.E-3*FSM(I,J)
    2 CONTINUE
      DO 3 J=1,JM
      DO 3 I=1,IM
    3 RHO(I,J,KB)=0.0  ! bottom
C      
      RETURN
      END
      

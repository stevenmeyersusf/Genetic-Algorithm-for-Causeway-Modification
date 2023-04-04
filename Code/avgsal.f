      SUBROUTINE AVGSAL(SAVG,SAVGTM,ISFLAG,IUSAL2)

c     calculate 25 hr average of salinity field
C     CHECK TIME. UNITS OF AVGTMS ARE DAYS FROM JAN 1, 1990

      include 'comdeck'

      PARAMETER(NAVGTM=19)
      DIMENSION SAVG(IM,JM,KB),AVGTMS(NAVGTM)
      INTEGER ISFLAG,CKFLAG,IUSAL2
      REAL SAVGTM

c     avgtms units are days
      DATA AVGTMS /80.,100.,125.,150.,175.,200.,250.,320.,360.,
     & 385.,400.,420.,450.,475.,515.,550.,567.,585.,620./   

      CKFLAG=0
      DO 100 M=1,NAVGTM
 100  IF (ABS(AVGTMS(M)-TIME)*24.*60.*60.LT.DTI/2.) CKFLAG=1

c      write(6,*) time,ckflag

C     NOT TIME TO DO AVERAGE
      IF ((CKFLAG.EQ.0).AND.(ISFLAG.EQ.0)) GOTO 999
 
C     START AVERAGING
      IF ((CKFLAG.EQ.1).AND.(ISFLAG.EQ.0)) THEN
        DO 200 I=1,IM
        DO 200 J=1,JM
        DO 200 K=1,KB
           SAVG(I,J,K)=S(I,J,K)
 200    CONTINUE
        ISFLAG=1
        write(6,*) '25hr ISFLAG=',ISFLAG
        write(6,*) 'TIME AVGTMS M=',TIME,AVGTMS(M),M        
      ENDIF

C     CONTINUE AVERAGING OVER TIME SAVGTM
      IF ((CKFLAG.EQ.0).AND.(ISFLAG.GT.0)) THEN
        DO 300 I=1,IM
        DO 300 J=1,JM
        DO 300 K=1,KB
           SAVG(I,J,K)=S(I,J,K)+SAVG(I,J,K)
 300    CONTINUE
        ISFLAG = ISFLAG+1

        IF (ABS(ISFLAG-SAVGTM*60.*60/DTI).LT.1) THEN 
           DO 400 I=1,IM
           DO 400 J=1,JM
           DO 400 K=1,KB
              SAVG(I,J,K)=SAVG(I,J,K)/ISFLAG
 400       CONTINUE
           WRITE(IUSAL2,*) 'TIME (DAYS)'
           WRITE(IUSAL2,*) TIME
           WRITE(IUSAL2,500) (((SAVG(I,J,K), I=1,IM),J=1,JM),K=1,KB)
           write(6,*) 'writing 25hr ',ISFLAG
           ISFLAG=0
        ENDIF

      ENDIF

 500  format(8(f9.3))

 999  CONTINUE

      END


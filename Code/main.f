c23456789012345678901234567890123456789012345678901234567890123456789012
C     VERSION(11/11/90)

      INCLUDE 'comdeck'
      include 'genetic'
C 
CMSV model_11-3.f  revise BC south of Egmont in SETDOM (11-3-95)
CMSV model_11-8.f  added driving of Bunces Pass and Passagrille (11-7-95)
CMSV               added options for real-time project i.e. seeded time and ramp          
CMSV               for a seeded cold start
CMSV               added ifpointfile and ifopenede and v  for tseries files 6-12-97
CMSV model_06_26_97nc.f revised TINE SETDOM to use 10 versus 11 obc and rewrote logic
Cdcb added drifter sub to avect particles with markovian diff. throughout
Cdcb the domain
Cdcb hardwired in time for first ccc print (arctccc)10/6/99

c    modified for genetic search of causeway cuts began 4/19/2022
c    by S. Meyers
c
C----------------------------------------------------------------------|
C                  GENERAL CIRCULATION MODEL                           |
C          ORTHOGONAL CURVILINEAR COORDINATE VERSION                   |
C                         ECOM3D                                       |
C                                                                      |
C                                                                      |
C     THIS IS A VERSION OF THE THREE DIMENSIONAL, TIME DEPENDENT,      |
C   PRIMITIVE EQUATION, CIRCULATION MODEL DEVELOPED BY GEORGE MELLOR   |
C   AND ALAN BLUMBERG WITH SUBSEQUENT CONTRIBUTIONS BY LEO OEY AND     |
C   BORIS GALPERIN.                                                    |
C                                                                      |
C                                                                      |
C     FOR DETAILS OF THE GOVERNING EQUATIONS AND SOLUTION              |
C     TECHNIQUES THE INTERESTED READER IS REFERRED TO:                 |
C                                                                      |
C     BLUMBERG, A.F. AND G.L. MELLOR, DIAGNOSTIC AND PROGNOSTIC        |
C     NUMERICAL CIRCULATION STUDIES OF THE SOUTH ATLANTIC BIGHT        |
C     J. GEOPHYS. RES. 88, 4579-4592, 1983.                            |
C                                                                      |
C     AND                                                              |
C                                                                      |
C     BLUMBERG, A.F. AND G.L. MELLOR, A DESCRIPTION OF A THREE         |
C     COASTAL OCEAN CIRCULATION MODEL, THREE DIMENSIONAL SHELF         |
C     MODELS, COASTAL AND ESTUARINE SCIENCES, 5, N.HEAPS, ED.,         |
C     AMERICAN GEOPHYSICAL UNION, 1987.                                |
C                                                                      |
C                                                                      |
C     IN SUBROUTINE PROFQ THE MODEL MAKES USE OF THE TURBULENCE        |
C     CLOSURE SUB-MODEL MOST RECENTLY DESCRIBED IN:                    |
C                                                                      |
C     MELLOR, G.L. AND T. YAMADA, DEVELOMENT OF A TURBULENCE CLOSURE   |
C     MODEL FOR GEOPHYSICAL FLUID PROBLEMS, REV. GEOPHYS. SPACE PHYS., |
C     851-879, 1982.                                                   |
C                                                                      |
C                                                                      |
C        THIS PROGRAM CONTAINS A NUMBER OF COMMENT CARDS THAT          |
C     SHOULD ENABLE AN ASTUTE OCEAN SCIENTIST TO MAKE SIMULATIONS.     |
C     PLEASE DIRECT CRITICISMS AND SUGGESTIONS TO  ME.                 |
C                                                                      |
C                                             ALAN BLUMBERG            |
C                                             HYDROQUAL, INC           |
C                                               MAHWAH, NJ             |
C                                               12/28/89               |
C                                                                      |
C----------------------------------------------------------------------|
C                                                                      |
C                   LOOP LIMITS                                        |
C                                                                      |
C              T,S,etc :  J=2,JMM1                                     |
C                         I=2,IMM1                                     |
C                                                                      |
C                 U    :  J=2,JMM1                                     |
C                         I=3,IMM1                                     |
C                                                                      |
C                 V    :  J=3,JMM1                                     |
C                         I=2,IMM1                                     |
C                                                                      |
C----------------------------------------------------------------------|
C
      DIMENSION UMONTH(IM,JM,KB),VMONTH(IM,JM,KB),SMONTH(IM,JM,KB)
      DIMENSION ADVUA(IM,JM),ADVVA(IM,JM),ADVUU(IM,JM),ADVVV(IM,JM)
      DIMENSION DRHOX(IM,JM,KB),DRHOY(IM,JM,KB),TRNU(IM,JM),TRNV(IM,JM)
      DIMENSION COM(80)
      DIMENSION IVAR(IM),PRT(IM,KB)
      EQUIVALENCE (IVAR,TPS),(PRT,A)
C
      CHARACTER*10 RESTAR,OPTION
      CHARACTER*40 PATHV
      CHARACTER*33 inpath
      CHARACTER*34 gpath,rpath
      CHARACTER*35 outpath
      character*55 gfile
      character*80 trackfile,root
      character*4 ysub  ! standard
      character*1 ab(2),case
      character*6 sfile
      character*19 cpcomm

C     NEW VARIABLES INVOLVED IN RUNNING-MEAN CALC OF SALINITY
      REAL SAVGTM   ! TIME OVER WHICH TO CALC AVERAGE OF SALINITY (HRS)
      REAL SAVG(IM,JM,KB),MONAVGSAL(IM,JM,KB)
      INTEGER ISFLAG,MONFLAG,MCNT,OPEN4,INTRSRT,negsal,igoodsal
      integer iab
      real fint,dint

c     must be consistent with genmain.f
      data ispan/0,0,0,0,1,0/
      data ncutx/0,0,0,0,2,0/
      data is/8,9,10,11,12,13,14,15,16/
      data js/33,32,32,32,32,32,32,32,32/
      data idxc/5,5,5,5,5,5,5,5,5/
      data ncellspan/0,0,0,0,9,0/

C 
      IURUN=5
      IUPRT=7
      IUGRD=8
      IUHYD=9
      IUPLT=10
      IUWRS=13
      IURRS=14
      IUFLW=15
      IUCCCRES=16
      IUTSR=17
c     full variables output files
      IUSAL1=80   ! SALINITY OUTPUT FILE
      IUVEL1=81  ! VELOCITY 
      IVVEL1=82  ! VELOCITY 
      IWVEL1=83  ! VELOCITY 
      IUELEV1=86
      IUTEMP1=87
      IUHVIS1=88
      IUVVIS1=89
      IUCCC = 79
      IEF = 70  ! LUN for control flags  Feb 2014  S.D. Meyers

      IUSAL=60   ! SPECIALIZED SALINITY OUTPUT FILE
      IUSAL2=61  ! 25 HR SALINITY FIELDS
      IUSAL3=63  ! MONTHLY SALINITY AVERAGE FIELD
      IUSAL4=64  ! SALINITY AT A FEW POINTS

      IUS1 = 61
      IUS2 = 62
      IUS3 = 63
      IUS4 = 64

      IELEV = 75

      IUT90=90
      IUT91=91
      IUT92=92
      IUT93=93
      IUT94=94

      IBH=48  ! model grid x-index to be output for BH input
      JBH=28  ! model grid y-index to be output for BH input
      IMB=55
      JMB=66
      IPR=60  ! number of time steps (model minutes) between writing output
      INTRSRT = 60*24*10 ! number of time steps (model minutes) between restarts
Crealtime_
      iuseedhr=95
      iulasthr=96
      iujtime=97
C
      hdry = 0.25  !  m, set ~0 so not active in 2001-03 runs
      hmin = 0.3+0.34443  ! min grid depth MSL
      doffset = 0.  ! SLR offset

c     read control flags
      open(IEF,file='eflags',status='old')
      READ(IEF,11)  (COM(I),I=1,80)
      READ(IEF,*)  ibath,iwds,idist,ireg,doffset,ibuoy
      close(ief)
c      write(6,*)  'eflags ',ibath,iwds,idist,ireg,doffset,ibuoy

       inpath = '/home/meyers/Sci/ECOM/GA/PT/Code/'
       rpath = '/home/meyers/Sci/ECOM/GA/Runfiles/'
       gpath = '/home/meyers/Sci/ECOM/GA/PT/Grids/'
       outpath = '/home/meyers/Sci/ECOM/GA/PT/Output/'

c     ------------------------------------------
c     determine scenario # and string id (sid)
c     ------------------------------------------
       open(10,file=inpath//'otb_config.dat',status='old')
       read(10,*) sidfile
       read(10,*) root   ! scenario 
       close(10)
 122   format(a10)

       ilen = len(root)
       ichrom = sidfile(5:6) ! assumes only 2 digits, see 'genetic'       
       open(10,file=inpath//sidfile,status='old')
       read(10,*) sid         ! gridid code, from mkrandgrid.f
       read(10,55) igen,ttest,tmin  ! generation #, (d) for ftest mkrandgrid.f
       close(10)
       if (igen.lt.10) then 
         write(sgen,'(i1)') igen
         sgen = "0"//sgen
       endif
       if ((igen.ge.10).and.(igen.lt.100)) then 
         write(sgen,'(i2)') igen
       endif
       open(11,file=inpath//sidfile(1:6)//'_'//sgen,status='replace')
       write(11,*) sid
       write(11,*) igen,ttest,tmin
       close(11)

       if (tmin.ge.ttest) then 
         write(6,*) 'ERROR tmin>ttest ',sid
         return
       endif
c       write(6,*) 'starting ',sid(1:nc),root(1:ilen)
 55    format(i4,2x,f10.2,2x,f10.2)
       

c     ------------------------------------------
c       read boundary file
c     ------------------------------------------
       OPEN (IURUN,FILE=rpath//'run_data2',status='old')

       if (ibath.eq.1) then 
          gfile = gpath//'otb_grid0.txt'
          iglen = 55
       endif
       if (ibath.eq.2) then 
          gfile = ''
          iglen = 1
       endif
       if (ibath.eq.3) then 
          gfile = ''
          iglen = 1
       endif
       open(iugrd,file = gfile(1:iglen),status='old')
c       write(6,*) 'Using grid file: ',gfile(1:iglen)

      OPEN (IUPRT,FILE='Output/gcmprt'//ichrom)
      OPEN (IUPLT,FILE='Output/gcmplt'//ichrom)
      OPEN (IUTSR,FILE='Output/gcmtsr'//ichrom)

      OPEN (IUFLW,FILE='gcm_flow'//ichrom,FORM='unformatted')
      OPEN (IUCCCRES,FILE='gcm_ccc'//ichrom,FORM='formatted')

      OPEN (IUT90,FILE='gcm_temp90_'//ichrom)
      OPEN (IUT91,FILE='gcm_temp91_'//ichrom)
      OPEN (IUT92,FILE='gcm_temp92_'//ichrom)
      OPEN (IUT93,FILE='gcm_temp93_'//ichrom)
      OPEN (IUT94,FILE='gcm_temp94_'//ichrom)
Crealtime
c      open (iuseedhr,file='seedhour.dat')
c      open (iulasthr,file='lasthour.dat')
      open (iujtime,file='jtime.dat')
      read (iujtime,603) rjtime 
c .....................................................     
C ... rjtime is the last time which data are prepared
c .....................................................       
C++++ 6-12-97 msv added option to open individual files for time series
C++++ open these in TINE ARCHIVE; reserve file numbers 20-40 for these
C++++ hardwired option to use ind. point files
c      ifpointfile=1
C++++ ifopenede and v will be reset to 1 after files are opened in SUB ARCHIVE
c      ifopenede=0
c      ifopenedv=0
C++++ 6-13-97 added path for individual tseries files
c      datapath='./ts_files/'      

       seedhour = 0.      
C
C++++ 8-16-97 aaded two new lines to run_data for depth scaling
c      write(6,*) "========= CONTROL PARAMS ",ibath,iwds,idist

      READ(IURUN,11)  (COM(I),I=1,80)
      READ(IURUN,99)  depthLB,depthMB,depthHB,depthOT
 99   FORMAT (4(f8.3))
C++++
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
c 20   format(a9)   ! must match size of sidfile in 'genetic'
 11   FORMAT(80A1)
 12   FORMAT(/1X,80A1/)
C
      READ(IURUN,11)  (COM(I),I=1,80)
      READ(IURUN,1) DTI,ISPLIT,IRAMP
 1    FORMAT(1F10.4,3I5)
      dti0 = dti  ! reference time step
C
      DTE=DTI/FLOAT(ISPLIT)
      DTE2=2.0*DTE
      DTI2=2.0*DTI
      ISPI=1.0/FLOAT(ISPLIT)
      ISP2I=0.5*ISPI
      DAYI=1.0/86400.0
      GRAV=9.806
C 
      WRITE(IUPRT,21) DTI,DTE,ISPLIT,IRAMP
 21   FORMAT(
     . ' BAROCLINIC TIME STEP IS                 ',F10.4,' SECONDS',//,
     . ' BAROTROPIC TIME STEP IS                 ',F10.4,' SECONDS',//,
     . ' INTERNAL/EXTERNAL MODE SPLITTING IS     ',I10//,
     . ' NUMBER OF RAMP TIME STEPS               ',I10//)
C
C-----------------------------------------------------------------------
C      TYPE OF RUN -
C      BAROTROPIC: 2-D CALCULATION (BOTTOM STRESS CALCULATED IN ADVAVE)
C      PROGNOSTIC: 3-D CALCULATION (BOTTOM STRESS CALCULATED IN PROFU,V)
C      DIAGNOSTIC: 3-D CALCULATION WITH T AND S HELD FIXED
C----------------------------------------------------------------------- 
C         3-D - TYPE OF MOMENTUM ADVECTION AND BOTTOM FRICTION
C      LINEAR    : ALL MOMENTUM ADVECTION NEGLECTED
C      NON-LINEAR: COMPLETE PHYSICS 
C----------------------------------------------------------------------  
C
      READ(IURUN,11)  (COM(I),I=1,80)
      READ(IURUN,10) NSTEPS,IPRINT,RESTAR,TOR,ADVECT
c      write(6,10) NSTEPS,IPRINT,RESTAR,TOR,ADVECT
      READ(IURUN,11)  (COM(I),I=1,80)
      READ(IURUN,3) BFRIC,Z0B,NU,ALPHA,TLAG
      READ(IURUN,11)  (COM(I),I=1,80)
      READ(IURUN,23) HORZMIX,HORCON,HPRNU
      horcon0=horcon  ! for adaptive time stepping
      READ(IURUN,11)  (COM(I),I=1,80)
      READ(IURUN,23) VERTMIX,UMOL,VPRNU
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,7) JHM,IAVGE
      WRITE(IUPRT,7) JHM,IAVGE 
      
C      WRITE(*,*) 'JHM,IAVGE',JHM,IAVGE
      IF(JHM.EQ.0 ) GO TO 200
      DEI=1.0/FLOAT(IAVGE)
      IF (MOD(JHM,10).EQ.0) THEN
      NJHM=JHM/10
      ELSE
      NJHM=JHM/10+1
      END IF    
      DO 17 NJH=1,NJHM
      NJHB=(NJH-1)*10+1
      NJHE=MIN(JHM,10*NJH)            
      READ(IURUN,6) (IHIST(I,2),I=NJHB,NJHE)      
      WRITE(IUPRT,6)(IHIST(I,2),I=NJHB,NJHE)  
c      WRITE(*,6)(IHIST(I,2),I=NJHB,NJHE)      
   17 CONTINUE      
   10 FORMAT(2I10,1X,A10,1X,A10,1X,A10)
   23 FORMAT(A10,2E10.3)
    3 FORMAT(5E10.3)
    5 FORMAT(16I5)
    6 FORMAT(10I8)
    7 FORMAT(8I10)
   49 FORMAT(4(2I5,1X,A4,I5))
  200 READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,7) ISKILL
      WRITE(IUPRT,7) ISKILL
      IF(ISKILL.EQ.0) THEN
       SKILLI=1.0
      ELSE
       SKILLI=1.0/FLOAT(ISKILL)
      ENDIF
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,5)  EPTS
      WRITE(IUPRT,5) EPTS
      IF(EPTS.EQ.0) GOTO 201
      READ(IURUN,5)  (INXIE(I),INXJE(I),I=1,EPTS)
      WRITE(IUPRT,5) (INXIE(I),INXJE(I),I=1,EPTS)
201   READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,5)  VPTS
      WRITE(IUPRT,5) VPTS
      IF(VPTS.EQ.0) GOTO 202
      READ(IURUN,5)  (INXIV(I),INXJV(I),I=1,VPTS)
      WRITE(IUPRT,5) (INXIV(I),INXJV(I),I=1,VPTS)
 202  READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,5)  FPTS
      WRITE(IUPRT,5) FPTS
      IF(FPTS.EQ.0) GOTO 203
      READ(IURUN,49) (ISFLX(N),JSFLX(N),DIRFLX(N),NFLXE(N),N=1,FPTS)
      WRITE(IUPRT,49) (ISFLX(N),JSFLX(N),DIRFLX(N),NFLXE(N),N=1,FPTS)
 203  READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,7) JTM,ITRACE
      WRITE(IUPRT,7) JTM,ITRACE
      IF(JTM.EQ.0 ) GO TO 204
      DEIT=1./FLOAT(ITRACE)
      READ(IURUN,6) (ITRAC(I,2),I=1,JTM)
      WRITE(IUPRT,6) (ITRAC(I,2),I=1,JTM)
 204  READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
c      WRITE(6,12) (COM(I),I=1,80)
      READ(IURUN,*) nccctime
      IF(nccctime.EQ.0 ) GO TO 205
      READ(IURUN,*) (iccctime(i),i=1,nccctime)
      do 199 i=1,nccctime
         iccctime(i) = iccctime(i)*60./dti
c         write(6,*) 'iccc',i,iccctime(i)
 199  enddo
      WRITE(IUPRT,*) (iccctime(i),i=1,nccctime)
c      WRITE(6,*) 'release times ',(iccctime(i),i=1,nccctime)
      do i=1,nccctime
        DEITCCC(i)=1./FLOAT(iccctime(i))
      enddo
c      READ(IURUN,*) iccctime
c      IF(iccctime.EQ.0 ) GO TO 205
c      WRITE(IUPRT,*) iccctime    
c      WRITE(6,*) 'iccctime=',iccctime    
c      DEITCCC=1./FLOAT(iccctime)
 205  READ(IURUN,11)  (COM(I),I=1,80)
      READ(IURUN,*)  SAVGTM
c      write(6,*) 'savgtm=',savgtm
 206  CONTINUE
C
      IF(TOR.NE.'BAROTROPIC'.AND.TOR.NE.'PROGNOSTIC'.AND.
     . TOR.NE.'DIAGNOSTIC') THEN
      WRITE(IUPRT,25) TOR
 25   FORMAT(//' TYPE OF RUN (TOR=',A10,') IS SPECIFIED INCORRECTLY! '
     .    ,'FIX AND RESUBMIT'//)
      STOP
      END IF
C
      IF (ADVECT.NE.'LINEAR    ' .AND.ADVECT.NE.'NON-LINEAR') THEN
      WRITE(IUPRT,26) ADVECT
  26  FORMAT(//'  TYPE OF ADVECTION (ADVECT=',A10,') IS SPECIFIED ',
     1 'INCORRECTLY!    FIX AND RESUBMIT'//)
      STOP
      END IF
C
      IF (HORZMIX.NE.'CLOSURE   ' .AND.HORZMIX.NE.'CONSTANT  ') THEN
      WRITE(IUPRT,27) HORZMIX
  27  FORMAT(//'  TYPE OF HORIZONTAL MIXING (HORZMIX=',A10,
     .  ') IS SPECIFIED INCORRECTLY!    FIX AND RESUBMIT'//)
      STOP
      END IF
C
      IF (VERTMIX.NE.'CLOSURE   ' .AND.VERTMIX.NE.'CONSTANT  ') THEN
      WRITE(IUPRT,28) VERTMIX
  28  FORMAT(//'  TYPE OF VERTICAL MIXING (VERTMIX=',A10,
     .  ') IS SPECIFIED INCORRECTLY!    FIX AND RESUBMIT'//)
      STOP
      END IF
C
      IF(TOR.EQ.'BAROTROPIC') THEN
      WRITE(IUPRT,14) TOR
      ELSE
      WRITE(IUPRT,13) TOR
      END IF
 13   FORMAT(/' THIS IS A THREE DIMENSIONAL MODEL RUN',2X,A10/)
 14   FORMAT(/' THIS IS A TWO DIMENSIONAL MODEL RUN',2X,A10/)
C
      WRITE(IUPRT,22) ADVECT
 22   FORMAT(/' THIS SIMULATION HAS ',A10,' MOMENTUN ADVECTION '/)
C
      IF(HORZMIX.EQ.'CLOSURE   ') THEN
      WRITE(IUPRT,29) HORZMIX,HORCON,HPRNU
      ELSE
      WRITE(IUPRT,31) HORZMIX,HORCON,HPRNU
      END IF
 29   FORMAT(/' THIS SIMULATION HAS ',A10,' HORIZONTAL MIXING ',
     . ' HORCON = ',1PE10.3,' HPRNU = ',1PE10.3/)
 31   FORMAT(/' THIS SIMULATION HAS ',A10,' HORIZONTAL MIXING ',
     . ' CONSTANT = ',1PE10.3,'m**2/s  HPRNU = ',1PE10.3/)
C
      IF(VERTMIX.EQ.'CLOSURE   ') THEN
      WRITE(IUPRT,32) VERTMIX,UMOL,VPRNU
      ELSE
      WRITE(IUPRT,33) VERTMIX,UMOL,VPRNU
      END IF
 32   FORMAT(/' THIS SIMULATION HAS ',A10,' VERTICAL MIXING ',
     . ' UMOL = ',1PE10.3,' HPRNU = ',1PE10.3/)
 33   FORMAT(/' THIS SIMULATION HAS ',A10,' VERTICAL MIXING ',
     . ' CONSTANT = ',1PE10.3,'m**2/s  HPRNU = ',1PE10.3/)
C
      CALL ZEROES(ADVUA,ADVVA,ADVUU,ADVVV,DRHOX,DRHOY,TRNU,TRNV)
      DO 46 I=1,IM
      DO 46 J=1,JM
      DO 46 K=1,KB
      UMONTH(I,J,K)=0.0
      VMONTH(I,J,K)=0.0
      SMONTH(I,J,K)=0.0

c     initialization 1/1/2001 obc -0.226
c     initialization 1/1/2004 obc +0.1925
      e0 = 0.
      EL(i,j) = e0
      ELF(I,J)= e0
      ELB(I,J)= e0
   46 CONTINUE
C
      
      IF(RESTAR.EQ.'COLD START') THEN
c         write(6,*) 'INITIATING COLD START'
         CALL SETDOM(ADVUA,ADVVA,ADVUU,ADVVV,DRHOX,DRHOY,TRNU,TRNV)
         CALL BCDATA2
cc      write(6,*) 'Hd 54 27 ',h(54,27),s(54,27,1)
      CALL TANDS(0)  ! 0-raw start  1-input file
      CALL DENS

Crealtime
c      read(iuseedhr,603) seedhour
      int=1   !seedhour * 60 * 60 * 1/DTI
      intcold=int
c      write(6,*) 'INITIAL S',S(36,16,1),s(36,17,1),s(37,15,1)

      ELSE
C
Crealtime
c         write(6,*) 'OPENING RESTART FILE'
         OPEN (IURRS,FORM='UNFORMATTED',FILE='restart_'//sid)
         READ (IURRS) 
     .     intcold,INT,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,AAM,AAM2D,
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,WUBOT,WVBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN
c     .     CCC,CCCB,CCCMEAN  ! commented for nowcast restart
c      int=intcold
c         write(6,*) 'dvm',dvm(47,37),dvm(48,37),dvm(48,38)
c         write(6,*) 'v',v(47,37,1),v(48,37,1),v(48,38,1)
c         write(6,*) 'restarting at int=',int,int/60.
c         do i=10,60
c         do j=10,90
c           write(6,*) i,j,el(i,j)
c         enddo
c         enddo
c         write(6,*) 'starting boundary el=',el(18,52),elb(18,52),
c     .        u(18,52,1),ub(18,52,1),v(18,52,1),vb(18,52,1)
      CLOSE (IURRS)
      CALL BCDATA2
      TIME=FLOAT(INT)*DAYI*DTI
C
      END IF
C
C -----  Temperature is disconnected from density and is used
C -----  as a passive scalar, 2-D versus 3-D; source is introduced
C -----  through the boundary conditions
C
cdcb      DO 34 K=1,KB
cdcb      DO 34 J=1,JM
cdcb      DO 34 I=1,IM
C SET INITIAL BAY TEMPERATURE
cdcb      T(I,J,K)=22.000
cdcb      TB(I,J,K)=22.0     
cdcb   34 CONTINUE
C
      DO 45 J=1,JM
      DO 45 I=1,IM
      Z0=Z0B
      CBCMIN=BFRIC
  45  CBC(I,J)=AMAX1(CBCMIN,.16/ALOG((ZZ(KBM1)-Z(KB))*H(I,J)/Z0)**2)*
     .    FSM(I,J)
C
      CALL PRINTS(DRHOX,DRHOY,TRNU,TRNV)
C
      ISTART=INT+1
      TNDAYS=FLOAT(NSTEPS)*DAYI*DTI
c      WRITE(6,*) 'TNDAYS=',TNDAYS
      IEND=INT+NSTEPS
      TPRT=FLOAT(IPRINT)*DAYI*DTI
      TAVG=FLOAT(IAVGE)*DAYI*DTI
      TSKILL=FLOAT(ISKILL)*DAYI*DTI
      TRACE=FLOAT(ITRACE)*DAYI*DTI
      DO 16 I=1,JHM
 16   IHIST(I,1)=IHIST(I,2)-IAVGE+1
      DO 18 I=1,JTM
 18    ITRAC(I,1)=ITRAC(I,2)-ITRACE+1
C 
      WRITE(IUPRT,15) ISTART,IEND 
c      WRITE(6,*)'ISTART,IEND =',ISTART,IEND 
 15   FORMAT(//30H MODEL STARTING UP...ISTART = ,I6,8H IEND = ,I6/)
      WRITE(IUPRT,20) TNDAYS
 20   FORMAT(//32H NUMBER OF DAYS IN SIMULATION = ,F6.2/) 
      WRITE(IUPRT,30) TPRT,IPRINT,TAVG,IAVGE,TSKILL,ISKILL
 30   FORMAT(//' TPRT =   ',F10.3,'  IPRINT =  ',I10,//
     3         ' TAVG =   ',F10.3,'  IAVGE =   ',I10,//
     4         ' TSKILL = ',F10.3,'  ISKILL =  ',I10//)
      WRITE(IUPRT,35)
      WRITE(IUPRT,40) (IHIST(I,1),IHIST(I,2),I=1,JHM)
 35   FORMAT(//' HISTORY TAKEN AT TIMESTEPS  START---STOP ')
 40   FORMAT(28X,I5,2X,I5)
      WRITE(IUPRT,47)
      WRITE(IUPRT,48) (ITRAC(I,1),ITRAC(I,2),I=1,JTM)
47    FORMAT(//' QUALITY PARAMETERS FOR AESOP INTEGRATED OVER TIMESTEPS
     .      START      STOP')
48    FORMAT(57X,I8,2X,I8)
C
      WRITE(IUPRT,4) BFRIC,Z0B,NU,ALPHA,TLAG
 4    FORMAT(//' BFRIC            =   ',F10.4,' nondimensional'/
     .         ' Z0B              =   ',F10.4,' m'/
     .         ' NU               =   ',F10.4,' nondimensional'/
     .         ' ALPHA            =   ',F10.4,' nondimensional'/
     .         ' INFLOW LAG TIME  =   ',F10.0,  ' seconds'/)
C
      TLAG=1./TLAG
      WRITE(IUPRT,70)
  70  FORMAT(/1X,' K',6X,'Z',10X,'ZZ',8X,'DZ',/)
      DO 90 K=1,KBM1
      WRITE(IUPRT,80) K,Z(K),ZZ(K),DZ(K)
  90  CONTINUE
      WRITE(IUPRT,80) K,Z(KB)
  80  FORMAT(I3,3F10.3)
C
      THOUR=FLOAT(INT)*DTI/3600.
c      write(6,*) 'Hr 54 27 ',h(54,27),s(54,27,1)
      CALL FIRST
c      write(6,*) 'Hs 54 27 ',h(54,27),s(54,27,1)
C 
      AREA=0.0
      VOLUME=0.0
      DO 280 J=1,JM
      DO 280 I=1,IM
      AREA=AREA+FSM(I,J)*H1(I,J)*H2(I,J)
      VOLUME=VOLUME+H(I,J)*FSM(I,J)*H1(I,J)*H2(I,J)
  280 CONTINUE
C
C-----------------------------------------------------------------------
C
C                 BEGIN NUMERICAL INTEGRATION
C
C-----------------------------------------------------------------------
C
      OPEN4=0
      MONFLAG=0
      ISFLAG=0
CDCB initialize the temp loop
      NOFI = 0

c----------------------------------------------
c     setup files for sampout
c----------------------------------------------


      ISMP(1) = 13
      JSMP(1) = 36  

      ISMP(2) = 55
      JSMP(2) = 71

      ISMP(3) = 54
      JSMP(3) = 58

      ISMP(4) = 55
      JSMP(4) = 66  

      ISMP(5) = 23
      JSMP(5) = 70

      ISMP(6) = 54
      JSMP(6) = 54

      ISMP(7) = 54
      JSMP(7) = 70

      ISMP(8) = 28
      JSMP(8) = 39

      ISMP(9) = 48
      JSMP(9) = 73 
      ISMP(10) = 48
      JSMP(10) = 72    

      ISMP(11) = 61
      JSMP(11) = 56    

      ISMP(12) = 63
      JSMP(12) = 73

c      write(ius1,*) nsmpnts
c      write(ius2,*) nsmpnts
c      write(ius3,*) nsmpnts
c      write(ius4,*) nsmpnts
c      do i=1,nsmpnts
c         write(ius1,*) ISMP(I),JSMP(I)
cc         write(ius2,*) ISMP(I),JSMP(I)
cc         write(ius3,*) ISMP(I),JSMP(I)
c         write(ius4,*) ISMP(I),JSMP(I)
c      enddo

c----------------------------------------------
Crealtime This checks to see if run out of bc data from rtmon
c----------------------------------------------
      if (IEND/60..gt.rjtime) then
       print *,'exceeded bc time',IEND/60.,rjtime
       goto 999
      endif

      iwrite=0
      ido_wd = 0  ! wet-dry flag
      ido_adj = 0 ! adaptive time step flag

c************************************************
c     Begin time loop
c************************************************
      cflmin = 5.
      cflmin0 = 5
      nadj = 1                  ! # time step
      dint = 1.                 ! increment per step for INT
      RAMP=TANH(FLOAT(INT-intcold)/FLOAT(IRAMP+1))     
      TIME=FLOAT(INT)*DAYI*DTI
      THOUR=FLOAT(INT)*DTI/3600.

      fint = int
 9000 continue
      if (int.gt.iend) goto 9001
c      DO 9000 INT=ISTART,IEND
c      if (mod(int,1000).eq.0)write(6,*)
c      if (mod(int,10000).lt.dint)write(6,*)
c     & 'int=',fint,iend,cflmin,cflmin0,time,thour,dti,horcon
c         write(6,*) 'int=',int,iend
         smin = 9999.
         smax = -9999.

c---------------------------------
c  determine if new half-year
c---------------------------------
      iopen = 0
      iclose=0
      ab(1) = 'a'
      ab(2) = 'b'

c --------
c   output files
c --------
      if (int.eq.1) then   ! open first years
c      if (abs(int-1).lt.dint) then   ! open first years
         iopen = 1
         iclose = 0
         ysub = 'test'
c         ysub = '1990'
         iab = 1
         iwrite=1
      endif
c      if (int.eq.262800*60./dti) then   
c      if (abs(int-262800*60./dti).lt.dint) then   
c         iopen = 1
c         iclose = 1
c         ysub = '2001'
cc         ysub = '1990'
c         iab = 2
c         iwrite=1
c      endif


      if (iclose.eq.1) then 
	CLOSE(IUSAL1)
	CLOSE(IUVEL1)
	CLOSE(IVVEL1)
	CLOSE(IWVEL1)
	CLOSE(IUELEV1)
	CLOSE(IUTEMP1)
	CLOSE(IUHVIS1)
	CLOSE(IUVVIS1)
      endif  ! iclose


c  release dye as ccc>1
c      if (int-iccctime(1).eq.0) then 
c         iopen=1
c         iwrite=1
c      endif
c      if (iopen.eq.0) iwrite=0

c     no output of prognostic variables
c      iwrite=0
c      iopen=0

c     each year split in half bc of memory limitations
      if (iopen.eq.1) then 
         iwrite=1  ! write at this time
         iwritedo=1  ! enable writing
c         write(6,*) 'AT INT=',int 
c      write(6,*) 'opening ',outpath//' file for '//
c     .        ysub//ab(iab)//root(1:ilen)//'.dat'
c      OPEN (IUSAL1,FILE=outpath//'salinity_'//
c     .     ysub//ab(iab)//root(1:ilen)//ichrom//'.dat',
c     & FORM='unformatted')
      OPEN (IUVEL1,FILE=outpath//'u_velocity_'//
     .     ysub//ab(iab)//root(1:ilen)//ichrom//'.dat',
     & FORM='unformatted')
      OPEN (IVVEL1,FILE=outpath//'v_velocity_'//
     .     ysub//ab(iab)//root(1:ilen)//ichrom//'.dat',
     & FORM='unformatted')
      OPEN (IWVEL1,FILE=outpath//'w_velocity_'//
     .     ysub//ab(iab)//root(1:ilen)//ichrom//'.dat',
     & FORM='unformatted')
      OPEN (IUELEV1,FILE=outpath//'elevation_'//
     .     ysub//ab(iab)//root(1:ilen)//ichrom//'.dat',
     & FORM='unformatted')
c      OPEN (IUTEMP1,FILE=outpath//'temperature_'//
c     .     ysub//ab(iab)//root(1:ilen)//'.dat',
c     & FORM='unformatted')
c      OPEN (IUHVIS1,FILE=outpath//'horvisc_'//
c     .     ysub//ab(iab)//root(1:ilen)//ichrom//'.dat',
c     & FORM='unformatted')
c      OPEN (IUVVIS1,FILE=outpath//'vertvisc_'//
c     .     ysub//ab(iab)//root(1:ilen)//ichrom//'.dat',
c     & FORM='unformatted')
c      OPEN (IUCCC,FILE=outpath//'ccc_'//root(1:ilen)//
c     &     '.dat',FORM='unformatted')
c      write(iuccc) time,ccc  ! write initial configuration
c
      iopen=0
      ENDIF  ! iopen

Crealtime

c      RAMP=TANH(FLOAT(INT-intcold)/FLOAT(IRAMP+1))     
c      TIME=FLOAT(INT)*DAYI*DTI
c      THOUR=FLOAT(INT)*DTI/3600.

c--------------------------------------------------------
c     check stability relative to CFL
c     adjust time step if needed
c--------------------------------------------------------
c      if (ido_adj.eq.1) then
         cflmin0 = cflmin
cc         write(6,*) 'checking cfl ',dint,dti
c         call testcfl(cflmin)   ! check CFL risk
c         if (int.gt.istart+iramp) then 
c            call tadjust(nadj,dint,cflmin,cflmin0) ! adjust time step if needed
c            DTE=DTI/FLOAT(ISPLIT)
c            DTE2=2.0*DTE
c            DTI2=2.0*DTI
cc         if (dti.ge.dti0)horcon = horcon0*(dti0/dti)**1
cc         write(6,*) 'DTI test ',dti,dti2,dte
cc         cflmin0 = cflmin
c         endif
c      endif  ! ido_adj

c--------------------------------------------------------
c     begin adjustable time loop
c     must be 2^n step
c--------------------------------------------------------

      do 6000 iadj=1,nadj
         RAMP=TANH((fINT-FLOAT(intcold))/FLOAT(IRAMP+1))     
         TIME=fINT*DAYI*DTI0
         THOUR=fINT*DTI0/3600.


CMSV
c      write(6,*)'CALL BCOND2, RAMP=',RAMP,fINT,intcold,IRAMP
      CALL BCOND2(7,DTI2)
C
      IF(TOR.NE.'BAROTROPIC') THEN
C
      CALL BAROPG(DRHOX,DRHOY,TRNU,TRNV)
C
      DO 50 J=1,JM
      DO 50 I=1,IM
      TRNU(I,J)=TRNU(I,J)+ADVUU(I,J)-ADVUA(I,J)
 50   TRNV(I,J)=TRNV(I,J)+ADVVV(I,J)-ADVVA(I,J)
C
      DO 120 J=1,JM
      DO 120 I=1,IM
 120  EGF(I,J)=EL(I,J)*ISPI   !****
      DO 400 J=2,JM
      DO 400 I=2,IM
      UTF(I,J)=UA(I,J)*(D(I,J)+D(I-1,J))*ISP2I
      VTF(I,J)=VA(I,J)*(D(I,J)+D(I,J-1))*ISP2I  
400   CONTINUE
C
      ENDIF
C
      IF(HORZMIX.EQ.'CLOSURE   ') CALL SMAG
C
C-----------------------------------------------------------------------
C         BEGIN EXTERNAL MODE
C-----------------------------------------------------------------------
C
                   DO 8000 IEXT=1,ISPLIT
c                      write(6,*) ' el(18,52) = ',iext,el(18,52)
C
c                      write(6,*) 'call extrnl ',iext,int
      CALL EXTRNL(ADVUA,ADVVA,TRNU,TRNV,DTE2,BFRIC,DTI2)
C
      IF(TOR.EQ.'BAROTROPIC') GO TO 440
      IF(IEXT.LT.(ISPLIT-2)) GO TO 440
      IF(IEXT.EQ.(ISPLIT-2)) THEN
      DO 402 J=1,JM
      DO 402 I=1,IM
 402  ETF(I,J)=0.25*NU*ELF(I,J)
      GO TO 440
      ENDIF
      IF(IEXT.EQ.(ISPLIT-1)) THEN
      DO 404 J=1,JM
         DO 404 I=1,IM 
 404  ETF(I,J)=ETF(I,J)+.5*(1.-.5*NU)*ELF(I,J)
      GO TO 440
      ENDIF
      IF(IEXT.EQ.(ISPLIT-0)) THEN
      DO 406 J=1,JM
      DO 406 I=1,IM
 406  ETF(I,J)=(ETF(I,J)+.5*ELF(I,J))*FSM(I,J)
      ENDIF

 440  CONTINUE

C-----------------------------------------------------------------------
C         APPLY FILTER TO REMOVE TIME SPLIT AND RESET TIME SEQUENCE
C-----------------------------------------------------------------------
C
      DO 150 J=1,JM 
      DO 150 I=1,IM 
      UA(I,J)=UA(I,J)+.5*NU*(UAB(I,J)-2.*UA(I,J)+UAF(I,J)) !****
      VA(I,J)=VA(I,J)+.5*NU*(VAB(I,J)-2.*VA(I,J)+VAF(I,J)) 
      EL(I,J)=EL(I,J)+.5*NU*(ELB(I,J)-2.*EL(I,J)+ELF(I,J)) 
c      if ((s(i,j,1).le.0).and.(fsm(i,j).eq.1))write(6,*)'s',
c     &     int,s(i,j,1),i,j,h(i,j),el(i,j)
c      dij = h(i,j)+el(i,j)
c      if ((fsm(i,j).eq.1).and.(dij.le.hdry)) then 
c         el(i,j) = -h(i,j)+hdry-0.0000001  ! extra is for roundoff
c         wd(i,j) = 0                    
c      endif
  150  continue

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
               write(6,*) 'DRY at ',i,j
               if (uaf(i,j).lt.0)  uaf(i,j)=0
               if (uaf(i+1,j).gt.0) uaf(i+1,j)=0
               if (vaf(i,j).lt.0) vaf(i,j)=0
               if (vaf(i,j+1).gt.0) vaf(i,j+1)=0
            endif       
 408     continue
c         if (ndry.gt.0) write(6,*)'ndry=',ndry,int*dti/3600,' hr'
      endif                     !ido

c      write(6,*) int,'el(18,52)',el(48,37)
       igoodsal = 0  ! flag for all s=0
       negsal = 0   ! flag for negative salinity
      DO 160 J=1,JM 
      DO 160 I=1,IM 
      ELB(I,J)=EL(I,J)
      EL(I,J)=ELF(I,J)
      D(I,J)=H(I,J)+EL(I,J) 
      if ((fsm(i,j)*s(i,j,1)).gt.0.) then 
         igoodsal = 1 
      endif
      IF ((fsm(i,j).eq.1).and.(D(I,J).LE.0.0)) THEN 
      WRITE(6,161) IEXT,I,J,H(I,J),EL(I,J),elb(i,j),d(i,j)
      WRITE(IUPRT,161) IEXT,I,J,H(I,J),EL(I,J),elb(i,j),d(i,j)
      write(6,*) 'final int=',INT,int*dti/3600.
      write(IUPRT,*) 'final int=',INT,int*dti/3600.
  161 FORMAT(//3X,'****************************' 
     1        /3X,'** NEGATIVE D **************' 
     2        /3X,'****************************' 
     3       //3X,'IEXT=',I10 
     4        /3X,'I=',I5,'  J=',I5,'  H=',F10.4,'  EL=',3(1x,F10.4))  
      write(6,*) 'ndry=''',ndry
      write(iuelev1) time,el
c      write(iusal1) time,s
c      write(iutemp1) time,t
c      write(iuvel1) time,u
c      write(ivvel1) time,v
      CALL PRINTS(DRHOX,DRHOY,TRNU,TRNV)
      STOP
      END IF
      UAB(I,J)=UA(I,J)
      UA(I,J)=UAF(I,J)  ! ****
      VAB(I,J)=VA(I,J)
 160  VA(I,J)=VAF(I,J)
C
      if (negsal.eq.1) then
      WRITE(6,162) negsal,igoodsal
      WRITE(IUPRT,162) negsal,igoodsal
 162  FORMAT(//3X,'****************************' 
     1        /3X,'** BAD SALINITY  **************' 
     2        /3X,'****************************'
     3        /3x,'negsal=',i2,'  igoodsal=',i2)
      write(6,*) 'final int=',INT,int/60.
      write(IUPRT,*) 'final int=',INT,int*dti/3600., ' hr'
      STOP
      endif
      



      IF(TOR.EQ.'BAROTROPIC') GO TO 8000
      IF(IEXT.EQ.ISPLIT) GO TO 8000
      DO 445 J=1,JM
      DO 445 I=1,IM
 445  EGF(I,J)=EGF(I,J)+EL(I,J)*ISPI
      DO 450 J=2,JM
      DO 450 I=2,IM
      UTF(I,J)=UTF(I,J)+UA(I,J)*(D(I,J)+D(I-1,J))*ISP2I
 450  VTF(I,J)=VTF(I,J)+VA(I,J)*(D(I,J)+D(I,J-1))*ISP2I
C
 8000              CONTINUE
c                   write(6,*) 'next el(32,18) = ',el(32,18)
C
C-----------------------------------------------------------------------
C         END EXTERNAL (2-D) MODE CALCULATION
C      AND CONTINUE WITH INTERNAL (3-D) MODE CALCULATION
C-----------------------------------------------------------------------
C 
Crealtime
      IF((INT-intcold).EQ.1) GO TO 8200
Crealtime      IF(INT.EQ.1) GO TO 8200 
      IF(TOR.EQ.'BAROTROPIC') GO TO 8300
C
C-----------------------------------------------------------------------
C         ADJUST U(Z) AND V(Z) SUCH THAT
C      VERTICAL AVERAGE OF (U,V) = (UA,VA)
C-----------------------------------------------------------------------
C
      DO 299 J=1,JM
      DO 299 I=1,IM
 299  TPS(I,J)=0.0
      DO 300 K=1,KBM1
      DO 300 J=2,JM
      DO 300 I=2,IM
 300  TPS(I,J)=TPS(I,J)+U(I,J,K)*DZ(K)
      DO 302 K=1,KBM1
      DO 302 J=2,JM
      DO 302 I=2,IM
 302  U(I,J,K)=(U(I,J,K)-TPS(I,J))+
     .     (UTB(I,J)+UTF(I,J))/(DT(I,J)+DT(I-1,J))
      DO 303 J=1,JM
      DO 303 I=1,IM
 303  TPS(I,J)=0.0
      DO 304 K=1,KBM1
      DO 304 J=2,JM
      DO 304 I=2,IM
 304  TPS(I,J)=TPS(I,J)+V(I,J,K)*DZ(K)
      DO 306 K=1,KBM1
      DO 306 J=2,JM
      DO 306 I=2,IM
 306  V(I,J,K)=(V(I,J,K)-TPS(I,J))+
     .     (VTB(I,J)+VTF(I,J))/(DT(I,J)+DT(I,J-1))
C
C-----------------------------------------------------------------------
C         CALCULATE HORIZONTAL MASS FLUXES, (H2*U*D) AND (H1*V*D)
C-----------------------------------------------------------------------
C
      DO 310 K=1,KBM1
      DO 311 J=2,JMM1
      DO 311 I=2,IM
      XMFL3D(I,J,K)=0.25*(H2(I-1,J)+H2(I,J))*(DT(I-1,J)+DT(I,J))*
     .   U(I,J,K)
  311 CONTINUE
      DO 312 J=2,JM
      DO 312 I=2,IMM1
      YMFL3D(I,J,K)=0.25*(H1(I,J-1)+H1(I,J))*(DT(I,J-1)+DT(I,J))*
     .   V(I,J,K)
  312 CONTINUE
  310 CONTINUE
C
C-----------------------------------------------------------------------
C         VERTVL INPUT = U,V,DT(=H+ET),ETF,ETB; OUTPUT = W
C---------------------------------------------------------------------- - 
C
      CALL VERTVL(DTI2)
      CALL BCOND2(5,DTI2)
c      CALL ARCHIVE(1)
C
C-----------------------------------------------------------------------
C         COMPUTE Q2F AND Q2LF USING UF AND VF AS TEMPORARY VARIABLES
C-----------------------------------------------------------------------
C
      IF(VERTMIX.EQ.'CLOSURE   ') THEN
      CALL ADVQ(Q2B,Q2,DTI2,UF)
      CALL ADVQ(Q2LB,Q2L,DTI2,VF)
      CALL PROFQ(DTI2)
      CALL BCOND2(6,DTI2)
      DO 325 K=1,KB 
      DO 325 J=1,JM 
      DO 325 I=1,IM 
      Q2 (I,J,K)=Q2 (I,J,K)+.5*NU*(UF(I,J,K)+Q2B(I,J,K)-2.*Q2(I,J,K)) 
      Q2B(I,J,K)=Q2(I,J,K)
      Q2(I,J,K)=UF(I,J,K) 
 325  CONTINUE
      DO 335 K=1,KB 
      DO 335 J=1,JM 
      DO 335 I=1,IM 
      Q2L(I,J,K)=Q2L(I,J,K)+.5*NU*(VF(I,J,K)+Q2LB(I,J,K)-
     2  2.*Q2L(I,J,K))
      Q2LB(I,J,K)=Q2L(I,J,K)
      Q2L(I,J,K)=VF(I,J,K)
 335  CONTINUE
      END IF
C
C-----------------------------------------------------------------------
C         COMPUTE TF AND SF USING UF AND VF AS TEMPORARY VARIABLES
C-----------------------------------------------------------------------
C
      IF(TOR.EQ.'PROGNOSTIC') THEN
      CALL ADVT(TB,T,TMEAN,DTI2,UF,TDIFF)
      CALL ADVT(SB,S,SMEAN,DTI2,VF,SDIFF)
      CALL ADVT(CCCB,CCC,CCCMEAN,DTI2,CCCF,CCCDIFF)

      CALL PROFT(UF,WTSURF,DTI2)
      CALL PROFT(VF,WSSURF,DTI2)
      CALL PROFT(CCCF,WCCCSURF,DTI2)

      CALL BCOND2(4,DTI2)
      DO 345 K=1,KB 
      DO 345 J=1,JM 
      DO 345 I=1,IM 
c     wd change by S Meyers April 2013
c      if ((fsm(i,j).eq.1).and.(wd(i,j).eq.0)) goto 345  ! skip update
      T(I,J,K)=T(I,J,K)+.5*NU*(UF(I,J,K)+TB(I,J,K)-2.*T(I,J,K))
      TB(I,J,K)=T(I,J,K)
      T(I,J,K)=UF(I,J,K)
 345  CONTINUE

      DO 355 K=1,KB 
      DO 355 J=1,JM 
      DO 355 I=1,IM 
c     change by S Meyers April 2013
c      if ((fsm(i,j).eq.1).and.(wd(i,j).eq.0)) goto 355  ! skip update
      S(I,J,K)=S(I,J,K)+.5*NU*(VF(I,J,K)+SB(I,J,K)-2.*S(I,J,K))
      SB(I,J,K)=S(I,J,K)
      S(I,J,K)=VF(I,J,K)
c      if ((fsm(i,j).eq.1).and.(s(i,j,k).le.0)) then 
c         write(6,*)'BAD Salinity ',i,j,k,s(i,j,k),wd(i,j),h(i,j),el(i,j)
c      endif
 355  CONTINUE

      DO 365 K=1,KB 
      DO 365 J=1,JM 
      DO 365 I=1,IM 
c     change by S Meyers April 2013
c      if ((fsm(i,j).eq.1).and.(wd(i,j).eq.0)) goto 365  ! skip update
      CCC(I,J,K)=CCC(I,J,K)+
     &.5*NU*(CCCF(I,J,K)+CCCB(I,J,K)-2.*CCC(I,J,K))
      CCCB(I,J,K)=CCC(I,J,K)
      CCC(I,J,K)=CCCF(I,J,K)
 365  CONTINUE
C 
      CALL DENS
C
      END IF
C
C-----------------------------------------------------------------------
C         COMPUTE UF AND VF
C-----------------------------------------------------------------------
C
      CALL ADVU(DRHOX,ADVUU,DTI2)
      CALL ADVV(DRHOY,ADVVV,DTI2)
      CALL PROFU(DTI2)
      CALL PROFV(DTI2)
      CALL BCOND2(3,DTI2)
C
      DO 369 J=1,JM 
      DO 369 I=1,IM 
 369  TPS(I,J)=0.0
      DO 370 K=1,KBM1
      DO 370 J=1,JM 
      DO 370 I=1,IM 
 370  TPS(I,J)=TPS(I,J)+(UF(I,J,K)+UB(I,J,K)-2.*U(I,J,K))*DZ(K)
      DO 372 K=1,KBM1
      DO 372 J=1,JM 
      DO 372 I=1,IM 
 372  U(I,J,K)=U(I,J,K)+.5*NU*(UF(I,J,K)+UB(I,J,K)-2.*U(I,J,K)
     .        -TPS(I,J))
      DO 373 J=1,JM 
      DO 373 I=1,IM 
 373  TPS(I,J)=0.0
      DO 374 K=1,KBM1
      DO 374 J=1,JM 
      DO 374 I=1,IM 
 374  TPS(I,J)=TPS(I,J)+(VF(I,J,K)+VB(I,J,K)-2.*V(I,J,K))*DZ(K)
      DO 376 K=1,KBM1
      DO 376 J=1,JM 
      DO 376 I=1,IM 
 376  V(I,J,K)=V(I,J,K)+.5*NU*(VF(I,J,K)+VB(I,J,K)-2.*V(I,J,K)
     .        -TPS(I,J))
      DO 377 K=1,KB 
      DO 377 J=1,JM 
      DO 377 I=1,IM 
      UB(I,J,K)=U(I,J,K)
      U(I,J,K)=UF(I,J,K)
      VB(I,J,K)=V(I,J,K)  
 377  V(I,J,K)=VF(I,J,K)
C
      CALL WREAL(DTI2)
C
 8200 CONTINUE

c      write(6,*) 'vr30=',int,v(35,18,1)

C***************drifter call here*************

      startpart = 0    ! default
      jwrite = 0
      irel = 0
      ntrackdays = 60  ! max tracking
      do mp=1,nccctime ! loop drifter start times
c         write(6,*) 'icc ',iccctime(mp),int,nsteps
         if ((int.ge.iccctime(mp)).and.
     &        (int.le.iccctime(mp)+ntrackdays*24*60)) then 
         irel=mp
         if (mod(int-iccctime(mp),iprint).eq.0)iwrite=1
         if (int.eq.iccctime(mp)) startpart=1
         CALL DRIFTER3(mp,root,trackfile)  
c         CALL DRIFTER4(jwrite)  
      endif
      enddo


c************************************************************
c     OUTPUT TO FILES
c************************************************************

c      iprint=1  ! print all time steps
      if ((iwrite.eq.1).and.(mod(int,iprint).eq.0).and.
     &     (iwritedo.eq.1)) then  
c         if (fint.gt.12000)write(6,*)'output at ',
c     &        fint,time,cflmin,cflmin0,horcon
            write(iuelev1) time,el
c            write(iusal1) time,s
c            write(iutemp1) time,t
c            write(iuvel1) time,u
c            write(ivvel1) time,v
c            write(iwvel1) time,w
c            write(iuhvis1) time,aam
c            write(iuvvis1) time,km

      endif
 77   format(9f10.4)

c*****************************************************************

      DO 380 J=1,JM 
      DO 380 I=1,IM 
      EGB(I,J)=EGF(I,J)
      ETB(I,J)=ET(I,J)
      ET(I,J)=ETF(I,J)
      DT(I,J)=H(I,J)+ET(I,J)
      UTB(I,J)=UTF(I,J)
 380  VTB(I,J)=VTF(I,J)
C
      DO 451 I=1,IM
      DO 451 J=1,JM
      DO 451 K=1,KB
      UMONTH(I,J,K)=UMONTH(I,J,K)+U(I,J,K)
      VMONTH(I,J,K)=VMONTH(I,J,K)+V(I,J,K)
      SMONTH(I,J,K)=SMONTH(I,J,K)+S(I,J,K)
  451 CONTINUE
C
 8300 CONTINUE

C
c      IF(MOD(INT,IPRINT).EQ.0)  THEN
c      CALL PRINTS(DRHOX,DRHOY,TRNU,TRNV)
c         CALL PRINTSAL(IUSAL)   ! SALINITY OUTPUT
c         write(6,*) int,time,' outputting limited salinity'
c      ENDIF
C
C     CHECK IF TIME TO CALCULATE AVERAGED SALINITY FIELDS
c      CALL AVGSAL(SAVG,SAVGTM,ISFLAG,IUSAL2)         
c      CALL MONSAL(IUSAL3,MONAVGSAL,MONFLAG,MCNT)
c      CALL SALPOINTS(IUSAL4,OPEN4)  

      
c       if (mod(int,60).eq.0) CALL SAMPOUT

c      IF (MOD(INT,60*6).EQ.0) THEN   
c        WRITE(6,*) 'OUTPUTING ELF ',time
c        WRITE(IUELV1,8998) TIME
c        WRITE(IUELV1,8999) ((EL(I,J), I=1,IM),J=1,JM)
c 8998   format(f10.3)
c 8999   format(10f7.3)
c      ENDIF

c     restart before particle release
c      IF (INT.EQ.ICCCTIME(1)-3) then 
c         WRITE(UNIT=case, FMT='(I1)') mp
c         if ((ibath.eq.1).and.(idist.eq.1))
c     &        OPEN (IUWRS,FORM='UNFORMATTED',FILE='restart_tracks_1')
c         WRITE(IUWRS) 
c     .        intcold,INT,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,AAM,AAM2D,
c     .        ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,WUBOT,WVBOT,
c     .        UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
c     .        VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,KM,KH,KQ,Q2,Q2B,
c     .        Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
c     .        CCC,CCCB,CCCMEAN
c         CLOSE(IUWRS)
c      ENDIF

c     general restart
c      if (mod(int,intrsrt).eq.0) then 
      if (abs(mod(int,intrsrt)-0).lt.dint) then 
c         write(6,*) 'WRITING RESTART AT INT=',int
         OPEN (IUWRS,FORM='UNFORMATTED',FILE='restart')
        WRITE(IUWRS) 
     .     intcold,INT,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,AAM,AAM2D,
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,WUBOT,WVBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .     CCC,CCCB,CCCMEAN
        CLOSE(IUWRS)
      ENDIF

c      CALL ARCHIVE(2)
c      CALL ARCHIVE(3)
C
c      if (mod(int,60).eq.0) 
c      if ((int.gt.100).and.(mod(int,60).eq.0)) call testcfl
c      call testcfl(cflmin)

      fint=fint+dint
c      write(6,*)'-------------------'
 6000 continue   ! iadj, adaptive time loop
      int = fint
      goto 9000
c 9000              CONTINUE
Crealtime
 9001 continue
C
      DO 500 I=1,IM
      DO 500 J=1,JM
      DO 500 K=1,KB
      UMONTH(I,J,K)=UMONTH(I,J,K)/FLOAT(NSTEPS)
      VMONTH(I,J,K)=VMONTH(I,J,K)/FLOAT(NSTEPS)
      SMONTH(I,J,K)=SMONTH(I,J,K)/FLOAT(NSTEPS)
  500 CONTINUE
      CALL PRTXY(UMONTH ,DUM,IM,JM,KBM1,IVAR,10000.,IUPRT)
      CALL PRTXY(VMONTH ,DVM,IM,JM,KBM1,IVAR,10000.,IUPRT)
      CALL PRTXY(SMONTH ,FSM,IM,JM,KBM1,IVAR,100.,IUPRT)
C
      CALL SLICEXZ(VMONTH,DUM,IM,JM,KB,6,PRT,100.,IUPRT)
      CALL SLICEXZ(VMONTH,DUM,IM,JM,KB,19,PRT,100.,IUPRT)
      CALL SLICEXZ(VMONTH,DUM,IM,JM,KB,26,PRT,100.,IUPRT)
      CALL SLICEXZ(VMONTH,DUM,IM,JM,KB,35,PRT,100.,IUPRT)
      CALL SLICEYZ(VMONTH,DVM,IM,JM,KB,11,PRT,100.,IUPRT)
      CALL SLICEYZ(VMONTH,DVM,IM,JM,KB,23,PRT,100.,IUPRT)
      CALL SLICEYZ(VMONTH,DVM,IM,JM,KB,29,PRT,100.,IUPRT)
      CALL SLICEXZ(SMONTH,FSM,IM,JM,KB,6,PRT,1.,IUPRT)
      CALL SLICEXZ(SMONTH,FSM,IM,JM,KB,19,PRT,1.,IUPRT)
      CALL SLICEXZ(SMONTH,FSM,IM,JM,KB,26,PRT,1.,IUPRT)
      CALL SLICEXZ(SMONTH,FSM,IM,JM,KB,35,PRT,1.,IUPRT)
      CALL SLICEYZ(SMONTH,FSM,IM,JM,KB,11,PRT,1.,IUPRT)
      CALL SLICEYZ(SMONTH,FSM,IM,JM,KB,23,PRT,1.,IUPRT)
      CALL SLICEYZ(SMONTH,FSM,IM,JM,KB,19,PRT,1.,IUPRT)
C
      WRITE(IUPRT,602) TIME
602   FORMAT(/2X,'JOB SUCCESSFULLY COMPLETED; TIME = ',1P1E10.2,' DAYS',
     . //)
C
c      WRITE(IUWRS) 
c     .     intcold,IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,AAM,AAM2D,
c     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,WUBOT,WVBOT,
c     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
c     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,KM,KH,KQ,Q2,Q2B,
c     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
c     .     CCC,CCCB,CCCMEAN
c      CLOSE (IUWRS)
C 
Crealtime
c      rlasthour=time*24.
c      write(iulasthr,603) rlasthour

c      print *,'time',time
c      print *,'thour',thour
c      print *,'rlasthour',rlasthour
c      print *,'time*24.',time*24.

603   format(f10.3)

 999  continue 
      close (iulasthr)
c      close (iuseedhr)
      close (iujtime)
c      close (iucccres)

      CLOSE (IUT90)
      CLOSE (IUT91)
      CLOSE (IUT92)
      CLOSE (IUT93)
      CLOSE (IUT94)
c      CALL SYSTEM ('rm gcm_temp*')

c------------------------------------------------------
c  fitness test based on drifter output
c------------------------------------------------------

      call ftest1(trackfile,fscr)

c     get rank file name created by 'modelgroup' (need to change
c     nomenclature from rank to fittest) 
c      open(10,file='rankfilename',status='old')
c      read(10,*) irlen
c      read(10,*) rankfile
c      close(10)      

c     write info to rank file
c      open(10,file=inpath//rankfile(1:irlen),status='old',
c     & access = 'append')
c      write(10,1500) sid,igen,fscr
c      close(10)
c 1500 format(a9,2x,i3,2x,f10.5)    ! 'a' format must = 'nc'

c--------------------------------------------------------
c  output this program has ended
c--------------------------------------------------------
c     get end file name created by 'modelgroup' 
      open(10,file='endfilename',status='old')
      read(10,*) irlen
      read(10,*) rankfile
      close(10)
      open(10,file=rankfile(1:irlen),status='old',access='append')
      write(10,*) sid
      close(10)

      write(6,*) 'NORMAL END OF PROGRAM'
c      STOP
      END
      


C23456789012345678901234567890123456789012345678901234567890123456789012   
c for present-day if n=100 then nump=2273000 if 1879 then nump=2360000
      SUBROUTINE DRIFTER3(MP,root,trackfile)
C     VERSION(06/21/99)
c     VERSION(07/5/2002) BY S.D. MEYERS
c     modified open boundary for MB model 9/15/2004 by Steven Meyers
C     examine model days 120, 220, and 540
C     ONLY INITIALIZES SPECIFIED SUBSET OF GRID
c     updated random number generator for PGI, S. Meyers 8/29/2013

      INCLUDE 'comdeck'  
      include 'genetic'

C
C-----------------------------------------------------------------------
C     THIS SUBROUTINE INTEGRATES PARTICALS FORWARD TO THE 
C            PRESENT INTERNAL TIME STEP
C
C  P here is dimentioned as particle number, parameter
C  where the first par. is I from the grid corner
C            second par. is J from the grid corner
C            third par. is K from the level bottom
C            fourth par. is DX of the particle
C            fifth par. is DY of the particle
C            sixth par. is DZ of the particle
C            seventh par. is age of the particle
C            eighth par. is density of the particle
C
C DX,DY,DZ are all 0 at the start of a grid (thier origin)
C and 1 at the boundary to the next grid, or level
C-----------------------------------------------------------------------
C
      DIMENSION BUP(NUMPX),BUP1(NUMPX),BUP2(NUMPX)
     .,TUP(NUMPX),TUP1(NUMPX),TUP2(NUMPX),UP(NUMPX)
     .,BVP(NUMPX),BVP1(NUMPX),BVP2(NUMPX)
     .,TVP(NUMPX),TVP1(NUMPX),TVP2(NUMPX),VP(NUMPX)
     .,BWP(NUMPX),BWP1(NUMPX),BWP2(NUMPX)
     .,TWP(NUMPX),TWP1(NUMPX),TWP2(NUMPX),WP(NUMPX)
     .,PH2(NUMPX),PH1(NUMPX),PH3(NUMPX)
     .,ut(im,jm,kb),vt(im,jm,kb),wt(im,jm,kb),IPCNT(IM,JM,KBM1),
     .IPCNT2(IM,JM,KBM1)
 
      REAL MOVX,MOVY,MOVZ,clev,dispx,dispy,dispz,dxp,dyp,dzp,px,py,pz
      INTEGER JP,ZKP,PK,IPF,JPF,ZKPF,IPR,INTFIL
      INTEGER II,JJ,MIJ,IX,IY,IZ
      INTEGER pong,stick,ob,pstart,ntemp1,ntemp2,ntemp3,slingshot
      INTEGER*4 LK,SEQ,npskip,ntrack,ip,np
      INTEGER mp,zkpf0  ! index for particle start time
      INTEGER isource  ! flag for initializing from file (1) or self-generate (2)
      integer ilen,ilen2
      character*80 case,cfile,tfile,trackfile
      character*2 si,sj
      character*1 s1,sreg
c      character*17  tfile,tfile0
      character*35 outpath
      character*80 root,root2

c      write(6,*) 'RUNNING DRFTRS INT=',INT,root,igen

c      npskip = 1 ! number of particles to skip when writing traj position file
      npskip = 1  ! for 2reg
      SEQ=INT
      outpath = '/home/meyers/Sci/ECOM/GA/PT/Output/'
      seedx = -1
      call random_seed(seedx)   ! initialize PGI random number generator

      if (mp.lt.10) WRITE(UNIT=case, FMT='(I1)') mp
      if (mp.ge.10) WRITE(UNIT=case, FMT='(I2)') mp

c      tfile0 = 'dummy'  ! dummy initial value

      IUTR=105
      IUTRS=104
      IURES=106
      INTFIL = 100

c   create output file names
       ilen = len(root)
c       write(6,*) 'ROOT=',ilen,trim(root)//'_'//trim(sid)
       root2 = trim(root)//'_'//trim(sid) 
       ilen2 = len(root2)
     
      if (igen.lt.10) then
        write(sgen,'(i1)') igen
        sgen = "0"//sgen
      endif
      if ((igen.ge.10).and.(igen.lt.100)) then
        write(sgen,'(i2)') igen
      endif
      root2 = trim(root2)//'_'//trim(sgen)
      ilen2 = len(root2)
c      write(6,*) 'root2=',root2(1:ilen2),' ',sgen(1:2)


      IF (startpart .EQ. 1) THEN
c        write(6,*) 'INITIALIZING DRIFTERS'       
        write(6,*)'partcl output:',outpath//'pcount_'//root2(1:ilen2)
        trackfile = outpath//'tracks_'//root2(1:ilen2) 
        OPEN(IUTR,FILE=outpath//'pcount_'//root2(1:ilen2))
        OPEN(IURES,FILE=outpath//'tracks_'//root2(1:ilen2))
c        write(6,*) 'starting particle tracking at ',int
      ELSE
        OPEN(IUTR,FILE=outpath//'pcount_'//root2(1:ilen2),
     &        access='append')
        OPEN(IURES,FILE=outpath//'tracks_'//root2(1:ilen2),
     &       access='append')
c        write(6,*) 'reopening  track file'
      ENDIF                     ! startpart=1

c the number of time steps between outputs IPR
      IPR=IPRINT
      PONG=1
      stick=0
      ob=3  ! otb baseline
c      slingshot = 1
      clev=1.0
      pstart=1

c      dispx=0.005  ! Burwell original as fixed value
c      dispy=0.005
      dispx=0.0087   ! Meyers April 2014
      dispy=dispx
      dispz=dispx*0.0002  ! vertical/horizontal aspect ratio

      do 3 i=1,im
      do 3 j=1,jm
      do 3 k=1,kb
         ut(i,j,k)=u(i,j,k)*dum(i,j)
         vt(i,j,k)=v(i,j,k)*dvm(i,j)
         wt(i,j,k)=w(i,j,k)
 3    continue
 
      DO 1 LK=1,NUMPX
         BUP(LK)=0.0
         BUP1(LK)=0.0
         BUP2(LK)=0.0
         TUP(LK)=0.0
         TUP1(LK)=0.0
         TUP2(LK)=0.0
         UP(LK)=0.0
         BVP(LK)=0.0
         BVP1(LK)=0.0
         BVP2(LK)=0.0
         TVP(LK)=0.0
         TVP1(LK)=0.0
         TVP2(LK)=0.0
         VP(LK)=0.0
         BWP(LK)=0.0
         BWP1(LK)=0.0
         BWP2(LK)=0.0
         TWP(LK)=0.0
         TWP1(LK)=0.0
         TWP2(LK)=0.0
         WP(LK)=0.0
         PH1(LK)=0.0
         PH2(LK)=0.0
         PH3(LK)=0.0

 1    CONTINUE
      MOVX=0.0
      MOVY=0.0
      MOVZ=0.0
      
c     always start count from zero
      do 100 i=1,im
      do 100 j=1,jm
      do 100 k=1,kbm1
         ipcnt(i,j,k)=0
         ipcnt2(i,j,k)=0
 100  continue

c      seedx=-1
c      seedy=-2
c      seedz=-3

c--------------------------------------------------------------------------
c  initialize particles based on isource: 1-read 2-self generate
c--------------------------------------------------------------------------
c     initialize drifter positions
      IF (startpart.eq.1) THEN 
c         write(6,*) 'INITIALIZING DRFTRS INT=',INT,ibuoy
         
c.............................................................................
c     read from file
         lk=0
         if (idist.eq.2) then 
            write(6,*) 'reading initial positions',ibath,ireg,ibuoy
            
c    define INPUT filename holding initial positions at 1m depth
c            ilen = 12
c            root(1:ilen) = 'init_tracks_'
c            ilen = 8
cc            root(1:ilen) = 'tracks0_'
c            if (ibath.eq.1) then 
c               root(1:ilen+5) = root(1:ilen)//'Pres_'
c               ilen = ilen+5
c            endif
c            if (ibath.eq.2) then 
c               root(1:ilen+5) = root(1:ilen)//'1879_'
c               ilen = ilen+5
c            endif
c            if (ireg.eq.2) then 
c               root(1:ilen+6) = root(1:ilen)//'100reg'
c               ilen = ilen+6
c            endif            
c            if (ireg.ne.2) then 
c               root(1:ilen+3) = root(1:ilen)//'reg'
c               ilen = ilen+3
c            endif
            
            write(6,*) 'reading initial positions ',root(1:ilen)
            OPEN(INTFIL,FILE=root(1:ilen),status='old')
            read(intfil,*) np
            if (np.gt.numpx) then 
              write(6,*) 'ERROR: number of particles in file .gt. NUMPX'
            endif
c            write(6,*) 'reading and creating particle positions',np
            DO 500  ip=1,np     ! read predefined (x,y) positions
               READ(INTFIL,*) px,py,pz !,MIJ  ! this pz ~ 1 m depth, 1reg* files
               anz=1.
               ip0=0
               ip1=0

               if (ireg.eq.2) then ! middle layers
                  ip0 = 2
                  ip1 = 7
               endif
               if (ireg.eq.3) then ! mid-depth
                  ip0 = 6
                  ip1 = 6
                  anz = 1
               endif

               if (ibuoy.eq.1) then  ! floating
                  ip0 = 1
                  ip1 = 1
                  anz = 1
               endif

c        Lines 500-550 distribute particles through layers  ! 2reg* files
               do 550 pz=ip0,ip1    !fine which layers to populate
                  if ((ireg.eq.2).and.((pz.eq.5).or.(pz.eq.6)))anz = 2
                  if ((ireg.eq.2).and.((pz.eq.4).or.(pz.eq.7)))anz = 1
                  if ((ireg.eq.2).and.((pz.le.3).or.(pz.ge.8)))anz = 0
                  if (ireg.eq.1) write(6,*) 'IREG ERROR @550 DRIFTER3'

               do 555 iz=1,anz  ! relative depth within layer

               lk = lk+1
c               write(6,*) lk,ip,np,px,py,pz
c               isactive(lk) = 1
               P(LK,1) = IFIX(px) ! cell number
               P(LK,2) = IFIX(py) ! cell number
               P(LK,3) = IFIX(pz) ! cell number
               P(LK,4) = px-P(LK,1)
               P(LK,5) = py-P(LK,2)
               P(LK,6) = pz-P(LK,3)
c     options for initial vertical disribution
               if (ireg.eq.3)p(lk,6) = 0.
               if (ireg.eq.2)P(LK,6) = (FLOAT(IZ-1))/(anz)+1./anz/2.

               P(LK,7) = 1.0
               P(LK,8) = 0.0
               P(LK,9) = 1.0
               P(LK,10) = INT
               n1=ifix(px)
               n2=ifix(py)
               n3=ifix(pz)
               ipcnt(n1,n2,n3)=ipcnt(n1,n2,n3)+1
c     write particle positions to track file
 555           continue  ! iz
 550           continue  ! k
 500        CONTINUE     ! ip
            close(INTFIL)
c            write(6,*) 'finished initial positions'
c     drop few particles to round to 10 
c            icut = lk-10*IFIX(lk/10.)
c            lk = lk-icut
c            write(6,*) 'cutting ',icut,lk

            nump=lk
            np = nump
            write(6,*)'NP=',np,numpx
            startpart=0    ! reset
         endif                  ! idist=2
            
c.............................................................................
c     generate  particle distribution as same # per cell
         if (idist.eq.1) then 
            lk=0
c            write(6,*) 'creating initial positions'
            anx = 2.            ! number within each cell
            any = 2.
            anz = 1.
            DO 400 II=1,IM-1  ! loop all cells in model
            DO 400 JJ=1,JM-1
c           proceed for sea-cells
               if (fsm(ii,jj).eq.0) goto 410
               do 11 ie=1,numebc  ! do not initialize drifters in open boundary
               if ((ii.eq.ieta(ie)).and.(jj.eq.jeta(ie)))goto 410
 11         enddo

            DO 420 k=1,8       !loop vertical cells
            DO 50 IX=1,anx     ! DISTRIBUTE WITHIN EACH CELL
            DO 50 IY=1,any
            DO 50 IZ=1,anz
               lk = lk+1
               P(LK,1) = FLOAT(II) ! cell number
               P(LK,2) = FLOAT(JJ)
               P(LK,3) = FLOAT(K)
               P(LK,4) = (FLOAT(IX)-1.)/(anx)!+1./anx/2. ! fraction of cell. 
               P(LK,5) = (FLOAT(IY)-1.)/(any)!+1./any/2.   
               P(LK,6) = 0.5  !(FLOAT(IZ)-1.)/(anz)
               P(LK,7) = 1.0
               P(LK,8) = 0.0
               P(LK,9) = 1.0
               P(LK,10) = INT
               ipcnt(ii,jj,k)=ipcnt(II,JJ,K)+1
 50         continue
 420        continue
 410        continue
 400        continue
            nump=lk
c............................................................................. 
c            WRITE(6,*) 'TOTAL # OF DRIFTERS=',LK,NUMP
            if (nump.gt.numpx) then 
              write(6,*) 'ERROR: number of particles in file .gt. NUMPX'
            endif
          endif                  ! idist=1
c.............................................................................
c  how many particles to print to track file
c  print initial configuration
c.............................................................................
         ntrack = nump/npskip
c         write(6,*) '==========>   NTRACK=',ntrack
         write(iures,*) ntrack
c         write(iures,*) 'Time'  ! used for easier searching of output file
         write(iures,*) int
         DO 401  LK=1,nump,NPSKIP ! horizontal number of cells
            WRITE(IURES,177) p(lk,1)+p(lk,4),p(lk,2)+p(lk,5),
     &           p(lk,3)+p(lk,6)
 401     CONTINUE

c         WRITE(IUTR,*)  'TIME STEP'
         WRITE(IUTR,*)  SEQ

         DO 12 k=1,kbm1
         do 12 j=1,jm
         do 12 i=1,im
            WRITE(IUTR,*)  ipcnt(I,J,K)
 12      continue

      ENDIF                     ! startpart=1
C          
c***********************************************************************
c  interpolate velocity field then integrate particle positions
C     For every particle               
c***********************************************************************
c      write(6,*) 'NUMP=',nump
      DO 10 LK=1,NUMP         
C     set old i,j,k
         IP=IFIX(P(LK,1))
         JP=IFIX(P(LK,2))
         ZKP=IFIX(P(LK,3))

         if (fsm(ip,jp).eq.0) goto 10  ! skip if particle on land point
         
C     linerly interp h1

         IF (P(LK,5) .GT. 0.5) THEN
            ph1(LK)=H1(IP,JP)+(P(LK,5)-0.5)
     &           *(H1(IP,JP+1)-H1(IP,JP))
         ELSE
            ph1(LK)=H1(IP,JP-1)+(P(LK,5)+0.5)
     &           *(H1(IP,JP)-H1(IP,JP-1))
         ENDIF
         
C     linerly interp h2
         
         IF (P(LK,4) .GT. 0.5) THEN
            ph2(LK)=H2(IP,JP)+(P(LK,4)-0.5)
     &           *(H2(IP+1,JP)-H2(IP,JP))
         ELSE
            ph2(LK)=H2(IP-1,JP)+(P(LK,4)+0.5)
     &           *(H2(IP,JP)-H2(IP-1,JP))
         ENDIF
         
C     linerly interp sigma
         
         ph3(LK)=DZ(ZKP)
     &        +2.0*(P(LK,6)*DZ(ZKP)+DZ(ZKP-1))
     &        *(DZ(ZKP)-DZ(ZKP-1))/(DZ(ZKP)+DZ(ZKP-1))
         
         
C     linerly interp U
         
C     if dy is further than half way through the grid
C     then the closest neighbors are IP+1,JP IP+1,JP+1
C     IP,JP IP,JP+1
         
         IF (P(LK,5) .GT. 0.5) THEN
            
C     if dz is lt .5  the lower level is ZKP+1 and the upper is ZKP
C     and the particle is in the ZKP level
            
            IF (P(LK,6) .LT. 0.5) THEN
C     calculate the lower two linear interp U's

               bup1(LK)=UT(IP+1,JP,ZKP+1)+(P(LK,5)-0.5)
     &              *(UT(IP+1,JP+1,ZKP+1)-UT(IP+1,JP,ZKP+1))
               bup2(LK)=UT(IP,JP,ZKP+1)+(P(LK,5)-0.5)
     &              *(UT(IP,JP+1,ZKP+1)-UT(IP,JP,ZKP+1))
               
C     combine by linearly interp those two u's
               
               bup(LK)=bup2(LK) + P(LK,4)
     &              *(bup1(LK)-bup2(LK))
               
C     do the same for the level above
               
               tup1(LK)=UT(IP+1,JP,ZKP)+(P(LK,5)-0.5)
     &              *(UT(IP+1,JP+1,ZKP)-UT(IP+1,JP,ZKP))
               tup2(LK)=UT(IP,JP,ZKP)+(P(LK,5)-0.5)
     &              *(UT(IP,JP+1,ZKP)-UT(IP,JP,ZKP))
               tup(LK)=tup2(LK) + P(LK,4)*(tup1(LK)-tup2(LK))
               
C     combine the upper and lower by linear interp
               
               up(LK)=bup(LK)+(2.0*P(LK,6)*DZ(ZKP)+DZ(ZKP+1))
     &              *(tup(LK)-bup(LK))/(DZ(ZKP)+DZ(ZKP+1))
               
C     if dz is ge .5 the lower grid is ZKP and the upper grid is ZKP-1
C     the particle is still in level ZKP
               
            ELSE
               
c      write(6,*)'calculate the lower two linear interp Us '
C     calculate the lower two linear interp U's 
               
               bup1(LK)=UT(IP+1,JP,ZKP)+(P(LK,5)-0.5)
     &              *(UT(IP+1,JP+1,ZKP)-UT(IP+1,JP,ZKP))
               bup2(LK)=UT(IP,JP,ZKP)+(P(LK,5)-0.5)
     &              *(UT(IP,JP+1,ZKP)-UT(IP,JP,ZKP))
               
C     combine by linearly interp those two u's
               
               bup(LK)=bup2(LK) + P(LK,4)
     &              *(bup1(LK)-bup2(LK))
               
C     do the same for the level above
               
               tup1(LK)=UT(IP+1,JP,ZKP-1)+(P(LK,5)-0.5)
     &              *(UT(IP+1,JP+1,ZKP-1)-UT(IP+1,JP,ZKP-1))
               tup2(LK)=UT(IP,JP,ZKP-1)+(P(LK,5)-0.5)
     &              *(UT(IP,JP+1,ZKP-1)-UT(IP,JP,ZKP-1))

               tup(LK)=tup2(LK) + P(LK,4)*(tup1(LK)-tup2(LK))
               
C     combine the upper and lower via linear interp.
               
               up(LK)=bup(LK)+(2.0*(P(LK,6)-0.5)*DZ(ZKP)+DZ(ZKP-1))
     &              *(tup(LK)-bup(LK))/(DZ(ZKP)+DZ(ZKP-1))
               
            ENDIF
            
C     if dy is less than or equal to half the grid dist
C     then the closet neighbors are JP-1
            
         ELSE
            
C     if dz is lt .5  the lower grid is ZKP+1 and the upper is ZKP
C     and the particle is in the ZKP level
            
            IF (P(LK,6) .LT. 0.5) THEN
               
C     calculate the lower two linear interp U's
               
               bup1(LK)=UT(IP+1,JP-1,ZKP+1)+(P(LK,5)+0.5)
     &              *(UT(IP+1,JP,ZKP+1)-UT(IP+1,JP-1,ZKP+1))
               bup2(LK)=UT(IP,JP-1,ZKP+1)+(P(LK,5)+0.5)
     &              *(UT(IP,JP,ZKP+1)-UT(IP,JP-1,ZKP+1))
               
C     combine by linearly interp those two u's
               
               bup(LK)=bup2(LK) + P(LK,4)
     &              *(bup1(LK)-bup2(LK))
               
C     do the same for the level above
               
               tup1(LK)=UT(IP+1,JP-1,ZKP)+(P(LK,5)+0.5)
     &              *(UT(IP+1,JP,ZKP)-UT(IP+1,JP-1,ZKP))
               tup2(LK)=UT(IP,JP-1,ZKP)+(P(LK,5)+0.5)
     &              *(UT(IP,JP,ZKP)-UT(IP,JP-1,ZKP))

               tup(LK)=tup2(LK) + P(LK,4)*(tup1(LK)-tup2(LK))
               
C     combine the upper and lower by linear interp
               
               up(LK)=bup(LK)+(2.0*P(LK,6)*DZ(ZKP)+DZ(ZKP+1))
     &              *(tup(LK)-bup(LK))/(DZ(ZKP)+DZ(ZKP+1))

c         if (mod(lk,100).eq.0) then 
c            write(6,*) 'UP ',LK,ip,jp,zkp,up(lk)
c            write(6,*) bup(lk),tup(lk)
c         endif
C     if dz is ge .5 the lower grid is ZKP and the upper grid is ZKP-1
C     the particle is still in level ZKP
               
            ELSE
               
c      write(6,*)'calculate the lower two linear interp Us'
C     calculate the lower two linear interp U's
               
               bup1(LK)=UT(IP+1,JP-1,ZKP)+(P(LK,5)+0.5)
     &              *(UT(IP+1,JP,ZKP)-UT(IP+1,JP-1,ZKP))
               bup2(LK)=UT(IP,JP-1,ZKP)+(P(LK,5)+0.5)
     &              *(UT(IP,JP,ZKP)-UT(IP,JP-1,ZKP))
               
C     combine by linearly interp those two u's
               
               bup(LK)=bup2(LK) + P(LK,4)
     &              *(bup1(LK)-bup2(LK))
               
C     do the same for the level above
               
               tup1(LK)=UT(IP+1,JP-1,ZKP-1)+(P(LK,5)+0.5)
     &              *(UT(IP+1,JP,ZKP-1)-UT(IP+1,JP-1,ZKP-1))
               tup2(LK)=UT(IP,JP-1,ZKP-1)+(P(LK,5)+0.5)
     &              *(UT(IP,JP,ZKP-1)-UT(IP,JP-1,ZKP-1))

               tup(LK)=tup2(LK) + P(LK,4)*(tup1(LK)-tup2(LK))
               
C     combine the upper and lower via linear interp.
               
               up(LK)=bup(LK)+(2.0*(P(LK,6)-0.5)*DZ(ZKP)+DZ(ZKP-1))
     &              *(tup(LK)-bup(LK))/(DZ(ZKP)+DZ(ZKP-1))
               
            ENDIF
            
         ENDIF
         
c      write(6,*)'linearly interp V'
C     linerly interp V
         
C     if dx is further than half way through the grid
C     then the closest neighbors are IP+1,JP IP+1,JP+1
C     IP,JP IP,JP+1
         
         IF (P(LK,4) .GT. 0.5) THEN
            
C     if dz is lt .5  the lower grid is ZKP+1 and the upper is ZKP
C     and the particle is in the ZKP level
            
            IF (P(LK,6) .LT. 0.5) THEN
               
C     calculate the lower two linear interp U's
               
               bvp1(LK)=VT(IP+1,JP,ZKP+1)+(P(LK,4)-0.5)
     &              *(VT(IP+1,JP+1,ZKP+1)-VT(IP+1,JP,ZKP+1))
               bvp2(LK)=VT(IP,JP,ZKP+1)+(P(LK,4)-0.5)
     &              *(VT(IP,JP+1,ZKP+1)-VT(IP,JP,ZKP+1))
               
C     combine by linearly interp those two u's
               
               bvp(LK)=bvp2(LK) + P(LK,5)
     &              *(bvp1(LK)-bvp2(LK))
               
C     do the same for the level above
               
               tvp1(LK)=VT(IP+1,JP,ZKP)+(P(LK,4)-0.5)
     &              *(VT(IP+1,JP+1,ZKP)-VT(IP+1,JP,ZKP))
               tvp2(LK)=VT(IP,JP,ZKP)+(P(LK,4)-0.5)
     &              *(VT(IP,JP+1,ZKP)-VT(IP,JP,ZKP))
               tvp(LK)=tvp2(LK) + P(LK,5)*(tvp1(LK)-tvp2(LK))
               
C     combine the upper and lower by linear interp
               
               vp(LK)=bvp(LK)+(2.0*P(LK,6)*DZ(ZKP)+DZ(ZKP+1))
     &              *(tvp(LK)-bvp(LK))/(DZ(ZKP)+DZ(ZKP+1))
               
C     if dz is ge .5 the lower grid is ZKP and the upper grid is ZKP-1
C     the particle is still in level ZKP
               
            ELSE
               
C     calculate the lower two linear interp U's 
               
               bvp1(LK)=VT(IP+1,JP,ZKP)+(P(LK,4)-0.5)
     &              *(VT(IP+1,JP+1,ZKP)-VT(IP+1,JP,ZKP))
               bvp2(LK)=VT(IP,JP,ZKP)+(P(LK,4)-0.5)
     &              *(VT(IP,JP+1,ZKP)-VT(IP,JP,ZKP))
               
C     combine by linearly interp those two u's
               
               bvp(LK)=bvp2(LK) + P(LK,5)
     &              *(bvp1(LK)-bvp2(LK))
               
C     do the same for the level above
               
               tvp1(LK)=VT(IP+1,JP,ZKP-1)+(P(LK,4)-0.5)
     &              *(VT(IP+1,JP+1,ZKP-1)-VT(IP+1,JP,ZKP-1))
               tvp2(LK)=VT(IP,JP,ZKP-1)+(P(LK,4)-0.5)
     &              *(VT(IP,JP+1,ZKP-1)-VT(IP,JP,ZKP-1))
               tvp(LK)=tvp2(LK) + P(LK,5)*(tvp1(LK)-tvp2(LK))
               
C     combine the upper and lower via linear interp.
               
               vp(LK)=bvp(LK)+(2.0*(P(LK,6)-0.5)*DZ(ZKP)+DZ(ZKP-1))
     &              *(tvp(LK)-bvp(LK))/(DZ(ZKP)+DZ(ZKP-1))
               
            ENDIF
            
C     if dy is less than or equal to half the grid dist
C     then the closet neighbors are IP,JP IP-1,JP
C     IP,JP+1 IP-1,JP+1
            
            
            
         ELSE
            
C     if dz is lt .5  the lower grid is ZKP+1 and the upper is ZKP
C     and the particle is in the ZKP level
            
            IF (P(LK,6) .LT. 0.5) THEN
               
C     calculate the lower two linear interp U's
               
               bvp1(LK)=VT(IP,JP,ZKP+1)+(P(LK,5)+0.5)
     &              *(VT(IP-1,JP,ZKP+1)-VT(IP,JP,ZKP+1))
               bvp2(LK)=VT(IP,JP+1,ZKP+1)+(P(LK,5)+0.5)
     &              *(VT(IP-1,JP+1,ZKP+1)-VT(IP,JP+1,ZKP+1))
               
C     combine by linearly interp those two u's
               
               bvp(LK)=bvp2(LK) + P(LK,4)
     &              *(bvp1(LK)-bvp2(LK))
               
C     do the same for the level above
               
               tvp1(LK)=VT(IP,JP,ZKP)+(P(LK,5)+0.5)
     &              *(VT(IP-1,JP,ZKP)-VT(IP,JP,ZKP))
               tvp2(LK)=VT(IP,JP+1,ZKP)+(P(LK,5)+0.5)
     &              *(VT(IP-1,JP+1,ZKP)-VT(IP,JP+1,ZKP))
               tvp(LK)=tvp2(LK) + P(LK,4)*(tvp1(LK)-tvp2(LK))
               
C     combine the upper and lower by linear interp
               
               vp(LK)=bvp(LK)+(2.0*P(LK,6)*DZ(ZKP)+DZ(ZKP+1))
     &              *(tvp(LK)-bvp(LK))/(DZ(ZKP)+DZ(ZKP+1))
               
c        if (mod(lk,100).eq.0) then 
c            write(6,*) 'VP ',LK,ip,jp,zkp,vp(lk)
c            write(6,*) bvp(lk),tvp(lk)
c         endif

C     if dz is ge .5 the lower grid is ZKP and the upper grid is ZKP-1
C     the particle is still in level ZKP
               
            ELSE
               
C     calculate the lower two linear interp U's
               
               bvp1(LK)=VT(IP,JP,ZKP)+(P(LK,5)+0.5)
     &              *(VT(IP-1,JP,ZKP)-VT(IP,JP,ZKP))
               bvp2(LK)=VT(IP,JP+1,ZKP)+(P(LK,5)+0.5)
     &              *(VT(IP-1,JP+1,ZKP)-VT(IP,JP+1,ZKP))
               
C     combine by linearly interp those two u's
               
               bvp(LK)=bvp2(LK) + P(LK,4)
     &              *(bvp1(LK)-bvp2(LK))
               
C     do the same for the level above
               
               tvp1(LK)=VT(IP,JP,ZKP-1)+(P(LK,5)+0.5)
     &              *(VT(IP-1,JP,ZKP-1)-VT(IP,JP,ZKP-1))
               tvp2(LK)=VT(IP,JP+1,ZKP-1)+(P(LK,5)+0.5)
     &              *(VT(IP-1,JP+1,ZKP-1)-VT(IP,JP+1,ZKP-1))
               tvp(LK)=tvp2(LK) + P(LK,4)*(tvp1(LK)-tvp2(LK))
               
C     combine the upper and lower via linear interp.
               
               vp(LK)=bvp(LK)+(2.0*(P(LK,6)-0.5)*DZ(ZKP)+DZ(ZKP-1))
     &              *(tvp(LK)-bvp(LK))/(DZ(ZKP)+DZ(ZKP-1))
               
c        if (mod(lk,100).eq.0) then 
c            write(6,*) 'VP ',LK,ip,jp,zkp,vp(lk)
c            write(6,*) bvp(lk),tvp(lk)
c         endif
            ENDIF
            
         ENDIF
         
         
c      write(6,*)'linearly interp W'
C     linerly interp W
         
         IF (P(LK,4) .GE. 0.5) THEN
            IF (P(LK,5) .GE. 0.5) THEN
               bwp1(LK)=WT(IP+1,JP,ZKP)+(P(LK,5)-0.5)
     &              *(WT(IP+1,JP+1,ZKP)-WT(IP+1,JP,ZKP))
               bwp2(LK)=WT(IP,JP,ZKP)+(P(LK,5)-0.5)
     &              *(WT(IP,JP+1,ZKP)-WT(IP,JP,ZKP))
               bwp(LK)=bwp2(LK) + P(LK,4)
     &              *(bwp1(LK)-bwp2(LK))
               twp1(LK)=WT(IP+1,JP,ZKP-1)+(P(LK,5)-0.5)
     &              *(WT(IP+1,JP+1,ZKP-1)-WT(IP+1,JP,ZKP-1))
               twp2(LK)=WT(IP,JP,ZKP-1)+(P(LK,5)-0.5)
     &              *(WT(IP,JP+1,ZKP-1)-WT(IP,JP,ZKP-1))
               twp(LK)=twp2(LK) + P(LK,4)*(twp1(LK)-twp2(LK))
               wp(LK)=bwp(LK)+(P(LK,6))*(twp(LK)-bwp(LK))
            ELSE
               bwp1(LK)=WT(IP+1,JP-1,ZKP)+(P(LK,5))
     &              *(WT(IP+1,JP,ZKP)-WT(IP+1,JP-1,ZKP))
               bwp2(LK)=WT(IP,JP-1,ZKP)+(P(LK,5))
     &              *(WT(IP,JP,ZKP)-WT(IP,JP-1,ZKP))
               bwp(LK)=bwp2(LK) + P(LK,4)
     &              *(bwp1(LK)-bwp2(LK))
               twp1(LK)=WT(IP+1,JP-1,ZKP-1)+(P(LK,5))
     &              *(WT(IP+1,JP,ZKP-1)-WT(IP+1,JP-1,ZKP-1))
               twp2(LK)=WT(IP,JP-1,ZKP-1)+(P(LK,5))
     &              *(WT(IP,JP,ZKP-1)-WT(IP,JP-1,ZKP-1))
               twp(LK)=twp2(LK) + P(LK,4)*(twp1(LK)-twp2(LK))
               wp(LK)=bwp(LK)+(P(LK,6))*(twp(LK)-bwp(LK))
               
            ENDIF
         ELSE
            
            IF (P(LK,5) .GE. 0.5) THEN
               bwp1(LK)=WT(IP,JP,ZKP)+(P(LK,5)-0.5)
     &              *(WT(IP,JP+1,ZKP)-WT(IP,JP,ZKP))
               bwp2(LK)=WT(IP-1,JP,ZKP)+(P(LK,5)-0.5)
     &              *(WT(IP-1,JP+1,ZKP)-WT(IP-1,JP,ZKP))
               bwp(LK)=bwp2(LK) + P(LK,4)
     &              *(bwp1(LK)-bwp2(LK))
               twp1(LK)=WT(IP,JP,ZKP-1)+(P(LK,5)-0.5)
     &              *(WT(IP,JP+1,ZKP-1)-WT(IP,JP,ZKP-1))
               twp2(LK)=WT(IP-1,JP,ZKP-1)+(P(LK,5)-0.5)
     &              *(WT(IP-1,JP+1,ZKP-1)-WT(IP-1,JP,ZKP-1))
               twp(LK)=twp2(LK) + P(LK,4)*(twp1(LK)-twp2(LK))
               wp(LK)=bwp(LK)+(P(LK,6))*(twp(LK)-bwp(LK))
            ELSE
               bwp1(LK)=WT(IP,JP-1,ZKP)+(P(LK,5))
     &              *(WT(IP,JP,ZKP)-WT(IP,JP-1,ZKP))
               bwp2(LK)=WT(IP-1,JP-1,ZKP)+(P(LK,5))
     &              *(WT(IP-1,JP,ZKP)-WT(IP-1,JP-1,ZKP))
               bwp(LK)=bwp2(LK) + P(LK,4)
     &              *(bwp1(LK)-bwp2(LK))
               twp1(LK)=WT(IP,JP-1,ZKP-1)+(P(LK,5))
     &              *(WT(IP,JP,ZKP-1)-WT(IP,JP-1,ZKP-1))
               twp2(LK)=WT(IP-1,JP-1,ZKP-1)+(P(LK,5))
     &              *(WT(IP-1,JP,ZKP-1)-WT(IP-1,JP-1,ZKP-1))
               twp(LK)=twp2(LK) + P(LK,4)*(twp1(LK)-twp2(LK))

               wp(LK)=bwp(LK)+(P(LK,6))*(twp(LK)-bwp(LK))
               
            ENDIF
            
         ENDIF
         
         
C     Integrate the vel. and add to old position
c         MOVX=UP(LK)*p(lk,7)*p(lk,9)*DTI/ph1(lk)
         MOVX=UP(LK)*p(lk,7)*p(lk,9)*DTI/h1(ip,jp)
         call random_number(sx)
         sx = (sx-0.5)*dispx*p(lk,7)*p(lk,9)
         MOVX=MOVX+sx

         MOVY=VP(LK)*p(lk,7)*p(lk,9)*DTI/h2(ip,jp)
         call random_number(sy)
         sy = (sy-0.5)*dispy*p(lk,7)*p(lk,9)
         MOVY = MOVY+sy
         
         MOVZ=WP(LK)*p(lk,7)*p(lk,9)*DTI/dz(zkp)
         call random_number(sz)
         sz = (sz-0.5)*dispz*p(lk,7)*p(lk,9)
         MOVZ=MOVZ+sz
         if (ibuoy.eq.1) then  ! floating
            MOVZ = 0.
         endif

C     Update the I, DX position
c         write(6,*) seq,lk,p(lk,1),p(lk,4),movx
         p04 = p(lk,4)
         p4 = p04+movx
         if (p4.gt.1) then  ! modified May 2014 SD Meyers
            p4 = (p4-1)*h1(ip,jp)/h1(ip+1,jp)
            p(lk,1) = p(lk,1)+1
         endif
         if (p4.lt.0) then  ! modified May 2014 SD Meyers
            p(lk,1) = p(lk,1)-1
            p4 = 1+p4*h1(ip,jp)/h1(ip-1,jp)
         endif
         p(lk,4)=p4
         
C     Update the J, DY position
         p05 = p(lk,5)
         p5 = p05+movy
         if (p5.ge.1) then  ! modified May 2014 SD Meyers
            p(lk,2) = p(lk,2)+1
            p5 = (p5-1.)*h2(ip,jp)/h2(ip,jp+1)
         endif
         if (p5.lt.0) then  ! modified May 2014 SD Meyers
            p(lk,2) = p(lk,2)-1
            p5 = 1+p5*h2(ip,jp)/h2(ip,jp-1)
         endif
         p(lk,5)=p5
         
C     Update the K, DZ position
         p06 = p(lk,6)
         p6 = p06-movz  !use -movz bc downward w is >0 wrt layers
         if ((p6.ge.1).and.(p(lk,3)+1.lt.10)) then  ! modified May 2014 SD Meyers
            p6 = (p6-1.)*dz(zkp)/dz(zkp+1)
            p(lk,3) = p(lk,3)+1
         endif
         if ((p6.lt.0).and.(p(lk,3)-1.gt.1)) then  ! modified May 2014 SD Meyers
            p6 = 1+p6*dz(zkp)/dz(zkp-1)
             p(lk,3) = p(lk,3)-1
         endif
         if (p6.lt.0) p6 = 0.
         if (p6.gt.1) p6 = 0.99
         p(lk,6)=p6

         if (ibuoy.eq.1) then  ! floating
            p(lk,3) = 1
            p(lk,6) = 0.5
         endif
         if (abs(p6).gt.1)write(6,*)'P6 ERROR:',int,lk,p(lk,3),p06,p6

         IPF=IFIX(P(LK,1))
         JPF=IFIX(P(LK,2))
         ZKPF=IFIX(P(LK,3))

c------------------------------------------------------------         
c  interpolate depth for change in column (i,j)
c  z(k) < 0  relative depths
c  dz(k) >0 relative thickness of layers 
c------------------------------------------------------------         
c         if (((ipf.ne.ip).or.(jpf.ne.jp)).and.
c     &        (fsm(ipf,jpf).eq.1).and.(fsm(ip,jp).eq.1)) then 
c            zkpf0 = zkpf
c            hold = h(ip,jp)    ! old cell depth
c            hnew = h(ipf,jpf)  ! new cell depth
c            d0 = -hold*z(zkpf0)+hold*dz(zkpf0)*p(lk,6)    ! depth of particle in meters
c
c            ! find which layer particle is now within in new column
c            do kk=1,kbm1  
c               if (abs(hnew*z(kk)).le.d0) zkpf=kk
cc               if (lk.eq.87971)write(6,*)d0,kk,abs(hnew*z(kk)),zkpf
c            enddo
c
c            ! traps particles if jump in depth not realistic
c            delz = (d0+hnew*z(zkpf))  ! depth within cell in meters
c            if (d0.lt.hnew) then   ! pong off side walls
c               p(lk,6) = delz/(hnew*dz(zkpf))
c               p(lk,3) = zkpf
c            else
c               ipf = ip  ! pong off side walls
c               jpf = jp
cc               p(lk,3) = 9  ! start at bottom
cc               p(lk,6) = 0.5
c            endif
c
c            if (ibuoy.eq.1) then ! floating
c               p(lk,3) = 1
c               p(lk,6) = 0.5
c            endif
c
cc     validate computation
c            d1 = -hnew*z(zkpf)+hnew*dz(zkpf)*p(lk,6)    ! depth of particle in meters
c       if ((p(lk,6).gt.1))then 
c          write(6,*)'K ERROR ',
c     &    seq,lk,hold,hnew,d0,d1,zkpf0,zkpf,z(zkpf)*hnew,z(zkpf+1)*hnew,
c     &           delz,delz/dz(zkpf)
cc          stop
c       endif
c         endif

c-------------------------------------------------------
c do OBC before pong
c-------------------------------------------------------
c         if (ob.eq.3) THEN  ! OTB baseline
cc         if (lk.eq.1) write(6,*) 'OTB OBC'
cc           if((jpf.le.4).and.(lk.eq.34887))
cc     &    write(6,*)'OTC OBC',lk,JPF,P(LK,5),MOVY
c          IF (JPF .lt. 3) THEN
c             P(LK,7)=0.0
c             P(LK,4)=0.0
c             P(LK,5)=0.0
c             IPF   = 1
c             JPF   = 1
c             ZKPF  = 1
c             isactive(lk) = 0
c             nactive=nactive-1
cc            if(lk.eq.34887) write(6,*)IPF,JPF
c          ENDIF
c         endif


c-------------------------------------------------------------
c  bounce particles off solid boundaries
c-------------------------------------------------------------
c         If ((pong .eq. 1).and.(isactive(lk).eq.1)) then
          if ((pong.eq.1).and.(JPF.ge.3)) then
C     pong x
c            write(6,*) 'pong x',ipf,jpf
            IF (FSM(IPF,JPF) .EQ. 0.0) THEN
c            write(6,*) 'PONGING ',lk,ip,jp,ipf,jpf
            IF (FSM(IPF,JP) .EQ. 0.0) THEN
c      	       print *,'pong x',lk,IP,IPF,P(lk,4)
               IPF=IP
               P(LK,4)=1.0-P(LK,4)
c               print *,'     to ',IPF,P(lk,4)
            ENDIF
C     pong y
c            write(6,*) 'pong y'
            IF (FSM(IP,JPF) .EQ. 0.0) THEN
c      	       print *,'pong y',lk,JP,JPF,P(lk,5)
               JPF=JP
               P(LK,5)=1.0-P(LK,5)
c               print *,'     to ',JPF,P(lk,5)
            ENDIF
            
            IF ((FSM(IP,JPF) .EQ. 1.0).AND.(FSM(IPF,JP).EQ.1.0)) THEN
c               write(6,*) 'mov pong'
               IF (MOVY .GE. MOVX) THEN
                  IPF=IP
                  P(LK,4)=1.0-P(LK,4)
               ELSE
                  JPF=JP
                  P(LK,5)=1.0-P(LK,5)
               ENDIF
            ENDIF
c    if still trapped then return to original position
c         if (fsm(ipf,jpf).eq.0) then 
c           write(6,*) 'saving ',lk,p10,p20
c           p(lk,1)=p10
c           p(lk,4)=p40
c           p(lk,2)=p20
c           p(lk,5)=p50
c         endif
         ENDIF

C     pong z
c            write(6,*) 'pong z'
            IF (ZKPF .ge. 9.0) THEN  ! randomize bounce to levels 
c              call random_number(sx)
c               sx = sx*2.+7.
c               ZKPF = IFIX(sx)
c               P(LK,6) = sx-zkpf  !1.0-P(LK,6)
               ZKPF=ZKP
               P(LK,6)=1.0-P(LK,6)

            ENDIF
            IF (ZKPF .lT. 2.0) THEN
c               ZKPF=2
c               P(LK,6)=0.75  !1.0-P(LK,6)
               ZKPF=ZKP
               P(LK,6)=1.0-P(LK,6)

cc     print *,'pong z top',ZKP,'zkp'   
            ENDIF
         endif  ! pong
         


c-------------------------------------------------------------
c  particles stick to solid boundaries
c-------------------------------------------------------------
         If (stick .eq. 1) then
C     stick x
            IF (FSM(IPF,JP) .EQ. 0.0) THEN
               IPF=IP
               P(LK,4)=0.0
c     print *,'stick x',IP,'ip'   
            ENDIF
C     stick y
            IF (FSM(IP,JPF) .EQ. 0.0) THEN
               JPF=JP
               P(LK,5)=0.0
c     print *,'stick y',JP,'jp'   
            ENDIF
C     stick z
            IF (ZKPF .GT. 10.0) THEN
               ZKPF=ZKP
               P(LK,6)=0.0
c     print *,'stick z bottom',ZKP,'zkp'   
            ENDIF
            IF (ZKPF .LT. 1.0) THEN
               ZKPF=ZKP
               P(LK,6)=0.0
c     print *,'stick z top',ZKP,'zkp'   
            ENDIF
c            write(6,*) 'end pong'
         endif 
         
              

c********TEST******** set k=4.0 for testing horizontal advection of particles
c         p(lk,3) = 4.
c         p(lk,6) = 0.
c*******************


c-----------------------------------------------------
c  no return condition at open boundary
c  *** should be before pong ***
c-----------------------------------------------------
c         write(6,*) 'LK OBC',lk,ob
         If (ob .eq. 1) then  ! tampa bay

         if (LK.eq.1) write(6,*)'TAMPA OBC'
c      baywide model
            IF (JPF .EQ. 20)  THEN
               IF ((IPF .LE. 19).and.(IPF .GE. 17))  THEN
                  P(LK,7)=0.0
                  P(LK,4)=0.0
                  P(LK,5)=0.0
                  IPF   = 1
                  JPF   = 1
                  ZKPF  = 1
               ENDIF
            ENDIF
            IF ((IPF .EQ. 23).and.(JPF .EQ. 19))  THEN
               P(LK,7)=0.0
               P(LK,4)=0.0
               P(LK,5)=0.0
               IPF   = 1
               JPF   = 1
               ZKPF  = 1
            ENDIF

            IF (IPF .EQ. 29)  THEN
               IF ((JPF .LE. 18).and.(JPF .GE. 17))  THEN
                  P(LK,7)=0.0
                  P(LK,4)=0.0
                  P(LK,5)=0.0
                  IPF   = 1
                  JPF   = 1
                  ZKPF  = 1
                  isactive(lk) = 0
                  nactive=nactive-1
               ENDIF
            ENDIF
            
            IF (JPF .EQ. 15)  THEN
               IF ((IPF .LE. 39).and.(IPF .GE. 33))  THEN
                  P(LK,7)=0.0
                  P(LK,4)=0.0
                  P(LK,5)=0.0
                  IPF   = 1
                  JPF   = 1
                  ZKPF  = 1
                  isactive(lk) = 0
                  nactive=nactive-1
               ENDIF
            ENDIF
         endif  ! ob

c  now moved before pong
         if (ob.eq.3) THEN  ! OTB baseline
c         if (lk.eq.1) write(6,*) 'OTB OBC'
c           if((jpf.le.4).and.(lk.eq.34887)) 
c     &    write(6,*)'OTC OBC',lk,JPF,P(LK,5),MOVY
          IF (JPF .lt. 3) THEN
             P(LK,7)=0.0
             P(LK,4)=0.0
             P(LK,5)=0.0
             IPF   = 1
             JPF   = 1
             ZKPF  = 1
             isactive(lk) = 0
             nactive=nactive-1
c             write(6,*)'OBC ',int,lk,jp,jpf
c            if(lk.eq.34887) write(6,*)IPF,JPF
          ENDIF
          endif  !ob

         P(LK,1)=IPF
         P(LK,2)=JPF
         P(LK,3)=ZKPF
         P(LK,8)=INT
         
c     count total particles
           n1=ifix(p(lk,1))
           n2=ifix(p(lk,2))
           n3=ifix(p(lk,3))
           ipcnt(n1,n2,n3)=ipcnt(n1,n2,n3)+1

 10   CONTINUE  ! lk=1,nump

c-----------------------------------------------------------------------
c     WRITE OUT DATA FOR TRACKING PARTICLE TRAJECTORIES
c-----------------------------------------------------------------------
      IF ((MOD(SEQ-iccctime(mp),IPR).EQ.0).and.
     &     (startpart.eq.0)) THEN
c         write(6,*) 'writing particles at time step ',seq
c         write(iures,*) 'Time'
         write(iures,*) int
         do 30 lk=1,nump,NPSKIP
            WRITE(IURES,177) p(lk,1)+p(lk,4),p(lk,2)+p(lk,5),
     &           p(lk,3)+p(lk,6)
 30      continue

c     write to pcount
c         WRITE(IUTR,*)  'TIME STEP'
         WRITE(IUTR,*)  SEQ
c         WRITE(IUTRS,*)  'TIME STEP'
c         WRITE(IUTRS,*)  SEQ
         DO 13 k=1,kbm1
         do 13 j=1,jm
         do 13 i=1,im
         WRITE(IUTR,*)  ipcnt(I,J,K)
c         WRITE(IUTRS,*)  ipcnt2(I,J,K)
c         if (ipcnt(i,j,k).gt.0) write(6,*) i,j,k,ipcnt(i,j,k)
 133        format(i4)
 13      continue

         ENDIF  ! MOD

 177  format(3(f8.3))
c 133     format(4(i6,2x),i8,2x,3(2x,f6.4))
c 134     format(i10,2x,i6,2x,3(2x,f7.4))



      close(IUTR)
c      close(IUTRS)
      close(IURES)
      RETURN
      END
      
  

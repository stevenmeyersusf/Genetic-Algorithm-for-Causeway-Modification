c23456789012345678901234567890123456789012345678901234567890123456789012   
c for present-day if n=100 then nump=2273000 if 1879 then nump=2360000
      SUBROUTINE DRIFTERc(mp)
C     VERSION(06/21/99)
c     VERSION(07/5/2002) BY S.D. MEYERS
c     modified for continuous release of particles by SD Meyers Feb 2013
c     nrel = # time steps between release
c     nact # active particles
c     nlast time step of last release

      INCLUDE 'comdeck'  

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
      DIMENSION BUP(NUMP),BUP1(NUMP),BUP2(NUMP)
     .,TUP(NUMP),TUP1(NUMP),TUP2(NUMP),UP(NUMP)
     .,BVP(NUMP),BVP1(NUMP),BVP2(NUMP)
     .,TVP(NUMP),TVP1(NUMP),TVP2(NUMP),VP(NUMP)
     .,BWP(NUMP),BWP1(NUMP),BWP2(NUMP)
     .,TWP(NUMP),TWP1(NUMP),TWP2(NUMP),WP(NUMP)
     .,PH2(NUMP),PH1(NUMP),PH3(NUMP)
     .,ut(im,jm,kb),vt(im,jm,kb),wt(im,jm,kb),IPCNT(IM,JM,KBM1),
     .IPCNT2(IM,JM,KBM1)
 
      REAL MOVX,MOVY,MOVZ,clev,dispx,dispy,dispz,dxp,dyp,dzp
      REAL x0,y0,z0  ! release grid position
      INTEGER IP,JP,ZKP,PK,IPF,JPF,ZKPF,SEQ,IPR,INTFIL
      INTEGER II,JJ,MIJ,IX,IY,IZ
      INTEGER pong,stick,ob,pstart,ntemp1,ntemp2,ntemp3,slingshot
      INTEGER LK
      INTEGER mp  ! index for particle start time
c      integer fileopen(nump)
      character*80 case
      character*2 si,sj
      character*1 s1
      character*13  root
      character*17  tfile,tfile0
      character*35 outpath

      SEQ=INT
      outpath = '/DATA3/meyers/Sci/ECOM/Test/Output/'

      WRITE(UNIT=case, FMT='(I1)') mp

c      root = 'Output/tracks'//case
      tfile0 = 'dummy'  ! dummy initial value

      nrel = 5
      IUTR=105
      IUTRS=104
      IURES=106
      INTFIL = 100
c	Write(6,*) startpart 
      IF (startpart.EQ. 1) THEN
        
c        OPEN(IUTR,FILE=outpath//'pcount_mm1_'//case)    ! number of particles in cell
c        OPEN(IUTRS,FILE=outpath//'pcount2_mm1_'//case)  ! num of orig part remaining
        OPEN(IURES,FILE=outpath//'tracks_nwlens_'//case,
     &       form='formatted')
        write(6,*) 'starting particle tracking at ',int
        nact = 0
        nlast = int-nrel
      ELSE
c        OPEN(IUTR,FILE=outpath//'pcount_mm1_'//case,
c     &        access='append')
c       OPEN(IUTRS,FILE=outpath//'pcount2_mm1_'//case,
c     &       access='append')
        OPEN(IURES,FILE=outpath//'tracks_nwlens_'//case,
     &       access='append')
      ENDIF                     ! startpart=1

c the number of time steps between outputs IPR
      IPR=60
      PONG=1
      stick=0
      ob=1
      slingshot = 1
      clev=1.0
      pstart=1

      dispx=0.005
      dispy=0.005
      dispz=0.0

      do 3 i=1,im
      do 3 j=1,jm
      do 3 k=1,kb
         ut(i,j,k)=u(i,j,k)*dum(i,j)
         vt(i,j,k)=v(i,j,k)*dvm(i,j)
         wt(i,j,k)=w(i,j,k)
 3    continue
 
      DO 1 LK=1,NUMP
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
c      do 100 i=1,im
c      do 100 j=1,jm
c      do 100 k=1,kbm1
c         ipcnt(i,j,k)=0
c         ipcnt2(i,j,k)=0
c 100  continue


      OPEN(INTFIL,FILE='initpos_lens.dat',status='old')
      read(intfil,*) x0,y0,z0
c      write(6,*) x0,y0,z0
      CLOSE(INTFIL)
         
c------------------------------------------------------------
c add particle
c------------------------------------------------------------
c      write(6,*) int,nlast,nrel,nact
      IF (((int-nlast).eq.nrel).and.(nact.lt.nump)) THEN 
         nact = nact+1
c         write(6,*)'INITIALIZING NEW DRIFTER INT=',INT
         seedx=-1
         seedy=-2
         seedz=-3
c         write(6,*) 'nact=',nact
         lk = nact
         P(LK,1) = FLOAT(IFIX(x0)) ! cell number
         P(LK,2) = FLOAT(IFIX(y0))
         P(LK,3) = FLOAT(IFIX(z0))
         P(LK,4) = x0-p(lk,1)
         P(LK,5) = y0-p(lk,2)
         P(LK,6) = z0-p(lk,3)
         P(LK,7) = 1.0
         P(LK,8) = 0.0
         P(LK,9) = 1.0
         P(LK,10) = INT
c         IORIG(LK) = II
c         JORIG(LK) = JJ
c         KORIG(LK) = K
         nlast = int
      ENDIF                     ! nact<nump
C 

c-------------------------------------------------------------
c     do particle counts (pcount)
c-------------------------------------------------------------
c      do 23 lk=1,nact 
c         ntemp1=ifix(p(lk,1))
c         ntemp2=ifix(p(lk,2))
c         ntemp3=ifix(p(lk,3))
c         ipcnt(ntemp1,ntemp2,ntemp3)=ipcnt(ntemp1,ntemp2,ntemp3)+1
c         write(6,*) p(lk,1),p(lk,2),p(lk,3)
c         write(6,*) ntemp1,ntemp2,ntemp3,ipcnt(ntemp1,ntemp2,ntemp3)
c 23   continue
c      DO 33 k=1,kbm1
c      do 33 j=1,jm
c      do 33 i=1,im
cc         if (ipcnt(i,j,k).gt.0) write(6,*) i,j,k,ipcnt(i,j,k)
c         ipcnt(i,j,k)=0  ! reset
c 33   continue


c-------------------------------------------------------------
c     increment position by velocity trilinear interpolation
c-------------------------------------------------------------
C     For every particle               
      DO 10 LK=1,nact         
         
         IP=IFIX(P(LK,1))
         JP=IFIX(P(LK,2))
         ZKP=IFIX(P(LK,3))
         
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

c      if (abs(ut(ip+1,jp,zkp+1)).gt.0) 
c     & write(6,*) 'getting 1',lk,UT(IP+1,JP,ZKP+1),ip,jp,zkp
c      if (abs(ut(ip,jp,zkp+1)).gt.0) 
c     & write(6,*) 'getting 2',lk,UT(IP,JP,ZKP+1),ip,jp,zkp
               
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
               
C     calculate the lower two linear interp U's 
               
               bup1(LK)=UT(IP+1,JP,ZKP)+(P(LK,5)-0.5)
     &              *(UT(IP+1,JP+1,ZKP)-UT(IP+1,JP,ZKP))
               bup2(LK)=UT(IP,JP,ZKP)+(P(LK,5)-0.5)
     &              *(UT(IP,JP+1,ZKP)-UT(IP,JP,ZKP))
               
C     combine by linearly interp those two u's
               
               bup(LK)=bup2(LK) + P(LK,4)
     &              *(bup1(LK)-bup2(LK))
               
c       if (mod(lk,50).eq.0) write(6,*) lk,bup1(lk),bup2(lk),bup(lk)

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
               
c        if (mod(lk,100).eq.0) then 
c            write(6,*) 'UP ',LK,ip,jp,zkp,up(lk)
c            write(6,*) bup(lk),tup(lk)
c         endif
            ENDIF
            
         ENDIF
         
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
         MOVX=UP(LK)*p(lk,7)*p(lk,9)*DTI/ph1(lk)
c         MOVX=MOVX+MOVX*dispx
         MOVX=MOVX+(rand(seedx)-0.5)*dispx*p(lk,7)*p(lk,9)
         
         MOVY=VP(LK)*p(lk,7)*p(lk,9)*DTI/ph2(lk)
c         MOVY=MOVY+MOVY*dispy
         MOVY=MOVY+(rand(seedy)-0.5)*dispy*p(lk,7)*p(lk,9)
         
         MOVZ=WP(LK)*p(lk,7)*p(lk,9)*DTI/ph3(lk)*clev
c         MOVZ=MOVZ+MOVZ*dispz
         MOVZ=MOVZ+(rand(seedz)-0.5)*dispz*p(lk,7)*p(lk,9)
         
c     print *,rand(seed),lk

C     Update the I, DX position
c         if (mod(lk,100).eq.0) then 
c            write(6,*) 'U ',LK,UP(LK),MOVX,p(lk,7),p(lk,9)
c         endif
         p10 = p(lk,1)
         p40 = p(lk,4)
         P(LK,4)=P(LK,4)+MOVX
         P(LK,1)=P(LK,1)+IFIX(P(LK,4))
         P(LK,4)=P(LK,4)-IFIX(P(LK,4))
         IF (P(LK,4) .LT. 0.0) THEN
            P(LK,4)=1.0+P(LK,4)
            P(LK,1)=P(LK,1)-1.0
         ENDIF 
         
C     Update the J, DY position
c         if (mod(lk,100).eq.0) then 
c            write(6,*) 'V ',LK,VP(LK),MOVY,p(lk,7),p(lk,9)
c         endif
         p20 = p(lk,2)
         p50 = p(lk,5)
         P(LK,5)=P(LK,5)+MOVY
         P(LK,2)=P(LK,2)+IFIX(P(LK,5))
         P(LK,5)=P(LK,5)-IFIX(P(LK,5))
         
         IF (P(LK,5) .LT. 0.0) THEN
            P(LK,5)=1.0+P(LK,5)
            P(LK,2)=P(LK,2)-1.0
         ENDIF 
         
C     Update the K, DZ position
         P(LK,6)=P(LK,6)+MOVZ
         P(LK,3)=P(LK,3)-IFIX(P(LK,6))
         P(LK,6)=P(LK,6)-IFIX(P(LK,6))
         
         IF (P(LK,6) .LT. 0.0) THEN
            P(LK,6)=1.0+P(LK,6)
            P(LK,3)=P(LK,3)+1.0
         ENDIF 
         IPF=IFIX(P(LK,1))
         JPF=IFIX(P(LK,2))
         ZKPF=IFIX(P(LK,3))
         
         If (pong .eq. 1) then
C     pong x
            IF (FSM(IPF,JPF) .EQ. 0.0) THEN
c            write(6,*) 'PONGING ',lk,ip,jp,ipf,jpf
            IF (FSM(IPF,JP) .EQ. 0.0) THEN
c      	       print *,'pong x',lk,IP,IPF,P(lk,4)
               IPF=IP
               P(LK,4)=1.0-P(LK,4)
c               print *,'     to ',IPF,P(lk,4)
            ENDIF
C     pong y
            IF (FSM(IP,JPF) .EQ. 0.0) THEN
c      	       print *,'pong y',lk,JP,JPF,P(lk,5)
               JPF=JP
               P(LK,5)=1.0-P(LK,5)
c               print *,'     to ',JPF,P(lk,5)
            ENDIF
            
c            if (fsm(ipf,jpf).eq.0) then 
c              write(6,*) 'STILL ON LAND'
c            endif

            IF ((FSM(IP,JPF) .EQ. 1.0).AND.(FSM(IPF,JP).EQ. 1.0)) THEN
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
            IF (ZKPF .GT. 10.0) THEN
c      print *,'pong z bottom',ZKP,ZKPF
               ZKPF=ZKP
               P(LK,6)=1.0-P(LK,6)
            ENDIF
            IF (ZKPF .LT. 1.0) THEN
               ZKPF=ZKP
               P(LK,6)=1.0-P(LK,6)
c     print *,'pong z top',ZKP,'zkp'   
            ENDIF
         endif  ! pong
         
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
         endif
         
         If (ob .eq. 1) then
C     open boundary index
c     East Bay Model
c            IF (IPF .EQ. 2)  THEN  
c               IF ((JPF .LE. 26).and.(IPF .GE. 25))  THEN

c      baywide model
            IF (JPF .EQ. 20)  THEN
               IF ((IPF .LE. 19).and.(IPF .GE. 17))  THEN
                  P(LK,7)=0.0
                  P(LK,4)=0.0
                  P(LK,5)=0.0
c     print *,'ob ',IPF,JPF   
               ENDIF
            ENDIF
            IF ((IPF .EQ. 23).and.(JPF .EQ. 19))  THEN
               P(LK,7)=0.0
               P(LK,4)=0.0
               P(LK,5)=0.0
c     print *,'ob ',IPF,JPF   
            ENDIF

            IF (IPF .EQ. 29)  THEN
               IF ((JPF .LE. 18).and.(JPF .GE. 17))  THEN
                  P(LK,7)=0.0
                  P(LK,4)=0.0
                  P(LK,5)=0.0
c     print *,'ob ',IPF,JPF   
               ENDIF
            ENDIF
            
            IF (JPF .EQ. 15)  THEN
               IF ((IPF .LE. 39).and.(IPF .GE. 33))  THEN
                  P(LK,7)=0.0
                  P(LK,4)=0.0
                  P(LK,5)=0.0
c     print *,'ob ',IPF,JPF   
               ENDIF
            ENDIF
         endif  ! ob
         
         P(LK,1)=IPF
         P(LK,2)=JPF
         P(LK,3)=ZKPF
         P(LK,8)=INT
         
c   count total particles
           ntemp1=ifix(p(lk,1))
           ntemp2=ifix(p(lk,2))
           ntemp3=ifix(p(lk,3))
           ipcnt(ntemp1,ntemp2,ntemp3)=ipcnt(ntemp1,ntemp2,ntemp3)+1

c   count original particles
c           if (ntemp1.eq.iorig(lk).and.ntemp2.eq.jorig(lk).and.
c     &          ntemp3.eq.korig(lk)) ipcnt2(ntemp1,ntemp2,ntemp3)=
c     &          ipcnt2(ntemp1,ntemp2,ntemp3)+1

c           write(6,*) 'lk',lk,seq,p(lk,1),p(lk,2),
c     &          p(lk,3)
c           write(6,*) ntemp1,ntemp2,ntemp3,ipcnt(ntemp1,ntemp2,ntemp3)
c           write(6,*) iorig(lk),jorig(lk),korig(lk)

 10   CONTINUE  ! lk=1,nact

c      write(6,*) 'finished initial time step'
c      nfill=0

c     WRITE OUT DATA FOR TRACKING PARTICLE TRAJECTORIES
      IF (MOD(SEQ,IPR).EQ.0) THEN
         WRITE(6,*)  'drifter time =',SEQ,nact
          write(iures,*) int,nact
         do 30 lk=1,nact
            I = IFIX(p(lk,1))
            J = IFIX(p(lk,2))
            K = IFIX(p(lk,3))
            dxp = p(lk,4)
            dyp = p(lk,5)
            dzp = p(lk,6)
          WRITE(IURES,177) lk,float(i)+dxp,float(j)+dyp,
     &           float(k)+dzp             
 30      continue
      ENDIF

 177  format(i8,3(f10.3))
c 133     format(4(i6,2x),i8,2x,3(2x,f6.4))
c 134     format(i10,2x,i6,2x,3(2x,f7.4))


c     tEST
c         write(6,*) '**************************'
c         write(6,*) 'CHECK drifter count t=',seq
c         DO 14 k=1,kbm1
c         do 14 j=1,jm
c         do 14 i=1,im
c 14         if (ipcnt(i,j,k).gt.0) write(6,*) seq,i,j,k,ipcnt(i,j,k)



c     WRITE OUT DATA FOR PARTICLE COUNT
c      IF (startpart.eq.1.or.(MOD(SEQ,IPR).EQ.0)) THEN
cc         write(6,*) '**************************'
c         write(6,*) 'writing drifter count t=',seq,ipcnt(30,20,1)
c         WRITE(IUTR,*)  'TIME STEP'
c         WRITE(IUTR,*)  SEQ
c         WRITE(IUTRS,*)  'TIME STEP'
c         WRITE(IUTRS,*)  SEQ
c         DO 13 k=1,kbm1
c         do 13 j=1,jm
c         do 13 i=1,im
c         WRITE(IUTR,*)  ipcnt(I,J,K)
c         WRITE(IUTRS,*)  ipcnt2(I,J,K)
cc         if (ipcnt(i,j,k).gt.0) write(6,*) i,j,k,ipcnt(i,j,k)
c 133        format(i4)
c 13      continue
cc         write(6,*) '**************************'
c      ENDIF

c      close(IUTR)
c      close(IUTRS)
      close(IURES)
      RETURN
      END
      
      
      


c Initial fitness test for OTB genetic algorithm May 14, 2022
c This test is # particles in remaining in 
c Feather Sound region after ttest days
c lower count is higher fitness. Local, counting only in region (ifs,jfs). 

c for main project make this part of galgorithm.f

       subroutine ftest1(trackfile,fscr)
c      file: particle file from drifter3.f
c      ttest: time to extract particle count (d), defined in genmain.f
c             should equal release time
c      fscr: fitness score

c       parameter(nfs = 55)   ! # gridcells in ALL test region
       parameter(nfs=18)     ! # gridcells in South test region
       include 'comdeck'
       include 'genetic'
       integer t0,i,j,k,ii,icnt,ifs(nfs),jfs(nfs),ifslist(numpx)
       integer ipart
       real xp,yp,zp,fscr
       character*80 trackfile


c    study region should be separate subroutine

c    define cells in which to obtain particle counts
c    ALL original FS cells (55)
c       data ifs/6,6,6,7,7,7,7,8,8,8,8,8,8,8,
c     & 9,9,9,9,9,9,9,9,
c     & 10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,
c     & 12,12,12,12,12,12,12,12,13,13,13,13,14/  ! feather sound array indices
c       data jfs/25,26,27,24,25,26,27,21,22,23,24,25,26,27,
c     & 20,21,22,23,24,25,26,27,
c     & 21,22,23,24,25,26,27,28,29,19,20,21,22,23,24,25,26,27,28,29,
c     & 20,21,22,23,24,25,26,27,20,21,22,23,21/

C    Southern portion of FS (18)
       data ifs/8,8,9,9,9,10,10,11,11,11,11,
     & 12,12,12,13,13,13,14/  ! feather sound array indices
       data jfs/21,22,20,21,22,21,22,19,20,21,22,
     & 20,21,22,20,21,22,21/    

C    SE portion of WCCC (12)
c      data ifs/12,12,12,13,13,13,14,14,14,15,15,15/
c      data jfs/29,30,31,29,30,31,29,30,31,29,30,31/


       fscr = -9999 ! default
       icnt = 0     ! initialize particle count
c-----------------------------------------------------
c loop tracks file, find test-time, count particles
c-----------------------------------------------------
       write(6,*) trackfile 
       open(12,file=trackfile,status='old')
       read(12,*) nump  ! # particles
c       write(6,*) 'ftest0 nump=',nump


c----------------------------------------------
c   initialize counters and flags
c----------------------------------------------
       itm=0   ! flag to control reading 'start' of particles
       icnt=0
       icnt2=0
       t0=0    ! particle time, initialize
       n0=0    ! # particles in 'start' time

c-----------------------------------------------
c  loop through tracks file and count # in FS
c  starting at min time and then at time ttest
c----------------------------------------------- 
       do while (t0.lt.ttest)   ! loop file until ttest
         read(12,*,END=100) t0
c         write(6,*) t0,tmin,ttest
         if (t0.lt.ttest) then  ! skip this data
           ipart = 0
           do ii=1,nump
             read(12,*) xp,yp,zp
             ipart = ipart+1
c            check if first time past min time
             if ((t0.ge.tmin).and.(itm.eq.0))then  ! if first, get # particles
c               write(6,*) 'finding n0...'
               do jj=1,nfs   ! loop FS cells
                if (ifix(xp).eq.ifs(jj).and.
     &            (ifix(yp).eq.jfs(jj))) then 
                    n0 = n0+1
                    ifslist(n0) = ipart
c                    write(6,*) n0,ifslist(n0),ipart
                endif
               enddo
             endif
           enddo !ii
           if ((t0.ge.tmin).and.(itm.eq.0)) then 
c              write(6,*)'n0=',n0
              itm=1  ! switch flag
           endif

         else ! time is >= ttest
           if (icnt2.eq.0) then !only one test
            do ii=1,nump
              read(12,*),xp,yp,zp
              isin = 0  ! flag if originally in FS
              do icheck=1,n0
                if (ifslist(icheck).eq.ii) isin=1
              enddo
              do ip=1,nfs
               if ((ifix(xp).eq.ifs(ip)).and.
     &         (ifix(yp).eq.jfs(ip)).and.
     &         (isin.eq.1))  icnt=icnt+1
              enddo
            enddo ! ii
            icnt2=1
          endif ! cnt>0
         endif  ! t0 < ttest

       enddo
 100   continue
      close(12)
 
 200  format(80a1)

c----------------------------------------------------
c   compute fitness score
c----------------------------------------------------
c      write(6,*) 'computing fitness ',sid
      fscr=real(n0-icnt)/n0  ! fraction flushed out
      write(6,*) 'results:',icnt,n0,fscr
   
      open(10,file='rankfilename',status='old')
      read(10,*) irlen
      read(10,*) rankfile
      close(10)
      open(10,file=rankfile(1:irlen),status='old',access='append')
      write(10,1500) sid,igen,fscr
      close(10)

 1500 format(a9,2x,i3,2x,f10.5)   ! 'a' format must = 'nc'

      return
      end


c Sitness test for OTB genetic algorithm June 10, 2022
c This test is # particles in remaining in 
c Feather Sound region after ttest days
c lower count is higher fitness. Local, counting only in region (ifs,jfs). 
c smaller FS region than ftest1.f

c for main project make this part of galgorithm.f

       subroutine ftest4(trackfile,fscr)
c      file: particle file from drifter3.f
c      ttest: time to extract particle count (d), defined in genmain.f
c             should equal release time
c      fscr: fitness score

       parameter(nfs = 6)   ! # gridcells in test region
       include 'comdeck'
       include 'genetic'
       integer t0,i,j,k,ii,icnt,ifs(nfs),jfs(nfs),ifslist(numpx)
       integer ipart
       real xp,yp,zp,fscr

c    define cells in which to obtain particle counts
       data ifs/14,14,15,15,16,16/ ! feather sound array indices
       data jfs/23,24,23,24,23,25/
       character*80 trackfile

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


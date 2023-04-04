c  mutates retained chromosomes and outputs new SID files
c  this is prototype version with hard-coded options
c  redefine npop value

c   BUGGY, replaced by mutate2.f

      subroutine mutate()

      include 'genetic'
      include 'comdeck'

      integer iopen(2)         ! # open cells in span
      integer iranksort(npopmax)
      integer iset(2,npopmax)  ! indices of open cells
      real rankid(npopmax),rarry(nc)
      real ranksort(npopmax)
      character*1 s1
      character*9 sid0,sidlist2(npopmax) 
      character*2 sp       ! index for new SID file


c     1.  get sid from last generation
c     2.  identify open cells, # open is icnt
c     3.  loop the open cells...
c     3a.  close cell
c     3b.  randomly open another
c     3c.  check for duplicates 
c     4.  save to new SID
c     5.  replace one with elite

      nkeep = 2   ! # of highest scoring SIDs to keep
      ncell = nc
c-----------------------------------------------
c     loop the retained SIDs and mutate
c-----------------------------------------------
      write(6,*) 'starting mutate',nkeep,nc,igen
      
      igen = igen+1       ! generation # to create
      write(6,*) 'CREATING GENERATION ',igen
      npop2 = 0          ! index for SID files to be created

c     keep was defined in roulette.f
      do k=1,nkeep        ! loop # sid to keep`
c   1.  get sid
        sid(1:nc) = sidlist(keep(k))   ! get SID string
        icnt=0     ! closed cells
        write(6,*) k,keep(k),sid(1:nc)
        write(6,*) sid

c   2.  identify which grid cells are open
        do i=1,nc
          s1 = trim(sid(i:i))
          if (trim(s1).eq.'1') then
            icnt = icnt+1
            iopen(icnt) = i  ! list currently open
          endif
        enddo
        write(6,*) 'partitioning ',icnt,sid(1:nc),iopen(1),iopen(2),igen

c    3. close each iopen in sequence and randomly open another
        do io=1,icnt  ! loop each open
          do idt=1,1  ! loop daughters of io
          sid0 = sid(1:nc)
          ialt = 2  ! alternative index, hardcoded
          if (io.eq.2) ialt=1

c   3a.   close this cell
          ii = iopen(io)
          
          sid0(ii:ii) = '0'       
c
c   3b.   randomly choose new cell to open
 550      call genran2(rankid,iranksort,ranksort,ncell)
          ir = 1 ! index for which sorted value to use
          inew = iranksort(ir)   
c         if this cells is already open, then change  
          do while ((inew.eq.iopen(ialt)).or.(inew.eq.ii))
            ir = ir+1 
            inew = iranksort(ir)  ! use next cell
          enddo
          sid0(inew:inew) = '1'         
          npop2 = npop2+1       ! new sid counter

c   3c.   check for duplicates in this generation
          iset(1,npop2) = ialt  ! record indices
          iset(2,npop2) = inew
          sidlist2(npop2) = sid0
          continue
          if (npop2.gt.1) then
            do ip2=1,npop2-1  ! loop previous and check
               isdiff=1       ! start by assuming is different

c         this method is specific to prototype (2 open cells)
              if (((iset(1,npop2).eq.iset(1,ip2)).and.
     &         (iset(2,npop2).eq.iset(2,ip2)))
     &        .or.
     &         ((iset(1,npop2).eq.iset(2,ip2)).and.
     &          (iset(2,npop2).eq.iset(1,ip2))))
     &          isdiff=0

              if (isdiff.eq.0) then
                sid0(inew:inew) = '0'
                ir = ir+1
                inew = iranksort(ir) !open new cell, then recheck            
                sid0(inew:inew) = '1'
                npop2 = npop2-1
                go to 550     ! start over with new random set 
              endif
            enddo
          endif  ! npop2>1
          if (npop2-1.lt.10) then 
            write(sp,'(i1)') npop2-1
            sp = '0'//sp
          endif
          if ((npop2-1.ge.10).and.(npop2-1.lt.100)) then
            write(sp,'(i2)') npop2-1
          endif
          if (npop2.ge.100) then 
            write(6,*) 'ERROR in mutate.f: npop too large, ~L54'
            return
          endif  

c     4.  output new SID file
          sidlist(npop2) = sid0
          write(6,*) sid,'->',sid0
          x =system('cd /home/meyers/Sci/ECOM/GA/PT/Code')
          open(10,file='SID_'//sp//'.dat')
          write(10,11) sid0
          write(10,21) igen,ttest,tmin
          close(10)
          enddo ! idt
        enddo  !  io
       enddo   ! k

c  5 . option for elitism, identified in rdfitness.f
       write(6,*) 'IELITE=',ielite
       if (ielite.gt.0) then   
c        check if elite already in population
         isincl = 0 ! assume not included
         do iid = 1,npop2
           sid0 = trim(sidlist(iid))
           if (trim(sid0).eq.trim(esid)) then 
             isincl=1
             write(6,*) 'Elite is included'
           endif
         enddo
c        randomly overwrite one of the SIDs, replace with elite
         if (isincl.eq.0) then ! proceed to replace randomly chosen
           write(6,*) 'Elite is NOT included'
           call genran2(rankid,iranksort,ranksort,npop2)
           i0 = iranksort(1)  ! take lowest scoring to replace
           sid0 = trim(sidlist(i0))  ! being replaced
           if (i0.lt.10) then
             write(sp,'(i1)') i0
             sp = '0'//sp
           endif
           if ((i0.ge.10).and.(i0.lt.100)) then
             write(sp,'(i2)') i0
           endif
           open(10,file='SID_'//sp//'.dat')
           write(10,11) esid(1:nc)
           write(10,21) igen,ttest,tmin
           close(10) 
         endif  ! isincl=0
       endif    ! ielite>0     

 11    FORMAT(a9)
 21    format(i4,2x,f10.2,2x,f10.2)
       return
       end

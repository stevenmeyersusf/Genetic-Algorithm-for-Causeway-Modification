c does npop loop to create unique cell configurations

       subroutine mkrandgrid(g)  !,root,suff)

       include 'comdeck'
       include 'genetic'
       character*2 sp
       real rset(nc,npop)
       integer iset(nc,npop)

c-----------------------------------------
c  first generation, stricly random
c-----------------------------------------
c       write(6,*) '1st generation grids...'
C      loop members of population for this generation
       do ip=1,npop  
         write(6,*) 'ip=',ip,npop

c        index for each population member
         ipm1=ip-1   ! use -1 for linux compatability
         if (ipm1.lt.10) then
           write(sp,'(i1)') ipm1
           sp = '0'//sp
         endif
         if (ipm1.ge.10.and.ipm1.lt.100) then
           write(sp,'(i2)') ipm1
         endif
         if (ip.ge.100) then 
            write(6,*) 'ERROR: mkranvals, # pop>npop'
            return
         endif
c         write(6,*) ip,sp  ! population identifier

c        generate random values for each cell
c        prototype r is a vector, full GA will be array
         call genran(r,isort,rsort)

c----------------------------------------------------------------
c     detect and replace duplicates 
c     look for duplicate cells in any order
c----------------------------------------------------------------
 550    continue
        if (ip.gt.1) then
        do ip2=1,ip-1  ! loop previous and check
          isdiff=1     ! assume is unique

c      this method is specific to prototype (2 open cells)
c      write(6,*) ip,ip2,isort(1),isort(2),iset(ip2,1),iset(ip2,2)
          if (((isort(1).eq.iset(1,ip2)).and.
     &         (isort(2).eq.iset(2,ip2)))
     &    .or.
     &       ((isort(1).eq.iset(2,ip2)).and.
     &        (isort(2).eq.iset(1,ip2))))
     &        isdiff=0


          if (isdiff.eq.0) then 
             write(6,*) 'duplicate:',r,isort(1),isort(2),
     &      iset(1,ip2),iset(2,ip2)
             call genran(r,isort,rsort)  ! replace
             go to 550
          endif
        enddo
        endif  ! ip>1

c       save to array for later comparison
        do ii=1,nc
          rset(ii,ip) = r(ii)
          iset(ii,ip) = isort(ii)
        enddo

c------------------------------------------------------
c     define open cells and sid 
c     prototype: only for 1 span, 2 cells in span
c------------------------------------------------------
         icnt=0  ! index for 'id'
         do i=1,nc
           id(i) = 0 ! default setting is closed cell
         enddo
         do n=1,nspan
c           write(6,*) 'n=',n
           if (ispan(n).gt.0) then !is active causeway
c          open two randomly selected cells
             id(isort(1)) = 1
             id(isort(2)) = 1
c             write(6,*) 'FINAL:',ip,isort(1),isort(2)
c
             call mksid()
c             write(6,11) (sid)
           endif  ! isspan>0
         enddo  ! nspan  

c--------------------------------------------
c  output sid to read file
c  bc each ECOM .exe must be run from shell
c--------------------------------------------
         open(10,file='SID_'//sp//'.dat',status='new')
         write(10,11) sid
c         write(6,*) 'igen=',igen,ttest
         write(10,21) igen,ttest,tmin
         close(10)

c        also save according to generation #
c         if (igen.lt.10) then
c          write(sgen,'(i1)') igen
c          sgen = "0"//sgen
c         endif
c         if ((igen.ge.10).and.(igen.lt.100)) then
c          write(sgen,'(i2)') igen
c         endif
c         write(6,*) 'SGEN=',sgen(1:2)
c         open(10,file='SID_'//sp//'_'//sgen(1:2))
c         write(10,*) sid
c         write(10,21)igen,ttest,tmin
c         close(10)

      enddo  ! npop

 11   FORMAT(a9)
 21   format(i4,2x,f10.2,2x,f10.2)
      return
      end


c     sorts rank of current generation
      subroutine rdfitness(igen0)

c     rankfile: file containing fitness test scores
c     igen0: generation # being analyzed
c     sidlist: list of SIDs in this generation
c     ifitsort: index of sids sorted by fitness score
c     icnt: # members in this generation

      include 'genetic'
c      include 'comdeck'

      character*80 COM  ! dummy variable
      integer igen0,igen1,n0
      real flist(npop)
      logical swapped
      escr = -99999   ! elite fscore
      iall = 0   ! index for all sids
c--------------------------------------------------
c     get fitness information
c--------------------------------------------------
      open(10,file='rankfilename',status='old')
      read(10,*) irlen   ! # characters in rank file name
      read(10,*) rankfile
      close(10)

c     write info to rank file
      ipopcnt = 0  ! initialize, # active population
      igen1 = 0    ! dummy variable
      write(6,*) 'opening ',rankfile(1:irlen),igen0
      open(10,file=rankfile(1:irlen),status='old')
      read(10,100) COM
      do while (igen1.le.igen0)
        read(10,1500,end=50) sid,igen1,fscr
        iall = iall+1
        sidall(iall) = sid(1:nc)    ! not used in prototype
        if (igen1.eq.igen0) then    ! examine this generation
          ipopcnt = ipopcnt+1      
          sidlist(ipopcnt) = sid(1:nc)
          fscore(ipopcnt) = fscr    ! raw score
c          write(6,*) 'ipop',ipopcnt,sidlist(ipopcnt),igen1,fscr

c         identify highest scoring sid
          if ((ielite.gt.0).and.(fscr.gt.escr)) then 
            esid(1:nc) = sid(1:nc)
            escr = fscr
            write(6,*) 'NEW elite ',fscr,' ',sid(1:nc)
          endif
        endif  ! igen1=igen0
      enddo
  50  close(10)
 1500 format(a9,2x,i3,2x,f10.5)   ! 'a' format must = 'nc'
  100 format(a80)

      write(6,*) 'ELITE SID=',esid(1:nc)
      if (ipopcnt.ne.npop) then 
c         write(6,*) 'RESETTING npop=',npop,ipopcnt
         npop = ipopcnt
      endif

      igen = igen0
      return
      end


c    weighted roulette selection of survivors
c    adds elite if noted in 'genetic'

      subroutine roulette(nkeep)

c     ielite: flag to use elitism (0-no, 1-yes)
c     slist: 
c     nkeep:  # of survivors 

      include 'genetic'
      include 'comdeck'
      integer ikeep(npopmax)
      real score2(npopmax),ssort2(npopmax),fx
      integer isort2(npopmax)

c---------------------------------------------
c   roulette
c   1. generate npop random numbers 
c   2. calculate weighted fitness
c   3. calculate roulette 
c   4. sort to find nkeep survivors
c   5. create list of survivors
c--------------------------------------------- 

      write(6,*) 'roulette ',ipopcnt,nkeep
c   1. generate list of random #s
      do i=1,ipopcnt
        call random_number(sx)
        r(i) = sx
      enddo
      
c   2. weighted fitness
c     sort to get highest value
      call bubsort(fscore,fsort,ifsort,ipopcnt)
      fx = fsort(ipopcnt)
      ifx = ifsort(ipopcnt)
      do i=1,ipopcnt
        score2(i) = r(i)*fscore(i)/fx
c        write(6,*) i,r(i),fscore(i),fsort(i),score2(i),fx
      enddo

c   3. true roulette normalizes score
c   not needed, just use sort    

c   4. sort by weighted score
c      write(6,*) 'sorting score2'
      call bubsort(score2,ssort2,isort2,ipopcnt)
c      write(6,*) 'size=',shape(score2)

c   5. list of survivors
       do i=1,nkeep
         keep(i) = isort2(ipopcnt-i+1)  ! keep  highest scoring
         write(6,*) 'keep ',ipopcnt-i+1,keep(i),sidlist(keep(i))
       enddo

c---------------------------------------------
c     option to include elite sid index
c---------------------------------------------
c      if (ielite.eq.1) then
c        isnew = 1  ! assume elite not already included
c        do i=1,nkeep
c          if (ifx.eq.keep(i)) isnew=0
c        enddo
c        if (isnew.eq.1) then
c          nkeep = nkeep+1
c          keep(nkeep) = ifx
c        endif
c        write(6,*) keep
c      endif

c      write(6,*) 'finished roulette ',nkeep

      return
      end

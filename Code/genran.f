C23456789012345678901234567890123456789012345678901234567890123456789012

      subroutine genran(r,irsort,rsort)
c     creates random [0,1] array for each active cell and sorted index

      include 'comdeck'
      include 'genetic'
c    r: array of random values for each causeway
c    isort: descending sort of r


      logical swapped

c      call random_seed()
c      write(6,*) 'creating random'
      do  n=1,nspan           ! loop spans
c      write(6,*) 'n=',n
      if (ispan(n).eq.1) then !check if this span editable
        np = ncellspan(n)     ! # points in this span
c        write(6,*) 'SPAN ',n,ncellspan(n),nspan,np 
c       generate random values for each cell in this span
        do j=1,np
          call random_number(sx)  ! call depends on compiler
          r(j) = sx
c          write(6,*) 'make r ',j,r(j)
        enddo
      endif
      enddo

c   use bubble sort 
      call bubsort(r,rsort,irsort,np)

c      write(6,*) 'finished genran'

      return
      end

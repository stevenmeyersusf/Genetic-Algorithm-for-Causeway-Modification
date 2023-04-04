C23456789012345678901234567890123456789012345678901234567890123456789012
       subroutine genran2(rarray,iasort,rasort,na)
c      include 'comdeck'
c      include 'genetic'
c    r: array of random values for each causeway
c    isort: descending sort of r

      include 'comdeck'
      include 'genetic'
c
      real rarray(na),rasort(na)
      integer iasort(na)

c     generate random values for each cell in this span
      do j=1,na
        call random_number(sx)  ! call depends on compiler
        rarray(j) = sx
c        write(6,*) 'make r ',j,rarray(j)
      enddo

c   use bubble sort 
      call bubsort(rarray,rasort,iasort,na)

      return
      end

c      does ascending bubble sort on array
       subroutine bubsort(array,sarray,iarray,n)
       logical swapped
       integer n
       real array(n),sarray(n)

       integer iarray(n)
       do i=1,n
         sarray(i) = array(i)
         iarray(i) = i
       enddo

c       write(6,*) 'sorting array size ',n
c      begin sorting
       do j = n-1,1,-1
        swapped = .false.
        do i = 1, j
          if (sarray(i).gt.sarray(i+1)) then
            temp = sarray(i)  ! sort value
            sarray(i) = sarray(i+1)
            sarray(i+1) = temp
            temp = iarray(i)           ! sort index
            iarray(i) = iarray(i+1)
            iarray(i+1) = temp
            swapped = .true.
           endif
         enddo
         if (.not. swapped) exit
       enddo
       return
       end

